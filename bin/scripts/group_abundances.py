#!/usr/bin/env python3

"""Create the FBMN cluster summary with bounded memory use."""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


BASE_COLUMNS = ["row ID", "row m/z", "row retention time"]
DEFAULT_CHUNK_SIZE = 250


def create_attribute_group_list(metadata_df):
    """Return ATTRIBUTE_ columns and their values in legacy output order."""
    all_attributes = [
        column
        for column in metadata_df.columns
        if column.upper().startswith("ATTRIBUTE_")
    ]
    all_attribute_groups = []
    for attribute in all_attributes:
        for group in metadata_df[attribute].unique().tolist():
            all_attribute_groups.append({"attribute": attribute, "group": group})
    return all_attributes, all_attribute_groups


def _clean_metadata(metadata_df):
    """Apply the filename normalization used by the original implementation."""
    if len(metadata_df) == 0:
        return metadata_df.copy()
    if "filename" not in metadata_df.columns:
        raise Exception("Metadata does not contain filename column")

    cleaned = metadata_df.loc[metadata_df["filename"].notnull()].copy()
    cleaned["filename"] = cleaned["filename"].map(
        lambda value: os.path.basename(str(value)).rstrip()
    )
    return cleaned


def _peak_area_columns(columns):
    return [column for column in columns if column.endswith("Peak area")]


def _column_to_filename(column):
    # Keep the legacy replacement behavior (not suffix removal only).
    return column.replace("Peak area", "").rstrip()


def _build_group_specs(feature_columns, metadata_df):
    """Map every legacy output group to its abundance input columns.

    A column can occur more than once in a spec when duplicate metadata rows
    existed. This intentionally retains the weighting produced by the legacy
    melt+merge implementation.
    """
    _, attribute_groups = create_attribute_group_list(metadata_df)
    abundance_columns = [
        (column, _column_to_filename(column))
        for column in _peak_area_columns(feature_columns)
    ]

    specs = []
    for item in attribute_groups:
        attribute = item["attribute"]
        group = item["group"]
        columns = []

        # pandas groupby drops NA grouping keys by default. The legacy output
        # nevertheless created an NA-named column and filled it with zero.
        if not pd.isna(group):
            filename_counts = metadata_df.loc[
                metadata_df[attribute] == group, "filename"
            ].value_counts(sort=False)
            # Legacy melt emitted abundance columns in feature-table order,
            # then merge repeated each row in matching metadata-row order.
            # Retaining that value order is required for byte-identical
            # floating-point group means.
            for column, filename in abundance_columns:
                columns.extend([column] * int(filename_counts.get(filename, 0)))

        specs.append((f"{attribute}:GNPSGROUP:{group}", columns))
    return specs


def _legacy_group_mean(feature_chunk, abundance_columns):
    """Match pandas groupby mean while retaining bounded memory.

    pandas' groupby mean uses compensated summation. A direct row-wise mean
    changes the last few floating-point bits, which in turn changes TSV bytes.
    Applying the same reduction down each feature row reproduces the legacy
    melt+merge+groupby result without materializing the tall table.
    """
    row_count = len(feature_chunk)
    totals = np.zeros(row_count, dtype=np.float64)
    compensation = np.zeros(row_count, dtype=np.float64)
    counts = np.zeros(row_count, dtype=np.int64)

    for column in abundance_columns:
        values = feature_chunk[column].to_numpy(dtype=np.float64, copy=False)
        valid = ~np.isnan(values)
        adjusted = np.where(valid, values - compensation, 0.0)
        updated = totals + adjusted
        compensation = np.where(valid, (updated - totals) - adjusted, compensation)
        totals = updated
        counts += valid

    means = np.full(row_count, np.nan, dtype=np.float64)
    np.divide(totals, counts, out=means, where=counts != 0)
    return pd.Series(means, index=feature_chunk.index)


def _calculate_chunk(feature_chunk, group_specs):
    output_columns = {
        "cluster index": feature_chunk["row ID"],
        "parent mass": feature_chunk["row m/z"],
        "RTMean": feature_chunk["row retention time"],
    }

    for output_column, abundance_columns in group_specs:
        if not abundance_columns:
            output_columns[output_column] = 0
        else:
            # Repeated names are intentional when metadata contains duplicates.
            output_columns[output_column] = _legacy_group_mean(
                feature_chunk, abundance_columns
            )
    return pd.DataFrame(output_columns, index=feature_chunk.index)


def calculate_groups_metadata(feature_table_df, metadata_df):
    """In-memory API retained for callers and small-data unit tests."""
    cleaned_metadata = _clean_metadata(metadata_df)
    specs = _build_group_specs(feature_table_df.columns, cleaned_metadata)
    return _calculate_chunk(feature_table_df, specs).reset_index(drop=True)


def calculate_groups_file(input_featuretable, input_metadata, output_path, chunk_size):
    """Stream a feature table and append completed cluster-summary chunks."""
    feature_columns = pd.read_csv(input_featuretable, sep=",", nrows=0).columns.tolist()
    missing = [column for column in BASE_COLUMNS if column not in feature_columns]
    if missing:
        raise ValueError(f"Feature table is missing required columns: {missing}")

    try:
        metadata_df = pd.read_csv(input_metadata, sep="\t")
    except Exception:
        metadata_df = pd.DataFrame()
    cleaned_metadata = _clean_metadata(metadata_df)
    group_specs = _build_group_specs(feature_columns, cleaned_metadata)

    required_columns = list(BASE_COLUMNS)
    required_set = set(required_columns)
    for _, columns in group_specs:
        for column in columns:
            if column not in required_set:
                required_columns.append(column)
                required_set.add(column)

    output_path = Path(output_path)
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    rows_written = 0
    try:
        for chunk_number, feature_chunk in enumerate(
            pd.read_csv(
                input_featuretable,
                sep=",",
                usecols=required_columns,
                chunksize=chunk_size,
            ),
            start=1,
        ):
            output_chunk = _calculate_chunk(feature_chunk, group_specs)
            output_chunk.to_csv(
                temporary_path,
                sep="\t",
                index=False,
                mode="w" if chunk_number == 1 else "a",
                header=chunk_number == 1,
            )
            rows_written += len(output_chunk)
            if chunk_number == 1 or chunk_number % 100 == 0:
                print(
                    f"[group_abundances] wrote {rows_written:,} feature rows",
                    flush=True,
                )

        # Preserve the legacy header-only behavior: converting an empty record
        # list back to a DataFrame produced a file containing just a newline.
        if rows_written == 0:
            pd.DataFrame().to_csv(temporary_path, sep="\t", index=False)
        os.replace(temporary_path, output_path)
        print(
            f"[group_abundances] complete: wrote {rows_written:,} feature rows",
            flush=True,
        )
    except BaseException:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass
        raise


def main():
    parser = argparse.ArgumentParser(description="Creates enriched cluster info summary")
    parser.add_argument("input_featuretable", help="input_featuretable")
    parser.add_argument("input_metadata", help="input_group_mapping_filename")
    parser.add_argument(
        "output_clusterinfosummary_filename", help="output cluster-info summary"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help=f"feature rows processed at once (default: {DEFAULT_CHUNK_SIZE})",
    )
    args = parser.parse_args()
    if args.chunk_size < 1:
        parser.error("--chunk-size must be at least 1")

    calculate_groups_file(
        args.input_featuretable,
        args.input_metadata,
        args.output_clusterinfosummary_filename,
        args.chunk_size,
    )


if __name__ == "__main__":
    main()
