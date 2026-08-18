#!/usr/bin/env python3

"""Convert an FBMN quantification table to tall form with bounded memory.

Rows are emitted deterministically in feature-chunk, abundance-column, then
feature-row order. Changing the chunk size can therefore change row order,
but never the output rows or values.
"""

import argparse
import glob
import os
from pathlib import Path

import pandas as pd


ID_COLUMNS = ["row ID", "row m/z", "row retention time"]
OUTPUT_COLUMNS = [
    "row ID",
    "filename",
    "area",
    "raw_data",
    "row m/z",
    "row retention time",
]
NO_ABUNDANCE_OUTPUT_COLUMNS = [
    "filename",
    "area",
    "raw_data",
    "row ID",
    "row m/z",
    "row retention time",
]
DEFAULT_CHUNK_SIZE = 100


def _peak_area_columns(columns):
    # Preserve the legacy substring match used by this script.
    return [column for column in columns if "Peak area" in column]


def _filename_from_column(column):
    return column.replace(" Peak area", "")


def _raw_filenames(raw_data):
    return {
        os.path.basename(path)
        for path in glob.glob(os.path.join(raw_data, "*"))
    }


def _make_tall_chunk(quant_chunk, abundance_columns, raw_filenames):
    tall = quant_chunk.melt(
        id_vars=ID_COLUMNS,
        value_vars=abundance_columns,
        var_name="filename",
        value_name="area",
    )
    tall["filename"] = tall["filename"].map(_filename_from_column)
    tall["raw_data"] = tall["filename"].isin(raw_filenames).astype(int)
    return tall.loc[:, OUTPUT_COLUMNS]


def create_tall_file(quant_file, raw_data, output_tall, chunk_size):
    columns = pd.read_csv(quant_file, sep=",", nrows=0).columns.tolist()
    missing = [column for column in ID_COLUMNS if column not in columns]
    if missing:
        raise ValueError(f"Quantification table is missing required columns: {missing}")

    abundance_columns = _peak_area_columns(columns)
    required_columns = ID_COLUMNS + abundance_columns
    raw_filenames = _raw_filenames(raw_data)
    matched_raw = {
        _filename_from_column(column)
        for column in abundance_columns
        if _filename_from_column(column) in raw_filenames
    }
    print(f"{len(matched_raw)} raw files found", flush=True)

    output_path = Path(output_tall)
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    input_rows = 0
    output_rows = 0
    try:
        for chunk_number, quant_chunk in enumerate(
            pd.read_csv(
                quant_file,
                sep=",",
                usecols=required_columns,
                chunksize=chunk_size,
            ),
            start=1,
        ):
            tall_chunk = _make_tall_chunk(
                quant_chunk, abundance_columns, matched_raw
            )
            tall_chunk.to_csv(
                temporary_path,
                index=False,
                sep="\t",
                mode="w" if chunk_number == 1 else "a",
                header=chunk_number == 1,
            )
            input_rows += len(quant_chunk)
            output_rows += len(tall_chunk)
            if chunk_number == 1 or chunk_number % 100 == 0:
                print(
                    f"[create_tall_quant] processed {input_rows:,} feature rows; "
                    f"wrote {output_rows:,} tall rows",
                    flush=True,
                )

        if input_rows == 0:
            pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(
                temporary_path, index=False, sep="\t"
            )
        elif not abundance_columns:
            # With no value columns, the legacy melt+merge path placed the
            # empty melt columns before the merge key.
            pd.DataFrame(columns=NO_ABUNDANCE_OUTPUT_COLUMNS).to_csv(
                temporary_path, index=False, sep="\t"
            )
        os.replace(temporary_path, output_path)
        print(
            f"[create_tall_quant] complete: processed {input_rows:,} feature rows; "
            f"wrote {output_rows:,} tall rows",
            flush=True,
        )
    except BaseException:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Create a tall quantification table for visualization and linkout"
    )
    parser.add_argument("quant_file", help="quant file")
    parser.add_argument("raw_data", help="raw data folder")
    parser.add_argument("output_tall", help="output file")
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help=(
            f"feature rows processed at once (default: {DEFAULT_CHUNK_SIZE}); "
            "changing it can change output row order"
        ),
    )
    args = parser.parse_args()
    if args.chunk_size < 1:
        parser.error("--chunk-size must be at least 1")

    create_tall_file(args.quant_file, args.raw_data, args.output_tall, args.chunk_size)


if __name__ == "__main__":
    main()
