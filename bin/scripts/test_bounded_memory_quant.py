import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import create_tall_quant
import group_abundances


SCRIPT_DIR = Path(__file__).parent
BASE = ["row ID", "row m/z", "row retention time"]


def _feature_table(style="minimal"):
    data = {
        "row ID": [11, 12, 13],
        "row m/z": [101.1, 202.2, 303.3],
        "row retention time": [1.5, 2.5, 3.5],
        "alpha.mzML Peak area": [0.0, 20.0, np.nan],
        "beta file.mzML Peak area": [4.0, 0.0, 8.0],
        "unmatched.mzML Peak area": [9.0, 10.0, 11.0],
    }
    frame = pd.DataFrame(data)
    if style == "rich":
        frame.insert(3, "row ion mobility", [1.0, np.nan, 3.0])
        frame.insert(4, "charge", [1, 2, 1])
        frame["source_file"] = ["one", "two", "three"]
        frame["alpha.mzML Peak height"] = [100.0, 200.0, 300.0]
    return frame


def _metadata():
    return pd.DataFrame(
        {
            "filename": [
                "/staged/alpha.mzML ",
                "beta file.mzML",
                "beta file.mzML",
                "not-in-features.mzML",
                "missing-name.mzML",
                None,
            ],
            "ATTRIBUTE_condition": ["case", "case", "control", "control", np.nan, "case"],
            "ATTRIBUTE_numeric": [1, 1, 2, 2, 3, 4],
            "ordinary_column": ["a", "b", "c", "d", "e", "f"],
        }
    )


def _legacy_group(feature_table_df, metadata_df):
    cluster_summary_df = pd.DataFrame()
    cluster_summary_df["cluster index"] = feature_table_df["row ID"]
    cluster_summary_df["parent mass"] = feature_table_df["row m/z"]
    cluster_summary_df["RTMean"] = feature_table_df["row retention time"]
    if len(metadata_df) == 0:
        return cluster_summary_df
    if "filename" not in metadata_df.columns:
        raise Exception("Metadata does not contain filename column")

    metadata_df = metadata_df[metadata_df["filename"].notnull()].copy()
    metadata_df["filename"] = metadata_df["filename"].map(str)
    metadata_df["filename"] = metadata_df["filename"].map(
        lambda value: os.path.basename(value).rstrip()
    )
    filename_columns = [
        column for column in feature_table_df.columns if column.endswith("Peak area")
    ]
    _, attribute_groups = group_abundances.create_attribute_group_list(metadata_df)
    tall = pd.melt(
        feature_table_df,
        id_vars=BASE,
        value_vars=filename_columns,
        var_name="filename",
        value_name="area",
    )
    tall["filename"] = tall["filename"].map(
        lambda value: value.replace("Peak area", "").rstrip()
    )
    tall = tall.merge(metadata_df, on="filename", how="inner")
    records = cluster_summary_df.to_dict(orient="records")
    means = {}
    for attribute in set(item["attribute"] for item in attribute_groups):
        means[attribute] = tall.groupby(
            ["row ID", attribute], observed=True
        )["area"].mean()
    for record in records:
        for item in attribute_groups:
            try:
                value = means[item["attribute"]].loc[
                    (record["cluster index"], item["group"])
                ]
            except KeyError:
                value = 0
            record[
                f'{item["attribute"]}:GNPSGROUP:{item["group"]}'
            ] = value
    return pd.DataFrame(records)


def _legacy_tall(quant_df, raw_dir):
    raw_files = list(Path(raw_dir).glob("*"))
    filenames = []
    for header in quant_df.columns:
        if "Peak area" in header:
            filenames.append(header.replace("Peak area", "").rstrip())
    intersection = {}
    for raw_file in raw_files:
        for filename in filenames:
            if raw_file.name == filename:
                intersection[filename] = str(raw_file)
                break

    columns = [f"{filename} Peak area" for filename in filenames]
    tall = quant_df.melt(
        id_vars=["row ID"],
        value_vars=columns,
        var_name="filename",
        value_name="area",
    )
    tall["filename"] = tall["filename"].str.replace(" Peak area", "", regex=False)
    tall["raw_data"] = 0
    for filename in intersection:
        tall.loc[tall["filename"] == filename, "raw_data"] = 1
    return pd.merge(
        tall,
        quant_df[BASE],
        on="row ID",
        how="left",
    )


def _sort_tall(frame):
    return frame.sort_values(["filename", "row ID"], kind="stable").reset_index(drop=True)


@pytest.mark.parametrize("style", ["minimal", "rich"])
@pytest.mark.parametrize("chunk_size", [1, 2, 100])
def test_group_file_matches_legacy_for_standardized_input_styles(tmp_path, style, chunk_size):
    features = _feature_table(style)
    metadata = _metadata()
    feature_path = tmp_path / "features.csv"
    metadata_path = tmp_path / "metadata.tsv"
    output_path = tmp_path / "groups.tsv"
    features.to_csv(feature_path, index=False)
    metadata.to_csv(metadata_path, sep="\t", index=False)

    group_abundances.calculate_groups_file(
        feature_path, metadata_path, output_path, chunk_size
    )
    actual = pd.read_csv(output_path, sep="\t")
    expected = _legacy_group(features, metadata)
    assert_frame_equal(actual, expected, check_dtype=False)


@pytest.mark.parametrize("metadata_kind", ["missing", "empty", "no_attributes"])
def test_group_file_metadata_variants_match_legacy(tmp_path, metadata_kind):
    features = _feature_table("rich")
    feature_path = tmp_path / "features.csv"
    output_path = tmp_path / "groups.tsv"
    features.to_csv(feature_path, index=False)

    if metadata_kind == "missing":
        metadata_path = tmp_path / "does-not-exist.tsv"
        expected_metadata = pd.DataFrame()
    elif metadata_kind == "empty":
        metadata_path = tmp_path / "metadata.tsv"
        metadata_path.write_text("")
        expected_metadata = pd.DataFrame()
    else:
        metadata_path = tmp_path / "metadata.tsv"
        expected_metadata = pd.DataFrame(
            {"filename": ["alpha.mzML"], "sample_type": ["case"]}
        )
        expected_metadata.to_csv(metadata_path, sep="\t", index=False)

    group_abundances.calculate_groups_file(
        feature_path, metadata_path, output_path, chunk_size=1
    )
    actual = pd.read_csv(output_path, sep="\t")
    expected = _legacy_group(features, expected_metadata)
    assert_frame_equal(actual, expected, check_dtype=False)


def test_group_file_rejects_metadata_without_filename_like_legacy(tmp_path):
    feature_path = tmp_path / "features.csv"
    metadata_path = tmp_path / "metadata.tsv"
    _feature_table().to_csv(feature_path, index=False)
    pd.DataFrame({"ATTRIBUTE_group": ["case"]}).to_csv(
        metadata_path, sep="\t", index=False
    )
    with pytest.raises(Exception, match="filename"):
        group_abundances.calculate_groups_file(
            feature_path, metadata_path, tmp_path / "out.tsv", 1
        )


def test_group_duplicate_metadata_rows_keep_legacy_weighting():
    features = _feature_table()
    metadata = pd.DataFrame(
        {
            "filename": ["alpha.mzML", "beta file.mzML", "beta file.mzML"],
            "ATTRIBUTE_group": ["same", "same", "same"],
        }
    )
    actual = group_abundances.calculate_groups_metadata(features, metadata)
    expected = _legacy_group(features, metadata)
    assert_frame_equal(actual, expected, check_dtype=False)
    # Row 11 is weighted as (alpha + beta + beta) / 3, exactly as merge did.
    assert actual.loc[0, "ATTRIBUTE_group:GNPSGROUP:same"] == pytest.approx(8 / 3)


@pytest.mark.parametrize("chunk_size", [1, 2, 100])
def test_group_file_is_byte_identical_to_legacy_float_reduction(
    tmp_path, chunk_size
):
    features = pd.DataFrame(
        {
            "row ID": [1, 2, 3, 4],
            "row m/z": [100.1, 200.2, 300.3, 400.4],
            "row retention time": [1.1, 2.2, 3.3, 4.4],
            "third.mzML Peak area": [1.0e16, 0.1, np.nan, 7.0],
            "first.mzML Peak area": [1.0, 0.2, np.nan, 11.0],
            "second.mzML Peak area": [-1.0e16, 0.3, np.nan, 13.0],
        }
    )
    # Metadata order intentionally differs from feature-column order, and a
    # duplicate row exercises the legacy merge weighting.
    metadata = pd.DataFrame(
        {
            "filename": [
                "second.mzML",
                "first.mzML",
                "third.mzML",
                "first.mzML",
            ],
            "ATTRIBUTE_group": ["mixed", "mixed", "mixed", "mixed"],
        }
    )
    feature_path = tmp_path / "features.csv"
    metadata_path = tmp_path / "metadata.tsv"
    output_path = tmp_path / "groups.tsv"
    expected_path = tmp_path / "expected.tsv"
    features.to_csv(feature_path, index=False)
    metadata.to_csv(metadata_path, sep="\t", index=False)
    _legacy_group(features, metadata).to_csv(expected_path, sep="\t", index=False)

    group_abundances.calculate_groups_file(
        feature_path, metadata_path, output_path, chunk_size
    )

    assert output_path.read_bytes() == expected_path.read_bytes()


def test_group_file_header_only_input(tmp_path):
    features = _feature_table().iloc[0:0]
    feature_path = tmp_path / "features.csv"
    metadata_path = tmp_path / "metadata.tsv"
    output_path = tmp_path / "groups.tsv"
    features.to_csv(feature_path, index=False)
    _metadata().to_csv(metadata_path, sep="\t", index=False)
    group_abundances.calculate_groups_file(
        feature_path, metadata_path, output_path, chunk_size=1
    )
    expected_path = tmp_path / "expected.tsv"
    _legacy_group(features, _metadata()).to_csv(expected_path, sep="\t", index=False)
    assert output_path.read_bytes() == expected_path.read_bytes()


@pytest.mark.parametrize("style", ["minimal", "rich"])
@pytest.mark.parametrize("raw_mode", ["none", "partial", "all"])
@pytest.mark.parametrize("chunk_size", [1, 2, 100])
def test_tall_file_matches_legacy_for_input_and_raw_variants(
    tmp_path, style, raw_mode, chunk_size
):
    features = _feature_table(style)
    feature_path = tmp_path / "features.csv"
    raw_dir = tmp_path / "raw"
    output_path = tmp_path / "tall.tsv"
    raw_dir.mkdir()
    features.to_csv(feature_path, index=False)
    if raw_mode in {"partial", "all"}:
        (raw_dir / "alpha.mzML").touch()
    if raw_mode == "all":
        (raw_dir / "beta file.mzML").touch()
        (raw_dir / "unmatched.mzML").touch()
    (raw_dir / "not-a-quant-column.mzML").touch()

    create_tall_quant.create_tall_file(
        feature_path, raw_dir, output_path, chunk_size
    )
    actual = pd.read_csv(output_path, sep="\t")
    expected = _legacy_tall(features, raw_dir)
    assert_frame_equal(
        _sort_tall(actual), _sort_tall(expected), check_dtype=False
    )
    assert list(actual.columns) == list(expected.columns)


@pytest.mark.parametrize("chunk_size", [1, 2, 100])
def test_tall_file_is_byte_stable_for_fixed_chunk_size(tmp_path, chunk_size):
    features = _feature_table("rich")
    feature_path = tmp_path / "features.csv"
    raw_dir = tmp_path / "raw"
    first_path = tmp_path / "first.tsv"
    second_path = tmp_path / "second.tsv"
    raw_dir.mkdir()
    (raw_dir / "alpha.mzML").touch()
    features.to_csv(feature_path, index=False)

    create_tall_quant.create_tall_file(
        feature_path, raw_dir, first_path, chunk_size
    )
    create_tall_quant.create_tall_file(
        feature_path, raw_dir, second_path, chunk_size
    )

    assert first_path.read_bytes() == second_path.read_bytes()


def test_tall_file_chunk_major_order_is_explicit(tmp_path):
    features = _feature_table()
    feature_path = tmp_path / "features.csv"
    output_path = tmp_path / "tall.tsv"
    features.to_csv(feature_path, index=False)

    create_tall_quant.create_tall_file(
        feature_path, tmp_path / "NO_FILE", output_path, chunk_size=2
    )
    actual = pd.read_csv(output_path, sep="\t")

    assert list(zip(actual["filename"], actual["row ID"])) == [
        ("alpha.mzML", 11),
        ("alpha.mzML", 12),
        ("beta file.mzML", 11),
        ("beta file.mzML", 12),
        ("unmatched.mzML", 11),
        ("unmatched.mzML", 12),
        ("alpha.mzML", 13),
        ("beta file.mzML", 13),
        ("unmatched.mzML", 13),
    ]


def test_tall_file_missing_raw_path_and_header_only_input(tmp_path):
    features = _feature_table().iloc[0:0]
    feature_path = tmp_path / "features.csv"
    output_path = tmp_path / "tall.tsv"
    features.to_csv(feature_path, index=False)
    create_tall_quant.create_tall_file(
        feature_path, tmp_path / "NO_FILE", output_path, chunk_size=1
    )
    actual = pd.read_csv(output_path, sep="\t")
    assert actual.empty
    assert list(actual.columns) == create_tall_quant.OUTPUT_COLUMNS


def test_tall_file_without_abundance_columns(tmp_path):
    features = _feature_table()[BASE]
    feature_path = tmp_path / "features.csv"
    output_path = tmp_path / "tall.tsv"
    features.to_csv(feature_path, index=False)
    create_tall_quant.create_tall_file(
        feature_path, tmp_path / "NO_FILE", output_path, chunk_size=1
    )
    actual = pd.read_csv(output_path, sep="\t")
    expected = _legacy_tall(features, tmp_path / "NO_FILE")
    assert_frame_equal(actual, expected, check_dtype=False)


@pytest.mark.parametrize("implementation", ["groups", "tall"])
def test_required_columns_are_validated_before_output(tmp_path, implementation):
    feature_path = tmp_path / "features.csv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"row ID": [1], "row m/z": [100.0]}).to_csv(
        feature_path, index=False
    )
    with pytest.raises(ValueError, match="row retention time"):
        if implementation == "groups":
            group_abundances.calculate_groups_file(
                feature_path, tmp_path / "NO_FILE", output_path, 1
            )
        else:
            create_tall_quant.create_tall_file(
                feature_path, tmp_path / "NO_FILE", output_path, 1
            )
    assert not output_path.exists()
    assert not Path(str(output_path) + ".tmp").exists()


@pytest.mark.parametrize(
    "script,extra_args",
    [
        ("group_abundances.py", ["metadata.tsv", "groups.tsv"]),
        ("create_tall_quant.py", ["raw", "tall.tsv"]),
    ],
)
def test_cli_full_path_and_chunk_option(tmp_path, script, extra_args):
    feature_path = tmp_path / "features.csv"
    _feature_table().to_csv(feature_path, index=False)
    _metadata().to_csv(tmp_path / "metadata.tsv", sep="\t", index=False)
    (tmp_path / "raw").mkdir()
    (tmp_path / "raw" / "alpha.mzML").touch()
    command = [
        sys.executable,
        str(SCRIPT_DIR / script),
        str(feature_path),
        *(str(tmp_path / value) for value in extra_args),
        "--chunk-size",
        "1",
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    assert (tmp_path / extra_args[-1]).is_file()


@pytest.mark.parametrize("script", ["group_abundances.py", "create_tall_quant.py"])
def test_cli_rejects_nonpositive_chunk_size(tmp_path, script):
    result = subprocess.run(
        [sys.executable, str(SCRIPT_DIR / script), "a", "b", "c", "--chunk-size", "0"],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "--chunk-size must be at least 1" in result.stderr
