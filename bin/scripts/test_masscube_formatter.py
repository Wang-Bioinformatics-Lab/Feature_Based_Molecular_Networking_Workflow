import os
import re
import pandas as pd
import masscube_formatter

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data", "MassCube_data")
FEATURE_TABLE = os.path.join(DATA_DIR, "aligned_feature_table_before_annotation.txt")
MSP = os.path.join(DATA_DIR, "features.msp")


def test_feature_csv_headers_and_rows(tmp_path):
    out = str(tmp_path / "ft.csv")
    masscube_formatter.convert_to_feature_csv(FEATURE_TABLE, out)
    df = pd.read_csv(out)

    # First three columns are the GNPS-standard ones, in order
    assert list(df.columns[:3]) == ["row ID", "row m/z", "row retention time"]

    # Every remaining column is a per-sample "Peak area" column (13 samples)
    peak_cols = [c for c in df.columns if c.endswith(" Peak area")]
    assert len(peak_cols) == 13
    assert len(df.columns) == 3 + 13

    # All features are kept and row IDs are the feature IDs 1..29851
    assert len(df) == 29851
    assert df["row ID"].min() == 1
    assert df["row ID"].max() == 29851
    assert df["row ID"].is_unique

    # Spot-check feature 1: m/z 306.207, RT 8.37
    row1 = df[df["row ID"] == 1].iloc[0]
    assert abs(row1["row m/z"] - 306.207) < 1e-3
    assert abs(row1["row retention time"] - 8.37) < 1e-3


def test_mgf_only_ms2_features(tmp_path):
    out = str(tmp_path / "specs.mgf")
    written = masscube_formatter.convert_mgf(MSP, out)

    # Only MS2-bearing features (Num Peaks > 0) become MGF blocks
    assert written == 3624

    content = open(out).read()
    assert content.count("BEGIN IONS") == 3624
    assert content.count("END IONS") == 3624

    # No empty spectra: every block carries at least one peak line
    blocks = content.split("BEGIN IONS")[1:]
    for block in blocks:
        body = block.split("END IONS")[0]
        peak_lines = [ln for ln in body.strip().splitlines()
                      if ln and not ln.startswith(("TITLE", "SCANS", "PEPMASS"))]
        assert len(peak_lines) >= 1


def test_mgf_scan_pepmass_and_linkage(tmp_path):
    mgf_out = str(tmp_path / "specs.mgf")
    csv_out = str(tmp_path / "ft.csv")
    masscube_formatter.convert_mgf(MSP, mgf_out)
    masscube_formatter.convert_to_feature_csv(FEATURE_TABLE, csv_out)

    content = open(mgf_out).read()

    # First feature: SCANS=1 with PEPMASS 306.207
    first_block = content.split("END IONS")[0]
    assert "SCANS=1" in first_block
    assert "TITLE=SCAN=1" in first_block
    assert "PEPMASS=306.207" in first_block

    # Every SCANS in the MGF exists as a row ID in the quant table
    scans = set(int(s) for s in re.findall(r"SCANS=(\d+)", content))
    row_ids = set(pd.read_csv(csv_out)["row ID"].tolist())
    assert scans.issubset(row_ids)
    assert len(scans) == 3624
