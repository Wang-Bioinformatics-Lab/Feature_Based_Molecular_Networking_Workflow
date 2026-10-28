import csv
from pathlib import Path
import subprocess
import sys

import pytest

from validate_feature_spectrum_alignment import AlignmentSkewError, validate_alignment


SCRIPT = Path(__file__).with_name("validate_feature_spectrum_alignment.py")


def inputs(tmp_path, rows, spectra):
    feature = tmp_path / "features.csv"
    with feature.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row ID", "row m/z", "row retention time"])
        writer.writerows(rows)
    mgf = tmp_path / "spectra.mgf"
    mgf.write_text("".join(
        "BEGIN IONS\nSCANS={}\nPEPMASS={}\n{}100 10\nEND IONS\n".format(identifier, mz, rt)
        for identifier, mz, rt in spectra
    ))
    return feature, mgf


@pytest.mark.parametrize("feature_rt,unit,mgf_rt", [
    (2.0, "minutes", "RTINSECONDS=132\n"),
    (120.0, "seconds", "RTINMINUTES=2.2\n"),
    (2.0, "minutes", "RTINMINUTES=1.8\n"),
])
def test_loose_tolerance_boundaries_and_rt_units(tmp_path, feature_rt, unit, mgf_rt):
    # Out-of-order IDs, integer-valued CSV IDs, and MS1-only feature rows are valid.
    feature, mgf = inputs(tmp_path, [(7, 300, 8), (12.0, 200, feature_rt)],
                          [(12, 201, mgf_rt)])
    assert validate_alignment(feature, mgf, unit) == (1, 1)


@pytest.mark.parametrize("mz,rt,reason", [
    (201.01, "RTINMINUTES=2\n", "precursor m/z"),
    (198.99, "RTINMINUTES=2\n", "precursor m/z"),
    (200, "RTINSECONDS=132.6\n", "spectrum RT"),
    (200, "RTINMINUTES=1.79\n", "spectrum RT"),
])
def test_same_ids_with_wrong_measurements_fail(tmp_path, mz, rt, reason):
    feature, mgf = inputs(tmp_path, [(33, 200, 2)], [(33, mz, rt)])
    with pytest.raises(AlignmentSkewError, match=reason):
        validate_alignment(feature, mgf, "minutes")


def test_missing_spectrum_rt_still_checks_precursor(tmp_path):
    feature, mgf = inputs(tmp_path, [(33, 150.12791, 2.839)], [(33, 154.990494, "")])
    with pytest.raises(AlignmentSkewError, match="precursor m/z"):
        validate_alignment(feature, mgf, "minutes")


def test_missing_rt_and_unknown_feature_rt_units_are_allowed(tmp_path):
    feature, mgf = inputs(tmp_path, [(33, 200, 120)], [(33, 200.9, "")])
    assert validate_alignment(feature, mgf, "seconds") == (1, 0)
    mgf.write_text(mgf.read_text().replace("PEPMASS=200.9", "PEPMASS=200.9 1000\nRTINMINUTES=2"))
    assert validate_alignment(feature, mgf) == (1, 0)


def test_cli_fails_with_readable_persisted_summary_and_clears_it_on_success(tmp_path):
    feature, mgf = inputs(tmp_path, [(33, 200, 2)], [(1, 200, "")])
    summary = tmp_path / "published/nf_output/workflow_summary.txt"
    command = [sys.executable, str(SCRIPT), str(feature), str(mgf),
               "--summary", str(summary), "--featurefindingtool", "MSDIAL5",
               "--feature-source", "original features.txt", "--spectra-source", "original spectra.mgf"]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 1
    message = summary.read_text()
    assert 'An ID alignment skew has been detected between feature table "original features.txt" and MS/MS spectra "original spectra.mgf".' in message
    assert "If you believe this is not a result of unrelated input data, please contact the GNPS2 developers." in message
    assert "Spectrum ID 1 has no matching feature-table row ID" in message
    assert len(message) < 1000
    assert message in result.stderr
    mgf.write_text(mgf.read_text().replace("SCANS=1\n", "SCANS=33\n"))
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0
    assert not summary.exists()
    summary.write_text("An unrelated workflow summary.\n")
    assert subprocess.run(command, capture_output=True).returncode == 0
    assert summary.read_text() == "An unrelated workflow summary.\n"


def test_repository_mzmine_example_passes():
    folder = Path(__file__).resolve().parents[2] / "data/mzmine2/gnps_featurefinding"
    checked, rt_checked = validate_alignment(
        folder / "features_quant.csv", folder / "spectra/features.mgf", "minutes"
    )
    assert checked == rt_checked == 2594
