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


@pytest.mark.parametrize("mz", [199, 201])
def test_loose_precursor_tolerance_boundaries(tmp_path, mz):
    # Out-of-order IDs, integer-valued CSV IDs, and MS1-only rows are valid.
    feature, mgf = inputs(tmp_path, [(7, 300, 8), (12.0, 200, 2)], [(12, mz, "")])
    assert validate_alignment(feature, mgf) == 1


@pytest.mark.parametrize("mz", [201.01, 198.99])
def test_same_ids_with_wrong_precursors_fail(tmp_path, mz):
    feature, mgf = inputs(tmp_path, [(33, 200, 2)], [(33, mz, "")])
    with pytest.raises(AlignmentSkewError, match="precursor m/z"):
        validate_alignment(feature, mgf)


@pytest.mark.parametrize("rt", ["", "RTINSECONDS=1800\n", "RTINMINUTES=30\n"])
def test_rt_is_not_used_to_reject_an_aligned_precursor(tmp_path, rt):
    feature, mgf = inputs(tmp_path, [(33, 200, 2)], [(33, 200.9, rt)])
    assert validate_alignment(feature, mgf) == 1


def test_cli_fails_with_readable_persisted_summary_and_clears_it_on_success(tmp_path):
    feature, mgf = inputs(tmp_path, [(33, 200, 2)], [(1, 200, "")])
    summary = tmp_path / "published/nf_output/workflow_summary.txt"
    command = [sys.executable, str(SCRIPT), str(feature), str(mgf),
               "--summary", str(summary),
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


def test_cli_ignores_rt_but_still_rejects_a_wrong_precursor(tmp_path):
    feature, mgf = inputs(tmp_path, [(33, 200, 2)], [(33, 200.5, "RTINSECONDS=1800\n")])
    summary = tmp_path / "nf_output/workflow_summary.txt"
    command = [sys.executable, str(SCRIPT), str(feature), str(mgf), "--summary", str(summary)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert not summary.exists()
    mgf.write_text(mgf.read_text().replace("PEPMASS=200.5", "PEPMASS=201.1"))
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 1
    assert "precursor m/z" in summary.read_text()


def test_repository_mzmine_example_passes():
    folder = Path(__file__).resolve().parents[2] / "data/mzmine2/gnps_featurefinding"
    checked = validate_alignment(
        folder / "features_quant.csv", folder / "spectra/features.mgf"
    )
    assert checked == 2594
