"""Offline precursor/ID contracts from actual feature-finder task excerpts."""
from pathlib import Path
import os
import re

import pytest

import everything_bagel_formatter
from validate_feature_spectrum_alignment import AlignmentSkewError, validate_alignment


FIXTURES = Path(__file__).parent / "fixtures"
TOOLS = ["asari", "everything_bagel", "masscube", "mzmine", "openms", "xcms"]


def converted_mgf(tool, tmp_path):
    source = FIXTURES / (tool + ".mgf")
    if tool != "everything_bagel":
        return source
    output = tmp_path / "eb.mgf"
    assert everything_bagel_formatter.convert_mgf(source, output) == (2, 1)
    return output


@pytest.mark.parametrize("tool", TOOLS)
def test_actual_export_excerpts_have_matching_ids_and_precursors(tmp_path, tool):
    assert validate_alignment(FIXTURES / (tool + ".csv"), converted_mgf(tool, tmp_path)) == 2


@pytest.mark.parametrize("tool", TOOLS)
def test_each_importer_rejects_a_shifted_spectrum_id(tmp_path, tool):
    source = converted_mgf(tool, tmp_path)
    bad = tmp_path / "shifted.mgf"
    bad.write_text(re.sub(r"(?m)^SCANS=.*$", "SCANS=999999", source.read_text(), count=1))
    with pytest.raises(AlignmentSkewError, match="no matching feature-table row ID"):
        validate_alignment(FIXTURES / (tool + ".csv"), bad)


@pytest.mark.parametrize("tool", TOOLS)
def test_each_importer_rejects_a_wrong_precursor(tmp_path, tool):
    source = converted_mgf(tool, tmp_path)
    bad = tmp_path / "wrong_precursor.mgf"
    bad.write_text(re.sub(r"(?m)^PEPMASS=.*$", "PEPMASS=9999", source.read_text(), count=1))
    with pytest.raises(AlignmentSkewError, match="precursor m/z"):
        validate_alignment(FIXTURES / (tool + ".csv"), bad)


@pytest.mark.parametrize("tool,count", [
    ("asari", 1657), ("everything_bagel", 2679), ("masscube", 3864),
    ("mzmine", 2565), ("openms", 2371), ("xcms", 3031), ("msdial", None),
])
def test_full_local_task_outputs(tmp_path, tool, count):
    root = Path(os.environ.get("GNPS_FEATURE_FINDER_TEST_DATA",
                               Path(__file__).resolve().parents[3] / "data/feature_finder_validation"))
    folder = root / tool / "nf_output/clustering"
    if not folder.exists():
        pytest.skip("Optional task files are not present locally; tests never download them")
    feature = folder / "featuretable_reformated.csv"
    spectra = folder / "specs_ms.mgf"
    if tool == "everything_bagel":
        clean = tmp_path / "eb.mgf"
        assert everything_bagel_formatter.convert_mgf(spectra, clean) == (2679, 3308)
        spectra = clean
    if tool == "msdial":
        # The original archived task is intentionally broken, not a golden output.
        with pytest.raises(AlignmentSkewError, match="no matching feature-table row ID"):
            validate_alignment(feature, spectra)
    else:
        assert validate_alignment(feature, spectra) == count
