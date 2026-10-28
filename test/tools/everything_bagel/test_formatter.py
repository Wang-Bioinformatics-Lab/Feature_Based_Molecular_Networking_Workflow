"""EB normalization retains real spectra and original IDs, including rescue MS2."""
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

import everything_bagel_formatter
from validate_feature_spectrum_alignment import AlignmentSkewError, mgf_headers, validate_alignment


ROOT = Path(__file__).resolve().parents[3]
FIXTURES = Path(__file__).parents[1] / "alignment/fixtures"


def test_removes_only_padding_and_preserves_real_blocks_verbatim(tmp_path):
    source = FIXTURES / "everything_bagel.mgf"
    output = tmp_path / "spectra.mgf"
    assert everything_bagel_formatter.convert_mgf(source, output) == (2, 1)
    source_blocks = source.read_text().split("BEGIN IONS")[1:]
    expected = "".join("BEGIN IONS" + block for block in source_blocks if "PEPMASS=0.0\n" not in block)
    assert output.read_text() == expected
    headers = list(mgf_headers(output))
    assert [int(h["SCANS"]) for h in headers] == [1, 3607]
    assert headers[1]["SOURCE_FEATURE_ID"] == "-1"
    assert validate_alignment(FIXTURES / "everything_bagel.csv", output) == 2


@pytest.mark.parametrize("header,peaks", [
    ("PEPMASS=0\nFEATURE_ID=0\n", ""),
    ("PEPMASS=0\n", "100 10\n"),
    ("PEPMASS=200\n", ""),
    ("PEPMASS=nan\n", ""),
    ("", ""),
    ("PEPMASS=0\nTITLE=actual_source_spectrum\n", ""),
])
def test_does_not_hide_real_or_invalid_spectra(tmp_path, header, peaks):
    source = tmp_path / "source.mgf"
    source.write_text("BEGIN IONS\nSCANS=0\n" + header + peaks + "END IONS\n")
    output = tmp_path / "output.mgf"
    assert everything_bagel_formatter.convert_mgf(source, output) == (1, 0)
    assert output.read_bytes() == source.read_bytes()


def test_sparse_and_zero_ids_are_never_renumbered(tmp_path):
    source = tmp_path / "source.mgf"
    source.write_text("".join(
        "BEGIN IONS\nSCANS={}\nPEPMASS=200\n100 10\nEND IONS\n".format(i)
        for i in [55, 0, 32]
    ))
    output = tmp_path / "output.mgf"
    assert everything_bagel_formatter.convert_mgf(source, output) == (3, 0)
    assert output.read_bytes() == source.read_bytes()


def test_duplicate_padding_id_cannot_be_hidden(tmp_path):
    source = tmp_path / "source.mgf"
    source.write_text("BEGIN IONS\nSCANS=11\nPEPMASS=0\nEND IONS\n" * 2)
    with pytest.raises(ValueError, match="Duplicate Everything Bagel spectrum ID"):
        everything_bagel_formatter.convert_mgf(source, tmp_path / "output.mgf")


@pytest.mark.parametrize("text", [
    "BEGIN IONS\nSCANS=1\nPEPMASS=200\n", "BEGIN IONS\nBEGIN IONS\nEND IONS\n",
    "BEGIN IONS\nPEPMASS=200\nEND IONS\n",
])
def test_truncated_nested_or_unidentified_blocks_fail(tmp_path, text):
    source = tmp_path / "source.mgf"
    source.write_text(text)
    with pytest.raises(ValueError):
        everything_bagel_formatter.convert_mgf(source, tmp_path / "output.mgf")


@pytest.mark.parametrize("tool", ["EVERYTHING_BAGEL", "MZMINE"])
def test_dispatch_uses_eb_normalization_without_changing_mzmine(tmp_path, tool):
    spectrum_dir = tmp_path / "input_spectra"
    spectrum_dir.mkdir()
    source = FIXTURES / "everything_bagel.mgf"
    (spectrum_dir / "input.mgf").write_bytes(source.read_bytes())
    feature, spectra = tmp_path / "features.csv", tmp_path / "spectra.mgf"
    result = subprocess.run([
        sys.executable, str(ROOT / "bin/scripts/reformat_quantification.py"), tool,
        str(FIXTURES / "everything_bagel.csv"), str(feature), str(spectrum_dir),
        str(spectra), str(tmp_path / "metadata.tsv"),
    ], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert pd.read_csv(feature)["row ID"].tolist() == [1, 3607]
    if tool == "EVERYTHING_BAGEL":
        assert validate_alignment(feature, spectra) == 2
        assert len(list(mgf_headers(spectra))) == 2
    else:
        assert spectra.read_bytes() == source.read_bytes()


def test_normalization_preserves_existing_filtered_outputs(tmp_path):
    source = FIXTURES / "everything_bagel.mgf"
    compact = tmp_path / "compact.mgf"
    everything_bagel_formatter.convert_mgf(source, compact)
    outputs = []
    for name, spectra in [("original", source), ("compact", compact)]:
        unfiltered, filtered = tmp_path / (name + ".mgf"), tmp_path / (name + "_filtered.mgf")
        result = subprocess.run([
            sys.executable, str(ROOT / "bin/scripts/filter_spectra.py"),
            str(spectra), str(unfiltered), str(filtered),
            "--FILTER_PRECURSOR_WINDOW", "0", "--WINDOW_FILTER", "0",
        ], cwd=tmp_path, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        outputs.append((unfiltered.read_bytes(), filtered.read_bytes()))
    assert outputs[0] == outputs[1]
