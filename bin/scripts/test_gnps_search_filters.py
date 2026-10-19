"""Exercise filter switches through the wrapper and bundled GNPS search engine."""

import csv
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
CLEAN = [(100 + i, 100 - i) for i in range(6)] + [(250, 80), (350, 70)]


def write_spectrum(path, peaks):
    path.write_text(
        "BEGIN IONS\nPEPMASS=1000\nCHARGE=1+\nSCANS=1\n"
        "SPECTRUMID=FILTER_TEST\nNAME=filter_test\n"
        + "".join(f"{mz} {intensity}\n" for mz, intensity in peaks)
        + "END IONS\n"
    )


def search(tmp_path, contaminated_side, switches):
    query, library = tmp_path / "query.mgf", tmp_path / "library.mgf"
    contaminated = CLEAN + [(106, 90), (995, 10000)]
    write_spectrum(query, contaminated if contaminated_side == "query" else CLEAN)
    write_spectrum(library, contaminated if contaminated_side == "library" else CLEAN)
    (tmp_path / "results").mkdir()
    command = [
        sys.executable, str(ROOT / "bin/scripts/library_search_wrapper.py"),
        str(query), str(library), "results", str(ROOT / "bin/binaries/convert"),
        str(ROOT / "bin/binaries/main_execmodule.allcandidates"),
        "--fragment_tolerance", "0.01", "--library_min_cosine", "0",
        "--library_min_matched_peaks", "1",
    ]
    if switches is not None:
        precursor, window = switches
        command += ["--filter_precursor", str(precursor), "--filter_window", str(window)]
    result = subprocess.run(
        command, cwd=tmp_path, env=dict(os.environ, OMP_NUM_THREADS="1"),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    assert result.returncode == 0, result.stdout
    with (tmp_path / "results/query.mgf_library.mgf.tsv").open() as handle:
        hits = list(csv.DictReader(handle, delimiter="\t"))
    assert len(hits) == 1, result.stdout
    return float(hits[0]["MQScore"])


@pytest.mark.parametrize("contaminated_side", ["query", "library"])
@pytest.mark.parametrize("precursor,window", [(0, 0), (0, 1), (1, 0), (1, 1)])
def test_filters_on_each_side(tmp_path, contaminated_side, precursor, window):
    score = search(tmp_path, contaminated_side, (precursor, window))
    # GNPS square-roots intensities before cosine normalization. The two extra
    # peaks cannot match; only the selected filters should remove their signal.
    clean_intensity = sum(intensity for _, intensity in CLEAN)
    remaining_extra = (0 if precursor else 10000) + (0 if window else 90)
    expected = math.sqrt(clean_intensity / (clean_intensity + remaining_extra))
    assert score == pytest.approx(expected, abs=1e-5)


def test_defaults_keep_both_filters_enabled(tmp_path):
    assert search(tmp_path, "library", None) == pytest.approx(1, abs=1e-5)
