"""Execute the production shell gate to test importer enforcement scope."""
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest


ROOT = Path(__file__).resolve().parents[3]
ENFORCED = {"MSDIAL5", "MSDIAL", "MZMINE", "EVERYTHING_BAGEL", "ASARI", "MASSCUBE", "OPENMS", "XCMS3"}


@pytest.mark.parametrize("tool", [
    "MSDIAL5", "MSDIAL", "MZMINE", "EVERYTHING_BAGEL", "ASARI", "OPENMS",
    "XCMS3", "MASSCUBE", "PROGENESIS", "METABOSCAPE", "OPTIMUS",
    "AGILENT", "AGILENT2", "AGILENT_EXPLORER2", "MZTABM",
])
def test_workflow_enforces_only_validated_importer_paths(tmp_path, tool):
    workflow = (ROOT / "nf_workflow.nf").read_text()
    start = workflow.index('    case "$params.featurefindingtool" in')
    end = workflow.index("\n    esac", start) + len("\n    esac")
    command = textwrap.dedent(workflow[start:end])
    command = command.replace("${file(params.publishdir).toAbsolutePath()}", str(tmp_path))
    command = command.replace("$params.featurefindingtool", tool)
    command = command.replace("$TOOL_FOLDER", str(ROOT / "bin"))
    command = command.replace("$input_features", "source_features.txt")
    command = command.replace("$input_spectra", "source_spectra.mgf")
    command = command.replace("python ", '"{}" '.format(sys.executable))
    (tmp_path / "featuretable_reformated.csv").write_text("row ID,row m/z\n33,200\n")
    (tmp_path / "specs_ms.mgf").write_text("BEGIN IONS\nSCANS=1\nPEPMASS=200\n100 10\nEND IONS\n")
    result = subprocess.run(["bash", "-e", "-c", command], cwd=tmp_path,
                            capture_output=True, text=True)
    summary = tmp_path / "nf_output/workflow_summary.txt"
    if tool in ENFORCED:
        assert result.returncode == 1
        assert "An ID alignment skew has been detected" in summary.read_text()
    else:
        assert result.returncode == 0, result.stderr
        assert not summary.exists()
