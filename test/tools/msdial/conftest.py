"""MS-DIAL fixtures and subprocess helpers, independent of the current cwd."""
from pathlib import Path
import os
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / 'bin/scripts'
sys.path.insert(0, str(SCRIPTS))
FIXTURES = Path(__file__).parent / 'fixtures'


@pytest.fixture
def scripts_dir():
    return SCRIPTS


@pytest.fixture
def fixture_dir():
    return FIXTURES


@pytest.fixture
def local_data():
    return Path(os.environ.get('GNPS_MSDIAL_TEST_DATA', ROOT / 'data/msdial'))


@pytest.fixture
def dispatch(tmp_path):
    """Run the actual reformat_quantification entry point, not a mock dispatch."""
    def run(tool, feature, spectra, normalization='None'):
        output = tmp_path / tool
        output.mkdir()
        spectrum_dir = output / 'input_spectra'
        spectrum_dir.mkdir()
        source = spectrum_dir / spectra.name
        source.write_bytes(spectra.read_bytes())
        feature_out = output / 'features.csv'
        mgf_out = output / 'spectra.mgf'
        result = subprocess.run(
            [sys.executable, str(SCRIPTS / 'reformat_quantification.py'),
             tool, str(feature), str(feature_out), str(spectrum_dir),
             str(mgf_out), str(output / 'embedded_metadata.tsv'),
             '--QUANT_FILE_NORM', normalization],
            cwd=output, capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return feature_out, mgf_out, source
    return run
