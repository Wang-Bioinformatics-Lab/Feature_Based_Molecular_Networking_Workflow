"""Known-bad historical outputs must be rejected, even after conversion is fixed."""
import subprocess
import sys

import pytest

from validate_feature_spectrum_alignment import AlignmentSkewError, mgf_headers, validate_alignment


@pytest.mark.parametrize('case', ['msdial5_interlab_original', 'msdial5_interlab_extensionless_metadata'])
def test_known_bad_published_outputs_are_rejected(local_data, case):
    folder = local_data / case / 'observed_outputs/nf_output/clustering'
    if not folder.exists():
        pytest.skip('Download fixtures: python data/download_msdial_test_data.py')
    with pytest.raises(AlignmentSkewError, match='Spectrum ID 1 has no matching feature-table row ID'):
        validate_alignment(folder / 'featuretable_reformated.csv', folder / 'spectra_reformatted.mgf', 'minutes')


def test_filtering_preserves_original_zero_id_and_its_precursor(tmp_path, fixture_dir, scripts_dir):
    output = tmp_path / 'spectra_reformatted.mgf'
    filtered = tmp_path / 'spectra_filtered.mgf'
    result = subprocess.run(
        [sys.executable, str(scripts_dir / 'filter_spectra.py'), str(fixture_dir / 'spectra.mgf'),
         str(output), str(filtered), '--FILTER_PRECURSOR_WINDOW', '0', '--WINDOW_FILTER', '0'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for path in [output, filtered]:
        actual = {int(s['SCANS']): float(s['PEPMASS']) for s in mgf_headers(path) if float(s['PEPMASS']) > 0}
        assert actual == {0: pytest.approx(300.12345), 32: pytest.approx(150.12791), 55: pytest.approx(451.0354)}
