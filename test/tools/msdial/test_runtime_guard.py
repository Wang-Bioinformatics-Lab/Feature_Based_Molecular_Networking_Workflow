"""Known-bad historical outputs must be rejected, even after conversion is fixed."""
import subprocess
import sys
import textwrap

import pytest

from validate_feature_spectrum_alignment import AlignmentSkewError, mgf_headers, validate_alignment


@pytest.mark.parametrize('case', ['msdial5_interlab_original', 'msdial5_interlab_extensionless_metadata'])
def test_known_bad_published_outputs_are_rejected(local_data, case):
    folder = local_data / case / 'observed_outputs/nf_output/clustering'
    if not folder.exists():
        pytest.skip('Download fixtures: python data/download_msdial_test_data.py')
    with pytest.raises(AlignmentSkewError, match='Spectrum ID 1 has no matching feature-table row ID'):
        validate_alignment(folder / 'featuretable_reformated.csv', folder / 'spectra_reformatted.mgf')


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


@pytest.mark.parametrize('tool', [
    'MSDIAL5', 'MSDIAL', 'MZMINE', 'EVERYTHING_BAGEL', 'ASARI', 'OPENMS',
    'XCMS3', 'MASSCUBE', 'PROGENESIS', 'METABOSCAPE', 'OPTIMUS',
    'AGILENT', 'AGILENT2', 'AGILENT_EXPLORER2', 'MZTABM',
])
def test_workflow_enforces_guard_only_for_msdial5(tmp_path, scripts_dir, tool):
    # Execute the real shell gate from the Nextflow process, with skewed
    # converted inputs. Other exporters have not established the same ID and
    # padding contract, so this PR must not introduce failures for those tools.
    workflow = (scripts_dir.parents[1] / 'nf_workflow.nf').read_text()
    start = workflow.index('    if [ "$params.featurefindingtool" = "MSDIAL5" ]; then')
    end = workflow.index('\n    fi', start) + len('\n    fi')
    command = textwrap.dedent(workflow[start:end])
    summary = tmp_path / 'nf_output/workflow_summary.txt'
    command = command.replace('${file(params.publishdir).toAbsolutePath()}', str(tmp_path))
    command = command.replace('$params.featurefindingtool', tool)
    command = command.replace('$TOOL_FOLDER', str(scripts_dir.parent))
    command = command.replace('$input_features', 'source_features.txt')
    command = command.replace('$input_spectra', 'source_spectra.mgf')
    # Use this test's interpreter without relying on the shell's Python setup.
    command = command.replace('python ', '"{}" '.format(sys.executable))
    (tmp_path / 'featuretable_reformated.csv').write_text(
        'row ID,row m/z,row retention time\n33,200,2\n'
    )
    (tmp_path / 'specs_ms.mgf').write_text(
        'BEGIN IONS\nSCANS=1\nPEPMASS=200\nRTINSECONDS=120\n100 10\nEND IONS\n'
    )
    result = subprocess.run(['bash', '-e', '-c', command], cwd=tmp_path,
                            capture_output=True, text=True)
    if tool == 'MSDIAL5':
        assert result.returncode == 1
        assert 'An ID alignment skew has been detected' in summary.read_text()
    else:
        assert result.returncode == 0, result.stderr
        assert not summary.exists()
