"""Expected MSDIAL5 behavior; these tests expose the original conversion bugs."""
import pandas as pd
import pytest
from pyteomics import mgf

import msdial5_formatter
from validate_feature_spectrum_alignment import validate_alignment


@pytest.fixture
def converted(tmp_path, fixture_dir, monkeypatch):
    # The old converter uses a cwd-local temp.mgf; keep tests isolated.
    monkeypatch.chdir(tmp_path)
    feature = tmp_path / 'features.csv'
    spectra = tmp_path / 'spectra.mgf'
    msdial5_formatter.convert_to_feature_csv(str(fixture_dir / 'v5_features.txt'), str(feature))
    msdial5_formatter.convert_mgf(str(fixture_dir / 'spectra.mgf'), str(spectra))
    return feature, spectra


def test_feature_ids_preserve_original_sparse_and_zero_alignment_ids(converted):
    feature, _ = converted
    assert pd.read_csv(feature)['row ID'].tolist() == [32, 0, 55, 99]


def test_mgf_ids_preserve_original_ids_despite_different_spectrum_order(converted):
    _, spectra = converted
    with mgf.read(str(spectra)) as reader:
        assert [int(s['params']['scans']) for s in reader] == [55, 32, 0]


def test_only_real_sample_columns_become_peak_areas(converted):
    feature, _ = converted
    table = pd.read_csv(feature)
    assert table.columns.tolist() == [
        'row ID', 'row m/z', 'row retention time',
        'sample_A Peak area', 'sample_B Peak area',
    ]
    assert table[['sample_A Peak area', 'sample_B Peak area']].values.tolist() == [
        [10, 20], [30, 40], [50, 60], [70, 80],
    ]


def test_precursor_rt_and_charge_remain_attached_to_the_original_id(converted):
    _, spectra = converted
    expected = {55: (451.0354, 7.617), 32: (150.12791, 2.839), 0: (300.12345, 4.25)}
    with mgf.read(str(spectra)) as reader:
        for spectrum in reader:
            params = spectrum['params']
            identifier = int(params['scans'])
            assert identifier in expected
            mz, rt = expected[identifier]
            assert params['pepmass'][0] == pytest.approx(mz, abs=1e-6)
            assert params['charge'][0] == 1
            # Read either legal MGF RT unit; do not require a serialization style.
            actual_rt = (float(params['rtinseconds']) / 60
                         if 'rtinseconds' in params else float(params['rtinminutes']))
            assert actual_rt == pytest.approx(rt, abs=1e-6)


def test_existing_one_percent_fragment_filter_stays_unchanged(converted):
    _, spectra = converted
    with mgf.read(str(spectra)) as reader:
        for spectrum in reader:
            assert spectrum['m/z array'].tolist() == [100.0, 110.0]
            assert spectrum['intensity array'].tolist() == [1000.0, 500.0]


def test_converted_outputs_pass_runtime_alignment_guard(converted):
    feature, spectra = converted
    assert validate_alignment(feature, spectra, 'minutes') == (3, 3)


def test_v5_dispatch_preserves_alignment_and_sample_values(fixture_dir, dispatch):
    feature, spectra, _ = dispatch(
        'MSDIAL5', fixture_dir / 'v5_features.txt', fixture_dir / 'spectra.mgf'
    )
    table = pd.read_csv(feature)
    assert table['row ID'].tolist() == [32, 0, 55, 99]
    assert table.columns.tolist() == [
        'row ID', 'row m/z', 'row retention time',
        'sample_A Peak area', 'sample_B Peak area',
    ]
    assert validate_alignment(feature, spectra, 'minutes') == (3, 3)


@pytest.mark.parametrize('case', ['msdial5_interlab_original', 'msdial5_interlab_extensionless_metadata'])
def test_real_v5_inputs_preserve_ids_and_align_after_conversion(tmp_path, local_data, case, monkeypatch):
    folder = local_data / case
    if not folder.exists():
        pytest.skip('Download fixtures: python data/download_msdial_test_data.py')
    monkeypatch.chdir(tmp_path)
    source = next((folder / 'inputfeatures').iterdir())
    source_table = pd.read_csv(source, sep='\t', skiprows=4)
    feature = tmp_path / 'features.csv'
    spectra = tmp_path / 'spectra.mgf'
    msdial5_formatter.convert_to_feature_csv(str(source), str(feature))
    msdial5_formatter.convert_mgf(str(next((folder / 'inputspectra').iterdir())), str(spectra))
    result = pd.read_csv(feature)
    assert result['row ID'].tolist() == source_table['Alignment ID'].tolist()
    assert result.columns.tolist() == [
        'row ID', 'row m/z', 'row retention time',
        *[c + ' Peak area' for c in source_table.columns[35:]],
    ]
    assert validate_alignment(feature, spectra, 'minutes') == (4653, 4653)


def test_compact_v5_export_retains_existing_sample_column_support(tmp_path, fixture_dir):
    output = tmp_path / 'compact.csv'
    msdial5_formatter.convert_to_feature_csv(str(fixture_dir / 'v5_compact_features.txt'), str(output))
    table = pd.read_csv(output)
    assert table.columns.tolist() == [
        'row ID', 'row m/z', 'row retention time',
        'sample_A Peak area', 'sample_B Peak area',
    ]
    assert table[['sample_A Peak area', 'sample_B Peak area']].values.tolist() == [
        [10, 20], [30, 40], [50, 60], [70, 80],
    ]
