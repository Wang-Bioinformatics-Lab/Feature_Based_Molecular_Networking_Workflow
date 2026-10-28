"""Regression contracts for the existing MSDIAL (pre-v5) conversion branch."""
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import msdial_formatter
from validate_feature_spectrum_alignment import validate_alignment


def test_legacy_dispatch_preserves_shared_ids_acquisition_and_mgf(fixture_dir, dispatch):
    feature_out, mgf_out, source = dispatch(
        'MSDIAL', fixture_dir / 'legacy_features.txt', fixture_dir / 'spectra.mgf'
    )
    table = pd.read_csv(feature_out)
    assert table['row ID'].tolist() == [32, 0, 55, 99]
    assert table['row m/z'].tolist() == [150.12791, 300.12345, 451.0354, 777.0]
    assert table['row retention time'].tolist() == [2.839, 4.25, 7.617, 9.0]
    assert table.columns.tolist() == [
        'row ID', 'row m/z', 'row retention time',
        'sample_A Peak area', 'sample_B Peak area',
    ]
    assert table[['sample_A Peak area', 'sample_B Peak area']].values.tolist() == [
        [10, 20], [30, 40], [50, 60], [70, 80],
    ]
    # The legacy path intentionally copies source MGF bytes, retaining charge,
    # source RT, sparse IDs, zero, and an order different from the feature CSV.
    assert mgf_out.read_bytes() == source.read_bytes()
    assert validate_alignment(feature_out, mgf_out, 'minutes') == (3, 3)


def test_legacy_pasef_columns_do_not_become_sample_abundances(tmp_path, fixture_dir):
    source = pd.read_csv(fixture_dir / 'legacy_features.txt', sep='\t', skiprows=3)
    source.insert(3, 'Average drift time', [1.1, 1.2, 1.3, 1.4])
    source.insert(4, 'Average CCS', [150, 160, 170, 180])
    input_path = tmp_path / 'pasef.txt'
    input_path.write_text('\n\n\n' + source.to_csv(sep='\t', index=False))
    output = tmp_path / 'features.csv'
    msdial_formatter.convert_to_feature_csv(str(input_path), str(output))
    result = pd.read_csv(output)
    assert result.columns.tolist() == [
        'row ID', 'row m/z', 'row retention time',
        'sample_A Peak area', 'sample_B Peak area',
    ]
    assert result['row ID'].tolist() == [32, 0, 55, 99]


@pytest.mark.parametrize('case', ['legacy_american_gut_2019', 'legacy_pasef_2019'])
def test_historical_inputs_reproduce_reference_output_and_align(tmp_path, local_data, case):
    folder = local_data / case
    if not folder.exists():
        pytest.skip('Download fixtures: python data/download_msdial_test_data.py')
    output = tmp_path / 'features.csv'
    msdial_formatter.convert_to_feature_csv(
        str(next((folder / 'inputfeatures').iterdir())), str(output)
    )
    assert_frame_equal(
        pd.read_csv(output), pd.read_csv(next((folder / 'reference_outputs').iterdir())),
        check_dtype=False, rtol=1e-10, atol=1e-8,
    )
    expected = 8689 if case == 'legacy_american_gut_2019' else 842
    assert validate_alignment(output, next((folder / 'inputspectra').iterdir()), 'minutes') == (expected, expected)
