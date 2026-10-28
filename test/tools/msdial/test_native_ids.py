"""Keep original MS-DIAL IDs across the positional IDs used by native tools.

Everything is local, including the tiny spectra and the bundled executables.
The mapping helper is imported inside tests so its absence is an ordinary
regression-test failure on the pre-fix revision, rather than a collection error.
"""
import csv
import importlib
import os
from pathlib import Path
import platform
import subprocess
import sys

import pytest

from validate_feature_spectrum_alignment import mgf_headers


PEAKS = '50 1000\n60 900\n70 800\n80 700\n90 600\n100 500\n'


def write_mgf(path, ids, empty_ids=()):
    path.write_text(''.join(
        'BEGIN IONS\nSCANS={}\nTITLE=source_{}\nPEPMASS={}\nCHARGE=1+\n'
        'RTINSECONDS=120\n{}END IONS\n\n'.format(
            identifier, identifier, 0 if identifier in empty_ids else 200,
            '' if identifier in empty_ids else PEAKS,
        ) for identifier in ids
    ))
    return path


def run(command, cwd):
    result = subprocess.run([str(v) for v in command], cwd=cwd,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


@pytest.fixture
def native_binary(scripts_dir):
    path = scripts_dir.parent / 'binaries/main_execmodule'
    if (platform.system() != 'Linux' or platform.machine() not in ('x86_64', 'amd64')
            or not path.exists() or not os.access(path, os.X_OK)):
        pytest.skip('Bundled GNPS native tools require an executable Linux x86_64 binary')
    return path


def test_indexed_mgf_preserves_sparse_zero_ids_and_all_other_metadata(tmp_path):
    source = write_mgf(tmp_path / 'source.mgf', [55, 0, 32])
    output = tmp_path / 'indexed.mgf'
    mapping = importlib.import_module('spectrum_id_mapping').index_mgf_ids(source, output)
    assert mapping == {1: 55, 2: 0, 3: 32}
    assert [int(s['SCANS']) for s in mgf_headers(output)] == [1, 2, 3]
    assert source.read_text() == output.read_text().replace('SCANS=1\n', 'SCANS=55\n').replace(
        'SCANS=2\n', 'SCANS=0\n').replace('SCANS=3\n', 'SCANS=32\n')


def test_restoring_network_ids_uses_mgf_order_and_preserves_other_columns(tmp_path):
    source = write_mgf(tmp_path / 'source.mgf', [55, 0, 32])
    pairs = tmp_path / 'pairs.tsv'
    pairs.write_text('CLUSTERID1\tCLUSTERID2\tCosine\tDeltaMZ\n1\t3\t0.95\t1.2\n2\t1\t0.85\t3.4\n')
    importlib.import_module('spectrum_id_mapping').restore_networking_ids(source, pairs)
    with pairs.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert [(int(r['CLUSTERID1']), int(r['CLUSTERID2'])) for r in rows] == [(55, 32), (0, 55)]
    assert [(float(r['Cosine']), float(r['DeltaMZ'])) for r in rows] == [(0.95, 1.2), (0.85, 3.4)]


@pytest.mark.parametrize('source_ids,empty_ids,expected', [
    ([0, 1, 2], [1], {0, 2}),
    ([1, 2, 3], [2], {1, 3}),
])
def test_native_network_pairs_restore_original_ids_with_padded_scans(
        tmp_path, scripts_dir, native_binary, source_ids, empty_ids, expected):
    spectra = write_mgf(tmp_path / 'query.mgf', source_ids, empty_ids)
    params = tmp_path / 'params'
    params.mkdir()
    run([sys.executable, scripts_dir / 'prep_molecular_networking_parameters.py',
         spectra, params, '--parallelism', '1', '--min_matched_peaks', '6',
         '--min_cosine', '0.7'], tmp_path)
    pairs = tmp_path / 'pairs.tsv'
    run([native_binary, 'ExecMolecularParallelPairs', params / '0.params',
         '-ccms_INPUT_SPECTRA_MS2', spectra, '-ccms_output_aligns', pairs], tmp_path)
    importlib.import_module('spectrum_id_mapping').restore_networking_ids(spectra, pairs)
    with pairs.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert rows, 'The local identical six-peak spectra must produce a network edge'
    actual = {int(r[c]) for r in rows for c in ('CLUSTERID1', 'CLUSTERID2')}
    assert actual == expected


@pytest.mark.parametrize('source_ids', [[0, 1, 55], [1, 2, 55]])
def test_native_library_search_restores_zero_without_colliding_with_one(
        tmp_path, scripts_dir, native_binary, source_ids):
    spectra = write_mgf(tmp_path / 'query.mgf', source_ids)
    library = write_mgf(tmp_path / 'library.mgf', [1])
    library.write_text(library.read_text().replace('TITLE=source_1',
                       'TITLE=local_reference\nNAME=local_reference\nSPECTRUMID=LOCAL_REFERENCE_1'))
    results = tmp_path / 'results'
    results.mkdir()
    run([sys.executable, scripts_dir / 'library_search_wrapper.py', spectra, library,
         results, scripts_dir.parent / 'binaries/convert', native_binary,
         '--pm_tolerance', '0.01', '--fragment_tolerance', '0.01',
         '--library_min_cosine', '0.7', '--library_min_matched_peaks', '6'], tmp_path)
    output = results / (spectra.name + '_' + library.name + '.tsv')
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == len(source_ids), 'Each query must match the identical local reference'
    assert sorted(int(row['#Scan#']) for row in rows) == sorted(source_ids)
    assert {row['SpectrumFile'] for row in rows} == {spectra.name}
    assert {row['FileScanUniqueID'] for row in rows} == {
        spectra.name + str(identifier) for identifier in source_ids
    }
