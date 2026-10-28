# Scoped MS-DIAL conversion tests

Run from the repository root:

```bash
python -m pytest -q test/tools/msdial
```

These tests run offline. The small, checked-in `fixtures/` files exercise the
legacy `MSDIAL` and separate `MSDIAL5` branches through the real formatter entry
point as well as direct conversion functions. The runtime dependencies in
`bin/conda_env.yml` and pytest are required; initialize the `GNPS_sharedcode`
submodule for the dispatch/filtering tests.

The synthetic files have realistic 22-annotation legacy/compact-v5 and
35-annotation newer-v5 layouts. They deliberately use sparse IDs, ID zero,
out-of-order spectra, distinct precursor/RT pairs, two samples, a feature without
MS/MS, and fragment peaks straddling the existing 1% intensity filter. They assert
original shared IDs, correct sample boundaries and values, acquisition metadata,
end-to-end alignment, and preservation of ID zero through spectra filtering.

The native-tool tests exercise the bundled networking and library-search
executables with tiny local spectra. They verify restoration of positional IDs,
zero/one collision handling, positive-only regressions, and stable per-query IDs.
Native tests skip on platforms where the bundled Linux x86_64 executable cannot
run; the local mapping and conversion tests still run.

The v5 preservation tests fail normally on the original broken implementation;
they are not marked `xfail` and do not treat wrong output as an expected result.
The legacy path, compact sample layout, and existing fragment filter have passing
regression tests. Generic runtime guard tolerance/summary tests remain in
`bin/scripts/test_validate_feature_spectrum_alignment.py`.

Optional real-data tests use the already-local files under `data/msdial/`, or a
path set with `GNPS_MSDIAL_TEST_DATA`. They skip if those files are absent and
**never contact a network service**. To populate them manually before testing:

```bash
python data/download_msdial_test_data.py
```

Those tests compare both 2019 legacy conversions with their historical reference
CSVs, convert both Interlab v5 task input pairs and verify source IDs/sample
headers/runtime alignment, and ensure the archived known-bad published outputs
are rejected. Archived buggy outputs are evidence inputs, not golden conversion
outputs. Large fixture data must remain uncommitted.
