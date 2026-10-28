# Cross-tool ID and precursor validation

Run offline from the repository root:

```bash
python -m pytest -q test/tools/alignment test/tools/everything_bagel
```

The checked-in fixtures are small excerpts from the six supplied Asari,
Everything Bagel, MassCube, MZmine, OpenMS and XCMS tasks. `fixtures/provenance.json`
records task IDs and selected feature IDs. Fragment arrays are truncated to six
peaks (and `NUM PEAKS` adjusted); these fixtures test ID/precursor contracts, not
spectral matching or metabolite identification. The EB fixture includes a rescued
representative with a 414 ppm difference, still within the chosen **1 Da** limit,
and one known empty padding block. No RT is inspected.

Tests verify matching inputs, missing IDs, mismatched precursors, and the real
workflow shell gate. The gate covers MS-DIAL5 and legacy MS-DIAL (whose formatters
have separate tests), plus the six importer paths above. Other supported input
formats are not enabled without their own verified contracts.

Optional full-data checks use `data/feature_finder_validation/` or the directory in
`GNPS_FEATURE_FINDER_TEST_DATA`. Each tool directory must contain the downloaded
task files at `nf_output/clustering/featuretable_reformated.csv` and
`nf_output/clustering/specs_ms.mgf`. Tests skip if these are absent, and never
download data. The original MS-DIAL5 task is an intentionally rejected negative
case; its buggy outputs are not a golden conversion result.

This check validates IDs and precursor consistency only. It does not certify
sample abundance completeness, filename-to-metadata matching, or assignments
between isobaric features. In particular, XCMS's known abundance-column issue
is a separate concern; passing the ID check does not bless that data loss.
