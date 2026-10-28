The MS-DIAL fixture manifest is tracked; the large downloaded files are local.

```bash
python data/download_msdial_test_data.py
python data/download_msdial_test_data.py --verify-only
python -m pytest bin/scripts/test_validate_feature_spectrum_alignment.py
```

The default destination is `data/msdial/` (approximately 168 MiB). Use `--case`
to select one case or `--destination` to use another directory. Files are checked
against their recorded sizes and SHA-256 checksums, and matching files are reused.
GNPS2 downloads are public; no authentication is needed. Historical source URLs
are recorded in `msdial_manifest.json`. Changed remote content fails verification
rather than silently replacing the baseline.

The two 2019 fixtures use the legacy `MSDIAL` conversion path and include correct
reference feature tables. The two `MSDIAL5` cases are the original Interlab task
and its extensionless-metadata rerun. Their `observed_outputs/` files contain the
reported bugs and are evidence of prior behavior, not correct expected outputs.
Source task parameters reference server paths and libraries not bundled here.

To keep the downloads out of commits without changing shared ignore rules:

```bash
printf '\n/data/msdial/\n' >> .git/info/exclude
```

The workflow enforces the runtime alignment check for **MSDIAL5, MSDIAL, MZMINE,
EVERYTHING_BAGEL, ASARI, MASSCUBE, OPENMS, and XCMS3**, after conversion and before
spectra filtering and networking. These importer contracts have local fixtures.
Other importer paths remain unchanged until validated with their own data.
The EB formatter removes recognizable empty padding without renumbering IDs;
real or malformed zero-precursor spectra are not silently discarded.
The generic checker remains available for explicit diagnostics. It compares IDs and precursor m/z
with a **1 Da** absolute tolerance. It does not inspect retention times.
The checker allows feature rows without MS/MS. It does not infer ID mappings, repair IDs, or change input
data. The failure message is written directly to `nf_output/workflow_summary.txt`
before the process exits nonzero, so it survives a failed Nextflow process.

Cross-tool tests and their source-task provenance are under `test/tools/alignment/`;
EB-specific normalization tests are under `test/tools/everything_bagel/`. They use
small checked-in local fixtures and never contact external services. Optional
full-task tests use existing local files selected by `GNPS_FEATURE_FINDER_TEST_DATA`.
