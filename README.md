# Nextflow Template

To run the workflow to test simply do

```
make run
```

To learn NextFlow checkout this documentation:

https://www.nextflow.io/docs/latest/index.html

## Library Annotation (offline mode)

Library matches are annotated by `bin/scripts/getGNPS_library_annotations.py`. The
`forceoffline` parameter (default `Yes`) controls how the compound metadata
(Compound_Name, Smiles, Adduct, …) is resolved for each hit:

- `forceoffline = No` — annotate GNPS (`CCMSLIB...`) hits by calling the live GNPS
  API per hit.
- `forceoffline = Yes` — skip the API and resolve every hit locally from a
  `library_summary.tsv` that is built from the library `.mgf` files by
  `bin/scripts/library_summary.py`. This is much faster (no per-hit network calls).

**Important — non-unique SCANS in the current GNPS2 libraries.** Many of the
libraries served from https://library.gnps2.org reuse the same `SCANS=` value
across spectra (SCANS is a per-source-file counter, not a global unique id — e.g.
`GNPS-MSMLS.mgf` has 863 spectra but only 56 distinct SCANS). `library_summary.py`
must therefore read with `mgf.read(..., use_index=False)` and key on the real
unique id (`SPECTRUMID`). Reading with `index_by_scans=True` collapses all spectra
that share a SCANS value down to one, dropping ~94% of the library from the summary,
so offline lookups miss and every field comes back `N/A`.

## Shared code / modules

The scripts under `bin/scripts/` are mostly symlinks into the `GNPS_sharedcode`
submodule. `GNPS_sharedcode` (github.com/mwang87/GNPS_sharedcode) is being
**deprecated** in favor of the consolidated `NextflowModules` submodule
(github.com/Wang-Bioinformatics-Lab/NextflowModules), which vendors a copy of the
shared code plus the workflow-specific, fixed scripts. `library_summary.py` was
removed from `GNPS_sharedcode` as part of that deprecation, so it is checked in
here as a **real file** (`bin/scripts/library_summary.py`) rather than a symlink.
It is kept byte-for-byte identical to the canonical version at
`NextflowModules/bin/library_search/library_summary.py`.

TODO: migrate this workflow onto the `NextflowModules` submodule (as
`LibrarySearch_Workflow` already has) and repoint the remaining `bin/scripts/`
symlinks to it, so shared code has a single source of truth.

## Installation

You will need to have conda, mamba, and nextflow installed. 

## Deployment to GNPS2

In order to deploy, we have a set of deployment tools that will enable deployment to the various gnps systems. To run the deployment, use the following commands from the deploy_gnps2 folder. 

```
make deploy-prod
```

You might need to checkout the module, do this by running

```
git submodule init
git submodule update
```
