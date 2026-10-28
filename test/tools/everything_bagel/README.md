# Everything Bagel import normalization

The EB formatter reuses MZmine feature-table normalization and streams the MGF,
removing only zero-precursor blocks with no fragment lines or real feature/source
metadata, and a generated `SCAN=...` title (or no title). It preserves all real
blocks verbatim, including zero/sparse IDs, rescued MS2, source references and
acquisition metadata. Duplicate, nested, unterminated or unidentified blocks fail
instead of disappearing. Invalid real precursors are retained for the strict
runtime validator to reject.

The tests use local fixtures. They exercise the actual dispatcher, verify MZmine
remains a byte-for-byte passthrough, and compare existing filtered output before
and after EB normalization. The legacy downstream padding in `filter_spectra.py`
is unchanged in this PR; removing it can be considered separately now that native
tool result IDs are explicitly mapped back.
