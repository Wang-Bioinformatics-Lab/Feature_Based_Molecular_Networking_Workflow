#!/usr/bin/env python3
"""
mztab_library_merge.py

Reads a GNPS library-match TSV and an mzTab file, then writes a new mzTab
file that is identical to the input except that the SME section is rebuilt
entirely from library-match data.

Linkage chain
─────────────
  library '#Scan#'  →  SMF table 'SMF_ID'
  SMF 'SME_ID_REFS' →  the SME_ID(s) to populate

Column population rules per the mzTab-M 2.0 spec
─────────────────────────────────────────────────
Column                         | Nullable | Source
─────────────────────────────────────────────────────────────────────────────
SME_ID                         | FALSE    | preserved from original SME table
evidence_input_id              | FALSE    | built from SMF: "{rt}_{mz}m/z"
database_identifier            | TRUE     | library: SpectrumID
chemical_formula               | TRUE     | library: molecular_formula
smiles                         | TRUE     | library: Smiles
inchi                          | TRUE     | library: INCHI
chemical_name                  | TRUE     | library: Compound_Name
uri                            | TRUE     | library: library_usi
derivatized_form               | TRUE     | null  (nullable, no source)
adduct_ion                     | TRUE     | library: Adduct
exp_mass_to_charge             | FALSE    | carried from SMF: exp_mass_to_charge
charge                         | FALSE    | carried from SMF: charge
theoretical_mass_to_charge     | FALSE    | library: Precursor_MZ
spectra_ref                    | FALSE    | carried verbatim from original SME row
identification_method          | FALSE    | fixed: "[,,GNPS2 Feature Based Molecular Networking,]"
ms_level                       | FALSE    | fixed: "[MS,MS:1000511,ms level,2]"
id_confidence_measure[1]       | TRUE     | library: MQScore
id_confidence_measure[2]       | TRUE     | library: TIC_Query
id_confidence_measure[3]       | TRUE     | library: MZErrorPPM
rank                           | FALSE    | fixed: 1
opt_global_*                   | TRUE     | null  (optional columns, dropped/nulled)

For SME rows with NO library match every nullable column → "null";
non-nullable columns (exp_mass_to_charge, charge, theoretical_mass_to_charge,
spectra_ref, identification_method, ms_level, rank) are still populated from
SMF / fixed values as above so the file stays spec-compliant.

Usage
─────
    python mztab_library_merge.py \\
        --library  <merged_results_with_gnps.tsv> \\
        --mztab    <quantification_table.mztab> \\
        --output   <updated_table.mztab> \\
        [--scan-col "#Scan#"] \\
        [--verbose]
"""

import argparse
import csv
import re
from pathlib import Path
from typing import Optional


# ─────────────────────────────────────────────────────────────────────────────
# Fixed / constant values
# ─────────────────────────────────────────────────────────────────────────────
FIXED_IDENTIFICATION_METHOD = "[,,GNPS2 Feature Based Molecular Networking,]"
FIXED_MS_LEVEL              = "[MS,MS:1000511,ms level,2]"
FIXED_RANK                  = "1"


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def parse_pipe_refs(s: str) -> list:
    """Split a '|'-delimited ref string; return stripped, non-empty tokens."""
    return [t.strip() for t in s.split("|") if t.strip()]


def lib_val(lib_row: dict, col: str) -> str:
    """Return library value for col, or '' if absent/empty."""
    v = lib_row.get(col, "")
    return v.strip() if v else ""


# ─────────────────────────────────────────────────────────────────────────────
# mzTab raw parser
# ─────────────────────────────────────────────────────────────────────────────

class MztabRaw:
    """
    Stores raw lines of an mzTab file and exposes parsed SMF and SME data
    so we can reconstruct the file with a replacement SME section.
    """

    def __init__(self, path: str) -> None:
        with open(path, encoding="utf-8") as fh:
            self.raw_lines: list = fh.readlines()

        self.total_cols: int = max(
            line.count("\t") + 1 for line in self.raw_lines
        )
        self._parse()

    # ------------------------------------------------------------------
    def _parse(self) -> None:
        header_to_prefix = {"SMH": "SML", "SFH": "SMF", "SEH": "SME"}
        current_prefix   = None
        current_headers  = []

        self.sme_headers:    list = []
        self.seh_line:       int  = -1
        self.sme_line_start: int  = -1
        self.sme_line_end:   int  = -1

        # SMF rows keyed by SMF_ID (int)
        self.smf_by_id:  dict = {}
        # SME rows: list of dicts (original, used only for spectra_ref carry-over)
        self.sme_rows:   list = []

        sfh_headers: list = []

        for i, line in enumerate(self.raw_lines):
            cols   = line.rstrip("\n").split("\t")
            prefix = cols[0].strip()

            if not prefix:
                if current_prefix == "SME" and self.sme_line_end == -1:
                    self.sme_line_end = i
                current_prefix  = None
                current_headers = []
                continue

            if prefix in header_to_prefix:
                current_prefix  = header_to_prefix[prefix]
                current_headers = [c.strip() for c in cols[1:]]
                # strip trailing empty strings from headers
                while current_headers and current_headers[-1] == "":
                    current_headers.pop()

                if current_prefix == "SME":
                    self.seh_line    = i
                    self.sme_headers = list(current_headers)
                elif current_prefix == "SMF":
                    sfh_headers = list(current_headers)

            elif prefix == current_prefix:
                values = [v.strip() for v in cols[1:]]
                while len(values) < len(current_headers):
                    values.append("")
                row = dict(zip(current_headers, values[: len(current_headers)]))

                if current_prefix == "SME":
                    if self.sme_line_start == -1:
                        self.sme_line_start = i
                    self.sme_line_end = i + 1
                    self.sme_rows.append(row)

                elif current_prefix == "SMF":
                    try:
                        smf_id = int(row.get("SMF_ID", ""))
                        self.smf_by_id[smf_id] = row
                    except ValueError:
                        pass

        if self.sme_line_end == -1 and self.sme_line_start != -1:
            self.sme_line_end = len(self.raw_lines)

    # ------------------------------------------------------------------
    def smf_to_sme_map(self) -> dict:
        """Return {smf_id (int): [sme_id_str, ...]}."""
        result = {}
        for smf_id, row in self.smf_by_id.items():
            refs = parse_pipe_refs(row.get("SME_ID_REFS", ""))
            if refs:
                result[smf_id] = refs
        return result

    # ------------------------------------------------------------------
    def original_sme_index(self) -> dict:
        """Return {sme_id_str: original_sme_row_dict} for carry-over fields."""
        return {row["SME_ID"]: row for row in self.sme_rows if row.get("SME_ID")}


# ─────────────────────────────────────────────────────────────────────────────
# Build evidence_input_id from SMF row
# ─────────────────────────────────────────────────────────────────────────────

def make_evidence_input_id(smf_row: dict) -> str:
    """
    Format: "{retention_time_in_seconds}_{exp_mass_to_charge}m/z"
    e.g.  "717.05_303.0674m/z"
    """
    rt  = smf_row.get("retention_time_in_seconds", "null")
    mz  = smf_row.get("exp_mass_to_charge", "null")
    return f"{rt}_{mz}m/z"


# ─────────────────────────────────────────────────────────────────────────────
# SME row renderer
# ─────────────────────────────────────────────────────────────────────────────

def render_sme_row(
    sme_id:      str,
    smf_row:     Optional[dict],   # SMF row linked to this SME (may be None)
    orig_sme:    Optional[dict],   # original SME row (for spectra_ref)
    lib_row:     Optional[dict],   # library match row (may be None)
    sme_headers: list,
    total_cols:  int,
) -> str:
    """
    Build one tab-delimited SME line.

    Non-nullable columns are always populated (from SMF / fixed values / lib).
    Nullable columns use library data when available, otherwise "null".
    """
    # ── per-column logic ──────────────────────────────────────────────────────
    def cell(col: str) -> str:
        matched = lib_row is not None

        # ── NON-NULLABLE ──────────────────────────────────────────────────────
        if col == "SME_ID":
            return sme_id

        if col == "evidence_input_id":
            # Build from SMF if available, else "null" would violate spec –
            # but SMF should always be present for every SME in this file.
            if smf_row:
                return make_evidence_input_id(smf_row)
            # Fall back to original value rather than "null" (non-nullable)
            if orig_sme:
                return orig_sme.get("evidence_input_id", "null")
            return "null"

        if col == "exp_mass_to_charge":
            if smf_row:
                return smf_row.get("exp_mass_to_charge", "null")
            if orig_sme:
                return orig_sme.get("exp_mass_to_charge", "null")
            return "null"

        if col == "charge":
            if smf_row:
                return smf_row.get("charge", "null")
            if orig_sme:
                return orig_sme.get("charge", "null")
            return "null"

        if col == "theoretical_mass_to_charge":
            if matched:
                v = lib_val(lib_row, "Precursor_MZ")
                return v if v else "null"
            # non-nullable: use exp_mass_to_charge as proxy when no match
            if smf_row:
                return smf_row.get("exp_mass_to_charge", "null")
            return "null"

        if col == "spectra_ref":
            # carry verbatim from original SME row
            if orig_sme:
                v = orig_sme.get("spectra_ref", "").strip()
                return v if v else "null"
            return "null"

        if col == "identification_method":
            return FIXED_IDENTIFICATION_METHOD

        if col == "ms_level":
            return FIXED_MS_LEVEL

        if col == "rank":
            return FIXED_RANK

        # ── NULLABLE ─────────────────────────────────────────────────────────
        if col == "database_identifier":
            if matched:
                v = lib_val(lib_row, "SpectrumID")
                return v if v else "null"
            return "null"

        if col == "chemical_formula":
            if matched:
                v = lib_val(lib_row, "molecular_formula")
                return v if v else "null"
            return "null"

        if col == "smiles":
            if matched:
                v = lib_val(lib_row, "Smiles")
                return v if v else "null"
            return "null"

        if col == "inchi":
            if matched:
                v = lib_val(lib_row, "INCHI")
                return v if v else "null"
            return "null"

        if col == "chemical_name":
            if matched:
                v = lib_val(lib_row, "Compound_Name")
                return v if v else "null"
            return "null"

        if col == "uri":
            if matched:
                v = lib_val(lib_row, "library_usi")
                return v if v else "null"
            return "null"

        if col == "derivatized_form":
            return "null"

        if col == "adduct_ion":
            if matched:
                v = lib_val(lib_row, "Adduct")
                return v if v else "null"
            return "null"

        if col == "id_confidence_measure[1]":
            if matched:
                v = lib_val(lib_row, "MQScore")
                return v if v else "null"
            return "null"

        if col == "id_confidence_measure[2]":
            if matched:
                v = lib_val(lib_row, "TIC_Query")
                return v if v else "null"
            return "null"

        if col == "id_confidence_measure[3]":
            if matched:
                v = lib_val(lib_row, "MZErrorPPM")
                return v if v else "null"
            return "null"

        # opt_* columns and anything else not explicitly handled → null
        return "null"

    # ── assemble values ───────────────────────────────────────────────────────
    values = [cell(col) for col in sme_headers]

    # Pad with empty strings to the file-wide column count
    while len(values) < total_cols - 1:
        values.append("")

    return "SME\t" + "\t".join(values) + "\n"


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main(args) -> None:
    verbose = args.verbose

    # ── load library matches ──────────────────────────────────────────────────
    if verbose:
        print(f"[*] Loading library results: {args.library}")
    with open(args.library, newline="", encoding="utf-8") as fh:
        lib_rows = list(csv.DictReader(fh, delimiter="\t"))
    if verbose:
        print(f"    {len(lib_rows)} library match row(s)")

    # scan# (= SMF_ID) → library row  (first match wins for duplicates)
    scan_to_lib: dict = {}
    for row in lib_rows:
        try:
            smf_id = int(row[args.scan_col])
            scan_to_lib.setdefault(smf_id, row)
        except (KeyError, ValueError):
            if verbose:
                print(f"    [!] Cannot parse scan from: {row}")

    # ── parse mzTab ───────────────────────────────────────────────────────────
    if verbose:
        print(f"[*] Parsing mzTab: {args.mztab}")
    mztab = MztabRaw(args.mztab)

    smf_to_sme  = mztab.smf_to_sme_map()    # {smf_id: [sme_id, ...]}
    orig_sme_ix = mztab.original_sme_index() # {sme_id: orig_sme_row}
    all_sme_ids = [row["SME_ID"] for row in mztab.sme_rows if row.get("SME_ID")]

    if verbose:
        print(f"    SMF rows: {len(mztab.smf_by_id)},  SME rows: {len(all_sme_ids)}")
        print(f"    SME columns ({len(mztab.sme_headers)}): {mztab.sme_headers}")
        print(f"    File-wide column count: {mztab.total_cols}")

    # ── build SME_ID → (smf_row, lib_row) ────────────────────────────────────
    # Reverse index: sme_id → smf_row (to get exp_mass_to_charge etc.)
    sme_to_smf: dict = {}
    for smf_id, smf_row in mztab.smf_by_id.items():
        for sme_id in smf_to_sme.get(smf_id, []):
            sme_to_smf.setdefault(sme_id, smf_row)

    # sme_id → library row (via SMF_ID match)
    sme_to_lib: dict = {}
    for smf_id, lib_row in scan_to_lib.items():
        for sme_id in smf_to_sme.get(smf_id, []):
            sme_to_lib.setdefault(sme_id, lib_row)

    if verbose:
        print(f"    SME IDs with SMF reference: {sorted(sme_to_smf.keys())}")
        print(f"    SME IDs matched to library: {sorted(sme_to_lib.keys())}")

    # ── build replacement SME lines ───────────────────────────────────────────
    new_sme_lines = []
    matched = 0
    for sme_id in all_sme_ids:
        smf_row  = sme_to_smf.get(sme_id)
        lib_row  = sme_to_lib.get(sme_id)
        orig_sme = orig_sme_ix.get(sme_id)
        if lib_row:
            matched += 1
        new_sme_lines.append(
            render_sme_row(
                sme_id      = sme_id,
                smf_row     = smf_row,
                orig_sme    = orig_sme,
                lib_row     = lib_row,
                sme_headers = mztab.sme_headers,
                total_cols  = mztab.total_cols,
            )
        )

    # ── assemble output ───────────────────────────────────────────────────────
    out_lines = (
        mztab.raw_lines[: mztab.seh_line]        # everything before SEH
        + [mztab.raw_lines[mztab.seh_line]]       # SEH header verbatim
        + new_sme_lines                           # new SME data
        + mztab.raw_lines[mztab.sme_line_end:]   # anything after last SME
    )

    out_path = Path(args.output)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.writelines(out_lines)

    if verbose:
        unmatched = len(all_sme_ids) - matched
        print(
            f"[✓] {len(all_sme_ids)} SME rows written "
            f"({matched} library-matched, {unmatched} null-filled) → {out_path}"
        )
    else:
        print(f"Done. Output written to {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser():
    p = argparse.ArgumentParser(
        description=(
            "Merge GNPS library-match results into an mzTab-M file. "
            "The output mzTab is identical to the input except the SME section "
            "is rebuilt entirely from library-match data per the mzTab-M 2.0 spec."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--library", "-l", required=True, metavar="FILE",
                   help="Path to the GNPS library-match TSV file")
    p.add_argument("--mztab",   "-m", required=True, metavar="FILE",
                   help="Path to the input mzTab file")
    p.add_argument("--output",  "-o", default="updated_table.mztab", metavar="FILE",
                   help="Path for the output mzTab file (default: updated_table.mztab)")
    p.add_argument("--scan-col", default="#Scan#", metavar="COL",
                   help="Library TSV column whose value equals SMF_ID (default: '#Scan#')")
    p.add_argument("--verbose", "-v", action="store_true",
                   help="Print progress information")
    return p


if __name__ == "__main__":
    parser = build_parser()
    args   = parser.parse_args()
    main(args)
