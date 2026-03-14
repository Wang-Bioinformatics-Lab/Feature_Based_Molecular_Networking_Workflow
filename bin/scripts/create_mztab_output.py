#!/usr/bin/env python3
"""
create_mztab_output.py

Produces a valid mzTab-M 2.0 (Summary mode) file from FBMN workflow outputs.
Runs for ALL featurefindingtool inputs (MZMINE, OpenMS, XCMS3, etc.).

Inputs
------
--cluster_summary   enriched cluster summary TSV  (from enrichClusterSummary)
--library_results   merged library results TSV    (from librarygetGNPSAnnotations)
--metadata          merged metadata TSV            (from createMetadataFile; optional / NO_FILE)
--feature_table     reformatted feature table CSV  (from quantification_table_reformatted)
--ms_run_file       path or name of the spectra file used during networking
--featurefindingtool  tool name (informational; goes into MTD description)
--output            output .mztab file path

Output structure (mzTab-M 2.0 Summary)
---------------------------------------
MTD  – metadata block
SML  – one row per feature  (abundances per assay)
SMF  – one row per feature  (same as SML in Summary mode)
SME  – one row per library hit (features without a hit still get a null-filled row)
"""

import argparse
import os
import re
import sys
import math
from datetime import date
from io import StringIO
from urllib.parse import quote, unquote

import pandas as pd


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _read_workflow_version(yaml_path: str) -> str:
    """Extract workflowversion from workflowinput.yaml without requiring PyYAML."""
    if not yaml_path or not os.path.exists(yaml_path):
        return "Unknown"
    with open(yaml_path, encoding="utf-8") as fh:
        for line in fh:
            m = re.match(r'^workflowversion\s*:\s*["\']?([^"\'#\n]+)["\']?', line)
            if m:
                return m.group(1).strip()
    return "Unknown"


def _null(v) -> str:
    """Return 'null' when a value is missing / NaN / empty."""
    if v is None:
        return "null"
    if isinstance(v, float) and math.isnan(v):
        return "null"
    s = str(v).strip()
    if s in ("", "nan", "NaN", "NA", "N/A", "None"):
        return "null"
    return s


def _tab(*cells) -> str:
    return "\t".join(str(c) for c in cells) + "\n"


def _to_uri(name: str) -> str:
    """Convert a sample/file name into a valid file:/// URI."""
    # Percent-encode characters that are not valid in a URI path
    encoded = quote(name, safe="/-._~")
    return "file:///{}".format(encoded)


def _format_adduct(raw: str) -> str:
    """Normalize adduct string to mzTab-M 2.0 format: [M+X]c.

    The spec requires adducts to match ^\\[\\d*M([+-][\\w\\d]+)*\\]\\d*[+-]$
    e.g. [M+H]+, [M-H]-, [M+Na]+, [2M+H]+
    """
    if not raw or raw == "null":
        return "null"
    s = raw.strip()
    # Already in correct format
    if re.match(r"^\[\d*M([+-][\w\d]+)*\]\d*[+-]$", s):
        return s
    # Common shorthand like "M+H" -> "[M+H]+"
    m = re.match(r"^(\d*M([+-][\w\d]+)*)$", s)
    if m:
        adduct_core = m.group(1)
        # Infer charge sign from the adduct
        if any(x in s for x in ["+H", "+Na", "+K", "+NH4", "+Li"]):
            return "[{}]+".format(adduct_core)
        elif any(x in s for x in ["-H", "-Cl", "-FA"]):
            return "[{}]-".format(adduct_core)
        else:
            return "[{}]+".format(adduct_core)
    # If it has brackets but missing charge sign
    m2 = re.match(r"^\[(\d*M([+-][\w\d]+)*)\]$", s)
    if m2:
        adduct_core = m2.group(1)
        if "-H" in s or "-Cl" in s:
            return "[{}]-".format(adduct_core)
        return "[{}]+".format(adduct_core)
    return "null"


def _infer_ms_run_format(filename: str) -> tuple:
    """
    Infer the correct CV terms for ms_run format and id_format from the
    file extension of the spectra file.

    Returns (format_cv, id_format_cv) as mzTab-M CV param strings.

    CV terms used:
      MS:1000584  mzML format          (mzML files)
      MS:1001062  Mascot MGF format    (MGF files — the standard FBMN input)
      MS:1001530  mzML unique identifier       (id_format for mzML)
      MS:1000774  multiple peak list nativeID format  (id_format for MGF/PKL)
    """
    ext = os.path.splitext(filename)[1].lower()
    if ext == ".mzml":
        fmt = "[MS, MS:1000584, mzML format, ]"
        id_fmt = "[MS, MS:1001530, mzML unique identifier, ]"
    elif ext in (".mgf", ".pkl"):
        # CHANGED from MS:1000564 (PSI mzData format — wrong) to:
        #   MS:1001062  Mascot MGF format
        #   MS:1000774  multiple peak list nativeID format
        fmt = "[MS, MS:1001062, Mascot MGF format, ]"
        id_fmt = "[MS, MS:1000774, multiple peak list nativeID format, ]"
    else:
        # Fallback: mzML is the most common vendor-neutral format
        fmt = "[MS, MS:1000584, mzML format, ]"
        id_fmt = "[MS, MS:1001530, mzML unique identifier, ]"
    return fmt, id_fmt


# ---------------------------------------------------------------------------
# data loading
# ---------------------------------------------------------------------------

def load_cluster_summary(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return df


def load_library_results(path: str) -> pd.DataFrame:
    if not path or path == "NO_FILE" or not os.path.exists(path):
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return pd.DataFrame()
        return df
    except Exception:
        return pd.DataFrame()


def load_metadata(path: str) -> pd.DataFrame:
    if not path or path == "NO_FILE" or not os.path.exists(path):
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty or "filename" not in df.columns:
            return pd.DataFrame()
        return df
    except Exception:
        return pd.DataFrame()


def load_feature_table(path: str) -> pd.DataFrame:
    if not path or path == "NO_FILE" or not os.path.exists(path):
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep=",")
        if df.empty:
            return pd.DataFrame()
        return df
    except Exception:
        return pd.DataFrame()


def _split_pipe_refs(value) -> list:
    if value is None:
        return []
    return [token.strip() for token in str(value).split("|") if token.strip()]


def _next_index(records: list, prefix: str) -> int:
    indices = []
    pattern = re.compile(r"^{}\[(\d+)\]".format(re.escape(prefix)))
    for key, _ in records:
        m = pattern.match(key)
        if m:
            indices.append(int(m.group(1)))
    return (max(indices) + 1) if indices else 1


def _basename_from_location(location: str) -> str:
    raw = str(location).strip()
    if raw.startswith("file://"):
        raw = raw[7:]
    raw = raw.replace("\\", "/")
    decoded = unquote(raw)
    return os.path.basename(decoded) or os.path.basename(str(location).strip())


def _parse_table_section(lines: list) -> pd.DataFrame:
    if not lines:
        return pd.DataFrame()
    try:
        df = pd.read_csv(StringIO("".join(lines)), sep="\t")
        if df.empty:
            return pd.DataFrame()
        df.columns = df.columns.str.strip()
        df = df.loc[:, ~df.columns.duplicated(keep="first")]
        return df
    except Exception:
        return pd.DataFrame()


def load_source_mztab_context(path: str) -> dict:
    if not path or path == "NO_FILE" or not os.path.exists(path):
        return {}

    try:
        with open(path, encoding="utf-8") as handle:
            raw_lines = handle.readlines()
    except Exception:
        return {}

    mtd_lines = []
    smf_lines = []
    sme_lines = []
    for line in raw_lines:
        prefix = line.split("\t")[0].strip()
        if prefix == "MTD":
            mtd_lines.append(line)
        elif prefix in ("SFH", "SMF"):
            smf_lines.append(line)
        elif prefix in ("SEH", "SME"):
            sme_lines.append(line)

    records = []
    if mtd_lines:
        try:
            mtd_df = pd.read_csv(StringIO("".join(mtd_lines)), sep="\t", header=None)[[1, 2]]
            for _, row in mtd_df.iterrows():
                key = str(row.iloc[0]).strip()
                value = str(row.iloc[1]).strip()
                records.append((key, value))
        except Exception:
            records = []

    ms_runs = {}
    sample_defs = {}
    assay_defs = {}
    study_var_defs = {}

    for key, value in records:
        ms_match = re.match(r"^ms_run\[(\d+)\]-(.+)$", key)
        sample_match = re.match(r"^sample\[(\d+)\](?:-(.+))?$", key)
        assay_match = re.match(r"^assay\[(\d+)\](?:-(.+))?$", key)
        sv_match = re.match(r"^study_variable\[(\d+)\](?:-(.+))?$", key)

        if ms_match:
            idx = int(ms_match.group(1))
            field = ms_match.group(2)
            ms_runs.setdefault(idx, {})[field] = value
        elif sample_match:
            idx = int(sample_match.group(1))
            field = sample_match.group(2)
            sample_defs.setdefault(idx, {"index": idx})
            if field:
                sample_defs[idx][field] = value
            else:
                sample_defs[idx]["name"] = value
        elif assay_match:
            idx = int(assay_match.group(1))
            field = assay_match.group(2)
            assay_defs.setdefault(idx, {"index": idx})
            if field:
                assay_defs[idx][field] = value
            else:
                assay_defs[idx]["name"] = value
        elif sv_match:
            idx = int(sv_match.group(1))
            field = sv_match.group(2)
            study_var_defs.setdefault(idx, {"index": idx})
            if field:
                study_var_defs[idx][field] = value
            else:
                study_var_defs[idx]["name"] = value

    assays = []
    for idx in sorted(assay_defs):
        assay = assay_defs[idx]
        ms_run_ref = assay.get("ms_run_ref", "")
        ms_run_idx = None
        ms_run_match = re.match(r"^ms_run\[(\d+)\]$", ms_run_ref)
        if ms_run_match:
            ms_run_idx = int(ms_run_match.group(1))
        location = ms_runs.get(ms_run_idx, {}).get("location", "") if ms_run_idx else ""
        basename = _basename_from_location(location) if location else ""
        candidates = []
        for candidate in [basename, os.path.splitext(basename)[0], assay.get("name", "")]:
            candidate = str(candidate).strip()
            if candidate and candidate not in candidates:
                candidates.append(candidate)
        assays.append({
            "index": idx,
            "name": assay.get("name", "assay_{}".format(idx)),
            "sample_ref": assay.get("sample_ref", ""),
            "ms_run_ref": ms_run_ref,
            "abundance_key": basename or assay.get("name", "assay_{}".format(idx)),
            "abundance_candidates": candidates or [assay.get("name", "assay_{}".format(idx))],
        })

    samples = []
    for idx in sorted(sample_defs):
        sample = sample_defs[idx]
        samples.append({
            "index": idx,
            "name": sample.get("name", "sample_{}".format(idx)),
            "fields": {
                key: value for key, value in sample.items()
                if key not in ("index", "name")
            },
        })

    study_variables = []
    for idx in sorted(study_var_defs):
        sv = study_var_defs[idx]
        assay_refs = []
        for ref in _split_pipe_refs(sv.get("assay_refs", "")):
            match = re.match(r"^assay\[(\d+)\]$", ref)
            if match:
                assay_refs.append(int(match.group(1)))
        study_variables.append({
            "index": idx,
            "name": sv.get("name", "study_variable_{}".format(idx)),
            "description": sv.get("description", sv.get("name", "study_variable_{}".format(idx))),
            "assay_refs": assay_refs,
            "fields": {
                key: value for key, value in sv.items()
                if key not in ("index", "name", "description", "assay_refs")
            },
        })

    smf_df = _parse_table_section(smf_lines)
    sme_df = _parse_table_section(sme_lines)
    smf_spectra_refs = {}
    if not smf_df.empty and not sme_df.empty and "SMF_ID" in smf_df.columns and "SME_ID" in sme_df.columns and "spectra_ref" in sme_df.columns:
        sme_df = sme_df.copy()
        sme_df["_SME_ID_STR"] = sme_df["SME_ID"].astype(str).str.strip()
        spectra_lookup = sme_df.set_index("_SME_ID_STR")["spectra_ref"].to_dict()
        for _, row in smf_df.iterrows():
            try:
                smf_id = int(str(row["SMF_ID"]).strip())
            except Exception:
                continue
            refs = []
            for sme_id in _split_pipe_refs(row.get("SME_ID_REFS", "")):
                spectra_ref = spectra_lookup.get(str(sme_id).strip())
                if spectra_ref and spectra_ref not in refs:
                    refs.append(spectra_ref)
            if refs:
                smf_spectra_refs[smf_id] = " | ".join(refs)

    return {
        "ms_runs": ms_runs,
        "samples": samples,
        "assays": assays,
        "study_variables": study_variables,
        "smf_spectra_refs": smf_spectra_refs,
    }


# ---------------------------------------------------------------------------
# sample / assay / study_variable helpers
# ---------------------------------------------------------------------------

def get_sample_columns(feature_df: pd.DataFrame, cluster_df: pd.DataFrame) -> list:
    """Return list of sample column basenames (without ' Peak area' suffix).

    Prefer the reformatted feature table (has per-file columns) over the
    cluster summary (may have group-level ATTRIBUTE_ columns instead).
    """
    # Try feature table first
    for df in [feature_df, cluster_df]:
        if df is not None and not df.empty:
            cols = [
                c.replace(" Peak area", "").strip()
                for c in df.columns
                if c.endswith(" Peak area")
            ]
            if cols:
                return cols

    return []


def build_assays(feature_df: pd.DataFrame, cluster_df: pd.DataFrame,
                 source_ctx: dict = None) -> list:
    source_ctx = source_ctx or {}
    if source_ctx.get("assays"):
        return source_ctx["assays"]

    samples = get_sample_columns(feature_df, cluster_df)
    return [
        {
            "index": i + 1,
            "name": sample,
            "sample_ref": "",
            "ms_run_ref": "ms_run[{}]".format(i + 1),
            "abundance_key": sample,
            "abundance_candidates": [sample],
        }
        for i, sample in enumerate(samples)
    ]


def build_study_variables(assays: list, metadata_df: pd.DataFrame,
                          source_ctx: dict = None) -> list:
    """
    Return a list of dicts:
    [
      {"name": "group1", "assay_refs": [1, 3, 5]},
      ...
    ]
    Groups come from ATTRIBUTE_* columns in metadata.
    If no metadata or no ATTRIBUTE_ columns, a single study_variable with all assays.
    """
    source_ctx = source_ctx or {}
    if source_ctx.get("study_variables"):
        return source_ctx["study_variables"]

    assay_map = {assay["abundance_key"]: assay["index"] for assay in assays}
    assay_count = len(assays)

    if metadata_df.empty:
        return [{"index": 1, "name": "all_samples", "description": "all samples",
                 "assay_refs": list(range(1, assay_count + 1))}]

    attr_cols = [c for c in metadata_df.columns if c.upper().startswith("ATTRIBUTE_")]
    if not attr_cols:
        return [{"index": 1, "name": "all_samples", "description": "all samples",
                 "assay_refs": list(range(1, assay_count + 1))}]

    # Use the first ATTRIBUTE_ column for grouping
    attr = attr_cols[0]

    # Clean metadata filenames to match sample column names
    meta = metadata_df.copy()
    meta["_basename"] = meta["filename"].apply(
        lambda x: os.path.basename(str(x)).strip() if pd.notna(x) else ""
    )

    groups = {}
    for _, row in meta.iterrows():
        fname = row["_basename"]
        group_val = str(row.get(attr, "unknown")).strip()
        if fname in assay_map:
            groups.setdefault(group_val, []).append(assay_map[fname])

    # Samples not in metadata → "unassigned" group
    meta_fnames = set(meta["_basename"].tolist())
    unassigned = [assay_map[s] for s in assay_map if s not in meta_fnames]
    if unassigned:
        groups.setdefault("unassigned", []).extend(unassigned)

    svs = []
    for idx, (g, refs) in enumerate(sorted(groups.items()), 1):
        svs.append({
            "index": idx,
            "name": re.sub(r"[^\w\-. ]", "_", str(g)),
            "description": str(g),
            "assay_refs": sorted(refs),
        })
    return svs


# ---------------------------------------------------------------------------
# MTD section builder
# ---------------------------------------------------------------------------

def build_mtd(assays: list, study_variables: list, ms_run_file: str,
              featurefindingtool: str, workflow_version: str = "Unknown",
              source_ctx: dict = None) -> list:
    source_ctx = source_ctx or {}

    lines = []

    def m(key, val):
        lines.append(_tab("MTD", key, val))

    m("mzTab-version", "2.0.0-M")
    m("mzTab-ID", "FBMN_" + date.today().strftime("%Y%m%d"))
    m("title", "FBMN Results - {}".format(featurefindingtool))
    m("description",
      "Feature-Based Molecular Networking results generated from {} input".format(
          featurefindingtool))

    # Infer correct format/id_format CV terms from file extension.
    # CHANGED: was always MS:1000564 (PSI mzData format) which is incorrect
    # for FBMN inputs (mzML / MGF). Now inferred per file:
    #   mzML → MS:1000584 + MS:1001530
    #   MGF  → MS:1001062 + MS:1000774
    ms_run_fmt, ms_run_id_fmt = _infer_ms_run_format(ms_run_file)

    if source_ctx.get("ms_runs"):
        for i in sorted(source_ctx["ms_runs"]):
            ms_run = source_ctx["ms_runs"][i]
            m("ms_run[{}]-location".format(i), ms_run.get("location", _to_uri("ms_run_{}".format(i))))
            m("ms_run[{}]-format".format(i), ms_run.get("format", ms_run_fmt))
            m("ms_run[{}]-id_format".format(i), ms_run.get("id_format", ms_run_id_fmt))
            m("ms_run[{}]-scan_polarity[1]".format(i),
              ms_run.get("scan_polarity[1]", "[MS, MS:1000130, positive scan, ]"))
    else:
        for assay in assays:
            i = assay["index"]
            loc = _to_uri(assay["abundance_key"])
            m("ms_run[{}]-location".format(i), loc)
            m("ms_run[{}]-format".format(i), ms_run_fmt)
            m("ms_run[{}]-id_format".format(i), ms_run_id_fmt)
            m("ms_run[{}]-scan_polarity[1]".format(i),
              "[MS, MS:1000130, positive scan, ]")

    if source_ctx.get("samples"):
        for sample in source_ctx["samples"]:
            m("sample[{}]".format(sample["index"]), sample["name"])
            for field, value in sample.get("fields", {}).items():
                m("sample[{}]-{}".format(sample["index"], field), value)

    for assay in assays:
        i = assay["index"]
        m("assay[{}]".format(i), assay["name"])
        if assay.get("sample_ref"):
            m("assay[{}]-sample_ref".format(i), assay["sample_ref"])
        m("assay[{}]-ms_run_ref".format(i), assay.get("ms_run_ref", "ms_run[{}]".format(i)))

    for sv in study_variables:
        refs = "|".join("assay[{}]".format(r) for r in sv["assay_refs"])
        m("study_variable[{}]".format(sv["index"]), sv["name"])
        m("study_variable[{}]-description".format(sv["index"]), sv["description"])
        m("study_variable[{}]-assay_refs".format(sv["index"]), refs)
        if source_ctx.get("study_variables"):
            for field, value in sv.get("fields", {}).items():
                if field in ("description", "assay_refs", "variation_function"):
                    continue
                m("study_variable[{}]-{}".format(sv["index"], field), value)
        else:
            # MS:1002884  mean  — the arithmetic mean is used for abundance_study_variable[n]
            m("study_variable[{}]-average_function".format(sv["index"]),
              "[MS, MS:1002884, mean, ]")
        # MS:1002886  standard deviation  — population stdev (divided by n) is used
        # for abundance_variation_study_variable[n].
        # Note: MS:1002885 (standard error = stdev/sqrt(n)) is a different statistic
        # and would be incorrect here.
        m("study_variable[{}]-variation_function".format(sv["index"]),
          "[MS, MS:1002886, standard deviation, ]")

    # MS:1000531  software  (correct parent term for a software entry)
    m("software[1]", "[MS, MS:1000531, GNPS2 FBMN, {}]".format(workflow_version))

    # MS:1001834  LC-MS label-free quantitation analysis
    # Correct usage: describes the *quantification workflow/method*, not a unit.
    m("quantification_method", "[MS, MS:1001834, LC-MS label-free quantitation analysis, ]")

    # For FBMN spectral library matches the appropriate MSI level is 2
    # (putatively annotated compounds). Using MS:1002896 with value "Level 2".
    m("small_molecule-identification_reliability",
      "[MS, MS:1002896, compound identification confidence level, Level 2]")

    # CV terms
    m("cv[1]-label", "MS")
    m("cv[1]-full_name", "PSI-MS controlled vocabulary")
    m("cv[1]-version", "4.1.238")
    m("cv[1]-uri",
      "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo")

    # database
    m("database[1]", "[, , GNPS spectral library, ]")
    m("database[1]-prefix", "gnps")
    m("database[1]-version", "Unknown")
    m("database[1]-uri", "https://gnps2.org/")

    # MS:1001844  MS1 feature area  — the correct term for LC-MS peak area
    # used as the quantification value in a label-free experiment.
    m("small_molecule-quantification_unit",
      "[MS, MS:1001844, MS1 feature area, ]")
    m("small_molecule_feature-quantification_unit",
      "[MS, MS:1001844, MS1 feature area, ]")

    # The correct term for the MQScore (cosine similarity) used in GNPS
    # library matching is:
    #   MS:1003304  spectral dot product  (synonym: cosine similarity)
    m("id_confidence_measure[1]",
      "[MS, MS:1003304, spectral dot product, ]")

    return lines


# ---------------------------------------------------------------------------
# SML section builder
# ---------------------------------------------------------------------------

SML_FIXED_COLS = [
    "SML_ID", "SMF_ID_REFS", "database_identifier", "chemical_formula",
    "smiles", "inchi", "chemical_name", "uri", "theoretical_neutral_mass",
    "adduct_ions", "reliability", "best_id_confidence_measure",
    "best_id_confidence_value",
]


def _get_abundance(row, assay, feature_row=None):
    """Get peak area abundance from the best available source."""
    for candidate in assay.get("abundance_candidates", [assay.get("abundance_key", assay.get("name", ""))]):
        col = candidate + " Peak area"
        # Try cluster summary row first
        v = row.get(col, None)
        if v is not None and not (isinstance(v, float) and math.isnan(v)):
            return _null(v)
        # Fallback to feature table row
        if feature_row is not None:
            v = feature_row.get(col, None)
            if v is not None and not (isinstance(v, float) and math.isnan(v)):
                return _null(v)
    return "null"


def build_sml(cluster_df: pd.DataFrame, lib_df: pd.DataFrame,
              feature_df: pd.DataFrame,
              assays: list, study_variables: list) -> list:
    lines = []

    # Build library lookup: cluster_index -> best library row
    lib_lookup = {}
    if not lib_df.empty and "#Scan#" in lib_df.columns:
        lb = lib_df.sort_values("MQScore", ascending=False) if "MQScore" in lib_df.columns else lib_df
        lb = lb.drop_duplicates(subset=["#Scan#"], keep="first")
        lib_lookup = lb.set_index("#Scan#").to_dict(orient="index")

    # Build feature table lookup by row ID
    feat_lookup = {}
    if not feature_df.empty and "row ID" in feature_df.columns:
        feat_lookup = feature_df.set_index("row ID").to_dict(orient="index")

    # Build SV abundance columns (mean per study_variable)
    sv_ab_cols = ["abundance_study_variable[{}]".format(j)
                  for j in range(1, len(study_variables) + 1)]
    sv_var_cols = ["abundance_variation_study_variable[{}]".format(j)
                   for j in range(1, len(study_variables) + 1)]

    # Assay abundance columns
    assay_ab_cols = ["abundance_assay[{}]".format(assay["index"]) for assay in assays]
    assay_by_index = {assay["index"]: assay for assay in assays}

    # Header
    header = SML_FIXED_COLS + assay_ab_cols + sv_ab_cols + sv_var_cols
    lines.append(_tab("SMH", *header))

    for _, row in cluster_df.iterrows():
        cluster_idx = int(row["cluster index"])
        sml_id = cluster_idx
        smf_ref = str(cluster_idx)

        lib = lib_lookup.get(cluster_idx, {})
        feat_row = feat_lookup.get(cluster_idx, None)

        db_id = _null(lib.get("SpectrumID", ""))
        formula = _null(lib.get("molecular_formula", ""))
        smiles = _null(lib.get("Smiles", ""))
        inchi = _null(lib.get("INCHI", ""))
        chem_name = _null(lib.get("Compound_Name", ""))
        uri = _null(lib.get("library_usi", ""))
        theor_neutral_mass = "null"
        adduct_ions = _format_adduct(_null(lib.get("Adduct", "")))

        # reliability: MSI level 2 for library-matched features, 3 for unmatched.
        # (mzTab-M reliability column is an integer 1–4 matching MSI levels)
        reliability = "2" if lib else "3"

        # CHANGED: best_id_confidence_measure now references MS:1003304
        # (spectral dot product / cosine similarity) instead of the wrong
        # MS:1002896 (compound identification confidence level).
        best_conf_measure = "[MS, MS:1003304, spectral dot product, ]" if lib else "null"
        best_conf_value = _null(lib.get("MQScore", "")) if lib else "null"

        # Per-assay abundances
        assay_vals = []
        for assay in assays:
            assay_vals.append(_get_abundance(row, assay, feat_row))

        # Per study_variable mean abundances
        sv_ab_vals = []
        sv_var_vals = []
        for sv in study_variables:
            vals = []
            for ai in sv["assay_refs"]:
                assay = assay_by_index.get(ai)
                if assay is None:
                    continue
                abundance = _get_abundance(row, assay, feat_row)
                if abundance != "null":
                    try:
                        vals.append(float(abundance))
                    except (ValueError, TypeError):
                        pass
            if vals:
                sv_ab_vals.append(str(sum(vals) / len(vals)))
                if len(vals) > 1:
                    mean = sum(vals) / len(vals)
                    var = math.sqrt(sum((x - mean) ** 2 for x in vals) / len(vals))
                    sv_var_vals.append(str(var))
                else:
                    sv_var_vals.append("null")
            else:
                sv_ab_vals.append("null")
                sv_var_vals.append("null")

        cells = [
            str(sml_id), smf_ref, db_id, formula, smiles, inchi, chem_name,
            uri, theor_neutral_mass, adduct_ions, reliability,
            best_conf_measure, best_conf_value,
        ] + assay_vals + sv_ab_vals + sv_var_vals

        lines.append(_tab("SML", *cells))

    return lines


# ---------------------------------------------------------------------------
# SMF section builder
# ---------------------------------------------------------------------------

SMF_FIXED_COLS = [
    "SMF_ID", "SME_ID_REFS", "SME_ID_REF_ambiguity_code",
    "adduct_ion", "isotopomer",
    "exp_mass_to_charge", "charge", "retention_time_in_seconds",
    "retention_time_in_seconds_start", "retention_time_in_seconds_end",
]


def build_smf(cluster_df: pd.DataFrame, lib_df: pd.DataFrame,
              feature_df: pd.DataFrame,
              assays: list) -> list:
    lines = []

    lib_lookup = {}
    if not lib_df.empty and "#Scan#" in lib_df.columns:
        lb = lib_df.sort_values("MQScore", ascending=False) if "MQScore" in lib_df.columns else lib_df
        lb = lb.drop_duplicates(subset=["#Scan#"], keep="first")
        lib_lookup = lb.set_index("#Scan#").to_dict(orient="index")

    feat_lookup = {}
    if not feature_df.empty and "row ID" in feature_df.columns:
        feat_lookup = feature_df.set_index("row ID").to_dict(orient="index")

    assay_ab_cols = ["abundance_assay[{}]".format(assay["index"]) for assay in assays]
    header = SMF_FIXED_COLS + assay_ab_cols
    lines.append(_tab("SFH", *header))

    for _, row in cluster_df.iterrows():
        cluster_idx = int(row["cluster index"])
        feat_row = feat_lookup.get(cluster_idx, None)

        # SME_ID_REFS: same as SMF_ID (one SME per feature)
        sme_ref = str(cluster_idx)

        lib = lib_lookup.get(cluster_idx, {})
        adduct = _format_adduct(_null(lib.get("Adduct", ""))) if lib else "null"

        mz = _null(row.get("parent mass", row.get("row m/z", None)))
        charge = "1"

        rt_raw = row.get("RTMean", row.get("row retention time", None))
        if rt_raw is not None and not (isinstance(rt_raw, float) and math.isnan(rt_raw)):
            try:
                rt_sec = float(rt_raw) * 60.0
                rt_str = str(rt_sec)
            except (ValueError, TypeError):
                rt_str = "null"
        else:
            rt_str = "null"

        assay_vals = []
        for assay in assays:
            assay_vals.append(_get_abundance(row, assay, feat_row))

        cells = [
            str(cluster_idx), sme_ref, "null",
            adduct, "null",
            mz, charge, rt_str,
            "null", "null",
        ] + assay_vals

        lines.append(_tab("SMF", *cells))

    return lines


# ---------------------------------------------------------------------------
# SME section builder
# ---------------------------------------------------------------------------

SME_FIXED_COLS = [
    "SME_ID", "evidence_input_id", "database_identifier", "chemical_formula",
    "smiles", "inchi", "chemical_name", "uri", "derivatized_form",
    "adduct_ion", "exp_mass_to_charge", "charge", "theoretical_mass_to_charge",
    "spectra_ref", "identification_method", "ms_level",
    "id_confidence_measure[1]", "rank",
]

# CHANGED: was MS:1003143 (mass array — a binary data array term, completely
# wrong context). The identification method should describe the *search
# approach*. MS:1001031 (spectral library search) is the correct term for
# GNPS library matching.
FIXED_IDENTIFICATION_METHOD = "[MS, MS:1001031, spectral library search, ]"

# MS:1000511  ms level  — correct, value 2 = MS2 (tandem MS spectra used
# for library matching)
FIXED_MS_LEVEL = "[MS, MS:1000511, ms level, 2]"

FIXED_RANK = "1"


def build_sme(cluster_df: pd.DataFrame, lib_df: pd.DataFrame,
              source_ctx: dict = None) -> list:
    lines = []
    source_ctx = source_ctx or {}

    lib_lookup = {}
    if not lib_df.empty and "#Scan#" in lib_df.columns:
        lb = lib_df.sort_values("MQScore", ascending=False) if "MQScore" in lib_df.columns else lib_df
        lb = lb.drop_duplicates(subset=["#Scan#"], keep="first")
        lib_lookup = lb.set_index("#Scan#").to_dict(orient="index")

    lines.append(_tab("SEH", *SME_FIXED_COLS))

    for _, row in cluster_df.iterrows():
        cluster_idx = int(row["cluster index"])
        lib = lib_lookup.get(cluster_idx, {})

        mz = _null(row.get("parent mass", row.get("row m/z", None)))

        rt_raw = row.get("RTMean", row.get("row retention time", None))
        rt_str = "null"
        if rt_raw is not None and not (isinstance(rt_raw, float) and math.isnan(rt_raw)):
            try:
                rt_str = str(float(rt_raw) * 60.0)
            except (ValueError, TypeError):
                rt_str = "null"

        evidence_input_id = "{}_{}m/z".format(rt_str, mz)

        # Preserve original spectra_ref when the input itself was mzTab-M.
        spectra_ref = source_ctx.get("smf_spectra_refs", {}).get(
            cluster_idx, "ms_run[1]:index={}".format(cluster_idx)
        )

        if lib:
            db_id = _null(lib.get("SpectrumID", ""))
            formula = _null(lib.get("molecular_formula", ""))
            smiles = _null(lib.get("Smiles", ""))
            inchi = _null(lib.get("INCHI", ""))
            chem_name = _null(lib.get("Compound_Name", ""))
            uri = _null(lib.get("library_usi", ""))
            adduct = _format_adduct(_null(lib.get("Adduct", "")))
            theor_mz = _null(lib.get("Precursor_MZ", ""))
            conf_val_1 = _null(lib.get("MQScore", ""))
        else:
            db_id = formula = smiles = inchi = chem_name = uri = "null"
            adduct = theor_mz = conf_val_1 = "null"

        cells = [
            str(cluster_idx),       # SME_ID
            evidence_input_id,      # evidence_input_id
            db_id,
            formula,
            smiles,
            inchi,
            chem_name,
            uri,
            "null",                 # derivatized_form (nullable)
            adduct,
            mz,                     # exp_mass_to_charge (non-nullable)
            "1",                    # charge (non-nullable)
            theor_mz if theor_mz != "null" else mz,  # theoretical_mass_to_charge (non-nullable proxy)
            spectra_ref,            # spectra_ref (non-nullable)
            FIXED_IDENTIFICATION_METHOD,
            FIXED_MS_LEVEL,
            conf_val_1,             # id_confidence_measure[1]
            FIXED_RANK,
        ]

        lines.append(_tab("SME", *cells))

    return lines


# ---------------------------------------------------------------------------
# main writer
# ---------------------------------------------------------------------------

def write_mztab(args):
    cluster_df = load_cluster_summary(args.cluster_summary)
    lib_df = load_library_results(args.library_results)
    metadata_df = load_metadata(args.metadata)
    feature_df = load_feature_table(args.feature_table)
    source_ctx = load_source_mztab_context(args.source_mztab)

    assays = build_assays(feature_df, cluster_df, source_ctx)
    if not assays:
        # Fallback: if no "Peak area" columns, use all non-standard columns
        skip = {"cluster index", "parent mass", "RTMean", "row m/z",
                "row retention time", "component", "Compound_Name", "#Scan#"}
        assays = [
            {
                "index": i + 1,
                "name": c,
                "sample_ref": "",
                "ms_run_ref": "ms_run[{}]".format(i + 1),
                "abundance_key": c,
                "abundance_candidates": [c],
            }
            for i, c in enumerate(cluster_df.columns)
            if c not in skip
        ]

    study_variables = build_study_variables(assays, metadata_df, source_ctx)

    workflow_version = _read_workflow_version(args.workflowinput_yaml)
    mtd_lines = build_mtd(assays, study_variables,
                          args.ms_run_file, args.featurefindingtool, workflow_version, source_ctx)
    sml_lines = build_sml(cluster_df, lib_df, feature_df, assays, study_variables)
    smf_lines = build_smf(cluster_df, lib_df, feature_df, assays)
    sme_lines = build_sme(cluster_df, lib_df, source_ctx)
    all_section_lines = mtd_lines + sml_lines + smf_lines + sme_lines
    max_cols = max((line.count("\t") + 1) for line in all_section_lines) if all_section_lines else 1
    section_separator = ("\t" * (max_cols - 1)) + "\n"

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write("COM\tmzTab-M file generated by FBMN workflow\n")
        fh.write("COM\tFeature finding tool: {}\n".format(args.featurefindingtool))
        fh.write("COM\tGenerated: {}\n".format(date.today().isoformat()))
        fh.write("\n")
        fh.writelines(mtd_lines)
        fh.write(section_separator)
        fh.writelines(sml_lines)
        fh.write(section_separator)
        fh.writelines(smf_lines)
        fh.write(section_separator)
        fh.writelines(sme_lines)

    print("[create_mztab_output] Written {} features to {}".format(
        len(cluster_df), args.output))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Create a valid mzTab-M 2.0 file from FBMN workflow outputs."
    )
    p.add_argument("--cluster_summary", required=True,
                   help="Enriched cluster summary TSV (from enrichClusterSummary)")
    p.add_argument("--library_results", default="NO_FILE",
                   help="Merged library results TSV (from librarygetGNPSAnnotations)")
    p.add_argument("--metadata", default="NO_FILE",
                   help="Merged metadata TSV (from createMetadataFile)")
    p.add_argument("--feature_table", default="NO_FILE",
                   help="Reformatted feature table CSV (from quantification_table_reformatted)")
    p.add_argument("--source_mztab", default="NO_FILE",
                   help="Original mzTab-M input file when featurefindingtool is MZTABM")
    p.add_argument("--ms_run_file", default="spectra.mgf",
                   help="Name / path of the spectra file used during networking")
    p.add_argument("--featurefindingtool", default="UNKNOWN",
                   help="Feature finding tool name")
    p.add_argument("--workflowinput_yaml", default="",
                   help="Path to workflowinput.yaml to read workflowversion")
    p.add_argument("--output", required=True,
                   help="Output mzTab file path")
    args = p.parse_args()
    write_mztab(args)


if __name__ == "__main__":
    main()
