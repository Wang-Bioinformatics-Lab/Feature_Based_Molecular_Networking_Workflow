#!/usr/bin/env python3
"""Run LibrarySearch_Rust and preserve FBMN's legacy result-table contract."""

import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile


LEGACY_COLUMNS = [
    "#Scan#",
    "SpectrumFile",
    "Annotation",
    "OrigAnnotation",
    "Protein",
    "dbIndex",
    "numMods",
    "matchOrientation",
    "startMass",
    "Charge",
    "MQScore",
    "p-value",
    "isDecoy",
    "StrictEnvelopeScore",
    "UnstrictEvelopeScore",
    "CompoundName",
    "Organism",
    "FileScanUniqueID",
    "FDR",
    "LibraryName",
    "mzErrorPPM",
    "LibMetaData",
    "Smiles",
    "Inchi",
    "LibSearchSharedPeaks",
    "Abundance",
    "ParentMassDiff",
    "SpecMZ",
    "ExactMass",
    "LibrarySpectrumID",
]


def _parse_charge(value):
    value = str(value or "").strip()
    if not value:
        return 0
    token = value.split()[0].rstrip("+-")
    try:
        charge = int(float(token))
    except ValueError:
        return 0
    return -abs(charge) if value.endswith("-") else charge


def read_query_metadata(mgf_path):
    """Return scan -> charge/RT/TIC without changing the query peak list."""
    spectra = {}
    headers = None
    tic = 0.0
    ordinal = 0
    with open(mgf_path, encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "BEGIN IONS":
                headers = {}
                tic = 0.0
                ordinal += 1
            elif line == "END IONS" and headers is not None:
                scan_text = headers.get("SCANS", headers.get("FEATURE_ID", ordinal))
                try:
                    scan = int(float(scan_text))
                except ValueError:
                    scan = ordinal
                rt = headers.get("RTINSECONDS", headers.get("RT", "0"))
                spectra[scan] = {
                    "charge": _parse_charge(headers.get("CHARGE", "0")),
                    "rt": rt,
                    "tic": tic,
                }
                headers = None
            elif headers is not None and "=" in line and not line[0:1].isdigit():
                key, value = line.split("=", 1)
                headers[key.upper()] = value
            elif headers is not None and line:
                fields = line.split()
                if len(fields) >= 2:
                    try:
                        tic += float(fields[1])
                    except ValueError:
                        pass
    return spectra


def _first(row, *keys, default=""):
    for key in keys:
        value = row.get(key)
        if value is not None and str(value) != "":
            return value
    return default


def convert_results(rust_tsv, output_tsv, query_mgf):
    query_metadata = read_query_metadata(query_mgf)
    query_basename = os.path.basename(query_mgf)

    with open(rust_tsv, newline="", encoding="utf-8") as src:
        rows = list(csv.DictReader(src, delimiter="\t"))

    with open(output_tsv, "w", newline="", encoding="utf-8") as dst:
        writer = csv.DictWriter(dst, fieldnames=LEGACY_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            scan = int(float(row["query_scan"]))
            query_mz = float(row["query_precursor_mz"])
            library_mz = float(row["library_precursor_mz"])
            delta_mz = float(row["delta_mz"])
            metadata = query_metadata.get(scan, {"charge": 0, "rt": "0", "tic": 0.0})
            spectrum_id = _first(row, "SPECTRUMID", "SpectrumID", "TITLE", "USI")
            if str(spectrum_id).startswith("mzspec:"):
                spectrum_id = str(spectrum_id).rsplit(":", 1)[-1]
            mz_error_ppm = abs(delta_mz) / library_mz * 1e6 if library_mz else math.nan

            writer.writerow(
                {
                    "#Scan#": scan,
                    "SpectrumFile": query_basename,
                    "Annotation": _first(row, "SEQ", default="*..*"),
                    "OrigAnnotation": "",
                    "Protein": f":{row['library_scan']}",
                    "dbIndex": row["library_idx"],
                    "numMods": -1,
                    "matchOrientation": 0,
                    "startMass": library_mz,
                    "Charge": metadata["charge"],
                    "MQScore": row["cosine"],
                    "p-value": metadata["rt"],
                    "isDecoy": 0,
                    "StrictEnvelopeScore": 0,
                    "UnstrictEvelopeScore": metadata["tic"],
                    "CompoundName": _first(row, "NAME", "COMPOUND_NAME"),
                    "Organism": _first(row, "ORGANISM"),
                    "FileScanUniqueID": f"{query_basename}{scan}",
                    "FDR": -1,
                    "LibraryName": os.path.basename(row["library_filename"]),
                    "mzErrorPPM": mz_error_ppm,
                    "LibMetaData": "",
                    "Smiles": _first(row, "SMILES"),
                    "Inchi": _first(row, "INCHI"),
                    "LibSearchSharedPeaks": row["matched_peaks"],
                    "Abundance": 0,
                    "ParentMassDiff": delta_mz,
                    "SpecMZ": query_mz,
                    "ExactMass": _first(row, "EXACTMASS", default=0),
                    "LibrarySpectrumID": spectrum_id,
                }
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("spectrum_file")
    parser.add_argument("library_file")
    parser.add_argument("result_folder")
    parser.add_argument("library_search_binary")
    parser.add_argument("--pm_tolerance", type=float, default=0.5)
    parser.add_argument("--fragment_tolerance", type=float, default=0.5)
    parser.add_argument("--library_min_cosine", type=float, default=0.7)
    parser.add_argument("--library_min_matched_peaks", type=int, default=6)
    parser.add_argument("--topk", type=int, default=1)
    parser.add_argument("--analog_search", type=int, default=0)
    parser.add_argument("--library_analog_max_shift", type=float, default=1999)
    parser.add_argument("--filter_precursor", type=int, choices=(0, 1), default=1)
    parser.add_argument("--filter_window", type=int, choices=(0, 1), default=1)
    args = parser.parse_args()

    output_dir = Path(args.result_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / (
        f"{Path(args.spectrum_file).name}_{Path(args.library_file).name}.tsv"
    )

    with tempfile.TemporaryDirectory(prefix="fbmn-library-search-") as temp_dir:
        rust_output = Path(temp_dir) / "rust_results.tsv"
        command = [
            args.library_search_binary,
            "-q",
            args.spectrum_file,
            "-l",
            args.library_file,
            "--pm-tolerance",
            str(args.pm_tolerance),
            "--tolerance",
            str(args.fragment_tolerance),
            "--threshold",
            str(args.library_min_cosine),
            "--min-matches",
            str(args.library_min_matched_peaks),
            "--top-k",
            str(args.topk),
            "--no-cache",
            "-o",
            str(rust_output),
        ]
        if args.analog_search:
            command.extend(["--analog", "--max-analog-delta", str(args.library_analog_max_shift)])
        if not args.filter_precursor:
            command.append("--no-filter-precursor")
        if not args.filter_window:
            command.append("--no-filter-window")

        subprocess.run(command, check=True)
        convert_results(rust_output, output_path, args.spectrum_file)


if __name__ == "__main__":
    main()
