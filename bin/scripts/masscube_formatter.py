#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026-06-01

@purpose: convert MassCube aligned feature table + MSP spectra into the
          GNPS-standard quantification CSV and MGF consumed by FBMN.
"""
import sys
import pandas as pd

# In a MassCube aligned feature table, columns up to and including
# "database" are metadata; every column after it is a per-sample intensity.
LAST_METADATA_COLUMN = "database"
DEFAULT_METADATA_COLUMN_COUNT = 30


def convert_to_feature_csv(input_filename, output_filename):
    df = pd.read_csv(input_filename, sep="\t")
    columns = list(df.columns)

    if LAST_METADATA_COLUMN in columns:
        boundary = columns.index(LAST_METADATA_COLUMN) + 1
    else:
        print("Warning: '{}' column not found; falling back to fixed metadata "
              "column count {}".format(LAST_METADATA_COLUMN,
                                       DEFAULT_METADATA_COLUMN_COUNT))
        boundary = DEFAULT_METADATA_COLUMN_COUNT

    sample_names = columns[boundary:]

    output_records = []
    for record in df.to_dict(orient="records"):
        output_record = {
            "row ID": str(int(record["feature_ID"])),
            "row m/z": str(record["m/z"]),
            "row retention time": str(record["RT"]),
        }
        for sample_name in sample_names:
            output_record[sample_name + " Peak area"] = record[sample_name]
        output_records.append(output_record)

    output_headers = ["row ID", "row m/z", "row retention time"]
    output_headers += [sample_name + " Peak area" for sample_name in sample_names]

    output_df = pd.DataFrame(output_records)
    output_df.to_csv(output_filename, sep=",", index=False, columns=output_headers)
    return


def _write_block(out, scan, precursor_mz, peaks):
    """Write one MGF block. Returns True if written, False if skipped
    (empty spectrum or missing precursor)."""
    if scan is None or precursor_mz is None or len(peaks) == 0:
        return False
    out.write("BEGIN IONS\n")
    out.write("TITLE=SCAN={}\n".format(scan))
    out.write("SCANS={}\n".format(scan))
    out.write("PEPMASS={}\n".format(precursor_mz))
    for mz, intensity in peaks:
        out.write("{} {}\n".format(mz, intensity))
    out.write("END IONS\n\n")
    return True


def convert_mgf(input_msp, output_mgf):
    """Convert a MassCube MSP file to MGF. One block per feature that has at
    least one MS2 peak (Num Peaks > 0); empty spectra are skipped. SCANS is the
    MSP ID, which equals the feature table's feature_ID / row ID."""
    written = 0
    scan = None
    precursor_mz = None
    peaks = []
    reading_peaks = False

    with open(output_mgf, "w") as out:
        with open(input_msp) as fh:
            for raw_line in fh:
                line = raw_line.strip()

                if line == "":
                    # Blank line terminates a block
                    if _write_block(out, scan, precursor_mz, peaks):
                        written += 1
                    scan = None
                    precursor_mz = None
                    peaks = []
                    reading_peaks = False
                    continue

                if reading_peaks:
                    parts = line.split()
                    if len(parts) >= 2:
                        peaks.append((parts[0], parts[1]))
                    continue

                lower = line.lower()
                if lower.startswith("id:"):
                    scan = line.split(":", 1)[1].strip()
                elif lower.startswith("precursormz:"):
                    precursor_mz = line.split(":", 1)[1].strip()
                elif lower.startswith("num peaks:"):
                    reading_peaks = True

        # Flush the final block if the file doesn't end with a blank line
        if _write_block(out, scan, precursor_mz, peaks):
            written += 1

    return written


if __name__ == "__main__":
    convert_to_feature_csv(sys.argv[1], sys.argv[2])
