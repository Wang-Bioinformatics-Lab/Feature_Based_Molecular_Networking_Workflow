#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert MS-DIAL5 GNPS exports without changing their shared alignment IDs."""

import sys
import tempfile

import pandas as pd
from pyteomics import mgf


def convert_to_feature_csv(input_filename, output_filename):
    input_format = pd.read_csv(input_filename, sep="\t", skiprows=4)
    headers = list(input_format.columns)
    # GNPS exports place sample abundances after the last spectrum annotation.
    # Annotation counts vary between MS-DIAL5 versions and mobility exports.
    if "MS/MS spectrum" not in headers:
        raise ValueError("MS-DIAL5 GNPS feature table is missing the MS/MS spectrum column")
    sample_names = headers[headers.index("MS/MS spectrum") + 1:]

    output_format = input_format[["Alignment ID", "Average Mz", "Average Rt(min)"] + sample_names].copy()
    output_format.rename(columns={
        "Alignment ID": "row ID",
        "Average Mz": "row m/z",
        "Average Rt(min)": "row retention time",
        **{name: name + " Peak area" for name in sample_names},
    }, inplace=True)
    output_format.to_csv(output_filename, index=False)


def _filtered_spectra(reader):
    for spectrum in reader:
        # SCANS and feature-table Alignment ID already identify the same feature.
        # Neither a +1 shift nor an independent sequential index is appropriate.
        if "scans" not in spectrum["params"]:
            raise ValueError("MS-DIAL5 GNPS spectrum is missing its SCANS alignment ID")
        intensities = spectrum["intensity array"]
        if len(intensities):
            keep = intensities > max(intensities) * 0.01
            spectrum["m/z array"] = spectrum["m/z array"][keep]
            spectrum["intensity array"] = intensities[keep]
        yield spectrum


def convert_mgf(input_mgf, output_mgf):
    # Some exports contain MSP-style Num Peaks lines. Remove these without
    # creating a fixed-name temp.mgf that concurrent conversions could share.
    with tempfile.TemporaryFile(mode="w+") as cleaned:
        with open(input_mgf) as source:
            for line in source:
                if not line.lstrip().startswith("Num Peaks"):
                    cleaned.write(line)
        cleaned.seek(0)
        with mgf.read(cleaned, use_index=False) as reader, open(output_mgf, "w") as output:
            # Retain SCANS, TITLE, precursor, RT, charge, and other source metadata.
            # Keep the existing >1% base-peak filter and four-decimal peak output.
            mgf.write(_filtered_spectra(reader), output=output,
                      fragment_format="{:.4f} {:.4f}", use_numpy=False,
                      write_charges=False)


if __name__ == "__main__":
    convert_to_feature_csv(sys.argv[1], sys.argv[2])
