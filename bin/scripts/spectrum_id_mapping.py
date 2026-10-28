#!/usr/bin/env python3
"""Translate native executable spectrum positions back to original MGF IDs."""

import argparse
import csv
from pathlib import Path

from validate_feature_spectrum_alignment import mgf_headers


def index_mgf_ids(source, output):
    """Use positive query IDs inside the native library search only.

    The executable substitutes a position for scan zero. Give it unambiguous
    positive indices and retain an explicit map for restoring the results.
    Source/published MGF files keep their original SCANS identifiers.
    """
    mapping = {}
    index = 0
    with open(source) as incoming, open(output, "w") as outgoing:
        for line in incoming:
            if line.strip().upper() == "BEGIN IONS":
                index += 1
            if line.strip().upper().startswith("SCANS="):
                mapping[index] = int(line.split("=", 1)[1].strip())
                line = "SCANS={}\n".format(index)
            outgoing.write(line)
    if len(mapping) != index:
        raise ValueError("Every query spectrum must have a SCANS identifier")
    return mapping


def restore_networking_ids(spectra, pairs_path):
    # Native networking reports spectrum positions (1-based), even when SCANS
    # is present. For positive, padded inputs this map is the identity; it also
    # handles zero-based and sparse IDs without shifting their annotations.
    mapping = {index: int(headers["SCANS"])
               for index, headers in enumerate(mgf_headers(spectra), start=1)}
    pairs_path = Path(pairs_path)
    temporary = pairs_path.with_name(pairs_path.name + ".tmp")
    try:
        with pairs_path.open(newline="") as source, temporary.open("w", newline="") as output:
            reader = csv.DictReader(source, delimiter="\t")
            if not {"CLUSTERID1", "CLUSTERID2"}.issubset(reader.fieldnames or []):
                raise ValueError("Networking results lack CLUSTERID1 or CLUSTERID2")
            writer = csv.DictWriter(output, fieldnames=reader.fieldnames,
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                for column in ["CLUSTERID1", "CLUSTERID2"]:
                    row[column] = mapping[int(row[column])]
                writer.writerow(row)
        temporary.replace(pairs_path)
    finally:
        if temporary.exists():
            temporary.unlink()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectra")
    parser.add_argument("pairs")
    args = parser.parse_args()
    restore_networking_ids(args.spectra, args.pairs)
