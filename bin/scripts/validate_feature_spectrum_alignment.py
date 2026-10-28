#!/usr/bin/env python3
"""Stop FBMN when converted feature IDs no longer identify their spectra."""

import argparse
import csv
from decimal import Decimal, InvalidOperation
import math
from pathlib import Path
import sys


MZ_TOLERANCE_DA = 1.0
ERROR_PREFIX = "An ID alignment skew has been detected between "


class AlignmentSkewError(ValueError):
    pass


def canonical_id(value):
    """Treat integer CSV values such as 12 and 12.0 as the same ID."""
    value = value.strip()
    try:
        number = Decimal(value)
        if number.is_finite() and number == number.to_integral_value():
            return str(int(number))
    except InvalidOperation:
        pass
    return value


def finite_number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (ValueError, TypeError):
        return None


def read_features(path):
    features = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"row ID", "row m/z"}
        if not required.issubset(reader.fieldnames or []):
            raise AlignmentSkewError("The converted feature table lacks row ID or row m/z.")
        for row in reader:
            identifier = canonical_id(row["row ID"])
            mz = finite_number(row["row m/z"])
            if not identifier or mz is None or mz <= 0:
                raise AlignmentSkewError("The converted feature table has an invalid ID or precursor m/z.")
            if identifier in features:
                raise AlignmentSkewError("Feature ID {} occurs more than once.".format(identifier))
            features[identifier] = mz
    return features


def mgf_headers(path):
    """Stream MGF metadata without materializing fragment peaks."""
    headers = None
    with open(path, encoding="utf-8-sig") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.upper() == "BEGIN IONS":
                headers = {}
            elif line.upper() == "END IONS" and headers is not None:
                yield headers
                headers = None
            elif headers is not None and "=" in line:
                key, value = line.split("=", 1)
                headers[key.strip().upper()] = value.strip()


def validate_alignment(feature_path, spectrum_path):
    """Check spectra by ID; feature rows without MS/MS are allowed.

    Run before filter_spectra.py injects dummy scans. Compare precursor m/z,
    not acquisition RT, which need not equal an aligned feature's RT center.
    """
    features = read_features(feature_path)
    checked = 0
    seen = set()
    for spectrum in mgf_headers(spectrum_path):
        identifier = canonical_id(spectrum.get("SCANS", ""))
        if not identifier:
            raise AlignmentSkewError("An MS/MS spectrum has no SCANS identifier.")
        if identifier in seen:
            raise AlignmentSkewError("Spectrum ID {} occurs more than once.".format(identifier))
        seen.add(identifier)
        if identifier not in features:
            raise AlignmentSkewError(
                "Spectrum ID {} has no matching feature-table row ID.".format(identifier)
            )
        mz_values = spectrum.get("PEPMASS", "").split()
        mz = finite_number(mz_values[0]) if mz_values else None
        if mz is None or mz <= 0:
            raise AlignmentSkewError("Spectrum ID {} has no valid precursor m/z.".format(identifier))
        feature_mz = features[identifier]
        if abs(mz - feature_mz) > MZ_TOLERANCE_DA + 1e-9:
            raise AlignmentSkewError(
                "ID {}: spectrum precursor m/z {:.6f}, feature m/z {:.6f} "
                "(allowed difference: 1 Da).".format(identifier, mz, feature_mz)
            )
        checked += 1
    if not checked:
        raise AlignmentSkewError("No MS/MS spectra were available to check against the feature table.")
    return checked


def source_name(path):
    source = Path(path)
    if source.is_dir():
        files = [p for p in source.iterdir() if p.is_file()]
        if len(files) == 1:
            source = files[0]
    # Keep both the requested user guidance and the diagnostic within the
    # GNPS2 status page's 1000-character summary limit.
    return source.name[:160]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("feature_table")
    parser.add_argument("spectra")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--feature-source")
    parser.add_argument("--spectra-source")
    args = parser.parse_args()
    summary = Path(args.summary)
    try:
        checked = validate_alignment(args.feature_table, args.spectra)
    except AlignmentSkewError as error:
        message = (
            ERROR_PREFIX + 'feature table "{}" and MS/MS spectra "{}". '
            "If you believe this is not a result of unrelated input data, "
            "please contact the GNPS2 developers.\n\n{}\n"
        ).format(
            source_name(args.feature_source or args.feature_table),
            source_name(args.spectra_source or args.spectra),
            str(error)[:350],
        )
        # Failed Nextflow processes do not publish their declared outputs.
        # Write directly to the task's summary path before exiting nonzero.
        summary.parent.mkdir(parents=True, exist_ok=True)
        temporary = summary.with_name(summary.name + ".tmp")
        temporary.write_text(message, encoding="utf-8")
        temporary.replace(summary)
        print(message, file=sys.stderr, end="")
        return 1
    # Do not leave this check's old failure visible after a successful rerun.
    if summary.exists() and summary.read_text(encoding="utf-8").startswith(ERROR_PREFIX):
        summary.unlink()
    print("ID alignment check passed: {} spectra.".format(checked))
    return 0


if __name__ == "__main__":
    sys.exit(main())
