#!/usr/bin/env python3
"""Normalize Everything Bagel FBMN exports without changing feature IDs."""

import mzmine_formatter
from validate_feature_spectrum_alignment import canonical_id, finite_number


def convert_to_feature_csv(input_filename, output_filename):
    # EB's feature table uses the MZmine2-style column contract.
    mzmine_formatter.convert_to_feature_csv(input_filename, output_filename)


def _is_padding(headers, has_fragment_lines):
    values = headers.get("PEPMASS", "").split()
    mz = finite_number(values[0]) if values else None
    title = headers.get("TITLE")
    return (
        mz == 0
        and not has_fragment_lines
        and not any(key in headers for key in (
            "FEATURE_ID", "SOURCE_FEATURE_ID", "SOURCE_FILE", "SOURCE_SCAN"
        ))
        and (title is None or title == "SCAN=" + headers.get("SCANS", ""))
    )


def convert_mgf(input_filename, output_filename):
    """Drop only empty zero-precursor padding; retain real blocks verbatim.

    Inspect one block at a time and preserve sparse/zero SCANS identifiers,
    rescue annotations, source references, and all fragment lines. An invalid
    real spectrum is left for the strict alignment check, rather than hidden.
    Returns (real_blocks_written, padding_blocks_removed).
    """
    seen = set()
    block = None
    written = removed = 0
    with open(input_filename, encoding="utf-8-sig", newline="") as source, \
            open(output_filename, "w", encoding="utf-8", newline="") as output:
        for line in source:
            text = line.strip()
            if text.upper() == "BEGIN IONS":
                if block is not None:
                    raise ValueError("Nested BEGIN IONS in Everything Bagel MGF")
                block = [line]
                headers = {}
                has_fragment_lines = False
            elif block is None:
                output.write(line)
            else:
                block.append(line)
                if text.upper() == "END IONS":
                    identifier = canonical_id(headers.get("SCANS", ""))
                    if not identifier:
                        raise ValueError("Everything Bagel spectrum has no SCANS identifier")
                    if identifier in seen:
                        raise ValueError("Duplicate Everything Bagel spectrum ID: " + identifier)
                    seen.add(identifier)
                    if _is_padding(headers, has_fragment_lines):
                        removed += 1
                    else:
                        output.writelines(block)
                        written += 1
                    block = None
                elif "=" in text:
                    key, value = text.split("=", 1)
                    headers[key.strip().upper()] = value.strip()
                elif text and not text.startswith(("#", ";", "!", "/")):
                    # Unexpected body text must not be mistaken for padding.
                    has_fragment_lines = True
        if block is not None:
            raise ValueError("Unterminated spectrum in Everything Bagel MGF")
    print("Everything Bagel: retained {} real spectra; removed {} empty padding spectra.".format(
        written, removed
    ))
    return written, removed
