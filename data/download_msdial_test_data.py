#!/usr/bin/env python3
"""Download or verify the locally excluded MS-DIAL regression data."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "test/data/msdial_manifest.json"


def matches(path, entry):
    if not path.is_file() or path.stat().st_size != entry["bytes"]:
        return False
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest() == entry["sha256"]


def fetch(entry, destination, verify_only=False):
    path = destination / entry["path"]
    if matches(path, entry):
        print("Verified {}".format(entry["path"]))
        return
    if verify_only:
        raise ValueError("Missing or changed fixture: {}".format(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".download")
    try:
        with urlopen(entry["source_url"], timeout=120) as response, temporary.open("wb") as output:
            for chunk in iter(lambda: response.read(1024 * 1024), b""):
                output.write(chunk)
        if not matches(temporary, entry):
            raise ValueError("Checksum/size mismatch downloading {}".format(entry["source_url"]))
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()
    print("Downloaded {}".format(entry["path"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=ROOT / "data/msdial")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--case", action="append", help="Download only this case (repeatable).")
    args = parser.parse_args()
    manifest = json.loads(DEFAULT_MANIFEST.read_text())
    names = {case["name"] for case in manifest["cases"]}
    if args.case and not set(args.case).issubset(names):
        parser.error("Unknown case; choose from: " + ", ".join(sorted(names)))
    try:
        for case in manifest["cases"]:
            if args.case and case["name"] not in args.case:
                continue
            for entry in case["files"]:
                fetch(entry, args.destination, args.verify_only)
    except (OSError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 1
    if not args.verify_only:
        args.destination.mkdir(parents=True, exist_ok=True)
        (args.destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
