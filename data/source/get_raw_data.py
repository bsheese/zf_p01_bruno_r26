#!/usr/bin/env python3
"""
get_raw_data.py — Download the raw distortion-corrected tracking data.

The raw frame-by-frame CSV (~153 MB, 4.85 M rows) is too large for GitHub and
is hosted externally. This script downloads it into data/source/ and verifies
the SHA-256 checksum, then the pipeline can run:

    python data/source/get_raw_data.py
    python src/01_ingest.py        # -> data/merged/tracking_merged.csv
    python src/03_features.py      # -> data/features/

The file is the Attempt-05 distortion-corrected Mathematica triallevel output
(columns: run, dayafternoon, tank, fishid, framenumber, x, y, deltax, deltay).
"""

import hashlib
import sys
import urllib.request
from pathlib import Path

# Zenodo deposit: https://zenodo.org/records/20707067  (CC BY 4.0)
RAW_DATA_URL = "https://zenodo.org/records/20707067/files/tracking_triallevel_distortion_corrected.csv?download=1"

DEST = Path(__file__).resolve().parent / "tracking_triallevel_distortion_corrected.csv"
SHA256 = "9f5c88414f86f93079fe100b166dc56ebe5c287dcefc839e4faddb1fb2747a4c"
SIZE_BYTES = 160_703_963


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if DEST.exists():
        print(f"Found existing {DEST.name}; verifying checksum …")
        if sha256_of(DEST) == SHA256:
            print("  ✓ checksum matches — nothing to do.")
            return 0
        print("  ✗ checksum mismatch — re-downloading.")

    if RAW_DATA_URL.startswith("REPLACE_ME"):
        print(
            "ERROR: RAW_DATA_URL is not set.\n"
            "Edit data/source/get_raw_data.py and paste the Zenodo/OSF direct\n"
            "download link, or obtain the file from the study authors.",
            file=sys.stderr,
        )
        return 1

    print(f"Downloading raw data ({SIZE_BYTES / 1e6:.0f} MB) …")
    print(f"  {RAW_DATA_URL}")
    tmp = DEST.with_suffix(".csv.part")
    urllib.request.urlretrieve(RAW_DATA_URL, tmp)

    print("Verifying checksum …")
    got = sha256_of(tmp)
    if got != SHA256:
        tmp.unlink(missing_ok=True)
        print(f"ERROR: checksum mismatch.\n  expected {SHA256}\n  got      {got}",
              file=sys.stderr)
        return 1

    tmp.rename(DEST)
    print(f"  ✓ saved → {DEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
