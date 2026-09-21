#!/usr/bin/env python3
"""Download the official NOAA IBTrACS v04r01 global NetCDF archive."""

from __future__ import annotations

import argparse
import os
import urllib.request
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_URL = (
    "https://www.ncei.noaa.gov/data/"
    "international-best-track-archive-for-climate-stewardship-ibtracs/"
    "v04r01/access/netcdf/IBTrACS.ALL.v04r01.nc"
)
DEFAULT_OUTPUT = REPOSITORY_ROOT / "external" / "tracks" / "IBTrACS.ALL.v04r01.nc"


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--url", default=DEFAULT_URL)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    return p


def download(url: str, output: Path, overwrite: bool) -> None:
    output = output.expanduser().resolve()
    part = output.with_suffix(output.suffix + ".part")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file() and output.stat().st_size > 0 and not overwrite:
        print(f"Existing file retained: {output} ({output.stat().st_size:,} bytes)")
        return
    if overwrite:
        part.unlink(missing_ok=True)

    offset = part.stat().st_size if part.exists() else 0
    request = urllib.request.Request(url, headers={"User-Agent": "GloCoFlood-data-downloader/1.0"})
    if offset:
        request.add_header("Range", f"bytes={offset}-")
    with urllib.request.urlopen(request, timeout=120) as response:
        resumed = offset > 0 and getattr(response, "status", None) == 206
        mode = "ab" if resumed else "wb"
        if not resumed:
            offset = 0
        total_header = response.headers.get("Content-Length")
        total = offset + int(total_header) if total_header else None
        written = offset
        with part.open(mode) as stream:
            while True:
                chunk = response.read(8 * 1024 * 1024)
                if not chunk:
                    break
                stream.write(chunk)
                written += len(chunk)
                if total:
                    print(f"\r{written / 1e6:,.1f}/{total / 1e6:,.1f} MB", end="", flush=True)
                else:
                    print(f"\r{written / 1e6:,.1f} MB", end="", flush=True)
    print()
    if not part.is_file() or part.stat().st_size == 0:
        raise RuntimeError("Download produced an empty file")
    os.replace(part, output)
    print(f"Saved: {output} ({output.stat().st_size:,} bytes)")


def main() -> int:
    args = parser().parse_args()
    print("Source:", args.url)
    print("Target:", args.output.expanduser().resolve())
    if args.dry_run:
        return 0
    download(args.url, args.output, args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
