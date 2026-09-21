#!/usr/bin/env python3
"""Report which external datasets are present in a GloCoFlood data root."""

from __future__ import annotations

import argparse
from pathlib import Path


EXPECTED = {
    "IBTrACS": "tracks/IBTrACS.ALL.v04r01.nc",
    "ERA5 TCR environment": "era5/tcr_environment",
    "ERA5 VIC forcing": "era5/vic_forcing",
    "ERA5-Land rainfall": "era5/era5_land_rainfall",
    "C15 lookup profiles": "C15_predata",
    "CaMa-Flood maps": "camaflood/map_v420/glb_15min",
    "HydroMT-SFINCS static inputs": "HydroMT-SFINCS/examples/data",
    "ADCIRC tide model": "tides",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_root", type=Path)
    args = parser.parse_args()
    root = args.data_root.expanduser().resolve()
    print("Data root:", root)
    present = 0
    for label, relative in EXPECTED.items():
        path = root / relative
        ok = path.exists() and (not path.is_file() or path.stat().st_size > 0)
        present += int(ok)
        print(f"[{'OK' if ok else 'MISSING':7}] {label}: {path}")
    print(f"Present: {present}/{len(EXPECTED)}")
    return 0 if present == len(EXPECTED) else 2


if __name__ == "__main__":
    raise SystemExit(main())
