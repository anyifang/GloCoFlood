#!/usr/bin/env python3
"""Create or execute chunked CDS requests for GloCoFlood ERA5 inputs."""

from __future__ import annotations

import argparse
import calendar
import json
import os
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPOSITORY_ROOT / "external" / "era5"
MONTHS = tuple(f"{month:02d}" for month in range(1, 13))
SIX_HOURLY = ("00:00", "06:00", "12:00", "18:00")


def parse_months(value: str) -> tuple[str, ...]:
    if value.lower() == "all":
        return MONTHS
    result = tuple(f"{int(item):02d}" for item in value.split(","))
    if not result or any(int(item) not in range(1, 13) for item in result):
        raise argparse.ArgumentTypeError("months must be 'all' or comma-separated values 1-12")
    return result


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "kind",
        choices=("tcr-environment", "vic-forcing", "era5-land-rainfall"),
        help="Dataset/request template to use.",
    )
    p.add_argument("--start-year", type=int, required=True)
    p.add_argument("--end-year", type=int, required=True)
    p.add_argument("--months", type=parse_months, default=MONTHS)
    p.add_argument(
        "--area",
        type=float,
        nargs=4,
        metavar=("NORTH", "WEST", "SOUTH", "EAST"),
        default=(90.0, -180.0, -90.0, 180.0),
    )
    p.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--execute", action="store_true", help="Submit requests; default is a request preview.")
    p.add_argument("--overwrite", action="store_true")
    return p


def monthly_environment(year: int, months: tuple[str, ...], area: list[float]):
    dataset = "reanalysis-era5-pressure-levels-monthly-means"
    request = {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": ["temperature", "u_component_of_wind", "v_component_of_wind"],
        "pressure_level": ["250", "600", "850"],
        "year": [str(year)],
        "month": list(months),
        "time": ["00:00"],
        "area": area,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }
    # This exact annual filename is consumed by the historical SFINCS/TCR
    # forcing builder. The file also retains 250/850-hPa winds for track
    # preprocessing and diagnostics.
    return dataset, request, Path("tcr_environment") / f"era5_t600_monthly_{year}.nc"


def hourly_request(kind: str, year: int, month: str, area: list[float]):
    days = [f"{day:02d}" for day in range(1, calendar.monthrange(year, int(month))[1] + 1)]
    if kind == "vic-forcing":
        dataset = "reanalysis-era5-single-levels"
        variables = [
            "2m_temperature", "2m_dewpoint_temperature", "total_precipitation",
            "surface_pressure", "surface_solar_radiation_downwards",
            "surface_thermal_radiation_downwards", "10m_u_component_of_wind",
            "10m_v_component_of_wind",
        ]
        folder = "vic_forcing"
        stem = "era5_vic"
    else:
        dataset = "reanalysis-era5-land"
        variables = ["total_precipitation"]
        folder = "era5_land_rainfall"
        stem = "era5_land_tp"
    request = {
        "variable": variables,
        "year": [str(year)],
        "month": [month],
        "day": days,
        "time": list(SIX_HOURLY),
        "area": area,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }
    return dataset, request, Path(folder) / f"{stem}_{year}{month}.nc"


def iter_requests(args: argparse.Namespace):
    area = [float(value) for value in args.area]
    for year in range(args.start_year, args.end_year + 1):
        if args.kind == "tcr-environment":
            yield monthly_environment(year, args.months, area)
        else:
            for month in args.months:
                yield hourly_request(args.kind, year, month, area)


def main() -> int:
    args = build_parser().parse_args()
    if args.end_year < args.start_year:
        raise SystemExit("--end-year must be >= --start-year")
    output_root = args.output_root.expanduser().resolve()
    jobs = list(iter_requests(args))
    print(f"Requests: {len(jobs)}")
    print(f"Output  : {output_root}")
    if not args.execute:
        for dataset, request, relative in jobs[:3]:
            print(json.dumps({"dataset": dataset, "target": str(relative), "request": request}, indent=2))
        if len(jobs) > 3:
            print(f"... {len(jobs) - 3} additional request(s); add --execute to submit")
        return 0

    try:
        import cdsapi
    except ImportError as exc:
        raise SystemExit("Install cdsapi first: python -m pip install -r data_download/requirements.txt") from exc
    client = cdsapi.Client()
    for index, (dataset, request, relative) in enumerate(jobs, start=1):
        target = output_root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.is_file() and target.stat().st_size > 0 and not args.overwrite:
            print(f"[{index}/{len(jobs)}] skip existing {target}")
            continue
        part = target.with_suffix(target.suffix + ".part")
        if part.exists():
            part.unlink()
        print(f"[{index}/{len(jobs)}] {dataset} -> {target}")
        client.retrieve(dataset, request, str(part))
        if not part.is_file() or part.stat().st_size == 0:
            raise RuntimeError(f"CDS returned an empty file: {part}")
        os.replace(part, target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
