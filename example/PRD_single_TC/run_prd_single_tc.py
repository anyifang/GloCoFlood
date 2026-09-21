#!/usr/bin/env python3
"""Build and verify one historical Pearl River Delta (PRD) TC case."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import xarray as xr
except ModuleNotFoundError as exc:
    raise SystemExit(
        "The PRD verifier requires xarray and netCDF4. Install them with: "
        "python -m pip install -r example/PRD_single_TC/requirements.txt"
    ) from exc


EXAMPLE_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = EXAMPLE_DIR.parents[1]
BUILD_DIR = REPOSITORY_ROOT / "SFINCS" / "global_build"
if str(BUILD_DIR) not in sys.path:
    sys.path.insert(0, str(BUILD_DIR))

DEFAULT_CONFIG = BUILD_DIR / "config" / "paths.toml"
BUILDER = BUILD_DIR / "P2_build_sfincs_historical_return_period_compound_forcing.py"
BUNDLED_CASE = EXAMPLE_DIR / "data" / "bundled_case" / "GTC_0009"
BUNDLED_MANIFEST = EXAMPLE_DIR / "data" / "MANIFEST.csv"
PRD_BLOCK = "ADC_WNP_04"
PRD_MODEL = "GTC_0009"
PRD_EVENT = "track_000261_tracks_GL_era5_197501_201412_000261"
EXPECTED_FILES = (
    "sfincs.inp",
    "sfincs.bnd",
    "sfincs.bzs",
    "sfincs.dep",
    "sfincs.dis",
    "sfincs.ind",
    "sfincs.man",
    "sfincs.msk",
    "sfincs.src",
    "sfincs_precipitation.nc",
    "sfincs_pressure.nc",
    "sfincs_wind.nc",
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate one reproducible historical TC compound-forcing case for "
            "the Pearl River Delta SFINCS domain."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-root", type=Path, default=EXAMPLE_DIR / "output")
    parser.add_argument("--return-period", type=int, default=None)
    parser.add_argument("--cama-quantile", type=float, default=None)
    parser.add_argument("--copy-mode", choices=("copy", "hardlink"), default="copy")
    parser.add_argument(
        "--full-build",
        action="store_true",
        help=(
            "Rebuild forcing from the external ADCIRC, CaMa-Flood, track, "
            "ERA5 and C15 datasets configured in paths.toml. By default the "
            "repository's compact, ready-to-run PRD case is installed."
        ),
    )
    parser.add_argument(
        "--sfincs-executable",
        type=Path,
        default=None,
        help="Optionally run this SFINCS executable after preparing the case.",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def resolve_paths(config_path: Path) -> tuple[dict, dict[str, Path]]:
    # Keep the lightweight bundled example independent of the production
    # configuration module and its full dependency chain.
    from run_global_build_from_config import PATH_ENV, load_config

    config = load_config(config_path)
    raw_paths = config.get("paths", {})
    required = set(PATH_ENV) | {"output_root"}
    missing = sorted(required - set(raw_paths))
    if missing:
        raise ValueError("Missing [paths] settings: " + ", ".join(missing))
    paths = {
        key: Path(os.path.expandvars(str(value))).expanduser().resolve()
        for key, value in raw_paths.items()
    }
    return config, paths


def bundled_case_dir(output_root: Path) -> Path:
    return output_root / "historical" / PRD_BLOCK / PRD_EVENT / PRD_MODEL


def verify_bundled_manifest() -> None:
    with BUNDLED_MANIFEST.open("r", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            path = EXAMPLE_DIR / "data" / row["file"]
            if not path.is_file():
                raise FileNotFoundError(f"Bundled example file is missing: {path}")
            if path.stat().st_size != int(row["size_bytes"]):
                raise RuntimeError(f"Bundled example size mismatch: {path}")
            digest = hashlib.sha256()
            with path.open("rb") as source:
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    digest.update(chunk)
            if digest.hexdigest() != row["sha256"]:
                raise RuntimeError(f"Bundled example checksum mismatch: {path}")


def install_bundled_case(output_root: Path, overwrite: bool, dry_run: bool) -> Path:
    """Copy the compact published case without requiring global source data."""
    verify_bundled_manifest()
    verify_case(BUNDLED_CASE)
    case_dir = bundled_case_dir(output_root)
    if dry_run:
        print("Bundled data verified:", BUNDLED_CASE)
        return case_dir
    if case_dir.exists():
        if not overwrite:
            verify_case(case_dir)
            print("Existing PRD example verified:", case_dir)
            return case_dir
        shutil.rmtree(case_dir)
    case_dir.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(BUNDLED_CASE, case_dir, copy_function=shutil.copy2)
    verify_case(case_dir)
    print("Bundled PRD example installed and verified:", case_dir)
    return case_dir


def run_sfincs(executable: Path, case_dir: Path) -> None:
    executable = executable.expanduser().resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"SFINCS executable not found: {executable}")
    print("Running SFINCS:", executable)
    completed = subprocess.run([str(executable)], cwd=case_dir, check=False)
    if completed.returncode:
        raise RuntimeError(f"SFINCS exited with code {completed.returncode}")


def preflight(paths: dict[str, Path], return_period: int) -> None:
    model_dir = paths["sfincs_model_root"] / PRD_MODEL
    result_dir = paths["adcirc_result_parent"] / f"{return_period}yr" / PRD_BLOCK / PRD_EVENT
    run_dir = paths["adcirc_run_parent"] / f"{return_period}yr" / PRD_BLOCK / PRD_EVENT

    checks = {
        "PRD base model": model_dir / "sfincs.inp",
        "ADCIRC water level": result_dir / "fort.63",
        "ADCIRC TC forcing": run_dir / "fort.22",
        "ADCIRC TC metadata": run_dir / "fort22_meta.txt",
    }
    missing = [f"{label}: {path}" for label, path in checks.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError("PRD example inputs are incomplete:\n  " + "\n  ".join(missing))

    membership = paths["membership_csv"]
    with membership.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = csv.DictReader(stream)
        mapped = any(
            row.get("block_id") == PRD_BLOCK and row.get("model_domain_id") == PRD_MODEL
            for row in rows
        )
    if not mapped:
        raise ValueError(f"{PRD_BLOCK}/{PRD_MODEL} is absent from {membership}")


def verify_case(case_dir: Path) -> None:
    missing = [name for name in EXPECTED_FILES if not (case_dir / name).is_file()]
    empty = [name for name in EXPECTED_FILES if (case_dir / name).is_file() and not (case_dir / name).stat().st_size]
    if missing or empty:
        raise RuntimeError(f"Incomplete output; missing={missing}, empty={empty}")

    expected_variables = {
        "sfincs_precipitation.nc": {"Precipitation"},
        "sfincs_wind.nc": {"eastward_wind", "northward_wind"},
        "sfincs_pressure.nc": {"barometric_pressure"},
    }
    for filename, expected in expected_variables.items():
        with xr.open_dataset(case_dir / filename) as dataset:
            absent = expected - set(dataset.data_vars)
            if absent or dataset.sizes.get("time", 0) < 2:
                raise RuntimeError(
                    f"Invalid {filename}: absent variables={sorted(absent)}, "
                    f"time steps={dataset.sizes.get('time', 0)}"
                )


def main() -> int:
    args = build_parser().parse_args()
    output_root = args.output_root.expanduser().resolve()

    if not args.full_build:
        print("Mode        : bundled ready-to-run case")
        print("PRD domain  :", PRD_MODEL)
        print("ADCIRC block:", PRD_BLOCK)
        print("TC event    :", PRD_EVENT)
        case_dir = install_bundled_case(output_root, args.overwrite, args.dry_run)
        if args.sfincs_executable is not None and not args.dry_run:
            run_sfincs(args.sfincs_executable, case_dir)
        return 0

    config_path = args.config.expanduser().resolve()
    config, paths = resolve_paths(config_path)
    from run_global_build_from_config import PATH_ENV

    build = config.get("build", {})
    return_period = args.return_period or int(build.get("return_period", 100))
    cama_quantile = (
        args.cama_quantile
        if args.cama_quantile is not None
        else float(build.get("cama_quantile", 0.50))
    )
    preflight(paths, return_period)
    env = os.environ.copy()
    for key, variable in PATH_ENV.items():
        env[variable] = str(paths[key])

    command = [
        sys.executable,
        str(BUILDER),
        "--return-period",
        str(return_period),
        "--cama-quantile",
        str(cama_quantile),
        "--cama-quantile-mode",
        "tc_inlet_p50",
        "--blocks",
        PRD_BLOCK,
        "--models",
        PRD_MODEL,
        "--events",
        PRD_EVENT,
        "--limit-events",
        "1",
        "--copy-mode",
        args.copy_mode,
        "--output-root",
        str(output_root),
        "--membership-csv",
        str(paths["membership_csv"]),
        "--adcirc-mesh-root",
        str(paths["adcirc_mesh_root"]),
        "--boundary-map-cache-dir",
        str(paths["boundary_map_cache"]),
        "--rain-model",
        "tcr",
        "--tcr-wind-model",
        "c15",
        "--c15-rmax-out-of-range",
        "rescale",
    ]
    if args.dry_run:
        command.append("--dry-run")
    if args.overwrite:
        command.append("--overwrite")

    print("PRD domain :", PRD_MODEL)
    print("ADCIRC block:", PRD_BLOCK)
    print("TC event    :", PRD_EVENT)
    print("Command     :", subprocess.list2cmdline(command))
    completed = subprocess.run(command, cwd=BUILD_DIR, env=env, check=False)
    if completed.returncode:
        return completed.returncode
    if args.dry_run:
        print("Dry-run completed; no case was written.")
        return 0

    case_dir = output_root / "historical" / PRD_BLOCK / PRD_EVENT / PRD_MODEL
    verify_case(case_dir)
    print("PRD example generated and verified:", case_dir)
    if args.sfincs_executable is not None:
        run_sfincs(args.sfincs_executable, case_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
