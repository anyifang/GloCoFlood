#!/usr/bin/env python3
"""Portable TOML-configured launcher for the four global SFINCS builders."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Sequence

try:
    import tomllib
except ModuleNotFoundError as exc:  # pragma: no cover - Python <3.11
    raise RuntimeError("Python 3.11 or newer is required for TOML configuration") from exc


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
RUNNER = SCRIPT_DIR / "run_four_compound_builds_parallel.py"

PATH_ENV = {
    "project_dir": "GLOCOFLOOD_PROJECT_DIR",
    "sfincs_model_root": "SFINCS_MODEL_ROOT",
    "adcirc_info_root": "ADCIRC_INFO_ROOT",
    "membership_csv": "SFINCS_MEMBERSHIP_CSV",
    "adcirc_mesh_root": "ADCIRC_MESH_ROOT",
    "adcirc_result_parent": "ADCIRC_RESULT_PARENT",
    "adcirc_run_parent": "ADCIRC_RUN_PARENT",
    "cama_historical_root": "CAMA_HISTORICAL_ROOT",
    "cama_future_root": "CAMA_FUTURE_ROOT",
    "cama_diagnostic_root": "CAMA_DIAGNOSTIC_ROOT",
    "cama_diagnostic_cache": "CAMA_DIAGNOSTIC_CACHE",
    "historical_track_nc": "TC_TRACK_HISTORICAL_NC",
    "era5_t600_root": "ERA5_T600_ROOT",
    "cmip6_t600_root": "CMIP6_T600_ROOT",
    "c15_predata_dir": "C15_PREDATA_DIR",
    "climada_tcr_data_dir": "CLIMADA_TCR_DATA_DIR",
    "vlm_adjusted_dep_root": "VLM_ADJUSTED_DEP_ROOT",
    "boundary_map_cache": "SFINCS_BOUNDARY_MAP_CACHE",
}

FILE_KEYS = {"membership_csv", "historical_track_nc"}
WRITABLE_KEYS = {"boundary_map_cache", "output_root"}


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run the global historical and SSP SFINCS build from TOML paths.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--check-only", action="store_true")
    p.add_argument("--return-period", type=int)
    p.add_argument("--cama-quantile", type=float)
    p.add_argument("--workers", type=int, choices=(4, 8))
    p.add_argument("--copy-mode", choices=("hardlink", "copy"))
    p.add_argument("--continue-on-error", action="store_true", default=None)
    p.add_argument("--overwrite", action="store_true", default=None)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--limit-events", type=int)
    p.add_argument("--blocks")
    p.add_argument("--models")
    p.add_argument("--events")
    return p


def load_config(path: Path) -> dict:
    with path.expanduser().resolve().open("rb") as stream:
        return tomllib.load(stream)


def configured(value, fallback):
    return fallback if value is None else value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    cfg = load_config(args.config)
    paths = cfg.get("paths", {})
    build = cfg.get("build", {})
    missing = sorted((set(PATH_ENV) | {"output_root"}) - set(paths))
    if missing:
        raise ValueError("Missing [paths] settings: " + ", ".join(missing))

    resolved: dict[str, Path] = {}
    for key, value in paths.items():
        path = Path(os.path.expandvars(str(value))).expanduser()
        # Keep the public configuration portable: relative entries are always
        # interpreted from the repository root, not from the caller's cwd.
        if not path.is_absolute():
            path = REPOSITORY_ROOT / path
        path = path.resolve()
        resolved[key] = path
        if key not in WRITABLE_KEYS:
            expected = path.is_file() if key in FILE_KEYS else path.is_dir()
            if not expected:
                kind = "file" if key in FILE_KEYS else "directory"
                raise FileNotFoundError(f"{key} {kind} does not exist: {path}")

    env = os.environ.copy()
    for key, variable in PATH_ENV.items():
        env[variable] = str(resolved[key])

    command = [sys.executable, str(RUNNER)]
    scalar_options = {
        "--return-period": configured(args.return_period, build.get("return_period", 100)),
        "--cama-quantile": configured(args.cama_quantile, build.get("cama_quantile", 0.50)),
        "--workers": configured(args.workers, build.get("workers", 8)),
        "--copy-mode": configured(args.copy_mode, build.get("copy_mode", "hardlink")),
        "--log-level": build.get("log_level", "INFO"),
        "--output-root": resolved["output_root"],
        "--membership-csv": resolved["membership_csv"],
        "--adcirc-mesh-root": resolved["adcirc_mesh_root"],
        "--boundary-map-cache-dir": resolved["boundary_map_cache"],
    }
    for option, value in scalar_options.items():
        command.extend([option, str(value)])
    for option, value in (
        ("--limit-events", args.limit_events),
        ("--blocks", args.blocks),
        ("--models", args.models),
        ("--events", args.events),
    ):
        if value is not None:
            command.extend([option, str(value)])

    if bool(configured(args.continue_on_error, build.get("continue_on_error", False))):
        command.append("--continue-on-error")
    if bool(configured(args.overwrite, build.get("overwrite", False))):
        command.append("--overwrite")
    if args.dry_run:
        command.append("--dry-run")

    print("Configuration is valid.")
    print("Output root:", resolved["output_root"])
    print("Command:", subprocess.list2cmdline(command))
    if args.check_only:
        return 0
    return subprocess.call(command, cwd=SCRIPT_DIR, env=env)


if __name__ == "__main__":
    raise SystemExit(main())
