#!/usr/bin/env python3
"""Build event-specific SFINCS forcing for ERA5 return-period TC events.

The script couples four existing data sources without rebuilding the SFINCS
grids:

* ADCIRC fort.14 + fort.63 -> sfincs.bzs
* Independent per-inlet TC-month-weighted CaMa-Flood P50 -> sfincs.dis
* ADCIRC fort.22 -> spatial wind and pressure NetCDF forcing
* CLIMADA-Petals TCR using a C15 profile generated directly from the
  translation-free gradient-wind peak (or R-CLIPER)
  -> spatial precipitation NetCDF forcing
* ERA5 monthly 600-hPa temperature -> TCR diagnoses 950-hPa saturation
  specific humidity at every track time using the track's total surface
  cyclone Vmax required by CLIMADA-Petals

The original model directories are never modified.  Each event/domain gets a
stand-alone runnable directory below ``--output-root``.  Static model files are
hard-linked by default and every mutable forcing/config file is newly written.

Water-level time convention
---------------------------
The SFINCS run spans WindowStart..WindowEnd (impact -3 d .. impact +2 d).
The first 2.5 d are an explicit zero-water-level buffer.  ADCIRC water levels
are then introduced with a six-hour linear ramp.  Native fort.63 timestamps
are used; no first/last-value extrapolation is allowed.

CaMa-Flood inflow convention
----------------------------
By default, every inlet cell is evaluated independently with the same observed
TC-landfall-month weighted P50 algorithm as ``P1_plot_sfincs_inflow_diagnostics.py``.
Multiple SFINCS source points forming one inlet section share that inlet flow
equally. GTC-total summary and other legacy modes remain available explicitly.

Examples
--------
Inventory only (no writes and no 53-GiB CaMa scan)::


Full production build (requires CLIMADA Core/Petals 6.0.x)::
    conda activate sfincs_tcr
    python P2_build_sfincs_historical_return_period_compound_forcing.py --return-period 200 --cama-quantile-mode tc_inlet_p90 --cama-quantile 0.90

"""

from __future__ import annotations

import argparse
import calendar
import hashlib
import importlib.util
import io
import json
import logging
import math
import os
import re
import shutil
import sys
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

# Windows conda does not always propagate these activation variables when a
# script is launched by absolute python.exe path. They must be set before
# importing GeoPandas/GDAL-backed packages; setting them later can leave GDAL
# partially initialized and cause a native process exit during TCR raster IO.
_CONDA_SHARE = Path(sys.prefix) / "Library" / "share"
_CONDA_DLL_HANDLES: list[object] = []
_CONDA_BIN = Path(sys.prefix) / "Library" / "bin"
if os.name == "nt" and _CONDA_BIN.is_dir():
    os.environ["PATH"] = str(_CONDA_BIN) + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        _CONDA_DLL_HANDLES.append(os.add_dll_directory(str(_CONDA_BIN)))
if not os.environ.get("GDAL_DATA") and (_CONDA_SHARE / "gdal").is_dir():
    os.environ["GDAL_DATA"] = str(_CONDA_SHARE / "gdal")
if not os.environ.get("PROJ_LIB") and (_CONDA_SHARE / "proj").is_dir():
    os.environ["PROJ_LIB"] = str(_CONDA_SHARE / "proj")

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS, Transformer
from scipy.spatial import cKDTree


LOG = logging.getLogger("sfincs_compound_forcing")
EARTH_RADIUS_KM = 6371.0088
FORT63_INVALID_ABS = 1.0e4
ADCIRC_LOCAL_SEARCH_EXTRA_KM = 1.0
CAMA_FILL_LOW = -9000.0
CAMA_FILL_HIGH = 1.0e19

SCRIPT_DIR = Path(__file__).resolve().parent


def configured_path(name: str, default: str | Path) -> Path:
    """Return an environment-configured path without embedding user paths."""
    return Path(os.environ.get(name, str(default))).expanduser()


PROJECT_DIR = Path(
    os.environ.get(
        "GLOCOFLOOD_PROJECT_DIR",
        os.environ.get("GLOBAL_FLOOD_PROJECT_DIR", SCRIPT_DIR.parent.parent),
    )
)
DEFAULT_MODEL_ROOT = configured_path(
    "SFINCS_MODEL_ROOT",
    PROJECT_DIR / "external" / "sfincs_models" / "global_sfincs_partition_models_cama_15min",
)
DEFAULT_ADCIRC_INFO_ROOT = configured_path(
    "ADCIRC_INFO_ROOT", PROJECT_DIR / "ADCIRC" / "global_build" / "catalog"
)
DEFAULT_MEMBERSHIP = configured_path(
    "SFINCS_MEMBERSHIP_CSV", DEFAULT_ADCIRC_INFO_ROOT / "global_tc_adcirc_block_members.csv"
)
DEFAULT_MESH_ROOT = configured_path(
    "ADCIRC_MESH_ROOT", PROJECT_DIR / "external" / "adcirc" / "adcirc_fort14_meshes"
)
DEFAULT_BOUNDARY_MAP_CACHE_DIR = configured_path(
    "SFINCS_BOUNDARY_MAP_CACHE", SCRIPT_DIR / "cache" / "adcirc_sfincs_boundary_maps"
)
DEFAULT_RETURN_PERIOD = 100
DEFAULT_ADCIRC_RESULT_PARENT = configured_path(
    "ADCIRC_RESULT_PARENT", PROJECT_DIR / "external" / "adcirc" / "return_out"
)
DEFAULT_ADCIRC_RUN_PARENT = configured_path(
    "ADCIRC_RUN_PARENT",
    PROJECT_DIR / "external" / "adcirc" / "return_period_tc_event_reruns",
)
DEFAULT_CAMA_ROOT = configured_path(
    "CAMA_HISTORICAL_ROOT",
    PROJECT_DIR / "external" / "camaflood" / "output" / "ERA5",
)
DEFAULT_CAMA_DIAGNOSTIC_ROOT = configured_path(
    "CAMA_DIAGNOSTIC_ROOT",
    PROJECT_DIR / "external" / "camaflood" / "output" / "gtc_cama_inflow_figures",
)
DEFAULT_CAMA_DIAGNOSTIC_CACHE_DIR = DEFAULT_CAMA_DIAGNOSTIC_ROOT / "cache"
DEFAULT_CAMA_MONTH_WEIGHTS_CSV = (
    DEFAULT_CAMA_DIAGNOSTIC_ROOT / "tables" / "gtc_tc_month_weights.csv"
)
DEFAULT_CAMA_MONTHLY_INLET_CSV = (
    DEFAULT_CAMA_DIAGNOSTIC_ROOT
    / "tables"
    / "gtc_monthly_climatology_by_inlet_cell.csv"
)
DEFAULT_CAMA_FLOW_SUMMARY_CSV = (
    DEFAULT_CAMA_DIAGNOSTIC_ROOT / "tables" / "gtc_flow_summary_by_case.csv"
)
DEFAULT_CAMA_INLET_CELLS_CSV = (
    DEFAULT_CAMA_DIAGNOSTIC_ROOT / "tables" / "gtc_cama_inlet_cells.csv"
)
DEFAULT_TRACK_NC = configured_path(
    "TC_TRACK_HISTORICAL_NC",
    PROJECT_DIR / "external" / "tracks" / "tracks_GL_era5_197501_201412.nc",
)
DEFAULT_ERA5_T600_ROOT = configured_path(
    "ERA5_T600_ROOT", PROJECT_DIR / "external" / "era5" / "tcr_environment"
)
DEFAULT_C15_PREDATA_DIR = configured_path(
    "C15_PREDATA_DIR", PROJECT_DIR / "external" / "C15_predata"
)
C15_TCR_MODEL_ID = 91015
C15_PROFILE_RE = re.compile(r"Wind_C15_data_Vmax(\d+)_Rmax(\d+)\.mat$", re.IGNORECASE)
C15_RESCALE_WARNED: set[tuple[int, int, int]] = set()
C15_MIN_LOOKUP_WIND_MS = 15.0
C15_MAX_LOOKUP_WIND_MS = 120.0
C15_MIN_RMAX_M = 30_000.0
C15_AIR_DENSITY_KG_M3 = 1.15
C15_ENVIRONMENTAL_PRESSURE_HPA = 1013.0
C15_HOLLAND_B_MIN = 1.0
C15_HOLLAND_B_MAX = 2.5
C15_BLEND_INNER_RADIUS_M = 500_000.0
C15_BLEND_OUTER_RADIUS_M = 700_000.0
C15_BETA_INNER_BASE_DEG = 10.0
C15_BETA_INNER_SLOPE_DEG = 10.0
C15_BETA_MID_BASE_DEG = 20.0
C15_BETA_MID_SLOPE_DEG = 25.0
C15_BETA_OUTER_DEG = 25.0
C15_BETA_MID_RADIUS_FACTOR = 1.2
# ``vmax_trks`` is the unambiguous near-surface total maximum wind supplied by
# the downscaled track. CLIMADA-Petals expects that surface value for both the
# public ``max_sustained_wind`` field and the saturation-humidity diagnosis.
# Its C15 dynamics, however, require a translation-free gradient wind. Follow
# CLIMADA's surface-to-gradient conversion exactly once.  The 0.9 factor is
# shared with the q950 thermodynamic diagnosis while Rmax and Holland B retain
# the ADCIRC/fort.22 surface-wind definitions:
#     v_surface_rotational = max(vmax_surface - v_translation, 0)
#     v_gradient = v_surface_rotational / 0.9
# The C15 field subsequently adds the radius-decaying translation vector once.
# Do not read ``v_trks`` here: its level/averaging semantics are not used by
# this build workflow.
TCR_GRADIENT_TO_SURFACE_WIND_FACTOR = 0.9
TCR_ENVIRONMENTAL_WIND_LEVELS_HPA = (250, 850)
TCR_REQUIRED_TRACK_WIND_VARIABLES = (
    "vmax_trks",
    "u250_trks",
    "v250_trks",
    "u850_trks",
    "v850_trks",
)


@dataclass(frozen=True)
class Event:
    block_id: str
    event_id: str
    result_dir: Path
    run_dir: Path
    fort63: Path
    fort14: Path
    fort22: Path
    meta_path: Path
    meta: dict[str, str]
    start: pd.Timestamp
    stop: pd.Timestamp
    impact: pd.Timestamp | None
    dt_seconds: int

    @property
    def duration_seconds(self) -> int:
        return int(round((self.stop - self.start).total_seconds()))


@dataclass
class ModelContext:
    model_id: str
    block_id: str
    base_dir: Path
    config: dict[str, str]
    crs: CRS
    bnd_xy: np.ndarray
    bnd_lonlat: np.ndarray
    candidate_node_ids: np.ndarray | None = None
    candidate_dist_km: np.ndarray | None = None
    src_xy: np.ndarray = field(default_factory=lambda: np.empty((0, 2), dtype=float))
    src_mapping: pd.DataFrame = field(default_factory=pd.DataFrame)


@dataclass
class AdcircMesh:
    node_ids: np.ndarray
    lon: np.ndarray
    lat: np.ndarray
    tree: cKDTree


@dataclass
class Fort22Data:
    times: pd.DatetimeIndex
    lon_min: float
    lat_max: float
    dlon: float
    dlat: float
    values: np.ndarray  # (time, lat, lon, [u, v, p])


def path_arg(value: str) -> Path:
    return Path(os.path.expandvars(value)).expanduser()


def return_period_tag(return_period: int) -> str:
    """Return the directory/name token used by the ADCIRC workflow."""
    if return_period <= 0:
        raise ValueError("return-period must be a positive integer number of years")
    return f"{return_period}yr"


def flow_quantile_tag(quantile: float) -> str:
    """Return a compact qXX token for automatically generated output names."""
    return f"q{int(round(float(quantile) * 100)):02d}"


def resolve_runtime_paths(args: argparse.Namespace) -> None:
    """Resolve return-period-dependent defaults while preserving explicit paths."""
    period_tag = return_period_tag(args.return_period)
    if args.adcirc_result_root is None:
        args.adcirc_result_root = args.adcirc_result_parent / period_tag
    if args.adcirc_run_root is None:
        args.adcirc_run_root = args.adcirc_run_parent / period_tag
    if args.output_root is None:
        args.output_root = (
            SCRIPT_DIR
            / f"global_sfincs_{period_tag}_{flow_quantile_tag(args.cama_quantile)}_tcr"
        )
    # The historical and future builders share the same common output root.
    # Future calls set climate_scenario and already append sspXXX before base.run().
    if (
        getattr(args, "climate_scenario", None) is None
        and args.output_root.name.casefold() != "historical"
    ):
        args.output_root = args.output_root / "historical"
    # Keep computational caches outside the runnable SFINCS output tree.  The
    # event/model directories should contain only files referenced by sfincs.inp.
    if args.q50_cache is None:
        args.q50_cache = args.cama_diagnostic_cache_dir / "builder_quantiles" / (
            f"cama_{args.cama_quantile_mode}_"
            f"q{int(round(args.cama_quantile * 100)):02d}_"
            f"{args.cama_start_year}_{args.cama_end_year}.csv"
        )


def csv_set(value: str | None) -> set[str] | None:
    if not value:
        return None
    result = {part.strip() for part in value.split(",") if part.strip()}
    return result or None


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate complete SFINCS compound-flood forcing for ERA5 return-period TC events.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--return-period",
        type=int,
        default=DEFAULT_RETURN_PERIOD,
        metavar="YEARS",
        help=(
            "TC/ADCIRC return period in years. It selects "
            "<adcirc-result-parent>/<YEARS>yr and is included in the default output name."
        ),
    )
    p.add_argument("--model-root", type=path_arg, default=DEFAULT_MODEL_ROOT)
    p.add_argument("--membership-csv", type=path_arg, default=DEFAULT_MEMBERSHIP)
    p.add_argument(
        "--adcirc-result-parent",
        type=path_arg,
        default=DEFAULT_ADCIRC_RESULT_PARENT,
        help="Parent directory containing <return-period>yr ADCIRC result folders.",
    )
    p.add_argument(
        "--adcirc-result-root",
        type=path_arg,
        default=None,
        help=(
            "Explicit ADCIRC result directory override. By default this is "
            "<adcirc-result-parent>/<return-period>yr."
        ),
    )
    p.add_argument(
        "--adcirc-run-parent",
        type=path_arg,
        default=DEFAULT_ADCIRC_RUN_PARENT,
        help=(
            "Parent directory containing <return-period>yr ADCIRC/TC rerun "
            "folders with fort.14, fort.22 and fort22_meta.txt."
        ),
    )
    p.add_argument(
        "--adcirc-run-root",
        type=path_arg,
        default=None,
        help=(
            "Explicit ADCIRC/TC run-directory override. By default this is "
            "<adcirc-run-parent>/<return-period>yr."
        ),
    )
    p.add_argument("--adcirc-mesh-root", type=path_arg, default=DEFAULT_MESH_ROOT)
    p.add_argument(
        "--boundary-map-cache-dir",
        type=path_arg,
        default=DEFAULT_BOUNDARY_MAP_CACHE_DIR,
        help=(
            "Persistent ADCIRC-node to SFINCS-boundary mapping cache. The mapping "
            "is independent of TC event and climate scenario."
        ),
    )
    p.add_argument(
        "--prepare-boundary-maps-only",
        action="store_true",
        help=(
            "Generate/validate all selected ADCIRC-to-SFINCS boundary maps, then exit "
            "before reading fort.63 or generating cases."
        ),
    )
    p.add_argument(
        "--prepare-flow-cache-only",
        action="store_true",
        help=(
            "Generate/validate the full selected CaMa-Flood flow cache, then exit. "
            "This prevents block-sharded workers from redundantly computing a partial cache."
        ),
    )
    p.add_argument(
        "--rebuild-boundary-map-cache",
        action="store_true",
        help="Force regeneration of the selected persistent boundary-map files.",
    )
    p.add_argument("--cama-root", type=path_arg, default=DEFAULT_CAMA_ROOT)
    p.add_argument("--track-nc", type=path_arg, default=DEFAULT_TRACK_NC)
    p.add_argument(
        "--t600-source",
        "--era5-monthly-root",
        dest="t600_source",
        type=path_arg,
        default=DEFAULT_ERA5_T600_ROOT,
        help=(
            "Historical mode: directory containing annual "
            "era5_t600_monthly_<year>.nc files. Future mode supplies one "
            "bias-corrected 12-month t600 NetCDF. The old "
            "--era5-monthly-root spelling remains as a compatibility alias."
        ),
    )
    p.add_argument(
        "--output-root",
        type=path_arg,
        default=None,
        help=(
            "Common historical/future output root override. Historical cases are written "
            "to its historical subdirectory. By default: "
            "global_sfincs_<return-period>yr_q<quantile>_tcr beside this script."
        ),
    )
    p.add_argument("--blocks", help="Comma-separated ADCIRC block IDs.")
    p.add_argument("--models", help="Comma-separated GTC model IDs.")
    p.add_argument("--events", help="Comma-separated exact event directory names.")
    p.add_argument("--limit-events", type=int, default=None)
    p.add_argument("--dry-run", action="store_true", help="Inventory and validate paths only.")
    p.add_argument(
        "--update-existing-inflow-only",
        action="store_true",
        help=(
            "Directly refresh sfincs.dis in already-built forcing directories from the "
            "selected diagnostic flow statistic; do not rebuild other forcing."
        ),
    )
    p.add_argument(
        "--allow-incomplete-models",
        action="store_true",
        help="Proceed with the runnable subset if a mapped base SFINCS model is incomplete.",
    )
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--continue-on-error", action="store_true")
    p.add_argument("--copy-mode", choices=("hardlink", "copy"), default="hardlink")

    p.add_argument("--buffer-days", type=float, default=2.5)
    p.add_argument("--hotstart-hours", type=float, default=6.0)
    p.add_argument("--forcing-dt-seconds", type=int, default=3600)
    p.add_argument(
        "--map-output-interval-seconds",
        type=int,
        default=None,
        help=(
            "Instantaneous sfincs_map.nc output interval; by default use the "
            "complete event duration so the NetCDF container is written only "
            "at the final time. Explicit 0 disables sfincs_map.nc on some "
            "SFINCS versions and is therefore not recommended."
        ),
    )
    p.add_argument(
        "--max-output-interval-seconds",
        type=int,
        default=None,
        help=(
            "Maximum-map aggregation interval; by default use the complete "
            "event duration so only one event-maximum map is retained."
        ),
    )
    p.add_argument(
        "--adcirc-neighbors",
        type=int,
        default=100,
        help=(
            "Number of nearest ADCIRC nodes examined for every SFINCS water-level "
            "boundary point. Candidates are restricted to no more than 1 km beyond "
            "the geometrically nearest ADCIRC node."
        ),
    )
    p.add_argument(
        "--adcirc-missing-waterlevel-fill-m",
        type=float,
        default=-1.0,
        help=(
            "Water-surface elevation assigned to dry/missing fort.63 records from the "
            "selected nearby node (default: -1 m)."
        ),
    )
    p.add_argument(
        "--min-water-valid-fraction",
        type=float,
        default=0.90,
        help=(
            "Deprecated compatibility option; node selection no longer rejects a nearby "
            "node based on a minimum valid fraction."
        ),
    )

    p.add_argument("--cama-start-year", type=int, default=1975)
    p.add_argument("--cama-end-year", type=int, default=2014)
    p.add_argument("--cama-quantile", type=float, default=0.50)
    p.add_argument(
        "--cama-quantile-mode",
        choices=(
            "tc_inlet_p50", "tc_inlet_p90", "tc_summary_p50", "tc_monthly_mean",
            "tc_weighted", "all_daily",
        ),
        default="tc_inlet_p50",
        help=(
            "tc_inlet_p50/tc_inlet_p90 independently compute the diagnostic TC-weighted "
            "percentile for each inlet cell (P50 is default); tc_summary_p50 uses the GTC-total P50; "
            "tc_monthly_mean derives a monthly mean; tc_weighted preserves the GTC-total "
            "daily P50 while allocating it to cells; all_daily uses a daily quantile."
        ),
    )
    p.add_argument(
        "--cama-diagnostic-cache-dir",
        type=path_arg,
        default=DEFAULT_CAMA_DIAGNOSTIC_CACHE_DIR,
        help="Directory containing cell_daily_<case>_*.npz diagnostic extraction caches.",
    )
    p.add_argument(
        "--cama-month-weights-csv",
        type=path_arg,
        default=DEFAULT_CAMA_MONTH_WEIGHTS_CSV,
        help="GTC-specific 12-month TC-landfall probability table.",
    )
    p.add_argument(
        "--cama-monthly-inlet-csv",
        type=path_arg,
        default=DEFAULT_CAMA_MONTHLY_INLET_CSV,
        help="Monthly inlet-cell climatology written by P1_plot_sfincs_inflow_diagnostics.py.",
    )
    p.add_argument(
        "--cama-flow-summary-csv",
        type=path_arg,
        default=DEFAULT_CAMA_FLOW_SUMMARY_CSV,
        help="GTC case summary containing tc_p50_flow_m3s.",
    )
    p.add_argument(
        "--cama-inlet-cells-csv",
        type=path_arg,
        default=DEFAULT_CAMA_INLET_CELLS_CSV,
        help="GTC inlet table whose q_reference_m3s values define inlet allocation fractions.",
    )
    p.add_argument("--cama-chunk-days", type=int, default=16)
    p.add_argument("--min-cama-valid-fraction", type=float, default=0.90)
    p.add_argument("--q50-cache", type=path_arg, default=None)
    p.add_argument("--recompute-q50", action="store_true")

    p.add_argument(
        "--rain-model",
        choices=("tcr", "rcliper", "none"),
        default="tcr",
        help="CLIMADA rain model; TCR is the physics-based option.",
    )
    p.add_argument(
        "--tcr-wind-model",
        choices=("c15", "er11"),
        default="c15",
        help=(
            "Gradient-wind workflow used internally by physics-based TCR rainfall. "
            "C15 selects the profile with surface_core/0.9 as its target peak; "
            "ER11 restores the CLIMADA default."
        ),
    )
    p.add_argument(
        "--c15-predata-dir",
        type=path_arg,
        default=DEFAULT_C15_PREDATA_DIR,
        help="Directory containing Wind_C15_data_Vmax*_Rmax*.mat lookup tables.",
    )
    p.add_argument(
        "--c15-rmax-out-of-range",
        choices=("zero", "rescale", "error"),
        default="rescale",
        help=(
            "How to handle a computed Rmax outside the C15 table range. "
            "Rescale is the validated default: it uses the nearest table profile "
            "while preserving the computed peak radius."
        ),
    )
    p.add_argument(
        "--met-forcing",
        choices=("fort22", "none"),
        default="fort22",
        help="Generate SFINCS spatial wind/pressure forcing from ADCIRC fort.22.",
    )
    p.add_argument("--forcing-resolution-m", type=float, default=5000.0)
    p.add_argument("--rain-max-eye-distance-km", type=float, default=700.0)
    p.add_argument("--rain-max-memory-gb", type=float, default=4.0)
    p.add_argument("--tcr-elevation-tif", type=path_arg, default=None)
    p.add_argument("--tcr-drag-tif", type=path_arg, default=None)
    p.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING"), default="INFO")
    return p


def parse_key_value_file(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith(("#", "!", "[", "-")) or "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.strip()] = value.strip()
    return result


def parse_sfincs_inp(path: Path) -> dict[str, str]:
    return {k.lower(): v for k, v in parse_key_value_file(path).items()}


def parse_timestamp(meta: dict[str, str], key: str) -> pd.Timestamp | None:
    value = meta.get(key)
    if value is None:
        return None
    try:
        return pd.Timestamp(value)
    except Exception as exc:
        raise ValueError(f"Invalid {key}={value!r}") from exc


def resolve_event_files(
    block_id: str,
    result_dir: Path,
    run_root: Path,
    mesh_root: Path,
) -> Event:
    event_id = result_dir.name
    run_dir = run_root / block_id / event_id
    meta_path = run_dir / "fort22_meta.txt"
    if not meta_path.is_file():
        raise FileNotFoundError(f"Missing event metadata: {meta_path}")
    meta = parse_key_value_file(meta_path)
    start = parse_timestamp(meta, "WindowStart")
    stop = parse_timestamp(meta, "WindowEnd")
    if start is None or stop is None or stop <= start:
        raise ValueError(f"Invalid WindowStart/WindowEnd in {meta_path}")
    impact = parse_timestamp(meta, "TC_aligned_impact_time")
    dt_seconds = int(round(float(meta.get("WTIMINC", "3600"))))
    fort14 = run_dir / "fort.14"
    if not fort14.is_file():
        fort14 = mesh_root / block_id / "fort.14"
    fort22 = run_dir / "fort.22"
    return Event(
        block_id=block_id,
        event_id=event_id,
        result_dir=result_dir,
        run_dir=run_dir,
        fort63=result_dir / "fort.63",
        fort14=fort14,
        fort22=fort22,
        meta_path=meta_path,
        meta=meta,
        start=start,
        stop=stop,
        impact=impact,
        dt_seconds=dt_seconds,
    )


def discover_events(args: argparse.Namespace, selected_blocks: set[str] | None) -> list[Event]:
    selected_events = csv_set(args.events)
    events: list[Event] = []
    for block_dir in sorted(args.adcirc_result_root.glob("ADC_*")):
        if not block_dir.is_dir():
            continue
        block_id = block_dir.name
        if selected_blocks and block_id not in selected_blocks:
            continue
        for result_dir in sorted(block_dir.iterdir()):
            if not result_dir.is_dir() or not (result_dir / "fort.63").is_file():
                continue
            if selected_events and result_dir.name not in selected_events:
                continue
            try:
                events.append(
                    resolve_event_files(
                        block_id,
                        result_dir,
                        args.adcirc_run_root,
                        args.adcirc_mesh_root,
                    )
                )
            except Exception as exc:
                if args.continue_on_error:
                    LOG.error("Skip invalid event %s/%s: %s", block_id, result_dir.name, exc)
                    continue
                raise
            if args.limit_events is not None and len(events) >= args.limit_events:
                return events
    return events


def load_membership(args: argparse.Namespace) -> pd.DataFrame:
    df = pd.read_csv(args.membership_csv, dtype=str)
    required = {"block_id", "model_domain_id"}
    if not required.issubset(df.columns):
        raise ValueError(f"{args.membership_csv} lacks {sorted(required - set(df.columns))}")
    df = df[["block_id", "model_domain_id"]].drop_duplicates().sort_values(
        ["block_id", "model_domain_id"]
    )
    selected_blocks = csv_set(args.blocks)
    selected_models = csv_set(args.models)
    if selected_blocks:
        df = df[df["block_id"].isin(selected_blocks)]
    if selected_models:
        df = df[df["model_domain_id"].isin(selected_models)]
    return df.reset_index(drop=True)


def read_xy(path: Path) -> np.ndarray:
    if not path.is_file() or path.stat().st_size == 0:
        return np.empty((0, 2), dtype=float)
    arr = np.loadtxt(path, dtype=float, usecols=(0, 1), ndmin=2)
    return np.asarray(arr, dtype=float)


def model_crs(config: dict[str, str]) -> CRS:
    if int(float(config.get("crsgeo", "0"))) == 1:
        return CRS.from_epsg(4326)
    epsg = config.get("epsg")
    if epsg is None:
        raise ValueError("sfincs.inp has neither crsgeo=1 nor epsg")
    return CRS.from_epsg(int(float(epsg)))


def validate_base_model(model_dir: Path) -> tuple[bool, str]:
    required = ("sfincs.inp", "sfincs.bnd", "sfincs.dep", "sfincs.msk", "sfincs.ind")
    missing = [name for name in required if not (model_dir / name).is_file()]
    if missing:
        return False, "missing " + ", ".join(missing)
    return True, "ok"


def load_src_mapping(model_dir: Path, crs: CRS) -> tuple[np.ndarray, pd.DataFrame]:
    src_path = model_dir / "sfincs.src"
    geojson = model_dir / "cama_upstream_rivers.geojson"
    if not src_path.is_file():
        return np.empty((0, 2), dtype=float), pd.DataFrame()
    src_xy = read_xy(src_path)
    if src_xy.size == 0:
        return src_xy, pd.DataFrame()
    if not geojson.is_file():
        raise FileNotFoundError(f"{model_dir.name} has sfincs.src but no {geojson.name}")
    gdf = gpd.read_file(geojson)
    needed = {"feature_type", "sfincs_src_col", "cama_row", "cama_col"}
    if not needed.issubset(gdf.columns):
        raise ValueError(f"{geojson} lacks {sorted(needed - set(gdf.columns))}")
    points = gdf[gdf["feature_type"].astype(str) == "SFINCS_Boundary_Point"].copy()
    points["src_col"] = pd.to_numeric(points["sfincs_src_col"], errors="coerce")
    points["cama_row"] = pd.to_numeric(points["cama_row"], errors="coerce")
    points["cama_col"] = pd.to_numeric(points["cama_col"], errors="coerce")
    points = points.dropna(subset=["src_col", "cama_row", "cama_col"])
    points["src_col"] = points["src_col"].round().astype(int)
    points["cama_row"] = points["cama_row"].round().astype(int)
    points["cama_col"] = points["cama_col"].round().astype(int)
    if points["src_col"].duplicated().any():
        dup = points.loc[points["src_col"].duplicated(False), "src_col"].tolist()
        raise ValueError(f"Duplicate sfincs_src_col in {geojson}: {dup}")
    points = points.sort_values("src_col")
    expected = list(range(1, len(src_xy) + 1))
    if points["src_col"].tolist() != expected:
        raise ValueError(
            f"{model_dir.name}: mapped src columns do not exactly cover 1..{len(src_xy)}"
        )
    if not (
        points["cama_row"].between(0, 719).all()
        and points["cama_col"].between(0, 1439).all()
    ):
        raise ValueError(f"{model_dir.name}: CaMa row/col must be zero-based within 720x1440")

    points["cama_cell_lon"] = -180.0 + (points["cama_col"] + 0.5) * 0.25
    points["cama_cell_lat"] = 90.0 - (points["cama_row"] + 0.5) * 0.25
    for source_name, calc_name in (("cama_lon", "cama_cell_lon"), ("cama_lat", "cama_cell_lat")):
        if source_name in points:
            source = pd.to_numeric(points[source_name], errors="coerce")
            finite = source.notna()
            if finite.any() and not np.allclose(
                source[finite], points.loc[finite, calc_name], atol=1.0e-6
            ):
                raise ValueError(f"{model_dir.name}: {source_name} disagrees with zero-based row/col")

    transformer = Transformer.from_crs(crs, 4326, always_xy=True)
    src_lon, src_lat = transformer.transform(src_xy[:, 0], src_xy[:, 1])
    keep_cols = [
        "src_col",
        "cama_row",
        "cama_col",
        "cama_cell_lon",
        "cama_cell_lat",
    ]
    for optional in ("inlet_id", "cama_boundary_id", "uparea_km2"):
        if optional in points.columns:
            keep_cols.append(optional)
    mapping = points[keep_cols].reset_index(drop=True)
    mapping.insert(1, "sfincs_x", src_xy[:, 0])
    mapping.insert(2, "sfincs_y", src_xy[:, 1])
    mapping.insert(3, "sfincs_lon", src_lon)
    mapping.insert(4, "sfincs_lat", src_lat)
    counts = mapping.groupby(["cama_row", "cama_col"])["src_col"].transform("count")
    mapping["same_cell_src_count"] = counts.astype(int)
    mapping["discharge_fraction"] = 1.0 / counts
    return src_xy, mapping


def load_model_contexts(args: argparse.Namespace, membership: pd.DataFrame) -> dict[str, list[ModelContext]]:
    contexts: dict[str, list[ModelContext]] = {}
    incomplete: list[str] = []
    for row in membership.itertuples(index=False):
        block_id = str(row.block_id)
        model_id = str(row.model_domain_id)
        model_dir = args.model_root / model_id
        ok, reason = validate_base_model(model_dir)
        if not ok:
            LOG.warning("Skip incomplete base model %s: %s", model_id, reason)
            incomplete.append(f"{model_id} ({reason})")
            continue
        config = parse_sfincs_inp(model_dir / "sfincs.inp")
        crs = model_crs(config)
        bnd_xy = read_xy(model_dir / "sfincs.bnd")
        if bnd_xy.size == 0:
            LOG.warning("Skip %s: empty sfincs.bnd", model_id)
            incomplete.append(f"{model_id} (empty sfincs.bnd)")
            continue
        to_ll = Transformer.from_crs(crs, 4326, always_xy=True)
        lon, lat = to_ll.transform(bnd_xy[:, 0], bnd_xy[:, 1])
        src_xy, src_mapping = load_src_mapping(model_dir, crs)
        contexts.setdefault(block_id, []).append(
            ModelContext(
                model_id=model_id,
                block_id=block_id,
                base_dir=model_dir,
                config=config,
                crs=crs,
                bnd_xy=bnd_xy,
                bnd_lonlat=np.column_stack([lon, lat]),
                src_xy=src_xy,
                src_mapping=src_mapping,
            )
        )
    for block_id in contexts:
        contexts[block_id].sort(key=lambda item: item.model_id)
    if incomplete and not args.dry_run and not args.allow_incomplete_models:
        raise RuntimeError(
            "Mapped base models are incomplete: "
            + "; ".join(incomplete)
            + ". Repair them or explicitly use --allow-incomplete-models."
        )
    return contexts


def lonlat_to_unit_xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lon_r = np.radians(np.asarray(lon, dtype=float))
    lat_r = np.radians(np.asarray(lat, dtype=float))
    cos_lat = np.cos(lat_r)
    return np.column_stack((cos_lat * np.cos(lon_r), cos_lat * np.sin(lon_r), np.sin(lat_r)))


def chord_to_km(chord: np.ndarray) -> np.ndarray:
    return 2.0 * EARTH_RADIUS_KM * np.arcsin(np.clip(np.asarray(chord) / 2.0, 0.0, 1.0))


def read_fort14_mesh(path: Path) -> AdcircMesh:
    LOG.info("Read ADCIRC mesh: %s", path)
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        f.readline()
        header = f.readline().split()
        if len(header) < 2:
            raise ValueError(f"Cannot parse fort.14 header: {path}")
        n_elements, n_nodes = int(header[0]), int(header[1])
        del n_elements
        ids = np.empty(n_nodes, dtype=np.int64)
        lon = np.empty(n_nodes, dtype=float)
        lat = np.empty(n_nodes, dtype=float)
        for i in range(n_nodes):
            parts = f.readline().split()
            if len(parts) < 4:
                raise EOFError(f"Unexpected EOF in fort.14 node table at row {i + 1}")
            ids[i] = int(parts[0])
            lon[i] = float(parts[1])
            lat[i] = float(parts[2])
    return AdcircMesh(ids, lon, lat, cKDTree(lonlat_to_unit_xyz(lon, lat)))


def assign_boundary_candidates(
    contexts: Sequence[ModelContext], mesh: AdcircMesh, k: int
) -> np.ndarray:
    wanted: list[np.ndarray] = []
    k_use = min(max(1, k), len(mesh.node_ids))
    for context in contexts:
        chord, idx = mesh.tree.query(
            lonlat_to_unit_xyz(context.bnd_lonlat[:, 0], context.bnd_lonlat[:, 1]),
            k=k_use,
        )
        if k_use == 1:
            chord = chord[:, None]
            idx = idx[:, None]
        context.candidate_node_ids = mesh.node_ids[idx]
        context.candidate_dist_km = chord_to_km(chord)
        wanted.append(context.candidate_node_ids.ravel())
    return np.unique(np.concatenate(wanted)) if wanted else np.empty(0, dtype=np.int64)


def _boundary_geometry_digest(context: ModelContext) -> str:
    digest = hashlib.sha256()
    digest.update(context.model_id.encode("utf-8"))
    digest.update(np.ascontiguousarray(context.bnd_lonlat, dtype="<f8").tobytes())
    return digest.hexdigest()


def _boundary_map_path(
    cache_dir: Path, block_id: str, model_id: str, neighbors: int
) -> Path:
    return cache_dir / f"k{neighbors}" / block_id / f"{model_id}.npz"


def _boundary_map_metadata(
    context: ModelContext, mesh_path: Path, neighbors: int
) -> dict[str, object]:
    stat = mesh_path.stat()
    return {
        "format_version": 2,
        "block_id": context.block_id,
        "model_id": context.model_id,
        "neighbors": int(neighbors),
        "boundary_geometry_sha256": _boundary_geometry_digest(context),
        "mesh_size_bytes": int(stat.st_size),
        "mesh_mtime_ns": int(stat.st_mtime_ns),
    }


def _load_boundary_map(
    context: ModelContext,
    mesh_path: Path,
    cache_path: Path,
    neighbors: int,
) -> bool:
    if not cache_path.is_file() or cache_path.stat().st_size <= 0:
        return False
    expected = _boundary_map_metadata(context, mesh_path, neighbors)
    try:
        with np.load(cache_path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata"].item()))
            node_ids = np.asarray(data["candidate_node_ids"], dtype=np.int64)
            distances = np.asarray(data["candidate_dist_km"], dtype=float)
        if any(metadata.get(key) != value for key, value in expected.items()):
            return False
        if (
            node_ids.shape != distances.shape
            or node_ids.ndim != 2
            or node_ids.shape[0] != len(context.bnd_xy)
            or node_ids.shape[1] != int(neighbors)
        ):
            return False
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return False
    context.candidate_node_ids = node_ids
    context.candidate_dist_km = distances
    return True


def _write_boundary_map(
    context: ModelContext,
    mesh_path: Path,
    cache_path: Path,
    neighbors: int,
) -> None:
    assert context.candidate_node_ids is not None
    assert context.candidate_dist_km is not None
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache_path.with_name(f".{cache_path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("wb") as stream:
            np.savez_compressed(
                stream,
                metadata=np.asarray(
                    json.dumps(
                        _boundary_map_metadata(context, mesh_path, neighbors),
                        sort_keys=True,
                    )
                ),
                candidate_node_ids=np.asarray(
                    context.candidate_node_ids, dtype=np.int32
                ),
                candidate_dist_km=np.asarray(
                    context.candidate_dist_km, dtype=np.float32
                ),
            )
        os.replace(temporary, cache_path)
    finally:
        if temporary.exists():
            temporary.unlink()


def load_or_build_boundary_candidates(
    args: argparse.Namespace,
    block_id: str,
    contexts: Sequence[ModelContext],
) -> np.ndarray:
    """Hydrate one block's reusable boundary mapping, rebuilding it if stale."""
    mesh_path = args.adcirc_mesh_root / block_id / "fort.14"
    if not mesh_path.is_file():
        raise FileNotFoundError(f"Missing static ADCIRC mesh for {block_id}: {mesh_path}")
    cache_paths = {
        context.model_id: _boundary_map_path(
            args.boundary_map_cache_dir,
            block_id,
            context.model_id,
            args.adcirc_neighbors,
        )
        for context in contexts
    }
    loaded = False
    if not args.rebuild_boundary_map_cache:
        loaded = all(
            _load_boundary_map(
                context,
                mesh_path,
                cache_paths[context.model_id],
                args.adcirc_neighbors,
            )
            for context in contexts
        )
    if not loaded:
        mesh = read_fort14_mesh(mesh_path)
        assign_boundary_candidates(contexts, mesh, args.adcirc_neighbors)
        for context in contexts:
            _write_boundary_map(
                context,
                mesh_path,
                cache_paths[context.model_id],
                args.adcirc_neighbors,
            )
        LOG.info(
            "Boundary map cached: %s (%s GTC model(s), k=%s)",
            block_id,
            len(contexts),
            args.adcirc_neighbors,
        )
    else:
        LOG.info(
            "Boundary map cache hit: %s (%s GTC model(s), k=%s)",
            block_id,
            len(contexts),
            args.adcirc_neighbors,
        )
    wanted = [
        context.candidate_node_ids.ravel()
        for context in contexts
        if context.candidate_node_ids is not None
    ]
    return np.unique(np.concatenate(wanted)) if wanted else np.empty(0, dtype=np.int64)


def prepare_boundary_maps(
    args: argparse.Namespace, contexts: dict[str, list[ModelContext]]
) -> None:
    for index, block_id in enumerate(sorted(contexts), start=1):
        LOG.info(
            "Prepare boundary map [%s/%s]: %s",
            index,
            len(contexts),
            block_id,
        )
        load_or_build_boundary_candidates(args, block_id, contexts[block_id])


def read_fort63_selected(path: Path, selected_node_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    selected = {int(node_id): i for i, node_id in enumerate(selected_node_ids)}
    LOG.info("Read %s selected ADCIRC nodes from %s", len(selected), path)
    with path.open("r", encoding="utf-8", errors="ignore", buffering=1024 * 1024) as f:
        f.readline()
        header = f.readline().split()
        if len(header) < 2:
            raise ValueError(f"Cannot parse fort.63 header: {path}")
        n_times, n_nodes = int(header[0]), int(header[1])
        times = np.empty(n_times, dtype=float)
        values = np.full((n_times, len(selected)), np.nan, dtype=np.float32)
        for it in range(n_times):
            record = f.readline()
            while record and not record.strip():
                record = f.readline()
            if not record:
                raise EOFError(f"Unexpected EOF before fort.63 record {it + 1}/{n_times}")
            parts = record.split()
            times[it] = float(parts[0])
            for _ in range(n_nodes):
                parts = f.readline().split()
                if len(parts) < 2:
                    raise EOFError(f"Unexpected EOF inside fort.63 record {it + 1}")
                col = selected.get(int(parts[0]))
                if col is not None:
                    value = float(parts[1])
                    if np.isfinite(value) and abs(value) < FORT63_INVALID_ABS:
                        values[it, col] = value
    if np.any(np.diff(times) <= 0):
        raise ValueError(f"fort.63 time is not strictly increasing: {path}")
    return times, values


def select_boundary_waterlevels(
    context: ModelContext,
    selected_node_ids: np.ndarray,
    native_times: np.ndarray,
    selected_values: np.ndarray,
    target_seconds: np.ndarray,
    buffer_seconds: float,
    ramp_seconds: float,
    missing_waterlevel_fill_m: float,
) -> tuple[np.ndarray, pd.DataFrame]:
    assert context.candidate_node_ids is not None
    assert context.candidate_dist_km is not None
    col_by_id = {int(node_id): i for i, node_id in enumerate(selected_node_ids)}
    n_bnd, n_candidates = context.candidate_node_ids.shape
    chosen_col = np.empty(n_bnd, dtype=int)
    chosen_rank = np.empty(n_bnd, dtype=int)
    chosen_id = np.empty(n_bnd, dtype=np.int64)
    chosen_dist = np.empty(n_bnd, dtype=float)
    valid_fraction = np.empty(n_bnd, dtype=float)
    for i in range(n_bnd):
        nearest_distance = float(context.candidate_dist_km[i, 0])
        local_limit = nearest_distance + ADCIRC_LOCAL_SEARCH_EXTRA_KM
        chosen: tuple[int, float, int, int] | None = None
        for rank in range(n_candidates):
            node_id = int(context.candidate_node_ids[i, rank])
            col = col_by_id[node_id]
            count = int(np.count_nonzero(np.isfinite(selected_values[:, col])))
            distance = float(context.candidate_dist_km[i, rank])
            if distance > local_limit + 1.0e-9:
                break
            score = (-count, distance, rank, col)
            if chosen is None or score < chosen:
                chosen = score
        assert chosen is not None
        count = -chosen[0]
        chosen_dist[i] = chosen[1]
        chosen_rank[i] = chosen[2] + 1
        chosen_col[i] = chosen[3]
        chosen_id[i] = int(context.candidate_node_ids[i, chosen[2]])
        valid_fraction[i] = count / len(native_times)
    if native_times[0] > target_seconds[-1] or native_times[-1] < target_seconds[-1]:
        raise ValueError(
            f"fort.63 coverage {native_times[0]:.0f}..{native_times[-1]:.0f} s does not reach "
            f"SFINCS stop {target_seconds[-1]:.0f} s"
        )

    out = np.zeros((len(target_seconds), n_bnd), dtype=np.float32)
    active = target_seconds > buffer_seconds
    for j in range(n_bnd):
        series = selected_values[:, chosen_col[j]].astype(float)
        # ADCIRC dry/missing records carry no usable water level. Retain the
        # nearby node and apply the configured dry-level convention rather than
        # migrating to a distant node. Valid negative water levels are unchanged.
        series = np.where(np.isfinite(series), series, missing_waterlevel_fill_m)
        use = active & (target_seconds >= native_times[0]) & (target_seconds <= native_times[-1])
        out[use, j] = np.interp(target_seconds[use], native_times, series)
    if ramp_seconds > 0:
        weights = np.clip((target_seconds - buffer_seconds) / ramp_seconds, 0.0, 1.0)
        out *= weights[:, None].astype(np.float32)
    out[target_seconds <= buffer_seconds, :] = 0.0
    match = pd.DataFrame(
        {
            "sfincs_bnd_id": np.arange(1, n_bnd + 1),
            "sfincs_x": context.bnd_xy[:, 0],
            "sfincs_y": context.bnd_xy[:, 1],
            "sfincs_lon": context.bnd_lonlat[:, 0],
            "sfincs_lat": context.bnd_lonlat[:, 1],
            "adcirc_node_id": chosen_id,
            "candidate_rank": chosen_rank,
            "match_distance_km": chosen_dist,
            "fort63_valid_fraction": valid_fraction,
        }
    )
    return out, match


def write_table_forcing(path: Path, seconds: np.ndarray, values: np.ndarray, decimals: int) -> None:
    if values.ndim != 2 or values.shape[0] != len(seconds):
        raise ValueError(f"Invalid forcing array for {path}: {values.shape}")
    value_fmt = f" {{:.{decimals}f}}"
    with path.open("w", encoding="ascii", newline="\n", buffering=1024 * 1024) as f:
        for sec, row in zip(seconds, values):
            f.write(f"{sec:.1f}" + "".join(value_fmt.format(float(v)) for v in row) + "\n")


def collect_cama_cells(contexts: dict[str, list[ModelContext]]) -> list[tuple[int, int]]:
    cells: set[tuple[int, int]] = set()
    for models in contexts.values():
        for context in models:
            if context.src_mapping.empty:
                continue
            cells.update(
                (int(row.cama_row), int(row.cama_col))
                for row in context.src_mapping.itertuples(index=False)
            )
    return sorted(cells)


def cama_year_files(root: Path, start_year: int, end_year: int) -> list[tuple[int, Path]]:
    files = [(year, root / f"outflw{year}.bin") for year in range(start_year, end_year + 1)]
    missing = [str(path) for _, path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing CaMa files: " + "; ".join(missing[:10]))
    return files


def expected_daily_dates(start_year: int, end_year: int) -> np.ndarray:
    return pd.date_range(
        f"{start_year:04d}-01-01", f"{end_year:04d}-12-31", freq="D"
    ).values.astype("datetime64[D]")


def load_cama_diagnostic_values(
    cache_dir: Path,
    case_name: str,
    cells: Sequence[tuple[int, int]],
    start_year: int,
    end_year: int,
) -> tuple[np.ndarray, np.ndarray, Path] | None:
    """Load the newest complete diagnostic cache covering all requested cells."""
    if not cache_dir.is_dir():
        return None
    candidates = sorted(
        cache_dir.glob(f"cell_daily_{case_name}_*.npz"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    wanted_dates = expected_daily_dates(start_year, end_year)
    for path in candidates:
        try:
            with np.load(path, allow_pickle=False) as npz:
                dates = np.asarray(npz["dates"]).astype("datetime64[D]")
                rows = np.asarray(npz["rows"], dtype=np.int64)
                cols = np.asarray(npz["cols"], dtype=np.int64)
                values = np.asarray(npz["values"], dtype=np.float32)
            if values.shape != (len(dates), len(rows)) or len(cols) != len(rows):
                continue
            if not np.array_equal(dates, wanted_dates):
                continue
            lookup = {
                (int(row), int(col)): index
                for index, (row, col) in enumerate(zip(rows, cols, strict=True))
            }
            if not set(cells).issubset(lookup):
                continue
            selected = values[:, [lookup[cell] for cell in cells]]
            return dates, selected, path
        except Exception as exc:
            LOG.debug("Ignore unusable CaMa diagnostic cache %s: %s", path, exc)
    return None


def load_gtc_month_weights(
    path: Path, model_ids: Sequence[str]
) -> dict[str, dict[int, float]]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing GTC TC-month weight table: {path}")
    table = pd.read_csv(path)
    required = {"gtc_id", "month", "month_weight"}
    if not required.issubset(table.columns):
        raise ValueError(f"{path} lacks {sorted(required - set(table.columns))}")
    table["gtc_id"] = table["gtc_id"].astype(str).str.strip()
    table["month"] = pd.to_numeric(table["month"], errors="coerce")
    table["month_weight"] = pd.to_numeric(table["month_weight"], errors="coerce")
    result: dict[str, dict[int, float]] = {}
    for model_id in sorted(set(model_ids)):
        subset = table[table["gtc_id"] == model_id]
        if len(subset) != 12 or set(subset["month"].dropna().astype(int)) != set(range(1, 13)):
            raise ValueError(f"{path}: {model_id} does not have exactly one row for every month")
        weights = {
            int(row.month): float(row.month_weight)
            for row in subset.itertuples(index=False)
        }
        if not np.all(np.isfinite(list(weights.values()))) or any(
            value < 0 for value in weights.values()
        ):
            raise ValueError(f"{path}: invalid month weights for {model_id}")
        if not np.isclose(sum(weights.values()), 1.0, atol=1.0e-8):
            raise ValueError(f"{path}: month weights for {model_id} do not sum to one")
        result[model_id] = weights
    return result


def tc_weighted_quantile(
    values: np.ndarray,
    months: np.ndarray,
    month_weights: dict[int, float],
    quantile: float,
) -> float:
    """Match ``weighted_quantile(tc_weighted_samples(...))`` in the diagnostic script."""
    sample_parts: list[np.ndarray] = []
    weight_parts: list[np.ndarray] = []
    values = np.asarray(values, dtype=float)
    months = np.asarray(months, dtype=int)
    for month, month_weight in sorted(month_weights.items()):
        if month_weight <= 0:
            continue
        mask = (months == month) & np.isfinite(values)
        count = int(np.count_nonzero(mask))
        if count == 0:
            continue
        sample_parts.append(values[mask])
        weight_parts.append(np.full(count, month_weight / count, dtype=float))
    if not sample_parts:
        return math.nan
    samples = np.concatenate(sample_parts)
    weights = np.concatenate(weight_parts)
    order = np.argsort(samples)
    samples = samples[order]
    cumulative = np.cumsum(weights[order])
    if cumulative[-1] <= 0:
        return math.nan
    return float(np.interp(float(quantile) * cumulative[-1], cumulative, samples))


def unique_river_contexts(
    contexts: dict[str, list[ModelContext]],
) -> dict[str, ModelContext]:
    result: dict[str, ModelContext] = {}
    for models in contexts.values():
        for context in models:
            if context.src_mapping.empty:
                continue
            previous = result.get(context.model_id)
            if previous is not None:
                left = previous.src_mapping[["cama_row", "cama_col"]].drop_duplicates()
                right = context.src_mapping[["cama_row", "cama_col"]].drop_duplicates()
                if set(map(tuple, left.to_numpy())) != set(map(tuple, right.to_numpy())):
                    raise ValueError(f"Inconsistent CaMa mappings for repeated {context.model_id}")
                continue
            result[context.model_id] = context
    return result


def diagnostic_month_weighted_cell_flows(
    monthly_csv: Path,
    weights_csv: Path,
    model_ids: Sequence[str],
    scenario: str,
    climate_models: Sequence[str],
) -> pd.DataFrame:
    """Return sum(month probability * monthly mean flow) for each inlet cell/case."""
    if not monthly_csv.is_file():
        raise FileNotFoundError(f"Missing diagnostic monthly inlet table: {monthly_csv}")
    monthly = pd.read_csv(monthly_csv)
    required = {
        "gtc_id", "scenario", "model", "cama_row", "cama_col", "month",
        "mean_flow_m3s",
    }
    if not required.issubset(monthly.columns):
        raise ValueError(f"{monthly_csv} lacks {sorted(required - set(monthly.columns))}")
    wanted_models = set(map(str, model_ids))
    wanted_climate = set(map(str, climate_models))
    monthly["gtc_id"] = monthly["gtc_id"].astype(str).str.strip()
    monthly["scenario"] = monthly["scenario"].astype(str).str.strip()
    monthly["model"] = monthly["model"].astype(str).str.strip()
    monthly = monthly[
        monthly["gtc_id"].isin(wanted_models)
        & (monthly["scenario"] == scenario)
        & monthly["model"].isin(wanted_climate)
    ].copy()
    if monthly.empty:
        raise ValueError(
            f"No diagnostic monthly flows for scenario={scenario}, models={sorted(wanted_climate)}"
        )
    for name in ("cama_row", "cama_col", "month", "mean_flow_m3s"):
        monthly[name] = pd.to_numeric(monthly[name], errors="coerce")
    if monthly[["cama_row", "cama_col", "month", "mean_flow_m3s"]].isna().any().any():
        raise ValueError(f"Non-numeric or missing flow fields in {monthly_csv}")

    weights = pd.read_csv(weights_csv, usecols=["gtc_id", "month", "month_weight"])
    weights["gtc_id"] = weights["gtc_id"].astype(str).str.strip()
    weights["month"] = pd.to_numeric(weights["month"], errors="coerce")
    weights["month_weight"] = pd.to_numeric(weights["month_weight"], errors="coerce")
    merged = monthly.merge(
        weights,
        on=["gtc_id", "month"],
        how="left",
        validate="many_to_one",
    )
    if merged["month_weight"].isna().any():
        raise ValueError("Diagnostic monthly flows lack matching GTC month weights")
    keys = ["gtc_id", "scenario", "model", "cama_row", "cama_col"]
    merged["weighted_flow_m3s"] = merged["mean_flow_m3s"] * merged["month_weight"]
    result = merged.groupby(keys, as_index=False).agg(
        q_m3s=("weighted_flow_m3s", "sum"),
        month_count=("month", "nunique"),
        month_weight_sum=("month_weight", "sum"),
    )
    if not result["month_count"].eq(12).all():
        bad = result.loc[~result["month_count"].eq(12), keys].head(10)
        raise ValueError(f"Incomplete 12-month diagnostic climatology:\n{bad.to_string(index=False)}")
    if not np.allclose(result["month_weight_sum"], 1.0, atol=1.0e-8):
        raise ValueError("Diagnostic GTC month weights do not sum to one")
    if not np.all(np.isfinite(result["q_m3s"])) or (result["q_m3s"] < 0).any():
        raise ValueError("Invalid month-weighted inlet flow")
    result["cama_row"] = result["cama_row"].astype(int)
    result["cama_col"] = result["cama_col"].astype(int)
    return result


def diagnostic_summary_p50_totals(
    summary_csv: Path,
    model_ids: Sequence[str],
    scenario: str,
    climate_models: Sequence[str],
) -> pd.DataFrame:
    """Read the final per-GTC ``tc_p50_flow_m3s`` values without recomputation."""
    if not summary_csv.is_file():
        raise FileNotFoundError(f"Missing diagnostic GTC flow summary: {summary_csv}")
    table = pd.read_csv(summary_csv)
    required = {"gtc_id", "scenario", "model", "tc_p50_flow_m3s"}
    if not required.issubset(table.columns):
        raise ValueError(f"{summary_csv} lacks {sorted(required - set(table.columns))}")
    table["gtc_id"] = table["gtc_id"].astype(str).str.strip()
    table["scenario"] = table["scenario"].astype(str).str.strip()
    table["model"] = table["model"].astype(str).str.strip()
    table["tc_p50_flow_m3s"] = pd.to_numeric(
        table["tc_p50_flow_m3s"], errors="coerce"
    )
    selected = table[
        table["gtc_id"].isin(set(map(str, model_ids)))
        & (table["scenario"] == scenario)
        & table["model"].isin(set(map(str, climate_models)))
    ][["gtc_id", "scenario", "model", "tc_p50_flow_m3s"]].copy()
    if selected.empty:
        raise ValueError(
            f"No diagnostic tc_p50 rows for scenario={scenario}, models={list(climate_models)}"
        )
    if selected.duplicated(["gtc_id", "scenario", "model"]).any():
        raise ValueError(f"Duplicate GTC/scenario/model rows in {summary_csv}")
    if selected["tc_p50_flow_m3s"].isna().any() or (
        selected["tc_p50_flow_m3s"] < 0
    ).any():
        raise ValueError(f"Invalid tc_p50_flow_m3s in {summary_csv}")
    return selected


def diagnostic_inlet_reference_fractions(
    inlet_csv: Path,
    needed: set[tuple[str, int, int]],
) -> pd.DataFrame:
    if not inlet_csv.is_file():
        raise FileNotFoundError(f"Missing diagnostic GTC inlet table: {inlet_csv}")
    table = pd.read_csv(inlet_csv)
    required = {"gtc_id", "cama_row", "cama_col", "q_reference_m3s"}
    if not required.issubset(table.columns):
        raise ValueError(f"{inlet_csv} lacks {sorted(required - set(table.columns))}")
    table["gtc_id"] = table["gtc_id"].astype(str).str.strip()
    for name in ("cama_row", "cama_col", "q_reference_m3s"):
        table[name] = pd.to_numeric(table[name], errors="coerce")
    table = table.dropna(subset=["cama_row", "cama_col", "q_reference_m3s"]).copy()
    table["cama_row"] = table["cama_row"].astype(int)
    table["cama_col"] = table["cama_col"].astype(int)
    table = table[
        table.apply(
            lambda row: (str(row.gtc_id), int(row.cama_row), int(row.cama_col))
            in needed,
            axis=1,
        )
    ].copy()
    keys = set(zip(table.gtc_id, table.cama_row, table.cama_col))
    missing = sorted(needed - keys)
    if missing:
        raise KeyError(f"Diagnostic inlet table lacks required cells: {missing[:20]}")
    if table.duplicated(["gtc_id", "cama_row", "cama_col"]).any():
        raise ValueError(f"Duplicate GTC/CaMa cell rows in {inlet_csv}")
    if (table["q_reference_m3s"] < 0).any():
        raise ValueError(f"Negative q_reference_m3s in {inlet_csv}")
    totals = table.groupby("gtc_id")["q_reference_m3s"].transform("sum")
    counts = table.groupby("gtc_id")["gtc_id"].transform("size")
    table["allocation_fraction"] = np.where(
        totals > 0, table["q_reference_m3s"] / totals, 1.0 / counts
    )
    return table[
        ["gtc_id", "cama_row", "cama_col", "q_reference_m3s", "allocation_fraction"]
    ]


def q50_cache_path(args: argparse.Namespace) -> Path:
    if args.q50_cache is not None:
        return args.q50_cache
    if args.cama_quantile_mode in {"tc_inlet_p50", "tc_inlet_p90"}:
        return args.output_root / "_cache" / (
            f"cama_independent_inlet_tc_weighted_p{int(round(args.cama_quantile * 100)):02d}_"
            f"{args.cama_start_year}_{args.cama_end_year}.csv"
        )
    if args.cama_quantile_mode == "tc_summary_p50":
        return args.output_root / "_cache" / "cama_tc_p50_from_diagnostic_summary.csv"
    if args.cama_quantile_mode == "tc_monthly_mean":
        return args.output_root / "_cache" / "cama_tc_landfall_month_weighted_mean.csv"
    if args.cama_quantile_mode == "tc_weighted":
        return args.output_root / "_cache" / (
            f"cama_tc_landfall_month_weighted_q{int(round(args.cama_quantile * 100)):02d}_"
            f"{args.cama_start_year}_{args.cama_end_year}.csv"
        )
    return args.output_root / "_cache" / (
        f"cama_daily_q{int(round(args.cama_quantile * 100)):02d}_"
        f"{args.cama_start_year}_{args.cama_end_year}.csv"
    )


def load_cama_summary_p50(
    args: argparse.Namespace,
    contexts: dict[str, list[ModelContext]] | None,
    allow_compute: bool,
) -> pd.DataFrame:
    if contexts is None:
        raise ValueError("Diagnostic summary P50 requires GTC model contexts")
    river_contexts = unique_river_contexts(contexts)
    needed = {
        (model_id, int(item.cama_row), int(item.cama_col))
        for model_id, context in river_contexts.items()
        for item in context.src_mapping[["cama_row", "cama_col"]]
        .drop_duplicates()
        .itertuples(index=False)
    }
    cache = q50_cache_path(args)
    definition = "direct_gtc_tc_p50_flow_m3s_from_diagnostic_summary"
    columns = [
        "model_id", "cama_row", "cama_col", "q_m3s", "gtc_total_q_m3s",
        "allocation_fraction", "q_reference_m3s", "definition",
        "flow_summary_csv", "inlet_cells_csv",
    ]
    if cache.is_file() and not args.recompute_q50:
        cached = pd.read_csv(cache)
        if set(columns).issubset(cached.columns):
            keys = set(
                zip(
                    cached.model_id.astype(str),
                    cached.cama_row.astype(int),
                    cached.cama_col.astype(int),
                )
            )
            valid = (
                needed.issubset(keys)
                and set(cached.definition.astype(str)) == {definition}
                and set(cached.flow_summary_csv.astype(str))
                == {str(args.cama_flow_summary_csv)}
                and set(cached.inlet_cells_csv.astype(str))
                == {str(args.cama_inlet_cells_csv)}
            )
            if valid:
                LOG.info("Use cached diagnostic-summary TC P50: %s", cache)
                return cached[
                    cached.apply(
                        lambda row: (
                            str(row.model_id), int(row.cama_row), int(row.cama_col)
                        ) in needed,
                        axis=1,
                    )
                ].copy()
    if not allow_compute or not needed:
        return pd.DataFrame(columns=columns)

    totals = diagnostic_summary_p50_totals(
        args.cama_flow_summary_csv,
        sorted(river_contexts),
        "historical",
        ["ERA5"],
    ).rename(columns={"gtc_id": "model_id", "tc_p50_flow_m3s": "gtc_total_q_m3s"})
    fractions = diagnostic_inlet_reference_fractions(
        args.cama_inlet_cells_csv, needed
    ).rename(columns={"gtc_id": "model_id"})
    result = fractions.merge(
        totals[["model_id", "gtc_total_q_m3s"]],
        on="model_id",
        how="left",
        validate="many_to_one",
    )
    if result["gtc_total_q_m3s"].isna().any():
        raise KeyError("Diagnostic flow summary lacks one or more required ERA5 GTCs")
    result["q_m3s"] = result["gtc_total_q_m3s"] * result["allocation_fraction"]
    result["cama_lon"] = -180.0 + (result["cama_col"] + 0.5) * 0.25
    result["cama_lat"] = 90.0 - (result["cama_row"] + 0.5) * 0.25
    result["definition"] = definition
    result["allocation_method"] = "q_reference_m3s_fraction"
    result["flow_summary_csv"] = str(args.cama_flow_summary_csv)
    result["inlet_cells_csv"] = str(args.cama_inlet_cells_csv)
    cache.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(cache, index=False, encoding="utf-8-sig")
    LOG.info("Wrote direct diagnostic-summary TC P50 cache: %s", cache)
    return result


def load_cama_month_weighted_mean(
    args: argparse.Namespace,
    contexts: dict[str, list[ModelContext]] | None,
    allow_compute: bool,
) -> pd.DataFrame:
    if contexts is None:
        raise ValueError("Month-weighted CaMa mean requires GTC model contexts")
    river_contexts = unique_river_contexts(contexts)
    needed = {
        (model_id, int(item.cama_row), int(item.cama_col))
        for model_id, context in river_contexts.items()
        for item in context.src_mapping[["cama_row", "cama_col"]]
        .drop_duplicates()
        .itertuples(index=False)
    }
    columns = [
        "model_id", "cama_row", "cama_col", "q_m3s", "gtc_total_q_m3s",
        "allocation_fraction", "definition", "month_weights_csv", "monthly_inlet_csv",
    ]
    cache = q50_cache_path(args)
    if cache.is_file() and not args.recompute_q50:
        cached = pd.read_csv(cache)
        if set(columns).issubset(cached.columns):
            cached_keys = set(
                zip(
                    cached.model_id.astype(str),
                    cached.cama_row.astype(int),
                    cached.cama_col.astype(int),
                )
            )
            valid_definition = set(cached.definition.astype(str)) == {
                "sum_of_gtc_landfall_month_probability_times_monthly_inlet_mean"
            }
            valid_sources = (
                set(cached.month_weights_csv.astype(str)) == {str(args.cama_month_weights_csv)}
                and set(cached.monthly_inlet_csv.astype(str))
                == {str(args.cama_monthly_inlet_csv)}
            )
            if needed.issubset(cached_keys) and valid_definition and valid_sources:
                LOG.info("Use cached TC-month-weighted mean CaMa flow: %s", cache)
                return cached[
                    cached.apply(
                        lambda row: (
                            str(row.model_id), int(row.cama_row), int(row.cama_col)
                        ) in needed,
                        axis=1,
                    )
                ].copy()
    if not allow_compute or not needed:
        return pd.DataFrame(columns=columns)

    flow = diagnostic_month_weighted_cell_flows(
        args.cama_monthly_inlet_csv,
        args.cama_month_weights_csv,
        sorted(river_contexts),
        "historical",
        ["ERA5"],
    ).rename(columns={"gtc_id": "model_id"})
    flow_keys = set(zip(flow.model_id, flow.cama_row, flow.cama_col))
    missing = sorted(needed - flow_keys)
    if missing:
        raise KeyError(f"Diagnostic monthly table lacks required ERA5 inlet cells: {missing[:20]}")
    flow = flow[
        flow.apply(
            lambda row: (str(row.model_id), int(row.cama_row), int(row.cama_col))
            in needed,
            axis=1,
        )
    ].copy()
    totals = flow.groupby("model_id")["q_m3s"].transform("sum")
    flow["gtc_total_q_m3s"] = totals
    counts = flow.groupby("model_id")["model_id"].transform("size")
    flow["allocation_fraction"] = np.where(
        totals > 0, flow["q_m3s"] / totals, 1.0 / counts
    )
    flow["cama_lon"] = -180.0 + (flow["cama_col"] + 0.5) * 0.25
    flow["cama_lat"] = 90.0 - (flow["cama_row"] + 0.5) * 0.25
    flow["definition"] = (
        "sum_of_gtc_landfall_month_probability_times_monthly_inlet_mean"
    )
    flow["month_weights_csv"] = str(args.cama_month_weights_csv)
    flow["monthly_inlet_csv"] = str(args.cama_monthly_inlet_csv)
    cache.parent.mkdir(parents=True, exist_ok=True)
    flow.to_csv(cache, index=False, encoding="utf-8-sig")
    LOG.info("Wrote TC-month-weighted mean CaMa flow cache: %s", cache)
    return flow


def load_or_compute_tc_weighted_cama_quantile(
    args: argparse.Namespace,
    cells: Sequence[tuple[int, int]],
    contexts: dict[str, list[ModelContext]] | None,
    allow_compute: bool,
) -> pd.DataFrame:
    if contexts is None:
        raise ValueError("TC-weighted CaMa q50 requires GTC model contexts")
    river_contexts = unique_river_contexts(contexts)
    model_ids = sorted(river_contexts)
    columns = [
        "model_id",
        "cama_row",
        "cama_col",
        "q_m3s",
        "gtc_total_q_m3s",
        "allocation_fraction",
        "valid_count",
        "expected_count",
    ]
    cache = q50_cache_path(args)
    independent_inlets = args.cama_quantile_mode in {"tc_inlet_p50", "tc_inlet_p90"}
    definition = (
        "independent_inlet_tc_landfall_month_probability_weighted_quantile"
        if independent_inlets
        else "gtc_total_tc_landfall_month_probability_weighted_quantile_allocated_by_cell_quantile"
    )
    needed = {
        (model_id, int(row.cama_row), int(row.cama_col))
        for model_id, context in river_contexts.items()
        for row in context.src_mapping[["cama_row", "cama_col"]]
        .drop_duplicates()
        .itertuples(index=False)
    }
    if cache.is_file() and not args.recompute_q50:
        cached = pd.read_csv(cache)
        required = set(columns) | {"definition", "month_weights_csv"}
        if required.issubset(cached.columns):
            cached_keys = set(
                zip(
                    cached.model_id.astype(str),
                    cached.cama_row.astype(int),
                    cached.cama_col.astype(int),
                )
            )
            definition_ok = set(cached.definition.astype(str)) == {definition}
            weights_ok = set(cached.month_weights_csv.astype(str)) == {
                str(args.cama_month_weights_csv)
            }
            if needed.issubset(cached_keys) and definition_ok and weights_ok:
                LOG.info("Use cached TC-month-weighted CaMa q50: %s", cache)
                return cached[
                    cached.apply(
                        lambda row: (
                            str(row.model_id), int(row.cama_row), int(row.cama_col)
                        )
                        in needed,
                        axis=1,
                    )
                ].copy()
    if not allow_compute:
        return pd.DataFrame(columns=columns)
    if not needed:
        return pd.DataFrame(columns=columns)

    loaded = load_cama_diagnostic_values(
        args.cama_diagnostic_cache_dir,
        "ERA5",
        cells,
        args.cama_start_year,
        args.cama_end_year,
    )
    if loaded is None:
        raise FileNotFoundError(
            "No complete ERA5 diagnostic NPZ cache covers the selected CaMa cells and "
            f"{args.cama_start_year}-{args.cama_end_year}: {args.cama_diagnostic_cache_dir}"
        )
    dates, values, source_cache = loaded
    LOG.info("Use ERA5 daily CaMa diagnostic values: %s", source_cache)
    weights_by_model = load_gtc_month_weights(args.cama_month_weights_csv, model_ids)
    cell_index = {cell: index for index, cell in enumerate(cells)}
    months = pd.DatetimeIndex(dates).month.to_numpy(dtype=int)
    rows: list[dict[str, object]] = []
    expected_count = len(dates)
    for model_id, context in river_contexts.items():
        weights = weights_by_model[model_id]
        unique_cells = context.src_mapping[["cama_row", "cama_col"]].drop_duplicates()
        model_cells = [
            (int(item.cama_row), int(item.cama_col))
            for item in unique_cells.itertuples(index=False)
        ]
        model_values = values[:, [cell_index[cell] for cell in model_cells]].astype(float)
        any_valid = np.isfinite(model_values).any(axis=1)
        total_series = np.nansum(model_values, axis=1)
        total_series[~any_valid] = np.nan
        total_q = tc_weighted_quantile(
            total_series, months, weights, args.cama_quantile
        )
        cell_q_values: list[float] = []
        cell_valid_counts: list[int] = []
        for cell_index_in_model, cell in enumerate(model_cells):
            series = model_values[:, cell_index_in_model]
            valid_count = int(np.count_nonzero(np.isfinite(series)))
            if (
                not independent_inlets
                and valid_count / expected_count < args.min_cama_valid_fraction
            ):
                raise ValueError(f"{model_id}/{cell}: CaMa valid fraction is too low")
            q_value = tc_weighted_quantile(
                series, months, weights, args.cama_quantile
            )
            if not np.isfinite(q_value):
                raise ValueError(f"{model_id}/{cell}: non-finite TC-weighted CaMa quantile")
            cell_q_values.append(q_value)
            cell_valid_counts.append(valid_count)
        if not np.isfinite(total_q):
            raise ValueError(f"{model_id}: non-finite total TC-weighted CaMa quantile")
        cell_q_array = np.asarray(cell_q_values, dtype=float)
        raw_sum = float(cell_q_array.sum())
        if raw_sum > 0:
            fractions = cell_q_array / raw_sum
        else:
            fractions = np.full(len(model_cells), 1.0 / len(model_cells))
        if independent_inlets:
            allocated = cell_q_array
            total_q = float(cell_q_array.sum())
            fractions = (
                cell_q_array / total_q
                if total_q > 0
                else np.full(len(model_cells), 1.0 / len(model_cells))
            )
        else:
            allocated = total_q * fractions
            if not np.isclose(float(allocated.sum()), total_q):
                raise RuntimeError(f"{model_id}: allocated CaMa discharge is not conservative")
        for cell, q_value, fraction, valid_count in zip(
            model_cells,
            allocated,
            fractions,
            cell_valid_counts,
            strict=True,
        ):
            rows.append(
                {
                    "model_id": model_id,
                    "cama_row": cell[0],
                    "cama_col": cell[1],
                    "cama_lon": -180.0 + (cell[1] + 0.5) * 0.25,
                    "cama_lat": 90.0 - (cell[0] + 0.5) * 0.25,
                    "q_m3s": q_value,
                    "gtc_total_q_m3s": total_q,
                    "allocation_fraction": fraction,
                    "valid_count": valid_count,
                    "expected_count": expected_count,
                    "quantile_nonexceedance": args.cama_quantile,
                    "start_year": args.cama_start_year,
                    "end_year": args.cama_end_year,
                    "definition": definition,
                    "month_weights_csv": str(args.cama_month_weights_csv),
                    "daily_values_cache": str(source_cache),
                }
            )
    result = pd.DataFrame(rows)
    cache.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(cache, index=False, encoding="utf-8-sig")
    LOG.info("Wrote TC-month-weighted CaMa quantile cache: %s", cache)
    return result


def load_or_compute_cama_quantile(
    args: argparse.Namespace,
    cells: Sequence[tuple[int, int]],
    allow_compute: bool = True,
    contexts: dict[str, list[ModelContext]] | None = None,
) -> pd.DataFrame:
    if args.cama_quantile_mode == "tc_summary_p50":
        return load_cama_summary_p50(args, contexts, allow_compute)
    if args.cama_quantile_mode == "tc_monthly_mean":
        return load_cama_month_weighted_mean(args, contexts, allow_compute)
    if args.cama_quantile_mode in {"tc_inlet_p50", "tc_inlet_p90", "tc_weighted"}:
        return load_or_compute_tc_weighted_cama_quantile(
            args, cells, contexts, allow_compute
        )
    columns = ["cama_row", "cama_col", "q_m3s", "valid_count", "expected_count"]
    cache = q50_cache_path(args)
    needed = set(cells)
    if cache.is_file() and not args.recompute_q50:
        cached = pd.read_csv(cache)
        if set(columns).issubset(cached.columns):
            cached_cells = set(zip(cached["cama_row"].astype(int), cached["cama_col"].astype(int)))
            if needed.issubset(cached_cells):
                LOG.info("Use cached CaMa P50: %s", cache)
                return cached[
                    cached.apply(lambda r: (int(r.cama_row), int(r.cama_col)) in needed, axis=1)
                ].copy()
    if not allow_compute:
        return pd.DataFrame(columns=columns)
    if not cells:
        return pd.DataFrame(columns=columns)

    year_files = cama_year_files(args.cama_root, args.cama_start_year, args.cama_end_year)
    nx, ny = 1440, 720
    ncell = nx * ny
    flat_indices = np.array([row * nx + col for row, col in cells], dtype=np.int64)
    expected_total = sum(366 if calendar.isleap(year) else 365 for year, _ in year_files)
    all_values = np.full((expected_total, len(cells)), np.nan, dtype=np.float32)
    cursor = 0
    LOG.info(
        "Compute daily CaMa Q%.0f for %s cells by one sequential scan of %s years",
        args.cama_quantile * 100,
        len(cells),
        len(year_files),
    )
    for year, path in year_files:
        expected_days = 366 if calendar.isleap(year) else 365
        expected_bytes = expected_days * ncell * np.dtype("<f4").itemsize
        if path.stat().st_size != expected_bytes:
            raise ValueError(
                f"Unexpected CaMa file size for {year}: {path.stat().st_size}, expected {expected_bytes}"
            )
        LOG.info("  CaMa %s", year)
        with path.open("rb") as f:
            days_read = 0
            while days_read < expected_days:
                n_days = min(args.cama_chunk_days, expected_days - days_read)
                chunk = np.fromfile(f, dtype="<f4", count=n_days * ncell)
                if chunk.size != n_days * ncell:
                    raise EOFError(f"Short read in {path} at day {days_read}")
                selected = chunk.reshape(n_days, ncell)[:, flat_indices]
                bad = (~np.isfinite(selected)) | (selected <= CAMA_FILL_LOW) | (selected >= CAMA_FILL_HIGH)
                selected = selected.astype(np.float32, copy=True)
                selected[bad] = np.nan
                selected[selected < 0] = 0.0
                all_values[cursor : cursor + n_days] = selected
                cursor += n_days
                days_read += n_days
    if cursor != expected_total:
        raise RuntimeError(f"CaMa sample count {cursor} != expected {expected_total}")
    valid_count = np.count_nonzero(np.isfinite(all_values), axis=0)
    valid_fraction = valid_count / expected_total
    if np.any(valid_fraction < args.min_cama_valid_fraction):
        bad_cells = [cells[i] for i in np.flatnonzero(valid_fraction < args.min_cama_valid_fraction)]
        raise ValueError(f"CaMa cells below valid-data threshold: {bad_cells}")
    q = np.nanquantile(all_values, args.cama_quantile, axis=0)
    result = pd.DataFrame(
        {
            "cama_row": [cell[0] for cell in cells],
            "cama_col": [cell[1] for cell in cells],
            "cama_lon": [-180.0 + (cell[1] + 0.5) * 0.25 for cell in cells],
            "cama_lat": [90.0 - (cell[0] + 0.5) * 0.25 for cell in cells],
            "q_m3s": q,
            "valid_count": valid_count,
            "expected_count": expected_total,
            "quantile_nonexceedance": args.cama_quantile,
            "start_year": args.cama_start_year,
            "end_year": args.cama_end_year,
            "definition": "all_daily_routed_outflow_nonexceedance_quantile",
        }
    )
    cache.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(cache, index=False, encoding="utf-8-sig")
    LOG.info("Wrote CaMa quantile cache: %s", cache)
    return result


def discharge_for_context(
    context: ModelContext,
    q_table: pd.DataFrame,
    target_seconds: np.ndarray,
) -> tuple[np.ndarray | None, pd.DataFrame]:
    if context.src_mapping.empty:
        return None, context.src_mapping
    context_q = q_table
    if "model_id" in q_table.columns:
        context_q = q_table[q_table["model_id"].astype(str) == context.model_id]
        if context_q.empty:
            raise KeyError(f"No TC-month-weighted CaMa discharge for {context.model_id}")
    q_lookup = {
        (int(row.cama_row), int(row.cama_col)): float(row.q_m3s)
        for row in context_q.itertuples(index=False)
    }
    mapping = context.src_mapping.copy()
    mapping_cells = [
        (int(row), int(col))
        for row, col in zip(mapping["cama_row"], mapping["cama_col"], strict=True)
    ]
    missing = sorted(set(mapping_cells) - set(q_lookup))
    if missing:
        raise KeyError(f"Missing CaMa discharge for {context.model_id}: {missing}")
    mapping["cell_q50_m3s"] = [q_lookup[cell] for cell in mapping_cells]
    mapping["assigned_q50_m3s"] = mapping["cell_q50_m3s"] * mapping["discharge_fraction"]
    mapping["cell_flow_m3s"] = mapping["cell_q50_m3s"]
    mapping["assigned_flow_m3s"] = mapping["assigned_q50_m3s"]
    q_src = mapping["assigned_q50_m3s"].to_numpy(dtype=np.float64)
    values = np.broadcast_to(q_src[None, :], (len(target_seconds), len(q_src))).copy()
    for _, group in mapping.groupby(["cama_row", "cama_col"]):
        if not np.isclose(group["assigned_q50_m3s"].sum(), group["cell_q50_m3s"].iloc[0]):
            raise RuntimeError(f"Discharge split is not conservative for {context.model_id}")
    return values, mapping


def update_existing_inflow_only(
    args: argparse.Namespace,
    contexts: dict[str, list[ModelContext]],
    q_table: pd.DataFrame,
) -> int:
    """Refresh only river forcing in existing complete SFINCS directories."""
    context_lookup = {
        (context.block_id, context.model_id): context
        for models in contexts.values()
        for context in models
        if not context.src_mapping.empty
    }
    updated = 0
    skipped = 0
    for (block_id, model_id), context in sorted(context_lookup.items()):
        block_dir = args.output_root / block_id
        if not block_dir.is_dir():
            continue
        for event_dir in sorted(path for path in block_dir.iterdir() if path.is_dir()):
            target = event_dir / model_id
            if not (target / "sfincs.inp").is_file():
                continue
            time_source = target / "sfincs.dis"
            if not time_source.is_file():
                time_source = target / "sfincs.bzs"
            if not time_source.is_file():
                LOG.warning(
                    "Skip %s: neither sfincs.dis nor sfincs.bzs supplies timestamps",
                    target,
                )
                skipped += 1
                continue
            table = np.loadtxt(time_source, ndmin=2)
            if table.ndim != 2 or table.shape[0] == 0:
                LOG.warning("Skip %s: invalid timestamp table %s", target, time_source)
                skipped += 1
                continue
            target_seconds = np.asarray(table[:, 0], dtype=float)
            discharge, _ = discharge_for_context(context, q_table, target_seconds)
            if discharge is None:
                continue
            write_table_forcing(
                target / "sfincs.dis", target_seconds, discharge, decimals=6
            )
            update_sfincs_inp(
                target / "sfincs.inp",
                {"srcfile": "sfincs.src", "disfile": "sfincs.dis"},
            )
            updated += 1
    LOG.info(
        "Updated sfincs.dis in %s existing directories; skipped=%s",
        updated,
        skipped,
    )
    return 0


MUTABLE_TOP_LEVEL = {
    "sfincs.inp",
    "sfincs.bzs",
    "sfincs.dis",
    "sfincs.wnd",
    "sfincs_precipitation.nc",
    "sfincs_wind.nc",
    "sfincs_pressure.nc",
    "forcing_complete.json",
    "forcing_failed.json",
    "adcirc_bnd_match.csv",
    "cama_src_match.csv",
}


def sfincs_file_references(inp_path: Path) -> dict[str, Path]:
    """Return relative files explicitly referenced by a SFINCS input file."""
    references: dict[str, Path] = {}
    for key, raw_value in parse_sfincs_inp(inp_path).items():
        if not key.endswith("file"):
            continue
        value = raw_value.strip().strip('"').strip("'")
        if not value or value.casefold() in {"none", "null"}:
            continue
        references[key] = Path(value)
    return references


def runtime_model_complete(target: Path) -> bool:
    """A model is complete when every file named by sfincs.inp exists."""
    inp_path = target / "sfincs.inp"
    if not inp_path.is_file() or inp_path.stat().st_size == 0:
        return False
    try:
        references = sfincs_file_references(inp_path)
    except (OSError, ValueError):
        return False
    for reference in references.values():
        path = reference if reference.is_absolute() else target / reference
        if not path.is_file() or path.stat().st_size == 0:
            return False
    return True


def copy_runtime_file(source: Path, destination: Path, copy_mode: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if copy_mode == "hardlink":
        try:
            os.link(source, destination)
            return
        except OSError:
            pass
    shutil.copy2(source, destination)


def install_shared_table_forcing(
    args: argparse.Namespace,
    context: ModelContext,
    filename: str,
    seconds: np.ndarray,
    values: np.ndarray,
    decimals: int,
    destination: Path,
) -> Path:
    """Write invariant table forcing once and link/copy it into each case."""
    digest = hashlib.sha256()
    digest.update(filename.encode("utf-8"))
    digest.update(str(decimals).encode("ascii"))
    digest.update(np.ascontiguousarray(seconds, dtype="<f8").tobytes())
    digest.update(np.ascontiguousarray(values, dtype="<f8").tobytes())
    shared = (
        args.output_root
        / "_shared_inputs"
        / context.model_id
        / filename.removeprefix("sfincs.")
        / f"{digest.hexdigest()[:24]}_{filename}"
    )
    if not shared.is_file() or shared.stat().st_size <= 0:
        shared.parent.mkdir(parents=True, exist_ok=True)
        temporary = shared.with_name(f".{shared.name}.{os.getpid()}.tmp")
        write_table_forcing(temporary, seconds, values, decimals=decimals)
        try:
            # Exclusive publication avoids two workers replacing the same cache inode.
            os.link(temporary, shared)
        except FileExistsError:
            pass
        except OSError:
            if not shared.exists():
                os.replace(temporary, shared)
        finally:
            if temporary.exists():
                temporary.unlink()
    if destination.exists():
        destination.unlink()
    copy_runtime_file(shared, destination, args.copy_mode)
    return shared


def is_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def prepare_output_model(
    context: ModelContext,
    target: Path,
    output_root: Path,
    overwrite: bool,
    copy_mode: str,
) -> bool:
    if runtime_model_complete(target) and not overwrite:
        LOG.info("Skip complete output: %s", target)
        return False
    if target.exists():
        if not overwrite:
            raise FileExistsError(f"Incomplete output exists; use --overwrite: {target}")
        if not is_within(target, output_root) or target.resolve() == output_root.resolve():
            raise ValueError(f"Refuse to remove unsafe output target: {target}")
        shutil.rmtree(target)
    target.mkdir(parents=True)
    references = sfincs_file_references(context.base_dir / "sfincs.inp")
    copied: set[Path] = set()
    for reference in references.values():
        if reference.is_absolute():
            continue
        source = context.base_dir / reference
        if source.name in MUTABLE_TOP_LEVEL or not source.is_file():
            continue
        destination = target / reference
        if destination in copied:
            continue
        copy_runtime_file(source, destination, copy_mode)
        copied.add(destination)
    shutil.copy2(context.base_dir / "sfincs.inp", target / "sfincs.inp")
    return True


def update_sfincs_inp(path: Path, updates: dict[str, str | int | float | None]) -> None:
    original = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    normalized = {key.lower(): value for key, value in updates.items()}
    seen: set[str] = set()
    output: list[str] = []
    key_re = re.compile(r"^\s*([A-Za-z0-9_]+)\s*=")
    for line in original:
        match = key_re.match(line)
        if not match:
            output.append(line)
            continue
        key = match.group(1).lower()
        if key not in normalized:
            output.append(line)
            continue
        seen.add(key)
        value = normalized[key]
        if value is not None:
            output.append(f"{key:<20s} = {value}")
    for key, value in normalized.items():
        if key not in seen and value is not None:
            output.append(f"{key:<20s} = {value}")
    path.write_text("\n".join(output) + "\n", encoding="utf-8")


def model_bounds(config: dict[str, str]) -> tuple[float, float, float, float]:
    x0 = float(config["x0"])
    y0 = float(config["y0"])
    width = int(float(config["mmax"])) * float(config["dx"])
    height = int(float(config["nmax"])) * float(config["dy"])
    angle = math.radians(float(config.get("rotation", "0")))
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    corners = []
    for local_x, local_y in ((0, 0), (width, 0), (width, height), (0, height)):
        corners.append(
            (
                x0 + local_x * cos_a - local_y * sin_a,
                y0 + local_x * sin_a + local_y * cos_a,
            )
        )
    xy = np.asarray(corners)
    return float(xy[:, 0].min()), float(xy[:, 1].min()), float(xy[:, 0].max()), float(xy[:, 1].max())


def forcing_grid(context: ModelContext, resolution: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if resolution <= 0:
        raise ValueError("forcing resolution must be positive")
    xmin, ymin, xmax, ymax = model_bounds(context.config)
    xmin = math.floor(xmin / resolution) * resolution
    ymin = math.floor(ymin / resolution) * resolution
    xmax = math.ceil(xmax / resolution) * resolution
    ymax = math.ceil(ymax / resolution) * resolution
    x = np.arange(xmin, xmax + 0.5 * resolution, resolution, dtype=float)
    y = np.arange(ymin, ymax + 0.5 * resolution, resolution, dtype=float)
    if len(x) < 2:
        x = np.array([xmin, xmin + resolution])
    if len(y) < 2:
        y = np.array([ymin, ymin + resolution])
    xx, yy = np.meshgrid(x, y)
    to_ll = Transformer.from_crs(context.crs, 4326, always_xy=True)
    lon, lat = to_ll.transform(xx, yy)
    return x, y, np.asarray(lon), np.asarray(lat)


def read_fort22(event: Event) -> Fort22Data:
    meta = event.meta
    required = ("NWLON", "NWLAT", "WLONMIN", "WLATMAX", "WLONINC", "WLATINC")
    missing = [key for key in required if key not in meta]
    if missing:
        raise ValueError(f"{event.meta_path} lacks {missing}")
    nlon = int(float(meta["NWLON"]))
    nlat = int(float(meta["NWLAT"]))
    lon_min = float(meta["WLONMIN"])
    lat_max = float(meta["WLATMAX"])
    dlon = abs(float(meta["WLONINC"]))
    dlat = abs(float(meta["WLATINC"]))
    n_steps = int(round(event.duration_seconds / event.dt_seconds)) + 1
    expected = n_steps * nlat * nlon * 3
    LOG.info("Read ADCIRC fort.22: %s", event.fort22)
    raw = np.fromfile(event.fort22, dtype=np.float32, sep=" ")
    if raw.size != expected:
        raise ValueError(f"fort.22 has {raw.size} values; expected {expected}")
    values = raw.reshape(n_steps, nlat, nlon, 3)
    times = pd.date_range(event.start, periods=n_steps, freq=pd.Timedelta(seconds=event.dt_seconds))
    return Fort22Data(times, lon_min, lat_max, dlon, dlat, values)


def sample_regular_latlon(
    data: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    lon_min: float,
    lat_max: float,
    dlon: float,
    dlat: float,
    point_chunk: int = 50000,
    warn_outside: bool = True,
) -> np.ndarray:
    if data.ndim != 3:
        raise ValueError(f"Expected data(time,lat,lon), got {data.shape}")
    n_time, n_lat, n_lon = data.shape
    lon_flat = np.asarray(lon, dtype=float).ravel()
    lat_flat = np.asarray(lat, dtype=float).ravel()
    grid_mid = lon_min + 0.5 * (n_lon - 1) * dlon
    lon_adjusted = lon_flat + 360.0 * np.round((grid_mid - lon_flat) / 360.0)
    fx = (lon_adjusted - lon_min) / dlon
    fy = (lat_max - lat_flat) / dlat
    tolerance = 1.0e-6
    outside = (
        (fx < -tolerance)
        | (fx > n_lon - 1 + tolerance)
        | (fy < -tolerance)
        | (fy > n_lat - 1 + tolerance)
    )
    if warn_outside and np.any(outside):
        excess_cells = np.maximum.reduce(
            (
                np.maximum(-fx, 0.0),
                np.maximum(fx - (n_lon - 1), 0.0),
                np.maximum(-fy, 0.0),
                np.maximum(fy - (n_lat - 1), 0.0),
            )
        )
        LOG.warning(
            "SFINCS forcing grid has %d/%d points outside fort.22; "
            "using the nearest fort.22 edge value (maximum %.3f source-grid cells outside)",
            int(np.count_nonzero(outside)),
            int(outside.size),
            float(np.max(excess_cells[outside])),
        )
    fx = np.clip(fx, 0, n_lon - 1)
    fy = np.clip(fy, 0, n_lat - 1)
    x0 = np.floor(fx).astype(int)
    y0 = np.floor(fy).astype(int)
    x1 = np.minimum(x0 + 1, n_lon - 1)
    y1 = np.minimum(y0 + 1, n_lat - 1)
    wx = (fx - x0).astype(np.float32)
    wy = (fy - y0).astype(np.float32)
    out = np.empty((n_time, len(fx)), dtype=np.float32)
    for start in range(0, len(fx), point_chunk):
        stop = min(start + point_chunk, len(fx))
        sl = slice(start, stop)
        v00 = data[:, y0[sl], x0[sl]]
        v01 = data[:, y0[sl], x1[sl]]
        v10 = data[:, y1[sl], x0[sl]]
        v11 = data[:, y1[sl], x1[sl]]
        out[:, sl] = (
            v00 * (1 - wx[sl]) * (1 - wy[sl])
            + v01 * wx[sl] * (1 - wy[sl])
            + v10 * (1 - wx[sl]) * wy[sl]
            + v11 * wx[sl] * wy[sl]
        )
    return out.reshape((n_time,) + lon.shape)


def write_grid_netcdf(
    path: Path,
    times: pd.DatetimeIndex,
    x: np.ndarray,
    y: np.ndarray,
    variables: dict[str, np.ndarray],
    units: dict[str, str],
    title: str,
) -> None:
    data_vars = {}
    for name, values in variables.items():
        if values.shape != (len(times), len(y), len(x)):
            raise ValueError(f"{name} shape {values.shape} does not match time/y/x")
        data_vars[name] = (("time", "y", "x"), np.asarray(values, dtype=np.float32))
    ds = xr.Dataset(data_vars, coords={"time": times, "y": y, "x": x})
    for name, unit in units.items():
        ds[name].attrs["units"] = unit
    ds.attrs.update({"title": title, "Conventions": "CF-1.8"})
    encoding: dict[str, dict] = {
        name: {"dtype": "float32", "zlib": True, "complevel": 4, "shuffle": True}
        for name in variables
    }
    encoding["time"] = {
        "dtype": "float64",
        "units": f"minutes since {times[0].strftime('%Y-%m-%d %H:%M:%S')}",
        "calendar": "proleptic_gregorian",
    }
    ds.to_netcdf(path, engine="netcdf4", encoding=encoding)


def write_fort22_meteorology(
    target: Path,
    fort22: Fort22Data,
    x: np.ndarray,
    y: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
) -> None:
    u = sample_regular_latlon(
        fort22.values[..., 0], lon, lat, fort22.lon_min, fort22.lat_max, fort22.dlon, fort22.dlat
    )
    v = sample_regular_latlon(
        fort22.values[..., 1],
        lon,
        lat,
        fort22.lon_min,
        fort22.lat_max,
        fort22.dlon,
        fort22.dlat,
        warn_outside=False,
    )
    pressure = sample_regular_latlon(
        fort22.values[..., 2],
        lon,
        lat,
        fort22.lon_min,
        fort22.lat_max,
        fort22.dlon,
        fort22.dlat,
        warn_outside=False,
    )
    write_grid_netcdf(
        target / "sfincs_wind.nc",
        fort22.times,
        x,
        y,
        {"eastward_wind": u, "northward_wind": v},
        {"eastward_wind": "m s-1", "northward_wind": "m s-1"},
        "ADCIRC fort.22 C15 wind forcing for SFINCS",
    )
    write_grid_netcdf(
        target / "sfincs_pressure.nc",
        fort22.times,
        x,
        y,
        {"barometric_pressure": pressure},
        {"barometric_pressure": "Pa"},
        "ADCIRC fort.22 pressure forcing for SFINCS",
    )


def read_track_points(meta_path: Path) -> pd.DataFrame:
    lines = meta_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    header_index = next(
        (i for i, line in enumerate(lines) if line.strip().startswith("index,original_time_index,time,")),
        None,
    )
    if header_index is None:
        raise ValueError(f"No [ForcingTrackPoints] CSV table in {meta_path}")
    frame = pd.read_csv(io.StringIO("\n".join(lines[header_index:])))
    required = {"original_time_index", "time", "lon180", "lat", "vmax_ms", "pressure_hPa"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Track table lacks {sorted(required - set(frame.columns))}")
    frame["time"] = pd.to_datetime(frame["time"])
    for name in required - {"time"}:
        frame[name] = pd.to_numeric(frame[name], errors="coerce")
    return frame.dropna(subset=list(required)).sort_values("time").reset_index(drop=True)


def matlab_track_step_motion(
    lon: np.ndarray, lat: np.ndarray, times: pd.DatetimeIndex
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Match the forward storm motion used by the ADCIRC MATLAB builder.

    Returns eastward velocity, northward velocity and speed in m/s. MATLAB
    evaluates the current-to-next motion at each forcing time and uses zero
    motion at the final active track point.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    if lon.shape != lat.shape or lon.ndim != 1 or len(lon) != len(times):
        raise ValueError("Track longitude, latitude and time arrays must be one-dimensional and equal")
    dlon = (np.diff(lon) + 180.0) % 360.0 - 180.0
    dx = dlon * np.cos(np.deg2rad(0.5 * (lat[:-1] + lat[1:]))) * 111_320.0
    dy = np.diff(lat) * 110_540.0
    dt = np.diff(times.values).astype("timedelta64[s]").astype(float)
    east = np.divide(dx, dt, out=np.zeros_like(dx), where=dt > 0)
    north = np.divide(dy, dt, out=np.zeros_like(dy), where=dt > 0)
    east = np.r_[east, 0.0]
    north = np.r_[north, 0.0]
    return east, north, np.hypot(east, north)


def matlab_track_step_speed(
    lon: np.ndarray, lat: np.ndarray, times: pd.DatetimeIndex
) -> np.ndarray:
    """Return the ADCIRC MATLAB forward storm-translation speed in m/s."""
    return matlab_track_step_motion(lon, lat, times)[2]


def c15_parameters_from_track(
    surface_vmax_ms: np.ndarray,
    translation_speed_ms: np.ndarray,
    central_pressure_hpa: np.ndarray,
    latitude_deg: np.ndarray,
    cv: float,
) -> dict[str, np.ndarray]:
    """Build C15 parameters while keeping surface and gradient winds separate.

    ``surface_vmax_ms`` is the track's near-surface total ``vmax_trks`` value.
    As in the ADCIRC fort.22 builder, the translation-free *surface* rotational
    wind controls weak-storm pressure scaling and the Cv/Rmax relationship,
    while Holland B uses total surface Vmax.  Only the TCR dynamic C15 wind is
    converted to gradient level using ``surface_rotational / 0.9``.  The custom
    spatial field adds the translation vector once at a later stage.
    """
    surface_vmax_ms = np.asarray(surface_vmax_ms, dtype=float)
    translation_speed_ms = np.asarray(translation_speed_ms, dtype=float)
    central_pressure_hpa = np.asarray(central_pressure_hpa, dtype=float)
    latitude_deg = np.asarray(latitude_deg, dtype=float)
    if not (
        surface_vmax_ms.shape
        == translation_speed_ms.shape
        == central_pressure_hpa.shape
        == latitude_deg.shape
    ):
        raise ValueError("C15 track parameter arrays must have identical shapes")
    if not np.isfinite(cv) or cv <= 0:
        raise ValueError(f"Invalid Cv_used={cv!r}; a positive finite value is required")

    surface_core_raw = surface_vmax_ms - np.maximum(translation_speed_ms, 0.0)
    surface_core = np.maximum(surface_core_raw, 0.0)
    gradient_core = surface_core / TCR_GRADIENT_TO_SURFACE_WIND_FACTOR
    weak_scale = np.clip(
        surface_core / C15_MIN_LOOKUP_WIND_MS, 0.0, 1.0
    )
    surface_vm_lookup = np.clip(
        np.maximum(surface_core, C15_MIN_LOOKUP_WIND_MS),
        C15_MIN_LOOKUP_WIND_MS,
        C15_MAX_LOOKUP_WIND_MS,
    )
    pressure_model_hpa = C15_ENVIRONMENTAL_PRESSURE_HPA - weak_scale * (
        C15_ENVIRONMENTAL_PRESSURE_HPA - central_pressure_hpa
    )
    rmax_m = np.maximum(
        cv
        * 51.6
        * np.exp(-0.0223 * surface_vm_lookup + 0.0281 * np.abs(latitude_deg))
        * 1000.0,
        C15_MIN_RMAX_M,
    )

    pressure_deficit_hpa = np.maximum(
        C15_ENVIRONMENTAL_PRESSURE_HPA - pressure_model_hpa, 1.0
    )
    holland_b_raw = (
        surface_vmax_ms**2
        * C15_AIR_DENSITY_KG_M3
        * math.e
        / pressure_deficit_hpa
        / 100.0
    )
    holland_b = np.clip(
        holland_b_raw, C15_HOLLAND_B_MIN, C15_HOLLAND_B_MAX
    )
    return {
        # Backward-compatible name used by the custom TCR C15 pathway.
        "vm_raw": gradient_core,
        "surface_core_raw": surface_core_raw,
        "surface_core": surface_core,
        "gradient_core": gradient_core,
        "weak_scale": weak_scale,
        # vm_lookup is the surface-wind lookup used by the Rmax relationship.
        "vm_lookup": surface_vm_lookup,
        "surface_vm_lookup": surface_vm_lookup,
        "pressure_model_hpa": pressure_model_hpa,
        "rmax_m": rmax_m,
        "holland_b_raw": holland_b_raw,
        "holland_b": holland_b,
    }


def basin_from_event(event: Event, raw_basin: str | None) -> str:
    if raw_basin and raw_basin.lower() not in {"nan", "none", ""}:
        return {"AU": "SP"}.get(raw_basin.upper(), raw_basin.upper())
    prefix = event.block_id.split("_")[1]
    return {"WNP": "WP", "NIO": "NI", "SIO": "SI", "AUSSP": "SP", "NA": "NA"}.get(
        prefix, "NA"
    )


def saffir_category(max_wind_ms: float) -> int:
    knots = max_wind_ms / 0.514444
    if knots < 34:
        return -1
    if knots < 64:
        return 0
    if knots < 83:
        return 1
    if knots < 96:
        return 2
    if knots < 113:
        return 3
    if knots < 137:
        return 4
    return 5


def _periodic_longitude_field(
    field: xr.DataArray, lon_name: str
) -> xr.DataArray:
    """Add cyclic end cells so linear sampling works across 0/360 degrees."""
    longitude = np.asarray(field[lon_name].values, dtype=float)
    if (
        longitude.ndim != 1
        or longitude.size < 2
        or np.any(np.diff(longitude) <= 0)
    ):
        raise ValueError("Longitude coordinate must be a strictly increasing vector")
    left = field.isel({lon_name: [-1]}).assign_coords(
        {lon_name: [longitude[-1] - 360.0]}
    )
    right = field.isel({lon_name: [0]}).assign_coords(
        {lon_name: [longitude[0] + 360.0]}
    )
    return xr.concat([left, field, right], dim=lon_name)


def sample_monthly_t600(
    source: Path,
    times: pd.DatetimeIndex,
    lon: np.ndarray,
    lat: np.ndarray,
) -> np.ndarray:
    """Sample ERA5 monthly 600-hPa temperature along a historical TC track."""
    if not source.is_dir():
        raise NotADirectoryError(
            f"Historical t600 source must be an annual-file directory: {source}"
        )
    result = np.full(len(times), np.nan, dtype=float)
    groups = pd.DataFrame(
        {"year": times.year, "month": times.month}
    ).groupby(["year", "month"]).groups
    for (year, month), positions in groups.items():
        positions = np.asarray(list(positions), dtype=int)
        path = source / f"era5_t600_monthly_{year}.nc"
        if not path.is_file():
            raise FileNotFoundError(
                f"TCR t600 input is missing for {year}: {path}"
            )
        with xr.open_dataset(path) as dataset:
            name = next(
                (
                    candidate
                    for candidate in ("t600", "t", "temperature")
                    if candidate in dataset
                ),
                None,
            )
            if name is None:
                raise KeyError(f"{path}: missing t600/temperature variable")
            temperature = dataset[name]
            level_name = next(
                (
                    candidate
                    for candidate in ("level", "pressure_level", "plev")
                    if candidate in temperature.dims
                ),
                None,
            )
            if level_name is not None:
                levels = np.asarray(temperature[level_name].values, dtype=float)
                target = 600.0 if np.nanmax(levels) < 2_000.0 else 60_000.0
                temperature = temperature.sel(
                    {level_name: target}, method="nearest"
                )
            lat_name = "latitude" if "latitude" in temperature.coords else "lat"
            lon_name = "longitude" if "longitude" in temperature.coords else "lon"
            if "time" not in temperature.dims:
                raise ValueError(f"{path}: t600 has no time dimension")
            if "time" in temperature.coords:
                month_matches = np.flatnonzero(
                    np.asarray(
                        temperature["time"].dt.month.values, dtype=int
                    )
                    == int(month)
                )
                if month_matches.size != 1:
                    raise ValueError(
                        f"{path}: expected one record for month {month}, "
                        f"found {month_matches.size}"
                    )
                temperature = temperature.isel(time=int(month_matches[0]))
            else:
                if temperature.sizes["time"] != 12:
                    raise ValueError(
                        f"{path}: undecodable time dimension has "
                        f"{temperature.sizes['time']} records, expected 12"
                    )
                temperature = temperature.isel(time=int(month) - 1)
            units = (
                str(temperature.attrs.get("units", "K"))
                .strip()
                .lower()
                .replace(" ", "")
            )
            if units in {"c", "degc", "degree_celsius", "degreescelsius"}:
                temperature = temperature + 273.15
            elif units not in {"", "k", "kelvin"}:
                raise ValueError(f"{path}: unsupported t600 units {units!r}")
            if np.any(np.diff(temperature[lat_name].values) < 0):
                temperature = temperature.sortby(lat_name)
            longitude = np.mod(
                np.asarray(temperature[lon_name].values, dtype=float), 360.0
            )
            temperature = temperature.assign_coords(
                {lon_name: (lon_name, longitude)}
            ).sortby(lon_name)
            temperature = _periodic_longitude_field(temperature, lon_name)
            source_lon = np.asarray(temperature[lon_name].values, dtype=float)
            center = 0.5 * (source_lon[1] + source_lon[-2])
            query_lon = (
                np.asarray(lon[positions], dtype=float) - center + 180.0
            ) % 360.0 + center - 180.0
            sampled = temperature.interp(
                {
                    lon_name: xr.DataArray(query_lon, dims="points"),
                    lat_name: xr.DataArray(
                        np.asarray(lat[positions], dtype=float), dims="points"
                    ),
                },
                method="linear",
            ).values
            result[positions] = np.asarray(sampled, dtype=float).reshape(-1)
    bad = ~np.isfinite(result) | (result < 150.0) | (result > 350.0)
    if bad.any():
        raise ValueError(
            f"Invalid ERA5 t600 at {bad.sum()} track positions; "
            f"valid range is 150..350 K"
        )
    return result


@lru_cache(maxsize=8)
def c15_available_keys(root_text: str) -> frozenset[tuple[int, int]]:
    """Return the available (Vmax m/s, Rmax km) C15 lookup-table keys."""
    root = Path(root_text)
    keys: set[tuple[int, int]] = set()
    if root.is_dir():
        for path in root.glob("Wind_C15_data_Vmax*_Rmax*.mat"):
            match = C15_PROFILE_RE.fullmatch(path.name)
            if match:
                keys.add((int(match.group(1)), int(match.group(2))))
    if not keys:
        raise FileNotFoundError(f"No C15 wind-profile lookup tables found under {root}")
    return frozenset(keys)


@lru_cache(maxsize=4096)
def load_c15_profile(
    root_text: str, vmax_key: int, rmax_key: int
) -> tuple[np.ndarray, np.ndarray]:
    """Load one C15 radial wind profile; radii are in metres and winds in m/s."""
    from scipy.io import loadmat

    path = (
        Path(root_text)
        / f"Wind_C15_data_Vmax{vmax_key}_Rmax{rmax_key}.mat"
    )
    raw = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "Wind_C15_data" not in raw:
        raise ValueError(f"Missing Wind_C15_data variable in {path}")
    profile = raw["Wind_C15_data"]
    rr = np.asarray(profile.rr, dtype=float).reshape(-1)
    c15_profile_wind_ms = np.asarray(profile.vg, dtype=float).reshape(-1)
    valid = np.isfinite(rr) & np.isfinite(c15_profile_wind_ms)
    rr, c15_profile_wind_ms = rr[valid], c15_profile_wind_ms[valid]
    order = np.argsort(rr)
    rr = rr[order]
    c15_profile_wind_ms = np.maximum(c15_profile_wind_ms[order], 0.0)
    if rr.size < 2 or np.any(np.diff(rr) <= 0):
        raise ValueError(f"Invalid C15 radial profile in {path}")
    return rr, c15_profile_wind_ms


def matlab_round_positive(value: float) -> int:
    """Match MATLAB round for the non-negative C15 lookup indices."""
    return int(math.floor(value + 0.5))


@contextmanager
def install_c15_tcr_wind_model(tcr_module: object, args: argparse.Namespace):
    """Temporarily make TCR use the C15 gradient-wind workflow.

    The translation-free surface core is divided by 0.9 first. That gradient
    peak selects the C15 profile, so the returned profile is already the
    gradient-wind profile required by TCR and must not be divided by 0.9 a
    second time. The patched field otherwise uses the same forward translation
    vector, radius-dependent inflow angle, asymmetric ``vmoc`` term and
    500--700 km outer blend as ``calc_c15_uvp_field_grid``. ADCIRC's separate
    0.893 conversion from 1-min to 10-min surface wind is not used by TCR.
    """
    # CLIMADA-Petals uses this module-level constant inside the t600-to-q950
    # thermodynamic diagnosis.  Override it only for this calculation context
    # so q950 and the custom C15 dynamic wind use exactly the same factor,
    # without modifying the installed package on disk.
    had_petals_surface_factor = hasattr(
        tcr_module, "GRADIENT_LEVEL_TO_SURFACE_WINDS"
    )
    original_petals_surface_factor = getattr(
        tcr_module, "GRADIENT_LEVEL_TO_SURFACE_WINDS", None
    )
    tcr_module.GRADIENT_LEVEL_TO_SURFACE_WINDS = (
        TCR_GRADIENT_TO_SURFACE_WIND_FACTOR
    )
    root_text = str(args.c15_predata_dir.resolve())
    keys = c15_available_keys(root_text)
    speeds = sorted({key[0] for key in keys})
    radii_by_speed = {
        speed: np.asarray(sorted(key[1] for key in keys if key[0] == speed), dtype=int)
        for speed in speeds
    }
    original_compute = tcr_module.compute_angular_windspeeds
    original_horizontal_winds = tcr_module._horizontal_winds
    original_w_topo = tcr_module._w_topo
    original_w_frict_stretch = tcr_module._w_frict_stretch
    original_tctrack_to_si = tcr_module.tctrack_to_si
    model_map = tcr_module.MODEL_VANG
    had_c15 = "C15" in model_map
    old_c15 = model_map.get("C15")

    def resolve_profile(
        gradient_core_speed: float,
        actual_rmax_m: float,
    ) -> tuple[np.ndarray, np.ndarray, float] | None:
        requested_speed = matlab_round_positive(
            float(
                np.clip(
                    max(gradient_core_speed, C15_MIN_LOOKUP_WIND_MS),
                    C15_MIN_LOOKUP_WIND_MS,
                    C15_MAX_LOOKUP_WIND_MS,
                )
            )
        )
        requested_rmax = matlab_round_positive(actual_rmax_m / 1000.0)
        requested_key = (requested_speed, requested_rmax)
        if requested_key in keys:
            rr, gradient_profile_ms = load_c15_profile(root_text, *requested_key)
            return rr, gradient_profile_ms, 1.0

        if args.c15_rmax_out_of_range == "error":
            raise ValueError(
                f"C15 table missing Vmax={requested_speed} m/s, "
                f"Rmax={requested_rmax} km"
            )
        if args.c15_rmax_out_of_range == "zero":
            warning_key = (requested_speed, requested_rmax, 0)
            if warning_key not in C15_RESCALE_WARNED:
                LOG.warning(
                    "C15 has no Vmax=%s/Rmax=%s-km table; using zero wind "
                    "to match the ADCIRC MATLAB builder",
                    requested_speed,
                    requested_rmax,
                )
                C15_RESCALE_WARNED.add(warning_key)
            return None

        lookup_speed = min(speeds, key=lambda value: abs(value - requested_speed))
        available_radii = radii_by_speed[lookup_speed]
        lookup_rmax = int(
            available_radii[np.argmin(np.abs(available_radii - requested_rmax))]
        )
        warning_key = (requested_speed, requested_rmax, lookup_rmax)
        if warning_key not in C15_RESCALE_WARNED:
            LOG.warning(
                "C15 has no Vmax=%s/Rmax=%s-km table; using Vmax=%s/Rmax=%s-km "
                "profile rescaled to the computed radius",
                requested_speed,
                requested_rmax,
                lookup_speed,
                lookup_rmax,
            )
            C15_RESCALE_WARNED.add(warning_key)
        rr, gradient_profile_ms = load_c15_profile(
            root_text, lookup_speed, lookup_rmax
        )
        return (
            rr,
            gradient_profile_ms,
            lookup_rmax * 1000.0 / actual_rmax_m,
        )

    def c15_gradient_profile(
        si_track: xr.Dataset,
        radius_m: np.ndarray,
        mask_centr_close: np.ndarray,
    ) -> np.ndarray:
        """Return the C15 profile selected by the target gradient-wind peak.

        The target peak is ``surface_core / 0.9``. Therefore the selected
        profile already represents TCR gradient wind; no second division by
        0.9 is applied to ``C15.vg``.
        """
        result = np.zeros_like(radius_m, dtype=float)
        if "vmax_gradient" not in si_track:
            raise ValueError("C15 TCR SI track lacks gradient-level Vmax")
        gradient_core_speeds = np.asarray(
            si_track["vmax_gradient"].values, dtype=float
        )
        radii_m = np.asarray(si_track["rad"].values, dtype=float)

        for time_i, (gradient_core_speed, actual_rmax_m) in enumerate(
            zip(gradient_core_speeds, radii_m, strict=True)
        ):
            close = np.asarray(mask_centr_close[time_i], dtype=bool)
            if (
                not close.any()
                or not np.isfinite(gradient_core_speed)
                or not np.isfinite(actual_rmax_m)
                or actual_rmax_m <= 0
            ):
                continue
            weak_scale = float(
                np.clip(
                    max(gradient_core_speed, 0.0) / C15_MIN_LOOKUP_WIND_MS,
                    0.0,
                    1.0,
                )
            )
            if weak_scale <= 0:
                continue
            profile = resolve_profile(
                float(gradient_core_speed), float(actual_rmax_m)
            )
            if profile is None:
                continue
            rr, gradient_profile_ms, radius_scale = profile
            query_radius = radius_m[time_i, close] * radius_scale
            result[time_i, close] = weak_scale * np.interp(
                query_radius,
                rr,
                gradient_profile_ms,
                left=0.0,
                right=0.0,
            )
        return result

    def matlab_total_wind_vectors(
        si_track: xr.Dataset,
        radius_m: np.ndarray,
        radial_direction: np.ndarray,
        mask_centr_close: np.ndarray,
    ) -> np.ndarray:
        """Return [northward, eastward] C15 gradient-level winds for TCR."""
        radius = np.maximum(np.asarray(radius_m, dtype=float), 1.0)
        gradient_wind_ms = c15_gradient_profile(
            si_track, radius, mask_centr_close
        )
        rmax = np.asarray(si_track["rad"].values, dtype=float)[:, None]
        gradient_core_speed = np.asarray(
            si_track["vmax_gradient"].values, dtype=float
        )[:, None]
        weak_scale = np.clip(
            np.maximum(gradient_core_speed, 0.0) / C15_MIN_LOOKUP_WIND_MS,
            0.0,
            1.0,
        )

        beta = np.full_like(radius, C15_BETA_OUTER_DEG, dtype=float)
        inner = radius < rmax
        middle = (radius >= rmax) & (radius < C15_BETA_MID_RADIUS_FACTOR * rmax)
        beta[inner] = C15_BETA_INNER_BASE_DEG + C15_BETA_INNER_SLOPE_DEG * (
            radius[inner] / np.broadcast_to(rmax, radius.shape)[inner]
        )
        beta[middle] = C15_BETA_MID_BASE_DEG + C15_BETA_MID_SLOPE_DEG * (
            radius[middle] / np.broadcast_to(rmax, radius.shape)[middle] - 1.0
        )

        hemisphere = np.sign(np.asarray(si_track["lat"].values, dtype=float))
        hemisphere[hemisphere == 0] = 1.0
        cta = np.arctan2(radial_direction[..., 0], radial_direction[..., 1])
        wind_angle = cta + np.deg2rad(hemisphere[:, None] * (90.0 + beta))
        rotation_direction = np.stack(
            [np.sin(wind_angle), np.cos(wind_angle)], axis=-1
        )

        vtrans = np.asarray(si_track["vtrans"].values, dtype=float)
        translation_speed = np.linalg.norm(vtrans, axis=1)
        translation_direction = np.divide(
            vtrans,
            translation_speed[:, None],
            out=np.zeros_like(vtrans),
            where=translation_speed[:, None] > 0,
        )
        vmoc = (
            weak_scale
            * translation_speed[:, None]
            * radius
            * rmax
            / (radius**2 + rmax**2)
        )
        # The profile was selected using surface_core/0.9 and is already at
        # gradient level. ADCIRC's independent 1-min-to-10-min factor (0.893)
        # must not enter TCR.
        total = (
            gradient_wind_ms[..., None] * rotation_direction
            + vmoc[..., None] * translation_direction[:, None, :]
        )

        blend = np.ones_like(radius)
        middle_blend = (
            (radius >= C15_BLEND_INNER_RADIUS_M)
            & (radius <= C15_BLEND_OUTER_RADIUS_M)
        )
        blend[middle_blend] = (
            C15_BLEND_OUTER_RADIUS_M - radius[middle_blend]
        ) / (C15_BLEND_OUTER_RADIUS_M - C15_BLEND_INNER_RADIUS_M)
        blend[radius > C15_BLEND_OUTER_RADIUS_M] = 0.0
        total *= blend[..., None]
        total[~np.asarray(mask_centr_close, dtype=bool)] = 0.0
        return total

    def compute_c15(
        si_track: xr.Dataset,
        d_centr: np.ndarray,
        mask_centr_close: np.ndarray,
        model: int,
        cyclostrophic: bool = False,
        model_kwargs: dict[str, object] | None = None,
    ) -> np.ndarray:
        if model != C15_TCR_MODEL_ID:
            return original_compute(
                si_track,
                d_centr,
                mask_centr_close,
                model,
                cyclostrophic=cyclostrophic,
                model_kwargs=model_kwargs,
            )
        return c15_gradient_profile(si_track, d_centr, mask_centr_close)

    def horizontal_winds_matlab(
        si_track: xr.Dataset,
        d_centr: dict[str, np.ndarray],
        mask_centr_close: np.ndarray,
        model: int,
        matlab_ref_mode: bool = False,
    ) -> dict[str, np.ndarray]:
        winds = original_horizontal_winds(
            si_track,
            d_centr,
            mask_centr_close,
            model,
            matlab_ref_mode=matlab_ref_mode,
        )
        if model != C15_TCR_MODEL_ID:
            return winds
        winds["matlab_total_vector"] = matlab_total_wind_vectors(
            si_track, d_centr[""], d_centr["dir"], mask_centr_close
        )
        for rstep in ["+h", "-h"]:
            winds[f"matlab_total_vector_r{rstep},t"] = matlab_total_wind_vectors(
                si_track, d_centr[rstep], d_centr["dir"], mask_centr_close
            )
        return winds

    def read_raster_sample_with_gradients_safe(
        path: str | Path,
        lat: np.ndarray,
        lon: np.ndarray,
        method: tuple[str, str] = ("linear", "linear"),
    ) -> tuple[np.ndarray, np.ndarray]:
        """Match CLIMADA raster sampling without ``rasterio.rowcol``.

        Rasterio 1.4.3 in the Windows ``sfincs_tcr`` environment can terminate
        the process inside ``rasterio.transform.rowcol``. This implementation
        reads the same padded raster window, then delegates interpolation and
        finite-difference gradients to CLIMADA's own routines. No physical or
        numerical TCR setting is changed.
        """
        import rasterio
        from affine import Affine
        from rasterio.windows import Window

        lat = np.asarray(lat, dtype=float)
        lon = np.asarray(lon, dtype=float)
        if lat.size == 0:
            return np.zeros(0), np.zeros((0, 2))
        with rasterio.open(path, "r") as src:
            transform = src.transform
            if not np.isclose(transform.b, 0.0) or not np.isclose(transform.d, 0.0):
                raise ValueError(f"Rotated rasters are unsupported in safe sampler: {path}")
            cols = (lon - transform.c) / transform.a
            rows = (lat - transform.f) / transform.e
            pad = 3
            col0 = max(int(np.floor(np.nanmin(cols))) - pad, 0)
            col1 = min(int(np.ceil(np.nanmax(cols))) + pad + 1, src.width)
            row0 = max(int(np.floor(np.nanmin(rows))) - pad, 0)
            row1 = min(int(np.ceil(np.nanmax(rows))) + pad + 1, src.height)
            if col1 <= col0 or row1 <= row0:
                return np.zeros(lat.size), np.zeros((lat.size, 2))
            window = Window(col0, row0, col1 - col0, row1 - row0)
            data = src.read(1, window=window).astype(float)
            local_transform = transform * Affine.translation(col0, row0)
            nodata = src.nodata
            crs = src.crs
        fill_value = float(nodata) if nodata is not None else 0.0
        data[~np.isfinite(data)] = fill_value
        values = tcr_module.u_coord.interp_raster_data(
            data,
            lat,
            lon,
            local_transform,
            method=method[0],
            fill_value=fill_value,
        )
        is_latlon = crs is not None and crs.to_epsg() == 4326
        gradient_data, gradient_transform = tcr_module.u_coord._raster_gradient(
            data, local_transform, latlon_to_m=is_latlon
        )
        gradients = tcr_module.u_coord.interp_raster_data(
            gradient_data,
            lat,
            lon,
            gradient_transform,
            method=method[1],
            fill_value=0.0,
        )
        return values, gradients

    def w_topo_matlab(
        si_track: xr.Dataset,
        d_centr: dict[str, np.ndarray],
        h_winds: dict[str, np.ndarray],
        centroids: np.ndarray,
        elevation_tif: str | Path | None = None,
    ) -> np.ndarray:
        if "matlab_total_vector" not in h_winds:
            return original_w_topo(
                si_track, d_centr, h_winds, centroids, elevation_tif=elevation_tif
            )
        if elevation_tif is None:
            elevation_tif = tcr_module.default_elevation_tif()
        h, h_grad = read_raster_sample_with_gradients_safe(
            elevation_tif,
            centroids[:, 0],
            centroids[:, 1],
            method=("linear", "linear"),
        )
        mask_onland = h > -1
        h_grad[~mask_onland, :] = 0.0
        h_grad_red = h_grad[None, :, :] * (
            np.clip(
                (150_000.0 - d_centr[""]) / 30_000.0, 0.2, 0.6
            )[:, :, None]
        )
        return (h_winds["matlab_total_vector"] * h_grad_red).sum(axis=-1)

    def w_frict_stretch_matlab(
        si_track: xr.Dataset,
        d_centr: dict[str, np.ndarray],
        h_winds: dict[str, np.ndarray],
        centroids: np.ndarray,
        res_radial_m: float = 2000.0,
        c_drag_tif: str | Path | None = None,
        min_c_drag: float = 0.001,
    ) -> np.ndarray:
        if "matlab_total_vector" not in h_winds:
            return original_w_frict_stretch(
                si_track,
                d_centr,
                h_winds,
                centroids,
                res_radial_m=res_radial_m,
                c_drag_tif=c_drag_tif,
                min_c_drag=min_c_drag,
            )
        if c_drag_tif is None:
            c_drag_tif = tcr_module.default_drag_tif()
        vnet = {
            f"r{rstep}h,t": np.linalg.norm(
                h_winds[f"matlab_total_vector_r{rstep}h,t"], axis=-1
            )
            for rstep in ["+", "-"]
        }
        cd, cd_grad = read_raster_sample_with_gradients_safe(
            c_drag_tif,
            centroids[:, 0],
            centroids[:, 1],
            method=("linear", "linear"),
        )
        mask_onland = cd >= min_c_drag
        cd[~mask_onland] = min_c_drag
        cd_grad[~mask_onland, :] = 0.0
        cd_hstep = (
            cd_grad[None] * (0.5 * res_radial_m * d_centr["dir"])
        ).sum(axis=-1)
        cd_values = {
            "": cd,
            "r+h": np.clip(cd + cd_hstep, 0.0, 0.01),
            "r-h": np.clip(cd - cd_hstep, 0.0, 0.01),
        }
        tau = {
            f"r{rstep}h,t": -cd_values[f"r{rstep}h"]
            * h_winds[f"r{rstep}h,t"]
            * vnet[f"r{rstep}h,t"]
            for rstep in ["+", "-"]
        }
        d_m_dr = {
            f"r{rstep}h,t": d_centr[f"{rstep}h"]
            * (
                si_track["cp"].values[:, None]
                + (1 if rstep == "+" else -1)
                * (h_winds[f"r{rstep},t"] - h_winds["r,t"])
                / res_radial_m
            )
            + h_winds[f"r{rstep}h,t"]
            for rstep in ["+", "-"]
        }
        pre_wf_wt = {
            f"r{rstep}h,t": d_centr[rstep] ** 2
            / np.fmax(10.0, d_m_dr[f"r{rstep}h,t"])
            * (
                np.fmin(
                    1.0,
                    -1.0
                    + 2.0
                    * (d_centr[rstep] / si_track["rad"].values[:, None]) ** 2,
                )
                * tcr_module.H_TROP
                * (h_winds[f"r{rstep},t+"] - h_winds[f"r{rstep},t-"])
                / (2.0 * si_track["tstep"].values[:, None])
                - tau[f"r{rstep}h,t"]
            )
            for rstep in ["+", "-"]
        }
        return (pre_wf_wt["r+h,t"] - pre_wf_wt["r-h,t"]) / (
            res_radial_m * d_centr[""]
        )

    def tctrack_to_si_matlab(track: xr.Dataset, *posargs, **kwargs) -> xr.Dataset:
        si_track = original_tctrack_to_si(track, *posargs, **kwargs)
        required = {
            "translation_velocity_east",
            "translation_velocity_north",
            "c15_vmax_gradient",
            "c15_vmax_surface_core",
            "vmax_total_surface",
            "ushear",
            "vshear",
        }
        if not required.issubset(track.variables):
            raise ValueError(
                "C15 TCR track lacks gradient Vmax, direct 250-850 hPa shear, "
                "or forward translation-vector variables"
            )
        gradient_vmax = np.asarray(
            track["c15_vmax_gradient"].values, dtype=float
        )
        surface_core_vmax = np.asarray(
            track["c15_vmax_surface_core"].values, dtype=float
        )
        expected_surface_vmax = np.asarray(
            track["vmax_total_surface"].values, dtype=float
        )
        if not np.allclose(
            np.asarray(si_track["vmax"].values, dtype=float),
            expected_surface_vmax,
            rtol=1.0e-12,
            atol=1.0e-12,
        ):
            raise ValueError(
                "C15 TCR wind-level mismatch: max_sustained_wind must equal "
                "the track's total surface vmax_trks for the q950 diagnosis"
            )
        # si_track['vmax'] remains the total surface input used by the q950
        # diagnosis. Custom C15 dynamics use the translation-free gradient
        # value derived once from that same surface input.
        si_track["vmax_gradient"] = ("time", gradient_vmax.copy())
        si_track["vmax_surface_core"] = ("time", surface_core_vmax.copy())
        vtrans = np.stack(
            [
                np.asarray(track["translation_velocity_north"].values, dtype=float),
                np.asarray(track["translation_velocity_east"].values, dtype=float),
            ],
            axis=1,
        )
        si_track["vtrans"].values[:] = vtrans
        si_track["vtrans_norm"].values[:] = np.linalg.norm(vtrans, axis=1)
        return si_track

    model_map["C15"] = C15_TCR_MODEL_ID
    tcr_module.compute_angular_windspeeds = compute_c15
    tcr_module._horizontal_winds = horizontal_winds_matlab
    tcr_module._w_topo = w_topo_matlab
    tcr_module._w_frict_stretch = w_frict_stretch_matlab
    tcr_module.tctrack_to_si = tctrack_to_si_matlab
    try:
        yield
    finally:
        tcr_module.compute_angular_windspeeds = original_compute
        tcr_module._horizontal_winds = original_horizontal_winds
        tcr_module._w_topo = original_w_topo
        tcr_module._w_frict_stretch = original_w_frict_stretch
        tcr_module.tctrack_to_si = original_tctrack_to_si
        if had_petals_surface_factor:
            tcr_module.GRADIENT_LEVEL_TO_SURFACE_WINDS = (
                original_petals_surface_factor
            )
        else:
            delattr(tcr_module, "GRADIENT_LEVEL_TO_SURFACE_WINDS")
        if had_c15:
            model_map["C15"] = old_c15
        else:
            model_map.pop("C15", None)


def build_climada_track(event: Event, args: argparse.Namespace) -> xr.Dataset | None:
    frame_all = read_track_points(event.meta_path)
    frame = frame_all[(frame_all["time"] >= event.start) & (frame_all["time"] <= event.stop)].copy()
    if len(frame) < 3:
        LOG.warning("%s has fewer than three active track points; rainfall will be zero", event.event_id)
        return None
    track_index = int(round(float(event.meta["track_index"])))
    original_idx = frame["original_time_index"].round().astype(int).to_numpy() - 1
    with xr.open_dataset(args.track_nc, decode_times=False) as raw:
        raw_track = raw.isel(n_trk=track_index - 1)
        missing_winds = [
            name
            for name in TCR_REQUIRED_TRACK_WIND_VARIABLES
            if name not in raw_track.variables
        ]
        if missing_winds:
            raise ValueError(
                f"{args.track_nc} lacks TCR environmental-wind variables: "
                + ", ".join(missing_winds)
            )
        raw_surface_vmax_historical = raw_track["vmax_trks"].isel(
            time=xr.DataArray(original_idx, dims="points")
        ).values
        u250 = raw_track["u250_trks"].isel(
            time=xr.DataArray(original_idx, dims="points")
        ).values
        v250 = raw_track["v250_trks"].isel(
            time=xr.DataArray(original_idx, dims="points")
        ).values
        u850 = raw_track["u850_trks"].isel(time=xr.DataArray(original_idx, dims="points")).values
        v850 = raw_track["v850_trks"].isel(time=xr.DataArray(original_idx, dims="points")).values
        raw_basin_value = raw_track["tc_basins"].values.item() if "tc_basins" in raw_track else None
    u250 = np.asarray(u250, dtype=float)
    v250 = np.asarray(v250, dtype=float)
    u850 = np.asarray(u850, dtype=float)
    v850 = np.asarray(v850, dtype=float)
    ushear = u250 - u850
    vshear = v250 - v850
    raw_surface_vmax_historical = np.asarray(
        raw_surface_vmax_historical, dtype=float
    )
    required_tcr_winds = np.stack(
        [u250, v250, u850, v850, ushear, vshear, raw_surface_vmax_historical]
    )
    if not np.all(np.isfinite(required_tcr_winds)):
        raise ValueError(
            f"Non-finite vmax_trks or 250/850-hPa environmental winds for {event.event_id}"
        )
    times = pd.DatetimeIndex(frame["time"])
    lon = frame["lon180"].to_numpy(dtype=float)
    lat = frame["lat"].to_numpy(dtype=float)
    vmax = frame["vmax_ms"].to_numpy(dtype=float)
    pressure = frame["pressure_hPa"].to_numpy(dtype=float)
    # Future P4 metadata retains both historical and future surface Vmax. The
    # historical column must trace back to the source NetCDF vmax_trks; the
    # current ``vmax`` is then used directly (historical or future) throughout
    # the TCR interface.
    intensity_ratio = np.ones_like(vmax)
    historical_surface_vmax = vmax.copy()
    if {"historical_vmax_ms", "future_vmax_ms"}.issubset(frame.columns):
        historical_vmax = frame["historical_vmax_ms"].to_numpy(dtype=float)
        future_vmax = frame["future_vmax_ms"].to_numpy(dtype=float)
        historical_surface_vmax = historical_vmax
        intensity_ratio = np.divide(
            future_vmax,
            historical_vmax,
            out=np.ones_like(future_vmax),
            where=np.isfinite(historical_vmax) & (historical_vmax > 0.0),
        )
        if not np.all(np.isfinite(intensity_ratio)) or np.any(intensity_ratio <= 0.0):
            raise ValueError(f"Invalid future intensity ratio for {event.event_id}")
    if not np.allclose(
        historical_surface_vmax,
        raw_surface_vmax_historical,
        rtol=1.0e-6,
        atol=1.0e-6,
    ):
        max_error = float(
            np.nanmax(np.abs(historical_surface_vmax - raw_surface_vmax_historical))
        )
        raise ValueError(
            f"Historical surface Vmax is inconsistent with source vmax_trks for "
            f"{event.event_id}; maximum absolute error={max_error:.6g} m/s"
        )
    t600 = sample_monthly_t600(
        args.t600_source, times, lon, lat
    )
    translation_east, translation_north, translation = matlab_track_step_motion(
        lon, lat, times
    )
    use_c15 = args.rain_model == "tcr" and args.tcr_wind_model == "c15"
    c15_params: dict[str, np.ndarray] | None = None
    cv: float | None = None
    if use_c15:
        if "Cv_used" not in event.meta:
            raise ValueError(
                f"Missing Cv_used in {event.meta_path}; refusing the previous "
                "silent Cv=1 fallback"
            )
        cv = float(event.meta["Cv_used"])
        c15_params = c15_parameters_from_track(
            vmax, translation, pressure, lat, cv
        )
        # q950 receives the unchanged surface vmax_trks. The C15 dynamics use
        # the translation-free gradient core derived from that same value.
        gradient_core_vmax = np.maximum(c15_params["vm_raw"], 0.0)
        rain_vmax = vmax
        pressure_for_tcr = c15_params["pressure_model_hpa"]
        rmax_m = c15_params["rmax_m"]
        LOG.debug(
            "%s C15: Cv=%.6f, raw B %.3f..%.3f, clipped B %.3f..%.3f",
            event.event_id,
            cv,
            float(np.nanmin(c15_params["holland_b_raw"])),
            float(np.nanmax(c15_params["holland_b_raw"])),
            float(np.nanmin(c15_params["holland_b"])),
            float(np.nanmax(c15_params["holland_b"])),
        )
    else:
        rain_vmax = vmax
        pressure_for_tcr = pressure
        fallback_cv = float(event.meta.get("Cv_used", "1.0"))
        fallback_gradient_core = np.maximum(
            vmax - np.maximum(translation, 0.0), 0.0
        ) / TCR_GRADIENT_TO_SURFACE_WIND_FACTOR
        fallback_vm = np.clip(
            np.maximum(fallback_gradient_core, C15_MIN_LOOKUP_WIND_MS),
            C15_MIN_LOOKUP_WIND_MS,
            C15_MAX_LOOKUP_WIND_MS,
        )
        rmax_m = np.maximum(
            fallback_cv
            * 51.6
            * np.exp(-0.0223 * fallback_vm + 0.0281 * np.abs(lat))
            * 1000.0,
            C15_MIN_RMAX_M,
        )
    dt_hours = np.diff(times.values).astype("timedelta64[s]").astype(float) / 3600.0
    dt_hours = np.r_[dt_hours, dt_hours[-1]]
    basin = basin_from_event(event, str(raw_basin_value) if raw_basin_value is not None else None)
    sid = event.meta.get("track_id", event.event_id)
    track_vars: dict[str, tuple[str, np.ndarray]] = {
        "lat": ("time", lat),
        "lon": ("time", lon),
        "time_step": ("time", dt_hours),
        "radius_max_wind": ("time", rmax_m / 1852.0),
        "radius_oci": ("time", np.full(len(frame), np.nan)),
        "max_sustained_wind": ("time", rain_vmax),
        "central_pressure": ("time", pressure_for_tcr),
        "environmental_pressure": (
            "time",
            np.full(len(frame), C15_ENVIRONMENTAL_PRESSURE_HPA),
        ),
        "basin": ("time", np.full(len(frame), basin, dtype="U2")),
        # Do not add q950 here. CLIMADA-Petals gives q950 precedence over
        # t600; omitting q950 is therefore required for TCR to diagnose the
        # 950-hPa saturation specific humidity from t600 and instantaneous
        # si_track["vmax"] at every track point.
        "t600": ("time", t600),
        # Direct environmental vertical shear. Presence of ushear/vshear
        # makes CLIMADA use u250-u850 and v250-v850 instead of BAM inference.
        "ushear": ("time", ushear),
        "vshear": ("time", vshear),
        "u250": ("time", u250),
        "v250": ("time", v250),
        "u850": ("time", u850),
        "v850": ("time", v850),
        "vmax_total_surface": ("time", vmax),
        "vmax_historical_surface": ("time", historical_surface_vmax),
        "vmax_source_track_surface": ("time", raw_surface_vmax_historical),
        "future_intensity_ratio": ("time", intensity_ratio),
        "translation_speed": ("time", translation),
        "translation_velocity_east": ("time", translation_east),
        "translation_velocity_north": ("time", translation_north),
        "central_pressure_wpr": ("time", pressure),
    }
    if c15_params is not None:
        track_vars.update(
            {
                "c15_vm_raw": ("time", c15_params["vm_raw"]),
                "c15_vmax_gradient": ("time", gradient_core_vmax),
                "c15_vm_lookup": ("time", c15_params["vm_lookup"]),
                "c15_vmax_surface_core": ("time", c15_params["surface_core"]),
                "c15_surface_vm_lookup": (
                    "time",
                    c15_params["surface_vm_lookup"],
                ),
                "c15_weak_scale": ("time", c15_params["weak_scale"]),
                "holland_b_raw": ("time", c15_params["holland_b_raw"]),
                "holland_b": ("time", c15_params["holland_b"]),
            }
        )
    track = xr.Dataset(
        data_vars=track_vars,
        coords={"time": times},
        attrs={
            "max_sustained_wind_unit": "m/s",
            "central_pressure_unit": "mb",
            "name": sid,
            "sid": sid,
            "orig_event_flag": True,
            "data_provider": "ERA5 synthetic track enriched from fort22_meta",
            "t600_source": str(args.t600_source),
            "tcr_humidity_pathway": (
                "CLIMADA _qs_from_t_diff_level(t600, total surface vmax_trks, "
                "600 hPa, 950 hPa); q950 intentionally absent"
            ),
            "tcr_wind_profile": (
                "C15_PROFILE_SELECTED_BY_SURFACE_CORE_DIVIDED_BY_0.9"
                if use_c15
                else args.tcr_wind_model.upper()
            ),
            "tcr_dynamic_wind_level": (
                "gradient" if use_c15 else "model_default"
            ),
            "tcr_q950_vmax_level": (
                "track_total_surface" if use_c15 else "track_default"
            ),
            "tcr_gradient_to_surface_wind_factor": (
                TCR_GRADIENT_TO_SURFACE_WIND_FACTOR if use_c15 else np.nan
            ),
            "tcr_vertical_shear": "direct_u250_minus_u850_v250_minus_v850",
            "tcr_vertical_shear_levels_hpa": "250,850",
            "c15_lookup_wind_level": (
                "gradient_peak_from_surface_core_divided_by_0.9"
                if use_c15
                else "not_used"
            ),
            "c15_second_profile_division_by_0.9": False,
            "c15_core_wind_rule": (
                "C15 lookup peak=max(total surface vmax_trks-translation "
                "speed,0)/"
                f"{TCR_GRADIENT_TO_SURFACE_WIND_FACTOR:g}; lookup profile "
                "is used directly as TCR gradient wind"
                if use_c15
                else "not_used"
            ),
            "c15_rmax_wind_level": (
                "translation-free surface rotational wind, matching fort.22"
                if use_c15
                else "not_used"
            ),
            "c15_holland_b_wind": (
                "total surface vmax_trks, matching fort.22"
                if use_c15
                else "not_used"
            ),
            "c15_translation_rule": (
                "forward current-to-next vector; final point zero" if use_c15 else "not_used"
            ),
            "c15_asymmetry_rule": (
                "translation removed before surface-to-gradient conversion, then "
                "added once as radius-decaying vmoc; radius-dependent beta inflow angle"
                if use_c15
                else "not_used"
            ),
            "c15_outer_blend_km": "500..700" if use_c15 else "not_used",
            "c15_cv_used": float(cv) if cv is not None else np.nan,
            "id_no": track_index,
            "category": saffir_category(float(np.nanmax(vmax))),
        },
    )
    required_values = track[
        [
            "lat",
            "lon",
            "max_sustained_wind",
            "central_pressure",
            "ushear",
            "vshear",
        ]
    ].to_array().values
    if not np.all(np.isfinite(required_values)):
        raise ValueError(f"Non-finite required CLIMADA track variables for {event.event_id}")
    return track


def load_climada_petals_tcr_module() -> object:
    """Load ``tc_rainfield`` without importing unrelated forecast plugins.

    Some CLIMADA-Petals releases import ``tc_tracks_forecast`` from the
    ``climada_petals.hazard`` package initializer. That optional forecast
    module requires the native ecCodes library even though TCR rainfall does
    not use GRIB data. Loading the TCR source module directly keeps the build
    independent of that unrelated optional dependency.
    """
    module_name = "_sfincs_climada_petals_tc_rainfield"
    if module_name in sys.modules:
        return sys.modules[module_name]
    petals_spec = importlib.util.find_spec("climada_petals")
    if petals_spec is None or not petals_spec.submodule_search_locations:
        raise ImportError("climada_petals is not installed")
    module_path = (
        Path(next(iter(petals_spec.submodule_search_locations)))
        / "hazard"
        / "tc_rainfield.py"
    )
    if not module_path.is_file():
        raise ImportError(f"Cannot locate CLIMADA-Petals TCR module: {module_path}")
    module_spec = importlib.util.spec_from_file_location(module_name, module_path)
    if module_spec is None or module_spec.loader is None:
        raise ImportError(f"Cannot create module specification for {module_path}")
    module = importlib.util.module_from_spec(module_spec)
    sys.modules[module_name] = module
    try:
        module_spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(module_name, None)
        raise
    return module


def configure_conda_geospatial_data_paths() -> None:
    """Supply GDAL/PROJ data paths when Windows conda activation omitted them."""
    conda_share = Path(sys.prefix) / "Library" / "share"
    candidates = {
        "GDAL_DATA": conda_share / "gdal",
        "PROJ_LIB": conda_share / "proj",
    }
    for variable, path in candidates.items():
        if not os.environ.get(variable) and path.is_dir():
            os.environ[variable] = str(path)


def compute_rainfall(
    track: xr.Dataset | None,
    model: str,
    target_times: pd.DatetimeIndex,
    lon: np.ndarray,
    lat: np.ndarray,
    args: argparse.Namespace,
) -> np.ndarray:
    out = np.zeros((len(target_times),) + lon.shape, dtype=np.float32)
    if track is None:
        return out
    configure_conda_geospatial_data_paths()
    try:
        from climada.hazard import Centroids, TCTracks
        tcr_module = load_climada_petals_tcr_module()
    except (ImportError, ModuleNotFoundError) as exc:
        raise RuntimeError(
            "CLIMADA rainfall requested but climada/climada_petals is not installed. "
            "Use a Python 3.11 conda environment with CLIMADA-Petals 6.0.x."
        ) from exc
    centroids = Centroids.from_lat_lon(lat.ravel(), lon.ravel())
    tracks = TCTracks(data=[track])
    model_kwargs: dict[str, object] = {}
    if args.tcr_elevation_tif is not None:
        model_kwargs["elevation_tif"] = str(args.tcr_elevation_tif)
    if args.tcr_drag_tif is not None:
        model_kwargs["c_drag_tif"] = str(args.tcr_drag_tif)
    if model == "tcr":
        model_kwargs["wind_model"] = (
            "C15" if args.tcr_wind_model == "c15" else "ER11"
        )
    wind_context = (
        install_c15_tcr_wind_model(tcr_module, args)
        if model == "tcr" and args.tcr_wind_model == "c15"
        else nullcontext()
    )
    with wind_context:
        rain = tcr_module.TCRain.from_tracks(
            tracks,
            centroids=centroids,
            model="TCR" if model == "tcr" else "R-CLIPER",
            model_kwargs=model_kwargs,
            ignore_distance_to_coast=True,
            store_rainrates=True,
            intensity_thres=0.0,
            max_dist_eye_km=args.rain_max_eye_distance_km,
            max_memory_gb=args.rain_max_memory_gb,
        )
    if not rain.rainrates:
        raise RuntimeError("CLIMADA returned no time-dependent rain rates")
    native = np.asarray(rain.rainrates[0].toarray(), dtype=np.float32).reshape(
        (track.sizes["time"],) + lon.shape
    )
    target_lookup = {pd.Timestamp(time): i for i, time in enumerate(target_times)}
    for source_i, time in enumerate(pd.DatetimeIndex(track.time.values)):
        target_i = target_lookup.get(pd.Timestamp(time))
        if target_i is not None:
            out[target_i] = native[source_i]
    return out


def preflight(args: argparse.Namespace, events: Sequence[Event], dry_run: bool) -> list[str]:
    issues: list[str] = []
    for path, label in (
        (args.model_root, "SFINCS model root"),
        (args.membership_csv, "ADCIRC membership CSV"),
        (args.adcirc_result_root, "ADCIRC return-period result root"),
        (args.adcirc_run_root, "ADCIRC run root"),
        (args.cama_root, "CaMa output root"),
    ):
        if not path.exists():
            issues.append(f"Missing {label}: {path}")
    if args.rain_model != "none":
        if not args.track_nc.is_file():
            issues.append(f"Missing synthetic-track NetCDF: {args.track_nc}")
        else:
            try:
                with xr.open_dataset(args.track_nc, decode_times=False) as track_ds:
                    missing_winds = [
                        name
                        for name in TCR_REQUIRED_TRACK_WIND_VARIABLES
                        if name not in track_ds.variables
                    ]
                if missing_winds:
                    issues.append(
                        "Synthetic-track NetCDF lacks direct TCR shear inputs: "
                        + ", ".join(missing_winds)
                    )
            except (OSError, ValueError) as exc:
                issues.append(
                    f"Cannot validate synthetic-track environmental winds: {exc}"
                )
        if args.rain_model == "tcr":
            if args.t600_source.is_dir():
                required_years: set[int] = set()
                for event in events:
                    required_years.update(
                        range(event.start.year, event.stop.year + 1)
                    )
                missing_t600 = [
                    args.t600_source / f"era5_t600_monthly_{year}.nc"
                    for year in sorted(required_years)
                    if not (
                        args.t600_source
                        / f"era5_t600_monthly_{year}.nc"
                    ).is_file()
                ]
                if missing_t600:
                    issues.append(
                        "Missing ERA5 t600 annual files: "
                        + "; ".join(str(path) for path in missing_t600[:8])
                    )
            elif not args.t600_source.is_file():
                issues.append(f"Missing TCR t600 source: {args.t600_source}")
        if args.rain_model == "tcr" and args.tcr_wind_model == "c15":
            if not args.c15_predata_dir.is_dir():
                issues.append(f"Missing C15 lookup-table directory: {args.c15_predata_dir}")
            else:
                try:
                    keys = c15_available_keys(str(args.c15_predata_dir.resolve()))
                    LOG.info(
                        "C15 TCR wind library: %s profiles from %s",
                        len(keys),
                        args.c15_predata_dir,
                    )
                except (FileNotFoundError, ValueError) as exc:
                    issues.append(str(exc))
        missing_modules = [
            name for name in ("climada", "climada_petals") if importlib.util.find_spec(name) is None
        ]
        if missing_modules:
            message = "Missing optional rainfall packages: " + ", ".join(missing_modules)
            if dry_run:
                issues.append(message + " (production build will stop)")
            else:
                raise RuntimeError(message)
    if args.met_forcing == "fort22":
        missing = [str(event.fort22) for event in events if not event.fort22.is_file()]
        if missing:
            issues.append("Missing fort.22 for events: " + "; ".join(missing[:5]))
    if args.cama_quantile_mode == "tc_summary_p50":
        if not args.cama_flow_summary_csv.is_file():
            issues.append(f"Missing GTC diagnostic flow summary: {args.cama_flow_summary_csv}")
        if not args.cama_inlet_cells_csv.is_file():
            issues.append(f"Missing GTC diagnostic inlet table: {args.cama_inlet_cells_csv}")
    if args.cama_quantile_mode in {
        "tc_inlet_p50", "tc_inlet_p90", "tc_monthly_mean", "tc_weighted"
    }:
        if not args.cama_month_weights_csv.is_file():
            issues.append(f"Missing GTC TC-month weight table: {args.cama_month_weights_csv}")
    if args.cama_quantile_mode == "tc_monthly_mean":
        if not args.cama_monthly_inlet_csv.is_file():
            issues.append(
                "Missing P1_plot_sfincs_inflow_diagnostics.py monthly inlet table: "
                f"{args.cama_monthly_inlet_csv}"
            )
    if args.cama_quantile_mode in {"tc_inlet_p50", "tc_inlet_p90", "tc_weighted"}:
        if not args.cama_diagnostic_cache_dir.is_dir():
            issues.append(
                "Missing P1_plot_sfincs_inflow_diagnostics.py daily cache directory: "
                f"{args.cama_diagnostic_cache_dir}"
            )
    return issues


def validate_options(args: argparse.Namespace) -> None:
    if args.dry_run and args.update_existing_inflow_only:
        raise ValueError("--dry-run and --update-existing-inflow-only are mutually exclusive")
    exclusive_modes = sum(
        bool(value)
        for value in (
            args.update_existing_inflow_only,
            args.prepare_boundary_maps_only,
            args.prepare_flow_cache_only,
        )
    )
    if exclusive_modes > 1:
        raise ValueError(
            "--update-existing-inflow-only, --prepare-boundary-maps-only and "
            "--prepare-flow-cache-only are mutually exclusive"
        )
    if args.buffer_days < 0 or args.hotstart_hours < 0:
        raise ValueError("buffer-days and hotstart-hours must be non-negative")
    if args.forcing_dt_seconds <= 0 or args.cama_chunk_days <= 0:
        raise ValueError("forcing-dt-seconds and cama-chunk-days must be positive")
    if (
        args.map_output_interval_seconds is not None
        and args.map_output_interval_seconds < 0
    ):
        raise ValueError("map-output-interval-seconds must be non-negative")
    if (
        args.max_output_interval_seconds is not None
        and args.max_output_interval_seconds <= 0
    ):
        raise ValueError("max-output-interval-seconds must be positive when specified")
    if args.adcirc_neighbors <= 0:
        raise ValueError("adcirc-neighbors must be positive")
    if not np.isfinite(args.adcirc_missing_waterlevel_fill_m):
        raise ValueError("adcirc-missing-waterlevel-fill-m must be finite")
    if not 0.0 <= args.cama_quantile <= 1.0:
        raise ValueError("cama-quantile must be in [0, 1]")
    if args.cama_quantile_mode == "tc_inlet_p50" and not np.isclose(
        args.cama_quantile, 0.5
    ):
        raise ValueError("tc_inlet_p50 requires --cama-quantile 0.5")
    if args.cama_quantile_mode == "tc_inlet_p90" and not np.isclose(
        args.cama_quantile, 0.9
    ):
        raise ValueError("tc_inlet_p90 requires --cama-quantile 0.9")
    if args.cama_start_year > args.cama_end_year:
        raise ValueError("cama-start-year must not exceed cama-end-year")
    if not 0.0 <= args.min_water_valid_fraction <= 1.0:
        raise ValueError("min-water-valid-fraction must be in [0, 1]")
    if not 0.0 <= args.min_cama_valid_fraction <= 1.0:
        raise ValueError("min-cama-valid-fraction must be in [0, 1]")


def event_target_seconds(event: Event, dt_seconds: int) -> np.ndarray:
    if event.duration_seconds % dt_seconds != 0:
        raise ValueError(
            f"Event duration {event.duration_seconds}s is not divisible by forcing dt {dt_seconds}s"
        )
    return np.arange(0, event.duration_seconds + dt_seconds, dt_seconds, dtype=float)


def map_output_intervals(
    args: argparse.Namespace, event: Event
) -> tuple[int, int]:
    dtout = (
        int(round(event.duration_seconds))
        if args.map_output_interval_seconds is None
        else int(args.map_output_interval_seconds)
    )
    dtmaxout = (
        int(round(event.duration_seconds))
        if args.max_output_interval_seconds is None
        else int(args.max_output_interval_seconds)
    )
    if dtmaxout > event.duration_seconds:
        raise ValueError(
            f"max-output-interval-seconds={dtmaxout} exceeds event duration "
            f"{event.duration_seconds:.0f}s for {event.event_id}"
        )
    return dtout, dtmaxout


def build_one_model(
    args: argparse.Namespace,
    event: Event,
    context: ModelContext,
    target_seconds: np.ndarray,
    selected_node_ids: np.ndarray,
    fort63_times: np.ndarray,
    fort63_values: np.ndarray,
    q_table: pd.DataFrame,
    fort22: Fort22Data | None,
    climada_track: xr.Dataset | None,
) -> dict[str, object]:
    target = args.output_root / event.block_id / event.event_id / context.model_id
    if not prepare_output_model(
        context, target, args.output_root, args.overwrite, args.copy_mode
    ):
        return {"status": "skipped", "target": str(target)}
    try:
        dtout, dtmaxout = map_output_intervals(args, event)
        water, _ = select_boundary_waterlevels(
            context,
            selected_node_ids,
            fort63_times,
            fort63_values,
            target_seconds,
            args.buffer_days * 86400.0,
            args.hotstart_hours * 3600.0,
            args.adcirc_missing_waterlevel_fill_m,
        )
        write_table_forcing(target / "sfincs.bzs", target_seconds, water, decimals=4)

        discharge, _ = discharge_for_context(context, q_table, target_seconds)
        if discharge is not None:
            install_shared_table_forcing(
                args,
                context,
                "sfincs.dis",
                target_seconds,
                discharge,
                decimals=6,
                destination=target / "sfincs.dis",
            )

        x, y, lon, lat = forcing_grid(context, args.forcing_resolution_m)
        if fort22 is not None:
            write_fort22_meteorology(target, fort22, x, y, lon, lat)
        if args.rain_model != "none":
            target_times = pd.DatetimeIndex(
                event.start + pd.to_timedelta(target_seconds, unit="s")
            )
            rainfall = compute_rainfall(
                climada_track, args.rain_model, target_times, lon, lat, args
            )
            write_grid_netcdf(
                target / "sfincs_precipitation.nc",
                target_times,
                x,
                y,
                {"Precipitation": rainfall},
                {"Precipitation": "mm hr-1"},
                (
                    f"CLIMADA TCR-{args.tcr_wind_model.upper()} TC rainfall for SFINCS"
                    if args.rain_model == "tcr"
                    else f"CLIMADA {args.rain_model.upper()} TC rainfall for SFINCS"
                ),
            )

        updates: dict[str, str | int | float | None] = {
            "tref": event.start.strftime("%Y%m%d %H%M%S"),
            "tstart": event.start.strftime("%Y%m%d %H%M%S"),
            "tstop": event.stop.strftime("%Y%m%d %H%M%S"),
            "tspinup": 0.0,
            "dtout": dtout,
            "dtmaxout": dtmaxout,
            "outputformat": "net",
            "bndfile": "sfincs.bnd",
            "bzsfile": "sfincs.bzs",
            "srcfile": "sfincs.src" if discharge is not None else None,
            "disfile": "sfincs.dis" if discharge is not None else None,
            "wndfile": None,
            "spwfile": None,
            "netamuamvfile": "sfincs_wind.nc" if fort22 is not None else None,
            "netampfile": "sfincs_pressure.nc" if fort22 is not None else None,
            "wind": 1 if fort22 is not None else 0,
            "baro": 1 if fort22 is not None else 0,
            "dtwnd": args.forcing_dt_seconds,
            "amprfile": None,
            "netamprfile": "sfincs_precipitation.nc" if args.rain_model != "none" else None,
            "ampr_block": 1 if args.rain_model != "none" else None,
        }
        update_sfincs_inp(target / "sfincs.inp", updates)
        marker = {
            "status": "complete",
            "return_period_years": args.return_period,
            "flow_quantile_nonexceedance": args.cama_quantile,
            "flow_quantile_mode": args.cama_quantile_mode,
            "block_id": event.block_id,
            "event_id": event.event_id,
            "model_id": context.model_id,
            "window_start": str(event.start),
            "window_end": str(event.stop),
            "map_output_interval_seconds": dtout,
            "max_output_interval_seconds": dtmaxout,
            "map_output_mode": (
                "event_maximum_plus_final_snapshot"
                if dtout == dtmaxout == int(round(event.duration_seconds))
                else "instantaneous_and_maximum"
            ),
            "water_buffer_hours": args.buffer_days * 24.0,
            "water_hotstart_ramp_hours": args.hotstart_hours,
            "fort63_first_seconds": float(fort63_times[0]),
            "fort63_last_seconds": float(fort63_times[-1]),
            "adcirc_valid_search_extra_km": ADCIRC_LOCAL_SEARCH_EXTRA_KM,
            "adcirc_missing_waterlevel_fill_m": args.adcirc_missing_waterlevel_fill_m,
            "cama_definition": getattr(
                args,
                "cama_definition",
                "CaMa-Flood routed-outflow quantile",
            ),
            "climate_scenario": getattr(args, "climate_scenario", None),
            "t600_source": str(args.t600_source),
            "tcr_humidity_pathway": (
                "t600_plus_surface_equivalent_vmax_to_saturation_q950"
                if args.rain_model == "tcr"
                else None
            ),
            "q950_track_variable_present": (
                bool(climada_track is not None and "q950" in climada_track)
                if args.rain_model == "tcr"
                else None
            ),
            "cama_models": getattr(args, "cama_models_used", None),
            "rain_model": args.rain_model,
            "tcr_wind_model": (
                args.tcr_wind_model if args.rain_model == "tcr" else None
            ),
            "c15_predata_dir": (
                str(args.c15_predata_dir)
                if args.rain_model == "tcr" and args.tcr_wind_model == "c15"
                else None
            ),
            "c15_lookup_wind_level": (
                "gradient_peak_from_surface_core_divided_by_0.9"
                if args.rain_model == "tcr" and args.tcr_wind_model == "c15"
                else None
            ),
            "c15_profile_conversion_inside_tcr": (
                "lookup_peak=surface_core/0.9; C15 profile used directly"
                if args.rain_model == "tcr" and args.tcr_wind_model == "c15"
                else None
            ),
            "tcr_gradient_to_surface_wind_factor": (
                TCR_GRADIENT_TO_SURFACE_WIND_FACTOR
                if args.rain_model == "tcr" and args.tcr_wind_model == "c15"
                else None
            ),
            "tcr_vertical_shear": (
                "direct_u250_minus_u850_v250_minus_v850"
                if args.rain_model == "tcr"
                else None
            ),
            "c15_rmax_out_of_range": (
                args.c15_rmax_out_of_range
                if args.rain_model == "tcr" and args.tcr_wind_model == "c15"
                else None
            ),
            "c15_cv_used": (
                float(climada_track.attrs["c15_cv_used"])
                if climada_track is not None
                and args.rain_model == "tcr"
                and args.tcr_wind_model == "c15"
                else None
            ),
            "c15_core_wind_rule": (
                str(climada_track.attrs["c15_core_wind_rule"])
                if climada_track is not None
                and args.rain_model == "tcr"
                and args.tcr_wind_model == "c15"
                else None
            ),
            "holland_b_raw_min": (
                float(climada_track["holland_b_raw"].min().item())
                if climada_track is not None and "holland_b_raw" in climada_track
                else None
            ),
            "holland_b_raw_max": (
                float(climada_track["holland_b_raw"].max().item())
                if climada_track is not None and "holland_b_raw" in climada_track
                else None
            ),
            "holland_b_clipped_min": (
                float(climada_track["holland_b"].min().item())
                if climada_track is not None and "holland_b" in climada_track
                else None
            ),
            "holland_b_clipped_max": (
                float(climada_track["holland_b"].max().item())
                if climada_track is not None and "holland_b" in climada_track
                else None
            ),
            "meteorology": args.met_forcing,
            "adcirc_neighbor_candidates": args.adcirc_neighbors,
            "source_fort63": str(event.fort63),
            "source_fort14": str(event.fort14),
            "source_fort22_meta": str(event.meta_path),
        }
        if not runtime_model_complete(target):
            missing = [
                str(reference)
                for reference in sfincs_file_references(target / "sfincs.inp").values()
                if not (
                    reference if reference.is_absolute() else target / reference
                ).is_file()
            ]
            raise RuntimeError(
                f"Generated SFINCS model still lacks referenced inputs: {missing}"
            )
        return {"status": "complete", "target": str(target), **marker}
    except Exception:
        # Do not leave diagnostic markers or half-built models in the runtime tree.
        if target.exists() and is_within(target, args.output_root):
            shutil.rmtree(target)
        raise


def run(args: argparse.Namespace) -> int:
    resolve_runtime_paths(args)
    validate_options(args)
    if not hasattr(args, "cama_definition"):
        if args.cama_quantile_mode in {"tc_inlet_p50", "tc_inlet_p90"}:
            args.cama_definition = (
                "Independent per-inlet TC-landfall-month-probability-weighted "
                f"P{int(round(args.cama_quantile * 100)):02d}; "
                "source points on the same inlet section are equally split"
            )
        elif args.cama_quantile_mode == "tc_summary_p50":
            args.cama_definition = (
                "Direct tc_p50_flow_m3s from P1_plot_sfincs_inflow_diagnostics.py "
                "summary; inlet allocation by q_reference_m3s fraction"
            )
        elif args.cama_quantile_mode == "tc_monthly_mean":
            args.cama_definition = (
                "GTC-specific sum of TC-landfall month probability times "
                "diagnostic monthly mean routed outflow"
            )
        elif args.cama_quantile_mode == "tc_weighted":
            args.cama_definition = (
                f"{args.cama_start_year}-{args.cama_end_year} GTC-specific "
                "TC-landfall-month-probability-weighted total routed-outflow "
                f"Q{int(round(args.cama_quantile * 100)):02d}"
            )
        else:
            args.cama_definition = (
                f"{args.cama_start_year}-{args.cama_end_year} all-daily routed-outflow "
                f"non-exceedance Q{int(round(args.cama_quantile * 100)):02d}"
            )
    membership = load_membership(args)
    selected_blocks = set(membership["block_id"])
    contexts = load_model_contexts(args, membership)
    model_count = sum(len(items) for items in contexts.values())
    src_model_count = sum(
        not context.src_mapping.empty for items in contexts.values() for context in items
    )
    source_count = sum(
        len(context.src_mapping) for items in contexts.values() for context in items
    )
    cells = collect_cama_cells(contexts)
    if args.prepare_boundary_maps_only:
        prepare_boundary_maps(args, contexts)
        LOG.info(
            "Prepared persistent ADCIRC-to-SFINCS maps for %s block(s); cache=%s",
            len(contexts),
            args.boundary_map_cache_dir,
        )
        return 0
    if args.prepare_flow_cache_only:
        args.output_root.mkdir(parents=True, exist_ok=True)
        q_table = load_or_compute_cama_quantile(
            args, cells, allow_compute=True, contexts=contexts
        )
        LOG.info(
            "Prepared complete CaMa-Flood cache for %s selected GTC model(s): %s rows",
            model_count,
            len(q_table),
        )
        return 0
    if args.update_existing_inflow_only:
        args.output_root.mkdir(parents=True, exist_ok=True)
        q_table = load_or_compute_cama_quantile(
            args, cells, allow_compute=True, contexts=contexts
        )
        return update_existing_inflow_only(args, contexts, q_table)

    events = discover_events(args, selected_blocks)
    LOG.info(
        "Inventory: %s mapped rows, %s runnable models in %s blocks, %s events, "
        "%s river models/%s src points/%s unique CaMa cells",
        len(membership),
        model_count,
        len(contexts),
        len(events),
        src_model_count,
        source_count,
        len(cells),
    )
    if args.max_output_interval_seconds is None:
        selected_durations = sorted(
            {int(round(event.duration_seconds)) for event in events}
        )
        max_output_description = (
            "full event"
            if not selected_durations
            else "full event (" + ", ".join(map(str, selected_durations)) + " s)"
        )
    else:
        max_output_description = f"{args.max_output_interval_seconds} s"
    if args.map_output_interval_seconds is None:
        map_output_description = max_output_description
    else:
        map_output_description = f"{args.map_output_interval_seconds} s"
    LOG.info(
        "SFINCS map output: dtout=%s; dtmaxout=%s",
        map_output_description,
        max_output_description,
    )
    issues = preflight(args, events, args.dry_run)
    for issue in issues:
        LOG.warning(issue)
    if args.dry_run:
        return 0
    if issues:
        raise RuntimeError("Preflight failed: " + "; ".join(issues))
    if not events:
        raise RuntimeError("No ADCIRC events selected")

    args.output_root.mkdir(parents=True, exist_ok=True)
    q_table = load_or_compute_cama_quantile(
        args, cells, allow_compute=True, contexts=contexts
    )
    summary: list[dict[str, object]] = []
    boundary_map_cache: dict[str, np.ndarray] = {}
    for event_index, event in enumerate(events, start=1):
        models = contexts.get(event.block_id, [])
        if not models:
            LOG.warning("No runnable SFINCS models for %s", event.block_id)
            continue
        pending = []
        for context in models:
            target = args.output_root / event.block_id / event.event_id / context.model_id
            if not runtime_model_complete(target) or args.overwrite:
                pending.append(context)
        if not pending:
            LOG.info("[%s/%s] Event already complete: %s", event_index, len(events), event.event_id)
            continue
        LOG.info(
            "[%s/%s] %s / %s -> %s model(s)",
            event_index,
            len(events),
            event.block_id,
            event.event_id,
            len(pending),
        )
        try:
            if event.block_id not in boundary_map_cache:
                boundary_map_cache[event.block_id] = load_or_build_boundary_candidates(
                    args,
                    event.block_id,
                    contexts[event.block_id],
                )
            wanted = boundary_map_cache[event.block_id]
            fort63_times, fort63_values = read_fort63_selected(event.fort63, wanted)
            target_seconds = event_target_seconds(event, args.forcing_dt_seconds)
            expected_buffer = args.buffer_days * 86400.0
            if abs(fort63_times[0] - expected_buffer) > 2 * args.forcing_dt_seconds:
                LOG.warning(
                    "fort.63 first time %.0fs differs from requested buffer %.0fs",
                    fort63_times[0],
                    expected_buffer,
                )
            fort22 = read_fort22(event) if args.met_forcing == "fort22" else None
            climada_track = (
                build_climada_track(event, args) if args.rain_model != "none" else None
            )
            for context in pending:
                try:
                    result = build_one_model(
                        args,
                        event,
                        context,
                        target_seconds,
                        wanted,
                        fort63_times,
                        fort63_values,
                        q_table,
                        fort22,
                        climada_track,
                    )
                    result.update(
                        {"block_id": event.block_id, "event_id": event.event_id, "model_id": context.model_id}
                    )
                    summary.append(result)
                except Exception as exc:
                    LOG.exception("Failed %s/%s/%s", event.block_id, event.event_id, context.model_id)
                    summary.append(
                        {
                            "status": "failed",
                            "block_id": event.block_id,
                            "event_id": event.event_id,
                            "model_id": context.model_id,
                            "error": repr(exc),
                        }
                    )
                    if not args.continue_on_error:
                        raise
        except Exception as exc:
            LOG.exception("Failed event %s/%s", event.block_id, event.event_id)
            if not args.continue_on_error:
                raise
            summary.append(
                {
                    "status": "event_failed",
                    "block_id": event.block_id,
                    "event_id": event.event_id,
                    "error": repr(exc),
                }
            )
    failed = sum(str(row.get("status", "")).endswith("failed") for row in summary)
    LOG.info("Finished: %s records, %s failed", len(summary), failed)
    return 1 if failed else 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
        force=True,
    )
    try:
        return run(args)
    except KeyboardInterrupt:
        LOG.error("Interrupted by user")
        return 130
    except Exception:
        LOG.exception("Build failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
