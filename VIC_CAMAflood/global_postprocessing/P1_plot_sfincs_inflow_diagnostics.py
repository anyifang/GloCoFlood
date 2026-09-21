#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Plot per-GTC SFINCS inflow diagnostics and the global TC-season flow map.

For each GTC partition that has matched upstream CaMa-Flood cells, this script:
  1. plots the SFINCS domain and matched upstream CaMa river network;
  2. extracts historical and future CaMa-Flood discharge at matched cells;
  3. plots total-inflow daily seasonality;
  4. plots inlet-cell monthly seasonality;
  5. plots TC-season cumulative distribution curves.

Default inputs are set for the local project layout:
  - CaMa output: <external-data-root>/camaflood/output
  - CaMa map:    <external-data-root>/camaflood/map_v420/glb_15min
  - SFINCS GTC:  <work-root>/global_sfincs_partition_models_cama_15min
                 global_flood/global_sfincs_partition_models_cama_15min

The complete workflow is implemented in this file. It first rereads the CaMa
``outflw*.bin`` files and redraws every matched GTC, then reads the major-outlet
tables produced by ``plot_output_flow_change_global.py`` and draws the global
TC-season distribution map.

Run the complete workflow:
  python P1_plot_sfincs_inflow_diagnostics.py

Fast test with one GTC and one year step subset:
  python P1_plot_sfincs_inflow_diagnostics.py --ids GTC_0007 --max_files_per_case 1 --max_steps_per_file 10
"""

from __future__ import annotations

import argparse
import calendar
import hashlib
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm

try:
    import geopandas as gpd
except Exception:  # pragma: no cover - script can still plot CaMa geojson without geopandas
    gpd = None

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except Exception as exc:  # pragma: no cover
    raise RuntimeError("This combined script requires cartopy for the global map.") from exc


# =========================
# User settings
# =========================
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(os.environ.get("GLOCOFLOOD_CAMA_ROOT", SCRIPT_DIR / "external"))
INPUT_DIR = ROOT_DIR / "output"
MAP_DIR = ROOT_DIR / "glb_15min"
MODELS_DIR = Path(os.environ.get("GLOCOFLOOD_SFINCS_MODEL_ROOT", SCRIPT_DIR / "external" / "sfincs_models"))
OUT_DIR = INPUT_DIR / "gtc_cama_inflow_figures"
GLOBAL_RESULT_DIR = INPUT_DIR / "flow_change_and_outlet_seasonality"
GLOBAL_TABLE_DIR = GLOBAL_RESULT_DIR / "tables"
GLOBAL_FIGURE_DIR = GLOBAL_RESULT_DIR / "figures"

HISTORICAL_CASE = "ERA5"
FUTURE_CASE_REGEX = r"^cmip6_.*_ssp(126|245|370)_2061_2100$"
SIM_GLOB = "outflw*.bin"
HIST_START_YEAR = 1975
HIST_END_YEAR = 2014
FUTURE_START_YEAR = 2061
FUTURE_END_YEAR = 2100
BINARY_ENDIAN = "little"
ARRAY_ORDER = "tyx"
INVALID_THRESHOLD = -9990.0
MIN_VALID_FLOW = 0.0

# TC-season months used when no IBTrACS event-date table is supplied.
# Northern Hemisphere: Jun-Nov. Southern Hemisphere: Nov-Apr.
TC_MONTHS_NH = (6, 7, 8, 9, 10, 11)
TC_MONTHS_SH = (11, 12, 1, 2, 3, 4)
IBTRACS_PATH = Path(os.environ.get("GLOCOFLOOD_IBTRACS_NC", SCRIPT_DIR / "external" / "IBTrACS.ALL.v04r01.nc"))
FIGS9_BASIN_EVENTS_CSV = Path(
    os.environ.get(
        "GLOCOFLOOD_TC_EVENTS_CSV",
        SCRIPT_DIR / "external" / "05_basin_landfall_events.csv",
    )
)
TC_ANALYSIS_START_YEAR = 1975
TC_ANALYSIS_END_YEAR = 2014
TC_LANDFALL_WIND_THRESHOLD_MS = 33.0

SCENARIO_COLORS = {
    "historical": "#111111",
    "ssp126": "#2c7fb8",
    "ssp245": "#f28e2b",
    "ssp370": "#d62728",
    "ssp585": "#7b3294",
}
MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
MONTH_LENGTHS_365 = np.asarray([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=int)
GLOBAL_SCENARIO_ORDER = ("ssp126", "ssp245", "ssp370", "ssp585")


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "axes.unicode_minus": False,
        "font.size": 11,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


@dataclass(frozen=True)
class GridInfo:
    nx: int
    ny: int
    dx: float
    west: float
    east: float
    south: float
    north: float


@dataclass(frozen=True)
class CaseInfo:
    name: str
    directory: Path
    files: Tuple[Path, ...]
    scenario: str
    model: str


@dataclass(frozen=True)
class InletCell:
    gtc_id: str
    inlet_id: str
    cell_id: str
    row: int
    col: int
    lon: float
    lat: float
    uparea_km2: float
    q_reference_m3s: float
    sfincs_src_cols: str
    source_feature_type: str


@dataclass
class GTCModel:
    gtc_id: str
    directory: Path
    cama_geojson: Optional[Path]
    active_domain_geojson: Optional[Path]
    cells: List[InletCell]
    centroid_lon: float
    centroid_lat: float
    basin_id: str = ""
    basin_label: str = ""


@dataclass
class CaseExtract:
    case: CaseInfo
    cache_path: Path
    dates: pd.DatetimeIndex
    values: np.ndarray


@dataclass(frozen=True)
class TCBasinDef:
    basin_id: str
    basin_label: str
    lon_min: float
    lon_max: float
    lat_min: float
    lat_max: float


@dataclass(frozen=True)
class TCWeightInfo:
    basin_id: str
    basin_label: str
    source: str
    months: Tuple[int, ...]
    month_weights: Dict[int, float]
    event_count: int
    analysis_start_year: int
    analysis_end_year: int


@dataclass
class GTCSeries:
    case: CaseInfo
    daily_clim: np.ndarray
    monthly_clim: np.ndarray
    inlet_monthly: np.ndarray
    tc_samples: np.ndarray
    tc_sample_weights: np.ndarray
    mean_flow_m3s: float
    p50_tc_m3s: float
    p90_tc_m3s: float
    inlet_mean_flow_m3s: np.ndarray
    inlet_p50_tc_m3s: np.ndarray
    inlet_p90_tc_m3s: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input_dir", type=Path, default=INPUT_DIR)
    parser.add_argument("--map_dir", type=Path, default=MAP_DIR)
    parser.add_argument("--models_dir", type=Path, default=MODELS_DIR)
    parser.add_argument("--out_dir", type=Path, default=OUT_DIR)
    parser.add_argument("--historical_case", default=HISTORICAL_CASE)
    parser.add_argument("--future_case_regex", default=FUTURE_CASE_REGEX)
    parser.add_argument("--sim_glob", default=SIM_GLOB)
    parser.add_argument("--ids", nargs="*", default=None, help="GTC ids to process, e.g. GTC_0007 12 59.")
    parser.add_argument("--limit", type=int, default=None, help="Only process the first N matched GTC partitions.")
    parser.add_argument(
        "--tables_only",
        action="store_true",
        help="Write diagnostic tables without redrawing per-GTC or global figures.",
    )
    parser.add_argument("--hist_start_year", type=int, default=HIST_START_YEAR)
    parser.add_argument("--hist_end_year", type=int, default=HIST_END_YEAR)
    parser.add_argument("--future_start_year", type=int, default=FUTURE_START_YEAR)
    parser.add_argument("--future_end_year", type=int, default=FUTURE_END_YEAR)
    parser.add_argument("--max_files_per_case", type=int, default=None, help="Quick-test limit for each case.")
    parser.add_argument("--max_steps_per_file", type=int, default=None, help="Quick-test time-step limit for each file.")
    parser.add_argument("--binary_endian", choices=["little", "big"], default=BINARY_ENDIAN)
    parser.add_argument("--array_order", choices=["tyx", "txy"], default=ARRAY_ORDER)
    parser.add_argument("--read_mode", choices=["sparse", "sequential"], default="sparse")
    parser.add_argument("--chunk_time", type=int, default=16, help="Sequential read chunk length.")
    parser.add_argument("--invalid_threshold", type=float, default=INVALID_THRESHOLD)
    parser.add_argument("--min_valid_flow", type=float, default=MIN_VALID_FLOW)
    parser.add_argument("--tc_months_nh", nargs="+", type=int, default=list(TC_MONTHS_NH))
    parser.add_argument("--tc_months_sh", nargs="+", type=int, default=list(TC_MONTHS_SH))
    parser.add_argument("--tc_months", nargs="+", type=int, default=None, help="Force one TC-season month list for all GTCs.")
    parser.add_argument(
        "--tc_weight_source",
        choices=["figs9_ibtracs", "fixed_months"],
        default="figs9_ibtracs",
        help="Use Figure S9 observed landfall events plus IBTrACS timestamps, or fixed hemisphere months.",
    )
    parser.add_argument("--figs9_basin_events_csv", type=Path, default=FIGS9_BASIN_EVENTS_CSV)
    parser.add_argument("--ibtracs_path", type=Path, default=IBTRACS_PATH)
    parser.add_argument("--tc_analysis_start_year", type=int, default=TC_ANALYSIS_START_YEAR)
    parser.add_argument("--tc_analysis_end_year", type=int, default=TC_ANALYSIS_END_YEAR)
    parser.add_argument("--tc_landfall_wind_ms", type=float, default=TC_LANDFALL_WIND_THRESHOLD_MS)
    parser.add_argument("--dpi", type=int, default=350, help="Per-GTC figure DPI.")
    cache_mode = parser.add_mutually_exclusive_group()
    cache_mode.add_argument("--force", dest="force", action="store_true", help="Reread all bin files (default).")
    cache_mode.add_argument(
        "--reuse_gtc_cache",
        dest="force",
        action="store_false",
        help="Reuse complete GTC extraction caches instead of rereading all bin files.",
    )
    parser.set_defaults(force=True)
    parser.add_argument("--global_table_dir", type=Path, default=GLOBAL_TABLE_DIR)
    parser.add_argument("--global_figure_dir", type=Path, default=GLOBAL_FIGURE_DIR)
    parser.add_argument("--global_monthly_table", default="outlet_monthly_climatology_all_cases.csv")
    parser.add_argument("--global_outlet_table", default="major_outlets.csv")
    parser.add_argument("--global_months_nh", nargs="+", type=int, default=list(TC_MONTHS_NH))
    parser.add_argument("--global_months_sh", nargs="+", type=int, default=list(TC_MONTHS_SH))
    parser.add_argument("--global_change_vlim", type=float, default=None)
    parser.add_argument("--global_size_min", type=float, default=45.0)
    parser.add_argument("--global_size_max", type=float, default=430.0)
    parser.add_argument("--global_dpi", type=int, default=450)
    parser.add_argument("--label_outlets", action="store_true")
    parser.add_argument("--global_output_name", default="global_major_outlet_tc_season_flow_change_map")
    return parser.parse_args()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_name(text: str) -> str:
    out = re.sub(r"[^\w\-.]+", "_", str(text).strip())
    return out.strip("._") or "case"


def normalize_gtc_id(value: str) -> str:
    text = str(value).strip()
    match = re.search(r"(\d+)", text)
    if match:
        return f"GTC_{int(match.group(1)):04d}"
    return text.upper()


def parse_cama_params(map_dir: Path) -> GridInfo:
    values: List[float] = []
    with (map_dir / "params.txt").open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            raw = raw.split("!!", 1)[0].strip()
            if not raw:
                continue
            match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", raw)
            if match:
                values.append(float(match.group(0)))
    if len(values) < 8:
        raise ValueError(f"Could not parse CaMa params from {map_dir / 'params.txt'}")
    return GridInfo(
        nx=int(values[0]),
        ny=int(values[1]),
        dx=float(values[3]),
        west=float(values[4]),
        east=float(values[5]),
        south=float(values[6]),
        north=float(values[7]),
    )


def row_col_to_lon_lat(row: int, col: int, grid: GridInfo) -> Tuple[float, float]:
    lon = grid.west + (col + 0.5) * grid.dx
    dy = abs((grid.north - grid.south) / grid.ny)
    lat = grid.north - (row + 0.5) * dy
    return float(lon), float(lat)


def infer_year(path: Path) -> Optional[int]:
    match = re.search(r"(19|20|21)\d{2}", path.stem)
    return int(match.group(0)) if match else None


def infer_dt_hours(nt: int) -> float:
    if nt in {365, 366}:
        return 24.0
    if nt in {365 * 4, 366 * 4}:
        return 6.0
    if nt in {365 * 8, 366 * 8}:
        return 3.0
    return 24.0


def build_dates(path: Path, nt: int) -> pd.DatetimeIndex:
    year = infer_year(path)
    if year is None:
        raise ValueError(f"Cannot infer year from {path.name}")
    return pd.date_range(
        start=pd.Timestamp(year=year, month=1, day=1),
        periods=nt,
        freq=pd.Timedelta(hours=infer_dt_hours(nt)),
    )


def endian_dtype(binary_endian: str) -> np.dtype:
    return np.dtype("<f4" if binary_endian == "little" else ">f4")


def infer_time_count(path: Path, grid: GridInfo, dtype: np.dtype) -> int:
    cells = grid.nx * grid.ny
    record_bytes = cells * dtype.itemsize
    file_size = path.stat().st_size
    if file_size % record_bytes != 0:
        raise ValueError(f"{path} size {file_size} is not divisible by one {grid.ny}x{grid.nx} record.")
    nt = file_size // record_bytes
    if nt <= 0:
        raise ValueError(f"No records found in {path}")
    return int(nt)


def reshape_time_major(raw: np.ndarray, nt: int, grid: GridInfo, array_order: str) -> np.ndarray:
    if array_order == "tyx":
        return raw.reshape(nt, grid.ny, grid.nx)
    if array_order == "txy":
        return np.transpose(raw.reshape(nt, grid.nx, grid.ny), (0, 2, 1))
    raise ValueError(array_order)


def filter_files(files: Sequence[Path], start_year: Optional[int], end_year: Optional[int], max_files: Optional[int]) -> Tuple[Path, ...]:
    selected: List[Path] = []
    for path in sorted(files, key=lambda p: (infer_year(p) or -1, p.name)):
        year = infer_year(path)
        if year is None:
            continue
        if start_year is not None and year < start_year:
            continue
        if end_year is not None and year > end_year:
            continue
        selected.append(path)
    if max_files is not None:
        selected = selected[: max(0, int(max_files))]
    return tuple(selected)


def parse_future_case_name(name: str) -> Tuple[str, str]:
    match = re.match(r"cmip6_(.+)_(ssp\d+)_(\d{4})_(\d{4})$", name)
    if match:
        return match.group(2), match.group(1)
    scenario_match = re.search(r"(ssp\d+)", name)
    scenario = scenario_match.group(1) if scenario_match else "future"
    model = re.sub(r"^cmip6_", "", name)
    model = re.sub(r"_ssp\d+.*$", "", model)
    return scenario, model or name


def discover_cases(args: argparse.Namespace) -> List[CaseInfo]:
    cases: List[CaseInfo] = []
    hist_dir = args.input_dir / args.historical_case
    hist_files = filter_files(
        tuple(hist_dir.glob(args.sim_glob)),
        args.hist_start_year,
        args.hist_end_year,
        args.max_files_per_case,
    )
    if not hist_files:
        raise FileNotFoundError(f"No historical files found: {hist_dir / args.sim_glob}")
    cases.append(
        CaseInfo(
            name=args.historical_case,
            directory=hist_dir,
            files=hist_files,
            scenario="historical",
            model=args.historical_case,
        )
    )

    pattern = re.compile(args.future_case_regex)
    for directory in sorted(p for p in args.input_dir.iterdir() if p.is_dir()):
        if directory.name == args.historical_case:
            continue
        if not pattern.search(directory.name):
            continue
        files = filter_files(
            tuple(directory.glob(args.sim_glob)),
            args.future_start_year,
            args.future_end_year,
            args.max_files_per_case,
        )
        if not files:
            print(f"[WARN] skip {directory.name}: no {args.sim_glob} files after year filter")
            continue
        scenario, model = parse_future_case_name(directory.name)
        cases.append(CaseInfo(directory.name, directory, files, scenario, model))
    if len(cases) == 1:
        print("[WARN] No future cases matched --future_case_regex; only historical will be plotted.")
    return cases


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def find_gtc_geojson(directory: Path) -> Optional[Path]:
    candidates = [
        directory / "cama_upstream_rivers.geojson",
        directory / "gis" / "cama_upstream_rivers.geojson",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def find_active_domain_geojson(directory: Path) -> Optional[Path]:
    candidates = [
        directory / "active_domain.geojson",
        directory / "region_grid.geojson",
        directory / "region.geojson",
        directory / "gis" / "active_domain.geojson",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def first_number(value, default: float = math.nan) -> float:
    if value is None:
        return default
    if isinstance(value, (int, float, np.number)):
        out = float(value)
        return out if math.isfinite(out) else default
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(value))
    if not match:
        return default
    out = float(match.group(0))
    return out if math.isfinite(out) else default


def normalize_row_col(row_raw, col_raw, lon_raw, lat_raw, grid: GridInfo) -> Optional[Tuple[int, int]]:
    if row_raw is None or col_raw is None:
        return None
    try:
        row0 = int(round(float(row_raw)))
        col0 = int(round(float(col_raw)))
    except Exception:
        return None

    candidates: List[Tuple[int, int, float]] = []
    lon_ref = first_number(lon_raw)
    lat_ref = first_number(lat_raw)
    for row, col in ((row0, col0), (row0 - 1, col0 - 1)):
        if 0 <= row < grid.ny and 0 <= col < grid.nx:
            if math.isfinite(lon_ref) and math.isfinite(lat_ref):
                lon, lat = row_col_to_lon_lat(row, col, grid)
                score = abs(lon - lon_ref) + abs(lat - lat_ref)
            else:
                score = 0.0 if (row, col) == (row0, col0) else 1.0
            candidates.append((row, col, score))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[2])
    return int(candidates[0][0]), int(candidates[0][1])


def feature_type(feature: dict) -> str:
    props = feature.get("properties") or {}
    return str(props.get("feature_type") or props.get("type") or "").strip()


def read_gtc_cells(gtc_id: str, cama_geojson: Path, grid: GridInfo) -> List[InletCell]:
    data = load_json(cama_geojson)
    features = data.get("features") or []
    preferred = [f for f in features if feature_type(f) == "CaMa_Inlet_Point"]
    if not preferred:
        preferred = [f for f in features if "cama_row" in (f.get("properties") or {}) and "cama_col" in (f.get("properties") or {})]

    cells: List[InletCell] = []
    seen: set[Tuple[int, int]] = set()
    for seq, feature in enumerate(preferred, start=1):
        props = feature.get("properties") or {}
        normalized = normalize_row_col(
            props.get("cama_row"),
            props.get("cama_col"),
            props.get("cama_lon"),
            props.get("cama_lat"),
            grid,
        )
        if normalized is None:
            continue
        row, col = normalized
        if (row, col) in seen:
            continue
        seen.add((row, col))
        lon, lat = row_col_to_lon_lat(row, col, grid)
        inlet_id = str(props.get("inlet_id") or props.get("boundary_id") or seq)
        cells.append(
            InletCell(
                gtc_id=gtc_id,
                inlet_id=inlet_id,
                cell_id=f"r{row:04d}_c{col:04d}",
                row=row,
                col=col,
                lon=first_number(props.get("cama_lon"), lon),
                lat=first_number(props.get("cama_lat"), lat),
                uparea_km2=first_number(props.get("uparea_km2")),
                q_reference_m3s=first_number(props.get("q_total_m3s")),
                sfincs_src_cols=str(props.get("sfincs_src_cols") or ""),
                source_feature_type=feature_type(feature),
            )
        )
    return cells


def geometry_bounds(geometry: dict) -> Optional[Tuple[float, float, float, float]]:
    xs: List[float] = []
    ys: List[float] = []

    def visit_coords(coords):
        if not coords:
            return
        first = coords[0]
        if isinstance(first, (int, float)):
            if len(coords) >= 2:
                xs.append(float(coords[0]))
                ys.append(float(coords[1]))
            return
        for item in coords:
            visit_coords(item)

    visit_coords(geometry.get("coordinates"))
    if not xs:
        return None
    return min(xs), min(ys), max(xs), max(ys)


def active_domain_metadata(path: Optional[Path]) -> Tuple[float, float, str, str]:
    if path is None or not path.exists():
        return math.nan, math.nan, "", ""
    try:
        data = load_json(path)
        feature = (data.get("features") or [{}])[0]
        props = feature.get("properties") or {}
        lon = first_number(props.get("lon"))
        lat = first_number(props.get("lat"))
        basin_id = str(props.get("basin_id") or props.get("quicklook_basin_id") or "")
        basin_label = str(props.get("basin_label") or props.get("basin_label_domain") or "")
        if not (math.isfinite(lon) and math.isfinite(lat)):
            bounds = geometry_bounds(feature.get("geometry") or {})
            if bounds:
                lon = 0.5 * (bounds[0] + bounds[2])
                lat = 0.5 * (bounds[1] + bounds[3])
        return lon, lat, basin_id, basin_label
    except Exception:
        return math.nan, math.nan, "", ""


def discover_gtc_models(args: argparse.Namespace, grid: GridInfo) -> Tuple[List[GTCModel], pd.DataFrame]:
    requested = {normalize_gtc_id(x) for x in args.ids} if args.ids else None
    status_rows: List[Dict[str, object]] = []
    models: List[GTCModel] = []

    for directory in sorted(p for p in args.models_dir.iterdir() if p.is_dir() and re.search(r"GTC_\d{4}", p.name)):
        gtc_id = normalize_gtc_id(directory.name)
        if requested and gtc_id not in requested:
            continue
        cama_geojson = find_gtc_geojson(directory)
        active_domain = find_active_domain_geojson(directory)
        cells = read_gtc_cells(gtc_id, cama_geojson, grid) if cama_geojson else []
        lon, lat, basin_id, basin_label = active_domain_metadata(active_domain)
        if not math.isfinite(lon) and cells:
            lon = float(np.nanmean([c.lon for c in cells]))
            lat = float(np.nanmean([c.lat for c in cells]))

        status = "matched" if cells else "skipped_no_cama_cells"
        status_rows.append(
            {
                "gtc_id": gtc_id,
                "directory": str(directory),
                "has_cama_upstream_geojson": cama_geojson is not None,
                "has_active_domain_geojson": active_domain is not None,
                "n_cama_cells": len(cells),
                "status": status,
                "centroid_lon": lon,
                "centroid_lat": lat,
                "basin_id": basin_id,
                "basin_label": basin_label,
            }
        )
        models.append(
            GTCModel(
                gtc_id=gtc_id,
                directory=directory,
                cama_geojson=cama_geojson,
                active_domain_geojson=active_domain,
                cells=cells,
                centroid_lon=lon,
                centroid_lat=lat,
                basin_id=basin_id,
                basin_label=basin_label,
            )
        )

    matched = [m for m in models if m.cells]
    if args.limit is not None:
        keep_ids = {m.gtc_id for m in matched[: max(0, int(args.limit))]}
        models = [m for m in models if (m.gtc_id in keep_ids or not m.cells)]
        matched = [m for m in models if m.cells]
    if not matched:
        raise RuntimeError("No GTC partition has matched CaMa cells after filtering.")
    return models, pd.DataFrame(status_rows)


def selected_models(models: Sequence[GTCModel]) -> List[GTCModel]:
    return [m for m in models if m.cells]


def unique_cells_from_models(models: Sequence[GTCModel]) -> List[Tuple[int, int]]:
    cells = sorted({(cell.row, cell.col) for model in models for cell in model.cells})
    return [(int(r), int(c)) for r, c in cells]


def write_inlet_table(models: Sequence[GTCModel], path: Path) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for model in models:
        for order, cell in enumerate(model.cells, start=1):
            rows.append(
                {
                    "gtc_id": model.gtc_id,
                    "inlet_order": order,
                    "inlet_id": cell.inlet_id,
                    "cell_id": cell.cell_id,
                    "cama_row": cell.row,
                    "cama_col": cell.col,
                    "cama_lon": cell.lon,
                    "cama_lat": cell.lat,
                    "uparea_km2": cell.uparea_km2,
                    "q_reference_m3s": cell.q_reference_m3s,
                    "sfincs_src_cols": cell.sfincs_src_cols,
                    "source_feature_type": cell.source_feature_type,
                    "gtc_directory": str(model.directory),
                    "cama_upstream_geojson": str(model.cama_geojson or ""),
                    "active_domain_geojson": str(model.active_domain_geojson or ""),
                    "basin_id": model.basin_id,
                    "basin_label": model.basin_label,
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return df


def extraction_hash(unique_cells: Sequence[Tuple[int, int]], args: argparse.Namespace) -> str:
    payload = {
        "cells": list(unique_cells),
        "array_order": args.array_order,
        "binary_endian": args.binary_endian,
        "hist_start_year": args.hist_start_year,
        "hist_end_year": args.hist_end_year,
        "future_start_year": args.future_start_year,
        "future_end_year": args.future_end_year,
        "max_files_per_case": args.max_files_per_case,
        "max_steps_per_file": args.max_steps_per_file,
        "invalid_threshold": args.invalid_threshold,
        "min_valid_flow": args.min_valid_flow,
    }
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:12]


def cache_path_for_case(cache_dir: Path, case: CaseInfo, key: str) -> Path:
    return cache_dir / f"cell_daily_{safe_name(case.name)}_{key}.npz"


def clean_values(values: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    out = values.astype(np.float32, copy=False)
    valid = np.isfinite(out)
    valid &= out > args.invalid_threshold
    valid &= out >= args.min_valid_flow
    return np.where(valid, out, np.nan).astype(np.float32, copy=False)


def collapse_to_daily(dates: pd.DatetimeIndex, values: np.ndarray) -> Tuple[pd.DatetimeIndex, np.ndarray]:
    day_index = pd.DatetimeIndex(dates).normalize()
    if len(day_index) == len(pd.unique(day_index)):
        return day_index, values
    df = pd.DataFrame(values)
    df["date"] = day_index
    grouped = df.groupby("date", sort=True).mean(numeric_only=True)
    return pd.DatetimeIndex(grouped.index), grouped.to_numpy(dtype=np.float32)


def read_file_cells_sparse(
    path: Path,
    rows: np.ndarray,
    cols: np.ndarray,
    grid: GridInfo,
    dtype: np.dtype,
    args: argparse.Namespace,
) -> Tuple[pd.DatetimeIndex, np.ndarray]:
    nt_total = infer_time_count(path, grid, dtype)
    nt = min(nt_total, int(args.max_steps_per_file)) if args.max_steps_per_file else nt_total
    if args.array_order == "tyx":
        flat_idx = rows * grid.nx + cols
        shape = (nt_total, grid.ny * grid.nx)
    else:
        flat_idx = cols * grid.ny + rows
        shape = (nt_total, grid.nx * grid.ny)
    mm = np.memmap(path, dtype=dtype, mode="r", shape=shape)
    values = np.asarray(mm[:nt, flat_idx], dtype=np.float32)
    del mm
    return build_dates(path, nt), clean_values(values, args)


def read_file_cells_sequential(
    path: Path,
    rows: np.ndarray,
    cols: np.ndarray,
    grid: GridInfo,
    dtype: np.dtype,
    args: argparse.Namespace,
) -> Tuple[pd.DatetimeIndex, np.ndarray]:
    nt_total = infer_time_count(path, grid, dtype)
    nt = min(nt_total, int(args.max_steps_per_file)) if args.max_steps_per_file else nt_total
    cells = grid.nx * grid.ny
    chunk_time = max(1, int(args.chunk_time))
    chunks: List[np.ndarray] = []
    with path.open("rb") as handle:
        for t0 in range(0, nt, chunk_time):
            t1 = min(t0 + chunk_time, nt)
            n_time = t1 - t0
            raw = np.fromfile(handle, dtype=dtype, count=n_time * cells)
            if raw.size != n_time * cells:
                raise IOError(f"Unexpected end of file while reading {path}")
            block = reshape_time_major(raw, n_time, grid, args.array_order)
            chunks.append(clean_values(block[:, rows, cols], args))
    return build_dates(path, nt), np.vstack(chunks)


def load_npz_extract(path: Path) -> Tuple[pd.DatetimeIndex, np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as npz:
        dates = pd.to_datetime(npz["dates"].astype("datetime64[D]"))
        values = np.asarray(npz["values"], dtype=np.float32)
        rows = np.asarray(npz["rows"], dtype=np.int64)
        cols = np.asarray(npz["cols"], dtype=np.int64)
    return pd.DatetimeIndex(dates), values, rows, cols


def cache_is_complete(path: Path, n_cells: int) -> bool:
    if not path.exists() or path.stat().st_size <= 0:
        return False
    try:
        dates, values, rows, cols = load_npz_extract(path)
    except Exception:
        return False
    return len(dates) == values.shape[0] and values.shape[1] == n_cells and len(rows) == n_cells and len(cols) == n_cells


def extract_case_cells(
    case: CaseInfo,
    unique_cells: Sequence[Tuple[int, int]],
    grid: GridInfo,
    args: argparse.Namespace,
    cache_dir: Path,
    key: str,
) -> CaseExtract:
    cache_path = cache_path_for_case(cache_dir, case, key)
    if not args.force and cache_is_complete(cache_path, len(unique_cells)):
        dates, values, _, _ = load_npz_extract(cache_path)
        print(f"[SKIP] {case.name}: existing complete cache -> {cache_path.name}")
        return CaseExtract(case, cache_path, dates, values)

    rows = np.asarray([r for r, _ in unique_cells], dtype=np.int64)
    cols = np.asarray([c for _, c in unique_cells], dtype=np.int64)
    dtype = endian_dtype(args.binary_endian)
    dates_all: List[pd.DatetimeIndex] = []
    values_all: List[np.ndarray] = []
    print(f"[READ] {case.name}: extracting {len(unique_cells)} CaMa cells from {len(case.files)} file(s)")
    for file_idx, path in enumerate(case.files, start=1):
        print(f"  [INFO] {case.name}: {file_idx}/{len(case.files)} {path.name}")
        if args.read_mode == "sparse":
            dates, values = read_file_cells_sparse(path, rows, cols, grid, dtype, args)
        else:
            dates, values = read_file_cells_sequential(path, rows, cols, grid, dtype, args)
        dates, values = collapse_to_daily(dates, values)
        dates_all.append(dates)
        values_all.append(values)

    dates_out = pd.DatetimeIndex(np.concatenate([d.values.astype("datetime64[D]") for d in dates_all]))
    values_out = np.vstack(values_all).astype(np.float32, copy=False)

    tmp_path = cache_path.with_suffix(".tmp.npz")
    with tmp_path.open("wb") as handle:
        np.savez(handle, dates=dates_out.values.astype("datetime64[D]"), values=values_out, rows=rows, cols=cols)
    tmp_path.replace(cache_path)
    return CaseExtract(case, cache_path, dates_out, values_out)


def clim_day_index(date: pd.Timestamp) -> int:
    if date.month == 2 and date.day == 29:
        return -1
    doy = int(date.dayofyear)
    if calendar.isleap(date.year) and doy > 60:
        doy -= 1
    return doy - 1


def daily_climatology(dates: pd.DatetimeIndex, values: np.ndarray) -> np.ndarray:
    total = np.zeros(365, dtype=np.float64)
    count = np.zeros(365, dtype=np.int64)
    for date, value in zip(pd.DatetimeIndex(dates), values):
        idx = clim_day_index(pd.Timestamp(date))
        if idx < 0 or not np.isfinite(value):
            continue
        total[idx] += float(value)
        count[idx] += 1
    out = np.full(365, np.nan, dtype=np.float64)
    good = count > 0
    out[good] = total[good] / count[good]
    return out


def monthly_climatology(dates: pd.DatetimeIndex, values: np.ndarray) -> np.ndarray:
    total = np.zeros(12, dtype=np.float64)
    count = np.zeros(12, dtype=np.int64)
    for date, value in zip(pd.DatetimeIndex(dates), values):
        if not np.isfinite(value):
            continue
        idx = int(pd.Timestamp(date).month) - 1
        total[idx] += float(value)
        count[idx] += 1
    out = np.full(12, np.nan, dtype=np.float64)
    good = count > 0
    out[good] = total[good] / count[good]
    return out


def normalize_lon_360(lon: float) -> float:
    value = float(lon)
    return value % 360.0


def tc_basin_definitions() -> Dict[str, TCBasinDef]:
    return {
        "BASIN_NATL": TCBasinDef("BASIN_NATL", "North Atlantic", -100.0, -25.0, 0.0, 52.0),
        "BASIN_WNP": TCBasinDef("BASIN_WNP", "Western North Pacific", 100.0, 180.0, 0.0, 52.0),
        "BASIN_NIO": TCBasinDef("BASIN_NIO", "North Indian Ocean", 35.0, 105.0, 0.0, 32.0),
        "BASIN_SIO": TCBasinDef("BASIN_SIO", "South Indian Ocean", 30.0, 105.0, -45.0, 6.0),
        "BASIN_AUSSP": TCBasinDef("BASIN_AUSSP", "South Pacific", 105.0, 240.0, -45.0, 8.0),
        "BASIN_ENP": TCBasinDef("BASIN_ENP", "Eastern North Pacific", -180.0, -75.0, 0.0, 52.0),
    }


def lon_in_basin(lon: float, basin: TCBasinDef) -> bool:
    lon360 = normalize_lon_360(lon)
    lon_min = normalize_lon_360(basin.lon_min)
    lon_max = normalize_lon_360(basin.lon_max)
    if lon_min <= lon_max:
        return lon_min <= lon360 <= lon_max
    return lon360 >= lon_min or lon360 <= lon_max


def basin_for_model(model: GTCModel, basin_defs: Dict[str, TCBasinDef]) -> Optional[TCBasinDef]:
    basin_id = str(model.basin_id or "").strip()
    if basin_id in basin_defs:
        return basin_defs[basin_id]
    lon = model.centroid_lon
    lat = model.centroid_lat
    if not (math.isfinite(lon) and math.isfinite(lat)):
        return None
    candidates = [b for b in basin_defs.values() if b.lat_min <= lat <= b.lat_max and lon_in_basin(lon, b)]
    if candidates:
        return candidates[0]
    return None


def tc_months_for_gtc(model: GTCModel, args: argparse.Namespace) -> Tuple[int, ...]:
    if args.tc_months:
        return tuple(sorted({int(m) for m in args.tc_months if 1 <= int(m) <= 12}))
    months = args.tc_months_nh if model.centroid_lat >= 0 else args.tc_months_sh
    return tuple(sorted({int(m) for m in months if 1 <= int(m) <= 12}))


def fixed_month_tc_weight_info(model: GTCModel, args: argparse.Namespace, source: str = "fixed_months") -> TCWeightInfo:
    months = tc_months_for_gtc(model, args)
    if not months:
        months = tuple(range(1, 13))
    weight = 1.0 / float(len(months))
    basin_id = str(model.basin_id or "")
    basin_label = str(model.basin_label or "")
    return TCWeightInfo(
        basin_id=basin_id,
        basin_label=basin_label,
        source=source,
        months=months,
        month_weights={int(m): weight for m in months},
        event_count=0,
        analysis_start_year=int(args.tc_analysis_start_year),
        analysis_end_year=int(args.tc_analysis_end_year),
    )


def decode_char_array(row: np.ndarray) -> str:
    items: List[str] = []
    for item in row:
        if isinstance(item, bytes):
            items.append(item.decode("ascii", errors="ignore"))
        else:
            items.append(str(item))
    return "".join(items).strip()


def load_figs9_observed_event_months(args: argparse.Namespace) -> Tuple[Dict[str, Dict[int, int]], pd.DataFrame]:
    event_path = Path(args.figs9_basin_events_csv)
    ibtracs_path = Path(args.ibtracs_path)
    if not event_path.exists():
        raise FileNotFoundError(f"Figure S9 basin event CSV not found: {event_path}")
    if not ibtracs_path.exists():
        raise FileNotFoundError(f"IBTrACS NetCDF not found: {ibtracs_path}")

    events = pd.read_csv(event_path)
    required = {"dataset_role", "basin_id", "track_index", "year", "landfall_vmax_ms", "max_time_index"}
    missing = required.difference(events.columns)
    if missing:
        raise ValueError(f"Figure S9 event CSV missing required columns: {sorted(missing)}")
    events = events.loc[events["dataset_role"].astype(str).str.lower().eq("observed")].copy()
    events = events.loc[
        (events["year"].astype(float) >= float(args.tc_analysis_start_year))
        & (events["year"].astype(float) <= float(args.tc_analysis_end_year))
        & (events["landfall_vmax_ms"].astype(float) > float(args.tc_landfall_wind_ms))
    ].copy()
    if events.empty:
        raise ValueError("No observed Figure S9 landfall events remain after year/intensity filtering.")

    try:
        from netCDF4 import Dataset, num2date
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("Missing python package: netCDF4") from exc

    with Dataset(ibtracs_path) as ds:
        time_var = ds.variables["time"]
        max_storm = int(ds.dimensions["storm"].size)
        months: List[float] = []
        for _, row in events.iterrows():
            storm_idx = int(row["track_index"])
            time_idx = int(row["max_time_index"])
            if storm_idx < 0 or storm_idx >= max_storm or time_idx < 0 or time_idx >= time_var.shape[1]:
                months.append(math.nan)
                continue
            value = time_var[storm_idx, time_idx]
            if np.ma.is_masked(value) or not np.isfinite(float(value)):
                months.append(math.nan)
                continue
            dt = num2date(
                float(value),
                units=time_var.units,
                calendar=getattr(time_var, "calendar", "standard"),
                only_use_cftime_datetimes=False,
                only_use_python_datetimes=True,
            )
            months.append(float(dt.month))

    events["tc_month"] = months
    events = events.loc[events["tc_month"].notna()].copy()
    events["tc_month"] = events["tc_month"].astype(int)
    month_counts: Dict[str, Dict[int, int]] = {}
    for basin_id, group in events.groupby("basin_id"):
        counts = group["tc_month"].value_counts().sort_index()
        month_counts[str(basin_id)] = {int(month): int(count) for month, count in counts.items()}
    return month_counts, events


def build_tc_weight_infos(
    models: Sequence[GTCModel],
    args: argparse.Namespace,
    table_dir: Path,
) -> Tuple[Dict[str, TCWeightInfo], pd.DataFrame, pd.DataFrame]:
    basin_defs = tc_basin_definitions()
    month_counts: Dict[str, Dict[int, int]] = {}
    event_df = pd.DataFrame()
    if args.tc_weight_source == "figs9_ibtracs":
        try:
            month_counts, event_df = load_figs9_observed_event_months(args)
            if not event_df.empty:
                event_df.to_csv(table_dir / "tc_observed_landfall_events_with_month.csv", index=False, encoding="utf-8-sig")
            print(f"[INFO] Loaded observed TC landfall month weights from {args.figs9_basin_events_csv}")
        except Exception as exc:
            print(f"[WARN] Could not load Figure S9/IBTrACS TC weights, falling back to fixed months: {exc}")

    infos: Dict[str, TCWeightInfo] = {}
    rows: List[Dict[str, object]] = []
    for model in models:
        basin_def = basin_for_model(model, basin_defs)
        counts = month_counts.get(basin_def.basin_id, {}) if basin_def is not None else {}
        if counts:
            total = float(sum(counts.values()))
            weights = {int(month): float(count) / total for month, count in sorted(counts.items())}
            info = TCWeightInfo(
                basin_id=basin_def.basin_id,
                basin_label=basin_def.basin_label,
                source="figs9_observed_landfall_ibtracs_month",
                months=tuple(sorted(weights)),
                month_weights=weights,
                event_count=int(total),
                analysis_start_year=int(args.tc_analysis_start_year),
                analysis_end_year=int(args.tc_analysis_end_year),
            )
        else:
            info = fixed_month_tc_weight_info(model, args, source="fixed_months_fallback")
        infos[model.gtc_id] = info
        for month in range(1, 13):
            rows.append(
                {
                    "gtc_id": model.gtc_id,
                    "tc_basin_id": info.basin_id,
                    "tc_basin_label": info.basin_label,
                    "tc_weight_source": info.source,
                    "tc_event_count": info.event_count,
                    "month": month,
                    "month_weight": float(info.month_weights.get(month, 0.0)),
                }
            )
    weights_df = pd.DataFrame(rows)
    weights_df.to_csv(table_dir / "gtc_tc_month_weights.csv", index=False, encoding="utf-8-sig")
    return infos, weights_df, event_df


def weighted_quantile(values: np.ndarray, weights: np.ndarray, percentile: float) -> float:
    finite = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not finite.any():
        return math.nan
    v = values[finite].astype(np.float64, copy=False)
    w = weights[finite].astype(np.float64, copy=False)
    order = np.argsort(v)
    v = v[order]
    w = w[order]
    cumulative = np.cumsum(w)
    total = cumulative[-1]
    if total <= 0:
        return math.nan
    target = float(percentile) / 100.0 * total
    return float(np.interp(target, cumulative, v))


def tc_weighted_samples(
    dates: pd.DatetimeIndex,
    values: np.ndarray,
    month_weights: Dict[int, float],
) -> Tuple[np.ndarray, np.ndarray]:
    date_index = pd.DatetimeIndex(dates)
    months = date_index.month.to_numpy(dtype=int)
    samples: List[np.ndarray] = []
    weights: List[np.ndarray] = []
    for month, month_weight in sorted(month_weights.items()):
        if month_weight <= 0:
            continue
        mask = (months == int(month)) & np.isfinite(values)
        if not mask.any():
            continue
        vals = values[mask].astype(np.float64, copy=False)
        sample_weight = np.full(vals.shape, float(month_weight) / float(vals.size), dtype=np.float64)
        samples.append(vals)
        weights.append(sample_weight)
    if not samples:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
    return np.concatenate(samples), np.concatenate(weights)


def gtc_series_from_case(
    model: GTCModel,
    case_extract: CaseExtract,
    cell_to_index: Dict[Tuple[int, int], int],
    tc_weight_info: TCWeightInfo,
) -> GTCSeries:
    col_indices = [cell_to_index[(cell.row, cell.col)] for cell in model.cells]
    values = case_extract.values[:, col_indices]
    any_valid = np.isfinite(values).any(axis=1)
    total = np.nansum(values, axis=1).astype(np.float64)
    total[~any_valid] = np.nan

    daily = daily_climatology(case_extract.dates, total)
    monthly = monthly_climatology(case_extract.dates, total)
    inlet_monthly = np.vstack([monthly_climatology(case_extract.dates, values[:, i]) for i in range(values.shape[1])])
    tc_samples, tc_sample_weights = tc_weighted_samples(case_extract.dates, total, tc_weight_info.month_weights)
    mean_flow = float(np.nanmean(total)) if np.isfinite(total).any() else math.nan
    p50 = weighted_quantile(tc_samples, tc_sample_weights, 50) if tc_samples.size else math.nan
    p90 = weighted_quantile(tc_samples, tc_sample_weights, 90) if tc_samples.size else math.nan
    inlet_mean = np.full(values.shape[1], np.nan, dtype=np.float64)
    inlet_p50 = np.full(values.shape[1], np.nan, dtype=np.float64)
    inlet_p90 = np.full(values.shape[1], np.nan, dtype=np.float64)
    for inlet_idx in range(values.shape[1]):
        inlet_values = values[:, inlet_idx]
        if np.isfinite(inlet_values).any():
            inlet_mean[inlet_idx] = float(np.nanmean(inlet_values))
        inlet_samples, inlet_weights = tc_weighted_samples(
            case_extract.dates,
            inlet_values,
            tc_weight_info.month_weights,
        )
        if inlet_samples.size:
            inlet_p50[inlet_idx] = weighted_quantile(inlet_samples, inlet_weights, 50)
            inlet_p90[inlet_idx] = weighted_quantile(inlet_samples, inlet_weights, 90)
    return GTCSeries(
        case=case_extract.case,
        daily_clim=daily,
        monthly_clim=monthly,
        inlet_monthly=inlet_monthly,
        tc_samples=tc_samples.astype(np.float64, copy=False),
        tc_sample_weights=tc_sample_weights.astype(np.float64, copy=False),
        mean_flow_m3s=mean_flow,
        p50_tc_m3s=p50,
        p90_tc_m3s=p90,
        inlet_mean_flow_m3s=inlet_mean,
        inlet_p50_tc_m3s=inlet_p50,
        inlet_p90_tc_m3s=inlet_p90,
    )


def group_series(series: Sequence[GTCSeries]) -> Dict[str, List[GTCSeries]]:
    grouped: Dict[str, List[GTCSeries]] = {}
    for item in series:
        grouped.setdefault(item.case.scenario, []).append(item)
    return grouped


def nanmean_stack(items: Sequence[np.ndarray]) -> np.ndarray:
    if not items:
        return np.array([])
    stack = np.vstack(items).astype(np.float64, copy=False)
    valid = np.isfinite(stack)
    count = valid.sum(axis=0)
    total = np.where(valid, stack, 0.0).sum(axis=0)
    out = np.full(stack.shape[1], np.nan, dtype=np.float64)
    good = count > 0
    out[good] = total[good] / count[good]
    return out


def nanpercentile_stack(items: Sequence[np.ndarray], percentile: float) -> np.ndarray:
    if not items:
        return np.array([])
    stack = np.vstack(items).astype(np.float64, copy=False)
    out = np.full(stack.shape[1], np.nan, dtype=np.float64)
    for col in range(stack.shape[1]):
        vals = stack[:, col]
        vals = vals[np.isfinite(vals)]
        if vals.size:
            out[col] = np.percentile(vals, percentile)
    return out


def empirical_cdf(samples: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    samples = np.sort(samples[np.isfinite(samples)])
    if samples.size == 0:
        return np.full_like(x_grid, np.nan, dtype=np.float64)
    return np.searchsorted(samples, x_grid, side="right") / samples.size * 100.0


def weighted_cdf(samples: np.ndarray, weights: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    finite = np.isfinite(samples) & np.isfinite(weights) & (weights > 0)
    if not finite.any():
        return np.full_like(x_grid, np.nan, dtype=np.float64)
    values = samples[finite].astype(np.float64, copy=False)
    sample_weights = weights[finite].astype(np.float64, copy=False)
    order = np.argsort(values)
    values = values[order]
    sample_weights = sample_weights[order]
    cumulative = np.cumsum(sample_weights)
    total = cumulative[-1]
    if total <= 0:
        return np.full_like(x_grid, np.nan, dtype=np.float64)
    indices = np.searchsorted(values, x_grid, side="right") - 1
    out = np.zeros_like(x_grid, dtype=np.float64)
    valid = indices >= 0
    out[valid] = cumulative[indices[valid]] / total * 100.0
    return out


def iter_coords(geometry: dict) -> Iterable[Tuple[np.ndarray, np.ndarray]]:
    geom_type = geometry.get("type")
    coords = geometry.get("coordinates")
    if not coords:
        return
    if geom_type == "Point":
        yield np.asarray([coords[0]], dtype=float), np.asarray([coords[1]], dtype=float)
    elif geom_type == "MultiPoint":
        xy = np.asarray(coords, dtype=float)
        yield xy[:, 0], xy[:, 1]
    elif geom_type == "LineString":
        xy = np.asarray(coords, dtype=float)
        if xy.ndim == 2 and xy.shape[1] >= 2:
            yield xy[:, 0], xy[:, 1]
    elif geom_type == "MultiLineString":
        for part in coords:
            xy = np.asarray(part, dtype=float)
            if xy.ndim == 2 and xy.shape[1] >= 2:
                yield xy[:, 0], xy[:, 1]
    elif geom_type == "Polygon":
        for ring in coords:
            xy = np.asarray(ring, dtype=float)
            if xy.ndim == 2 and xy.shape[1] >= 2:
                yield xy[:, 0], xy[:, 1]
    elif geom_type == "MultiPolygon":
        for polygon in coords:
            for ring in polygon:
                xy = np.asarray(ring, dtype=float)
                if xy.ndim == 2 and xy.shape[1] >= 2:
                    yield xy[:, 0], xy[:, 1]


def collect_geojson_bounds(path: Optional[Path]) -> Optional[Tuple[float, float, float, float]]:
    if path is None or not path.exists():
        return None
    xs: List[float] = []
    ys: List[float] = []
    try:
        data = load_json(path)
    except Exception:
        return None
    for feature in data.get("features") or []:
        for x, y in iter_coords(feature.get("geometry") or {}):
            xs.extend(x.tolist())
            ys.extend(y.tolist())
    if not xs:
        return None
    return min(xs), min(ys), max(xs), max(ys)


def plot_geojson_map(ax: plt.Axes, model: GTCModel) -> None:
    ax.set_facecolor("#eef7fb")

    bounds_list: List[Tuple[float, float, float, float]] = []
    for path in (model.active_domain_geojson, model.cama_geojson):
        bounds = collect_geojson_bounds(path)
        if bounds is not None:
            bounds_list.append(bounds)

    if model.active_domain_geojson and model.active_domain_geojson.exists():
        try:
            data = load_json(model.active_domain_geojson)
            for feature in data.get("features") or []:
                for x, y in iter_coords(feature.get("geometry") or {}):
                    ax.fill(x, y, facecolor="#d9d9d9", edgecolor="#111111", linewidth=1.0, alpha=0.35, zorder=1)
                    ax.plot(x, y, color="#111111", linewidth=1.0, zorder=3)
        except Exception as exc:
            print(f"[WARN] {model.gtc_id}: could not plot active domain: {exc}")

    if model.cama_geojson and model.cama_geojson.exists():
        try:
            data = load_json(model.cama_geojson)
            for feature in data.get("features") or []:
                ftype = feature_type(feature)
                geom = feature.get("geometry") or {}
                if ftype == "River_Line":
                    for x, y in iter_coords(geom):
                        ax.plot(x, y, color="#306b8c", linewidth=0.75, alpha=0.85, zorder=2)
                elif ftype == "Snap_Link":
                    for x, y in iter_coords(geom):
                        ax.plot(x, y, color="#8a8a8a", linewidth=0.6, alpha=0.55, zorder=2)
            inlet_x = [cell.lon for cell in model.cells]
            inlet_y = [cell.lat for cell in model.cells]
            if inlet_x:
                ax.scatter(inlet_x, inlet_y, s=38, c="#d62728", edgecolors="white", linewidths=0.7, zorder=5)
        except Exception as exc:
            print(f"[WARN] {model.gtc_id}: could not plot upstream rivers: {exc}")

    if bounds_list:
        xmin = min(b[0] for b in bounds_list)
        ymin = min(b[1] for b in bounds_list)
        xmax = max(b[2] for b in bounds_list)
        ymax = max(b[3] for b in bounds_list)
        pad_x = max(0.25, 0.08 * (xmax - xmin))
        pad_y = max(0.25, 0.08 * (ymax - ymin))
        ax.set_xlim(xmin - pad_x, xmax + pad_x)
        ax.set_ylim(ymin - pad_y, ymax + pad_y)

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, color="#d0d0d0", linewidth=0.45, alpha=0.6)
    ax.tick_params(direction="in", top=True, right=True)
    handles = [
        Line2D([0], [0], color="#111111", lw=1.1, label="SFINCS domain"),
        Line2D([0], [0], color="#306b8c", lw=1.1, label="Upstream CaMa rivers"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#d62728", markeredgecolor="white", markersize=6, label="CaMa inlet cells"),
    ]
    ax.legend(handles=handles, loc="best", frameon=True, framealpha=0.92)
    ax.text(0.02, 0.98, f"a  {model.gtc_id}", transform=ax.transAxes, ha="left", va="top", fontweight="bold", fontsize=15)


def month_spans_365(months: Sequence[int]) -> List[Tuple[int, int, int]]:
    starts = np.r_[1, 1 + np.cumsum(MONTH_LENGTHS_365[:-1])]
    ends = np.cumsum(MONTH_LENGTHS_365)
    return [(int(starts[m - 1]), int(ends[m - 1]), int(m)) for m in months if 1 <= int(m) <= 12]


def format_flow_axis(ax: plt.Axes) -> None:
    ax.tick_params(direction="in", top=True, right=True)
    ax.grid(True, color="#dcdcdc", linewidth=0.55, alpha=0.7)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def plot_tc_month_bars_on_twin_axis(
    ax: plt.Axes,
    tc_weight_info: TCWeightInfo,
) -> Optional[plt.Axes]:
    weights = {int(m): float(w) for m, w in tc_weight_info.month_weights.items() if 1 <= int(m) <= 12 and float(w) > 0}
    if not weights:
        return None

    starts = np.r_[1, 1 + np.cumsum(MONTH_LENGTHS_365[:-1])]
    centers = starts + MONTH_LENGTHS_365 / 2
    probabilities = np.array([weights.get(month, 0.0) for month in range(1, 13)], dtype=np.float64)

    tc_ax = ax.twinx()
    tc_ax.set_zorder(ax.get_zorder() - 1)
    ax.patch.set_alpha(0.0)
    tc_ax.patch.set_alpha(0.0)
    tc_ax.bar(
        centers,
        probabilities,
        width=MONTH_LENGTHS_365 * 0.78,
        color="#9b59b6",
        edgecolor="#6c3483",
        linewidth=0.65,
        alpha=0.28,
        label="TC landfall probability",
        zorder=0,
    )

    ymax = float(np.nanmax(probabilities)) if np.isfinite(probabilities).any() else 0.0
    tc_ax.set_ylim(0.0, ymax * 1.18 if ymax > 0 else 1.0)
    tc_ax.set_ylabel("TC landfall probability", color="#6c3483")
    tc_ax.tick_params(axis="y", direction="in", left=False, right=True, colors="#6c3483")
    tc_ax.spines["right"].set_color("#6c3483")
    tc_ax.spines["right"].set_linewidth(1.0)
    tc_ax.grid(False)
    return tc_ax


def plot_daily_seasonality(ax: plt.Axes, series: Sequence[GTCSeries], tc_weight_info: TCWeightInfo) -> None:
    grouped = group_series(series)
    x = np.arange(1, 366)
    tc_ax = plot_tc_month_bars_on_twin_axis(ax, tc_weight_info)

    hist = grouped.get("historical", [])
    if hist:
        ax.plot(x, hist[0].daily_clim, color=SCENARIO_COLORS["historical"], linewidth=2.0, label="Historical", zorder=3)

    for scenario in sorted(k for k in grouped if k != "historical"):
        items = grouped[scenario]
        color = SCENARIO_COLORS.get(scenario, "#666666")
        mean = nanmean_stack([item.daily_clim for item in items])
        lo = nanpercentile_stack([item.daily_clim for item in items], 10)
        hi = nanpercentile_stack([item.daily_clim for item in items], 90)
        if items and len(items) > 1:
            ax.fill_between(x, lo, hi, color=color, alpha=0.14, linewidth=0, zorder=1)
        ax.plot(x, mean, color=color, linewidth=1.8, label=scenario, zorder=2)

    starts = np.r_[1, 1 + np.cumsum(MONTH_LENGTHS_365[:-1])]
    centers = starts + MONTH_LENGTHS_365 / 2
    ax.set_xlim(1, 365)
    ax.set_xticks(centers[::2])
    ax.set_xticklabels([MONTH_LABELS[i] for i in range(0, 12, 2)])
    ax.set_ylabel("Discharge (m$^3$ s$^{-1}$)")
    ax.set_xlabel("Day of year")
    ax.text(0.02, 0.98, "b", transform=ax.transAxes, ha="left", va="top", fontweight="bold", fontsize=15)
    handles, labels = ax.get_legend_handles_labels()
    if tc_ax is not None:
        tc_handles, tc_labels = tc_ax.get_legend_handles_labels()
        handles.extend(tc_handles)
        labels.extend(tc_labels)
    if handles:
        ax.legend(handles, labels, loc="best", frameon=False, ncol=2)
    format_flow_axis(ax)
    if tc_ax is not None:
        ax.tick_params(axis="y", right=False)


def plot_monthly_inlets(ax: plt.Axes, model: GTCModel, series: Sequence[GTCSeries]) -> None:
    grouped = group_series(series)
    x = np.arange(1, 13)
    hist = grouped.get("historical", [])
    if hist:
        hist_item = hist[0]
        if hist_item.inlet_monthly.shape[0] > 1:
            for inlet_values in hist_item.inlet_monthly:
                ax.plot(x, inlet_values, color="#9a9a9a", linewidth=0.8, alpha=0.65)
        ax.plot(x, hist_item.monthly_clim, color=SCENARIO_COLORS["historical"], linewidth=2.0, label="Historical total")

    for scenario in sorted(k for k in grouped if k != "historical"):
        items = grouped[scenario]
        color = SCENARIO_COLORS.get(scenario, "#666666")
        mean = nanmean_stack([item.monthly_clim for item in items])
        ax.plot(x, mean, color=color, linewidth=1.8, label=f"{scenario} total")

    ax.set_xlim(1, 12)
    ax.set_xticks(x)
    ax.set_xticklabels(MONTH_LABELS, rotation=35, ha="right")
    ax.set_ylabel("Monthly discharge (m$^3$ s$^{-1}$)")
    ax.set_xlabel("Month")
    ax.text(0.02, 0.98, "c", transform=ax.transAxes, ha="left", va="top", fontweight="bold", fontsize=15)
    ax.legend(loc="best", frameon=False, fontsize=8)
    format_flow_axis(ax)
    if len(model.cells) > 1:
        ax.text(
            0.98,
            0.05,
            f"grey lines: {len(model.cells)} inlet cells",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color="#555555",
        )


def plot_tc_cdf(ax: plt.Axes, series: Sequence[GTCSeries]) -> None:
    grouped = group_series(series)
    sample_arrays = [item.tc_samples for item in series if item.tc_samples.size]
    if not sample_arrays:
        ax.text(0.5, 0.5, "No TC-season samples", transform=ax.transAxes, ha="center", va="center")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 100)
    else:
        all_samples = np.concatenate(sample_arrays)
        x_max = float(np.nanpercentile(all_samples, 99.5))
        if not math.isfinite(x_max) or x_max <= 0:
            x_max = max(1.0, float(np.nanmax(all_samples)))
        x_grid = np.linspace(0.0, x_max * 1.05, 300)
        for scenario in ["historical"] + sorted(k for k in grouped if k != "historical"):
            items = grouped.get(scenario, [])
            sample_list = [item.tc_samples for item in items if item.tc_samples.size]
            weight_list = [item.tc_sample_weights for item in items if item.tc_samples.size]
            if not sample_list:
                continue
            samples = np.concatenate(sample_list)
            weights = np.concatenate(weight_list)
            color = SCENARIO_COLORS.get(scenario, "#666666")
            cdf = weighted_cdf(samples, weights, x_grid)
            label = "Historical" if scenario == "historical" else scenario
            ax.plot(x_grid, cdf, color=color, linewidth=1.9, label=label)
            p50 = weighted_quantile(samples, weights, 50)
            p90 = weighted_quantile(samples, weights, 90)
            ax.plot(p50, weighted_cdf(samples, weights, np.asarray([p50]))[0], "s", ms=5, mfc="white", mec=color, mew=1.2)
            ax.plot(p90, weighted_cdf(samples, weights, np.asarray([p90]))[0], "^", ms=6, mfc="white", mec=color, mew=1.2)
        ax.set_xlim(0, x_grid[-1])
        ax.set_ylim(0, 100)
    ax.set_ylabel("TC-weighted CDF (%)")
    ax.set_xlabel("Discharge (m$^3$ s$^{-1}$)")
    ax.text(0.02, 0.98, "d", transform=ax.transAxes, ha="left", va="top", fontweight="bold", fontsize=15)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="best", frameon=False, ncol=2)
    format_flow_axis(ax)


def save_gtc_figure(
    model: GTCModel,
    series: Sequence[GTCSeries],
    tc_weight_info: TCWeightInfo,
    figure_dir: Path,
    dpi: int,
) -> Path:
    fig = plt.figure(figsize=(13.2, 8.7), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, left=0.065, right=0.985, bottom=0.075, top=0.965, wspace=0.22, hspace=0.30)
    ax_map = fig.add_subplot(gs[0, 0])
    ax_daily = fig.add_subplot(gs[0, 1])
    ax_monthly = fig.add_subplot(gs[1, 0])
    ax_cdf = fig.add_subplot(gs[1, 1])

    plot_geojson_map(ax_map, model)
    plot_daily_seasonality(ax_daily, series, tc_weight_info)
    plot_monthly_inlets(ax_monthly, model, series)
    plot_tc_cdf(ax_cdf, series)

    png_path = figure_dir / f"{model.gtc_id}_cama_inflow_seasonality_tc_cdf.png"
    pdf_path = figure_dir / f"{model.gtc_id}_cama_inflow_seasonality_tc_cdf.pdf"
    fig.savefig(png_path, dpi=dpi)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path


def rows_for_tables(
    model: GTCModel,
    series: Sequence[GTCSeries],
    tc_weight_info: TCWeightInfo,
) -> Tuple[
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
]:
    daily_rows: List[Dict[str, object]] = []
    monthly_rows: List[Dict[str, object]] = []
    inlet_monthly_rows: List[Dict[str, object]] = []
    inlet_summary_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    for item in series:
        common = {
            "gtc_id": model.gtc_id,
            "case": item.case.name,
            "scenario": item.case.scenario,
            "model": item.case.model,
        }
        for day, value in enumerate(item.daily_clim, start=1):
            daily_rows.append({**common, "clim_day": day, "mean_flow_m3s": value})
        for month, value in enumerate(item.monthly_clim, start=1):
            monthly_rows.append({**common, "month": month, "mean_flow_m3s": value})
        for inlet_idx, cell in enumerate(model.cells):
            for month, value in enumerate(item.inlet_monthly[inlet_idx, :], start=1):
                inlet_monthly_rows.append(
                    {
                        **common,
                        "inlet_id": cell.inlet_id,
                        "cell_id": cell.cell_id,
                        "cama_row": cell.row,
                        "cama_col": cell.col,
                        "month": month,
                        "mean_flow_m3s": value,
                    }
                )
            inlet_summary_rows.append(
                {
                    **common,
                    "inlet_order": inlet_idx + 1,
                    "inlet_id": cell.inlet_id,
                    "cell_id": cell.cell_id,
                    "cama_row": cell.row,
                    "cama_col": cell.col,
                    "cama_lon": cell.lon,
                    "cama_lat": cell.lat,
                    "uparea_km2": cell.uparea_km2,
                    "q_reference_m3s": cell.q_reference_m3s,
                    "tc_weight_source": tc_weight_info.source,
                    "tc_event_count": tc_weight_info.event_count,
                    "mean_flow_m3s": float(item.inlet_mean_flow_m3s[inlet_idx]),
                    "tc_p50_flow_m3s": float(item.inlet_p50_tc_m3s[inlet_idx]),
                    "tc_p90_flow_m3s": float(item.inlet_p90_tc_m3s[inlet_idx]),
                    "basin_id": model.basin_id,
                    "basin_label": model.basin_label,
                }
            )
        summary_rows.append(
            {
                **common,
                "n_cama_cells": len(model.cells),
                "tc_basin_id": tc_weight_info.basin_id,
                "tc_basin_label": tc_weight_info.basin_label,
                "tc_weight_source": tc_weight_info.source,
                "tc_event_count": tc_weight_info.event_count,
                "tc_months": ",".join(str(m) for m in tc_weight_info.months),
                "tc_month_weights": ";".join(f"{m}:{tc_weight_info.month_weights.get(m, 0.0):.6f}" for m in range(1, 13)),
                "mean_flow_m3s": item.mean_flow_m3s,
                "tc_p50_flow_m3s": item.p50_tc_m3s,
                "tc_p90_flow_m3s": item.p90_tc_m3s,
                "basin_id": model.basin_id,
                "basin_label": model.basin_label,
            }
        )
    return daily_rows, monthly_rows, inlet_monthly_rows, inlet_summary_rows, summary_rows


def weighted_inlet_statistics(
    dates: pd.DatetimeIndex,
    values: np.ndarray,
    month_weights: Dict[int, float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=np.float64)
    if values.ndim == 1:
        values = values[:, None]
    finite = np.isfinite(values)
    counts = finite.sum(axis=0)
    means = np.divide(
        np.nansum(values, axis=0),
        counts,
        out=np.full(values.shape[1], np.nan, dtype=np.float64),
        where=counts > 0,
    )

    months = pd.DatetimeIndex(dates).month.to_numpy(dtype=np.int8)
    sample_weights = np.zeros(values.shape, dtype=np.float64)
    for month, month_weight in month_weights.items():
        if month_weight <= 0:
            continue
        row_mask = months == int(month)
        if not row_mask.any():
            continue
        valid_month = finite[row_mask, :]
        valid_counts = valid_month.sum(axis=0)
        scale = np.divide(
            float(month_weight),
            valid_counts,
            out=np.zeros(values.shape[1], dtype=np.float64),
            where=valid_counts > 0,
        )
        sample_weights[row_mask, :] = valid_month * scale[None, :]

    sortable = np.where(finite, values, np.inf)
    order = np.argsort(sortable, axis=0)
    sorted_values = np.take_along_axis(sortable, order, axis=0)
    sorted_weights = np.take_along_axis(sample_weights, order, axis=0)
    cumulative = np.cumsum(sorted_weights, axis=0)
    totals = cumulative[-1, :]
    p50 = np.full(values.shape[1], np.nan, dtype=np.float64)
    p90 = np.full(values.shape[1], np.nan, dtype=np.float64)
    for col in range(values.shape[1]):
        valid_sorted = np.isfinite(sorted_values[:, col]) & (sorted_weights[:, col] > 0)
        if not valid_sorted.any() or totals[col] <= 0:
            continue
        cumulative_col = cumulative[valid_sorted, col]
        values_col = sorted_values[valid_sorted, col]
        p50[col] = float(np.interp(0.50 * totals[col], cumulative_col, values_col))
        p90[col] = float(np.interp(0.90 * totals[col], cumulative_col, values_col))
    return means, p50, p90


def write_inlet_summary_only(
    models: Sequence[GTCModel],
    extracts: Sequence[CaseExtract],
    cell_to_index: Dict[Tuple[int, int], int],
    tc_weight_infos: Dict[str, TCWeightInfo],
    path: Path,
) -> pd.DataFrame:
    grouped: Dict[Tuple[Tuple[int, float], ...], List[Tuple[GTCModel, int, InletCell, int, TCWeightInfo]]] = {}
    for model in models:
        tc_weight_info = tc_weight_infos.get(model.gtc_id)
        if tc_weight_info is None:
            raise KeyError(f"Missing TC month weights for {model.gtc_id}")
        signature = tuple(sorted((int(month), float(weight)) for month, weight in tc_weight_info.month_weights.items()))
        for inlet_idx, cell in enumerate(model.cells):
            col_index = cell_to_index[(cell.row, cell.col)]
            grouped.setdefault(signature, []).append(
                (model, inlet_idx, cell, col_index, tc_weight_info)
            )

    rows: List[Dict[str, object]] = []
    for case_idx, case_extract in enumerate(extracts, start=1):
        print(f"[INLET] case {case_idx}/{len(extracts)} {case_extract.case.name}")
        for signature, entries in grouped.items():
            month_weights = {month: weight for month, weight in signature}
            col_indices = [entry[3] for entry in entries]
            values = case_extract.values[:, col_indices]
            mean_flow, p50, p90 = weighted_inlet_statistics(
                case_extract.dates,
                values,
                month_weights,
            )
            for entry_idx, (model, inlet_idx, cell, _, tc_weight_info) in enumerate(entries):
                rows.append(
                    {
                        "gtc_id": model.gtc_id,
                        "case": case_extract.case.name,
                        "scenario": case_extract.case.scenario,
                        "model": case_extract.case.model,
                        "inlet_order": inlet_idx + 1,
                        "inlet_id": cell.inlet_id,
                        "cell_id": cell.cell_id,
                        "cama_row": cell.row,
                        "cama_col": cell.col,
                        "cama_lon": cell.lon,
                        "cama_lat": cell.lat,
                        "uparea_km2": cell.uparea_km2,
                        "q_reference_m3s": cell.q_reference_m3s,
                        "tc_weight_source": tc_weight_info.source,
                        "tc_event_count": tc_weight_info.event_count,
                        "mean_flow_m3s": float(mean_flow[entry_idx]),
                        "tc_p50_flow_m3s": float(p50[entry_idx]),
                        "tc_p90_flow_m3s": float(p90[entry_idx]),
                        "basin_id": model.basin_id,
                        "basin_label": model.basin_label,
                    }
                )
    result = pd.DataFrame(rows)
    result.to_csv(path, index=False, encoding="utf-8-sig")
    return result


def run_gtc_diagnostics(args: argparse.Namespace) -> None:
    table_dir = ensure_dir(args.out_dir / "tables")
    cache_dir = ensure_dir(args.out_dir / "cache")
    figure_dir = ensure_dir(args.out_dir / "figures")

    grid = parse_cama_params(args.map_dir)
    cases = discover_cases(args)
    models_all, status_df = discover_gtc_models(args, grid)
    models = selected_models(models_all)
    unique_cells = unique_cells_from_models(models)
    cell_to_index = {cell: i for i, cell in enumerate(unique_cells)}
    cache_key = extraction_hash(unique_cells, args)
    tc_weight_infos, tc_weight_df, _ = build_tc_weight_infos(models, args, table_dir)

    status_df.to_csv(table_dir / "gtc_status.csv", index=False, encoding="utf-8-sig")
    inlet_df = write_inlet_table(models, table_dir / "gtc_cama_inlet_cells.csv")
    print(f"[INFO] Matched GTC partitions: {len(models)}")
    print(f"[INFO] Unique CaMa cells to extract: {len(unique_cells)}")
    print(f"[INFO] Case count: {len(cases)}")
    print(f"[INFO] Inlet table -> {table_dir / 'gtc_cama_inlet_cells.csv'}")
    print(f"[INFO] TC month weights -> {table_dir / 'gtc_tc_month_weights.csv'} ({len(tc_weight_df)} rows)")

    extracts: List[CaseExtract] = []
    for case in cases:
        extracts.append(extract_case_cells(case, unique_cells, grid, args, cache_dir, cache_key))

    if args.tables_only:
        inlet_summary_path = table_dir / "gtc_flow_summary_by_inlet_cell.csv"
        inlet_summary = write_inlet_summary_only(
            models,
            extracts,
            cell_to_index,
            tc_weight_infos,
            inlet_summary_path,
        )
        print(f"[DONE] Independent inlet summary -> {inlet_summary_path}")
        print(f"[DONE] Independent inlet rows: {len(inlet_summary)}")
        return

    daily_rows: List[Dict[str, object]] = []
    monthly_rows: List[Dict[str, object]] = []
    inlet_monthly_rows: List[Dict[str, object]] = []
    inlet_summary_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    figure_rows: List[Dict[str, object]] = []

    for idx, model in enumerate(models, start=1):
        print(f"[PLOT] {idx}/{len(models)} {model.gtc_id}: {len(model.cells)} CaMa cell(s)")
        tc_weight_info = tc_weight_infos.get(model.gtc_id, fixed_month_tc_weight_info(model, args, source="fixed_months_fallback"))
        series = [gtc_series_from_case(model, case_extract, cell_to_index, tc_weight_info) for case_extract in extracts]
        d_rows, m_rows, im_rows, is_rows, s_rows = rows_for_tables(
            model, series, tc_weight_info
        )
        daily_rows.extend(d_rows)
        monthly_rows.extend(m_rows)
        inlet_monthly_rows.extend(im_rows)
        inlet_summary_rows.extend(is_rows)
        summary_rows.extend(s_rows)
        if not args.tables_only:
            png_path = save_gtc_figure(model, series, tc_weight_info, figure_dir, args.dpi)
            figure_rows.append(
                {
                    "gtc_id": model.gtc_id,
                    "figure_png": str(png_path),
                    "n_cama_cells": len(model.cells),
                    "tc_basin_id": tc_weight_info.basin_id,
                    "tc_weight_source": tc_weight_info.source,
                    "tc_event_count": tc_weight_info.event_count,
                }
            )

    pd.DataFrame(daily_rows).to_csv(table_dir / "gtc_daily_climatology_total.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(monthly_rows).to_csv(table_dir / "gtc_monthly_climatology_total.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(inlet_monthly_rows).to_csv(table_dir / "gtc_monthly_climatology_by_inlet_cell.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(inlet_summary_rows).to_csv(
        table_dir / "gtc_flow_summary_by_inlet_cell.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(summary_rows).to_csv(table_dir / "gtc_flow_summary_by_case.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(figure_rows).to_csv(table_dir / "gtc_figure_index.csv", index=False, encoding="utf-8-sig")

    print("[DONE] Outputs:")
    print(f"  figures: {figure_dir}")
    print(f"  tables : {table_dir}")
    print(f"  cache  : {cache_dir}")
    print(f"  matched inlet rows: {len(inlet_df)}")


def clean_global_months(months: Sequence[int]) -> Tuple[int, ...]:
    out = sorted({int(month) for month in months if 1 <= int(month) <= 12})
    if not out:
        raise ValueError("TC-season month list is empty after validation.")
    return tuple(out)


def global_tc_months_for_lat(
    lat: float,
    months_nh: Sequence[int],
    months_sh: Sequence[int],
) -> Tuple[int, ...]:
    return tuple(months_nh if float(lat) >= 0.0 else months_sh)


def global_weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    good = values.notna() & weights.notna() & (weights > 0)
    if not bool(good.any()):
        return math.nan
    return float(
        np.average(
            values.loc[good].astype(float).to_numpy(),
            weights=weights.loc[good].astype(float).to_numpy(),
        )
    )


def global_scenario_sort_key(scenario: str) -> Tuple[int, str]:
    scenario = str(scenario)
    if scenario in GLOBAL_SCENARIO_ORDER:
        return GLOBAL_SCENARIO_ORDER.index(scenario), scenario
    return len(GLOBAL_SCENARIO_ORDER), scenario


def build_global_tc_season_table(
    monthly: pd.DataFrame,
    outlets: pd.DataFrame,
    months_nh: Sequence[int],
    months_sh: Sequence[int],
) -> pd.DataFrame:
    monthly_required = {
        "case",
        "scenario",
        "model",
        "outlet_id",
        "month",
        "mean_flow_m3s",
        "n_steps",
    }
    outlet_required = {"outlet_id", "rank", "lon", "lat", "uparea_km2"}
    missing = monthly_required.difference(monthly.columns)
    if missing:
        raise ValueError(f"Missing columns in monthly table: {sorted(missing)}")
    missing = outlet_required.difference(outlets.columns)
    if missing:
        raise ValueError(f"Missing columns in outlet table: {sorted(missing)}")

    outlet_meta = outlets[
        ["outlet_id", "rank", "lon", "lat", "uparea_km2"]
    ].drop_duplicates("outlet_id")
    rows: List[Dict[str, object]] = []
    for _, outlet in outlet_meta.iterrows():
        outlet_id = outlet["outlet_id"]
        lat = float(outlet["lat"])
        months = global_tc_months_for_lat(lat, months_nh, months_sh)
        subset = monthly[
            (monthly["outlet_id"] == outlet_id) & monthly["month"].isin(months)
        ].copy()
        if subset.empty:
            continue

        historical = subset[subset["scenario"] == "historical"]
        historical_flow = global_weighted_mean(
            historical["mean_flow_m3s"], historical["n_steps"]
        )
        if not np.isfinite(historical_flow):
            continue

        scenarios = subset.loc[
            subset["scenario"] != "historical", "scenario"
        ].dropna().unique()
        for scenario in sorted(scenarios, key=global_scenario_sort_key):
            future = subset[subset["scenario"] == scenario]
            future_flow = global_weighted_mean(
                future["mean_flow_m3s"], future["n_steps"]
            )
            if not np.isfinite(future_flow):
                continue
            change_m3s = future_flow - historical_flow
            change_percent = (
                change_m3s / historical_flow * 100.0
                if historical_flow != 0
                else math.nan
            )
            rows.append(
                {
                    "scenario": scenario,
                    "outlet_id": outlet_id,
                    "rank": int(outlet["rank"]),
                    "lon": float(outlet["lon"]),
                    "lat": lat,
                    "uparea_km2": float(outlet["uparea_km2"]),
                    "hemisphere": "NH" if lat >= 0 else "SH",
                    "tc_months": ",".join(str(month) for month in months),
                    "historical_tc_flow_m3s": historical_flow,
                    "future_tc_flow_m3s": future_flow,
                    "change_m3s": change_m3s,
                    "change_percent": change_percent,
                    "n_future_models": int(future["case"].nunique()),
                }
            )

    if not rows:
        raise RuntimeError("No TC-season outlet rows could be built from the monthly table.")
    return pd.DataFrame(rows).sort_values(["scenario", "rank"]).reset_index(drop=True)


def global_marker_sizes(
    flows: Iterable[float],
    size_min: float,
    size_max: float,
) -> np.ndarray:
    flow = np.asarray(list(flows), dtype=float)
    finite = np.isfinite(flow) & (flow > 0)
    out = np.full(flow.shape, size_min, dtype=float)
    if not finite.any():
        return out
    logged = np.log10(flow[finite])
    low = float(np.nanmin(logged))
    high = float(np.nanmax(logged))
    if high <= low:
        out[finite] = 0.5 * (size_min + size_max)
    else:
        out[finite] = size_min + (logged - low) / (high - low) * (size_max - size_min)
    return out


def global_size_for_legend(
    value: float,
    flow_min: float,
    flow_max: float,
    size_min: float,
    size_max: float,
) -> float:
    return float(global_marker_sizes([value, flow_min, flow_max], size_min, size_max)[0])


def global_pretty_flow(value: float) -> str:
    if value >= 1000:
        return f"{value:,.0f}"
    return f"{value:g}"


def choose_global_size_legend_values(flows: np.ndarray) -> List[float]:
    finite = flows[np.isfinite(flows) & (flows > 0)]
    if finite.size == 0:
        return [1000.0, 10000.0, 100000.0]
    candidates = [1000.0, 5000.0, 10000.0, 50000.0, 100000.0, 300000.0]
    selected = [
        value
        for value in candidates
        if finite.min() * 0.6 <= value <= finite.max() * 1.4
    ]
    if len(selected) >= 3:
        return selected[:3] if selected[0] >= np.nanmedian(finite) else selected[-3:]
    return [float(value) for value in np.nanpercentile(finite, [25, 50, 90])]


def add_global_basemap(ax) -> None:
    ax.set_global()
    ax.set_facecolor("#dff1f8")
    ax.add_feature(
        cfeature.OCEAN.with_scale("110m"),
        facecolor="#dff1f8",
        edgecolor="none",
        zorder=0,
    )
    ax.add_feature(
        cfeature.LAND.with_scale("110m"),
        facecolor="#f1efe6",
        edgecolor="#b7b7b7",
        linewidth=0.35,
        zorder=1,
    )
    ax.add_feature(
        cfeature.COASTLINE.with_scale("110m"),
        edgecolor="#7c7c7c",
        linewidth=0.35,
        zorder=2,
    )
    ax.add_feature(
        cfeature.BORDERS.with_scale("110m"),
        edgecolor="#c6c6c6",
        linewidth=0.25,
        zorder=2,
    )
    gridlines = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.35,
        color="#b8c7cf",
        alpha=0.7,
        linestyle="-",
    )
    gridlines.top_labels = False
    gridlines.right_labels = False
    gridlines.xlabel_style = {"size": 9}
    gridlines.ylabel_style = {"size": 9}


def plot_global_tc_map(
    tc_table: pd.DataFrame,
    args: argparse.Namespace,
) -> Tuple[Path, Path]:
    scenarios = sorted(
        tc_table["scenario"].dropna().unique(),
        key=global_scenario_sort_key,
    )
    figure = plt.figure(figsize=(6.3 * len(scenarios), 4.7), constrained_layout=True)
    axes = [
        figure.add_subplot(
            1,
            len(scenarios),
            index + 1,
            projection=ccrs.Robinson(central_longitude=0),
        )
        for index in range(len(scenarios))
    ]

    changes = tc_table["change_percent"].astype(float).to_numpy()
    if args.global_change_vlim is None:
        finite = np.abs(changes[np.isfinite(changes)])
        vmax = float(np.nanpercentile(finite, 95)) if finite.size else 10.0
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = 10.0
        vmax = max(5.0, vmax)
    else:
        vmax = float(args.global_change_vlim)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    all_flows = tc_table["future_tc_flow_m3s"].astype(float).to_numpy()
    positive_flows = all_flows[np.isfinite(all_flows) & (all_flows > 0)]
    flow_min = float(np.nanmin(positive_flows))
    flow_max = float(np.nanmax(positive_flows))
    last_scatter = None

    for ax, scenario in zip(axes, scenarios):
        add_global_basemap(ax)
        data = tc_table[tc_table["scenario"] == scenario].copy()
        sizes = global_marker_sizes(
            data["future_tc_flow_m3s"],
            args.global_size_min,
            args.global_size_max,
        )
        last_scatter = ax.scatter(
            data["lon"],
            data["lat"],
            s=sizes,
            c=data["change_percent"],
            cmap="RdBu_r",
            norm=norm,
            edgecolor="#222222",
            linewidth=0.35,
            alpha=0.88,
            transform=ccrs.PlateCarree(),
            zorder=5,
        )
        ax.text(
            0.02,
            0.98,
            scenario,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=14,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=2.5),
        )
        if args.label_outlets:
            for _, row in data.iterrows():
                ax.text(
                    row["lon"] + 2.0,
                    row["lat"] + 1.0,
                    str(row["outlet_id"]),
                    transform=ccrs.PlateCarree(),
                    fontsize=7,
                    color="#111111",
                    zorder=6,
                )

    if last_scatter is not None:
        colorbar = figure.colorbar(
            last_scatter,
            ax=axes,
            orientation="horizontal",
            pad=0.055,
            shrink=0.48,
            aspect=34,
        )
        colorbar.set_label(
            "Future change relative to historical TC-season discharge (%)"
        )
        colorbar.ax.tick_params(direction="in", length=3)

    legend_values = choose_global_size_legend_values(all_flows)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="white",
            markeredgecolor="#333333",
            markersize=math.sqrt(
                global_size_for_legend(
                    value,
                    flow_min,
                    flow_max,
                    args.global_size_min,
                    args.global_size_max,
                )
            ),
            label=global_pretty_flow(value),
        )
        for value in legend_values
    ]
    axes[0].legend(
        handles=handles,
        title="TC-season discharge\n(m$^3$ s$^{-1}$)",
        loc="lower left",
        frameon=True,
        framealpha=0.88,
        borderpad=0.8,
        handletextpad=1.5,
    )

    figure_dir = ensure_dir(args.global_figure_dir)
    png_path = figure_dir / f"{args.global_output_name}.png"
    pdf_path = figure_dir / f"{args.global_output_name}.pdf"
    figure.savefig(png_path, dpi=args.global_dpi)
    figure.savefig(pdf_path)
    plt.close(figure)
    return png_path, pdf_path


def run_global_outlet_map(args: argparse.Namespace) -> None:
    table_dir = ensure_dir(args.global_table_dir)
    monthly_path = table_dir / args.global_monthly_table
    outlet_path = table_dir / args.global_outlet_table
    if not monthly_path.exists():
        raise FileNotFoundError(monthly_path)
    if not outlet_path.exists():
        raise FileNotFoundError(outlet_path)

    monthly = pd.read_csv(monthly_path)
    outlets = pd.read_csv(outlet_path)
    tc_table = build_global_tc_season_table(
        monthly,
        outlets,
        clean_global_months(args.global_months_nh),
        clean_global_months(args.global_months_sh),
    )

    output_table = table_dir / "major_outlet_tc_season_flow_change_by_scenario.csv"
    tc_table.to_csv(output_table, index=False, encoding="utf-8-sig")
    png_path, pdf_path = plot_global_tc_map(tc_table, args)
    print("[DONE] TC-season outlet flow table:")
    print(f"  {output_table}")
    print("[DONE] TC-season outlet flow map:")
    print(f"  {png_path}")
    print(f"  {pdf_path}")


def main() -> None:
    args = parse_args()
    if args.hist_start_year > args.hist_end_year:
        raise ValueError("hist_start_year must not be later than hist_end_year")
    if args.future_start_year > args.future_end_year:
        raise ValueError("future_start_year must not be later than future_end_year")
    print(
        f"[PERIOD] ERA5 historical: {args.hist_start_year}-{args.hist_end_year}; "
        f"CMIP6 future: {args.future_start_year}-{args.future_end_year}"
    )
    print("[STEP 1/2] Extracting bin data and plotting every matched GTC...")
    run_gtc_diagnostics(args)
    if not args.tables_only:
        print("[STEP 2/2] Plotting the global TC-season outlet distribution...")
        run_global_outlet_map(args)
    print("[DONE] Combined GTC and global-map workflow completed.")


if __name__ == "__main__":
    main()
