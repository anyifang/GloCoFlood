# -*- coding: utf-8 -*-
"""
Batch coupling for future synthetic-TC SFINCS cases.

This script is based on sfincs_couple_adcirc_cama_era5_PRD.py, but replaces the
ERA5 precipitation step with an R-CLIPER tropical-cyclone rain field generated
from fort22_meta.txt. It loops over one or more estuaries, scenarios and return
periods from the ADCIRC future_adirc_fort directory.

Typical use:
    python sfincs_couple_adcirc_cama_tc_rain_future_batch.py
    python sfincs_couple_adcirc_cama_tc_rain_future_batch.py --site PRD
    python sfincs_couple_adcirc_cama_tc_rain_future_batch.py --site PRD --scenarios ssp585_mean --return-periods rp100yr
    python sfincs_couple_adcirc_cama_tc_rain_future_batch.py --all-sites

Default output layout:
    examples/PRD_future/PRD__ssp585_mean_flowq50_rp100yr/
    examples/YRD_future/YRD__ssp585_mean_flowq50_rp100yr/
    examples/BoB_future/BoB__ssp585_mean_flowq50_rp100yr/
    examples/Misp_future/Mississippi__ssp585_mean_flowq50_rp100yr/

Notes:
    - Prefer ADCIRC fort.63 / fort.63.nc time-series files for sfincs.bzs.
    - If no fort.63 exists, maxele.63 can be used as a constant boundary only
      when ALLOW_STATIC_MAXELE63_BZS or --allow-static-maxele is enabled.
    - VLM subsidence is applied to sfincs.dep once per output model. By default
      the uniform lowering uses the Nature supplementary Table 1
      area-weighted VLM rate instead of averaging the gridVLM tif.
      vlm_subsidence_applied.txt marks completed lowering to avoid repeats.
"""

import argparse
import csv
import math
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray
import xarray as xr
from pyproj import Transformer
from scipy.spatial import cKDTree

import sfincs_couple_adcirc_cama_era5_PRD as base


# =============================================================================
# 0. User configuration
# =============================================================================

EXAMPLES_ROOT = Path(__file__).resolve().parent

LOCAL_FUTURE_ADCIRC_ROOT = Path(
    r"H:\Global_compoundflood\ADCIRC\OceanMesh2D-Projection\Examples"
    r"\cases_downscaled_ERA5_globalMeanCv\future_adirc_fort"
)
SERVER_FUTURE_ADCIRC_ROOT = Path(
    "/publicfs01/fs1-m8/home/m8s001451/zayf/ADCIRC/global_model/future_adirc_fort"
)

FUTURE_ADCIRC_ROOT = (
    LOCAL_FUTURE_ADCIRC_ROOT
    if LOCAL_FUTURE_ADCIRC_ROOT.exists()
    else SERVER_FUTURE_ADCIRC_ROOT
)
ADCIRC_RESULT_ROOT = FUTURE_ADCIRC_ROOT / "maxele63_results"

TARGET_SITES = ["PRD", "YRD", "BoB", "Mississippi"]
SCENARIO_TAGS = [
    "historical_era5",
    "ssp126_mean",
    "ssp245_mean",
    "ssp370_mean",
    "ssp585_mean",
]
RETURN_PERIOD_TAGS = ["rp200yr"]

OUTPUT_PARENT_ROOT = None
OVERWRITE_EXISTING_OUTPUT_MODEL = True
STOP_ON_ERROR = False

ENABLE_CAMA_DISCHARGE = True
SKIP_CAMA_IF_MISSING = False
FLOW_FORCING_MODE = "stats"  # "stats", "cama", or "none"; "stats_p90" is still accepted
FLOW_STATS_ROOT = Path(r"H:\Global_compoundflood\Camaflood\result")
FLOW_STATS_FILENAME_TEMPLATE = "sfincs_boundary_tc_flow_statistics__{site}.csv"
FLOW_STATS_PERCENTILE = 90  # 50 -> allocatedP50, 90 -> allocatedP90, 95 -> allocatedP95
FLOW_STATS_SCENARIO_MAP = {
    "historical_era5": "historical_gcm",
    "ssp126_mean": "ssp126",
    "ssp245_mean": "ssp245",
    "ssp370_mean": "ssp370",
    "ssp585_mean": "ssp585",
}
CAMA_RUN_SCRIPT_TEMPLATE = (
    "/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/cmf_v420_pkg/gosh/"
    "auto_cama_vic_coupled_module_{region_tag}.sh"
)
CAMA_OUTFLOW_BIN_TEMPLATE = (
    "/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/cmf_v420_pkg/gosh/"
    "work_{region_name}_vic_camaflood/camaout/outflw{year}.bin"
)
CAMA_OUTPUT_DT_HOURS = 24
CAMA_DTYPE = "<f4"
CAMA_ARRAY_ORDER = "tyx"
CAMA_TIME_START = None
CAMA_MAP_DIR = None
SPLIT_CAMA_DISCHARGE_FOR_DUPLICATE_SRC = True
CAMA_MATCH_MAX_DIST_DEG = 0.5

TARGET_FREQ = "1H"
ADCIRC_MATCH_MAX_DIST_DEG = 0.5
ALLOW_STATIC_MAXELE63_BZS = False

TC_RAIN_OUT_NC_NAME = "sfincs_tc_rain_rcliper.nc"
TC_RAIN_SOURCE_NC_NAME = "tc_rain_rcliper_source.nc"
RAIN_MAX_DIST_KM = 300.0
TC_RAIN_COMPRESSION_LEVEL = 9
TC_RAIN_LEAST_SIGNIFICANT_DIGIT = 2  # about 0.01 mm/hr precision, stored as float32
TC_RAIN_ZERO_BELOW_MMHR = 0.001
KEEP_TC_RAIN_SOURCE_NC = False
CLEAN_TEMP_RAIN_NC = True

APPLY_VLM_SUBSIDENCE = True
APPLY_VLM_TO_HISTORICAL = False
VLM_YEARS = 85.0
VLM_SOURCE = "paper_table"  # "paper_table" or "grid_tif"
PAPER_VLM_RATE_STAT = "mean"  # "mean", "median", or "area_weighted"
PAPER_VLM_RATE_UNIT = "mm/yr"
PAPER_VLM_TABLE_PATH = EXAMPLES_ROOT / "paper_delta_subsidence_supp_tables.xlsx"
VLM_RATE_UNIT = "cm/yr"  # Used only when VLM_SOURCE == "grid_tif"
VLM_ROOT = EXAMPLES_ROOT / "data" / "gridVLM"

PAPER_VLM_RATES_MM_YR: Dict[str, Dict[str, float]] = {
    # Nature supplementary Table 1: land subsidence / VLM rates for 40 deltas.
    # Negative VLM means downward land motion; SFINCS applies the magnitude.
    "Pearl": {
        "mean": -3.462510,
        "median": -2.703990,
        "area_weighted": -3.144292,
    },
    "Yangtze": {
        "mean": -1.815526,
        "median": -1.391256,
        "area_weighted": -1.509529,
    },
    "Ganges-Brahmaputra": {
        "mean": -3.943874,
        "median": -3.262842,
        "area_weighted": -3.602161,
    },
    "Mississippi": {
        "mean": -3.330013,
        "median": -3.196811,
        "area_weighted": -3.252499,
    },
}


@dataclass(frozen=True)
class SiteConfig:
    site_name: str
    model_tag: str
    cama_region_tag: str
    vlm_tif: str
    flow_stats_tag: str
    paper_delta_name: str


SITE_CONFIGS: Dict[str, SiteConfig] = {
    "PRD": SiteConfig("PRD", "PRD", "PRD", "pearl_vlm.tif", "PRD", "Pearl"),
    "YRD": SiteConfig("YRD", "YRD", "YRD", "yangtze_vlm.tif", "YRD", "Yangtze"),
    "BoB": SiteConfig("BoB", "BoB", "BoB", "ganges_vlm.tif", "BoB", "Ganges-Brahmaputra"),
    "Mississippi": SiteConfig(
        "Mississippi",
        "Misp",
        "Misp",
        "mississippi_vlm.tif",
        "MRD",
        "Mississippi",
    ),
}

SITE_ALIASES: Dict[str, str] = {
    "prd": "PRD",
    "pearl": "PRD",
    "yrd": "YRD",
    "yangtze": "YRD",
    "bob": "BoB",
    "ganges": "BoB",
    "misp": "Mississippi",
    "mississippi": "Mississippi",
}


def log(msg: str):
    print(msg, flush=True)


def normalize_site_name(site: str) -> str:
    key = str(site).strip()
    if key in SITE_CONFIGS:
        return key
    alias = SITE_ALIASES.get(key.lower())
    if alias:
        return alias
    raise KeyError(
        f"Unknown site '{site}'. Known sites: {list(SITE_CONFIGS)}; "
        f"aliases: {sorted(SITE_ALIASES)}"
    )


def normalize_lon_to_180(lon_vals):
    lon_vals = np.asarray(lon_vals, dtype=np.float64)
    return ((lon_vals + 180.0) % 360.0) - 180.0


def safe_adcirc_case_name(site: str, scenario: str, rp: str) -> str:
    return f"{site}__{scenario}__{rp}"


def flow_case_tag() -> str:
    if FLOW_FORCING_MODE in ("stats", "stats_p90"):
        return f"flowq{normalize_flow_percentile(FLOW_STATS_PERCENTILE)}"
    if FLOW_FORCING_MODE == "cama":
        return "flowcama"
    if FLOW_FORCING_MODE == "none":
        return "flownone"
    return f"flow{FLOW_FORCING_MODE}"


def safe_output_case_name(site: str, scenario: str, rp: str) -> str:
    return f"{site}__{scenario}_{flow_case_tag()}_{rp}"


def make_output_root(base_model_root: Path, case_name: str) -> str:
    if OUTPUT_PARENT_ROOT:
        return str(Path(OUTPUT_PARENT_ROOT) / case_name)
    base_name = base_model_root.name
    if base_name.endswith("_SFINCS"):
        site_tag = base_name[: -len("_SFINCS")]
    else:
        site_tag = base_name
    return str(base_model_root.parent / f"{site_tag}_future" / case_name)


def case_iter(
    sites: Sequence[str],
    scenarios: Sequence[str],
    return_periods: Sequence[str],
) -> Iterable[Tuple[SiteConfig, str, str, str, str, Path]]:
    for site in sites:
        site = normalize_site_name(site)
        cfg = SITE_CONFIGS[site]
        for scenario in scenarios:
            for rp in return_periods:
                output_case_name = safe_output_case_name(cfg.site_name, scenario, rp)
                adcirc_case_name = safe_adcirc_case_name(cfg.site_name, scenario, rp)
                yield (
                    cfg,
                    scenario,
                    rp,
                    output_case_name,
                    adcirc_case_name,
                    FUTURE_ADCIRC_ROOT / adcirc_case_name,
                )


# =============================================================================
# 1. fort22_meta track parsing and R-CLIPER rain
# =============================================================================

def parse_fort22_track_points(meta_path: Path) -> pd.DataFrame:
    lines = meta_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip() == "[ForcingTrackPoints]":
            start = i + 1
            break
    if start is None:
        raise ValueError(f"{meta_path} has no [ForcingTrackPoints] table")

    table_lines = []
    for line in lines[start:]:
        s = line.strip()
        if not s:
            if table_lines:
                break
            continue
        if s.startswith("[") and table_lines:
            break
        table_lines.append(s)

    if len(table_lines) < 2:
        raise ValueError(f"{meta_path} has an empty [ForcingTrackPoints] table")

    reader = csv.DictReader(table_lines)
    rows = list(reader)
    if not rows:
        raise ValueError(f"{meta_path} has no forcing track rows")

    df = pd.DataFrame(rows)
    rename = {}
    for col in df.columns:
        lc = col.strip().lower()
        if lc in {"time", "datetime", "date"}:
            rename[col] = "time"
        elif lc in {"lon180", "lon", "longitude"}:
            rename[col] = "lon180"
        elif lc in {"lat", "latitude"}:
            rename[col] = "lat"
        elif lc in {"vmax_ms", "vmax", "wind_ms", "max_wind_ms"}:
            rename[col] = "vmax_ms"
        elif lc in {"pressure_hpa", "pc_hpa", "pressure"}:
            rename[col] = "pressure_hPa"
    df = df.rename(columns=rename)

    required = {"time", "lon180", "lat", "vmax_ms"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{meta_path} forcing track table missing columns: {sorted(missing)}")

    df["time"] = pd.to_datetime(df["time"])
    for col in ["lon180", "lat", "vmax_ms", "pressure_hPa"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = (
        df.dropna(subset=["time", "lon180", "lat", "vmax_ms"])
        .sort_values("time")
        .drop_duplicates(subset=["time"], keep="last")
        .reset_index(drop=True)
    )
    if df.empty:
        raise ValueError(f"{meta_path} has no valid forcing track rows")
    return df


def interpolate_track_to_times(track: pd.DataFrame, target_time: pd.DatetimeIndex) -> pd.DataFrame:
    track = track.copy().set_index("time").sort_index()
    union_time = track.index.union(target_time).sort_values()

    out = pd.DataFrame(index=target_time)
    lon_rad = np.unwrap(np.deg2rad(track["lon180"].astype(float).values))
    lon_interp = (
        pd.Series(lon_rad, index=track.index)
        .reindex(union_time)
        .interpolate("time")
        .ffill()
        .bfill()
        .reindex(target_time)
    )
    out["lon180"] = normalize_lon_to_180(np.rad2deg(lon_interp.values))

    for col in ["lat", "vmax_ms", "pressure_hPa"]:
        if col not in track.columns:
            continue
        out[col] = (
            track[col].astype(float)
            .reindex(union_time)
            .interpolate("time")
            .ffill()
            .bfill()
            .reindex(target_time)
            .values
        )

    out["time"] = target_time
    return out.reset_index(drop=True)


def rcliper_rainrate_mmhr(
    vmax_ms: float,
    center_lon: float,
    center_lat: float,
    grid_lon: np.ndarray,
    grid_lat: np.ndarray,
    max_dist_km: float = RAIN_MAX_DIST_KM,
) -> np.ndarray:
    """R-CLIPER rain rate in mm/h after CLIMADA petals' Tuleya et al. setup."""
    rain = np.zeros(grid_lon.shape, dtype=np.float32)
    if not np.isfinite(vmax_ms) or vmax_ms <= 0:
        return rain

    dlon = normalize_lon_to_180(grid_lon - center_lon)
    dx_km = dlon * np.cos(np.deg2rad(grid_lat)) * 111.12
    dy_km = (grid_lat - center_lat) * 111.12
    dist_km = np.sqrt(dx_km * dx_km + dy_km * dy_km)
    close = dist_km <= max_dist_km
    if not np.any(close):
        return rain

    kn_to_ms = 0.514444
    inch_to_mm = 25.4
    hours_per_day = 24.0

    a1, a2, a3, a4 = -1.10, -1.60, 64.5, 150.0
    b1, b2, b3, b4 = 3.96, 4.80, -13.0, -16.0

    u_norm = 1.0 + (float(vmax_ms) / kn_to_ms - 35.0) / 33.0
    rainr_0 = a1 + b1 * u_norm
    rainr_m = a2 + b2 * u_norm
    rad_m = max(a3 + b3 * u_norm, 1.0)
    rad_e = max(a4 + b4 * u_norm, 1.0)

    d = dist_km[close]
    rain_close = np.zeros(d.shape, dtype=np.float64)
    inner = d <= rad_m
    rain_close[inner] = rainr_0 + (rainr_m - rainr_0) * (d[inner] / rad_m)
    rain_close[~inner] = rainr_m * np.exp(-(d[~inner] - rad_m) / rad_e)
    rain_close *= inch_to_mm / hours_per_day
    rain_close[~np.isfinite(rain_close)] = 0.0
    rain_close[rain_close < 0.0] = 0.0

    rain[close] = rain_close.astype(np.float32)
    return rain


def build_tc_rain_dataset(
    track: pd.DataFrame,
    target_time: pd.DatetimeIndex,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    model_lon2d: np.ndarray,
    model_lat2d: np.ndarray,
    model_crs,
    case_name: str,
) -> xr.Dataset:
    track_i = interpolate_track_to_times(track, target_time)
    rain_stack = np.zeros(
        (len(target_time), len(y_grid), len(x_grid)),
        dtype=np.float32,
    )

    for it, row in track_i.iterrows():
        rain_stack[it, :, :] = rcliper_rainrate_mmhr(
            vmax_ms=float(row["vmax_ms"]),
            center_lon=float(row["lon180"]),
            center_lat=float(row["lat"]),
            grid_lon=model_lon2d,
            grid_lat=model_lat2d,
        )

    pr = xr.DataArray(
        rain_stack,
        dims=("time", "y", "x"),
        coords={
            "time": target_time,
            "y": y_grid.astype(np.float64),
            "x": x_grid.astype(np.float64),
        },
        name="precip",
        attrs={
            "units": "mm/hr",
            "model": "R-CLIPER",
            "source": "fort22_meta.txt [ForcingTrackPoints]",
        },
    )
    ds = xr.Dataset(
        {"precip": pr},
        attrs={
            "title": f"R-CLIPER tropical cyclone rain forcing for {case_name}",
            "history": (
                "Generated from fort22_meta.txt track points using the "
                "Tuleya et al. R-CLIPER parameterization as implemented in "
                "CLIMADA petals tc_rainfield."
            ),
        },
    )
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    ds = ds.rio.write_crs(model_crs, inplace=False)
    ds["time"].encoding["calendar"] = "proleptic_gregorian"
    ds["time"].encoding["units"] = (
        f"hours since {pd.Timestamp(target_time[0]).strftime('%Y-%m-%d %H:%M:%S')}"
    )
    return ds


def prepare_tc_rain_for_compression(ds: xr.Dataset) -> xr.Dataset:
    """
    Keep rain as SFINCS-readable float32, but reduce entropy before compression.

    The final NetCDF variable remains a normal float32 array. Tiny near-zero
    rain rates are set to zero and values are rounded to a small number of
    decimal places. This makes zlib/shuffle compression much more effective
    without relying on non-portable half-precision NetCDF storage.
    """
    attrs = dict(ds["precip"].attrs)
    pr = ds["precip"].astype(np.float32)
    if TC_RAIN_ZERO_BELOW_MMHR is not None and float(TC_RAIN_ZERO_BELOW_MMHR) > 0.0:
        pr = pr.where(pr >= float(TC_RAIN_ZERO_BELOW_MMHR), 0.0)
    if TC_RAIN_LEAST_SIGNIFICANT_DIGIT is not None:
        decimals = int(TC_RAIN_LEAST_SIGNIFICANT_DIGIT)
        pr = pr.round(decimals).astype(np.float32)
        attrs["least_significant_digit"] = decimals
        attrs["value_precision"] = "rounded to {} decimal places in mm/hr, stored as float32".format(decimals)
    ds["precip"] = pr
    ds["precip"].attrs.update(attrs)
    return ds


def write_tc_rain_netcdf_compressed(ds: xr.Dataset, out_nc: Path):
    """Write SFINCS rain forcing as compressed, standard float32 NetCDF4."""
    encoding = {}
    for name, da in ds.data_vars.items():
        if np.issubdtype(da.dtype, np.floating):
            chunks = []
            for dim in da.dims:
                dim_len = int(da.sizes[dim])
                if dim == "time":
                    chunks.append(1)
                else:
                    chunks.append(min(512, dim_len))
            encoding[name] = {
                "zlib": True,
                "complevel": int(TC_RAIN_COMPRESSION_LEVEL),
                "shuffle": True,
                "dtype": "float32",
                "chunksizes": tuple(chunks),
            }
            if name == "precip" and TC_RAIN_LEAST_SIGNIFICANT_DIGIT is not None:
                encoding[name]["least_significant_digit"] = int(TC_RAIN_LEAST_SIGNIFICANT_DIGIT)
        else:
            encoding[name] = {
                "zlib": True,
                "complevel": int(TC_RAIN_COMPRESSION_LEVEL),
                "shuffle": True,
            }

    if "time" in ds.coords:
        encoding["time"] = {
            "zlib": True,
            "complevel": int(TC_RAIN_COMPRESSION_LEVEL),
            "shuffle": True,
        }

    ds.to_netcdf(out_nc, engine="netcdf4", encoding=encoding)


# =============================================================================
# 2. ADCIRC water-level source helpers
# =============================================================================

def find_adcirc_elevation_source(case_dir: Path, case_name: str) -> Tuple[Path, str]:
    fort63_candidates = [
        case_dir / "fort.63",
        case_dir / "fort.63.nc",
        ADCIRC_RESULT_ROOT / case_name / "fort.63",
        ADCIRC_RESULT_ROOT / case_name / "fort.63.nc",
    ]
    for path in fort63_candidates:
        if path.exists():
            return path, "fort63"

    maxele_candidates = [
        case_dir / "maxele.63",
        ADCIRC_RESULT_ROOT / case_name / "maxele.63",
    ]
    for path in maxele_candidates:
        if path.exists() and ALLOW_STATIC_MAXELE63_BZS:
            return path, "maxele63_static"

    checked = fort63_candidates + maxele_candidates
    raise FileNotFoundError(
        "No ADCIRC elevation source found. Checked: "
        + "; ".join(str(p) for p in checked)
    )


def read_maxele63_first_dataset(maxele_path: Path) -> pd.DataFrame:
    with maxele_path.open("r", encoding="utf-8", errors="ignore") as f:
        _ = f.readline()
        hdr = f.readline().split()
        if len(hdr) < 2:
            raise ValueError(f"Bad maxele.63 header in {maxele_path}")
        nnodes = int(float(hdr[1]))
        _ = f.readline()

        node_ids = np.zeros(nnodes, dtype=np.int64)
        values = np.full(nnodes, np.nan, dtype=np.float64)
        for inode in range(1, nnodes + 1):
            vals = f.readline().split()
            if not vals:
                continue
            if len(vals) >= 2:
                try:
                    node_ids[inode - 1] = int(float(vals[0]))
                    values[inode - 1] = float(vals[1])
                except Exception:
                    node_ids[inode - 1] = inode
                    values[inode - 1] = np.nan
            else:
                node_ids[inode - 1] = inode
                values[inode - 1] = float(vals[0])

    values[np.isin(values, [-99999.0, -9999.0, -1.0e30])] = np.nan
    return pd.DataFrame({"node_id": node_ids, "maxele_m": values})


def match_bnd_to_static_maxele(
    cand: pd.DataFrame,
    bnd_lon: np.ndarray,
    bnd_lat: np.ndarray,
    maxele_path: Path,
    target_time: pd.DatetimeIndex,
    k_nearest: int = 50,
):
    maxele = read_maxele63_first_dataset(maxele_path)
    cand2 = cand.merge(maxele, on="node_id", how="left")
    valid = cand2[np.isfinite(cand2["maxele_m"].values)].copy()
    if valid.empty:
        raise RuntimeError(f"No valid maxele values in {maxele_path}")

    tree = cKDTree(np.c_[cand2["lon"].values, cand2["lat"].values])
    valid_tree = cKDTree(np.c_[valid["lon"].values, valid["lat"].values])
    k_use = min(k_nearest, len(cand2))
    dist_k, idx_k = tree.query(np.c_[bnd_lon, bnd_lat], k=k_use)
    if k_use == 1:
        dist_k = dist_k[:, None]
        idx_k = idx_k[:, None]

    chosen_ids = []
    chosen_lons = []
    chosen_lats = []
    chosen_dists = []
    chosen_vals = []

    for i in range(len(bnd_lon)):
        best = None
        for j in range(idx_k.shape[1]):
            rr = cand2.iloc[int(idx_k[i, j])]
            val = float(rr["maxele_m"]) if np.isfinite(rr["maxele_m"]) else np.nan
            if np.isfinite(val):
                best = (
                    int(rr["node_id"]),
                    float(rr["lon"]),
                    float(rr["lat"]),
                    float(dist_k[i, j]),
                    val,
                )
                break
        if best is None:
            d_one, idx_one = valid_tree.query([[bnd_lon[i], bnd_lat[i]]], k=1)
            rr = valid.iloc[int(idx_one[0])]
            best = (
                int(rr["node_id"]),
                float(rr["lon"]),
                float(rr["lat"]),
                float(d_one[0]),
                float(rr["maxele_m"]),
            )
        chosen_ids.append(best[0])
        chosen_lons.append(best[1])
        chosen_lats.append(best[2])
        chosen_dists.append(best[3])
        chosen_vals.append(best[4])

    arr = np.tile(np.asarray(chosen_vals, dtype=float), (len(target_time), 1))
    wl_df = pd.DataFrame(arr, index=target_time, columns=np.arange(1, len(chosen_vals) + 1))
    return (
        wl_df,
        chosen_ids,
        np.asarray(chosen_lons),
        np.asarray(chosen_lats),
        np.asarray(chosen_dists),
    )


def build_waterlevel_forcing(
    case_name: str,
    case_dir: Path,
    fort14_path: Path,
    bnd_lon: np.ndarray,
    bnd_lat: np.ndarray,
    tstart: pd.Timestamp,
    tstop: pd.Timestamp,
    dt_adcirc: int,
    target_time: pd.DatetimeIndex,
):
    node_df, _open_bnd_nodes = base.parse_fort14_nodes_and_open_boundary(str(fort14_path))
    cand = node_df.copy()
    elev_path, source_kind = find_adcirc_elevation_source(case_dir, case_name)
    log(f"   ADCIRC elevation source: {elev_path} ({source_kind})")

    if source_kind == "fort63":
        wl_df, matched_nodes, matched_lons, matched_lats, matched_dists = (
            base.match_bnd_to_valid_adcirc_nodes(
                cand=cand,
                bnd_lon=bnd_lon,
                bnd_lat=bnd_lat,
                fort63_path=str(elev_path),
                start_time=tstart,
                dt_seconds=dt_adcirc,
                k_nearest=20,
            )
        )
        wl_df.index = pd.to_datetime(wl_df.index)
        wl_df = (
            wl_df.reindex(wl_df.index.union(target_time))
            .sort_index()
            .interpolate("time")
            .reindex(target_time)
            .interpolate(limit_direction="both")
            .ffill()
            .bfill()
            .fillna(0.0)
        )
    else:
        log(
            "   [WARNING] Using maxele.63 as a constant boundary. "
            "Use fort.63 when a transient downstream water-level boundary is available."
        )
        wl_df, matched_nodes, matched_lons, matched_lats, matched_dists = (
            match_bnd_to_static_maxele(
                cand=cand,
                bnd_lon=bnd_lon,
                bnd_lat=bnd_lat,
                maxele_path=elev_path,
                target_time=target_time,
            )
        )

    return wl_df, matched_nodes, matched_lons, matched_lats, matched_dists, elev_path, source_kind


# =============================================================================
# 3. VLM subsidence
# =============================================================================

def should_apply_vlm(scenario: str) -> bool:
    if not APPLY_VLM_SUBSIDENCE:
        return False
    if scenario == "historical_era5" and not APPLY_VLM_TO_HISTORICAL:
        return False
    return True


def read_site_subsidence_from_paper_table(cfg: SiteConfig) -> Tuple[Path, float, str, float]:
    delta_name = cfg.paper_delta_name
    if delta_name not in PAPER_VLM_RATES_MM_YR:
        raise KeyError(
            f"No paper-table VLM rate configured for {cfg.site_name} "
            f"(delta name: {delta_name})"
        )
    stat = PAPER_VLM_RATE_STAT.lower()
    if stat not in PAPER_VLM_RATES_MM_YR[delta_name]:
        raise KeyError(
            f"Unsupported PAPER_VLM_RATE_STAT={PAPER_VLM_RATE_STAT}. "
            "Use 'mean', 'median', or 'area_weighted'."
        )
    rate_mm_yr = float(PAPER_VLM_RATES_MM_YR[delta_name][stat])
    subsidence_m = abs(rate_mm_yr) * VLM_YEARS / 1000.0
    return PAPER_VLM_TABLE_PATH, rate_mm_yr, PAPER_VLM_RATE_UNIT, subsidence_m


def read_site_subsidence_from_grid_tif(cfg: SiteConfig) -> Tuple[Path, float, str, float]:
    vlm_path = VLM_ROOT / cfg.vlm_tif
    if not vlm_path.exists():
        raise FileNotFoundError(f"Cannot find VLM raster: {vlm_path}")

    da = rioxarray.open_rasterio(vlm_path, masked=True).squeeze(drop=True)
    arr = np.asarray(da.values, dtype=np.float64)
    valid = np.isfinite(arr)
    if not np.any(valid):
        raise ValueError(f"No valid VLM values in {vlm_path}")

    mean_rate = float(np.nanmean(arr[valid]))
    # VLM rasters commonly use negative values for downward land motion.
    # Here the user request is to lower the whole SFINCS domain by the
    # subsidence magnitude, so keep the applied amount positive.
    if VLM_RATE_UNIT.lower() in ("cm/yr", "cm/year"):
        subsidence_m = abs(mean_rate) * VLM_YEARS / 100.0
    elif VLM_RATE_UNIT.lower() in ("mm/yr", "mm/year"):
        subsidence_m = abs(mean_rate) * VLM_YEARS / 1000.0
    else:
        raise ValueError(f"Unsupported VLM_RATE_UNIT: {VLM_RATE_UNIT}")
    return vlm_path, mean_rate, VLM_RATE_UNIT, subsidence_m


def read_site_subsidence(cfg: SiteConfig) -> Tuple[Path, float, str, float]:
    source = VLM_SOURCE.lower()
    if source == "paper_table":
        return read_site_subsidence_from_paper_table(cfg)
    if source == "grid_tif":
        return read_site_subsidence_from_grid_tif(cfg)
    raise ValueError(f"Unsupported VLM_SOURCE: {VLM_SOURCE}")


def apply_uniform_subsidence_to_dep(
    model_root: Path,
    subsidence_m: float,
    vlm_path: Path,
    vlm_rate: float,
    vlm_rate_unit: str,
    scenario: str,
    case_name: str,
):
    dep_path = model_root / "sfincs.dep"
    marker_path = model_root / "vlm_subsidence_applied.txt"
    backup_path = model_root / "sfincs.dep.before_vlm"

    if not dep_path.exists():
        raise FileNotFoundError(dep_path)
    if abs(subsidence_m) < 1.0e-9:
        log("   VLM mean is near zero; sfincs.dep unchanged.")
        return

    if marker_path.exists():
        txt = marker_path.read_text(encoding="utf-8", errors="ignore")
        old = None
        for line in txt.splitlines():
            if line.startswith("subsidence_m="):
                old = float(line.split("=", 1)[1])
                break
        if old is not None and abs(old - subsidence_m) < 1.0e-6:
            log("   VLM subsidence marker exists; sfincs.dep was already lowered for this amount.")
            return
        raise RuntimeError(
            f"{marker_path} exists with a different subsidence amount. "
            "Use --overwrite to rebuild this output model from the base model."
        )

    shutil.copy2(dep_path, backup_path)
    vals = np.fromfile(dep_path, dtype="<f4")
    vals = vals.astype(np.float32, copy=True)
    finite = np.isfinite(vals)
    vals[finite] = vals[finite] - np.float32(subsidence_m)
    vals.astype("<f4").tofile(dep_path)

    marker = "\n".join(
        [
            f"case_name={case_name}",
            f"scenario={scenario}",
            f"vlm_source={VLM_SOURCE}",
            f"vlm_source_path={vlm_path}",
            f"vlm_rate_stat={PAPER_VLM_RATE_STAT if VLM_SOURCE.lower() == 'paper_table' else 'grid_mean'}",
            f"vlm_rate={vlm_rate:.8f}",
            f"vlm_rate_unit={vlm_rate_unit}",
            f"years={VLM_YEARS:.3f}",
            f"subsidence_m={subsidence_m:.8f}",
            "operation=sfincs.dep values minus subsidence_m",
            f"backup={backup_path.name}",
        ]
    )
    marker_path.write_text(marker + "\n", encoding="utf-8")
    log(
        f"   VLM applied: source={VLM_SOURCE}, rate={vlm_rate:.4f} {vlm_rate_unit}, "
        f"years={VLM_YEARS:g}, lowering={subsidence_m:.4f} m"
    )


# =============================================================================
# 4. Flow statistics discharge
# =============================================================================

def scenario_to_flow_key(scenario: str) -> str:
    if scenario in FLOW_STATS_SCENARIO_MAP:
        return FLOW_STATS_SCENARIO_MAP[scenario]
    if scenario.startswith("ssp"):
        return scenario.split("_", 1)[0]
    if scenario.startswith("historical"):
        return "historical_gcm"
    return scenario


def flow_stats_path(cfg: SiteConfig) -> Path:
    return FLOW_STATS_ROOT / FLOW_STATS_FILENAME_TEMPLATE.format(site=cfg.flow_stats_tag)


def normalize_flow_percentile(value) -> int:
    text = str(value).strip().upper()
    if text.startswith("P"):
        text = text[1:]
    if text.endswith("%"):
        text = text[:-1]
    percentile = float(text)
    if not percentile.is_integer():
        raise ValueError(
            f"Flow statistics percentile must be an integer for allocatedP* columns, got {value!r}"
        )
    percentile = int(percentile)
    if percentile < 0 or percentile > 100:
        raise ValueError(f"Flow statistics percentile must be between 0 and 100, got {value!r}")
    return percentile


def flow_stats_quantile_suffix() -> str:
    return f"P{normalize_flow_percentile(FLOW_STATS_PERCENTILE)}"


def parse_src_cols(value) -> List[int]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    out = []
    for part in text.replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        out.append(int(float(part)))
    return out


def add_flow_stats_discharge(
    sf,
    cfg: SiteConfig,
    scenario: str,
    output_model_root: Path,
    src_xy: np.ndarray,
    model_crs,
    target_time: pd.DatetimeIndex,
):
    if src_xy.shape[0] == 0:
        log("   No sfincs.src found; flow-statistics discharge skipped.")
        return None, None

    stats_path = flow_stats_path(cfg)
    if not stats_path.exists():
        raise FileNotFoundError(f"Cannot find flow statistics CSV: {stats_path}")

    flow_key = scenario_to_flow_key(scenario)
    flow_quantile = flow_stats_quantile_suffix()
    value_col = f"allocated{flow_quantile}"
    df = pd.read_csv(stats_path)
    if "scenarioKey" not in df.columns:
        raise ValueError(f"{stats_path} is missing scenarioKey")
    if value_col not in df.columns:
        raise ValueError(f"{stats_path} is missing {value_col}")
    if "sfincsSrcCols" not in df.columns:
        raise ValueError(f"{stats_path} is missing sfincsSrcCols")

    rows = df[df["scenarioKey"].astype(str).str.lower() == flow_key.lower()].copy()
    if rows.empty:
        raise ValueError(
            f"No rows for scenarioKey={flow_key} in {stats_path}. "
            f"Available keys: {sorted(df['scenarioKey'].astype(str).unique())}"
        )

    q_per_src = np.zeros(src_xy.shape[0], dtype=np.float64)
    diag_rows = []
    for _, row in rows.iterrows():
        cols = parse_src_cols(row["sfincsSrcCols"])
        if not cols:
            continue
        q_total = float(row[value_col])
        q_each = q_total / float(len(cols))
        for col in cols:
            if col < 1 or col > len(q_per_src):
                raise IndexError(
                    f"{stats_path}: sfincsSrcCol {col} is outside 1..{len(q_per_src)}"
                )
            q_per_src[col - 1] += q_each
        diag_rows.append({
            "scenario": scenario,
            "scenarioKey": flow_key,
            "inletId": row.get("inletId", np.nan),
            "inletLabel": row.get("inletLabel", ""),
            "boundarySide": row.get("boundarySide", ""),
            "sfincsSrcCols": ",".join(str(c) for c in cols),
            "sfincsSrcNpts": len(cols),
            "flow_stat_column": value_col,
            "allocated_total_m3s": q_total,
            "assigned_each_src_m3s": q_each,
        })

    q = np.tile(q_per_src.reshape(1, -1), (len(target_time), 1))
    q_df = pd.DataFrame(q, index=target_time, columns=np.arange(1, len(q_per_src) + 1))

    src_gdf = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(src_xy[:, 0], src_xy[:, 1]),
        crs=model_crs,
    )
    src_gdf.index = np.arange(1, len(src_gdf) + 1)
    sf.setup_discharge_forcing(timeseries=q_df, locations=src_gdf)

    out_csv = output_model_root / "flow_statistics_src_discharge.csv"
    per_src = pd.DataFrame({
        "sfincs_src_id": np.arange(1, len(q_per_src) + 1),
        "sfincs_x": src_xy[:, 0],
        "sfincs_y": src_xy[:, 1],
        "flow_m3s": q_per_src,
        "source_csv": str(stats_path),
        "scenarioKey": flow_key,
        "flow_stat_column": value_col,
        "flow_stat_percentile": normalize_flow_percentile(FLOW_STATS_PERCENTILE),
    })
    per_src.to_csv(out_csv, index=False, encoding="utf-8-sig")

    inlet_csv = output_model_root / "flow_statistics_inlet_allocation.csv"
    pd.DataFrame(diag_rows).to_csv(inlet_csv, index=False, encoding="utf-8-sig")

    log(
        f"   Flow statistics discharge: {stats_path.name}, scenarioKey={flow_key}, "
        f"{value_col}, total={float(np.sum(q_per_src)):.3f} m3/s."
    )
    return q_df, out_csv


# =============================================================================
# 5. CaMa discharge
# =============================================================================

def resolve_cama_geojson(
    output_model_root: Path,
    cama_run_script: str,
    sh_cfg: Dict[str, str],
) -> str:
    candidates = [
        sh_cfg.get("GEOJSON_PATH"),
        str(output_model_root / "cama_upstream_rivers.geojson"),
    ]
    for cand in candidates:
        if cand and os.path.exists(cand):
            return cand
    raise FileNotFoundError(
        "Cannot find cama_upstream_rivers.geojson. Checked: "
        + "; ".join(str(c) for c in candidates if c)
        + f"; CAMA script: {cama_run_script}"
    )


def add_cama_discharge(
    sf,
    cfg: SiteConfig,
    output_model_root: Path,
    src_xy: np.ndarray,
    src_lon: np.ndarray,
    src_lat: np.ndarray,
    model_crs,
    target_time: pd.DatetimeIndex,
    event_year: int,
):
    if src_xy.shape[0] == 0:
        log("   No sfincs.src found; CaMa discharge skipped.")
        return None, None
    if not ENABLE_CAMA_DISCHARGE:
        log("   ENABLE_CAMA_DISCHARGE=False; CaMa discharge skipped.")
        return None, None

    cama_run_script = CAMA_RUN_SCRIPT_TEMPLATE.format(region_tag=cfg.cama_region_tag)
    sh_cfg = base.parse_shell_assignments(cama_run_script) if os.path.exists(cama_run_script) else {}
    cama_res = sh_cfg.get("CAMA_RESOLUTION", "15min")

    try:
        cama_geojson = resolve_cama_geojson(output_model_root, cama_run_script, sh_cfg)
        cama_outflow_bin = base.resolve_cama_outflow_path(
            cama_run_script,
            CAMA_OUTFLOW_BIN_TEMPLATE,
            event_year=event_year,
            region_tag=cfg.cama_region_tag,
        )
        if not os.path.exists(cama_outflow_bin):
            raise FileNotFoundError(f"Cannot find CaMa outflw file: {cama_outflow_bin}")
    except Exception:
        if SKIP_CAMA_IF_MISSING:
            log("   [WARNING] CaMa input missing; discharge skipped.")
            return None, None
        raise

    log(f"   CaMa outflw file: {cama_outflow_bin}")
    cama_map_dir = base.resolve_cama_map_dir(
        cama_run_script,
        cama_outflow_bin=cama_outflow_bin,
        cama_res=cama_res,
        region_tag=cfg.cama_region_tag,
        explicit_map_dir=CAMA_MAP_DIR,
    )
    cama_grid = base.read_cama_map_params(cama_map_dir, cama_res=cama_res)

    west = float(cama_grid["west"])
    east = float(cama_grid["east"])
    south = float(cama_grid["south"])
    north = float(cama_grid["north"])
    cama_res_deg = float(cama_grid["res_deg"])
    nx_c = int(cama_grid["nx"])
    ny_c = int(cama_grid["ny"])

    cama_lons = west + (np.arange(nx_c) + 0.5) * cama_res_deg
    cama_lats = north - (np.arange(ny_c) + 0.5) * cama_res_deg
    log(
        "   CaMa grid: "
        f"nx={nx_c}, ny={ny_c}, res={cama_res_deg:.6f}, "
        f"bounds=({west:.3f}, {south:.3f}, {east:.3f}, {north:.3f})"
    )

    src_cama_mapping = base.match_sfincs_src_to_cama_boundaries(
        cama_geojson,
        src_lon=src_lon,
        src_lat=src_lat,
    )
    matched_cama_lonlat = src_cama_mapping[["cama_lon", "cama_lat"]].astype(float).values
    dist_src = src_cama_mapping["match_dist_src_deg"].astype(float).values

    outflw = base.read_cama_outflow_binary(
        cama_outflow_bin,
        ny=ny_c,
        nx=nx_c,
        dtype=CAMA_DTYPE,
        array_order=CAMA_ARRAY_ORDER,
    )
    nt_c = outflw.shape[0]
    cama_time_start = (
        pd.Timestamp(CAMA_TIME_START)
        if CAMA_TIME_START is not None
        else pd.Timestamp(f"{event_year}-01-01 00:00:00")
    )
    cama_time = pd.date_range(cama_time_start, periods=nt_c, freq=f"{int(CAMA_OUTPUT_DT_HOURS)}H")

    lon2, lat2 = np.meshgrid(cama_lons, cama_lats)
    tree_cama = cKDTree(np.c_[lon2.ravel(), lat2.ravel()])
    dist_cell, idx_cell = tree_cama.query(matched_cama_lonlat, k=1)
    iy = idx_cell // nx_c
    ix = idx_cell % nx_c

    q = np.full((nt_c, len(src_xy)), np.nan, dtype=float)
    cell_keys = [(int(iy[j]), int(ix[j])) for j in range(len(src_xy))]
    cell_counts = {}
    for key in cell_keys:
        cell_counts[key] = cell_counts.get(key, 0) + 1
    split_counts = np.array([cell_counts[key] for key in cell_keys], dtype=int)
    split_fractions = 1.0 / split_counts.astype(float)

    for j in range(len(src_xy)):
        q_raw = outflw[:, iy[j], ix[j]]
        if SPLIT_CAMA_DISCHARGE_FOR_DUPLICATE_SRC:
            q[:, j] = q_raw / float(split_counts[j])
        else:
            q[:, j] = q_raw

    q_df = pd.DataFrame(q, index=cama_time, columns=np.arange(1, len(src_xy) + 1))
    q_df = q_df.interpolate(limit_direction="both")
    q_df = (
        q_df.reindex(q_df.index.union(target_time))
        .sort_index()
        .interpolate("time")
        .reindex(target_time)
        .interpolate(limit_direction="both")
        .ffill()
        .bfill()
        .fillna(0.0)
    )

    src_gdf = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(src_xy[:, 0], src_xy[:, 1]),
        crs=model_crs,
    )
    src_gdf.index = np.arange(1, len(src_gdf) + 1)
    sf.setup_discharge_forcing(timeseries=q_df, locations=src_gdf)

    out_src_match_csv = output_model_root / "cama_src_match.csv"
    df_src_match = pd.DataFrame({
        "sfincs_src_id": np.arange(1, len(src_xy) + 1),
        "sfincs_x": src_xy[:, 0],
        "sfincs_y": src_xy[:, 1],
        "sfincs_lon": src_lon,
        "sfincs_lat": src_lat,
        "mapping_method": src_cama_mapping["mapping_method"].values,
        "cama_boundary_id": src_cama_mapping["cama_boundary_id"].values,
        "inlet_id": src_cama_mapping["inlet_id"].values,
        "geojson_cama_row": src_cama_mapping["cama_row"].values,
        "geojson_cama_col": src_cama_mapping["cama_col"].values,
        "geojson_uparea_km2": src_cama_mapping["uparea_km2"].values,
        "snap_start_lon": src_cama_mapping["source_lon"].values,
        "snap_start_lat": src_cama_mapping["source_lat"].values,
        "cama_lon": matched_cama_lonlat[:, 0],
        "cama_lat": matched_cama_lonlat[:, 1],
        "cama_cell_lon": lon2.ravel()[idx_cell],
        "cama_cell_lat": lat2.ravel()[idx_cell],
        "cama_cell_iy": iy,
        "cama_cell_ix": ix,
        "shared_src_count_for_same_cama_cell": split_counts,
        "discharge_fraction_assigned": split_fractions,
        "match_dist_snap_deg": dist_src,
        "match_dist_cama_cell_deg": dist_cell,
    })
    df_src_match.to_csv(out_src_match_csv, index=False, encoding="utf-8-sig")

    if np.nanmax(dist_cell) > CAMA_MATCH_MAX_DIST_DEG:
        log(f"   [WARNING] CaMa match distance exceeds {CAMA_MATCH_MAX_DIST_DEG} deg.")

    log(
        f"   CaMa discharge written for {len(src_gdf)} source points "
        f"(min={float(np.nanmin(q_df.values)):.3f}, max={float(np.nanmax(q_df.values)):.3f})."
    )
    return q_df, out_src_match_csv


# =============================================================================
# 5. One-case workflow
# =============================================================================

def process_one_case(
    cfg: SiteConfig,
    scenario: str,
    rp: str,
    case_name: str,
    adcirc_case_name: str,
    case_dir: Path,
):
    log("\n=======================================================")
    log(f"Case: {case_name}")
    log(f"ADCIRC source case: {adcirc_case_name}")
    log("=======================================================")

    if not case_dir.exists():
        raise FileNotFoundError(f"ADCIRC case directory not found: {case_dir}")

    meta_path = case_dir / "fort22_meta.txt"
    if not meta_path.exists():
        raise FileNotFoundError(meta_path)
    meta = base.parse_simple_key_value_file(str(meta_path))
    tstart, tstop, event_year = base.get_event_time_info(meta)
    dt_adcirc = int(float(meta.get("WTIMINC", 3600)))
    target_time = pd.date_range(tstart, tstop, freq=TARGET_FREQ)
    tc_id = meta.get("tc_id", case_name)

    base_model_root = EXAMPLES_ROOT / f"{cfg.model_tag}_SFINCS"
    output_model_root = Path(
        base.prepare_output_model_root(
            base_model_root=str(base_model_root),
            tc_id=tc_id,
            explicit_output_root=make_output_root(base_model_root, case_name),
            overwrite_existing=OVERWRITE_EXISTING_OUTPUT_MODEL,
        )
    )

    log(f"   event window: {tstart} -> {tstop}")
    log(f"   event year  : {event_year}")
    log(f"   output model: {output_model_root}")

    if should_apply_vlm(scenario):
        vlm_path, vlm_rate, vlm_rate_unit, subsidence_m = read_site_subsidence(cfg)
        apply_uniform_subsidence_to_dep(
            output_model_root,
            subsidence_m=subsidence_m,
            vlm_path=vlm_path,
            vlm_rate=vlm_rate,
            vlm_rate_unit=vlm_rate_unit,
            scenario=scenario,
            case_name=case_name,
        )
    else:
        log("   VLM subsidence skipped for this scenario.")

    sf = base.open_existing_sfincs_model(str(output_model_root))
    model_crs = base.get_model_crs(sf)
    x_grid, y_grid = base.get_model_xy(sf)
    transformer_to_ll = Transformer.from_crs(model_crs, "EPSG:4326", always_xy=True)

    bnd_path = output_model_root / "sfincs.bnd"
    src_path = output_model_root / "sfincs.src"
    inp_path = output_model_root / "sfincs.inp"
    if not bnd_path.exists():
        raise FileNotFoundError(bnd_path)

    sf.setup_config(
        tref=tstart.strftime("%Y%m%d %H%M%S"),
        tstart=tstart.strftime("%Y%m%d %H%M%S"),
        tstop=tstop.strftime("%Y%m%d %H%M%S"),
    )

    bnd_xy = base.read_xy_file(str(bnd_path))
    if bnd_xy.shape[0] == 0:
        raise RuntimeError(f"{bnd_path} is empty")
    bnd_lon, bnd_lat = transformer_to_ll.transform(bnd_xy[:, 0], bnd_xy[:, 1])

    if src_path.exists():
        src_xy = base.read_xy_file(str(src_path))
        if src_xy.shape[0] > 0:
            src_lon, src_lat = transformer_to_ll.transform(src_xy[:, 0], src_xy[:, 1])
        else:
            src_lon = np.array([])
            src_lat = np.array([])
    else:
        src_xy = np.empty((0, 2), dtype=float)
        src_lon = np.array([])
        src_lat = np.array([])

    fort14_path = case_dir / "fort.14"
    if not fort14_path.exists() and "fort14" in meta:
        fort14_path = Path(meta["fort14"])
    if not fort14_path.exists():
        raise FileNotFoundError(f"Cannot find fort.14 for ADCIRC case {adcirc_case_name}")

    log("   [1/3] ADCIRC downstream water level -> sfincs.bzs")
    (
        wl_df,
        matched_nodes,
        matched_lons,
        matched_lats,
        matched_dists,
        elev_path,
        elev_source_kind,
    ) = build_waterlevel_forcing(
        case_name=adcirc_case_name,
        case_dir=case_dir,
        fort14_path=fort14_path,
        bnd_lon=bnd_lon,
        bnd_lat=bnd_lat,
        tstart=tstart,
        tstop=tstop,
        dt_adcirc=dt_adcirc,
        target_time=target_time,
    )

    bzs_path = output_model_root / "sfincs.bzs"
    base.write_sfincs_bzs_file(str(bzs_path), wl_df)

    out_bnd_match_csv = output_model_root / "adcirc_bnd_match.csv"
    df_bnd_match = pd.DataFrame({
        "sfincs_bnd_id": np.arange(1, len(bnd_lon) + 1),
        "sfincs_x": bnd_xy[:, 0],
        "sfincs_y": bnd_xy[:, 1],
        "sfincs_lon": bnd_lon,
        "sfincs_lat": bnd_lat,
        "adcirc_node_id": matched_nodes,
        "adcirc_lon": matched_lons,
        "adcirc_lat": matched_lats,
        "match_dist_deg": matched_dists,
        "elevation_source": str(elev_path),
        "elevation_source_kind": elev_source_kind,
    })
    df_bnd_match.to_csv(out_bnd_match_csv, index=False, encoding="utf-8-sig")
    if np.nanmax(matched_dists) > ADCIRC_MATCH_MAX_DIST_DEG:
        log(f"   [WARNING] ADCIRC match distance exceeds {ADCIRC_MATCH_MAX_DIST_DEG} deg.")
    log(
        f"   bzs written: {bzs_path} "
        f"(min={float(np.nanmin(wl_df.values)):.4f}, max={float(np.nanmax(wl_df.values)):.4f})"
    )

    log("   [2/3] River discharge -> sfincs.dis")
    out_src_match_csv = None
    if FLOW_FORCING_MODE in ("stats", "stats_p90"):
        q_df_final, out_src_match_csv = add_flow_stats_discharge(
            sf=sf,
            cfg=cfg,
            scenario=scenario,
            output_model_root=output_model_root,
            src_xy=src_xy,
            model_crs=model_crs,
            target_time=target_time,
        )
    elif FLOW_FORCING_MODE == "cama":
        q_df_final, out_src_match_csv = add_cama_discharge(
            sf=sf,
            cfg=cfg,
            output_model_root=output_model_root,
            src_xy=src_xy,
            src_lon=src_lon,
            src_lat=src_lat,
            model_crs=model_crs,
            target_time=target_time,
            event_year=event_year,
        )
    elif FLOW_FORCING_MODE == "none":
        q_df_final = None
        log("   FLOW_FORCING_MODE='none'; river discharge skipped.")
    else:
        raise ValueError(f"Unknown FLOW_FORCING_MODE: {FLOW_FORCING_MODE}")

    log("   [3/3] fort22_meta TC track -> R-CLIPER rain netamprfile")
    x2d, y2d = np.meshgrid(x_grid, y_grid)
    model_lon2d, model_lat2d = transformer_to_ll.transform(x2d, y2d)
    track = parse_fort22_track_points(meta_path)
    ds_pr = build_tc_rain_dataset(
        track=track,
        target_time=target_time,
        x_grid=x_grid,
        y_grid=y_grid,
        model_lon2d=model_lon2d,
        model_lat2d=model_lat2d,
        model_crs=model_crs,
        case_name=case_name,
    )
    ds_pr = prepare_tc_rain_for_compression(ds_pr)

    tc_rain_src_nc = output_model_root / TC_RAIN_SOURCE_NC_NAME
    tc_rain_out_nc = output_model_root / TC_RAIN_OUT_NC_NAME
    if KEEP_TC_RAIN_SOURCE_NC:
        write_tc_rain_netcdf_compressed(ds_pr, tc_rain_src_nc)

    rain_temp_candidates = [
        output_model_root / "sfincs.nc",
        output_model_root / "sfincs_netampr.nc",
        output_model_root / "precip.nc",
        output_model_root / "precip_2d.nc",
    ]
    preexisting_temp_nc = {p.resolve() for p in rain_temp_candidates if p.exists()}

    sf.setup_precip_forcing_from_grid(
        precip=ds_pr["precip"],
        aggregate=False,
        cumulative_input=False,
    )
    sf.write_forcing()

    generated_netcdf = None
    for cand_fn in [
        "sfincs.nc",
        "sfincs_netampr.nc",
        tc_rain_out_nc.name,
        "precip.nc",
    ]:
        cand_path = output_model_root / cand_fn
        if cand_path.exists():
            generated_netcdf = cand_path
            break

    if generated_netcdf is not None and generated_netcdf.resolve() != tc_rain_out_nc.resolve():
        with xr.open_dataset(generated_netcdf, engine="netcdf4") as ds_tmp:
            write_tc_rain_netcdf_compressed(ds_tmp.load(), tc_rain_out_nc)
    elif generated_netcdf is None:
        write_tc_rain_netcdf_compressed(ds_pr, tc_rain_out_nc)

    if not KEEP_TC_RAIN_SOURCE_NC and tc_rain_src_nc.exists():
        tc_rain_src_nc.unlink()

    if CLEAN_TEMP_RAIN_NC:
        for temp_nc in rain_temp_candidates:
            if not temp_nc.exists():
                continue
            if temp_nc.resolve() == tc_rain_out_nc.resolve():
                continue
            if temp_nc.resolve() in preexisting_temp_nc:
                continue
            temp_nc.unlink()

    rain_max = float(ds_pr["precip"].max().values)
    log(
        f"   rain written: {tc_rain_out_nc} "
        f"(stored=float32, zlib level={TC_RAIN_COMPRESSION_LEVEL}, "
        f"least_significant_digit={TC_RAIN_LEAST_SIGNIFICANT_DIGIT}, "
        f"zero_below={TC_RAIN_ZERO_BELOW_MMHR} mm/hr, max={rain_max:.3f} mm/hr)"
    )

    sf.setup_config(
        tref=tstart.strftime("%Y%m%d %H%M%S"),
        tstart=tstart.strftime("%Y%m%d %H%M%S"),
        tstop=tstop.strftime("%Y%m%d %H%M%S"),
    )
    sf.set_config("bndfile", "sfincs.bnd")
    sf.set_config("bzsfile", "sfincs.bzs")
    if src_xy.shape[0] > 0 and q_df_final is not None:
        sf.set_config("srcfile", "sfincs.src")
        sf.set_config("disfile", "sfincs.dis")
    sf.set_config("netamprfile", tc_rain_out_nc.name)
    sf.set_config("ampr_block", 1)
    sf.write()

    updates = {
        "tref": tstart.strftime("%Y%m%d %H%M%S"),
        "tstart": tstart.strftime("%Y%m%d %H%M%S"),
        "tstop": tstop.strftime("%Y%m%d %H%M%S"),
        "bndfile": "sfincs.bnd",
        "bzsfile": "sfincs.bzs",
        "netamprfile": tc_rain_out_nc.name,
        "ampr_block": "1",
    }
    if src_xy.shape[0] > 0 and q_df_final is not None:
        updates.update({"srcfile": "sfincs.src", "disfile": "sfincs.dis"})
    base.update_sfincs_inp(str(inp_path), updates)

    base.write_sfincs_bzs_file(str(bzs_path), wl_df)
    if src_xy.shape[0] > 0 and q_df_final is not None:
        dis_path = output_model_root / "sfincs.dis"
        base.write_sfincs_dis_file(str(dis_path), q_df_final)

    log(f"   [OK] completed: {output_model_root}")
    return {
        "case": case_name,
        "adcirc_case": adcirc_case_name,
        "output_model_root": str(output_model_root),
        "bzs": str(bzs_path),
        "rain": str(tc_rain_out_nc),
        "bnd_match": str(out_bnd_match_csv),
        "src_match": str(out_src_match_csv) if out_src_match_csv is not None else "",
    }


# =============================================================================
# 6. CLI
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build SFINCS forcing for future synthetic-TC ADCIRC cases."
    )
    parser.add_argument(
        "--site",
        action="append",
        help="Site name. Can be repeated. Accepts PRD, YRD, BoB, Mississippi, or Misp.",
    )
    parser.add_argument(
        "--all-sites",
        action="store_true",
        help="Run the four target estuaries: PRD, YRD, BoB and Mississippi.",
    )
    parser.add_argument("--scenarios", nargs="+", help="Scenario tags to run.")
    parser.add_argument("--return-periods", nargs="+", help="Return-period tags to run.")
    parser.add_argument("--future-root", help="Override future_adirc_fort root.")
    parser.add_argument("--output-parent", help="Put all output model folders under this directory.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-copy output models from the base model. This is already the default.",
    )
    parser.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Do not re-copy existing output model folders.",
    )
    parser.add_argument("--stop-on-error", action="store_true", help="Stop at the first failed case.")
    parser.add_argument(
        "--flow-mode",
        choices=["stats", "stats_p90", "cama", "none"],
        help="River discharge source. Default is stats.",
    )
    parser.add_argument("--flow-stats-root", help="Directory containing sfincs_boundary_tc_flow_statistics__*.csv.")
    parser.add_argument(
        "--flow-quantile",
        default=None,
        help="Flow quantile suffix used with allocated*, e.g. P50, P90, or P95. Overrides FLOW_STATS_PERCENTILE.",
    )
    parser.add_argument(
        "--flow-percentile",
        default=None,
        help="Flow percentile used with allocatedP*, e.g. 50, 90, or 95. Overrides FLOW_STATS_PERCENTILE.",
    )
    parser.add_argument("--skip-cama", action="store_true", help="Do not write CaMa discharge forcing.")
    parser.add_argument("--skip-flow", action="store_true", help="Do not write river discharge forcing.")
    parser.add_argument("--skip-missing-cama", action="store_true", help="Continue if CaMa inputs are missing.")
    parser.add_argument("--no-vlm", action="store_true", help="Disable VLM subsidence lowering.")
    parser.add_argument("--vlm-historical", action="store_true", help="Also apply VLM to historical_era5.")
    parser.add_argument(
        "--vlm-source",
        choices=["paper_table", "grid_tif"],
        help="VLM source for uniform SFINCS lowering. Default is paper_table.",
    )
    parser.add_argument(
        "--vlm-stat",
        choices=["mean", "median", "area_weighted"],
        help="Paper Table 1 VLM statistic to use when --vlm-source paper_table. Default is area_weighted.",
    )
    parser.add_argument(
        "--rain-compression-level",
        type=int,
        help="NetCDF zlib compression level for TC rain forcing, 0..9. Default is 9.",
    )
    parser.add_argument(
        "--rain-least-significant-digit",
        type=int,
        help="Round TC rain to this many decimal places in mm/hr before compression. Default is 2.",
    )
    parser.add_argument(
        "--rain-zero-below",
        type=float,
        help="Set rain rates below this mm/hr threshold to zero before compression. Default is 0.001.",
    )
    parser.add_argument(
        "--keep-rain-source-nc",
        action="store_true",
        help="Keep the diagnostic tc_rain_rcliper_source.nc file.",
    )
    parser.add_argument(
        "--keep-temp-rain-nc",
        action="store_true",
        help="Keep new temporary HydroMT rain NetCDF files such as sfincs.nc.",
    )
    parser.add_argument(
        "--allow-static-maxele",
        action="store_true",
        help="Allow maxele.63 fallback as a constant bzs boundary.",
    )
    parser.add_argument(
        "--no-static-maxele",
        action="store_true",
        help="Require transient fort.63 / fort.63.nc and reject maxele.63 fallback.",
    )
    return parser.parse_args()


def main():
    global FUTURE_ADCIRC_ROOT, ADCIRC_RESULT_ROOT
    global OUTPUT_PARENT_ROOT, OVERWRITE_EXISTING_OUTPUT_MODEL, STOP_ON_ERROR
    global ENABLE_CAMA_DISCHARGE, SKIP_CAMA_IF_MISSING
    global APPLY_VLM_SUBSIDENCE, APPLY_VLM_TO_HISTORICAL, ALLOW_STATIC_MAXELE63_BZS
    global VLM_SOURCE, PAPER_VLM_RATE_STAT
    global FLOW_FORCING_MODE, FLOW_STATS_ROOT, FLOW_STATS_PERCENTILE
    global TC_RAIN_COMPRESSION_LEVEL, KEEP_TC_RAIN_SOURCE_NC, CLEAN_TEMP_RAIN_NC
    global TC_RAIN_LEAST_SIGNIFICANT_DIGIT, TC_RAIN_ZERO_BELOW_MMHR

    args = parse_args()
    if args.future_root:
        FUTURE_ADCIRC_ROOT = Path(args.future_root)
        ADCIRC_RESULT_ROOT = FUTURE_ADCIRC_ROOT / "maxele63_results"
    if args.output_parent:
        OUTPUT_PARENT_ROOT = args.output_parent
    if args.overwrite:
        OVERWRITE_EXISTING_OUTPUT_MODEL = True
    if args.no_overwrite:
        OVERWRITE_EXISTING_OUTPUT_MODEL = False
    if args.stop_on_error:
        STOP_ON_ERROR = True
    if args.flow_mode:
        if args.flow_mode == "stats_p90":
            FLOW_FORCING_MODE = "stats"
            if not args.flow_quantile and not args.flow_percentile:
                FLOW_STATS_PERCENTILE = 90
        else:
            FLOW_FORCING_MODE = args.flow_mode
    if args.flow_stats_root:
        FLOW_STATS_ROOT = Path(args.flow_stats_root)
    if args.flow_quantile:
        FLOW_STATS_PERCENTILE = normalize_flow_percentile(args.flow_quantile)
    if args.flow_percentile:
        FLOW_STATS_PERCENTILE = normalize_flow_percentile(args.flow_percentile)
    if args.skip_cama:
        ENABLE_CAMA_DISCHARGE = False
        if FLOW_FORCING_MODE == "cama":
            FLOW_FORCING_MODE = "none"
    if args.skip_flow:
        FLOW_FORCING_MODE = "none"
    if args.skip_missing_cama:
        SKIP_CAMA_IF_MISSING = True
    if args.no_vlm:
        APPLY_VLM_SUBSIDENCE = False
    if args.vlm_historical:
        APPLY_VLM_TO_HISTORICAL = True
    if args.vlm_source:
        VLM_SOURCE = args.vlm_source
    if args.vlm_stat:
        PAPER_VLM_RATE_STAT = args.vlm_stat
    if args.rain_compression_level is not None:
        if args.rain_compression_level < 0 or args.rain_compression_level > 9:
            raise ValueError("--rain-compression-level must be between 0 and 9")
        TC_RAIN_COMPRESSION_LEVEL = args.rain_compression_level
    if args.rain_least_significant_digit is not None:
        if args.rain_least_significant_digit < 0:
            raise ValueError("--rain-least-significant-digit must be >= 0")
        TC_RAIN_LEAST_SIGNIFICANT_DIGIT = args.rain_least_significant_digit
    if args.rain_zero_below is not None:
        if args.rain_zero_below < 0.0:
            raise ValueError("--rain-zero-below must be >= 0")
        TC_RAIN_ZERO_BELOW_MMHR = args.rain_zero_below
    if args.keep_rain_source_nc:
        KEEP_TC_RAIN_SOURCE_NC = True
    if args.keep_temp_rain_nc:
        CLEAN_TEMP_RAIN_NC = False
    if args.allow_static_maxele:
        ALLOW_STATIC_MAXELE63_BZS = True
    if args.no_static_maxele:
        ALLOW_STATIC_MAXELE63_BZS = False

    if args.all_sites:
        sites = list(TARGET_SITES)
    elif args.site:
        sites = [normalize_site_name(site) for site in args.site]
    else:
        sites = TARGET_SITES

    scenarios = args.scenarios or SCENARIO_TAGS
    return_periods = args.return_periods or RETURN_PERIOD_TAGS

    log("=======================================================")
    log(" SFINCS future synthetic-TC coupling batch")
    log("=======================================================")
    log(f"future ADCIRC root: {FUTURE_ADCIRC_ROOT}")
    log(f"sites             : {sites}")
    log(f"scenarios         : {scenarios}")
    log(f"return periods    : {return_periods}")
    log(f"overwrite output  : {OVERWRITE_EXISTING_OUTPUT_MODEL}")
    log(f"VLM enabled       : {APPLY_VLM_SUBSIDENCE}")
    log(f"VLM source        : {VLM_SOURCE}")
    if VLM_SOURCE == "paper_table":
        log(f"VLM paper stat    : {PAPER_VLM_RATE_STAT} ({PAPER_VLM_RATE_UNIT})")
    log(f"flow mode         : {FLOW_FORCING_MODE}")
    if FLOW_FORCING_MODE in ("stats", "stats_p90"):
        log(f"flow statistics   : {FLOW_STATS_ROOT} ({flow_stats_quantile_suffix()})")
    log(f"flow case tag     : {flow_case_tag()}")
    log(f"CaMa enabled      : {ENABLE_CAMA_DISCHARGE}")
    log(
        "rain precision    : float32, rounded with least_significant_digit={}, "
        "zero_below={} mm/hr".format(TC_RAIN_LEAST_SIGNIFICANT_DIGIT, TC_RAIN_ZERO_BELOW_MMHR)
    )
    log(f"rain compression  : zlib level {TC_RAIN_COMPRESSION_LEVEL}")
    log(f"keep rain source  : {KEEP_TC_RAIN_SOURCE_NC}")
    log(f"static maxele ok  : {ALLOW_STATIC_MAXELE63_BZS}")

    results = []
    failures = []
    for cfg, scenario, rp, case_name, adcirc_case_name, case_dir in case_iter(sites, scenarios, return_periods):
        try:
            results.append(process_one_case(cfg, scenario, rp, case_name, adcirc_case_name, case_dir))
        except Exception as exc:
            failures.append((case_name, repr(exc)))
            log(f"   [FAILED] {case_name}: {exc!r}")
            if STOP_ON_ERROR:
                raise

    log("\n=======================================================")
    log("Batch summary")
    log("=======================================================")
    log(f"completed: {len(results)}")
    for item in results:
        log(f"   OK {item['case']} -> {item['output_model_root']}")
    log(f"failed   : {len(failures)}")
    for case_name, err in failures:
        log(f"   FAIL {case_name}: {err}")


if __name__ == "__main__":
    main()
