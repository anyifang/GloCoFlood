#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build SFINCS models from global TC coastal partition polygons.

This script contains the full embedded BoB SFINCS modelling workflow, but
replaces the fixed MODEL_BBOX with one vector boundary per partition.

Output per partition:
    global_sfincs_partition_models_cama_<tag>min/No_<P2 label>/
        sfincs.inp, sfincs.dep, sfincs.msk, sfincs.bnd, sfincs.bzs, ...
        region_grid.geojson
        active_domain.geojson
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_HYDROMT_EXAMPLES = Path(
    os.environ.get(
        "HYDROMT_SFINCS_EXAMPLES",
        REPO_ROOT / "external" / "HydroMT-SFINCS" / "examples",
    )
).expanduser()
DEFAULT_CAMA_BASE_DIR = os.environ.get("CAMA_MAP_ROOT")

# Change this tag at the top of the script to switch the CaMa-Flood map
# resolution used for matching. Supported examples: "03", "15".
DEFAULT_CAMA_RES_TAG = "15"


def normalize_cama_res_tag(cama_res_tag: str) -> str:
    tag = str(cama_res_tag).replace("min", "").strip().zfill(2)
    if not tag.isdigit():
        raise ValueError(f"CAMA_RES_TAG must look like '03' or '15', got: {cama_res_tag}")
    return tag


def default_output_dir_for_cama_tag(cama_res_tag: str) -> Path:
    tag = normalize_cama_res_tag(cama_res_tag)
    return SCRIPT_DIR / f"global_sfincs_partition_models_cama_{tag}min"


@dataclass
class Config:
    partition_geojson: Path = (
        SCRIPT_DIR
        / "global_tc_coastal_model_partitions"
        / "global_tc_coastal_model_domain_boundaries.geojson"
    )
    domain_csv: Path = (
        SCRIPT_DIR
        / "global_tc_coastal_model_partitions"
        / "global_tc_coastal_model_domains.csv"
    )
    cama_res_tag: str = DEFAULT_CAMA_RES_TAG
    output_dir: Path | None = None
    status_csv: Path | None = None
    res_m: float = 200.0
    overwrite: bool = True
    stop_on_error: bool = False
    # GTC identifiers are the stable keys consumed by the forcing builders.
    folder_naming: str = "model-id"
    hydromt_examples: Path = DEFAULT_HYDROMT_EXAMPLES
    cama_base_dir: Path | None = (
        Path(DEFAULT_CAMA_BASE_DIR).expanduser() if DEFAULT_CAMA_BASE_DIR else None
    )

    def __post_init__(self):
        self.cama_res_tag = normalize_cama_res_tag(self.cama_res_tag)
        if self.output_dir is None:
            self.output_dir = default_output_dir_for_cama_tag(self.cama_res_tag)
        if self.status_csv is None:
            self.status_csv = self.output_dir / "sfincs_partition_build_status.csv"



EMBEDDED_BOB_TEMPLATE = r'''# -*- coding: utf-8 -*-
import xarray as xr
import rioxarray
import numpy as np
import pandas as pd
import geopandas as gpd
import os
import glob
import rasterio.features
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from shapely.geometry import Point, LineString
from shapely.ops import linemerge, unary_union, nearest_points
from pyproj import Transformer
from scipy.ndimage import distance_transform_edt, binary_dilation, label, binary_fill_holes
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from hydromt_sfincs import SfincsModel
from collections import defaultdict

print("=======================================================")
print("   SFINCS GLOBAL ADAPTIVE MEGA-MODEL")
print("   (Pre-filled 1m River Depth + CFD Dilation Rendering)")
print("=======================================================")

# =========================================================
# 0. 全局配置 (Global Configurations)
# =========================================================
FABDEM_PATH = "data/FABDEM_Coastal_200km/fabdem.vrt"
LAND_POLY_SHP = "data/coastaline/land_polygons.shp"
RIVER_DIR = "data/river/"
LANDCOVER_DIR = "data/Landcover"
MODEL_FOLDER = "BoB_SFINCS"  

MODEL_BBOX = [87.4, 21.5, 91.7, 23.7]
TIME_START = "20220914 000000"
TIME_STOP  = "20220918 000000"

# ★ 自适应 UTM 投影带推算逻辑 ★
center_lon = (MODEL_BBOX[0] + MODEL_BBOX[2]) / 2.0
center_lat = (MODEL_BBOX[1] + MODEL_BBOX[3]) / 2.0
utm_zone = int((center_lon + 180) / 6) + 1
# 纬度 > 0 为北半球 (326xx)，否则为南半球 (327xx)
AUTO_EPSG = 32600 + utm_zone if center_lat > 0 else 32700 + utm_zone

print(f"\n   >>> Auto-calculated EPSG Code: {AUTO_EPSG} (UTM Zone {utm_zone})")

RIVER_SLOPE = 1e-5  
MAX_RIVER_UPSTREAM_M = 500000  
BOX_EDGE_CHECK_CELLS = 15
CAMA_SFINCS_MATCH_MAX_M = 50000
CAMA_BOX_EDGE_PRIORITY_M = 50000
CAMA_WIDE_SECTION_EXTRA_DIST_M = 60000
CAMA_WIDTH_MATCH_MAX_RATIO = 3.0

USE_LOCAL_MSL_CORRECTION = True
LOCAL_MSL_OFFSET_M = 0.0
USE_MDT_GRID_MSL_CORRECTION = True
MDT_PATH = "data/cnes_obs-sl_glo_phy-mdt_my_0.125deg_P20Y_1776683935070.nc"
MDT_VARIABLE = "mdt"

# =========================================================
# CaMa-Flood 分辨率配置：只需要改这里
# 可填 "03" 或 "15"
# =========================================================
CAMA_RES_TAG = "15"
CAMA_BASE_DIR_CANDIDATES = [
    path for path in [
        os.environ.get("CAMA_MAP_ROOT"),
        os.path.join("data", "camaflood", "map_v420"),
    ] if path
]


def resolve_cama_base_dir(cama_res_tag, candidates):
    tag = str(cama_res_tag).replace("min", "").strip().zfill(2)
    target_name = f"glb_{tag}min"
    for path in candidates:
        base_norm = os.path.normpath(path)
        if os.path.basename(base_norm).lower() == target_name and os.path.isdir(base_norm):
            return base_norm
        if os.path.isdir(os.path.join(base_norm, target_name)):
            return base_norm
    return candidates[0]


CAMA_BASE_DIR = resolve_cama_base_dir(CAMA_RES_TAG, CAMA_BASE_DIR_CANDIDATES)

UPAREA_THRESHOLD_KM2 = 10000
SPECIFIC_YIELD = 0.03

# 核心修正：分离淘汰阈值与渲染阈值
MIN_PHYSICAL_WIDTH_M = 50.0
CFD_RESOLVE_CELLS = 1.0
# 核心修正：分离淘汰阈值与渲染阈值
SEA_DIKE_HEIGHT_M = 0.0
RIVER_DIKE_COAST_HEIGHT_M = SEA_DIKE_HEIGHT_M
RIVER_DIKE_INLAND_HEIGHT_M = 0.0
RIVER_DIKE_GRADIENT_DISTANCE_M = 150000.0
BELOW_SEA_LEVEL_MANNING = 0.01
DEFAULT_LAND_MANNING = 0.04
LANDCOVER_NODATA_VALUE = 255
MANNING_MAP = {
    111: 0.10, 112: 0.10, 113: 0.10, 114: 0.10, 115: 0.10, 116: 0.10,  # closed forest
    121: 0.08, 122: 0.08, 123: 0.08, 124: 0.08, 125: 0.08, 126: 0.08,  # open forest
    20: 0.06,   # shrubland
    30: 0.04,   # herbaceous vegetation
    40: 0.04,   # cropland
    50: 0.10,   # built-up
    60: 0.03,   # bare / sparse vegetation
    70: 0.03,   # snow & ice
    80: 0.025,  # permanent water bodies
    90: 0.08,   # herbaceous wetland
    100: 0.04,  # moss & lichen
    200: 0.01  # open sea
}

# 预先定义，避免后面未赋值时报错
q_vals_dynamic = []
valid_src_points = []
CAMA_INLET_TARGETS = {}


def fill_internal_holes_with_elevation_interp(
    sf,
    original_land_mask,
    valid_dem_mask_2d,
    final_river_mask_2d,
    original_dep=None,
    ring_iterations=3,
):
    """
    回填由于河道宽度筛选导致的内部空洞：
    1) 识别被 active mask 包围的内部 hole
    2) 将这些 hole 的 msk 恢复为 1（普通陆地区域）
    3) 用周边陆地高程点插值回填
       - 优先 linear
       - 再用 nearest 补齐
    """
    current_msk = sf.grid["msk"].values.copy()
    current_dep = sf.grid["dep"].values.copy()

    active_mask = (current_msk > 0) & valid_dem_mask_2d
    land_footprint = (np.asarray(original_land_mask.values) > 0) & valid_dem_mask_2d

    filled_mask = binary_fill_holes(active_mask)
    internal_hole_mask = filled_mask & (~active_mask) & land_footprint

    if not np.any(internal_hole_mask):
        print("   >>> No internal holes detected after river screening.")
        return internal_hole_mask

    print(f"   >>> Detected {int(internal_hole_mask.sum())} internal hole cells. Filling them...")

    current_msk[internal_hole_mask] = 1

    if original_dep is None:
        original_dep_values = current_dep
    else:
        original_dep_values = np.asarray(
            original_dep.values if hasattr(original_dep, "values") else original_dep
        )

    z_fill = original_dep_values[internal_hole_mask]

    if np.any(~np.isfinite(z_fill)):
        known_original_mask = np.isfinite(original_dep_values)
        if np.any(known_original_mask):
            _, indices = distance_transform_edt(~known_original_mask, return_indices=True)
            fallback_dep = original_dep_values[indices[0], indices[1]]
            z_fill = np.where(
                np.isfinite(z_fill),
                z_fill,
                fallback_dep[internal_hole_mask]
            )

    if np.any(~np.isfinite(z_fill)):
        known_current_mask = (current_msk > 0) & np.isfinite(current_dep)
        _, indices = distance_transform_edt(~known_current_mask, return_indices=True)
        fallback_dep = current_dep[indices[0], indices[1]]
        z_fill = np.where(
            np.isfinite(z_fill),
            z_fill,
            fallback_dep[internal_hole_mask]
        )

    current_dep[internal_hole_mask] = z_fill

    sf.grid["msk"] = xr.DataArray(
        current_msk,
        coords=sf.grid["msk"].coords,
        dims=sf.grid["msk"].dims,
    )
    sf.grid["dep"] = xr.DataArray(
        current_dep,
        coords=sf.grid["dep"].coords,
        dims=sf.grid["dep"].dims,
    )

    print("   >>> Internal holes filled using original FABDEM elevations.")
    return internal_hole_mask


def get_cama_config(cama_res_tag, cama_base_dir):
    """
    根据 CaMa-Flood 分辨率标签自动生成：
    - 地图目录
    - uparea / nextxy 文件路径
    - 网格分辨率（度）
    - 全局网格大小 nx, ny
    - 匹配搜索半径 search_radius
    """
    tag = str(cama_res_tag).replace("min", "").strip().zfill(2)

    if not tag.isdigit():
        raise ValueError(f"CAMA_RES_TAG 必须类似 '03' 或 '15'，当前为: {cama_res_tag}")

    minutes = int(tag)
    res_deg = minutes / 60.0

    nx = int(round(360.0 / res_deg))
    ny = int(round(180.0 / res_deg))

    base_norm = os.path.normpath(cama_base_dir)
    if os.path.basename(base_norm).lower() == f"glb_{tag}min":
        cama_dir = base_norm
    else:
        cama_dir = os.path.join(base_norm, f"glb_{tag}min")
    uparea_bin = os.path.join(cama_dir, "uparea.bin")
    nextxy_bin = os.path.join(cama_dir, "nextxy.bin")
    width_bin = os.path.join(cama_dir, "width.bin")

    if minutes <= 3:
        search_radius = 0.05
    elif minutes <= 6:
        search_radius = 0.1
    else:
        search_radius = 0.3

    if not os.path.isdir(cama_dir):
        raise FileNotFoundError(f"CaMa 地图目录不存在: {cama_dir}")
    if not os.path.exists(uparea_bin):
        raise FileNotFoundError(f"找不到 uparea.bin: {uparea_bin}")
    if not os.path.exists(nextxy_bin):
        raise FileNotFoundError(f"找不到 nextxy.bin: {nextxy_bin}")

    if not os.path.exists(width_bin):
        width_bin = None

    print("\n=======================================================")
    print(f"   >>> CaMa-Flood resolution tag : {tag} min")
    print(f"   >>> CaMa-Flood grid size      : nx={nx}, ny={ny}")
    print(f"   >>> CaMa-Flood resolution     : {res_deg:.6f} degree")
    print(f"   >>> CaMa-Flood search radius  : {search_radius} cell(s)")
    print(f"   >>> CaMa-Flood map dir        : {cama_dir}")
    print("=======================================================")

    return {
        "tag": tag,
        "minutes": minutes,
        "res_deg": res_deg,
        "nx": nx,
        "ny": ny,
        "search_radius": search_radius,
        "dir": cama_dir,
        "uparea_bin": uparea_bin,
        "nextxy_bin": nextxy_bin,
        "width_bin": width_bin,
    }


def find_best_cama_pixel_for_point(
    lon, lat,
    uparea_global,
    cama_lons, cama_lats,
    cama_res,
    cama_nx, cama_ny,
    search_radius_cells,
    uparea_threshold_km2
):
    """
    给一个 SFINCS 边界点（lon, lat），在 CaMa 网格中寻找最合适的上游面积像元。
    返回:
        best_idx   : (iy, ix) or None
        best_up    : m2
        best_dist  : degree distance
    """
    center_ix = int(np.floor((lon - (-180)) / cama_res))
    center_iy = int(np.floor((90 - lat) / cama_res))

    search_radius = max(1, int(np.ceil(search_radius_cells)))

    best_dist = np.inf
    best_idx = None
    best_up = -1.0

    for dy in range(-search_radius, search_radius + 1):
        for dx in range(-search_radius, search_radius + 1):
            iy = center_iy + dy
            ix = center_ix + dx

            if 0 <= iy < cama_ny and 0 <= ix < cama_nx:
                current_up = uparea_global[iy, ix]
                if current_up >= uparea_threshold_km2 * 1e6:
                    c_lon = cama_lons[ix]
                    c_lat = cama_lats[iy]
                    dist = np.sqrt((lon - c_lon) ** 2 + (lat - c_lat) ** 2)

                    if dist < best_dist:
                        best_dist = dist
                        best_idx = (iy, ix)
                        best_up = current_up

    return best_idx, best_up, best_dist


def choose_best_cama_for_component(component_records):
    """
    一个连续边界连通域内部，可能有多个点分别匹配到不同的 CaMa 像元。
    这里强制整个连通域只选一个 CaMa 像元。

    选择规则：
    1) 命中点数最多
    2) 若并列，则平均距离最小
    3) 若再并列，则 uparea 最大
    """
    stats = defaultdict(list)
    up_map = {}

    for rec in component_records:
        cama_idx = rec["cama_idx"]
        stats[cama_idx].append(rec["dist"])
        up_map[cama_idx] = rec["uparea"]

    best_key = None
    best_tuple = None

    for cama_idx, dists in stats.items():
        score = (
            len(dists),
            -np.mean(dists),
            up_map[cama_idx]
        )
        if (best_tuple is None) or (score > best_tuple):
            best_tuple = score
            best_key = cama_idx

    return best_key


def infer_bbox_side_for_lonlat(lon, lat, model_bbox):
    west, south, east, north = model_bbox
    candidates = [
        ("west", abs(float(lon) - west)),
        ("east", abs(float(lon) - east)),
        ("south", abs(float(lat) - south)),
        ("north", abs(float(lat) - north)),
    ]
    return min(candidates, key=lambda rec: rec[1])[0]


def geometry_points_for_target(geom):
    if geom is None or geom.is_empty:
        return []
    geom_type = geom.geom_type
    if geom_type == "Point":
        return [geom]
    if geom_type == "MultiPoint":
        return list(geom.geoms)
    if geom_type in ("LineString", "LinearRing"):
        coords = list(geom.coords)
        if not coords:
            return []
        return [Point(coords[0]), Point(coords[-1])]
    if hasattr(geom, "geoms"):
        points = []
        for part in geom.geoms:
            points.extend(geometry_points_for_target(part))
        return points
    return []


def cama_active_domain_crossing_target(
    upstream_lon,
    upstream_lat,
    inside_lon,
    inside_lat,
    model_bbox,
):
    if "ACTIVE_DOMAIN_GEOM_WGS" not in globals():
        return project_lonlat_to_bbox_edge(inside_lon, inside_lat, model_bbox)

    line = LineString([
        (float(upstream_lon), float(upstream_lat)),
        (float(inside_lon), float(inside_lat)),
    ])
    inside_pt = Point(float(inside_lon), float(inside_lat))
    try:
        intersection = line.intersection(ACTIVE_DOMAIN_GEOM_WGS.boundary)
        target_points = geometry_points_for_target(intersection)
        if target_points:
            target_pt = min(target_points, key=lambda pt: pt.distance(inside_pt))
        else:
            target_pt = nearest_points(inside_pt, ACTIVE_DOMAIN_GEOM_WGS.boundary)[1]
        match_lon = float(target_pt.x)
        match_lat = float(target_pt.y)
        side = infer_bbox_side_for_lonlat(match_lon, match_lat, model_bbox)
        return side, match_lon, match_lat
    except Exception:
        return project_lonlat_to_bbox_edge(inside_lon, inside_lat, model_bbox)


def find_cama_domain_inlets_for_debug(
    uparea_global,
    nextx,
    nexty,
    cama_lons,
    cama_lats,
    model_bbox,
    uparea_threshold_km2,
):
    global CAMA_INLET_TARGETS
    CAMA_INLET_TARGETS = {}

    west, south, east, north = model_bbox
    up_threshold = uparea_threshold_km2 * 1e6

    valid = uparea_global >= up_threshold

    dlon = abs(cama_lons[1] - cama_lons[0])
    dlat = abs(cama_lats[1] - cama_lats[0])
    pad = 2
    lon_search = (cama_lons >= west - pad * dlon) & (cama_lons <= east + pad * dlon)
    lat_search = (cama_lats >= south - pad * dlat) & (cama_lats <= north + pad * dlat)
    search = lat_search[:, None] & lon_search[None, :]

    inside = np.zeros_like(search, dtype=bool)
    active_geom = globals().get("ACTIVE_DOMAIN_GEOM_WGS", None)
    if active_geom is not None and not active_geom.is_empty:
        y_search_all, x_search_all = np.where(search)
        for yy, xx in zip(y_search_all, x_search_all):
            inside[yy, xx] = bool(
                active_geom.covers(Point(float(cama_lons[xx]), float(cama_lats[yy])))
            )
    else:
        lon_inside = (cama_lons >= west) & (cama_lons <= east)
        lat_inside = (cama_lats >= south) & (cama_lats <= north)
        inside = lat_inside[:, None] & lon_inside[None, :]

    inlet_keys = set()
    has_inside_parent = set()
    outside_parents_for_inside = defaultdict(list)

    def register_target(inlet_key, upstream_key=None, target_method="active_domain_boundary_crossing"):
        in_y, in_x = int(inlet_key[0]), int(inlet_key[1])
        in_lon = float(cama_lons[in_x])
        in_lat = float(cama_lats[in_y])
        if upstream_key is not None:
            up_y, up_x = int(upstream_key[0]), int(upstream_key[1])
            up_lon = float(cama_lons[up_x])
            up_lat = float(cama_lats[up_y])
            target_side, match_lon, match_lat = cama_active_domain_crossing_target(
                up_lon,
                up_lat,
                in_lon,
                in_lat,
                model_bbox,
            )
            CAMA_INLET_TARGETS[(in_y, in_x)] = {
                "target_bbox_side": target_side,
                "match_lon": float(match_lon),
                "match_lat": float(match_lat),
                "target_method": target_method,
                "entry_from_row": int(up_y),
                "entry_from_col": int(up_x),
                "entry_from_lon": float(up_lon),
                "entry_from_lat": float(up_lat),
                "entry_to_row": int(in_y),
                "entry_to_col": int(in_x),
                "entry_to_lon": float(in_lon),
                "entry_to_lat": float(in_lat),
            }
        else:
            target_side, match_lon, match_lat = project_lonlat_to_bbox_edge(
                in_lon,
                in_lat,
                model_bbox,
            )
            CAMA_INLET_TARGETS[(in_y, in_x)] = {
                "target_bbox_side": target_side,
                "match_lon": float(match_lon),
                "match_lat": float(match_lat),
                "target_method": "nearest_active_domain_boundary",
                "entry_to_row": int(in_y),
                "entry_to_col": int(in_x),
                "entry_to_lon": float(in_lon),
                "entry_to_lat": float(in_lat),
            }

    def downstream_inside_valid(y0, x0):
        y1 = int(nexty[y0, x0])
        x1 = int(nextx[y0, x0])
        return (
            0 <= y1 < uparea_global.shape[0]
            and 0 <= x1 < uparea_global.shape[1]
            and inside[y1, x1]
            and valid[y1, x1]
        )

    y_search, x_search = np.where(search & valid)
    for y0, x0 in zip(y_search, x_search):
        y1 = int(nexty[y0, x0])
        x1 = int(nextx[y0, x0])
        if not (0 <= y1 < uparea_global.shape[0] and 0 <= x1 < uparea_global.shape[1]):
            continue

        if inside[y0, x0] and inside[y1, x1]:
            has_inside_parent.add((y1, x1))
        elif (not inside[y0, x0]) and inside[y1, x1]:
            if downstream_inside_valid(y1, x1):
                inlet_keys.add((y1, x1))
                outside_parents_for_inside[(y1, x1)].append((y0, x0))

    for inlet_key, upstream_keys in outside_parents_for_inside.items():
        upstream_key = max(
            upstream_keys,
            key=lambda key: float(uparea_global[int(key[0]), int(key[1])]),
        )
        register_target(inlet_key, upstream_key=upstream_key)

    y_inside, x_inside = np.where(inside & valid)
    for y0, x0 in zip(y_inside, x_inside):
        y1 = int(nexty[y0, x0])
        x1 = int(nextx[y0, x0])
        downstream_inside = (
            0 <= y1 < uparea_global.shape[0]
            and 0 <= x1 < uparea_global.shape[1]
            and inside[y1, x1]
            and valid[y1, x1]
        )
        if active_geom is not None and not active_geom.is_empty:
            near_edge = Point(float(cama_lons[x0]), float(cama_lats[y0])).distance(
                active_geom.boundary
            ) <= 0.75 * max(dlon, dlat)
        else:
            near_edge = (
                abs(cama_lons[x0] - west) <= 0.5 * dlon
                or abs(cama_lons[x0] - east) <= 0.5 * dlon
                or abs(cama_lats[y0] - south) <= 0.5 * dlat
                or abs(cama_lats[y0] - north) <= 0.5 * dlat
            )
        if near_edge and downstream_inside and (y0, x0) not in has_inside_parent:
            inlet_keys.add((y0, x0))
            if (y0, x0) not in CAMA_INLET_TARGETS:
                register_target((y0, x0), upstream_key=None)

    return sorted(inlet_keys, key=lambda k: uparea_global[k[0], k[1]], reverse=True)


def _ordered_lonlat_from_xy(points_xy, transformer_to_latlon):
    coords = np.array(points_xy, dtype=float)
    if len(coords) == 0:
        return []
    if len(coords) == 1:
        ordered = coords
    else:
        centered = coords - coords.mean(axis=0)
        try:
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
            order = np.argsort(centered @ vh[0])
        except Exception:
            order = np.lexsort((coords[:, 1], coords[:, 0]))
        ordered = coords[order]

    return [
        transformer_to_latlon.transform(float(x_val), float(y_val))
        for x_val, y_val in ordered
    ]


def _edge_band_mask(shape, iterations):
    edge_mask = np.zeros(shape, dtype=bool)
    edge_mask[0, :] = True
    edge_mask[-1, :] = True
    edge_mask[:, 0] = True
    edge_mask[:, -1] = True
    if iterations <= 0:
        return edge_mask
    return binary_dilation(edge_mask, iterations=int(iterations))


def write_pre_match_boundary_geojsons(
    model_folder,
    upstream_boundary_mask,
    boundary_labels,
    boundary_source_map,
    valid_river_mask,
    dist_2d,
    sf,
    transformer_to_latlon,
    cama_inlets,
    uparea_global,
    nextx,
    nexty,
    cama_lons,
    cama_lats,
):
    os.makedirs(model_folder, exist_ok=True)
    for old_name in [
        "pre_match_sfincs_boundaries.geojson",
        "pre_match_sfincs_flow_boundaries.geojson",
        "pre_match_camaflood_boundaries.geojson",
        "pre_match_camaflood_flow_boundaries.geojson",
    ]:
        old_path = os.path.join(model_folder, old_name)
        if os.path.exists(old_path):
            os.remove(old_path)

    sfincs_features = []
    source_labels = {
        1: "max_distance",
        2: "box_edge",
        3: "max_distance_and_box_edge",
        4: "low_elevation_edge",
        5: "partition_edge",
    }
    y_src, x_src = np.where(upstream_boundary_mask)
    for iy_src, ix_src in zip(y_src, x_src):
        x_coord = float(sf.grid["x"].values[ix_src])
        y_coord = float(sf.grid["y"].values[iy_src])
        lon_pt, lat_pt = transformer_to_latlon.transform(x_coord, y_coord)
        comp_id = int(boundary_labels[iy_src, ix_src])
        source_code = int(boundary_source_map[iy_src, ix_src])
        sfincs_features.append({
            "geometry": Point(lon_pt, lat_pt),
            "feature_type": "sfincs_boundary_point",
            "boundary_component_id": comp_id,
            "boundary_source": source_labels.get(source_code, "unknown"),
            "source_code": source_code,
            "grid_row": int(iy_src),
            "grid_col": int(ix_src),
            "sfincs_x": x_coord,
            "sfincs_y": y_coord,
            "dist_km": float(dist_2d[iy_src, ix_src] / 1000.0),
            "is_valid_river": bool(valid_river_mask[iy_src, ix_src]),
        })

    component_ids = np.unique(boundary_labels[upstream_boundary_mask])
    component_ids = component_ids[component_ids > 0]
    for comp_id in component_ids:
        comp_mask = boundary_labels == comp_id
        yy_comp, xx_comp = np.where(comp_mask)
        if len(yy_comp) == 0:
            continue

        points_xy = [
            (
                float(sf.grid["x"].values[ix_src]),
                float(sf.grid["y"].values[iy_src]),
            )
            for iy_src, ix_src in zip(yy_comp, xx_comp)
        ]
        lonlat = _ordered_lonlat_from_xy(points_xy, transformer_to_latlon)
        geom = LineString(lonlat) if len(lonlat) > 1 else Point(lonlat[0])

        sfincs_features.append({
            "geometry": geom,
            "feature_type": "sfincs_boundary_component",
            "boundary_component_id": int(comp_id),
            "n_points": int(len(yy_comp)),
            "boundary_source": source_labels.get(
                int(np.max(boundary_source_map[yy_comp, xx_comp])),
                "unknown",
            ),
            "source_code": int(np.max(boundary_source_map[yy_comp, xx_comp])),
            "mean_dist_km": float(np.nanmean(dist_2d[yy_comp, xx_comp]) / 1000.0),
            "max_dist_km": float(np.nanmax(dist_2d[yy_comp, xx_comp]) / 1000.0),
        })

    if sfincs_features:
        out_sfincs = os.path.join(model_folder, "pre_match_sfincs_boundaries.geojson")
        sfincs_pre_match_gdf = gpd.GeoDataFrame(sfincs_features, crs="EPSG:4326")
        sfincs_pre_match_gdf.to_file(out_sfincs, driver="GeoJSON")
        print(f"   >>> Pre-match SFINCS boundaries written: {out_sfincs}")
        out_sfincs_flow = os.path.join(model_folder, "pre_match_sfincs_flow_boundaries.geojson")
        sfincs_pre_match_gdf.to_file(out_sfincs_flow, driver="GeoJSON")
        print(f"   >>> Pre-match SFINCS flow boundaries written: {out_sfincs_flow}")

    cama_features = []
    cama_flow_features = []
    for rank, (cama_y, cama_x) in enumerate(cama_inlets, start=1):
        c_lon = float(cama_lons[cama_x])
        c_lat = float(cama_lats[cama_y])
        target_info = CAMA_INLET_TARGETS.get((cama_y, cama_x), {})
        if target_info:
            target_side = target_info["target_bbox_side"]
            match_lon = float(target_info["match_lon"])
            match_lat = float(target_info["match_lat"])
        else:
            target_side, match_lon, match_lat = project_lonlat_to_bbox_edge(
                c_lon,
                c_lat,
                MODEL_BBOX,
            )
        nxt_x = int(nextx[cama_y, cama_x])
        nxt_y = int(nexty[cama_y, cama_x])

        props = {
            "feature_type": "camaflood_boundary_point",
            "rank": int(rank),
            "cama_row": int(cama_y),
            "cama_col": int(cama_x),
            "cama_lon": c_lon,
            "cama_lat": c_lat,
            "target_bbox_side": target_side,
            "match_lon": float(match_lon),
            "match_lat": float(match_lat),
            "uparea_km2": float(uparea_global[cama_y, cama_x] / 1e6),
            "threshold_km2": float(UPAREA_THRESHOLD_KM2),
        }
        for target_key in [
            "target_method",
            "entry_from_row",
            "entry_from_col",
            "entry_from_lon",
            "entry_from_lat",
            "entry_to_row",
            "entry_to_col",
            "entry_to_lon",
            "entry_to_lat",
        ]:
            if target_key in target_info:
                props[target_key] = target_info[target_key]

        if 0 <= nxt_y < uparea_global.shape[0] and 0 <= nxt_x < uparea_global.shape[1]:
            props.update({
                "next_row": int(nxt_y),
                "next_col": int(nxt_x),
                "next_lon": float(cama_lons[nxt_x]),
                "next_lat": float(cama_lats[nxt_y]),
                "next_uparea_km2": float(uparea_global[nxt_y, nxt_x] / 1e6),
            })

        cama_features.append({
            "geometry": Point(c_lon, c_lat),
            **props,
        })
        cama_flow_features.append({
            "geometry": Point(c_lon, c_lat),
            **{
                **props,
                "feature_type": "camaflood_flow_boundary_point",
            },
        })
        cama_flow_features.append({
            "geometry": LineString([(c_lon, c_lat), (match_lon, match_lat)]),
            **{
                **props,
                "feature_type": "camaflood_bbox_match_link",
            },
        })

        if "next_lon" in props:
            cama_features.append({
                "geometry": LineString([(c_lon, c_lat), (props["next_lon"], props["next_lat"])]),
                **{
                    **props,
                    "feature_type": "camaflood_nextxy_link",
                },
            })
            cama_flow_features.append({
                "geometry": LineString([(c_lon, c_lat), (props["next_lon"], props["next_lat"])]),
                **{
                    **props,
                    "feature_type": "camaflood_flow_nextxy_link",
                },
            })

    if cama_features:
        out_cama = os.path.join(model_folder, "pre_match_camaflood_boundaries.geojson")
        cama_pre_match_gdf = gpd.GeoDataFrame(cama_features, crs="EPSG:4326")
        cama_pre_match_gdf.to_file(out_cama, driver="GeoJSON")
        print(f"   >>> Pre-match CaMa-Flood boundaries written: {out_cama}")
        out_cama_flow = os.path.join(model_folder, "pre_match_camaflood_flow_boundaries.geojson")
        gpd.GeoDataFrame(cama_flow_features, crs="EPSG:4326").to_file(out_cama_flow, driver="GeoJSON")
        print(f"   >>> Pre-match CaMa-Flood flow boundaries written: {out_cama_flow}")


def boundary_source_label_from_codes(source_codes):
    codes = set(int(code) for code in np.asarray(source_codes).ravel() if int(code) > 0)
    if 5 in codes:
        return "partition_edge"
    if 3 in codes:
        return "max_distance_and_box_edge"
    if 2 in codes:
        return "box_edge"
    if 1 in codes:
        return "max_distance"
    if 4 in codes:
        return "low_elevation_edge"
    return "unknown"


def cama_flows_to_any(start_key, target_keys, nextx, nexty, max_steps=2000):
    target_keys = set(target_keys)
    trace_y, trace_x = int(start_key[0]), int(start_key[1])
    visited = set()

    for _ in range(max_steps):
        nxt_x = int(nextx[trace_y, trace_x])
        nxt_y = int(nexty[trace_y, trace_x])

        if nxt_x < 0 or nxt_y < 0:
            return None
        if (nxt_y, nxt_x) == (trace_y, trace_x):
            return None
        if (nxt_y, nxt_x) in visited:
            return None
        if (nxt_y, nxt_x) in target_keys:
            return (nxt_y, nxt_x)

        visited.add((nxt_y, nxt_x))
        trace_y, trace_x = nxt_y, nxt_x

    return None


def project_lonlat_to_bbox_edge(lon, lat, model_bbox):
    west, south, east, north = model_bbox
    if "ACTIVE_DOMAIN_GEOM_WGS" in globals():
        try:
            point = Point(float(lon), float(lat))
            nearest_on_boundary = nearest_points(point, ACTIVE_DOMAIN_GEOM_WGS.boundary)[1]
            edge_lon = float(nearest_on_boundary.x)
            edge_lat = float(nearest_on_boundary.y)
            edge_candidates = [
                ("west", abs(edge_lon - west)),
                ("east", abs(edge_lon - east)),
                ("south", abs(edge_lat - south)),
                ("north", abs(edge_lat - north)),
            ]
            side, _ = min(edge_candidates, key=lambda rec: rec[1])
            return side, edge_lon, edge_lat
        except Exception:
            pass

    lon_clamped = min(max(float(lon), west), east)
    lat_clamped = min(max(float(lat), south), north)

    edge_candidates = [
        ("west", west, lat_clamped, abs(float(lon) - west)),
        ("east", east, lat_clamped, abs(float(lon) - east)),
        ("south", lon_clamped, south, abs(float(lat) - south)),
        ("north", lon_clamped, north, abs(float(lat) - north)),
    ]
    side, edge_lon, edge_lat, _ = min(edge_candidates, key=lambda rec: rec[3])
    return side, edge_lon, edge_lat


def component_box_edge_sides(yy_comp, xx_comp, shape, edge_cells):
    nrows, ncols = shape
    sides = set()
    edge_cells = max(1, int(edge_cells))

    if np.any(xx_comp <= edge_cells):
        sides.add("west")
    if np.any(xx_comp >= ncols - 1 - edge_cells):
        sides.add("east")
    if np.any(yy_comp <= edge_cells):
        sides.add("north")
    if np.any(yy_comp >= nrows - 1 - edge_cells):
        sides.add("south")

    return sides


def select_single_row_boundary_cells(grid_cells, target_side, shape):
    """
    Reduce a box-edge boundary component from a band/patch to one cell-thick row.
    The row is selected on the side where the CaMa inlet enters the SFINCS box.
    """
    if len(grid_cells) == 0:
        return []

    cells = np.array(grid_cells, dtype=int)
    yy_cells = cells[:, 0]
    xx_cells = cells[:, 1]
    nrows, ncols = shape

    if target_side == "north":
        selected = cells[yy_cells == np.min(yy_cells)]
        order = np.argsort(selected[:, 1])
    elif target_side == "south":
        selected = cells[yy_cells == np.max(yy_cells)]
        order = np.argsort(selected[:, 1])
    elif target_side == "west":
        selected = cells[xx_cells == np.min(xx_cells)]
        order = np.argsort(selected[:, 0])
    elif target_side == "east":
        selected = cells[xx_cells == np.max(xx_cells)]
        order = np.argsort(selected[:, 0])
    else:
        edge_dists = {
            "north": int(np.min(yy_cells)),
            "south": int((nrows - 1) - np.max(yy_cells)),
            "west": int(np.min(xx_cells)),
            "east": int((ncols - 1) - np.max(xx_cells)),
        }
        inferred_side = min(edge_dists, key=edge_dists.get)
        return select_single_row_boundary_cells(grid_cells, inferred_side, shape)

    if len(selected) < 2 and len(cells) <= 3:
        if target_side in ("north", "south"):
            order = np.argsort(cells[:, 1])
        elif target_side in ("west", "east"):
            order = np.argsort(cells[:, 0])
        else:
            order = np.lexsort((cells[:, 1], cells[:, 0]))
        selected = cells[order]

    selected = selected[order]
    return [(int(iy), int(ix)) for iy, ix in selected]


def select_nearest_contiguous_run(line_cells, target_x, target_y, sf, target_side):
    if len(line_cells) <= 1:
        return line_cells

    sorted_cells = list(line_cells)
    runs = []
    current = [sorted_cells[0]]
    gap_cells = 2

    for cell in sorted_cells[1:]:
        prev = current[-1]
        if target_side in ("north", "south"):
            is_next = abs(cell[1] - prev[1]) <= gap_cells
        else:
            is_next = abs(cell[0] - prev[0]) <= gap_cells

        if is_next:
            current.append(cell)
        else:
            runs.append(current)
            current = [cell]
    runs.append(current)

    if len(runs) == 1:
        return runs[0]

    def run_score(run):
        coords = np.array([
            (
                float(sf.grid["x"].values[ix_src]),
                float(sf.grid["y"].values[iy_src]),
            )
            for iy_src, ix_src in run
        ])
        dists = np.sqrt((coords[:, 0] - target_x) ** 2 + (coords[:, 1] - target_y) ** 2)
        return (float(np.min(dists)), -len(run))

    return min(runs, key=run_score)


def select_complete_boundary_front_cells(grid_cells, target_side, shape):
    """
    Build a one-cell-thick front for the already matched SFINCS boundary component.
    This only uses cells inside the matched component, so supplementation cannot
    leak into another river.
    """
    if len(grid_cells) == 0:
        return []

    cells = np.array(grid_cells, dtype=int)
    yy_cells = cells[:, 0]
    xx_cells = cells[:, 1]
    nrows, ncols = shape

    if target_side == "north":
        selected = []
        for ix_val in np.unique(xx_cells):
            col_cells = cells[xx_cells == ix_val]
            selected.append(col_cells[np.argmin(col_cells[:, 0])])
        selected = np.array(selected, dtype=int)
        order = np.argsort(selected[:, 1])
    elif target_side == "south":
        selected = []
        for ix_val in np.unique(xx_cells):
            col_cells = cells[xx_cells == ix_val]
            selected.append(col_cells[np.argmax(col_cells[:, 0])])
        selected = np.array(selected, dtype=int)
        order = np.argsort(selected[:, 1])
    elif target_side == "west":
        selected = []
        for iy_val in np.unique(yy_cells):
            row_cells = cells[yy_cells == iy_val]
            selected.append(row_cells[np.argmin(row_cells[:, 1])])
        selected = np.array(selected, dtype=int)
        order = np.argsort(selected[:, 0])
    elif target_side == "east":
        selected = []
        for iy_val in np.unique(yy_cells):
            row_cells = cells[yy_cells == iy_val]
            selected.append(row_cells[np.argmax(row_cells[:, 1])])
        selected = np.array(selected, dtype=int)
        order = np.argsort(selected[:, 0])
    else:
        edge_dists = {
            "north": int(np.min(yy_cells)),
            "south": int((nrows - 1) - np.max(yy_cells)),
            "west": int(np.min(xx_cells)),
            "east": int((ncols - 1) - np.max(xx_cells)),
        }
        inferred_side = min(edge_dists, key=edge_dists.get)
        return select_complete_boundary_front_cells(grid_cells, inferred_side, shape)

    selected = selected[order]
    return [(int(iy), int(ix)) for iy, ix in selected]


def complete_inflow_cells_with_matched_boundary_front(
    component_grid_cells,
    initial_line_cells,
    target_x,
    target_y,
    sf,
    target_side,
    shape,
):
    """
    Check whether the selected SFINCS inflow cells cover the matched boundary
    component's full front. If not, supplement only within this matched component.
    """
    complete_front = select_complete_boundary_front_cells(
        component_grid_cells,
        target_side,
        shape,
    )
    if not complete_front:
        return initial_line_cells, {
            "front_npts": 0,
            "initial_npts": int(len(initial_line_cells)),
            "final_npts": int(len(initial_line_cells)),
            "completed": False,
        }

    front_run = select_nearest_contiguous_run(
        complete_front,
        target_x,
        target_y,
        sf,
        target_side,
    )
    if not front_run:
        front_run = complete_front

    initial_set = set((int(iy), int(ix)) for iy, ix in initial_line_cells)
    front_set = set((int(iy), int(ix)) for iy, ix in front_run)
    missing_cells = front_set - initial_set
    completed = len(missing_cells) > 0

    final_cells = front_run if completed else initial_line_cells
    return final_cells, {
        "front_npts": int(len(front_run)),
        "initial_npts": int(len(initial_line_cells)),
        "final_npts": int(len(final_cells)),
        "completed": bool(completed),
    }


def build_matched_downstream_river_mask(
    valid_river_mask,
    matched_boundary_grid_cells,
    dist_2d,
    pixel_size_m,
):
    """
    Preserve every matched SFINCS river downstream of its inflow section.

    Important deltas can have multiple valid upstream inlets in different
    distributaries. Keeping only the longest/upstream section lets the normal
    width screen remove other already-matched inlet sections.
    """
    protected_mask = np.zeros(valid_river_mask.shape, dtype=bool)
    if not matched_boundary_grid_cells or not np.any(valid_river_mask):
        return protected_mask

    conn_structure = np.ones((3, 3), dtype=np.uint8)
    river_component_labels, _ = label(
        valid_river_mask.astype(np.uint8),
        structure=conn_structure,
    )

    finite_dist = dist_2d[np.isfinite(dist_2d)]
    if finite_dist.size == 0:
        return protected_mask

    min_drop = 0.05 * pixel_size_m
    max_steps = max(100, int(np.nanmax(finite_dist) / max(pixel_size_m, 1.0)) + 20)

    def trace_downstream_from_section(grid_cells):
        if len(grid_cells) == 0:
            return

        cells = np.array(grid_cells, dtype=int)
        yy_cells = cells[:, 0]
        xx_cells = cells[:, 1]
        inside = (
            (yy_cells >= 0)
            & (yy_cells < valid_river_mask.shape[0])
            & (xx_cells >= 0)
            & (xx_cells < valid_river_mask.shape[1])
        )
        yy_cells = yy_cells[inside]
        xx_cells = xx_cells[inside]
        if len(yy_cells) == 0:
            return

        river_ids = np.unique(river_component_labels[yy_cells, xx_cells])
        river_ids = river_ids[river_ids > 0]
        if len(river_ids) == 0:
            return

        river_id_counts = [
            (int(river_id), int(np.sum(river_component_labels[yy_cells, xx_cells] == river_id)))
            for river_id in river_ids
        ]
        section_river_id = max(river_id_counts, key=lambda rec: rec[1])[0]
        same_section_river = river_component_labels[yy_cells, xx_cells] == section_river_id
        yy_cells = yy_cells[same_section_river]
        xx_cells = xx_cells[same_section_river]
        if len(yy_cells) == 0:
            return

        section_dist = dist_2d[yy_cells, xx_cells]
        finite = np.isfinite(section_dist)
        if not np.any(finite):
            return

        start_cells = [
            (int(iy), int(ix))
            for iy, ix in zip(yy_cells[finite], xx_cells[finite])
        ]
        current_front = set(start_cells)
        visited = set()

        for _ in range(max_steps):
            current_front = {
                cell for cell in current_front
                if cell not in visited
            }
            if not current_front:
                break

            next_front = set()
            for iy, ix in current_front:
                if not (
                    0 <= iy < valid_river_mask.shape[0]
                    and 0 <= ix < valid_river_mask.shape[1]
                    and river_component_labels[iy, ix] == section_river_id
                    and np.isfinite(dist_2d[iy, ix])
                ):
                    continue

                visited.add((iy, ix))
                protected_mask[iy, ix] = True
                current_dist = float(dist_2d[iy, ix])
                downstream_candidates = []

                for dy in (-1, 0, 1):
                    for dx in (-1, 0, 1):
                        if dy == 0 and dx == 0:
                            continue
                        iy_n = iy + dy
                        ix_n = ix + dx
                        if not (
                            0 <= iy_n < valid_river_mask.shape[0]
                            and 0 <= ix_n < valid_river_mask.shape[1]
                            and river_component_labels[iy_n, ix_n] == section_river_id
                            and np.isfinite(dist_2d[iy_n, ix_n])
                        ):
                            continue

                        neighbor_dist = float(dist_2d[iy_n, ix_n])
                        if neighbor_dist < current_dist - min_drop:
                            downstream_candidates.append((neighbor_dist, iy_n, ix_n))

                if not downstream_candidates:
                    continue

                best_next_dist = max(rec[0] for rec in downstream_candidates)
                for neighbor_dist, iy_n, ix_n in downstream_candidates:
                    if neighbor_dist >= best_next_dist - 0.75 * pixel_size_m:
                        next_front.add((int(iy_n), int(ix_n)))

            current_front = next_front

    for grid_cells in matched_boundary_grid_cells:
        trace_downstream_from_section(grid_cells)

    return protected_mask


def _grid_shared_edge_segment(ix, iy, ix_n, iy_n, x_values, y_values, dx, dy):
    x0 = float(x_values[ix])
    y0 = float(y_values[iy])
    x1 = float(x_values[ix_n])
    y1 = float(y_values[iy_n])
    hx = 0.5 * abs(float(dx))
    hy = 0.5 * abs(float(dy))

    if iy_n != iy and ix_n == ix:
        y_edge = 0.5 * (y0 + y1)
        return LineString([(x0 - hx, y_edge), (x0 + hx, y_edge)])
    if ix_n != ix and iy_n == iy:
        x_edge = 0.5 * (x0 + x1)
        return LineString([(x_edge, y0 - hy), (x_edge, y0 + hy)])
    raise ValueError("Only 4-neighbour grid edges are supported.")


def _split_long_lines(geoms, max_points=4000):
    out_geoms = []
    for geom in geoms:
        if geom.is_empty or geom.geom_type != "LineString":
            continue
        coords = list(geom.coords)
        if len(coords) <= max_points:
            out_geoms.append(geom)
            continue
        start = 0
        while start < len(coords) - 1:
            end = min(start + max_points, len(coords))
            if end - start >= 2:
                out_geoms.append(LineString(coords[start:end]))
            if end == len(coords):
                break
            start = end - 1
    return out_geoms


def _line_point_count(geom):
    if geom is None or geom.is_empty or geom.geom_type != "LineString":
        return 0
    return len(geom.coords)


def _resample_line_by_point_count(geom, n_points):
    n_points = max(2, int(n_points))
    if geom.length <= 0 or _line_point_count(geom) <= 2:
        return geom
    distances = np.linspace(0.0, geom.length, n_points)
    coords = [geom.interpolate(float(dist)).coords[0] for dist in distances]
    return LineString(coords)


def _limit_weir_gdf_points(gdf, max_points, pixel_size_m):
    if gdf.empty:
        return gdf, 0, 0

    initial_points = int(sum(_line_point_count(geom) for geom in gdf.geometry))
    max_points = int(max_points)
    if initial_points <= max_points:
        return gdf, initial_points, initial_points

    min_length_m = 2.0 * pixel_size_m
    tolerance_m = 1.0 * pixel_size_m
    work = gdf.copy()

    for _ in range(8):
        work["geometry"] = work.geometry.apply(
            lambda geom: geom.simplify(tolerance_m, preserve_topology=False)
        )
        work = work[
            work.geometry.notnull()
            & (~work.geometry.is_empty)
            & (work.geometry.geom_type == "LineString")
            & (work.geometry.length >= min_length_m)
            & (work.geometry.apply(_line_point_count) >= 2)
        ].copy()
        current_points = int(sum(_line_point_count(geom) for geom in work.geometry))
        if current_points <= max_points:
            return work.reset_index(drop=True), initial_points, current_points
        tolerance_m *= 1.5

    work["_length_m"] = work.geometry.length
    work = work.sort_values(by="_length_m", ascending=False).drop(columns="_length_m").copy()
    total_length = float(work.geometry.length.sum())
    if total_length <= 0:
        return work.iloc[0:0].copy(), initial_points, 0

    lengths = work.geometry.length.values.astype(float)
    raw_alloc = np.maximum(2, np.floor(max_points * lengths / total_length).astype(int))
    raw_alloc = np.minimum(raw_alloc, [max(2, _line_point_count(geom)) for geom in work.geometry])

    while raw_alloc.sum() > max_points:
        reducible = np.where(raw_alloc > 2)[0]
        if len(reducible) == 0:
            break
        idx = reducible[np.argmax(raw_alloc[reducible])]
        raw_alloc[idx] -= 1

    if raw_alloc.sum() > max_points:
        keep_count = max(1, max_points // 2)
        work = work.iloc[:keep_count].copy()
        raw_alloc = np.full(len(work), 2, dtype=int)

    work["geometry"] = [
        _resample_line_by_point_count(geom, npts)
        for geom, npts in zip(work.geometry, raw_alloc)
    ]
    final_points = int(sum(_line_point_count(geom) for geom in work.geometry))
    return work.reset_index(drop=True), initial_points, final_points


def build_interface_weir_gdf(
    sf,
    source_mask,
    target_land_mask,
    crest_height_m,
    name_prefix,
    par1=0.6,
):
    source_mask = np.asarray(source_mask, dtype=bool)
    target_land_mask = np.asarray(target_land_mask, dtype=bool)
    if not np.any(source_mask) or not np.any(target_land_mask):
        return gpd.GeoDataFrame(columns=["name", "z", "par1", "geometry"], crs=sf.crs)

    x_values = sf.grid["x"].values
    y_values = sf.grid["y"].values
    dx, dy = sf.grid.raster.res
    nrows, ncols = source_mask.shape

    segments = []
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    yy_land, xx_land = np.where(target_land_mask)
    for iy, ix in zip(yy_land, xx_land):
        for dy_idx, dx_idx in directions:
            iy_n = iy + dy_idx
            ix_n = ix + dx_idx
            if not (0 <= iy_n < nrows and 0 <= ix_n < ncols):
                continue
            if not source_mask[iy_n, ix_n]:
                continue
            segments.append(
                _grid_shared_edge_segment(
                    ix=ix,
                    iy=iy,
                    ix_n=ix_n,
                    iy_n=iy_n,
                    x_values=x_values,
                    y_values=y_values,
                    dx=dx,
                    dy=dy,
                )
            )

    if not segments:
        return gpd.GeoDataFrame(columns=["name", "z", "par1", "geometry"], crs=sf.crs)

    merged = linemerge(unary_union(segments))
    if merged.geom_type == "LineString":
        geoms = [merged]
    else:
        geoms = [
            geom
            for geom in getattr(merged, "geoms", [])
            if geom.geom_type == "LineString" and not geom.is_empty
        ]
    geoms = _split_long_lines(geoms)

    records = []
    for i, geom in enumerate(geoms, start=1):
        records.append({
            "name": f"{name_prefix}_{i:04d}",
            "z": float(crest_height_m),
            "par1": float(par1),
            "geometry": geom,
        })

    return gpd.GeoDataFrame(records, crs=sf.crs)


def apply_sea_dike_as_terrain(
    sf,
    final_river_mask_2d,
    sea_dike_height_m,
):
    current_dep = sf.grid["dep"].values.copy()
    current_msk = sf.grid["msk"].values
    final_river_mask = np.asarray(final_river_mask_2d, dtype=bool)
    open_boundary_mask = current_msk == 2
    land_mask = (current_msk == 1) & (~final_river_mask)

    sea_dike_terrain_mask = (
        binary_dilation(open_boundary_mask, structure=np.ones((3, 3), dtype=bool))
        & land_mask
        & np.isfinite(current_dep)
    )

    if np.any(sea_dike_terrain_mask):
        current_dep[sea_dike_terrain_mask] = np.maximum(
            current_dep[sea_dike_terrain_mask],
            float(sea_dike_height_m),
        )
        sf.grid["dep"] = xr.DataArray(
            current_dep,
            coords=sf.grid["dep"].coords,
            dims=sf.grid["dep"].dims,
        )

    print(
        "   >>> Sea dike applied as terrain: "
        f"{int(sea_dike_terrain_mask.sum())} land-side cell(s), "
        f"crest={float(sea_dike_height_m):.2f} m"
    )
    return sea_dike_terrain_mask


def apply_terrain_dikes_after_terrain(
    sf,
    final_river_mask_2d,
    model_folder,
    sea_dike_height_m,
    river_dike_coast_height_m,
    river_dike_inland_height_m,
    river_dike_gradient_distance_m,
):
    for stale_name in ["sfincs.weir", "dike_weir_lines.geojson"]:
        stale_path = os.path.join(model_folder, stale_name)
        if os.path.exists(stale_path):
            os.remove(stale_path)
    try:
        sf.config.pop("weirfile", None)
    except Exception:
        pass

    current_msk = sf.grid["msk"].values
    current_dep = sf.grid["dep"].values.copy()
    final_river_mask = np.asarray(final_river_mask_2d, dtype=bool)
    land_mask = (current_msk == 1) & (~final_river_mask)
    open_boundary_mask = current_msk == 2

    sea_dike_mask = (
        binary_dilation(open_boundary_mask, structure=np.ones((3, 3), dtype=bool))
        & land_mask
        & np.isfinite(current_dep)
    )
    river_dike_mask = (
        binary_dilation(final_river_mask, structure=np.ones((3, 3), dtype=bool))
        & land_mask
        & np.isfinite(current_dep)
    )

    pixel_size_m = abs(float(sf.grid.raster.res[0]))
    dist_to_open_boundary_m = distance_transform_edt(
        ~open_boundary_mask,
        sampling=(pixel_size_m, pixel_size_m),
    )
    gradient_distance_m = max(float(river_dike_gradient_distance_m), pixel_size_m)
    river_gradient_factor = np.clip(
        dist_to_open_boundary_m / gradient_distance_m,
        0.0,
        1.0,
    )
    river_dike_raise = (
        float(river_dike_coast_height_m)
        + river_gradient_factor
        * (float(river_dike_inland_height_m) - float(river_dike_coast_height_m))
    )

    dike_raise = np.zeros_like(current_dep, dtype=float)
    dike_raise[sea_dike_mask] = np.maximum(
        dike_raise[sea_dike_mask],
        float(sea_dike_height_m),
    )
    dike_raise[river_dike_mask] = np.maximum(
        dike_raise[river_dike_mask],
        river_dike_raise[river_dike_mask],
    )

    dike_mask = dike_raise > 0
    current_dep[dike_mask] = current_dep[dike_mask] + dike_raise[dike_mask]
    sf.grid["dep"] = xr.DataArray(
        current_dep,
        coords=sf.grid["dep"].coords,
        dims=sf.grid["dep"].dims,
    )

    print(
        "   >>> Terrain dikes applied: "
        f"sea_cells={int(sea_dike_mask.sum())}, "
        f"river_cells={int(river_dike_mask.sum())}, "
        f"total_cells={int(dike_mask.sum())}"
    )
    print(
        "   >>> Dike terrain raise: "
        f"sea=+{float(sea_dike_height_m):.2f} m, "
        f"river=+{float(river_dike_coast_height_m):.2f}-"
        f"{float(river_dike_inland_height_m):.2f} m over "
        f"{float(river_dike_gradient_distance_m) / 1000.0:.1f} km"
    )
    if np.any(river_dike_mask):
        river_raise_values = river_dike_raise[river_dike_mask]
        print(
            "   >>> River dike raise stats: "
            f"min=+{float(np.nanmin(river_raise_values)):.2f} m, "
            f"mean=+{float(np.nanmean(river_raise_values)):.2f} m, "
            f"max=+{float(np.nanmax(river_raise_values)):.2f} m"
    )
    return sea_dike_mask, river_dike_mask


def refresh_coastline_waterlevel_boundaries(
    sf,
    coast_da,
    valid_dem_mask_2d,
    low_elevation_land_mask,
    final_river_mask_2d,
):
    """
    Re-apply msk=2 after river processing.
    This keeps river mouths and every coastline-shapefile boundary cell in the
    water-level boundary set used by sfincs.bnd/sfincs.bzs.
    """
    current_msk = sf.grid["msk"].values.copy()
    coast_mask_2d = np.asarray(coast_da.values) > 0
    low_land_mask_2d = np.asarray(low_elevation_land_mask.values) > 0
    final_river_mask = np.asarray(final_river_mask_2d, dtype=bool)

    coastline_boundary_mask = (
        coast_mask_2d
        & np.asarray(valid_dem_mask_2d, dtype=bool)
        & (
            (current_msk > 0)
            | low_land_mask_2d
            | final_river_mask
        )
    )

    previous_boundary_count = int(np.sum(current_msk == 2))
    river_mouth_boundary_count = int(np.sum(coastline_boundary_mask & final_river_mask))
    current_msk[coastline_boundary_mask] = 2

    sf.grid["msk"] = xr.DataArray(
        current_msk,
        coords=sf.grid["msk"].coords,
        dims=sf.grid["msk"].dims,
    )

    yy_bnd, xx_bnd = np.where(current_msk == 2)
    x_bnd = sf.grid["x"].values[xx_bnd]
    y_bnd = sf.grid["y"].values[yy_bnd]
    bnd_points = [Point(x, y) for x, y in zip(x_bnd, y_bnd)]
    bnd_gdf = gpd.GeoDataFrame(geometry=bnd_points, crs=sf.crs)
    bnd_gdf.index = range(1, len(bnd_gdf) + 1)
    sf.geoms["bnd"] = bnd_gdf

    print(
        "   >>> Coastline water-level boundaries refreshed: "
        f"previous={previous_boundary_count}, final={len(bnd_gdf)}, "
        f"river_mouth_cells={river_mouth_boundary_count}"
    )

    return bnd_gdf, x_bnd, y_bnd, coastline_boundary_mask


def setup_landcover_manning_after_terrain(
    sf,
    landcover_dir,
    manning_map,
    below_sea_level_manning,
    default_land_manning,
    landcover_nodata_value,
):
    landcover_files = sorted(glob.glob(os.path.join(landcover_dir, "*.tif")))
    if not landcover_files:
        raise FileNotFoundError(f"No landcover tif found in: {landcover_dir}")
    landcover_path = landcover_files[0]

    out_shape = (sf.grid.sizes["y"], sf.grid.sizes["x"])
    lc_on_grid = np.full(
        out_shape,
        int(landcover_nodata_value),
        dtype=np.uint16,
    )

    with rasterio.open(landcover_path) as src:
        src_nodata = src.nodata
        if src_nodata is None:
            src_nodata = landcover_nodata_value
        reproject(
            source=rasterio.band(src, 1),
            destination=lc_on_grid,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src_nodata,
            dst_transform=sf.grid.raster.transform,
            dst_crs=sf.crs,
            dst_nodata=landcover_nodata_value,
            resampling=Resampling.nearest,
        )

    current_dep = sf.grid["dep"].values
    current_msk = sf.grid["msk"].values
    active_mask = current_msk > 0

    manning = np.full(
        out_shape,
        -9999.0,
        dtype=np.float32,
    )
    manning[active_mask] = float(default_land_manning)

    for lc_code, manning_value in manning_map.items():
        code_mask = active_mask & (lc_on_grid == int(lc_code))
        if np.any(code_mask):
            manning[code_mask] = float(manning_value)

    below_sea_mask = active_mask & np.isfinite(current_dep) & (current_dep < 0.0)
    manning[below_sea_mask] = float(below_sea_level_manning)

    manning_da = xr.DataArray(
        manning,
        coords=sf.grid["dep"].coords,
        dims=sf.grid["dep"].dims,
        name="manning",
    )
    manning_da.attrs.update({
        "standard_name": "manning roughness",
        "unit": "s.m-1/3",
    })
    sf.grid["manning"] = manning_da
    sf.set_config("manningfile", "sfincs.man")
    for key in ["manning", "manning_land", "manning_sea", "rgh_lev_land"]:
        sf.config.pop(key, None)

    used_codes, used_counts = np.unique(lc_on_grid[active_mask], return_counts=True)
    mapped_codes = set(int(code) for code in manning_map.keys())
    unmapped_codes = [
        int(code)
        for code in used_codes
        if int(code) not in mapped_codes and int(code) != int(landcover_nodata_value)
    ]

    print(
        "   >>> Spatial Manning roughness generated from landcover: "
        f"{landcover_path}"
    )
    print(
        "   >>> Manning active cells: "
        f"active={int(active_mask.sum())}, "
        f"below_sea={int(below_sea_mask.sum())}, "
        f"unmapped_lc_codes={unmapped_codes[:20]}"
    )
    print(
        "   >>> Manning stats: "
        f"min={float(np.nanmin(manning[active_mask])):.3f}, "
        f"mean={float(np.nanmean(manning[active_mask])):.3f}, "
        f"max={float(np.nanmax(manning[active_mask])):.3f}"
    )
    return manning_da


def reproject_mdt_to_sfincs_grid(
    sf,
    mdt_path,
    mdt_variable,
    fallback_offset_m,
):
    if not os.path.exists(mdt_path):
        raise FileNotFoundError(f"MDT file not found: {mdt_path}")

    ds_mdt = xr.open_dataset(mdt_path)
    if mdt_variable not in ds_mdt:
        raise KeyError(
            f"Variable '{mdt_variable}' not found in MDT file. "
            f"Available variables: {list(ds_mdt.data_vars)}"
        )

    da_mdt = ds_mdt[mdt_variable]
    if "time" in da_mdt.dims:
        da_mdt = da_mdt.isel(time=0)

    lon_name = "longitude" if "longitude" in da_mdt.coords else "lon"
    lat_name = "latitude" if "latitude" in da_mdt.coords else "lat"
    lons = np.asarray(da_mdt[lon_name].values, dtype=float)
    lats = np.asarray(da_mdt[lat_name].values, dtype=float)

    if lons[0] > lons[-1]:
        da_mdt = da_mdt.sortby(lon_name)
        lons = np.asarray(da_mdt[lon_name].values, dtype=float)
    if lats[0] > lats[-1]:
        da_mdt = da_mdt.sortby(lat_name)
        lats = np.asarray(da_mdt[lat_name].values, dtype=float)

    dx = float(np.nanmedian(np.diff(lons)))
    dy = float(np.nanmedian(np.diff(lats)))
    left = float(np.nanmin(lons) - 0.5 * abs(dx))
    top = float(np.nanmax(lats) + 0.5 * abs(dy))
    src_transform = from_origin(left, top, abs(dx), abs(dy))

    src_data = np.asarray(da_mdt.values, dtype=np.float32)
    # rasterio expects row 0 at the north side.
    src_data = src_data[::-1, :]

    out_shape = (sf.grid.sizes["y"], sf.grid.sizes["x"])
    mdt_on_grid = np.full(out_shape, np.nan, dtype=np.float32)

    reproject(
        source=src_data,
        destination=mdt_on_grid,
        src_transform=src_transform,
        src_crs="EPSG:4326",
        src_nodata=np.nan,
        dst_transform=sf.grid.raster.transform,
        dst_crs=sf.crs,
        dst_nodata=np.nan,
        resampling=Resampling.bilinear,
    )

    valid_mdt = np.isfinite(mdt_on_grid)
    if np.any(valid_mdt):
        _, nearest_idx = distance_transform_edt(
            ~valid_mdt,
            return_indices=True,
        )
        mdt_on_grid = np.where(
            valid_mdt,
            mdt_on_grid,
            mdt_on_grid[nearest_idx[0], nearest_idx[1]],
        )
    else:
        mdt_on_grid[:, :] = float(fallback_offset_m)
        print(
            "   [WARNING] MDT has no valid values on the SFINCS grid; "
            f"fallback constant offset {float(fallback_offset_m):.3f} m is used."
        )

    return xr.DataArray(
        mdt_on_grid.astype(np.float32),
        coords=sf.grid["dep"].coords,
        dims=sf.grid["dep"].dims,
        name="local_msl_offset",
    )


def apply_fabdem_to_local_msl_correction(
    sf,
    use_mdt_grid,
    mdt_path,
    mdt_variable,
    constant_offset_m,
):
    if use_mdt_grid:
        local_msl_offset = reproject_mdt_to_sfincs_grid(
            sf=sf,
            mdt_path=mdt_path,
            mdt_variable=mdt_variable,
            fallback_offset_m=constant_offset_m,
        )
    else:
        local_msl_offset = xr.full_like(
            sf.grid["dep"],
            float(constant_offset_m),
            dtype=np.float32,
        )
        local_msl_offset.name = "local_msl_offset"

    dep_values = sf.grid["dep"].values.copy()
    offset_values = local_msl_offset.values
    valid_dep = np.isfinite(dep_values)
    dep_values[valid_dep] = dep_values[valid_dep] - offset_values[valid_dep]

    sf.grid["dep"] = xr.DataArray(
        dep_values,
        coords=sf.grid["dep"].coords,
        dims=sf.grid["dep"].dims,
    )
    sf.grid["local_msl_offset"] = local_msl_offset

    valid_offset = np.isfinite(offset_values[valid_dep])
    if np.any(valid_offset):
        vals = offset_values[valid_dep][valid_offset]
        print(
            "   >>> FABDEM vertical datum corrected to local MSL: "
            f"dep = dep - MDT, offset min={float(np.nanmin(vals)):.3f} m, "
            f"mean={float(np.nanmean(vals)):.3f} m, "
            f"max={float(np.nanmax(vals)):.3f} m"
        )
    else:
        print(
            "   >>> FABDEM vertical datum corrected to local MSL "
            f"with constant offset {float(constant_offset_m):.3f} m"
        )

    return local_msl_offset


CAMA_CFG = get_cama_config(CAMA_RES_TAG, CAMA_BASE_DIR)
CAMA_DIR = CAMA_CFG["dir"]
CAMA_UPAREA_BIN = CAMA_CFG["uparea_bin"]
CAMA_NEXTXY_BIN = CAMA_CFG["nextxy_bin"]
CAMA_WIDTH_BIN = CAMA_CFG.get("width_bin")

# =========================================================
# 1-3. 初始化网格与加载地形
# =========================================================
print("\n[STEP 1-3] Initializing Grid, Topo and Coastal Boundaries...")
sf = SfincsModel(root=MODEL_FOLDER, mode="w+")
sf.setup_config(tref=TIME_START, tstart=TIME_START, tstop=TIME_STOP, dtmaxout=3600, zsini=0.0)

sf.setup_grid_from_region(region={"bbox": MODEL_BBOX}, res=200, crs=AUTO_EPSG)

gdf = gpd.read_file(LAND_POLY_SHP, bbox=tuple(MODEL_BBOX)).to_crs(sf.crs)
gdf["geometry"] = gdf.geometry.buffer(0)
land_gdf = gdf[gdf.geometry.type.isin(["Polygon", "MultiPolygon"])]

fabdem_full = rioxarray.open_rasterio(FABDEM_PATH)
minx, miny, maxx, maxy = MODEL_BBOX
fabdem_local = fabdem_full.rio.clip_box(minx=minx - 0.1, miny=miny - 0.1, maxx=maxx + 0.1, maxy=maxy + 0.1)
if fabdem_local.ndim == 3:
    fabdem_local = fabdem_local.isel(band=0)

fabdem_clean = fabdem_local.where((fabdem_local > -50.0) & (fabdem_local != 0))
fabdem_clean.name = "elevtn"
ds_topo = fabdem_clean.to_dataset()
ds_topo.rio.write_crs("EPSG:4326", inplace=True)

sf.setup_dep(datasets_dep=[ds_topo])
if USE_LOCAL_MSL_CORRECTION:
    apply_fabdem_to_local_msl_correction(
        sf=sf,
        use_mdt_grid=USE_MDT_GRID_MSL_CORRECTION,
        mdt_path=MDT_PATH,
        mdt_variable=MDT_VARIABLE,
        constant_offset_m=LOCAL_MSL_OFFSET_M,
    )
original_fabdem_dep = sf.grid["dep"].copy()

# =========================================================
# 4. 生成计算掩膜与边界
# =========================================================
print("\n[STEP 4/9] Generating Boundaries & Eradicating Void Artifacts...")
sf.setup_mask_active(mask=land_gdf, reset_mask=True)

valid_dem_mask_2d = ~np.isnan(sf.grid["dep"].values)

sf.grid["msk"] = xr.where(valid_dem_mask_2d, sf.grid["msk"], 0)

original_land_mask = sf.grid["msk"].copy()
sf.grid["msk"] = xr.where(sf.grid["dep"] > 15.0, 0, sf.grid["msk"])
low_elevation_land_mask = sf.grid["msk"].copy()

coastline_geom = land_gdf.unary_union.boundary
transform = sf.grid.raster.transform
out_shape = (sf.grid.sizes["y"], sf.grid.sizes["x"])

active_domain_gdf = ACTIVE_DOMAIN_GDF_WGS.to_crs(sf.crs)
active_domain_geom = active_domain_gdf.unary_union
if active_domain_geom is None or active_domain_geom.is_empty:
    raise RuntimeError(f"Empty active partition geometry: {ACTIVE_DOMAIN_GEOJSON}")
active_domain_mask = rasterio.features.rasterize(
    [(active_domain_geom, 1)],
    out_shape=out_shape,
    transform=transform,
    fill=0,
    dtype=np.uint8,
    all_touched=True,
)
active_domain_da = xr.DataArray(
    active_domain_mask,
    coords=sf.grid["msk"].coords,
    dims=sf.grid["msk"].dims,
)
original_land_mask = xr.where(active_domain_da == 1, original_land_mask, 0)
low_elevation_land_mask = xr.where(active_domain_da == 1, low_elevation_land_mask, 0)
sf.grid["msk"] = xr.where(active_domain_da == 1, sf.grid["msk"], 0)
partition_boundary_mask = rasterio.features.rasterize(
    [(active_domain_geom.boundary, 1)],
    out_shape=out_shape,
    transform=transform,
    fill=0,
    dtype=np.uint8,
    all_touched=True,
)
partition_boundary_band = binary_dilation(
    partition_boundary_mask.astype(bool),
    iterations=BOX_EDGE_CHECK_CELLS,
)
print(f"   >>> Active partition mask: {ACTIVE_DOMAIN_GEOJSON}")

coast_mask = rasterio.features.rasterize(
    [(coastline_geom, 1)],
    out_shape=out_shape,
    transform=transform,
    fill=0,
    dtype=np.uint8,
    all_touched=True,
)
coast_da = xr.DataArray(coast_mask, coords=sf.grid["msk"].coords, dims=sf.grid["msk"].dims)

is_boundary = (coast_da == 1) & (sf.grid["msk"] == 1)
sf.grid["msk"] = xr.where(is_boundary, 2, sf.grid["msk"])

y_idx, x_idx = np.where(is_boundary.values)
x_coords = sf.grid["x"].values[x_idx]
y_coords = sf.grid["y"].values[y_idx]

points = [Point(x, y) for x, y in zip(x_coords, y_coords)]
bnd_gdf = gpd.GeoDataFrame(geometry=points, crs=sf.crs)
bnd_gdf.index += 1
sf.geoms["bnd"] = bnd_gdf
print(f"   [OK] Forcefully created {len(bnd_gdf)} boundary points within valid DEM extent.")

# =========================================================
# 5. 跨洲河网检索 + 拓扑合并 + 智能扩宽
# =========================================================
print("\n[STEP 5/9] Global River Database Scanning & Topological Parsing...")

grit_files = glob.glob(os.path.join(RIVER_DIR, "GRITv06_reaches_*_EPSG4326.gpkg"))
final_river_mask_2d = np.zeros(out_shape, dtype=bool)
main_river_weir_mask_2d = np.zeros(out_shape, dtype=bool)

if grit_files:
    gdf_list = []
    for gpkg_file in grit_files:
        try:
            chunk_gdf = gpd.read_file(gpkg_file, layer="lines", bbox=tuple(MODEL_BBOX))
            if not chunk_gdf.empty:
                gdf_list.append(chunk_gdf)
        except Exception:
            pass

    if gdf_list:
        rivers_gdf = pd.concat(gdf_list, ignore_index=True).to_crs(sf.crs)
        rivers_gdf = rivers_gdf.explode(index_parts=False).reset_index(drop=True)

        pixel_size_m = abs(sf.grid.raster.res[0])
        cfd_draw_width_m = CFD_RESOLVE_CELLS * pixel_size_m
        print(f"   >>> SFINCS Grid Resolution: {pixel_size_m:.1f} m")
        print(f"   >>> Physical Drop Threshold after CaMa matching: < {MIN_PHYSICAL_WIDTH_M:.1f} m")
        print(f"   >>> CFD Dilation Enforced Width: {cfd_draw_width_m:.1f} m")

        print("   >>> Rebuilding complete physical reaches from topological nodes...")
        node_to_segs = defaultdict(list)
        seg_to_nodes = {}

        for idx, row in rivers_gdf.iterrows():
            geom = row.geometry
            if geom is None or geom.is_empty or geom.geom_type != "LineString":
                continue
            coords = list(geom.coords)
            start_node = (round(coords[0][0], 0), round(coords[0][1], 0))
            end_node   = (round(coords[-1][0], 0), round(coords[-1][1], 0))
            node_to_segs[start_node].append(idx)
            node_to_segs[end_node].append(idx)
            seg_to_nodes[idx] = (start_node, end_node)

        split_nodes = {node for node, segs in node_to_segs.items() if len(segs) != 2}

        visited = set()
        complete_reaches = []

        for split_node in split_nodes:
            for start_seg in node_to_segs[split_node]:
                if start_seg not in visited:
                    reach_segs = []
                    curr_seg = start_seg
                    curr_node = split_node

                    while True:
                        reach_segs.append(curr_seg)
                        visited.add(curr_seg)
                        nodes = seg_to_nodes.get(curr_seg)
                        if not nodes:
                            break
                        next_node = nodes[0] if nodes[1] == curr_node else nodes[1]
                        if next_node in split_nodes:
                            break
                        next_segs = node_to_segs[next_node]
                        curr_seg = next_segs[0] if next_segs[0] != curr_seg else next_segs[1]
                        curr_node = next_node

                    complete_reaches.append(reach_segs)

        for seg in seg_to_nodes.keys():
            if seg not in visited:
                reach_segs = []
                curr_seg = seg
                curr_node = seg_to_nodes[curr_seg][0]
                while curr_seg not in visited:
                    reach_segs.append(curr_seg)
                    visited.add(curr_seg)
                    nodes = seg_to_nodes[curr_seg]
                    next_node = nodes[0] if nodes[1] == curr_node else nodes[1]
                    next_segs = node_to_segs[next_node]
                    curr_seg = next_segs[0] if next_segs[0] != curr_seg else next_segs[1]
                    curr_node = next_node
                complete_reaches.append(reach_segs)

        # ---------------------------------------------------------
        # 评估去留与强制 CFD 拓宽
        # ---------------------------------------------------------
        candidate_seg_indices = []
        rivers_gdf["draw_width"] = np.maximum(
            rivers_gdf["grwl_width_median"].fillna(cfd_draw_width_m),
            cfd_draw_width_m,
        )
        rivers_gdf["width_filter_keep"] = False
        rivers_gdf["reach_avg_width_m"] = np.nan
        kept_count = 0
        dropped_count = 0

        for reach_segs in complete_reaches:
            reach_gdf = rivers_gdf.loc[reach_segs]

            lengths = reach_gdf.geometry.length
            total_length = lengths.sum()
            if total_length == 0:
                continue

            weights = lengths / total_length
            avg_width = (reach_gdf["grwl_width_median"] * weights).sum()

            candidate_seg_indices.extend(reach_segs)
            rivers_gdf.loc[reach_segs, "reach_avg_width_m"] = avg_width
            if avg_width >= MIN_PHYSICAL_WIDTH_M:
                rivers_gdf.loc[reach_segs, "width_filter_keep"] = True
                kept_count += 1
            else:
                dropped_count += 1

        print(
            "   >>> Topology Evaluation before CaMa matching: "
            f"{kept_count} reaches pass width screening, "
            f"{dropped_count} narrow reaches deferred."
        )

        if len(candidate_seg_indices) > 0:
            rivers_gdf = rivers_gdf.loc[candidate_seg_indices].copy()
            rivers_gdf = rivers_gdf.sort_values(by="draw_width", ascending=True)

            raw_depth = 0.083 * (rivers_gdf["grwl_width_median"] ** 0.6)
            rivers_gdf["est_depth"] = np.maximum(raw_depth, 3.0)
            rivers_gdf["match_width"] = (
                rivers_gdf["grwl_width_median"]
                .fillna(rivers_gdf["draw_width"])
                .fillna(cfd_draw_width_m)
                .astype(float)
            )

            rivers_gdf["geometry"] = rivers_gdf.geometry.buffer(rivers_gdf["draw_width"] / 2.0)

            river_depth_raster = rasterio.features.rasterize(
                ((geom, depth) for geom, depth in zip(rivers_gdf.geometry, rivers_gdf["est_depth"])),
                out_shape=out_shape,
                transform=transform,
                fill=0,
                dtype=np.float32,
                all_touched=True,
            )
            river_width_raster = rasterio.features.rasterize(
                ((geom, width) for geom, width in zip(rivers_gdf.geometry, rivers_gdf["match_width"])),
                out_shape=out_shape,
                transform=transform,
                fill=0,
                dtype=np.float32,
                all_touched=True,
            )
            width_keep_shapes = [
                (geom, 1)
                for geom, keep in zip(rivers_gdf.geometry, rivers_gdf["width_filter_keep"])
                if bool(keep)
            ]
            if width_keep_shapes:
                river_width_keep_raster = rasterio.features.rasterize(
                    width_keep_shapes,
                    out_shape=out_shape,
                    transform=transform,
                    fill=0,
                    dtype=np.uint8,
                    all_touched=True,
                )
            else:
                river_width_keep_raster = np.zeros(out_shape, dtype=np.uint8)
            river_depth_da = xr.DataArray(
                river_depth_raster, coords=sf.grid["dep"].coords, dims=sf.grid["dep"].dims
            )

            sf.grid["msk"] = xr.where((river_depth_da > 0) & (original_land_mask == 1), 1, sf.grid["msk"])

            river_mask = (river_depth_raster > 0) & (sf.grid["msk"].values == 1)
            coast_mask_bin = (sf.grid["msk"].values == 2)
            coast_dilated = binary_dilation(coast_mask_bin, iterations=3)
            outlet_mask = river_mask & coast_dilated

            print("   >>> Building river topological routing graph")
            y_idx_r, x_idx_r = np.where(river_mask)
            node_map = np.full(river_mask.shape, -1, dtype=int)
            node_map[y_idx_r, x_idx_r] = np.arange(len(y_idx_r))

            edges, weights = [], []
            dirs = [
                (0, 1, pixel_size_m),
                (1, 0, pixel_size_m),
                (1, 1, pixel_size_m * 1.4142),
                (1, -1, pixel_size_m * 1.4142),
            ]

            for dy, dx, w in dirs:
                y_neighbor = y_idx_r + dy
                x_neighbor = x_idx_r + dx
                valid = (
                    (y_neighbor >= 0)
                    & (y_neighbor < river_mask.shape[0])
                    & (x_neighbor >= 0)
                    & (x_neighbor < river_mask.shape[1])
                )
                valid_neighbors = valid.copy()
                valid_neighbors[valid] = river_mask[y_neighbor[valid], x_neighbor[valid]]

                u = node_map[y_idx_r[valid_neighbors], x_idx_r[valid_neighbors]]
                v = node_map[y_neighbor[valid_neighbors], x_neighbor[valid_neighbors]]

                edges.append(np.column_stack((u, v)))
                edges.append(np.column_stack((v, u)))
                weights.extend([w] * len(u))
                weights.extend([w] * len(u))

            edges = np.vstack(edges)
            weights = np.array(weights)
            graph = csr_matrix((weights, (edges[:, 0], edges[:, 1])), shape=(len(y_idx_r), len(y_idx_r)))
            outlet_nodes = node_map[outlet_mask]

            if len(outlet_nodes) > 0:
                dist_matrix = shortest_path(graph, directed=False, indices=outlet_nodes)
                min_dist = dist_matrix.min(axis=0)
                min_dist[np.isinf(min_dist)] = np.inf

                dist_2d = np.full(river_mask.shape, np.inf, dtype=float)
                dist_2d[y_idx_r, x_idx_r] = min_dist

                valid_river_mask = river_mask & (dist_2d < MAX_RIVER_UPSTREAM_M)
                final_river_mask_2d = valid_river_mask.copy()

                beyond_river = river_mask & (dist_2d >= MAX_RIVER_UPSTREAM_M)
                max_distance_boundary = valid_river_mask & binary_dilation(beyond_river)

                box_edge_band = _edge_band_mask(
                    valid_river_mask.shape,
                    iterations=BOX_EDGE_CHECK_CELLS,
                )
                box_edge_cut_boundary = (
                    valid_river_mask
                    & box_edge_band
                    & (~coast_dilated)
                )

                low_elevation_keep_mask = np.asarray(low_elevation_land_mask.values) > 0
                high_elevation_river = valid_river_mask & (~low_elevation_keep_mask)
                low_elevation_edge_boundary = (
                    valid_river_mask
                    & low_elevation_keep_mask
                    & binary_dilation(
                        high_elevation_river,
                        structure=np.ones((3, 3), dtype=bool),
                    )
                    & (~coast_dilated)
                )
                partition_edge_cut_boundary = (
                    valid_river_mask
                    & partition_boundary_band
                    & (~coast_dilated)
                )

                upstream_boundary_mask = (
                    max_distance_boundary
                    | box_edge_cut_boundary
                    | low_elevation_edge_boundary
                    | partition_edge_cut_boundary
                )
                boundary_source_map = np.zeros(valid_river_mask.shape, dtype=np.uint8)
                boundary_source_map[max_distance_boundary] = 1
                boundary_source_map[box_edge_cut_boundary] = np.where(
                    boundary_source_map[box_edge_cut_boundary] == 1,
                    3,
                    2,
                )
                boundary_source_map[low_elevation_edge_boundary] = np.where(
                    boundary_source_map[low_elevation_edge_boundary] > 0,
                    boundary_source_map[low_elevation_edge_boundary],
                    4,
                )
                boundary_source_map[partition_edge_cut_boundary] = 5
                print(
                    "   >>> SFINCS upstream boundary candidates: "
                    f"max-distance={int(max_distance_boundary.sum())}, "
                    f"box-edge={int(box_edge_cut_boundary.sum())}, "
                    f"low-elevation-edge={int(low_elevation_edge_boundary.sum())}, "
                    f"partition-edge={int(partition_edge_cut_boundary.sum())}, "
                    f"total={int(upstream_boundary_mask.sum())}"
                )

                sf.grid["msk"] = xr.where(river_mask & ~valid_river_mask, 0, sf.grid["msk"])

                dep_before_river_excavation = sf.grid["dep"].copy()
                dist_2d[~valid_river_mask] = 0.0
                dist_da = xr.DataArray(dist_2d, coords=sf.grid["dep"].coords, dims=sf.grid["dep"].dims)
                synthetic_wse_da = dist_da * RIVER_SLOPE

                sf.grid["dep"] = xr.where(
                    valid_river_mask,
                    synthetic_wse_da - river_depth_da,
                    sf.grid["dep"],
                )
                print(f"   >>> Riverbeds excavated with forced {cfd_draw_width_m}m Dilation.")

                # =========================================================
                # SFINCS 连续入流边界 与 CaMa-Flood 匹配逻辑（修正版）
                # 规则：
                # 1) 一个连续边界只能匹配一个 CaMa 像元
                # 2) 一个 CaMa 像元如果匹配到多个连续边界，只保留点数最多的那个
                # =========================================================
                y_src, x_src = np.where(upstream_boundary_mask)
                if len(y_src) > 0:
                    transformer_to_latlon = Transformer.from_crs(sf.crs, "EPSG:4326", always_xy=True)

                    cama_nx = CAMA_CFG["nx"]
                    cama_ny = CAMA_CFG["ny"]
                    cama_res = CAMA_CFG["res_deg"]
                    search_radius_cells = CAMA_CFG["search_radius"]

                    cama_lons = np.linspace(-180 + cama_res / 2, 180 - cama_res / 2, cama_nx)
                    cama_lats = np.linspace(90 - cama_res / 2, -90 + cama_res / 2, cama_ny)

                    if os.path.exists(CAMA_UPAREA_BIN) and os.path.exists(CAMA_NEXTXY_BIN):
                        uparea_global = np.fromfile(CAMA_UPAREA_BIN, dtype="<f4").reshape(cama_ny, cama_nx)
                        nextxy_global = np.fromfile(CAMA_NEXTXY_BIN, dtype="<i4").reshape(2, cama_ny, cama_nx)
                        width_global = None
                        if CAMA_WIDTH_BIN is not None and os.path.exists(CAMA_WIDTH_BIN):
                            width_global = np.fromfile(CAMA_WIDTH_BIN, dtype="<f4").reshape(cama_ny, cama_nx)
                        nextx = nextxy_global[0, :, :] - 1
                        nexty = nextxy_global[1, :, :] - 1

                        # 1) 连通域标记：每个 component 就是一条连续的 SFINCS 入流边界
                        conn_structure = np.ones((3, 3), dtype=np.uint8)
                        boundary_labels, n_components = label(
                            upstream_boundary_mask.astype(np.uint8),
                            structure=conn_structure
                        )

                        print(f"   >>> Found {n_components} continuous upstream boundary component(s).")

                        cama_pre_match_inlets = find_cama_domain_inlets_for_debug(
                            uparea_global=uparea_global,
                            nextx=nextx,
                            nexty=nexty,
                            cama_lons=cama_lons,
                            cama_lats=cama_lats,
                            model_bbox=MODEL_BBOX,
                            uparea_threshold_km2=UPAREA_THRESHOLD_KM2,
                        )
                        write_pre_match_boundary_geojsons(
                            model_folder=MODEL_FOLDER,
                            upstream_boundary_mask=upstream_boundary_mask,
                            boundary_labels=boundary_labels,
                            boundary_source_map=boundary_source_map,
                            valid_river_mask=valid_river_mask,
                            dist_2d=dist_2d,
                            sf=sf,
                            transformer_to_latlon=transformer_to_latlon,
                            cama_inlets=cama_pre_match_inlets,
                            uparea_global=uparea_global,
                            nextx=nextx,
                            nexty=nexty,
                            cama_lons=cama_lons,
                            cama_lats=cama_lats,
                        )

                        component_to_best_cama = {}
                        component_to_points = {}
                        component_to_mean_dist = {}
                        component_to_uparea = {}
                        lonlat_mapping = {}

                        # 2) 每个连通域内部逐点找候选 CaMa，但最后整个连通域只保留一个 CaMa
                        for comp_id in range(1, n_components + 1):
                            comp_mask = (boundary_labels == comp_id)
                            yy_comp, xx_comp = np.where(comp_mask)

                            if len(yy_comp) == 0:
                                continue

                            component_records = []
                            component_points_xy = []

                            for iy_src, ix_src in zip(yy_comp, xx_comp):
                                x_coord = sf.grid["x"].values[ix_src]
                                y_coord = sf.grid["y"].values[iy_src]
                                component_points_xy.append((x_coord, y_coord))

                                lon_pt, lat_pt = transformer_to_latlon.transform(x_coord, y_coord)

                                best_idx, best_up, best_dist = find_best_cama_pixel_for_point(
                                    lon=lon_pt,
                                    lat=lat_pt,
                                    uparea_global=uparea_global,
                                    cama_lons=cama_lons,
                                    cama_lats=cama_lats,
                                    cama_res=cama_res,
                                    cama_nx=cama_nx,
                                    cama_ny=cama_ny,
                                    search_radius_cells=search_radius_cells,
                                    uparea_threshold_km2=UPAREA_THRESHOLD_KM2,
                                )

                                if best_idx is not None:
                                    component_records.append({
                                        "cama_idx": best_idx,
                                        "uparea": best_up,
                                        "dist": best_dist,
                                        "x": x_coord,
                                        "y": y_coord,
                                    })
                                    lonlat_mapping[best_idx] = (cama_lons[best_idx[1]], cama_lats[best_idx[0]])

                            # 这个连续边界完全找不到合适 CaMa，就跳过
                            if len(component_records) == 0:
                                continue

                            best_cama_idx = choose_best_cama_for_component(component_records)

                            dists_best = [rec["dist"] for rec in component_records if rec["cama_idx"] == best_cama_idx]
                            up_best = [rec["uparea"] for rec in component_records if rec["cama_idx"] == best_cama_idx][0]

                            component_to_best_cama[comp_id] = best_cama_idx
                            component_to_points[comp_id] = component_points_xy
                            component_to_mean_dist[comp_id] = float(np.mean(dists_best))
                            component_to_uparea[comp_id] = up_best

                        # 3) 如果同一个 CaMa 像元匹配到多个连续边界，只保留“最宽”的那个（点数最多）
                        cama_to_components = defaultdict(list)
                        for comp_id, cama_idx in component_to_best_cama.items():
                            cama_to_components[cama_idx].append(comp_id)

                        selected_component_for_cama = {}
                        for cama_idx, comp_list in cama_to_components.items():
                            best_comp = None
                            best_score = None
                            for comp_id in comp_list:
                                score = (
                                    len(component_to_points[comp_id]),
                                    -component_to_mean_dist[comp_id],
                                    component_to_uparea[comp_id],
                                )
                                if (best_score is None) or (score > best_score):
                                    best_score = score
                                    best_comp = comp_id
                            selected_component_for_cama[cama_idx] = best_comp

                        # 4) CaMa 拓扑去重：若某个已选 CaMa 在另一个已选 CaMa 下游，则删除
                        nextxy_global = np.fromfile(CAMA_NEXTXY_BIN, dtype="<i4").reshape(2, cama_ny, cama_nx)
                        nextx = nextxy_global[0, :, :] - 1
                        nexty = nextxy_global[1, :, :] - 1

                        selected_cama_keys = list(selected_component_for_cama.keys())
                        sorted_cama_keys = sorted(
                            selected_cama_keys,
                            key=lambda k: component_to_uparea[selected_component_for_cama[k]]
                        )

                        final_valid_cama_keys = []
                        covered_downstream_pixels = set()

                        for curr_y, curr_x in sorted_cama_keys:
                            if (curr_y, curr_x) in covered_downstream_pixels:
                                continue

                            final_valid_cama_keys.append((curr_y, curr_x))

                            trace_y, trace_x = curr_y, curr_x
                            visited_nodes = set()
                            while True:
                                nxt_x = nextx[trace_y, trace_x]
                                nxt_y = nexty[trace_y, trace_x]

                                if nxt_x < 0 or nxt_y < 0:
                                    break
                                if (nxt_y, nxt_x) == (trace_y, trace_x):
                                    break
                                if (nxt_y, nxt_x) in visited_nodes:
                                    break

                                visited_nodes.add((nxt_y, nxt_x))
                                covered_downstream_pixels.add((nxt_y, nxt_x))
                                trace_y, trace_x = nxt_y, nxt_x

                        # Override the old SFINCS-driven matching with a CaMa-driven pass.
                        # Only CaMa cells detected as domain inlets are allowed to create src forcing.
                        component_to_points = {}
                        component_to_grid_cells = {}
                        component_to_source_codes = {}
                        component_to_source_label = {}
                        component_to_edge_sides = {}
                        component_to_width_m = {}
                        component_to_mean_dist_m = {}
                        component_to_max_dist_m = {}
                        for comp_id in range(1, n_components + 1):
                            comp_mask = boundary_labels == comp_id
                            yy_comp, xx_comp = np.where(comp_mask)
                            if len(yy_comp) == 0:
                                continue

                            source_codes_all = boundary_source_map[yy_comp, xx_comp]
                            structural_local = np.isin(source_codes_all, [1, 2, 3, 5])
                            if np.any(structural_local):
                                yy_use = yy_comp[structural_local]
                                xx_use = xx_comp[structural_local]
                                source_codes = source_codes_all[structural_local]
                            else:
                                yy_use = yy_comp
                                xx_use = xx_comp
                                source_codes = source_codes_all

                            component_to_points[comp_id] = [
                                (
                                    float(sf.grid["x"].values[ix_src]),
                                    float(sf.grid["y"].values[iy_src]),
                                )
                                for iy_src, ix_src in zip(yy_use, xx_use)
                            ]
                            component_to_grid_cells[comp_id] = [
                                (int(iy_src), int(ix_src))
                                for iy_src, ix_src in zip(yy_use, xx_use)
                            ]
                            component_to_source_codes[comp_id] = source_codes
                            component_to_source_label[comp_id] = boundary_source_label_from_codes(source_codes)
                            component_to_edge_sides[comp_id] = component_box_edge_sides(
                                yy_use,
                                xx_use,
                                valid_river_mask.shape,
                                BOX_EDGE_CHECK_CELLS,
                            )
                            component_to_width_m[comp_id] = float(np.nanmax(river_width_raster[yy_use, xx_use]))
                            component_to_mean_dist_m[comp_id] = float(np.nanmean(dist_2d[yy_use, xx_use]))
                            component_to_max_dist_m[comp_id] = float(np.nanmax(dist_2d[yy_use, xx_use]))

                        transformer_from_latlon = Transformer.from_crs(
                            "EPSG:4326",
                            sf.crs,
                            always_xy=True,
                        )
                        selected_component_for_cama = {}
                        component_to_uparea = {}
                        component_to_match_dist_m = {}
                        lonlat_mapping = {}
                        cama_match_lonlat_mapping = {}
                        assigned_components = set()
                        cama_match_rows = []
                        cama_match_rows_by_key = {}

                        for candidate_rank, cama_idx in enumerate(cama_pre_match_inlets, start=1):
                            cama_y, cama_x = cama_idx
                            c_lon = float(cama_lons[cama_x])
                            c_lat = float(cama_lats[cama_y])
                            target_info = CAMA_INLET_TARGETS.get(cama_idx, {})
                            if target_info:
                                target_side = target_info["target_bbox_side"]
                                match_lon = float(target_info["match_lon"])
                                match_lat = float(target_info["match_lat"])
                            else:
                                target_side, match_lon, match_lat = project_lonlat_to_bbox_edge(
                                    c_lon,
                                    c_lat,
                                    MODEL_BBOX,
                                )
                            target_x, target_y = transformer_from_latlon.transform(match_lon, match_lat)
                            uparea_km2 = float(uparea_global[cama_y, cama_x] / 1e6)
                            cama_width_m = np.nan
                            if width_global is not None:
                                cama_width_m = float(width_global[cama_y, cama_x])
                                if not np.isfinite(cama_width_m) or cama_width_m <= 0:
                                    cama_width_m = np.nan
                            match_radius_m = float(CAMA_SFINCS_MATCH_MAX_M)
                            if uparea_km2 >= 100000.0:
                                match_radius_m = max(match_radius_m, 80000.0)

                            candidates = []
                            for comp_id, points_xy in component_to_points.items():
                                if comp_id in assigned_components:
                                    continue

                                coords = np.array(points_xy, dtype=float)
                                dist_m = float(np.min(np.sqrt(
                                    (coords[:, 0] - target_x) ** 2
                                    + (coords[:, 1] - target_y) ** 2
                                )))
                                if dist_m > match_radius_m:
                                    continue

                                source_codes = set(
                                    int(code)
                                    for code in np.asarray(component_to_source_codes[comp_id]).ravel()
                                )
                                has_box_edge = (2 in source_codes) or (3 in source_codes)
                                has_structural_cut = any(
                                    code in source_codes
                                    for code in (1, 2, 3, 5)
                                )
                                if not has_structural_cut:
                                    continue
                                sfincs_width_m = float(component_to_width_m.get(comp_id, 0.0))
                                width_ratio = np.inf
                                width_log_error = np.inf
                                if np.isfinite(cama_width_m) and cama_width_m > 0 and sfincs_width_m > 0:
                                    width_ratio = max(sfincs_width_m, cama_width_m) / max(
                                        min(sfincs_width_m, cama_width_m),
                                        1.0,
                                    )
                                    width_log_error = abs(np.log(max(sfincs_width_m, 1.0) / max(cama_width_m, 1.0)))
                                candidates.append({
                                    "component_id": int(comp_id),
                                    "dist_m": float(dist_m),
                                    "npts": int(len(points_xy)),
                                    "width_m": sfincs_width_m,
                                    "cama_width_m": cama_width_m,
                                    "width_ratio": float(width_ratio),
                                    "width_log_error": float(width_log_error),
                                    "has_box_edge": bool(has_box_edge),
                                    "source_rank": 0 if has_structural_cut else 1,
                                    "same_edge_side": target_side in component_to_edge_sides.get(comp_id, set()),
                                })

                            if not candidates:
                                row = {
                                    "status": "unmatched_no_sfincs_section_in_range",
                                    "candidate_rank": int(candidate_rank),
                                    "cama_row": int(cama_y),
                                    "cama_col": int(cama_x),
                                    "cama_lon": c_lon,
                                    "cama_lat": c_lat,
                                    "target_bbox_side": target_side,
                                    "match_lon": float(match_lon),
                                    "match_lat": float(match_lat),
                                    "uparea_km2": uparea_km2,
                                    "cama_width_m": float(cama_width_m) if np.isfinite(cama_width_m) else np.nan,
                                    "match_radius_km": float(match_radius_m / 1000.0),
                                }
                                for target_key in [
                                    "target_method",
                                    "entry_from_row",
                                    "entry_from_col",
                                    "entry_from_lon",
                                    "entry_from_lat",
                                    "entry_to_row",
                                    "entry_to_col",
                                    "entry_to_lon",
                                    "entry_to_lat",
                                ]:
                                    if target_key in target_info:
                                        row[target_key] = target_info[target_key]
                                cama_match_rows.append(row)
                                cama_match_rows_by_key[cama_idx] = row
                                continue

                            # Large rivers may need a wider search radius near irregular
                            # partition boundaries, but the extended ring is only a
                            # fallback. If a section is already available within the
                            # normal radius, keep the original closer-match behavior.
                            normal_radius_candidates = [
                                rec for rec in candidates
                                if rec["dist_m"] <= CAMA_SFINCS_MATCH_MAX_M
                            ]
                            if normal_radius_candidates:
                                candidates = normal_radius_candidates

                            best_source_rank = min(
                                rec["source_rank"]
                                for rec in candidates
                            )
                            candidates = [
                                rec for rec in candidates
                                if rec["source_rank"] == best_source_rank
                            ]

                            nearest_dist_m = min(rec["dist_m"] for rec in candidates)
                            spatial_pool = [
                                rec for rec in candidates
                                if rec["dist_m"] <= nearest_dist_m + CAMA_WIDE_SECTION_EXTRA_DIST_M
                            ]
                            width_similar_pool = [
                                rec for rec in spatial_pool
                                if rec["width_ratio"] <= CAMA_WIDTH_MATCH_MAX_RATIO
                            ]
                            if width_similar_pool:
                                best_match = min(
                                    width_similar_pool,
                                    key=lambda rec: (
                                        rec["width_log_error"],
                                        not rec["same_edge_side"],
                                        rec["dist_m"],
                                        -rec["npts"],
                                    ),
                                )
                            else:
                                same_side_candidates = [
                                    rec for rec in candidates
                                    if rec["same_edge_side"]
                                ]
                                if same_side_candidates:
                                    candidates_for_choice = same_side_candidates
                                else:
                                    box_edge_candidates = [
                                        rec for rec in candidates
                                        if rec["has_box_edge"]
                                    ]
                                    candidates_for_choice = box_edge_candidates or candidates

                                nearest_dist_m = min(
                                    rec["dist_m"]
                                    for rec in candidates_for_choice
                                )
                                width_pool = [
                                    rec for rec in candidates_for_choice
                                    if rec["dist_m"] <= nearest_dist_m + CAMA_WIDE_SECTION_EXTRA_DIST_M
                                ]
                                best_match = max(
                                    width_pool,
                                    key=lambda rec: (
                                        rec["width_m"],
                                        rec["npts"],
                                        -rec["dist_m"],
                                    ),
                                )

                            comp_id = int(best_match["component_id"])
                            match_dist_m = float(best_match["dist_m"])
                            assigned_components.add(comp_id)
                            selected_component_for_cama[cama_idx] = comp_id
                            component_to_uparea[comp_id] = float(uparea_global[cama_y, cama_x])
                            component_to_match_dist_m[comp_id] = float(match_dist_m)
                            lonlat_mapping[cama_idx] = (c_lon, c_lat)
                            cama_match_lonlat_mapping[cama_idx] = (
                                float(match_lon),
                                float(match_lat),
                                target_side,
                            )

                            row = {
                                "status": "matched",
                                "candidate_rank": int(candidate_rank),
                                "cama_row": int(cama_y),
                                "cama_col": int(cama_x),
                                "cama_lon": c_lon,
                                "cama_lat": c_lat,
                                "target_bbox_side": target_side,
                                "match_lon": float(match_lon),
                                "match_lat": float(match_lat),
                                "uparea_km2": uparea_km2,
                                "cama_width_m": float(cama_width_m) if np.isfinite(cama_width_m) else np.nan,
                                "component_id": int(comp_id),
                                "match_dist_km": float(match_dist_m / 1000.0),
                                "boundary_npts": int(len(component_to_points[comp_id])),
                                "boundary_width_m": float(component_to_width_m.get(comp_id, np.nan)),
                                "width_ratio": float(best_match.get("width_ratio", np.nan)),
                                "width_log_error": float(best_match.get("width_log_error", np.nan)),
                                "boundary_source": component_to_source_label[comp_id],
                                "same_edge_side": bool(best_match["same_edge_side"]),
                                "boundary_mean_dist_km": float(component_to_mean_dist_m[comp_id] / 1000.0),
                                "boundary_max_dist_km": float(component_to_max_dist_m[comp_id] / 1000.0),
                            }
                            for target_key in [
                                "target_method",
                                "entry_from_row",
                                "entry_from_col",
                                "entry_from_lon",
                                "entry_from_lat",
                                "entry_to_row",
                                "entry_to_col",
                                "entry_to_lon",
                                "entry_to_lat",
                            ]:
                                if target_key in target_info:
                                    row[target_key] = target_info[target_key]
                            cama_match_rows.append(row)
                            cama_match_rows_by_key[cama_idx] = row

                        sorted_cama_keys = sorted(
                            selected_component_for_cama.keys(),
                            key=lambda key: uparea_global[key[0], key[1]],
                            reverse=True,
                        )
                        final_valid_cama_keys = []
                        pruned_cama_keys = set()
                        for cama_idx in sorted_cama_keys:
                            downstream_kept = cama_flows_to_any(
                                cama_idx,
                                final_valid_cama_keys,
                                nextx,
                                nexty,
                            )
                            if downstream_kept is not None:
                                pruned_cama_keys.add(cama_idx)
                                row = cama_match_rows_by_key.get(cama_idx)
                                if row is not None:
                                    row["status"] = "pruned_by_cama_downstream_match"
                                    row["downstream_kept_row"] = int(downstream_kept[0])
                                    row["downstream_kept_col"] = int(downstream_kept[1])
                                continue

                            final_valid_cama_keys.append(cama_idx)

                        if cama_match_rows:
                            old_match_summary = os.path.join(
                                MODEL_FOLDER,
                                "cama_sfincs_match_summary.csv",
                            )
                            if os.path.exists(old_match_summary):
                                os.remove(old_match_summary)

                        print(
                            "   >>> CaMa-driven inlet selection: "
                            f"{len(cama_pre_match_inlets)} CaMa candidate(s), "
                            f"{len(final_valid_cama_keys)} matched after CaMa-topology pruning, "
                            f"{len(pruned_cama_keys)} pruned."
                        )

                        # 5) 组织最终 SFINCS src 点和动态流量
                        matched_boundary_grid_cells = [
                            component_to_grid_cells[selected_component_for_cama[cama_idx]]
                            for cama_idx in final_valid_cama_keys
                            if (
                                cama_idx in selected_component_for_cama
                                and selected_component_for_cama[cama_idx] in component_to_grid_cells
                            )
                        ]
                        matched_downstream_river_mask = build_matched_downstream_river_mask(
                            valid_river_mask=valid_river_mask,
                            matched_boundary_grid_cells=matched_boundary_grid_cells,
                            dist_2d=dist_2d,
                            pixel_size_m=pixel_size_m,
                        )
                        main_river_weir_mask_2d = matched_downstream_river_mask.copy()
                        width_screen_keep_mask = river_width_keep_raster > 0
                        post_match_keep_mask = width_screen_keep_mask | matched_downstream_river_mask
                        post_match_removed_mask = valid_river_mask & (~post_match_keep_mask)

                        if np.any(post_match_removed_mask):
                            low_elevation_keep_mask = np.asarray(low_elevation_land_mask.values) > 0
                            post_match_restore_mask = post_match_removed_mask & low_elevation_keep_mask
                            post_match_outside_original_mask = (
                                post_match_removed_mask
                                & (~low_elevation_keep_mask)
                            )
                            valid_river_mask = valid_river_mask & post_match_keep_mask
                            final_river_mask_2d = valid_river_mask.copy()
                            dist_2d[post_match_removed_mask] = 0.0
                            sf.grid["dep"] = xr.where(
                                post_match_restore_mask,
                                dep_before_river_excavation,
                                sf.grid["dep"],
                            )
                            sf.grid["msk"] = xr.where(
                                post_match_outside_original_mask,
                                0,
                                sf.grid["msk"],
                            )
                        else:
                            post_match_restore_mask = post_match_removed_mask
                            post_match_outside_original_mask = post_match_removed_mask

                        screening_pruned_cama_keys = []
                        screened_final_valid_cama_keys = []
                        for cama_idx in final_valid_cama_keys:
                            comp_id = selected_component_for_cama.get(cama_idx)
                            grid_cells = component_to_grid_cells.get(comp_id, [])
                            keeps_boundary = any(
                                valid_river_mask[iy_src, ix_src]
                                for iy_src, ix_src in grid_cells
                            )
                            if keeps_boundary:
                                screened_final_valid_cama_keys.append(cama_idx)
                            else:
                                screening_pruned_cama_keys.append(cama_idx)
                                row = cama_match_rows_by_key.get(cama_idx)
                                if row is not None:
                                    row["status"] = "pruned_by_post_match_width_screening"

                        final_valid_cama_keys = screened_final_valid_cama_keys

                        print(
                            "   >>> Post-match river width screening: "
                            f"{int(matched_downstream_river_mask.sum())} main matched river cell(s) protected, "
                            f"{int(post_match_restore_mask.sum())} unmatched narrow river cell(s) restored to original elevation, "
                            f"{int(post_match_outside_original_mask.sum())} outside original <=15m domain removed from mask, "
                            f"{len(screening_pruned_cama_keys)} matched inlet(s) pruned by normal width screening."
                        )

                        valid_src_points = []
                        q_vals_dynamic = []
                        river_features = []
                        final_src_features = []
                        inlet_source_counts = []
                        boundary_coverage_records = []
                        upstream_graph = defaultdict(list)

                        valid_trace_mask = uparea_global > 100 * 1e2
                        y_valid, x_valid = np.where(valid_trace_mask)
                        for y_v, x_v in zip(y_valid, x_valid):
                            nx_idx, ny_idx = nextx[y_v, x_v], nexty[y_v, x_v]
                            if nx_idx >= 0 and ny_idx >= 0:
                                upstream_graph[(ny_idx, nx_idx)].append((y_v, x_v))

                        for inlet_id, cama_idx in enumerate(final_valid_cama_keys, start=1):
                            comp_id = selected_component_for_cama[cama_idx]
                            uparea_km2 = component_to_uparea[comp_id] / 1e6
                            c_lon, c_lat = lonlat_mapping[cama_idx]
                            match_lon, match_lat, target_side = cama_match_lonlat_mapping.get(
                                cama_idx,
                                (c_lon, c_lat, "unknown"),
                            )
                            boundary_source = component_to_source_label.get(comp_id, "unknown")
                            boundary_width_m = component_to_width_m.get(comp_id, np.nan)
                            line_cells = select_single_row_boundary_cells(
                                component_to_grid_cells[comp_id],
                                target_side,
                                valid_river_mask.shape,
                            )
                            target_x, target_y = transformer_from_latlon.transform(match_lon, match_lat)
                            line_cells = select_nearest_contiguous_run(
                                line_cells,
                                target_x,
                                target_y,
                                sf,
                                target_side,
                            )
                            if not line_cells:
                                line_cells = component_to_grid_cells[comp_id]
                            if boundary_source in {"low_elevation_edge", "partition_edge"}:
                                width_for_limit_m = boundary_width_m
                                if not np.isfinite(width_for_limit_m) or width_for_limit_m <= 0:
                                    width_for_limit_m = 2.0 * pixel_size_m
                                max_section_npts = max(
                                    2,
                                    int(np.ceil(float(width_for_limit_m) / max(float(pixel_size_m), 1.0))) + 4,
                                )
                                if len(line_cells) > max_section_npts:
                                    cell_coords = np.array([
                                        (
                                            float(sf.grid["x"].values[ix_src]),
                                            float(sf.grid["y"].values[iy_src]),
                                        )
                                        for iy_src, ix_src in line_cells
                                    ])
                                    dists_to_target = np.sqrt(
                                        (cell_coords[:, 0] - target_x) ** 2
                                        + (cell_coords[:, 1] - target_y) ** 2
                                    )
                                    keep_idx = np.argsort(dists_to_target)[:max_section_npts]
                                    line_cells = [line_cells[int(idx)] for idx in keep_idx]
                                    if target_side in ("north", "south"):
                                        line_cells = sorted(line_cells, key=lambda cell: cell[1])
                                    elif target_side in ("west", "east"):
                                        line_cells = sorted(line_cells, key=lambda cell: cell[0])
                                    else:
                                        line_cells = sorted(line_cells)
                                coverage_record = {
                                    "front_npts": int(len(component_to_grid_cells[comp_id])),
                                    "initial_npts": int(len(line_cells)),
                                    "final_npts": int(len(line_cells)),
                                    "completed": False,
                                }
                            else:
                                line_cells, coverage_record = complete_inflow_cells_with_matched_boundary_front(
                                    component_grid_cells=component_to_grid_cells[comp_id],
                                    initial_line_cells=line_cells,
                                    target_x=target_x,
                                    target_y=target_y,
                                    sf=sf,
                                    target_side=target_side,
                                    shape=valid_river_mask.shape,
                                )
                            if not line_cells:
                                line_cells = component_to_grid_cells[comp_id]
                                coverage_record["final_npts"] = int(len(line_cells))
                            coverage_record.update({
                                "inlet_id": int(inlet_id),
                                "component_id": int(comp_id),
                            })
                            boundary_coverage_records.append(coverage_record)

                            pts_xy = [
                                (
                                    float(sf.grid["x"].values[ix_src]),
                                    float(sf.grid["y"].values[iy_src]),
                                )
                                for iy_src, ix_src in line_cells
                            ]
                            inlet_source_counts.append(len(pts_xy))

                            cama_cell_lon = float(cama_lons[cama_idx[1]])
                            cama_cell_lat = float(cama_lats[cama_idx[0]])
                            match_dist_km = component_to_match_dist_m.get(comp_id, np.nan) / 1000.0
                            candidate_rank = int(cama_pre_match_inlets.index(cama_idx) + 1)

                            total_q = max(uparea_km2 * SPECIFIC_YIELD, 100.0)
                            q_per_point = total_q / len(pts_xy)

                            pts_list = [Point(xy[0], xy[1]) for xy in pts_xy]
                            src_col_start = len(valid_src_points) + 1
                            src_col_end = src_col_start + len(pts_list) - 1
                            src_cols = ",".join(str(col) for col in range(src_col_start, src_col_end + 1))
                            coverage_props = {
                                "sfincs_boundary_front_npts": int(coverage_record["front_npts"]),
                                "sfincs_initial_run_npts": int(coverage_record["initial_npts"]),
                                "sfincs_boundary_completed": bool(coverage_record["completed"]),
                            }

                            river_features.append({
                                "geometry": Point(c_lon, c_lat),
                                "feature_type": "CaMa_Inlet_Point",
                                "geometry_role": "camaflood_downstream_inflow_point",
                                "inlet_id": int(inlet_id),
                                "cama_boundary_id": int(candidate_rank),
                                "candidate_rank": int(candidate_rank),
                                "cama_row": int(cama_idx[0]),
                                "cama_col": int(cama_idx[1]),
                                "cama_lon": float(cama_cell_lon),
                                "cama_lat": float(cama_cell_lat),
                                "match_lon": float(match_lon),
                                "match_lat": float(match_lat),
                                "target_bbox_side": target_side,
                                "uparea_km2": uparea_km2,
                                "component_id": int(comp_id),
                                "boundary_npts": int(len(pts_xy)),
                                "boundary_source": boundary_source,
                                "boundary_width_m": float(boundary_width_m),
                                "sfincs_src_col_start": int(src_col_start),
                                "sfincs_src_col_end": int(src_col_end),
                                "sfincs_src_cols": src_cols,
                                "sfincs_src_npts": int(len(pts_list)),
                                **coverage_props,
                                "match_dist_km": float(match_dist_km),
                                "q_total_m3s": float(total_q),
                            })

                            section_lonlat = _ordered_lonlat_from_xy(pts_xy, transformer_to_latlon)
                            river_features.append({
                                "geometry": (
                                    LineString(section_lonlat)
                                    if len(section_lonlat) > 1
                                    else Point(section_lonlat[0])
                                ),
                                "feature_type": "SFINCS_Boundary_Section",
                                "geometry_role": "sfincs_inflow_boundary_line",
                                "inlet_id": int(inlet_id),
                                "cama_boundary_id": int(candidate_rank),
                                "candidate_rank": int(candidate_rank),
                                "cama_row": int(cama_idx[0]),
                                "cama_col": int(cama_idx[1]),
                                "cama_lon": float(cama_cell_lon),
                                "cama_lat": float(cama_cell_lat),
                                "uparea_km2": uparea_km2,
                                "component_id": int(comp_id),
                                "boundary_npts": int(len(pts_xy)),
                                "boundary_source": boundary_source,
                                "boundary_width_m": float(boundary_width_m),
                                "sfincs_src_col_start": int(src_col_start),
                                "sfincs_src_col_end": int(src_col_end),
                                "sfincs_src_cols": src_cols,
                                "sfincs_src_npts": int(len(pts_list)),
                                **coverage_props,
                                "match_dist_km": float(match_dist_km),
                                "q_total_m3s": float(total_q),
                            })

                            section_center_lon = float(np.mean([pt[0] for pt in section_lonlat]))
                            section_center_lat = float(np.mean([pt[1] for pt in section_lonlat]))
                            river_features.append({
                                "geometry": LineString([(section_center_lon, section_center_lat), (c_lon, c_lat)]),
                                "feature_type": "Section_Snap_Link",
                                "geometry_role": "section_to_camaflood_link",
                                "inlet_id": int(inlet_id),
                                "cama_boundary_id": int(candidate_rank),
                                "candidate_rank": int(candidate_rank),
                                "cama_row": int(cama_idx[0]),
                                "cama_col": int(cama_idx[1]),
                                "uparea_km2": uparea_km2,
                                "component_id": int(comp_id),
                                "boundary_npts": int(len(pts_xy)),
                                "boundary_source": boundary_source,
                                "boundary_width_m": float(boundary_width_m),
                                "sfincs_src_col_start": int(src_col_start),
                                "sfincs_src_col_end": int(src_col_end),
                                "sfincs_src_cols": src_cols,
                                "sfincs_src_npts": int(len(pts_list)),
                                **coverage_props,
                                "match_dist_km": float(match_dist_km),
                                "q_total_m3s": float(total_q),
                            })

                            for src_point_id, pt in enumerate(pts_list, start=1):
                                sfincs_src_col = src_col_start + src_point_id - 1
                                valid_src_points.append(pt)
                                q_vals_dynamic.append(q_per_point)
                                sfincs_lon, sfincs_lat = transformer_to_latlon.transform(pt.x, pt.y)

                                final_src_features.append({
                                    "geometry": Point(sfincs_lon, sfincs_lat),
                                    "inlet_id": int(inlet_id),
                                    "cama_boundary_id": int(candidate_rank),
                                    "src_point_id": int(src_point_id),
                                    "sfincs_src_col": int(sfincs_src_col),
                                    "candidate_rank": int(candidate_rank),
                                    "cama_row": int(cama_idx[0]),
                                    "cama_col": int(cama_idx[1]),
                                    "cama_lon": float(cama_cell_lon),
                                    "cama_lat": float(cama_cell_lat),
                                    "match_lon": float(match_lon),
                                    "match_lat": float(match_lat),
                                    "target_bbox_side": target_side,
                                    "uparea_km2": uparea_km2,
                                    "component_id": int(comp_id),
                                    "boundary_npts": int(len(pts_xy)),
                                    "boundary_source": boundary_source,
                                    "boundary_width_m": float(boundary_width_m),
                                    **coverage_props,
                                    "match_dist_km": float(match_dist_km),
                                    "q_point_m3s": float(q_per_point),
                                    "q_total_m3s": float(total_q),
                                })

                                river_features.append({
                                    "geometry": Point(sfincs_lon, sfincs_lat),
                                    "feature_type": "SFINCS_Boundary_Point",
                                    "geometry_role": "sfincs_inflow_point",
                                    "inlet_id": int(inlet_id),
                                    "cama_boundary_id": int(candidate_rank),
                                    "src_point_id": int(src_point_id),
                                    "sfincs_src_col": int(sfincs_src_col),
                                    "candidate_rank": int(candidate_rank),
                                    "cama_row": int(cama_idx[0]),
                                    "cama_col": int(cama_idx[1]),
                                    "cama_lon": float(cama_cell_lon),
                                    "cama_lat": float(cama_cell_lat),
                                    "match_lon": float(match_lon),
                                    "match_lat": float(match_lat),
                                    "target_bbox_side": target_side,
                                    "uparea_km2": uparea_km2,
                                    "component_id": int(comp_id),
                                    "boundary_npts": int(len(pts_xy)),
                                    "boundary_source": boundary_source,
                                    "boundary_width_m": float(boundary_width_m),
                                    "sfincs_src_col_start": int(src_col_start),
                                    "sfincs_src_col_end": int(src_col_end),
                                    "sfincs_src_cols": src_cols,
                                    "sfincs_src_npts": int(len(pts_list)),
                                    **coverage_props,
                                    "match_dist_km": float(match_dist_km),
                                    "q_point_m3s": float(q_per_point),
                                    "q_total_m3s": float(total_q),
                                })

                                river_features.append({
                                    "geometry": LineString([(sfincs_lon, sfincs_lat), (c_lon, c_lat)]),
                                    "feature_type": "Snap_Link",
                                    "geometry_role": "sfincs_point_to_camaflood_link",
                                    "inlet_id": int(inlet_id),
                                    "cama_boundary_id": int(candidate_rank),
                                    "src_point_id": int(src_point_id),
                                    "sfincs_src_col": int(sfincs_src_col),
                                    "candidate_rank": int(candidate_rank),
                                    "cama_row": int(cama_idx[0]),
                                    "cama_col": int(cama_idx[1]),
                                    "cama_lon": float(cama_cell_lon),
                                    "cama_lat": float(cama_cell_lat),
                                    "match_lon": float(match_lon),
                                    "match_lat": float(match_lat),
                                    "target_bbox_side": target_side,
                                    "uparea_km2": uparea_km2,
                                    "component_id": int(comp_id),
                                    "boundary_npts": int(len(pts_xy)),
                                    "boundary_source": boundary_source,
                                    "boundary_width_m": float(boundary_width_m),
                                    "sfincs_src_col_start": int(src_col_start),
                                    "sfincs_src_col_end": int(src_col_end),
                                    "sfincs_src_cols": src_cols,
                                    "sfincs_src_npts": int(len(pts_list)),
                                    **coverage_props,
                                    "match_dist_km": float(match_dist_km),
                                    "q_point_m3s": float(q_per_point),
                                    "q_total_m3s": float(total_q),
                                })

                            queue = [cama_idx]
                            visited_cama_trace = set()
                            while queue:
                                curr_y, curr_x = queue.pop(0)
                                if (curr_y, curr_x) in visited_cama_trace:
                                    continue
                                visited_cama_trace.add((curr_y, curr_x))
                                children = upstream_graph.get((curr_y, curr_x), [])
                                for child_y, child_x in children:
                                    lon1, lat1 = cama_lons[child_x], cama_lats[child_y]
                                    lon2, lat2 = cama_lons[curr_x], cama_lats[curr_y]
                                    line = LineString([(lon1, lat1), (lon2, lat2)])
                                    river_features.append({
                                        "geometry": line,
                                        "feature_type": "River_Line",
                                        "geometry_role": "camaflood_upstream_river_line",
                                        "inlet_id": int(inlet_id),
                                        "cama_boundary_id": int(candidate_rank),
                                        "candidate_rank": int(candidate_rank),
                                        "cama_row": int(cama_idx[0]),
                                        "cama_col": int(cama_idx[1]),
                                        "uparea_km2": uparea_global[child_y, child_x] / 1e6,
                                        "component_id": int(comp_id),
                                        "boundary_npts": int(len(pts_xy)),
                                        "boundary_source": boundary_source,
                                        "boundary_width_m": float(boundary_width_m),
                                        "sfincs_src_col_start": int(src_col_start),
                                        "sfincs_src_col_end": int(src_col_end),
                                        "sfincs_src_cols": src_cols,
                                        "sfincs_src_npts": int(len(pts_list)),
                                        **coverage_props,
                                        "match_dist_km": float(match_dist_km),
                                    })
                                    if (child_y, child_x) not in visited_cama_trace:
                                        queue.append((child_y, child_x))

                        if river_features:
                            rivers_gdf_gis = gpd.GeoDataFrame(river_features, crs="EPSG:4326")
                            out_gis_file = os.path.join(MODEL_FOLDER, "cama_upstream_rivers.geojson")
                            rivers_gdf_gis.to_file(out_gis_file, driver="GeoJSON")

                        old_final_src_file = os.path.join(MODEL_FOLDER, "sfincs_final_inflow_points.geojson")
                        if os.path.exists(old_final_src_file):
                            os.remove(old_final_src_file)

                        if len(valid_src_points) > 0:
                            src_gdf = gpd.GeoDataFrame(geometry=valid_src_points, crs=sf.crs)
                            src_gdf.index = range(1, len(src_gdf) + 1)

                            q_vals_arr = np.array(q_vals_dynamic)
                            times = pd.date_range(start=TIME_START, end=TIME_STOP, freq="H")
                            q_data = np.tile(q_vals_arr, (len(times), 1))
                            df_dis = pd.DataFrame(q_data, index=times, columns=src_gdf.index)

                            sf.setup_discharge_forcing(timeseries=df_dis, locations=src_gdf)

                            print(f"   >>> Final matched CaMa inlets: {len(final_valid_cama_keys)}")
                            print(f"   >>> Final continuous SFINCS inflow points: {len(valid_src_points)}")
                            if inlet_source_counts:
                                print(
                                    "   >>> SFINCS source points per inlet: "
                                    f"min={min(inlet_source_counts)}, "
                                    f"mean={np.mean(inlet_source_counts):.1f}, "
                                    f"max={max(inlet_source_counts)}"
                                )
                            if boundary_coverage_records:
                                completed_count = sum(
                                    1
                                    for rec in boundary_coverage_records
                                    if rec["completed"]
                                )
                                print(
                                    "   >>> SFINCS boundary coverage check: "
                                    f"{completed_count}/{len(boundary_coverage_records)} "
                                    "matched boundary section(s) supplemented within their own component."
                                )

        else:
            print("   [WARNING] All rivers dropped due to strict resolution limit.")

# =========================================================
# 5.5 拓扑连通性清理
# =========================================================
print("\n[STEP 5.5] Cleaning up isolated fragmented domains...")
active_cells = sf.grid["msk"].values > 0
labeled_array, num_features = label(active_cells)
boundary_cells = sf.grid["msk"].values == 2
connected_labels = np.unique(labeled_array[boundary_cells])
connected_labels = connected_labels[connected_labels > 0]

protected_labels = []
matched_river_cleanup_mask = active_cells & np.asarray(main_river_weir_mask_2d, dtype=bool)
if np.any(matched_river_cleanup_mask):
    protected_labels.extend(
        int(label_id)
        for label_id in np.unique(labeled_array[matched_river_cleanup_mask])
        if int(label_id) > 0
    )

if valid_src_points:
    x_vals_cleanup = np.asarray(sf.grid["x"].values, dtype=float)
    y_vals_cleanup = np.asarray(sf.grid["y"].values, dtype=float)
    for src_point in valid_src_points:
        ix_src = None
        iy_src = None
        if x_vals_cleanup.size > 0 and y_vals_cleanup.size > 0:
            ix_src = int(np.argmin(np.abs(x_vals_cleanup - float(src_point.x))))
            iy_src = int(np.argmin(np.abs(y_vals_cleanup - float(src_point.y))))
        if (
            ix_src is not None
            and iy_src is not None
            and 0 <= iy_src < labeled_array.shape[0]
            and 0 <= ix_src < labeled_array.shape[1]
        ):
            label_id = int(labeled_array[iy_src, ix_src])
            if label_id > 0:
                protected_labels.append(label_id)

if protected_labels:
    connected_labels = np.unique(
        np.concatenate([
            connected_labels,
            np.asarray(protected_labels, dtype=connected_labels.dtype if connected_labels.size else int),
        ])
    )
    connected_labels = connected_labels[connected_labels > 0]
    print(
        "   >>> Isolated-domain cleanup protection: "
        f"{len(set(protected_labels))} src/matched-river component(s) retained."
    )

if len(connected_labels) > 0:
    is_connected = np.isin(labeled_array, connected_labels)
    isolated_mask = active_cells & ~is_connected
    sf.grid["msk"] = xr.where(isolated_mask, 0, sf.grid["msk"])

# =========================================================
# 5.6 回填内部空洞 + 高程点插值
# =========================================================
print("\n[STEP 5.6] Filling internal voids caused by river-width screening...")
internal_hole_mask = fill_internal_holes_with_elevation_interp(
    sf=sf,
    original_land_mask=original_land_mask,
    valid_dem_mask_2d=valid_dem_mask_2d,
    final_river_mask_2d=final_river_mask_2d,
    original_dep=original_fabdem_dep,
    ring_iterations=3,
)

# =========================================================
# 5.8 基于修补后的地形，进行边界高程继承
# =========================================================
print("\n[STEP 5.8] Post-hole-fill boundary elevation inheritance...")
current_dep = sf.grid["dep"].values
current_msk = sf.grid["msk"].values

# 这里不再只限于 msk==1，凡是 active 区域都允许继承补值
valid_active_mask = (current_msk > 0) & (~np.isnan(current_dep))

if np.any(valid_active_mask):
    invalid_mask = ~valid_active_mask
    dist, indices = distance_transform_edt(invalid_mask, return_indices=True)
    filled_dep = current_dep[indices[0], indices[1]]
    fill_condition = (current_msk > 0) & invalid_mask
    sf.grid["dep"] = xr.where(fill_condition, filled_dep, sf.grid["dep"])

sf.grid["dep"] = sf.grid["dep"].where(sf.grid["msk"] > 0)

# =========================================================
# 5.9 强制抬升所有非河道陆地高程 >= 1.0 米
# =========================================================
print("\n[STEP 5.9] Enforcing minimum land elevation (1.0m) outside rivers...")
is_non_river_land = (sf.grid["msk"].values == 1) & (~final_river_mask_2d)
sf.grid["dep"] = xr.where(is_non_river_land & (sf.grid["dep"] < 1.0), 1.0, sf.grid["dep"])

# =========================================================
# 5.95 强制修改河道地形最大值为 -1.0m
# =========================================================
print("\n[STEP 5.95] Forcing maximum river bed elevation to -1.0m...")
sf.grid["dep"] = xr.where(final_river_mask_2d & (sf.grid["dep"] > -3.0), -3.0, sf.grid["dep"])

# =========================================================
# 5.96 Re-apply coastline water-level boundaries after river processing
# =========================================================
print("\n[STEP 5.96] Refreshing coastline water-level boundaries...")
bnd_gdf, x_coords, y_coords, coastline_boundary_mask = refresh_coastline_waterlevel_boundaries(
    sf=sf,
    coast_da=coast_da,
    valid_dem_mask_2d=valid_dem_mask_2d,
    low_elevation_land_mask=low_elevation_land_mask,
    final_river_mask_2d=final_river_mask_2d,
)

# =========================================================
# 5.97 Add sea and river dikes as terrain after terrain is finalized
# =========================================================
print("\n[STEP 5.97] Applying sea and river dikes as terrain raise...")
apply_terrain_dikes_after_terrain(
    sf=sf,
    final_river_mask_2d=final_river_mask_2d,
    model_folder=MODEL_FOLDER,
    sea_dike_height_m=SEA_DIKE_HEIGHT_M,
    river_dike_coast_height_m=RIVER_DIKE_COAST_HEIGHT_M,
    river_dike_inland_height_m=RIVER_DIKE_INLAND_HEIGHT_M,
    river_dike_gradient_distance_m=RIVER_DIKE_GRADIENT_DISTANCE_M,
)

# =========================================================
# 5.98 Setup spatial Manning roughness from landcover
# =========================================================
print("\n[STEP 5.98] Setting up spatial Manning roughness from landcover...")
setup_landcover_manning_after_terrain(
    sf=sf,
    landcover_dir=LANDCOVER_DIR,
    manning_map=MANNING_MAP,
    below_sea_level_manning=BELOW_SEA_LEVEL_MANNING,
    default_land_manning=DEFAULT_LAND_MANNING,
    landcover_nodata_value=LANDCOVER_NODATA_VALUE,
)

zs_init = np.zeros_like(sf.grid["dep"].values)
sf.grid["zsini"] = xr.DataArray(zs_init, coords=sf.grid["dep"].coords, dims=sf.grid["dep"].dims)

# =========================================================
# 6. 配置风场与挂载原生文件
# =========================================================
print("\n[STEP 6/9] Setting up Wind and Coupling configs...")
dates_1h = pd.date_range(start=TIME_START, end=TIME_STOP, freq="H")
wind_speed = 0 + 0 * np.sin(np.linspace(0, np.pi, len(dates_1h)))
wind_df = pd.DataFrame({"mag": wind_speed, "dir": np.full(len(dates_1h), 90)}, index=dates_1h)
sf.setup_wind_forcing(timeseries=wind_df)
sf.set_config("wind", 1)

sf.set_config("bndfile", "sfincs.bnd")
sf.set_config("bzsfile", "sfincs.bzs")
if "src" in sf.geoms:
    sf.set_config("srcfile", "sfincs.src")
    sf.set_config("disfile", "sfincs.dis")
sf.write()

# =========================================================
# 7. 手动生成 ADCIRC 耦合原生 bnd 和 bzs 文件
# =========================================================
print("\n[STEP 7/9] Generating Native SFINCS Boundary Files...")
num_points = len(x_coords)
bnd_file_path = os.path.join(MODEL_FOLDER, "sfincs.bnd")
with open(bnd_file_path, "w") as f:
    for x, y in zip(x_coords, y_coords):
        f.write(f"{x:.2f} {y:.2f}\n")

dates_10m = pd.date_range(start=TIME_START, end=TIME_STOP, freq="10min")
time_seconds = (dates_10m - dates_10m[0]).total_seconds().values

bzs_file_path = os.path.join(MODEL_FOLDER, "sfincs.bzs")
with open(bzs_file_path, "w") as f:
    for t_sec in time_seconds:
        hours = t_sec / 3600.0
        wl = 0.0 + 0.0 * np.sin(2 * np.pi * hours / 12.42)
        wl_row_string = "    ".join([f"{wl:.3f}"] * num_points)
        f.write(f"{t_sec:8.1f}    {wl_row_string}\n")

if "src" in sf.geoms:
    src_gdf = sf.geoms["src"]
    with open(os.path.join(MODEL_FOLDER, "sfincs.src"), "w") as f:
        for geom in src_gdf.geometry:
            f.write(f"{geom.x:.2f} {geom.y:.2f}\n")

    with open(os.path.join(MODEL_FOLDER, "sfincs.dis"), "w") as f:
        for t_sec in time_seconds:
            q_row_string = "    ".join([f"{q:.1f}" for q in q_vals_dynamic])
            f.write(f"{t_sec:8.1f}    {q_row_string}\n")

print("\n=======================================================")
print(f"GLOBAL ADAPTIVE SUCCESS! Check folder: {MODEL_FOLDER}")
print("=======================================================")
'''


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", nargs="*", default=None, help="Partition numbers to build, e.g. No_001 No_067.")
    parser.add_argument("--limit", type=int, default=None, help="Only build the first N selected partitions.")
    parser.add_argument("--start-rank", type=int, default=None, help="Start at P2 plot_rank >= this value.")
    parser.add_argument("--end-rank", type=int, default=None, help="Stop at P2 plot_rank <= this value.")
    parser.add_argument("--res", type=float, default=None, help="SFINCS grid resolution in metres.")
    parser.add_argument("--cama-res", default=None, help="CaMa-Flood map resolution tag, e.g. 03, 06, or 15.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Directory for generated partition SFINCS models.")
    parser.add_argument("--partition-geojson", type=Path, default=None, help="Partition polygon GeoJSON.")
    parser.add_argument("--domain-csv", type=Path, default=None, help="Partition metadata CSV.")
    parser.add_argument(
        "--hydromt-examples",
        type=Path,
        default=None,
        help="HydroMT-SFINCS examples directory containing the data/ inputs.",
    )
    parser.add_argument(
        "--cama-base-dir",
        type=Path,
        default=None,
        help="CaMa-Flood map root containing glb_<resolution>min.",
    )
    parser.add_argument(
        "--folder-naming",
        choices=("plot-rank", "model-id"),
        default=None,
        help="Name model folders as No_001... or GTC_0001....",
    )
    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=None, help="Overwrite existing model folders.")
    parser.add_argument("--no-overwrite", dest="overwrite", action="store_false", help="Skip partitions with an existing sfincs.inp.")
    parser.add_argument("--legacy-fast", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--dry-run", action="store_true", help="List selected partitions without building.")
    parser.add_argument("--stop-on-error", action="store_true", help="Stop after the first failed partition.")
    parser.add_argument(
        "--plot",
        dest="plot_after_build",
        action="store_true",
        default=True,
        help="Run plot_sfincs_partition_models.py after building selected models.",
    )
    parser.add_argument(
        "--no-plot",
        dest="plot_after_build",
        action="store_false",
        help="Do not run plot_sfincs_partition_models.py after building.",
    )
    return parser.parse_args()


def auto_utm_epsg(lon: float, lat: float) -> int:
    lon = wrap_lon(float(lon))
    zone = int((lon + 180.0) / 6.0) + 1
    zone = max(1, min(60, zone))
    return (32600 if lat >= 0 else 32700) + zone


def wrap_lon(lon: float) -> float:
    return ((lon + 180.0) % 360.0) - 180.0


def clean_geometry(geom):
    if geom is None or geom.is_empty:
        return geom
    try:
        geom = geom.buffer(0)
    except Exception:
        pass
    return geom


def coerce_rank(value) -> int | None:
    try:
        rank = float(value)
    except Exception:
        return None
    if not np.isfinite(rank):
        return None
    return int(round(rank))


def partition_model_folder_name(
    model_id: str, plot_rank=None, folder_naming: str = "plot-rank"
) -> str:
    if folder_naming == "model-id":
        return str(model_id)
    if folder_naming != "plot-rank":
        raise ValueError(f"Unsupported folder naming mode: {folder_naming}")
    rank = coerce_rank(plot_rank)
    if rank is None:
        return str(model_id)
    return f"No_{rank:03d}"


def partition_display_label(model_id: str, plot_rank=None) -> str:
    rank = coerce_rank(plot_rank)
    if rank is None:
        return str(model_id)
    return f"No.{rank}"


def model_id_from_folder_name(name: str) -> str:
    match = re.search(r"(GTC_\d{4})", str(name))
    if match:
        return match.group(1)
    return str(name)


def plot_rank_from_folder_name(name: str) -> int | None:
    match = re.match(r"No_(\d+)(?:_GTC_\d{4})?$", str(name))
    if not match:
        return None
    return coerce_rank(match.group(1))


def require_saved_plot_rank(gdf: gpd.GeoDataFrame, source: Path) -> gpd.GeoDataFrame:
    if "plot_rank" not in gdf.columns:
        raise RuntimeError(
            "Missing 'plot_rank' in partition outputs. "
            "Re-run P2_partition_for_SFINCS_domains.m first so the domain CSV/GeoJSON "
            f"are saved in quicklook order with plot_rank. Checked: {source}"
        )
    out = gdf.copy()
    out["plot_rank"] = pd.to_numeric(out["plot_rank"], errors="coerce")
    if out["plot_rank"].isna().any():
        bad = out.loc[out["plot_rank"].isna(), "model_domain_id"].astype(str).head(5).tolist()
        raise RuntimeError(f"Invalid empty plot_rank values in partition outputs, e.g. {bad}")
    return out


def load_partitions(cfg: Config) -> gpd.GeoDataFrame:
    if not cfg.partition_geojson.exists():
        raise FileNotFoundError(cfg.partition_geojson)
    gdf = gpd.read_file(cfg.partition_geojson).to_crs(4326)
    gdf["geometry"] = gdf.geometry.apply(clean_geometry)
    gdf = gdf[~gdf.geometry.is_empty].copy()
    if cfg.domain_csv.exists():
        meta = pd.read_csv(cfg.domain_csv)
        if "model_domain_id" in meta.columns:
            gdf = gdf.merge(
                meta.drop(columns=[c for c in ["geometry"] if c in meta.columns]),
                on="model_domain_id",
                how="left",
                suffixes=("", "_domain"),
            )
    gdf = require_saved_plot_rank(gdf, cfg.domain_csv if cfg.domain_csv.exists() else cfg.partition_geojson)
    gdf = gdf.sort_values(["plot_rank", "model_domain_id"], kind="mergesort")
    return gdf.reset_index(drop=True)


def select_partitions(gdf: gpd.GeoDataFrame, args: argparse.Namespace) -> gpd.GeoDataFrame:
    out = gdf
    if args.ids:
        ids = {model_id_from_folder_name(v) for v in args.ids if model_id_from_folder_name(v).startswith("GTC_")}
        ranks = {plot_rank_from_folder_name(v) for v in args.ids}
        ranks.discard(None)
        mask = out["model_domain_id"].isin(ids)
        if ranks and "plot_rank" in out.columns:
            mask = mask | pd.to_numeric(out["plot_rank"], errors="coerce").isin(ranks)
        out = out[mask]
    if args.start_rank is not None:
        out = out[pd.to_numeric(out["plot_rank"], errors="coerce") >= args.start_rank]
    if args.end_rank is not None:
        out = out[pd.to_numeric(out["plot_rank"], errors="coerce") <= args.end_rank]
    if args.limit is not None:
        out = out.head(args.limit)
    return out.reset_index(drop=True)


def prepare_partition_domain_files(row, cfg: Config) -> dict:
    model_id = str(row["model_domain_id"])
    plot_rank = row.get("plot_rank", np.nan)
    model_dir = cfg.output_dir / partition_model_folder_name(
        model_id, plot_rank, cfg.folder_naming
    )
    domain_wgs = gpd.GeoDataFrame([row.drop(labels="geometry")], geometry=[row.geometry], crs=4326)
    domain_wgs["geometry"] = domain_wgs.geometry.apply(clean_geometry)
    inp = model_dir / "sfincs.inp"
    if inp.exists() and not cfg.overwrite:
        return {
            "model_domain_id": model_id,
            "plot_rank": plot_rank,
            "status": "skipped_exists",
            "message": "",
            "model_dir": str(model_dir),
            "skip": True,
        }

    if cfg.overwrite and model_dir.exists():
        resolved_out = cfg.output_dir.resolve()
        resolved_model = model_dir.resolve()
        if resolved_out not in resolved_model.parents and resolved_model != resolved_out:
            raise RuntimeError(f"Refusing to overwrite outside output_dir: {model_dir}")
        shutil.rmtree(model_dir)

    model_dir.mkdir(parents=True, exist_ok=True)

    centroid = domain_wgs.geometry.iloc[0].centroid
    epsg = auto_utm_epsg(centroid.x, centroid.y)
    domain_utm = domain_wgs.to_crs(epsg)
    region_wgs = gpd.GeoDataFrame(
        geometry=[clean_geometry(domain_wgs.geometry.iloc[0].envelope)],
        crs=4326,
    )

    active_domain_path = model_dir / "active_domain.geojson"
    region_path = model_dir / "region_grid.geojson"
    domain_wgs.to_file(active_domain_path, driver="GeoJSON")
    region_wgs.to_file(region_path, driver="GeoJSON")
    bbox = tuple(float(v) for v in domain_wgs.total_bounds)

    return {
        "model_domain_id": model_id,
        "plot_rank": plot_rank,
        "model_dir": model_dir,
        "active_domain_path": active_domain_path,
        "region_path": region_path,
        "bbox": bbox,
        "epsg": epsg,
        "skip": False,
    }


def render_embedded_bob_for_partition(cfg: Config, prep: dict) -> str:
    """
    Render the user-provided BoB workflow for one global partition.

    This version intentionally keeps the BoB script unchanged except for the
    study-range/output substitutions required by the partition builder:
      - MODEL_BBOX is set from the current partition extent;
      - MODEL_FOLDER is set to the current partition output directory;
      - the SFINCS grid resolution is taken from --res / Config.res_m.

    No custom river-port, inlet-boundary, CaMa matching, width-scoring, or
    active-domain clipping logic is injected here; those parts remain exactly
    as implemented in BoB_building_15min_fix.py.
    """
    source = EMBEDDED_BOB_TEMPLATE
    bbox = [float(v) for v in prep["bbox"]]
    bbox_literal = "[" + ", ".join(f"{v:.10f}" for v in bbox) + "]"
    model_dir_literal = repr(str(Path(prep["model_dir"]).resolve()).replace("\\", "/"))
    active_domain_literal = repr(str(Path(prep["active_domain_path"]).resolve()).replace("\\", "/"))

    source = re.sub(
        r'^CAMA_RES_TAG\s*=.*$',
        f'CAMA_RES_TAG = "{cfg.cama_res_tag}"',
        source,
        count=1,
        flags=re.MULTILINE,
    )
    if cfg.cama_base_dir is not None:
        cama_literal = repr(str(cfg.cama_base_dir.resolve()).replace("\\", "/"))
        source = re.sub(
            r'^CAMA_BASE_DIR\s*=\s*resolve_cama_base_dir\(.*$',
            f"CAMA_BASE_DIR = {cama_literal}",
            source,
            count=1,
            flags=re.MULTILINE,
        )
    source = re.sub(
        r'^MODEL_FOLDER\s*=.*$',
        f"MODEL_FOLDER = {model_dir_literal}",
        source,
        count=1,
        flags=re.MULTILINE,
    )
    source = re.sub(
        r'^MODEL_BBOX\s*=.*$',
        (
            f"MODEL_BBOX = {bbox_literal}\n"
            f"ACTIVE_DOMAIN_GEOJSON = {active_domain_literal}\n"
            "ACTIVE_DOMAIN_GDF_WGS = gpd.read_file(ACTIVE_DOMAIN_GEOJSON).to_crs('EPSG:4326')\n"
            "ACTIVE_DOMAIN_GDF_WGS['geometry'] = ACTIVE_DOMAIN_GDF_WGS.geometry.buffer(0)\n"
            "ACTIVE_DOMAIN_GDF_WGS = ACTIVE_DOMAIN_GDF_WGS[~ACTIVE_DOMAIN_GDF_WGS.geometry.is_empty]\n"
            "if ACTIVE_DOMAIN_GDF_WGS.empty:\n"
            "    raise RuntimeError(f'Empty active partition geometry: {ACTIVE_DOMAIN_GEOJSON}')\n"
            "ACTIVE_DOMAIN_GEOM_WGS = ACTIVE_DOMAIN_GDF_WGS.unary_union"
        ),
        source,
        count=1,
        flags=re.MULTILINE,
    )

    grid_needle = 'sf.setup_grid_from_region(region={"bbox": MODEL_BBOX}, res=200, crs=AUTO_EPSG)'
    grid_patch = f'sf.setup_grid_from_region(region={{"bbox": MODEL_BBOX}}, res={float(cfg.res_m):.12g}, crs=AUTO_EPSG)'
    if grid_needle not in source:
        raise RuntimeError("BoB template changed: cannot find setup_grid_from_region line.")
    source = source.replace(grid_needle, grid_patch, 1)

    return source


def build_one_partition_with_embedded_bob(row, cfg: Config) -> dict:
    prep = prepare_partition_domain_files(row, cfg)
    if prep.get("skip"):
        return prep

    model_id = prep["model_domain_id"]
    model_dir = Path(prep["model_dir"])
    script_path = model_dir / "_run_embedded_bob_for_partition.py"
    log_path = model_dir / "embedded_bob_build.log"
    rendered = render_embedded_bob_for_partition(cfg, prep)
    script_path.write_text(rendered, encoding="utf-8")

    cmd = [sys.executable, str(script_path.resolve())]
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("Command: " + " ".join(cmd) + "\n")
        log.write(f"Working directory: {cfg.hydromt_examples}\n\n")
        log.flush()
        proc = subprocess.run(
            cmd,
            cwd=str(cfg.hydromt_examples),
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )

    if proc.returncode != 0:
        raise RuntimeError(f"BoB template workflow failed for {model_id}; see {log_path}")

    inp = model_dir / "sfincs.inp"
    if not inp.exists():
        raise RuntimeError(f"BoB template finished but sfincs.inp was not written for {model_id}")

    return {
        "model_domain_id": model_id,
        "plot_rank": row.get("plot_rank", np.nan),
        "status": "ok",
        "message": "",
        "model_dir": str(model_dir),
        "epsg": prep["epsg"],
        "bbox_west": prep["bbox"][0],
        "bbox_south": prep["bbox"][1],
        "bbox_east": prep["bbox"][2],
        "bbox_north": prep["bbox"][3],
        "active_domain": str(prep["active_domain_path"]),
        "embedded_script": str(script_path),
        "embedded_log": str(log_path),
        "rank_domain": row.get("rank_domain", np.nan),
        "basin_id": row.get("basin_id", ""),
    }


def build_one_partition(row, cfg: Config) -> dict:
    return build_one_partition_with_embedded_bob(row, cfg)


def write_status(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8")


def run_plot_after_build(args: argparse.Namespace, cfg: Config, status_rows: list[dict]):
    plot_script = SCRIPT_DIR / "plot_sfincs_partition_models.py"
    if not plot_script.exists():
        raise FileNotFoundError(f"Cannot find plotting script: {plot_script}")

    model_ids: list[str] = []
    for row in status_rows:
        model_id = str(row.get("model_domain_id", ""))
        model_dir = Path(str(row.get("model_dir", cfg.output_dir / model_id)))
        status = str(row.get("status", ""))
        if status in {"ok", "skipped_exists"} and model_id and (model_dir / "sfincs.inp").exists():
            model_ids.append(model_id)

    if not model_ids:
        print("\nNo successfully built or existing SFINCS models to plot.")
        return

    cmd = [
        sys.executable,
        str(plot_script),
        "--models-dir",
        str(cfg.output_dir),
        "--partition-geojson",
        str(cfg.partition_geojson),
        "--domain-csv",
        str(cfg.domain_csv),
        "--land-polygons",
        str(cfg.hydromt_examples / "data" / "coastaline" / "land_polygons.shp"),
        "--membership-csv",
        str(REPO_ROOT / "ADCIRC" / "global_build" / "catalog" / "global_tc_adcirc_block_members.csv"),
        "--ids",
        *model_ids,
    ]
    if args.stop_on_error:
        cmd.append("--stop-on-error")

    print("\n============================================================")
    print("Running plot_sfincs_partition_models.py after model build")
    print(f"Models to plot: {len(model_ids)}")
    print("Command       : " + " ".join(cmd))
    print("============================================================")

    subprocess.run(cmd, cwd=str(SCRIPT_DIR), check=True)


def main():
    args = parse_args()
    cfg = Config()
    if args.partition_geojson is not None:
        cfg.partition_geojson = args.partition_geojson.expanduser().resolve()
    if args.domain_csv is not None:
        cfg.domain_csv = args.domain_csv.expanduser().resolve()
    if args.hydromt_examples is not None:
        cfg.hydromt_examples = args.hydromt_examples.expanduser().resolve()
    if args.cama_base_dir is not None:
        cfg.cama_base_dir = args.cama_base_dir.expanduser().resolve()
    if args.cama_res is not None:
        cfg.cama_res_tag = normalize_cama_res_tag(args.cama_res)
        if args.output_dir is None:
            cfg.output_dir = default_output_dir_for_cama_tag(cfg.cama_res_tag)
            cfg.status_csv = cfg.output_dir / "sfincs_partition_build_status.csv"
    if args.res is not None:
        cfg.res_m = float(args.res)
    if args.output_dir is not None:
        cfg.output_dir = args.output_dir
        cfg.status_csv = cfg.output_dir / "sfincs_partition_build_status.csv"
    if args.folder_naming is not None:
        cfg.folder_naming = args.folder_naming
    if args.overwrite is not None:
        cfg.overwrite = bool(args.overwrite)
    cfg.stop_on_error = bool(args.stop_on_error)
    if args.legacy_fast:
        raise RuntimeError("--legacy-fast is disabled; this script now always uses the full embedded BoB workflow.")

    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    gdf = load_partitions(cfg)
    selected = select_partitions(gdf, args)

    print("============================================================")
    print("Build SFINCS models from TC coastal partition polygons")
    print(f"Partition polygons : {cfg.partition_geojson}")
    print(f"CaMa map tag       : {cfg.cama_res_tag}min")
    print(f"Output dir         : {cfg.output_dir}")
    print(f"Selected partitions: {len(selected)} / {len(gdf)}")
    print(f"Resolution         : {cfg.res_m:.1f} m")
    print(f"Folder naming      : {cfg.folder_naming}")
    print("Workflow           : full embedded BoB workflow")
    print(f"Mode               : {'preview only' if args.dry_run else 'build models'}")
    print("============================================================")

    if args.dry_run:
        cols = [
            c
            for c in ["plot_rank", "rank_domain", "quicklook_basin_id", "basin_id", "grid_cell_count"]
            if c in selected.columns
        ]
        preview = selected[cols].copy()
        preview.insert(
            0,
            "output_folder",
            [
                partition_model_folder_name(model_id, rank, cfg.folder_naming)
                for model_id, rank in zip(selected["model_domain_id"], selected["plot_rank"])
            ],
        )
        print(preview.to_string(index=False))
        return

    if not cfg.hydromt_examples.is_dir():
        raise FileNotFoundError(
            "HydroMT-SFINCS examples directory is missing: "
            f"{cfg.hydromt_examples}. Set --hydromt-examples or "
            "HYDROMT_SFINCS_EXAMPLES."
        )
    if cfg.cama_base_dir is not None and not cfg.cama_base_dir.is_dir():
        raise FileNotFoundError(f"CaMa-Flood map root is missing: {cfg.cama_base_dir}")

    status_rows: list[dict] = []
    for i, (_, row) in enumerate(selected.iterrows(), start=1):
        model_id = str(row["model_domain_id"])
        plot_rank = row.get("plot_rank", np.nan)
        print(f"\n===== [{i}/{len(selected)}] {partition_display_label(model_id, plot_rank)} =====")
        try:
            rec = build_one_partition(row, cfg)
            print(f"  -> {rec['status']}: {rec.get('active_cells', '')} active cells")
        except Exception as exc:
            rec = {
                "model_domain_id": model_id,
                "plot_rank": plot_rank,
                "status": "failed",
                "message": str(exc),
                "traceback": traceback.format_exc(),
                "model_dir": str(
                    cfg.output_dir
                    / partition_model_folder_name(model_id, plot_rank, cfg.folder_naming)
                ),
                "rank_domain": row.get("rank_domain", np.nan),
                "basin_id": row.get("basin_id", ""),
            }
            print(f"  !! failed: {exc}")
            if cfg.stop_on_error:
                status_rows.append(rec)
                write_status(cfg.status_csv, status_rows)
                raise
        status_rows.append(rec)
        write_status(cfg.status_csv, status_rows)

    print(f"\nDone. Status CSV:\n  {cfg.status_csv}")
    if args.plot_after_build:
        run_plot_after_build(args, cfg, status_rows)
    else:
        print("\nSkipping plot step because --no-plot was specified.")


if __name__ == "__main__":
    main()
