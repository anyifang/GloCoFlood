# -*- coding: utf-8 -*-
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
from shapely.ops import linemerge, unary_union
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
MODEL_FOLDER = "Misp_SFINCS"  

MODEL_BBOX = [-93.1,28.7, -88.7, 31.5]
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
MAX_RIVER_UPSTREAM_M = 250000  
BOX_EDGE_CHECK_CELLS = 15
CAMA_SFINCS_MATCH_MAX_M = 50000
CAMA_BOX_EDGE_PRIORITY_M = 50000
CAMA_WIDE_SECTION_EXTRA_DIST_M = 60000

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
    r"H:\Global_compoundflood\Camaflood\cmf_v420_pkg\map",
    r"H:\Global_compoundflood\Camaflood\cmf_v420_pkg\map\glb_15min",
    "/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/cmf_v420_pkg/map",
    "/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/cmf_v420_pkg/map/glb_15min",
]
CAMA_BASE_DIR = next(
    (path for path in CAMA_BASE_DIR_CANDIDATES if os.path.isdir(path)),
    CAMA_BASE_DIR_CANDIDATES[0],
)

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


def find_cama_domain_inlets_for_debug(
    uparea_global,
    nextx,
    nexty,
    cama_lons,
    cama_lats,
    model_bbox,
    uparea_threshold_km2,
):
    west, south, east, north = model_bbox
    up_threshold = uparea_threshold_km2 * 1e6

    lon_inside = (cama_lons >= west) & (cama_lons <= east)
    lat_inside = (cama_lats >= south) & (cama_lats <= north)
    inside = lat_inside[:, None] & lon_inside[None, :]
    valid = uparea_global >= up_threshold

    dlon = abs(cama_lons[1] - cama_lons[0])
    dlat = abs(cama_lats[1] - cama_lats[0])
    pad = 2
    lon_search = (cama_lons >= west - pad * dlon) & (cama_lons <= east + pad * dlon)
    lat_search = (cama_lats >= south - pad * dlat) & (cama_lats <= north + pad * dlat)
    search = lat_search[:, None] & lon_search[None, :]

    inlet_keys = set()
    has_inside_parent = set()

    y_search, x_search = np.where(search & valid)
    for y0, x0 in zip(y_search, x_search):
        y1 = int(nexty[y0, x0])
        x1 = int(nextx[y0, x0])
        if not (0 <= y1 < uparea_global.shape[0] and 0 <= x1 < uparea_global.shape[1]):
            continue

        if inside[y0, x0] and inside[y1, x1]:
            has_inside_parent.add((y1, x1))
        elif (not inside[y0, x0]) and inside[y1, x1]:
            inlet_keys.add((y1, x1))

    y_inside, x_inside = np.where(inside & valid)
    for y0, x0 in zip(y_inside, x_inside):
        y1 = int(nexty[y0, x0])
        x1 = int(nextx[y0, x0])
        downstream_inside = (
            0 <= y1 < uparea_global.shape[0]
            and 0 <= x1 < uparea_global.shape[1]
            and inside[y1, x1]
        )
        near_edge = (
            abs(cama_lons[x0] - west) <= 0.5 * dlon
            or abs(cama_lons[x0] - east) <= 0.5 * dlon
            or abs(cama_lats[y0] - south) <= 0.5 * dlat
            or abs(cama_lats[y0] - north) <= 0.5 * dlat
        )
        if near_edge and downstream_inside and (y0, x0) not in has_inside_parent:
            inlet_keys.add((y0, x0))

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
    return

    sfincs_features = []
    source_labels = {
        1: "max_distance",
        2: "box_edge",
        3: "max_distance_and_box_edge",
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
    if 3 in codes:
        return "max_distance_and_box_edge"
    if 2 in codes:
        return "box_edge"
    if 1 in codes:
        return "max_distance"
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
    Preserve only the main matched SFINCS river.
    The main river is defined by the matched inflow section with the longest
    distance to the outlet, then traced downstream along decreasing dist_2d.
    """
    protected_mask = np.zeros(valid_river_mask.shape, dtype=bool)
    if not matched_boundary_grid_cells or not np.any(valid_river_mask):
        return protected_mask

    conn_structure = np.ones((3, 3), dtype=np.uint8)
    river_component_labels, _ = label(
        valid_river_mask.astype(np.uint8),
        structure=conn_structure,
    )

    main_section = None
    main_score = None
    for grid_cells in matched_boundary_grid_cells:
        if len(grid_cells) == 0:
            continue

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
            continue

        river_ids = np.unique(river_component_labels[yy_cells, xx_cells])
        river_ids = river_ids[river_ids > 0]
        if len(river_ids) == 0:
            continue

        river_id_counts = [
            (int(river_id), int(np.sum(river_component_labels[yy_cells, xx_cells] == river_id)))
            for river_id in river_ids
        ]
        main_river_id = max(river_id_counts, key=lambda rec: rec[1])[0]
        same_section_river = river_component_labels[yy_cells, xx_cells] == main_river_id
        yy_cells = yy_cells[same_section_river]
        xx_cells = xx_cells[same_section_river]
        if len(yy_cells) == 0:
            continue

        section_dist = dist_2d[yy_cells, xx_cells]
        finite = np.isfinite(section_dist)
        if not np.any(finite):
            continue

        section_score = float(np.nanmax(section_dist[finite]))
        if main_score is None or section_score > main_score:
            main_score = section_score
            main_section = (
                [(int(iy), int(ix)) for iy, ix in zip(yy_cells[finite], xx_cells[finite])],
                int(main_river_id),
            )

    if main_section is None:
        return protected_mask

    start_cells, main_river_id = main_section
    current_front = set(start_cells)
    visited = set()
    min_drop = 0.05 * pixel_size_m
    max_steps = max(100, int(np.nanmax(dist_2d[np.isfinite(dist_2d)]) / max(pixel_size_m, 1.0)) + 20)

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
                and river_component_labels[iy, ix] == main_river_id
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
                        and river_component_labels[iy_n, ix_n] == main_river_id
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

                upstream_boundary_mask = max_distance_boundary | box_edge_cut_boundary
                boundary_source_map = np.zeros(valid_river_mask.shape, dtype=np.uint8)
                boundary_source_map[max_distance_boundary] = 1
                boundary_source_map[box_edge_cut_boundary] = np.where(
                    boundary_source_map[box_edge_cut_boundary] == 1,
                    3,
                    2,
                )
                print(
                    "   >>> SFINCS upstream boundary candidates: "
                    f"max-distance={int(max_distance_boundary.sum())}, "
                    f"box-edge={int(box_edge_cut_boundary.sum())}, "
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

                            component_to_points[comp_id] = [
                                (
                                    float(sf.grid["x"].values[ix_src]),
                                    float(sf.grid["y"].values[iy_src]),
                                )
                                for iy_src, ix_src in zip(yy_comp, xx_comp)
                            ]
                            component_to_grid_cells[comp_id] = [
                                (int(iy_src), int(ix_src))
                                for iy_src, ix_src in zip(yy_comp, xx_comp)
                            ]
                            source_codes = boundary_source_map[yy_comp, xx_comp]
                            component_to_source_codes[comp_id] = source_codes
                            component_to_source_label[comp_id] = boundary_source_label_from_codes(source_codes)
                            component_to_edge_sides[comp_id] = component_box_edge_sides(
                                yy_comp,
                                xx_comp,
                                valid_river_mask.shape,
                                BOX_EDGE_CHECK_CELLS,
                            )
                            component_to_width_m[comp_id] = float(np.nanmax(river_width_raster[yy_comp, xx_comp]))
                            component_to_mean_dist_m[comp_id] = float(np.nanmean(dist_2d[yy_comp, xx_comp]))
                            component_to_max_dist_m[comp_id] = float(np.nanmax(dist_2d[yy_comp, xx_comp]))

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
                            target_side, match_lon, match_lat = project_lonlat_to_bbox_edge(
                                c_lon,
                                c_lat,
                                MODEL_BBOX,
                            )
                            target_x, target_y = transformer_from_latlon.transform(match_lon, match_lat)
                            uparea_km2 = float(uparea_global[cama_y, cama_x] / 1e6)

                            candidates = []
                            for comp_id, points_xy in component_to_points.items():
                                if comp_id in assigned_components:
                                    continue

                                coords = np.array(points_xy, dtype=float)
                                dist_m = float(np.min(np.sqrt(
                                    (coords[:, 0] - target_x) ** 2
                                    + (coords[:, 1] - target_y) ** 2
                                )))
                                if dist_m > CAMA_SFINCS_MATCH_MAX_M:
                                    continue

                                source_codes = set(
                                    int(code)
                                    for code in np.asarray(component_to_source_codes[comp_id]).ravel()
                                )
                                has_box_edge = (2 in source_codes) or (3 in source_codes)
                                candidates.append({
                                    "component_id": int(comp_id),
                                    "dist_m": float(dist_m),
                                    "npts": int(len(points_xy)),
                                    "width_m": float(component_to_width_m.get(comp_id, 0.0)),
                                    "has_box_edge": bool(has_box_edge),
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
                                    "match_radius_km": float(CAMA_SFINCS_MATCH_MAX_M / 1000.0),
                                }
                                cama_match_rows.append(row)
                                cama_match_rows_by_key[cama_idx] = row
                                continue

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
                                "component_id": int(comp_id),
                                "match_dist_km": float(match_dist_m / 1000.0),
                                "boundary_npts": int(len(component_to_points[comp_id])),
                                "boundary_width_m": float(component_to_width_m.get(comp_id, np.nan)),
                                "boundary_source": component_to_source_label[comp_id],
                                "same_edge_side": bool(best_match["same_edge_side"]),
                                "boundary_mean_dist_km": float(component_to_mean_dist_m[comp_id] / 1000.0),
                                "boundary_max_dist_km": float(component_to_max_dist_m[comp_id] / 1000.0),
                            }
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
                            boundary_source = component_to_source_label.get(comp_id, "unknown")
                            boundary_width_m = component_to_width_m.get(comp_id, np.nan)
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
