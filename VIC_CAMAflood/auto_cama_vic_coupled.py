# -*- coding: utf-8 -*-
"""
CaMa-Flood + VIC Image Driver 自动耦合建模脚本（NetCDF 版）
================================================================
功能：
1. 建立 CaMa-Flood 区域计算域（reg map）
2. 基于 CaMa 区域图，构建 6min (0.1°) VIC Image Driver 计算域
3. 生成 VIC Image Driver 的 domain.nc / params.nc / global_parameter.txt
4. 支持两种 forcing 模式：
   - READY_VIC_NC: 直接使用已经整理好的 VIC forcing 年度 NetCDF
   - ERA5_DAILY_TO_VIC: 读取你当前按“天”为单位下载的 ERA5 NetCDF（每年一个文件夹），拼接成年尺度 VIC forcing
5. 可选运行 VIC Image Driver
6. 从 VIC 输出 NetCDF 中提取 OUT_RUNOFF + OUT_BASEFLOW，整理成 CaMa-Flood 可读的年度 runoff NetCDF
7. 重新生成 CaMa-Flood 针对 VIC runoff 的 diminfo / inpmat
8. 生成 CaMa-Flood go script（NetCDF runoff 输入）
9. 可选运行 CaMa-Flood

说明：
- 这版改成 VIC Image Driver，核心原因是 Image Driver 原生使用 NetCDF I/O，forcing / domain / params / outputs 都是 NetCDF。
- 你现在已经在 Python 下载 ERA5，并且 CaMa-Flood 也支持 NetCDF runoff，因此这版统一走 NetCDF 链路更顺。
- 当前默认：RUN_VIC = False, RUN_CAMA = False。
- 这版假设你的 ERA5 年度文件已经包含 time/lat/lon 和所需气象变量；如果变量名不同，改 ERA5_VARMAP 即可。
- 植被参数这里先用“单一草地 tile”占位版，以保证框架先跑通。后续你可再接入 MODIS/ESA land cover 做多植被类型参数化。
"""

import os
import re
import glob
import shutil
import subprocess
from typing import Dict, Tuple, Optional, List

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import xarray as xr
from scipy.spatial import cKDTree

print("=======================================================")
print("  🌊 CaMa-Flood + VIC Image Driver AutoPipeline 🌊")
print("=======================================================")

# =========================================================
# 1. 用户配置区
# =========================================================
CAMA_DIR = "../../cmf_v420_pkg"
GEOJSON_PATH = "../../../HydroMT-SFINCS/examples/PRD_SFINCS/cama_upstream_rivers.geojson"

REGION_NAME = "PRD_Auto"
CAMA_RESOLUTION = "03min"   # CaMa 支撑域分辨率
VIC_RESOLUTION = "06min"    # 6min = 0.1°

YEAR_START = 2000
YEAR_END = 2000
BUFFER_DEG = 1.0

# ---------- CaMa-Flood ----------
PREPARE_CAMA_RUN_SCRIPT = True
RUN_CAMA = True

# ---------- VIC Image Driver ----------
VIC_EXEC = "/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/VIC-master/vic/drivers/image/vic_image.exe"
RUN_VIC = True

# VIC Image Driver 使用 NetCDF forcing
# 模式 1: READY_VIC_NC -> 你已经有 VIC 格式 forcing，只需要给 prefix
#   文件应为: <VIC_READY_FORCING_PREFIX><year>.nc
# 模式 2: ERA5_DAILY_TO_VIC -> 读取你当前按“天”为单位下载的 ERA5 文件，拼接成年尺度 VIC forcing
VIC_FORCING_MODE = "ERA5_DAILY_TO_VIC"   # READY_VIC_NC or ERA5_DAILY_TO_VIC or ERA5_TO_VIC

# 如果你已经有 VIC 格式 forcing，填这个前缀
VIC_READY_FORCING_PREFIX = "/publicfs01/fs1-m8/home/m8s001451/zayf/VIC_Data/vic_forcing/vic_forcing_"

# 你当前 ERA5 下载结果的根目录（按年分子目录）
ERA5_DAILY_ROOT = "/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/cmf_v420_pkg/VIC_Data/ERA5"
ERA5_DAILY_FILE_TEMPLATE = "ERA5_VIC_global_{year}_{month}_{day}.nc"
ERA5_EXPECTED_HOURS = [0, 6, 12, 18]
ERA5_REQUIRE_COMPLETE_YEAR = True   # 缺天则报错；设为 False 则尽量拼接已有文件

# 如果你仍然想兼容“每年一个综合 nc”的旧方式，可以继续保留这个模板
ERA5_INPUT_TEMPLATE = "/publicfs01/fs1-m8/home/m8s001451/zayf/VIC_Data/ERA5/era5_{year}.nc"

# 你的 ERA5 文件中的变量名映射（已按 downloadERA5.py 的长变量名适配）
ERA5_VARMAP = {
    "t2m": "2m_temperature",                     # [K]
    "tp": "total_precipitation",                # [m per saved step] or [m/day]
    "sp": "surface_pressure",                   # [Pa]
    "ssrd": "surface_solar_radiation_downwards",# [J/m2 over saved step] or [W/m2]
    "strd": "surface_thermal_radiation_downwards",# [J/m2 over saved step] or [W/m2]
    "d2m": "2m_dewpoint_temperature",           # [K]
    "u10": "10m_u_component_of_wind",           # [m/s]
    "v10": "10m_v_component_of_wind",           # [m/s]
}

# 坐标名自动识别失败时，可手动指定
ERA5_TIME_NAME = None
ERA5_LAT_NAME = None
ERA5_LON_NAME = None

# ERA5 单位/格式控制
ERA5_PREC_IN_M = True                 # tp 是米 -> 乘 1000 变 mm
ERA5_RAD_IS_ACCUM_J = True            # ssrd/strd 是每步累计 J/m2 -> 除以每步秒数变 W/m2
ERA5_FORCE_TO_DAILY = False           # True: 把 subdaily 强制聚合为 daily；False: 保持原时间步

# VIC 时间步设置
# 你当前下载脚本每天保留 4 个时次：00/06/12/18 UTC，因此默认设为 4。
# 如果后面你决定把 ERA5 先聚合成 daily forcing，再改回 1。
VIC_MODEL_STEPS_PER_DAY = 4
VIC_SNOW_STEPS_PER_DAY = 4
VIC_RUNOFF_STEPS_PER_DAY = 4

# ---------- VIC 参数数据 ----------
DIR_OPENLANDMAP = "/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/cmf_v420_pkg/VIC_Data"
LAND_DEM_NC = "/publicfs01/fs1-m8/home/m8s001451/zayf/HydroMT-SFINCS/examples/data/gebco/GEBCO_2025_sub_ice.nc"
GLOBAL_SOIL_TXT = "/publicfs01/fs1-m8/home/m8s001451/zayf/VIC_Data/global/global_soil_param_new.txt"
ANNUAL_PRECIP_TIF = "/publicfs01/fs1-m8/home/m8s001451/zayf/VIC_Data/ERA5/annual_mean_precip.tif"

# 如果你有更好的陆地 DEM（如 FABDEM / MERIT DEM），把 LAND_DEM_NC 改掉更合理

# 工作目录
WORK_ROOT = os.path.join(CAMA_DIR, "gosh", f"work_{REGION_NAME}_image")
VIC_WORKDIR = os.path.join(WORK_ROOT, "vic_image")
VIC_FORCING_DIR = os.path.join(VIC_WORKDIR, "forcing")
VIC_PARAM_DIR = os.path.join(VIC_WORKDIR, "params")
VIC_RESULT_DIR = os.path.join(VIC_WORKDIR, "result")
VIC_LOG_DIR = os.path.join(VIC_WORKDIR, "log")
CAMA_RUNOFF_NC_DIR = os.path.join(WORK_ROOT, "cama_runoff_from_vic")

# USDA soil texture -> VIC 参数查找表
SOIL_LUT = {
    1:  [0.36, 0.17, 127.057, 27.56, 0.098],
    2:  [0.37, 0.25, 96.161,  22.52, 0.103],
    3:  [0.31, 0.23, 208.449, 29.00, 0.098],
    4:  [0.34, 0.21, 113.240, 19.04, 0.081],
    5:  [0.36, 0.21, 116.950, 17.96, 0.086],
    6:  [0.27, 0.17, 146.218, 20.32, 0.060],
    7:  [0.29, 0.14, 123.880, 13.60, 0.088],
    8:  [0.32, 0.12, 64.863, 10.58, 0.140],
    9:  [0.21, 0.09, 291.743, 12.68, 0.057],
    10: [0.28, 0.08, 306.902, 9.10, 0.071],
    11: [0.15, 0.06, 895.365, 10.98, 0.050],
    12: [0.08, 0.03, 5794.287, 11.20, 0.050]
}

# 简化植被参数（单一草地 tile 占位）
DEFAULT_VEG_CLASS = 1
DEFAULT_ROOT_DEPTH = np.array([0.10, 0.50, 1.00], dtype=np.float64)
DEFAULT_ROOT_FRACT = np.array([0.10, 0.30, 0.60], dtype=np.float64)
DEFAULT_LAI = np.array([1.0] * 12, dtype=np.float64)
DEFAULT_ALBEDO = np.array([0.20] * 12, dtype=np.float64)
DEFAULT_VEG_ROUGH = np.array([0.037] * 12, dtype=np.float64)      # ~0.123 * 0.3 m
DEFAULT_DISPLACEMENT = np.array([0.201] * 12, dtype=np.float64)   # ~0.67 * 0.3 m

# =========================================================
# 2. 通用工具函数
# =========================================================
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def ensure_file_exists(path: str, desc: str = "file") -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"❌ Missing {desc}: {path}")


def ensure_nonempty_file(path: str, desc: str = "file") -> None:
    ensure_file_exists(path, desc)
    if os.path.getsize(path) == 0:
        raise RuntimeError(f"❌ Empty {desc}: {path}")


def run_cmd(cmd: str, cwd: Optional[str] = None) -> None:
    print(f"   >>> Executing: {cmd}")
    result = subprocess.run(cmd, shell=True, cwd=cwd, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"❌ Command failed: {cmd}\n   cwd={cwd}")


def robust_modify_shell_vars(filepath: str, replacements: Dict[str, str]) -> None:
    ensure_file_exists(filepath, "script")
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    for key, val in replacements.items():
        val = str(val)
        pattern = rf"(?m)^(\s*(?:export\s+)?{re.escape(key)}\s*=\s*).*$"
        if re.search(pattern, content):
            content = re.sub(pattern, lambda m, v=val: f"{m.group(1)}{v}", content)
            print(f"      ✅ Updated shell var: {key}={val}")
        else:
            content += f"\n{key}={val}\n"
            print(f"      ➕ Added shell var: {key}={val}")

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)
    subprocess.run(["chmod", "+x", filepath], check=False)


def get_bbox_from_geojson(geojson_file: str, buffer: float = 1.0) -> Tuple[int, int, int, int]:
    print(f"\n[STEP] Analyzing spatial extent from: {geojson_file}")
    ensure_file_exists(geojson_file, "GeoJSON")
    gdf = gpd.read_file(geojson_file)
    minx, miny, maxx, maxy = gdf.total_bounds
    west = int(np.floor(minx - buffer))
    east = int(np.ceil(maxx + buffer))
    south = int(np.floor(miny - buffer))
    north = int(np.ceil(maxy + buffer))
    print(f"   ✅ BBOX = [{west}, {east}, {south}, {north}]")
    return west, east, south, north


def get_resolution_deg(resolution: str) -> float:
    res_map = {
        "15min": 15.0 / 60.0,
        "06min": 6.0 / 60.0,
        "05min": 5.0 / 60.0,
        "03min": 3.0 / 60.0,
        "01min": 1.0 / 60.0,
    }
    if resolution not in res_map:
        raise ValueError(f"❌ Unsupported resolution: {resolution}")
    return res_map[resolution]


def infer_coord_name(ds: xr.Dataset, candidates: List[str]):
    lower_map = {name.lower(): name for name in list(ds.coords) + list(ds.dims)}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def sample_raster(tif_path: str, lons, lats, multiplier: float = 1.0, nodata_val=None) -> np.ndarray:
    if not os.path.exists(tif_path):
        print(f"      ⚠️ Warning: Missing raster {os.path.basename(tif_path)}. Filling NaN.")
        return np.full(len(lons), np.nan, dtype=float)

    with rasterio.open(tif_path) as src:
        pts = [(float(lon), float(lat)) for lon, lat in zip(lons, lats)]
        vals = np.array([v[0] for v in src.sample(pts)], dtype=float)
        src_nodata = src.nodata

    if src_nodata is not None:
        vals[vals == src_nodata] = np.nan
    if nodata_val is not None:
        vals[vals == nodata_val] = np.nan
    return vals * multiplier


def sample_netcdf_dem(nc_path: str, lons, lats, var_name: Optional[str] = None, method: str = "nearest") -> np.ndarray:
    if not os.path.exists(nc_path):
        print("      ⚠️ Warning: DEM NetCDF missing. Filling NaN.")
        return np.full(len(lons), np.nan, dtype=float)

    try:
        with xr.open_dataset(nc_path) as ds:
            if var_name is None:
                var_name = "elevation" if "elevation" in ds.data_vars else list(ds.data_vars)[0]
            lon_name = infer_coord_name(ds, ["lon", "longitude", "x"])
            lat_name = infer_coord_name(ds, ["lat", "latitude", "y"])
            if lon_name is None or lat_name is None:
                raise ValueError("Cannot identify lon/lat coordinate names in DEM NetCDF.")

            lon_da = xr.DataArray(np.asarray(lons, dtype=float), dims="points")
            lat_da = xr.DataArray(np.asarray(lats, dtype=float), dims="points")
            vals = ds[var_name].sel({lon_name: lon_da, lat_name: lat_da}, method=method).values
            return np.asarray(vals, dtype=float)
    except Exception as e:
        print(f"      ⚠️ Warning: DEM extraction failed: {e}")
        return np.full(len(lons), np.nan, dtype=float)


def fill_nan_by_nearest_xy(values: np.ndarray, lons: np.ndarray, lats: np.ndarray, default_value: float) -> np.ndarray:
    vals = np.asarray(values, dtype=float).copy()
    valid = np.isfinite(vals)
    if valid.all():
        return vals
    if not valid.any():
        vals[:] = default_value
        return vals
    tree = cKDTree(np.c_[lons[valid], lats[valid]])
    _, idx = tree.query(np.c_[lons[~valid], lats[~valid]], k=1)
    vals[~valid] = vals[valid][idx]
    return vals


def compute_regular_latlon_cell_area(lat_center_deg: np.ndarray, dlat_deg: float, dlon_deg: float) -> np.ndarray:
    """规则经纬网格单元面积 [m2]，返回 shape=(nlat,)"""
    R = 6371000.0
    lat1 = np.deg2rad(lat_center_deg - dlat_deg / 2.0)
    lat2 = np.deg2rad(lat_center_deg + dlat_deg / 2.0)
    dlon = np.deg2rad(dlon_deg)
    area = (R ** 2) * dlon * (np.sin(lat2) - np.sin(lat1))
    return np.abs(area)


def open_year_dataset(path: str) -> xr.Dataset:
    ensure_file_exists(path, "NetCDF forcing file")
    ds = xr.open_dataset(path)

    time_name = ERA5_TIME_NAME or infer_coord_name(ds, ["time"])
    lat_name = ERA5_LAT_NAME or infer_coord_name(ds, ["lat", "latitude", "y"])
    lon_name = ERA5_LON_NAME or infer_coord_name(ds, ["lon", "longitude", "x"])
    if time_name is None or lat_name is None or lon_name is None:
        raise ValueError(f"❌ Cannot infer time/lat/lon coordinates from {path}")

    rename_map = {}
    if time_name != "time":
        rename_map[time_name] = "time"
    if lat_name != "lat":
        rename_map[lat_name] = "lat"
    if lon_name != "lon":
        rename_map[lon_name] = "lon"
    if rename_map:
        ds = ds.rename(rename_map)

    if np.any(np.diff(ds["lat"].values) < 0):
        ds = ds.sortby("lat")
    if np.any(np.diff(ds["lon"].values) < 0):
        ds = ds.sortby("lon")
    return ds


def get_era5_daily_file(year: int, month: int, day: int) -> str:
    year_dir = os.path.join(ERA5_DAILY_ROOT, f"{year:04d}")
    return os.path.join(
        year_dir,
        ERA5_DAILY_FILE_TEMPLATE.format(year=year, month=f"{month:02d}", day=f"{day:02d}")
    )


def open_era5_daily_dataset_for_year(year: int) -> xr.Dataset:
    """
    读取你当前按天下载的 ERA5 文件，并拼接成该年的一个 Dataset。
    文件结构应为：
      ERA5_DAILY_ROOT/YYYY/ERA5_VIC_global_YYYY_MM_DD.nc
    """
    dates = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
    files = []
    missing = []
    for d in dates:
        fp = get_era5_daily_file(d.year, d.month, d.day)
        if os.path.exists(fp):
            files.append(fp)
        else:
            missing.append(fp)

    if ERA5_REQUIRE_COMPLETE_YEAR and missing:
        raise FileNotFoundError(
            f"❌ ERA5 daily files missing for year {year}. Missing count = {len(missing)}.\n"
            f"First missing file: {missing[0]}"
        )

    if not files:
        raise FileNotFoundError(f"❌ No ERA5 daily NetCDF files found for year {year} under {ERA5_DAILY_ROOT}")

    print(f"   >>> Found {len(files)} daily ERA5 files for {year}")
    ds = xr.open_mfdataset(files, combine="by_coords")
    ds = open_year_dataset_from_ds(ds)

    if ERA5_EXPECTED_HOURS is not None and len(ERA5_EXPECTED_HOURS) > 0:
        hours = pd.DatetimeIndex(ds["time"].values).hour
        unique_hours = sorted(set(int(h) for h in hours))
        if unique_hours != sorted(ERA5_EXPECTED_HOURS):
            print(f"      ⚠️ Warning: expected hours {ERA5_EXPECTED_HOURS}, but found {unique_hours}")

    return ds


def open_year_dataset_from_ds(ds: xr.Dataset) -> xr.Dataset:
    time_name = ERA5_TIME_NAME or infer_coord_name(ds, ["time"])
    lat_name = ERA5_LAT_NAME or infer_coord_name(ds, ["lat", "latitude", "y"])
    lon_name = ERA5_LON_NAME or infer_coord_name(ds, ["lon", "longitude", "x"])
    if time_name is None or lat_name is None or lon_name is None:
        raise ValueError("❌ Cannot infer time/lat/lon coordinates from in-memory ERA5 dataset")

    rename_map = {}
    if time_name != "time":
        rename_map[time_name] = "time"
    if lat_name != "lat":
        rename_map[lat_name] = "lat"
    if lon_name != "lon":
        rename_map[lon_name] = "lon"
    if rename_map:
        ds = ds.rename(rename_map)

    if np.any(np.diff(ds["lat"].values) < 0):
        ds = ds.sortby("lat")
    if np.any(np.diff(ds["lon"].values) < 0):
        ds = ds.sortby("lon")
    return ds


def infer_steps_per_day(time_index: pd.DatetimeIndex) -> int:
    if len(time_index) < 2:
        return 1
    dt_sec = int((time_index[1] - time_index[0]).total_seconds())
    if dt_sec <= 0:
        return 1
    steps = int(round(86400 / dt_sec))
    return max(1, steps)


def dewpoint_to_vapor_pressure_kpa(td_c: xr.DataArray) -> xr.DataArray:
    """Tetens 公式，输出 kPa"""
    return 0.6108 * np.exp((17.27 * td_c) / (td_c + 237.3))

# =========================================================
# 3. CaMa-Flood 区域域构建
# =========================================================
def build_cama_regional_map() -> Tuple[str, Tuple[int, int, int, int], float]:
    print("\n=======================================================")
    print("[PHASE 1] Building CaMa-Flood regional domain")
    print("=======================================================")

    reg_map_dir = os.path.join(CAMA_DIR, "map", f"reg_{REGION_NAME}_{CAMA_RESOLUTION}")
    src_region_dir = os.path.join(reg_map_dir, "src_region")
    src_param_dir = os.path.join(reg_map_dir, "src_param")

    bbox = get_bbox_from_geojson(GEOJSON_PATH, buffer=BUFFER_DEG)
    cama_res_deg = get_resolution_deg(CAMA_RESOLUTION)

    print(f"\n[STEP 1/4] Creating regional map for {REGION_NAME}")
    ensure_dir(reg_map_dir)
    if not os.path.exists(src_region_dir):
        shutil.copytree(os.path.join(CAMA_DIR, "map", "src", "src_region"), src_region_dir)

    run_cmd("make all", cwd=src_region_dir)
    script_s01_reg = os.path.join(src_region_dir, "s01-regional_map.sh")
    path_str = f'"../../glb_{CAMA_RESOLUTION}/"'
    robust_modify_shell_vars(script_s01_reg, {
        "SOURCE": path_str,
        "MAPGLB": path_str,
        "CDIR": path_str,
        "WEST": str(bbox[0]),
        "EAST": str(bbox[1]),
        "SOUTH": str(bbox[2]),
        "NORTH": str(bbox[3]),
    })
    run_cmd("./s01-regional_map.sh", cwd=src_region_dir)
    ensure_nonempty_file(os.path.join(reg_map_dir, "params.txt"), "params.txt")
    ensure_nonempty_file(os.path.join(reg_map_dir, "nextxy.bin"), "nextxy.bin")
    print("   ✅ Regional map created and validated")

    print("\n[STEP 2/4] Generating channel parameters")
    if not os.path.exists(src_param_dir):
        shutil.copytree(os.path.join(CAMA_DIR, "map", "src", "src_param"), src_param_dir)
    run_cmd("make all", cwd=src_param_dir)
    script_s01_param = os.path.join(src_param_dir, "s01-channel_params.sh")
    subprocess.run(["chmod", "+x", script_s01_param], check=False)
    run_cmd("./s01-channel_params.sh", cwd=src_param_dir)
    ensure_nonempty_file(os.path.join(reg_map_dir, "rivwth.bin"), "rivwth.bin")
    print("   ✅ Channel parameters generated")

    print("\n[STEP 3/4] Generating runoff input matrix for default sample forcing")
    script_s02_inp = os.path.join(src_param_dir, "s02-generate_inpmat.sh")
    subprocess.run(["chmod", "+x", script_s02_inp], check=False)
    run_cmd("./s02-generate_inpmat.sh", cwd=src_param_dir)

    print("\n[STEP 4/4] Regional map ready")
    return reg_map_dir, bbox, cama_res_deg

# =========================================================
# 4. 6min VIC 网格与 mask / frac
# =========================================================
def build_vic_support_grid(
    reg_map_dir: str,
    cama_bbox: Tuple[int, int, int, int],
    cama_res_deg: float,
    vic_bbox: Tuple[int, int, int, int],
    vic_res_deg: float,
) -> Dict[str, np.ndarray]:
    print("\n=======================================================")
    print("[PHASE 2] Building 6min VIC support grid from CaMa domain")
    print("=======================================================")

    west_c, east_c, south_c, north_c = cama_bbox
    nx_c = int(round((east_c - west_c) / cama_res_deg))
    ny_c = int(round((north_c - south_c) / cama_res_deg))

    nextxy_path = os.path.join(reg_map_dir, "nextxy.bin")
    ensure_nonempty_file(nextxy_path, "nextxy.bin")
    raw = np.fromfile(nextxy_path, dtype="<i4")
    expected = 2 * ny_c * nx_c
    if raw.size != expected:
        raise RuntimeError(
            f"❌ nextxy.bin size mismatch. expected={expected}, got={raw.size}."
        )
    nextxy = raw.reshape(2, ny_c, nx_c)
    cama_active = (nextxy[0] != -9999)

    west_v, east_v, south_v, north_v = vic_bbox
    nx_v = int(round((east_v - west_v) / vic_res_deg))
    ny_v = int(round((north_v - south_v) / vic_res_deg))

    lon_v = west_v + (np.arange(nx_v) + 0.5) * vic_res_deg
    lat_v = north_v - (np.arange(ny_v) + 0.5) * vic_res_deg

    lon_c = west_c + (np.arange(nx_c) + 0.5) * cama_res_deg
    lat_c = north_c - (np.arange(ny_c) + 0.5) * cama_res_deg
    lon2_c, lat2_c = np.meshgrid(lon_c, lat_c)

    counts = np.zeros((ny_v, nx_v), dtype=np.int32)
    active_counts = np.zeros((ny_v, nx_v), dtype=np.int32)

    ix = np.floor((lon2_c - west_v) / vic_res_deg).astype(int)
    iy = np.floor((north_v - lat2_c) / vic_res_deg).astype(int)
    inside = (ix >= 0) & (ix < nx_v) & (iy >= 0) & (iy < ny_v)

    np.add.at(counts, (iy[inside], ix[inside]), 1)
    np.add.at(active_counts, (iy[inside & cama_active], ix[inside & cama_active]), 1)

    frac = np.zeros((ny_v, nx_v), dtype=np.float64)
    valid = counts > 0
    frac[valid] = active_counts[valid] / counts[valid]
    mask = (frac > 0).astype(np.int32)

    area_lat = compute_regular_latlon_cell_area(lat_v, vic_res_deg, vic_res_deg)
    area = np.repeat(area_lat[:, None], nx_v, axis=1)

    print(f"   ✅ VIC 6min grid: ny={ny_v}, nx={nx_v}")
    print(f"   ✅ Active VIC cells: {int(mask.sum())}")

    return {
        "lat": lat_v,
        "lon": lon_v,
        "mask": mask,
        "frac": frac,
        "area": area,
    }

# =========================================================
# 5. VIC Image Driver 参数文件（domain.nc / params.nc / global.txt）
# =========================================================
def build_vic_image_domain_nc(vic_grid: Dict[str, np.ndarray], out_path: str) -> str:
    print("\n[STEP] Writing VIC Image domain.nc")
    ensure_dir(os.path.dirname(out_path))

    ds = xr.Dataset(
        data_vars={
            "mask": (("lat", "lon"), vic_grid["mask"].astype(np.int32)),
            "area": (("lat", "lon"), vic_grid["area"].astype(np.float64)),
            "frac": (("lat", "lon"), vic_grid["frac"].astype(np.float64)),
        },
        coords={
            "lat": vic_grid["lat"].astype(np.float64),
            "lon": vic_grid["lon"].astype(np.float64),
        },
        attrs={"description": "VIC Image Driver domain file"}
    )
    ds["lat"].attrs.update({"units": "degrees_north", "long_name": "latitude of grid cell center"})
    ds["lon"].attrs.update({"units": "degrees_east", "long_name": "longitude of grid cell center"})
    ds["mask"].attrs.update({"long_name": "mask"})
    ds["area"].attrs.update({"units": "m2", "long_name": "area"})
    ds["frac"].attrs.update({"long_name": "frac"})

    ds.to_netcdf(out_path)
    print(f"   ✅ {out_path}")
    return out_path


def generate_vic_soil_dataframe(vic_grid: Dict[str, np.ndarray]) -> pd.DataFrame:
    print("\n[STEP] Sampling soil/topography for VIC 6min parameter grid")
    lat_vals = vic_grid["lat"]
    lon_vals = vic_grid["lon"]
    mask = vic_grid["mask"]

    lon2, lat2 = np.meshgrid(lon_vals, lat_vals)
    run_idx = mask.astype(bool)

    df = pd.DataFrame({
        "lat": lat2[run_idx],
        "lon": lon2[run_idx],
        "gridcell": np.arange(1, int(run_idx.sum()) + 1, dtype=np.int32),
    })

    # USDA texture classes
    for i, depth in enumerate([0, 30, 100], start=1):
        tif_name = f"sol_texture.class_usda.tt_m_250m_b{depth}..{depth}cm_1950..2017_v0.2.tif"
        tif = os.path.join(DIR_OPENLANDMAP, tif_name)
        classes = sample_raster(tif, df["lon"].values, df["lat"].values, nodata_val=255)
        classes = fill_nan_by_nearest_xy(classes, df["lon"].values, df["lat"].values, default_value=1.0)
        df[f"class_{i}"] = np.clip(np.rint(classes), 1, 12).astype(int)

    # bulk density
    for i, depth in enumerate([0, 30, 100], start=1):
        tif_name = f"sol_bulkdens.fineearth_usda.4a1h_m_250m_b{depth}..{depth}cm_1950..2017_v0.2.tif"
        tif = os.path.join(DIR_OPENLANDMAP, tif_name)
        bd = sample_raster(tif, df["lon"].values, df["lat"].values, multiplier=10.0, nodata_val=255)
        bd[bd < 800] = np.nan
        bd = fill_nan_by_nearest_xy(bd, df["lon"].values, df["lat"].values, default_value=1300.0)
        df[f"bulk_density_{i}"] = bd

    # DEM / annual precip
    elev = sample_netcdf_dem(LAND_DEM_NC, df["lon"].values, df["lat"].values, method="nearest")
    elev = fill_nan_by_nearest_xy(elev, df["lon"].values, df["lat"].values, default_value=10.0)
    df["elev"] = elev

    annual_prec = sample_raster(ANNUAL_PRECIP_TIF, df["lon"].values, df["lat"].values)
    annual_prec = fill_nan_by_nearest_xy(annual_prec, df["lon"].values, df["lat"].values, default_value=1000.0)
    df["annual_prec"] = annual_prec

    # Global soil txt (IDW)
    soil_cols = [
        "init_moist_1", "init_moist_2", "init_moist_3", "avg_T",
        "bubble_1", "bubble_2", "bubble_3",
        "quartz_1", "quartz_2", "quartz_3"
    ]
    if os.path.exists(GLOBAL_SOIL_TXT):
        global_df = pd.read_csv(GLOBAL_SOIL_TXT, sep=r"\s+", header=None, low_memory=False)
        gl_lons = global_df[3].values
        gl_lats = global_df[2].values
        tree = cKDTree(np.c_[gl_lons, gl_lats])
        distances, indices = tree.query(np.c_[df["lon"].values, df["lat"].values], k=4)
        distances = np.maximum(distances, 1e-6)
        weights = 1.0 / (distances ** 2)
        weights /= weights.sum(axis=1)[:, np.newaxis]

        def apply_idw(col_idx: int) -> np.ndarray:
            vals = global_df[col_idx].values[indices]
            return np.sum(vals * weights, axis=1)

        df["init_moist_1"] = apply_idw(18)
        df["init_moist_2"] = apply_idw(19)
        df["init_moist_3"] = apply_idw(20)
        df["avg_T"] = apply_idw(25)
        df["bubble_1"] = apply_idw(27)
        df["bubble_2"] = apply_idw(28)
        df["bubble_3"] = apply_idw(29)
        df["quartz_1"] = apply_idw(30)
        df["quartz_2"] = apply_idw(31)
        df["quartz_3"] = apply_idw(32)
        for col in soil_cols:
            df[col] = fill_nan_by_nearest_xy(df[col].values, df["lon"].values, df["lat"].values, default_value=10.0)
    else:
        print("      ⚠️ Warning: GLOBAL_SOIL_TXT not found. Using placeholders.")
        for col in soil_cols:
            df[col] = 10.0

    # Map lookup table
    for i in [1, 2, 3]:
        cls = df[f"class_{i}"]
        df[f"Wcr_{i}"] = cls.map(lambda c: SOIL_LUT.get(int(c), SOIL_LUT[1])[0])
        df[f"Wpwp_{i}"] = cls.map(lambda c: SOIL_LUT.get(int(c), SOIL_LUT[1])[1])
        df[f"Ksat_{i}"] = cls.map(lambda c: SOIL_LUT.get(int(c), SOIL_LUT[1])[2])
        df[f"expt_{i}"] = cls.map(lambda c: SOIL_LUT.get(int(c), SOIL_LUT[1])[3])
        df[f"resid_moist_{i}"] = cls.map(lambda c: SOIL_LUT.get(int(c), SOIL_LUT[1])[4])
        df[f"soil_density_{i}"] = 2685.0
        porosity = 1.0 - df[f"bulk_density_{i}"] / df[f"soil_density_{i}"]
        max_resid = df[f"Wpwp_{i}"] * porosity * 0.95
        df[f"resid_moist_{i}"] = np.minimum(df[f"resid_moist_{i}"], max_resid)
        df[f"phi_s_{i}"] = -999.0

    # fixed params
    df["infilt"] = 0.5
    df["Ds"] = 0.01
    df["Dsmax"] = 8.0
    df["Ws"] = 0.8
    df["c"] = 2.0
    df["depth_1"] = 0.1
    df["depth_2"] = 0.5
    df["depth_3"] = 2.0
    df["dp"] = 4.0
    df["off_gmt"] = df["lon"] * 24.0 / 360.0
    df["rough"] = 0.01
    df["snow_rough"] = 0.03
    df["fs_active"] = 1
    df["July_Tavg"] = 1.0

    print(f"   ✅ Sampled {len(df)} active VIC cells")
    return df


def build_vic_image_params_nc(vic_grid: Dict[str, np.ndarray], out_path: str) -> str:
    print("\n[STEP] Writing VIC Image params.nc")
    ensure_dir(os.path.dirname(out_path))

    df = generate_vic_soil_dataframe(vic_grid)
    lat_vals = vic_grid["lat"]
    lon_vals = vic_grid["lon"]
    ny, nx = len(lat_vals), len(lon_vals)
    mask = vic_grid["mask"].astype(np.int32)

    # Row/col index for active cells
    lat_to_idx = {float(v): i for i, v in enumerate(lat_vals)}
    lon_to_idx = {float(v): j for j, v in enumerate(lon_vals)}

    # helper arrays
    fill = np.nan
    int_fill = -2147483647
    run_cell = np.zeros((ny, nx), dtype=np.int32)
    gridcell = np.zeros((ny, nx), dtype=np.int32)
    lats2 = np.repeat(lat_vals[:, None], nx, axis=1)
    lons2 = np.repeat(lon_vals[None, :], ny, axis=0)

    nlayer = 3
    nveg = 1
    nroot = 3
    nmonth = 12

    # allocate
    infilt = np.full((ny, nx), fill)
    Ds = np.full((ny, nx), fill)
    Dsmax = np.full((ny, nx), fill)
    Ws = np.full((ny, nx), fill)
    c = np.full((ny, nx), fill)
    expt = np.full((nlayer, ny, nx), fill)
    Ksat = np.full((nlayer, ny, nx), fill)
    phi_s = np.full((nlayer, ny, nx), fill)
    init_moist = np.full((nlayer, ny, nx), fill)
    elev = np.full((ny, nx), fill)
    depth = np.full((nlayer, ny, nx), fill)
    avg_T = np.full((ny, nx), fill)
    dp = np.full((ny, nx), fill)
    bubble = np.full((nlayer, ny, nx), fill)
    quartz = np.full((nlayer, ny, nx), fill)
    bulk_density = np.full((nlayer, ny, nx), fill)
    soil_density = np.full((nlayer, ny, nx), fill)
    off_gmt = np.full((ny, nx), fill)
    Wcr_FRACT = np.full((nlayer, ny, nx), fill)
    Wpwp_FRACT = np.full((nlayer, ny, nx), fill)
    rough = np.full((ny, nx), fill)
    snow_rough = np.full((ny, nx), fill)
    annual_prec = np.full((ny, nx), fill)
    resid_moist = np.full((nlayer, ny, nx), fill)
    fs_active = np.zeros((ny, nx), dtype=np.int32)
    July_Tavg = np.full((ny, nx), fill)

    Nveg_arr = np.zeros((ny, nx), dtype=np.int32)
    Cv = np.zeros((nveg, ny, nx), dtype=np.float64)
    root_depth = np.zeros((nveg, nroot, ny, nx), dtype=np.float64)
    root_fract = np.zeros((nveg, nroot, ny, nx), dtype=np.float64)
    LAI = np.zeros((nveg, nmonth, ny, nx), dtype=np.float64)
    overstory = np.zeros((nveg, ny, nx), dtype=np.int32)
    rarc = np.zeros((nveg, ny, nx), dtype=np.float64)
    rmin = np.zeros((nveg, ny, nx), dtype=np.float64)
    wind_h = np.zeros((nveg, ny, nx), dtype=np.float64)
    RGL = np.zeros((nveg, ny, nx), dtype=np.float64)
    rad_atten = np.zeros((nveg, ny, nx), dtype=np.float64)
    wind_atten = np.zeros((nveg, ny, nx), dtype=np.float64)
    trunk_ratio = np.zeros((nveg, ny, nx), dtype=np.float64)
    albedo = np.zeros((nveg, nmonth, ny, nx), dtype=np.float64)
    veg_rough = np.zeros((nveg, nmonth, ny, nx), dtype=np.float64)
    displacement = np.zeros((nveg, nmonth, ny, nx), dtype=np.float64)

    # fill active cells
    for _, row in df.iterrows():
        i = lat_to_idx[float(row["lat"])]
        j = lon_to_idx[float(row["lon"])]
        run_cell[i, j] = 1
        gridcell[i, j] = int(row["gridcell"])
        infilt[i, j] = row["infilt"]
        Ds[i, j] = row["Ds"]
        Dsmax[i, j] = row["Dsmax"]
        Ws[i, j] = row["Ws"]
        c[i, j] = row["c"]
        elev[i, j] = row["elev"]
        avg_T[i, j] = row["avg_T"]
        dp[i, j] = row["dp"]
        off_gmt[i, j] = row["off_gmt"]
        rough[i, j] = row["rough"]
        snow_rough[i, j] = row["snow_rough"]
        annual_prec[i, j] = row["annual_prec"]
        fs_active[i, j] = int(row["fs_active"])
        July_Tavg[i, j] = row["July_Tavg"]

        expt[:, i, j] = [row["expt_1"], row["expt_2"], row["expt_3"]]
        Ksat[:, i, j] = [row["Ksat_1"], row["Ksat_2"], row["Ksat_3"]]
        phi_s[:, i, j] = [row["phi_s_1"], row["phi_s_2"], row["phi_s_3"]]
        init_moist[:, i, j] = [row["init_moist_1"], row["init_moist_2"], row["init_moist_3"]]
        depth[:, i, j] = [row["depth_1"], row["depth_2"], row["depth_3"]]
        bubble[:, i, j] = [row["bubble_1"], row["bubble_2"], row["bubble_3"]]
        quartz[:, i, j] = [row["quartz_1"], row["quartz_2"], row["quartz_3"]]
        bulk_density[:, i, j] = [row["bulk_density_1"], row["bulk_density_2"], row["bulk_density_3"]]
        soil_density[:, i, j] = [row["soil_density_1"], row["soil_density_2"], row["soil_density_3"]]
        Wcr_FRACT[:, i, j] = [row["Wcr_1"], row["Wcr_2"], row["Wcr_3"]]
        Wpwp_FRACT[:, i, j] = [row["Wpwp_1"], row["Wpwp_2"], row["Wpwp_3"]]
        resid_moist[:, i, j] = [row["resid_moist_1"], row["resid_moist_2"], row["resid_moist_3"]]

        # single grass tile
        Nveg_arr[i, j] = 1
        Cv[0, i, j] = 1.0
        root_depth[0, :, i, j] = DEFAULT_ROOT_DEPTH
        root_fract[0, :, i, j] = DEFAULT_ROOT_FRACT
        LAI[0, :, i, j] = DEFAULT_LAI
        overstory[0, i, j] = 0
        rarc[0, i, j] = 2.0
        rmin[0, i, j] = 100.0
        wind_h[0, i, j] = 10.0
        RGL[0, i, j] = 100.0
        rad_atten[0, i, j] = 0.5
        wind_atten[0, i, j] = 0.5
        trunk_ratio[0, i, j] = 0.2
        albedo[0, :, i, j] = DEFAULT_ALBEDO
        veg_rough[0, :, i, j] = DEFAULT_VEG_ROUGH
        displacement[0, :, i, j] = DEFAULT_DISPLACEMENT

    # outside active cells but inside domain arrays
    run_cell = np.where(mask > 0, run_cell, 0)
    Nveg_arr = np.where(mask > 0, Nveg_arr, 0)

    ds = xr.Dataset(
        coords={
            "lat": lat_vals.astype(np.float64),
            "lon": lon_vals.astype(np.float64),
            "nlayer": np.arange(1, nlayer + 1, dtype=np.int32),
            "veg_class": np.array([DEFAULT_VEG_CLASS], dtype=np.int32),
            "root_zone": np.arange(1, nroot + 1, dtype=np.int32),
            "month": np.arange(1, nmonth + 1, dtype=np.int32),
        },
        data_vars={
            "run_cell": (("lat", "lon"), run_cell),
            "gridcell": (("lat", "lon"), gridcell),
            "gridcel": (("lat", "lon"), gridcell),
            "lats": (("lat", "lon"), lats2),
            "lons": (("lat", "lon"), lons2),
            "infilt": (("lat", "lon"), infilt),
            "Ds": (("lat", "lon"), Ds),
            "Dsmax": (("lat", "lon"), Dsmax),
            "Ws": (("lat", "lon"), Ws),
            "c": (("lat", "lon"), c),
            "expt": (("nlayer", "lat", "lon"), expt),
            "Ksat": (("nlayer", "lat", "lon"), Ksat),
            "phi_s": (("nlayer", "lat", "lon"), phi_s),
            "init_moist": (("nlayer", "lat", "lon"), init_moist),
            "elev": (("lat", "lon"), elev),
            "depth": (("nlayer", "lat", "lon"), depth),
            "avg_T": (("lat", "lon"), avg_T),
            "dp": (("lat", "lon"), dp),
            "bubble": (("nlayer", "lat", "lon"), bubble),
            "quartz": (("nlayer", "lat", "lon"), quartz),
            "bulk_density": (("nlayer", "lat", "lon"), bulk_density),
            "soil_density": (("nlayer", "lat", "lon"), soil_density),
            "off_gmt": (("lat", "lon"), off_gmt),
            "Wcr_FRACT": (("nlayer", "lat", "lon"), Wcr_FRACT),
            "Wpwp_FRACT": (("nlayer", "lat", "lon"), Wpwp_FRACT),
            "rough": (("lat", "lon"), rough),
            "snow_rough": (("lat", "lon"), snow_rough),
            "annual_prec": (("lat", "lon"), annual_prec),
            "resid_moist": (("nlayer", "lat", "lon"), resid_moist),
            "fs_active": (("lat", "lon"), fs_active),
            "July_Tavg": (("lat", "lon"), July_Tavg),
            "Nveg": (("lat", "lon"), Nveg_arr),
            "Cv": (("veg_class", "lat", "lon"), Cv),
            "root_depth": (("veg_class", "root_zone", "lat", "lon"), root_depth),
            "root_fract": (("veg_class", "root_zone", "lat", "lon"), root_fract),
            "LAI": (("veg_class", "month", "lat", "lon"), LAI),
            "overstory": (("veg_class", "lat", "lon"), overstory),
            "rarc": (("veg_class", "lat", "lon"), rarc),
            "rmin": (("veg_class", "lat", "lon"), rmin),
            "wind_h": (("veg_class", "lat", "lon"), wind_h),
            "RGL": (("veg_class", "lat", "lon"), RGL),
            "rad_atten": (("veg_class", "lat", "lon"), rad_atten),
            "wind_atten": (("veg_class", "lat", "lon"), wind_atten),
            "trunk_ratio": (("veg_class", "lat", "lon"), trunk_ratio),
            "albedo": (("veg_class", "month", "lat", "lon"), albedo),
            "veg_rough": (("veg_class", "month", "lat", "lon"), veg_rough),
            "displacement": (("veg_class", "month", "lat", "lon"), displacement),
            "veg_descr": (("veg_class",), np.array(["grass_default"], dtype=object)),
        },
        attrs={"description": "VIC parameter file"}
    )

    ds["lat"].attrs.update({"units": "degrees_north", "long_name": "latitude of grid cell center"})
    ds["lon"].attrs.update({"units": "degrees_east", "long_name": "longitude of grid cell center"})

    ds.to_netcdf(out_path)
    print(f"   ✅ {out_path}")
    return out_path


def write_vic_image_global_param(
    out_path: str,
    domain_nc: str,
    params_nc: str,
    forcing_prefix: str,
    result_dir: str,
    log_dir: str,
) -> str:
    print("\n[STEP] Writing VIC Image global parameter file")
    ensure_dir(os.path.dirname(out_path))
    ensure_dir(result_dir)
    ensure_dir(log_dir)

    text = f"""
#######################################################################
# VIC Image Driver global parameter file
#######################################################################
NLAYER              3
NODES               3
MODEL_STEPS_PER_DAY {VIC_MODEL_STEPS_PER_DAY}
SNOW_STEPS_PER_DAY  {VIC_SNOW_STEPS_PER_DAY}
RUNOFF_STEPS_PER_DAY {VIC_RUNOFF_STEPS_PER_DAY}
STARTYEAR           {YEAR_START}
STARTMONTH          1
STARTDAY            1
STARTSEC            0
ENDYEAR             {YEAR_END}
ENDMONTH            12
ENDDAY              31
FULL_ENERGY         FALSE
FROZEN_SOIL         FALSE
QUICK_FLUX          TRUE
QUICK_SOLVE         TRUE
NOFLUX              FALSE
IMPLICIT            FALSE
EXP_TRANS           TRUE
SNOW_DENSITY        DENS_BRAS
BLOWING             FALSE
COMPUTE_TREELINE    FALSE
CARBON              FALSE
CONTINUEONERROR     TRUE
AERO_RESIST_CANSNOW AR_406_FULL

#######################################################################
# Forcing Files and Parameters
#######################################################################
FORCING1            {forcing_prefix}
FORCE_TYPE          AIR_TEMP  tas
FORCE_TYPE          PREC      prcp
FORCE_TYPE          PRESSURE  pres
FORCE_TYPE          SWDOWN    dswrf
FORCE_TYPE          LWDOWN    dlwrf
FORCE_TYPE          VP        vp
FORCE_TYPE          WIND      wind
CANOPY_LAYERS       3
WIND_H              10.0

#######################################################################
# DOMAIN INFO
#######################################################################
DOMAIN              {domain_nc}
DOMAIN_TYPE         LAT      lat
DOMAIN_TYPE         LON      lon
DOMAIN_TYPE         MASK     mask
DOMAIN_TYPE         AREA     area
DOMAIN_TYPE         FRAC     frac
DOMAIN_TYPE         YDIM     lat
DOMAIN_TYPE         XDIM     lon

#######################################################################
# Land Surface Files and Parameters
#######################################################################
PARAMETERS          {params_nc}
BASEFLOW            ARNO
JULY_TAVG_SUPPLIED  TRUE
ORGANIC_FRACT       FALSE
ALB_SRC             FROM_VEGLIB
LAI_SRC             FROM_VEGLIB
FCAN_SRC            FROM_DEFAULT
SNOW_BAND           FALSE

#######################################################################
# Output Files and Parameters
#######################################################################
LOG_DIR             {log_dir}
RESULT_DIR          {result_dir}

OUTFILE             hydrology
AGGFREQ             NDAYS 1
HISTFREQ            NYEARS 1
OUT_FORMAT          NETCDF4_CLASSIC
OUTVAR              OUT_RUNOFF    *  OUT_TYPE_FLOAT   1  AGG_TYPE_SUM
OUTVAR              OUT_BASEFLOW  *  OUT_TYPE_FLOAT   1  AGG_TYPE_SUM
OUTVAR              OUT_PREC      *  OUT_TYPE_FLOAT   1  AGG_TYPE_SUM
OUTVAR              OUT_EVAP      *  OUT_TYPE_FLOAT   1  AGG_TYPE_SUM
""".strip() + "\n"

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"   ✅ {out_path}")
    return out_path

# =========================================================
# 6. ERA5 -> VIC Image forcing NetCDF
# =========================================================
def prepare_vic_forcing_from_era5(
    vic_grid: Dict[str, np.ndarray],
    out_dir: str,
    prefix: str = "vic_forcing_",
    source_mode: str = "annual"
) -> str:
    print("\n=======================================================")
    if source_mode == "daily":
        print("[PHASE 3] Preparing VIC forcing NetCDF from daily ERA5 files")
    else:
        print("[PHASE 3] Preparing VIC forcing NetCDF from ERA5")
    print("=======================================================")
    ensure_dir(out_dir)

    target_lat = xr.DataArray(vic_grid["lat"], dims=("lat",), coords={"lat": vic_grid["lat"]})
    target_lon = xr.DataArray(vic_grid["lon"], dims=("lon",), coords={"lon": vic_grid["lon"]})

    for year in range(YEAR_START, YEAR_END + 1):
        if source_mode == "daily":
            print(f"   >>> ERA5 daily source root: {ERA5_DAILY_ROOT}/{year:04d}")
            ds = open_era5_daily_dataset_for_year(year)
        else:
            in_nc = ERA5_INPUT_TEMPLATE.format(year=year)
            print(f"   >>> ERA5 source: {in_nc}")
            ds = open_year_dataset(in_nc)

        required = list(ERA5_VARMAP.values())
        for v in required:
            if v not in ds.data_vars:
                raise KeyError(f"❌ ERA5 variable not found for year {year}: {v}")

        # subset bbox first
        lat_min, lat_max = vic_grid["lat"].min() - 0.5, vic_grid["lat"].max() + 0.5
        lon_min, lon_max = vic_grid["lon"].min() - 0.5, vic_grid["lon"].max() + 0.5
        ds = ds.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))

        # optional resample to daily
        if ERA5_FORCE_TO_DAILY:
            ds_daily = xr.Dataset()
            for key, varname in ERA5_VARMAP.items():
                da = ds[varname]
                if key in ["tp", "ssrd", "strd"]:
                    ds_daily[varname] = da.resample(time="1D").sum(keep_attrs=True)
                else:
                    ds_daily[varname] = da.resample(time="1D").mean(keep_attrs=True)
            ds = ds_daily

        # interp to VIC grid
        ds_i = ds.interp(lat=target_lat, lon=target_lon, method="linear")

        # infer timestep
        steps_per_day = infer_steps_per_day(pd.DatetimeIndex(ds_i["time"].values))
        if steps_per_day != VIC_MODEL_STEPS_PER_DAY:
            raise ValueError(
                f"❌ ERA5/VIC timestep mismatch: forcing={steps_per_day} steps/day, "
                f"but VIC_MODEL_STEPS_PER_DAY={VIC_MODEL_STEPS_PER_DAY}. "
                f"请修改 VIC_MODEL_STEPS_PER_DAY，或者修改 ERA5_FORCE_TO_DAILY。"
            )
        seconds_per_step = int(round(86400 / steps_per_day))

        tas = ds_i[ERA5_VARMAP["t2m"]] - 273.15
        prcp = ds_i[ERA5_VARMAP["tp"]] * (1000.0 if ERA5_PREC_IN_M else 1.0)
        pres = ds_i[ERA5_VARMAP["sp"]] / 1000.0   # Pa -> kPa
        td_c = ds_i[ERA5_VARMAP["d2m"]] - 273.15
        vp = dewpoint_to_vapor_pressure_kpa(td_c)
        wind = np.sqrt(ds_i[ERA5_VARMAP["u10"]] ** 2 + ds_i[ERA5_VARMAP["v10"]] ** 2)

        if ERA5_RAD_IS_ACCUM_J:
            dswrf = ds_i[ERA5_VARMAP["ssrd"]] / seconds_per_step
            dlwrf = ds_i[ERA5_VARMAP["strd"]] / seconds_per_step
        else:
            dswrf = ds_i[ERA5_VARMAP["ssrd"]]
            dlwrf = ds_i[ERA5_VARMAP["strd"]]

        out = xr.Dataset(
            data_vars={
                "tas": tas.astype(np.float32),
                "prcp": prcp.astype(np.float32),
                "pres": pres.astype(np.float32),
                "dswrf": dswrf.astype(np.float32),
                "dlwrf": dlwrf.astype(np.float32),
                "vp": vp.astype(np.float32),
                "wind": wind.astype(np.float32),
            },
            coords={
                "time": ds_i["time"],
                "lat": target_lat,
                "lon": target_lon,
            },
            attrs={
                "title": f"VIC Image forcing from ERA5, {year}",
                "Conventions": "CF-1.6"
            }
        )
        out["prcp"].attrs.update({"long_name": "PREC", "units": "mm"})
        out["tas"].attrs.update({"long_name": "AIR_TEMP", "units": "C"})
        out["dswrf"].attrs.update({"long_name": "SWDOWN", "units": "W/m2"})
        out["dlwrf"].attrs.update({"long_name": "LWDOWN", "units": "W/m2"})
        out["pres"].attrs.update({"long_name": "PRESSURE", "units": "kPa"})
        out["vp"].attrs.update({"long_name": "VP", "units": "kPa"})
        out["wind"].attrs.update({"long_name": "WIND", "units": "m/s"})

        out_nc = os.path.join(out_dir, f"{prefix}{year}.nc")
        out.to_netcdf(out_nc)
        print(f"   ✅ Wrote VIC forcing: {out_nc}")

    return os.path.join(out_dir, prefix)


def resolve_vic_forcing_prefix(vic_grid: Dict[str, np.ndarray]) -> str:
    ensure_dir(VIC_FORCING_DIR)
    if VIC_FORCING_MODE == "READY_VIC_NC":
        print("\n[INFO] Using prebuilt VIC NetCDF forcing prefix")
        test_file = f"{VIC_READY_FORCING_PREFIX}{YEAR_START}.nc"
        ensure_file_exists(test_file, "VIC forcing NetCDF")
        return VIC_READY_FORCING_PREFIX
    elif VIC_FORCING_MODE == "ERA5_TO_VIC":
        return prepare_vic_forcing_from_era5(vic_grid, VIC_FORCING_DIR, prefix="vic_forcing_", source_mode="annual")
    elif VIC_FORCING_MODE == "ERA5_DAILY_TO_VIC":
        return prepare_vic_forcing_from_era5(vic_grid, VIC_FORCING_DIR, prefix="vic_forcing_", source_mode="daily")
    else:
        raise ValueError(f"❌ Unsupported VIC_FORCING_MODE: {VIC_FORCING_MODE}")

# =========================================================
# 7. 运行 VIC Image Driver
# =========================================================
def run_vic_image(global_param_path: str) -> None:
    print("\n=======================================================")
    print("[PHASE 4] Running VIC Image Driver")
    print("=======================================================")
    ensure_file_exists(VIC_EXEC, "VIC Image executable")
    run_cmd(f'{VIC_EXEC} -g {global_param_path}', cwd=VIC_WORKDIR)
    print("   ✅ VIC Image finished")

# =========================================================
# 8. VIC 输出 -> CaMa-Flood 年 runoff NetCDF
# =========================================================
def convert_vic_output_to_cama_runoff_nc(
    vic_result_dir: str,
    out_dir: str,
    prefix: str = "vic_runoff_",
    runoff_var_name: str = "Runoff",
) -> List[str]:
    print("\n=======================================================")
    print("[PHASE 5] Converting VIC output to CaMa-Flood runoff NetCDF")
    print("=======================================================")
    ensure_dir(out_dir)

    nc_files = sorted(glob.glob(os.path.join(vic_result_dir, "*.nc")))
    if not nc_files:
        raise FileNotFoundError(f"❌ No VIC NetCDF outputs found in {vic_result_dir}")

    target_files = []
    for f in nc_files:
        try:
            with xr.open_dataset(f) as ds0:
                if "OUT_RUNOFF" in ds0.data_vars and "OUT_BASEFLOW" in ds0.data_vars:
                    target_files.append(f)
        except Exception:
            pass

    if not target_files:
        raise RuntimeError("❌ Could not find VIC output files containing OUT_RUNOFF and OUT_BASEFLOW")

    ds = xr.open_mfdataset(target_files, combine="by_coords")
    runoff = ds["OUT_RUNOFF"] + ds["OUT_BASEFLOW"]

    # 如有必要，聚合为 daily mm/day
    steps_per_day = infer_steps_per_day(pd.DatetimeIndex(ds["time"].values))
    if steps_per_day > 1:
        runoff = runoff.resample(time="1D").sum()
    else:
        runoff = runoff.resample(time="1D").sum()

    out_files = []
    for year in range(YEAR_START, YEAR_END + 1):
        sel = runoff.sel(time=slice(f"{year}-01-01", f"{year}-12-31"))
        out_ds = xr.Dataset(
            data_vars={
                runoff_var_name: (sel.dims, sel.values.astype(np.float32), {
                    "long_name": "VIC total runoff for CaMa-Flood (OUT_RUNOFF + OUT_BASEFLOW)",
                    "units": "mm/day"
                })
            },
            coords={k: sel.coords[k] for k in sel.coords},
            attrs={
                "title": f"VIC runoff for CaMa-Flood, {REGION_NAME}, {year}",
                "note": "Runoff = OUT_RUNOFF + OUT_BASEFLOW"
            }
        )
        out_nc = os.path.join(out_dir, f"{prefix}{year}.nc")
        out_ds.to_netcdf(out_nc)
        out_files.append(out_nc)
        print(f"   ✅ Wrote {out_nc}")

    ds.close()
    return out_files

# =========================================================
# 9. 为 VIC runoff 重新生成 CaMa input matrix
# =========================================================
def prepare_cama_inpmat_for_vic(reg_map_dir: str, bbox: Tuple[int, int, int, int], vic_res_deg: float) -> Tuple[str, str]:
    print("\n=======================================================")
    print("[PHASE 6] Preparing CaMa input matrix for VIC runoff")
    print("=======================================================")

    src_param_dir = os.path.join(reg_map_dir, "src_param")
    if not os.path.exists(src_param_dir):
        shutil.copytree(os.path.join(CAMA_DIR, "map", "src", "src_param"), src_param_dir)
    run_cmd("make all", cwd=src_param_dir)

    script_s02_inp = os.path.join(src_param_dir, "s02-generate_inpmat.sh")
    subprocess.run(["chmod", "+x", script_s02_inp], check=False)

    west, east, south, north = bbox
    nxin = int(round((east - west) / vic_res_deg))
    nyin = int(round((north - south) / vic_res_deg))

    diminfo_name = os.path.join(reg_map_dir, f"diminfo_vic_{VIC_RESOLUTION}.txt")
    inpmat_name = os.path.join(reg_map_dir, f"inpmat_vic_{VIC_RESOLUTION}.bin")

    replacements = {
        "WEST": str(west),
        "EAST": str(east),
        "SOUTH": str(south),
        "NORTH": str(north),
        "GSIZE": str(vic_res_deg),
        "DLON": str(vic_res_deg),
        "DLAT": str(vic_res_deg),
        "RES": str(vic_res_deg),
        "LONRES": str(vic_res_deg),
        "LATRES": str(vic_res_deg),
        "NXIN": str(nxin),
        "NYIN": str(nyin),
        "IX": str(nxin),
        "IY": str(nyin),
        # VIC forcing lat 按升序写入 netCDF
        "LATREV": "FALSE",
        "LAT_REVERSE": "FALSE",
        "LREVERSE": "FALSE",
        "DIMINFO": f'"{diminfo_name}"',
        "CDIMINFO": f'"{diminfo_name}"',
        "INPMAT": f'"{inpmat_name}"',
        "CINPMAT": f'"{inpmat_name}"',
    }
    robust_modify_shell_vars(script_s02_inp, replacements)
    run_cmd("./s02-generate_inpmat.sh", cwd=src_param_dir)

    diminfo_path = diminfo_name if os.path.exists(diminfo_name) else None
    inpmat_path = inpmat_name if os.path.exists(inpmat_name) else None
    if diminfo_path is None:
        candidates = sorted(glob.glob(os.path.join(reg_map_dir, "*diminfo*.txt")), key=os.path.getmtime)
        if candidates:
            diminfo_path = candidates[-1]
    if inpmat_path is None:
        candidates = sorted(glob.glob(os.path.join(reg_map_dir, "*inpmat*.bin")), key=os.path.getmtime)
        if candidates:
            inpmat_path = candidates[-1]

    if diminfo_path is None or inpmat_path is None:
        raise RuntimeError("❌ Failed to generate CaMa diminfo/inpmat for VIC runoff.")

    ensure_nonempty_file(diminfo_path, "diminfo")
    ensure_nonempty_file(inpmat_path, "inpmat")
    print(f"   ✅ diminfo: {diminfo_path}")
    print(f"   ✅ inpmat : {inpmat_path}")
    return diminfo_path, inpmat_path

# =========================================================
# 10. 生成 CaMa go script
# =========================================================
def prepare_cama_run_script_for_vic(
    diminfo_file_abs: str,
    inpmat_file_abs: str,
    runoff_nc_dir: str,
    runoff_prefix: str = "vic_runoff_",
    runoff_var: str = "Runoff",
) -> str:
    print("\n=======================================================")
    print("[PHASE 7] Preparing CaMa go script driven by VIC runoff")
    print("=======================================================")
    gosh_dir = os.path.join(CAMA_DIR, "gosh")
    ensure_file_exists(gosh_dir, "gosh directory")

    template_scripts = glob.glob(os.path.join(gosh_dir, "test*.sh"))
    if not template_scripts:
        raise FileNotFoundError("❌ Missing CaMaFlood template script: gosh/test*.sh")

    template_script = template_scripts[0]
    target_script = os.path.join(gosh_dir, f"run_{REGION_NAME}_from_vic_image.sh")
    shutil.copy(template_script, target_script)

    replacements = {
        "EXP": f"{REGION_NAME}_from_vic_image",
        "FMAP": f'"${{BASE}}/map/reg_{REGION_NAME}_{CAMA_RESOLUTION}"',
        "YSTA": str(YEAR_START),
        "YEND": str(YEAR_END),
        "LINPCDF": ".TRUE.",
        "CROFDIR": f'"{runoff_nc_dir}"',
        "CRUNOFFDIR": f'"{runoff_nc_dir}"',
        "CROFCDF": f'"{runoff_prefix}"',
        "CVNROF": f'"{runoff_var}"',
        "CROFVAR": f'"{runoff_var}"',
        "CDIMINFO": f'"{diminfo_file_abs}"',
        "CINPMAT": f'"{inpmat_file_abs}"',
        "SYEARIN": str(YEAR_START),
        "SMONIN": "1",
        "SDAYIN": "1",
        "SHOURIN": "0",
        "IFRQ_INP": "24",
        "DROFUNIT": "86400000",
        "LOUTCDF": ".TRUE.",
        "LFLDOUT": ".TRUE.",
    }
    robust_modify_shell_vars(target_script, replacements)
    print(f"   ✅ CaMa go script: {target_script}")
    return target_script


def run_cama(script_path: str) -> None:
    print("\n=======================================================")
    print("[PHASE 8] Running CaMa-Flood")
    print("=======================================================")
    gosh_dir = os.path.dirname(script_path)
    run_cmd(f"./{os.path.basename(script_path)}", cwd=gosh_dir)
    print("   ✅ CaMa-Flood finished")

# =========================================================
# 主流程
# =========================================================
def main() -> None:
    ensure_dir(WORK_ROOT)
    ensure_dir(VIC_WORKDIR)
    ensure_dir(VIC_PARAM_DIR)
    ensure_dir(VIC_RESULT_DIR)
    ensure_dir(VIC_LOG_DIR)
    ensure_dir(CAMA_RUNOFF_NC_DIR)

    # 1. CaMa regional map
    reg_map_dir, bbox, cama_res_deg = build_cama_regional_map()

    # 2. VIC 6min grid based on CaMa map extent
    vic_res_deg = get_resolution_deg(VIC_RESOLUTION)
    vic_grid = build_vic_support_grid(
        reg_map_dir=reg_map_dir,
        cama_bbox=bbox,
        cama_res_deg=cama_res_deg,
        vic_bbox=bbox,
        vic_res_deg=vic_res_deg,
    )

    # 3. VIC Image domain / params / forcing / global param
    domain_nc = build_vic_image_domain_nc(vic_grid, os.path.join(VIC_PARAM_DIR, f"{REGION_NAME}.domain.nc"))
    params_nc = build_vic_image_params_nc(vic_grid, os.path.join(VIC_PARAM_DIR, f"{REGION_NAME}.params.nc"))
    forcing_prefix = resolve_vic_forcing_prefix(vic_grid)
    global_param = write_vic_image_global_param(
        out_path=os.path.join(VIC_WORKDIR, f"global_{REGION_NAME}_image.txt"),
        domain_nc=domain_nc,
        params_nc=params_nc,
        forcing_prefix=forcing_prefix,
        result_dir=VIC_RESULT_DIR,
        log_dir=VIC_LOG_DIR,
    )

    # 4. run VIC image if requested
    if RUN_VIC:
        run_vic_image(global_param)
    else:
        print("\n[INFO] RUN_VIC = False -> 当前仅准备 VIC Image Driver 输入，不运行 VIC")

    # 5. convert VIC output to CaMa runoff nc
    # 如果 RUN_VIC=False，这里要求你之前已经在 VIC_RESULT_DIR 放好 VIC 输出 nc
    runoff_files = convert_vic_output_to_cama_runoff_nc(
        vic_result_dir=VIC_RESULT_DIR,
        out_dir=CAMA_RUNOFF_NC_DIR,
        prefix="vic_runoff_",
        runoff_var_name="Runoff",
    )

    # 6. regenerate CaMa inpmat for 0.1° runoff
    diminfo_file_abs, inpmat_file_abs = prepare_cama_inpmat_for_vic(reg_map_dir, bbox, vic_res_deg)

    # 7. prepare CaMa run script
    cama_script = None
    if PREPARE_CAMA_RUN_SCRIPT:
        cama_script = prepare_cama_run_script_for_vic(
            diminfo_file_abs=diminfo_file_abs,
            inpmat_file_abs=inpmat_file_abs,
            runoff_nc_dir=CAMA_RUNOFF_NC_DIR,
            runoff_prefix="vic_runoff_",
            runoff_var="Runoff",
        )
    else:
        print("\n[INFO] PREPARE_CAMA_RUN_SCRIPT = False -> Skip preparing CaMa go script")

    # 8. run CaMa if requested
    if RUN_CAMA:
        if cama_script is None:
            raise RuntimeError("❌ RUN_CAMA=True 但未生成 CaMa go script")
        run_cama(cama_script)
    else:
        print("\n[INFO] RUN_CAMA = False -> 当前仅准备 CaMa 输入，不运行 CaMa-Flood")

    print("\n=======================================================")
    print("🎉 NetCDF VIC-Image + CaMa pipeline prepared")
    print(f"🔹 CaMa reg map  : {reg_map_dir}")
    print(f"🔹 VIC domain nc : {domain_nc}")
    print(f"🔹 VIC params nc : {params_nc}")
    print(f"🔹 VIC global    : {global_param}")
    print(f"🔹 VIC result dir: {VIC_RESULT_DIR}")
    print(f"🔹 CaMa runoff nc: {runoff_files[0] if runoff_files else 'N/A'}")
    if cama_script:
        print(f"🔹 CaMa go script: {cama_script}")
    print("=======================================================")


if __name__ == "__main__":
    main()
