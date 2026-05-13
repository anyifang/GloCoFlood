#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# CaMa-Flood + VIC Image Driver 自动耦合建模脚本（Shell 主控版）
# ---------------------------------------------------------------------------
# 说明：
# 1. Shell 负责总控与原生模型调用；NetCDF / GeoJSON / 栅格处理使用内嵌 Python。
# 2. Python 步骤进入独立子 shell：先 module load miniforge，再调用指定解释器；
#    Python 结束后自动退出，不污染后续 VIC / CaMa-Flood 的运行环境。
# 3. VIC 使用 Image Driver + NetCDF forcing / params / outputs。
###############################################################################

############################
# 0. 用户配置区
############################

# 模块环境
MODULE_INIT_SH="/public1/soft/modules/module.sh"
PYTHON_MODULE="miniforge/25.3"
PYTHON_BIN="python"
VIC_RUNTIME_MODULE="netcdf/4.4.1-icc-kd"

# 基础编译/运行库
source "$MODULE_INIT_SH"
module load netcdf/4.4.1-icc-kd
module load mpi/openmpi/4.1.0-gcc8.5.0

# 运行 VIC 时需要的动态库路径
NETCDF_LIB_DIR="/miniforge/25.3/lib"
export HDF5_USE_FILE_LOCKING=FALSE

# 根目录与区域
CAMA_DIR="/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/cmf_v420_pkg"
GEOJSON_PATH="/publicfs01/fs1-m8/home/m8s001451/zayf/HydroMT-SFINCS/examples/PRD_SFINCS/cama_upstream_rivers.geojson"
REGION_NAME="PRD_Auto"
CAMA_RESOLUTION="15min"
VIC_RESOLUTION="15min"

YEAR_START=2000
YEAR_END=2015
BUFFER_DEG=1.0

PREPARE_CAMA_RUN_SCRIPT=1
RUN_VIC=1
RUN_CAMA=1

VIC_EXEC="/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/VIC-master/vic/drivers/image/vic_image.exe"

# ERA5 -> VIC forcing
VIC_FORCING_MODE="ERA5_DAILY_TO_VIC"   # READY_VIC_NC or ERA5_DAILY_TO_VIC; 调试阶段建议切到 READY_VIC_NC 复用已生成 forcing
VIC_READY_FORCING_PREFIX="/publicfs01/fs1-m8/home/m8s001451/zayf/VIC_Data/vic_forcing/vic_forcing_"
ERA5_DAILY_ROOT="/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/cmf_v420_pkg/VIC_Data/ERA5"
ERA5_DAILY_FILE_TEMPLATE="ERA5_VIC_global_{year}_{month}_{day}.nc"
ERA5_EXPECTED_HOURS="0,6,12,18"
ERA5_REQUIRE_COMPLETE_YEAR=1
ERA5_FORCE_TO_DAILY=0
ERA5_PREC_IN_M=1
ERA5_RAD_IS_ACCUM_J=1

# VIC 步长（当前对应每日 4 个 6 小时文件）
VIC_MODEL_STEPS_PER_DAY=4
VIC_SNOW_STEPS_PER_DAY=4
VIC_RUNOFF_STEPS_PER_DAY=4
VIC_RUNOFF_HOURS=$((24 / VIC_RUNOFF_STEPS_PER_DAY))

# 参数数据源
DIR_OPENLANDMAP="/publicfs01/fs1-m8/home/m8s001451/zayf/Camaflood/cmf_v420_pkg/VIC_Data"
LAND_DEM_NC="/publicfs01/fs1-m8/home/m8s001451/zayf/HydroMT-SFINCS/examples/data/gebco/GEBCO_2025_sub_ice.nc"
GLOBAL_SOIL_TXT="${DIR_OPENLANDMAP}/global_soil_param_new.txt"
ANNUAL_PRECIP_TIF="${DIR_OPENLANDMAP}/annual_mean_precip.tif"
GLOBCOVER_TIF="${DIR_OPENLANDMAP}/GlobCover2009_merged.tif"

# 工作目录
WORKNAME="vic_camaflood"
WORK_ROOT="${CAMA_DIR}/gosh/work_${REGION_NAME}_${WORKNAME}"
VIC_WORKDIR="${WORK_ROOT}/vic_image"
VIC_FORCING_DIR="${VIC_WORKDIR}/forcing"
VIC_PARAM_DIR="${VIC_WORKDIR}/params"
VIC_RESULT_DIR="${VIC_WORKDIR}/result"
VIC_LOG_DIR="${VIC_WORKDIR}/log"
CAMA_RUNOFF_NC_DIR="${WORK_ROOT}/cama_runoff_from_vic"

############################
# 1. 通用函数
############################
log() {
  echo "$@"
}

ensure_dir() {
  mkdir -p "$1"
}

ensure_file_exists() {
  local path="$1"
  local desc="${2:-file}"
  if [[ ! -f "$path" ]]; then
    echo "❌ Missing ${desc}: ${path}" >&2
    exit 1
  fi
}

ensure_nonempty_file() {
  local path="$1"
  local desc="${2:-file}"
  ensure_file_exists "$path" "$desc"
  if [[ ! -s "$path" ]]; then
    echo "❌ Empty ${desc}: ${path}" >&2
    exit 1
  fi
}

run_cmd() {
  local cmd="$1"
  local cwd="${2:-}"
  log "   >>> Executing: ${cmd}"
  if [[ -n "$cwd" ]]; then
    (cd "$cwd" && bash -lc "$cmd")
  else
    bash -lc "$cmd"
  fi
}

replace_shell_var() {
  local file="$1"
  local key="$2"
  local value="$3"
  ensure_file_exists "$file" "script"

  if grep -Eq "^[[:space:]]*(export[[:space:]]+)?${key}[[:space:]]*=" "$file"; then
    sed -i -E "s|^([[:space:]]*(export[[:space:]]+)?${key}[[:space:]]*=[[:space:]]*).*$|\\1${value}|" "$file"
    log "      ✅ Updated shell var: ${key}=${value}"
  else
    printf '\n%s=%s\n' "$key" "$value" >> "$file"
    log "      ➕ Added shell var: ${key}=${value}"
  fi
}

# Python 步骤：在子 shell 中加载 miniforge，但固定使用指定解释器
run_python_stdin() {
  (
    source "$MODULE_INIT_SH"
    module load "$PYTHON_MODULE" >/dev/null 2>&1 || true
    "$PYTHON_BIN" -
  )
}

# 运行 VIC / CaMa 前，清理 Python 模块残留
ensure_native_runtime_env() {
  source "$MODULE_INIT_SH"
  module unload "$PYTHON_MODULE" >/dev/null 2>&1 || true
  hash -r
}

get_resolution_deg() {
  case "$1" in
    15min) echo "0.25" ;;
    06min) echo "0.1" ;;
    05min) echo "0.0833333333333333" ;;
    03min) echo "0.05" ;;
    01min) echo "0.0166666666666667" ;;
    *)
      echo "❌ Unsupported resolution: $1" >&2
      exit 1
      ;;
  esac
}

############################
# 2. GeoJSON -> bbox
############################
get_bbox_from_geojson() {
  local geojson_path="$1"
  local buffer_deg="$2"
  run_python_stdin <<PY
import geopandas as gpd
import numpy as np

gdf = gpd.read_file(r"${geojson_path}")
minx, miny, maxx, maxy = gdf.total_bounds
west = int(np.floor(minx - float(${buffer_deg})))
east = int(np.ceil(maxx + float(${buffer_deg})))
south = int(np.floor(miny - float(${buffer_deg})))
north = int(np.ceil(maxy + float(${buffer_deg})))
print(f"{west} {east} {south} {north}")
PY
}

############################
# 3. PHASE 1: 建立 CaMa 区域域
############################

build_cama_regional_map() {
  log ""
  log "======================================================="
  log "[PHASE 1] Building CaMa-Flood regional domain"
  log "======================================================="

  read -r WEST EAST SOUTH NORTH <<< "$(get_bbox_from_geojson "$GEOJSON_PATH" "$BUFFER_DEG")"
  export WEST EAST SOUTH NORTH
  log "   ✅ BBOX = [${WEST}, ${EAST}, ${SOUTH}, ${NORTH}]"

  REG_MAP_DIR="${WORK_ROOT}/cama_map_${REGION_NAME}_${CAMA_RESOLUTION}"
  SRC_REGION_DIR="${REG_MAP_DIR}/src_region"
  SRC_PARAM_DIR="${REG_MAP_DIR}/src_param"
  CAMA_RES_DEG="$(get_resolution_deg "$CAMA_RESOLUTION")"
  export REG_MAP_DIR SRC_REGION_DIR SRC_PARAM_DIR CAMA_RES_DEG

  ensure_dir "$REG_MAP_DIR"

  log ""
  log "[STEP 1/3] Creating regional map for ${REGION_NAME}"
  run_cmd "make" "${CAMA_DIR}/map/src/src_region"
  if [[ ! -d "$SRC_REGION_DIR" ]]; then
    cp -r "${CAMA_DIR}/map/src/src_region" "$SRC_REGION_DIR"
  fi


  local script_s01_reg="${SRC_REGION_DIR}/s01-regional_map.sh"
  local path_str="${CAMA_DIR}/map/glb_${CAMA_RESOLUTION}/"
  replace_shell_var "$script_s01_reg" "SOURCE" "$path_str"
  replace_shell_var "$script_s01_reg" "MAPGLB" "$path_str"
  replace_shell_var "$script_s01_reg" "CDIR" "$path_str"
  replace_shell_var "$script_s01_reg" "WEST" "$WEST"
  replace_shell_var "$script_s01_reg" "EAST" "$EAST"
  replace_shell_var "$script_s01_reg" "SOUTH" "$SOUTH"
  replace_shell_var "$script_s01_reg" "NORTH" "$NORTH"

  run_cmd "./s01-regional_map.sh" "$SRC_REGION_DIR"
  ensure_nonempty_file "${REG_MAP_DIR}/params.txt" "params.txt"
  ensure_nonempty_file "${REG_MAP_DIR}/nextxy.bin" "nextxy.bin"
  log "   ✅ Regional map created and validated"

  log ""
  log "[STEP 2/3] Generating channel parameters"
  run_cmd "make" "${CAMA_DIR}/map/src/src_param"
  if [[ ! -d "$SRC_PARAM_DIR" ]]; then
    cp -r "${CAMA_DIR}/map/src/src_param" "$SRC_PARAM_DIR"
  fi
  replace_shell_var "${SRC_PARAM_DIR}/s01-channel_params.sh" "CROFBIN" "${CAMA_DIR}/map/data/ELSE_GPCC_coastmod_dayclm-1981-2010.one"

  chmod +x "${SRC_PARAM_DIR}/s01-channel_params.sh"
  run_cmd "./s01-channel_params.sh" "$SRC_PARAM_DIR"
  ensure_nonempty_file "${REG_MAP_DIR}/rivwth.bin" "rivwth.bin"
  log "   ✅ Channel parameters generated"


  log ""
  log "[STEP 3/3] Regional map ready"
}

############################
# 4. PHASE 2/3: Python 生成 VIC domain/params/forcing/global
############################
prepare_vic_inputs() {
  log ""
  log "======================================================="
  log "[PHASE 2-3] Building VIC grid/domain/params/forcing/global"
  log "======================================================="

  ensure_dir "$WORK_ROOT"
  ensure_dir "$VIC_WORKDIR"
  ensure_dir "$VIC_FORCING_DIR"
  ensure_dir "$VIC_PARAM_DIR"
  ensure_dir "$VIC_RESULT_DIR"
  ensure_dir "$VIC_LOG_DIR"
  ensure_dir "$CAMA_RUNOFF_NC_DIR"

  export WEST EAST SOUTH NORTH
  export CAMA_DIR REGION_NAME CAMA_RESOLUTION VIC_RESOLUTION WORKNAME
  export YEAR_START YEAR_END
  export REG_MAP_DIR VIC_WORKDIR VIC_FORCING_DIR VIC_PARAM_DIR VIC_RESULT_DIR VIC_LOG_DIR CAMA_RUNOFF_NC_DIR
  export VIC_READY_FORCING_PREFIX VIC_FORCING_MODE
  export ERA5_DAILY_ROOT ERA5_DAILY_FILE_TEMPLATE ERA5_EXPECTED_HOURS ERA5_REQUIRE_COMPLETE_YEAR ERA5_FORCE_TO_DAILY ERA5_PREC_IN_M ERA5_RAD_IS_ACCUM_J
  export VIC_MODEL_STEPS_PER_DAY VIC_SNOW_STEPS_PER_DAY VIC_RUNOFF_STEPS_PER_DAY
  export DIR_OPENLANDMAP LAND_DEM_NC GLOBAL_SOIL_TXT ANNUAL_PRECIP_TIF GLOBCOVER_TIF GEOJSON_PATH

  run_python_stdin <<'PY'
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import xarray as xr
from rasterio.features import geometry_mask
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from rasterio.windows import from_bounds
from scipy.spatial import cKDTree

WEST = int(os.environ["WEST"])
EAST = int(os.environ["EAST"])
SOUTH = int(os.environ["SOUTH"])
NORTH = int(os.environ["NORTH"])
REGION_NAME = os.environ["REGION_NAME"]
WORKNAME = os.environ["WORKNAME"]
REG_MAP_DIR = os.environ["REG_MAP_DIR"]
CAMA_RESOLUTION = os.environ["CAMA_RESOLUTION"]
VIC_RESOLUTION = os.environ["VIC_RESOLUTION"]
YEAR_START = int(os.environ["YEAR_START"])
YEAR_END = int(os.environ["YEAR_END"])
VIC_WORKDIR = os.environ["VIC_WORKDIR"]
VIC_FORCING_DIR = os.environ["VIC_FORCING_DIR"]
VIC_PARAM_DIR = os.environ["VIC_PARAM_DIR"]
VIC_RESULT_DIR = os.environ["VIC_RESULT_DIR"]
VIC_LOG_DIR = os.environ["VIC_LOG_DIR"]
VIC_FORCING_MODE = os.environ["VIC_FORCING_MODE"]
VIC_READY_FORCING_PREFIX = os.environ["VIC_READY_FORCING_PREFIX"]
ERA5_DAILY_ROOT = os.environ["ERA5_DAILY_ROOT"]
ERA5_DAILY_FILE_TEMPLATE = os.environ["ERA5_DAILY_FILE_TEMPLATE"]
ERA5_REQUIRE_COMPLETE_YEAR = bool(int(os.environ["ERA5_REQUIRE_COMPLETE_YEAR"]))
ERA5_FORCE_TO_DAILY = bool(int(os.environ["ERA5_FORCE_TO_DAILY"]))
ERA5_PREC_IN_M = bool(int(os.environ["ERA5_PREC_IN_M"]))
ERA5_RAD_IS_ACCUM_J = bool(int(os.environ["ERA5_RAD_IS_ACCUM_J"]))
VIC_MODEL_STEPS_PER_DAY = int(os.environ["VIC_MODEL_STEPS_PER_DAY"])
VIC_SNOW_STEPS_PER_DAY = int(os.environ["VIC_SNOW_STEPS_PER_DAY"])
VIC_RUNOFF_STEPS_PER_DAY = int(os.environ["VIC_RUNOFF_STEPS_PER_DAY"])
DIR_OPENLANDMAP = os.environ["DIR_OPENLANDMAP"]
LAND_DEM_NC = os.environ["LAND_DEM_NC"]
GLOBAL_SOIL_TXT = os.environ["GLOBAL_SOIL_TXT"]
ANNUAL_PRECIP_TIF = os.environ["ANNUAL_PRECIP_TIF"]
GLOBCOVER_TIF = os.environ["GLOBCOVER_TIF"]
GEOJSON_PATH = os.environ["GEOJSON_PATH"]

ERA5_VARMAP = {
    "t2m": "2m_temperature",
    "tp": "total_precipitation",
    "sp": "surface_pressure",
    "ssrd": "surface_solar_radiation_downwards",
    "strd": "surface_thermal_radiation_downwards",
    "d2m": "2m_dewpoint_temperature",
    "u10": "10m_u_component_of_wind",
    "v10": "10m_v_component_of_wind",
}

SOIL_LUT = {
    1:  [0.36, 0.17, 127.057, 27.56, 0.098],
    2:  [0.37, 0.25, 96.161,  22.52, 0.103],
    3:  [0.31, 0.23, 208.449, 29.00, 0.098],
    4:  [0.34, 0.21, 113.240, 19.04, 0.081],
    5:  [0.36, 0.21, 116.950, 17.96, 0.086],
    6:  [0.27, 0.17, 146.218, 20.32, 0.060],
    7:  [0.29, 0.14, 123.880, 13.60, 0.088],
    8:  [0.32, 0.12,  64.863, 10.58, 0.140],
    9:  [0.21, 0.09, 291.743, 12.68, 0.057],
    10: [0.28, 0.08, 306.902,  9.10, 0.071],
    11: [0.15, 0.06, 895.365, 10.98, 0.050],
    12: [0.08, 0.03, 5794.287, 11.20, 0.050],
}
VIC_VEG_CLASS_IDS = np.arange(1, 12, dtype=np.int32)
VIC_VEG_DESCR = np.array([
    "Needleleaf",
    "Broadleaf",
    "Deciduous_Needleleaf",
    "Deciduous_Broadleaf",
    "Mixed_Cover",
    "Woodland",
    "Wooded_Grasslands",
    "Closed_Shrublands",
    "Open_Shrublands",
    "Grasslands",
    "Cropland",
], dtype=object)
VIC_VEG_LAI = np.array([
    [3.400, 3.400, 3.500, 3.700, 4.000, 4.400, 4.400, 4.300, 4.200, 3.700, 3.500, 3.400],
    [3.400, 3.400, 3.500, 3.700, 4.000, 4.400, 4.400, 4.300, 4.200, 3.700, 3.500, 3.400],
    [1.680, 1.520, 1.680, 2.900, 4.900, 5.000, 5.000, 4.600, 3.440, 3.040, 2.160, 2.000],
    [1.680, 1.520, 1.680, 2.900, 4.900, 5.000, 5.000, 4.600, 3.440, 3.040, 2.160, 2.000],
    [1.680, 1.520, 1.680, 2.900, 4.900, 5.000, 5.000, 4.600, 3.440, 3.040, 2.160, 2.000],
    [1.680, 1.520, 1.680, 2.900, 4.900, 5.000, 5.000, 4.600, 3.440, 3.040, 2.160, 2.000],
    [2.000, 2.250, 2.950, 3.850, 3.750, 3.500, 3.550, 3.200, 3.300, 2.850, 2.600, 2.200],
    [2.000, 2.250, 2.950, 3.850, 3.750, 3.500, 3.550, 3.200, 3.300, 2.850, 2.600, 2.200],
    [2.000, 2.250, 2.950, 3.850, 3.750, 3.500, 3.550, 3.200, 3.300, 2.850, 2.600, 2.200],
    [2.000, 2.250, 2.950, 3.850, 3.750, 3.500, 3.550, 3.200, 3.300, 2.850, 2.600, 2.200],
    [0.500, 0.500, 0.500, 0.500, 1.500, 3.000, 4.500, 5.000, 2.500, 0.500, 0.500, 0.020],
], dtype=np.float64)
VIC_VEG_ROOT_DEPTH = np.array([
    [0.10, 1.00, 5.00],
    [0.10, 1.00, 5.00],
    [0.10, 1.00, 5.00],
    [0.10, 1.00, 5.00],
    [0.10, 1.00, 5.00],
    [0.10, 1.00, 1.00],
    [0.10, 1.00, 1.00],
    [0.10, 1.00, 0.50],
    [0.10, 1.00, 0.50],
    [0.10, 1.00, 0.50],
    [0.10, 0.75, 0.50],
], dtype=np.float64)
VIC_VEG_ROOT_FRACT = np.array([
    [0.05, 0.45, 0.50],
    [0.05, 0.45, 0.50],
    [0.05, 0.45, 0.50],
    [0.05, 0.45, 0.50],
    [0.05, 0.45, 0.50],
    [0.10, 0.65, 0.25],
    [0.10, 0.65, 0.25],
    [0.10, 0.65, 0.25],
    [0.10, 0.65, 0.25],
    [0.10, 0.70, 0.20],
    [0.10, 0.60, 0.30],
], dtype=np.float64)
VIC_VEG_OVERSTORY = np.array([1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=np.int32)
DEFAULT_ALBEDO = np.array([0.20] * 12, dtype=np.float64)
# 按 VIC 文档 rough≈0.123*height, displacement≈0.67*height，以下采用工程经验高度
VIC_VEG_HEIGHT = np.array([20.0, 20.0, 18.0, 18.0, 18.0, 10.0, 8.0, 1.5, 1.0, 0.5, 1.0], dtype=np.float64)
VIC_VEG_ROUGH_MONTHLY = 0.123 * VIC_VEG_HEIGHT[:, None] * np.ones((1, 12), dtype=np.float64)
VIC_VEG_DISPLACEMENT_MONTHLY = 0.67 * VIC_VEG_HEIGHT[:, None] * np.ones((1, 12), dtype=np.float64)
VIC_VEG_RARC = np.full(len(VIC_VEG_CLASS_IDS), 2.0, dtype=np.float64)
VIC_VEG_RMIN = np.full(len(VIC_VEG_CLASS_IDS), 100.0, dtype=np.float64)
VIC_VEG_WIND_H = np.maximum(10.0, VIC_VEG_HEIGHT + 2.0)
# VIC 文档：trees约30 W/m2，crops约100 W/m2；此处将 woody classes 设为30，其余设为100
VIC_VEG_RGL = np.array([30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 100.0, 100.0, 100.0, 100.0], dtype=np.float64)
VIC_VEG_RAD_ATTEN = np.full(len(VIC_VEG_CLASS_IDS), 0.5, dtype=np.float64)
VIC_VEG_WIND_ATTEN = np.full(len(VIC_VEG_CLASS_IDS), 0.5, dtype=np.float64)
VIC_VEG_TRUNK_RATIO = np.array([0.3, 0.3, 0.3, 0.3, 0.3, 0.25, 0.2, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
VIC_FALLBACK_CLASS = 10
GLOBCOVER_TO_VIC = np.zeros(256, dtype=np.int32)
GLOBCOVER_TO_VIC[11] = 11
GLOBCOVER_TO_VIC[14] = 11
GLOBCOVER_TO_VIC[20] = 11
GLOBCOVER_TO_VIC[30] = 7
GLOBCOVER_TO_VIC[40] = 2
GLOBCOVER_TO_VIC[50] = 4
GLOBCOVER_TO_VIC[60] = 4
GLOBCOVER_TO_VIC[70] = 1
GLOBCOVER_TO_VIC[90] = 1
GLOBCOVER_TO_VIC[100] = 5
GLOBCOVER_TO_VIC[110] = 6
GLOBCOVER_TO_VIC[120] = 7
GLOBCOVER_TO_VIC[130] = 8
GLOBCOVER_TO_VIC[140] = 10
GLOBCOVER_TO_VIC[150] = 9
GLOBCOVER_TO_VIC[160] = 2
GLOBCOVER_TO_VIC[170] = 2
GLOBCOVER_TO_VIC[180] = 7


def get_resolution_deg(resolution: str) -> float:
    return {
        "15min": 15.0 / 60.0,
        "06min": 6.0 / 60.0,
        "05min": 5.0 / 60.0,
        "03min": 3.0 / 60.0,
        "01min": 1.0 / 60.0,
    }[resolution]


def infer_coord_name(ds: xr.Dataset, candidates):
    lower_map = {name.lower(): name for name in list(ds.coords) + list(ds.dims)}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def sample_raster(tif_path: str, lons, lats, multiplier: float = 1.0, nodata_val=None) -> np.ndarray:
    if not os.path.exists(tif_path):
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


def sample_netcdf_dem(nc_path: str, lons, lats, var_name=None, method="nearest") -> np.ndarray:
    if not os.path.exists(nc_path):
        return np.full(len(lons), np.nan, dtype=float)
    with xr.open_dataset(nc_path, engine="netcdf4") as ds:
        if var_name is None:
            var_name = "elevation" if "elevation" in ds.data_vars else list(ds.data_vars)[0]
        lon_name = infer_coord_name(ds, ["lon", "longitude", "x"])
        lat_name = infer_coord_name(ds, ["lat", "latitude", "y"])
        lon_da = xr.DataArray(np.asarray(lons, dtype=float), dims="points")
        lat_da = xr.DataArray(np.asarray(lats, dtype=float), dims="points")
        vals = ds[var_name].sel({lon_name: lon_da, lat_name: lat_da}, method=method).values
        return np.asarray(vals, dtype=float)


def fill_nan_by_nearest_xy(values, lons, lats, default_value):
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
    R = 6371000.0
    lat1 = np.deg2rad(lat_center_deg - dlat_deg / 2.0)
    lat2 = np.deg2rad(lat_center_deg + dlat_deg / 2.0)
    dlon = np.deg2rad(dlon_deg)
    area = (R ** 2) * dlon * (np.sin(lat2) - np.sin(lat1))
    return np.abs(area)


def build_globcover_fractional_cv(globcover_path: str, west: float, east: float, south: float, north: float, vic_res_deg: float, run_idx: np.ndarray) -> np.ndarray:
    """
    只读取模型范围内的 GlobCover，并按最终 VIC 11 类直接聚合为 fractional Cv。
    这样比逐个 GlobCover code 重投影快很多，并且最后强制每个活动格点 sum(Cv)=1.0。
    """
    nveg = len(VIC_VEG_CLASS_IDS)
    ny, nx = run_idx.shape
    cv = np.zeros((nveg, ny, nx), dtype=np.float64)

    if not os.path.exists(globcover_path):
        raise FileNotFoundError(f"Missing GlobCover raster: {globcover_path}")

    dst_transform = from_origin(west, north, vic_res_deg, vic_res_deg)
    dst_shape = (ny, nx)

    with rasterio.open(globcover_path) as src:
        # 只裁剪模型范围，避免整张全球 tif 参与后续重投影
        window = from_bounds(west, south, east, north, src.transform)
        window = window.round_offsets().round_lengths()

        src_arr = src.read(1, window=window, masked=True)
        src_transform = src.window_transform(window)

        src_mask = np.ma.getmaskarray(src_arr)
        src_data = np.asarray(src_arr.filled(0), dtype=np.int16)
        if src.nodata is not None:
            src_mask |= (src_data == src.nodata)

        # 直接按最终 VIC 植被类聚合，最多只做 11 次 reproject
        for vic_cls in VIC_VEG_CLASS_IDS:
            gc_codes = np.where(GLOBCOVER_TO_VIC == vic_cls)[0]
            if len(gc_codes) == 0:
                continue

            src_binary = np.where(np.isin(src_data, gc_codes) & (~src_mask), 1.0, 0.0).astype(np.float32)
            if np.all(src_binary == 0.0):
                continue

            dst_frac = np.zeros(dst_shape, dtype=np.float32)
            reproject(
                source=src_binary,
                destination=dst_frac,
                src_transform=src_transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs="EPSG:4326",
                resampling=Resampling.average,
                src_nodata=0.0,
                dst_nodata=0.0,
            )
            cv[vic_cls - 1, :, :] = dst_frac

    cv[:, ~run_idx] = 0.0

    # 去掉极小碎片覆盖
    cv[cv < 0.01] = 0.0

    # 常规归一化
    totals = cv.sum(axis=0)
    valid = run_idx & (totals > 0)
    if np.any(valid):
        cv[:, valid] = cv[:, valid] / totals[valid]

    # 全空格点回退到 Grasslands
    fallback = run_idx & (~valid)
    if np.any(fallback):
        cv[VIC_FALLBACK_CLASS - 1, fallback] = 1.0

    # 再做一次精确闭合，保证每个活动格点 sum(Cv)=1.0
    iy_ix = np.argwhere(run_idx)
    for iy, ix in iy_ix:
        cell = cv[:, iy, ix].copy()

        if np.all(cell <= 0.0):
            cell[:] = 0.0
            cell[VIC_FALLBACK_CLASS - 1] = 1.0
            cv[:, iy, ix] = cell
            continue

        s = cell.sum()
        if s > 0.0:
            cell = cell / s

        pos = np.where(cell > 0.0)[0]
        if len(pos) == 0:
            cell[:] = 0.0
            cell[VIC_FALLBACK_CLASS - 1] = 1.0
        else:
            last = pos[-1]
            others = cell.sum() - cell[last]
            cell[last] = max(0.0, 1.0 - others)

        s2 = cell.sum()
        if s2 > 0.0 and abs(s2 - 1.0) > 1e-14:
            cell = cell / s2
            pos = np.where(cell > 0.0)[0]
            last = pos[-1]
            others = cell.sum() - cell[last]
            cell[last] = max(0.0, 1.0 - others)

        cv[:, iy, ix] = cell

    cv[:, ~run_idx] = 0.0
    return cv


def open_year_dataset_from_ds(ds: xr.Dataset) -> xr.Dataset:
    time_name = infer_coord_name(ds, ["time"])
    lat_name = infer_coord_name(ds, ["lat", "latitude", "y"])
    lon_name = infer_coord_name(ds, ["lon", "longitude", "x"])
    rename_map = {}
    if time_name != "time":
        rename_map[time_name] = "time"
    if lat_name != "lat":
        rename_map[lat_name] = "lat"
    if lon_name != "lon":
        rename_map[lon_name] = "lon"
    if rename_map:
        ds = ds.rename(rename_map)

    # 保持原有读取逻辑，只在读取后统一坐标
    if np.any(np.diff(ds["lat"].values) < 0):
        ds = ds.sortby("lat")

    # ERA5 若为 0~360 经度，这里统一转为 -180~180
    lon = np.asarray(ds["lon"].values, dtype=np.float64)
    lon = ((lon + 180.0) % 360.0) - 180.0
    ds = ds.assign_coords(lon=lon)

    # 经度转换后，连同数据矩阵一起按经度重新排序
    if np.any(np.diff(ds["lon"].values) < 0):
        ds = ds.sortby("lon")

    # 防止归一化后出现重复经度（如 0/360）
    lon_round = np.round(np.asarray(ds["lon"].values, dtype=np.float64), 10)
    _, unique_idx = np.unique(lon_round, return_index=True)
    if len(unique_idx) < ds["lon"].size:
        ds = ds.isel(lon=np.sort(unique_idx))

    return ds


def get_era5_daily_file(year: int, month: int, day: int) -> str:
    year_dir = os.path.join(ERA5_DAILY_ROOT, f"{year:04d}")
    return os.path.join(year_dir, ERA5_DAILY_FILE_TEMPLATE.format(year=year, month=f"{month:02d}", day=f"{day:02d}"))


def open_era5_daily_dataset_for_year(
    year: int,
    lat_min: float,
    lat_max: float,
    lon_min: float,
    lon_max: float,
) -> xr.Dataset:
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
        raise FileNotFoundError(f"ERA5 daily files missing for year {year}. First missing: {missing[0]}")
    if not files:
        raise FileNotFoundError(f"No ERA5 daily files found for year {year} under {ERA5_DAILY_ROOT}")
    print(f"   >>> Found {len(files)} daily ERA5 files for {year}")
    print(
        f"   >>> Cropping ERA5 on read: "
        f"lat=[{lat_min:.3f}, {lat_max:.3f}], lon=[{lon_min:.3f}, {lon_max:.3f}]"
    )

    keep_vars = list(dict.fromkeys(ERA5_VARMAP.values()))
    datasets = []
    try:
        for fp in files:
            with xr.open_dataset(fp, engine="netcdf4") as ds_one:
                ds_one = open_year_dataset_from_ds(ds_one)
                ds_one = ds_one[keep_vars]
                ds_one = ds_one.sel(
                    lat=slice(lat_min, lat_max),
                    lon=slice(lon_min, lon_max),
                )
                ds_one = ds_one.load()
                datasets.append(ds_one)
        ds = xr.concat(datasets, dim="time")
        ds = ds.sortby("time")
        ds = open_year_dataset_from_ds(ds)
        return ds
    finally:
        datasets.clear()


def infer_steps_per_day(time_index: pd.DatetimeIndex) -> int:
    if len(time_index) < 2:
        return 1
    dt_sec = int((time_index[1] - time_index[0]).total_seconds())
    if dt_sec <= 0:
        return 1
    return max(1, int(round(86400 / dt_sec)))


def dewpoint_to_vapor_pressure_kpa(td_c: xr.DataArray) -> xr.DataArray:
    return 0.6108 * np.exp((17.27 * td_c) / (td_c + 237.3))


# ---------- build vic support grid ----------
cama_res_deg = get_resolution_deg(CAMA_RESOLUTION)
vic_res_deg = get_resolution_deg(VIC_RESOLUTION)
nx_c = int(round((EAST - WEST) / cama_res_deg))
ny_c = int(round((NORTH - SOUTH) / cama_res_deg))
nextxy_path = os.path.join(REG_MAP_DIR, "nextxy.bin")
raw = np.fromfile(nextxy_path, dtype="<i4")
nextxy = raw.reshape(2, ny_c, nx_c)
cama_active = (nextxy[0] != -9999)

nx_v = int(round((EAST - WEST) / vic_res_deg))
ny_v = int(round((NORTH - SOUTH) / vic_res_deg))
lon_v = WEST + (np.arange(nx_v) + 0.5) * vic_res_deg
lat_v = NORTH - (np.arange(ny_v) + 0.5) * vic_res_deg

# ---------- basin polygon mask on VIC grid ----------
gdf_basin = gpd.read_file(GEOJSON_PATH)
if gdf_basin.crs is not None:
    try:
        if gdf_basin.crs.to_epsg() != 4326:
            gdf_basin = gdf_basin.to_crs(4326)
    except Exception:
        gdf_basin = gdf_basin.to_crs(4326)
geoms = [geom for geom in gdf_basin.geometry if geom is not None and not geom.is_empty]
if len(geoms) == 0:
    raise ValueError(f"No valid geometries found in basin polygon file: {GEOJSON_PATH}")
transform_vic = from_origin(WEST, NORTH, vic_res_deg, vic_res_deg)
basin_mask = geometry_mask(
    geoms,
    out_shape=(ny_v, nx_v),
    transform=transform_vic,
    invert=True,
    all_touched=False,
).astype(bool)
lon_c = WEST + (np.arange(nx_c) + 0.5) * cama_res_deg
lat_c = NORTH - (np.arange(ny_c) + 0.5) * cama_res_deg
lon2_c, lat2_c = np.meshgrid(lon_c, lat_c)
counts = np.zeros((ny_v, nx_v), dtype=np.int32)
active_counts = np.zeros((ny_v, nx_v), dtype=np.int32)
ix = np.floor((lon2_c - WEST) / vic_res_deg).astype(int)
iy = np.floor((NORTH - lat2_c) / vic_res_deg).astype(int)
inside = (ix >= 0) & (ix < nx_v) & (iy >= 0) & (iy < ny_v)
np.add.at(counts, (iy[inside], ix[inside]), 1)
np.add.at(active_counts, (iy[inside & cama_active], ix[inside & cama_active]), 1)
frac = np.zeros((ny_v, nx_v), dtype=np.float64)
valid = counts > 0
frac[valid] = active_counts[valid] / counts[valid]
# Only keep VIC cells whose centers fall inside the basin polygon
mask = basin_mask.astype(np.int32)
frac = np.where(mask > 0, frac, 0.0)
area_lat = compute_regular_latlon_cell_area(lat_v, vic_res_deg, vic_res_deg)
area = np.repeat(area_lat[:, None], nx_v, axis=1)
print(f"   ✅ VIC 6min grid: ny={ny_v}, nx={nx_v}")
print(f"   ✅ Basin-masked VIC cells: {int(mask.sum())}")

# ---------- domain.nc ----------
domain_nc = os.path.join(VIC_PARAM_DIR, f"{REGION_NAME}.domain.nc")
ds_domain = xr.Dataset(
    data_vars={
        "mask": (("lat", "lon"), mask.astype(np.int32)),
        "area": (("lat", "lon"), area.astype(np.float64)),
        "frac": (("lat", "lon"), frac.astype(np.float64)),
    },
    coords={"lat": lat_v.astype(np.float64), "lon": lon_v.astype(np.float64)},
    attrs={"description": "VIC Image Driver domain file"},
)
ds_domain.to_netcdf(domain_nc, engine="netcdf4")
print(f"   ✅ {domain_nc}")

# ---------- soil dataframe ----------
lon2, lat2 = np.meshgrid(lon_v, lat_v)
run_idx = mask.astype(bool)
df = pd.DataFrame({
    "lat": lat2[run_idx],
    "lon": lon2[run_idx],
    "gridcell": np.arange(1, int(run_idx.sum()) + 1, dtype=np.int32),
})

for i, depth_cm in enumerate([0, 30, 100], start=1):
    tif_name = f"sol_texture.class_usda.tt_m_250m_b{depth_cm}..{depth_cm}cm_1950..2017_v0.2.tif"
    tif = os.path.join(DIR_OPENLANDMAP, tif_name)
    classes = sample_raster(tif, df["lon"].values, df["lat"].values, nodata_val=255)
    classes = fill_nan_by_nearest_xy(classes, df["lon"].values, df["lat"].values, 1.0)
    df[f"class_{i}"] = np.clip(np.rint(classes), 1, 12).astype(int)

for i, depth_cm in enumerate([0, 30, 100], start=1):
    tif_name = f"sol_bulkdens.fineearth_usda.4a1h_m_250m_b{depth_cm}..{depth_cm}cm_1950..2017_v0.2.tif"
    tif = os.path.join(DIR_OPENLANDMAP, tif_name)
    bd = sample_raster(tif, df["lon"].values, df["lat"].values, multiplier=10.0, nodata_val=255)
    bd[bd < 800] = np.nan
    bd = fill_nan_by_nearest_xy(bd, df["lon"].values, df["lat"].values, 1300.0)
    df[f"bulk_density_{i}"] = bd

elev = sample_netcdf_dem(LAND_DEM_NC, df["lon"].values, df["lat"].values, method="nearest")
elev = fill_nan_by_nearest_xy(elev, df["lon"].values, df["lat"].values, 10.0)
df["elev"] = elev

annual_prec = sample_raster(ANNUAL_PRECIP_TIF, df["lon"].values, df["lat"].values)
annual_prec = fill_nan_by_nearest_xy(annual_prec, df["lon"].values, df["lat"].values, 1000.0)
df["annual_prec"] = annual_prec

soil_cols = ["init_moist_1", "init_moist_2", "init_moist_3", "avg_T", "bubble_1", "bubble_2", "bubble_3", "quartz_1", "quartz_2", "quartz_3"]
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
        df[col] = fill_nan_by_nearest_xy(df[col].values, df["lon"].values, df["lat"].values, 10.0)
else:
    for col in soil_cols:
        df[col] = 10.0

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
df["July_Tavg"] = 1
print(f"   ✅ Sampled {len(df)} active VIC cells")

# ---------- fractional vegetation cover from GlobCover ----------
Cv_frac_all = build_globcover_fractional_cv(
    globcover_path=GLOBCOVER_TIF,
    west=WEST,
    east=EAST,
    south=SOUTH,
    north=NORTH,
    vic_res_deg=vic_res_deg,
    run_idx=run_idx,
)
print("   ✅ Built fractional Cv from GlobCover")

# ---------- params.nc ----------
lat_vals = lat_v.astype(np.float64)
lon_vals = lon_v.astype(np.float64)
ny, nx = len(lat_vals), len(lon_vals)
lat_to_idx = {float(v): i for i, v in enumerate(lat_vals)}
lon_to_idx = {float(v): j for j, v in enumerate(lon_vals)}
fill = np.nan
nlayer, nveg, nroot, nmonth = 3, len(VIC_VEG_CLASS_IDS), 3, 12
run_cell = np.zeros((ny, nx), dtype=np.int32)
gridcell = np.zeros((ny, nx), dtype=np.int32)
lats2 = np.repeat(lat_vals[:, None], nx, axis=1)
lons2 = np.repeat(lon_vals[None, :], ny, axis=0)
def alloc2(): return np.full((ny, nx), fill)
def alloc3(): return np.full((nlayer, ny, nx), fill)
infilt = alloc2(); Ds = alloc2(); Dsmax = alloc2(); Ws = alloc2(); c = alloc2(); elev2 = alloc2(); avg_T = alloc2(); dp2 = alloc2(); off_gmt = alloc2(); rough = alloc2(); snow_rough = alloc2(); annual_prec2 = alloc2(); July_Tavg = alloc2()
expt = alloc3(); Ksat = alloc3(); phi_s = alloc3(); init_moist = alloc3(); depth = alloc3(); bubble = alloc3(); quartz = alloc3(); bulk_density = alloc3(); soil_density = alloc3(); Wcr_FRACT = alloc3(); Wpwp_FRACT = alloc3(); resid_moist = alloc3()
fs_active = np.zeros((ny, nx), dtype=np.int32)
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

for k in range(nveg):
    root_depth[k, :, :, :] = VIC_VEG_ROOT_DEPTH[k][:, None, None]
    root_fract[k, :, :, :] = VIC_VEG_ROOT_FRACT[k][:, None, None]
    LAI[k, :, :, :] = VIC_VEG_LAI[k][:, None, None]
    overstory[k, :, :] = VIC_VEG_OVERSTORY[k]
    rarc[k, :, :] = VIC_VEG_RARC[k]
    rmin[k, :, :] = VIC_VEG_RMIN[k]
    wind_h[k, :, :] = VIC_VEG_WIND_H[k]
    RGL[k, :, :] = VIC_VEG_RGL[k]
    rad_atten[k, :, :] = VIC_VEG_RAD_ATTEN[k]
    wind_atten[k, :, :] = VIC_VEG_WIND_ATTEN[k]
    trunk_ratio[k, :, :] = VIC_VEG_TRUNK_RATIO[k]
    albedo[k, :, :, :] = DEFAULT_ALBEDO[:, None, None]
    veg_rough[k, :, :, :] = VIC_VEG_ROUGH_MONTHLY[k][:, None, None]
    displacement[k, :, :, :] = VIC_VEG_DISPLACEMENT_MONTHLY[k][:, None, None]

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
    elev2[i, j] = row["elev"]
    avg_T[i, j] = row["avg_T"]
    dp2[i, j] = row["dp"]
    off_gmt[i, j] = row["off_gmt"]
    rough[i, j] = row["rough"]
    snow_rough[i, j] = row["snow_rough"]
    annual_prec2[i, j] = row["annual_prec"]
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

    cell_cv = Cv_frac_all[:, i, j].copy()
    if np.all(cell_cv <= 0.0):
        cell_cv[:] = 0.0
        cell_cv[VIC_FALLBACK_CLASS - 1] = 1.0

    # 写入 params.nc 前再做一次精确闭合，避免 VIC 报 Sum of veg tile area fractions != 1.0
    s = cell_cv.sum()
    if s > 0.0:
        cell_cv = cell_cv / s
        pos = np.where(cell_cv > 0.0)[0]
        last = pos[-1]
        others = cell_cv.sum() - cell_cv[last]
        cell_cv[last] = max(0.0, 1.0 - others)

    Cv[:, i, j] = cell_cv
    Nveg_arr[i, j] = int(np.count_nonzero(cell_cv > 1e-12))

params_nc = os.path.join(VIC_PARAM_DIR, f"{REGION_NAME}.params.nc")
ds_params = xr.Dataset(
    coords={
        "lat": lat_vals.astype(np.float64),
        "lon": lon_vals.astype(np.float64),
        "nlayer": np.arange(1, nlayer + 1, dtype=np.int32),
        "veg_class": VIC_VEG_CLASS_IDS.astype(np.int32),
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
        "elev": (("lat", "lon"), elev2),
        "depth": (("nlayer", "lat", "lon"), depth),
        "avg_T": (("lat", "lon"), avg_T),
        "dp": (("lat", "lon"), dp2),
        "bubble": (("nlayer", "lat", "lon"), bubble),
        "quartz": (("nlayer", "lat", "lon"), quartz),
        "bulk_density": (("nlayer", "lat", "lon"), bulk_density),
        "soil_density": (("nlayer", "lat", "lon"), soil_density),
        "off_gmt": (("lat", "lon"), off_gmt),
        "Wcr_FRACT": (("nlayer", "lat", "lon"), Wcr_FRACT),
        "Wpwp_FRACT": (("nlayer", "lat", "lon"), Wpwp_FRACT),
        "rough": (("lat", "lon"), rough),
        "snow_rough": (("lat", "lon"), snow_rough),
        "annual_prec": (("lat", "lon"), annual_prec2),
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
        "veg_descr": (("veg_class",), VIC_VEG_DESCR),
    },
    attrs={"description": "VIC parameter file with fractional Cv from GlobCover"},
)
ds_params.to_netcdf(params_nc, engine="netcdf4")
print(f"   ✅ {params_nc}")

# ---------- forcing ----------
def prepare_vic_forcing_from_era5():
    target_lat = xr.DataArray(lat_v, dims=("lat",), coords={"lat": lat_v})
    target_lon = xr.DataArray(lon_v, dims=("lon",), coords={"lon": lon_v})
    lat_min, lat_max = float(lat_v.min() - 0.5), float(lat_v.max() + 0.5)
    lon_min, lon_max = float(lon_v.min() - 0.5), float(lon_v.max() + 0.5)
    for year in range(YEAR_START, YEAR_END + 1):
        if VIC_FORCING_MODE == "READY_VIC_NC":
            continue
        ds = open_era5_daily_dataset_for_year(
            year=year,
            lat_min=lat_min,
            lat_max=lat_max,
            lon_min=lon_min,
            lon_max=lon_max,
        )
        if ERA5_FORCE_TO_DAILY:
            ds_daily = xr.Dataset()
            for key, varname in ERA5_VARMAP.items():
                da = ds[varname]
                if key in ["tp", "ssrd", "strd"]:
                    ds_daily[varname] = da.resample(time="1D").sum(keep_attrs=True)
                else:
                    ds_daily[varname] = da.resample(time="1D").mean(keep_attrs=True)
            ds = ds_daily
        ds_i = ds.interp(lat=target_lat, lon=target_lon, method="linear")
        steps_per_day = infer_steps_per_day(pd.DatetimeIndex(ds_i["time"].values))
        if steps_per_day != VIC_MODEL_STEPS_PER_DAY:
            raise ValueError(f"ERA5/VIC timestep mismatch: forcing={steps_per_day}, VIC={VIC_MODEL_STEPS_PER_DAY}")
        seconds_per_step = int(round(86400 / steps_per_day))

        tas = ds_i[ERA5_VARMAP["t2m"]] - 273.15
        if ERA5_PREC_IN_M:
            prcp = ds_i[ERA5_VARMAP["tp"]] * 1000.0 * 24 / steps_per_day
        else:
            prcp = ds_i[ERA5_VARMAP["tp"]] * 24 / steps_per_day
        pres = ds_i[ERA5_VARMAP["sp"]] /1000
        td_c = ds_i[ERA5_VARMAP["d2m"]] - 273.15
        vp = dewpoint_to_vapor_pressure_kpa(td_c)
        wind = np.sqrt(ds_i[ERA5_VARMAP["u10"]] ** 2 + ds_i[ERA5_VARMAP["v10"]] ** 2)
        if ERA5_RAD_IS_ACCUM_J:
            dswrf = ds_i[ERA5_VARMAP["ssrd"]] / 3600
            dlwrf = ds_i[ERA5_VARMAP["strd"]] / 3600
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
            coords={"time": ds_i["time"], "lat": target_lat, "lon": target_lon},
            attrs={"title": f"VIC Image forcing from ERA5, {year}", "Conventions": "CF-1.6"},
        )
        out["time"].encoding["calendar"] = "proleptic_gregorian"
        out["time"].encoding["units"] = f"hours since {year}-01-01 00:00:00"
        out_nc = os.path.join(VIC_FORCING_DIR, f"vic_forcing_{year}.nc")
        out.to_netcdf(out_nc, engine="netcdf4")
        print(f"   ✅ Wrote VIC forcing: {out_nc}")

if VIC_FORCING_MODE == "READY_VIC_NC":
    test_file = f"{VIC_READY_FORCING_PREFIX}{YEAR_START}.nc"
    if not os.path.exists(test_file):
        raise FileNotFoundError(f"Missing VIC forcing NetCDF: {test_file}")
    forcing_prefix = VIC_READY_FORCING_PREFIX
else:
    prepare_vic_forcing_from_era5()
    forcing_prefix = os.path.join(VIC_FORCING_DIR, "vic_forcing_")

# ---------- global param ----------
global_param = os.path.join(VIC_WORKDIR, f"global_{REGION_NAME}_image.txt")
text = f"""
#######################################################################
# VIC Image Driver global parameter file
#######################################################################
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
CALENDAR            PROLEPTIC_GREGORIAN
OUT_TIME_UNITS      DAYS

FULL_ENERGY         FALSE
FROZEN_SOIL         FALSE

SNOW_DENSITY        DENS_BRAS
BLOWING             FALSE
COMPUTE_TREELINE    FALSE
CARBON              FALSE
AERO_RESIST_CANSNOW AR_406_FULL

#######################################################################
# Forcing Files and Parameters
#######################################################################
FORCING1            {forcing_prefix}
FORCE_TYPE          AIR_TEMP    tas
FORCE_TYPE          PREC        prcp
FORCE_TYPE          PRESSURE    pres
FORCE_TYPE          SWDOWN      dswrf
FORCE_TYPE          LWDOWN      dlwrf
FORCE_TYPE          VP          vp
FORCE_TYPE          WIND        wind
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
ALB_SRC             FROM_VEGPARAM
LAI_SRC             FROM_VEGPARAM
FCAN_SRC            FROM_DEFAULT
SNOW_BAND           FALSE

#######################################################################
# Output Files and Parameters
#######################################################################
LOG_DIR             {VIC_LOG_DIR}
RESULT_DIR          {VIC_RESULT_DIR}

OUTFILE             hydrology
AGGFREQ             NSTEPS 1
HISTFREQ            NYEARS 1
OUT_FORMAT          NETCDF4_CLASSIC
OUTVAR              OUT_RUNOFF    *  OUT_TYPE_FLOAT   1  AGG_TYPE_SUM
OUTVAR              OUT_BASEFLOW  *  OUT_TYPE_FLOAT   1  AGG_TYPE_SUM
OUTVAR              OUT_PREC      *  OUT_TYPE_FLOAT   1  AGG_TYPE_SUM
OUTVAR              OUT_EVAP      *  OUT_TYPE_FLOAT   1  AGG_TYPE_SUM
""".strip() + "\n"

with open(global_param, "w", encoding="utf-8") as f:
    f.write(text)

print("\n[STEP] Writing VIC Image global parameter file")
print(f"   ✅ {global_param}")
PY
}

############################
# 5. PHASE 4: 运行 VIC
############################
run_vic_image() {
  log ""
  log "======================================================="
  log "[PHASE 4] Running VIC Image Driver"
  log "======================================================="

  ensure_file_exists "$VIC_EXEC" "VIC Image executable"
  ensure_native_runtime_env
  export LD_LIBRARY_PATH="${NETCDF_LIB_DIR}:${LD_LIBRARY_PATH:-}"

  local global_param="${VIC_WORKDIR}/global_${REGION_NAME}_image.txt"
  ensure_file_exists "$global_param" "VIC global parameter file"

  local vic_exec_abs
  local global_param_abs
  vic_exec_abs="$(readlink -f "$VIC_EXEC")"
  global_param_abs="$(readlink -f "$global_param")"

  run_cmd "\"${vic_exec_abs}\" -g \"${global_param_abs}\"" "$VIC_WORKDIR"
  log "   ✅ VIC Image finished"
}

############################
# 6. PHASE 5/6/7: Python 转 runoff + 生成 CaMa 脚本
############################
postprocess_vic_and_prepare_cama() {
  log ""
  log "======================================================="
  log "[PHASE 5-7] Convert VIC runoff and prepare CaMa inputs"
  log "======================================================="

  export CAMA_DIR REGION_NAME CAMA_RESOLUTION VIC_RESOLUTION WORKNAME
  export YEAR_START YEAR_END VIC_RUNOFF_STEPS_PER_DAY
  export REG_MAP_DIR VIC_RESULT_DIR CAMA_RUNOFF_NC_DIR WORK_ROOT
  export PREPARE_CAMA_RUN_SCRIPT RUN_CAMA
  export WEST EAST SOUTH NORTH

  run_python_stdin <<'PY'
import os
import glob
import shutil
import subprocess
import numpy as np
import xarray as xr

CAMA_DIR = os.environ["CAMA_DIR"]
REGION_NAME = os.environ["REGION_NAME"]
WORKNAME = os.environ["WORKNAME"]
CAMA_RESOLUTION = os.environ["CAMA_RESOLUTION"]
VIC_RESOLUTION = os.environ["VIC_RESOLUTION"]
YEAR_START = int(os.environ["YEAR_START"])
YEAR_END = int(os.environ["YEAR_END"])
VIC_RUNOFF_STEPS_PER_DAY = int(os.environ["VIC_RUNOFF_STEPS_PER_DAY"])   # 关键修复
REG_MAP_DIR = os.environ["REG_MAP_DIR"]
VIC_RESULT_DIR = os.environ["VIC_RESULT_DIR"]
CAMA_RUNOFF_NC_DIR = os.environ["CAMA_RUNOFF_NC_DIR"]
WORK_ROOT = os.environ["WORK_ROOT"]
PREPARE_CAMA_RUN_SCRIPT = bool(int(os.environ["PREPARE_CAMA_RUN_SCRIPT"]))
RUN_CAMA = bool(int(os.environ["RUN_CAMA"]))
WEST = int(os.environ["WEST"])
EAST = int(os.environ["EAST"])
SOUTH = int(os.environ["SOUTH"])
NORTH = int(os.environ["NORTH"])


def run_cmd(cmd, cwd=None):
    print(f"   >>> Executing: {cmd}")
    result = subprocess.run(cmd, shell=True, cwd=cwd, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\n cwd={cwd}")


def robust_modify_shell_vars(filepath, replacements):
    import re
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()
    for key, val in replacements.items():
        val = str(val)
        pattern = rf"(?m)^(\s*(?:export\s+)?{re.escape(key)}\s*=\s*).*$"
        if re.search(pattern, content):
            content = re.sub(pattern, lambda m, v=val: f"{m.group(1)}{v}", content)
        else:
            content += f"\n{key}={val}\n"
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)
    subprocess.run(["chmod", "+x", filepath], check=False)


def get_resolution_deg(resolution: str) -> float:
    return {
        "15min": 15.0 / 60.0,
        "06min": 6.0 / 60.0,
        "05min": 5.0 / 60.0,
        "03min": 3.0 / 60.0,
        "01min": 1.0 / 60.0,
    }[resolution]


print("\n=======================================================")
print("[PHASE 5] Converting VIC output to CaMa-Flood runoff NetCDF")
print("=======================================================")

nc_files = sorted(glob.glob(os.path.join(VIC_RESULT_DIR, "*.nc")))
if not nc_files:
    raise FileNotFoundError(f"No VIC NetCDF outputs found in {VIC_RESULT_DIR}")

target_files = []
for f in nc_files:
    try:
        with xr.open_dataset(f, engine="netcdf4") as ds0:
            if "OUT_RUNOFF" in ds0.data_vars and "OUT_BASEFLOW" in ds0.data_vars:
                target_files.append(f)
    except Exception:
        pass

if not target_files:
    raise RuntimeError("Could not find VIC output files containing OUT_RUNOFF and OUT_BASEFLOW")

datasets = []
try:
    for fp in target_files:
        ds_one = xr.open_dataset(fp, engine="netcdf4")
        datasets.append(ds_one)

    ds = xr.concat(datasets, dim="time")
    ds = ds.sortby("time")

    hours_per_step = 24.0 / float(VIC_RUNOFF_STEPS_PER_DAY)
    runoff = (ds["OUT_RUNOFF"] + ds["OUT_BASEFLOW"]) / hours_per_step
finally:
    for ds_one in datasets:
        try:
            ds_one.close()
        except Exception:
            pass

os.makedirs(CAMA_RUNOFF_NC_DIR, exist_ok=True)
out_files = []
for year in range(YEAR_START, YEAR_END + 1):
    sel = runoff.sel(time=slice(f"{year}-01-01", f"{year}-12-31"))
    out_ds = xr.Dataset(
        data_vars={
            "Runoff": (sel.dims, sel.values.astype(np.float32), {
                "long_name": "VIC total runoff for CaMa-Flood (OUT_RUNOFF + OUT_BASEFLOW)",
                "units": "mm/hour"
            })
        },
        coords={k: sel.coords[k] for k in sel.coords},
        attrs={"title": f"VIC runoff for CaMa-Flood, {REGION_NAME}, {year}"},
    )
    out_nc = os.path.join(CAMA_RUNOFF_NC_DIR, f"vic_runoff_{year}.nc")
    out_ds.to_netcdf(out_nc, engine="netcdf4")
    out_files.append(out_nc)
    print(f"   ✅ Wrote {out_nc}")

print("\n=======================================================")
print("[PHASE 6] Preparing CaMa input matrix for VIC runoff")
print("=======================================================")

src_param_dir = os.path.join(REG_MAP_DIR, "src_param")
#run_cmd("make all", cwd=src_param_dir)

script_s02_inp = os.path.join(src_param_dir, "s02-generate_inpmat.sh")
subprocess.run(["chmod", "+x", script_s02_inp], check=False)

vic_res_deg = get_resolution_deg(VIC_RESOLUTION)

diminfo_name = f"diminfo_vic_{VIC_RESOLUTION}.txt"
inpmat_name = f"inpmat_vic_{VIC_RESOLUTION}.bin"

replacements = {
    "GRSIZEIN": str(vic_res_deg),
    "WESTIN": str(WEST),
    "EASTIN": str(EAST),
    "NORTHIN": str(NORTH),
    "SOUTHIN": str(SOUTH),
    "OLAT": '"NtoS"',
    "DIMINFO": f'"{diminfo_name}"',
    "INPMAT": f'"{inpmat_name}"',
}

robust_modify_shell_vars(script_s02_inp, replacements)
run_cmd("./s02-generate_inpmat.sh", cwd=src_param_dir)

if PREPARE_CAMA_RUN_SCRIPT:
    print("\n=======================================================")
    print("[PHASE 7] Preparing CaMa go script driven by VIC runoff")
    print("=======================================================")

    gosh_dir = os.path.join(CAMA_DIR, "gosh")
    template_scripts = glob.glob(os.path.join(gosh_dir, "Autorun_Camaflood_exp.sh"))
    if not template_scripts:
        raise FileNotFoundError("Cannot find template script: Autorun_Camaflood_exp.sh")

    template_script = template_scripts[0]
    target_script = os.path.join(gosh_dir, f"run_{REGION_NAME}_from_vic_image.sh")
    shutil.copy(template_script, target_script)

    input_interval_hours = int(round(24 / VIC_RUNOFF_STEPS_PER_DAY))

    replacements = {
        "EXP": f"work_{REGION_NAME}_{WORKNAME}",
        "FMAP": f'"{REG_MAP_DIR}"',
        "YSTA": str(YEAR_START),
        "YEND": str(YEAR_END),
        "LINPCDF": ".TRUE.",
        "CROFDIR": f'"{CAMA_RUNOFF_NC_DIR}"',
        "CRUNOFFDIR": f'"{CAMA_RUNOFF_NC_DIR}"',
        "CROFPRE": '"vic_runoff_"',
        "CVNROF": '"Runoff"',
   	"CDIMINFO": f'"${{FMAP}}/{diminfo_name}"',
    	"CINPMAT": f'"${{FMAP}}/{inpmat_name}"',
        "SMONIN": "1",
        "SDAYIN": "1",
        "SHOURIN": "0",
        "IFRQ_INP": str(input_interval_hours),   # 关键修复
        "DROFUNIT": "3600000",
    }

    robust_modify_shell_vars(target_script, replacements)
    print(f"   ✅ CaMa go script: {target_script}")

    if RUN_CAMA:
        print("\n=======================================================")
        print("[PHASE 8] Running CaMa-Flood")
        print("=======================================================")
        run_cmd(f"./{os.path.basename(target_script)}", cwd=gosh_dir)
        print("   ✅ CaMa-Flood finished")
PY
}

############################
# 7. 主流程
############################
main() {
  build_cama_regional_map
  prepare_vic_inputs

  if [[ "$RUN_VIC" -eq 1 ]]; then
    run_vic_image
  else
    log ""
    log "[INFO] RUN_VIC = False -> 当前仅准备 VIC Image Driver 输入，不运行 VIC"
    log "[INFO] 如果 RUN_VIC = False，则不会继续做 VIC 输出转 CaMa runoff。"
    exit 0
  fi

  postprocess_vic_and_prepare_cama

  log ""
  log "======================================================="
  log "🎉 Shell pipeline finished"
  log "🔹 CaMa reg map  : ${REG_MAP_DIR}"
  log "🔹 VIC domain nc : ${VIC_PARAM_DIR}/${REGION_NAME}.domain.nc"
  log "🔹 VIC params nc : ${VIC_PARAM_DIR}/${REGION_NAME}.params.nc"
  log "🔹 VIC global    : ${VIC_WORKDIR}/global_${REGION_NAME}_image.txt"
  log "🔹 VIC result dir: ${VIC_RESULT_DIR}"
  log "🔹 CaMa runoff nc: ${CAMA_RUNOFF_NC_DIR}/vic_runoff_${YEAR_START}.nc"
  if [[ "$PREPARE_CAMA_RUN_SCRIPT" -eq 1 ]]; then
    log "🔹 CaMa go script: ${CAMA_DIR}/gosh/run_${REGION_NAME}_from_vic_image.sh"
  fi
  log "======================================================="
}

main "$@"
