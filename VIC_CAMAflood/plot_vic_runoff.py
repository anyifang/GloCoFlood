#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# =========================================================
# 0. 固定配置：直接在当前 hydrology 文件夹运行
# =========================================================
FILE_GLOB = "hydrology.*.nc"
OUTDIR = os.path.join(os.getcwd(), "vic_annual_check")

VAR_PREC = "OUT_PREC"
VAR_RUNOFF = "OUT_RUNOFF"
VAR_BASEFLOW = "OUT_BASEFLOW"
VAR_EVAP = "OUT_EVAP"

# 降雨等值线，重点突出中国常用的 400/800 mm
PREC_CONTOURS_MAIN = [400, 800]
PREC_CONTOURS_ALL = [200, 400, 800, 1200, 1600, 2000, 2400]

RUNOFF_CONTOURS = [50, 100, 200, 400, 800, 1200]
EVAP_CONTOURS = [200, 400, 600, 800, 1000, 1200]
RESID_CONTOURS = [-400, -200, -100, 100, 200, 400]


# =========================================================
# 1. 工具函数
# =========================================================
def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def infer_coord_name(ds, candidates):
    lower_map = {name.lower(): name for name in list(ds.coords) + list(ds.dims)}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def standardize_dataset(ds: xr.Dataset) -> xr.Dataset:
    time_name = infer_coord_name(ds, ["time"])
    lat_name = infer_coord_name(ds, ["lat", "latitude", "y"])
    lon_name = infer_coord_name(ds, ["lon", "longitude", "x"])

    rename_map = {}
    if time_name and time_name != "time":
        rename_map[time_name] = "time"
    if lat_name and lat_name != "lat":
        rename_map[lat_name] = "lat"
    if lon_name and lon_name != "lon":
        rename_map[lon_name] = "lon"

    if rename_map:
        ds = ds.rename(rename_map)

    if "lat" not in ds.coords or "lon" not in ds.coords:
        raise ValueError("无法识别 lat/lon 坐标。")

    if np.any(np.diff(ds["lat"].values) < 0):
        ds = ds.sortby("lat")
    if np.any(np.diff(ds["lon"].values) < 0):
        ds = ds.sortby("lon")

    return ds


def collect_files():
    files = sorted(glob.glob(FILE_GLOB))
    files = [f for f in files if os.path.isfile(f)]
    if not files:
        raise FileNotFoundError(f"当前目录 {os.getcwd()} 下没有找到 {FILE_GLOB}")
    return files


def year_from_filename(fp: str):
    base = os.path.basename(fp)
    # hydrology.2000-01-01-00000.nc
    parts = base.split(".")
    if len(parts) >= 2:
        try:
            return int(parts[1][:4])
        except Exception:
            return None
    return None


def open_one_file(fp: str) -> xr.Dataset:
    ds = xr.open_dataset(fp, engine="netcdf4")
    ds = standardize_dataset(ds)

    needed = [VAR_PREC, VAR_RUNOFF, VAR_BASEFLOW, VAR_EVAP]
    missing = [v for v in needed if v not in ds.data_vars]
    if missing:
        raise KeyError(f"{fp} 缺少变量: {missing}. 可用变量: {list(ds.data_vars)}")

    return ds[needed]


def calc_annual_maps(ds: xr.Dataset):
    """
    OUT_PREC / OUT_RUNOFF / OUT_BASEFLOW / OUT_EVAP
    在当前输出设置下是每个时间步的累计量，沿 time 求和就是年总量（mm/year）。
    """
    P = ds[VAR_PREC].sum(dim="time", skipna=True)
    R = (ds[VAR_RUNOFF] + ds[VAR_BASEFLOW]).sum(dim="time", skipna=True)
    E = ds[VAR_EVAP].sum(dim="time", skipna=True)

    # 关键修复：把流域外的 0 值区域剔除，不参与后续平均
    valid_mask = (P > 0) | (R > 0) | (E > 0)

    P = P.where(valid_mask)
    R = R.where(valid_mask)
    E = E.where(valid_mask)
    RES = (P - R - E).where(valid_mask)

    return P, R, E, RES, valid_mask


def plot_contour_map(
    da: xr.DataArray,
    title: str,
    out_png: str,
    cmap="viridis",
    contour_levels=None,
    main_levels=None,
    cbar_label="mm",
    center_zero=False,
):
    lon = da["lon"].values
    lat = da["lat"].values
    xx, yy = np.meshgrid(lon, lat)
    zz = da.values

    # 跳过全 NaN
    if np.all(~np.isfinite(zz)):
        print(f"[WARN] {title} 全部为 NaN，跳过绘图。")
        return

    if center_zero:
        vmax0 = np.nanmax(np.abs(zz))
        vmin, vmax = -vmax0, vmax0
        cmap = "RdBu_r"
    else:
        vmin = np.nanmin(zz)
        vmax = np.nanmax(zz)

    fig, ax = plt.subplots(figsize=(11, 7))

    cf = ax.contourf(xx, yy, zz, levels=20, cmap=cmap, vmin=vmin, vmax=vmax)
    cb = plt.colorbar(cf, ax=ax)
    cb.set_label(cbar_label)

    if contour_levels is not None and len(contour_levels) > 0:
        valid_levels = [lv for lv in contour_levels if np.nanmin(zz) <= lv <= np.nanmax(zz)]
        if len(valid_levels) > 0:
            cs = ax.contour(xx, yy, zz, levels=valid_levels, colors="k", linewidths=0.8)
            ax.clabel(cs, fmt="%d", fontsize=9)

    if main_levels is not None and len(main_levels) > 0:
        valid_main = [lv for lv in main_levels if np.nanmin(zz) <= lv <= np.nanmax(zz)]
        if len(valid_main) > 0:
            cs2 = ax.contour(xx, yy, zz, levels=valid_main, colors="red", linewidths=1.8)
            ax.clabel(cs2, fmt="%d", fontsize=10, colors="red")

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


# =========================================================
# 2. 主程序
# =========================================================
def main():
    ensure_dir(OUTDIR)

    files = collect_files()
    print("=" * 80)
    print("检测到 hydrology 文件：")
    for f in files:
        print("  ", os.path.basename(f))
    print("=" * 80)

    yearly_records = []

    P_list = []
    R_list = []
    E_list = []
    RES_list = []
    MASK_list = []

    for fp in files:
        year = year_from_filename(fp)
        print(f"\n处理: {os.path.basename(fp)}")

        with open_one_file(fp) as ds:
            P, R, E, RES, valid_mask = calc_annual_maps(ds)

            P_list.append(P)
            R_list.append(R)
            E_list.append(E)
            RES_list.append(RES)
            MASK_list.append(valid_mask)

            p_mean = float(P.mean(skipna=True).item())
            r_mean = float(R.mean(skipna=True).item())
            e_mean = float(E.mean(skipna=True).item())
            res_mean = float(RES.mean(skipna=True).item())
            rc = r_mean / p_mean if abs(p_mean) > 1e-12 else np.nan

            yearly_records.append({
                "year": year,
                "domain_mean_precip_mm": p_mean,
                "domain_mean_runoff_mm": r_mean,
                "domain_mean_evap_mm": e_mean,
                "domain_mean_residual_mm": res_mean,
                "runoff_coefficient": rc,
                "n_valid_cells": int(valid_mask.sum().item()),
            })

            # 每年单独出图
            plot_contour_map(
                P,
                title=f"{year} Annual Precipitation",
                out_png=os.path.join(OUTDIR, f"precip_annual_{year}.png"),
                cmap="YlGnBu",
                contour_levels=PREC_CONTOURS_ALL,
                main_levels=PREC_CONTOURS_MAIN,
                cbar_label="mm/year",
            )

            plot_contour_map(
                R,
                title=f"{year} Annual Runoff",
                out_png=os.path.join(OUTDIR, f"runoff_annual_{year}.png"),
                cmap="viridis",
                contour_levels=RUNOFF_CONTOURS,
                main_levels=None,
                cbar_label="mm/year",
            )

            plot_contour_map(
                E,
                title=f"{year} Annual Evaporation",
                out_png=os.path.join(OUTDIR, f"evap_annual_{year}.png"),
                cmap="YlOrBr",
                contour_levels=EVAP_CONTOURS,
                main_levels=None,
                cbar_label="mm/year",
            )

    # 多年平均年总量
    P_mean = xr.concat(P_list, dim="year").mean(dim="year", skipna=True)
    R_mean = xr.concat(R_list, dim="year").mean(dim="year", skipna=True)
    E_mean = xr.concat(E_list, dim="year").mean(dim="year", skipna=True)
    RES_mean = xr.concat(RES_list, dim="year").mean(dim="year", skipna=True)

    # 多年联合有效掩膜：只保留至少某一年有有效值的区域
    joint_mask = xr.concat(MASK_list, dim="year").any(dim="year")
    P_mean = P_mean.where(joint_mask)
    R_mean = R_mean.where(joint_mask)
    E_mean = E_mean.where(joint_mask)
    RES_mean = RES_mean.where(joint_mask)

    # 多年平均图
    plot_contour_map(
        P_mean,
        title="Multi-year Mean Annual Precipitation",
        out_png=os.path.join(OUTDIR, "precip_annual_mean.png"),
        cmap="YlGnBu",
        contour_levels=PREC_CONTOURS_ALL,
        main_levels=PREC_CONTOURS_MAIN,
        cbar_label="mm/year",
    )

    plot_contour_map(
        R_mean,
        title="Multi-year Mean Annual Runoff",
        out_png=os.path.join(OUTDIR, "runoff_annual_mean.png"),
        cmap="viridis",
        contour_levels=RUNOFF_CONTOURS,
        main_levels=None,
        cbar_label="mm/year",
    )

    plot_contour_map(
        E_mean,
        title="Multi-year Mean Annual Evaporation",
        out_png=os.path.join(OUTDIR, "evap_annual_mean.png"),
        cmap="YlOrBr",
        contour_levels=EVAP_CONTOURS,
        main_levels=None,
        cbar_label="mm/year",
    )

    plot_contour_map(
        RES_mean,
        title="Multi-year Mean Annual Residual (P - R - E)",
        out_png=os.path.join(OUTDIR, "residual_annual_mean.png"),
        contour_levels=RESID_CONTOURS,
        main_levels=None,
        cbar_label="mm/year",
        center_zero=True,
    )

    # 导出 nc
    out_nc = xr.Dataset(
        data_vars={
            "precip_annual_mean": P_mean.astype(np.float32),
            "runoff_annual_mean": R_mean.astype(np.float32),
            "evap_annual_mean": E_mean.astype(np.float32),
            "residual_annual_mean": RES_mean.astype(np.float32),
            "valid_mask": joint_mask.astype(np.int32),
        },
        coords={
            "lat": P_mean["lat"],
            "lon": P_mean["lon"],
        },
        attrs={
            "description": "Multi-year mean annual VIC diagnostics from hydrology.*.nc, excluding zero-value outside cells"
        }
    )
    out_nc_path = os.path.join(OUTDIR, "vic_annual_mean_diagnostics.nc")
    out_nc.to_netcdf(out_nc_path, engine="netcdf4")

    # 导出表格
    df = pd.DataFrame(yearly_records).sort_values("year")
    csv_path = os.path.join(OUTDIR, "yearly_water_balance_summary.csv")
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")

    # 文本摘要
    p_multi = float(P_mean.mean(skipna=True).item())
    r_multi = float(R_mean.mean(skipna=True).item())
    e_multi = float(E_mean.mean(skipna=True).item())
    res_multi = float(RES_mean.mean(skipna=True).item())
    runoff_coef = r_multi / p_multi if abs(p_multi) > 1e-12 else np.nan

    txt_path = os.path.join(OUTDIR, "summary.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("VIC annual diagnostic summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Working directory: {os.getcwd()}\n")
        f.write(f"Files found      : {len(files)}\n")
        f.write(f"Years            : {df['year'].min()} - {df['year'].max()}\n")
        f.write(f"Valid cells      : {int(joint_mask.sum().item())}\n")
        f.write("\n")
        f.write("Domain mean multi-year annual values (excluding zero-value outside cells):\n")
        f.write(f"  Precipitation : {p_multi:.3f} mm/year\n")
        f.write(f"  Runoff        : {r_multi:.3f} mm/year\n")
        f.write(f"  Evaporation   : {e_multi:.3f} mm/year\n")
        f.write(f"  Residual      : {res_multi:.3f} mm/year\n")
        f.write(f"  Runoff coeff  : {runoff_coef:.3f}\n")

    print("\n完成，输出目录：")
    print(OUTDIR)
    print("\n修正后重点图件：")
    print("  precip_annual_mean.png")
    print("  runoff_annual_mean.png")
    print("  evap_annual_mean.png")
    print("  residual_annual_mean.png")
    print("  yearly_water_balance_summary.csv")
    print("  summary.txt")


if __name__ == "__main__":
    main()