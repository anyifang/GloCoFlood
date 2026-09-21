#!/usr/bin/env python
"""Check a SFINCS sfincs_map.nc file and plot maximum-inundation diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from pyproj import Transformer


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "sfincs_map_diagnostics"


def scalar(value) -> float:
    return float(np.asarray(value).reshape(-1)[0])


def quantiles(values: np.ndarray) -> dict[str, float]:
    values = values[np.isfinite(values)]
    levels = (0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0)
    return {f"q{int(q * 100):02d}": float(np.quantile(values, q)) for q in levels}


def inspect_map(map_path: Path, output_dir: Path) -> tuple[dict, Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    with xr.open_dataset(map_path) as ds:
        required = {"msk", "zb", "h", "zs", "hmax", "zsmax", "status"}
        missing = sorted(required.difference(ds.variables))
        if missing:
            raise ValueError(f"Missing required variables: {', '.join(missing)}")

        msk = ds["msk"].values
        zb = ds["zb"].values
        hmax = ds["hmax"].isel(timemax=0).values
        zsmax = ds["zsmax"].isel(timemax=0).values
        h_initial = ds["h"].isel(time=0).values
        h_final = ds["h"].isel(time=-1).values
        x = ds["x"].values
        y = ds["y"].values

        active = np.isfinite(zb) & (msk > 0)
        wet = active & np.isfinite(hmax)
        # A simple topographic land mask. River-channel cells below datum are
        # deliberately excluded from the land-inundation area.
        land_wet = wet & (zb >= 0.0)
        wet_count = int(np.count_nonzero(wet))
        if wet_count == 0:
            raise ValueError("hmax contains no finite active cells")

        input_attrs = ds["inp"].attrs if "inp" in ds.variables else ds.attrs
        dx = abs(float(input_attrs.get("dx", 1.0)))
        dy = abs(float(input_attrs.get("dy", 1.0)))
        cell_area_km2 = dx * dy / 1.0e6
        thresholds = (0.05, 0.10, 0.30, 0.50, 1.0, 2.0, 5.0)
        area_all = {
            str(depth): float(np.count_nonzero(wet & (hmax >= depth)) * cell_area_km2)
            for depth in thresholds
        }
        area_land = {
            str(depth): float(
                np.count_nonzero(land_wet & (hmax >= depth)) * cell_area_km2
            )
            for depth in thresholds
        }

        residual_mask = wet & np.isfinite(zsmax)
        residual = hmax[residual_mask] - (
            zsmax[residual_mask] - zb[residual_mask]
        )

        maximum_index = np.unravel_index(
            np.nanargmax(np.where(wet, hmax, np.nan)), hmax.shape
        )
        maximum_x = float(x[maximum_index])
        maximum_y = float(y[maximum_index])
        if "crs" in ds.variables:
            crs_attrs = ds["crs"].attrs
            crs = crs_attrs.get(
                "epsg_code",
                crs_attrs.get("EPSG", f"EPSG:{int(scalar(ds['crs'].values))}"),
            )
        else:
            crs = input_attrs.get("epsg", "EPSG:4326")
        crs = str(crs)
        if crs.isdigit():
            crs = f"EPSG:{crs}"
        to_geographic = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        lon, lat = to_geographic.transform(maximum_x, maximum_y)

        time_values = np.asarray(ds["time"].values).astype("datetime64[s]")
        timemax_values = np.asarray(ds["timemax"].values).astype("datetime64[s]")
        status = int(scalar(ds["status"].values))
        report = {
            "file": str(map_path.resolve()),
            "file_size_bytes": map_path.stat().st_size,
            "sfincs_build_revision": ds.attrs.get(
                "Build-Revision", ds.attrs.get("Build_revision")
            ),
            "status": status,
            "status_interpretation": (
                "normal completion" if status == 0 else "non-zero SFINCS status"
            ),
            "grid": {
                "shape": list(zb.shape),
                "crs": str(crs),
                "dx_m": dx,
                "dy_m": dy,
                "cell_area_km2": cell_area_km2,
                "active_cells": int(np.count_nonzero(active)),
                "wet_hmax_cells": wet_count,
                "land_wet_cells_zb_ge_0": int(np.count_nonzero(land_wet)),
            },
            "time": {
                "records": [str(t) for t in time_values],
                "maximum_records": [str(t) for t in timemax_values],
                "input_tstart": input_attrs.get("tstart"),
                "input_tstop": input_attrs.get("tstop"),
                "input_dtout_seconds": float(input_attrs.get("dtout", np.nan)),
                "input_dtmaxout_seconds": float(
                    input_attrs.get("dtmaxout", np.nan)
                ),
            },
            "hmax_m_quantiles_all_wet": quantiles(hmax[wet]),
            "hmax_m_quantiles_land_zb_ge_0": quantiles(hmax[land_wet]),
            "inundated_area_km2_all_active": area_all,
            "inundated_area_km2_land_zb_ge_0": area_land,
            "maximum_hmax": {
                "depth_m": float(hmax[maximum_index]),
                "zsmax_m": float(zsmax[maximum_index]),
                "bed_elevation_m": float(zb[maximum_index]),
                "x_m": maximum_x,
                "y_m": maximum_y,
                "longitude": float(lon),
                "latitude": float(lat),
            },
            "consistency_hmax_equals_zsmax_minus_zb": {
                "mean_absolute_residual_m": float(np.mean(np.abs(residual))),
                "p99_absolute_residual_m": float(
                    np.quantile(np.abs(residual), 0.99)
                ),
                "maximum_absolute_residual_m": float(np.max(np.abs(residual))),
            },
            "instantaneous_h": {
                "initial_quantiles_active_m": quantiles(h_initial[active]),
                "final_quantiles_active_m": quantiles(h_final[active]),
            },
        }

        plt.rcParams.update({"font.size": 9})
        model_projection = ccrs.epsg(int(crs.split(":")[-1]))
        geographic_projection = ccrs.PlateCarree()
        model_extent = [
            float(np.nanmin(x) - dx / 2.0),
            float(np.nanmax(x) + dx / 2.0),
            float(np.nanmin(y) - dy / 2.0),
            float(np.nanmax(y) + dy / 2.0),
        ]
        corner_x = [
            model_extent[0],
            model_extent[1],
            model_extent[0],
            model_extent[1],
        ]
        corner_y = [
            model_extent[2],
            model_extent[2],
            model_extent[3],
            model_extent[3],
        ]
        corner_lon, corner_lat = to_geographic.transform(corner_x, corner_y)
        geographic_extent = [
            float(np.min(corner_lon)),
            float(np.max(corner_lon)),
            float(np.min(corner_lat)),
            float(np.max(corner_lat)),
        ]

        sand_cmap = LinearSegmentedColormap.from_list(
            "sfincs_sand", ["#f3edda", "#d8c79c", "#af9867"]
        )
        terrain_raster = np.ma.masked_where(
            ~(active & (zb >= 0.0)), np.clip(zb, 0.0, 60.0)
        )
        minimum_mapped_depth = 0.05
        inundation_raster = np.ma.masked_where(
            ~(land_wet & (hmax >= minimum_mapped_depth)), hmax
        )
        coastal_zs_raster = np.ma.masked_where(
            ~(wet & (zb < 5.0)), zsmax
        )
        land_depths = hmax[land_wet]

        depth_bounds = (0.05, 0.10, 0.30, 0.50, 1.0, 2.0, 3.5)
        depth_cmap = ListedColormap(
            [
                "#d7ebf7",
                "#b7d8ed",
                "#8fc2df",
                "#5aa6d1",
                "#287ab8",
                "#08519c",
                "#08306b",
            ],
            name="sfincs_inundation_blues",
        )
        depth_cmap.set_over("#08306b")
        depth_norm = BoundaryNorm(
            depth_bounds, depth_cmap.N, clip=False, extend="max"
        )

        def style_geographic_axis(axis) -> None:
            axis.set_extent(geographic_extent, crs=geographic_projection)
            axis.add_feature(
                cfeature.OCEAN.with_scale("10m"),
                facecolor="#fbfbf8",
                edgecolor="none",
                zorder=0,
            )
            axis.add_feature(
                cfeature.LAND.with_scale("10m"),
                facecolor="#d8c79c",
                edgecolor="none",
                zorder=0,
            )
            axis.coastlines(
                resolution="10m", color="#7e7258", linewidth=0.55, zorder=4
            )
            gridlines = axis.gridlines(
                crs=geographic_projection,
                draw_labels=True,
                linewidth=0.55,
                color="#777777",
                alpha=0.55,
                linestyle="--",
                zorder=5,
            )
            gridlines.top_labels = False
            gridlines.right_labels = False
            gridlines.xformatter = LongitudeFormatter()
            gridlines.yformatter = LatitudeFormatter()
            gridlines.xlabel_style = {"size": 8}
            gridlines.ylabel_style = {"size": 8}

        def draw_terrain(axis) -> None:
            axis.imshow(
                terrain_raster,
                origin="lower",
                extent=model_extent,
                transform=model_projection,
                interpolation="nearest",
                cmap=sand_cmap,
                vmin=0.0,
                vmax=60.0,
                alpha=0.72,
                zorder=1,
            )

        figure = plt.figure(figsize=(13, 10), constrained_layout=True)
        axes = np.empty((2, 2), dtype=object)
        axes[0, 0] = figure.add_subplot(
            2, 2, 1, projection=geographic_projection
        )
        axes[0, 1] = figure.add_subplot(
            2, 2, 2, projection=geographic_projection
        )
        axes[1, 0] = figure.add_subplot(2, 2, 3)
        axes[1, 1] = figure.add_subplot(2, 2, 4)

        style_geographic_axis(axes[0, 0])
        draw_terrain(axes[0, 0])
        image = axes[0, 0].imshow(
            inundation_raster,
            origin="lower",
            extent=model_extent,
            transform=model_projection,
            interpolation="nearest",
            cmap=depth_cmap,
            norm=depth_norm,
            alpha=0.92,
            zorder=3,
        )
        axes[0, 0].set_title(
            "(a) Maximum land inundation depth (hmax ≥ 0.05 m)"
        )
        figure.colorbar(
            image,
            ax=axes[0, 0],
            label="Maximum depth (m)",
            ticks=depth_bounds,
            extend="max",
            shrink=0.9,
        )

        coastal = wet & (zb < 5.0)
        coastal_zs = zsmax[coastal]
        if coastal_zs.size:
            vmin, vmax = np.nanquantile(coastal_zs, [0.01, 0.99])
            style_geographic_axis(axes[0, 1])
            draw_terrain(axes[0, 1])
            image = axes[0, 1].imshow(
                coastal_zs_raster,
                origin="lower",
                extent=model_extent,
                transform=model_projection,
                interpolation="nearest",
                cmap="coolwarm",
                vmin=vmin,
                vmax=vmax,
                alpha=0.92,
                zorder=3,
            )
            figure.colorbar(
                image,
                ax=axes[0, 1],
                label="Maximum water level (m datum)",
                shrink=0.9,
            )
        axes[0, 1].set_title("(b) zsmax where bed elevation < 5 m")

        axes[1, 0].hist(
            land_depths,
            bins=np.geomspace(0.01, max(0.011, float(np.nanmax(land_depths))), 55),
            color="#277da1",
            edgecolor="none",
        )
        axes[1, 0].set_xscale("log")
        axes[1, 0].set_yscale("log")
        axes[1, 0].set_xlabel("Maximum depth on cells with zb >= 0 (m)")
        axes[1, 0].set_ylabel("Cell count")
        axes[1, 0].set_title("(c) Maximum-depth distribution")
        axes[1, 0].grid(alpha=0.25, which="both")

        labels = [f"{depth:g}" for depth in thresholds]
        all_values = [area_all[str(depth)] for depth in thresholds]
        land_values = [area_land[str(depth)] for depth in thresholds]
        positions = np.arange(len(thresholds))
        width = 0.38
        axes[1, 1].bar(
            positions - width / 2,
            all_values,
            width,
            color="#4d908e",
            label="All active cells",
        )
        axes[1, 1].bar(
            positions + width / 2,
            land_values,
            width,
            color="#f9844a",
            label="Land proxy: zb >= 0",
        )
        axes[1, 1].set_xticks(positions, labels)
        axes[1, 1].set_xlabel("Depth threshold (m)")
        axes[1, 1].set_ylabel("Area above threshold (km²)")
        axes[1, 1].set_title("(d) Thresholded maximum-inundation area")
        axes[1, 1].legend()
        axes[1, 1].grid(alpha=0.25, axis="y")

        figure.suptitle(
            f"SFINCS map-result check: {map_path.parent.name}\n"
            f"{time_values[0]} to {time_values[-1]} | status={status}",
            fontsize=13,
        )

        stem = f"{map_path.parent.name}_sfincs_map_result_check"
        png_path = output_dir / f"{stem}.png"
        pdf_path = output_dir / f"{stem}.pdf"
        json_path = output_dir / f"{stem}.json"
        figure.savefig(png_path, dpi=240)
        figure.savefig(pdf_path)
        plt.close(figure)

        map_stem = f"{map_path.parent.name}_maximum_inundation_map"
        inundation_png_path = output_dir / f"{map_stem}.png"
        inundation_pdf_path = output_dir / f"{map_stem}.pdf"
        map_figure = plt.figure(figsize=(7.2, 7.8), constrained_layout=True)
        map_axis = map_figure.add_subplot(1, 1, 1, projection=geographic_projection)
        style_geographic_axis(map_axis)
        draw_terrain(map_axis)
        image = map_axis.imshow(
            inundation_raster,
            origin="lower",
            extent=model_extent,
            transform=model_projection,
            interpolation="nearest",
            cmap=depth_cmap,
            norm=depth_norm,
            alpha=0.94,
            zorder=3,
        )
        map_axis.text(
            0.018,
            0.978,
            "a",
            transform=map_axis.transAxes,
            ha="left",
            va="top",
            fontsize=18,
            fontweight="bold",
            zorder=6,
        )
        colorbar = map_figure.colorbar(
            image,
            ax=map_axis,
            orientation="horizontal",
            pad=0.055,
            fraction=0.045,
            ticks=depth_bounds,
            extend="max",
        )
        colorbar.set_label(
            "Maximum water depth on land (m; cells with hmax ≥ 0.05 m)"
        )
        map_figure.savefig(inundation_png_path, dpi=300)
        map_figure.savefig(inundation_pdf_path, dpi=300)
        plt.close(map_figure)

        json_path.write_text(
            json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
        )
    return report, png_path, inundation_png_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "map_file",
        type=Path,
        help="Path to sfincs_map.nc",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=f"Diagnostic output folder (default: {DEFAULT_OUTPUT_DIR})",
    )
    args = parser.parse_args()
    map_file = args.map_file.resolve()
    if not map_file.is_file():
        raise FileNotFoundError(map_file)
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else DEFAULT_OUTPUT_DIR
    )
    report, png_path, inundation_png_path = inspect_map(map_file, output_dir)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print(f"\nDiagnostic plot: {png_path}")
    print(f"Maximum-inundation map: {inundation_png_path}")


if __name__ == "__main__":
    main()
