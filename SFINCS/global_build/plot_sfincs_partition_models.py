#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""Plot existing SFINCS partition models.

This script reads written SFINCS native files and creates one PNG per model.
The left panel shows the broader coastal-shelf context with the matched
CaMa-Flood river network; the right panel shows the SFINCS computational
domain and forcing boundaries.
  python .\plot_sfincs_partition_models.py --models-dir .\global_sfincs_partition_models_cama_15min
"""

from __future__ import annotations

import argparse
import math
import os
import re
import traceback
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
from pyproj import CRS, Transformer
from shapely.geometry import GeometryCollection, LineString, MultiPolygon, Polygon, box
from shapely.ops import transform as shapely_transform


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
HYDROMT_EXAMPLES = Path(
    os.environ.get(
        "HYDROMT_SFINCS_EXAMPLES",
        REPO_ROOT / "external" / "HydroMT-SFINCS" / "examples",
    )
).expanduser()
LAND_POLY_SHP = HYDROMT_EXAMPLES / "data" / "coastaline" / "land_polygons.shp"
SEA_FACE_COLOR = "#d9eef7"
DEFAULT_PARTITION_GEOJSON = (
    SCRIPT_DIR
    / "global_tc_coastal_model_partitions"
    / "global_tc_coastal_model_domain_boundaries.geojson"
)
DEFAULT_DOMAIN_CSV = (
    SCRIPT_DIR
    / "global_tc_coastal_model_partitions"
    / "global_tc_coastal_model_domains.csv"
)
DEFAULT_MEMBERSHIP_CSV = (
    REPO_ROOT
    / "ADCIRC"
    / "global_build"
    / "catalog"
    / "global_tc_adcirc_block_members.csv"
)
DEFAULT_MODELS_DIR = SCRIPT_DIR / "global_sfincs_partition_models"
FONT_FAMILY = "Arial"
BASE_FONT_SIZE = 14
AXIS_LABEL_FONT_SIZE = 14
TICK_LABEL_FONT_SIZE = 14
LEGEND_FONT_SIZE = 14
LEGEND_TITLE_FONT_SIZE = 14
DOMAIN_LABEL_FONT_SIZE = 14
PANEL_LABEL_FONT_SIZE = 20
COLORBAR_LABEL_FONT_SIZE = 14
COLORBAR_TICK_FONT_SIZE = 14

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": [FONT_FAMILY],
        "font.size": BASE_FONT_SIZE,
        "axes.labelsize": AXIS_LABEL_FONT_SIZE,
        "axes.titlesize": AXIS_LABEL_FONT_SIZE,
        "xtick.labelsize": TICK_LABEL_FONT_SIZE,
        "ytick.labelsize": TICK_LABEL_FONT_SIZE,
        "legend.fontsize": LEGEND_FONT_SIZE,
        "legend.title_fontsize": LEGEND_TITLE_FONT_SIZE,
        "figure.titlesize": AXIS_LABEL_FONT_SIZE,
        "axes.unicode_minus": False,
    }
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", nargs="*", default=None, help="Partition numbers to plot, e.g. No_001 No_067.")
    parser.add_argument("--limit", type=int, default=None, help="Only plot the first N selected models.")
    parser.add_argument("--start-rank", type=int, default=None, help="Start at P2 plot_rank >= this value.")
    parser.add_argument("--end-rank", type=int, default=None, help="Stop at P2 plot_rank <= this value.")
    parser.add_argument("--models-dir", type=Path, default=DEFAULT_MODELS_DIR, help="Directory containing model folders.")
    parser.add_argument("--plot-dir", type=Path, default=None, help="Directory for generated PNG files.")
    parser.add_argument("--partition-geojson", type=Path, default=DEFAULT_PARTITION_GEOJSON)
    parser.add_argument("--domain-csv", type=Path, default=DEFAULT_DOMAIN_CSV)
    parser.add_argument(
        "--land-polygons",
        type=Path,
        default=LAND_POLY_SHP,
        help="Land polygon shapefile used for the regional context panel.",
    )
    parser.add_argument(
        "--membership-csv",
        type=Path,
        default=DEFAULT_MEMBERSHIP_CSV,
        help="GTC-to-ADCIRC membership table used in plot labels.",
    )
    parser.add_argument("--stop-on-error", action="store_true", help="Stop after the first failed plot.")
    return parser.parse_args()


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


def partition_model_folder_name(model_id: str, plot_rank=None) -> str:
    rank = coerce_rank(plot_rank)
    if rank is None:
        return str(model_id)
    return f"No_{rank:03d}"


def partition_plot_stem(
    model_id: str,
    plot_rank=None,
    adc_block_id: str | None = None,
) -> str:
    rank = coerce_rank(plot_rank)
    parts = [f"No.{rank}" if rank is not None else "No.unknown", str(model_id)]
    if adc_block_id:
        parts.append(str(adc_block_id))
    return "_".join(parts)


def partition_legend_title(model_id: str, plot_rank=None) -> str:
    rank = coerce_rank(plot_rank)
    return f"No.{rank}" if rank is not None else str(model_id)


def partition_display_label(model_id: str, plot_rank=None, adc_block_id: str | None = None) -> str:
    rank = coerce_rank(plot_rank)
    parts = [str(model_id)]
    if rank is not None:
        parts.append(f"No.{rank}")
    if adc_block_id:
        parts.append(str(adc_block_id))
    return " | ".join(parts)


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


def find_partition_model_dir(models_dir: Path, model_id: str, plot_rank=None) -> Path:
    models_dir = Path(models_dir)
    candidates: list[Path] = []
    rank = coerce_rank(plot_rank)
    candidates.append(models_dir / partition_model_folder_name(model_id, plot_rank))
    if rank is not None:
        candidates.append(models_dir / f"No_{rank}")
        candidates.append(models_dir / f"No_{rank:03d}_{model_id}")
        candidates.append(models_dir / f"No_{rank}_{model_id}")
        candidates.extend(sorted(models_dir.glob(f"No_{rank:03d}*")))
        candidates.extend(sorted(models_dir.glob(f"No_{rank}*")))
    candidates.append(models_dir / str(model_id))
    candidates.extend(sorted(models_dir.glob(f"No_*_{model_id}")))

    seen: set[Path] = set()
    unique_candidates: list[Path] = []
    for path in candidates:
        if path in seen:
            continue
        seen.add(path)
        unique_candidates.append(path)

    for path in unique_candidates:
        if (path / "sfincs.inp").exists():
            return path
    for path in unique_candidates:
        if path.exists():
            return path
    return unique_candidates[0]


def load_partitions(partition_geojson: Path, domain_csv: Path) -> gpd.GeoDataFrame:
    if not partition_geojson.exists():
        raise FileNotFoundError(partition_geojson)
    gdf = gpd.read_file(partition_geojson).to_crs(4326)
    gdf["geometry"] = gdf.geometry.apply(clean_geometry)
    gdf = gdf[~gdf.geometry.is_empty].copy()
    if domain_csv.exists():
        meta = pd.read_csv(domain_csv)
        if "model_domain_id" in meta.columns:
            gdf = gdf.merge(
                meta.drop(columns=[c for c in ["geometry"] if c in meta.columns]),
                on="model_domain_id",
                how="left",
                suffixes=("", "_domain"),
            )
    gdf = require_saved_plot_rank(gdf, domain_csv if domain_csv.exists() else partition_geojson)
    gdf = gdf.sort_values(["plot_rank", "model_domain_id"], kind="mergesort")
    return gdf.reset_index(drop=True)


def select_model_ids(args: argparse.Namespace) -> list[str]:
    if args.partition_geojson.exists():
        gdf = load_partitions(args.partition_geojson, args.domain_csv)
        selected = gdf
        if args.ids:
            ids = {model_id_from_folder_name(v) for v in args.ids if model_id_from_folder_name(v).startswith("GTC_")}
            ranks = {plot_rank_from_folder_name(v) for v in args.ids}
            ranks.discard(None)
            mask = selected["model_domain_id"].isin(ids)
            if ranks and "plot_rank" in selected.columns:
                mask = mask | pd.to_numeric(selected["plot_rank"], errors="coerce").isin(ranks)
            selected = selected[mask]
        if args.start_rank is not None:
            selected = selected[pd.to_numeric(selected["plot_rank"], errors="coerce") >= args.start_rank]
        if args.end_rank is not None:
            selected = selected[pd.to_numeric(selected["plot_rank"], errors="coerce") <= args.end_rank]
        if args.limit is not None:
            selected = selected.head(args.limit)
        return [str(v) for v in selected["model_domain_id"]]

    if args.ids:
        ids = [model_id_from_folder_name(v) for v in args.ids]
    else:
        ids = sorted(set(
            model_id_from_folder_name(p.name)
            for p in args.models_dir.iterdir()
            if (p / "sfincs.inp").exists()
        ))
    if args.limit is not None:
        ids = ids[: args.limit]
    return ids


def read_sfincs_inp(path: Path) -> dict:
    cfg: dict[str, object] = {}
    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, val = [part.strip() for part in line.split("=", 1)]
        val = val.split("#", 1)[0].strip()
        try:
            if "." in val or "e" in val.lower():
                cfg[key] = float(val)
            else:
                if val and all(ch in "+-0123456789" for ch in val):
                    cfg[key] = int(val)
                else:
                    cfg[key] = val
        except ValueError:
            cfg[key] = val
    return cfg


def sfincs_active_indices_to_world(inp: dict, active_ind: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nmax = int(inp["nmax"])
    x0 = float(inp.get("x0", 0.0))
    y0 = float(inp.get("y0", 0.0))
    dx = float(inp.get("dx", 1.0))
    dy = float(inp.get("dy", dx))
    rot = math.radians(float(inp.get("rotation", 0.0)))

    idx0 = active_ind.astype(np.int64) - 1
    ii = idx0 // nmax
    jj = idx0 % nmax
    xl = (ii.astype(float) + 0.5) * dx
    yl = (jj.astype(float) + 0.5) * dy
    xw = x0 + xl * math.cos(rot) - yl * math.sin(rot)
    yw = y0 + xl * math.sin(rot) + yl * math.cos(rot)
    return xw, yw


def sfincs_grid_edges_to_world(
    inp: dict, col0: int, col1: int, row0: int, row1: int
) -> tuple[np.ndarray, np.ndarray]:
    x0 = float(inp.get("x0", 0.0))
    y0 = float(inp.get("y0", 0.0))
    dx = float(inp.get("dx", 1.0))
    dy = float(inp.get("dy", dx))
    rot = math.radians(float(inp.get("rotation", 0.0)))

    cols = np.arange(col0, col1 + 2, dtype=float)
    rows = np.arange(row0, row1 + 2, dtype=float)
    xl, yl = np.meshgrid(cols * dx, rows * dy)
    xw = x0 + xl * math.cos(rot) - yl * math.sin(rot)
    yw = y0 + xl * math.sin(rot) + yl * math.cos(rot)
    return xw, yw


def read_xy_file(path: Path) -> np.ndarray:
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 2), dtype=float)
    data = np.loadtxt(path, dtype=float)
    if data.size == 0:
        return np.empty((0, 2), dtype=float)
    data = np.atleast_2d(data)
    if data.shape[1] < 2:
        return np.empty((0, 2), dtype=float)
    return data[:, :2]


def read_domain_gdf(model_dir: Path) -> gpd.GeoDataFrame | None:
    for name in ["active_domain.geojson", "region_grid.geojson"]:
        path = model_dir / name
        if path.exists():
            gdf = gpd.read_file(path).to_crs(4326)
            gdf["geometry"] = gdf.geometry.apply(clean_geometry)
            gdf = gdf[~gdf.geometry.is_empty].copy()
            if not gdf.empty:
                return gdf
    return None


def lon_to_domain_frame(lon, bounds):
    minx, _, maxx, _ = bounds
    arr = np.asarray(lon, dtype=float).copy()
    if np.isfinite(minx) and np.isfinite(maxx):
        if maxx > 180.0:
            arr[arr < minx] += 360.0
        elif minx < -180.0:
            arr[arr > maxx] -= 360.0
    return arr


def geometry_to_domain_frame(geom, bounds):
    minx, _, maxx, _ = bounds
    if not (np.isfinite(minx) and np.isfinite(maxx)):
        return geom

    def _shift(x, y, z=None):
        x_arr = np.asarray(x, dtype=float)
        if maxx > 180.0:
            x_out = np.where(x_arr < minx, x_arr + 360.0, x_arr)
        elif minx < -180.0:
            x_out = np.where(x_arr > maxx, x_arr - 360.0, x_arr)
        else:
            x_out = x_arr
        if z is None:
            return x_out, y
        return x_out, y, z

    return shapely_transform(_shift, geom)


def shift_geometry_lon(geom, delta: float):
    return shapely_transform(
        lambda x, y, z=None: (np.asarray(x, dtype=float) + delta, y)
        if z is None
        else (np.asarray(x, dtype=float) + delta, y, z),
        geom,
    )


def format_lon_label(value, _pos=None) -> str:
    value = float(value)
    hemi = "E" if value >= 0 else "W"
    return f"{abs(value):g}°{hemi}"


def format_lat_label(value, _pos=None) -> str:
    value = float(value)
    hemi = "N" if value >= 0 else "S"
    return f"{abs(value):g}°{hemi}"


def apply_lonlat_formatters(ax):
    ax.xaxis.set_major_formatter(FuncFormatter(format_lon_label))
    ax.yaxis.set_major_formatter(FuncFormatter(format_lat_label))
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        top=True,
        right=True,
        labelsize=TICK_LABEL_FONT_SIZE,
    )


def domain_label_text(domain_gdf: gpd.GeoDataFrame | None, model_id: str, plot_rank=None) -> str:
    for value in [plot_rank]:
        try:
            if pd.notna(value):
                return str(int(round(float(value))))
        except Exception:
            pass
    if domain_gdf is not None and not domain_gdf.empty:
        for col in ["plot_rank", "position_connected_rank", "rank_domain"]:
            if col not in domain_gdf.columns:
                continue
            value = domain_gdf[col].iloc[0]
            try:
                if pd.notna(value):
                    return str(int(round(float(value))))
            except Exception:
                pass
    try:
        return str(int(str(model_id).split("_")[-1]))
    except Exception:
        return str(model_id)


def annotate_domain_label(ax, domain_gdf: gpd.GeoDataFrame | None, model_id: str, ctx=None, plot_rank=None):
    if domain_gdf is None or domain_gdf.empty:
        return
    geom = domain_gdf.geometry.iloc[0]
    if geom is None or geom.is_empty:
        return
    minx, miny, maxx, maxy = [float(v) for v in geom.bounds]
    cx = 0.5 * (minx + maxx)
    cy = 0.5 * (miny + maxy)
    if ctx is not None:
        west, south, east, north = [float(v) for v in ctx]
        distances = {
            "west": abs(minx - west),
            "east": abs(east - maxx),
            "south": abs(miny - south),
            "north": abs(north - maxy),
        }
        side = max(distances, key=distances.get)
        offset = max(0.06 * max(maxx - minx, maxy - miny, 1.0), 0.18)
        if side == "west":
            x, y, ha, va = minx - offset, cy, "right", "center"
        elif side == "east":
            x, y, ha, va = maxx + offset, cy, "left", "center"
        elif side == "south":
            x, y, ha, va = cx, miny - offset, "center", "top"
        else:
            x, y, ha, va = cx, maxy + offset, "center", "bottom"
    else:
        point = geom.representative_point()
        x, y, ha, va = float(point.x), float(point.y), "center", "center"
    ax.text(
        float(x),
        float(y),
        domain_label_text(domain_gdf, model_id, plot_rank),
        ha=ha,
        va=va,
        fontsize=DOMAIN_LABEL_FONT_SIZE,
        fontfamily=FONT_FAMILY,
        fontweight="bold",
        color="#7f0000",
        zorder=11,
    )


def annotate_panel_label(ax, label: str):
    ax.text(
        0.02,
        0.98,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=PANEL_LABEL_FONT_SIZE,
        fontfamily=FONT_FAMILY,
        fontweight="bold",
        color="black",
        bbox={"facecolor": "none", "edgecolor": "none", "alpha": 0.75, "pad": 2.0},
        zorder=30,
    )


def read_land_for_context(ctx: tuple[float, float, float, float]) -> gpd.GeoDataFrame:
    if not LAND_POLY_SHP.exists():
        return gpd.GeoDataFrame(geometry=[], crs=4326)

    west, south, east, north = [float(v) for v in ctx]
    parts = []

    def _read_one(bbox, shift=0.0):
        gdf = gpd.read_file(LAND_POLY_SHP, bbox=tuple(bbox)).to_crs(4326)
        gdf = gdf[gdf.geometry.type.isin(["Polygon", "MultiPolygon"])].copy()
        if gdf.empty:
            return None
        if shift:
            gdf["geometry"] = gdf.geometry.apply(lambda geom: shift_geometry_lon(geom, shift))
        return gdf

    if east > 180.0:
        requests = [
            ((max(-180.0, west), south, 180.0, north), 0.0),
            ((-180.0, south, east - 360.0, north), 360.0),
        ]
    elif west < -180.0:
        requests = [
            ((west + 360.0, south, 180.0, north), -360.0),
            ((-180.0, south, min(180.0, east), north), 0.0),
        ]
    else:
        requests = [((west, south, east, north), 0.0)]

    for bbox, shift in requests:
        try:
            part = _read_one(bbox, shift)
            if part is not None:
                parts.append(part)
        except Exception:
            pass

    if not parts:
        return gpd.GeoDataFrame(geometry=[], crs=4326)
    out = pd.concat(parts, ignore_index=True)
    return gpd.GeoDataFrame(out, geometry="geometry", crs=parts[0].crs)


def polygon_exterior_lines(geom) -> list[LineString]:
    geom = clean_geometry(geom)
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [LineString(geom.exterior.coords)] if geom.exterior is not None else []
    if isinstance(geom, MultiPolygon):
        lines = []
        for part in geom.geoms:
            lines.extend(polygon_exterior_lines(part))
        return lines
    if isinstance(geom, GeometryCollection):
        lines = []
        for part in geom.geoms:
            lines.extend(polygon_exterior_lines(part))
        return lines
    return []


def read_coastline_for_context(ctx: tuple[float, float, float, float]) -> gpd.GeoDataFrame:
    land = read_land_for_context(ctx)
    if land.empty:
        return gpd.GeoDataFrame(geometry=[], crs=4326)

    coast_geoms = []
    for geom in land.geometry:
        coast_geoms.extend(polygon_exterior_lines(geom))
    if not coast_geoms:
        return gpd.GeoDataFrame(geometry=[], crs=land.crs)
    return gpd.GeoDataFrame(geometry=coast_geoms, crs=land.crs)


def split_elevation_cmap() -> LinearSegmentedColormap:
    negative = plt.cm.Blues(np.linspace(0.95, 0.28, 128))
    positive = plt.cm.YlGn(np.linspace(0.22, 0.95, 128))
    cmap = LinearSegmentedColormap.from_list("sfincs_split_topobathy", np.vstack([negative, positive]))
    cmap.set_bad((1.0, 1.0, 1.0, 0.0))
    return cmap


def read_cama_match_gdf(model_dir: Path, domain_gdf: gpd.GeoDataFrame | None) -> gpd.GeoDataFrame:
    path = model_dir / "cama_upstream_rivers.geojson"
    if not path.exists() or path.stat().st_size == 0:
        return gpd.GeoDataFrame(geometry=[], crs=4326)

    cama = gpd.read_file(path).to_crs(4326)
    if domain_gdf is not None and not domain_gdf.empty:
        domain_bounds = tuple(float(v) for v in domain_gdf.total_bounds)
        cama = cama.copy()
        cama["geometry"] = cama.geometry.apply(lambda geom: geometry_to_domain_frame(geom, domain_bounds))
    return cama


def plot_cama_match_layers(ax, cama: gpd.GeoDataFrame, *, compact: bool):
    if cama.empty:
        return

    role = cama["geometry_role"] if "geometry_role" in cama.columns else pd.Series([""] * len(cama), index=cama.index)
    geom_type = cama.geometry.geom_type
    line_geom = geom_type.isin(["LineString", "MultiLineString"])
    point_geom = geom_type == "Point"
    upstream = cama[(role == "camaflood_upstream_river_line") & line_geom]
    links = cama[role.isin(["section_to_camaflood_link", "sfincs_point_to_camaflood_link"]) & line_geom]
    sfincs_lines = cama[(role == "sfincs_inflow_boundary_line") & line_geom]
    src_pts = cama[(role == "sfincs_inflow_point") & point_geom]

    if not upstream.empty:
        if "uparea_km2" in upstream.columns:
            widths = np.asarray(upstream["uparea_km2"], dtype=float)
            widths = np.where(np.isfinite(widths) & (widths > 0), widths, np.nan)
            if np.any(np.isfinite(widths)):
                log_width = np.log10(np.maximum(widths, 1.0))
                log_min = float(np.nanmin(log_width))
                log_max = float(np.nanmax(log_width))
                denom = max(log_max - log_min, 1e-6)
                line_widths = 0.25 + 1.75 * (log_width - log_min) / denom
            else:
                line_widths = np.full(len(upstream), 0.8)
        else:
            line_widths = np.full(len(upstream), 0.8)
        for geom, lw in zip(upstream.geometry, line_widths):
            if geom is not None and not geom.is_empty:
                gpd.GeoSeries([geom], crs=upstream.crs).plot(
                    ax=ax,
                    color="#2166ac",
                    linewidth=float(lw),
                    alpha=0.78,
                    zorder=3,
                )
    if not links.empty:
        links.plot(ax=ax, color="#969696", linewidth=0.55, alpha=0.65, zorder=4, label="_nolegend_")
    if not sfincs_lines.empty:
        sfincs_lines.plot(ax=ax, color="#7b3294", linewidth=1.3, zorder=5, label="SFINCS inflow section")
    if not src_pts.empty:
        ax.scatter(
            src_pts.geometry.x,
            src_pts.geometry.y,
            s=36 if compact else 46,
            c="#6a1b9a",
            marker="^",
            edgecolors="none",
            linewidths=0,
            zorder=7,
            label="SFINCS river-discharge boundary",
        )


def plot_land_context(
    ax,
    domain_gdf: gpd.GeoDataFrame | None,
    cama: gpd.GeoDataFrame,
    model_id: str,
    plot_rank=None,
    adc_block_id: str | None = None,
):
    ax.set_facecolor(SEA_FACE_COLOR)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, color="0.88", linewidth=0.5)

    if domain_gdf is None or domain_gdf.empty:
        ax.text(
            0.5,
            0.5,
            "No active_domain.geojson",
            transform=ax.transAxes,
            ha="center",
            fontsize=BASE_FONT_SIZE,
            fontfamily=FONT_FAMILY,
        )
        return None

    minx, miny, maxx, maxy = [float(v) for v in domain_gdf.total_bounds]
    if not cama.empty:
        cminx, cminy, cmaxx, cmaxy = [float(v) for v in cama.total_bounds]
        minx = min(minx, cminx)
        miny = min(miny, cminy)
        maxx = max(maxx, cmaxx)
        maxy = max(maxy, cmaxy)

    span = max(float(maxx - minx), float(maxy - miny), 1.0)
    pad = max(2.0, span * 0.35)
    ctx = (
        float(minx) - pad,
        max(-60.0, float(miny) - pad),
        float(maxx) + pad,
        min(85.0, float(maxy) + pad),
    )

    land = read_land_for_context(ctx)
    if not land.empty:
        land.plot(ax=ax, facecolor="#e5e0d4", edgecolor="#9c968a", linewidth=0.35, zorder=1)

    plot_cama_match_layers(ax, cama, compact=True)
    domain_gdf.plot(ax=ax, facecolor="#d7191c", edgecolor="none", alpha=0.14, zorder=7)
    domain_gdf.boundary.plot(ax=ax, color="#d7191c", linewidth=1.4, zorder=8, label="SFINCS domain")
    ax.set_xlim(ctx[0], ctx[2])
    ax.set_ylim(ctx[1], ctx[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_anchor("E")
    apply_lonlat_formatters(ax)
    custom_handles = [
        Line2D([0], [0], color="#2166ac", linewidth=1.2, label="CaMa-flood river network"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="#6a1b9a", markeredgecolor="#6a1b9a", markersize=7, label="SFINCS river discharge boundary"),
        Line2D([0], [0], color="#fdae61", linewidth=1.8, label="SFINCS water level boundary"),
        Line2D([0], [0], color="#d7191c", linewidth=1.4, label="SFINCS domain"),
    ]
    legend = ax.legend(
        handles=custom_handles,
        loc="lower left",
        fontsize=LEGEND_FONT_SIZE,
        frameon=True,
        title=partition_legend_title(model_id, plot_rank),
        title_fontsize=LEGEND_TITLE_FONT_SIZE,
    )
    legend.get_title().set_fontweight("bold")
    legend.get_title().set_fontfamily(FONT_FAMILY)
    for text in legend.get_texts():
        text.set_fontfamily(FONT_FAMILY)
    return ctx


def equal_aspect_axis_width(
    ctx: tuple[float, float, float, float],
    fig,
    height_frac: float,
    *,
    min_width: float,
    max_width: float,
) -> float:
    west, south, east, north = [float(v) for v in ctx]
    dx = max(abs(east - west), 1e-6)
    dy = max(abs(north - south), 1e-6)
    fig_w, fig_h = fig.get_size_inches()
    width = height_frac * (dx / dy) * (fig_h / fig_w)
    return float(np.clip(width, min_width, max_width))


def apply_compact_equal_aspect_layout(fig, ax_left, ax_right, cax, left_ctx, right_ctx):
    bottom = 0.11
    height = 0.82
    left = 0.055
    gap = 0.070
    cbar_gap = 0.020
    cbar_width = 0.018

    left_width = equal_aspect_axis_width(
        left_ctx,
        fig,
        height,
        min_width=0.34,
        max_width=0.47,
    )
    right_width = equal_aspect_axis_width(
        right_ctx,
        fig,
        height,
        min_width=0.17,
        max_width=0.34,
    )

    max_total_right = 0.965
    total = left + left_width + gap + right_width + cbar_gap + cbar_width
    if total > max_total_right:
        overflow = total - max_total_right
        shrink_left = min(overflow * 0.65, max(0.0, left_width - 0.30))
        left_width -= shrink_left
        overflow -= shrink_left
        shrink_right = min(overflow, max(0.0, right_width - 0.15))
        right_width -= shrink_right

    right_x = left + left_width + gap
    cax_x = right_x + right_width + cbar_gap
    ax_left.set_position([left, bottom, left_width, height])
    ax_right.set_position([right_x, bottom, right_width, height])
    cax.set_position([cax_x, bottom + 0.05, cbar_width, height - 0.10])
    ax_left.set_anchor("E")
    ax_right.set_anchor("W")


def create_sfincs_domain_forcing_plot(
    model_dir: Path,
    model_id: str,
    plot_dir: Path,
    partition_geom=None,
    plot_rank=None,
    adc_block_id: str | None = None,
) -> Path:
    inp_path = model_dir / "sfincs.inp"
    ind_path = model_dir / "sfincs.ind"
    dep_path = model_dir / "sfincs.dep"
    msk_path = model_dir / "sfincs.msk"
    if not all(p.exists() for p in [inp_path, ind_path, dep_path, msk_path]):
        raise FileNotFoundError(f"Missing one or more SFINCS files in {model_dir}")

    inp = read_sfincs_inp(inp_path)
    epsg = int(inp["epsg"])
    ind = np.fromfile(ind_path, dtype="<i4")
    if ind.size < 2:
        raise RuntimeError(f"Invalid sfincs.ind: {ind_path}")
    n_active = int(ind[0])
    active_ind = ind[1:]
    if active_ind.size != n_active:
        active_ind = active_ind[:n_active]

    dep = np.fromfile(dep_path, dtype="<f4")
    msk = np.fromfile(msk_path, dtype=np.uint8)
    n = min(n_active, dep.size, msk.size, active_ind.size)
    active_ind = active_ind[:n]
    dep = dep[:n]
    msk = msk[:n]

    valid = np.isfinite(dep) & (msk > 0)
    active_ind = active_ind[valid]
    dep = dep[valid]

    xw, yw = sfincs_active_indices_to_world(inp, active_ind)
    transformer = Transformer.from_crs(CRS.from_epsg(epsg), CRS.from_epsg(4326), always_xy=True)
    lon, lat = transformer.transform(xw, yw)

    mmax = int(inp["mmax"])
    nmax = int(inp["nmax"])
    idx0 = active_ind.astype(np.int64) - 1
    cols = idx0 // nmax
    rows = idx0 % nmax
    inside = (cols >= 0) & (cols < mmax) & (rows >= 0) & (rows < nmax)
    if not np.any(inside):
        raise RuntimeError(f"No valid active SFINCS cells for plotting: {model_dir}")
    cols = cols[inside]
    rows = rows[inside]
    dep_for_grid = dep[inside]

    dep_grid = np.full((nmax, mmax), np.nan, dtype=np.float32)
    dep_grid[rows, cols] = dep_for_grid.astype(np.float32, copy=False)
    finite_grid = np.isfinite(dep_grid)
    active_rows = np.where(np.any(finite_grid, axis=1))[0]
    active_cols = np.where(np.any(finite_grid, axis=0))[0]
    if active_rows.size == 0 or active_cols.size == 0:
        raise RuntimeError(f"No finite SFINCS depth values for plotting: {model_dir}")

    pad_cells = 2
    row0 = max(0, int(active_rows.min()) - pad_cells)
    row1 = min(nmax - 1, int(active_rows.max()) + pad_cells)
    col0 = max(0, int(active_cols.min()) - pad_cells)
    col1 = min(mmax - 1, int(active_cols.max()) + pad_cells)
    dep_crop = dep_grid[row0 : row1 + 1, col0 : col1 + 1]

    xedge, yedge = sfincs_grid_edges_to_world(inp, col0, col1, row0, row1)
    lon_edge, lat_edge = transformer.transform(xedge, yedge)

    bnd_xy = read_xy_file(model_dir / "sfincs.bnd")
    if bnd_xy.size:
        bnd_lon, bnd_lat = transformer.transform(bnd_xy[:, 0], bnd_xy[:, 1])
    else:
        bnd_lon = np.array([])
        bnd_lat = np.array([])

    src_xy = read_xy_file(model_dir / "sfincs.src")
    if src_xy.size:
        src_lon, src_lat = transformer.transform(src_xy[:, 0], src_xy[:, 1])
    else:
        src_lon = np.array([])
        src_lat = np.array([])

    domain_gdf = read_domain_gdf(model_dir)
    if domain_gdf is not None and not domain_gdf.empty and plot_rank is not None:
        domain_gdf = domain_gdf.copy()
        domain_gdf["plot_rank"] = plot_rank
    if domain_gdf is not None and not domain_gdf.empty:
        domain_bounds = tuple(float(v) for v in domain_gdf.total_bounds)
        lon = lon_to_domain_frame(lon, domain_bounds)
        lon_edge = lon_to_domain_frame(lon_edge, domain_bounds)
        if bnd_lon.size:
            bnd_lon = lon_to_domain_frame(bnd_lon, domain_bounds)
        if src_lon.size:
            src_lon = lon_to_domain_frame(src_lon, domain_bounds)

    cama = read_cama_match_gdf(model_dir, domain_gdf)

    if partition_geom is not None and domain_gdf is not None and not domain_gdf.empty:
        domain_bounds = tuple(float(v) for v in domain_gdf.total_bounds)
        partition_geom = geometry_to_domain_frame(clean_geometry(partition_geom), domain_bounds)

    plot_dir.mkdir(parents=True, exist_ok=True)
    out_png = plot_dir / f"{partition_plot_stem(model_id, plot_rank, adc_block_id)}.png"
    fig = plt.figure(figsize=(16.5, 6.4))
    ax_context = fig.add_axes([0.055, 0.11, 0.47, 0.82])
    ax = fig.add_axes([0.55, 0.11, 0.30, 0.82])
    cax = fig.add_axes([0.89, 0.16, 0.018, 0.72])
    left_ctx = plot_land_context(
        ax_context,
        domain_gdf,
        cama,
        model_id,
        plot_rank,
        adc_block_id,
    )

    vmax = 15.0
    vmin = -15.0
    cmap = split_elevation_cmap()
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    ax.set_facecolor(SEA_FACE_COLOR)
    sc = ax.pcolormesh(
        lon_edge,
        lat_edge,
        np.clip(dep_crop, vmin, vmax),
        cmap=cmap,
        norm=norm,
        shading="auto",
        rasterized=True,
    )

    finite_lon = lon[np.isfinite(lon)]
    finite_lat = lat[np.isfinite(lat)]
    extent_lons = [finite_lon]
    extent_lats = [finite_lat]
    if bnd_lon.size:
        extent_lons.append(bnd_lon[np.isfinite(bnd_lon)])
        extent_lats.append(bnd_lat[np.isfinite(bnd_lat)])
    if src_lon.size:
        extent_lons.append(src_lon[np.isfinite(src_lon)])
        extent_lats.append(src_lat[np.isfinite(src_lat)])
    if partition_geom is not None and not partition_geom.is_empty:
        pminx, pminy, pmaxx, pmaxy = [float(v) for v in partition_geom.bounds]
        extent_lons.append(np.array([pminx, pmaxx], dtype=float))
        extent_lats.append(np.array([pminy, pmaxy], dtype=float))
    extent_lons = [arr for arr in extent_lons if arr.size]
    extent_lats = [arr for arr in extent_lats if arr.size]

    if extent_lons and extent_lats:
        all_extent_lon = np.concatenate(extent_lons)
        all_extent_lat = np.concatenate(extent_lats)
        right_minx = float(np.nanmin(all_extent_lon))
        right_maxx = float(np.nanmax(all_extent_lon))
        right_miny = float(np.nanmin(all_extent_lat))
        right_maxy = float(np.nanmax(all_extent_lat))
    elif domain_gdf is not None and not domain_gdf.empty:
        right_minx, right_miny, right_maxx, right_maxy = [float(v) for v in domain_gdf.total_bounds]
    else:
        right_minx, right_miny, right_maxx, right_maxy = (-0.5, -0.5, 0.5, 0.5)

    right_span = max(
        float(right_maxx - right_minx) if np.isfinite(right_maxx - right_minx) else 0.0,
        float(right_maxy - right_miny) if np.isfinite(right_maxy - right_miny) else 0.0,
        0.05,
    )
    right_pad = max(0.02, 0.035 * right_span)
    right_ctx = (
        float(right_minx) - right_pad,
        float(right_miny) - right_pad,
        float(right_maxx) + right_pad,
        float(right_maxy) + right_pad,
    )
    if left_ctx is not None:
        apply_compact_equal_aspect_layout(fig, ax_context, ax, cax, left_ctx, right_ctx)
    coastline = read_coastline_for_context(right_ctx)
    right_land = read_land_for_context(right_ctx)
    if not right_land.empty and domain_gdf is not None and not domain_gdf.empty:
        domain_bounds = tuple(float(v) for v in domain_gdf.total_bounds)
        right_land = right_land.copy()
        right_land["geometry"] = right_land.geometry.apply(lambda geom: geometry_to_domain_frame(geom, domain_bounds))
    if not right_land.empty:
        right_land.plot(
            ax=ax,
            facecolor="#e5e0d4",
            edgecolor="none",
            linewidth=0.0,
            zorder=0.6,
        )

    if not coastline.empty and domain_gdf is not None and not domain_gdf.empty:
        domain_bounds = tuple(float(v) for v in domain_gdf.total_bounds)
        coastline = coastline.copy()
        coastline["geometry"] = coastline.geometry.apply(lambda geom: geometry_to_domain_frame(geom, domain_bounds))

    if not coastline.empty:
        coastline.plot(ax=ax, color="#4d4d4d", linewidth=0.75, alpha=0.9, zorder=3, label="Coastline")

    if partition_geom is not None and not partition_geom.is_empty:
        gpd.GeoSeries([partition_geom], crs=4326).boundary.plot(
            ax=ax,
            color="#d7191c",
            linewidth=1.35,
            linestyle="-",
            zorder=5,
            label="Partition boundary",
        )

    if bnd_lon.size:
        b_step = max(1, int(math.ceil(bnd_lon.size / 50000)))
        ax.scatter(
            bnd_lon[::b_step],
            bnd_lat[::b_step],
            s=8,
            c="#fdae61",
            marker=".",
            label="SFINCS coupling boundary",
            linewidths=0,
            zorder=7,
        )
    else:
        ax.scatter([], [], s=8, c="#fdae61", marker=".", label="SFINCS coupling boundary")

    if src_lon.size:
        ax.scatter(
            src_lon,
            src_lat,
            s=70,
            c="#6a1b9a",
            marker="^",
            label="SFINCS river-discharge boundary",
            linewidths=0,
            zorder=8,
        )
    else:
        ax.scatter([], [], s=70, c="#6a1b9a", marker="^", label="SFINCS river-discharge boundary")

    if finite_lon.size and finite_lat.size:
        ax.set_xlim(right_ctx[0], right_ctx[2])
        ax.set_ylim(right_ctx[1], right_ctx[3])

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, color="0.85", linewidth=0.5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_anchor("W")
    apply_lonlat_formatters(ax)
    cb = fig.colorbar(sc, cax=cax)
    cb.set_ticks([-15, -10, -5, 0, 5, 10, 15])
    cb.ax.tick_params(labelsize=COLORBAR_TICK_FONT_SIZE)
    cb.set_label("Elevation relative to MSL (m)", fontsize=COLORBAR_LABEL_FONT_SIZE, fontfamily=FONT_FAMILY)

    annotate_panel_label(ax_context, "a")
    annotate_panel_label(ax, "b")

    fig.savefig(out_png, dpi=600, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    return out_png


def write_status(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8")


def main():
    global LAND_POLY_SHP
    args = parse_args()
    LAND_POLY_SHP = args.land_polygons.expanduser().resolve()
    models_dir = args.models_dir
    plot_dir = args.plot_dir or (models_dir / "model_plots")
    model_ids = select_model_ids(args)
    partition_lookup = {}
    adc_block_lookup: dict[str, str] = {}
    if args.partition_geojson.exists():
        try:
            partition_gdf = load_partitions(args.partition_geojson, args.domain_csv)
            partition_lookup = {
                str(row["model_domain_id"]): {
                    "geometry": clean_geometry(row.geometry),
                    "plot_rank": row.get("plot_rank", None),
                }
                for _, row in partition_gdf.iterrows()
            }
        except Exception:
            partition_lookup = {}
    if args.membership_csv.is_file():
        membership = pd.read_csv(args.membership_csv, dtype=str)
        required = {"model_domain_id", "block_id"}
        if not required.issubset(membership.columns):
            raise ValueError(
                f"{args.membership_csv} lacks {sorted(required - set(membership.columns))}"
            )
        if membership["model_domain_id"].duplicated().any():
            raise ValueError(f"Duplicate model_domain_id in {args.membership_csv}")
        adc_block_lookup = dict(
            zip(membership["model_domain_id"], membership["block_id"])
        )

    print("============================================================")
    print("Plot existing SFINCS partition models")
    print(f"Models dir : {models_dir}")
    print(f"Plot dir   : {plot_dir}")
    print(f"Models     : {len(model_ids)}")
    print("============================================================")

    status_rows: list[dict] = []
    for i, model_id in enumerate(model_ids, start=1):
        meta = partition_lookup.get(model_id, {})
        plot_rank = meta.get("plot_rank", None) if isinstance(meta, dict) else None
        partition_geom = meta.get("geometry", None) if isinstance(meta, dict) else None
        adc_block_id = adc_block_lookup.get(model_id)
        model_dir = find_partition_model_dir(models_dir, model_id, plot_rank)
        print(
            f"\n===== [{i}/{len(model_ids)}] "
            f"{partition_display_label(model_id, plot_rank, adc_block_id)} plot ====="
        )
        try:
            plot_path = create_sfincs_domain_forcing_plot(
                model_dir,
                model_id,
                plot_dir,
                partition_geom=partition_geom,
                plot_rank=plot_rank,
                adc_block_id=adc_block_id,
            )
            rec = {
                "model_domain_id": model_id,
                "plot_rank": plot_rank,
                "adc_block_id": adc_block_id,
                "status": "ok",
                "message": "",
                "model_dir": str(model_dir),
                "plot_path": str(plot_path),
            }
            print(f"  -> plot: {plot_path}")
        except Exception as exc:
            rec = {
                "model_domain_id": model_id,
                "plot_rank": plot_rank,
                "adc_block_id": adc_block_id,
                "status": "failed",
                "message": str(exc),
                "traceback": traceback.format_exc(),
                "model_dir": str(model_dir),
            }
            print(f"  !! plot failed: {exc}")
            if args.stop_on_error:
                status_rows.append(rec)
                write_status(models_dir / "sfincs_partition_plot_status.csv", status_rows)
                raise
        status_rows.append(rec)
        write_status(models_dir / "sfincs_partition_plot_status.csv", status_rows)

    print(f"\nDone. Plot status CSV:\n  {models_dir / 'sfincs_partition_plot_status.csv'}")


if __name__ == "__main__":
    main()
