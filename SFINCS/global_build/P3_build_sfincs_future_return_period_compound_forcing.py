#!/usr/bin/env python3
"""Build future return-period SFINCS compound-flood forcing for three SSPs.

This is a future-climate front end for
``P2_build_sfincs_historical_return_period_compound_forcing.py``. It deliberately reuses the
historical builder's SFINCS boundary matching, fort.63 dry-value convention,
fort.22 interpolation, MATLAB-equivalent C15/TCR wind field, and model-copy
workflow. Future-specific inputs are:

* P4-scaled ADCIRC cases under ``<return-period>yr_ssp126/245/370``;
* bias-corrected eight-model monthly t600 climatologies; CLIMADA TCR combines
  t600 with each track point's total surface Vmax to diagnose saturation
  q950, while C15 derives its gradient core from that same surface Vmax;
* 2061--2100 CaMa-Flood routed outflow for the same eight GCMs.

By default, every inlet cell is evaluated independently with the diagnostic
TC-landfall-month weighted P50 algorithm, then its eight GCM P50 values are
averaged with equal weights. Multiple SFINCS source points forming one inlet
section share that inlet flow equally. GTC-total and monthly modes remain
available explicitly.
conda activate sfincs_tcr
 python P3_build_sfincs_future_return_period_compound_forcing.py --return-period 200 --scenarios ssp126 ssp245 ssp370 --cama-quantile-mode tc_inlet_p90 --cama-quantile 0.90
"""

from __future__ import annotations

import argparse
import calendar
import copy
import logging
import os
import re
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import xarray as xr

import P2_build_sfincs_historical_return_period_compound_forcing as base


LOG = logging.getLogger("sfincs_future_compound_forcing")
SCRIPT_DIR = Path(__file__).resolve().parent


def configured_path(name: str, default: str | Path) -> Path:
    return Path(os.environ.get(name, str(default))).expanduser()

SCENARIOS = ("ssp126", "ssp245", "ssp370")
GCM_MODELS = (
    "CanESM5",
    "CMCC-CM2-SR5",
    "EC-Earth3",
    "GFDL-ESM4",
    "INM-CM4-8",
    "MPI-ESM1-2-LR",
    "MRI-ESM2-0",
    "NorESM2-LM",
)

DEFAULT_ADCIRC_RESULT_PARENT = configured_path(
    "ADCIRC_RESULT_PARENT", base.PROJECT_DIR / "external" / "adcirc" / "return_out"
)
DEFAULT_ADCIRC_RUN_PARENT = configured_path(
    "ADCIRC_RUN_PARENT",
    base.PROJECT_DIR / "external" / "adcirc" / "return_period_tc_event_reruns",
)
DEFAULT_FUTURE_CAMA_ROOT = configured_path(
    "CAMA_FUTURE_ROOT", base.PROJECT_DIR / "external" / "camaflood" / "output"
)
DEFAULT_DIAGNOSTIC_CACHE = configured_path(
    "CAMA_DIAGNOSTIC_CACHE", DEFAULT_FUTURE_CAMA_ROOT / "gtc_cama_inflow_figures" / "cache"
)
DEFAULT_FUTURE_T600_ROOT = configured_path(
    "CMIP6_T600_ROOT", base.PROJECT_DIR / "external" / "cmip6" / "tcr_environment"
)
DEFAULT_VLM_ADJUSTED_DEP_ROOT = configured_path(
    "VLM_ADJUSTED_DEP_ROOT",
    base.PROJECT_DIR / "external" / "vlm_data" / "gtc_terrain_overlay" / "adjusted_dep",
)

CAMA_NX = 1440
CAMA_NY = 720
CAMA_NCELL = CAMA_NX * CAMA_NY
CAMA_DTYPE = np.dtype("<f4")
CAMA_INVALID_LOW = -9990.0
CAMA_INVALID_HIGH = 1.0e19


ORIGINAL_PREFLIGHT = base.preflight
ORIGINAL_BUILD_CLIMADA_TRACK = base.build_climada_track
ORIGINAL_PREPARE_OUTPUT_MODEL = base.prepare_output_model

# Set for each scenario before base.run() enters the shared case builder.
ACTIVE_FUTURE_ARGS: argparse.Namespace | None = None


def build_parser() -> argparse.ArgumentParser:
    parser = base.build_parser()
    parser.description = (
        "Generate future return-period SFINCS compound-flood forcing for SSP126, "
        "SSP245 and SSP370."
    )
    parser.set_defaults(
        output_root=None,
        cama_root=DEFAULT_FUTURE_CAMA_ROOT,
        cama_start_year=2061,
        cama_end_year=2100,
        rain_model="tcr",
        tcr_wind_model="c15",
        cama_quantile_mode="tc_inlet_p50",
        cama_diagnostic_cache_dir=DEFAULT_DIAGNOSTIC_CACHE,
    )
    next(action for action in parser._actions if action.dest == "output_root").help = (
        "Common historical/future output root override. By default: "
        "global_sfincs_<return-period>yr_q<quantile>_tcr beside this script; "
        "each selected scenario is written below it."
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        choices=SCENARIOS,
        default=list(SCENARIOS),
        help="Future emissions scenarios to build.",
    )
    parser.add_argument(
        "--future-adcirc-result-parent",
        type=base.path_arg,
        default=DEFAULT_ADCIRC_RESULT_PARENT,
        help=(
            "Parent containing <return-period>yr_<scenario> ADCIRC result directories."
        ),
    )
    parser.add_argument(
        "--future-adcirc-run-parent",
        type=base.path_arg,
        default=DEFAULT_ADCIRC_RUN_PARENT,
        help=(
            "Parent containing P4 <return-period>yr_<scenario> ADCIRC run directories."
        ),
    )
    parser.add_argument(
        "--future-t600-root",
        "--future-q950-root",
        dest="future_t600_root",
        type=base.path_arg,
        default=DEFAULT_FUTURE_T600_ROOT,
        help=(
            "Directory containing the three bias-corrected eight-GCM t600 "
            "climatology files. --future-q950-root remains only as a "
            "backward-compatible option spelling."
        ),
    )
    parser.add_argument(
        "--future-cama-root",
        type=base.path_arg,
        default=DEFAULT_FUTURE_CAMA_ROOT,
        help="Parent containing cmip6_<GCM>_<scenario>_2061_2100 directories.",
    )
    parser.add_argument(
        "--vlm-adjusted-dep-root",
        type=base.path_arg,
        default=DEFAULT_VLM_ADJUSTED_DEP_ROOT,
        help=(
            "Directory containing GTC_XXXX/sfincs.dep produced by "
            "build_and_plot_continuous_global_vlm.py --write-adjusted-dep. "
            "Each newly built future case reuses its matching immutable terrain "
            "through the selected copy mode (hard link by default)."
        ),
    )
    parser.add_argument(
        "--no-vlm-adjustment",
        action="store_true",
        help=(
            "Explicitly retain the original base-model sfincs.dep. By default, "
            "future cases require and use the adjusted VLM terrain."
        ),
    )
    parser.add_argument(
        "--gcm-models",
        nargs="+",
        default=list(GCM_MODELS),
        help="GCMs used in the equal-weight future inflow ensemble.",
    )
    parser.add_argument(
        "--diagnostic-cache-dir",
        dest="cama_diagnostic_cache_dir",
        type=base.path_arg,
        default=DEFAULT_DIAGNOSTIC_CACHE,
        help="NPZ cache directory written by P1_plot_sfincs_inflow_diagnostics.py.",
    )
    parser.add_argument(
        "--no-diagnostic-cache",
        action="store_true",
        help="Ignore diagnostic NPZ files and read the original CaMa binaries.",
    )
    parser.add_argument(
        "--require-diagnostic-cache",
        dest="require_cama_diagnostic_cache",
        action="store_true",
        help=(
            "Do not fall back to the much slower raw CaMa scan when a complete "
            "P1_plot_sfincs_inflow_diagnostics.py cache is unavailable."
        ),
    )
    parser.add_argument(
        "--require-complete-adcirc",
        action="store_true",
        help=(
            "Stop unless every selected P4 run case has a non-empty transferred fort.63. "
            "By default, currently available results are built and missing transfers are reported."
        ),
    )
    return parser


def t600_path(root: Path, scenario: str, start_year: int, end_year: int) -> Path:
    return root / (
        f"t600_bias_corrected_multimodel_mean_{scenario}_"
        f"{start_year:04d}-{end_year:04d}.nc"
    )


def validate_t600_file(path: Path, scenario: str, models: Sequence[str]) -> list[str]:
    issues: list[str] = []
    if not path.is_file():
        return [f"Missing future t600 climatology: {path}"]
    try:
        with xr.open_dataset(path) as ds:
            if "t600" not in ds:
                issues.append(f"{path}: missing t600 variable")
            elif (
                "month" not in ds["t600"].dims
                or ds.sizes.get("month") != 12
            ):
                issues.append(
                    f"{path}: t600 must contain exactly 12 calendar months"
                )
            else:
                units = str(ds["t600"].attrs.get("units", "")).strip().lower()
                if units not in {"k", "kelvin"}:
                    issues.append(f"{path}: t600 units must be K, found {units!r}")
                values = ds["t600"]
                minimum = float(values.min())
                maximum = float(values.max())
                if minimum < 150.0 or maximum > 350.0:
                    issues.append(
                        f"{path}: implausible t600 range "
                        f"{minimum:.3f}..{maximum:.3f} K"
                    )
            file_scenario = str(ds.attrs.get("scenario", "")).lower()
            if file_scenario != scenario:
                issues.append(
                    f"{path}: scenario attribute {file_scenario!r} does not match {scenario}"
                )
            model_count = int(ds.attrs.get("model_count", -1))
            if model_count != len(models):
                issues.append(
                    f"{path}: model_count={model_count}, expected {len(models)}"
                )
            used = str(ds.attrs.get("models_used", "")).split()
            if set(used) != set(models):
                issues.append(f"{path}: models_used does not match requested 8-GCM ensemble")
    except Exception as exc:
        issues.append(f"Cannot validate future t600 {path}: {exc}")
    return issues


def _periodic_longitude(field: xr.DataArray, lon_name: str) -> xr.DataArray:
    lon = np.asarray(field[lon_name].values, dtype=float)
    if lon.ndim != 1 or lon.size < 2 or np.any(np.diff(lon) <= 0):
        raise ValueError("t600 longitude coordinate must be a strictly increasing vector")
    left = field.isel({lon_name: [-1]}).assign_coords({lon_name: [lon[-1] - 360.0]})
    right = field.isel({lon_name: [0]}).assign_coords({lon_name: [lon[0] + 360.0]})
    return xr.concat([left, field, right], dim=lon_name)


def sample_future_monthly_t600(
    path: Path,
    times: pd.DatetimeIndex,
    lon: np.ndarray,
    lat: np.ndarray,
) -> np.ndarray:
    """Sample the bias-corrected 12-month future t600 along a TC track."""
    if not path.is_file():
        raise FileNotFoundError(f"Missing future t600 climatology: {path}")
    result = np.full(len(times), np.nan, dtype=float)
    with xr.open_dataset(path) as ds:
        if "t600" not in ds:
            raise KeyError(f"{path}: missing t600 variable")
        temperature = ds["t600"]
        level_name = next(
            (
                name
                for name in ("level", "pressure_level", "plev")
                if name in temperature.dims
            ),
            None,
        )
        if level_name is not None:
            levels = np.asarray(temperature[level_name].values, dtype=float)
            target = 600.0 if np.nanmax(levels) < 2_000.0 else 60_000.0
            temperature = temperature.sel(
                {level_name: target}, method="nearest"
            )
        units = str(temperature.attrs.get("units", "K")).strip().lower()
        if units in {"c", "degc", "degree_celsius", "degreescelsius"}:
            temperature = temperature + 273.15
        elif units not in {"", "k", "kelvin"}:
            raise ValueError(f"{path}: unsupported t600 units {units!r}")
        lat_name = "latitude" if "latitude" in temperature.coords else "lat"
        lon_name = "longitude" if "longitude" in temperature.coords else "lon"
        if "month" not in temperature.dims or temperature.sizes["month"] != 12:
            raise ValueError(f"{path}: t600 must have a 12-record month dimension")
        if np.any(np.diff(temperature[lat_name].values) < 0):
            temperature = temperature.sortby(lat_name)
        source_longitude = np.mod(
            np.asarray(temperature[lon_name].values, dtype=float), 360.0
        )
        temperature = temperature.assign_coords(
            {lon_name: (lon_name, source_longitude)}
        ).sortby(lon_name)
        temperature = _periodic_longitude(temperature, lon_name)
        source_lon = np.asarray(temperature[lon_name].values, dtype=float)
        query_lon = np.asarray(lon, dtype=float)
        center = 0.5 * (source_lon[1] + source_lon[-2])
        query_lon = (query_lon - center + 180.0) % 360.0 + center - 180.0
        for month in np.unique(times.month):
            positions = np.flatnonzero(times.month == month)
            field = temperature.sel(month=int(month))
            sampled = field.interp(
                {
                    lon_name: xr.DataArray(query_lon[positions], dims="points"),
                    lat_name: xr.DataArray(np.asarray(lat)[positions], dims="points"),
                },
                method="linear",
            ).values
            result[positions] = np.asarray(sampled, dtype=float).reshape(-1)
    bad = ~np.isfinite(result) | (result < 150.0) | (result > 350.0)
    if bad.any():
        raise ValueError(
            f"Invalid future t600 at {bad.sum()} track positions from {path}"
        )
    return result


def read_future_track_points(meta_path: Path) -> pd.DataFrame:
    """Read P4 FutureForcingTrackPoints and expose the fields expected by TCR."""
    lines = meta_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    header_index = next(
        (
            i
            for i, line in enumerate(lines)
            if line.strip().startswith("index,original_time_index,time,lon180,lat,")
        ),
        None,
    )
    if header_index is None:
        raise ValueError(f"No [FutureForcingTrackPoints] CSV table in {meta_path}")
    frame = pd.read_csv(base.io.StringIO("\n".join(lines[header_index:])))
    required = {
        "original_time_index",
        "time",
        "lon180",
        "lat",
        "future_vmax_ms",
        "future_pressure_hPa",
    }
    if not required.issubset(frame.columns):
        raise ValueError(f"{meta_path} lacks {sorted(required - set(frame.columns))}")
    # Preserve both historical and future surface-intensity columns so the
    # shared builder can verify provenance against source vmax_trks. The
    # aliases retain the historical-builder interface.
    frame["vmax_ms"] = frame["future_vmax_ms"]
    frame["pressure_hPa"] = frame["future_pressure_hPa"]
    frame["time"] = pd.to_datetime(frame["time"])
    numeric = required - {"time"}
    for name in numeric | {"vmax_ms", "pressure_hPa"}:
        if name in frame:
            frame[name] = pd.to_numeric(frame[name], errors="coerce")
    needed = [
        "original_time_index",
        "time",
        "lon180",
        "lat",
        "historical_vmax_ms",
        "future_vmax_ms",
        "vmax_ms",
        "pressure_hPa",
    ]
    return frame.dropna(subset=needed).sort_values("time").reset_index(drop=True)


def resolve_future_event(
    block_id: str,
    result_dir: Path,
    run_root: Path,
    mesh_root: Path,
) -> base.Event:
    event_id = result_dir.name
    run_dir = run_root / block_id / event_id
    future_meta_path = run_dir / "fort22_meta.txt"
    if not future_meta_path.is_file():
        raise FileNotFoundError(f"Missing future event metadata: {future_meta_path}")
    future_meta = base.parse_key_value_file(future_meta_path)
    historical_meta_path = Path(
        future_meta.get(
            "source_historical_fort22_meta",
            str(Path(future_meta.get("source_case_dir", "")) / "fort22_meta.txt"),
        )
    )
    historical_meta = (
        base.parse_key_value_file(historical_meta_path)
        if historical_meta_path.is_file()
        else {}
    )
    meta = {**historical_meta, **future_meta}
    start = base.parse_timestamp(meta, "WindowStart")
    stop = base.parse_timestamp(meta, "WindowEnd")
    if start is None or stop is None or stop <= start:
        raise ValueError(f"Invalid WindowStart/WindowEnd in {future_meta_path}")
    impact = base.parse_timestamp(meta, "TC_aligned_impact_time")
    dt_seconds = int(round(float(meta.get("WTIMINC", "3600"))))
    fort14 = run_dir / "fort.14"
    if not fort14.is_file():
        fort14 = mesh_root / block_id / "fort.14"
    return base.Event(
        block_id=block_id,
        event_id=event_id,
        result_dir=result_dir,
        run_dir=run_dir,
        fort63=result_dir / "fort.63",
        fort14=fort14,
        fort22=run_dir / "fort.22",
        meta_path=future_meta_path,
        meta=meta,
        start=start,
        stop=stop,
        impact=impact,
        dt_seconds=dt_seconds,
    )


def discover_future_events(
    args: argparse.Namespace, selected_blocks: set[str] | None
) -> list[base.Event]:
    selected_events = base.csv_set(args.events)
    expected: set[tuple[str, str]] = set()
    available: set[tuple[str, str]] = set()
    available_results: list[tuple[str, Path]] = []
    for block_dir in sorted(args.adcirc_run_root.glob("ADC_*")):
        if not block_dir.is_dir():
            continue
        block_id = block_dir.name
        if selected_blocks and block_id not in selected_blocks:
            continue
        for run_case in sorted(block_dir.iterdir()):
            if not run_case.is_dir() or not (run_case / "fort22_meta.txt").is_file():
                continue
            if selected_events and run_case.name not in selected_events:
                continue
            key = (block_id, run_case.name)
            expected.add(key)
            result_dir = args.adcirc_result_root / block_id / run_case.name
            fort63 = result_dir / "fort.63"
            if not fort63.is_file() or fort63.stat().st_size < 16:
                continue
            available.add(key)
            available_results.append((block_id, result_dir))
    missing = sorted(expected - available)
    args.future_expected_event_count = len(expected)
    args.future_available_event_count = len(available)
    args.future_missing_event_count = len(missing)
    LOG.info(
        "%s ADCIRC transfer inventory: expected=%s, available non-empty fort.63=%s, missing=%s",
        args.climate_scenario,
        len(expected),
        len(available),
        len(missing),
    )
    if missing:
        LOG.warning(
            "%s ADCIRC transfer is incomplete; first missing cases: %s",
            args.climate_scenario,
            "; ".join(f"{block}/{event}" for block, event in missing[:8]),
        )
    if args.require_complete_adcirc and missing:
        raise RuntimeError(
            f"{args.climate_scenario}: {len(missing)} selected ADCIRC fort.63 files are not transferred"
        )
    if args.limit_events is not None:
        available_results = available_results[: args.limit_events]
    events: list[base.Event] = []
    for block_id, result_dir in available_results:
        try:
            events.append(
                resolve_future_event(
                    block_id, result_dir, args.adcirc_run_root, args.adcirc_mesh_root
                )
            )
        except Exception as exc:
            if args.continue_on_error:
                LOG.error("Skip invalid future event %s/%s: %s", block_id, result_dir.name, exc)
                continue
            raise
    return events


def future_climada_track(
    event: base.Event, args: argparse.Namespace
) -> xr.Dataset | None:
    track = ORIGINAL_BUILD_CLIMADA_TRACK(event, args)
    if track is not None:
        if args.rain_model == "tcr":
            required_shear = {"u250", "v250", "u850", "v850", "ushear", "vshear"}
            missing_shear = sorted(required_shear - set(track.variables))
            if missing_shear:
                raise ValueError(
                    "Future TCR track lacks direct 250-850 hPa shear variables: "
                    + ", ".join(missing_shear)
                )
            np.testing.assert_allclose(
                track["ushear"].values,
                track["u250"].values - track["u850"].values,
                rtol=0.0,
                atol=0.0,
                err_msg="Future TCR zonal shear is not u250-u850",
            )
            np.testing.assert_allclose(
                track["vshear"].values,
                track["v250"].values - track["v850"].values,
                rtol=0.0,
                atol=0.0,
                err_msg="Future TCR meridional shear is not v250-v850",
            )
        if args.rain_model == "tcr" and args.tcr_wind_model == "c15":
            required_c15 = {
                "c15_vmax_gradient",
                "c15_vmax_surface_core",
                "vmax_total_surface",
                "vmax_historical_surface",
                "vmax_source_track_surface",
                "translation_speed",
                "future_intensity_ratio",
            }
            missing_c15 = sorted(required_c15 - set(track.variables))
            if missing_c15:
                raise ValueError(
                    "Future C15 TCR track lacks wind-level diagnostics: "
                    + ", ".join(missing_c15)
                )
            expected_surface_core = np.maximum(
                track["vmax_total_surface"].values
                - track["translation_speed"].values,
                0.0,
            )
            np.testing.assert_allclose(
                track["c15_vmax_surface_core"].values,
                expected_surface_core,
                rtol=1.0e-12,
                atol=1.0e-12,
                err_msg=(
                    "Future C15 surface core is inconsistent with surface "
                    "vmax_trks minus translation"
                ),
            )
            np.testing.assert_allclose(
                track["c15_vmax_gradient"].values,
                expected_surface_core
                / base.TCR_GRADIENT_TO_SURFACE_WIND_FACTOR,
                rtol=1.0e-12,
                atol=1.0e-12,
                err_msg=(
                    "Future C15 gradient wind is inconsistent with the surface "
                    "vmax_trks to gradient conversion"
                ),
            )
            np.testing.assert_allclose(
                track["max_sustained_wind"].values,
                track["vmax_total_surface"].values,
                rtol=1.0e-12,
                atol=1.0e-12,
                err_msg=(
                    "Future C15 TCR q950 interface did not receive future "
                    "total surface vmax_trks"
                ),
            )
        track.attrs["climate_scenario"] = args.climate_scenario
        track.attrs["future_intensity_scale_factor"] = float(
            event.meta.get("intensity_scale_factor", "nan")
        )
        track.attrs["future_t600_source"] = str(args.t600_source)
        track.attrs["tcr_humidity_pathway"] = (
            "CLIMADA saturation q950 diagnosed from bias-corrected t600 and "
            "future total surface vmax_trks; q950 intentionally absent"
        )
        track.attrs["future_tcr_wind_and_shear_contract"] = (
            "total surface vmax_trks at q950 interface; C15 profile selected "
            "with gradient_core=max(vmax_trks-translation,0)/0.9 and used "
            "directly as TCR gradient wind; direct u250-u850/v250-v850 shear"
        )
    return track


def future_case_dir(args: argparse.Namespace, model: str, scenario: str) -> Path:
    return args.future_cama_root / (
        f"cmip6_{model}_{scenario}_{args.cama_start_year:04d}_"
        f"{args.cama_end_year:04d}"
    )


def future_year_files(args: argparse.Namespace, model: str, scenario: str) -> list[Path]:
    directory = future_case_dir(args, model, scenario)
    return [directory / f"outflw{year}.bin" for year in range(args.cama_start_year, args.cama_end_year + 1)]


def validate_future_flow_inputs(args: argparse.Namespace, scenario: str) -> list[str]:
    issues: list[str] = []
    for model in args.gcm_models:
        for year, path in zip(
            range(args.cama_start_year, args.cama_end_year + 1),
            future_year_files(args, model, scenario),
            strict=True,
        ):
            if not path.is_file() or path.stat().st_size == 0:
                issues.append(f"Missing future CaMa file: {path}")
                continue
            record_bytes = CAMA_NCELL * CAMA_DTYPE.itemsize
            if path.stat().st_size % record_bytes:
                issues.append(f"Invalid CaMa file size: {path} ({path.stat().st_size} bytes)")
                continue
            nstep = path.stat().st_size // record_bytes
            days = 366 if calendar.isleap(year) else 365
            if nstep not in (days, 4 * days, 8 * days):
                issues.append(f"Unexpected CaMa time count in {path}: {nstep}")
    return issues


def q50_cache_path(args: argparse.Namespace, scenario: str) -> Path:
    if args.q50_cache is None:
        if args.cama_quantile_mode in {"tc_inlet_p50", "tc_inlet_p90"}:
            return args.output_root / "_cache" / (
                "cama_8gcm_independent_inlet_tc_weighted_"
                f"p{int(round(args.cama_quantile * 100)):02d}_{scenario}.csv"
            )
        if args.cama_quantile_mode == "tc_summary_p50":
            return args.output_root / "_cache" / (
                f"cama_8gcm_tc_p50_from_diagnostic_summary_{scenario}.csv"
            )
        if args.cama_quantile_mode == "tc_monthly_mean":
            return args.output_root / "_cache" / (
                f"cama_8gcm_tc_landfall_month_weighted_mean_{scenario}.csv"
            )
        return args.output_root / "_cache" / (
            f"cama_8gcm_tc_landfall_month_weighted_q50_{scenario}_"
            f"{args.cama_start_year}_{args.cama_end_year}.csv"
        )
    text = str(args.q50_cache)
    if "{scenario}" in text:
        return Path(text.format(scenario=scenario))
    # base.resolve_runtime_paths supplies one generic cache name before this
    # future hook runs. Always make it scenario-specific, including when this
    # process was launched for only one scenario; otherwise two block-sharded
    # workers (or separate SSP jobs) race on and overwrite the same CSV.
    if args.q50_cache.stem.casefold().endswith(f"_{scenario}".casefold()):
        return args.q50_cache
    return args.q50_cache.with_name(
        f"{args.q50_cache.stem}_{scenario}{args.q50_cache.suffix or '.csv'}"
    )


def _expected_dates(start_year: int, end_year: int) -> np.ndarray:
    return pd.date_range(
        f"{start_year:04d}-01-01", f"{end_year:04d}-12-31", freq="D"
    ).values.astype("datetime64[D]")


def load_diagnostic_values(
    args: argparse.Namespace,
    model: str,
    scenario: str,
    cells: Sequence[tuple[int, int]],
) -> np.ndarray | None:
    if args.no_diagnostic_cache:
        return None
    case_name = (
        f"cmip6_{model}_{scenario}_{args.cama_start_year:04d}_{args.cama_end_year:04d}"
    )
    loaded = base.load_cama_diagnostic_values(
        args.cama_diagnostic_cache_dir,
        case_name,
        cells,
        args.cama_start_year,
        args.cama_end_year,
    )
    if loaded is None:
        return None
    _, values, path = loaded
    LOG.info("Use complete diagnostic CaMa cache: %s", path)
    return values


def read_future_cama_cells(
    args: argparse.Namespace,
    model: str,
    scenario: str,
    cells: Sequence[tuple[int, int]],
) -> np.ndarray:
    cached = load_diagnostic_values(args, model, scenario, cells)
    if cached is not None:
        return cached
    if args.require_cama_diagnostic_cache:
        raise FileNotFoundError(
            f"No complete diagnostic cache covering all requested cells for "
            f"cmip6_{model}_{scenario}_{args.cama_start_year}_{args.cama_end_year}"
        )
    flat = np.asarray([row * CAMA_NX + col for row, col in cells], dtype=np.int64)
    yearly: list[np.ndarray] = []
    for year, path in zip(
        range(args.cama_start_year, args.cama_end_year + 1),
        future_year_files(args, model, scenario),
        strict=True,
    ):
        days = 366 if calendar.isleap(year) else 365
        record_bytes = CAMA_NCELL * CAMA_DTYPE.itemsize
        if not path.is_file() or path.stat().st_size % record_bytes:
            raise FileNotFoundError(f"Missing or invalid future CaMa file: {path}")
        nstep = path.stat().st_size // record_bytes
        if nstep not in (days, 4 * days, 8 * days):
            raise ValueError(f"Unexpected time count {nstep} in {path}")
        LOG.info("  %s %s %s", scenario, model, path.name)
        mm = np.memmap(path, dtype=CAMA_DTYPE, mode="r", shape=(nstep, CAMA_NCELL))
        values = np.asarray(mm[:, flat], dtype=np.float32)
        del mm
        bad = (
            ~np.isfinite(values)
            | (values <= CAMA_INVALID_LOW)
            | (values >= CAMA_INVALID_HIGH)
            | (values < 0)
        )
        values[bad] = np.nan
        steps_per_day = nstep // days
        if steps_per_day > 1:
            values = np.nanmean(
                values.reshape(days, steps_per_day, len(cells)), axis=1
            ).astype(np.float32)
        yearly.append(values)
    return np.vstack(yearly)


def load_or_compute_future_q50(
    args: argparse.Namespace,
    cells: Sequence[tuple[int, int]],
    allow_compute: bool = True,
    contexts: dict[str, list[base.ModelContext]] | None = None,
) -> pd.DataFrame:
    if contexts is None:
        raise ValueError("Future TC-weighted q50 requires GTC model contexts")
    if args.cama_quantile_mode == "tc_summary_p50":
        return load_future_summary_p50(args, contexts, allow_compute=allow_compute)
    if args.cama_quantile_mode == "tc_monthly_mean":
        return load_future_month_weighted_mean(
            args, contexts, allow_compute=allow_compute
        )
    scenario = args.climate_scenario
    cache = q50_cache_path(args, scenario)
    independent_inlets = args.cama_quantile_mode in {"tc_inlet_p50", "tc_inlet_p90"}
    percentile_label = f"p{int(round(args.cama_quantile * 100)):02d}"
    ensemble_method = (
        "equal_mean_of_eight_gcm_independent_inlet_tc_landfall_month_weighted_"
        + percentile_label
        if independent_inlets
        else "equal_mean_of_eight_gcm_per_gtc_tc_landfall_month_weighted_total_q50_allocated_by_cell_quantile"
    )
    river_contexts = base.unique_river_contexts(contexts)
    weights_by_model = base.load_gtc_month_weights(
        args.cama_month_weights_csv, sorted(river_contexts)
    )
    wanted = {
        (model_id, int(row.cama_row), int(row.cama_col))
        for model_id, context in river_contexts.items()
        for row in context.src_mapping[["cama_row", "cama_col"]]
        .drop_duplicates()
        .itertuples(index=False)
    }
    required = {
        "model_id", "cama_row", "cama_col", "q_m3s", "scenario",
        "gtc_total_q_m3s", "allocation_fraction", "model_count", "models_used",
        "ensemble_method", "month_weights_csv",
    }
    if cache.is_file() and not args.recompute_q50:
        table = pd.read_csv(cache)
        if required.issubset(table.columns):
            cached_keys = set(
                zip(
                    table.model_id.astype(str),
                    table.cama_row.astype(int),
                    table.cama_col.astype(int),
                )
            )
            correct_scenario = set(table.scenario.astype(str)) == {scenario}
            correct_count = set(pd.to_numeric(table.model_count, errors="coerce")) == {
                len(args.gcm_models)
            }
            correct_models = set(table.models_used.astype(str)) == {" ".join(args.gcm_models)}
            correct_method = set(table.ensemble_method.astype(str)) == {ensemble_method}
            correct_weights = set(table.month_weights_csv.astype(str)) == {
                str(args.cama_month_weights_csv)
            }
            if (
                wanted.issubset(cached_keys)
                and correct_scenario
                and correct_count
                and correct_models
                and correct_method
                and correct_weights
            ):
                LOG.info("Use cached future 8-GCM TC-month-weighted q50: %s", cache)
                return table[
                    table.apply(
                        lambda row: (
                            str(row.model_id), int(row.cama_row), int(row.cama_col)
                        ) in wanted,
                        axis=1,
                    )
                ].copy()
    columns = [
        "model_id", "cama_row", "cama_col", "q_m3s",
        "gtc_total_q_m3s", "allocation_fraction", "valid_count", "expected_count",
    ]
    if not allow_compute or not wanted:
        return pd.DataFrame(columns=columns)

    cell_index = {cell: index for index, cell in enumerate(cells)}
    row_keys = sorted(wanted)
    row_key_index = {key: index for index, key in enumerate(row_keys)}
    months = pd.DatetimeIndex(
        _expected_dates(args.cama_start_year, args.cama_end_year)
    ).month.to_numpy(dtype=int)
    model_q50: list[np.ndarray] = []
    model_counts: list[np.ndarray] = []
    expected_count = len(months)
    for index, model in enumerate(args.gcm_models, start=1):
        LOG.info(
            "%s future CaMa q50 model %s/%s: %s",
            scenario,
            index,
            len(args.gcm_models),
            model,
        )
        values = read_future_cama_cells(args, model, scenario, cells)
        if values.shape != (expected_count, len(cells)):
            raise ValueError(
                f"{scenario}/{model}: extracted shape {values.shape}, "
                f"expected {(expected_count, len(cells))}"
            )
        q_values = np.full(len(row_keys), np.nan, dtype=float)
        valid_counts = np.zeros(len(row_keys), dtype=np.int64)
        for model_id, context in river_contexts.items():
            model_cells = [
                (int(item.cama_row), int(item.cama_col))
                for item in context.src_mapping[["cama_row", "cama_col"]]
                .drop_duplicates()
                .itertuples(index=False)
            ]
            model_values = values[
                :, [cell_index[cell] for cell in model_cells]
            ].astype(float)
            any_valid = np.isfinite(model_values).any(axis=1)
            total_series = np.nansum(model_values, axis=1)
            total_series[~any_valid] = np.nan
            weights = weights_by_model[model_id]
            total_q = base.tc_weighted_quantile(
                total_series, months, weights, args.cama_quantile
            )
            cell_q: list[float] = []
            cell_counts: list[int] = []
            for cell_i, cell in enumerate(model_cells):
                series = model_values[:, cell_i]
                valid_count = int(np.count_nonzero(np.isfinite(series)))
                if (
                    not independent_inlets
                    and valid_count / expected_count < args.min_cama_valid_fraction
                ):
                    raise ValueError(
                        f"{scenario}/{model}/{model_id}/{cell}: "
                        "CaMa valid fraction is too low"
                    )
                cell_q.append(
                    base.tc_weighted_quantile(
                        series, months, weights, args.cama_quantile
                    )
                )
                cell_counts.append(valid_count)
            cell_q_array = np.asarray(cell_q, dtype=float)
            if not np.isfinite(total_q) or not np.all(np.isfinite(cell_q_array)):
                raise ValueError(
                    f"{scenario}/{model}/{model_id}: non-finite TC-month-weighted q50"
                )
            raw_sum = float(cell_q_array.sum())
            fractions = (
                cell_q_array / raw_sum
                if raw_sum > 0
                else np.full(len(model_cells), 1.0 / len(model_cells))
            )
            if independent_inlets:
                allocated = cell_q_array
            else:
                allocated = total_q * fractions
                if not np.isclose(float(allocated.sum()), total_q):
                    raise RuntimeError(
                        f"{scenario}/{model}/{model_id}: discharge allocation is not conservative"
                    )
            for cell, q_value, valid_count in zip(
                model_cells, allocated, cell_counts, strict=True
            ):
                output_i = row_key_index[(model_id, cell[0], cell[1])]
                q_values[output_i] = q_value
                valid_counts[output_i] = valid_count
        if not np.all(np.isfinite(q_values)):
            raise ValueError(f"{scenario}/{model}: non-finite TC-month-weighted q50")
        model_q50.append(q_values)
        model_counts.append(valid_counts)
        del values

    q_stack = np.stack(model_q50)
    if q_stack.shape[0] != len(args.gcm_models) or not np.all(np.isfinite(q_stack)):
        raise ValueError(f"{scenario}: incomplete or non-finite GCM q50 ensemble")
    ensemble = np.mean(q_stack, axis=0)
    count_stack = np.stack(model_counts)
    ensemble_totals: dict[str, float] = {}
    for row_index, key in enumerate(row_keys):
        ensemble_totals[key[0]] = ensemble_totals.get(key[0], 0.0) + float(
            ensemble[row_index]
        )
    allocation_fractions = np.asarray(
        [
            float(ensemble[index]) / ensemble_totals[key[0]]
            if ensemble_totals[key[0]] > 0
            else 1.0 / sum(candidate[0] == key[0] for candidate in row_keys)
            for index, key in enumerate(row_keys)
        ],
        dtype=float,
    )
    result = pd.DataFrame(
        {
            "model_id": [key[0] for key in row_keys],
            "cama_row": [key[1] for key in row_keys],
            "cama_col": [key[2] for key in row_keys],
            "cama_lon": [-180.0 + (key[2] + 0.5) * 0.25 for key in row_keys],
            "cama_lat": [90.0 - (key[1] + 0.5) * 0.25 for key in row_keys],
            "q_m3s": ensemble,
            "gtc_total_q_m3s": [ensemble_totals[key[0]] for key in row_keys],
            "allocation_fraction": allocation_fractions,
            "valid_count": np.min(count_stack, axis=0),
            "expected_count": expected_count,
            "quantile_nonexceedance": args.cama_quantile,
            "start_year": args.cama_start_year,
            "end_year": args.cama_end_year,
            "scenario": scenario,
            "model_count": len(args.gcm_models),
            "models_used": " ".join(args.gcm_models),
            "ensemble_method": ensemble_method,
            "month_weights_csv": str(args.cama_month_weights_csv),
        }
    )
    for model, q_values in zip(args.gcm_models, model_q50, strict=True):
        safe_model = re.sub(r"[^A-Za-z0-9]+", "_", model).strip("_")
        result[f"q{int(round(args.cama_quantile * 100)):02d}_m3s__{safe_model}"] = q_values
    cache.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(cache, index=False, encoding="utf-8-sig")
    LOG.info("Wrote future 8-GCM TC-month-weighted q50 cache: %s", cache)
    return result


def load_future_summary_p50(
    args: argparse.Namespace,
    contexts: dict[str, list[base.ModelContext]],
    allow_compute: bool = True,
) -> pd.DataFrame:
    scenario = args.climate_scenario
    river_contexts = base.unique_river_contexts(contexts)
    needed = {
        (model_id, int(item.cama_row), int(item.cama_col))
        for model_id, context in river_contexts.items()
        for item in context.src_mapping[["cama_row", "cama_col"]]
        .drop_duplicates()
        .itertuples(index=False)
    }
    cache = q50_cache_path(args, scenario)
    definition = "equal_mean_of_eight_gcm_direct_gtc_tc_p50_from_diagnostic_summary"
    required = {
        "model_id", "cama_row", "cama_col", "q_m3s", "gtc_total_q_m3s",
        "allocation_fraction", "q_reference_m3s", "scenario", "model_count",
        "models_used", "definition", "flow_summary_csv", "inlet_cells_csv",
    }
    if cache.is_file() and not args.recompute_q50:
        cached = pd.read_csv(cache)
        if required.issubset(cached.columns):
            keys = set(
                zip(
                    cached.model_id.astype(str),
                    cached.cama_row.astype(int),
                    cached.cama_col.astype(int),
                )
            )
            valid = (
                needed.issubset(keys)
                and set(cached.scenario.astype(str)) == {scenario}
                and set(cached.definition.astype(str)) == {definition}
                and set(pd.to_numeric(cached.model_count, errors="coerce"))
                == {len(args.gcm_models)}
                and set(cached.models_used.astype(str)) == {" ".join(args.gcm_models)}
                and set(cached.flow_summary_csv.astype(str))
                == {str(args.cama_flow_summary_csv)}
                and set(cached.inlet_cells_csv.astype(str))
                == {str(args.cama_inlet_cells_csv)}
            )
            if valid:
                LOG.info("Use cached future direct-summary TC P50: %s", cache)
                return cached[
                    cached.apply(
                        lambda row: (
                            str(row.model_id), int(row.cama_row), int(row.cama_col)
                        ) in needed,
                        axis=1,
                    )
                ].copy()
    columns = [
        "model_id", "cama_row", "cama_col", "q_m3s", "gtc_total_q_m3s",
        "allocation_fraction",
    ]
    if not allow_compute or not needed:
        return pd.DataFrame(columns=columns)

    per_gcm = base.diagnostic_summary_p50_totals(
        args.cama_flow_summary_csv,
        sorted(river_contexts),
        scenario,
        args.gcm_models,
    ).rename(
        columns={
            "gtc_id": "model_id",
            "model": "climate_model",
            "tc_p50_flow_m3s": "gcm_total_p50_m3s",
        }
    )
    model_counts = per_gcm.groupby("model_id")["climate_model"].nunique()
    if not model_counts.eq(len(args.gcm_models)).all():
        bad = model_counts[~model_counts.eq(len(args.gcm_models))]
        raise ValueError(f"Future GTC summaries do not contain all eight GCMs:\n{bad}")
    ensemble_totals = per_gcm.groupby("model_id", as_index=False).agg(
        gtc_total_q_m3s=("gcm_total_p50_m3s", "mean")
    )
    fractions = base.diagnostic_inlet_reference_fractions(
        args.cama_inlet_cells_csv, needed
    ).rename(columns={"gtc_id": "model_id"})
    result = fractions.merge(
        ensemble_totals,
        on="model_id",
        how="left",
        validate="many_to_one",
    )
    if result["gtc_total_q_m3s"].isna().any():
        raise KeyError("Future diagnostic summary lacks one or more required GTCs")
    result["q_m3s"] = result["gtc_total_q_m3s"] * result["allocation_fraction"]
    result["cama_lon"] = -180.0 + (result["cama_col"] + 0.5) * 0.25
    result["cama_lat"] = 90.0 - (result["cama_row"] + 0.5) * 0.25
    result["scenario"] = scenario
    result["model_count"] = len(args.gcm_models)
    result["models_used"] = " ".join(args.gcm_models)
    result["definition"] = definition
    result["allocation_method"] = "q_reference_m3s_fraction"
    result["flow_summary_csv"] = str(args.cama_flow_summary_csv)
    result["inlet_cells_csv"] = str(args.cama_inlet_cells_csv)
    pivot = per_gcm.pivot(
        index="model_id", columns="climate_model", values="gcm_total_p50_m3s"
    )
    for model in args.gcm_models:
        safe_model = re.sub(r"[^A-Za-z0-9]+", "_", model).strip("_")
        result[f"gtc_total_p50_m3s__{safe_model}"] = result["model_id"].map(
            pivot[model]
        )
    cache.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(cache, index=False, encoding="utf-8-sig")
    LOG.info("Wrote future direct-summary TC P50 cache: %s", cache)
    return result


def load_future_month_weighted_mean(
    args: argparse.Namespace,
    contexts: dict[str, list[base.ModelContext]],
    allow_compute: bool = True,
) -> pd.DataFrame:
    scenario = args.climate_scenario
    river_contexts = base.unique_river_contexts(contexts)
    needed = {
        (model_id, int(item.cama_row), int(item.cama_col))
        for model_id, context in river_contexts.items()
        for item in context.src_mapping[["cama_row", "cama_col"]]
        .drop_duplicates()
        .itertuples(index=False)
    }
    cache = q50_cache_path(args, scenario)
    definition = (
        "equal_mean_of_eight_gcm_sum_of_gtc_landfall_month_probability_"
        "times_monthly_inlet_mean"
    )
    required = {
        "model_id", "cama_row", "cama_col", "q_m3s", "gtc_total_q_m3s",
        "allocation_fraction", "scenario", "model_count", "models_used",
        "definition", "month_weights_csv", "monthly_inlet_csv",
    }
    if cache.is_file() and not args.recompute_q50:
        cached = pd.read_csv(cache)
        if required.issubset(cached.columns):
            keys = set(
                zip(
                    cached.model_id.astype(str),
                    cached.cama_row.astype(int),
                    cached.cama_col.astype(int),
                )
            )
            valid = (
                needed.issubset(keys)
                and set(cached.scenario.astype(str)) == {scenario}
                and set(cached.definition.astype(str)) == {definition}
                and set(pd.to_numeric(cached.model_count, errors="coerce"))
                == {len(args.gcm_models)}
                and set(cached.models_used.astype(str)) == {" ".join(args.gcm_models)}
                and set(cached.month_weights_csv.astype(str))
                == {str(args.cama_month_weights_csv)}
                and set(cached.monthly_inlet_csv.astype(str))
                == {str(args.cama_monthly_inlet_csv)}
            )
            if valid:
                LOG.info("Use cached future TC-month-weighted mean flow: %s", cache)
                return cached[
                    cached.apply(
                        lambda row: (
                            str(row.model_id), int(row.cama_row), int(row.cama_col)
                        ) in needed,
                        axis=1,
                    )
                ].copy()
    columns = [
        "model_id", "cama_row", "cama_col", "q_m3s", "gtc_total_q_m3s",
        "allocation_fraction",
    ]
    if not allow_compute or not needed:
        return pd.DataFrame(columns=columns)

    per_gcm = base.diagnostic_month_weighted_cell_flows(
        args.cama_monthly_inlet_csv,
        args.cama_month_weights_csv,
        sorted(river_contexts),
        scenario,
        args.gcm_models,
    ).rename(columns={"gtc_id": "model_id", "model": "climate_model"})
    per_gcm_keys = set(zip(per_gcm.model_id, per_gcm.cama_row, per_gcm.cama_col))
    missing = sorted(needed - per_gcm_keys)
    if missing:
        raise KeyError(f"Diagnostic monthly table lacks future inlet cells: {missing[:20]}")
    per_gcm = per_gcm[
        per_gcm.apply(
            lambda row: (str(row.model_id), int(row.cama_row), int(row.cama_col))
            in needed,
            axis=1,
        )
    ].copy()
    group_keys = ["model_id", "cama_row", "cama_col"]
    model_counts = per_gcm.groupby(group_keys)["climate_model"].nunique()
    if not model_counts.eq(len(args.gcm_models)).all():
        bad = model_counts[~model_counts.eq(len(args.gcm_models))].head(20)
        raise ValueError(f"Future inlet cells do not contain all eight GCMs:\n{bad}")
    ensemble = per_gcm.groupby(group_keys, as_index=False).agg(q_m3s=("q_m3s", "mean"))
    totals = ensemble.groupby("model_id")["q_m3s"].transform("sum")
    counts = ensemble.groupby("model_id")["model_id"].transform("size")
    ensemble["gtc_total_q_m3s"] = totals
    ensemble["allocation_fraction"] = np.where(
        totals > 0, ensemble["q_m3s"] / totals, 1.0 / counts
    )
    ensemble["cama_lon"] = -180.0 + (ensemble["cama_col"] + 0.5) * 0.25
    ensemble["cama_lat"] = 90.0 - (ensemble["cama_row"] + 0.5) * 0.25
    ensemble["scenario"] = scenario
    ensemble["model_count"] = len(args.gcm_models)
    ensemble["models_used"] = " ".join(args.gcm_models)
    ensemble["definition"] = definition
    ensemble["month_weights_csv"] = str(args.cama_month_weights_csv)
    ensemble["monthly_inlet_csv"] = str(args.cama_monthly_inlet_csv)
    pivot = per_gcm.pivot(
        index=group_keys, columns="climate_model", values="q_m3s"
    )
    for model in args.gcm_models:
        safe_model = re.sub(r"[^A-Za-z0-9]+", "_", model).strip("_")
        lookup = pivot[model].to_dict()
        ensemble[f"weighted_mean_m3s__{safe_model}"] = [
            lookup[(row.model_id, row.cama_row, row.cama_col)]
            for row in ensemble.itertuples(index=False)
        ]
    cache.parent.mkdir(parents=True, exist_ok=True)
    ensemble.to_csv(cache, index=False, encoding="utf-8-sig")
    LOG.info("Wrote future TC-month-weighted mean flow cache: %s", cache)
    return ensemble


def future_preflight(
    args: argparse.Namespace,
    events: Sequence[base.Event],
    dry_run: bool,
) -> list[str]:
    issues = ORIGINAL_PREFLIGHT(args, events, dry_run)
    if not args.no_vlm_adjustment:
        membership = base.load_membership(args)
        model_ids = sorted(set(membership["model_domain_id"].astype(str)))
        if not args.vlm_adjusted_dep_root.is_dir():
            issues.append(
                "Missing adjusted VLM terrain root: "
                f"{args.vlm_adjusted_dep_root}. First run "
                "build_and_plot_continuous_global_vlm.py --gtc-only "
                "--write-adjusted-dep."
            )
        else:
            for model_id in model_ids:
                source = args.vlm_adjusted_dep_root / model_id / "sfincs.dep"
                base_dep = args.model_root / model_id / "sfincs.dep"
                if not source.is_file():
                    issues.append(f"Missing adjusted VLM terrain for {model_id}: {source}")
                    continue
                if base_dep.is_file() and source.stat().st_size != base_dep.stat().st_size:
                    issues.append(
                        f"Adjusted VLM terrain size mismatch for {model_id}: "
                        f"{source.stat().st_size} bytes versus base "
                        f"{base_dep.stat().st_size} bytes"
                    )
    issues.extend(
        validate_t600_file(
            args.t600_source, args.climate_scenario, args.gcm_models
        )
    )
    if args.cama_quantile_mode not in {"tc_summary_p50", "tc_monthly_mean"}:
        flow_issues = validate_future_flow_inputs(args, args.climate_scenario)
        issues.extend(flow_issues[:20])
        if len(flow_issues) > 20:
            issues.append(f"... and {len(flow_issues) - 20} additional future CaMa issues")
    return issues


def prepare_future_output_model(
    context: base.ModelContext,
    target: Path,
    output_root: Path,
    overwrite: bool,
    copy_mode: str,
) -> bool:
    """Create a future case and link/copy its shared VLM-adjusted terrain."""
    created = ORIGINAL_PREPARE_OUTPUT_MODEL(
        context, target, output_root, overwrite, copy_mode
    )
    if not created:
        return False
    args = ACTIVE_FUTURE_ARGS
    if args is None:
        raise RuntimeError("Future VLM configuration was not initialized")

    destination = target / "sfincs.dep"
    if args.no_vlm_adjustment:
        source = context.base_dir / "sfincs.dep"
        method = "original_base_model_terrain"
    else:
        source = args.vlm_adjusted_dep_root / context.model_id / "sfincs.dep"
        method = "vlm_adjusted_terrain"
    if not source.is_file():
        raise FileNotFoundError(f"Missing {method} for {context.model_id}: {source}")
    base_dep = context.base_dir / "sfincs.dep"
    if base_dep.is_file() and source.stat().st_size != base_dep.stat().st_size:
        raise ValueError(
            f"{context.model_id}: terrain byte count differs from the base grid "
            f"({source.stat().st_size} != {base_dep.stat().st_size})"
        )

    # sfincs.dep is immutable during a SFINCS run. Link it from the per-GTC VLM
    # source by default; never overwrite a hardlink in place.
    if destination.exists():
        destination.unlink()
    base.copy_runtime_file(source, destination, copy_mode)
    LOG.info(
        "Installed %s for %s: %s", method, context.model_id, destination
    )
    return True


def install_future_hooks() -> None:
    base.discover_events = discover_future_events
    base.read_track_points = read_future_track_points
    base.sample_monthly_t600 = sample_future_monthly_t600
    base.build_climada_track = future_climada_track
    base.load_or_compute_cama_quantile = load_or_compute_future_q50
    base.preflight = future_preflight
    # Prevent sfincs.dep from coming from the unadjusted base model;
    # prepare_future_output_model installs the selected shared terrain.
    base.MUTABLE_TOP_LEVEL.add("sfincs.dep")
    base.prepare_output_model = prepare_future_output_model


def resolve_future_output_root(args: argparse.Namespace, period_tag: str) -> None:
    """Resolve the future return-period-dependent output default."""
    if args.output_root is None:
        args.output_root = (
            SCRIPT_DIR
            / f"global_sfincs_{period_tag}_{base.flow_quantile_tag(args.cama_quantile)}_tcr"
        )


def run(args: argparse.Namespace) -> int:
    global ACTIVE_FUTURE_ARGS
    period_tag = base.return_period_tag(args.return_period)
    resolve_future_output_root(args, period_tag)
    scenarios = list(dict.fromkeys(args.scenarios))
    if args.no_diagnostic_cache and args.require_cama_diagnostic_cache:
        raise ValueError(
            "--no-diagnostic-cache and --require-diagnostic-cache cannot be used together"
        )
    if args.cama_quantile_mode == "tc_inlet_p90" and not np.isclose(
        args.cama_quantile, 0.9
    ):
        raise ValueError("tc_inlet_p90 requires --cama-quantile 0.9")
    if args.cama_quantile_mode not in {
        "tc_summary_p50", "tc_monthly_mean", "tc_inlet_p90"
    } and not np.isclose(args.cama_quantile, 0.5):
        raise ValueError("Future compound forcing is defined with q50; use --cama-quantile 0.5")
    if len(args.gcm_models) != 8 or set(args.gcm_models) != set(GCM_MODELS):
        raise ValueError(
            "Future river forcing requires exactly the configured eight GCMs: "
            + ", ".join(GCM_MODELS)
        )
    if (args.cama_start_year, args.cama_end_year) != (2061, 2100):
        LOG.warning(
            "Using nonstandard future CaMa period %s-%s",
            args.cama_start_year,
            args.cama_end_year,
        )
    install_future_hooks()
    status = 0
    for scenario in scenarios:
        scenario_args = copy.copy(args)
        scenario_args.climate_scenario = scenario
        scenario_args.adcirc_result_root = (
            args.future_adcirc_result_parent / f"{period_tag}_{scenario}"
        )
        scenario_args.adcirc_run_root = (
            args.future_adcirc_run_parent / f"{period_tag}_{scenario}"
        )
        scenario_args.output_root = args.output_root / scenario
        scenario_args.cama_root = args.future_cama_root
        scenario_args.t600_source = t600_path(
            args.future_t600_root,
            scenario,
            args.cama_start_year,
            args.cama_end_year,
        )
        scenario_args.vlm_adjusted_dep_root = args.vlm_adjusted_dep_root.resolve()
        scenario_args.cama_models_used = " ".join(args.gcm_models)
        if args.cama_quantile_mode in {"tc_inlet_p50", "tc_inlet_p90"}:
            scenario_args.cama_definition = (
                "Equal mean of eight GCM independent per-inlet TC-landfall-month-"
                f"probability-weighted P{int(round(args.cama_quantile * 100)):02d} values; "
                "source points on the same inlet section are equally split "
                f"({scenario})"
            )
        elif args.cama_quantile_mode == "tc_summary_p50":
            scenario_args.cama_definition = (
                "Equal mean of eight GCM tc_p50_flow_m3s values read directly from "
                "P1_plot_sfincs_inflow_diagnostics.py summary; inlet allocation by "
                f"q_reference_m3s fraction ({scenario})"
            )
        else:
            scenario_args.cama_definition = (
                f"{args.cama_start_year}-{args.cama_end_year} equal mean of eight GCM "
                "per-inlet sum of GTC TC-landfall month probability times diagnostic "
                f"monthly mean routed outflow ({scenario})"
            )
        LOG.info("=" * 72)
        LOG.info("Build future SFINCS scenario: %s", scenario)
        LOG.info("ADCIRC results: %s", scenario_args.adcirc_result_root)
        LOG.info("ADCIRC P4 runs: %s", scenario_args.adcirc_run_root)
        LOG.info("Future t600: %s", scenario_args.t600_source)
        LOG.info(
            "Terrain: %s",
            "original base model (--no-vlm-adjustment)"
            if scenario_args.no_vlm_adjustment
            else f"VLM adjusted from {scenario_args.vlm_adjusted_dep_root}",
        )
        LOG.info("Output: %s", scenario_args.output_root)
        ACTIVE_FUTURE_ARGS = scenario_args
        try:
            status = max(status, base.run(scenario_args))
        finally:
            ACTIVE_FUTURE_ARGS = None
    return status


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(message)s",
        force=True,
    )
    try:
        return run(args)
    except Exception:
        LOG.exception("Future SFINCS compound-forcing build failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
