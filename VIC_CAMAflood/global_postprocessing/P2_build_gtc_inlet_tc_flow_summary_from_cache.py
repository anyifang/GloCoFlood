#!/usr/bin/env python3
"""Build TC-month-weighted flow statistics for every independent CaMa inlet."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_ROOT = Path(
    os.environ.get(
        "GLOCOFLOOD_CAMA_DIAGNOSTIC_ROOT",
        Path(__file__).resolve().parent / "output" / "gtc_cama_inflow_figures",
    )
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--cache-key", default="17b34a8f6c68")
    return parser.parse_args()


def parse_case(path: Path, cache_key: str) -> tuple[str, str, str]:
    suffix = f"_{cache_key}.npz"
    if not path.name.startswith("cell_daily_") or not path.name.endswith(suffix):
        raise ValueError(f"Unexpected cache filename: {path.name}")
    case = path.name[len("cell_daily_") : -len(suffix)]
    if case == "ERA5":
        return case, "historical", "ERA5"
    match = re.fullmatch(r"cmip6_(.+)_(ssp\d+)_2061_2100", case)
    if match is None:
        raise ValueError(f"Cannot parse future case: {case}")
    return case, match.group(2), match.group(1)


def weighted_inlet_statistics(
    dates: np.ndarray,
    values: np.ndarray,
    month_weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return mean, TC-month-weighted P50 and P90 for each values column."""
    values = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(values)
    counts = finite.sum(axis=0)
    means = np.divide(
        np.nansum(values, axis=0),
        counts,
        out=np.full(values.shape[1], np.nan),
        where=counts > 0,
    )

    months = pd.DatetimeIndex(dates).month.to_numpy(dtype=np.int8)
    sample_weights = np.zeros(values.shape, dtype=np.float64)
    for month in range(1, 13):
        row_mask = months == month
        if not row_mask.any():
            continue
        valid_month = finite[row_mask, :]
        valid_counts = valid_month.sum(axis=0)
        scale = np.divide(
            month_weights[month - 1, :],
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

    p50 = np.full(values.shape[1], np.nan)
    p90 = np.full(values.shape[1], np.nan)
    for column in range(values.shape[1]):
        keep = np.isfinite(sorted_values[:, column]) & (sorted_weights[:, column] > 0)
        if not keep.any() or totals[column] <= 0:
            continue
        cumulative_column = cumulative[keep, column]
        values_column = sorted_values[keep, column]
        p50[column] = np.interp(0.50 * totals[column], cumulative_column, values_column)
        p90[column] = np.interp(0.90 * totals[column], cumulative_column, values_column)
    return means, p50, p90


def load_month_weights(weights_path: Path, inlets: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    weights = pd.read_csv(weights_path)
    required = {
        "gtc_id", "tc_weight_source", "tc_event_count", "month", "month_weight"
    }
    missing = required.difference(weights.columns)
    if missing:
        raise ValueError(f"Month-weight table lacks columns: {sorted(missing)}")

    counts = weights.groupby("gtc_id").size()
    if not (counts == 12).all():
        raise ValueError("Every GTC must have exactly 12 month-weight rows")
    sums = weights.groupby("gtc_id")["month_weight"].sum()
    if not np.allclose(sums.to_numpy(), 1.0, atol=1e-6):
        raise ValueError("TC month weights do not sum to one")

    pivot = weights.pivot(index="gtc_id", columns="month", values="month_weight")
    pivot = pivot.reindex(columns=range(1, 13))
    if pivot.isna().any().any():
        raise ValueError("TC month-weight matrix contains missing values")
    matrix = pivot.loc[inlets["gtc_id"]].to_numpy(dtype=np.float64).T

    metadata = (
        weights[["gtc_id", "tc_weight_source", "tc_event_count"]]
        .drop_duplicates("gtc_id")
        .set_index("gtc_id")
    )
    return matrix, metadata


def main() -> None:
    args = parse_args()
    table_dir = args.root / "tables"
    cache_dir = args.root / "cache"
    inlet_path = table_dir / "gtc_cama_inlet_cells.csv"
    weight_path = table_dir / "gtc_tc_month_weights.csv"
    output_path = table_dir / "gtc_flow_summary_by_inlet_cell.csv"

    inlets = pd.read_csv(inlet_path)
    inlet_required = {
        "gtc_id", "inlet_order", "inlet_id", "cell_id", "cama_row", "cama_col",
        "cama_lon", "cama_lat", "uparea_km2", "q_reference_m3s", "basin_id",
        "basin_label",
    }
    missing = inlet_required.difference(inlets.columns)
    if missing:
        raise ValueError(f"Inlet table lacks columns: {sorted(missing)}")
    inlet_keys = ["gtc_id", "inlet_id", "cell_id"]
    if inlets.duplicated(inlet_keys).any():
        raise ValueError("Independent inlet keys are not unique")

    month_weights, weight_metadata = load_month_weights(weight_path, inlets)
    cache_files = sorted(cache_dir.glob(f"cell_daily_*_{args.cache_key}.npz"))
    if len(cache_files) != 25:
        raise ValueError(f"Expected 25 cache files, found {len(cache_files)}")

    output_rows: list[pd.DataFrame] = []
    for file_index, cache_path in enumerate(cache_files, start=1):
        case, scenario, model = parse_case(cache_path, args.cache_key)
        print(f"[{file_index:02d}/{len(cache_files):02d}] {case}", flush=True)
        with np.load(cache_path, allow_pickle=False) as cache:
            dates = cache["dates"]
            values = cache["values"]
            rows = cache["rows"].astype(np.int64)
            cols = cache["cols"].astype(np.int64)

        cell_to_column = {(int(row), int(col)): idx for idx, (row, col) in enumerate(zip(rows, cols))}
        requested_cells = list(zip(inlets["cama_row"].astype(int), inlets["cama_col"].astype(int)))
        missing_cells = sorted(set(requested_cells).difference(cell_to_column))
        if missing_cells:
            raise ValueError(f"{cache_path.name} lacks {len(missing_cells)} requested inlet cells")
        inlet_columns = np.array([cell_to_column[cell] for cell in requested_cells], dtype=np.int64)
        inlet_values = values[:, inlet_columns]
        mean_flow, p50_flow, p90_flow = weighted_inlet_statistics(
            dates, inlet_values, month_weights
        )

        frame = inlets[
            [
                "gtc_id", "inlet_order", "inlet_id", "cell_id", "cama_row", "cama_col",
                "cama_lon", "cama_lat", "uparea_km2", "q_reference_m3s", "basin_id",
                "basin_label",
            ]
        ].copy()
        frame.insert(1, "case", case)
        frame.insert(2, "scenario", scenario)
        frame.insert(3, "model", model)
        frame["tc_weight_source"] = frame["gtc_id"].map(weight_metadata["tc_weight_source"])
        frame["tc_event_count"] = frame["gtc_id"].map(weight_metadata["tc_event_count"])
        frame["mean_flow_m3s"] = mean_flow
        frame["tc_p50_flow_m3s"] = p50_flow
        frame["tc_p90_flow_m3s"] = p90_flow
        output_rows.append(frame)

    result = pd.concat(output_rows, ignore_index=True)
    expected_rows = len(inlets) * len(cache_files)
    if len(result) != expected_rows:
        raise ValueError(f"Expected {expected_rows} rows, built {len(result)}")
    historical = result[result["scenario"] == "historical"]
    if len(historical) != len(inlets) or historical.duplicated(inlet_keys).any():
        raise ValueError("Historical output does not contain one row per independent inlet")
    for scenario in ("ssp126", "ssp245", "ssp370"):
        future = result[result["scenario"] == scenario]
        model_counts = future.groupby(inlet_keys)["model"].nunique()
        if len(future) != len(inlets) * 8 or not (model_counts == 8).all():
            raise ValueError(f"{scenario} does not contain eight models per independent inlet")
    if not np.isfinite(result[["tc_p50_flow_m3s", "tc_p90_flow_m3s"]]).all().all():
        raise ValueError("Inlet flow output contains non-finite weighted quantiles")

    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    result.to_csv(temporary_path, index=False, encoding="utf-8-sig")
    temporary_path.replace(output_path)
    print(f"Wrote {len(result)} rows for {len(inlets)} independent inlets: {output_path}")


if __name__ == "__main__":
    main()
