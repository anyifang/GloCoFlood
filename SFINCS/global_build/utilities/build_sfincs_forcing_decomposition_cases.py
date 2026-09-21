#!/usr/bin/env python3
"""Build three leave-one-driver-out SFINCS experiment sets.

The existing ``global_sfincs_100yr_q50_tcr`` tree is the all-driver baseline.
For every baseline case this script writes only ``sfincs.inp`` into three new
case directories: no_river, no_rain and no_boundary.  All immutable model and
forcing inputs are referenced through portable relative paths, so the source
and decomposition roots must remain sibling directories after transfer.

Open-boundary cells (msk=2) cannot safely be left without a boundary series.
Therefore no_boundary keeps the original bnd locations and points bzsfile to a
shared, constant-zero series.  Shared zero series are stored once per boundary
point count/time window outside the individual case directories.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SOURCE = SCRIPT_DIR / "global_sfincs_100yr_q50_tcr"
DEFAULT_OUTPUT = SCRIPT_DIR / "global_sfincs_100yr_q50_tcr_decomposition"
SCENARIOS = ("historical", "ssp126", "ssp245", "ssp370")
VARIANTS = ("no_river", "no_rain", "no_boundary")

FILE_KEYS = {
    "depfile", "mskfile", "indexfile", "manningfile", "inifile", "rstfile",
    "bndfile", "bzsfile", "bzifile", "netbndbzsbzifile",
    "srcfile", "disfile", "netsrcdisfile", "obsfile", "crsfile",
    "thdfile", "thsfile", "weirfile", "drnfile", "sbgfile", "qtrfile",
    "spwfile", "spw2file", "wndfile", "netamuamvfile", "netampfile",
    "netamprfile", "precipfile", "amprfile", "qinffile", "netqinf_file",
    "psifile", "sigmafile", "ksfile", "f0file", "fcfile", "kdfile",
}
RIVER_KEYS = {"srcfile", "disfile", "netsrcdisfile"}
RAIN_KEYS = {"netamprfile", "precipfile", "amprfile"}
BOUNDARY_SERIES_KEYS = {"bzsfile", "netbndbzsbzifile"}


@dataclass(frozen=True)
class InputLine:
    raw: str
    key_original: str | None = None
    key: str | None = None
    value: str | None = None


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Create no-river, no-rain and no-boundary SFINCS cases."
    )
    p.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    p.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--scenarios", nargs="+", choices=SCENARIOS, default=list(SCENARIOS))
    p.add_argument("--limit-cases", type=int, default=0,
                   help="Limit baseline cases per scenario for QA; 0 means all.")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    return p


def clean_value(value: str) -> str:
    value = re.split(r"[#!]", value, maxsplit=1)[0].strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        value = value[1:-1]
    return value


def read_inp(path: Path) -> list[InputLine]:
    result: list[InputLine] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith(("#", "!")) or "=" not in raw:
            result.append(InputLine(raw=raw))
            continue
        left, right = raw.split("=", 1)
        original = left.strip()
        result.append(InputLine(raw=raw, key_original=original,
                                key=original.lower(), value=clean_value(right)))
    return result


def inp_values(lines: Iterable[InputLine]) -> dict[str, str]:
    return {line.key: line.value for line in lines if line.key is not None}


def parse_time(value: str) -> datetime:
    return datetime.strptime(value, "%Y%m%d %H%M%S")


def relative_reference(source: Path, case_dir: Path) -> str:
    return os.path.relpath(source, case_dir).replace(os.sep, "/")


def referenced_source_file(source_case: Path, value: str, key: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        raise ValueError(f"Absolute {key} is not portable: {path}")
    resolved = source_case / path
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise FileNotFoundError(f"Missing/nonempty {key}: {resolved}")
    return resolved


def boundary_point_count(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="ignore") as stream:
        for line in stream:
            text = line.strip()
            if text and not text.startswith(("#", "!")):
                count += 1
    if count <= 0:
        raise ValueError(f"No boundary points in {path}")
    return count


def zero_boundary_file(
    shared_root: Path,
    count: int,
    start_seconds: int,
    stop_seconds: int,
    dry_run: bool,
) -> Path:
    name = f"zero_bzs_np{count:06d}_t{start_seconds}_{stop_seconds}.bzs"
    path = shared_root / name
    if path.is_file() and path.stat().st_size > 0:
        return path
    if not dry_run:
        shared_root.mkdir(parents=True, exist_ok=True)
        zeros = " 0" * count
        path.write_text(
            f"{start_seconds}{zeros}\n{stop_seconds}{zeros}\n",
            encoding="ascii",
            newline="\n",
        )
    return path


def render_inp(
    lines: list[InputLine],
    source_case: Path,
    target_case: Path,
    variant: str,
    zero_bzs: Path | None,
) -> str:
    output: list[str] = []
    seen_zero_boundary = False
    for line in lines:
        if line.key is None:
            output.append(line.raw)
            continue
        key = line.key
        assert line.key_original is not None and line.value is not None
        if variant == "no_river" and key in RIVER_KEYS:
            continue
        if variant == "no_rain" and (key in RAIN_KEYS or key == "ampr_block"):
            continue
        if variant == "no_boundary" and key in BOUNDARY_SERIES_KEYS:
            if key != "bzsfile":
                raise ValueError(
                    "Only ASCII bndfile/bzsfile boundary forcing is supported "
                    f"for no_boundary, found {key} in {source_case}"
                )
            assert zero_bzs is not None
            value = relative_reference(zero_bzs, target_case)
            seen_zero_boundary = True
        elif key in FILE_KEYS and line.value not in ("", "-999"):
            value = relative_reference(
                referenced_source_file(source_case, line.value, key), target_case
            )
        else:
            value = line.value
        output.append(f"{line.key_original:<20} = {value}")
    if variant == "no_boundary" and not seen_zero_boundary:
        raise ValueError(f"Missing bzsfile in {source_case / 'sfincs.inp'}")
    return "\n".join(output).rstrip() + "\n"


def discover(source_root: Path, scenario: str, limit: int) -> list[Path]:
    root = source_root / scenario
    if not root.is_dir():
        raise FileNotFoundError(root)
    files = sorted(root.rglob("sfincs.inp"))
    if limit > 0:
        files = files[:limit]
    return files


def main() -> int:
    args = parser().parse_args()
    source_root = args.source_root.resolve()
    output_root = args.output_root.resolve()
    if source_root == output_root:
        raise ValueError("Source and output roots must differ")
    if args.limit_cases < 0:
        raise ValueError("--limit-cases must be non-negative")

    shared_root = output_root / "_shared_zero_boundary"
    boundary_count_cache: dict[tuple[str, int], int] = {}
    rows: list[dict[str, object]] = []
    written = skipped = 0

    for scenario in args.scenarios:
        cases = discover(source_root, scenario, args.limit_cases)
        print(f"{scenario}: {len(cases)} baseline cases")
        for number, inp_path in enumerate(cases, 1):
            source_case = inp_path.parent
            relative_case = source_case.relative_to(source_root / scenario)
            lines = read_inp(inp_path)
            values = inp_values(lines)
            required = {"depfile", "mskfile", "indexfile", "manningfile",
                        "bndfile", "bzsfile", "netamprfile", "tref", "tstart", "tstop"}
            missing = sorted(required - values.keys())
            if missing:
                raise ValueError(f"{inp_path} lacks required keys: {missing}")
            if abs(float(values.get("zsini", "0"))) > 1.0e-12:
                raise ValueError(
                    f"no_boundary is defined at 0 m, but zsini is not zero: {inp_path}"
                )
            tref = parse_time(values["tref"])
            tstart = parse_time(values["tstart"])
            tstop = parse_time(values["tstop"])
            start_seconds = int((tstart - tref).total_seconds())
            stop_seconds = int((tstop - tref).total_seconds())
            if stop_seconds <= start_seconds:
                raise ValueError(f"Invalid model time window: {inp_path}")

            bnd_path = referenced_source_file(source_case, values["bndfile"], "bndfile")
            cache_key = (relative_case.parts[-1], bnd_path.stat().st_size)
            if cache_key not in boundary_count_cache:
                boundary_count_cache[cache_key] = boundary_point_count(bnd_path)
            point_count = boundary_count_cache[cache_key]
            zero_bzs = zero_boundary_file(
                shared_root, point_count, start_seconds, stop_seconds, args.dry_run
            )

            river_present = all(k in values for k in ("srcfile", "disfile"))
            for variant in VARIANTS:
                target_case = output_root / scenario / variant / relative_case
                target_inp = target_case / "sfincs.inp"
                text = render_inp(
                    lines, source_case, target_case, variant,
                    zero_bzs if variant == "no_boundary" else None,
                )
                status = "would_write" if args.dry_run else "written"
                if target_inp.exists() and not args.overwrite:
                    current = target_inp.read_text(encoding="utf-8", errors="ignore")
                    if current != text:
                        raise FileExistsError(
                            f"Existing input differs; use --overwrite: {target_inp}"
                        )
                    status = "unchanged"
                    skipped += 1
                elif not args.dry_run:
                    target_case.mkdir(parents=True, exist_ok=True)
                    target_inp.write_text(text, encoding="utf-8", newline="\n")
                    written += 1
                rows.append({
                    "scenario": scenario,
                    "variant": variant,
                    "river_enabled": int(variant != "no_river" and river_present),
                    "rain_enabled": int(variant != "no_rain"),
                    "boundary_enabled": int(variant != "no_boundary"),
                    "river_present_in_baseline": int(river_present),
                    "source_case": source_case.as_posix(),
                    "source_full_result": (source_case / "sfincs_map.nc").as_posix(),
                    "decomposition_case": target_case.as_posix(),
                    "sfincs_inp": target_inp.as_posix(),
                    "status": status,
                })
            if number % 250 == 0:
                print(f"  {number}/{len(cases)} baseline cases processed")

    if not args.dry_run:
        output_root.mkdir(parents=True, exist_ok=True)
        manifest = output_root / "decomposition_case_manifest.csv"
        with manifest.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        readme = output_root / "README.txt"
        readme.write_text(
            "Three leave-one-driver-out experiment sets.\n"
            "Full baseline: ../global_sfincs_100yr_q50_tcr\n"
            "no_river: src/dis forcing omitted.\n"
            "no_rain: spatial precipitation forcing omitted.\n"
            "no_boundary: original boundary locations forced at constant 0 m.\n"
            "Wind and pressure are retained in every experiment.\n"
            "If river_present_in_baseline=0 in the manifest, no_river is identical\n"
            "to the full baseline and may be skipped during submission.\n"
            "Each case directory initially contains only sfincs.inp; immutable inputs\n"
            "are read through relative paths. Keep source/output roots as siblings.\n",
            encoding="utf-8",
        )
        print(f"Manifest: {manifest}")
    print(f"Generated rows={len(rows)}; written={written}; unchanged={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
