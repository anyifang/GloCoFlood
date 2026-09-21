#!/usr/bin/env python3
"""Run historical, SSP126, SSP245 and SSP370 SFINCS builders in parallel.

Eight independent child processes are started by default: every scenario is
split into two disjoint ADCIRC-block shards. Each process is limited to one
numerical-library thread and pinned to a different available logical CPU, so
the batch uses approximately eight CPU cores without nested BLAS/OpenMP
oversubscription. Each builder writes to its own scenario folder:

    <output-root>/historical
    <output-root>/ssp126
    <output-root>/ssp245
    <output-root>/ssp370

Existing complete cases are skipped by the underlying builders unless
``--overwrite`` is explicitly supplied.

conda activate sfincs_tcr
python run_four_compound_builds_parallel.py --return-period 100 --cama-quantile 0.50 --continue-on-error --overwrite
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
HISTORICAL_BUILDER = SCRIPT_DIR / "P2_build_sfincs_historical_return_period_compound_forcing.py"
FUTURE_BUILDER = SCRIPT_DIR / "P3_build_sfincs_future_return_period_compound_forcing.py"
SCENARIOS = ("historical", "ssp126", "ssp245", "ssp370")
DEFAULT_GLOBAL_FLOOD_PROJECT_DIR = Path(
    os.environ.get(
        "GLOCOFLOOD_PROJECT_DIR",
        os.environ.get("GLOBAL_FLOOD_PROJECT_DIR", SCRIPT_DIR.parent.parent),
    )
)
DEFAULT_MEMBERSHIP = (
    DEFAULT_GLOBAL_FLOOD_PROJECT_DIR
    / "ADCIRC"
    / "global_build"
    / "catalog"
    / "global_tc_adcirc_block_members.csv"
)
DEFAULT_ADCIRC_MESH_ROOT = (
    DEFAULT_GLOBAL_FLOOD_PROJECT_DIR
    / "external"
    / "adcirc"
    / "adcirc_fort14_meshes"
)
DEFAULT_BOUNDARY_MAP_CACHE_DIR = Path(
    os.environ.get(
        "SFINCS_BOUNDARY_MAP_CACHE",
        SCRIPT_DIR / "cache" / "adcirc_sfincs_boundary_maps",
    )
)
CLIMADA_TCR_DATA_DIR = Path(
    os.environ.get(
        "CLIMADA_TCR_DATA_DIR",
        DEFAULT_GLOBAL_FLOOD_PROJECT_DIR
        / "external/climada/data/hazard/tc_rainfield",
    )
)
DEFAULT_TCR_ELEVATION_TIF = (
    CLIMADA_TCR_DATA_DIR
    / "topography_land_360as/v1/topography_land_360as.tif"
)
DEFAULT_TCR_DRAG_TIF = CLIMADA_TCR_DATA_DIR / "c_drag_500/v1/c_drag_500.tif"


@dataclass
class RunningBuild:
    name: str
    command: list[str]
    process: subprocess.Popen[str]
    log_path: Path
    log_handle: object
    reader: threading.Thread
    cpu: int | None


@dataclass(frozen=True)
class BuildJob:
    name: str
    scenario: str
    blocks: tuple[str, ...]
    command: list[str]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run historical and three future SFINCS compound-forcing builders "
            "as block-sharded one-core processes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--return-period", type=int, default=100, metavar="YEARS")
    parser.add_argument(
        "--cama-quantile",
        type=float,
        default=0.50,
        help="CaMa-Flood non-exceedance probability; 0.50 or 0.90.",
    )
    parser.add_argument(
        "--cama-quantile-mode",
        choices=("tc_inlet_p50", "tc_inlet_p90"),
        default=None,
        help="If omitted, inferred from --cama-quantile.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help=(
            "Common output root. Default: "
            "global_sfincs_<return-period>yr_q<quantile>_tcr beside this script."
        ),
    )
    parser.add_argument("--blocks", help="Comma-separated ADCIRC block IDs.")
    parser.add_argument("--models", help="Comma-separated GTC model IDs.")
    parser.add_argument("--events", help="Comma-separated exact event directory names.")
    parser.add_argument("--limit-events", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--workers",
        type=int,
        choices=(4, 8),
        default=8,
        help=(
            "Parallel one-core builder processes. Four uses one process per scenario; "
            "eight splits every scenario into two ADCIRC-block shards."
        ),
    )
    parser.add_argument(
        "--copy-mode",
        choices=("hardlink", "copy"),
        default="hardlink",
        help="Reuse immutable/shared files by hard link when possible.",
    )
    parser.add_argument(
        "--membership-csv", type=Path, default=DEFAULT_MEMBERSHIP
    )
    parser.add_argument(
        "--adcirc-mesh-root", type=Path, default=DEFAULT_ADCIRC_MESH_ROOT
    )
    parser.add_argument(
        "--boundary-map-cache-dir",
        type=Path,
        default=DEFAULT_BOUNDARY_MAP_CACHE_DIR,
        help="Persistent ADCIRC-to-SFINCS boundary-node mapping cache.",
    )
    parser.add_argument(
        "--rebuild-boundary-map-cache",
        action="store_true",
        help="Force the one-time boundary-map preparation step to rebuild its cache.",
    )
    parser.add_argument(
        "--tcr-elevation-tif",
        type=Path,
        default=DEFAULT_TCR_ELEVATION_TIF,
        help=(
            "Local TCR topography raster passed to all four builders. Providing "
            "this explicitly prevents CLIMADA from querying its online API."
        ),
    )
    parser.add_argument(
        "--tcr-drag-tif",
        type=Path,
        default=DEFAULT_TCR_DRAG_TIF,
        help=(
            "Local TCR drag-coefficient raster passed to all four builders. "
            "Providing this explicitly prevents CLIMADA API access."
        ),
    )
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING"), default="INFO"
    )
    parser.add_argument(
        "--no-pin-cores",
        action="store_true",
        help=(
            "Do not pin the four child processes to distinct CPUs. Numerical "
            "thread counts are still limited to one per process."
        ),
    )
    parser.add_argument(
        "--python",
        type=Path,
        default=None,
        help=(
            "Python executable used to launch each builder. If omitted, the "
            "sfincs_tcr Conda environment is located automatically."
        ),
    )
    return parser


def _python_executable(prefix: Path) -> Path:
    """Return the platform-specific Python executable below an environment."""
    return prefix / ("python.exe" if os.name == "nt" else "bin/python")


def _conda_environment_prefixes() -> list[Path]:
    """Read Conda's environment registry without requiring shell activation."""
    conda = shutil.which("conda")
    if not conda:
        return []
    try:
        result = subprocess.run(
            [conda, "info", "--envs", "--json"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        payload = json.loads(result.stdout)
        return [Path(value) for value in payload.get("envs", [])]
    except (OSError, subprocess.SubprocessError, ValueError, TypeError):
        return []


def resolve_builder_python(requested: Path | None) -> Path:
    """Select the requested Python or automatically find ``sfincs_tcr``."""
    if requested is not None:
        candidate = requested.expanduser()
        if not candidate.is_absolute():
            candidate = (Path.cwd() / candidate).resolve()
        if not candidate.is_file():
            raise FileNotFoundError(f"Python executable does not exist: {candidate}")
        return candidate

    candidates: list[Path] = []
    active_prefix = os.environ.get("CONDA_PREFIX")
    if active_prefix and Path(active_prefix).name.lower() == "sfincs_tcr":
        candidates.append(_python_executable(Path(active_prefix)))
    if Path(sys.prefix).name.lower() == "sfincs_tcr":
        candidates.append(Path(sys.executable))
    candidates.extend(
        _python_executable(prefix)
        for prefix in _conda_environment_prefixes()
        if prefix.name.lower() == "sfincs_tcr"
    )
    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate).lower()
        if key not in seen and candidate.is_file():
            return candidate.resolve()
        seen.add(key)
    raise RuntimeError(
        "Cannot locate the sfincs_tcr Conda environment. Activate it first, "
        "or pass --python <path-to-sfincs_tcr-python.exe>."
    )


def validate_tcr_python(python: Path, env: dict[str, str]) -> None:
    """Fail once, before launching four jobs, if TCR dependencies are absent."""
    probe = subprocess.run(
        [
            str(python),
            "-c",
            "import climada, climada_petals; print('climada and climada_petals OK')",
        ],
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    if probe.returncode != 0:
        details = (probe.stderr or probe.stdout).strip()
        raise RuntimeError(
            f"Selected Python cannot import the TCR dependencies: {python}\n{details}"
        )


def configure_conda_data_paths(env: dict[str, str], python: Path) -> None:
    """Supply data directories normally set by Conda activation on Windows."""
    prefix = python.parent
    env["CONDA_PREFIX"] = str(prefix)
    env["CONDA_DEFAULT_ENV"] = prefix.name
    if os.name == "nt":
        path_entries = [
            prefix,
            prefix / "Scripts",
            prefix / "Library/bin",
            prefix / "Library/usr/bin",
            prefix / "Library/mingw-w64/bin",
        ]
        current_path = env.get("PATH", "")
        env["PATH"] = os.pathsep.join(
            [str(path) for path in path_entries if path.is_dir()] + [current_path]
        )
    for variable, relative in (
        ("GDAL_DATA", Path("Library/share/gdal")),
        ("PROJ_DATA", Path("Library/share/proj")),
        ("PROJ_LIB", Path("Library/share/proj")),
    ):
        candidate = prefix / relative
        if candidate.is_dir():
            env[variable] = str(candidate)


def quantile_mode(quantile: float, requested: str | None) -> str:
    if requested is not None:
        expected = 0.50 if requested == "tc_inlet_p50" else 0.90
        if abs(float(quantile) - expected) > 1.0e-9:
            raise ValueError(
                f"{requested} requires --cama-quantile {expected:.2f}, "
                f"not {quantile}"
            )
        return requested
    if abs(float(quantile) - 0.50) <= 1.0e-9:
        return "tc_inlet_p50"
    if abs(float(quantile) - 0.90) <= 1.0e-9:
        return "tc_inlet_p90"
    raise ValueError(
        "Automatic mode selection supports --cama-quantile 0.50 or 0.90; "
        "specify a supported quantile."
    )


def output_root(args: argparse.Namespace) -> Path:
    if args.output_root is not None:
        path = args.output_root.expanduser()
        return path.resolve() if path.is_absolute() else (Path.cwd() / path).resolve()
    qtag = f"q{int(round(float(args.cama_quantile) * 100)):02d}"
    return (
        SCRIPT_DIR
        / f"global_sfincs_{int(args.return_period)}yr_{qtag}_tcr"
    ).resolve()


def resolve_required_local_file(path: Path, option: str) -> Path:
    """Resolve and validate a local file before any child builder is started."""
    candidate = path.expanduser()
    if not candidate.is_absolute():
        candidate = (Path.cwd() / candidate).resolve()
    else:
        candidate = candidate.resolve()
    if not candidate.is_file():
        raise FileNotFoundError(f"{option} file does not exist: {candidate}")
    if candidate.stat().st_size <= 0:
        raise ValueError(f"{option} file is empty: {candidate}")
    return candidate


def common_builder_args(
    args: argparse.Namespace, root: Path, mode: str
) -> list[str]:
    result = [
        "--return-period",
        str(int(args.return_period)),
        "--cama-quantile-mode",
        mode,
        "--cama-quantile",
        f"{float(args.cama_quantile):.12g}",
        "--output-root",
        str(root),
        "--c15-rmax-out-of-range",
        "rescale",
        "--copy-mode",
        args.copy_mode,
        "--membership-csv",
        str(args.membership_csv),
        "--adcirc-mesh-root",
        str(args.adcirc_mesh_root),
        "--boundary-map-cache-dir",
        str(args.boundary_map_cache_dir),
        "--log-level",
        args.log_level,
        "--tcr-elevation-tif",
        str(args.tcr_elevation_tif),
        "--tcr-drag-tif",
        str(args.tcr_drag_tif),
    ]
    for option, value in (
        ("--models", args.models),
        ("--events", args.events),
    ):
        if value:
            result.extend([option, str(value)])
    if args.limit_events is not None:
        result.extend(["--limit-events", str(int(args.limit_events))])
    if args.overwrite:
        result.append("--overwrite")
    if args.continue_on_error:
        result.append("--continue-on-error")
    if args.dry_run:
        result.append("--dry-run")
    return result


def _csv_values(value: str | None) -> set[str] | None:
    if not value:
        return None
    values = {part.strip() for part in value.split(",") if part.strip()}
    return values or None


def selected_blocks(args: argparse.Namespace) -> list[str]:
    """Return blocks left after the runner's block/model filters."""
    requested_blocks = _csv_values(args.blocks)
    requested_models = _csv_values(args.models)
    found_blocks: set[str] = set()
    with args.membership_csv.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"block_id", "model_domain_id"}.issubset(reader.fieldnames or []):
            raise ValueError(
                f"{args.membership_csv} must contain block_id and model_domain_id"
            )
        for row in reader:
            block_id = str(row["block_id"]).strip()
            model_id = str(row["model_domain_id"]).strip()
            if requested_blocks is not None and block_id not in requested_blocks:
                continue
            if requested_models is not None and model_id not in requested_models:
                continue
            found_blocks.add(block_id)
    if requested_blocks is not None:
        missing = sorted(requested_blocks - found_blocks)
        if missing:
            raise ValueError(
                "Selected block(s) have no matching GTC models: " + ", ".join(missing)
            )
    if not found_blocks:
        raise ValueError("No ADCIRC blocks remain after applying --blocks/--models")
    return sorted(found_blocks)


def split_blocks(blocks: Sequence[str], parts: int) -> list[tuple[str, ...]]:
    shards: list[list[str]] = [[] for _ in range(parts)]
    for index, block_id in enumerate(blocks):
        shards[index % parts].append(block_id)
    return [tuple(shard) for shard in shards if shard]


def make_jobs(
    builder_python: Path,
    common: Sequence[str],
    blocks: Sequence[str],
    workers: int,
) -> list[BuildJob]:
    shards_per_scenario = workers // len(SCENARIOS)
    shards = split_blocks(blocks, shards_per_scenario)
    jobs: list[BuildJob] = []
    for scenario in SCENARIOS:
        for shard_index, shard in enumerate(shards, start=1):
            name = (
                scenario
                if len(shards) == 1
                else f"{scenario}_part{shard_index}"
            )
            if scenario == "historical":
                command = [str(builder_python), str(HISTORICAL_BUILDER), *common]
            else:
                command = [
                    str(builder_python),
                    str(FUTURE_BUILDER),
                    "--scenarios",
                    scenario,
                    *common,
                ]
            command.extend(["--blocks", ",".join(shard)])
            jobs.append(BuildJob(name, scenario, shard, command))
    return jobs


def prepare_boundary_maps(
    builder_python: Path,
    common: Sequence[str],
    blocks: Sequence[str],
    rebuild: bool,
    env: dict[str, str],
    log_path: Path,
) -> None:
    command = [
        str(builder_python),
        str(HISTORICAL_BUILDER),
        *common,
        "--blocks",
        ",".join(blocks),
        "--prepare-boundary-maps-only",
    ]
    if rebuild:
        command.append("--rebuild-boundary-map-cache")
    print("Preparing persistent ADCIRC-to-SFINCS boundary maps first...", flush=True)
    with log_path.open("w", encoding="utf-8", buffering=1) as log_handle:
        log_handle.write("COMMAND: " + subprocess.list2cmdline(command) + "\n")
        process = subprocess.Popen(
            command,
            cwd=SCRIPT_DIR,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            log_handle.write(line)
            print(f"[boundary-map] {line}", end="", flush=True)
        code = int(process.wait())
    if code != 0:
        raise RuntimeError(
            f"Boundary-map preparation failed with exit code {code}; log={log_path}"
        )


def prepare_flow_caches(
    builder_python: Path,
    common: Sequence[str],
    blocks: Sequence[str],
    env: dict[str, str],
    log_dir: Path,
) -> None:
    """Build each scenario's complete flow cache before block sharding."""
    commands: list[tuple[str, list[str]]] = []
    for scenario in SCENARIOS:
        if scenario == "historical":
            command = [str(builder_python), str(HISTORICAL_BUILDER), *common]
        else:
            command = [
                str(builder_python),
                str(FUTURE_BUILDER),
                "--scenarios",
                scenario,
                *common,
            ]
        command.extend(
            ["--blocks", ",".join(blocks), "--prepare-flow-cache-only"]
        )
        commands.append((scenario, command))
    print("Preparing four complete scenario flow caches...", flush=True)
    processes: list[tuple[str, subprocess.Popen[str], object, Path]] = []
    try:
        for scenario, command in commands:
            log_path = log_dir / f"flow_cache_{scenario}.log"
            handle = log_path.open("w", encoding="utf-8", buffering=1)
            handle.write("COMMAND: " + subprocess.list2cmdline(command) + "\n")
            process = subprocess.Popen(
                command,
                cwd=SCRIPT_DIR,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            processes.append((scenario, process, handle, log_path))
        failures: list[str] = []
        for scenario, process, _, log_path in processes:
            code = int(process.wait())
            state = "ready" if code == 0 else f"failed(exit={code})"
            print(f"  flow cache {scenario}: {state}; log={log_path}", flush=True)
            if code != 0:
                failures.append(f"{scenario}: exit={code}, log={log_path}")
        if failures:
            raise RuntimeError("Flow-cache preparation failed: " + "; ".join(failures))
    except KeyboardInterrupt:
        for _, process, _, _ in processes:
            if process.poll() is None:
                process.terminate()
        raise
    finally:
        for _, _, handle, _ in processes:
            handle.close()  # type: ignore[attr-defined]


def available_cpus() -> list[int]:
    try:
        import psutil  # type: ignore

        cpus = list(psutil.Process().cpu_affinity())
        if cpus:
            return [int(cpu) for cpu in cpus]
    except Exception:
        pass
    count = os.cpu_count() or 1
    return list(range(int(count)))


def pin_process(process: subprocess.Popen[str], cpu: int) -> bool:
    try:
        import psutil  # type: ignore

        psutil.Process(process.pid).cpu_affinity([int(cpu)])
        return True
    except Exception:
        pass
    if hasattr(os, "sched_setaffinity"):
        try:
            os.sched_setaffinity(process.pid, {int(cpu)})  # type: ignore[attr-defined]
            return True
        except Exception:
            pass
    return False


def stream_output(
    name: str,
    process: subprocess.Popen[str],
    log_handle: object,
) -> None:
    assert process.stdout is not None
    for line in process.stdout:
        text = line.rstrip("\r\n")
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        rendered = f"{timestamp} [{name}] {text}\n"
        log_handle.write(rendered)  # type: ignore[attr-defined]
        log_handle.flush()  # type: ignore[attr-defined]
        print(rendered, end="", flush=True)


def terminate_all(builds: Sequence[RunningBuild]) -> None:
    for build in builds:
        if build.process.poll() is None:
            build.process.terminate()
    deadline = time.monotonic() + 10.0
    for build in builds:
        remaining = max(0.0, deadline - time.monotonic())
        if build.process.poll() is None:
            try:
                build.process.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                build.process.kill()


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.return_period <= 0:
        raise ValueError("--return-period must be positive")
    builder_python = resolve_builder_python(args.python)
    args.tcr_elevation_tif = resolve_required_local_file(
        args.tcr_elevation_tif, "--tcr-elevation-tif"
    )
    args.tcr_drag_tif = resolve_required_local_file(
        args.tcr_drag_tif, "--tcr-drag-tif"
    )
    mode = quantile_mode(args.cama_quantile, args.cama_quantile_mode)
    for attribute in ("membership_csv", "adcirc_mesh_root", "boundary_map_cache_dir"):
        value = getattr(args, attribute).expanduser()
        setattr(
            args,
            attribute,
            value.resolve() if value.is_absolute() else (Path.cwd() / value).resolve(),
        )
    if not args.membership_csv.is_file():
        raise FileNotFoundError(f"Membership CSV does not exist: {args.membership_csv}")
    if not args.adcirc_mesh_root.is_dir():
        raise FileNotFoundError(f"ADCIRC mesh root does not exist: {args.adcirc_mesh_root}")
    block_ids = selected_blocks(args)
    root = output_root(args)
    root.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = root / "_parallel_build_logs" / timestamp
    log_dir.mkdir(parents=True, exist_ok=False)

    common = common_builder_args(args, root, mode)
    jobs = make_jobs(builder_python, common, block_ids, args.workers)

    child_env = os.environ.copy()
    configure_conda_data_paths(child_env, builder_python)
    for variable in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "GDAL_NUM_THREADS",
    ):
        child_env[variable] = "1"
    child_env["PYTHONUNBUFFERED"] = "1"
    validate_tcr_python(builder_python, child_env)

    if not args.dry_run:
        prepare_boundary_maps(
            builder_python,
            common,
            block_ids,
            args.rebuild_boundary_map_cache,
            child_env,
            log_dir / "boundary_map_precompute.log",
        )
        prepare_flow_caches(builder_python, common, block_ids, child_env, log_dir)

    cpus = available_cpus()
    selected_cpus = [cpus[index % len(cpus)] for index in range(len(jobs))]
    print(f"Output root: {root}", flush=True)
    print(f"Logs: {log_dir}", flush=True)
    print(f"Builder Python: {builder_python}", flush=True)
    print(f"TCR elevation: {args.tcr_elevation_tif}", flush=True)
    print(f"TCR drag: {args.tcr_drag_tif}", flush=True)
    print("C15 out-of-range Rmax: rescale", flush=True)
    print(f"Boundary map cache: {args.boundary_map_cache_dir}", flush=True)
    print(
        f"Starting {len(jobs)} one-core builders: "
        + ", ".join(
            f"{job.name}@CPU{selected_cpus[index]}[{len(job.blocks)} blocks]"
            for index, job in enumerate(jobs)
        ),
        flush=True,
    )

    builds: list[RunningBuild] = []
    try:
        for index, job in enumerate(jobs):
            name = job.name
            log_path = log_dir / f"{name}.log"
            log_handle = log_path.open("w", encoding="utf-8", buffering=1)
            command = job.command
            log_handle.write("COMMAND: " + subprocess.list2cmdline(command) + "\n")
            process = subprocess.Popen(
                command,
                cwd=SCRIPT_DIR,
                env=child_env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
            )
            cpu = None
            if not args.no_pin_cores and pin_process(process, selected_cpus[index]):
                cpu = selected_cpus[index]
            reader = threading.Thread(
                target=stream_output,
                args=(name, process, log_handle),
                name=f"log-{name}",
                daemon=True,
            )
            reader.start()
            builds.append(
                RunningBuild(
                    name=name,
                    command=command,
                    process=process,
                    log_path=log_path,
                    log_handle=log_handle,
                    reader=reader,
                    cpu=cpu,
                )
            )

        exit_codes: dict[str, int] = {}
        for build in builds:
            exit_codes[build.name] = int(build.process.wait())
        for build in builds:
            build.reader.join(timeout=5.0)

        print("\nParallel build summary:", flush=True)
        for build in builds:
            code = exit_codes[build.name]
            state = "SUCCESS" if code == 0 else "FAILED"
            affinity = f"CPU{build.cpu}" if build.cpu is not None else "not pinned"
            print(
                f"  {build.name:10s} {state:7s} exit={code} "
                f"({affinity}) log={build.log_path}",
                flush=True,
            )
        return 0 if all(code == 0 for code in exit_codes.values()) else 1
    except KeyboardInterrupt:
        print("\nInterrupted; terminating child builders...", file=sys.stderr)
        terminate_all(builds)
        return 130
    finally:
        for build in builds:
            try:
                build.log_handle.close()  # type: ignore[attr-defined]
            except Exception:
                pass


if __name__ == "__main__":
    raise SystemExit(main())
