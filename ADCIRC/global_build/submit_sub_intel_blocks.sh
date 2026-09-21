#!/usr/bin/env bash
# Batch submit ADCIRC block runs from relative ADC_* block directories.
#
# Examples:
#   bash submit_sub_intel_blocks.sh
#   bash submit_sub_intel_blocks.sh --dry-run
#   bash submit_sub_intel_blocks.sh --pattern 'ADC_NA_*'
#   bash submit_sub_intel_blocks.sh --pattern 'ADC_NA_01' --force

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ -d "adcirc_fort14_meshes" ]]; then
    MESH_ROOT="adcirc_fort14_meshes"
else
    MESH_ROOT="."
fi
PATTERN="ADC_*"
DRY_RUN=0
FORCE=0
MAX_SUBMIT=0
SLEEP_SECONDS=0.2

usage() {
    cat <<EOF
Usage: bash $(basename "$0") [options]

Options:
  --mesh-root DIR      Mesh block root. Default: ./adcirc_fort14_meshes if it exists, otherwise current script directory.
  --pattern GLOB      Block directory name glob. Default: ADC_*
  --dry-run           Print jobs that would be submitted, but do not call sbatch.
  --force             Submit even when slurm-*.out or fort.63 already exists.
  --max N             Submit at most N jobs. Default: 0 means no limit.
  --sleep SEC         Sleep seconds between sbatch calls. Default: 0.2
  -h, --help          Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mesh-root)
            MESH_ROOT="${2%/}"
            shift 2
            ;;
        --pattern)
            PATTERN="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --max)
            MAX_SUBMIT="$2"
            shift 2
            ;;
        --sleep)
            SLEEP_SECONDS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ ! -d "$MESH_ROOT" ]]; then
    echo "ERROR: mesh root does not exist: $MESH_ROOT" >&2
    exit 1
fi

if [[ "$DRY_RUN" -eq 0 ]] && ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch was not found. Run this script on the cluster login node, or use --dry-run." >&2
    exit 1
fi

submitted=0
skipped=0
failed=0

shopt -s nullglob
block_dirs=( "$MESH_ROOT"/$PATTERN )
shopt -u nullglob

if [[ ${#block_dirs[@]} -eq 0 ]]; then
    echo "No block directories matched: $MESH_ROOT/$PATTERN"
    exit 0
fi

printf 'Mesh root: %s\n' "$MESH_ROOT"
printf 'Pattern:   %s\n' "$PATTERN"
printf 'Dry run:   %s\n' "$DRY_RUN"
printf 'Force:     %s\n\n' "$FORCE"

for block_dir in "${block_dirs[@]}"; do
    [[ -d "$block_dir" ]] || continue
    block_name="$(basename "$block_dir")"
    job_script="$block_dir/sub_intel.sh"

    if [[ ! -f "$job_script" ]]; then
        echo "SKIP $block_name: missing sub_intel.sh"
        ((skipped+=1))
        continue
    fi

    missing_inputs=()
    for f in fort.13 fort.14 fort.15 fort.19 fort.22; do
        [[ -f "$block_dir/$f" ]] || missing_inputs+=( "$f" )
    done
    if [[ ${#missing_inputs[@]} -gt 0 ]]; then
        echo "SKIP $block_name: missing inputs: ${missing_inputs[*]}"
        ((skipped+=1))
        continue
    fi

    if [[ "$FORCE" -eq 0 ]]; then
        if compgen -G "$block_dir/slurm-*.out" >/dev/null || [[ -f "$block_dir/fort.63" ]]; then
            echo "SKIP $block_name: existing slurm output or fort.63 found; use --force to resubmit"
            ((skipped+=1))
            continue
        fi
    fi

    if [[ "$MAX_SUBMIT" -gt 0 && "$submitted" -ge "$MAX_SUBMIT" ]]; then
        echo "Reached --max $MAX_SUBMIT; stopping."
        break
    fi

    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "DRY-RUN submit $block_name: (cd '$block_dir' && sbatch sub_intel.sh)"
        ((submitted+=1))
        continue
    fi

    echo "SUBMIT $block_name"
    if (cd "$block_dir" && sbatch sub_intel.sh); then
        ((submitted+=1))
    else
        echo "FAILED $block_name" >&2
        ((failed+=1))
    fi

    sleep "$SLEEP_SECONDS"
done

printf '\nDone. submitted=%d skipped=%d failed=%d\n' "$submitted" "$skipped" "$failed"
if [[ "$failed" -gt 0 ]]; then
    exit 1
fi
