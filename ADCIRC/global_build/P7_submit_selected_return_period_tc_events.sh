#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "$BASE_DIR/../submit_sub_intel_jobs.sh" ]]; then
  WORK_DIR="$(cd "$BASE_DIR/.." && pwd)"
elif [[ -f "$BASE_DIR/../../submit_sub_intel_jobs.sh" ]]; then
  WORK_DIR="$(cd "$BASE_DIR/../.." && pwd)"
else
  echo "Could not find submit_sub_intel_jobs.sh from $BASE_DIR" >&2
  exit 2
fi
BLOCK_PATTERN="${1:-ADC_*}"
MAX_JOBS="${MAX_JOBS:-490}"
MAX_SUBMIT="${MAX_SUBMIT:-0}"
DRY_RUN_ARG="${DRY_RUN_ARG:-}"

if [[ -d "$BASE_DIR/100yr" ]]; then
  bash "$WORK_DIR/submit_sub_intel_jobs.sh" --root "$BASE_DIR/100yr" --block-pattern "$BLOCK_PATTERN" --max-jobs "$MAX_JOBS" --max-submit "$MAX_SUBMIT" $DRY_RUN_ARG
fi
if [[ -d "$BASE_DIR/200yr" ]]; then
  bash "$WORK_DIR/submit_sub_intel_jobs.sh" --root "$BASE_DIR/200yr" --block-pattern "$BLOCK_PATTERN" --max-jobs "$MAX_JOBS" --max-submit "$MAX_SUBMIT" $DRY_RUN_ARG
fi
if [[ -d "$BASE_DIR/500yr" ]]; then
  bash "$WORK_DIR/submit_sub_intel_jobs.sh" --root "$BASE_DIR/500yr" --block-pattern "$BLOCK_PATTERN" --max-jobs "$MAX_JOBS" --max-submit "$MAX_SUBMIT" $DRY_RUN_ARG
fi
