#!/usr/bin/env bash
# Collect completed return-period rerun outputs.
#
# Default layout:
#   100yr/<block>/<case>/fort.63
#   100yr/<block>/<case>/maxele.63
#   500yr/<block>/<case>/fort.63
#   500yr/<block>/<case>/maxele.63
#
# Output layout:
#   return_out/<rp>/<block>/<case>/fort.63
#   return_out/<rp>/<block>/<case>/maxele.63
#
# The script moves fort.63 and maxele.63 by default, but never removes the
# original case directories.

#SBATCH -J collect_return_out
#SBATCH -o collect_return_out_%j.out
#SBATCH -e collect_return_out_%j.err

set -euo pipefail
shopt -s nullglob

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RP_DIRS="${RP_DIRS:-100yr 500yr}"
BLOCK_PATTERN="${BLOCK_PATTERN:-ADC_*}"
CASE_PATTERN="${CASE_PATTERN:-track_*}"
OUT_ROOT="${OUT_ROOT:-$BASE_DIR/return_out}"
DRY_RUN="${DRY_RUN:-0}"
COPY_MODE="${COPY_MODE:-0}"
OVERWRITE="${OVERWRITE:-0}"
MIN_AGE_SECONDS="${MIN_AGE_SECONDS:-0}"
LOG_STAMP="$(date '+%Y%m%d_%H%M%S')"
LOG_CSV="${LOG_CSV:-$BASE_DIR/collect_return_out_log_${LOG_STAMP}.csv}"
MISSING_CSV="${MISSING_CSV:-$BASE_DIR/collect_return_out_missing_cases.csv}"

usage() {
    cat <<EOF
Usage: bash $(basename "$0") [options]

Move fort.63 and maxele.63 from completed selected return-period reruns into
return_out, preserving return-period/block/case subfolders.

Options:
  --rp-dirs "LIST"       Return-period folders to scan. Default: "100yr 500yr"
  --block-pattern GLOB   Block directory glob. Default: ADC_*
  --case-pattern GLOB    Case directory glob. Default: track_*
  --out-root DIR         Output root. Default: ./return_out
  --copy                 Copy files instead of moving them.
  --overwrite            Overwrite existing collected files.
  --min-age SEC          Only move source files older than SEC seconds. Default: 0
  --dry-run              Print/log actions without changing files.
  -h, --help             Show this help.

Examples:
  bash $(basename "$0") --dry-run
  bash $(basename "$0")
  sbatch $(basename "$0")
  OUT_ROOT=./work/return_out bash $(basename "$0") --copy

This script does not delete the original case directories.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --rp-dirs)
            RP_DIRS="$2"; shift 2 ;;
        --block-pattern)
            BLOCK_PATTERN="$2"; shift 2 ;;
        --case-pattern)
            CASE_PATTERN="$2"; shift 2 ;;
        --out-root)
            OUT_ROOT="$2"; shift 2 ;;
        --copy)
            COPY_MODE=1; shift ;;
        --overwrite)
            OVERWRITE=1; shift ;;
        --min-age)
            MIN_AGE_SECONDS="$2"; shift 2 ;;
        --dry-run)
            DRY_RUN=1; shift ;;
        -h|--help)
            usage; exit 0 ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2 ;;
    esac
done

csv_escape() {
    local s="${1:-}"
    s="${s//\"/\"\"}"
    printf '"%s"' "$s"
}

csv_row() {
    local first=1
    local field
    for field in "$@"; do
        if [[ "$first" -eq 0 ]]; then
            printf ','
        fi
        csv_escape "$field"
        first=0
    done
    printf '\n'
}

file_size_bytes() {
    local f="$1"
    if [[ -f "$f" ]]; then
        stat -c '%s' "$f" 2>/dev/null || wc -c < "$f" | awk '{print $1}'
    else
        printf ''
    fi
}

file_age_seconds() {
    local f="$1"
    local now mtime
    if [[ ! -e "$f" ]]; then
        printf ''
        return
    fi
    now="$(date '+%s')"
    mtime="$(stat -c '%Y' "$f" 2>/dev/null || printf '%s' "$now")"
    printf '%s' "$((now - mtime))"
}

is_done_status() {
    case "$1" in
        moved|copied|dry_run_move|dry_run_copy|already_collected)
            return 0 ;;
        *)
            return 1 ;;
    esac
}

transfer_one_file() {
    local src="$1"
    local dst="$2"
    local dst_dir age
    dst_dir="$(dirname "$dst")"

    if [[ -f "$src" ]]; then
        if [[ ! -s "$src" ]]; then
            printf 'empty_source'
            return
        fi
        if [[ "$MIN_AGE_SECONDS" -gt 0 ]]; then
            age="$(file_age_seconds "$src")"
            if [[ -n "$age" && "$age" -lt "$MIN_AGE_SECONDS" ]]; then
                printf 'too_new'
                return
            fi
        fi
        if [[ -f "$dst" && "$OVERWRITE" -eq 0 ]]; then
            printf 'dest_exists'
            return
        fi
        if [[ "$DRY_RUN" -eq 1 ]]; then
            if [[ "$COPY_MODE" -eq 1 ]]; then
                printf 'dry_run_copy'
            else
                printf 'dry_run_move'
            fi
            return
        fi
        mkdir -p "$dst_dir"
        if [[ "$COPY_MODE" -eq 1 ]]; then
            cp -p -f "$src" "$dst"
            printf 'copied'
        else
            mv -f "$src" "$dst"
            printf 'moved'
        fi
        return
    fi

    if [[ -s "$dst" ]]; then
        printf 'already_collected'
    elif [[ -f "$dst" ]]; then
        printf 'empty_collected'
    else
        printf 'missing'
    fi
}

printf 'Base dir        : %s\n' "$BASE_DIR"
printf 'RP dirs         : %s\n' "$RP_DIRS"
printf 'Block pattern   : %s\n' "$BLOCK_PATTERN"
printf 'Case pattern    : %s\n' "$CASE_PATTERN"
printf 'Output root     : %s\n' "$OUT_ROOT"
printf 'Mode            : %s\n' "$([[ "$COPY_MODE" -eq 1 ]] && printf copy || printf move)"
printf 'Overwrite       : %s\n' "$OVERWRITE"
printf 'Min age seconds : %s\n' "$MIN_AGE_SECONDS"
printf 'Dry run         : %s\n' "$DRY_RUN"
printf 'Log CSV         : %s\n' "$LOG_CSV"
printf 'Missing CSV     : %s\n' "$MISSING_CSV"

if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "$OUT_ROOT"
fi

csv_row timestamp rp_tag block_id track_case case_dir dest_dir \
    fort63_status fort63_src fort63_dst fort63_size_bytes \
    maxele63_status maxele63_src maxele63_dst maxele63_size_bytes > "$LOG_CSV"
csv_row timestamp rp_tag block_id track_case case_dir missing_or_not_collected > "$MISSING_CSV"

cases_scanned=0
cases_complete=0
cases_partial=0
cases_missing_all=0
fort63_done=0
maxele63_done=0
dest_exists_count=0
empty_source_count=0
too_new_count=0
missing_file_count=0

for rp_dir in $RP_DIRS; do
    root="$BASE_DIR/$rp_dir"
    if [[ ! -d "$root" ]]; then
        echo "WARN: missing return-period folder, skipped: $root" >&2
        continue
    fi

    for block_dir in "$root"/$BLOCK_PATTERN; do
        [[ -d "$block_dir" ]] || continue
        block_id="$(basename "$block_dir")"

        for case_dir in "$block_dir"/$CASE_PATTERN; do
            [[ -d "$case_dir" ]] || continue
            track_case="$(basename "$case_dir")"
            dest_dir="$OUT_ROOT/$rp_dir/$block_id/$track_case"
            fort_src="$case_dir/fort.63"
            maxele_src="$case_dir/maxele.63"
            fort_dst="$dest_dir/fort.63"
            maxele_dst="$dest_dir/maxele.63"

            cases_scanned=$((cases_scanned + 1))
            fort_size="$(file_size_bytes "$fort_src")"
            maxele_size="$(file_size_bytes "$maxele_src")"
            fort_status="$(transfer_one_file "$fort_src" "$fort_dst")"
            maxele_status="$(transfer_one_file "$maxele_src" "$maxele_dst")"

            csv_row "$(date '+%F %T')" "$rp_dir" "$block_id" "$track_case" "$case_dir" "$dest_dir" \
                "$fort_status" "$fort_src" "$fort_dst" "$fort_size" \
                "$maxele_status" "$maxele_src" "$maxele_dst" "$maxele_size" >> "$LOG_CSV"

            if is_done_status "$fort_status"; then
                fort63_done=$((fort63_done + 1))
            fi
            if is_done_status "$maxele_status"; then
                maxele63_done=$((maxele63_done + 1))
            fi
            if [[ "$fort_status" == "dest_exists" || "$maxele_status" == "dest_exists" ]]; then
                dest_exists_count=$((dest_exists_count + 1))
            fi
            if [[ "$fort_status" == "empty_source" || "$maxele_status" == "empty_source" ]]; then
                empty_source_count=$((empty_source_count + 1))
            fi
            if [[ "$fort_status" == "too_new" || "$maxele_status" == "too_new" ]]; then
                too_new_count=$((too_new_count + 1))
            fi

            missing=()
            is_done_status "$fort_status" || missing+=("fort.63=$fort_status")
            is_done_status "$maxele_status" || missing+=("maxele.63=$maxele_status")

            if [[ ${#missing[@]} -eq 0 ]]; then
                cases_complete=$((cases_complete + 1))
            elif [[ ${#missing[@]} -eq 2 ]]; then
                cases_missing_all=$((cases_missing_all + 1))
                missing_file_count=$((missing_file_count + 2))
                csv_row "$(date '+%F %T')" "$rp_dir" "$block_id" "$track_case" "$case_dir" "${missing[*]}" >> "$MISSING_CSV"
            else
                cases_partial=$((cases_partial + 1))
                missing_file_count=$((missing_file_count + 1))
                csv_row "$(date '+%F %T')" "$rp_dir" "$block_id" "$track_case" "$case_dir" "${missing[*]}" >> "$MISSING_CSV"
            fi

            printf '[%s] %s/%s/%s fort.63=%s maxele.63=%s\n' \
                "$(date '+%F %T')" "$rp_dir" "$block_id" "$track_case" "$fort_status" "$maxele_status"
        done
    done
done

printf '\nSummary\n'
printf '  cases scanned       : %d\n' "$cases_scanned"
printf '  cases complete      : %d\n' "$cases_complete"
printf '  cases partial       : %d\n' "$cases_partial"
printf '  cases missing all   : %d\n' "$cases_missing_all"
printf '  fort.63 collected   : %d\n' "$fort63_done"
printf '  maxele.63 collected : %d\n' "$maxele63_done"
printf '  dest-exists cases   : %d\n' "$dest_exists_count"
printf '  empty-source cases  : %d\n' "$empty_source_count"
printf '  too-new cases       : %d\n' "$too_new_count"
printf '  not-collected files : %d\n' "$missing_file_count"
printf '  output root         : %s\n' "$OUT_ROOT"
printf '  log CSV             : %s\n' "$LOG_CSV"
printf '  missing CSV         : %s\n' "$MISSING_CSV"
