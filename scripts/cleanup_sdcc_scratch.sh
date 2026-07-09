#!/bin/bash
#
# cleanup_sdcc_scratch.sh -- lightflavorspectra_OO200
#
# Clears out accumulated batch-job scratch: scheduler stdout/stderr logs and
# Generator auxiliary (schedFiles) output that pile up from repeated
# submissions of xml/runPicoBinner_OO200_SDCC_template.xml AND
# xml/haddPicoBinner_SDCC_template.xml. Each star-submit call writes a fresh
# .condor class-ad file per process into its Generator <Location> dir and
# never cleans up the old ones -- across enough submissions (especially the
# hadd job's many-worker-process fan-out) this is what fills the scratch
# quota and makes star-submit itself fail with "Disk quota exceeded" before a
# job ever runs.
#
# Does NOT touch permanent physics output under
# /star/data03/pwg/liggett/PicoBinner_OO200/ -- only scratch-area
# log/bookkeeping clutter.
#
# Usage:
#   ./cleanup_sdcc_scratch.sh            # shows sizes, prompts for confirmation
#   ./cleanup_sdcc_scratch.sh -y         # skip confirmation
#   ./cleanup_sdcc_scratch.sh --dry-run  # show what would be removed, delete nothing

set -euo pipefail

LOG_DIR="/star/data03/scratch/liggett/picoBinner_OO200"
SCHED_DIR="/star/data03/scratch/liggett/schedFiles/picoBinner_OO200"
HADD_LOG_DIR="/star/data03/scratch/liggett/haddPicoBinner_OO200"
HADD_SCHED_DIR="/star/data03/scratch/liggett/schedFiles/haddPicoBinner_OO200"

DRY_RUN=false
ASSUME_YES=false

for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=true ;;
    -y|--yes) ASSUME_YES=true ;;
    -h|--help)
      echo "Usage: $0 [-y|--yes] [--dry-run]"
      exit 0
      ;;
    *)
      echo "Unknown option: $arg" >&2
      echo "Usage: $0 [-y|--yes] [--dry-run]" >&2
      exit 1
      ;;
  esac
done

dir_size() {
  local dir="$1"
  if [ -d "$dir" ]; then
    du -sh "$dir" 2>/dev/null | cut -f1
  else
    echo "(missing)"
  fi
}

echo "PicoBinner job log dir:        $LOG_DIR        [$(dir_size "$LOG_DIR")]"
echo "PicoBinner scheduler files dir: $SCHED_DIR  [$(dir_size "$SCHED_DIR")]"
echo "hadd job log dir:              $HADD_LOG_DIR        [$(dir_size "$HADD_LOG_DIR")]"
echo "hadd scheduler files dir:       $HADD_SCHED_DIR  [$(dir_size "$HADD_SCHED_DIR")]"
echo

if $DRY_RUN; then
  echo "[dry run] Would remove all contents of:"
  echo "  $LOG_DIR/*"
  echo "  $SCHED_DIR/*"
  echo "  $HADD_LOG_DIR/*"
  echo "  $HADD_SCHED_DIR/*"
  echo "Nothing was deleted."
  exit 0
fi

if ! $ASSUME_YES; then
  read -r -p "Delete all job logs and scheduler files listed above? [y/N] " reply
  case "$reply" in
    [yY][eE][sS]|[yY]) ;;
    *) echo "Aborted -- nothing deleted."; exit 0 ;;
  esac
fi

for dir in "$LOG_DIR" "$SCHED_DIR" "$HADD_LOG_DIR" "$HADD_SCHED_DIR"; do
  if [ -d "$dir" ]; then
    echo "Clearing $dir ..."
    find "$dir" -mindepth 1 -delete
  else
    echo "$dir does not exist yet -- skipping."
  fi
done

echo "Done. Permanent output under /star/data03/pwg/liggett/PicoBinner_OO200/ was not touched."
