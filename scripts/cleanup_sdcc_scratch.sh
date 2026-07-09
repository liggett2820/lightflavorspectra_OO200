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
# job ever runs (errno 122 / EDQUOT -- this is a quota on the scratch
# fileset itself, separate from the `quota -s` home-directory quota, and
# apparently small enough that a few GB of stale scratch clutter is enough
# to trip it).
#
# Also cleans up two things discovered during a 2026-07-09 investigation:
#
#   1. A stray lowercase "schedfiles" directory alongside the real
#      "schedFiles" one -- an orphaned duplicate (nothing in this repo's
#      XML templates or this script ever reference the lowercase name) that
#      was quietly eating multiple GB with old .condor class-ad files.
#
#   2. Dated top-level scratch dirs (YYYY_MM_DD, e.g. "2026_07_01") that
#      match Condor's per-job ${SCRATCH} variable. The job XML is supposed
#      to move its output to permanent storage and self-clean this
#      (see the "mv .../*_Processed.root ... ; rm ..." lines in the XML
#      templates), but a job that gets held or fails before reaching that
#      line leaves the whole dated SCRATCH dir behind, full of
#      stdout_stderr logs (and possibly un-migrated *_Processed.root
#      output). This script only removes dated dirs that (a) aren't today's
#      date -- so it never touches a currently-running job's live SCRATCH --
#      and (b) contain no *.root files. Any dated dir that DOES contain
#      *.root files is left alone and flagged, since that likely means real
#      output never got moved to permanent storage and needs to be rescued
#      by hand first (see README's permanent-output paths).
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

SCRATCH_ROOT="/star/data03/scratch/liggett"

LOG_DIR="$SCRATCH_ROOT/picoBinner_OO200"
SCHED_DIR="$SCRATCH_ROOT/schedFiles/picoBinner_OO200"
HADD_LOG_DIR="$SCRATCH_ROOT/haddPicoBinner_OO200"
HADD_SCHED_DIR="$SCRATCH_ROOT/schedFiles/haddPicoBinner_OO200"
TYPO_SCHED_DIR="$SCRATCH_ROOT/schedfiles"

TODAY=$(date +%Y_%m_%d)

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

# Find dated per-job ${SCRATCH} leftovers (YYYY_MM_DD) other than today's,
# splitting them into ones that are safe to remove (no *.root inside) vs.
# ones that need to be rescued by hand first (contain *.root output that
# apparently never got moved to permanent storage).
DATED_SAFE=()
DATED_UNSAFE=()
if [ -d "$SCRATCH_ROOT" ]; then
  for dir in "$SCRATCH_ROOT"/[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9]; do
    [ -d "$dir" ] || continue
    base=$(basename "$dir")
    [ "$base" = "$TODAY" ] && continue
    if find "$dir" -name '*.root' -print -quit | grep -q .; then
      DATED_UNSAFE+=("$dir")
    else
      DATED_SAFE+=("$dir")
    fi
  done
fi

echo "PicoBinner job log dir:        $LOG_DIR        [$(dir_size "$LOG_DIR")]"
echo "PicoBinner scheduler files dir: $SCHED_DIR  [$(dir_size "$SCHED_DIR")]"
echo "hadd job log dir:              $HADD_LOG_DIR        [$(dir_size "$HADD_LOG_DIR")]"
echo "hadd scheduler files dir:       $HADD_SCHED_DIR  [$(dir_size "$HADD_SCHED_DIR")]"
echo "stray lowercase schedfiles dir: $TYPO_SCHED_DIR  [$(dir_size "$TYPO_SCHED_DIR")]"
if [ ${#DATED_SAFE[@]} -gt 0 ]; then
  echo "Stale dated SCRATCH dirs (no *.root inside, safe to remove):"
  for dir in "${DATED_SAFE[@]}"; do
    echo "  $dir  [$(dir_size "$dir")]"
  done
fi
if [ ${#DATED_UNSAFE[@]} -gt 0 ]; then
  echo
  echo "!!! Dated SCRATCH dirs containing *.root files -- NOT touched, rescue by hand:"
  for dir in "${DATED_UNSAFE[@]}"; do
    echo "  $dir  [$(dir_size "$dir")]"
  done
fi
echo

if $DRY_RUN; then
  echo "[dry run] Would remove all contents of:"
  echo "  $LOG_DIR/*"
  echo "  $SCHED_DIR/*"
  echo "  $HADD_LOG_DIR/*"
  echo "  $HADD_SCHED_DIR/*"
  echo "  $TYPO_SCHED_DIR (entirely)"
  for dir in "${DATED_SAFE[@]}"; do
    echo "  $dir (entirely)"
  done
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

if [ -d "$TYPO_SCHED_DIR" ]; then
  echo "Removing stray lowercase $TYPO_SCHED_DIR ..."
  rm -rf "$TYPO_SCHED_DIR"
fi

for dir in "${DATED_SAFE[@]}"; do
  echo "Removing stale dated SCRATCH dir $dir ..."
  rm -rf "$dir"
done

echo "Done. Permanent output under /star/data03/pwg/liggett/PicoBinner_OO200/ was not touched."
if [ ${#DATED_UNSAFE[@]} -gt 0 ]; then
  echo "Reminder: ${#DATED_UNSAFE[@]} dated SCRATCH dir(s) with *.root files were left in place -- see warning above."
fi
