#!/bin/bash
# hadd_final_merge.sh -- lightflavorspectra_OO200
# (formerly hadd_picobinner_sdcc_output.sh -- renamed 2026-07-11 so the STAGE 1/STAGE 2
# split with hadd_stage1_chunk_worker.sh is obvious from the filenames alone.)
#
# WHY THIS EXISTS:
# The PicoBinner SDCC batch job (xml/runPicoBinner_OO200_SDCC_template.xml) splits the
# catalog query's input files across many worker-node processes (maxFilesPerProcess="12"),
# each of which writes its own separate *_Processed.root output file into the job's
# permanent output area (currently /star/data03/pwg/liggett/PicoBinner_OO200/SL24y_<PARTNAME>/,
# per that XML's <output> tag). Every downstream stage (ZFitter's RunZFitter.C, in
# particular) expects ONE yield file per species, not dozens of per-job fragments -- this
# script merges a directory of ROOT files with `hadd`, which correctly SUMS histograms of
# the same name across files, exactly equivalent to having run PicoBinner once over the
# combined input. This is unrelated to run_picobinner_resilient.sh's chunk-merge step --
# that script chunks and merges a LOCAL run; this one merges an SDCC batch job's
# already-produced output.
#
# This script is generic (just <inputDir>/<pattern> -> one hadd'd <outputFile>) and is
# actually used in TWO different places in the pipeline -- pick whichever matches what
# you're merging:
#
#   1. Directly on PicoBinner's raw per-worker fragments (fine for a modest fragment
#      count -- a single serial hadd over everything):
#        ./scripts/hadd_final_merge.sh \
#          /star/data03/pwg/liggett/PicoBinner_OO200/SL24y_PION/ \
#          /star/data03/pwg/liggett/PicoBinner_OO200/yieldHistos_OO200_pion.root
#        (pattern defaults to *_Processed.root, matching the SDCC XML's <output> tag)
#
#   2. As STAGE 2 of the two-stage "chunk-then-merge" workflow in
#      xml/haddPicoBinner_SDCC_template.xml (used when there are hundreds of fragments --
#      STAGE 1 of that XML fans a parallel hadd out across worker processes first,
#      producing a handful of chunk_*.root files instead of hundreds of fragments):
#        ./scripts/hadd_final_merge.sh \
#          /star/data03/pwg/liggett/PicoBinner_OO200/SL24y_PION_chunks/ \
#          /star/data03/pwg/liggett/PicoBinner_OO200/yieldHistos_OO200_pion.root \
#          "chunk_*.root"
#        (pattern must be passed explicitly here -- the default won't match chunk files)
#
# USAGE:
#   ./scripts/hadd_final_merge.sh <inputDir> <outputFile> [pattern] [--force]
#
#   <inputDir>    directory containing the ROOT files to merge -- either PicoBinner's raw
#                 per-job output (usage 1 above) or STAGE 1's chunk_*.root output
#                 (usage 2 above)
#   <outputFile>  path to write the single merged ROOT file to (e.g.
#                 /star/data03/pwg/liggett/PicoBinner_OO200/yieldHistos_OO200_pion.root)
#   [pattern]     glob pattern (relative to inputDir) matching the files to merge
#                 (default: *_Processed.root -- override to "chunk_*.root" for usage 2)
#
# Refuses to run (exits nonzero without touching anything) if no files match the pattern,
# or if <outputFile> already exists -- pass --force as the last argument to overwrite.

set -u

FORCE=0
args=()
for a in "$@"; do
  if [ "$a" == "--force" ]; then
    FORCE=1
  else
    args+=("$a")
  fi
done
set -- "${args[@]}"

if [ $# -lt 2 ]; then
  echo "Usage: $0 <inputDir> <outputFile> [pattern=*_Processed.root] [--force]" >&2
  exit 1
fi

INPUT_DIR="$1"
OUTPUT_FILE="$2"
PATTERN="${3:-*_Processed.root}"

if [ ! -d "$INPUT_DIR" ]; then
  echo "!!! Input directory does not exist: $INPUT_DIR" >&2
  exit 1
fi
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"

# Use find (not a bare glob) so this fails loudly and clearly if nothing matches,
# instead of hadd silently receiving a literal unexpanded "*_Processed.root" string
# as a single (nonexistent) filename argument.
mapfile -t MATCHED_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -name "$PATTERN" | sort)

if [ ${#MATCHED_FILES[@]} -eq 0 ]; then
  echo "!!! No files matching '$PATTERN' found in $INPUT_DIR -- nothing to merge." >&2
  exit 1
fi

echo "Found ${#MATCHED_FILES[@]} file(s) matching '$PATTERN' in $INPUT_DIR:"
for f in "${MATCHED_FILES[@]}"; do
  echo "  $f"
done

if [ -e "$OUTPUT_FILE" ] && [ $FORCE -ne 1 ]; then
  echo "!!! Output file already exists: $OUTPUT_FILE (pass --force to overwrite)" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_FILE")"

echo "Merging ${#MATCHED_FILES[@]} file(s) into $OUTPUT_FILE ..."
hadd -f "$OUTPUT_FILE" "${MATCHED_FILES[@]}"
status=$?

if [ $status -ne 0 ]; then
  echo "!!! hadd FAILED (exit $status)" >&2
  exit $status
fi

echo "Done: $OUTPUT_FILE"
