#!/bin/bash
# hadd_picobinner_sdcc_output.sh -- lightflavorspectra_OO200
#
# WHY THIS EXISTS:
# The PicoBinner SDCC batch job (xml/runPicoBinner_OO200_SDCC_template.xml) splits the
# catalog query's input files across many worker-node processes (maxFilesPerProcess="12"),
# each of which writes its own separate *_Processed.root output file into the job's
# permanent output area (currently /star/data03/pwg/liggett/PicoBinner_OO200/SL24y_<PARTNAME>/,
# per that XML's <output> tag). Every downstream stage (ZFitter's RunZFitter.C, in
# particular) expects ONE yield file per species, not dozens of per-job fragments -- this
# script merges them with `hadd`, which correctly SUMS histograms of the same name across
# files, exactly equivalent to having run PicoBinner once over the combined input. This is
# unrelated to run_picobinner_resilient.sh's chunk-merge step -- that script chunks and
# merges a LOCAL run; this one merges an SDCC batch job's already-produced output.
#
# USAGE:
#   ./scripts/hadd_picobinner_sdcc_output.sh <inputDir> <outputFile> [pattern]
#
#   <inputDir>    directory containing the per-job PicoBinner outputs (e.g.
#                 /star/data03/pwg/liggett/PicoBinner_OO200/SL24y_PION/)
#   <outputFile>  path to write the single merged ROOT file to (e.g.
#                 /star/data03/pwg/liggett/PicoBinner_OO200/yieldHistos_OO200_pion.root)
#   [pattern]     glob pattern (relative to inputDir) matching the files to merge
#                 (default: *_Processed.root, matching the SDCC XML's <output> tag)
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
