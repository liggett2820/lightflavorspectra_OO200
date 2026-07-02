#!/bin/bash
# run_picobinner_sdcc.bash -- lightflavorspectra_OO200
#
# Thin wrapper invoked from inside xml/runPicoBinner_OO200_SDCC_template.xml's
# <command> block, on an SDCC/RCF worker node. It exists only to (a) resolve paths
# relative to this repo's actual layout the way LoadPicoBinnerLibs.C/RunPicoBinner.C
# expect (they're written to be run from inside macros/, per the main README) and
# (b) name the output file the way the XML's <output fromScratch="*_Processed.root">
# glob expects -- mirroring lightflavorspectra_etof's scripts/runPicoBinner.bash
# calling convention, but pointed at THIS repo's actual macro pair (RunPicoBinner.C +
# LoadPicoBinnerLibs.C) instead of etof's RunPicoBinner_Liggett.cxx.
#
# Not meant to be run standalone on a laptop -- for local runs use
# scripts/run_picobinner_resilient.sh instead (chunking + retry logic that a single
# SDCC job doesn't need, since the scheduler itself already parallelizes across many
# jobs/chunks of the filelist).
#
# Args (all positional, matching RunPicoBinner.C's signature):
#   $1  fileList                       path to the (already-awk'd, one-file-per-line) filelist
#   $2  outputDir                      directory to write the _Processed.root file into (the XML passes ${SCRATCH}/)
#   $3  partName                       PION | KAON | PROTON | ALL
#   $4  makeBasicHistos                true | false
#   $5  correctT0                      true | false
#   $6  ignoreBTOFOfNonParticleOfInterest   true | false
#   $7  doMassSquared                  true | false

set -u

if [ $# -lt 7 ]; then
  echo "Usage: $0 <fileList> <outputDir> <PION|KAON|PROTON|ALL> <makeBasicHistos> <correctT0> <ignoreBTOFNonPOI> <doMassSquared>"
  exit 1
fi

FILE_LIST="$1"
OUTPUT_DIR="$2"
PART_NAME="$3"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MACROS_DIR="$REPO_ROOT/macros"

OUT_BASE=$(basename "$FILE_LIST" .list)
OUT_FILE="${OUTPUT_DIR%/}/${OUT_BASE}_${PART_NAME}_Processed.root"

echo "run_picobinner_sdcc.bash: fileList=$FILE_LIST partName=$PART_NAME outFile=$OUT_FILE"

(
  cd "$MACROS_DIR" &&
  root -l -b -q LoadPicoBinnerLibs.C \
    "RunPicoBinner.C(\"$FILE_LIST\",\"$OUT_FILE\",9,\"$PART_NAME\",\"$4\",\"$5\",\"$6\",\"$7\")"
)
status=$?

echo "run_picobinner_sdcc.bash: root exited $status"
echo "End of run_picobinner_sdcc.bash"
exit $status
