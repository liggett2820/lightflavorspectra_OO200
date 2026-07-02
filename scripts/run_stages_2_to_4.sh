#!/bin/bash
# run_stages_2_to_4.sh -- lightflavorspectra_OO200
#
# WHY THIS EXISTS:
# Stage 1 (PicoBinner) has its own dedicated driver, scripts/run_picobinner_resilient.sh,
# because of that stage's chunking/crash-retry needs. Stages 2-4 (ZFitter, RawSpectraModifier,
# SpectraFitter) have no such complication -- they're just a fixed sequence of `root`
# invocations that has to happen in the right order with the right files handed from one
# stage to the next (see the top-level README's "Pipeline stages" section for the full
# explanation of each stage and why the arguments below are what they are). This script
# is that fixed sequence, so you don't have to retype the six ZFitter invocations (one per
# particle x detector combination) and the hadd merge step by hand every time.
#
# USAGE:
#   ./scripts/run_stages_2_to_4.sh <pionYield.root> <kaonYield.root> <protonYield.root> \
#       <etaBetaGammaDeDxFile.root> <efficiencyFile.root> <workDir>
#
#   <pionYield.root>/<kaonYield.root>/<protonYield.root>  stage 1 (PicoBinner) outputs
#   <etaBetaGammaDeDxFile.root>   ZFitter's eta/betaGamma/dEdx lookup input (see RunZFitter.C)
#   <efficiencyFile.root>         embedding-derived TPC/BTOF efficiency + energy-loss
#                                 corrections for THIS dataset -- not shipped with this
#                                 repo, you supply your own (see RunRawSpectraModifier.C)
#   <workDir>                     directory for all intermediate + final outputs (created
#                                 if it doesn't exist)
#
# OUTPUT (all written into <workDir>):
#   zfitter_out/spectra_OO200_*.root   -- six per-(particle,detector) ZFitter outputs
#   spectra_combined.root              -- hadd of the six above (stage 3's input)
#   spectra_modified.root              -- RawSpectraModifier output (stage 4's input)
#   fit_output.root                    -- SpectraFitter output (final result)
#   run_log.txt                        -- combined stdout/stderr from every step
#
# This script uses RawSpectraModifier's and SpectraFitter's DEFAULT arguments for
# everything not passed in here (momentum ranges, rapidity cuts, fit seeds, etc. -- see
# each Run*.C's own header comment for what those defaults are and which of them you
# should sanity-check/override for your specific dataset before trusting the final
# output). Edit this script directly if you need to pass non-default values through.

set -u

if [ $# -lt 6 ]; then
  echo "Usage: $0 <pionYield.root> <kaonYield.root> <protonYield.root> <etaBetaGammaDeDxFile.root> <efficiencyFile.root> <workDir>"
  exit 1
fi

PION_YIELD="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
KAON_YIELD="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
PROTON_YIELD="$(cd "$(dirname "$3")" && pwd)/$(basename "$3")"
DEDX_FILE="$(cd "$(dirname "$4")" && pwd)/$(basename "$4")"
EFFICIENCY_FILE="$(cd "$(dirname "$5")" && pwd)/$(basename "$5")"
WORK_DIR="$6"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MACROS_DIR="$REPO_ROOT/macros"

mkdir -p "$WORK_DIR"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
ZFITTER_OUT_DIR="$WORK_DIR/zfitter_out/"
mkdir -p "$ZFITTER_OUT_DIR"

RUN_LOG="$WORK_DIR/run_log.txt"
: > "$RUN_LOG"

run_root() {
  echo "+++ root -l -q -b $* +++" | tee -a "$RUN_LOG"
  ( cd "$MACROS_DIR" && root -l -q -b "$@" ) >> "$RUN_LOG" 2>&1
  status=$?
  if [ $status -ne 0 ]; then
    echo "!!! FAILED (exit $status): root -l -q -b $* -- see $RUN_LOG" >&2
    exit $status
  fi
}

echo "===== Stage 2: ZFitter (6 invocations: pion/kaon/proton x TPC/BTOF) ====="
for partIndex in 0 1 2; do
  for detIndex in 0 1; do
    run_root LoadZFitterLibs.C "RunZFitter.C($partIndex,$detIndex,\"$PION_YIELD\",\"$KAON_YIELD\",\"$PROTON_YIELD\",\"$ZFITTER_OUT_DIR\",\"$DEDX_FILE\")"
  done
done

echo "===== Merging the six ZFitter spectra files with hadd ====="
SPECTRA_COMBINED="$WORK_DIR/spectra_combined.root"
rm -f "$SPECTRA_COMBINED"
hadd "$SPECTRA_COMBINED" "$ZFITTER_OUT_DIR"/spectra_OO200_*.root >> "$RUN_LOG" 2>&1
if [ $? -ne 0 ]; then
  echo "!!! hadd FAILED -- see $RUN_LOG" >&2
  exit 1
fi

echo "===== Stage 3: RawSpectraModifier ====="
SPECTRA_MODIFIED="$WORK_DIR/spectra_modified.root"
run_root LoadRawSpectraModifierLibs.C "RunRawSpectraModifier.C(\"$PION_YIELD\",\"$KAON_YIELD\",\"$PROTON_YIELD\",\"$SPECTRA_COMBINED\",\"$EFFICIENCY_FILE\",\"$SPECTRA_MODIFIED\")"

echo "===== Stage 4: SpectraFitter ====="
FIT_OUTPUT="$WORK_DIR/fit_output.root"
run_root LoadSpectraFitterLibs.C "RunSpectraFitter.C(\"$SPECTRA_MODIFIED\",\"$FIT_OUTPUT\")"

echo "===== Done. Final output: $FIT_OUTPUT ====="
echo "Full log: $RUN_LOG"
