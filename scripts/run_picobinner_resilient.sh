#!/bin/bash
# run_picobinner_resilient.sh -- lightflavorspectra_OO200
#
# WHY THIS EXISTS:
# Running PicoBinner as one long-lived process over the full fastoffline filelist
# (3430 files, ~3M events) hits a non-deterministic segfault somewhere inside ROOT's
# TBranchElement/TClonesArray reading machinery (confirmed: identical binary + data
# crashed at different event indices across repeated runs -- see chat history / repo
# notes). It is NOT fixable from our code; it's effectively random. So instead of one
# process reading everything, this script runs PicoBinner in small chunks, each as its
# own `root` subprocess:
#   - a chunk that crashes only loses that chunk's events, not the whole run
#   - failed chunks are logged so they can be identified/retried (optionally at a
#     smaller chunk size to further isolate)
#   - successful chunks' output histograms get merged at the end with `hadd` (which
#     correctly SUMS histograms of the same name -- exactly equivalent to having
#     accumulated all those events in one PicoBinner run)
#
# USAGE:
#   ./scripts/run_picobinner_resilient.sh <master_filelist> <output_dir> [particle] [chunk_size] [max_retries]
#
#   <master_filelist>  path to the full .list file (e.g. filelist_OO_fastoffline.list)
#   <output_dir>       directory to write per-chunk outputs + logs + final merged file
#   [particle]         PION | KAON | PROTON | ALL   (default: PION)
#   [chunk_size]        number of input files per subprocess (default: 1 -- safest;
#                       raise this for speed once you've seen how often chunks fail,
#                       at the cost of losing more events per crash)
#   [max_retries]       times to retry a chunk that crashed before giving up on it
#                       (default: 5). Since the crash is confirmed non-deterministic
#                       (same binary+data crashes at different points on separate runs),
#                       a fresh subprocess retry has a real chance of succeeding even
#                       when the previous attempt on the exact same chunk crashed.
#
# OUTPUT:
#   <output_dir>/chunk_XXXXX.root       -- one per successful chunk
#   <output_dir>/failed_chunks.log      -- filelists of chunks that still failed after
#                                          all retries were exhausted
#   <output_dir>/run_log.txt            -- full stdout/stderr from every attempt
#   <output_dir>/merged.root            -- hadd of all successful chunk outputs
#
# Retrying failures further: rerun this script with a filelist built from
# failed_chunks.log (optionally with a smaller chunk_size, down to 1, and/or more
# max_retries) and hadd the new merged.root into the existing one.

set -u

if [ $# -lt 2 ]; then
  echo "Usage: $0 <master_filelist> <output_dir> [particle=PION] [chunk_size=1] [max_retries=5]"
  exit 1
fi

MASTER_FILELIST="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
OUTPUT_DIR="$2"
PARTICLE="${3:-PION}"
CHUNK_SIZE="${4:-1}"
MAX_RETRIES="${5:-5}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MACROS_DIR="$REPO_ROOT/macros"
BIN_DIR="$REPO_ROOT/bin"

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

FAILED_LOG="$OUTPUT_DIR/failed_chunks.log"
RUN_LOG="$OUTPUT_DIR/run_log.txt"
: > "$FAILED_LOG"
: > "$RUN_LOG"

export DYLD_LIBRARY_PATH="$BIN_DIR:${DYLD_LIBRARY_PATH:-}"

TMP_DIR="$(mktemp -d /tmp/picobinner_chunks_XXXXXX)"
trap 'rm -rf "$TMP_DIR"' EXIT

# Split the master filelist into chunks of CHUNK_SIZE lines each.
split -l "$CHUNK_SIZE" "$MASTER_FILELIST" "$TMP_DIR/chunk_"

NUM_CHUNKS=$(ls "$TMP_DIR"/chunk_* 2>/dev/null | wc -l | tr -d ' ')
echo "Split $MASTER_FILELIST into $NUM_CHUNKS chunk(s) of up to $CHUNK_SIZE file(s) each."
echo "Output dir: $OUTPUT_DIR"

i=0
n_ok=0
n_failed=0
total_attempts=0

# Checks whether a_1 is a complete, readable ROOT file with at least one key in it.
#
# WHY THIS EXISTS (found by inspecting run_log.txt for a chunk that needed all 5
# retries): source/PicoBinner.cxx calls outFile->Close() and prints "Wrote Histograms"
# BEFORE returning (see the lines right around "Wrote Histograms" in PicoBinner.cxx) --
# but on several attempts, the `root` process still exited non-zero AFTER printing that
# and "Finished PicoBinner Macro", i.e. the crash happens during ROOT's OWN process
# teardown (a known ROOT/cling shutdown-crash pattern), strictly after the output file
# was already fully written and closed. Trusting only the process exit status (the
# original check here) meant perfectly good, complete output files were being deleted
# and the chunk retried for no reason related to data integrity -- almost certainly the
# real explanation for the unexpectedly high permanent-failure rate seen in practice.
# Checking the file itself (not the exit code) avoids throwing away good output.
is_valid_root_output() {
  local f="$1"
  [ -f "$f" ] || return 1
  root -l -b -q -e "TFile* rf = TFile::Open(\"$f\"); bool ok = rf && !rf->IsZombie() && rf->GetNkeys() > 0; gSystem->Exit(ok ? 0 : 1);" > /dev/null 2>&1
  return $?
}

for chunk_list in "$TMP_DIR"/chunk_*; do
  i=$((i+1))
  chunk_name=$(printf "chunk_%05d" "$i")
  chunk_out="$OUTPUT_DIR/${chunk_name}.root"

  attempt=0
  success=0
  while [ $attempt -lt $MAX_RETRIES ]; do
    attempt=$((attempt+1))
    total_attempts=$((total_attempts+1))
    echo "=== [$i/$NUM_CHUNKS] attempt $attempt/$MAX_RETRIES: $(wc -l < "$chunk_list" | tr -d ' ') file(s) -> $chunk_out ===" | tee -a "$RUN_LOG"

    (
      cd "$MACROS_DIR" &&
      root -l -q -b LoadPicoBinnerLibs.C \
        "RunPicoBinner.C(\"$chunk_list\",\"$chunk_out\",9,\"$PARTICLE\",\"true\",\"false\",\"false\",\"false\")"
    ) >> "$RUN_LOG" 2>&1
    status=$?

    if [ $status -ne 0 ]; then
      echo "    (process exited $status -- checking if $chunk_out is a valid, complete file anyway before giving up on this attempt)" | tee -a "$RUN_LOG"
    fi

    if is_valid_root_output "$chunk_out"; then
      success=1
      break
    fi
    rm -f "$chunk_out"  # don't let a genuinely incomplete/corrupt file linger before retry
  done

  if [ $success -eq 1 ]; then
    n_ok=$((n_ok+1))
    if [ $attempt -gt 1 ]; then
      echo "    -> succeeded on retry $attempt/$MAX_RETRIES" | tee -a "$RUN_LOG"
    fi
  else
    n_failed=$((n_failed+1))
    echo "!!! Chunk $i FAILED after $MAX_RETRIES attempts. Files:" >> "$FAILED_LOG"
    cat "$chunk_list" >> "$FAILED_LOG"
    echo "---" >> "$FAILED_LOG"
  fi
done

echo ""
echo "Done: $n_ok/$NUM_CHUNKS chunks succeeded, $n_failed failed permanently (after up to $MAX_RETRIES attempts each). Total subprocess launches: $total_attempts."
if [ $n_failed -gt 0 ]; then
  echo "Failed chunk file lists logged in: $FAILED_LOG"
fi

CHUNK_OUTPUTS=("$OUTPUT_DIR"/chunk_*.root)
if [ -e "${CHUNK_OUTPUTS[0]}" ]; then
  echo "Merging $n_ok successful chunk output(s) into $OUTPUT_DIR/merged.root ..."
  hadd -f "$OUTPUT_DIR/merged.root" "${CHUNK_OUTPUTS[@]}"
else
  echo "No successful chunk outputs to merge."
fi
