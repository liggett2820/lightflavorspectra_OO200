#!/bin/bash
# hadd_chunk_from_filelist.sh -- lightflavorspectra_OO200
#
# WHY THIS EXISTS:
# Companion to xml/haddPicoBinner_SDCC_template.xml's parallel chunk-then-merge strategy
# for combining PicoBinner's many per-job SDCC output fragments (*_Processed.root)
# faster than a single serial `hadd` over all of them. SUMS splits the full file list
# across N worker processes (maxFilesPerProcess in that XML); each process calls this
# script with ITS OWN slice of the file list (${FILELIST}, populated automatically by
# the scheduler) and hadds just that slice into one chunk file. Once every process's
# chunk has landed in the permanent output directory, run
#   scripts/hadd_picobinner_sdcc_output.sh <chunkDir> <finalOutput> "chunk_*.root"
# to do the final merge -- much cheaper than the original merge, since it's now
# combining only a handful of chunk files instead of hundreds of small fragments.
#
# USAGE:
#   ./scripts/hadd_chunk_from_filelist.sh <fileList> <outputChunkFile>
#
#   <fileList>         text file, one input *_Processed.root path per line (this is
#                      ${FILELIST}, provided automatically by the SUMS scheduler)
#   <outputChunkFile>  path to write this process's merged chunk to

set -u

if [ $# -lt 2 ]; then
  echo "Usage: $0 <fileList> <outputChunkFile>" >&2
  exit 1
fi

FILE_LIST="$1"
OUTPUT_CHUNK="$2"

if [ ! -f "$FILE_LIST" ]; then
  echo "!!! File list does not exist: $FILE_LIST" >&2
  exit 1
fi

# awk '{print $1}' in case the list has a second column (e.g. event counts) -- matches
# the same defensive habit used for catalog-derived FILELISTs elsewhere in this repo.
mapfile -t INPUT_FILES < <(awk '{print $1}' "$FILE_LIST")

if [ ${#INPUT_FILES[@]} -eq 0 ]; then
  echo "!!! File list $FILE_LIST is empty -- nothing to merge." >&2
  exit 1
fi

echo "Merging ${#INPUT_FILES[@]} file(s) from $FILE_LIST into $OUTPUT_CHUNK ..."
hadd -f "$OUTPUT_CHUNK" "${INPUT_FILES[@]}"
status=$?

if [ $status -ne 0 ]; then
  echo "!!! hadd FAILED (exit $status)" >&2
  exit $status
fi

echo "Done: $OUTPUT_CHUNK"
