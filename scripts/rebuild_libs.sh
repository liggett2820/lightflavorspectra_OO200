#!/bin/bash
# rebuild_libs.sh -- lightflavorspectra_OO200
#
# WHY THIS EXISTS:
# `root -l -q -b macros/makeLibs.C` compiles most of this repo via ACLiC's
# gSystem->CompileMacro(), which only recompiles a file if it thinks the source is newer
# than the .so already sitting in bin/. If an updated source file (e.g. SetCutClass.C)
# reaches SDCC in a way that preserves an older mtime (rsync/scp -p, etc.), that check
# silently passes and the OLD compiled binary keeps running even though the source text
# on disk is correct -- this bit us once already with a stale centrality-cut scheme (see
# README's "Centrality binning" section, "UPDATE 2026-07-07").
#
# This script forces a real rebuild every time by default: it wipes bin/* via
# `makeLibs.C("clean")` before recompiling, so there's no stale-.so ambiguity to debug.
# It also handles the one build step makeLibs.C does NOT do itself -- the PicoDst reader
# submodule has its own Makefile and its output (libStPicoDst.so) just gets
# gSystem->Load()'ed from bin/, so it has to be rebuilt and copied back in every time
# bin/ gets wiped, or makeLibs.C fails immediately with "libStPicoDst.so does not exist".
#
# That submodule's own Makefile only tracks libStPicoDst.so as a target -- if it's
# already sitting in submodule/PicoDstReader_SL24y/ from a previous build, plain `make`
# reports "Nothing to be done" and skips the rootcint step entirely. `make clean` first
# forces the dictionary and library to actually regenerate.
#
# StPicoDst_Dict_rdict.pcm is NOT part of this: xml/runPicoBinner_OO200_SDCC_template.xml
# (lines ~94-102) documents, from an actual SDCC build directory listing, that SDCC's
# ROOT/rootcling never produces that file for this submodule at all -- it's a macOS-local
# build sidecar only, and isn't needed on SDCC (confirmed by lightflavorspectra_etof's own
# working production job, which packages only libStPicoDst.so). So this script copies it
# opportunistically if present and just skips it otherwise, instead of treating it as
# required.
#
# USAGE (run from anywhere -- it cd's to the repo root itself):
#   ./scripts/rebuild_libs.sh              # full clean rebuild (default, recommended)
#   ./scripts/rebuild_libs.sh --no-clean    # skip the bin/ wipe, just recompile changed
#                                           # files (faster, but re-introduces the same
#                                           # stale-.so risk this script exists to avoid --
#                                           # only use this if you're sure bin/ is already
#                                           # in sync with source)

set -eu

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

CLEAN=1
if [ "${1:-}" = "--no-clean" ]; then
  CLEAN=0
fi

echo "=== [1/4] Building PicoDstReader_SL24y submodule ==="
( cd submodule/PicoDstReader_SL24y && make clean && make )

if [ $CLEAN -eq 1 ]; then
  echo "=== [2/4] Clearing bin/ (root -l -q -b 'macros/makeLibs.C(\"clean\")') ==="
  root -l -q -b 'macros/makeLibs.C("clean")'
else
  echo "=== [2/4] Skipped (--no-clean) ==="
fi

echo "=== [3/4] Copying libStPicoDst.so into bin/ ==="
mkdir -p bin
cp submodule/PicoDstReader_SL24y/libStPicoDst.so bin/
if [ -f submodule/PicoDstReader_SL24y/StPicoDst_Dict_rdict.pcm ]; then
  cp submodule/PicoDstReader_SL24y/StPicoDst_Dict_rdict.pcm bin/
  echo "    (also copied StPicoDst_Dict_rdict.pcm -- present on this build)"
else
  echo "    (no StPicoDst_Dict_rdict.pcm produced here -- expected on SDCC, not needed)"
fi

echo "=== [4/4] Compiling everything else (root -l -q -b macros/makeLibs.C) ==="
root -l -q -b macros/makeLibs.C

echo ""
echo "Done. Contents of bin/ (check that .so timestamps all match this run):"
ls -la bin/*.so 2>/dev/null
