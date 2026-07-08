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
# submodule has its own Makefile and its outputs (libStPicoDst.so,
# StPicoDst_Dict_rdict.pcm) just get gSystem->Load()'ed from bin/, so they have to be
# rebuilt and copied back in every time bin/ gets wiped, or makeLibs.C fails immediately
# with "libStPicoDst.so does not exist".
#
# That submodule's own Makefile only tracks libStPicoDst.so as a target -- if it's
# already sitting in submodule/PicoDstReader_SL24y/ from a previous build, plain `make`
# reports "Nothing to be done" and skips the rootcint step entirely, so
# StPicoDst_Dict_rdict.pcm (a side-effect file, not its own tracked target) never gets
# (re)created and the later `cp` fails with "No such file or directory". `make clean`
# first forces both the dictionary and the library to actually regenerate.
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

echo "=== [3/4] Copying libStPicoDst.so + StPicoDst_Dict_rdict.pcm into bin/ ==="
mkdir -p bin
cp submodule/PicoDstReader_SL24y/libStPicoDst.so bin/
cp submodule/PicoDstReader_SL24y/StPicoDst_Dict_rdict.pcm bin/

echo "=== [4/4] Compiling everything else (root -l -q -b macros/makeLibs.C) ==="
root -l -q -b macros/makeLibs.C

echo ""
echo "Done. Contents of bin/ (check that .so timestamps all match this run):"
ls -la bin/*.so 2>/dev/null
