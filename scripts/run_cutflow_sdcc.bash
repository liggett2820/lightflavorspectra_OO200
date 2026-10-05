#!/bin/bash
# run_cutflow_sdcc.bash -- lightflavorspectra_OO200
# Called from xml/runCutFlow_OO200_SDCC.xml on a worker node, from the unpacked SandBox root.
# Args: $1 file list (one picoDst per line), $2 output ROOT file
set -u
if [ $# -lt 2 ]; then echo "Usage: $0 <fileList> <outFile>"; exit 1; fi
# The SandBox packs submodule/PicoDstReader_SL24y at the top level; LoadCutFlowLibs.C
# expects it under submodule/.
if [ -d PicoDstReader_SL24y ] && [ ! -d submodule/PicoDstReader_SL24y ]; then
  mkdir -p submodule && ln -s ../PicoDstReader_SL24y submodule/PicoDstReader_SL24y
fi
root -l -b -q macros/LoadCutFlowLibs.C "macros/MakeEventCutFlow_OO200.C+(\"$1\",\"$2\")"
status=$?
echo "run_cutflow_sdcc.bash: root exited $status"
exit $status
