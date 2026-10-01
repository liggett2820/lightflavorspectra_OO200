#!/bin/csh -f
#----------------------------------------------------------------------------------------------------
# Build the STAR Glauber/centrality libraries (with O+O) with cons.
# cons compiles packages under ./StRoot, so the packages in ./Makers are synced there first.
# Run from this directory on an RCF interactive node:  ./build_OO200.csh
#----------------------------------------------------------------------------------------------------
if ( ! $?GLAUBER_STARVER ) then
  echo "Set the STAR library version first, e.g.  setenv GLAUBER_STARVER SL19b"
  exit 1
endif
starver $GLAUBER_STARVER

mkdir -p StRoot
foreach pkg ( StCentralityMaker StFastGlauberMcMaker StGlauberAnalysisMaker StGlauberTree StGlauberUtilities )
  rsync -a --delete Makers/$pkg/ StRoot/$pkg/
end
cons
