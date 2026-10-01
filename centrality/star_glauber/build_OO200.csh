#!/bin/csh -f
#----------------------------------------------------------------------------------------------------
# Build the STAR Glauber/centrality libraries (with O+O) with cons.
# cons compiles packages under ./StRoot, so the packages in ./Makers are synced there first.
# Run from this directory on an RCF interactive node:  ./build_OO200.csh
#----------------------------------------------------------------------------------------------------
# Run INSIDE the SDCC SL7 container (same one PicoBinner uses). Default SL24y, as for PicoBinner.
if ( ! $?GLAUBER_STARVER ) setenv GLAUBER_STARVER SL24y
starver $GLAUBER_STARVER

mkdir -p StRoot
foreach pkg ( StCentralityMaker StFastGlauberMcMaker StGlauberAnalysisMaker StGlauberTree StGlauberUtilities )
  rsync -a --delete Makers/$pkg/ StRoot/$pkg/
end
cons
