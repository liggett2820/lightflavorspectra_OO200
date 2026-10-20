#!/bin/csh
#----------------------------------------------------------------------------------------------------
# One NBD scan job for O+O 200 GeV (called by submit_doScan_OO200.pl through condor)
# Arguments: npp k eff xbin xmin xmax [multCut [dataFile]]
#   x is scanned inside the job; multCut defaults to 15; dataFile defaults to hRefMult_OO200.root.
#   For the EPD totnMIP scan, dataFile = hTotnMIPCorr_OO200.root (its "hRefMult" is the corrected totnMIP,
#   1000 unit bins 0-1000, from macros/FitEpdCent_OO200.C).
#----------------------------------------------------------------------------------------------------
if ( $?GLAUBER_STARVER ) then
  starver $GLAUBER_STARVER
else
  starver SL24y
endif

set npp     = $1
set k       = $2
set eff     = $3
set xbin    = $4
set xmin    = $5
set xmax    = $6

set nevents = "1000000" # 1M per (npp, k, x, eff) point
set real    = "hRefMult_OO200.root"   # raw refMult after event cuts, from make_hRefMult_OO200.C
set mc      = "ncoll_npart.root"      # from addNcollVsNpart.C on the default O+O trees
set multCut = "15"                    # fit refMult >= 15 only (trigger/vertex inefficiency at low refMult); vary 10/20 as a check
if ( $#argv >= 7 ) set multCut = "$7"      # optional 7th argument overrides it (submit_doScan_OO200.pl -multcut N)
if ( $#argv >= 8 ) set real    = "$8"      # optional 8th argument: data file (submit_doScan_OO200.pl -totnmip)
if ( ! -e $real ) then
  echo "doScanX_OO200.csh: data file $real not found in `pwd`"
  exit 1
endif

root4star -b <<EOF2
  .L doNbdFitMaker_OO200.C
  scan($nevents, "$real", "$mc", $multCut, 1, $npp, $npp, 1, $k, $k, $xbin, $xmin, $xmax, $eff, 1.00, kFALSE);
  .q
EOF2
