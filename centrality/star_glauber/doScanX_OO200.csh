#!/bin/csh
#----------------------------------------------------------------------------------------------------
# One NBD scan job for O+O 200 GeV (called by submit_doScan_OO200.pl through condor)
# Arguments: npp k eff xbin xmin xmax   (x is scanned inside the job)
#----------------------------------------------------------------------------------------------------
if ( $?GLAUBER_STARVER ) then
  starver $GLAUBER_STARVER
else
  starver SL19b
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

root4star -b <<EOF2
  .L doNbdFitMaker_OO200.C
  scan($nevents, "$real", "$mc", $multCut, 1, $npp, $npp, 1, $k, $k, $xbin, $xmin, $xmax, $eff, 1.00, kFALSE);
  .q
EOF2
