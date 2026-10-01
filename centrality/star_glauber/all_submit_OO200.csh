#!/bin/csh -f
#----------------------------------------------------------------------------------------------------
# Step 1 for O+O 200 GeV: submit StFastGlauberMcMaker jobs for every systematic type
# Same as all_submit_doFastGlauberMcMaker.csh, with the system fixed to OO at 200 GeV and
# spherical nuclei (-def kFALSE; 16O has no deformation parameters in InitOO).
#
# Usage: ./all_submit_OO200.csh [Begin jobIndex] [End jobIndex] [0=test, 1=submit]
#   e.g. ./all_submit_OO200.csh 0 10 1   -> 10 jobs x 200k = 2M accepted events per type
#----------------------------------------------------------------------------------------------------
if ( $#argv != 3 ) then
  echo ""
  echo "  Usage : $0 [Begin jobIndex] [End jobIndex] [0=test, 1=submit]"
  echo ""
  exit
endif

set nevents = 200000
set begin   = "$1"
set end     = "$2"
set flag    = "$3"
set types   = ( "default" "large" "small" "largeXsec" "smallXsec" "gauss" "smallNpp" "largeNpp" )

set run = ""
if ( $flag == 1 ) then
  set run = "-run"
endif

foreach type ($types)
  ./submit_glauber.pl -sleep 0 -v $run -n $nevents -sys OO -e 200 -t $type -def kFALSE $begin $end
end
