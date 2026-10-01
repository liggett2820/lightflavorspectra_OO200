#!/bin/csh

if ( $#argv != 6 ) then
  echo ""
  echo " Usage: $0 [System] [Input file list] [Output filename] [type] [re-weighting correction (true or false)] [Unit weight (true or false)]"
  echo ""
  echo ""
  exit
endif

set system = "$1"
set input  = "$2"
set output = "$3"
set type   = "$4"
set table  = "./table"
set reweighting = "$5"
set unitweight = "$6"

# Use one STAR library version for cons and every step: setenv GLAUBER_STARVER <version> before running
if ( $?GLAUBER_STARVER ) then
  starver $GLAUBER_STARVER
else
  starver SL16d
endif

root4star -b <<EOF
  .L doAnalysisMaker.C
  doAnalysisMaker("$system", "$input", "$output", "$type", "$table", $reweighting, $unitweight);
  .q
EOF

