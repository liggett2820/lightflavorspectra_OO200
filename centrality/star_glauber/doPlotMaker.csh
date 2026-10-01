#!/bin/csh

if ( $#argv != 2 ) then
  echo ""
  echo " Usage : $0 [name] [energy]"
  echo "    Current available name's are:"
  echo "  ImpactParameter, Npart, Ncoll, Multiplicity, AreaRP, AreaPP, EccRP, EccPP"
  echo ""
  echo "  NOTE: name is case insensitive"
  echo ""
  exit
endif

# Use one STAR library version for cons and every step: setenv GLAUBER_STARVER <version> before running
if ( $?GLAUBER_STARVER ) then
  starver $GLAUBER_STARVER
else
  starver SL16d
endif

set name = "$1"
set energy = "$2"
#set mode = "$3"
set mode = "1"
set option = "-b"

root4star $option <<EOF
  .L doPlotMaker.C
  doPlotMaker("$name", "$energy", $mode);
  .q
EOF

