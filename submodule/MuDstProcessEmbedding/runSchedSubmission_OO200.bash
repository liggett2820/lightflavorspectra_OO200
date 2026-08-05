#!/bin/bash
# runSchedSubmission_OO200.bash -- lightflavorspectra_OO200/submodule/MuDstProcessEmbedding
#
# Ported from lightflavorspectra_etof/submodule/MuDstProcessEmbedding/runSchedSubmission_Collider.bash.
# Usage: runSchedSubmission_OO200.bash <partIndex 0=pi/1=K/2=p> <charge +1/-1> <requestID> <energy tag, e.g. 200> <starver, e.g. SL23c>
#
# WHAT CHANGED FROM THE ORIGINAL, AND WHY:
#   - $OutputDir and FileListDir below are the ORIGINAL AUTHOR's (harasty's) personal
#     RCF directories. TODO: replace both with your own RCF paths -- OutputDir must
#     match the /star/data03/pwg/... path you used in muEmbed_OO200.xml's <output> URLs.
#   - Trimmed the part1==4/5/6/7 (deuteron/triton/helion/alpha) and part1==9..14
#     (Lambda/Xi/Omega/K0S/Phi) branches -- this repo's scope is pion/kaon/proton only
#     (matches ParticleInfo's 3-species convention used everywhere else in this repo).
#     Add them back the same way as the pi/K/p branches below if you need them later.
#   - XMLSCRIPTNAME points at muEmbed_OO200.xml (this repo's copy) instead of
#     muEmbed_Collider.xml.

NumFiles=all
FilePerProcess=16
XMLSCRIPTNAME=muEmbed_OO200.xml
# TODO: replace REPLACE_WITH_YOUR_RCF_USERNAME with your actual RCF username/path
FileListDir=/star/u/REPLACE_WITH_YOUR_RCF_USERNAME/lightflavorspectra_OO200/submodule/MuDstProcessEmbedding/fileLists/${4}GeV_OO200/full/
OutputDir=/star/data03/pwg/REPLACE_WITH_YOUR_RCF_USERNAME/embedding/

mkdir -p /star/data03/scratch/REPLACE_WITH_YOUR_RCF_USERNAME/tmp/muEmb

mkdir -p $OutputDir/pip/logs/
mkdir -p $OutputDir/pim/logs/
mkdir -p $OutputDir/Kp/logs/
mkdir -p $OutputDir/Km/logs/
mkdir -p $OutputDir/pro/logs/
mkdir -p $OutputDir/pbar/logs/

echo --------------------------------------------------
echo ___________  REMEMBER TO RECOMPILE with compile.bash  ________

    if [ $1 == 0 ]; then
	if [ $2 == 1 ]; then
            mkdir -p $OutputDir/pip/$4_$3
	    star-submit-template -template $XMLSCRIPTNAME -entities filelistdir=$FileListDir,part=pip,starver=$5,requestID=$3,nFiles=$NumFiles,energy=$4,filesPerProcess=$FilePerProcess,motherID=0
	elif [ $2 == -1 ]; then
            mkdir -p $OutputDir/pim/$4_$3
	    star-submit-template -template $XMLSCRIPTNAME -entities filelistdir=$FileListDir,part=pim,starver=$5,requestID=$3,nFiles=$NumFiles,energy=$4,filesPerProcess=$FilePerProcess,motherID=0
	fi
    elif [ $1 == 1 ]; then
	if [ $2 == 1 ]; then
            mkdir -p $OutputDir/Kp/$4_$3
	    star-submit-template -template $XMLSCRIPTNAME -entities filelistdir=$FileListDir,part=Kp,starver=$5,requestID=$3,nFiles=$NumFiles,energy=$4,filesPerProcess=$FilePerProcess,motherID=0
	elif [ $2 == -1 ]; then
            mkdir -p $OutputDir/Km/$4_$3
	    star-submit-template -template $XMLSCRIPTNAME -entities filelistdir=$FileListDir,part=Km,starver=$5,requestID=$3,nFiles=$NumFiles,energy=$4,filesPerProcess=$FilePerProcess,motherID=0
	fi
    elif [ $1 == 2 ]; then
	if [ $2 == 1 ]; then
            mkdir -p $OutputDir/pro/$4_$3
	    star-submit-template -template $XMLSCRIPTNAME -entities filelistdir=$FileListDir,part=pro,starver=$5,requestID=$3,nFiles=$NumFiles,energy=$4,filesPerProcess=$FilePerProcess,motherID=0
	elif [ $2 == -1 ]; then
            mkdir -p $OutputDir/pbar/$4_$3
	    star-submit-template -template $XMLSCRIPTNAME -entities filelistdir=$FileListDir,part=pbar,starver=$5,requestID=$3,nFiles=$NumFiles,energy=$4,filesPerProcess=$FilePerProcess,motherID=0
	fi
    fi

exit
