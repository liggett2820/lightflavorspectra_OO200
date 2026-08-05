#!/bin/bash
echo Preparing Directory For New Compilation using cleanup.bash
rm ./*.session.xml
rm ./*.dataset
rm ./schedFiles/*
rm -r ./*.package
rm ./*.zip
rm ./bin_temp/*
rm ./headers/CutClass.h
rm ./headers/PhysMath.h
rm ./headers/HistogramUtilities.h
rm ./headers/makefile_toggles.h
rm ./makefile_toggles.h
rm -r ./submodule/*
rm -r ./.sl73_gcc485/

mkdir -p ./headers
mkdir -p ./bin_temp
mkdir -p ./submodule
mkdir -p ./submodule/ParticleInfo
mkdir -p ./submodule/ParticleInfo/ParticleInfo/
