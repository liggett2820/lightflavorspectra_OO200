# This Makefile Structure Was Copied from https://medium.com/@petrstepanov/cern-root-roofit-makefile-structure-on-macos-and-linux-3305c822e1c4

# Select compiler command depending on environment
OS:=$(shell uname)
ifeq ($(OS),Darwin)
	CXX=clang++   # for mac / root 6
else
	CXX=g++       # for other systems
endif

# Variables for compiling and linking
CXXFLAGS = $(shell root-config --cflags) -fPIC -Wno-deprecated-declarations
LDFLAGS  = $(shell root-config --ldflags)
LIBS     = $(shell root-config --glibs) -lHtml -lMinuit -lFumili -lGenVector
INCDIR   = $(shell root-config --incdir)

# Define variables for directories
BASEDIR  = $(LIGHTFLAVORSPECTRA)
HEAD_DIR = $(BASEDIR)./headers
SRC_DIR  = $(BASEDIR)./source
OBJ_DIR  = $(BASEDIR)./build
BIN_DIR  = $(BASEDIR)./bin

# Set the name of your executable
APP_NAME=lightflavorspectra

# Set variables with dictionary .cxx and .pcm filenames
DICT_CXX_FILENAME=$(BIN_DIR)/$(APP_NAME)_Dict.cxx
DICT_PCM_FILENAME=$(BIN_DIR)/$(APP_NAME)_Dict_rdict.pcm

# Shared library and its debug symbols filenames
SHARED_LIBRARY=$(BIN_DIR)/$(APP_NAME)_Lib.so
SHARED_LIBRARY_DS=$(BIN_DIR)/$(APP_NAME)_Lib.so.dSYM

# Binary executable file names with path
EXECUTABLE=$(BIN_DIR)/$(APP_NAME)
EXECUTABLE_LOCAL=$(BIN_DIR)/$(APP_NAME)_local

# Define list of header files (.h), source files (.cpp) and object files (.o)
HEADERS_SIMPLE = CutClass.h SplineFitter.h  OffsetVertexAnalysis.h  ZFitter.h SpectraFitter.h RawSpectraModifier.h HybridFeeddownMaker.h StyleSettings.h ZFitterEditor.h LinEquMinimizer.h GausMixModel.h FitUtilities.h ExtraMath.h SpectraStructure.h KnockoutProtonAnalyzer.h  GraphFitter.h FitStorageHelper.h DeDxCalibration.h EfficiencyFitter.h DnDyAnalysis.h MattDstStudies.h SystematicErrorStudies.h MattMcTrack.h MattMcEvent.h MattEvent.h MattTrack.h EmbeddingBinner.h ThermusResult.h Histo2D.h
HEADERS_SIMPLE += Helix.h PhysMath.h HistogramUtilities.h PhysMathHistogramUtilities.h MattMcFDEvent.h MattMcFDTrack.h EmbeddingFeedDownBinner.h
HEADERS_SIMPLE_FULL = $(shell printf "$(HEAD_DIR)/%s " $(HEADERS_SIMPLE))
HEADERS_SIMPLE_FULL += $(BASEDIR)./submodule/ParticleInfo/ParticleInfo/ParticleInfo.h

SOURCES_SIMPLE = SpectraStructure.cxx StyleSettings.cxx CalibrationClass.cxx EfficiencyFitter.cxx CutClass.cxx  GraphFitter.cxx utilityFunctions.cxx LinEquMinimizer.cxx GausMixModel.cxx DrawSimpleAllRap.cxx SplineFitter.cxx SpectraFitter.cxx particleYields_UrQMD_root.cxx DnDyAnalysisFXT.cxx ZFitterEditor.cxx FitUtilities.cxx RawSpectraModifier.cxx ZFitter.cxx ExtraMath.cxx MattDstStudies.cxx SystematicErrorStudies.cxx HybridFeeddownMaker.cxx PionRatioFit.cxx FitStorageHelper.cxx OffsetVertexAnalysis.cxx KnockoutProtonAnalyzer.cxx DnDyAnalysis.cxx DeDxCalibration.cxx DrawEmbedQA.cxx MattEvent.cxx MattTrack.cxx MattMcTrack.cxx MattMcEvent.cxx EmbeddingBinner.cxx ThermusResult.cxx Histo2D.cxx
SOURCES_SIMPLE += Helix.cxx PhysMath.cxx HistogramUtilities.cxx PhysMathHistogramUtilities.cxx MattMcFDEvent.cxx MattMcFDTrack.cxx EmbeddingFeedDownBinner.cxx
SOURCES_SIMPLE_FULL = $(shell printf "$(SRC_DIR)/%s " $(SOURCES_SIMPLE))
SOURCES_SIMPLE_FULL += $(BASEDIR)./submodule/ParticleInfo/ParticleInfo/ParticleInfo.cxx

GLAUBER_SO = DetermineCentralityBin_cxx.so MakeResultsTable_cxx.so FindBestFitNegativeBinomialParameters_cxx.so PrintResults_cxx.so FindCentralityBinCuts_cxx.so ReturnStructs_cxx.so FindNpartNcollDistributions_cxx.so SystematicErrors_cxx.so GlauberClass_cxx.so TriggerEfficiency_cxx.so GlauberModel_cxx.so
GLAUBER_SO_FULL = $(shell printf "bin/%s " $(GLAUBER_SO))

# Define object files from source files
#OBJECTS = $(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(SOURCES_SIMPLE_FULL))
#OBJECTS = $(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(SOURCES_SIMPLE))
#OBJECTS += $(OBJ_DIR)/ParticleInfo.o

OBJECTS_SIMPLE = SpectraStructure.o StyleSettings.o CalibrationClass.o EfficiencyFitter.o CutClass.o  GraphFitter.o utilityFunctions.o LinEquMinimizer.o GausMixModel.o SplineFitter.o SpectraFitter.o DnDyAnalysisFXT.o  ZFitterEditor.o FitUtilities.o RawSpectraModifier.o ZFitter.o ExtraMath.o MattDstStudies.o SystematicErrorStudies.o HybridFeeddownMaker.o PionRatioFit.o FitStorageHelper.o OffsetVertexAnalysis.o KnockoutProtonAnalyzer.o DnDyAnalysis.o DeDxCalibration.o DrawEmbedQA.o MattEvent.o MattTrack.o MattMcTrack.o MattMcEvent.o EmbeddingBinner.o ThermusResult.o Histo2D.o ParticleInfo.o DrawSimpleAllRap.o particleYields_UrQMD_root.o
OBJECTS_SIMPLE += Helix.o PhysMath.o HistogramUtilities.o PhysMathHistogramUtilities.o MattMcFDEvent.o MattMcFDTrack.o EmbeddingFeedDownBinner.o
OBJECTS = $(shell printf "$(OBJ_DIR)/%s " $(OBJECTS_SIMPLE))


# Convenience variable for creating new directories
dir_guard=@mkdir -p $(@D)

SO_DEP=makefile

# Path for installing the shared library and dictionary .pcm
DYNAMIC_LIBRARY_PATH = ./

# Path for installing the executable binary
ifeq ($(PREFIX),)
	PREFIX := ./
endif

# Default target builds production executable
all: basic

# This is basic
basic: setup ./bin/StRefMultCorr_SL21a.so
	$(MAKE) libStPicoDst_WO_C11.so --directory=submodule/PicoDstReader/
	cp -f ./submodule/PicoDstReader/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader/libStPicoDst.so ./macros
	root -l -q -b macros/makeLibs_RCF.C

basic_SL22dev: setup ./bin/StRefMultCorr_SL22dev.so
	$(MAKE) libStPicoDst_WO_C11.so --directory=submodule/PicoDstReader/
	cp -f ./submodule/PicoDstReader/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader/libStPicoDst.so ./macros
	root -l -q -b macros/makeLibs_RCF.C

basic_SL22c: setup ./bin/StRefMultCorr_SL22dev.so
	$(MAKE) libStPicoDst.so --directory=submodule/PicoDstReader_SL22c/
	cp -f ./submodule/PicoDstReader_SL22c/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader_SL22c/libStPicoDst.so ./macros
	root -l -q -b macros/makeLibs_RCF.C

# COLLIDER options for RCF
col27: setup ./bin/StRefMultCorr_SL19b.so
	starver SL19b
	$(MAKE) libStPicoDst.so --directory=submodule/PicoDstReader_SL19b/
	cp -f ./submodule/PicoDstReader_SL19b/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader_SL19b/libStPicoDst.so ./macros/
	root -l -q -b macros/makeLibs_RCF.C

col19p6: setup ./bin/StRefMultCorr_SL21a.so
	$(MAKE) libStPicoDst_WO_C11.so --directory=submodule/PicoDstReader_SL23d/
	cp -f ./submodule/PicoDstReader_SL23d/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader_SL23d/libStPicoDst.so ./macros
	root -l -q -b macros/makeLibs_RCF.C

col17p3: setup ./bin/StRefMultCorr_9p2_11p5_17p3GeV.so
	$(MAKE) libStPicoDst_WO_C11.so --directory=submodule/PicoDstReader_SL23d/
	cp -f ./submodule/PicoDstReader_SL23d/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader_SL23d/libStPicoDst.so ./macros
	root -l -q -b macros/makeLibs_RCF.C

col14p6: setup ./bin/StRefMultCorr_SL22dev.so
	$(MAKE) libStPicoDst_WO_C11.so --directory=submodule/PicoDstReader_SL23d/
	cp -f ./submodule/PicoDstReader_SL23d/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader_SL23d/libStPicoDst.so ./macros
	root -l -q -b macros/makeLibs_RCF.C

col11p5: setup ./bin/StRefMultCorr_9p2_11p5_17p3GeV.so
	$(MAKE) libStPicoDst_WO_C11.so --directory=submodule/PicoDstReader_SL23a/
	cp -f ./submodule/PicoDstReader_SL23a/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader_SL23a/libStPicoDst.so ./macros
	root -l -q -b macros/makeLibs_RCF.C

col9p2: setup ./bin/StRefMultCorr_9p2_11p5_17p3GeV.so
	$(MAKE) libStPicoDst_WO_C11.so --directory=submodule/PicoDstReader_SL23a/
	cp -f ./submodule/PicoDstReader_SL23a/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader_SL23a/libStPicoDst.so ./macros
	root -l -q -b macros/makeLibs_RCF.C

col7p7: setup ./bin/StRefMultCorr_7p7GeV.so
	$(MAKE) libStPicoDst_WO_C11.so --directory=submodule/PicoDstReader_SL22b/
	cp -f ./submodule/PicoDstReader_SL22b/libStPicoDst.so ./bin/
	cp -f ./submodule/PicoDstReader_SL22b/libStPicoDst.so ./macros
	root -l -q -b macros/makeLibs_RCF.C

./bin/StRefMultCorr_SL19b.so:
	$(CXX) -c $(CXXFLAGS) submodule/StRoot/StRefMultCorr_SL19b/StRefMultCorr.cxx -o ./bin/StRefMultCorr.o
	$(CXX) -c $(CXXFLAGS) submodule/StRoot/StRefMultCorr_SL19b/CentralityMaker.cxx -o ./bin/CentralityMaker.o
	@echo "Compiling StRefMultCorr_SL19b.so "$@
	$(CXX) -shared -std=c++11 $(CXXFLAGS) ./bin/StRefMultCorr.o ./bin/CentralityMaker.o -o $@

./bin/StRefMultCorr_SL21a.so:
	$(CXX) -c $(CXXFLAGS) -std=c++11 submodule/StRoot/StRefMultCorr_SL21a/StRefMultCorr.cxx -o ./bin/StRefMultCorr.o
	$(CXX) -c $(CXXFLAGS) -std=c++11 submodule/StRoot/StRefMultCorr_SL21a/CentralityMaker.cxx -o ./bin/CentralityMaker.o
	@echo "Compiling StRefMultCorr_SL21a.so "$@
	$(CXX) -shared -std=c++11 $(CXXFLAGS) ./bin/StRefMultCorr.o ./bin/CentralityMaker.o -o $@

./bin/StRefMultCorr_SL22dev.so:
	$(CXX) -c $(CXXFLAGS) -std=c++11 submodule/StRoot/StRefMultCorr_SL22dev/StRefMultCorr.cxx -o ./bin/StRefMultCorr.o
	$(CXX) -c $(CXXFLAGS) -std=c++11 submodule/StRoot/StRefMultCorr_SL22dev/CentralityMaker.cxx -o ./bin/CentralityMaker.o
	@echo "Compiling StRefMultCorr_SL22dev.so "$@
	$(CXX) -shared -std=c++11 $(CXXFLAGS) ./bin/StRefMultCorr.o ./bin/CentralityMaker.o -o $@

./bin/StRefMultCorr_7p7GeV.so:
	$(CXX) -c $(CXXFLAGS) -std=c++11 submodule/StRoot/StRefMultCorr_7p7GeV/StRefMultCorr.cxx -o ./bin/StRefMultCorr.o
	$(CXX) -c $(CXXFLAGS) -std=c++11 submodule/StRoot/StRefMultCorr_7p7GeV/CentralityMaker.cxx -o ./bin/CentralityMaker.o
	@echo "Compiling StRefMultCorr_7p7GeV.so "$@
	$(CXX) -shared -std=c++11 $(CXXFLAGS) ./bin/StRefMultCorr.o ./bin/CentralityMaker.o -o $@

./bin/StRefMultCorr_9p2_11p5_17p3GeV.so:
	$(CXX) -c $(CXXFLAGS) -std=c++11 submodule/StRoot/StRefMultCorr_9p2_11p5_17p3GeV/StRefMultCorr.cxx -o ./bin/StRefMultCorr.o
	$(CXX) -c $(CXXFLAGS) -std=c++11 submodule/StRoot/StRefMultCorr_9p2_11p5_17p3GeV/CentralityMaker.cxx -o ./bin/CentralityMaker.o
	@echo "Compiling StRefMultCorr_9p2_11p5_17p3GeV.so "$@
	$(CXX) -shared -std=c++11 $(CXXFLAGS) ./bin/StRefMultCorr.o ./bin/CentralityMaker.o -o $@


setup:
	mkdir -p $(BIN_DIR)

# THE REST IS FOR ROOT6

# Add -O3 optimization level for the production release
production: CXXFLAGS+=-O3
production: directories $(DICT_CXX_FILENAME) $(SHARED_LIBRARY) $(OBJECTS) $(EXECUTABLE) move_files

# Add -O3 optimization level for the test release
release: CXXFLAGS+=-O3
release: directories $(DICT_CXX_FILENAME) $(SHARED_LIBRARY) $(OBJECTS) $(EXECUTABLE_LOCAL) move_files

# Add flags for debug optimizations: -Og -ggdb -DDEBUG
debug: CXXFLAGS+=-g -O0 -ggdb -DDEBUG
debug: directories $(DICT_CXX_FILENAME) $(OBJECTS) $(SHARED_LIBRARY)

# Target for dictionary generation (creates .cxx and .pcm)
$(DICT_CXX_FILENAME): echoheaders echo_objects
	rootcint -f $@ -c -D_VANILLA_ROOT_ -DROOT_CINT -DCLING_DEBUG=1 -D__ROOT__ -I. -I$(INCDIR) $(HEADERS_SIMPLE_FULL) ./Linkdef.h

# Target for compiling the shared library
#$(SHARED_LIBRARY): $(DICT_CXX_FILENAME) makefile
#	$(CXX) -shared -o $@ $(LDFLAGS) $(CXXFLAGS) $(LIBS) $(DICT_CXX_FILENAME) $(SOURCES_SIMPLE_FULL)
$(SHARED_LIBRARY):  $(DICT_CXX_FILENAME) $(OBJECTS) makefile
	$(CXX) -shared -o $@ $(LDFLAGS) $(CXXFLAGS) $(LIBS) $(OBJECTS) $(DICT_CXX_FILENAME)

$(OBJ_DIR)/ZFitter.o: $(SRC_DIR)/ZFitter.cxx $(SRC_DIR)/ZFitter_2DFits.cxx $(SRC_DIR)/ZFitter_BTOF_ColliderCenter.cxx $(SRC_DIR)/ZFitter_BTOF_ColliderCenter_SimulCent.cxx $(SRC_DIR)/ZFitter_bTOFFuncts.cxx $(SRC_DIR)/ZFitter_dEdxFuncts.cxx $(SRC_DIR)/ZFitter_eTOFFuncts.cxx $(SRC_DIR)/ZFitter_PionKaonTPC.cxx $(SRC_DIR)/ZFitter_PionKaonTPC_SimulCent.cxx $(SRC_DIR)/ZFitter_PionKaonTPC_SimulCent_SingleLoop.cxx $(SRC_DIR)/ZFitter_ProtonTPC.cxx $(SRC_DIR)/ZFitter_ProtonTPC_SimulCent.cxx $(SRC_DIR)/ZFitter_ProtonTPC_SimulCent_SingleLoop.cxx $(SRC_DIR)/ZFitter_TPC_FXT_SimulCent.cxx $(SRC_DIR)/ZFitter_TPCFuncts.cxx
	$(dir_guard)
	@echo "Compiling Specific "$@
	$(CXX) -c $(CXXFLAGS)  $< -o $@

$(OBJ_DIR)/ParticleInfo.o: $(SRC_DIR)/../submodule/ParticleInfo/ParticleInfo/ParticleInfo.cxx
	$(dir_guard)
	@echo "Compiling Specific "$@
	$(CXX) -c $(CXXFLAGS)  $< -o $@



# Target for compiling the object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx
	$(dir_guard)
	@echo "Compiling "$@
	$(CXX) -c $(CXXFLAGS) $< -o $@

# Target for linking the production executable
$(EXECUTABLE): $(SHARED_LIBRARY)
	@echo "Linking "$@
	$(CXX) -o $@ $(SHARED_LIBRARY) $(LIBS) $(GLAUBER_SO_FULL) ./submodule/PicoDstReader/libStPicoDst.so

# Target for linking the test release and development executable
$(EXECUTABLE_LOCAL):
ifeq ($(OS),Darwin)
	$(CXX) -o $@ $(OBJECTS) $(SHARED_LIBRARY) $(LIBS) $(GLAUBER_SO_FULL) ./submodule/PicoDstReader/libStPicoDst.so
	install_name_tool -change $(SHARED_LIBRARY) @executable_path/$(SHARED_LIBRARY) $(EXECUTABLE_LOCAL)
else
	$(CXX) -o $@ $(OBJECTS) $(SHARED_LIBRARY) $(LIBS) -Wl,-rpath,'$$ORIGIN'
endif

glauber:
	$(MAKE) --directory=submodule/glaubermcmodel/
stpicodst:
	$(MAKE) $(STPICODST) --directory=submodule/PicoDstReader/
	ln -fs $(BASEDIR)./submodule/PicoDstReader/libStPicoDst.so $(BASEDIR)/bin/
	ln -fs $(BASEDIR)./submodule/PicoDstReader/StPicoDst_Dict_rdict.pcm $(BASEDIR)/bin/
	ln -fs $(BASEDIR)./submodule/PicoDstReader/libStPicoDst.so $(BASEDIR)/
	ln -fs $(BASEDIR)./submodule/PicoDstReader/StPicoDst_Dict_rdict.pcm $(BASEDIR)/

# Target moves .so and .pcm files to ./dist
move_files:
	mv $(SHARED_LIBRARY) $(BIN_DIR)/
	mv $(DICT_PCM_FILENAME) $(BIN_DIR)/
	rm $(DICT_CXX_FILENAME)

# Target moves .so.dSYM debug symbols to /dist folder
move_debug_symbols:
ifeq ($(OS),Darwin)
	rm -fr $(BIN_DIR)/$(SHARED_LIBRARY_DS)
	mv $(SHARED_LIBRARY_DS) $(BIN_DIR)/
endif

# Target installs the executable and moves shared library to $ROOTSYS/lib
install:
	sudo install -m 755 $(EXECUTABLE) $(DESTDIR)$(PREFIX)/bin/
	sudo install -m 755 $(BIN_DIR)/$(SHARED_LIBRARY) $(DYNAMIC_LIBRARY_PATH)/
	sudo install -m 755 $(BIN_DIR)/$(DICT_PCM_FILENAME) $(DYNAMIC_LIBRARY_PATH)/

# Target to copy the application icon and create the launcher on Linux
install-linux-launcher:
	xdg-icon-resource install --context apps --size 128 ./resources/$(APP_NAME).png $(APP_NAME)
	xdg-desktop-menu install ./resources/$(APP_NAME).desktop

# Target that cleans the build
clean:
	rm -f -r $(OBJ_DIR)
	rm -f -r $(BIN_DIR)
	rm -f ./libStPicoDst.so
	rm -f ./StPicoDst_Dict_rdict.pcm
	root -l -q -b macros/makeLibs.C\(\"clean\"\)
	$(MAKE) clean --directory=submodule/PicoDstReader/
	$(MAKE) clean --directory=submodule/PicoDstReader_SL19b/
	$(MAKE) clean --directory=submodule/PicoDstReader_SL22b/
	$(MAKE) clean --directory=submodule/PicoDstReader_SL22c/
	$(MAKE) clean --directory=submodule/PicoDstReader_SL23a/
	$(MAKE) clean --directory=submodule/PicoDstReader_SL23d/
	$(MAKE) clean --directory=submodule/glaubermcmodel/

# Target for creating the ./build and ./dist directories
directories:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(BIN_DIR)

# Target for test output
echo:
	$(info SOURCES: $(SOURCES_SIMPLE_FULL))
	$(info HEADERS: $(HEADERS_SIMPLE_FULL))

echoheaders:
	$(info HEADERS: $(HEADERS_SIMPLE_FULL))

echo_objects:
	$(info OBJECTS: $(OBJECTS))

# List of special targets that do not generate files
.PHONY: clean directories move_files move_debug_symbols echo
