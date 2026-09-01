# Define compiler
CXX = g++

# Define flags. -D_VANILLA_ROOT_ is needed to avoid StMessMgr confusion
#CFLAGS = $(shell root-config --cflags) -O2 -fPIC -Wall -W -Woverloaded-virtual -Wno-deprecated-declarations
#CFLAGS += -pipe -std=c++11 -D_VANILLA_ROOT_ -I.

CFLAGS_WO_C11 = $(shell root-config --cflags) -O2 -fPIC -Wall -W -Woverloaded-virtual -Wno-deprecated-declarations  -pipe -D_VANILLA_ROOT_ -I.
# NOTE: this used to force "-std=c++11" here, appended AFTER root-config's own
# --cflags (which already specifies whatever standard your ROOT build actually
# requires, e.g. -std=c++17). Since the LAST -std= flag on the command line wins,
# that forced c++11 back on top of a ROOT install that hard-requires c++17,
# producing "ROOT requires support for C++17 or higher." Dropped it so root-config
# decides the standard, matching your actual ROOT install.
# FIXED 2026-09-01 (RCF/SL23c only -- do NOT copy this back into the SL24y sibling
# Makefile, whose whole point is letting root-config decide): SDCC's ROOT 5.34.38
# root-config --cflags supplies NO -std= flag at all (unlike a modern ROOT, which
# hard-requires e.g. c++17 and supplies it itself -- that's the case the NOTE above
# is about). Without any -std= flag, GCC 4.8.5 defaults to gnu++98, where `auto` is
# the old C++98 storage-class keyword instead of C++11 type inference -- breaking
# StPicoBEmcPidTraits.cxx's `auto getConstrainedShort = [](float x) {...}` lambda
# with "'auto' changes meaning in C++11" / "does not name a type". Re-added
# -std=c++11 here specifically, appended after root-config's own (empty-for-this-
# ROOT) --cflags, so there's nothing for it to override -- confirmed by a clean
# `make distclean && make` at RCF under `starver SL23c` after this change.
CFLAGS = $(CFLAGS_WO_C11) -std=c++11

LIBS = $(shell root-config --libs)
INCS = $(shell root-config --incdir)

# Define output library
STPICODST = libStPicoDst.so
STPICODST_WO_C11 = libStPicoDst_WO_C11.so

# Compile all *.cxx classes in the directory
SRC = $(shell find . -name "*.cxx")

all: $(STPICODST)

rcf: $(STPICODST_WO_C11)

# $(SRC:.cc=.o)
$(STPICODST): $(SRC:.cxx=.o) StPicoDst_Dict.C
	$(CXX) $(CFLAGS) -shared $^ -o $(STPICODST) $(LIBS)

$(STPICODST_WO_C11): $(SRC:.cxx=.o) StPicoDst_Dict.C
	$(CXX) $(CFLAGS_WO_C11) -shared $^ -o $(STPICODST_WO_C11) $(LIBS)
	mv $(STPICODST_WO_C11) $(STPICODST)


%.o: %.cxx
	$(CXX) -fPIC $(CFLAGS) -c -o $@ $<

StPicoDst_Dict.C: $(shell find . -name "*.h" ! -name "*LinkDef*")
	rootcint -f $@ -c -D_VANILLA_ROOT_ -DROOT_CINT -D__ROOT__ -I. -I$(INCS) $^ StPicoDstLinkDef.h

.PHONY: clean distclean

clean:
	rm -vf *.o StarRoot/*.o StPicoDst_Dict*

distclean:
	rm -vf *.o StarRoot/*.o StPicoDst_Dict* $(STPICODST)
