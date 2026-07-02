// RunPicoBinner.C -- lightflavorspectra_OO200
//
// Ported from the original repo's macros/RunPicoBinner.C. The only real change: since
// _USE_CALIB_CLASS_ is never defined in this repo's makefile_toggles.h, the
// #ifdef _USE_CALIB_CLASS_ / #else branch collapses to just the #else path (calib =
// NULL) -- that's the branch the original actually took for this build too, so nothing
// about the actual behavior changes, this just removes dead branches.
//
// Run from within macros/, loading LoadPicoBinnerLibs.C FIRST in the same invocation:
//   root -l -q -b LoadPicoBinnerLibs.C 'RunPicoBinner.C("myFileList.list","out.root",9,"PION","true","false","false","false")'

#include "../makefile_toggles.h"
// Needed so cling knows these types when this macro is interpreted in a fresh ROOT
// session -- gSystem->Load() below loads the compiled implementations, but doesn't
// by itself guarantee the class declarations are known to the interpreter.
#include "../headers/CutClass.h"
#include "../headers/CalibrationClass.h"

// NOTE: this file used to also carry forward declarations for SetCutClass() and
// PicoBinner() here. Turns out they're unnecessary (and, worse, conflicting) once
// the libraries are loaded via a separate LoadPicoBinnerLibs.C step first -- cling
// picks up the genuine declarations from the loaded dictionaries automatically, and
// a second, hand-written declaration for the same function becomes an ambiguous
// duplicate candidate at the call site. Removed.

// NOTE: the gSystem->Load() calls used to be right here, but that caused two
// different problems in a row: (1) inside RunPicoBinner()'s own body, cling needs
// the whole function link-resolved before any statement in it runs, so Load()
// couldn't satisfy a call later in the SAME function; (2) moved to file scope here,
// mixing top-level executable statements with this function in one file made cling
// eagerly JIT everything (including this function's heavy std::string usage) as one
// big unit for "runStaticInitializersOnce", which tripped a libc++ ABI-tag symbol
// mismatch on this machine's toolchain. Splitting them into two separate files
// avoids both: run LoadPicoBinnerLibs.C first, THEN this file, as two arguments to
// the same `root` invocation (see the usage note at the top of this file).

void RunPicoBinner(string inFile,
                   string outFile,
                   int numParticles,           // pass 9 -- see PicoBinner.cxx header note
                   string partName,            // "PION", "KAON", "PROTON", or "ALL"
                   string makeBasicHistos,
                   string correct_T0,
                   string useBTOFToIgnoreNonParticleOfInterest,
                   string doMassSquared){
  std::cerr << "RunPicoBinner start" << endl;

  CutClass *cuts = new CutClass();
  std::cerr << " SetCutClass in RunPicoBinner" << endl;
  SetCutClass(cuts);
  std::cerr << " Finished SetCutClass in RunPicoBinner" << endl;

  cout << "Before" << endl;
  CalibrationClass* calib = NULL;
  cout << "After Calib" << endl;

  cout << "Finished Loading Calibration Class: " << calib << endl;
  int partIndex = 0;
  if(partName == "KAON"){ partIndex = 1;}
  if(partName == "PROTON"){ partIndex = 2;}
  if(partName == "ALL"){ partIndex = -1;}

  bool makeBasicHistos_bool = true;
  if(makeBasicHistos == "false"){ makeBasicHistos_bool = false;}
  bool ignoreBTOFOFNONPOI = false;
  if(useBTOFToIgnoreNonParticleOfInterest == "true"){ ignoreBTOFOFNONPOI = true;}
  bool correctT0 = false;
  if(correct_T0 == "true"){correctT0 = true;}
  bool doMassSqu = false;
  if(doMassSquared == "true"){doMassSqu = true;}

  std::cout << " Starting Call to PicoBinner in RunPicoBinner" << endl;
  PicoBinner(inFile, outFile, numParticles, cuts, partIndex, makeBasicHistos_bool, correctT0, ignoreBTOFOFNONPOI, doMassSqu, calib);
  std::cout << " Finished PicoBinner in RunPicoBinner" << endl;
  cout << "Finished PicoBinner Macro" << endl;
  gSystem->ProcessEvents();
  gSystem->Sleep(1);

  return;
}
