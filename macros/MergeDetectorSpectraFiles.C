// macros/MergeDetectorSpectraFiles.C
//
// WHY THIS EXISTS (2026-07-21):
// RunZFitter.C's a_detectorIndex gating (0=TPC, 1=BTOF -- see its header comment
// and the fitPionTPC/fitPionBTOF/etc split around line 132-142) means each ZFitter
// run only fits ONE detector, and ZFitter::makeSpectra() writes only ONE output
// file per run (fresh TFile("RECREATE") every time -- ZFitter.cxx line ~3733).
// So Andrew's separate TPC and BTOF ZFitter runs produced two SEPARATE spectra
// files, each with only its own detector's data actually populated.
//
// But RawSpectraModifier::loadSpectraFile() (RawSpectraModifier.cxx line ~635)
// needs ONE file with BOTH present -- it reads DeDx_FitData/<species>/... AND
// BTOF_FitData/<species>/... (and ETOF_FitData/<species>/... if present) from the
// SAME a_filename (confirmed by direct read of loadSpectraFile(), lines 678-689).
//
// WHY NOT JUST `hadd`:
// Traced ZFitter::makeSpectra() (ZFitter.cxx ~3729-4073) top to bottom against
// RunZFitter.C's call sequence. Confirmed:
//   - makeFitDataHistos(a_partIndex, a_detectorIndex) IS gated by a_detectorIndex
//     internally (ZFitter.cxx ~2008: "if(a_detectorIndex < 0 || a_detectorIndex
//     == 0)" for TPC structures, "== 1" for BTOF) -- so DeDx_FitData/BTOF_FitData
//     end up genuinely non-overlapping between a TPC-only file and a BTOF-only
//     file: each has real histograms in exactly one of the two directories, and
//     just an empty mkdir'd directory (no objects -- HistogramUtilities::
//     ConditionalWrite() skips null pointers, confirmed via headers/
//     HistogramUtilities.h line 462-468) for the other. A plain hadd would be
//     fine for JUST these two directories.
//   - BUT fitter->fitEtaBetaGammaDeDx(a_etaBetaGammaDeDxFile,true,0.4,false) is
//     called UNCONDITIONALLY in RunZFitter.C (line ~503), with NO a_detectorIndex
//     gate at all -- so BOTH a TPC-only run and a BTOF-only run populate the
//     betagamma_FitData/ section (ZFitter.cxx ~3930) from the SAME input yield
//     data. A plain `hadd` sums histograms found at the same object path in both
//     input files -- so hadd-ing the two spectra files together would silently
//     DOUBLE every histogram under betagamma_FitData/ in the merged file.
//   - Confirmed RawSpectraModifier::loadSpectraFile() never reads anything under
//     betagamma_FitData/ (or DeDx_mom_FitData/, the other detector-agnostic
//     section) -- so this duplication wouldn't actually corrupt anything Andrew's
//     pipeline consumes. But it's silent, wrong data sitting in the merged file
//     that would mislead anyone who later looks at or tries to use that section
//     (e.g. for a future dE/dx-vs-beta*gamma cross-check), so this macro
//     sidesteps the whole question by copying ONLY the three directories
//     loadSpectraFile() actually reads (DeDx_FitData/BTOF_FitData/ETOF_FitData),
//     verbatim, from whichever input file has them.
//
// The recursive directory-copy helper below is the standard ROOT idiom for this
// (same structure as ROOT's own tutorials/io/copyFiles.C CopyDir()), adapted so
// the copied contents land directly in a caller-named destination directory
// instead of being wrapped in an extra subdirectory named after the source.
//
// Usage:
//   root -l -q 'MergeDetectorSpectraFiles.C("spectra_TPC.root","spectra_BTOF.root","","spectra_merged.root")'
// Pass "" for a_etofSpectraFile if you don't have one (ETOF fit drivers were
// never ported into this repo -- see headers/ZFitter.h and RunZFitter.C's header
// comment: "a_detectorIndex: 0=TPC, 1=BTOF (ETOF fit drivers were not ported)").
//
// STATUS: written 2026-07-21, UNTESTED (no ROOT available to Claude to run this).
// After running, sanity-check the output: TFile::ls() should show both
// DeDx_FitData/ and BTOF_FitData/ populated with real (non-empty) histograms,
// and a spot-check histogram (e.g. dEdxSpectra_PionPlus_Cent00) should have the
// same integral/entries in the merged file as in the original TPC-only file.
//
// FIX (2026-07-21, same day): the first real run of this macro (bare
// `root -l -q 'MergeDetectorSpectraFiles.C(...)'`, no lib-loading prefix) printed
// "Warning in <TClass::Init>: no dictionary for class Histo2D is available" followed
// by hundreds of "Error in <TBufferFile::CheckByteCount>: object of class TNamed
// read too few bytes" errors. Root cause: DeDx_FitData/<species>/Cent/
// mTm0_ZTPC_Proj2D/ (written from ZFitter::m_mTm0_ZTPC_by_Rapidity, ZFitter.cxx
// ~3752-3756/3789-3793) is genuinely populated with Histo2D-typed objects --
// Histo2D (headers/Histo2D.h) is a custom class (`class Histo2D : public TNamed`,
// ClassDef(Histo2D,1)) compiled+dictionary-generated via ACLiC into
// bin/Histo2D_cxx.so. Without that dictionary loaded, TClass::Init() falls back to
// an emulated class built from the file's embedded TStreamerInfo, which failed to
// deserialize Histo2D correctly here -- hence the cascade. copyDirContents()'s
// `TObject* obj = key->ReadObj();` needs the real compiled dictionary present to
// read these objects back correctly, so it's loaded explicitly below rather than
// relying on the caller to remember a LoadZFitterLibs.C prefix (this macro was
// meant to be self-contained/standalone per the Usage note above -- that's still
// true, it just needs to load what it actually depends on).
//
// Practical impact on the file already produced by the un-fixed run: the corrupted
// objects are confined to the diagnostic-only mTm0_ZTPC_Proj2D subdirectory --
// RawSpectraModifier::loadSpectraFile() never reads anything under that path (it
// reads DeDx_FitData/<species>/dEdxSpectra_..., which are plain TH2D writes, always
// available in bare ROOT, and not affected by this). Even so, regenerate the merged
// file with this fixed macro before using it -- the regeneration is cheap (a file
// copy, not a re-run of the fit) and removes any doubt.

void copyDirContents(TDirectory* a_source, TDirectory* a_dest){
  TDirectory* savdir = gDirectory;
  a_dest->cd();
  TKey* key;
  TIter nextkey(a_source->GetListOfKeys());
  while((key = (TKey*) nextkey())){
    const char* classname = key->GetClassName();
    TClass* cl = gROOT->GetClass(classname);
    if(!cl) continue;
    if(cl->InheritsFrom(TDirectory::Class())){
      a_source->cd(key->GetName());
      TDirectory* subSource = gDirectory;
      TDirectory* subDest = a_dest->mkdir(key->GetName());
      copyDirContents(subSource, subDest);
      a_dest->cd();
    } else if(cl->InheritsFrom(TTree::Class())){
      TTree* T = (TTree*) a_source->Get(key->GetName());
      a_dest->cd();
      TTree* newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      a_source->cd();
      TObject* obj = key->ReadObj();
      a_dest->cd();
      obj->Write();
      delete obj;
    }
  }
  a_dest->SaveSelf(kTRUE);
  savdir->cd();
}

void copyTopLevelDir(string a_inFileName, string a_dirName, TFile* a_outFile){
  if(a_inFileName == "") return;
  TFile* inFile = new TFile(a_inFileName.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "ERROR: MergeDetectorSpectraFiles: couldn't open " << a_inFileName << endl;
    if(inFile) delete inFile;
    return;
  }
  TDirectory* srcDir = (TDirectory*) inFile->Get(a_dirName.c_str());
  if(!srcDir){
    cout << "WARNING: MergeDetectorSpectraFiles: no " << a_dirName << " directory found in " << a_inFileName << " -- skipping" << endl;
    inFile->Close();
    delete inFile;
    return;
  }
  a_outFile->cd();
  TDirectory* destDir = a_outFile->mkdir(a_dirName.c_str());
  copyDirContents(srcDir, destDir);
  cout << "Copied " << a_dirName << " from " << a_inFileName << endl;
  inFile->Close();
  delete inFile;
}

void MergeDetectorSpectraFiles(string a_tpcSpectraFile, string a_btofSpectraFile, string a_etofSpectraFile, string a_outFileName){

  // Load the dictionary for the custom Histo2D class (see the FIX note above) so
  // copyDirContents() can correctly read+write Histo2D-typed objects (found under
  // DeDx_FitData/<species>/Cent/mTm0_ZTPC_Proj2D/) instead of falling back to a
  // broken StreamerInfo-emulated class. namespaces_cxx.so is loaded first since
  // it's the base dependency nearly everything else in this repo is built against
  // (same load order as macros/LoadZFitterLibs.C). gSystem->Load() is a harmless
  // no-op if a library is already loaded, so this is safe to run even when this
  // macro is invoked after a LoadZFitterLibs.C prefix too.
  gSystem->Load("../bin/namespaces_cxx.so");
  gSystem->Load("../bin/Histo2D_cxx.so");

  TFile* outFile = new TFile(a_outFileName.c_str(),"RECREATE");

  copyTopLevelDir(a_tpcSpectraFile,  "DeDx_FitData",  outFile);
  copyTopLevelDir(a_btofSpectraFile, "BTOF_FitData",  outFile);
  copyTopLevelDir(a_etofSpectraFile, "ETOF_FitData",  outFile);

  outFile->Write();
  outFile->Close();
  cout << "Wrote merged spectra file: " << a_outFileName << endl;
}
