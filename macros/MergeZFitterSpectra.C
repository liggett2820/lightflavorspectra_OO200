// MergeZFitterSpectra.C -- lightflavorspectra_OO200
//
// Drop-in replacement for the bare `hadd` call used to merge the six ZFitter (stage 2)
// output files into stage 3's input (see scripts/run_stages_2_to_4.sh, line ~84, and any
// manual re-run of that merge step).
//
// WHY THIS EXISTS: a bare `hadd spectra_combined.root spectra_OO200_*.root` runs in a
// plain ROOT session that has never loaded this repo's compiled Histo2D dictionary.
// Histo2D (source/Histo2D.cxx) inherits from TNamed and defines no Merge(TCollection*)
// method of its own. Without its dictionary loaded, TFile can't deserialize the
// DeDx_FitData/<Species>/Cent/mTm0_ZTPC_Proj2D objects as Histo2D at all -- it falls back
// to reading them with the generic TNamed streamer, which reads too few bytes and
// silently discards the actual histogram content. That's exactly what shows up as:
//   Error in <TBufferFile::CheckByteCount>: object of class TNamed read too few bytes: 18 instead of <N>
//   Error in <TFileMerger::MergeRecursive>: TKey and object retrieve disagree on type (Histo2D vs TNamed). Continuing with TNamed.
// during a bare `hadd` run -- i.e. those three DeDx_FitData histograms end up corrupted
// (reduced to an empty TNamed) in the merged output, for every one of the six input files.
//
// FIX: load the same two libraries makeLibs_SL24y.C compiles Histo2D from --
// GausMixModel_cxx.so (headers/Histo2D.h's dependency chain via headers/ZFitter.h needs
// this first) then Histo2D_cxx.so itself -- before running the merge. Once the dictionary
// is present, TFileMerger reads Histo2D objects correctly. Histo2D has no Merge() method,
// so for a key path written identically into all six input files (as DeDx_FitData appears
// to be -- it's not species/detector-specific despite the six separate ZFitter runs),
// TFileMerger just keeps that one copy in the output rather than summing/duplicating it --
// which is the correct behavior here, not a compromise.
//
// USAGE (run from the macros/ directory, same convention as LoadRawSpectraModifierLibs.C
// and friends -- library paths below are relative to macros/):
//
//   cd /path/to/lightflavorspectra_OO200/macros
//   root -l -q -b 'MergeZFitterSpectra.C("/path/to/zfitter_output", "spectra_OO200_*_RapBin_*.root", "/path/to/spectra_combined_OO200.root")'
//
// After it runs, search the output for "Error"/"TBufferFile"/"TFileMerger" -- there
// should be none. If there still are, the object triggering them is a DIFFERENT custom
// class whose library also needs loading below.

void MergeZFitterSpectra(string a_inputDir, string a_inputGlob, string a_outputFile){

  // Histo2D's own dependency chain (see makeLibs_SL24y.C: GausMixModel.cxx is compiled
  // immediately before Histo2D.cxx because headers/ZFitter.h's #include order requires
  // GausMixModel.h to already be known before Histo2D.h is parsed).
  gSystem->Load("../bin/GausMixModel_cxx.so");
  gSystem->Load("../bin/Histo2D_cxx.so");

  TString listing = gSystem->GetFromPipe(TString::Format("ls %s/%s", a_inputDir.c_str(), a_inputGlob.c_str()));

  TFileMerger merger;
  // Force full object-level merging (not hadd's fast byte-copy path) so the Histo2D
  // dictionary we just loaded is actually exercised for every key, every file.
  merger.SetFastMethod(kFALSE);

  if (!merger.OutputFile(a_outputFile.c_str(), "RECREATE")){
    cout << "ERROR: could not open output file " << a_outputFile << endl;
    return;
  }

  int nFiles = 0;
  TObjArray* lines = listing.Tokenize("\n");
  for (int i = 0; i < lines->GetEntries(); i++){
    TString line = ((TObjString*)lines->At(i))->GetString();
    line = line.Strip(TString::kBoth);
    if (line.Length() == 0) continue;
    cout << "Adding: " << line << endl;
    if (!merger.AddFile(line)){
      cout << "ERROR: could not add file " << line << endl;
      return;
    }
    nFiles++;
  }
  delete lines;

  if (nFiles == 0){
    cout << "ERROR: no files matched " << a_inputDir << "/" << a_inputGlob << endl;
    return;
  }
  if (nFiles != 6){
    cout << "WARNING: expected 6 ZFitter output files (3 species x 2 detectors), found " << nFiles << " -- double check the glob before trusting the result." << endl;
  }

  cout << "Merging " << nFiles << " files into " << a_outputFile << " ..." << endl;
  if (!merger.Merge()){
    cout << "ERROR: TFileMerger::Merge() failed" << endl;
    return;
  }

  cout << "Done. Search the output above for ERROR/TBufferFile/TFileMerger -- there should be none now." << endl;
}
