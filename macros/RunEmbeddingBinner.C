// RunEmbeddingBinner.C -- lightflavorspectra_OO200
//
// Ported from the original repo's macros/RunEmbeddingBinner_Collider.C. This is stage 2
// of the embedding pipeline: it reprocesses the TTree(s) StMuAnalysisMaker produced (at
// RCF, over your P23ic MuDst embedding sample, one .root per RequestID after hadd'ing
// each RequestID's per-job outputs together) into the TH2D-style efficiency/energy-loss
// histograms that EfficiencyFitter::loadEmbeddingFiles() reads. Runs LOCALLY (not RCF).
//
// PREREQUISITES: for each species/charge, you need the hadd'd StMuAnalysisMaker output
// file(s) for every RequestID your P23ic embedding request produced for that species/
// charge (see submodule/MuDstProcessEmbedding/README and the runBatchSpecificRequest_
// OO200GeV.bash / makeEmbedFileLists_OO200.bash TODOs for how those get produced).
//
// WHAT CHANGED FROM THE ORIGINAL, AND WHY:
//   - The original is a 7-energy-scan macro with a combineAllEnergiesExcept27 mode and
//     per-energy dca_override/z_range_override/baseDir_Eff arrays. This repo only has
//     ONE energy (O+O 200 GeV), so all of that is collapsed to a single configuration --
//     no energy index, no combine-across-energies branch.
//   - Kept the original's "combine" loop shape (one EmbeddingBinner object per species/
//     charge, looping loadEmbeddingOutput()+fillWeightedHistograms()+fillHistograms()+
//     closeInputFile() over every RequestID file for that species/charge, then ONE
//     write() at the end) rather than its "separate" per-file mode -- EfficiencyFitter::
//     loadEmbeddingFiles() takes at most two input files per species/charge (see its own
//     doc comment: "if a_inFileNameTwo is not empty, it adds both plus and minus
//     [request-ID files] together"), so if your request has more than 2 RequestIDs per
//     species/charge, combining them HERE (not at the EfficiencyFitter step) is the only
//     way to fold them all in.
//   - inputEmbeddingFiles below are PLACEHOLDERS -- I don't know your actual hadd'd
//     output filenames/paths. Fill in the real ones (however many RequestIDs your P23ic
//     request has per species/charge) before running.
//   - inputMattDstAnalysisFile (for loadDataWeightHistograms, which reweights the
//     embedding sample to match your real data's y/mTm0/DCA distribution) is left empty
//     ("") by default, meaning that reweighting step is SKIPPED -- matching the
//     original's own behavior for energies that didn't have a MattDstStudies output
//     ready yet. Point it at your real MattDstStudies output once you have one, if you
//     want the weighted efficiency (fillWeightedHistograms) to actually be weighted
//     rather than falling back to unweighted.
//   - nCentBins matches this repo's SetCutClass.C (6 bins), same convention already used
//     in RunRawSpectraModifier.C/RunEfficiencyFitter.C.

#include "../makefile_toggles.h"
#include "../headers/EmbeddingBinner.h"
#include "../macros/SetCutClass.C"

void RunEmbeddingBinner(string a_mattDstAnalysisFile = ""){

  // TODO: fill in your actual hadd'd StMuAnalysisMaker output files, one entry per
  // RequestID, for each species/charge you've processed at RCF so far. Species/charge
  // with an empty vector are skipped entirely below.
  vector<string> inputFiles_pip  = { /* "muEmbed_200_pip_REPLACE_WITH_REQUESTID_f6.root", */ };
  vector<string> inputFiles_pim  = { /* "muEmbed_200_pim_REPLACE_WITH_REQUESTID_f6.root", */ };
  vector<string> inputFiles_Kp   = { /* "muEmbed_200_Kp_REPLACE_WITH_REQUESTID_f6.root", */ };
  vector<string> inputFiles_Km   = { /* "muEmbed_200_Km_REPLACE_WITH_REQUESTID_f6.root", */ };
  vector<string> inputFiles_pro  = { /* "muEmbed_200_pro_REPLACE_WITH_REQUESTID_f6.root", */ };
  vector<string> inputFiles_pbar = { /* "muEmbed_200_pbar_REPLACE_WITH_REQUESTID_f6.root", */ };

  string outDir = "./"; // TODO: point this at wherever you want the reprocessed output written

  vector<string> partNames  = {"pip","pim","Kp","Km","pro","pbar"};
  vector< vector<string> > inputFilesBySpecies = {
    inputFiles_pip, inputFiles_pim, inputFiles_Kp, inputFiles_Km, inputFiles_pro, inputFiles_pbar
  };

  CutClass* cuts = new CutClass();
  SetCutClass(cuts);

  for(unsigned int specIndex = 0; specIndex < partNames.size(); specIndex++){
    if(inputFilesBySpecies[specIndex].empty()){
      cout << "Skipping " << partNames[specIndex] << " -- no input files given yet (see TODO above)." << endl;
      continue;
    }
    cout << "-------------------------------   PROCESSING " << partNames[specIndex] << "    ------------------------------------" << endl;
    string outputFileName = Form("%s/embed_OO200_%s_Reprocessed.root", outDir.c_str(), partNames[specIndex].c_str());

    EmbeddingBinner* obj = new EmbeddingBinner(partNames[specIndex], cuts);
    obj->makeHistograms();
    for(unsigned int fileIndex = 0; fileIndex < inputFilesBySpecies[specIndex].size(); fileIndex++){
      cout << "  Loading: " << inputFilesBySpecies[specIndex][fileIndex] << endl;
      obj->loadEmbeddingOutput(inputFilesBySpecies[specIndex][fileIndex]);
      if(a_mattDstAnalysisFile != "") obj->loadDataWeightHistograms(a_mattDstAnalysisFile, 3);
      obj->fillWeightedHistograms();
      obj->fillHistograms();
      obj->closeInputFile();
    }
    obj->write(outputFileName);
    cout << "  Wrote: " << outputFileName << endl;
    delete obj;
  }

  cout << "-------------  End of RunEmbeddingBinner  -------------" << endl;
}
