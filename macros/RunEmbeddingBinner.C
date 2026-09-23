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
void SetCutClass(CutClass* a_cuts); // impl loaded via LoadEmbeddingBinnerLibs.C

void RunEmbeddingBinner(string a_mattDstAnalysisFile = ""){

  // TODO: fill in your actual hadd'd StMuAnalysisMaker output files, one entry per
  // RequestID, for each species/charge you've processed at RCF so far. Species/charge
  // with an empty vector are skipped entirely below.
  vector<string> inputFiles_pip;
  inputFiles_pip.push_back("/star/data03/pwg/liggett2820/embedding/pip/200_20232003/muEmbed_pip_200_20232003_merged.root");
  vector<string> inputFiles_pim;
  inputFiles_pim.push_back("/star/data03/pwg/liggett2820/embedding/pim/200_20232003/muEmbed_pim_200_20232003_merged.root");
  vector<string> inputFiles_Kp;
  inputFiles_Kp.push_back("/star/data03/pwg/liggett2820/embedding/Kp/200_20232003/muEmbed_Kp_200_20232003_merged.root");
  vector<string> inputFiles_Km;
  inputFiles_Km.push_back("/star/data03/pwg/liggett2820/embedding/Km/200_20232003/muEmbed_Km_200_20232003_merged.root");
  vector<string> inputFiles_pro;
  inputFiles_pro.push_back("/star/data03/pwg/liggett2820/embedding/pro/200_20232003/muEmbed_pro_200_20232003_merged.root");
  vector<string> inputFiles_pbar;
  inputFiles_pbar.push_back("/star/data03/pwg/liggett2820/embedding/pbar/200_20232003/muEmbed_pbar_200_20232003_merged.root");

  string outDir = "/star/u/liggett2820/lightflavorspectra_OO200/embedding_reprocessed";

  vector<string> partNames;
    partNames.push_back("pip");
    partNames.push_back("pim");
    partNames.push_back("Kp");
    partNames.push_back("Km");
    partNames.push_back("pro");
    partNames.push_back("pbar");
  vector< vector<string> > inputFilesBySpecies;
    inputFilesBySpecies.push_back(inputFiles_pip);
    inputFilesBySpecies.push_back(inputFiles_pim);
    inputFilesBySpecies.push_back(inputFiles_Kp);
    inputFilesBySpecies.push_back(inputFiles_Km);
    inputFilesBySpecies.push_back(inputFiles_pro);
    inputFilesBySpecies.push_back(inputFiles_pbar);

  CutClass* cuts = new CutClass();
  SetCutClass(cuts);
  const char* envSpeciesC = gSystem->Getenv("EMBED_SPECIES");
  const char* envInputFileC = gSystem->Getenv("EMBED_INPUTFILE");
  string envSpecies = (envSpeciesC != 0) ? envSpeciesC : "";
  string envInputFile = (envInputFileC != 0) ? envInputFileC : "";
  if(envSpecies != "" && envInputFile != ""){
    string outPartialDir = outDir + "/partials";
    string base = envInputFile.substr(envInputFile.find_last_of("/") + 1);
    size_t dotPos = base.find_last_of(".");
    if(dotPos != string::npos) base = base.substr(0, dotPos);
    string outputFileName = Form("%s/embed_OO200_%s_%s_partial.root", outPartialDir.c_str(), envSpecies.c_str(), base.c_str());
    cout << "-------------------------------   PROCESSING (per-file) " << envSpecies << "    ------------------------------------" << endl;
    EmbeddingBinner* obj = new EmbeddingBinner(envSpecies, cuts);
    obj->makeHistograms();
    cout << "  Loading: " << envInputFile << endl;
    obj->loadEmbeddingOutput(envInputFile);
    if(a_mattDstAnalysisFile != "") obj->loadDataWeightHistograms(a_mattDstAnalysisFile, 3);
    obj->fillHistograms();
    obj->fillWeightedHistograms();
    obj->closeInputFile();
    obj->write(outputFileName);
    cout << "  Wrote: " << outputFileName << endl;
    delete obj;
    cout << "-------------  End of RunEmbeddingBinner (per-file)  -------------" << endl;
    return;
  }

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
      obj->fillHistograms();
      obj->fillWeightedHistograms();
      obj->closeInputFile();
    }
    obj->write(outputFileName);
    cout << "  Wrote: " << outputFileName << endl;
    delete obj;
  }

  cout << "-------------  End of RunEmbeddingBinner  -------------" << endl;
}
