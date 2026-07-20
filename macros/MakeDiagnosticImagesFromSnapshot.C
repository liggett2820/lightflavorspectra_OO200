// MakeDiagnosticImagesFromSnapshot.C -- lightflavorspectra_OO200, added 2026-07
//
// Turns a diagnostic snapshot file (written by a ZFitter run with
// setSaveDiagnosticSnapshots(true, ...) on -- see RunZFitter.C's saveDiagnosticSnapshots
// toggle) into the same per-bin diagnostic PNGs a live setSaveDiagnosticImages(true) run
// would have produced, WITHOUT re-running the fit. This is the deferred half of the
// "run fast now, make pictures later" split -- see source/ZFitter_DiagnosticSnapshot.cxx
// for the shared rendering code both this macro and the live fit drivers call, and that
// file's header comment for the one thing this does NOT cover (ProtonTPC's SimulCent
// all-centralities overlay plot -- stays live-only, see that comment for why).
//
// USAGE (run from within macros/, loading LoadZFitterLibs.C first in the same
// invocation, matching RunZFitter.C's own pattern):
//
//   root -l -b -q LoadZFitterLibs.C 'MakeDiagnosticImagesFromSnapshot.C("/path/to/diagnosticSnapshots_yieldFits_OO200_Pion_TPC_AllRap.root","/path/to/outputDir/","yieldFits_OO200_Pion_TPC_AllRap")'
//
// a_snapshotFile: the diagnosticSnapshots_*.root file a prior run wrote (see
//                 RunZFitter.C's diagnosticSnapshotFileName)
// a_outputDir:    same a_outputDir you gave that RunZFitter.C call
// a_imageDirName: same imageDirName that run used -- RunZFitter.C prints this
//                 ("image dir name: ...") when it starts, or reconstruct it from
//                 partNames[a_partIndex]/detNames[a_detectorIndex] the same way
//                 RunZFitter.C does

#include <set>
#include "../makefile_toggles.h"
#include "../headers/CutClass.h"
#include "../headers/ZFitter.h"

void MakeDiagnosticImagesFromSnapshot(string a_snapshotFile, string a_outputDir, string a_imageDirName){

  TFile* snapshotFile = new TFile(a_snapshotFile.c_str(),"READ");
  if(snapshotFile == nullptr || snapshotFile->IsZombie()){
    cout << "ERROR: could not open snapshot file " << a_snapshotFile << endl;
    return;
  }
  TTree* tree = (TTree*) snapshotFile->Get("DiagnosticSnapshots");
  if(tree == nullptr){
    cout << "ERROR: no DiagnosticSnapshots TTree found in " << a_snapshotFile << endl;
    return;
  }

  // Same CutClass setup RunZFitter.C does before constructing a ZFitter -- cheap
  // insurance against relying on any not-otherwise-initialized cut/centrality-binning
  // state inside the ZFitter constructor or setImageDir() below, even though the
  // rendering methods themselves only read m_imagePreDir/m_imgDirName (set via
  // setImageDir just below) and nothing from any loaded yield data, so no
  // loadDataHistograms(...)/loadQuickStartFile(...) call is needed here at all.
  CutClass* cuts = new CutClass();
  SetCutClass(cuts);
  ZFitter* fitter = new ZFitter();
  fitter->setCuts(cuts);
  // a_detectorIndex left at its <0 default here on purpose: this macro doesn't know (or
  // need to know) which detector the snapshot is for ahead of time -- driverType on each
  // row tells us that already, and setImageDir()'s <0 fallback just creates the
  // superset of subdirectories, which is harmless (see that method's own comment for
  // the fuller "only creates what's live" explanation this predates).
  fitter->setImageDir(a_outputDir, a_imageDirName, -1);

  int driverType, charge, centIndex, rapBin, mtm0Bin, fittingRound;
  double lowRange, highRange, currentLowFitRange, currentHighFitRange, sigmaGuess;
  bool drawInitialSeeds, saveNoLogImages;
  double N[3], N_err[3], mu[3], mu_err[3], sigma[3], sigma_err[3], gamma[3], gamma_err[3];
  double seedN[3], seedMu[3], seedSigma[3], seedGamma[3];
  int doFit[3], color[3];
  string* particleName = new string();
  string* roundTag = new string();
  vector<string>* topTextLines = new vector<string>();
  vector<string>* bottomTextLines = new vector<string>();
  TH1D* hist = nullptr;

  tree->SetBranchAddress("driverType", &driverType);
  tree->SetBranchAddress("particleName", &particleName);
  tree->SetBranchAddress("charge", &charge);
  tree->SetBranchAddress("centIndex", &centIndex);
  tree->SetBranchAddress("rapBin", &rapBin);
  tree->SetBranchAddress("mtm0Bin", &mtm0Bin);
  tree->SetBranchAddress("fittingRound", &fittingRound);
  tree->SetBranchAddress("roundTag", &roundTag);
  tree->SetBranchAddress("lowRange", &lowRange);
  tree->SetBranchAddress("highRange", &highRange);
  tree->SetBranchAddress("currentLowFitRange", &currentLowFitRange);
  tree->SetBranchAddress("currentHighFitRange", &currentHighFitRange);
  tree->SetBranchAddress("sigmaGuess", &sigmaGuess);
  tree->SetBranchAddress("drawInitialSeeds", &drawInitialSeeds);
  tree->SetBranchAddress("saveNoLogImages", &saveNoLogImages);
  tree->SetBranchAddress("N", N);
  tree->SetBranchAddress("N_err", N_err);
  tree->SetBranchAddress("mu", mu);
  tree->SetBranchAddress("mu_err", mu_err);
  tree->SetBranchAddress("sigma", sigma);
  tree->SetBranchAddress("sigma_err", sigma_err);
  tree->SetBranchAddress("gamma", gamma);
  tree->SetBranchAddress("gamma_err", gamma_err);
  tree->SetBranchAddress("seedN", seedN);
  tree->SetBranchAddress("seedMu", seedMu);
  tree->SetBranchAddress("seedSigma", seedSigma);
  tree->SetBranchAddress("seedGamma", seedGamma);
  tree->SetBranchAddress("doFit", doFit);
  tree->SetBranchAddress("color", color);
  tree->SetBranchAddress("topTextLines", &topTextLines);
  tree->SetBranchAddress("bottomTextLines", &bottomTextLines);
  tree->SetBranchAddress("hist", &hist);

  Long64_t nEntries = tree->GetEntries();
  cout << "MakeDiagnosticImagesFromSnapshot: " << nEntries << " snapshot rows in " << a_snapshotFile << endl;

  // Directories are created up front per (particleName,detector-flavor) combination
  // actually present in the file, rather than relying on ZFitter::setImageDir()'s
  // m_particleLoaded gating (which this macro never populates, since it never loads any
  // yield data) -- straightforward gSystem->mkdir(...,kTRUE) instead.
  set<string> dirsMade;

  for(Long64_t iii = 0; iii < nEntries; iii++){
    tree->GetEntry(iii);

    ZFitterDiagnosticSnapshot snapshot;
    snapshot.driverType   = driverType;
    snapshot.particleName = *particleName;
    snapshot.charge       = charge;
    snapshot.centIndex    = centIndex;
    snapshot.rapBin       = rapBin;
    snapshot.mtm0Bin      = mtm0Bin;
    snapshot.fittingRound = fittingRound;
    snapshot.roundTag     = *roundTag;
    snapshot.lowRange            = lowRange;
    snapshot.highRange           = highRange;
    snapshot.currentLowFitRange  = currentLowFitRange;
    snapshot.currentHighFitRange = currentHighFitRange;
    snapshot.sigmaGuess          = sigmaGuess;
    snapshot.drawInitialSeeds    = drawInitialSeeds;
    snapshot.saveNoLogImages     = saveNoLogImages;
    for(int sp = 0; sp < 3; sp++){
      snapshot.N[sp] = N[sp]; snapshot.N_err[sp] = N_err[sp];
      snapshot.mu[sp] = mu[sp]; snapshot.mu_err[sp] = mu_err[sp];
      snapshot.sigma[sp] = sigma[sp]; snapshot.sigma_err[sp] = sigma_err[sp];
      snapshot.gamma[sp] = gamma[sp]; snapshot.gamma_err[sp] = gamma_err[sp];
      snapshot.seedN[sp] = seedN[sp]; snapshot.seedMu[sp] = seedMu[sp];
      snapshot.seedSigma[sp] = seedSigma[sp]; snapshot.seedGamma[sp] = seedGamma[sp];
      snapshot.doFit[sp] = (doFit[sp] != 0);
      snapshot.color[sp] = color[sp];
    }
    snapshot.topTextLines    = *topTextLines;
    snapshot.bottomTextLines = *bottomTextLines;
    snapshot.hist = hist;

    string dirKey = Form("%s_%d", snapshot.particleName.c_str(), snapshot.centIndex);
    if(dirsMade.find(dirKey) == dirsMade.end()){
      gSystem->mkdir(Form("%s/%s/%s/dEdxFits_Cent%02d",a_outputDir.c_str(), a_imageDirName.c_str(), snapshot.particleName.c_str(), snapshot.centIndex), kTRUE);
      dirsMade.insert(dirKey);
    }

    if(snapshot.driverType == 0)      fitter->renderPionKaonTPCDiagnosticPlot(snapshot);
    else if(snapshot.driverType == 1) fitter->renderProtonTPCDiagnosticPlot(snapshot);
    else cout << "WARNING: unrecognized driverType " << snapshot.driverType << " on row " << iii << " -- skipped" << endl;

    if(iii % 500 == 0) cout << "  " << iii << " / " << nEntries << endl;
  }

  snapshotFile->Close();
  cout << "MakeDiagnosticImagesFromSnapshot: done." << endl;
}
