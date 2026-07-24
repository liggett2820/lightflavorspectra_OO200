// ZFitter_DiagnosticSnapshot.cxx -- lightflavorspectra_OO200, added 2026-07
//
// WHY THIS FILE EXISTS:
//
// The per-bin diagnostic fit-overlay PNGs (gated by m_saveDiagnosticImages) are by far
// the most expensive part of a ZFitter run -- a full TCanvas render + PNG encode for
// every rapidity/mTm0/centrality/round combination. Andrew asked to be able to run the
// fit fast (skip image rendering) and generate the diagnostic images afterward, in a
// second pass, from a separate macro -- not just move this code to a new file while
// still running it inline every time.
//
// That requires two things this file provides:
//
//   1. A snapshot format (see ZFitterDiagnosticSnapshot in headers/ZFitter.h) that
//      captures everything the plotting code reads out of live fit state -- the raw
//      per-bin histogram, the final fit parameters/errors, the initial seed values used
//      for the optional "_init" overlay curves, and the exact pre-formatted text-box
//      contents (captured while the minimizer/paramNames/limits are still valid, rather
//      than re-derived later from raw numbers, to avoid any risk of drifting from what a
//      live run would have shown) -- plus the machinery to write/read that as a TTree
//      (ensureDiagnosticSnapshotTreeOpen / fillDiagnosticSnapshotTreeRow / closeDiagnosticSnapshotFile,
//      all below).
//
//   2. Two shared rendering functions (renderPionKaonTPCDiagnosticPlot /
//      renderProtonTPCDiagnosticPlot) that do the ACTUAL canvas/TF1/text-box/pull-
//      histogram/SaveAs work, parameterized off a ZFitterDiagnosticSnapshot instead of
//      the enclosing fit function's locals. These are called from both places:
//        - live, from within ZFitter_PionKaonTPC_SimulCent_SingleLoop.cxx /
//          ZFitter_ProtonTPC_SimulCent_SingleLoop.cxx, when m_saveDiagnosticImages is on
//          (same visual output as before this change -- the code was moved here nearly
//          verbatim, just reading off the snapshot struct instead of local variables)
//        - deferred, from macros/MakeDiagnosticImagesFromSnapshot.C, reading back a
//          snapshot file written by a previous m_saveDiagnosticSnapshots=true run
//      so the two driver files' plotting logic can't silently diverge from what the
//      replay macro produces -- there is exactly one copy of this code, not two.
//
// SCOPE LIMITATION (flagged to Andrew when this was designed, not silently dropped):
// ProtonTPC's one-per-bin "SimulCent_RapBin_...png" overlay (all 6 centralities drawn
// together on one canvas -- see ZFitter_ProtonTPC_SimulCent_SingleLoop.cxx, the
// "DRAW IMAGE WITH ALL CENTRALITIES" block) is NOT covered here. It draws from ALL
// centralities' histograms/fits at once, which doesn't fit this snapshot's one-row-per-
// (bin x centrality) shape -- and it's 6x cheaper anyway (one per bin, not one per bin x
// centrality), so it's a small slice of the actual cost problem. That plot stays gated on
// m_saveDiagnosticImages only, unchanged, in its original location.
//
// #included directly into ZFitter.cxx's translation unit, same pattern as the other
// ZFitter_*.cxx driver files (see that file's #include list) -- so this compiles into
// the same library, with no ABI-mismatch risk between separately-built .so's.

//############################################################################
// TREE WRITE SIDE -- called from the live fit drivers
//############################################################################

void ZFitter::ensureDiagnosticSnapshotTreeOpen(){
  if(m_diagnosticSnapshotTree != nullptr) return; // already open

  m_diagnosticSnapshotFile = new TFile(m_diagnosticSnapshotFilePath.c_str(),"RECREATE");
  m_diagnosticSnapshotTree = new TTree("DiagnosticSnapshots","ZFitter per-bin diagnostic plot snapshots (2026-07)");

  ZFitterDiagnosticSnapshot& row = m_diagnosticSnapshotRow;
  m_diagnosticSnapshotTree->Branch("driverType",    &row.driverType,    "driverType/I");
  m_diagnosticSnapshotTree->Branch("particleName",  &row.particleName);
  m_diagnosticSnapshotTree->Branch("charge",        &row.charge,        "charge/I");
  m_diagnosticSnapshotTree->Branch("centIndex",     &row.centIndex,     "centIndex/I");
  m_diagnosticSnapshotTree->Branch("rapBin",        &row.rapBin,        "rapBin/I");
  m_diagnosticSnapshotTree->Branch("mtm0Bin",       &row.mtm0Bin,       "mtm0Bin/I");
  m_diagnosticSnapshotTree->Branch("fittingRound",  &row.fittingRound,  "fittingRound/I");
  m_diagnosticSnapshotTree->Branch("roundTag",      &row.roundTag);

  m_diagnosticSnapshotTree->Branch("lowRange",             &row.lowRange,             "lowRange/D");
  m_diagnosticSnapshotTree->Branch("highRange",             &row.highRange,             "highRange/D");
  m_diagnosticSnapshotTree->Branch("currentLowFitRange",   &row.currentLowFitRange,   "currentLowFitRange/D");
  m_diagnosticSnapshotTree->Branch("currentHighFitRange",  &row.currentHighFitRange,  "currentHighFitRange/D");
  m_diagnosticSnapshotTree->Branch("sigmaGuess",            &row.sigmaGuess,            "sigmaGuess/D");
  m_diagnosticSnapshotTree->Branch("drawInitialSeeds",      &row.drawInitialSeeds,      "drawInitialSeeds/O");
  m_diagnosticSnapshotTree->Branch("saveNoLogImages",       &row.saveNoLogImages,       "saveNoLogImages/O");

  m_diagnosticSnapshotTree->Branch("N",         row.N,         "N[3]/D");
  m_diagnosticSnapshotTree->Branch("N_err",     row.N_err,     "N_err[3]/D");
  m_diagnosticSnapshotTree->Branch("mu",        row.mu,        "mu[3]/D");
  m_diagnosticSnapshotTree->Branch("mu_err",    row.mu_err,    "mu_err[3]/D");
  m_diagnosticSnapshotTree->Branch("sigma",     row.sigma,     "sigma[3]/D");
  m_diagnosticSnapshotTree->Branch("sigma_err", row.sigma_err, "sigma_err[3]/D");
  m_diagnosticSnapshotTree->Branch("gamma",     row.gamma,     "gamma[3]/D");
  m_diagnosticSnapshotTree->Branch("gamma_err", row.gamma_err, "gamma_err[3]/D");

  m_diagnosticSnapshotTree->Branch("seedN",     row.seedN,     "seedN[3]/D");
  m_diagnosticSnapshotTree->Branch("seedMu",    row.seedMu,    "seedMu[3]/D");
  m_diagnosticSnapshotTree->Branch("seedSigma", row.seedSigma, "seedSigma[3]/D");
  m_diagnosticSnapshotTree->Branch("seedGamma", row.seedGamma, "seedGamma[3]/D");
  // bool[3] doesn't get a clean ROOT leaflist type -- bind to the stable-address int[3]
  // shadow instead (kept in sync with row.doFit in fillDiagnosticSnapshotTreeRow below).
  m_diagnosticSnapshotTree->Branch("doFit", m_diagnosticSnapshotDoFitInt, "doFit[3]/I");
  m_diagnosticSnapshotTree->Branch("color", row.color, "color[3]/I");

  m_diagnosticSnapshotTree->Branch("topTextLines",    &row.topTextLines);
  m_diagnosticSnapshotTree->Branch("bottomTextLines", &row.bottomTextLines);

  m_diagnosticSnapshotHistPtr = nullptr;
  m_diagnosticSnapshotTree->Branch("hist","TH1D",&m_diagnosticSnapshotHistPtr);
}

void ZFitter::fillDiagnosticSnapshotTreeRow(const ZFitterDiagnosticSnapshot& a_snapshot){
  ensureDiagnosticSnapshotTreeOpen();

  // Copy every field except the bool[3] arrays (handled via the static int[3] shadow
  // branches booked in ensureDiagnosticSnapshotTreeOpen() above -- TTree doesn't have a
  // clean native bool-array leaflist) and hist (handled separately below, since it needs
  // a uniquely-named clone rather than a straight pointer copy).
  m_diagnosticSnapshotRow = a_snapshot;
  m_diagnosticSnapshotRow.hist = nullptr; // don't let the row hold a dangling reference to the caller's live histogram
  for(int iii = 0; iii < 3; iii++) m_diagnosticSnapshotDoFitInt[iii] = a_snapshot.doFit[iii] ? 1 : 0;

  // Unique clone name so ROOT doesn't warn about replacing an existing object of the
  // same name across the many Fill() calls in one run.
  TString uniqueHistName = Form("snapHist_d%d_%s_c%d_cent%d_r%03d_m%04d_s%d",
    a_snapshot.driverType, a_snapshot.particleName.c_str(), a_snapshot.charge,
    a_snapshot.centIndex, a_snapshot.rapBin, a_snapshot.mtm0Bin, a_snapshot.fittingRound);
  m_diagnosticSnapshotHistPtr = (a_snapshot.hist != nullptr) ? (TH1D*) a_snapshot.hist->Clone(uniqueHistName) : nullptr;

  m_diagnosticSnapshotTree->Fill();

  if(m_diagnosticSnapshotHistPtr != nullptr){
    delete m_diagnosticSnapshotHistPtr;
    m_diagnosticSnapshotHistPtr = nullptr;
  }

  // Durability: this codebase's fit drivers have segfaulted mid-run more than once this
  // analysis cycle (see the project log). AutoSave periodically so a crash doesn't lose
  // every snapshot written so far, not just the ones from a clean exit.
  static int fillsSinceAutoSave = 0;
  fillsSinceAutoSave++;
  if(fillsSinceAutoSave >= 50){
    m_diagnosticSnapshotTree->AutoSave("SaveSelf");
    fillsSinceAutoSave = 0;
  }
}

void ZFitter::closeDiagnosticSnapshotFile(){
  if(m_diagnosticSnapshotTree != nullptr){
    // BUGFIX (2026-07-21): TObject::Write() writes into gDirectory, NOT necessarily the
    // tree's own file -- by the time this runs, makeSpectra() has already closed the main
    // spectra outFile, which resets gDirectory back to gROOT ("Rint"). Writing without
    // first cd()-ing into m_diagnosticSnapshotFile silently failed with
    // "Error in <TROOT::WriteTObject>: The current directory (Rint) is not associated
    // with a file. The object () has not been written." -- confirmed via the 2026-07-21
    // proton TPC run's log. AutoSave() during the fit loop (see fillDiagnosticSnapshotTreeRow
    // above) still correctly targeted the tree's own file, since gDirectory hadn't been
    // clobbered yet at that point -- so this bug only risked losing (a) entries filled
    // since the last AutoSave threshold and (b) a clean final tree header, not the whole
    // snapshot file. cd() into the right file first so this Write() actually lands there.
    if(m_diagnosticSnapshotFile != nullptr) m_diagnosticSnapshotFile->cd();
    m_diagnosticSnapshotTree->Write("", TObject::kOverwrite);
  }
  if(m_diagnosticSnapshotFile != nullptr){
    m_diagnosticSnapshotFile->Close();
    delete m_diagnosticSnapshotFile;
  }
  m_diagnosticSnapshotFile = nullptr;
  m_diagnosticSnapshotTree = nullptr;
}

//############################################################################
// RENDER SIDE -- called live (m_saveDiagnosticImages) and from the replay macro
//############################################################################

void ZFitter::renderPionKaonTPCDiagnosticPlot(const ZFitterDiagnosticSnapshot& a_snapshot){
  if(a_snapshot.hist == nullptr) return;

  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  fittingCanvas->Draw();
  TPad* topPad    = new TPad("topPad_PiKTPC_snap","topPad_PiKTPC_snap",0.0,0.3,1.0,1.0);
  TPad* bottomPad = new TPad("bottomPad_PiKTPC_snap","bottomPad_PiKTPC_snap",0.0,0.0,1.0,0.3);
  topPad->Draw();
  bottomPad->Draw();
  fittingCanvas->cd();
  topPad->cd();
  topPad->SetLogy(true);
  topPad->SetRightMargin(0.4);
  bottomPad->SetLogy(false);
  bottomPad->SetRightMargin(0.4);
  gStyle->SetOptFit(0111);

  TH1D* hist = a_snapshot.hist;
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerSize(0.5);
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlack);

  TString nameBase = Form("snap_%s_c%d_cent%d_r%03d_m%04d_s%d",
    a_snapshot.particleName.c_str(), a_snapshot.charge, a_snapshot.centIndex,
    a_snapshot.rapBin, a_snapshot.mtm0Bin, a_snapshot.fittingRound);

  #ifndef _ZFITTER_PHYSMATH_BYPASS_
    TF1* combinedFunct = PhysMath::getSkewAndTwoGausFunct(Form("%s_combined",nameBase.Data()),
      a_snapshot.N[0], a_snapshot.mu[0], a_snapshot.sigma[0], a_snapshot.gamma[0],
      a_snapshot.N[1], a_snapshot.mu[1], a_snapshot.sigma[1],
      a_snapshot.N[2], a_snapshot.mu[2], a_snapshot.sigma[2]);
    TF1* pionFunct  = PhysMath::skewNormalFunct(a_snapshot.N[0], a_snapshot.mu[0], a_snapshot.sigma[0], a_snapshot.gamma[0]);
    TF1* kaonFunct  = PhysMath::getGaussianFunction(Form("%s_kaon",nameBase.Data()), a_snapshot.N[1], a_snapshot.mu[1], a_snapshot.sigma[1]);
    TF1* electronFunct = PhysMath::getGaussianFunction(Form("%s_electron",nameBase.Data()), a_snapshot.N[2], a_snapshot.mu[2], a_snapshot.sigma[2]);
    pionFunct->SetLineColor(a_snapshot.color[0]);
    pionFunct->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
    electronFunct->SetLineColor(a_snapshot.color[2]);
    electronFunct->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
    kaonFunct->SetLineColor(a_snapshot.color[1]);
    kaonFunct->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
    TF1* pionFunct_init = NULL;
    TF1* koanFunct_init = NULL;
    TF1* electronFunct_init = NULL;
    if(a_snapshot.drawInitialSeeds){
      pionFunct_init  = PhysMath::skewNormalFunct(a_snapshot.seedN[0], a_snapshot.seedMu[0], a_snapshot.seedSigma[0], a_snapshot.seedGamma[0]);
      koanFunct_init  = PhysMath::getGaussianFunction(Form("%s_kaon_init",nameBase.Data()), a_snapshot.seedN[1], a_snapshot.seedMu[1], a_snapshot.seedSigma[1]);
      electronFunct_init = PhysMath::getGaussianFunction(Form("%s_electron_init",nameBase.Data()), a_snapshot.seedN[2], a_snapshot.seedMu[2], a_snapshot.seedSigma[2]);
      pionFunct_init->SetLineColor(kYellow+2);
      pionFunct_init->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
      electronFunct_init->SetLineColor(kOrange);
      electronFunct_init->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
      koanFunct_init->SetLineColor(kYellow+1);
      koanFunct_init->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
    }
  #endif

  double highestValue  = HistogramUtilities::getMaxInRange(hist, a_snapshot.currentLowFitRange, a_snapshot.currentHighFitRange);
  double lowestValue   = HistogramUtilities::getMinInRange(hist, a_snapshot.currentLowFitRange, a_snapshot.currentHighFitRange);
  TLine* lowRangeLine  = new TLine(a_snapshot.currentLowFitRange, lowestValue, a_snapshot.currentLowFitRange, highestValue);
  TLine* highRangeLine = new TLine(a_snapshot.currentHighFitRange, lowestValue, a_snapshot.currentHighFitRange, highestValue);
  lowRangeLine->SetLineColor(kViolet);
  highRangeLine->SetLineColor(kViolet);

  topPad->cd();
  HistogramUtilities::setAxisRanges(hist, a_snapshot.lowRange - 0.5*a_snapshot.sigmaGuess, a_snapshot.highRange + 0.5*a_snapshot.sigmaGuess, true);
  topPad->SetLogy(true);
  hist->SetStats(kFALSE);
  hist->Draw("E");
  #ifndef _ZFITTER_PHYSMATH_BYPASS_
    if(pionFunct_init && a_snapshot.doFit[0]) pionFunct_init->Draw("same");
    if(koanFunct_init && a_snapshot.doFit[1]) koanFunct_init->Draw("same");
    if(electronFunct_init && a_snapshot.doFit[2]) electronFunct_init->Draw("same");
    combinedFunct->Draw("same");
    pionFunct->Draw("same");
    electronFunct->Draw("same");
    kaonFunct->Draw("same");
  #endif
  lowRangeLine->Draw("same");
  highRangeLine->Draw("same");

  TPaveText* fitParsTxt = new TPaveText(0.65,0.05,0.95,0.90,"NDC");
  for(unsigned int iii = 0; iii < a_snapshot.topTextLines.size(); iii++) fitParsTxt->AddText(a_snapshot.topTextLines[iii].c_str());
  fitParsTxt->SetFillStyle(0);
  fitParsTxt->SetBorderSize(5);
  fitParsTxt->Draw("same");

  bottomPad->cd();
  #ifndef _ZFITTER_PHYSMATH_BYPASS_
  TH1D* pullHisto = HistogramUtilities::makeFitPullHistogram(hist, combinedFunct, true);
  pullHisto->SetMarkerStyle(20);
  pullHisto->SetMarkerSize(0.5);
  pullHisto->SetMarkerColor(kBlack);
  pullHisto->SetLineColor(kBlack);
  #endif

  TPaveText* fitTxt = new TPaveText(0.65,0.05,0.95,0.95,"NDC");
  for(unsigned int iii = 0; iii < a_snapshot.bottomTextLines.size(); iii++) fitTxt->AddText(a_snapshot.bottomTextLines[iii].c_str());
  fitTxt->SetFillStyle(0);
  fitTxt->SetBorderSize(5);

  TH1F* frame_pull = bottomPad->DrawFrame(hist->GetBinLowEdge(hist->FindBin(a_snapshot.lowRange - 0.5*a_snapshot.sigmaGuess)), -15.0,
         hist->GetBinLowEdge(hist->FindBin(a_snapshot.highRange + 0.5*a_snapshot.sigmaGuess)), 15.0);
  frame_pull->SetTitle(Form("; Z_{TPC}; Standardized Pull"));
  #ifndef _ZFITTER_PHYSMATH_BYPASS_
  pullHisto->Draw("Same");
  #endif
  TLine* oneLine = new TLine(a_snapshot.lowRange - 0.5*a_snapshot.sigmaGuess, 0.0, a_snapshot.highRange + 0.5*a_snapshot.sigmaGuess, 0.0);
  oneLine->SetLineStyle(7);
  oneLine->SetLineColor(kBlack);
  oneLine->Draw("same");
  fitTxt->Draw("same");

  gPad->Update();
  gSystem->ProcessEvents();
  fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/RapBin_%02d_STEP_%02d_mTm0Bin_%02d.png",
    m_imagePreDir.c_str(), m_imgDirName.c_str(), a_snapshot.particleName.c_str(),
    a_snapshot.centIndex, a_snapshot.rapBin, a_snapshot.fittingRound, a_snapshot.mtm0Bin));

  if(a_snapshot.saveNoLogImages){
    topPad->SetLogy(false);
    topPad->Update();
    gSystem->ProcessEvents();
    fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/NoLog_RapBin_%02d_STEP_%02d_mTm0Bin_%02d.png",
      m_imagePreDir.c_str(), m_imgDirName.c_str(), a_snapshot.particleName.c_str(),
      a_snapshot.centIndex, a_snapshot.rapBin, a_snapshot.fittingRound, a_snapshot.mtm0Bin));
  }

  #ifndef _SAVE_POINTERS_
    #ifndef _ZFITTER_PHYSMATH_BYPASS_
      delete combinedFunct;
      delete pionFunct;
      delete electronFunct;
      delete kaonFunct;
      if(pionFunct_init) delete pionFunct_init;
      if(koanFunct_init) delete koanFunct_init;
      if(electronFunct_init) delete electronFunct_init;
      delete pullHisto;
    #endif
    delete lowRangeLine;
    delete highRangeLine;
    delete oneLine;
    delete fitParsTxt;
    delete fitTxt;
    delete topPad;
    delete bottomPad;
    delete fittingCanvas;
  #endif
}

void ZFitter::renderProtonTPCDiagnosticPlot(const ZFitterDiagnosticSnapshot& a_snapshot){
  if(a_snapshot.hist == nullptr) return;

  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  fittingCanvas->Draw();
  TPad* topPad    = new TPad("topPad_P_TPC_snap","topPad_P_TPC_snap",0.0,0.3,1.0,1.0);
  TPad* bottomPad = new TPad("bottomPad_P_TPC_snap","bottomPad_P_TPC_snap",0.0,0.0,1.0,0.3);
  topPad->Draw();
  bottomPad->Draw();
  fittingCanvas->cd();
  topPad->cd();
  topPad->SetLogy(true);
  topPad->SetRightMargin(0.4);
  bottomPad->SetLogy(false);
  bottomPad->SetRightMargin(0.4);
  gStyle->SetOptFit(0111);

  TH1D* hist = a_snapshot.hist;
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerSize(0.5);
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlack);

  TString nameBase = Form("snap_%s_c%d_cent%d_r%03d_m%04d_s%d",
    a_snapshot.particleName.c_str(), a_snapshot.charge, a_snapshot.centIndex,
    a_snapshot.rapBin, a_snapshot.mtm0Bin, a_snapshot.fittingRound);

  #ifndef _ZFITTER_PHYSMATH_BYPASS_
    TF1* combinedFunct = PhysMath::getThreeSkewFunct(Form("%s_combined",nameBase.Data()),
      a_snapshot.N[0], a_snapshot.mu[0], a_snapshot.sigma[0], a_snapshot.gamma[0],
      a_snapshot.N[1], a_snapshot.mu[1], a_snapshot.sigma[1], a_snapshot.gamma[1],
      a_snapshot.N[2], a_snapshot.mu[2], a_snapshot.sigma[2], a_snapshot.gamma[2]);
    TF1* pionFunct    = PhysMath::skewNormalFunct(a_snapshot.N[0], a_snapshot.mu[0], a_snapshot.sigma[0], a_snapshot.gamma[0]);
    TF1* kaonFunct    = PhysMath::skewNormalFunct(a_snapshot.N[1], a_snapshot.mu[1], a_snapshot.sigma[1], a_snapshot.gamma[1]);
    TF1* protonFunct  = PhysMath::skewNormalFunct(a_snapshot.N[2], a_snapshot.mu[2], a_snapshot.sigma[2], a_snapshot.gamma[2]);
    pionFunct->SetLineColor(a_snapshot.color[0]);
    pionFunct->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
    kaonFunct->SetLineColor(a_snapshot.color[1]);
    kaonFunct->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
    protonFunct->SetLineColor(a_snapshot.color[2]);
    protonFunct->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
    TF1* pionFunct_init = NULL;
    TF1* kaonFunct_init = NULL;
    TF1* protonFunct_init = NULL;
    if(a_snapshot.drawInitialSeeds){
      protonFunct_init = PhysMath::skewNormalFunct(a_snapshot.seedN[2], a_snapshot.seedMu[2], a_snapshot.seedSigma[2], a_snapshot.seedGamma[2]);
      kaonFunct_init   = PhysMath::skewNormalFunct(a_snapshot.seedN[1], a_snapshot.seedMu[1], a_snapshot.seedSigma[1], a_snapshot.seedGamma[1]);
      pionFunct_init   = PhysMath::skewNormalFunct(a_snapshot.seedN[0], a_snapshot.seedMu[0], a_snapshot.seedSigma[0], a_snapshot.seedGamma[0]);
      pionFunct_init->SetLineColor(kYellow+2);
      pionFunct_init->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
      kaonFunct_init->SetLineColor(kYellow+1);
      kaonFunct_init->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
      protonFunct_init->SetLineColor(kOrange);
      protonFunct_init->SetRange(a_snapshot.lowRange, a_snapshot.highRange);
      if(!a_snapshot.doFit[0]) pionFunct_init->SetParameter(0,0.0);
      if(!a_snapshot.doFit[1]) kaonFunct_init->SetParameter(0,0.0);
      if(!a_snapshot.doFit[2]) protonFunct_init->SetParameter(0,0.0);
    }
  #endif

  double highestValue  = HistogramUtilities::getMaxInRange(hist, a_snapshot.currentLowFitRange, a_snapshot.currentHighFitRange);
  double lowestValue   = HistogramUtilities::getMinInRange(hist, a_snapshot.currentLowFitRange, a_snapshot.currentHighFitRange);
  TLine* lowRangeLine  = new TLine(a_snapshot.currentLowFitRange, lowestValue, a_snapshot.currentLowFitRange, highestValue);
  TLine* highRangeLine = new TLine(a_snapshot.currentHighFitRange, lowestValue, a_snapshot.currentHighFitRange, highestValue);
  lowRangeLine->SetLineColor(kViolet);
  highRangeLine->SetLineColor(kViolet);

  topPad->cd();
  HistogramUtilities::setAxisRanges(hist, a_snapshot.lowRange - 0.5*a_snapshot.sigmaGuess, a_snapshot.highRange + 0.5*a_snapshot.sigmaGuess, true);
  topPad->SetLogy(true);
  hist->SetStats(kFALSE);
  hist->Draw("E");
  #ifndef _ZFITTER_PHYSMATH_BYPASS_
    if(pionFunct_init) pionFunct_init->Draw("same");
    if(kaonFunct_init) kaonFunct_init->Draw("same");
    if(protonFunct_init) protonFunct_init->Draw("same");
    combinedFunct->Draw("same");
    pionFunct->Draw("same");
    protonFunct->Draw("same");
    kaonFunct->Draw("same");
  #endif
  lowRangeLine->Draw("same");
  highRangeLine->Draw("same");

  TPaveText* fitParsTxt = new TPaveText(0.65,0.05,0.95,0.95,"NDC");
  for(unsigned int iii = 0; iii < a_snapshot.topTextLines.size(); iii++) fitParsTxt->AddText(a_snapshot.topTextLines[iii].c_str());
  fitParsTxt->SetFillStyle(0);
  fitParsTxt->SetBorderSize(5);
  fitParsTxt->Draw("same");

  bottomPad->cd();
  TH1D* pullHisto = nullptr;
  #ifndef _ZFITTER_PHYSMATH_BYPASS_
  pullHisto = HistogramUtilities::makeFitPullHistogram(hist, combinedFunct, true);
  pullHisto->SetMarkerStyle(20);
  pullHisto->SetMarkerSize(0.5);
  pullHisto->SetMarkerColor(kBlack);
  pullHisto->SetLineColor(kBlack);
  pullHisto->GetXaxis()->SetTitle(Form("Z_{TPC}"));
  pullHisto->GetYaxis()->SetTitle(Form("Standardized Pull"));
  pullHisto->GetXaxis()->SetRangeUser(hist->GetBinLowEdge(hist->FindBin(a_snapshot.lowRange - 0.5*a_snapshot.sigmaGuess)), hist->GetBinLowEdge(hist->FindBin(a_snapshot.highRange + 0.5*a_snapshot.sigmaGuess)+1));
  pullHisto->GetYaxis()->SetRangeUser(-15,15);
  #endif

  TPaveText* fitTxt = new TPaveText(0.65,0.03,0.95,0.97,"NDC");
  for(unsigned int iii = 0; iii < a_snapshot.bottomTextLines.size(); iii++) fitTxt->AddText(a_snapshot.bottomTextLines[iii].c_str());
  fitTxt->SetFillStyle(0);
  fitTxt->SetBorderSize(5);

  if(pullHisto) pullHisto->Draw("E");
  TLine* oneLine = new TLine(a_snapshot.lowRange - 0.5*a_snapshot.sigmaGuess, 0.0, a_snapshot.highRange + 0.5*a_snapshot.sigmaGuess, 0.0);
  oneLine->SetLineStyle(7);
  oneLine->SetLineColor(kBlack);
  oneLine->Draw("same");
  fitTxt->Draw("same");
  gPad->Update();
  gSystem->ProcessEvents();
  fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/RapBin_%02d_STEP_%02d_mTm0Bin_%02d_%s.png",
    m_imagePreDir.c_str(), m_imgDirName.c_str(), a_snapshot.particleName.c_str(),
    a_snapshot.centIndex, a_snapshot.rapBin, a_snapshot.fittingRound, a_snapshot.mtm0Bin, a_snapshot.roundTag.c_str()));

  if(a_snapshot.saveNoLogImages){
    topPad->SetLogy(false);
    topPad->Update();
    gSystem->ProcessEvents();
    fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/NoLog_RapBin_%02d_STEP_%02d_mTm0Bin_%02d_%s.png",
      m_imagePreDir.c_str(), m_imgDirName.c_str(), a_snapshot.particleName.c_str(),
      a_snapshot.centIndex, a_snapshot.rapBin, a_snapshot.fittingRound, a_snapshot.mtm0Bin, a_snapshot.roundTag.c_str()));
  }

  #ifndef _SAVE_POINTERS_
    #ifndef _ZFITTER_PHYSMATH_BYPASS_
      delete combinedFunct;
      delete pionFunct;
      delete protonFunct;
      delete kaonFunct;
      if(pionFunct_init)   delete pionFunct_init;
      if(kaonFunct_init)   delete kaonFunct_init;
      if(protonFunct_init) delete protonFunct_init;
      delete pullHisto;
    #endif
    delete lowRangeLine;
    delete highRangeLine;
    delete oneLine;
    delete fitParsTxt;
    delete fitTxt;
    delete topPad;
    delete bottomPad;
    delete fittingCanvas;
  #endif
}
