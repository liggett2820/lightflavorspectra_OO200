// TestEnergyLossFix.C -- lightflavorspectra_OO200
//
// One-off sanity check for the fitEnergyLoss() fix in source/EfficiencyFitter.cxx
// (mTm0LossGraph now gets populated via HistogramUtilities::makeGraphAsymm() when
// m_useEnergyLossTProfiles is true -- see that commit's message for the full
// diagnosis). Restricts to ONE species and ONE rapidity bin via setRapBinsToTest()
// so this runs in seconds instead of the full multi-species job.
//
// Fill in embFile below with the same reprocessed-proton file path you already use
// in RunEfficiencyFitter_Embedding.C (RunEmbeddingBinner.C's output), then run from
// the repo root:
//   root -l -q -b macros/LoadEfficiencyFitterLibs.C macros/TestEnergyLossFix.C
//
// What to check afterward:
//   1. No more "ERROR: Something went really wrong... no energy loss graph" lines
//      in the console output for this bin.
//   2. A new "combined formula: ..." line DOES appear -- that only prints once the
//      residual/Legendre-fit code actually runs, which is the real proof the fix
//      is working (not just that the error message went away).
//   3. Open EfficiencyFitterImages_Test/proton/EnergyLoss_AllCent/Residules_yIndex_XX.png
//      and eyeball it -- real scattered points near zero with a sane fit line
//      through them, not empty or wildly diverging.

void TestEnergyLossFix(){

  EfficiencyFitter* fitter = new EfficiencyFitter();
  fitter->setImageDir("./EfficiencyFitterImages_Test/");
  fitter->setUseIntegralChiSqr(false);
  fitter->setNoFitNumBins(25, 120);

  fitter->setRapBinsToTest({0.0});   // mid-rapidity only -- seconds instead of minutes

  const int PROTON = 2;   // pion (0) is skipped entirely by fitEnergyLoss -- use kaon(1) or proton(2)

  string embFile = "FILL_IN_YOUR_REPROCESSED_PROTON_FILE_PATH_HERE.root"; // TODO: match RunEfficiencyFitter_Embedding.C's embFile_pro

  cout << "Loading embedding for proton..." << endl;
  fitter->loadEmbeddingFiles(embFile, "", PROTON, 1);   // charge +1, a_inFileNameTwo=""

  cout << "Fitting AllCent energy loss for proton..." << endl;
  fitter->fitEnergyLoss(PROTON, 1, -1);   // AllCent

  cout << "Done. Check console above for 'combined formula:' and look at the Residules PNG." << endl;
}
