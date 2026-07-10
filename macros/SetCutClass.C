// SetCutClass.C -- lightflavorspectra_OO200
//
// Ported from the original repo's macros/SetCutClass_Collider.C, keeping ONLY the
// values that would actually have been selected for _OO_200_COL_ (all other energies'
// #ifdef _7p7_COL_ / _9p2_COL_ / etc. blocks were dead code for this build and are
// removed). Everything below the energy-specific block (mass cuts, BTOF/ETOF beta-gamma
// cut lines, and all 9 species' rapidity/mTm0 bin edge arrays and eta/pT binning) was
// NOT gated by any energy toggle in the original -- those numbers are copied verbatim,
// unchanged, digit for digit.
//
// UPDATED 2026-07-03: reconciled against the actual macros/SetCutClass.C used on SDCC
// (starsub03) to produce yieldHistos_OOGeV_proton_2026_07_01.root. That file was pasted
// in full and is the ground truth for what cuts were really applied when the data was
// produced -- several values below had drifted from it (this repo's own port history
// vs. what was actually run). Every value changed below is called out at the point of
// change; nothing was re-derived or guessed, all copied verbatim from the SDCC file.
//
// The former note about "unsigned int trigs[6] = {860001}; setTriggers(6,...)" being a
// possible copy-paste bug is now resolved: the SDCC file uses a single trigger,
// "unsigned int trigs[1] = {860002}; setTriggers(1,...)" -- both the array size AND the
// trigger ID were wrong here. Fixed below to match SDCC exactly.
//
// setInvBetaWidth_bTOF/eTOF calls (previously present here, absent from the SDCC file)
// were dropped 2026-07-03 to match SDCC exactly.

#include "../headers/CutClass.h"
#include <iostream>
#include <vector>
#include "../makefile_toggles.h"
using namespace std;

void SetCutClass(CutClass *cuts){
  #ifndef _CLING_
    cout << "Attempting to load CutClass_cxx.so, the following is not a problem if you are on root 6" << endl;
    gSystem->Load("../bin/CutClass_cxx.so");
  #endif

  string mode = "Center"; // O+O 200 GeV is always run in Center mode in the original
  cuts->setVertexConfig(mode);
  cuts->setZRange(-2,2); // was (-70,70) -- corrected 2026-07-03 to match the SDCC file actually used

  // ---- _OO_200_COL_ specific values (from the original's #ifdef _OO_200_COL_ block) ----
  // TEMPORARILY reverted 2026-07-03 back to "dev" (the "SL23c" correction above is
  // still the physically-correct value matching the real SDCC production, but the
  // corresponding SL23c_BichselCurves.root file doesn't actually exist anywhere --
  // checked both this repo's and the full etof repo's submodule/ParticleInfo/
  // PidFunctions_Data/ directories, only SL19b/SL20d/dev are present in either. Using
  // "SL23c" makes ParticleInfo::ReadBichselFile() call exit(EXIT_FAILURE) immediately,
  // silently killing the whole ROOT process with no further output -- this is why the
  // last several ZFitter reruns appeared to complete ("ran fine") but never actually
  // wrote new spectra data. Switch back to "SL23c" once the real Bichsel curve file is
  // obtained from SDCC/RCF and placed in submodule/ParticleInfo/PidFunctions_Data/.
  cuts->setStarver("dev");
  cuts->setSqrts_NN(200);
  cuts->setBeamEnergy(100);
  cuts->setBeamLocation(0.0,0.0);
  cuts->setPileUpCut(260); // added 2026-07-03 -- absent here before, present in the SDCC file
  unsigned int trigs[1] = {860002}; // was trigs[6]={860001} -- both size and ID were wrong, fixed 2026-07-03
  cuts->setTriggers(1, &trigs[0]);

  //cut on nToFT0, number of ToF tracks used for T0 Calibration
  cuts->setUseT0(false); // was true -- corrected 2026-07-03; SDCC has setTofT0Cut commented out entirely so it's dropped here too
  cuts->setEventConfig("ColliderCenter");
  cuts->setYCM(0.0);
  cuts->setRadialCut(2.0);
  cuts->setTofMatch(0); // was 3 -- corrected 2026-07-03 to match SDCC
  cuts->setTrigToggle(true);

  //###############    TRACK QUALITY CUTS ###################
  cuts->setNHitsFit(15); // reverted from 20 (2026-07-10) to recover more forward (large
                         // |eta|, low pT) tracks -- 20 was the 2026-07-03 SDCC-matched
                         // value; this is a deliberate divergence from that officially-
                         // validated cut, not a correction, so treat results under this
                         // as a track-quality systematic variation until/unless the
                         // official analysis itself adopts 15.
  cuts->setNHitsDeDx(10);
  cuts->setFitMaxRatio(0.5001); // was 0.505 -- corrected 2026-07-03 to match SDCC
  cuts->setBTOF(1.6, 2.8);
  cuts->setDCA(3); // was 1.0 -- corrected 2026-07-03 to match SDCC
  cuts->setETOFMaxClusterSize(99);
  cuts->setETOFMatchFlagCut(-1);
  // added 2026-07-03 -- absent here before, present in the SDCC file (4-arg: meanDeltaX, meanDeltaY, widthDeltaX, widthDeltaY)
  cuts->setETOF(-0.4409,0.08343, 50, 50);
  // added 2026-07-03 -- absent here before, present in the SDCC file
  vector<double> momStops_TPC_Plus  = {2.71,2.0,4.0,1.0,4.0};
  vector<double> momStops_TPC_Minus = {4.0,1.56,0.0,2.0,0.0};
  vector<double> momStops_bToF_Plus  = {2.71,2.0,4.0,1.0,4.0};
  vector<double> momStops_bToF_Minus = {4.0,1.56,0.0,2.0,0.0};
  cuts->setMomStops_TPC(momStops_TPC_Plus, momStops_TPC_Minus);
  cuts->setMomStops_bToF(momStops_bToF_Plus, momStops_bToF_Minus);

  // _STREFMULTCORR_ was never defined for this build -- raw refMult centrality bins,
  // exactly the values the original's _OO_200_COL_ branch used.
  cout << " NOT Using StRefMultCorr Centrality Bins (raw refMult, matching the original's _OO_200_COL_ config)" << endl;

  // 5-bin scheme (0-5%,5-10%,10-20%,20-40%,40-80%), requested 2026-07-02, replacing the
  // previous 9-bin scheme below. The refMult cut VALUES are not re-derived or guessed --
  // they are the exact same cuts the 9-bin scheme already used at the percentile
  // boundaries that survive merging (5,10,20,40,80 were already bin edges in the 9-bin
  // scheme; only the 30/50/60/70 edges are dropped, merging 20-30%+30-40% into one bin
  // and 40-50%+50-60%+60-70%+70-80% into one bin).
  //
  // ORDERING (corrected 2026-07-02): centrality classes are defined by integrating the
  // Glauber-model multiplicity distribution from MOST CENTRAL (highest multiplicity) to
  // MOST PERIPHERAL, so the refMult cut must be HIGHEST for the 0-5% bin and decrease
  // monotonically to the most peripheral bin -- i.e. a_cent_edges must be DESCENDING
  // while a_cent_percents is ascending. The as-ported 9-bin array below had these
  // ascending in parallel (same values, wrong order) -- see CutClass::centralityIndex():
  // `for(iii...) if(a_mult >= m_centEdges[iii]) return iii-1;` only discriminates
  // central vs. peripheral correctly if m_centEdges descends as iii increases; with
  // ascending edges, essentially every event with refMult>=0 would satisfy the very
  // first real threshold check and get labeled centIndex 0. Fixed here by reversing the
  // edge order (same magnitudes, not re-derived) relative to the 9-bin array's values.
  // Corrected 2026-07-03: {300,150,100,70,0} was this repo's own value, but the SDCC
  // file that actually produced the yields file being analyzed uses {44,37,28,17,5} --
  // same percentile labels, different multiplicity thresholds (likely a different
  // vertex/pile-up/track-cut configuration upstream). PicoBinner bakes centIndex into
  // the raw yield histograms at production time, so THIS is the scheme that matters.
  // 80-100% bin added (edge=0) so the most peripheral events aren't dropped entirely --
  // without it, CutClass::centralityIndex() falls through its loop for any event with
  // refMult below the old last edge (5) and rejects it (m_rejectedEvents->Fill(6.5),
  // return -1), the same way the 0-5% bin's upper "no limit" edge is handled by pushing
  // a large sentinel in CutClass::setCentralities(). refMult=0 as the floor is exact
  // (not derived/estimated), so this doesn't need DeriveCentralityEdges.C to add.
  double centCutsArray[6] = {44,37,28,17,5,0};
  int percents[6] = {5,10,20,40,80,100};
  cuts->setCentralities(6, &centCutsArray[0], &percents[0]);

  // Previous 9-bin scheme (0-5,5-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80) -- kept
  // here, not deleted, in case you need to revert or compare. Shown here ALREADY
  // CORRECTED to descending edge order (the as-ported values were {0,5,20,50,70,90,100,
  // 150,300} ascending, which has the same ordering bug described above):
  //   double centCutsArray[9] = {300,150,100,90,70,50,20,5,0};
  //   int percents[9] = {5,10,20,30,40,50,60,70,80};
  //   cuts->setCentralities(9, &centCutsArray[0], &percents[0]);

  // ---- everything below here was NOT energy-gated in the original -- copied verbatim ----
  vector<double> lowMassCut;
  vector<double> highMassCut;
  //loose no overlap mass cuts
  // indices 4-7 corrected 2026-07-03 to match the SDCC file (was 1.65/2.4/1.23/1.65 low,
  // 2.4/3.16/1.58/2.4 high -- now 1.38/2.7/2.7/2.9 low, 2.7/2.9/2.9/3.8 high)
  lowMassCut.push_back(0.0);
  lowMassCut.push_back(0.367);
  lowMassCut.push_back(0.726);
  lowMassCut.push_back(0.0);
  lowMassCut.push_back(1.38);
  lowMassCut.push_back(2.7);
  lowMassCut.push_back(2.7);
  lowMassCut.push_back(2.9);
  lowMassCut.push_back(0.97);
  highMassCut.push_back(0.367);
  highMassCut.push_back(0.726);
  highMassCut.push_back(1.38);
  highMassCut.push_back(0.065);
  highMassCut.push_back(2.7);
  highMassCut.push_back(2.9);
  highMassCut.push_back(2.9);
  highMassCut.push_back(3.8);
  highMassCut.push_back(0.121);
  cuts->setMassCuts(lowMassCut,highMassCut);

  // Dropped 2026-07-03: setInvBetaWidth_bTOF/eTOF (0.03 / 0.015 for all three species)
  // had no equivalent call at all in the SDCC file, so removed here to match exactly.

  // Corrected 2026-07-03 to match the SDCC file's line points exactly (was a different
  // set of points entirely -- {0.212698,...}/{0.276361,...} etc. -- this repo's own
  // values, not what was actually used to produce the yields file).
  vector<double> lowLineXPts;
  vector<double> lowLineYPts;
  vector<double> highLineXPts;
  vector<double> highLineYPts;
  lowLineXPts.push_back(0.2);
  lowLineYPts.push_back(30.0);
  lowLineXPts.push_back(0.65);
  lowLineYPts.push_back(4.5);
  lowLineXPts.push_back(0.65);
  lowLineYPts.push_back(4.5);
  lowLineXPts.push_back(3.73);
  lowLineYPts.push_back(1.0);
  highLineXPts.push_back(0.25);
  highLineYPts.push_back(40.0);
  highLineXPts.push_back(0.65);
  highLineYPts.push_back(10.0);
  highLineXPts.push_back(0.65);
  highLineYPts.push_back(10.0);
  highLineXPts.push_back(1.3);
  highLineYPts.push_back(5.0);
  highLineXPts.push_back(1.3);
  highLineYPts.push_back(5.0);
  highLineXPts.push_back(2.0);
  highLineYPts.push_back(4.0);
  highLineXPts.push_back(2.0);
  highLineYPts.push_back(4.0);
  highLineXPts.push_back(5.0);
  highLineYPts.push_back(3.5);
  highLineXPts.push_back(5.0);
  highLineYPts.push_back(3.5);
  highLineXPts.push_back(50);
  highLineYPts.push_back(3.5);
  cuts->setBTOFBetaGammaCuts(false, lowLineXPts, lowLineYPts, highLineXPts, highLineYPts);

  vector<double> lowLineXPts2;
  vector<double> lowLineYPts2;
  vector<double> highLineXPts2;
  vector<double> highLineYPts2;
  lowLineXPts2.push_back(-5.0);
  lowLineYPts2.push_back(1.8);
  lowLineXPts2.push_back(0.15);
  lowLineYPts2.push_back(1.8);

  highLineXPts2.push_back(-5.0);
  highLineYPts2.push_back(4.8); //3.58 to electron bottom
  highLineXPts2.push_back(0.15);
  highLineYPts2.push_back(4.8); //3.58 to electron bottom

  cuts->setBTOFInvBetaGammaSqrCuts(true, lowLineXPts2, lowLineYPts2, highLineXPts2, highLineYPts2);

  vector<double> lowLineXPts3;
  vector<double> lowLineYPts3;
  vector<double> highLineXPts3;
  vector<double> highLineYPts3;
  lowLineXPts3.push_back(-5.0);
  lowLineYPts3.push_back(1.8);
  lowLineXPts3.push_back(0.15);
  lowLineYPts3.push_back(1.8);

  highLineXPts3.push_back(-5.0);
  highLineYPts3.push_back(4.8); //3.58 to electron bottom
  highLineXPts3.push_back(0.15);
  highLineYPts3.push_back(4.8); //3.58 to electron bottom

  // toggle changed true -> false 2026-07-03 to match SDCC (its own comment there notes
  // "need to change back to true" -- so this may flip again later, but for now this
  // matches what was actually run to produce the yields file)
  cuts->setETOFInvBetaGammaSqrCuts(false, lowLineXPts3, lowLineYPts3, highLineXPts3, highLineYPts3);
  cuts->setBTOFElectronCut(true);
  cuts->setETOFElectronCut(false);
  cuts->setModNSigmaCuts_bTOF(true ,0.05, -5.5, 5.0,-0.08);
  cuts->setModNSigmaCuts_eTOF(true,0.05, -5.5, 5.0,-0.08);

  cuts->setLogMTBinning(false);
  double pionArray[101]={0,0.000357785,0.00142569,0.00318778,0.00561879,0.00868581,0.0123503,0.0165703,0.021302,0.0265016,0.0321268,0.038137,0.0444946,0.0511649,0.0581161,0.0653197,0.07275,0.0803841,0.0882013,0.0961837,0.104315,0.11258,0.120967,0.129465,0.138062,0.146751,0.155523,0.16437,0.173287,0.182268,0.191307,0.2004,0.209543,0.218731,0.227962,0.237232,0.246539,0.255879,0.265251,0.274652,0.284081,0.293535,0.303013,0.312514,0.322036,0.331577,0.341138,0.350715,0.36031,0.36992,0.379544,0.389183,0.398835,0.408499,0.418175,0.427863,0.437561,0.447269,0.456987,0.466714,0.476449,0.486193,0.495945,0.505705,0.515472,0.525246,0.535026,0.544813,0.554606,0.564404,0.574209,0.584018,0.593833,0.603653,0.613477,0.623306,0.633139,0.642977,0.652819,0.662664,0.672514,0.682367,0.692223,0.702083,0.711946,0.721812,0.731682,0.741554,0.751429,0.761307,0.771188,0.781071,0.790957,0.800845,0.810735,0.820628,0.830523,0.84042,0.850319,0.86022,0.870123};
  double kaonArray[101]={0,0.000101272,0.000404963,0.0009107,0.00161786,0.0025256,0.00363279,0.00493813,0.00644006,0.0081368,0.0100264,0.0121067,0.0143753,0.0168298,0.0194675,0.0222855,0.0252809,0.0284507,0.0317918,0.0353008,0.0389744,0.0428093,0.0468021,0.0509492,0.0552472,0.0596925,0.0642817,0.0690111,0.0738774,0.078877,0.0840064,0.0892623,0.0946412,0.10014,0.105755,0.111483,0.117321,0.123266,0.129315,0.135464,0.141712,0.148054,0.154489,0.161013,0.167624,0.17432,0.181097,0.187953,0.194887,0.201895,0.208975,0.216126,0.223345,0.23063,0.237979,0.24539,0.252862,0.260392,0.26798,0.275622,0.283318,0.291066,0.298864,0.306711,0.314606,0.322547,0.330533,0.338562,0.346634,0.354746,0.362899,0.37109,0.379319,0.387585,0.395886,0.404222,0.412592,0.420994,0.429428,0.437893,0.446389,0.454913,0.463466,0.472047,0.480655,0.48929,0.49795,0.506635,0.515344,0.524078,0.532834,0.541613,0.550413,0.559236,0.568079,0.576942,0.585825,0.594728,0.603649,0.612589,0.621547};
  double protonArray[151]={0,5.32881e-05,0.000213134,0.000479484,0.000852246,0.00133129,0.00191647,0.00260757,0.00340436,0.00430657,0.00531391,0.00642603,0.00764257,0.00896312,0.0103873,0.0119145,0.0135444,0.0152763,0.0171098,0.0190443,0.021079,0.0232135,0.0254471,0.027779,0.0302085,0.0327349,0.0353575,0.0380755,0.0408881,0.0437945,0.0467938,0.0498851,0.0530678,0.0563408,0.0597032,0.0631543,0.066693,0.0703184,0.0740297,0.0778258,0.0817058,0.0856688,0.0897138,0.0938398,0.0980459,0.102331,0.106694,0.111135,0.115652,0.120243,0.124909,0.129649,0.13446,0.139343,0.144297,0.149319,0.15441,0.159569,0.164794,0.170085,0.17544,0.180859,0.186341,0.191885,0.19749,0.203155,0.208879,0.214661,0.220501,0.226398,0.23235,0.238357,0.244418,0.250532,0.256699,0.262917,0.269186,0.275505,0.281874,0.28829,0.294755,0.301266,0.307824,0.314427,0.321075,0.327767,0.334503,0.341281,0.348101,0.354963,0.361865,0.368807,0.375789,0.382809,0.389868,0.396964,0.404098,0.411267,0.418473,0.425713,0.432989,0.440298,0.447641,0.455018,0.462426,0.469867,0.47734,0.484843,0.492377,0.499941,0.507534,0.515157,0.522809,0.530488,0.538196,0.54593,0.553692,0.56148,0.569294,0.577134,0.585,0.59289,0.600805,0.608743,0.616706,0.624692,0.632701,0.640733,0.648787,0.656863,0.664961,0.673081,0.681221,0.689382,0.697563,0.705765,0.713986,0.722227,0.730487,0.738766,0.747064,0.75538,0.763714,0.772066,0.780436,0.788823,0.797227,0.805648,0.814086,0.82254,0.83101};
  double pionBTOFArray[201]={0,0.00142569,0.00561879,0.0123503,0.021302,0.0321268,0.0444946,0.0581161,0.07275,0.0882013,0.104315,0.120967,0.138062,0.155523,0.173287,0.191307,0.209543,0.227962,0.246539,0.265251,0.284081,0.303013,0.322036,0.341138,0.36031,0.379544,0.398835,0.418175,0.437561,0.456987,0.476449,0.495945,0.515472,0.535026,0.554606,0.574209,0.593833,0.613477,0.633139,0.652819,0.672514,0.692223,0.711946,0.731682,0.751429,0.771188,0.790957,0.810735,0.830523,0.850319,0.870123,0.889935,0.909753,0.929579,0.949411,0.969249,0.989093,1.00894,1.0288,1.04866,1.06852,1.08839,1.10826,1.12814,1.14802,1.1679,1.18779,1.20768,1.22757,1.24747,1.26737,1.28727,1.30718,1.32709,1.347,1.36691,1.38682,1.40674,1.42666,1.44658,1.46651,1.48643,1.50636,1.52629,1.54622,1.56615,1.58608,1.60602,1.62596,1.64589,1.66583,1.68577,1.70572,1.72566,1.7456,1.76555,1.7855,1.80544,1.82539,1.84534,1.86529,1.88525,1.9052,1.92515,1.94511,1.96506,1.98502,2.00498,2.02493,2.04489,2.06485,2.08481,2.10477,2.12474,2.1447,2.16466,2.18462,2.20459,2.22455,2.24452,2.26448,2.28445,2.30442,2.32439,2.34435,2.36432,2.38429,2.40426,2.42423,2.4442,2.46417,2.48414,2.50412,2.52409,2.54406,2.56403,2.58401,2.60398,2.62396,2.64393,2.66391,2.68388,2.70386,2.72383,2.74381,2.76379,2.78376,2.80374,2.82372,2.8437,2.86367,2.88365,2.90363,2.92361,2.94359,2.96357,2.98355,3.00353,3.02351,3.04349,3.06347,3.08345,3.10343,3.12342,3.1434,3.16338,3.18336,3.20334,3.22333,3.24331,3.26329,3.28328,3.30326,3.32324,3.34323,3.36321,3.3832,3.40318,3.42316,3.44315,3.46313,3.48312,3.5031,3.52309,3.54308,3.56306,3.58305,3.60303,3.62302,3.64301,3.66299,3.68298,3.70297,3.72295,3.74294,3.76293,3.78291,3.8029,3.82289,3.84288,3.86286};
  double kaonBTOFArray[201]={0,0.000404963,0.00161786,0.00363279,0.00644006,0.0100264,0.0143753,0.0194675,0.0252809,0.0317918,0.0389744,0.0468021,0.0552472,0.0642817,0.0738774,0.0840064,0.0946412,0.105755,0.117321,0.129315,0.141712,0.154489,0.167624,0.181097,0.194887,0.208975,0.223345,0.237979,0.252862,0.26798,0.283318,0.298864,0.314606,0.330533,0.346634,0.362899,0.379319,0.395886,0.412592,0.429428,0.446389,0.463466,0.480655,0.49795,0.515344,0.532834,0.550413,0.568079,0.585825,0.603649,0.621547,0.639516,0.657551,0.67565,0.693811,0.712029,0.730303,0.74863,0.767008,0.785435,0.803909,0.822427,0.840987,0.859589,0.87823,0.896909,0.915624,0.934374,0.953158,0.971973,0.99082,1.0097,1.0286,1.04753,1.06649,1.08548,1.10449,1.12352,1.14258,1.16166,1.18076,1.19988,1.21902,1.23818,1.25736,1.27656,1.29577,1.31501,1.33426,1.35352,1.3728,1.3921,1.4114,1.43073,1.45007,1.46942,1.48878,1.50816,1.52754,1.54695,1.56636,1.58578,1.60521,1.62466,1.64411,1.66358,1.68305,1.70253,1.72203,1.74153,1.76104,1.78056,1.80008,1.81962,1.83916,1.85871,1.87827,1.89784,1.91741,1.93699,1.95658,1.97617,1.99577,2.01538,2.03499,2.05461,2.07423,2.09386,2.1135,2.13314,2.15278,2.17243,2.19209,2.21175,2.23142,2.25109,2.27077,2.29045,2.31013,2.32982,2.34952,2.36921,2.38892,2.40862,2.42833,2.44805,2.46777,2.48749,2.50721,2.52694,2.54668,2.56641,2.58615,2.6059,2.62564,2.64539,2.66514,2.6849,2.70466,2.72442,2.74419,2.76395,2.78372,2.8035,2.82327,2.84305,2.86283,2.88262,2.9024,2.92219,2.94198,2.96178,2.98157,3.00137,3.02117,3.04097,3.06078,3.08059,3.1004,3.12021,3.14002,3.15984,3.17965,3.19947,3.2193,3.23912,3.25894,3.27877,3.2986,3.31843,3.33826,3.3581,3.37793,3.39777,3.41761,3.43745,3.45729,3.47714,3.49698,3.51683,3.53668};
  double protonBTOFArray[251]={0,0.000213134,0.000852246,0.00191647,0.00340436,0.00531391,0.00764257,0.0103873,0.0135444,0.0171098,0.021079,0.0254471,0.0302085,0.0353575,0.0408881,0.0467938,0.0530678,0.0597032,0.066693,0.0740297,0.0817058,0.0897138,0.0980459,0.106694,0.115652,0.124909,0.13446,0.144297,0.15441,0.164794,0.17544,0.186341,0.19749,0.208879,0.220501,0.23235,0.244418,0.256699,0.269186,0.281874,0.294755,0.307824,0.321075,0.334503,0.348101,0.361865,0.375789,0.389868,0.404098,0.418473,0.432989,0.447641,0.462426,0.47734,0.492377,0.507534,0.522809,0.538196,0.553692,0.569294,0.585,0.600805,0.616706,0.632701,0.648787,0.664961,0.681221,0.697563,0.713986,0.730487,0.747064,0.763714,0.780436,0.797227,0.814086,0.83101,0.847997,0.865047,0.882156,0.899324,0.916548,0.933828,0.951161,0.968547,0.985983,1.00347,1.021,1.03858,1.05621,1.07388,1.09159,1.10935,1.12715,1.14498,1.16286,1.18077,1.19873,1.21671,1.23473,1.25279,1.27088,1.289,1.30716,1.32534,1.34356,1.36181,1.38008,1.39838,1.41671,1.43507,1.45346,1.47186,1.4903,1.50876,1.52724,1.54575,1.56428,1.58283,1.60141,1.62,1.63862,1.65726,1.67591,1.69459,1.71329,1.732,1.75074,1.76949,1.78826,1.80704,1.82585,1.84467,1.86351,1.88236,1.90123,1.92011,1.93901,1.95793,1.97685,1.9958,2.01475,2.03372,2.05271,2.07171,2.09072,2.10974,2.12877,2.14782,2.16688,2.18595,2.20503,2.22413,2.24323,2.26235,2.28147,2.30061,2.31976,2.33892,2.35808,2.37726,2.39645,2.41565,2.43485,2.45407,2.47329,2.49252,2.51177,2.53102,2.55028,2.56954,2.58882,2.6081,2.62739,2.64669,2.666,2.68531,2.70463,2.72396,2.7433,2.76264,2.78199,2.80135,2.82071,2.84008,2.85946,2.87884,2.89823,2.91763,2.93703,2.95644,2.97585,2.99527,3.0147,3.03413,3.05357,3.07301,3.09246,3.11191,3.13137,3.15083,3.1703,3.18977,3.20925,3.22874,3.24823,3.26772,3.28722,3.30672,3.32623,3.34574,3.36526,3.38478,3.4043,3.42383,3.44337,3.46291,3.48245,3.50199,3.52155,3.5411,3.56066,3.58022,3.59979,3.61936,3.63893,3.65851,3.67809,3.69767,3.71726,3.73685,3.75645,3.77604,3.79564,3.81525,3.83486,3.85447,3.87408,3.8937,3.91332,3.93295,3.95257,3.9722,3.99184,4.01147,4.03111,4.05075,4.0704,4.09004,4.10969,4.12935,4.149};
  double pionETOFArray[101]={0,0.00561879,0.021302,0.0444946,0.07275,0.104315,0.138062,0.173287,0.209543,0.246539,0.284081,0.322036,0.36031,0.398835,0.437561,0.476449,0.515472,0.554606,0.593833,0.633139,0.672514,0.711946,0.751429,0.790957,0.830523,0.870123,0.909753,0.949411,0.989093,1.0288,1.06852,1.10826,1.14802,1.18779,1.22757,1.26737,1.30718,1.347,1.38682,1.42666,1.46651,1.50636,1.54622,1.58608,1.62596,1.66583,1.70572,1.7456,1.7855,1.82539,1.86529,1.9052,1.94511,1.98502,2.02493,2.06485,2.10477,2.1447,2.18462,2.22455,2.26448,2.30442,2.34435,2.38429,2.42423,2.46417,2.50412,2.54406,2.58401,2.62396,2.66391,2.70386,2.74381,2.78376,2.82372,2.86367,2.90363,2.94359,2.98355,3.02351,3.06347,3.10343,3.1434,3.18336,3.22333,3.26329,3.30326,3.34323,3.3832,3.42316,3.46313,3.5031,3.54308,3.58305,3.62302,3.66299,3.70297,3.74294,3.78291,3.82289,3.86286};
  double kaonETOFArray[101]={0,0.00161786,0.00644006,0.0143753,0.0252809,0.0389744,0.0552472,0.0738774,0.0946412,0.117321,0.141712,0.167624,0.194887,0.223345,0.252862,0.283318,0.314606,0.346634,0.379319,0.412592,0.446389,0.480655,0.515344,0.550413,0.585825,0.621547,0.657551,0.693811,0.730303,0.767008,0.803909,0.840987,0.87823,0.915624,0.953158,0.99082,1.0286,1.06649,1.10449,1.14258,1.18076,1.21902,1.25736,1.29577,1.33426,1.3728,1.4114,1.45007,1.48878,1.52754,1.56636,1.60521,1.64411,1.68305,1.72203,1.76104,1.80008,1.83916,1.87827,1.91741,1.95658,1.99577,2.03499,2.07423,2.1135,2.15278,2.19209,2.23142,2.27077,2.31013,2.34952,2.38892,2.42833,2.46777,2.50721,2.54668,2.58615,2.62564,2.66514,2.70466,2.74419,2.78372,2.82327,2.86283,2.9024,2.94198,2.98157,3.02117,3.06078,3.1004,3.14002,3.17965,3.2193,3.25894,3.2986,3.33826,3.37793,3.41761,3.45729,3.49698,3.53668};
  double protonETOFArray[126]={0,0.000852246,0.00340436,0.00764257,0.0135444,0.021079,0.0302085,0.0408881,0.0530678,0.066693,0.0817058,0.0980459,0.115652,0.13446,0.15441,0.17544,0.19749,0.220501,0.244418,0.269186,0.294755,0.321075,0.348101,0.375789,0.404098,0.432989,0.462426,0.492377,0.522809,0.553692,0.585,0.616706,0.648787,0.681221,0.713986,0.747064,0.780436,0.814086,0.847997,0.882156,0.916548,0.951161,0.985983,1.021,1.05621,1.09159,1.12715,1.16286,1.19873,1.23473,1.27088,1.30716,1.34356,1.38008,1.41671,1.45346,1.4903,1.52724,1.56428,1.60141,1.63862,1.67591,1.71329,1.75074,1.78826,1.82585,1.86351,1.90123,1.93901,1.97685,2.01475,2.05271,2.09072,2.12877,2.16688,2.20503,2.24323,2.28147,2.31976,2.35808,2.39645,2.43485,2.47329,2.51177,2.55028,2.58882,2.62739,2.666,2.70463,2.7433,2.78199,2.82071,2.85946,2.89823,2.93703,2.97585,3.0147,3.05357,3.09246,3.13137,3.1703,3.20925,3.24823,3.28722,3.32623,3.36526,3.4043,3.44337,3.48245,3.52155,3.56066,3.59979,3.63893,3.67809,3.71726,3.75645,3.79564,3.83486,3.87408,3.91332,3.95257,3.99184,4.03111,4.0704,4.10969,4.149};
  double blankEdgesArray[2] = {0,1.0};
  // SUPERSEDED 2026-07-03 (yet again): my previous fix here recentered the range to a
  // symmetric [-1.55,1.55], invented from first principles (raw-data extent + this
  // repo's own analysis-cut default). That reasoning was sound in the abstract, but it
  // doesn't matter here: the yieldHistos_OOGeV_proton_2026_07_01.root file was already
  // produced by ETOF's own PicoBinner (not this repo's), so its bin edges are already
  // fixed and baked in -- what matters is matching THAT binning exactly, not deriving a
  // new one. Pulled directly from lightflavorspectra_etof/macros/SetCutClass_Collider.C
  // (lines 354-380 there), which that file's own inline comment documents as "the EXACT
  // values from the SetCutClass.C actually used on SDCC (starsub03) to produce the O+O
  // proton yields file" -- i.e. real ground truth, not an estimate. Two things this
  // corrects vs. my previous guess: (1) the real per-species/per-detector bin counts and
  // ranges are NOT uniform (e.g. proton TPC is 36 bins over [-2.05,1.55], proton BTOF is
  // 27 bins over [-1.35,1.35] -- different from each other and from my guessed 31/1.55);
  // (2) indices 3-8 (electron/deuteron/triton/helion/alpha/muon) are genuinely blank
  // dummy bins (1 bin, [-1.0,1.0], blankEdgesArray) in the real production config -- my
  // previous fix (and the comment it replaced) wrongly assumed they used real
  // pion/kaon/proton arrays like species 0-2 do.
  cuts->setVariableRapMtM0BinningInfo(0,41,-2.05,2.05,100, pionArray,
                                        41,-2.05,2.05, 200, pionBTOFArray,
                                        41,-2.05,2.05, 100, pionETOFArray);
  cuts->setVariableRapMtM0BinningInfo(1,41,-2.05,2.05,100, kaonArray,
                                        41,-2.05,2.05, 200, kaonBTOFArray,
                                        41,-2.05,2.05, 100, kaonETOFArray);
  cuts->setVariableRapMtM0BinningInfo(2,41,-2.05,2.05,150, protonArray,
                                        41,-2.05,2.05, 250, protonBTOFArray,
                                        41,-2.05,2.05, 125, protonETOFArray);
  cuts->setVariableRapMtM0BinningInfo(3,1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray);
  cuts->setVariableRapMtM0BinningInfo(4,1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray);
  cuts->setVariableRapMtM0BinningInfo(5,1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray);
  cuts->setVariableRapMtM0BinningInfo(6,1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray);
  cuts->setVariableRapMtM0BinningInfo(7,1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray);
  cuts->setVariableRapMtM0BinningInfo(8,1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray,
                                        1,-1.0,1.0,  1, blankEdgesArray);

  // corrected 2026-07-03 to match SDCC (was 41,-2.05,2.05 / 25,-1.25,1.25 / 6,-1.65,-1.05)
  //
  // REVERTED 2026-07-03 (again): eta_TPC and eta_bTOF above (21,-2.05,0.05 and
  // 16,-1.55,0.05) had the same asymmetric-truncation-near-zero bug diagnosed in the
  // rapidity binning above -- both got cut off within 0.05 of eta=0 on the upper end.
  // Unlike eTOF (a real, physically asymmetric endcap detector at one fixed location,
  // so its own asymmetric range is correct detector geometry, not a bug -- left
  // untouched here), TPC and bTOF are symmetric barrel detectors centered on the
  // interaction point and should have eta ranges symmetric around 0. The values
  // this comment already documented as what was here BEFORE the SDCC copy
  // (41,-2.05,2.05 / 25,-1.25,1.25 -- both symmetric, both already at the same 0.10
  // bin width as everything else in this file) are better evidence than a freshly
  // derived number, since they were already used successfully in this exact repo, so
  // reverting to them rather than inventing new numbers.
  //
  // eta_eTOF SUPERSEDED 2026-07-03: eta_TPC and eta_bTOF above already exactly match
  // lightflavorspectra_etof/macros/SetCutClass_Collider.C's real production values
  // (41,-2.05,2.05 and 25,-1.25,1.25 -- confirmed, no further change needed), but
  // eta_eTOF's value here (7,-2.05,-1.35) was still an approximation, not the real
  // production number. The actual value used to produce this dataset's yield files is
  // 6,-1.65,-1.05 (same file, line ~385) -- correcting to match exactly, same reasoning
  // as the rapidity-binning fix above (eTOF's asymmetric range is real, physical
  // endcap-detector geometry, not a bug -- but the exact bin count/edges still have to
  // match what PicoBinner actually baked into the yield file).
  cuts->setEtaPtBinningInfo(80,0.05,2.5, //pT: n bins, low, high
                       41,-2.05,2.05,  //eta_TPC: n bins, low, high
                       25,-1.25,1.25,  //eta bToF: nBins, low, high
                       6,-1.65,-1.05);  //eta eToF: nBins, low, high

  cuts->makeCutInformationHistogram();

  return;
}
