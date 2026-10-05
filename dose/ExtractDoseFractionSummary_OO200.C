// ExtractDoseFractionSummary_OO200.C -- lightflavorspectra_OO200
//
// PURPOSE: turns Andrew's three-source dose-fraction framing into actual numbers per
// centrality bin:
//   fireball-produced particles (pions, kaons, antiprotons)  -> "fireball" bucket
//   stopped net protons (net-proton central/y=0 Gaussian)    -> "participant" bucket
//   forward-rapidity net protons (net-proton side Gaussians) -> "spectator" bucket
// This is the BULK (integrated-over-rapidity, integrated-over-pT) normalization step:
// it produces one number per bucket per centrality bin, in both raw per-event yield
// (dN/dy, summed over the relevant species) and actual cross section (mb, via the
// Woods-Saxon Glauber sigma_reaction already computed in
// GenerateGlauberMC_OO200.C, the C++/ROOT macro that superseded the original
// EstimateNpartNspecFractions_OO200.py this constant was first established
// in). It is the direct input for PlotDoseFractionComparison_OO200.C's
// stacked-bar comparison chart (the Python version of that plotting macro is
// likewise retired).
//
// This is DELIBERATELY SEPARATE from ExtractDoubleDifferentialCrossSections_OO200.C,
// which does the much more involved (mT-m0, y) -> (E, theta) Jacobian conversion
// needed for an actual d^2 sigma/dE dOmega handoff to a transport code. That full
// differential treatment is only possible for the fireball species and the raw
// proton/antiproton spectra (they have a measured pT shape at each rapidity bin to
// convert). Net-proton's participant/spectator split does NOT have that -- it comes
// from a rapidity-only Gaussian decomposition with no associated pT information, so
// this macro reports participant/spectator as BULK numbers only, never claims a
// double-differential for them. See that macro's header for the Fermi-motion-based
// alternative treatment those two components would need.
//
// THE PAIR-PRODUCTION ACCOUNTING (READ BEFORE TRUSTING THE FIREBALL NUMBER): naively
// summing raw yields as fireball = pi+ + pi- + K+ + K- + pbar double-counts net-zero
// against the participant bucket. Net-proton = p - pbar is exactly what removes
// pair-produced (fireball) protons from the "stopped participant" bucket -- so if raw
// antiprotons are ALSO placed in the fireball bucket, antiproton dose cancels out of
// the combined total: fireball + participant = (pi+K+pbar) + (p-pbar) = pi+K+p, and
// the pbar term vanishes entirely. Antiprotons are real, dose-depositing particles
// (they eventually annihilate, typically depositing MORE energy per particle than a
// proton) -- they should not net to zero. The internally consistent fix, used below:
// assign BOTH members of each produced p-pbar pair to the fireball bucket. Since pair
// production makes roughly equal numbers of protons and antiprotons, the pair-produced
// proton content is approximated by the measured antiproton yield, so:
//     fireball_dNdy = dNdy(pi+) + dNdy(pi-) + dNdy(K+) + dNdy(K-) + 2 * dNdy(pbar)
// and participant stays net-proton = p - pbar as already computed. This 2x factor is
// the whole reason this macro exists rather than just summing
// ExtractNetProtonStoppingFraction_OO200.C's CSV with a naive raw-antiproton term.
//
// A SECOND CAVEAT -- ACCEPTANCE MISMATCH BETWEEN BUCKETS: the participant/spectator
// numbers come from ExtractNetProtonStoppingFraction_OO200.C's fit, EXTRAPOLATED all
// the way to y_beam ~ 5.36 (that macro's whole point). The fireball number here is
// NOT extrapolated -- it is the raw dNdy_*_Cent%02d_Nominal histogram integral over
// only the MEASURED acceptance, |y| < 1.55 (see globalDefinitions.h). Particle
// PRODUCTION is expected to fall off faster in rapidity than beam-rapidity-peaked
// spectator baryons do, so this under-count is likely smaller in relative terms for
// pions/kaons/antiprotons than it would be for an unextrapolated baryon measurement --
// but it is still an apples-to-oranges comparison as computed here, not a rigorous
// one. Extrapolating the fireball species the same way (their own rapidity-shape fit,
// analogous to net-proton's Gaussian but for particle PRODUCTION, which typically
// follows a single midrapidity plateau/Gaussian rather than a double-humped shape) is
// a natural next step, not yet done here. Read the "fireball" fraction below as a
// measured-acceptance first-pass number, not a beam-rapidity-integrated one.
//
// INPUTS:
//   a_fitOutputFile : SpectraFitter::writeOutputs() ROOT file (same as
//                      ExtractNetProtonStoppingFraction_OO200.C's input). Reads:
//        ParticleYields/PionPlus/dNdy_PionPlus_Cent%02d_Nominal
//        ParticleYields/PionMinus/dNdy_PionMinus_Cent%02d_Nominal
//        ParticleYields/KaonPlus/dNdy_KaonPlus_Cent%02d_Nominal
//        ParticleYields/KaonMinus/dNdy_KaonMinus_Cent%02d_Nominal
//        ParticleYields/ProtonMinus/dNdy_ProtonMinus_Cent%02d_Nominal  (antiproton)
//     (naming confirmed the same way as the net-proton macro: ParticleInfo::
//      GetParticleName(0,+-1)="PionPlus"/"PionMinus", (1,+-1)="KaonPlus"/"KaonMinus",
//      particleSpeciesName[] confirmed directly in submodule/ParticleInfo/ParticleInfo/
//      ParticleInfo.cxx: index 0="Pion", 1="Kaon", 2="Proton".)
//   a_netProtonCSV  : the CSV written by ExtractNetProtonStoppingFraction_OO200.C
//                      (NetProtonStopping_Summary.csv) -- run that macro FIRST.
//   a_sigmaReactionMb : Woods-Saxon Glauber sigma_reaction, mb. Default below (1173.0)
//                      matches EstimateNpartNspecFractions_OO200.py's Woods-Saxon
//                      result at the time this macro was written -- re-check that
//                      script's printed value hasn't changed (e.g. different RNG seed
//                      or event count) before trusting the default.
//
// OUTPUT: a_outputCSV (default DoseFractionSummary.csv), one row per centrality bin:
//   centIndex, centLabel, sigma_bin_mb,
//   fireball_dNdy(+-err), participant_dNdy(+-err), spectator_dNdy(+-err), total_dNdy,
//   fireball_frac(+-err), participant_frac(+-err), spectator_frac(+-err),
//   fireball_sigma_mb, participant_sigma_mb, spectator_sigma_mb
// Fractions are normalization-independent (sigma_bin cancels in the ratio) -- they are
// meaningful even if sigma_reaction itself later needs revision. Only the *_sigma_mb
// absolute columns depend on getting a_sigmaReactionMb right.
//
// USAGE:
//   root -l -b -q 'ExtractDoseFractionSummary_OO200.C("fit_output.root","NetProtonStopping_Summary.csv")'
//
// UNTESTED, SAME CAVEAT AS THE OTHER TWO SCRIPTS IN THIS DIRECTORY: no ROOT in this
// sandbox, so this could not be run against real data. Checked for brace/paren
// balance and cross-referenced every histogram/column name against source this
// session (see comments above); still worth a careful read before trusting output.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

using namespace std;

const int N_CENT_BINS = 6;
const char* CENT_LABELS[N_CENT_BINS] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
// Matches EstimateNpartNspecFractions_OO200.py's CENT_EDGES = [0,.05,.10,.20,.40,.80,1.00]
const double CENT_BIN_WIDTH_FRAC[N_CENT_BINS] = {0.05, 0.05, 0.10, 0.20, 0.40, 0.20};

struct NetProtonRow {
  int centIndex;
  double integralTotal, integralTotalErr;
  double fracParticipant, fracParticipantErr;
  double fracSpectator, fracSpectatorErr;
  bool found;
};

// Minimal CSV parser for NetProtonStopping_Summary.csv's fixed column layout:
// centIndex,centLabel,model,chi2ndfDouble,chi2ndfTriple,pValueFTest,
// integralTotal,integralTotalErr,fracParticipant,fracParticipantErr,
// fracSpectator,fracSpectatorErr
vector<NetProtonRow> readNetProtonCSV(const string& a_path){
  vector<NetProtonRow> rows;
  ifstream in(a_path.c_str());
  if(!in.is_open()){
    cerr << "ERROR: could not open net-proton CSV '" << a_path << "'" << endl;
    return rows;
  }
  string line;
  bool firstLine = true;
  while(getline(in, line)){
    if(firstLine){ firstLine = false; continue; } // skip header
    if(line.empty()) continue;
    vector<string> fields;
    stringstream ss(line);
    string field;
    while(getline(ss, field, ',')) fields.push_back(field);
    if(fields.size() < 12) continue; // malformed/short row -- skip
    NetProtonRow row;
    row.centIndex          = atoi(fields[0].c_str());
    row.integralTotal       = atof(fields[6].c_str());
    row.integralTotalErr    = atof(fields[7].c_str());
    row.fracParticipant     = atof(fields[8].c_str());
    row.fracParticipantErr  = atof(fields[9].c_str());
    row.fracSpectator       = atof(fields[10].c_str());
    row.fracSpectatorErr    = atof(fields[11].c_str());
    row.found = true;
    rows.push_back(row);
  }
  return rows;
}

// Helper: integral + error over a histogram's full range (measured acceptance only,
// |y| < 1.55 -- see header caveat on the fireball/participant-spectator acceptance
// mismatch).
bool integralAndError(TFile* a_file, const TString& a_histName, double& a_val, double& a_err){
  TH1D* h = (TH1D*)a_file->Get(a_histName);
  if(!h){
    cerr << "WARNING: missing histogram '" << a_histName << "'" << endl;
    a_val = 0; a_err = 0;
    return false;
  }
  a_val = h->IntegralAndError(1, h->GetNbinsX(), a_err);
  return true;
}

void ExtractDoseFractionSummary_OO200(string a_fitOutputFile,
                                       string a_netProtonCSV,
                                       double a_sigmaReactionMb = 1173.0, // Woods-Saxon Glauber default -- see header
                                       string a_outputCSV = "DoseFractionSummary.csv"){

  cout << "==================================================================" << endl;
  cout << " ExtractDoseFractionSummary_OO200" << endl;
  cout << " sigma_reaction (input)  = " << a_sigmaReactionMb << " mb" << endl;
  cout << "==================================================================" << endl;

  TFile* inFile = TFile::Open(a_fitOutputFile.c_str(), "READ");
  if(!inFile || inFile->IsZombie()){
    cerr << "ERROR: could not open input file '" << a_fitOutputFile << "'" << endl;
    return;
  }

  vector<NetProtonRow> netProtonRows = readNetProtonCSV(a_netProtonCSV);
  if(netProtonRows.empty()){
    cerr << "ERROR: no rows read from net-proton CSV -- run "
         << "ExtractNetProtonStoppingFraction_OO200.C first and pass its CSV here." << endl;
    inFile->Close();
    return;
  }

  ofstream csv(a_outputCSV.c_str());
  csv << "centIndex,centLabel,sigma_bin_mb,"
      << "fireball_dNdy,fireball_dNdy_err,participant_dNdy,participant_dNdy_err,"
      << "spectator_dNdy,spectator_dNdy_err,total_dNdy,"
      << "fireball_frac,fireball_frac_err,participant_frac,participant_frac_err,"
      << "spectator_frac,spectator_frac_err,"
      << "fireball_sigma_mb,participant_sigma_mb,spectator_sigma_mb\n";

  for(int centIndex = 0; centIndex < N_CENT_BINS; centIndex++){

    // ---- Fireball species: pi+, pi-, K+, K-, and 2x antiproton (see pair-production caveat) ----
    double pionPlus, pionPlusErr, pionMinus, pionMinusErr;
    double kaonPlus, kaonPlusErr, kaonMinus, kaonMinusErr;
    double antiproton, antiprotonErr;

    bool ok = true;
    ok &= integralAndError(inFile, Form("ParticleYields/PionPlus/dNdy_PionPlus_Cent%02d_Nominal", centIndex), pionPlus, pionPlusErr);
    ok &= integralAndError(inFile, Form("ParticleYields/PionMinus/dNdy_PionMinus_Cent%02d_Nominal", centIndex), pionMinus, pionMinusErr);
    ok &= integralAndError(inFile, Form("ParticleYields/KaonPlus/dNdy_KaonPlus_Cent%02d_Nominal", centIndex), kaonPlus, kaonPlusErr);
    ok &= integralAndError(inFile, Form("ParticleYields/KaonMinus/dNdy_KaonMinus_Cent%02d_Nominal", centIndex), kaonMinus, kaonMinusErr);
    ok &= integralAndError(inFile, Form("ParticleYields/ProtonMinus/dNdy_ProtonMinus_Cent%02d_Nominal", centIndex), antiproton, antiprotonErr);

    if(!ok){
      cerr << "WARNING: centIndex " << centIndex << " (" << CENT_LABELS[centIndex]
           << ") -- missing one or more fireball species histograms, skipping this bin." << endl;
      continue;
    }

    double fireball_dNdy = pionPlus + pionMinus + kaonPlus + kaonMinus + 2.0 * antiproton;
    double fireball_dNdy_err = sqrt(pionPlusErr*pionPlusErr + pionMinusErr*pionMinusErr
                                   + kaonPlusErr*kaonPlusErr + kaonMinusErr*kaonMinusErr
                                   + 4.0*antiprotonErr*antiprotonErr); // quadrature; treats species as uncorrelated (approximate)

    // ---- Participant/spectator: from the net-proton stopping-fraction fit ----
    NetProtonRow* netProtonRow = NULL;
    for(size_t i = 0; i < netProtonRows.size(); i++){
      if(netProtonRows[i].centIndex == centIndex){ netProtonRow = &netProtonRows[i]; break; }
    }
    if(!netProtonRow){
      cerr << "WARNING: centIndex " << centIndex << " (" << CENT_LABELS[centIndex]
           << ") -- not found in net-proton CSV, skipping this bin." << endl;
      continue;
    }

    double participant_dNdy    = netProtonRow->integralTotal * netProtonRow->fracParticipant;
    double spectator_dNdy       = netProtonRow->integralTotal * netProtonRow->fracSpectator;
    // Propagate errors treating integralTotal and the fraction as independent (approximate --
    // in reality they come from the same fit and are correlated; a fully rigorous propagation
    // would need the fit's covariance matrix, not just the two marginal errors already in the CSV).
    double participant_dNdy_err = fabs(participant_dNdy) * sqrt(
        pow(netProtonRow->integralTotalErr / TMath::Max(fabs(netProtonRow->integralTotal),1e-9), 2)
      + pow(netProtonRow->fracParticipantErr / TMath::Max(netProtonRow->fracParticipant,1e-9), 2));
    double spectator_dNdy_err = fabs(spectator_dNdy) * sqrt(
        pow(netProtonRow->integralTotalErr / TMath::Max(fabs(netProtonRow->integralTotal),1e-9), 2)
      + pow(netProtonRow->fracSpectatorErr / TMath::Max(netProtonRow->fracSpectator,1e-9), 2));

    double total_dNdy = fireball_dNdy + participant_dNdy + spectator_dNdy;
    if(total_dNdy <= 0){
      cerr << "WARNING: centIndex " << centIndex << " -- total dNdy <= 0, skipping." << endl;
      continue;
    }

    double fireball_frac     = fireball_dNdy / total_dNdy;
    double participant_frac  = participant_dNdy / total_dNdy;
    double spectator_frac     = spectator_dNdy / total_dNdy;
    // Fractional errors approximated from each bucket's own error relative to the total
    // (ignores anti-correlation between buckets induced by the shared denominator -- a
    // reasonable first-pass estimate, not an exact propagation).
    double fireball_frac_err     = fireball_dNdy_err / total_dNdy;
    double participant_frac_err  = participant_dNdy_err / total_dNdy;
    double spectator_frac_err     = spectator_dNdy_err / total_dNdy;

    double sigma_bin_mb = CENT_BIN_WIDTH_FRAC[centIndex] * a_sigmaReactionMb;
    double fireball_sigma_mb     = sigma_bin_mb * fireball_dNdy;
    double participant_sigma_mb  = sigma_bin_mb * participant_dNdy;
    double spectator_sigma_mb     = sigma_bin_mb * spectator_dNdy;

    cout << "------------------------------------------------------------------" << endl;
    cout << " Cent " << CENT_LABELS[centIndex] << ": fireball=" << 100.0*fireball_frac
         << "%  participant=" << 100.0*participant_frac << "%  spectator=" << 100.0*spectator_frac << "%" << endl;

    csv << centIndex << "," << CENT_LABELS[centIndex] << "," << sigma_bin_mb << ","
        << fireball_dNdy << "," << fireball_dNdy_err << ","
        << participant_dNdy << "," << participant_dNdy_err << ","
        << spectator_dNdy << "," << spectator_dNdy_err << "," << total_dNdy << ","
        << fireball_frac << "," << fireball_frac_err << ","
        << participant_frac << "," << participant_frac_err << ","
        << spectator_frac << "," << spectator_frac_err << ","
        << fireball_sigma_mb << "," << participant_sigma_mb << "," << spectator_sigma_mb << "\n";
  }

  csv.close();
  inFile->Close();
  cout << "==================================================================" << endl;
  cout << "Wrote " << a_outputCSV << endl;
  cout << "Feed this into PlotDoseFractionComparison_OO200.C for the stacked-bar comparison "
       << "(the Python version of that macro is retired)." << endl;
}
