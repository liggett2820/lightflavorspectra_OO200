// EstimateEnergyDoseFraction_OO200_LiteratureStopping.C -- lightflavorspectra_OO200
//
// PURPOSE: the fireball/participant/spectator split Andrew originally asked
// for ("fractional energy dosage of spectators, participants, and produced
// particles by centrality"), computed by ENERGY for the first time this
// session -- everything before this macro (EstimateFireballFractionFromRefMult_OO200.C,
// PlotFireballFractionFromRefMult_OO200.C) was an explicitly-labeled
// particle-COUNT proxy, one unit per track/nucleon regardless of energy.
//
// THE ENERGY MODEL (Andrew's, this session): a straightforward energy-
// conservation split per nucleon in the collision.
//   - a SPECTATOR nucleon never interacts -> keeps its full beam energy.
//   - a PARTICIPANT nucleon interacts -> keeps (beam energy - stopping loss);
//     the stopping loss it doesn't keep has to go somewhere by energy
//     conservation, and the standard reading of a "nuclear stopping"
//     measurement like the one cited below IS that this lost energy converts
//     into newly PRODUCED particles (the fireball) -- so:
//   - PRODUCED-PARTICLE (fireball) energy = Npart * (per-participant stopping
//     loss), i.e. exactly the energy participants gave up.
// This means E_participant_total + E_produced_total = Npart * E_beam always
// (the split within a participant's energy budget, not a separate pool), and
// the whole-event energy budget E_spectator+E_participant+E_produced =
// (Nspec+Npart)*E_beam = 2*A*E_beam is CONSTANT across centrality bins (a
// built-in consistency check printed below -- should read 3200 GeV for every
// bin, since 2*A=32 nucleons * 100 GeV each for O+O at 200 GeV).
//
// INPUTS:
//   Npart, Nspec per centrality bin: read live from a_glauberRoot (default
//   glauber_oo200_results_NBDfit.root), same reader/schema as
//   EstimateFireballFractionFromRefMult_OO200.C -- real, NBD-Glauber-fit-
//   calibrated O+O numbers, not literature.
//
//   E_beam per nucleon = sqrt(s_NN)/2 = 100 GeV -- standard heavy-ion
//   kinematic definition, same sqrt(s_NN)=200 GeV this whole analysis uses.
//
//   E_loss per participant = 73 +/- 6 GeV -- BRAHMS, "Nuclear Stopping in
//   Au+Au Collisions at sqrt(s_NN)=200 GeV," Phys. Rev. Lett. 93, 102301
//   (2004) (via arXiv:nucl-ex/0312023, https://arxiv.org/abs/nucl-ex/0312023),
//   same paper already cited in EstimateBaryonStoppingFromLiterature_AuAu200.C
//   for delta_y. CENTRAL (0-5%) Au+Au only in the paper -- applied here to
//   ALL SIX O+O centrality bins under Andrew's explicit instruction (this
//   session, following EstimateBaryonStoppingFromLiterature_AuAu200.C) to
//   assume stopping is centrality-independent at this energy.
//
// CAVEATS (same two as EstimateBaryonStoppingFromLiterature_AuAu200.C, worth
// repeating here since they now bias the ENERGY split, not a particle-count
// split):
//   - SYSTEM SIZE: 73 GeV is measured in Au+Au (A=197), not O+O (A=16).
//     Lighter/thinner systems are expected to show LESS stopping, i.e. a
//     participant nucleon likely keeps MORE than (100-73)=27 GeV in O+O, not
//     less. That means this estimate likely OVERSTATES fireball_frac (too
//     much energy assumed converted to produced particles) and UNDERSTATES
//     participant_frac (too little energy assumed retained by participants).
//     spectator_frac is comparatively robust -- it's just Nspec*100 GeV,
//     independent of the stopping-loss assumption entirely.
//   - CENTRALITY: only 0-5% has a direct literature value; the other five
//     bins carry it under the explicit centrality-independence assumption,
//     not an independent measurement.
//
// OUTPUT:
//   Console table: Npart, Nspec, E_spectator, E_participant, E_produced (GeV)
//   and the three fractions, per centrality bin, plus the constant-budget
//   consistency check.
//   OO200_EnergyDoseFractionFromLiteratureStopping.png -- stacked bar, same
//   color convention as PlotFireballFractionFromRefMult_OO200.C /
//   PlotDoseFractionComparison_OO200.C (brick red fireball, navy participant,
//   amber spectator) for visual continuity across this repo's figures. Now
//   carries a 7th "Min-Bias" bar -- see CROSS-SECTION-WEIGHTED MIN-BIAS
//   AVERAGE below.
//
// CROSS-SECTION-WEIGHTED MIN-BIAS AVERAGE (added after Andrew flagged that a
// single combined number needs to weight each bin by how much of the total
// cross section it actually covers, not just average the six printed
// percentages straight across):
//   The six centrality bins are NOT equal width -- 0-5%/5-10% are 5 points
//   wide, 10-20% is 10, 20-40% is 20, 40-80% is 40, 80-100% is 20 (same
//   CENT_EDGES_FRAC boundaries FitGlauberNBDToRefMult_OO200.C uses to derive
//   Npart/Nspec per bin, duplicated here -- keep the two in sync if either
//   changes). An unweighted mean of the six fireball/participant/spectator
//   percentages would silently let the 40-80% bin (40% of all collisions)
//   count the same as 0-5% (5% of all collisions) -- wrong for a "what does
//   a random min-bias O+O collision look like" summary number. Instead each
//   bin's fraction is weighted by its own width in percentile terms
//   (binWidth[i] = CENT_EDGES_FRAC[i+1]-CENT_EDGES_FRAC[i], summing to 1.0
//   over all six bins), which is the same event-share assumption the
//   percentile-cut centrality definition itself is built on. This is a
//   cross-section-share weighting of the model's own per-bin outputs, not an
//   independent measurement -- it inherits every caveat already attached to
//   fireballFrac/participantFrac/spectatorFrac above (the Au+Au stopping
//   number applied to O+O, foremost).
//
// USAGE:
//   root -l -b -q 'EstimateEnergyDoseFraction_OO200_LiteratureStopping.C()'
// or with an explicit Glauber ROOT file path:
//   root -l -b -q 'EstimateEnergyDoseFraction_OO200_LiteratureStopping.C("glauber_oo200_results_NBDfit.root")'
//
// STATUS: UNTESTED -- no ROOT in the sandbox that wrote this (same caveat as
// every other C++ macro in this session). Andrew HAS run an earlier
// (JSON-reading) version of this macro (the currently-embedded
// OO200_EnergyDoseFractionFromLiteratureStopping.png / oo200_energy_dose_fraction.png
// reflects real numbers) -- but TWO things below have NOT been exercised
// against real data yet and need a re-run before trusting them: the min-bias
// weighted-average addition (7th bar + console block), and the Glauber-table
// reader, just switched from hand-rolled JSON line-scanning to a TTree read
// (Python/JSON retired from this pipeline, per Andrew's instruction).

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TStyle.h"
#include "TSystem.h"

using namespace std;

//============================================================================
// Palette -- identical hex values to PlotDoseFractionComparison_OO200.C /
// PlotFireballFractionFromRefMult_OO200.C.
//============================================================================
const int COLOR_FIREBALL      = TColor::GetColor("#C0392B"); // brick red
const int COLOR_PARTICIPANT   = TColor::GetColor("#1E2761"); // navy
const int COLOR_SPECTATOR     = TColor::GetColor("#E8A33D"); // amber
const int COLOR_TEXT_ON_DARK  = kWhite;
const int COLOR_TEXT_ON_LIGHT = TColor::GetColor("#1A1A1A");
const int COLOR_CAVEAT        = TColor::GetColor("#7A7A7A");

const int N_CENT_BINS = 6;
const char* CENT_LABELS[N_CENT_BINS] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
// Percentile bin edges -- IDENTICAL to FitGlauberNBDToRefMult_OO200.C's
// CENT_EDGES_FRAC (kept in sync manually, same policy as the species-weight
// table shared with GCRSpectrum.hh). Used only to weight the min-bias
// average below by each bin's true share of the total cross section.
const double CENT_EDGES_FRAC[N_CENT_BINS + 1] = {0.0, 0.05, 0.10, 0.20, 0.40, 0.80, 1.00};

// Cited inputs -- see header.
const double E_BEAM_PER_NUCLEON_GEV = 100.0; // sqrt(s_NN)/2, sqrt(s_NN)=200 GeV
const double E_LOSS_PER_PARTICIPANT_GEV = 73.0; // BRAHMS, central Au+Au 200 GeV
const double E_LOSS_PER_PARTICIPANT_ERR_GEV = 6.0;

//============================================================================
// Glauber centrality-bin reader -- reads the shared "GlauberCentralityBins"
// TTree schema (see GenerateGlauberMC_OO200.C's header), trimmed to just
// cent/mean_Npart/mean_Nspec (no refMult edge needed here). Replaces the
// former hand-rolled JSON line-scanning reader now that Python/JSON are
// retired from this pipeline (same idiom/scope as
// EstimateFireballFractionFromRefMult_OO200.C's ReadGlauberCentralityTable()).
//============================================================================
struct GlauberBinRow {
  string cent;
  double meanNpart = 0.0, meanNspec = 0.0;
  bool foundNpart = false, foundNspec = false;
};

vector<GlauberBinRow> ReadGlauberCentralityTable(const string& a_path){
  vector<GlauberBinRow> out;
  TFile* f = TFile::Open(a_path.c_str(), "READ");
  if(!f || f->IsZombie()) return out;

  TTree* t = (TTree*) f->Get("GlauberCentralityBins");
  if(!t){ f->Close(); return out; }

  Char_t cent_in[16];
  Double_t meanNpart_in, meanNspec_in;
  t->SetBranchAddress("cent", cent_in);
  t->SetBranchAddress("meanNpart", &meanNpart_in);
  t->SetBranchAddress("meanNspec", &meanNspec_in);

  Long64_t nEntries = t->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    t->GetEntry(i);
    GlauberBinRow row;
    row.cent = string(cent_in);
    row.meanNpart = meanNpart_in;
    row.meanNspec = meanNspec_in;
    row.foundNpart = true;
    row.foundNspec = true;
    out.push_back(row);
  }
  f->Close();
  return out;
}

void DrawPercentLabel(double x, double y, double val, int color, double textSize = 0.028){
  if(val < 3.0) return;
  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextSize(textSize);
  lat.SetTextAlign(22);
  lat.SetTextColor(color);
  lat.DrawLatex(x, y, Form("%.0f%%", val));
}

//============================================================================
// Main macro
//============================================================================
void EstimateEnergyDoseFraction_OO200_LiteratureStopping(
    string a_glauberRoot = "glauber_oo200_results_NBDfit.root",
    string a_outputDir   = "."
){
  gStyle->SetOptStat(0);
  gSystem->mkdir(a_outputDir.c_str(), kTRUE);

  //--------------------------------------------------------------------------
  // Load Npart/Nspec per bin -- required, no fabrication if missing.
  //--------------------------------------------------------------------------
  vector<GlauberBinRow> rows = ReadGlauberCentralityTable(a_glauberRoot);
  if(rows.empty()){
    cerr << "ERROR: could not read any centrality-bin rows from '" << a_glauberRoot << "'."
         << " Run FitGlauberNBDToRefMult_OO200.C (or GenerateGlauberMC_OO200.C for the purely"
         << " geometric file) first, or pass the right path." << endl;
    return;
  }

  double meanNpart[N_CENT_BINS], meanNspec[N_CENT_BINS];
  for(int i = 0; i < N_CENT_BINS; i++){
    const GlauberBinRow* match = NULL;
    for(size_t j = 0; j < rows.size(); j++){
      if(rows[j].cent == CENT_LABELS[i]){ match = &rows[j]; break; }
    }
    if(!match || !match->foundNpart || !match->foundNspec){
      cerr << "ERROR: '" << a_glauberRoot << "' has no mean_Npart/mean_Nspec for cent bin '"
           << CENT_LABELS[i] << "' -- aborting rather than guessing." << endl;
      return;
    }
    meanNpart[i] = match->meanNpart;
    meanNspec[i] = match->meanNspec;
  }

  cout << "==================================================================" << endl;
  cout << "EstimateEnergyDoseFraction_OO200_LiteratureStopping" << endl;
  cout << "O+O sqrt(s_NN)=200 GeV -- real Npart/Nspec (NBD-Glauber-fit), literature stopping" << endl;
  cout << "E_beam=" << E_BEAM_PER_NUCLEON_GEV << " GeV/nucleon   E_loss(participant)="
       << E_LOSS_PER_PARTICIPANT_GEV << " +/- " << E_LOSS_PER_PARTICIPANT_ERR_GEV
       << " GeV (BRAHMS, central Au+Au -- applied to all bins, see header)" << endl;
  cout << "==================================================================" << endl;

  //--------------------------------------------------------------------------
  // Energy accounting per bin.
  //--------------------------------------------------------------------------
  double eSpectator[N_CENT_BINS], eParticipant[N_CENT_BINS], eProduced[N_CENT_BINS], eTotal[N_CENT_BINS];
  double fireballFrac[N_CENT_BINS], participantFrac[N_CENT_BINS], spectatorFrac[N_CENT_BINS];

  printf("%-10s%10s%10s%14s%14s%14s%12s\n", "Cent", "Npart", "Nspec", "E_spec[GeV]", "E_part[GeV]", "E_prod[GeV]", "E_tot[GeV]");
  for(int i = 0; i < N_CENT_BINS; i++){
    eSpectator[i]   = meanNspec[i] * E_BEAM_PER_NUCLEON_GEV;
    eParticipant[i] = meanNpart[i] * (E_BEAM_PER_NUCLEON_GEV - E_LOSS_PER_PARTICIPANT_GEV);
    eProduced[i]    = meanNpart[i] * E_LOSS_PER_PARTICIPANT_GEV;
    eTotal[i]       = eSpectator[i] + eParticipant[i] + eProduced[i];

    fireballFrac[i]    = eProduced[i]    / eTotal[i];
    participantFrac[i] = eParticipant[i] / eTotal[i];
    spectatorFrac[i]   = eSpectator[i]   / eTotal[i];

    printf("%-10s%10.2f%10.2f%14.1f%14.1f%14.1f%12.1f\n", CENT_LABELS[i], meanNpart[i], meanNspec[i],
           eSpectator[i], eParticipant[i], eProduced[i], eTotal[i]);
  }

  cout << endl << "Consistency check: E_tot should equal 2*A*E_beam = 32*100 = 3200 GeV in EVERY"
       << " bin (Npart+Nspec=32 always) -- if any row above deviates from 3200, something's off." << endl;

  cout << endl;
  printf("%-10s%16s%16s%16s\n", "Cent", "fireball %", "participant %", "spectator %");
  for(int i = 0; i < N_CENT_BINS; i++){
    printf("%-10s%16.1f%16.1f%16.1f\n", CENT_LABELS[i], 100.0*fireballFrac[i], 100.0*participantFrac[i], 100.0*spectatorFrac[i]);
  }

  //--------------------------------------------------------------------------
  // Cross-section-weighted min-bias average -- see header note above. Each
  // bin's fraction is weighted by its own percentile width, not counted
  // once-per-bin, so the wide 40-80% bin properly outweighs the narrow 0-5%
  // bin in the combined number.
  //--------------------------------------------------------------------------
  double binWidth[N_CENT_BINS], sumBinWidth = 0.0;
  for(int i = 0; i < N_CENT_BINS; i++){
    binWidth[i] = CENT_EDGES_FRAC[i + 1] - CENT_EDGES_FRAC[i];
    sumBinWidth += binWidth[i];
  }
  if(TMath::Abs(sumBinWidth - 1.0) > 1e-9){
    cout << "WARNING: CENT_EDGES_FRAC bin widths sum to " << sumBinWidth
         << ", not 1.0 -- check the edges array before trusting the min-bias average." << endl;
  }

  double minBiasFireballFrac = 0.0, minBiasParticipantFrac = 0.0, minBiasSpectatorFrac = 0.0;
  for(int i = 0; i < N_CENT_BINS; i++){
    minBiasFireballFrac    += binWidth[i] * fireballFrac[i];
    minBiasParticipantFrac += binWidth[i] * participantFrac[i];
    minBiasSpectatorFrac   += binWidth[i] * spectatorFrac[i];
  }

  cout << endl << "Cross-section-weighted min-bias (0-100%) average, weighted by each bin's"
       << " percentile width (5,5,10,20,40,20):" << endl;
  printf("%-10s%16.1f%16.1f%16.1f\n", "0-100%", 100.0*minBiasFireballFrac,
         100.0*minBiasParticipantFrac, 100.0*minBiasSpectatorFrac);
  cout << "(For comparison, the UNWEIGHTED mean of the six bins above would be "
       << Form("%.1f/%.1f/%.1f", 100.0*(fireballFrac[0]+fireballFrac[1]+fireballFrac[2]+fireballFrac[3]+fireballFrac[4]+fireballFrac[5])/6.0,
               100.0*(participantFrac[0]+participantFrac[1]+participantFrac[2]+participantFrac[3]+participantFrac[4]+participantFrac[5])/6.0,
               100.0*(spectatorFrac[0]+spectatorFrac[1]+spectatorFrac[2]+spectatorFrac[3]+spectatorFrac[4]+spectatorFrac[5])/6.0)
       << " -- wrong for a min-bias summary since it treats all six bins as equally likely.)" << endl;

  //--------------------------------------------------------------------------
  // Plot: stacked bar, same style/colors as the other dose-fraction figures.
  // A 7th "Min-Bias" bar (cross-section-weighted average, see above) is
  // appended after the six centrality bins, visually separated by extra bar
  // spacing so it isn't mistaken for a 7th centrality class.
  //--------------------------------------------------------------------------
  // Layout: 6 centrality bins (bin positions 1-6), then the Min-Bias bar
  // immediately after (position 7, labeled "0-100%") -- no empty spacer
  // bin. The dashed TLine drawn at the position-6/7 boundary below is what
  // marks it as a separate summary rather than a 7th centrality class, so
  // an empty bin isn't needed to make that visual distinction.
  const int N_PLOT_BINS = N_CENT_BINS + 1;
  const int MINBIAS_BIN = N_CENT_BINS + 1; // position 7
  TCanvas* c = new TCanvas("cEnergyDoseFraction", "Energy Dose Fraction (Literature Stopping)", 950, 700);
  gPad->SetBottomMargin(0.24);
  gPad->SetGridy();

  TH1D* hStackTop = new TH1D("hStackTop", "", N_PLOT_BINS, 0, N_PLOT_BINS);
  TH1D* hStackMid = new TH1D("hStackMid", "", N_PLOT_BINS, 0, N_PLOT_BINS);
  TH1D* hStackBot = new TH1D("hStackBot", "", N_PLOT_BINS, 0, N_PLOT_BINS);

  for(int i = 0; i < N_CENT_BINS; i++){
    double fb = 100.0 * fireballFrac[i];
    double pt = 100.0 * participantFrac[i];
    double sp = 100.0 * spectatorFrac[i];
    hStackTop->SetBinContent(i + 1, fb + pt + sp);
    hStackMid->SetBinContent(i + 1, fb + pt);
    hStackBot->SetBinContent(i + 1, fb);
    hStackTop->GetXaxis()->SetBinLabel(i + 1, CENT_LABELS[i]);
  }
  {
    double fb = 100.0 * minBiasFireballFrac;
    double pt = 100.0 * minBiasParticipantFrac;
    double sp = 100.0 * minBiasSpectatorFrac;
    hStackTop->SetBinContent(MINBIAS_BIN, fb + pt + sp);
    hStackMid->SetBinContent(MINBIAS_BIN, fb + pt);
    hStackBot->SetBinContent(MINBIAS_BIN, fb);
    hStackTop->GetXaxis()->SetBinLabel(MINBIAS_BIN, "0-100%");
  }

  double barWidth = 0.6, barOffset = (1.0 - barWidth) / 2.0;
  hStackTop->SetBarWidth(barWidth); hStackTop->SetBarOffset(barOffset);
  hStackMid->SetBarWidth(barWidth); hStackMid->SetBarOffset(barOffset);
  hStackBot->SetBarWidth(barWidth); hStackBot->SetBarOffset(barOffset);

  hStackTop->SetFillColor(COLOR_SPECTATOR);
  hStackMid->SetFillColor(COLOR_PARTICIPANT);
  hStackBot->SetFillColor(COLOR_FIREBALL);

  hStackTop->SetMaximum(100.0);
  hStackTop->SetMinimum(0.0);
  hStackTop->GetXaxis()->SetLabelSize(0.040);
  hStackTop->GetYaxis()->SetTitle("Fraction of Total Kinetic Energy [%]");
  hStackTop->GetXaxis()->SetTitle("Centrality");
  // ROOT right-aligns axis titles by default, which put "Centrality" at the
  // far right of the axis -- exactly where the Min-Bias bin's label lives.
  // That's the real cause of the overlap (independent of how that label is
  // formatted); center the title instead so it sits under the middle of the
  // six centrality bins, clear of the Min-Bias bin entirely.
  hStackTop->GetXaxis()->CenterTitle(kTRUE);
  hStackTop->SetTitle("O+O #sqrt{s_{NN}}=200 GeV -- Toy Energy Split vs. Centrality (Literature Stopping Loss)");

  hStackTop->Draw("bar");
  hStackMid->Draw("bar same");
  hStackBot->Draw("bar same");

  // Separator line between the six centrality bins and the Min-Bias summary
  // bar, at the position-6/7 boundary.
  TLine* sepLine = new TLine((double) N_CENT_BINS, 0.0, (double) N_CENT_BINS, 100.0);
  sepLine->SetLineColor(TColor::GetColor("#B8B8B8"));
  sepLine->SetLineStyle(2);
  sepLine->SetLineWidth(1);
  sepLine->Draw();

  for(int i = 0; i < N_CENT_BINS; i++){
    double x = hStackTop->GetXaxis()->GetBinCenter(i + 1);
    double fb = 100.0 * fireballFrac[i];
    double pt = 100.0 * participantFrac[i];
    double sp = 100.0 * spectatorFrac[i];
    DrawPercentLabel(x, fb / 2.0, fb, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt / 2.0, pt, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt + sp / 2.0, sp, COLOR_TEXT_ON_LIGHT);
  }
  {
    double x = hStackTop->GetXaxis()->GetBinCenter(MINBIAS_BIN);
    double fb = 100.0 * minBiasFireballFrac;
    double pt = 100.0 * minBiasParticipantFrac;
    double sp = 100.0 * minBiasSpectatorFrac;
    DrawPercentLabel(x, fb / 2.0, fb, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt / 2.0, pt, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt + sp / 2.0, sp, COLOR_TEXT_ON_LIGHT);
  }

  TLegend* leg = new TLegend(0.12, 0.08, 0.90, 0.16);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hStackBot, "Fireball-produced (energy lost by participants)", "f");
  leg->AddEntry(hStackMid, "Participant-retained (E_beam - stopping loss)", "f");
  leg->AddEntry(hStackTop, "Spectator (full E_beam, untouched)", "f");
  leg->Draw();

  TLatex caveat;
  caveat.SetNDC();
  caveat.SetTextFont(42);
  caveat.SetTextSize(0.022);
  caveat.SetTextColor(COLOR_CAVEAT);
  caveat.SetTextAlign(21);
  caveat.DrawLatex(0.5, 0.02,
    "Stopping loss (73#pm6 GeV/participant) from BRAHMS Au+Au 200 GeV, applied to all O+O bins -- likely overstates fireball_frac for a lighter system.");

  string outPath = a_outputDir + "/OO200_EnergyDoseFractionFromLiteratureStopping.png";
  c->SaveAs(outPath.c_str());
  cout << endl << "Wrote " << outPath << endl;

  cout << endl << "Reminder: E_loss=73 GeV is an Au+Au (not O+O) measurement applied under an"
       << " explicit centrality-independence assumption -- treat fireball_frac here as likely"
       << " an OVERESTIMATE and participant_frac as likely an UNDERESTIMATE for the real O+O"
       << " system (spectator_frac is not sensitive to this assumption). See header." << endl;
}
