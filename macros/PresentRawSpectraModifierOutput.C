// macros/PresentRawSpectraModifierOutput.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS (2026-07-21):
// Andrew asked to see the corrected pion spectra written by RunRawSpectraModifier.C
// (see macros/RunRawSpectraModifier.C and the project log's "RunRawSpectraModifier.C
// command for corrected pion spectra" section). There was no macro in this repo that
// reads RawSpectraModifier's OWN output format directly -- macros/PresentZFitterSpectra.C
// reads a ZFitter spectra file (TH2-typed DeDx_FitData/BTOF_FitData objects, one 2D
// histogram per species/charge/centrality, rapidity vs mT-m0), and
// macros/DrawSpectraAndDnDy.C reads SpectraFitter's output (a separate, not-yet-run
// downstream step) -- neither reads RawSpectraModifier::writeRapiditySpectra()'s actual
// output format, which is a TGraphErrors PER (species, charge, centrality, rapidity bin)
// under a "Spectra_<Species><Plus/Minus>/" directory (confirmed via direct read of
// RawSpectraModifier.cxx's writeRapiditySpectra(), ~line 4750-4783).
//
// WHAT'S ACTUALLY IN THE FILE (confirmed via direct read of source/RawSpectraModifier.cxx):
//   - Directory naming: "Spectra_<Species><Plus/Minus>" (e.g. "Spectra_PionPlus"), from
//     ParticleInfo::GetParticleName(partIndex,+-1) -- species name and Plus/Minus are
//     CONCATENATED with no underscore (confirmed via direct read of
//     submodule/ParticleInfo/ParticleInfo/ParticleInfo.cxx's GetParticleName()).
//   - Object naming INSIDE that directory uses a DIFFERENT convention -- species name
//     and Plus/Minus ARE underscore-separated there, because the underscore is
//     hardcoded into the Form(...) string itself rather than coming from
//     GetParticleName(): "corrected_dEdxSpectra_<Species>_<Plus/Minus>_Cent<CC>_yIndex_<YY>"
//     (TPC dE/dx-derived), "corrected_InvBetaBTOFSpectra_<Species>_<Plus/Minus>_Cent<CC>_
//     yIndex_<YY>" (BTOF 1/beta-derived), and "corrected_InvBetaETOFSpectra_..." (ETOF,
//     only present if this repo's ETOF fit drivers were ever exercised -- they weren't,
//     per RunZFitter.C's own header comment, so this will simply be absent and is
//     handled as such below, not treated as an error).
//   - There are also "raw_..." (pre-correction) and (if populateInclusiveSpectra() was
//     called, which RunRawSpectraModifier.C's live call sequence always does unless
//     a_rawOnly=true) "inclusive_..." (corrected minus feed-down) variants of the same
//     objects. This macro only reads the "corrected_..." ones -- the actual final
//     product of a normal (a_rawOnly=false) RunRawSpectraModifier.C run.
//   - Each TGraphErrors' NAME was overwritten to the "corrected_..." convention above,
//     but its TITLE was never touched -- it's still whatever
//     HistogramUtilities::VerticalProjectionGraphs() set it to when the graph was first
//     created: "<original TH2 title> with y=[<low>,<high>]" (confirmed via direct read,
//     source/HistogramUtilities.cxx ~line 4134-4157). This is the ONLY place the actual
//     rapidity value for a given yIndex survives into the output file, so
//     parseRapidityCenterFromTitle() below recovers it from there rather than assuming
//     a fixed bin width/offset that could silently drift out of sync with whatever
//     rapidity binning a given run actually used.
//
// STYLE (revised 2026-07-21, see LATEST REVISION below): originally matched
// macros/PresentZFitterSpectra.C's color-by-|y| convention (single-hue blue ramp,
// hand-drawn color-scale bar). Andrew asked for this macro's plots to instead look
// like a reference STAR Au+Au fixed-target deuteron spectra plot: plain black
// markers with NO color-by-rapidity at all, no color-scale bar, and the
// decade-offset trick expressed as a compact on-plot formula rather than a color
// legend or prose. Marker styles 20/24/25 for TPC/BTOF/ETOF and the log-y DrawFrame
// convention are unchanged from PresentZFitterSpectra.C. The
// centralityLabel()/particleSymbolWithCharge() helper functions below are still
// copied near-verbatim from that macro (this is a separate interpreted
// macro/translation unit, so they can't just be #included from there) -- see that
// macro's own header comments for the full rationale behind each design choice if
// it's not repeated here.
//
// UNLIKE PresentZFitterSpectra.C, this macro (through REVISION 7) did NOT run any of
// its own suppression cuts (significance / chi^2 / local-trend-outlier) -- those were
// ZFitter-output-stage QA cuts applied to raw yields before fitting; RawSpectraModifier's
// output was already downstream of that (and of cleanSpectra()/trimSpectra(),
// RawSpectraModifier's own cleanup step), so this macro just drew whatever points
// survived that pipeline, as-is. REVISION 8 below changes this for the significance cut
// specifically -- see that entry for why.
//
// REVISION 8 (2026-07-23): Andrew asked to bring PresentZFitterSpectra.C's statistical
// significance cut (content/error, i.e. here point-Y/point-EY, below a_minSignificance
// gets dropped -- see that macro's suppressInsignificantBins() and its STATISTICAL
// SIGNIFICANCE CUT header comment for the full tuning history, settled at 5.0 sigma)
// into this macro too, as the first of the mT-m0-cut methodologies discussed for
// finalizing the "Determine mT-m0 cuts" Next Steps item. Added suppressInsignificantPoints()
// below -- functionally the same cut, adapted for this macro's TGraphErrors (rather than
// TH1D) data representation -- applied to each rapidity bin's TPC/BTOF/ETOF graphs right
// after they're fetched from file, in drawOneCorrectedSpectraPanel()'s yIndex loop, so
// both the numeric N_evt/BTOF checks and the drawn plot itself see the same
// already-cut points. New a_minSignificance parameter on PresentRawSpectraModifierOutput()
// defaults to 5.0, matching PresentZFitterSpectra.C's tuned value; pass 0 (or negative)
// to disable, same convention as that macro.
//
// REVISION 9 (2026-07-23, same day): Andrew asked for the chi^2/ndf fit-quality cut too
// (the other half of PresentZFitterSpectra.C's suppression pipeline -- see that macro's
// suppressBadFitQualityBins()). Checked first whether the underlying data even supports
// it: RawSpectraModifier::loadSpectraFile() (source/RawSpectraModifier.cxx ~line 635)
// never loads the ZFitter file's ZTPC_ChiSqr_.../ZbTOF_ChiSqr_... histograms at all, so
// chi^2/ndf never enters RawSpectraModifier's output file -- there's nothing for this
// macro to read there, populated or not. Also confirmed (by reading
// convertSpectraToInvariant() and applyBTOFEff()) that RawSpectraModifier only rescales
// each bin's already-fitted content/error by a scalar (N_evt/efficiency normalization) --
// it never re-fits anything, so a given (rapidity, mT-m0) bin's chi^2/ndf, wherever it's
// read from, is numerically identical to whatever ZFitter originally computed for it.
// Andrew then ran macros/CheckChiSqrValues.C directly against the Proton ZFitter TPC/BTOF
// files (Cent00) to check whether the values are even populated for this species -- UNLIKE
// pions (see PresentZFitterSpectra.C's a_maxChiSqrNdf header note: "this dataset's
// chi^2/ndf histograms turned out to be all zero"), protons come back with real,
// non-degenerate values: TPC Plus median 2.06 (p95=4.46), BTOF Plus median 1.63,
// BTOF Minus median 1.09, all reasonably well-behaved. TPC Minus (antiproton) is the
// exception -- median already 8.16, with the same "confidently wrong fit, deceptively
// small quoted error" failure mode this cut exists to catch: rapBin=9/mtm0Bin=17 has
// yield=0.499 but chi^2/ndf=1.03e13, a fit that clears even a very tight significance cut
// (565947 sigma) despite being obvious garbage. Rather than duplicate the chi^2/ndf
// histograms into RawSpectraModifier's output (a compiled-source change + rebuild, given
// the finding above that it wouldn't be any DIFFERENT data, just a copy), this macro
// instead opens the ORIGINAL ZFitter TPC/BTOF file(s) directly -- same files
// PresentZFitterSpectra.C reads -- and looks up each drawn point's chi^2/ndf by
// COORDINATE (FindBin() on rapidity and mT-m0) rather than by a passed-in bin index, since
// this macro's data is TGraphErrors points, not TH1D bins on a fixed axis. This is valid
// specifically because RawSpectraModifier never rebins these TH2s (see above) -- the
// ZFitter file's axis binning is exactly what the corrected TGraphErrors points were
// projected from. New a_tpcZFitterFile/a_btofZFitterFile parameters (both "" by default,
// i.e. chi^2 cut off) point at those files; a_maxChiSqrNdf (TPC Plus/BTOF Plus/BTOF
// Minus/ETOF) and a_maxChiSqrNdfTPCMinus (TPC Minus only, kept separate given how much
// higher that category's distribution runs) default to 0 (disabled), same convention as
// a_minSignificance and as PresentZFitterSpectra.C's own a_maxChiSqrNdf.
//
// PREREQUISITE: a_correctedSpectraFile is RunRawSpectraModifier.C's a_outputFile
// argument (e.g. "spectra_OO200_Pion_Corrected.root").
//
// USAGE (run from within macros/ or anywhere -- this macro has no ../bin/ dependency,
// it only reads a pre-built ROOT file with plain TGraphErrors, no custom-class
// dictionary needed):
//   root -l -q 'PresentRawSpectraModifierOutput.C("/Users/aliggett/data/OO/zfitter_output/spectra_OO200_Pion_Corrected.root")'
//   root -l -q 'PresentRawSpectraModifierOutput.C("spectra_OO200_Pion_Corrected.root",0,"./corrected_pion_plots")'
//   root -l -q 'PresentRawSpectraModifierOutput.C("spectra_OO200_Proton_Corrected.root",2,"./corrected_proton_plots",true,3.0)'  // e.g. a looser 3-sigma significance cut
//   root -l -q 'PresentRawSpectraModifierOutput.C("spectra_OO200_Proton_Corrected.root",2,"./corrected_proton_plots",true,5.0,"spectra_OO200_Proton_TPC_RapBin_-1_-1_-1_-1_-1_-1.root","spectra_OO200_Proton_BTOF_RapBin_-1_-1_-1_-1_-1_-1.root",7.0,20.0)'  // + chi^2/ndf cut: 7.0 everywhere except TPC Minus (antiprotons), which gets a looser 20.0 -- see REVISION 9 header comment for why
//
// STATUS: written 2026-07-21, validated end-to-end by Andrew (all 12 Pion panels
// render; the N_evt and TPC-point numeric checks print correctly). A subsequent
// revision fixed a misleading claim in printNormalizationCheck()'s own comment (see
// that function for the full story) and added a genuine BTOF-point correction check
// (printBTOFCorrectionCheck()) as a stronger confirmatory test -- both since VALIDATED
// (Andrew's kaon-run check landed in the expected 2.2-2.27 range, matching pion).
//
// REVISION 2 (2026-07-21, same day): a_partIndex already made this macro
// species-agnostic from the start (nothing kaon-specific needed for Andrew's first
// kaon run) -- this revision was purely visual polish, per Andrew's request. Found
// that gStyle->SetOptTitle(0) (set below) silently suppressed EVERY plot's DrawFrame
// title string, including the "rapidity bins offset by decades for readability" note
// that used to live inside it -- so that scaling info, which isn't captured anywhere
// in the y-axis label, was never actually visible on any rendered plot. Fixed by
// drawing the title and a scaling/annotation subtitle explicitly via TLatex instead.
// Andrew's follow-up feedback: "there's too many words" (the 2-line prose subtitle
// was too verbose) and asked for the rapidity coloring to not be symmetric-by-|y|.
//
// REVISION 3 (2026-07-21, same day): Andrew then sent a reference plot (STAR Au+Au
// fixed-target deuteron spectra, sqrt(s_NN)=3.2 GeV) and asked for this macro's
// output to look like it. That reference uses plain BLACK markers with NO
// color-by-rapidity scheme at all (not even an asymmetric one), no color-scale bar,
// and expresses the decade-offset trick as one compact on-plot formula line (e.g.
// "#bullet #times3^{#pm10y}") rather than prose or a color legend. This directly
// answers both pieces of feedback above at once -- removing color-by-rapidity
// entirely is a stronger form of "not symmetric," and the compact formula line
// replaces the verbose subtitle. Implemented here by:
//   (1) removing color-by-|y| entirely -- drawColorScaleBar(), the color-ramp setup
//       in PresentRawSpectraModifierOutput(), and the a_nColors/a_kRapColorMax
//       parameters are all gone; every marker is now plain black (kBlack).
//   (2) replacing the rank-index decade offset, mult = 10^b where b was the bin's
//       position in the bins vector (0..~40 now that the rapidity cut was widened
//       earlier this session -- this could reach mult=10^40, an unintended
//       regression from that widening, never actually exercised because Andrew
//       hadn't re-rendered since), with a bounded, y-VALUE-based exponential --
//       mult = kOffsetBase^(kOffsetScaleFactor * y). This naturally stays bounded
//       since y itself is bounded (roughly +-2.05 for this analysis), fixing the
//       latent bug and mirroring the reference plot's convention. kOffsetScaleFactor
//       was originally tuned down to 2 (smaller than the reference's implied "10")
//       to keep the total offset span modest -- but Andrew found that left too little
//       vertical separation between rapidity-bin curves after the rapidity-cut
//       widening added many more bins, reading as extra background overlap in the
//       lower-right (high mT-m0) region. REVISION 6 (2026-07-22) raised
//       kOffsetScaleFactor to 10, matching the reference plot's own convention
//       exactly -- see drawOneCorrectedSpectraPanel(). REVISION 7 (2026-07-22, same
//       day) then reverted this whole y-value-based approach: Andrew tried the
//       REVISION 6 render and said he preferred the ORIGINAL rank-index formula
//       (mult=10^b) for less error-bar overlap ("less noisy errorbars") -- see
//       drawOneCorrectedSpectraPanel() for the reinstated formula and the tradeoff
//       this reintroduces.
//   (3) replacing the verbose 2-line TLatex title/subtitle block with a compact
//       3-line corner annotation (title+centrality, system/energy, and one bullet
//       formula line folding in the offset formula and N_evt) -- matching the
//       reference's stacked-corner-annotation style directly.
// Untested -- please check a few rendered panels and confirm the new annotation is
// legible, markers read clearly as plain black, and the offset separation still
// looks reasonable (not too cramped, not too spread out) across the widest-rapidity
// panels.
//
// REVISION 4 (2026-07-21, same day): Andrew asked to bring color-by-rapidity back,
// using "the default coloring scheme" -- clarified as ROOT's own built-in default
// palette (kBird, a blue->green->yellow perceptually-uniform ramp), NOT a hand-rolled
// gradient table like the single-hue blue ramp REVISION 1-2 copied from
// PresentZFitterSpectra.C. Implemented by simply NOT calling gStyle->SetPalette(...)
// or TColor::CreateGradientColorTable(...) at all -- TColor::GetColorPalette(index)
// then samples whatever ROOT's active default palette already is, with no custom
// color-table setup needed. The rest of REVISION 3's changes are kept: the bounded
// y-value-based decade-offset formula (fixing the old 10^40 latent bug) and the
// compact 2-line corner annotation both stay as-is -- only the "plain black, no
// color" decision from REVISION 3 is reverted. The color-scale bar is back too
// (needed to make the restored color-by-|y| scheme interpretable, same right-side
// placement as PresentZFitterSpectra.C), and the corner annotation moved from
// top-right (REVISION 3, where the bar now sits again) to top-left so the two
// don't collide; the TPC/BTOF marker-style legend moved below the annotation, still
// on the left.
//
// REVISION 5 (2026-07-22): Andrew asked to "undo the symmetric rapidity coloring" --
// clarified this means REVISION 4's color-by-|y| scheme, where a +y bin and its
// mirror -y bin land on the exact same color (symmetric about y=0), not "remove
// color-by-rapidity again." Fixed by coloring by SIGNED y instead of |y|: colorFrac
// now maps linearly over the full [-a_kRapColorMax,+a_kRapColorMax] range instead of
// [0,a_kRapColorMax], so the active palette (still ROOT's default, unchanged from
// REVISION 4 -- no gStyle->SetPalette(...) call here either) is spread across the
// whole rapidity range and +y/-y bins get visibly different colors. drawColorScaleBar()
// was updated to match: its ticks/labels now span the signed range and its axis title
// changed from "|y|" to "y". Nothing else from REVISION 4 changed -- same color bar
// placement, same top-left annotation/legend positions, same offset formula.
//
// REVISION 6 (2026-07-22, same day): Andrew asked why the plot looked visually
// noisier/more background-overlapped in the lower-right (high mT-m0) region than
// before REVISION 3. Traced to REVISION 3's offset formula being tuned to a modest
// kOffsetScaleFactor=2 at the same time the rapidity-cut widening (a separate,
// same-session change) roughly doubled the number of visible rapidity-bin curves --
// packing far more curves into a narrower offset band than before REVISION 3. Fixed
// by raising kOffsetScaleFactor from 2 to 10, matching the reference plot's own
// "x3^{+-10y}" convention exactly.
//
// REVISION 7 (2026-07-22, same day): after seeing the REVISION 6 render, Andrew said
// he preferred the ORIGINAL rank-index decade offset (mult=10^b, b = the curve's
// position in the bins vector) over the REVISION 3-6 bounded y-value-based formula --
// less error-bar overlap between curves. Reverted the offset formula to that original
// unbounded form (see drawOneCorrectedSpectraPanel()); this reintroduces the latent
// property REVISION 3 had originally set out to fix (mult can reach ~10^40 for the
// widest-rapidity panels, since the rapidity-cut widening earlier this session means
// a panel can now have ~41 rapidity bins) -- but that large multiplier is exactly what
// gives the wide separation Andrew is asking for, and ROOT's log-y axis handles the
// resulting huge dynamic range without difficulty (it's a display multiplier applied
// per curve, not a physical quantity). Everything else (signed-y color mapping from
// REVISION 5, the compact corner annotation, plain-black-vs-colored choice from
// REVISION 3/4) is unchanged.

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TColor.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>

using namespace std;

//_______________________________________________________________________________
// RESTORED in REVISION 4 (2026-07-21) -- removed in REVISION 3, brought back now that
// color-by-|y| is back too. Copied near-verbatim from macros/PresentZFitterSpectra.C
// (see that macro's own header comment, ~line 583, for the full rationale behind every
// design choice here -- the TGaxis-didn't-render workaround, the barPad sizing, the
// top/bottom label clipping fix, etc.). a_nColors is passed in from the caller, which
// as of REVISION 4 is TColor::GetNumberOfColors() (ROOT's active default palette's
// color count) rather than a custom gradient table's fixed size -- the bar-drawing
// logic itself is otherwise unchanged.
//
// REVISION 5 (2026-07-22): Andrew asked to "undo the symmetric rapidity coloring" --
// clarified as: keep color-by-rapidity, but stop mapping +y and -y to the SAME color
// via |y| (REVISION 4's scheme). a_maxAbsY here is still the same fixed scale constant
// (kRapColorMax=1.3) but is now interpreted as the RANGE MAGNITUDE of a SIGNED y axis
// running from -a_maxAbsY (bottom of bar) to +a_maxAbsY (top of bar), matching the new
// signed colorFrac formula in drawOneCorrectedSpectraPanel() below -- so a +0.5 bin and
// a -0.5 bin now land on visibly different colors instead of the same one. The color-
// segment-drawing loop itself is unchanged (fracLo already ran bottom-to-top 0->1,
// which now just means "sampled color for signed y = -max + fracLo*2*max" instead of
// "sampled color for |y| = fracLo*max"); only the tick/label values and the axis title
// changed to reflect the new signed range.
void drawColorScaleBar(TPad* a_parentPad, double a_maxAbsY, int a_nColors, string a_namePrefix){
  if(a_maxAbsY <= 0) return;
  a_parentPad->cd();
  double barX1 = 0.86, barX2 = 0.99, barY1 = 0.12, barY2 = 0.88;
  TPad* barPad = new TPad((a_namePrefix+"_colorBarPad").c_str(),"",barX1,barY1,barX2,barY2);
  barPad->SetFillStyle(0);
  barPad->SetLeftMargin(0); barPad->SetRightMargin(0);
  barPad->SetTopMargin(0);  barPad->SetBottomMargin(0);
  barPad->Draw();
  barPad->cd();
  barPad->Range(0,0,1,1);
  double barFrac = 0.30;
  int nSegments = 100;
  for(int s = 0; s < nSegments; s++){
    double fracLo = (double)s/nSegments;
    double fracHi = (double)(s+1)/nSegments;
    int colorIndex = TColor::GetColorPalette((int)(fracLo*(a_nColors-1)));
    TBox* seg = new TBox(0,fracLo,barFrac,fracHi);
    seg->SetFillColor(colorIndex);
    seg->SetLineWidth(0);
    seg->Draw();
  }
  TBox* borderBox = new TBox(0,0,barFrac,1);
  borderBox->SetFillStyle(0);
  borderBox->SetLineColor(kBlack);
  borderBox->Draw();
  // REVISION 5: ticks/labels now span the SIGNED range [-a_maxAbsY,+a_maxAbsY] instead
  // of [0,a_maxAbsY] -- tick spacing widened from 0.1 to 0.2 (in y) since the covered
  // span doubled, keeping a similar total tick count to before.
  double tickStep = 0.2;
  int nTicks = TMath::Max(1,(int)TMath::Nint((2.0*a_maxAbsY)/tickStep));
  double tickLen  = 0.06;
  double labelGap = 0.02;
  for(int i = 0; i <= nTicks; i++){
    double yFrac = (double)i/nTicks;
    double v = -a_maxAbsY + yFrac*(2.0*a_maxAbsY);
    TLine* tick = new TLine(barFrac,yFrac,barFrac+tickLen,yFrac);
    tick->SetLineColor(kBlack);
    tick->Draw();
    TLatex* label = new TLatex(barFrac+tickLen+labelGap,yFrac,Form("%.1f",v));
    label->SetTextFont(43);
    label->SetTextSize(16);
    int vAlign = (i == 0) ? 1 : (i == nTicks) ? 3 : 2;
    label->SetTextAlign(10+vAlign);
    label->Draw();
  }
  TLatex* title = new TLatex(0.90,0.5,"y");
  title->SetTextFont(43);
  title->SetTextSize(20);
  title->SetTextAlign(22);
  title->SetTextAngle(90);
  title->Draw();
  a_parentPad->cd();
}

//_______________________________________________________________________________
// Copied verbatim from macros/PresentZFitterSpectra.C, ~line 702 -- see that macro's
// own header comment for where these percent ranges come from (macros/SetCutClass.C's
// cuts->setCentralities(6, ...) call).
string centralityLabel(int a_centIndex){
  static const char* kCentralityRanges[6] = {"0-5%","5-10%","10-20%","20-40%","40-80%","80-100%"};
  if(a_centIndex < 0 || a_centIndex >= 6) return Form("Cent%02d",a_centIndex);
  return kCentralityRanges[a_centIndex];
}

//_______________________________________________________________________________
// Copied verbatim from macros/PresentZFitterSpectra.C, ~line 717.
string particleSymbolWithCharge(string a_speciesName, string a_chargeLabel){
  bool isMinus = (!a_chargeLabel.empty() && a_chargeLabel[0] == 'M');
  if(a_speciesName == "Proton" && isMinus) return "#bar{p}";
  string symbol = "?";
  if(a_speciesName == "Pion")        symbol = "#pi";
  else if(a_speciesName == "Kaon")   symbol = "K";
  else if(a_speciesName == "Proton") symbol = "p";
  string chargeSign = isMinus ? "-" : "+";
  return symbol + "^{" + chargeSign + "}";
}

//_______________________________________________________________________________
// Recover the rapidity bin CENTER for a TGraphErrors written by
// RawSpectraModifier::makeRapiditySpectraTGraphs() -- see the file-format header
// comment above for why this is the only place the actual rapidity value survives
// into the output file. Returns -999 (and prints a warning) if the title doesn't
// contain a parseable "y=[lo,hi]" substring -- this should only happen if a future
// change to HistogramUtilities::VerticalProjectionGraphs()'s title format isn't
// reflected here too.
double parseRapidityCenterFromTitle(TGraphErrors* a_graph){
  if(!a_graph) return -999;
  double lo = 0, hi = 0;
  const char* title = a_graph->GetTitle();
  const char* marker = strstr(title, "y=[");
  if(!marker){
    cout << "WARNING: parseRapidityCenterFromTitle: could not find 'y=[' in title: " << title << endl;
    return -999;
  }
  if(sscanf(marker, "y=[%lf,%lf]", &lo, &hi) != 2){
    cout << "WARNING: parseRapidityCenterFromTitle: could not parse rapidity range out of title: " << title << endl;
    return -999;
  }
  return 0.5*(lo+hi);
}

//_______________________________________________________________________________
// STATISTICAL SIGNIFICANCE CUT (added 2026-07-23) -- see REVISION 8 in this file's
// header comment for the full rationale. Adapted from PresentZFitterSpectra.C's
// suppressInsignificantBins() for TGraphErrors rather than TH1D: removes any point
// whose significance (Y/EY) falls below a_minSignificance. Walks from the LAST point
// backward so RemovePoint() (which shifts every later index down by one) never
// invalidates an index this loop hasn't visited yet. A point with EY<=0 is treated as
// insignificant (significance forced to 0), matching suppressInsignificantBins()'s
// error<=0 convention -- not skipped/kept by default. a_minSignificance<=0 disables the
// cut entirely and every point survives, same convention as that macro.
int suppressInsignificantPoints(TGraphErrors* a_graph, double a_minSignificance){
  if(!a_graph || a_minSignificance <= 0) return 0;
  int nSuppressed = 0;
  for(int p = a_graph->GetN()-1; p >= 0; p--){
    double y = a_graph->GetY()[p];
    if(y == 0) continue; // already-empty point -- nothing to suppress
    double ey = a_graph->GetEY()[p];
    double significance = (ey > 0) ? (y/ey) : 0.0;
    if(significance < a_minSignificance){
      a_graph->RemovePoint(p);
      nSuppressed++;
    }
  }
  return nSuppressed;
}

//_______________________________________________________________________________
// FIT-QUALITY (chi^2/ndf) CUT (added 2026-07-23) -- see REVISION 9 in this file's header
// comment for the full rationale, including why this looks the value up by COORDINATE in
// a separately-opened ZFitter file rather than reading it out of a_correctedSpectraFile
// (RawSpectraModifier's output never carries chi^2/ndf at all). a_chiSqrHist is one of
// ZFitter's own per-bin chi^2/ndf TH2s (DeDx_FitData/.../ZTPC_ChiSqr_... or
// BTOF_FitData/.../ZbTOF_ChiSqr_...), X-axis=rapidity, Y-axis=mT-m0 -- same axes as the
// yield TH2 each point in a_graph was originally projected from (see
// PresentZFitterSpectra.C's own file-format comment). a_rapidityY is this graph's
// (already-parsed, single, shared) rapidity center -- passed in rather than re-derived
// per point since every point in a_graph comes from the same rapidity bin. Matches
// suppressBadFitQualityBins()'s badFit convention exactly: non-finite, <=0, or above
// a_maxChiSqrNdf all count as a bad fit. a_chiSqrHist==NULL (ZFitter file/histogram not
// found or not supplied) or a_maxChiSqrNdf<=0 disables the cut entirely, same convention
// as the other two cuts in this file.
int suppressBadFitQualityPoints(TGraphErrors* a_graph, double a_rapidityY, TH2* a_chiSqrHist, double a_maxChiSqrNdf){
  if(!a_graph || !a_chiSqrHist || a_maxChiSqrNdf <= 0) return 0;
  int rapBinIdx = a_chiSqrHist->GetXaxis()->FindBin(a_rapidityY);
  int nSuppressed = 0;
  for(int p = a_graph->GetN()-1; p >= 0; p--){
    double mtm0 = a_graph->GetX()[p];
    int mtm0BinIdx = a_chiSqrHist->GetYaxis()->FindBin(mtm0);
    double chiSqrNdf = a_chiSqrHist->GetBinContent(rapBinIdx, mtm0BinIdx);
    bool badFit = !TMath::Finite(chiSqrNdf) || chiSqrNdf <= 0 || chiSqrNdf > a_maxChiSqrNdf;
    if(badFit){
      a_graph->RemovePoint(p);
      nSuppressed++;
    }
  }
  return nSuppressed;
}

//_______________________________________________________________________________
// Andrew asked (2026-07-21) to double-check that the corrected spectra are actually
// divided by N_evt -- reasonable, since the offset-by-decade plot alone doesn't make
// that visible (each successive rapidity bin's curve is multiplied by a growing power
// of 10 purely for readability, which dominates the apparent y-axis scale and has
// nothing to do with whether N_evt division happened). This prints an explicit,
// numeric before/after comparison instead of asking anyone to eyeball a plot.
//
// Confirmed via direct read of RawSpectraModifier::convertSpectraToInvariant()
// (source/RawSpectraModifier.cxx ~line 2724-2766, called by RunRawSpectraModifier.C
// right after loadSpectraFile() and before makeRapiditySpectraTGraphs()) that EVERY
// bin is scaled by
//   1 / (dy * d(mT-m0) * 2*pi * mT * N_evt * efficiency)
// where efficiency is the a_efficiency argument RunRawSpectraModifier.C passes
// (overallEfficiency=1.0 in the command given to Andrew -- a placeholder for a future
// per-run efficiency scalar, NOT the real BTOF efficiency correction, which is a
// separate, later step -- see below). N_evt comes from m_centEvents[partIndex], loaded
// straight from the raw PicoBinner yield file's centEvents<Species> histogram.
//
// IMPORTANT CORRECTION (2026-07-21, after Andrew's first real rerun showed raw==
// corrected exactly, ratio(actual/predicted) in the thousands, for every single TPC
// point): this function's ORIGINAL comment here claimed the printed "ratio" should
// come out near 1.0 for a TPC point. That was wrong, and it's worth spelling out why,
// since it's a real gotcha in this file's structure, not just a units slip. Direct
// re-read of RawSpectraModifier::makeRapiditySpectraTGraphs() (~line 4360-4453) shows
// BOTH the "raw_..." and "corrected_..." TGraphErrors for a given point are built from
// projections of the SAME m_dEdxSpectraPlus TH2 -- and convertSpectraToInvariant()
// (the function this formula is quoting) already ran on that TH2 BEFORE either graph
// was created. So "raw_" here is NOT a truly-raw (uncorrected-for-N_evt) yield -- it's
// already fully N_evt/mT/dy/dmT-normalized, same as "corrected_". The raw/corrected
// split in THIS FILE is about a *later* step
// (loadAndApplyTPCEffAndEnergyLossAndBTOFEffFile(), which multiplies "corrected_"
// in place by TPC efficiency/energy-loss and/or BTOF efficiency, but leaves "raw_"
// untouched) -- not about N_evt normalization at all. A TPC point in a run with no TPC
// efficiency/energy-loss file (this run) legitimately gets NO further multiplication
// in that later step, so actualScale=correctedY/rawY==1.0 exactly is the CORRECT,
// expected result for it -- it says nothing whatsoever about whether N_evt division
// happened, because both graphs already have it baked in equally. The "ratio (actual/
// predicted)" number below is consequently not a pass/fail check either -- since
// actualScale is 1 by construction here, that ratio just reprints 1/predictedScale,
// i.e. the reciprocal of the assumed-still-raw-counts scale factor. It's kept below
// for transparency but should not be read as a normalization check.
//
// The actually meaningful numeric check, added below, is a BACK-CALCULATION: treat
// correctedY as "some true raw bin count, times predictedScaleNoEfficiency" and solve
// for that implied raw count. If N_evt normalization were genuinely happening, this
// should come out as a physically sensible raw yield count (order 1-few thousand for
// a high-statistics near-midrapidity TPC bin in this analysis) that trends sensibly
// with centrality (fewer particles per event in peripheral collisions despite more
// events). If N_evt were silently missing from the real formula, this back-calculation
// would instead imply an absurd sub-1 raw count. See printNormalizationCheck()'s body
// for the actual numbers from Andrew's first real rerun, which came out physically
// sensible (hundreds to ~2300, falling off toward peripheral) -- solid evidence FOR
// correct N_evt normalization even though the raw/corrected ratio itself couldn't
// show it directly for the reason above.
//
// This function ALSO now runs the same back-calculation-style check on a genuine BTOF
// point (see printBTOFCorrectionCheck() below) where a real, non-trivial raw/corrected
// difference IS expected (the ~44-50% BTOF efficiency correction from earlier in this
// analysis) -- that's a much stronger confirmatory test than the TPC point above,
// since it directly demonstrates the "corrected_" graph gets modified by a real,
// independently-known correction factor, rather than trivially equaling "raw_".
void printNormalizationCheck(TDirectory* a_dir, string a_speciesName, string a_chargeLabel,
                              int a_centIndex, int a_yIndex, double a_y, TGraphErrors* a_correctedTpc,
                              TH1D* a_centEventsHist, double a_speciesMassGeV){
  cout << "  [checked rapidity bin: yIndex " << a_yIndex << ", y~" << a_y << "]" << endl;
  if(!a_centEventsHist){
    cout << "NOTE: no centEvents" << a_speciesName << " (or bare centEvents) histogram found in this file "
         << "-- can't print a numeric N_evt normalization check for " << a_speciesName << a_chargeLabel
         << " Cent" << a_centIndex << "." << endl;
    return;
  }
  double nEvt = a_centEventsHist->GetBinContent(a_centIndex+1);
  cout << "N_evt check: " << a_speciesName << a_chargeLabel << " Cent" << a_centIndex << " -- N_evt = " << nEvt << endl;
  // Print explicitly WHY this bails, rather than returning silently -- the caller
  // picks whichever bin has a non-null TPC pointer, but "non-null" doesn't guarantee
  // the graph actually has points in it (VerticalProjectionGraphs() can still create
  // and write an empty TGraphErrors for a yIndex whose bin content was all zero).
  if(nEvt <= 0){
    cout << "NOTE: N_evt <= 0 for this centrality -- skipping the numeric scale check." << endl;
    return;
  }
  if(!a_correctedTpc || a_correctedTpc->GetN() == 0){
    cout << "NOTE: the TPC graph passed in for this check has zero points (yIndex " << a_yIndex
         << ") -- skipping. This can happen for an edge rapidity bin with a nearly-empty projection; "
         << "not a sign anything is wrong." << endl;
    return;
  }

  TGraphErrors* rawTpc = (TGraphErrors*) a_dir->Get(Form("raw_dEdxSpectra_%s_%s_Cent%02d_yIndex_%02d",
                            a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,a_yIndex));
  if(!rawTpc || rawTpc->GetN() == 0){
    cout << "NOTE: no raw_dEdxSpectra counterpart found for yIndex " << a_yIndex
         << " -- skipping the numeric before/after check for this panel (drawing continues normally)." << endl;
    return;
  }

  // Match the corrected graph's first point to its nearest-in-X raw counterpart --
  // cleanSpectra()/trimSpectra() can drop points between raw and corrected, so point
  // INDEX isn't guaranteed to line up even though both graphs came from the same
  // underlying TH2 bin grid.
  double correctedX = a_correctedTpc->GetX()[0];
  double correctedY = a_correctedTpc->GetY()[0];
  int bestRawPoint = -1;
  double bestDist = 1e300;
  for(int p = 0; p < rawTpc->GetN(); p++){
    double dist = TMath::Abs(rawTpc->GetX()[p] - correctedX);
    if(dist < bestDist){ bestDist = dist; bestRawPoint = p; }
  }
  if(bestRawPoint < 0) return;
  double rawX = rawTpc->GetX()[bestRawPoint];
  double rawY = rawTpc->GetY()[bestRawPoint];
  if(rawY == 0) return;

  double dy = 0, dmTm0 = 2.0*a_correctedTpc->GetEX()[0]; // GetEX() was set to half the mT-m0 bin width when this graph was first made (see file-format header comment)
  const char* title = a_correctedTpc->GetTitle();
  double lo=0, hi=0;
  if(sscanf(strstr(title,"y=["), "y=[%lf,%lf]", &lo, &hi) == 2) dy = hi - lo;
  double mT = correctedX + a_speciesMassGeV;

  double predictedScaleNoEfficiency = 1.0/(dy * dmTm0 * 2.0*TMath::Pi() * mT * nEvt);
  double actualScale = correctedY / rawY;
  // See this function's header comment for why actualScale==1 here is EXPECTED, not
  // a red flag: raw_/corrected_ in this file are both already N_evt-normalized, and
  // this TPC point receives no further correction in a run with no TPC eff/eloss file.
  // The genuinely informative number is the back-calculation below.
  double impliedRawCount = (predictedScaleNoEfficiency != 0) ? correctedY/predictedScaleNoEfficiency : 0;

  cout << "  TPC point check (mT-m0~" << rawX << " matched to " << correctedX << "): raw=" << rawY
       << " corrected=" << correctedY << " actual scale=" << actualScale
       << "  predicted scale (1/(dy*dmT*2pi*mT*N_evt), no efficiency/energy-loss term)=" << predictedScaleNoEfficiency
       << "  [note: actualScale==1 here is the EXPECTED result for a TPC point in this run -- NOT a normalization check; see header comment]" << endl;
  cout << "  Back-calculated implied raw bin count (correctedY / predictedScaleNoEfficiency) = " << impliedRawCount
       << " -- this IS the meaningful check: does this look like a plausible raw particle count for a "
       << "high-statistics near-midrapidity TPC bin (order 1-few thousand)? If N_evt normalization were "
       << "silently missing from the real formula, this would instead come out far below 1 (absurd)." << endl;
}

//_______________________________________________________________________________
// A STRONGER confirmatory test than printNormalizationCheck()'s TPC point (added
// 2026-07-21, same investigation): unlike a TPC point, a BTOF point's "corrected_..."
// graph DOES get modified relative to "raw_..." in this run, by the real, independently
// established BTOF efficiency (~44-50% this analysis, per earlier work this session --
// see the project log). That means correctedY/rawY for a BTOF point should come out
// close to 1/efficiency, i.e. roughly 2.0-2.3, NOT 1.0. Seeing a genuine, non-trivial
// ratio here (rather than the trivial exact-1.0 match the TPC point necessarily gives
// in this run) is much more convincing evidence that "corrected_" graphs really do get
// multiplicatively modified relative to "raw_" by a real physical correction factor --
// which in turn makes it far more credible that the earlier (invariant-normalization,
// N_evt-included) step is doing its job too, since both live in the same code path.
void printBTOFCorrectionCheck(TDirectory* a_dir, string a_speciesName, string a_chargeLabel,
                               int a_centIndex, int a_yIndex, double a_y, TGraphErrors* a_correctedBtof){
  cout << "  [checked rapidity bin: yIndex " << a_yIndex << ", y~" << a_y << "]" << endl;
  if(!a_correctedBtof || a_correctedBtof->GetN() == 0){
    cout << "NOTE: the BTOF graph passed in for the BTOF correction check has zero points (yIndex "
         << a_yIndex << ") -- skipping. Not a sign anything is wrong." << endl;
    return;
  }
  TGraphErrors* rawBtof = (TGraphErrors*) a_dir->Get(Form("raw_InvBetaBTOFSpectra_%s_%s_Cent%02d_yIndex_%02d",
                            a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,a_yIndex));
  if(!rawBtof || rawBtof->GetN() == 0){
    cout << "NOTE: no raw_InvBetaBTOFSpectra counterpart found for yIndex " << a_yIndex
         << " -- skipping the BTOF correction check for this panel (drawing continues normally)." << endl;
    return;
  }

  double correctedX = a_correctedBtof->GetX()[0];
  double correctedY = a_correctedBtof->GetY()[0];
  int bestRawPoint = -1;
  double bestDist = 1e300;
  for(int p = 0; p < rawBtof->GetN(); p++){
    double dist = TMath::Abs(rawBtof->GetX()[p] - correctedX);
    if(dist < bestDist){ bestDist = dist; bestRawPoint = p; }
  }
  if(bestRawPoint < 0) return;
  double rawX = rawBtof->GetX()[bestRawPoint];
  double rawY = rawBtof->GetY()[bestRawPoint];
  if(rawY == 0){
    cout << "NOTE: matched raw_InvBetaBTOFSpectra point (yIndex " << a_yIndex
         << ") has Y=0 -- can't form a ratio, skipping." << endl;
    return;
  }

  double actualScale = correctedY / rawY;
  cout << "  BTOF point check (mT-m0~" << rawX << " matched to " << correctedX << "): raw=" << rawY
       << " corrected=" << correctedY << " actual scale (corrected/raw)=" << actualScale
       << "  -- expect roughly 2.0-2.3 here (1/BTOF-efficiency, ~44-50% this analysis), and in any case "
       << "clearly NOT 1.0 -- if this comes out exactly 1.0 like the TPC check, that WOULD be a red flag "
       << "(it would mean the BTOF efficiency correction isn't actually being applied to this point)." << endl;
  if(actualScale == 1.0){
    cout << "  *** actual scale == 1.0 exactly -- this IS the red-flag case described above. If this persists "
         << "at a genuinely mid-rapidity bin (see the [checked rapidity bin: y~...] line just above -- this "
         << "selection now prefers the smallest |y| available, specifically to rule out an edge-rapidity-bin "
         << "explanation), grep the RunRawSpectraModifier.C run's own log for the line "
         << "\"WARNING: Not applying BTOF Eff for Part: ... y: <that y value>\" at this species/centrality to "
         << "confirm whether RawSpectraModifier itself skipped the BTOF correction for this exact bin (missing "
         << "efficiency object for this (Cent,yIndex) in efficiency_OO200_BTOF_v3.root) versus a deeper bug." << endl;
  }
}

//_______________________________________________________________________________
// One (species, charge, centrality) panel: gather every populated rapidity bin's
// corrected TPC/BTOF/ETOF TGraphErrors from a_dir, then draw them all overlaid, offset
// (optional) purely for visual separation between curves, colored by signed y on the
// fixed [-a_kRapColorMax,+a_kRapColorMax] scale using ROOT's active default palette
// (REVISION 5, 2026-07-22 -- see header comment). The offset multiplier is the
// ORIGINAL rank-index-based mult = 10^b, where b is the curve's position in the bins
// vector (REVISION 7, 2026-07-22 -- reverted from REVISION 3-6's bounded, y-value-based
// formula per Andrew's explicit preference for less error-bar overlap; see header
// comment for the full history, including the fact this can reach ~10^40 for the
// widest-rapidity panels now that the rapidity cut has been widened -- intentional,
// not a bug).
void drawOneCorrectedSpectraPanel(TDirectory* a_dir, string a_speciesName, string a_chargeLabel,
                                   int a_centIndex, double a_maxMtM0, bool a_offsetByDecade,
                                   string a_outDir, string a_systemLabel, int a_nColors,
                                   double a_kRapColorMax,
                                   TH1D* a_centEventsHist, double a_speciesMassGeV,
                                   double a_minSignificance,
                                   TH2* a_tpcChiSqrHist, TH2* a_btofChiSqrHist, TH2* a_etofChiSqrHist,
                                   double a_maxChiSqrNdf, double a_maxChiSqrNdfTPCMinus){

  // Offset formula history: REVISION 3 (2026-07-21) replaced the original rank-index
  // formula (mult=10^b) with a bounded, y-value-based one (mult=base^(scaleFactor*y))
  // to fix a latent unboundedness bug (b could reach ~40 after the rapidity cut was
  // widened, so mult could reach 10^40). REVISION 6 (2026-07-22) raised that formula's
  // scaleFactor from 2 to 10 to widen the curve separation, matching the reference
  // plot's own "x3^{+-10y}" convention. REVISION 7 (2026-07-22, same day): after
  // seeing the REVISION 6 render, Andrew said he preferred the ORIGINAL rank-index
  // formula (mult=10^b) over the bounded y-value-based one -- less error-bar overlap
  // between curves. Reverted to that original formula below (see the mult=... lines
  // in the pass-1/pass-2 loops over `bins`). This intentionally reintroduces the
  // unbounded-multiplier property REVISION 3 had fixed (mult can reach ~10^40 for the
  // widest-rapidity panels) -- accepted as the tradeoff for the wider separation
  // Andrew is asking for; ROOT's log-y axis handles the resulting large dynamic range
  // fine, since it's a per-curve display multiplier, not a physical quantity.

  struct RapBin { double y; int yIndex; TGraphErrors* tpc; TGraphErrors* btof; TGraphErrors* etof; };
  vector<RapBin> bins;

  // STATISTICAL SIGNIFICANCE CUT (REVISION 8, 2026-07-23) -- see suppressInsignificantPoints()
  // and this file's header comment for the rationale. Tallied per-panel and reported once
  // below, same convention as PresentZFitterSpectra.C's per-panel suppression summary.
  int nSignificanceSuppressed = 0;

  // FIT-QUALITY (chi^2/ndf) CUT (REVISION 9, 2026-07-23) -- see suppressBadFitQualityPoints()
  // and this file's header comment for the rationale, including why TPC gets its own
  // threshold when a_chargeLabel=="Minus" (antiprotons): that category's chi^2/ndf
  // distribution runs roughly 4x higher than every other detector/charge combination.
  int nChiSqrSuppressed = 0;
  double tpcMaxChiSqrNdf = (a_chargeLabel == "Minus") ? a_maxChiSqrNdfTPCMinus : a_maxChiSqrNdf;

  for(int yIndex = 0; yIndex < 120; yIndex++){
    TGraphErrors* tpc  = (TGraphErrors*) a_dir->Get(Form("corrected_dEdxSpectra_%s_%s_Cent%02d_yIndex_%02d",
                              a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,yIndex));
    TGraphErrors* btof = (TGraphErrors*) a_dir->Get(Form("corrected_InvBetaBTOFSpectra_%s_%s_Cent%02d_yIndex_%02d",
                              a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,yIndex));
    TGraphErrors* etof = (TGraphErrors*) a_dir->Get(Form("corrected_InvBetaETOFSpectra_%s_%s_Cent%02d_yIndex_%02d",
                              a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,yIndex));
    if(!tpc && !btof && !etof) continue;

    // Applied right after fetch, before this bin's rapidity is parsed or it's stored in
    // `bins` -- so every downstream consumer (the N_evt/BTOF numeric checks below and the
    // Pass 1/Pass 2 drawing loops further down) sees the same already-cut points, rather
    // than risking the checks and the plot disagreeing about which points survived.
    if(tpc)  nSignificanceSuppressed += suppressInsignificantPoints(tpc,  a_minSignificance);
    if(btof) nSignificanceSuppressed += suppressInsignificantPoints(btof, a_minSignificance);
    if(etof) nSignificanceSuppressed += suppressInsignificantPoints(etof, a_minSignificance);

    double y = -999;
    if(tpc)       y = parseRapidityCenterFromTitle(tpc);
    else if(btof) y = parseRapidityCenterFromTitle(btof);
    else          y = parseRapidityCenterFromTitle(etof);
    if(y == -999){
      cout << "WARNING: skipping " << a_speciesName << a_chargeLabel << " Cent" << a_centIndex
           << " yIndex " << yIndex << " -- could not determine its rapidity, see warning above." << endl;
      continue;
    }

    // Applied after the significance cut and after y is known (needed to look up the
    // right row of a_*ChiSqrHist) but still before this bin is stored -- same
    // already-cut-before-anyone-reads-it guarantee as the significance cut above.
    if(tpc)  nChiSqrSuppressed += suppressBadFitQualityPoints(tpc,  y, a_tpcChiSqrHist,  tpcMaxChiSqrNdf);
    if(btof) nChiSqrSuppressed += suppressBadFitQualityPoints(btof, y, a_btofChiSqrHist, a_maxChiSqrNdf);
    if(etof) nChiSqrSuppressed += suppressBadFitQualityPoints(etof, y, a_etofChiSqrHist, a_maxChiSqrNdf);

    RapBin b;
    b.y = y; b.yIndex = yIndex; b.tpc = tpc; b.btof = btof; b.etof = etof;
    bins.push_back(b);
  }

  if(bins.empty()){
    cout << "No corrected spectra found for " << a_speciesName << a_chargeLabel << " Cent" << a_centIndex
         << " -- skipping (expected for centralities beyond what this run actually populated)." << endl;
    return;
  }

  cout << "Drawing " << a_speciesName << a_chargeLabel << " Cent" << a_centIndex
       << ": " << bins.size() << " rapidity bin(s)" << endl;
  if(a_minSignificance > 0 && nSignificanceSuppressed > 0){
    cout << "  Suppressed (significance < " << a_minSignificance << " sigma): "
         << nSignificanceSuppressed << " point(s) across TPC/BTOF/ETOF, all rapidity bins, "
         << a_speciesName << a_chargeLabel << " Cent" << a_centIndex << endl;
  }
  if((a_maxChiSqrNdf > 0 || a_maxChiSqrNdfTPCMinus > 0) && nChiSqrSuppressed > 0){
    cout << "  Suppressed (chi^2/ndf > " << tpcMaxChiSqrNdf << " TPC, " << a_maxChiSqrNdf
         << " BTOF/ETOF): " << nChiSqrSuppressed << " point(s) across TPC/BTOF/ETOF, all rapidity bins, "
         << a_speciesName << a_chargeLabel << " Cent" << a_centIndex << endl;
  }

  // Explicit numeric N_evt/normalization check (Andrew, 2026-07-21) -- run once per
  // panel. IMPORTANT (revised same day): originally this picked the FIRST bin (in
  // yIndex order, i.e. lowest yIndex) with a non-empty point. That turned out to bias
  // toward an EDGE rapidity bin (RawSpectraModifier's own output covers a wider
  // rapidity range than the BTOF efficiency was actually fit across -- the
  // EfficiencyFitter diagnostic-image survey earlier this session only sampled
  // mid-rapidity y~0 and edge |y|~0.9, not necessarily this file's full range), which
  // is exactly the kind of bin least likely to have real efficiency-correction data
  // available for it. Now picks the bin with the SMALLEST |y| among those with a
  // non-empty point instead -- a mid-rapidity bin is far more likely to sit inside
  // whatever rapidity range the efficiency measurement actually covered, so this is a
  // fairer, more representative test than "whichever bin happened to be first."
  int bestTpcBin = -1; double bestTpcAbsY = 1e300;
  for(size_t b = 0; b < bins.size(); b++){
    if(bins[b].tpc && bins[b].tpc->GetN() > 0 && TMath::Abs(bins[b].y) < bestTpcAbsY){
      bestTpcBin = (int) b; bestTpcAbsY = TMath::Abs(bins[b].y);
    }
  }
  if(bestTpcBin >= 0){
    printNormalizationCheck(a_dir, a_speciesName, a_chargeLabel, a_centIndex, bins[bestTpcBin].yIndex,
                             bins[bestTpcBin].y, bins[bestTpcBin].tpc, a_centEventsHist, a_speciesMassGeV);
  }else{
    cout << "NOTE: no rapidity bin with a non-empty TPC point was found for " << a_speciesName << a_chargeLabel
         << " Cent" << a_centIndex << " -- skipping the numeric N_evt/scale check for this panel "
         << "(the plot itself is unaffected)." << endl;
  }

  // Same idea, but on a BTOF point (see printBTOFCorrectionCheck()'s header comment
  // for why this is the STRONGER confirmatory test -- a real, non-trivial correction
  // is expected here, unlike the TPC point above which trivially shows scale=1 in a
  // run with no TPC efficiency/energy-loss file). Same mid-rapidity-preferred
  // selection as the TPC check above, for the same reason.
  int bestBtofBin = -1; double bestBtofAbsY = 1e300;
  for(size_t b = 0; b < bins.size(); b++){
    if(bins[b].btof && bins[b].btof->GetN() > 0 && TMath::Abs(bins[b].y) < bestBtofAbsY){
      bestBtofBin = (int) b; bestBtofAbsY = TMath::Abs(bins[b].y);
    }
  }
  if(bestBtofBin >= 0){
    printBTOFCorrectionCheck(a_dir, a_speciesName, a_chargeLabel, a_centIndex, bins[bestBtofBin].yIndex,
                              bins[bestBtofBin].y, bins[bestBtofBin].btof);
  }else{
    cout << "NOTE: no rapidity bin with a non-empty BTOF point was found for " << a_speciesName << a_chargeLabel
         << " Cent" << a_centIndex << " -- skipping the BTOF correction check for this panel "
         << "(the plot itself is unaffected; expected if this bin's mT-m0 range never reaches the BTOF "
         << "momentum handoff, or if BTOF spectra weren't written for this species/charge/centrality)." << endl;
  }

  // Pass 1: figure out the y-axis range from the actual (post-offset) point values.
  double globalMinPositiveY = 1e300;
  double globalMaxY = 0;
  for(size_t b = 0; b < bins.size(); b++){
    double mult = a_offsetByDecade ? TMath::Power(10.0,(double)b) : 1.0;
    TGraphErrors* graphs[3] = {bins[b].tpc, bins[b].btof, bins[b].etof};
    for(int gi = 0; gi < 3; gi++){
      TGraphErrors* g = graphs[gi];
      if(!g) continue;
      for(int p = 0; p < g->GetN(); p++){
        double val = g->GetY()[p]*mult;
        if(val <= 0) continue;
        if(val < globalMinPositiveY) globalMinPositiveY = val;
        if(val > globalMaxY) globalMaxY = val;
      }
    }
  }

  if(globalMaxY <= 0){
    cout << "WARNING: every point for " << a_speciesName << a_chargeLabel << " Cent" << a_centIndex
         << " is zero/negative -- skipping plot." << endl;
    return;
  }

  double yFloor = TMath::Power(10.0, TMath::Floor(TMath::Log10(globalMinPositiveY)));

  TCanvas* canv = new TCanvas(Form("canv_%s%s_Cent%02d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex),"",900,700);
  canv->SetLogy();
  // REVISION 4 (2026-07-21): color-scale bar is back (see header comment), so restore
  // the right margin that makes room for it -- same convention as PresentZFitterSpectra.C.
  canv->SetRightMargin(0.16);
  // Extra headroom for the compact corner annotation block drawn below (replacing the
  // DrawFrame title string, which gStyle->SetOptTitle(0) silently suppresses -- see the
  // comment on the DrawFrame call itself). Must be set BEFORE DrawFrame() so the frame
  // is actually sized to leave this space, not just drawn under it.
  canv->SetTopMargin(0.14);

  // 2026-07-21: gStyle->SetOptTitle(0) (set in PresentRawSpectraModifierOutput() below)
  // suppresses ANY title string passed to DrawFrame -- it was silently never actually
  // rendered on any plot, including the "rapidity bins offset by decades for
  // readability" note that used to live inside it. Andrew asked for exactly that kind
  // of scaling information (anything not already captured by the y-axis label) to be
  // made visible. Pass an empty title here and draw the real labels explicitly via
  // TLatex further below instead, so they're guaranteed to actually render.
  TH1F* frame = canv->DrawFrame(0, yFloor, a_maxMtM0, globalMaxY*3, "");
  frame->GetXaxis()->SetTitle("m_{T}-m_{0} [GeV/c^{2}]");
  // 2026-07-21: fixed from a placeholder "Corrected Yield [arb. norm.]" after Andrew
  // asked about the normalization -- confirmed via direct read of
  // RawSpectraModifier::convertSpectraToInvariant() that every point here already has
  // 1/(2*pi*mT*N_evt*dy*d(mT-m0)) baked in (see printNormalizationCheck()'s header
  // comment for the full derivation + a runtime numeric check), PLUS the real BTOF
  // efficiency correction for BTOF points specifically. This is a proper superset of
  // PresentZFitterSpectra.C's own "1/(2pi N_evt) dN/dy d(mT-m0)" label -- that plot's
  // simpler label was always a stand-in for its earlier, pre-efficiency pipeline
  // stage (see that macro's own NORMALIZATION header comment) -- so the two axis
  // labels are intentionally NOT the same string; this one is the more complete,
  // standard invariant-yield form.
  frame->GetYaxis()->SetTitle("1/(2#pi m_{T} N_{evt}) d^{2}N/dy d(m_{T}-m_{0}) [(GeV/c^{2})^{-2}]");
  frame->GetYaxis()->CenterTitle(true);
  frame->GetYaxis()->SetTitleOffset(1.4);

  // Compact corner annotation block (kept from REVISION 3, per Andrew: "there's too
  // many words" -- title line, a system/energy line, and one short bullet line for
  // the offset formula + N_evt). REVISION 4 (2026-07-21): moved from top-right to
  // top-left, since the color-scale bar is back and occupies the top-right corner
  // (the right-margin strip) again.
  double annX = 0.13;
  TLatex* titleLine = new TLatex(annX,0.965,
    Form("%s Spectra, %s, %s",
         particleSymbolWithCharge(a_speciesName,a_chargeLabel).c_str(),
         a_systemLabel.c_str(),
         centralityLabel(a_centIndex).c_str()));
  titleLine->SetNDC();
  titleLine->SetTextFont(43);
  titleLine->SetTextSize(20);
  titleLine->SetTextAlign(13); // left, top
  titleLine->Draw();

  //double nEvtForLabel = (a_centEventsHist) ? a_centEventsHist->GetBinContent(a_centIndex+1) : -1;
  //if(a_offsetByDecade){
  //  string bulletLabel = "#bullet #times10^{b} offset (b = curve rank)";
  //  if(nEvtForLabel > 0) bulletLabel += Form(", N_{evt}=%.0f", nEvtForLabel);
  //  TLatex* bulletLine = new TLatex(annX,0.925,bulletLabel.c_str());
  //  bulletLine->SetNDC();
  //  bulletLine->SetTextFont(43);
  //  bulletLine->SetTextSize(15);
  //  bulletLine->SetTextAlign(13); // left, top
  //  bulletLine->Draw();
  //}else if(nEvtForLabel > 0){
  //  TLatex* nEvtLine = new TLatex(annX,0.925,Form("N_{evt} = %.0f",nEvtForLabel));
  //  nEvtLine->SetNDC();
  //  nEvtLine->SetTextFont(43);
  //  nEvtLine->SetTextSize(15);
  //  nEvtLine->SetTextAlign(13); // left, top
  //  nEvtLine->Draw();
  //}*/

  // Pass 2: draw every rapidity bin's TPC/BTOF/ETOF points, offset per the y-value-based
  // formula above, colored by SIGNED y (REVISION 5, 2026-07-22 -- see below) using
  // ROOT's active default palette. Cloned before scaling so the offset multiplier never
  // mutates the actual objects loaded live out of the file.
  //
  // REVISION 5 (2026-07-22): Andrew asked to "undo the symmetric rapidity coloring" --
  // REVISION 4's colorFrac used |y|, so a +0.5 and a -0.5 bin mapped to the exact same
  // color (symmetric about y=0), which is exactly the "symmetric coloring" complaint.
  // Fixed by mapping colorFrac from SIGNED y over the full [-a_kRapColorMax,
  // +a_kRapColorMax] range instead of |y| over [0,a_kRapColorMax] -- so the full active
  // palette (e.g. kBird's blue->green->yellow) is now spread across the whole rapidity
  // range, with -y and +y bins landing on visibly different colors. drawColorScaleBar()
  // above was updated to match (ticks/labels/title now reflect the signed range).
  for(size_t b = 0; b < bins.size(); b++){
    double mult = a_offsetByDecade ? TMath::Power(10.0,(double)b) : 1.0;
    double colorFrac = TMath::Max(0.0, TMath::Min(1.0, (bins[b].y + a_kRapColorMax) / (2.0*a_kRapColorMax)));
    int colorIndex = TColor::GetColorPalette((int)(colorFrac*(a_nColors-1)));

    if(bins[b].tpc){
      TGraphErrors* g = (TGraphErrors*) bins[b].tpc->Clone();
      for(int p = 0; p < g->GetN(); p++) g->SetPoint(p, g->GetX()[p], g->GetY()[p]*mult);
      g->SetMarkerStyle(20); // filled circle = TPC, same convention as PresentZFitterSpectra.C
      g->SetMarkerColor(colorIndex);
      g->SetLineColor(colorIndex);
      g->SetMarkerSize(0.7);
      g->Draw("PE SAME");
    }
    if(bins[b].btof){
      TGraphErrors* g = (TGraphErrors*) bins[b].btof->Clone();
      for(int p = 0; p < g->GetN(); p++) g->SetPoint(p, g->GetX()[p], g->GetY()[p]*mult);
      g->SetMarkerStyle(24); // open circle = BTOF
      g->SetMarkerColor(colorIndex);
      g->SetLineColor(colorIndex);
      g->SetMarkerSize(0.7);
      g->Draw("PE SAME");
    }
    if(bins[b].etof){
      TGraphErrors* g = (TGraphErrors*) bins[b].etof->Clone();
      for(int p = 0; p < g->GetN(); p++) g->SetPoint(p, g->GetX()[p], g->GetY()[p]*mult);
      g->SetMarkerStyle(25); // open square = ETOF
      g->SetMarkerColor(colorIndex);
      g->SetLineColor(colorIndex);
      g->SetMarkerSize(0.7);
      g->Draw("PE SAME");
    }
  }

  // REVISION 4 (2026-07-21): drawColorScaleBar() call restored, same right-side
  // placement as PresentZFitterSpectra.C.
  //
  // REVISION 9 (2026-07-22, same day as REVISION 8): per Andrew's explicit request,
  // moved the TPC/BTOF marker-style legend from the top-left (below the corner
  // annotation) to the upper-right portion of the plot instead. The color-scale bar
  // itself only occupies x in [0.86,0.99] (see drawColorScaleBar()'s barX1/barX2), so
  // there's genuinely open space between the corner annotation (x in [0.13,0.5], y in
  // [0.89,0.975]) and the bar -- the legend now sits at x in [0.55,0.85], y in
  // [0.78,0.885], just left of the bar and below the annotation's row, so it collides
  // with neither.
  drawColorScaleBar((TPad*)gPad, a_kRapColorMax, a_nColors,
                     Form("bar_%s%s_Cent%02d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex));

  TLegend* leg = new TLegend(0.55,0.74,0.72,0.855);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  TGraphErrors* dummyTPC  = new TGraphErrors(); dummyTPC->SetMarkerStyle(20);  dummyTPC->SetMarkerColor(kBlack);
  TGraphErrors* dummyBTOF = new TGraphErrors(); dummyBTOF->SetMarkerStyle(24); dummyBTOF->SetMarkerColor(kBlack);
  leg->AddEntry(dummyTPC,Form("#times10^{b} offset (b = curve rank)"),"p");
  leg->AddEntry(dummyTPC,"TPC","p");
  leg->AddEntry(dummyBTOF,"BTOF","p");
  leg->Draw();

  string outName = Form("%s/CorrectedSpectra_%s%s_Cent%02d.png",
                         a_outDir.c_str(), a_speciesName.c_str(), a_chargeLabel.c_str(), a_centIndex);
  canv->SaveAs(outName.c_str());
  delete canv;
}

//_______________________________________________________________________________
void PresentRawSpectraModifierOutput(string a_correctedSpectraFile,
                                      int a_partIndex = 0, // 0=Pion, 1=Kaon, 2=Proton -- matches loadDataFile's convention throughout this repo
                                      string a_outDir = "./corrected_spectra_plots",
                                      bool a_offsetByDecade = true,
                                      double a_minSignificance = 5.0, // content/error cut on drawn points -- see REVISION 8 header comment; pass 0 (or negative) to disable
                                      string a_tpcZFitterFile = "", // ZFitter TPC output file (e.g. RunZFitter.C's TPC spectra_OO200_<Species>_TPC_RapBin_...root) -- see REVISION 9 header comment. "" = chi^2/ndf cut off
                                      string a_btofZFitterFile = "", // same, BTOF output file (also supplies ETOF's chi^2 hist, per PresentZFitterSpectra.C's naming quirk)
                                      double a_maxChiSqrNdf = 0, // chi^2/ndf ceiling for TPC Plus/BTOF Plus/BTOF Minus/ETOF -- pass 0 (or negative) to disable
                                      double a_maxChiSqrNdfTPCMinus = 0){ // separate, looser TPC ceiling for antiprotons -- see REVISION 9 header comment for why; 0 (or negative) disables

  gSystem->mkdir(a_outDir.c_str(), true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  static const char* kSpeciesNames[3] = {"Pion","Kaon","Proton"};
  if(a_partIndex < 0 || a_partIndex > 2){
    cout << "ERROR: a_partIndex must be 0 (Pion), 1 (Kaon), or 2 (Proton)." << endl;
    return;
  }
  string speciesName = kSpeciesNames[a_partIndex];

  // Per-species mT-m0 axis ceiling -- same values/reasoning as
  // PresentZFitterSpectra.C's speciesMaxMtM0GeV[] (chosen so BTOF/ETOF curves are
  // never clipped; see that macro's own comment for the exact numbers' provenance).
  const double kSpeciesMaxMtM0GeV[3] = {4.0, 4.0, 5.0};
  double maxMtM0 = kSpeciesMaxMtM0GeV[a_partIndex];

  // Standard PDG charged-particle masses (GeV/c^2) -- used only for the mT = (mT-m0)+m0
  // conversion in printNormalizationCheck()'s numeric check below. Deliberately
  // hardcoded here rather than instantiating ParticleInfo (which RawSpectraModifier
  // itself uses internally) so this macro stays free of any ../bin/ compiled-library
  // dependency -- these are the standard PDG values ParticleInfo::GetParticleMass()
  // is expected to match.
  const double kSpeciesMassGeV[3] = {0.13957, 0.49368, 0.93827}; // Pion+-, Kaon+-, Proton/antiproton
  double speciesMassGeV = kSpeciesMassGeV[a_partIndex];

  TFile* inFile = TFile::Open(a_correctedSpectraFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "ERROR: could not open " << a_correctedSpectraFile << endl;
    return;
  }

  // FIT-QUALITY (chi^2/ndf) CUT (REVISION 9, 2026-07-23) -- see this file's header
  // comment for the full rationale. These are opened once here (not per panel) and the
  // relevant TH2 looked up fresh per (species,charge,centIndex) below, same pattern
  // PresentZFitterSpectra.C uses for its own chi^2 histogram reads. Left NULL (cut off
  // for that detector) if the corresponding filename argument is empty or fails to open --
  // this macro still runs and draws everything else normally either way.
  TFile* tpcZFitterFile = (a_tpcZFitterFile != "") ? TFile::Open(a_tpcZFitterFile.c_str(),"READ") : NULL;
  if(a_tpcZFitterFile != "" && (!tpcZFitterFile || tpcZFitterFile->IsZombie())){
    cout << "WARNING: could not open a_tpcZFitterFile (" << a_tpcZFitterFile
         << ") -- the TPC chi^2/ndf cut will be skipped for every panel." << endl;
    tpcZFitterFile = NULL;
  }
  TFile* btofZFitterFile = (a_btofZFitterFile != "") ? TFile::Open(a_btofZFitterFile.c_str(),"READ") : NULL;
  if(a_btofZFitterFile != "" && (!btofZFitterFile || btofZFitterFile->IsZombie())){
    cout << "WARNING: could not open a_btofZFitterFile (" << a_btofZFitterFile
         << ") -- the BTOF/ETOF chi^2/ndf cut will be skipped for every panel." << endl;
    btofZFitterFile = NULL;
  }

  // N_evt lives in centEvents<Species> (per-species PicoBinner yield file, loaded by
  // RawSpectraModifier::loadDataFile() and re-written into this output file by
  // writeMemberData() -- confirmed via direct read, RawSpectraModifier.cxx ~line
  // 3649), falling back to the bare "centEvents" name used when every species was
  // binned into one shared file (same two-name fallback convention already used by
  // PresentZFitterSpectra.C's getNumEvents()).
  TH1D* centEventsHist = (TH1D*) inFile->Get(Form("centEvents%s", speciesName.c_str()));
  if(!centEventsHist) centEventsHist = (TH1D*) inFile->Get("centEvents");
  if(!centEventsHist){
    cout << "WARNING: no centEvents" << speciesName << " (or bare centEvents) histogram found in "
         << a_correctedSpectraFile << " -- the numeric N_evt normalization check below will be skipped "
         << "for every panel (the plots themselves are unaffected -- N_evt normalization already happened "
         << "inside RunRawSpectraModifier.C regardless of whether this histogram survived into this file)." << endl;
  }

  // REVISION 4 (2026-07-21): color-by-|y| restored, per Andrew's request to use
  // "the default coloring scheme" -- clarified as ROOT's own built-in default palette
  // (kBird), not the single-hue blue gradient table REVISION 1-2 had copied from
  // PresentZFitterSpectra.C. Originally deliberately did NOT call gStyle->SetPalette(...)
  // or TColor::CreateGradientColorTable(...) at all -- TColor::GetColorPalette(idx)
  // below simply samples whatever ROOT's already-active default palette is.
  //
  // REVISION 8 (2026-07-22, same day as REVISION 7): Andrew asked for the color scale
  // to reach red. Asked via AskUserQuestion which approach he wanted -- a dedicated
  // diverging blue-gray-red ramp (drafted earlier this session via the dataviz skill,
  // matched to signed-y/polarity data, but never implemented) vs. ROOT's built-in
  // kRainBow palette (reaches red, but is the same non-perceptually-monotonic rainbow
  // palette already replaced once this session in PresentZFitterSpectra.C, since nearby
  // values there could look similar). Andrew explicitly chose kRainBow. Implemented by
  // adding a single `gStyle->SetPalette(kRainBow);` call right here, so
  // TColor::GetColorPalette(idx) below now samples kRainBow (blue->cyan->green->yellow
  // ->red) instead of the previously-active kBird (blue->green->yellow, no red) --
  // every other REVISION 4/5 choice (signed-y colorFrac mapping, no custom gradient
  // table beyond this one explicit SetPalette call, drawColorScaleBar()'s tick/label/
  // axis-title logic) is unchanged. This is a deliberate step back from the dataviz
  // skill's "never a rainbow" guidance, made knowingly per Andrew's explicit choice --
  // not an oversight.
  gStyle->SetPalette(kRainBow);
  int nColors = TColor::GetNumberOfColors();
  const double kRapColorMax = 1.3; // same fixed scale as PresentZFitterSpectra.C

  // Compact system/energy label for the corner annotation drawn in
  // drawOneCorrectedSpectraPanel() -- hardcoded to this analysis's actual collision
  // system, matching the reference image's own "AuAu"/"sqrt(s_NN)=3.2 GeV" convention.
  string systemLabel = "O+O, #sqrt{s_{NN}} = 200 GeV";

  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    string chargeLabel = (pmIndex == 0) ? "Plus" : "Minus";
    // Directory naming uses ParticleInfo::GetParticleName(partIndex,+-1)'s convention:
    // species name and Plus/Minus CONCATENATED with no underscore (see file-format
    // header comment above) -- e.g. "Spectra_PionPlus", not "Spectra_Pion_Plus".
    string dirName = Form("Spectra_%s%s", speciesName.c_str(), chargeLabel.c_str());
    TDirectory* dir = (TDirectory*) inFile->Get(dirName.c_str());
    if(!dir){
      cout << "WARNING: directory " << dirName << " not found in " << a_correctedSpectraFile << " -- skipping this charge sign entirely." << endl;
      continue;
    }

    for(int centIndex = 0; centIndex < 16; centIndex++){
      // Object naming/paths match PresentZFitterSpectra.C's own chi^2 histogram reads
      // exactly (see that macro, ~line 1235-1240) -- including the ETOF quirk where the
      // object still lives under "ETOF_FitData/.../ZbTOF_ChiSqr_..." (not "ZeTOF_..."),
      // sourced from the SAME btofZFitterFile as the BTOF histogram, not a separate file.
      TH2* tpcChiSqrHist  = tpcZFitterFile  ? (TH2*) tpcZFitterFile->Get(Form("DeDx_FitData/%s/ZTPC_ChiSqr_%s%s_Cent%02d",
                                speciesName.c_str(),speciesName.c_str(),chargeLabel.c_str(),centIndex)) : NULL;
      TH2* btofChiSqrHist = btofZFitterFile ? (TH2*) btofZFitterFile->Get(Form("BTOF_FitData/%s/ZbTOF_ChiSqr_%s%s_Cent%02d",
                                speciesName.c_str(),speciesName.c_str(),chargeLabel.c_str(),centIndex)) : NULL;
      TH2* etofChiSqrHist = btofZFitterFile ? (TH2*) btofZFitterFile->Get(Form("ETOF_FitData/%s/ZbTOF_ChiSqr_%s%s_Cent%02d",
                                speciesName.c_str(),speciesName.c_str(),chargeLabel.c_str(),centIndex)) : NULL;

      drawOneCorrectedSpectraPanel(dir, speciesName, chargeLabel, centIndex, maxMtM0,
                                    a_offsetByDecade, a_outDir, systemLabel, nColors,
                                    kRapColorMax, centEventsHist, speciesMassGeV,
                                    a_minSignificance,
                                    tpcChiSqrHist, btofChiSqrHist, etofChiSqrHist,
                                    a_maxChiSqrNdf, a_maxChiSqrNdfTPCMinus);
    }
  }

  cout << "Done. Plots written to " << a_outDir << endl;
}
