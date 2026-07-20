// PresentZFitterSpectra.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: the original request (before the pion TPC crash/Bichsel-curve
// detour took over) was a macro to plot the pion spectra straight out of ZFitter's
// output -- raw yield vs mT-m0, overlaying all rapidity bins, combining TPC and BTOF
// on the same plot. This is that macro, generalized to any of the three species
// (pion/kaon/proton) this analysis fits, since the plotting logic doesn't depend on
// species beyond the name string.
//
// WHAT IT READS (names/paths confirmed 2026-07-16 directly against a real ZFitter
// output file's TFile::ls(), not just read from source -- see the naming caveat below
// for why that direct confirmation mattered):
//   TPC:  DeDx_FitData/<Species>/dEdxSpectra_<Species>Plus_Cent%02d
//         DeDx_FitData/<Species>/dEdxSpectra_<Species>Minus_Cent%02d
//   BTOF: BTOF_FitData/<Species>/InvBetaBTOFSpectra_<Species>Plus_Cent%02d
//         BTOF_FitData/<Species>/InvBetaBTOFSpectra_<Species>Minus_Cent%02d
//   ETOF: ETOF_FitData/<Species>/InvBetaETOFSpectra_<Species>Plus_Cent%02d
//         ETOF_FitData/<Species>/InvBetaETOFSpectra_<Species>Minus_Cent%02d
//         (only drawn if present -- ETOF booking in ZFitter.cxx is gated on
//         a_detectorIndex==2, which RunZFitter.C never actually passes -- see its own
//         header comment, "ETOF fit drivers were not ported" -- so in practice these
//         are expected to stay absent regardless of which file is passed)
// All three are TH2D with X-axis = rapidity bin (title "y_{...}"), Y-axis = mT-m0
// [GeV/c^2] (title "m_{T}-m_{...}"), bin content = fitted yield/amplitude.
// "Overlaying all rapidity bins" means, for a fixed species/charge/centrality: take a
// Y-axis projection (ProjectionY) at each non-empty X bin (rapidity bin) to get a 1D
// yield-vs-(mT-m0) spectrum, then draw all of those on one pad. As of 2026-07-16
// (Round 9), each +y bin is further combined (averaged) with its -y mirror bin before
// drawing -- see COMBINING +-y RAPIDITY BINS below -- so what's actually drawn is one
// curve per |y| bin, not one curve per raw rapidity bin.
//
// NAMING CAVEAT (2026-07-16, corrected from an earlier wrong "fix"): there is NO
// underscore between species and Plus/Minus. ZFitter.cxx's makeFitDataHistos() has TWO
// SetName() calls touching the same m_Spectra_Cent_ZTPC element in sequence -- an
// earlier one (~line 2513) using Form("dEdxSpectra_%s_Plus_Cent%02d", bare species
// name, ...) which DOES have an underscore, and a later one in the same function
// (~line 2269) using Form("dEdxSpectra_%s_Cent%02d", GetParticleName(idx,chargePM),
// ...) where GetParticleName with a nonzero charge already returns "PionPlus" as one
// string -- no separate underscore. The LATER call is what's actually on disk
// (confirmed via TFile::ls() against a real output file: "dEdxSpectra_PionPlus_Cent00").
// An earlier pass at this macro (and at RawSpectraModifier.cxx, see that file's own
// history) added an underscore based on the EARLIER, non-final SetName() call and got
// it backwards. If you ever see NULL histograms reading a ZFitter output file with a
// name that looks *almost* right, check TFile::ls() against the real file directly --
// don't trust a single grep hit in a large function that may set the same object's
// name more than once.
//
// LEGEND CHOICE: a fit typically has dozens of rapidity bins (the pion TPC crash log
// referenced RapBin 35, i.e. at least 36 bins) -- one legend entry per bin would be
// unreadable, and even after the Round 9 +y/-y combination below roughly halves that
// count, ~21 curves is still too many for a clean legend box. Instead, each combined
// |y|-bin curve is colored along a continuous rainbow gradient (dark = y~0/center,
// bright = largest |y|). Rounds 9-10 tried inline "y=+-<bincenter>" text labels next to
// each curve's last point (every other curve, to reduce clutter), but Round 11 replaced
// that with a proper COLZ-style color-scale bar in the pad's right margin (see
// drawColorScaleBar() and its own Round 11 header comment) per Andrew's request --
// a continuous bar with a real tick-labeled axis reads |y| for any curve directly from
// its color, without needing text scattered across the plot at all. Detector is still
// distinguished by marker style (filled circle = TPC, open circle = BTOF, open square =
// ETOF) via a small 3-entry legend.
//
// TWO SEPARATE INPUT FILES, NOT ONE (updated 2026-07-16): RunZFitter.C's
// a_detectorIndex (0=TPC, 1=BTOF) makes each ZFitter invocation fit only ONE detector
// -- confirmed in RunZFitter.C, where a_detectorIndex directly sets fit<Species>TPC vs
// fit<Species>BTOF. So a completed TPC run's output file has a populated
// DeDx_FitData/ folder but an empty/absent BTOF_FitData/ (no BTOF fit ever ran in that
// process), and vice versa for a BTOF run's output file. To actually combine TPC+BTOF
// on one plot, this macro reads them from the two separate output files ZFitter
// produced. ETOF is read from the BTOF file too, since ETOF_FitData (when non-empty)
// is written by the same makeSpectra() call as BTOF_FitData.
//
// Usage:
//   root -l -q 'PresentZFitterSpectra.C("/path/to/tpc_output.root","/path/to/btof_output.root",0,0,".","/path/to/yieldHistos_OO200_pion.root")'
//     -> pion (partIndex 0), centrality bin 0, both charges, output to ".", N_evt read
//        from the raw PicoBinner yield file (required for N_evt normalization -- see below)
//   root -l -q 'PresentZFitterSpectra.C("/path/to/tpc_output.root","",1,2,"/path/to/outDir","/path/to/yieldHistos_OO200_kaon.root")'
//     -> kaon (partIndex 1), centrality bin 2, TPC only (no BTOF file yet), output to /path/to/outDir
//
// a_tpcZFitterFile:  path to the TPC-detector ZFitter output file (a_detectorIndex=0 run).
//                    Pass "" if no TPC run exists yet -- BTOF-only will still draw.
// a_btofZFitterFile: path to the BTOF-detector ZFitter output file (a_detectorIndex=1 run).
//                    Pass "" (or omit) if no BTOF run exists yet -- draws TPC only.
// a_partIndex: 0=Pion, 1=Kaon, 2=Proton (matches ParticleInfo's speciesIndex ordering)
// a_centIndex: centrality bin, matches ZFitter's nCentToRun indexing (0 = most central
//              in this analysis's 6-bin scheme, per README's centrality section)
// a_yieldFile: path to the RAW PicoBinner yield file for this species (e.g. what
//              RunZFitter.C was given as a_pionYieldFile/a_kaonYieldFile/a_protonYieldFile)
//              -- NOT a ZFitter output file. This is where N_evt actually lives (see
//              NORMALIZATION below); pass "" to skip N_evt normalization entirely.
// a_minSignificance: minimum content/error ("sigma") a bin must have to be drawn at all
//              -- see STATISTICAL SIGNIFICANCE CUT below. Default 5.0 as of the tail-
//              cleanup follow-ups (was 2.0 since Round 12, loosened from an original 3.0;
//              then tightened back to 4.0, then 5.0, after Andrew found the high-mT-m0
//              tail still visibly scattered); pass 0 (or negative) to disable the cut
//              entirely and draw every nonzero bin again, e.g. for debugging what's
//              actually in the raw fit output.
// a_maxChiSqrNdf: maximum per-bin chi^2/ndf (from ZFitter's own ZTPC_ChiSqr_.../
//              ZbTOF_ChiSqr_... histograms) a bin's underlying fit may have and still be
//              drawn -- see suppressBadFitQualityBins() below. Independent of, and in
//              addition to, the significance cut: a bad fit can converge to a
//              "confidently wrong" value with a small quoted error, which no amount of
//              tightening the significance threshold alone can catch (this is what
//              motivated adding this cut -- two isolated outlier points survived even a
//              5-sigma significance cut). DEFAULT 0 (disabled) for now -- an initial
//              3.0 guess turned out to suppress EVERY bin (reported by Andrew as "all
//              rapidity bins empty"), meaning either this analysis's real chi^2/ndf
//              distribution runs well above the textbook "good fit" range of ~1-3, or
//              there's an axis-alignment bug in how suppressBadFitQualityBins() reads
//              the chi^2 histogram -- CheckChiSqrValues.C (new diagnostic macro) prints
//              real (yield, chi^2/ndf) pairs so the right threshold (and correctness of
//              the axis alignment) can be confirmed empirically instead of guessed
//              again. Pass a positive value once that's known to enable the cut.
//
// NORMALIZATION (added 2026-07-16, per Andrew's request to normalize/sum-over-2pi/
// normalize-by-rapidity-bin/rescale axes): each rapidity-bin projection is converted
// from a raw fit yield to
//   1/(2*pi*N_evt*dy) * dN/d(mT-m0)   [per event, per unit rapidity]
// i.e. divided by the number of events in this centrality bin, by 2*pi (since the
// TPC/BTOF acceptance is already summed over the full 2*pi azimuth -- this is the
// standard convention for writing an invariant-ish differential yield), by dy (that
// rapidity bin's width, read off the source TH2's X-axis, so bins are compared per
// unit rapidity rather than per whatever-width-this-bin-happens-to-be), AND by the
// mT-m0 bin width (ROOT's Scale(1.0,"width")) so bins of different widths are
// comparable. This mirrors two pieces of precedent already in this codebase:
//   - ZFitter::makeSingleFitHisto() (~line 1732) already divides an individual
//     rapidity-bin/mT-m0-bin projection by N_evt (with the same relative-error
//     quadrature-sum convention used here) and by bin width via Scale(1.0,"width")
//     (but NOT by dy, since it's just extracting one fixed rapidity bin, not
//     presenting multiple rapidity bins side by side the way this macro does).
//   - RawSpectraModifier::convertSpectraToInvariant() (~line 2686) does the fuller
//     invariant-yield conversion later in the pipeline: divides by rapidity_width *
//     mTm0Width * 2*pi * mT * N_evt * efficiency -- this macro now matches that
//     formula except for the mT and efficiency factors (see below).
// This macro deliberately stops short of that fuller RawSpectraModifier conversion --
// no division by mT, no efficiency correction -- since efficiency depends on an
// embedding-derived file this macro doesn't have and isn't needed for a quick
// per-detector sanity-check plot straight out of ZFitter. If N_evt can't be found
// (missing "centEvents<Species>"/"centEvents" histogram in the given a_yieldFile, or no
// a_yieldFile given at all), the projections are left un-normalized by N_evt (2*pi, dy,
// and bin-width scaling still applied) and a warning is printed once.
//
// WHERE N_EVT ACTUALLY LIVES (found 2026-07-16, after the first version of this
// normalization "looked multiplied instead of divided"): centEvents is NEVER written
// into a ZFitter output file -- traced ZFitter::makeSpectra() (the function that
// actually writes DeDx_FitData/BTOF_FitData/ETOF_FitData) end to end and confirmed it.
// centEvents<Species> is written by PicoBinner.cxx into the RAW YIELD FILE (e.g.
// yieldHistos_OO200_pion.root); ZFitter only ever reads it back out of that file for
// its own internal use (loadDataHistograms(), ~line 904) and never re-writes it. So
// N_evt has to be read from a_yieldFile below, not from a_tpcZFitterFile/
// a_btofZFitterFile -- looking for it in the ZFitter output files (the first version of
// this macro's mistake) will always silently come back empty.
//
// STATISTICAL SIGNIFICANCE CUT (added 2026-07-16, Round 8, after Andrew looked at the
// actual output PNG and asked how to clean up the raw spectra): past roughly mT-m0 ~0.5
// GeV/c^2 the plotted points stopped looking like real physics and turned into a scatter
// of huge vertical error bars spanning many orders of magnitude down to ~1e-17/1e-18 --
// i.e. statistically empty or near-empty fit bins whose tiny positive fluctuation still
// survives normalization and gets plotted as if it were a real data point. Two things
// followed from that: (1) those bins are visually meaningless and should not be drawn at
// all, and (2) they were exactly what was setting the Round 7 auto y-floor down near
// 1e-17/1e-18, squeezing the ~3-4 real decades of physics into a sliver at the top of an
// 18-decade axis. Fix: any bin whose significance (content/error) falls below
// a_minSignificance (default 5.0 as of the second tail-cleanup follow-up -- was 4.0,
// then 2.0 since Round 12, and 3.0 before that; tightened further per Andrew's
// request, roughly a 5-sigma requirement now) is zeroed out
// (SetBinContent/SetBinError to 0) right after normalizeProjection() runs on it, via the
// new suppressInsignificantBins() helper below -- zero-content bins are already skipped
// by ROOT on a log-y draw, so this both stops them from being drawn AND stops them from
// setting globalMinPositiveY/globalMaxY in the scan that sizes the frame, fixing the
// y-floor problem as a side effect rather than needing separate logic for it. content/
// error (the significance) is unchanged by the N_evt division and Scale(1.0,"width") in
// normalizeProjection() -- both divide content and error by the same positive factor per
// bin, which cancels in the ratio -- so applying the cut after normalization gives the
// identical result as applying it to the raw ZFitter fit yield. Suppressed-bin counts are
// summed per detector across all rapidity bins and printed once per panel (see the
// N_evt-availability prints below) so cleanup is visible in the log, not silent.
//
// AXIS RANGES: x-axis (mT-m0) upper bound is now PER-SPECIES (changed 2026-07-16 from a
// flat [0,3] for every species, per Andrew's request to "extend the max range to be
// whatever is in SetCutClass rounded up to the nearest GeV"). The real per-species
// ceiling is set in macros/SetCutClass.C's setVariableRapMtM0BinningInfo() calls
// (lines 284-292): the TPC array's own max mT-m0 is much smaller (pion 0.870123, kaon
// 0.621547, proton 0.83101 GeV/c^2 -- this is exactly why TPC curves were "stopping at
// 0.88" before the DrawFrame fix below), but BTOF and ETOF share a common, larger max
// per species (pion 3.86286, kaon 3.53668, proton 4.149 GeV/c^2 -- pulled straight from
// pionBTOFArray/kaonBTOFArray/protonBTOFArray's last entries, which match
// pionETOFArray/kaonETOFArray/protonETOFArray's last entries exactly). The frame is
// sized to the largest of these three detectors' real ceilings for that species,
// rounded up to the next whole GeV, so BTOF/ETOF curves are never clipped: Pion -> 4,
// Kaon -> 4, Proton -> 5 (see speciesMaxMtM0GeV[] in PresentZFitterSpectra() below). The
// x-axis lower bound stays 0. y-axis (yield) floor was briefly a single flat 1e-3 for
// every panel (2026-07-16), then changed same day (per Andrew's follow-up request to
// "do the same for the y floor" as the x-max fix) to auto-derive PER PANEL from that
// panel's own data: the smallest positive normalized-yield bin actually visible in the
// [0,a_maxMtM0] window, rounded DOWN to the nearest power of 10 (globalMinPositiveY
// below, floor = 10^floor(log10(globalMinPositiveY))) -- so each panel's floor sits just
// under its own real data instead of one guessed constant that might sit awkwardly
// above or far below any given panel's actual smallest values. The top of the range
// still auto-scales to 3x the tallest visible bin (globalMaxY*3) so nothing gets
// clipped off the top -- the same auto-from-data philosophy, just at the other end.
//
// COMBINING +-y RAPIDITY BINS (added 2026-07-16, Round 9, per Andrew's request): a
// symmetric collision system (O+O, like most A+A at a collider) is expected to
// produce a rapidity spectrum symmetric under y -> -y, so the +y and -y bins at the
// same |y| are, physically, two independent measurements of the same underlying
// value -- averaging them exploits that symmetry to improve effective statistics
// per |y| slice instead of drawing two separate (and separately noisy) curves.
// Implementation: SetCutClass.C books all three species' rapidity axis with an ODD
// bin count (41) uniformly spanning a range symmetric about y=0 (-2.05 to 2.05 for
// all of TPC/BTOF/ETOF, all species -- confirmed directly in setVariableRapMtM0
// BinningInfo(), lines 284-292), which guarantees bin i and bin nRapBins+1-i are
// exact mirror images of each other (bin (nRapBins+1)/2, when nRapBins is odd, sits
// exactly at y=0 and has no distinct partner). combineMirroredProjections() does the
// actual per-mT-m0-bin averaging (see its own comment for exactly what "average"
// means when one side is empty/insignificant), called once per mirror pair
// (i, nRapBins+1-i) for i=1..(nRapBins+1)/2 in drawOneChargeSign(). This roughly
// halves the number of drawn curves (e.g. 41 raw rapidity bins -> ~21 combined |y|
// bins), which is also why the color gradient and inline labels changed meaning:
// they now index |y| bins (0=center to max=edge) rather than raw rapidity bin index
// (most-negative to most-positive). A runtime check (comparing bin i's and bin
// nRapBins+1-i's actual bin centers) verifies the assumed symmetric binning rather
// than blindly trusting it, and warns once if some future SetCutClass.C change ever
// breaks that assumption for a given species.
//
// LABEL PLACEMENT FIX (Round 10, 2026-07-16, superseded by Round 11 below): the first
// Round 9 rendering piled every label on top of each other around mT-m0~0.6-0.9
// GeV/c^2, unreadable. Cause: the label-placement helper was only checked against
// tpcProj first, falling through to btofProj/etofProj only if TPC had no positive bins
// at all -- but nearly every curve DOES have TPC content (TPC's own booked range just
// ends much earlier than BTOF's, see AXIS RANGES above), so almost every label landed
// at TPC's shared ceiling regardless of how much further that same curve's BTOF tail
// actually extended. Fixed at the time by checking all three detectors and keeping
// whichever result had the largest x. This entire inline-labeling approach (Rounds
// 9-10) was then replaced by Round 11's color-scale bar (see below) -- kept here for
// history since it explains why a naive "just check TPC first" label placement is a
// trap for ANY future per-curve annotation in this macro, not just labels.
//
// COLOR-SCALE BAR (Round 11, 2026-07-16, replaces Rounds 9-10's inline labels):
// Andrew asked to replace the per-curve "y=+-X.XX" labels with "a color code on the
// right side of the plots similar to what colz does for 2d histograms." Implemented
// as drawColorScaleBar() -- a small nested TPad in the main pad's (now-widened, see
// gPad->SetRightMargin() in drawOneChargeSign()) right margin, filled with a stack of
// TBox gradient segments plus hand-drawn tick labels (see Round 12 follow-up comment
// in drawColorScaleBar() -- originally a TGaxis, replaced after it rendered nothing),
// rather than an actual TH2
// COLZ draw (there's no real 2D histogram to hang a ROOT palette off of here). The
// bar's linear 0-to-max|y| scale is exactly (not just approximately) consistent with
// how curves are colored: given SetCutClass.C's uniform rapidity binning (see above),
// a mirror pair's shared |y| magnitude is an exact affine function of the group index
// g, and so is Pass 4's colorFrac = (nGroups-g)/(nGroups-1) -- algebraically,
// colorFrac(g) == yMagnitude(g)/yMagnitude(g=1) for every g, not merely close. This
// removed the need for findCurveEndpoint() (deleted) and the per-curve TLatex labels
// entirely -- see the updated LEGEND CHOICE header comment above.
//
// FIXED COLOR RANGE (Round 12, 2026-07-16): the Round 11 scheme above normalized color
// to each PANEL's own largest |y| (yMagnitude(g=1)), which meant the same color could
// mean different |y| values on different panels (e.g. if one species/centrality's
// combined groups only reached |y|~1.5 while another reached |y|~2.0, "brightest" meant
// two different things). Andrew asked to fix the color range instead, so the same
// color always means the same |y| everywhere, making panels comparable at a glance.
// colorFrac is now TMath::Min(1.0, yMagnitude(g)/kRapColorMax) (see kRapColorMax in
// drawOneChargeSign()) -- groups whose |y| exceeds kRapColorMax clip to the brightest
// color rather than extending the scale. drawColorScaleBar() is passed this same fixed
// kRapColorMax rather than the panel's own max, so the bar's printed axis range
// matches what the colors actually mean. Range was set to [0,1.7] initially, then
// narrowed to [0,1.3] (still Round 12, same day) once CheckLastRapBin.C showed the
// real pion fit output never reaches past y=1.3 (TPC) / y=0.9 (BTOF) -- see the
// kRapColorMax declaration in drawOneChargeSign() for the full rationale.

#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TPad.h"
#include "TBox.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TMath.h"
#include <iostream>
#include <vector>
using namespace std;

//_______________________________________________________________________________
// Divide a single rapidity-bin projection by N_evt (with the same relative-error
// quadrature-sum convention ZFitter::makeSingleFitHisto() uses at ~line 1732-1751,
// reimplemented inline here so this macro doesn't need to link against PhysMath), then
// fold in 1/(2*pi), 1/(rapidity bin width), and the mT-m0 bin-width division via ROOT's
// Scale(1.0,"width"). a_rapidityBinWidth is that one rapidity bin's width (dy) taken
// from the source TH2's X-axis at this projection's rapBin -- added 2026-07-16 per
// Andrew's request to also normalize per unit rapidity, matching the dy factor
// RawSpectraModifier::convertSpectraToInvariant() (~line 2719) applies in its fuller
// invariant-yield conversion.
// If a_nEvt <= 0 (N_evt histogram wasn't found for this file), skip the N_evt division
// but still apply the 2*pi + rapidity-width + bin-width scaling, and let the caller
// know via the return value so it can warn once instead of per-projection.
bool normalizeProjection(TH1D* a_proj, double a_nEvt, double a_rapidityBinWidth){
  bool didDivideByNEvt = (a_nEvt > 0);
  if(didDivideByNEvt){
    for(int binX = 1; binX <= a_proj->GetNbinsX(); binX++){
      double content = a_proj->GetBinContent(binX);
      if(content == 0.0) continue;
      double error   = a_proj->GetBinError(binX);
      double nEvtErr = TMath::Sqrt(a_nEvt);
      double newContent = content / a_nEvt;
      double relErrContent = error / content;
      double relErrNEvt    = nEvtErr / a_nEvt;
      double newErr = newContent * TMath::Sqrt(relErrContent*relErrContent + relErrNEvt*relErrNEvt);
      if(TMath::Finite(newContent) && TMath::Finite(newErr)){
        a_proj->SetBinContent(binX,newContent);
        a_proj->SetBinError(binX,newErr);
      }else{
        a_proj->SetBinContent(binX,0.0);
        a_proj->SetBinError(binX,0.0);
      }
    }
  }
  double rapidityWidth = (a_rapidityBinWidth > 0) ? a_rapidityBinWidth : 1.0;
  a_proj->Scale(1.0/(2.0*TMath::Pi()*rapidityWidth),"width");
  return didDivideByNEvt;
}

//_______________________________________________________________________________
// Zero out any bin whose significance (content/error) is below a_minSignificance, so it
// won't be drawn (ROOT skips zero-content bins on a log-y axis) and won't influence the
// globalMaxY/globalMinPositiveY scan that sizes the frame -- see the STATISTICAL
// SIGNIFICANCE CUT header comment for the real-plot symptom (huge scattered error bars
// down to ~1e-17/1e-18) that motivated this. A bin with content>0 but error<=0 is
// treated as insignificant (degenerate/undefined significance) rather than divided by
// zero. a_minSignificance<=0 disables the cut (every already-nonzero bin passes).
// Returns how many bins were suppressed, so the caller can print a per-detector,
// per-panel count instead of silently dropping points.
int suppressInsignificantBins(TH1D* a_proj, double a_minSignificance){
  if(a_minSignificance <= 0) return 0;
  int nSuppressed = 0;
  for(int binX = 1; binX <= a_proj->GetNbinsX(); binX++){
    double content = a_proj->GetBinContent(binX);
    if(content <= 0) continue; // already empty -- nothing to suppress
    double error = a_proj->GetBinError(binX);
    double significance = (error > 0) ? (content/error) : 0.0;
    if(significance < a_minSignificance){
      a_proj->SetBinContent(binX,0.0);
      a_proj->SetBinError(binX,0.0);
      nSuppressed++;
    }
  }
  return nSuppressed;
}

//_______________________________________________________________________________
// Zero out any bin whose underlying fit had chi^2/ndf above a_maxChiSqrNdf, so a
// "confidently wrong" fit (small quoted error, but a genuinely bad fit underneath)
// doesn't survive suppressInsignificantBins() just because its error happens to be
// small -- added after Andrew spotted two isolated outlier points in the high-mT-m0
// tail that had tight, low-error markers but sat 1-1.5 decades off the trend of their
// neighbors, even after significance was already tightened to 5-sigma. Significance
// (content/error) and fit quality (chi^2/ndf) are independent things: a fit can
// converge to the wrong peak with high confidence, which raising the significance
// threshold alone can never catch.
//
// a_chiSqrHist is ZFitter's own per-bin chi^2/ndf histogram (chi^2/ndf already computed
// and stored there directly, NOT raw chi^2 -- see source/ZFitter_TPCFuncts.cxx /
// ZFitter_bTOFFuncts.cxx, which compute fitResult->GetChisquare()/fitResult->GetNDF()
// right after each individual per-(rapidity,mTm0)-bin fit), on the SAME (rapidity bin,
// mTm0 bin) axes as the yield TH2s this analysis already reads -- see
// ZFitter::makeFitDataHistos() in source/ZFitter.cxx, which clones the chi2 histograms
// from the exact same binning template as the yield spectra. a_rapBin must be the same
// rapidity-bin index already used to make a_proj (a_tpcHist/a_btofHist/a_etofHist's
// ProjectionY(...,rapBin,rapBin) call), so a_proj's bin N lines up with
// a_chiSqrHist->GetBinContent(a_rapBin,N) directly, no reindexing needed.
//
// a_chiSqrHist==NULL (e.g. an older ZFitter run predating this histogram, or this
// detector has no fits at all) or a_maxChiSqrNdf<=0 disables the cut entirely rather
// than crashing or suppressing everything. A chi^2/ndf that comes back <=0 or
// non-finite (NaN/Inf, e.g. an NDF=0 degenerate fit) is treated as a failed fit and
// suppressed, same as exceeding the threshold.
int suppressBadFitQualityBins(TH1D* a_proj, TH2* a_chiSqrHist, int a_rapBin, double a_maxChiSqrNdf){
  if(!a_proj || !a_chiSqrHist || a_maxChiSqrNdf <= 0) return 0;
  int nSuppressed = 0;
  for(int binY = 1; binY <= a_proj->GetNbinsX(); binY++){
    if(a_proj->GetBinContent(binY) <= 0) continue; // already empty -- nothing to suppress
    double chiSqrNdf = a_chiSqrHist->GetBinContent(a_rapBin,binY);
    bool badFit = !TMath::Finite(chiSqrNdf) || chiSqrNdf <= 0 || chiSqrNdf > a_maxChiSqrNdf;
    if(badFit){
      a_proj->SetBinContent(binY,0.0);
      a_proj->SetBinError(binY,0.0);
      nSuppressed++;
    }
  }
  return nSuppressed;
}

//_______________________________________________________________________________
// Combine two rapidity-bin projections that are mirror images of each other around
// y=0 into one averaged curve -- see the COMBINING +-y RAPIDITY BINS header comment
// for why this exploits the forward-backward symmetry expected in this collision
// system. a_projI/a_projJ order doesn't matter -- both are treated symmetrically.
// Per mT-m0 bin: if both sides have a positive (already-significant, post-Round-8-cut)
// entry, average them, propagating the error for an average of two independent
// measurements (err_avg = sqrt(errI^2+errJ^2)/2); if only one side has an entry (the
// other was empty, or suppressed by the significance cut), use that one value AS-IS
// rather than averaging it with an implicit zero -- doing that would artificially
// halve a real yield just because its mirror partner happened to be statistically
// insignificant, which is not what "average the two" is supposed to mean. If neither
// side has an entry, the combined bin stays empty. Returns NULL only if both inputs
// are NULL; if exactly one is NULL, returns a plain clone of the other (nothing to
// combine).
TH1D* combineMirroredProjections(TH1D* a_projI, TH1D* a_projJ, string a_newName){
  if(!a_projI && !a_projJ) return NULL;
  if(!a_projI) return (TH1D*) a_projJ->Clone(a_newName.c_str());
  if(!a_projJ) return (TH1D*) a_projI->Clone(a_newName.c_str());
  TH1D* combined = (TH1D*) a_projI->Clone(a_newName.c_str());
  for(int b = 1; b <= combined->GetNbinsX(); b++){
    double contentI = a_projI->GetBinContent(b);
    double contentJ = a_projJ->GetBinContent(b);
    double errorI = a_projI->GetBinError(b);
    double errorJ = a_projJ->GetBinError(b);
    if(contentI > 0 && contentJ > 0){
      combined->SetBinContent(b, 0.5*(contentI+contentJ));
      combined->SetBinError(b, 0.5*TMath::Sqrt(errorI*errorI + errorJ*errorJ));
    }else if(contentI > 0){
      combined->SetBinContent(b, contentI);
      combined->SetBinError(b, errorI);
    }else if(contentJ > 0){
      combined->SetBinContent(b, contentJ);
      combined->SetBinError(b, errorJ);
    }else{
      combined->SetBinContent(b, 0.0);
      combined->SetBinError(b, 0.0);
    }
  }
  return combined;
}

//_______________________________________________________________________________
// Draw a COLZ-style vertical color-scale bar in the current pad's right margin,
// mapping the same linear |y| -> color-palette fraction used for the combined curves
// (see drawOneChargeSign()'s Pass 4) onto a labeled gradient strip, in place of the
// Round 9/10 inline per-curve "y=+-X.XX" labels (Round 11, per Andrew's request:
// "instead of labeling the rapidity bins, do a color code on the right side of the
// plots similar to what colz does for 2d histograms"). Built manually (a small nested
// TPad holding a stack of TBox segments for the gradient, plus hand-drawn TLine tick
// marks and TLatex labels/title -- see the Round 12 follow-up comment inside the
// function for why this is TLine/TLatex rather than a TGaxis) rather than via an
// actual TH2 COLZ draw, since there's no real 2D histogram here to hang a palette off
// of -- this reproduces the same visual (gradient bar + labeled ticks) without needing
// one. a_parentPad is the pad this bar's own nested pad gets placed inside (one of the
// two Plus/Minus pads from PresentZFitterSpectra()'s c->Divide(2,1)); a_maxAbsY is the
// top of the |y| scale the bar's axis should show -- as of Round 12 this is always the
// fixed kRapColorMax constant (see drawOneChargeSign()), NOT this panel's own largest
// |y|, so the bar (and the colors it explains) mean the same thing on every panel
// rather than each panel silently rescaling to its own data. a_nColors must match the
// nColors used for the curves' TColor::GetColorPalette() calls so the bar's gradient
// matches their actual colors.
void drawColorScaleBar(TPad* a_parentPad, double a_maxAbsY, int a_nColors, string a_namePrefix){
  if(a_maxAbsY <= 0) return; // nothing meaningful to show a scale for (e.g. only 1 group)
  a_parentPad->cd();
  // barPad needs to be wide enough to hold BOTH the gradient strip AND the tick-label
  // text/title next to it -- an earlier version sized barPad to just the strip itself
  // (0.90-0.94 NDC), which clipped every label and the "|y|" title entirely (nothing to
  // paint them on past the pad's own right edge). Now barPad spans the bulk of the
  // margin gPad->SetRightMargin() reserved, and only the LEFT ~30% of barPad's own
  // local coordinate space (barFrac below) is the actual colored strip -- the remaining
  // ~70% is where the hand-drawn tick labels/title physically render.
  double barX1 = 0.86, barX2 = 0.99, barY1 = 0.12, barY2 = 0.88; // NDC, within a_parentPad
  TPad* barPad = new TPad((a_namePrefix+"_colorBarPad").c_str(),"",barX1,barY1,barX2,barY2);
  barPad->SetFillStyle(0);
  barPad->SetLeftMargin(0); barPad->SetRightMargin(0);
  barPad->SetTopMargin(0);  barPad->SetBottomMargin(0);
  barPad->Draw();
  barPad->cd();
  barPad->Range(0,0,1,1);
  double barFrac = 0.30; // fraction of barPad's own width that's the actual gradient strip
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
  // Tick marks + numeric labels + "|y|" title, drawn MANUALLY with TLine/TLatex
  // instead of TGaxis (Round 12 follow-up, 2026-07-18): the TGaxis-based axis that used
  // to be here rendered NOTHING in the actual rerun -- no ticks, no numbers, no title,
  // just the gradient strip and its TBox border (confirmed via the delivered PNG: the
  // border box, itself a TBox, painted fine in this exact pad, but the axis right next
  // to it painted nothing at all). This is the first time this code path had actually
  // executed -- the previous attempt was blocked by the "too many arguments" compile
  // error on SetNdivisions(), so the bug was never caught before now. TGaxis has its
  // own single-argument SetNdivisions() (confirmed via that compile error), and passing
  // it a negative value evidently does not paint anything here, unlike the documented
  // TAttAxis/TAxis "negative forces exact division count" convention it was modeled on.
  // Rather than keep guessing at TGaxis's behavior with no local ROOT available to test
  // against, this draws the same visual by hand out of primitives already proven to
  // render correctly in this pad: a short tick (TLine) plus a left-aligned numeric
  // label (TLatex) at every 0.1 units of |y|, and a separately-drawn rotated title.
  // nTicks/tickStep computed the same way the old nDiv was, so ticks still land exactly
  // every 0.1 units up to a_maxAbsY (assumed to be a multiple of 0.1).
  int nTicks = TMath::Max(1,(int)TMath::Nint(a_maxAbsY/0.1));
  double tickLen  = 0.06; // in barPad's own 0-1 local coordinate space (set by Range() above)
  double labelGap = 0.02;
  for(int i = 0; i <= nTicks; i++){
    double v = i*(a_maxAbsY/nTicks); // integer step count avoids float drift from repeated += 0.1
    double yFrac = v/a_maxAbsY;
    TLine* tick = new TLine(barFrac,yFrac,barFrac+tickLen,yFrac);
    tick->SetLineColor(kBlack);
    tick->Draw();
    // Font codes ending in "3" (precision 3) size text in FIXED PIXELS rather than as a
    // fraction of the current pad's height -- needed here because barPad is a small
    // nested sub-pad (only ~0.13 NDC wide), so a fractional size comes out tiny once
    // scaled to barPad's own small dimensions. Fixed-pixel sizing keeps labels legible
    // regardless of how narrow barPad is (same lesson as the old TGaxis version).
    TLatex* label = new TLatex(barFrac+tickLen+labelGap,yFrac,Form("%.1f",v));
    label->SetTextFont(43);
    label->SetTextSize(16);
    // Vertically-centered alignment (align 12) is fine for every INTERIOR tick, but the
    // bottom (yFrac=0) and top (yFrac=1) ticks sit exactly ON barPad's own Range(0,0,1,1)
    // boundary -- a vertically-centered label there has half its height fall outside
    // [0,1] in local coordinates, and ROOT clips anything outside the pad, so those two
    // labels came out chopped in half ("0.0" and the top value both showed as cut-off
    // digits in the actual rerun). Fix: anchor the bottom label's BOTTOM edge at yFrac=0
    // (align *1, text grows upward, never crosses below 0) and the top label's TOP edge
    // at yFrac=1 (align *3, text grows downward, never crosses above 1); every other
    // label keeps center alignment (align *2) since it isn't sitting on a hard boundary.
    int vAlign = (i == 0) ? 1 : (i == nTicks) ? 3 : 2;
    label->SetTextAlign(10+vAlign); // left-adjusted horizontally (1x), vAlign vertically
    label->Draw();
  }
  // "|y|" title, rotated 90 degrees, vertically centered along the axis and clear of
  // every numeric label (all of which are 3 characters, "0.0".."1.3", well short of
  // x=0.90 in barPad's local coordinates).
  TLatex* title = new TLatex(0.90,0.5,"|y|");
  title->SetTextFont(43);
  title->SetTextSize(20);
  title->SetTextAlign(22); // centered horizontally and vertically
  title->SetTextAngle(90);
  title->Draw();
  a_parentPad->cd(); // leave the caller back in the main (parent) pad, not the bar's
}

//_______________________________________________________________________________
// Read the per-centrality event-count histogram out of a raw PicoBinner yield file
// (NOT a ZFitter output file -- see the NORMALIZATION/"WHERE N_EVT ACTUALLY LIVES"
// header comments for why), trying the per-species name first (written when species
// were binned separately) and falling back to the bare "centEvents" name (written when
// every species was binned into one file) -- same two-name fallback ZFitter.cxx/
// RawSpectraModifier.cxx both use when they read this same histogram back out.
// Returns -1 if neither is found.
double getNumEvents(TFile* a_file, string a_speciesName, int a_centIndex){
  if(!a_file) return -1;
  TH1* centEvents = (TH1*) a_file->Get(Form("centEvents%s",a_speciesName.c_str()));
  if(!centEvents) centEvents = (TH1*) a_file->Get("centEvents");
  if(!centEvents) return -1;
  if(a_centIndex+1 < 1 || a_centIndex+1 > centEvents->GetNbinsX()) return -1;
  return centEvents->GetBinContent(a_centIndex+1);
}

//_______________________________________________________________________________
// Map a centIndex (0 = most central) to its actual percent range, for display on the
// plot title instead of the opaque "CentNN" bin index. These ranges come directly from
// macros/SetCutClass.C's cuts->setCentralities(6, centCutsArray, percents) call --
// percents = {5,10,20,40,80,100} are the UPPER edge of each of the 6 bins in centIndex
// order (see that macro's ORDERING comment for why edges must be read this way: index 0
// is the highest-refMult/most-central bin, 0-5%). Hardcoded here rather than reading
// SetCutClass.C at runtime since this is purely a display label -- if the centrality
// scheme ever changes, both places need updating together (see README's "Centrality
// binning" section, which documents this exact 6-bin scheme as of 2026-07-02/07-09).
string centralityLabel(int a_centIndex){
  static const char* kCentralityRanges[6] = {"0-5%","5-10%","10-20%","20-40%","40-80%","80-100%"};
  if(a_centIndex < 0 || a_centIndex >= 6) return Form("Cent%02d",a_centIndex); // out-of-scheme fallback
  return kCentralityRanges[a_centIndex];
}

//_______________________________________________________________________________
// Map a species name ("Pion"/"Kaon"/"Proton") + charge label ("Plus"/"Minus") to the
// actual physics symbol with a +/- superscript (e.g. "#pi^{+}"), for display in the
// legend instead of spelling the species out in English -- per Andrew's request.
// a_chargeLabel is matched by its first character ('P'/'M') since that's exactly what
// drawOneChargeSign()'s callers already pass (see PresentZFitterSpectra() below).
// Proton Minus is a special case: rather than "p^{-}", use the standard antiproton
// notation "#bar{p}" (p-bar, no separate charge sign -- the bar itself already conveys
// it) -- per Andrew's follow-up request.
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
// One pad's worth of work: given the Plus (or Minus) TPC/BTOF/ETOF TH2Ds for one
// species+centrality, project out every rapidity bin, combine each +y/-y mirror pair
// into one averaged curve (see combineMirroredProjections()/COMBINING +-y RAPIDITY
// BINS header comment), and draw all of those overlaid -- TPC/BTOF/ETOF distinguished
// by marker style, |y| distinguished by color along a rainbow gradient (dark near
// y=0, bright at the largest |y|), read off a color-scale bar in the pad's right
// margin (see drawColorScaleBar()/Round 11 header comment) rather than per-curve text.
// a_*ChiSqrHist (may be individually NULL) are ZFitter's own per-bin chi^2/ndf TH2s,
// same (rapidity bin, mTm0 bin) axes as a_*Hist -- see suppressBadFitQualityBins().
void drawOneChargeSign(TH2* a_tpcHist, TH2* a_btofHist, TH2* a_etofHist,
                        TH2* a_tpcChiSqrHist, TH2* a_btofChiSqrHist, TH2* a_etofChiSqrHist,
                        string a_speciesName, string a_chargeLabel, int a_centIndex,
                        double a_tpcNumEvents, double a_btofNumEvents, double a_maxMtM0,
                        double a_minSignificance, double a_maxChiSqrNdf){

  gPad->SetLogy();
  gPad->SetRightMargin(0.16); // widened from 0.05 (Round 11) to make room for the |y|
                              // color-scale bar drawn in place of the old inline labels
                              // -- see drawColorScaleBar() / Round 11 header comment.
  gPad->SetTopMargin(0.08);
  // Widened from ROOT's default (~0.10) to give the Y-axis title enough horizontal
  // clearance from the numeric tick labels ("10", "1", "10^-1", ...) -- confirmed via
  // the actual rerun that after CenterTitle() moved the title into the middle of the
  // axis (see below), it visually overlapped/ran through those numbers, since ROOT's
  // Y-axis TitleOffset is calibrated as a multiple of PAD HEIGHT, not width. This pad
  // is much taller than it is wide (two side-by-side panels from c->Divide(2,1)), so a
  // default offset that looks fine on a normal square pad ends up placing the title too
  // close to the axis here -- not far enough left to clear the number labels.
  gPad->SetLeftMargin(0.14);

  if(!a_tpcHist && !a_btofHist && !a_etofHist){
    TLatex missingLabel;
    missingLabel.SetTextSize(0.05);
    missingLabel.DrawLatexNDC(0.15,0.5,Form("No TPC/BTOF/ETOF spectra found for %s %s, %s",
                                             a_speciesName.c_str(),a_chargeLabel.c_str(),centralityLabel(a_centIndex).c_str()));
    return;
  }

  // Use whichever histogram is available to set up the frame/axis ranges and to
  // determine how many rapidity bins exist (TPC, BTOF, and ETOF share the same
  // rapidity binning by construction -- all cloned from the same structHisto-family
  // template in ZFitter.cxx -- so any of them is a valid reference).
  TH2* referenceHist = a_tpcHist ? a_tpcHist : (a_btofHist ? a_btofHist : a_etofHist);
  int nRapBins = referenceHist->GetNbinsX();

  // Track whether N_evt normalization actually happened (vs. was skipped because no
  // centEvents histogram was found) so we can warn once at the end instead of once per
  // rapidity bin.
  bool anyTpcProj = false, anyBtofProj = false, anyEtofProj = false;
  bool tpcNormalizedByNEvt = true, btofNormalizedByNEvt = true, etofNormalizedByNEvt = true;

  // Running totals of bins suppressed by the significance cut (see
  // suppressInsignificantBins() / STATISTICAL SIGNIFICANCE CUT header comment),
  // summed across every rapidity bin so one line can report the total per detector
  // instead of printing once per bin. Chi^2/ndf-suppressed counts (see
  // suppressBadFitQualityBins()) are tracked the same way, separately, so the two cuts'
  // effects can be reported and reasoned about independently.
  int nTpcSuppressed = 0, nBtofSuppressed = 0, nEtofSuppressed = 0;
  int nTpcChiSqrSuppressed = 0, nBtofChiSqrSuppressed = 0, nEtofChiSqrSuppressed = 0;

  // Pass 1: build the per-rapidity-bin projections (one entry per rapBin, possibly
  // NULL), normalized and significance-cut exactly as before Round 9. These get
  // combined into +-y pairs in pass 2 below -- the y-axis range (pass 3) is sized from
  // the COMBINED curves, not these raw per-bin ones, since the combined curves are
  // what's actually drawn.
  vector<TH1D*> tpcProjections, btofProjections, etofProjections;
  for(int rapBin = 1; rapBin <= nRapBins; rapBin++){
    TH1D* tpcProj  = a_tpcHist  ? a_tpcHist->ProjectionY(Form("tpcProj_%s_%s_Cent%02d_rap%03d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,rapBin),rapBin,rapBin) : NULL;
    TH1D* btofProj = a_btofHist ? a_btofHist->ProjectionY(Form("btofProj_%s_%s_Cent%02d_rap%03d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,rapBin),rapBin,rapBin) : NULL;
    TH1D* etofProj = a_etofHist ? a_etofHist->ProjectionY(Form("etofProj_%s_%s_Cent%02d_rap%03d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,rapBin),rapBin,rapBin) : NULL;

    if(tpcProj  && tpcProj->Integral()  <= 0){ delete tpcProj;  tpcProj  = NULL; }
    if(btofProj && btofProj->Integral() <= 0){ delete btofProj; btofProj = NULL; }
    if(etofProj && etofProj->Integral() <= 0){ delete etofProj; etofProj = NULL; }

    // Normalize each surviving projection: 1/(2*pi*N_evt*dy) dN/d(mT-m0). ETOF shares
    // the BTOF file's N_evt (see header comment -- ETOF is read out of the same file as
    // BTOF). Rapidity bin width (dy) is taken from each detector's own source TH2 at
    // this rapBin (they're expected to share the same rapidity binning by construction,
    // but read each one's own axis rather than assuming that). See normalizeProjection()'s
    // own comment for exactly what this does.
    if(tpcProj){  anyTpcProj  = true; tpcNormalizedByNEvt  = normalizeProjection(tpcProj,  a_tpcNumEvents,  a_tpcHist->GetXaxis()->GetBinWidth(rapBin)); }
    if(btofProj){ anyBtofProj = true; btofNormalizedByNEvt = normalizeProjection(btofProj, a_btofNumEvents, a_btofHist->GetXaxis()->GetBinWidth(rapBin)); }
    if(etofProj){ anyEtofProj = true; etofNormalizedByNEvt = normalizeProjection(etofProj, a_btofNumEvents, a_etofHist->GetXaxis()->GetBinWidth(rapBin)); }

    // Suppress statistically-insignificant bins (see STATISTICAL SIGNIFICANCE CUT
    // header comment) right after normalization -- significance (content/error) is
    // unchanged by the divisions/scaling normalizeProjection() just did, so this gives
    // the same result as cutting on the raw fit yield, just simpler to place here.
    if(tpcProj)  nTpcSuppressed  += suppressInsignificantBins(tpcProj,  a_minSignificance);
    if(btofProj) nBtofSuppressed += suppressInsignificantBins(btofProj, a_minSignificance);
    if(etofProj) nEtofSuppressed += suppressInsignificantBins(etofProj, a_minSignificance);

    // Suppress bad-fit-quality bins (see suppressBadFitQualityBins() header comment) --
    // an independent cut from significance above: a fit can converge to the wrong value
    // with a small, "confident" error, which significance alone can't catch. Order
    // relative to the significance cut doesn't matter (each only looks at bins the
    // OTHER cut hasn't already zeroed, and does nothing to an already-zero bin), so this
    // just runs right after it for locality.
    if(tpcProj)  nTpcChiSqrSuppressed  += suppressBadFitQualityBins(tpcProj,  a_tpcChiSqrHist,  rapBin, a_maxChiSqrNdf);
    if(btofProj) nBtofChiSqrSuppressed += suppressBadFitQualityBins(btofProj, a_btofChiSqrHist, rapBin, a_maxChiSqrNdf);
    if(etofProj) nEtofChiSqrSuppressed += suppressBadFitQualityBins(etofProj, a_etofChiSqrHist, rapBin, a_maxChiSqrNdf);

    tpcProjections.push_back(tpcProj);
    btofProjections.push_back(btofProj);
    etofProjections.push_back(etofProj);
  }

  if(!tpcNormalizedByNEvt && anyTpcProj)   cout << "Note: N_evt unavailable (no a_yieldFile given, or no centEvents found in it) -- "  << a_speciesName << " " << a_chargeLabel << " TPC curves are only scaled by 1/(2*pi*dy*binWidth), not divided by N_evt." << endl;
  if(!btofNormalizedByNEvt && anyBtofProj) cout << "Note: N_evt unavailable (no a_yieldFile given, or no centEvents found in it) -- " << a_speciesName << " " << a_chargeLabel << " BTOF curves are only scaled by 1/(2*pi*dy*binWidth), not divided by N_evt." << endl;
  if(!etofNormalizedByNEvt && anyEtofProj) cout << "Note: N_evt unavailable (no a_yieldFile given, or no centEvents found in it) -- " << a_speciesName << " " << a_chargeLabel << " ETOF curves are only scaled by 1/(2*pi*dy*binWidth), not divided by N_evt." << endl;

  if(a_minSignificance > 0 && (nTpcSuppressed || nBtofSuppressed || nEtofSuppressed)){
    cout << "Suppressed (significance < " << a_minSignificance << " sigma) for " << a_speciesName << " " << a_chargeLabel
         << ", Cent" << Form("%02d",a_centIndex) << ": " << nTpcSuppressed << " TPC bin(s), "
         << nBtofSuppressed << " BTOF bin(s), " << nEtofSuppressed << " ETOF bin(s) -- not drawn, "
         << "not counted toward the y-axis range." << endl;
  }
  if(a_maxChiSqrNdf > 0 && (nTpcChiSqrSuppressed || nBtofChiSqrSuppressed || nEtofChiSqrSuppressed)){
    cout << "Suppressed (chi^2/ndf > " << a_maxChiSqrNdf << ") for " << a_speciesName << " " << a_chargeLabel
         << ", Cent" << Form("%02d",a_centIndex) << ": " << nTpcChiSqrSuppressed << " TPC bin(s), "
         << nBtofChiSqrSuppressed << " BTOF bin(s), " << nEtofChiSqrSuppressed << " ETOF bin(s) -- not drawn, "
         << "not counted toward the y-axis range." << endl;
  }

  // Pass 2: combine each +y/-y mirror pair (bin i <-> bin nRapBins+1-i) into one
  // averaged curve per detector -- see COMBINING +-y RAPIDITY BINS header comment for
  // why this exact index pairing is guaranteed correct (not approximate) given how
  // SetCutClass.C books the rapidity axis. When nRapBins is odd there's one unpaired
  // central bin (i == nRapBins+1-i) sitting exactly at y=0 -- carried through as-is,
  // labeled "y=0" below rather than "y=+-0.00". A runtime check against the real bin
  // centers (not just the index arithmetic) verifies the assumed symmetry and warns
  // once if it's ever violated.
  int nGroups = (nRapBins+1)/2;
  vector<TH1D*> tpcCombined(nGroups,(TH1D*)NULL), btofCombined(nGroups,(TH1D*)NULL), etofCombined(nGroups,(TH1D*)NULL);
  vector<double> yMagnitudes(nGroups,0.0);
  bool warnedAsymmetry = false;

  for(int g = 1; g <= nGroups; g++){
    int i = g;
    int j = nRapBins+1-i;

    double centerI = referenceHist->GetXaxis()->GetBinCenter(i);
    double centerJ = referenceHist->GetXaxis()->GetBinCenter(j);
    if(i != j && !warnedAsymmetry){
      double binWidth = referenceHist->GetXaxis()->GetBinWidth(i);
      if(TMath::Abs(centerI+centerJ) > 0.1*binWidth){
        cout << "!!! Rapidity bin " << i << " (y=" << centerI << ") and its assumed mirror bin "
             << j << " (y=" << centerJ << ") for " << a_speciesName << " " << a_chargeLabel
             << " are NOT symmetric about y=0 as expected -- check SetCutClass.C's rapidity "
             << "binning for this species; combined +-y curves may not mean what they're "
             << "labeled to mean." << endl;
        warnedAsymmetry = true;
      }
    }
    yMagnitudes[g-1] = (i==j) ? 0.0 : 0.5*(TMath::Abs(centerI)+TMath::Abs(centerJ));

    string tag = Form("%s_%s_Cent%02d_grp%03d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,g);
    if(i == j){
      // Central bin at y=0 (only when nRapBins is odd) -- nothing to combine with,
      // just carry the single projection through as-is.
      tpcCombined[g-1]  = tpcProjections[i-1]  ? (TH1D*) tpcProjections[i-1]->Clone(("tpcCmb_"+tag).c_str())   : NULL;
      btofCombined[g-1] = btofProjections[i-1] ? (TH1D*) btofProjections[i-1]->Clone(("btofCmb_"+tag).c_str()) : NULL;
      etofCombined[g-1] = etofProjections[i-1] ? (TH1D*) etofProjections[i-1]->Clone(("etofCmb_"+tag).c_str()) : NULL;
    }else{
      tpcCombined[g-1]  = combineMirroredProjections(tpcProjections[i-1],  tpcProjections[j-1],  "tpcCmb_"+tag);
      btofCombined[g-1] = combineMirroredProjections(btofProjections[i-1], btofProjections[j-1], "btofCmb_"+tag);
      etofCombined[g-1] = combineMirroredProjections(etofProjections[i-1], etofProjections[j-1], "etofCmb_"+tag);
    }
  }

  // Pass 3: size the frame from the COMBINED curves (not the pre-combine per-bin
  // ones), since those are what's actually drawn below.
  double globalMaxY = 0;
  double globalMinPositiveY = 1e30;
  for(int g = 1; g <= nGroups; g++){
    for(TH1D* proj : {tpcCombined[g-1],btofCombined[g-1],etofCombined[g-1]}){
      if(!proj) continue;
      for(int b = 1; b <= proj->GetNbinsX(); b++){
        double binCenter = proj->GetXaxis()->GetBinCenter(b);
        if(binCenter < 0 || binCenter > a_maxMtM0) continue;
        double content = proj->GetBinContent(b);
        if(content <= 0) continue;
        if(content > globalMaxY) globalMaxY = content;
        if(content < globalMinPositiveY) globalMinPositiveY = content;
      }
    }
  }

  if(globalMaxY <= 0){
    TLatex emptyLabel;
    emptyLabel.SetTextSize(0.05);
    emptyLabel.DrawLatexNDC(0.15,0.5,Form("All rapidity bins empty for %s %s, %s",
                                          a_speciesName.c_str(),a_chargeLabel.c_str(),centralityLabel(a_centIndex).c_str()));
    return;
  }

  int nColors = 200;
  gStyle->SetPalette(kRainBow);

  // Fixed |y| color-scale range (Round 12, per Andrew's request), replacing the Round
  // 11 behavior of normalizing color to whatever this particular species/panel's
  // largest |y| happened to be. A fixed [0,kRapColorMax] range means the same color always means
  // the same |y| across every species/charge/centrality panel -- important for actually
  // comparing panels against each other, which a per-panel-normalized scale can't do
  // (e.g. pion's largest |y| combined group is ~2.0, so under the old scheme pion's
  // brightest color was |y|~2.0 while some other panel's brightest color could mean a
  // completely different |y|). Values above kRapColorMax clip to the brightest color
  // rather than overflowing the scale. Originally set to 1.7 (this analysis's nominal
  // rapidity acceptance edge is ~2.05, see SetCutClass.C), then lowered to 1.3 once
  // CheckLastRapBin.C showed the real fit output never actually reaches that far --
  // pion's last filled bin is TPC y=1.3, BTOF y=0.9 (out of 41 nominally booked bins
  // reaching y=2.0). 1.7 left the top ~25% of the bar as dead space no curve ever hit,
  // compressing the color contrast between real curves into the lower part of the
  // range; 1.3 matches TPC's actual reach so the full bar is used.
  const double kRapColorMax = 1.3;

  // y-axis floor: the smallest positive visible bin (globalMinPositiveY, tracked in the
  // scan above), rounded DOWN to the nearest power of 10 -- see AXIS RANGES header
  // comment (Round 7) for the reasoning. Unchanged in spirit by Round 9 -- just now
  // computed from the combined curves instead of the pre-combine ones.
  double yFloor = TMath::Power(10.0, TMath::Floor(TMath::Log10(globalMinPositiveY)));

  // Draw an explicit frame FIRST, fixed to [0,a_maxMtM0] x [yFloor, globalMaxY*3], and
  // draw every combined curve "SAME" on top of it -- see AXIS RANGES / Round 5 header
  // comment for why SetRangeUser() can't be used here instead.
  // Y-axis title changed from "1/(2#pi N_evt #Deltay) dN/d(mT-m0)" to the double-
  // differential "dN/dy d(mT-m0)" notation Andrew asked for -- physically the same
  // quantity (the "#Deltay" division is still happening in the Scale(...,"width") call
  // above, it's just no longer spelled out as a discrete bin-width divisor in the
  // label; a per-rapidity-bin-width-normalized spectrum IS the discretized dN/dy).
  // Differential is over (m_T-m_0), matching the actual x-axis variable, not just m_T
  // (fixed per Andrew's follow-up: "dN/dydmt-m0", i.e. d(m_T-m_0) not dm_T).
  // Title shows the actual centrality PERCENT RANGE (e.g. "0-5%") via centralityLabel()
  // instead of the opaque internal "CentNN" bin index -- per Andrew's request, since the
  // index alone doesn't mean anything without cross-referencing SetCutClass.C. The "--"
  // that used to separate the centrality range from the y-axis-units clause was dropped
  // per Andrew's follow-up request (just a plain comma-separated list now). Species+charge
  // now shown via particleSymbolWithCharge() (e.g. "#pi^{+}") instead of the spelled-out
  // "Pion Plus" -- per Andrew's follow-up ("the actual word pion is still there"), ROOT's
  // histogram/frame titles render the same LaTeX-like markup as axis titles and the
  // legend header, so this is consistent with where the symbol was first introduced.
  TH1F* frame = gPad->DrawFrame(0,yFloor,a_maxMtM0,globalMaxY*3,
                                 Form("%s, %s, 1/(2#pi N_{evt}) dN/dy d(m_{T}-m_{0}) vs m_{T}-m_{0}, +y/-y rapidity bins combined",
                                      particleSymbolWithCharge(a_speciesName,a_chargeLabel).c_str(),centralityLabel(a_centIndex).c_str()));
  frame->GetXaxis()->SetTitle("m_{T}-m_{0} [GeV/c^{2}]");
  frame->GetYaxis()->SetTitle("1/(2#pi N_{evt}) dN/dy d(m_{T}-m_{0}) [(GeV/c^{2})^{-1}]");
  // CenterTitle(true): WITHOUT this, ROOT's default Y-axis title placement anchors the
  // title flush against the axis's HIGH end (top of the frame) with no buffer above it
  // -- confirmed via the actual rerun's pixel data, where the title's topmost character
  // landed within ~2px of the frame's own top border. A long title string anchored
  // there with zero headroom is exactly the kind of thing that clips at the top against
  // the pad boundary. CenterTitle() moves the anchor to the middle of the axis instead,
  // giving equal buffer above and below -- the same fix already used for the color-scale
  // bar's "|y|" title earlier in this file (see drawColorScaleBar()'s CenterTitle note).
  frame->GetYaxis()->CenterTitle(true);
  // Push the title further from the axis (default offset ~1) so it clears the numeric
  // tick labels instead of running through them -- see the SetLeftMargin() comment
  // above in drawOneChargeSign() for why the default offset isn't enough clearance in
  // this particular (tall, narrow) pad. Paired with the widened left margin so there's
  // actually room for the title to move into.
  frame->GetYaxis()->SetTitleOffset(1.4);

  // Pass 4: draw every combined curve, colored by |y| (dark near y=0, bright at the
  // largest |y|) -- as of Round 11, a color-scale bar drawn after this loop (see
  // drawColorScaleBar()) communicates |y| instead of the Round 9/10 inline per-curve
  // labels, which is why there's no labeling logic left in this loop.
  for(int g = 1; g <= nGroups; g++){
    // Color by ACTUAL |y| value against the fixed kRapColorMax scale above (Round 12)
    // -- not by group index/nGroups as before Round 12, which normalized to whatever
    // this panel's own largest |y| happened to be rather than a common scale. Clamp to
    // 1.0 so a group whose |y| exceeds kRapColorMax (shouldn't normally happen given
    // this analysis's ~2.05 rapidity acceptance edge, but just in case) still gets a
    // valid color instead of an out-of-range palette index.
    double colorFrac = TMath::Min(1.0, yMagnitudes[g-1]/kRapColorMax);
    int colorIndex = TColor::GetColorPalette((int)(colorFrac*(nColors-1)));

    TH1D* tpcProj  = tpcCombined[g-1];
    TH1D* btofProj = btofCombined[g-1];
    TH1D* etofProj = etofCombined[g-1];

    if(tpcProj){
      tpcProj->SetMarkerStyle(20); // filled circle = TPC
      tpcProj->SetMarkerColor(colorIndex);
      tpcProj->SetLineColor(colorIndex);
      tpcProj->SetMarkerSize(0.7);
      // 2026-07-20: "PE SAME" (was "P SAME") -- explicitly requests the y-error bars
      // this histogram already carries. Added after Andrew flagged residual
      // "fluctuations" at the high-mT-m0/low-yield tail even after tightening the
      // significance cut. Investigated at length (see project log): in THIS ROOT build,
      // plain "P" already painted error bars identically to "PE" (confirmed by a
      // byte-for-byte pixel diff), so this line is a portability/clarity improvement,
      // not a functional fix -- don't expect a rerun to look any different. The actual
      // finding was that real per-bin uncertainties ARE present and correctly computed
      // throughout (verified via debug prints + a 5-test isolation harness), but on this
      // plot's ~6-decade log-y axis, a bin with a "normal" few-percent relative error
      // (the significance cut already requires >=10sigma, i.e. <=10% relative error) is
      // only a couple of pixels tall -- visually indistinguishable from no bar at all.
      // The residual tail scatter is real, small-but-nonzero-uncertainty statistical
      // noise rendered at the resolution limit of this axis, not a missing-error-bar bug.
      tpcProj->Draw("PE SAME"); // frame is already drawn -- always SAME now
    }
    if(btofProj){
      btofProj->SetMarkerStyle(24); // open circle = BTOF
      btofProj->SetMarkerColor(colorIndex);
      btofProj->SetLineColor(colorIndex);
      btofProj->SetMarkerSize(0.7);
      btofProj->Draw("PE SAME"); // see tpcProj's "PE SAME" comment just above
    }
    if(etofProj){
      etofProj->SetMarkerStyle(25); // open square = ETOF
      etofProj->SetMarkerColor(colorIndex);
      etofProj->SetLineColor(colorIndex);
      etofProj->SetMarkerSize(0.7);
      etofProj->Draw("PE SAME"); // see tpcProj's "PE SAME" comment just above
    }
  }

  // Color-scale bar (Round 11, fixed range since Round 12) in place of the old inline
  // "y=+-X.XX" labels -- see drawColorScaleBar() header comment. Passes the same fixed
  // kRapColorMax used to color the curves above, not this panel's own largest |y|, so
  // the bar always reflects what the colors actually mean.
  drawColorScaleBar((TPad*)gPad, kRapColorMax, nColors,
                     Form("bar_%s_%s_Cent%02d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex));

  // Caption + legend moved to the TOP-RIGHT corner, caption above legend (Round 12
  // follow-up, per Andrew's request) -- the caption used to sit at x=0.15,y=0.85 (upper
  // LEFT), which was fine back when it only had to clear the Round 11 inline per-curve
  // labels, but the actual spectra curves peak right there (near mT-m0~0, top-left of
  // the frame -- see AXIS RANGES header comment), so the caption was overlapping real
  // data points once those inline labels were removed. Curves are already low (well
  // below the top of frame) by the time you get out to the x-range this corner covers,
  // since the spectrum falls off steeply with mT-m0 -- the same reason the TPC/BTOF/ETOF
  // legend has safely lived in the top-right since Round 11 without colliding with data.
  // Right-aligned (align 31: horizontal=right, vertical=bottom) so the caption's anchor
  // point is its own bottom-right corner, hugging the same right edge (x=0.82) as the
  // legend below it, clear of the color-scale bar's reserved margin (x>0.84, see
  // gPad->SetRightMargin() above and drawColorScaleBar()'s barX1/barX2). Baseline kept
  // at y=0.88, not higher, so the text (height ~0.03 NDC, growing upward from a
  // bottom-anchored baseline) stays clear of the pad's TopMargin(0.08) boundary at
  // y=0.92 -- above that is where the panel's own title text lives (see DrawFrame()
  // above), and this caption shouldn't run into it.
  TLatex collabLabel;
  collabLabel.SetNDC();
  collabLabel.SetTextSize(0.03);
  collabLabel.SetTextAlign(31);
  // Species dropped entirely from the caption (was shown via particleSymbolOnly()) --
  // per Andrew's request, since it's redundant with the legend header (charge-specific,
  // e.g. "#pi^{+}") and the panel title above the frame, both of which already show it.
  collabLabel.DrawLatex(0.82,0.88,"STAR, O+O #sqrt{s_{NN}} = 200 GeV");

  // Marker-style legend (TPC vs BTOF vs ETOF) -- NOT one entry per curve, see header
  // comment. Built from dummy histograms so the legend swatches show the correct
  // marker style independent of which bins actually had data. Moved down slightly
  // (top edge 0.91->0.85) to leave room for the caption directly above it, same right
  // edge (x=0.82) as before so it's still clear of the color-scale bar's margin.
  // y1 lowered further (0.72->0.68) to fit the new SetHeader() row below without
  // shrinking the existing entries.
  TLegend* styleLegend = new TLegend(0.58,0.68,0.82,0.85);
  styleLegend->SetBorderSize(0);
  styleLegend->SetFillStyle(0);
  styleLegend->SetTextSize(0.035);
  // Header shows the actual particle symbol + charge (e.g. "#pi^{+}") instead of
  // relying on the panel title/caption text alone -- per Andrew's request ("use the pi
  // symbol for pion and label as + or - as needed"). See particleSymbolWithCharge().
  styleLegend->SetHeader(particleSymbolWithCharge(a_speciesName,a_chargeLabel).c_str());
  TH1D* tpcDummy = new TH1D(Form("tpcDummy_%s_%s_%d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex),"",1,0,1);
  TH1D* btofDummy = new TH1D(Form("btofDummy_%s_%s_%d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex),"",1,0,1);
  TH1D* etofDummy = new TH1D(Form("etofDummy_%s_%s_%d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex),"",1,0,1);
  tpcDummy->SetMarkerStyle(20); tpcDummy->SetMarkerColor(kBlack);
  btofDummy->SetMarkerStyle(24); btofDummy->SetMarkerColor(kBlack);
  etofDummy->SetMarkerStyle(25); etofDummy->SetMarkerColor(kBlack);
  // Detector-only labels now -- "(dE/dx)"/"(1/#beta)" dropped per Andrew's request;
  // that method-of-measurement detail isn't needed to read the plot.
  if(a_tpcHist)  styleLegend->AddEntry(tpcDummy,"TPC","P");
  if(a_btofHist) styleLegend->AddEntry(btofDummy,"BTOF","P");
  if(a_etofHist) styleLegend->AddEntry(etofDummy,"ETOF","P");
  styleLegend->Draw();
}

//_______________________________________________________________________________
void PresentZFitterSpectra(string a_tpcZFitterFile, string a_btofZFitterFile = "", int a_partIndex = 0, int a_centIndex = 0, string a_outDir = ".", string a_yieldFile = "", double a_minSignificance = 5.0, double a_maxChiSqrNdf = 0){

  const string speciesNames[3] = {"Pion","Kaon","Proton"}; // matches ParticleInfo's speciesIndex ordering for the 3 species this analysis fits
  if(a_partIndex < 0 || a_partIndex > 2){
    cout << "!!! a_partIndex must be 0 (Pion), 1 (Kaon), or 2 (Proton). Got " << a_partIndex << "." << endl;
    return;
  }
  string speciesName = speciesNames[a_partIndex];

  // Per-species x-axis (mT-m0) upper bound, rounded up to the nearest whole GeV, taken
  // from macros/SetCutClass.C's setVariableRapMtM0BinningInfo() calls (lines 284-292,
  // as of 2026-07-16) -- specifically the larger of the BTOF/ETOF arrays' last entries
  // for that species (TPC's own array tops out much lower -- see AXIS RANGES header
  // comment above for why that distinction matters). If SetCutClass.C's binning is ever
  // changed again, these three numbers need to be re-derived from it by hand -- this
  // macro has no way to read CutClass's binning directly out of the ZFitter/yield files.
  const double speciesMaxMtM0GeV[3] = {4.0, 4.0, 5.0}; // Pion (3.86286->4), Kaon (3.53668->4), Proton (4.149->5)
  double maxMtM0 = speciesMaxMtM0GeV[a_partIndex];

  TFile* tpcFile = NULL;
  if(a_tpcZFitterFile != ""){
    tpcFile = new TFile(a_tpcZFitterFile.c_str(),"READ");
    if(!tpcFile || tpcFile->IsZombie()){
      cout << "!!! Could not open TPC file " << a_tpcZFitterFile << endl;
      return;
    }
  }
  TFile* btofFile = NULL;
  if(a_btofZFitterFile != ""){
    btofFile = new TFile(a_btofZFitterFile.c_str(),"READ");
    if(!btofFile || btofFile->IsZombie()){
      cout << "!!! Could not open BTOF file " << a_btofZFitterFile << endl;
      return;
    }
  }
  if(!tpcFile && !btofFile){
    cout << "!!! Neither a TPC nor a BTOF file was given -- nothing to plot." << endl;
    return;
  }

  TH2* tpcPlus   = tpcFile  ? (TH2*) tpcFile->Get(Form("DeDx_FitData/%s/dEdxSpectra_%sPlus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* tpcMinus  = tpcFile  ? (TH2*) tpcFile->Get(Form("DeDx_FitData/%s/dEdxSpectra_%sMinus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* btofPlus  = btofFile ? (TH2*) btofFile->Get(Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%sPlus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* btofMinus = btofFile ? (TH2*) btofFile->Get(Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%sMinus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* etofPlus  = btofFile ? (TH2*) btofFile->Get(Form("ETOF_FitData/%s/InvBetaETOFSpectra_%sPlus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* etofMinus = btofFile ? (TH2*) btofFile->Get(Form("ETOF_FitData/%s/InvBetaETOFSpectra_%sMinus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;

  // Per-bin chi^2/ndf histograms, same (rapidity bin, mTm0 bin) axes as the yield TH2s
  // just above -- see ZFitter::makeFitDataHistos()/makeSpectra() in source/ZFitter.cxx
  // and suppressBadFitQualityBins()'s header comment for the object-name derivation.
  // NOTE the ETOF path's object name still says "ZbTOF_ChiSqr_..." (not "ZeTOF_..."),
  // even though it lives under ETOF_FitData -- that's a real copy-paste naming bug
  // upstream in ZFitter.cxx (SetName() reuses the BTOF chi2 histogram's name string
  // when building the ETOF one), confirmed by reading that source directly. Matching it
  // here isn't a mistake -- getting this "wrong" (i.e. using "ZeTOF") would just return
  // NULL and silently disable the ETOF chi^2 cut. Missing histograms (NULL, e.g. an
  // older ZFitter run predating this) degrade gracefully -- see suppressBadFitQualityBins().
  TH2* tpcPlusChiSqr   = tpcFile  ? (TH2*) tpcFile->Get(Form("DeDx_FitData/%s/ZTPC_ChiSqr_%sPlus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* tpcMinusChiSqr  = tpcFile  ? (TH2*) tpcFile->Get(Form("DeDx_FitData/%s/ZTPC_ChiSqr_%sMinus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* btofPlusChiSqr  = btofFile ? (TH2*) btofFile->Get(Form("BTOF_FitData/%s/ZbTOF_ChiSqr_%sPlus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* btofMinusChiSqr = btofFile ? (TH2*) btofFile->Get(Form("BTOF_FitData/%s/ZbTOF_ChiSqr_%sMinus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* etofPlusChiSqr  = btofFile ? (TH2*) btofFile->Get(Form("ETOF_FitData/%s/ZbTOF_ChiSqr_%sPlus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;
  TH2* etofMinusChiSqr = btofFile ? (TH2*) btofFile->Get(Form("ETOF_FitData/%s/ZbTOF_ChiSqr_%sMinus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex)) : NULL;

  if(!tpcPlus && !tpcMinus && !btofPlus && !btofMinus && !etofPlus && !etofMinus){
    cout << "!!! No TPC, BTOF, or ETOF spectra found for " << speciesName << ", Cent" << Form("%02d",a_centIndex)
         << " in the given file(s). Either that centrality/species wasn't fit yet in this ZFitter "
         << "run, or the file(s) predate a completed run for it (check the ZFitter log for errors)." << endl;
    return;
  }
  if(a_tpcZFitterFile == "") cout << "Note: no TPC file given -- drawing BTOF/ETOF only." << endl;
  if(a_btofZFitterFile == "") cout << "Note: no BTOF file given -- drawing TPC only." << endl;
  if(!etofPlus && !etofMinus) cout << "Note: no ETOF spectra found -- this analysis has ETOF available (m_hasETOF), but no ETOF ZFitter run has completed yet, or none exists in the BTOF file." << endl;

  // N_evt for this centrality bin. IMPORTANT (found 2026-07-16, after Andrew reported
  // the normalized plot "looks multiplied instead of divided"): centEvents is NOT
  // written into the ZFitter fit-output file at all -- traced ZFitter::makeSpectra()
  // (the function that actually writes DeDx_FitData/BTOF_FitData/ETOF_FitData, i.e. the
  // files a_tpcZFitterFile/a_btofZFitterFile point to) end to end and it never writes
  // centEvents anywhere. The only place ZFitter.cxx ever Write()s centEvents is inside
  // createQuickStartFile() (~line 1397), a completely different cache file RunZFitter.C
  // doesn't produce in the normal pipeline. centEvents<Species> is actually written by
  // PicoBinner.cxx (~line 1933, HistogramUtilities::ConditionalWrite(centEvents)) into
  // the RAW YIELD FILE (e.g. yieldHistos_OO200_pion.root) -- ZFitter only reads it back
  // out of THAT file (loadDataHistograms()/loadBasicDataHistograms(), ~line 904) for its
  // own internal per-bin normalization; it never re-writes it into its own output. So
  // N_evt has to come from a_yieldFile (the same PicoBinner output RunZFitter.C was
  // given as a_pionYieldFile/a_kaonYieldFile/a_protonYieldFile), not from tpcFile/
  // btofFile -- both of those will always come back with no centEvents histogram,
  // which is exactly the "un-normalized-looking" symptom that surfaced this bug (the
  // 1/(2*pi*dy*binWidth) factors alone, with no N_evt division to compensate, can
  // easily push values up past the raw counts instead of down). One yield file's worth
  // of centEvents covers both TPC and BTOF here, since they're the same events.
  double numEvents = -1;
  TFile* yieldFile = NULL;
  if(a_yieldFile != ""){
    yieldFile = new TFile(a_yieldFile.c_str(),"READ");
    if(!yieldFile || yieldFile->IsZombie()){
      cout << "!!! Could not open yield file " << a_yieldFile << " -- proceeding without N_evt normalization." << endl;
      yieldFile = NULL;
    }else{
      numEvents = getNumEvents(yieldFile, speciesName, a_centIndex);
    }
  }else{
    cout << "Note: no a_yieldFile given -- can't find N_evt (it only lives in the raw "
         << "PicoBinner yield file, not the ZFitter output). Curves will only be scaled "
         << "by 1/(2*pi*dy*binWidth), not divided by N_evt." << endl;
  }
  cout << "N_evt for " << speciesName << " Cent" << Form("%02d",a_centIndex) << ": " << numEvents
       << (a_yieldFile != "" ? (" (from " + a_yieldFile + ")") : "") << endl;
  double tpcNumEvents = numEvents, btofNumEvents = numEvents;

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(Form("PresentZFitterSpectra_%s_Cent%02d",speciesName.c_str(),a_centIndex), // internal object name, kept as CentNN (not user-visible)
                            Form("%s ZFitter spectra, %s",speciesName.c_str(),centralityLabel(a_centIndex).c_str()),2000,900); // window title, shown in ROOT's GUI title bar if run interactively
  c->Divide(2,1);

  c->cd(1);
  drawOneChargeSign(tpcPlus,btofPlus,etofPlus,tpcPlusChiSqr,btofPlusChiSqr,etofPlusChiSqr,speciesName,"Plus",a_centIndex,tpcNumEvents,btofNumEvents,maxMtM0,a_minSignificance,a_maxChiSqrNdf);

  c->cd(2);
  drawOneChargeSign(tpcMinus,btofMinus,etofMinus,tpcMinusChiSqr,btofMinusChiSqr,etofMinusChiSqr,speciesName,"Minus",a_centIndex,tpcNumEvents,btofNumEvents,maxMtM0,a_minSignificance,a_maxChiSqrNdf);

  string outPath = a_outDir + Form("/PresentZFitterSpectra_%s_Cent%02d.png",speciesName.c_str(),a_centIndex);
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
