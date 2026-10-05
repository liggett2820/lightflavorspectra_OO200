// ExampleGausVsStudentFit.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: Andrew asked whether the deck had a plot actually showing the
// Gaussian and Student's-t fit CURVES overlaid on a dE/dx-type distribution (as
// opposed to the "Cent%02d_RapBin_%02d.png" Overview plots, which only show the
// fitted mean+-sigma trend line over the raw 2D histogram -- fit RESULTS, not fit
// SHAPES). From reading ZFitter.cxx directly: no such plot is ever produced by the
// real pipeline. ZFitter.cxx defines three functions built for exactly this --
//   drawIndividualFits_Gaus()     (line ~3018)
//   drawIndividualFits_Student()  (line ~3043)
//   drawIndividualFits_LogNorm()  (line ~3068)
// -- but none of the three is ever CALLED anywhere in the codebase (dead code).
//
// This macro fits REAL data out of the pion yield histos file -- the same file
// RunZFitter.C takes as its a_pionYieldFile argument (see that macro's header, line
// ~55, and ZFitter::loadDataHistograms(), called on that file at RunZFitter.C ~line
// 476). Paths read: /yields/Pion_space/zTPCPlus_Pion_Cent<N> and
// /yields/Pion_space/zBTOFPlus_Pion_Cent<N>, confirmed against ZFitter::
// loadDataHistograms() (source/ZFitter.cxx ~line 951, ~line 1014-1022). A single
// (rapidity bin, mT-m0 window) slice comes out of each via TH3::ProjectionZ(name,
// rapBin, rapBin, mtLow, mtHigh) -- the same call ZFitter::grabSpecificHistogram()
// uses internally (~line 1733).
//
// FUNCTIONAL FORMS -- copied verbatim from ZFitter.cxx:
//   TPC  (ZFitter::gausAmp / ZFitter::ZTPCFunct, used because
//         m_useStudentTDistributionsForTPC = false):
//           amp * 0.3989422804 * exp(-0.5*(x-mean)^2/sigma^2) / sigma
//   BTOF (ZFitter::studentTAmp / ZFitter::ZTOFFunct_student, used because
//         m_useStudentTDistributionsForTOF = true):
//           amp * TMath::Student((x-mean)/sigma, ndf) / sigma
//         ("ndf" here -- matching ZFitter.cxx's own a_ndf and ROOT's own
//         TMath::Student(T,ndf) parameter name -- is really the Student's-t
//         distribution's continuous tail-shape parameter (usually called nu/ν in the
//         statistics literature), not a classical statistical degrees-of-freedom
//         count; it isn't restricted to integers. See ROUND 15 below. This macro's own
//         variables/labels use "nu" to avoid that ambiguity, even though the upstream
//         ROOT/ZFitter API this formula is copied from calls it ndf.)
//
// ROUND 2 (this version) -- what changed and why:
// Andrew ran round 1 against the real file and sent back the log. Two real facts came
// out of that log that round 1 had no way to know (no ROOT reachable in the session
// that wrote it): the actual axis ranges are Z_TPC=[-2,2] (600 bins) and
// Z_BTOF=[-0.5,0.5] (1000 bins), not the values guessed before. Both slices had
// healthy statistics (Integral ~985k for TPC, ~231k for BTOF) -- so the bad round-1
// result (TPC fit status 3, BTOF status 1, and a "Pion+" component parked at
// mean=-1.98654, i.e. sitting right on the zMin=-2 boundary) was NOT a statistics
// problem. It was a seeding bug: round 1's mean seed was
// `zMin + meanFracOfRange[s]*range`, which for pion (frac=0.0) evaluates to zMin
// itself (-2) -- not to the physically-motivated Z=0 the header comment actually
// argued for. That -2 seed sat at the SetParLimits() boundary Migrad was given, so
// the pion component never had anywhere to go and the fit never found the real peak
// structure.
//
// Rather than just fix that one arithmetic slip and re-guess new seed offsets blindly
// (which would just be trading one unverified guess for another), this version seeds
// from the data itself using TSpectrum::Search() -- ROOT's standard peak-finder for
// exactly this situation (seeding a multi-peak fit without knowing peak positions in
// advance). It finds up to 3 local maxima in the projected slice, sorts them by
// height, and seeds the fit's 3 species components from the found (position, height)
// pairs directly -- initial sigma for each comes from half the spacing to its nearest
// neighboring peak (or a fallback fraction of the axis range if peaks weren't found /
// there's no neighbor). The tallest found peak is assumed to be Pion+ (~70% of the
// yield by design, so it should be the dominant peak almost regardless of where it
// sits), and SetParLimits() now bounds each mean to a window around ITS OWN seed
// (not the full axis range) so a component can't wander off to a boundary the way
// pion did in round 1. The fit is also run twice (Migrad seeded from its own round-1
// result the second time) since a second pass from near the minimum often mops up
// residual non-convergence -- standard practice, not a guarantee.
//
// ROUND 3 (this version) -- what changed and why:
// Round 2 converged cleanly (fit status 0 on both panels) but Andrew sent back the
// actual PNG, and it showed the "clean convergence" was hollow: on the TPC panel, the
// Pion+ fit is excellent, but "Kaon+" was just tracing the pion peak's own falling
// shoulder, and "Proton" caught a tiny blip -- meanwhile the real data has an obvious,
// visually unmistakable secondary bump around Z~0.9-1.0 that none of the three curves
// came anywhere near. Round 2's TSpectrum threshold (0.02, i.e. a candidate peak had
// to clear 2% of the tallest peak's height) was the cause: the pion peak is ~37600
// counts tall, 2% of that is ~750, and the real secondary bump is only ~1% of the
// pion height -- so TSpectrum discarded it as noise, FindTopPeaks() fell through to
// its synthetic padding fallback for both remaining slots, and THAT'S what produced
// the "Kaon+"/"Proton" curves sitting in the wrong place with suspiciously round
// numbers (exactly 10% of the pion peak's height, printed and manually confirmed to
// match the padding formula). Fit status 0 just means Migrad found a stable minimum
// for the parameters it was given -- it says nothing about whether those parameters
// correspond to anything physically real, which they didn't here.
// Fix: threshold dropped from 0.02 to 0.003 (with slightly wider search-sigma
// smoothing to compensate) so a real ~1%-of-peak bump survives. FindTopPeaks() now
// also prints every raw peak TSpectrum finds, not just the top 3 kept for seeding, so
// the full peak landscape (including whatever's happening around Z~1.7-1.8 in the TPC
// panel, which may or may not be a real fourth feature -- this simple macro still only
// fits 3 species, matching the deck's pion/kaon/proton scope) is visible in the log
// even when it's not one of the 3 components actually fit.
//
// ROUND 6 (this version) -- what changed and why:
// Round 5's higher-mT-m0 test (mtm0Frac 0.45-0.60) was a real success for TPC: three
// genuinely separated, well-fit peaks, all three curves visibly tracking real data
// bumps. But the species LABELS came out backwards -- "Proton" landed on the peak
// closer to pion (Z~0.45), "Kaon+" on the farther one (Z~1.15). Cause: species slots
// 1 and 2 were assigned by height rank in the FindTopPeaks() output, and in that
// particular run the closer, real peak (0.463333, height 1041) fell ~0.007 short of
// minSeparationTPC and got excluded from the "kept" seed list, leaving only the
// farther, smaller peak (1.15, height 351) as the height-rank-2 slot ("Kaon+"); the
// padded fallback slot ("Proton") then drifted, purely because its fit bounds
// happened to be wide enough, onto the real closer peak during optimization. The fit
// itself was fine (both curves land on real Gaussians) -- only the identity labels
// were wrong, because pion/kaon/proton is a MASS ordering, and mass ordering should
// map onto Z-distance-from-pion ordering (kaon, being closer in mass to pion, reads
// less dE/dx separation than proton), not onto whichever peak happened to survive a
// height-based selection process. OrderRemainingByDistanceFromPion() now re-sorts
// every non-pion slot by ascending |position - pionPosition| after FindTopPeaks()
// returns, so slot 1 (kaon) is always whichever real peak sits closer to pion and
// slot 2 (proton) is always farther out, regardless of which one happened to be
// taller or which one survived the separation filter.
//
// Round 6's fix confirmed clean: re-run with the combined window (TPC 0.45-0.60,
// BTOF 0.20-0.35) gave correctly-labeled Kaon+/Proton curves on genuinely real peaks
// in both panels. The one open problem left was BTOF's total-fit-curve-rises-at-both-
// edges symptom, present since round 3 and only ever partially masked (never fixed)
// by raising the Student's-t ndf floor in rounds 4 and 6 -- both non-pion components
// were still pinned exactly at whatever floor was set, which is what "the fit wants
// something a peaked model alone can't supply" looks like across three straight
// rounds. Andrew asked explicitly to fix that properly before touching the deck.
//
// ROUND 7 (this version) -- what changed and why:
// Added a genuine flat background term to BOTH fits (not just BTOF -- TPC's own real
// PNGs have shown a similar, if less dramatic, elevated tail past the outermost peak
// in every round, so the same gap plausibly exists there too, just smaller relative
// to TPC's much taller peaks). EstimateFlatBackground() reads the level directly off
// each histogram (median of all bin contents -- robust because peaks are narrow
// relative to the full axis, so the large majority of bins sit on the flat baseline)
// rather than guessing a number. Both total-fit TF1s gained a 10th/13th flat
// parameter, seeded from that estimate, with generous but bounded limits. A dashed
// orange horizontal line for the fitted background level is now drawn and legended on
// both panels so its scale is visible directly against the peaks, not just buried in
// the fit parameters.
//
// Round 7 confirmed: BTOF's total fit curve now levels off at both edges instead of
// diverging (background fitted to a plausible value on both panels, matching the
// visible flat shelf in the real data). One new wrinkle: BTOF's Kaon+ came off the
// ndf=2 floor as hoped, but landed pinned at the ndf=30 CEILING instead -- the
// opposite direction from every prior round, i.e. it wants to be even more
// Gaussian-like than the ceiling allows now that it isn't being asked to fake tails to
// cover background.
//
// ROUND 8 (this version) -- what changed and why:
// Andrew asked to switch from a combined mT-m0 WINDOW (summing several mT-m0 bins in
// one ProjectionZ call, via a_mtm0FracLowX/a_mtm0FracHighX) to a single mT-m0 BIN
// (a_mtm0FracX), matching the already-single rapidity bin. This is worth testing
// directly against the ndf=30-ceiling finding above: dE/dx and time-of-flight
// resolution are both momentum-dependent, so each species' true Z-distribution
// (position AND width) can shift somewhat across an mT-m0 window. Summing multiple
// mT-m0 bins therefore convolves each species' distribution with its own
// momentum-dependent drift across that window -- which broadens/smears the peaks in a
// way that can look exactly like the "background shelf" EstimateFlatBackground() was
// built to soak up in round 7, or like extra heavy-tailedness (low ndf) before that.
// A single mT-m0 bin removes that smearing source entirely, so this round isolates
// whether real per-bin peaks are actually narrower / more Gaussian-like (supporting
// "combined-window smearing" as part of what round 7's background term and the ndf=30
// pinning were both partially compensating for) or whether the shape holds regardless
// (pointing to genuine background/contamination instead). Trade-off: a single bin has
// far less integrated statistics than a 15-bin-wide window, so expect noisier peaks
// and possibly a status-nonzero fit if the chosen bin is too thin -- that's expected,
// not necessarily a bug; compare the peak counts/shapes to prior rounds' logs rather
// than expecting identical statistics.
//
// This still hasn't been run by me -- there's still no ROOT reachable in this
// session; every real-data image so far (rounds 2 through 7) came from Andrew running
// it and reading the PNG back over the device bridge to his Mac, not from anything
// generated here. If TSpectrum needs an explicit `gSystem->Load("libSpectrum")` on
// your ROOT build (older ROOT5-style installs sometimes do; ROOT6 usually
// autoloads it), un-comment the Load() line right below the includes.
//
// NOTE -- signature change: the four a_mtm0FracLowTPC/a_mtm0FracHighTPC/
// a_mtm0FracLowBTOF/a_mtm0FracHighBTOF arguments are now just two:
// a_mtm0FracTPC, a_mtm0FracBTOF. Any old invocation passing all four positionally
// will now bind those values to the wrong parameters (mtm0FracBTOF and outPngPath) --
// update the call, e.g.:
//   root -l -b -q 'macros/ExampleGausVsStudentFit.C("<path>", 0, -1, -1, 0.275, 0.275)'
//
// Round 8 confirmed: both fits converge cleanly at a single mT-m0 bin, and all three
// peaks in both panels are real and correctly labeled. BTOF's Kaon+ ndf=30 ceiling
// pin PERSISTED even with cross-bin smearing removed -- evidence against
// "combined-window smearing was the cause," and more likely just a symptom of low
// per-bin statistics (Kaon+ was only ~37 counts here) leaving the tail-shape
// parameter poorly constrained. The background terms also fitted to values (~1 on
// both panels) suspiciously close to the bare Poisson floor of a mostly-empty
// histogram at this statistics level, rather than a clearly resolved shelf. Andrew's
// call: the peak fits themselves look fine as-is, so remove the background term
// entirely rather than keep chasing a number that single-bin statistics can't
// reliably resolve anyway.
//
// ROUND 9 (this version) -- what changed and why:
// Removed the flat background term added in round 7 (EstimateFlatBackground(), the
// 10th/13th fit parameter on each TF1, the fitted background printout, and the dashed
// orange background line + legend entry on both panels) per the above. Both TF1s are
// back to 9/12 pure-peak parameters, matching every round through round 6. The
// earlier edge-blowup problem (rounds 3-6) was on the combined-window default; at a
// single mT-m0 bin (round 8 on) statistics are low enough that it's not clear the
// background term was resolving anything real rather than just soaking up noise, so
// this isn't a regression back to the round-6 problem -- it's dropping a term that
// wasn't earning its keep at the current (single-bin) operating point.
//
// Round 9 confirmed: with no background term, all three BTOF ndf values settled
// interior to their [2,30] bounds (5.47, 4.38, 3.30) -- no pinning at either boundary,
// for the first time across all 9 rounds. TPC unchanged (still 3 clean peaks).
//
// ROUND 10 (this version) -- what changed and why:
// Andrew asked to double-check the plot titles carry the EXACT rapidity and mT-m0
// ranges, not just the bin index (the titles/log had said "RapBin 21, mT-m0 bin 27",
// which doesn't say what y or mT-m0 range that bin actually spans). Both TH3Is are
// booked with variable-width bin edges, not fixed ones (see PicoBinner.cxx's
// tpcRapEdges/mtm0Edges/btofRapEdges/btofMtM0Edges, built via
// HistogramUtilities::getBinEdges() from getRapMtM0BinStructure()) -- so a bin's
// physical width isn't inferable from its index alone, and has to come from
// GetXaxis()/GetYaxis()->GetBinLowEdge()/GetBinUpEdge() on the actual histogram.
// X-axis = rapidity (dimensionless "y"), Y-axis = mT-m0 (GeV/c^2), confirmed against
// PicoBinner.cxx's TH3I constructor argument order and the "m_{T}-m_{0} [GeV/c^{2}]"
// axis label convention used elsewhere in the codebase (EfficiencyFitter.cxx). Titles
// and the log now print the exact [low,high] edges for both axes instead of a bare
// bin number.
//
// ROUND 11 (this version) -- what changed and why:
// Three cosmetic changes, all per Andrew's request now that the fit itself (round 9)
// and the physical ranges in the titles (round 10) are settled:
// - Dropped "-- sum-of-Gaussians fit" / "-- sum-of-Student's-t fit" off the end of
//   both plot titles.
// - Dropped the parenthetical roles ("(tallest peak)", "(closer to pion)",
//   "(farther from pion)") from the LEGEND only -- added speciesLegendLabels[] for
//   that ("Pion+"/"Kaon+"/"Proton"), while keeping the original, more descriptive
//   speciesNames[] for the console log printout, where the extra context is still
//   useful for debugging peak assignment.
// - Replaced the bare "Cent 0" in both titles with the actual centrality percent
//   range ("0-5%"), via a new CentralityLabel() helper -- same 6-bin scheme
//   ({"0-5%","5-10%","10-20%","20-40%","40-80%","80-100%"}, index 0 = most
//   central) already used by macros/PresentZFitterSpectra.C's own
//   centralityLabel() for the deck's other centrality-binned plots, copied here
//   rather than shared via #include since these are independent ROOT macros.
//
// ROUND 12 (this version) -- what changed and why:
// Andrew asked for the fitted parameters themselves to appear on the PNG, not just in
// the console log (which isn't part of the deliverable image). Added a TPaveText per
// panel -- amp/mean/sigma per species for TPC, plus ndf for BTOF (the parameter this
// whole BTOF-iteration thread has been about) -- positioned upper-left in NDC
// (pad-relative) coordinates to clear both the legend (upper-right) and the pion peak,
// which has sat near the center of the x-axis in every real run so far. No explicit
// text size set, so ROOT auto-fits the text to the box for however many lines are
// added; transparent fill/no border to match the legend's own style.
//
// ROUND 13 (this version) -- what changed and why:
// Round 12's fit-parameters box came back with the font too small -- letting ROOT
// auto-fit the text to the box height (no explicit SetTextSize) shrank it hard to
// cram 4 lines in. Fixed with an explicit SetTextSize(0.032) on both panels' boxes
// (same rough scale as the legend text), and widened the box (x2: 0.42 -> 0.50) so
// BTOF's longer lines (which add "ndf=...") have room without wrapping or overflowing
// past the box edge.
//
// ROUND 14 (this version) -- what changed and why:
// Round 13 was still broken -- worse, actually: Andrew reported the text now
// overlapping the plot. Two compounding causes, found by actually working out the
// NDC-to-data coordinate mapping instead of eyeballing it: (1) round 13's box right
// edge (NDC x=0.50) maps to roughly data x=0.0 under ROOT's default pad margins,
// which is almost exactly the pion peak's mean -- so the box itself was sitting on
// top of the peak, not just near it; (2) each species' one long combined line ("name:
// amp=... mean=... sigma=...") was wide enough at the explicit 0.032 text size to
// overflow past even that box's right edge, into the legend. Fixed by actually
// computing a safe box x-range from the axis limits and default margins (x2=0.40 now
// maps to roughly data x=-0.4/-0.13 for TPC/BTOF respectively, clear of every real
// peak seen in any round) and by never emitting a single line long enough to overflow
// in the first place -- one value per line (name / A / mu / sigma[/ ndf]) instead of
// one combined line per species. This trades line count for line width, which is the
// right trade here: the box's x-range no longer overlaps the peak, so there's no cost
// to making it taller to fit more short lines.
//
// ROUND 15 (this version) -- what changed and why:
// Andrew asked why the BTOF tail-shape parameter came out non-integer (e.g. 4.38) --
// answer: "ndf" is a misleading name here. It only means an integer count in the
// classical derivation of the Student's-t distribution (sample-mean/sample-stddev of
// n normal samples, degrees of freedom = n-1). But TMath::Student(T,ndf) is just the
// resulting PDF formula, a smooth function of its shape parameter for any real value >
// 0 (the Gamma function it's built from is defined for non-integer arguments). ZFitter
// uses it purely as a flexible peak shape -- Migrad fits that parameter continuously
// like any other, with no knowledge of (or obligation to) its historical name. To stop
// implying a phantom integer constraint, renamed the local variable (ndfSeed ->
// nuSeed), the console log label ("ndf=" -> "nu="), and the on-plot label ("ndf =" ->
// "#nu =", matching the #mu/#sigma already used) to "nu", the standard statistics
// symbol for this exact continuous shape parameter. Did NOT rename ROOT's own
// TMath::Student(T,ndf) call or the upstream ZFitter.cxx source this formula is copied
// from -- those aren't this macro's to rename, and the header's FUNCTIONAL FORMS note
// above now flags the naming mismatch explicitly instead.
//
// ROUND 16 (this version) -- what changed and why:
// Andrew noticed the printed "A" values don't match the peak heights visible on the
// curves (e.g. TPC pion A=401.3, but the curve clearly peaks around 2400-2500). This
// is because A is the AREA-normalization constant of the functional form -- for the
// Gaussian, integral = A, and the value AT the peak (x=mean) works out to
// A/(sigma*sqrt(2*pi)); for the Student's-t form, peak height = A*TMath::Student(0,
// nu)/sigma (TMath::Student(0,nu) plays the same "density at the mode" role, just
// nu-dependent instead of the Gaussian's fixed 1/sqrt(2pi)). Checking the pion TPC
// numbers above: 401.293/(0.0674024*sqrt(2*pi)) = ~2375, which matches the curve's
// visible peak. This isn't a bug in the fit -- A is exactly the parameter ZFitter's
// own gausAmp()/studentTAmp() use, and it's what's fit by Migrad -- it's just not the
// same number as "how tall does the curve look." Added a "peak = ..." line (computed
// from A/sigma[/nu]) to both the on-plot parameter box and the console log, right
// after A, so both quantities are visible and it's unambiguous which is which.

#include <vector>
#include <algorithm>
#include <utility>
#include <string>

// gSystem->Load("libSpectrum"); // uncomment if TSpectrum isn't found (see note above)

//=========================================================================
// Finds up to a_nPeaksWanted local maxima in a_hIn via TSpectrum, sorted by
// DESCENDING peak height (tallest first). Returns (position, height) pairs. If
// TSpectrum finds fewer than a_nPeaksWanted peaks, the list is padded out by
// offsetting from the tallest found peak (or the histogram's own max bin, if
// TSpectrum found nothing at all) so callers always get exactly a_nPeaksWanted
// entries to seed with -- the padded entries are marked with a small padded
// amplitude so a fit component that isn't really there can shrink toward zero
// instead of chasing real yield away from a genuine peak.
//=========================================================================
std::vector< std::pair<double,double> > FindTopPeaks(TH1D* a_hIn, int a_nPeaksWanted,
    double a_searchSigmaBins, double a_threshold, double a_minSeparation){

  TSpectrum spectrum(4*a_nPeaksWanted + 4);
  int nFound = spectrum.Search(a_hIn, a_searchSigmaBins, "nobackground", a_threshold);
  double* peakX = spectrum.GetPositionX();
  double* peakY = spectrum.GetPositionY();

  std::vector< std::pair<double,double> > peaks; // (height, position)
  for(int i = 0; i < nFound; i++){
    peaks.push_back(std::make_pair(peakY[i], peakX[i]));
  }
  std::sort(peaks.begin(), peaks.end(), std::greater< std::pair<double,double> >());

  cout << "  [" << a_hIn->GetName() << "] TSpectrum found " << nFound
       << " raw peaks (minSeparation=" << a_minSeparation << ") (position, height):" << endl;
  for(size_t i = 0; i < peaks.size(); i++){
    cout << "    " << peaks[i].second << ", " << peaks[i].first << endl;
  }

  // ROUND 4: greedy-select tallest-first, but SKIP a candidate if it's within
  // a_minSeparation of a peak already selected. Without this, a tall shoulder/
  // inflection artifact sitting right next to the dominant peak can outrank a
  // smaller but genuinely separated real peak further out -- which is exactly what
  // happened in round 3 (a "Proton" component that just traced the pion's own
  // shoulder, while a real bump around Z~1.7-1.8 in the TPC panel was never even
  // considered because it ranked below that shoulder artifact in raw height).
  std::vector< std::pair<double,double> > result; // (position, height), tallest first
  for(size_t i = 0; i < peaks.size() && (int)result.size() < a_nPeaksWanted; i++){
    double candX = peaks[i].second;
    double candY = peaks[i].first;
    bool tooClose = false;
    for(size_t j = 0; j < result.size(); j++){
      if(TMath::Abs(candX - result[j].first) < a_minSeparation){ tooClose = true; break; }
    }
    if(tooClose) continue;
    result.push_back(std::make_pair(candX, candY));
  }
  cout << "  [" << a_hIn->GetName() << "] kept after separation filter:" << endl;
  for(size_t i = 0; i < result.size(); i++){
    cout << "    " << result[i].first << ", " << result[i].second << endl;
  }

  double rangeLow  = a_hIn->GetXaxis()->GetXmin();
  double rangeHigh = a_hIn->GetXaxis()->GetXmax();
  double fallbackX = a_hIn->GetXaxis()->GetBinCenter(a_hIn->GetMaximumBin());
  double fallbackY = a_hIn->GetMaximum();
  if(result.empty()){
    result.push_back(std::make_pair(fallbackX, fallbackY));
  }
  int padIndex = 1;
  while((int)result.size() < a_nPeaksWanted){
    double anchorX = result[0].first;
    double anchorY = result[0].second;
    double offset = 0.08*(rangeHigh-rangeLow)*padIndex;
    double padX = anchorX + offset;
    if(padX > rangeHigh) padX = anchorX - offset;
    result.push_back(std::make_pair(padX, 0.1*anchorY));
    padIndex++;
  }
  return result;
}

//=========================================================================
// ROUND 6: reorders a FindTopPeaks() result so that species identity is assigned by
// PHYSICAL proximity to the pion peak, not by height rank. Round 5's higher-mT-m0
// TPC result exposed why height rank isn't good enough: the separation filter
// happened to exclude the real, closer-to-pion peak from the "kept" seed list (it
// fell ~0.007 short of minSeparationTPC), leaving only the farther real peak in the
// height-ranked list -- so "Kaon+" (slot 1, by height rank) ended up on the FARTHER
// peak while the padded "Proton" slot (slot 2) drifted, during fitting, onto the
// CLOSER real peak purely because its wide parameter limits happened to reach it.
// The fit itself converged fine either way (both are real Gaussians on real data) --
// only the LABELING was backwards, because pion-kaon-proton is a mass ordering, and
// mass ordering should map to Z-distance-from-pion ordering (kaon, being closer in
// mass to pion, should show less dE/dx separation than proton). Element 0 (tallest =
// pion, essentially guaranteed both by abundance and by construction: Z is defined
// relative to the pion hypothesis, so pion sits at Z~0 almost regardless of
// statistics) is left in place; every other element is re-sorted by ascending
// distance from element 0's position, so slot 1 (kaon) is always whichever peak is
// physically closest to pion and slot 2 (proton) is always farther out.
//=========================================================================
void OrderRemainingByDistanceFromPion(std::vector< std::pair<double,double> >& a_peaks){
  if(a_peaks.size() < 2) return;
  double pionX = a_peaks[0].first;
  std::sort(a_peaks.begin()+1, a_peaks.end(),
    [pionX](const std::pair<double,double>& a, const std::pair<double,double>& b){
      return TMath::Abs(a.first-pionX) < TMath::Abs(b.first-pionX);
    });
}

//=========================================================================
// ROUND 11: maps a centIndex (0 = most central) to its actual percent range, for
// display instead of the opaque "Cent0" bin index -- same helper (and same 6-bin
// scheme) as PresentZFitterSpectra.C's centralityLabel(), copied here rather than
// shared via #include since these are independent ROOT macros. Ranges come from
// macros/SetCutClass.C's cuts->setCentralities(6, centCutsArray, percents) call:
// percents = {5,10,20,40,80,100} are the UPPER edge of each of the 6 bins in centIndex
// order, index 0 being the highest-refMult/most-central bin (0-5%). If the centrality
// scheme ever changes, this and PresentZFitterSpectra.C's copy both need updating.
//=========================================================================
std::string CentralityLabel(int a_centIndex){
  static const char* kCentralityRanges[6] = {"0-5%","5-10%","10-20%","20-40%","40-80%","80-100%"};
  if(a_centIndex < 0 || a_centIndex >= 6) return Form("Cent%02d", a_centIndex); // out-of-scheme fallback
  return kCentralityRanges[a_centIndex];
}

void ExampleGausVsStudentFit(
    const char* a_pionYieldFile = "/Users/aliggett/data/OO/yieldHistos_OO200_pion.root",
    int a_centIndex = 0,
    int a_rapBinTPC = -1,   // -1 = auto-pick the middle rapidity bin
    int a_rapBinBTOF = -1,  // -1 = auto-pick the middle rapidity bin
    // ROUND 5 exposed these as parameters instead of hardcoded fractions. ROUND 8
    // collapsed the old Low/High WINDOW pair into a single BIN fraction each, per
    // Andrew's request -- see the ROUND 8 header note above for why (isolating
    // momentum-dependent smearing from combining multiple mT-m0 bins). 0.275 sits in
    // the middle of the old default window (0.20-0.35); try higher (e.g. 0.45-0.60)
    // for a higher-momentum slice, same as before.
    double a_mtm0FracTPC = 0.275,
    double a_mtm0FracBTOF = 0.275,
    const char* a_outPngPath = "ExampleGausVsStudentFit.png"){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameBorderMode(0);

  const int nSpecies = 3;
  // ROUND 11: kept the parenthetical roles here (used only in the log printout, where
  // they're useful diagnostic context) but dropped them from the legend -- see
  // speciesLegendLabels below, used on the plot itself per Andrew's request.
  const char* speciesNames[nSpecies] = {"Pion+ (tallest peak)", "Kaon+ (closer to pion)", "Proton (farther from pion)"};
  const char* speciesLegendLabels[nSpecies] = {"Pion+", "Kaon+", "Proton"};
  int speciesColors[nSpecies] = {kRed+1, kBlue+1, kGreen+2};

  TFile* inFile = new TFile(a_pionYieldFile, "READ");
  if(!inFile || inFile->IsZombie()){
    cout << "ERROR: could not open " << a_pionYieldFile << endl;
    return;
  }

  TString tpcPath  = Form("/yields/Pion_space/zTPCPlus_Pion_Cent%d",  a_centIndex);
  TString btofPath = Form("/yields/Pion_space/zBTOFPlus_Pion_Cent%d", a_centIndex);

  TH3I* h3TPC  = (TH3I*) inFile->Get(tpcPath);
  TH3I* h3BTOF = (TH3I*) inFile->Get(btofPath);

  if(!h3TPC){
    cout << "ERROR: couldn't find " << tpcPath << " in " << a_pionYieldFile << endl;
    return;
  }
  if(!h3BTOF){
    cout << "ERROR: couldn't find " << btofPath << " in " << a_pionYieldFile << endl;
    return;
  }

  cout << "TPC  histo: " << tpcPath  << "  nRapBins=" << h3TPC->GetNbinsX()
       << " nMtm0Bins=" << h3TPC->GetNbinsY()  << " nZBins=" << h3TPC->GetNbinsZ()
       << " Zrange=[" << h3TPC->GetZaxis()->GetXmin() << "," << h3TPC->GetZaxis()->GetXmax() << "]" << endl;
  cout << "BTOF histo: " << btofPath << "  nRapBins=" << h3BTOF->GetNbinsX()
       << " nMtm0Bins=" << h3BTOF->GetNbinsY() << " nZBins=" << h3BTOF->GetNbinsZ()
       << " Zrange=[" << h3BTOF->GetZaxis()->GetXmin() << "," << h3BTOF->GetZaxis()->GetXmax() << "]" << endl;

  //=========================================================================
  // Pull one (rapidity, mT-m0-window) slice out of each 3D histogram.
  //=========================================================================

  // ROUND 8: single mT-m0 BIN, not a summed window -- mtLow/mtHigh kept as separate
  // variables (both equal to the same bin) purely so the rest of the code below
  // (ProjectionZ calls, histogram naming, titles) didn't need touching.
  int rapBinTPC = (a_rapBinTPC > 0) ? a_rapBinTPC : (h3TPC->GetNbinsX()+1)/2;
  int mtBinTPC  = TMath::Max(1, (int)(a_mtm0FracTPC*h3TPC->GetNbinsY()));
  int mtLowTPC  = mtBinTPC;
  int mtHighTPC = mtBinTPC;

  int rapBinBTOF = (a_rapBinBTOF > 0) ? a_rapBinBTOF : (h3BTOF->GetNbinsX()+1)/2;
  int mtBinBTOF  = TMath::Max(1, (int)(a_mtm0FracBTOF*h3BTOF->GetNbinsY()));
  int mtLowBTOF  = mtBinBTOF;
  int mtHighBTOF = mtBinBTOF;

  // ROUND 10: exact physical (rapidity, mT-m0) bin edges, for the titles and log --
  // Andrew asked to double-check these against the real bin edges rather than just
  // printing bin INDICES. Both TH3Is are booked with variable-width bin edges (see
  // PicoBinner.cxx's tpcRapEdges/mtm0Edges/btofRapEdges/btofMtM0Edges vectors, read
  // via HistogramUtilities::getBinEdges()), so a bin's width isn't a fixed fraction of
  // the axis -- GetBinLowEdge()/GetBinUpEdge() read the actual edges ROOT stored for
  // that histogram, which is correct regardless of uniform vs variable binning. X-axis
  // = rapidity (dimensionless, "y"), Y-axis = mT-m0 (GeV/c^2) -- confirmed against
  // PicoBinner.cxx's TH3I booking order (nBinsX/xEdges=rapidity, nBinsY/yEdges=mT-m0,
  // nBinsZ/zEdges=Z) and the "m_{T}-m_{0} [GeV/c^{2}]" axis label used elsewhere in
  // the codebase (e.g. EfficiencyFitter.cxx).
  double rapLowTPC   = h3TPC->GetXaxis()->GetBinLowEdge(rapBinTPC);
  double rapHighTPC  = h3TPC->GetXaxis()->GetBinUpEdge(rapBinTPC);
  double mtLowValTPC  = h3TPC->GetYaxis()->GetBinLowEdge(mtBinTPC);
  double mtHighValTPC = h3TPC->GetYaxis()->GetBinUpEdge(mtBinTPC);

  double rapLowBTOF   = h3BTOF->GetXaxis()->GetBinLowEdge(rapBinBTOF);
  double rapHighBTOF  = h3BTOF->GetXaxis()->GetBinUpEdge(rapBinBTOF);
  double mtLowValBTOF  = h3BTOF->GetYaxis()->GetBinLowEdge(mtBinBTOF);
  double mtHighValBTOF = h3BTOF->GetYaxis()->GetBinUpEdge(mtBinBTOF);

  TH1D* hTPC = h3TPC->ProjectionZ(
    Form("hTPC_Pion_Cent%d_RapBin%d_mTm0_%d_%d", a_centIndex, rapBinTPC, mtLowTPC, mtHighTPC),
    rapBinTPC, rapBinTPC, mtLowTPC, mtHighTPC);
  hTPC->SetDirectory(0);

  TH1D* hBTOF = h3BTOF->ProjectionZ(
    Form("hBTOF_Pion_Cent%d_RapBin%d_mTm0_%d_%d", a_centIndex, rapBinBTOF, mtLowBTOF, mtHighBTOF),
    rapBinBTOF, rapBinBTOF, mtLowBTOF, mtHighBTOF);
  hBTOF->SetDirectory(0);

  cout << "TPC  slice: rapBin=" << rapBinTPC  << " y=[" << rapLowTPC  << "," << rapHighTPC
       << "]  mT-m0 bin=" << mtBinTPC  << " mT-m0=[" << mtLowValTPC  << "," << mtHighValTPC
       << "] GeV/c^2  Integral=" << hTPC->Integral()  << endl;
  cout << "BTOF slice: rapBin=" << rapBinBTOF << " y=[" << rapLowBTOF << "," << rapHighBTOF
       << "]  mT-m0 bin=" << mtBinBTOF << " mT-m0=[" << mtLowValBTOF << "," << mtHighValBTOF
       << "] GeV/c^2  Integral=" << hBTOF->Integral() << endl;

  if(hTPC->Integral() <= 0){
    cout << "ERROR: TPC slice is empty -- try a different a_centIndex/a_rapBinTPC." << endl;
    return;
  }
  if(hBTOF->Integral() <= 0){
    cout << "ERROR: BTOF slice is empty -- try a different a_centIndex/a_rapBinBTOF." << endl;
    return;
  }

  //=========================================================================
  // TPC panel: sum-of-Gaussians fit to the real Z_TPC slice, seeded from
  // TSpectrum-found peaks instead of guessed fractions of the axis range.
  //=========================================================================

  double zMinTPC = hTPC->GetXaxis()->GetXmin();
  double zMaxTPC = hTPC->GetXaxis()->GetXmax();
  double binWidthTPC = hTPC->GetBinWidth(1);
  double rangeTPC = zMaxTPC - zMinTPC;

  // search sigma in bins: guess ~2% of the axis width's worth of bins as a typical
  // peak width -- 600 bins over range 4 means ~12 bins; adjust if the printed peak
  // list below looks wrong (too many spurious peaks = raise this; peaks merged into
  // one = lower it). ROUND 3: threshold dropped from 0.02 to 0.003 -- round 2's real
  // plot showed a genuine secondary bump at ~1% of the pion peak's height (missed
  // entirely at the old 2% threshold, which is what forced the padding fallback to
  // kick in and produce a "Kaon+" curve that was really just tracing the pion's own
  // shoulder). 0.003 leaves a safety margin below that ~1% real feature while still
  // being well above single-bin statistical noise, which TSpectrum's deconvolution
  // already suppresses independently of this threshold.
  double searchSigmaBinsTPC = 0.02*hTPC->GetNbinsX();
  // ROUND 4: candidates within minSeparationTPC of an already-picked (taller) peak
  // are skipped -- see FindTopPeaks' header comment. 6x the search-sigma width, in Z
  // units, is enough to clear a shoulder/inflection artifact right next to the pion
  // peak while still letting a genuinely separated peak (the real Z~0.9-1.0 kaon
  // bump, ~13x pion's own fitted sigma away in round 2's log) through.
  double minSeparationTPC = 6.0*searchSigmaBinsTPC*binWidthTPC;
  std::vector< std::pair<double,double> > peaksTPC = FindTopPeaks(hTPC, nSpecies, searchSigmaBinsTPC, 0.003, minSeparationTPC);
  OrderRemainingByDistanceFromPion(peaksTPC);
  cout << "TPC peaks after distance-from-pion reordering (position, height):" << endl;
  for(size_t i = 0; i < peaksTPC.size(); i++){
    cout << "  " << peaksTPC[i].first << ", " << peaksTPC[i].second << endl;
  }
  std::vector<double> peakPosSortedTPC;
  for(size_t i = 0; i < peaksTPC.size(); i++) peakPosSortedTPC.push_back(peaksTPC[i].first);
  std::sort(peakPosSortedTPC.begin(), peakPosSortedTPC.end());

  TF1* totalGausFit = new TF1("totalGausFit",
    "[0]*0.3989422804*TMath::Exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))/[2]"
    "+[3]*0.3989422804*TMath::Exp(-0.5*(x-[4])*(x-[4])/([5]*[5]))/[5]"
    "+[6]*0.3989422804*TMath::Exp(-0.5*(x-[7])*(x-[7])/([8]*[8]))/[8]",
    // ROUND 9: flat background term (round 7-8) removed -- see header note above
    zMinTPC, zMaxTPC);
  for(int s = 0; s < nSpecies; s++){
    double meanSeed = peaksTPC[s].first;
    double heightSeed = peaksTPC[s].second;
    // nearest-neighbor spacing among ALL found peaks, for a data-driven sigma seed
    double nearestSpacing = rangeTPC;
    for(size_t j = 0; j < peakPosSortedTPC.size(); j++){
      if(peakPosSortedTPC[j] == meanSeed) continue;
      double d = TMath::Abs(peakPosSortedTPC[j]-meanSeed);
      if(d < nearestSpacing) nearestSpacing = d;
    }
    double sigmaSeed = TMath::Max(2.0*binWidthTPC, TMath::Min(0.25*nearestSpacing, 0.1*rangeTPC));
    double ampSeed = heightSeed * TMath::Sqrt(2.0*TMath::Pi()) * sigmaSeed; // invert peak-height formula
    totalGausFit->SetParameter(3*s+0, ampSeed);
    totalGausFit->SetParameter(3*s+1, meanSeed);
    totalGausFit->SetParameter(3*s+2, sigmaSeed);
    totalGausFit->SetParLimits(3*s+0, 0.0, 20.0*ampSeed + 1.0);
    totalGausFit->SetParLimits(3*s+1, meanSeed - 0.5*nearestSpacing, meanSeed + 0.5*nearestSpacing);
    totalGausFit->SetParLimits(3*s+2, 0.5*binWidthTPC, TMath::Max(0.3*nearestSpacing, 3.0*binWidthTPC));
  }
  TFitResultPtr fitResultTPC = hTPC->Fit(totalGausFit, "RQS0");
  fitResultTPC = hTPC->Fit(totalGausFit, "RQS0"); // second pass from the first minimum
  cout << "TPC fit status: " << fitResultTPC->Status()
       << " (0 = converged cleanly; nonzero = check peak list / searchSigmaBinsTPC above)" << endl;

  TF1* gausComponent[nSpecies];
  for(int s = 0; s < nSpecies; s++){
    double fittedMean  = totalGausFit->GetParameter(3*s+1);
    double fittedSigma = totalGausFit->GetParameter(3*s+2);
    gausComponent[s] = new TF1(Form("gausComponent%d", s),
      "0.3989422804*[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))/[2]",
      fittedMean - 6.0*fittedSigma, fittedMean + 6.0*fittedSigma);
    gausComponent[s]->SetParameter(0, totalGausFit->GetParameter(3*s+0));
    gausComponent[s]->SetParameter(1, fittedMean);
    gausComponent[s]->SetParameter(2, fittedSigma);
    gausComponent[s]->SetNpx(300);
    gausComponent[s]->SetLineColor(speciesColors[s]);
    gausComponent[s]->SetLineWidth(2);
    // ROUND 16: peak height (curve value at x=mean) added alongside amp -- amp is the
    // AREA-normalization constant (matches ZFitter's gausAmp()), not the height you'd
    // read off the plot; see header note. peak = amp/(sigma*sqrt(2pi)).
    double peakHeightLogTPC = totalGausFit->GetParameter(3*s+0) * 0.3989422804 / fittedSigma;
    cout << "  " << speciesNames[s] << " (TPC): amp=" << totalGausFit->GetParameter(3*s+0)
         << " mean=" << fittedMean << " sigma=" << fittedSigma
         << " peakHeight=" << peakHeightLogTPC << endl;
  }
  totalGausFit->SetLineColor(kBlack);
  totalGausFit->SetLineStyle(2);
  totalGausFit->SetLineWidth(2);
  totalGausFit->SetNpx(500);

  //=========================================================================
  // BTOF panel: sum-of-Student's-t fit to the real Z_BTOF slice, same
  // TSpectrum-seeding approach.
  //=========================================================================

  double zMinBTOF = hBTOF->GetXaxis()->GetXmin();
  double zMaxBTOF = hBTOF->GetXaxis()->GetXmax();
  double binWidthBTOF = hBTOF->GetBinWidth(1);
  double rangeBTOF = zMaxBTOF - zMinBTOF;

  double searchSigmaBinsBTOF = 0.02*hBTOF->GetNbinsX();
  // BTOF peaks sit much closer together than TPC's do (round 3's log had kaon only
  // ~0.1 from pion, out of a total 1.0-wide axis), so a 6x multiplier here would
  // reject the real kaon peak. 3x is enough to clear a shoulder artifact without
  // excluding that close-in real structure.
  double minSeparationBTOF = 3.0*searchSigmaBinsBTOF*binWidthBTOF;
  std::vector< std::pair<double,double> > peaksBTOF = FindTopPeaks(hBTOF, nSpecies, searchSigmaBinsBTOF, 0.003, minSeparationBTOF);
  OrderRemainingByDistanceFromPion(peaksBTOF);
  cout << "BTOF peaks after distance-from-pion reordering (position, height):" << endl;
  for(size_t i = 0; i < peaksBTOF.size(); i++){
    cout << "  " << peaksBTOF[i].first << ", " << peaksBTOF[i].second << endl;
  }
  std::vector<double> peakPosSortedBTOF;
  for(size_t i = 0; i < peaksBTOF.size(); i++) peakPosSortedBTOF.push_back(peaksBTOF[i].first);
  std::sort(peakPosSortedBTOF.begin(), peakPosSortedBTOF.end());

  TF1* totalStudentFit = new TF1("totalStudentFit",
    "[0]*TMath::Student((x-[1])/[2],[3])/[2]"
    "+[4]*TMath::Student((x-[5])/[6],[7])/[6]"
    "+[8]*TMath::Student((x-[9])/[10],[11])/[10]",
    // ROUND 9: flat background term (round 7-8) removed -- see header note above
    zMinBTOF, zMaxBTOF);
  for(int s = 0; s < nSpecies; s++){
    double meanSeed = peaksBTOF[s].first;
    double heightSeed = peaksBTOF[s].second;
    double nearestSpacing = rangeBTOF;
    for(size_t j = 0; j < peakPosSortedBTOF.size(); j++){
      if(peakPosSortedBTOF[j] == meanSeed) continue;
      double d = TMath::Abs(peakPosSortedBTOF[j]-meanSeed);
      if(d < nearestSpacing) nearestSpacing = d;
    }
    double sigmaSeed = TMath::Max(2.0*binWidthBTOF, TMath::Min(0.25*nearestSpacing, 0.1*rangeBTOF));
    // ROUND 15: renamed from ndfSeed to nuSeed -- see header note. This is the
    // Student's-t distribution's continuous tail-shape parameter (conventionally
    // called nu/ndf in its own PDF formula), NOT a classical statistical "degrees of
    // freedom" count -- low nu = heavy tails, see header.
    double nuSeed = 3.0;
    // Student's-t peak density at x=mean is TMath::Student(0,nu)/sigma; invert that
    // (rather than the Gaussian 1/sqrt(2pi) constant) to get a consistent amp seed
    double peakDensityAtZero = TMath::Student(0.0, nuSeed);
    double ampSeed = heightSeed * sigmaSeed / peakDensityAtZero;
    totalStudentFit->SetParameter(4*s+0, ampSeed);
    totalStudentFit->SetParameter(4*s+1, meanSeed);
    totalStudentFit->SetParameter(4*s+2, sigmaSeed);
    totalStudentFit->SetParameter(4*s+3, nuSeed);
    totalStudentFit->SetParLimits(4*s+0, 0.0, 20.0*ampSeed + 1.0);
    totalStudentFit->SetParLimits(4*s+1, meanSeed - 0.5*nearestSpacing, meanSeed + 0.5*nearestSpacing);
    totalStudentFit->SetParLimits(4*s+2, 0.5*binWidthBTOF, TMath::Max(0.3*nearestSpacing, 3.0*binWidthBTOF));
    // ROUND 4 raised this floor from 1.0 to 2.0 as a stopgap for the combined-window
    // edge-blowup problem (rounds 3-6). ROUND 7-8 added and then tested a background
    // term for that; ROUND 9 removed it again (see header) since it wasn't clearly
    // resolving anything at single-mT-m0-bin statistics. This floor/ceiling stays as a
    // sanity bound on the tail-shape parameter (nu), not as background-avoidance.
    totalStudentFit->SetParLimits(4*s+3, 2.0, 30.0);
  }
  TFitResultPtr fitResultBTOF = hBTOF->Fit(totalStudentFit, "RQS0");
  fitResultBTOF = hBTOF->Fit(totalStudentFit, "RQS0"); // second pass from the first minimum
  cout << "BTOF fit status: " << fitResultBTOF->Status()
       << " (0 = converged cleanly; nonzero = check peak list / searchSigmaBinsBTOF above)" << endl;

  TF1* studentComponent[nSpecies];
  for(int s = 0; s < nSpecies; s++){
    studentComponent[s] = new TF1(Form("studentComponent%d", s),
      "[0]*TMath::Student((x-[1])/[2],[3])/[2]", zMinBTOF, zMaxBTOF);
    studentComponent[s]->SetParameter(0, totalStudentFit->GetParameter(4*s+0));
    studentComponent[s]->SetParameter(1, totalStudentFit->GetParameter(4*s+1));
    studentComponent[s]->SetParameter(2, totalStudentFit->GetParameter(4*s+2));
    studentComponent[s]->SetParameter(3, totalStudentFit->GetParameter(4*s+3));
    studentComponent[s]->SetNpx(300);
    studentComponent[s]->SetLineColor(speciesColors[s]);
    studentComponent[s]->SetLineWidth(2);
    // ROUND 16: peak height added alongside amp, same reasoning as the TPC panel --
    // peak = amp*TMath::Student(0,nu)/sigma for the Student's-t form.
    double sigmaLogBTOF = totalStudentFit->GetParameter(4*s+2);
    double nuLogBTOF    = totalStudentFit->GetParameter(4*s+3);
    double peakHeightLogBTOF = totalStudentFit->GetParameter(4*s+0) * TMath::Student(0.0, nuLogBTOF) / sigmaLogBTOF;
    cout << "  " << speciesNames[s] << " (BTOF): amp=" << totalStudentFit->GetParameter(4*s+0)
         << " mean=" << totalStudentFit->GetParameter(4*s+1)
         << " sigma=" << sigmaLogBTOF
         << " nu=" << nuLogBTOF
         << " peakHeight=" << peakHeightLogBTOF << endl;
  }
  totalStudentFit->SetLineColor(kBlack);
  totalStudentFit->SetLineStyle(2);
  totalStudentFit->SetLineWidth(2);
  totalStudentFit->SetNpx(500);

  //=========================================================================
  // Draw: two panels side by side, same 1600x800 canvas convention the rest of
  // the deck's dE/dx diagnostic plots use.
  //=========================================================================

  TCanvas* canv = new TCanvas("canv", "canv", 1600, 800);
  canv->Divide(2, 1);

  canv->cd(1);
  gPad->SetLogy();
  hTPC->SetMarkerStyle(20);
  hTPC->SetMarkerSize(0.6);
  hTPC->SetMarkerColor(kGray+2);
  hTPC->SetLineColor(kGray+2);
  hTPC->SetTitle(Form("TPC: Pion+ Z_{TPC}, Cent %s, y=[%1.3f,%1.3f], m_{T}-m_{0}=[%1.3f,%1.3f] GeV/c^{2}",
    CentralityLabel(a_centIndex).c_str(), rapLowTPC, rapHighTPC, mtLowValTPC, mtHighValTPC));
  hTPC->GetXaxis()->SetTitle("Z_{TPC}");
  hTPC->GetYaxis()->SetRangeUser(0.5, hTPC->GetMaximum()*3.0);
  hTPC->Draw("E1");
  totalGausFit->Draw("same");
  for(int s = 0; s < nSpecies; s++) gausComponent[s]->Draw("same");

  // Fitted parameters printed directly on the plot (not just the console log), per
  // Andrew's request (round 12). ROUND 14: round 13's wide box (x2=0.50) and one long
  // line per species actually overlapped the PEAK ITSELF, not just the legend --
  // ROOT's default pad margins put the frame at roughly NDC x=[0.1,0.9] mapped onto
  // the data's [-2,2] range, so NDC x=0.50 lands at data x~0.0, right where the pion
  // peak (mean~-0.03) sits; the long combined "amp=... mean=... sigma=..." line then
  // overflowed rightward into the legend on top of that. Fixed both problems: (1)
  // narrowed+shifted the box to NDC x=[0.13,0.40], which maps to roughly data
  // x=[-1.85,-0.4] -- comfortably left of every real peak seen so far in any round --
  // and (2) split each species onto separate one-value-per-line entries (name, A, mu,
  // sigma) instead of one long combined line, so no single line is ever wide enough to
  // overflow regardless of exact font metrics. Box height extended (down to y=0.36)
  // to fit the extra lines -- safe now that the box's x-range no longer overlaps the
  // peak, so extending it vertically only crosses the sparse, near-empty tail region
  // to its left, not the fitted curve. Species name lines are colored to match their
  // curve/legend color for a quick visual match.
  TPaveText* paramsTPC = new TPaveText(0.13, 0.36, 0.40, 0.90, "NDC");
  paramsTPC->SetBorderSize(0);
  paramsTPC->SetFillStyle(0);
  paramsTPC->SetTextAlign(12);
  paramsTPC->SetTextFont(42);
  paramsTPC->SetTextSize(0.022);
  paramsTPC->AddText("Fit parameters:");
  for(int s = 0; s < nSpecies; s++){
    TText* nameLineTPC = paramsTPC->AddText(speciesLegendLabels[s]);
    nameLineTPC->SetTextColor(speciesColors[s]);
    double ampValTPC   = totalGausFit->GetParameter(3*s+0);
    double sigmaValTPC = totalGausFit->GetParameter(3*s+2);
    // ROUND 16: A is the Gaussian's AREA-normalization constant (matching ZFitter's
    // own gausAmp() convention), not the visual peak height Andrew was comparing it
    // against -- see header note. Added the actual peak height (curve value at x=mean)
    // as its own line so both are on the plot and it's unambiguous which is which.
    double peakHeightTPC = ampValTPC * 0.3989422804 / sigmaValTPC;
    paramsTPC->AddText(Form("  A = %.4g", ampValTPC));
    paramsTPC->AddText(Form("  #mu = %.4g", totalGausFit->GetParameter(3*s+1)));
    paramsTPC->AddText(Form("  #sigma = %.4g", sigmaValTPC));
    paramsTPC->AddText(Form("  peak = %.4g", peakHeightTPC));
  }
  paramsTPC->Draw();

  TLegend* legTPC = new TLegend(0.58, 0.62, 0.88, 0.88);
  legTPC->SetBorderSize(0);
  legTPC->SetFillStyle(0);
  legTPC->AddEntry(hTPC, "Data", "lep");
  legTPC->AddEntry(totalGausFit, "Total fit", "l");
  for(int s = 0; s < nSpecies; s++) legTPC->AddEntry(gausComponent[s], speciesLegendLabels[s], "l");
  legTPC->Draw();

  canv->cd(2);
  gPad->SetLogy();
  hBTOF->SetMarkerStyle(20);
  hBTOF->SetMarkerSize(0.6);
  hBTOF->SetMarkerColor(kGray+2);
  hBTOF->SetLineColor(kGray+2);
  hBTOF->SetTitle(Form("BTOF: Pion+ Z_{BTOF}, Cent %s, y=[%1.3f,%1.3f], m_{T}-m_{0}=[%1.3f,%1.3f] GeV/c^{2}",
    CentralityLabel(a_centIndex).c_str(), rapLowBTOF, rapHighBTOF, mtLowValBTOF, mtHighValBTOF));
  hBTOF->GetXaxis()->SetTitle("Z_{BTOF}");
  hBTOF->GetYaxis()->SetRangeUser(0.5, hBTOF->GetMaximum()*3.0);
  hBTOF->Draw("E1");
  totalStudentFit->Draw("same");
  for(int s = 0; s < nSpecies; s++) studentComponent[s]->Draw("same");

  // Same fitted-parameters box as the TPC panel, including nu (BTOF's extra
  // Student's-t tail-shape parameter) since that's exactly the number that's been the
  // whole point of the last several rounds of iteration. ROUND 14: same narrower box +
  // one-value-per-line fix as the TPC panel above -- BTOF's peaks sit even closer to
  // the box's right edge in data terms (narrower [-0.5,0.5] axis), so the same NDC box
  // (x=[0.13,0.40] maps to roughly data x=[-0.46,-0.13] here) clears them with more
  // margin than on the TPC panel, not less. ROUND 15: label changed from "ndf" to the
  // TLatex Greek "#nu" -- see header note on why "ndf" is the wrong name for this
  // (it's the Student's-t distribution's own continuous shape parameter, not a
  // classical statistical degrees-of-freedom count, so it isn't restricted to
  // integers -- consistent with #mu/#sigma already used for the other parameters).
  TPaveText* paramsBTOF = new TPaveText(0.13, 0.30, 0.40, 0.90, "NDC");
  paramsBTOF->SetBorderSize(0);
  paramsBTOF->SetFillStyle(0);
  paramsBTOF->SetTextAlign(12);
  paramsBTOF->SetTextFont(42);
  paramsBTOF->SetTextSize(0.022);
  paramsBTOF->AddText("Fit parameters:");
  for(int s = 0; s < nSpecies; s++){
    TText* nameLineBTOF = paramsBTOF->AddText(speciesLegendLabels[s]);
    nameLineBTOF->SetTextColor(speciesColors[s]);
    double ampValBTOF   = totalStudentFit->GetParameter(4*s+0);
    double sigmaValBTOF = totalStudentFit->GetParameter(4*s+2);
    double nuValBTOF    = totalStudentFit->GetParameter(4*s+3);
    // ROUND 16: same A-vs-peak-height distinction as the TPC panel -- for the
    // Student's-t form, peak height = A*TMath::Student(0,nu)/sigma (TMath::Student(0,
    // nu) plays the same "density at the mode" role the Gaussian's 0.3989.../sigma
    // constant does, just nu-dependent instead of fixed).
    double peakHeightBTOF = ampValBTOF * TMath::Student(0.0, nuValBTOF) / sigmaValBTOF;
    paramsBTOF->AddText(Form("  A = %.4g", ampValBTOF));
    paramsBTOF->AddText(Form("  #mu = %.4g", totalStudentFit->GetParameter(4*s+1)));
    paramsBTOF->AddText(Form("  #sigma = %.4g", sigmaValBTOF));
    paramsBTOF->AddText(Form("  #nu = %.3g", nuValBTOF));
    paramsBTOF->AddText(Form("  peak = %.4g", peakHeightBTOF));
  }
  paramsBTOF->Draw();

  TLegend* legBTOF = new TLegend(0.58, 0.62, 0.88, 0.88);
  legBTOF->SetBorderSize(0);
  legBTOF->SetFillStyle(0);
  legBTOF->AddEntry(hBTOF, "Data", "lep");
  legBTOF->AddEntry(totalStudentFit, "Total fit", "l");
  for(int s = 0; s < nSpecies; s++) legBTOF->AddEntry(studentComponent[s], speciesLegendLabels[s], "l");
  legBTOF->Draw();

  canv->SaveAs(a_outPngPath);
  cout << "Wrote " << a_outPngPath << endl;
}
