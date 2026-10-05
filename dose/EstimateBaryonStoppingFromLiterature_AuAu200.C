// EstimateBaryonStoppingFromLiterature_AuAu200.C -- lightflavorspectra_OO200
//
// PURPOSE: a LITERATURE-BASED (not O+O-data-driven) estimate of the
// participant("stopped")/spectator("un-stopped") net-proton split at
// sqrt(s_NN)=200 GeV, using the SAME double/triple-Gaussian convention
// dose/ExtractNetProtonStoppingFraction_OO200.C defines for the eventual real
// O+O fit -- so this can be read as a physics-informed PRIOR for what that
// macro should find once fit_output.root exists, not as a substitute result.
//
// NAMING: deliberately NOT "_OO200" -- this is built entirely from PUBLISHED
// Au+Au measurements (a much heavier, thicker system than O-16+O-16), not
// from this analysis's own data. Do not treat its output as an O+O number.
//
// SOURCES (fetched and read this session -- see the discussion this macro
// followed for the fetch details):
//   - delta_y = 2.05 +/- 0.17 : rapidity loss (y_beam - <y_netproton>) at
//     sqrt(s_NN)=200 GeV Au+Au, BRAHMS, as compiled in F. Videbaek,
//     "Stopping and Baryon Transport in Heavy Ion Reactions"
//     (https://www4.rcf.bnl.gov/~videbaek/public/icpaqgp_videbaek.pdf).
//   - dN/dy(net-proton) ~ 7, flat for |y|<1, rising to ~12 by y~2-3, in
//     CENTRAL (0-5%) Au+Au at 200 GeV, still rising (no turnover) at the edge
//     of the measured window -- BRAHMS, "Nuclear Stopping in Au+Au Collisions
//     at sqrt(s_NN)=200 GeV," Phys. Rev. Lett. 93, 102301 (2004)
//     (abstract via arXiv:nucl-ex/0312023, https://arxiv.org/abs/nucl-ex/0312023).
//   - y_beam = ln(sqrt(s_NN)/m_p) = ln(200/0.938) = 5.36 -- same formula/value
//     already used in ExtractNetProtonStoppingFraction_OO200.C's own header.
//   - Au+Au net-proton sum rule = 2*Z_Au = 158 (79 protons/nucleus x 2 nuclei),
//     used ONLY as a consistency check on the crude reconstruction below, the
//     same role the analogous "==16" check plays for O+O elsewhere in dose/.
//
// EXTRACTION CAVEAT -- READ BEFORE TRUSTING THE NUMBERS: dN/dy(0)~7 and
// dN/dy(y~2.5)~12 were read off a conference-talk PDF and an arXiv abstract
// via automated fetch, NOT transcribed from the paper's own published data
// table or fit-parameter values (that full text wasn't retrievable from this
// sandbox). The double-Gaussian (y0, sigma, A) reconstructed below from just
// these two representative points is therefore a rough back-of-envelope
// stand-in for BRAHMS's own (unseen, in this session) fit -- the sum-rule
// check below is the honesty check on how rough: expect it to land somewhere
// around 80-90% of the target 158, not exact, consistent with "close, not
// exact" for a 2-point reconstruction of a real published shape. Refine
// NET_PROTON_DNDY_MIDRAP/NET_PROTON_DNDY_Y2P5 below against the actual PRL
// 93,102301 data table if better precision is ever needed.
//
// WHY THE STRICT (Andrew's chosen) CONVENTION GIVES participant_frac ~ 0:
// ExtractNetProtonStoppingFraction_OO200.C defines "participant/stopped" as
// the amplitude of a SEPARATE third Gaussian centered at y=0, on top of the
// mirror side-pair -- not "whatever fraction of a wide single Gaussian happens
// to land near y=0." The measured BRAHMS shape (flat then monotonically
// RISING toward higher |y|, no turnover/bump at y=0) is exactly that macro's
// own documented signature for preferring the double-Gaussian (side-pair-only)
// model over the triple. Under that same convention, applied honestly here,
// there is no separately-resolvable participant component at 200 GeV -- so
// participant_frac ~ 0, spectator_frac ~ 1. This is NOT this macro claiming to
// have fit a third Gaussian and found it exactly zero -- with only two
// representative published points there are too few constraints to even fit a
// 3-parameter (y0, sigma, A_center) triple-Gaussian (3 unknowns, 2 data
// points). It means: the two published numbers in hand do not REQUIRE a
// central component to explain them, and the paper's own qualitative
// description (flat-then-rising, "high degree of transparency") points the
// same way. A real quantitative upper LIMIT on a hidden central component
// would need the paper's actual chi2(double) vs chi2(triple) comparison,
// which this session could not retrieve.
//
// SYSTEM-SIZE CAVEAT: this is Au+Au (A=197), not O+O (A=16). Smaller/lighter
// systems are expected to show LESS stopping (more transparency) than Au+Au
// at the same energy -- thinner nuclear path length, fewer average NN
// re-scatterings per participant. So if this Au+Au-based participant_frac~0
// estimate is wrong in a specific direction for O+O, it is more likely an
// OVERESTIMATE of O+O's true participant/stopped fraction, not an
// underestimate -- treat this as a plausible UPPER BOUND going into O+O, not
// a central value.
//
// CENTRALITY CAVEAT: the two dN/dy reference numbers above are for CENTRAL
// (0-5%) Au+Au only -- this session could not retrieve BRAHMS's own
// centrality-differential net-proton dN/dy(y) values for the other bins. The
// literature does describe central vs. peripheral net-proton distributions
// differing (more peripheral generally trends toward LESS stopping, i.e. an
// even smaller participant fraction, not a larger one -- so applying the
// central-bin value to all bins is, if anything, another upper-bound-leaning
// choice, same direction as the Au+Au-vs-O+O caveat above).
//
// EXPLICIT ASSUMPTION (Andrew's instruction, applied 2026-08-21): stopping is
// treated as CENTRALITY-INDEPENDENT at this energy, so the 0-5% bin's
// participant_frac/spectator_frac split is now applied to all six bins below
// -- an assumption this macro applies because it was explicitly requested,
// not because it was derived or independently found in the literature. If
// real centrality-differential numbers turn up later, replace this uniform
// fill with them.
//
// USAGE:
//   root -l -b -q 'EstimateBaryonStoppingFromLiterature_AuAu200.C()'
//
// STATUS: UNTESTED -- no ROOT in the sandbox that wrote this (same caveat as
// every other C++ macro in this session). This one is pure arithmetic (no
// histograms, no fits, no external file I/O) so the syntax-risk surface is
// small, but the numbers it's built on carry the extraction caveat above.

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const int N_CENT_BINS = 6;
const char* CENT_LABELS[N_CENT_BINS] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};

//============================================================================
// Cited inputs -- see header for sources. CENTRAL (0-5%) Au+Au only.
//============================================================================
const double Y_BEAM_200GEV          = 5.36;  // ln(200/0.938), same value ExtractNetProtonStoppingFraction_OO200.C uses
const double DELTA_Y                = 2.05;  // BRAHMS/Videbaek, sqrt(s_NN)=200 GeV Au+Au
const double DELTA_Y_ERR            = 0.17;
const double NET_PROTON_DNDY_MIDRAP = 7.0;   // dN/dy at |y|<1, central 0-5% Au+Au (BRAHMS PRL 93,102301)
const double Y_REF2                 = 2.5;   // representative y in the quoted "y~2-3" range
const double NET_PROTON_DNDY_Y2P5   = 12.0;  // dN/dy at y~2.5, central 0-5% Au+Au (same source)
const double AU_SUM_RULE_TARGET     = 158.0; // 2*79, Au+Au net-proton sum rule -- consistency check only

//============================================================================
// Closed-form solve for a mirror-symmetric double Gaussian,
//   dNdy(y) = A*[exp(-(y-y0)^2/2sigma^2) + exp(-(y+y0)^2/2sigma^2)],
// given y0 = y_beam - delta_y and the two reference dN/dy points above.
// Derivation: at y=0, dNdy(0)=2*A*exp(-y0^2/2sigma^2). At y=Y_REF2 (with
// Y_REF2 << y0 so the far mirror term is negligible), dNdy(Y_REF2) is
// approx A*exp(-(Y_REF2-y0)^2/2sigma^2). Dividing these two conditions and
// solving for sigma^2 gives the closed form below -- computed here, not
// hardcoded, so it updates automatically if the cited inputs are refined.
//============================================================================
void SolveDoubleGaussian(double deltaY, double& y0Out, double& sigmaOut, double& aOut, double& sumRuleOut){
  double y0 = Y_BEAM_200GEV - deltaY;
  double dTerm = Y_REF2 * (2.0 * y0 - Y_REF2); // = y0^2 - (Y_REF2-y0)^2
  double lnTerm = log(2.0 * NET_PROTON_DNDY_Y2P5 / NET_PROTON_DNDY_MIDRAP);
  double sigma2 = dTerm / (2.0 * lnTerm);
  double sigma = sqrt(sigma2);
  double a = (NET_PROTON_DNDY_MIDRAP / 2.0) * exp(y0 * y0 / (2.0 * sigma2));
  double sumRule = 2.0 * a * sigma * sqrt(2.0 * M_PI); // integral of both side-Gaussians, -inf..inf

  y0Out = y0; sigmaOut = sigma; aOut = a; sumRuleOut = sumRule;
}

void EstimateBaryonStoppingFromLiterature_AuAu200(){

  cout << "==================================================================" << endl;
  cout << "EstimateBaryonStoppingFromLiterature_AuAu200" << endl;
  cout << "Literature-based double-Gaussian reconstruction, Au+Au sqrt(s_NN)=200 GeV" << endl;
  cout << "NOT an O+O result, NOT a fit to this analysis's own data -- see header." << endl;
  cout << "==================================================================" << endl;
  cout << "Cited inputs: y_beam=" << Y_BEAM_200GEV << "  delta_y=" << DELTA_Y << " +/- " << DELTA_Y_ERR
       << "  dN/dy(|y|<1)=" << NET_PROTON_DNDY_MIDRAP << "  dN/dy(y~" << Y_REF2 << ")=" << NET_PROTON_DNDY_Y2P5 << endl;

  //--------------------------------------------------------------------------
  // Central value plus a crude +/-1sigma envelope from delta_y's own quoted
  // uncertainty (simple symmetric shift, not full error propagation -- an
  // order-of-magnitude sanity envelope, not a rigorous uncertainty).
  //--------------------------------------------------------------------------
  double y0, sigma, a, sumRule;
  SolveDoubleGaussian(DELTA_Y, y0, sigma, a, sumRule);

  double y0Lo, sigmaLo, aLo, sumRuleLo;
  SolveDoubleGaussian(DELTA_Y + DELTA_Y_ERR, y0Lo, sigmaLo, aLo, sumRuleLo); // larger delta_y -> smaller y0

  double y0Hi, sigmaHi, aHi, sumRuleHi;
  SolveDoubleGaussian(DELTA_Y - DELTA_Y_ERR, y0Hi, sigmaHi, aHi, sumRuleHi); // smaller delta_y -> larger y0

  cout << endl << "Reconstructed double-Gaussian (central value):" << endl;
  cout << "  y0 (side-Gaussian centroid)    = " << y0    << "   [range " << y0Lo    << " - " << y0Hi    << " from delta_y +/- " << DELTA_Y_ERR << "]" << endl;
  cout << "  sigma (side-Gaussian width)    = " << sigma << "   [range " << sigmaLo << " - " << sigmaHi << "]" << endl;
  cout << "  A (per-side amplitude)         = " << a     << "   [range " << aLo     << " - " << aHi     << "]" << endl;

  cout << endl << "Sum-rule consistency check (Au+Au target = " << AU_SUM_RULE_TARGET << " net protons):" << endl;
  cout << "  Integrated double-Gaussian     = " << sumRule
       << "  (" << 100.0 * sumRule / AU_SUM_RULE_TARGET << "% of target)" << endl;
  cout << "  Range from delta_y uncertainty = " << sumRuleLo << " - " << sumRuleHi
       << "  (" << 100.0*sumRuleLo/AU_SUM_RULE_TARGET << "% - " << 100.0*sumRuleHi/AU_SUM_RULE_TARGET << "% of target)" << endl;
  cout << "  A large deviation from 100% here is expected -- this is a 2-point" << endl;
  cout << "  back-of-envelope reconstruction, not BRAHMS's own published fit. See header." << endl;

  //--------------------------------------------------------------------------
  // Participant/spectator split under the STRICT convention (Andrew's chosen
  // option): participant = amplitude of a SEPARATE central Gaussian, which
  // this double-Gaussian-only model does not have by construction. See header
  // for why this is a real reading of the data, not an artifact of this
  // macro's own choices.
  //--------------------------------------------------------------------------
  cout << endl << "==================================================================" << endl;
  cout << "Participant/spectator split (STRICT convention, matches" << endl;
  cout << "ExtractNetProtonStoppingFraction_OO200.C's own definition):" << endl;
  cout << "==================================================================" << endl;
  const double PARTICIPANT_FRAC_0TO5 = 0.0;
  const double SPECTATOR_FRAC_0TO5   = 1.0;

  printf("%-10s%18s%18s%12s\n", "Cent", "participant_frac", "spectator_frac", "source");
  printf("%-10s%18.3f%18.3f%12s\n", CENT_LABELS[0], PARTICIPANT_FRAC_0TO5, SPECTATOR_FRAC_0TO5, "BRAHMS");
  for(int i = 1; i < N_CENT_BINS; i++){
    printf("%-10s%18.3f%18.3f%12s\n", CENT_LABELS[i], PARTICIPANT_FRAC_0TO5, SPECTATOR_FRAC_0TO5, "assumed");
  }
  cout << endl << "Only 0-5% has a direct literature basis (BRAHMS central Au+Au data). The other" << endl;
  cout << "five bins carry the SAME split under Andrew's explicit assumption that stopping" << endl;
  cout << "is centrality-independent at this energy -- not an independent literature result" << endl;
  cout << "for those bins (see header). Since real stopping is expected to WEAKEN, not" << endl;
  cout << "strengthen, toward peripheral collisions, this uniform fill likely still reads as" << endl;
  cout << "an upper bound on participant_frac for bins beyond 0-5%, same direction as the" << endl;
  cout << "Au+Au-vs-O+O caveat below." << endl;
  cout << endl << "Reminder: this is Au+Au, a much heavier/thicker system than O+O -- read as a" << endl;
  cout << "plausible UPPER BOUND on O+O's participant fraction at 200 GeV, not a direct" << endl;
  cout << "estimate of it, until the real O+O fit_output.root-based extraction exists." << endl;
}
