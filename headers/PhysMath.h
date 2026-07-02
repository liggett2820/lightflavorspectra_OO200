#ifndef _PHYS_MATH_
#define _PHYS_MATH_

#include "../makefile_toggles.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Math/Boost.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/LorentzVector.h"
#include "Math/GenVector/Boost.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TSystem.h"
#include "Math/IFunction.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/LorentzVector.h"
#include "TGraphAsymmErrors.h"

#include "Math/DisplacementVector3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#ifdef _MATHMORE_
  #include "Math/SpecFuncMathMore.h"
#endif
#include "TVirtualFitter.h"
#ifdef _BOOST_
  #include <boost/math/distributions/skew_normal.hpp>
  #include <boost/math/special_functions/owens_t.hpp>
#endif

#include "TCanvas.h"
#include <iostream>
#include <complex>

//#ifndef _NO_HELIX_
//#include "../headers/Helix.h"
//#endif

//#include "../headers/HelixPhysical.h"
//#include "../headers/HistogramUtilities.h"


using namespace std;

//typedef double Scalar;
//typedef DisplacementVector3D< Cartesian3D<Scalar> > BetaVector;


const double physMath_Li2_P[8] = {
      1.0706105563309304277e+0,
     -4.5353562730201404017e+0,
      7.4819657596286408905e+0,
     -6.0516124315132409155e+0,
      2.4733515209909815443e+0,
     -4.6937565143754629578e-1,
      3.1608910440687221695e-2,
     -2.4630612614645039828e-4
  };


const double physMath_Li2_Q[8] = {
      1.0000000000000000000e+0,
     -4.5355682121856044935e+0,
      8.1790029773247428573e+0,
     -7.4634190853767468810e+0,
      3.6245392503925290187e+0,
     -8.9936784740041174897e-1,
      9.8554565816757007266e-2,
     -3.2116618742475189569e-3
   };

const double physMath_sqrtTwo = sqrt(2.0);


namespace PhysMath{
  // a_m is mass   |   a_p is total momentum   |  a_y is rapidity
  double dmtde (double a_y);
  double dthetadomega2pi(double a_theta);
  double dydtheta(double a_m, double a_pt, double a_theta);
  double invBeta(double a_m, double a_p);
  double invBeta_from_MassSquaredSigned(double a_signed_m2, double a_p);
  double invBetaWidth_from_MassSquaredWidth(double a_signed_m2, double a_width, double a_p);
  double beta_from_MassSquaredSigned(double a_signed_m2, double a_p);
  double rapidity(double a_m,double a_pT,double a_pZ);
  void   rapidity_mTm0(ROOT::Math::PxPyPzEVector* a_vect, double &a_y, double &a_mTm0);
  double rapFromEtaPt(double a_m, double a_eta, double a_pt);
  double rapCM(double y_lab, double yCM);
  double mass_fromBeta(double a_p, double a_beta);//always positive
  double massSquared_signed_fromBeta(double a_p, double a_beta); // negative for beta > 1
  double mass_signed_fromBeta(double a_p, double a_beta);
  double massSquared_signed_fromInvBeta(double a_p, double a_invbeta); // negative for invbeta < 1
  double massSquaredWidth_from_invBetaWidth(double a_p, double a_invBeta, double a_width);
  double mTm0(double a_m, double a_pt);
  double mTm0_from_y_p(double a_y, double a_p, double a_m);// mT-m0 from rapidity(y) and momentum(p) and mass (m)
  double mT_from_mTm0_m(double a_mTm0, double a_m);
  double pT(double a_px, double a_py);
  double pT_from_mTm0(double a_mTm0, double a_m);
  double pZ_from_y_mTm0(double a_y, double a_mTm0, double a_m);
  double pZ_from_eta_pT(double a_eta, double a_pT);
  double pTot(double a_px, double a_py, double a_pz);
  double pTot_from_y_mTm0(double a_y, double a_mTm0, double a_m);
  double pTot_from_eta_pT(double a_eta, double a_pT);
  double eta_from_p(double a_px, double a_py, double a_pz);
  double eta_from_y_mTm0(double a_y, double a_mTm0, double a_m);
  double theta_from_eta(double a_eta); //note theta = pi/2 when eta = 0
  double eta_from_theta(double a_theta);
  double theta_from_y_mTm0(double a_y, double a_mTm0, double a_m);
  double invBeta_width_prediction(double a_timeRes,double a_y, double a_mTm0, double a_m, double a_bTOFRadius = 2.0);
  double m2_width_prediction(double a_timeRes,double a_y, double a_mTm0, double a_m, double a_bTOFRadius = 2.0);
  double m2_width_prediction(double a_InvBetaRes,double a_p, double a_m);


  bool tofSelectionCut(double a_p,double a_m_0,double a_inverseBeta,double a_cut = 0.03 );


  vector<double> calcMassSquaredEdges(int a_numBins, double a_timingResolution, double a_momToSample, double a_lowestParticleMass, double a_highestParticleMass);


  //this is the converstion factor between dN/dy*dpT = fact * dN/dn*dpT
  double eta_to_rap_factor(double a_m, double a_eta, double a_pt);

  // converstion factor from dN/dy*dpT = fact
  double pt_to_mT_factor(double a_m, double a_pt);

  double beta_to_betagamma(double a_beta);
  double betagamma(double a_p, double a_m);
  double gamma(double a_p, double a_m);
  double beta(double a_p, double a_m);
  double invBetaGammaSqr_signed_fromInvBeta(double a_invBeta);

  double integralWithMtScale(TF1* a_funct, double a_mass, double a_stepSize = 0.001, double a_lowbound = 0.0, double a_highbound = 35.0, bool a_extratrapolateToInf = true);
  double interpolate(double a_x, double a_x1, double a_y1, double a_x2,double a_y2);
  void   interpolate(double a_x, double a_x1, double a_y1, double a_yErr1, double a_x2, double a_y2, double a_yErr2, double &a_yValue, double &a_yError);
  void   interpolateAsymmError(double a_x, double a_x1, double a_y1, double a_yErrLow1, double a_yErrHigh1, double a_x2, double a_y2, double a_yErrLow2, double a_yErrHigh2, double &a_yValue, double &a_yErrLow, double &a_yErrHigh);

  // Interpolation, x2 > x1 and y2 > y1 and x,y must be inside. Z1 is bottom left  Z2 is bottom right, Z3 is top left, Z4 if top right in xy plane
  double interpolateSquare2D(double a_x, double a_y, double a_x1, double a_y1, double a_x2,double a_y2, double a_z1, double a_z2, double a_z3, double a_z4);
  double signed_sqrt(double a_number);
  double signed_square(double a_number);
  double addInQuadriture(double a_num1, double a_num2);

  // These calculate the remaining fraction of particles with mean lifetime a_tau. Use the third one if you want it to have a seam between bTof and eTOF
  double fractionRemainingAfterHorizontalDistance(double a_rapidity, double a_mTm0, double a_mass, double a_tau, double a_horizontalDistance); // exponentials are crazy, this is both the fraction after/total and the prob at that distance
  double fractionRemainingAfterRadialDistance(    double a_rapidity, double a_mTm0, double a_mass, double a_tau, double a_radialDistance);
  double fractionRemainingAfterDistance(          double a_rapidity, double a_mTm0, double a_mass, double a_tau, double a_radialDistance, double a_horizontalDistance);

  double fractionRemainingAfterHorizontalDistance_EtaPt(double a_eta, double a_pT, double a_mass, double a_tau, double a_horizontalDistance);
  double fractionRemainingAfterRadialDistance_EtaPt(    double a_eta, double a_pT, double a_mass, double a_tau, double a_radialDistance);
  double fractionRemainingAfterDistance_EtaPt(          double a_eta, double a_pT, double a_mass, double a_tau, double a_radialDistance, double a_horizontalDistance);

  double fractionRemainingAfterDistance(double a_pTot, double a_mass, double a_c_tau, double a_distance);


  double calc_pct_kaon_left(double a_rapidity, double a_mTm0, double a_mass);


  //Decays Particle in Rest Frame, with a random 3D angle. Must setRandGenSeed before using this function // the lorentz vector needs to be px py pz E
  void centerOfMassTwoBodyDecay(TRandom3* a_randObj, double a_parent_mass, double a_first_daughter_mass, double a_second_daughter_mass, TLorentzVector* &a_first_daughter_vect , TLorentzVector* &a_second_daughter_vect); // the vectors passed are filled
  void centerOfMassTwoBodyDecay(TRandom3* a_randObj, double a_parent_mass, double a_first_daughter_mass, double a_second_daughter_mass, ROOT::Math::PxPyPzEVector* &a_first_daughter_vect , ROOT::Math::PxPyPzEVector* &a_second_daughter_vect); // the vectors passed are filled

  void boostedTwoBodyDecay_Vect(TRandom3* a_randObj, double a_parent_mass, TLorentzVector* &a_parent_vect, double a_first_daughter_mass, double a_second_daughter_mass,   TLorentzVector* &a_first_daughter_vect , TLorentzVector* &a_second_daughter_vect);
  void boostedTwoBodyDecay(TRandom3* a_randObj, double a_parent_mass, double a_parent_rapidity, double a_parent_mTm0,
                           double a_first_daughter_mass, double &a_first_daughter_rapidity, double &a_first_daughter_mTm0,
                           double a_second_daughter_mass, double &a_second_daughter_rapidity, double &a_second_daughter_mTm0);

  void boostedTwoBodyDecay(ROOT::Math::Boost* a_boostObj,TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex, double a_parent_mass, ROOT::Math::PxPyPzEVector* a_parent_vect,
                                double a_first_daughter_mass,  ROOT::Math::PxPyPzEVector* &a_first_daughter_vect,
                                double a_second_daughter_mass, ROOT::Math::PxPyPzEVector* &a_second_daughter_vect,
                                ROOT::Math::XYZVector* &a_SecondaryVertex, double a_parent_MeanLifetime);

  // These compute a DCA With an approximation of the daughter vector going in a straight line
  double dca_from_decay(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect);
  void dca_from_decay(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect, double &a_dca_x, double &a_dca_y, double &a_dca_z);
  // Charges should be -1 or 1, distance vectors should be in cm, dca returned in cm
  //double dca_from_decay_helix(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect, int a_daughter_charge);
  //void dca_from_decay_helix(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect, int a_daughter_charge, double &a_dca_x, double &a_dca_y, double &a_dca_z);



  void boostedThreeBodyDecay_Vect(double a_parent_mass,TLorentzVector* &a_parent_vect, double a_first_daughter_mass, double a_second_daughter_mass,double a_third_daughter_mass, TLorentzVector* &a_first_daughter_vect , TLorentzVector* &a_second_daughter_vect,TLorentzVector* &a_third_daughter_vect); // the vectors passed are filled
  void boostedThreeBodyDecay(TRandom3* a_randObj, double a_parent_mass, double a_parent_rapidity, double a_parent_mTm0, double a_first_daughter_mass, double &a_first_daughter_rapidity, double &a_first_daughter_mTm0, double a_second_daughter_mass, double &a_second_daughter_rapidity, double &a_second_daughter_mTm0, double a_third_daughter_mass, double &a_third_daughter_rapidity, double &a_third_daughter_mTm0);
  void boostedThreeBodyDecay(TGenPhaseSpace* a_generator, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                                double a_parent_mass,          ROOT::Math::PxPyPzEVector* a_parent_vect,
                                double a_first_daughter_mass,  ROOT::Math::PxPyPzEVector* &a_first_daughter_vect,
                                double a_second_daughter_mass, ROOT::Math::PxPyPzEVector* &a_second_daughter_vect,
                                double a_third_daughter_mass,  ROOT::Math::PxPyPzEVector* &a_third_daughter_vect,
                                ROOT::Math::XYZVector* a_secondaryVertex, double a_parent_MeanLifetime);
  double decayLength(TRandom3* a_randObj,double a_momentum, double a_mass, double a_meanLifetimeInSeconds); // returns random length in centimeters

  //specific decays  // OLD METHOD Without DCA Calculation
  void lambdaDecay(TRandom3* a_randObj, double a_lambda_rapidity,double a_lambda_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pion_rapidity, double &a_pion_mTm0); // Lambda -> p + piMinus
  void sigmaPlusDecay(TRandom3* a_randObj, double a_sigma_rapidity,double a_sigma_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pion_rapidity, double &a_pion_mTm0); // SigmaPlus -> p + pi0
  void xiMinusDecay(TRandom3* a_randObj, double a_xi_rapidity,double a_xi_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pion_rapidity, double &a_pion_mTm0, double &a_second_pion_rapidity, double &a_second_pion_mTm0); // XiMinus -> p + piMinus + piMinus
  void xiZeroDecay(TRandom3* a_randObj, double a_xi_rapidity,double a_xi_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pionZero_rapidity, double &a_pionZero_mTm0, double &a_pionMinus_rapidity, double &a_pionMinus_mTm0); // XiZero -> p + pi0 + piMinus
  //Branching ratios built in to omega Minus decay
  void omegaMinusDecay(TRandom3* a_randObj, double a_omega_rapidity,double a_omega_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_kaon_rapidity, double &a_kaon_mTm0, double &a_pionMinus_rapidity, double &a_pionMinus_mTm0, double &a_second_pion_rapidity, double &a_second_pion_mTm0);
  void sigmaZeroDecay(TRandom3* a_randObj, double a_sigma_rapidity,double a_sigma_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pion_rapidity, double &a_pion_mTm0); // SigmaZero -> p + pi- + gamma
  void kaonZeroShortDecay(TRandom3* a_randObj, double a_kaon_rapidity,double a_kaon_mTm0, double &a_piPlus_rapidity, double &a_piPlus_mTm0, double &a_piMinus_rapidity, double &a_piMinus_mTm0); // K0S -> pi+ + pi-
  void etaDecay(TRandom3* a_randObj, double a_eta_rapidity,double a_eta_mTm0, double &a_piPlus_rapidity, double &a_piPlus_mTm0, double &a_piMinus_rapidity, double &a_piMinus_mTm0);

  void lambdaDecay(    ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_lambda_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pion_vect, ROOT::Math::XYZVector* &a_SecondaryVertex); // Lambda -> p + piMinus
  void sigmaPlusDecay( ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_sigma_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pion_vect, ROOT::Math::XYZVector* &a_SecondaryVertex); // SigmaPlus -> p + pi0
  void xiMinusDecay(    ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_xi_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pion_vect, ROOT::Math::PxPyPzEVector* &a_second_pion_vect, ROOT::Math::XYZVector* &a_SecondaryVertex, ROOT::Math::XYZVector* &a_ThirdVertex); // XiMinus -> p + piMinus + piMinus
  void xiZeroDecay(  ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_xi_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pionZero_vect, ROOT::Math::PxPyPzEVector* &a_pionMinus_vect, ROOT::Math::XYZVector* &a_SecondaryVertex, ROOT::Math::XYZVector* &a_ThirdVertex); // XiZero -> p + pi0 + piMinus
  //Branching ratios built in to omega Minus decay
  void omegaMinusDecay( ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
          ROOT::Math::PxPyPzEVector* a_omega_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pionZero_vect,
          ROOT::Math::PxPyPzEVector* &a_kaon_vect, ROOT::Math::PxPyPzEVector* &a_pionMinus_vect, ROOT::Math::PxPyPzEVector* &a_second_pion_vect,
          ROOT::Math::XYZVector* &a_SecondaryVertex, ROOT::Math::XYZVector* &a_ThirdVertex, ROOT::Math::XYZVector* &a_FourthVertex,
          int &a_proton_vertex_index, int &a_kaon_vertex_index, int &a_pion_vertex_index, int &a_second_pion_vertex_index);
  void sigmaZeroDecay( ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_sigma_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pion_vect, ROOT::Math::XYZVector* &a_SecondaryVertex, ROOT::Math::XYZVector* &a_ThirdVertex); // SigmaZero -> p + pi- + gamma
  void kaonZeroShortDecay(ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_kaon_vect, ROOT::Math::PxPyPzEVector* &a_piPlus_vect, ROOT::Math::PxPyPzEVector* &a_piMinus_vect, ROOT::Math::XYZVector* &a_SecondaryVertex); // K0S -> pi+ + pi-
  void etaDecay(     TGenPhaseSpace* a_generator, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_eta_vect, ROOT::Math::PxPyPzEVector* &a_piPlus_vect, ROOT::Math::PxPyPzEVector* &a_piMinus_vect, ROOT::Math::XYZVector* &a_SecondaryVertex);



  //##################   MATH FUNCTIONS       ######################
  complex< double > ComplexGamma( complex<double> z);
  double  Li2(double a_x); //Dilogarithm x=[0,1]
  double  LogGammaR2(double a_x, double a_y); // returns log(abs(gamma(x+iy)/gamma(x))^2)
  double* choleskyDecomp(double *A, int n, bool a_print);
  double  crossProduct3D(double* a_vect_1,double* a_vect_2);
  void    crossProduct3D(double* a_vect_1,double* a_vect_2,double &a_cross_x, double &a_cross_y, double &a_cross_z);

  bool        isTMatrixDSymEmpty(TMatrixDSym a_matrix);
  TMatrixDSym reduceCovMatrix(TMatrixDSym a_covMatrix);
  TMatrixD    reduceCovMatrix(TMatrixD a_covMatrix);
  TMatrixDSym reduceCovMatrix_keepPassedVector(TMatrixDSym a_covMatrix, vector< int > a_rowsToKeep); //a_rowsToKeep should contain indecies of rows from covariance matrix
  TMatrixDSym addZerosToCovMatrix(TMatrixDSym a_covMatrix, int a_numExtraZeroRows); // adds extra rows and coloumns
  TH2D* makeCorrelationHistogram(string a_name, string a_title, TMatrixDSym a_covMatrix, vector< string > a_paramNames = {});
  TH2D* makeCovarianceHistogram(string a_name, string a_title, TMatrixDSym a_covMatrix, vector< string > a_paramNames = {});


  // when you take the ln of a histogram, you can fit with parabolas, this function changes those fit parameters into gaussian parameters
  // returns {N,mu,sig,N_error,mu_error,sig_error}
  vector< double > logParabToGaus(double a_scale, double a_mean,
                                  double a_lambda, double a_scaleError = 0.0,
                                  double a_meanError = 0.0, double a_lambdaError = 0.0);
  //returns {scale, mean, lambda}
  vector< double > logGausToParab(double a_N, double a_mu, double a_sig);
  double studentT_var_to_stdDev(double a_sigma_studentT, double a_nu_studentT);
  double studentT_var_from_stdDev(double a_stdDev_gaus, double a_nu_studentT);

  double gaus(double a_mean, double a_sig, double a_x);
  double gaus_by_deriviatives(double a_x,double a_yLoc,double a_ySlope,double a_ySecDer,double a_xLoc);
  double gaus_finite_integral(double a_N, double a_mu, double a_sig, double a_x_low, double a_x_high);
  void   get3PointGausParams(double &a_N, double &a_mu, double &a_sig, double a_x1,double a_y1,double a_x2,double a_y2,double a_x3,double a_y3);
  void   get3PointParabParams(double a_x0, double &a_y0, double &a_ySlope, double &a_ySecDer, double a_x1,double a_y1,double a_x2,double a_y2,double a_x3,double a_y3);
  void   gausParams_from_deriviatives(double a_yLoc,double a_ySlope,double a_ySecDer,double a_xLoc, double &a_N, double &a_mu, double &a_sig);
  void   derivatives_from_gausParams(double &a_yLoc,double &a_ySlope,double &a_ySecDer,double a_xLoc, double a_N, double a_mu, double a_sig);
  double lineByDerivatives(double a_x, double a_yLoc, double a_ySlope,double a_xLoc);
  double parabolaByDerivatives(double a_x, double a_yLoc, double a_ySlope, double a_ySecDer,double a_xLoc);

  double logNorm(double a_mu, double a_sig, double a_x);
  void   gausToLogNorm(double a_mean, double a_stdDev, double &a_mu, double &a_sig);
  void   logNormToGaus(double a_mu, double a_sig, double &a_mean, double &a_stdDev);
  double student_T(double a_mu, double a_sig, double a_nu, double a_x);
  double student_T_BinContent(double a_mu, double a_sig, double a_nu, double a_x_low, double a_x_high); //already divided by bin width (just 10 sum samplings)
  #ifdef _MATHMORE_
    double student_T_finite_integral(double a_N,double a_mu, double a_sig, double a_nu, double a_x_low, double a_x_high);
  #endif

  double m2_student_T(double a_mu, double a_sig, double a_nu, double a_x, double a_momentum); // normalized in invBeta space
  double m2_student_T_BinContent(double a_mu, double a_sig, double a_nu, double a_x_low, double a_x_high, double a_momentum);



  double skewNormal(double a_mu, double a_sig, double a_skew, double a_x);
  double skewNormalBinContent(double a_mu, double a_sig, double a_skew, double a_xLow,double a_xHigh); //already divided by bin width (just 10 sum samplings)
  #ifdef _BOOST_
    //alpha is a_skew/sqrt(1.0 - a_skew*a_skew) for speed reasons
    double skewNormal_finite_integral(double a_N,double a_mu, double a_sig, double a_alpha, double a_x_low, double a_x_high);
  #endif
  double generalNormalV2(double a_mode, double a_sigma, double a_skew, double a_x);

  double pearsonIV(double a_lambda, //location parameter [anything]
                  double a_alpha, //scale parameter (stdDev?)  [> 0]
                  double a_m, // skew    [> 1/2]
                  double a_nu, // tailness  [anything (0 is type VII)]
                  double a_x, // where to evaluate the function
                  bool a_normalize = true // normalization constant applied (false will still scale so f(mode)=1.0)
                );
  double pearsonIVOperator(double* a_x,const double* a_params);  // lambda ~ mode     alpha ~ stdDev      m ~ kurtosis       nu ~ skewness   // none of these are exact, but just main contributions
  double pearsonIVAmpToCounts(double a_amplitude,double a_lambda,double a_alpha,double a_m,double a_nu);
  double piecewiseTurnOn(double a_A, double a_B, double a_C, double a_xCut, double a_ySlopeCut, double a_m, double a_x);
  double piecewiseTurnOnTwoLines(double a_A, double a_B, double a_C,double a_x0_cut, double a_x1_cut, double a_y1_cut,  double a_ySlope1, double a_m, double a_x);

  TF1* getPearsonAndTwoGausFunct(string a_functName, double a_N = 1.0, double a_lambda = 0, double a_alpha = 1.0, double a_m = 1.0, double a_nu = 0.0, // PearsonIV
                                 double a_N_g1 = 1.0, double a_mu_g1 = 0.0, double a_sig_g1 = 1.0, //Gaussian 1
                                 double a_N_g2 = 1.0, double a_mu_g2 = 0.0, double a_sig_g2 = 1.0, bool a_normalizePearson = false); // Gaussian 2
  TF1* getPearsonAndTwoGausFunct_sigmaLinked(string a_functName, double a_N , double a_lambda, double a_alpha , double a_m , double a_nu ,
                                 double a_N_g1 , double a_mu_g1 , double a_sig_g1,
                                 double a_N_g2 , double a_mu_g2 ,bool a_normalizePearson);
  TF1* getPearsonAndTwoGausFunct_ratioLinked(string a_functName, double a_N , double a_lambda, double a_alpha , double a_m , double a_nu ,
                                 double a_fixedRatio , double a_mu_g1 , double a_sig_g1,
                                 double a_N_g2 , double a_mu_g2 ,double a_sig_g2, bool a_normalizePearson);

  TF1* getSkewAndTwoGausFunct(string a_functName, double a_N , double a_mu, double a_sig , double a_skew,
                                         double a_N_g1 , double a_mu_g1 , double a_sig_g1,
                                         double a_N_g2 , double a_mu_g2 , double a_sig_g2);
  TF1* getSkewAndTwoGausFunct_sigmaLinked(string a_functName, double a_N , double a_mu, double a_sig , double a_skew,
                                         double a_N_g1 , double a_mu_g1 , double a_sig_g1,
                                         double a_N_g2 , double a_mu_g2);
  TF1* getSkewAndTwoGausFunct_ratioLinked(string a_functName, double a_N , double a_mu, double a_sig , double a_skew,
                                 double a_N2ToN1or3Ratio, double a_mu_g1 , double a_sig_g1,
                                 double a_N_g2 , double a_mu_g2 ,double a_sig_g2, bool a_Ratio12InsteadOf23 = false);
  TF1* getTwoSkewFunct(string a_functName, double a_N1 , double a_mu1, double a_sig1 , double a_skew1, double a_N2, double a_mu2, double a_sig2, double a_skew2);
  TF1* getTwoSkewFunct_SafeParam4Amp(string a_functName, double a_N1 , double a_mu1, double a_sig1 , double a_skew1, double a_N2, double a_mu2, double a_sig2, double a_skew2);
  //These are for the simulatneous ZFitter Fitting (see there for parameter meanings)
  TF1* getThreeGausFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3);
  TF1* getThreeSkewFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3);
  TF1* getThreeSkewFunct(string a_functName, double a_N1, double a_mu1, double a_sig1, double a_gam1, double a_N2,  double a_mu2, double a_sig2, double a_gam2, double a_N3,  double a_mu3, double a_sig3, double a_gam3);
  TF1* getThreeStudentTFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3);
  TF1* getThreeM2StudentTFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3, double a_momentum);
  TF1* getThreePearsonIVFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3);
  TF1* getThreeDoubleGaussFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3);


  double blastWaveModelFit(double *x, double *par);
  double pTblastWaveModelFit(double *x, double *par);
  double legendre(int a_order,double a_x); // up to order 29
  double semiHermite(int a_order, double a_x);
  double tsallis(double a_dNdy,double a_T,double a_q,double a_m, double a_mTm0); //mTm0 is the variable // see https://indico.cern.ch/event/232225/contributions/1541626/attachments/386231/537252/tsallis_cern.pdf




  //Functions
  TF1* getSimpleExponentialFunction(string a_name, double a_scale, double a_exponent);
  TF1* getExponentialFunction(string a_name, double a_scale, double a_exponent, double a_leftRightShift, double a_y0);
  TF1* getFitSafeExponentialFunction(string a_name, double a_scale, double a_exponent);
  TF1* getFitSafeDoubleExponentialFunction(string a_name, double a_scale, double a_exponent, double a_scale_2, double a_exponent_2);
  TF1* getFitSafeTripleExponentialFunction(string a_name, double a_scale, double a_exponent, double a_scale_2, double a_exponent_2,double a_scale_3, double a_exponent_3);
  TF1* getExpoPoly4Function(string a_name,double a_scale, double a_p0, double a_p1, double a_p2, double a_p3, double a_p4);


  TF1* getStudentTFunction(string a_name, double a_N, double a_mu, double a_sig, double a_nu );
  TF1* getM2StudentTFunction(string a_name, double a_N, double a_mu, double a_sig, double a_nu, double a_momentum);
  TF1* getGaussianFunction(string a_name,double a_N, double a_mu, double a_sig);
  TF1* getGaussianFunction_byDerivatives_GausPars(string a_name,double a_N, double a_mu, double a_sig, double a_xLoc);
  TF1* getGaussianFunction_byDerivatives(string a_name,double a_yLoc, double a_ySlope, double a_ySecDer, double a_xLoc);
  TF1* getTwoGaussianFunction(string a_name,double a_N1, double a_mu1, double a_sig1,double a_N2, double a_mu2, double a_sig2);
  TF1* getGaussianConvolutionWithStep(string a_name, double a_N, double a_mu, double a_sig, double a_width);


  string getGaussianFunctionString(int a_startParam = 0);
  double doubleGauss(double a_N, double a_mu, double a_sig1, double a_fracN, double a_sig2,double a_x); //this function will never be negative (N looses meaning if r != [0,1])
  TF1* getDoubleGaussianFunction(string a_name="funct",double a_N=1, double a_mu=0, double a_sig1=1, double a_fracN=0.5, double a_sig2=2);
  string getDoubleGaussianString(int a_startParam = 0);
  double gaussAndStudentT(double a_N, double a_mu, double a_sig1, double a_nu, double a_fracN, double a_sig2, double a_x);
  TF1* getGaussAndStudentTFunction(string a_name="funct",double a_N=1, double a_mu=0, double a_sig1=1, double a_nu=5, double a_fracN=0.5, double a_sig2=1.1);
  string getGaussAndStudentTFunctionString(int a_startParam = 0);

  TF1* getFourierTrigFunction(string a_name, int a_order, double a_xScale = 1.0); // order is sin(wx) cos(wx) sin(2wx) cos(2wx) ...
  //void initializeFourierTrigFunction(TF1* a_funct, TGraphAsymmErrors* a_graph, double a_lowRange = -999.999, double a_highRange = -999.999); //sets parameter limits as well, leave range blank if you want to use whole range of graph.
  TF1* getPolynomialFunction(string a_name, int a_order);
  TF1* getLegendreExpansion(string a_name, int a_order,bool a_fixRange = false, double a_lowRange = -1.0, double a_highRange = 1.0);// the ranges are the last two parameters that you can set later if you want
  //TF1* getSemiHermiteExpansion(string a_name, int a_order,bool a_fixRange = false, double a_lowRange = -1.0, double a_highRange = 1.0);// the ranges are the last two parameters that you can set later if you want
  TF1* getTsallisFunction(string a_name,double a_dNdy,double a_T,double a_q,double a_m); // this if for 1/2pimT * dN/dmtdy   q>1 always
  TF1* getTsallisFunction_NotInvariant(string a_name,double a_dNdy,double a_T,double a_q,double a_m);

  TF1* skewNormalFunct(double a_N, double a_mu, double a_sig, double a_skew);
  TF1* generalNormalV2Funct(double a_N, double a_mu, double a_sig, double a_skew);
  TF1* pearsonIVFunct(double a_N = 1.0, double a_lambda = 0,double a_alpha = 1.0, double a_m = 1.0,double a_nu = 0.0, bool a_normalize = true); // lambda ~ mode     alpha ~ stdDev      m ~ kurtosis       nu ~ skewness
  string getLegendreExpansionString(int a_order);



  double blastWaveOperator(double *x, double *par);
  TF1* blastWaveFunct(string a_name, double a_A, double a_T, double a_BetaS, double a_n, double a_R, double a_mass);
  double blastWave(double a_mTm0, double a_A, double a_T, double a_BetaS, double a_n, double a_R, double a_mass);
  // This will return a normal blast wave and multiplies by 2pimT
  double blastWave_NotInvariant(double a_mTm0, double a_A, double a_T, double a_BetaS, double a_n, double a_R, double a_mass);
  TF1* blastWaveFunct_NotInvariant(string a_name, double a_A, double a_T, double a_BetaS, double a_n, double a_R, double a_mass);
  double blastWave_betaConstr_NotInvariant(double a_mTm0, double a_A, double a_T, double a_BetaS_p4, double a_BetaS_shift, double a_BetaS_p2, double a_BetaS_const, double a_n, double a_R, double a_mass, double a_rapidity);
  TF1* blastWaveFunct_betaConstr_NotInvariant(string a_name, double a_A, double a_T, double a_BetaS_p4, double a_BetaS_shift, double a_BetaS_p2, double a_BetaS_const, double a_n, double a_R, double a_mass, double a_rapidity);
  TF1* blastWaveFunctPlusFitSafeExponentialFunct(string a_name,double a_A, double a_T, double a_BetaS, double a_n, double a_mass, double a_scale, double a_exponent);
  TF1* parabolaFunct(string a_name, double a_yLoc, double a_ySlope, double a_ySecDer, double a_xLoc);
  TF1* parabolaFunct(string a_name, double a_x1, double a_y1, double a_x2, double a_y2, double a_x3, double a_y3);

  TF1* linearFunct(string a_name, double a_yLoc, double a_ySlope, double a_xLoc);
  TF1* linearFunct(string a_name, double a_x1, double a_y1, double a_x2, double a_y2); // will convert to xy loc and slope like above funciton
  TF1* piecewiseTurnOnFunct(string a_name, double a_A, double a_B, double a_C, double a_xCut, double a_ySlopeCut, double a_m);
  TF1* piecewiseTurnOnTwoLinesFunct(string a_name,double a_A, double a_B, double a_C,double a_x0_cut, double a_x1_cut, double a_y1_cut,  double a_ySlope1, double a_m);


  TF1* getInverseCoshFunction(string a_name, double a_Amp,double a_mu, double a_width, bool a_fixMu = true);

  double doubleBoseEinstein_CoulombCorrected_NotInvariant_Plus(double a_mTm0, double a_ALow, double a_TALow, double a_TmLow, double a_TsLow, double a_VC, double a_AHigh, double a_TAHigh, double a_TsHigh, double a_rapVal, double a_centIndex);
  double doubleBoseEinstein_CoulombCorrected_NotInvariant_Minus(double a_mTm0, double a_ALow, double a_TALow, double a_TmLow, double a_TsLow, double a_VC, double a_AHigh, double a_TAHigh, double a_TsHigh, double a_rapVal, double a_centIndex);
  TF1* doubleBoseEinsteinFunct_CoulombCorrected_NotInvariant(string a_name, double a_ALow, double a_TALow, double a_TmLow, double a_TsLow, double a_VC, double a_AHigh, double a_TAHigh, double a_TsHigh, double a_rapVal, int a_charge, double a_centIndex);

  double doubleBoseEinstein_CoulombCorrected_FreeTemps_NotInvariant_Plus(double a_mTm0, double a_ALow, double a_TLow,  double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_rapVal, double a_centIndex, double a_yCM);
  double doubleBoseEinstein_CoulombCorrected_FreeTemps_NotInvariant_Minus(double a_mTm0, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_rapVal, double a_centIndex, double a_yCM);
  TF1* doubleBoseEinsteinFunct_CoulombCorrected_FreeTemps_NotInvariant(string a_name, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_rapVal, int a_charge, double a_centIndex, double a_yCM);


  double doubleBoseEinstein_CoulombCorrected_FreeTemps_noMem_NotInvariant_Plus(double a_mTm0, double a_ALow, double a_TLow,  double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_proTemp, double a_rapCM);
  double doubleBoseEinstein_CoulombCorrected_FreeTemps_noMem_NotInvariant_Minus(double a_mTm0, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_proTemp, double a_rapCM);
  TF1* doubleBoseEinsteinFunct_CoulombCorrected_FreeTemps_noMem_NotInvariant(string a_name, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_proTemp, double a_rapCM, int a_charge);


  double boseEinstein_CoulombCorrected_NotInvariant_Plus(double a_mTm0, double a_A, double a_T, double a_mass, double a_VC, double a_rapVal);
  double boseEinstein_CoulombCorrected_NotInvariant_Minus(double a_mTm0, double a_A, double a_T, double a_mass, double a_VC, double a_rapVal);
  TF1*   boseEinsteinFunct_CoulombCorrected_NotInvariant(string a_name, double a_A, double a_T, double a_mass, double a_VC, double a_rapVal, int a_charge);

  TF1* getCeresDnDyFunction(string a_name, double a_dNdy_0, double a_sqrt_sNN, double a_mass);

  //Sigmoid Functions : Turn-on Functions
  TF1* getLogisticFunction();
  TF1* getGeneralizedLogisticFunction();
  TF1* getHyperbolicTangentFunction();
  TF1* getArcTangentFunction();
  TF1* getGundimannianFunction();
  TF1* getErrorFunction();
  TF1* getAlgebraicTurnOnFunction();
  TF1* getAlgebraicTurnOnFunction_2();

  //most probable dE/dx as a function of momentum
  TF1* getBichselFunction(string a_functionName,double a_mass);
  TF1* getLnBichselFunction(string a_functionName,double a_mass); // ln( bichsel ) for ln(dE/dx) Prediction

  //most probable dE/dx
  TF1* getBichselFunctionBetaGamma(string a_functionName);
  TF1* getLnBichselFunctionBetaGamma(string a_functionName); // ln( bichsel ) for ln(dE/dx) Prediction
  TF1* getBichselFunctionInvBetaGammaSqr(string a_functionName);
  TF1* getLnBichselFunctionInvBetaGammaSqr(string a_functionName); // ln( bichsel ) for ln(dE/dx) Prediction
  TF1* getBichselFunctionBetaGamma_ChrisFlores(string a_functionName); // This adds an inverse beta term to the log term
  TF1* getLnBichselFunctionBetaGamma_ChrisFlores(string a_functionName); // This adds an inverse beta term to the log term
  TF1* getBichselFunctionBetaGamma_Harasty(string a_functionName);
  TF1* getLnBichselFunctionBetaGamma_Harasty(string a_functionName);

  TF1* getInvBetaCalibrationFunction_byMtM0(string a_functionName, double a_mass, double a_mass_parameter, double a_rapidity, double a_timing_offset = 0.0);//,  double a_stopping = 0.001, double a_stopping_width = 0.001, double a_scaling = 1.0);

  #ifdef _CPP11_
    void SimultaneousSortDoubleInt( vector< double > a_values, vector< double > a_other);
  #endif

  //##############     Testers     ###################
  void testTurnOnFuncts();



  void testLegendreFunctions();

  void testTsallisNormalization();
  void testGausLikeFunctions();
  void testCrossProduct();



}


#endif
