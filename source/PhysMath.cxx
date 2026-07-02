#ifndef _PHYS_MATH_CXX_
#define _PHYS_MATH_CXX_

#include "../headers/PhysMath.h"
//#include "../headers/HistogramUtilities.h"
#include "Math/GenVector/Boost.h"

double PhysMath:: dmtde (double a_y){
    return pow(TMath::CosH(a_y),-1.0);
}

double PhysMath:: dthetadomega2pi(double a_theta){
    return pow(TMath::Sin(a_theta),-1.0);
}

double PhysMath:: dydtheta(double a_m, double a_pt, double a_theta){
    double numerator = -a_pt*pow(pow(TMath::Sin(a_theta),-1.0), 2.0);
    double denominator = pow(pow(a_m,2.0)+pow(a_pt*pow(TMath::Sin(a_theta),-1.0), 2.0),0.5);
    return numerator/denominator;
}

double PhysMath::invBeta(double a_m, double a_p){
  if(a_m < 0.0){
    return sqrt(-a_m*a_m + a_p*a_p)/a_p;
  }else{
    return sqrt(a_m*a_m + a_p*a_p)/a_p;
  }
}

double PhysMath::invBeta_from_MassSquaredSigned(double a_signed_m2, double a_p){
  double pSqr = a_p*a_p;
  if(a_signed_m2 + pSqr < 0.0) return 0.0;
  else return sqrt((a_signed_m2 + pSqr)/pSqr);
}

double PhysMath::invBetaWidth_from_MassSquaredWidth(double a_signed_m2, double a_width, double a_p){
  double pSqr = a_p*a_p;
  double invBeta_low = sqrt((a_signed_m2 - a_width + pSqr)/pSqr);
  double invBeta_high = sqrt((a_signed_m2 + a_width + pSqr)/pSqr);
  return 0.5*(invBeta_high - invBeta_low);
}


double PhysMath::beta_from_MassSquaredSigned(double a_signed_m2, double a_p){
  double pSqr = a_p*a_p;
  return sqrt(pSqr/(a_signed_m2 + pSqr));
}

double PhysMath::rapidity(double a_m,double a_pT,double a_pZ){
 double energy = TMath::Sqrt(a_pT*a_pT + a_pZ*a_pZ +a_m*a_m);
 return 0.5*TMath::Log((energy + a_pZ)/(energy - a_pZ));
}


void PhysMath::rapidity_mTm0(ROOT::Math::PxPyPzEVector* a_vect, double &a_y, double &a_mTm0){
  if(!a_vect){
    cout << "ERROR: PhysMath::rapidity_mTm0 handed invalid pointer: " << a_vect << endl;
    a_y = -999;
    a_mTm0 = -999;
  }
  a_y         = a_vect->Rapidity();
  double pTot = a_vect->P();
  double mass = sqrt(a_vect->E()*a_vect->E() - pTot*pTot );
  a_mTm0      = PhysMath::mTm0(mass,a_vect->Pt());
}

double PhysMath::rapFromEtaPt(double a_m, double a_eta, double a_pt){
  double numerator = sqrt(a_m*a_m + a_pt*a_pt*TMath::CosH(a_eta)*TMath::CosH(a_eta))
                     + a_pt * TMath::SinH(a_eta);
  double denominator = TMath::Sqrt(a_m*a_m + a_pt*a_pt);
  return TMath::Log(numerator/denominator);
}

double PhysMath::rapCM(double y_lab, double yCM){
  return (y_lab - yCM);
}

double PhysMath::mass_fromBeta(double a_p, double a_beta){
  return sqrt(a_p*a_p*(1.0/(a_beta*a_beta)-1.0));
}

double PhysMath::massSquared_signed_fromBeta(double a_p, double a_beta){
  return a_p*a_p*(1.0/(a_beta*a_beta) - 1.0);
}

double PhysMath::mass_signed_fromBeta(double a_p, double a_beta){
  double m2 = a_p*a_p*(1.0/(a_beta*a_beta) - 1.0);
  if(m2 > 0) return sqrt(m2);
  else       return -sqrt(-m2);
}

double PhysMath::invBetaGammaSqr_signed_fromInvBeta(double a_invBeta){
  return a_invBeta*a_invBeta - 1.0;
}

double PhysMath::massSquared_signed_fromInvBeta(double a_p, double a_invBeta){
  return a_p*a_p*(a_invBeta*a_invBeta - 1.0);
}

double PhysMath::massSquaredWidth_from_invBetaWidth(double a_p, double a_invBeta, double a_width){
  double pSqr = a_p*a_p;
  double m2_low = pSqr*((a_invBeta - a_width)*(a_invBeta - a_width) - 1.0);
  double m2_high = pSqr*((a_invBeta + a_width)*(a_invBeta + a_width) - 1.0);
  return 0.5*(m2_high - m2_low);
}


double PhysMath::mTm0(double a_m, double a_pt){
  return sqrt(a_m*a_m + a_pt*a_pt) - a_m;
}
double PhysMath::pT(double a_px, double a_py){
  return sqrt(a_px*a_px + a_py*a_py);
}

double PhysMath::pT_from_mTm0(double a_mTm0, double a_m){
  return sqrt((a_mTm0+a_m)*(a_mTm0+a_m)-a_m*a_m);
}

double PhysMath::pZ_from_y_mTm0(double a_y, double a_mTm0, double a_m){
  return (a_mTm0+a_m)*TMath::SinH(a_y);
}

double PhysMath::pZ_from_eta_pT(double a_eta, double a_pT){
  return a_pT*TMath::SinH(a_eta);
}

double PhysMath::pTot(double a_px, double a_py, double a_pz){
  return sqrt(a_px*a_px + a_py*a_py + a_pz*a_pz);
}


double PhysMath::pTot_from_y_mTm0(double a_y, double a_mTm0, double a_m){
  double pT = pT_from_mTm0(a_mTm0,a_m);
  double pZ = pZ_from_y_mTm0(a_y,a_mTm0,a_m);
  return TMath::Sqrt(pT*pT + pZ*pZ);
  //return sqrt(pow((a_mTm0 + a_m)*TMath::CosH(a_y),2) - a_m*a_m);
}

double PhysMath::pTot_from_eta_pT(double a_eta, double a_pT){
  return a_pT*TMath::CosH(a_eta);
}



double PhysMath::eta_from_p(double a_px, double a_py, double a_pz){
  double pTot = sqrt(a_px*a_px + a_py*a_py + a_pz*a_pz);
  double pT = sqrt(a_px*a_px + a_py*a_py);
  double eta = TMath::ACosH(pTot/pT);
  if(a_pz < 0.0) eta = -eta;
  return eta;
}

double PhysMath::theta_from_eta(double a_eta){
  return 2.0*atan(exp(-a_eta)); // exp always positive, so no need for atan2. Theta is from positive axis, so eta = 0 at pi/2
}

double PhysMath::eta_from_theta(double a_theta){
  return -log(tan(0.5*a_theta));
}

double PhysMath::theta_from_y_mTm0(double a_y, double a_mTm0, double a_m){
  return PhysMath::theta_from_eta(PhysMath::eta_from_y_mTm0(a_y,a_mTm0,a_m));
}


double PhysMath::invBeta_width_prediction(double a_timeRes,double a_y, double a_mTm0, double a_m, double a_bTOFRadius){
  return sin(PhysMath::theta_from_y_mTm0(a_y,a_mTm0,a_m))*a_timeRes*299792458.0/a_bTOFRadius;
}


double PhysMath::m2_width_prediction(double a_timeRes,double a_y, double a_mTm0, double a_m, double a_bTOFRadius){
  //double sigmaInvBeta = sin(PhysMath::theta_from_y_mTm0(a_y,a_mTm0,a_m))*a_timeRes*299792458.0/a_bTOFRadius;
  double mom = PhysMath::pTot_from_y_mTm0( a_y,  a_mTm0,  a_m);
  //double invBeta = PhysMath::invBeta(a_m,mom);
  //return 2.0*pow(mom,2.0)*invBeta*sigmaInvBeta;
  return 2.0*mom*sin(PhysMath::theta_from_y_mTm0(a_y,a_mTm0,a_m))*sqrt(a_m*a_m + mom*mom)*a_timeRes*299792458.0/a_bTOFRadius;
}

double PhysMath::m2_width_prediction(double a_InvBetaRes,double a_p, double a_m){
  double invBeta   = PhysMath::invBeta(a_m,a_p);
  double lowMass2  = PhysMath::massSquared_signed_fromInvBeta(a_p,invBeta - a_InvBetaRes);
  double highMass2 = PhysMath::massSquared_signed_fromInvBeta(a_p,invBeta + a_InvBetaRes);
  #ifdef _PHYS_MATH_DEBUG_
    cout << "Call to PhysMath::m2_width_prediction(a_InvBetaRes=" << a_InvBetaRes << ",a_p=" << a_p << ",a_m=" << a_m << ")" << endl;
    cout << "    invBeta -> " << invBeta - a_InvBetaRes << " < " << invBeta << " < " << invBeta + a_InvBetaRes << endl;
    cout << "    m2 ->  " << lowMass2 << " < " << a_m*a_m << " < " << highMass2 << endl;
  #endif
  return 0.5*(highMass2-lowMass2);
}

bool PhysMath::tofSelectionCut(double a_p,double a_m_0,double a_inverseBeta,double a_cut){
  return fabs(a_inverseBeta - sqrt((a_m_0*a_m_0)/(a_p*a_p) + 1.0)) < a_cut;
}


double PhysMath::eta_from_y_mTm0(double a_y, double a_mTm0, double a_m){
  double pT = pT_from_mTm0(a_mTm0,a_m);
  double pTot = pTot_from_y_mTm0(a_y,a_mTm0,a_m);
  double eta = TMath::ACosH(pTot/pT);
  if(a_y < 0.0) eta = -eta;
  return eta;
}

// mT-m0 from rapidity(y) and momentum(p) and mass (m)
double PhysMath::mTm0_from_y_p(double a_y, double a_p, double a_m){
  double sinhY = TMath::SinH(a_y);
  return TMath::Sqrt( (a_p*a_p + a_m*a_m)
   /(1.0 + sinhY*sinhY)
   ) - a_m;
}

double PhysMath::mT_from_mTm0_m(double a_mTm0, double a_m){
  return a_mTm0 + a_m;
}


//this is the converstion factor between dN/dy*dpT = fact * dN/dn*dpT
double PhysMath::eta_to_rap_factor(double a_m, double a_eta, double a_pt){
  double rap = rapFromEtaPt(a_m, a_eta, a_pt);
  double inside = a_m/(TMath::Sqrt(a_pt*a_pt + a_m*a_m)* TMath::CosH(rap));
  return  1.0/TMath::Sqrt(1.0 - inside*inside);
}

// converstion factor from dN/dy*dpT = fact
double PhysMath::pt_to_mT_factor(double a_m, double a_pt){
  double mT = mTm0(a_m,a_pt) + a_m;
  return mT/TMath::Sqrt(mT*mT - a_m*a_m);
}

double PhysMath::gaus(double a_mean, double a_sig, double a_x){
  double val = 0.3989422804 * TMath::Exp(-0.5*( a_x - a_mean )*( a_x - a_mean )/(a_sig*a_sig)) /(a_sig);
  if(std::isfinite(val)) return val;
  return 0.0;
}

double PhysMath::gaus_finite_integral(double a_N, double a_mu, double a_sig, double a_x_low, double a_x_high){
  double sqrtTwoSigma = physMath_sqrtTwo*a_sig;
  return a_N*(erf((a_x_high - a_mu)/sqrtTwoSigma) - erf((a_x_low - a_mu)/sqrtTwoSigma));
}

double PhysMath::logNorm(double a_mu, double a_sig, double a_x){
  return 0.3989422804 * TMath::Exp(
   -0.5*( TMath::Log(a_x) - a_mu )*( TMath::Log(a_x) - a_mu )/(a_sig*a_sig)
   ) /(a_sig * a_x);
}

void PhysMath::gausToLogNorm(double a_mean, double a_stdDev, double &a_mu, double &a_sig){
  a_mu = 2.0*log(a_mean) - 0.5*log(a_stdDev*a_stdDev + a_mean*a_mean);
  a_sig = sqrt( log(1.0 + a_stdDev*a_stdDev/(a_mean*a_mean)) );
}


void PhysMath::logNormToGaus(double a_mu, double a_sig, double &a_mean, double &a_stdDev){
  a_mean = exp(a_mu+0.5*a_sig*a_sig);
  a_stdDev = a_mean*sqrt(exp(a_sig*a_sig) - 1.0);
}

double PhysMath::student_T(double a_mu, double a_sig, double a_nu, double a_x){
  return TMath::Student((a_x - a_mu)/a_sig, a_nu)/a_sig;
}

double PhysMath::student_T_BinContent(double a_mu, double a_sig, double a_nu, double a_x_low, double a_x_high){
  double width = a_x_high - a_x_low;
  double step = 0.1*width;
  double content = 0.0;
  for(double point = a_x_low + 0.5*step; point < a_x_high; point += step){
    content += student_T(a_mu, a_sig, a_nu, point);
  }
  return 0.1*content;
}

#ifdef _MATHMORE_
double PhysMath::student_T_finite_integral(double a_N,double a_mu, double a_sig, double a_nu, double a_x_low, double a_x_high){
  double k = 0.5*(a_nu+1.0);
  return exp(ROOT::Math::lgamma(k) - ROOT::Math::lgamma(0.5*a_nu))*(ROOT::Math::hyperg(0.5,k,1.5,-a_x_high/a_nu) - ROOT::Math::hyperg(0.5,k,1.5,-a_x_high/a_nu))/sqrt(TMath::Pi()*a_nu);
}
#endif

double PhysMath::m2_student_T(double a_mu, double a_sig, double a_nu, double a_x, double a_momentum){
  double invBeta = sqrt(a_x/(a_momentum*a_momentum) + 1.0);
  return TMath::Student((invBeta - a_mu)/a_sig, a_nu)/(2.0*a_sig*a_momentum*a_momentum*invBeta);
}

double PhysMath::m2_student_T_BinContent(double a_mu, double a_sig, double a_nu, double a_x_low, double a_x_high, double a_momentum){
  double width = a_x_high - a_x_low;
  double step = 0.1*width;
  double content = 0.0;
  for(double point = a_x_low + 0.5*step; point < a_x_high; point += step){
    content += m2_student_T(a_mu, a_sig, a_nu, point, a_momentum);
  }
  return 0.1*content;
}





double PhysMath::skewNormal(double a_mu, double a_sig, double a_skew, double a_x){
  double skew = a_skew;
  if(skew == 1.0) skew = 0.99999999999;
  else if(skew == -1.0) skew = -0.99999999999;
  double alpha = skew/sqrt(fabs(1.0 - skew*skew));
  double loc;
  if(alpha > 0.0)        loc = a_mu - a_sig/sqrt(1.0 - 2.0*skew*skew/TMath::Pi())*(skew*sqrt(2.0/TMath::Pi()) - (4.0-TMath::Pi())/4.0*pow(skew*sqrt(2.0/TMath::Pi()),3)/(1.0-2.0*skew*skew/TMath::Pi()) - 0.5*exp(-2.0*TMath::Pi()/fabs(alpha)));
  else if(alpha < 0.0)   loc = a_mu - a_sig/sqrt(1.0 - 2.0*skew*skew/TMath::Pi())*(skew*sqrt(2.0/TMath::Pi()) - (4.0-TMath::Pi())/4.0*pow(skew*sqrt(2.0/TMath::Pi()),3)/(1.0-2.0*skew*skew/TMath::Pi()) + 0.5*exp(-2.0*TMath::Pi()/fabs(alpha)));
  else if(alpha == 0.0)  loc = a_mu - a_sig/sqrt(1.0 - 2.0*skew*skew/TMath::Pi())*(skew*sqrt(2.0/TMath::Pi()) - (4.0-TMath::Pi())/4.0*pow(skew*sqrt(2.0/TMath::Pi()),3)/(1.0-2.0*skew*skew/TMath::Pi()));
  double w = a_sig/sqrt(1.0 - 2.0*skew*skew/TMath::Pi());
  double val = 1.0/(w*sqrt(2.0*TMath::Pi()))*exp(-0.5*( a_x - loc)*( a_x - loc )/(w*w)) * (erf(alpha*(a_x-loc)/(w*sqrt(2.0)))+1.0);
  if(std::isfinite(val)) return val;
  return 0.0;
}


double PhysMath::skewNormalBinContent(double a_mu, double a_sig, double a_skew, double a_xLow,double a_xHigh){
  double width = a_xHigh - a_xLow;
  double step = 0.1*width;
  double content = 0.0;
  for(double point = a_xLow + 0.5*step; point < a_xHigh; point += step){
    content += skewNormal(a_mu, a_sig, a_skew, point);
  }
  return 0.1*content;
}



#ifdef _BOOST_
double PhysMath::skewNormal_finite_integral(double a_N,double a_mu, double a_sig, double a_alpha, double a_x_low, double a_x_high){
  double sqrtTwoSigma = physMath_sqrtTwo*a_sig;
  double arg_low  = (a_x_low - a_mu)/sqrtTwoSigma;
  double arg_high = (a_x_high - a_mu)/sqrtTwoSigma;
  return a_N*(erf(arg_high) - 2.0*boost::math::owens_t(arg_high*physMath_sqrtTwo,a_alpha) - erf(arg_low) + 2.0*boost::math::owens_t(arg_low*physMath_sqrtTwo,a_alpha));
}
#endif


/*double PhysMath::generalNormalV2(double a_mode, double a_sigma, double a_skew, double a_x){
    double hold = pow(a_skew,2)+sqrt(pow(a_skew,4)+4.0*pow(a_skew,2))+2.0;
    double k = 0.0;
    if     (a_skew >  0.0)  k = -sqrt(log(pow(hold/2.0,1.0/3.0) + pow(2.0/hold,1.0/3.0) - 1.0));
    else if(a_skew <  0.0)  k =  sqrt(log(pow(hold/2.0,1.0/3.0) + pow(2.0/hold,1.0/3.0) - 1.0));
    double alpha = a_sigma*k*exp(-0.5*k*k)/sqrt(exp(k*k) - 1.0);
    double xi    = (a_mode*k*exp(k*k) + alpha - alpha*exp(k*k))*exp(-k*k)/k;
    //if(mirr == 1.0){
    //   k = -k;
    //    xi = -xi;
    //}* /
    //if(fabs(a_x) > fabs(xi + alpha/k)) return 0.0;
    double gausFunct = exp(-0.5*pow(-log(1.0-k*(a_x - xi)/alpha)/k,2));
    double denom     = (alpha - k*(a_x - xi))*sqrt(2.0*TMath::Pi());
    return gausFunct/denom;
}*/

double PhysMath::generalNormalV2(double a_mode, double a_sigma, double a_skew, double a_x){
    double hold = pow(a_skew,2)+sqrt(pow(a_skew,4)+4.0*pow(a_skew,2))+2.0;
    double k = 0.0;
    double oneThird = 1.0/3.0;
    if     (a_skew >  0.0)  k = -sqrt(log(pow(0.5*hold,oneThird) + pow(2.0/hold,oneThird) - 1.0));
    else if(a_skew <  0.0)  k =  sqrt(log(pow(0.5*hold,oneThird) + pow(2.0/hold,oneThird) - 1.0));
    double k_squared = k*k;
    double alpha = a_sigma*k*exp(-0.5*k_squared)/sqrt(exp(k_squared) - 1.0);
    double xi    = (a_mode*k + alpha*(exp(-k_squared)-1.0))/k;
    //double xi    = (a_mode*k*exp(k*k) + alpha - alpha*exp(k*k))*exp(-k*k)/k;
    double gausFunct = exp(-0.5*pow(-log(1.0-k*(a_x - xi)/alpha)/k,2));
    double denom     = (alpha - k*(a_x - xi))*sqrt(2.0*TMath::Pi());
    return gausFunct/denom;
}


double PhysMath::LogGammaR2(double a_x, double a_y){
  //CDF/MEMO/STATISTICS/PUBLIC/6820    A Guide to the Pearson Type IV Distribution    Joel Heinrich—University of Pennsylvania    December 21, 2004
  //https://cran.r-project.org
  /* returns log(abs(gamma(x+iy)/gamma(x))^2) */
  const double y2 = a_y*a_y;
  const double xmin = (2.0*y2 > 10.0) ? 2.0*y2 : 10.0;
  double r = 0.0;
  double s = 1.0;
  double p = 1.0;
  double f = 0.0;
  while(a_x < xmin) {
    const double t = a_y/a_x++;
    r += log(1.0 + t*t);
  }
  while (p > s*DBL_EPSILON) { // DBL_EPSILON = 2.2204460e-16 is defined somewhere in c++ or root...
    p *= y2 + f*f;
    p /= a_x++ * ++f;
    s += p;
  }
  return -r-log(s);
}

/*
//Translated from R to C++, but it doesn't normalize to 1, but 1.1 or something close... cant figure out why
double PhysMath::pearsonIV(double a_lambda,double a_alpha,double a_m,double a_nu,double  a_x, bool a_normalize){
  //CDF/MEMO/STATISTICS/PUBLIC/6820    A Guide to the Pearson Type IV Distribution    Joel Heinrich—University of Pennsylvania    December 21, 2004
  //https://cran.r-project.org
  if(a_alpha <= 0 || a_m <= 0.5) return 0.0;
  int logNormalization = 0.0;
  if(a_normalize) logNormalization = -0.5723649429247000870717136756765293558236474064576557857568115357 + LogGammaR2(a_m,0.5*a_nu) + ROOT::Math::lgamma(a_m) - ROOT::Math::lgamma(a_m-0.5) - log(a_alpha); // -0.57 is -ln(sqrt(Pi))  Their way
  double location = (a_x-a_lambda)/a_alpha;
  return exp(logNormalization - a_m*log(1.0+(location*location)) - a_nu*TMath::ATan(location));
}*/

double PhysMath::pearsonIV(double a_lambda,double a_alpha,double a_m,double a_nu, double a_x, bool a_normalize){
  double loc = (a_x - a_lambda)/a_alpha;
  //double A = 1.0 + loc*loc;
  //double B = -TMath::ATan(loc);
  //double lnShape =  0.5*(a_nu-a_m)*(log(A)+B) + 0.5*(a_nu+a_m)*(B - log(A)); // one blows up the other is stable?
  //double lnShape =  a_nu*B - a_m*log(A); // one blows up the other is stable?
  double lnShape =  -a_nu*TMath::ATan(loc) - a_m*log(1.0 + loc*loc); // one blows up the other is stable?

  if(a_normalize){ // can take a long time for certian parameters
    double lnNorm = -0.5723649429247000870717136756765293558236474064576557857568115357 - ROOT::Math::lgamma(a_m-0.5) + ROOT::Math::lgamma(a_m);
    if(a_nu != 0.0) lnNorm += LogGammaR2(a_m,0.5*a_nu);
    return  exp(lnShape + lnNorm)/a_alpha;
  }else{
    double mode = -0.5*a_nu/a_m; // ~mode if lambda = 0 and alpha = 1
    lnShape += a_nu*TMath::ATan(mode) + a_m*log(1.0 + mode*mode);
    return exp(lnShape);
  }

  /*  Original, not safe for high m...
  double safe_m = a_m;
  double safe_nu = a_nu;
  double loc = (a_x - a_lambda)/a_alpha;
  double A = 1.0 + loc*loc;
  double B = exp(-TMath::ATan(loc));
  double shape =  pow(A*B,0.5*(safe_nu-safe_m)) * pow(B/A,0.5*(safe_nu+safe_m)); // one blows up the other is stable
  //double shape = pow(1.0 + pow(loc, 2), -safe_m) * exp( -safe_nu * TMath::ATan(loc) );
  if(a_normalize){
    double normalization = 1.0;
    //if(a_m > 120.0) a_m = 120; // sets a max limit on m
    if(a_nu != 0.0){
      complex<double> comNum(safe_m,0.5*safe_nu);
      complex<double> gammaComNum = PhysMath::ComplexGamma( comNum );
      double gammaM = TMath::Gamma(safe_m);
      complex<double> comFrac(gammaComNum.real()/gammaM, gammaComNum.imag()/gammaM);
      normalization = (comFrac.real()*comFrac.real() + comFrac.imag()*comFrac.imag() );
    }
    return  shape * normalization / (a_alpha * TMath::Beta(safe_m - 0.5,0.5));
  }
  return shape;
  */

}


double PhysMath::pearsonIVAmpToCounts(double a_amplitude,double a_lambda,double a_alpha,double a_m,double a_nu){
  double lnNorm = -0.5723649429247000870717136756765293558236474064576557857568115357 - ROOT::Math::lgamma(a_m-0.5) + ROOT::Math::lgamma(a_m);
  if(a_nu != 0.0) lnNorm += LogGammaR2(a_m,0.5*a_nu);
  return a_amplitude*a_alpha/exp(lnNorm);
}

double PhysMath::piecewiseTurnOn(double a_A, double a_B, double a_C, double a_xCut,double a_ySlopeCut, double a_m, double a_x){
  // A*exp(-b*pT^-c)    x = mT-m0
  if(a_x < a_xCut){ // use the turn on function
    return a_A*exp(-a_B*pow(fabs(pow(a_x+a_m,2.0)-a_m*a_m),-a_C*0.5));
  } // use a line
  double y_0 = a_A*exp(-a_B*pow(fabs(pow(a_xCut+a_m,2.0)-a_m*a_m),-a_C*0.5));
  return a_ySlopeCut*(a_x - a_xCut) + a_A*exp(-a_B*pow((pow(a_xCut+a_m,2.0)-a_m*a_m),-a_C*0.5));

}

double PhysMath::piecewiseTurnOnTwoLines(double a_A, double a_B, double a_C,double a_x0_cut, double a_x1_cut, double a_y1_cut,  double a_ySlope1, double a_m, double a_x){
  // A*exp(-b*pT^-c)    x = mT-m0  put m = 0 for a mTm0 instead of pT
  double x_low = a_x0_cut;
  double x_high = a_x1_cut;
  if(x_low > x_high){
    x_low = a_x1_cut;
    x_high = a_x0_cut;
  }
  if(a_x < x_low){ // use the turn on function
    return a_A*exp(-a_B*pow(fabs(pow(a_x+a_m,2.0)-a_m*a_m),-a_C*0.5));
  }else if(a_x < x_high){ // use a line
    double y_0 = a_A*exp(-a_B*pow(fabs(pow(x_low+a_m,2.0)-a_m*a_m),-a_C*0.5));
    return (a_y1_cut - y_0)*(a_x - x_high)/(x_high - x_low) + a_y1_cut;
  }// use a different line
  return a_ySlope1*(a_x - x_high) + a_y1_cut;

}


double PhysMath::pearsonIVOperator(double* a_x,const double* a_params){
  return a_params[0]*PhysMath::pearsonIV(a_params[1], a_params[2], a_params[3], a_params[4], a_x[0] );
}

TF1* PhysMath::skewNormalFunct(double a_N, double a_mu, double a_sig, double a_skew){
  TF1* funct = new TF1("skewNormal","[0]*PhysMath::skewNormal([1],[2],[3],x)",-10,10);
  funct->SetParNames("N","#mu","#sigma","#gamma");
  funct->SetParameter(0, a_N);
  funct->SetParameter(1, a_mu);
  funct->SetParameter(2, a_sig);
  funct->SetParameter(3, a_skew);
  funct->SetNpx(300);
  //funct->SetParLimits(0,0.0,10.0*a_N);
  //funct->SetParLimits(1,-10.0,10.0);
  //funct->SetParLimits(2,0.0, 300);
  //funct->SetParLimits(3,-20, 20);
  return funct;
}


TF1* PhysMath::generalNormalV2Funct(double a_N, double a_mode, double a_sig, double a_skew){
  TF1* funct = new TF1("generalNormalV2Funct","[0]*PhysMath::generalNormalV2([1],[2],[3],x)",-10,10);
  funct->SetParNames("N","mode","#sigma","#gamma");
  funct->SetParameter(0, a_N);
  funct->SetParameter(1, a_mode);
  funct->SetParameter(2, a_sig);
  funct->SetParameter(3, a_skew);
  funct->SetNpx(300);
  funct->SetParLimits(0,0.0,10.0);
  funct->SetParLimits(1,-10.0,10.0);
  funct->SetParLimits(2,0.0, 300);
  funct->SetParLimits(3,-20, 20);
  return funct;
}



TF1* PhysMath::pearsonIVFunct(double a_N,  double a_lambda,double a_alpha, double a_m,double a_nu, bool a_normalize){
  TF1* funct = NULL;
  if(a_normalize){
    funct = new TF1("pearsonIV","[0]*PhysMath::pearsonIV([1],[2],[3],[4],x,true)",-10,10);
  }else{
    funct = new TF1("pearsonIV","[0]*PhysMath::pearsonIV([1],[2],[3],[4],x,false)",-10,10);
  }

  funct->SetParNames("N","#lambda","#alpha","m","#nu");
  funct->SetParameter(0, a_N);
  funct->SetParameter(1, a_lambda);
  funct->SetParameter(2, a_alpha);
  funct->SetParameter(3, a_m);
  funct->SetParameter(4, a_nu);
  funct->SetNpx(300);

  funct->SetParLimits(2,0.0, 300);
  funct->SetParLimits(3,0.01, 3.0);
  funct->SetParLimits(4,-3.0, 3.0);


  return funct;
}


TF1* PhysMath::getPearsonAndTwoGausFunct(string a_functName, double a_N , double a_lambda, double a_alpha , double a_m , double a_nu ,
                                         double a_N_g1 , double a_mu_g1 , double a_sig_g1,
                                         double a_N_g2 , double a_mu_g2 , double a_sig_g2,bool a_normalizePearson){
  TF1* funct = new TF1(a_functName.c_str(),Form(
    "[0]*PhysMath::pearsonIV([1],[2],[3],[4],x,%s) + [5]*0.3989422804 * TMath::Exp(-0.5*( x - [6] )*(x - [6])/([7]*[7])) /[7] + [8]*0.3989422804 * TMath::Exp(-0.5*( x - [9] )*(x - [9])/([10]*[10])) /[10]", a_normalizePearson ? "true" : "false"),-10,10);


  funct->SetParNames("N_{1}","#lambda","#alpha","m","#nu","N_{2}","#mu_{2}","#sigma_{2}","N_{3}","#mu_{3}","#sigma_{3}");
  funct->SetParameter(0, a_N);
  funct->SetParameter(1, a_lambda);
  funct->SetParameter(2, a_alpha);
  funct->SetParameter(3, a_m);
  funct->SetParameter(4, a_nu);
  funct->SetParameter(5, a_N_g1);
  funct->SetParameter(6, a_mu_g1);
  funct->SetParameter(7, a_sig_g1);
  funct->SetParameter(8, a_N_g2);
  funct->SetParameter(9, a_mu_g2);
  funct->SetParameter(10, a_sig_g2);
  funct->SetNpx(300);

  funct->SetParLimits(2,0.0, 300);
  funct->SetParLimits(3,0.01, 3.0);
  funct->SetParLimits(4,-3.0, 3.0);

  return funct;
}


TF1* PhysMath::getPearsonAndTwoGausFunct_sigmaLinked(string a_functName, double a_N , double a_lambda, double a_alpha , double a_m , double a_nu ,
                                         double a_N_g1 , double a_mu_g1 , double a_sig_g1,
                                         double a_N_g2 , double a_mu_g2 ,bool a_normalizePearson){
  TF1* funct = new TF1(a_functName.c_str(),Form(
    "[0]*PhysMath::pearsonIV([1],[2],[3],[4],x,%s) + [5]*0.3989422804 * TMath::Exp(-0.5*( x - [6] )*(x - [6])/([7]*[7])) /[7] + [8]*0.3989422804 * TMath::Exp(-0.5*( x - [9] )*(x - [9])/([7]*[7])) /[7]", a_normalizePearson ? "true" : "false"),-10,10);
  funct->SetParNames("N_{1}","#lambda","#alpha","m","#nu","N_{2}","#mu_{2}","#sigma_{2&3}","N_{3}","#mu_{3}");
  funct->SetParameter(0, a_N);
  funct->SetParameter(1, a_lambda);
  funct->SetParameter(2, a_alpha);
  funct->SetParameter(3, a_m);
  funct->SetParameter(4, a_nu);
  funct->SetParameter(5, a_N_g1);
  funct->SetParameter(6, a_mu_g1);
  funct->SetParameter(7, a_sig_g1);
  funct->SetParameter(8, a_N_g2);
  funct->SetParameter(9, a_mu_g2);
  funct->SetNpx(300);

  funct->SetParLimits(2,0.0, 300);
  funct->SetParLimits(3,0.01, 3.0);
  funct->SetParLimits(4,-3.0, 3.0);

  return funct;
}

TF1* PhysMath::getPearsonAndTwoGausFunct_ratioLinked(string a_functName, double a_N , double a_lambda, double a_alpha , double a_m , double a_nu ,
                                 double a_fixedRatio, double a_mu_g1 , double a_sig_g1,
                                 double a_N_g2 , double a_mu_g2 ,double a_sig_g2, bool a_normalizePearson){
  TF1* funct = new TF1(a_functName.c_str(),Form(
    "[0]*PhysMath::pearsonIV([1],[2],[3],[4],x,%s) + [5]*[8]*PhysMath::gaus([6],[7],x)  + [8]*PhysMath::gaus([9],[10],x) ", a_normalizePearson ? "true" : "false"),-10,10);
  funct->SetParNames("N_{1}","#lambda","#alpha","m","#nu","N_{2}/N_{3}","#mu_{2}","#sigma_{3}","N_{3}","#mu_{3}","#sigma_{3}");
  funct->SetParameter(0, a_N);
  funct->SetParameter(1, a_lambda);
  funct->SetParameter(2, a_alpha);
  funct->SetParameter(3, a_m);
  funct->SetParameter(4, a_nu);
  funct->SetParameter(5, a_fixedRatio);
  funct->SetParameter(6, a_mu_g1);
  funct->SetParameter(7, a_sig_g1);
  funct->SetParameter(8, a_N_g2);
  funct->SetParameter(9, a_mu_g2);
  funct->SetParameter(10, a_sig_g2);

  funct->SetNpx(300);

  funct->SetParLimits(2,0.0, 300);
  funct->SetParLimits(3,0.01, 3.0);
  funct->SetParLimits(4,-3.0, 3.0);
  funct->FixParameter(5,a_fixedRatio);

  return funct;
}



TF1* PhysMath::getSkewAndTwoGausFunct(string a_functName, double a_N , double a_mu, double a_sig , double a_skew,
                                         double a_N_g1 , double a_mu_g1 , double a_sig_g1,
                                         double a_N_g2 , double a_mu_g2 , double a_sig_g2){
  TF1* funct = new TF1(a_functName.c_str(),
    "[0] * PhysMath::skewNormal([1],[2],[3],x) + [4]*PhysMath::gaus([5],[6],x) + [7]*PhysMath::gaus([8],[9],x)",-10,10);
    //"[0] * PhysMath::skewNormal([1],[2],[3],x) + [4]*0.3989422804*exp(-0.5*( x - [5] )*( x - [5] )/([6]*[6]))/[6] + [7]*0.3989422804*exp(-0.5*(x - [8])*( x - [8] )/([9]*[9]))/[9]",-10,10);
  funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","#gamma_{1}","N_{2}","#mu_{2}","#sigma_{2}","N_{3}","#mu_{3}","#sigma_{3}");
  funct->SetParameter(0, a_N);
  funct->SetParameter(1, a_mu);
  funct->SetParameter(2, a_sig);
  funct->SetParameter(3, a_skew);
  funct->SetParameter(4, a_N_g1);
  funct->SetParameter(5, a_mu_g1);
  funct->SetParameter(6, a_sig_g1);
  funct->SetParameter(7, a_N_g2);
  funct->SetParameter(8, a_mu_g2);
  funct->SetParameter(9, a_sig_g2);
  funct->SetNpx(300);

  funct->SetParLimits(0,0.0, 10.0*a_N);
  funct->SetParLimits(1, a_mu - 5.0*a_sig, a_mu + 5.0*a_sig);
  funct->SetParLimits(2, 0.0, 10.0*a_sig);
  funct->SetParLimits(3, -10.0*a_skew,10.0*a_skew);
  funct->SetParLimits(4, 0.0,10.0*a_N_g1);
  funct->SetParLimits(5, a_mu_g1 - 5.0*a_sig_g1,a_mu_g1 + 5.0*a_sig_g1);
  funct->SetParLimits(6, 0.0,10.0*a_sig_g1);
  funct->SetParLimits(7, 0.0,10.0*a_N_g2);
  funct->SetParLimits(8, a_mu_g2 - 5.0*a_sig_g2,a_mu_g2 + 5.0*a_sig_g2);
  funct->SetParLimits(9, 0.0,10.0*a_sig_g2);

  return funct;
}


TF1* PhysMath::getSkewAndTwoGausFunct_sigmaLinked(string a_functName, double a_N , double a_mu, double a_sig , double a_skew,
                                         double a_N_g1 , double a_mu_g1 , double a_sig_g1,
                                         double a_N_g2 , double a_mu_g2){
  TF1* funct = new TF1(a_functName.c_str(),
  "[0] * PhysMath::skewNormal([1],[2],[3],x) + [4]*PhysMath::gaus([5],[6],x) + [7]*PhysMath::gaus([8],[6],x)",-10,10);

  //"[0] * PhysMath::skewNormal([1],[2],[3],x) + [4]*0.3989422804*exp(-0.5*( x - [5] )*( x - [5] )/([6]*[6]))/[6] + [7]*0.3989422804*exp(-0.5*(x - [8])*( x - [8] )/([6]*[6]))/[6]",-10,10);
  funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","#gamma_{1}","N_{2}","#mu_{2}","#sigma_{2&3}","N_{3}","#mu_{3}");
  funct->SetParameter(0, a_N);
  funct->SetParameter(1, a_mu);
  funct->SetParameter(2, a_sig);
  funct->SetParameter(3, a_skew);
  funct->SetParameter(4, a_N_g1);
  funct->SetParameter(5, a_mu_g1);
  funct->SetParameter(6, a_sig_g1);
  funct->SetParameter(7, a_N_g2);
  funct->SetParameter(8, a_mu_g2);
  funct->SetNpx(300);

  funct->SetParLimits(0,0.0, 10.0*a_N);
  funct->SetParLimits(1, a_mu - 5.0*a_sig, a_mu + 5.0*a_sig);
  funct->SetParLimits(2, 0.0, 10.0*a_sig);
  funct->SetParLimits(3, -10.0*a_skew,10.0*a_skew);
  funct->SetParLimits(4, 0.0,10.0*a_N_g1);
  funct->SetParLimits(5, a_mu_g1 - 5.0*a_sig_g1,a_mu_g1 + 5.0*a_sig_g1);
  funct->SetParLimits(6, 0.0,10.0*a_sig_g1);
  funct->SetParLimits(7, 0.0,10.0*a_N_g2);
  funct->SetParLimits(8, a_mu_g2 - 5.0*a_sig_g1,a_mu_g2 + 5.0*a_sig_g1);
  //funct->SetParLimits(9, 0.0,10.0*a_sig_g1);

  return funct;
}

TF1* PhysMath::getSkewAndTwoGausFunct_ratioLinked(string a_functName, double a_N , double a_mu, double a_sig , double a_skew,
                                 double a_N2ToN1or3Ratio, double a_mu_g1 , double a_sig_g1,
                                 double a_N_g2 , double a_mu_g2 ,double a_sig_g2, bool a_Ratio12InsteadOf23){
  TF1* funct = NULL;
  if(a_Ratio12InsteadOf23){
    funct = new TF1(a_functName.c_str(),
      "[0]*PhysMath::skewNormal([1],[2],[3],x) + [4]*[0]*PhysMath::gaus([5],[6],x) + [7]*PhysMath::gaus([8],[9],x)",-10,10);
      //"[0]*PhysMath::skewNormal([1],[2],[3],x) + [4]*[0]*0.3989422804*exp(-0.5*(x-[5])*(x-[5])/([6]*[6]))/[6] + [7]*0.3989422804*exp(-0.5*(x - [8])*(x-[8])/([9]*[9]))/[9]",-10,10);
  }else{
    funct = new TF1(a_functName.c_str(),
      "[0]*PhysMath::skewNormal([1],[2],[3],x) + [4]*[7]*PhysMath::gaus([5],[6],x) + [7]*PhysMath::gaus([8],[9],x)",-10,10);
      //"[0]*PhysMath::skewNormal([1],[2],[3],x)  + [4]*[7]*0.3989422804*exp(-0.5*(x-[5])*(x-[5])/([6]*[6]))/[6] + [7]*0.3989422804*exp(-0.5*(x - [8])*(x-[8])/([9]*[9]))/[9]",-10,10);
  }
  string ratioName = "N_{2}/N_{3}";
  if(a_Ratio12InsteadOf23) ratioName = "N_{2}/N_{1}";
  funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","#gamma_{1}",ratioName.c_str(),"#mu_{2}","#sigma_{2}","N_{3}","#mu_{3}","#sigma_{3}");
  funct->SetParameter(0, a_N);
  funct->SetParameter(1, a_mu);
  funct->SetParameter(2, a_sig);
  funct->SetParameter(3, a_skew);
  funct->SetParameter(4, a_N2ToN1or3Ratio);
  funct->SetParameter(5, a_mu_g1);
  funct->SetParameter(6, a_sig_g1);
  funct->SetParameter(7, a_N_g2);
  funct->SetParameter(8, a_mu_g2);
  funct->SetParameter(9, a_sig_g2);

  funct->SetNpx(300);

  funct->SetParLimits(0,0.0, 10.0*a_N);
  funct->SetParLimits(1, a_mu - 5.0*a_sig, a_mu + 5.0*a_sig);
  funct->SetParLimits(2, 0.0, 10.0*a_sig);
  funct->SetParLimits(3, -10.0*a_skew,10.0*a_skew);
  funct->SetParameter(4,a_N2ToN1or3Ratio);
  funct->SetParLimits(4,0.0,10.0*a_N2ToN1or3Ratio);
  funct->SetParLimits(5, a_mu_g1 - 5.0*a_sig_g1,a_mu_g1 + 5.0*a_sig_g1);
  funct->SetParLimits(6, 0.0,10.0*a_sig_g1);
  funct->SetParLimits(7, 0.0,10.0*a_N_g2);
  funct->SetParLimits(8, a_mu_g2 - 5.0*a_sig_g2,a_mu_g2 + 5.0*a_sig_g2);
  funct->SetParLimits(9, 0.0,10.0*a_sig_g2);
  return funct;

}


TF1* PhysMath::getTwoSkewFunct(string a_functName, double a_N1 , double a_mu1, double a_sig1 , double a_skew1, double a_N2, double a_mu2, double a_sig2, double a_skew2){
  TF1* funct = new TF1(a_functName.c_str(), "[0]*PhysMath::skewNormal([1],[2],[3],x) + [4]*PhysMath::skewNormal([5],[6],[7],x)",-10,10);
  funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","#gamma_{1}","N_{2}","#mu_{2}","#sigma_{2}","#gamma_{2}");
  funct->SetParameter(0, a_N1);
  funct->SetParameter(1, a_mu1);
  funct->SetParameter(2, a_sig1);
  funct->SetParameter(3, a_skew1);
  funct->SetParameter(4, a_N2);
  funct->SetParameter(5, a_mu2);
  funct->SetParameter(6, a_sig2);
  funct->SetParameter(7, a_skew2);

  funct->SetNpx(300);

  funct->SetParLimits(0,0.0, 10.0*a_N1);
  funct->SetParLimits(1, a_mu1 - 5.0*a_sig1, a_mu1 + 5.0*a_sig1);
  funct->SetParLimits(2, 0.0, 10.0*a_sig1);
  funct->SetParLimits(3, -10.0*a_skew1,10.0*a_skew1);
  funct->SetParLimits(4,0.0, 10.0*a_N2);
  funct->SetParLimits(5, a_mu2 - 5.0*a_sig2, a_mu2 + 5.0*a_sig2);
  funct->SetParLimits(6, 0.0, 10.0*a_sig2);
  funct->SetParLimits(7, -10.0*a_skew2,10.0*a_skew2);

  return funct;

}

TF1* PhysMath::getTwoSkewFunct_SafeParam4Amp(string a_functName, double a_N1 , double a_mu1, double a_sig1 , double a_skew1, double a_N2, double a_mu2, double a_sig2, double a_skew2){
  TF1* funct = new TF1(a_functName.c_str(), "[0]*PhysMath::skewNormal([1],[2],[3],x) + tan([4]*0.5*TMath::Pi())*PhysMath::skewNormal([5],[6],[7],x)",-10,10);
  funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","#gamma_{1}","N_{2} safe","#mu_{2}","#sigma_{2}","#gamma_{2}");
  funct->SetParameter(0, a_N1);
  funct->SetParameter(1, a_mu1);
  funct->SetParameter(2, a_sig1);
  funct->SetParameter(3, a_skew1);
  funct->SetParameter(4, 2.0*atan(a_N2)/TMath::Pi());
  funct->SetParameter(5, a_mu2);
  funct->SetParameter(6, a_sig2);
  funct->SetParameter(7, a_skew2);

  funct->SetNpx(300);

  funct->SetParLimits(0,0.0, 10.0*a_N1);
  funct->SetParLimits(1, a_mu1 - 5.0*a_sig1, a_mu1 + 5.0*a_sig1);
  funct->SetParLimits(2, 0.0, 10.0*a_sig1);
  funct->SetParLimits(3, -10.0*a_skew1,10.0*a_skew1);
  funct->SetParLimits(4,0.0, 0.9999999999);
  funct->SetParLimits(5, a_mu2 - 5.0*a_sig2, a_mu2 + 5.0*a_sig2);
  funct->SetParLimits(6, 0.0, 10.0*a_sig2);
  funct->SetParLimits(7, -10.0*a_skew2,10.0*a_skew2);

  return funct;

}


TF1* PhysMath::getThreeGausFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3){
  TF1* funct = new TF1(a_functName.c_str(),
    "[0]*PhysMath::gausn([1],[2],x) + [3]*PhysMath::gausn([4],[5],x) + [6]*PhysMath::gausn([7],[8],x)",-10,10);
  funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","N_{2}","#mu_{2}","#sigma_{2}","N_{3}","#mu_{3}","#sigma_{3}");
  funct->SetParameter(0, a_N1);
  funct->SetParameter(1, a_params[0]);
  funct->SetParameter(2, a_params[1]);
  funct->SetParameter(3, a_N2);
  funct->SetParameter(4, a_params[4]);
  funct->SetParameter(5, a_params[5]);
  funct->SetParameter(6, a_N3);
  funct->SetParameter(7, a_params[8]);
  funct->SetParameter(8, a_params[9]);
  funct->SetNpx(300);
  return funct;
}

TF1* PhysMath::getThreeSkewFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3){
  TF1* funct = new TF1(a_functName.c_str(),
    "[0]*PhysMath::skewNormal([1],[2],[3],x) + [4]*PhysMath::skewNormal([5],[6],[7],x) + [8]*PhysMath::skewNormal([9],[10],[11],x)",-10,10);
  //funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","#gamma_{1}","N_{2}","#mu_{2}","#sigma_{2}","#gamma_{2}","N_{3}","#mu_{3}","#sigma_{3}","#gamma_{3}");
  funct->SetParameter(0,  a_N1);
  funct->SetParameter(1,  a_params[0]);
  funct->SetParameter(2,  a_params[1]);
  funct->SetParameter(3,  a_params[2]);
  funct->SetParameter(4,  a_N2);
  funct->SetParameter(5,  a_params[4]);
  funct->SetParameter(6,  a_params[5]);
  funct->SetParameter(7,  a_params[6]);
  funct->SetParameter(8,  a_N3);
  funct->SetParameter(9,  a_params[8]);
  funct->SetParameter(10, a_params[9]);
  funct->SetParameter(11, a_params[10]);
  funct->SetNpx(300);
  return funct;
}

TF1* PhysMath::getThreeSkewFunct(string a_functName, double a_N1, double a_mu1, double a_sig1, double a_gam1, double a_N2,  double a_mu2, double a_sig2, double a_gam2, double a_N3,  double a_mu3, double a_sig3, double a_gam3){
  TF1* funct = new TF1(a_functName.c_str(),
    "[0]*PhysMath::skewNormal([1],[2],[3],x) + [4]*PhysMath::skewNormal([5],[6],[7],x) + [8]*PhysMath::skewNormal([9],[10],[11],x)",-10,10);
  //funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","#gamma_{1}","N_{2}","#mu_{2}","#sigma_{2}","#gamma_{2}","N_{3}","#mu_{3}","#sigma_{3}","#gamma_{3}");
  funct->SetParameter(0,  a_N1);
  funct->SetParameter(1,  a_mu1);
  funct->SetParameter(2,  a_sig1);
  funct->SetParameter(3,  a_gam1);
  funct->SetParameter(4,  a_N2);
  funct->SetParameter(5,  a_mu2);
  funct->SetParameter(6,  a_sig2);
  funct->SetParameter(7,  a_gam2);
  funct->SetParameter(8,  a_N3);
  funct->SetParameter(9,  a_mu3);
  funct->SetParameter(10, a_sig3);
  funct->SetParameter(11, a_gam3);
  funct->SetNpx(300);
  return funct;
}


TF1* PhysMath::getThreeStudentTFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3){
  TF1* funct = new TF1(a_functName.c_str(),
    "[0]*PhysMath::student_T([1],[2],[3],x) + [4]*PhysMath::student_T([5],[6],[7],x) + [8]*PhysMath::student_T([9],[10],[11],x)",-10,10);
  //funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","#nu_{1}","N_{2}","#mu_{2}","#sigma_{2}","#nu_{2}","N_{3}","#mu_{3}","#sigma_{3}","#nu_{3}");
  funct->SetParameter(0,  a_N1);
  funct->SetParameter(1,  a_params[0]);
  funct->SetParameter(2,  a_params[1]);
  funct->SetParameter(3,  a_params[2]);
  funct->SetParameter(4,  a_N2);
  funct->SetParameter(5,  a_params[4]);
  funct->SetParameter(6,  a_params[5]);
  funct->SetParameter(7,  a_params[6]);
  funct->SetParameter(8,  a_N3);
  funct->SetParameter(9,  a_params[8]);
  funct->SetParameter(10, a_params[9]);
  funct->SetParameter(11, a_params[10]);
  funct->SetNpx(300);
  return funct;
}

TF1* PhysMath::getThreeM2StudentTFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3, double a_momentum){

  //TF1* F1 = PhysMath::getM2StudentTFunction("", a_N1, a_params[0], a_params[1], a_params[2],  a_momentum);
  //TF1* F2 = PhysMath::getM2StudentTFunction("", a_N2, a_params[4], a_params[5], a_params[6],  a_momentum);
  //TF1* F3 = PhysMath::getM2StudentTFunction("", a_N3, a_params[8], a_params[9], a_params[10], a_momentum);
  TF1* funct = new TF1(a_functName.c_str(),
    "[0]*PhysMath::m2_student_T([1],[2],[3],x,[12]) + [4]*PhysMath::m2_student_T([5],[6],[7],x,[12]) + [8]*PhysMath::m2_student_T([9],[10],[11],x,[12])",-10,10);
  //funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","#nu_{1}","N_{2}","#mu_{2}","#sigma_{2}","#nu_{2}","N_{3}","#mu_{3}","#sigma_{3}","#nu_{3}");
  funct->SetParameter(0,  a_N1);
  funct->SetParameter(1,  a_params[0]);
  funct->SetParameter(2,  a_params[1]);
  funct->SetParameter(3,  a_params[2]);
  funct->SetParameter(4,  a_N2);
  funct->SetParameter(5,  a_params[4]);
  funct->SetParameter(6,  a_params[5]);
  funct->SetParameter(7,  a_params[6]);
  funct->SetParameter(8,  a_N3);
  funct->SetParameter(9,  a_params[8]);
  funct->SetParameter(10, a_params[9]);
  funct->SetParameter(11, a_params[10]);
  funct->SetParameter(12, a_momentum);
  //funct->SetParameter(13, F1->GetParameter(5));
  //funct->SetParameter(14, F2->GetParameter(5));
  //funct->SetParameter(15, F3->GetParameter(5));
  funct->SetNpx(300);

  //delete F1;
  //delete F2;
  //delete F3;
  return funct;
}


TF1* PhysMath::getThreePearsonIVFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3){
  TF1* funct = new TF1(a_functName.c_str(),
    "[0]*PhysMath::pearsonIV([1],[2],[3],[4],x,true) + [5]*PhysMath::pearsonIV([6],[7],[8],[9],x,true) + [10]*PhysMath::pearsonIV([11],[12],[13],[14],x,true)",-10,10);
  //funct->SetParNames("N_{1}","#lambda_{1}","#alpha_{1}","m_{1}","#nu_{1}","N_{2}","#lambda_{2}","#alpha_{2}","m_{2}","#nu_{2}","N_{3}","#lambda_{3}","#alpha_{3}","m_{3}","#nu_{3}");
  funct->SetParameter(0,  a_N1);
  funct->SetParameter(1,  a_params[0]);
  funct->SetParameter(2,  a_params[1]);
  funct->SetParameter(3,  a_params[2]);
  funct->SetParameter(4,  a_params[3]);
  funct->SetParameter(5,  a_N2);
  funct->SetParameter(6,  a_params[4]);
  funct->SetParameter(7,  a_params[5]);
  funct->SetParameter(8,  a_params[6]);
  funct->SetParameter(9,  a_params[7]);
  funct->SetParameter(10, a_N3);
  funct->SetParameter(11, a_params[8]);
  funct->SetParameter(12, a_params[9]);
  funct->SetParameter(13, a_params[10]);
  funct->SetParameter(14, a_params[11]);
  funct->SetNpx(300);
  return funct;
}

TF1* PhysMath::getThreeDoubleGaussFunct(string a_functName, double* a_params, double a_N1, double a_N2, double a_N3){
  TF1* funct = new TF1(a_functName.c_str(),
    "PhysMath::doubleGauss([0],[1],[2],[3],[4],x) + PhysMath::doubleGauss([5],[6],[7],[8],[9],x) + PhysMath::doubleGauss([10],[11],[12],[13],[14],x)",-10,10);
  //funct->SetParNames("N_{1}","#lambda_{1}","#alpha_{1}","m_{1}","#nu_{1}","N_{2}","#lambda_{2}","#alpha_{2}","m_{2}","#nu_{2}","N_{3}","#lambda_{3}","#alpha_{3}","m_{3}","#nu_{3}");
  funct->SetParameter(0,  a_N1);
  funct->SetParameter(1,  a_params[0]);
  funct->SetParameter(2,  a_params[1]);
  funct->SetParameter(3,  a_params[2]);
  funct->SetParameter(4,  a_params[3]);
  funct->SetParameter(5,  a_N2);
  funct->SetParameter(6,  a_params[4]);
  funct->SetParameter(7,  a_params[5]);
  funct->SetParameter(8,  a_params[6]);
  funct->SetParameter(9,  a_params[7]);
  funct->SetParameter(10, a_N3);
  funct->SetParameter(11, a_params[8]);
  funct->SetParameter(12, a_params[9]);
  funct->SetParameter(13, a_params[10]);
  funct->SetParameter(14, a_params[11]);
  funct->SetNpx(300);
  return funct;
}





vector< double > PhysMath::logParabToGaus(double a_scale, double a_mean, double a_lambda, double a_scaleError, double a_meanError, double a_lambdaError){
  vector< double > params(6,0.0);
  params[0] = pow(a_lambda,-0.5)*2.506628274631*exp(a_scale);
  params[1] = a_mean; // same in parabola and gaussian
  params[2] = 1.0/sqrt(a_lambda);
  params[3] = sqrt(  params[0]*params[0]*a_scaleError *a_scaleError
                     +  0.25*params[0]*params[0]/(a_lambda*a_lambda) * a_lambdaError * a_lambdaError ); //first order only
                     //+ 0.5*a_scaleError*a_scaleError
  params[4] = a_meanError; // same in parabola and gaussian
  params[5] = 0.5*pow(a_lambda,-1.5)*a_lambdaError;

  return params;

}

vector< double > PhysMath::logGausToParab(double a_N, double a_mu, double a_sig){
  vector< double > params(3,0.0);
  params[0] = log(a_N) - log(a_sig) - 0.918938533204672741;
  params[1] = a_mu; // same in parabola and gaussian
  params[2] = 1.0/(a_sig*a_sig);
  return params;
}


double PhysMath::studentT_var_to_stdDev(double a_sigma_studentT, double a_nu_studentT){
  return sqrt(a_sigma_studentT*a_sigma_studentT*a_nu_studentT/(a_nu_studentT - 2.0));
}

double PhysMath::studentT_var_from_stdDev(double a_stdDev_gaus, double a_nu_studentT){
  return sqrt(a_stdDev_gaus*a_stdDev_gaus*(a_nu_studentT - 2.0)/a_nu_studentT);
}

double PhysMath::beta_to_betagamma(double a_beta){
  return a_beta*pow(1.0-a_beta*a_beta,-0.5);
}
double PhysMath::betagamma(double a_p, double a_m){
  return a_p/a_m;
  //return beta_to_betagamma( a_p/sqrt(a_m*a_m + a_p*a_p) );
}
double PhysMath::gamma(double a_p, double a_m){
  return sqrt(a_m*a_m + a_p*a_p)/a_m;
}
double PhysMath::beta(double a_p, double a_m){
  return a_p/sqrt(a_m*a_m + a_p*a_p);
}

//___________________________________________________________________
double PhysMath::blastWaveModelFit(double *x, double *par){

  double mTm0=x[0]; //mT-m0
  double Amp  = par[0]; //Amplitude
  double T     = par[1]; //Tkin
  double betaS = par[2]; //Surface Velocity
  double n     = par[3]; //Power in the Tranverse Velocity factor (betaT)
  double R     = par[4]; //Upper Limit of Radial Integral (SHOULD BE FIXED)
  double mass  = par[5]; //Should be FIXED


  //######  Radial Integral for the Blast Wave Model  ########
  //Integration Variable
  double r = 0.0; //The lower limit of integration is 0.0

  //Integration Steps
  double sum = 0.0;
  double stepWidth = 0.01;
  int    nSteps = (int)(fabs(R)/stepWidth);

  //Compute mT and pT
  double mT = mTm0 + mass;
  double pT = PhysMath::pT_from_mTm0(mTm0, mass);

  //Do the numerical integral
  double rho = 0;
  double BesselI0 = 0;
  double BesselK1 = 0;
  for (int iStep=0; iStep<nSteps; iStep++){
    rho   = TMath::ATanH(betaS * pow(r/R,n));
    BesselI0 = TMath::BesselI0((pT/T)*TMath::SinH(rho));
    BesselK1 = TMath::BesselK1((mT/T)*TMath::CosH(rho));
    sum += r * mT * BesselI0 * BesselK1 * stepWidth;
    r += stepWidth;
  }//End Loop Over iStep

  return (1.0/(2.0*TMath::Pi())) * Amp * sum / T;
}

//___________________________________________________________________
double PhysMath::pTblastWaveModelFit(double *x, double *par){

  double pT=x[0]; //mT-m0
  double Amp  = par[0]; //Amplitude
  double T     = par[1]; //Tkin
  double betaS = par[2]; //Surface Velocity
  double n     = par[3]; //Power in the Tranverse Velocity factor (betaT)
  double R     = par[4]; //Upper Limit of Radial Integral (SHOULD BE FIXED)
  double mass  = par[5]; //Should be FIXED


  //######  Radial Integral for the Blast Wave Model  ########
  //Integration Variable
  double r = 0.0; //The lower limit of integration is 0.0

  //Integration Steps
  double sum = 0.0;
  double stepWidth = 0.01;
  int    nSteps = (int)(fabs(R)/stepWidth);

  //Compute mT and pT
  double mT = PhysMath::mTm0(mass,pT) + mass;
  //double mT = mTm0 + mass;
  //double pT = PhysMath::pT_from_mTm0(mTm0, mass);

  //Do the numerical integral
  double rho = 0;
  double BesselI0 = 0;
  double BesselK1 = 0;
  for (int iStep=0; iStep<nSteps; iStep++){
    rho   = TMath::ATanH(betaS * pow(r/R,n));
    BesselI0 = TMath::BesselI0((pT/T)*TMath::SinH(rho));
    BesselK1 = TMath::BesselK1((mT/T)*TMath::CosH(rho));
    sum += r * mT * BesselI0 * BesselK1 * stepWidth;
    r += stepWidth;
  }//End Loop Over iStep

  return (1.0/(2.0*TMath::Pi())) * Amp * sum / T;
}
//________________________________________________


double PhysMath::legendre(int a_order,double a_x){
  if(a_x < -1.0) return 0.0;
  if(a_x > 1.0)  return 0.0;
  //if(a_x < -1.0) a_x = -1.0;
  //if(a_x > 1.0)  a_x = 1.0;

  switch(a_order){
      case 0 : return 1.0;
      case 1 : return a_x;
      case 2 : return 0.5*(3.0*a_x*a_x-1.0);
      case 3 : return 0.5*(5.0*a_x*a_x*a_x-3.0*a_x);
      case 4 : return (35.0*a_x*a_x*a_x*a_x-30.0*a_x*a_x+3.0)/8.0;
      case 5 : return (63.0*pow(a_x,5)-70.0*a_x*a_x*a_x+15.0*a_x)/8.0;
      case 6 : return (231.0*pow(a_x,6)-315.0*pow(a_x,4)+105.0*a_x*a_x-5.0)/16.0;
      case 7 : return (429.0*pow(a_x,7)-693.0*pow(a_x,5)+315.0*a_x*a_x*a_x-35.0*a_x)/16.0;
      case 8 : return (35.0-1260.0*a_x*a_x+6930.0*a_x*a_x*a_x*a_x-12012.0*pow(a_x,6)+6435.0*pow(a_x,8))/128.0;
      case 9 : return (315.0*a_x-4620.0*a_x*a_x*a_x+18018.0*pow(a_x,5)-25740.0*pow(a_x,7)+12155.0*pow(a_x,9))/128.0;
      case 10 : return (-63.0+3465.0*a_x*a_x-30030.0*pow(a_x,4)+90090.0*pow(a_x,6)-109395.0*pow(a_x,8)+46189.0*pow(a_x,10))/256.0;
      case 11 : return (-693.0*a_x+15015.0*a_x*a_x*a_x-90090.0*pow(a_x,5)+218790.0*pow(a_x,7)-230945.0*pow(a_x,9)+88179.0*pow(a_x,11))/256.0;
      case 12 : return (231.0-18018.0*a_x*a_x+225225.0*pow(a_x,4)-1021020.0*pow(a_x,6)+2078505.0*pow(a_x,8)-1939938.0*pow(a_x,10)+676039.0*pow(a_x,12))/1024.0;
      case 13 : return (3003.0*a_x-90090.0*a_x*a_x*a_x+765765.0*pow(a_x,5)-2771340.0*pow(a_x,7)+4849845.0*pow(a_x,9)-4056234.0*pow(a_x,11)+1300075.0*pow(a_x,13))/1024.0;
      case 14 : return (-429.0+45045.0*a_x*a_x-765765.0*pow(a_x,4)+4849845.0*pow(a_x,6)-14549535.0*pow(a_x,8)+22309287.0*pow(a_x,10)-16900975.0*pow(a_x,12)+5014575.0*pow(a_x,14))/2048.0;
      case 15 : return (-6435.0*a_x+255255.0*a_x*a_x*a_x-2909907.0*pow(a_x,5)+14549535.0*pow(a_x,7)-37182145.0*pow(a_x,9)+50702925.0*pow(a_x,11)-35102025.0*pow(a_x,13)+9694845.0*pow(a_x,15))/2048.0;
      case 16 : return (6435.0-875160.0*a_x*a_x+19399380.0*pow(a_x,4)-162954792.0*pow(a_x,6)+669278610.0*pow(a_x,8)-1487285800.0*pow(a_x,10)+1825305300.0*pow(a_x,12)-1163381400.0*pow(a_x,14)+300540195.0*pow(a_x,16))/32768.0;
      case 17 : return (109395.0*a_x-5542680.0*a_x*a_x*a_x+81477396.0*pow(a_x,5)-535422888.0*pow(a_x,7)+1859107250.0*pow(a_x,9)-3650610600.0*pow(a_x,11)+4071834900.0*pow(a_x,13)-2404321560.0*pow(a_x,15)+583401555.0*pow(a_x,17))/32768.0;
      case 18 : return (-12155.0+2078505.0*a_x*a_x-58198140.0*pow(a_x,4)+624660036.0*pow(a_x,6)-3346393050.0*pow(a_x,8)+10039179150.0*pow(a_x,10)-17644617900.0*pow(a_x,12)+18032411700.0*pow(a_x,14)-9917826435.0*pow(a_x,16)+2268783825.0*pow(a_x,18))/65536.0;
      case 19 : return (-230945.0*a_x+14549535.0*a_x*a_x*a_x-267711444.0*pow(a_x,5)+2230928700.0*pow(a_x,7)-10039179150.0*pow(a_x,9)+26466926850.0*pow(a_x,11)-42075627300.0*pow(a_x,13)+39671305740.0*pow(a_x,15)-20419054425.0*pow(a_x,17)+4418157975.0*pow(a_x,19))/65536.0;
      case 20 : return (46189.0-9699690.0*a_x*a_x+334639305.0*pow(a_x,4)-4461857400.0*pow(a_x,6 )+ 30117537450.0*pow(a_x,8 )- 116454478140.0*pow(a_x,10) + 273491577450.0*pow(a_x,12) - 396713057400.0*pow(a_x,14) + 347123925225.0*pow(a_x,16) - 167890003050.0*pow(a_x,18) + 34461632205.0*pow(a_x,20))/262144.0;
      case 21 : return (969969.0*a_x-74364290.0*pow(a_x,3)+1673196525.0*pow(a_x,5) - 17210021400.0*pow(a_x,7) + 97045398450.0*pow(a_x,9) - 328189892940.0*pow(a_x,11) + 694247850450.0*pow(a_x,13) - 925663800600.0*pow(a_x,15) + 755505013725.0*pow(a_x,17) - 344616322050.0*pow(a_x,19) + 67282234305.0*pow(a_x,21))/262144.0;
      case 22 : return (-88179.0+22309287.0*a_x*a_x-929553625.0*pow(a_x,4)+15058768725.0*pow(a_x,6) - 124772655150.0*pow(a_x,8) + 601681470390.0*pow(a_x,10) - 1805044411170.0*pow(a_x,12) + 3471239252250.0*pow(a_x,14) - 4281195077775.0*pow(a_x,16) + 3273855059475.0*pow(a_x,18) - 1412926920405.0*pow(a_x,20) + 263012370465.0*pow(a_x,22))/524288.0;
      case 23 : return (-2028117.0*a_x+185910725.0*pow(a_x,3)-5019589575.0*pow(a_x,5 )+ 62386327575.0*pow(a_x,7 )- 429772478850.0*pow(a_x,9 )+ 1805044411170.0*pow(a_x,11) - 4859734953150.0*pow(a_x,13) + 8562390155550.0*pow(a_x,15) - 9821565178425.0*pow(a_x,17) + 7064634602025.0*pow(a_x,19) - 2893136075115.0*pow(a_x,21) + 514589420475.0*pow(a_x,23))/524288.0;
      case 24 : return (676039.0-202811700.0*a_x*a_x+10039179150.0*pow(a_x,4)-194090796900.0*pow(a_x,6) + 1933976154825.0*pow(a_x,8 )- 11345993441640.0*pow(a_x,10) + 42117702927300.0*pow(a_x,12) - 102748681866600.0*pow(a_x,14) + 166966608033225.0*pow(a_x,16) - 178970743251300.0*pow(a_x,18) + 121511715154830.0*pow(a_x,20) - 47342226683700.0*pow(a_x,22) + 8061900920775.0*pow(a_x,24))/4194304.0;
      case 25 : return (16900975.0*a_x-1825305300.0*pow(a_x,3)+58227239070.0*pow(a_x,5 )- 859544957700.0*pow(a_x,7 )+ 7091245901025.0*pow(a_x,9 )- 36100888223400.0*pow(a_x,11) + 119873462177700.0*pow(a_x,13) - 267146572853160.0*pow(a_x,15) + 402684172315425.0*pow(a_x,17) - 405039050516100.0*pow(a_x,19) + 260382246760350.0*pow(a_x,21) - 96742811049300.0*pow(a_x,23) + 15801325804719.0*pow(a_x,25))/4194304.0;
      case 26 : return (-1300075.0+456326325.0*a_x*a_x-26466926850.0*pow(a_x,4)+601681470390.0*pow(a_x,6 )- 7091245901025.0*pow(a_x,8 )+ 49638721307175.0*pow(a_x,10) - 222622144044300.0*pow(a_x,12) + 667866432132900.0*pow(a_x,14) - 1369126185872445.0*pow(a_x,16) + 1923935489951475.0*pow(a_x,18) - 1822675727322450.0*pow(a_x,20) + 1112542327066950.0*pow(a_x,22) - 395033145117975.0*pow(a_x,24) + 61989816618513.0*pow(a_x,26))/8388608.0;
      case 27 : return (-35102025.0*a_x+4411154475.0*pow(a_x,3)-164094946470.0*pow(a_x,5 )+ 2836498360410.0*pow(a_x,7 )- 27577067392875.0*pow(a_x,9 )+ 166966608033225.0*pow(a_x,11) - 667866432132900.0*pow(a_x,13) + 1825501581163260.0*pow(a_x,15) - 3463083881912655.0*pow(a_x,17) + 4556689318306125.0*pow(a_x,19) - 4079321865912150.0*pow(a_x,21) + 2370198870707850.0*pow(a_x,23) - 805867616040669.0*pow(a_x,25) + 121683714103007.0*pow(a_x,27))/8388608.0;
      case 28 : return (5014575.0-2035917450.0*a_x*a_x+136745788725.0*pow(a_x,4)-3610088822340.0*pow(a_x,6 )+ 49638721307175.0*pow(a_x,8 )- 408140597414550.0*pow(a_x,10) + 2170565904431925.0*pow(a_x,12) - 7823578204985400.0*pow(a_x,14) + 19624141997505045.0*pow(a_x,16) - 34630838819126550.0*pow(a_x,18) + 42832879592077575.0*pow(a_x,20) - 36343049350853700.0*pow(a_x,22) + 20146690401016725.0*pow(a_x,24) - 6570920561562378.0*pow(a_x,26) + 956086325095055.0*pow(a_x,28))/33554432.0;
      case 29 : return (145422675.0*a_x-21037813650.0*pow(a_x,3)+902522205585.0*pow(a_x,5 )- 18050444111700.0*pow(a_x,7 )+ 204070298707275.0*pow(a_x,9 )- 1447043936287950.0*pow(a_x,11) + 6845630929362225.0*pow(a_x,13) - 22427590854291480.0*pow(a_x,15) + 51946258228689825.0*pow(a_x,17) - 85665759184155150.0*pow(a_x,19) + 99943385714847675.0*pow(a_x,21) - 80586761604066900.0*pow(a_x,23) + 42710983650155457.0*pow(a_x,25) - 13385208551330770.0*pow(a_x,27) + 1879204156221315.0*pow(a_x,29))/33554432.0;


      default :
         cout << "ERROR: Invalid Order in PhysMath::legendre >29 or <0 : " << a_order  << endl;
   }
   return 0.0;
}


double PhysMath::semiHermite(int a_order, double a_x){
  switch(a_order){
      case 0  : return 1.0;
      case 1  : return 1.0;
      case 2  : return 1.0;
      case 3  : return 1.0;
      case 4  : return 1.0;
      case 5  : return 1.0;
      case 6  : return 1.0;
      case 7  : return 1.0;
      case 8  : return 1.0;
      case 9  : return 1.0;
      case 10 : return 1.0;
      case 11 : return 1.0;
      case 12 : return 1.0;
      case 13 : return 1.0;
      case 14 : return 1.0;
      case 15 : return 1.0;
      case 16 : return 1.0;
      case 17 : return 1.0;
      case 18 : return 1.0;
      case 19 : return 1.0;
      case 20 : return 1.0;
      case 21 : return 1.0;
      case 22 : return 1.0;
      case 23 : return 1.0;
      case 24 : return 1.0;
      case 25 : return 1.0;
      case 26 : return 1.0;
      case 27 : return 1.0;
      case 28 : return 1.0;
      case 29 : return 1.0;


      default :
         cout << "ERROR: Invalid Order in PhysMath::legendre >19 or <0 : " << a_order  << endl;
   }
   return 0.0;
}
//________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

double PhysMath::tsallis(double a_dNdy,double a_T,double a_q_minus_one,double a_m, double a_mTm0){
  return a_dNdy*(1.0-a_q_minus_one)*pow(1.0 + a_mTm0*a_q_minus_one/(a_T+a_m*a_q_minus_one),-(a_q_minus_one+1.0)/a_q_minus_one)/((a_T+a_m)*(a_T+a_m*a_q_minus_one)*6.2831853072);
}

//________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________




double PhysMath::integralWithMtScale(TF1* a_funct,double a_mass,double a_stepSize, double a_lowbound,double a_highbound, bool a_extrapolateToInf){
  double integral = 0.0;
  for(double evalPt = a_lowbound; evalPt < a_highbound - a_stepSize; evalPt += a_stepSize){
    integral += 0.5*6.2831853072*(a_funct->Eval(evalPt) * (evalPt + a_mass) + a_funct->Eval(evalPt + a_stepSize) * (evalPt + a_stepSize + a_mass));
  }
  integral = integral*a_stepSize;
  if(a_extrapolateToInf){
    double x1 = a_highbound - 5.0*a_stepSize;
    double x2 = a_highbound;
    double y1 = a_funct->Eval(x1) * 6.2831853072 * (x1 + a_mass);
    double y2 = a_funct->Eval(x2) * 6.2831853072 * (x2 + a_mass);
    double slope = (y2 - y1)/(x2 - x1);
    double intercept = y2 - slope*x2;
    double xCrossing = -intercept/slope;
    double area = 0.5*(xCrossing - x2)*y2;
    if(area > 0.0) integral += area;
  }
  return integral;
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


bool PhysMath::isTMatrixDSymEmpty(TMatrixDSym a_matrix){
  for(int row = 0; row < a_matrix.GetNrows(); row++){
    for(int col = row; col < a_matrix.GetNcols(); col++){
      if(a_matrix[row][col] != 0.0){
        return false;
      }
    }
  }
  return true;
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


TMatrixDSym PhysMath::reduceCovMatrix(TMatrixDSym a_covMatrix){
  int newSize = 0;
  for(int iii = 0; iii < a_covMatrix.GetNrows(); iii++){
    if(a_covMatrix[iii][iii] != 0.0){
      newSize++;
    }
  }

  TMatrixDSym redMatrix(newSize);

  int counterRow = 0;
  for(int row = 0; row < a_covMatrix.GetNrows(); row++){
    if(a_covMatrix[row][row] == 0.0) continue;
    int counterCol = 0;
    for(int col = 0; col < a_covMatrix.GetNcols(); col++){
      if(a_covMatrix[col][col] == 0.0) continue;

      redMatrix[counterRow][counterCol] = a_covMatrix[row][col];
      counterCol++;
    }
    counterRow++;
  }

  return redMatrix;

}


//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


TMatrixD PhysMath::reduceCovMatrix(TMatrixD a_covMatrix){
  int newSize = 0;
  for(int iii = 0; iii < a_covMatrix.GetNrows(); iii++){
    if(a_covMatrix[iii][iii] != 0.0){
      newSize++;
    }
  }

  TMatrixD redMatrix(newSize,newSize);

  int counterRow = 0;
  for(int row = 0; row < a_covMatrix.GetNrows(); row++){
    if(a_covMatrix[row][row] == 0.0) continue;
    int counterCol = 0;
    for(int col = 0; col < a_covMatrix.GetNcols(); col++){
      if(a_covMatrix[col][col] == 0.0) continue;

      redMatrix[counterRow][counterCol] = a_covMatrix[row][col];
      counterCol++;
    }
    counterRow++;
  }

  return redMatrix;

}


//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


TMatrixDSym PhysMath::reduceCovMatrix_keepPassedVector(TMatrixDSym a_covMatrix, vector< int > a_rowsToKeep){
  int newSize = (int) a_rowsToKeep.size();

  cout << "Reducing covMatrix" << endl;

  TMatrixDSym nullMatrix;

  if(newSize == 0){
    cout << "no rows to keep were passed" << endl;
    return nullMatrix;
  }

  TMatrixDSym redMatrix(newSize);

  for(int row = 0; row < newSize; row++){
    int rowVal = a_rowsToKeep[row];
    for(int col = 0; col < newSize; col++){
      int colVal = a_rowsToKeep[col];
      if(rowVal == -999 || colVal == -999) redMatrix[row][col] = 0.0;
      else redMatrix[row][col] = a_covMatrix[rowVal][colVal];
    }
  }

  return redMatrix;

}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


TMatrixDSym PhysMath::addZerosToCovMatrix(TMatrixDSym a_covMatrix, int a_numExtraZeroRows){
  int originalSize = a_covMatrix.GetNrows();
  int newSize = originalSize + a_numExtraZeroRows;
  TMatrixDSym bigMatrix(newSize);
  for(int row = 0; row < newSize; row++){
    for(int col = 0; col < newSize; col++){
      if(row >= originalSize || col >= originalSize) bigMatrix[row][col] = 0.0;
      else bigMatrix[row][col] = a_covMatrix[row][col];
    }
  }
  #ifdef _PHYS_MATH_DEBUG_
    cout << "PhysMath::addZerosToCovMatrix with " << a_numExtraZeroRows << " extra zero rows." << endl;
    cout << "     Original Size: " << originalSize << "   New Size: " << newSize << "   actual size of new matrix: " << bigMatrix.GetNrows() << endl;
  #endif
  return bigMatrix;
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

TH2D* PhysMath::makeCorrelationHistogram(string a_name, string a_title, TMatrixDSym a_covMatrix, vector< string > a_paramNames){
  int numParams = a_covMatrix.GetNrows();
  TH2D* correlation_histo = new TH2D(a_name.c_str(),a_title.c_str(),numParams,0.0,numParams,numParams,0.0,numParams);
  correlation_histo->SetStats(kFALSE);
  #ifdef _PHYS_MATH_DEBUG_
    cout << " PhysMath::makeCorrelationHistogram has " << numParams << " parameters from the covMatrix, but " << a_paramNames.size() << " names " << endl;
    for(int paramIndex = 0; paramIndex < a_paramNames.size(); paramIndex++){
      cout << " param: " << paramIndex << "  " << a_paramNames[paramIndex] << endl;
    }
  #endif

  //label the bins if the user supplied the paramter names
  if(numParams <= a_paramNames.size()){
    for(int paramIndex = 0; paramIndex < numParams; paramIndex++){
      correlation_histo->GetXaxis()->SetBinLabel(paramIndex + 1,a_paramNames[paramIndex].c_str());
      correlation_histo->GetYaxis()->SetBinLabel(paramIndex + 1,a_paramNames[paramIndex].c_str());
    }
  }
  for(int binX = 1; binX <= correlation_histo->GetNbinsX(); binX++){
    for(int binY = 1; binY <= correlation_histo->GetNbinsY(); binY++){
      if(binX == binY){
        correlation_histo->SetBinContent(binX,binY,1.0);
      }else{
        correlation_histo->SetBinContent(binX,binY, a_covMatrix[binX - 1][binY - 1]/(sqrt(a_covMatrix[binX - 1][binX - 1])*sqrt(a_covMatrix[binY - 1][binY - 1])));
      }
      if(!std::isfinite(correlation_histo->GetBinContent(binX,binY))){
        correlation_histo->SetBinContent(binX,binY,0.0);
      }
    }
  }
  return correlation_histo;
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

TH2D* PhysMath::makeCovarianceHistogram(string a_name, string a_title, TMatrixDSym a_covMatrix, vector< string > a_paramNames){
  int numParams = a_covMatrix.GetNrows();
  TH2D* covarianceHisto = new TH2D(a_name.c_str(),a_title.c_str(),numParams,0.0,numParams,numParams,0.0,numParams);
  covarianceHisto->SetStats(kFALSE);

  //label the bins if the user supplied the paramter names
  if(numParams <= a_paramNames.size()){
    for(int paramIndex = 0; paramIndex < numParams; paramIndex++){
      covarianceHisto->GetXaxis()->SetBinLabel(paramIndex + 1,a_paramNames[paramIndex].c_str());
      covarianceHisto->GetYaxis()->SetBinLabel(paramIndex + 1,a_paramNames[paramIndex].c_str());
    }
  }
  for(int binX = 1; binX <= covarianceHisto->GetNbinsX(); binX++){
    for(int binY = 1; binY <= covarianceHisto->GetNbinsY(); binY++){
      covarianceHisto->SetBinContent(binX,binY, a_covMatrix[binX - 1][binY - 1]);
    }
  }
  return covarianceHisto;
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
double PhysMath::interpolate(double a_x, double a_x1, double a_y1, double a_x2,double a_y2){
  return (a_y2 - a_y1)*(a_x - a_x1)/(a_x2 - a_x1) + a_y1;
}

//__________________________________________________________________________________________________________________________________________________
void PhysMath::interpolate(double a_x, double a_x1, double a_y1, double a_yErr1, double a_x2, double a_y2, double a_yErr2, double &a_yValue, double &a_yError){
  double A = (a_x - a_x1)/(a_x2 - a_x1);
  double B = (a_x2 - a_x)/(a_x2 - a_x1);
  a_yValue = A*(a_y2 - a_y1) + a_y1;
  a_yError = sqrt(A*A*a_yErr2*a_yErr2 + B*B*a_yErr1*a_yErr1);
}

void PhysMath::interpolateAsymmError(double a_x, double a_x1, double a_y1, double a_yErrLow1, double a_yErrHigh1, double a_x2, double a_y2, double a_yErrLow2, double a_yErrHigh2, double &a_yValue, double &a_yErrLow, double &a_yErrHigh){
  double A = (a_x - a_x1)/(a_x2 - a_x1);
  double B = (a_x2 - a_x)/(a_x2 - a_x1);
  a_yValue = A*(a_y2 - a_y1) + a_y1;
  a_yErrLow  = sqrt(A*A*a_yErrLow2*a_yErrLow2 + B*B*a_yErrLow1*a_yErrLow1);
  a_yErrHigh = sqrt(A*A*a_yErrHigh2*a_yErrHigh2 + B*B*a_yErrHigh1*a_yErrHigh1);
  // maybe doulbe check math here
}







//__________________________________________________________________________________________________________________________________________________

double PhysMath::interpolateSquare2D(double a_x, double a_y, double a_x1, double a_y1, double a_x2,double a_y2, double a_z1, double a_z2, double a_z3, double a_z4){
  #ifdef _PHYSMATH_DEBUG_
    cout << "PhysMath::interpolateSquare2D(" << a_x << "," << a_y << "," << a_x1 << "," << a_y1 << "," << a_x2 << "," << a_y2 << "," << a_z1 << "," << a_z2 << "," << a_z3 << "," << a_z4 << ")" << endl;
  #endif
  double x_width = fabs(a_x2 - a_x1);
  double y_width = fabs(a_y2 - a_y1);
  double z_xLow  = (1.0 - fabs(a_y - a_y1)/y_width)*a_z1   + (1.0 - fabs(a_y2 - a_y)/y_width)*a_z3;
  double z_xHigh = (1.0 - fabs(a_y - a_y1)/y_width)*a_z2   + (1.0 - fabs(a_y2 - a_y)/y_width)*a_z4;
  double z_x     = (1.0 - fabs(a_x - a_x1)/x_width)*z_xLow + (1.0 - fabs(a_x2 - a_x)/x_width)*z_xHigh;
  return z_x;
  //double z_yLow  = (1.0 - (a_x - a_x1)/x_width)*a_z1   + (1.0 - (a_x2 - a_x)/x_width)*a_z2;
  //double z_yHigh = (1.0 - (a_x - a_x1)/x_width)*a_z3   + (1.0 - (a_x2 - a_x)/x_width)*a_z4;
  //double z_y     = (1.0 - (a_y - a_y1)/y_width)*z_yLow + (1.0 - (a_y2 - a_y)/y_width)*z_yHigh;
  //return 0.5*(z_x + z_y);
}

double PhysMath::signed_sqrt(double a_number){
  if(a_number < 0.0) return -sqrt(-a_number);
  else               return  sqrt(a_number);
}

double PhysMath::signed_square(double a_number){
  if(a_number < 0.0) return -a_number*a_number;
  else               return  a_number*a_number;
}

double PhysMath::addInQuadriture(double a_num1, double a_num2){return sqrt(a_num1*a_num1 + a_num2*a_num2);};




vector<double> PhysMath::calcMassSquaredEdges(int a_numBins, double a_timingResolution, double a_momToSample, double a_lowestParticleMass, double a_highestParticleMass){
        double timingLightDistance = a_timingResolution * 3.0E8;

	//get particle m^2 - 3sigma at mom = 2.0 GeV/c
	double momForMaxSigma = 2.0;
	double lowestMassSigma = timingLightDistance*sqrt(a_lowestParticleMass*a_lowestParticleMass + momForMaxSigma*momForMaxSigma)*momForMaxSigma;
        double lowM2 = a_lowestParticleMass*a_lowestParticleMass - 3.0*lowestMassSigma;

        double highestMassSigma = timingLightDistance*sqrt(a_highestParticleMass*a_highestParticleMass + momForMaxSigma*momForMaxSigma)*momForMaxSigma;
        double highM2 = a_highestParticleMass*a_highestParticleMass + 3.0*highestMassSigma;

        double bottomSampleVal = timingLightDistance * a_momToSample * sqrt(lowM2 + a_momToSample*a_momToSample);
        if(lowM2 + a_momToSample*a_momToSample < 0.0) bottomSampleVal = 0.0;
        double topSampleVal = timingLightDistance * a_momToSample * sqrt(highM2 + a_momToSample*a_momToSample);
        double sampleWidth = (topSampleVal - bottomSampleVal)/((double) a_numBins);
        vector<double> edges;
        for(int edgeIndex = 0; edgeIndex <= a_numBins; edgeIndex++){
          double edgeValue = pow((bottomSampleVal + sampleWidth*((double) edgeIndex) )/(timingLightDistance * a_momToSample),2.0) - a_momToSample*a_momToSample;
          edges.push_back(edgeValue);
        }
        return edges;
}


//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double PhysMath::fractionRemainingAfterHorizontalDistance(double a_rapidity, double a_mTm0, double a_mass, double a_tau, double a_horizontalDistance){
  double theta    = PhysMath::theta_from_eta(PhysMath::eta_from_y_mTm0(a_rapidity, a_mTm0, a_mass));
  double momentum = PhysMath::pTot_from_y_mTm0(a_rapidity,a_mTm0,a_mass);
  double distance = fabs(a_horizontalDistance/cos(theta));
  return exp(-a_mass*distance/(momentum*2.99792e8*a_tau));
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double PhysMath::fractionRemainingAfterRadialDistance(double a_rapidity, double a_mTm0, double a_mass, double a_tau, double a_radialDistance){
  double theta    = PhysMath::theta_from_eta(PhysMath::eta_from_y_mTm0(a_rapidity, a_mTm0, a_mass));
  double momentum = PhysMath::pTot_from_y_mTm0(a_rapidity,a_mTm0,a_mass);
  double distance = fabs(a_radialDistance/sin(theta));
  return exp(-a_mass*distance/(momentum*2.99792e8*a_tau));
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double PhysMath::fractionRemainingAfterDistance(double a_rapidity, double a_mTm0, double a_mass, double a_tau, double a_radialDistance, double a_horizontalDistance){
  double theta        = fabs(PhysMath::theta_from_eta(PhysMath::eta_from_y_mTm0(a_rapidity, a_mTm0, a_mass)));
  double theta_corner = fabs(atan(a_radialDistance/a_horizontalDistance));
  double momentum = PhysMath::pTot_from_y_mTm0(a_rapidity,a_mTm0,a_mass);
  double distance = 0;
  if(theta < theta_corner) distance = fabs(a_horizontalDistance/cos(theta));
  else                     distance = fabs(a_radialDistance/sin(theta));
  return exp(-a_mass*distance/(momentum*2.99792e8*a_tau));
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double PhysMath::fractionRemainingAfterHorizontalDistance_EtaPt(double a_eta, double a_pT, double a_mass, double a_tau, double a_horizontalDistance){
  double momentum = PhysMath::pTot_from_eta_pT(a_eta,a_pT);
  double distance = fabs(a_horizontalDistance/cos(2.0*atan(exp(-a_eta))));
  return exp(-a_mass*distance/(momentum*2.99792e8*a_tau));
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double PhysMath::fractionRemainingAfterRadialDistance_EtaPt(double a_eta, double a_pT, double a_mass, double a_tau, double a_radialDistance){
  double momentum = PhysMath::pTot_from_eta_pT(a_eta,a_pT);
  double theta    = fabs(PhysMath::theta_from_eta(a_eta));
  double distance = fabs(a_radialDistance/sin(theta));
  return exp(-a_mass*distance/(momentum*2.99792e8*a_tau));
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double PhysMath::fractionRemainingAfterDistance_EtaPt(double a_eta, double a_pT, double a_mass, double a_tau, double a_radialDistance, double a_horizontalDistance){
  double theta        = fabs(PhysMath::theta_from_eta(a_eta));
  double theta_corner = fabs(atan(a_radialDistance/a_horizontalDistance));
  double momentum = PhysMath::pTot_from_eta_pT(a_eta,a_pT);
  double distance = 0;
  if(theta < theta_corner) distance = fabs(a_horizontalDistance/cos(theta));
  else                     distance = fabs(a_radialDistance/sin(theta));
  //return exp(-a_mass*distance/(momentum*2.99792e8*a_tau));
  return fractionRemainingAfterDistance(momentum, a_mass, 2.99792e8*a_tau, distance);

}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double PhysMath::fractionRemainingAfterDistance(double a_pTot, double a_mass, double a_c_tau, double a_distance){
  return exp(-a_mass*a_distance/(a_pTot*a_c_tau));
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double PhysMath::calc_pct_kaon_left(double a_rapidity, double a_mTm0, double a_mass){
  double kaonEta = PhysMath::eta_from_y_mTm0(a_rapidity, a_mTm0, a_mass);
  double kaonMom = PhysMath::pTot_from_y_mTm0(a_rapidity, a_mTm0, a_mass);
  double kaonBetaGamma = kaonMom/a_mass;
  //double kaonBeta = sqrt((kaonBetaGamma*kaonBetaGamma)/((kaonBetaGamma*kaonBetaGamma)+1.0));
  //double kaonGamma = 1.0/sqrt(1.0-(kaonBeta*kaonBeta));
  //double kaonTheta = 2.0*atan(exp(-kaonEta));
  double kaonTheta = PhysMath::theta_from_eta(kaonEta);
  double kaonBeta  = PhysMath::beta(kaonMom, a_mass);
  double kaonGamma = kaonBetaGamma/kaonBeta;

  double kaonV = (2.99792e8)*kaonBeta;
  double kaonDist = fabs(0.10/sin(kaonTheta));
  double kaonTime = kaonDist/kaonV;
  double kaonTau = 1.238e-8;
  double pctKaonLeft = exp(-kaonTime/(kaonGamma*kaonTau));

  #ifdef _PHYS_MATH_DEBUG_
    cout << "kaon eta: " << kaonEta << endl;
    cout << "kaon momentum: " << kaonMom << endl;
    cout << "kaon beta*gamma: " << kaonBetaGamma << endl;
    cout << "kaon beta: " << kaonBeta << endl;
    cout << "kaon velocity: " << kaonV << endl;
    cout << "kaon gamma: " << kaonGamma << endl;
    cout << "kaon Theta: " << kaonTheta << endl;
    cout << "kaon dist: " << kaonDist << endl;
    cout << "kaon time: " << kaonTime << endl;
    cout << "percent of kaons left: " << pctKaonLeft << endl;
  #endif

  return pctKaonLeft;

}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double PhysMath::dca_from_decay(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect){
  if(!a_daughter_vect) return -999;
  double displacement[3] = {a_secondaryVertex->X() - a_primaryVertex->X(), a_secondaryVertex->Y() - a_primaryVertex->Y(), a_secondaryVertex->Z() - a_primaryVertex->Z()};
  double distance = sqrt(displacement[0]*displacement[0] + displacement[1]*displacement[1] + displacement[2]*displacement[2]);
  double daughter_pTot = a_daughter_vect->P();
  double parent_direction_unit[3]   = {displacement[0]/distance,displacement[1]/distance,displacement[2]/distance};
  double daughter_direction_unit[3] = {a_daughter_vect->Px()/daughter_pTot,a_daughter_vect->Py()/daughter_pTot,a_daughter_vect->Pz()/daughter_pTot};
  return fabs(distance * crossProduct3D(daughter_direction_unit,parent_direction_unit));
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


void PhysMath::dca_from_decay(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect, double &a_dca_x, double &a_dca_y, double &a_dca_z){
  if(!a_daughter_vect) return;
  double displacement[3] = {a_secondaryVertex->X() - a_primaryVertex->X(), a_secondaryVertex->Y() - a_primaryVertex->Y(), a_secondaryVertex->Z() - a_primaryVertex->Z()};
  double distance = sqrt(displacement[0]*displacement[0] + displacement[1]*displacement[1] + displacement[2]*displacement[2]);
  double daughter_pTot = a_daughter_vect->P();
  //cout << "daughter mom check : " << daughter_pTot << "   =?  " << sqrt(a_daughter_vect->Px()*a_daughter_vect->Px() + a_daughter_vect->Py()*a_daughter_vect->Py() + a_daughter_vect->Pz()*a_daughter_vect->Pz()) << endl;
  //cout << "Proton Momentum: " << daughter_pTot <<  "   displacement of parent: " << distance << endl;
  double parent_direction_unit[3]   = {displacement[0]/distance,displacement[1]/distance,displacement[2]/distance};
  //cout << "  unit check = " << sqrt(parent_direction_unit[0]*parent_direction_unit[0] + parent_direction_unit[1]*parent_direction_unit[1] + parent_direction_unit[2]*parent_direction_unit[2]) << endl;
  //cout << "Parent to Unit: (" << displacement[0] << "," << displacement[1] << "," << displacement[2] << ")  to (" << parent_direction_unit[0] << "," << parent_direction_unit[1] << "," << parent_direction_unit[2] << ")" << endl;
  double daughter_direction_unit[3] = {a_daughter_vect->Px()/daughter_pTot,a_daughter_vect->Py()/daughter_pTot,a_daughter_vect->Pz()/daughter_pTot};
  //cout << "  unit check = " << sqrt(daughter_direction_unit[0]*daughter_direction_unit[0] + daughter_direction_unit[1]*daughter_direction_unit[1] + daughter_direction_unit[2]*daughter_direction_unit[2]) << endl;

  //cout << "Daughter to Unit: (" << a_daughter_vect->Px() << "," << a_daughter_vect->Py() << "," << a_daughter_vect->Pz() << ")  to (" << daughter_direction_unit[0] << "," << daughter_direction_unit[1] << "," << daughter_direction_unit[2] << ")" << endl;

  double outOfPlaneVect[3];
  //cout << "#####################initial cross product " << endl;
  crossProduct3D(daughter_direction_unit,parent_direction_unit,outOfPlaneVect[0],outOfPlaneVect[1],outOfPlaneVect[2]);

  //cout << "  Should be zero: " << daughter_direction_unit[0]*outOfPlaneVect[0] + daughter_direction_unit[1]*outOfPlaneVect[1]  + daughter_direction_unit[2]*outOfPlaneVect[2]  << endl;
  //cout << "  Should be zero: " << parent_direction_unit[0]*outOfPlaneVect[0] + parent_direction_unit[1]*outOfPlaneVect[1]  + parent_direction_unit[2]*outOfPlaneVect[2]  << endl;



  //double crossProductSize = sqrt(outOfPlaneVect[0]*outOfPlaneVect[0] + outOfPlaneVect[1]*outOfPlaneVect[1] + outOfPlaneVect[2]*outOfPlaneVect[2]);
  //double dca = crossProductSize*distance;
  //cout << "cross product size = " << crossProductSize << endl;
  //double crossProduct_unit[3] = {outOfPlaneVect[0]/crossProductSize, outOfPlaneVect[1]/crossProductSize, outOfPlaneVect[2]/crossProductSize};
  double dca[3];
  crossProduct3D(outOfPlaneVect,daughter_direction_unit,dca[0],dca[1],dca[2]);
  //cout << "  Should be zero: " << daughter_direction_unit[0]*dca[0] + daughter_direction_unit[1]*dca[1]  + daughter_direction_unit[2]*dca[2]  << endl;
  //cout << "  Should be zero: " << dca[0]*outOfPlaneVect[0] + dca[1]*outOfPlaneVect[1]  + dca[2]*outOfPlaneVect[2]  << endl;

  a_dca_x = dca[0] * distance;
  a_dca_y = dca[1] * distance;
  a_dca_z = dca[2] * distance;

}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
/*double PhysMath::dca_from_decay_helix(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect, int a_daughter_charge){
  if(a_daughter_charge == 0) return dca_from_decay(a_primaryVertex,a_secondaryVertex,a_daughter_vect);
  #ifndef _NO_HELIX_2_
  TVector3 momentum_inGeV(a_daughter_vect->Px(),a_daughter_vect->Py(),a_daughter_vect->Pz());
  TVector3 decay_vertex(0.01*a_primaryVertex->x(),0.01*a_primaryVertex->y(),0.01*a_primaryVertex->z());
  TVector3 event_vertex(0.01*a_secondaryVertex->x(),0.01*a_secondaryVertex->y(),0.01*a_secondaryVertex->z());
  Helix helix(momentum_inGeV,decay_vertex,-0.5*1E-11,a_daughter_charge);
  return 100.0*fabs(helix.geometricSignedDistance(event_vertex));
  #else
  return dca_from_decay(a_primaryVertex,a_secondaryVertex,a_daughter_vect);
  #endif
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


void PhysMath::dca_from_decay_helix(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect, int a_daughter_charge, double &a_dca_x, double &a_dca_y, double &a_dca_z){
  #ifndef _NO_HELIX_2_
  TVector3 momentum_inGeV(a_daughter_vect->Px(),a_daughter_vect->Py(),a_daughter_vect->Pz());
  TVector3 decay_vertex(0.01*a_primaryVertex->x(),0.01*a_primaryVertex->y(),0.01*a_primaryVertex->z());
  TVector3 event_vertex(0.01*a_secondaryVertex->x(),0.01*a_secondaryVertex->y(),0.01*a_secondaryVertex->z());
  Helix helix(momentum_inGeV,decay_vertex,-0.5*1E-11,a_daughter_charge);
  double dca = 0.0;
  helix.dca(event_vertex,dca,a_dca_x,a_dca_y,a_dca_z);
  dca = 100.0* fabs(dca);
  a_dca_x = 100.0* a_dca_x;
  a_dca_y = 100.0* a_dca_y;
  a_dca_z = 100.0* a_dca_z;
  #else
    dca_from_decay(a_primaryVertex,a_secondaryVertex,a_daughter_vect,a_dca_x,a_dca_y,a_dca_z);
  #endif
}*/









//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


void PhysMath::centerOfMassTwoBodyDecay(TRandom3* a_randObj, double a_parent_mass, double a_first_daughter_mass, double a_second_daughter_mass, TLorentzVector* &a_first_daughter_vect , TLorentzVector* &a_second_daughter_vect){
  double p_x = 0;
  double p_y = 0;
  double p_z = 0;
  a_randObj->Sphere(p_x, p_y, p_z, 1.0);
  double first_energy  = (a_parent_mass*a_parent_mass + a_first_daughter_mass*a_first_daughter_mass   - a_second_daughter_mass*a_second_daughter_mass)/(2.0*a_parent_mass);
  double second_energy = (a_parent_mass*a_parent_mass + a_second_daughter_mass*a_second_daughter_mass - a_first_daughter_mass*a_first_daughter_mass )/(2.0*a_parent_mass);
  double first_momentum  = sqrt(first_energy*first_energy - a_first_daughter_mass*a_first_daughter_mass);
  double second_momentum = sqrt(second_energy*second_energy - a_second_daughter_mass*a_second_daughter_mass);
  a_first_daughter_vect = new TLorentzVector();
  a_second_daughter_vect = new TLorentzVector();
  a_first_daughter_vect->SetPxPyPzE(   first_momentum*p_x,  first_momentum*p_y,    first_momentum*p_z, first_energy);
  a_second_daughter_vect->SetPxPyPzE(-second_momentum*p_x,-second_momentum*p_y,  -second_momentum*p_z, second_energy);
  return;
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void PhysMath::centerOfMassTwoBodyDecay(TRandom3* a_randObj, double a_parent_mass, double a_first_daughter_mass, double a_second_daughter_mass, ROOT::Math::PxPyPzEVector* &a_first_daughter_vect , ROOT::Math::PxPyPzEVector* &a_second_daughter_vect){
  double p_x = 0;
  double p_y = 0;
  double p_z = 0;
  a_randObj->Sphere(p_x, p_y, p_z, 1.0);
  double first_energy    = (a_parent_mass*a_parent_mass + a_first_daughter_mass*a_first_daughter_mass   - a_second_daughter_mass*a_second_daughter_mass)/(2.0*a_parent_mass);
  double second_energy   = (a_parent_mass*a_parent_mass + a_second_daughter_mass*a_second_daughter_mass - a_first_daughter_mass*a_first_daughter_mass )/(2.0*a_parent_mass);
  double first_momentum  = sqrt(first_energy*first_energy    -  a_first_daughter_mass*a_first_daughter_mass);
  double second_momentum = sqrt(second_energy*second_energy  -  a_second_daughter_mass*a_second_daughter_mass);
  a_first_daughter_vect  = new ROOT::Math::PxPyPzEVector();
  a_second_daughter_vect = new ROOT::Math::PxPyPzEVector();
  a_first_daughter_vect->SetPxPyPzE(   first_momentum*p_x,  first_momentum*p_y,    first_momentum*p_z, first_energy);
  a_second_daughter_vect->SetPxPyPzE(-second_momentum*p_x,-second_momentum*p_y,  -second_momentum*p_z, second_energy);
  #ifdef _PHYS_MATH_DEBUG_
    cout << "PhysMath::centerOfMassTwoBodyDecay: " << endl;
    cout << "  Daughter 1: " << "  E= " << first_energy << " p = " << first_momentum
        << " (" << a_first_daughter_vect->Px() << ", " << a_first_daughter_vect->Py() << ", " << a_first_daughter_vect->Pz() << ") " << endl;
    cout << "  Daughter 2: " << "  E= " << second_energy << " p = " << second_momentum
        << " (" << a_second_daughter_vect->Px() << ", " << a_second_daughter_vect->Py() << ", " << a_second_daughter_vect->Pz() << ") " << endl;
    cout << "  Mass Check: " << a_first_daughter_mass << "=?" << sqrt(first_energy*first_energy - first_momentum*first_momentum) << "=?" << sqrt(pow(a_first_daughter_vect->E(),2) - pow(a_first_daughter_vect->P(),2)) << endl;
    cout << "  Mass Check: " << a_second_daughter_mass << "=?" << sqrt(second_energy*second_energy - second_momentum*second_momentum) << "=?" << sqrt(pow(a_second_daughter_vect->E(),2) - pow(a_second_daughter_vect->P(),2)) << endl;
  #endif

  return;
}



//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void PhysMath::boostedThreeBodyDecay_Vect(double a_parent_mass,TLorentzVector* &a_parent_vect, double a_first_daughter_mass, double a_second_daughter_mass,double a_third_daughter_mass, TLorentzVector* &a_first_daughter_vect , TLorentzVector* &a_second_daughter_vect,TLorentzVector* &a_third_daughter_vect){
  TGenPhaseSpace generator;
  double masses[3];
  masses[0] = a_first_daughter_mass;
  masses[1] = a_second_daughter_mass;
  masses[2] = a_third_daughter_mass;
  generator.SetDecay(*a_parent_vect,3,masses);
  generator.Generate();

  a_first_daughter_vect = new TLorentzVector();
  a_second_daughter_vect = new TLorentzVector();
  a_third_daughter_vect = new TLorentzVector();
  a_first_daughter_vect->SetPxPyPzE(generator.GetDecay(0)->Px(),generator.GetDecay(0)->Py(),generator.GetDecay(0)->Pz(),generator.GetDecay(0)->E());
  a_second_daughter_vect->SetPxPyPzE(generator.GetDecay(1)->Px(),generator.GetDecay(1)->Py(),generator.GetDecay(1)->Pz(),generator.GetDecay(1)->E());
  a_third_daughter_vect->SetPxPyPzE(generator.GetDecay(2)->Px(),generator.GetDecay(2)->Py(),generator.GetDecay(2)->Pz(),generator.GetDecay(2)->E());

}


//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void PhysMath::boostedTwoBodyDecay_Vect(TRandom3* a_randObj, double a_parent_mass, TLorentzVector* &a_parent_vect, double a_first_daughter_mass, double a_second_daughter_mass,   TLorentzVector* &a_first_daughter_vect , TLorentzVector* &a_second_daughter_vect){
  #ifdef _PHYS_MATH_DEBUG_
    cout << " PhysMath::boostedTwoBodyDecay_Vect  with masses " << a_parent_mass << "  " << a_first_daughter_mass << "  " << a_second_daughter_mass << endl;
    cout << " Parent Rapidity: " << a_parent_vect->Rapidity() << " mTm0: " << a_parent_vect->Mt() - a_parent_mass << endl;
  #endif
  centerOfMassTwoBodyDecay(a_randObj, a_parent_mass, a_first_daughter_mass, a_second_daughter_mass, a_first_daughter_vect , a_second_daughter_vect);
  TVector3 boost = a_parent_vect->BoostVector();
  #ifdef _PHYS_MATH_DEBUG_
    cout << " Boost to Parent CM: " << boost.X() << " " << boost.Y() << "  " << boost.Z() << endl;
    cout << " Daugter 1 before boost: " << a_first_daughter_vect->Px() << " " << a_first_daughter_vect->Py() << " " << a_first_daughter_vect->Pz() << " " << a_first_daughter_vect->E() << "   y: " << a_first_daughter_vect->Rapidity() << "  mTm0: " << a_first_daughter_vect->Mt() - a_first_daughter_mass << endl;
    cout << " Daugter 2 before boost: " << a_second_daughter_vect->Px() << " " << a_second_daughter_vect->Py() << " " << a_second_daughter_vect->Pz() << " " << a_second_daughter_vect->E() << "   y: " << a_second_daughter_vect->Rapidity() << "  mTm0: " << a_second_daughter_vect->Mt() - a_second_daughter_mass  << endl;
  #endif
  a_first_daughter_vect->Boost(boost);
  a_second_daughter_vect->Boost(boost);
  #ifdef _PHYS_MATH_DEBUG_
    cout << " Daugter 1 after boost: " << a_first_daughter_vect->Px() << " " << a_first_daughter_vect->Py() << " " << a_first_daughter_vect->Pz() << " " << a_first_daughter_vect->E() << "   y: " << a_first_daughter_vect->Rapidity() << "  mTm0: " << a_first_daughter_vect->Mt() - a_first_daughter_mass  << endl;
    cout << " Daugter 2 after boost: " << a_second_daughter_vect->Px() << " " << a_second_daughter_vect->Py() << " " << a_second_daughter_vect->Pz() << " " << a_second_daughter_vect->E() << "   y: " << a_second_daughter_vect->Rapidity() << "  mTm0: " << a_second_daughter_vect->Mt() - a_second_daughter_mass  << endl;
  #endif
  return;
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void PhysMath::boostedTwoBodyDecay(ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                                double a_parent_mass, ROOT::Math::PxPyPzEVector* a_parent_vect,
                                double a_first_daughter_mass,  ROOT::Math::PxPyPzEVector* &a_first_daughter_vect,
                                double a_second_daughter_mass, ROOT::Math::PxPyPzEVector* &a_second_daughter_vect,
                                ROOT::Math::XYZVector* &a_SecondaryVertex, double a_parent_MeanLifetime){
  #ifdef _PHYS_MATH_DEBUG_
    cout << " PhysMath::boostedTwoBodyDecay with masses " << a_parent_mass << "  " << a_first_daughter_mass << "  " << a_second_daughter_mass << endl;
    cout << " Parent Rapidity: " << a_parent_vect->Rapidity() << " mTm0: " << a_parent_vect->Mt() - a_parent_mass << endl;
  #endif
  ROOT::Math::PxPyPzEVector* v_cm_d_1 = NULL;
  ROOT::Math::PxPyPzEVector* v_cm_d_2 = NULL;
  centerOfMassTwoBodyDecay(a_randObj, a_parent_mass, a_first_daughter_mass, a_second_daughter_mass, v_cm_d_1, v_cm_d_2);
  //a_boostObj->SetComponents(a_parent_vect->BoostToCM());

  double parent_momentum = a_parent_vect->P();
  double parent_beta     = a_parent_vect->Beta();
  a_boostObj->SetComponents(parent_beta*a_parent_vect->Px()/parent_momentum, parent_beta*a_parent_vect->Py()/parent_momentum, parent_beta*a_parent_vect->Pz()/parent_momentum);

  #ifdef _PHYS_MATH_DEBUG_
    cout << " Beta Check : " << parent_beta << "=?" << beta(parent_momentum, a_parent_mass) << endl;
  #endif

  ROOT::Math::PxPyPzEVector d1 = a_boostObj->operator()(*v_cm_d_1); //ridiculous, I know...
  ROOT::Math::PxPyPzEVector d2 = a_boostObj->operator()(*v_cm_d_2);
  a_first_daughter_vect  = new ROOT::Math::PxPyPzEVector(d1.Px(),d1.Py(),d1.Pz(),d1.E());
  a_second_daughter_vect = new ROOT::Math::PxPyPzEVector(d2.Px(),d2.Py(),d2.Pz(),d2.E());


  double decayLength = PhysMath::decayLength(a_randObj,parent_momentum, a_parent_mass, a_parent_MeanLifetime);
  if(a_parent_MeanLifetime == 0.0){
    decayLength = 0.0;
  }
  a_SecondaryVertex = new ROOT::Math::XYZVector(a_primaryVertex->x() + a_parent_vect->Px()*decayLength/parent_momentum,
                                                a_primaryVertex->y() + a_parent_vect->Py()*decayLength/parent_momentum,
                                                a_primaryVertex->z() + a_parent_vect->Pz()*decayLength/parent_momentum);

  #ifdef _PHYS_MATH_DEBUG_
    cout << "Boost Decay: " << endl;
    cout << "  Daughter 1: " << "  E= " << a_first_daughter_vect->E() << " p = " << a_first_daughter_vect->P()
        << " (" << a_first_daughter_vect->Px() << ", " << a_first_daughter_vect->Py() << ", " << a_first_daughter_vect->Pz() << ") " << endl;
    cout << "  Daughter 2: " << "  E= " << a_second_daughter_vect->E() << " p = " << a_second_daughter_vect->P()
        << " (" << a_second_daughter_vect->Px() << ", " << a_second_daughter_vect->Py() << ", " << a_second_daughter_vect->Pz() << ") " << endl;
    cout << "  Mass Check: " << a_first_daughter_mass << "=?" << sqrt(pow(a_first_daughter_vect->E(),2) - pow(a_first_daughter_vect->P(),2)) << endl;
    cout << "  Mass Check: " << a_second_daughter_mass << "=?" << sqrt(pow(a_second_daughter_vect->E(),2) - pow(a_second_daughter_vect->P(),2)) << endl;
    cout << "  Parent P=" << parent_momentum << "   decay length= " << decayLength << "     Second Vertex:   " << a_SecondaryVertex->X() << " " << a_SecondaryVertex->Y() << " " << a_SecondaryVertex->Z() << endl;
  #endif

  if(v_cm_d_1) delete v_cm_d_1;
  if(v_cm_d_2) delete v_cm_d_2;


  return;
}

//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void PhysMath::boostedThreeBodyDecay(TGenPhaseSpace* a_generator, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                                double a_parent_mass,          ROOT::Math::PxPyPzEVector* a_parent_vect,
                                double a_first_daughter_mass,  ROOT::Math::PxPyPzEVector* &a_first_daughter_vect,
                                double a_second_daughter_mass, ROOT::Math::PxPyPzEVector* &a_second_daughter_vect,
                                double a_third_daughter_mass,  ROOT::Math::PxPyPzEVector* &a_third_daughter_vect,
                                ROOT::Math::XYZVector* a_secondaryVertex, double a_parent_MeanLifetime){
  #ifdef _PHYS_MATH_DEBUG_
    cout << " PhysMath::boostedThreeBodyDecay" << endl;
  #endif

  double masses[3];
  masses[0] = a_first_daughter_mass;
  masses[1] = a_second_daughter_mass;
  masses[2] = a_third_daughter_mass;
  TLorentzVector parent_vect_oldClass;
  parent_vect_oldClass.SetPxPyPzE(a_parent_vect->Px(),a_parent_vect->Py(),a_parent_vect->Pz(),a_parent_vect->E());
  a_generator->SetDecay(parent_vect_oldClass,3,masses);
  a_generator->Generate();

  a_first_daughter_vect  = new ROOT::Math::PxPyPzEVector();
  a_second_daughter_vect = new ROOT::Math::PxPyPzEVector();
  a_third_daughter_vect  = new ROOT::Math::PxPyPzEVector();
  a_first_daughter_vect->SetPxPyPzE(a_generator->GetDecay(0)->Px(),a_generator->GetDecay(0)->Py(),a_generator->GetDecay(0)->Pz(),a_generator->GetDecay(0)->E());
  a_second_daughter_vect->SetPxPyPzE(a_generator->GetDecay(1)->Px(),a_generator->GetDecay(1)->Py(),a_generator->GetDecay(1)->Pz(),a_generator->GetDecay(1)->E());
  a_third_daughter_vect->SetPxPyPzE(a_generator->GetDecay(2)->Px(),a_generator->GetDecay(2)->Py(),a_generator->GetDecay(2)->Pz(),a_generator->GetDecay(2)->E());

  double parent_momentum = a_parent_vect->P();
  double decayLength = PhysMath::decayLength(a_randObj, parent_momentum, a_parent_mass, a_parent_MeanLifetime);
  if(a_parent_MeanLifetime == 0.0){
    decayLength = 0.0;
  }
  a_secondaryVertex = new ROOT::Math::XYZVector(a_primaryVertex->x() + a_parent_vect->Px()*decayLength/parent_momentum,
                                                a_primaryVertex->y() + a_parent_vect->Py()*decayLength/parent_momentum,
                                                a_primaryVertex->z() + a_parent_vect->Pz()*decayLength/parent_momentum);
}





//____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void PhysMath::boostedTwoBodyDecay(TRandom3* a_randObj, double a_parent_mass, double a_parent_rapidity,double a_parent_mTm0,
                                   double a_first_daughter_mass,  double &a_first_daughter_rapidity, double &a_first_daughter_mTm0,
                                   double a_second_daughter_mass, double &a_second_daughter_rapidity,double &a_second_daughter_mTm0){
  #ifdef _PHYS_MATH_DEBUG_
    cout << " PhysMath::boostedTwoBodyDecay" << endl;
    cout << " Parent Rapidity: " << a_parent_rapidity << " mTm0: " << a_parent_mTm0 << endl;
  #endif

  double parent_momentum = pTot_from_y_mTm0(a_parent_rapidity, a_parent_mTm0, a_parent_mass);
  double parent_energy   = sqrt(a_parent_mass*a_parent_mass + parent_momentum*parent_momentum);
  double parent_pT       = pT_from_mTm0(a_parent_mTm0, a_parent_mass);
  double parent_pz       = sqrt(parent_momentum*parent_momentum - parent_pT*parent_pT); // down beam line
  if(a_parent_rapidity < 0.0) parent_pz = -parent_pz;
  double parent_px = 0;
  double parent_py = 0;
  a_randObj->Circle(parent_px,parent_py,parent_pT); // initializes px and py of parrent

  #ifdef _PHYS_MATH_DEBUG_
    cout << "  E: " << parent_energy << " p: " << parent_momentum << " pT: " << parent_pT << " pZ: " << parent_pz << endl;
  #endif

  TLorentzVector* parent_vect = new TLorentzVector();
  parent_vect->SetPxPyPzE(parent_px,parent_py,parent_pz,parent_energy);
  TLorentzVector* first_daughter_vect = NULL;
  TLorentzVector* second_daughter_vect = NULL;
  boostedTwoBodyDecay_Vect(a_randObj,a_parent_mass, parent_vect, a_first_daughter_mass, a_second_daughter_mass,
                              first_daughter_vect,second_daughter_vect);
  a_first_daughter_rapidity  = first_daughter_vect->Rapidity();
  a_second_daughter_rapidity = second_daughter_vect->Rapidity();
  a_first_daughter_mTm0      = mTm0(a_first_daughter_mass, first_daughter_vect->Pt());
  a_second_daughter_mTm0     = mTm0(a_second_daughter_mass, second_daughter_vect->Pt());
  //a_first_daughter_pTot      = first_daughter_vect->P();
  //a_second_daughter_pTot     = second_daughter_vect->P();
  //double parent_direction_unit[3]          = {0,parent_pT/parent_momentum,parent_pz/parent_momentum};
  //double first_daughter_direction_unit[3]  = {first_daughter_vect->Px()/a_first_daughter_pTot,first_daughter_vect->Py()/a_first_daughter_pTot,first_daughter_vect->Pz()/a_first_daughter_pTot};
  //double second_daughter_direction_unit[3] = {second_daughter_vect->Px()/a_second_daughter_pTot,second_daughter_vect->Py()/a_second_daughter_pTot,second_daughter_vect->Pz()/a_second_daughter_pTot};
  //a_first_daughter_dca  = decayLength * crossProduct3D(first_daughter_direction_unit,parent_direction_unit);
  //a_second_daughter_dca = decayLength * crossProduct3D(second_daughter_direction_unit,parent_direction_unit);

  #ifdef _PHYS_MATH_DEBUG_
    cout << "  y y m m :" << a_first_daughter_rapidity << " " << a_second_daughter_rapidity << " " << a_first_daughter_mTm0 << " " << a_second_daughter_mTm0 << endl;
  #endif
  delete parent_vect;
  delete first_daughter_vect;
  delete second_daughter_vect;
  return;
}

void PhysMath::boostedThreeBodyDecay(TRandom3* a_randObj, double a_parent_mass, double a_parent_rapidity, double a_parent_mTm0, double a_first_daughter_mass, double &a_first_daughter_rapidity, double &a_first_daughter_mTm0, double a_second_daughter_mass, double &a_second_daughter_rapidity, double &a_second_daughter_mTm0, double a_third_daughter_mass, double &a_third_daughter_rapidity, double &a_third_daughter_mTm0){
  #ifdef _PHYS_MATH_DEBUG_
    cout << " PhysMath::boostedThreeBodyDecay" << endl;
    cout << " Parent Rapidity: " << a_parent_rapidity << " mTm0: " << a_parent_mTm0 << endl;
  #endif

  double parent_momentum = pTot_from_y_mTm0(a_parent_rapidity, a_parent_mTm0, a_parent_mass);
  double parent_energy   = sqrt(a_parent_mass*a_parent_mass + parent_momentum*parent_momentum);
  double parent_pT       = pT_from_mTm0(a_parent_mTm0, a_parent_mass);
  double parent_pz       = sqrt(parent_momentum*parent_momentum - parent_pT*parent_pT); // down beam line
  if(a_parent_rapidity < 0.0) parent_pz = -parent_pz;

  #ifdef _PHYS_MATH_DEBUG_
    cout << "  E: " << parent_energy << " p: " << parent_momentum << " pT: " << parent_pT << " pZ: " << parent_pz << endl;
  #endif

  TLorentzVector* parent_vect = new TLorentzVector();
  parent_vect->SetPxPyPzE(0.0,parent_pT,parent_pz,parent_energy);
  TLorentzVector* first_daughter_vect = NULL;
  TLorentzVector* second_daughter_vect = NULL;
  TLorentzVector* third_daughter_vect = NULL;
  boostedThreeBodyDecay_Vect(a_parent_mass, parent_vect, a_first_daughter_mass, a_second_daughter_mass,a_third_daughter_mass,first_daughter_vect,second_daughter_vect,third_daughter_vect);
  a_first_daughter_rapidity = first_daughter_vect->Rapidity();
  a_second_daughter_rapidity = second_daughter_vect->Rapidity();
  a_third_daughter_rapidity = third_daughter_vect->Rapidity();
  a_first_daughter_mTm0 = mTm0(a_first_daughter_mass, first_daughter_vect->Pt());
  a_second_daughter_mTm0  = mTm0(a_second_daughter_mass, second_daughter_vect->Pt());
  a_third_daughter_mTm0  = mTm0(a_second_daughter_mass, second_daughter_vect->Pt());
  #ifdef _PHYS_MATH_DEBUG_
    cout << "  y y y m m m:" << a_first_daughter_rapidity << " " << a_second_daughter_rapidity <<" " << a_third_daughter_rapidity <<  " " << a_first_daughter_mTm0 << " " << a_second_daughter_mTm0<< " " << a_third_daughter_mTm0 << endl;
  #endif
  delete parent_vect;
  delete first_daughter_vect;
  delete second_daughter_vect;
  delete third_daughter_vect;
  return;
}


double PhysMath::decayLength(TRandom3* a_randObj, double a_momentum, double a_mass, double a_meanLifetimeInSeconds){
  return 2.99792e10*betagamma(a_momentum, a_mass)*a_randObj->Exp(a_meanLifetimeInSeconds);
}



void PhysMath::lambdaDecay(TRandom3* a_randObj, double a_lambda_rapidity,double a_lambda_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pion_rapidity, double &a_pion_mTm0){
  boostedTwoBodyDecay(a_randObj, 1.115683,a_lambda_rapidity,a_lambda_mTm0, 0.938272081,a_proton_rapidity, a_proton_mTm0, 0.13957039, a_pion_rapidity, a_pion_mTm0);//,2.632E-10);
}

void PhysMath::lambdaDecay(    ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_lambda_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pion_vect, ROOT::Math::XYZVector* &a_SecondaryVertex){
  boostedTwoBodyDecay(a_boostObj,a_randObj,a_primaryVertex, 1.115683, a_lambda_vect,0.938272081,a_proton_vect, 0.13957039,
                      a_pion_vect,a_SecondaryVertex, 2.632E-10);
}

void PhysMath::kaonZeroShortDecay(TRandom3* a_randObj, double a_kaon_rapidity,double a_kaon_mTm0, double &a_piPlus_rapidity, double &a_piPlus_mTm0, double &a_piMinus_rapidity, double &a_piMinus_mTm0){
  boostedTwoBodyDecay(a_randObj, 0.497614,a_kaon_rapidity,a_kaon_mTm0, 0.13957039,a_piPlus_rapidity, a_piPlus_mTm0, 0.13957039, a_piMinus_rapidity, a_piMinus_mTm0);//,0.8954E-10);
}

void PhysMath::kaonZeroShortDecay(    ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_kaon_vect, ROOT::Math::PxPyPzEVector* &a_piPlus_vect, ROOT::Math::PxPyPzEVector* &a_piMinus_vect, ROOT::Math::XYZVector* &a_SecondaryVertex){
  boostedTwoBodyDecay(a_boostObj,a_randObj,a_primaryVertex, 0.497614, a_kaon_vect,0.13957039,a_piPlus_vect, 0.13957039,
                      a_piMinus_vect,a_SecondaryVertex, 0.8954E-10);
}

void PhysMath::sigmaPlusDecay(TRandom3* a_randObj, double a_sigma_rapidity,double a_sigma_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pion_rapidity, double &a_pion_mTm0){
  boostedTwoBodyDecay(a_randObj, 1.18937,a_sigma_rapidity,a_sigma_mTm0, 0.938272081,a_proton_rapidity, a_proton_mTm0, 0.1349768, a_pion_rapidity, a_pion_mTm0);//,0.8018E-10);
}

void PhysMath::sigmaPlusDecay(    ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_sigma_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pion_vect, ROOT::Math::XYZVector* &a_SecondaryVertex){
  boostedTwoBodyDecay(a_boostObj,a_randObj,a_primaryVertex, 1.18937, a_sigma_vect,0.938272081,a_proton_vect, 0.13957039,
                      a_pion_vect,a_SecondaryVertex, 0.8018E-10);
}

void PhysMath::xiMinusDecay(TRandom3* a_randObj, double a_xi_rapidity,double a_xi_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pion_rapidity, double &a_pion_mTm0, double &a_second_pion_rapidity, double &a_second_pion_mTm0){
  double lambda_rapidity = 0;
  double lambda_mTm0     = 0;
  double lambda_dca      = 0;
  boostedTwoBodyDecay(a_randObj, 1.32171,a_xi_rapidity,a_xi_mTm0, 1.115683,lambda_rapidity,lambda_mTm0, 0.13957039, a_pion_rapidity, a_pion_mTm0);//,1.639E-10);
  lambdaDecay(a_randObj, lambda_rapidity, lambda_mTm0, a_proton_rapidity, a_proton_mTm0, a_second_pion_rapidity, a_second_pion_mTm0);
}

void PhysMath::xiMinusDecay(    ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_xi_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pion_vect, ROOT::Math::PxPyPzEVector* &a_second_pion_vect, ROOT::Math::XYZVector* &a_SecondaryVertex, ROOT::Math::XYZVector* &a_ThirdVertex){
  ROOT::Math::PxPyPzEVector* lambda_vect = NULL;
  boostedTwoBodyDecay(a_boostObj,a_randObj,a_primaryVertex, 1.32171, a_xi_vect,1.115683,lambda_vect, 0.13957039,
                      a_pion_vect,a_SecondaryVertex, 1.639E-10);
  lambdaDecay(a_boostObj,a_randObj, a_SecondaryVertex, lambda_vect, a_proton_vect,a_second_pion_vect, a_ThirdVertex);
}


void PhysMath::xiZeroDecay(TRandom3* a_randObj, double a_xi_rapidity,double a_xi_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pionZero_rapidity, double &a_pionZero_mTm0, double &a_pionMinus_rapidity, double &a_pionMinus_mTm0){
  double lambda_rapidity = 0;
  double lambda_mTm0     = 0;
  double lambda_dca      = 0;
  boostedTwoBodyDecay(a_randObj, 1.31486,a_xi_rapidity,a_xi_mTm0, 1.115683,lambda_rapidity,lambda_mTm0, 0.1349768, a_pionZero_rapidity, a_pionZero_mTm0);//,2.90E-10,lambda_dca,a_pionZero_dca); //pi0
  lambdaDecay(a_randObj, lambda_rapidity, lambda_mTm0, a_proton_rapidity, a_proton_mTm0, a_pionMinus_rapidity, a_pionMinus_mTm0);
}

void PhysMath::xiZeroDecay( ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_xi_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pionZero_vect, ROOT::Math::PxPyPzEVector* &a_pionMinus_vect, ROOT::Math::XYZVector* &a_SecondaryVertex, ROOT::Math::XYZVector* &a_ThirdVertex){
  ROOT::Math::PxPyPzEVector* lambda_vect = NULL;
  boostedTwoBodyDecay(a_boostObj,a_randObj,a_primaryVertex, 1.31486, a_xi_vect,1.115683,lambda_vect, 0.1349768,
                      a_pionZero_vect,a_SecondaryVertex, 2.90E-10);
  lambdaDecay(a_boostObj,a_randObj, a_SecondaryVertex, lambda_vect, a_proton_vect,a_pionMinus_vect, a_ThirdVertex);
}


void PhysMath::omegaMinusDecay(TRandom3* a_randObj, double a_omega_rapidity,double a_omega_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_kaon_rapidity, double &a_kaon_mTm0, double &a_pionMinus_rapidity, double &a_pionMinus_mTm0, double &a_second_pion_rapidity, double &a_second_pion_mTm0){

  double branch = a_randObj->Uniform();
  if(branch < 0.678){ // LAMBDA DECAY
    double lambda_rapidity = 0;
    double lambda_mTm0 = 0;
    boostedTwoBodyDecay(a_randObj, 1.67245,a_omega_rapidity,a_omega_mTm0, 1.115683,lambda_rapidity,lambda_mTm0, 0.493677, a_kaon_rapidity, a_kaon_mTm0);
    if(a_randObj->Uniform() > 0.639) return;
    lambdaDecay(a_randObj, lambda_rapidity, lambda_mTm0, a_proton_rapidity, a_proton_mTm0, a_pionMinus_rapidity, a_pionMinus_mTm0);
  }else if(branch < 0.914){ // XI ZERO DECAY

    double xiZero_rapidity = 0;
    double xiZero_mTm0 = 0;
    double pionZero_rapidity = 0.0;
    double pionZero_mTm0 = 0.0;
    boostedTwoBodyDecay(a_randObj, 1.67245,a_omega_rapidity,a_omega_mTm0, 1.31486,xiZero_rapidity,xiZero_mTm0, 0.13957039, a_pionMinus_rapidity, a_pionMinus_mTm0);
    if(a_randObj->Uniform() > 0.99525) return;
    xiZeroDecay(a_randObj, xiZero_rapidity, xiZero_mTm0, a_proton_rapidity, a_proton_mTm0, pionZero_rapidity, pionZero_mTm0, a_second_pion_rapidity, a_second_pion_mTm0);
  }else{ // XI MINUS DECAY
    if(a_randObj->Uniform() > 0.99525) return;
    double xi_rapidity = 0;
    double xi_mTm0 = 0;
    double pionZero_rapidity = 0.0;
    double pionZero_mTm0 = 0.0;
    boostedTwoBodyDecay(a_randObj, 1.67245,a_omega_rapidity,a_omega_mTm0, 1.32171,xi_rapidity,xi_mTm0, 0.1349768, pionZero_rapidity, pionZero_mTm0);
    if(a_randObj->Uniform() > 0.99887) return;
    xiMinusDecay(a_randObj, xi_rapidity, xi_mTm0, a_proton_rapidity, a_proton_mTm0, a_pionMinus_rapidity, a_pionMinus_mTm0, a_second_pion_rapidity, a_second_pion_mTm0);
  }


}



void PhysMath::omegaMinusDecay( ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
          ROOT::Math::PxPyPzEVector* a_omega_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pionZero_vect,ROOT::Math::PxPyPzEVector* &a_kaon_vect, ROOT::Math::PxPyPzEVector* &a_pionMinus_vect, ROOT::Math::PxPyPzEVector* &a_second_pion_vect, ROOT::Math::XYZVector* &a_SecondaryVertex, ROOT::Math::XYZVector* &a_ThirdVertex, ROOT::Math::XYZVector* &a_FourthVertex,int &a_proton_vertex_index, int &a_kaon_vertex_index, int &a_pion_vertex_index, int &a_second_pion_vertex_index){

  a_proton_vertex_index = -1;
  a_kaon_vertex_index = -1;
  a_pion_vertex_index = -1;
  a_second_pion_vertex_index = -1;

  double branch = a_randObj->Uniform();
  if(branch < 0.678){ // LAMBDA DECAY
    ROOT::Math::PxPyPzEVector* lambda_vect = NULL;
    boostedTwoBodyDecay(a_boostObj,a_randObj,a_primaryVertex, 1.67245, a_omega_vect,1.115683,lambda_vect, 0.493677,
                      a_kaon_vect,a_SecondaryVertex, 0.821E-10);
    a_kaon_vertex_index = 2;
    if(a_randObj->Uniform() > 0.639) return;
    lambdaDecay(a_boostObj,a_randObj, a_SecondaryVertex, lambda_vect, a_proton_vect,a_pionMinus_vect, a_ThirdVertex);
    a_proton_vertex_index = 3;
    a_pion_vertex_index = 3;
  }else if(branch < 0.914){ // XI ZERO DECAY

    ROOT::Math::PxPyPzEVector* xiZero_vect   = NULL;
    ROOT::Math::PxPyPzEVector* pionZero_vect = NULL;
    boostedTwoBodyDecay(a_boostObj,a_randObj,a_primaryVertex, 1.67245, a_omega_vect,1.31486, xiZero_vect, 0.13957039,
                      a_pionMinus_vect,a_SecondaryVertex, 0.821E-10);
    a_pion_vertex_index = 2;
    if(a_randObj->Uniform() > 0.99525) return;
    xiZeroDecay(a_boostObj,a_randObj, a_SecondaryVertex, xiZero_vect, a_proton_vect, pionZero_vect,a_second_pion_vect, a_ThirdVertex, a_FourthVertex);
    a_proton_vertex_index = 4;
    a_second_pion_vertex_index = 4;
  }else{ // XI MINUS DECAY
    ROOT::Math::PxPyPzEVector* xi_vect   = NULL;
    ROOT::Math::PxPyPzEVector* pionZero_vect = NULL;
    boostedTwoBodyDecay(a_boostObj,a_randObj,a_primaryVertex, 1.67245, a_omega_vect,1.32171,xi_vect, 0.1349768,
            pionZero_vect,a_SecondaryVertex, 0.821E-10);
    if(a_randObj->Uniform() > 0.99887) return;
    xiMinusDecay(a_boostObj,a_randObj, a_SecondaryVertex, xi_vect, a_proton_vect, a_pionMinus_vect, a_second_pion_vect, a_ThirdVertex, a_FourthVertex);
    a_pion_vertex_index   = 3;
    a_proton_vertex_index = 4;
    a_second_pion_vertex_index = 4;
  }

}








void PhysMath::sigmaZeroDecay(TRandom3* a_randObj, double a_sigma_rapidity,double a_sigma_mTm0, double &a_proton_rapidity, double &a_proton_mTm0, double &a_pion_rapidity, double &a_pion_mTm0){
  double lambda_rapidity = 0;
  double lambda_mTm0 = 0;
  double gammaRapidity = 0;
  double gammamTm0 = 0;
  boostedTwoBodyDecay(a_randObj, 1.192642,a_sigma_rapidity,a_sigma_mTm0, 1.115683,lambda_rapidity,lambda_mTm0, 0.0, gammaRapidity, gammamTm0);
  lambdaDecay(a_randObj, lambda_rapidity, lambda_mTm0, a_proton_rapidity, a_proton_mTm0, a_pion_rapidity, a_pion_mTm0);
}


void PhysMath::sigmaZeroDecay( ROOT::Math::Boost* a_boostObj, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_sigma_vect, ROOT::Math::PxPyPzEVector* &a_proton_vect, ROOT::Math::PxPyPzEVector* &a_pion_vect, ROOT::Math::XYZVector* &a_SecondaryVertex, ROOT::Math::XYZVector* &a_ThirdVertex){
  ROOT::Math::PxPyPzEVector* lambda_vect = NULL;
  ROOT::Math::PxPyPzEVector* gamma_vect  = NULL;
  boostedTwoBodyDecay(a_boostObj,a_randObj,a_primaryVertex,1.192642, a_sigma_vect,1.115683,lambda_vect, 0.0,
                      gamma_vect,a_SecondaryVertex, 7.4E-20);
  lambdaDecay(a_boostObj,a_randObj, a_SecondaryVertex, lambda_vect, a_proton_vect,a_pion_vect, a_ThirdVertex);
}




void PhysMath::etaDecay(TRandom3* a_randObj, double a_eta_rapidity,double a_eta_mTm0, double &a_piPlus_rapidity, double &a_piPlus_mTm0, double &a_piMinus_rapidity, double &a_piMinus_mTm0){
  double branch = a_randObj->Uniform();
  if(branch < 0.2292){
    double piZeroRap = 0;
    double piZeroMtm0 = 0;
    boostedThreeBodyDecay(a_randObj, 0.547862, a_eta_rapidity, a_eta_mTm0, 0.13957039, a_piPlus_rapidity, a_piPlus_mTm0, 0.13957039, a_piMinus_rapidity, a_piMinus_mTm0, 0.1349768, piZeroRap, piZeroMtm0);
  }else if(branch < 0.2714){ //0.2292+0.0422
    double gamma_rapidity = 0.0;
    double gamma_mTm0 = 0.0;
    boostedThreeBodyDecay(a_randObj, 0.547862, a_eta_rapidity, a_eta_mTm0, 0.13957039, a_piPlus_rapidity, a_piPlus_mTm0, 0.13957039, a_piMinus_rapidity, a_piMinus_mTm0, 0.0, gamma_rapidity, gamma_mTm0);
  }
}


void PhysMath::etaDecay( TGenPhaseSpace* a_generator, TRandom3* a_randObj, ROOT::Math::XYZVector* a_primaryVertex,
                       ROOT::Math::PxPyPzEVector* a_eta_vect, ROOT::Math::PxPyPzEVector* &a_piPlus_vect, ROOT::Math::PxPyPzEVector* &a_piMinus_vect, ROOT::Math::XYZVector* &a_SecondaryVertex){
  double branch = a_randObj->Uniform();
  if(branch < 0.2292){
    ROOT::Math::PxPyPzEVector* piZero_vect = NULL;
    boostedThreeBodyDecay(a_generator,a_randObj,a_primaryVertex,0.547862, a_eta_vect,0.13957039,a_piPlus_vect, 0.13957039, a_piMinus_vect,0.1349768, piZero_vect,a_SecondaryVertex, 5.0E-19);
  }else if(branch < 0.2714){ //0.2292+0.0422
    ROOT::Math::PxPyPzEVector* gamma_vect = NULL;
    boostedThreeBodyDecay(a_generator,a_randObj,a_primaryVertex,0.547862, a_eta_vect,0.13957039,a_piPlus_vect, 0.13957039, a_piMinus_vect,0.0, gamma_vect,a_SecondaryVertex, 5.0E-19);
  }
}







TF1* PhysMath::getSimpleExponentialFunction(string a_name, double a_scale, double a_exponent){
  TF1* funct = new TF1(a_name.c_str(),
    "[0]*exp([1]*x)",
    -5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlue);
  funct->SetParNames("A","b");
  funct->SetParameter(0, a_scale);
  //funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, a_exponent);
  //funct->SetParLimits(1,0.0,1000.0);
  return funct;
}



TF1* PhysMath::getExponentialFunction(string a_name, double a_scale, double a_exponent, double a_leftRightShift, double a_y0){
  TF1* funct = new TF1(a_name.c_str(),
    "[0]*exp([1]*(x - [2])) + [3]",
    -5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlue);
  funct->SetParNames("A","b","x_{0}","y_{0}");
  funct->SetParameter(0, a_scale);
  //funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, a_exponent);
  //funct->SetParLimits(1,0.0,1000.0);
  funct->SetParameter(2, a_leftRightShift);
  //funct->SetParLimits(2, -2.0, 2.0);
  funct->SetParameter(3, a_y0);
  //funct->SetParLimits(3, -2.0, 2.0);
  return funct;
}

TF1* PhysMath::getFitSafeExponentialFunction(string a_name, double a_scale, double a_exponent){
  TF1* funct = new TF1(a_name.c_str(),"[0]*exp(-x/[1])/[1]",-5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlue);
  funct->SetParNames("N","T");
  funct->SetParameter(0, a_scale);
  funct->SetParLimits(0, 0.01*a_scale, 100.0*a_scale);
  funct->SetParameter(1, a_exponent);
  funct->SetParLimits(1,0.001*a_exponent,1000.0*a_exponent);
  return funct;
}

TF1* PhysMath::getFitSafeDoubleExponentialFunction(string a_name, double a_scale, double a_exponent, double a_scale_2, double a_exponent_2){
  TF1* funct = new TF1(a_name.c_str(),"[0]*exp(-x/[1])/[1] + [2]*exp(-x/[3])/[3]",-10.0,10.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlue);
  funct->SetParNames("N_1","T_1","N_2","T_2");
  funct->SetParameter(0, a_scale);
  funct->SetParLimits(0, 0.01*a_scale, 100.0*a_scale);
  funct->SetParameter(1, a_exponent);
  funct->SetParLimits(1,0.001*a_exponent,1000.0*a_exponent);
  funct->SetParameter(2, a_scale_2);
  funct->SetParLimits(2, 0.01*a_scale_2, 100.0*a_scale_2);
  funct->SetParameter(3, a_exponent_2);
  funct->SetParLimits(3,0.001*a_exponent_2,1000.0*a_exponent_2);
  return funct;
}
TF1* PhysMath::getFitSafeTripleExponentialFunction(string a_name, double a_scale, double a_exponent, double a_scale_2, double a_exponent_2,double a_scale_3, double a_exponent_3){
  TF1* funct = new TF1(a_name.c_str(),"[0]*exp(-x/[1])/[1] + [2]*exp(-x/[3])/[3] + [4]*exp(-x/[5])/[5]",-10.0,10.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlue);
  funct->SetParNames("N_1","T_1","N_2","T_2","N_3","T_3");
  funct->SetParameter(0, a_scale);
  funct->SetParLimits(0, 0.01*a_scale, 100.0*a_scale);
  funct->SetParameter(1, a_exponent);
  funct->SetParLimits(1,0.001*a_exponent,1000.0*a_exponent);
  funct->SetParameter(2, a_scale_2);
  funct->SetParLimits(2, 0.01*a_scale_2, 100.0*a_scale_2);
  funct->SetParameter(3, a_exponent_2);
  funct->SetParLimits(3,0.001*a_exponent_2,1000.0*a_exponent_2);
  funct->SetParameter(4, a_scale_3);
  funct->SetParLimits(4, 0.01*a_scale_3, 100.0*a_scale_3);
  funct->SetParameter(5, a_exponent_3);
  funct->SetParLimits(5,0.001*a_exponent_3,1000.0*a_exponent_3);
  return funct;
}

TF1* PhysMath::getExpoPoly4Function(string a_name,double a_scale, double a_p0, double a_p1, double a_p2, double a_p3, double a_p4){
  TF1* funct = new TF1(a_name.c_str(),"[0]*exp([1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x)",-10.0,10.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlue);
  funct->SetParNames("N","P_0","P_1","P_2","P_3","P_4");
  funct->SetParameter(0, a_scale);
  //funct->SetParLimits(0, 0.01*a_scale, 100.0*a_scale);
  funct->SetParameter(1, a_p0);
  //funct->SetParLimits(1,0.001*a_exponent,1000.0*a_exponent);
  funct->SetParameter(2, a_p1);
  //funct->SetParLimits(2, 0.01*a_scale_2, 100.0*a_scale_2);
  funct->SetParameter(3, a_p2);
  //funct->SetParLimits(3,0.001*a_exponent_2,1000.0*a_exponent_2);
  funct->SetParameter(4, a_p3);
  //funct->SetParLimits(4, 0.01*a_scale_3, 100.0*a_scale_3);
  funct->SetParameter(5, a_p4);
  //funct->SetParLimits(5,0.001*a_exponent_3,1000.0*a_exponent_3);
  return funct;
}



TF1* PhysMath::getStudentTFunction(string a_name, double a_N, double a_mu, double a_sig, double a_nu ){
  TF1* funct = new TF1(a_name.c_str(),
    "[0]*TMath::Student((x - [1])/[2], [3])/[2]",
    -10,10);
  funct->SetParNames("N","#mu","#sigma","#nu");
  funct->SetParameter(0,a_N);
  funct->SetParameter(1,a_mu);
  funct->SetParameter(2,a_sig);
  funct->SetParameter(3,a_nu);
  return funct;
}

TF1* PhysMath::getM2StudentTFunction(string a_name, double a_N, double a_mu, double a_sig, double a_nu, double a_momentum){
  TF1* funct = new TF1(a_name.c_str(),
    "[0]*TMath::Student((sqrt(x/([4]*[4]) + 1.0) - [1])/[2], [3])/(2.0*[2]*[4]*[4]*sqrt(x/([4]*[4]) + 1.0))",
    -10,10);
  funct->SetParNames("N","#mu","#sigma","#nu","p_Tot");
  funct->SetParameter(0,1.0);
  funct->SetParameter(1,a_mu);
  funct->SetParameter(2,a_sig);
  funct->SetParameter(3,a_nu);
  funct->FixParameter(4,a_momentum);
  //funct->SetParameter(5,1.0);
  //double lowRange  = PhysMath::massSquared_signed_fromInvBeta(a_momentum,a_mu - 6.0*a_sig);
  //double highRange = PhysMath::massSquared_signed_fromInvBeta(a_momentum,a_mu + 6.0*a_sig);
  //double normalization = funct->Integral(lowRange,highRange)/(highRange - lowRange);
  //funct->SetParameter(5,1.0/normalization);
  funct->SetParameter(0,a_N);
  return funct;
}



TF1* PhysMath::getGaussianFunction(string a_name, double a_N, double a_mu, double a_sig){
  TF1* funct = new TF1(a_name.c_str(),
    "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",
    -200,200);
  funct->SetParNames("N","#mu","#sigma");
  funct->SetParameter(0,a_N);
  funct->SetParameter(1,a_mu);
  funct->SetParameter(2,a_sig);
  return funct;
}

TF1* PhysMath::getTwoGaussianFunction(string a_name,double a_N1, double a_mu1, double a_sig1,double a_N2, double a_mu2, double a_sig2){
  TF1* funct = new TF1(a_name.c_str(),
    "0.3989422804 * ([0]*exp(-0.5*( x - [1] )*(x - [1])/([2]*[2]))/[2] + [3]*exp(-0.5*( x - [4] )*(x - [4])/([5]*[5]))/[5])",
    -200,200);
  funct->SetParNames("N_{1}","#mu_{1}","#sigma_{1}","N_{2}","#mu_{2}","#sigma_{2}");
  funct->SetParameter(0,a_N1);
  funct->SetParameter(1,a_mu1);
  funct->SetParameter(2,a_sig1);
  funct->SetParameter(3,a_N2);
  funct->SetParameter(4,a_mu2);
  funct->SetParameter(5,a_sig2);
  funct->SetNpx(300);
  return funct;
}

TF1* PhysMath::getGaussianConvolutionWithStep(string a_name, double a_N, double a_mu, double a_sig, double a_width){
  TF1* funct = new TF1(a_name.c_str(),"[0]*0.5*(erf(([1] - x + 0.5*[3])/(sqrt(2.0)*[2])) - erf(([1] - x - 0.5*[3])/(sqrt(2.0)*[2])))/[3]",a_mu - 0.5*a_width - 15.0*a_sig,a_mu + 0.5*a_width + 15.0*a_sig);
  funct->SetParNames("N","#mu","#sigma","width");
  funct->SetParameter(0,a_N);
  funct->SetParameter(1,a_mu);
  funct->SetParameter(2,a_sig);
  funct->SetParameter(3,a_width);
  funct->SetNpx(300);
  return funct;
}





TF1* PhysMath::getGaussianFunction_byDerivatives_GausPars(string a_name,double a_N, double a_mu, double a_sig, double a_x0){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::gaus_by_deriviatives(x,[0],[1],[2],[3])",-200,200);
  funct->SetParNames("y_0","y_slope","y_secDer","x_0");
  double par[3] = {0,0,0};
  PhysMath::derivatives_from_gausParams(par[0],par[1],par[2],a_x0,a_N, a_mu,a_sig);
  funct->SetParameter(0,par[0]);
  funct->SetParameter(1,par[1]);
  funct->SetParameter(2,par[2]);
  funct->SetParameter(3,a_x0);
  return funct;
}
TF1* PhysMath::getGaussianFunction_byDerivatives(string a_name,double a_yLoc, double a_ySlope, double a_ySecDer, double a_xLoc){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::gaus_by_deriviatives(x,[0],[1],[2],[3])",-200,200);
  funct->SetParNames("y_0","y_slope","y_secDer","x_0");
  funct->SetParameter(0,a_yLoc);
  funct->SetParameter(1,a_ySlope);
  funct->SetParameter(2,a_ySecDer);
  funct->SetParameter(3,a_xLoc);
  return funct;
}

TF1* PhysMath::getDoubleGaussianFunction(string a_name,double a_N, double a_mu, double a_sig1, double a_fracN, double a_sig2){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::doubleGauss([0],[1],[2],[3],[4],x)",
    //"[0]*0.3989422804*((1.0-[3])*exp(-0.5*pow((x-[1])/[2],2))/[2] + [3]*exp(-0.5*pow((x-[1])/[2]))/[4])",
    -10,10);
  funct->SetParNames("N","#mu","#sigma_{1}","r","#sigma_{2}");
  funct->SetParameter(0,a_N);
  funct->SetParameter(1,a_mu);
  funct->SetParameter(2,a_sig1);
  funct->SetParameter(3,a_fracN);
  funct->SetParameter(4,a_sig2);
  funct->SetParLimits(0,0.001*a_N,1000.0*a_N);
  funct->SetParLimits(1,a_mu - 1.5*a_sig1,a_mu + 1.5*a_sig1);
  funct->SetParLimits(2,0.01*a_sig1,10.0*a_sig1);
  funct->SetParLimits(3,0.0,1.0);
  funct->SetParLimits(4,0.02*a_sig2,20.0*a_sig2);
  return funct;
}
TF1* PhysMath::getGaussAndStudentTFunction(string a_name,double a_N, double a_mu, double a_sig1, double a_nu, double a_fracN, double a_sig2){
  TF1* funct = new TF1(a_name.c_str(),
    "[0]*((1.0 - [4])*TMath::Student((x - [1])/[2], [3])/[2] + [4]*0.3989422804*exp(-0.5*(x-[1])*(x-[1])/([5]*[5]))/[5])",
    -10,10);
  funct->SetParNames("N","#mu","#sigma_{1}","#nu","r","#sigma_{2}");
  funct->SetParameter(0,a_N);
  funct->SetParameter(1,a_mu);
  funct->SetParameter(2,a_sig1);
  funct->SetParameter(3,a_nu);
  funct->SetParameter(4,a_fracN);
  funct->SetParameter(5,a_sig2);
  funct->SetParLimits(0,0.001*a_N,1000.0*a_N);
  funct->SetParLimits(1,a_mu - 1.5*a_sig1,a_mu + 1.5*a_sig1);
  funct->SetParLimits(2,0.01*a_sig1,10.0*a_sig1);
  funct->SetParLimits(3,1.0,10000.0);
  funct->SetParLimits(4,0.0,1.0);
  funct->SetParLimits(5,0.02*a_sig2,20.0*a_sig2);

  return funct;
}


string PhysMath::getGaussianFunctionString(int a_startParam){
  string text = Form("[%d]*0.3989422804*TMath::Exp(-0.5*(x-[%d])*(x-[%d])/([%d]*[%d]))/[%d]",
    a_startParam,a_startParam+1,a_startParam+1,a_startParam+2,a_startParam+2,a_startParam+2);
  return text;
};

string PhysMath::getDoubleGaussianString(int a_startParam){
  string text = Form("PhysMath::doubleGauss([%d],[%d],[%d],[%d],[%d],x)",
    a_startParam,a_startParam+1,a_startParam+2,a_startParam+3,a_startParam+4);
  //string text = Form("[%d]*0.3989422804*((1.0-[%d])*exp(-0.5*(x-[%d])*(x-[%d])/([%d]*[%d]))/[%d] + [%d]*exp(-0.5*(x-[%d])*(x-[%d])/([%d]*[%d]))/[%d])",
  //  a_startParam,a_startParam+3,a_startParam+1,a_startParam+1,a_startParam+2,a_startParam+2,a_startParam+2,
  //  a_startParam+3,a_startParam+1,a_startParam+1,a_startParam+4,a_startParam+4,a_startParam+4);
  return text;
};

string PhysMath::getGaussAndStudentTFunctionString(int a_startParam){
  string text = Form("[%d]*((1.0-[%d])*TMath::Student((x-[%d])/[%d], [%d])/[%d] + [%d]*0.3989422804*exp(-0.5*(x-[%d])*(x-[%d])/([%d]*[%d]))/[%d])",
    a_startParam,a_startParam+4,a_startParam+1,a_startParam+2,a_startParam+3,a_startParam+2,
    a_startParam+4,a_startParam+1,a_startParam+1,a_startParam+5,a_startParam+5,a_startParam+5);
  return text;
};

double PhysMath::doubleGauss(double a_N, double a_mu, double a_sig1, double a_fracN, double a_sig2,double a_x){
  double loc1 = (a_x-a_mu)/a_sig1;
  double loc2 = (a_x-a_mu)/a_sig2;
  double value = a_N*0.3989422804*((1.0-a_fracN)*exp(-0.5*loc1*loc1)/a_sig1 + a_fracN*exp(-0.5*loc2*loc2)/a_sig2);
  if(value < 0.0) return 0.0;
  return value;
}

double PhysMath::gaussAndStudentT(double a_N, double a_mu, double a_sig1, double a_nu, double a_fracN, double a_sig2, double a_x){
  double loc = a_x - a_mu;
  return a_N*((1.0 - a_fracN)*TMath::Student(loc/a_sig1,a_nu)/a_sig1 + a_fracN*0.3989422804*exp(-0.5*loc*loc/(a_sig2*a_sig2))/a_sig2);

}




TF1* PhysMath::getFourierTrigFunction(string a_name, int a_order, double a_xScale){
  //Parameters: const sinAmp cosAmp sinAmp cosAmp... xScale
  if(a_order < 0){
    cout << "ERROR: PhysMath::getFourierTrigFunction: you requested an order < 0. I don't know what to do."<< endl;
    return new TF1(a_name.c_str(),"0.0",-1.0,1.0);
  }

  string functString = "[0]";
  int indexForXScale = 2*a_order + 1;
  for(int iii = 1; iii <= a_order; iii++){
    double freq = (double) iii;
    functString += Form("+[%d]*sin(%1.1f*[%d]*x)+[%d]*cos(%1.1f*[%d]*x)",2*iii - 1,freq,indexForXScale,2*iii,freq,indexForXScale);
  }

  TF1* funct = new TF1(a_name.c_str(),functString.c_str(),-1,1);
  funct->SetParName(0,"A_0");
  funct->SetParameter(0,0.0);
  funct->SetParName(indexForXScale,"x_scale");
  funct->SetParameter(indexForXScale,1.0);
  for(int iii = 1; iii <= a_order; iii++){
    funct->SetParName(2*iii-1,Form("A_%d",iii));
    funct->SetParName(2*iii,Form("B_%d",iii));
    funct->SetParameter(2*iii-1,1.0);
    funct->SetParameter(2*iii,1.0);
  }
  return funct;

}

TF1* PhysMath::getLegendreExpansion(string a_name, int a_order, bool a_fixRange, double a_lowRange, double a_highRange){
  //Parameters: A_0 A_1 ... xlowRange xhighRange // as you would expect
  if(a_order < 0){
    cout << "ERROR: PhysMath::getLegendreExpansion: you requested an order < 0. I don't know what to do."<< endl;
    return new TF1(a_name.c_str(),"0.0",-1.0,1.0);
  }
  if(a_order > 29){
    cout << "ERROR: PhysMath::getLegendreExpansion: you requested an order > 29. I don't know what to do. I only implement up to order 29. You can add them if you need."<< endl;
    return new TF1(a_name.c_str(),"0.0",-1.0,1.0);
  }

  int indexForLowRange = a_order;
  int indexForHighRange = a_order+1;
  string functString = "[0]";
  for(int iii = 1; iii < a_order; iii++){
    functString += Form("+[%d]*PhysMath::legendre(%d,(2.0*x-[%d]-[%d])/fabs([%d]-[%d]))", iii,iii,indexForLowRange,indexForHighRange,indexForHighRange,indexForLowRange);
  }

  TF1* funct = new TF1(a_name.c_str(),functString.c_str(),a_lowRange,a_highRange);
  for(int iii = 0; iii <= a_order; iii++){
    funct->SetParName(iii,Form("A_%d",iii));
    funct->SetParameter(iii,1.0);
  }
  funct->SetParName(indexForLowRange,"x_lowRange");
  funct->SetParName(indexForHighRange,"x_highRange");

  if(a_fixRange){
    funct->FixParameter(indexForLowRange,a_lowRange);
    funct->FixParameter(indexForHighRange,a_highRange);
  }else{
    funct->SetParameter(indexForLowRange, a_lowRange);
    funct->SetParameter(indexForHighRange,a_highRange);
  }
  funct->SetRange(a_lowRange, a_highRange);

  return funct;
}



TF1* PhysMath::getTsallisFunction(string a_name,double a_dNdy,double a_T,double a_q_minus_one,double a_m){
  //  0 a_dNdy,
  //  1 a_T,
  //  2 a_q_minus_one,
  //  3 a_m
  TF1* funct = new TF1(a_name.c_str(),"[0]*(1.0-[2])*pow(1.0 + x*[2]/([1]+[3]*[2]),-([2]+1.0)/[2])/(([1]+[3])*([1] + [3]*[2]))/6.2831853072",0.0,8.0);

  funct->SetParName(0,"dNdy");
  funct->SetParName(1,"T");
  funct->SetParName(2,"q-1"); // temp
  funct->SetParName(3,"m");

  funct->SetParameter(0,a_dNdy);
  funct->SetParameter(1,a_T);
  funct->SetParameter(2,a_q_minus_one); // temp
  funct->FixParameter(3,a_m);

  funct->SetParLimits(0,0.001*a_dNdy,1000.0*a_dNdy);
  funct->SetParLimits(1,0.0,1000.0*a_T); // temp
  funct->SetParLimits(2,0.0,1.0); //this number is going to be close to 0.0 bc that is the exponential function
  funct->SetNpx(300);
  return funct;
}

TF1* PhysMath::getTsallisFunction_NotInvariant(string a_name,double a_dNdy,double a_T,double a_q_minus_one,double a_m){
  //  0 a_dNdy,
  //  1 a_T,
  //  2 a_q_minus_one,
  //  3 a_m
  TF1* funct = new TF1(a_name.c_str(),"(x+[3])*[0]*(1.0-[2])*pow(1.0 + x*[2]/([1]+[3]*[2]),-([2]+1.0)/[2])/(([1]+[3])*([1]+[3]*[2]))",0.0,8.0);
  funct->SetParName(0,"dNdy");
  funct->SetParName(1,"T");
  funct->SetParName(2,"q-1"); // temp
  funct->SetParName(3,"m");

  funct->SetParameter(0,a_dNdy);
  funct->SetParameter(1,a_T);
  funct->SetParameter(2,a_q_minus_one); // temp
  funct->FixParameter(3,a_m);

  funct->SetParLimits(0,0.001*a_dNdy,1000.0*a_dNdy);
  funct->SetParLimits(1,0.0,1000.0*a_T); // temp
  funct->SetParLimits(2,0.0,1.0); //this number is going to be close to 0.0 bc that is the exponential function
  funct->SetNpx(300);
  return funct;
}





string PhysMath::getLegendreExpansionString(int a_order){
  int indexForLowRange = a_order;
  int indexForHighRange = a_order+1;
  //string functString = "[0]";
  string functString = Form("[0]*PhysMath::legendre(0,(2.0*x-[%d]-[%d])/fabs([%d]-[%d]))",indexForLowRange,indexForHighRange,indexForHighRange,indexForLowRange);
  for(int iii = 1; iii < a_order; iii++){
    functString += Form("+[%d]*PhysMath::legendre(%d,(2.0*x-[%d]-[%d])/fabs([%d]-[%d]))", iii,iii,indexForLowRange,indexForHighRange,indexForHighRange,indexForLowRange);
  }
  return functString;
}






/*
void PhysMath::initializeFourierTrigFunction(TF1* a_funct, TGraphAsymmErrors* a_graph, double a_lowRange, double a_highRange){
  double lowRange = a_lowRange;
  double highRange = a_highRange;
  if(a_lowRange != -999.999) lowRange = a_lowRange;
  if(a_highRange != -999.999) highRange = a_highRange;
  int order = (a_funct->GetNpar()-2)/2;
  int average = 0.0;

  //find initialization to parameters

}*/







TF1* PhysMath::getPolynomialFunction(string a_name, int a_order){
  //Parameters: A0 A1*x A2*x^2 ...
  if(a_order < 0){
    cout << "ERROR: PhysMath::getPolynomialFunction: you requested an order < 0. I don't know what to do."<< endl;
    return new TF1(a_name.c_str(),"0.0",-1.0,1.0);
  }

  string functString = "[0]";
  int indexForXScale = 2*a_order + 1;
  for(int iii = 1; iii <= a_order; iii++){
    functString += Form("+[%d]",iii);
    for(int jjj = 1; jjj <= iii; jjj++){
      functString += "*x";
    }
  }

  TF1* funct = new TF1(a_name.c_str(),functString.c_str(),-1,1);
  for(int iii = 0; iii <= a_order; iii++){
    funct->SetParName(iii,Form("A_%d",iii));
    funct->SetParameter(iii,1.0);
  }
  return funct;
}


double PhysMath::blastWaveOperator(double *x, double *par){

  double mTm0  = x[0]; //mT-m0
  double Amp   = par[0]; //Amplitude
  double T     = par[1]; //Tkin
  double betaS = par[2]; //Surface Velocity
  double n     = par[3]; //Power in the Tranverse Velocity factor (betaT)
  double R     = par[4]; //Upper Limit of Radial Integral (SHOULD BE FIXED)
  double mass  = par[5]; //Should be FIXED


  //######  Radial Integral for the Blast Wave Model  ########
  //Integration Variable
  double r = 0.0; //The lower limit of integration is 0.0

  //Integration Steps
  double sum = 0.0;
  double stepWidth = 0.01;
  int    nSteps = (int)(fabs(R)/stepWidth);

  //Compute mT and pT
  double mT = mTm0 + mass;
  double pT = PhysMath::pT_from_mTm0(mTm0, mass);

  //Do the numerical integral
  double rho = 0;
  double BesselI0 = 0;
  double BesselK1 = 0;
  for (int iStep=0; iStep<nSteps; iStep++){
    rho   = TMath::ATanH(betaS * pow(r/R,n));
    BesselI0 = TMath::BesselI0((pT/T)*TMath::SinH(rho));
    BesselK1 = TMath::BesselK1((mT/T)*TMath::CosH(rho));
    sum += r * mT * BesselI0 * BesselK1 * stepWidth;
    r += stepWidth;
  }//End Loop Over iStep

  return (1.0/(2.0*TMath::Pi())) * Amp * sum / T;
}

double PhysMath::blastWave(double a_mTm0, double a_A, double a_T, double a_BetaS, double a_n, double a_R, double a_mass){

  //######  Radial Integral for the Blast Wave Model  ########
  //Integration Variable
  double r = 0.0; //The lower limit of integration is 0.0

  //Integration Steps
  double sum = 0.0;
  double stepWidth = 0.01;
  int    nSteps = (int)(fabs(a_R)/stepWidth);

  //Compute mT and pT
  double mT = a_mTm0 + a_mass;
  double pT = PhysMath::pT_from_mTm0(a_mTm0, a_mass);

  //Do the numerical integral
  double rho = 0;
  double BesselI0 = 0;
  double BesselK1 = 0;
  for (int iStep=0; iStep<nSteps; iStep++){
    rho      = TMath::ATanH(a_BetaS * pow(r/a_R,a_n));
    BesselI0 = TMath::BesselI0((pT/a_T)*TMath::SinH(rho));
    BesselK1 = TMath::BesselK1((mT/a_T)*TMath::CosH(rho));
    sum += r * mT * BesselI0 * BesselK1 * stepWidth;
    r += stepWidth;
  }//End Loop Over iStep

  return (1.0/(2.0*TMath::Pi())) * a_A * sum / a_T;
}

double PhysMath::blastWave_NotInvariant(double a_mTm0, double a_A, double a_T, double a_BetaS, double a_n, double a_R, double a_mass){
  //######  Radial Integral for the Blast Wave Model  ########
  //Integration Variable
  double r = 0.0; //The lower limit of integration is 0.0

  //Integration Steps
  double sum = 0.0;
  double stepWidth = 0.01;
  int    nSteps = (int)(fabs(a_R)/stepWidth);

  //Compute mT and pT
  double mT = a_mTm0 + a_mass;
  double pT = PhysMath::pT_from_mTm0(a_mTm0, a_mass);

  //Do the numerical integral
  double rho = 0;
  double BesselI0 = 0;
  double BesselK1 = 0;
  for (int iStep=0; iStep<nSteps; iStep++){
    rho   = TMath::ATanH(a_BetaS * pow(r/a_R,a_n));
    BesselI0 = TMath::BesselI0((pT/a_T)*TMath::SinH(rho));
    BesselK1 = TMath::BesselK1((mT/a_T)*TMath::CosH(rho));
    sum += r * mT * BesselI0 * BesselK1 * stepWidth;
    r += stepWidth;
  }//End Loop Over iStep

  return mT * a_A * sum / a_T;
}

double PhysMath::blastWave_betaConstr_NotInvariant(double a_mTm0, double a_A, double a_T, double a_BetaS_p4, double a_BetaS_shift, double a_BetaS_p2, double a_BetaS_const, double a_n, double a_R, double a_mass, double a_rapidity){
  //######  Radial Integral for the Blast Wave Model  ########
  //Integration Variable
  double r = 0.0; //The lower limit of integration is 0.0

  //Integration Steps
  double sum = 0.0;
  double stepWidth = 0.01;
  int    nSteps = (int)(fabs(a_R)/stepWidth);

  //Compute mT and pT
  double mT = a_mTm0 + a_mass;
  double pT = PhysMath::pT_from_mTm0(a_mTm0, a_mass);

  double betaS = a_BetaS_p4*pow(a_rapidity - a_BetaS_shift,4) + a_BetaS_p2*pow(a_rapidity - a_BetaS_shift,2) + a_BetaS_const;

  //Do the numerical integral
  double rho = 0;
  double BesselI0 = 0;
  double BesselK1 = 0;
  for (int iStep=0; iStep<nSteps; iStep++){
    rho   = TMath::ATanH(betaS * pow(r/a_R,a_n));
    BesselI0 = TMath::BesselI0((pT/a_T)*TMath::SinH(rho));
    BesselK1 = TMath::BesselK1((mT/a_T)*TMath::CosH(rho));
    sum += r * mT * BesselI0 * BesselK1 * stepWidth;
    r += stepWidth;
  }//End Loop Over iStep

  return mT * a_A * sum / a_T;
}

TF1* PhysMath::blastWaveFunct(string a_name, double a_A, double a_T, double a_BetaS, double a_n, double a_R, double a_mass){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::blastWave(x,[0],[1],[2],[3],[4],[5])",0,20);
  funct->SetParameter(0,a_A);
  funct->SetParameter(1,a_T);
  funct->SetParameter(2,a_BetaS);
  funct->SetParameter(3,a_n);
  funct->SetParameter(4,a_R);
  funct->SetParameter(5,a_mass);
  funct->SetParNames("A","T_kin","#beta_S","n","R","m_p");
  return funct;
}

TF1* PhysMath::blastWaveFunct_NotInvariant(string a_name, double a_A, double a_T, double a_BetaS, double a_n, double a_R, double a_mass){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::blastWave_NotInvariant(x,[0],[1],[2],[3],[4],[5])",0,20);
  funct->SetParameter(0,a_A);
  funct->SetParameter(1,a_T);
  funct->SetParameter(2,a_BetaS);
  funct->SetParameter(3,a_n);
  funct->SetParameter(4,a_R);
  funct->SetParameter(5,a_mass);
  funct->SetParNames("A","T_kin","#beta_S","W","R","m_p");
  return funct;
}

TF1* PhysMath::blastWaveFunct_betaConstr_NotInvariant(string a_name, double a_A, double a_T, double a_BetaS_p4, double a_BetaS_shift, double a_BetaS_p2, double a_BetaS_const, double a_n, double a_R, double a_mass, double a_rapidity){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::blastWave_betaConstr_NotInvariant(x,[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])",0,20);
  funct->SetParameter(0,a_A);
  funct->SetParameter(1,a_T);
  funct->SetParameter(2,a_BetaS_p4);
  funct->SetParameter(3,a_BetaS_shift);
  funct->SetParameter(4,a_BetaS_p2);
  funct->SetParameter(5,a_BetaS_const);
  funct->SetParameter(6,a_n);
  funct->SetParameter(7,a_R);
  funct->SetParameter(8,a_mass);
  funct->SetParameter(9,a_rapidity);
  funct->SetParNames("A","T_kin","p4 #beta_S","shift #beta_S","p2 #beta_S","const #beta_S","W","R","m_p","y");
  return funct;
}


TF1* PhysMath::blastWaveFunctPlusFitSafeExponentialFunct(string a_name,double a_A, double a_T, double a_BetaS, double a_n, double a_mass, double a_scale, double a_exponent){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::blastWave(x,[0],[1],[2],[3],1.0,[4]) + [5]*exp(-x/[6])/[6]",0,20);
  funct->SetParameter(0,a_A);
  funct->SetParameter(1,a_T);
  funct->SetParameter(2,a_BetaS);
  funct->SetParameter(3,a_n);
  funct->SetParameter(4,a_mass);
  funct->SetParameter(5,a_scale);
  funct->SetParameter(6,a_exponent);
  funct->SetParLimits(0,0.0,1.0);
  funct->SetParLimits(1,0.0,10.0);
  funct->SetParLimits(2,0.0,1.0);
  funct->SetParLimits(3,0.0001,15.0);
  funct->SetParLimits(4,0.0,350.0);
  funct->SetParLimits(5,0.0,1.0);
  funct->SetParLimits(6,0.0,50.0);
  funct->SetParNames("A","T","#beta_S","n","m_p","N_2","T_2");
  return funct;
}





double PhysMath::boseEinstein_CoulombCorrected_NotInvariant_Plus(double a_mTm0, double a_A, double a_T, double a_mass, double a_VC, double a_rapVal){


  //Compute mT and pT
  double mT = a_mTm0 + a_mass;
  double mom = PhysMath::pTot_from_y_mTm0(a_rapVal,a_mTm0,a_mass);
  double Ef = mT*TMath::CosH(a_rapVal);
  double Tp = 0.12;
  double mPro = 0.9382723128;
  double Emax = sqrt(pow(mPro*mom/a_mass,2)+pow(mPro,2))-mPro;
  double Veff = a_VC*(1-exp(-Emax/Tp));

  return 2*TMath::Pi()*mT*(mom*(Ef-Veff)/Ef/sqrt(pow(Ef-Veff,2)-a_mass*a_mass))*a_A/(exp((Ef-Veff)/a_T)-1.0);

}

double PhysMath::boseEinstein_CoulombCorrected_NotInvariant_Minus(double a_mTm0, double a_A, double a_T, double a_mass, double a_VC, double a_rapVal){


  //Compute mT and pT
  double mT = a_mTm0 + a_mass;
  double mom = PhysMath::pTot_from_y_mTm0(a_rapVal,a_mTm0,a_mass);
  double Ef = mT*TMath::CosH(a_rapVal);
  double Tp = 0.12;
  double mPro = 0.9382723128;
  double Emax = sqrt(pow(mPro*mom/a_mass,2)+pow(mPro,2))-mPro;
  double Veff = a_VC*(1-exp(-Emax/Tp));

  return 2*TMath::Pi()*mT*(mom*(Ef+Veff)/Ef/sqrt(pow(Ef+Veff,2)-a_mass*a_mass))*a_A/(exp((Ef+Veff)/a_T)-1.0);

}

TF1* PhysMath::boseEinsteinFunct_CoulombCorrected_NotInvariant(string a_name, double a_A, double a_T, double a_mass, double a_VC, double a_rapVal, int a_charge){
  TF1* funct = NULL;
  if(a_charge == 1) funct = new TF1(a_name.c_str(),"PhysMath::boseEinstein_CoulombCorrected_NotInvariant_Plus(x,[0],[1],[2],[3],[4])",0,20);
  else if(a_charge == -1) funct = new TF1(a_name.c_str(),"PhysMath::boseEinstein_CoulombCorrected_NotInvariant_Minus(x,[0],[1],[2],[3],[4])",0,20);
  funct->SetParameter(0,a_A);
  funct->SetParameter(1,a_T);
  funct->SetParameter(2,a_mass);
  funct->SetParameter(3,a_VC);
  funct->SetParameter(4,a_rapVal);
  funct->SetParNames("A_{Low}","T_{Low}","m_{#pi}","V_{C}","y_{lab}");
  return funct;
}


TF1* PhysMath::getCeresDnDyFunction(string a_name, double a_dNdy_0, double a_sqrt_sNN, double a_mass){
  TF1* funct = new TF1(a_name.c_str(), "[0]*pow(cosh(0.75*x/(sqrt(log(0.5*[1]/[2])*(1.0-0.5*x*x*[2]/[1])))),-2)", -20,20 );
  funct->SetParNames("dNdy_0","#sqrt{s_{NN}}","m");
  funct->SetParameters(a_dNdy_0,a_sqrt_sNN,a_mass);
  return funct;
}




double PhysMath::doubleBoseEinstein_CoulombCorrected_NotInvariant_Plus(double a_mTm0, double a_ALow, double a_TALow, double a_TmLow, double a_TsLow, double a_VC, double a_AHigh, double a_TAHigh, double a_TsHigh, double a_rapVal, double a_centIndex){

  int centIndex = (int) a_centIndex;

  double Tpro_amp[7] = {0.252297,0.23596,0.20449,0.188534,0.163917,0.13915,0.127722};
  double Tpro_sig[7] = {1.57028,1.44414,1.37874,1.43276,1.37272,1.35113,1.46965};

  int tmp = 0;
  if(tmp == 0){
    //cout << "pFit" << endl;
    //cout << "a_rapVal: " << a_rapVal << endl;
    //cout << "a_centIndex: " << a_centIndex << endl;
    tmp = 1;
  }

  //Compute mT and pT
  double mPi = 0.13957;
  double mT = a_mTm0 + mPi;
  double mom = PhysMath::pTot_from_y_mTm0(a_rapVal - a_TmLow,a_mTm0,mPi);
  double Ef = mT*TMath::CosH(a_rapVal - a_TmLow);
  double Tp = Tpro_amp[centIndex]/TMath::CosH(Tpro_sig[centIndex]*(a_rapVal - a_TmLow));
  double mPro = 0.9382723128;
  double Emax = sqrt(pow(mPro*mom/mPi,2)+pow(mPro,2))-mPro;
  double Veff = a_VC/1000.0*(1-exp(-Emax/Tp));

  double TLow = a_TALow/TMath::CosH(a_TsLow*(a_rapVal - a_TmLow));
  double THigh = a_TAHigh/TMath::CosH(a_TsHigh*(a_rapVal - a_TmLow));

  return 2*TMath::Pi()*mT*(Ef*sqrt(pow(Ef-Veff,2)-mPi*mPi)*(Ef-Veff)/(mom*Ef))*(a_ALow/(exp((Ef-Veff)/TLow)-1.0)/(6.283185307179586 * TLow) + a_AHigh/(exp((Ef-Veff)/THigh)-1.0)/(6.283185307179586 * THigh));

}

double PhysMath::doubleBoseEinstein_CoulombCorrected_NotInvariant_Minus(double a_mTm0, double a_ALow, double a_TALow, double a_TmLow, double a_TsLow, double a_VC, double a_AHigh, double a_TAHigh, double a_TsHigh, double a_rapVal, double a_centIndex){

  int centIndex = (int) a_centIndex;

  double Tpro_amp[7] = {0.252297,0.23596,0.20449,0.188534,0.163917,0.13915,0.127722};
  double Tpro_sig[7] = {1.57028,1.44414,1.37874,1.43276,1.37272,1.35113,1.46965};

  int tmp = 0;
  if(tmp == 0){
    //cout << "mFit" << endl;
    //cout << "a_rapVal: " << a_rapVal << endl;
    //cout << "a_centIndex: " << a_centIndex << endl;
    tmp = 1;
  }

  //Compute mT and pT
  double mPi = 0.13957;
  double mT = a_mTm0 + mPi;
  double mom = PhysMath::pTot_from_y_mTm0(a_rapVal - a_TmLow,a_mTm0,mPi);
  double Ef = mT*TMath::CosH(a_rapVal - a_TmLow);
  double Tp = Tpro_amp[centIndex]/TMath::CosH(Tpro_sig[centIndex]*(a_rapVal - a_TmLow));
  double mPro = 0.9382723128;
  double Emax = sqrt(pow(mPro*mom/mPi,2)+pow(mPro,2))-mPro;
  double Veff = a_VC/1000.0*(1-exp(-Emax/Tp));

  double TLow = a_TALow/TMath::CosH(a_TsLow*(a_rapVal - a_TmLow));
  double THigh = a_TAHigh/TMath::CosH(a_TsHigh*(a_rapVal - a_TmLow));

  return 2*TMath::Pi()*mT*(Ef*sqrt(pow(Ef+Veff,2)-mPi*mPi)*(Ef+Veff)/(mom*Ef))*(a_ALow/(exp((Ef+Veff)/TLow)-1.0)/(6.283185307179586 * TLow) + a_AHigh/(exp((Ef+Veff)/THigh)-1.0)/(6.283185307179586 * THigh));

}

TF1* PhysMath::doubleBoseEinsteinFunct_CoulombCorrected_NotInvariant(string a_name, double a_ALow, double a_TALow, double a_TmLow, double a_TsLow, double a_VC, double a_AHigh, double a_TAHigh, double a_TsHigh, double a_rapVal, int a_charge, double a_centIndex){
  TF1* funct = NULL;

  int tmp = 0;
  if(tmp == 0){
    //cout << "fit Init" << endl;
    //cout << "a_rapVal: " << a_rapVal << endl;
    // cout << "a_centIndex: " << a_centIndex << endl;
    tmp = 1;
  }
  if(a_charge == 1) funct = new TF1(a_name.c_str(),"PhysMath::doubleBoseEinstein_CoulombCorrected_NotInvariant_Plus(x,[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])",0,20);
  else if(a_charge == -1) funct = new TF1(a_name.c_str(),"PhysMath::doubleBoseEinstein_CoulombCorrected_NotInvariant_Minus(x,[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])",0,20);
  funct->SetParameter(0,a_ALow);
  funct->SetParameter(1,a_TALow);
  funct->SetParameter(2,a_TmLow);
  funct->SetParameter(3,a_TsLow);
  funct->SetParameter(4,a_VC);
  funct->SetParameter(5,a_AHigh);
  funct->SetParameter(6,a_TAHigh);
  funct->SetParameter(7,a_TsHigh);
  funct->FixParameter(8,a_rapVal);
  funct->SetParLimits(8,a_rapVal,a_rapVal);
  funct->FixParameter(9,a_centIndex);
  funct->SetParLimits(9,a_centIndex,a_centIndex);
  funct->SetParNames("A_{Low}","A T_{Low}","#mu T_{Low}","#sigma T_{Low}","V_{C} (MeV)","A_{High}","A T_{High}","#sigma T_{High}","y_{lab}","centIndex");
  funct->SetNpx(300);
  return funct;
}



double PhysMath::doubleBoseEinstein_CoulombCorrected_FreeTemps_NotInvariant_Plus(double a_mTm0, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_rapVal, double a_centIndex, double a_yCM){

  int centIndex = (int) a_centIndex;

  double Tpro_amp[7] = {0.252297,0.23596,0.20449,0.188534,0.163917,0.13915,0.127722};
  double Tpro_sig[7] = {1.57028,1.44414,1.37874,1.43276,1.37272,1.35113,1.46965};

  int tmp = 0;
  if(tmp == 0){
    //cout << "pFit" << endl;
    //cout << "a_rapVal: " << a_rapVal << endl;
    //cout << "a_centIndex: " << a_centIndex << endl;
    tmp = 1;
  }

  //Compute mT and pT
  double mPi = 0.13957;
  double mT = a_mTm0 + mPi;
  double mom = PhysMath::pTot_from_y_mTm0(a_rapVal - a_yCM,a_mTm0,mPi);
  double Ef = mT*TMath::CosH(a_rapVal - a_yCM);
  double Tp = Tpro_amp[centIndex]/TMath::CosH(Tpro_sig[centIndex]*(a_rapVal - a_yCM));
  double mPro = 0.9382723128;
  double Emax = sqrt(pow(mPro*mom/mPi,2)+pow(mPro,2))-mPro;
  double Veff = a_VC*(1-exp(-Emax/Tp));

  double THigh = a_Amp_THigh/TMath::CosH(a_Sig_THigh*(a_rapVal - a_yCM));

  return 2*TMath::Pi()*mT*(Ef*sqrt(pow(Ef-Veff,2)-mPi*mPi)*(Ef-Veff)/(mom*Ef))*(a_ALow/(exp((Ef-Veff)/a_TLow)-1.0)/(6.283185307179586 * a_TLow) + a_AHigh/(exp((Ef-Veff)/THigh)-1.0)/(6.283185307179586 * THigh));

}

double PhysMath::doubleBoseEinstein_CoulombCorrected_FreeTemps_NotInvariant_Minus(double a_mTm0, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_rapVal, double a_centIndex, double a_yCM){

  int centIndex = (int) a_centIndex;

  double Tpro_amp[7] = {0.252297,0.23596,0.20449,0.188534,0.163917,0.13915,0.127722};
  double Tpro_sig[7] = {1.57028,1.44414,1.37874,1.43276,1.37272,1.35113,1.46965};

  int tmp = 0;
  if(tmp == 0){
    //cout << "mFit" << endl;
    //cout << "a_rapVal: " << a_rapVal << endl;
    //cout << "a_centIndex: " << a_centIndex << endl;
    tmp = 1;
  }

  //Compute mT and pT
  double mPi = 0.13957;
  double mT = a_mTm0 + mPi;
  double mom = PhysMath::pTot_from_y_mTm0(a_rapVal - a_yCM,a_mTm0,mPi);
  double Ef = mT*TMath::CosH(a_rapVal - a_yCM);
  double Tp = Tpro_amp[centIndex]/TMath::CosH(Tpro_sig[centIndex]*(a_rapVal - a_yCM));
  double mPro = 0.9382723128;
  double Emax = sqrt(pow(mPro*mom/mPi,2)+pow(mPro,2))-mPro;
  double Veff = a_VC*(1-exp(-Emax/Tp));

  double THigh = a_Amp_THigh/TMath::CosH(a_Sig_THigh*(a_rapVal - a_yCM));

  return 2*TMath::Pi()*mT*(Ef*sqrt(pow(Ef+Veff,2)-mPi*mPi)*(Ef+Veff)/(mom*Ef))*(a_ALow/(exp((Ef+Veff)/a_TLow)-1.0)/(6.283185307179586 * a_TLow) + a_AHigh/(exp((Ef+Veff)/THigh)-1.0)/(6.283185307179586 * THigh));

}

TF1* PhysMath::doubleBoseEinsteinFunct_CoulombCorrected_FreeTemps_NotInvariant(string a_name, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_rapVal, int a_charge, double a_centIndex, double a_yCM){
  TF1* funct = NULL;

  int tmp = 0;
  if(tmp == 0){
    //cout << "fit Init" << endl;
    //cout << "a_rapVal: " << a_rapVal << endl;
    // cout << "a_centIndex: " << a_centIndex << endl;
    tmp = 1;
  }
  if(a_charge == 1) funct = new TF1(a_name.c_str(),"PhysMath::doubleBoseEinstein_CoulombCorrected_FreeTemps_NotInvariant_Plus(x,[0],[1],[2],[3],[4],[5],[6],[7],[8])",0,20);
  else if(a_charge == -1) funct = new TF1(a_name.c_str(),"PhysMath::doubleBoseEinstein_CoulombCorrected_FreeTemps_NotInvariant_Minus(x,[0],[1],[2],[3],[4],[5],[6],[7],[8])",0,20);
  funct->SetParameter(0,a_ALow);
  funct->SetParameter(1,a_TLow);
  funct->SetParameter(2,a_VC);
  funct->SetParameter(3,a_AHigh);
  funct->SetParameter(4,a_Amp_THigh);
  funct->SetParameter(5,a_Sig_THigh);
  funct->FixParameter(6,a_rapVal);
  funct->SetParLimits(6,a_rapVal,a_rapVal);
  funct->FixParameter(7,a_centIndex);
  funct->SetParLimits(7,a_centIndex,a_centIndex);
  funct->FixParameter(8,a_yCM);
  funct->SetParLimits(8,a_yCM,a_yCM);
  funct->SetParNames("A_{Low}","T_{Low}","V_{C} (GeV)","A_{High}","Amp T_{High}","Sig T_{High}","y_{lab}","centIndex","y_{CM}");
  funct->SetNpx(300);
  return funct;
}



double PhysMath::doubleBoseEinstein_CoulombCorrected_FreeTemps_noMem_NotInvariant_Plus(double a_mTm0, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_proTemp, double a_rapCM){


  //Compute mT and pT
  double mPi = 0.13957;
  double mT = a_mTm0 + mPi;
  double mom = PhysMath::pTot_from_y_mTm0(a_rapCM,a_mTm0,mPi);
  double Ef = mT;
  //double Tp = Tpro_amp[centIndex]/TMath::CosH(Tpro_sig[centIndex]*(a_rapVal - a_yCM));
  double mPro = 0.9382723128;
  double Emax = sqrt(pow(mPro*mom/mPi,2)+pow(mPro,2))-mPro;
  double Veff = a_VC*(1-exp(-Emax/a_proTemp));

  double THigh = a_Amp_THigh/TMath::CosH(a_Sig_THigh*(a_rapCM));

  return 2*TMath::Pi()*mT*(Ef*sqrt(pow(Ef-Veff,2)-mPi*mPi)*(Ef-Veff)/(mom*Ef))*(a_ALow/(exp((Ef-Veff)/a_TLow)-1.0)/(6.283185307179586 * a_TLow) + a_AHigh/(exp((Ef-Veff)/THigh)-1.0)/(6.283185307179586 * THigh));

}

double PhysMath::doubleBoseEinstein_CoulombCorrected_FreeTemps_noMem_NotInvariant_Minus(double a_mTm0, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_proTemp, double a_rapCM){

  //Compute mT and pT
  double mPi = 0.13957;
  double mT = a_mTm0 + mPi;
  double mom = PhysMath::pTot_from_y_mTm0(a_rapCM,a_mTm0,mPi);
  double Ef = mT;
  //double Tp = Tpro_amp[centIndex]/TMath::CosH(Tpro_sig[centIndex]*(a_rapVal - a_yCM));
  double mPro = 0.9382723128;
  double Emax = sqrt(pow(mPro*mom/mPi,2)+pow(mPro,2))-mPro;
  double Veff = a_VC*(1-exp(-Emax/a_proTemp));

  double THigh = a_Amp_THigh/TMath::CosH(a_Sig_THigh*(a_rapCM));

  return 2*TMath::Pi()*mT*(Ef*sqrt(pow(Ef+Veff,2)-mPi*mPi)*(Ef+Veff)/(mom*Ef))*(a_ALow/(exp((Ef+Veff)/a_TLow)-1.0)/(6.283185307179586 * a_TLow) + a_AHigh/(exp((Ef+Veff)/THigh)-1.0)/(6.283185307179586 * THigh));

}

TF1* PhysMath::doubleBoseEinsteinFunct_CoulombCorrected_FreeTemps_noMem_NotInvariant(string a_name, double a_ALow, double a_TLow, double a_VC, double a_AHigh, double a_Amp_THigh, double a_Sig_THigh, double a_proTemp, double a_rapCM, int a_charge){
  TF1* funct = NULL;

  if(a_charge == 1) funct = new TF1(a_name.c_str(),"PhysMath::doubleBoseEinstein_CoulombCorrected_FreeTemps_noMem_NotInvariant_Plus(x,[0],[1],[2],[3],[4],[5],[6],[7])",0,20);
  else if(a_charge == -1) funct = new TF1(a_name.c_str(),"PhysMath::doubleBoseEinstein_CoulombCorrected_FreeTemps_noMem_NotInvariant_Minus(x,[0],[1],[2],[3],[4],[5],[6],[7])",0,20);
  funct->SetParameter(0,a_ALow);
  funct->SetParameter(1,a_TLow);
  funct->SetParameter(2,a_VC);
  funct->SetParameter(3,a_AHigh);
  funct->SetParameter(4,a_Amp_THigh);
  funct->SetParameter(5,a_Sig_THigh);
  funct->FixParameter(6,a_proTemp);
  funct->SetParLimits(6,a_proTemp,a_proTemp);
  funct->FixParameter(7,a_rapCM);
  funct->SetParLimits(7,a_rapCM,a_rapCM);
  funct->SetParNames("A_{Low}","T_{Low}","V_{C} (GeV)","A_{High}","Amp T_{High}","Sig T_{High}","T_{pro}","y_{CM}");
  funct->SetNpx(300);
  return funct;
}


TF1* PhysMath::parabolaFunct(string a_name, double a_yLoc, double a_ySlope, double a_ySecDer, double a_xLoc){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::parabolaByDerivatives(x,[0],[1],[2],[3])",-30,30);
  funct->SetParameter(0,a_yLoc);
  funct->SetParameter(1,a_ySlope);
  funct->SetParameter(2,a_ySecDer);
  funct->FixParameter(3,a_xLoc);
  funct->SetParNames("y_0","y'","y''","x_0");
  return funct;

}

TF1* PhysMath::parabolaFunct(string a_name, double a_x1, double a_y1, double a_x2, double a_y2, double a_x3, double a_y3){
  double x0 = a_x2;
  double y0 = -1;
  double ySlope = -1;
  double ySecDer = -1;
  PhysMath::get3PointParabParams( x0, y0, ySlope, ySecDer,  a_x1, a_y1, a_x2, a_y2, a_x3, a_y3);
  return PhysMath::parabolaFunct(a_name,  y0,  ySlope,  ySecDer,  x0);
}




TF1* PhysMath::linearFunct(string a_name, double a_yLoc, double a_ySlope, double a_xLoc){
  //TF1* funct = new TF1(a_name.c_str(),"PhysMath::lineByDerivatives(x,[0],[1],[2])",-30,30);
  TF1* funct = new TF1(a_name.c_str(),"[1]*(x-[2])+[0]",-30,30);
  funct->SetParameter(0,a_yLoc);
  funct->SetParameter(1,a_ySlope);
  funct->FixParameter(2,a_xLoc);
  funct->SetParNames("y_0","y'","x_0");
  return funct;

}

TF1* PhysMath::linearFunct(string a_name, double a_x1, double a_y1, double a_x2, double a_y2){
  double slope = (a_y1-a_y2)/(a_x1-a_x2);
  return PhysMath::linearFunct(a_name, a_y1, slope, a_x1);
}


TF1* PhysMath::piecewiseTurnOnFunct(string a_name, double a_A, double a_B, double a_C, double a_xCut, double a_ySlopeCut, double a_m){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::piecewiseTurnOn([0],[1],[2],[3],[4],[5],x)",0,6);
  funct->SetParNames("A","b","c","x_c","dy/dx_c","m");
  funct->SetParameter(0,a_A);
  funct->SetParameter(1,a_B);
  funct->SetParameter(2,a_C);
  funct->SetParameter(3,a_xCut);
  funct->SetParameter(4,a_ySlopeCut);
  funct->SetParameter(5,a_m);
  return funct;
}

TF1* PhysMath::piecewiseTurnOnTwoLinesFunct(string a_name,double a_A, double a_B, double a_C,double a_x0_cut, double a_x1_cut, double a_y1_cut,  double a_ySlope1, double a_m){
  TF1* funct = new TF1(a_name.c_str(),"PhysMath::piecewiseTurnOnTwoLines([0],[1],[2],[3],[4],[5],[6],[7],x)",0,6);
  funct->SetParNames("A","b","c","x_0","x_1","y_1","dy/dx_1","m");
  funct->SetParameter(0,a_A);
  funct->SetParameter(1,a_B);
  funct->SetParameter(2,a_C);
  funct->SetParameter(3,a_x0_cut);
  funct->SetParameter(4,a_x1_cut);
  funct->SetParameter(5,a_y1_cut);
  funct->SetParameter(6,a_ySlope1);
  funct->SetParameter(7,a_m);
  return funct;
}

TF1* PhysMath::getInverseCoshFunction(string a_name, double a_Amp, double a_mu, double a_width, bool a_fixMu){
  TF1* funct = new TF1(a_name.c_str(),"[0]/TMath::CosH([2]*(x-[1]))",-3.0,3.0);
  funct->SetParameter(0,a_Amp);
  funct->SetParLimits(0,0.0,10.0*a_Amp);
  if(a_fixMu){
    funct->FixParameter(1,a_mu);
  }else{
    funct->SetParameter(1,a_mu);
    funct->SetParLimits(1,a_mu-0.1*a_width,a_mu+0.1*a_width);
  }
  funct->SetParameter(2,a_width);
  funct->SetParLimits(2,0.1*a_width,15.0*a_width);
  funct->SetNpx(300);
  return funct;
}




TF1* PhysMath::getLogisticFunction(){
  TF1* funct = new TF1("LogisticFunction",
    "[0]*( 2.0 / (1.0 + exp(-[1]*(x - [2]))) - 1 ) + [3]",
    -5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlue);
  funct->SetParNames("A","m","x_{0}","y_{0}");
  funct->SetParameter(0, 1.0);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 1.0);
  funct->SetParLimits(1,0.0,1000.0);
  funct->SetParameter(2, 0.0);
  funct->SetParLimits(2, -2.0, 2.0);
  funct->SetParameter(3, 0.0);
  funct->SetParLimits(3, -2.0, 2.0);
  return funct;
}

TF1* PhysMath::getGeneralizedLogisticFunction(){
  /*TF1* funct = new TF1("GeneralizedLogisticFunction",
    "[0]*(2.0*pow( 1.0 + exp(-[1]*(x - [2] + log(pow(2.0,1.0/[4]) - 1.0)/[1])) ,-[4]) - 1.0) - [3]",
    -5.0,5.0);

  funct->SetNpx(500);
  funct->SetLineColor(kRed);
  funct->SetParNames("A","m","x_{0}","y_{0}","#alpha");
  funct->SetParameter(0, 1.0);
  funct->SetParameter(1, 1.0);
  funct->SetParLimits(1,0.0,1000.0);
  funct->SetParameter(2, 0.0);
  funct->SetParameter(3, 0.0);
  funct->SetParameter(4, 1.0);
  funct->SetParLimits(4,0.0,1000.0); */

  TF1* funct = new TF1("GeneralizedLogisticFunction",
    "[0]*( pow(1.0 + [5]*exp(-[1]*(x-[2])) ,-1.0/[4] ) - pow(1.0+[5],-1.0/[4]) )  /  (1.0 - pow(1.0+[5],-1.0/[4])) + [3]",
    -5.0,5.0);

  funct->SetNpx(500);
  funct->SetLineColor(kRed);
  funct->SetParNames("A","B","x_{0}","y_{0}","#nu","Q");
  funct->SetParameter(0, 1.0);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 1.0);
  funct->SetParLimits(1,0.0001,100.0);
  funct->SetParameter(2, 0.0);
  funct->SetParLimits(2, -2.0, 2.0);
  funct->SetParameter(3, 0.0);
  funct->SetParLimits(3, -2.0, 2.0);
  funct->SetParameter(4, 1.0);
  funct->SetParLimits(4,0.0,100.0);
  funct->SetParameter(5, 1.0);
  funct->SetParLimits(5,1E-6,100.0);
  return funct;
}

TF1* PhysMath::getHyperbolicTangentFunction(){
  TF1* funct = new TF1("HyperbolicTangentFunction",
    "[0]*TMath::TanH([1]*(x - [2])) + [3]",
    -5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kGreen);
  funct->SetParNames("A","m","x_{0}","y_{0}");
  funct->SetParameter(0, 1.0);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 1.0);
  funct->SetParLimits(1,0.0,1000.0);
  funct->SetParameter(2, 0.0);
  funct->SetParLimits(2, -2.0, 2.0);
  funct->SetParameter(3, 0.0);
  funct->SetParLimits(3, -2.0, 2.0);
  return funct;
}

TF1* PhysMath::getArcTangentFunction(){
  TF1* funct = new TF1("ArcTangentFunciton",
    "[0]*2.0/TMath::Pi()*TMath::ATan(TMath::Pi() / 2.0 * [1]*(x - [2])) + [3]",
    -5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kViolet);
  funct->SetParNames("A","m","x_{0}","y_{0}");
  funct->SetParameter(0, 1.0);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 1.0);
  funct->SetParLimits(1,0.0,1000.0);
  funct->SetParameter(2, 0.0);
  funct->SetParLimits(2, -2.0, 2.0);
  funct->SetParameter(3, 0.0);
  funct->SetParLimits(3, -2.0, 2.0);
  return funct;
}

TF1* PhysMath::getGundimannianFunction(){
  TF1* funct = new TF1("GundimannianFunction",
    "4.0/TMath::Pi()*[0]*TMath::ATan( TMath::TanH( TMath::Pi()/4.0*[1]*(x-[2]) ) ) + [3]",
    -5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kOrange);
  funct->SetParNames("A","m","x_{0}","y_{0}");
  funct->SetParameter(0, 1.0);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 1.0);
  funct->SetParLimits(1,0.0,1000.0);
  funct->SetParameter(2, 0.0);
  funct->SetParLimits(2, -2.0, 2.0);
  funct->SetParameter(3, 0.0);
  funct->SetParLimits(3, -2.0, 2.0);
  return funct;
}

TF1* PhysMath::getErrorFunction(){
  TF1* funct = new TF1("ErrorFunction",
    "[0] * TMath::Erf( sqrt(TMath::Pi())/2.0 * [1]*( x - [2]) ) + [3]",
    -5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlack);
  funct->SetParNames("A","m","x_{0}","y_{0}");
  funct->SetParameter(0, 1.0);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 1.0);
  funct->SetParLimits(1,0.0,1000.0);
  funct->SetParameter(2, 0.0);
  funct->SetParLimits(2, -2.0, 2.0);
  funct->SetParameter(3, 0.0);
  funct->SetParLimits(3, -2.0, 2.0);
  return funct;
}

TF1* PhysMath::getAlgebraicTurnOnFunction(){
  TF1* funct = new TF1("AlgebraicTurnOn",
    "[0]*[1]*(x - [2])*pow(1.0 + pow([1]*(x - [2]) , 2.0),-0.5) + [3]",
    -5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlue + 1);
  funct->SetParNames("A","m","x_{0}","y_{0}");
  funct->SetParameter(0, 1.0);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 1.0);
  funct->SetParLimits(1,0.0,1000.0);
  funct->SetParameter(2, 0.0);
  funct->SetParLimits(2, -2.0, 2.0);
  funct->SetParameter(3, 0.0);
  funct->SetParLimits(3, -2.0, 2.0);
  return funct;
}

TF1* PhysMath::getAlgebraicTurnOnFunction_2(){
  TF1* funct = new TF1("AlgebraicTurnOn_2",
    "[0]*[1]*(x - [2]) / ( 1.0 + fabs([1]*(x - [2])) ) + [3]",
    -5.0,5.0);
  funct->SetNpx(500);
  funct->SetLineColor(kBlue + 2);
  funct->SetParNames("A","m","x_{0}","y_{0}");
  funct->SetParameter(0, 1.0);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 1.0);
  funct->SetParLimits(1,0.0,1000.0);
  funct->SetParameter(2, 0.0);
  funct->SetParLimits(2, -2.0, 2.0);
  funct->SetParameter(3, 0.0);
  funct->SetParLimits(3, -2.0, 2.0);
  return funct;
}


TF1* PhysMath::getBichselFunction(string a_functionName, double a_mass){
  //H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195
  // x is momentum
  TF1* funct = new TF1(a_functionName.c_str(),
    "pow(sqrt([5]*[5] + x*x)/x,[3]) * ([0] - [1]*log([2] + pow([5]/x,[4]))) - [1]",
    0.0,20);
  //(beta)^-D * (A - B ln(C + betagamma^-E)) - B
  funct->SetParNames("A","B","C","D","E","m");
  funct->SetParameter(0,0.76);// A
  funct->SetParameter(1,0.076);// B
  funct->SetParameter(2,1.3E-5); //C
  funct->SetParameter(3,1.863); //D
  funct->SetParameter(4,2.25868);// E
  funct->SetParameter(5,a_mass); // m

  funct->SetNpx(300);
  return funct;

}

TF1* PhysMath::getLnBichselFunction(string a_functionName,double a_mass){
  //H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195
  TF1* funct = new TF1(a_functionName.c_str(),
    "log(pow(sqrt([5]*[5] + x*x)/x,[3]) * ([0] - [1]*log([2] + pow([5]/x,[4]))) - [1])",
    0.0,20);
  //   ln(   (beta)^-D * (A - B ln(C + betagamma^-E)) - B   )
  funct->SetParNames("A","B","C","D","E","m");
  funct->SetParameter(0,0.76);// A
  funct->SetParameter(1,0.076);// B
  funct->SetParameter(2,1.3E-5); //C
  funct->SetParameter(3,1.863); //D
  funct->SetParameter(4,2.25868);// E
  funct->SetParameter(5,a_mass); // m

  funct->SetNpx(300);
  return funct;


}


TF1* PhysMath::getBichselFunctionBetaGamma_ChrisFlores(string a_functionName){
  TF1* funct = new TF1(a_functionName.c_str(),
    "[0]*pow(1.0 + 1.0/(x*x),[1]) + [2]*pow(1.0 + 1.0/(x*x),[3])*pow(fabs([4]*x*x),[5]) - [6]",
    0.0,30);
  funct->SetParNames("A","B","C","D","E","F","G");
  funct->SetParameter(0,0.808032);// A
  funct->SetParameter(1,0.76);    // B
  funct->SetParameter(2,1.34E-5); // C
  funct->SetParameter(3,0.932);   // D
  funct->SetParameter(4,0.974);   // E
  funct->SetParameter(5,0.974);   // F
  funct->SetParameter(6,0.974);   // G

  funct->SetParLimits(0,0.0,5.0); // A
  funct->SetParLimits(1,0.0,5.0); // B
  funct->SetParLimits(2,0.0,5.0); // C
  funct->SetParLimits(3,0.0,5.0); // D
  funct->SetParLimits(4,0.0,5.0); // E
  funct->SetParLimits(5,0.0,5.0); // E
  funct->SetParLimits(6,0.0,5.0); // E
  funct->SetNpx(300);
  return funct;
}

TF1* PhysMath::getLnBichselFunctionBetaGamma_ChrisFlores(string a_functionName){
  TF1* funct = new TF1(a_functionName.c_str(),
    "log([0]*pow(1.0 + 1.0/(x*x),[1]) + [2]*pow(1.0 + 1.0/(x*x),[3])*pow(fabs([4]*x*x),[5]) - [6])",
    0.0,30);
  funct->SetParNames("A","B","C","D","E","F","G");
  funct->SetParameter(0,0.808032);// A
  funct->SetParameter(1,0.76);    // B
  funct->SetParameter(2,1.34E-5); // C
  funct->SetParameter(3,0.932);   // D
  funct->SetParameter(4,0.974);   // E
  funct->SetParameter(5,0.974);   // F
  funct->SetParameter(6,0.974);   // G

  funct->SetParLimits(0,0.0,5.0); // A
  funct->SetParLimits(1,0.0,5.0); // B
  funct->SetParLimits(2,0.0,5.0); // C
  funct->SetParLimits(3,0.0,5.0); // D
  funct->SetParLimits(4,0.0,5.0); // E
  funct->SetParLimits(5,0.0,5.0); // E
  funct->SetParLimits(6,0.0,5.0); // E


  funct->SetNpx(300);
  return funct;

}



TF1* PhysMath::getBichselFunctionBetaGamma_Harasty(string a_functionName){
  TF1* funct = new TF1(a_functionName.c_str(),
    "[9]*((1.0-[11])*(1.0-exp(-x/[10])) + [11])*pow([7]*x+[8],-[6]) + [0]*( pow(1.0 + [5]*exp(-[1]*(x-[2])) ,-1.0/[4] ) - pow(1.0+[5],-1.0/[4]) )  /  (1.0 - pow(1.0+[5],-1.0/[4])) + [3] + [12]*x",
    0.0,30);

  funct->SetNpx(500);
  funct->SetLineColor(kRed);
  funct->SetParNames("A","B","x_{0}","y_{0}","#nu","Q","a","b");
  funct->SetParameter(0, 0.743);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 0.4);
  funct->SetParLimits(1,0.0001,100.0);
  funct->SetParameter(2, -2);
  funct->SetParLimits(2, -20.0, 12.0);
  funct->SetParameter(3, 2.0);
  funct->SetParLimits(3, -2.0, 2.0);
  funct->SetParameter(4, 1.43E3);
  funct->SetParLimits(4,0.0,1E5);
  funct->SetParameter(5, 14.6);
  funct->SetParLimits(5,1E-6,100.0);
  funct->SetParameter(6, 2.37);
  funct->SetParLimits(6, 0.0, 100.0);
  funct->SetParameter(7, 0.21);
  funct->SetParLimits(7, 0.0, 10.0);
  funct->SetParameter(8, 1E-12);
  funct->SetParLimits(8,0.0,1E-6);
  funct->SetParameter(9, 0.05);
  funct->SetParLimits(9,0.0,1.0);
  funct->SetParameter(10, 2E-4);
  funct->SetParLimits(10,0.0,1.0);
  funct->SetParameter(11, 0.9);
  funct->SetParLimits(11,0.0,1.5);
  funct->SetParameter(12, 0.0);
  funct->SetParLimits(12,-0.1,5.0);

  return funct;
}


TF1* PhysMath::getLnBichselFunctionBetaGamma_Harasty(string a_functionName){
  TF1* funct = new TF1(a_functionName.c_str(),
    "log([9]*((1.0-[11])*(1.0-exp(-x/[10])) + [11])*pow([7]*x+[8],-[6]) + [0]*( pow(1.0 + [5]*exp(-[1]*(x-[2])) ,-1.0/[4] ) - pow(1.0+[5],-1.0/[4]) )  /  (1.0 - pow(1.0+[5],-1.0/[4])) + [3] + [12]*x)",
    0.0,30);

  funct->SetNpx(500);
  funct->SetLineColor(kRed);
  funct->SetParNames("A","B","x_{0}","y_{0}","#nu","Q","a","b");
  funct->SetParameter(0, 0.743);
  funct->SetParLimits(0, 0.0, 10.0);
  funct->SetParameter(1, 0.4);
  funct->SetParLimits(1,0.0001,100.0);
  funct->SetParameter(2, -2);
  funct->SetParLimits(2, -20.0, 12.0);
  funct->SetParameter(3, 2.0);
  funct->SetParLimits(3, -2.0, 2.0);
  funct->SetParameter(4, 1.43E3);
  funct->SetParLimits(4,0.0,1E5);
  funct->SetParameter(5, 14.6);
  funct->SetParLimits(5,1E-6,100.0);
  funct->SetParameter(6, 2.37);
  funct->SetParLimits(6, 0.0, 100.0);
  funct->SetParameter(7, 0.21);
  funct->SetParLimits(7, 0.0, 10.0);
  funct->SetParameter(8, 1E-12);
  funct->SetParLimits(8,0.0,1E-6);
  funct->SetParameter(9, 0.05);
  funct->SetParLimits(9,0.0,1.0);
  funct->SetParameter(10, 2E-4);
  funct->SetParLimits(10,0.0,1.0);
  funct->SetParameter(11, 0.9);
  funct->SetParLimits(11,0.0,1.5);
  funct->SetParameter(12, 0.0);
  funct->SetParLimits(12,-0.1,5.0);

  return funct;
}

TF1* PhysMath::getBichselFunctionBetaGamma(string a_functionName){
  //H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195
  /*TF1* funct = new TF1(a_functionName.c_str(),
    "pow(sqrt(1.0 + x*x),[3]) * ([0] - [1]*log([2] + pow(x,-[4]))) - [1]",
    0.0,20);

  //(beta)^-D * (A - B ln(C + betagamma^-E)) - B
  funct->SetParNames("A","B","C","D","E");
  funct->SetParameter(0,0.76);// A
  funct->SetParameter(1,0.076);// B
  funct->SetParameter(2,1.3E-5); //C
  funct->SetParameter(3,1.863); //D
  funct->SetParameter(4,2.25868);// E
  */

  /*
  TF1* funct = new TF1(a_functionName.c_str(),
    "[0]*pow(1.0 + 1.0/(x*x),[3])*([1] - pow(1.0 + 1.0/(x*x),-[3]) - log([2]+pow(x,-[4])))",
    0.0,20);
  funct->SetParNames("P1","P2","P3","P4","P5");
  funct->SetParameter(0,0.076);// A
  funct->SetParameter(1,10.632);// B
  funct->SetParameter(2,1.34E-5); //C
  funct->SetParameter(3,1.863); //D
  funct->SetParameter(4,1.948);// E
  */

  TF1* funct = new TF1(a_functionName.c_str(),
    "pow(1.0 + 1.0/(x*x),[3])*([0] - [1]*log([2] + pow(1.0/(x*x),[4]))) - [1]",
    0.0,30);
  funct->SetParNames("A","B","C","D","E","F");
  funct->SetParameter(0,0.808032);// A = P1*P2
  funct->SetParameter(1,0.76);// B = P1
  funct->SetParameter(2,1.34E-5); // C = P3
  funct->SetParameter(3,0.932); // D = P4/2
  funct->SetParameter(4,0.974); // E = P5/2
  //funct->SetParameter(5,0.76); // F = was P1, now new param

  funct->SetParLimits(0,0.0,5.0);// A = P1*P2
  funct->SetParLimits(1,0.0,5.0);// B = P1
  funct->SetParLimits(2,0.0,1.0); // C = P3
  funct->SetParLimits(3,0.0,5.0); // D = P4/2
  funct->SetParLimits(4,0.0,5.0); // E = P5/2
  //funct->SetParLimits(5,0.0,5.0); // F = was P1, now new param


  funct->SetNpx(300);
  return funct;

}

TF1* PhysMath::getLnBichselFunctionBetaGamma(string a_functionName){
  //H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195
  /*TF1* funct = new TF1(a_functionName.c_str(),
    "log(pow(sqrt(1.0 + x*x),[3]) * ([0] - [1]*log([2] + pow(x,-[4]))) - [1])",
    0.0,20);
  //   ln(   (beta)^-D * (A - B ln(C + betagamma^-E)) - B   )
  funct->SetParNames("A","B","C","D","E");
  funct->SetParameter(0,0.76);// A
  funct->SetParameter(1,0.076);// B
  funct->SetParameter(2,1.3E-5); //C
  funct->SetParameter(3,1.863); //D
  funct->SetParameter(4,2.25868);// E
  */

  TF1* funct = new TF1(a_functionName.c_str(),
    "log(pow(1.0 + 1.0/(x*x),[3])*([0] - [1]*log([2] + pow(1.0/(x*x),[4]))) - [1])",
    0.0,30);
  funct->SetParNames("A","B","C","D","E","F");
  funct->SetParameter(0,0.808032);// A = P1*P2
  funct->SetParameter(1,0.76);// B = P1
  funct->SetParameter(2,1.34E-5); // C = P3
  funct->SetParameter(3,0.932); // D = P4/2
  funct->SetParameter(4,0.974); // E = P5/2
  //funct->SetParameter(5,0.76); // F = was P1, now new param

  //funct->SetParLimits(0,0.0,5.0);// A = P1*P2
  //funct->SetParLimits(1,0.0,5.0);// B = P1
  //funct->SetParLimits(2,0.0,1.0); // C = P3
  //funct->SetParLimits(3,0.0,5.0); // D = P4/2
  //funct->SetParLimits(4,0.0,5.0); // E = P5/2
  //funct->SetParLimits(5,0.0,5.0); // F = was P1, now new param

  funct->SetNpx(300);
  funct->SetRange(0.0,20);
  return funct;


}

TF1* PhysMath::getBichselFunctionInvBetaGammaSqr(string a_functionName){
  //H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195
  TF1* funct = new TF1(a_functionName.c_str(),
    "pow(1.0 + x,[3])*([0] - [1]*log([2] + pow(x,[4]))) - [1]",
    0.0,20);
  funct->SetParNames("A","B","C","D","E","F");
  funct->SetParameter(0,0.808032);// A = P1*P2
  funct->SetParameter(1,0.76);// B = P1
  funct->SetParameter(2,1.34E-5); // C = P3
  funct->SetParameter(3,0.932); // D = P4/2
  funct->SetParameter(4,0.974); // E = P5/2
  //funct->SetParameter(5,0.76); // F = was P1, now new param

  funct->SetParLimits(0,0.0,5.0);// A = P1*P2
  funct->SetParLimits(1,0.0,5.0);// B = P1
  funct->SetParLimits(2,0.0,1.0); // C = P3
  funct->SetParLimits(3,0.0,5.0); // D = P4/2
  funct->SetParLimits(4,0.0,5.0); // E = P5/2
  //funct->SetParLimits(5,0.0,5.0); // F = was P1, now new param

  funct->SetNpx(300);
  return funct;
}

TF1* PhysMath::getLnBichselFunctionInvBetaGammaSqr(string a_functionName){
  //H. Bichsel / Nuclear Instruments and Methods in Physics Research A 562 (2006) page 195
  TF1* funct = new TF1(a_functionName.c_str(),
    "log(pow(1.0 + x,[3])*([0] - [1]*log([2] + pow(x,[4]))) - [1])",
    0.0,20);
  funct->SetParNames("A","B","C","D","E","F");
  funct->SetParameter(0,0.808032);// A = P1*P2
  funct->SetParameter(1,0.76);// B = P1
  funct->SetParameter(2,1.34E-5); // C = P3
  funct->SetParameter(3,0.932); // D = P4/2
  funct->SetParameter(4,0.974); // E = P5/2
  //funct->SetParameter(5,0.76); // F = was P1, now new param

  funct->SetParLimits(0,0.0,5.0);// A = P1*P2
  funct->SetParLimits(1,0.0,5.0);// B = P1
  funct->SetParLimits(2,0.0,1.0); // C = P3
  funct->SetParLimits(3,0.0,5.0); // D = P4/2
  funct->SetParLimits(4,0.0,5.0); // E = P5/2
  //funct->SetParLimits(5,0.0,5.0); // F = was P1, now new param

  funct->SetNpx(300);
  return funct;
}




TF1* PhysMath::getInvBetaCalibrationFunction_byMtM0(string a_functionName, double a_mass, double a_mass_parameter, double a_rapidity, double a_timing_offset){//, double a_stopping, double a_stopping_width, double a_scaling){
  TF1* funct = new TF1(a_functionName.c_str(), "sqrt( [0]*[0]/(pow((x + [3])*TMath::CosH([1]),2) - [3]*[3]) + 1.0 ) + [2]",0.0,20);
  //TF1* funct = new TF1(a_functionName.c_str(), "[6]*(1.0 - [4]*exp(-x/[5]))*sqrt( [0]*[0]/(pow((x + [3])*TMath::CosH([1]),2) - [3]*[3]) + 1.0 ) + [2]",0.0,20);
  funct->SetParNames("m","rapidity","#beta^{-1}_{0}","m_{fixed}");
  //funct->SetParNames("m","rapidity","#beta^{-1}_{0}","m_{fixed}","a","b","c");
  funct->SetParameter(0,a_mass_parameter);
  funct->FixParameter(1,a_rapidity);
  funct->SetParameter(2,a_timing_offset);
  funct->FixParameter(3,a_mass);
  funct->SetParLimits(0, 0.0, 10.0);
  //funct->SetParLimits(1,a_rapidity-0.05,a_rapidity+0.05);
  funct->SetParLimits(2, -3.0, 3.0);
  //funct->SetParLimits(3,0.0,10.0);
  /*
  funct->SetParameter(4,a_stopping);
  funct->SetParameter(5,a_stopping_width);
  funct->SetParameter(6,a_scaling);
  funct->SetParLimits(4, 0, 10.0);
  funct->SetParLimits(5, 0, 30.0);
  funct->SetParLimits(6,-10.0, 10.0);
  */

  return funct;
}






void PhysMath::get3PointGausParams(double &a_N, double &a_mu, double &a_sig, double a_x1,double a_y1,double a_x2,double a_y2,double a_x3,double a_y3){

  //log transforms gaussian to parabola
  double p1 = log(a_y1);
  double p2 = log(a_y2);
  double p3 = log(a_y3);

  //solve parabola system of equations: y = A*x^2 + B*x + C
  double denom = (a_x1 - a_x2)*(a_x1 - a_x3)*(a_x2 - a_x3);
  double A = (a_x1*(p3 - p2) + p1*(a_x2 - a_x3) - a_x2*p3 + p2*a_x3)/denom;
  double B = (a_x1*a_x1*(p2 - p3) + p1*(a_x3*a_x3 - a_x2*a_x2) + a_x2*a_x2*p3 - p2*a_x3*a_x3)/denom;
  double C = (a_x1*(a_x1*a_x2*p3 - a_x1*p2*a_x3 - a_x2*a_x2*p3 + p2*a_x3*a_x3) + p1*a_x2*a_x3*(a_x2 - a_x3))/denom;

  //cout << "  check1: " << exp(A*a_x1*a_x1 + B*a_x1 + C) << " =? " << a_y1 << endl;

  // solve for parab params    y = a*(x-x0)^2 + y0
  double a = A;
  double x0 = -0.5*B/A;
  double y0 = C - A*x0*x0;

  if(a > 0.0){
    cout << "ERROR: PhysMath::get3PointGausParams : The three points supplied cannot give a gaussian. The middle point is lower than the outside points." << endl;
    cout << " Arguments : " <<  " p1 " << a_x1 << "  " << a_y1 << "    p2 " << a_x2 << "  " << a_y2 << "     p3 " << a_x3 << "  " << a_y3 << endl;
    a_N = -999;
    a_mu = -999;
    a_sig = -999;
    return;
  }

  if(a_y1 <= 0.0 || a_y2 <= 0.0 || a_y3 <= 0.0){
    cout << "ERROR: PhysMath::get3PointGausParams : The three points supplied cannot give a gaussian. You supplied a y value <= 0" << endl;
    a_N = -999;
    a_mu = -999;
    a_sig = -999;
    return;
  }


  //cout << "  check2: " << exp(a*(a_x1-x0)*(a_x1-x0) + y0) << " =? " << a_y1 << endl;

  // trasform back to gaussian parameters
  a_mu = x0;
  a_sig = sqrt(-0.5/a);
  a_N = 2.5066282746*a_sig*exp(y0);

  //cout << "  check3: " << a_N*exp(-0.5*(a_x1-a_mu)*(a_x1-a_mu)/(a_sig*a_sig))/(sqrt(2.0*3.1415)*a_sig) << " =? " << a_y1 << endl;
  #ifdef _PHYS_MATH_DEBUG_
    cout << "PhysMath::get3PointGausParams" << endl;
    cout << " points (x,y) " << a_x1 << " " << a_y1 << "    " << a_x2 << " " << a_y2 << "    " << a_x3 << " " << a_y3 << endl;
    cout << " logs " << p1 << "  " << p2 << "  " << p3 << endl;
    cout << " denom " << denom << endl;
    cout << "  A: " << A << "  B: " << B << "    C: " << C << endl;
    cout << "  a: " << a << "  x0: " << x0 << "  y0: " << y0 << endl;
    cout << "  N: " << a_N << "  mu: " << a_mu << "  sig: " << a_sig << endl;
  #endif


  return;
}

void PhysMath::get3PointParabParams(double a_x0, double &a_y0, double &a_ySlope, double &a_ySecDer, double a_x1,double a_y1,double a_x2,double a_y2,double a_x3,double a_y3){

  //solve parabola system of equations: y = A*x^2 + B*x + C
  double denom = (a_x1 - a_x2)*(a_x1 - a_x3)*(a_x2 - a_x3);
  double A = (a_x1*(a_y3 - a_y2) + a_y1*(a_x2 - a_x3) - a_x2*a_y3 + a_y2*a_x3)/denom;
  double B = (a_x1*a_x1*(a_y2 - a_y3) + a_y1*(a_x3*a_x3 - a_x2*a_x2) + a_x2*a_x2*a_y3 - a_y2*a_x3*a_x3)/denom;
  double C = (a_x1*(a_x1*a_x2*a_y3 - a_x1*a_y2*a_x3 - a_x2*a_x2*a_y3 + a_y2*a_x3*a_x3) + a_y1*a_x2*a_x3*(a_x2 - a_x3))/denom;

  // solve for parab params    y = A*(x-mu)^2 + b
  a_ySecDer = 2.0*A;
  //double mu = -0.5*B/A;
  //a_ySlope = (a_x0 - mu)*a_ySecDer;
  a_ySlope = (2.0*A*a_x0) + B;
  //double b = C - A*a_x0*a_x0;
  //a_y0 = A*(a_x0 - mu)*(a_x0 - mu) + b;
  a_y0 = A*a_x0*a_x0 + B*a_x0 + C;

  #ifdef _PHYS_MATH_DEBUG_
    cout << "PhysMath::get3PointParabParams" << endl;
    cout << "  check1: " << A*a_x1*a_x1 + B*a_x1 + C << " =? " << a_y1 << endl;
    cout << "  check2: " << 0.5*a_ySecDer*a_x3*a_x3 + (a_ySlope - a_ySecDer*a_x0)*a_x3 + a_y0 + 0.5*a_ySecDer*a_x0*a_x0 - a_ySlope*a_x0  << " =? " << a_y3 << endl;
    cout << "  points (x,y) : (" << a_x1 << " , " << a_y1 << ")    (" << a_x2 << " , " << a_y2 << ")    (" << a_x3 << " , " << a_y3 << ")" << endl;
    cout << "  denom " << denom << endl;
    cout << "  A: " << A << "  B: " << B << "    C: " << C << endl;
    cout << "  y\": " << a_ySecDer  << " y': " << a_ySlope  << "  y_0: " << a_y0 << endl;
  #endif

  return;
}

void PhysMath::gausParams_from_deriviatives(double a_yLoc,double a_ySlope,double a_ySecDer,double a_xLoc, double &a_N, double &a_mu, double &a_sig){
  a_sig = sqrt(-1.0/a_ySecDer);
  a_mu  = a_xLoc + (a_sig*a_sig)*a_ySlope;
  a_N   = sqrt(2.0*TMath::Pi()) * a_sig * exp(a_yLoc + 0.5*(a_xLoc - a_mu)*(a_xLoc - a_mu)/(a_sig*a_sig));
}

void PhysMath::derivatives_from_gausParams(double &a_yLoc,double &a_ySlope,double &a_ySecDer,double a_xLoc, double a_N, double a_mu, double a_sig){
  a_ySecDer = -1.0/(a_sig*a_sig);
  a_ySlope  = -(a_xLoc - a_mu)/(a_sig*a_sig);
  a_yLoc    = log(a_N/(sqrt(2.0*TMath::Pi()) * a_sig)) - 0.5*(a_xLoc - a_mu)*(a_xLoc - a_mu)/(a_sig*a_sig);
}


double PhysMath::gaus_by_deriviatives(double a_x,double a_yLoc,double a_ySlope,double a_ySecDer,double a_xLoc){
  return exp(0.5*a_ySecDer*a_x*a_x + (a_ySlope - a_ySecDer*a_xLoc)*a_x + a_yLoc + 0.5*a_ySecDer*a_xLoc*a_xLoc - a_ySlope*a_xLoc);
}


double PhysMath::parabolaByDerivatives(double a_x, double a_yLoc, double a_ySlope, double a_ySecDer,double a_xLoc){
    //return 0.5*a_ySecDer*a_x*a_x + (a_ySlope - a_ySecDer*a_xLoc)*a_x + a_yLoc + 0.5*a_ySecDer*a_xLoc*a_xLoc - a_ySlope*a_xLoc;

    return 0.5*a_ySecDer*a_x*a_x + (a_ySlope - a_ySecDer*a_xLoc)*a_x + a_yLoc + 0.5*a_ySecDer*a_xLoc*a_xLoc - a_ySlope*a_xLoc;
}

double PhysMath::lineByDerivatives(double a_x, double a_yLoc, double a_ySlope,double a_xLoc){
   return  a_ySlope*(a_x - a_xLoc) + a_yLoc;
}








void PhysMath::testTurnOnFuncts(){

  TCanvas* canv = new TCanvas("TurnOnCanvas");
  TF1* funct1 = getLogisticFunction();
  TF1* funct2 = getGeneralizedLogisticFunction();
  TF1* funct3 = getHyperbolicTangentFunction();
  TF1* funct4 = getArcTangentFunction();
  TF1* funct5 = getGundimannianFunction();
  TF1* funct6 = getErrorFunction();
  TF1* funct7 = getAlgebraicTurnOnFunction();
  TF1* funct8 = getAlgebraicTurnOnFunction_2();

  funct1->Draw();
  funct2->Draw("same");
  funct3->Draw("same");
  funct4->Draw("same");
  funct5->Draw("same");
  funct6->Draw("same");
  funct7->Draw("same");
  funct8->Draw("same");

  return;

}


complex< double > PhysMath::ComplexGamma( complex<double> z){
  // copied from Loic from "https://bytes.com/topic/c/answers/576697-c-routine-complex-gamma-function"
  double p_complexGamma[9] = {
    0.99999999999980993, 676.5203681218851,
    -1259.1392167224028, 771.32342877765313, -176.61502916214059,
    12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
    1.5056327351493116e-7};
  if( real(z) < 0.5 ){
    return TMath::Pi() / (sin(TMath::Pi() * z) * PhysMath::ComplexGamma(1.0 - z));
  }
  z -= 1.0;
  complex<double> x = p_complexGamma[0];
  for(int iii = 1; iii < 9; iii++){
    x += p_complexGamma[iii]/(z + complex<double>(iii,0));
  }
  complex<double> t = z + (7.5);
  return sqrt(2*TMath::Pi()) * pow(t,z + 0.5) * exp(-t) * x;
}



#ifdef _CPP11_
void PhysMath::SimultaneousSortDoubleInt( vector< double > a_values, vector< double > a_other){
  std::sort(a_other.begin(), a_other.end(),[&a_values](size_t i, size_t j){return a_values[i] > a_values[j];});
  std::sort(a_values.begin(), a_values.end(),std::greater<double>());
}
#endif






/*
 Modified by Matt Harasty from
 https://github.com/Expander/polylogarithm/blob/master/src/cpp/Li2.cpp
 Original Author: Alexander Voigt

 Implemented as an economized Pade approximation with a
 maximum error of 4.16e-18.
*/
double PhysMath::Li2(double a_x){
   double y = 0.0, r = 0.0, s = 1.0;

   // transform to [0, 1/2)
   /*if (a_x < -1) {
      const double l = log(1 - a_x);
      y = 1/(1 - a_x);
      r = -PI*PI/6 + l*(0.5*l - log(-a_x));
      s = 1;
   } else if (a_x == -1) {
      return -PI*PI/12;
   } else if (a_x < 0) {
      const double l = log1p(-a_x);
      y = a_x/(a_x - 1);
      r = -0.5*l*l;
      s = -1;
   } else */ if (a_x == 0.0) {
      return 0.0;
   } else if (a_x < 0.5) {
      y = a_x;
      r = 0.0;
      s = 1.0;
   } else if (a_x < 1) {
      y = 1.0 - a_x;
      r = 1.6449340668482262 - log(a_x)*log(1.0 - a_x);
      s = -1.0;
   } else if (a_x == 1.0) {
      return 1.6449340668482262;
   } /*else if (a_x < 2) {
      const double l = log(a_x);
      y = 1 - 1/a_x;
      r = PI*PI/6 - l*(log(1 - 1/a_x) + 0.5*l);
      s = 1;
   } else {
      const double l = log(a_x);
      y = 1/a_x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
   }*/

   double z = y - 0.25;
   double z2 = z*z;
   double z4 = z2*z2;
   double p = physMath_Li2_P[0] + z * physMath_Li2_P[1] + z2 * (physMath_Li2_P[2] + z * physMath_Li2_P[3]) +
              z4 * (physMath_Li2_P[4] + z * physMath_Li2_P[5] + z2 * (physMath_Li2_P[6] + z * physMath_Li2_P[7]));
   double q = physMath_Li2_Q[0] + z * physMath_Li2_Q[1] + z2 * (physMath_Li2_Q[2] + z * physMath_Li2_Q[3]) +
              z4 * (physMath_Li2_Q[4] + z * physMath_Li2_Q[5] + z2 * (physMath_Li2_Q[6] + z * physMath_Li2_Q[7]));
   return r + s*y*p/q;
}



double* PhysMath::choleskyDecomp(double *A, int n, bool a_print){
    double *L = (double*) calloc(n * n, sizeof(double));
    if(L == NULL) exit(EXIT_FAILURE);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < (i+1); j++){
            double s = 0;
            for(int k = 0; k < j; k++) s += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (i == j) ?
                           sqrt(A[i * n + i] - s) :
                           (1.0 / L[j * n + j] * (A[i * n + j] - s));
        }
    if(a_print){
      for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) printf("%3.15f ", L[i * n + j]);
        printf("\n");
      }
    }
    return L;
}


double PhysMath::crossProduct3D(double* a_vect_1,double* a_vect_2){
  double f_1 = (a_vect_1[1]*a_vect_2[2]) - (a_vect_1[2]*a_vect_2[1]);
  double f_2 = (a_vect_1[2]*a_vect_2[0]) - (a_vect_1[0]*a_vect_2[2]);
  double f_3 = (a_vect_1[0]*a_vect_2[1]) - (a_vect_1[1]*a_vect_2[0]);
  return sqrt( (f_1*f_1) + (f_2*f_2) + (f_3*f_3) );
}

void PhysMath::crossProduct3D(double* a_vect_1,double* a_vect_2, double &a_cross_x, double &a_cross_y, double &a_cross_z){
  a_cross_x = (a_vect_1[1]*a_vect_2[2]) - (a_vect_1[2]*a_vect_2[1]);
  a_cross_y = (a_vect_1[2]*a_vect_2[0]) - (a_vect_1[0]*a_vect_2[2]);
  a_cross_z = (a_vect_1[0]*a_vect_2[1]) - (a_vect_1[1]*a_vect_2[0]);
  //cout << "Cross Product: " << endl;
  //cout << "     " << a_vect_1[0] << "  " << a_vect_1[1] << "  " << a_vect_1[2] << endl;
  //cout << "     " << a_vect_2[0] << "  " << a_vect_2[1] << "  " << a_vect_2[2] << endl;
  //cout << "     " << a_cross_x << "  " << a_cross_y << "  " << a_cross_z << endl;


}



void PhysMath::testLegendreFunctions(){
  TF1* Leg[20];
  Leg[0] = new TF1("L0","PhysMath::legendre(0,x)",-1,1);
  Leg[1] = new TF1("L1","PhysMath::legendre(1,x)",-1,1);
  Leg[2] = new TF1("L2","PhysMath::legendre(2,x)",-1,1);
  Leg[3] = new TF1("L3","PhysMath::legendre(3,x)",-1,1);
  Leg[4] = new TF1("L4","PhysMath::legendre(4,x)",-1,1);
  Leg[5] = new TF1("L5","PhysMath::legendre(5,x)",-1,1);
  Leg[6] = new TF1("L6","PhysMath::legendre(6,x)",-1,1);
  Leg[7] = new TF1("L7","PhysMath::legendre(7,x)",-1,1);
  Leg[8] = new TF1("L8","PhysMath::legendre(8,x)",-1,1);
  Leg[9] = new TF1("L9","PhysMath::legendre(9,x)",-1,1);
  Leg[10] = new TF1("L10","PhysMath::legendre(10,x)",-1,1);
  Leg[11] = new TF1("L11","PhysMath::legendre(11,x)",-1,1);
  Leg[12] = new TF1("L12","PhysMath::legendre(12,x)",-1,1);
  Leg[13] = new TF1("L13","PhysMath::legendre(13,x)",-1,1);
  Leg[14] = new TF1("L14","PhysMath::legendre(14,x)",-1,1);
  Leg[15] = new TF1("L15","PhysMath::legendre(15,x)",-1,1);
  Leg[16] = new TF1("L16","PhysMath::legendre(16,x)",-1,1);
  Leg[17] = new TF1("L17","PhysMath::legendre(17,x)",-1,1);
  Leg[18] = new TF1("L18","PhysMath::legendre(18,x)",-1,1);
  Leg[19] = new TF1("L19","PhysMath::legendre(19,x)",-1,1);

  TCanvas* canv = new TCanvas();

  for(int iii = 0; iii < 20; iii++){
    Leg[iii]->SetLineColor(kBlue);
    Leg[iii]->SetNpx(200);
  }
  TH1I* frame = new TH1I("asdf","",2,-1,1);
  frame->SetBinContent(1,-999);
  frame->SetBinContent(2,-999);

  frame->Draw();
  frame->GetYaxis()->SetRangeUser(-1.0,1.0);
  for(int iii = 0; iii < 20; iii++){
    Leg[iii]->SetLineColor(kBlue);
    Leg[iii]->SetNpx(200);
    Leg[iii]->Draw("same");
  }
}





void PhysMath::testTsallisNormalization(){
  TF1* funct = PhysMath::getTsallisFunction_NotInvariant("tsallis",1.0,0.35,0.05,0.15);
  double integral = funct->Integral(0.0,10.0);
  cout << " Tsallis dNdy Initialized to " << funct->GetParameter(0) << endl;
  cout << " Does this match the numeric integration =? " << integral << endl;
}

void PhysMath::testGausLikeFunctions(){
  TF1* student = getStudentTFunction("StudentT", 1.0, 0.0, 1.0, 3.5);
  TF1* gaus = getGaussianFunction("gaus", 1.0, 0.0, 1.0);
  TF1* doubleGaus = getDoubleGaussianFunction("doubleGaus",1.0, 0.0, 1.0, 0.2, 1.2);
  TF1* guasStudent = getGaussAndStudentTFunction("gausStudentT",1.0, 0.0, 1.0, 4.5, 0.2, 1.2);
  TF1* skew = skewNormalFunct( 1.0, 0.0, 1.0,5.0);
  TF1* genNorm = generalNormalV2Funct(1.0, 0.0, 1.0,0.5);
  TF1* pearson = pearsonIVFunct(1.0,0,3.0, 10.0,1.1,true); // lambda ~ mode     alpha ~ stdDev      m ~ kurtosis       nu ~ skewness

  student->SetLineColor(kRed);
  student->SetRange(-5,5);
  student->SetNpx(300);
  gaus->SetLineColor(kBlue);
  gaus->SetRange(-5,5);
  gaus->SetNpx(300);
  doubleGaus->SetLineColor(kGreen+1);
  doubleGaus->SetRange(-5,5);
  doubleGaus->SetNpx(300);
  guasStudent->SetLineColor(kOrange+1);
  guasStudent->SetRange(-5,5);
  guasStudent->SetNpx(300);
  skew->SetLineColor(kBlack);
  skew->SetRange(-5,5);
  skew->SetNpx(300);
  genNorm->SetLineColor(kViolet);
  genNorm->SetRange(-5,5);
  genNorm->SetNpx(300);
  pearson->SetLineColor(kYellow+1);
  pearson->SetRange(-5,5);
  pearson->SetNpx(300);


  TCanvas* canv = new TCanvas();
  student->Draw();
  gaus->Draw("same");
  doubleGaus->Draw("same");
  guasStudent->Draw("same");
  skew->Draw("same");
  genNorm->Draw("same");
  pearson->Draw("same");

}

void PhysMath::testCrossProduct(){

  TRandom3 randomThing;
  randomThing.SetSeed(time(NULL));

  double vec1[3];
  double vec2[3];
  double r1 = randomThing.Uniform(0,1);
  double r2 = randomThing.Uniform(0,1);
  randomThing.Sphere(vec1[0],vec1[1],vec1[2],r1);
  randomThing.Sphere(vec2[0],vec2[1],vec2[2],r2);
  cout << "Vect1: " << vec1[0] << "  " << vec1[1] << "  " << vec1[2] << endl;
  cout << "Vect2: " << vec2[0] << "  " << vec2[1] << "  " << vec2[2] << endl;

  double basicCross = crossProduct3D(vec1,vec2);
  double cross[3];
  crossProduct3D(vec1,vec2,cross[0],cross[1],cross[2]);
  //cross[0] = (vec1[2]*vec2[3]) - (vec1[3]*vec2[2]);
  //cross[1] = (vec1[3]*vec2[1]) - (vec1[1]*vec2[3]);
  //cross[2] = (vec1[1]*vec2[2]) - (vec1[2]*vec2[1]);

  double crossSize =  sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

  cout << "CROSS: " << cross[0] << "  " << cross[1] << "  " << cross[2] << "     mag = " << crossSize << endl;

  double dotProduct = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
  double theta   = acos(dotProduct/(r1*r2));
  double theta_2 = asin(crossSize/(r1*r2));
  double theta_3 = asin(basicCross/(r1*r2));

  cout << "Theta: " << theta << " =? " << theta_2 << " =? " << theta_3 <<  endl;

}



#endif
