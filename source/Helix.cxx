#ifndef _HELIX_CXX_
#define _HELIX_CXX_

#include "../headers/Helix.h"


//ClassImpT(Helix,double);

const Double_t Helix::NoSolution = 3.e+33;

//_________________
Helix::Helix() : mSingularity(false), mOrigin(0, 0, 0),
			     mDipAngle(0), mCurvature(0), mPhase(0),
			     mH(0), mCosDipAngle(0), mSinDipAngle(0),
			     mCosPhase(0), mSinPhase(0) {
  /*no-op*/
}

//_________________
Helix::Helix(Double_t c, Double_t d, Double_t phase,
			 const TVector3& o, Int_t h) {
  setParameters(c, d, phase, o, h);
}

Helix::Helix(const TVector3& p,
					 const TVector3& o,
					 Double_t B, Double_t q) {
  mH = (q*B <= 0) ? 1 : -1;
  if(p.y() == 0 && p.x() == 0) {
    setPhase((M_PI/4)*(1-2.*mH));
  }
  else {
    setPhase(atan2(p.y(),p.x())-mH*M_PI/2);
  }
  setDipAngle(atan2(p.z(),p.Perp()));
  mOrigin = o;

//#ifndef ST_NO_NAMESPACES
//  {
//    using namespace units;
//#endif
    setCurvature( ::fabs( (c_light*nanosecond/meter*q*B/tesla) /
			  ( p.Mag()/GeV*mCosDipAngle) / meter) );
//#ifndef ST_NO_NAMESPACES
//  }
//#endif
}

//_________________
Helix::Helix(const Helix &h) {
  mSingularity = h.mSingularity;
  mOrigin = h.mOrigin;
  mDipAngle = h.mDipAngle;
  mCurvature = h.mCurvature;
  mPhase = h.mPhase;
  mH = h.mH;
  mCosDipAngle = h.mCosDipAngle;
  mSinDipAngle = h.mSinDipAngle;
  mCosPhase = h.mCosPhase;
  mSinPhase = h.mSinPhase;
}

//_________________
Helix::~Helix() { /* noop */ };


Int_t Helix::h() const {return mH;}
Double_t Helix::dipAngle() const {return mDipAngle;}
Double_t Helix::curvature() const {return mCurvature;}
Double_t Helix::phase() const {return mPhase;}

Double_t Helix::x(Double_t s) const {
  if (mSingularity)
    return mOrigin.x() - s*mCosDipAngle*mSinPhase;
  else
    return mOrigin.x() + (cos(mPhase + s*mH*mCurvature*mCosDipAngle)-mCosPhase)/mCurvature;
}

Double_t Helix::y(Double_t s) const {
  if (mSingularity)
    return mOrigin.y() + s*mCosDipAngle*mCosPhase;
  else
    return mOrigin.y() + (sin(mPhase + s*mH*mCurvature*mCosDipAngle)-mSinPhase)/mCurvature;
}

Double_t Helix::z(Double_t s) const {
  return mOrigin.z() + s*mSinDipAngle;
}

Double_t Helix::cx(Double_t s) const {
  if (mSingularity)
    return -mCosDipAngle*mSinPhase;
  else
    return -sin(mPhase + s*mH*mCurvature*mCosDipAngle)*mH*mCosDipAngle;
}

Double_t Helix::cy(Double_t s) const {
  if (mSingularity)
    return mCosDipAngle*mCosPhase;
  else
    return cos(mPhase + s*mH*mCurvature*mCosDipAngle)*mH*mCosDipAngle;
}

Double_t Helix::cz(Double_t /* s */)  const { return mSinDipAngle; }
const TVector3& Helix::origin() const { return mOrigin; }
TVector3 Helix::at(Double_t s) const { return TVector3(x(s), y(s), z(s)); }
TVector3 Helix::cat(Double_t s) const { return TVector3(cx(s), cy(s), cz(s)); }
Double_t Helix::pathLength(Double_t X, Double_t Y) const { return fudgePathLength(TVector3(X, Y, 0)); }

Int_t Helix::bad(Double_t WorldSize) const {

  Int_t ierr = 0;
  if ( !std::isfinite(mDipAngle) ) {
    return 11;
  }
  if ( !std::isfinite(mCurvature) ) {
    return 12;
  }

  //ierr = mOrigin.bad(WorldSize);

  // The line above is commented and the StThreeVector::bad(double)
  // is rewritten here
  for(Int_t iIter=0; iIter<3; iIter++) {

    Double_t tmpVal;
    // Value StThreeVector.mX1[iter] ???
    switch(iIter) {
    case 0: tmpVal = mOrigin.X(); break;
    case 1: tmpVal = mOrigin.Y(); break;
    case 2: tmpVal = mOrigin.Z(); break;
    default: tmpVal = NAN;
    };

    if ( !std::isfinite( tmpVal ) ) {
      ierr = 10 + iIter;
    }
    if ( ::fabs( tmpVal ) > WorldSize ) {
      ierr = 20 + iIter;
    }
  } //for(Int_t iIter=0; iIter<3; iIter+)

  if (ierr) {
    return (3+ierr*100);
  }
  if ( ::fabs(mDipAngle) > 1.58 ) {
    return 21;
  }

  Double_t qwe = ::fabs( ::fabs(mDipAngle) - M_PI/2 );
  if ( qwe < 1./WorldSize ) {
    return 31;
  }

  if ( ::fabs(mCurvature) > WorldSize )	{
    return 22;
  }
  if ( mCurvature < 0 )	{
    return 32;
  }
  if (abs(mH) != 1 ) {
    return 24;
  }

  return 0;
}









//_________________
void Helix::setParameters(Double_t c, Double_t dip, Double_t phase,
				const TVector3& o, Int_t h) {

  //
  //  The order in which the parameters are set is important
  //  since setCurvature might have to adjust the others.
  //
  mH = (h>=0) ? 1 : -1;    // Default is: positive particle
                           //             positive field
  mOrigin   = o;
  setDipAngle(dip);
  setPhase(phase);

  //
  // Check for singularity and correct for negative curvature.
  // May change mH and mPhase. Must therefore be set last.
  //
  setCurvature(c);

  //
  // For the case B=0, h is ill defined. In the following we
  // always assume h = +1. Since phase = psi - h * pi/2
  // we have to correct the phase in case h = -1.
  // This assumes that the user uses the same h for phase
  // as the one he passed to the constructor.
  //
  if (mSingularity && mH == -1) {
    mH = +1;
    setPhase(mPhase-M_PI);
  }
}

//_________________
void Helix::setCurvature(Double_t val) {

  if (val < 0) {
    mCurvature = -val;
    mH = -mH;
    setPhase( mPhase+M_PI );
  }
  else {
    mCurvature = val;
  }

#ifndef ST_NO_NUMERIC_LIMITS
  if ( ::fabs(mCurvature) <= std::numeric_limits<Double_t>::epsilon() ) {
#else
  if ( ::fabs(mCurvature) <= static_cast<Double_t>(0) ) {
#endif
    mSingularity = true;  // straight line
  }
  else {
    mSingularity = false; // curved
  }
}

//_________________
void Helix::setPhase(Double_t val) {

  mPhase       = val;
  mCosPhase    = cos(mPhase);
  mSinPhase    = sin(mPhase);
  if ( ::fabs(mPhase) > M_PI) {
    mPhase = atan2(mSinPhase, mCosPhase);  // force range [-pi,pi]
  }
}

//_________________
void Helix::setDipAngle(Double_t val) {

  mDipAngle    = val;
  mCosDipAngle = cos(mDipAngle);
  mSinDipAngle = sin(mDipAngle);
}

//_________________
Double_t Helix::xcenter() const {

  if (mSingularity) {
    return 0;
  }
  else {
    return mOrigin.x()-mCosPhase/mCurvature;
  }
}

//_________________
Double_t Helix::ycenter() const {
  if (mSingularity) {
    return 0;
  }
  else {
    return mOrigin.y()-mSinPhase/mCurvature;
  }
}

//_________________
Double_t Helix::fudgePathLength(const TVector3& p) const {

  Double_t s;
  Double_t dx = p.x()-mOrigin.x();
  Double_t dy = p.y()-mOrigin.y();

  if (mSingularity) {
    s = (dy*mCosPhase - dx*mSinPhase)/mCosDipAngle;
  }
  else {
    s = atan2(dy*mCosPhase - dx*mSinPhase,
	      1/mCurvature + dx*mCosPhase+dy*mSinPhase)/
      (mH*mCurvature*mCosDipAngle);
  }
  return s;
}

//_________________
Double_t Helix::distance(const TVector3& p, Bool_t scanPeriods) const {
  return ( this->at( pathLength(p,scanPeriods) )-p ).Mag();
}

//_________________
Double_t Helix::pathLength(const TVector3& p, Bool_t scanPeriods) const {

  //
  //  Returns the path length at the distance of closest
  //  approach between the Helix and point p.
  //  For the case of B=0 (straight line) the path length
  //  can be calculated analytically. For B>0 there is
  //  unfortunately no easy solution to the problem.
  //  Here we use the Newton method to find the root of the
  //  referring equation. The 'fudgePathLength' serves
  //  as a starting value.
  //

  Double_t s;
  Double_t dx = p.x()-mOrigin.x();
  Double_t dy = p.y()-mOrigin.y();
  Double_t dz = p.z()-mOrigin.z();

  if (mSingularity) {
    s = mCosDipAngle*(mCosPhase*dy-mSinPhase*dx) +
      mSinDipAngle*dz;
  }
  else {
//#ifndef ST_NO_NUMERIC_LIMITS
//    {
//      using namespace units;
//#endif
  const Double_t MaxPrecisionNeeded = micrometer;
  const Int_t    MaxIterations      = 100;

  //
  // The math is taken from Maple with C(expr,optimized) and
  // some hand-editing. It is not very nice but efficient.
  //
  Double_t t34 = mCurvature*mCosDipAngle*mCosDipAngle;
  Double_t t41 = mSinDipAngle*mSinDipAngle;
  Double_t t6, t7, t11, t12, t19;

  //
  // Get a first guess by using the dca in 2D. Since
  // in some extreme cases we might be off by n periods
  // we add (subtract) periods in case we get any closer.
  //
  s = fudgePathLength(p);

  if (scanPeriods) {
	  Double_t ds = period();
	  Int_t    j;
	  Int_t jmin = 0;
	  Double_t d;
	  Double_t dmin = ( at(s) - p).Mag() ;

	  for(j=1; j<MaxIterations; j++) {
			d = ( at(s+j*ds) - p ).Mag();
			if ( d  < dmin) {
			dmin = d;
			jmin = j;
			}
			else {
			break;
			}
	  } //for(j=1; j<MaxIterations; j++)

	  for(j=-1; -j<MaxIterations; j--) {
			d = ( at(s+j*ds ) - p ).Mag() ;
			if ( d < dmin) {
			dmin = d;
			jmin = j;
			}
			else {
			break;
			}
	  } //for(j=-1; -j<MaxIterations; j--)

	  if (jmin) {
			s += jmin*ds;
	  }
  } //if (scanPeriods)

    //
    // Newtons method:
    // Stops after MaxIterations iterations or if the required
    // precision is obtained. Whatever comes first.
    //
    Double_t sOld = s;
    for (Int_t i=0; i<MaxIterations; i++) {
      t6  = mPhase+s*mH*mCurvature*mCosDipAngle;
      t7  = cos(t6);
      t11 = dx-(1/mCurvature)*(t7-mCosPhase);
      t12 = sin(t6);
      t19 = dy-(1/mCurvature)*(t12-mSinPhase);
      s  -= (t11*t12*mH*mCosDipAngle-t19*t7*mH*mCosDipAngle -
        (dz-s*mSinDipAngle)*mSinDipAngle)/
      (t12*t12*mCosDipAngle*mCosDipAngle+t11*t7*t34 +
        t7*t7*mCosDipAngle*mCosDipAngle +
        t19*t12*t34+t41);
      if (fabs(sOld-s) < MaxPrecisionNeeded) break;
      sOld = s;
    } //for (Int_t i=0; i<MaxIterations; i++)
//#ifndef ST_NO_NUMERIC_LIMITS
//    }
//#endif
  } //else
  //std::cout << "Path length to position (" << p.x() << ", " << p.y() << ", " << p.z() << "): " << s << std::endl;

  return s;
}

//_________________
Double_t Helix::period() const {
  if (mSingularity) {
//#ifndef ST_NO_NUMERIC_LIMITS
    return std::numeric_limits<Double_t>::max();
//#else
//    return DBL_MAX;
//#endif
  }
  else {
    return fabs(2*M_PI/(mH*mCurvature*mCosDipAngle));
  }
}

//_________________
std::pair<Double_t, Double_t> Helix::pathLength(Double_t r) const {

  std::pair<Double_t,Double_t> value;
  std::pair<Double_t,Double_t> VALUE(999999999.,999999999.);
  //
  // The math is taken from Maple with C(expr,optimized) and
  // some hand-editing. It is not very nice but efficient.
  // 'first' is the smallest of the two solutions (may be negative)
  // 'second' is the other.
  //
  if (mSingularity) {
    Double_t t1 = mCosDipAngle*(mOrigin.x()*mSinPhase-mOrigin.y()*mCosPhase);
    Double_t t12 = mOrigin.y()*mOrigin.y();
    Double_t t13 = mCosPhase*mCosPhase;
    Double_t t15 = r*r;
    Double_t t16 = mOrigin.x()*mOrigin.x();
    Double_t t20 = -mCosDipAngle*mCosDipAngle*(2.0*mOrigin.x()*mSinPhase*mOrigin.y()*mCosPhase +
					       t12-t12*t13-t15+t13*t16);
    if (t20<0.) {
      return VALUE;
    }
    t20 = ::sqrt(t20);
    value.first  = (t1-t20)/(mCosDipAngle*mCosDipAngle);
    value.second = (t1+t20)/(mCosDipAngle*mCosDipAngle);
  }
  else {
    Double_t t1 = mOrigin.y()*mCurvature;
    Double_t t2 = mSinPhase;
    Double_t t3 = mCurvature*mCurvature;
    Double_t t4 = mOrigin.y()*t2;
    Double_t t5 = mCosPhase;
    Double_t t6 = mOrigin.x()*t5;
    Double_t t8 = mOrigin.x()*mOrigin.x();
    Double_t t11 = mOrigin.y()*mOrigin.y();
    Double_t t14 = r*r;
    Double_t t15 = t14*mCurvature;
    Double_t t17 = t8*t8;
    Double_t t19 = t11*t11;
    Double_t t21 = t11*t3;
    Double_t t23 = t5*t5;
    Double_t t32 = t14*t14;
    Double_t t35 = t14*t3;
    Double_t t38 = 8.0*t4*t6 - 4.0*t1*t2*t8 - 4.0*t11*mCurvature*t6 +
      4.0*t15*t6 + t17*t3 + t19*t3 + 2.0*t21*t8 + 4.0*t8*t23 -
      4.0*t8*mOrigin.x()*mCurvature*t5 - 4.0*t11*t23 -
      4.0*t11*mOrigin.y()*mCurvature*t2 + 4.0*t11 - 4.0*t14 +
      t32*t3 + 4.0*t15*t4 - 2.0*t35*t11 - 2.0*t35*t8;
    Double_t t40 = (-t3*t38);
    if (t40<0.) {
      return VALUE;
    }
    t40 = ::sqrt(t40);

    Double_t t43 = mOrigin.x()*mCurvature;
    Double_t t45 = 2.0*t5 - t35 + t21 + 2.0 - 2.0*t1*t2 -2.0*t43 - 2.0*t43*t5 + t8*t3;
    Double_t t46 = mH*mCosDipAngle*mCurvature;

    value.first = (-mPhase + 2.0*atan((-2.0*t1 + 2.0*t2 + t40)/t45))/t46;
    value.second = -(mPhase + 2.0*atan((2.0*t1 - 2.0*t2 + t40)/t45))/t46;

    //
    //   Solution can be off by +/- one period, select smallest
    //
    Double_t p = period();

    if ( ! std::isnan(value.first) ) {
      if ( ::fabs(value.first-p) < ::fabs(value.first) ) {
	value.first = value.first-p;
      }
      else if ( ::fabs(value.first+p) < ::fabs(value.first) ) {
	value.first = value.first+p;
      }
    } //if ( ! std::isnan(value.first) )

    if (! std::isnan(value.second)) {
      if ( ::fabs(value.second-p) < ::fabs(value.second) ) {
	value.second = value.second-p;
      }
      else if ( ::fabs(value.second+p) < ::fabs(value.second) ) {
	value.second = value.second+p;
      }
    } //if (! std::isnan(value.second))
  } //else

  if (value.first > value.second) {
    std::swap(value.first,value.second);
  }

  return(value);
}

//_________________
std::pair<Double_t, Double_t> Helix::pathLength(Double_t r, Double_t x, Double_t y) {
  Double_t x0 = mOrigin.x();
  Double_t y0 = mOrigin.y();
  mOrigin.SetX(x0-x);
  mOrigin.SetY(y0-y);
  std::pair<Double_t, Double_t> result = this->pathLength(r);
  mOrigin.SetX(x0);
  mOrigin.SetY(y0);
  return result;
}

//________________
Double_t Helix::pathLength(const TVector3& r,
				 const TVector3& n) const {
  //
  // Vector 'r' defines the position of the center and
  // vector 'n' the normal vector of the plane.
  // For a straight line there is a simple analytical
  // solution. For curvatures > 0 the root is determined
  // by Newton method. In case no valid s can be found
  // the max. largest value for s is returned.
  //
  Double_t s;

  if (mSingularity) {
    Double_t t = n.z()*mSinDipAngle +
      n.y()*mCosDipAngle*mCosPhase -
      n.x()*mCosDipAngle*mSinPhase;
    if (t == 0) {
      s = NoSolution;
    }
    else {
      s = ((r - mOrigin)*n)/t;
    }
  }
  else {
    const Double_t MaxPrecisionNeeded = micrometer;
    const Int_t    MaxIterations      = 20;

    Double_t A = mCurvature*((mOrigin - r)*n) -
      n.x()*mCosPhase -
      n.y()*mSinPhase;
    Double_t t = mH*mCurvature*mCosDipAngle;
    Double_t u = n.z()*mCurvature*mSinDipAngle;

    Double_t a, f, fp;
    Double_t sOld = s = 0;
    //Double_t shiftOld = 0;
    Double_t shift;
    //  (cos(angMax)-1)/angMax = 0.1
    const Double_t angMax = 0.21;
    Double_t deltas = fabs(angMax/(mCurvature*mCosDipAngle));
    //  dampingFactor = exp(-0.5);
    //	Double_t dampingFactor = 0.60653;
    Int_t i;

    for (i=0; i<MaxIterations; i++) {
      a  = t*s+mPhase;
      Double_t sina = sin(a);
      Double_t cosa = cos(a);
      f = A + n.x()*cosa + n.y()*sina +	u*s;
      fp = -n.x()*sina*t + n.y()*cosa*t + u;

      if ( fabs(fp)*deltas <= fabs(f) ) { //too big step
        Int_t sgn = 1;
        if (fp<0.) sgn = -sgn;
        if (f <0.) sgn = -sgn;
        shift = sgn*deltas;
        if (shift<0) shift*=0.9;  // don't get stuck shifting +/-deltas
      }
      else {
        shift = f/fp;
      }
      s -= shift;
      //shiftOld = shift;
      if ( ::fabs(sOld-s) < MaxPrecisionNeeded ) {
        break;
      }
      sOld = s;
    } //for (i=0; i<MaxIterations; i++)

    if (i == MaxIterations) {
      return NoSolution;
    }
  } //else
  return s;
}

//_________________
std::pair<Double_t, Double_t> Helix::pathLengths(const Helix& h, Double_t minStepSize,
						  Double_t minRange) const {
  //
  //	Cannot handle case where one is a Helix
  //  and the other one is a straight line.
  //
  if (mSingularity != h.mSingularity) {
    return std::pair<Double_t, Double_t>(NoSolution, NoSolution);
  }

  Double_t s1, s2;

  if (mSingularity) {
    //
    //  Analytic solution
    //
    TVector3 dv = h.mOrigin - mOrigin;
    TVector3 a(-mCosDipAngle*mSinPhase,
	       mCosDipAngle*mCosPhase,
	       mSinDipAngle);
    TVector3 b(-h.mCosDipAngle*h.mSinPhase,
	       h.mCosDipAngle*h.mCosPhase,
	       h.mSinDipAngle);
    Double_t ab = a*b;
    Double_t g  = dv*a;
    Double_t k  = dv*b;
    s2 = (k-ab*g)/(ab*ab-1.);
    s1 = g+s2*ab;
    return std::pair<Double_t, Double_t>(s1, s2);
  } //if (mSingularity)
  else {
    //
    //  First step: get dca in the xy-plane as start value
    //
    Double_t dx = h.xcenter() - xcenter();
    Double_t dy = h.ycenter() - ycenter();
    Double_t dd = ::sqrt(dx*dx + dy*dy);
    Double_t r1 = 1/curvature();
    Double_t r2 = 1/h.curvature();

    Double_t cosAlpha = (r1*r1 + dd*dd - r2*r2)/(2*r1*dd);

    Double_t s;
    Double_t x, y;
    if ( ::fabs(cosAlpha) < 1) {           // two solutions
      Double_t sinAlpha = sin(acos(cosAlpha));
      x = xcenter() + r1*(cosAlpha*dx - sinAlpha*dy)/dd;
      y = ycenter() + r1*(sinAlpha*dx + cosAlpha*dy)/dd;
      s = pathLength(x, y);
      x = xcenter() + r1*(cosAlpha*dx + sinAlpha*dy)/dd;
      y = ycenter() + r1*(cosAlpha*dy - sinAlpha*dx)/dd;
      Double_t a = pathLength(x, y);
      if ( h.distance(at(a)) < h.distance(at(s)) ) {
        s = a;
      }
    } //if ( ::fabs(cosAlpha) < 1)
    else {                              // no intersection (or exactly one)
      Int_t rsign = ((r2-r1) > dd ? -1 : 1); // set -1 when *this* Helix is
      // completely contained in the other
      x = xcenter() + rsign*r1*dx/dd;
      y = ycenter() + rsign*r1*dy/dd;
      s = pathLength(x, y);
    } //else

    //
    //   Second step: scan in decreasing intervals around seed 's'
    //   minRange and minStepSize are passed as arguments to the method.
    //   They have default values defined in the header file.
    //
    Double_t dmin              = h.distance(at(s));
    Double_t range             = std::max(2*dmin, minRange);
    Double_t ds                = range/10;
    Double_t slast=-999999, ss, d;
    s1 = s - range/2.;
    s2 = s + range/2.;

    while (ds > minStepSize) {

      for ( ss=s1; ss<s2+ds; ss+=ds ) {
        d = h.distance(at(ss));
        if (d < dmin) {
          dmin = d;
          s = ss;
        }
        slast = ss;
      } //for ( ss=s1; ss<s2+ds; ss+=ds )

      //
      //  In the rare cases where the minimum is at the
      //  the border of the current range we shift the range
      //  and start all over, i.e we do not decrease 'ds'.
      //  Else we decrease the search intervall around the
      //  current minimum and redo the scan in smaller steps.
      //
      if (s == s1) {
        d = 0.8*(s2-s1);
        s1 -= d;
        s2 -= d;
      }
      else if (s == slast) {
        d = 0.8*(s2-s1);
        s1 += d;
        s2 += d;
      }
      else {
        s1 = s-ds;
        s2 = s+ds;
        ds /= 10;
      }
    } //while (ds > minStepSize)
    return std::pair<Double_t, Double_t>(s, h.pathLength(at(s)));
  } //else
}

//_________________
void Helix::moveOrigin(Double_t s) {

  if (mSingularity) {
    mOrigin = at(s);
  }
  else {
    TVector3 newOrigin = at(s);
    Double_t newPhase = atan2(newOrigin.y() - ycenter(),
			      newOrigin.x() - xcenter());
    mOrigin = newOrigin;
    setPhase(newPhase);
  }
}




//_________________
TVector3 Helix::momentum(Double_t B) const {

  if (mSingularity) {
    return(TVector3(0,0,0));
  }
  else {
//#ifndef ST_NO_NAMESPACES
//    {
//      using namespace units;
//#endif
      Double_t pt = GeV*fabs(c_light*nanosecond/meter*B/tesla)/(fabs(mCurvature)*meter);

      return ( TVector3( pt*cos(mPhase+mH*M_PI/2),   // pos part pos field
			 pt*sin(mPhase+mH*M_PI/2),
			 pt*tan(mDipAngle) ) );
//#ifndef ST_NO_NAMESPACES
//    }
//#endif
  } //else
}

//_________________
TVector3 Helix::momentumAt(Double_t S, Double_t B) const {
  // Obtain phase-shifted momentum from phase-shift of origin
  Helix tmp(*this);
  tmp.moveOrigin(S);
  return tmp.momentum(B);
}

//_________________
Int_t Helix::charge(Double_t B) const {
  return (B > 0 ? -mH : mH);
}

//_________________
Double_t Helix::geometricSignedDistance(Double_t x, Double_t y) {
  // Geometric signed distance
  Double_t thePath = this->pathLength(x,y);
  TVector3 DCA2dPosition = this->at(thePath);
  DCA2dPosition.SetZ(0);
  TVector3 position(x,y,0);
  TVector3 DCAVec = (DCA2dPosition-position);
  TVector3 momVec;
  // Deal with straight tracks
  if (this->mSingularity) {
    momVec = this->at(1)- this->at(0);
    momVec.SetZ(0);
  }
  else {
    momVec = this->momentumAt(thePath,1./tesla); // Don't care about Bmag.  Helicity is what matters.
    momVec.SetZ(0);
  }

  Double_t cross = DCAVec.x()*momVec.y() - DCAVec.y()*momVec.x();
  Double_t theSign = (cross>=0) ? 1. : -1.;
  return theSign*DCAVec.Perp();
}

//_________________
Double_t Helix::curvatureSignedDistance(Double_t x, Double_t y) {
  // Protect against mH = 0 or zero field
  if (this->mSingularity || abs(this->mH)<=0) {
    return (this->geometricSignedDistance(x,y));
  }
  else {
    return (this->geometricSignedDistance(x,y))/(this->mH);
  }
}

//_________________
Double_t Helix::geometricSignedDistance(const TVector3& pos) {
  Double_t sdca2d = this->geometricSignedDistance(pos.x(),pos.y());
  Double_t theSign = (sdca2d>=0) ? 1. : -1.;
  return (this->distance(pos))*theSign;
}

//_________________
Double_t Helix::curvatureSignedDistance(const TVector3& pos) {
  Double_t sdca2d = this->curvatureSignedDistance(pos.x(),pos.y());
  Double_t theSign = (sdca2d>=0) ? 1. : -1.;
  return (this->distance(pos))*theSign;
}


void Helix::dca(const TVector3& a_vertex, double &a_dca, double &a_dcaX, double &a_dcaY, double &a_dcaZ){
  TVector3 dca_vect = this->at( pathLength(a_vertex,true) )- a_vertex;
  a_dca = dca_vect.Mag();
  a_dcaX = dca_vect.x();
  a_dcaY = dca_vect.y();
  a_dcaZ = dca_vect.z();
}



TGraph2D* Helix::visualizeHelixTrajectory(double pathLengthStart, double pathLengthEnd, int numPoints) {
  TGraph2D* helixGraph = new TGraph2D();
  helixGraph->SetTitle("Helix Trajectory;X;Y;Z");
  helixGraph->SetMarkerStyle(20);
  helixGraph->SetMarkerSize(0.5);
  helixGraph->SetMarkerColor(kRed+1);

  double pathLengthStep = (pathLengthEnd - pathLengthStart) / numPoints;

  for (int i = 0; i <= numPoints; ++i) {
    double s = pathLengthStart + i * pathLengthStep;
    TVector3 position = this->at(s);
    helixGraph->SetPoint(i, position.x(), position.y(), position.z());
  }
  return helixGraph;
}






void testHelix(){

  // Test Helix
  TVector3 p(1,1,1);
  TVector3 o(0,0,0);
  Double_t B = 1;
  Double_t q = 1;
  Helix hp(p,o,B,q);
  TVector3 p2 = hp.momentum(B);
  std::cout << "Momentum: " << p2.x() << " " << p2.y() << " " << p2.z() << std::endl;
  std::cout << "Charge: " << hp.charge(B) << std::endl;
  std::cout << "Curvature: " << hp.curvature() << std::endl;
  std::cout << "Dip Angle: " << hp.dipAngle() << std::endl;
  std::cout << "Phase: " << hp.phase() << std::endl;
  std::cout << "Origin: " << hp.origin().x() << " " << hp.origin().y() << " " << hp.origin().z() << std::endl;
  std::cout << "Path Length: " << hp.pathLength(1,1) << std::endl;
  std::cout << "At: " << hp.at(1).x() << " " << hp.at(1).y() << " " << hp.at(1).z() << std::endl;
  std::cout << "Distance: " << hp.distance(TVector3(1,1,1)) << std::endl;
  std::cout << "Geometric Signed Distance: " << hp.geometricSignedDistance(1,1) << std::endl;
  std::cout << "Curvature Signed Distance: " << hp.curvatureSignedDistance(1,1) << std::endl;
  std::cout << "Geometric Signed Distance: " << hp.geometricSignedDistance(TVector3(1,1,1)) << std::endl;
  std::cout << "Curvature Signed Distance: " << hp.curvatureSignedDistance(TVector3(1,1,1)) << std::endl;
  std::cout << "Momentum At: " << hp.momentumAt(1,B).x() << " " << hp.momentumAt(1,B).y() << " " << hp.momentumAt(1,B).z() << std::endl;

  std::cout << "\n\n\n\n\n" << std::endl;

  std::cout << "Testing Straight Track:" << std::endl;
  TVector3 p_straight(100,100,0);
  TVector3 o_straight(1,1,0);
  Helix hp_straight(p_straight,o_straight,B,q);
  TVector3 p2_straight = hp_straight.momentum(B);
  std::cout << "Momentum: " << p_straight.x() << " " << p_straight.y() << " " << p_straight.z() << std::endl;
  std::cout << "Momentum: " << p2_straight.x() << " " << p2_straight.y() << " " << p2_straight.z() << std::endl;
  std::cout << "DCA from (1,1,0) to (0,0,0): " << hp_straight.geometricSignedDistance(TVector3(0,0,0)) << std::endl;


  TGraph* dca_graph = new TGraph();

  for(int index = 0; index < 1000; index++){
    double z_pos = -50 + ((double)index)*0.1;
    double phase = ((double)index)*2.0*TMath::Pi()/1000.0;
    TVector3 p_straight(0.3,0.3,0.3);
    TVector3 o_straight(1.0+cos(phase),1.0+sin(phase),0);
    Helix hp_straight(p_straight,o_straight,B,q);
    double dca =  hp_straight.geometricSignedDistance(TVector3(0,0,0));
    dca_graph->SetPoint(index,phase,dca);
  }
  TH1D* frame = new TH1D("frame","DCA from Changing Phase; Phase; DCA",1000,0,2*TMath::Pi());
  frame->SetMaximum(3);
  frame->SetMinimum(0.0);
  frame->SetStats(0);
  frame->Draw();
  dca_graph->Draw("LP");



  // TEST HELIX PATH WITH MATH

  double mass = 1.67262192e-27;
  double light = 2.99792458e8;
  double GeV = 1.60217663e-10;

  /*DCA Debug:    B:-5e+09   q:1
  Mom: (-0.280106, -0.879413, -0.0406705)
  Decay: (-0.0462766, -0.1185, 0.450629)
  Event: (0.00086155, -0.00081602, 0.446654)
  DCA: 0.126835
  DCAs: 1.32511  :  12.6835   12.6835   12.6835   12.6835
  */

  TVector3 momentum_inGeV(0.535903, 0.204477, -0.168263);
  TVector3 decay_vertex(-0.0462766, -0.1185, 0.450629);
  TVector3 event_vertex(0.00086155, -0.00081602, 0.446654);
  //TVector3 momentum(0.5*GeV/light,0.5*GeV/light,0.2*GeV/light);
  TVector3 momentum(momentum_inGeV.x()*GeV/light,momentum_inGeV.y()*GeV/light,momentum_inGeV.z()*GeV/light);
  double pTot = momentum.Mag();
  double betaTot = light*pTot/sqrt(light*light*pTot*pTot+mass*mass*pow(light,4));
  double gamma = 1.0/sqrt(1.0-betaTot*betaTot);
  TVector3 betaVect(momentum.x()*betaTot/pTot,momentum.y()*betaTot/pTot,momentum.z()*betaTot/pTot);
  double beta_T = sqrt(betaVect.x()*betaVect.x()+betaVect.y()*betaVect.y());
  Helix helix(momentum_inGeV,decay_vertex,0.5*1E-11,1);

  double omega = 1.60217663e-19*0.5/(gamma*mass);
  //double gm_qb = gamma*mass/(1.60217663e-19*0.5);
  double radius = gamma*mass*light*beta_T/(1.60217663e-19*0.5);
  double phase = atan2(betaVect.x(),betaVect.y());

  double c1 = decay_vertex.x() + radius*cos(phase);
  double c2 = decay_vertex.y() - radius*sin(phase);

  std::cout << "momentum: " << momentum.x() << " " << momentum.y() << " " << momentum.z() << std::endl;
  std::cout << "betaVect: " << betaVect.x() << " " << betaVect.y() << " " << betaVect.z() << std::endl;
  std::cout << "decay_vertex: " << decay_vertex.x() << " " << decay_vertex.y() << " " << decay_vertex.z() << std::endl;
  std::cout << "pTot: " << pTot << std::endl;
  std::cout << "Gamma: " << gamma << std::endl;
  std::cout << "Beta_tot: " << betaTot << std::endl;
  std::cout << "Beta_T: " << beta_T << std::endl;
  std::cout << "Omega: " << omega << std::endl;
  std::cout << "radius: " << radius << std::endl;
  std::cout << "Phase: " << phase << std::endl;
  std::cout << "C1: " << c1 << std::endl;
  std::cout << "C2: " << c2 << std::endl;

  TGraph* position_xy_graph = new TGraph();
  TGraph* position_xz_graph = new TGraph();
  TGraph* position_yz_graph = new TGraph();
  TGraph2D* position_xyz_graph = new TGraph2D();
  position_xy_graph->SetMarkerStyle(20);
  position_xz_graph->SetMarkerStyle(20);
  position_yz_graph->SetMarkerStyle(20);
  position_xyz_graph->SetMarkerStyle(20);
  position_xy_graph->SetMarkerColor(kGreen + 2);
  position_xz_graph->SetMarkerColor(kBlue + 1);
  position_yz_graph->SetMarkerColor(kRed + 1);

  int index = 0;
  double time_scale = 1E-6;
  double dca_check = 1E200;
  TVector3 dca_point(0,0,0);

  double dist_to_origin = 1E99;
  TVector3 closest_point(0,0,0);
  for(double time = -time_scale; time < time_scale; time+=1e-5*time_scale, index++){
    double x = c1 - radius*cos(omega*time + phase);
    double y = c2 + radius*sin(omega*time + phase);
    double z = decay_vertex.z() + light*betaVect.z()*time;
    position_xy_graph->SetPoint(index,x,y);
    position_xz_graph->SetPoint(index,x,z);
    position_yz_graph->SetPoint(index,y,z);
    position_xyz_graph->SetPoint(index,x,y,z);
    std::cout << "Time: " << time << " X: " << x << " Y: " << y << " Z: " << z << std::endl;

    //double dca = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    double dca = sqrt(pow(x-event_vertex.x(),2)+pow(y-event_vertex.y(),2)+pow(z-event_vertex.z(),2));
    if(dca < dca_check){
      dca_check = dca;
      dca_point.SetXYZ(x,y,z);
    }
    double dist = sqrt(pow(x-decay_vertex.x(),2)+pow(y-decay_vertex.y(),2)+pow(z-decay_vertex.z(),2));
    if(dist < dist_to_origin){
      dist_to_origin = dist;
      closest_point.SetXYZ(x,y,z);
    }
  }

  TGraph2D* helix_graph = helix.visualizeHelixTrajectory(-100, 100, 100000);

  std::cout << "DCA: " << dca_check << " Point: " << dca_point.x() << " " << dca_point.y() << " " << dca_point.z() << std::endl;
  std::cout << "Dist to Origin: " << dist_to_origin << " Point: " << closest_point.x() << " " << closest_point.y() << " " << closest_point.z() << std::endl;
  std::cout << "DCA Helix: " << helix.geometricSignedDistance(event_vertex) << std::endl;

  double scale = 15;
  TCanvas* canv2 = new TCanvas("canv2","canv2",800,800);
  TH1D* frame2 = new TH1D("frame2","Helix Paths",100,-scale,scale);
  frame2->SetMaximum(scale);
  frame2->SetMinimum(-scale);
  frame2->SetStats(0);
  frame2->Draw();
  position_xy_graph->Draw("P");
  position_xz_graph->Draw("P");
  position_yz_graph->Draw("P");

  TCanvas* canv3 = new TCanvas("canv3","canv3",800,800);
  position_xyz_graph->Draw("P0");
  helix_graph->Draw("PLSAME");





}


#endif
