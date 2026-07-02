#ifndef _HELIX_
#define _HELIX_

// Adaptation of StPicoHelix class from STAR software package
// C++ headers
#include <cmath>
#include <utility>
#include <algorithm>
#include <limits>
#include <iostream>
#include <float.h>
#include <ostream>

// ROOT headers
#include "TVector3.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"

/*
#ifdef _VANILLA_ROOT_
  #include "SystemOfUnits.h"
#else
  #include "StarClassLibrary/SystemOfUnits.h"
#endif
*/

// Declare C++ namespaces
//#if !defined(ST_NO_NAMESPACES)
//  using namespace std;
//  using std::pair;
//  using std::swap;
//  using std::max;
//  using std::numeric_limits;
//#endif


// Originally an inline "system of units copy" duplicating centimeter/micrometer/
// GeV/etc as global static const doubles, guarded by "#ifndef HEP_SYSTEM_OF_UNITS_H"
// (that guard is never #define'd here, so it was only ever a no-op protection
// against a SystemOfUnits.h that happened to be included earlier in the SAME
// translation unit -- in the original repo, Helix.h/.cxx were compiled inside
// namespaces.cxx as their own ACLiC unit, separate from PicoBinner.cxx's unit
// (which pulls in the PicoDstReader submodule's SystemOfUnits.h via
// StPicoHelix.h), so the two copies never collided there.
//
// In this repo, ROOT's interpreter (cling) carries symbols across ACLiC'd units
// compiled in the same session, so once namespaces.cxx's copy of "micrometer" (etc)
// was interpreted, PicoBinner.cxx's later pull-in of the submodule's
// "using namespace units;" (in submodule/PicoDstReader_SL23a/SystemOfUnits.h) made
// "micrometer" ambiguous. Since the values are identical either way, this now just
// includes that same SystemOfUnits.h (guarded by its own HEP_SYSTEM_OF_UNITS_H,
// so it's a no-op if already included) instead of redefining the constants --
// same numbers, one definition.
#include "../submodule/PicoDstReader_SL23a/SystemOfUnits.h"

// c_light isn't defined in SystemOfUnits.h's units namespace, but Helix.cxx uses it
// (setCurvature() and momentum()), so it stays here on its own -- no collision risk
// since nothing else defines "c_light".
static const double c_light = 2.99792458e+8 * meter/second;


//_________________
class Helix {

 public:
  /// Default constructor
  Helix();

  /// Constructor that takes next arguments:
  /// curvature, dip angle, phase, origin, h
  Helix(Double_t c, Double_t dip, Double_t phase,
	      const TVector3& o, Int_t h=-1);

  /// Copy constructor
  Helix(const Helix&);

  /// Constructor with momentum (in GeV/c), origin (in meters), signed Magnetic Field (-0.5E-11 for ReverseFullField)
  /// and Charge of particle (+/- 1)
  Helix(const TVector3&,
		      const TVector3&,
		      Double_t, Double_t);




  // Assignment operator (will use the one, provided by compiler)
  //Helix& operator=(const Helix&);

  /// Destructor
  ~Helix();

  /// Return dip angle
  Double_t dipAngle()   const;
  /// Return curvature: 1/R in xy-plane
  Double_t curvature()  const;
  /// Return phase: aziumth in xy-plane measured from ring center
  Double_t phase()      const;
  /// Return x-center of circle in xy-plane
  Double_t xcenter()    const;
  /// Return y-center of circle in xy-plane
  Double_t ycenter()    const;
  /// Return -sign(q*B);
  Int_t    h()          const;

  /// Return origin of the Helix = starting point
  const TVector3& origin() const;

  /// Set Helix parameters
  void setParameters(Double_t c, Double_t dip, Double_t phase, const TVector3& o, Int_t h);

  /// coordinates of Helix at point s
  Double_t x(Double_t s)  const;
  Double_t y(Double_t s)  const;
  Double_t z(Double_t s)  const;
  TVector3 at(Double_t s) const;

  /// pointing vector of Helix at point s
  Double_t cx(Double_t s)  const;
  Double_t cy(Double_t s)  const;
  Double_t cz(Double_t s = 0)  const;
  TVector3 cat(Double_t s) const;

  /// returns period length of Helix
  Double_t period()       const;

  /// path length at given r (cylindrical r)
  std::pair<Double_t, Double_t> pathLength(Double_t r)   const;

  /// path length at given r (cylindrical r, cylinder axis at x,y)
  std::pair<Double_t, Double_t> pathLength(Double_t r, Double_t x, Double_t y);

  /// path length at distance of closest approach to a given point
  Double_t pathLength(const TVector3& p, Bool_t scanPeriods = true) const;

  /// path length at intersection with plane
  Double_t pathLength(const TVector3& r,
		      const TVector3& n) const;

  /// path length at distance of closest approach in the xy-plane to a given point
  Double_t pathLength(Double_t x, Double_t y) const;

  /// path lengths at dca between two helices
  std::pair<Double_t, Double_t> pathLengths(const Helix&,
				       Double_t minStepSize = 10*micrometer,
				       Double_t minRange = 10*centimeter) const;

  /// minimal distance between point and Helix
  Double_t distance(const TVector3& p, Bool_t scanPeriods = true) const;

  /// checks for valid parametrization
  Bool_t valid(Double_t world = 1.e+5) const { return !bad(world); }
  Int_t  bad(Double_t world = 1.e+5) const;

  /// Move the origin along the Helix to s which becomes then s=0
  virtual void moveOrigin(Double_t s);

  static const Double_t NoSolution;

  /// Return the momentum at origin
  ///    \param bField magnetic field
  TVector3 momentum(Double_t) const;
  /// Return momemtum at S
  TVector3 momentumAt(Double_t, Double_t) const;
  /// Return charge of a particle
  Int_t charge(Double_t)   const;
  /// 2d DCA to x,y point signed relative to curvature
  Double_t curvatureSignedDistance(Double_t x, Double_t y) ;
  /// 2d DCA to x,y point signed relative to rotation
  Double_t geometricSignedDistance(Double_t x, Double_t y) ;
  /// 3d DCA to 3d point signed relative to curvature
  Double_t curvatureSignedDistance(const TVector3&) ;
  /// 3d DCA to 3d point signed relative to rotation
  // Give the vertex argument in meters. Output is in meters, so multiply by 100 to get cm.
  Double_t geometricSignedDistance(const TVector3&) ;

  //Following added by Matt Harasty
  void dca(const TVector3& a_vertex, double &a_dca, double &a_dcaX, double &a_dcaY, double &a_dcaZ);
  TGraph2D* visualizeHelixTrajectory(double pathLengthStart, double pathLengthEnd, int numPoints);


 protected:

  /// Set curvature of the Helix
  void setCurvature(Double_t);	/// performs also various checks
  /// Set phase of the Helix
  void setPhase(Double_t);
  /// Set dip angle of the Helix
  void setDipAngle(Double_t);
  /// Value of S where distance in x-y plane is minimal
  Double_t fudgePathLength(const TVector3&) const;

 protected:
  /// true for straight line case (B=0)
  Bool_t    mSingularity;
  /// starting point of a Helix
  TVector3  mOrigin;
  /// Dip angle
  Double_t  mDipAngle;
  /// Curvature = 1/R
  Double_t  mCurvature;
  /// Phase
  Double_t  mPhase;
  /// -sign(q*B);
  Int_t     mH;

  /// Cos of dip angle
  Double_t mCosDipAngle;
  /// Sin of dip angle
  Double_t mSinDipAngle;
  /// Cos of phase
  Double_t mCosPhase;
  //// Sin of phase
  Double_t mSinPhase;

  //ClassDef(Helix,1)
};

// tests helix class against straight math
// draws a test helix and computes some example DCAs
//void testHelix();



#endif
