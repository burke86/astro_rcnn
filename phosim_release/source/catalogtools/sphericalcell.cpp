///
/// @package phosim
/// @file sphericalCell.cpp
/// @brief SphericalCell funnctions which generate ra/dec (or 
///  theta/phi) cells for phosimcatgen star(and later galaxy) generation.
///
/// @brief Created by:
/// @author Glenn Sembroski (Purdue)
///
/// @brief Modified by:
/// @author 
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "sphericalcell.h"
const double PAL__DPI  = 3.1415926535897932384626433832795028841971693993751;
// Constructor
SphericalCell::SphericalCell()
{
  //Nothing here yet
}
// ***********************************************************************


void SphericalCell::GetGalacticFromRaDec(double RADeg, double DecDeg,
                                         double& longDeg, double& latDeg)
// *****************************************************************
// Convert from a supplied RA/Dec to specified Galactic Longitude/Latitude
// Using functions suggested by John Peterson who has the references.
// ****************************************************************** 
{
  // Longitude first
  double arg1Rad = (192.25-RADeg) * ( PAL__DPI / 180.);
  double arg2Rad = 27.4 * (PAL__DPI / 180.);
  double decRad  = DecDeg * ( PAL__DPI / 180.);

  double longitudeRad = atan( sin(arg1Rad) / 
          ( (cos(arg1Rad)*sin(arg2Rad) )  - (tan(decRad) * cos(arg2Rad) ) ) ); 
                              
  longDeg = longitudeRad * ( 180. / PAL__DPI );

  double latitudeRad = asin(  ( sin(decRad)*sin(arg2Rad) ) + 
                              ( cos(decRad) * cos(arg2Rad) * cos(arg1Rad) ) ) ;
  latDeg = latitudeRad *  ( 180. / PAL__DPI );

  return;
} 
// ********************************************************************

void SphericalCell::SetGalacticFromRaDec()
// *****************************************************************
// Convert from RA/Dec to Galactic Longitude/Latitude
// Using functions suggested by John Peterson who has the references.
// ****************************************************************** 
{
  GetGalacticFromRaDec(fRADeg,fDecDeg, fGalacticLongitudeDeg, 
                       fGalacticLatitudeDeg);
  return;
} 
// ********************************************************************


double SphericalCell::CalculateStarDensity(double M)
{
  //Assuming M is low edge of the bin.
  double a = -.0004 -.00006*std::abs(fGalacticLatitudeDeg);
  if (a < -.0040) {
    a = -.0040;
  }
  double lambda= -4. + (.5 * M) +  (a * pow(M,2)) - (.00017 * pow(M,3) );
  double density = pow(10., lambda);
  // cout<<" a,M, Theta1Deg,density: " << a << " " << M << " " << latitudeDeg 
  //    << density << endl;
  return density;
}
// **********************************************************************

void SphericalCell::SetBoundries(double ThetaDeg, double ThetaSpacingDeg, 
                                 double PhiDeg, double PhiSpacingDeg)
// *******************************************************
// Set min max of cell boudries in spherical coords.
// *******************************************************
{ 
  fThetaMinDeg=ThetaDeg - ThetaSpacingDeg/ 2.;
  if (fThetaMinDeg<0.0) {
    fThetaMinDeg = 0.0;
  }
  fThetaMaxDeg = ThetaDeg + ThetaSpacingDeg / 2.;
  if (fThetaMaxDeg> 180.0) {
    fThetaMaxDeg = 180.0;
  }
  fPhiMinDeg = PhiDeg  - PhiSpacingDeg / 2.;
  if (fPhiMinDeg<0.0) {
    fPhiMinDeg = 0.0;
  }
  fPhiMaxDeg = PhiDeg + PhiSpacingDeg / 2.;
  if (fPhiMaxDeg> 360.0) {
    fPhiMaxDeg = 360.0;
  }

  return;
}
// ************************************************************************

void SphericalCell::SetAreaDeg2()
{
  //Sr is for steradians
  double capArea1Sr = 2*PAL__DPI*(1 - cos(fThetaMinDeg * (PAL__DPI/180) ) );
  double capArea2Sr = 2*PAL__DPI*(1 - cos(fThetaMaxDeg * (PAL__DPI/180) ) );
  double binPhiFraction = ( (fPhiMaxDeg-fPhiMinDeg)/360.);
  double celAreaSr = (capArea2Sr-capArea1Sr) * binPhiFraction;

  fAreaDeg2 = celAreaSr * ( (180.*180)/ (PAL__DPI* PAL__DPI));

  return;
}
// *********************************************************************

void SphericalCell::GenerateCellRandomRADec(double& RADeg, double& DecDeg, 
                                        Random& random)
{
  double u = random.uniform();
  double v = random.uniform();
  RADeg = fPhiMinDeg + u*(fPhiMaxDeg-fPhiMinDeg);  //phi and RA same
  double r1 = ( cos(fThetaMinDeg * (PAL__DPI/180) ) + 1.)  / 2.;    
  double r2 = ( cos(fThetaMaxDeg * (PAL__DPI/180) ) + 1.)  / 2.;
  double r = r1 + v * (r2-r1);
  double thetaDeg = acos(2*r -1) * (180/PAL__DPI);
  DecDeg = 90.-thetaDeg;
  return;
}
// ***********************************************************************

