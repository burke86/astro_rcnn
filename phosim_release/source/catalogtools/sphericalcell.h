///
/// @package phosim
/// @file sphericalCell.h
/// @brief Class definition for SphericalCell class which genrates ra/dec (or 
/// theta/phi) cells for phosimcatgen star(and later galaxy) generation.
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
#ifndef SPHERICALCELL_H
#define SPHERICALCELL_H

#include <cmath>
#include "../ancillary/random.h"

// **********************************
// Cell width /height definitions.
// **********************************
const double kThetaHeightDeg = 1.0;
const double kPhiWidthDeg    = 1.0;

// **********************************


class SphericalCell
{
public:
  SphericalCell();

  void   GetGalacticFromRaDec(double RADeg, double DecDeg, double& longDeg, 
                            double& latDeg);
  void   SetGalacticFromRaDec();
  double CalculateStarDensity(double M);
  void   SetBoundries(double Theta, double ThetaSpacng, double PhiDeg,
                      double PhiSpacingDeg);
  void   SetAreaDeg2();
  void   GenerateCellRandomRADec(double& RADeg, double& decDeg, 
                                 Random& random);

  double fRADeg;
  double fDecDeg;
  int    fThetaIndex;
  int    fPhiIndex;
  double fThetaMinDeg;  //Filled by SetBoundries (spherical coords)
  double fThetaMaxDeg;
  double fPhiMinDeg;
  double fPhiMaxDeg;
  double fAreaDeg2;
  double fGalacticLongitudeDeg;  //Not used yet
  double fGalacticLatitudeDeg;   //Determines star density
  //Random fRandom;    //seed for this will be a funciton of base seed and 
  //                  //thetaIndex and phiIndex 
};

#endif
