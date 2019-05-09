///
/// @package phosim
/// @file GalaxyComplex.h
/// @brief Generates phosim sersicComplex galaxy parameters
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

#ifndef GALAXYCOMPLEX_H
#define GALAXYCOMPLEX_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include <algorithm>    //std::min_element, std::max_element, std:lower_bound
#include <sstream>

#include "../ancillary/random.h"
const double kPI = 3.1415926535897932384626433832795028841971693993751;
const int kNumValues = 10000;   //Number of Mabsolute,Mapparent and z 
                                //divisions of range.

class GalacticComplexObject
{
  // **********************************************************************
  //Structure to hold sersicComplex values for phosimcatgen GalaxyComplex class
  // **********************************************************************
 public:
  std::string ID;
  double RADeg;
  double DecDeg;
  double M;                                         //magnitude
  std::string SEDFileName;
  double z;                                         //Red Shift

  //sersicComplex parameters
  double a;
  double b;
  double c;
  double alpha;
  double beta;
  int    index;
  double clumpFrac;
  double clump;
  double clumpWidth;
  double spiralFrac;
  double alphaSpiral;
  double bar;
  double spiralWidth;
  double phi0;
};
// *****************************************************************

class GalaxyComplex
// ***************************************************************
// This class uses IDL code from John Peterson (2018-02-15), converted to C++
// by Glenn Sembroski. (Not sure of definitions or souces of any of this) 
// It produces a random galaxy, both  the disk and bulge magnitudes, 
// the red shift z
// and it fills the sersicComplex class for each (disk and bulge) in order for 
// the calling program to generate a galaxy object(2 objects, disk and bulge
//  actually) for a phosim catalog.
// ****************************************************************
{
 public:
  GalaxyComplex();
  void Init(); 
  bool GenerateGalaxyForCells(GalacticComplexObject& Bulge, 
                              GalacticComplexObject& Disk, double MinimumM, 
                              double MaximumM, Uint64 galaxySeed);
  double GetGalaxyDensityDeg2(){return fGalaxyDensityDeg2;};
  std::string  GenerateSersicComplexObjectString(
											   GalacticComplexObject& Object);
  void WriteObject(GalacticComplexObject& Object, std::ofstream* ofs);

private:
  double fGalaxyDensityDeg2;
  std::vector<double >::iterator fVectorIteratorUpper;

 public:
  std::vector < double > fMabsoluteValues;
  std::vector < double > fCumulativeMabsoluteProb;
  std::vector < double > fZValues;
  std::vector < double > fCumulativeZProb;
};
// *********************************************************************

#endif
