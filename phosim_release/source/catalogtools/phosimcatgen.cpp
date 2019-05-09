///
/// @package phosim
/// @file phosimcatgen.cpp
/// @brief Generates phosim Object catalog and files with 
///        random or grid of stars and/or galaxies.
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

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <algorithm>    // std::min_element, std::max_element
#include <sstream>

#include "sphericalcell.h"
#include "galaxycomplex.h"
#include "../ancillary/readtext.h"
#include "../ancillary/random.h"


// ***********************************************************************
// Utility function definitions for Sofa routines (see refereneces below)
// ***********************************************************************
double eraAnp(double a);
double eraPdp(double a[3], double b[3]);
double eraPm(double p[3]);
void   eraPxp(double a[3], double b[3], double axb[3]);
void   eraS2c(double theta, double phi, double c[3]);
double eraSepp(double a[3], double b[3]);
double eraSeps(double al, double ap, double bl, double bp);
void   palDtp2s ( double xi, double eta, double raz, double decz,
                double *ra, double *dec );
void   palDh2e ( double az, double el, double phi, double *ha, double *dec);

//Constants for pal* and era* functions (and also used as needed otherwise)
const double ERFA_D2PI = 6.283185307179586476925287;
const double PAL__DPI  = 3.1415926535897932384626433832795028841971693993751;
// ***********************************************************************


const int kPrintFreq = 10000;

// const double kBulgeToTotalFlux      = .5541;
// const double kBulgeToTotalFluxSigma = .0198;


// **************************************************************************
class  CGParams
{
public:
  CGParams(int argc, char* argv[]);
  void SetDefaults();
  void PrintOptionValues();
  bool makeStars(){return fMakeStars;};
  bool makeStarGrid(){return fMakeStarGrid;};
  bool makeGalaxies(){return fMakeGalaxies;};

  std::string fObservationID;
  uint32_t fRandomSeed;
  double   fRADegOrigin;
  double   fDecDegOrigin;


  //command arguments;
  double fStarsMinimumMagnitude;
  double fStarsMaximumMagnitude;
  double fStarsDiameterFOVDeg;

  double fStarGridSpacingDeg;
  double fStarGridMagnitude;
  double fStarGridWidthDeg;

  double fGalaxiesMinimumMagnitude;
  double fGalaxiesMaximumMagnitude;
  double fGalaxiesDiameterFOVDeg;


private:
  bool fMakeStars;
  bool fMakeStarGrid;
  bool fMakeGalaxies;

public:
  int           fCatalogObjectID;
  int           fCatalogObjectCount;
  double        fStarDensityPerDeg2;
  double        fGalaxyDensityPerDeg2;
  std::string   fStarGridSEDFileName;
  std::string   fCatFileName;
  std::ofstream fCatFile;

  //double fMinimumMagnitude;
  //double fMaximumMagnitude;
  //std::vector< double > fStarMagDensity;
  //std::vector< double > fGalaxyMagDensity;
  //std::vector< double > fStarMagCumulativeDensity;
  //std::vector< double > fGalaxyMagCumulativeDensity;
  //std::vector< double > fStarMagCumulativeProb;
  //std::vector< double > fGalaxyMagCumulativeProb;

};
// ***********************************************************************
  
CGParams::CGParams(int argc, char* argv[])
{
  // *****************************************************************
  // This code has little to no error checking. Assumed this program is run 
  // by another program (phosim.py usually)  which knows how things should
  // be setup
  // *****************************************************************

  SetDefaults();

  // This constructer parses the parametersa form the input (obs_*.pars) 
  // file comming in though stdin.

  std::string cmdArguments;

  // Read parameters from stdin.
  readText pars(std::cin);
  for (size_t t(0); t < pars.getSize(); t++) {
	std::string line(pars[t]);
	readText::get(line, "obshistid",   fObservationID);
	readText::get(line, "obsseed",     fRandomSeed);
	readText::get(line, "pointingra",  fRADegOrigin);
	readText::get(line, "pointingdec", fDecDegOrigin);

	std::istringstream iss(line);
	std::string keyName;
	iss >> keyName;

	if ( keyName == "stars") {
	  fMakeStars = true;
	  iss >> fStarsMinimumMagnitude >> fStarsMaximumMagnitude 
		  >> fStarsDiameterFOVDeg;
	}

  	if ( keyName == "stargrid") {
	  fMakeStarGrid = true;
	  iss >> fStarGridSpacingDeg >> fStarGridWidthDeg >> fStarGridMagnitude;
	}

  	if ( keyName == "galaxies") {
	  fMakeGalaxies = true;
	  iss >> fGalaxiesMinimumMagnitude >> fGalaxiesMaximumMagnitude 
		  >> fGalaxiesDiameterFOVDeg;
	}
  }
 }
// **********************************************************************

void CGParams::PrintOptionValues()
{
  std::cout<< " obervation ID: " << fObservationID << std::endl;
  std::cout<< " ra(deg): " << fRADegOrigin << std::endl;
  std::cout<< " dec(deg):  "<< fDecDegOrigin << std::endl;
  std::cout<< " random number generator seed: " << fRandomSeed << std::endl;
  if (fMakeStars) {
	std::cout << " make stars: True"      << std::endl;
	std::cout << "  min stars magnitude: " << fStarsMinimumMagnitude 
			  << std::endl;
	std::cout << "  max stars magnitude: " << fStarsMaximumMagnitude 
			  << std::endl;
	std::cout << "  stars FOV Diameter(deg): " << fStarsDiameterFOVDeg 
			  << std::endl;
    std::cout << "  stars FOV Area(deg^2): " << PAL__DPI * 
                                  pow(fStarsDiameterFOVDeg/2.,2) << std::endl;
    
  }
  else{
	std::cout << " make stars: False" << std::endl;
  }

  if (fMakeStarGrid) {
	std::cout<< " make stargrid: True" << std::endl;
	std::cout<< "  stargrid Spacing(deg): " << fStarGridSpacingDeg 
			 << std::endl;
	std::cout<< "  stargrid Magnitude(deg): " << fStarGridMagnitude 
			 << std::endl;
	std::cout<< "  stargrid FOV Diameter(deg): " << fStarGridWidthDeg 
			 << std::endl;
  }
  else{
	std::cout<< " make stargrid: False"<< std::endl;
  }
  
  if (fMakeGalaxies) {
	std::cout<< " make galaxies: True" << std::endl;
	std::cout<< "  min galaxies magnitude: " << fGalaxiesMinimumMagnitude 
			 << std::endl;
	std::cout<< "  max galaxies magnitude: " << fGalaxiesMaximumMagnitude 
			 << std::endl;
	std::cout<< "  galaxies FOV Diameter(deg): " << fGalaxiesDiameterFOVDeg 
			 << std::endl;
    std::cout<< "  galaxies FOV Area(deg^2): " << PAL__DPI * 
                              pow(fGalaxiesDiameterFOVDeg/2.,2) << std::endl;
  }
  else{
	std::cout<< " make galaxies: False" << std::endl;
  }
  
  std::cout << "-----------------------------------------------------------"
	"-------------------------------" << std::endl;
  return;
}
// **********************************************************************

void CGParams::SetDefaults()
{
  fMakeStars    = false;
  fMakeStarGrid = false;
  fMakeGalaxies = false;
  fStarGridSEDFileName = "../sky/sed_flat.txt";
  return;
}
// ************************************************************************



class CatalogObject
{
public:
  //CatalogObject(CGParams* pCGPar , Random* pRandom);
  CatalogObject(CGParams* pCGPar );
  void   GenerateSingleRandomRaDec();
  void   GenerateGridRaDec(double RARad,  double DecRad, 
				  double FOVDeg, double StepSizeDeg);
  double GetRADeg(){return fRARad*180.0/PAL__DPI;};
  double GetDecDeg(){return fDecRad*180.0/PAL__DPI;};
  void   SetMagnitude(double mag ){fMagnitude=mag;return;};
  double GetMagnitude(){return fMagnitude;};
  void   SetStarsSEDFileName()
              {fSEDFileName=pCGParams->fStarGridSEDFileName;return;};
  void   SetStarGridSEDFileName()
              {fSEDFileName=pCGParams->fStarGridSEDFileName;return;};
  void   SetGalaxySEDFileName()
              {fSEDFileName=pCGParams->fStarGridSEDFileName;return;};
  std::string GetSEDFileName(){return fSEDFileName;};
  void   GetSphericalCellsInFOV(double FOVDeg); 

                                          
  bool        fDebugPrint;

protected:
  CGParams* pCGParams;
  //Random* pRandom;
  Random fRandom;
  std::vector < double > fRightAscension;
  std::vector < double > fDeclination;

  //std::vector < double > fStarMagDensity;
  //std::vector < double > fGalaxyMagDensity;
  //std::vector < double > fStarMagCumulativeDensity;
  //std::vector < double > fGalaxyMagCumulativeDensity;
  //std::vector < double > fStarMagCumulativeProb;
  //std::vector < double > fGalaxyMagCumulativeProb;

  std::vector < SphericalCell > fCells;
  int fNumThetaCells;
  int fNumPhiCells;
  //Object vaiables
  double fRARad;
  double fDecRad;
  double fMagnitude;
  double fZenithMaxRange;
  std::string fSEDFileName;
  std::string fStarSEDFileName;
  std::string fStarGridSEDFileName;
  std::string fGalaxySEDFileName;
 

}; 
// ***************************************************************************

//CatalogObject::CatalogObject(CGParams* pCGPar, Random* prandom)
CatalogObject::CatalogObject(CGParams* pCGPar)
{
  pCGParams =  pCGPar;
  //pRandom   =  prandom;
  // Need to set random seed  before we use random numbers
  fRandom.setSeed32( pCGParams->fRandomSeed); 

}
// ********************************************************************** 

void CatalogObject::GenerateSingleRandomRaDec()
{
  // ****************************************************************
  // We have to be a little clever in generating the random 
  // palcement. 
  // From:   // http://mathworld.wolfram.com/SpherePointPicking.html
  // 1: Pick randomly in an alt/az over a spherical cap centered at
  //  zenith. alt=90,az=0,0.
  // 2:Use the conversion program to convert the az/el to the 
  // correct ra/dec
  // ****************************************************************
  // Gen random az,el in cap with fZenithMaxRange. Set in Constructor
  // ***************************************
  //double azimuthRad=pRandom->uniform()  * 2 * PAL__DPI;  //Uniform in azimuth
  double azimuthRad=fRandom.uniform()  * 2 * PAL__DPI;  //Uniform in azimuth
                                          //zenith 0 to zenithMaxRange in cap
  //double zenithRad=acos( 1. - 2. * pRandom->uniform() * fZenithMaxRange);
  double zenithRad=acos( 1. - 2. * fRandom.uniform() * fZenithMaxRange);
  double elevationRad = PAL__DPI / 2. - zenithRad;
  
  double longRad = pCGParams->fRADegOrigin * PAL__DPI / 180.0;
  double latRad  = pCGParams->fDecDegOrigin * PAL__DPI / 180.0;
  double haRad;
 
  // *********************************
  // Use generic conversion function from Slalib (From pallib actually for 
  // licensing issues)
  // ******************************** 

  palDh2e(azimuthRad,elevationRad ,latRad ,&haRad, &fDecRad );
  fRARad=longRad-haRad;
  if (fRARad<0.0) {
	fRARad= 2.0*PAL__DPI+fRARad;
  }
  return;
}
// ****************************************************************


void CatalogObject::GenerateGridRaDec(double RARad,  double DecRad, 
									  double FOVDeg, double StepSizeDeg)
// *****************************************************************
// Generates a square grid of focalplane locations(in deg) and 
// converts to equivalent ( tangent plane or gnomonic) projection to
// the celestial sphere in ra/dec centered at the specified ra/dec
// *****************************************************************
// Signs of directions will be determined later. Just go with 
// "whatever" for now.
// ****************************************************************
// Make up a vector of vectors grid of focal plane locations. 
// Outermost is "Y" which is along altitude coord.
// Results in Catalog members: "rightAscension" and "declination" vectors.
// ****************************************************************
{
  bool debugPrint = false;
  fRightAscension.clear();
  fDeclination.clear();


  // ****************************************************************
  // Determin xPintDeg,yPointDeg  (in focal plane);
  // *******************************
  
  int npoints =  FOVDeg / StepSizeDeg;

  // Since we want a point at the very center this should be odd
  if (npoints%2 == 0) {
    npoints=npoints+1;
  }

  std::cout << " #Star Grid spacing(deg):" << StepSizeDeg 
			<< std::endl;
  std::cout << " #Number of points across grid: " << npoints << std::endl;
  std::cout << " #Number of points in grid:     " << npoints*npoints 
			<< std::endl;
  
  if (npoints*npoints > kPrintFreq) {
	std::cout << " This may take a few minutes. Be patinet!" << std::endl;
	std::cout <<" Grid Ra/Dec creation progress: each # = " << kPrintFreq 
			  <<" grid points completed" << std::endl;
	std::cout<< " " << std::flush; //Start off a little offset to lineup with
	                               // above 
  }

  double halfWidth=( ( npoints - 1 ) /2 ) * StepSizeDeg;

  for ( int j = 0; j < npoints; j++ ) {

    double yPointDeg= (- halfWidth) + j * StepSizeDeg;

    for ( int i = 0 ;i < npoints ; i++ ) {

      double xPointDeg= (- halfWidth) + i * StepSizeDeg;

 	  if (debugPrint) {
		std::cout << "xPointDeg, yPointDeg: "<<  xPointDeg << " " << yPointDeg
				  << std::endl;
	  }

	  // *******************************************************
	  // Project xPointDeg and yPointDeg back to celestial sphere and 
	  // get Ra/dec for this grid point.
	  // ********************************************************
      double raGridRad;
      double decGridRad;

	  // **********************************************************
	  // Use palDtp2s. (Tangent Plane to Spherical)  gnomonic projection
	  // function which we coppied to here. (Allowed under pallib/erfam
	  // license)
	  // **********************************************************

	  palDtp2s(xPointDeg * PAL__DPI / 180., yPointDeg *  PAL__DPI / 180., 
			   RARad,  DecRad, &raGridRad,  &decGridRad);

      // *****************************
	  // Save the ra/dec of the grid point.
	  // ****************************
	  fRightAscension.push_back(raGridRad);
      fDeclination.push_back(decGridRad);
 
	  if (debugPrint) {
		std::cout<<xPointDeg<< " " <<  yPointDeg << " " << raGridRad << " "
			   << decGridRad<<std::endl;
		std::cout << "at point: " <<  fRightAscension.size() 
				  << ": raGridRad,decGridRad: "<<  raGridRad << " " 
				  << decGridRad << std::endl;
	  }

	  if ( fRightAscension.size()%kPrintFreq == 0) {
		std::cout<< "#" << std::flush;
	    if  ( fRightAscension.size()%(kPrintFreq*10) == 0) {
		  std::cout<< "|" << std::flush;
		}
	  }
    }
  }
  return;
}
// ****************************************************************


void CatalogObject::GetSphericalCellsInFOV(double FOVDiameterDeg) 
{
  // ******************************************************************
  // Find all the ra/dec cells that have any area in the FOV (centered at 
  // pCGParams->fRADegOrigin,pCGParams->fDecDegOrigin). As we find each cell 
  // place it in the vector "Cells". 
  // ******************************************************************
  // Basic approcch is to divide the celestial sphere (using theta,phi coords)
  // into cells. The height of each cell is ~kThetaHeightDeg. The width is 
  // ~kPhiWidthDeg (though that may change in the future to being a width 
  // inversly propotional to cos(Dec).  But, keep it simple for now)
  // Note: we make fine adjustments to the cell spacing so that they exactly 
  // cover all of the clestial sphere.
  // ******************************************************************
  // Use a brute force search. First increase the FOV radius by the larger of 
  // Phi angular Spacing or the Theta angular Spacing. Then, test each and 
  // every cell to see if it's center is within this expanded FOV. We can use 
  // the algorithum used for the SlaLib function SLA_DSEP to determine if the 
  // center of a cell is within the expanded FOV.
  // ********************************************************************

  //Theta,phi of FOV origin
  double originThetaDeg = 90.-pCGParams->fDecDegOrigin;
  double originPhiDeg =  pCGParams->fRADegOrigin;


  // make a grid on the sphere. Make sure we cover it all
  fNumThetaCells = 180/kThetaHeightDeg;  // note psossible round down
  fNumPhiCells   = 360/kPhiWidthDeg;
  
  // Adjust bin spacing to exactly fit(in case we had round down).
  // Note for phi the actutal spacing in distance on the sphere, the angular 
  // distance, is largest at dec=0 (theta=90.)
  double thetaCellSpacingDeg= 180./(double)fNumThetaCells;
  double phiCellSpacingDeg = 360./(double)fNumPhiCells;

  // ***********************************************************************
  // We need to find all those Theta,Phi cells that have any area in the FOV. 
  // To do this we add the larger of Theta angular distance spacing or phi 
  // angular distance spacing at dec=0 (theta=90.) to the FOV angular radius. 
  // This is to insure we find all possible cells, This is obviously way 
  // overkill (phi angular distance spacing is largest at the equator (dec=0)) 
  // I know,  but I think it results in us not missing anything.
  // ***********************************************************************
  double modifiedFOVRadiusDeg = FOVDiameterDeg/2;
  if (thetaCellSpacingDeg > phiCellSpacingDeg) {
    //Add another 10% for good luck
    modifiedFOVRadiusDeg += 1.1*thetaCellSpacingDeg; 
  }
  else{
    modifiedFOVRadiusDeg += 1.1*phiCellSpacingDeg;
  }   

  // Now iterate through ALL the cells, testing each one to see if it might
  // be in the FOV.
  // Following  could be more efficent(could make some obvious cuts on theta 
  // range for example) but this is clearer.
  for (int i=0; i < fNumThetaCells; i++ ) {
    double thetaDeg =   (i * thetaCellSpacingDeg ) + thetaCellSpacingDeg / 2.;
    for (int j=0; j <  fNumPhiCells; j++ ) {
      // Find Phi of the center of this cell
      double phiDeg = (j * phiCellSpacingDeg) + phiCellSpacingDeg / 2.;

      //See if cells center is within expanded FOV. Use Sofa eraSeps function
      double angDistDeg =  eraSeps(phiDeg * (PAL__DPI/180.), 
                                   (90.-thetaDeg) * (PAL__DPI/180.), 
                                   originPhiDeg *(PAL__DPI/180.), 
                                   (90.-originThetaDeg) *(PAL__DPI/180.) ) * 
                                                               (180./PAL__DPI);
      if (angDistDeg <= modifiedFOVRadiusDeg) {
        //This cell is in the field of view (FOV). Fill its values

        SphericalCell cell;
        cell.fRADeg      = phiDeg;            //Phi same as RA in deg
        cell.fDecDeg     = 90.-thetaDeg;     //Theta = 90-dec
        cell.fThetaIndex = i;
        cell.fPhiIndex   = j;
        // **********************
        // Now find area of this cell in deg**2
        // Get cell boundries
        // **********************
        cell.SetBoundries(thetaDeg, thetaCellSpacingDeg, phiDeg, 
                                                           phiCellSpacingDeg);
        cell.SetAreaDeg2();
        cell.SetGalacticFromRaDec();
 
        // And save this "hit" cell in the Catalog object vector fCells
        fCells.push_back(cell);
      }
    }
  }
  return;
}
// ***********************************************************************


class StarObject :public CatalogObject
// **********************************************************
// Note this class derived from CatalogObject and thus has 
// CatalogObject's methods
// **************************************************
{
 public:
  //StarObject(CGParams* pCGPar, std::ofstream& objFile, Random* prandom);
  StarObject(CGParams* pCGPar, std::ofstream& objFile);
  void GenerateGridStars();
  void GenerateRaDecFromStarXYMFile(std::string fileName);
  void GenerateGalacticStars(); 
private:
  std::ofstream* fOfs;
};
// *****************************************************************


StarObject::StarObject(CGParams* pCGPar, std::ofstream& objFile): 
                         CatalogObject(pCGPar)
  //StarObject::StarObject(CGParams* pCGPar, std::ofstream& objFile, 
  //Random* prandom):
  //                            CatalogObject(pCGPar,prandom)  
{
  //fZenithMaxRange used to generate  ra/dec for random stars.
  fZenithMaxRange = ( 1. - cos( ( pCGParams->fStarsDiameterFOVDeg / 2. )
							   * PAL__DPI / 180.0 ) ) / 2.0;

  fOfs = &objFile;
}


void StarObject::GenerateGridStars()
{
  // *******************************************************
  //Using a single default (for now) SED for all grid stars
  // *******************************************************
  SetStarGridSEDFileName();              //CatalogObject Method

  // ***************************************************************
  //Get Ra/Dec of where the tel is pointing
  // ***************************************************************
  double raRad = pCGParams->fRADegOrigin * PAL__DPI / 180.0;
  double decRad  = pCGParams->fDecDegOrigin * PAL__DPI / 180.0;


  // **************************************************************
  //Generate the grid locations, centered on raRad,decRad.
  // Places ra/dec values for grid in rightAscension and declination
  // vectors.
  // **************************************************************
  GenerateGridRaDec(raRad, decRad, pCGParams->fStarGridWidthDeg, 
					       pCGParams->fStarGridSpacingDeg);   
                                                 //CatalogObject method
  SetMagnitude(pCGParams->fStarGridMagnitude);
  
  int numRaDecPoint=fRightAscension.size();
  std::cout << "Making " << numRaDecPoint << "Star Grid Objects" << std::endl;

  // ************************************
  // Write star object to catalog file
  // ************************************
  for ( int i = 0; i < numRaDecPoint; i++ ) {
    pCGParams->fCatalogObjectID++;
    *fOfs<< "object " << pCGParams->fCatalogObjectID << " " << std::fixed 
       << std::setprecision(8) << fRightAscension.at(i) * 180.0 / PAL__DPI 
       << " " << fDeclination.at(i) * 180.0 / PAL__DPI 
       << "  " << std::setprecision(6) << GetMagnitude()
       << " "  <<  GetSEDFileName() 
       << " 0.0 0.0 0.0 0.0 0.0 0.0 star none none" 
       << std::endl;
  }
  pCGParams->fCatalogObjectCount += numRaDecPoint;
  return;
}
// ****************************************************************

void StarObject::GenerateGalacticStars()
// ****************************************************************
// Generate stars randomly (but repeatably) with a density that is
// galactic latitude dependent (may include longitude dependency later)
// Do this by dividing celestial sphere into ra and dec cells whithin which 
// density is aproximatly constant. Use density and area of cell to get 
// mean number in cell, fluxtuate it with approprate poison distribution and 
// create and wirte object stars to catalog file.
// ***************************************************************
{
  SetStarsSEDFileName();//This just sed_flat.txt until we find a better way.

  // ****************************************************
  // Find all the cells in our FOV (and then some for safty's sake)
  // Places cells in vector fCells
  // ****************************************************
  GetSphericalCellsInFOV(pCGParams->fStarsDiameterFOVDeg);

  // *****************************************************
  // Iterate through the cells
  int numCells = fCells.size();
  for (int i = 0; i < numCells; i++) {
    // We will now need to iterate though the Magnitude range since density is 
    // a function of magnitude. Use deltaM=1 steps
    // Since our densitiy is for magnitude bins with integer low edge 
    // and high edge, handle partial lowest and highest bins.

    int lowMEdge=int(pCGParams->fStarsMinimumMagnitude);  //Round DOWN!

    Random random;
    // Determine random seed for this cell. Using a different but cell 
    // index dependent seed for each cell insures we can reproduce same 
    // stars independent of origin of FOV
    // Set the seed in this cell random number generator
    //Note this way every cell will always have the same stars in it.
    int cellSeed = fCells.at(i).fPhiIndex +
                                      fCells.at(i).fThetaIndex * fNumPhiCells; 
    random.setSeed32(cellSeed);
    //std::cout << "i, cellSeed: " << i << " " << cellSeed <<std::endl;

    while ( lowMEdge < pCGParams->fStarsMaximumMagnitude ) {
      // ***************************************************
      // Mean number of stars in this Magnitude bin for this cell
      // Using low edge of bin here. 
      // Possible lowest and hists bins are fractional. Handle low edge first
      // We only reset m here if lowMEdge is not the minimum magnitude in this
      // bin. Otherwise we use lowMEdge for m.
      // ****************************************************
      double mBinFraction = 1;
      double m = (double) lowMEdge;
      if (double(lowMEdge) < pCGParams->fStarsMinimumMagnitude ) {
        m=pCGParams->fStarsMinimumMagnitude;
        mBinFraction = (lowMEdge+1)-pCGParams->fStarsMinimumMagnitude;
      }
      else  if (lowMEdge+1. >pCGParams->fStarsMaximumMagnitude ) {
        //Partial top Magnitude cell. Scale things as needed
        mBinFraction = (pCGParams->fStarsMaximumMagnitude-lowMEdge);
        //use defulat low edge (above)
      }

      double density = fCells.at(i).CalculateStarDensity(m);
      double meanNumStars = density * fCells.at(i).fAreaDeg2 * mBinFraction;

      //Possion fluxtuate this number to get correct statistical number
      int numStars = random.poisson(meanNumStars);
      
      //Now generate these stars and place in the catalog
      if (numStars >0 ) {
        for (int j = 0; j < numStars;  j++ ) {
          double RADeg=0.;
          double decDeg=0.;
          //fCells.at(i).GenerateRandomRADec(RADeg, decDeg);
          fCells.at(i).GenerateCellRandomRADec(RADeg, decDeg, random);
          
          // At this point we need to check that this particular star is within 
          // the exact FOV. (not the modified). This saves time and space.
          double angDistDeg = eraSeps( 
                                 RADeg * (PAL__DPI/180.), 
                                 decDeg * (PAL__DPI/180.), 
                                 pCGParams->fRADegOrigin *(PAL__DPI/180.), 
                                 pCGParams->fDecDegOrigin *(PAL__DPI/180.) ) * 
                                                             (180./PAL__DPI);

          // *************************************************************
          // In order for the random number stream to remain consistent for
          // this cell we need to pre calculate the random nnumber we would 
          // normaly use after the  following test
          double u = random.uniform();
          if ( angDistDeg <= pCGParams->fStarsDiameterFOVDeg/2.) { 
            //Do an even distribution within magnitude bin for M.
            //double magnitude= m * fCells.at(i).fRandom.uniform() * 
            //Handle special cases of first and last bin by setting m and 
            // mbinFraction as appropriate above.
            double magnitude;
            magnitude = m + (u * mBinFraction);

            //Debug print
            //double longDeg;
            //double latDeg;
            //fCells.at(i).GetGalacticFromRaDec(RADeg, decDeg, longDeg, latDeg);
            //std::cout << RADeg << " " << decDeg << " " << magnitude << " " 
            //          << mBinFraction << " " << latDeg << " " << longDeg 
            //          << std::endl;

            // Write star to catalog 
            pCGParams->fCatalogObjectID++;
            *fOfs<< "object " << pCGParams->fCatalogObjectID << " " 
                 << std::fixed << std::setprecision(8) << RADeg << " " 
                 << decDeg 
                 << "  " << std::setprecision(6) << magnitude << " "  
                 <<  GetSEDFileName() 
                 << " 0.0 0.0 0.0 0.0 0.0 0.0 star none none" << std::endl;
            pCGParams->fCatalogObjectCount++;        
          }
        }
      }
      //Bump to next magnitude bin.
      lowMEdge++;
    }
  }
  return;
}  
// **************************************************************************

class GalaxyObject :public CatalogObject
// **********************************************************
// Note this class derived from CatalogObject and thus has 
// CatalogObject's methods
// **************************************************
{
 public:
  //GalaxyObject(CGParams* pCGPar, std::ofstream& objFile, Random* pRandom);
  GalaxyObject(CGParams* pCGPar, std::ofstream& objFile);
  void GenerateComplexGalaxiesForCells();

private:
  std::ofstream* fOfs;
};
// *****************************************************************

GalaxyObject::GalaxyObject(CGParams* pCGPar, std::ofstream& objFile) :
                CatalogObject(pCGPar)  
//GalaxyObject::GalaxyObject(CGParams* pCGPar, std::ofstream& objFile, 
//                           Random* pRandom):
//  CatalogObject(pCGPar, pRandom)  
{

  //fZenithMaxRange used to generate  ra/dec
  fZenithMaxRange = ( 1. - cos( ( pCGParams->fGalaxiesDiameterFOVDeg / 2. )
							   * PAL__DPI / 180.0 ) ) / 2.0;
  fOfs = &objFile;
}
// **************************************************************************

void GalaxyObject::GenerateComplexGalaxiesForCells()
// ****************************************************************
// Generate sersicComplex galaxies randomly (but repeatably) using idl code 
// from John Peterson converted to C++
// **************************************************************
// Do this by dividing celestial sphere into ra and dec cells  Use density and
// area of cell to get  mean number of galaxies in a cell, fluxtuate it with 
// approprate poison distribution and create and write object disk and bulge 
// of galaxies to catalog file.
// ***************************************************************
{
  SetGalaxySEDFileName();//This just flat for now until we find a better way. 

  // ****************************************************
  // Find all the cells in our FOV (and then some for safty's sake)
  // Places cells in vector fCells
  // ****************************************************
  GetSphericalCellsInFOV(pCGParams->fGalaxiesDiameterFOVDeg);

  // Define a GalaxyComplex object which sets up all the cumulative 
  // probability arrays for galaxy density (and other important things). 

  GalaxyComplex complexGalaxy;
  double densityDeg2 =  complexGalaxy.GetGalaxyDensityDeg2();

  //std::cout << "densityDeg2:" << densityDeg2 <<std::endl;


  GalacticComplexObject disk;
  GalacticComplexObject bulge;

  int numGalaxiesInRange=0;

  // *****************************************************
  // Iterate through the ra/dec cells that were found in the FOV
  int numCells = fCells.size();
  for (int i = 0; i < numCells; i++) {
    // Generate a unique but repeatable random number sequence generator fot
    // this cell
    Random random;
    
    // Determine random seed for this cell. using a different,
    // but cell index dependent seed for each galaxy. This will 
    // insure we can reproduce the same galaxies independent of (tracking) 
    // origin of the FOV.
    // Set the seed in this cell's  random number generator
    // Note: You find fNumPhiCells in CatalogObject.
    Uint64 cellSeed =   fCells.at(i).fPhiIndex 
                                   + fCells.at(i).fThetaIndex * fNumPhiCells;
    random.setSeed64(cellSeed);
    
    // Determine for this cell the number of galaxies over the complete 
    // Absolute Magnitude range. We will cut on the Apparent Magnitude later.

    double meanNumGalaxies = densityDeg2 * fCells.at(i).fAreaDeg2;
    //std::cout << i << " " << fCells.at(i).fAreaDeg2 << std::endl;

    //Possion fluxtuate this number to get correct statistical number
    int numGalaxies = random.poisson(meanNumGalaxies);
 
    // Now generate these galaxies and if they are within the specified 
    // magnitude range, we write to the disk and bulge objects to the catalog 
    // file
    if (numGalaxies >0 ) {
      double RADeg=0.;
      double decDeg=0.;
      for (int j = 0; j < numGalaxies;  j++ ) {
          fCells.at(i).GenerateCellRandomRADec(RADeg, decDeg, random);
          
          // At this point we need to check that this particular Galaxy is 
          // within the exact FOV. (not the modified). 
          // This saves time and space.
          double angDistDeg = eraSeps( 
                                 RADeg * (PAL__DPI/180.), 
                                 decDeg * (PAL__DPI/180.), 
                                 pCGParams->fRADegOrigin *(PAL__DPI/180.), 
                                 pCGParams->fDecDegOrigin *(PAL__DPI/180.) ) * 
                                                              (180./PAL__DPI);
          
          if ( angDistDeg <= pCGParams->fGalaxiesDiameterFOVDeg/2.) { 

            //Now try to generate the galaxy objects (2 of them, a disk and 
            // a bulge). Reject if the galaxy Mapparent magnitude is not in
            // the range. GenerateGalaxyFor Cells does this.
            // ************************

            // **************************************************************
            // This gets a little triky. To insure repoducability within cells.
            // We will produce a new random function with a seed dependent on 
            // cell indices and galaxy index (j).
            // Do this so that the random number stream remains consistent.

            Uint64 galaxySeedC = 
                  ( (Uint64) j * (Uint64) fNumPhiCells * (Uint64) fNumThetaCells);
            Uint64 galaxySeed = cellSeed + galaxySeedC;
            //std::cout << galaxySeed << " " << galaxySeedC << " " << cellSeed 
            //          << " " << j << " " 
            //          << fCells.at(i).fPhiIndex   << " " << fNumPhiCells << " " 
            //          << fCells.at(i).fThetaIndex << " " << fNumThetaCells 
            //          << std::endl;

            //This here for debug print in GenerateGalaxyForCells method.
            bulge.RADeg  = RADeg;   
            bulge.DecDeg = decDeg;
            disk.RADeg   = RADeg;
            disk.DecDeg  = decDeg;
            bool weHaveAGalaxy = complexGalaxy.GenerateGalaxyForCells( 
                                   bulge, disk,
                                   pCGParams->fGalaxiesMinimumMagnitude,
                                   pCGParams->fGalaxiesMaximumMagnitude,
                                   galaxySeed);
            if( weHaveAGalaxy) {

              pCGParams->fCatalogObjectID++;   //Bump up the object ID
 
	
              std::ostringstream osBulgeID;
              osBulgeID << pCGParams->fCatalogObjectID << ".1";
              bulge.ID=osBulgeID.str();
              bulge.SEDFileName = GetSEDFileName();
              complexGalaxy.WriteObject( bulge, fOfs);

              std::ostringstream osDiskID;
              osDiskID << pCGParams->fCatalogObjectID << ".2";
              disk.ID = osDiskID.str();
              disk.SEDFileName = GetSEDFileName();
              complexGalaxy.WriteObject( disk,  fOfs);
              
              numGalaxiesInRange++;
            }   //end of magnitude test
          }  //end of ra/dec in FOV test
	  } //end of galaxies in cell loop
    }
  }   //end of cell loop

  //count the objects we add to catalog.
  //Each galaxy requirtes 2 objects(disk and bulge)
  std::cout<< "Number of galactic objects (Disk and Bulge) within requested "
              "magnitude range added to catalog: " << numGalaxiesInRange
           << std::endl;
  pCGParams->fCatalogObjectCount += 2*numGalaxiesInRange;
  return;
}
// ***********************************************************************

// **************************************************************************
// All of the following pal* and era* functions are used under the following 
// copyright and license.
// Included explicitly here to minimize dependencies.
// ***************************************************************************

/*
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  All rights reserved.
**  
**  This library is derived, with permission, from the International
**  Astronomical Union's "Standards of Fundamental Astronomy" library,
**  available from http://www.iausofa.org.
**  
**  The ERFA version is intended to retain identical functionality to
**  the SOFA library, but made distinct through different function and
**  file names, as set out in the SOFA license conditions.  The SOFA
**  original has a role as a reference standard for the IAU and IERS,
**  and consequently redistribution is permitted only in its unaltered
**  state.  The ERFA version is not subject to this restriction and
**  therefore can be included in distributions which do not support the
**  concept of "read only" software.
**  
**  Although the intent is to replicate the SOFA API (other than
**  replacement of prefix names) and results (with the exception of
**  bugs;  any that are discovered will be fixed), SOFA is not
**  responsible for any errors found in this version of the library.
**  
**  If you wish to acknowledge the SOFA heritage, please acknowledge
**  that you are using a library derived from SOFA, rather than SOFA
**  itself.
**  
**  
**  TERMS AND CONDITIONS
**  
**  Redistribution and use in source and binary forms, with or without
**  modification, are permitted provided that the following conditions
**  are met:
**  
**  1 Redistributions of source code must retain the above copyright
**    notice, this list of conditions and the following disclaimer.
**  
**  2 Redistributions in binary form must reproduce the above copyright
**    notice, this list of conditions and the following disclaimer in
**    the documentation and/or other materials provided with the
**    distribution.
**  
**  3 Neither the name of the Standards Of Fundamental Astronomy Board,
**    the International Astronomical Union nor the names of its
**    contributors may be used to endorse or promote products derived
**    from this software without specific prior written permission.
**  
**  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
**  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
**  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
**  FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
**  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
**  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
**  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
**  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
**  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
**  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
**  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
**  POSSIBILITY OF SUCH DAMAGE.
**  
*/
//#include "erfa.h"

double eraAnp(double a)
/*
**  - - - - - - -
**   e r a A n p
**  - - - - - - -
**
**  Normalize angle into the range 0 <= a < 2pi.
**
**  Given:
**     a        double     angle (radians)
**
**  Returned (function value):
**              double     angle in range 0-2pi
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double w;

   w = fmod(a, ERFA_D2PI);
   if (w < 0) w += ERFA_D2PI;

   return w;

}

double eraPdp(double a[3], double b[3])
/*
**  - - - - - - -
**   e r a P d p
**  - - - - - - -
**
**  p-vector inner (=scalar=dot) product.
**
**  Given:
**     a      double[3]     first p-vector
**     b      double[3]     second p-vector
**
**  Returned (function value):
**            double        a . b
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double w;

   w  = a[0] * b[0]
      + a[1] * b[1]
      + a[2] * b[2];

   return w;

}
// ***************************************************************************

double eraPm(double p[3])
/*
**  - - - - - -
**   e r a P m
**  - - - - - -
**
**  Modulus of p-vector.
**
**  Given:
**     p      double[3]     p-vector
**
**  Returned (function value):
**            double        modulus
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   return sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );

}
// ***************************************************************************

void eraPxp(double a[3], double b[3], double axb[3])
/*
**  - - - - - - -
**   e r a P x p
**  - - - - - - -
**
**  p-vector outer (=vector=cross) product.
**
**  Given:
**     a        double[3]      first p-vector
**     b        double[3]      second p-vector
**
**  Returned:
**     axb      double[3]      a x b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double xa, ya, za, xb, yb, zb;

   xa = a[0];
   ya = a[1];
   za = a[2];
   xb = b[0];
   yb = b[1];
   zb = b[2];
   axb[0] = ya*zb - za*yb;
   axb[1] = za*xb - xa*zb;
   axb[2] = xa*yb - ya*xb;

   return;

}
// ***************************************************************************

void eraS2c(double theta, double phi, double c[3])
/*
**  - - - - - - -
**   e r a S 2 c
**  - - - - - - -
**
**  Convert spherical coordinates to Cartesian.
**
**  Given:
**     theta    double       longitude angle (radians)
**     phi      double       latitude angle (radians)
**
**  Returned:
**     c        double[3]    direction cosines
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double cp;

   cp = cos(phi);
   c[0] = cos(theta) * cp;
   c[1] = sin(theta) * cp;
   c[2] = sin(phi);

   return;

}
// ***************************************************************************

double eraSepp(double a[3], double b[3])
/*
**  - - - - - - - -
**   e r a S e p p
**  - - - - - - - -
**
**  Angular separation between two p-vectors.
**
**  Given:
**     a      double[3]    first p-vector (not necessarily unit length)
**     b      double[3]    second p-vector (not necessarily unit length)
**
**  Returned (function value):
**            double       angular separation (radians, always positive)
**
**  Notes:
**
**  1) If either vector is null, a zero result is returned.
**
**  2) The angular separation is most simply formulated in terms of
**     scalar product.  However, this gives poor accuracy for angles
**     near zero and pi.  The present algorithm uses both cross product
**     and dot product, to deliver full accuracy whatever the size of
**     the angle.
**
**  Called:
**     eraPxp       vector product of two p-vectors
**     eraPm        modulus of p-vector
**     eraPdp       scalar product of two p-vectors
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double axb[3], ss, cs, s;

/* Sine of angle between the vectors, multiplied by the two moduli. */
   eraPxp(a, b, axb);
   ss = eraPm(axb);

/* Cosine of the angle, multiplied by the two moduli. */
   cs = eraPdp(a, b);

/* The angle. */
   s = ((ss != 0.0) || (cs != 0.0)) ? atan2(ss, cs) : 0.0;

   return s;

}
// ****************************************************************************

double eraSeps(double al, double ap, double bl, double bp)
/*
**  - - - - - - - -
**   e r a S e p s
**  - - - - - - - -
**
**  Angular separation between two sets of spherical coordinates.
**
**  Given:
**     al     double       first longitude (radians)
**     ap     double       first latitude (radians)
**     bl     double       second longitude (radians)
**     bp     double       second latitude (radians)
**
**  Returned (function value):
**            double       angular separation (radians)
**
**  Called:
**     eraS2c       spherical coordinates to unit vector
**     eraSepp      angular separation between two p-vectors
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double ac[3], bc[3], s;

/* Spherical to Cartesian. */
   eraS2c(al, ap, ac);
   eraS2c(bl, bp, bc);

/* Angle between the vectors. */
   s = eraSepp(ac, bc);

   return s;

}
// ****************************************************************************

/*
*+
*  Name:
*     palDh2e

*  Purpose:
*     Horizon to equatorial coordinates: Az,El to HA,Dec

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palDh2e( double az, double el, double phi, double * ha, double * dec );

*  Arguments:
*     az = double (Given)
*        Azimuth (radians)
*     el = double (Given)
*        Elevation (radians)
*     phi = double (Given)
*        Observatory latitude (radians)
*     ha = double * (Returned)
*        Hour angle (radians)
*     dec = double * (Returned)
*        Declination (radians)

*  Description:
*     Convert horizon to equatorial coordinates.

*  Authors:
*     PTW: Pat Wallace (STFC)
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - All the arguments are angles in radians.
*     - The sign convention for azimuth is north zero, east +pi/2.
*     - HA is returned in the range +/-pi.  Declination is returned
*       in the range +/-pi/2.
*     - The latitude is (in principle) geodetic.  In critical
*       applications, corrections for polar motion should be applied.
*     - In some applications it will be important to specify the
*       correct type of elevation in order to produce the required
*       type of HA,Dec.  In particular, it may be important to
*       distinguish between the elevation as affected by refraction,
*       which will yield the "observed" HA,Dec, and the elevation
*       in vacuo, which will yield the "topocentric" HA,Dec.  If the
*       effects of diurnal aberration can be neglected, the
*       topocentric HA,Dec may be used as an approximation to the
*       "apparent" HA,Dec.
*     - No range checking of arguments is done.
*     - In applications which involve many such calculations, rather
*       than calling the present routine it will be more efficient to
*       use inline code, having previously computed fixed terms such
*       as sine and cosine of latitude.

*  History:
*     2012-02-08 (TIMJ):
*        Initial version with documentation taken from Fortran SLA
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 1996 Rutherford Appleton Laboratory
*     Copyright (C) 2012 Science and Technology Facilities Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software: you can redistribute it and/or
*     modify it under the terms of the GNU Lesser General Public
*     License as published by the Free Software Foundation, either
*     version 3 of the License, or (at your option) any later
*     version.
*
*     This program is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU Lesser General Public License for more details.
*
*     You should have received a copy of the GNU Lesser General
*     License along with this program.  If not, see
*     <http://www.gnu.org/licenses/>.

*  Bugs:
*     {note_any_bugs_here}
*-
*/

//#include "pal.h"
#include <math.h>

void
palDh2e ( double az, double el, double phi, double *ha, double *dec) {

  double sa;
  double ca;
  double se;
  double ce;
  double sp;
  double cp;

  double x;
  double y;
  double z;
  double r;

  /*  Useful trig functions */
  sa = sin(az);
  ca = cos(az);
  se = sin(el);
  ce = cos(el);
  sp = sin(phi);
  cp = cos(phi);

  /*  HA,Dec as x,y,z */
  x = -ca * ce * sp + se * cp;
  y = -sa * ce;
  z = ca * ce * cp + se * sp;

  /*  To HA,Dec */
  r = sqrt(x * x + y * y);
  if (r == 0.) {
    *ha = 0.;
  } else {
    *ha = atan2(y, x);
  }
  *dec = atan2(z, r);

  return;
}
// ************************************************************************
/*
*+
*  Name:
*     palDtp2s

*  Purpose:
*     Tangent plane to spherical coordinates

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palDtp2s( double xi, double eta, double raz, double decz,
*               double *ra, double *dec);

*  Arguments:
*     xi = double (Given)
*        First rectangular coordinate on tangent plane (radians)
*     eta = double (Given)
*        Second rectangular coordinate on tangent plane (radians)
*     raz = double (Given)
*        RA spherical coordinate of tangent point (radians)
*     decz = double (Given)
*        Dec spherical coordinate of tangent point (radians)
*     ra = double * (Returned)
*        RA spherical coordinate of point to be projected (radians)
*     dec = double * (Returned)
*        Dec spherical coordinate of point to be projected (radians)

*  Description:
*     Transform tangent plane coordinates into spherical.

*  Authors:
*     PTW: Pat Wallace (STFC)
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  History:
*     2012-02-08 (TIMJ):
*        Initial version with documentation taken from Fortran SLA
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 1995 Rutherford Appleton Laboratory
*     Copyright (C) 2012 Science and Technology Facilities Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software: you can redistribute it and/or
*     modify it under the terms of the GNU Lesser General Public
*     License as published by the Free Software Foundation, either
*     version 3 of the License, or (at your option) any later
*     version.
*
*     This program is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU Lesser General Public License for more details.
*
*     You should have received a copy of the GNU Lesser General
*     License along with this program.  If not, see
*     <http://www.gnu.org/licenses/>.

*  Bugs:
*     {note_any_bugs_here}
*-
*/

//#include "pal.h"
//#include "pal1sofa.h"

#include <math.h>

void
palDtp2s ( double xi, double eta, double raz, double decz,
           double *ra, double *dec ) {

  double cdecz;
  double denom;
  double sdecz;
  double d;

  sdecz = sin(decz);
  cdecz = cos(decz);
  denom = cdecz - eta * sdecz;
  d = atan2(xi, denom) + raz;
  *ra = eraAnp(d);
  *dec = atan2(sdecz + eta * cdecz, sqrt(xi * xi + denom * denom));

  return;
}
// ***************************************************************************




int main(int argc, char* argv[])
{
  //std::cout << "-----------------------------------------------------------"
  //             "-------------------------------" << std::endl;
  std::cout << "Phosim Catalog Generator" << std::endl;
  std::cout << "-----------------------------------------------------------"
	           "-------------------------------" << std::endl;

  // First handle the input options.
  CGParams* pCGParams = new CGParams(argc,argv);
  pCGParams->PrintOptionValues();

  //Random*  pRandom = new Random();

  // Need to set random seed  before we use random numbers
  //pRandom->setSeed32( pCGParams->fRandomSeed); 
  
  // **********************************************************************
  // Initalize object IDs to start at zero. While this may duplicate object 
  // ids in other catalogs thats ok.
  // **********************************************************************
  pCGParams->fCatalogObjectID=0;    
  pCGParams->fCatalogObjectCount=0;    


  // *************************************************************************
  // Now we start to produce objects for the catalog.  We may have a combined 
  // catalog of stars, galaixies and a stargrid as chosen by user.
  // *********************************************************************
  // Set up the catalog file name and open the catalog
  if (pCGParams->makeStars() || pCGParams->makeStarGrid() || 
	  pCGParams->makeGalaxies()) {
	

	pCGParams->fCatFileName = "catgen_" + pCGParams->fObservationID +  ".cat";
	pCGParams->fCatFile.open(pCGParams->fCatFileName.c_str(), std::ios::out);
  }
  else{
	std::cout<<" No catalog request found. No catalog made"<< std::endl;
	pCGParams->fCatFile.close();
	exit(1);
  }
  // *************************************************************************
  // First are random stars.All use default SED.
  // *************************************************************************
  if ( pCGParams->makeStars() ) {
	StarObject starObj(pCGParams,pCGParams->fCatFile);
    // *********************
    // Generate stars with galactic latitude density distribution and 
    // magnitude. No longitude distribution as of yet.
    // This method also used the sphericalCell alogrithum for reproducibility.
    // *********************
    starObj.GenerateGalacticStars();
}

  // ***************************
  // Grid of stars. All have same magnitiude and defualt SED
  // ***************************
  if( pCGParams->makeStarGrid() ) {
	StarObject starGridObj(pCGParams,pCGParams->fCatFile);
	starGridObj.GenerateGridStars();
  }

  // ************************************************************************
  // Random Complex Galaxies (in ra/dec cells)
  // ************************************************************************
  if ( pCGParams->makeGalaxies() ) {
	GalaxyObject galaxyObj(pCGParams,pCGParams->fCatFile);
	galaxyObj.GenerateComplexGalaxiesForCells();
  }

  // Display number of objects created (stars + grid + galaxies) in file.
  std::cout<< "Total Number Simulated objects: " 
		   << pCGParams->fCatalogObjectCount << std::endl; 
 
  return 0;
}			  

