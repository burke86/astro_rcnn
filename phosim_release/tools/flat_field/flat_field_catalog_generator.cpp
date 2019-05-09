///
/// @package phosim
/// @file ISphereCatalogGenerator
/// @brief Generates phosim Object catalog for integrating sphere simulation
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
#include <vector>
#include <cstring>
#include <cstdlib>
#include <iomanip>

using namespace std;
// ***********************************************************************
// Following taken from sla library)
const double kPi  (3.1415926535897932384626433832795028841971693993751);
const double k2Pi (6.2831853071795864769252867665590057683943387987502);
/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                                                         :(A)+(B)*floor(-(A)/(B))):(A))
/* dsign(A,B) - magnitude of A with sign of B (double) */
#define dsign(A,B) ((B)<0.0?-(A):(A))


void AzEl2RADec(double fAz, double fEl, double& fRAApparent,
		double& fDecApparent);
void slaDh2e ( double az, double el, double phi, double *ha, double *dec );
double slaDranrm ( double angle );
double slaDrange ( double angle );

class FibinochiSpiralLattice
{
public:
  FibinochiSpiralLattice(){};
  ~FibinochiSpiralLattice(){};
  void genFibinochiLattice(int numTotalPointsSphere, double maxZnDeg,
			   vector <double>& ZnDeg, vector <double>& AzDeg);
};
    
void FibinochiSpiralLattice::genFibinochiLattice(int numTotalPointsSphere, 
						 double maxZnDeg,
						 vector <double>& ZnDeg, 
						 vector <double>& AzDeg)
{
  // ********************************************************************
  // from:
  //  "Measurement of areas on a sphere using Fibonacci and latitude-longitude 
  //  lattices", Alvaro Gonz'alez', arXiv:0912.4540v1 and Mathematical 
  //  Geosciences, in press.
  // *********************************************************************
  int number=(numTotalPointsSphere-1)/2;
  //std::cout<<"Number: "<<number<<std::endl;
  double num=(double)number;


  double phi =0;
  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;  //golden ratio
  //std::cout<<"phi(Golden ratio): "<<phi<<std::endl;
  //std::cout<<"num: "<<num<<std::endl;

  ZnDeg.clear();
  AzDeg.clear();

  int icount=0;
  for(int i=-number;i<=number ;i++){
    double point=(double)i;
    double latDeg=asin(2.*point/(2.*num+1.))*180./kPi;

    //Mod divide of i/phi
    int result = static_cast<int>(point/phi);
    double mod = point-static_cast<double>(result)*phi;
    double lonDeg=mod*360.0/phi;   

    if(lonDeg<-180){
      lonDeg=lonDeg+360.;
    }
    if (lonDeg>180){
      lonDeg=lonDeg-360.;
    }
    double zenithDeg=90.-latDeg;
    if (zenithDeg<maxZnDeg && zenithDeg>=0.0){
      ZnDeg.push_back(zenithDeg);
      //cout<<lonDeg<<endl;
      if(lonDeg<0){
	lonDeg=lonDeg+360.0;
      }
      AzDeg.push_back(lonDeg);
      // std::cout<<i<<" "<<zenithDeg<<" "<<lonDeg<<std::endl;
      icount++;
    }
  }
  std::cout<<" Number of points within "<<maxZnDeg<<" deg cap: "<<icount<<endl;
  return;
}
// ****************************************************************************


// ***********************************************************************
// Generates a catalog of a matrix of overlapping gaussian objects in the sky 
// whose combined light will will produce a distribution on an aperature with 
// an even lambertian distribution, simulating the light intensity and 
// distribution from the exit aperature of an integrating sphere.
// ***********************************************************************
// May be expanded later to simulate outputs with a gradiant etc by varing 
// relative strengths of the source object matrix. For now its all 1.0
// ***********************************************************************
//
//Use:  ./ISphereCatalogGenerator -h option to see a list of input parameters.
//
void usage()
{
  cout<<" ISphereCatalogGenerator usage:"<<endl;
  cout<<"-h : Help (usage)"<<endl;
  cout<<"-o : Output Object Catalog file name (Defaults to ISphere.cat)"<<endl;
  cout<<"-s : Sed File Name(as needed by phosim, with path!)  (Mandatory)"<<endl;
  cout<<"-M : Combined Source magnitude (defaults to 20.0)"<<endl;
  cout<<"-t : Max angle of light in Deg (defaults to ~45 Deg) "<<endl;
  // following are for testing usually
  cout<<"-N : Number of points in spherical cap (defaults to 10,000)"<<endl;
  cout<<"-G : sigma of gausian sources(arcsec)(defaults to 5 times spacing)"
      <<endl;
  exit(0);
}

int main(int argc, char** argv)
{
  // *********************************************************************
  //Read in all parameters and open files as needed
  // *********************************************************************
  string objectCatalogName="ISphere.cat"; 
  string sedFileName;
  bool   gotSed=false;
  double sourceMagnitude=20;    //defulat: combined object effective magnitude
  double maxZnDeg=45.0;         //Default: Approximate maximum theta of light 
                                //emitted:Limits source placement from center
  double sigmaSourceArcSec;     //Gaussion width of grid sources
  bool gotSigma=false;
  int    numSphericalCap=10000; //Number of points to aim for in spherical cap
  // *****************************************************
  // get input params
  // *****************************************************
  for (int i = 1; i < argc; i++) {

    if (strcmp(argv[i],"-h")==0) {
      usage();
    } 
    else if (strcmp(argv[i],"-o")==0) {
      objectCatalogName=argv[i+1];
      cout<<"Object Catalog Name: "<<objectCatalogName<<endl;
    } 
    else if (strcmp(argv[i],"-s")==0) {
      sedFileName=argv[i+1];
      gotSed=true;
      cout<<" Sed File Name: "<<sedFileName<<endl;
    } 
    else if (strcmp(argv[i],"-M")==0) {
      sourceMagnitude = atof(argv[i + 1]);
      cout<<" Effective Source Magnitude: "<< sourceMagnitude<<endl;
    }
    else if (strcmp(argv[i],"-t")==0) {
      maxZnDeg = atof(argv[i + 1]);
      cout<<" Max Zenith of emitted photons(Deg): "<<maxZnDeg<<endl;
    }
    else if (strcmp(argv[i],"-N")==0) {
      numSphericalCap = atoi(argv[i + 1]);
      cout<<" Number of points in spherical Cap: "<<numSphericalCap<<endl;
    }
    else if (strcmp(argv[i],"-G")==0) {
      sigmaSourceArcSec = atof(argv[i + 1]);
      cout<<" Sigma of gaussian sources: "<<sigmaSourceArcSec<<endl;
      gotSigma=true;
    }
  }
  
  // ****************************************************************
  // Check mandatory got set
  // ****************************************************************
  if (!gotSed) {
    cout<<" Fatal--Missing Sed file name option (-s) Its mandatory!"<<endl;
    cout<<" Use -h for usage help"<<endl;
    exit(1);
  }
  
  // *******************************************************************
  // Create output object catalog file
  // *******************************************************************
  ofstream os(objectCatalogName.c_str());

  // ********************************************************************
  // Idea is to generate a grid of gaussian objects that overlap  
  // the space of the requested spherical cap.
  // We use a Fibinochi spiral matrix to space the points
  // *********************************************************************
  // Find all the Ra/Dec locations of the gaussian sources. The center will be 
  // at Ra=0.0 and Dec=0.0
  // ********************************************************************
  // See functions below for references to the source of the fibinochi 
  // distribution of points
  // *********************************************************************
  // Find approx number of points we need on the sphere to give us the 
  // requested number in our Spherical cap.  
  // Below all sets R==1.0 Example=>area sphere=4*Pi*(R^2)=4*Pi
  // See: http://en.wikipedia.org/wiki/Spherical_cap
  // *********************************************************************
  double R=1.0;  //Just included for clarity
  double sphCapH=R*(1.0-cos(maxZnDeg*kPi/180.));
  double sphCapArea=2*kPi*R*sphCapH;
  double AreaRatio=sphCapArea/(4*kPi*R*R);
  int  numTotalPointsHalfSphere= (numSphericalCap/AreaRatio)/2;
  // We want numTotalPointsSphere odd
  int  numTotalPointsSphere=2*numTotalPointsHalfSphere;
  numTotalPointsSphere=numTotalPointsSphere+1;
  cout<<numTotalPointsSphere<<"  "<<numTotalPointsHalfSphere<<endl;

  vector <double> ZnDeg;
  vector <double> AzDeg;

  // ***********************************************
  // Generate fibiunochi spherical lattic
  // ***********************************************
  FibinochiSpiralLattice FibSpiral;
  FibSpiral.genFibinochiLattice(numTotalPointsSphere,maxZnDeg,ZnDeg,AzDeg);
  
  int numPointsCap=ZnDeg.size();

  if(!gotSigma){
    // ********************************************************************
    // Find average spacing so that we can find approximate sigma to use
    // We do this by finding number of points in our set that it takes
    // to girdle the equator of the sphere.
    // (I keep R in these equations for clairity)
    // ********************************************************************
    double areaCell=4.*kPi*R*R/ numTotalPointsSphere;
    //Area each point represents
    double cellR=sqrt(areaCell/kPi);// This is r
    double approxNumCellsEquator=2*kPi*R*R/cellR;  //Number cells around
    //equator
    double approximateSpacingArcSec=360.*60.*60./approxNumCellsEquator;
    sigmaSourceArcSec=5*approximateSpacingArcSec/2.; //5 is arbitray.
    cout<<sigmaSourceArcSec<<" "<<approximateSpacingArcSec<<" "
	<<approxNumCellsEquator<<" "<<cellR<<" "<<areaCell<<endl;
  }
  
  // **********************************************************************
  // Find the magnitude each source needs so thatt all the points sum up to 
  // the requested magnitude. Remember magnitudes are in logs of flux.
  // **********************************************************************
  double pointMagnitude=sourceMagnitude+2.5*log10(numPointsCap);

  // *********************************************************************
  // Set time to MJD2000.0 (I think, we may want this as a setrtable input
  // parameter some day 
  //  Find RA,Dec of center of cap
  // *********************************************************************
  double DecRad;
  double RaRad;
  AzEl2RADec(0.0,kPi/2.0, RaRad, DecRad);
  cout<<"Direction of center of cap (should be tracking direction!)RA: "
      <<RaRad*180./kPi<<"(deg.)   Dec: "<<DecRad*180./kPi<<endl;



  // ********************************************************************
  // Now we generate the object catalog
  // ********************************************************************
for (int i=0;i<numPointsCap; i++){
    // **************************************************************
    // Convert Zn/Az of this point to an RA/Dec.
    // Assume Ra=0,Dec=0 is at zenith (ZN=0)
    // **************************************************************
    AzEl2RADec(AzDeg.at(i)*kPi/180.,kPi/2.0-ZnDeg.at(i)*kPi/180, RaRad, DecRad);
    //std::cout<<AzDeg.at(i)<<" "<<ZnDeg.at(i)<<" "<<RaRad*180./kPi<<" "
    //	     <<DecRad*180.0/kPi<<std::endl;


    
    // ***************************************************************
    // Add this point as an object to the catalog file. Make it a gaussian 
    // source with a width several times (I use 5) the half spacing spacing
    // ****************************************************************
    os<<"object "<<i<<" "<<std::fixed<<std::setprecision(13)<<RaRad*180./kPi
      <<" "<<DecRad*180.0/kPi<<"  "<<std::setprecision(4)<<pointMagnitude
      <<" "<<sedFileName<<" 0.0 0.0 0.0 0.0 0.0 0.0 gaussian "
      <<sigmaSourceArcSec<<" none  none"<<endl;
  }
  

return 0;
}
// *************************************************************************

void AzEl2RADec(double fAz, double fEl, double& fRA, 
		double& fDec)
// ************************************************************************
// Convert Az and El to RaDec epoch 2000 where AZ= 0 => RA=0.0
// ************************************************************************
{
  //double kLatitude               = 0.552828;
  //double kEastLongitude          = -1.93649;
  double kLatitude               = 0.0;
  double fHourangle;
  
  slaDh2e(fAz, fEl, kLatitude, &fHourangle, &fDec);

  fRA = - fHourangle;
       
  fRA = slaDranrm(fRA);
  fDec = slaDrange(fDec);
  return;
}
// ************************************************************************

// *********************************************************************
// Following from sla Library
// *********************************************************************
void slaDh2e ( double az, double el, double phi, double *ha, double *dec )
  /*
  **  - - - - - - - -
  **   s l a D h 2 e
  **  - - - - - - - -
  **
  **  Horizon to equatorial coordinates:  Az,El to HA,Dec
  **
  **  (double precision)
  **
  **  Given:
  **     az          double       azimuth
  **     el          double       elevation
  **     phi         double       observatory latitude
  **
  **  Returned:
  **     *ha         double       hour angle
  **     *dec        double       declination
  **
  **  Notes:
  **
  **  1)  All the arguments are angles in radians.
  **
  **  2)  The sign convention for azimuth is north zero, east +pi/2.
  **
  **  3)  HA is returned in the range +/-pi.  Declination is returned
  **      in the range +/-pi/2.
  **
  **  4)  The is latitude is (in principle) geodetic.  In critical
  **      applications, corrections for polar motion should be applied.
  **
  **  5)  In some applications it will be important to specify the
  **      correct type of elevation in order to produce the required
  **      type of HA,Dec.  In particular, it may be important to
  **      distinguish between the elevation as affected by refraction,
  **      which will yield the "observed" HA,Dec, and the elevation
  **      in vacuo, which will yield the "topocentric" HA,Dec.  If the
  **      effects of diurnal aberration can be neglected, the
  **      topocentric HA,Dec may be used as an approximation to the
  **      "apparent" HA,Dec.
  **
  **  6)  No range checking of arguments is done.
  **
  **  7)  In applications which involve many such calculations, rather
  **      than calling the present routine it will be more efficient to
  **      use inline code, having previously computed fixed terms such
  **      as sine and cosine of latitude.
  **
  **  Last revision:   30 November 2000
  **
  **  Copyright P.T.Wallace.  All rights reserved.
  */
{
  double sa, ca, se, ce, sp, cp, x, y, z, r;
  
  /* Useful trig functions */
  sa = sin ( az );
  ca = cos ( az );
  se = sin ( el );
  ce = cos ( el );
  sp = sin ( phi );
  cp = cos ( phi );
  
  /* HA,Dec as x,y,z */
  x = - ca * ce * sp + se * cp;
  y = - sa * ce;
  z = ca * ce * cp + se * sp;
  
  /* To spherical */
  r = sqrt ( x * x + y * y );
  *ha = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
  *dec = atan2 ( z, r );
  return;
}
// **********************************************************************


double slaDranrm ( double angle )
  /*
  **  - - - - - - - - - -
  **   s l a D r a n r m
  **  - - - - - - - - - -
  **
  **  Normalize angle into range 0-2 pi.
  **
  **  (double precision)
  **
  **  Given:
  **     angle     double      the angle in radians
  **
  **  The result is angle expressed in the range 0-2 pi (double).
  **
  **  Defined in slamac.h:  k2Pi, dmod
  **
  **  Last revision:   19 March 1996
  **
  **  Copyright P.T.Wallace.  All rights reserved.
  */
{
  double w;
  
  w = dmod ( angle, k2Pi );
  return ( w >= 0.0 ) ? w : w + k2Pi;
}
// ********************************************************************

double slaDrange ( double angle )
/*
**  - - - - - - - - - -
**   s l a D r a n g e
**  - - - - - - - - - -
**
**  Normalize angle into range +/- pi.
**
**  (double precision)
**
**  Given:
**     angle     double      the angle in radians
**
**  The result is angle expressed in the +/- pi (double precision).
**
**  Defined in slamac.h:  kPi, k2Pi, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
  double w;

  w = dmod ( angle, k2Pi );
  return ( fabs ( w ) < kPi ) ? w : w - dsign ( k2Pi, angle );
}

