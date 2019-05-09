///
/// @package phosim
/// @file GalaxyComplex.cpp
/// @brief Generates phosim sersicComplex galaxy parameters randomly (but 
///   repeatably) using IDL code provided by John Peterson converted to C++
// **************************************************************
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
/// Original algorithm from IDL code provided by John Peterson.

#include "galaxycomplex.h"


// *********************************************************************
GalaxyComplex::GalaxyComplex()
{
  Init();
}

void GalaxyComplex::Init()
{
  // MagnitueValues and zValues are lookup tables for more efficency.
  fMabsoluteValues.resize( kNumValues, 0.0);
  fCumulativeMabsoluteProb.resize( kNumValues, 0.0);

  std:: vector < double > densityMabsolute;
  densityMabsolute.resize( kNumValues, 0.0);

  std:: vector < double > MabsoluteProb;
  MabsoluteProb.resize( kNumValues, 0.0);


  double mmin = -17.0;   //This is an lowest absolute magnitude we are 
                         //interested in
  double densityMabsoluteTotal = 0;
  for(int i = 0; i < kNumValues; i++) {
    // Create the Absolute Magnitude array. This has as M values starting at 
    // M=-25.  Has steps of .00001 M. Last value is M~=-17. Milky way is -20.8
    // Mabsolute. A look up table, sort of.
	fMabsoluteValues.at(i) = -25.0 + ((double)i / (double) kNumValues ) *
	                                                          (mmin + 25.0);
	
	// determine densityMabsolute.  It is distributed as a function of Mabsolute
    double ppower =  0.4 * ( -20.3 - fMabsoluteValues.at(i));
	densityMabsolute.at(i)= 0.4 * log(10.0) * 0.005 * 
			 pow(10., (ppower * ( -1.2 + 1.0))  ) * exp(- pow(10.,  ppower )) ;
    //Total density over Mabsolute range.
    densityMabsoluteTotal += densityMabsolute.at(i); 
  }

  // Cumulative Mabsolute Probabilty distribution
  // **********************************************************
  // Normalize densityMabsolute. Makes it into a probability. Create 
  // cumulative densityMabsolute prob. This will be used to find a 
  // galaxy with a random Mabsolute. 
  // **********************************************************
  fCumulativeMabsoluteProb.at(0) = MabsoluteProb.at(0);  //Initial value
  for( int i = 1; i < kNumValues; i++) { 
	MabsoluteProb.at(i) = densityMabsolute.at(i)/densityMabsoluteTotal;
    fCumulativeMabsoluteProb.at(i) = fCumulativeMabsoluteProb.at(i-1) + 
                                                             MabsoluteProb.at(i);
  }
  // ************************************************************

  //Redshift Distribution
  fZValues.resize(kNumValues,0.0);
  fCumulativeZProb.resize(kNumValues,0.0);

  std::vector < double > zProb;
  zProb.resize(kNumValues,0.0);


  double zc = 1.0;  //Some type of z scale?

  double totalZProb = 0.0;
  for (int i=0; i<kNumValues; i++) {
	//z=5 is max z we are interested in.
	fZValues.at(i) = ((double)i/(double)kNumValues)*5.0;
	zProb.at(i) = pow( fZValues.at(i), 2) * exp(-fZValues.at(i) / zc); 
	totalZProb += zProb.at(i);
  }

  //Cumulative Redshift Probability Distribution
  fCumulativeZProb.at(0) = zProb.at(0);
  for (int i=1; i< kNumValues ; i ++ ) {
	zProb.at(i) = zProb.at(i) /	totalZProb;
    fCumulativeZProb.at(i) = fCumulativeZProb.at(i-1) + zProb.at(i);
  }

  // Density of galaxies
  // ***********************************************************************
  // This version is for the spherical cell algorithum (cells may be of 
  // differeing areas ( in deg**2). So leave inclusion of areas to later)
  // To get a number in a cell ee will apply area in deg2 later. So now add 
  // conversions to get denGalaxyDeg2 the correct units(the 
  // "pow(60., 2)" deg**2 to arcmnin**2 term in denGalaxyDeg2).
  // Don't forget John's fudge factor which he says we will address later.
  // Note: We use a double for the densityMabsolute here. We calculate 
  // complete number later (we will mutilpy by an area in deg2 and fluxtuate 
  // appropriatly (poisson).).
  // ***********************************************************************
  // Breaking up density calculation into terms for computation clarity.
  // But where the breaks should be physics wise (ie what terms mean what) 
  // is unknown to me

  double term1 = ( mmin + 25.0 ) / kNumValues;
  double term2 = 4. * kPI * pow( (3e5 / 70.0), 3) ;
  double term3 = ( ( 4. * kPI) / pow( (kPI/180.), 2) )  * 60.0 * 60.0;
  double denGalaxyDeg2 = densityMabsoluteTotal * term1 * (term2 / term3) * 
                                      pow(60., 2);
  
  // I think  this is densityMabsolute of galaxies for the complete possible 
  // Mabsolute range!
  fGalaxyDensityDeg2 = ( denGalaxyDeg2 * (totalZProb/kNumValues) * 5.0 );  

  return;
}
// ***************************************************************************

std::string GalaxyComplex::GenerateSersicComplexObjectString(
												 GalacticComplexObject& Object)
{
  // **********************************************************************
  // Generate string that goes into phosim object catalog for the 
  // sersicComplex part
  // **********************************************************************
  std::ostringstream os;
  os << "sersicComplex" 
	 << std::fixed << std::setw(8) << std::setprecision(4)
	 << Object.a << " " << Object.b << " " << Object.c << " "
	 << Object.alpha << " " << Object.beta << " " << Object.index << " "
	 << Object.clumpFrac << " " << Object.clump << " " 
	 << Object.clumpWidth << " " << Object.spiralFrac << " " 
	 << Object.alphaSpiral << " " << Object.bar << " "
	 << Object.spiralWidth << " " << Object.phi0 << " ";

  return os.str();
}
// ************************************************************************

void GalaxyComplex::WriteObject( GalacticComplexObject& Object, 
								 std::ofstream* ofs)
// ************************************************************************
// Write the object line to the catalog file
// ************************************************************************
{
  // *********************************
  // Get the sersic part of the line
  // *********************************
  std::string sersicCmplxStr=GenerateSersicComplexObjectString( Object);

   *ofs<< "object " << Object.ID << " " << std::fixed 
       << std::setprecision(8) << Object.RADeg  << " " <<Object.DecDeg
       << " " << std::setprecision(2) << Object.M << " " 
       << Object.SEDFileName << " " << std::setprecision(3) << Object.z 
       << " 0.0 0.0 0.0 0.0 0.0 " << sersicCmplxStr << "none none " 
       << std::endl;
   return;
}
// *************************************************************************

bool GalaxyComplex::GenerateGalaxyForCells(GalacticComplexObject& Bulge, 
                                           GalacticComplexObject& Disk,
                                           double MinimumM, double MaximumM,
                                           Uint64 galaxySeed)
// ****************************************************************************
// This method generates a Mapparent for a galaxy and if within range of 
// Mminimum -Mmaximum it creates the Magnitudes, sizes and sersicComplex  
// paramteres for a random galaxy disk and bulge following simple 
// distributions of such.
// ****************************************************************************
// Note: a  random sequence seed is used here for cell/galaxy 
// repeatability. This seed in cell and galaxy index dependent
// ****************************************************************************
// GalacticComplexObject defined in GalaxyObject class
// ****************************************************************************
// Original algorithm from IDL code provided by John Peterson.
{
  Random random;
  random.setSeed64(galaxySeed);
   
  // ***************************************************************
  // To determin the Mapparent of a galaxy with this methiod we need at the 
  // minimum Mabsolute and the distance.
  // **************************************************************
  // For distance first pick a z from our z distribution.
  // Find index to z in the z arrays. Find index where the 
  // cumulatrive probability exceeds the random selection  
  
  double u = random.uniform();

  // To find the bin in a Cumulative probability distribution vector
  // (naturarly sorted with first bin >= 0 (and maybe a few more as 0) and with
  // the last bin < 1.0) which holds the random 
  // probability u we will use the C++ standard library algorithum 
  // "std::upper_bound". upper_bound does a binary search:
  // upper_bound(ForwardIterator first, ForwardIterator last, const T& val);
  // "Returns an iterator pointing to the first element in the range 
  // (first,last) which compares greater than val" ie. its the first bin who's 
  // value is   > val. This works for all values of u even u == 0.0.
  // Note: We start range at bin 1 not 0. Why? Think about it and you will see.

  fVectorIteratorUpper = std::upper_bound(fCumulativeZProb.begin()+1, 
                                          fCumulativeZProb.end(), u);
  //Find bin index for this value of u. Back up one bin
  int g = fVectorIteratorUpper - fCumulativeZProb.begin() - 1; 
  
  // Thus:
  double redShift = fZValues.at(g);     //redShift in range: 0 -> 5.0
  double distance=3e5/(70.0)*redShift;  //Units?


  // Now pick a random Mabsolute magnitude from its distribution
  u = random.uniform();

  fVectorIteratorUpper = std::upper_bound(fCumulativeMabsoluteProb.begin()+1, 
                                          fCumulativeMabsoluteProb.end(), u);
  //back up 1 bin
  g = fVectorIteratorUpper - fCumulativeMabsoluteProb.begin() - 1;

  //Get absolute magnitude of galaxy
  double MAbsolute = fMabsoluteValues.at(g);

  // Correct for distance. Generate total Apparent magnitude
  double MApparent = MAbsolute + (5.0 * log10( distance * 1e6/10.0 ) ); 

  //Now is the soonest we can check that this MApparent is within the range
  // we specified

  if (MApparent >MaximumM || MApparent < MinimumM) {
    return false;
  }
  else {
    //Debug print to check M distribution
    //std::cout << redShift << " " << distance << "  " << MApparent << " " 
    //          << MAbsolute << " " << Bulge.RADeg << " " << Bulge.DecDeg 
    //          << std::endl; 

    //OK its a good one! Make the bulge and disk objects!
    // Determine total flux and then convert to disk and bulge magnitude
    // Now the angular size of the galaxy
    double ssize =  pow( ( pow( 10.0, (-0.4*(MAbsolute+21.0)))  ), .333) *0.01 /
      ( distance * (kPI /(180.*3600.0) ) );  //Units? (arcsec?)
    
    
    //Divide up MApparent between disk and bulge randomly
    double bulgefrac =  random.uniform();
    
    Bulge.M = -2.5 * log10(    bulgefrac   * pow(10., (-0.4*MApparent)) );
    Disk.M  = -2.5 * log10( (1.-bulgefrac) * pow(10., (-0.4*MApparent)) );
    Bulge.z = redShift;
    Disk.z  = redShift;
    
    
    // Now fill in the sersic values for sersicComplex.
    Bulge.a=ssize*(0.8+0.2*random.uniform())*0.75;
    Bulge.b=ssize*(0.8+0.2*random.uniform())*0.75;
    Bulge.c=ssize*(0.8+0.2*random.uniform())*0.75;
    Bulge.alpha = 0.0;
    Bulge.beta  = 0.0;
    Bulge.index = 4;
    Bulge.clumpFrac   = random.uniform() * 0.5; 
    Bulge.clump       = 300;
    Bulge.clumpWidth  = ssize*0.1;
    Bulge.spiralFrac  = 0;
    Bulge.alphaSpiral = 0;
    Bulge.bar         = 0;
    Bulge.spiralWidth = 0;
    Bulge.phi0        = 0;
    
    Disk.a=ssize*(0.8+0.2*random.uniform())*1.25;
    Disk.b=ssize*(0.8+0.2*random.uniform())*1.25;
    Disk.c=ssize*(0.2*random.uniform())*1.25;
    Disk.alpha = random.uniform()*360.0;
    Disk.beta  = ( acos( 2.0 * random.uniform() - 1.0) 
                   - (kPI/2.0) ) * 180. / kPI;
    Disk.index = 1;
    Disk.clumpFrac   = random.uniform() * 0.5; 
    Disk.clump       = 300;
    Disk.clumpWidth  = ssize*0.1;
    Disk.spiralFrac  = random.uniform();
    Disk.alphaSpiral = random.uniform() * 60.0 - 30.0;
    Disk.bar         = (ssize*0.5);
    Disk.spiralWidth = (ssize*0.1);
    Disk.phi0        = random.uniform() * 360.0 ;
    
    return true; //return total Magnitude of galaxy. Used to slect M range
  }
}
// ************************************************************************

