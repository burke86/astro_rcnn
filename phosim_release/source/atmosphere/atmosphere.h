///
/// @package phosim
/// @file atmosphere.h
/// @brief header for atmosphere class
///
/// @brief Created by
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sstream>
#include <vector>
#include <stdexcept>

#include "raytrace/helpers.h"
#include "ancillary/random.h"
#include "ancillary/fft.h"
#include "ancillary/fits.h"
#include "ancillary/readtext.h"

class Atmosphere {

 public:

    int numlevel;
    float groundlevel;
    double constrainclouds;
    std::string datadir;
    std::string instrdir;
    std::vector<float> osests;
    std::vector<float> altitudes;
    std::vector<float> jests;
    Random random;

    void createAtmosphere(float monthnum, float constrainseeing, const std::string & outputfilename, const std::vector<int> & cloudscreen, long seed, double tai, double *tseeing);

    void turb2d(long seed, double see5, double outerx, double outers,
            double zenith, double wavelength, const std::string & name,
            long N_size=1024);

    void cloud(long seed, double cloheight, double pixsz,
               const std::string & name, long N_size = 1024);

    void airglow(long seed, const std::string & name, long screenSize=1024);

    void ccalc(double tai, int call);

    void magcalc(float monthnum, float altitude, float & magest, float & direst, double r1, double r2);
    float outerscale(float altitude, double tai, int n);

};
