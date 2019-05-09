///
/// @package phosim
/// @file air.h
/// @brief header file for air class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

class Air {

 public:
    double *tauWavelength;
    double **tau;
    double *airmassLayer;
    double **tauMoon;
    double *airmassLayerMoon;
    double air_refraction_adc;

    void opacitySetup(double zenith, double moonalt, std::vector<double> height, double groundlevel, double raynorm, double o2norm, double h2onorm, double o3norm, double aerosoltau, double aerosolindex, long layers, std::string dir, double *airmass);

};
