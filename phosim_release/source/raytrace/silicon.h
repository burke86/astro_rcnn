///
/// @package phosim
/// @file silicon.h
/// @brief silicon header file
///
/// @brief Created by:
/// @author Andy Rasmussen (SLAC)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "ancillary/random.h"

class Silicon {

public:

    double *meanFreePath;
    float *sigma, *sigmaX, *sigmaY, *fsigma, *gsigma, *hsigma, *isigma;
    float *gammaX, *gammaY;
    float *deltaX, *deltaY;
    float *nbulkmap, *deadLayer, *periodmap;
    double *indexRefraction;
    double *temperatureGrid;
    double *rho;
    double *wavelengthGrid;
    double *thicknessGrid;
    double *dopantGrid;
    double spaceChargeShield;
    double spaceChargeSpreadX;
    double spaceChargeSpreadY;
    double chargeStopCharge;
    double x0;
    double y0;
    long numWavelength;
    long numTemperature;
    long numThickness;
    long numDopant;
    Random random;

    double absorptionCoeffMCT(double lambda, double temperature, double x);
    double absorptionCoeffSi(double lambda, double temperature);
    double indexRefractionMCT (double lambda, double temperature, double x);
    double indexRefractionSi(double lambda);
    void setup(std::string devmaterial, double ccdtemp, double nBulk, double nF, double nB, double sF,
               double sB, double tSi, double overdepBias, std::string instrdir,
               long nampx, long nampy, double pixsize, long seedchip,
               float *impurityX, float *impurityY, int impurityvariation,
               double minwavelength, double maxwavelength);
    double dopeProfile(double z, double nbulk, double nf, double nb, double sf, double sb, double tsi);
    double muMCT(double x, double temperature);
    double muSi (double efield, double temperature, int polarity);
    double epsilonMCT(double x);

};
