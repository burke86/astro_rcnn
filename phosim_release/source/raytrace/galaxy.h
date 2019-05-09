///
/// @package phosim
/// @file galaxy.h
/// @brief galaxy header
///
/// @brief Created by:
/// @author Suzanne Lorenz (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <vector>
#include "ancillary/random.h"

class Galaxy {

public:

    double **fAll;
    double **cAll;
    double *dAll;
    double **fAll2d;
    double **cAll2d;
    double *dAll2d;
    double *clumpX;
    double *clumpY;
    double *clumpZ;
    double *phiArr;
    double *thetaArr;
    double **sphericalHarmonic;
    Random *random;

    void sersic(double a, double b, double c, double alpha, double beta,
                double n, double *xOut, double *yOut, int thread);

    void sersicComplex(double a, double b, double c, double alpha, double beta,
                       double n, double clumpFrac, double clump, double clumpWidth, double spiralFrac,
                       double alphaSpiral, double bar, double spiralWidth, double phi0, int thread, int *init,
                       double *xOut, double *yOut);

    void sersicDiskComplex(double a, double b, double c, double alpha, double beta,
                           double n, double clumpFrac, double clump, double clumpWidth, double spiralFrac,
                           double alphaSpiral, double bar, double spiralWidth, double phi0, int thread, int *init,
                           double *xOut, double *yOut);

    void sampleSersic(char*);

    void sersic2d(double a, double b, double beta, double n,
                  double *xOut, double *yOut, int thread);

    void sersicDisk(double a, double b, double c, double alpha, double beta, double n,
                    double *xOut, double *yOut, int thread);

    void sampleSersic2d(char*, int numthread);

    void distortedSphere(double alpha, double coeff0, double coeff1, double coeff2, double coeff3, double coeff4, double coeff5, double coeff6, double coeff7, double coeff8,
                         double *x_out, double *y_out, int thread);

};
