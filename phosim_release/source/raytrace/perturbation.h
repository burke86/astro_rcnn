///
/// @package phosim
/// @file perturbation.h
/// @brief header file for perturbation class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef Perturbation_H
#define Perturbation_H

class Perturbation {

 public:
    std::vector<double> eulerPhi;
    std::vector<double> eulerPsi;
    std::vector<double> eulerTheta;
    std::vector<double> decenterX;
    std::vector<double> decenterY;
    std::vector<double> defocus;
    double *zernike_r_grid;
    double *zernike_phi_grid;
    double *zernike_coeff;
    double *zernike_summed;
    double *zernike_summed_nr_p;
    double *zernike_summed_np_r;
    double *miescatter;
    double *rotationmatrix;
    double *inverserotationmatrix;
    double *jitterrot;
    double *jitterele;
    double *jitterazi;
    double *jittertime;
    double *windshake;
    std::vector<int> zernikeflag;
    std::vector<double> rmax;

};

#endif
