/// @brief Grating Class
///
/// @brief Created by:
/// @author Glenn Sembroski (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef GRATING_H
#define GRATING_H

#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "constants.h"

int const kMaxNumIntervals = 54000; // Number of intervals to approximate
double const kStepRad = static_cast<double>(PI)/static_cast<double>(kMaxNumIntervals);
double const kBlazeAngleDeg = 7.09;
int const kNumSlits = 300;
double const kDeltaNm = 3333.33;

class Grating {
// Class to determine diffraction from a diffraction grating

public:
    Grating();
    Grating(double blazeAngleDeg, int numSlits, double deltaNm);
    ~Grating();

    void diffract(double vxIn, double vyIn, double vzIn,
                  double vxGratingNormal, double vyGratingNormal,
                  double vzGratingNormal,
                  double& vxOut, double& vyOut, double& vzOut,
                  double wavelengthNm);

    void setAngleBlazeRad(double angRadians) {
        _angleBlazeRad = angRadians;
        return;
    };
    void setNumSlits(int number) {
        _numSlits = number;
        return;
    };
    void setDeltaNm(double deltaNm) {
        _deltaNm = deltaNm;
        return;
    };

private:
    double _angleBlazeRad;
    double _deltaNm;
    int    _numSlits;
    double _vxGratingNormal;
    double _vyGratingNormal;
    double _vzGratingNormal;
    std::vector < double > _integral;

    //Methods
    double _calculateFunction(double angleInRad, double angleOutRad,
                             double wavelengthNm);
    void _setGratingNormal(double vxGratingNormal, double vyGratingNormal,
                          double vzGratingNormal) {
        _vxGratingNormal = vxGratingNormal;
        _vyGratingNormal = vyGratingNormal;
        _vzGratingNormal = vzGratingNormal;
        return;
    };

    void   _makeTable(double angleInRad, double wavelengthNm);
    int    _binarySearch(double goal);
    double _calculateAngle(double angleInRad, double wavelengthNm);
};

#endif
