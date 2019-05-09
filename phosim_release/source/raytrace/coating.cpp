///
/// @package phosim
/// @file coating.cpp
/// @brief coating class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "coating.h"

void Coating::setup (long totalSurfaces) {

    wavelengthNumber = new long[totalSurfaces]();
    angleNumber = new long[totalSurfaces]();

    reflection = new double*[totalSurfaces]();
    transmission = new double*[totalSurfaces]();
    wavelength = new double*[totalSurfaces]();
    angle = new double*[totalSurfaces]();

}

void Coating::allocate (long surfaceIndex, long lines) {

    transmission[surfaceIndex] = new double[lines]();
    reflection[surfaceIndex] = new double[lines]();
    wavelength[surfaceIndex] = new double[lines]();
    angle[surfaceIndex] = new double[lines]();

}
