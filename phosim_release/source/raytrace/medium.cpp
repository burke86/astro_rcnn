///
/// @package phosim
/// @file medium.cpp
/// @brief medium class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "raytrace/parameters.h"
#include "ancillary/readtext.h"
#include "raytrace/medium.h"

void Medium::setup (int surfaceIndex, std::string mediumFileName) {

    readText mediumPars(mediumFileName);
    size_t nline = mediumPars.getSize();

    indexRefractionNumber[surfaceIndex] = static_cast<long>(nline);
    indexRefraction[surfaceIndex] = static_cast<double*>(calloc(nline, sizeof(double)));
    indexRefractionWavelength[surfaceIndex] = static_cast<double*>(calloc(nline, sizeof(double)));

    for (size_t tt(0); tt < nline; tt++) {
        std::istringstream isst(mediumPars[tt]);
        isst >> *(indexRefractionWavelength[surfaceIndex] + tt);
        isst >> *(indexRefraction[surfaceIndex] + tt);
    }

}
