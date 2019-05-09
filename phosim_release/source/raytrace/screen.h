///
/// @package phosim
/// @file screen.h
/// @brief header file for screen class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef Screen_H
#define Screen_H

#include <fftw3.h>
#include <fitsio.h>
#include <fitsio2.h>
#include "parameters.h"

class Screen {

public:
    double large_sizeperpixel;
    double coarse_sizeperpixel;
    double medium_sizeperpixel;
    double fine_sizeperpixel;
    double *hffunc;
    double *hffunc_n;
    float *turbulenceCoarseX;
    float *turbulenceCoarseY;
    float *turbulenceLargeX;
    float *turbulenceLargeY;
    float *turbulenceMediumX;
    float *turbulenceMediumY;
    float *phaseLarge;
    float *phaseCoarse;
    float *phaseMedium;
    float *phaseFine;
    float *phaseMediumH;
    float *phaseFineH;
    float *cloud[MAX_LAYER];
    float *see_norm, *phase_norm;
    float secondKickSize;
    double *phasescreen;
    double *focalscreen;
    double *tfocalscreen;
    fftw_complex *inscreen;
    fftw_complex *outscreen;
    double wavelengthfactor_nom;
    double *jitterwind;
    double *focalscreencum;
    float *pupil_values;
    double pupilscreenscale;
    double paddingfactor;

    double readScreen(int keynum, float *array, char *tempstring);

};


#endif
