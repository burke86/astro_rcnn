///
/// @package phosim
/// @file dust.cpp
/// @brief dust absorption
///
/// @brief Created by:
/// @author James Pizagno (UW)
///
/// @brief Modified by:
/// @author John Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "dust.h"
#include "helpers.h"

void Dust::setup (double maxwavelength) {

    int i;
    double a, b;

    wavelengthGrid = static_cast<double*>(calloc(maxwavelength, sizeof(double)));
    aGrid = static_cast<double*>(calloc(maxwavelength, sizeof(double)));
    bGrid = static_cast<double*>(calloc(maxwavelength, sizeof(double)));
    kGrid = static_cast<double*>(calloc(maxwavelength, sizeof(double)));
    for (i = 0; i < maxwavelength; i++) {
        wavelengthGrid[i] = (static_cast<double>(i))/1000.0;
        ccmSetup(wavelengthGrid[i], &a, &b);
        aGrid[i] = a;
        bGrid[i] = b;
        kGrid[i] = calzettiSetup(wavelengthGrid[i]);
    }

}


double Dust::ccm (double wavelength, double maxwavelength, double av, double rv) {

    long index;
    double a, b;
    double rindex;

    index = find_linear(wavelengthGrid, maxwavelength, wavelength, &rindex);
    a = interpolate_linear(aGrid, index, rindex);
    b = interpolate_linear(bGrid, index, rindex);
    return pow(10.0, -0.4*av*(a + b/rv));

}

double Dust::calzetti (double wavelength, double maxwavelength, double av, double rv) {

    long index;
    double rindex;
    double k;

    index = find_linear(wavelengthGrid, maxwavelength, wavelength, &rindex);
    k = interpolate_linear(kGrid, index, rindex);
    return pow(10.0, -0.4*av*(1 + k/rv));

}

// Cardelli, Clayton, Mathis Model
//
// Created 15 Jan 2010
// I changed "av" to A_V for clarity.
// This is the Cardelli, Clayton, Mathis 1989 ApJ, 345, 245 (CCM89) extinction curve.
// It returns the extinction in magnitudes defined by:
// equation1-CCM:  A(l)/A_V = a(x) + b(x)/R_V,
// where x=1/l
//   and av = V-band extinction in magnitudes.
//   and R_V  = av / E(B-V), is the selective extinction.
// R_V provides the shape of the extinction curve, see figs 3,4 in CCM89.
// This is for the Milky Way extinction curve.  There is a different
// curve for the Magellanic Clouds.
//
// Input: R_V = 3.1 for Diffuse ISM,  R_V=5.2 for Molecular Clouds,
//        av = V-band extinction (magnitudes),
//      and l=wavelength (microns)

void Dust::ccmSetup(double l, double *a, double *b) {

    double x, y;
    double fa, fb;

    if (l != 0.0) {
        x = 1.f/l;
    } else {
        x = 1000.0;
    }
    y = x - 1.82f;
    if (x >= 8.0f) {      // the Far-UV
        *a = (-1.073 - 0.628*(x - 8.0) + 0.137*pow(x - 8.0, 2) - 0.070*pow(x - 8.0, 3));
        *b = (13.670 + 4.257*(x - 8.0) - 0.420*pow(x - 8.0, 2) + 0.374*pow(x - 8.0, 3));
    } else if (x >= 3.3 && x < 8.0) {      //  UV
        if (x >= 5.9 && x <= 8.0) {
            fa = (-0.04473*pow(x - 5.9, 2) - 0.009779*pow(x - 5.9, 3));
            fb = (0.2130*pow(x - 5.9, 2) + 0.1207*pow(x - 5.9, 3));
        } else {
            fa = 0.0;
            fb = 0.0;
        }
        *a = ( 1.752 - 0.316*x - 0.104 / (pow(x - 4.67, 2) + 0.341) + fa);
        *b = (-3.090 + 1.825*x + 1.206 / (pow(x - 4.62, 2) + 0.263) + fb);
    } else if (x <= 1.1) {    //  the infrared
        *a = (0.574*pow(x, 1.61));
        *b = (-0.527*pow(x, 1.61));
    } else if (x > 1.1 && x < 3.3) {      // Optical/NIR
        *a = (1.0 + 0.17699*y - 0.50447*pow(y, 2) - 0.02427*pow(y, 3) + 0.72085*pow(y, 4) + 0.01979*pow(y, 5));
        *a = (*a - 0.77530*pow(y, 6) + 0.32999*pow(y, 7));
        *b = (1.41338*y + 2.28305*pow(y, 2) + 1.07233*pow(y, 3) - 5.38434*pow(y, 4) - 0.62251*pow(y, 5));
        *b = (*b + 5.30260*pow(y, 6) - 2.09002*pow(y, 7));
    }
}

// Calzetti et al. Model
//
// written by Jim Pizagno, 7 August 2009
// from Calzetti et al. 2000 ApJ, 533, 683
// equations 2-4.
// call calext(microns, R_V, A_V).
// returns exponential opacity defined by:
//    I_observed = I_true * exp(-tau).
// l=wavelength in microns
// RV = AV/E(B-V) typically 3.1 to 5.0
// AV = V-band extinction

double Dust::calzettiSetup(double l) {

    double k = 100.0;

    if (l < 0.12) {
        k = 100.0;
    } else if (l >= 0.12 && l<0.63) {
        k = (2.659*(-2.156 + 1.509/l - 0.198/pow(l, 2) + 0.011/pow(l, 3)));
    } else if (l>=0.63 && l < 2.2) {
        k = (2.659*(-1.857 + 1.04/l));
    }
    return k;

}
