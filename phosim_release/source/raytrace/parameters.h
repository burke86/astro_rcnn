///
/// @package phosim
/// @file parameters.h
/// @brief parameters
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef Par_H
#define Par_H

// general numerical accuracy
int const SURFACE_POINTS = 102400;
int const PERTURBATION_POINTS = 1024;
int const SCREEN_SIZE = 1024;
int const OPD_SCREEN_SIZE = 255;
int const OPD_SAMPLING = 4;
int const bin_number = 8000;
double const scale = 10.0;
int const MAX_SOURCE = 1000000;
int const MAX_SURF = 30;
int const MAX_LAYER = 100;
int const MAX_IMAGE = 100;
int const MAX_BOUNCE = 100;
int const NZERN = 28;
int const NCHEB = 21;

int const RAYTRACE_MAX_ITER = 10;
double const RAYTRACE_TOLERANCE = 1e-7;      // mm
double const RAYTRACE_MAX_LENGTH = 40000.0;  // mm
double const RAYTRACE_ERROR = 10.0;          // mm

int const maxChip = 400;
int const maxCatalog = 20;
int const SILICON_STEPS = 1024;
int const SILICON_SUB_STEPS = 256;

#endif
