///
/// @package phosim
/// @file constants.h
/// @brief fundamental constants
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

// MATH

#ifdef M_PI
    #undef M_PI
#endif

#define PI        (3.141592653589793238462643)
#define ARCSEC    (PI/180.0/3600.0)
#define DEGREE    (PI/180.0)
#define HALFSQ5M1 (0.5*(sqrt(5)-1.0))
#define HALF3MSQ5 (0.5*(3.0-sqrt(5)))

// BASE PHYSICS

#define H_PLANCK      6.62607015e-34   // J s              (exact)
#define C_LIGHT       2.99792458e8     // m s^-1           (exact)
#define E_CHARGE      1.6022176634e-19 // Coloumbs         (exact)
#define K_BOLTZMANN   1.380649e-23     // Joules K^-1      (exact)
#define G_NEWTON      6.67408e-11      // Newton kg^-2 m^2 (+/- 3.1e-15)

// DERIVED/MEASURED PHYSICS

#define H_CGS (H_PLANCK*1e7)                            // ergs s  (exact)
#define C_CGS (C_LIGHT*100.0)                           // cm s^-1 (exact)
#define EPSILON_0 (1e7/C_LIGHT/C_LIGHT/4.0/PI/100.0)    // Farads cm^-1
#define EPSILON_SI 11.7                                 // Relative permittivity for Si

// ASTRONOMY

#define RADIUS_EARTH 6371.0      // km
#define EARTH_SUN 149597870.700  // km
#define EARTH_MOON 384403.0      // km
