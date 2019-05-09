///
/// @package phosim
/// @file rng_mwc.cpp
/// @brief random number generator functions
///
/// @brief Created by:
/// @author Kreso Cosic (SLIP)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef RANDOM_H
#define RANDOM_H

#ifdef WIN32 //For microsoft's C++ compiler
#include <climits>

typedef __int64             Int64;
typedef unsigned __int64   Uint64;
typedef __int32             Int32;
typedef unsigned __int32   Uint32;
typedef __int16             Int16;
typedef unsigned __int16   Uint16;
typedef __int8              Int8 ;
typedef unsigned __int8    Uint8 ;

const Uint32 uint32_MAX = _UI32_MAX;
#elif __GNUC__  //for g++

#define __STDC_LIMIT_MACROS
#include <stdint.h>

typedef int64_t  Int64;
typedef uint64_t Uint64;
typedef int32_t  Int32;
typedef uint32_t Uint32;
typedef int16_t  Int16;
typedef uint16_t Uint16;
typedef int8_t   Int8;
typedef uint8_t  Uint8;

const Uint32 uint32_MAX = 0xFFFFFFFFU;

#else // for all other compilers, standard-compilant includes
#include <cstdint.h>

typedef std:: int64_t  Int64;
typedef std::uint64_t Uint64;
typedef std:: int32_t  Int32;
typedef std::uint32_t Uint32;
typedef std:: int16_t  Int16;
typedef std::uint16_t Uint16;
typedef std:: int8_t   Int8;
typedef std::uint8_t  Uint8;

const Uint32 uint32_MAX = std::uint32_MAX;
#endif

#include <time.h>

#define PHOSIM_RANDOM_ALIGN 256

class Random {

public:

    Uint32 m_z;
    Uint32 m_w;
    Uint32 m_z_reseed;
    Uint32 m_w_reseed;
    Uint32 m_z_correl;
    Uint32 m_w_correl;
    double intToDouble;
    double plusBiasDouble;

    char filler1[24];
    char filler2[64];
    char filler3[128];

    Random();

    double uniform();
    double normal();
    double exponential();
    long long poisson(double lambda);
    Uint32 uniformUint32();
    void unwind(long count);
    void setSeed(Uint32 z, Uint32 w);
    void setSeed64(Uint64 seed);
    void setSeed32(Uint32 seed);

    void setSeedFromTime();
    time_t timeSec1970();
    void getSeed(Uint32 *z, Uint32 *w);

    double uniformFixed();
    double normalFixed();
    double exponentialFixed();
    long long poissonFixed(double lambda);
    Uint32 uniformUint32Fixed();
    void unwindFixed(long count);
    void setSeedFixed(Uint32 z, Uint32 w);
    void setSeed64Fixed(Uint64 seed);
    void setSeed32Fixed(Uint32 seed);

    double uniformCorrel(double time, double t0, int noseed);
    double uniformCorrelWrap();
    double normalCorrel(double time, double t0);
    double exponentialCorrel(double time, double t0);
    long long poissonCorrel(double lambda, double time, double t0);
    Uint32 uniformUint32Correl();
    void unwindCorrel(long count);
    void setSeedCorrel(Uint32 z, Uint32 w);
    void setSeed64Correl(Uint64 seed);
    void setSeed32Correl(Uint32 seed);


} __attribute__ ((aligned (PHOSIM_RANDOM_ALIGN)));

#endif
