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

#include <assert.h>
#include <math.h>

#include "random.h"

Random::Random() {

    m_z = 1234;
    m_w = 42;
    m_z_reseed = 1234;
    m_w_reseed = 42;
    m_z_correl = 1234;
    m_w_correl = 42;
    intToDouble = 1.0/(uint32_MAX + 1.0);
    plusBiasDouble = intToDouble/2;

}

double Random::normal() {
    /// @fn double Random::normal()
    /// @brief Returns a normally distributed random number
    /// with mean of 0 and variance of 1.   Uses the
    /// Marsaglia & Bray 1964, SIAM Vol 6, 260 algorithm
    /// (Superior to Box-Mueller transforms).

    double u1=0.0;
    double u2=0.0;
    double v1=0.0;
    double v2=0.0;
    double s=2.0;
    while (s >= 1.0 || s == 0.0) {
        u1 = uniform();
        u2 = uniform();
        v1 = 2.0*u1 - 1.0;
        v2 = 2.0*u2 - 1.0;
        s = pow(v1, 2) + pow(v2, 2);
    }
    return(v1*sqrt((-2.0*log(s))/s));
}

double Random::exponential() {
    /// @fn double Random::exponential()
    /// @brief Returns a exponential distributed random number
    /// with rate of 1.  Uses an analytic form.

    return -log(uniform());
}

long long Random::poisson(double lambda) {
    /// @fn long long Random::poisson(long long lambda)
    /// @brief Returns a poisson distributed random number
    /// with rate of lambda.  Uses normal distribution approximation
    /// for lambda greater than 10.  For lambda less than 10 uses
    /// standard binomial distribution.

    long long count;
    double f;
    double g;
    long long d;

    if (lambda < 10) {
        f = uniform();
        d = 0;
        while (f >= exp(-lambda)) {
            g = uniform();
            f = f*g;
            d++;
        }
        count = d;
    } else {
        count = static_cast<long long>(lambda + sqrt(lambda)*normal());
        if (count < 0) count = 0;
    }
    return(count);
}


double Random::uniform() {
    /// @fn double Random::uniform()
    /// @brief Generates uniform random number in the interval
    /// from 0 to 1.  Converts uniformly generated 32 bit
    /// integer to double precision deviate.

    return uniformUint32()*intToDouble + plusBiasDouble;
}


Uint32 Random::uniformUint32() {
    /// @fn Uint32 Random::randomUint32()
    /// @brief Generates uniform random number for a 32-bit
    /// integer.  Based on multiply with carry methodology by George Marsaglia
    /// See for example Marsaglia 2003, JOASM Vol 2.
    /// An extremely fast algorithm which has an astronomically
    /// long repeat time (>2^60).

    m_z = 36969 * (m_z & 65535) + (m_z >> 16);
    m_w = 18000 * (m_w & 65535) + (m_w >> 16);
    return ((m_z << 16) + m_w);
}

void Random::unwind(long count) {
    /// @fn void Random::unwind(long count)
    /// @brief Unwinds random number generator by count values
    /// to remove any initial bias.

    for (long i = 0; i < count; i++) {
        uniformUint32();
    }

}

void Random::setSeed(Uint32 z, Uint32 w) {
    /// @fn void Random::setSeed(Uint32 z, Uint32 w)
    /// @brief Sets the pair of seed values.

    assert(w != 0 && z != 0);
    m_z = z;
    m_w = w;

}

void Random::getSeed(Uint32 *z, Uint32 *w) {
    /// @fn void Random::setSeed(Uint32 z, Uint32 w)
    /// @brief Gets the pair of seed values.

    *z = m_z;
    *w = m_w;

}

// set seed from 64-bit value
void Random::setSeed64(Uint64 seed) {
    /// @fn void Random::setSeed64(Uint64 seed)
    /// @brief Sets seed using 64 bit value.

    Uint32 z = Uint32(seed >> 32);
    Uint32 w = Uint32(seed & 0xFFFFFFFF);
    if (z == 0) z = 4232482;
    if (w == 0) w = 1234628;

    setSeed(z, w);
}

void Random::setSeed32(Uint32 seed) {
    /// @fn void Random::setSeed32(Uint32 seed)
    /// @brief Sets 64 bit seed from 32 mixed up value.

    Uint32 lobits = seed;
    Uint32 hibits = seed*91*53;
    Uint64 seed64 = (Uint64(hibits) << 32) + lobits;

    setSeed64(seed64);
}

time_t Random::timeSec1970() {
    /// @fn time_t Random::timeSec1970()
    /// @brief Returns time.
    return time(NULL);
}

void Random::setSeedFromTime() {
    /// @fn void Random::setSeedFromTime()
    /// @brief Sets seed from time.

    Uint32 lobits = Uint32(clock()) + Uint32(timeSec1970());
    Uint32 hibits = Uint32(clock())*53 + Uint32(timeSec1970())*91;
    Uint64 seed64 = ((Uint64)hibits << 32) + lobits;

    setSeed64(seed64);
}

// SECOND SET OF RANDOM NUMBER FUNCTIONS WHERE WE MESS WITH SEED

double Random::normalFixed() {
    /// @fn double Random::normalFixed()
    /// @brief Returns a normally distributed random number
    /// with mean of 0 and variance of 1.   Uses the
    /// Marsaglia & Bray 1964, SIAM Vol 6, 260 algorithm
    /// (Superior to Box-Mueller transforms).

    double u1 = 0.0;
    double u2 = 0.0;
    double v1 = 0.0;
    double v2 = 0.0;
    double s = 2.0;
    while (s >= 1.0 || s == 0.0) {
        u1 = uniformFixed();
        u2 = uniformFixed();
        v1 = 2.0*u1 - 1.0;
        v2 = 2.0*u2 - 1.0;
        s = pow(v1, 2) + pow(v2, 2);
    }
    return(v1*sqrt((-2.0*log(s))/s));
}

double Random::exponentialFixed() {
    /// @fn double Random::exponentialFixed()
    /// @brief Returns a exponential distributed random number
    /// with rate of 1.  Uses an analytic form.

    return -log(uniformFixed());
}

long long Random::poissonFixed(double lambda) {
    /// @fn long long Random::poissonFixed(long long lambda)
    /// @brief Returns a poisson distributed random number
    /// with rate of lambda.  Uses normal distribution approximation
    /// for lambda greater than 10.  For lambda less than 10 uses
    /// standard binomial distribution.

    long long count;
    double f;
    double g;
    long long d;

    if (lambda < 10) {
        f = uniformFixed();
        d = 0;
        while (f >= exp(-lambda)) {
            g = uniformFixed();
            f = f*g;
            d++;
        }
        count = d;
    } else {
        count = static_cast<long long>(lambda + sqrt(lambda)*normalFixed());
        if (count < 0) count = 0;
    }
    return(count);
}


double Random::uniformFixed() {
    /// @fn double Random::uniformFixed()
    /// @brief Generates uniform random number in the interval
    /// from 0 to 1.  Converts uniformly generated 32 bit
    /// integer to double precision deviate.

    return uniformUint32Fixed()*intToDouble + plusBiasDouble;
}


Uint32 Random::uniformUint32Fixed() {
    /// @fn Uint32 Random::randomUint32Fixed()
    /// @brief Generates uniform random number for a 32-bit
    /// integer.  Based on multiply with carry methodology by George Marsaglia
    /// See for example Marsaglia 2003, JOASM Vol 2.
    /// An extremely fast algorithm which has an astronomically
    /// long repeat time (>2^60).

    m_z_reseed = 36969 * (m_z_reseed & 65535) + (m_z_reseed >> 16);
    m_w_reseed = 18000 * (m_w_reseed & 65535) + (m_w_reseed >> 16);
    return ((m_z_reseed << 16) + m_w_reseed);
}

void Random::unwindFixed(long count) {
    /// @fn void Random::unwindFixed(long count)
    /// @brief Unwinds random number generator by count values
    /// to remove any initial bias.

    for (long i = 0; i < count; i++) {
        uniformUint32Fixed();
    }

}

void Random::setSeedFixed(Uint32 z, Uint32 w) {
    /// @fn void Random::setSeedFixed(Uint32 z, Uint32 w)
    /// @brief Sets the pair of seed values.

    assert(w != 0 && z != 0);
    m_z_reseed = z;
    m_w_reseed = w;

}

// set seed from 64-bit value
void Random::setSeed64Fixed(Uint64 seed) {
    /// @fn void Random::setSeed64Fixed(Uint64 seed)
    /// @brief Sets seed using 64 bit value.

    Uint32 z = Uint32(seed >> 32);
    Uint32 w = Uint32(seed & 0xFFFFFFFF);
    if (z == 0) z = 4232482;
    if (w == 0) w = 1234628;

    setSeedFixed(z, w);
}

void Random::setSeed32Fixed(Uint32 seed) {
    /// @fn void Random::setSeed32Fixed(Uint32 seed)
    /// @brief Sets 64 bit seed from 32 mixed up value.

    Uint32 lobits = seed;
    Uint32 hibits = seed*91*53;
    Uint64 seed64 = (Uint64(hibits) << 32) + lobits;

    setSeed64Fixed(seed64);
}

/// The time correlation random functions

double Random::normalCorrel(double time, double t0) {
    /// @fn double Random::normalFixed()
    /// @brief Returns a normally distributed random number
    /// with mean of 0 and variance of 1.   Uses the
    /// Marsaglia & Bray 1964, SIAM Vol 6, 260 algorithm
    /// (Superior to Box-Mueller transforms).

    double u1 = 0.0;
    double u2 = 0.0;
    double v1 = 0.0;
    double v2 = 0.0;
    double s = 2.0;
    int noseed = 0;
    while (s >= 1.0 || s == 0.0) {
        u1 = uniformCorrel(time, t0, noseed);
        noseed = 1;
        u2 = uniformCorrel(time, t0, noseed);
        v1 = 2.0*u1 - 1.0;
        v2 = 2.0*u2 - 1.0;
        s = pow(v1, 2) + pow(v2, 2);
    }
    return(v1*sqrt((-2.0*log(s))/s));
}

double Random::exponentialCorrel(double time, double t0) {
    /// @fn double Random::exponentialFixed()
    /// @brief Returns a exponential distributed random number
    /// with rate of 1.  Uses an analytic form.

    return -log(uniformCorrel(time, t0, 0));
}

long long Random::poissonCorrel(double lambda, double time, double t0) {
    /// @fn long long Random::poissonFixed(long long lambda)
    /// @brief Returns a poisson distributed random number
    /// with rate of lambda.  Uses normal distribution approximation
    /// for lambda greater than 10.  For lambda less than 10 uses
    /// standard binomial distribution.

    long long count;
    double f;
    double g;
    long long d;

    if (lambda < 10) {
        f = uniformCorrel(time, t0, 0);
        d = 0;
        while (f >= exp(-lambda)) {
            g = uniformCorrel(time, t0, 1);
            f = f*g;
            d++;
        }
        count = d;
    } else {
        count = static_cast<long long>(lambda + sqrt(lambda)*normalCorrel(time, t0));
        if (count < 0) count = 0;
    }
    return(count);
}


double Random::uniformCorrelWrap() {
    /// @fn double Random::uniformFixed()
    /// @brief Generates uniform random number in the interval
    /// from 0 to 1.  Converts uniformly generated 32 bit
    /// integer to double precision deviate.

    return uniformUint32Correl()*intToDouble + plusBiasDouble;
}


Uint32 Random::uniformUint32Correl() {
    /// @fn Uint32 Random::randomUint32Fixed()
    /// @brief Generates uniform random number for a 32-bit
    /// integer.  Based on multiply with carry methodology by George Marsaglia
    /// See for example Marsaglia 2003, JOASM Vol 2.
    /// An extremely fast algorithm which has an astronomically
    /// long repeat time (>2^60).

    m_z_correl = 36969 * (m_z_correl & 65535) + (m_z_correl >> 16);
    m_w_correl = 18000 * (m_w_correl & 65535) + (m_w_correl >> 16);
    return ((m_z_correl << 16) + m_w_correl);
}

void Random::unwindCorrel(long count) {
    /// @fn void Random::unwindFixed(long count)
    /// @brief Unwinds random number generator by count values
    /// to remove any initial bias.

    for (long i = 0; i < count; i++) {
        uniformUint32Correl();
    }

}

void Random::setSeedCorrel(Uint32 z, Uint32 w) {
    /// @fn void Random::setSeedFixed(Uint32 z, Uint32 w)
    /// @brief Sets the pair of seed values.

    assert(w != 0 && z != 0);
    m_z_correl = z;
    m_w_correl = w;

}

// set seed from 64-bit value
void Random::setSeed64Correl(Uint64 seed) {
    /// @fn void Random::setSeed64Fixed(Uint64 seed)
    /// @brief Sets seed using 64 bit value.

    Uint32 z = Uint32(seed >> 32);
    Uint32 w = Uint32(seed & 0xFFFFFFFF);
    if (z == 0) z = 4232482;
    if (w == 0) w = 1234628;

    setSeedCorrel(z, w);
}

void Random::setSeed32Correl(Uint32 seed) {
    /// @fn void Random::setSeed32Fixed(Uint32 seed)
    /// @brief Sets 64 bit seed from 32 mixed up value.

    Uint32 lobits = seed;
    Uint32 hibits = seed*91*53;
    Uint64 seed64 = (Uint64(hibits) << 32) + lobits;

    setSeed64Correl(seed64);
}


double Random::uniformCorrel(double time, double t0, int noseed) {

    Uint32 seed0 = floor(time/t0);
    Uint32 seed1 = seed0 + 1;
    double dt = time/t0 - seed0;
    if (noseed == 0) setSeed32Correl(seed0);
    double a = uniformCorrelWrap();
    if (noseed == 0) setSeed32Correl(seed1);
    double b = uniformCorrelWrap();
    return a*(1-dt)+b*dt;

}
