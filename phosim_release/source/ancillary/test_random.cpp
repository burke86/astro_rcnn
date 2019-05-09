///
/// @package phosim
/// @file test_random.cpp
/// @brief unit tests for random class.
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <math.h>
#include <stdio.h>
#include "ancillary/random.h"
#include "validation/unittest.h"

int main() {

    Random random;
    int flag;
    int flag2;
    double f;
    double g = 0.0;
    double h = 0.0;
    long n = 10000000;
    double m = static_cast<double>(n);


    random.setSeedFromTime();
    random.setSeed32Fixed(0);

    // Uniform unit test
    flag = 1;
    flag2 = 1;
    g = 0.0;
    h = 0.0;
    for (long i = 0; i < n; i++) {
        f = random.uniform();
        if (f < 0.0) flag = 0;
        if (f >= 1.0) flag2 = 0;
        g += f;
        h += f*f;
    }
    unitTestOutput(static_cast<double>(flag2), 1.0, "random::uniform", "less than 1", "none");
    unitTestOutput(static_cast<double>(flag), 1.0, "random::uniform", "greater than or equal to 0", "none");
    unitTestOutput(g/m, 0.5, 1e-2, "random::uniform", "mean", "none");
    unitTestOutput(sqrt(h/m - g*g/m/m), 1.0/sqrt(12.0), 1e-2, "random::uniform", "standard deviation", "none");

    // Normal unit test
    flag = 1;
    flag2 = 1;
    g = 0.0;
    h = 0.0;
    for (long i = 0; i < n; i++) {
        f = random.normal();
        g += f;
        h += f*f;
    }
    unitTestOutput(g/m, 0.0, 1e-2, "random::normal", "mean", "none");
    unitTestOutput(sqrt(h/m - g*g/m/m), 1.0, 1e-2, "random::normal", "standard deviation", "none");


    // Exponential unit test
    flag = 1;
    flag2 = 1;
    g = 0.0;
    h = 0.0;
    for (long i = 0; i < n; i++) {
        f = random.exponential();
        if (f < 0.0) flag = 0;
        g += f;
        h += f*f;
    }
    unitTestOutput(static_cast<double>(flag), 1.0, "random::exponential", "greater than or equal to 0", "none");
    unitTestOutput(g/m, 1.0, 1e-2, "random::exponential", "mean", "none");
    unitTestOutput(sqrt(h/m - g*g/m/m), 1.0, 1e-2, "random::exponential", "standard deviation", "none");

    // Poisson unit test
    flag = 1;
    flag2 = 1;
    g = 0.0;
    h = 0.0;
    for (long i = 0; i < n; i++) {
        f = random.poisson(5.0);
        if (f < 0.0) flag = 0;
        g += f;
        h += f*f;
    }
    unitTestOutput(static_cast<double>(flag), 1.0, "random::poisson", "greater than or equal to 0", "none");
    unitTestOutput(g/m, 5.0, 1e-2, "random::poisson", "mean", "none");
    unitTestOutput(sqrt(h/m - g*g/m/m), sqrt(5.0), 1e-2, "random::poisson", "standard deviation", "none");


    // Uniform unit test
    flag = 1;
    flag2 = 1;
    g = 0.0;
    h = 0.0;
    for (long i = 0; i < n; i++) {
        f = random.uniformFixed();
        if (f < 0.0) flag = 0;
        if (f >= 1.0) flag2 = 0;
        g += f;
        h += f*f;
    }
    unitTestOutput(static_cast<double>(flag2), 1.0, "random::uniformFixed", "less than 1", "none");
    unitTestOutput(static_cast<double>(flag), 1.0, "random::uniformFixed", "greater than or equal to 0", "none");
    unitTestOutput(g/m, 0.5, 1e-2, "random::uniformFixed", "mean", "none");
    unitTestOutput(sqrt(h/m - g*g/m/m), 1.0/sqrt(12.0), 1e-2, "random::uniformFixed", "standard deviation", "none");

    // Normal unit test
    flag = 1;
    flag2 = 1;
    g = 0.0;
    h = 0.0;
    for (long i = 0; i < n; i++) {
        f = random.normalFixed();
        g += f;
        h += f*f;
    }
    unitTestOutput(g/m, 0.0, 1e-2, "random::normalFixed", "mean", "none");
    unitTestOutput(sqrt(h/m - g*g/m/m), 1.0, 1e-2, "random::normalFixed", "standard deviation", "none");


    // Exponential unit test
    flag = 1;
    flag2 = 1;
    g = 0.0;
    h = 0.0;
    for (long i = 0; i < n; i++) {
        f = random.exponentialFixed();
        if (f < 0.0) flag = 0;
        g += f;
        h += f*f;
    }
    unitTestOutput(static_cast<double>(flag), 1.0, "random::exponentialFixed", "greater than or equal to 0", "none");
    unitTestOutput(g/m, 1.0, 1e-2, "random::exponentialFixed", "mean", "none");
    unitTestOutput(sqrt(h/m - g*g/m/m), 1.0, 1e-2, "random::exponentialFixed", "standard deviation", "none");

    // Poisson unit test
    flag = 1;
    flag2 = 1;
    g = 0.0;
    h = 0.0;
    for (long i = 0; i < n; i++) {
        f = random.poissonFixed(5.0);
        if (f < 0.0) flag = 0;
        g += f;
        h += f*f;
    }
    unitTestOutput(static_cast<double>(flag), 1.0, "random::poissonFixed", "greater than or equal to 0", "none");
    unitTestOutput(g/m, 5.0, 1e-2, "random::poissonFixed", "mean", "none");
    unitTestOutput(sqrt(h/m - g*g/m/m), sqrt(5.0), 1e-2, "random::poissonFixed", "standard deviation", "none");

}
