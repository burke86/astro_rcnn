///
/// @package phosim
/// @file surface.cpp
/// @brief surface class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "surface.h"

void Surface::setup (long surfaceTotal, long points) {

    radiusCurvature = new double[surfaceTotal]();
    conic = new double[surfaceTotal]();
    height = new double[surfaceTotal]();
    outerRadius = new double[surfaceTotal]();
    innerRadius = new double[surfaceTotal]();
    innerRadius0 = new double[surfaceTotal]();
    three = new double[surfaceTotal]();
    four = new double[surfaceTotal]();
    five = new double[surfaceTotal]();
    six = new double[surfaceTotal]();
    seven = new double[surfaceTotal]();
    eight = new double[surfaceTotal]();
    nine = new double[surfaceTotal]();
    ten = new double[surfaceTotal]();
    centerx = new double[surfaceTotal]();
    centery = new double[surfaceTotal]();
    rmax = new double[surfaceTotal]();

    surfacemed = new int[surfaceTotal]();
    surfacecoating = new int[surfaceTotal]();
    surfacetype = new int[surfaceTotal]();

    profile = new double[surfaceTotal*points]();
    radius = new double[surfaceTotal*points]();
    radiusArea = new double[surfaceTotal*points]();
    normal = new double[surfaceTotal*points]();

}

void Surface::asphere (long surfaceIndex, long points) {
    // Create surface profile and normal map

    double R = radiusCurvature[surfaceIndex]; // radius of curvature (mm)
    double rInner = innerRadius[surfaceIndex]; // inner radius (mm)
    double rOutter = outerRadius[surfaceIndex]; // outter radius (mm)
    double k = conic[surfaceIndex]; // conic constant (unitless)
    double radiusFraction, r, asphere, asphereDerivative; // unitless, mm, mm, mm
    double third = -three[surfaceIndex]*1e3; // 3rd order asphere coefficient (mm)
    double fourth = -four[surfaceIndex]*1e3; // 4th order asphere coefficient (mm)
    double fifth = -five[surfaceIndex]*1e3; // 5th order asphere coefficient (mm)
    double sixth = -six[surfaceIndex]*1e3; // 6th order asphere coefficient (mm)
    double seventh = -seven[surfaceIndex]*1e3; // 7th order asphere coefficient (mm)
    double eighth = -eight[surfaceIndex]*1e3; // 8th order asphere coefficient (mm)
    double ninth = -nine[surfaceIndex]*1e3; // 9th order asphere coefficient (mm)
    double tenth = -ten[surfaceIndex]*1e3; // 10th order asphere coefficient (mm)

    for (long i = 0; i < points; i++) {

        radiusFraction = (static_cast<double>(i))/(static_cast<double>(points) - 1);

        radius[points*surfaceIndex + i] = rInner + (rOutter - rInner)*radiusFraction;
        radiusArea[points*surfaceIndex + i] = sqrt((rOutter*rOutter - rInner*rInner)*radiusFraction + rInner*rInner);

        r = radius[points*surfaceIndex + i];

        asphere = third*pow(r, 3.0) + fourth*pow(r, 4.0) + fifth*pow(r, 5.0) + sixth*pow(r, 6.0) + seventh*pow(r, 7.0) + eighth*pow(r, 8.0) + ninth*pow(r, 9.0) + tenth*pow(r, 10.0);
        asphereDerivative = third*pow(r, 2.0)*3.0 + fourth*pow(r, 3.0)*4.0 + fifth*pow(r, 4.0)*5.0 + sixth*pow(r, 5.0)*6.0 + seventh*pow(r, 6.0)*7.0 + eighth*pow(r, 7.0)*8.0 + ninth*pow(r, 8.0)*9.0 + tenth*pow(r, 9.0)*10.0;

        if (R != 0) {
            // sagitta equation
            profile[points*surfaceIndex + i] = height[surfaceIndex] + r*r/(R*(1.0 + sqrt(1.0 - (1.0 + k)*r*r/(R*R)))) + asphere;
            normal[points*surfaceIndex + i] = r/(R*sqrt(1.0 - (1.0 + k)*r*r/(R*R))) + asphereDerivative;
        } else {
            profile[points*surfaceIndex + i] = height[surfaceIndex] + asphere;
            normal[points*surfaceIndex + i] = asphereDerivative;
        }
        
    }

}
