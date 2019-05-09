///
/// @package phosim
/// @file surface.h
/// @brief header file for surface class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///


class Coating {

public:
    double** reflection;
    double** transmission;
    double** wavelength;
    double** angle;

    long* wavelengthNumber;
    long* angleNumber;

    void setup(long totalSurface);
    void allocate(long surfaceIndex, long lines);

};
