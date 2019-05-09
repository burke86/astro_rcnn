///
/// @package phosim
/// @file medium.h
/// @brief header file for medium class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

class Medium {

 public:
    long indexRefractionNumber[MAX_SURF];
    double *indexRefraction[MAX_SURF];
    double *indexRefractionWavelength[MAX_SURF];

    void setup(int surfaceIndex, std::string mediumFileName);

};
