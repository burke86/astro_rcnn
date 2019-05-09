///
/// @package phosim
/// @file obstruction.h
/// @brief header file for obstruction class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

class Obstruction {

public:
    double par[MAX_SURF];
    double center[MAX_SURF];
    double width[MAX_SURF];
    double depth[MAX_SURF];
    double angle[MAX_SURF];
    double reference[MAX_SURF];
    double height[MAX_SURF];
    int type[MAX_SURF];
    int nspid;
    int pupil;

};
