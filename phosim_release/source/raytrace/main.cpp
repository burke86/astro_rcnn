///
/// @package phosim
/// @file main.cpp
/// @brief main for raytrace code
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "raytrace/image.h"
int main() {
    Image im;

    im.parser();
    im.background();
    im.atmSetup();
    im.telSetup();
    im.sourceLoop();
    im.cleanup();

    return(0);
}
