///
/// @package phosim
/// @file main.cpp
/// @brief main for trim program
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
#include "trim/trim.h"

int main() {

    Trim trim;

    trim.setup();
    if (trim.opdMode == 0) {
        for (int d = 0; d < trim.nChip; d++) {
            trim.getDetectorProperties(d);
        }
        trim.readCatalog();
    } else {
        trim.readOpdCatalog();
    }
    return(0);

}
