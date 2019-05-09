///
/// @package phosim
/// @file test_trim.cpp
/// @brief trim unit test
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

#include <math.h>

#include "trim/trim.h"
#include "validation/unittest.h"
#include "ancillary/readtext.h"

int main() {

    Trim trim;

    trim.setup();
    unitTestOutput(1, trim.nCatalog, "trim::setup", "Number of Catalogs", "none");
    unitTestOutput(1, trim.nChip, "trim::setup", "Number of Chips", "none");
    for (int d = 0; d < trim.nChip; d++) {
        trim.getDetectorProperties(d);
    }
    unitTestOutput(20480.0, trim.xDimension[0], "trim::getDetectorProperties", "x dimension", "microns");
    unitTestOutput(20480.0, trim.yDimension[0], "trim::getDetectorProperties", "y dimension", "microns");
    unitTestOutput(10.0, trim.pixelSize[0], "trim::getDetectorProperties", "pixel size", "microns");
    trim.readCatalog();
    readText catalog("trimcatalog_0_testchip.pars");
    for (size_t t(0); t < catalog.getSize(); t++) {
        std::string line(catalog[t]);
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;
        if (keyName == "object") {
            int i;
            iss >> i;
            unitTestOutput(1.0, i, "trim::readCatalog", "object ID", "none");
            continue;
        }
    }
    return(0);

}
