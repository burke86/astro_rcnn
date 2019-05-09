///
/// @package phosim
/// @file instrumentfiles.h
/// @brief Class to make instrument files for program insturment
///
/// @brief Created by:
/// @author Glenn Sembroski
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef INSTRUMENTFILES_H
#define INSTRUMENTFILES_H

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>

#include "ancillary/readtext.h"
#include "ancillary/random.h"
#include "raytrace/parameters.h"

class InstrumentFiles {

public:

    Random random;
    InstrumentFiles();
    ~InstrumentFiles();

    void makeTrackingFile(std::string trackingFileName, double vistime,
                          double jittertime);
    void readoutPars(readText& focalPlaneLayoutPars,readText& segmentationPars,
                     std::string readoutString, int camConfig);
    void focalPlanePars(readText& focalPlaneLayoutPars,
                        std::string outChipString, int camConfig, int perturbationMode, std::vector<int> pertSurf, std::vector<int> pertSecondSurf, int pertFlag);
    void makeSurfaceMap(std::string opticsFile);

private:
    std::map< std::string, int > fSurfaceMap;
    std::map< std::string, int >::iterator fSurfaceMapPos;
    int _lastSurface;
};
#endif
