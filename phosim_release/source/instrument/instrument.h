///
/// @package phosim
/// @file instrument.cpp
/// @brief instrument  application.
///
/// @brief Created by
/// @author En-Hsing Peng (Purdue)
///
/// @brief Modified by
/// @author Glenn Sembroski (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///
#ifndef INSTRUMENT_H
#define INSTRUMENT_H

#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "instrumentfiles.h"
#include "distortion.h"
#include "ancillary/readtext.h"
#include "ancillary/random.h"
#include "raytrace/constants.h"


InstrumentFiles instrumentFiles;

const int gDevCol = 10;
const double kPi (3.141592653589793238462643);
const int gTotalNumTags = 27;

#endif
