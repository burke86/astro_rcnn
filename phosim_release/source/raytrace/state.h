///
/// @package phosim
/// @file state.h
/// @brief header for state structure
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "event.h"
#include "counter.h"
#include <atomic>

struct State {

    std::atomic<int> *satupmap;
    std::atomic<int> *satdownmap;
    double *opd;
    double *opdcount;
    std::atomic<double> *dynamicTransmission;
    std::atomic<double> *dynamicTransmissionLow;
    std::atomic<unsigned long> *focal_plane;
    float *focal_plane_fl;
    double *cx;
    double *cy;
    double *cz;
    double *r0;
    double *epR;

    EventFile* pEventLogging;
    Clog counterLog;
    Tlog throughputLog;
    Clog globalLog;

};
