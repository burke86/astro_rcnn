///
/// @package phosim
/// @file photon.h
/// @brief header for photon class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

struct Photon {
    int direction;
    int polarization;
    double shiftedAngle;
    double airRefraction;
    double airRefractionPrevious;
    double airRefractionADC;
    double airRefractionPreviousADC;
    double adcx, adcy, adcz;
    double nprev, ncurr;
    double wavelength;
    double wavelengthFactor;
    long indexlx0, indexly0, indexlx1, indexly1;
    long indexcx0, indexcy0, indexcx1, indexcy1;
    long indexmx0, indexmy0, indexmx1, indexmy1;
    long indexfx0, indexfy0, indexfx1, indexfy1;
    double dlx, dly, dcx, dcy, dmx, dmy, dfx, dfy;
    long uuint, vvint, wwint;
    double windx, windy;
    long lindex;
    double xp, yp;
    double xporig, yporig, zporig;
    double opdx, opdy;
    long xPos, yPos;
    double xPosR, yPosR;
    double xpos, ypos;
    double time, prtime;
    double absoluteTime;
    double dvr;
    double op;
    long oindex;
    long counter;
    long maxcounter;
    int ghostFlag;
    long sourceOver_m;
    long sourceSaturationRadius;
    int saturationFlag;
    double z0;
    double collect_z;
    long xindex;
    double rxindex;
    long location;
    double saveRand[MAX_BOUNCE];
    double refractionInvariant;
    double refractionInvariantCenter;
    long thread;
};
