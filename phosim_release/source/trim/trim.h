///
/// @package phosim
/// @file trim.h
/// @brief trim header file
///
/// @brief Created by:
/// @author Alan Meert (Purdue)
///
/// @brief Modified by:
/// @author Justin Bankert (Purdue)
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "raytrace/parameters.h"
#include <vector>
#include <string>

class Trim {

public:
    std::string instrdir;
    std::vector<std::string> catalog;
    std::vector<std::string> chipid;
    std::string catalogOpd;
    char outputFilename[maxChip][4096];
    char outputOpdFilename[4096];
    int nChip;
    int nCatalog;
    int buffer;
    int filter;
    int strayLight;
    int minSource;
    int opdMode;
    int flipX;
    int flipY;
    int swap;
    std::string obshistid;
    long flatDirectory;
    double pointingRa;
    double pointingDec;
    double rotatez;
    double focalLength;
    double plateScale;
    double scale;
    double extendedBuffer;
    double xPosition[maxChip];
    double yPosition[maxChip];
    double xDimension[maxChip];
    double yDimension[maxChip];
    double pixelSize[maxChip];
    double angle[maxChip];
    double deltaX[maxChip];
    double deltaY[maxChip];

    void xyPosition(double alpha, double delta, double *x, double *y);
    void readCatalog();
    void readOpdCatalog();
    void getDetectorProperties(int d);
    void setup();

};
