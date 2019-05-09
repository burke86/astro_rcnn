///
/// @package phosim
/// @file event.h
/// @brief event logger header
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Glenn Sembroski (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <fitsio.h>
#include <iostream>
#include <string>
#include <vector>

struct PhotonLog {

    double* pX;
    double* pY;
    double* pZ;
    int* pSurface;
    long long counter;
    long long lastRowWritten;
    long long maxNumPhotons;
    long long bufferSize;
    std::string outputDir;
    std::string eventFitsFileName;
    fitsfile* pFitsFile;
 };

class EventFile {

 private:
    PhotonLog* pfPhotonLog;

 public:
    EventFile(int maxNumPhotons, std::string outputDir);
    EventFile(int maxNumPhotons, std::string outputDir, std::string eventFileName);
    ~EventFile();

    void eventFileInit(int maxNumPhotons, std::string outputDir,std::string eventFileName);
    void logPhoton(double a, double b, double c, int d);
    void eventFileCreate();
    void eventFileWrite(long long firstRow,long long numRowsToWrite);
    void eventFileClose();
};
