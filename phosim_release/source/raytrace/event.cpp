///
/// @package phosim
/// @file event.cpp
/// @brief event logger
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

#include "event.h"

EventFile::EventFile(int maxNumPhotons, std::string outputDir) {

    std::string eventFitsFileName = "output.fits";
    eventFileInit(maxNumPhotons, outputDir, eventFitsFileName );

}

EventFile::EventFile(int maxNumPhotons, std::string outputDir,
                     std::string eventFitsFileName) {
    eventFileInit(maxNumPhotons, outputDir, eventFitsFileName);
}
EventFile::~EventFile() {

    delete pfPhotonLog->pX;
    delete pfPhotonLog->pY;
    delete pfPhotonLog->pZ;
    delete pfPhotonLog->pSurface;
}
void EventFile::eventFileInit(int maxNumPhotons, std::string outputDir,
                              std::string eventFitsFileName) {

// Initalize the photonlog

    pfPhotonLog = new PhotonLog();

    pfPhotonLog->maxNumPhotons = maxNumPhotons;
    // 200000 is about as big as we can make these arrays and not seg fault

    if (pfPhotonLog->maxNumPhotons > 200000){
        pfPhotonLog->bufferSize = 200000;
    } else {
        pfPhotonLog->bufferSize = pfPhotonLog->maxNumPhotons;
    }
    //We use arrays here since thats what the c fits library call want.
    // We use the "p" pre-lable to indicate these are pointers to arrays.
    pfPhotonLog->pX = new double[pfPhotonLog->bufferSize];
    pfPhotonLog->pY = new double[pfPhotonLog->bufferSize];
    pfPhotonLog->pZ = new double[pfPhotonLog->bufferSize];
    pfPhotonLog->pSurface = new int[pfPhotonLog->bufferSize];
    pfPhotonLog->counter = 0;
    pfPhotonLog->lastRowWritten = 0;
    pfPhotonLog->pFitsFile = NULL;

    pfPhotonLog->outputDir = outputDir;
    pfPhotonLog->eventFitsFileName = eventFitsFileName;
    return;
}

void EventFile::logPhoton (double a, double b, double c, int d) {
// Save photons into the columun arrays. If arrays are full write them to
// the eventsfile. Create file when/(if) its  needed.
// This extention to > 200000 events is needed for calibration modeling
    //  See if we have already closed the file
    if (pfPhotonLog->counter < 0){
        return;
    }

    // Save this photon

    int logIndex = pfPhotonLog->counter-pfPhotonLog->lastRowWritten;

    pfPhotonLog->pX[logIndex]       = a;
    pfPhotonLog->pY[logIndex]       = b;
    pfPhotonLog->pZ[logIndex]       = c;
    pfPhotonLog->pSurface[logIndex] = d;
    pfPhotonLog->counter++;


    // See if this event will be the last we want.
    // Or is the array full?
    // setting photons.maxNumPhotonsons to -1 will also cause file to be closed*
    // Note pfPhotonLog->counter has been bumped by 1. Counting rows now

    long long numRowsToWrite = 0;
    numRowsToWrite = pfPhotonLog->counter - pfPhotonLog->lastRowWritten;

    long long fStartRow = 0;

    if (pfPhotonLog->counter >= pfPhotonLog->maxNumPhotons ||
        numRowsToWrite%pfPhotonLog->bufferSize == 0){
        // Write this group out.
        if (pfPhotonLog->pFitsFile == NULL){
            // file not created yet. do that
            eventFileCreate();
            fStartRow = 1;
        } else {
            // File already exists. Pick up where we left off
            fStartRow = pfPhotonLog->lastRowWritten + 1;
        }

        // Now write out the photons
        if (numRowsToWrite>0 && pfPhotonLog->maxNumPhotons>0){
            eventFileWrite(fStartRow, numRowsToWrite);
            pfPhotonLog->lastRowWritten = pfPhotonLog->counter;
        }
        // close file ?
        if (pfPhotonLog->counter >= pfPhotonLog->maxNumPhotons){
            eventFileClose();
        }
        return;
    }
}


void EventFile::eventFileCreate() {
// Create the output fits events file and the events table using cfitsio*
// library calls

    int status;
    std::string tempstring;

    tempstring = "!" + pfPhotonLog->outputDir + "/" + pfPhotonLog->eventFitsFileName;

    status = 0;
    fits_create_file(&pfPhotonLog->pFitsFile, tempstring.c_str(), &status);

    if (status != 0) {
        std::cout << "eventFileCreate failed to create event file: " << tempstring << " Status: " << status << std::endl;
        pfPhotonLog->pFitsFile = NULL;
        return;
    }

    // Now that the file is made, define an extention table for the events
    char* ttype[4];
    char* tform[4];
    ttype[0] = (char*)"x";
    ttype[1] = (char*)"y";
    ttype[2] = (char*)"z";
    ttype[3] = (char*)"surface";
    tform[0] = (char*)"1D";
    tform[1] = (char*)"1D";
    tform[2] = (char*)"1D";
    tform[3] = (char*)"1J";

    status = 0; //Must always init this to 0=ok before calling cfitsio routine*
    fits_create_tbl(pfPhotonLog->pFitsFile, BINARY_TBL, 0, 4, ttype, tform, NULL, NULL, &status);
    if (status != 0){
        std::cout << "eventFileCreate failed to create event file table. Status:" << status << std::endl;
        pfPhotonLog->pFitsFile = NULL;
        return;
    }
    return;
}

void EventFile::eventFileWrite(long long firstRow, long long numRowsToWrite) {
// Write whatever phots need writing to the rows specified from the
// buffers.

    int status;
    status = 0; //Must always init this to 0=ok before calling cfitsio routine
    fits_write_col(pfPhotonLog->pFitsFile,TDOUBLE, 1, firstRow, 1, numRowsToWrite,
                   pfPhotonLog->pX, &status);
    fits_write_col(pfPhotonLog->pFitsFile,TDOUBLE, 2, firstRow, 1, numRowsToWrite,
                   pfPhotonLog->pY, &status);
    fits_write_col(pfPhotonLog->pFitsFile, TDOUBLE, 3, firstRow, 1, numRowsToWrite,
                   pfPhotonLog->pZ, &status);
    fits_write_col(pfPhotonLog->pFitsFile, TINT, 4, firstRow, 1, numRowsToWrite,
                   pfPhotonLog->pSurface, &status);
    if (status != 0){
        std::cout << "eventFileWrite failed to write data to event file. Status: " << std::endl;
        pfPhotonLog->pFitsFile = NULL;
        return;
    }
    return;
}
void EventFile::eventFileClose() {
// Clean up (write partial buffers if we need to) and close up the events
// fits file

    int status;
    long long fStartRow;
    long long numRowsToWrite = 0;

    // File already closed?
    if (pfPhotonLog->counter == -1){
        return;
    }

    // Check to see if we need to write out some last events
    if (pfPhotonLog->counter>pfPhotonLog->lastRowWritten) {
        if (pfPhotonLog->pFitsFile == NULL) {
            // file not created yet. do that
            eventFileCreate();
            fStartRow = 1;
            numRowsToWrite = pfPhotonLog->counter;
        } else {
            // File already exists. Pick up where we left off
            fStartRow = pfPhotonLog->lastRowWritten + 1;
            numRowsToWrite = pfPhotonLog->counter-pfPhotonLog->lastRowWritten;
        }
        // Now write out the photons
        if (numRowsToWrite > 0){
            eventFileWrite(fStartRow, numRowsToWrite);
        }
    }
    // Check that the file was ever opened
    if (pfPhotonLog->pFitsFile != NULL){
        status = 0;
        fits_close_file(pfPhotonLog->pFitsFile, &status);
        if (status != 0) {
            std::cout << "eventFileClose failed to close event file. Status: " << std::endl;
        }
    }
    pfPhotonLog->counter = -1;  // Flag that file is closed
    return;
}
