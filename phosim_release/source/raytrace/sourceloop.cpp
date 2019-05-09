///
/// @package phosim
/// @file sourceloop.cpp
/// @brief source loop
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Colin Burke (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///
#include <pthread.h>
#include <unistd.h>

int Image::sourceLoop() {

    char tempstring[4096];
    long sourceCounter;
    long surfaceLimit;
    pthread_t *thread;
    thread_args *args;

    //    OPTIMIZATION INITIALIZATION
    thread = (pthread_t*)malloc(numthread*sizeof(pthread_t));
    args = (thread_args*)malloc(numthread*sizeof(thread_args));
    openthread = static_cast<int*>(malloc(numthread*sizeof(int)));
    surfaceLimit = natmospherefile*2 + nsurf*2 + 2;
    state.dynamicTransmission = static_cast<std::atomic<double>*>
        (malloc(surfaceLimit*(maxwavelength - minwavelength + 1)*sizeof(std::atomic<double>)));
    state.dynamicTransmissionLow = static_cast<std::atomic<double>*>
        (malloc(surfaceLimit*(maxwavelength - minwavelength + 1)*sizeof(std::atomic<double>)));
    for (int i = 0; i < maxwavelength - minwavelength + 1; i++) {
        for (int j = 0; j < surfaceLimit; j++) {
            state.dynamicTransmission[i*surfaceLimit + j] = -1.0;
            state.dynamicTransmissionLow[i*surfaceLimit + j] = 1.0;
        }
    }
    if (checkpointcount != 0) readCheckpoint(checkpointcount);
    counterInit(&state.counterLog);
    counterClear(&state.globalLog);
    if (opdfile == 1) pthread_mutex_init(&lock.lock4, NULL);
    if (opdfile == 1) pthread_mutex_init(&lock.lock5, NULL);
    if (opdfile == 1) pthread_mutex_init(&lock.lock6, NULL);
    if (opdfile == 1) pthread_mutex_init(&lock.lock7, NULL);
    pthread_mutex_init(&lock.lock8, NULL);
    if (opdfile == 1) pthread_cond_init(&lock.cond, NULL);
    if (opdfile == 1) {
        remain = 0;
        sourceperthread = 1;
    }
    // integration sequence bookkeeping
    int nsamples = nframes + nskip; // total number of samples in sequence
    int nframei = 0; // current frame number within group
    int ngroupi = 0; // current group number
    std::string outputfilenameold = outputfilename;

    //    INTEGRATION SEQUENCE LOOP
    for (int nframe = 0; nframe < ngroups*nsamples - nskip; nframe++) {

        for (int i = 0; i < numthread; i++) {
            openthread[i] = 0;
        }
        openthreads = 0;
        long subsource = 0;

        std::ostringstream outfile;
        if (ngroups > 1 || nsamples > 1) {
            std::cout << "------------------------------------------------------------------------------------------" << std::endl;
            std::cout << "Integration Sequence: Group " << (ngroupi + 1) << " Frame " << (nframei + 1) << std::endl;
            std::cout << "------------------------------------------------------------------------------------------" << std::endl;
            outfile << outputfilenameold << "_R" << std::setfill('0') << std::setw(3) << nframe;
        } else {
            outfile << outputfilenameold;
        }
        outputfilename = outfile.str();

        nframei++;
        // do not simulate frame skips
        if (nframei % nframes == 0) {
            nframe += nskip;
            nframei = 0;
            ngroupi++;
            tai += nskip*devvalue/86400.0;
        }
        if (nframe > 0) tai += devvalue/86400.0;

        //    LOGGING INITIALIZATION
        if (eventfile) {
            state.pEventLogging = new EventFile((int)(nphot*20), workdir, eventFitsFileName+".fits");
        } else {
            state.pEventLogging = NULL;
        }
        pGrating = new Grating();
        if (throughputfile) initThroughput(&state.throughputLog, minwavelength, maxwavelength, nsurf);
        
        //    SOURCE TYPE LOOP
        for (int sourceType = 646; sourceType >= 0; sourceType--) {

            sourceCounter = 0;

            if (static_cast<int>(round(sourceType*checkpointtotal/546)) == checkpointcount) {

                //    SOURCE LOOP
                // long subsource = 0;

                for (long source = 0; source < nsource; source++) {

                    if ((sourceType >= 0   && sourceType <=  99 && sources.type[source] == 0 && (source%100) == sourceType) ||
                        (sourceType >= 100 && sourceType <= 199 && sources.type[source] == 1 && (source%100) == (sourceType - 100)) ||
                        (sourceType >= 200 && sourceType <= 299 && sources.type[source] == 2 && (source%100) == (sourceType - 200)) ||
                        (sourceType >= 300 && sourceType <= 399 && sources.type[source] == 3 && (source%100) == (sourceType - 300)) ||
                        (sourceType >= 400 && sourceType <= 499 && sources.type[source] == 4 && (source%100) == (sourceType - 400)) ||
                        (sourceType >= 500 && sourceType <= 599 && sources.type[source] == 5 && (source%100) == (sourceType - 500)) ||
                        (sourceType == 600 && sources.type[source] == 6 && sources.mag[source] > 32.25) ||
                        (sourceType >= 601 && sourceType <= 645 && sources.type[source] == 6 &&
                         sources.mag[source] >  (332.25 - static_cast<double>(sourceType)/2.0) &&
                         sources.mag[source] <= (332.75 - static_cast<double>(sourceType)/2.0))
                        || (sourceType == 646 && sources.type[source] == 6 && sources.mag[source] <= 9.75)) {

                    while (openthreads >= numthread*sourceperthread) {
                        usleep(100*sourceperthread);
                    }
                        pthread_mutex_lock(&lock.lock8);
                        long bestthread = sourceperthread + 1;
                        for (int i = 0; i < numthread; i++) {
                            if ((openthread[i] < sourceperthread) && (openthread[i] < bestthread)) {
                                subsource = i;
                                bestthread = openthread[i];
                            }
                        }
                        openthread[subsource] += 1;
                        openthreads++;
                        sourceCounter++;
                        args[subsource].ssource[openthread[subsource] - 1] = source;
                        if (openthread[subsource] == sourceperthread) {
                            args[subsource].instance = this;
                            args[subsource].thread = subsource;
                            args[subsource].runthread = sourceperthread;
                            if (opdfile == 1) {
                                if (((sourceCounter-1) % numthread) == 0) {
                                    if (source < static_cast<long>(floor(nsource/numthread))*numthread) {
                                        remain = numthread;
                                    } else {
                                        remain = nsource - static_cast<long>(floor(nsource/numthread))*numthread;
                                    }
                                }
                            }
                            pthread_create(&thread[subsource], NULL, &Image::threadFunction, &args[subsource]);
                            if (opdfile == 0) pthread_detach(thread[subsource]);
                        }
                        pthread_mutex_unlock(&lock.lock8);
                        if (opdfile == 1) {
                            if ((((sourceCounter-1) % numthread) == (numthread - 1)) ||
                                (source == nsource - 1)) {
                                for (long ss = 0; ss < numthread; ss++) {
                                    if (openthread[ss] > 0) {
                                        pthread_join(thread[ss], NULL);
                                    }
                                }
                            }
                        }

                    }

                }
            }
            if (sourceType >= 0  && sourceType < 100)  sprintf(tempstring, "Dome Light     %3d%% ", (100 - sourceType));
            if (sourceType >= 100 && sourceType < 200) sprintf(tempstring, "Airglow-Coll   %3d%% ", (200 - sourceType));
            if (sourceType >= 200 && sourceType < 300) sprintf(tempstring, "Airglow-Phot   %3d%% ", (300 - sourceType));
            if (sourceType >= 300 && sourceType < 400) sprintf(tempstring, "Moon-Rayleigh  %3d%% ", (400 - sourceType));
            if (sourceType >= 400 && sourceType < 500) sprintf(tempstring, "Moon-Mie       %3d%% ", (500 - sourceType));
            if (sourceType >= 500 && sourceType < 600) sprintf(tempstring, "Zodiacal       %3d%% ", (600 - sourceType));
            if (sourceType == 600)                     sprintf(tempstring, "Astro Object m>32.0 ");
            if (sourceType >= 601 && sourceType < 646) sprintf(tempstring, "Astro Object m=%4.1f ", 332.5 - static_cast<double>(sourceType)/2.0);
            if (sourceType == 646)                     sprintf(tempstring, "Astro Object m<10.0 ");
            if (sourceCounter > 0) counterCheck(&state.counterLog, sourceCounter, tempstring);
        }
        for (int i=0; i<numthread; i++) {
            if (openthread[i] != sourceperthread && openthread[i] != 0) {
                args[i].instance = this;
                args[i].thread = i;
                args[i].runthread = openthread[i];
                pthread_create(&thread[i], NULL, &Image::threadFunction, &args[i]);
                pthread_detach(thread[i]);
            }
        }
        while (openthreads != 0) {
            usleep(100*sourceperthread);
        }

        // COSMIC RAYS
        long long ray;
        long long detRay;
        if (checkpointcount == checkpointtotal) {
            if (backgroundMode > 0) {
                detRay = 0;
                ray = 0;
                cosmicRays(&detRay);
                if (detRay > 0) {
                    for (long i = 0; i < detRay; i++) {
                        countGood(&state.counterLog, 1, &ray);
                    }
                    sprintf(tempstring, "Cosmic Rays         ");
                    counterCheck(&state.counterLog, detRay, tempstring);
                }
            }
        }

        // OUTPUT DATA
        if (checkpointcount == checkpointtotal) writeImageFile();
        if (opdfile) writeOPD();
        if (checkpointcount !=  checkpointtotal) writeCheckpoint(checkpointcount);
        if (centroidfile) writeCentroidFile(workdir, outputfilename, sourcePhoton, sourceXpos, sourceYpos, sources.id, nsource);
        if (throughputfile) writeThroughputFile(workdir, outputfilename, &state.throughputLog, minwavelength, maxwavelength, nsurf);
        if (eventfile) state.pEventLogging->eventFileClose();
    }

    return(0);

}
