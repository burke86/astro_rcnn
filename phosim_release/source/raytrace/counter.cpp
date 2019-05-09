///
/// @package phosim
/// @file counter.cpp
/// @brief various logging functions
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <stdio.h>
#include <stdlib.h>

#include "sys/time.h"
#include "constants.h"
#include "counter.h"

void counterInit(Clog *counterLog) {


    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Type                Sources         Photons  (Sat,Rem,Rej,Acc)%%  Time (s)       Photons/s\n");
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");

    counterLog->rejected = 0;
    counterLog->accepted = 0;
    counterLog->removed = 0;
    counterLog->removed_dt = 0;
    counterLog->totalPhoton = 0;
    struct timeval tim;
    gettimeofday(&tim, NULL);
    counterLog->previousWallTime = tim.tv_sec + (tim.tv_usec/1000000.0);
    counterLog->previousCPUTime = TicksToSec(GetNewTick());

}

void counterClear(Clog *counterLog) {

    counterLog->rejected = 0;
    counterLog->accepted = 0;
    counterLog->removed = 0;
    counterLog->removed_dt = 0;
    counterLog->totalPhoton = 0;

}

void counterAdd(Clog *localLog, Clog *counterLog) {

    counterLog->rejected += localLog->rejected;
    counterLog->accepted += localLog->accepted;
    counterLog->removed += localLog->removed;
    counterLog->removed_dt += localLog->removed_dt;
    counterLog->totalPhoton += localLog->totalPhoton;

}

void counterCheck(Clog *counterLog, long sourcecounter, char *name) {

    double newCpuTime, newWallTime;
    long rate;
    char sourceString[4096];
    char photonString[4096];
    char rateString[4096];
    long long tp;

    newCpuTime = TicksToSec(clock());
    struct timeval tim;
    gettimeofday(&tim, NULL);
    newWallTime = tim.tv_sec + (tim.tv_usec/1000000.0);
    rate = static_cast<long>(counterLog->totalPhoton/(newWallTime-counterLog->previousWallTime + 1e-2));

    if (sourcecounter < 1000) sprintf(sourceString, "%7ld", sourcecounter);
    if (sourcecounter >= 1000) sprintf(sourceString, "%3ld,%03ld", sourcecounter/1000, sourcecounter%1000);

    tp = counterLog->totalPhoton;
    if (tp < 1000) sprintf(photonString, "%15lld", tp);
    if (tp >= 1000 && tp < 1000000) {
        sprintf(photonString, "%11lld,%03lld", tp/1000, tp%1000);
    }
    if (tp >= 1000000 && tp < 1000000000) {
        sprintf(photonString, "%7lld,%03lld,%03lld", tp/1000000,
                (tp/1000)%1000, tp%1000);
    }
    if (tp >= 1000000000) {
        sprintf(photonString,"%3lld,%03lld,%03lld,%03lld", (tp/1000000000),
                (tp/1000000)%1000, (tp/1000)%1000,
                tp%1000);
    }

    if (rate < 1000) sprintf(rateString, "%15ld", rate);
    if (rate >= 1000 && rate < 1000000) sprintf(rateString, "%11ld,%03ld", rate/1000, rate%1000);
    if (rate >= 1000000 && rate < 1000000000) sprintf(rateString, "%7ld,%03ld,%03ld", rate/1000000,
                                                      (rate/1000)%1000, rate%1000);
    if (rate >= 1000000000) sprintf(rateString, "%3ld,%03ld,%03ld,%03ld", (rate/1000000000),
                                    (rate/1000000)%1000, (rate/1000)%1000, rate%1000);

    if (counterLog->totalPhoton > 1) {
        fprintf(stdout,"%s %s %s  (%3.0f,%3.0f,%3.0f,%3.0f)  %9.2f %s\n",name, sourceString, photonString,
                static_cast<double>(counterLog->rejected)/static_cast<double>(counterLog->totalPhoton)*100,
                static_cast<double>(counterLog->removed_dt)/static_cast<double>(counterLog->totalPhoton)*100,
                static_cast<double>(counterLog->removed)/static_cast<double>(counterLog->totalPhoton)*100,
                static_cast<double>(counterLog->accepted)/static_cast<double>(counterLog->totalPhoton)*100,
                newWallTime - counterLog->previousWallTime, rateString);
    }

    counterLog->previousWallTime = newWallTime;
    counterLog->previousCPUTime = newCpuTime;
    counterLog->rejected = 0;
    counterLog->accepted = 0;
    counterLog->removed = 0;
    counterLog->removed_dt = 0;
    counterLog->totalPhoton = 0;

}


void countGood(Clog *counterLog, long long photons, long long *ray) {

    counterLog->accepted += 1;
    counterLog->rejected += photons - 1;
    counterLog->totalPhoton += photons;
    *ray += photons;

}

void countBad(Clog *counterLog, long long photons, long long *ray) {

    counterLog->removed += photons;
    counterLog->totalPhoton += photons;
    *ray += photons;

}

void countBad_dt(Clog *counterLog, long long photons, long long *ray) {

    counterLog->removed_dt += photons;
    counterLog->totalPhoton += photons;
    *ray += photons;

}

void addThroughput (Tlog *throughputlog, double minwavelength, double maxwavelength, long surf, long waveindex, long long sourceover) {

    throughputlog->throughput[(surf + 1)*(static_cast<int>(maxwavelength - minwavelength + 1)) + waveindex] += static_cast<double>(sourceover);
}

void initThroughput (Tlog *throughputlog, double minwavelength, double maxwavelength, long nsurf) {

    throughputlog->throughput = static_cast<double*>(calloc((nsurf + 2)*(maxwavelength - minwavelength + 1), sizeof(double)));

}

void writeThroughputFile (const std::string & outputdir, const std::string & outputfilename, Tlog *throughputlog, double minwavelength, double maxwavelength, long nsurf) {

    FILE *outdafile;
    long i, k;
    char tempstring[4096];

        sprintf(tempstring, "%s/throughput_%s.txt", outputdir.c_str(), outputfilename.c_str());
        outdafile = fopen(tempstring, "w");
        for (k = 0; k < maxwavelength - minwavelength + 1; k++) {
            fprintf(outdafile, "%ld ", k);
            for (i = 0; i < nsurf + 2; i++) {
                fprintf(outdafile, "%lf ", throughputlog->throughput[i*(static_cast<int>(maxwavelength - minwavelength + 1)) + k]);
            }
            fprintf(outdafile, "\n");
        }
        fclose(outdafile);


}

void writeCentroidFile (const std::string & outputdir, const std::string & outputfilename,
                        long long *sourceSaturation, long long *sourceXpos, long long *sourceYpos,
                        std::vector<std::string> source_id, long nsource) {

    FILE *outdafile;
    long k;
    char tempstring[4096];


        sprintf(tempstring, "%s/centroid_%s.txt", outputdir.c_str(), outputfilename.c_str());
        outdafile = fopen(tempstring, "w");
        fprintf(outdafile, "SourceID Photons AvgX AvgY\n");
        for (k = 0; k < nsource; k++) {
            fprintf(outdafile, "%s %lld %lf %lf\n", source_id[k].c_str(), sourceSaturation[k],
                    (static_cast<double>(sourceXpos[k]))/(static_cast<double>(sourceSaturation[k])),
                    (static_cast<double>(sourceYpos[k]))/(static_cast<double>(sourceSaturation[k])));
        }
        fclose(outdafile);

}
