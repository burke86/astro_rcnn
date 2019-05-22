///
/// @package phosim
/// @file counter.h
/// @brief header for counter
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <time.h>
#include <string>
#include <vector>
#include <cstdatomic>

struct Clog {
    std::atomic<long long> rejected;
    std::atomic<long long> removed;
    std::atomic<long long> accepted;
    std::atomic<long long> removed_dt;
    std::atomic<long long> totalPhoton;
    double previousCPUTime;
    double previousWallTime;
};

struct Tlog {
    double *throughput;
};

void counterClear(Clog *counterLog);
void counterAdd(Clog *localLog, Clog *counterLog);
void counterInit (Clog *counterLog);
void counterCheck (Clog *counterLog, long sourcecounter, char *name);
void countBad (Clog *counterLog, long long photons, long long *ray);
void countBad_dt (Clog *counterLog, long long photons, long long *ray);
void countGood (Clog *counterLog, long long photons, long long *ray);

void writeThroughputFile (const std::string & outputdir, const std::string & outputfilename, Tlog *throughpuTlog, double minwavelength, double maxwavelength, long nsurf);
void addThroughput (Tlog *throughpuTlog, double minwavelength, double maxwavelength, long surf, long waveindex, long long sourceover);
void initThroughput (Tlog *throughpuTlog, double minwavelength, double maxwavelength, long nsurf);
void writeCentroidFile (const std::string & outputdir, const std::string & outputfilename,
                        long long *source_saturation, long long *source_xpos, long long *source_ypos,
                        std::vector<std::string> source_id, long nsource);

// Waits for a new tick. Returns number of ticks since program start.
inline clock_t GetNewTick() {
    clock_t prev = clock();
    clock_t cur = clock();
    while (cur == prev) {
        cur = clock();
    }
    return cur;
}

// Converts ticks to seconds
inline double TicksToSec(clock_t ticks) {
    return static_cast<double>(ticks/CLOCKS_PER_SEC);
}
