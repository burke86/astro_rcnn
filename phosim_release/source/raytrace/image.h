///
/// @package phosim
/// @file image.h
/// @brief header for image class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <fitsio.h>
#include <fitsio2.h>
#include <fftw3.h>
#include <vector>

#include "ancillary/random.h"
#include "galaxy.h"
#include "dust.h"
#include "observation.h"
#include "raytrace.h"
#include "surface.h"
#include "coating.h"
#include "silicon.h"
#include "perturbation.h"
#include "screen.h"
#include "air.h"
#include "medium.h"
#include "obstruction.h"
#include "chip.h"
#include "contamination.h"
#include "grating.h"
#include "photon.h"
#include "state.h"
#include "lock.h"

class Image : public Observation {

    public:

    // objects and structures
    Galaxy galaxy;
    Dust dust;
    Surface surface;
    Coating coating;
    Silicon silicon;
    Air air;
    Perturbation perturbation;
    Screen screen;
    Medium medium;
    Obstruction obstruction;
    Contamination contamination;
    Chip chip;
    Grating* pGrating;
    State state;
    Random* random;

    // setup and loop methods
    int atmSetup();
    int telSetup();
    int sourceLoop();
    void photonLoop(long ssource, long thread, int finish);

    // thread and optimization
    int dynamicTransmissionOptimization(long k, long *lastSurface, long *preGhost, long waveSurfaceIndex, long straylightcurrent, Photon *aph);
    static void* threadFunction(void *voidArgs);
    Lock lock;
    long remain;
    long openthreads;
    int *openthread;

    // physics methods
    int getWavelengthTime(Photon *aph, long source);
    int domeSeeing(Vector *angle, Photon *aph);
    int tracking(Vector *angle, double time);
    int atmosphericDispersion(Vector *angle, Photon *aph, long layer);
    int largeAngleScattering(Vector *largeAngle, Photon *aph);
    int secondKick(Vector *largeAngle, Photon *aph);
    int diffraction(Vector *position, Vector angle, Vector *largeAngle, Photon *aph);
    int samplePupil(Vector *position, long long ray, Photon *aph);
    int transmissionCheck(double transmission, long surfaceIndex, long waveSurfaceIndex, Photon *aph);
    int transmissionPreCheck(long surfaceIndex, long waveSurfaceIndex, Photon *aph);
    int chooseSurface(long *newSurface, long *oldSurface, Photon *aph);
    int findSurface(Vector angle, Vector position, double *distance, long surfaceIndex, Photon *aph);
    int goldenBisectSurface(double a, double b, double c, double *z, Vector angle,
                            Vector position, double *distance, long surfaceIndex, Photon *aph);
    int getIntercept(double x, double y, double *z, long surfaceIndex, Photon *aph);
    int getDeltaIntercept(double x, double y, double *zv, long surfaceIndex, Photon *aph);
    int bloom(int saturatedFlag, Photon *aph);
    void saturate(long source, Vector *largeAngle, Photon *aph, double shiftX, double shiftY);
    int photonSiliconPropagate(Vector *angle, Vector *position, double lambda, Vector normal,
                               double dh, long waveSurfaceIndex, Photon *aph);
    int electronSiliconPropagate(Vector *angle, Vector *position, Photon *aph);
    int contaminationSurfaceCheck(Vector position, Vector *angle, long surfaceIndex, Photon *aph);
    double airIndexRefraction(Photon *aph, long layer);
    double surfaceCoating(double wavelength, Vector angle,
                          Vector normal, long surfaceIndex, double *reflection, Photon *aph);
    double cloudOpacity(long layer, Photon *aph);
    double cloudOpacityMoon(long layer, Photon *aph);
    double atmosphereOpacity(Vector angle, long layer, Photon *aph);
    double atmosphereOpacityMoon(Vector angle, long layer, Photon *aph);
    double fringing (Vector angle, Vector normal, double wavelength, double nSi, double thickness, double meanFreePath, int polarization);
    void getAngle(Vector *angle, double time, long source);
    void getDeltaAngle(Vector *angle, Vector *position, long source, double *shiftX, double *shiftY, int thread, int *initGal, Photon *aph);
    void newRefractionIndex(long surfaceIndex, Photon *aph);
    void atmospherePropagate(Vector *position, Vector angle, long layer, int mode, Photon *aph);
    void atmosphereIntercept(Vector *position, long layer, Photon *aph);
    void atmosphereRefraction(Vector *angle, long layer, int mode, Photon *aph);
    void atmosphereDiffraction(Vector *angle, Photon *aph);
    void transform(Vector *angle, Vector *position, long surfaceIndex, int focusFlag, Photon *aph);
    void transformInverse(Vector *angle, Vector *position, long surfaceIndex);
    void interceptDerivatives(Vector *normal, Vector position, long surfaceIndex);
    void cosmicRays(long long *raynumber);

    // output
    void writeImageFile();
    void writeOPD();
    void readCheckpoint(int checkpointcount);
    void writeCheckpoint(int checkpointcount);
    void cleanup();

};

struct thread_args {
    Image *instance;
    long runthread;
    long ssource[4096];
    long thread;
};
