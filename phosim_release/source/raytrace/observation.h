///
/// @package phosim
/// @file observation.h
/// @brief observation header file
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <string>
#include <vector>
#include "parameters.h"
#include "vector_type.h"
#include "source.h"

class Observation {

public:
    double rotationjitter;
    double zenith;
    double azimuth;
    double rotationrate;
    double elevationjitter;
    double azimuthjitter;
    double windjitter;
    double windshake;
    double groundlevel;
    double xtelloc;
    double ytelloc;
    double latitude;
    double longitude;
    std::vector<double> seefactor;
    std::vector<double> wind;
    std::vector<double> winddir;
    std::vector<double> outerscale;
    std::vector<double> height;
    std::vector<double> cloudmean;
    std::vector<double> cloudvary;
    std::vector<std::string> extraCommandString;
    double *dtau;
    double pressure;
    double waterPressure;
    double temperature;
    double tempvar;
    double raynorm;
    double o2norm;
    double o3norm;
    double h2onorm;
    double raygradient;
    double o2gradient;
    double o3gradient;
    double h2ogradient;
    double aerosolgradient;
    double rayAngle;
    double o2Angle;
    double o3Angle;
    double h2oAngle;
    double aerosolAngle;
    double aerosoltau;
    double aerosolindex;
    std::vector<std::string> pertType;
    int NTERM;
    std::vector<std::vector<double> > izernike;
    std::vector<int> surfaceLink;
    std::vector<std::vector<double> > body;
    double minr;
    double shuttererror;
    double zodiacalLightMagnitude;
    double overrideZodiacalLightMagnitude;
    double transtol, backAlpha, backBeta, backRadius, backGamma, backDelta,
        backBuffer, backEpsilon, np, finiteDistance;
    long satbuffer;
    long activeBuffer;
    double screentol;
    double maxr;
    double domeseeing;
    double toypsf;
    double pixsize;
    double pra;
    double pdec;
    double spiderangle;
    double platescale;
    double centerx;
    double centery;
    double tai;
    double day;
    double laScatterProb;
    double exptime;
    double vistime;
    double timeoffset;
    double rotatex;
    double rotatey;
    double rotatez;
    double *sedCorr;
    double *sedDwdp;
    double *sedW;
    double *sedC;
    double largeScale;
    double coarseScale;
    double mediumScale;
    double fineScale;
    double totalnorm;
    double totalseeing;
    double moonalt;
    double moondist;
    double moonra;
    double moondec;
    double phaseang;
    double solarzen;
    double solaralt;
    double lst;
    double hourangle;
    double airglowcintensity;
    double airglowpintensity;
    double watervar;
    double minwavelength;
    double maxwavelength;
    double centralwavelength;
    double domelight;
    double domewave;
    double chipangle;
    double decenterx;
    double decentery;
    int telconfig;
    int checkpointcount;
    int checkpointtotal;
    int impurityvariation;
    int fieldanisotropy;
    int fringeflag;
    int detectorcollimate;
    int deadlayer;
    int chargesharing;
    int pixelerror;
    int chargediffusion;
    int photoelectric;
    int airrefraction;
    int aberration;
    int nutation;
    int precession;
    int onlyadjacent;
    int splitsources;
    int focusSteps;
    int focusFlag;
    double raydensity;
    double scalenumber;
    double devvalue;
    double airglowvariation;
    double focusStepZ;
    double focusStepX;
    float nbulk;
    float nf;
    float nb;
    float sf;
    float sb;
    float sensorthickness;
    float impurityX;
    float impurityY;
    float overdepbias;
    float siliconreflectivity;
    float sensorTempNominal;
    float sensorTempDelta;
    float qevariation;
    float *airglow;
    long long nphot;
    long long *sourceXpos;
    long long *sourceYpos;
    long long *sourcePhoton;
    long airglowScreenSize;
    long telescopeMode;
    long backgroundMode;
    long coatingmode;
    long poissonMode;
    long contaminationmode;
    long trackingMode;
    long windShakeMode;
    long ranseed;
    long obsseed;
    long zernikemode;
    long atmospheric_dispersion;
    long atmosphericdispcenter;
    long natmospherefile;
    long straylight;
    double straylightcut;
    long detectorMode;
    long opticsonlymode;
    long additivemode;
    long diffractionMode;
    long spiderMode;
    long pupilscreenMode;
    long spaceMode;
    long sequenceMode;
    long aperturemode;
    long areaExposureOverride;
    long filter;
    std::string filterName;
    long saturation;
    long eventfile;
    long opdfile;
    long opdsize;
    long opdsampling;
    long centroidfile;
    long throughputfile;
    long pixelsx, pixelsy, minx, maxx, miny, maxy;
    std::string obshistid;
    long pairid;
    long blooming;
    long well_depth;
    long *sedN;
    long *sedPtr;
    long nsedptr;
    long sedptr;
    long nsource;
    long nimage;
    long nreallocs;
    long sedMax;
    long ghostonly;
    long atmdebug;
    long largeGrid;
    long coarseGrid;
    long mediumGrid;
    long fineGrid;
    long nsnap;
    int nframes;
    int nskip;
    int ngroups; // not exposed
    std::vector<int> ghost;
    int flatdir;
    int tarfile;
    int date;
    int numthread;
    std::string devmaterial;
    std::string devtype;
    std::string devmode;
    std::vector<int> feaflag;
    std::vector<std::string> feafile;
    std::vector<double> feascaling;
    double normwave;
    int atmospheremode;
    int opacitymode;
    int sourceperthread;

    // should be part of image but have conflict with settings.c
    long nsurf;
    long nmirror;
    double airmass;

    // remainder should be ok
    std::vector<std::string> atmospherefile;
    std::vector<std::string> cloudfile;
    std::string trackingfile;
    long trackinglines;
    std::string outputfilename;
    std::string chipid;
    std::string focalplanefile;
    std::string eventFitsFileName;
    std::string outputdir, workdir, seddir, imagedir, datadir, instrdir, bindir;
    Vector tpx, tpy, tpz;

    long naxesb[MAX_IMAGE][2];
    float *tempptr[MAX_IMAGE];
    float *cumulativex[MAX_IMAGE];
    float cumulative[MAX_IMAGE];
    float cumulativef[MAX_IMAGE];

    Source sources;

    // functions
    int parser();
    int background();
    int addSource(const std::string & object, int sourcetype);
    int addOpd(const std::string & opd);
    int header(fitsfile *faptr);
    int settings();
    int filterTruncateSources();
    int splitSources();
    void readSed(const std::string & filename, int mode);

};
