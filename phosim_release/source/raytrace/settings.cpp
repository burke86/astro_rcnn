///
/// @package phosim
/// @file settings.cpp
/// @brief settings (part of observation class)
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Alan Meert (Purdue)
/// @author Colin Burke (Purdue)
/// @author Caleb Remocaldo (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "ancillary/fits.h"

int Observation::parser () {

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Photon Raytrace" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Installing Universe." << std::endl;

    double altitude;
    int extraCommandFlag = 0;
    sourceperthread = 1;
    nphot = 0;
    totalnorm = 0.0;
    nsedptr = 0;
    sedptr = 0;
    nsource = 0;
    nimage = 0;
    nsurf = 0;
    maxr = 4180; //mm
    minr = 2558; //mm
    exptime = 15.0;
    vistime = 33.0;
    nsnap = 2;
    nframes = 1;
    nskip = 0;
    shuttererror = 0.0;
    timeoffset = 0.0;
    pra = 0.0;
    pdec = 0.0;
    rotatex = -1;
    rotatey = -1;
    rotatez = 0.0;
    spiderangle = 0.0;
    zenith = 0.0;
    airmass = 1.0;
    azimuth = 0.0;
    windjitter = 2.5;
    rotationjitter = 1;
    elevationjitter = 0.02;
    azimuthjitter = 0.02;
    windshake = 0.05;
    impurityvariation = 1;
    fieldanisotropy = 1;
    fringeflag = 1;
    deadlayer = 1;
    chargesharing = 1;
    pixelerror = 1;
    poissonMode = 1;
    telescopeMode = 1;
    backgroundMode = 1;
    coatingmode = 1;
    chargediffusion = 1;
    photoelectric = 1;
    detectorcollimate = 0;
    contaminationmode = 1;
    trackingMode = 1;
    windShakeMode = 1;
    detectorMode = 1;
    diffractionMode = 1;
    spiderMode = 1;
    pupilscreenMode = 0;
    spaceMode = -1;
    pressure = 520;
    waterPressure = 8;
    temperature = 5;
    tempvar = 0.0;
    airrefraction = 1;
    raynorm = 1;
    o2norm = 1;
    o3norm = 1;
    h2onorm = 1;
    aerosoltau = 0.02;
    aerosolindex = -1.28;
    o2gradient = 0.0;
    o3gradient = 0.0;
    h2ogradient = 0.0;
    aerosolgradient = 0.0;
    raygradient = 0.0;
    rayAngle = 0.0;
    o2Angle = 0.0;
    o3Angle = 0.0;
    h2oAngle = 0.0;
    aerosolAngle = 0.0;
    ranseed = -1;
    obsseed = 0;
    zernikemode = 1;
    onlyadjacent = 1;
    splitsources = 1;
    atmospheric_dispersion = 1;
    atmosphericdispcenter = 1;
    aberration = 0;
    precession = 0;
    nutation = 0;
    outputdir = "../output";
    workdir = ".";
    seddir = "../data/SEDs";
    imagedir = "../data/images";
    datadir = "../data";
    instrdir = "../data/lsst";
    bindir = "../data/lsst";
    outputfilename = "focalplane";
    chipid = "R22_S11";
    trackingfile = ".";
    natmospherefile = 0;
    straylight = 1;
    straylightcut = 10.0;
    aperturemode = 0;
    ghostonly = 0;
    saturation = 1;
    eventfile = 0;
    opdfile = 0;
    opdsize = OPD_SCREEN_SIZE;
    opdsampling = OPD_SAMPLING;
    eventFitsFileName = "output.fits.gz";
    centroidfile = 0;
    throughputfile = 0;
    filter = 0;
    blooming = 1;
    obshistid = "0";
    pairid = 0;
    tai = 0.0;
    toypsf = 0.0;
    finiteDistance = 0.0;
    transtol = 0.0;
    backAlpha = 0.1;
    backGamma = 1.0;
    backDelta = 80.0;
    screentol = 0.01;
    overrideZodiacalLightMagnitude = 0;
    backBeta = 4.0;
    backRadius = 10.0; // arcseconds
    backBuffer = 100.0;
    activeBuffer = 10;
    backEpsilon = 1.0;
    np = 0.90;
    satbuffer = 5;
    date = 1;
    overdepbias = -45.0;
    sensorTempNominal = -1; // K
    sensorTempDelta = 0; // K
    qevariation = 0.0;
    airglowvariation = 1.0;
    airglowScreenSize = 1024;
    laScatterProb = 0.135;
    totalseeing = 0.67;
    flatdir = 0;
    tarfile = 0;
    atmdebug = 0;
    largeScale = 1.0;
    coarseScale = 1.0;
    mediumScale = 1.0;
    fineScale = 1.0;
    largeGrid = 1;
    coarseGrid = 1;
    mediumGrid = 1;
    fineGrid = 1;
    moonalt = -1.0*PI/2.0 + 0.00001;
    moondist = PI - 0.00001;
    phaseang = PI - 0.00001;
    solarzen = PI - 0.00001;
    solaralt = 0;
    airglowcintensity = 22.08;
    airglowpintensity = 22.08;
    moonra = 0; // degrees
    moondec = 0; // degrees
    domelight = 1000.0;
    domewave = 0.0;
    raydensity = 0.6;
    scalenumber = 8.0;
    checkpointtotal = 0;
    checkpointcount = 0;
    areaExposureOverride = 0;
    opticsonlymode = 0;
    additivemode = 1;
    focusFlag = 0;

    //update these values when instrdir is set
    //variables in focalplanelayout.txt
    centerx = -1; // microns
    centery = -1; // microns
    pixsize = -1; // microns
    pixelsx = -1;
    pixelsy = -1;
    minx = -1;
    miny = -1;
    maxx = -1;
    maxy = -1;
    sensorthickness = -1; // microns
    // variables in location.txt
    groundlevel = -1;
    xtelloc = -1; // m
    ytelloc = -1; // m
    latitude = -1; // degrees
    longitude = -1; // degrees
    // variables in central_wavelengths.txt
    filterName = "";
    minwavelength = -1;
    maxwavelength = -1;
    centralwavelength = -1;
    platescale = -1;

    // variables in sensor.txt
    well_depth = -1; // electrons
    nbulk = -1;
    nf = -1;
    nb = -1;
    sf = -1;
    sb = -1;
    // variables in tracking.txt
    windjitter = -1;
    windshake = -1;
    rotationjitter = -1;
    elevationjitter = -1;
    azimuthjitter = -1;
    // bookkeeping for atmosphere-related control structure
    opacitymode = 0;
    atmospheremode = 2; // 2 = everything on for the purposes of the code (use Kolmogorov diffraction)
                        // 1 = no turbulence (use pupil diffraction)
                        // 0 = atmosphere all off (skip some setup)

    /* atmosphere parameter arrays */
    atmospherefile.resize(MAX_LAYER);
    cloudfile.resize(MAX_LAYER);
    seefactor.resize(MAX_LAYER, 0);
    wind.resize(MAX_LAYER, 0);
    winddir.resize(MAX_LAYER, 0);
    outerscale.resize(MAX_LAYER, 0);
    height.resize(MAX_LAYER, 0);
    cloudmean.resize(MAX_LAYER, 0);
    cloudvary.resize(MAX_LAYER, 0);
    dtau = static_cast<double*>(calloc(MAX_LAYER, sizeof(double)));

    /* telescope parameter arrays */
    std::vector<std::vector<double> > tbody;
    std::vector<std::vector<double> > tizernike;
    izernike.resize(MAX_SURF);
    tizernike.resize(MAX_SURF);
    pertType.resize(MAX_SURF);
    NTERM = NZERN;
    if ( NZERN < NCHEB ) NTERM = NCHEB;
    for (int i = 0; i < MAX_SURF; i++) {
        izernike[i].resize(NTERM, 0);
        tizernike[i].resize(NTERM, 0);
    }
    for (int i = 0; i < MAX_SURF; i++) {
        surfaceLink.push_back(i);
    }
    body.resize(MAX_SURF);
    tbody.resize(MAX_SURF);
    for (int i = 0; i < MAX_SURF; i++) {
        body[i].resize(6, 0);
        tbody[i].resize(6, 0);
    }
    ghost.resize(MAX_SURF, 0);
    feaflag.resize(MAX_SURF, 0);
    feafile.resize(MAX_SURF*MAX_SURF*2);
    feascaling.resize(MAX_SURF*MAX_SURF*2);

    /* sky parameter arrays */
    sources.vx = static_cast<double*>(calloc(MAX_SOURCE, sizeof(double)));
    sources.vy = static_cast<double*>(calloc(MAX_SOURCE, sizeof(double)));
    sources.vz = static_cast<double*>(calloc(MAX_SOURCE, sizeof(double)));
    sources.norm = static_cast<double*>(calloc(MAX_SOURCE, sizeof(double)));
    sources.mag = static_cast<double*>(calloc(MAX_SOURCE, sizeof(double)));
    sources.type = static_cast<int*>(calloc(MAX_SOURCE, sizeof(int)));
    sedCorr = static_cast<double*>(calloc(MAX_SOURCE, sizeof(double)));
    sedDwdp = static_cast<double*>(calloc(MAX_SOURCE, sizeof(double)));
    sources.spatialtype = static_cast<int*>(calloc(MAX_SOURCE, sizeof(int)));
    sources.dusttype = static_cast<int*>(calloc(MAX_SOURCE, sizeof(int)));
    sources.dusttypez = static_cast<int*>(calloc(MAX_SOURCE, sizeof(int)));
    sources.sedptr = static_cast<long*>(calloc(MAX_SOURCE, sizeof(long)));
    // sources.id = (double*)calloc(MAX_SOURCE, sizeof(double));
    sources.skysameas = static_cast<long*>(calloc(MAX_SOURCE, sizeof(long)));
    sedN = static_cast<long*>(calloc(MAX_SOURCE, sizeof(long)));
    sedPtr = static_cast<long*>(calloc(MAX_SOURCE, sizeof(long)));
    sourceXpos = static_cast<long long*>(calloc(MAX_SOURCE, sizeof(long long)));
    sourceYpos = static_cast<long long*>(calloc(MAX_SOURCE, sizeof(long long)));
    sourcePhoton = static_cast<long long*>(calloc(MAX_SOURCE, sizeof(long long)));
    sources.spatialpar = (double**)calloc(MAX_SOURCE, sizeof(double*));
    for (int i = 0; i < MAX_SOURCE; i++) {
        sources.spatialpar[i] = static_cast<double*>(calloc(14, sizeof(double)));
    }
    sources.dustpar = static_cast<double**>(calloc(MAX_SOURCE, sizeof(double*)));
    for (int i = 0; i < MAX_SOURCE; i++) {
        sources.dustpar[i] = static_cast<double*>(calloc(2, sizeof(double)));
    }
    sources.dustparz = static_cast<double**>(calloc(MAX_SOURCE, sizeof(double*)));
    for (int i = 0; i < MAX_SOURCE; i++) {
        sources.dustparz[i] = static_cast<double*>(calloc(2, sizeof(double)));
    }
    sources.sedfilename.resize(MAX_SOURCE);
    sources.spatialname.resize(MAX_SOURCE);
    sources.dustname.resize(MAX_SOURCE);
    sources.dustnamez.resize(MAX_SOURCE);

    // get instrdir first
    readText pars(std::cin);
    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        readText::get(line, "instrdir", instrdir);
    }

    // read central wavelengths
    std::istringstream wavelengthPars(readText::get(instrdir + "/central_wavelengths.txt", filter));
    std::string filterNameT;
    double minwavelengthT, maxwavelengthT, centralwavelengthT, platescaleT;
    wavelengthPars >> filterNameT >> minwavelengthT >> maxwavelengthT >> centralwavelengthT >> platescaleT;
    if (filterName == "") filterName = filterNameT;
    if (minwavelength == -1) minwavelength = minwavelengthT*1000; // nm
    if (maxwavelength == -1) maxwavelength = maxwavelengthT*1000; // nm
    if (centralwavelength == -1) centralwavelength = centralwavelengthT; // um
    if (platescale == -1) platescale = platescaleT; // um deg-1

    // read tracking file
    // need to know optional telescope rotation before we look for sources
    std::string sss;
    sss = instrdir + "/tracking.txt";
    std::ifstream inStream3(sss.c_str());
    if (inStream3) {
        readText trackingPars(instrdir + "/tracking.txt");
        for (size_t t(0); t < trackingPars.getSize(); t++) {
            std::string line(trackingPars[t]);
            if (windjitter == -1) readText::get(line, "windjitter", windjitter);
            if (windshake == -1) readText::get(line, "windshake", windshake);
            if (rotationjitter == -1) readText::get(line, "rotationjitter", rotationjitter);
            if (elevationjitter == -1) readText::get(line, "elevationjitter", elevationjitter);
            if (azimuthjitter == -1) readText::get(line, "azimuthjitter", azimuthjitter);
            // fiducial field angles (equivalent to rotating entire telescope)
            if (rotatex == -1) {
                if (readText::getKey(line, "rotatex", rotatex)) rotatex *= (-DEGREE);
            }
            if (rotatey == -1) {
                if (readText::getKey(line, "rotatey", rotatey)) rotatey *= (-DEGREE);
            }
        }
        if (windshake == -1) windshake=0.05;
    } else {
        windjitter = 2.5;
        rotationjitter = 0.0;
        elevationjitter = 0.0;
        azimuthjitter = 0.0;
        windshake = 0.0;
    }
    if (rotatex == -1) rotatex = 0.0;
    if (rotatey == -1) rotatey = 0.0;

    // read parameters file
    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;

        if (keyName == "object") {
            std::string object;
            std::getline(iss, object);
            addSource(object,  6);
            continue;
        }
        if (keyName == "opd") {
            std::string opd;
            std::getline(iss, opd);
            addOpd(opd);
            continue;
        }

        readText::get(line, "outputdir", outputdir);
        readText::get(line, "workdir", workdir);
        readText::get(line, "seddir", seddir);
        readText::get(line, "imagedir", imagedir);
        readText::get(line, "datadir", datadir);
        readText::get(line, "bindir", bindir);
        readText::get(line, "sourceperthread", sourceperthread);
        readText::get(line, "thread", numthread);
        readText::get(line, "atmospheremode", atmospheremode);
        readText::get(line, "telescopemode", telescopeMode);
        readText::get(line, "poissonmode", poissonMode);
        readText::get(line, "backgroundmode", backgroundMode);
        readText::get(line, "impurityvariation", impurityvariation);
        readText::get(line, "fieldanisotropy", fieldanisotropy);
        readText::get(line, "fringing", fringeflag);
        readText::get(line, "deadlayer", deadlayer);
        readText::get(line, "chargesharing", chargesharing);
        readText::get(line, "photoelectric", photoelectric);
        readText::get(line, "detectorcollimate", detectorcollimate);
        readText::get(line, "pixelerror", pixelerror);
        readText::get(line, "chargediffusion", chargediffusion);
        readText::get(line, "coatingmode", coatingmode);
        readText::get(line, "contaminationmode", contaminationmode);
        readText::get(line, "trackingmode", trackingMode);
        readText::get(line, "windshakemode", windShakeMode);
        readText::get(line, "detectormode", detectorMode);
        readText::get(line, "diffractionmode", diffractionMode);
        readText::get(line, "spidermode", spiderMode);
        readText::get(line, "pupilscreenmode", pupilscreenMode);
        readText::get(line, "spacemode", spaceMode);
        readText::get(line, "zernikemode", zernikemode);
        readText::get(line, "straylight", straylight);
        readText::get(line, "straylightcut", straylightcut);// not exposed
        readText::get(line, "ghost", ghost);// no exposed
        readText::get(line, "ghostonly", ghostonly);
        readText::get(line, "aperturemode", aperturemode);
        readText::get(line, "onlyadjacent", onlyadjacent);
        readText::get(line, "splitsources", splitsources);
        readText::get(line, "areaexposureoverride", areaExposureOverride);
        readText::get(line, "opticsonlymode", opticsonlymode);
        readText::get(line, "additivemode", additivemode);
        readText::get(line, "minr", minr);
        readText::get(line, "maxr", maxr);
        readText::get(line, "exptime", exptime);
        if (readText::getKey(line, "nsnap", nsnap)) extraCommandFlag = 1;
        readText::get(line, "nframes", nframes);
        readText::get(line, "nskip", nskip);
        readText::get(line, "shuttererror", shuttererror);
        readText::get(line, "timeoffset", timeoffset);
        readText::get(line, "finitedistance", finiteDistance);
        readText::get(line, "transtol", transtol); // not exposed
        readText::get(line, "np", np); // not exposed
        readText::get(line, "satbuffer", satbuffer); // not exposed
        readText::get(line, "activebuffer", activeBuffer);// not exposed
        readText::get(line, "backalpha", backAlpha);// not exposed
        readText::get(line, "backgamma", backGamma);// not exposed
        readText::get(line, "backdelta", backDelta);// not exposed
        readText::get(line, "backepsilon", backEpsilon);// not exposed
        readText::get(line, "screentol", screentol);// not exposed
        if (readText::getKey(line, "zodiacallightmag", zodiacalLightMagnitude)) overrideZodiacalLightMagnitude = 1;
        readText::get(line, "backbeta", backBeta);//  not exposed
        readText::get(line, "backradius", backRadius);// not exposed
        readText::get(line, "backbuffer", backBuffer);// not exposed
        readText::get(line, "date", date);// not exposed
        readText::get(line, "dayofyear", day);// not exposed
        readText::get(line, "flatdir", flatdir);// not exposed
        readText::get(line, "tarfile", tarfile);// not exposed
        readText::get(line, "atmdebug", atmdebug);// not exposed
        readText::get(line, "large_grid", largeGrid);// not exposed
        readText::get(line, "coarse_grid", coarseGrid);// not exposed
        readText::get(line, "medium_grid", mediumGrid);// not exposed
        readText::get(line, "fine_grid", fineGrid);// not exposed
        readText::get(line, "large_scale", largeScale);// not exposed
        readText::get(line, "coarse_scale", coarseScale);// not exposed
        readText::get(line, "medium_scale", mediumScale);// not exposed
        readText::get(line, "fine_scale", fineScale);// not exposed
        readText::get(line, "opdfile", opdfile);// not exposed
        readText::get(line, "opdsampling", opdsampling);
        readText::get(line, "opdsize", opdsize);
        if (readText::getKey(line, "filter", filter)) {
            // read central wavelengths again
            std::istringstream wavelengthParsQ(readText::get(instrdir + "/central_wavelengths.txt", filter));
            std::string filterNameQ;
            double minwavelengthQ, maxwavelengthQ, centralwavelengthQ, platescaleQ;
            wavelengthParsQ >> filterNameQ >> minwavelengthQ >> maxwavelengthQ >> centralwavelengthQ >> platescaleQ;
            filterName = filterNameQ;
            minwavelength = minwavelengthQ*1000; // nm
            maxwavelength = maxwavelengthQ*1000; // nm
            centralwavelength = centralwavelengthQ; // um
            platescale = platescaleQ; // um deg-1
        }
        readText::get(line, "saturation", saturation);
        readText::get(line, "aberration", aberration);
        readText::get(line, "nutation", nutation);
        readText::get(line, "precession", precession);
        readText::get(line, "blooming", blooming);
        readText::get(line, "eventfile", eventfile);
        readText::get(line, "eventFitsFileName", eventFitsFileName);
        readText::get(line, "centroidfile", centroidfile);
        readText::get(line, "throughputfile", throughputfile);
        readText::get(line, "well_depth", well_depth);
        readText::get(line, "nbulk", nbulk);
        readText::get(line, "nf", nf);
        readText::get(line, "nb", nb);
        readText::get(line, "sf", sf);
        readText::get(line, "sb", sb);
        readText::get(line, "sensortempnominal", sensorTempNominal);
        readText::get(line, "sensorthickness", sensorthickness);
        readText::get(line, "overdepbias", overdepbias);
        readText::get(line, "sensortempdelta", sensorTempDelta);
        readText::get(line, "qevariation", qevariation);
        readText::get(line, "obshistid", obshistid); //
        readText::get(line, "exposureid", pairid);//
        readText::get(line, "tai", tai);//
        readText::get(line, "windjitter", windjitter);
        readText::get(line, "windshake", windshake);
        readText::get(line, "rotationjitter", rotationjitter);
        readText::get(line, "elevationjitter", elevationjitter);
        readText::get(line, "azimuthjitter", azimuthjitter);

        if (keyName == "izernike" ) {
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < NTERM; j++) {
                    tizernike[i][j] = izernike[i][j];
                    if (additivemode == 1) izernike[i][j] = 0.0;
                }
            }
            long surfaceIndex;
            readText::get(line, "izernike", izernike);
            iss >> surfaceIndex;
            pertType[surfaceIndex].assign("zern");
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < NTERM; j++) {
                    if (additivemode == 1) izernike[i][j] += tizernike[i][j];
                }
            }
        }
        if (keyName == "ichebyshev" ) {
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < NTERM; j++) {
                    tizernike[i][j] = izernike[i][j];
                    if (additivemode == 1) izernike[i][j] = 0.0;
                }
            }
            long surfaceIndex;
            readText::get(line, "ichebyshev", izernike);
            iss >> surfaceIndex;
            pertType[surfaceIndex].assign("chebyshev");
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < NTERM; j++) {
                    if (additivemode == 1) izernike[i][j] += tizernike[i][j];
                }
            }
        }
        readText::get(line, "surfacelink", surfaceLink);
        if (keyName == "body") {
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < 6; j++) {
                    tbody[i][j] = body[i][j];
                    if (additivemode == 1) body[i][j] = 0.0;
                }
            }
            readText::get(line, "body", body);
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < 6; j++) {
                    if (additivemode == 1) body[i][j] += tbody[i][j];
                }
            }

        }
        readText::get(line, "natmospherefile", natmospherefile);
        readText::get(line, "atmospherefile", atmospherefile);
        readText::get(line, "cloudfile", cloudfile);
        readText::get(line, "trackingfile", trackingfile);
        readText::get(line, "chipid", chipid);
        readText::get(line, "seeing", seefactor);
        readText::get(line, "wind", wind);
        readText::get(line, "winddir", winddir);
        readText::get(line, "outerscale", outerscale);
        readText::get(line, "height", height);
        readText::get(line, "rayangle", rayAngle);
        readText::get(line, "o2angle", o2Angle);
        readText::get(line, "o3angle", o3Angle);
        readText::get(line, "h2oangle", h2oAngle);
        readText::get(line, "aerosolangle", aerosolAngle);
        readText::get(line, "raygradient", raygradient);
        readText::get(line, "o2gradient", o2gradient);
        readText::get(line, "o3gradient", o3gradient);
        readText::get(line, "h2ogradient", h2ogradient);
        readText::get(line, "aerosolgradient", aerosolgradient);
        readText::get(line, "cloudmean", cloudmean);
        readText::get(line, "cloudvary", cloudvary);
        readText::get(line, "atmosphericdispersion", atmospheric_dispersion);
        readText::get(line, "atmosphericdispcenter", atmosphericdispcenter);
        readText::get(line, "seed", ranseed);
        readText::get(line, "obsseed", obsseed);//
        readText::get(line, "vistime", vistime);//
        readText::get(line, "pressure", pressure);
        readText::get(line, "waterpressure", waterPressure);
        readText::get(line, "temperature", temperature);
        readText::get(line, "tempvar", tempvar);
        readText::get(line, "airrefraction", airrefraction);
        readText::get(line, "reldensity", raynorm);
        readText::get(line, "relo2", o2norm);
        readText::get(line, "relo3", o3norm);
        readText::get(line, "relh2o", h2onorm);
        readText::get(line, "aerosoltau", aerosoltau);
        readText::get(line, "aerosolindex", aerosolindex);
        readText::get(line, "lascatprob", laScatterProb);
        readText::get(line, "domeseeing", domeseeing);
        readText::get(line, "toypsf", toypsf);
        readText::get(line, "airglowvariation", airglowvariation);//
        readText::get(line, "totalseeing", totalseeing);//
        readText::get(line, "airglowcintensity", airglowcintensity);//
        readText::get(line, "airglowpintensity", airglowpintensity);//
        readText::get(line, "domelight", domelight);//
        readText::get(line, "telconfig", telconfig);//
        readText::get(line, "checkpointcount", checkpointcount);//
        readText::get(line, "checkpointtotal", checkpointtotal);//
        readText::get(line, "domewave", domewave);//
        readText::get(line, "raydensity", raydensity);//
        readText::get(line, "scalenumber", scalenumber);//
        readText::get(line, "centerx", centerx);
        readText::get(line, "centery", centery);
        readText::get(line, "pixelsize", pixsize);
        readText::get(line, "pixelsx", pixelsx);
        readText::get(line, "pixelsy", pixelsy);
        readText::get(line, "minx", minx);
        readText::get(line, "miny", miny);
        readText::get(line, "maxx", maxx);
        readText::get(line, "maxy", maxy);
        readText::get(line, "wavelength", centralwavelength);
        readText::get(line, "platescale", platescale);
        readText::get(line, "groundlevel", groundlevel);
        readText::get(line, "xtellocation", xtelloc);
        readText::get(line, "ytellocation", ytelloc);// up to here
        if (readText::getKey(line, "latitude", latitude)) latitude *= DEGREE;
        if (readText::getKey(line, "longitude", longitude)) longitude *= DEGREE;
        if (readText::getKey(line, "pointingra", pra)) pra *= DEGREE;
        if (readText::getKey(line, "pointingdec", pdec)) pdec *= DEGREE;
        if (readText::getKey(line, "rotatex", rotatex)) rotatex *= (-DEGREE);
        if (readText::getKey(line, "rotatey", rotatey)) rotatey *= (-DEGREE);
        if (readText::getKey(line, "rotationangle", rotatez)) rotatez *= (-DEGREE);
        if (readText::getKey(line, "spiderangle", spiderangle)) spiderangle *= DEGREE;
        if (readText::getKey(line, "altitude", altitude)) {
            altitude *= DEGREE;
            zenith = PI/2.0 - altitude;
        }
        if (readText::getKey(line, "zenith", zenith)) zenith *= DEGREE;
        if (readText::getKey(line, "azimuth", azimuth)) azimuth *= DEGREE;
        if (readText::getKey(line, "moonra", moonra)) moonra *= DEGREE;
        if (readText::getKey(line, "moondec", moondec)) moondec *= DEGREE;
        if (readText::getKey(line, "moonalt", moonalt)) moonalt *= DEGREE;
        if (readText::getKey(line, "moondist", moondist)) moondist *= DEGREE;
        readText::getKey(line, "hourangle", hourangle);
        readText::getKey(line, "lst", lst);
        if (readText::getKey(line, "phaseang", phaseang)) phaseang = PI - phaseang*PI/100.0;
        if (readText::getKey(line, "solaralt", solaralt)) {
            solaralt *= DEGREE;
            solarzen = PI/2.0 - solaralt;
        }
        if (readText::getKey(line, "solarzen", solarzen)) solarzen *= DEGREE;



        if (keyName == "clearperturbations") {
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < 6; j++) {
                    body[i][j] = 0.0;
                }
            }
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < NTERM; j++) {
                    izernike[i][j] = 0.0;
                }
            }
            for (int i = 0; i < MAX_SURF; i++) {
                feaflag[i]=0;
            }
        }
        if (keyName == "cleartracking") {
            rotationjitter = 0.0;
            elevationjitter = 0.0;
            azimuthjitter = 0.0;
        }
        if (keyName == "clearclouds") {
            for (int i = 0; i < MAX_LAYER; i++) {
                cloudmean[i] = 0.0;
                cloudvary[i] = 0.0;
            }
        }
        if (keyName == "clearopacity") {
            h2onorm = 0.0;
            raynorm = 0.0;
            o2norm = 0.0;
            o3norm = 0.0;
            aerosoltau = 0.0;
            h2ogradient = 0.0;
            raygradient = 0.0;
            o2gradient = 0.0;
            o3gradient = 0.0;
            aerosolgradient = 0.0;
            for (int i = 0; i < MAX_LAYER; i++) {
                cloudmean[i] = 0.0;
                cloudvary[i] = 0.0;
            }
            opacitymode = 0;
        }
        if (keyName == "clearturbulence") {
            for (int i = 0; i < MAX_LAYER; i++) {
                seefactor[i] = 0.0;
            }
            domeseeing = 0.0;
            if (atmospheremode > 1) atmospheremode = 1;
        }
        if (keyName == "cleardefects") {
            impurityvariation = 0;
            fieldanisotropy = 0;
            deadlayer = 0;
            chargesharing = 0;
            pixelerror = 0;
        }
        if (keyName == "cleareverything") {
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < 6; j++) {
                    body[i][j] = 0.0;
                }
            }
            for (int i = 0; i < MAX_SURF; i++) {
                for (int j = 0; j < NTERM; j++) {
                    izernike[i][j] = 0.0;
                }
            }
            rotationjitter = 0.0;
            elevationjitter = 0.0;
            azimuthjitter = 0.0;
            for (int i = 0; i < MAX_LAYER; i++) {
                cloudmean[i] = 0.0;
                cloudvary[i] = 0.0;
            }
            h2onorm = 0.0;
            raynorm = 0.0;
            o2norm = 0.0;
            o3norm = 0.0;
            aerosoltau = 0.0;
            h2ogradient = 0.0;
            raygradient = 0.0;
            o2gradient = 0.0;
            o3gradient = 0.0;
            aerosolgradient = 0.0;
            for (int i = 0; i < MAX_LAYER; i++) {
                cloudmean[i] = 0.0;
                cloudvary[i] = 0.0;
            }
            for (int i = 0; i < MAX_LAYER; i++) {
                seefactor[i] = 0.0;
            }
            domeseeing = 0.0;
            impurityvariation = 0;
            fieldanisotropy = 0;
            deadlayer = 0;
            chargesharing = 0;
            photoelectric = 0;
            pixelerror = 0;
            detectorMode = 0;
            telescopeMode = 0;
            atmospheric_dispersion = 0;
            diffractionMode = 0;
            laScatterProb = 0.0;
            contaminationmode = 0;
            atmospheremode = 0;
            airrefraction = 0;
        }
        if (keyName == "surfacemap") {
            long surfaceIndex;
            iss >> surfaceIndex;
            iss >> feafile[surfaceIndex*MAX_SURF*2 + feaflag[surfaceIndex]];
            iss >> feascaling[surfaceIndex*MAX_SURF*2 + feaflag[surfaceIndex]];
            feaflag[surfaceIndex] += 1;
        }
        if (keyName == "focus") {
            focusFlag = 1;
            iss >> focusSteps;
            iss >> focusStepZ;
            focusStepZ = focusStepZ*1e-3;
            iss >> focusStepX; 
        }
        if (keyName == "distortion") {
            long surfaceIndex;
            std::ostringstream fileName1, fileName2;
            iss >> surfaceIndex;
            fileName1 << "fea_"  << obshistid << "_" << surfaceIndex << ".txt";
            feafile[surfaceIndex*MAX_SURF*2].assign(fileName1.str());
            feaflag[surfaceIndex] = 1;
            feascaling[surfaceIndex*MAX_SURF*2]=1.0;
            iss >> surfaceIndex;
            if (surfaceIndex >= 0) {
                fileName2 << "fea_"  << obshistid << "_" << surfaceIndex << ".txt";
                feafile[surfaceIndex*MAX_SURF*2].assign(fileName2.str());
                feaflag[surfaceIndex] = 1;
                feascaling[surfaceIndex*MAX_SURF*2]=1.0;
            }
        }
        if (extraCommandFlag > 1) {
            extraCommandString.push_back(line);
            extraCommandFlag++;
        }
        if (extraCommandFlag == 1) extraCommandFlag++;

    }

    // read location file
    sss = instrdir + "/location.txt";
    std::ifstream inStream(sss.c_str());
    if (inStream) {
        readText locationPars(instrdir + "/location.txt");
        for (size_t t(0); t < locationPars.getSize(); t++) {
            std::string line(locationPars[t]);
            if (groundlevel == -1) readText::get(line, "groundlevel", groundlevel);
            if (xtelloc == -1) readText::get(line, "xtellocation", xtelloc);
            if (ytelloc == -1) readText::get(line, "ytellocation", ytelloc);
            if (latitude == -1) if (readText::getKey(line, "latitude", latitude)) latitude *= DEGREE;
            if (longitude == -1) if (readText::getKey(line, "longitude", longitude)) longitude *= DEGREE;
            if (spaceMode == -1) readText::get(line, "spacemode", spaceMode);
        }
    } else {
        groundlevel = 0.0;
        xtelloc = 0.0;
        ytelloc = 0.0; 
        latitude = 0.0;
        longitude = 0.0;
    }
    if (spaceMode == -1) spaceMode = 0;

    if (spaceMode > 0) {
        atmospheric_dispersion = 0;
        //clear turbulence
        for (int i = 0; i < MAX_LAYER; i++) {
            seefactor[i] = 0.0;
        }
        domeseeing = 0.0;
        //clear opacity
        h2onorm = 0.0;
        raynorm = 0.0;
        o2norm = 0.0;
        o3norm = 0.0;
        aerosoltau = 0.0;
        h2ogradient = 0.0;
        raygradient = 0.0;
        o2gradient = 0.0;
        o3gradient = 0.0;
        aerosolgradient = 0.0;
        for (int i = 0; i < MAX_LAYER; i++) {
            cloudmean[i] = 0.0;
            cloudvary[i] = 0.0;
        }
        opacitymode = 0;
        windjitter = 0.0;
        windshake = 0.0;
        atmospheremode = 0;
    }

    if (opticsonlymode == 1) {
        detectorMode = 0;
        diffractionMode = 0;
        contaminationmode = 0;
        laScatterProb = 0.0;
        atmospheric_dispersion = 0;
        straylight = 0;
        rotationjitter = 0.0;
        elevationjitter = 0.0;
        azimuthjitter = 0.0;
        for (int i = 0; i < MAX_LAYER; i++) {
            seefactor[i] = 0.0;
        }
        domeseeing = 0.0;
        h2onorm = 0.0;
        raynorm = 0.0;
        o2norm = 0.0;
        o3norm = 0.0;
        aerosoltau = 0.0;
        h2ogradient = 0.0;
        raygradient = 0.0;
        o2gradient = 0.0;
        o3gradient = 0.0;
        aerosolgradient = 0.0;
        for (int i = 0; i < MAX_LAYER; i++) {
            cloudmean[i] = 0.0;
            cloudvary[i] = 0.0;
        }
        airrefraction = 0;
        atmospheremode = 0;
    }

    // this will signal that the atmosphere is completely off
    if (atmospheremode == 1 && atmospheric_dispersion == 0.0 && opacitymode == 0) atmospheremode = 0;

    if (atmospheremode == 0) {
        natmospherefile = 0;
        h2onorm = 0.0;
        raynorm = 0.0;
        o2norm = 0.0;
        o3norm = 0.0;
        aerosoltau = 0.0;
        h2ogradient = 0.0;
        raygradient = 0.0;
        o2gradient = 0.0;
        o3gradient = 0.0;
        aerosolgradient = 0.0;
        airrefraction = 0;
        for (int i = 0; i < MAX_LAYER; i++) {
            cloudmean[i] = 0.0;
            cloudvary[i] = 0.0;
        }
        for (int i = 0; i < MAX_LAYER; i++) {
            seefactor[i] = 0.0;
        }
    }

    if (telconfig != 2 && telconfig != 3) domelight = 1000.0;

    if (opdfile) aperturemode = 2;

    std::ostringstream outfile;
    std::ostringstream outfileevent;
    unsigned pos = instrdir.rfind("/") + 1;
    for (unsigned i = pos; i < instrdir.length(); i++) {
        outfile << instrdir[i];
        outfileevent << instrdir[i];
    }
    outfile << "_e_"  << obshistid << "_f"<< filter << "_" << chipid << "_E" << std::setfill('0') << std::setw(3) << pairid;
    outputfilename = outfile.str();
    outfileevent << "_r_"  << obshistid << "_f"<< filter << "_" << chipid << "_E" << std::setfill('0') << std::setw(3) << pairid;
    eventFitsFileName = outfileevent.str();

    if (flatdir == 1) {
        instrdir  =  ".";
        bindir  =  ".";
    }

    if (tarfile == 1) {
        std::ostringstream tarName;
        tarName << "raytrace_" << obshistid << ".tar";
        std::ifstream tarFile(tarName.str().c_str());
        if (tarFile.good()) {
            std::cout << "Untarring " << tarName.str()<<std::endl;
            std::string tarCommand = "tar xf " + tarName.str();
            system(tarCommand.c_str());
        }
    }

    // OPD
    focalplanefile = instrdir + "/focalplanelayout.txt";
    readText focalplanePar(focalplanefile);
    std::string tchipid;
    if (chipid == "opd") {
        std::string line(focalplanePar[focalplanePar.getSize()/2]);
        std::istringstream iss(line);
        iss >> tchipid;
        opdfile = 1;
        transtol = 1.0;
        aperturemode = 2;
        detectorMode = 0;
        diffractionMode = 0;
        coatingmode = 0;
        contaminationmode = 0;
        laScatterProb = 0.0;
        atmospheric_dispersion = 0;
        straylight = 0;
        airrefraction = 0;
        rotationjitter = 0.0;
        elevationjitter = 0.0;
        azimuthjitter = 0.0;
        backgroundMode = 0;
        for (int i = 0; i < MAX_LAYER; i++) {
            seefactor[i] = 0.0;
        }
        domeseeing = 0.0;
        h2onorm = 0.0;
        raynorm = 0.0;
        o2norm = 0.0;
        o3norm = 0.0;
        aerosoltau = 0.0;
        h2ogradient = 0.0;
        raygradient = 0.0;
        o2gradient = 0.0;
        o3gradient = 0.0;
        aerosolgradient = 0.0;
        for (int i = 0; i < MAX_LAYER; i++) {
            cloudmean[i] = 0.0;
            cloudvary[i] = 0.0;
        }
    } else {
        tchipid = chipid;
    }
    std::istringstream focalplanePars(readText::get(focalplanefile, tchipid));
    double centerxT, centeryT, pixsizeT;
    long pixelsxT, pixelsyT;
    double angle1, angle2;
    float sensorthicknessT;
    std::string grouptype;
    focalplanePars >> centerxT >> centeryT >> pixsizeT >> pixelsxT >> pixelsyT >>
        devmaterial >> devtype >> devmode >> devvalue >> sensorthicknessT >> grouptype >> chipangle >> angle1 >> angle2 >> decenterx >> decentery;
    decenterx *= 1000.0;
    decentery *= 1000.0;
    chipangle *= PI/180.0;
    if (centerx == -1) centerx = centerxT;
    if (centery == -1) centery = centeryT;
    if (pixsize == -1) pixsize = pixsizeT;
    if (pixelsx == -1) pixelsx = pixelsxT;
    if (pixelsy == -1) pixelsy = pixelsyT;
    if (minx == -1) minx = 0;
    if (miny == -1) miny = 0;
    if (maxx == -1) maxx = pixelsx - 1;
    if (maxy == -1) maxy = pixelsy - 1;
    if (sensorthickness == -1) sensorthickness = sensorthicknessT;
    int nsamples = nframes + nskip;
    if (nframes == 1 && nskip == 0) devmode = "frame"; //override
    if (devmode == "frame") {
        ngroups = 1;
        nframes = 1;
        nskip = 0;
        nsamples = 1;
    } else {
        ngroups = ceil((vistime + devvalue*nskip)/(devvalue*nsamples));
    }

    // read sensor file
    sss = instrdir + "/sensor.txt";
    std::ifstream inStream2(sss.c_str());
    if (inStream2) {
        readText siliconPars(instrdir + "/sensor.txt");
        for (size_t t(0); t < siliconPars.getSize(); t++) {
            std::string line(siliconPars[t]);
            if (well_depth == -1) readText::get(line, "wellDepth", well_depth);
            if (nbulk == -1) readText::get(line, "nbulk", nbulk);
            if (nf == -1) readText::get(line, "nf", nf);
            if (nb == -1) readText::get(line, "nb", nb);
            if (sf == -1) readText::get(line, "sf", sf);
            if (sb == -1) readText::get(line, "sb", sb);
            if (sensorTempNominal == -1) readText::get(line, "sensorTempNominal", sensorTempNominal);
        }
    } else {
        well_depth = 1e5; 
        nbulk = 1e12;
        nf = 0.0;
        nb = 0.0;
        sf = 0.0;
        sb = 0.0;
        sensorTempNominal = 173;
    }

<<<<<<< HEAD
    // scaling background parameters
    backRadius = backRadius*180000.0/platescale*pixsize/10.0;
    backBuffer = backBuffer*platescale/180000.0*10.0/pixsize;
    activeBuffer = activeBuffer*platescale/180000.0*10.0/pixsize;
    backBeta = backBeta*180000.0/platescale*pixsize/10.0;
    if (backGamma > 1) backGamma = backGamma*180000.0/platescale*pixsize/10.0;

    windjitter = windjitter*pow(vistime/60, 0.25);
    focusStepX = focusStepX * platescale / 3600.0 * 1e-3;


=======
    sss = instrdir + "/tracking.txt";
    std::ifstream inStream3(sss.c_str());
    if (inStream3) {
        readText trackingPars(instrdir + "/tracking.txt");
        for (size_t t(0); t < trackingPars.getSize(); t++) {
            std::string line(trackingPars[t]);
            if (windjitter == -1) readText::get(line, "windjitter", windjitter);
            if (rotationjitter == -1) readText::get(line, "rotationjitter", rotationjitter);
            if (elevationjitter == -1) readText::get(line, "elevationjitter", elevationjitter);
            if (azimuthjitter == -1) readText::get(line, "azimuthjitter", azimuthjitter);
        }
    } else {
        windjitter = 2.5;
        rotationjitter = 0.0;
        elevationjitter = 0.0;
        azimuthjitter = 0.0;
    }
    windjitter = windjitter * pow(vistime/60, 0.25);
>>>>>>> 83f06c65... # This is a combination of 5 commits.
    return(0);

}

int Observation::settings() {

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Basic Setup" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "[outputdir] Output directory:                              " << outputdir << std::endl;
    std::cout << "[outputfilename] Output file name:                         " << outputfilename << std::endl;
    std::cout << "[seddir] SED directory:                                    " << seddir << std::endl;
    std::cout << "[imagedir] Image directory:                                " << imagedir << std::endl;
    std::cout << "[centroidfile] Output centroid file (0=no/1=yes):          " << centroidfile << std::endl;
    std::cout << "[throughputfile] Output throughput file (0=no/1=yes):      " << throughputfile << std::endl;
    std::cout << "[eventfile] Output event file (0=no/1=yes):                " << eventfile << std::endl;
    std::cout << "[eventFitsFileName] Output event Fits file name:           " << eventFitsFileName << std::endl;
    std::cout << "[bindir] Binary directory:                                 " << bindir << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Module Switches" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "[backgroundmode] Background mode (0=off/1=on):             " << backgroundMode << std::endl;
    std::cout << "[telescopemode] Telescope mode (0=off/1=on):               " << telescopeMode << std::endl;
    std::cout << "[trackingmode] Tracking mode (0=off/1=on):                 " << trackingMode << std::endl;
    std::cout << "[windshakemode] Wind shake mode (0=off/1=on):              " << windShakeMode << std::endl;
    std::cout << "[detectormode] Detector mode (0=off/1=on):                 " << detectorMode << std::endl;
    std::cout << "[diffractionmode] Diffraction mode (0=off/1=on):           " << diffractionMode << std::endl;
    std::cout << "[spacemode] Space mode (0=off/1=LEO/2=L1/etc.):            " << spaceMode << std::endl;
    std::cout << "[zernikemode] Zernike mode (0=off/1=on):                   " << zernikemode << std::endl;
    std::cout << "[straylight] Straylight mode (0=off/1=on):                 " << straylight << std::endl;
    std::cout << "[aperturemode] Aperture mode (0=normal/1=on):              " << aperturemode << std::endl;
    std::cout << "[ghostonly] Ghost-only mode (0=normal/1=on):               " << ghostonly << std::endl;
    std::cout << "[saturation] Saturation mode (0=off/1=on):                 " << saturation << std::endl;
    std::cout << "[blooming] Blooming mode (0=off/1=on):                     " << blooming << std::endl;
    std::cout << "[atmosphericdispersion] Atmos. dispersion (0=off/1=on):    " << atmospheric_dispersion << std::endl;
    std::cout << "[atmosphericdispcenter] Atmos. disp. ctr. corr.:           " << atmosphericdispcenter << std::endl;
    std::cout << "[impurityvariation] Impurity variation (0=off/1=on):       " << impurityvariation << std::endl;
    std::cout << "[fieldanisotropy] Field anisotropy (0=off/1=on):           " << fieldanisotropy << std::endl;
    std::cout << "[fringing] Fringing (0=off/1=on):                          " << fringeflag << std::endl;
    std::cout << "[deadlayer] Dead layer (0=off/1=on):                       " << deadlayer << std::endl;
    std::cout << "[chargediffusion] Charge diffusion (0=off/1=on):           " << chargediffusion << std::endl;
    std::cout << "[photoelectric] Photoelectric (0=off/1=on):                " << photoelectric << std::endl;
    std::cout << "[chargesharing] Charge sharing (0=off/1=on):               " << chargesharing << std::endl;
    std::cout << "[pixelerror] Pixel error (0=off/1=on):                     " << pixelerror << std::endl;
    std::cout << "[coatingmode] Coating mode (0=off/1=on):                   " << coatingmode << std::endl;
    std::cout << "[contaminationmode] Contamination mode (0=off/1=on):       " << contaminationmode << std::endl;
    std::cout << "[aberration] Aberration mode (0=off/1=on):                 " << aberration << std::endl;
    std::cout << "[precession] Precession mode (0=off/1=on):                 " << precession << std::endl;
    std::cout << "[nutation] Nutation mode (0=off/1=on):                     " << nutation << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Telescope Operator and Bookkeeping" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Pointing right ascension (degrees):                        " << pra/DEGREE << std::endl;
    std::cout << "Pointing declination (degrees):                            " << pdec/DEGREE << std::endl;
    std::cout << "[rotationangle] Rotation angle (rotSkyPos) (degrees):      " << -rotatez/DEGREE << std::endl;
    std::cout << "Angle of spider (rotTelPos) (degrees):                     " << spiderangle/DEGREE << std::endl;
    if (spaceMode == 0) {
        std::cout << "Zenith angle (degrees):                                    " << zenith/DEGREE << std::endl;
        std::cout << "Azimuthal angle (degrees):                                 " << azimuth/DEGREE << std::endl;
    }
    std::cout << "Filter (number starting with 0):                           " << filter << std::endl;
    std::cout << "Filter Name:                                               " << filterName << std::endl;
    std::cout << "Random seed:                                               " << ranseed << std::endl;
    std::cout << "Sensor temperature (K):                                    " << (sensorTempNominal + sensorTempDelta) << std::endl;
    std::cout << "[sensortempdelta] Delta sensor temperature (K):            " << sensorTempDelta << std::endl;


    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Instantaneous Instrument and Site Characteristics" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "[instrdir] Instrument & site directory:                    " << instrdir << std::endl;

    //instrument
    std::cout << "[platescale] Plate scale:                                  " << platescale << std::endl;
    std::cout << "[minr] Minimum aperture radius:                            " << minr << std::endl;
    std::cout << "[maxr] Maximum aperture radius:                            " << maxr << std::endl;
    std::cout << "[chipid] Chip/Amplifier ID:                                " << chipid << std::endl;
    std::cout << "[centerx] Chip center x (microns):                         " << centerx << std::endl;
    std::cout << "[centery] Chip center y (microns):                         " << centery << std::endl;
    std::cout << "[pixelsx] Chip x pixels:                                   " << pixelsx << std::endl;
    std::cout << "[pixelsy] Chip y pixels:                                   " << pixelsy << std::endl;
    std::cout << "[minx] Minimum x pixel of amplifier:                       " << minx << std::endl;
    std::cout << "[maxx] Maximum x pixel of amplifier:                       " << maxx << std::endl;
    std::cout << "[miny] Minimum y pixel of amplifier:                       " << miny << std::endl;
    std::cout << "[maxy] Maximum y pixel of amplifier:                       " << maxy << std::endl;
    std::cout << "[pixelsize] Pixel Size (microns):                          " << pixsize << std::endl;
    std::cout << "[welldepth] Full well depth:                               " << well_depth << std::endl;
    std::cout << "[nbulk] Bulk doping density:                               " << std::scientific << nbulk << std::endl;
    std::cout << "[nf] Front side doping density:                            " << std::scientific << nf << std::endl;
    std::cout << "[nb] Back side doping density:                             " << std::scientific << nb << std::endl;
    std::cout << "[sf] Front side doping scale:                              " << std::scientific << sf << std::endl;
    std::cout << "[sb] Back side doping scale:                               " << std::scientific << sb << std::endl;
    std::cout << "[sensorthickness] Sensor Thickness (microns):              " << std::fixed << sensorthickness << std::endl;
    std::cout << "[overdepbias] Over depletion bias (volts):                 " << overdepbias << std::endl;
    std::cout << "[sensorTempNominal] Nominal sensor temperature (K):        " << sensorTempNominal << std::endl;
    std::cout << "[qevariation] QE variation:                                " << qevariation << std::endl;
    std::cout << "[exptime] Exposure time (s):                               " << exptime << std::endl;
    std::cout << "[nsnap] Number of snaps:                                   " << nsnap << std::endl;
    std::cout << "Number of groups:                                          " << ngroups << std::endl;
    std::cout << "[nframes] Number of frames:                                " << nframes << std::endl;
    std::cout << "[nskip] Number of frames to skip:                          " << nskip << std::endl;
    std::cout << "[shuttererror] Shutter error (s):                          " << shuttererror << std::endl;
    std::cout << "[timeoffset] Time offset (s):                              " << timeoffset << std::endl;
    std::cout << "[windjitter] Wind jitter (degrees):                        " << windjitter << std::endl;
    std::cout << "[rotationjitter] Rotation jitter (arcseconds):             " << rotationjitter << std::endl;
    std::cout << "[elevationjitter] Elevation jitter (arcseconds):           " << elevationjitter << std::endl;
    std::cout << "[azimuthjitter] Azimuthal jitter (arcseconds):             " << azimuthjitter << std::endl;
    std::cout << "[rotatex] Fiducial field angle X (degrees):                " << -rotatex/DEGREE << std::endl;
    std::cout << "[rotatey] Fiducial field angle Y (degrees):                " << -rotatey/DEGREE << std::endl;
    std::cout << "[izernike optic# zernike#] Zernike amplitude:              " << std::endl;
    for (long i = 0; i < nsurf + 1; i++) {
        for (long j = 0; j < NTERM/2; j++) {
            std::cout << izernike[i][j] << " ";
        }
        std::cout << std::endl;
        std::cout << " ";
        for (long j = NTERM/2; j < NTERM; j++) {
            std::cout << izernike[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "[body optic# dof#] Body motion of optics:                  " << std::endl;
    for (long i = 0; i < nsurf + 1; i++) {
        for (long j = 0; j < 6; j++) {
            std::cout <<  body[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "[lascatprob] Large angle scattering probability:           " << laScatterProb << std::endl;
    std::cout << "[toypsf] Toy PSF:                                          " << toypsf << std::endl;

    //site
    if (spaceMode == 0) {
        std::cout << "[domeseeing] Dome seeing:                                  " << domeseeing << std::endl;
        std::cout << "[groundlevel] Ground level (m):                            " << groundlevel << std::endl;
        std::cout << "[xtellocation] X Telescope location (m):                   " << xtelloc << std::endl;
        std::cout << "[ytellocation] Y Telescope location (m):                   " << ytelloc << std::endl;
        std::cout << "[latitude] Telescope latitude (degrees):                   " << latitude/DEGREE << std::endl;
        std::cout << "[longitude] Telescope longitude (degrees):                 " << longitude/DEGREE << std::endl;
    } else if (spaceMode == 1) {
        std::cout << "Location:                                                  low Earth orbit" << std::endl;
    } else if (spaceMode == 2) {
        std::cout << "Location:                                                  Lagrange point 1" << std::endl;
    } else if (spaceMode == 3) {
        std::cout << "Location:                                                  Lagrange point 2" << std::endl;
    } else if (spaceMode == 4) {
        std::cout << "Location:                                                  Lagrange point 3" << std::endl;
    } else if (spaceMode == 5) {
        std::cout << "Location:                                                  Lagrange point 4" << std::endl;
    } else if (spaceMode == 6) {
        std::cout << "Location:                                                  Lagrange point 5" << std::endl;
    }
    if (natmospherefile > 0) {
        std::cout << "[pressure] Air pressure (mmHg):                            " << pressure << std::endl;
        std::cout << "[waterpressure] Water vapor pressure (mmHg):               " << waterPressure << std::endl;
        std::cout << "[temperature] Ground temperature (degrees C):              " << temperature << std::endl;
        std::cout << "[tempvar] Temperature derivative (degrees per hour):       " << tempvar << std::endl;
        std::cout << "[reldensity] Relative density:                             " << raynorm << std::endl;
        std::cout << "[relo2] Relative O2 fraction:                              " << o2norm << std::endl;
        std::cout << "[relh2o] Relative H2O fraction:                            " << h2onorm << std::endl;
        std::cout << "[aerosoltau] Aerosol optical depth:                        " << aerosoltau << std::endl;
        std::cout << "[aerosolindex] Aerosol index:                              " << aerosolindex << std::endl;
        std::cout << "[relo3] Relative O3 fraction:                              " << o3norm << std::endl;
        std::cout << "[raygradient] Density gradient (fraction/km):              " << raygradient << std::endl;
        std::cout << "[o2gradient] O2 gradient (fraction/km):                    " << o2gradient << std::endl;
        std::cout << "[o3gradient] O3 gradient (fraction/km):                    " << o3gradient << std::endl;
        std::cout << "[h2ogradient] H2O gradient (fraction/km):                  " << h2ogradient << std::endl;
        std::cout << "[aerosolgradient] Aerosol gradient:                        " << aerosolgradient << std::endl;
        std::cout << "[rayangle] Density angle:                                  " << rayAngle << std::endl;
        std::cout << "[o2angle] O2 angle:                                        " << o2Angle << std::endl;
        std::cout << "[o3angle] O3 angle:                                        " << o3Angle << std::endl;
        std::cout << "[h2oangle] H2O angle:                                      " << h2oAngle << std::endl;
        std::cout << "[aerosolangle] Aerosol angle:                              " << aerosolAngle << std::endl;
    }
    std::cout << "[natmospherefile] Number of atmosphere layers:             " << natmospherefile << std::endl;
    if (natmospherefile > 0) {
        std::cout << "[seeing layer#] Seeing at 5000 angstroms (arcsec):" << std::endl;
        for (long i = 0; i < natmospherefile; i++) {
            std::cout << seefactor[i]/pow(1/cos(zenith), 0.6) << " ";
        }
        std::cout << std::endl;
        std::cout << "[wind layer#] Wind speed (m/s):" << std::endl;
        for (long i = 0; i < natmospherefile; i++) {
            std::cout << wind[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "[winddir layer#] Wind direction (degrees):" << std::endl;
        for (long i = 0; i < natmospherefile; i++) {
            std::cout << winddir[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "[height layer#] Layer height (km):" << std::endl;
        for (long i = 0; i < natmospherefile; i++) {
            std::cout << height[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "[outerscale layer#] Outer scale (m):" << std::endl;
        for (long i = 0; i < natmospherefile; i++) {
            std::cout << outerscale[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "[cloudmean layer#] Mean cloud extinction (km):" << std::endl;
        for (long i = 0; i < natmospherefile; i++) {
            std::cout << cloudmean[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "[cloudvary layer#] Variation of cloud extinction (km):" << std::endl;
        for (long i = 0; i < natmospherefile; i++) {
            std::cout << cloudvary[i] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;

    return(1);

}

int Observation::header(fitsfile *faptr) {

    long i;
    int status = 0;
    char tempstring[4096];
    char tempstring2[4096];
    double tempf1;

    for (long i = 0; i < static_cast<long>(extraCommandString.size()); i++) {
        sprintf(tempstring, "PHOV%ld", i);
        sprintf(tempstring2, "Physics Override Command %ld", i);
        fits_write_key(faptr, TSTRING, (char*)tempstring, (char*)extraCommandString[i].c_str(), (char*)tempstring2, &status);
    }
    fitsWriteKey(faptr, "OBSID", obshistid, "Observation ID");
    fitsWriteKey(faptr, "TAI", convertTai(tai + (timeoffset - exptime/2.0)/24.0/3600.0), "International Atomic Time scale");
    fitsWriteKey(faptr, "MJD-OBS", tai + (timeoffset - exptime/2.0)/24.0/3600.0, "Modified Julian date");
    fitsWriteKey(faptr, "HOURANG", hourangle, "Hour angle");
    fitsWriteKey(faptr, "LST", lst, "Local sidereal time");
    fitsWriteKey(faptr, "OUTPDIR", outputdir, "Output directory");
    fitsWriteKey(faptr, "OUTFILE", outputfilename, "Output file name");
    fitsWriteKey(faptr, "SEDDIR", seddir, "SED directory");
    fitsWriteKey(faptr, "IMGDIR", imagedir, "Image directory");
    fitsWriteKey(faptr, "BACMODE", backgroundMode, "Background mode (0=off/1=on)");
    fitsWriteKey(faptr, "TELMODE", telescopeMode, "Telescope mode (0=off/1=on)");
    fitsWriteKey(faptr, "SPCMODE", spaceMode, "Space mode (0=off)");
    fitsWriteKey(faptr, "TRKMODE", trackingMode, "Tracking mode (0=off/1=on)");
    fitsWriteKey(faptr, "WSKMODE", windShakeMode, "Wind shake mode (0=off/1=on)");
    fitsWriteKey(faptr, "DIFMODE", diffractionMode, "Diffraction mode (0=off/1=on)");
    fitsWriteKey(faptr, "DETMODE", detectorMode, "Detector mode (0=off/1=on)");
    fitsWriteKey(faptr, "ZERMODE", zernikemode, "Zernike mode (0=off/1=on)");
    fitsWriteKey(faptr, "STRLGHT", straylight, "Straylight mode (0=off/1=on)");
    fitsWriteKey(faptr, "APRMODE", aperturemode, "Aperture mode (0=normal/1=on)");
    fitsWriteKey(faptr, "GHOMODE", ghostonly, "Ghost-only mode (0=normal/1=on)");
    fitsWriteKey(faptr, "SATMODE", saturation, "Saturation mode (0=off/1=on)");
    fitsWriteKey(faptr, "BLOOMNG", blooming, "Blooming mode (0=off/1=on)");
    fitsWriteKey(faptr, "EVTFILE", eventfile, "Output event file (0=no/1=yes)");
    fitsWriteKey(faptr, "MOONDIS", moondist, "Distance to moon");
    fitsWriteKey(faptr, "MOONPHS", phaseang, "Phase of the Moon");
    fitsWriteKey(faptr, "MOONALT", moonalt, "Moon altitude (radians)");
    fitsWriteKey(faptr, "SUNALT", solaralt, "Sun altitude (radians)");

    fits_write_key(faptr, TLONG, (char*)"THRFILE", &throughputfile, (char*)"Output throughput file (0=no/1=yes)", &status);
    tempf1 = (3600.0*1000.0)/platescale;
    fits_write_key(faptr, TDOUBLE, (char*)"PLTSCAL", &tempf1, (char*)"Approx. Plate scale (arcsec/mm)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"MINR", &minr, (char*)"Minimum aperture radius", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"MAXR", &maxr, (char*)"Maximum aperture radius", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"EXPTIME", &exptime, (char*)"Exposure time", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DARKTIME", &exptime, (char*)"Actual exposed time", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SHUTERR", &shuttererror, (char*)"Shutter error", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TIMEOFF", &timeoffset, (char*)"Time offset", &status);
    fits_write_key(faptr, TLONG, (char*)"FILTNM", &filter, (char*)"Filter/optics configuration number", &status);
    fits_write_key(faptr, TSTRING, (char*)"FILTER", &filterName, (char*)"Filter", &status);
    fits_write_key(faptr, TLONG, (char*)"SEED", &ranseed, (char*)"Random seed", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRA", &pra, (char*)"Pointing RA (radians)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PDEC", &pdec, (char*)"Pointing Dec (radians)", &status);
    tempf1 = pra*180/PI;
    fits_write_key(faptr, TDOUBLE, (char*)"RA_DEG", &tempf1, (char*)"Pointing RA (decimal degrees)", &status);
    tempf1 = pdec*180/PI;
    fits_write_key(faptr, TDOUBLE, (char*)"DEC_DEG", &tempf1, (char*)"Pointing Dec (decimal degrees)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AIRMASS", &airmass, (char*)"Airmass", &status);
    tempf1 = -rotatex*180/PI;
    fits_write_key(faptr, TDOUBLE, (char*)"FLDANGX", &tempf1, (char*)"Fiducial field angle X (degrees)", &status);
    tempf1 = -rotatey*180/PI;
    fits_write_key(faptr, TDOUBLE, (char*)"FLDANGY", &tempf1, (char*)"Fiducial field angle Y (degrees)", &status);
    tempf1 = -rotatez*180/PI;
    fits_write_key(faptr, TDOUBLE, (char*)"ROTANG", &tempf1, (char*)"Rotation angle (rotSkyPos) (degrees)", &status);
    tempf1 = spiderangle*180/PI;
    fits_write_key(faptr, TDOUBLE, (char*)"SPIDANG", &tempf1, (char*)"Angle of spider (rotTelPos)", &status);
    tempf1 = zenith*180/PI;
    fits_write_key(faptr, TDOUBLE, (char*)"ZENITH", &tempf1, (char*)"Zenith angle (degrees)", &status);
    tempf1 = azimuth*180/PI;
    fits_write_key(faptr, TDOUBLE, (char*)"AZIMUTH", &tempf1, (char*)"Azimuthal angle (degrees)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"ROTJITT", &rotationjitter, (char*)"Rotation jitter (arcseconds)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"ELEJITT", &elevationjitter, (char*)"Elevation jitter (arcseconds)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AZIJITT", &azimuthjitter, (char*)"Azimuthal jitter (arcseconds)", &status);
    fits_write_key(faptr, TSTRING, (char*)"TRKFILE", (char*)trackingfile.c_str(), (char*)"Tracking file", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"GNDLEVL", &groundlevel, (char*)"Ground level (m)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"XTELLOC", &xtelloc, (char*)"X telescope location (m)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"YTELLOC", &ytelloc, (char*)"Y telescope location (m)", &status);

    for (long i = 0; i < nsurf + 1; i++) {
        for (long j = 0; j < NZERN; j++) {
            sprintf(tempstring, "ZER%2ld%2ld", i, j);
            fits_write_key(faptr, TDOUBLE, tempstring, &izernike[i][j], (char*)"Zernike amplitude", &status);
        }
    }
    for (long i = 0; i < nsurf + 1; i++) {
        for (long j = 0; j < 6; j++) {
            sprintf(tempstring, "BOD%2ld%2ld", i, j);
            fits_write_key(faptr, TDOUBLE, tempstring, &izernike[i][j], (char*)"Body motion misalignment", &status);
        }
    }

    fits_write_key(faptr, TLONG, (char*)"ATMFILE", &natmospherefile, (char*)"Number of atmosphere files", &status);
    for (i = 0; i < natmospherefile; i++) {
        sprintf(tempstring, "AFILE%ld", i);
        fits_write_key(faptr, TSTRING, tempstring, (char*)atmospherefile[i].c_str(), (char*)"Atmosphere file", &status);
        sprintf(tempstring, "CFILE%ld", i);
        fits_write_key(faptr, TSTRING, tempstring, (char*)cloudfile[i].c_str(), (char*)"Cloud file", &status);
        sprintf(tempstring, "SEE%ld", i);
        tempf1 = seefactor[i]/(pow(1/cos(zenith), 0.6));
        fits_write_key(faptr, TDOUBLE, tempstring, &tempf1, (char*)"Seeing at 5000 angstrom (sigma)", &status);
        sprintf(tempstring, "WIND%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &wind[i], (char*)"Wind speed (m/s)", &status);
        sprintf(tempstring, "WDIR%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &winddir[i], (char*)"Wind direction (degrees)", &status);
        sprintf(tempstring, "OSCL%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &outerscale[i], (char*)"Outer scale (m)", &status);
        sprintf(tempstring, "HGT%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &height[i], (char*)"Height (km)", &status);
        sprintf(tempstring, "CMEAN%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &cloudmean[i], (char*)"Mean cloud extinction (mag)", &status);
        sprintf(tempstring, "CVARY%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &cloudvary[i], (char*)"Variation of cloud ext. (mag)", &status);
    }
    fits_write_key(faptr, TDOUBLE, (char*)"RAYGRAD", &raygradient, (char*)"Density gradient", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"O2GRAD", &o2gradient, (char*)"O2 gradient", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"O3GRAD", &o3gradient, (char*)"O3 gradient", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"H2OGRAD", &h2ogradient, (char*)"H2O gradient", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AERGRAD", &aerosolgradient, (char*)"Aerosol gradient", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"RAYANG", &rayAngle, (char*)"Density angle", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"O2ANG", &o2Angle, (char*)"O2 angle", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"O3ANG", &o3Angle, (char*)"O3 angle", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"H2OANG", &h2oAngle, (char*)"H2O angle", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AERANG", &aerosolAngle, (char*)"Aerosol angle", &status);
    fits_write_key(faptr, TLONG, (char*)"ATMDISP", &atmospheric_dispersion, (char*)"Atmos. dispersion (0=off/1=on)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRESS", &pressure, (char*)"Air pressure (mmHg)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"H2OPRESS", &waterPressure, (char*)"Water vapor pressure (mmHg)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TEMPERA", &temperature, (char*)"Ground temperature (degrees C)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TEMPVAR", &tempvar, (char*)"Temperature derivative (degrees per hour)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"RELDENS", &raynorm, (char*)"Relative density", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"RELO2", &o2norm, (char*)"Relative O2 fraction", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"RELH2O", &h2onorm, (char*)"Relative H2O fraction", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AERTAU", &aerosoltau, (char*)"Aerosol optical depth", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AERIND", &aerosolindex, (char*)"Aerosol index", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"RELO3", &o3norm, (char*)"Relative O3 fraction", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"LASCPR", &laScatterProb, (char*)"Large angle scattering probability", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DOMESEE", &domeseeing, (char*)"Dome Seeing (arcseconds FWHM)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TOYPSF", &toypsf, (char*)"Toy PSF (arcseconds FWHM)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PIXSIZE", &pixsize, (char*)"Pixel Size (microns)", &status);
    fits_write_key(faptr, TSTRING, (char*)"CHIPID", (char*)chipid.c_str(), (char*)"Chip/Amplifier ID", &status);
    fits_write_key(faptr, TLONG, (char*)"PAIRID", &pairid, (char*)"Pair ID", &status);
    fits_write_key(faptr, TSTRING, (char*)"FPFILE", (char*)focalplanefile.c_str(), (char*)"Focal plane file name", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"CENTX", &centerx, (char*)"Chip center x (microns)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"CENTY", &centery, (char*)"Chip center y (microns)", &status);
    fits_write_key(faptr, TLONG, (char*)"PIXX", &pixelsx, (char*)"Chip x pixels", &status);
    fits_write_key(faptr, TLONG, (char*)"PIXY", &pixelsy, (char*)"Chip y pixels", &status);
    fits_write_key(faptr, TLONG, (char*)"MINX", &minx, (char*)"Minimum x pixel of amplifier", &status);
    fits_write_key(faptr, TLONG, (char*)"MAXX", &maxx, (char*)"Maximum x pixel of amplifier", &status);
    fits_write_key(faptr, TLONG, (char*)"MINY", &miny, (char*)"Minimum y pixel of amplifier", &status);
    fits_write_key(faptr, TLONG, (char*)"MAXY", &maxy, (char*)"Maximum y pixel of amplifier", &status);
    fits_write_key(faptr, TLONG, (char*)"WELDPT", &well_depth, (char*)"Full well depth (electrons)", &status);
    fits_write_key(faptr, TFLOAT, (char*)"NBULK", &nbulk, (char*)"Bulk doping density", &status);
    fits_write_key(faptr, TFLOAT, (char*)"NF", &nf, (char*)"Front side doping density", &status);
    fits_write_key(faptr, TFLOAT, (char*)"NB", &nb, (char*)"Back side doping density", &status);
    fits_write_key(faptr, TFLOAT, (char*)"SF", &sf, (char*)"Front side doping scale", &status);
    fits_write_key(faptr, TFLOAT, (char*)"SB", &sb, (char*)"Back side doping scale", &status);
    fits_write_key(faptr, TFLOAT, (char*)"SETHICK", &sensorthickness, (char*)"Sensor thickness (microns)", &status);
    fits_write_key(faptr, TFLOAT, (char*)"OVRDEP", &overdepbias, (char*)"Over depletion bias (volts)", &status);
    tempf1 = sensorTempNominal + sensorTempDelta;
    fits_write_key(faptr, TDOUBLE, (char*)"CCDTEMP", &tempf1, (char*)"Sensor temperature (K):", &status);
    fits_write_key(faptr, TFLOAT, (char*)"TRX0", &impurityX, (char*)"Tree ring center", &status);
    fits_write_key(faptr, TFLOAT, (char*)"TRY0", &impurityY, (char*)"Tree ring center", &status);

    return(0);

}
