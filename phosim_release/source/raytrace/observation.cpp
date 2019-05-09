///
/// @package phosim
/// @file observation.cpp
/// @brief observation class for raytrace code
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Alan Meert (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <iostream>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <fitsio2.h>

#include "raytrace.h"
#include "observation.h"
#include "constants.h"
#include "helpers.h"
#include "constants.h"
#include "ancillary/readtext.h"
#include "coordinate.cpp"
#include "settings.cpp"
#include "sed.cpp"

int Observation::addSource (const std::string & object, int sourcetype) {

    char tempstring[4096];
    int nspatialpar, ndustpar, ndustparz;
    long i, j;
    int status;
    char *ffptr;
    fitsfile *faptr;
    long naxes[2];
    int nfound;
    int anynull;
    float nullval;
    double x, y;
    double nn;
    double mag;
    double ra, dec, redshift, gamma1, gamma2, kappa, magnification,
        deltara, deltadec;
    std::string id;

    if (nsource >= MAX_SOURCE) {
        std::perror("Too many sources in catalog");
    }

    nspatialpar = 0;

    std::istringstream iss(object);
    iss >> id >> ra >> dec;
    sources.id.push_back(id);
    sources.split.push_back(0);
    sources.ra.push_back(ra);
    sources.dec.push_back(dec);
    sources.ra[nsource] *= DEGREE;
    sources.dec[nsource] *= DEGREE;
    iss >> mag >> sources.sedfilename[nsource] >> redshift >> gamma1 >> gamma2 >> kappa >> deltara >> deltadec;
    sources.redshift.push_back(redshift);
    sources.gamma1.push_back(gamma1);
    sources.gamma2.push_back(gamma2);
    sources.kappa.push_back(kappa);
    sources.deltara.push_back(deltara);
    sources.deltadec.push_back(deltadec);
    sources.deltara[nsource] *= DEGREE;
    sources.deltadec[nsource] *= DEGREE;
    sources.type[nsource] = sourcetype;
    magnification = 2.5*log10((1 - sources.kappa[nsource])*(1 - sources.kappa[nsource])-
                              sources.gamma1[nsource]*sources.gamma1[nsource]-
                              sources.gamma2[nsource]*sources.gamma2[nsource]);
    sources.norm[nsource] = pow(10.0, ((mag + magnification + 48.6)/(-2.5)));
    sources.mag[nsource] = mag + magnification;

    // read SED file
    if (nsource > 0) {
        for (i = 0; i < nsource; i++) {
            if (sources.sedfilename[i] == sources.sedfilename[nsource]) {
                sources.sedptr[nsource] = sources.sedptr[i];
                goto skipsedread;
            }
        }
    }
    sources.sedptr[nsource] = nsedptr;
    readSed(sources.sedfilename[nsource], 0);

skipsedread:;
    sources.norm[nsource] = sources.norm[nsource]/(normwave)*(1 + sources.redshift[nsource])*sedDwdp[sources.sedptr[nsource]];

    iss >> sources.spatialname[nsource];

    if (sources.spatialname[nsource] == "point") {
        sources.spatialtype[nsource] = POINT;
        nspatialpar = 0;
    }
    if (sources.spatialname[nsource].find("fit") != std::string::npos) {
        sources.spatialtype[nsource] = IMAGE;
        nspatialpar = 2;
    }
    if (sources.spatialname[nsource] == "gauss") {
        sources.spatialtype[nsource] = GAUSSIAN;
        nspatialpar = 1;
    }
    if (sources.spatialname[nsource] == "sersic") {
        sources.spatialtype[nsource] = SERSIC;
        nspatialpar = 6;
    }
    if (sources.spatialname[nsource] == "sersic2d") {
        sources.spatialtype[nsource] = SERSIC2D;
        nspatialpar = 4;
    }
    if (sources.spatialname[nsource] == "sersic2D") {
        sources.spatialtype[nsource] = SERSIC2D;
        nspatialpar = 4;
    }
    if (sources.spatialname[nsource] == "sersicComplex") {
        sources.spatialtype[nsource] = SERSICCOMPLEX;
        nspatialpar = 14;
    }
    if (sources.spatialname[nsource] == "sersiccomplex") {
        sources.spatialtype[nsource] = SERSICCOMPLEX;
        nspatialpar = 14;
    }
    if (sources.spatialname[nsource] == "sersicDiskComplex") {
        sources.spatialtype[nsource] = SERSICCOMPLEX;
        nspatialpar = 14;
    }
    if (sources.spatialname[nsource] == "sersicdiskcomplex") {
        sources.spatialtype[nsource] = SERSICCOMPLEX;
        nspatialpar = 14;
    }
    if (sources.spatialname[nsource] == "sersicDiscComplex") {
        sources.spatialtype[nsource] = SERSICCOMPLEX;
        nspatialpar = 14;
    }
    if (sources.spatialname[nsource] == "sersicdisccomplex") {
        sources.spatialtype[nsource] = SERSICCOMPLEX;
        nspatialpar = 14;
    }
    if (sources.spatialname[nsource] == "sersicDisk") {
        sources.spatialtype[nsource] = SERSICDISK;
        nspatialpar = 6;
    }
    if (sources.spatialname[nsource] == "sersicdisc") {
        sources.spatialtype[nsource] = SERSICDISK;
        nspatialpar = 6;
    }
    if (sources.spatialname[nsource] == "sersicDisc") {
        sources.spatialtype[nsource] = SERSICDISK;
        nspatialpar = 6;
    }
    if (sources.spatialname[nsource] == "sersicdisk") {
        sources.spatialtype[nsource] = SERSICDISK;
        nspatialpar = 6;
    }
    if (sources.spatialname[nsource] == "distortedSphere") {
        sources.spatialtype[nsource] = DISTORTEDSPHERE;
        nspatialpar = 10;
    }
    if (sources.spatialname[nsource] == "distortedSphere") {
        sources.spatialtype[nsource] = DISTORTEDSPHERE;
        nspatialpar = 10;
    }
    if (sources.spatialname[nsource] == "movingpoint") {
        sources.spatialtype[nsource] = MOVINGPOINT;
        nspatialpar = 2;
    }
    if (sources.spatialname[nsource] == "pinhole") {
        sources.spatialtype[nsource] = PINHOLE;
        nspatialpar = 4;
    }

    if (sources.spatialtype[nsource] == 1) {
        if (nsource > 0) {
            sources.skysameas[nsource] = -1;
            for (i = 0; i < nsource; i++) {
                if (sources.spatialname[nsource] == sources.spatialname[i]) sources.skysameas[nsource] = i;
            }
        } else {
            sources.skysameas[nsource] = -1;
        }


        // read image file

        if (sources.skysameas[nsource] == -1) {
            sprintf(tempstring, "%s/%s+0", imagedir.c_str(), sources.spatialname[nsource].c_str());

            ffptr = tempstring;
            status = 0;

            if (fits_open_file(&faptr, ffptr, READONLY, &status)) printf("Error opening %s\n", ffptr);
            fits_read_keys_lng(faptr, (char*)"NAXIS", 1, 2, naxes, &nfound, &status);
            tempptr[nimage] = static_cast<float*>(malloc(naxes[0]*naxes[1]*sizeof(float)));
            naxesb[nimage][0] = naxes[1];
            naxesb[nimage][1] = naxes[0];
            fits_read_img(faptr, TFLOAT, 1, naxes[0]*naxes[1], &nullval,
                          tempptr[nimage], &anynull, &status);
            fits_close_file(faptr, &status);


            cumulative[nimage] = 0;
            cumulativex[nimage] = static_cast<float*>(malloc(naxesb[nimage][0]*sizeof(float)));
            for (i = 0; i < naxesb[nimage][0]; i++) {
                cumulativex[nimage][i] = cumulative[nimage];
                for (j = 0; j < naxesb[nimage][1]; j++) {
                    if (*(tempptr[nimage] + i*naxesb[nimage][1] + j) < 0) {
                        *(tempptr[nimage] + i*naxesb[nimage][1] + j) = 0;
                    }
                    cumulative[nimage] += *(tempptr[nimage] + i*naxesb[nimage][1] + j);
                }
            }

            sources.spatialpar[nsource][2] = nimage;
            nimage++;
        }

    }



    if (nspatialpar > 0) {
        for (i = 0; i < nspatialpar; i++) {
            iss >> sources.spatialpar[nsource][i];
        }
    }

    iss >> sources.dustnamez[nsource];

    sources.dusttypez[nsource] = 0;
    ndustparz = 0;
    if (sources.dustnamez[nsource] == "ccm") {
        sources.dusttypez[nsource] = 1;
        ndustparz = 2;
    } else if (sources.dustnamez[nsource] == "calzetti") {
        sources.dusttypez[nsource] = 2;
        ndustparz = 2;
    } else if (sources.dustnamez[nsource] == "CCM") {
        sources.dusttypez[nsource] = 1;
        ndustparz = 2;
    } else if (sources.dustnamez[nsource] == "CALZETTI") {
        sources.dusttypez[nsource] = 2;
        ndustparz = 2;
    }

    if (ndustparz > 0) {
        for (i = 0; i < ndustparz; i++) {
            iss >> sources.dustparz[nsource][i];
        }
    }


    iss >> sources.dustname[nsource];

    sources.dusttype[nsource] = 0;
    ndustpar = 0;
    if (sources.dustname[nsource] == "ccm") {
        sources.dusttype[nsource] = 1;
        ndustpar = 2;
    } else if (sources.dustname[nsource] == "calzetti") {
        sources.dusttype[nsource] = 2;
        ndustpar = 2;
    } else if (sources.dustname[nsource] == "CCM") {
        sources.dusttype[nsource] = 1;
        ndustpar = 2;
    } else if (sources.dustname[nsource] == "CALZETTI") {
        sources.dusttype[nsource] = 2;
        ndustpar = 2;
    }

    if (ndustpar > 0) {
        for (i = 0; i < ndustpar; i++) {
            iss >> sources.dustpar[nsource][i];
        }
    }

    double spra, spdec;
    spra=pra;
    spdec=pdec;
    if (aberration == 1) {
        double xprime, yprime, zprime;
        aberrationAxes(tai, &xprime, &yprime, &zprime);
        aberrationShift(&spra, &spdec, xprime, yprime, zprime);
    }
    if (precession == 1) {
        precessionShift(&spra, &spdec, tai, 51544.5);
    }
    if (nutation ==1) {
        double dlong, dobl, obl;
        nutationValues(tai,&dlong,&dobl);
        obl=obliquity(tai);
        nutationShift(&spra,&spdec,obl,dobl,dlong);
    }
    setup_tangent(spra, spdec, &tpx, &tpy, &tpz);

    ra = sources.ra[nsource] + sources.deltara[nsource];
    dec = sources.dec[nsource] + sources.deltadec[nsource];
    if (aberration == 1) {
        double xprime, yprime, zprime;
        aberrationAxes(tai, &xprime, &yprime, &zprime);
        aberrationShift(&ra, &dec, xprime, yprime, zprime);
    }
    if (precession == 1) {
        precessionShift(&ra, &dec, tai, 51544.5);
    }
    if (nutation ==1) {
        double dlong, dobl, obl;
        nutationValues(tai,&dlong,&dobl);
        obl=obliquity(tai);
        nutationShift(&ra,&dec,obl,dobl,dlong);
    }
    tangent(ra, dec, &x, &y, &tpx, &tpy, &tpz);

    sources.vx[nsource] = x*cos(rotatez) - y*sin(rotatez) + rotatex;
    sources.vy[nsource] = x*sin(rotatez) + y*cos(rotatez) + rotatey;
    sources.vz[nsource] = -1.0;
    nn = sqrt((sources.vx[nsource])*(sources.vx[nsource]) +
              (sources.vy[nsource])*(sources.vy[nsource]) + 1);
    sources.vx[nsource] = sources.vx[nsource]/nn;
    sources.vy[nsource] = sources.vy[nsource]/nn;
    sources.vz[nsource] = sources.vz[nsource]/nn;

    nsource++;
    return(0);
}

int Observation::addOpd (const std::string & opd) {

    double x, y;
    double nn;
    double ra, dec;
    double wave;
    std::string id;

    if (nsource >= MAX_SOURCE) {
        std::perror("Too many OPD sources");
    }

    std::istringstream iss(opd);
    iss >> id >> ra >> dec;
    sources.id.push_back(id);
    sources.split.push_back(0);
    sources.ra.push_back(ra);
    sources.dec.push_back(dec);
    sources.ra[nsource] *= DEGREE;
    sources.dec[nsource] *= DEGREE;
    iss >> wave;
    sources.redshift.push_back(0.0);
    sources.gamma1.push_back(wave);
    sources.gamma2.push_back(0.0);
    sources.kappa.push_back(0.0);
    sources.deltara.push_back(0.0);
    sources.deltadec.push_back(0.0);
    sources.type[nsource] = 6;
    sources.norm[nsource] = opdsize*opdsize*opdsampling*opdsampling;
    sources.mag[nsource] = opdsize*opdsize*opdsampling*opdsampling;
    sources.spatialtype[nsource] = OPD;
    sources.dusttypez[nsource] = 0;
    sources.dusttype[nsource] = 0;

    setup_tangent(0.0, 0.0, &tpx, &tpy, &tpz);

    tangent(sources.ra[nsource] + sources.deltara[nsource], sources.dec[nsource] + sources.deltadec[nsource], &x, &y, &tpx, &tpy, &tpz);

    sources.vx[nsource] = x;
    sources.vy[nsource] = y;
    sources.vz[nsource] = -1.0;
    nn = sqrt((sources.vx[nsource])*(sources.vx[nsource]) +
              (sources.vy[nsource])*(sources.vy[nsource]) + 1);
    sources.vx[nsource] = sources.vx[nsource]/nn;
    sources.vy[nsource] = sources.vy[nsource]/nn;
    sources.vz[nsource] = sources.vz[nsource]/nn;

    nsource++;
    return(0);
}


int Observation::background () {

    // Bookkeeping setup
    char tempstring[4096];
    double focalLength = platescale/DEGREE;
    double aa, f, xp, yp, ra, dec, xv, yv, maxDistanceY, currbuffer, sourceDistanceX, sourceDistanceY, maxDistanceX;
    double dx , dy;
    long ndeci,  nrai, deci, rai;
    long over, i, j;
    double dra, dis, cosdis;
    int nspatialpar, ii;
    char line[4096];
    std::string dir;
    std::string dir2;
    double x, y;
    double nn;
    double mag = 100;
    int diffusetype;
    double angularSeparationDeg, angularSepRadians, moonMagnitude, moonApparentMagnitude,
        scatter_function, moon_illuminance, lunar_illuminance, darkskyc_magnitude, darkskyp_magnitude;

    if (flatdir == 0) {
        dir = datadir + "/sky";
        dir2 = "../sky";
    } else if (flatdir == 1) {
        dir = ".";
        dir2 = ".";
    }

    if (atmospheremode) {
        airglow = static_cast<float*>(calloc((airglowScreenSize)*(airglowScreenSize), sizeof(float)));
        sprintf(tempstring, "airglowscreen_%s.fits.gz", obshistid.c_str());
        fitsReadImage(tempstring, airglow);
    }

    // Calculations setup
    moonApparentMagnitude = -12.73 + 0.026 * (fabs(phaseang/DEGREE)) + (4E-9) * pow(phaseang/DEGREE,  4);
    moon_illuminance = pow(10, -0.4 * (moonApparentMagnitude + 16.57));
    lunar_illuminance = 26.3311157 - 1.08572918 * log(moon_illuminance);

    // Twighlight
    float temp;
    float backgroundMagnitude = -2.5*log10(0.5*pow(10.,-0.4*airglowpintensity) + 0.5*pow(10.,-0.4*airglowcintensity));
    float background_brightness = 34.08 * exp(20.72333 - 0.92104 * backgroundMagnitude);
    float darksky[6];
    float darksky_data[2511][2];
    char darksky_sedfile[4096];
    FILE *darkskyF;

    sprintf(darksky_sedfile, "%s/darksky_sed.txt", dir.c_str());

    darkskyF = fopen(darksky_sedfile, "r" );
    for (i = 0; i < 2511; i++ ){
        fgets(line, 4096, darkskyF );
        sscanf(line, "%f %f\n", &darksky_data[i][0], &temp );
        darksky_data[i][1] = temp * background_brightness / 6.299537E-18;
    }
    fclose(darkskyF);
    darksky[0] = darksky_data[39][1];
    darksky[1] = darksky_data[491][1];
    darksky[2] = darksky_data[1019][1];
    darksky[3] = darksky_data[1547][1];
    darksky[4] = darksky_data[1961][1];
    darksky[5] = darksky_data[2250][1];

    float lunar[6];
    float lunarData[7500][2];
    char lunarSedfile[4096];
    FILE *lunarFp;

    sprintf(lunarSedfile, "%s/lunar_sed.txt", dir.c_str());

    lunarFp = fopen (lunarSedfile, "r" );
    for (i = 0; i < 1441; i++){
        fgets( line, 200, lunarFp );
        sscanf( line, "%f %f\n", &lunarData[i][0], &temp );
        lunarData[i][1] = temp * background_brightness / 3.882815E-16;
    }
    fclose(lunarFp);

    lunar[0] = lunarData[600][1];
    lunar[1] = lunarData[1800][1];
    lunar[2] = lunarData[3200][1];
    lunar[3] = lunarData[4600][1];
    lunar[4] = lunarData[5700][1];
    lunar[5] = lunarData[6700][1];

    float a[6][3] = {{
            11.78, 1.376, -0.039
        }, {
            11.84, 1.411, -0.041
        }, {
            11.84, 1.518, -0.057
        }, {
            11.40, 1.567, -0.064
        }, {
            10.93, 1.470, -0.062
        }, {
            10.43, 1.420, -0.052
        }
    };
    float color[6] = { 0.67, 1.03, 0, -0.74, -1.90, -2.20 };

    float magnitude = 100.0, brightness = 0.0, angle = solarzen / PI * 180.0;
    j = filter;
    float alpha = 0.0, beta = 0.0;

    if ( angle <= 106 ) {
        alpha = 1.0 - ( angle - 95.0 ) / 11.0;
        beta = 1.0 - alpha;
        magnitude = a[j][0] + a[j][1] * ( angle - 95.0 ) + a[j][2] * ( angle - 95.0 ) * ( angle - 95.0 );
    } else if ( angle >= 254 ){
        alpha = 1.0 - ( 265 - angle ) / 11.0;
        beta = 1.0 - alpha;
        magnitude = a[j][0] + a[j][1] * ( 265.0 - angle ) + a[j][2] * ( 265.0 - angle ) * ( 265.0 - angle );
    } else if ((angle > 106) && (angle < 254)) {
        alpha = 0.0;
        beta = 1.0;
        magnitude = a[j][0] + a[j][1] * ( 11.0 ) + a[j][2] * ( 11.0 ) * ( 11.0 );
    }

    brightness = 34.08 * exp( 20.72333 - 0.92104 * ( magnitude - color[j] ) );

    float lunarMagnitude, darkskyMagnitude;
    float lunarBrightness, darkskyBrightness;
    float lunarPhotonCount, darkskyPhotonCount;

    lunarBrightness = 0.5 * ( brightness - 2.0 * alpha * lunar[j] );
    darkskyBrightness = 0.5 * ( brightness - 2.0 * beta * darksky[j] );

    if (lunarBrightness < 0) lunarBrightness = 0.0;
    if (darkskyBrightness < 0) darkskyBrightness = 0.0;

    lunarMagnitude = 26.33111 - 1.08573 * log(lunarBrightness) + color[j];
    darkskyMagnitude = 26.33111 - 1.08573 * log(darkskyBrightness) + color[j];

    lunarPhotonCount = expf(-0.4 * lunarMagnitude * 2.30258509);
    if (lunarPhotonCount < 0) lunarPhotonCount = 0;

    darkskyPhotonCount = expf(-0.4 * darkskyMagnitude * 2.30258509);
    if (darkskyPhotonCount < 0) darkskyPhotonCount = 0;

    if ((angle > 106) && (angle < 130)) darkskyPhotonCount *= exp( 1 - 24.0 / fabs( angle - 130.0 ) );
    if ((angle > 230) && (angle < 254)) darkskyPhotonCount *= exp( 1 - 24.0 / fabs( angle - 230.0 ) );
    if ((angle >= 130) && (angle <= 230)) darkskyPhotonCount = 0;

    if ((angle > 106) && (angle < 130)) lunarPhotonCount *= exp( 1 - 24.0 / fabs( angle - 130.0 ) );
    if ((angle > 230) && (angle < 254)) lunarPhotonCount *= exp( 1 - 24.0 / fabs( angle - 230.0 ) );
    if ((angle >= 130) && (angle <= 230)) lunarPhotonCount = 0;

    // Zodiacal light
    double zl; // Zodiacal light luminosity
    double xx;
    double yy;
    double zz;
    double xEcliptic;
    double yEcliptic;
    double zEcliptic;
    double lambdaEcliptic; // Ecliptic longitude
    double betaEcliptic; // Ecliptic latitude
    double eclipticObliquity = 23.4*DEGREE;
    double lambdaSun; // Sun longitude
    double deltaLambda; // Ecliptic longitude - Sun longitude 

    // Source type loop
    for (diffusetype = 0; diffusetype < 6; diffusetype++){
        // Conditional checks for observational parameters and modes
        if ((diffusetype == 0 && domelight < 100 && spaceMode == 0) ||
            (diffusetype == 1 && airglowcintensity < 100 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
            (diffusetype == 2 && airglowpintensity < 100 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
            (diffusetype == 3 && moonalt > 0 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
            (diffusetype == 4 && moonalt > 0 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
            (diffusetype == 5 && backgroundMode == 1 && telconfig == 0)) {

            // Populate sources
            ndeci = (long)(180*3600/backRadius);
            over = (long)(60/backRadius);
            if (over == 0) over=1;
            for (deci = 0; deci < ndeci; deci += over) {
                dec = ((deci + 0.5 - ndeci/2)/(ndeci))*PI;
                nrai = ((long)(ndeci*cos(dec)*2));

                // Coarse placement
                for (rai = 0; rai < nrai; rai += over) {
                    ra = 2*PI*((rai + 0.5 - nrai/2)/(nrai));

                    aa = cos(dec)*cos(ra - pra);
                    f = focalLength/(sin(pdec)*sin(dec) + aa*cos(pdec));
                    yp = f*(cos(pdec)*sin(dec) - aa*sin(pdec));
                    xp = f*cos(dec)*sin(ra - pra);
                    xv = (xp*cos(-rotatez) + yp*sin(-rotatez));
                    yv = (-xp*sin(-rotatez) + yp*cos(-rotatez));
                    currbuffer = backBuffer + (3*backBeta*backRadius + 60)*ARCSEC*focalLength/pixsize;
                    maxDistanceX = pixelsx * pixsize/2.0 + currbuffer * pixsize;
                    maxDistanceY = pixelsy * pixsize/2.0 + currbuffer * pixsize;
                    dx = xv - centerx - decenterx;
                    dy = yv - centery - decentery;
                    sourceDistanceX = fabs(cos(chipangle)*dx + sin(chipangle)*dy);
                    sourceDistanceY = fabs(-sin(chipangle)*dx + cos(chipangle)*dy);

                    dra = fabs(ra - pra);
                    if (dra > PI) dra = 2*PI - dra;
                    cosdis = sin(dec)*sin(pdec) + cos(dec)*cos(pdec)*cos(dra);
                    if (cosdis > 1) cosdis = 1.0;
                    if (cosdis < -1) cosdis = -1.0;
                    dis = acos(cosdis);

                    if ((sourceDistanceX <= maxDistanceX) && (sourceDistanceY <= maxDistanceY) &&
                        (dis < PI/2) ){

                        // Fine placement
                        for (i = 0; i < over; i++) {
                            for (j = 0; j < over; j++) {
                                dec = ((deci + i + 0.5 - ndeci/2)/ndeci)*PI;
                                ra = 2*PI*((rai + j + 0.5 - nrai/2)/(nrai));
                                aa = cos(dec)*cos(ra - pra);
                                f = focalLength/(sin(pdec)*sin(dec) + aa*cos(pdec));
                                yp = f*(cos(pdec)*sin(dec) - aa*sin(pdec));
                                xp = f*cos(dec)*sin(ra - pra);
                                xv = (xp * cos(-rotatez) + yp * sin(-rotatez));
                                yv = (-xp * sin(-rotatez) + yp * cos(-rotatez));
                                currbuffer = backBuffer + 3*backBeta*backRadius*ARCSEC*focalLength/pixsize;
                                maxDistanceX = pixelsx * pixsize/2.0 + currbuffer * pixsize;
                                maxDistanceY = pixelsy * pixsize/2.0 + currbuffer * pixsize;
                                dx = xv - centerx - decenterx;
                                dy = yv - centery - decentery;
                                sourceDistanceX = fabs(cos(chipangle)*dx + sin(chipangle)*dy);
                                sourceDistanceY = fabs(-sin(chipangle)*dx + cos(chipangle)*dy);

                                dra = fabs(ra - pra);
                                if (dra > PI) dra = 2*PI - dra;
                                cosdis = sin(dec)*sin(pdec) + cos(dec)*cos(pdec)*cos(dra);
                                if (cosdis > 1) cosdis = 1.0;
                                if (cosdis < -1) cosdis = -1.0;
                                dis = acos(cosdis);

                                if ((sourceDistanceX <= maxDistanceX) && (sourceDistanceY <= maxDistanceY) &&
                                    (dis < PI/2)){

                                    dra = fabs(ra - moonra);
                                    if (dra > PI) dra = 2*PI - dra;
                                    cosdis = sin(dec)*sin(moondec) + cos(dec)*cos(moondec)*cos(dra);
                                    if (cosdis > 1) cosdis = 1.0;
                                    if (cosdis < -1) cosdis = -1.0;
                                    angularSepRadians = acos(cosdis);
                                    angularSeparationDeg  =  angularSepRadians * 180.0 / PI;

                                    if (diffusetype == 3) {
                                        scatter_function = 16.57 + 26.3311157 - 2.5*log10(4.0*PI) +
                                            2.5*log10(PI/180.0/3600.0*PI/180.0/3600.0) +
                                            2.5*log10(2.0/2.78666 * (1.06 + cos(angularSepRadians) * cos(angularSepRadians)));
                                    } else {
                                        scatter_function = 16.57 + 26.3311157 - 2.5*log10(4.0*PI) +
                                            2.5*log10(PI/180.0/3600.0*PI/180.0/3600.0) +
                                            2.5*log10(2.0/0.284181 * (exp((- angularSeparationDeg / 40.0)*2.30258509) + 1e-1));
                                    }
                                    if (moonalt < 0) moonMagnitude = 10000;
                                    else moonMagnitude = lunar_illuminance - scatter_function;

                                    moonMagnitude = -2.5*log10(lunarPhotonCount + pow(10.0, -0.4*moonMagnitude));
                                    darkskyc_magnitude = -2.5*log10(darkskyPhotonCount + pow(10.0, -0.4*airglowcintensity));
                                    darkskyp_magnitude = -2.5*log10(darkskyPhotonCount + pow(10.0, -0.4*airglowpintensity));
 
                                    // Airglow variation
                                    long ax0, ax1, ay0, ay1;
                                    double dx, dy;
                                    find_linear_wrap(xv, platescale*15.0/3600, airglowScreenSize, &ax0, &ax1, &dx);
                                    find_linear_wrap(yv, platescale*15.0/3600, airglowScreenSize, &ay0, &ay1, &dy);

                                    if (diffusetype == 0) mag = domelight - 2.5*log10(backRadius*backRadius);
                                    if (diffusetype == 1) {
                                        double airglowv;
                                        airglowv = airglowvariation*(static_cast<double>(interpolate_bilinear_float_wrap(airglow,
                                            airglowScreenSize, ax0, ax1, dx, ay0, ay1, dy))); 
                                        mag = darkskyc_magnitude + airglowv - 2.5*log10(backRadius*backRadius);
                                    }
                                    if (diffusetype == 2) mag = darkskyp_magnitude - 2.5*log10(backRadius*backRadius);
                                    if (diffusetype == 3) mag = moonMagnitude - 2.5*log10(backRadius*backRadius);
                                    if (diffusetype == 4) mag = moonMagnitude - 2.5*log10(backRadius*backRadius);
                                    if (diffusetype == 5) {

                                        if (overrideZodiacalLightMagnitude == 0) {
         
                                            // Calculate ecliptic longitude and latitude
                                            // Source position in Cartesian coordinates 
                                            xx = cos(ra) * cos(dec);
                                            yy = sin(ra) * cos(dec);
                                            zz = sin(dec);
                                            // Transform to ecliptic coordinates
                                            xEcliptic = xx;
                                            yEcliptic = yy*cos(eclipticObliquity) + zz*sin(eclipticObliquity);
                                            zEcliptic = -yy*sin(eclipticObliquity) + zz*cos(eclipticObliquity);
                                            // Ecliptic longitue and latitude
                                            lambdaEcliptic = atan2(yEcliptic, xEcliptic)/DEGREE; // Degrees
                                            betaEcliptic = atan(zEcliptic/sqrt(pow(xEcliptic, 2) + pow(yEcliptic, 2)))/DEGREE; // Degrees

                                            // Zodiacal light variation
                                            // Follows Kwon et al. New Astronomy 10-2 (2004) pp91-107
                                            if (lambdaEcliptic < 0) {
                                                lambdaEcliptic += 360.0;
                                            }
                                            lambdaSun = fmod((640.466 - 0.985607*(day - 51544.0)), 360.0); // Degrees
                                            if (lambdaSun < 0) {
                                                lambdaSun += 360.0;
                                            }
                                            deltaLambda = fmod(lambdaEcliptic - lambdaSun, 360.0); // Degrees
                                            if (deltaLambda < 0) {
                                                deltaLambda += 360.0;
                                            }
                                            zl = (140.0 + 20.0/(1 + pow(fabs(180.0 - deltaLambda)/20.0, 2.0)))*
                                                (5.0/14.0 + 9.0/14.0/(1 + pow((betaEcliptic)/20.0, 2.0)));

                                            zodiacalLightMagnitude = -2.5*(log10(zl)) + 27.78;

                                        }

                                        mag = zodiacalLightMagnitude - 2.5*log10(backRadius*backRadius);

                                    }

                                    if (mag < 100) {

                                        nspatialpar = 0;

                                        if (nsource >= MAX_SOURCE) {
                                            std::perror("Too many background sources");
                                        }

                                        sources.id.push_back("0.0");
                                        sources.split.push_back(0);
                                        sources.ra.push_back(ra);
                                        sources.dec.push_back(dec);
                                        if (diffusetype == 0 && domewave == 0.0) sources.sedfilename[nsource] = dir + "/sed_dome.txt";
                                        if (diffusetype == 0 && domewave != 0.0) sources.sedfilename[nsource] = "laser";
                                        if (diffusetype == 1) sources.sedfilename[nsource] = dir + "/airglowc_sed.txt";
                                        if (diffusetype == 2) sources.sedfilename[nsource] = dir + "/airglowp_sed.txt";
                                        if (diffusetype == 3) sources.sedfilename[nsource] = dir + "/lunar_sed.txt";
                                        if (diffusetype == 4) sources.sedfilename[nsource] = dir + "/lunar_sed.txt";
                                        if (diffusetype == 5) sources.sedfilename[nsource] = dir + "/zodiacal_sed.txt";
                                        sources.redshift.push_back(0.0);
                                        sources.gamma1.push_back(0.0);
                                        sources.gamma2.push_back(0.0);
                                        sources.kappa.push_back(0.0);
                                        sources.deltara.push_back(0.0);
                                        sources.deltadec.push_back(0.0);
                                        sources.type[nsource] = diffusetype;
                                        sources.norm[nsource] = pow(10.0, ((mag + 48.6)/(-2.5)));
                                        sources.mag[nsource] = mag;

                                        // Read SED file
                                        if (nsource > 0) {
                                            for (ii = 0; ii < nsource; ii++) {
                                                if (sources.sedfilename[ii] == sources.sedfilename[nsource]) {
                                                    sources.sedptr[nsource] = sources.sedptr[ii];
                                                    goto skipsedread;
                                                }
                                            }
                                        }

                                        sources.sedptr[nsource] = nsedptr;
                                        readSed(sources.sedfilename[nsource], 1);

                                    skipsedread:;
                                        sources.norm[nsource] = sources.norm[nsource]/(normwave)*(1 + sources.redshift[nsource])*
                                            sedDwdp[sources.sedptr[nsource]];

                                        sources.spatialname[nsource] = "gauss";
                                        if (sources.spatialname[nsource] == "gauss") {
                                            sources.spatialtype[nsource] = 2;
                                            nspatialpar = 1;
                                        }
                                        sources.spatialpar[nsource][0] = backRadius*backEpsilon;
                                        sources.dustnamez[nsource] = "none";
                                        sources.dustname[nsource] = "none";

                                        setup_tangent(pra, pdec, &tpx, &tpy, &tpz);

                                        tangent(sources.ra[nsource] + sources.deltara[nsource], sources.dec[nsource] + sources.deltadec[nsource],
                                                &x, &y, &tpx, &tpy, &tpz);

                                        sources.vx[nsource] = x*cos(rotatez) - y*sin(rotatez) + rotatex;
                                        sources.vy[nsource] = x*sin(rotatez) + y*cos(rotatez) + rotatey;
                                        sources.vz[nsource] = -1.0;
                                        nn = sqrt((sources.vx[nsource])*(sources.vx[nsource]) +
                                                  (sources.vy[nsource])*(sources.vy[nsource]) + 1);
                                        sources.vx[nsource] = sources.vx[nsource]/nn;
                                        sources.vy[nsource] = sources.vy[nsource]/nn;
                                        sources.vz[nsource] = sources.vz[nsource]/nn;
                                        nsource++;

                                    }

                                }


                            }
                        }
                    }

                }
            }
        }
    }

    return(0);
}


int Observation::filterTruncateSources () {

    double filterlow,  filterhigh;
    long lowptr, highptr;
    double lowvalue, highvalue;

    filterlow = 0;
    filterhigh = maxwavelength;

    for (long i = 0; i < nsedptr; i++) {
        highvalue = 0;
        lowvalue = 0;
        lowptr = 0;
        highptr = 0;
        for (long j = 0; j < sedN[i]; j++) {

            if (sedW[sedPtr[i] + j] < filterlow) {
                lowptr = j;
                lowvalue = sedC[sedPtr[i] + j];
            }
            if (sedW[sedPtr[i] + j] < filterhigh) {
                highptr = j;
                highvalue = sedC[sedPtr[i] + j];
            }

        }

        sedCorr[i]=(sedC[sedPtr[i] + highptr] - sedC[sedPtr[i] + lowptr]);

        for (long j = lowptr; j < highptr + 1; j++) {
            sedC[sedPtr[i] + j]=(sedC[sedPtr[i] + j] - lowvalue)/
                (highvalue - lowvalue);
        }
        sedPtr[i] = sedPtr[i] + lowptr;
        sedN[i] = highptr - lowptr;
    }

    for (long i = 0; i < nsource; i++) {
        if (sources.spatialtype[i] != OPD) {
            sources.norm[i] = sources.norm[i]*sedCorr[sources.sedptr[i]];
        }
    }

    return(0);

}


int Observation::splitSources () {

    long count = 0;

    if (numthread > 1) {

        while (1) {
            long brightest = -1;
            double bright = 0.0;
            for (long i = 0; i < nsource; i++) {
                if ((sources.norm[i] > bright) && (sources.type[i] == 6) && (sources.spatialtype[i] != OPD)  && (sources.split[i]==0)) {
                    bright = sources.norm[i];
                    brightest = i;
                }
            }

            if (brightest == -1) break;
            long totalbrightsplit = 0;
            for (long i = 0; i < nsource; i++) {
                if ((sources.norm[i] > bright) && (sources.split[i]==1)) totalbrightsplit++;
            }
            if (totalbrightsplit >= numthread*2) break;

            sources.split[brightest]=1;
            sources.norm[brightest]=sources.norm[brightest]/(static_cast<double>(numthread));
            sources.mag[brightest]=sources.mag[brightest]+2.5*log10(numthread);

            for (long j = 0; j < (numthread-1); j++) sources.ra.push_back(sources.ra[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.dec.push_back(sources.dec[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.redshift.push_back(sources.redshift[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.gamma1.push_back(sources.gamma1[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.gamma2.push_back(sources.gamma2[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.kappa.push_back(sources.kappa[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.deltara.push_back(sources.deltara[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.deltadec.push_back(sources.deltadec[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.split.push_back(sources.split[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.sedfilename.push_back(sources.sedfilename[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.spatialname.push_back(sources.spatialname[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.dustname.push_back(sources.dustname[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.dustnamez.push_back(sources.dustnamez[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.id.push_back(sources.id[brightest]);
            for (long j = 0; j < (numthread-1); j++) sources.vx[nsource+j]=sources.vx[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.vy[nsource+j]=sources.vy[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.vz[nsource+j]=sources.vz[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.mag[nsource+j]=sources.mag[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.norm[nsource+j]=sources.norm[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.spatialtype[nsource+j]=sources.spatialtype[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.dusttype[nsource+j]=sources.dusttype[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.dusttypez[nsource+j]=sources.dusttypez[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.type[nsource+j]=sources.type[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.skysameas[nsource+j]=sources.skysameas[brightest];
            for (long j = 0; j < (numthread-1); j++) sources.sedptr[nsource+j]=sources.sedptr[brightest];
            for (long j = 0; j < (numthread-1); j++) for (long k = 0; k< 14; k++) sources.spatialpar[nsource+j][k]=sources.spatialpar[brightest][k];
            for (long j = 0; j < (numthread-1); j++) for (long k = 0; k< 2; k++) sources.dustpar[nsource+j][k]=sources.dustpar[brightest][k];
            for (long j = 0; j < (numthread-1); j++) for (long k = 0; k< 2; k++) sources.dustparz[nsource+j][k]=sources.dustparz[brightest][k];

            for (long j = 0; j < (numthread-1); j++) nsource++;

            count++;
            if (count > 1000) {
                std::cout << "Error with too many split sources.  This should not occur." << std::endl;
                exit(1);
            }
        }

    }
    return(0);

}
