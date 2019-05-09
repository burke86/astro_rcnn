///
/// @package phosim
/// @file atmospheresetup.cpp
/// @brief setup for atmosphere (part of image class)
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

#include <stdlib.h>
#include "ancillary/readtext.h"
#include "constants.h"

int Image::atmSetup () {

    char tempstring[4096];
    std::string dir;

    if (flatdir == 0) dir = datadir + "/atmosphere";
    else if (flatdir == 1) dir = ".";

    sources.sedfilename.clear();
    sources.spatialname.clear();
    sources.dustname.clear();
    sources.dustnamez.clear();

    ranseed = obsseed + pairid;

    posix_memalign((void **)&random, PHOSIM_RANDOM_ALIGN, numthread*sizeof(Random));
    posix_memalign((void **)&galaxy.random, PHOSIM_RANDOM_ALIGN, numthread*sizeof(Random));

    for (int i = 0; i < numthread; i++) {
        new (&random[i]) Random();
        new (&galaxy.random[i]) Random();
    }
    if (obsseed == -1) {
        for (int i = 0; i < numthread; i++) {
            random[i].setSeedFromTime();
            galaxy.random[i].setSeedFromTime();
        }
    } else {
        for (int i = 0; i < numthread; i++) {
            random[i].setSeed32(ranseed+i*100);
            galaxy.random[i].setSeed32(ranseed+9999+i*100);
        }
    }

    for (int i = 0; i < numthread; i++) {
        random[i].unwind(100);
        galaxy.random[i].unwind(100);
    }

    if (devtype == "CMOS") {
        exptime = vistime/nsnap;
        timeoffset = 0.5*exptime + pairid*(exptime) - 0.5*vistime;
    } else {
        exptime = (vistime - (nsnap - 1)*devvalue)/nsnap;
        timeoffset = 0.5*exptime + pairid*(devvalue + exptime) - 0.5*vistime;
    }

    exptime = exptime + shuttererror*random[0].normal();

    // SKY SETUP
    char sersicdata[4096];
    if (flatdir == 0) sprintf(sersicdata, "%s/sky/sersic_const.txt", datadir.c_str());
    else if (flatdir == 1) sprintf(sersicdata, "sersic_const.txt");
    galaxy.sampleSersic(sersicdata);
    galaxy.sampleSersic2d(sersicdata, numthread);
    dust.setup(maxwavelength);
    filterTruncateSources();
    if (splitsources == 1) splitSources();
    totalnorm = 0.0;
    for (long i = 0; i < nsource; i++) {
        totalnorm += sources.norm[i];
    }

    // ATMOSPHERIC SETUP
    for (long i = 0; i < natmospherefile; i++) {
        seefactor[i] = pow(1/cos(zenith), 0.6)*seefactor[i];
    }

    screen.large_sizeperpixel = 5120.0;
    screen.coarse_sizeperpixel = 640.0;
    screen.medium_sizeperpixel = 80.0;
    screen.fine_sizeperpixel = 10.0;

    screen.wavelengthfactor_nom = pow(0.5, -0.2);

    // Air Setup
    if (atmospheremode) {
        fprintf(stdout, "Creating Air.\n");
        air.opacitySetup(zenith, moonalt, height, groundlevel, raynorm,
                     o2norm, h2onorm, o3norm, aerosoltau, aerosolindex,
                     natmospherefile, dir, &airmass);
    }

    // Diffraction Screens
    screen.phasescreen = new double[SCREEN_SIZE*SCREEN_SIZE]();
    screen.focalscreen = new double[SCREEN_SIZE*SCREEN_SIZE]();
    screen.tfocalscreen = new double[SCREEN_SIZE*SCREEN_SIZE]();
    screen.outscreen = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*SCREEN_SIZE*SCREEN_SIZE);
    screen.inscreen = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*SCREEN_SIZE*SCREEN_SIZE);
    screen.pupil_values = new float[SCREEN_SIZE*SCREEN_SIZE]();
    screen.focalscreencum = new double[SCREEN_SIZE*SCREEN_SIZE]();

    // Atmospheric Turbulence & Clouds
    if (atmospheremode) {
        fprintf(stdout, "Generating Turbulence.\n");
        screen.turbulenceLargeX = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceLargeY = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceCoarseX = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceCoarseY = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceMediumX = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceMediumY = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();

        screen.phaseLarge = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.phaseCoarse = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.phaseMedium = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.phaseFine = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();

        screen.phaseMediumH = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.phaseFineH = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();

        screen.see_norm = new float[natmospherefile]();
        screen.phase_norm = new float[natmospherefile]();

        for (long layer = 0; layer < natmospherefile; layer++) {

            if (cloudvary[layer] != 0 || cloudmean[layer] != 0) {
                screen.cloud[layer] = static_cast<float*>(calloc(SCREEN_SIZE*SCREEN_SIZE, sizeof(float)));
                sprintf(tempstring, "%s.fits.gz", cloudfile[layer].c_str());
                screen.readScreen(-1, screen.cloud[layer], tempstring);
             }

            sprintf(tempstring, "%s_largep.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseLarge + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            sprintf(tempstring, "%s_coarsep.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseCoarse + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            sprintf(tempstring, "%s_mediump.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseMedium + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            sprintf(tempstring, "%s_finep.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseFine + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);

            sprintf(tempstring, "%s_mediumh.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseMediumH + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            sprintf(tempstring, "%s_fineh.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseFineH + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);

            sprintf(tempstring, "%s_largex.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceLargeX + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            sprintf(tempstring, "%s_largey.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceLargeY + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            sprintf(tempstring, "%s_coarsex.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceCoarseX + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            sprintf(tempstring, "%s_coarsey.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceCoarseY + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            sprintf(tempstring, "%s_mediumx.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceMediumX + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            sprintf(tempstring, "%s_mediumy.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceMediumY + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);

            for (long i = 0; i < SCREEN_SIZE; i++) {
                for (long j = 0; j < SCREEN_SIZE; j++) {
                    *(screen.turbulenceLargeX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceLargeX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceLargeY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceLargeY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceCoarseX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceCoarseX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceCoarseY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceCoarseY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceMediumX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceMediumX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceMediumY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceMediumY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                }
            }

        }
    }

    double sigma=1.0/centralwavelength;
    double temp=temperature+273.15;
    double ps=pressure/760.00*1013.25;
    double pw=waterPressure/760.00*1013.25;
    double dw=(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/temp-710.792/temp/temp+7.75141e4/temp/temp/temp))*pw/temp;
    double ds=(1+ps*(57.90e-8-9.325e-4/temp+0.25844/temp/temp))*ps/temp;
    double n=(2371.34+683939.7/(130.0-pow(sigma,2))+4547.3/(38.9-pow(sigma,2)))*ds;
    n=n+(6478.31-58.058*pow(sigma,2)-0.71150*pow(sigma,4)+0.08851*pow(sigma,6))*dw;
    air.air_refraction_adc = 1e-8*n;

    // air.air_refraction_adc = 64.328 + 29498.1/(146 - 1/centralwavelength/centralwavelength) + 255.4/(41 - 1/centralwavelength/centralwavelength);
    // air.air_refraction_adc = air.air_refraction_adc*pressure*(1 + (1.049 - 0.0157*temperature)*1e-6*pressure)/720.883/(1 + 0.003661*temperature);
    // air.air_refraction_adc = air.air_refraction_adc - ((0.0624 - 0.000680/centralwavelength/centralwavelength)/
    //                                                    (1 + 0.003661*temperature)*waterPressure);
    // air.air_refraction_adc = air.air_refraction_adc/1e6;

    return(0);

}
