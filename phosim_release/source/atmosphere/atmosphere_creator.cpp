///
/// @package phosim
/// @file turb2d.cpp
/// @brief turbulence generator
///
/// @brief Created by:
/// @author Mallory Young (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.

float mean(float nums[], int howmany) {
    // calculates the mean of a given array of values
    float sum = 0.0;
    int i = 0;
    for (i = 0; i < howmany; i++) {
        sum = sum + nums[i];
    }
    float mean = 0.0;
    mean = sum/howmany;
    return mean;
}

float stdev(float nums[], int howmany) {
    // calculates the standard deviation of a given array of values
    int n = 0;
    float stdv = 0;
    float ave = 0.0;
    float sum = 0.0;
    float sumsqr = 0.0;
    for (n = 0; n < howmany; n++) {
        sum = sum + nums[n];
    }
    ave = sum/howmany;
    for (n = 0; n < howmany; n++) {
        sumsqr = sumsqr + pow(nums[n] - ave, 2);
    }
    stdv = sumsqr/howmany;
    stdv = sqrt(stdv);
    return stdv;
}

float interpol(float xx, float xs[], float ys[], int hm) {
    // linearly interpolates between two known data arrays
    float x1 = 0, x2 = 0, h1 = 0, h2 = 0, h = 0;
    int i = 0, j = 0;
    for (i = 0; i < hm; i++) {
        if (xs[i] >= xx) {
            x1 = xs[i];
            j = i;
        }
    }
    x2 = xs[j + 1];
    h1 = ys[j];
    h2 = ys[j + 1];
    h = h1 + ((xx - x1)*(h2 - h1))/(x2 - x1);
    return h;
}

void Atmosphere::createAtmosphere(float monthnum, float constrainseeing,
                                  const std::string & outputfilename, const std::vector<int> & cloudscreen, long seed, double tai, double *tseeing) {

    /// Generates estimates for the following atmospheric parameters: wind speed, wind
    /// direction, Cn^2, outer scale as a function of height
    /// Method: Please see the corresponding paper entitled "Atmospheric Parameters for the Image
    /// Simulator Based on the Site of the Large Synoptic Survey Telescope" by Young for further
    /// discussion.

    random.setSeed32(seed);
    random.unwind(10000);

    float totalj;
    float totalseeing;
    float bestseeing = 0.0;
    std::vector<float> relseeing(numlevel, 0);

    altitudes[numlevel - 1] = 10.0 + groundlevel;
    altitudes[numlevel - 2] = 20.0 + groundlevel;
    for (int i = 0; i < numlevel - 2; i++) {
        altitudes[i] = pow(2.0, (numlevel - 2) - i/((static_cast<float>(numlevel) - 2)/6.0))
            *(16000.0/pow(2.0, numlevel - 2)) + groundlevel;
    }

    std::vector<float> magests(numlevel, 0);
    std::vector<float> dirests(numlevel, 0);
    double r1 = random.uniformCorrel(tai, 600.0/(24.0*3600.0), 0);
    double r2 = random.normalCorrel(tai + 1000, 600.0/(24.0*3600.0));
    for (int i = 0; i < numlevel; i++) {
        magcalc(monthnum, altitudes[i], magests[i], dirests[i], r1, r2);
        if (i < (numlevel-1)) osests[i] = outerscale(altitudes[i], tai, i);
        if (i == (numlevel-1)) osests[i] = 1.0;
        if (i == 0 || i == 1) {
            osests[i] *= 1.5;
        }
        if (osests[i] > (altitudes[i] - groundlevel)/2) {
            osests[i] = (altitudes[i] - groundlevel)/2;
        }
    }

    if (constrainseeing < 0.0) {
        ccalc(tai, 0);
        totalj = 0.0;
        for (int b = 0; b < numlevel-1; b++) {
            totalj += jests[b];
        }
        totalseeing = 5.25/(pow(5000e-10, 0.2))/(2.0*sqrt(2.0*log(2.0)))
            /(M_PI/180/3600)*pow(totalj, 0.6);
        for (int b = 0; b < numlevel-1; b++) {
            relseeing[b] = pow(jests[b]/totalj, 0.6);
        }

        float renorm = 0.0;
        for (int b = 0; b < numlevel-1; b++) {
            renorm += relseeing[b]*relseeing[b];
        }
        renorm = sqrt(renorm);
        for (int b=0; b < numlevel-1; b++) {
            relseeing[b] = relseeing[b]*totalseeing/renorm;
        }
    } else if (constrainseeing > 0.0) {
        int bestindex = 0;
        float seeingdiff = 0.0, minseeingdiff = 0.0;
        std::vector<std::vector<float> > jestarray(100);
        float besttotalj = 0.0, totalj = 0.0;
        for (int c = 0; c < 100; c++) {
            // 100 different versions of cests are created
            jestarray[c].resize(numlevel-1, 0);
            ccalc(tai, c);
            for (int b = 0; b < numlevel-1; b++){
                jestarray[c][b] = jests[b];
            }
            totalj = 0.0;
            for (int b = 0; b < numlevel-1; b++) {
                totalj += jestarray[c][b];
            }
            totalseeing = 5.25/(pow(5000e-10, 0.2))/(2.0*sqrt(2.0*log(2.0)))
                /(M_PI/180/3600)*pow(totalj, 0.6);
            // here the totalseeing value closest to constrainseeing is
            // found by comparing differences between the two
            seeingdiff = fabs(constrainseeing - totalseeing);
            if (c == 0 || seeingdiff < minseeingdiff) {
                minseeingdiff = seeingdiff;
                bestseeing = totalseeing;
                bestindex = c;
                besttotalj = totalj;
            }
        }
        totalseeing = constrainseeing;
        // the totalseeing value that is closest to constrainseeing is
        // used to generate the relseeing at each level
        for (int b = 0; b < numlevel-1; b++){
            relseeing[b] = pow(jestarray[bestindex][b]/besttotalj, 0.6);
        }
        float renorm = 0.0;
        for (int b = 0; b < numlevel-1; b++) {
            renorm += relseeing[b]*relseeing[b];
        }
        renorm = sqrt(renorm);
        for (int b = 0; b < numlevel-1; b++) {
            relseeing[b] = relseeing[b]*constrainseeing/renorm;
        }
    } else {
        for (int c = 0; c < numlevel-1; c++) {
            relseeing[c] = 0.0;
            jests[c] = 0.0;
        }
        totalseeing = 0.0;
    }

    double cloudmean1, cloudsigma1, cloudmean2, cloudsigma2, cloudvariationscale, domeseeingmean=0.0, domeseeingmedian=0.0;
    readText pars(instrdir + "/site.txt");
    for (size_t t(0); t < pars.getSize(); t++) {
        readText::get(pars[t], "cloudmean1", cloudmean1);
        readText::get(pars[t], "cloudsigma1", cloudsigma1);
        readText::get(pars[t], "cloudmean2", cloudmean2);
        readText::get(pars[t], "cloudsigma2", cloudsigma2);
        readText::get(pars[t], "cloudvariationscale", cloudvariationscale);
        readText::get(pars[t], "domeseeingmean", domeseeingmean);
        readText::get(pars[t], "domeseeingmedian", domeseeingmedian);
    }

    //domeseeing
    double domeseeing=0.0;
    if ((domeseeingmean!=0.0) && (domeseeingmedian != 0.0)) {
        float mu = log(domeseeingmedian);
        float sigma = pow(2.0*(log(domeseeingmean) - mu), 0.5);
        mu = log(domeseeingmean) - (1/2)*(pow(sigma, 2.0));
        domeseeing = exp(mu + sigma*random.normalCorrel(tai + 1234, 600.0/(24.0*3600.0)));
    }
    relseeing[numlevel -1]=domeseeing/2.35;

    FILE *afile;
    afile = fopen(outputfilename.c_str(), "wt");
    float rvalue[5];
    rvalue[0] = 1.0 + 0.01*random.normal();
    rvalue[1] = 1.0 + 0.0002*random.normal();
    rvalue[2] = exp(0.18*random.normal());
    rvalue[3] = exp(0.20*random.normal());
    rvalue[4] = 0.02 + 0.01*random.normal();
    if (rvalue[4] < 0) rvalue[4] = 0.0;
    fprintf(afile, "natmospherefile %i\n", numlevel);
    fprintf(afile, "totalseeing %f\n", totalseeing*(2.0*sqrt(2.0*log(2.0))));
    *tseeing=totalseeing*(2.0*sqrt(2.0*log(2.0)));
    fprintf(afile, "reldensity %f\n", rvalue[0]);
    fprintf(afile, "relo2 %f\n", rvalue[1]);
    fprintf(afile, "relh2o %f\n", rvalue[2]);
    fprintf(afile, "relo3 %f\n\n", rvalue[3]);
    fprintf(afile, "aerosoltau %f\n\n", rvalue[4]);
    rvalue[0] = 0.01/4.0*random.uniform();
    rvalue[1] = 0.0002/4.0*random.uniform();
    rvalue[2] = 0.18/4.0*random.uniform();
    rvalue[3] = 0.20/4.0*random.uniform();
    rvalue[4] = 0.01/4.0*random.uniform();
    fprintf(afile, "raygradient %f\n\n", rvalue[0]);
    fprintf(afile, "o2gradient %f\n\n", rvalue[1]);
    fprintf(afile, "o3gradient %f\n\n", rvalue[2]);
    fprintf(afile, "h2ogradient %f\n\n", rvalue[3]);
    fprintf(afile, "aerosolgradient %f\n\n", rvalue[4]);
    fprintf(afile, "rayangle %f\n\n", random.uniform()*2.0*M_PI);
    fprintf(afile, "o2angle %f\n\n", random.uniform()*2.0*M_PI);
    fprintf(afile, "o3angle %f\n\n", random.uniform()*2.0*M_PI);
    fprintf(afile, "h2oangle %f\n\n", random.uniform()*2.0*M_PI);
    fprintf(afile, "aerosolangle %f\n\n", random.uniform()*2.0*M_PI);
    for (int i = 0; i < numlevel; i++) {
        fprintf(afile, "height %d %f\n", i, (altitudes[i] - groundlevel)/1000.0);
        fprintf(afile, "wind %d %f\n", i, magests[i]);
        fprintf(afile, "winddir %d %f\n", i, dirests[i]*(180.0/M_PI));
        fprintf(afile, "outerscale %d %f\n", i, osests[i]);
        fprintf(afile, "seeing %d %f\n", i, relseeing[i]);

        if (cloudscreen[i]) {
            float cloudmean = 0.0, cloudvary;
            int counter = 0;
        redo:;
            if (i == 1) cloudmean = exp(log(cloudmean1) + cloudsigma1*random.normalCorrel(tai, 8.0/24.0));
            if (i == 2) cloudmean = exp(log(cloudmean2) + cloudsigma2*random.normalCorrel(tai + 1000, 8.0/24.0));
            if (constrainclouds > 0.0) {
                if ((cloudmean >= constrainclouds/2.0) && (counter<10000)) {
                        counter++;
                        goto redo;
                    }
            }
            cloudvary = random.uniformCorrel(tai, 8.0/24.0, 0)*(cloudmean/cloudvariationscale + 0.04);
            fprintf(afile, "cloudmean %d %f\n\n", i, cloudmean);
            fprintf(afile, "cloudvary %d %f\n\n", i, cloudvary);
        }
    }
    fclose(afile);

}

void Atmosphere::ccalc(double tai, int call) {

    char filename[4096];
    float data[7][11];
    FILE *file;
    sprintf(filename, "%s/site.txt", instrdir.c_str());
    file = fopen(filename, "r");
    int i, j;
    for (i = 0; i < 7; i++) {
        fscanf(file, "%f %f %f %f %f %f %f %f %f %f %f\n", &data[i][0],
               &data[i][1], &data[i][2], &data[i][3], &data[i][4], &data[i][5],
               &data[i][6], &data[i][7], &data[i][8], &data[i][9], &data[i][10]);
    }
    fclose(file);
    float j1grid[7];
    float j2grid[7];
    float jtgrid[7];
    std::vector<float> hlow(numlevel-1, 0);
    std::vector<float> hhigh(numlevel-1, 0);
    float r1, r2, xl, xh, overlap;

    for (i = 0; i < 7; i++) {

        r1 = random.normalCorrel(tai + 20000 + 1000*call + 50*i, 600.0/(3600.0*24.0));
        r2 = random.normalCorrel(tai + 20000 + 1000*call + 500 + 50*i, 600.0/(3600.0*24.0));

        if (data[i][4] != 0.0) {
            j1grid[i] = exp(log(data[i][4]) + data[i][6]*r1);
        } else {
            j1grid[i] = 0.0;
        }
        if (j1grid[i] < 0) {
            j1grid[i] = 0.0;
        }
        if (data[i][8] != 0.0) {
            j2grid[i] = exp(log(data[i][8]) + data[i][10]*r2);
        } else {
            j2grid[i] = 0.0;
        }
        if (j2grid[i] < 0) {
            j2grid[i] = 0.0;
        }
        jtgrid[i] = (j1grid[i] + j2grid[i])*1e-13;
    }
    for (i = 0; i < numlevel-1; i++) {
        if (i > 0) {
            hhigh[i] = 0.5*(altitudes[i - 1] + altitudes[i]);
        }
        if (i < numlevel - 2) {
            hlow[i] = 0.5*(altitudes[i] + altitudes[i + 1]);
        }
    }
    hlow[numlevel - 2] = 0.0;
    hhigh[0] = 20000 + groundlevel;

    for (i=0; i < numlevel-1; i++) {
        jests[i] = 0.0;
        for (j = 0; j < 7; j++) {
            if (data[j][0] > hlow[i]) {
                xl = data[j][0];
            } else {
                xl = hlow[i];
            }
            if (data[j][2] < hhigh[i]) {
                xh = data[j][2];
            } else {
                xh = hhigh[i];
            }
            overlap = (xh - xl)/(data[j][2] - data[j][0]);
            if (data[j][2] < hlow[i]) {
                overlap = 0.0;
            }
            if (data[j][0] > hhigh[i]) {
                overlap = 0.0;
            }
            jests[i] = jests[i] + overlap*jtgrid[j];
        }
    }
}

void Atmosphere::magcalc(float monthnum, float altitude, float & magest, float & direst, double r1, double r2) {

    /// @fn void AtmosphereCreator::magcalc(float monthnum, float altitude, float & magest, float & direst)
    /// @brief creates wind speed and direction estimates.  Wind speed is based
    /// on a Rayleigh distribution and interpolation of known data values.
    /// The speed in given in m/s.  The wind direction is based on a Gaussian distribution and
    /// is ultimately converted to degrees, where 0 degrees points directly east.
    /// "pressure_altitude_relationship.txt" - contains corresponding pressures for each altitude
    /// "vmonthly10.txt" through "vmonthly1000.txt" - these 17 data files contain velocity data
    /// values for the v-wind componet for each pressure level.  This data is taken from the website
    /// www.noaa.gov using their NCAR Reanalysis data.
    /// "umonthly10.txt" through "umonthly1000.txt" - these files contain data for the u-component
    /// of the wind.


    char filename[4096];
    static int nsig(17);
    float sigmas[nsig], udata[13][61], mags[61], sqmags[61], meanmag[nsig];
    float uvel[61], vvel[61], stddevs[nsig], avedirs[nsig], dirs[61],
        vdata[13][61];
    float sum = 0.0, sigma = 0.0, u = 0.0, v = 0.0, dir = 0.0, uave = 0.0, vave = 0.0;
    char* pres[] = {(char*)"10",(char*)"20",(char*)"30",(char*)"50",
                    (char*)"70",(char*)"100",(char*)"150",(char*)"200",
                    (char*)"250",(char*)"300",(char*)"400",(char*)"500",
                    (char*)"600", (char*)"700",(char*)"850",(char*)"925",
                    (char*)"1000"};
    int i = 0;
    for (i = 0; i < nsig; i++) {
        sprintf(filename, "%s/atmosphere/vmonthly%s.txt", datadir.c_str(), pres[i]);
        FILE *f = fopen(filename, "r");
        int k = 0;

        for (k = 0; k < 61; k++) {/*reads in v velocities for a given pressure level*/
            fscanf(f,"%f %f %f %f %f %f %f %f %f %f %f %f %f", &vdata[0][k],
                   &vdata[1][k], &vdata[2][k], &vdata[3][k], &vdata[4][k],
                   &vdata[5][k], &vdata[6][k], &vdata[7][k], &vdata[8][k],
                   &vdata[9][k], &vdata[10][k], &vdata[11][k], &vdata[12][k]);
        }
        fclose(f);

        sprintf(filename, "%s/atmosphere/umonthly%s.txt", datadir.c_str(), pres[i]);
        FILE *f2 = fopen(filename, "r");
        for (k = 0; k < 61; k++) {/*reads in u velocities*/
            fscanf(f2, "%f %f %f %f %f %f %f %f %f %f %f %f %f", &udata[0][k],
                   &udata[1][k], &udata[2][k], &udata[3][k], &udata[4][k],
                   &udata[5][k], &udata[6][k], &udata[7][k], &udata[8][k],
                   &udata[9][k], &udata[10][k], &udata[11][k], &udata[12][k]);
        }
        fclose(f2);
        for (k = 0; k < 61; k++){
            // finds the appropriate velocities for the chosen month
            uvel[k] = udata[static_cast<int>(monthnum)][k];
            vvel[k] = vdata[static_cast<int>(monthnum)][k];
        }
        for (k = 0; k < 61; k++) {
            // finds the magnitudes of those velocity components
            mags[k] = sqrt(pow(uvel[k], 2) + pow(vvel[k], 2));
        }
        uave = mean(uvel, 61);
        vave = mean(vvel, 61);
        int j = 0;
        for (j = 0; j < 61; j++) {
            // finds the direction in radians of each velocity component pair*/
            u = uvel[j];
            v = vvel[j];
            if (u > 0.0) {
                if (v > 0.0) {
                    dirs[j] = atan(v/u);
                } else {
                    dirs[j] = atan(v/u) + 2*M_PI;
                }
            } else {
                if (v > 0.0) {
                    dirs[j] = atan(v/u) + M_PI;
                } else {
                    dirs[j] = atan(v/u) + M_PI;
                }
            }
        }
        // finds the standard deviation of the directions
        stddevs[i] = stdev(dirs, 61);
        float dirs2[61];
        float stddevs2[nsig];
        int n = 0;
        for (n = 0; n < 61; n++) {
            dirs2[n] = dirs[n];
        }
        for (n = 0; n < 61; n++) {
            if (dirs[n] <= M_PI) {
                dirs2[n] = dirs2[n] + 2*M_PI;
            }
        }
        stddevs2[i] = stdev(dirs2, 61);
        if (stddevs2[i] < stddevs[i]) {
            // finds the true standard deviation by taking the lower of
            // the two ie. when taken from 0 to 2pi or pi to 3pi
            stddevs[i] = stddevs2[i];
        }
        if (uave > 0.0) {/*finds the direction of the average magnitude*/
            if (vave > 0.0) {
                dir = atan(vave/uave);
            } else {
                dir = atan(vave/uave) + 2*M_PI;
            }
        } else {
            if (vave > 0) {
                dir = atan(vave/uave) + M_PI;
            } else {
                dir = atan(vave/uave) + M_PI;
            }
        }
        avedirs[i] = dir;
        meanmag[i] = mean(mags, 61);/*creates an array of mean magnitudes for each height*/
        int b = 0;
        sum = 0.0;
        for (b = 0; b < 61; b++) {
            sqmags[b] = pow(mags[b], 2);
            sum = sum + sqmags[b];
        }
        sigma = sqrt(sum/122.0);
        sigmas[i] = sigma;
    }
    float hdata[2][nsig];
    float alts[nsig];
    FILE *hfile;
    sprintf(filename, "%s/atmosphere/pressure_altitude_relationship.txt", datadir.c_str());
    hfile = fopen(filename, "r");
    for (i = 0; i < nsig; i++) { /*gets the altitudes for each set pressure level*/
        fscanf(hfile, "%f %f ", &hdata[0][i], &hdata[1][i]);
    }
    fclose(hfile);
    for (i = 0; i < nsig; i++) {
        alts[i] = hdata[1][i];
    }

    // interpolates a sigma value for the desired height
    float mean = interpol(altitude, alts, sigmas, nsig);
    float squareroot = sqrt(-2*log(r1));
    magest = mean*squareroot;

    // interpolates a direction for the desired altitude
    float intdir = interpol(altitude, alts, avedirs, nsig);
    // interpolates a standard deviation for the desired altitude
    float intstd = interpol(altitude, alts, stddevs, nsig);
    direst = intdir + r2*intstd;
    if (direst < 0.0) {
        direst = direst + 2*M_PI;
    }
}

float Atmosphere::outerscale(float altitude, double tai, int n) {
    /// @fn float AtmosphereCreator::outerscale(float altitude)
    //  @brief creates an outerscale estimate based on a lognormal distribution.

    float mean, median;

    readText pars(instrdir + "/site.txt");
    for (size_t t(0); t < pars.getSize(); t++) {
        readText::get(pars[t], "mean", mean);
        readText::get(pars[t], "median", median);
    }

    float mu = log(median);
    float sigma = pow(2.0*(log(mean) - mu), 0.5);
    mu = log(mean) - (1/2)*(pow(sigma, 2.0));
    float osest = exp(mu + sigma*random.normalCorrel(tai + n*1000 + 2000, 600.0/(24.0*3600.0)));
    return osest;
}
