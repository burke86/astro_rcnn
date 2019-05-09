///
/// @package phosim
/// @file silicon.cpp
/// @brief silicon class
///
/// @brief Created by:
/// @author Andy Rasmussen (SLAC)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <istream>
#include <vector>
#include <fitsio.h>
#include <fitsio2.h>

#include "silicon.h"
#include "constants.h"
#include "parameters.h"
#include "ancillary/readtext.h"

double Silicon::absorptionCoeffMCT (double lambda, double temperature, double x) {

    // Hg(1-x)Cd(x)Te follows Chu et al. J. Appl. Phys. 75 1234-5. (Kane)
    // and Finkman and Schacham J. Appl. Phys. 56 2896-900 (Urbach).
    double energy = 12398/lambda; // eV
    double alpha_o = exp(53.61*x - 18.88);
    double E_o = -0.3424 + 1.838*x + 0.148*pow(x, 2);
    double T_o = 81.9; // K
    double sigma = 3.267e4*(1 + x); // K eV-1
    double alpha_T = 100 + 5000*x;
    double E_T = (T_o + temperature)/sigma*log(alpha_T/alpha_o) + E_o; // eV
    double E_g = -0.302 + 1.93*x + 5.35e-4 * temperature * (1 - 2*x) - 0.81 * pow(x, 2) + 0.832*pow(x, 3); // eV  Hanson et al.
    double beta = alpha_T/pow(E_T - E_g, 0.5);
    double alpha = 0.0;

    if (energy < E_g) { // Urbach
        alpha = alpha_o*exp(sigma*(energy - E_o)/(temperature + T_o));
    } else { // Kane region (unsure about best model)
        alpha = beta*pow(energy - E_g, 0.5);
    }

    // alpha is in cm-1	
    return(alpha);

}


double Silicon::absorptionCoeffSi (double lambda, double temperature) {
    
    // follows Rajkanan et al. Solid-state electronics 22, pp793-795.
    double egT[2], eg0[] = {1.1557f, 2.5f}; // eV
    double edgT = 0.0, egd0 = 3.2f;         // eV
    double ep[] = {1.827e-2f, 5.773e-2f};   // eV
    double C[] = {5.5f, 4.0f};              // no dimension
    double A[] = {3.231e2f, 7.237e3f};      // cm-1 eV-2
    double ad = 1.052e6f;                   // cm-1 eV-2
    double kBoltzmann = 8.617e-5f;          // eV K-1
    double eV;                              // eV
    double beta = 7.021e-4f;                // eV K-1
    double gamma = 1108;                    // K
    double alpha = 0.0;                     // cm-1
    double deltaE0, deltaE1[2][2][2];       // eV

    edgT   = egd0   - (beta*temperature*temperature/(temperature + gamma));
    egT[0] = eg0[0] - (beta*temperature*temperature/(temperature + gamma));
    egT[1] = eg0[1] - (beta*temperature*temperature/(temperature + gamma));

    if (lambda > 3100) {
        eV = 12398/lambda;
        deltaE0 = eV - edgT;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            deltaE1[i][j][0] = eV - egT[j] + ep[i];
            deltaE1[i][j][1] = eV - egT[j] - ep[i];
        }
    }

    alpha = ad*sqrt((deltaE0>0)?deltaE0:0);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            alpha += C[i]*A[j]*
                (((deltaE1[i][j][0]>0)?pow(deltaE1[i][j][0], 2)/(exp(ep[i]/(kBoltzmann*temperature)) - 1):0) +
                 ((deltaE1[i][j][1]>0)?pow(deltaE1[i][j][1], 2)/(1 - exp(-ep[i]/(kBoltzmann*temperature))):0));
        }
    }
    // alpha is in cm-1.
    } else {
    // wavelength range is between 413 & 3100 Angstroms - not well characterized. use a constant.
        alpha = 1.35e6;
    }
    return(alpha);
}


double Silicon::indexRefractionMCT (double lambda, double temperature, double x) {
    
    // follows Liu et al. J. Appl. Phys. 75 4176.
    double A = 13.173 - 9.852*x + 2.909*pow(x, 2) + 0.001*(300 - temperature);
    double B = 0.83 - 0.246*x - 0.0961*pow(x, 2) + 8e-4 *(300 - temperature);
    double C = 6.706 - 14.437*x + 8.531*pow(x, 2) + 7e-4*(300 - temperature);
    double D = 1.953e-4 - 0.00128*x + 1.853e-4*pow(x, 2);
    double n = pow(A + B/(1 - pow(C/lambda, 2)) + D * pow(lambda, 2), 0.5);
    
    return(n);
}


double Silicon::indexRefractionSi (double lambda) {

    double energy = 12398/lambda;
    // data points were digitized from Philipp & Taft (1960) and fit between 0 and 3.4 eV
    return(double)(3.364967 + 0.2119184*energy + 2.78878*exp((energy - 3.30168)/0.397862));

}


double Silicon::dopeProfile (double z, double nbulk, double nf, double nb, double sf, double sb, double tsi) {

    return(nbulk + nb*exp(-(tsi - z)/sb) + nf*exp(-(z)/sf));

}

double Silicon::muMCT(double x, double temperature) {

    // follows Rosbeck et al. (1982)
    double s = pow(0.2/x, 7.5);
    double r = pow(0.2/x, 0.6);
    double mu_e = 9e8*s/pow(temperature ,2*r);
    return(mu_e/100); //mu_h

}


double Silicon::muSi (double efield, double temperature, int polarity) {

    // Jacobini et al. (1977) equation (9)
    double vm, ec, beta;
    if (polarity == -1) {
        vm = 1.53e9 * pow(temperature, -0.87); // cm/s
        ec = 1.01 * pow(temperature, 1.55); // V/cm
        beta = 2.57e-2 * pow(temperature, 0.66); // index
    } else {
        vm = 1.62e8 * pow(temperature, -0.52); // cm/s
        ec = 1.24 * pow(temperature, 1.68); // V/cm
        beta = 0.46 * pow(temperature, 0.17); // index
    }
    return((vm/ec)/pow(1 + pow(fabs(efield)/ec, beta), 1/beta));

}


double Silicon::epsilonMCT(double x) {

    // follows Dornhaus et al. (1983)
    // high frequency approximation
    double epsilon_inf = 15.2 - 15.6*x + 8.2*pow(x,2);
    return(epsilon_inf);

}


void Silicon::setup (std::string devmaterial, double sensorTemp, double nbulk, double nf, double nb,
                     double sf, double sb, double tsi, double overdepBias,
                     std::string instrdir,
                     long nampx, long nampy, double pixsize, long seedchip,
                     float *impurityX, float *impurityY, int impurityvariation,
                     double minwavelength, double maxwavelength) {

    sf = sf/1e4f;
    sb = sb/1e4f;
    tsi = tsi/1e4f;

    numWavelength = 1024;
    numTemperature = 16;
    numDopant = 16;
    numThickness = SILICON_STEPS;

    wavelengthGrid = static_cast<double*>(calloc(numWavelength, sizeof(double)));
    if (wavelengthGrid == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    temperatureGrid = static_cast<double*>(calloc(numTemperature, sizeof(double)));
    if (temperatureGrid == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    rho = static_cast<double*>(calloc(numTemperature, sizeof(double)));
    if (rho == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    dopantGrid = static_cast<double*>(calloc(numDopant, sizeof(double)));
    if (dopantGrid == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    thicknessGrid = static_cast<double*>(calloc(numThickness, sizeof(double)));
    if (thicknessGrid == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    meanFreePath = static_cast<double*>(calloc(numWavelength*numTemperature, sizeof(double)));
    if (meanFreePath == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    indexRefraction = static_cast<double*>(calloc(numWavelength, sizeof(double)));
    if (indexRefraction == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    sigma = static_cast<float*>(calloc(numTemperature*numThickness*numDopant, sizeof(float)));
    if (sigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    hsigma = static_cast<float*>(calloc(numTemperature*numThickness*numDopant, sizeof(float)));
    if (hsigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    isigma = static_cast<float*>(calloc(numTemperature*numThickness*numDopant, sizeof(float)));
    if (isigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    fsigma = static_cast<float*>(calloc(numThickness*numDopant, sizeof(float)));
    if (fsigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    gsigma = static_cast<float*>(calloc(numThickness*numDopant, sizeof(float)));
    if (gsigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    sigmaX = static_cast<float*>(calloc(nampx*nampy, sizeof(float)));
    if (sigmaX == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    sigmaY = static_cast<float*>(calloc(nampx*nampy, sizeof(float)));
    if (sigmaY == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    gammaX = static_cast<float*>(calloc(nampx*nampy, sizeof(float)));
    if (gammaX == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    gammaY = static_cast<float*>(calloc(nampx*nampy, sizeof(float)));
    if (gammaY == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    deltaX = static_cast<float*>(calloc(nampx*nampy, sizeof(float)));
    if (deltaX == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    deltaY = static_cast<float*>(calloc(nampx*nampy, sizeof(float)));
    if (deltaY == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    nbulkmap = static_cast<float*>(calloc(nampx*nampy, sizeof(float)));
    if (nbulkmap == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    periodmap = static_cast<float*>(calloc(nampx*nampy, sizeof(float)));
    if (periodmap == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    deadLayer = static_cast<float*>(calloc(nampx*nampy, sizeof(float)));
    if (deadLayer == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    
    double compositionHgCd = 0.5;
    double period = 100.0;
    double amplitude = 0.0;
    double ratio = 10.0;
    double surfaceCharge;
    double depth = 0.0;
    double width = 3000.0;
    double height = 500.0;
    double xoverlap = 0.0;
    double yoverlap = 0.0;
    double pixelVarX = 0.0;
    double pixelVarY = 0.0;
    double epsilon;
    double periodVariation = 0.0;
    double periodDecay = 100000.0;
    int polarity = -1;
    std::vector<double> xinit, yinit, angle;
    spaceChargeShield = 0.0;
    spaceChargeSpreadX = 0.0;
    spaceChargeSpreadY = 0.0;
    polarity = -1;
    chargeStopCharge = 0.0;

    std::string sss;
    sss = instrdir + "/sensor.txt";
    std::ifstream inStream(sss.c_str());
    if (inStream) {
        readText pars(instrdir + "/sensor.txt");
        for (size_t t(0); t < pars.getSize(); t++) {
            std::string line(pars[t]);
            readText::get(line, "compositionHgCd", compositionHgCd);
            readText::get(line, "treeRingPeriod", period);
            readText::get(line, "treeRingVariation", periodVariation);
            readText::get(line, "treeRingDecay", periodDecay);
            readText::get(line, "treeRingAmplitude", amplitude);
            readText::get(line, "treeRingRatio", ratio);
            readText::get(line, "edgeSurfaceCharge", surfaceCharge);
            readText::get(line, "deadLayerDepth", depth);
            readText::get(line, "deadLayerWidth", width);
            readText::get(line, "deadLayerHeight", height);
            readText::get(line, "deadLayerXoverlap", xoverlap);
            readText::get(line, "deadLayerYoverlap", yoverlap);
            readText::get(line, "pixelVarX", pixelVarX);
            readText::get(line, "pixelVarY", pixelVarY);
            readText::get(line, "deadLayerXinit", xinit);
            readText::get(line, "deadLayerYinit", yinit);
            readText::get(line, "deadLayerAngle", angle);
            readText::get(line, "spaceChargeShield", spaceChargeShield);
            readText::get(line, "spaceChargeSpreadX", spaceChargeSpreadX);
            readText::get(line, "spaceChargeSpreadY", spaceChargeSpreadY);
            readText::get(line, "chargeStopCharge", chargeStopCharge);
            readText::get(line, "siliconType", polarity);
        }
    }

    spaceChargeSpreadX = (spaceChargeSpreadX/pixsize)*(spaceChargeSpreadX/pixsize);
    spaceChargeSpreadY = (spaceChargeSpreadY/pixsize)*(spaceChargeSpreadY/pixsize);
    period = period/pixsize;
    periodDecay = periodDecay/pixsize;
    periodVariation = periodVariation/pixsize;
    
    // tree ring map
    random.setSeed32Fixed(1000 + seedchip);
    *impurityX = random.uniformFixed()*(ratio*nampx) - (ratio - 1)/2*nampx;
    *impurityY = random.uniformFixed()*(ratio*nampy) - (ratio - 1)/2*nampy;
    for (long k = 0; k < 10; k++) {
        double phase = random.uniformFixed()*2.0*PI;
        double localVariation = periodVariation*random.normalFixed();
        for (long i = 0; i < nampx; i++) {
            for (long j = 0; j < nampy; j++) {
                double radius = sqrt(pow(i - *impurityX, 2) + pow(j - *impurityY, 2));
                double localPeriod = (period + localVariation)*(exp(-radius/periodDecay));
                periodmap[nampx*j + i] += localPeriod/10.0;
                if (k == 0) nbulkmap[nampx*j + i] = 1.0 + amplitude/sqrt(10.0)*sin(2*PI*radius/localPeriod + phase);
                if (k != 0) nbulkmap[nampx*j + i] += amplitude/sqrt(10.0)*sin(2*PI*radius/localPeriod + phase);

            }
        }
    }

    if (devmaterial == "MCT") epsilon = epsilonMCT(compositionHgCd);
    else epsilon = EPSILON_SI;

    // lateral deflection map
    double scale = fabs(surfaceCharge)/nbulk;
    if (scale == 0.0) scale = 1e-4*pixsize;
    for (long i = 0; i < nampx; i++) {
        for (long j = 0; j < nampy; j++) {
            deltaX[nampx*j + i] = surfaceCharge/EPSILON_0/epsilon*exp(-(i*pixsize*1e-4)/scale)*E_CHARGE-
                surfaceCharge/EPSILON_0/EPSILON_SI*exp(-((nampx - 1 - i)*pixsize*1e-4)/scale)*E_CHARGE;
            deltaY[nampx*j + i] = surfaceCharge/EPSILON_0/epsilon*exp(-(j*pixsize*1e-4)/scale)*E_CHARGE-
                surfaceCharge/EPSILON_0/EPSILON_SI*exp(-((nampy - 1 - j)*pixsize*1e-4)/scale)*E_CHARGE;
        }
    }

    // add tree rings
    for (long i = 1; i < (nampx - 1); i++) {
        for (long j = 1; j < (nampy - 1); j++) {
            sigmaX[nampx*j + i] = (nbulkmap[nampx*j + i + 1] - nbulkmap[nampx*j + i - 1])*
                nbulk*(2.0*tsi)*(periodmap[nampx*j + i]*pixsize*1e-4)/(2.0*PI)/(2.0*PI)*E_CHARGE/
                epsilon/EPSILON_0/2.0/(pixsize*1e-4)/2.0;
            sigmaY[nampx*j + i] = (nbulkmap[nampx*(j + 1) + i] - nbulkmap[nampx*(j - 1) + i])*
                nbulk*(2.0*tsi)*(periodmap[nampx*j + i]*pixsize*1e-4)/(2.0*PI)/(2.0*PI)*E_CHARGE/
                epsilon/EPSILON_0/2.0/(pixsize*1e-4)/2.0;
        }
    }

    if (impurityvariation == 0) {
        for (long i = 0; i < nampx; i++) {
            for (long j = 0; j < nampy; j++) {
                nbulkmap[nampx*j + i] = 1.0;
            }
        }
    }

    // int status;
    // long naxesa[2];
    // fitsfile *faptr;
    // status = 0;
    // fits_create_file(&faptr,"!defl.fits",&status);
    // naxesa[0] = nampx; naxesa[1] = nampy;
    // fits_create_img(faptr,FLOAT_IMG,2,naxesa,&status);
    // fits_write_img(faptr,TFLOAT,1,nampx*nampy,deltaX,&status);
    // fits_close_file(faptr,&status);

    // nonuniform pixels
    for (long i = 0; i < nampx; i++) {
        for (long j = 0; j < nampy; j++) {
            gammaX[nampx*j + i] = random.normalFixed()*pixelVarX*pixsize;
            gammaY[nampx*j + i] = random.normalFixed()*pixelVarY*pixsize;
        }
    }

    // dead layer
    for (size_t k = 0; k < xinit.size(); k++) {
        for (long i = 0; i < nampx; i++) {
            for (long j = 0; j < nampy; j++) {
                double ip = i*cos(angle[k]*DEGREE) - j*sin(angle[k]*DEGREE);
                double jp = i*sin(angle[k]*DEGREE) + j*cos(angle[k]*DEGREE);
                double di = fmod(ip*pixsize - xinit[k], width - xoverlap);
                double dj = fmod(jp*pixsize - yinit[k], height - yoverlap);
                while (di<0) di += width - xoverlap;
                while (dj<0) dj += height - yoverlap;
                if (di < width && dj < height) {
                    deadLayer[nampx*j + i] += depth;
                    if (di < xoverlap) deadLayer[nampx*j + i] += depth;
                    if (dj < yoverlap) deadLayer[nampx*j + i] += depth;
                    if (di < xoverlap && dj < yoverlap) deadLayer[nampx*j + i] += depth;
                }
            }
        }
    }

    if (devmaterial == "MCT") { // HgCdTe
        for (long i = 0; i < numTemperature; i++) {
            for (long j = 0; j < numWavelength; j++) {
                temperatureGrid[i] = sensorTemp + (static_cast<double>(i))/(static_cast<double>(numTemperature))*2.0 - 1.0;
                rho[i] = tsi*(static_cast<double>(i + 1))/(static_cast<double>(numTemperature))/5.0;
                wavelengthGrid[j] = minwavelength/1000.0 +
                    (static_cast<double>(j))/(static_cast<double>(numWavelength - 1))*(maxwavelength - minwavelength)/1000.0;
                meanFreePath[i*numWavelength + j] = 10.0/absorptionCoeffMCT(wavelengthGrid[j]*10000.0, temperatureGrid[i], compositionHgCd);
            }
        }
        for (long i = 0; i < numWavelength; i++) {
            indexRefraction[i] = indexRefractionMCT(wavelengthGrid[i], temperatureGrid[i], compositionHgCd);
        }
    } else { // Silicon
        for (long i = 0; i < numTemperature; i++) {
            for (long j = 0; j < numWavelength; j++) {
                temperatureGrid[i] = sensorTemp + (static_cast<double>(i))/(static_cast<double>(numTemperature))*2.0 - 1.0;
                rho[i] = tsi*(static_cast<double>(i + 1))/(static_cast<double>(numTemperature))/5.0;
                wavelengthGrid[j] = minwavelength/1000.0 +
                    (static_cast<double>(j))/(static_cast<double>(numWavelength - 1))*(maxwavelength - minwavelength)/1000.0;
                meanFreePath[i*numWavelength + j] = 10.0/absorptionCoeffSi(wavelengthGrid[j]*10000.0, temperatureGrid[i]);
            }
        }
        for (long i = 0; i < numWavelength; i++) {
            indexRefraction[i] = indexRefractionSi(wavelengthGrid[i]*10000.0);
        }
    }
    for (long i = 0; i < numDopant; i++) {
        dopantGrid[i] = (static_cast<double>(i))/(static_cast<double>(numDopant))*2*amplitude*nbulk +
            (1.0 - amplitude)*nbulk;
    }

    double z[SILICON_STEPS], E[SILICON_STEPS], v[SILICON_STEPS],
        tcol[SILICON_STEPS], tp[SILICON_STEPS], gp[SILICON_STEPS], hp[SILICON_STEPS], ip[SILICON_STEPS],
        tcoli[SILICON_STEPS], tpi[SILICON_STEPS], gpi[SILICON_STEPS], hpi[SILICON_STEPS], ipi[SILICON_STEPS];

    double minE, dz;
    int i, minEindex;
    long n = SILICON_STEPS;

    if (devmaterial == "MCT") { // MCT

        double factor = E_CHARGE/EPSILON_0/epsilonMCT(compositionHgCd)/(4*PI)*spaceChargeShield;

        for (long l = 0; l < numDopant; l++) {
            for (long k = 0; k < numTemperature; k++) {
                dz = (-tsi)/static_cast<double>(n - 1);
                for (i = 0; i < n; i++) {
                    z[i] = tsi + i*dz;
                    thicknessGrid[i] = tsi + i*dz;
                    if (i == 0) {
                        E[0] = overdepBias/tsi;
                        v[0] = 0;
                    } else {
                        // sample the doping profile to evaluate E[i]
                        int nSample = 100;
                        double dp = 0;
                        int j;
                        dp = 0;
                        for (j = 0; j < nSample; j++) {
                            dp += dopeProfile(tsi + dz*(i + (j)/(static_cast<double>(nSample))), dopantGrid[l], nf, nb, sf, sb, tsi);
                        }
                        dp /= nSample;
                        E[i] = E[i - 1] + dp*E_CHARGE/(EPSILON_0*epsilonMCT(compositionHgCd))*dz;
                        v[i] = v[i - 1] - 0.5f*(E[i] + E[i - 1])*dz;
                    }
                }
                minEindex = 0;
                minE = 0;
                for (i = 0; i < n; i++) {
                    // find minimum value for E and count from that index
                    if (E[i] < minE) {
                        minEindex = i;
                        minE = E[i];
                    }
                }
                i = n - 1;
                tcoli[i] = (dz/(minE*muMCT(compositionHgCd, temperatureGrid[k])));
                tpi[i] = dz/minE*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                    (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                gpi[i] = dz/minE;
                hpi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                ipi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                while (i--) {
                    if (i > minEindex) {
                        tcoli[i] = (dz/(minE*muMCT(compositionHgCd, temperatureGrid[k])));
                        tpi[i] = dz/minE*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                            (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                        gpi[i] = dz/minE;
                        hpi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        ipi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                    } else {
                        if (E[i] + E[i + 1] > 0) {
                            tcoli[i] = 0.1;
                            tpi[i] = 0.1*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                                (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                            gpi[i] = 0.1;
                            hpi[i] = 0.1*factor/(z[i]*z[i] + rho[k]*rho[k]);
                            ipi[i] = 0.1*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        } else {
                            tcoli[i] = 2*dz/(E[i]*muMCT(compositionHgCd, temperatureGrid[k]) +
                                            E[i + 1]*muMCT(compositionHgCd, temperatureGrid[k]));
                            tpi[i] = 2*dz/(E[i] + E[i + 1])*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                                (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                            gpi[i] = 2*dz/(E[i] + E[i + 1]);
                            hpi[i] = 2*dz/(E[i] + E[i + 1])*factor/(z[i]*z[i] + rho[k]*rho[k]);
                            ipi[i] = 2*dz/(E[i] + E[i + 1])*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        }
                    }
                }
                i = n;
                while (i--) {
                    tcol[i] = 0.0;
                    tp[i] = 0.0;
                    gp[i] = 0.0;
                    hp[i] = 0.0;
                    ip[i] = 0.0;
                    for (long j = i; j < i + SILICON_SUB_STEPS; j++) {
                        if (j < n) {
                            tcol[i] += tcoli[j];
                            tp[i] += tpi[j];
                            gp[i] += gpi[j];
                            hp[i] += hpi[j];
                            ip[i] += ipi[j];
                        }
                    }
                }
                i = n;
                while (i--) {
                    sigma[l*numThickness*numTemperature + k*numThickness + i] = sqrt(2*K_BOLTZMANN*temperatureGrid[k]
                                                                                     /E_CHARGE*muMCT(compositionHgCd, temperatureGrid[k]) * tcol[i])*1e4f;
                    fsigma[l*numThickness + i] = tp[i]*1e4f;
                    gsigma[l*numThickness + i] = gp[i]*1e4f;
                    hsigma[l*numThickness*numTemperature + k*numThickness + i] = hp[i]*1e4f;
                    isigma[l*numThickness*numTemperature + k*numThickness + i] = ip[i]*1e4f;
                }

            }
        }
    } else { // Silicon
        
        double factor = E_CHARGE/EPSILON_0/EPSILON_SI/(4*PI)*spaceChargeShield;

        for (long l = 0; l < numDopant; l++) {
            for (long k = 0; k < numTemperature; k++) {
                dz = (-tsi)/static_cast<double>(n - 1);
                for (i = 0; i < n; i++) {
                    z[i] = tsi + i*dz;
                    thicknessGrid[i] = tsi + i*dz;
                    if (i == 0) {
                        E[0] = overdepBias/tsi;
                        v[0] = 0;
                    } else {
                        // sample the doping profile to evaluate E[i]
                        int nSample = 100;
                        double dp = 0;
                        int j;
                        dp = 0;
                        for (j = 0; j < nSample; j++) {
                            dp += dopeProfile(tsi + dz*(i + (j)/(static_cast<double>(nSample))), dopantGrid[l], nf, nb, sf, sb, tsi);
                        }
                        dp /= nSample;
                        E[i] = E[i - 1] + dp*E_CHARGE/(EPSILON_0*EPSILON_SI)*dz;
                        v[i] = v[i - 1] - 0.5f*(E[i] + E[i - 1])*dz;
                    }
                }
                minEindex = 0;
                minE = 0;
                for (i = 0; i < n; i++) {
                    // find minimum value for E and count from that index
                    if (E[i] < minE) {
                        minEindex = i;
                        minE = E[i];
                    }
                }
                i = n - 1;
                tcoli[i] = (dz/(minE*muSi(minE, temperatureGrid[k], polarity)));
                tpi[i] = dz/minE*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                    (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                gpi[i] = dz/minE;
                hpi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                ipi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                while (i--) {
                    if (i > minEindex) {
                        tcoli[i] = (dz/(minE*muSi(minE, temperatureGrid[k], polarity)));
                        tpi[i] = dz/minE*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                            (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                        gpi[i] = dz/minE;
                        hpi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        ipi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                    } else {
                        if (E[i] + E[i + 1] > 0) {
                            tcoli[i] = 0.1;
                            tpi[i] = 0.1*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                                (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                            gpi[i] = 0.1;
                            hpi[i] = 0.1*factor/(z[i]*z[i] + rho[k]*rho[k]);
                            ipi[i] = 0.1*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        } else {
                            tcoli[i] = 2*dz/(E[i]*muSi(E[i], temperatureGrid[k], polarity) +
                                            E[i + 1]*muSi(E[i + 1], temperatureGrid[k], polarity));
                            tpi[i] = 2*dz/(E[i] + E[i + 1])*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                                (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                            gpi[i] = 2*dz/(E[i] + E[i + 1]);
                            hpi[i] = 2*dz/(E[i] + E[i + 1])*factor/(z[i]*z[i] + rho[k]*rho[k]);
                            ipi[i] = 2*dz/(E[i] + E[i + 1])*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        }
                    }
                }
                i = n;
                while (i--) {
                    tcol[i] = 0.0;
                    tp[i] = 0.0;
                    gp[i] = 0.0;
                    hp[i] = 0.0;
                    ip[i] = 0.0;
                    for (long j = i; j < i + SILICON_SUB_STEPS; j++) {
                        if (j < n) {
                            tcol[i] += tcoli[j];
                            tp[i] += tpi[j];
                            gp[i] += gpi[j];
                            hp[i] += hpi[j];
                            ip[i] += ipi[j];
                        }
                    }
                }
                i = n;
                while (i--) {
                    sigma[l*numThickness*numTemperature + k*numThickness + i] = sqrt(2*K_BOLTZMANN*temperatureGrid[k]
                                                                                     /E_CHARGE*muSi(0, temperatureGrid[k], polarity) * tcol[i])*1e4f;
                    fsigma[l*numThickness + i] = tp[i]*1e4f;
                    gsigma[l*numThickness + i] = gp[i]*1e4f;
                    if (isnan(tp[i]) || isnan(gp[i])) printf("A %ld %d %e %e ",l,i,tp[i],gp[i]);
                    hsigma[l*numThickness*numTemperature + k*numThickness + i] = hp[i]*1e4f;
                    isigma[l*numThickness*numTemperature + k*numThickness + i] = ip[i]*1e4f;
                }

            }
        }
    }

}
