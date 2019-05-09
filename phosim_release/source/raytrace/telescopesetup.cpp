///
/// @package phosim
/// @file telescopesetup.cpp
/// @brief setup for telescope (part of image class)
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

#include <vector>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include "ancillary/readtext.h"
#include "fea.h"

using namespace fea;

int Image::telSetup () {

    double runningz;

    long seedchip = 0;
    for (size_t m = 0; m < chipid.size(); m++) {
        seedchip += static_cast<long>((static_cast<int>(chipid.c_str()[m]%10))*pow(10, m));
    }

    chip.nampx = (maxx - minx + 1);
    chip.nampy = (maxy - miny + 1);
    chip.buffer = 400;
    chip.midpoint = pixelsx/2;
    if (saturation == 1) {

        state.satupmap = static_cast<std::atomic <int>*>(calloc(chip.nampx*chip.nampy, sizeof(int)));
        state.satdownmap = static_cast<std::atomic <int>*>(calloc(chip.nampx*chip.nampy, sizeof(int)));
        for (long i = minx; i <= maxx; i++) {
            for (long j = miny; j <= maxy; j++) {
                *(state.satdownmap + chip.nampx*(j - miny) + (i - minx)) = i - 1;
                *(state.satupmap + chip.nampx*(j - miny) + (i - minx)) = i + 1;
            }
        }

    }
    state.focal_plane = static_cast<std::atomic<unsigned long>*>(calloc(chip.nampx*chip.nampy, sizeof(std::atomic<unsigned long>)));
    state.focal_plane_fl = static_cast<float*>(calloc(chip.nampx*chip.nampy, sizeof(float)));
    if (opdfile) {
        state.opd = static_cast<double*>(calloc(nsource*OPD_SCREEN_SIZE*OPD_SCREEN_SIZE, sizeof(double)));
        state.opdcount = static_cast<double*>(calloc(nsource*OPD_SCREEN_SIZE*OPD_SCREEN_SIZE, sizeof(double)));
        state.cx = static_cast<double*>(calloc(nsource, sizeof(double)));
        state.cy = static_cast<double*>(calloc(nsource, sizeof(double)));
        state.cz = static_cast<double*>(calloc(nsource, sizeof(double)));
        state.r0 = static_cast<double*>(calloc(nsource, sizeof(double)));
        state.epR = static_cast<double*>(calloc(nsource, sizeof(double)));
    }

    // OPTICS AND COATINGS
    fprintf(stdout, "Building Optics.\n");
    std::ostringstream opticsFile;
    opticsFile << instrdir << "/optics_" << filter << ".txt";
    readText opticsPars(opticsFile.str());
    long totSurf = opticsPars.getSize() + 1;  //detector, focalplane, exit pupil

    surface.setup(totSurf, SURFACE_POINTS);
    coating.setup(totSurf);

    nsurf = 0;
    runningz = 0.0;
    nmirror = 0;
    for (long t = 0; t < totSurf - 1; t++){
        std::istringstream iss(opticsPars[t]);
        std::string surfaceName, surfacetype, coatingFile, mediumFile;
        iss >> surfaceName >> surfacetype;

        if (surfacetype == "mirror") surface.surfacetype[nsurf] = MIRROR;
        else if (surfacetype == "lens") surface.surfacetype[nsurf] = LENS;
        else if (surfacetype == "filter") surface.surfacetype[nsurf] = FILTER;
        else if (surfacetype == "det") surface.surfacetype[nsurf] = DETECTOR;
        else if (surfacetype == "grating") surface.surfacetype[nsurf] = GRATING;
        else if (surfacetype == "exitpupil") {
            nsurf++;
            surface.surfacetype[nsurf] = EXITPUPIL;
        }
        if (surfacetype == "mirror") nmirror++;

        iss >> surface.radiusCurvature[nsurf];
        double dz;
        iss >> dz;
        runningz += dz;
        surface.height[nsurf] = runningz;
        iss >> surface.outerRadius[nsurf];
        iss >> surface.innerRadius[nsurf];
        iss >> surface.conic[nsurf];
        iss >> surface.three[nsurf];
        iss >> surface.four[nsurf];
        iss >> surface.five[nsurf];
        iss >> surface.six[nsurf];
        iss >> surface.seven[nsurf];
        iss >> surface.eight[nsurf];
        iss >> surface.nine[nsurf];
        iss >> surface.ten[nsurf];
        iss >> coatingFile;
        iss >> mediumFile;
        surface.centerx[nsurf] = 0.0;
        surface.centery[nsurf] = 0.0;
        surface.rmax[nsurf] = surface.outerRadius[nsurf];
        surface.innerRadius0[nsurf] = surface.innerRadius[nsurf];
        if (opdfile) surface.innerRadius[nsurf] = 0.0;
        if (surfacetype != "none") {
            surface.asphere(nsurf, SURFACE_POINTS);
            if (surface.surfacetype[nsurf] == MIRROR || surface.surfacetype[nsurf] == DETECTOR ||
                surface.surfacetype[nsurf] == LENS || surface.surfacetype[nsurf] == FILTER) {
                perturbation.zernikeflag.push_back(zernikemode);
            } else if (surface.surfacetype[nsurf] == EXITPUPIL) {
                perturbation.zernikeflag.push_back(0);
                perturbation.zernikeflag.push_back(0);
            } else {
                perturbation.zernikeflag.push_back(0);
            }

            // COATINGS
            surface.surfacecoating[nsurf] = 0;
            if (coatingFile != "none") {
                readText coatingPars(instrdir + "/" + coatingFile);
                size_t nline = coatingPars.getSize();
                coating.allocate(nsurf, nline);
                long j = 0;
                long i = 0;
                double angle, angle0 = 0.0;
                for (size_t tt(0); tt<nline; tt++) {
                    std::istringstream isst(coatingPars[tt]);
                    isst >> angle;
                    if (i == 0 && j == 0) angle0 = angle;
                    if (angle > angle0) {
                        i++;
                        coating.wavelengthNumber[nsurf] = j;
                        if (j != coating.wavelengthNumber[nsurf] && i > 1) {
                            std::cout << "Error in format of " << coatingFile << std::endl;
                            exit(1);
                        }
                        j = 0;
                    }
                    *(coating.angle[nsurf] + i) = angle;
                    isst >> *(coating.wavelength[nsurf] + j);
                    isst >> *(coating.transmission[nsurf] + tt);
                    isst >> *(coating.reflection[nsurf] + tt);
                    j++;
                    angle0 = angle;
                }
                i++;
                coating.wavelengthNumber[nsurf] = j;
                coating.angleNumber[nsurf] = i;
                coating.wavelength[nsurf] = static_cast<double*>(realloc(coating.wavelength[nsurf], coating.wavelengthNumber[nsurf]*sizeof(double)));
                coating.angle[nsurf] = static_cast<double*>(realloc(coating.angle[nsurf], coating.angleNumber[nsurf]*sizeof(double)));
                surface.surfacecoating[nsurf] = 1;
            }

            // MEDIUM
            surface.surfacemed[nsurf] = 0;
            if (mediumFile == "air") surface.surfacemed[nsurf] = 2;
            if (mediumFile != "vacuum" && mediumFile != "air") {
                medium.setup(nsurf, instrdir + "/" + mediumFile);
                surface.surfacemed[nsurf] = 1;
            }
            nsurf++;

        }
    }
    // if (opdfile) nsurf -= 2; //exit pupil and focalplane

    if (telescopeMode == 1) {
        maxr = surface.outerRadius[0];
        minr = surface.innerRadius0[0];
    }


    // OBSTRUCTIONS
    obstruction.nspid = 0;
    if (diffractionMode >= 1) {
        std::string sss;
        sss = instrdir + "/spider.txt";
        std::ifstream inStream(sss.c_str());
        if (inStream) {
            readText spiderPars(instrdir + "/spider.txt");
            std::cout << "Placing Obstructions." << std::endl;
            for (size_t t(0); t < spiderPars.getSize(); t++) {
                std::istringstream iss(spiderPars[t]);
                std::string spiderType;
                iss >> spiderType;
                if (spiderType == "pupilscreen"){
                    pupilscreenMode = 1;
                } else {
                    iss >> obstruction.height[obstruction.nspid];
                    iss >> obstruction.width[obstruction.nspid];
                    iss >> obstruction.center[obstruction.nspid];
                    double tempf4;
                    iss >> tempf4;
                    obstruction.depth[obstruction.nspid] = tempf4;
                    obstruction.angle[obstruction.nspid] = tempf4;
                    obstruction.reference[obstruction.nspid] = tempf4;
                    if (spiderType == "crossx") obstruction.type[obstruction.nspid] = 1;
                    if (spiderType == "crossy") obstruction.type[obstruction.nspid] = 2;
                    obstruction.nspid++;
                }
            }
        }
    }

    // PERTURBATIONS
    fprintf(stdout, "Perturbing Design.\n");
    perturbation.zernike_coeff = static_cast<double*>(calloc(NTERM*nsurf, sizeof(double)));
    // chuck's old model: 0.82312e-3*pow(zv, -1.2447)

    for (long i = 0; i < nsurf; i++) {
        for (long j = 0; j < NTERM; j++) {
            *(perturbation.zernike_coeff + i*NTERM + j) = izernike[i][j];
        }
    }


    perturbation.eulerPhi.reserve(totSurf);
    perturbation.eulerPsi.reserve(totSurf);
    perturbation.eulerTheta.reserve(totSurf);
    perturbation.decenterX.reserve(totSurf);
    perturbation.decenterY.reserve(totSurf);
    perturbation.defocus.reserve(totSurf);
    perturbation.rmax.reserve(totSurf);

    double *zernikeR, *zernikePhi, *zernikeNormalR, *zernikeNormalPhi;
    double *chebyshev, *chebyshevR, *chebyshevPhi;
    zernikeR = static_cast<double*>(calloc(NZERN*PERTURBATION_POINTS, sizeof(double)));
    zernikePhi = static_cast<double*>(calloc(NZERN*PERTURBATION_POINTS, sizeof(double)));
    zernikeNormalR = static_cast<double*>(calloc(NZERN*PERTURBATION_POINTS, sizeof(double)));
    zernikeNormalPhi = static_cast<double*>(calloc(NZERN*PERTURBATION_POINTS, sizeof(double)));
    chebyshev = static_cast<double*>(calloc(NCHEB*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    chebyshevR = static_cast<double*>(calloc(NCHEB*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    chebyshevPhi = static_cast<double*>(calloc(NCHEB*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_r_grid = static_cast<double*>(calloc(PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_phi_grid = static_cast<double*>(calloc(PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_summed = static_cast<double*>(calloc(nsurf*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_summed_nr_p = static_cast<double*>(calloc(nsurf*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_summed_np_r = static_cast<double*>(calloc(nsurf*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));

    zernikes(zernikeR, zernikePhi, perturbation.zernike_r_grid,
             perturbation.zernike_phi_grid, zernikeNormalR, zernikeNormalPhi,
             PERTURBATION_POINTS, NZERN);
    chebyshevs (perturbation.zernike_r_grid, perturbation.zernike_phi_grid,
                chebyshev, chebyshevR, chebyshevPhi, PERTURBATION_POINTS, NCHEB);

    for (long k = 0; k < nsurf; k++) {
        if (pertType[k] == "zern") {
            for (long j = 0; j < PERTURBATION_POINTS; j++) {
                for (long l = 0; l < PERTURBATION_POINTS; l++) {
                    *(perturbation.zernike_summed + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) = 0;
                    *(perturbation.zernike_summed_nr_p + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) = 0;
                    *(perturbation.zernike_summed_np_r + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) = 0;
                    for (long m = 0; m < NZERN; m++) {
                        *(perturbation.zernike_summed + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) +=
                            *(perturbation.zernike_coeff + k*NTERM + m)*(*(zernikeR + m*PERTURBATION_POINTS + j))*
                            (*(zernikePhi + m*PERTURBATION_POINTS + l));
                        *(perturbation.zernike_summed_nr_p + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) +=
                            *(perturbation.zernike_coeff + k*NTERM + m)*(*(zernikeNormalR + m*PERTURBATION_POINTS + j))*
                            (*(zernikePhi + m*PERTURBATION_POINTS + l));

                        *(perturbation.zernike_summed_np_r + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) +=
                            *(perturbation.zernike_coeff + k*NTERM + m)*(*(zernikeR + m*PERTURBATION_POINTS + j))*
                            (*(zernikeNormalPhi + m*PERTURBATION_POINTS + l));

                    }
                }
            }
        } else if (pertType[k] == "chebyshev") {
            for (long j = 0; j < PERTURBATION_POINTS; j++) {
                for (long l = 0; l < PERTURBATION_POINTS; l++) {
                    long idx = k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j;
                    perturbation.zernike_summed[idx] = 0;
                    perturbation.zernike_summed_nr_p[idx] = 0;
                    perturbation.zernike_summed_np_r[idx] = 0;
                    for (long m = 0; m < NCHEB; m++) {
                        long idxi = m*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j;
                        perturbation.zernike_summed[idx] += perturbation.zernike_coeff[k*NTERM + m]*chebyshev[idxi];
                        perturbation.zernike_summed_nr_p[idx] += perturbation.zernike_coeff[k*NTERM + m]*chebyshevR[idxi];
                        perturbation.zernike_summed_np_r[idx] += perturbation.zernike_coeff[k*NTERM + m]*chebyshevPhi[idxi];
                    }
                }
            }
        }
    }

    perturbation.rotationmatrix = static_cast<double*>(calloc(totSurf*3*3, sizeof(double)));
    perturbation.inverserotationmatrix = static_cast<double*>(calloc(totSurf*3*3, sizeof(double)));
    for (long i = 0; i< nsurf + 1; i++) {
        perturbation.eulerPhi[i] = body[i][0];
        perturbation.eulerPsi[i] = body[i][1];
        perturbation.eulerTheta[i] = body[i][2];
        perturbation.decenterX[i] = body[i][3];
        perturbation.decenterY[i] = body[i][4];
        perturbation.defocus[i] = body[i][5];
    }
    satbuffer += static_cast<long>(2.0*fabs(perturbation.defocus[nsurf]/(pixsize*1e-3)));
    perturbation.decenterX[nsurf] += centerx*1e-3;
    perturbation.decenterY[nsurf] += centery*1e-3;

    surface.centerx[nsurf] = 0.0;
    surface.centery[nsurf] = 0.0;
    surface.rmax[nsurf] = sqrt(pixelsx*pixelsx + pixelsy*pixelsy)*sqrt(2.0)/2.0*pixsize*1e-3;
    surface.height[nsurf] =  surface.height[nsurf - 1];
    for (long i = 0; i < totSurf; i++) {
        *(perturbation.rotationmatrix + 9*i + 0*3 + 0) = cos(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 0*3 + 1) = cos(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 0*3 + 2) = sin(perturbation.eulerPsi[i])*sin(perturbation.eulerTheta[i]);
        *(perturbation.rotationmatrix + 9*i + 1*3 + 0) = -sin(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 1*3 + 1) = -sin(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 1*3 + 2) = cos(perturbation.eulerPsi[i])*sin(perturbation.eulerTheta[i]);
        *(perturbation.rotationmatrix + 9*i + 2*3 + 0) = sin(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i]);
        *(perturbation.rotationmatrix + 9*i + 2*3 + 1) = -sin(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i]);
        *(perturbation.rotationmatrix + 9*i + 2*3 + 2) = cos(perturbation.eulerTheta[i]);
        *(perturbation.inverserotationmatrix + 9*i + 0*3 + 0) = cos(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 0*3 + 1) = -sin(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 0*3 + 2) = sin(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 1*3 + 0) = cos(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 1*3 + 1) = -sin(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 1*3 + 2) = -sin(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 2*3 + 0) = sin(perturbation.eulerTheta[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 2*3 + 1) = sin(perturbation.eulerTheta[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 2*3 + 2) = cos(perturbation.eulerTheta[i]);
    }
    for (long i = 0; i < totSurf; i++) {
        perturbation.rmax[i] = surface.rmax[surfaceLink[i]];
    }


    /* Update perturbation model from FEA file */
    size_t nNear = 4;
    // int degree=1;
    int leafSize = 10;
    for (long k = 0; k < nsurf; k++) {
        for (long m = 0; m < feaflag[k]; m++) {
            // std::cout<<"Loading "<<feafile[k*MAX_SURF*2 + m]<<std::endl;
            std::vector<double> ix(PERTURBATION_POINTS*PERTURBATION_POINTS), iy(PERTURBATION_POINTS*PERTURBATION_POINTS);
            for (long j = 0; j < PERTURBATION_POINTS; j++) {
                for (long l = 0; l < PERTURBATION_POINTS; l++) {
                    ix[l*PERTURBATION_POINTS + j] = perturbation.rmax[k]*perturbation.zernike_r_grid[j]*cos(perturbation.zernike_phi_grid[l]);
                    iy[l*PERTURBATION_POINTS + j] = perturbation.rmax[k]*perturbation.zernike_r_grid[j]*sin(perturbation.zernike_phi_grid[l]);
                }
            }
            long sizeFea;
            Fea feaTree(feafile[k*MAX_SURF*2 + m], leafSize, surface, perturbation, k, feascaling[k*MAX_SURF*2 + m], &sizeFea);
            nNear = round(0.5*sqrt(static_cast<double>(sizeFea)));
            if (nNear < 4) nNear = 4;
            feaTree.knnQueryFitDegree1(ix, iy, &perturbation.zernike_summed[k*PERTURBATION_POINTS*PERTURBATION_POINTS],
                    &perturbation.zernike_summed_nr_p[k*PERTURBATION_POINTS*PERTURBATION_POINTS],
                    &perturbation.zernike_summed_np_r[k*PERTURBATION_POINTS*PERTURBATION_POINTS], nNear);
            // for (long ii = 0; ii < PERTURBATION_POINTS; ii++) {
            //     for (long jj = 0; jj < PERTURBATION_POINTS; jj++) {
            //         if (perturbation.zernike_summed[k*PERTURBATION_POINTS*PERTURBATION_POINTS + ii*PERTURBATION_POINTS + jj] < 10e-3 &&
            //             perturbation.zernike_summed[k*PERTURBATION_POINTS*PERTURBATION_POINTS + ii*PERTURBATION_POINTS + jj] > -10e-3) {

            //         }
            //     }
            // }

        }
    }


    /// Large Angle Scattering
    if (laScatterProb > 0) {
        perturbation.miescatter = static_cast<double*>(calloc(10000, sizeof(double)));
        for (long j = 0; j < 10000; j++) {
            *(perturbation.miescatter + j) = 1.0/(1.0 + pow(((static_cast<float>(j))/10000.0*0.1*PI/180.0)/
                                                        (1.0*PI/180.0/3600.0), 3.5))*2.0*PI*
            (static_cast<float>(j));
        }
        double tempf1 = 0.0;
        for (long j = 0; j < 10000; j++) {
            tempf1 = tempf1 + *(perturbation.miescatter + j);
        }
        for (long j = 0; j < 10000; j++) {
            *(perturbation.miescatter + j) = *(perturbation.miescatter + j)/tempf1;
        }
        for (long j = 1; j < 10000; j++) {
            *(perturbation.miescatter + j) = *(perturbation.miescatter + j) + *(perturbation.miescatter + j - 1);
        }
    }

    // Tracking
    double angle = PI/2 - zenith;
    if (fabs(angle) < 0.5*PI/180.0) angle=0.5*PI/180.0;
    rotationrate = 15.04*cos(latitude)*cos(azimuth)/cos(angle);
    if (trackingfile == ".") trackingMode = 0;

    std::ifstream inStream(trackingfile.c_str());
    if (inStream) {
        readText trackingPars(trackingfile);
        size_t j = trackingPars.getSize();
        perturbation.jitterrot = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.jitterele = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.jitterazi = static_cast<double*>(calloc(j, sizeof(double)));
        screen.jitterwind = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.jittertime = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.windshake = static_cast<double*>(calloc(j, sizeof(double)));
        trackinglines = j;
        for (size_t t(0); t < j; t++) {
            std::istringstream iss(trackingPars[t]);
            iss >> perturbation.jittertime[t];
            double f1, f2, f3, f4, f5;
            iss >> f1 >> f2 >> f3 >> f4 >> f5;
            perturbation.jitterele[t] = f1*elevationjitter;
            perturbation.jitterazi[t] = f2*azimuthjitter;
            perturbation.jitterrot[t] = f3*rotationjitter;
            perturbation.windshake[t] = f5*windshake;
            // printf("%e\n",perturbation.windshake[t]);
            if (trackingMode == 0) {
                perturbation.jitterele[t] = 0.0;
                perturbation.jitterazi[t] = 0.0;
                perturbation.jitterrot[t] = 0.0;
            }
            if (windShakeMode == 0) {
                perturbation.windshake[t] = 0.0;
            }
            screen.jitterwind[t] = f4*windjitter;
        }
    }

    /// Setup chip material
    fprintf(stdout, "Electrifying Devices.\n");
    if (detectorMode != 0) {
        silicon.setup(devmaterial, sensorTempNominal + sensorTempDelta, nbulk, nf, nb, sf, sb, sensorthickness, overdepbias,
                  instrdir, chip.nampx, chip.nampy, pixsize, seedchip,
                  &impurityX, &impurityY, impurityvariation, minwavelength, maxwavelength);
    }

    if (contaminationmode == 1) {
        fprintf(stdout, "Contaminating Surfaces.\n");
        contamination.setup(surface.innerRadius, surface.outerRadius, nsurf - 1, PERTURBATION_POINTS,
                        pixsize, chip.nampx, chip.nampy, qevariation, seedchip);
    }

    fprintf(stdout, "Diffracting.\n");
   // 2nd kick
    // {

    //     Photon photon;
    //     Vector position;
    //     Vector angle;
    //     double shiftedAngle;
    //     double r,  phi;
    //     double rindex;
    //     long index;
    //     double radius;
    //     long atmtempdebug = 0;
    //     double bestvalue;
    //     double bestscale = 1.0;
    //     Vector largeAngle;
    //     double stdx = 0.0, stdy = 0.0;
    //     double radx = 0.0, rady = 0.0;
    //     double ncount = 0.0;
    //     screen.hffunc = static_cast<double*>(calloc(10000, sizeof(double)));
    //     screen.hffunc_n = static_cast<double*>(calloc(10000, sizeof(double)));

    //     atmtempdebug = atmdebug;

    //     for (long l = 0; l <= atmtempdebug; l++) {
    //         atmdebug = l;

    //         for (long zz = 0 ; zz < 100; zz++) {

    //         for (long i = 0; i < SCREEN_SIZE; i++) {
    //             for (long j = 0; j < SCREEN_SIZE; j++) {
    //                 screen.tfocalscreen[i*SCREEN_SIZE + j] = 0;
    //             }
    //         }

    //         photon.prtime = -1.0;

    //         for (long k = 0; k < 10; k++) {

    //             photon.wavelength = 0.3 + static_cast<double>(zz)/100.0*0.9;
    //             photon.wavelengthFactor = pow(photon.wavelength, -0.2)/screen.wavelengthfactor_nom;
    //             photon.time = random.uniform()*exptime;
    //             photon.absoluteTime = photon.time - exptime/2.0 + timeoffset;
    //             shiftedAngle = spiderangle + photon.time*rotationrate*ARCSEC;
    //             r = sqrt(random.uniform()*(maxr*maxr - minr*minr) + minr*minr);
    //             phi = random.uniform()*2*PI;
    //             position.x = r*cos(phi);
    //             position.y = r*sin(phi);
    //             index = find_linear(const_cast<double *>(&surface.radius[0]), SURFACE_POINTS, r, &rindex);
    //             position.z = interpolate_linear(const_cast<double *>(&surface.profile[0]), index, rindex);
    //             photon.xp = position.x;
    //             photon.yp = position.y;
    //             angle.x = 0.0;
    //             angle.y = 0.0;
    //             angle.z = -1.0;

    //             atmospherePropagate(&position, angle, -1, 2, &photon);
    //             for (long layer = 0; layer < natmospherefile; layer++) {
    //                 atmospherePropagate(&position, angle, layer, 2, &photon);
    //                 atmosphereIntercept(&position, layer, &photon);
    //                 atmosphereRefraction(&angle, layer, 2, &photon);
    //             }
    //             atmosphereDiffraction(&angle, &photon);
    //             for (long i = 0; i < SCREEN_SIZE; i++) {
    //                 for (long j = 0; j < SCREEN_SIZE; j++) {
    //                     screen.tfocalscreen[i*SCREEN_SIZE + j] += screen.focalscreen[i*SCREEN_SIZE + j];
    //                 }
    //             }


    //         }



    //         for (long i = 0; i < 10000; i++) {
    //             *(screen.hffunc + i) = 0.0;
    //         }
    //         for (long i = 0; i < 10000; i++) {
    //             *(screen.hffunc_n + i) = 0.0;
    //         }
    //         for (long i = 0; i < SCREEN_SIZE; i++) {
    //             for (long j = 0; j < SCREEN_SIZE; j++) {
    //                 radius = sqrt(pow((i - (SCREEN_SIZE/2 + 1)), 2) + pow((j - (SCREEN_SIZE/2 + 1)), 2));
    //                 *(screen.hffunc + ((long)(radius))) += screen.tfocalscreen[i*SCREEN_SIZE + j];
    //                 *(screen.hffunc_n + ((long)(radius))) += 1;
    //             }
    //         }

    //         for (long j = 0; j < 10000; j++) {
    //             if (screen.hffunc_n[j] < 1) {
    //                 *(screen.hffunc + j) = 0;
    //             } else {
    //                 *(screen.hffunc + j) = *(screen.hffunc + j)/(screen.hffunc_n[j])*((double)(j));
    //             }
    //         }
    //         double tempf1 = 0.0;
    //         for (long j = 0; j < 10000; j++) {
    //             tempf1 = tempf1 + *(screen.hffunc + j);
    //         }
    //         for (long j = 0; j < 10000; j++) {
    //             *(screen.hffunc + j) = *(screen.hffunc + j)/tempf1;
    //         }
    //         for (long j = 1; j < 10000; j++) {
    //             *(screen.hffunc + j) = *(screen.hffunc + j) + *(screen.hffunc + j - 1);
    //         }
    //         for (long j = 0; j < 10000; j++) {
    //             printf("xx %ld %ld %lf\n",zz,j,*(screen.hffunc + j));
    //         }

    //         }
    //     }
    // }


    // 2nd kick
    {

        Photon photon;
        Vector position;
        Vector angle;
        double shiftedAngle;
        double r,  phi;
        double rindex;
        long index;
        double radius;
        long atmtempdebug = 0;
        double bestvalue;
        double bestscale = 1.0;
        Vector largeAngle;
        double stdx = 0.0, stdy = 0.0;
        double radx = 0.0, rady = 0.0;
        double ncount = 0.0;
        screen.hffunc = static_cast<double*>(calloc(10000, sizeof(double)));
        screen.hffunc_n = static_cast<double*>(calloc(10000, sizeof(double)));
        photon.thread = 0;
        atmtempdebug = atmdebug;

        for (long l = 0; l <= atmtempdebug; l++) {
            atmdebug = l;

            for (long i = 0; i < SCREEN_SIZE; i++) {
                for (long j = 0; j < SCREEN_SIZE; j++) {
                    screen.tfocalscreen[i*SCREEN_SIZE + j] = 0;
                }
            }

            photon.prtime = -1.0;

            for (long k = 0; k < 10; k++) {

                photon.wavelength = 0.5;
                photon.wavelengthFactor = pow(photon.wavelength, -0.2)/screen.wavelengthfactor_nom;
                photon.time = random[0].uniform()*exptime;
                photon.absoluteTime = photon.time - exptime/2.0 + timeoffset;
                shiftedAngle = spiderangle + photon.time*rotationrate*ARCSEC;
                r = sqrt(random[0].uniform()*(maxr*maxr - minr*minr) + minr*minr);
                phi = random[0].uniform()*2*PI;
                position.x = r*cos(phi);
                position.y = r*sin(phi);
                index = find_linear(const_cast<double *>(&surface.radius[0]), SURFACE_POINTS, r, &rindex);
                position.z = interpolate_linear(const_cast<double *>(&surface.profile[0]), index, rindex);
                photon.xp = position.x;
                photon.yp = position.y;
                angle.x = 0.0;
                angle.y = 0.0;
                angle.z = -1.0;

                atmospherePropagate(&position, angle, -1, 2, &photon);
                for (long layer = 0; layer < natmospherefile; layer++) {
                    atmospherePropagate(&position, angle, layer, 2, &photon);
                    atmosphereIntercept(&position, layer, &photon);
                    atmosphereRefraction(&angle, layer, 3, &photon);
                }
                atmosphereDiffraction(&angle, &photon);
                for (long i = 0; i < SCREEN_SIZE; i++) {
                    for (long j = 0; j < SCREEN_SIZE; j++) {
                        screen.tfocalscreen[i*SCREEN_SIZE + j] += screen.focalscreen[i*SCREEN_SIZE + j];
                    }
                }


            }



            for (long i = 0; i < 10000; i++) {
                *(screen.hffunc + i) = 0.0;
            }
            for (long i = 0; i < 10000; i++) {
                *(screen.hffunc_n + i) = 0.0;
            }
            for (long i = 0; i < SCREEN_SIZE; i++) {
                for (long j = 0; j < SCREEN_SIZE; j++) {
                    radius = sqrt(pow((i - (SCREEN_SIZE/2 + 1)), 2) + pow((j - (SCREEN_SIZE/2 + 1)), 2));
                    *(screen.hffunc + ((long)(radius))) += screen.tfocalscreen[i*SCREEN_SIZE + j];
                    *(screen.hffunc_n + ((long)(radius))) += 1;
                }
            }

            for (long j = 0; j < 10000; j++) {
                if (screen.hffunc_n[j] < 1) {
                    *(screen.hffunc + j) = 0;
                } else {
                    *(screen.hffunc + j) = *(screen.hffunc + j)/(screen.hffunc_n[j])*((double)(j));
                }
            }
            double tempf1 = 0.0;
            for (long j = 0; j < 10000; j++) {
                tempf1 = tempf1 + *(screen.hffunc + j);
            }
            for (long j = 0; j < 10000; j++) {
                *(screen.hffunc + j) = *(screen.hffunc + j)/tempf1;
            }
            for (long j = 1; j < 10000; j++) {
                *(screen.hffunc + j) = *(screen.hffunc + j) + *(screen.hffunc + j - 1);
            }

            if (l == 0) {

                bestvalue = 1e10;
                bestscale = 1.0;
                for (double scales = 0.0; scales <= 2.01; scales += 0.01) {
                    stdx = 0.0;
                    stdy = 0.0;
                    ncount = 0.0;

                    for (long k = 0; k < 10000; k++) {

                        screen.secondKickSize = scales;
                        photon.wavelength = 0.5;
                        photon.wavelengthFactor = pow(photon.wavelength, -0.2)/screen.wavelengthfactor_nom;
                        photon.time = random[0].uniform()*1000.0;
                        photon.absoluteTime = photon.time - 1000.0/2.0 + timeoffset;
                        shiftedAngle = spiderangle + photon.time*rotationrate*ARCSEC;
                        r = sqrt(random[0].uniform()*(maxr*maxr - minr*minr) + minr*minr);
                        phi = random[0].uniform()*2*PI;
                        position.x = r*cos(phi);
                        position.y = r*sin(phi);
                        index = find_linear(const_cast<double *>(&surface.radius[0]), SURFACE_POINTS, r, &rindex);
                        position.z = interpolate_linear(const_cast<double *>(&surface.profile[0]), index, rindex);
                        photon.xp = position.x;
                        photon.yp = position.y;
                        angle.x = 0.0;
                        angle.y = 0.0;
                        angle.z = -1.0;
                        largeAngle.x = 0;
                        largeAngle.y = 0;

                        secondKick(&largeAngle, &photon);
                        atmospherePropagate(&position, angle, -1, 1, &photon);
                        for (long layer = 0; layer < natmospherefile; layer++) {
                            atmospherePropagate(&position, angle, layer, 1, &photon);
                            atmosphereIntercept(&position, layer, &photon);
                            atmosphereRefraction(&angle, layer, 1, &photon);
                        }
                        angle.x = angle.x + largeAngle.x;
                        angle.y = angle.y + largeAngle.y;

                        radx = angle.x/((totalseeing + 1e-6)/2.35482*ARCSEC*pow(1/cos(zenith), 0.6)*photon.wavelengthFactor);
                        rady = angle.y/((totalseeing + 1e-6)/2.35482*ARCSEC*pow(1/cos(zenith), 0.6)*photon.wavelengthFactor);
                        stdx += radx*radx*exp(-(radx*radx + rady*rady)/2.0/1.0);
                        stdy += rady*rady*exp(-(radx*radx + rady*rady)/2.0/1.0);
                        ncount += exp(-(radx*radx + rady*rady)/2.0/1.0);

                    }

                    stdx /= ncount;
                    stdy /= ncount;
                    screen.secondKickSize = sqrt(stdx + stdy);
                    if (fabs(screen.secondKickSize - 1.0) < bestvalue) {
                        bestvalue = fabs(screen.secondKickSize - 1.0);
                        bestscale = scales;
                    }

                }
                screen.secondKickSize = bestscale;
            }

        }

        atmdebug = atmtempdebug;


    }


    if (areaExposureOverride == 0) {
        nphot = static_cast<long long>(PI*(maxr/10.0*maxr/10.0 - minr/10.0*minr/10.0)*exptime*(totalnorm/H_CGS));
    } else {
        nphot = static_cast<long long>(1e6*(totalnorm/H_CGS));
    }
    if (opdfile == 1) {
        nphot = totalnorm;
    }

    settings();

    fprintf(stdout, "Number of Sources: %ld\n", nsource);
    fprintf(stdout, "Photons: %6.2e  Photons/cm^2/s: %6.2e\n", static_cast<double>(nphot), totalnorm/H_CGS);

    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Photon Raytrace\n");
    readText versionPars(bindir + "/version");
    fprintf(stdout, "%s\n", versionPars[0].c_str());
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");

    return(0);

}
