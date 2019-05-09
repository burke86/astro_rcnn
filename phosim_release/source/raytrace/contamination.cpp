///
/// @package phosim
/// @file contamination.cpp
/// @brief contamination class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <math.h>
#include <fitsio.h>
#include <fitsio2.h>

#include "contamination.h"
#include "constants.h"
#include "ancillary/random.h"

void Contamination::setup (double *innerRadius, double *outerRadius, long totalSurface, long points,
                           double pixsize, long nx, long ny, float qevariation, long seedchip) {

    double mmd = 0.3;
    double st = 0.0025;
    double stw = 0.0025;
    double sd, sdw, sdd;
    double sigma = 1.0;
    double mmdw = 100.0;
    double sigmaw = 1.0;
    double sizetol = 0.1;
    double gridsize;
    double smallgridsize;
    double area;
    long long particles, particlesd;
    double size;
    double fraction;
    long r0, phi0, x0, y0;
    double deltar, deltaphi, deltax, deltay;
    long rd, pd, xd, yd;
    long rc, pc, xc, yc;
    long dx, dy, dr, dphi;
    double x, y, r, phi;
    long number;
    double oversim;
    double totalloss;

    random.setSeed32Fixed(1000);
    oversim = 100.0;
    sd = st/(PI*exp(2.0*log(mmd) + 0.5*2.0*2.0*sigma*sigma)*1e-3*1e-3);
    sdd = st/(PI*exp(2.0*log(mmd) + 0.5*2.0*2.0*sigma*sigma)*1e-3*1e-3);
    sdw = stw/(PI*exp(2.0*log(mmdw) + 0.5*2.0*2.0*sigmaw*sigmaw)*1e-3*1e-3);

    //optics
    surfacenumber = new long[totalSurface]();
    for (long i = 0; i < totalSurface; i++) {
        surfacenumber[i] = 0;
    }

    transmission = new float[totalSurface*points*points]();
    for (long i = 0; i < totalSurface; i++) {
        for (long j = 0; j < points; j++) {
            for (long k = 0; k < points; k++) {
                transmission[i*points*points + j*points + k] = 1.0;
            }
        }
    }

    chiptransmission = new float[nx*ny]();
    for (long j = 0; j < nx; j++) {
        for (long k = 0; k < ny; k++) {
            chiptransmission[k*nx + j] = 1.0 - qevariation*fabs(random.normalFixed());
        }
    }

    surfacelistmap = new long[totalSurface*points*points]();
    for (long i = 0; i < totalSurface; i++) {
        for (long j = 0; j < points; j++) {
            for (long k = 0; k < points; k++) {
                surfacelistmap[i*points*points + j*points + k] = -1;
            }
        }
    }

    surfacelistx.resize(totalSurface);
    surfacelisty.resize(totalSurface);
    surfacelists.resize(totalSurface);
    surfacelistt.resize(totalSurface);

    for (long i = 0; i<totalSurface; i++) {
        deltar = outerRadius[i] - sqrt((outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i])*(points - 2)/(points - 1)+
                                       innerRadius[i]*innerRadius[i]);
        deltaphi = 2*PI/(points - 1);
        totalloss = 0.0;
        area = PI*(outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i]);
        particles = static_cast<long>(ceil(sd*area + sdw*area));
        particlesd = static_cast<long>(ceil(sd*area));
        gridsize = sizetol*sqrt(area/PI)/(static_cast<double>(points));
        smallgridsize = sizetol*gridsize;
        long long j = 0;
        long long c1 = 0;
        long long c2 = 0;
        long long c3 = 0;
        while (j < particles) {
            if (j < particlesd) {
                size = exp(log(mmd) + sigma*random.normalFixed())*1e-3;
            } else {
                size = exp(log(mmdw) + sigmaw*random.normalFixed())*1e-3;
            }
            number = oversim;
            if (number < 1) number = 1;
            if (size > gridsize) {
                // large contaminant
                for (long k = 0; k < number; k++) {
                    c1++;
                    r = sqrt(random.uniformFixed()*(outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i]) +
                             innerRadius[i]*innerRadius[i]);
                    phi = random.uniformFixed()*2*PI;
                    r0 = static_cast<long>((r*r - innerRadius[i]*innerRadius[i])/(outerRadius[i]*outerRadius[i]
                                                                                  - innerRadius[i]*innerRadius[i])*(static_cast<double>(points)-1));
                    phi0 = static_cast<long>((phi)/(2*PI)*(static_cast<double>(points) - 1));
                    x = r*cos(phi);
                    y = r*sin(phi);
                    surfacelistx[i].push_back(x);
                    surfacelisty[i].push_back(y);
                    surfacelists[i].push_back(size);
                    if (j < particlesd) {
                        surfacelistt[i].push_back(0);
                    } else {
                        surfacelistt[i].push_back(1);
                    }
                    dr = ceil((size/deltar) + 1);
                    dphi = ceil((size/(deltaphi*sqrt((outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i])*
                                                     (static_cast<double>(r0))/(static_cast<double>(points) - 1) +
                                                     innerRadius[i]*innerRadius[i]))) + 1);
                    for (rd = r0 - dr; rd <= r0 + dr; rd++) {
                        for (pd = phi0 - dphi; pd <= phi0 + dphi; pd++) {
                            rc = rd;
                            pc = pd;
                            if (rc < 0) rc = 0;
                            if (rc > points - 1) rc = points - 1;
                            if (pc < 0) pc = pc + points;
                            if (pc > points - 1) pc = pc - points;
                            if (surfacelistmap[i*points*points + rc*points + pc] != -1) {
                                if (surfacelists[i][surfacenumber[i]] > surfacelists[i][surfacelistmap[i*points*points + rc*points + pc]]) {
                                    surfacelistmap[i*points*points + rc*points + pc] = surfacenumber[i];
                                }
                            } else {
                                    surfacelistmap[i*points*points + rc*points + pc] = surfacenumber[i];
                            }
                        }
                    }
                    surfacenumber[i]++;
                }
            } else if (size > smallgridsize) {
                // small contaminant
                for (long k = 0; k < number; k++) {
                    c2++;
                    r = sqrt(random.uniformFixed()*(outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i]) +
                             innerRadius[i]*innerRadius[i]);
                    phi = random.uniformFixed()*2*PI;
                    r0 = static_cast<long>((r*r - innerRadius[i]*innerRadius[i])/(outerRadius[i]*outerRadius[i]
                                                                                  - innerRadius[i]*innerRadius[i])*(static_cast<double>(points)-1));
                    phi0 = static_cast<long>((phi)/(2*PI)*(static_cast<double>(points) - 1));
                    if (r0 < 0) r0 = 0;
                    if (r0 > points - 1) r0 = points - 1;
                    if (phi0 < 0) phi0 = phi0 + points;
                    if (phi0 > points - 1) phi0 = phi0 - points;
                    fraction = PI*size*size/(area/(static_cast<double>(points))/(static_cast<double>(points)));
                    transmission[i*points*points + r0*points + phi0] -= fraction;
                    if (transmission[i*points*points + r0*points + phi0] < 0) {
                        transmission[i*points*points + r0*points + phi0] = 0.0;
                    }
                }
            } else {
                // super small contaminant
                totalloss += number*PI*size*size;
                c3 += number;
            }

            j += number;
        }
        double avg = 0.0;
        for (long k = 0; k < points; k++) {
            for (long l = 0; l < points; l++) {
                transmission[i*points*points + k*points + l] -= totalloss/area;
                if (transmission[i*points*points + k*points + l] < 0) {
                    transmission[i*points*points + k*points + l] = 0.0;
                }
                avg += transmission[i*points*points + k*points + l];
            }
        }

    }

    // int status;
    // long naxesa[2];
    // fitsfile *faptr;
    // status=0;
    // fits_create_file(&faptr,"!cont_m.fits",&status);
    // naxesa[0]=points; naxesa[1]=points;
    // fits_create_img(faptr,FLOAT_IMG,2,naxesa,&status);
    // fits_write_img(faptr,TFLOAT,1,points*points,transmission,&status);
    // fits_close_file(faptr,&status);

    // status=0;
    // fits_create_file(&faptr,"!cont_p.fits",&status);
    // naxesa[0]=points; naxesa[1]=points;
    // fits_create_img(faptr,LONG_IMG,2,naxesa,&status);
    // fits_write_img(faptr,TLONG,1,points*points,surfacelistmap,&status);
    // fits_close_file(faptr,&status);

    //chip
    random.setSeed32Fixed(1000 + seedchip);
    chiplistmap = new long[nx*ny]();
    for (long i = 0; i < nx; i++) {
        for (long j = 0; j < ny; j++) {
            chiplistmap[j*nx + i] = -1;
        }
    }
    chipnumber = 0;
    area = (nx*pixsize*ny*pixsize*1e-3*1e-3);
    particles = static_cast<long>(ceil(sdd*area));
    gridsize = sizetol*pixsize*1e-3/sqrt(PI);
    smallgridsize = sizetol*gridsize;
    deltax = pixsize*1e-3;
    deltay = pixsize*1e-3;
    for (long j = 0; j < particles; j++) {
        size = exp(log(mmd) + sigma*random.normalFixed())*1e-3;
        x = random.uniformFixed()*nx*pixsize*1e-3;
        y = random.uniformFixed()*ny*pixsize*1e-3;
        x0 = round((long)(x/pixsize/1e-3));
        y0 = round((long)(y/pixsize/1e-3));

        if (size > gridsize) {
            chiplistx.push_back(x);
            chiplisty.push_back(y);
            chiplists.push_back(size);
            dx = ceil((size/deltax) + 2);
            dy = ceil((size/deltay) + 2);
            for (xd = x0 - dx; xd <= x0 + dx; xd++) {
                for (yd = y0 - dy; yd <= y0 + dy; yd++) {
                    xc = xd;
                    yc = yd;
                    if (xc < 0) xc = 0;
                    if (xc >= nx) xc = nx - 1;
                    if (yc < 0) yc = 0;
                    if (yc >= ny) yc = ny - 1;
                    if (chiplistmap[nx*yc + xc] != -1) {
                        if (chiplists[chipnumber] > chiplists[chiplistmap[nx*yc + xc]]) {
                                chiplistmap[nx*yc + xc] = chipnumber;
                            }
                    } else {
                        chiplistmap[nx*yc + xc] = chipnumber;
                    }
                }
            }
            chipnumber++;
        } else {
            fraction = PI*size*size/(pixsize*pixsize*1e-3*1e-3);
            if (x0 < 0) x0 = 0;
            if (y0 < 0) y0 = 0;
            if (x0 > nx - 1) x0 = nx - 1;
            if (y0 > ny - 1) y0 = ny - 1;
            chiptransmission[y0*nx + x0] -= fraction;
            if (chiptransmission[y0*nx + x0] < 0) {
                chiptransmission[y0*nx + x0] = 0.0;
            }
        }

    }

    elements = 100000;
    henyey_greenstein = new double[elements]();
    henyey_greenstein_mu = new double[elements]();
    double g = 0.6;
    for (long i = 0; i < elements; i++) {
        henyey_greenstein_mu[i] = -1 + 2.0*(static_cast<double>(i))/(elements - 1);
        henyey_greenstein[i] = (1 - g*g)/pow(1 + g*g - 2*g*henyey_greenstein_mu[i], 1.5);
    }
    double total = 0.0;
    for (long i = 0; i < elements; i++) {
        total += henyey_greenstein[i];
    }
    for (long i = 0; i < elements; i++) {
        henyey_greenstein[i] = henyey_greenstein[i]/total;
    }
    for (long i = 1; i < elements; i++) {
        henyey_greenstein[i] = henyey_greenstein[i] + henyey_greenstein[i - 1];
    }
    absorptionLength = 10.0;

    henyey_greenstein_w = new double[elements]();
    henyey_greenstein_mu_w = new double[elements]();
    g = 0.9;
    for (long i = 0; i < elements; i++) {
        henyey_greenstein_mu_w[i] = -1 + 2.0*(static_cast<double>(i))/(elements - 1);
        henyey_greenstein_w[i] = (1 - g*g)/pow(1 + g*g - 2*g*henyey_greenstein_mu_w[i], 1.5);
    }
    total = 0.0;
    for (long i = 0; i < elements; i++) {
        total += henyey_greenstein_w[i];
    }
    for (long i = 0; i < elements; i++) {
        henyey_greenstein_w[i] = henyey_greenstein_w[i]/total;
    }
    for (long i = 1; i < elements; i++) {
        henyey_greenstein_w[i] = henyey_greenstein_w[i] + henyey_greenstein_w[i - 1];
    }

    // status = 0;
    // fits_create_file(&faptr, "!cont.fits", &status);
    // naxesa[0] = nx; naxesa[1] = ny;
    // fits_create_img(faptr, LONG_IMG, 2, naxesa, &status);
    // fits_write_img(faptr, TLONG, 1, nx*ny, chiplistmap, &status);
    // fits_close_file(faptr, &status);

    // status = 0;
    // fits_create_file(&faptr, "!cont_t.fits", &status);
    // naxesa[0] = nx; naxesa[1] = ny;
    // fits_create_img(faptr, FLOAT_IMG, 2, naxesa, &status);
    // fits_write_img(faptr, TFLOAT, 1, nx*ny, chiptransmission, &status);
    // fits_close_file(faptr, &status);

}
