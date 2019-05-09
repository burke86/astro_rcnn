///
/// @package phosim
/// @file airglow.cpp
/// @brief airglow generator
///
/// @brief Created by:
/// @author En-Hsin Peng (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

void Atmosphere::airglow(long seed, const std::string & name, long screenSize) {

    long nx = screenSize;
    long ny = screenSize;
    long n = nx*ny;

    std::vector<double> in1(n, 0);
    std::vector<double> in(n, 0);

    static int ns(3);
    std::vector<double> kmax(ns, 0);
    std::vector<double> kmin(ns, 0);
    std::vector<double> pixs(ns, 0);
    pixs[0] = 15.0/3600;  //in degrees
    kmin[0] = 4.0/screenSize/pixs[0];
    kmax[0] = 1.0/pixs[0];
    for (int ii = 1; ii < ns; ii++) {
        pixs[ii] = pixs[ii - 1]*16;
        kmax[ii] = kmin[ii - 1];
        kmin[ii] = kmin[ii - 1]/16;
    }

    fftw_complex *out = NULL;
    fftMalloc(&out, nx*ny);
    fftw_complex *incom = NULL;
    fftMalloc(&incom, nx*ny);

    /* initialize random seed */
    random.setSeed32(seed);
    random.unwind(10000);
    double idxa = 1.3, idxb = 3.0;
    double n0a = 6.6e-3;
    double n0b = 1e-5;

    for (int ii = 0; ii < ns; ii++) {
        std::vector<double> xgrid(screenSize, 0);
        for (int ix = 0; ix < nx; ix++) {
            xgrid[ix] = ix*pixs[ii];
            double kx = (ix < (nx/2)) ? static_cast<double>(ix): static_cast<double>(nx - ix);
            for (int iy = 0; iy < ny; iy++) {
                double ky = (iy < (ny/2)) ? static_cast<double>(iy): static_cast<double>(ny - iy);
                long i = iy + ny*ix;
                double dkk = sqrt(kx*kx + ky*ky)/screenSize/pixs[ii];
                out[i][0] = 0.0;
                out[i][1] = 0.0;
                if (dkk >= kmax[ii] || dkk < kmin[ii]) {
                    out[i][0] = 0.0;
                } else {
                    double dkkp = (n0a*pow(dkk, -idxa) + n0b*pow(dkk, -idxb))/pixs[ii];
                    out[i][0] = dkkp*random.normal();
                }
            }
        }
        out[0][0] = 0.0;
        inverseFFT(nx, ny, out, incom);

        for (long i=0; i < n; i++) {
            in1[i] = incom[i][1]/screenSize;
        }

        double xpos = random.uniform()*(nx*pixs[ii] - nx*pixs[0]);
        double ypos = random.uniform()*(ny*pixs[ii] - ny*pixs[0]);
        for (int ix = 0; ix < nx; ix++) {
            for (int iy=0; iy < ny; iy++) {
                long i = iy + ny*ix;
                double dx, dy;
                long idx0 = find_linear(const_cast<double *>(&xgrid[0]), screenSize,
                                        xpos + ix*pixs[0], &dx);
                long idy0 = find_linear(const_cast<double *>(&xgrid[0]), screenSize,
                                        ypos + iy*pixs[0], &dy);
                in[i] += interpolate_bilinear(&in1[0], screenSize, idx0, dx, idy0, dy);
            }
        }
    }

    double avg = 0.0;
    double rms = 0.0;
    double imax = 0.0;
    double imin = 0.0;
    for (long i = 0; i < n; i++) {
        avg += in[i];
        imax = (in[i] > imax) ? in[i] : imax;
        imin = (in[i] < imin) ? in[i] : imin;
    }
    avg /= n;
    for (long i = 0; i < n; i++) {
        rms += (in[i] - avg)*(in[i] - avg);
    }
    rms /= n;
    rms = sqrt(rms);

    for (long i = 0; i < n; i++) {
        in[i] -= avg;
    }

    fitsfile *fptr;
    fitsCreateImage(&fptr, "!" + name + ".fits.gz");
    fitsWriteImage(fptr, nx, ny, const_cast<double*>(&in[0]));

}
