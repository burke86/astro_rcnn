///
/// @package phosim
/// @file turb2d.cpp
/// @brief turbulence generator
///
/// @brief Created by:
/// @author J. Garrett Jernigan (UCB)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
/// @author J. Chiang (SLAC)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

void Atmosphere::turb2d(long seed, double see5, double outerx, double outers,
                        double zenith, double wavelength,
                        const std::string & name, long N_size) {

    double lambda;
    lambda = 2.0 + (5.0/3.0);

    double r0;
    double wavelengthfactor;
    double wavelengthfactor_nom;
    double newsee;

    wavelengthfactor_nom = pow(0.5, -0.2);
    wavelengthfactor = pow(wavelength, -0.2)/wavelengthfactor_nom;

    newsee = see5;
    if (newsee < 1e-1) newsee = 1e-1;

    r0 = (0.98 * (wavelength*1e4*1.0e-8)/ ( (3.1425926/(180.0*60.0*60.0))
                                            * (newsee)*pow(1/cos(zenith*M_PI/180.),0.6)*wavelengthfactor) );

    long nx = N_size;
    long ny = N_size;
    int N = nx*ny;
    int ixp;
    int iyp;
    int ixm;
    int iym;

    double norm = 0.0;
    double finenorm = 0.0;
    double dennorm = 0.0;

    double xrms, yrms, xavg, yavg, xyrms;
    double rms, avg;

    fftw_complex *out_high = NULL;
    fftMalloc(&out_high, nx*ny);
    fftw_complex *out_low = NULL;
    fftMalloc(&out_low, nx*ny);
    fftw_complex *out_full = NULL;
    fftMalloc(&out_full, nx*ny);

    fftw_complex *inc_high = NULL;
    fftMalloc(&inc_high, nx*ny);
    fftw_complex *inc_low = NULL;
    fftMalloc(&inc_low, nx*ny);
    fftw_complex *inc_full = NULL;
    fftMalloc(&inc_full, nx*ny);


    double *in_high = NULL;
    fftMalloc(&in_high, nx*ny);
    double *in_low = NULL;
    fftMalloc(&in_low, nx*ny);
    double *in_full = NULL;
    fftMalloc(&in_full, nx*ny);


    std::vector<float> derx(nx*ny, 0);
    std::vector<float> dery(nx*ny, 0);

    std::vector<float> largex(nx*ny, 0);
    std::vector<float> largey(nx*ny, 0);
    std::vector<float> coarsex(nx*ny, 0);
    std::vector<float> coarsey(nx*ny, 0);
    std::vector<float> mediumx(nx*ny, 0);
    std::vector<float> mediumy(nx*ny, 0);
    std::vector<float> finex(nx*ny, 0);
    std::vector<float> finey(nx*ny, 0);
    std::vector<float> largep(nx*ny, 0);
    std::vector<float> largeh(nx*ny, 0);
    std::vector<float> coarsep(nx*ny, 0);
    std::vector<float> coarseh(nx*ny, 0);
    std::vector<float> mediump(nx*ny, 0);
    std::vector<float> mediumh(nx*ny, 0);
    std::vector<float> finep(nx*ny, 0);
    std::vector<float> fineh(nx*ny, 0);
    // std::vector<float> largel(nx*ny, 0);
    // std::vector<float> coarsel(nx*ny, 0);
    // std::vector<float> mediuml(nx*ny, 0);
    // std::vector<float> finel(nx*ny, 0);

    float w;
    float rcrit;
    float powerbuffer;
    int derstep;

    derstep = 1;
    powerbuffer = 8.0*16.0/8.0;
    rcrit = 1.5;
    w = 4.0; // Birkett

    // set the seed for random numbers
    random.setSeed32(seed);
    random.unwind(10000);

    double kx;
    double ky;
    double dkk;
    double dkksqr;
    double dkkp;
    double dklow;
    double scale = 0.0, scale_next = 0.0;

    for (int k = 0; k < 4; k++) {

        if (k == 0) {
            scale = 512.0;
            scale_next = 1000000.0;
        }
        if (k == 1) {
            scale = 64.0;
            scale_next = 512.0;
        }
        if (k == 2) {
            scale = 8.0;
            scale_next = 64.0;
        }
        if (k == 3) {
            scale = 1.0;
            scale_next = 8.0;
        }

        dkksqr = 0.5*sqrt((double)(nx*nx) + (double)(ny*ny));
        dkksqr /= scale;
        dklow = (1024.0/outers);

        for (int ix = 0; ix < nx; ix++) {
            if (ix < (nx/2)) {
                kx = static_cast<double>(ix);
            } else {
                kx = static_cast<double>(nx - ix);
            }
            for (int iy=0; iy < ny; iy++) {
                int i = iy + ny*ix;
                if (iy < (ny/2)) {
                    ky = static_cast<double>(iy);
                } else {
                    ky = static_cast<double>(ny - iy);
                }

                dkk = sqrt(kx*kx + ky*ky);
                dkk /= scale;

                dkkp = sqrt(1.0/pow((dkk*dkk + dklow*dklow), 0.5*lambda))/scale;
                if (i==0) dkkp = 0.0;

                // screen power division
                if (dkk > (1024.0/scale)/sqrt(2.0)/sqrt(powerbuffer) ||
                    dkk <= (1024.0/scale_next)/sqrt(2.0)/sqrt(powerbuffer)) {
                    dkkp = 0.0;
                }

                if (dkk > dkksqr) dkkp = 0.0;

                out_full[i][0] = dkkp*random.normal();
                out_full[i][1] = dkkp*random.normal();
                out_low[i][0] = out_full[i][0];
                out_low[i][1] = out_full[i][1];
                out_high[i][0] = out_full[i][0];
                out_high[i][1] = out_full[i][1];


                // second kick division
                if (dkk > (1024.0/(r0*rcrit)) ) {
                    out_low[i][0] *= 0.0;
                    out_low[i][1] *= 0.0;
                }
                if (dkk <= (1024.0/(r0*rcrit)) ) {
                    out_high[i][0] *= 0.0;
                    out_high[i][1] *= 0.0;
                }

                // outer outer scale
                if ( dkk < (1024.0/outerx) ) {
                    out_full[i][0] *= 0.0;
                    out_full[i][1] *= 0.0;
                    out_low[i][0] *= 0.0;
                    out_low[i][1] *= 0.0;
                    out_high[i][0] *= 0.0;
                    out_high[i][1] *= 0.0;
                }
            }
        }

        // inverse fft
        inverseFFT(nx, ny, out_full, inc_full);
        inverseFFT(nx, ny, out_high, inc_high);
        inverseFFT(nx, ny, out_low, inc_low);


        for (int i=0; i < N; i++) {
            in_full[i] = inc_full[i][0]/static_cast<double>(N);
            in_high[i] = inc_high[i][0]/static_cast<double>(N);
            in_low[i] = inc_low[i][0]/static_cast<double>(N);
        }

        // compute derivatives
        for (int ix=0; ix < nx; ix++) {
            for (int iy=0; iy < ny; iy++) {

                ixp = ix + derstep;
                if (ixp > nx - 1) ixp -= nx;
                iyp = iy + derstep;
                if (iyp > ny - 1) iyp -= ny;
                ixm = ix - derstep;
                if (ixm < 0) ixm += nx;
                iym = iy - derstep;
                if (iym < 0) iym += ny;

                derx[iy + ny*ix] = ((in_low[iy  + ny*ixp] - in_low[iy  + ny*ixm])*w+
                                    (in_low[iyp + ny*ixp] - in_low[iyp + ny*ixm])+
                                    (in_low[iym + ny*ixp] - in_low[iym + ny*ixm]))
                    /2.0/(static_cast<float>(derstep))/(scale)/(2.0 + w);

                dery[iy + ny*ix] = ((in_low[iyp + ny*ix ] - in_low[iym + ny*ix ])*w+
                                    (in_low[iyp + ny*ixp] - in_low[iym + ny*ixp])+
                                    (in_low[iyp + ny*ixm] - in_low[iym + ny*ixm]))
                    /2.0/(static_cast<float>(derstep))/(scale)/(2.0 + w);

            }
        }

        // compute rms of derivatives
        xrms = 0.0;
        yrms = 0.0;
        xavg = 0.0;
        yavg = 0.0;
        xyrms = 0.0;
        for (int ix=0; ix < nx; ix++) {
            for (int iy=0; iy < ny; iy++) {
                xrms += derx[iy + ny*ix]*derx[iy + ny*ix];
                yrms += dery[iy + ny*ix]*dery[iy + ny*ix];
                xavg += derx[iy + ny*ix];
                yavg += dery[iy + ny*ix];
            }
        }
        xavg = xavg/((double)(nx*ny));
        yavg = yavg/((double)(nx*ny));
        xrms = sqrt(xrms/((double)(nx*ny)) - xavg*xavg);
        yrms = sqrt(yrms/((double)(nx*ny)) - yavg*yavg);
        xyrms = sqrt((xrms*xrms + yrms*yrms)/2.0);
        norm += xyrms*xyrms;
        if (k == 3 || k == 2) finenorm += xyrms*xyrms;

        // compute rms of phase
        rms = 0.0;
        avg = 0.0;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                rms += in_full[iy + ny*ix]*in_full[iy + ny*ix];
                avg += in_full[iy + ny*ix];
            }
        }
        avg = avg/((double)(nx*ny));
        rms = sqrt(rms/((double)(nx*ny)) - avg*avg);
        dennorm += rms*rms;

        // copy to arrays
        double totalx = 0.0;
        double totaly = 0.0;
        double totalp = 0.0;
        double totalh = 0.0;
        double totall = 0.0;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                if (k == 0) {
                    largex[iy + ny*ix] = derx[iy + ny*ix];
                    largey[iy + ny*ix] = dery[iy + ny*ix];
                    largep[iy + ny*ix] = in_full[iy + ny*ix];
                    largeh[iy + ny*ix] = in_high[iy + ny*ix];
                    // largel[iy + ny*ix] = in_low[iy + ny*ix];
                }
                if (k == 1) {
                    coarsex[iy + ny*ix] = derx[iy + ny*ix];
                    coarsey[iy + ny*ix] = dery[iy + ny*ix];
                    coarsep[iy + ny*ix] = in_full[iy + ny*ix];
                    coarseh[iy + ny*ix] = in_high[iy + ny*ix];
                    // coarsel[iy + ny*ix] = in_low[iy + ny*ix];
                }
                if (k == 2) {
                    mediumx[iy + ny*ix] = derx[iy + ny*ix];
                    mediumy[iy + ny*ix] = dery[iy + ny*ix];
                    mediump[iy + ny*ix] = in_full[iy + ny*ix];
                    mediumh[iy + ny*ix] = in_high[iy + ny*ix];
                    // mediuml[iy + ny*ix]=in_low[iy + ny*ix];
                }
                if (k == 3) {
                    finex[iy + ny*ix] = derx[iy + ny*ix];
                    finey[iy + ny*ix] = dery[iy + ny*ix];
                    finep[iy + ny*ix] = in_full[iy + ny*ix];
                    fineh[iy + ny*ix] = in_high[iy + ny*ix];
                    // finel[iy + ny*ix]=in_low[iy + ny*ix];
                }

                totalx += derx[iy + ny*ix];
                totaly += derx[iy + ny*ix];
                totalp += in_full[iy + ny*ix];
                totalh += in_high[iy + ny*ix];
                totall += in_low[iy + ny*ix];

            }
        }

    }

    // co-normalize
    norm = sqrt(norm);
    for (int ix=0; ix < nx; ix++) {
        for (int iy=0; iy < ny; iy++) {
            largex[iy + ny*ix] /= norm;
            largey[iy + ny*ix] /= norm;
            coarsex[iy + ny*ix] /= norm;
            coarsey[iy + ny*ix] /= norm;
            mediumx[iy + ny*ix] /= norm;
            mediumy[iy + ny*ix] /= norm;
            finex[iy + ny*ix]   /= norm;
            finey[iy + ny*ix]   /= norm;
        }
    }

    // output fits files
    fitsfile *fptr = NULL;

    fitsCreateImage(&fptr, "!" + name + "_largex.fits.gz");
    fitsWriteKey(fptr, "NORM", norm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&largex[0]));

    fitsCreateImage(&fptr, "!" + name + "_largey.fits.gz");
    fitsWriteKey(fptr, "NORM", norm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&largey[0]));

    fitsCreateImage(&fptr, "!" + name + "_coarsex.fits.gz");
    fitsWriteKey(fptr, "NORM", norm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&coarsex[0]));

    fitsCreateImage(&fptr, "!" + name + "_coarsey.fits.gz");
    fitsWriteKey(fptr, "NORM", norm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&coarsey[0]));

    fitsCreateImage(&fptr, "!" + name + "_mediumx.fits.gz");
    fitsWriteKey(fptr, "NORM", norm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&mediumx[0]));

    fitsCreateImage(&fptr, "!" + name + "_mediumy.fits.gz");
    fitsWriteKey(fptr, "NORM", norm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&mediumy[0]));

    fitsCreateImage(&fptr, "!" + name + "_finex.fits.gz");
    fitsWriteKey(fptr, "NORM", norm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&finex[0]));

    fitsCreateImage(&fptr, "!" + name + "_finey.fits.gz");
    fitsWriteKey(fptr, "NORM", norm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&finey[0]));

    dennorm = sqrt(dennorm);
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            largep[iy + ny*ix] /= dennorm;
            largeh[iy + ny*ix] /= dennorm;
            coarsep[iy + ny*ix] /= dennorm;
            coarseh[iy + ny*ix] /= dennorm;
            mediump[iy + ny*ix] /= dennorm;
            mediumh[iy + ny*ix] /= dennorm;
            finep[iy + ny*ix]   /= dennorm;
            fineh[iy + ny*ix]   /= dennorm;
        }
    }

    fitsCreateImage(&fptr, "!" + name + "_largep.fits.gz");
    fitsWriteKey(fptr, "NORM", dennorm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&largep[0]));

    fitsCreateImage(&fptr, "!" + name + "_coarsep.fits.gz");
    fitsWriteKey(fptr, "NORM", dennorm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&coarsep[0]));

    fitsCreateImage(&fptr, "!" + name + "_mediump.fits.gz");
    fitsWriteKey(fptr, "NORM", dennorm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&mediump[0]));

    fitsCreateImage(&fptr, "!" + name + "_mediumh.fits.gz");
    fitsWriteKey(fptr, "NORM", dennorm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&mediumh[0]));

    fitsCreateImage(&fptr, "!" + name + "_finep.fits.gz");
    fitsWriteKey(fptr, "NORM", dennorm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&finep[0]));

    fitsCreateImage(&fptr, "!" + name + "_fineh.fits.gz");
    fitsWriteKey(fptr, "NORM", dennorm, "");
    fitsWriteImage(fptr, nx, ny, const_cast<float*>(&fineh[0]));

    // fitsCreateImage(&fptr, "!" + name + "_largeh.fits.gz");
    // fitsWriteKey(fptr, "NORM", dennorm, "");
    // fitsWriteImage(fptr, nx, ny, const_cast<float*>(&largeh[0]));

    // fitsCreateImage(&fptr, "!" + name + "_coarseh.fits.gz");
    // fitsWriteKey(fptr, "NORM", dennorm, "");
    // fitsWriteImage(fptr, nx, ny, const_cast<float*>(&coarseh[0]));

    // fitsCreateImage(&fptr, "!" + name + "largel.fits.gz");
    // fitsWriteKey(fptr, "NORM", dennorm, "");
    // fitsWriteImage(fptr, nx, ny, const_cast<float*>(&largel[0]));

    // fitsCreateImage(&fptr, "!" + name + "_coarsel.fits.gz");
    // fitsWriteKey(fptr, "NORM", dennorm, "");
    // fitsWriteImage(fptr, nx, ny, const_cast<float*>(&coarsel[0]));

    // fitsCreateImage(&fptr, "!" + name + "_mediuml.fits.gz");
    // fitsWriteKey(fptr, "NORM", dennorm, "");
    // fitsWriteImage(fptr, nx, ny, const_cast<float*>(&mediuml[0]));

    // fitsCreateImage(&fptr, "!" + name + "_finel.fits.gz");
    // fitsWriteKey(fptr, "NORM", dennorm, "");
    // fitsWriteImage(fptr, nx, ny, const_cast<float*>(&finel[0]));

}
