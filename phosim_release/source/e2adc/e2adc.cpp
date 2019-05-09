///
/// @package phosim
/// @file e2adc.cpp
/// @brief electron digitization code
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author En-Hsin Peng (Purdue)
/// @author Colin Burke (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>

#include "e2adc.h"
#include "ancillary/readtext.h"

void E2adc::setup() {

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Electron to ADC Image Converter" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;

    chipid = "R22_S11";
    instrdir = "../data/lsst";
    filter = 0;
    flatdir = 0;
    tarfile = 0;
    readorder = 1;
    serialcte = 0.999995;
    parallelcte = 1.0;
    nonlinear = 0.0;
    vistime = 33.0;
    nsnap = 2;
    ngroups = 1;
    nframes = 1;
    nskip = 0;
    wellDepth = 100000;
    minx = 0;
    miny = 0;
    sequenceMode = -1;
    
    // Read parameters from stdin.
    readText pars(std::cin);

    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        readText::get(line, "instrdir", instrdir);
        readText::get(line, "flatdir", flatdir);
        readText::get(line, "tarfile", tarfile);
        readText::get(line, "filter", filter);
        readText::get(line, "nonlinear", nonlinear);
        readText::get(line, "welldepth", wellDepth);
        readText::get(line, "parallelcte", parallelcte);
        readText::get(line, "serialcte", serialcte);
        readText::get(line, "readorder", readorder);
        readText::get(line, "adcerr", adcerr);
        readText::get(line, "obshistid", obshistid);
        readText::get(line, "exposureid", exposureid);
        readText::get(line, "chipid", chipid);
        readText::get(line, "obsseed", obsseed);
        readText::get(line, "vistime", vistime);
        readText::get(line, "nsnap", nsnap);
        readText::get(line, "nframes", nframes);
        readText::get(line, "nskip", nskip);
        readText::get(line, "parallelread", parallelread);
        readText::get(line, "serialread", serialread);
        readText::get(line, "parallelprescan", parallelPrescan);
        readText::get(line, "serialoverscan", serialOverscan);
        readText::get(line, "serialprescan", serialPrescan);
        readText::get(line, "paralleloverscan", parallelOverscan);
        readText::get(line, "bias", bias);
        readText::get(line, "gain", gain);
        readText::get(line, "readnoise", readnoise);
        readText::get(line, "darkcurrent", darkcurrent);
        readText::get(line, "hotpixelrate", hotpixelrate);
        readText::get(line, "hotcolumnrate", hotcolumnrate);
    }

    instr = "";
    unsigned pos = instrdir.rfind("/") + 1;
    for (unsigned i = pos; i < instrdir.length(); i++) {
        instr += instrdir[i];
    }

    if (flatdir == 1) instrdir = ".";

    if (adcerr.size() == 0) adcerr.resize(16, 0);
    
    std::istringstream focalplanePars(readText::get(instrdir + "/focalplanelayout.txt", chipid));
    double centerxt, centeryt, pixsizet;
    long pixelsxt, pixelsyt;
    focalplanePars >> centerxt>> centeryt >> pixsizet >> pixelsxt >> pixelsyt >> devmaterial >> devtype >> devmode >> devvalue;
    // allow devmode override
    focalplanefile = instrdir + "/segmentation.txt";
    std::vector<std::string> amplifiers;
    readText::readSegmentation(focalplanefile, chipid, amplifiers);
    namp = static_cast<long>(amplifiers.size());
    outchipid.resize(namp);
    outminx.resize(namp);
    outmaxx.resize(namp);
    outminy.resize(namp);
    outmaxy.resize(namp);
    crosstalk.resize(namp);
    dataBit.resize(namp);
    for (long j = 0; j < namp; j++) {
        std::istringstream segmentationPars(amplifiers[j]);
        segmentationPars >> outchipid[j] >> outminx[j] >> outmaxx[j] >> outminy[j] >> outmaxy[j];
        double v;
        for (int i = 0; i < 16; i++) {
            segmentationPars >> v;
        }
        segmentationPars >> dataBit[j];
        crosstalk[j].resize(namp);
        for (long k = 0; k < namp; k++) {
            segmentationPars >> crosstalk[j][k];
        }
    }

    long seed = obsseed + exposureid;
    if (obsseed == -1) {
        random.setSeedFromTime();
    } else {
        random.setSeed64(seed);
    }
    random.unwind(10000);

    if (devtype == "CMOS") {
        exptime = vistime/nsnap;
    } else {
        exptime = (vistime - (nsnap - 1)*devvalue)/nsnap;
    }

    int nsamples = nframes + nskip;
    if (nframes == 1 && nskip == 0) devmode = "frame";
    if (devmode == "sequence") {
        exptime = devvalue;
        ngroups = ceil((vistime + devvalue*nskip)/(devvalue*nsamples));
        emap.resize(ngroups + 1); //number of groups plus reference frame
    } else {
        ngroups = 1;
        nframes = 1;
        nskip = 0;
        nsamples = 1;
        emap.resize(1);
    }

    int nframei = 0; // current frame number within group
    int ngroupi = 0; // current group number

    // loop through each frame in integration sequence
    for (int nframe = 0; nframe < ngroups*nsamples - nskip; nframe++) {

        infile.clear();
        infile.str("");
        if (devmode == "sequence") {
            infile << instr << "_e_"  << obshistid << "_f" << filter << "_" << chipid << "_E" << std::setfill('0') << std::setw(3) << exposureid << "_R" << std::setfill('0') << std::setw(3) << nframe << ".fits.gz";
        } else {
            infile << instr << "_e_"  << obshistid << "_f" << filter << "_" << chipid << "_E" << std::setfill('0') << std::setw(3) << exposureid << ".fits.gz";
        }
        nframei++;
        // do not simulate frame skips
        if (nframei % nframes == 0) {
            nframe += nskip;
            nframei = 0;
            ngroupi++;
            if (ngroupi > 1) {
                emap[ngroupi].reserve(onaxes[0]*onaxes[1]);
                // continue "up the ramp" integration from previous group
                for (int i = 0; i < onaxes[0]; i++) {
                    for (int j = 0; j < onaxes[1]; j++) {
                        emap[ngroupi][onaxes[0]*j + i] = emap[ngroupi - 1][onaxes[0]*j + i];
                    }
                }
            }
        }
        int status = 0;
        std::string ss2;
        float nullval;
        int anynull;
        int nfound;
        ss2 = infile.str();
        std::vector<float> emaptemp;
        fits_open_file(&foptr, ss2.c_str(), READONLY, &status);
        fits_read_keys_lng(foptr, (char*)"NAXIS", 1, 2, onaxes, &nfound, &status);
        emaptemp.reserve(onaxes[0]*onaxes[1]);
        if ((nframe == 0 || nframei == 0) && devmode == "sequence") emap[ngroupi].assign(onaxes[0]*onaxes[1], 0);
        if (nframe == 0) emap[0].reserve(onaxes[0]*onaxes[1]);
        fits_read_img(foptr, TFLOAT, 1, onaxes[0]*onaxes[1], &nullval, const_cast<float *>(&emaptemp[0]), &anynull, &status);

        // average-combine group
        for (int i = 0; i < onaxes[0]; i++) {
            for (int j = 0; j < onaxes[1]; j++) {
                if (nframe == 0) { // reference frame (not actually a group frame)
                    emap[0][onaxes[0]*j + i] = emaptemp[onaxes[0]*j + i];
                } 
                if (devmode == "sequence") {
                    emap[ngroupi][onaxes[0]*j + i] += emaptemp[onaxes[0]*j + i]/nframes;
                }
            }
        }

    }
    
    adcmap.reserve(onaxes[0]*onaxes[1]);

    Uint32 seedchip = 0;
    for (size_t m(0); m < chipid.size(); m++) {
        seedchip += static_cast<Uint32>(((static_cast<int>(chipid.c_str()[m])%10)*pow(10, m)));
    }
    random.setSeed32Fixed(seedchip);

}

void E2adc::setHotpixels() {

    hotpixelmap.resize(onaxes[0]*onaxes[1], 0);

    // amplifier loop
    for (long l = 0; l < namp; l++) {
        for (long j = (outminy[l] - miny); j <= (outmaxy[l] - miny); j++) {
            for (long i = (outminx[l] - minx); i<= (outmaxx[l] - minx); i++) {
                if (random.uniformFixed() < hotpixelrate[l]) {
                    hotpixelmap[onaxes[0]*j + i] = 1;
                }
            }
        }
    }
    // amplifier loop (hot columns)
    for (long l = 0; l < namp; l++) {
        for (long j = (outminy[l] - miny); j <= (outmaxy[l] - miny); j++) {
            if (parallelread[l]==-1) {
                for (long i = (outminx[l] - minx); i <= (outmaxx[l] - minx); i++) {
                    if (random.uniformFixed() < hotcolumnrate[l]/(outmaxx[l] - outminx[l] + 1)*2) {
                        for (long ii = i; ii <= (outmaxx[l] - minx); ii++) {
                            hotpixelmap[onaxes[0]*j + ii] = 1;
                        }
                        break;
                    }
                }
            } else {
                for (long i = (outmaxx[l] - minx); i >= (outminx[l] - minx); i--) {
                    if (random.uniformFixed() < hotcolumnrate[l]/(outmaxx[l] - outminx[l] + 1)*2) {
                        for (long ii = i; ii >= (outminx[l] - minx); ii--) {
                            hotpixelmap[onaxes[0]*j + ii] = 1;
                        }
                        break;
                    }
                }
            }
        }
    }
}

void E2adc::convertADC() {

    int nstop;
    if (devmode == "sequence") {
        fullReadoutMap.resize(ngroups + 1);
        fullReadoutMapL.resize(ngroups + 1);
        nstop = ngroups + 1;
    } else {
        fullReadoutMap.resize(1);
        fullReadoutMapL.resize(1);
        nstop = 1;
    }

    std::vector<std::vector<float> > readoutmap;
    readoutmap.resize(namp);
    nx.resize(namp);
    ny.resize(namp);
    
    // loop through integration sequence groups
    for (int ngroup = 0; ngroup < nstop; ngroup++) {

        // amplifier loop
        for (long l = 0; l < namp; l++) {
            std::cout << "Reading out chip " << chipid << " with amplifier chain " << outchipid[l] << std::endl;
            double nbackground = darkcurrent[l]*exptime;
            for (long i = (outminx[l] - minx); i <= (outmaxx[l] - minx); i++) {
                for (long j = (outminy[l] - miny); j <= (outmaxy[l] - miny); j++) {
                    // readnoise and dark current
                    long factor = roundl(random.normal()*sqrt(nbackground) + nbackground);
                    if (factor < 0) factor = 0;
                    adcmap[onaxes[0]*j + i] = emap[ngroup][onaxes[0]*j + i] + factor;
                    if (adcmap[onaxes[0]*j + i] > wellDepth) adcmap[onaxes[0]*j + i] = wellDepth;
                }
            }

            for (long i = (outminx[l] - minx); i <= (outmaxx[l] - minx); i++) {
                for (long j = (outminy[l] - miny); j <= (outmaxy[l] - miny); j++) {
                    if (hotpixelmap[onaxes[0]*j + i] == 1) {
                        adcmap[onaxes[0]*j + i] = wellDepth;
                    }
                }
            }

            // small adc map
            long naxes[2];
            naxes[0] = outmaxx[l] - outminx[l] + 1;
            naxes[1] = outmaxy[l] - outminy[l] + 1;
            std::vector<float> smalladcmap(naxes[0]*naxes[1]);
            for (long i = 0; i < naxes[0]; i++) {
                for (long j = 0; j < naxes[1]; j++) {
                    smalladcmap[naxes[0]*j + i] = adcmap[onaxes[0]*(j + (outminy[l] - miny))+(i + (outminx[l] - minx))];
                }
            }

            // Charge transfer inefficiency   a*b+(1-a)*b+(1-b)=a*b+b-a*b+1-b=1
            if (parallelread[l] == 1 && serialread[l] == -1) {
                for (long i = naxes[0] - 1; i >= 0; i--) {
                    for (long j = 0; j < naxes[1]; j++) {
                        float origct = smalladcmap[naxes[0]*j + i];
                        for (long k = 0; k < (static_cast<long>(origct)); k++) {
                            if (random.uniform() < parallelcte) {
                                if (random.uniform() >= serialcte) {
                                    if (j != naxes[1] - 1) {
                                        smalladcmap[naxes[0]*(j + 1) + i] += 1.0;
                                        smalladcmap[naxes[0]*j + i] -= 1.0;
                                    }
                                }
                            } else {
                                if (i != 0) {
                                    smalladcmap[naxes[0]*j + (i - 1)] += 1.0;
                                    smalladcmap[naxes[0]*j + i] -= 1.0;
                                }
                            }
                        }
                    }
                }
            }

            if (parallelread[l] == 1 && serialread[l] == 1) {
                for (long i = naxes[0] - 1; i >= 0; i--) {
                    for (long j = naxes[1] - 1; j >= 0; j--) {
                        float origct = smalladcmap[naxes[0]*j + i];
                        for (long k = 0; k < (static_cast<long>(origct)); k++) {
                            if (random.uniform() < parallelcte) {
                                if (random.uniform() >= serialcte) {
                                    if (j != 0) {
                                        smalladcmap[naxes[0]*(j - 1) + i] += 1.0;
                                        smalladcmap[naxes[0]*j + i] -= 1.0;
                                    }
                                }
                            } else {
                                if (i != 0) {
                                    smalladcmap[naxes[0]*j + (i - 1)] += 1.0;
                                    smalladcmap[naxes[0]*j + i] -= 1.0;
                                }
                            }
                        }
                    }
                }
            }

            if (parallelread[l] == -1 && serialread[l] == -1) {
                for (long i = 0; i < naxes[0]; i++) {
                    for (long j = 0; j < naxes[1]; j++) {
                        float origct = smalladcmap[naxes[0]*j + i];
                        for (long k = 0; k < (static_cast<long>(origct)); k++) {
                            if (random.uniform() < parallelcte) {
                                if (random.uniform() >= serialcte) {
                                    if (j != naxes[1] - 1) {
                                        smalladcmap[naxes[0]*(j + 1) + i] += 1.0;
                                        smalladcmap[naxes[0]*j + i] -= 1.0;
                                    }
                                }
                            } else {
                                if (i != naxes[0] - 1) {
                                    smalladcmap[naxes[0]*j + (i + 1)] += 1.0;
                                    smalladcmap[naxes[0]*j + i] -= 1.0;
                                }
                            }
                        }
                    }
                }
            }

            if (parallelread[l] == -1 && serialread[l] == 1) {
                for (long i = 0; i < naxes[0]; i++) {
                    for (long j = naxes[1] - 1; j >= 0; j--) {
                        float origct = smalladcmap[naxes[0]*j + i];
                        for (long k = 0; k < (static_cast<long>(origct)); k++) {
                            if (random.uniform() < parallelcte) {
                                if (random.uniform() >= serialcte) {
                                    if (j != 0) {
                                        smalladcmap[naxes[0]*(j - 1) + i] += 1.0;
                                        smalladcmap[naxes[0]*j + i] -= 1.0;
                                    }
                                }
                            } else {
                                if (i != naxes[0] - 1) {
                                    smalladcmap[naxes[0]*j + (i + 1)] += 1.0;
                                    smalladcmap[naxes[0]*j + i] -= 1.0;
                                }
                            }
                        }
                    }
                }
            }

            // electron to ADC conversion
            for (long i = 0; i < naxes[0]; i++) {
                for (long j = 0; j < naxes[1]; j++) {
                    smalladcmap[naxes[0]*j + i] = smalladcmap[naxes[0]*j + i]/(gain[l])/
                        (1 + nonlinear*(smalladcmap[naxes[0]*j + i]/wellDepth)) + bias[l] + roundl(random.normal()*readnoise[l]);
                }
            }

            nx[l] = naxes[0] + serialOverscan[l] + serialPrescan[l];
            ny[l] = naxes[1] + parallelPrescan[l] + parallelOverscan[l];
            std::vector<float> fulladcmap(nx[l]*ny[l]);

            if (parallelread[l] == 1 && serialread[l] == -1) {
                for (long i = 0; i < naxes[0]; i++) {
                    for (long j = 0; j < naxes[1]; j++) {
                        fulladcmap[nx[l]*(j + parallelPrescan[l]) + i + serialPrescan[l]] = smalladcmap[naxes[0]*j + i];
                    }
                }
                for (long i = (nx[l] - serialOverscan[l]); i < nx[l]; i++) {
                    for (long j = 0; j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < serialPrescan[l]; i++) {
                    for (long j = 0; j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < nx[l]; i++) {
                    for (long j = 0; j < parallelPrescan[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < nx[l]; i++) {
                    for (long j = (ny[l] - parallelOverscan[l]); j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
            } else if (parallelread[l] == -1 && serialread[l] == -1) {
                for (long i = 0; i < naxes[0]; i++) {
                    for (long j = 0; j < naxes[1]; j++) {
                        fulladcmap[nx[l]*(j + parallelPrescan[l]) + i + serialOverscan[l]] = smalladcmap[naxes[0]*j + i];
                    }
                }
                for (long i = 0; i < serialOverscan[l]; i++) {
                    for (long j = 0; j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = (nx[l] - serialPrescan[l]); i < nx[l]; i++) {
                    for (long j = 0; j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < nx[l]; i++) {
                    for (long j = 0; j < parallelPrescan[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < nx[l]; i++) {
                    for (long j = (ny[l] - parallelOverscan[l]); j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
            } else if (parallelread[l] == 1 && serialread[l] == 1) {
                for (long i = 0; i < naxes[0]; i++) {
                    for (long j = 0; j < naxes[1]; j++) {
                        fulladcmap[nx[l]*(j + parallelOverscan[l]) + i + serialPrescan[l]] = smalladcmap[naxes[0]*j + i];
                    }
                }
                for (long i = (nx[l] - serialOverscan[l]); i < nx[l]; i++) {
                    for (long j = 0; j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < serialPrescan[l]; i++) {
                    for (long j = 0; j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < nx[l]; i++) {
                    for (long j = 0; j < parallelOverscan[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < nx[l]; i++) {
                    for (long j = (ny[l] - parallelPrescan[l]); j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
            } else if (parallelread[l] == -1 && serialread[l] == 1) {
                for (long i = 0; i < naxes[0]; i++) {
                    for (long j = 0; j < naxes[1]; j++) {
                        fulladcmap[nx[l]*(j + parallelOverscan[l]) + i + serialOverscan[l]] = smalladcmap[naxes[0]*j + i];
                    }
                }
                for (long i = 0; i < serialOverscan[l]; i++) {
                    for (long j = 0; j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = (nx[l] - serialPrescan[l]); i < nx[l]; i++) {
                    for (long j = 0; j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < nx[l]; i++) {
                    for (long j = 0; j < parallelOverscan[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
                for (long i = 0; i < nx[l]; i++) {
                    for (long j = (ny[l] - parallelPrescan[l]); j < ny[l]; j++) {
                        fulladcmap[nx[l]*j + i] = bias[l] + roundl(random.normal()*readnoise[l]);
                    }
                }
            }

            // ADC digitization
            std::vector<float> readoutmap_orig(nx[l]*ny[l]);
            for (long i = 0; i < nx[l]; i++) {
                for (long j = 0; j < ny[l]; j++) {
                    readoutmap_orig[nx[l]*j + i] = 0.0;
                    for (int k = 0; k < 16; k++) {
                        long bit = ((static_cast<int>(fulladcmap[nx[l]*j + i]/(pow(2, k) + adcerr[k]))) % 2);
                        readoutmap_orig[nx[l]*j + i] += bit*pow(2, k);
                    }
                }
            }

            // change to readout order
            readoutmap[l].resize(nx[l]*ny[l]);
            if (readorder == 1) {
                if (parallelread[l] == 1 && serialread[l] == -1) {
                    for (long i = 0; i < nx[l]; i++) {
                        for (long j = 0; j < ny[l]; j++) {
                            readoutmap[l][ny[l]*i + j] = readoutmap_orig[nx[l]*j + (nx[l] - 1 - i)];
                        }
                    }
                } else if (parallelread[l] == -1 && serialread[l] == -1) {
                    for (long i = 0; i < nx[l]; i++) {
                        for (long j = 0; j < ny[l]; j++) {
                            readoutmap[l][ny[l]*i + j] = readoutmap_orig[nx[l]*j + i];
                        }
                    }
                } else if (parallelread[l] == 1 && serialread[l] == 1) {
                    for (long i = 0; i < nx[l]; i++) {
                        for (long j = 0; j < ny[l]; j++) {
                            readoutmap[l][ny[l]*i + j] = readoutmap_orig[nx[l]*(ny[l] - 1 - j)+(nx[l] - 1 - i)];
                        }
                    }
                } else if (parallelread[l] == -1 && serialread[l] == 1) {
                    for (long i = 0; i < nx[l]; i++) {
                        for (long j = 0; j < ny[l]; j++) {
                            readoutmap[l][ny[l]*i + j] = readoutmap_orig[nx[l]*(ny[l] - 1 - j) + i];
                        }
                    }
                }
            } else {
                for (long i = 0; i < nx[l]*ny[l]; i++) {
                    readoutmap[l][i] = readoutmap_orig[i];
                }
            }
        }

        // crosstalk
        fullReadoutMap[ngroup].resize(namp);
        fullReadoutMapL[ngroup].resize(namp);
        for (long l = 0; l < namp; l++) {
            size_t mapSize = readoutmap[l].size();
            if (dataBit[l] <= 16) {
                fullReadoutMap[ngroup][l].resize(mapSize);
            } else {
                fullReadoutMapL[ngroup][l].resize(mapSize);
            }
            for (size_t i = 0; i < mapSize; i++) {
                float sum = 0.0;
                for (long k = 0; k < namp; k++) {
                    if (i < readoutmap[k].size()) {
                        sum += crosstalk[l][k]*readoutmap[k][i];
                    }
                }
                if (dataBit[l] <= 16) {
                    if (sum >= 0 && sum < pow(2, 16)) fullReadoutMap[ngroup][l][i] = static_cast<unsigned short>(sum);
                    if (sum < 0) fullReadoutMap[ngroup][l][i] = static_cast<unsigned short>(0);
                    if (sum >= pow(2, 16)) fullReadoutMap[ngroup][l][i] = static_cast<unsigned short>(pow(2, 16) - 1);
                } else {
                    if (sum >= 0 && sum < pow(2, 32)) fullReadoutMapL[ngroup][l][i] = static_cast<unsigned long>(sum);
                    if (sum < 0) fullReadoutMapL[ngroup][l][i] = static_cast<unsigned long>(0);
                    if (sum >= pow(2, 32)) fullReadoutMapL[ngroup][l][i] = static_cast<unsigned long>(pow(2, 32) - 1);
                }
            }
        }
    }

}

void E2adc::writeFitsImage() {

    char filename[4096];
    char comment[4096];
    fitsfile *fptr = NULL;

    int nstop;
    if (devmode == "sequence") {
        nstop = ngroups + 1;
    } else {
        nstop = 1;
    }

    // loop through integration sequence groups
    for (int ngroup = 0; ngroup < nstop; ngroup++) {
        std::string tarFiles(infile.str());
        // amplifier loop
        for (long l = 0; l < namp; l++) {
            // create file
            std::ostringstream outfile;
            outfile << "!" << instr << "_a_" << obshistid << "_f"<< filter << "_" << outchipid[l] << "_E";
            outfile << std::setfill('0') << std::setw(3) << exposureid;
            if (devmode == "sequence") {
                outfile << "_G" << std::setfill('0') << std::setw(3) << ngroup << ".fits.gz";
            } else {
                outfile << ".fits.gz";
            }
            std::string sss;
            sss = outfile.str();
            sprintf(filename, "%s", sss.c_str());
            int status = 0;
            fits_create_file(&fptr, filename, &status);
            // write header keys
            fits_copy_header(foptr, fptr, &status);
            if (dataBit[l] <= 16) {
                long temp = 16;
                fitsUpdateKey(fptr, "BITPIX", temp, "number of bits per data pixel");
                unsigned long temp2 = 32768;
                fitsUpdateKey(fptr, "BZERO", temp2, "offset data range to that of unsigned short");
            } else {
                long temp = 32;
                fitsUpdateKey(fptr, "BITPIX", temp, "number of bits per data pixel");
                unsigned long temp2 = 2147483648;
                fitsUpdateKey(fptr, "BZERO", temp2, "offset data range to that of unsigned long");
            }
            fitsUpdateKey(fptr, "BSCALE", static_cast<long>(1), "default scaling factor");
            fitsWriteKey(fptr, "BIAS", bias[l], "Bias");
            fitsWriteKey(fptr, "GAIN", gain[l], "Gain");
            fitsWriteKey(fptr, "SCTE", serialcte, "Serial CTE");
            fitsWriteKey(fptr, "PCTE", parallelcte, "Parallel CTE");
            fitsWriteKey(fptr, "NONLIN", nonlinear, "Non-linear gain");
            fitsWriteKey(fptr, "E2AWLDP", wellDepth, "E2adc well depth");
            fitsWriteKey(fptr, "SATURATE", wellDepth/gain[l] + bias[l], "Saturation estimate");
            fitsWriteKey(fptr, "PREAD", parallelread[l], "Parallel read out direction");
            fitsWriteKey(fptr, "SREAD", serialread[l], "Serial read out direction");
            fitsWriteKey(fptr, "PSCANP", parallelPrescan[l], "Pre-scan parallel");
            fitsWriteKey(fptr, "OSCANS", serialOverscan[l], "Over-scan serial");
            fitsWriteKey(fptr, "PSCANS", serialPrescan[l], "Pre-scan serial");
            fitsWriteKey(fptr, "OSCANP", parallelOverscan[l], "Over-scan parallel");
            for (int i = 0; i < 16; i++) {
                std::ostringstream ss;
                ss << "ADCER" << i;
                std::string ss2;
                ss2 = ss.str();
                fitsWriteKey(fptr, ss2.c_str(), adcerr[i], "ADC error bit");
            }
            fitsWriteKey(fptr, "E2AICHI", chipid.c_str(), "E2adc input chip ID");
            fitsWriteKey(fptr, "E2AOCHI", outchipid[l].c_str(), "E2adc output chip ID");
            unsigned pos = outchipid[l].find_last_of("_");
            fitsWriteKey(fptr, "CCDID", outchipid[l].substr(0, pos).c_str(), "CCD ID");
            fitsWriteKey(fptr, "AMPID", outchipid[l].substr(pos + 1).c_str(), "Amplifier ID");
            fitsWriteKey(fptr, "RDORDER", readorder, "0=CCS; 1=readorder");
            fitsWriteKey(fptr, "RDNOISE", readnoise[l]/gain[l], "Readout noise (ADU/pixel)");
            fitsWriteKey(fptr, "DRKCURR", darkcurrent[l], "Dark Current (e-/pixel/s)");
            if (readorder == 1) {
                std::ostringstream ss;
                std::string ss2;
                ss.str("");
                ss << "[" << std::setw(4) << 1 << ":" << std::setw(4) << ny[l];
                ss << "," << std::setw(4) << 1 << ":" << std::setw(4) << serialOverscan[l] << "]";
                ss2 = ss.str();
                fitsWriteKey(fptr, "BIASSEC", ss2.c_str(), "Scan section of amplifier");
                ss.str("");
                ss << "[" << std::setw(4) << parallelPrescan[l] + 1 << ":" << std::setw(4) << ny[l] - parallelOverscan[l];
                ss << "," << std::setw(4) << serialOverscan[l] + 1 << ":" << std::setw(4) << nx[l] - serialPrescan[l] << "]";
                ss2 = ss.str();
                fitsWriteKey(fptr, "TRIMSEC", ss2.c_str(), "Trimmed section of amplifier");
                fitsWriteKey(fptr, "DATASEC", ss2.c_str(), "Data section of amplifier");
            }
            double value = 0.0;
            fitsReadKey(fptr, "CRPIX1", &value, comment);
            fitsUpdateKey(fptr, "CRPIX1", value - outminx[l], comment);
            fitsReadKey(fptr, "CRPIX2", &value, comment);
            fitsUpdateKey(fptr, "CRPIX2", value - outminy[l], comment);
            if (readorder == 1) {
                double key[8];
                double newKey[8];
                fitsReadKey(fptr, "CRPIX1", &key[0], comment);
                fitsReadKey(fptr, "CRVAL1", &key[1], comment);
                fitsReadKey(fptr, "CRPIX2", &key[2], comment);
                fitsReadKey(fptr, "CRVAL2", &key[3], comment);
                fitsReadKey(fptr, "CD1_1", &key[4], comment);
                fitsReadKey(fptr, "CD1_2", &key[5], comment);
                fitsReadKey(fptr, "CD2_1", &key[6], comment);
                fitsReadKey(fptr, "CD2_2", &key[7], comment);
                // x->y
                newKey[0] = key[2];
                newKey[1] = key[1];
                newKey[2] = key[0];
                newKey[3] = key[3];
                newKey[4] = key[5];
                newKey[5] = key[4];
                newKey[6] = key[7];
                newKey[7] = key[6];
                // flip x
                if (serialread[l] == 1) {
                    newKey[0] = -newKey[0] + ny[l] - 1;
                    newKey[4] = -newKey[4];
                    newKey[6] = -newKey[6];
                }
                // flip y
                if (parallelread[l] == 1) {
                    newKey[2] = -newKey[2] + nx[l] - 1;
                    newKey[5] = -newKey[5];
                    newKey[7] = -newKey[7];
                }
                fitsUpdateKey(fptr, "CRPIX1", newKey[0], comment);
                fitsUpdateKey(fptr, "CRVAL1", newKey[1], comment);
                fitsUpdateKey(fptr, "CRPIX2", newKey[2], comment);
                fitsUpdateKey(fptr, "CRVAL2", newKey[3], comment);
                fitsUpdateKey(fptr, "CD1_1", newKey[4], comment);
                fitsUpdateKey(fptr, "CD1_2", newKey[5], comment);
                fitsUpdateKey(fptr, "CD2_1", newKey[6], comment);
                fitsUpdateKey(fptr, "CD2_2", newKey[7], comment);
            }
            // write image data
            if (dataBit[l] <= 16) {
                fitsWriteImage(fptr, ny[l], nx[l], const_cast<unsigned short *>(&fullReadoutMap[ngroup][l][0]));
            } else {
                fitsWriteImage(fptr, ny[l], nx[l], const_cast<unsigned long *>(&fullReadoutMapL[ngroup][l][0]));
            }
            tarFiles += " " + outfile.str();
        }

        if (tarfile == 1) {
            std::ostringstream tarName;
            tarName << instr << "_" << obshistid << "_f" << filter << "_" << chipid << "_E" << std::setfill('0') << std::setw(3) << exposureid << std::setfill('0') << std::setw(3) << ngroup << ".tar";
            std::cout << "Tarring " << tarName.str() << std::endl;
            std::string tarCommand = "tar cf " + tarName.str() + " " + tarFiles + " --remove-files";
            system(tarCommand.c_str());
        }
    }

}


int main(void) {

    E2adc e2adc;

    e2adc.setup();
    e2adc.setHotpixels();
    e2adc.convertADC();
    e2adc.writeFitsImage();

    return 0;

}
