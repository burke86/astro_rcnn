///
/// @package phosim
/// @file photonloop.cpp
/// @brief main photon loop
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///


void* Image::threadFunction(void *voidArgs) {
    thread_args *args = (thread_args*)voidArgs;
    for (int i=0; i<args->runthread-1; i++) {
        args->instance->photonLoop(args->ssource[i], args->thread,0);
    }
    args->instance->photonLoop(args->ssource[args->runthread-1], args->thread,args->runthread);
    return NULL;
}

void Image::photonLoop(long ssource, long thread, int finish) {

    Vector angle, largeAngle, position, positionDiffraction;
    Vector angleDiffraction, angleDeltaDiffraction, angleDeltaTotalDiffraction;
    Vector positionPrevious, anglePrevious, normal;
    double transmission, moonTransmission, reflection, distance;
    double shiftX = 0.0, shiftY = 0.0;
    long long ray;
    long long photonSource;
    long long detRay;
    long long backSourceOver;
    double sourceSaturation;
    double backSigma = 0.0;
    long preGhost = 0;
    long sourceOver;
    long newSurf = 0;
    long oldSurf;
    long waveIndex;
    long waveSurfaceIndex;
    long straylightcurrent;
    long surfaceLimit;
    int miss;
    int initFlag;
    int initGal;
    Photon photon;
    Clog localLog;
    double rate = 0.1;
    double photonSourceDouble;
    int lock4=0;

    //    SETUP PHOTONS
    photon.thread = thread;
    counterClear(&localLog);
    photon.maxcounter = 0;
    surfaceLimit = natmospherefile*2 + nsurf*2 + 2;
    photonSourceDouble = nphot*sources.norm[ssource]/totalnorm;
    if (nsource == 1) photonSourceDouble = nphot;
    if (poissonMode == 1){
        photonSource = random[thread].poisson(photonSourceDouble);
    } else {
        photonSource = static_cast<long long>(photonSourceDouble);
    }
    if (telconfig != 0 && sources.type[ssource] != 0) photonSource = 0;
    // if (aperturemode == 2) photonSource = OPD_SCREEN_SIZE*OPD_SCREEN_SIZE;

    if (opdfile == 1) {
        pthread_mutex_lock(&lock.lock4);
        lock4 = 1;
        for (int surfIdx = 0; surfIdx < nsurf; surfIdx++) {
            surface.innerRadius[surfIdx] = 0.0;
            surface.asphere(surfIdx, SURFACE_POINTS);
        }
    }

    ray = 0;
    initGal = 0;
    sourceOver = 1;
    photon.sourceOver_m = 1;
    sourceSaturation = 1.0;
    photon.sourceSaturationRadius = 0;
    obstruction.pupil = 0;
    photon.prtime = -1.0;
    if (sources.type[ssource] < 6) {
        rate = (static_cast<double>(state.globalLog.accepted) + 10)/(static_cast<double>(state.globalLog.totalPhoton) + 1000);
        if (rate < 0.01) rate = 0.01;
        if (backGamma < 1.0) backSourceOver = (long long)(backAlpha*sqrt(rate*photonSource));
        else backSourceOver = (long long)(backAlpha/backGamma*sqrt(rate*photonSource));
        if (backAlpha <= 0.0) backSourceOver = static_cast<long long>(backGamma);
        if (backSourceOver < 1) backSourceOver = 1;
        if (backSourceOver > photonSource) backSourceOver = photonSource;
    } else {
        backSourceOver = 1;
    }
    if (sources.type[ssource] >= 6) {
        backSigma = 0.0;
    } else {
        backSigma = backBeta*backRadius/3600.0*platescale/1000.0;
    }
    if (sources.mag[ssource] < straylightcut && straylight == 1) {
        straylightcurrent = 1;
    } else {
        straylightcurrent = 0;
    }
    //   MAIN PHOTON LOOP
    photonloop:
    while (ray < photonSource) {
        sourceOver = 1;
        photon.sourceOver_m = 1;

        //   Get Wavelength and Time
        miss = getWavelengthTime(&photon, ssource);
        photon.absoluteTime = photon.time - exptime/2.0 + timeoffset;

        // Get Angle
        getAngle(&angle, photon.time, ssource);

        if (eventfile) {
            state.pEventLogging->logPhoton(angle.x, angle.y, photon.wavelength, 0);
            state.pEventLogging->logPhoton(photon.time, stof(sources.id[ssource]), 0, 1);
        }

        //   Saturation Optimization
        if (saturation) {
            if (sourceSaturation > 1.0 && photon.sourceSaturationRadius > 0 && sources.type[ssource] >= 6) {
                sourceOver = floor(sourceSaturation);
                if (random[thread].uniform() < (sourceSaturation - floor(sourceSaturation))) sourceOver++;
                if (sourceOver > well_depth) sourceOver = well_depth;
                if (sourceOver < 1) sourceOver = 1;
                photon.sourceOver_m = floor(((sourceSaturation - np)/(1 - np))*1.00);
                if (random[thread].uniform() < ((sourceSaturation - np)/(1 - np)*1.00) - floor((sourceSaturation - np)/(1 - np))*1.00) photon.sourceOver_m++;
                if (photon.sourceOver_m < 1) {
                    photon.sourceOver_m = 1;
                    sourceOver = 1;
                }
                // if ((ray%10000000) ==0) printf("%lf %ld %ld %lf %lld %lld\n",sourceSaturation,sourceOver,photon.sourceOver_m,photon.sourceSaturationRadius,ray,photonSource);
            }
        }
        // if (((ray%10000000)==0) && (ray!=0)) printf("%lf |%s| %lld %ld %ld %ld\n",sources.mag[ssource],sources.id[ssource].c_str(),ray,sourceOver,photon.sourceOver_m,photon.sourceSaturationRadius);
        if (sources.type[ssource] < 6 && backGamma > 1.0) {
            sourceOver = static_cast<long long>(backGamma);
            photon.sourceOver_m = sourceOver;
            if (backSourceOver*sourceOver > photonSource) {
                sourceOver = static_cast<long long>(photonSource/backSourceOver);
                photon.sourceOver_m = static_cast<long long>(photonSource/backSourceOver);
            }
            if (sourceOver < 1) {
                sourceOver = 1;
                photon.sourceOver_m = 1;
            }
        }
        if (miss) {
            countBad(&localLog, sourceOver*backSourceOver, &ray);
            goto photonloop;
        }
        waveIndex = static_cast<long>(photon.wavelength*1000 - minwavelength);
        waveSurfaceIndex = waveIndex*surfaceLimit;
        initFlag = 0;
        if (throughputfile) addThroughput(&state.throughputLog, minwavelength, maxwavelength, -1, waveIndex, sourceOver*backSourceOver);

        // Dynamic Transmission Optimization
        long kmax;
        if (sources.type[ssource] < 6 && backGamma > 1.0) {
            kmax = 1;
        } else {
            kmax = sourceOver;
        }
        // if ((ray % 1000000)==0) printf("yy %lld %lf %ld %ld %ld\n",ray,sourceSaturation,sourceOver,photon.sourceOver_m,photon.sourceSaturationRadius);
        for (long k = 0; k < kmax; k++) {

            if ((k == 0) || (k > 0 && straylightcurrent == 1)) {
                long lastSurface = -1;
                miss = dynamicTransmissionOptimization(k, &lastSurface, &preGhost, waveSurfaceIndex, straylightcurrent, &photon);
                if (miss == 1) {
                    if (throughputfile && lastSurface >= 0) {
                        for (long j = 0; j <= lastSurface; j++) {
                            addThroughput(&state.throughputLog, minwavelength, maxwavelength, j, waveIndex, sourceOver*backSourceOver);
                        }
                    }
                    countBad_dt(&localLog, sourceOver*backSourceOver, &ray);

                    goto photonloop;
                }
            }

        redodiff:;
            vectorInit(&largeAngle);
            shiftX = 0.0;
            shiftY = 0.0;
            //  Sample Pupil

            miss = samplePupil(&positionDiffraction, ray, &photon);
            if (miss) {
                if (k > 0) goto redodiff;
                countBad(&localLog, sourceOver*backSourceOver, &ray);

                goto photonloop;
            }

            miss = samplePupil(&position, ray, &photon);
            if (miss) {
                if (k > 0) goto redodiff;
                countBad(&localLog, sourceOver*backSourceOver, &ray);

                goto photonloop;
            }
            if (finiteDistance != 0.0) {
                double finiteDistanceR = sqrt(pow(finiteDistance, 2) + pow(position.x, 2) +
                                              pow(position.y, 2));
                angle.x += position.x/finiteDistanceR;
                angle.y += position.y/finiteDistanceR;
                angle.z = smallAnglePupilNormalize(angle.x, angle.y);
            }
            if (opdfile && ray == 0) {
                angle.x = 0.0;
                angle.y = 1e-5;
                angle.z = smallAnglePupilNormalize(angle.x, angle.y);
            }
            if (diffractionMode == 5) vectorCopy(position, &positionDiffraction);
            photon.xp = position.x;
            photon.yp = position.y;
            vectorInit(&largeAngle);
            if (initFlag == 0) {
                photon.shiftedAngle = spiderangle + photon.time*rotationrate*ARCSEC;
                photon.wavelengthFactor = pow(photon.wavelength, -0.2)/screen.wavelengthfactor_nom;
                initFlag = 1;
            }

            //  Diffraction
            if (diffractionMode >= 1 && spiderMode == 1) {
                miss = diffraction(&positionDiffraction, angle, &largeAngle, &photon);
                if (miss) {
                    if (k > 0) goto redodiff;
                    countBad(&localLog, sourceOver*backSourceOver, &ray);

                    goto photonloop;
                }
            }
            //   Large Angle Scattering
            largeAngleScattering(&largeAngle, &photon);

            getDeltaAngle(&angle, &position, ssource, &shiftX, &shiftY, thread, &initGal, &photon);

            //   Second Kick
            if (diffractionMode == 1 && sources.type[ssource] !=  0) secondKick(&largeAngle, &photon);

            // Saturation Optimization Calculation
            if (photon.sourceSaturationRadius > 0) {
                // if ((ray % 1000000)==0) printf("xx %lld %lf %ld %ld %ld %lf %ld\n",ray,sourceSaturation,sourceOver,photon.sourceOver_m,photon.sourceSaturationRadius,modulus(&largeAngle)/DEGREE*platescale/pixsize,preGhost);
                // if (modulus(&largeAngle)/DEGREE*platescale/pixsize
                //     > photon.sourceSaturationRadius || preGhost>= 2) {
                if (sqrt(pow(largeAngle.x+shiftX,2.0)+pow(largeAngle.y+shiftY,2.0))/DEGREE*platescale/pixsize
                     > photon.sourceSaturationRadius || preGhost>= 2) {
                    photon.sourceOver_m = 1;
                    sourceSaturation -= ((1.0 - np)/np + 1.0)*0.01;
                    if (sourceSaturation < 1) sourceSaturation = 1;
                    break;
                }
            }

        }
        if (photon.sourceSaturationRadius > 0) sourceSaturation += 0.01;
        if (np <= 0.0) sourceSaturation = 1;
        // if (photon.sourceSaturationRadius > 0.0) printf("xx %lf %ld %ld %lf\n",sourceSaturation,sourceOver,photon.sourceOver_m,photon.sourceSaturationRadius);
        photon.counter = -1;
        if (opdfile) photon.op = 0.0;

        // Get Delta Angle
        // getDeltaAngle(&angle, &position, ssource, &shiftX, &shiftY);

       // ATMOSPHERE
        photon.airRefraction = 0.0;
        photon.airRefractionPrevious = 0.0;
        photon.ncurr = 1.0;
        // Loop over Atmosphere Layers
        moonTransmission = 1.0;
        vectorInit(&angleDeltaTotalDiffraction);
        for (int layer = -1; layer < natmospherefile; layer++) {

            // Atmospheric Dispersion
            photon.airRefractionPrevious = photon.airRefraction;
            photon.airRefraction = airIndexRefraction(&photon, layer);
            photon.ncurr = 1.0 + photon.airRefraction;
            if (sources.type[ssource] !=  0) {
                vectorCopy(angle,&angleDiffraction);
                atmosphericDispersion(&angleDiffraction, &photon, layer);
                vectorSubtract(angleDiffraction,angle,&angleDeltaDiffraction);
                vectorAdd(angleDeltaTotalDiffraction,angleDeltaDiffraction,&angleDeltaTotalDiffraction);
            }

             // Atmosphere Propagate
            atmospherePropagate(&position, angle, layer, diffractionMode, &photon);

            if (sources.type[ssource] != 0) {
                if (atmospheremode == 2) {
                    if (layer >= 0) {

                        // Atmosphere Intercept
                        atmosphereIntercept(&position, layer, &photon);

                        // Atmosphere Refraction
                        atmosphereRefraction(&angle, layer, diffractionMode, &photon);

                        // Clouds
                        transmission = cloudOpacity(layer, &photon);
                        if (sources.type[ssource] == 4) {
                            moonTransmission *= transmission;
                            transmission = cloudOpacityMoon(layer, &photon);
                        }
                        if (sources.type[ssource] == 3) {
                            transmission = cloudOpacityMoon(layer, &photon);
                        }
                        if (transmissionCheck(transmission, 1 + layer*2, waveSurfaceIndex, &photon)) {
                            countBad(&localLog, sourceOver*backSourceOver, &ray);
                            goto photonloop;
                        }
                    }

                    // Atmosphere Opacity
                    transmission = atmosphereOpacity(angle, layer, &photon);
                    if (sources.type[ssource] == 3) {
                        moonTransmission *= transmission;
                        transmission = atmosphereOpacityMoon(angle, layer, &photon);
                    }
                    if (sources.type[ssource] == 4) {
                        transmission = atmosphereOpacityMoon(angle, layer, &photon);
                    }
                    if (transmissionCheck(transmission, 2 + layer*2, waveSurfaceIndex, &photon)) {
                        countBad(&localLog, sourceOver*backSourceOver, &ray);
                        goto photonloop;
                    }

                }

                if (eventfile) {
                    state.pEventLogging->logPhoton(position.x, position.y, position.z, layer + 100);
                }

            }
        }
        if (sources.type[ssource] == 3 || sources.type[ssource] == 4) {
            moonTransmission = (1.0 - moonTransmission);
            double randNumber = random[thread].uniform();
            if (randNumber > moonTransmission) {
                countBad(&localLog, sourceOver*backSourceOver, &ray);
                goto photonloop;
            }
        }
        vectorAdd(angle,angleDeltaTotalDiffraction,&angle);

        // Atmosphere Diffraction
        if (diffractionMode == 2 && sources.type[ssource] != 0) atmosphereDiffraction(&angle, &photon);

        // Dome Seeing
        if (domeseeing > 0.0 || toypsf > 0.0) domeSeeing(&angle, &photon);

        // Tracking
        if (trackingMode) tracking(&angle, photon.absoluteTime);

        // Large Angle
        angle.x += largeAngle.x + shiftX;
        angle.y += largeAngle.y + shiftY;
        angle.z = smallAnglePupilNormalize(angle.x, angle.y);

        if (telescopeMode == 0) {
            newSurf = nsurf - 2;
            double p0 = platescale/DEGREE/1000/fabs(angle.z);
            if (nmirror % 2 == 0) {
                propagate(&position, angle, ((surface.height[nsurf - 1] + p0) - position.z)/angle.z);
                angle.x -= (position.x/p0);
                angle.y -= (position.y/p0);
                angle.z = -1.0;
                normalize(&angle);
            } else {
                propagate(&position, angle, ((surface.height[nsurf - 1] - p0) - position.z)/angle.z);
                angle.x -= (position.x/p0);
                angle.y -= (position.y/p0);
                angle.z = 1.0;
                normalize(&angle);
            }
            if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, 200);
        } else {
            newSurf = -1;
        }


        // OPTICS AND DETECTOR
        photon.direction = 1;
        photon.ghostFlag = 0;
        photon.saturationFlag = 0;
    surfaceloop: while (1) {
            oldSurf = newSurf;
            if (photon.direction == 1) {
                newSurf++;
            } else {
                newSurf--;
            }

        redostraylight:;

            // Find intercept
            if (newSurf>= 0 && newSurf<nsurf) {
                transform(&position, &angle, newSurf, 0, &photon);
                if (surface.surfacetype[newSurf] == DETECTOR) {
                    transform(&position, &angle, newSurf + 1, focusFlag, &photon);
                }
                miss = findSurface(angle, position, &distance, newSurf, &photon);
            } else {
                miss = 1;
            }

            //   Missed surface or wrong direction
            if (miss == 1 || (distance < 0)) {
                photon.ghostFlag = 1;
                if (straylightcurrent == 0) {
                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;
                } else {
                    if (chooseSurface(&newSurf, &oldSurf, &photon) == 1) {
                        countBad(&localLog, sourceOver*backSourceOver, &ray);
                       goto photonloop;
                    } else {
                        goto redostraylight;
                    }
                }
            }
            if (onlyadjacent == 1) {
                if ((newSurf > (oldSurf+1)) || (newSurf < (oldSurf-1))) {
                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;
                }
            }

            propagate(&position, angle, distance);
            if (opdfile) {

                if (newSurf == 0) {
                    photon.opdx = -position.z*angle.x/angle.z + position.x;
                    photon.opdy = -position.z*angle.y/angle.z + position.y;
                    distance = position.x*angle.x + position.y*angle.y + (position.z - 20000)*angle.z;
                }
                photon.op -= distance*photon.ncurr;
            }

            if (throughputfile == 1 && photon.direction == 1) {
                addThroughput(&state.throughputLog, minwavelength, maxwavelength, newSurf, waveIndex, sourceOver*backSourceOver);
            }

            if (photon.direction == -1) photon.ghostFlag = 1;

            //   DERIVATIVES
            interceptDerivatives(&normal, position, newSurf);

            //   CONTAMINATION
            if (surface.surfacetype[newSurf] !=  DETECTOR && contaminationmode == 1) {
                miss = contaminationSurfaceCheck(position, &angle, newSurf, &photon);
                if (miss) {
                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;
                }
 
            }



            //   SURFACE COATINGS
            transmission = surfaceCoating(photon.wavelength, angle, normal, newSurf, &reflection, &photon);

            if (transmissionCheck(transmission, natmospherefile*2 + 1 + newSurf*2, waveSurfaceIndex, &photon)) {
                if (straylightcurrent == 1 && ghost[newSurf] == 0) {
                    if (transmissionCheck(reflection + transmission, natmospherefile*2 + 1 + newSurf*2 + 1, waveSurfaceIndex, &photon)) {
                        countBad(&localLog, sourceOver*backSourceOver, &ray);
                        goto photonloop;
                    } else {
                        photon.direction = -photon.direction;
                        reflect(&angle, normal);
                        transformInverse(&position, &angle, newSurf);
                        if (surface.surfacetype[newSurf] == DETECTOR) {
                            transformInverse(&position, &angle, newSurf + 1);
                        }
                        if (surface.surfacetype[newSurf] == MIRROR) {
                            countBad(&localLog, sourceOver*backSourceOver, &ray);
                            goto photonloop;
                        }
                        goto surfaceloop;
                    }
                } else {
                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;
                }
            }
            //   INTERACTIONS
            if (surface.surfacetype[newSurf] == MIRROR) {

                //   MIRROR
                reflect(&angle, normal);
                transformInverse(&position, &angle, newSurf);
                if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);

            } else if (surface.surfacetype[newSurf] == LENS || surface.surfacetype[newSurf] == FILTER) {

                //   LENS/FILTER
                newRefractionIndex(newSurf, &photon);
                refract(&angle, normal, photon.nprev, photon.ncurr);
                transformInverse(&position, &angle, newSurf);
                if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);

            } else if (surface.surfacetype[newSurf] == GRATING) {

                //   GRATING
                double wavelengthNm = photon.wavelength*1000.0;
                Vector angleOut;
                pGrating->diffract(angle.x, angle.y, angle.z, normal.x, normal.y, normal.z,
                                   angleOut.x, angleOut.y, angleOut.z, wavelengthNm);
                vectorCopy(angleOut, &angle);
                transformInverse(&position, &angle, newSurf);
                if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);

            } else if (surface.surfacetype[newSurf] == DETECTOR) {

                if (eventfile || opdfile) {
                    transformInverse(&position, &angle, newSurf + 1);
                    transformInverse(&position, &angle, newSurf);
                    if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);
                    if (opdfile) {
                        if (ray == 0) {
                            state.epR[ssource] = -position.y*angle.z/angle.y;
                            // state.epR = 2730.1296937159;
                            printf("EPR: %.9f wavelength %.9f\n", state.epR[ssource], photon.wavelength);
                            //double epR[] = {2720.6523592541, 2730.1296937159,
                            //    2734.7054985315, 2738.3089059062,
                            //    2740.5318819465, 2741.7425507257};
                            //epR[1] = 2730.0677516303;
                            //epR[1] = 2730.1916350528;
                        } else if (ray == 1) {
                            state.cx[ssource] = position.x;
                            state.cy[ssource] = position.y;
                            state.cz[ssource] = position.z;
                            state.r0[ssource] = sqrt(pow(state.epR[ssource], 2) + pow(state.cx[ssource], 2) + pow(state.cy[ssource], 2));
                            printf("chief ray position: %.10f %.10f %.10f %.10f\n", state.cx[ssource], state.cy[ssource], state.cz[ssource], state.r0[ssource]);
                            pthread_mutex_unlock(&lock.lock4);
                            lock4 = 0;
                            pthread_mutex_lock(&lock.lock6);
                            remain--;
                            if (remain == 0) {
                                pthread_cond_broadcast(&lock.cond);
                            } else {
                                while (remain != 0) {
                                    pthread_cond_wait(&lock.cond,&lock.lock6);
                                }
                            }
                            pthread_mutex_unlock(&lock.lock6);
                            pthread_mutex_lock(&lock.lock7);
                            for (int surfIdx = 0; surfIdx < nsurf; surfIdx++) {
                                surface.innerRadius[surfIdx] = surface.innerRadius0[surfIdx];
                                surface.asphere(surfIdx, SURFACE_POINTS);
                            }
                            pthread_mutex_unlock(&lock.lock7);
                        }
                        //solve line-sphere intersection analytically
                        double ocx = position.x - state.cx[ssource];
                        double ocy = position.y - state.cy[ssource];
                        double ocz = position.z - state.cz[ssource];
                        double ocsqr = ocx*ocx + ocy*ocy + ocz*ocz;
                        double ocproj = ocx*angle.x + ocy*angle.y + ocz*angle.z;
                        distance = - ocproj + sqrt(ocproj*ocproj - ocsqr + state.r0[ssource]*state.r0[ssource]);
                        if (ray == 0) distance = state.epR[ssource];
                        photon.op -= distance;
                    }
                    transform(&position, &angle, newSurf, 0, &photon);
                    transform(&position, &angle, newSurf + 1, 0, &photon);
                }
                vectorCopy(position, &positionPrevious);
                vectorCopy(angle, &anglePrevious);
                detRay = 0;
                
            detectorloop: while (detRay < backSourceOver) {

                    position.x = positionPrevious.x + random[thread].normal()*backSigma;
                    position.y = positionPrevious.y + random[thread].normal()*backSigma;
                    position.z = positionPrevious.z;
                    vectorCopy(anglePrevious, &angle);

                    //   SILICON
                    if (detectorMode) {

                        photon.xPos = static_cast<long>(floor(position.x*1000/pixsize + pixelsx/2));
                        photon.yPos = static_cast<long>(floor(position.y*1000/pixsize + pixelsy/2));
                        if (sources.type[ssource] < 6) {
                            if (photon.xPos <= minx - activeBuffer || photon.xPos >= maxx + activeBuffer ||
                                photon.yPos <= miny - activeBuffer || photon.yPos >= maxy + activeBuffer) {
                                if (saturation) saturate(ssource, &largeAngle, &photon, shiftX, shiftY);
                                countBad(&localLog, sourceOver, &ray);
                                detRay++;
                                goto detectorloop;
                            }
                        }

                        double dh;
                        miss = getDeltaIntercept(position.x, position.y, &dh, newSurf, &photon);

                        if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny &&
                            photon.yPos <= maxy && contaminationmode == 1) {
                            if (random[thread].uniform()>(double)(*(contamination.chiptransmission +
                                                            chip.nampx*(photon.yPos - miny) + (photon.xPos - minx)))) {
                                countBad(&localLog, sourceOver, &ray);
                                detRay++;
                                goto detectorloop;
                            }
                            if (*(contamination.chiplistmap + chip.nampx*(photon.yPos - miny) + (photon.xPos - minx)) != -1) {
                                double xx = ((position.x*1000/pixsize) + pixelsx/2)*1e-3*pixsize;
                                double yy = ((position.y*1000/pixsize) + pixelsy/2)*1e-3*pixsize;
                                long cc = *(contamination.chiplistmap + chip.nampx*(photon.yPos - miny) + (photon.xPos - minx));
                                if (sqrt(pow(xx - contamination.chiplistx[cc], 2.0) +
                                         pow(yy - contamination.chiplisty[cc], 2.0)) <
                                    contamination.chiplists[cc]) {
                                    if (random[thread].uniform() > exp(-contamination.absorptionLength*contamination.chiplists[cc])) {
                                        countBad(&localLog, sourceOver, &ray);
                                        detRay++;
                                        goto detectorloop;
                                    }
                                    long index;
                                    find(contamination.henyey_greenstein, contamination.elements, random[thread].uniform(), &index);
                                    double mu = contamination.henyey_greenstein_mu[index];
                                    double phi = random[thread].uniform()*2*PI;
                                    shift_mu(&angle, mu, phi);
                                    if (mu < 0) {
                                        countBad(&localLog, sourceOver, &ray);
                                        detRay++;
                                        goto detectorloop;
                                    }
                                    photon.ghostFlag=1;
                                }
                            }
                        }


                        miss = photonSiliconPropagate(&angle, &position, photon.wavelength, normal, dh, waveSurfaceIndex, &photon);

                        if (miss == 0) {
                            if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, 300);
                        } else {
                            if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, 304);
                            countBad(&localLog, sourceOver, &ray);
                            detRay++;
                            goto detectorloop;
                        }

                        int notFinished = 1;
                        int eCounter = 0;
                        while (notFinished == 1) {
                            photon.xPos = static_cast<long>(floor(position.x*1000/pixsize + pixelsx/2));
                            photon.yPos = static_cast<long>(floor(position.y*1000/pixsize + pixelsy/2));
                            photon.xPosR = position.x*1000/pixsize - floor(position.x*1000/pixsize) - 0.5;
                            photon.yPosR = position.y*1000/pixsize - floor(position.y*1000/pixsize) - 0.5;

                            miss = electronSiliconPropagate(&angle, &position, &photon);
                            if (miss == 1) notFinished = 0;
                            if (photon.z0 > photon.collect_z) {
                                if (position.z <= photon.collect_z) {
                                    notFinished = 0;
                                    position.z = photon.collect_z;
                                }
                            }
                            if (photon.z0 <= photon.collect_z) {
                                if (position.z >= photon.collect_z) {
                                    notFinished = 0;
                                    position.z = photon.collect_z;
                                }
                            }
                            eCounter++;
                            if (eCounter > (static_cast<double>(SILICON_STEPS)/static_cast<double>(SILICON_SUB_STEPS))) notFinished = 0;
                            if (eventfile == 1 && notFinished == 0) state.pEventLogging->logPhoton(position.x, position.y, position.z, 302);
                            if (eventfile == 1 && notFinished == 1) state.pEventLogging->logPhoton(position.x, position.y, position.z, 301);

                        }

                    }

                    photon.xPos = (long)(floor(position.x*1000/pixsize + pixelsx/2));
                    photon.yPos = (long)(floor(position.y*1000/pixsize + pixelsy/2));

                    if (eventfile) {
                        if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {
                            state.pEventLogging->logPhoton(static_cast<double>(photon.xPos),
                                                           static_cast<double>(photon.yPos), 0.0, 303);
                        }
                    }

                    if (centroidfile) {
                        if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {
                            sourceXpos[ssource] += photon.xPos*sourceOver;
                            sourceYpos[ssource] += photon.yPos*sourceOver;
                            sourcePhoton[ssource] += sourceOver;
                        }
                    }


                    if (opdfile && ray > 0) {
                        // if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {
                        long xx = floor(photon.opdx/maxr/2*(OPD_SCREEN_SIZE - 1.0) + OPD_SCREEN_SIZE/2.0);
                        long yy = floor(photon.opdy/maxr/2*(OPD_SCREEN_SIZE - 1.0) + OPD_SCREEN_SIZE/2.0);
                            if (xx >= 0 && xx < OPD_SCREEN_SIZE && yy >= 0 && yy < OPD_SCREEN_SIZE) {
                                pthread_mutex_lock(&lock.lock5);
                                *(state.opd + ssource*OPD_SCREEN_SIZE*OPD_SCREEN_SIZE + OPD_SCREEN_SIZE*yy + xx) += photon.op;
                                *(state.opdcount + ssource*OPD_SCREEN_SIZE*OPD_SCREEN_SIZE + OPD_SCREEN_SIZE*yy + xx) += 1;
                                pthread_mutex_unlock(&lock.lock5);
                            }
                        }
                    // }

                    if (sources.type[ssource] < 6 && backGamma > 1.0) {
                        Vector newPosition;
                        vectorCopy(position, &newPosition);
                        long long lmax;
                        lmax = photon.sourceOver_m;
                        photon.sourceOver_m = 1;
                        for (long long l = 0; l < lmax; l++) {
                            if (l > 0) {
                                position.x = newPosition.x + random[thread].normal()*backSigma/backDelta;
                                position.y = newPosition.y + random[thread].normal()*backSigma/backDelta;
                            }
                            photon.xPos = static_cast<long>(floor(position.x*1000/pixsize + pixelsx/2));
                            photon.yPos = static_cast<long>(floor(position.y*1000/pixsize + pixelsy/2));
                              if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {

                                // pthread_mutex_lock(&lock.lock1);
                                if (saturation) {
                                    saturate(ssource, &largeAngle, &photon, shiftX, shiftY);
                                } else {
                                    *(state.focal_plane + chip.nampx*(photon.yPos - miny) +
                                      (photon.xPos - minx)) += photon.sourceOver_m;
                                }
                                // pthread_mutex_unlock(&lock.lock1);
                                countGood(&localLog, photon.sourceOver_m, &ray);
                            } else {

                                countBad(&localLog, photon.sourceOver_m, &ray);
                            }
                        }

                        photon.sourceOver_m = lmax;
                    } else {
                        if (photon.xPos >= minx && photon.xPos <= maxx && photon.yPos >= miny && photon.yPos <= maxy) {
                            // pthread_mutex_lock(&lock.lock1);
                            if (saturation) {
                                saturate(ssource, &largeAngle, &photon, shiftX, shiftY);
                            } else {
                                *(state.focal_plane + chip.nampx*(photon.yPos - miny) + (photon.xPos - minx)) += photon.sourceOver_m;
                            }
                           // pthread_mutex_unlock(&lock.lock1);

                        } else {
                            if (saturation) saturate(ssource, &largeAngle, &photon, shiftX, shiftY);
                            // countBad(&localLog, sourceOver, &ray);
                            countBad(&localLog, photon.sourceOver_m, &ray);
                            detRay++;
                            goto detectorloop;
                        }
                    }

                    if (throughputfile) addThroughput(&state.throughputLog, minwavelength, maxwavelength, nsurf, waveIndex, photon.sourceOver_m);
                    detRay++;
                    if (sources.type[ssource] >= 6 || (sources.type[ssource] < 6 && backGamma <= 1.0)) {
                        countGood(&localLog, photon.sourceOver_m, &ray);
                    }
                }
                break;

            }
        }
    }

    // pthread_mutex_lock(&lock.lock2);
    counterAdd(&localLog, &state.counterLog);
    counterAdd(&localLog, &state.globalLog);
    // pthread_mutex_unlock(&lock.lock2);

    if (finish != 0) {
        pthread_mutex_lock(&lock.lock8);
        openthreads -= finish;
        openthread[thread] = 0;
        pthread_mutex_unlock(&lock.lock8);
    }
    if (lock4 == 1) pthread_mutex_unlock(&lock.lock4);

    // pthread_exit(NULL);

}
