///
/// @package phosim
/// @file photonmanipulate.cpp
/// @brief photon manipulation routines (part of Image class)
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

int Image::domeSeeing(Vector *angle, Photon *aph) {

    double phi, r;

    phi = 2*PI*random[aph->thread].uniform();
    r = sqrt(domeseeing*domeseeing + toypsf*toypsf)*ARCSEC/2.35482*random[aph->thread].normal();

    angle->x = angle->x + r*cos(phi);
    angle->y = angle->y + r*sin(phi);
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    return(0);

}

int Image::tracking(Vector *angle, double time) {

    double rindex, vxp, vyp;
    long index;

    index = find_linear(perturbation.jittertime, trackinglines, time, &rindex);
    vxp = angle->x - (cos(spiderangle)*perturbation.jitterele[index] +
                      sin(spiderangle)*perturbation.jitterazi[index])*ARCSEC;
    vyp = angle->y + ((-1.0)*(sin(spiderangle)*perturbation.jitterele[index]) +
                      cos(spiderangle)*perturbation.jitterazi[index])*ARCSEC;
    angle->y = vyp*cos(perturbation.jitterrot[index]*ARCSEC) + (-1.0)*vxp*sin(perturbation.jitterrot[index]*ARCSEC);
    angle->x = vyp*sin(perturbation.jitterrot[index]*ARCSEC) + (vxp)*cos(perturbation.jitterrot[index]*ARCSEC);
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    //wind shake
    double wx = perturbation.windshake[index]*ARCSEC*cos(winddir[7]*DEGREE - azimuth);
    double wy = perturbation.windshake[index]*ARCSEC*sin(winddir[7]*DEGREE - azimuth);
    // printf("%e %e %e %e\n",perturbation.windshake[index],screen.jitterwind[index],wx - wy*screen.jitterwind[index]*DEGREE,wy + wx*screen.jitterwind[index]*DEGREE);

    angle->x += wx - wy*screen.jitterwind[index]*DEGREE;
    angle->y += wy + wx*screen.jitterwind[index]*DEGREE;
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    return(0);
}

int Image::largeAngleScattering(Vector *largeAngle, Photon *aph) {

    long index;
    double r, phi;

    r = random[aph->thread].uniform();
    if (r < laScatterProb) {
        r = random[aph->thread].uniform();
        find(perturbation.miescatter, 10000, r, &index);
        r = ((double)(index))/10000.0*0.1*DEGREE;
        phi = 2.0*PI*random[aph->thread].uniform();
        largeAngle->x += r*cos(phi);
        largeAngle->y += r*sin(phi);
    }

    return(0);

}

int Image::secondKick(Vector *largeAngle, Photon *aph) {

    long index;
    double r, phi;
    if (pupilscreenMode == 1 && atmospheremode < 2) {
        int i;
        int j;
        r = random[aph->thread].uniform();
        // 2D sample from normalized cumulative FFT distribution
        for (i = 0; i < SCREEN_SIZE; i++) {
            if (screen.focalscreencum[i*SCREEN_SIZE] > r) {
                i--;
                goto jumpi;
            }
        }
        i--;
        jumpi:;
        for (j = 0; j < SCREEN_SIZE; j++) {
            if (screen.focalscreencum[i*SCREEN_SIZE + j] > r) {
                j--;
                goto jumpj;
            }
        }
        j--;
        jumpj:;
        i -= SCREEN_SIZE/2;
        j -= SCREEN_SIZE/2;
        r = (static_cast<double>(sqrt(i*i + j*j)))*screen.pupilscreenscale*1e-3/
            (SCREEN_SIZE*screen.fine_sizeperpixel*screen.paddingfactor)*aph->wavelength;
        phi = static_cast<double>(atan2(j, i));
    }
    else {
        r = random[aph->thread].uniform();
        find(screen.hffunc, 10000, r, &index);
        // If turbulence off, do pupil diffraction. Otherwise, do Kolmogorov diffraction.
        if (atmospheremode < 2) {
            r = ((double)(index))*0.5*1e-3/(surface.outerRadius[0]*screen.paddingfactor)*aph->wavelength;
        }
        else
            r = ((double)(index))*0.5*1e-3/(SCREEN_SIZE*screen.fine_sizeperpixel)*aph->wavelengthFactor;
        phi = 2*PI*random[aph->thread].uniform();
    }
    largeAngle->x += r*cos(phi);
    largeAngle->y += r*sin(phi);

    return(0);

}


int Image::diffraction(Vector *position, Vector angle, Vector *largeAngle, Photon *aph) {

    double distance;
    double xl, yl, zorig, r;
    int i;
    double diffdist, adiffdist;
    int signv = 0;
    double diffx = 0.0, diffy = 0.0, mindiffdist;

    mindiffdist = 1e30;
    zorig = position->z;
    double cosShiftedAngle = cos(aph->shiftedAngle);
    double sinShiftedAngle = sin(aph->shiftedAngle);
    double mindist = 1e-3;

    for (i = 0; i < obstruction.nspid; i++) {

        distance = (obstruction.height[i] - (position->z))/angle.z;
        propagate(position, angle, distance);
        xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
        yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;

        if (obstruction.type[i] == 1) {

            if (obstruction.angle[i] != 0) {
                distance = (obstruction.height[i] +
                            (fabs(yl) - obstruction.reference[i])*sin(DEGREE*obstruction.angle[i]) -
                            (position->z))/angle.z;
                propagate(position, angle, distance);
                xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
                yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;
            }

            diffdist = xl - obstruction.center[i];
            if (diffdist > 0) {
                signv = 1;
            } else {
                signv = -1;
            }
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.width[i] + mindist) return(1);
            if (adiffdist - obstruction.width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.width[i];
                diffx = cosShiftedAngle*diffdist/adiffdist;
                diffy = -sinShiftedAngle*diffdist/adiffdist;
            }

            distance = (obstruction.height[i] - obstruction.depth[i] +
                        (fabs(yl) - obstruction.reference[i])*sin(DEGREE*obstruction.angle[i]) -
                        (position->z))/angle.z;
            propagate(position, angle, distance);
            xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
            yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;

            diffdist = xl - obstruction.center[i];
            if (diffdist > 0 && signv == -1) return(1);
            if (diffdist < 0 && signv == 1) return(1);
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.width[i] + mindist) return(1);
            if (adiffdist - obstruction.width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.width[i];
                diffx = cosShiftedAngle*diffdist/adiffdist;
                diffy = -sinShiftedAngle*diffdist/adiffdist;
            }

        }

        if (obstruction.type[i] == 2) {

            if (obstruction.angle[i] != 0) {
                distance = (obstruction.height[i] +
                            (fabs(xl)-obstruction.reference[i])*sin(DEGREE*obstruction.angle[i]) -
                            (position->z))/angle.z;
                propagate(position, angle, distance);
                xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
                yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;
            }

            diffdist = yl - obstruction.center[i];
            if (diffdist > 0) {
                signv = 1;
            } else {
                signv = -1;
            }
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.width[i] + mindist) return(1);
            if (adiffdist - obstruction.width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.width[i];
                diffx = sinShiftedAngle*diffdist/adiffdist;
                diffy = cosShiftedAngle*diffdist/adiffdist;
            }

            distance = (obstruction.height[i] - obstruction.depth[i] +
                        (fabs(xl)-obstruction.reference[i])*sin(DEGREE*obstruction.angle[i]) -
                        (position->z))/angle.z;
            propagate(position, angle, distance);
            xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
            yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;

            diffdist = yl - obstruction.center[i];
            if (diffdist > 0 && signv == -1) return(1);
            if (diffdist < 0 && signv == 1) return(1);
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.width[i] + mindist) return(1);
            if (adiffdist - obstruction.width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.width[i];
                diffx = sinShiftedAngle*diffdist/adiffdist;
                diffy = cosShiftedAngle*diffdist/adiffdist;
            }
        }
    }

    distance = (zorig - (position->z))/angle.z;
    propagate(position, angle, distance);

    if (obstruction.pupil == 1) {
        r = sqrt((position->x)*(position->x) + (position->y)*(position->y));
        diffdist = surface.outerRadius[0] - r;
        if (diffdist < mindist) return(1);
        if (diffdist < mindiffdist) {
            mindiffdist = diffdist;
            diffx = (position->x)/r;
            diffy = (position->y)/r;
        }
        diffdist = r - surface.innerRadius[0];
        if (diffdist < mindist) return(1);
        if (diffdist < mindiffdist) {
            mindiffdist = diffdist;
            diffx = (position->x)/r;
            diffy = (position->y)/r;
        }
    }

    if (mindiffdist < surface.outerRadius[0] - surface.innerRadius[0]) {
        r = (aph->wavelength/1000.0)/(4*PI*mindiffdist);
        if (diffractionMode != 5) {
            largeAngle->x += r*diffx;
            largeAngle->y += r*diffy;
        }
    }
    return(0);

}

double Image::airIndexRefraction(Photon *aph, long layer) {

    if (airrefraction) {
        if (layer == -1) {
            aph->airRefractionPreviousADC=0.0;
            aph->airRefractionADC = 0.0;
            return(0.0);
        } else {

            double currentheight=height[layer]+groundlevel/1e3;
            double temp=(temperature+273.15)*(1-currentheight/44329.)/(1-groundlevel/1e3/44329.);
            double press=pressure*pow(1-currentheight/44329.,5.255876)/pow(1-groundlevel/1e3/44329.,5.255876);
            // double temp=(temperature+273.15);

            double sigma=1.0/aph->wavelength;
            double ps=press/760.00*1013.25;
            double pw=waterPressure/760.00*1013.25;
            double dw=(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/temp-710.792/temp/temp+7.75141e4/temp/temp/temp))*pw/temp;
            double ds=(1+ps*(57.90e-8-9.325e-4/temp+0.25844/temp/temp))*ps/temp;
            double n=(2371.34+683939.7/(130.0-pow(sigma,2))+4547.3/(38.9-pow(sigma,2)))*ds;
            n=n+(6478.31-58.058*pow(sigma,2)-0.71150*pow(sigma,4)+0.08851*pow(sigma,6))*dw;
            n = 1e-8*n;

            aph->airRefractionPreviousADC=aph->airRefractionADC;
            sigma=1.0/centralwavelength;
            ps=press/760.00*1013.25;
            pw=waterPressure/760.00*1013.25;
            dw=(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/temp-710.792/temp/temp+7.75141e4/temp/temp/temp))*pw/temp;
            ds=(1+ps*(57.90e-8-9.325e-4/temp+0.25844/temp/temp))*ps/temp;
            aph->airRefractionADC=(2371.34+683939.7/(130.0-pow(sigma,2))+4547.3/(38.9-pow(sigma,2)))*ds;
            aph->airRefractionADC=aph->airRefractionADC+(6478.31-58.058*pow(sigma,2)-0.71150*pow(sigma,4)+0.08851*pow(sigma,6))*dw;
            aph->airRefractionADC = 1e-8*aph->airRefractionADC;

            return(n);
        }

        // double airRefraction = 64.328 + 29498.1/(146 - 1/aph->wavelength/aph->wavelength) + 255.4/(41 - 1/aph->wavelength/aph->wavelength);
        // airRefraction = airRefraction*pressure*(1 + (1.049 - 0.0157*temperature)*1e-6*pressure)/720.883/(1 + 0.003661*temperature);
        // airRefraction = airRefraction - ((0.0624 - 0.000680/aph->wavelength/aph->wavelength)/(1 + 0.003661*temperature)*waterPressure);
        // airRefraction = airRefraction/1e6;
        // return(airRefraction);
    } else {
        return(0.0);
    }

}

int Image::atmosphericDispersion (Vector *angle, Photon *aph, long layer) {

    double dx, dy, phi, theta, psi;
    //aph->shiftedAngle=0.0;
    double signa=1.0;
    double signb=1.0;

    if (atmospheric_dispersion) {
    if (layer == -1) {
        phi = aph->shiftedAngle*signa;
        theta = zenith;
        psi = -aph->shiftedAngle*signb;
        dx = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*(-angle->x) +
            (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*(-angle->y) +
            (sin(psi)*sin(theta))*fabs(angle->z);
        dy = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*(-angle->x) +
            (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*(-angle->y) +
            (cos(psi)*sin(theta))*fabs(angle->z);
        aph->refractionInvariant=(height[0]+groundlevel/1e3+RADIUS_EARTH)*(sqrt(dx*dx + dy*dy));
        if (atmosphericdispcenter) {
            phi = aph->shiftedAngle*signa;
            theta = zenith;
            psi = -aph->shiftedAngle*signb;
            dx = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*(0.0) +
                (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*(0.0) +
                (sin(psi)*sin(theta))*fabs(1.0);
            dy = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*(0.0) +
                (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*(0.0) +
                (cos(psi)*sin(theta))*fabs(1.0);
            aph->refractionInvariantCenter=(height[0]+groundlevel/1e3+RADIUS_EARTH)*(sqrt(dx*dx + dy*dy));
            aph->adcx = 0.0;
            aph->adcy = 0.0;
            aph->adcz = smallAnglePupilNormalize(aph->adcx, aph->adcy);
        }
    } else {

        double currentheight=height[layer]+groundlevel/1e3+RADIUS_EARTH;
        double correctionCenter=0.0;
        if (atmosphericdispcenter) {
        phi = aph->shiftedAngle*signa;
        theta = zenith;
        psi = -aph->shiftedAngle*signb;
        dx = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*(0.0) +
            (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*(0.0) +
            (sin(psi)*sin(theta))*fabs(1.0);
        dy = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*(0.0) +
            (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*(0.0) +
            (cos(psi)*sin(theta))*fabs(1.0);
        correctionCenter=aph->refractionInvariantCenter/sqrt(pow(currentheight*(1+aph->airRefractionADC),2.0) -
                                                             pow(aph->refractionInvariantCenter,2.0))*
            (aph->airRefractionADC - aph->airRefractionPreviousADC)/(1 + aph->airRefractionADC);
            if ((dx == 0.0) && (dy == 0.0)) {
                aph->adcx = correctionCenter;
                aph->adcy = 0.0;
            } else {
                aph->adcx = correctionCenter*dx/sqrt(dx*dx + dy*dy);
                aph->adcy = correctionCenter*dy/sqrt(dx*dx + dy*dy);
            }
            aph->adcz = smallAnglePupilNormalize(aph->adcx, aph->adcy);
        }
        double correction=aph->refractionInvariant/sqrt(pow(currentheight*(1+aph->airRefraction),2.0) -
                                                        pow(aph->refractionInvariant,2.0))*
            (aph->airRefraction - aph->airRefractionPrevious)/(1 + aph->airRefraction);
        phi = aph->shiftedAngle*signa;
        theta = zenith;
        psi = -aph->shiftedAngle*signb;
        dx = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*(-angle->x) +
            (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*(-angle->y) +
            (sin(psi)*sin(theta))*fabs(angle->z);
        dy = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*(-angle->x) +
            (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*(-angle->y) +
            (cos(psi)*sin(theta))*fabs(angle->z);
        if ((dx == 0.0) && (dy == 0.0)) {
            angle->x += correction;
        } else {
            angle->x += correction*dx/sqrt(dx*dx + dy*dy);
            angle->y += correction*dy/sqrt(dx*dx + dy*dy);
        }
        // printf("a %ld %e\n",layer,(correction-correctionCenter)/ARCSEC);
        if (atmosphericdispcenter) {
            angle->x -= aph->adcx;
            angle->y -= aph->adcy;
        }
        // printf("b %ld %e %e %e\n",layer,correctionCenter,aph->refractionInvariantCenter,aph->airRefraction);
        // printf("%e %e %e %e\n",correction,correctionCenter,aph->airRefraction,aph->airRefractionADC);
        angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    }
    }
    // if (layer == 0) {

    //     if (atmosphericdispcenter) {
    //         dx = zenith*sin(aph->shiftedAngle);
    //         dy = zenith*cos(aph->shiftedAngle);
    //         adcx = tan(sqrt(dx*dx + dy*dy))*dx/sqrt(dx*dx + dy*dy)*air.air_refraction_adc;
    //         adcy = tan(sqrt(dx*dx + dy*dy))*dy/sqrt(dx*dx + dy*dy)*air.air_refraction_adc;
    //         if (zenith == 0.0) {
    //             adcx = 0.0;
    //             adcy = 0.0;
    //         }
    //     }
    //     if (atmospheric_dispersion) {
    //         dx = -angle->x + zenith*sin(aph->shiftedAngle);
    //         dy = -angle->y + zenith*cos(aph->shiftedAngle);
    //         angle->x = angle->x + tan(sqrt(dx*dx + dy*dy))*dx/sqrt(dx*dx + dy*dy)*aph->airRefraction;
    //         angle->y = angle->y + tan(sqrt(dx*dx + dy*dy))*dy/sqrt(dx*dx + dy*dy)*aph->airRefraction;
    //         if (atmosphericdispcenter) {
    //             angle->x = angle->x - adcx;
    //             angle->y = angle->y - adcy;
    //         }
    //         angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    //     }

    // }
    return(0);

}


int Image::samplePupil (Vector *position, long long ray, Photon *aph) {

    double r, phi; // mm, radians

    if (aperturemode == 1) {
        double x, y;
        x = ((static_cast<double>((ray % 40000) % 200))/200.0)*2.0*maxr - maxr;
        y = ((static_cast<double>((ray % 40000) / 200))/200.0)*2.0*maxr - maxr;
        r = sqrt(x*x + y*y);
        phi = atan2(y, x);
        if (r < minr || r > maxr) return(1);
    } else if (aperturemode == 2) {
        if (ray < 2) {
            r = 1e-14;
            phi = 0.0;
        } else {
            double x, y;
            long opdsize2 = opdsize*opdsize*opdsampling*opdsampling;
            x = ((static_cast<double>((ray % opdsize2) % (opdsampling*opdsize)))/(opdsampling*opdsize - 1))*2.0*maxr - maxr;
            y = ((static_cast<double>((ray % opdsize2) / (opdsampling*opdsize)))/(opdsampling*opdsize - 1))*2.0*maxr - maxr;
            r = sqrt(x*x + y*y);
            phi = atan2(y, x);
            if (r < minr || r > maxr) return(1);
        }
    } else {
        r = sqrt(random[aph->thread].uniform()*(maxr*maxr - minr*minr) + minr*minr);
        phi = random[aph->thread].uniform()*2*PI;
    }

    int surfaceIndex = 0; // start at first surface in optics file

    double R = surface.radiusCurvature[surfaceIndex]; // radius of curvature (mm)
    double asphere; // sum of asphere coefficients (mm)
    double third = -surface.three[surfaceIndex]*1e3; // 3rd order asphere coefficient (mm)
    double fourth = -surface.four[surfaceIndex]*1e3; // 4th order asphere coefficient (mm)
    double fifth = -surface.five[surfaceIndex]*1e3; // 5th order asphere coefficient (mm)
    double sixth = -surface.six[surfaceIndex]*1e3; // 6th order asphere coefficient (mm)
    double seventh = -surface.seven[surfaceIndex]*1e3; // 7th order asphere coefficient (mm)
    double eighth = -surface.eight[surfaceIndex]*1e3; // 8th order asphere coefficient (mm)
    double ninth = -surface.nine[surfaceIndex]*1e3; // 9th order asphere coefficient (mm)
    double tenth = -surface.ten[surfaceIndex]*1e3; // 10th order asphere coefficient (mm)

    asphere = third*pow(r, 3.0) + fourth*pow(r, 4.0) + fifth*pow(r, 5.0) + sixth*pow(r, 6.0) + seventh*pow(r, 7.0) + eighth*pow(r, 8.0) + ninth*pow(r,9.0) + tenth*pow(r, 10.0);

    position->x = r*cos(phi);
    position->y = r*sin(phi);
    if (R != 0.0) {
        double k = surface.conic[surfaceIndex]; // conic constant (unitless)
        // sagitta equation
        position->z = surface.height[surfaceIndex] + r*r/(R*(1.0 + sqrt(1.0 - (1.0 + k)*r*r/(R*R)))) + asphere;
    } else {
        position->z = surface.height[surfaceIndex] + asphere;
    }

    return(0);

}


int Image::transmissionCheck(double transmission, long surfaceIndex, long waveSurfaceIndex, Photon *aph) {

    double randNum;

    (aph->counter)++;
    if (transmission > state.dynamicTransmission[waveSurfaceIndex + surfaceIndex]) {
        state.dynamicTransmission[waveSurfaceIndex + surfaceIndex] = transmission;
    }
    if (transmission < state.dynamicTransmissionLow[waveSurfaceIndex + surfaceIndex]) {
        state.dynamicTransmissionLow[waveSurfaceIndex + surfaceIndex] = transmission;
    }
    if (aph->counter <= aph->maxcounter) {
        randNum = aph->saveRand[aph->counter];
    } else {
        randNum = random[aph->thread].uniform();
    }
    if (randNum > transmission) {
        return(1);
    } else {
        return(0);
    }

}

int Image::transmissionPreCheck(long surfaceIndex, long waveSurfaceIndex, Photon *aph) {

    (aph->counter)++;
    aph->saveRand[aph->counter] = random[aph->thread].uniform();
    if (aph->saveRand[aph->counter] > (fabs(state.dynamicTransmission[waveSurfaceIndex + surfaceIndex]) + transtol)) {
        return(1);
    } else {
        if (aph->saveRand[aph->counter] > (fabs(state.dynamicTransmissionLow[waveSurfaceIndex + surfaceIndex]) - transtol)) {
            // double f1=state.dynamicTransmission[waveSurfaceIndex + surfaceIndex];
            // double f2=state.dynamicTransmissionLow[waveSurfaceIndex + surfaceIndex];
            // double f3=aph->saveRand[aph->counter];
            // printf("%ld %ld %lf %lf %lf\n",surfaceIndex,waveSurfaceIndex,f1,f2,f3);
            return(2);
        } else {
            return(0);
        }
    }

}

double Image::surfaceCoating (double wavelength, Vector angle, Vector normal, long newSurf, double *reflection, Photon *aph) {

    double filterAngle, rindex, crindex;
    long index, cindex;

    if (surface.surfacecoating[newSurf] != 0 && coatingmode == 1) {
        if (coating.angleNumber[newSurf] > 1) {
            Vector tempangle;
            vectorCopy(angle, &tempangle);
            refract(&tempangle, normal, aph->ncurr, 1.0);
            double arg = fabs(normal.x*tempangle.x + normal.y*tempangle.y + normal.z*tempangle.z);
            if (arg > 1) arg = 1.0;
            filterAngle = acos(arg)/DEGREE;
            index = find_linear(coating.wavelength[newSurf], coating.wavelengthNumber[newSurf], wavelength, &rindex);
            cindex = find_linear(coating.angle[newSurf], coating.angleNumber[newSurf], filterAngle, &crindex);
            if (crindex - cindex < 0.0) crindex = static_cast<double>(cindex);
            if (crindex - cindex > 1.0) crindex = static_cast<double>(cindex + 1);
            if (rindex - index < 0.0) rindex = static_cast<double>(index);
            if (rindex - index > 1.0) rindex = static_cast<double>(index + 1);
            *reflection = interpolate_bilinear(coating.reflection[newSurf], coating.wavelengthNumber[newSurf], cindex, crindex, index, rindex);
            return(interpolate_bilinear(coating.transmission[newSurf], coating.wavelengthNumber[newSurf], cindex, crindex, index, rindex));
        } else {
            index = find_linear(coating.wavelength[newSurf], coating.wavelengthNumber[newSurf], wavelength, &rindex);
            if (rindex - index < 0.0) rindex = static_cast<double>(index);
            if (rindex - index > 1.0) rindex = static_cast<double>(index + 1);
            *reflection = interpolate_linear(coating.reflection[newSurf], index, rindex);
            return(interpolate_linear(coating.transmission[newSurf], index, rindex));
        }
    } else {
        return(1.0);
    }

}

void Image::newRefractionIndex(long surfaceIndex, Photon *aph) {

    long index, newMedium;
    double rindex;

    aph->nprev = aph->ncurr;
    if (aph->direction == 1) {
        newMedium = surfaceIndex;
    } else {
        newMedium = surfaceIndex - 1;
    }
    if (newMedium >= 0) {
        if (surface.surfacemed[newMedium] == 0) {
            aph->ncurr = 1.0;
        } else if (surface.surfacemed[newMedium] == 2) {
            aph->ncurr = 1.0 + aph->airRefraction;
        } else {
            index = find_linear(medium.indexRefractionWavelength[newMedium], medium.indexRefractionNumber[newMedium], aph->wavelength, &rindex);
            aph->ncurr = interpolate_linear(medium.indexRefraction[newMedium], index, rindex);
        }
    } else {
        aph->ncurr = 1.0 + aph->airRefraction;
    }

}


void Image::atmospherePropagate(Vector *position, Vector angle, long layer, int mode, Photon *aph) {

    if (layer == -1) {

        propagate(position, angle, (1e6*70.0 - position->z)/(angle.z));
        aph->xporig = position->x;
        aph->yporig = position->y;
        aph->zporig = position->z;

        if (mode == 2) {
            if (fabs(aph->time - aph->prtime) > screentol) {
                for (long k = 0; k < SCREEN_SIZE; k++) {
                    for (long l = 0; l < SCREEN_SIZE; l++) {
                        *(screen.phasescreen + k*SCREEN_SIZE + l) = 0;
                    }
                }
            }
        }

    } else {

        double distance = (1e6*height[layer] - position->z)/(angle.z);
        propagate(position, angle, distance);

    }

}

void Image::atmosphereIntercept(Vector *position, long layer, Photon *aph) {

    double rindex;
    long index;
    double wx = wind[layer]*1.0e3*(aph->absoluteTime)*cos(winddir[layer]*DEGREE - azimuth);
    double wy = wind[layer]*1.0e3*(aph->absoluteTime)*sin(winddir[layer]*DEGREE - azimuth);


    // wind blur
    index = find_linear(perturbation.jittertime, trackinglines, aph->absoluteTime, &rindex);

    aph->windx = wx - wy*screen.jitterwind[index]*DEGREE;
    aph->windy = wy + wx*screen.jitterwind[index]*DEGREE;

    aph->lindex = SCREEN_SIZE*SCREEN_SIZE*layer;

    aph->xpos = position->x + xtelloc + aph->windx;
    aph->ypos = position->y + ytelloc + aph->windy;

    find_linear_wrap(aph->xpos, screen.large_sizeperpixel, SCREEN_SIZE, &(aph->indexlx0), &(aph->indexlx1), &(aph->dlx));
    find_linear_wrap(aph->ypos, screen.large_sizeperpixel, SCREEN_SIZE, &(aph->indexly0), &(aph->indexly1), &(aph->dly));
    find_linear_wrap(aph->xpos, screen.coarse_sizeperpixel, SCREEN_SIZE, &(aph->indexcx0), &(aph->indexcx1), &(aph->dcx));
    find_linear_wrap(aph->ypos, screen.coarse_sizeperpixel, SCREEN_SIZE, &(aph->indexcy0), &(aph->indexcy1), &(aph->dcy));
    find_linear_wrap(aph->xpos, screen.medium_sizeperpixel, SCREEN_SIZE, &(aph->indexmx0), &(aph->indexmx1), &(aph->dmx));
    find_linear_wrap(aph->ypos, screen.medium_sizeperpixel, SCREEN_SIZE, &(aph->indexmy0), &(aph->indexmy1), &(aph->dmy));
    find_linear_wrap(aph->xpos, screen.fine_sizeperpixel, SCREEN_SIZE, &(aph->indexfx0), &(aph->indexfx1), &(aph->dfx));
    find_linear_wrap(aph->ypos, screen.fine_sizeperpixel, SCREEN_SIZE, &(aph->indexfy0), &(aph->indexfy1), &(aph->dfy));


}

void Image::atmosphereRefraction (Vector *angle, long layer, int mode, Photon *aph) {

    double scaleOuter;

    if (mode == 1 || mode == 5) {
        scaleOuter = aph->wavelengthFactor*ARCSEC*screen.secondKickSize;
    } else {
        scaleOuter = aph->wavelengthFactor*ARCSEC;
    }

    if (mode <= 1 || mode >= 4) {

        if (atmdebug == 0) {

            (angle->x) += (interpolate_bilinear_float_wrap(screen.turbulenceLargeX + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                         aph->dlx, aph->indexly0, aph->indexly1, aph->dly) +
                         interpolate_bilinear_float_wrap(screen.turbulenceCoarseX + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                         aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy) +
                         interpolate_bilinear_float_wrap(screen.turbulenceMediumX + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                         aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy))*scaleOuter;

            (angle->y) += (interpolate_bilinear_float_wrap(screen.turbulenceLargeY + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                         aph->dlx, aph->indexly0, aph->indexly1, aph->dly) +
                         interpolate_bilinear_float_wrap(screen.turbulenceCoarseY + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                         aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy) +
                         interpolate_bilinear_float_wrap(screen.turbulenceMediumY + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                         aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy))*scaleOuter;

            angle->z = smallAnglePupilNormalize(angle->x, angle->y);

        } else {

            aph->indexlx0 = static_cast<long>(floor((aph->indexlx0)/largeGrid)*largeGrid);
            aph->indexly0 = static_cast<long>(floor((aph->indexly0)/largeGrid)*largeGrid);
            aph->indexcx0 = static_cast<long>(floor((aph->indexcx0)/coarseGrid)*coarseGrid);
            aph->indexcy0 = static_cast<long>(floor((aph->indexcy0)/coarseGrid)*coarseGrid);
            aph->indexmx0 = static_cast<long>(floor((aph->indexmx0)/mediumGrid)*mediumGrid);
            aph->indexmy0 = static_cast<long>(floor((aph->indexmy0)/mediumGrid)*mediumGrid);
            aph->indexfx0 = static_cast<long>(floor((aph->indexfx0)/fineGrid)*fineGrid);
            aph->indexfy0 = static_cast<long>(floor((aph->indexfy0)/fineGrid)*fineGrid);
            aph->indexlx1 = static_cast<long>(floor((aph->indexlx1)/largeGrid)*largeGrid);
            aph->indexly1 = static_cast<long>(floor((aph->indexly1)/largeGrid)*largeGrid);
            aph->indexcx1 = static_cast<long>(floor((aph->indexcx1)/coarseGrid)*coarseGrid);
            aph->indexcy1 = static_cast<long>(floor((aph->indexcy1)/coarseGrid)*coarseGrid);
            aph->indexmx1 = static_cast<long>(floor((aph->indexmx1)/mediumGrid)*mediumGrid);
            aph->indexmy1 = static_cast<long>(floor((aph->indexmy1)/mediumGrid)*mediumGrid);
            aph->indexfx1 = static_cast<long>(floor((aph->indexfx1)/fineGrid)*fineGrid);
            aph->indexfy1 = static_cast<long>(floor((aph->indexfy1)/fineGrid)*fineGrid);

            (angle->x) += (interpolate_bilinear_float_wrap(screen.turbulenceLargeX + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                         aph->dlx, aph->indexly0, aph->indexly1, aph->dly)*largeScale +
                         interpolate_bilinear_float_wrap(screen.turbulenceCoarseX + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                         aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy)*coarseScale +
                         interpolate_bilinear_float_wrap(screen.turbulenceMediumX + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                         aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)*mediumScale)*scaleOuter;

            (angle->y) += (interpolate_bilinear_float_wrap(screen.turbulenceLargeY + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                         aph->dlx, aph->indexly0, aph->indexly1, aph->dly)*largeScale +
                         interpolate_bilinear_float_wrap(screen.turbulenceCoarseY + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                         aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy)*coarseScale +
                         interpolate_bilinear_float_wrap(screen.turbulenceMediumY + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                         aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)*mediumScale)*scaleOuter;

            angle->z = smallAnglePupilNormalize(angle->x, angle->y);


        }

    } else {
        if (fabs(aph->time - aph->prtime) > screentol) {
            double randomi = random[aph->thread].uniform();
            double randomj = random[aph->thread].uniform();
            for (long i = 0; i < SCREEN_SIZE; i++) {
                for (long j = 0; j < SCREEN_SIZE; j++) {

                    find_linear_wrap(aph->xpos - aph->xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.large_sizeperpixel, SCREEN_SIZE, &(aph->indexlx0), &(aph->indexlx1), &(aph->dlx));
                    find_linear_wrap(aph->xpos - aph->xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.coarse_sizeperpixel, SCREEN_SIZE, &(aph->indexcx0), &(aph->indexcx1), &(aph->dcx));
                    find_linear_wrap(aph->xpos - aph->xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.medium_sizeperpixel, SCREEN_SIZE, &(aph->indexmx0), &(aph->indexmx1), &(aph->dmx));
                    find_linear_wrap(aph->xpos - aph->xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.fine_sizeperpixel,  SCREEN_SIZE, &(aph->indexfx0), &(aph->indexfx1), &(aph->dfx));
                    find_linear_wrap(aph->ypos - aph->yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.large_sizeperpixel, SCREEN_SIZE, &(aph->indexly0), &(aph->indexly1), &(aph->dly));
                    find_linear_wrap(aph->ypos - aph->yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.coarse_sizeperpixel, SCREEN_SIZE, &(aph->indexcy0), &(aph->indexcy1), &(aph->dcy));
                    find_linear_wrap(aph->ypos - aph->yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.medium_sizeperpixel, SCREEN_SIZE, &(aph->indexmy0), &(aph->indexmy1), &(aph->dmy));
                    find_linear_wrap(aph->ypos - aph->yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.fine_sizeperpixel,  SCREEN_SIZE,  &(aph->indexfy0),  &(aph->indexfy1), &(aph->dfy));

                    if (atmdebug == 0) {
                        if (mode == 2) {

                            *(screen.phasescreen + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phaseLarge + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                                 aph->dlx, aph->indexly0, aph->indexly1, aph->dly)+
                                 interpolate_bilinear_float_wrap(screen.phaseCoarse + aph->lindex,SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                                 aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy)+
                                 interpolate_bilinear_float_wrap(screen.phaseMedium + aph->lindex,SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                                 aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)+
                                 interpolate_bilinear_float_wrap(screen.phaseFine + aph->lindex  ,SCREEN_SIZE, aph->indexfx0, aph->indexfx1,
                                                                 aph->dfx, aph->indexfy0, aph->indexfy1, aph->dfy))*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        } else {

                            *(screen.phasescreen + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phaseMediumH + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                                 aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)+
                                 interpolate_bilinear_float_wrap(screen.phaseFineH + aph->lindex, SCREEN_SIZE, aph->indexfx0, aph->indexfx1,
                                                                 aph->dfx,  aph->indexfy0, aph->indexfy1, aph->dfy))*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        }

                    } else {

                        aph->indexlx0 = static_cast<long>(floor((aph->indexlx0)/largeGrid)*largeGrid);
                        aph->indexly0 = static_cast<long>(floor((aph->indexly0)/largeGrid)*largeGrid);
                        aph->indexcx0 = static_cast<long>(floor((aph->indexcx0)/coarseGrid)*coarseGrid);
                        aph->indexcy0 = static_cast<long>(floor((aph->indexcy0)/coarseGrid)*coarseGrid);
                        aph->indexmx0 = static_cast<long>(floor((aph->indexmx0)/mediumGrid)*mediumGrid);
                        aph->indexmy0 = static_cast<long>(floor((aph->indexmy0)/mediumGrid)*mediumGrid);
                        aph->indexfx0 = static_cast<long>(floor((aph->indexfx0)/fineGrid)*fineGrid);
                        aph->indexfy0 = static_cast<long>(floor((aph->indexfy0)/fineGrid)*fineGrid);
                        aph->indexlx1 = static_cast<long>(floor((aph->indexlx1)/largeGrid)*largeGrid);
                        aph->indexly1 = static_cast<long>(floor((aph->indexly1)/largeGrid)*largeGrid);
                        aph->indexcx1 = static_cast<long>(floor((aph->indexcx1)/coarseGrid)*coarseGrid);
                        aph->indexcy1 = static_cast<long>(floor((aph->indexcy1)/coarseGrid)*coarseGrid);
                        aph->indexmx1 = static_cast<long>(floor((aph->indexmx1)/mediumGrid)*mediumGrid);
                        aph->indexmy1 = static_cast<long>(floor((aph->indexmy1)/mediumGrid)*mediumGrid);
                        aph->indexfx1 = static_cast<long>(floor((aph->indexfx1)/fineGrid)*fineGrid);
                        aph->indexfy1 = static_cast<long>(floor((aph->indexfy1)/fineGrid)*fineGrid);

                        if (mode == 2) {

                            *(screen.phasescreen + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phaseLarge + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                                 aph->dlx, aph->indexly0, aph->indexly1, aph->dly)*largeScale +
                                 interpolate_bilinear_float_wrap(screen.phaseCoarse + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                                 aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy)*coarseScale +
                                 interpolate_bilinear_float_wrap(screen.phaseMedium + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                                 aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)*mediumScale +
                                 interpolate_bilinear_float_wrap(screen.phaseFine + aph->lindex, SCREEN_SIZE, aph->indexfx0, aph->indexfx1,
                                                                 aph->dfx, aph->indexfy0, aph->indexfy1, aph->dfy)*fineScale)*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        } else {
                            *(screen.phasescreen + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phaseMediumH + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                                 aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)*mediumScale +
                                 interpolate_bilinear_float_wrap(screen.phaseFineH + aph->lindex, SCREEN_SIZE, aph->indexfx0, aph->indexfx1,
                                                                 aph->dfx, aph->indexfy0, aph->indexfy1, aph->dfy)*fineScale)*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        }


                    }
                }
            }
        }

    }

}


double Image::cloudOpacity (long layer, Photon *aph) {

    double transmission, transmissionC, airmassl, densityScale5;

    transmission = 1.0;
    transmissionC = 1.0;
    if (aph->dvr == 0.0) {
        airmassl = 1.0;
    } else {
        airmassl = RADIUS_EARTH*(sqrt(pow(1 + height[layer]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1) -
                                 sqrt(pow(1 + height[layer + 1]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1))*
            sin(aph->dvr)/(height[layer] - height[layer + 1]);
    }
    densityScale5 = 1.0 + aerosolgradient*(aph->xpos*cos(aerosolAngle) + aph->ypos*sin(aerosolAngle))/(1000000.0);
    transmission *= exp(-((*(air.tau[5*(layer + 1) + 4] + aph->oindex)*densityScale5)*airmassl/air.airmassLayer[layer + 1]));

    if (cloudmean[layer] != 0 || cloudvary[layer] != 0) {
        transmissionC = pow(10.0, -0.4*(cloudmean[layer] +
                                       cloudvary[layer]*((double)(*(screen.cloud[layer] + aph->indexcx0*SCREEN_SIZE + aph->indexcy0)))));
        if (transmissionC > 1.0) transmissionC = 1.0;
    }
    transmission *= transmissionC;
    return(transmission);
}

double Image::cloudOpacityMoon (long layer, Photon *aph) {

    double transmission, transmissionC;

    transmission = 1.0;
    transmissionC = 1.0;
    transmission *= exp(-((*(air.tauMoon[5*(layer + 1) + 4] + aph->oindex))));

    if (cloudmean[layer] != 0 || cloudvary[layer] != 0) {
        transmissionC = pow(10.0, -0.4*(cloudmean[layer] +
                                       cloudvary[layer]*((double)(*(screen.cloud[layer] + aph->indexcx0*SCREEN_SIZE + aph->indexcy0)))));
        if (transmissionC > 1.0) transmissionC = 1.0;
    }
    transmission *= transmissionC;
    return(transmission);
}


double Image::atmosphereOpacity (Vector angle, long layer, Photon *aph) {

    double dvx, dvy, airmassl;
    double rindex;
    double densityScale, densityScale2, densityScale3, densityScale4, densityScale5;

    if (layer == -1) {

        dvx = -angle.x + zenith*sin(aph->shiftedAngle);
        dvy = -angle.y + zenith*cos(aph->shiftedAngle);
        aph->dvr = sqrt(dvx*dvx + dvy*dvy);
        if (aph->dvr == 0.0) {
            airmassl = 1.0;
        } else {
            airmassl = RADIUS_EARTH*(sqrt(pow(1 + 70.0/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1) -
                                     sqrt(pow(1 + height[layer + 1]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1))*
                sin(aph->dvr)/(70.0 - height[layer + 1]);
        }

        aph->oindex = find_linear(air.tauWavelength, 180001, aph->wavelength, &rindex);
        return(exp(-(*(air.tau[0] + aph->oindex))*airmassl/air.airmassLayer[layer + 1]));


    } else {

        if (aph->dvr == 0.0) {
            airmassl = 1.0;
        } else {
            airmassl = RADIUS_EARTH*(sqrt(pow(1 + height[layer]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1) -
                                     sqrt(pow(1 + height[layer + 1]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1))*
                sin(aph->dvr)/(height[layer] - height[layer + 1]);
        }

        densityScale = 1.0 + raygradient*(aph->xpos*cos(rayAngle) + aph->ypos*sin(rayAngle))/(1000000.0);
        densityScale2 = 1.0 + o3gradient*(aph->xpos*cos(o3Angle) + aph->ypos*sin(o3Angle))/(1000000.0);
        densityScale3 = 1.0 + o2gradient*(aph->xpos*cos(o2Angle) + aph->ypos*sin(o2Angle))/(1000000.0);
        densityScale4 = 1.0 + h2ogradient*(aph->xpos*cos(h2oAngle) + aph->ypos*sin(h2oAngle))/(1000000.0);
        densityScale5 = 1.0 + aerosolgradient*(aph->xpos*cos(aerosolAngle) + aph->ypos*sin(aerosolAngle))/(1000000.0);

        return(exp(-(*(air.tau[5*(layer + 1) + 0] + aph->oindex)*densityScale +
                     *(air.tau[5*(layer + 1) + 1] + aph->oindex)*densityScale2 +
                     *(air.tau[5*(layer + 1) + 2] + aph->oindex)*densityScale3 +
                     *(air.tau[5*(layer + 1) + 3] + aph->oindex)*densityScale4)*airmassl/air.airmassLayer[layer + 1]));

    }

}

double Image::atmosphereOpacityMoon (Vector angle, long layer, Photon *aph) {

    double rindex;

    if (layer == -1) {

        aph->oindex = find_linear(air.tauWavelength, 180001, aph->wavelength, &rindex);
        return(exp(-(*(air.tau[0] + aph->oindex))));

    } else {

        return(exp(-(*(air.tauMoon[5*(layer + 1) + 0] + aph->oindex) +
                     *(air.tauMoon[5*(layer + 1) + 1] + aph->oindex) +
                     *(air.tauMoon[5*(layer + 1) + 2] + aph->oindex) +
                     *(air.tauMoon[5*(layer + 1) + 3] + aph->oindex))));

    }

}


// transform to optic frame
void Image::transform( Vector *position, Vector *angle, long surfaceIndex, int focusFlag, Photon *aph) {

    double vprime1, vprime2, vprime3;

    // defocus
    position->z  =  position->z - perturbation.defocus[surfaceIndex] - surface.height[surfaceIndex];
    // decenter
    position->x = position->x - perturbation.decenterX[surfaceIndex] - surface.centerx[surfaceIndex];
    position->y = position->y - perturbation.decenterY[surfaceIndex] - surface.centery[surfaceIndex];

    if (focusFlag == 1) {
        int step;
        double shiftedValue;
        step = floor(aph->time/exptime*focusSteps);
        shiftedValue = step - (focusSteps - 1)/2.0;
        position->z += shiftedValue*focusStepZ;
        if (step == 0) shiftedValue -= 1.0;
        position->x += shiftedValue*focusStepX*cos(aph->shiftedAngle);
        position->y += shiftedValue*focusStepX*sin(aph->shiftedAngle);
    }

    // angles
    vprime1 = (*(perturbation.rotationmatrix + 9*surfaceIndex+0*3 + 0))*(angle->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 1))*(angle->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 2))*(angle->z);
    vprime2 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 0))*(angle->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 1))*(angle->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 2))*(angle->z);
    vprime3 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 0))*(angle->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 1))*(angle->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 2))*(angle->z);

    angle->x = vprime1;
    angle->y = vprime2;
    angle->z = vprime3;

    // position
    vprime1 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 0))*(position->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 1))*(position->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 2))*(position->z);
    vprime2 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 0))*(position->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 1))*(position->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 2))*(position->z);
    vprime3 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 0))*(position->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 1))*(position->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 2))*(position->z);

    position->x = vprime1;
    position->y = vprime2;
    position->z = vprime3;

    position->z = position->z + surface.height[surfaceIndex];
    position->x = position->x + surface.centerx[surfaceIndex];
    position->y = position->y + surface.centery[surfaceIndex];

}

// transform back to lab frame
void Image::transformInverse( Vector *position, Vector *angle, long surfaceIndex) {

    double vprime1, vprime2, vprime3;


    position->z = position->z - surface.height[surfaceIndex];
    position->x = position->x - surface.centerx[surfaceIndex];
    position->y = position->y - surface.centery[surfaceIndex];

    // angles
    vprime1 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 0))*(angle->x)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 1))*(angle->y)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 2))*(angle->z);
    vprime2 = (*(perturbation.inverserotationmatrix+9*surfaceIndex + 1*3 + 0))*(angle->x)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 1))*(angle->y)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 2))*(angle->z);
    vprime3 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 0))*(angle->x)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 1))*(angle->y)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 2))*(angle->z);

    angle->x = vprime1;
    angle->y = vprime2;
    angle->z = vprime3;

    // position
    vprime1 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 0))*(position->x) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 1))*(position->y) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 2))*(position->z);
    vprime2 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 0))*(position->x) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 1))*(position->y) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 2))*(position->z);
    vprime3 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 0))*(position->x) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 1))*(position->y) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 2))*(position->z);

    position->x = vprime1;
    position->y = vprime2;
    position->z = vprime3;


    // defocus
    position->z = position->z + perturbation.defocus[surfaceIndex] + surface.height[surfaceIndex];

    // decenter
    position->x = position->x + perturbation.decenterX[surfaceIndex] + surface.centerx[surfaceIndex];
    position->y = position->y + perturbation.decenterY[surfaceIndex] + surface.centery[surfaceIndex];

}

void Image::interceptDerivatives(Vector *normal, Vector position, long surfaceIndex) {

    double rxd, udmin = 0.0, wdmin = 0.0;
    long umin = 0, rx, wmin = 0;
    double normal3, normal2;
    double phi, r;
    double dx, dy;
    double normalpr;

    dx = position.x - surface.centerx[surfaceIndex];
    dy = position.y - surface.centery[surfaceIndex];

    r = sqrt(dx*dx + dy*dy);
    rx = find_linear(&surface.radius[SURFACE_POINTS*surfaceIndex], SURFACE_POINTS, r, &rxd);
    if (perturbation.zernikeflag[surfaceIndex] == 1) {
        if (r < RAYTRACE_TOLERANCE) {
            phi = atan2(dy, dx);
            if (phi < 0) phi += 2*PI;
            dx += RAYTRACE_TOLERANCE*cos(phi);
            dy += RAYTRACE_TOLERANCE*sin(phi);
            r = sqrt(dx*dx + dy*dy);
            rx = find_linear(&surface.radius[SURFACE_POINTS*surfaceIndex], SURFACE_POINTS, r, &rxd);
        }
        wmin = find_linear(perturbation.zernike_r_grid, PERTURBATION_POINTS, r/perturbation.rmax[surfaceIndex], &wdmin);
        phi = atan2(dy, dx);
        if (phi < 0) phi += 2*PI;
        umin = find_linear(perturbation.zernike_phi_grid, PERTURBATION_POINTS, phi, &udmin);
    }

    normal3 = interpolate_linear(&surface.normal[SURFACE_POINTS*surfaceIndex], rx, rxd);
    if (perturbation.zernikeflag[surfaceIndex] == 0) {
        normal->x = -normal3*dx/r;
        normal->y = -normal3*dy/r;
        normal->z = 1.0;
    } else {
        normal2 = 0;
        normalpr = normal3;
        normal3 += (interpolate_bilinear(perturbation.zernike_summed_nr_p +
                                         surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS, PERTURBATION_POINTS, umin, udmin, wmin, wdmin))
            /perturbation.rmax[surfaceIndex];
        normal2 += interpolate_bilinear(perturbation.zernike_summed_np_r +
                                        surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS, PERTURBATION_POINTS, umin, udmin, wmin, wdmin);
        normal->x = -normal3*dx/r + dy*normal2/(r*r);
        normal->y = -normal3*dy/r - dx*normal2/(r*r);
        normal->z = 1.0;
        // printf("%e %e %e %e %e %e %ld %ld %e %e %e %e %ld\n",normal2,normal3,dx,
               // dy,r,perturbation.rmax[surfaceIndex],umin,wmin,
// *(perturbation.zernike_summed_np_r + surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS + umin*PERTURBATION_POINTS + wmin),
//                *(perturbation.zernike_summed_np_r + surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS + (umin+1)*PERTURBATION_POINTS + wmin),
// *(perturbation.zernike_summed_np_r + surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS + umin*PERTURBATION_POINTS + wmin+1),
//                *(perturbation.zernike_summed_np_r + surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS + (umin+1)*PERTURBATION_POINTS + wmin +1),
// surfaceIndex);
    }
    normalize(normal);

}

int Image::chooseSurface (long *newSurf, long *oldSurf, Photon *aph) {

    if (aph->direction == 1) {
    trynextsurface:;
        if ((*newSurf) <= (*oldSurf) + 1) {
            (*newSurf)--;
            if ((*newSurf) == (*oldSurf)) (*newSurf) = (*oldSurf) - 1;
            if ((*newSurf) <= -1) (*newSurf) = (*oldSurf) + 2;
        } else {
            (*newSurf)++;
        }
        if ((*newSurf) > nsurf - 1) {
            return(1);
        } else {
            if (ghost[(*newSurf)] == 1) goto trynextsurface;
            return(0);
        }
    } else {
    trynextsurfaceback:;
        if ((*newSurf) >= (*oldSurf) - 1) {
            (*newSurf)++;
            if ((*newSurf) == (*oldSurf)) (*newSurf) = (*oldSurf) + 1;
            if ((*newSurf) > nsurf - 1) (*newSurf) = (*oldSurf) - 2;
        } else {
            (*newSurf)--;
        }
        if ((*newSurf) <= -1) {
            return(1);
        } else {
            if (ghost[(*newSurf)] == 1) goto trynextsurfaceback;
            return(0);
        }
    }

}


void Image::atmosphereDiffraction (Vector *angle, Photon *aph) {

    double radius, cr, cc, tf;
    fftw_plan pb;
    long i, j = 0, ix, jx;
    double dx, dy;
    double norm;
    double seeing;
    double rinner, router;
    double datamax;
    char* tempstring;

    // seeing = (totalseeing + 1e-6)*pow(1/cos(zenith), 0.6)*aph->wavelengthFactor;
    // norm = sqrt(0.0229*pow(0.98*(1e-4*aph->wavelength)/(seeing*ARCSEC),-5./3.))*3.75e8;

    // printf("a %e %e %e\n",seeing,norm,aph->wavelength);

    // seeing = (totalseeing + 1e-6)*pow(1/cos(zenith), 0.6)*aph->wavelengthFactor;
    // norm = sqrt(0.0229*pow(0.98*(1e-4*0.5)/(seeing*ARCSEC),-5./3.))*3.14e8;

    // printf("b %e %e %e\n",seeing,norm,aph->wavelength);

    if (atmospheremode < 2) { //turbulence off or atmosphere completely off
        seeing = 0.0;
        norm = 0.0;
        if (pupilscreenMode == 1) {
            tempstring = (char*)((instrdir + "/pupilscreen.fits").c_str());
            datamax = 255.0;
            datamax = screen.readScreen(9, screen.pupil_values, tempstring);
            tempstring = (char*)((instrdir + "/pupilscreen.fits").c_str());
            screen.paddingfactor = 8;
            screen.paddingfactor = screen.readScreen(10, screen.pupil_values, tempstring);
            tempstring = (char*)((instrdir + "/pupilscreen.fits").c_str());
            screen.pupilscreenscale = 1.0;
            screen.pupilscreenscale = screen.readScreen(11, screen.pupil_values, tempstring);
        } else {
            screen.paddingfactor = 32;
            rinner = screen.fine_sizeperpixel*SCREEN_SIZE*surface.innerRadius[0]/(2*surface.outerRadius[0]*screen.paddingfactor);
            router = screen.fine_sizeperpixel*SCREEN_SIZE/(2*screen.paddingfactor);
        }
    } else {
        seeing = (totalseeing + 1e-6)*pow(1/cos(zenith), 0.6)*aph->wavelengthFactor;
        norm = sqrt(0.0229*pow(0.98*(1e-4*0.5)/(seeing*ARCSEC), -5.0/3.0))*3.25e8;
        screen.paddingfactor = 1;
        rinner = surface.innerRadius[0]/screen.paddingfactor;
        router = surface.outerRadius[0]/screen.paddingfactor;
    }

    for (i = 0; i < SCREEN_SIZE; i++) {
        for (j = 0; j < SCREEN_SIZE; j++) {
            radius = sqrt((i - SCREEN_SIZE/2 + 0.5)*(i - SCREEN_SIZE/2 + 0.5)+
                          (j - SCREEN_SIZE/2 + 0.5)*(j - SCREEN_SIZE/2 + 0.5))*screen.fine_sizeperpixel;
            dx = (i - SCREEN_SIZE/2 + 0.5)*screen.fine_sizeperpixel*cos(aph->shiftedAngle) +
                (j - SCREEN_SIZE/2 + 0.5)*screen.fine_sizeperpixel*sin(aph->shiftedAngle);
            dy = -(i - SCREEN_SIZE/2 + 0.5)*screen.fine_sizeperpixel*sin(aph->shiftedAngle) +
                (j - SCREEN_SIZE/2 + 0.5)*screen.fine_sizeperpixel*cos(aph->shiftedAngle);
            if (atmospheremode < 2 && pupilscreenMode == 1) {
                screen.inscreen[SCREEN_SIZE*i + j][0] = screen.pupil_values[j*SCREEN_SIZE + i]/datamax;
                screen.inscreen[SCREEN_SIZE*i + j][1] = 0.0;
            } else {
                if (radius > rinner && radius < router) {
                    screen.inscreen[SCREEN_SIZE*i + j][0] = cos(screen.phasescreen[i*SCREEN_SIZE + j]*norm);
                    screen.inscreen[SCREEN_SIZE*i + j][1] = sin(screen.phasescreen[i*SCREEN_SIZE + j]*norm);
                    // for (k = 0; k<nspid; k++) {
                    //     if (obstruction.type[k] = =1) {
                    //         if (fabs(dx-obstruction.center[k]) < obstruction.width[k]) {
                    //             screen.inscreen[SCREEN_SIZE*i+j][0] = 0.0;
                    //             screen.inscreen[SCREEN_SIZE*i+j][1] = 0.0;
                    //         }
                    //     }
                    //     if (obstruction.type[k]==2) {
                    //         if (fabs(dy-obstruction.center[k]) < obstruction.width[k]) {
                    //             screen.inscreen[SCREEN_SIZE*i+j][0] = 0.0;
                    //             screen.inscreen[SCREEN_SIZE*i+j][1] = 0.0;
                    //         }
                    //     }
                    // }
                } else {
                    screen.inscreen[SCREEN_SIZE*i + j][0] = 0.0;
                    screen.inscreen[SCREEN_SIZE*i + j][1] = 0.0;
                }
            }
        }
    }


    pb = fftw_plan_dft_2d(SCREEN_SIZE, SCREEN_SIZE, screen.inscreen, screen.outscreen, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(pb);
    fftw_destroy_plan(pb);

    double total = 0.0;
    double value;
    for (i = 0; i<SCREEN_SIZE; i++) {
        for (j = 0; j<SCREEN_SIZE; j++) {
            ix = i + SCREEN_SIZE/2;
            jx = j + SCREEN_SIZE/2;
            ix = ix % (SCREEN_SIZE);
            jx = jx % (SCREEN_SIZE);
            if (atmospheremode < 2 && pupilscreenMode == 1) {
                value = pow(screen.outscreen[i*SCREEN_SIZE + j][0], 2) + pow(screen.outscreen[i*SCREEN_SIZE + j][1], 2);
                total += value;
                screen.focalscreen[ix*SCREEN_SIZE + jx] = value;
            } else {
                screen.focalscreen[ix*SCREEN_SIZE + jx] = pow(screen.outscreen[i*SCREEN_SIZE + j][0], 2)+
                    pow(screen.outscreen[i*SCREEN_SIZE + j][1], 2);
            }
        }
    }

    double cump = 0.0;
    tf = 0.0;
    for (i = 0; i < SCREEN_SIZE; i++) {
        for (j = 0; j < SCREEN_SIZE; j++) {
            tf += screen.focalscreen[i*SCREEN_SIZE + j];
            if (atmospheremode < 2 && pupilscreenMode == 1) { // normalize
                screen.focalscreencum[i*SCREEN_SIZE + j] = cump;
                cump += screen.focalscreen[i*SCREEN_SIZE + j]/total;
            }
        }
    }
    cr = random[aph->thread].uniform()*tf;
    cc = 0.0;
    for (i = 0; i<SCREEN_SIZE; i++) {
        for (j = 0; j<SCREEN_SIZE; j++) {
            cc += screen.focalscreen[i*SCREEN_SIZE + j];
            if (cr < cc) goto breakphoton;
        }
    }

breakphoton:;

    if (atmospheremode > 1) { // else, do properly in secondKick instead
        angle->x = angle->x + (i-SCREEN_SIZE/2 + random[aph->thread].uniform()-0.5)*aph->wavelength*1e-3/(SCREEN_SIZE*screen.fine_sizeperpixel);
        angle->y = angle->y + (j-SCREEN_SIZE/2 + random[aph->thread].uniform()-0.5)*aph->wavelength*1e-3/(SCREEN_SIZE*screen.fine_sizeperpixel);
        angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    }

    aph->prtime = aph->time;

}


int Image::bloom (int saturatedFlag, Photon *aph) {

    long newxpos, oldxpos, stepadd, iii;
    long startstep = 0;
    long location;

    location = chip.nampx*(aph->yPos - miny) + (aph->xPos - minx);

    if (aph->xPos >= chip.midpoint) {
        stepadd = 1;
    } else {
        stepadd = 0;
    }

    if (*(state.satupmap + location) >= 0 || *(state.satdownmap + location) >= 0) {

        oldxpos = aph->xPos;
        if (*(state.satupmap + location) < 0) startstep = aph->xPos - *(state.satdownmap + location);
        if (*(state.satdownmap + location) < 0) startstep = *(state.satupmap + location) - aph->xPos;
        if (*(state.satupmap + location) >= 0 && *(state.satdownmap + location) >= 0) {
            if (aph->xPos - *(state.satdownmap + location) < *(state.satupmap + location) - aph->xPos)
                startstep = aph->xPos - *(state.satdownmap + location); else
                startstep = *(state.satupmap + location) - aph->xPos;
        }

        for (iii = startstep; iii < chip.midpoint; iii++) {

            if (*(state.satupmap + location) >= 0) {
                newxpos = oldxpos + iii;
                if (newxpos < (chip.midpoint + stepadd*chip.midpoint)) {
                    if (*(state.focal_plane + chip.nampx*(aph->yPos - miny) + (newxpos - minx)) < well_depth) {
                        aph->xPos = newxpos;
                        if (saturatedFlag == 1) *(state.satupmap + location) = newxpos;
                        return(0);
                    } else {
                        if (newxpos > *(state.satupmap + location) && saturatedFlag == 1) *(state.satupmap + location) = newxpos;
                    }
                } else {
                    if (saturatedFlag == 1) *(state.satupmap + location) = -1;
                }
            }

            if (*(state.satdownmap + location) >= 0) {
                newxpos = oldxpos - iii;
                if (newxpos >= stepadd*chip.midpoint) {
                    if (*(state.focal_plane + chip.nampx*(aph->yPos-miny) + (newxpos - minx)) < well_depth) {
                        aph->xPos = newxpos;
                        if (saturatedFlag == 1) *(state.satdownmap + location) = newxpos;
                        return(0);
                    } else {
                        if (newxpos < *(state.satdownmap + location) && saturatedFlag == 1) *(state.satdownmap + location) = newxpos;
                    }
                } else {
                    if (saturatedFlag == 1) *(state.satdownmap + location) = -1;
                }
            }

        }
    }


    return(1);

}



void Image::saturate (long source, Vector* largeAngle, Photon *aph, double shiftX, double shiftY)

{
    long location, origlocation;
    std::atomic<unsigned long> leftover;
    std::atomic<unsigned long> wellDepth;
    long minrad;

    leftover = aph->sourceOver_m;
    origlocation = chip.nampx*(aph->yPos - miny) + (aph->xPos - minx);
    location = origlocation;

    if (aph->xPos >= minx && aph->xPos <= maxx && aph->yPos >= miny && aph->yPos <= maxy) {
        wellDepth = static_cast<unsigned long>(well_depth);

    rebloom:;
        // pthread_mutex_lock(&lock.lock1);
        *(state.focal_plane + location) += leftover;
        // pthread_mutex_unlock(&lock.lock1);
        if (*(state.focal_plane + location) > well_depth) {
            leftover = *(state.focal_plane + location) - wellDepth;
            *(state.focal_plane + location) = well_depth;
            if (blooming == 1) {
                if (bloom(1, aph)) goto fullysat;
                location = chip.nampx*(aph->yPos - miny) + (aph->xPos - minx);
                goto rebloom;
            }
        }

    fullysat:;


        if (*(state.focal_plane + origlocation) >= well_depth) {
            aph->saturationFlag = 1;
            if (aph->ghostFlag == 0 && sources.spatialtype[source] != 4 && sources.spatialtype[source] != 1  && sources.type[source] >= 6) {
                minrad = static_cast<long>(fabs((largeAngle->y + shiftY)/DEGREE*platescale/pixsize)) - satbuffer;
                if (minrad < 0) minrad = 0;
                if (minrad == (aph->sourceSaturationRadius + 1)) {
                    // if (sources.id[source]=="3333362765") printf("A %ld %ld %lf %lf %lf %lf\n",minrad,satbuffer,((largeAngle->y)/DEGREE*platescale/pixsize),((largeAngle->x)/DEGREE*platescale/pixsize),((shiftY)/DEGREE*platescale/pixsize),((shiftX)/DEGREE*platescale/pixsize));
                    aph->sourceSaturationRadius = minrad;
               }
            }
        }


    } else {

        if (aph->ghostFlag == 0 && sources.spatialtype[source] != 4 && sources.spatialtype[source] != 1
            && sources.type[source] >= 6 && aph->saturationFlag == 0) {
            if (aph->xPos > 1e6 || aph->xPos < -1e6 || aph->yPos > 1e6 || aph->yPos < -1e6) {
                // printf("Error %ld %ld\n",aph->xPos,aph->yPos);
            } else {
            minrad = 0;
            double deltaX, deltaY;
            deltaX = (largeAngle->x + shiftX)/DEGREE*platescale/pixsize;
            deltaY = (largeAngle->y + shiftY)/DEGREE*platescale/pixsize;
            if (aph->yPos - deltaY > maxy && aph->yPos > maxy) {
                if (-deltaY - satbuffer > minrad) minrad = static_cast<long>(-deltaY - satbuffer);
            }
            if (aph->yPos - deltaY < miny && aph->yPos < miny) {
                if (deltaY - satbuffer > minrad) minrad = static_cast<long>(deltaY - satbuffer);
            }
            if (aph->xPos - deltaX > maxx && aph->xPos > maxx) {
                if (-deltaX - satbuffer > minrad) minrad = static_cast<long>(-deltaX - satbuffer);
            }
            if (aph->xPos - deltaX < minx && aph->xPos < minx) {
                if (deltaX - satbuffer > minrad) minrad = static_cast<long>(deltaX - satbuffer);
            }
            if (minrad < 0) minrad = 0;
            if (minrad == (aph->sourceSaturationRadius + 1)) {
                // if (sources.id[source]=="3333362765") printf("B %ld %ld\n",minrad,satbuffer);
                aph->sourceSaturationRadius = minrad;
            }
            }
        }


    }

}

int Image::findSurface (Vector angle, Vector position, double *distance, long surfaceIndex, Photon *aph) {


    double zmin;
    long miss;

    double a, b, c, l;
    double disc, l1, l2;
    double tol;
    double discrepancy;
    double zv;
    double radcurv;
    int iter = 0;

    if (surface.radiusCurvature[surfaceIndex] == 0.0) {

        l = (surface.height[surfaceIndex] - position.z)/angle.z;
        miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
        discrepancy = fabs(position.z + angle.z*l - zv);
        while (discrepancy > RAYTRACE_TOLERANCE) {
            l = l + (zv - (position.z + angle.z*l))/angle.z;
            miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
            discrepancy = fabs(position.z + angle.z*l - zv);
            iter++;
            if (iter > RAYTRACE_MAX_ITER) goto iterjump;
        }
        zmin = discrepancy;
        *distance = l;
    } else {
        radcurv = surface.radiusCurvature[surfaceIndex];
        a = angle.x*angle.x + angle.y*angle.y + (1 + surface.conic[surfaceIndex])/2.0*angle.z*angle.z;
        b = (1 + surface.conic[surfaceIndex])*angle.z*(position.z - surface.height[surfaceIndex]) +
            2.0*angle.x*position.x + 2.0*angle.y*position.y - 2.0*radcurv*angle.z;
        c = (1 + surface.conic[surfaceIndex])/2.0*(position.z*position.z - 2.0*position.z*surface.height[surfaceIndex] +
                                                   surface.height[surfaceIndex]*surface.height[surfaceIndex]) +
            position.x*position.x + position.y*position.y - 2.0*radcurv*(position.z - surface.height[surfaceIndex]);
        disc = b*b - 4*a*c;
        if (a == 0.0) {
            l = -c/b;
            if (l <= 0.0 || b == 0.0) {
                tol = RAYTRACE_MAX_LENGTH;
                miss = goldenBisectSurface(l - tol, l, l + tol, &zmin, angle, position, distance, surfaceIndex, aph);
            } else {
                miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
                discrepancy = fabs(position.z + angle.z*l - zv);
                while (discrepancy > RAYTRACE_TOLERANCE) {
                    l = l + (zv - (position.z + angle.z*l))/angle.z;
                    miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
                    discrepancy = fabs(position.z + angle.z*l - zv);
                    iter++;
                    if (iter > RAYTRACE_MAX_ITER) goto iterjump;
                }
                zmin = discrepancy;
                *distance = l;
            }
        } else {
            if (disc > 0) {
                l1 = ((-b + sqrt(disc))/2.0/a);
                l2 = ((-b - sqrt(disc))/2.0/a);
                l = 0.0;
                if ((l2 > 0) && (l1 > 0)) {
                    if (l2 < l1) {
                        l = l2;
                    } else {
                        l = l1;
                    }
                }
                if ((l2 > 0) && (l1 < 0)) l = l2;
                if ((l2 < 0) && (l1 > 0)) l = l1;
                if (l == 0.0) {
                    tol = RAYTRACE_MAX_LENGTH;
                    miss = goldenBisectSurface(l - tol, l, l + tol, &zmin, angle, position, distance, surfaceIndex, aph);
                } else {
                    miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
                    discrepancy = fabs(position.z + angle.z*l - zv);
                    while (discrepancy > RAYTRACE_TOLERANCE) {
                        l = l + (zv - (position.z + angle.z*l))/angle.z;
                        miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
                        discrepancy = fabs(position.z + angle.z*l - zv);
                        iter++;
                        if (iter > RAYTRACE_MAX_ITER || miss == 1) goto iterjump;
                    }
                    zmin = discrepancy;
                    *distance = l;
                }
            } else {
            iterjump:;
                l = 0.0;
                tol = RAYTRACE_MAX_LENGTH;
                miss = goldenBisectSurface(l - tol, l, l + tol, &zmin, angle, position, distance, surfaceIndex, aph);
            }
        }
    }

    if (zmin > RAYTRACE_ERROR || miss == 1) {
        return(1);
    } else {
        return(0);
    }

}


int Image::goldenBisectSurface(double a, double b, double c, double *z, Vector angle, Vector position,
                               double *distance, long surfaceIndex, Photon *aph) {

    double f0, f1, f2, f3;
    double x0, x1, x2, x3;
    double zv;
    long miss;

    x0 = a;
    x3 = c;
    if (fabs(c - b) > fabs(b - a)) {
        x1 = b;
        x2 = (b + HALFSQ5M1*(c - b));
    } else {
        x2 = b;
        x1 = (b - HALFSQ5M1*(b - a));
    }

    miss = getIntercept(position.x + angle.x*x1, position.y + angle.y*x1, &zv, surfaceIndex, aph);
    f1 = fabs(position.z + angle.z*x1 - zv);
    miss = getIntercept(position.x + angle.x*x2, position.y + angle.y*x2, &zv, surfaceIndex, aph);
    f2 = fabs(position.z + angle.z*x2 - zv);

    while (fabs(x2 - x1) > RAYTRACE_TOLERANCE) {

        if (f2 < f1) {
            x0 = x1;
            x1 = x2;
            x2 = HALF3MSQ5*x1 + HALFSQ5M1*x3;
            miss = getIntercept(position.x + angle.x*x2, position.y + angle.y*x2, &zv, surfaceIndex, aph);
            f0 = f1;
            f1 = f2;
            f2 = fabs(position.z + angle.z*x2 - zv);
        } else {
            x3 = x2;
            x2 = x1;
            x1 = HALF3MSQ5*x2 + HALFSQ5M1*x0;
            miss = getIntercept(position.x + angle.x*x1, position.y + angle.y*x1, &zv, surfaceIndex, aph);
            f3 = f2;
            f2 = f1;
            f1 = fabs(position.z + angle.z*x1 - zv);
        }
    }

    if (f1 < f2) {
        *z = f1;
        *distance = x1;
    } else {
        *z = f2;
        *distance = x2;
    }

    return(miss);

}

int Image::getIntercept(double x, double y, double *z, long surfaceIndex, Photon *aph) {

    double r, phi;
    double uu, ww, tt;
    double dx, dy;
    long ttint;
    long lSurfaceIndex;

    dx = x - surface.centerx[surfaceIndex];
    dy = y - surface.centery[surfaceIndex];
    r = sqrt(dx*dx + dy*dy);
    lSurfaceIndex = SURFACE_POINTS*surfaceIndex;
    if ((r > surface.radius[lSurfaceIndex + SURFACE_POINTS - 1]) ||
        (r < surface.radius[lSurfaceIndex + 0])) {
        return(1);
    }

    ttint = find_linear(&surface.radius[lSurfaceIndex], SURFACE_POINTS, r, &tt);
    *z = interpolate_linear(&surface.profile[lSurfaceIndex], ttint, tt);

    find(&surface.radiusArea[lSurfaceIndex], PERTURBATION_POINTS, r, &(aph->vvint));

    aph->wwint = find_linear(perturbation.zernike_r_grid, PERTURBATION_POINTS, r/perturbation.rmax[surfaceIndex], &ww);
    phi = atan2(dy, dx);
    if (phi < 0) phi += 2*PI;
    aph->uuint = find_linear(perturbation.zernike_phi_grid, PERTURBATION_POINTS, phi, &uu);
    if (perturbation.zernikeflag[surfaceIndex] == 1) {
        *z += interpolate_bilinear(perturbation.zernike_summed + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex,
                                   PERTURBATION_POINTS, aph->uuint, uu, aph->wwint, ww);
    }

    return(0);

}

int Image::getDeltaIntercept(double x, double y, double *z, long surfaceIndex, Photon *aph) {

    double r, phi;
    double uu, ww;
    long uuint, wwint;
    double dx, dy;
    int miss;

    miss = 0;

    dx = x - surface.centerx[surfaceIndex];
    dy = y - surface.centery[surfaceIndex];
    r = sqrt(dx*dx + dy*dy);


    find(&surface.radiusArea[SURFACE_POINTS*surfaceIndex], SURFACE_POINTS, r, &(aph->vvint));
    if ((r > surface.radius[SURFACE_POINTS*surfaceIndex + SURFACE_POINTS - 1]) ||
        (r < surface.radius[SURFACE_POINTS*surfaceIndex + 0])) miss = 1;
    *z = 0.0;

    if (perturbation.zernikeflag[surfaceIndex] == 1) {
        wwint = find_linear(perturbation.zernike_r_grid, PERTURBATION_POINTS, r/perturbation.rmax[surfaceIndex], &ww);
        phi = atan2(dy, dx);
        if (phi < 0) phi += 2*PI;
        uuint = find_linear(perturbation.zernike_phi_grid, PERTURBATION_POINTS, phi, &uu);
        *z += interpolate_bilinear(perturbation.zernike_summed + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex,
                                   PERTURBATION_POINTS, uuint, uu, wwint, ww);
    }

    return(miss);
}


int Image::getWavelengthTime (Photon *aph, long source) {

    double tempf1;
    long index;
    double dustvalue = 1.0;

    // select wavelength
    if (sources.spatialtype[source] != OPD) {
        tempf1 = random[aph->thread].uniform();
        find(sedC + sedPtr[sources.sedptr[source]], sedN[sources.sedptr[source]], tempf1, &index);
        if (sources.type[source] == 1 || (sources.type[source] == 0 && domewave != 0.0)) {
            aph->wavelength = *(sedW + sedPtr[sources.sedptr[source]] + index);
        } else {
            aph->wavelength = interpolate(sedW + sedPtr[sources.sedptr[source]], sedC + sedPtr[sources.sedptr[source]], tempf1, index);
        }
        aph->wavelength = aph->wavelength/1000.0;
    } else {
        aph->wavelength = sources.gamma1[source]/1000.0;
    }

    // dust at source
    if (sources.dusttypez[source] != 0) {
        if (sources.dusttypez[source] == 1) dustvalue = dust.ccm(aph->wavelength, maxwavelength, sources.dustparz[source][0], sources.dustparz[source][1]);
        if (sources.dusttypez[source] == 2) dustvalue = dust.calzetti(aph->wavelength, maxwavelength, sources.dustparz[source][0], sources.dustparz[source][1]);
        if (random[aph->thread].uniform() > dustvalue) { return(1);}
    }

    // redshift aph
    aph->wavelength = (aph->wavelength)*(1 + sources.redshift[source]);
    if (aph->wavelength < minwavelength/1000.0 || aph->wavelength >= maxwavelength/1000.0) return(1);

    // dust at z=0
    if (sources.dusttype[source] != 0) {
        if (sources.dusttype[source] == 1) dustvalue = dust.ccm(aph->wavelength, maxwavelength, sources.dustpar[source][0], sources.dustpar[source][1]);
        if (sources.dusttype[source] == 2) dustvalue = dust.calzetti(aph->wavelength, maxwavelength, sources.dustpar[source][0], sources.dustpar[source][1]);
        if (random[aph->thread].uniform() > dustvalue) {return(1);}
    }

    // choose time
    aph->time = static_cast<double>((random[aph->thread].uniform())*exptime);

    // choose polarization
    if (random[aph->thread].uniform() > 0.5) aph->polarization=1; else aph->polarization=-1;

    return(0);

}

void Image::getAngle (Vector *angle, double time, long source) {

    if (sources.spatialtype[source] == MOVINGPOINT) {
        angle->x = sources.vx[source] + (sources.spatialpar[source][0]*cos(rotatez)+
                                         sources.spatialpar[source][1]*(-sin(rotatez)))*
            (time - exptime/2.0 + timeoffset)*ARCSEC;
        angle->y = sources.vy[source] + (sources.spatialpar[source][0]*sin(rotatez)+
                                         sources.spatialpar[source][1]*cos(rotatez))*
            (time - exptime/2.0 + timeoffset)*ARCSEC;
    } else {
        angle->x = sources.vx[source];
        angle->y = sources.vy[source];
    }
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

}

void Image::getDeltaAngle(Vector *angle, Vector *position, long source, double *shiftX, double *shiftY, int thread, int *initGal, Photon *aph) {

    double dx = 0.0;
    double dy = 0.0;

    if (sources.spatialtype[source] != POINT && sources.spatialtype[source] != MOVINGPOINT
        && sources.spatialtype[source] != OPD) {
        if (sources.spatialtype[source] == SERSIC2D) {
            galaxy.sersic2d(sources.spatialpar[source][0], sources.spatialpar[source][1],
                            sources.spatialpar[source][2]*DEGREE, sources.spatialpar[source][3],
                            &dx, &dy, thread);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == GAUSSIAN) {
            dx = random[aph->thread].normal()*ARCSEC*sources.spatialpar[source][0];
            dy = random[aph->thread].normal()*ARCSEC*sources.spatialpar[source][0];
       } else if (sources.spatialtype[source] == PINHOLE) {
            double r = sqrt(random[aph->thread].uniform()*sources.spatialpar[source][3]*sources.spatialpar[source][3]);
            double phi = random[aph->thread].uniform()*2*PI;
            double xt = position->x - r*cos(phi) - sources.spatialpar[source][0];
            double yt = position->y - r*sin(phi) - sources.spatialpar[source][1];
            double zt = position->z - surface.height[0] - sources.spatialpar[source][2];
            double finiteDistanceR = sqrt(xt*xt + yt*yt + zt*zt);
            dx = -xt/finiteDistanceR;
            dy = yt/finiteDistanceR;
        } else if (sources.spatialtype[source] == SERSIC) {
            galaxy.sersic(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3]*DEGREE,
                          sources.spatialpar[source][4]*DEGREE, sources.spatialpar[source][5],
                          &dx, &dy, thread);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == SERSICDISK) {
            galaxy.sersicDisk(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3]*DEGREE,
                          sources.spatialpar[source][4]*DEGREE, sources.spatialpar[source][5],
                              &dx, &dy, thread);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        }else if (sources.spatialtype[source] == DISTORTEDSPHERE) {
            galaxy.distortedSphere(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3],
                          sources.spatialpar[source][4], sources.spatialpar[source][5],
                                   sources.spatialpar[source][6], sources.spatialpar[source][7],
                                   sources.spatialpar[source][8], sources.spatialpar[source][9],
                                   &dx, &dy, thread);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == SERSICCOMPLEX) {
            galaxy.sersicComplex(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3]*DEGREE,
                          sources.spatialpar[source][4]*DEGREE, sources.spatialpar[source][5],
                          sources.spatialpar[source][6], sources.spatialpar[source][7],
                          sources.spatialpar[source][8], sources.spatialpar[source][9],
                          sources.spatialpar[source][10]*DEGREE, sources.spatialpar[source][11],
                          sources.spatialpar[source][12], sources.spatialpar[source][13]*DEGREE,
                                 thread, initGal, &dx, &dy);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        }else if (sources.spatialtype[source] == SERSICDISKCOMPLEX) {
            galaxy.sersicComplex(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3]*DEGREE,
                          sources.spatialpar[source][4]*DEGREE, sources.spatialpar[source][5],
                          sources.spatialpar[source][6], sources.spatialpar[source][7],
                          sources.spatialpar[source][8], sources.spatialpar[source][9],
                          sources.spatialpar[source][10]*DEGREE, sources.spatialpar[source][11],
                          sources.spatialpar[source][12], sources.spatialpar[source][13]*DEGREE,
                                 thread, initGal, &dx, &dy);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == IMAGE) {
            float imagepixel, localcumulative;
            long ii, jj, selectedsky;
            selectedsky = static_cast<long>(sources.spatialpar[source][2]);
            imagepixel = (random[aph->thread].uniform()*cumulative[selectedsky]);
            localcumulative = 0;
            for (ii = 0; ii < naxesb[selectedsky][0]; ii++) {
                if (cumulativex[selectedsky][ii] > imagepixel) {
                    ii--;
                    goto jumpa;
                }
            }
            ii--;
        jumpa:;
            localcumulative = cumulativex[selectedsky][ii];
            for (jj = 0; jj < naxesb[selectedsky][1]; jj++) {
                localcumulative += (*(tempptr[selectedsky] + (ii)*naxesb[selectedsky][1] + jj));
                if (localcumulative > imagepixel) {
                    jj--;
                    goto jumpb;
                }
            }
            jj--;
        jumpb:;
            dx = ((jj - (naxesb[selectedsky][1]/2.0))*cos(sources.spatialpar[source][1]*DEGREE)+
                (-(ii - (naxesb[selectedsky][0]/2.0)))*sin(sources.spatialpar[source][1]*DEGREE))*
                sources.spatialpar[source][0]*ARCSEC;
            dy = -((jj - (naxesb[selectedsky][1]/2.0))*(-sin(sources.spatialpar[source][1]*DEGREE))+
                (-(ii - (naxesb[selectedsky][0]/2.0)))*cos(sources.spatialpar[source][1]*DEGREE))*
                sources.spatialpar[source][0]*ARCSEC;
        }
        if (sources.gamma1[source] != 0.0 || sources.gamma2[source] != 0.0) {
            double dxp = dx*(1 + sources.gamma1[source] - sources.kappa[source]) - dy*(sources.gamma2[source]);
            double dyp = dy*(1 - sources.gamma1[source] - sources.kappa[source]) - dx*(sources.gamma2[source]);
            dx = dxp;
            dy = dyp;
        }
        *shiftX = - dx*cos(rotatez) - dy*(sin(rotatez));
        *shiftY = - dx*(sin(rotatez)) + dy*cos(rotatez);
        // angle->x += *shiftX;
        // angle->y += *shiftY;
        // angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    }


}



int Image::photonSiliconPropagate(Vector *angle, Vector *position, double lambda, Vector normal, double dh, long waveSurfaceIndex, Photon *aph) {

    double travel, dtravel;
    long yindex;
    double ryindex;
    Vector origAngle;
    double dead;

    aph->z0 = position->z;
    if (detectorcollimate == 1) {
        angle->z = angle->z/fabs(angle->z);
        angle->x = 0.0;
        angle->y = 0.0;
    }
    vectorCopy(*angle, &origAngle);

    // silicon refraction
    yindex = find_linear(silicon.wavelengthGrid, silicon.numWavelength, aph->wavelength, &ryindex);
    double nSi = interpolate_linear(silicon.indexRefraction, yindex, ryindex);
    refract(angle, normal, 1, nSi);

    // photo-electron conversion
    aph->xindex = find_linear(silicon.temperatureGrid, silicon.numTemperature, sensorTempNominal + sensorTempDelta, &(aph->rxindex));
    double mfp = interpolate_bilinear(silicon.meanFreePath, silicon.numWavelength, aph->xindex, aph->rxindex, yindex, ryindex);
    if (std::isinf(mfp) || std::isnan(mfp)) return(1);
    double randNum;
    double conversion;
    if (photoelectric == 1) {
        conversion = 1.0 - exp(-2*sensorthickness/1e3/fabs(angle->z)/mfp);
    } else {
        conversion = 1.0;
        mfp = 0.0;
    }
    (aph->counter)++;
    if (conversion > state.dynamicTransmission[2*natmospherefile + 2*nsurf + 1 + waveSurfaceIndex]) {
        state.dynamicTransmission[2*natmospherefile + 2*nsurf + 1 + waveSurfaceIndex] = conversion;
    }
    if (aph->counter <= aph->maxcounter) {
        randNum = aph->saveRand[aph->counter];
    } else {
        randNum = random[aph->thread].uniform();
    }
    if (randNum > conversion) return(1);
    travel = mfp*(-log(1.0 - randNum));
    aph->location = chip.nampx*(aph->yPos - miny) + (aph->xPos - minx);
    if (aph->xPos < minx || aph->xPos > maxx || aph->yPos < miny || aph->yPos > maxy) {
        long xL, yL;
        xL = aph->xPos;
        yL = aph->yPos;
        if (aph->xPos < minx) xL = minx;
        if (aph->xPos > maxx) xL = maxx;
        if (aph->yPos < miny) yL = miny;
        if (aph->yPos > maxy) yL = maxy;
        aph->location = chip.nampx*(yL - miny) + (xL - minx);
    }
    if (deadlayer == 1) {
        dead = silicon.deadLayer[aph->location];
    } else {
        dead = 0.0;
    }
    if (fabs(travel*angle->z) < dead/1e6) return(1);
    aph->collect_z = aph->z0 + angle->z/fabs(angle->z)*sensorthickness/1e3;
    if (fabs(travel*angle->z) >= sensorthickness/1e3) {
        double conversion = fringing(origAngle, normal, aph->wavelength, nSi, static_cast<double>(sensorthickness) + dh*1000.0*fringeflag, mfp, aph->polarization);
        if (random[aph->thread].uniform() > conversion) return(1);
        aph->ghostFlag = 1;
        while (travel > 0) {
            if (fabs(travel*angle->z) >= sensorthickness/1e3) {
                dtravel = fabs(sensorthickness/1e3/angle->z);
                propagate(position, *angle, dtravel);
                reflect(angle, normal);
                travel -= dtravel;
            } else {
                propagate(position, *angle, travel);
                travel = 0;
            }
        }
    } else {
        propagate(position, *angle, travel);
    }
    return(0);
}

int Image::electronSiliconPropagate(Vector *angle, Vector *position, Photon *aph) {

    long windex, zindex, uindex;
    double rwindex, rzindex, ruindex;
    double dopant;

    // charge diffusion
    zindex = find_linear(silicon.thicknessGrid, silicon.numThickness, sensorthickness/1e4 - fabs((position->z - aph->z0)/10.0), &rzindex);
    aph->location = chip.nampx*(aph->yPos - miny) + (aph->xPos - minx);
    if (aph->xPos < minx || aph->xPos > maxx || aph->yPos < miny || aph->yPos > maxy) {
        long xL, yL;
        xL = aph->xPos;
        yL = aph->yPos;
        if (aph->xPos < minx) xL = minx;
        if (aph->xPos > maxx) xL = maxx;
        if (aph->yPos < miny) yL = miny;
        if (aph->yPos > maxy) yL = maxy;
        aph->location = chip.nampx*(yL - miny) + (xL - minx);
    }
    if (impurityvariation == 1) {
        dopant = silicon.nbulkmap[aph->location]*nbulk;
    } else {
        dopant = nbulk;
    }
    windex = find_linear(silicon.dopantGrid, silicon.numDopant, dopant, &rwindex);
    double sg = interpolate_trilinear(silicon.sigma, silicon.numTemperature, silicon.numThickness,
                                      windex, rwindex, aph->xindex, aph->rxindex, zindex, rzindex);
    if (chargediffusion == 0) sg = 0.0;

    // complications to charge diffusion (effect of lateral fields)
    double fsg = interpolate_bilinear_float(silicon.fsigma, silicon.numThickness, windex, rwindex, zindex, rzindex);
    double gsg = interpolate_bilinear_float(silicon.gsigma, silicon.numThickness, windex, rwindex, zindex, rzindex);
    double sa, sb, ga, gb, da, db;
    double chsx, chsy, chv, chp;
    sa = silicon.sigmaX[aph->location];
    sb = silicon.sigmaY[aph->location];
    ga = silicon.gammaX[aph->location];
    gb = silicon.gammaY[aph->location];
    da = silicon.deltaX[aph->location];
    db = silicon.deltaY[aph->location];
    if (chargesharing == 1) {
        chsx = 0.0;
        chsy = 0.0;
        double rho, rhoprime, cost, sint;
        if (aph->xPos <= maxx && aph->xPos >= minx && aph->yPos >= miny && aph->yPos <= maxy) {
            rho = sqrt((aph->xPosR)*(aph->xPosR) + (aph->yPosR)*(aph->yPosR));
            cost = aph->xPosR/rho;
            sint = aph->yPosR/rho;
            rhoprime = sqrt(rho*rho + cost*cost*silicon.spaceChargeSpreadX + sint*sint*silicon.spaceChargeSpreadY);
            uindex = find_linear(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-4, &ruindex);
            chp = *(state.focal_plane + chip.nampx*(aph->yPos - miny) + (aph->xPos - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*cost;
            chsy += chp*chv*sint;
        }
        if (aph->xPos + 1 <= maxx && aph->xPos + 1 >= minx && aph->yPos >= miny && aph->yPos <= maxy) {
            rho = sqrt((aph->xPosR - 1)*(aph->xPosR - 1) + (aph->yPosR)*(aph->yPosR));
            cost = (aph->xPosR - 1)/rho;
            sint = aph->yPosR/rho;
            rhoprime = sqrt(rho*rho + cost*cost*silicon.spaceChargeSpreadX + sint*sint*silicon.spaceChargeSpreadY);
            uindex = find_linear(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-4, &ruindex);
            chp = *(state.focal_plane + chip.nampx*(aph->yPos - miny) + (aph->xPos + 1 - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*cost;
            chsy += chp*chv*sint;
        }
        if (aph->xPos - 1 <= maxx && aph->xPos - 1 >= minx && aph->yPos >= miny && aph->yPos <= maxy) {
            rho = sqrt((aph->xPosR + 1)*(aph->xPosR + 1) + (aph->yPosR)*(aph->yPosR));
            cost = (aph->xPosR + 1)/rho;
            sint = aph->yPosR/rho;
            rhoprime = sqrt(rho*rho + cost*cost*silicon.spaceChargeSpreadX + sint*sint*silicon.spaceChargeSpreadY);
            uindex = find_linear(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-4, &ruindex);
            chp = *(state.focal_plane + chip.nampx*(aph->yPos - miny) + (aph->xPos - 1 - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*cost;
            chsy += chp*chv*sint;
        }
        if (aph->xPos <= maxx && aph->xPos >= minx && aph->yPos + 1 >= miny && aph->yPos + 1 <= maxy) {
            rho = sqrt((aph->xPosR)*(aph->xPosR) + (aph->yPosR - 1)*(aph->yPosR - 1));
            cost = aph->xPosR/rho;
            sint = (aph->yPosR - 1)/rho;
            rhoprime = sqrt(rho*rho + cost*cost*silicon.spaceChargeSpreadX + sint*sint*silicon.spaceChargeSpreadY);
            uindex = find_linear(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-4, &ruindex);
            chp = *(state.focal_plane + chip.nampx*(aph->yPos + 1 - miny) + (aph->xPos - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*cost;
            chsy += chp*chv*sint;
        }
        if (aph->xPos <= maxx && aph->xPos >= minx && aph->yPos - 1 >= miny && aph->yPos - 1 <= maxy) {
            rho = sqrt((aph->xPosR)*(aph->xPosR) + (aph->yPosR + 1)*(aph->yPosR + 1));
            cost = aph->xPosR/rho;
            sint = (aph->yPosR + 1)/rho;
            rhoprime = sqrt(rho*rho + cost*cost*silicon.spaceChargeSpreadX + sint*sint*silicon.spaceChargeSpreadY);
            uindex = find_linear(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-4, &ruindex);
            chp = *(state.focal_plane + chip.nampx*(aph->yPos - 1 - miny) + (aph->xPos - minx));
            chv = interpolate_trilinear(silicon.hsigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsx += chp*chv*cost;
            chsy += chp*chv*sint;
        }
        if (aph->xPos <= maxx && aph->xPos >= minx && aph->yPos >= miny && aph->yPos - 1 <= maxy) {
            rho = fabs(aph->yPosR + 0.5);
            uindex = find_linear(silicon.rho, silicon.numTemperature, rho*pixsize*1e-4, &ruindex);
            chp = silicon.chargeStopCharge;
            chv = interpolate_trilinear(silicon.isigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsy += chp*chv*(aph->yPosR + 0.5)/rho;
        }
        if (aph->xPos <= maxx && aph->xPos >= minx && aph->yPos + 1 >= miny && aph->yPos <= maxy) {
            rho = fabs(aph->yPosR - 0.5);
            uindex = find_linear(silicon.rho, silicon.numTemperature, rho*pixsize*1e-4, &ruindex);
            chp = silicon.chargeStopCharge;
            chv = interpolate_trilinear(silicon.isigma, silicon.numTemperature, silicon.numThickness,
                                        windex, rwindex, uindex, ruindex, zindex, rzindex);
            chsy += chp*chv*(aph->yPosR - 0.5)/rho;
        }
    } else {
        chsx = 0.0;
        chsy = 0.0;
    }
    if (fieldanisotropy == 0) {
        sa = 0.0;
        sb = 0.0;
        da = 0.0;
        db = 0.0;
    }
    if (pixelerror == 0) {
        ga = 0.0;
        gb = 0.0;
    }
    double nr = static_cast<double>(SILICON_STEPS)/static_cast<double>(SILICON_SUB_STEPS);
    position->x += (sg*random[aph->thread].normal() + fsg*sa + gsg*da + ga + chsx)/1e3;
    position->y += (sg*random[aph->thread].normal() + fsg*sb + gsg*db + gb + chsy)/1e3;
    if (position->z <= aph->collect_z) {
        position->z += sensorthickness/1e3/nr;
    } else {
        position->z -= sensorthickness/1e3/nr;
    }
    return(0);

}

double Image::fringing (Vector angle, Vector normal, double wavelength, double nSi, double thickness, double meanFreePath, int polarization) {

    // double arg = fabs(normal.x*angle.x + normal.y*angle.y + normal.z*angle.z);
    // if (arg > 1) arg = 1.0;
    // double airAngle = acos(arg);
    // double siliconAngle = asin(sin(airAngle)/nSi);
    // double airAngleOut = -airAngle;
    // double delta = 2*PI/wavelength*thickness*nSi*sqrt(1.0 - sin(airAngle)*sin(airAngle)/(nSi*nSi));
    // std::complex<double> rhot1 ((pow(cos(airAngle), polarization) - nSi*pow(cos(siliconAngle), polarization))/
    //                             (pow(cos(airAngle), polarization) + nSi*pow(cos(siliconAngle), polarization)), 0.0);
    // std::complex<double> rhot2 ((nSi*pow(cos(siliconAngle), polarization) - pow(cos(airAngleOut), polarization))/
    //                             (pow(cos(airAngleOut), polarization) + nSi*pow(cos(siliconAngle), polarization)), 0.0);
    // std::complex<double> phase (cos(-2.0*delta), sin(-2.0*delta));
    // std::complex<double> gamma;
    // std::complex<double> one (1, 0);
    // gamma = (rhot1 + rhot2*phase)/(one + rhot1*rhot2*phase);
    // double reflection = real(gamma*conj(gamma));
    // return(reflection);

    std::complex<double> nSilicon (nSi, -wavelength/meanFreePath/1e3/4/PI);
    std::complex<double> nOx (1.5, 0);
    std::complex<double> one (1, 0);
    std::complex<double> imaginary (0, 1);
    double thicknessOx = 2.0/1e3;
    double thicknessSi = thickness;

    double arg = fabs(normal.x*angle.x + normal.y*angle.y + normal.z*angle.z);
    if (arg > 1) arg = 1.0;
    double airAngle = acos(arg);
    double siliconAngle = asin(sin(airAngle)/real(nSilicon));
    double oxideAngle = asin(sin(airAngle)/real(nOx));
    double airAngleOut = -airAngle;

    std::complex<double> factor1 (sqrt(1.0 - sin(airAngle)*sin(airAngle)/(real(nSilicon)*real(nSilicon))), 0);
    std::complex<double> factor2 (sqrt(1.0 - sin(airAngle)*sin(airAngle)/(real(nOx)*real(nOx))), 0);
    std::complex<double> delta1;
    std::complex<double> delta2;

    delta1 = 2*PI/wavelength*thicknessSi*nSilicon*factor1;
    delta2 = 2*PI/wavelength*thicknessOx*nOx*factor2;

    std::complex<double> rhot1;
    std::complex<double> rhot2;
    std::complex<double> rhot3;
    rhot1 = (pow(cos(airAngle), polarization) - nSilicon*pow(cos(siliconAngle), polarization))/
            (pow(cos(airAngle), polarization) + nSilicon*pow(cos(siliconAngle), polarization));
    rhot2 = (nSilicon*pow(cos(siliconAngle), polarization) - nOx*pow(cos(oxideAngle), polarization))/
            (nSilicon*pow(cos(siliconAngle), polarization) + nOx*pow(cos(oxideAngle), polarization));
    rhot3 = (nOx*pow(cos(oxideAngle), polarization) - pow(cos(airAngleOut), polarization))/
            (nOx*pow(cos(oxideAngle), polarization) + pow(cos(airAngleOut), polarization));

    std::complex<double> phase1;
    std::complex<double> phase2;
    phase1 = exp(-2.0*imaginary*delta1);
    phase2 = exp(-2.0*imaginary*delta2);

    std::complex<double> gamma1;
    std::complex<double> gamma2;
    gamma2 = (rhot2 + rhot3*phase2)/(one + rhot2*rhot3*phase2);
    gamma1 = (rhot1 + gamma2*phase1)/(one + rhot1*gamma2*phase1);

    double convert = real(gamma1*conj(gamma1)) + real(phase1*conj(phase1))*real(phase2*conj(phase2))*(1.0 - real(rhot3*conj(rhot3)));
    return(convert);


}


int Image::contaminationSurfaceCheck(Vector position, Vector *angle, long surfaceIndex, Photon *aph) {

    if (random[aph->thread].uniform() > (double)(*(contamination.transmission +
                                      PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
                                      aph->vvint*PERTURBATION_POINTS + aph->uuint))) {
        return(1);
    }

    if (*(contamination.surfacelistmap + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
          aph->vvint*PERTURBATION_POINTS + aph->uuint) != -1) {
        long cc = *(contamination.surfacelistmap + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
                    aph->vvint*PERTURBATION_POINTS + aph->uuint);

        if (sqrt(pow(position.x - contamination.surfacelistx[surfaceIndex][cc], 2.0)+
                 pow(position.y - contamination.surfacelisty[surfaceIndex][cc], 2.0)) <
            contamination.surfacelists[surfaceIndex][cc]) {

            if (random[aph->thread].uniform() > exp(-contamination.absorptionLength*contamination.surfacelists[surfaceIndex][cc])) return(1);
            long index;
            double mu;

            if (contamination.surfacelistt[surfaceIndex][cc] == 0) {
                find(contamination.henyey_greenstein, contamination.elements, random[aph->thread].uniform(), &index);
                mu = contamination.henyey_greenstein_mu[index];
            } else {
                find(contamination.henyey_greenstein_w, contamination.elements, random[aph->thread].uniform(), &index);
                mu = contamination.henyey_greenstein_mu_w[index];
            }
            double phi = random[aph->thread].uniform()*2.0*PI;
            shift_mu(angle, mu, phi);
            if (mu < 0) {
                return(1);
            }
            aph->ghostFlag = 1;
        }
    }
    return(0);


}
