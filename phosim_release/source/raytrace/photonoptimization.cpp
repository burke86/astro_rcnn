///
/// @package phosim
/// @file photonoptimization.cpp
/// @brief photon optimization routines (part of Image class)
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

int Image::dynamicTransmissionOptimization(long k, long *lastSurface, long *preGhost, long waveSurfaceIndex, long straylightcurrent, Photon *aph) {

    long newSurf;
    int status;

 redoghlp:;
    newSurf = -1;
    aph->direction = 1;
    aph->counter = -1;
    *preGhost = 0;

    for (int i = -1; i < natmospherefile; i++) {

        if (i >= 0) {
            if (aph->counter >= (MAX_BOUNCE - 1)) goto maxbounce;
            status = transmissionPreCheck(i*2 + 1, waveSurfaceIndex, aph);
            if (status == 1) {
                if (k > 0) goto redoghlp;
                return(1);
            }
            // if (status == 2) {
            //     (*preGhost)++;
            //     goto maxbounce;
            // }
        }

        if (aph->counter >= (MAX_BOUNCE - 1)) goto maxbounce;
        status = transmissionPreCheck(i*2 + 2, waveSurfaceIndex, aph);
        if (status == 1) {
            if (k > 0) goto redoghlp;
            return(1);
        }
        // if (status == 2) {
        //     (*preGhost)++;
        //     goto maxbounce;
        // }
    }
    while (1) {
        if (aph->direction == 1) {
            newSurf++;
        } else {
            newSurf--;
        }
        if (newSurf == -1) {
            if (k > 0) goto redoghlp;
            return(1);
        }

        if (throughputfile == 1 && aph->direction == 1) {
            *lastSurface = newSurf;
        }
        if (aph->counter >= (MAX_BOUNCE - 1)) goto maxbounce;
        status = transmissionPreCheck(natmospherefile*2 + 1 + 2*newSurf, waveSurfaceIndex, aph);
        if (status == 1) {
            if (straylightcurrent == 1) {
                if (aph->counter >= (MAX_BOUNCE - 1)) goto maxbounce;
                status = transmissionPreCheck(natmospherefile*2 + 1 + 2*newSurf + 1, waveSurfaceIndex, aph);
                if (status == 1) {
                    if (k > 0) goto redoghlp;
                    return(1);
                } else if (status == 0) {
                    aph->direction = -aph->direction;
                    (*preGhost)++;
                } else if (status == 2) {
                    (*preGhost)++;
                    goto maxbounce;
                }
            } else {
                if (k > 0) goto redoghlp;
                return(1);
            }
        }
        if (aph->direction == 1 && newSurf == nsurf - 1) {
            if (aph->counter >= (MAX_BOUNCE - 1)) goto maxbounce;
            status = transmissionPreCheck(natmospherefile*2 + 1 + 2*nsurf, waveSurfaceIndex, aph);
            if (status == 1) {
                if (k > 0) goto redoghlp;
                return(1);
            } else if (status == 0) {
                goto maxbounce;
            } else if (status == 2) {
                (*preGhost)++;
                goto maxbounce;
            }
        }

        if (aph->direction == -1 && newSurf == 0) {
            if (k > 0) goto redoghlp;
            return(1);
        }
    }
 maxbounce:;
    aph->maxcounter = aph->counter;

    return(0);

}
