///
/// @package phosim
/// @file screen.cpp
/// @brief screen class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "screen.h"

double Screen::readScreen (int keynum, float *array, char *tempstring) {


        char *ffptr;
        fitsfile *faptr;
        long naxes[2];
        int nfound;
        int anynull;
        float nullval;
        char keyname[4096];
        char comment[4096];
        char value[4096];
        int status;

        ffptr=tempstring;
        status=0;
        if (fits_open_file(&faptr,ffptr,READONLY,&status)) {
            printf("Error opening %s\n",ffptr);
            exit(1);
        }
        fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
        if (keynum >= 0) fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
        fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,array,&anynull,&status);
        fits_close_file(faptr,&status);

        return(atof(value));

}
