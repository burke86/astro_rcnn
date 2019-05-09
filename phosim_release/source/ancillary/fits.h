///
/// @package phosim
/// @file fits.h
/// @brief fits helpers
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <fitsio.h>
#include <fitsio2.h>


inline void fitsWriteKey(fitsfile *fptr, const char *keyname, float value, const char *comment) {
    int status = 0;
    fits_write_key(fptr, TFLOAT, keyname, &value, comment, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsWriteKey(fitsfile *fptr, const char *keyname, double value, const char *comment) {
    int status = 0;
    fits_write_key(fptr, TDOUBLE, keyname, &value, comment, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsWriteKey(fitsfile *fptr, const char *keyname, long value, const char *comment) {
    int status = 0;
    fits_write_key(fptr, TLONG, keyname, &value, comment, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}


inline void fitsWriteKey(fitsfile *fptr, const char *keyname, int value, const char *comment) {
    int status = 0;
    fits_write_key(fptr, TINT, keyname, &value, comment, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsWriteKey(fitsfile *fptr, const char *keyname, const char *value, const char *comment) {
    int status = 0;
    char tempstring[4096];
    sprintf(tempstring, "%s", value);
    if (strlen(tempstring) < 69) {
        fits_write_key(fptr, TSTRING, keyname, tempstring, comment, &status);
    } else {
        fits_write_key_longstr(fptr, keyname, tempstring, comment, &status);
    }
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsWriteKey(fitsfile *fptr, const char *keyname, std::string value, const char *comment) {
    int status = 0;
    char tempstring[4096];
    sprintf(tempstring, "%s", value.c_str());
    if (strlen(tempstring) < 69) {
        fits_write_key(fptr, TSTRING, keyname, tempstring, comment, &status);
    } else {
        fits_write_key_longstr(fptr, keyname, tempstring, comment, &status);
    }
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsUpdateKey(fitsfile *fptr, const char *keyname, float value, char *comment) {
    int status = 0;
    char tempstring[4096];
    sprintf(tempstring, "%s", comment);
    fits_update_key(fptr, TFLOAT, keyname, &value, tempstring, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsUpdateKey(fitsfile *fptr, const char *keyname, double value, char *comment) {
    int status = 0;
    char tempstring[4096];
    sprintf(tempstring, "%s", comment);
    fits_update_key(fptr, TDOUBLE, keyname, &value, tempstring, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsUpdateKey(fitsfile *fptr, const char *keyname, int value, const char *comment) {
    int status = 0;
    char tempstring[4096];
    sprintf(tempstring, "%s", comment);
    fits_update_key(fptr, TSHORT, keyname, &value, tempstring, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}


inline void fitsUpdateKey(fitsfile *fptr, const char *keyname, long value, const char *comment) {
    int status = 0;
    char tempstring[4096];
    sprintf(tempstring, "%s", comment);
    fits_update_key(fptr, TLONG, keyname, &value, tempstring, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsUpdateKey(fitsfile *fptr, const char *keyname, unsigned long value, const char *comment) {
    int status = 0;
    char tempstring[4096];
    sprintf(tempstring, "%s", comment);
    fits_update_key(fptr, TULONG, keyname, &value, tempstring, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsUpdateKey(fitsfile *fptr, const char *keyname, unsigned short value, char *comment) {
    int status = 0;
    char tempstring[4096];
    sprintf(tempstring, "%s", comment);
    fits_update_key(fptr, TUSHORT, keyname, &value, tempstring, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsUpdateKey(fitsfile *fptr, const char *keyname, const char *value, char *comment) {
    int status = 0;
    char tempstring[4096];
    sprintf(tempstring, "%s", value);
    if (strlen(tempstring) < 69) {
        fits_update_key(fptr, TSTRING, keyname, tempstring, comment, &status);
    } else {
        fits_update_key_longstr(fptr, keyname, tempstring, comment, &status);
    }
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsReadKey(fitsfile *fptr, const char *keyname, float *value, char *comment) {
    int status = 0;
    fits_read_key(fptr, TFLOAT, keyname, value, comment, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsReadKey(fitsfile *fptr, const char *keyname, double *value, char *comment) {
    int status = 0;
    fits_read_key(fptr, TDOUBLE, keyname, value, comment, &status);
    if (status != 0) {
        printf("FITS Header Key Error: %d\n", status);
        exit(1);
    }
}

inline void fitsCreateImage(fitsfile **fptr, std::string filename) {

    int status = 0;
    long naxes[2];

    naxes[0] = 1;
    naxes[1] = 1;
    fits_create_file(fptr, filename.c_str(), &status);
    if (status != 0) {
        printf("FITS Create FILE Error: %d\n", status);
        exit(1);
    }
    fits_create_img(*fptr, FLOAT_IMG, 2, naxes, &status);
    if (status != 0) {
        printf("FITS Create Image Error: %d\n", status);
        exit(1);
    }

}


inline void fitsCreateImage(fitsfile **fptr, const char *filename) {

    int status = 0;
    long naxes[2];

    naxes[0] = 1;
    naxes[1] = 1;
    fits_create_file(fptr, filename, &status);
    if (status != 0) {
        printf("FITS Create FILE Error: %d\n", status);
        exit(1);
    }
    fits_create_img(*fptr, FLOAT_IMG, 2, naxes, &status);
    if (status != 0) {
        printf("FITS Create Image Error: %d\n", status);
        exit(1);
    }

}
inline void fitsCreateImageLong(fitsfile **fptr, const char *filename) {

    int status = 0;
    long naxes[2];

    naxes[0] = 1;
    naxes[1] = 1;
    fits_create_file(fptr, filename, &status);
    if (status != 0) {
        printf("FITS Create FILE Error: %d\n", status);
        exit(1);
    }
    fits_create_img(*fptr, LONG_IMG, 2, naxes, &status);
    if (status != 0) {
        printf("FITS Create Image Error: %d\n", status);
        exit(1);
    }

}
inline void fitsCreateImageShort(fitsfile **fptr, const char *filename) {

    int status = 0;
    long naxes[2];

    naxes[0] = 1;
    naxes[1] = 1;
    fits_create_file(fptr, filename, &status);
    if (status != 0) {
        printf("FITS Create FILE Error: %d\n", status);
        exit(1);
    }
    fits_create_img(*fptr, SHORT_IMG, 2, naxes, &status);
    if (status != 0) {
        printf("FITS Create Image Error: %d\n", status);
        exit(1);
    }

}

inline void fitsWriteImage(fitsfile *fptr, long x, long y, float *image) {

    int status = 0;
    long naxes[2];

    naxes[0] = x;
    naxes[1] = y;
    fits_update_key(fptr, TLONG, (char*)"NAXIS1", &naxes[0], NULL, &status);
    if (status != 0) {
        printf("FITS Create Image Error: %d\n", status);
        exit(1);
    }
    fits_update_key(fptr, TLONG, (char*)"NAXIS2", &naxes[1], NULL, &status);
    if (status != 0) {
        printf("FITS Update Key Error: %d\n", status);
        exit(1);
    }
    fits_write_img(fptr, TFLOAT, 1, x*y, image, &status);
    if (status != 0) {
        printf("FITS Update Key Error: %d\n", status);
        exit(1);
    }
    fits_close_file(fptr, &status);
    if (status != 0) {
        printf("FITS Close Image Error: %d\n", status);
        exit(1);
    }


}

inline void fitsWriteImage(fitsfile *fptr, long x, long y, double *image) {

    int status = 0;
    long naxes[2];

    naxes[0] = x;
    naxes[1] = y;
    fits_update_key(fptr, TLONG, (char*)"NAXIS1", &naxes[0], NULL, &status);
    if (status != 0) {
        printf("FITS Create Image Error: %d\n", status);
        exit(1);
    }
    fits_update_key(fptr, TLONG, (char*)"NAXIS2", &naxes[1], NULL, &status);
    if (status != 0) {
        printf("FITS Update Key Error: %d\n", status);
        exit(1);
    }
    fits_write_img(fptr, TDOUBLE, 1, x*y, image, &status);
    if (status != 0) {
        printf("FITS Update Key Error: %d\n", status);
        exit(1);
    }
    fits_close_file(fptr, &status);
    if (status != 0) {
        printf("FITS Close Image Error: %d\n", status);
        exit(1);
    }


}

inline void fitsWriteImage(fitsfile *fptr, long x, long y, unsigned long *image) {

    int status = 0;
    long naxes[2];

    naxes[0] = x;
    naxes[1] = y;
    fits_update_key(fptr, TLONG, (char*)"NAXIS1", &naxes[0], NULL, &status);
    if (status != 0) {
        printf("FITS Update Key Error: %d\n", status);
        exit(1);
    }
    fits_update_key(fptr, TLONG, (char*)"NAXIS2", &naxes[1], NULL, &status);
    if (status != 0) {
        printf("FITS Update Key Error: %d\n", status);
        exit(1);
    }
    fits_write_img(fptr, TULONG, 1, x*y, image, &status);
    if (status != 0) {
        printf("FITS Write Image Error: %d\n", status);
        status = 0;
        /* exit(1); */
    }
    fits_close_file(fptr, &status);
    if (status != 0) {
        printf("FITS Close Image Error: %d\n", status);
        exit(1);
    }

}

inline void fitsWriteImage(fitsfile *fptr, long x, long y, unsigned short *image) {

    int status = 0;
    long naxes[2];

    naxes[0] = x;
    naxes[1] = y;
    fits_update_key(fptr, TLONG, (char*)"NAXIS1", &naxes[0], NULL, &status);
    if (status != 0) {
        printf("FITS Update Key Error: %d\n", status);
        exit(1);
    }
    fits_update_key(fptr, TLONG, (char*)"NAXIS2", &naxes[1], NULL, &status);
    if (status != 0) {
        printf("FITS Update Key Error: %d\n", status);
        exit(1);
    }
    fits_write_img(fptr, TUSHORT, 1, x*y, image, &status);
    if (status != 0) {
        printf("FITS Write Image Error: %d\n", status);
        exit(1);
    }
    fits_close_file(fptr, &status);
    if (status != 0) {
        printf("FITS Close Image Error: %d\n", status);
        exit(1);
    }


}


inline void fitsReadImage(std::string filename, float *image) {

    fitsfile *faptr;
    long naxes[2];
    int nfound = 0;
    int anynull = 0;
    float nullval = 0.0;
    int status = 0;

    if (fits_open_file(&faptr, filename.c_str(), READONLY, &status)) {
        printf("Error opening %s\n", filename.c_str());
        exit(1);
    }
    fits_read_keys_lng(faptr, (char*)"NAXIS", 1, 2, naxes, &nfound, &status);
    fits_read_img(faptr, TFLOAT, 1, naxes[0]*naxes[1], &nullval, image, &anynull, &status);
    fits_close_file(faptr, &status);

}
