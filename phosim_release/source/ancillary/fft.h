///
/// @package phosim
/// @file fft.h
/// @brief helper functions for fftw3
///
/// @brief Created by
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "fftw3.h"

inline void fftMalloc(fftw_complex **array, long n) {

    *array = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n);
    if (*array == NULL) {
        printf("fftw3 malloc failed.\n");
        exit(1);
    }
    for (long i = 0; i < n; i++) {
        (*array)[i][0] = 0.0;
        (*array)[i][1] = 0.0;
    }

}

inline void fftMalloc(double **array, long n) {

    *array = static_cast<double*>(fftw_malloc(sizeof(double)*n));
    if (*array == NULL) {
        printf("fftw3 malloc failed.\n");
        exit(1);
    }
    for (long i = 0; i < n; i++) {
        (*array)[i] = 0.0;
    }

}


inline void inverseFFT (long nx, long ny, fftw_complex *out, fftw_complex *in) {

    fftw_plan plan;

    plan = fftw_plan_dft_2d(nx, ny, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

}
