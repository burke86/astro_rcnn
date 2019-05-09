///
/// @package phosim
/// @file helpers.h
/// @brief helpers
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <vector>

inline void find (double *y, long n, double x, long *index) {
    long high, mid, low;

    low = 0;
    high = n - 1;
    while (high - low > 1) {
        mid = (high + low) >> 1;
        if (x > y[mid]) {
            low = mid;
        } else {
            high = mid;
        }
    }
     *index = low;
 }
inline long find_linear (double *xx, long n, double x, double *dj) {
    long j;
    *dj = (x - xx[0])/(xx[n - 1] - xx[0])*(n - 1);
    j = (long)(*dj);
    if (j > n - 2) j = n - 2;
    if (j < 0) j = 0;
    return j;
 }
inline long find_linear_v (std::vector<double> & xx, double x, double *dj) {
    long j;
    long n = xx.size();

    *dj = (x - xx[0])/(xx[n - 1] - xx[0])*(n - 1);
    j = (long)(*dj);
    if (j > n - 2) j = n - 2;
    if (j < 0) j = 0;

    return j;
}

inline void find_linear_wrap (double x, double dx, long n, long *i, long *j, double *di) {

    double xx, fxx;

    xx = x/dx + (static_cast<double>(n)/2 - 0.5);
    fxx = floor(xx);

    *i = static_cast<long>(fxx);
    *j = (*i) + 1;

    *i = (*i) % n;
    if ((*i) < 0) (*i) += n;

    *j = (*j) % n;
    if ((*j) < 0) (*j) += n;

    *di = xx - fxx;

}


inline void normalize (double *vx, double *vy, double *vz) {

    double r;

    r = sqrt(((*vx)*(*vx)) + ((*vy)*(*vy)) + ((*vz)*(*vz)));

    *vx = *vx/r;
    *vy = *vy/r;
    *vz = *vz/r;

}

inline double interpolate (double *yy, double *xx, double x, long index) {

    return (yy[index]*(xx[index + 1] - x) + yy[index + 1]*(x - xx[index]))/(xx[index + 1] - xx[index]);

}

inline double interpolate_linear (double *yy, long index,  double rindex) {

    double dr = rindex - index;
    return (yy[index]*(1.0 - dr) + yy[index + 1]*dr);

}


inline double interpolate_bilinear (double *zz, long nelem, long xindex, double rxindex, long yindex, double ryindex) {

    double dx, dy;

    dx = rxindex - xindex;
    dy = ryindex - yindex;

    return (zz[xindex*nelem + yindex]*(1 - dx)*(1 - dy)+
            zz[(xindex + 1)*nelem + yindex]*(dx)*(1 - dy)+
            zz[xindex*nelem + yindex + 1]*(1 - dx)*(dy)+
            zz[(xindex + 1)*nelem + yindex + 1]*(dx)*(dy));


}

inline double interpolate_bilinear_float (float *zz, long nelem, long xindex, double rxindex, long yindex, double ryindex) {

    double dx, dy;

    dx = rxindex - xindex;
    dy = ryindex - yindex;

    return (zz[xindex*nelem + yindex]*(1 - dx)*(1 - dy)+
            zz[(xindex + 1)*nelem + yindex]*(dx)*(1 - dy)+
            zz[xindex*nelem + yindex + 1]*(1 - dx)*(dy) +
            zz[(xindex + 1)*nelem + yindex + 1]*(dx)*(dy));


}

inline float interpolate_trilinear (float *zz, long nelemy, long nelemz, long xindex, double rxindex,
                                     long yindex, double ryindex, long zindex, double rzindex ) {

    double dx, dy, dz;

    dx = rxindex - xindex;
    dy = ryindex - yindex;
    dz = rzindex - zindex;

    return (zz[xindex*nelemy*nelemz + yindex*nelemz + zindex]*(1 - dx)*(1 - dy)*(1 - dz) +
            zz[(xindex + 1)*nelemy*nelemz + yindex*nelemz + zindex]*(dx)*(1 - dy)*(1 - dz) +
            zz[xindex*nelemy*nelemz + (yindex + 1)*nelemz + zindex]*(1 - dx)*(dy)*(1 - dz) +
            zz[(xindex + 1)*nelemy*nelemz + (yindex + 1)*nelemz + zindex]*(dx)*(dy)*(1 - dz) +
            zz[xindex*nelemy*nelemz + yindex*nelemz + (zindex + 1)]*(1 - dx)*(1 - dy)*(dz) +
            zz[(xindex + 1)*nelemy*nelemz + yindex*nelemz + (zindex + 1)]*(dx)*(1 - dy)*(dz) +
            zz[xindex*nelemy*nelemz + (yindex + 1)*nelemz + (zindex + 1)]*(1 - dx)*(dy)*(dz) +
            zz[(xindex + 1)*nelemy*nelemz + (yindex + 1)*nelemz + (zindex + 1)]*(dx)*(dy)*(dz));

}

inline double interpolate_bilinear_float_wrap (float *zz, long nelem, long xindex, long xindex1, double dx, long yindex, long yindex1, double dy ) {

    double dxdy = dx*dy;

    return (zz[xindex*nelem  + yindex ]*(1 - dx - dy + dxdy) +
            zz[xindex1*nelem + yindex ]*(dx - dxdy) +
            zz[xindex*nelem  + yindex1]*(dy - dxdy) +
            zz[xindex1*nelem + yindex1]*(dxdy));


}



// Interpolation search
// asserts that: sorted_array[0] <= value
// returns position of first element less or equal to parameter value
inline int interpolationSearch(double sorted_array[], long n, double value) {
    int beg = 0;
    int end = n - 1;

    double begVal = sorted_array[beg];
    double endVal = sorted_array[end];

    if (endVal <= value) {
        return end;
    }
    //invariant for following loop: endVal>value

    while (begVal < endVal) {
        int    testPos = beg + int((value - begVal)*(end - beg)/(endVal - begVal));
        double testVal = sorted_array[testPos];

        if (value < testVal) {
            end = testPos;
            endVal = testVal;
        } else {
            if (testPos > beg){
                beg = testPos;
                begVal = testVal;
            } else {
                beg++;
                begVal = sorted_array[beg];
                if (begVal > value) {
                    return beg - 1;
                }
            }
        }
    }
    return beg;
}

inline double smallAnglePupilNormalize (double x, double y) {

    return(-sqrt(1.0 - x*x - y*y));
}


#ifdef ISO_LIBRARIES_ONLY
    inline bool isnan(double d) { return d != d; }
#endif
