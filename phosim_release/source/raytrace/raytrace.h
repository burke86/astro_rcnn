#include "vector_type.h"

void vectorCopy(Vector vectorIn, Vector *vectorOut);
void setup_tangent (double ra, double dec, Vector *tpx, Vector *tpy, Vector *tpz);
void tangent (double ra, double dec, double *x, double *y, Vector *tpx, Vector *tpy, Vector *tpz);
double vectorDot( Vector *vector1, Vector *vector2);
void vectorAdd(Vector vectorA, Vector vectorB, Vector *vectorOut);
void vectorInit (Vector *vector);
void vectorSubtract( Vector vectorA, Vector vectorB, Vector *vectorOut);
void normalize (Vector *vector);
double modulus (Vector *vector);
void propagate(Vector *position, Vector angle, double distance);
void reflect(Vector *vectorIn, Vector normal);
void refract (Vector *vectorIn, Vector normal, double n_1, double n_2);
double rotateInverseX(Vector *vector, double angle);
double rotateInverseY(Vector *vector, double angle);
void shift_mu(Vector *vector, double mu, double phi);
void zernikes(double *zernike_r, double *zernike_phi, double *zernike_r_grid,
              double *zernike_phi_grid, double *zernike_normal_r, double *zernike_normal_phi,
              long numelements, long nzernikes);
double chebyshevT(int n, double x);
void chebyshevT2D(int n, double x, double y, double *t);
void chebyshevs(double *r_grid, double *phi_grid, double *chebyshev, double *chebyshev_r, double *chebyshev_phi, long nPoint, long nTerm);
