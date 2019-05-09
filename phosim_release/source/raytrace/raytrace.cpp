///
/// @package phosim
/// @file raytrace.cpp
/// @brief raytrace functions
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author J. Garrett Jernigan (UCB)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <math.h>
#include "raytrace.h"
#include "constants.h"
#include "helpers.h"

void vectorCopy (Vector vectorIn, Vector *vectorOut) {

    vectorOut->x = vectorIn.x;
    vectorOut->y = vectorIn.y;
    vectorOut->z = vectorIn.z;

}

void vectorAdd (Vector vectorA, Vector vectorB, Vector *vectorOut) {

    vectorOut->x = vectorA.x + vectorB.x;
    vectorOut->y = vectorA.y + vectorB.y;
    vectorOut->z = vectorA.z + vectorB.z;

}

void vectorSubtract (Vector vectorA, Vector vectorB, Vector *vectorOut) {

    vectorOut->x = vectorA.x-vectorB.x;
    vectorOut->y = vectorA.y - vectorB.y;
    vectorOut->z = vectorA.z - vectorB.z;

}

void vectorInit (Vector *vector) {

    vector->x = 0.0;
    vector->y = 0.0;
    vector->z = 0.0;

}

void dcrsph (Vector *vector, double longitude, double latitude) {

    vector->x = cos(longitude)*cos(latitude);
    vector->y = sin(longitude)*cos(latitude);
    vector->z = sin(latitude);

}

void dcarsph (Vector *vector, double *longitude, double *latitude) {

    *latitude = asin(vector->z);
    *longitude = atan2(vector->y, vector->x);
    if (*longitude<0.0) *longitude += 2.0*PI;

}

double vectorDot(Vector *vectorA, Vector *vectorB) {

    return(vectorA->x*vectorB->x + vectorA->y*vectorB->y + vectorA->z*vectorB->z);
}

void vectorCross(Vector *vectorA, Vector *vectorB, Vector *vectorC) {

    vectorC->x = vectorA->y*vectorB->z - vectorA->z*vectorB->y;
    vectorC->y = vectorA->z*vectorB->x - vectorA->x*vectorB->z;
    vectorC->z = vectorA->x*vectorB->y - vectorA->y*vectorB->x;
}


void setup_tangent ( double ra, double dec, Vector *tpx, Vector *tpy, Vector *tpz) {

    dcrsph(tpz, ra, dec);
    dcrsph(tpx, ra + (PI/2.0), 0.0);
    vectorCross(tpz, tpx, tpy);
}

void tangent (double ra, double dec, double *x, double *y, Vector *tpx, Vector *tpy, Vector *tpz) {

    double dotz;
    Vector s;

    dcrsph(&s, ra, dec);
    dotz = vectorDot(&s, tpz);
    *x = vectorDot(&s, tpx)/dotz;
    *y = vectorDot(&s, tpy)/dotz;

    // s.x = tpz->x  +   (*x)*tpx->x   +  (*y)*tpy->x;
    // s.y = tpz->y  +   (*x)*tpx->y   +  (*y)*tpy->y;
    // s.z = tpz->z  +   (*x)*tpx->z   +  (*y)*tpy->z;
    // normalize(&S);

}

void shift_mu (Vector *vector, double mu, double phi) {

    double vx = vector->x;
    double vy = vector->y;
    double vz = vector->z;
    if (vx==0.0 && vy==0.0) vy = vy + 1e-9;
    double vr = sqrt(vx*vx + vy*vy);

    double x = sqrt(1 - mu*mu)*cos(phi);
    double y = sqrt(1 - mu*mu)*sin(phi);

    vector->x = (-vy)/vr*x + vz*vx/vr*y + vx*mu;
    vector->y = vx/vr*x + vz*vy/vr*y + vy*mu;
    vector->z = (-vr)*y + vz*mu;

}

void normalize (Vector *vector) {

    double r;

    r = modulus(vector);
    vector->x /= r;
    vector->y /= r;
    vector->z /= r;

}

double rotateInverseX(Vector *vector, double angle) {

    return(vector->x*cos(angle) + vector->y*sin(angle));

}

double rotateInverseY(Vector *vector, double angle) {

    return(vector->y*cos(angle) - vector->x*sin(angle));

}


double modulus (Vector *vector) {

    return(sqrt((vector->x)*(vector->x) + (vector->y)*(vector->y) + (vector->z)*(vector->z)));

}



void propagate(Vector *position, Vector angle, double distance) {

    position->x = position->x + angle.x*distance;
    position->y = position->y + angle.y*distance;
    position->z = position->z + angle.z*distance;

}

void reflect( Vector *vectorIn, Vector normal) {

    double twicevidotvn;

    twicevidotvn = 2*((vectorIn->x)*normal.x + (vectorIn->y)*normal.y + (vectorIn->z)*normal.z);

    vectorIn->x = vectorIn->x - twicevidotvn*normal.x;
    vectorIn->y = vectorIn->y - twicevidotvn*normal.y;
    vectorIn->z = vectorIn->z - twicevidotvn*normal.z;

}

void refract (Vector *vectorIn, Vector normal, double n1, double n2) {

    double vidotvn;
    double vidotvn2;
    double twicevidotvn;
    double a, b, b2, d2;

    vidotvn = (vectorIn->x)*normal.x + (vectorIn->y)*normal.y + (vectorIn->z)*normal.z;

    if (vidotvn>=0.0) {
        normal.x = -normal.x;
        normal.y = -normal.y;
        normal.z = -normal.z;
        vidotvn = -vidotvn;
    }

    vidotvn2 = vidotvn*vidotvn;
    b = n1/n2;
    b2 = b*b;
    d2 = 1.0 - b2*(1.0 - vidotvn2);

    if (d2 >= 0.0) {
        if (vidotvn >= 0.0) {
            a = -b*vidotvn + sqrt(d2);
            vectorIn->x = a*normal.x + b*(vectorIn->x);
            vectorIn->y = a*normal.y + b*(vectorIn->y);
            vectorIn->z = a*normal.z + b*(vectorIn->z);
        } else {
            a = -b*vidotvn-sqrt(d2);
            vectorIn->x = a*normal.x + b*(vectorIn->x);
            vectorIn->y = a*normal.y + b*(vectorIn->y);
            vectorIn->z = a*normal.z + b*(vectorIn->z);
        }
    } else {
        // Total Internal Reflection
        twicevidotvn = vidotvn + vidotvn;
        vectorIn->x = -twicevidotvn*normal.x + (vectorIn->x);
        vectorIn->y = -twicevidotvn*normal.y + (vectorIn->y);
        vectorIn->z = -twicevidotvn*normal.z + (vectorIn->z);
    }

}

void zernikes(double *zernikeR, double *zernikePhi, double *zernikeRGrid, double *zernikePhiGrid,
              double *zernikeNormalR, double *zernikeNormalPhi, long numelements, long nzernikes) {

    double r, phi;

    for (long i = 0; i < numelements; i++) {

        r = (static_cast<double>(i))/(static_cast<double>(numelements) - 1);
        phi = (static_cast<double>(i))/(static_cast<double>(numelements))*2*PI;

        *(zernikeRGrid + i) = r;
        *(zernikePhiGrid + i) = phi;

        // defocus 2/0
        *(zernikeR + 0*numelements + i) = 1.0;
        *(zernikePhi + 0*numelements + i) = 1.0;
        *(zernikeNormalR + 0*numelements + i) = 0.0;
        *(zernikeNormalPhi + 0*numelements + i) = 0.0;

        // defocus 2/0
        *(zernikeR + 1*numelements + i) = 2*r;
        *(zernikePhi + 1*numelements + i) = cos(phi);
        *(zernikeNormalR + 1*numelements + i) = 2.0;
        *(zernikeNormalPhi + 1*numelements + i) = -sin(phi);

        // defocus 2/0
        *(zernikeR + 2*numelements + i) = 2*r;
        *(zernikePhi + 2*numelements + i) = sin(phi);
        *(zernikeNormalR + 2*numelements + i) = 2.0;
        *(zernikeNormalPhi + 2*numelements + i) = cos(phi);

        // defocus 2/0
        *(zernikeR + 3*numelements + i) = sqrt(3.0)*(2*r*r - 1);
        *(zernikePhi + 3*numelements + i) = 1.0;
        *(zernikeNormalR + 3*numelements + i) = sqrt(3.0)*(4*r);
        *(zernikeNormalPhi + 3*numelements + i) = 0.0;

        // astigmatism 2/2
        *(zernikeR + 4*numelements + i) = sqrt(6.0)*r*r;
        *(zernikePhi + 4*numelements + i) = sin(2*phi);
        *(zernikeNormalR + 4*numelements + i) = sqrt(6.0)*2*r;
        *(zernikeNormalPhi + 4*numelements + i) = cos(2*phi)*2;

        *(zernikeR + 5*numelements + i) = sqrt(6.0)*r*r;
        *(zernikePhi + 5*numelements + i) = cos(2*phi);
        *(zernikeNormalR + 5*numelements + i) = sqrt(6.0)*2*r;
        *(zernikeNormalPhi + 5*numelements + i) = -sin(2*phi)*2;

        // coma 3/1
        *(zernikeR + 6*numelements + i) = sqrt(8.0)*(3*r*r*r - 2*r);
        *(zernikePhi + 6*numelements + i) = sin(phi);
        *(zernikeNormalR + 6*numelements + i) = sqrt(8.0)*(9*r*r - 2);
        *(zernikeNormalPhi + 6*numelements + i) = cos(phi);

        *(zernikeR + 7*numelements + i) = sqrt(8.0)*(3*r*r*r - 2*r);
        *(zernikePhi + 7*numelements + i) = cos(phi);
        *(zernikeNormalR + 7*numelements + i) = sqrt(8.0)*(9*r*r - 2);
        *(zernikeNormalPhi + 7*numelements + i) = -sin(phi);

        // trefoil 3/3
        *(zernikeR + 8*numelements + i) = sqrt(8.0)*(r*r*r);
        *(zernikePhi + 8*numelements + i) = sin(3*phi);
        *(zernikeNormalR + 8*numelements + i) = sqrt(8.0)*(3*r*r);
        *(zernikeNormalPhi + 8*numelements + i) = cos(3*phi)*3;

        *(zernikeR + 9*numelements + i) = sqrt(8.0)*(r*r*r);
        *(zernikePhi + 9*numelements + i) = cos(3*phi);
        *(zernikeNormalR + 9*numelements + i) = sqrt(8.0)*(3*r*r);
        *(zernikeNormalPhi + 9*numelements + i) = -sin(3*phi)*3;

        // spherical aberration 4/0
        *(zernikeR + 10*numelements + i) = sqrt(5.0)*(6*r*r*r*r - 6*r*r + 1);
        *(zernikePhi + 10*numelements + i) = 1.0;
        *(zernikeNormalR + 10*numelements + i) = sqrt(5.0)*(24*r*r*r - 12*r);
        *(zernikeNormalPhi + 10*numelements + i) = 0.0;

        // 2nd astigmatism 4/2
        *(zernikeR + 11*numelements + i) = sqrt(10.0)*(4*r*r*r*r - 3*r*r);
        *(zernikePhi + 11*numelements + i) = cos(2*phi);
        *(zernikeNormalR + 11*numelements + i) = sqrt(10.0)*(16*r*r*r - 6*r);
        *(zernikeNormalPhi + 11*numelements + i) = -sin(2*phi)*2;

        *(zernikeR + 12*numelements + i) = sqrt(10.0)*(4*r*r*r*r - 3*r*r);
        *(zernikePhi + 12*numelements + i) = sin(2*phi);
        *(zernikeNormalR + 12*numelements + i) = sqrt(10.0)*(16*r*r*r - 6*r);
        *(zernikeNormalPhi + 12*numelements + i) = cos(2*phi)*2;

        // quadfoil 4/4
        *(zernikeR + 13*numelements + i) = sqrt(10.0)*(r*r*r*r);
        *(zernikePhi + 13*numelements + i) = cos(4*phi);
        *(zernikeNormalR + 13*numelements + i) = sqrt(10.0)*(4*r*r*r);
        *(zernikeNormalPhi + 13*numelements + i) = -sin(4*phi)*4;

        *(zernikeR + 14*numelements + i) = sqrt(10.0)*(r*r*r*r);
        *(zernikePhi + 14*numelements + i) = sin(4*phi);
        *(zernikeNormalR + 14*numelements + i) = sqrt(10.0)*(4*r*r*r);
        *(zernikeNormalPhi + 14*numelements + i) = cos(4*phi)*4;

        // 2nd coma 5/1
        *(zernikeR + 15*numelements + i) = sqrt(12.0)*(10*r*r*r*r*r - 12*r*r*r + 3*r);
        *(zernikePhi + 15*numelements + i) = cos(phi);
        *(zernikeNormalR + 15*numelements + i) = sqrt(12.0)*(50*r*r*r*r - 36*r*r + 3);
        *(zernikeNormalPhi + 15*numelements + i) = -sin(phi);

        *(zernikeR + 16*numelements + i) = sqrt(12.0)*(10*r*r*r*r*r - 12*r*r*r + 3*r);
        *(zernikePhi + 16*numelements + i) = sin(phi);
        *(zernikeNormalR + 16*numelements + i) = sqrt(12.0)*(50*r*r*r*r - 36*r*r + 3);
        *(zernikeNormalPhi + 16*numelements + i) = cos(phi);

        // 2nd trefoil 5/3
        *(zernikeR + 17*numelements + i) = sqrt(12.0)*(5*r*r*r*r*r - 4*r*r*r);
        *(zernikePhi + 17*numelements + i) = cos(3*phi);
        *(zernikeNormalR + 17*numelements + i) = sqrt(12.0)*(25*r*r*r*r - 12*r*r);
        *(zernikeNormalPhi + 17*numelements + i) = -sin(3*phi)*3;

        *(zernikeR + 18*numelements + i) = sqrt(12.0)*(5*r*r*r*r*r - 4*r*r*r);
        *(zernikePhi + 18*numelements + i) = sin(3*phi);
        *(zernikeNormalR + 18*numelements + i) = sqrt(12.0)*(25*r*r*r*r - 12*r*r);
        *(zernikeNormalPhi + 18*numelements + i) = cos(3*phi)*3;

        // pentafoil 5/5
        *(zernikeR + 19*numelements + i) = sqrt(12.0)*(r*r*r*r*r);
        *(zernikePhi + 19*numelements + i) = cos(5*phi);
        *(zernikeNormalR + 19*numelements + i) = sqrt(12.0)*(5*r*r*r*r);
        *(zernikeNormalPhi + 19*numelements + i) = -sin(5*phi)*5;

        *(zernikeR + 20*numelements + i) = sqrt(12.0)*(r*r*r*r*r);
        *(zernikePhi + 20*numelements + i) = sin(5*phi);
        *(zernikeNormalR + 20*numelements + i) = sqrt(12.0)*(5*r*r*r*r);
        *(zernikeNormalPhi + 20*numelements + i) = cos(5*phi)*5;

        *(zernikeR + 21*numelements + i) = sqrt(7.0)*(20.0*r*r*r*r*r*r - 30.0*r*r*r*r + 12.0*r*r - 1.0);
        *(zernikePhi + 21*numelements + i) = 1.0;
        *(zernikeNormalR + 21*numelements + i) = sqrt(7.0)*(120.0*r*r*r*r*r - 120.0*r*r*r + 24.0*r);
        *(zernikeNormalPhi + 21*numelements + i) = 0.0;

        *(zernikeR + 22*numelements + i) = sqrt(14.0)*(15.0*r*r*r*r*r*r - 20.0*r*r*r*r + 6*r*r);
        *(zernikePhi + 22*numelements + i) = sin(2*phi);
        *(zernikeNormalR + 22*numelements + i) = sqrt(14.0)*(90.0*r*r*r*r*r - 80.0*r*r*r + 12.0*r);
        *(zernikeNormalPhi + 22*numelements + i) = cos(2*phi)*2;

        *(zernikeR + 23*numelements + i) = sqrt(14.0)*(15.0*r*r*r*r*r*r - 20.0*r*r*r*r + 6*r*r);
        *(zernikePhi + 23*numelements + i) = cos(2*phi);
        *(zernikeNormalR + 23*numelements + i) = sqrt(14.0)*(90.0*r*r*r*r*r - 80.0*r*r*r + 12.0*r);
        *(zernikeNormalPhi + 23*numelements + i) = -sin(2*phi)*2;

        *(zernikeR + 24*numelements + i) = sqrt(14.0)*(6.0*r*r*r*r*r*r - 5.0*r*r*r*r);
        *(zernikePhi + 24*numelements + i) = sin(4*phi);
        *(zernikeNormalR + 24*numelements + i) = sqrt(14.0)*(36.0*r*r*r*r*r - 20.0*r*r*r);
        *(zernikeNormalPhi + 24*numelements + i) = cos(4*phi)*4;

        *(zernikeR + 25*numelements + i) = sqrt(14.0)*(6.0*r*r*r*r*r*r - 5.0*r*r*r*r);
        *(zernikePhi + 25*numelements + i) = cos(4*phi);
        *(zernikeNormalR + 25*numelements + i) = sqrt(14.0)*(36.0*r*r*r*r*r - 20.0*r*r*r);
        *(zernikeNormalPhi + 25*numelements + i) = -sin(4*phi)*4;

        *(zernikeR + 26*numelements + i) = sqrt(14.0)*(r*r*r*r*r*r);
        *(zernikePhi + 26*numelements + i) = sin(6*phi);
        *(zernikeNormalR + 26*numelements + i) = sqrt(14.0)*(6.0*r*r*r*r*r);
        *(zernikeNormalPhi + 26*numelements + i) = cos(6*phi)*6;

        *(zernikeR + 27*numelements + i) = sqrt(14.0)*(r*r*r*r*r*r);
        *(zernikePhi + 27*numelements + i) = cos(6*phi);
        *(zernikeNormalR + 27*numelements + i) = sqrt(14.0)*(6.0*r*r*r*r*r);
        *(zernikeNormalPhi + 27*numelements + i) = -sin(6*phi)*6;


    }


}

//Chebyshev polynomial of the first kind
double chebyshevT (int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    return 2*x*chebyshevT(n - 1, x) - chebyshevT(n - 2, x);
}

//differentiation of Chebyshev polynomial of the first kind
double chebyshevT_x (int n, double x) {
    if (n == 0) return 0.0;
    if (n == 1) return 1.0;
    return 2*chebyshevT(n - 1, x) + 2*chebyshevT_x(n - 1, x) - chebyshevT_x(n - 2, x);
}

void chebyshevT2D (int n, double x, double y, double *t) {
    int degree, d(n + 1);
    for (degree = 0; degree <= n; degree++) {
        d -= (degree + 1);
        if (d <= 0) {
            d += degree;
            break;
        }
    }
    double tx = chebyshevT(degree - d, x);
    double ty = chebyshevT(d, y);
    double tx_x = chebyshevT_x(degree - d, x);
    double ty_y = chebyshevT_x(d, y);
    t[0] = tx*ty;
    t[1] = tx_x*ty;
    t[2] = tx*ty_y;
}

void chebyshevs (double *r_grid, double *phi_grid, double *chebyshev, double *chebyshev_r, double *chebyshev_phi, long nPoint, long nTerm) {
    for (long j = 0; j < nPoint; j++) {
        double r = r_grid[j];
        for (long l = 0; l < nPoint; l++) {
            double x = r*cos(phi_grid[l]);
            double y = r*sin(phi_grid[l]);
            for (long i = 0; i < nTerm; i++) {
                double t[3], t_r, t_phi;
                chebyshevT2D(i, x, y, &t[0]);
                t_r = t[1]*x/r + t[2]*y/r;         // dz/dr = dz/dx cos(theta) + dz/dy sin(theta)
                t_phi = r*(-t[1]*y/r + t[2]*x/r);  // dz/dtheta = r [ -dz/dx sin(theta) + dz/dy cos(theta) ]
                chebyshev[i*nPoint*nPoint + l*nPoint + j] = t[0];
                chebyshev_r[i*nPoint*nPoint + l*nPoint + j] = t_r;
                chebyshev_phi[i*nPoint*nPoint + l*nPoint + j] = t_phi;
            }
        }
    }
}
