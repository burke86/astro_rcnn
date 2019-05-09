///
/// @package phosim
/// @file fea.cpp
/// @brief class to load and interpolate FEA data
///
/// @brief Created by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "fea.h"

namespace fea {

    void Fea::loadData(const std::string & feaFile) {
    readText feaPars(feaFile);
    size_t N = feaPars.getSize();
    m_node.pts.resize(N);
    m_z.resize(N);
    minR = 1e10;
    maxR = 0.0;
    meanTx = 0.0;
    meanTy = 0.0;
    meanTz = 0.0;
    meanRx = 0.0;
    meanRy = 0.0;
    meanRz = 0.0;
    double sx, cx, sy, cy, sz, cz;
    double a11, a12, a13, a21, a22, a23, a31, a32, a33;
    for (size_t t(0); t < N; t++){

        double x0, y0, z0;
        double tx, ty, tz, rx, ry, rz;

        std::istringstream iss(feaPars[t]);
        iss >> x0 >> y0 >> z0 >> tx >> ty >> tz >> rx >> ry >> rz;
        rx *= (M_PI/180.0);
        ry *= (M_PI/180.0);
        rz *= (M_PI/180.0);
        meanTx += tx;
        meanTy += ty;
        meanTz += tz;
        meanRx += rx;
        meanRy += ry;
        meanRz += rz;
        sx = sin(rx);
        cx = cos(rx);
        sy = sin(ry);
        cy = cos(ry);
        sz = sin(rz);
        cz = cos(rz);
        a11 = cz*cy;
        a12 = cz*sx*sy + sz*cx;
        a13 = -cz*cx*sy + sz*sx;
        a21 = -sz*cy;
        a22 = -sz*sx*sy + cz*cx;
        a23 = sz*cx*sy + cz*sx;
        a31 = sy;
        a32 = -sx*cy;
        a33 = cx*cy;
        //Assuming tilts are done first, in the order x, y, and z, and then the decenters are done.
        m_node.pts[t].x = a11*x0 + a12*y0 + a13*z0;
        m_node.pts[t].y = a21*x0 + a22*y0 + a23*z0;
        m_z[t]          = a31*x0 + a32*y0 + a33*z0;
        m_node.pts[t].x +=  tx;
        m_node.pts[t].y += ty;
        m_z[t] += tz;
    }
    meanTx /= N;
    meanTy /= N;
    meanTz /= N;
    meanRx /= N;
    meanRy /= N;
    meanRz /= N;
    sx = sin(meanRx);
    cx = cos(meanRx);
    sy = sin(meanRy);
    cy = cos(meanRy);
    sz = sin(meanRz);
    cz = cos(meanRz);
    a11 = cz*cy;
    a21 = cz*sx*sy + sz*cx;
    a31 = -cz*cx*sy + sz*sx;
    a12 = -sz*cy;
    a22 = -sz*sx*sy + cz*cx;
    a32 = sz*cx*sy + cz*sx;
    a13 = sy;
    a23 = -sx*cy;
    a33 = cx*cy;
    rotationMatrix[0][0] = a11;
    rotationMatrix[1][0] = a21;
    rotationMatrix[2][0] = a31;
    rotationMatrix[0][1] = a12;
    rotationMatrix[1][1] = a22;
    rotationMatrix[2][1] = a32;
    rotationMatrix[0][2] = a13;
    rotationMatrix[1][2] = a23;
    rotationMatrix[2][2] = a33;

    for (size_t t(0); t < N; t++){
        double x, y, z, xp, yp, zp;
        x = m_node.pts[t].x - meanTx;
        y = m_node.pts[t].y - meanTy;
        z = m_z[t] - meanTz;
        xp = a11*x + a12*y + a13*z;
        yp = a21*x + a22*y + a23*z;
        zp = a31*x + a32*y + a33*z;
        m_node.pts[t].x = xp;
        m_node.pts[t].y = yp;
        m_z[t]          = zp;
        double r = xp*xp + yp*yp;
        if (r>maxR) maxR = r;
        if (r<minR) minR = r;
    }
    maxR = sqrt(maxR);
    minR = sqrt(minR);

    m_kdTree = new KDTreeAdaptor(2, m_node, KDTreeSingleIndexAdaptorParams(m_leaf));
    m_kdTree->buildIndex();
}

    void Fea::loadZmxData(const std::string & feaFile, double scaling) {
    readText feaPars(feaFile);
    size_t N = feaPars.getSize() - 1;
    m_node.pts.resize(N*9);
    m_z.resize(N*9);
    std::istringstream iss(feaPars[0]);
    size_t nx, ny;
    double dx, dy;
    iss >> nx >> ny >> dx >> dy;
    long nn = 0;
    for (size_t t(0); t < N; t++) {
        std::istringstream iss(feaPars[t + 1]);
        double z1, z2, z3, z4;
        iss >> z1 >> z2 >> z3 >> z4;
        z1 *= scaling;
        z2 *= scaling;
        z3 *= scaling;
        z4 *= scaling;
        // m_z[4*t] = -z1;
        size_t i = t % nx;
        size_t j = t / nx;
        double x0, y0, x, y;
        x0 = (i - nx/2.0 + 0.5)*dx;
        y0 = -(j - ny/2.0 + 0.5)*dy;

        for (long ii=-1; ii<=1; ii++) {
            for (long jj = -1; jj <= 1; jj++) {
                x = x0 + ii*dx/3.0;
                y = y0 + jj*dy/3.0;
                m_z[nn] = (z1 + ii*z2*dx/3.0 + jj*z3*dy/3.0 + ii*jj*z4*dx/3.0*dy/3.0);
                m_node.pts[nn].x = x;
                m_node.pts[nn].y = y;
                if (m_z[nn] < 10e-3 && m_z[nn] > -10e-3) {
                } else {
                    printf("%e %ld\n", m_z[nn], nn);
                }
                nn++;
            }
        }
    }
    m_node.pts.resize(nn);
    m_z.resize(nn);
    meanTx = 0.0;
    meanTy = 0.0;
    meanTz = 0.0;
    rotationMatrix[0][0] = 1.0;
    rotationMatrix[1][0] = 0.0;
    rotationMatrix[2][0] = 0.0;
    rotationMatrix[0][1] = 0.0;
    rotationMatrix[1][1] = 1.0;
    rotationMatrix[2][1] = 0.0;
    rotationMatrix[0][2] = 0.0;
    rotationMatrix[1][2] = 0.0;
    rotationMatrix[2][2] = 1.0;

    m_kdTree = new KDTreeAdaptor(2, m_node, KDTreeSingleIndexAdaptorParams(m_leaf));
    m_kdTree->buildIndex();
}


void Fea::radiusQuery(std::vector<double> & x, std::vector<double> & y, std::vector<double> & z, double radius, double p) {
    double radiusSqr = radius*radius;
    for (size_t t(0); t < x.size(); t++) {
        std::vector<std::pair<size_t, double> > indices_dists;
        nanoflann::RadiusResultSet<double, size_t> resultSet(radiusSqr, indices_dists);
        double query_pt[2] = {x[t], y[t]};
        m_kdTree->findNeighbors(resultSet, query_pt, nanoflann::SearchParams(m_leaf));
        double wsum = 0.0, wz = 0.0;
        for  (size_t n(0); n < indices_dists.size(); n++) {
            if (indices_dists[n].second < 1e-20) {
                wz = m_z[indices_dists[n].first];
                wsum = 1.0;
                break;
            }
            double dist = sqrt(indices_dists[n].second);
            double w = pow((radius - dist)/radius/dist, p);
            wz += m_z[indices_dists[n].first]*w;
            wsum += w;
        }
        z[t] = wz/wsum;
    }
}

void Fea::knnQuery(std::vector<double> & x, std::vector<double> & y, std::vector<double> & z, size_t nNear, double p) {
    for (size_t t(0); t < x.size(); t++) {
        std::vector<size_t> ret_indexes(nNear);
        std::vector<double> out_dists_sqr(nNear);
        nanoflann::KNNResultSet<double> resultSet(nNear);
        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        double query_pt[2] = {x[t], y[t]};
        m_kdTree->findNeighbors(resultSet, query_pt, nanoflann::SearchParams(m_leaf));
        double wsum = 0.0, wz = 0.0;
        for  (size_t n(0); n < nNear; n++) {
            if (out_dists_sqr[n] < 1e-20) {
                wz = m_z[ret_indexes[n]];
                wsum = 1.0;
                break;
            }
            double w = pow(out_dists_sqr[n], -p/2.0);
            wz += m_z[ret_indexes[n]]*w;
            wsum += w;
        }
        z[t] = wz/wsum;
    }
}
/*
void Fea::knnQueryFit(std::vector<double> & x, std::vector<double> & y, double *z, double *z_r, double *z_phi, size_t nNear, int degree) {
    int coordType=2; //1: Cartesian  2: polar
    int sameNeighbors=0;
    double minDistanceSqr=0.;
    size_t t0(0);
    int nTerm(0);
    double r, phi;
    for (int i=0;i<=degree;i++) nTerm += (i+1);
    VectorXd c0(nTerm), c(nTerm);
    for (size_t t(0);t<x.size();t++) {
        //if (t>0) {
        //    double rSqr=pow(x[t]-x[t0],2)+pow(y[t]-y[t0],2);
        //    sameNeighbors=(rSqr<minDistanceSqr)? 1:0;
        //}
        if (sameNeighbors == 0) {
            std::vector<size_t> ret_indexes(nNear);
            std::vector<double> out_dists_sqr(nNear);
            nanoflann::KNNResultSet<double> resultSet(nNear);
            resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
            double query_pt[2]={x[t],y[t]};
            m_kdTree->findNeighbors(resultSet, query_pt, nanoflann::SearchParams(m_leaf));
            minDistanceSqr=out_dists_sqr[0];
            MatrixXd A(nNear,nTerm);
            VectorXd b(nNear);
            for  (size_t n(0);n<nNear;n++) {
                size_t idx=ret_indexes[n];
                A(n,0)=1.;
                int mm(1);
                if (coordType == 2) {
                    r=sqrt(m_node.pts[idx].x*m_node.pts[idx].x+m_node.pts[idx].y*m_node.pts[idx].y);
                    phi=atan2(m_node.pts[idx].y,m_node.pts[idx].x);
                }
                for (int d(1);d<=degree;d++) {
                    for (int dp(0);dp<=d;dp++) {
                        if (coordType == 1) A(n,mm)=pow(m_node.pts[idx].x,dp)*pow(m_node.pts[idx].y,d-dp);
                        else if (coordType == 2) A(n,mm)=zernike(d,dp,r,phi);
                        mm++;
                    }
                }
                b(n)=m_z[idx];
            }
            //c = A.fullPivLu().solve(b);  //solve A*c=b
            c = A.colPivHouseholderQr().solve(b);  //solve A*c=b
            c0=c;
            t0=t;

        } else {
            c=c0;
        }
        z[t]=c(0);
        z_r[t]=0.;
        z_phi[t]=0.;
        int m(1);
        double z_x=0., z_y=0.;
        r=sqrt(x[t]*x[t]+y[t]*y[t]);
        phi=atan2(y[t],x[t]);
        for (int d(1);d<=degree;d++) {
            for (int dp(0);dp<=d;dp++) {
                if (coordType == 1) {
                    z[t] += c(m)*pow(x[t],dp)*pow(y[t],d-dp);
                    if(dp>=1) z_x += c(m)*dp*pow(x[t],dp-1)*pow(y[t],d-dp);
                    if(d-dp>=1) z_y += c(m)*(d-dp)*pow(x[t],dp)*pow(y[t],d-dp-1);
                } else if (coordType == 2) {
                    z[t] += c(m)*zernike(d,dp,r,phi);
                    z_r[t] += c(m)*zernike_r(d,dp,r,phi)*rmax;
                    z_phi[t] += c(m)*zernike_phi(d,dp,r,phi);
                }
                m++;
            }
        }
        if (coordType == 1) {
            if (r>0) {
                z_r[t]=(z_x*x[t]/r+z_y*y[t]/r)*rmax;     // dz/dr = dz/dx cos(theta) + dz/dy sin(theta)
                z_phi[t]=r*(-z_x*y[t]/r+z_y*x[t]/r);     // dz/dtheta = r [ -dz/dx sin(theta) + dz/dy cos(theta) ]
            }
        }
    }
}*/

void Fea::knnQueryFitDegree1(std::vector<double> & x, std::vector<double> & y, double *z, double *z_r, double *z_phi, size_t nNear) {
    int coordType = 2; //1: Cartesian  2: polar
    int sameNeighbors = 0;
    double minDistanceSqr = 0.0;
    size_t t0(0);
    int degree(1);
    double r, phi;
    double c0[3], c[3], H[3][3];
    std::vector<double> b(nNear);
    for (size_t t(0); t < x.size(); t++) {
        //if (t>0) {
        //    double rSqr = pow(x[t]-x[t0],2) + pow(y[t]-y[t0],2);
        //    sameNeighbors = (rSqr<minDistanceSqr)? 1:0;
        //}
        if (sameNeighbors == 0) {
            std::vector<size_t> ret_indexes(nNear);
            std::vector<double> out_dists_sqr(nNear);
            nanoflann::KNNResultSet<double> resultSet(nNear);
            resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
            double query_pt[2] = {x[t], y[t]};
            m_kdTree->findNeighbors(resultSet, query_pt, nanoflann::SearchParams(m_leaf));
            minDistanceSqr = out_dists_sqr[0];
            /* c =  inverse(transpose(A) A) transpose(A) b = H transpose(A) b */
            /* ct = transpose(A) b */
            double s0(0.0), sx(0.0), sy(0.0), vx(0.0), vy(0.0), xy(0.0), a1, a2;
            double ct[3] = {0.0, 0.0, 0.0};
            for  (size_t n(0); n < nNear; n++) {
                size_t idx = ret_indexes[n];
                if (coordType == 1) {
                    a1 = m_node.pts[idx].y;
                    a2 = m_node.pts[idx].x;
                } else if (coordType == 2) {
                    r = sqrt(m_node.pts[idx].x*m_node.pts[idx].x + m_node.pts[idx].y*m_node.pts[idx].y);
                    phi = atan2(m_node.pts[idx].y, m_node.pts[idx].x);
                    a1 = zernike(1, 0, r, phi);
                    a2 = zernike(1, 1, r, phi);
                }
                b[n] = m_z[idx];
                s0 += 1;
                sx += a1;
                sy += a2;
                vx += a1*a1;
                vy += a2*a2;
                xy += a1*a2;
                ct[0] += b[n];
                ct[1] += a1*b[n];
                ct[2] += a2*b[n];
            }
            double h0 = (s0*xy*xy - 2*xy*sx*sy + vy*sx*sx + vx*sy*sy - s0*vy*vx);
            H[0][0] = (xy*xy - vx*vy)/h0;
            H[1][0] = (-xy*sy + sx*vy)/h0;
            H[2][0] = (-xy*sx + sy*vx)/h0;
            H[0][1] = (-xy*sy + sx*vy)/h0;
            H[1][1] = (sy*sy - s0*vy)/h0;
            H[2][1] = (s0*xy - sy*sx)/h0;
            H[0][2] = (-xy*sx + sy*vx)/h0;
            H[1][2] = (s0*xy - sy*sx)/h0;
            H[2][2] = (sx*sx - s0*vx)/h0;
            for (size_t i(0); i < 3; i++) {
                c[i] = 0.0;
                for (size_t j(0); j < 3; j++)
                    c[i] += H[i][j]*ct[j];
            }
            for (size_t i(0); i < 3; i++) {
                c0[i] = c[i];
            }
            t0 = t;
        } else {
            for (size_t i(0); i < 3; i++) {
                c[i] = c0[i];
            }
        }
        z[t] += c[0];
        z_r[t] += 0.0;
        z_phi[t] += 0.0;
        int m(1);
        double z_x = 0.0;
        double z_y = 0.0;
        r = sqrt(x[t]*x[t] + y[t]*y[t]);
        phi = atan2(y[t], x[t]);
        for (int d(1); d <= degree; d++) {
            for (int dp(0); dp <= d; dp++) {
                if (coordType == 1) {
                    z[t] += c[m]*pow(x[t], dp)*pow(y[t], d - dp);
                    if (dp >= 1) z_x += c[m]*dp*pow(x[t], dp - 1)*pow(y[t], d - dp);
                    if (d - dp >= 1) z_y += c[m]*(d - dp)*pow(x[t], dp)*pow(y[t], d - dp - 1);
                } else if (coordType == 2) {
                    z[t] += c[m]*zernike(d, dp, r, phi);
                    z_r[t] += c[m]*zernike_r(d, dp, r, phi)*rmax;
                    z_phi[t] += c[m]*zernike_phi(d, dp, r, phi);
                }
                m++;
            }
        }
        if (coordType == 1) {
            if (r>0) {
                z_r[t] += (z_x*x[t]/r + z_y*y[t]/r)*rmax;     // dz/dr = dz/dx cos(theta) + dz/dy sin(theta)
                z_phi[t] += r*(-z_x*y[t]/r + z_y*x[t]/r);     // dz/dtheta = r [ -dz/dx sin(theta) + dz/dy cos(theta) ]
            }
        }
    }
}

    void Fea::subtractSurface(Surface & surface, Perturbation & perturbation, long surfaceIndex, double scaling) {
    rmax = perturbation.rmax[surfaceIndex];
    double radiusofcurv = -surface.radiusCurvature[surfaceIndex];
    double third = surface.three[surfaceIndex]*1e3;
    double fourth = surface.four[surfaceIndex]*1e3;
    double fifth = surface.five[surfaceIndex]*1e3;
    double sixth = surface.six[surfaceIndex]*1e3;
    double seventh = surface.seven[surfaceIndex]*1e3;
    double eighth = surface.eight[surfaceIndex]*1e3;
    double ninth = surface.nine[surfaceIndex]*1e3;
    double tenth = surface.ten[surfaceIndex]*1e3;
    double conic = surface.conic[surfaceIndex];

    if (radiusofcurv != 0) {
        for (size_t t(0); t < m_z.size(); t++) {
            double r = sqrt(m_node.pts[t].x*m_node.pts[t].x + m_node.pts[t].y*m_node.pts[t].y);
            double h = surface.height[surfaceIndex];
            h -= (pow(r, 2.0)/radiusofcurv/(1.0 + sqrt(1.0 - (conic + 1.0)*pow(r/radiusofcurv, 2.0))) +
                  third*pow(r, 3.0) + fourth*pow(r, 4.0) +
                  fifth*pow(r, 5.0) + sixth*pow(r, 6.0) +
                  seventh*pow(r, 7.0) + eighth*pow(r, 8.0) +
                  ninth*pow(r, 9.0) + tenth*pow(r, 10.0));
            m_z[t] -= h;
        }
    } else {
        for (size_t t(0); t < m_z.size(); t++) {
            m_z[t] -= surface.height[surfaceIndex];
        }
    }

    for (size_t t(0); t < m_z.size(); t++) {
        m_z[t] *= scaling;
    }

}
void Fea::getTransformation(Perturbation & perturbation, long surfaceIndex) {
    //phosim convention
    //x_optics = Rotation (x_lab - dx)
    //x_lab = InverseRotation (x_optics) + dx
    perturbation.defocus[surfaceIndex] = meanTz;
    perturbation.decenterX[surfaceIndex] = meanTx;
    perturbation.decenterY[surfaceIndex] = meanTy;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            perturbation.rotationmatrix[surfaceIndex*9 + i*3 + j] = rotationMatrix[i][j];
            perturbation.inverserotationmatrix[surfaceIndex*9 + i*3 + j] = rotationMatrix[j][i];
        }
    }
}

double Fea::zernike(int n, int m, double r, double phi) {

    double s(0);
    if (n == 1) {
        if (m == 0) s = 2*r*sin(phi);
        else s = 2*r*cos(phi);
    } else if (n == 2) {
        if (m == 0) s = sqrt(3.0)*(2*r*r - 1);
        else if (m == 1) s = sqrt(6.0)*r*r*cos(2*phi);
        else s = sqrt(6.0)*r*r*sin(2*phi);
    } else if (n == 3) {
        if (m == 0) s = sqrt(8.0)*(3*r*r*r - 2*r)*cos(phi);
        else if (m == 1) s = sqrt(8.0)*(3*r*r*r - 2*r)*sin(phi);
        else if (m == 2) s = sqrt(8.0)*(r*r*r)*cos(3*phi);
        else s = sqrt(8.0)*(r*r*r)*sin(3*phi);
    } else if (n == 4) {
        if (m == 0) s = sqrt(5.0)*(6*r*r*r*r - 6*r*r + 1);
        else if (m == 1) s = sqrt(10.0)*(4*r*r*r*r - 3*r*r)*cos(2*phi);
        else if (m == 2) s = sqrt(10.0)*(4*r*r*r*r - 3*r*r)*sin(2*phi);
        else if (m == 3) s = sqrt(10.0)*(r*r*r*r)*cos(4*phi);
        else s = sqrt(10.0)*(r*r*r*r)*sin(4*phi);
    }
    return s;
}

double Fea::zernike_r(int n, int m, double r, double phi) {
    double s(0);
    if (n == 1) {
        if (m == 0) s = 2*sin(phi);
        else s = 2*cos(phi);
    } else if (n == 2) {
        if (m == 0) s = sqrt(3.0)*(4*r);
        else if (m == 1) s = sqrt(6.0)*2*r*cos(2*phi);
        else s = sqrt(6.0)*2*r*sin(2*phi);
    } else if (n == 3) {
        if (m == 0) s = sqrt(8.0)*(9*r*r - 2)*cos(phi);
        else if (m == 1) s = sqrt(8.0)*(9*r*r - 2)*sin(phi);
        else if (m == 2) s = sqrt(8.0)*(3*r*r)*cos(3*phi);
        else s = sqrt(8.0)*(3*r*r)*sin(3*phi);
    } else if (n == 4) {
        if (m == 0) s = sqrt(5.0)*(24*r*r*r - 12*r);
        else if (m == 1) s = sqrt(10.0)*(16*r*r*r - 6*r)*cos(2*phi);
        else if (m == 2) s = sqrt(10.0)*(16*r*r*r - 6*r)*sin(2*phi);
        else if (m == 3) s = sqrt(10.0)*(4*r*r*r)*cos(4*phi);
        else s = sqrt(10.0)*(4*r*r*r)*sin(4*phi);
    }
    return s;
}

double Fea::zernike_phi(int n, int m, double r, double phi) {

    double s(0);
    if (n == 1) {
        if (m == 0) s = 2*r*cos(phi);
        else s = -2*r*sin(phi);
    } else if (n == 2) {
        if (m == 0) s = 0.0;
        else if (m == 1) s = -2*sqrt(6.0)*r*r*sin(2*phi);
        else s = 2*sqrt(6.0)*r*r*cos(2*phi);
    } else if (n == 3) {
        if (m == 0) s = -sqrt(8.0)*(3*r*r*r - 2*r)*sin(phi);
        else if (m == 1) s = sqrt(8.0)*(3*r*r*r - 2*r)*cos(phi);
        else if (m == 2) s = -3*sqrt(8.0)*(r*r*r)*sin(3*phi);
        else s = 3*sqrt(8.0)*(r*r*r)*cos(3*phi);
    } else if (n == 4) {
        if (m == 0) s = 0.0;
        else if (m == 1) s = -2*sqrt(10.0)*(4*r*r*r*r - 3*r*r)*sin(2*phi);
        else if (m == 2) s = 2*sqrt(10.0)*(4*r*r*r*r - 3*r*r)*cos(2*phi);
        else if (m == 3) s = -4*sqrt(10.0)*(r*r*r*r)*sin(4*phi);
        else s = 4*sqrt(10.0)*(r*r*r*r)*cos(4*phi);
    }
    return s;
}

}
