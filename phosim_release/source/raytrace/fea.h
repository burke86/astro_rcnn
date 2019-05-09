///
/// @package phosim
/// @file fea.h
/// @brief class to load and interpolate FEA data
///
/// @brief Created by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///


#ifndef FEA_H
#define FEA_H

#include <math.h>
#include <cstdlib>
#include <iostream>
#include "nanoflann.hpp"
//#include <Eigen/Dense>
#include "perturbation.h"
#include "surface.h"
#include "ancillary/readtext.h"

using namespace nanoflann;
//using namespace Eigen;

namespace fea {

struct Node {
    struct Point {
        double x, y;
    };
    std::vector<Point>  pts;
    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
    inline double kdtree_distance(const double *p1, const size_t idx_p2, size_t size) const {
        double d0 = p1[0] - pts[idx_p2].x;
        double d1 = p1[1] - pts[idx_p2].y;
        return d0*d0 + d1*d1;
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline double kdtree_get_pt(const size_t idx, int dim) const {
        if (dim == 0) return pts[idx].x;
        else return pts[idx].y;
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const { return false; }

};

class Fea {
    public:
        Fea();
        Fea(const std::string & feaFile, int leafSize) {
            m_leaf = leafSize;
            loadData(feaFile);
        }
        Fea(const std::string & feaFile, int leafSize, Surface & surface, Perturbation & perturbation, long surfaceIndex, double scaling, long *sizeFea) {
            m_leaf = leafSize;
            readText test(feaFile);
            std::istringstream iss(test[0]);
            size_t nx, ny;
            iss >> nx >> ny;
            size_t N = test.getSize() - 1;
            if (nx*ny == N) {  //zemax grid data
                rmax = perturbation.rmax[surfaceIndex];
                loadZmxData(feaFile, scaling);
            } else {
                loadData(feaFile);
                subtractSurface(surface, perturbation, surfaceIndex, scaling);
            }
            *sizeFea = N;
        }
        void loadData(const std::string & feaFile);
        void loadZmxData(const std::string & feaFile, double scaling);
        void radiusQuery(std::vector<double> & x, std::vector<double> & y, std::vector<double> & z, double radius, double p);
        void knnQuery(std::vector<double> & x, std::vector<double> & y, std::vector<double> & z, size_t nNear, double p);
        //void knnQueryFit(std::vector<double> & x, std::vector<double> & y, double *z, double *z_r, double *z_phi, size_t nNear, int degree);
        void knnQueryFitDegree1(std::vector<double> & x, std::vector<double> & y, double *z, double *z_r, double *z_phi, size_t nNear);

        ~Fea() {
            delete m_kdTree;
        }
        double getMinR() {
            return minR;
        }
        double getMaxR() {
            return maxR;
        }
        void getTransformation(Perturbation & perturbation, long surfaceIndex);
        void subtractSurface(Surface & surface, Perturbation & perturbation, long surfaceIndex, double scaling);
        double zernike(int n, int m, double r, double phi);
        double zernike_r(int n, int m, double r, double phi);
        double zernike_phi(int n, int m, double r, double phi);

    private:
        Node m_node;
        std::vector<double> m_z;
        int m_leaf;
        double minR, maxR, rmax;
        double meanTx, meanTy, meanTz, meanRx, meanRy, meanRz;
        double meanZ;
        double rotationMatrix[3][3];

        typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, Node >, Node, 2 /* dim */> KDTreeAdaptor;
        KDTreeAdaptor* m_kdTree;
};

} // namespace fea

#endif
