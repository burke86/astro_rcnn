///
/// @package phosim
/// @file source.h
/// @brief source structure
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

enum SourceTypes {POINT=0, IMAGE=1, GAUSSIAN=2, MOVINGPOINT=4, SERSIC=5, SERSIC2D=6, PINHOLE=7, OPD=8, SERSICCOMPLEX=9, SERSICDISK=10, SERSICDISKCOMPLEX=11, DISTORTEDSPHERE=12};

struct Source {
    std::vector<double> ra;
    std::vector<double> redshift;
    std::vector<double> gamma1;
    std::vector<double> gamma2;
    std::vector<double> kappa;
    std::vector<double> deltara;
    std::vector<double> deltadec;
    std::vector<double> dec;
    std::vector<int> split;
    double *vx;
    double *vy;
    double *vz;
    double *norm;
    double *mag;
    double **spatialpar;
    double **dustpar;
    double **dustparz;
    std::vector<std::string> id;
    int *spatialtype;
    int *dusttype;
    int *dusttypez;
    int *type;
    std::vector<std::string> sedfilename;
    std::vector<std::string> spatialname;
    std::vector<std::string> dustname;
    std::vector<std::string> dustnamez;
    long *skysameas;
    long *sedptr;
};
