#include "grating.h"

Grating::Grating() {

    _angleBlazeRad = kBlazeAngleDeg * DEGREE;
    _integral.resize(kMaxNumIntervals);
    _numSlits = kNumSlits;
    _deltaNm = kDeltaNm;

}

Grating::Grating(double blazeAngleDeg, int numSlits, double deltaNm)  {

    _angleBlazeRad =  blazeAngleDeg * DEGREE;
    _numSlits = numSlits;
    _deltaNm = deltaNm;
    _integral.resize(kMaxNumIntervals);

}

Grating::~Grating() { }

double Grating::_calculateFunction(double angleInRad, double angleOutRad,
                  double wavelengthNm) {

//    calculate the value of the intensity function
    double vPrime;
    double v;
    vPrime = PI*_deltaNm * (sin(angleOutRad) + sin(angleInRad)) / wavelengthNm;
    v = PI*_deltaNm*cos(_angleBlazeRad) * (sin(angleOutRad - _angleBlazeRad) + sin(angleInRad - _angleBlazeRad) ) / wavelengthNm;
    return pow(sin(_numSlits * vPrime) / (_numSlits * sin(vPrime)), 2.0) * pow(sin(v)/v, 2.0);

}

void Grating::_makeTable(double angleInRad, double wavelengthNm) {
//    make a table of integration based on the angleOutRad.
//    so far I don't know if there's a more efficient way
//    to do the integration. For the numerical method of
//    integration, I think the bottleneck is actually from
//    some of the operations like sin or square. Some other
//    numerical method may have less kStepRad for each wavelengthNm
//    but may result in more inefficient operations. Anyway,
//    I'll try other methods later.

    double angleOutRad;
    int i;

    _integral.clear();
    _integral.resize(kMaxNumIntervals);

    angleOutRad = -(PI/2.0);
    for (i = 0; i < kMaxNumIntervals; i++){
        angleOutRad += kStepRad;
        if (i == 0){
            _integral.at(i) = kStepRad *
                _calculateFunction(angleInRad, angleOutRad, wavelengthNm);
        } else {
            _integral.at(i) = _integral.at(i - 1) + kStepRad *
                _calculateFunction(angleInRad, angleOutRad, wavelengthNm);
        }
    }
    return;
}

int Grating::_binarySearch(double goal) {
//  return the index of the interval which correspond to the upper_limit of
//  the integration.

    int l = 0;
    int r = kMaxNumIntervals - 1;
    int mid;

    while (l <= r){
        mid = l + ((r - l) >> 1);
        if (goal > _integral.at(mid) ){
            l = mid + 1;
        } else {
            r = mid - 1;
        }
    }
    return mid;
}



double Grating::_calculateAngle(double angleInRad, double wavelengthNm) {
// calculate the angleOutRad of each photon after grating

    double angleOutRad;
    double upperLimit;

    _makeTable(angleInRad, wavelengthNm);
    upperLimit = rand()/(RAND_MAX + 1.0) * _integral.at(kMaxNumIntervals - 1);
    angleOutRad = -(PI/2.0) + kStepRad * _binarySearch(upperLimit);
    return angleOutRad;
}


void Grating::diffract(double vxIn, double vyIn, double vzIn,
                       double vxGratingNormal, double vyGratingNormal,
                       double vzGratingNormal, double& vxOut,
                       double& vyOut, double& vzOut, double wavelengthNm) {
    // Main method: Calculates the direction out of a photon after it interacts
    // with the grating.
    // vxIn,vyIn,vzIn is input direction of the photon
    // vxOut,vyOut,vzOut is output direction of the photon.

    _setGratingNormal(vxGratingNormal, vyGratingNormal, vzGratingNormal);

    // To use the grating, we need at first transform the phothons from the lab
    // frame to the optic frame, and then call THIS function. After that,
    // we need to transform the coordinates back to the lab frame.
    // (call transform_inverse)
    //     o_hat = a * i_hat + b * n_hat;
    //     vnDotvi = n_hat * i_hat = cos(angleInRad)
    //     I'm not quite sure about these geometries because so far
    //   I didn't test it in the whole LSST simulation. Please help
    //   me double check that.
    // GHS: We need negative of vin vector to get angleInRad
    // Put a test in here later to see if photon is hitting bottom of grating

    double vnDotvi = -(vxIn*_vxGratingNormal + vyIn*_vyGratingNormal +
                       vzIn*_vzGratingNormal); //This should be just -vzIn if N=0,0,1
    double angleInRad = acos(vnDotvi);

    // Change sign on angleout for our calculation. Put in test later.
    //If it comes back positive there was a problem

    double angleOutRad = _calculateAngle(angleInRad, wavelengthNm);

    // Now generate the outgoing vector
    double a = 0;
    double b = 0;

    a = sin(-angleOutRad) / sin(angleInRad);
    b = cos(-angleOutRad) + a*cos(angleInRad);

    vxOut = a*vxIn + b*_vxGratingNormal;
    vyOut = a*vyIn + b*_vyGratingNormal;
    vzOut = a*vzIn + b*_vzGratingNormal;

    return;
}
