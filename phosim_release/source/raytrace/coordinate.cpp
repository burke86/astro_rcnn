void aberrationAxes (double mjd, double *xprime, double *yprime, double *zprime) {

    // Ron-Vondrak Aberration Model in Chap 23 Meeus
    double T = (mjd - 51544.5)/36525.0;
    double L2 = 3.1761467 + 1021.3285546*T;
    double L3 = 1.7534703 + 628.3075849*T;
    double L4 = 6.2034809 + 334.0612431*T;
    double L5 = 0.5995465 + 52.9690965*T;
    double L6 = 0.8740168 + 21.3299095*T;
    double L7 = 5.4812939 + 7.4781599*T;
    double L8 = 5.311863 + 3.813306*T;
    double Lp = 3.8103444 + 8399.6847337*T;
    double D = 5.1984667 + 7771.3771486*T;
    double Mp = 2.3555559 + 8328.6914289*T;
    double F = 1.6279052 + 8433.4661601*T;
    double c = 17314463350.0;

    *xprime = (-1719914 - 2*T)*sin(L3) - 25*cos(L3)
        + (6434 + 141*T)*sin(2*L3) + (28807 - 107*T)*cos(2*L3)
        + 715*sin(L5)
        + 715*sin(Lp)
        + (486 - 5*T)*sin(3*L3) + (-236 - 4*T)*cos(3*L3)
        + 159*sin(L6)
        + 0
        + 39*sin(Lp + Mp)
        + 33*sin(2*L5) - 10*cos(2*L5)
        + 31*sin(2*L3-L5) + 1*cos(2*L3-L5)
        + 8*sin(3*L3-8*L4+3*L5) - 28*cos(3*L3-8*L4+3*L5)
        + 8*sin(5*L3-8*L4+3*L5) - 28*cos(3*L3-8*L4*3*L5)
        + 21*sin(2*L2-L3)
        - 19*sin(L2)
        + 17*sin(L7)
        + 16*sin(L3-2*L5)
        + 16*sin(L8)
        + 11*sin(L3+L5) - 1*cos(L3+L5)
        + 0 - 11*cos(2*L2-2*L3)
        - 11*sin(L3-L5) - 2*cos(L3-L5)
        -7*sin(4*L3) -8*cos(4*L3)
        -10*sin(3*L3-2*L5)
        -9*sin(L2-2*L3)
        -9*sin(2*L2-3*L3)
        -9*cos(2*L6)
        -9*cos(2*L2-4*L3)
        +8*sin(3*L3-2*L4)
        +8*sin(Lp + 2*D - Mp)
        -4*sin(8*L2 - 12*L3) - 7*cos(8*L2-12*L3)
        -4*sin(8*L2 - 14*L3) - 7*cos(8*L2 - 14*L3)
        -6*sin(2*L4) - 5*cos(2*L4)
        -sin(3*L2-4*L3) - cos(3*L2-4*L3)
        +4*sin(2*L3-2*L5) -6*cos(2*L3-2*L5)
        -7*cos(3*L2-3*L3)
        +5*sin(2*L3-2*L4)-5*cos(2*L3-2*L4)
        +5*sin(Lp-2*D);

    *yprime = (25-13*T)*sin(L3) + (1578089+156*T)*cos(L3)
        + (25697-95*T)*sin(2*L3) + (-5904-130*T)*cos(2*L3)
        + 6*sin(L5) - 657*sin(Lp)
        - 656*cos(Lp)
        + (-216-4*T)*sin(3*L3) + (-446+5*T)*cos(3*L3)
        + 2*sin(L6) - 147*cos(L6)
        + 26*cos(F)
        - 36*cos(Lp + Mp)
        - 9*sin(2*L5) - 30*cos(2*L5)
        + 1*sin(2*L3-L5) - 28*cos(2*L3-L5)
        + 25*sin(3*L3-8*L4+3*L5) + 8*cos(3*L3-8*L4+3*L5)
        - 25*sin(5*L3-8*L4+3*L5) - 8*cos(3*L3-8*L4*3*L5)
        - 19*cos(2*L2-L3)
        + 17*cos(L2)
        - 16*cos(L7)
        + 15*cos(L3-2*L5)
        + 1*sin(L8) - 15*cos(L8)
        - 1*sin(L3+L5) - 10*cos(L3+L5)
        - 10*sin(2*L2-2*L3)
        - 2*sin(L3-L5) +9*cos(L3-L5)
        - 8*sin(4*L3) +6*cos(4*L3)
        + 9*cos(3*L3-2*L5)
        - 9*cos(L2-2*L3)
        - 8*cos(2*L2-3*L3)
        - 8*sin(2*L6)
        + 8*sin(2*L2-4*L3)
        - 8*cos(3*L3-2*L4)
        - 7*cos(Lp + 2*D - Mp)
        - 6*sin(8*L2 - 12*L3) + 4*cos(8*L2-12*L3)
        + 6*sin(8*L2 - 14*L3) - 4*cos(8*L2 - 14*L3)
        -4*sin(2*L4) + 5*cos(2*L4)
        -2*sin(3*L2-4*L3) - 7*cos(3*L2-4*L3)
        -5*sin(2*L3-2*L5) -4*cos(2*L3-2*L5)
        -6*sin(3*L2-3*L3)
        -4*sin(2*L3-2*L4)-5*cos(2*L3-2*L4)
        -5*cos(Lp-2*D);

        *zprime = (10+32*T)*sin(L3) + (684185-358*T)*cos(L3)
        + (11141-48*T)*sin(2*L3) + (-2559-55*T)*cos(2*L3)
        - 15*sin(L5) - 282*sin(Lp)
        - 285*cos(Lp)
        + (-94)*sin(3*L3) + (-193)*cos(3*L3)
        - 6*sin(L6) - 61*cos(L6)
        - 59*cos(F)
        - 16*cos(Lp + Mp)
        - 5*sin(2*L5) - 13*cos(2*L5)
        + 0*sin(2*L3-L5) - 12*cos(2*L3-L5)
        + 11*sin(3*L3-8*L4+3*L5) + 3*cos(3*L3-8*L4+3*L5)
        - 11*sin(5*L3-8*L4+3*L5) - 3*cos(3*L3-8*L4*3*L5)
        - 8*cos(2*L2-L3)
        + 8*cos(L2)
        - 7*cos(L7)
        + 1*sin(L3-2*L5) + 7*cos(L3-2*L5)
        - 3*sin(L8) - 6*cos(L8)
        - 1*sin(L3+L5) - 5*cos(L3+L5)
        - 4*sin(2*L2-2*L3)
        - 1*sin(L3-L5) +4*cos(L3-L5)
        - 3*sin(4*L3) +3*cos(4*L3)
        + 4*cos(3*L3-2*L5)
        - 4*cos(L2-2*L3)
        - 4*cos(2*L2-3*L3)
        - 3*sin(2*L6)
        + 3*sin(2*L2-4*L3)
        - 3*cos(3*L3-2*L4)
        - 3*cos(Lp + 2*D - Mp)
        - 3*sin(8*L2 - 12*L3) + 2*cos(8*L2-12*L3)
        + 3*sin(8*L2 - 14*L3) - 2*cos(8*L2 - 14*L3)
        -2*sin(2*L4) + 2*cos(2*L4)
        +1*sin(3*L2-4*L3) - 4*cos(3*L2-4*L3)
        -2*sin(2*L3-2*L5) -2*cos(2*L3-2*L5)
        -3*sin(3*L2-3*L3)
        -2*sin(2*L3-2*L4)-2*cos(2*L3-2*L4)
        -2*cos(Lp-2*D);


    *xprime = *xprime/c;
    *yprime = *yprime/c;
    *zprime = *zprime/c;
}

void aberrationShift (double *ra, double *dec, double xprime, double yprime, double zprime) {

    //Chapter 23 Meeus

    double ra0 = *ra;
    double dec0 = *dec;
    *ra += (yprime*cos(ra0) - xprime*sin(ra0))/(cos(dec0));
    *dec += -((xprime*cos(ra0) + yprime*sin(ra0))*sin(dec0) - zprime*cos(dec0));

}

void precessionShift (double *ra, double *dec, double mjd, double mjd0) {

    //Chapter 21 Meeus
    double T = (mjd0 - 51544.5)/36525.0;
    double t = (mjd - 51544.5)/36525.0;

    double zeta = ((2306.2181 + 1.39656*T - 0.000139*T*T)*t
                   + (0.30188 - 0.000344*T)*t*t + 0.017998*t*t*t)*ARCSEC;
    double z = ((2306.2181 + 1.39656*T - 0.000139*T*T)*t
                + (1.09468 + 0.000066*T)*t*t + 0.018203*t*t*t)*ARCSEC;
    double theta = ((2004.3019 - 0.86330*T - 0.000217*T*T)*t
                    - (0.42665 + 0.000217*T)*t*t - 0.041833*t*t*t)*ARCSEC;

    double A = cos((*dec))*sin((*ra)+zeta);
    double B = cos(theta)*cos((*dec))*cos((*ra) + zeta) - sin(theta)*sin((*dec));
    double C = sin(theta)*cos((*dec))*cos((*ra) + zeta) + cos(theta)*sin((*dec));

    *ra = atan2(A, B) + z;
    if ((*dec > PI/2.0-DEGREE) || (*dec < -PI/2.0+DEGREE)) {
        *dec = acos(sqrt(A*A+B*B));
    } else {
        *dec = asin(C);
    }
}

void nutationShift (double *ra, double *dec, double obliquity, double dobliquity, double dlongitude) {

    double ra0 = *ra;
    double dec0 = *dec;

    *ra += (cos(obliquity*DEGREE) + sin(obliquity*DEGREE)*sin(ra0)*tan(dec0))*dlongitude*ARCSEC
        - (cos(ra0)*tan(dec0))*dobliquity*ARCSEC;
    *dec += (sin(obliquity*DEGREE)*cos(ra0))*dlongitude*ARCSEC + sin(ra0)*dobliquity*ARCSEC;

}

void nutationValues (double mjd, double *dlongitude, double *dobliquity) {
    //Meeus Chapter 22
    //IAU Theory of Nutation (1980)
    //output is in arcseconds

    double t = (mjd - 51544.5)/36525.0;
    double d = 297.85036 + 455267.111480*t - 0.0019142*pow(t,2) + pow(t,3)/189474.0;
    double m = 357.52772 + 35999.050340*t - 0.0001603*pow(t,2) - pow(t,3)/300000.0;
    double mp = 134.96298 + 477198.867398*t + 0.0086972*pow(t,2) + pow(t,3)/56250.0;
    double f = 93.27191 + 483202.017538*t - 0.0036825*pow(t,2) + pow(t,3)/327270.0;
    double omega = 125.04452 - 1934.136261*t + 0.0020708*pow(t,2) + pow(t,3)/450000.0;
    *dlongitude = 0.0001*((-171996-174.2*t)*sin(0*d+0*m+0*mp+0*f+1*omega) +
                             (-13187-1.6*t)*sin(-2*d+0*m+0*mp+2*f+2*omega) +
                              (-2274-0.2*t)*sin(0*d+0*m+0*mp+2*f+2*omega) +
                               (2062+0.2*t)*sin(0*d+0*m+0*mp+0*f+2*omega) +
                               (1426-3.4*t)*sin(0*d+1*m+0*mp+0*f+0*omega) +
                                (712+0.1*t)*sin(0*d+0*m+1*mp+0*f+0*omega) +
                               (-517+1.2*t)*sin(-2*d+1*m+0*mp+2*f+2*omega) +
                               (-386-0.4*t)*sin(0*d+0*m+0*mp+2*f+1*omega) +
                                     (-301)*sin(0*d+0*m+1*mp+2*f+2*omega) +
                                (217-0.5*t)*sin(-2*d-1*m+0*mp+2*f+2*omega) +
                                     (-158)*sin(-2*d+0*m+1*mp+0*f+0*omega) +
                                (129+0.1*t)*sin(-2*d+0*m+0*mp+2*f+1*omega) +
                                      (123)*sin(0*d+0*m-1*mp+2*f+2*omega) +
                                       (63)*sin(2*d+0*m+0*mp+0*f+0*omega) +
                                 (63+0.1*t)*sin(0*d+0*m+1*mp+0*f+1*omega) +
                                      (-59)*sin(2*d+0*m-1*mp+2*f+2*omega) +
                                (-58-0.1*t)*sin(0*d+0*m-1*mp+0*f+1*omega) +
                                      (-51)*sin(0*d+0*m+1*mp+2*f+1*omega) +
                                       (48)*sin(-2*d+0*m+2*mp+0*f+0*omega) +
                                       (46)*sin(0*d+0*m+-2*mp+2*f+1*omega) +
                                      (-38)*sin(0*d+0*m+0*mp+2*f+2*omega) +
                                      (-31)*sin(0*d+0*m+2*mp+2*f+2*omega) +
                                       (29)*sin(0*d+0*m+2*mp+0*f+0*omega) +
                                       (29)*sin(0*d+0*m+1*mp+2*f+2*omega) +
                                       (26)*sin(0*d+0*m+0*mp+2*f+0*omega) +
                                      (-22)*sin(-2*d+0*m+0*mp+2*f+0*omega) +
                                 (17-0.1*t)*sin(0*d+2*m+0*mp+0*f+0*omega) +
                                       (16)*sin(0*d+0*m-1*mp+0*f+1*omega) +
                                (-16+0.1*t)*sin(-2*d+2*m+0*mp+2*f+2*omega) +
                                      (-15)*sin(0*d+1*m+0*mp+0*f+1*omega) +
                                      (-13)*sin(-2*d+0*m+1*mp+0*f+1*omega) +
                                      (-12)*sin(0*d-1*m+0*mp+0*f+1*omega) +
                                       (11)*sin(0*d+0*m+2*mp-2*f+0*omega) +
                                      (-10)*sin(2*d+0*m-1*mp+2*f+1*omega) +
                                       (-8)*sin(2*d+0*m+1*mp+2*f+2*omega) +
                                        (7)*sin(0*d+1*m+0*mp+2*f+2*omega) +
                                       (-7)*sin(-2*d+1*m+1*mp+0*f+0*omega) +
                                       (-7)*sin(0*d-1*m+0*mp+2*f+2*omega) +
                                       (-7)*sin(0*d+1*m+0*mp+2*f+2*omega) +
                                        (6)*sin(2*d+0*m+1*mp+0*f+0*omega) +
                                        (6)*sin(-2*d+0*m+2*mp+2*f+2*omega) +
                                        (6)*sin(-2*d+0*m+1*mp+2*f+1*omega) +
                                       (-6)*sin(2*d+0*m-2*mp+0*f+1*omega) +
                                       (-6)*sin(2*d+0*m+0*mp+0*f+1*omega) +
                                        (5)*sin(0*d-1*m+1*mp+0*f+0*omega) +
                                       (-5)*sin(-2*d-1*m+0*mp+2*f+1*omega) +
                                       (-5)*sin(-2*d+0*m+0*mp+0*f+1*omega) +
                                       (-5)*sin(0*d+0*m+2*mp+2*f+1*omega) +
                                        (4)*sin(-2*d+0*m+2*mp+0*f+1*omega) +
                                        (4)*sin(-2*d+1*m+0*mp+2*f+1*omega) +
                                        (4)*sin(0*d+0*m+1*mp-2*f+0*omega) +
                                       (-4)*sin(-1*d+0*m+1*mp+0*f+0*omega) +
                                       (-4)*sin(-2*d+1*m+0*mp+0*f+0*omega) +
                                       (-4)*sin(1*d+0*m+0*mp+0*f+0*omega) +
                                        (3)*sin(0*d+0*m+1*mp+2*f+0*omega) +
                                       (-3)*sin(0*d+0*m-2*mp+2*f+2*omega) +
                                       (-3)*sin(-1*d-1*m+1*mp+0*f+0*omega) +
                                       (-3)*sin(0*d+1*m+1*mp+0*f+0*omega) +
                                       (-3)*sin(0*d-1*m+1*mp+2*f+2*omega) +
                                       (-3)*sin(2*d-1*m-1*mp+2*f+2*omega) +
                                       (-3)*sin(0*d+0*m+3*mp+2*f+2*omega) +
                                       (-3)*sin(2*d-1*m+0*mp+2*f+2*omega));
    *dobliquity = 0.0001*((92025 + 8.9*t)*cos(0*d+0*m+0*mp+0*f+1*omega) +
                           (5736 - 3.1*t)*cos(-2*d+0*m+0*mp+2*f+2*omega) +
                            (977 - 0.5*t)*cos(0*d+0*m+0*mp+2*f+2*omega) +
                           (-895 + 0.5*t)*cos(0*d+0*m+0*mp+0*f+2*omega) +
                             (54 - 0.1*t)*cos(0*d+1*m+0*mp+0*f+0*omega) +
                                     (-7)*cos(0*d+0*m+1*mp+0*f+0*omega) +
                            (224 - 0.6*t)*cos(-2*d+1*m+0*mp+2*f+2*omega) +
                                    (200)*cos(0*d+0*m+0*mp+2*f+1*omega) +
                            (129 - 0.1*t)*cos(0*d+0*m+1*mp+2*f+2*omega) +
                            (-95 + 0.3*t)*cos(-2*d-1*m+0*mp+2*f+2*omega) +
                                    (-70)*cos(-2*d+0*m+0*mp+2*f+1*omega) +
                                    (-53)*cos(0*d+0*m-1*mp+2*f+2*omega) +
                                    (-33)*cos(0*d+0*m+1*mp+0*f+1*omega) +
                                     (26)*cos(2*d+0*m-1*mp+2*f+2*omega) +
                                     (32)*cos(0*d+0*m-1*mp+0*f+1*omega) +
                                     (27)*cos(0*d+0*m+1*mp+2*f+1*omega) +
                                    (-24)*cos(0*d+0*m-2*mp+2*f+1*omega) +
                                     (16)*cos(2*d+0*m+0*mp+2*f+2*omega) +
                                     (13)*cos(0*d+0*m+2*mp+2*f+2*omega) +
                                    (-12)*cos(-2*d+0*m+1*mp+2*f+2*omega) +
                                    (-10)*cos(0*d+0*m-1*mp+2*f+1*omega) +
                                     (-8)*cos(2*d+0*m-1*mp+0*f+1*omega) +
                                      (7)*cos(-2*d+2*m+0*mp+2*f+2*omega) +
                                      (9)*cos(0*d+1*m+0*mp+0*f+1*omega) +
                                      (7)*cos(-2*d+0*m+1*mp+0*f+1*omega) +
                                      (6)*cos(0*d-1*m+0*mp+0*f+1*omega) +
                                      (5)*cos(2*d+0*m-1*mp+2*f+1*omega) +
                                      (3)*cos(2*d+0*m+1*mp+2*f+2*omega) +
                                      (-3)*cos(0*d+1*m+0*mp+2*f+2*omega) +
                                      (3)*cos(0*d-1*m+0*mp+2*f+2*omega) +
                                      (3)*cos(2*d+0*m+0*mp+2*f+1*omega) +
                                      (-3)*cos(-2*d+0*m+2*mp+2*f+2*omega) +
                                      (-3)*cos(-2*d+0*m+1*mp+2*f+1*omega) +
                                      (3)*cos(2*d+0*m-2*mp+0*f+1*omega) +
                                      (3)*cos(2*d+0*m+0*mp+0*f+1*omega) +
                                      (3)*cos(-2*d-1*m+0*mp+2*f+1*omega) +
                                      (3)*cos(-2*d+0*m+0*mp+0*f+1*omega) +
                                      (3)*cos(0*d+0*m+2*mp+2*f+1*omega));

}

double obliquity (double mjd) {
    // Laskar in Meeus Chapter 22
    // valid +/- 10,000 years around 2000.0
    double t= (mjd - 51544.5)/36525.0;
    double u = t/100;
    return(  (23.0 + 26.0/60.0 + 21.448/3600.0)
             + 1.0/3600.0*(-4680.93*u - 1.55*pow(u,2)
                           +1999.25*pow(u,3) - 51.38*pow(u,4)
                           -249.67*pow(u,5) - 39.05*pow(u,6)
                           +7.12*pow(u,7) + 27.87*pow(u,8)
                           +5.79*pow(u,9) + 2.45*pow(u,10)));
}


double convertTai(double mjd) {
    //http://maia.usno.navy.mil/ser7/tai-utc.dat

    double dt = 0.0;

    if (mjd < 37300.0) { dt=0.0; }
    else if (mjd < 37512.0) { dt = 1.4228180 + (mjd-37300)*0.001296; }
    else if (mjd < 37665.0) { dt = 1.3728180 + (mjd-37300)*0.001296; }
    else if (mjd < 38334.0) { dt = 1.8458580 + (mjd-37665)*0.0011232; }
    else if (mjd < 38395.0) { dt = 1.9458580 + (mjd-37665)*0.0011232; }
    else if (mjd < 38486.0) { dt = 3.2401300 + (mjd-38761)*0.001296; }
    else if (mjd < 38639.0) { dt = 3.3401300 + (mjd-38761)*0.001296; }
    else if (mjd < 38761.0) { dt = 3.4401300 + (mjd-38761)*0.001296; }
    else if (mjd < 38820.0) { dt = 3.5401300 + (mjd-38761)*0.001296; }
    else if (mjd < 38942.0) { dt = 3.6401300 + (mjd-38761)*0.001296; }
    else if (mjd < 39004.0) { dt = 3.7401300 + (mjd-38761)*0.001296; }
    else if (mjd < 39126.0) { dt = 3.8401300 + (mjd-38761)*0.001296; }
    else if (mjd < 39887.0) { dt = 4.3131700 + (mjd-39126)*0.002592; }
    else if (mjd < 41317.0) { dt = 4.2131700 + (mjd-39126)*0.002592; }
    else if (mjd < 41499.0) { dt = 10.0; }
    else if (mjd < 41683.0) { dt = 11.0; }
    else if (mjd < 42048.0) { dt = 12.0; }
    else if (mjd < 42413.0) { dt = 13.0; }
    else if (mjd < 42778.0) { dt = 14.0; }
    else if (mjd < 43144.0) { dt = 15.0; }
    else if (mjd < 43509.0) { dt = 16.0; }
    else if (mjd < 43874.0) { dt = 17.0; }
    else if (mjd < 44239.0) { dt = 18.0; }
    else if (mjd < 44786.0) { dt = 19.0; }
    else if (mjd < 45151.0) { dt = 20.0; }
    else if (mjd < 45516.0) { dt = 21.0; }
    else if (mjd < 46247.0) { dt = 22.0; }
    else if (mjd < 47161.0) { dt = 23.0; }
    else if (mjd < 47892.0) { dt = 24.0; }
    else if (mjd < 48257.0) { dt = 25.0; }
    else if (mjd < 48804.0) { dt = 26.0; }
    else if (mjd < 49169.0) { dt = 27.0; }
    else if (mjd < 49534.0) { dt = 28.0; }
    else if (mjd < 50083.0) { dt = 29.0; }
    else if (mjd < 50630.0) { dt = 30.0; }
    else if (mjd < 51179.0) { dt = 31.0; }
    else if (mjd < 53736.0) { dt = 32.0; }
    else if (mjd < 54832.0) { dt = 33.0; }
    else if (mjd < 56109.0) { dt = 34.0; }
    else if (mjd < 57204.0) { dt = 35.0; }
    else if (mjd < 57754.0) { dt = 36.0; }
    else { dt = 37.0; }

    return (mjd + dt/86400.0);

}
