///
/// @package phosim
/// @file operator.cpp
/// @brief operator class
///
/// @brief Created by
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "../raytrace/constants.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

double localSiderealTime(double mjd, double longitude) {

    // Local Sidereal Time
    // from Mees Equation 12.4  the sidereal time at Greenwich at 0 UT
    // according to 1982 IAU definition expressed in degrees

    double t = (mjd - 51544.5)/36525;
    double theta = 280.46061837 + (360.98564736629 * (mjd - 51544.5)) + pow(t,2)*0.000387933 - pow(t,3)/38710000.0;
    double lst = (theta - longitude)/15.0;
    lst = fmod(lst, 24.0);
    if (lst < 0.0) lst += 24.0;
    return(lst);

}


void altAztoRaDec(double lst, double latitude, double alt, double az, double *ra, double *dec) {
    // alt/az -> ra/dec

    double ha = atan2( -sin(az*DEGREE)*cos(alt*DEGREE),
                       -cos(az*DEGREE)*sin(latitude*DEGREE)*cos(alt*DEGREE)+
                       sin(alt*DEGREE)*cos(latitude*DEGREE))/DEGREE;
    if (ha < 0) ha += 360.0;
    ha = fmod(ha, 360.0);
    *dec = asin(sin(latitude*DEGREE)*sin(alt*DEGREE) + cos(latitude*DEGREE)*cos(alt*DEGREE)*cos(az*DEGREE))/DEGREE;
    *ra = lst*15.0 - ha;
    *ra = fmod(*ra, 360.0);
    if (*ra < 0) *ra += 360.0;
}



void raDectoAltAz (double lst, double latitude, double ra, double dec, double *alt, double *az) {
    // rad/dec -> alt/az
    double ha = lst*15.0 - ra;
    ha = fmod(ha, 360.0);
    if (ha < 0) ha += 360.0;
    *az = atan2( sin(ha*DEGREE), cos(ha*DEGREE)*sin(latitude*DEGREE) -
                 tan(dec*DEGREE)*cos(latitude*DEGREE))/DEGREE + 180.0;
    *alt = asin(sin(latitude*DEGREE)*sin(dec*DEGREE) + cos(latitude*DEGREE)*cos(dec*DEGREE)*cos(ha*DEGREE))/DEGREE;
}





void jdToCalendar (double mjd, long *month, long *day, long *year, double *frac) {
    // Meeus Chapter 7 algorithm to convert JD to Calendar Date
    double jd = mjd + 2400000.5;
    long a, b, c, d, e, f;
    double tjd = jd + 0.5;
    long z = floor (tjd);
    f = tjd - z;
    if (z < 2299161) {
        a = z;
    } else {
        long alpha = floor ((z - 1867216.25)/36524.25);
        a = z + 1 + alpha - floor(alpha/4.0);
    }
    b = a + 1524;
    c = floor((b-122.1)/365.25);
    d = floor(365.25*c);
    e = floor((b-d)/30.6001);
    *day = b - d - floor(30.6001*e) + f;
    if (e < 14) {
        *month = e - 1;
    } else {
        *month = e - 13;
    }
    if (*month > 2) {
        *year = c - 4716;
    } else {
        *year = c - 4715;
    }
    *frac = tjd - z;
    if (*frac > 1.0) {
        *frac = *frac - 1.0;
    }

}


double calculateDayOfYear (long month, long day, long year, double frac) {
    // Meuss Chapter 7 algorithm to convert calander date to day of the year
    int isLeapYear = 0;
    int k;

    if (fmod(year, 100.0) == 0) {
        if (fmod(year, 400.0) == 0) {
            isLeapYear = 1;
        }
    } else if (fmod(year, 4.0) == 0) {
      isLeapYear = 1;
    }

    if (isLeapYear == 1) {
        k = 1;
    } else {
        k = 2;
    }

    return(floor((275*month)/9.0) - k*floor((month + 9)/12.0) + day - 30 + frac);

}

double calculateRotSkyPos (double ha, double latitude, double declination, double rotTelPos) {

    return(atan2(sin(ha*DEGREE)*cos(latitude*DEGREE),sin(latitude*DEGREE)*cos(declination*DEGREE)-sin(declination*DEGREE)*cos(ha*DEGREE)*cos(latitude*DEGREE))/DEGREE-rotTelPos);

}

double calculateRotTelPos (double ha, double latitude, double declination, double rotSkyPos) {

    return(atan2(sin(ha*DEGREE)*cos(latitude*DEGREE),sin(latitude*DEGREE)*cos(declination*DEGREE)-sin(declination*DEGREE)*cos(ha*DEGREE)*cos(latitude*DEGREE))/DEGREE-rotSkyPos);

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




void sunPos (double mjd, double *ra, double *dec) {
    // Sun Position Chapter 25 Meeus
    double t = (mjd -51544.5)/36525.0;
    double l0 = 280.46646 + 36000.76983*t + 0.0003032*pow(t,2);
    double m = 357.52911 + 35999.05029*t - 0.0001537*pow(t,2);
    double c = (1.914602 - 0.004817*t - 0.000015*pow(t,2))*sin(m*DEGREE) +
        (0.019993 - 0.000101*t)*sin(2*m*DEGREE) +
        0.000289*sin(3*m*DEGREE);
    double omega = 125.04 - 1934.136*t;
    double sun = l0 + c;
    //for now will correct for nutation/aberration but maybe not later
    double lambda = sun - 0.00569 - 0.00478*sin(omega*DEGREE);
    double epsilon = obliquity(mjd);
    epsilon += 0.00256*cos(omega);
    *ra = atan2( cos(epsilon*DEGREE)*sin(lambda*DEGREE), cos(lambda*DEGREE))/DEGREE;
    *dec = asin( sin(epsilon*DEGREE)*sin(lambda*DEGREE))/DEGREE;
    if (*ra < 0) *ra += 360.0;
}

void moonPos (double mjd, double *ra, double *dec) {
    double t = (mjd -51544.5)/36525.0;
    double lp = 218.3164477 + 481267.88123421*t
        -0.0015786*pow(t,2)+pow(t,3)/538841.0 - pow(t,4)/65194000.0;
    double d = 297.8501921 + 445267.1114034*t
        -0.0018819*pow(t,2) + pow(t,3)/545868.0 - pow(t,4)/113065000.0;
    double m = 357.5291092 + 35999.0502909*t
        -0.0001536*pow(t,2) + pow(t,3)/24490000.0;
    double mp = 134.9633964 + 477198.8675055*t
        +0.0087414*pow(t,2) + pow(t,3)/69699 - pow(t,4)/14712000.0;
    double f = 93.2720950 + 483202.0175233*t
        -0.0036539*pow(t,2) - pow(t,3)/3526000.0 + pow(t,4)/863310000.0;
    d=d*DEGREE;
    m=m*DEGREE;
    mp=mp*DEGREE;
    f=f*DEGREE;
    double a1 = 119.75 + 131.849*t;
    double a2 = 53.09 + 479264.290*t;
    double a3 = 313.45 + 481266.484*t;
    double e = 1.0 - 0.002516*t - 0.0000074*pow(t,2);
    double lambda = lp
        +1/1000000.0*(6288774*sin(0*d+0*m+1*mp+0*f)+
                      1274027*sin(2*d+0*m-1*mp+0*f)+
                      658314*sin(2*d+0*m+0*mp+0*f)+
                      213618*sin(0*d+0*m+2*mp+0*f)+
                      -185116*e*sin(0*d+1*m+0*mp+0*f)+
                      -114332*sin(0*d+0*m+0*mp+2*f)+
                      58793*sin(2*d+0*m-2*mp+0*f)+
                      57066*e*sin(2*d-1*m-1*mp+0*f)+
                      53322*sin(2*d+0*m+1*mp+0*f)+
                      45758*e*sin(2*d-1*m+0*mp+0*f)+
                      -40923*e*sin(0*d+1*m-1*mp+0*f)+
                      -34720*sin(1*d+0*m+0*mp+0*f)+
                      -30383*e*sin(0*d+1*m+1*mp+0*f)+
                      15327*sin(2*d+0*m+0*mp-2*f)+
                      -12528*sin(0*d+0*m+1*mp+2*f)+
                      10980*sin(0*d+0*m+1*mp-2*f)+
                      10675*sin(4*d+0*m-1*mp+0*f)+
                      10034*sin(0*d+0*m+3*mp+0*f)+
                      8548*sin(4*d+0*m-2*mp+0*f)+
                      -7888*e*sin(2*d+1*m-1*mp+0*f)+
                      -6766*e*sin(2*d+1*m+0*mp+0*f)+
                      -5163*sin(1*d+0*m-1*mp+0*f)+
                      4987*e*sin(1*d+1*m+0*mp+0*f)+
                      4036*e*sin(2*d-1*m+1*mp+0*f)+
                      3994*sin(2*d+0*m+2*mp+0*f)+
                      3861*sin(4*d+0*m+0*mp+0*f)+
                      3665*sin(2*d+0*m-3*mp+0*f)+
                      -2689*e*sin(0*d+1*m-2*mp+0*f)+
                      -2602*sin(2*d+0*m-1*mp+2*f)+
                      2390*e*sin(2*d-1*m-2*mp+0*f)+
                      -2348*sin(1*d+0*m+1*mp+0*f)+
                      2236*e*e*sin(2*d-2*m+0*mp+0*f)+
                      -2120*e*sin(0*d+1*m+2*mp+0*f)+
                      -2069*e*e*sin(0*d+2*m+0*mp+0*f)+
                      2048*e*e*sin(2*d-2*m-1*mp+0*f)+
                      -1773*sin(2*d+0*m+1*mp-2*f)+
                      -1595*sin(2*d+0*m+0*mp+2*f)+
                      1215*e*sin(4*d-1*m-1*mp+0*f)+
                      -1110*sin(0*d+0*m+2*mp+2*f)+
                      -892*sin(3*d+0*m-1*mp+0*f)+
                      -810*e*sin(2*d+1*m+1*mp+0*f)+
                      759*e*sin(4*d-1*m-2*mp+0*f)+
                      -713*e*e*sin(0*d+2*m-1*mp+0*f)+
                      -700*e*e*sin(2*d+2*m-1*mp+0*f)+
                      691*e*sin(2*d+1*m-2*mp+0*f)+
                      596*e*sin(2*d-1*m+0*mp-2*f)+
                      549*sin(4*d+0*m+1*mp+0*f)+
                      537*sin(0*d+0*m+4*mp+0*f)+
                      520*e*sin(4*d-1*m+0*mp+0*f)+
                      -487*sin(1*d+0*m-2*mp+0*f)+
                      -399*e*sin(2*d+1*m+0*mp-2*f)+
                      -381*sin(0*d+0*m+2*mp-2*f)+
                      351*e*sin(1*d+1*m+1*mp+0*f)+
                      -340*sin(3*d+0*m-2*mp+0*f)+
                      330*sin(4*d+0*m-3*mp+0*f)+
                      327*e*sin(2*d-1*m+2*mp+0*f)+
                      -323*e*e*sin(0*d+2*m+1*mp+0*f)+
                      299*e*sin(1*d+1*m-1*mp+0*f)+
                      294*sin(2*d+0*m+3*mp+0*f));
    double beta = 1.0/1000000.0*(5128122*sin(0*d+0*m+0*mp+1*f)+
                               280602*sin(0*d+0*m+1*mp+1*f)+
                               277693*sin(0*d+0*m+1*mp-1*f)+
                               173237*sin(2*d+0*m+0*mp-1*f)+
                               55413*sin(2*d+0*m-1*mp+1*f)+
                               46271*sin(2*d+0*m-1*mp-1*f)+
                               32573*sin(2*d+0*m+0*mp+1*f)+
                               17198*sin(0*d+0*m+2*mp+1*f)+
                               9266*sin(2*d+0*m+1*mp-1*f)+
                               8822*sin(0*d+0*m+2*mp-1*f)+
                               8216*e*sin(2*d-1*m+0*mp-1*f)+
                               4324*sin(2*d+0*m-2*mp-1*f)+
                               4200*sin(2*d+0*m+1*mp+1*f)+
                               -3359*e*sin(2*d+1*m+0*mp-1*f)+
                               2463*e*sin(2*d-1*m-1*mp+1*f)+
                               2211*e*sin(2*d-1*m+0*mp+1*f)+
                               2065*e*sin(2*d-1*m-1*mp-1*f)+
                               -1870*e*sin(0*d+1*m-1*mp-1*f)+
                               1828*sin(4*d+0*m-1*mp-1*f)+
                               -1794*e*sin(0*d+1*m+0*mp+1*f)+
                               -1749*sin(0*d+0*m+0*mp+3*f)+
                               -1565*e*sin(0*d+1*m-1*mp+1*f)+
                               -1491*sin(1*d+0*m+0*mp+1*f)+
                               -1475*e*sin(0*d+1*m+1*mp+1*f)+
                               -1410*e*sin(0*d+1*m+1*mp-1*f)+
                               -1344*e*sin(0*d+1*m+0*mp-1*f)+
                               -1335*sin(1*d+0*m+0*mp-1*f)+
                               1107*sin(0*d+0*m+3*mp+1*f)+
                               1021*sin(4*d+0*m+0*mp-1*f)+
                               833*sin(4*d+0*m-1*mp+1*f)+
                               777*sin(0*d+0*m+1*mp-3*f)+
                               671*sin(4*d+0*m-2*mp+1*f)+
                               607*sin(2*d+0*m+0*mp-3*f)+
                               596*sin(2*d+0*m+2*mp-1*f)+
                               491*e*sin(2*d-1*m+1*mp-1*f)+
                               -451*sin(2*d+0*m-2*mp+1*f)+
                               439*sin(0*d+0*m+3*mp-1*f)+
                               422*sin(2*d+0*m+2*mp+1*f)+
                               421*sin(2*d+0*m-3*mp-1*f)+
                               -366*e*sin(2*d+1*m-1*mp+1*f)+
                               -351*e*sin(2*d+1*m+0*mp+1*f)+
                               331*sin(4*d+0*m+0*mp+1*f)+
                               315*e*sin(2*d-1*m+1*mp+1*f)+
                               302*e*e*sin(2*d-2*m+0*mp-1*f)+
                               -283*sin(0*d+0*m+1*mp+3*f)+
                               -229*e*sin(2*d+1*m+1*mp-1*f)+
                               223*e*sin(1*d+1*m+0*mp-1*f)+
                               223*e*sin(1*d+1*m+0*mp+1*f)+
                               -220*e*sin(0*d+1*m-2*mp-1*f)+
                               -220*e*sin(2*d+1*m-1*mp-1*f)+
                               -185*sin(1*d+0*m+1*mp+1*f)+
                               181*e*sin(2*d-1*m-2*mp-1*f)+
                               -177*e*sin(0*d+1*m+2*mp+1*f)+
                               176*sin(4*d+0*m-2*mp-1*f)+
                               166*e*sin(4*d-1*m-1*mp-1*f)+
                               -164*sin(1*d+0*m+1*mp-1*f)+
                               132*sin(4*d+0*m+1*mp-1*f)+
                               -119*sin(1*d+0*m-1*mp-1*f)+
                               -115*e*sin(4*d-1*m+0*mp-1*f)+
                               107*e*e*sin(2*d-2*m+0*mp+1*f));
    d=d/DEGREE;
    m=m/DEGREE;
    mp=mp/DEGREE;
    f=f/DEGREE;
    lambda += 1.0/1000000.0*(3958*sin(a1*DEGREE) + 1962*sin((lp - f)*DEGREE) + 318*sin(a2*DEGREE));
    beta += 1.0/1000000.0*(-2235*sin(lp*DEGREE) + 382*sin(a3*DEGREE) + 175*sin((a1-f)*DEGREE) +
                           175*sin((a1+f)*DEGREE) + 127*sin((lp-mp)*DEGREE) - 115*sin((lp+mp)*DEGREE));
    double epsilon=obliquity(mjd);
    *ra = atan2(sin(lambda*DEGREE)*cos(epsilon*DEGREE)-tan(beta*DEGREE)*sin(epsilon*DEGREE),
                cos(lambda*DEGREE))/DEGREE;
    *dec = asin(sin(beta*DEGREE)*cos(epsilon*DEGREE) + cos(beta*DEGREE)*sin(epsilon*DEGREE)*sin(lambda*DEGREE))/DEGREE;
    if (*ra < 0) *ra += 360.0;

}

double distSphere(double ra1, double dec1, double ra2, double dec2) {

    double d = sin(dec1*DEGREE)*sin(dec2*DEGREE) + cos(dec1*DEGREE)*cos(dec2*DEGREE)*cos((ra1-ra2)*DEGREE);
    return (acos(d)/DEGREE);

}

double moonPhase(double raSun, double decSun, double raMoon, double decMoon) {
    //Chapter 48 Meeus

    double psi = acos(sin(decSun*DEGREE)*sin(decMoon*DEGREE) + cos(decSun*DEGREE)*cos(decMoon*DEGREE)*cos(raSun*DEGREE - raMoon*DEGREE));
    double i = atan2(EARTH_SUN*sin(psi), EARTH_MOON - EARTH_SUN*cos(psi));
    return(100*(1+cos(i))/2.0);

}
