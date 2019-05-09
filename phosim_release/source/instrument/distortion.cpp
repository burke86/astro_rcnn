///
/// @package phosim
/// @file distortion.cpp
/// @brief distortion
///
/// @brief Created by
/// @author John Peterson (Purdue)
///
/// @brief Modified by
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "distortion.h"
#include "raytrace/constants.h"
#include "ancillary/random.h"
#include <unistd.h>
#include <pthread.h>


double asphere(double r, double radiusofcurv, double height, double conic, double third, double fourth,
               double fifth, double sixth, double seventh, double eighth, double ninth, double tenth) {
    double h;
    if (radiusofcurv != 0) {
        h = (height - (pow(r, 2)/radiusofcurv/(1.0 + sqrt(1.0 - (conic + 1.0)*(pow(r/radiusofcurv, 2))))
         +third*1e3*pow(r, 3) + fourth*1e3*pow(r, 4) + fifth*1e3*pow(r, 5)+
                       sixth*1e3*pow(r, 6) + seventh*1e3*pow(r, 7) + eighth*1e3*pow(r, 8) + ninth*1e3*pow(r, 9) + tenth*1e3*pow(r, 10)))/1e3;
    } else {
        h = height/1e3;
    }
    return h;
}

double asphereDerivative(double r, double radiusofcurv, double height, double conic, double third, double fourth,
               double fifth, double sixth, double seventh, double eighth, double ninth, double tenth) {
    double dhdr;
    if (radiusofcurv != 0) {
        dhdr = (- r/radiusofcurv/(sqrt(1.0 - (conic + 1.0)*(pow(r/radiusofcurv, 2))))
         +3.0*third*1e3*pow(r, 2) + 4.0*fourth*1e3*pow(r, 3) + 5.0*fifth*1e3*pow(r, 4)+
                       6.0*sixth*1e3*pow(r, 5) + 7.0*seventh*1e3*pow(r, 6) +
                8.0*eighth*1e3*pow(r, 7) + 9.0*ninth*1e3*pow(r, 8) + 10.0*tenth*1e3*pow(r, 9));
    } else {
        dhdr = 0.0;
    }
    return dhdr;
}

double annularZernike (double r, double phi, double e, int n) {

    if (n==0) return(1.0);

    if (n==1) return(2.0*r*cos(phi)/sqrt(1+e*e));
    if (n==2) return(2.0*r*sin(phi)/sqrt(1+e*e));

    if (n==3) return(sqrt(3.0)*(2.0*r*r - 1.0 - e*e)/(1.0 - e*e));
    if (n==4) return(sqrt(6.0)*r*r/sqrt(1+e*e+pow(e,4.0))*sin(2.0*phi));
    if (n==5) return(sqrt(6.0)*r*r/sqrt(1+e*e+pow(e,4.0))*cos(2.0*phi));

    return(0.0);
}

double angular (double phi, int n) {

    if (n==0) return(sin(1*phi));
    if (n==1) return(cos(1*phi));
    if (n==2) return(sin(2*phi));
    if (n==3) return(cos(2*phi));
    if (n==4) return(sin(3*phi));
    if (n==5) return(cos(3*phi));
    if (n==6) return(sin(4*phi));
    if (n==7) return(cos(4*phi));
    if (n==8) return(sin(5*phi));
    if (n==9) return(cos(5*phi));
    if (n==10) return(sin(6*phi));
    if (n==11) return(cos(6*phi));
    if (n==12) return(sin(7*phi));
    if (n==13) return(cos(7*phi));
    if (n==14) return(sin(8*phi));
    if (n==15) return(cos(8*phi));
    if (n==16) return(sin(9*phi));
    if (n==17) return(cos(9*phi));
    if (n==18) return(sin(10*phi));
    if (n==19) return(cos(10*phi));
    if (n==20) return(sin(11*phi));
    if (n==21) return(cos(11*phi));
    if (n==22) return(sin(12*phi));
    if (n==23) return(cos(12*phi));

    return(0.0);

}

double legendre (double r, int n) {

    if (n==0) return(sqrt(0.5)*1.0);
    if (n==1) return(sqrt(1.5)*r);
    if (n==2) return(sqrt(2.5)*0.5*(3*r*r-1.0));
    if (n==3) return(sqrt(3.5)*0.5*(5*r*r*r-3.0*r));
    if (n==4) return(sqrt(4.5)*0.125*(35*pow(r,4.0)-30.0*pow(r,2.0)+3.0));
    if (n==5) return(sqrt(5.5)*0.125*(63*pow(r,5)-70.0*pow(r,3.0)+15.0));
    if (n==6) return(sqrt(6.5)*0.0625*(231*pow(r,6)-315*pow(r,4)+105*pow(r,2)-5));
    if (n==7) return(sqrt(7.5)*0.0625*(429*pow(r,7)-693*pow(r,5)+315*pow(r,3)-35*r));
    if (n==8) return(sqrt(8.5)/128*(6435*pow(r,8)-12012*pow(r,6)+6930*pow(r,4)-1260*pow(r,2)+35));
    if (n==9) return(sqrt(9.5)/128*(12155*pow(r,9)-25740*pow(r,7)+18018*pow(r,5)-4620*pow(r,3)+315*r));
    if (n==10) return(sqrt(10.5)/256*(46189*pow(r,10)-109395*pow(r,8)+90090*pow(r,6)-30030*pow(r,4)+3465*pow(r,2)-63));
    if (n==11) return(sqrt(11.5)/256*(88179*pow(r,11)-230945*pow(r,9)+218790*pow(r,7)-90090*pow(r,5)+15015*pow(r,3)-693*r));
    if (n==12) return(sqrt(12.5)/1024*(676039*pow(r,12)-1939938*pow(r,10)+2078505*pow(r,8)-1021020*pow(r,6)+225225*pow(r,4)-18018*pow(r,2)+231));
    if (n==13) return(sqrt(13.5)/1024*(1300075*pow(r,13)-4056234*pow(r,11)+4849845*pow(r,9)-2771340*pow(r,7)+765765*pow(r,5)-90090*pow(r,3)+3003*r));
    if (n==14) return(sqrt(14.5)/2048*(5014575*pow(r,14)-16900975*pow(r,12)+22309287*pow(r,10)-14549535*pow(r,8)+4849845*pow(r,6)-765765*pow(r,4)+
                                       45045*pow(r,2)-429));
    if (n==15) return(sqrt(15.5)/2048*(9694845*pow(r,15)-35102025*pow(r,13)+50702925*pow(r,11)-37182145*pow(r,9)+14549535*pow(r,7)-2909907*pow(r,5)+
                                       255255*pow(r,3)-6435*r));
    if (n==16) return(sqrt(16.5)/32768*(300540195*pow(r,16)-1163381400*pow(r,14)+1825305300*pow(r,12)-1487285800*pow(r,10)+669278610*pow(r,8)-
                                        162954792*pow(r,6)+19399380*pow(r,4)-875160*pow(r,2)+6435));
    if (n==17) return(sqrt(17.5)/32768*(583401555*pow(r,17)-2404321560*pow(r,15)+4071834900*pow(r,13)-3650610600*pow(r,11)+1859107250*pow(r,9)-
                                        535422888*pow(r,7)+81477396*pow(r,5)-5542680*pow(r,3)+109395*r));
    if (n==18) return(sqrt(18.5)/65536*(2268783825*pow(r,18)-9917826435*pow(r,16)+18032411700*pow(r,14)-17644617900*pow(r,12)+
                                        10039179150*pow(r,10)-3346393050*pow(r,8)+624660036*pow(r,6)-58198140*pow(r,4)+
                                        2078505*pow(r,2)-12155));
    if (n==19) return(sqrt(19.5)/65536*(4418157975*pow(r,19)-20419054425*pow(r,17)+39671305740*pow(r,15)-42075627300*pow(r,13)+
                                        26466926850*pow(r,11)-10039179150*pow(r,9)+2230928700*pow(r,7)-267711444*pow(r,5)+14549535*pow(r,3)-
                                        230945*r));
    if (n==20) return(sqrt(20.5)/262144*(34461632205*pow(r,20)-167890003050*pow(r,18)+347123925225*pow(r,16)-396713057400*pow(r,14)+
                                         273491577450*pow(r,12)-116454478140*pow(r,10)+30117537450*pow(r,8)-4461857400*pow(r,6)+
                                         334639305*pow(r,4)-9699690*pow(r,2)+46189));
    if (n==21) return(sqrt(21.5)/262144*(67282234305*pow(r,21) - 344616322050*pow(r,19) + 755505013725*pow(r,17) - 925663800600*pow(r,15) +
                                         694247850450*pow(r,13) - 328189892940*pow(r,11) + 97045398450*pow(r,9) - 17210021400*pow(r,7)+
                                         1673196525*pow(r,5) - 74364290*pow(r,3) + 969969*r));
    if (n==22) return(sqrt(22.5)/524288*(263012370465*pow(r,22) - 1412926920405*pow(r,20) + 3273855059475*pow(r,18) - 4281195077775*pow(r,16) +
                                         3471239252250*pow(r,14) - 1805044411170*pow(r,12) + 601681470390*pow(r,10) - 124772655150*pow(r,8) +
                                         15058768725*pow(r,6) - 929553625*pow(r,4) + 22309287*pow(r,2) - 88179));
    if (n==23) return(sqrt(23.5)/524288*(514589420475*pow(r,23) - 2893136075115*pow(r,21) + 7064634602025*pow(r,19) - 9821565178425*pow(r,17) +
                                         8562390155550*pow(r,15) - 4859734953150*pow(r,13) + 1805044411170*pow(r,11) - 429772478850*pow(r,9) +
                                         62386327575*pow(r,7) - 5019589575*pow(r,5) + 185910725*pow(r,3) - 2028117*r));
    return(0.0);
}

double zernike(double r, double phi, int n) {

    if (n==0) return(1.0);

    if (n==1) return(2.0*r*cos(phi));
    if (n==2) return(2.0*r*sin(phi));

    if (n==3) return(sqrt(3.0)*(2.0*r*r - 1.0));
    if (n==4) return(sqrt(6.0)*r*r*sin(2.0*phi));
    if (n==5) return(sqrt(6.0)*r*r*cos(2.0*phi));

    if (n==6) return(sqrt(8.0)*(3.0*r*r*r - 2.0*r)*sin(phi));
    if (n==7) return(sqrt(8.0)*(3.0*r*r*r - 2.0*r)*cos(phi));
    if (n==8) return(sqrt(8.0)*r*r*r*sin(3.0*phi));
    if (n==9) return(sqrt(8.0)*r*r*r*cos(3.0*phi));

    if (n==10) return(sqrt(5.0)*(6.0*pow(r,4.0)-6.0*pow(r,2.0)+1.0));
    if (n==11) return(sqrt(10.0)*(4.0*pow(r,4.0)-3.0*pow(r,2.0))*cos(2.0*phi));
    if (n==12) return(sqrt(10.0)*(4.0*pow(r,4.0)-3.0*pow(r,2.0))*sin(2.0*phi));
    if (n==13) return(sqrt(10.0)*pow(r,4.0)*cos(4.0*phi));
    if (n==14) return(sqrt(10.0)*pow(r,4.0)*sin(4.0*phi));

    if (n==15) return(sqrt(12.0)*(10.0*pow(r,5.0)-12.0*pow(r,3.0)+3*r)*sin(phi));
    if (n==16) return(sqrt(12.0)*(10.0*pow(r,5.0)-12.0*pow(r,3.0)+3*r)*cos(phi));
    if (n==17) return(sqrt(12.0)*(5.0*pow(r,5.0)-4.0*pow(r,3.0))*sin(3.0*phi));
    if (n==18) return(sqrt(12.0)*(5.0*pow(r,5.0)-4.0*pow(r,3.0))*cos(3.0*phi));
    if (n==19) return(sqrt(12.0)*pow(r,5.0)*sin(5.0*phi));
    if (n==20) return(sqrt(12.0)*pow(r,5.0)*cos(5.0*phi));

    if (n==21) return(sqrt(7.0)*(20.0*pow(r,6.0)-30.0*pow(r,4.0)+12.0*pow(r,2.0)-1.0));
    if (n==22) return(sqrt(14.0)*(15.0*pow(r,6.0)-20.0*pow(r,4.0)+6.0*pow(r,2.0))*cos(2.0*phi));
    if (n==23) return(sqrt(14.0)*(15.0*pow(r,6.0)-20.0*pow(r,4.0)+6.0*pow(r,2.0))*sin(2.0*phi));
    if (n==24) return(sqrt(14.0)*(6.0*pow(r,6.0)-5.0*pow(r,4.0))*cos(4.0*phi));
    if (n==25) return(sqrt(14.0)*(6.0*pow(r,6.0)-5.0*pow(r,4.0))*sin(4.0*phi));
    if (n==26) return(sqrt(14.0)*pow(r,6.0)*cos(6.0*phi));
    if (n==27) return(sqrt(14.0)*pow(r,6.0)*sin(6.0*phi));

    if (n==28) return(4.0*(35.0*pow(r,7.0)-60.0*pow(r,5)+30.0*pow(r,3)-4.0*r)*sin(phi));
    if (n==29) return(4.0*(35.0*pow(r,7.0)-60.0*pow(r,5)+30.0*pow(r,3)-4.0*r)*cos(phi));
    if (n==30) return(4.0*(21.0*pow(r,7.0)-30.0*pow(r,5)+10.0*pow(r,3))*sin(3.0*phi));
    if (n==31) return(4.0*(21.0*pow(r,7.0)-30.0*pow(r,5)+10.0*pow(r,3))*cos(3.0*phi));
    if (n==32) return(4.0*(7.0*pow(r,7.0)-6.0*pow(r,5.0))*sin(5.0*phi));
    if (n==33) return(4.0*(7.0*pow(r,7.0)-6.0*pow(r,5.0))*cos(5.0*phi));
    if (n==34) return(4.0*(pow(r,7.0))*sin(7.0*phi));
    if (n==35) return(4.0*(pow(r,7.0))*cos(7.0*phi));

    if (n==36) return(sqrt(9.0)*(70.0*pow(r,8.0)-140.0*pow(r,6.0)+90.0*pow(r,4.0)-20.0*pow(r,2.0)+1));
    if (n==37) return(sqrt(18.0)*(56.0*pow(r,8.0)-105.0*pow(r,6.0)+60.0*pow(r,4.0)-10.0*pow(r,2.0))*cos(2.0*phi));
    if (n==38) return(sqrt(18.0)*(56.0*pow(r,8.0)-105.0*pow(r,6.0)+60.0*pow(r,4.0)-10.0*pow(r,2.0))*sin(2.0*phi));
    if (n==39) return(sqrt(18.0)*(28.0*pow(r,8.0)-42.0*pow(r,4.0)+15.0*pow(r,4.0))*cos(4.0*phi));
    if (n==40) return(sqrt(18.0)*(28.0*pow(r,8.0)-42.0*pow(r,4.0)+15.0*pow(r,4.0))*sin(4.0*phi));
    if (n==41) return(sqrt(18.0)*(8.0*pow(r,8.0)-7.0*pow(r,6.0))*cos(6.0*phi));
    if (n==42) return(sqrt(18.0)*(8.0*pow(r,8.0)-7.0*pow(r,6.0))*sin(6.0*phi));
    if (n==43) return(sqrt(18.0)*(pow(r,8.0))*cos(8.0*phi));
    if (n==44) return(sqrt(18.0)*(pow(r,8.0))*sin(8.0*phi));

    if (n==45) return(sqrt(20.0)*(126.0*pow(r,9.0)-280.0*pow(r,7)+210.0*pow(r,5)-60.0*pow(r,3)-5.0*r)*sin(phi));
    if (n==46) return(sqrt(20.0)*(126.0*pow(r,9.0)-280.0*pow(r,7)+210.0*pow(r,5)-60.0*pow(r,3)-5.0*r)*cos(phi));
    if (n==47) return(sqrt(20.0)*(84.0*pow(r,9.0)-168.0*pow(r,7)+105.0*pow(r,5)-20.0*pow(r,3))*sin(3.0*phi));
    if (n==48) return(sqrt(20.0)*(84.0*pow(r,9.0)-168.0*pow(r,7)+105.0*pow(r,5)-20.0*pow(r,3))*cos(3.0*phi));
    if (n==49) return(sqrt(20.0)*(36.0*pow(r,9.0)-56.0*pow(r,7)+21.0*pow(r,5))*sin(5.0*phi));
    if (n==50) return(sqrt(20.0)*(36.0*pow(r,9.0)-56.0*pow(r,7)+21.0*pow(r,5))*cos(5.0*phi));
    if (n==51) return(sqrt(20.0)*(9.0*pow(r,9.0)-8.0*pow(r,7))*sin(7.0*phi));
    if (n==52) return(sqrt(20.0)*(9.0*pow(r,9.0)-8.0*pow(r,7))*cos(7.0*phi));
    if (n==53) return(sqrt(20.0)*(pow(r,9.0))*sin(9.0*phi));
    if (n==54) return(sqrt(20.0)*(pow(r,9.0))*cos(9.0*phi));

    if (n==55) return(sqrt(11.0)*(252.0*pow(r,10.0)-630.0*pow(r,8.0)+560.0*pow(r,6)-210.0*pow(r,4)+30.0*pow(r,2)-1.0));
    if (n==56) return(sqrt(22.0)*(210.0*pow(r,10.0)-504.0*pow(r,8.0)+420.0*pow(r,6)-140.0*pow(r,2)+15.0)*cos(2.0*phi));
    if (n==57) return(sqrt(22.0)*(210.0*pow(r,10.0)-504.0*pow(r,8.0)+420.0*pow(r,6)-140.0*pow(r,2)+15.0)*sin(2.0*phi));
    if (n==58) return(sqrt(22.0)*(120.0*pow(r,10.0)-252.0*pow(r,8.0)+168.0*pow(r,6)-35.0*pow(r,4))*cos(4.0*phi));
    if (n==59) return(sqrt(22.0)*(120.0*pow(r,10.0)-252.0*pow(r,8.0)+168.0*pow(r,6)-35.0*pow(r,4))*sin(4.0*phi));
    if (n==60) return(sqrt(22.0)*(45.0*pow(r,10.0)-72.0*pow(r,8.0)+28.0*pow(r,6))*cos(6.0*phi));
    if (n==61) return(sqrt(22.0)*(45.0*pow(r,10.0)-72.0*pow(r,8.0)+28.0*pow(r,6))*sin(6.0*phi));
    if (n==62) return(sqrt(22.0)*(10.0*pow(r,10.0)-9.0*pow(r,8.0))*cos(8.0*phi));
    if (n==63) return(sqrt(22.0)*(10.0*pow(r,10.0)-9.0*pow(r,8.0))*sin(8.0*phi));
    if (n==64) return(sqrt(22.0)*(pow(r,10.0)*cos(10.0*phi)));
    if (n==65) return(sqrt(22.0)*(pow(r,10.0)*sin(10.0*phi)));

    if (n==66) return(sqrt(24.0)*(462*pow(r,11)-1260*pow(r,8)+1260*pow(r,6)-560*pow(r,4)+105*pow(r,2)-6)*sin(phi));
    if (n==67) return(sqrt(24.0)*(462*pow(r,11)-1260*pow(r,8)+1260*pow(r,6)-560*pow(r,4)+105*pow(r,2)-6)*cos(phi));
    if (n==68) return(sqrt(24.0)*(330*pow(r,11)-840*pow(r,9)+756*pow(r,7)-280*pow(r,5)+35*pow(r,3))*sin(3*phi));
    if (n==69) return(sqrt(24.0)*(330*pow(r,11)-840*pow(r,9)+756*pow(r,7)-280*pow(r,5)+35*pow(r,3))*cos(3*phi));
    if (n==70) return(sqrt(24.0)*(165*pow(r,11)-360*pow(r,9)+252*pow(r,7)-56*pow(r,5))*sin(5*phi));
    if (n==71) return(sqrt(24.0)*(165*pow(r,11)-360*pow(r,9)+252*pow(r,7)-56*pow(r,5))*cos(5*phi));
    if (n==72) return(sqrt(24.0)*(55*pow(r,11)-90*pow(r,9)+36*pow(r,7))*sin(7*phi));
    if (n==73) return(sqrt(24.0)*(55*pow(r,11)-90*pow(r,9)+36*pow(r,7))*cos(7*phi));
    if (n==74) return(sqrt(24.0)*(11*pow(r,11)-10*pow(r,9))*sin(9*phi));
    if (n==75) return(sqrt(24.0)*(11*pow(r,11)-10*pow(r,9))*cos(9*phi));
    if (n==76) return(sqrt(24.0)*(pow(r,11))*sin(11*phi));
    if (n==77) return(sqrt(24.0)*(pow(r,11))*cos(11*phi));

    if (n==78) return(sqrt(13.0)*(924*pow(r,12)-2772*pow(r,10)+3150*pow(r,8)-1680*pow(r,6)+420*pow(r,4)-42*pow(r,2)+1));
    if (n==79) return(sqrt(26.0)*(792*pow(r,12)-2310*pow(r,10)+2520*pow(r,8)-1260*pow(r,6)+280*pow(r,4)-21*pow(r,2))*cos(2.0*phi));
    if (n==80) return(sqrt(26.0)*(792*pow(r,12)-2310*pow(r,10)+2520*pow(r,8)-1260*pow(r,6)+280*pow(r,4)-21*pow(r,2))*sin(2.0*phi));
    if (n==81) return(sqrt(26.0)*(495*pow(r,12)-1320*pow(r,10)+1260*pow(r,8)-504*pow(r,6)+70*pow(r,4))*cos(4.0*phi));
    if (n==82) return(sqrt(26.0)*(495*pow(r,12)-1320*pow(r,10)+1260*pow(r,8)-504*pow(r,6)+70*pow(r,4))*sin(4.0*phi));
    if (n==83) return(sqrt(26.0)*(220*pow(r,12)-495*pow(r,10)+360*pow(r,8)-84*pow(r,6))*cos(6.0*phi));
    if (n==84) return(sqrt(26.0)*(220*pow(r,12)-495*pow(r,10)+360*pow(r,8)-84*pow(r,6))*sin(6.0*phi));
    if (n==85) return(sqrt(26.0)*(66*pow(r,12)-110*pow(r,10)+45*pow(r,8))*cos(8.0*phi));
    if (n==86) return(sqrt(26.0)*(66*pow(r,12)-110*pow(r,10)+45*pow(r,8))*sin(8.0*phi));
    if (n==87) return(sqrt(26.0)*(12*pow(r,12)-11*pow(r,10))*cos(10.0*phi));
    if (n==88) return(sqrt(26.0)*(12*pow(r,12)-11*pow(r,10))*sin(10.0*phi));
    if (n==89) return(sqrt(26.0)*(pow(r,12))*cos(12.0*phi));
    if (n==90) return(sqrt(26.0)*(pow(r,12))*sin(12.0*phi));

    return(0.0);

}


int nodeperthread = 256;
int openthread[100];
int openthreads=0;
pthread_mutex_t lock1;

struct dsstruct {
    double *exx;
    double *exy;
    double *exz;
    double *eyy;
    double *eyz;
    double *ezz;
    double *connection;
    double *volume;
    double a3d;
    double kn;
    double ks;
    double x1d;
    double y1d;
    double z1d;
    double *xo, *yo, *zo;
    double *x, *y, *z;
    double *xoo, *yoo, *zoo;
    double *vx, *vy, *vz;
    int *cx, *cy, *cz;
    int *e;
    double dt;
    double grav0, theta, ophi, m;
    double *fraction;
    int *surface;
    double *surfaceError;
    int *surfaceErrorPoint, *surfaceErrorZ;
    double length;
    std::atomic<double> firstDisp, secondDisp;
    double *dhdx, *dhdy;
    double *wfe;
    int maxzern;
    int controlState;
    int control;
    int *actuator;
    std::vector<double> actuatorX;
    std::vector<double> actuatorY;
    std::vector<double> actuatorR;
    double *pert;
    int nleg;
    int nphi;
    int nzern;
    int npert;
    std::atomic<double> vc, vmag;
    int final;
    int actRedo;
    int *actuatorClosestI, *actuatorClosestJ;
    double *acte;
    int mountType, zernikeMode;
    double dzmean, dzcount;
    double alpha;
    double dtdx, tphi,ttheta;
    double averageHeight;
    double temperature0, dtemp;
    double *moveOffset;
    int actuatorStart;
    int zernikeStart;
 };

dsstruct ds;

struct arstruct {
    long i[2048];
    long j[2048];
    long k[2048];
    long N;
    int cc;
    int nt;
};

void* elastic(void *voidArgs) {

    arstruct *ar = (arstruct*)voidArgs;
    double q2 = 0.8;
    double grav = 9.8;

    for (int c=0;c<ar->cc;c++) {

    long index = ar->i[c]*ar->N*ar->N + ar->j[c]*ar->N + ar->k[c];

    if (ds.e[index] == 1) {

        long ip = ar->i[c] + 1;
        if (ip > ar->N - 1) ip = ar->N - 1;
        if (ds.e[ip*ar->N*ar->N + ar->j[c]*ar->N + ar->k[c]] != 1) ip = ar->i[c];
        long im = ar->i[c] - 1;
        if (im < 0) im = 0;
        if (ds.e[im*ar->N*ar->N + ar->j[c]*ar->N + ar->k[c]] != 1) im = ar->i[c];
        long jp = ar->j[c] + 1;
        if (jp > ar->N - 1) jp = ar->N - 1;
        if (ds.e[ar->i[c]*ar->N*ar->N + jp*ar->N + ar->k[c]] != 1) jp = ar->j[c];
        long jm = ar->j[c] - 1;
        if (jm < 0) jm = 0;
        if (ds.e[ar->i[c]*ar->N*ar->N + jm*ar->N + ar->k[c]] != 1) jm = ar->j[c];
        long kp = ar->k[c] + 1;
        if (kp > ar->N - 1) kp = ar->N - 1;
        if (ds.e[ar->i[c]*ar->N*ar->N + ar->j[c]*ar->N + kp] != 1) kp = ar->k[c];
        long km = ar->k[c] - 1;
        if (km < 0) km = 0;
        if (ds.e[ar->i[c]*ar->N*ar->N + ar->j[c]*ar->N + km] != 1) km = ar->k[c];

        ds.exx[index] = 0.0;
        ds.exy[index] = 0.0;
        ds.exz[index] = 0.0;
        ds.eyy[index] = 0.0;
        ds.eyz[index] = 0.0;
        ds.ezz[index] = 0.0;
        long indexip = ip*ar->N*ar->N + ar->j[c]*ar->N + ar->k[c];
        long indexim = im*ar->N*ar->N + ar->j[c]*ar->N + ar->k[c];
        long indexjp = ar->i[c]*ar->N*ar->N + jp*ar->N + ar->k[c];
        long indexjm = ar->i[c]*ar->N*ar->N + jm*ar->N + ar->k[c];
        long indexkp = ar->i[c]*ar->N*ar->N + ar->j[c]*ar->N + kp;
        long indexkm = ar->i[c]*ar->N*ar->N + ar->j[c]*ar->N + km;

        if (im != ip) {
            ds.exx[index] += (((ds.x[indexip] - ds.xo[indexip])-(ds.x[indexim] - ds.xo[indexim]))/(ds.xo[indexip] - ds.xo[indexim]));
            ds.exy[index] += 0.5*(((ds.y[indexip] - ds.yo[indexip])-(ds.y[indexim] - ds.yo[indexim]))/(ds.xo[indexip] - ds.xo[indexim]));
            ds.exz[index] += 0.5*(((ds.z[indexip] - ds.zo[indexip])-(ds.z[indexim] - ds.zo[indexim]))/(ds.xo[indexip] - ds.xo[indexim]));
        }
        if (jm != jp) {
            ds.eyy[index] += (((ds.y[indexjp] - ds.yo[indexjp]) -(ds.y[indexjm] - ds.yo[indexjm]))/(ds.yo[indexjp] - ds.yo[indexjm]));
            ds.exy[index] += 0.5*(((ds.x[indexjp] - ds.xo[indexjp]) -(ds.x[indexjm] - ds.xo[indexjm]))/(ds.yo[indexjp] - ds.yo[indexjm]));
            ds.eyz[index] += 0.5*(((ds.z[indexjp] - ds.zo[indexjp]) -(ds.z[indexjm] - ds.zo[indexjm]))/(ds.yo[indexjp] - ds.yo[indexjm]));
        }
        if (km != kp) {
            ds.ezz[index] += (((ds.z[indexkp] - ds.zo[indexkp]) -(ds.z[indexkm] - ds.zo[indexkm]))/(ds.zo[indexkp] - ds.zo[indexkm]));
            ds.exz[index] += 0.5*(((ds.x[indexkp] - ds.xo[indexkp])-(ds.x[indexkm] - ds.xo[indexkm]))/(ds.zo[indexkp] - ds.zo[indexkm]));
            ds.eyz[index] += 0.5*(((ds.y[indexkp] - ds.yo[indexkp]) -(ds.y[indexkm] - ds.yo[indexkm]))/(ds.zo[indexkp] - ds.zo[indexkm]));
        }

 
        double forcex = 0.0;
        double forcey = 0.0;
        double forcez = 0.0;

        ds.x1d = ds.x[index] - ds.xo[index];
        ds.y1d = ds.y[index] - ds.yo[index];
        ds.z1d = ds.z[index] - ds.zo[index];

        for (int l = 0; l < 18; l++) {

            long ia, ja, ka;

            if (l == 0) {
                ia = ar->i[c] - 1;
                ja = ar->j[c];
                ka = ar->k[c];
            }
            if (l == 1) {
                ia = ar->i[c] + 1;
                ja = ar->j[c];
                ka = ar->k[c];
            }
            if (l == 2) {
                ia = ar->i[c];
                ja = ar->j[c] - 1;
                ka = ar->k[c];
            }
            if (l == 3) {
                ia = ar->i[c];
                ja = ar->j[c] + 1;
                ka = ar->k[c];
            }
            if (l == 4) {
                ia = ar->i[c];
                ja = ar->j[c];
                ka = ar->k[c] - 1;
            }
            if (l == 5) {
                ia = ar->i[c];
                ja = ar->j[c];
                ka = ar->k[c] + 1;
            }
            if (l == 6) {
                ia = ar->i[c] - 1;
                ja = ar->j[c] + 1;
                ka = ar->k[c];
            }
            if (l == 7) {
                ia = ar->i[c] - 1;
                ja = ar->j[c] - 1;
                ka = ar->k[c];
            }
            if (l == 8) {
                ia = ar->i[c] - 1;
                ja = ar->j[c];
                ka = ar->k[c] + 1;
            }
            if (l == 9) {
                ia = ar->i[c] - 1;
                ja = ar->j[c];
                ka = ar->k[c] - 1;
            }
            if (l == 10) {
                ia = ar->i[c];
                ja = ar->j[c] - 1;
                ka = ar->k[c] - 1;
            }
            if (l == 11) {
                ia = ar->i[c];
                ja = ar->j[c] - 1;
                ka = ar->k[c] + 1;
            }
            if (l == 12) {
                ia = ar->i[c];
                ja = ar->j[c] + 1;
                ka = ar->k[c] - 1;
            }
            if (l == 13) {
                ia = ar->i[c];
                ja = ar->j[c] + 1;
                ka = ar->k[c] + 1;
            }
            if (l == 14) {
                ia = ar->i[c] + 1;
                ja = ar->j[c] + 1;
                ka = ar->k[c];
            }
            if (l == 15) {
                ia = ar->i[c] + 1;
                ja = ar->j[c] - 1;
                ka = ar->k[c];
            }
            if (l == 16) {
                ia = ar->i[c] + 1;
                ja = ar->j[c];
                ka = ar->k[c] + 1;
            }
            if (l == 17) {
                ia = ar->i[c] + 1;
                ja = ar->j[c];
                ka = ar->k[c] - 1;
            }

            long indexa = ia*ar->N*ar->N + ja*ar->N + ka;
            if (ia >= 0 && ja >= 0 && ka >= 0 && ia <= ar->N - 1 && ja <= ar->N - 1 && ka <= ar->N - 1) {

                if (ds.e[indexa] == 1) {

                    double r = sqrt(pow(ds.x[index] - ds.x[indexa], 2) +
                                    pow(ds.y[index] - ds.y[indexa], 2) +
                                    pow(ds.z[index] - ds.z[indexa], 2));
                    double ro = sqrt(pow(ds.xo[index] - ds.xo[indexa], 2) +
                                     pow(ds.yo[index] - ds.yo[indexa], 2) +
                                     pow(ds.zo[index] - ds.zo[indexa], 2));

                    double nx = (ds.x[indexa] - ds.x[index])/r;
                    double ny = (ds.y[indexa] - ds.y[index])/r;
                    double nz = (ds.z[indexa] - ds.z[index])/r;

                    double x2d = ds.x[indexa] - ds.xo[indexa];
                    double y2d = ds.y[indexa] - ds.yo[indexa];
                    double z2d = ds.z[indexa] - ds.zo[indexa];

                    double ebxx = 0.5*(ds.exx[index] + ds.exx[indexa]);
                    double ebxy = 0.5*(ds.exy[index] + ds.exy[indexa]);
                    double ebxz = 0.5*(ds.exz[index] + ds.exz[indexa]);
                    double ebyy = 0.5*(ds.eyy[index] + ds.eyy[indexa]);
                    double ebyz = 0.5*(ds.eyz[index] + ds.eyz[indexa]);
                    double ebzz = 0.5*(ds.ezz[index] + ds.ezz[indexa]);

                    double ednx = (ebxx*nx + ebxy*ny + ebxz*nz)*ro;
                    double edny = (ebxy*nx + ebyy*ny + ebyz*nz)*ro;
                    double ednz = (ebxz*nx + ebyz*ny + ebzz*nz)*ro;

                    double edndn = ednx*nx + edny*ny + ednz*nz;

                    double ux = x2d - ds.x1d;
                    double uy = y2d - ds.y1d;
                    double uz = z2d - ds.z1d;
                    double udn = ux*nx + uy*ny + uz*nz;

                    double factor = ds.a3d/ds.connection[index]*ds.volume[index]*3.0/2.0;

                    forcex += ds.kn*udn*nx*factor + ds.ks*(ednx - edndn*nx)*factor;
                    forcey += ds.kn*udn*ny*factor + ds.ks*(edny - edndn*ny)*factor;
                    forcez += ds.kn*udn*nz*factor + ds.ks*(ednz - edndn*nz)*factor;

                }
            }
        }


        double mlocal = ds.m*ds.fraction[index];

        forcex += - mlocal*grav*sin(ds.theta)*cos(ds.ophi);
        forcey += - mlocal*grav*sin(ds.theta)*sin(ds.ophi);
        forcez += - mlocal*grav*cos(ds.theta) + mlocal*ds.grav0*grav;

        ds.vx[index] += forcex/mlocal*ds.dt - q2*fabs(forcex)/mlocal*ds.dt*ds.vx[index]/(fabs(ds.vx[index]) + 1e-10);
        ds.vy[index] += forcey/mlocal*ds.dt - q2*fabs(forcey)/mlocal*ds.dt*ds.vy[index]/(fabs(ds.vy[index]) + 1e-10);
        ds.vz[index] += forcez/mlocal*ds.dt - q2*fabs(forcez)/mlocal*ds.dt*ds.vz[index]/(fabs(ds.vz[index]) + 1e-10);

        if (ds.cx[index] == 0) ds.x[index] += ds.vx[index]*ds.dt;
        if (ds.cy[index] == 0) ds.y[index] += ds.vy[index]*ds.dt;
        if (ds.cz[index] == 0) ds.z[index] += ds.vz[index]*ds.dt;

        if (ds.surface[index] == 1) {
            double rr = sqrt(pow(ds.xoo[index], 2) +
                             pow(ds.yoo[index], 2));
            double pp = atan2(ds.yoo[ar->i[c]*ar->N*ar->N + ar->j[c]*ar->N +ar->k[c]], ds.xoo[ar->i[c]*ar->N*ar->N +ar->j[c]*ar->N + ar->k[c]]);
            double r = rr/ds.length;

            double dz = 0.0;
            for (int l=0; l < ds.maxzern; l++) dz+=ds.wfe[l]*zernike(r,pp,l);

            double disp = pow(ds.z[index] + dz - ds.zoo[index] - ds.dhdx[index]*(ds.x[index]-ds.xoo[index]) -
                              ds.dhdy[index]*(ds.y[index]-ds.yoo[index]), 2.0);

            ds.surfaceError[ar->i[c]*ar->N+ar->j[c]]=ds.z[index] + dz - ds.zoo[index] - ds.dhdx[index]*(ds.x[index]-ds.xoo[index]) -
                ds.dhdy[index]*(ds.y[index]-ds.yoo[index]);
            ds.surfaceErrorPoint[ar->i[c]*ar->N+ar->j[c]]=1;
            ds.surfaceErrorZ[ar->i[c]*ar->N+ar->j[c]]=ar->k[c];
            if (ds.controlState == 1) {
                ds.firstDisp = ds.firstDisp + disp;
            } else {
                ds.secondDisp = ds.firstDisp + disp;
            }
        }
        if (((ds.cz[index] == 1) || (ds.cx[index] == 1) || (ds.cy[index] == 1)) && ((ds.control==1) || (ds.control==3))) {
            double rr = sqrt(pow(ds.xoo[index], 2) +
                             pow(ds.yoo[index], 2));
            double pp = atan2(ds.yoo[index], ds.xoo[index]);
            double sr = rr/ds.length;

            rr=ds.actuatorR[ds.actuator[index]];
            pp=atan2(ds.actuatorY[ds.actuator[index]],ds.actuatorX[ds.actuator[index]]);
            sr = rr/ds.length;

            double dz = 0.0;
            double dx=0.0, dy=0.0;
            long indexu = index;
            long ii=ar->i[c];
            long jj=ar->j[c];

            ii=ds.actuatorClosestI[ds.actuator[index]];
            jj=ds.actuatorClosestJ[ds.actuator[index]];
            if (ds.mountType==0) for (int l = 0; l < ar->N ; l++) if (ds.e[ii*ar->N*ar->N + jj*ar->N + l] == 1) indexu = ii*ar->N*ar->N + jj*ar->N + l;
            if (ds.mountType==1) for (int l = ar->N-1; l >=0 ; l--) if (ds.e[ii*ar->N*ar->N + jj*ar->N + l] == 1) indexu = ii*ar->N*ar->N + jj*ar->N + l;

            dz=ds.dzmean;

            if (ds.zernikeMode == 1) {
                for (int l=0;l<ds.npert-ds.nleg*1-ds.nphi;l++) dz += ds.pert[l]*zernike(sr,pp,l);

                for (int l=ds.npert-ds.nleg*1-ds.nphi;l<ds.npert-ds.nleg*0-ds.nphi;l++) {
                    dz += ds.pert[l]*legendre(sr,l-(ds.npert-ds.nleg*1-ds.nphi));
                }
                for (int l=ds.npert-ds.nphi;l<ds.npert;l++) {
                    dz += ds.pert[l]*angular(pp,l-(ds.npert-ds.nphi));
                }
            }
            if (ds.zernikeMode == 0) {
                if (ds.actuator[index] >= 0) dz += ds.pert[ds.actuator[index]];
                for (int l=ds.npert-ds.nleg*1-ds.nphi;l<ds.npert-ds.nleg*0-ds.nphi;l++) {
                    dz += ds.pert[l]*legendre(sr,l-(ds.npert-ds.nleg*1-ds.nphi));
                }
                for (int l=ds.npert-ds.nphi;l<ds.npert;l++) {
                    dz += ds.pert[l]*angular(pp,l-(ds.npert-ds.nphi));
                }
            }
            if ((ds.final == 1) || (ds.actRedo == 1)) dz += ds.acte[ds.actuator[index]];

            ds.z[index] = ds.zo[index] - dz;
            ds.x[index] = ds.xo[index] - dx;
            ds.y[index] = ds.yo[index] - dy;
        }
        if (((ds.cz[index] == 1) || (ds.cx[index] == 1) || (ds.cy[index] == 1)) && ((ds.control==0) || (ds.control==2))) {
            double rr = sqrt(pow(ds.xoo[index], 2) +
                             pow(ds.yoo[index], 2));
            double pp = atan2(ds.yoo[index], ds.xoo[index]);
            double sr = rr/ds.length;

            rr=ds.actuatorR[ds.actuator[index]];
            pp=atan2(ds.actuatorY[ds.actuator[index]],ds.actuatorX[ds.actuator[index]]);
            sr = rr/ds.length;

            ds.z[index] = ds.zo[index];
            for (int l=0;l<ds.npert-ds.nleg*1-ds.nphi;l++) ds.z[index] += ds.moveOffset[l+ds.zernikeStart]*1e-6*zernike(sr,pp,l);
            if (ds.actuator[index] >= 0) ds.z[index] += ds.moveOffset[ds.actuator[index]+ds.actuatorStart]*1e-6;
            if ((ds.final == 1) || (ds.actRedo == 1)) ds.z[index] += ds.acte[ds.actuator[index]];
            ds.x[index] = ds.xo[index];
            ds.y[index] = ds.yo[index];
        }
        if (ds.cx[index] == 0 && ds.cy[index] == 0 && ds.cz[index] == 0) {
            ds.vmag = ds.vmag + sqrt(pow(ds.vx[index], 2) + pow(ds.vy[index], 2) + pow(ds.vz[index], 2));
            ds.vc = ds.vc + 1.0;
        }
    }
    }

    pthread_mutex_lock(&lock1);
    openthreads--;
    openthread[ar->nt]=0;
    pthread_mutex_unlock(&lock1);
    return NULL;

}


void distortion(int surfaceIndex, int secondSurfaceIndex, long N, double tol, int control, double wfeError, double actError, long seed, const std::string & obsID, int filter,
                const std::string & instrdir, double zenith, double temperature, double azimuth, double temperatureChange, double angTol, double initStep, int zernikeStart, int actuatorStart, double *moveOffset, int pertDebug) {

    double q1 = 0.1;
    double tol2 = 1.0;
    ds.actRedo = 0;

    int warn = 0;
    int *done;
    int *numericalThickness, *numericalLow, *numericalHigh;
    double *actualThickness;
    double *numericalError;
    double *cstep;
    double *dhdr;
    double *actuatorDistance;
    double  *surfaceDerivativeX, *surfaceDerivativeY;
    int  *goodDerivative;
    double totalDerivativePoints = 0.0;
    double avgr0, avgr0n;
    double avgr1, avgr1n;

    // dsstruct ds;

    ds.control=control;
    pthread_t *thread;
    arstruct *args;
    long numthread=8;
    thread = (pthread_t*)malloc(numthread*sizeof(pthread_t));
    args = (arstruct*)malloc(numthread*sizeof(arstruct));
    openthreads=0;
    for (int i=0;i<numthread;i++) openthread[i]=0;
    pthread_mutex_init(&lock1, NULL);

    ds.x = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.y = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.z = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.xo = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.yo = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.zo = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.xoo = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.yoo = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.zoo = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.vx = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.vy = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.vz = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.e = static_cast<int*>(malloc(N*N*N*sizeof(int)));
    ds.cx = static_cast<int*>(malloc(N*N*N*sizeof(int)));
    ds.cy = static_cast<int*>(malloc(N*N*N*sizeof(int)));
    ds.cz = static_cast<int*>(malloc(N*N*N*sizeof(int)));
    done = static_cast<int*>(malloc(N*N*N*sizeof(int)));
    ds.connection = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.volume = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.fraction = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.surface = static_cast<int*>(malloc(N*N*N*sizeof(int)));
    ds.actuator = static_cast<int*>(malloc(N*N*N*sizeof(int)));
    actuatorDistance = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    numericalError = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.surfaceError = static_cast<double*>(malloc(N*N*sizeof(double)));
    ds.surfaceErrorPoint = static_cast<int*>(malloc(N*N*sizeof(int)));
    goodDerivative = static_cast<int*>(malloc(N*N*sizeof(int)));
    ds.surfaceErrorZ = static_cast<int*>(malloc(N*N*sizeof(int)));
    surfaceDerivativeX = static_cast<double*>(malloc(N*N*sizeof(double)));
    surfaceDerivativeY = static_cast<double*>(malloc(N*N*sizeof(double)));
    numericalThickness = static_cast<int*>(malloc(N*N*sizeof(int)));
    actualThickness = static_cast<double*>(malloc(N*N*sizeof(double)));
    numericalLow = static_cast<int*>(malloc(N*N*sizeof(int)));
    numericalHigh = static_cast<int*>(malloc(N*N*sizeof(int)));
    ds.dhdx = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.dhdy = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    dhdr = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.exx = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.exy = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.exz = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.eyy = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.eyz = static_cast<double*>(malloc(N*N*N*sizeof(double)));
    ds.ezz = static_cast<double*>(malloc(N*N*N*sizeof(double)));

    double radiusofcurv[2];
    double height[2];
    double nextHeight[2];
    double nextInner[2];
    double nextOuter[2];
    double averageRadius[2];
    double averageNextRadius[2];
    double outerradius[2];
    double innerradius[2];
    double conic[2];
    double third[2];
    double fourth[2];
    double fifth[2];
    double sixth[2];
    double seventh[2];
    double eighth[2];
    double ninth[2];
    double tenth[2];
    double focus0=1.0;
    double focus1=1.0;
    double focus[2];

    Random random;

    if (seed == -1) random.setSeedFromTime();
    else random.setSeed32(seed);
    random.unwind(10000);


    std::ostringstream opticsFile;
    opticsFile << instrdir << "/optics_" << filter << ".txt";
    readText opticsPars(opticsFile.str());
    int nsurf(0), surfIdx(0);
    double runningz(0);
    std::string materialFile;
    std::string surfacetype;
    int nextFlag=0;
    for (size_t t(0); t < opticsPars.getSize(); t++){
        std::istringstream iss(opticsPars[t]);
        std::string surfaceName, coatingFile, mediumFile;
        double dz, radiusCurvature, inner, outer;
        iss >> surfaceName >> surfacetype >> radiusCurvature >> dz >> outer >> inner;
        runningz += dz;
        //only calculate a pair of surfaces
        if (surfacetype != "none") {
            if (nextFlag==1) {
                nextHeight[surfIdx-1]=runningz;
                nextOuter[surfIdx-1]=outer;
                nextInner[surfIdx-1]=inner;
                nextFlag=0;
            }
            if (nsurf == surfaceIndex || nsurf == secondSurfaceIndex) {
                radiusofcurv[surfIdx] = - radiusCurvature;
                height[surfIdx] = runningz;
                outerradius[surfIdx] = outer;
                innerradius[surfIdx] = inner;
                iss >> conic[surfIdx]
                    >> third[surfIdx]
                    >> fourth[surfIdx]
                    >> fifth[surfIdx]
                    >> sixth[surfIdx]
                    >> seventh[surfIdx]
                    >> eighth[surfIdx]
                    >> ninth[surfIdx]
                    >> tenth[surfIdx];
                iss >> coatingFile >> mediumFile >> materialFile;
                surfIdx++;
                nextFlag=1;
            }
            nsurf++;
        }
    }
    for (int i = 0 ; i < 2 ; i++) {
        averageRadius[i]=1.0/3.0*(pow(outerradius[i],3.0)-pow(innerradius[i],3.0))/
            (1.0/2.0*(pow(outerradius[i],2.0)-pow(innerradius[i],2.0)));
        averageNextRadius[i]=1.0/3.0*(pow(nextOuter[i],3.0)-pow(nextInner[i],3.0))/
            (1.0/2.0*(pow(nextOuter[i],2.0)-pow(nextInner[i],2.0)));
    }

    // read material file
    double density = 0.0;
    double fillFactor = 1.0;
    double youngs = 0.0;
    double nu = 0.0;
    ds.alpha = 0.0;
    ds.temperature0 = 20.0;
    ds.grav0 = 1.0;
    double heatCapacity = 0.0;
    double thickness = 0.0;
    double thermalConductivity = 0.0;
    int mountPoints = 0;
    ds.mountType = 0;
    int mountRim = 0;
    int mountVertical = 0;
    double constraintSize = 0.1;
    double constraintRim =0.0;
    double coolingRate = 0.0;
    int blankType = 0;
    int mountRing =0;
    ds.moveOffset=moveOffset;
    ds.zernikeStart=zernikeStart;
    ds.actuatorStart=actuatorStart;

    std::string sss;
    sss = instrdir + "/" + materialFile;
    std::ifstream inStream(sss.c_str());
    if (inStream) {
        readText materialPars(instrdir + "/" + materialFile);
        for (size_t t(0); t < materialPars.getSize(); t++) {
            std::string line(materialPars[t]);
            readText::get(line, "mountPoints", mountPoints);
            readText::get(line, "mountVertical", mountVertical);
            readText::get(line, "mountRim", mountRim);
            readText::get(line, "mountRing", mountRing);
            readText::get(line, "mountType", ds.mountType);
            readText::get(line, "density", density);
            readText::get(line, "heatCapacity", heatCapacity);
            readText::get(line, "thickness", thickness);
            readText::get(line, "fillFactor", fillFactor);
            readText::get(line, "thermalConductivity", thermalConductivity);
            readText::get(line, "youngsModulus", youngs);
            readText::get(line, "poissonRatio", nu);
            readText::get(line, "nominalTemperature", ds.temperature0);
            readText::get(line, "nominalGravity", ds.grav0);
            readText::get(line, "thermalExpansionCoefficient", ds.alpha);
            readText::get(line, "constraintSize", constraintSize);
            readText::get(line, "constraintRim", constraintRim);
            readText::get(line, "blankType", blankType);
            readText::get(line, "firstFocus", focus0);
            readText::get(line, "secondFocus", focus1);
            readText::get(line, "coolingPerformanceRate", coolingRate);
        }
    }
    if (mountVertical==0) mountVertical=mountPoints;
    if (constraintRim==0.0) constraintRim=constraintSize;
    if (secondSurfaceIndex < 0) height[1]=height[0]-thickness;
    double height0 = (height[0] + height[1]) * 0.5;
    height[0] -= height0;
    height[1] -= height0;
    nextHeight[0] -= height0;
    nextHeight[1] -= height0;
    printf("Height:   %lf %lf\n",height[0],height[1]);
    focus[0]=focus0;
    focus[1]=focus1;

    int fused=0;
    int single=0;
    int doub=0;
    if (secondSurfaceIndex < 0) {
        single=1;
    } else if (ds.mountType==2) {
        doub=1;
    } else {
        fused=1;
    }

    density = density*fillFactor;
    ds.dtemp = temperature;
    ds.theta = zenith;
    ds.ophi = azimuth;

    //have to randomly choose these with timescale
    ds.ttheta = acos(2.0*drand48()-1.0)*180./PI;
    ds.ttheta = PI/2.0;
    ds.tphi = 2*PI*drand48();
    double minzz=1e30;
    double maxzz=-1e30;
    ds.length=0.0;
    double midpoint=0.0;
    double minlengthhigh=0.0;
    double lengthlow=0.0;
    double weightedmidpoint;
    double minlength=0.0;
    int outsurf=0;
    int insurf=0;
    if ((fused == 1) || (doub==1)) {
        if (outerradius[0] > outerradius[1]) ds.length=outerradius[0]/1e3; else ds.length=outerradius[1]/1e3;
        if (outerradius[0] > outerradius[1]) lengthlow=innerradius[0]/1e3; else lengthlow=innerradius[1]/1e3;
        if (outerradius[0] > outerradius[1]) minlength=innerradius[1]/1e3; else minlength=innerradius[0]/1e3;
        if (outerradius[0] > outerradius[1]) minlengthhigh=outerradius[1]/1e3; else minlengthhigh=outerradius[0]/1e3;
        if (outerradius[0] > outerradius[1]) midpoint=0.5*(outerradius[1]/1e3+innerradius[0]/1e3); else midpoint=0.5*(outerradius[0]/1e3+innerradius[1]/1e3);
        if (outerradius[0] > outerradius[1]) { outsurf=0; insurf=1;} else { outsurf=1; insurf=0;}
    } else {
        ds.length=outerradius[0]/1e3;
        minlength=innerradius[0]/1e3;
    }
    weightedmidpoint = sqrt(0.5*(ds.length*ds.length+minlength*minlength));
    printf("Length: %lf %lf %lf\n",ds.length,minlength,midpoint);
    //derived
    double d = 2.0*(ds.length)/float(N - 1); //in m
    ds.a3d = 15.0;
    ds.kn = 3.0*youngs/(1.0 - 2.0*nu)*d/(ds.a3d);
    ds.ks = 3.0*(1.0 - 4.0*nu)*youngs/(1.0 + nu)/(1.0 - 2.0*nu)*d/(ds.a3d);
    ds.m = density*d*d*d;
    ds.dt = sqrt(ds.m/(youngs*d/3.0/(1 - 2.0*nu) + 4.0/3.0*youngs*d/2.0/(1 + nu)))*q1;

    double meanDerivative=0.0;
    double meanRDerivative=0.0;
    double meanRDerivative1=0.0;
    double meanDerivativeA=0.0;
    double meanDerivativeB=0.0;
    double minDerivative=1e30;
    long kk1 = -N;
    long kk2 = N;
    long jj1 = -N;
    long jj2 = N;
    long ii1 = -N;
    long ii2 = N;
    long ia = 0, ja = 0, ka = 0;
    long tpoints=0, apoints=0;
    double hh[2];
    double rh[2];
    double rr=0.0;
    double hmin=0.0, hmax=0.0;
    double xx, yy, zz;
    double avgRadius = 0.0;
    ds.averageHeight = 0.0;
    double hx=0.0, hy=0.0;
    for (long i = 0; i < N; i++) {
        for (long j = 0; j < N; j++) {
            for (long k = 0; k < N; k++) {
                long index=i*N*N + j*N + k;
                ds.e[index] = 0;
                ds.cx[index] = 0;
                ds.cy[index] = 0;
                ds.cz[index] = 0;
                ds.vx[index] = 0.0;
                ds.vy[index] = 0.0;
                ds.vz[index] = 0.0;
                ds.surfaceError[i*N + j] = 0.0;
                ds.surfaceErrorPoint[i*N + j]=0;
                goodDerivative[i*N + j]=0;
                ds.surfaceErrorZ[i*N + j]=0;
                surfaceDerivativeX[i*N + j]=0.0;
                surfaceDerivativeY[i*N + j]=0.0;
                xx = -(ds.length) + 2*(ds.length)*(static_cast<double>(i))/
                    (static_cast<double>(N - 1));
                yy = -(ds.length) + 2*(ds.length)*(static_cast<double>(j))/
                    (static_cast<double>(N - 1));
                zz = -(ds.length) + 2*(ds.length)*(static_cast<double>(k))/
                    (static_cast<double>(N - 1));
                rr = sqrt(pow(xx, 2) + pow(yy, 2))*1e3;
                hh[0]=0.0;
                hh[1]=0.0;
                rh[0]=0.0;
                rh[1]=0.0;
                hx = 0.0;
                hy = 0.0;
                if (fused == 0) {
                    hh[0] = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]);
                    hx= asphereDerivative(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                          sixth[0], seventh[0], eighth[0], ninth[0], tenth[0])*xx/(rr/1e3);
                    hy= asphereDerivative(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                          sixth[0], seventh[0], eighth[0], ninth[0], tenth[0])*yy/(rr/1e3);
                    rh[0]=innerradius[0]/1e3;
                    rh[1]=outerradius[0]/1e3;
                }
                if (fused == 1) {
                    rh[0]=minlength;
                    rh[1]=ds.length;
                    if ((rr/1e3 <= ds.length) && (rr/1e3 > lengthlow)) {
                        hh[0] = asphere(rr, radiusofcurv[outsurf], height[outsurf], conic[outsurf], third[outsurf],
                                        fourth[outsurf], fifth[outsurf], sixth[outsurf], seventh[outsurf],
                                        eighth[outsurf], ninth[outsurf], tenth[outsurf]);
                        hx= asphereDerivative(rr, radiusofcurv[outsurf], height[outsurf], conic[outsurf], third[outsurf],
                                              fourth[outsurf], fifth[outsurf], sixth[outsurf], seventh[outsurf],
                                              eighth[outsurf], ninth[outsurf], tenth[outsurf])*xx/(rr/1e3);
                        hy= asphereDerivative(rr, radiusofcurv[outsurf], height[outsurf], conic[outsurf], third[outsurf],
                                              fourth[outsurf], fifth[outsurf], sixth[outsurf], seventh[outsurf],
                                              eighth[outsurf], ninth[outsurf], tenth[outsurf])*yy/(rr/1e3);
                    }
                    if ((rr/1e3 <= minlengthhigh) && (rr/1e3 > minlength)) {
                        hh[0] = asphere(rr, radiusofcurv[insurf], height[insurf], conic[insurf], third[insurf], fourth[insurf], fifth[insurf],
                                        sixth[insurf], seventh[insurf], eighth[insurf], ninth[insurf], tenth[insurf]);
                        hx= asphereDerivative(rr, radiusofcurv[insurf], height[insurf], conic[insurf], third[insurf],
                                              fourth[insurf], fifth[insurf], sixth[insurf], seventh[insurf],
                                              eighth[insurf], ninth[insurf], tenth[insurf])*xx/(rr/1e3);
                        hy= asphereDerivative(rr, radiusofcurv[insurf], height[insurf], conic[insurf], third[insurf],
                                              fourth[insurf], fifth[insurf], sixth[insurf], seventh[insurf],
                                              eighth[insurf], ninth[insurf], tenth[insurf])*yy/(rr/1e3);
                    }
                    if ((rr/1e3 <= lengthlow) && (rr/1e3 > minlengthhigh)) {
                        hh[0] = 0.5*(asphere(rr, radiusofcurv[insurf], height[insurf], conic[insurf], third[insurf], fourth[insurf], fifth[insurf],
                                        sixth[insurf], seventh[insurf], eighth[insurf], ninth[insurf], tenth[insurf])+
                                     asphere(rr, radiusofcurv[outsurf], height[outsurf], conic[outsurf], third[outsurf],
                                             fourth[outsurf], fifth[outsurf], sixth[outsurf], seventh[outsurf],
                                             eighth[outsurf], ninth[outsurf], tenth[outsurf]));
                        hx= 0.5*(asphereDerivative(rr, radiusofcurv[insurf], height[insurf], conic[insurf], third[insurf],
                                              fourth[insurf], fifth[insurf], sixth[insurf], seventh[insurf],
                                              eighth[insurf], ninth[insurf], tenth[insurf])*xx/(rr/1e3)+
                                 asphereDerivative(rr, radiusofcurv[outsurf], height[outsurf], conic[outsurf], third[outsurf],
                                              fourth[outsurf], fifth[outsurf], sixth[outsurf], seventh[outsurf],
                                                   eighth[outsurf], ninth[outsurf], tenth[outsurf])*xx/(rr/1e3));
                        hy= 0.5*(asphereDerivative(rr, radiusofcurv[insurf], height[insurf], conic[insurf], third[insurf],
                                              fourth[insurf], fifth[insurf], sixth[insurf], seventh[insurf],
                                              eighth[insurf], ninth[insurf], tenth[insurf])*yy/(rr/1e3)+
                                 asphereDerivative(rr, radiusofcurv[outsurf], height[outsurf], conic[outsurf], third[outsurf],
                                              fourth[outsurf], fifth[outsurf], sixth[outsurf], seventh[outsurf],
                                              eighth[outsurf], ninth[outsurf], tenth[outsurf])*yy/(rr/1e3));
                    }
                }

                if (doub==1) {
                    hh[1] = asphere(rr, radiusofcurv[1], height[1], conic[1], third[1], fourth[1], fifth[1],
                                    sixth[1], seventh[1], eighth[1], ninth[1], tenth[1]);
                } else {
                    if (fused==0) {
                    if (blankType == 0) hh[1] = asphere(rr, 0.0, height[1], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                    if (blankType == 1) hh[1] = asphere(rr, radiusofcurv[0], height[1], conic[0], third[0],
                                                        fourth[0], fifth[0], sixth[0], seventh[0], eighth[0],
                                                        ninth[0], tenth[0]);
                    } else {
                    if (blankType == 0) hh[1] = asphere(rr, 0.0, height[0]-thickness, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                    if (blankType == 1) hh[1] = asphere(rr, radiusofcurv[0], height[0]-thickness, conic[0], third[0],
                                                        fourth[0], fifth[0], sixth[0], seventh[0], eighth[0],
                                                        ninth[0], tenth[0]);


                    }
                }
                hmin=0.0;
                hmax=0.0;
                if (hh[1] > hh[0]) { hmax=hh[1]; hmin=hh[0];}
                if (hh[1] <= hh[0]) { hmax=hh[0]; hmin=hh[1];}
                if (rr/1e3 <= rh[1] && rr/1e3 > rh[0] && zz > hmin && zz <= hmax) {
                    ds.fraction[index]=1.0;
                    if (fabs(zz-hmax) < d) {
                         ds.fraction[index]=(hmax-zz)/d + 0.5;
                         zz=hmax;
                    }
                    if (fabs(zz-hmin) < d) {
                         ds.fraction[index]=(zz-hmin)/d + 0.5;
                         zz=hmin;
                    }
                    ds.e[index] = 1;
                    ds.x[index] = xx;
                    ds.y[index] = yy;
                    ds.z[index] = zz;
                    if (ds.z[i*N*N+j*N+k] < minzz) minzz=ds.z[i*N*N+j*N+k];
                    if (ds.z[i*N*N+j*N+k] > maxzz) maxzz=ds.z[i*N*N+j*N+k];
                    ds.xo[index] = xx;
                    ds.yo[index] = yy;
                    ds.zo[index] = zz;
                    ds.xoo[index] = xx;
                    ds.yoo[index] = yy;
                    ds.zoo[index] = zz;
                    ds.dhdx[index] = hx;
                    ds.dhdy[index] = hy;
                    dhdr[index]= sqrt(hx*hx+hy*hy);
                    if (dhdr[index] < minDerivative) minDerivative=dhdr[index];
                    ds.surface[index] = 0;
                    numericalError[index] = hmax - zz;
                    ds.actuator[index] = -1;
                    actuatorDistance[index] = 1e10;
                    avgRadius += sqrt(xx*xx+yy*yy);
                    ds.averageHeight += zz;
                    tpoints++;
                    if (k > kk1) kk1 = k;
                    if (k < kk2) kk2 = k;
                    if (j > jj1) jj1 = j;
                    if (j < jj2) jj2 = j;
                    if (i > ii1) ii1 = i;
                    if (i < ii2) ii2 = i;
                }
                apoints++;
            }
        }
    }
    long ir = ii1 - ii2 + 1;
    long jr = jj1 - jj2 + 1;
    long kr = kk1 - kk2 + 1;
    long *kc, *jc, *ic, *nc, *indc;
    ic = static_cast<long*>(malloc(ir*jr*kr*sizeof(long)));
    jc = static_cast<long*>(malloc(ir*jr*kr*sizeof(long)));
    kc = static_cast<long*>(malloc(ir*jr*kr*sizeof(long)));
    nc = static_cast<long*>(malloc(ir*jr*kr*sizeof(long)));
    indc = static_cast<long*>(malloc(ir*jr*kr*sizeof(long)));
    printf("All points:  %ld\n",apoints);
    printf("Total points:  %ld\n",tpoints);
    printf("Points in reduced lattice:  %ld\n",ir*jr*kr);
    printf("Spacing:  %e\n",d);
    printf("Point mass:  %e\n",ds.m);
    double totalMass=0.0;
    for (long i = ii2; i <= ii1; i++) {
        for (long j = jj2; j <= jj1; j++) {
            for (long k = kk2; k <= kk1; k++) {
                if (ds.e[i*N*N + j*N +k]==1) {
                    totalMass += ds.m*ds.fraction[i*N*N + j*N +k];
                }
            }
        }
    }
    printf("Mass:  %e\n",totalMass);
    avgRadius=avgRadius/(static_cast<double>(tpoints));
    ds.averageHeight=ds.averageHeight/(static_cast<double>(tpoints));
    printf("Average radius:  %lf\n",avgRadius);
    printf("Average height:  %lf\n",ds.averageHeight);
    printf("Min derivative:   %e\n",minDerivative);

    // connecting points
    for (long i = ii2; i <= ii1; i++) {
        for (long j = jj2; j <= jj1; j++) {
            for (long k = kk2; k <= kk1; k++) {
                if (ds.e[i*N*N + j*N + k] == 1) {
                    ds.connection[i*N*N + j*N + k] = 0;
                    ds.volume[i*N*N + j*N + k] = 0;
                    double sidea=0.0,sideb=0.0,sidec=0.0;
                    for (int l = 0; l < 18; l++) {

                        if (l == 0) {
                            ia = i - 1;
                            ja = j;
                            ka = k;
                        }
                        if (l == 1) {
                            ia = i + 1;
                            ja = j;
                            ka = k;
                        }
                        if (l == 2) {
                            ia = i;
                            ja = j - 1;
                            ka = k;
                        }
                        if (l == 3) {
                            ia = i;
                            ja = j + 1;
                            ka = k;
                        }
                        if (l == 4) {
                            ia = i;
                            ja = j;
                            ka = k - 1;
                        }
                        if (l == 5) {
                            ia = i;
                            ja = j;
                            ka = k + 1;
                        }
                        if (l == 6) {
                            ia = i - 1;
                            ja = j + 1;
                            ka = k;
                        }
                        if (l == 7) {
                            ia = i - 1;
                            ja = j - 1;
                            ka = k;
                        }
                        if (l == 8) {
                            ia = i - 1;
                            ja = j;
                            ka = k + 1;
                        }
                        if (l == 9) {
                            ia = i - 1;
                            ja = j;
                            ka = k - 1;
                        }
                        if (l == 10) {
                            ia = i;
                            ja = j - 1;
                            ka = k - 1;
                        }
                        if (l == 11) {
                            ia = i;
                            ja = j - 1;
                            ka = k + 1;
                        }
                        if (l == 12) {
                            ia = i;
                            ja = j + 1;
                            ka = k - 1;
                        }
                        if (l == 13) {
                            ia = i;
                            ja = j + 1;
                            ka = k + 1;
                        }
                        if (l == 14) {
                            ia = i + 1;
                            ja = j + 1;
                            ka = k;
                        }
                        if (l == 15) {
                            ia = i + 1;
                            ja = j - 1;
                            ka = k;
                        }
                        if (l == 16) {
                            ia = i + 1;
                            ja = j;
                            ka = k + 1;
                        }
                        if (l == 17) {
                            ia = i + 1;
                            ja = j;
                            ka = k - 1;
                        }
                        if (ia >= 0 && ja >= 0 && ka >= 0 && ia <= N - 1 && ja <= N - 1 && ka <= N - 1) {
                            if (ds.e[ia*N*N + ja*N +ka]==1) {
                                double distance = sqrt((ia-i)*(ia-i) + (ja - j)*(ja - j) + (ka - k)*(ka - k));
                                ds.connection[i*N*N + j*N + k] += 0.5*distance*distance;
                                if (l==0) sidea+=0.5;
                                if (l==1) sidea+=0.5;
                                if (l==2) sideb+=0.5;
                                if (l==3) sideb+=0.5;
                                if (l==4) sidec+=ds.fraction[i*N*N + j*N + k]/2.0;
                                if (l==5) sidec+=ds.fraction[i*N*N + j*N + k]/2.0;
                            }
                        }
                    }
                    ds.volume[i*N*N + j*N + k] = sidea*sideb*sidec;
                }
            }
        }
    }


    // surface points
    long surfacePoints = 0;
    int top=0, bottom=0;
    for (long i = ii2; i <= ii1; i++) {
        for (long j = jj2; j <= jj1; j++) {
            numericalLow[i*N+j]=-1;
            numericalHigh[i*N+j]=-1;
            for (long k = kk2; k <= kk1; k++) {
                long idx = i*N*N + j*N + k;
                if (ds.e[idx] == 1) {
                    top = 0;
                    bottom = 0;
                    if (k + 1 < N) {
                        if (ds.e[idx + 1] == 0) top = 1;
                    } else top = 1;
                    if (k > 0) {
                        if (ds.e[idx - 1] == 0) bottom = 1;
                    } else bottom = 1;
                    if (top == 1 && bottom == 1) {
                        if (warn==0) std::cout << "Warning: single layer" << std::endl;
                        warn++;
                    }
                    if (bottom == 1) numericalLow[i*N+j]=k;
                    if (top == 1) numericalHigh[i*N+j]=k;
                    if (top == 1 || bottom == 1) {

                        double x0 = ds.xoo[idx] * 1e3; // in mm
                        double y0 = ds.yoo[idx] * 1e3;

                        if (single == 1) {
                            if ((ds.mountType == 0) && (top == 1)) {
                                ds.surface[idx] = 1;
                                meanDerivative += dhdr[idx];
                                surfacePoints++;
                            } else if ((ds.mountType == 1) && (bottom == 1)) {
                                ds.surface[idx] = 1;
                                meanDerivative += dhdr[idx];
                                surfacePoints++;
                            }
                        }
                        if (doub == 1) {
                            if (bottom == 1) {
                                ds.surface[idx] = 1;
                                meanDerivative += dhdr[idx];
                                surfacePoints++;
                            }
                            if (top == 1) {
                                ds.surface[idx] = 1;
                                meanDerivative += dhdr[idx];
                                surfacePoints++;
                            }

                        }
                        if (fused == 1) {
                            if ((ds.mountType == 0) && (top == 1)) {
                                double rr = sqrt(pow(x0, 2) + pow(y0, 2));
                                if ((rr <= outerradius[0]) && (rr >= innerradius[0])) {
                                    ds.surface[idx] = 1;
                                    meanDerivative += dhdr[idx];
                                    meanDerivativeA += dhdr[idx];
                                    surfacePoints++;
                                }
                                if ((rr <= outerradius[1]) && (rr >= innerradius[1])) {
                                    ds.surface[idx] = 1;
                                    meanDerivative += dhdr[idx];
                                    meanDerivativeB += dhdr[idx];
                                    surfacePoints++;
                                }
                            }
                            if ((ds.mountType == 1) && (bottom == 1)) {
                                double rr = sqrt(pow(x0, 2) + pow(y0, 2));
                                if ((rr <= outerradius[0]) && (rr >= innerradius[0])) {
                                    ds.surface[idx] = 1;
                                    meanDerivative += dhdr[idx];
                                    surfacePoints++;
                                }
                                if ((rr <= outerradius[1]) && (rr >= innerradius[1])) {
                                    ds.surface[idx] = 1;
                                    meanDerivative += dhdr[idx];
                                    surfacePoints++;
                                }
                            }
                        }
                    }
                }
            }
            numericalThickness[i*N+j]=numericalHigh[i*N+j]-numericalLow[i*N+j]+1;
            actualThickness[i*N+j]=ds.zoo[i*N*N+j*N+numericalHigh[i*N+j]]-ds.zoo[i*N*N+j*N+numericalLow[i*N+j]];
            if ((numericalLow[i*N+j] < 0) || (numericalHigh[i*N+j] < 0)) numericalThickness[i*N+j]=-1;
        }
    }

    meanDerivative /= static_cast<double>(surfacePoints);
    meanDerivativeA /= static_cast<double>(surfacePoints);
    meanDerivativeB /= static_cast<double>(surfacePoints);
    printf("Mean derivative:  %e\n",meanDerivative);
    // printf("Mean derivative A:  %e\n",meanDerivativeA);
    // printf("Mean derivative B:  %e\n",meanDerivativeB);

    int aa[N];
    double maxa[N];
    double mina[N];
    for (long i=0;i<N;i++) {
        aa[i]=0;
        maxa[i]=0.0;
        mina[i]=1e10;
    }
    double meanThickness=0.0;
    double meanThicknessCounter=0.0;
    for (long i = ii2; i <= ii1; i++) {
        for (long j = jj2; j <= jj1; j++) {
            if (numericalThickness[i*N+j]-1 > 0) {
                meanThickness += actualThickness[i*N+j];
                meanThicknessCounter += 1.0;
                aa[numericalThickness[i*N+j]-1]++;
                double radius = sqrt((i - (N/2-0.5))*(i - (N/2-0.5)) + (j - (N/2-0.5))*(j - (N/2-0.5)));
                if (radius > maxa[numericalThickness[i*N+j]-1]) maxa[numericalThickness[i*N+j]-1]=radius;
                if (radius < mina[numericalThickness[i*N+j]-1]) mina[numericalThickness[i*N+j]-1]=radius;
            }
        }
    }
    meanThickness /= meanThicknessCounter;
    printf("Thickness:   %e %e\n",thickness,meanThickness);
    printf("Cooling effective conductivity:   %e\n",coolingRate*meanThickness/ds.length/ds.length/PI);
    double timescale = meanThickness*meanThickness/(thermalConductivity+coolingRate*meanThickness/ds.length/ds.length/PI)
        *density*heatCapacity/PI/PI;
    double dTdt = temperatureChange;
    ds.dtdx = timescale*dTdt/3600.0/ds.length;
    printf("Thermal Timescale = %lf hours\n",timescale/3600.0);
    printf("Thermal Gradient = %lf deg/m\n",ds.dtdx);

    int thicknessPoint=0;
    int thicknessValue[N];
    long j=0;
    for (long i=0;i<N;i++) {
        if (aa[i]!=0) {
            std::cout << "Thickness " << i+1 << " " << aa[i] << " " << mina[i] << " " << maxa[i] << std::endl;
            thicknessValue[j]=i+1;
            thicknessPoint++;
            j++;
        }
    }
    std::cout << "Surface points:  " << surfacePoints << std::endl;


    //constraints
    long nactuator = 0;
    int *actuatorClosestK;
    ds.actuatorClosestI = static_cast<int*>(malloc(((mountPoints+1)*(mountPoints+1)+mountRim+mountRing*mountRing*7)*sizeof(int)));
    ds.actuatorClosestJ = static_cast<int*>(malloc(((mountPoints+1)*(mountPoints+1)+mountRim+mountRing*mountRing*7)*sizeof(int)));
    actuatorClosestK = static_cast<int*>(malloc(((mountPoints+1)*(mountPoints+1)+mountRim+mountRing*mountRing*7)*sizeof(int)));
    double *actuatorClosestDistance;
    actuatorClosestDistance = static_cast<double*>(malloc(((mountPoints+1)*(mountPoints+1)+mountRim+mountRing*mountRing*7)*sizeof(double)));
    for (long i=0;i<((mountPoints+1)*(mountPoints+1)+mountRim+mountRing*mountRing*7);i++) actuatorClosestDistance[i]=1e30;
    long constraint =0;
    
    if ((ds.mountType==0) || (ds.mountType==1)) {
        // pushed from bottom
        long NC = mountPoints;
        long NCa = mountVertical;
        for (long ja=0; ja<NCa; ja++) {
            long ncc = 0;
            if ((ja % 2)==1) ncc=NC; else ncc=NC-1;
            for (long ia=0; ia<ncc; ia++) {
                double dx =0.0;
                if ((ja % 2)==1) dx=0.0; else dx=0.5;
                double xx=-(ds.length)+2*(ds.length)*((double)ia+dx)/((double)(NC-1));
                double yy=-(ds.length)+2*(ds.length)*((double)ja)/((double)(NCa-1));
                double rr=sqrt(xx*xx+yy*yy);
                if ((rr < ds.length) && (rr >minlength)) {
                    double mindist=1e30;
                    long ib, jb, kb;
                    int found=0;
                    int foundT=0;
                    for (long i=0; i<N; i++) {
                        for (long j=0; j<N; j++) {
                            for (long k=0; k<N; k++) {
                                if (ds.e[i*N*N+j*N+k]==1) {
                                    double dist;
                                    if (ds.mountType==0) {
                                        dist=sqrt(pow(xx-ds.x[i*N*N+j*N+k],2.0)+pow(yy-ds.y[i*N*N+j*N+k],2.0)+
                                                     pow(minzz-ds.z[i*N*N+j*N+k],2.0));
                                    } else {
                                        dist=sqrt(pow(xx-ds.x[i*N*N+j*N+k],2.0)+pow(yy-ds.y[i*N*N+j*N+k],2.0)+
                                                     pow(maxzz-ds.z[i*N*N+j*N+k],2.0));
                                    }
                                    if (dist < mindist && dist >= constraintSize) {
                                        mindist=dist;
                                        ib=i;
                                        jb=j;
                                        kb=k;
                                        found=1;
                                    }
                                    top = 0;
                                    bottom = 0;
                                    int ok = 0;
                                    long idx = i*N*N + j*N + k;
                                    if (k + 1 < N) {
                                        if (ds.e[idx + 1] == 0) top = 1;
                                    } else top = 1;
                                    if (k > 0) {
                                        if (ds.e[idx - 1] == 0) bottom = 1;
                                    } else bottom = 1;
                                    if (ds.mountType == 0 && bottom == 1) ok = 1;
                                    if (ds.mountType == 1 && top == 1) ok = 1;
                                    if (dist < constraintSize && ok == 1) {
                                        foundT=1;
                                        mindist=dist;
                                        if (ds.cx[i*N*N+j*N+k]==0) {
                                            ds.cx[i*N*N+j*N+k]=1;
                                            ds.cy[i*N*N+j*N+k]=1;
                                            ds.cz[i*N*N+j*N+k]=1;
                                            if (dist < actuatorDistance[i*N*N+j*N+k]) {
                                                ds.actuator[i*N*N+j*N+k] = nactuator;
                                                actuatorDistance[i*N*N+j*N+k] = dist;
                                            }
                                            if (dist < actuatorClosestDistance[nactuator]) {
                                                actuatorClosestDistance[nactuator]=dist;
                                                ds.actuatorClosestI[nactuator]=i;
                                                ds.actuatorClosestJ[nactuator]=j;
                                                actuatorClosestK[nactuator]=k;
                                            }
                                            constraint++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (found==0) {
                        printf("Error:  couldn't find constraint.\n");
                    } else {
                        if (foundT == 0) {
                            if (ds.cx[ib*N*N+jb*N+kb]==0) {
                                ds.cx[ib*N*N+jb*N+kb]=1;
                                ds.cy[ib*N*N+jb*N+kb]=1;
                                ds.cz[ib*N*N+jb*N+kb]=1;
                                if (mindist < actuatorDistance[ib*N*N+jb*N+kb]) {
                                    ds.actuator[ib*N*N+jb*N+kb] = nactuator;
                                    actuatorDistance[ib*N*N+jb*N+kb] = mindist;
                                }
                                actuatorClosestDistance[nactuator]=mindist;
                                ds.actuatorClosestI[nactuator]=ib;
                                ds.actuatorClosestJ[nactuator]=jb;
                                actuatorClosestK[nactuator]=kb;
                                    constraint++;
                            }
                        }
                    }
                    ds.actuatorX.push_back(xx);
                    ds.actuatorY.push_back(yy);
                    ds.actuatorR.push_back(rr);
                    nactuator++;
                }
            }
        }
    }
    // ring mount
    if (((ds.mountType==0) || (ds.mountType==1)) && (mountRing>0)) {
        // pushed from bottom
        long NC = mountRing;
        for (long ja=0; ja<NC; ja++) {
            long ncc = ceill(2*PI*((double)ja)/((double)(NC-1))*NC);
            for (long ia=0; ia<ncc; ia++) {
                double rr=(ds.length)*((double)ja)/((double)(NC-1));
                double phi=2*PI*((double)ia)/((double)(ncc));
                double xx=rr*cos(phi);
                double yy=rr*sin(phi);
                if (rr >=minlength) {
                   double mindist=1e30;
                    long ib, jb, kb;
                    int found=0;
                    int foundT=0;
                    for (long i=0; i<N; i++) {
                        for (long j=0; j<N; j++) {
                            for (long k=0; k<N; k++) {
                                if (ds.e[i*N*N+j*N+k]==1) {
                                    double dist;
                                    if (ds.mountType==0) {
                                        dist=sqrt(pow(xx-ds.x[i*N*N+j*N+k],2.0)+pow(yy-ds.y[i*N*N+j*N+k],2.0)+
                                                     pow(minzz-ds.z[i*N*N+j*N+k],2.0));
                                    } else {
                                        dist=sqrt(pow(xx-ds.x[i*N*N+j*N+k],2.0)+pow(yy-ds.y[i*N*N+j*N+k],2.0)+
                                                     pow(maxzz-ds.z[i*N*N+j*N+k],2.0));
                                    }
                                    if (dist < mindist && dist >= constraintSize) {
                                        mindist=dist;
                                        ib=i;
                                        jb=j;
                                        kb=k;
                                        found=1;
                                    }
                                    top = 0;
                                    bottom = 0;
                                    int ok = 0;
                                    long idx = i*N*N + j*N + k;
                                    if (k + 1 < N) {
                                        if (ds.e[idx + 1] == 0) top = 1;
                                    } else top = 1;
                                    if (k > 0) {
                                        if (ds.e[idx - 1] == 0) bottom = 1;
                                    } else bottom = 1;
                                    if (ds.mountType == 0 && bottom == 1) ok = 1;
                                    if (ds.mountType == 1 && top == 1) ok = 1;
                                    if (dist < constraintSize && ok == 1) {
                                        foundT=1;
                                        mindist=dist;
                                        if (ds.cx[i*N*N+j*N+k]==0) {
                                            ds.cx[i*N*N+j*N+k]=1;
                                            ds.cy[i*N*N+j*N+k]=1;
                                            ds.cz[i*N*N+j*N+k]=1;
                                            if (dist < actuatorDistance[i*N*N+j*N+k]) {
                                                ds.actuator[i*N*N+j*N+k] = nactuator;
                                                actuatorDistance[i*N*N+j*N+k] = dist;
                                            }
                                            if (dist < actuatorClosestDistance[nactuator]) {
                                                actuatorClosestDistance[nactuator]=dist;
                                                ds.actuatorClosestI[nactuator]=i;
                                                ds.actuatorClosestJ[nactuator]=j;
                                                actuatorClosestK[nactuator]=k;
                                            }
                                            constraint++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (found==0) {
                        printf("Error:  couldn't find constraint.\n");
                    } else {
                        if (foundT == 0) {
                            if (ds.cx[ib*N*N+jb*N+kb]==0) {
                                ds.cx[ib*N*N+jb*N+kb]=1;
                                ds.cy[ib*N*N+jb*N+kb]=1;
                                ds.cz[ib*N*N+jb*N+kb]=1;
                                if (mindist < actuatorDistance[ib*N*N+jb*N+kb]) {
                                    ds.actuator[ib*N*N+jb*N+kb] = nactuator;
                                    actuatorDistance[ib*N*N+jb*N+kb] = mindist;
                                }
                                actuatorClosestDistance[nactuator]=mindist;
                                ds.actuatorClosestI[nactuator]=ib;
                                ds.actuatorClosestJ[nactuator]=jb;
                                actuatorClosestK[nactuator]=kb;
                                    constraint++;
                            }
                        }
                    }
                    ds.actuatorX.push_back(xx);
                    ds.actuatorY.push_back(yy);
                    ds.actuatorR.push_back(rr);
                    nactuator++;
                }
            }
        }
    }


    if (ds.mountType==2 || (mountRim>0)) {
        // contrained from side
        long mountSide;
        if (ds.mountType==2) mountSide=mountPoints;
        if (mountRim>0) mountSide=mountRim;
        if (mountRim>0) constraintSize=constraintRim;
         for (long l=0; l<mountSide; l++) {
             double phi = 2*PI*static_cast<double>(l)/(static_cast<double>(mountSide));
             double xx = ds.length*cos(phi);
             double yy = ds.length*sin(phi);
            double mindist=1e30;
            long ib, jb, kb;
            int found=0;
            int foundT=0;
            for (long i=0; i<N; i++) {
                for (long j=0; j<N; j++) {
                    for (long k=0; k<N; k++) {
                        double midzz = 0.5*(maxzz+minzz);
                        double dist;
                        if (ds.mountType==0) {
                            dist=sqrt(pow(xx-ds.x[i*N*N+j*N+k],2.0)+pow(yy-ds.y[i*N*N+j*N+k],2.0)+
                                      pow(minzz-ds.z[i*N*N+j*N+k],2.0));
                        } else if (ds.mountType==1) {
                            dist=sqrt(pow(xx-ds.x[i*N*N+j*N+k],2.0)+pow(yy-ds.y[i*N*N+j*N+k],2.0)+
                                      pow(maxzz-ds.z[i*N*N+j*N+k],2.0));
                        } else {
                            dist=sqrt(pow(xx-ds.x[i*N*N+j*N+k],2.0)+pow(yy-ds.y[i*N*N+j*N+k],2.0)+
                                      pow(midzz-ds.z[i*N*N+j*N+k],2.0));
                        }
                        if (dist < mindist && dist > constraintSize) {
                            mindist=dist;
                            ib=i;
                            jb=j;
                            kb=k;
                            found=1;
                        }
                        int ok =0;
                        double rr = sqrt(pow(ds.x[i*N*N+j*N+k],2.0)+pow(ds.y[i*N*N+j*N+k],2.0));
                        if ((ds.mountType == 2) && (rr > (ds.length - d))) ok = 1;
                        top = 0;
                        bottom = 0;
                        long idx = i*N*N + j*N + k;
                        if (k + 1 < N) {
                            if (ds.e[idx + 1] == 0) top = 1;
                        } else top = 1;
                        if (k > 0) {
                            if (ds.e[idx - 1] == 0) bottom = 1;
                        } else bottom = 1;
                        if (ds.mountType == 0 && bottom == 1) ok = 1;
                        if (ds.mountType == 1 && top == 1) ok = 1;
                        if (dist < constraintSize && ok==1) {
                            foundT=1;
                            mindist=dist;
                            if (ds.cx[i*N*N+j*N+k]==0) {
                                ds.cx[i*N*N+j*N+k]=1;
                                ds.cy[i*N*N+j*N+k]=1;
                                ds.cz[i*N*N+j*N+k]=1;
                                if (dist < actuatorDistance[i*N*N+j*N+k]) {
                                    ds.actuator[i*N*N+j*N+k] = nactuator;
                                    actuatorDistance[i*N*N+j*N+k] = dist;
                                }
                                if (dist < actuatorClosestDistance[nactuator]) {
                                    actuatorClosestDistance[nactuator]=dist;
                                    ds.actuatorClosestI[nactuator]=i;
                                    ds.actuatorClosestJ[nactuator]=j;
                                    actuatorClosestK[nactuator]=k;
                                }
                                constraint++;
                            }
                        }
                    }
                }
            }
            if (found==0) {
                printf("Error:  couldn't find constraint.\n");
            } else {
                if (foundT == 0) {
                    ds.cx[ib*N*N+jb*N+kb]=1;
                    ds.cy[ib*N*N+jb*N+kb]=1;
                    ds.cz[ib*N*N+jb*N+kb]=1;
                    if (mindist < actuatorDistance[ib*N*N+jb*N+kb]) {
                        ds.actuator[ib*N*N+jb*N+kb] = nactuator;
                        actuatorDistance[ib*N*N+jb*N+kb] = mindist;
                    }
                    actuatorClosestDistance[nactuator]=mindist;
                    ds.actuatorClosestI[nactuator]=ib;
                    ds.actuatorClosestJ[nactuator]=jb;
                    actuatorClosestK[nactuator]=kb;
                    constraint++;
                }
            }
            ds.actuatorX.push_back(xx);
            ds.actuatorY.push_back(yy);
            ds.actuatorR.push_back(ds.length);
            nactuator++;
         }
    }
    printf("Actuators/Supports:  %ld\n",nactuator);
    printf("Constraints: %ld\n",constraint);

    // temporary
    // thicknessPoint = 0;
    ds.zernikeMode = 1;
    printf("Angular Tolerance:  %lf\n",angTol);
    printf("Initial Step:  %e\n",initStep);
    double angularTolerance = angTol*angTol*(static_cast<double>(N))*(static_cast<double>(N))/64.0/64.0;
    double surfaceTolerance = 1e-9;
    ds.nleg=24*0;
    ds.nphi=24*0;
    ds.nzern=(91);
    ds.npert=ds.nzern+ds.nleg*1+ds.nphi;

    if (ds.zernikeMode == 0) ds.npert=nactuator+ds.nleg+ds.nphi;

    int nhistory=100;
    double chisigma=1e6;
    double dispsigma=1.0;
    double *csigma, *perthist, *chihist, *disphist;
    ds.pert = static_cast<double*>(malloc(ds.npert*sizeof(double)));
    cstep = static_cast<double*>(malloc(ds.npert*sizeof(double)));
    csigma = static_cast<double*>(malloc(ds.npert*sizeof(double)));
    perthist =  static_cast<double*>(malloc(ds.npert*nhistory*sizeof(double)));
    chihist = static_cast<double*>(malloc(nhistory*sizeof(double)));
    disphist = static_cast<double*>(malloc(nhistory*sizeof(double)));
    for (int l=0; l<ds.npert; l++) ds.pert[l]=0.0;
    double minstep = 1e-10;
    int modCounter = 10000;
    for (int l=0; l<ds.npert; l++) ds.pert[l]=0.0;
    for (int l=0; l<ds.npert; l++) cstep[l]=0.0;
    for (int l=0; l<ds.npert; l++) csigma[l]=0.0;
    for (int l=0; l<ds.npert-ds.nleg-ds.nphi; l++) csigma[l]=initStep;
    for (int l=ds.npert-ds.nleg-ds.nphi; l<ds.npert; l++) csigma[l]=initStep;
    for (int l=0; l<ds.npert; l++) for (int k=0; k<nhistory;k++) perthist[k*ds.npert+l]=0.0;
    for (int l=0; l<nhistory; l++) chihist[l]=0.0;
    for (int l=0; l<nhistory; l++) disphist[l]=0.0;


    int histcounter=0;
    int thistcounter=0;
    int *ord;
    double *ordw;
    ds.maxzern=91;
    ds.acte = static_cast<double*>(malloc(nactuator*sizeof(double)));
    ds.wfe = static_cast<double*>(malloc(ds.maxzern*sizeof(double)));
    ord = static_cast<int*>(malloc(ds.maxzern*sizeof(int)));
    ordw = static_cast<double*>(malloc(ds.maxzern*sizeof(double)));
    for (int l=0;l<ds.maxzern;l++) ord[l]=0;
    for (int l=0;l<=0;l++) ord[l]=1;
    for (int l=1;l<=2;l++) ord[l]=2;
    for (int l=3;l<=5;l++) ord[l]=3;
    for (int l=6;l<=9;l++) ord[l]=4;
    for (int l=10;l<=14;l++) ord[l]=5;
    for (int l=15;l<=20;l++) ord[l]=6;
    for (int l=21;l<=27;l++) ord[l]=7;
    for (int l=28;l<=35;l++) ord[l]=8;
    for (int l=36;l<=44;l++) ord[l]=9;
    for (int l=45;l<=54;l++) ord[l]=10;
    for (int l=55;l<=65;l++) ord[l]=11;
    for (int l=66;l<=77;l++) ord[l]=12;
    for (int l=78;l<=90;l++) ord[l]=13;
    if (ds.zernikeMode == 1) for (int l=0; l<ds.npert-ds.nleg-ds.nphi; l++) csigma[l]=csigma[l]/(ord[l])/(ord[l]);
    if (ds.zernikeMode == 1) csigma[0]=0.0;
    long rej=0, acc=0;
    double firstChi=0.0, secondChi=0.0, bestChi=0.0;
    ds.firstDisp=0.0;
    ds.secondDisp=0.0;
    double bestDisp=0.0;
    int firstTime = 0;
    int secondTime =0;
    int newFirst = 0;
    ds.controlState =0;
    int pause = 1;
    double tact=0.0, twfe=0.0;
    for (int l=0;l<ds.maxzern;l++) ordw[l]=0.0;
    // normalized gamma(n+1-11/6)/gamma(n+2+11/6)*(n+1)
    // excluding first two terms
    ordw[1]=0.0;
    ordw[2]=0.0;
    ordw[3]=0.6667/3.0;
    ordw[4]=0.1780/4.0;
    ordw[5]=0.0705/5.0;
    ordw[6]=0.0342/6.0;
    ordw[7]=0.0188/7.0;
    ordw[8]=0.0113/8.0;
    ordw[9]=0.0072/9.0;
    ordw[10]=0.0049/10.0;
    ordw[11]=0.0034/11.0;
    ordw[12]=0.0024/12.0;
    ordw[13]=0.0018/13.0;
    for (int l=0;l<nactuator;l++) {
        ds.acte[l] = 0.0;
    }
    for (int l=0;l<ds.maxzern;l++) {
        ds.wfe[l] = 0.0;
    }
    for (int l=0;l<ds.maxzern;l++) {
        if ((control == 1) || (control == 3)) {
            ds.wfe[l] = random.normal()*wfeError*ordw[ord[l]]*1e-9;
        } else {
            ds.wfe[l] = 0.0;
        }
        twfe+= ds.wfe[l]*ds.wfe[l];
    }
    for (int l=0;l<nactuator;l++) {
        ds.acte[l] = random.normal()*actError*1e-9;
        tact+= ds.acte[l]*ds.acte[l];
    }
    printf("Actuator Error: %lf %e %e\n",sqrt(tact/static_cast<double>(nactuator))*1e9,tact,static_cast<double>(nactuator));
    printf("Wavefront Estimate Error: %lf\n",sqrt(twfe)*1e9);


    ds.dzmean = 0.0;
    ds.dzcount = 0.0;
    for (long i = ii2; i <= ii1; i++) {
        for (long j = jj2; j <= jj1; j++) {
            for (long k = kk2; k <= kk1; k++) {
                long index = i*N*N + j*N + k;
                double tlocal = (ds.dtemp + ds.xoo[index]*ds.dtdx*sin(ds.ttheta)*cos(ds.tphi) +
                                 ds.yoo[index]*ds.dtdx*sin(ds.ttheta)*sin(ds.tphi) +
                                 (ds.zoo[index]-ds.averageHeight)*ds.dtdx*cos(ds.ttheta)) - ds.temperature0;
                ds.xo[index]=ds.xoo[index]*(1+ds.alpha*tlocal);
                ds.yo[index]=ds.yoo[index]*(1+ds.alpha*tlocal);
                ds.zo[index]=(ds.zoo[index]-ds.averageHeight)*(1+ds.alpha*tlocal)+ds.averageHeight;
                if (ds.surface[index]==1) {
                    ds.dzmean+=(ds.zo[index]-ds.zoo[index]);
                    ds.dzcount+=1.0;
                }
            }
        }
    }
    ds.dzmean = ds.dzmean/ds.dzcount;

    long counter = 0;
    double vpeak = 1e-10;
    ds.vmag = 1.0;
    ds.final = 0;
    int skipresimflag=0;

    while (ds.vmag > tol*surfaceTolerance/ds.dt || ds.vmag > tol2*tol*sqrt((dispsigma+1e-14)/(static_cast<double>(surfacePoints)))/ds.dt ||
           (chisigma >= sqrt(2.0*static_cast<double>(ds.npert))*angularTolerance) || (ds.final==0)) {


        // move constraint
        if (ds.vmag > tol*surfaceTolerance/ds.dt || ds.vmag > tol2*tol*sqrt((dispsigma+1e-14)/(static_cast<double>(surfacePoints)))/ds.dt) {
            pause=1;
            if ((control == 1) || (control == 3)) {
                if (ds.controlState == 1) {
                    firstChi = 0.0;
                    ds.firstDisp = 0.0;
                } else {
                    secondChi = 0.0;
                    ds.secondDisp = 0.0;
                }
            } else {
                chisigma=0.0;
                bestChi = secondChi;
                bestDisp = ds.secondDisp;
                firstChi = 0;
                ds.firstDisp = 0.0;
                secondChi = 0;
                ds.secondDisp = 0.0;
            }
        } else {
            pause=0;
            if ((control == 1) || (control==3)) {
                if (chisigma < sqrt(2.0*static_cast<double>(ds.npert))*angularTolerance) {
                    if (ds.final == 0) {
                        // for (int l = 0; l < nactuator; l++) ds.pert[l]+=ds.acte[l];
                    }
                    ds.final=1;
                    chisigma=0.0;
                    ds.controlState = 0;
                } else {
                skipresim:;
                    if (ds.controlState == 0) {
                        if (firstTime != 0) {
                            double burn = static_cast<double>(thistcounter)/(1.0*static_cast<double>(nhistory));
                            if (burn > 1.0) burn=1.0;
                            if (burn < 0.01) burn=0.01;
                            // printf("%d %lf %lf\n",newFirst,sqrt(firstChi/totalDerivativePoints),sqrt(secondChi/totalDerivativePoints));
                            if (newFirst==0) {
                                if ((secondChi <= firstChi) ||
                                    ((secondChi > firstChi) && (random.uniform() < exp(-burn*fabs(firstChi-secondChi)/2.0/angularTolerance)))) {
                                    bestChi = secondChi;
                                    bestDisp = ds.secondDisp;
                                    acc++;
                                } else {
                                    bestChi = firstChi;
                                    bestDisp = ds.firstDisp;
                                    for (int l=0; l<ds.npert; l++) ds.pert[l] -= cstep[l];
                                    rej++;
                                }
                            } else {
                                if ((firstChi <= secondChi) ||
                                    ((firstChi > secondChi) && (random.uniform() < exp(-burn*fabs(firstChi-secondChi)/2.0/angularTolerance)))) {

                                    for (int l=0; l<ds.npert; l++) ds.pert[l] += cstep[l];
                                    bestChi = firstChi;
                                    bestDisp = ds.firstDisp;
                                    acc++;
                                } else {
                                    bestChi = secondChi;
                                    bestDisp = ds.secondDisp;
                                    rej++;
                                }
                            }
                        // reconfigure steps
                            for (int l=0; l<ds.npert; l++) perthist[histcounter*ds.npert+l]=ds.pert[l];
                            chihist[histcounter]=bestChi;
                            disphist[histcounter]=bestDisp;
                            histcounter++;
                            thistcounter++;
                            if (histcounter>=nhistory) histcounter=0;
                            if (thistcounter >= nhistory) {
                                for (int l=0; l<ds.npert; l++) {
                                    double mean=0.0,dev=0.0;
                                    for (int k=0; k<nhistory; k++) {
                                        mean+=perthist[k*ds.npert+l];
                                        dev+=perthist[k*ds.npert+l]*perthist[k*ds.npert+l];
                                    }
                                    double num=static_cast<double>(nhistory);
                                    csigma[l] = sqrt(1.0/num*(dev-mean*mean/num));
                                    if (csigma[l] < minstep) csigma[l]=minstep;
                                }
                                double mean=0.0,dev=0.0;
                                for (int k=0; k<nhistory; k++) {
                                    mean+=chihist[k];
                                    dev+=chihist[k]*chihist[k];
                                }
                                double num=static_cast<double>(nhistory);
                                chisigma=sqrt(1.0/num*(dev-mean*mean/num));
                                mean=0.0,dev=0.0;
                                for (int k=0; k<nhistory; k++) {
                                    mean+=disphist[k];
                                    dev+=disphist[k]*disphist[k];
                                }
                                num=static_cast<double>(nhistory);
                                if ((dev-mean*mean/num) > 0) dispsigma=sqrt(1.0/num*(dev-mean*mean/num)); else dispsigma=0.0;
                            } else {
                                chisigma=1e6;
                                dispsigma=1;
                                if (totalDerivativePoints>0) chisigma=chisigma*totalDerivativePoints;
                                if (surfacePoints>0) dispsigma=dispsigma*surfacePoints;
                            }
                        }
                        if (random.uniform() < 0.5) newFirst=1; else newFirst=0;
                        if (firstTime == 0) newFirst=1;
                        // initialize
                        firstChi = 0;
                        ds.firstDisp = 0.0;
                        secondChi = 0;
                        ds.secondDisp = 0.0;
                        if (newFirst == 1) {
                            for (int l=0; l < ds.npert; l++) cstep[l] = random.normal()*csigma[l]*2.38/sqrt(static_cast<double>(ds.npert));
                            for (int l=0; l < ds.npert; l++) ds.pert[l] += cstep[l];
                        }
                        for (int l=0;l<nactuator;l++) {
                            if ((control == 1) || (control==3)) {
                                ds.acte[l] = random.normal()*actError*1e-9;
                            } else {
                                ds.acte[l] = 0.0;
                            }
                        }
                    }
                    if (ds.controlState == 1) {
                        if (newFirst == 0) {
                            for (int l=0; l < ds.npert; l++) cstep[l] = random.normal()*csigma[l]*2.38/sqrt(static_cast<double>(ds.npert));
                            for (int l=0; l < ds.npert; l++) ds.pert[l] += cstep[l];
                        } else {
                            for (int l=0; l < ds.npert; l++) ds.pert[l] -= cstep[l];
                        }
                    }
                    ds.controlState++;
                    if (ds.controlState == 2) ds.controlState = 0;
                    //
                    // if controlstate==1 && newfirst==1 will be trying out new (need to continue)
                    // if controlstate==0 && newfirst==1 will be keeping old so secondChi=bestChi & can skip
                    // if controlstate==1 && newfirst==0 will be keeping old so firstChi=bestChi & can skip
                    // if controlstate==0 && newfirst==0 now trying out new (need to continue)
                    skipresimflag=0;
                    if ((ds.controlState==0) && (newFirst==1) && (secondTime==1)) {
                        secondChi=bestChi;
                        ds.secondDisp=bestDisp;
                        skipresimflag=1;
                    }
                    if ((ds.controlState==1) && (newFirst==0) && (secondTime==1)) {
                        firstChi=bestChi;
                        ds.firstDisp=bestDisp;
                        skipresimflag=1;
                    }
                    firstTime = 1;
                    if (ds.controlState == 0) secondTime=1;

                }
            } else {
                ds.final = 1;
            }
        }

        if ((((control==1) || (control==3)) && (((pause == 0) && (ds.controlState == 0)))) ||
            (((control==0) || (control==2)) && ((counter % 100)==0)) ||
            (((control==1) || (control==3)) && ((counter % 1000)==0))){
            printf("%6d %8ld %8.4lf%% %4.0lf%% %8.4lf arcsec %8.4lf arcsec %8.4lf nm ", thistcounter,counter, ds.vmag/(surfaceTolerance/ds.dt)*100, static_cast<double>(rej)/static_cast<double>(rej+acc)*100,sqrt(bestChi/totalDerivativePoints),sqrt(chisigma/totalDerivativePoints),sqrt((dispsigma+1e-14)/(static_cast<double>(surfacePoints)))*1e9);
            for (int l=0; l<ds.npert; l++) if (fabs(ds.pert[l]*1e6) > 0.01) printf("%d: %4.0lf   ",l,ds.pert[l]*1e9);
            printf("\n");
        }
        if (skipresimflag==1) {
            skipresimflag=0;
            goto skipresim;
        }


        long nn=0;
        for (long i=ii2; i<=ii1; i++) {
            for (long j=jj2; j<=jj1; j++) {
                for (long k=kk2; k<=kk1; k++) {
                    long index=i*N*N+j*N+k;
                    done[index]=0;
                    if (ds.e[index]==1) {
                        ic[nn]=i;
                        jc[nn]=j;
                        kc[nn]=k;
                        indc[nn]=index;
                        nc[nn]=nn;
                        nn++;
                    }
                }
            }
        }
        long nnorig=nn;

        ds.vc = 0.0;
        ds.vmag = 0.0;
        for (long iii = 0; iii < nnorig; iii=iii+nodeperthread) {

            long ln2 = nodeperthread;
            if (ln2 > (nnorig - iii)) {
                ln2 = (nnorig - iii);
            }

            int bestthread=0;
            for (int i = 0; i < numthread; i++) {
                if (openthread[i]==0) bestthread=i;
            }
            pthread_mutex_lock(&lock1);
            openthread[bestthread]=1;
            openthreads++;
            pthread_mutex_unlock(&lock1);
            for (long kkk = 0; kkk < ln2; kkk++) {


                // memmove(&nc[vv],&nc[vv+1],(nn-(vv+1))*sizeof(long));
                // nn--;

                long vv=0;
                do {
                    vv = floor(random.uniform()*nn);
                } while (done[indc[nc[vv]]] == 1);

                done[indc[nc[vv]]] = 1;
                args[bestthread].i[kkk]=ic[nc[vv]];
                args[bestthread].j[kkk]=jc[nc[vv]];
                args[bestthread].k[kkk]=kc[nc[vv]];
            }
            args[bestthread].N=N;
            args[bestthread].cc=ln2;
            args[bestthread].nt=bestthread;
            pthread_create(&thread[bestthread],NULL,elastic,&args[bestthread]);
            pthread_detach(thread[bestthread]);
            while (openthreads >= numthread) {
                usleep(10);
            }

        }

        while (openthreads != 0) {
            usleep(10);
        }

        ds.vmag = ds.vmag/ds.vc;
        if (ds.vmag > vpeak) vpeak = ds.vmag;
        if (vpeak < 1e-10) vpeak=1e-10;

        counter++;

        //derivative of surface
        if (ds.vmag > tol*surfaceTolerance/ds.dt || ds.vmag > tol2*tol*sqrt((dispsigma+1e-14)/(static_cast<double>(surfacePoints)))/ds.dt) {

        } else {

        avgr0=0.0;
        avgr0n=0.0;
        avgr1=0.0;
        avgr1n=0.0;
        for (long i = ii2; i <= ii1; i++) {
            for (long j = jj2; j <= jj1; j++) {
                if (ds.surfaceErrorPoint[i*N+j]==1) {
                    goodDerivative[i*N+j]=1;
                    long steppp=1;
                    long ilow = i-steppp;
                    if (ilow < ii2) {
                        ilow=ii2;
                        goodDerivative[i*N+j]=0;
                    }
                    long ihigh = i+steppp;
                    if (ihigh >= ii1) {
                        ihigh=ii1;
                        goodDerivative[i*N+j]=0;
                    }
                    long jlow = j-steppp;
                    if (jlow < jj2) {
                        jlow=jj2;
                        goodDerivative[i*N+j]=0;
                    }
                    long jhigh = j+steppp;
                    if (jhigh >= jj1) {
                        jhigh=jj1;
                        goodDerivative[i*N+j]=0;
                    }

                long s11 = i*N + j;
                long s12 = i*N + jhigh;
                long s10 = i*N + jlow;
                long s01 = ilow*N + j;
                long s02 = ilow*N + jhigh;
                long s00 = ilow*N + jlow;
                long s21 = ihigh*N + j;
                long s22 = ihigh*N + jhigh;
                long s20 = ihigh*N + jlow;

                double factor=1.0/(2.0*ds.length/(static_cast<double>(N-1)))/ARCSEC/static_cast<double>(steppp);

                surfaceDerivativeY[s11] = ((ds.surfaceError[s12]-ds.surfaceError[s10])*4.0 +
                                           (ds.surfaceError[s02]-ds.surfaceError[s00])*1.0 +
                                           (ds.surfaceError[s22]-ds.surfaceError[s20])*1.0)/2.0/6.0*factor;
                surfaceDerivativeX[s11] = ((ds.surfaceError[s21]-ds.surfaceError[s01])*4.0 +
                                           (ds.surfaceError[s20]-ds.surfaceError[s00])*1.0 +
                                           (ds.surfaceError[s22]-ds.surfaceError[s02])*1.0)/2.0/6.0*factor;

                if (ds.surfaceErrorPoint[s11]==0) goodDerivative[i*N+j]=0;
                if (ds.surfaceErrorPoint[s12]==0) goodDerivative[i*N+j]=0;
                if (ds.surfaceErrorPoint[s10]==0) goodDerivative[i*N+j]=0;
                if (ds.surfaceErrorPoint[s01]==0) goodDerivative[i*N+j]=0;
                if (ds.surfaceErrorPoint[s02]==0) goodDerivative[i*N+j]=0;
                if (ds.surfaceErrorPoint[s00]==0) goodDerivative[i*N+j]=0;
                if (ds.surfaceErrorPoint[s21]==0) goodDerivative[i*N+j]=0;
                if (ds.surfaceErrorPoint[s22]==0) goodDerivative[i*N+j]=0;
                if (ds.surfaceErrorPoint[s20]==0) goodDerivative[i*N+j]=0;
                if (fused==1) {
                    double x0,y0,r0;
                    int part=0;
                    x0 = ds.xoo[s11*N + ds.surfaceErrorZ[s11]] * 1e3;
                    y0 = ds.yoo[s11*N + ds.surfaceErrorZ[s11]] * 1e3;
                    r0= sqrt(x0*x0 + y0*y0);
                    if (((r0 < innerradius[0]) || (r0 > outerradius[0])) &&
                        ((r0 < innerradius[1]) || (r0 > outerradius[1]))) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[0]) && (r0 <= outerradius[0])) part=0;
                    if ((r0 >= innerradius[1]) && (r0 <= outerradius[1])) part=1;
                    x0 = ds.xoo[s12*N + ds.surfaceErrorZ[s12]] * 1e3;
                    y0 = ds.yoo[s12*N + ds.surfaceErrorZ[s12]] * 1e3;
                    r0= sqrt(x0*x0 + y0*y0);
                    if (((r0 < innerradius[0]) || (r0 > outerradius[0])) &&
                        ((r0 < innerradius[1]) || (r0 > outerradius[1]))) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[0]) && (r0 <= outerradius[0]) && (part==1)) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[1]) && (r0 <= outerradius[1]) && (part==0)) goodDerivative[s11]=0;
                    x0 = ds.xoo[s10*N + ds.surfaceErrorZ[s10]] * 1e3;
                    y0 = ds.yoo[s10*N + ds.surfaceErrorZ[s10]] * 1e3;
                    r0= sqrt(x0*x0 + y0*y0);
                    if (((r0 < innerradius[0]) || (r0 > outerradius[0])) &&
                        ((r0 < innerradius[1]) || (r0 > outerradius[1]))) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[0]) && (r0 <= outerradius[0]) && (part==1)) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[1]) && (r0 <= outerradius[1]) && (part==0)) goodDerivative[s11]=0;
                    x0 = ds.xoo[s21*N + ds.surfaceErrorZ[s21]] * 1e3;
                    y0 = ds.yoo[s21*N + ds.surfaceErrorZ[s21]] * 1e3;
                    r0= sqrt(x0*x0 + y0*y0);
                    if (((r0 < innerradius[0]) || (r0 > outerradius[0])) &&
                        ((r0 < innerradius[1]) || (r0 > outerradius[1]))) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[0]) && (r0 <= outerradius[0]) && (part==1)) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[1]) && (r0 <= outerradius[1]) && (part==0)) goodDerivative[s11]=0;
                    x0 = ds.xoo[s22*N + ds.surfaceErrorZ[s22]] * 1e3;
                    y0 = ds.yoo[s22*N + ds.surfaceErrorZ[s22]] * 1e3;
                    r0= sqrt(x0*x0 + y0*y0);
                    if (((r0 < innerradius[0]) || (r0 > outerradius[0])) &&
                        ((r0 < innerradius[1]) || (r0 > outerradius[1]))) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[0]) && (r0 <= outerradius[0]) && (part==1)) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[1]) && (r0 <= outerradius[1]) && (part==0)) goodDerivative[s11]=0;
                    x0 = ds.xoo[s20*N + ds.surfaceErrorZ[s20]] * 1e3;
                    y0 = ds.yoo[s20*N + ds.surfaceErrorZ[s20]] * 1e3;
                    r0= sqrt(x0*x0 + y0*y0);
                    if (((r0 < innerradius[0]) || (r0 > outerradius[0])) &&
                        ((r0 < innerradius[1]) || (r0 > outerradius[1]))) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[0]) && (r0 <= outerradius[0]) && (part==1)) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[1]) && (r0 <= outerradius[1]) && (part==0)) goodDerivative[s11]=0;
                    x0 = ds.xoo[s01*N + ds.surfaceErrorZ[s01]] * 1e3;
                    y0 = ds.yoo[s01*N + ds.surfaceErrorZ[s01]] * 1e3;
                    r0= sqrt(x0*x0 + y0*y0);
                    if (((r0 < innerradius[0]) || (r0 > outerradius[0])) &&
                        ((r0 < innerradius[1]) || (r0 > outerradius[1]))) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[0]) && (r0 <= outerradius[0]) && (part==1)) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[1]) && (r0 <= outerradius[1]) && (part==0)) goodDerivative[s11]=0;
                    x0 = ds.xoo[s02*N + ds.surfaceErrorZ[s02]] * 1e3;
                    y0 = ds.yoo[s02*N + ds.surfaceErrorZ[s02]] * 1e3;
                    r0= sqrt(x0*x0 + y0*y0);
                    if (((r0 < innerradius[0]) || (r0 > outerradius[0])) &&
                        ((r0 < innerradius[1]) || (r0 > outerradius[1]))) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[0]) && (r0 <= outerradius[0]) && (part==1)) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[1]) && (r0 <= outerradius[1]) && (part==0)) goodDerivative[s11]=0;
                    x0 = ds.xoo[s00*N + ds.surfaceErrorZ[s00]] * 1e3;
                    y0 = ds.yoo[s00*N + ds.surfaceErrorZ[s00]] * 1e3;
                    r0= sqrt(x0*x0 + y0*y0);
                    if (((r0 < innerradius[0]) || (r0 > outerradius[0])) &&
                        ((r0 < innerradius[1]) || (r0 > outerradius[1]))) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[0]) && (r0 <= outerradius[0]) && (part==1)) goodDerivative[s11]=0;
                    if ((r0 >= innerradius[1]) && (r0 <= outerradius[1]) && (part==0)) goodDerivative[s11]=0;
                }
                if (goodDerivative[s11]==1) {
                    double x0 = ds.xoo[s11*N + ds.surfaceErrorZ[s11]] * 1e3;
                    double y0 = ds.yoo[s11*N + ds.surfaceErrorZ[s11]] * 1e3;
                    double r0= sqrt(x0*x0 + y0*y0);
                    if (fused==1) {
                        if ((r0 > innerradius[0]) && (r0 < outerradius[0])) {
                            avgr0+=r0;
                            avgr0n+=1.0;
                        }
                        if ((r0 > innerradius[1]) && (r0 < outerradius[1])) {
                            avgr1+=r0;
                            avgr1n+=1.0;
                        }
                    } else {
                        avgr0+=r0;
                        avgr0n+=1.0;
                    }
                }
                }
            }
        }
        avgr0 /= avgr0n;
        avgr1 /= avgr1n;
        totalDerivativePoints = 0.0;
        meanRDerivative=0.0;
        meanRDerivative1=0.0;
        double totalRDerivativePoints = 0.0;
        double totalRDerivativePoints1 = 0.0;
        for (long i = ii2; i <= ii1; i++) {
            for (long j = jj2; j <= jj1; j++) {
                if (goodDerivative[i*N+j]==1) {
                    long s11 = i*N + j;
                    double x0 = ds.xoo[i*N*N+j*N + ds.surfaceErrorZ[s11]] * 1e3;
                    double y0 = ds.yoo[i*N*N+j*N + ds.surfaceErrorZ[s11]] * 1e3;
                    double r0= sqrt(x0*x0 + y0*y0);
                    if (fused == 1) {
                        if ((r0 <= outerradius[0]) && (r0 >= innerradius[0])) {
                            meanRDerivative += surfaceDerivativeX[s11]*x0/r0 + surfaceDerivativeY[s11]*y0/r0;
                            totalRDerivativePoints+=1.0;
                        }
                        if ((r0 <= outerradius[1]) && (r0 >= innerradius[1])) {
                            meanRDerivative1 += surfaceDerivativeX[s11]*x0/r0 + surfaceDerivativeY[s11]*y0/r0;
                            totalRDerivativePoints1+=1.0;
                        }
                    } else {
                        meanRDerivative += surfaceDerivativeX[s11]*x0/r0 + surfaceDerivativeY[s11]*y0/r0;
                        totalRDerivativePoints+=1.0;
                    }
                }
            }
        }
        totalDerivativePoints = totalRDerivativePoints + totalRDerivativePoints1;
        meanRDerivative /= totalRDerivativePoints;
        meanRDerivative1 /= totalRDerivativePoints1;
        meanRDerivative = -1.0/(-radiusofcurv[0])*ds.alpha*(ds.dtemp-ds.temperature0)/ARCSEC;
        meanRDerivative1 = -1.0/(-radiusofcurv[1])*ds.alpha*(ds.dtemp-ds.temperature0)/ARCSEC;
        // could add correction (1-(1+kappa)*r^2/r0^2)^(-3/2.) for non-parabolic but not as simple as change of radius of curvature for these
        double totalDerivative = 0.0;
        for (long i = ii2; i <= ii1; i++) {
            for (long j = jj2; j <= jj1; j++) {
                if (goodDerivative[i*N+j]==1) {
                    long s11 = i*N + j;
                    double x0 = ds.xoo[i*N*N+j*N + ds.surfaceErrorZ[s11]] * 1e3;
                    double y0 = ds.yoo[i*N*N+j*N + ds.surfaceErrorZ[s11]] * 1e3;
                    double r0= sqrt(x0*x0 + y0*y0);
                    if (fused == 1) {
                        if ((r0 <= outerradius[0]) && (r0 >= innerradius[0])) {
                            totalDerivative += pow(surfaceDerivativeX[s11]-meanRDerivative*x0,2.0) +
                                pow(surfaceDerivativeY[s11]-meanRDerivative*y0,2.0);
                        }
                        if ((r0 <= outerradius[1]) && (r0 >= innerradius[1])) {
                            totalDerivative += pow(surfaceDerivativeX[s11]-meanRDerivative1*x0,2.0) +
                                pow(surfaceDerivativeY[s11]-meanRDerivative1*y0,2.0);
                       }
                    } else {
                        totalDerivative += pow(surfaceDerivativeX[s11]-meanRDerivative*x0,2.0) +
                                pow(surfaceDerivativeY[s11]-meanRDerivative*y0,2.0);

                    }
                }
            }
        }



        if (ds.controlState == 1) {
            firstChi = totalDerivative;
        } else {
            secondChi = totalDerivative;
        }
        }



        if (((((counter-1) % modCounter) == 0) || (modCounter==1)) && (pertDebug==1)) {

            std::ostringstream surface1, surface2;
            surface1 << "feac_"  << obsID << "_" << surfaceIndex << "_" << floor((counter-1)/static_cast<double>(modCounter)) << ".txt";
            surface2 << "feac_"  << obsID << "_" << secondSurfaceIndex << "_" << floor((counter-1)/static_cast<double>(modCounter)) << ".txt";

            std::ofstream output1(surface1.str().c_str());
            std::ofstream output2(surface2.str().c_str());
            top=0;
            bottom=0;
            for (long i = ii2; i <= ii1; i++) {
                for (long j = jj2; j <= jj1; j++) {
                    for (long k = kk2; k <= kk1; k++) {
                        long idx = i*N*N + j*N + k;
                        if (ds.e[idx] == 1) {
                            top = 0;
                            bottom = 0;
                            if (k + 1 < N) {
                                if (ds.e[idx + 1] == 0) top = 1;
                            } else top = 1;
                            if (k > 0) {
                                if (ds.e[idx - 1] == 0) bottom = 1;
                            } else bottom = 1;
                            if (top == 1 && bottom == 1) {
                                if (warn==0) std::cout << "Warning single layer" << std::endl;
                                warn++;
                            }
                            if (top == 1 || bottom == 1) {
                                double x0 = ds.xoo[idx] * 1e3; // in mm
                                double y0 = ds.yoo[idx] * 1e3;
                                double rr = sqrt(pow(x0, 2) + pow(y0, 2));
                                double dx = (ds.x[idx] - ds.xoo[idx]) * 1e3;
                                double dy = (ds.y[idx] - ds.yoo[idx]) * 1e3;
                                double dz = (ds.z[idx] - ds.zoo[idx]) * 1e3;
                                if (single == 1) {
                                    if ((ds.mountType == 0) && (top == 1)) {
                                        double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                output1 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 << " "
                                        << dx << " " << dy << " " << dz << " 0 0 0\n";
                            } else if ((ds.mountType == 1) && (bottom == 1)) {
                                double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                output1 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 << " "
                                        << dx << " " << dy << " " << dz << " 0 0 0\n";
                            }
                        }
                        if (doub == 1) {
                            if (bottom == 1) {
                                double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                output1 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 << " "
                                        << dx << " " << dy << " " << dz << " 0 0 0\n";
                            }
                            if (top == 1) {
                                double z0 = asphere(rr, radiusofcurv[1], height[1], conic[1], third[1], fourth[1], fifth[1],
                                                    sixth[1], seventh[1], eighth[1], ninth[1], tenth[1]) * 1e3 + height0;
                                output2 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 << " "
                                        << dx << " " << dy << " " << dz << " 0 0 0\n";
                            }

                        }
                        if (fused == 1) {
                            if ((ds.mountType == 0) && (top == 1)) {
                                double rr = sqrt(pow(x0, 2) + pow(y0, 2));
                                if ((rr <= outerradius[0]) && (rr >= innerradius[0])) {
                                double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                output1 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 << " "
                                        << dx << " " << dy << " " << dz << " 0 0 0\n";
                                }
                                if ((rr <= outerradius[1]) && (rr >= innerradius[1])) {
                                double z0 = asphere(rr, radiusofcurv[1], height[1], conic[1], third[1], fourth[1], fifth[1],
                                                    sixth[1], seventh[1], eighth[1], ninth[1], tenth[1]) * 1e3 + height0;
                                output2 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 << " "
                                        << dx << " " << dy << " " << dz << " 0 0 0\n";
                                }
                            }
                            if ((ds.mountType == 1) && (bottom == 1)) {
                                double rr = sqrt(pow(x0, 2) + pow(y0, 2));
                                if ((rr <= outerradius[0]) && (rr >= innerradius[0])) {
                                double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                output1 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 << " "
                                        << dx << " " << dy << " " << dz << " 0 0 0\n";
                                }
                                if ((rr <= outerradius[1]) && (rr >= innerradius[1])) {
                                double z0 = asphere(rr, radiusofcurv[1], height[1], conic[1], third[1], fourth[1], fifth[1],
                                                    sixth[1], seventh[1], eighth[1], ninth[1], tenth[1]) * 1e3 + height0;
                                output2 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 << " "
                                        << dx << " " << dy << " " << dz << " 0 0 0\n";
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    output1.close();
    output2.close();

        }
    }


    int outFlag1=1;
    int outFlag2=1;
    std::ostringstream surface1, surface2;
    surface1 << "fea_"  << obsID << "_" << surfaceIndex << ".txt";
    surface2 << "fea_"  << obsID << "_" << secondSurfaceIndex << ".txt";

    std::ofstream output1(surface1.str().c_str());
    std::ofstream output2(surface2.str().c_str());
    top=0;
    bottom=0;
    for (long i = ii2; i <= ii1; i++) {
        for (long j = jj2; j <= jj1; j++) {
            for (long k = kk2; k <= kk1; k++) {
                long idx = i*N*N + j*N + k;
                if (ds.e[idx] == 1) {
                    top = 0;
                    bottom = 0;
                    if (k + 1 < N) {
                        if (ds.e[idx + 1] == 0) top = 1;
                    } else top = 1;
                    if (k > 0) {
                        if (ds.e[idx - 1] == 0) bottom = 1;
                    } else bottom = 1;
                    if (top == 1 && bottom == 1) {
                        if (warn==0) std::cout << "Warning single layer" << std::endl;
                        warn++;
                    }
                    if (top == 1 || bottom == 1) {

                        double x0 = ds.xoo[idx] * 1e3; // in mm
                        double y0 = ds.yoo[idx] * 1e3;
                        double rr = sqrt(pow(x0, 2) + pow(y0, 2));
                        double dx = (ds.x[idx] - ds.xoo[idx]) * 1e3;
                        double dy = (ds.y[idx] - ds.yoo[idx]) * 1e3;
                        double dz = (ds.z[idx] - ds.zoo[idx]) * 1e3;
                        if (single == 1) {
                            if ((ds.mountType == 0) && (top == 1)) {
                                double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                double defocus = ds.alpha*(ds.dtemp-ds.temperature0)*radiusofcurv[0]/2.0*focus[0];
                                if ((control==2) || (control==3)) defocus=0.0;
                                if (outFlag1) output1 << std::scientific << std::setprecision(17)
                                         << x0 << " " << y0 << " " << z0 + defocus << " "
                                                      << dx << " " << dy << " " << dz << " 0 0 0\n";
                                    else output1 << std::scientific << std::setprecision(17)
                                         << x0 << " " << y0 << " " << z0 << " "
                                                 <<  " 0 0 0 0 0 0\n";
                            } else if ((ds.mountType == 1) && (bottom == 1)) {
                                double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                double defocus = ds.alpha*(ds.dtemp-ds.temperature0)*radiusofcurv[0]/2.0*focus[0];
                                if ((control==2) || (control==3)) defocus=0.0;
                                if (outFlag1) output1 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 + defocus << " "
                                                      << dx << " " << dy << " " << dz << " 0 0 0\n";
                                    else output1 << std::scientific << std::setprecision(17)
                                         << x0 << " " << y0 << " " << z0 << " "
                                                 <<  " 0 0 0 0 0 0\n";
                            }
                        }
                        if (doub == 1) {
                            if (bottom == 1) {
                                double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                double defocus = 0.0;
                                if ((control==2) || (control==3)) defocus=0.0;
                                if (outFlag1) output1 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 + defocus << " "
                                                      << dx << " " << dy << " " << dz << " 0 0 0\n";
                                    else output1 << std::scientific << std::setprecision(17)
                                         << x0 << " " << y0 << " " << z0 << " "
                                                 <<  " 0 0 0 0 0 0\n";
                            }
                            if (top == 1) {
                                double z0 = asphere(rr, radiusofcurv[1], height[1], conic[1], third[1], fourth[1], fifth[1],
                                                    sixth[1], seventh[1], eighth[1], ninth[1], tenth[1]) * 1e3 + height0;
                                double defocus = 0.0;
                                if ((control==2) || (control==3)) defocus=0.0;
                                 if (outFlag2) output2 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 + defocus << " "
                                                      << dx << " " << dy << " " << dz << " 0 0 0\n";
                                    else output2 << std::scientific << std::setprecision(17)
                                         << x0 << " " << y0 << " " << z0 << " "
                                                 <<  " 0 0 0 0 0 0\n";
                            }

                        }
                        if (fused == 1) {
                            if ((ds.mountType == 0) && (top == 1)) {
                                double rr = sqrt(pow(x0, 2) + pow(y0, 2));
                                if ((rr <= outerradius[0]) && (rr >= innerradius[0])) {
                                double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                double defocus = ds.alpha*(ds.dtemp-ds.temperature0)*radiusofcurv[0]/2.0*focus[0];
                                if ((control==2) || (control==3)) defocus=0.0;
                                 if (outFlag1) output1 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 + defocus << " "
                                        << dx << " " << dy << " " << dz << " 0 0 0\n";
                                else output1 << std::scientific << std::setprecision(17)
                                         << x0 << " " << y0 << " " << z0 << " "
                                             <<  " 0 0 0 0 0 0\n";
                                }
                                if ((rr <= outerradius[1]) && (rr >= innerradius[1])) {
                                double z0 = asphere(rr, radiusofcurv[1], height[1], conic[1], third[1], fourth[1], fifth[1],
                                                    sixth[1], seventh[1], eighth[1], ninth[1], tenth[1]) * 1e3 + height0;
                                double defocus = ds.alpha*(ds.dtemp-ds.temperature0)*radiusofcurv[1]/2.0*focus[1];
                                if ((control==2) || (control==3)) defocus=0.0;
                                if (outFlag2) output2 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 + defocus << " "
                                                      << dx << " " << dy << " " << dz << " 0 0 0\n";
                                else output2 << std::scientific << std::setprecision(17)
                                         << x0 << " " << y0 << " " << z0 << " "
                                             <<  " 0 0 0 0 0 0\n";
                                }
                            }
                            if ((ds.mountType == 1) && (bottom == 1)) {
                                double rr = sqrt(pow(x0, 2) + pow(y0, 2));
                                if ((rr <= outerradius[0]) && (rr >= innerradius[0])) {
                                double z0 = asphere(rr, radiusofcurv[0], height[0], conic[0], third[0], fourth[0], fifth[0],
                                                    sixth[0], seventh[0], eighth[0], ninth[0], tenth[0]) * 1e3 + height0;
                                double defocus = ds.alpha*(ds.dtemp-ds.temperature0)*radiusofcurv[0]/2.0*focus[0];
                                if ((control==2) || (control==3)) defocus=0.0;
                                if (outFlag1) output1 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 + defocus << " "
                                                      << dx << " " << dy << " " << dz << " 0 0 0\n";
                                    else output1 << std::scientific << std::setprecision(17)
                                         << x0 << " " << y0 << " " << z0 << " "
                                                 <<  " 0 0 0 0 0 0\n";
                                }
                                if ((rr <= outerradius[1]) && (rr >= innerradius[1])) {
                                double z0 = asphere(rr, radiusofcurv[1], height[1], conic[1], third[1], fourth[1], fifth[1],
                                                    sixth[1], seventh[1], eighth[1], ninth[1], tenth[1]) * 1e3 + height0;
                                double defocus = ds.alpha*(ds.dtemp-ds.temperature0)*radiusofcurv[1]/2.0*focus[1];
                                if ((control==2) || (control==3)) defocus=0.0;
                                if (outFlag2) output2 << std::scientific << std::setprecision(17)
                                        << x0 << " " << y0 << " " << z0 + defocus << " "
                                                      << dx << " " << dy << " " << dz << " 0 0 0\n";
                                else output2 << std::scientific << std::setprecision(17)
                                         << x0 << " " << y0 << " " << z0 << " "
                                             <<  " 0 0 0 0 0 0\n";
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    output1.close();
    output2.close();

    // for (int l=0; l<ds.npert; l++) printf("%d: %6.1lf %7.4lf %7.4lf \n",l,ds.pert[l]*1e6,ds.pert[l]*1e6,csigma[l]*1e6);

    if (pertDebug==1) {
    std::ostringstream vol;
    vol << "vol_"  << obsID << ".txt";
    std::ofstream vol1(vol.str().c_str());
    for (long i = ii2; i <= ii1; i++) {
        for (long j = jj2; j <= jj1; j++) {
            for (long k = kk2; k <= kk1; k++) {
                long idx = i*N*N + j*N + k;
                if (ds.e[idx] == 1) {
                    vol1 << std::scientific << std::setprecision(17) <<
                        ds.x[idx] << " " << ds.y[idx] << " " << ds.z[idx] << " "
                         << ds.xoo[idx] << " " << ds.yoo[idx] << " " << ds.zoo[idx] << " "
                         << ds.cx[idx] << " " << i << " " << j << " " << k << std::endl;
                }
            }
        }
    }
    vol1.close();
    }

}
