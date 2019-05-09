
#include <math.h>
#include <stdio.h>

#include "../source/raytrace/helpers.h"
#include "../source/raytrace/raytrace.h"

#include "unittest.h"

int main() {

    double xx[100],yy[100];
    float ff[100];
    double x,y,z;
    float f;
    long index;
    long i,j;
    Vector v, normal, position, angle;

    printf("\n");
    printf("-------------------------------------\n");
    printf("    PhoSim Unit Testing Framework\n");
    printf("-------------------------------------\n\n");
    for (i=0;i<100;i++) xx[i]=((double)i)/100.0;
    for (i=0;i<100;i++) yy[i]=((double)i)/100.0;
    for (i=0;i<100;i++) ff[i]=((double)i)/100.0;


    find(xx,100,0.601,&index);
    unitTestOutput(index,60,"find","Index finding","none");

    index=find_linear(xx,100,0.601,&x);
    unitTestOutput(index,60,"findLinear","Index finding","none");

    find_linear_wrap(1.0,0.1,100,&i,&j,&x);
    unitTestOutput(j,60,"findLinearWrap","Index finding","none");

    x=1.0; y=1.0; z=0.0;
    normalize(&x,&y,&z);
    unitTestOutput(x,1.0/sqrt(2.0),"normalize","Normalize vector","none");

    x=interpolate(yy,xx,0.605,60);
    unitTestOutput(x,0.605,"interpolate","Interpolate array","none");

    x=interpolate_linear(yy,60,60.5);
    unitTestOutput(x,0.605,"interpolateLinear","Interpolate array","none");

    x=interpolate_bilinear(yy,0,0,0.0,60,60.5);
    unitTestOutput(x,0.605,"interpolateBilinear","Interpolate 2-D array","none");

    f=interpolate_bilinear_float_wrap(ff,0,0,1,0.0,60,61,0.5);
    unitTestOutput((double)f,0.605,"interpolateBilinearFloatWrap","Interpolate 2-D array","none");

    v.x=0.0; v.y=1.0/sqrt(2.0); v.z=-1.0/sqrt(2.0);
    normal.x=0.0; normal.y=0.0; normal.z=1.0;
    reflect(&v,normal);
    unitTestOutput(v.z,1/sqrt(2.0),"reflect","Ray component after reflection","none");

    v.x=0.0; v.y=1.0/sqrt(2.0); v.z=-1.0/sqrt(2.0);
    normal.x=0.0; normal.y=0.0; normal.z=1.0;
    refract(&v,normal,1.0,2.0);
    unitTestOutput(v.y,1/sqrt(2.0)/2.0,"refract","Ray component after snell's law","none");

    position.x=0.; position.y=0.; position.z=0.;
    angle.x=0; angle.y=1/sqrt(2.0); angle.z=1/sqrt(2.0);
    propagate(&position,angle,sqrt(2.0));
    unitTestOutput(position.y,1.0,"propagate","Ray position after propagation","none");

    return(0);
}
