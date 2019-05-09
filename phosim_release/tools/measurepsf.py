##
## @package phosim
## @file phosim
## @brief main script to measurepsf.py
##
## @brief Created by:
## @author Jun Cheng (Purdue)
##
##
## @warning This code is not fully validated
## and not ready for full release.  Please
## treat results with caution.
##
## Requirement:
## python library:  astropy.
##
## Example:
## python measurepsf.py cut_test.fits
##
## This code will return: e1, e2, ellip, pa, flux,rms



from astropy.io import fits
import math

def readFitsImage(imageName):
    hdulist = fits.open(imageName)
    data =  hdulist[0].data
    hdulist.close()
    return data

def getEllipticity(data):
    size = data.shape
    t1 = 0
    t2 = 0
    t3 = 0
    t4 = 0
    t5 = 0
    t6 = 0

    f=open("data.txt", 'w')
    f.write("[")
    for i in range(size[0]):
        f.write("[")
        for j in range(size[1]):
            f.write(str(data[i][j]))
            if j!=size[1]-1:
                f.write(",")
            t1 += i*j*data[i][j]
            t2 += data[i][j]
            t3 += i*data[i][j]
            t4 += j*data[i][j]
            t5 += i*i*data[i][j]
            t6 += j*j*data[i][j]
        f.write("],")
    f.write("]\n")
    covxy = t1/t2-t3*t4/t2/t2
    resultx = t5/t2-t3*t3/t2/t2
    resulty = t6/t2-t4*t4/t2/t2

    medx = size[0]/2.0
    medy = size[1]/2.0
    alphax = math.sqrt(resultx*2.0)
    alphay = math.sqrt(resulty*2.0)
    alphaxy = covxy*math.sqrt(2.0)

    for trials in range(100):
        t1 = 0
        t2 = 0
        t3 = 0
        t4 = 0
        t5 = 0
        t6 = 0
        t7 = 0
        for i in range(size[0]):
            for j in range(size[1]):
                weight = (math.exp(-((i-medx)**2/alphax/alphax-2.0*alphaxy/alphax/alphax/alphay/alphay*(i-medx)*(j-medy)+(j-medy)**2/alphay/alphay)/2.0/(1-(alphaxy/alphax/alphay)**2)))/(2*math.pi*alphax*alphay*math.sqrt(1-(alphaxy/alphax/alphay)**2))

                t1 += i*j*data[i][j]*weight
                t2 += data[i][j]*weight
                t3 += i*data[i][j]*weight
                t4 += j*data[i][j]*weight
                t5 += i*i*data[i][j]*weight
                t6 += j*j*data[i][j]*weight
                t7 += data[i][j]

        covxy=(t1/t2-t3*t4/t2/t2)
        resultx=(t5/t2-t3*t3/t2/t2)
        resulty=(t6/t2-t4*t4/t2/t2)
        medx=(t3/t2)
        medy=(t4/t2)
        flux=t7

        rms=math.sqrt(resultx+resulty)

        e1=(resulty-resultx)/(resultx+resulty)
        e2=(2.0*covxy)/(resultx+resulty)

        ellip= math.sqrt(e1**2+e2**2)
        pa=0.5*math.degrees(math.atan(e2/e1))-90

        if (abs(alphax-math.sqrt(2.0*resultx)) < 1e-6 and abs(alphay-math.sqrt(2.0*resulty)) < 1e-6 and abs(alphaxy-2.0*covxy) < 1e-6) :
            break
        alphax=math.sqrt(resultx*2.0)
        alphay=math.sqrt(resulty*2.0)
        alphaxy=covxy*2.0

    return e1, e2, ellip, pa, flux,rms


def main():
    imageName = 'cut_test.fits'
    data = readFitsImage(imageName)
    e1, e2, ellip, pa, flux, rms = getEllipticity(data)
    print "e1 = ", e1
    print "e2 = ", e2
    print "e  = ", ellip
    print "Position angle = ", pa
    print "Flux = ", flux
    print "RMS = ", rms


if __name__=='__main__':
    main()
