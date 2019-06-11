import numpy as np
import os
import sys

""" 
To download images in the NOAO archive 
The image will be named as 'Jxxxx+xxxx.fits' based on ra and dec
"""

def _make_name(ra,dec):

    # generate the name for the objects (e.g. Jxxxx+xxxx)    

    ra = ra/15.
    rah=int(ra)
    ram=int(60*(ra-rah))
    ras=3600.*(ra-rah-(ram/60.))

    if (dec >= 0.0):
        decd=int(dec)
        decm=int(60.*(dec-decd))
        decs=3600.*(dec-decd-(decm/60.))
        return "J%02d%02d%05.2f+%02d%02d%04.1f" % (rah,ram,ras,decd,decm,decs)
    else:
        decd=-int(dec)
        dec=-dec
        decm=int(60.*(dec-decd))
        decs=3600.*(dec-decd-(decm/60.))
        dec=-dec        
        return "J%02d%02d%05.2f-%02d%02d%04.1f" % (rah,ram,ras,decd,decm,decs)

def get_multi_images(ra,dec,size,output_dir):

    make_dir(output_dir)
    bands = ["g","r","z"]
    for band in bands:
        output_path = "%s/img_%s.fits" % (output_dir,band)
        get_image(ra,dec,size,band,output_path)

def get_image(ra,dec,size,band,output_path):

    # set up some values
    '''
    band = "g"
    name = _make_name(float(ra),float(dec)) 
    '''
    query = "wget 'http://legacysurvey.org/viewer/fits-cutout/?ra=%s&dec=%s&" % (ra,dec)+\
            "size=%s&layer=decals-dr7&pixscale=0.262&bands=%s' -O out.fits" % (size,band)

    # download the image
    os.system(query)

    # change the name
    os.system("mv out.fits %s" % output_path)

def check_set_number(testset_dir):

    number = len(os.listdir(testset_dir))

    return number

def make_dir(dir_name):

    if not os.path.exists(dir_name):
        print("-- Making dir: %s" % dir_name)
        os.system("mkdir %s" % dir_name)
    

if __name__ == "__main__":

    testset_dir = "testset"

    make_dir(testset_dir)

    if len(sys.argv) != 4 :
        print("Incorrect number of parameters !!")
        print("Input Format:python download_image.py [ra(deg)] [dec(deg)] "+\
              "[size(arcsec)] ")

    else:
        ra = sys.argv[1]
        dec = sys.argv[2]
        size = sys.argv[3]
        number = check_set_number(testset_dir)
        get_multi_images(ra,dec,size,"%s/set_%s" %(testset_dir,number))
