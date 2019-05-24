import numpy as np
import os
import sys

""" 
To download images in the NOAO archive 
The image will be named as 'Jxxxx+xxxx.fits' based on ra and dec
"""

def _make_name(ra,dec):

    # generate the name for the objects (e.g. Jxxxx+xxxx)

    ################
    # Fill in code #
    ################

    return name


def get_image(ra,dec,size):

    # set up some values
    band = "g"
    name = _make_name(ra,dec) 

    query = "" % (ra,dec,band)

    # download the image
    os.system(query)

    # change the name
    os.system("mv out.fits %s.fits" % name)
    
    ################
    # Fill in code #
    ################

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Incorrect number of parameters !!")
        print("Input Format:python download_image.py [ra(deg)] [dec(deg)] "+\
              "[size(arcsec)]")
    else:
        ra = sys.argv[1]
        dec = sys.argv[2]
        size = sys.argv[3]
        get_image(ra,dec,size)
