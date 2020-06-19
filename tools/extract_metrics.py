import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt



def get_fits_intensity(fp):
    data = 0
    with fits.open(fp) as hdul:
        data = hdul[0].data
    intensity = np.sum(data)
    return intensity

def intensity_hist(fp):
    data = 0
    with fits.open(fp) as hdul:
        data = hdul[0].data
    _ = plt.hist(data[0:256,0:256].flatten(), bins = "auto")
    plt.show()
    return None
    













if __name__ == "__main__":
    ROOT_PATH = "/Users/anshul/Desktop/astro_rcnn/hsc_test/set_0/"
    print(intensity_hist(os.path.join(ROOT_PATH, 'img_g_l.fits')))
    #print(intensity_hist(os.path.join(ROOT_PATH, 'img_r.fits')))
    #print(intensity_hist(os.path.join(ROOT_PATH, 'img_g.fits')))
    #print(intensity_hist(os.path.join(ROOT_PATH, 'img_z.fits')))
    pass
