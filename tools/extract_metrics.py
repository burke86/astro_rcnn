import os
import numpy as np
from astropy.io import fits



def get_fits_intensity(fp):
    data = 0
    with fits.open(fp) as hdul:
        data = hdul[0].data
    intensity = np.sum(data)
    return intensity












if __name__ == "__main__":
    ROOT_PATH = "../../example/set_0/"
    print(get_fits_intensity(os.path.join(ROOT_PATH, 'img_r.fits')))
    print(get_fits_intensity(os.path.join(ROOT_PATH, 'img_g.fits')))
    print(get_fits_intensity(os.path.join(ROOT_PATH, 'img_z.fits')))
    pass
