import os
import numpy as np
from astropy.io import fits


def save_znorm(in_fp, out_fp):
    data = 0
    with fits.open(in_fp) as hdul:
        data = hdul[0].data
    zscore_data = 1e4 * ((data - np.mean(data))/np.std(data))
    hdu = fits.PrimaryHDU(zscore_data)
    hdul_new = fits.HDUList([hdu])
    hdul_new.writeto(out_fp)
    return None

def save_logscale(in_fp, out_fp):
    data = 0
    with fits.open(in_fp) as hdul:
        data = hdul[0].data
    log_scale_data = np.log(np.abs((1e3) * data) + 1) / np.log(1e4)
    hdu = fits.PrimaryHDU(log_scale_data)
    hdul_new = fits.HDUList([hdu])
    hdul_new.writeto(out_fp)
    return None
    

if __name__ == "__main__":
   #save_znorm('img_g.fits','img_g_z.fits') 
   save_logscale('img_g.fits', 'img_g_l.fits')
