import os
import numpy as np
from astropy.io import fits


def getKey(filename):
    key = 0
    try:
        key = int(filename[0:3])
    except:
        try:
            key = int(filename[0:2])
        except:
            key = int(filename[0])
    return key


def reformat_filepaths(root_path):
    filepaths = os.listdir(root_path)
    os.chdir(ROOT_PATH)
    files = []
    for elem in filepaths:
        if elem.endswith(".fits"):
            files.append(elem)

    files.sort(key = lambda x: getKey(x))
    print(files)
    for i in range(0, len(files)):
        if (i % 3 == 0):
            os.system("mkdir set_" + str(int(i/3)))
            move_cmd = "mv " + files[i] + " set_" + str(int(i/3)) + "/img_r.fits"
            os.system(move_cmd)

        elif (i % 3 == 1):
            move_cmd = "mv " + files[i] + " set_" + str(int(i/3)) + "/img_g.fits"
            os.system(move_cmd)
        else:
            move_cmd = "mv " + files[i] + " set_" + str(int(i/3)) + "/img_z.fits"
            os.system(move_cmd)

def compress_to_layer(root_path, l):
    dir_names = os.listdir(root_path)
    for set_dir in dir_names:
        curr_files = os.listdir(os.path.join(root_path, set_dir))
        for filename in curr_files:
            filepath = os.path.join(os.path.join(root_path, set_dir), filename)
            data = 0
            with fits.open(filepath) as hdul:
                data = hdul[l].data

            remove_cmd = "rm " + filepath
            os.system(remove_cmd)
            hdu = fits.PrimaryHDU(data)
            hdul_new = fits.HDUList([hdu])
            hdul_new.writeto(filepath)


if __name__ == "__main__":
    ROOT_PATH = "/Users/anshul/Desktop/astro_rcnn/hsc_test"
    #reformat_filepaths(ROOT_PATH)
    compress_to_layer(ROOT_PATH, 1)

