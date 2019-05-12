# train mask-RCNN on phosim images

# this assumes ./phosim_release/generate_training_set.py has been run already
# and fits images and masks are sitting in ./phosim_release/output directory

# setup
import os
import sys
import random
import math
import re
import time
import numpy as np
import cv2
import matplotlib
import matplotlib.pyplot as plt
from astropy.io.fits import getdata

# Root directory of the project
ROOT_DIR = os.path.abspath("./Mask_RCNN")

# Import Mask RCNN
sys.path.append(ROOT_DIR)  # To find local version of the library
from mrcnn.config import Config
from mrcnn import utils
import mrcnn.model as modellib
from mrcnn import visualize
from mrcnn.model import log

%matplotlib inline

# Directory to save logs and trained model
MODEL_DIR = os.path.join(ROOT_DIR, "logs")

# Local path to trained weights file
COCO_MODEL_PATH = os.path.join(ROOT_DIR, "mask_rcnn_coco.h5")
# Download COCO trained weights from Releases if needed
if not os.path.exists(COCO_MODEL_PATH):
    utils.download_trained_weights(COCO_MODEL_PATH)

class PhoSimDataset(utils.Dataset):

    def load_sources(self):
        # load specifications for image Dataset
        # follows load_shapes example
        black = (0,0,0)
        height = 128
        width = 128
        # add DES classes
        self.add_class("des", 1, "star")
        self.add_class("des", 2, "galaxy")

        # add image ids and specs from phosim output Directory
        i = 0
        for setdir in os.path.listdir('./phosim_release/output/'):
            # image loop
            sources = 0
            for image in os.path.listdir(setdir):
                # find masks
                if 'set' in set:
                    if image.endswith('.fits') and not image.contains('img'):
                        sources +=1
                    self.add_image("des", image_id=i, path=None,
                            width=width, height=height,
                            bg_color=black, sources=sources)
                    i += 1


    def load_image(self, image_id):
        # load image set via image_id from phosim output directory
        # each set directory contains seperate files for images and masks
        for setdir in os.path.listdir('./phosim_release/output/'):
            # image loop
            search_str = 'set%d' % image_id
            if search_str in setdir:
                for image in os.path.listdir(setdir):
                    # find image id
                    if image.endswith('.fits'):
                        # image
                        if image.contains('img'):
                            data = getdata(image)
                            break
        return data

    def load_mask(self, image_id):
        info = self.image_info[image_id]

        # load image set via image_id from phosim output directory
        red = (255,0,0) # star mask
        blue = (0,0,255) # galaxy mask
        threshold = 0.01 # pixel values above this % of the max value in the
        sources = info['sources'] # number of sources in image
        count = len(sources)
        mask = np.zeros([info['height'], info['width'], count], dtype=np.uint8)

        # load image set via image_id from phosim output directory
        # each set directory contains seperate files for images and masks
        for setdir in os.path.listdir('./phosim_release/output/'):
            # image loop
            search_str = 'set%d' % image_id
            if search_str in setdir:
                for image in os.path.listdir(setdir):
                    # find image id
                    if image.endswith('.fits') and not image.contains('img'):
                        data = getdata(image)
                        mask = np.where(data/np.max(data) > threshold)
                        break
        return mask
