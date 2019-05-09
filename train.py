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

    def load_images(self):
        # load simulated phosim load_images
        width = 128
        height = 128

        # add DES classes
        self.add_class("des", 1, "star")
        self.add_class("des", 2, "galaxy")

        # Add images
        # each set directory contains seperate files for images and masks
        for setdir in os.path.listdir('./phosim_release/output/'):
            # image loop
            for image in os.path.listdir(setdir):
                if image.endswith('.fits'):
                    # image
                    if image.contains('img'):
                        data = getdata(image)
                    # mask
                    else:
                        data = getdata(image)
                self.add_image("des", image_id=i, path=None,
                           width=width, height=height,
                           bg_color=bg_color, shapes=shapes)
