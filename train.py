# train mask-RCNN on phosim images

# this assumes ./phosim_release/generate_training_set.py has been run already
# and fits images and masks are sitting in ./phosim_release/output directory

# setup
import os
import sys
import gzip
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

#%matplotlib inline

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
        output = './phosim_release/output/'
        for setdir in os.listdir(output):
            sources = 0
            # set_X
            if 'set' in setdir:
                # add tranining image set
                self.add_image("des", image_id=i, path=None,
                        width=width, height=height,
                        bg_color=black, sources=sources)
                i += 1


    def load_image(self, image_id):
        # load image set via image_id from phosim output directory
        # each set directory contains seperate files for images and masks
        info = self.image_info[image_id]
        output = './phosim_release/output/'
        for setdir in os.listdir(output):
            # image loop
            search_str = 'set_%d' % image_id
            if search_str in setdir:
                # image loop
                for image in os.listdir(os.path.join(output,setdir)):
                    if image.endswith('.fits.gz') and 'img' in image:
                        print(image)
                        data = getdata(os.path.join(output,setdir,image))
                        break
        # convert format
        image = np.array(data).reshape([1, 1, 3]).astype(np.uint8)
        return image

    def image_reference(self, image_id):
        # return the shapes data of the image
        info = self.image_info[image_id]
        if info["source"] == "sources":
            return info["sources"]
        else:
            super(self.__class__).image_reference(self, image_id)

    def load_mask(self, image_id):
        info = self.image_info[image_id]

        # load image set via image_id from phosim output directory
        red = (255,0,0) # star mask color
        blue = (0,0,255) # galaxy mask color
        threshold = 0.01 # pixel values above this % of the max value in the
        sources = info['sources'] # number of sources in image
        count = len(sources)
        mask = np.zeros([info['height'], info['width'], count], dtype=np.uint8)

        # load image set via image_id from phosim output directory
        # each set directory contains seperate files for images and masks
        output = './phosim_release/output/'
        for setdir in os.listdir(output):
            # image loop
            search_str = 'set_%d' % image_id
            if search_str in setdir:
                # image loop
                for image in os.listdir(os.path.join(output,setdir)):
                    if image.endswith('.fits.gz') and not 'img' in image:
                        print(image)
                        data = data = getdata(os.path.join(output,setdir,image))
                        mask *= np.where(data/np.max(data) > threshold)
                        break
        # occulsions? colors?
        return mask


# Training dataset
dataset_train = PhoSimDataset()
dataset_train.load_sources()
dataset_train.prepare()

# Load and display random samples
image_ids = np.random.choice(dataset_train.image_ids, 1)
for image_id in image_ids:
    image = dataset_train.load_image(image_id)
    mask, class_ids = dataset_train.load_mask(image_id)
    #visualize.display_top_masks(image, mask, class_ids, dataset_train.class_names)
