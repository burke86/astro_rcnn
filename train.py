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

# number of training image (and associtaed masks) expected
N_TRAINING = 1

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

# get images from phosim output directory
def load_images():
    i = 0
    for image_name in os.path.listdir('./phosim_release/output/')
        if image_name.endswith('training_img') and image_name.endswith('.fits'):
            i = image_name[-6]
            data = getdata(image_name)



    pass
