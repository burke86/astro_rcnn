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
sys.path.append("./")
sys.path.append(ROOT_DIR)  # To find local version of the library
from mrcnn.config import Config
from mrcnn import utils
import mrcnn.model as modellib
from mrcnn import visualize
from mrcnn.model import log

# Directory to save logs and trained model
MODEL_DIR = os.path.join(ROOT_DIR, "logs")

from train import DESConfig

class InferenceConfig(DESConfig):
    GPU_COUNT = 1
    IMAGES_PER_GPU = 1

def detect():

    inference_config = InferenceConfig()

    # Recreate the model in inference mode
    model = modellib.MaskRCNN(mode="inference",
                              config=inference_config,
                              model_dir=MODEL_DIR)

    # Get path to saved weights
    model_path = os.path.join(MODEL_DIR, "mask_rcnn_des.h5")

    # Load trained weights
    print("Loading weights from ", model_path)
    model.load_weights(model_path, by_name=True)

    # Load DES Y1 image as test
    data = get_data('./c4d_130830_034536_ooi_r_d1.fits.fz',1)
    data /= np.max(des_image)
    data *= 255
    # convert format
    image = np.ones([data.shape[0], data.shape[1], 3], dtype=np.uint8)
    image[:,:,0] = data
    image[:,:,1] = data
    image[:,:,2] = data
    image = np.flip(image,0)

    # Test on a full DES image
    results = model.detect([image], verbose=1)

    r = results[0]
    visualize.display_instances(image, r['rois'], r['masks'], r['class_ids'],
                                dataset_val.class_names, r['scores'], ax=get_ax())

if __name__ == "__main__":
    detect()
