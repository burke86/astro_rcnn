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
from skimage.morphology import watershed

# Root directory of the project
ROOT_DIR = os.path.abspath("./Mask_RCNN")
OUT_DIR = os.path.abspath("./phosim_release/output/")

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


## CONFIG

class SourcesConfig(Config):
    """Configuration for training on the toy shapes dataset.
    Derives from the base Config class and overrides values specific
    to the toy shapes dataset.
    """
    # Give the configuration a recognizable name
    NAME = "sources"

    # Train on 1 GPU and 8 images per GPU. We can put multiple images on each
    # GPU because the images are small. Batch size is 8 (GPUs * images/GPU).
    GPU_COUNT = 1
    IMAGES_PER_GPU = 8

    # Number of classes (including background)
    NUM_CLASSES = 1 + 2  # background + star and galaxy

    # Use small images for faster training. Set the limits of the small side
    # the large side, and that determines the image shape.
    IMAGE_MIN_DIM = 128
    IMAGE_MAX_DIM = 128

    # Use smaller anchors because our image and objects are small
    RPN_ANCHOR_SCALES = (8, 16, 32, 64, 128)  # anchor side in pixels

    # Reduce training ROIs per image because the images are small and have
    # few objects. Aim to allow ROI sampling to pick 33% positive ROIs.
    TRAIN_ROIS_PER_IMAGE = 32

    # Use a small epoch since the data is simple
    STEPS_PER_EPOCH = 100

    # use small validation steps since the epoch is small
    VALIDATION_STEPS = 5

config = SourcesConfig()
config.display()

## DATASET

class PhoSimDataset(utils.Dataset):

    def load_sources(self, max_num=None):
        # load specifications for image Dataset
        # follows load_shapes example
        black = (0,0,0)
        height = 128
        width = 128
        # add DES classes
        self.add_class("des", 1, "star")
        self.add_class("des", 2, "galaxy")

        # find number of sets:
        num_sets = 0
        for setdir in os.listdir(OUT_DIR):
            if 'set_' in setdir:
                num_sets += 1
                if max_num is not None and num_sets > max_num:
                    break

        # add image ids and specs from phosim output directory in order
        for i in range(num_sets):
            for setdir in os.listdir(OUT_DIR):
                search_str = 'set_%d' % i
                if search_str in setdir:
                    # count sources
                    sources = 0
                    for image in os.listdir(os.path.join(OUT_DIR,setdir)):
                        if image.endswith('.fits.gz') and not 'img' in image:
                            sources += 1
                    # add tranining image set
                    self.add_image("des", image_id=i, path=None,
                            width=width, height=height,
                            bg_color=black, sources=sources)

    def load_image(self, image_id):
        # load image set via image_id from phosim output directory
        # each set directory contains seperate files for images and masks
        info = self.image_info[image_id]
        for setdir in os.listdir(OUT_DIR):
            search_str = 'set_%d' % image_id
            if search_str in setdir:
                # image loop
                for image in os.listdir(os.path.join(OUT_DIR,setdir)):
                    if image.endswith('.fits.gz') and 'img' in image:
                        data = getdata(os.path.join(OUT_DIR,setdir,image))
                        data /= np.max(data)
                        data *= 255
                        break
        # convert format
        image = np.ones([info['height'], info['width'], 3], dtype=np.uint8)
        image[:,:,0] = data
        image[:,:,1] = data
        image[:,:,2] = data
        image = np.flip(image,0)
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
        threshold = 0.01 # pixel values above this % of the max value in the
        sources = info['sources'] # number of sources in image
        mask = np.zeros([info['height'], info['width'], sources], dtype=np.uint8)
        markers = np.zeros([info['height'], info['width']], dtype=np.uint8)
        # load image set via image_id from phosim output directory
        # each set directory contains seperate files for images and masks
        class_ids = np.zeros(sources,dtype=np.uint8)
        for setdir in os.listdir(OUT_DIR):
            search_str = 'set_%d' % image_id
            if search_str in setdir:
                # image loop
                i = 0
                for image in os.listdir(os.path.join(OUT_DIR,setdir)):
                    if image.endswith('.fits.gz') and not 'img' in image:
                        data = getdata(os.path.join(OUT_DIR,setdir,image))
                        data /= np.max(data)
                        args = np.argwhere(data > threshold)
                        for arg in args:
                            mask[arg[0],arg[1],i] = 1
                        # Gaussian blur
                        mask[:,:,i] = cv2.GaussianBlur(mask[:,:,i],(25,25),2)
                        # re-apply threshold
                        args = np.argwhere(data > threshold)
                        for arg in args:
                            mask[arg[0],arg[1],i] = 1
                        if 'star' in image:
                            class_ids[i] = 1
                        elif 'gal' in image:
                            class_ids[i] = 2
                            x,y = np.unravel_index(np.argmax(data),data.shape)
                            markers[x,y] = 1
                        i += 1
        # galaxy-galaxy occlusions
        #for j in range(sources):
        #    if class_ids[j] == 1: continue
        #    image = self.load_image(image_id)[:,:,0]
        #    label = watershed(-image, markers, mask=np.sum(mask,2))
        #    if not j % 2 == 0:
        #        mask[:,:,j] = label
        mask = np.flip(mask,0)
        return mask, class_ids


# Training dataset
dataset_train = PhoSimDataset()
dataset_train.load_sources()
dataset_train.prepare()

# Validation dataset
dataset_val = PhoSimDataset()
dataset_val.load_sources(3)
dataset_val.prepare()

# Load and display random samples
image_ids = np.random.choice(dataset_train.image_ids, 4)
for image_id in image_ids:
    image = dataset_train.load_image(image_id)
    mask, class_ids = dataset_train.load_mask(image_id)
    visualize.display_top_masks(image, mask, class_ids, dataset_train.class_names)

## CREATE MODEL

# Create model in training mode
model = modellib.MaskRCNN(mode="training", config=config,
                          model_dir=MODEL_DIR)

# Which weights to start with?
init_with = "coco"  # imagenet, coco, or last

if init_with == "imagenet":
    model.load_weights(model.get_imagenet_weights(), by_name=True)
elif init_with == "coco":
    # Load weights trained on MS COCO, but skip layers that
    # are different due to the different number of classes
    # See README for instructions to download the COCO weights
    model.load_weights(COCO_MODEL_PATH, by_name=True,
                       exclude=["mrcnn_class_logits", "mrcnn_bbox_fc",
                                "mrcnn_bbox", "mrcnn_mask"])
elif init_with == "last":
    # Load the last model you trained and continue training
    model.load_weights(model.find_last(), by_name=True)

# Train the head branches
# Passing layers="heads" freezes all layers except the head
# layers. You can also pass a regular expression to select
# which layers to train by name pattern.
model.train(dataset_train, dataset_val,
            learning_rate=config.LEARNING_RATE,
            epochs=1,
            layers='heads')

# Fine tune all layers
# Passing layers="all" trains all layers. You can also
# pass a regular expression to select which layers to
# train by name pattern.
model.train(dataset_train, dataset_val,
            learning_rate=config.LEARNING_RATE / 10,
            epochs=2,
            layers="all")

# Save weights
# Typically not needed because callbacks save after every epoch
# Uncomment to save manually
# model_path = os.path.join(MODEL_DIR, "mask_rcnn_shapes.h5")
# model.keras_model.save_weights(model_path)
