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
from multiprocessing.dummy import Pool as ThreadPool

# Root directory of the project
ROOT_DIR = os.path.abspath("./Mask_RCNN")
OUT_DIR = os.path.abspath("./trainingset")

# Import Mask RCNN
sys.path.append(ROOT_DIR)  # To find local version of the library
from mrcnn.config import Config
from mrcnn import utils
import mrcnn.model as modellib
from mrcnn import visualize
from mrcnn.model import log

# Directory to save logs and trained model
MODEL_DIR = os.path.join(ROOT_DIR, "logs")

# Local path to trained weights file
COCO_MODEL_PATH = os.path.join(ROOT_DIR, "mask_rcnn_coco.h5")
# Download COCO trained weights from Releases if needed
if not os.path.exists(COCO_MODEL_PATH):
    utils.download_trained_weights(COCO_MODEL_PATH)


## CONFIG

class DESConfig(Config):

    # Give the configuration a recognizable name
    NAME = "DES"

    # Train on 1 GPU and 8 images per GPU. We can put multiple images on each
    # GPU because the images are small. Batch size is 8 (GPUs * images/GPU).
    GPU_COUNT = 1
    IMAGES_PER_GPU = 8

    # Number of classes (including background)
    NUM_CLASSES = 1 + 2  # background + star and galaxy

    # Use small images for faster training. Set the limits of the small side
    # the large side, and that determines the image shape.
    IMAGE_MIN_DIM = 512
    IMAGE_MAX_DIM = 512

    # Use smaller anchors because our image and objects are small
    RPN_ANCHOR_SCALES = (8, 16, 32, 64, 128)  # anchor side in pixels

    # Reduce training ROIs per image because the images are small and have
    # few objects. Aim to allow ROI sampling to pick 33% positive ROIs.
    TRAIN_ROIS_PER_IMAGE = 500

    # Use a small epoch since the data is simple
    STEPS_PER_EPOCH = 100

    # use small validation steps since the epoch is small
    VALIDATION_STEPS = 5

    # Store masks inside the bounding boxes (looses some accuracy but speeds up training)
    USE_MINI_MASK = True

class PhoSimDataset(utils.Dataset):

    def load_sources(self, max_num=None):
        # load specifications for image Dataset
        # follows load_shapes example
        black = (0,0,0)
        height = 512
        width = 512
        # add DES classes
        self.add_class("des", 1, "star")
        self.add_class("des", 2, "galaxy")

        # find number of sets:
        num_sets = 0
        for setdir in os.listdir(OUT_DIR):
            if 'set_' in setdir:
                num_sets += 4
                if max_num is not None and num_sets > max_num:
                    break

        # add image ids and specs from phosim output directory in order
        for i in range(num_sets):
            setdir = 'set_%d' % (i // 4) # read same set dir for all 4
            sources = 0 # count sources
            for image in os.listdir(os.path.join(OUT_DIR,setdir)):
                if image.endswith('.fits.gz') and not 'img' in image:
                    #if none on chip
                    # rename masks with source id number
                    if 'star' in image and not 'star_' in image:
                        image = os.path.join(OUT_DIR,setdir,image)
                        os.rename(image, image.split('star')[0]+"star_"+str(sources)+".fits.gz")
                    elif 'gal' in image and not 'gal_' in image:
                        image = os.path.join(OUT_DIR,setdir,image)
                        os.rename(image, image.split('gal')[0]+"gal_"+str(sources)+".fits.gz")
                    sources += 1
            # add tranining image set
            self.add_image("des", image_id=i, path=None,
                    width=width, height=height,
                    bg_color=black, sources=sources)

    def load_image(self, image_id):
        # load image set via image_id from phosim output directory
        # each set directory contains seperate files for images and masks
        info = self.image_info[image_id]
        if image_id % 4 == 0:
            for setdir in os.listdir(OUT_DIR):
                print(setdir)
                if setdir == 'set_%d' % image_id:
                    # image loop
                    for image in os.listdir(os.path.join(OUT_DIR,setdir)):
                        if image.endswith('.fits.gz') and 'img_g' in image:
                            g = getdata(os.path.join(OUT_DIR,setdir,image))
                            g /= np.max(g)
                            g *= 255
                        elif image.endswith('.fits.gz') and 'img_r' in image:
                            r = getdata(os.path.join(OUT_DIR,setdir,image))
                            r /= np.max(r)
                            r *= 255
                        elif image.endswith('.fits.gz') and 'img_i' in image:
                            i = getdata(os.path.join(OUT_DIR,setdir,image))
                            i /= np.max(i)
                            i *= 255
                else:
                    print('Warning: set_%d not found!' % image_id:)
                    return
            # convert format
            image = np.zeros([info['height'], info['width'], 3], dtype=np.uint8)
            image[:,:,0] = np.swapaxes(g,0,1) # b
            image[:,:,1] = np.swapaxes(r,0,1) # g
            image[:,:,2] = np.swapaxes(i,0,1) # r
        else: # get other 2 mirrors
            image = self.load_image(image_id//4)
        image = np.rot90(image,1+image_id%4)
        return image

    def read_mask(self,image):
        thresh = 0.00001
        if image.endswith('.fits.gz') and not 'img' in image:
            try: data = getdata(image)
            except: return
            mask_temp = np.zeros([data.shape[0],data.shape[1]], dtype=np.uint8)
            # Normalize
            max_data = np.max(data)
            # Check for empty mask (falls off chip)
            if max_data == 0: return
            data /= max_data
            # Apply threshold
            inds = np.argwhere(data>thresh)
            for ind in inds:
                mask_temp[ind[0],ind[1]] = 1
            # Gaussian blur
            mask_temp = cv2.GaussianBlur(mask_temp,(9,9),2)
            if 'star' in image:
                i = int(image.split('star_')[1].split('.')[0])
                self.class_ids[i] = 1
            elif 'gal' in image:
                i = int(image.split('gal_')[1].split('.')[0])
                self.class_ids[i] = 2
            self.mask[:,:,i] = mask_temp.astype(np.bool)
            return

    def load_mask(self, image_id):
        info = self.image_info[image_id]
        if image_id % 4 == 0:
            # load image set via image_id from phosim output directory
            sources = info['sources'] # number of sources in image
            self.mask = np.full([info['height'], info['width'], sources], False)
            # load image set via image_id from phosim output directory
            # each set directory contains seperate files for images and masks
            self.class_ids = np.zeros(sources,dtype=np.uint8)
            for setdir in os.listdir(OUT_DIR):
                if setdir == 'set_%d' % image_id:
                     image = os.listdir(os.path.join(OUT_DIR,setdir))
                     image_full = [os.path.join(OUT_DIR,setdir,f) for f in image]
                     pool = ThreadPool(128)
                     out = pool.map(self.read_mask, image_full)
                     pool.close()
                     pool.join()
        else: # get other 3 rotations
            self.mask,self.class_ids = self.load_mask(image_id//4)
        self.mask = np.rot90(self.mask,1+image_id%4,axes=(0,1))
        return self.mask.astype(np.bool), self.class_ids.astype(np.int32)

def train():

    #TODO: add option for small medium or large training set

    config = DESConfig()
    config.display()

    ## DATASET

    # Training dataset
    dataset_train = PhoSimDataset()
    dataset_train.load_sources()
    dataset_train.prepare()

    # Validation dataset
    dataset_val = PhoSimDataset()
    dataset_val.load_sources()
    dataset_val.prepare()

    # Load and display random samples
    image_ids = np.random.choice(dataset_train.image_ids, 1)
    for image_id in image_ids:
        print("loading image to visualize. Be patient...")
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
    model_path = os.path.join(MODEL_DIR, "mask_rcnn_des.h5")
    model.keras_model.save_weights(model_path)

if __name__ == "__main__":
    train()
