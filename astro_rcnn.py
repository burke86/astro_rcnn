# Run mask-RCNN on phosim images
#
# Written by Colin J. Burke (Illinois)
# Adapted from Mask_RCNN/samples/nucleus/nucleus.py

# this assumes simulate.py has been run already
# or the training data is in the ./trainingset directory
# and validation data is in the ./validationset directory

# Type ./astro_rcnn -h for usage

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
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool
from imgaug import augmenters as iaa

# Root directory of the project
ROOT_DIR = os.path.abspath("./Mask_RCNN")
TRAIN_DIR = os.path.abspath("./trainingset")
VAL_DIR = os.path.abspath("./validationset")
TEST_DIR = os.path.abspath("./testset")
REAL_DIR = os.path.abspath("./realset")

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

class DESConfig(Config):

    # Give the configuration a recognizable name
    NAME = "DES"

    # Batch size (images/step) is (GPUs * images/GPU).
    GPU_COUNT = 4
    IMAGES_PER_GPU = 6

    # Number of classes (including background)
    NUM_CLASSES = 1 + 2  # background + star and galaxy

    # Use small images for faster training. Set the limits of the small side
    # the large side, and that determines the image shape.
    IMAGE_MIN_DIM = 512
    IMAGE_MAX_DIM = 512

    # Use smaller anchors because our image and objects are small
    RPN_ANCHOR_SCALES = (8, 16, 32, 64, 128)  # anchor side in pixels
    
    # How many anchors per image to use for RPN training
    RPN_TRAIN_ANCHORS_PER_IMAGE = 512

    # Reduce training ROIs per image because the images are small and have
    # few objects. Aim to allow ROI sampling to pick 33% positive ROIs.
    TRAIN_ROIS_PER_IMAGE = 250

    # Maximum number of ground truth instances (objects) in one image
    MAX_GT_INSTANCES = 300
    DETECTION_MAX_INSTANCES = 300

    # Mean pixel values (RGB)
    MEAN_PIXEL = np.array([444, 289, 161])

    # Note the images per epoch = steps/epoch * images/GPU * GPUs 
    # So the training time is porportional to the batch size
    # Use a small epoch since the batch size is large
    STEPS_PER_EPOCH = max(1, 1000 // (IMAGES_PER_GPU * GPU_COUNT))

    # Use small validation steps since the epoch is small
    VALIDATION_STEPS = max(1, 50 // (IMAGES_PER_GPU * GPU_COUNT))

    # Store masks inside the bounding boxes (looses some accuracy but speeds up training)
    USE_MINI_MASK = True
   


class InferenceConfig(DESConfig):
    GPU_COUNT = 1
    IMAGES_PER_GPU = 1
    DETECTION_MIN_CONFIDENCE = 0.5

class PhoSimDataset(utils.Dataset):

    def load_sources(self, set_dir, dataset="validation"):
        # dataset should be "validation", "training", "test", or "real"
        self.dataset = dataset
        self.out_dir = set_dir
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
        for setdir in os.listdir(self.out_dir):
            if 'set_' in setdir:
                num_sets += 1

        # add image ids and specs from phosim output directory in order
        for i in range(num_sets):
            setdir = 'set_%d' % i # read same set dir for all 4
            sources = 0 # count sources
            for image in os.listdir(os.path.join(self.out_dir,setdir)):
                if (image.endswith('.fits.gz') or image.endswith('.fits')) and not 'img' in image:
                    # rename masks with source id number
                    if 'star' in image and not 'star_' in image:
                        image = os.path.join(self.out_dir,setdir,image)
                        os.rename(image, image.split('star')[0]+"star_"+str(sources)+".fits.gz")
                    elif 'gal' in image and not 'gal_' in image:
                        image = os.path.join(self.out_dir,setdir,image)
                        os.rename(image, image.split('gal')[0]+"gal_"+str(sources)+".fits.gz")
                    sources += 1
            # add tranining image set
            self.add_image("des", image_id=i, path=None,
                    width=width, height=height,
                    bg_color=black, sources=sources)

    def load_image(self, image_id, scale=True):
        # load image set via image_id from phosim output directory
        # each set directory contains seperate files for images and masks
        info = self.image_info[image_id]
        setdir = 'set_%d' % image_id
        # adjust scale (do not use for training)
        if self.dataset == "real":
            norm_max = 50
            zero = .005
        else:
            norm_max = 180000
            zero = 0.0
        # image loop
        found = 0
        for image in os.listdir(os.path.join(self.out_dir,setdir)):
            if image == "img_g.fits" or image == "img_g.fits.gz":
                g = getdata(os.path.join(self.out_dir,setdir,image))
                if scale:
                    g = 65535*(g - np.mean(g))/np.std(g)
                else:
                    g = np.add(g,zero)
                    g = g*65535/norm_max
                found += 1
            elif image == "img_r.fits" or image == "img_r.fits.gz":
                r = getdata(os.path.join(self.out_dir,setdir,image))
                if scale:
                    r = 65535*(r - np.mean(r))/np.std(r)
                else:
                    r = np.add(r,zero)
                    r = r*65535/norm_max
                found += 1
            elif image == "img_z.fits" or image == "img_z.fits.gz":
                z = getdata(os.path.join(self.out_dir,setdir,image))
                if scale:
                    z = 65535*(z - np.mean(z))/np.std(z)
                else:
                    z = np.add(z,zero)
                    z = z*65535/norm_max
                found += 1
            # found all 3 bands
            if found == 3:
                break
        if found != 3:
            print("WARNING: Did not find files for all 3 bands.")
        # convert format
        image = np.zeros([info['height'], info['width'], 3], dtype=np.uint16)
        image[:,:,0] = z # red
        image[:,:,1] = r # green
        image[:,:,2] = g # blue
        return image

    def read_mask(self,image):
        # read and create mask
        if (image.endswith('.fits.gz') or image.endswith('.fits')) and not 'img' in image:
            data = getdata(image)
            mask_temp = np.zeros([data.shape[0],data.shape[1]], dtype=np.uint8)
            # Normalize
            max_data = np.max(data)
            # Check for empty mask (falls off chip)
            if max_data == 0: return
            data /= max_data
            if 'star' in image:
                i = int(image.split('star_')[1].split('.')[0])
                self.class_ids[i] = 1
                thresh = 0.005
            elif 'gal' in image:
                i = int(image.split('gal_')[1].split('.')[0])
                self.class_ids[i] = 2
                thresh = 0.08
            # Apply threshold (can make multiple contours)
            inds = np.argwhere(data>thresh)
            for ind in inds:
                mask_temp[ind[0],ind[1]] = 1
            mask_temp = cv2.GaussianBlur(mask_temp,(9,9),2)
            self.mask[:,:,i] = mask_temp.astype(np.bool)
            return

    def load_mask(self, image_id):
        info = self.image_info[image_id]
        # load image set via image_id from phosim output directory
        sources = info['sources'] # number of sources in image
        self.mask = np.full([info['height'], info['width'], sources], False)
        # load image set via image_id from phosim output directory
        # each set directory contains seperate files for images and masks
        self.class_ids = np.zeros(sources,dtype=np.uint8)
        setdir = 'set_%d' % image_id
        image = os.listdir(os.path.join(self.out_dir,setdir))
        # use all threads to load masks
        image_full = [os.path.join(self.out_dir,setdir,f) for f in image]
        threads = np.clip(mp.cpu_count()-2,0,None)
        pool = ThreadPool(threads)
        out = pool.map(self.read_mask, image_full)
        pool.close()
        pool.join()
        return self.mask.astype(np.bool), self.class_ids.astype(np.int32)

def train():

    start_time = time.time()

    config = DESConfig()
    config.display()

    ## DATASET

    # Training dataset
    dataset_train = PhoSimDataset()
    dataset_train.load_sources(TRAIN_DIR)
    dataset_train.prepare()

    # Validation dataset
    dataset_val = PhoSimDataset()
    dataset_val.load_sources(VAL_DIR)
    dataset_val.prepare()

    # Image augmentation
    augmentation = iaa.SomeOf((0, 4), [
        iaa.Fliplr(0.5),
        iaa.Flipud(0.5),
        iaa.OneOf([iaa.Affine(rotate=90),
                   iaa.Affine(rotate=180),
                   iaa.Affine(rotate=270)]),
        iaa.GaussianBlur(sigma=(0.0, 3.0)),
        iaa.AddElementwise((-25, 25))
    ])
    
    # Create model in training mode
    model = modellib.MaskRCNN(mode="training", config=config,
                              model_dir=MODEL_DIR)

    # Which weights to start with?
    init_with = "last"  # imagenet, coco, or last

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
    #model.train(dataset_train, dataset_val,
    #            learning_rate=config.LEARNING_RATE,
    #            augmentation=augmentation,
    #            epochs=15,
    #            layers='heads')

    # Fine tune all layers
    # Passing layers="all" trains all layers. You can also
    # pass a regular expression to select which layers to
    # train by name pattern.
    model.train(dataset_train, dataset_val,
                learning_rate=config.LEARNING_RATE / 10,
                augmentation=augmentation,
                epochs=40,
                layers="all")

    # Save weights
    model_path = os.path.join(MODEL_DIR, "astro_rcnn_des.h5")
    model.keras_model.save_weights(model_path)
    
    print("Done in %.2f hours." % float((time.time() - start_time)/3600))

    return

def detect(mode="simulated"):

    print("Model in inference mode.")
    inference_config = InferenceConfig()

    # Use most recent weight file
    # Add code to download it if it does not exist
    model_path = os.path.join(MODEL_DIR, "astro_rcnn_des.h5")

    # Recreate the model in inference mode
    model = modellib.MaskRCNN(mode="inference", 
                          config=inference_config,
                          model_dir=MODEL_DIR)

    # Load trained weights
    print("Loading weights from ", model_path)
    model.load_weights(model_path, by_name=True)

    dataset = PhoSimDataset()

    if mode == "real":
        # Load real DES image of abell cluster
        dataset.load_sources(REAL_DIR)
    elif mode == "simulated":
        # Load images from test set
        dataset.load_sources(TEST_DIR)
    else:
        print("Inference mode not recongnized.")  

    dataset.prepare()
    
    # Plot precision-recall curve over range of IOU thresholds
    mean_APs_star = []
    mean_ps_star = []
    mean_rs_star = []
    mean_APs_gal = []
    mean_ps_gal = []
    mean_rs_gal = []
    iou_thresholds = np.arange(0.5, 1.0, 0.05)
    for i in iou_thresholds:
        # star
        APs,ps,rs = utils.compute_performance(dataset,model,inference_config,1,i)
        mean_APs_star.append(APs)
        mean_ps_star.append(ps)
        mean_rs_star.append(rs)
        # galaxy
        APs,ps,rs = utils.compute_performance(dataset,model,inference_config,2,i)
        mean_APs_gal.append(APs)
        mean_ps_gal.append(ps)
        mean_rs_gal.append(rs)
    # Plot precision-recall
    visualize.plot_precision_recall_range(mean_APs_star,iou_thresholds,mean_ps_star,mean_rs_star,save_fig=True,title="star")
    visualize.plot_precision_recall_range(mean_APs_gal,iou_thresholds,mean_ps_gal,mean_rs_gal,save_fig=True,title="galaxy")

    # Loop over images
    submission = []
    for image_id in dataset.image_ids:
        # Load image and run detection
        image = dataset.load_image(image_id)
        # Detect objects
        r = model.detect([image], verbose=0)[0]
        # Encode image to RLE. Returns a string of multiple lines
        source_id = dataset.image_info[image_id]["id"]
        # Save image predictions with masks
        visualize.display_instances(
            image, r['rois'], r['masks'], r['class_ids'],
            dataset.class_names, r['scores'],
            show_bbox=False, show_mask=False,
            title="Predictions")
        plt.savefig("{}/{}.png".format(MODEL_DIR, dataset.image_info[image_id]["id"]))
        # Save mask as sextractor-style .fits file

    return

if __name__ == "__main__":
    import argparse
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Mask R-CNN for star/galaxy detection, classification, and deblending')
    parser.add_argument("command",metavar="<command>",help="'train', 'detect', or 'detect_real'")
    args = parser.parse_args()

    # Train or evaluate
    if args.command == "train":
        train()
    elif args.command == "detect":
        detect("simulated")
    elif args.command == "detect_real":
        detect("real")
    else:
        print("'{}' is not recognized. "
              "Use 'train', 'detect', or 'detect_real'".format(args.command))

