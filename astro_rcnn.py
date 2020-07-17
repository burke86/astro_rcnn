"""
Astro-RCNN code written on top of Mask R-CNN

Written by Colin J. Burke, Anshul Shah (UIUC)
Adapted from Mask_RCNN/samples/nucleus/nucleus.py

"""

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
from astropy.io import fits
from astropy.io.fits import getdata
from astropy.visualization import make_lupton_rgb
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool
from imgaug import augmenters as iaa
#from photutils.isophote import Ellipse, EllipseGeometry
#from photutils.aperture import EllipticalAperture


# Root directory of the project
ROOT_DIR = os.path.abspath("./Mask_RCNN")

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
    IMAGES_PER_GPU = 4

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
    MEAN_PIXEL = np.array([-200, -200, -200])

    # Note the images per epoch = steps/epoch * images/GPU * GPUs
    # So the training time is porportional to the batch size
    # Use a small epoch since the batch size is large
    STEPS_PER_EPOCH = max(1, 1000 // (IMAGES_PER_GPU * GPU_COUNT))

    # Use small validation steps since the epoch is small
    VALIDATION_STEPS = max(1, 250 // (IMAGES_PER_GPU * GPU_COUNT))

    # Store masks inside the bounding boxes (looses some accuracy but speeds up training)
    USE_MINI_MASK = True



class InferenceConfig(DESConfig):
    GPU_COUNT = 1
    IMAGES_PER_GPU = 1
    DETECTION_MIN_CONFIDENCE = 0.5


class PhoSimDataset(utils.Dataset):

    def __init__(self, height=512, width=512, stretch=0.005, Q=10):
        self.height = height
        self.width = width
        self.stretch = stretch
        self.Q = Q
        super(PhoSimDataset, self).__init__()


    def load_sources(self, set_dir, dataset="validation", normalize="zscore", store_raw=False):
        # Load sources in dataset with proper id
        # This happens once, upon calling dataset.prepare()
        self.dataset = dataset
        self.out_dir = set_dir
        # load specifications for image Dataset
        # follows load_shapes example
        black = (0,0,0)
        # add DES classes
        self.add_class("des", 1, "star")
        self.add_class("des", 2, "galaxy")

        # find number of sets:
        num_sets = 0
        for setdir in os.listdir(self.out_dir):
            if 'set_' in setdir:
                # add tranining image set
                self.add_image("des", image_id=num_sets,path=os.path.join(self.out_dir,set_dir),
                    width=self.width,height=self.height,bg_color=black)
                num_sets += 1

        # store data in memory
        self.images = [None]*(num_sets)
        if store_raw:
            self.raws = [None]*(num_sets)

        self.masks = [None]*num_sets
        self.class_ids_mem = [None]*num_sets
        threads = np.clip(mp.cpu_count(),1,num_sets)
        print("Loading images from disk.")
        pool = ThreadPool(threads)
        pool.starmap(self.load_image_disk, [(i, normalize, store_raw) for i in range(num_sets)])
        if dataset == "training" or dataset == "validation":
            print("Loading masks from disk (this may take several minutes).")
            pool.map(self.load_mask_disk, range(num_sets))
        pool.close()
        pool.join()
        return

    def load_image(self, image_id, raw = False):
        if raw:
            return self.raws[image_id]
        else:
            return self.images[image_id]

    def load_image_disk(self, image_id, normalize='zscore', store_raw = False):
        # load from disk -- each set directory contains seperate files for images and masks
        info = self.image_info[image_id]
        setdir = 'set_%d' % image_id
        # read images
        g = getdata(os.path.join(self.out_dir,setdir,"img_g.fits"),memmap=False)
        r = getdata(os.path.join(self.out_dir,setdir,"img_r.fits"),memmap=False)
        z = getdata(os.path.join(self.out_dir,setdir,"img_z.fits"),memmap=False)

        image = np.zeros([info['height'], info['width'], 3], dtype=np.int16)

        # store raw image
        if store_raw:
            image_raw = np.zeros([info['height'], info['width'], 3], dtype=np.float64)
            
            image_raw[:,:,0] = z # red
            image_raw[:,:,1] = r # green
            image_raw[:,:,2] = g # blue
            self.raws[image_id] = image_raw

        I = (z+r+g)/3.0
        stretch = self.stretch
        Q = self.Q

        if normalize=='lupton':
            z = z*np.arcsinh(stretch*Q*(I))/(Q*I)
            r = r*np.arcsinh(stretch*Q*(I))/(Q*I)
            g = g*np.arcsinh(stretch*Q*(I))/(Q*I)
        elif normalize=='zscore':
            I = I*np.mean([np.std(g),np.std(r),np.std(z)])
            z = (z - np.mean(z))/I
            r = (r - np.mean(r))/I
            g = (g - np.mean(g))/I

        max_RGB = np.percentile([z,r,g], 99.95)
        # avoid saturation
        r = r/max_RGB; g = g/max_RGB; z = z/max_RGB

        # Rescale to 16-bit int
        int16_max = np.iinfo(np.int16).max
        r = r * int16_max
        g = g * int16_max
        z = z * int16_max

        image[:,:,0] = z # red
        image[:,:,1] = r # green
        image[:,:,2] = g # blue
        self.images[image_id] = image
        return image

    def load_mask(self, image_id):
        return self.masks[image_id], self.class_ids_mem[image_id]

    def load_mask_disk(self, image_id):
        # Load from disk
        info = self.image_info[image_id]
        # load image set via image_id from phosim output directory
        setdir = 'set_%d' % image_id
        maskdir = os.path.join(self.out_dir,setdir,"masks.fits")
        with fits.open(maskdir,memmap=False,lazy_load_hdus=False) as hdul:
            sources = len(hdul)
            print(hdul[0].header["CLASS_ID"])
            data = [hdu.data/np.max(hdu.data) for hdu in hdul]
            class_ids = [hdu.header["CLASS_ID"] for hdu in hdul]
        # make mask from threshold
        thresh = [0.005 if i == 1 else 0.08 for i in class_ids]
        masks = np.zeros([info['height'], info['width'], sources],dtype=np.uint8)
        for i in range(sources):
            """
            # inital guess
            x0, y0 = np.unravel_index(np.argmax(data[i]), masks.shape)
            sma = 10 # semi-major axis
            eps = 0 # ellipticity
            g = EllipseGeometry(x0, y0, sma, eps, pa)
            ellipse = Ellipse(data, geometry=g)
            isolist = ellipse.fit_image()
            # convert Petrosian isophot to mask
            position = [isolist.x0, isolist.y0]
            sma = isolist.sma
            b = sma*np.sqrt(1-isolist.eps**2)
            aper = EllipticalAperture(position, sma, b, isolist.pa)
            # create mask
            masks[:,:,i] = aper.to_mask(method='subpixel')
            """
            masks[:,:,i][data[i]>thresh[i]] = 1
            masks[:,:,i] = cv2.GaussianBlur(masks[:,:,i],(9,9),2)
        self.class_ids_mem[image_id] = np.array(class_ids,dtype=np.uint8)
        self.masks[image_id] = np.array(masks,dtype=np.bool)
        return self.masks[image_id], self.class_ids_mem[image_id]

def train(train_dir,val_dir):

    start_time = time.time()

    config = DESConfig()
    config.display()

    ## DATASET

    # Training dataset
    dataset_train = PhoSimDataset()
    dataset_train.load_sources(train_dir,dataset="training")
    dataset_train.prepare()

    # Validation dataset
    dataset_val = PhoSimDataset()
    dataset_val.load_sources(val_dir,dataset="validation")
    dataset_val.prepare()

    # Image augmentation
    augmentation = iaa.SomeOf((0, 4), [
        iaa.Fliplr(0.5),
        iaa.Flipud(0.5),
        iaa.OneOf([iaa.Affine(rotate=90),
                   iaa.Affine(rotate=180),
                   iaa.Affine(rotate=270)]),
        iaa.GaussianBlur(sigma=(0.0, np.random.random_sample()*4+2)),
        iaa.AddElementwise((-25, 25))
    ])

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
                augmentation=augmentation,
                epochs=15,
                layers='heads')

    # Fine tune all layers
    # Passing layers="all" trains all layers. You can also
    # pass a regular expression to select which layers to
    # train by name pattern.
    model.train(dataset_train, dataset_val,
                learning_rate=config.LEARNING_RATE / 10,
                augmentation=augmentation,
                epochs=25,
                layers="all")

    # Do one more with an even lower learning rate
    model.train(dataset_train, dataset_val,
                learning_rate=config.LEARNING_RATE / 100,
                augmentation=augmentation,
                epochs=35,
                layers="all")

    # Final stage
    model.train(dataset_train, dataset_val,
                learning_rate=config.LEARNING_RATE / 1000,
                augmentation=augmentation,
                epochs=50,
                layers="all")

    # Save weights
    model_path = os.path.join(MODEL_DIR, "astro_rcnn_decam.h5")
    model.keras_model.save_weights(model_path)

    print("Done in %.2f hours." % float((time.time() - start_time)/3600))

    return

def detect(directory, mode="detect", outdir = ".", normalize="zscore", plot_instances=False):

    print("Model in inference mode.")
    inference_config = InferenceConfig()

    # Use most recent weight file
    # Add code to download it if it does not exist
    model_path = os.path.join(MODEL_DIR, "astro_rcnn_decam.h5")

    # Recreate the model in inference mode
    model = modellib.MaskRCNN(mode="inference",
                          config=inference_config,
                          model_dir=MODEL_DIR)

    # Load trained weights
    print("Loading weights from ", model_path)
    model.load_weights(model_path, by_name=True)

    dataset = PhoSimDataset()

    # Assess performance
    if mode == "assess":

        # Load images
        dataset.load_sources(directory, dataset="validation")

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
            print(ps)
            print(rs)
            # galaxy
            APs,ps,rs = utils.compute_performance(dataset,model,inference_config,2,i)
            mean_APs_gal.append(APs)
            mean_ps_gal.append(ps)
            mean_rs_gal.append(rs)
            print(ps)
            print(rs)
        # Plot precision-recall
        visualize.plot_precision_recall_range(mean_APs_star,iou_thresholds,mean_ps_star,mean_rs_star,save_fig=True,title="star")
        visualize.plot_precision_recall_range(mean_APs_gal,iou_thresholds,mean_ps_gal,mean_rs_gal,save_fig=True,title="galaxy")

    # Detect
    else:

        # Load images
        dataset.load_sources(directory, dataset="test", normalize=normalize, store_raw=plot_instances)

        dataset.prepare()

        start_time = time.time()

        # Loop over batch of images (NOTE: assume batch size of one for now)
        results = []
        # Loop over images in batch
        for image_id in range(len(dataset.image_info)):
            # Load image and run detection
            image = dataset.load_image(image_id)
            # Detect objects
            r = np.array(model.detect([image],verbose=0))
            results.append(r[0])
            # Visualize as it steps through
            if plot_instances:
                image_raw = dataset.load_image(image_id, raw=True)
                im_disp = make_lupton_rgb(image_raw[:,:,0], image_raw[:,:,1], image_raw[:,:,2], minimum=np.percentile(image_raw, 50), stretch=dataset.stretch, Q=dataset.Q)
                visualize.display_instances(im_disp, r[0]['rois'], r[0]['masks'], r[0]['class_ids'], dataset.class_names, r[0]['scores'],save_fig=True)

        print("Detected %d images in %.2f seconds with batch size of 1." % (len(dataset.image_info), float(time.time() - start_time)))

        # save masks as fits file
        for j,r in enumerate(results):
            hdul = fits.HDUList()
            for i,mask in enumerate(r["class_ids"]):
                hdr = fits.Header()
                hdr["BITPIX"] = 8
                hdr["CLASS_ID"] = r["class_ids"][i]
                hdr["SCORE"] = round(r["scores"][i],3)
                hdr["BBOX"] = str(r["rois"][i])
                x0 = r["rois"][i][1]
                y0 = r["rois"][i][0]
                x1 = r["rois"][i][3]
                y1 = r["rois"][i][2]
                hdr["WEIGHTS"] = os.path.basename(model_path)
                mask_i = r["masks"][y0:y1,x0:x1,i].astype(dtype=np.uint8)
                hdul.append(fits.ImageHDU(mask_i,header=hdr))

            print("Writing to output_%d.fits" % j)
            hdul.writeto(os.path.join(outdir, ("output_%d.fits" % j)) ,overwrite=True)

        print("Success!")

    return

if __name__ == "__main__":
    import argparse
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Mask R-CNN for star/galaxy detection, classification, and deblending.')
    parser.add_argument("command",metavar="<command>",help="'train', 'detect', or 'assess'")
    parser.add_argument("datapath",metavar="<datapath>",default="none",help="path to set of FITS images e.g. 'example' example directory")
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--normalize", default="zscore")
    args = parser.parse_args()
    datapath = os.path.abspath(args.datapath.split(",")[0])

    # Train or evaluate
    if args.command == "train":
        validationpath = os.path.abspath(args.datapath.split(",")[1])
        train(datapath,validationpath)
    elif args.command == "detect":
        detect(datapath, outdir = args.outdir, normalize=args.normalize)
    elif args.command == "assess":
        detect(datapath,mode="assess")
    else:
        print("'{}' is not recognized. "
              "Use 'train', 'detect', or 'detect_assess'".format(args.command))

