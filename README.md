# Astro R-CNN

Detect, classify, and deblend sources in astronomical images using [Mask R-CNN](https://github.com/matterport/Mask_RCNN).

*Reference Paper:* [Burke et al. (submitted), MNRAS.](http://adsabs.harvard.edu/doi/10.1093/mnras/stz2845)

*Corresponding Author:* 
[Colin J. Burke](https://astro.illinois.edu/directory/profile/colinjb2), University of Illinois at Urbana-Champaign

*Additional Repository Authors:* Patrick D. Aleo, Yu-Ching Chen, Joshua Yao-Yu Lin.

## Description:

Astro R-CNN is a deep learning method for efficiently performing all tasks of source detection, classification, and deblending on astronomical images.

Setup:
```
pip install -r requirements.txt
```

Usage:
```
./astro_rcnn detect example
```
This will run the model in inference mode with pre-trained DECam weights (use GPU for best performance). The result will be a multi-extension FITS file ```output_0.fits``` with a segmentation mask cutout in each extension corresponding to an object detection (extension number=SOURCE_ID). In each header, you will find the CLASS_ID (star=1,galaxy=2), bounding box (BBOX: y1,x1,y2,x2), and detection confidence (SCORE).

![infrence](https://user-images.githubusercontent.com/13906989/61251399-f3588400-a71f-11e9-896d-e73008a4e0e3.png)
Example of Astro R-CNN detection on a real DECam image. See [demo.ipynb](https://github.com/burke86/deblend_maskrcnn/blob/master/demo.ipynb) for an interactive demonstration, including how to train on your own images. 

<img src="https://user-images.githubusercontent.com/13906989/61023273-e1b55c00-a36e-11e9-85df-cf7471a44aa9.png" alt="deblending" width="512"/>

Examples of Astro R-CNN deblending on a real DECam image.

This is a simple repository intended for demonstration purposes. For use with full-scale images or surveys, please contact the authors.

## Training:

To train your own model, first download PhoSim training data (or [make your own](https://bitbucket.org/phosim/phosim_release)) into the project root directory: [training set (1,000 images)](https://uofi.box.com/s/svlkblkh5o4a3q3qwu7iks6r21cmmu64) [validation set (250 images)](https://uofi.box.com/s/m22q747nawtxq8e5iihjulpapwlvucr5).

Then, try:
```
./astro_rcnn train trainingset,validationset
```
Depending on your setup, you should adjust the configuration settings and decide which weights to initialize with in ```astro_rcnn.py```.

If you would like a simulated test dataset beyond ```example``` (1 image) to assess the network's performance: [test set (50 images)](https://uofi.box.com/s/bmtkjrj9g832w9qybjd1yc4l6cyqx6cs).

```
./astro_rcnn assess testset
```
This will generate mean AP score plots for stars and galaxies in ```testset```.

Also available are real DECam datasets of clusters of galaxies: [ACO 1689 (50 images)](https://uofi.box.com/s/7cy1yuahmaiucq857wgo3exln8wvc825) [ACO 3266 (50 images)](https://uofi.box.com/s/jaug7rjvpkotqevfnmm3pajo7011f48v).
