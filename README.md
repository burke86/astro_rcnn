# Astro R-CNN

Detect, classify, and deblend sources in astronomical images using [Mask R-CNN](https://github.com/matterport/Mask_RCNN).

*Authors:* 
[Colin J. Burke](https://astro.illinois.edu/directory/profile/colinjb2), [Patrick D. Aleo](https://astro.illinois.edu/directory/profile/paleo2), [Yu-Ching Chen](https://astro.illinois.edu/directory/profile/ycchen)

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
The result will be a multi-extension FITS file ```output_0.fits``` with a segmentation mask cutout in each exension corresponding to an object detection (extension number=SOURCE_ID). A table with the CLASS_ID (star=1,galaxy=2), bounding box (BBOX: y1,x1,y2,x2), and detection confidence (SCORE) are in the header.

![infrence](https://user-images.githubusercontent.com/13906989/61251399-f3588400-a71f-11e9-896d-e73008a4e0e3.png)
Example of Astro R-CNN detection on a real DECam image. See [demo.ipynb](https://github.com/burke86/deblend_maskrcnn/blob/master/demo.ipynb) for an interactive demonstration, including how to train on your own images. 

<img src="https://user-images.githubusercontent.com/13906989/61023273-e1b55c00-a36e-11e9-85df-cf7471a44aa9.png" alt="deblending" width="512"/>

Examples of Astro R-CNN deblending on a real DECam image.

This is a simple repository intended for demonstration purposes. For use with full-scale images or surveys, please contact the authors.

## Cite

```
@article{Burke2019,
  title={Classifying and Deblending Astronomical Sources with Mask R-CNN Deep Learning},
  author={Burke et al (in prep.)},
  year={2019},
  journal={MNRAS},
}
```
