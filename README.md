## Astro R-CNN

Detect, classify, and deblend sources in astronomical images using [Mask R-CNN](https://github.com/matterport/Mask_RCNN).

*Authors:* 
[Colin J. Burke](https://astro.illinois.edu/directory/profile/colinjb2), [Patrick D. Aleo](https://astro.illinois.edu/directory/profile/paleo2), [Yu-Ching Chen](https://astro.illinois.edu/directory/profile/ycchen)

# Description:

Astro R-CNN is a deep learning method for efficiently performing all tasks of source detection, classification, and deblending on astronomical images.

Usage:
```
astro_rcnn detect example.fits
```
The result will be a FITS file ```output.fits``` with the segmentation mask in the image exension. A table with the mask_id, class_id (star=1, galaxy=2), and confidence are in the header.

![infrence](https://user-images.githubusercontent.com/13906989/61251399-f3588400-a71f-11e9-896d-e73008a4e0e3.png)
Example of Astro R-CNN detection on a real DECam image. See [train.ipynb](https://github.com/burke86/deblend_maskrcnn/blob/master/train.ipynb) for an interactive demonstration. 

<img src="https://user-images.githubusercontent.com/13906989/61023273-e1b55c00-a36e-11e9-85df-cf7471a44aa9.png" alt="deblending" width="512"/>

Examples of Astro R-CNN deblending on a real DECam image.

This is a simple repository intended for demonstration purposes. For use with full-scale images or surveys, please contact the authors.

# Cite

```
@article{Burke2019,
  title={Classifying and Deblending Astronomical Sources with Mask R-CNN Deep Learning},
  author={Burke et al (in prep.)},
  year={2019},
  journal={MNRAS},
}
```
