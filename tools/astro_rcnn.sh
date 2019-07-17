#!/bin/sh
cd ..
module load powerai
python astro_rcnn.py train trainingset,validationset
