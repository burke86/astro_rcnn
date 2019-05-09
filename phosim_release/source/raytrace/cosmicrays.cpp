///
/// @package phosim
/// @file cosmicrays.cpp
/// @brief main for raytrace code
///
/// @brief Created by:
/// @author Mallory Young (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///


// Name:
//      create_rays
//  Purpose:
//      Randomly puts cosmic rays onto entire focal plane one chip at a time.
//  Inputs:
//      needs the focalplane fits files
//  Outputs:
//      fits files with superimposed cosmic rays
//  Written by:
//      Mallory Young. NOVEMBER 2009
//  For:
//      John Peterson. Purdue University
//  Method:
//     (3) determines an estimated number of rays per chip according to the input values for number of rays per
//      cm^2 per second, pixel size, and exposure time.  The established ray density according to Kirk
//      Gilmore is .6/cm^2/15 seconds.  The number of rays per chip is adjusted accordingly.
//      (4) The probability per pixel of having a ray present is calculated by dividing the estimated number of rays
//          by the total number of pixels.
//     (5)    loops through each pixel in the array.  For each pixel a random number between 0 and 1 is generated; if the
//           number falls below the probability, then a ray is placed there.
//     (6) randomly transforms ray image to one of 8 possible orientations formed by 90, 180, 270 degree rotation and
//         vertical flip around horizontal axis.
//     (7)    The transformed cosmic ray array to be superimposed onto the existing fits file is randomly chosen between the
//          130 ray text files.  Each ray file contains a different sized ray.
//     (8)    The pixel values of the ray file are added to the pixel values of the existing fits file, thus adding the ray,
//          until the desired number of rays are added.

int *transformArray(int* matrix, int cols, int rows, int ro, int f);

void Image::cosmicRays(long long *raynumber) {

    float numrays = 0.0;
    /* (3)  calculates estimated number of rays to be superimposed upon the fits file image and the
       probability per pixel of there being a ray there*/
    numrays = pixsize*pow(10 , -6);
    numrays = pow(numrays, 2)*(exptime/15.0)*pow(10, 4)*raydensity*pixelsx*pixelsy;
    double probability = 0.0;
    probability = (double)(numrays/pixelsx/pixelsy);
    int i = 0 , l = 0;
    double nm = 0.0;
    int filenum = 0;
    float rnd = 0.0;
    long xPos, yPos;

    for (i = 0; i < pixelsx; i++) {
        for (l = 0; l < pixelsy; l++) {
            nm = random[0].uniform();
            if (nm < probability)    {

                rnd = random[0].uniform();
                if (rnd != 1.0) {
                    filenum = 1 + (int)(rnd*130.0);/*chooses random ray file*/
                } else {
                    filenum = 130;
                }
                char str[4];
                char rayfile[4096]; /*indicate directory/location of ray files*/
                if (flatdir == 0) sprintf(rayfile , "%s/cosmic_rays/iray" , datadir.c_str());
                if (flatdir == 1) sprintf(rayfile , "./iray");
                sprintf(str , "%i" , filenum);
                strcat(rayfile , str);
                strcat(rayfile , ".txt");
                int rows = 0, cols = 0;
                FILE *f = fopen(rayfile , "r");
                int k = 0 , j = 0;
                fscanf(f, "%i %i" , &cols , &rows);/*reads in the number of columns and rows in the ray array*/
                /*makes sure the ray array will fit in the fits array; if so, it adds the ray*/
                int raywidth = rows;
                if (cols > rows) raywidth = cols;
                *raynumber = *raynumber + 1;
                int data[cols][rows];
                int temp[500];
                int z = 0;

                /* loop through 2D ray array & store the values in a temp 1D array */
                while (!feof(f)) {
                    fscanf(f, "%d", &temp[z]);
                    z++;
                }
                fclose(f);

                /* go thru temp 1D array and put values into a 2D array */
                for (j = 0; j < rows; j++) {
                    for (k = 0; k < cols; k++) {
                        data[k][j] = temp[cols*j + k];
                    }
                }
                int fl = 0 , ro = 0;
                fl = (int) (random[0].uniform()*2);
                ro = (int) (random[0].uniform()*4);
                int rows2 = 0, cols2 = 0;
                if ((ro == 1) || (ro == 3)) {
                    rows2 = cols;
                    cols2 = rows;
                } else {
                    rows2 = rows;
                    cols2 = cols;
                }
                int *newmatrix = NULL;
                int value;
                newmatrix = static_cast<int *>(malloc(sizeof(int)*cols*rows));
                newmatrix = transformArray(*data, cols, rows , ro , fl );

                /* adding in the ray array to the fits file array at a random location
                   determined by corner_pix */
                for (k = 0; k < rows2; k++) {
                    for (j = 0; j < cols2; j++) {
                        xPos = i + k + minx;
                        yPos = l + j + miny;
                        if (xPos >= minx && xPos <= maxx && yPos >= miny && yPos <= maxy) {
                            value = newmatrix[j*rows2 + k];
                            if (value < 0) value = 0;
                            *(state.focal_plane + chip.nampx*(yPos - miny) + (xPos - minx)) += static_cast<unsigned long>(scalenumber*value);
                            if (*(state.focal_plane + chip.nampx*(yPos - miny) + (xPos - minx)) > well_depth)
                                *(state.focal_plane + chip.nampx*(yPos - miny) + (xPos - minx)) = well_depth;
                        }

                    }
                }

            }
        }
    }

}

int *transformArray(int* matrix, int cols, int rows, int ro, int fl) {

    int *newmatrix = NULL;
    int *tempmatrix = NULL;
    tempmatrix = static_cast<int *>(malloc(sizeof(int)*cols*rows));
    newmatrix = static_cast<int *>(malloc(sizeof(int)*cols*rows));
    int x = 0, y = 0, i = 0;
    //flip vertically
    if (fl == 0) {
        for (y = 0; y < cols; y++) {
            for (x = 0; x < rows; x++) {
                tempmatrix[rows - x - 1 + rows*y] = matrix[rows*y + x];
            }
        }
    }
    //don't flip
    if (fl == 1) {
        for (y = 0; y < rows; y++) {
            for (x = 0; x < cols; x++) {
                tempmatrix[cols*y + x] = matrix[y*cols + x];
            }
        }
    }
    //rotate 90 degrees
    if (ro == 1) {
        for (y = 0; y < cols; y++) {
            for (x = 0; x < rows; x++) {
                newmatrix[y + cols*(rows - 1 - x)] = tempmatrix[i];
                i++;
            }
        }
    }
    //rotate 180 degrees
    if (ro == 2) {
        for (x = 0; x < cols*rows; x++) {
            newmatrix[rows*cols - 1 - x] = tempmatrix[x];
        }
    }
    //rotate 270 degrees
    if (ro == 3) {
        int i = 0;
        for (y = 0; y < cols; y++) {
            for (x = 0; x < rows; x++) {
                newmatrix[(x*cols) + (cols - 1) - y] = tempmatrix[i];
                i++;
            }
        }
    }
    if (ro == 0) {
        for (x = 0; x < cols*rows; x++) {
            newmatrix[x] = tempmatrix[x];
        }
    }
    return newmatrix;
}
