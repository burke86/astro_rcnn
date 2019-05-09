///
/// @package phosim
/// @file sed.cpp
/// @brief sed reader
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

void Observation::readSed(const std::string & filename, int mode) {

    long lsedptr;
    char *sptr;
    char *sptr2;
    char *sptr3;
    char tempstring[4096];
    double tempf1;
    double tempf2;
    double dw;
    double cdw;
    long closestw = 0;
    long j;
    FILE *indafile;
    gzFile ingzfile;
    char line[4096];
    char tempstring1[512];
    char tempstring2[512];
    double monochromaticFactor;
    double standardwavelength = 500.0; // nm 

    normwave = standardwavelength;
    monochromaticFactor = 1.0;

    if (sedptr == 0) {
        nreallocs = 0;
        sedMax = 4194304; // 4M elements, 32MB memory
        sedW = static_cast<double*>(malloc((long)(sedMax*sizeof(double))));
        sedC = static_cast<double*>(malloc((long)(sedMax*sizeof(double))));
    } else {
        if (sedptr > (sedMax - 100000)) {
            nreallocs++;
            sedMax = 2*sedMax;
            sedW = static_cast<double*>(realloc(sedW, (long)((sedMax)*sizeof(double))));
            sedC = static_cast<double*>(realloc(sedC, (long)((sedMax)*sizeof(double))));
        }
    }

    lsedptr = 0;
    if (mode == 0) sprintf(tempstring, "%s/%s", seddir.c_str(), filename.c_str());
    if (mode == 1) sprintf(tempstring, "%s", filename.c_str());
    if (flatdir == 1) {
        sptr = strtok_r(tempstring, "/", &sptr2);
        do {
            sptr3 = sptr;
            sptr = strtok_r(NULL, "/", &sptr2);
        } while (sptr != NULL);
        sprintf(tempstring, "%s", sptr3);
    }

    if (strstr(tempstring, "laser") != NULL) {
        normwave = domewave;
        closestw = 0;
        cdw = 1e30;

        tempf1 = 0.0;
        tempf2 = 0.0;
        for (int k = 0; k < 5; k++) {

            if (k == 0) {
                tempf1 = minwavelength;
                tempf2 = 0.0;
            }
            if (k == 1) {
                tempf1 = domewave - 0.5e-6;
                tempf2 = 0.0;
            }
            if (k == 2) {
                tempf1 = domewave;
                tempf2 = 1.0/tempf1;
            }
            if (k == 3) {
                tempf1 = domewave + 0.5e-6;
                tempf2 = 0.0;
            }
            if (k == 4) {
                tempf1 = maxwavelength;
                tempf2 = 0.0;
            }
			
            sedW[sedptr + lsedptr] = tempf1;
            sedC[sedptr + lsedptr] = tempf2*tempf1;

            dw = fabs(tempf1 - domewave);

            if (dw < cdw) {
                cdw = dw;
                closestw = lsedptr;
            }
            lsedptr++;
        }
        monochromaticFactor = normwave/0.5e-6*(0.01*H_CGS/pow(10.0, (20.0+48.6)/(-2.5)));

    } else {


        if (strstr(tempstring, ".gz") == NULL) {

            indafile = fopen(tempstring, "r");
            if (indafile == NULL) {
                fprintf(stderr, "Can't find SED file: %s\n", tempstring);
                exit(1);
            }

            closestw = 0;
            cdw = 1e30;
            while (fgets(line, 4096, indafile)) {
                sscanf(line, "%s %s", tempstring1, tempstring2);
                tempf1 = strtod(tempstring1, NULL);
                tempf2 = strtod(tempstring2, NULL);
                sedW[sedptr + lsedptr] = tempf1;
                sedC[sedptr + lsedptr] = tempf2*tempf1;
                if (tempstring1[0] != '#' && (tempstring2[0] == 'n' || tempstring2[0] == 'N')) {
                    fprintf(stderr, "Warning:   SED file: %s contains a NaN!\n", tempstring);
                    sedC[sedptr + lsedptr] = 0.0;
                }

                dw = fabs(tempf1 - standardwavelength);

                if (dw < cdw && tempf2 > 0.0) {
                    cdw = dw;
                    closestw = lsedptr;
                }
                lsedptr++;
                if (lsedptr >= 100000) {
                    fprintf(stderr, "Error:  Too many lines in SED file: %s\n", tempstring);
                    exit(1);
                }
            }
            fclose(indafile);

        } else {

            ingzfile = gzopen(tempstring, "r");
            if (ingzfile == NULL) {
                fprintf(stderr, "Can't find SED file: %s\n", tempstring);
                exit(1);
            }

            closestw = 0;
            cdw = 1e30;
            while (gzgets(ingzfile, line, 4096)) {
                sscanf(line, "%s %s", tempstring1, tempstring2);
                tempf1 = strtod(tempstring1, NULL);
                tempf2 = strtod(tempstring2, NULL);

                sedW[sedptr + lsedptr] = tempf1;
                sedC[sedptr + lsedptr] = tempf2*tempf1;
                if (tempstring1[0] != '#' && (tempstring2[0] == 'n' || tempstring2[0] == 'N')) {
                    fprintf(stderr, "Warning:   SED file: %s contains a NaN!\n", tempstring);
                    sedC[sedptr + lsedptr] = 0.0;
                }

                dw = fabs(tempf1 - standardwavelength);

                if (dw < cdw && tempf2 > 0.0) {
                    cdw = dw;
                    closestw = lsedptr;
                }
                lsedptr = lsedptr + 1;
                if (lsedptr >= 100000) {
                    fprintf(stderr, "Error:  Too many lines in SED file: %s\n", tempstring);
                    exit(1);
                }
            }
            gzclose(ingzfile);

        }

        // monochromatic exception
        if (lsedptr == 1) {
            lsedptr--;
            normwave = tempf1;
            closestw = 0;
            cdw = 1e30;

            tempf1 = 0.0;
            tempf2 = 0.0;
            for (int k = 0; k < 5; k++) {

                if (k == 0) {
                    tempf1 = minwavelength;
                    tempf2 = 0.0;
                }
                if (k == 1) {
                    tempf1 = normwave - 0.5e-6;
                    tempf2 = 0.0;
                }
                if (k == 2) {
                    tempf1 = normwave;
                    tempf2 = 1.0/normwave;
                }
                if (k == 3) {
                    tempf1 = normwave + 0.5e-6;
                    tempf2 = 0.0;
                }
                if (k == 4) {
                    tempf1 = maxwavelength;
                    tempf2 = 0.0;
                }

                sedW[sedptr + lsedptr] = tempf1;
                sedC[sedptr + lsedptr] = tempf2*tempf1;

                dw = fabs(tempf1 - normwave);

                if (dw < cdw) {
                    cdw = dw;
                    closestw = lsedptr;
                }
                lsedptr++;
            }
            monochromaticFactor = normwave/0.5e-6*(0.01*H_CGS/pow(10.0, (20.0+48.6)/(-2.5)));
        }


    }

    for (j = 0; j < lsedptr; j++) {
        if (j != 0 && j != (lsedptr - 1)) sedC[sedptr + j] = sedC[sedptr + j]*(sedW[sedptr + j + 1] - sedW[sedptr + j - 1])/2.0;
        if (j == 0) sedC[sedptr + j] = sedC[sedptr + j]*(sedW[sedptr + j + 1] - sedW[sedptr + j]);
        if (j == (lsedptr - 1)) sedC[sedptr + j] = sedC[sedptr + j]*(sedW[sedptr + j] - sedW[sedptr + j - 1]);
    }

    tempf1 = 0;
    for (j = 0; j < lsedptr; j++) {
        tempf1 += sedC[sedptr + j];
    }
    for (j = 0; j < lsedptr; j++) {
        sedC[sedptr + j] = sedC[sedptr + j]/tempf1;
    }

    if (closestw == 0) {
        sedDwdp[nsedptr] = (sedW[sedptr + closestw + 1] - sedW[sedptr + closestw])/sedC[sedptr + closestw];
    } else {
        sedDwdp[nsedptr] = (sedW[sedptr + closestw + 1] - sedW[sedptr + closestw - 1])/sedC[sedptr + closestw]/2.0;
    }

    if (sedC[sedptr + closestw] <= 0.0) {
        printf("Error in SED file; 0 value at normalizing wavelength\n");
        sedDwdp[nsedptr] = 0.0;
    }
    sedDwdp[nsedptr] = sedDwdp[nsedptr]*monochromaticFactor;

    for (j = 1; j < lsedptr; j++) {
        sedC[sedptr + j] += sedC[sedptr + j - 1];
    }

    sedN[nsedptr] = lsedptr;
    sedPtr[nsedptr] = sedptr;
    sedptr += lsedptr;
    nsedptr++;

    if (nsedptr >= 10000) {
        printf("Error:   Too many SED files\n");
        exit(1);
    }

};
