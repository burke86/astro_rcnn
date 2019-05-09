///
/// @package phosim
/// @file air.cpp
/// @brief air class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Emily Grace (Purdue)
/// @author En-Hsin Peng (Purdue)
/// @author Michael Wood-Vasey (Pitt)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <math.h>
#include <vector>

#include "ancillary/readtext.h"
#include "air.h"
#include "constants.h"
#include "helpers.h"

void Air::opacitySetup(double zenith, double moonalt, std::vector<double> height, double groundlevel, double raynorm, double o2norm, double h2onorm,
                       double o3norm, double aerosoltau, double aerosolindex, long layers, std::string dir, double *airmass) {

    std::vector<double> rayprofile;
    std::vector<double> rayprofile_wavelength;
    std::vector<double> nprofile;
    std::vector<double> nprofile_height;
    std::vector<double> o3profile;
    std::vector<double> o3profile_height;
    std::vector<double> o3cs;
    std::vector<double> o3cs_wavelength;
    std::vector<std::vector<double> > o2cs;
    std::vector<double> o2cs_wavelength;
    std::vector<std::vector<double> > h2ocs;
    std::vector<double> h2ocs_wavelength;
    std::vector<double> h2oprofile;
    std::vector<double> h2oprofile_height;
    std::vector<double> airmass_height;
    std::vector<double> airmass_heightMoon;

    double chi, minchi, cheight, cheightnext, prheight;
    int closestlayer;
    long index, index2, index3, index4;
    double rindex;
    double transmission, transmission2, transmission3, transmission4, transmission5;
    double transmissionMoon, transmissionMoon2, transmissionMoon3, transmissionMoon4, transmissionMoon5;
    double wavelength;
    double moonZenith;

    moonZenith = PI/2 - moonalt;

    readText::readCol(dir + "/rayprofile.txt", rayprofile_wavelength, rayprofile);
    readText::readCol(dir + "/nprofile.txt", nprofile_height, nprofile);
    readText::readCol(dir + "/o3cs.txt", o3cs_wavelength, o3cs);
    readText::readCol(dir + "/o3profile.txt", o3profile_height, o3profile);
    readText::readMultiCol(dir + "/o2cs.txt", o2cs_wavelength, o2cs);
    readText::readMultiCol(dir + "/h2ocs.txt", h2ocs_wavelength, h2ocs);
    readText::readCol(dir + "/h2oprofile.txt", h2oprofile_height, h2oprofile);

    airmassLayer = new double[layers + 1]();
    airmassLayerMoon = new double[layers + 1]();
    tauWavelength = new double[180001]();
    tau = new double*[5*(layers + 1)]();
    tauMoon = new double*[5*(layers + 1)]();
    for (long i = 0; i < (5*(layers + 1)); i++) {
        tau[i] = new double[180001]();
        tauMoon[i] = new double[180001]();
    }
    for (long i = 0; i < 180001; i++) {
        tauWavelength[i] = 0.3 + 0.9*((static_cast<double>(i))/180000.0);
    }

    prheight = 70.0; // km

    // Airmass calculations

    // approximate airmass (for fits header)
    *airmass = 1.0/ (cos(zenith) + 0.025*exp(-11.0*cos(zenith)));

    // Approximation for ozone at 20km for Radius_Earth = 6378km
    // See, e.g., B. E. Schaefer (1993) VA, 35, 4, 311
    // airmass_ozone=1.0/ sqrt(1.0 - pow(sin(zenith) / (1+ (20.0/6378.0)),2));
    //  But instead can do even better by calculating this at each height as part of the tabulation below
    // Aerosol airmass
    // Eq. 3a from Schaefer 1993
    // Gives the airmass for a component with a given scale height in the atmosphere.
    // Schaefer93 Suggested 1.5 km for aerosols.  Referred to Hayes and Latham 1975.
    // airmass_aerosol=1./(cos(zenith)+0.01*sqrt(1.5)*exp((-30*cos(zenith))/sqrt(1.5)));
    // Previous was simple treatment of overall average airmass
    // But we can do better since we have the height information directly
    // And, in fact, this should actually be used with all of the absorption terms, not just ozone
    // This is the airmass term for a given height
    // but do we actually need the differential instead
    // Let's calculate this outside of the wavelength loop to save 180,000 calculations.

    airmass_height.resize(401, 0);
    airmass_heightMoon.resize(401, 0);
    for (long k = 0; k < 401; k++) {
        cheight = pow(10.0, 1.845 - (static_cast<double>(k))/100.0);
        cheightnext = pow(10.0, 1.845 - (static_cast<double>(k) + 1)/100.0);
        if (zenith == 0.0) {
            airmass_height[k] = 1.0;
        } else {
            airmass_height[k] = RADIUS_EARTH*(sqrt(pow(1 + cheight/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1)-
                                              sqrt(pow(1 + cheightnext/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1))*
                sin(zenith)/(cheight - cheightnext);
        }
        if (moonZenith == 0.0) {
            airmass_heightMoon[k] = 1.0;
        } else {
            airmass_heightMoon[k] = RADIUS_EARTH*(sqrt(pow(1 + cheight/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1)-
                                              sqrt(pow(1 + cheightnext/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1))*
                sin(moonZenith)/(cheight - cheightnext);
        }
    }

    for (long k = 0; k < layers + 1; k++) {
        if (zenith == 0.0) {
            airmassLayer[k] = 1.0;
        } else {
            if (k == 0) {
                airmassLayer[k] = RADIUS_EARTH*(sqrt(pow(1 + 70.0/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1)-
                                                sqrt(pow(1 + height[k]/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1))*
                    sin(zenith)/(70.0 - height[k]);
            } else {
                airmassLayer[k] = RADIUS_EARTH*(sqrt(pow(1 + height[k - 1]/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1)-
                                                sqrt(pow(1 + height[k]/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1))*
                    sin(zenith)/(height[k - 1] - height[k]);
            }
        }
        if (moonZenith == 0.0) {
            airmassLayerMoon[k] = 1.0;
        } else {
            if (k == 0) {
               airmassLayerMoon[k] = RADIUS_EARTH*(sqrt(pow(1 + 70.0/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1)-
                                                sqrt(pow(1 + height[k]/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1))*
                    sin(moonZenith)/(70.0 - height[k]);
            } else {
               airmassLayerMoon[k] = RADIUS_EARTH*(sqrt(pow(1 + height[k - 1]/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1)-
                                                sqrt(pow(1 + height[k]/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1))*
                    sin(moonZenith)/(height[k - 1] - height[k]);
            }
        }
    }

    for (long i = 0; i < 180001; i++) {
        for (long k = 0; k < 401; k++) {
            wavelength = tauWavelength[i];
            cheight = pow(10.0, 1.845 - (static_cast<double>(k))/100.0);

            transmission = 0.0;
            transmission2 = 0.0;
            transmission3 = 0.0;
            transmission4 = 0.0;
            transmission5 = 0.0;

            // rayleigh scattering
            index = find_linear_v(rayprofile_wavelength, wavelength, &rindex);
            index3 = find_linear_v(nprofile_height, cheight*1000.0 + groundlevel, &rindex);
            transmission += raynorm*airmass_height[k]*rayprofile[index]*nprofile[index3]*(prheight - cheight)*1e5;

            // O3 absorption
            index = find_linear_v(o3cs_wavelength, wavelength, &rindex);
            index2 = find_linear_v(o3profile_height, cheight*1000.0 + groundlevel, &rindex);
            transmission2 += o3norm*airmass_height[k]*o3cs[index]*o3profile[index2]*(prheight - cheight)*1e5;


            // O2 absorption
            index = find_linear_v(o2cs_wavelength, wavelength, &rindex);
            index2 = ((long)((cheight + groundlevel/1000.0)/7.0));
            if (index2 <= 8) {
                if (index2 < 8) {
                    double tempf1 = ((groundlevel/1000.0 + cheight) - index2*7.0)/7.0;
                    transmission3 += o2norm*0.2096*airmass_height[k]*(1 - tempf1)*o2cs[index][index2]*nprofile[index3]*(prheight - cheight)*1e5;
                    transmission3 += o2norm*0.2096*airmass_height[k]*tempf1*o2cs[index][index2 + 1]*nprofile[index3]*(prheight - cheight)*1e5;
                    // H2O absorption
                    index = find_linear_v(h2ocs_wavelength, wavelength, &rindex);
                    index4 = find_linear_v(h2oprofile_height, cheight*1000.0 + groundlevel, &rindex);
                    transmission4 += h2onorm*h2oprofile[index4]*airmass_height[k]*(1 - tempf1)*
                        h2ocs[index][index2]*nprofile[index3]*(prheight - cheight)*1e5;
                    transmission4 += h2onorm*h2oprofile[index4]*airmass_height[k]*tempf1*
                        h2ocs[index][index2 + 1]*nprofile[index3]*(prheight - cheight)*1e5;
                } else {
                    transmission3 += o2norm*0.2096*airmass_height[k]*o2cs[index][index2]*nprofile[index3]*(prheight - cheight)*1e5;
                    // H2O absorption
                    index = find_linear_v(h2ocs_wavelength, wavelength, &rindex);
                    index4 = find_linear_v(h2oprofile_height, cheight*1000.0 + groundlevel, &rindex);
                    transmission4 += h2onorm*h2oprofile[index4]*airmass_height[k]*h2ocs[index][index2]*nprofile[index3]*(prheight - cheight)*1e5;
                }
            }

            // putting at 1.5 km (hack)
            if (k == 167) transmission5 += aerosoltau*airmass_height[k]*pow(wavelength/0.5, aerosolindex);


            // Moon

            transmissionMoon = 0.0;
            transmissionMoon2 = 0.0;
            transmissionMoon3 = 0.0;
            transmissionMoon4 = 0.0;
            transmissionMoon5 = 0.0;

            // rayleigh scattering
            index = find_linear_v(rayprofile_wavelength, wavelength, &rindex);
            index3 = find_linear_v(nprofile_height, cheight*1000.0 + groundlevel, &rindex);
            transmissionMoon += raynorm*airmass_heightMoon[k]*rayprofile[index]*nprofile[index3]*(prheight - cheight)*1e5;

            // O3 absorption
            index = find_linear_v(o3cs_wavelength, wavelength, &rindex);
            index2 = find_linear_v(o3profile_height, cheight*1000.0 + groundlevel, &rindex);
            transmissionMoon2 += o3norm*airmass_heightMoon[k]*o3cs[index]*o3profile[index2]*(prheight - cheight)*1e5;


            // O2 absorption
            index = find_linear_v(o2cs_wavelength, wavelength, &rindex);
            index2 = ((long)((cheight + groundlevel/1000.0)/7.0));
            if (index2 <= 8) {
                if (index2 < 8) {
                    double tempf1 = ((groundlevel/1000.0 + cheight) - index2*7.0)/7.0;
                    transmissionMoon3 += o2norm*0.2096*airmass_heightMoon[k]*(1 - tempf1)*o2cs[index][index2]*nprofile[index3]*(prheight - cheight)*1e5;
                    transmissionMoon3 += o2norm*0.2096*airmass_heightMoon[k]*tempf1*o2cs[index][index2 + 1]*nprofile[index3]*(prheight - cheight)*1e5;
                    // H2O absorption
                    index = find_linear_v(h2ocs_wavelength, wavelength, &rindex);
                    index4 = find_linear_v(h2oprofile_height, cheight*1000.0 + groundlevel, &rindex);
                    transmissionMoon4 += h2onorm*h2oprofile[index4]*airmass_heightMoon[k]*(1 - tempf1)*
                        h2ocs[index][index2]*nprofile[index3]*(prheight - cheight)*1e5;
                    transmissionMoon4 += h2onorm*h2oprofile[index4]*airmass_heightMoon[k]*tempf1*
                        h2ocs[index][index2 + 1]*nprofile[index3]*(prheight - cheight)*1e5;
                } else {
                    transmissionMoon3 += o2norm*0.2096*airmass_heightMoon[k]*o2cs[index][index2]*nprofile[index3]*(prheight - cheight)*1e5;
                    // H2O absorption
                    index = find_linear_v(h2ocs_wavelength, wavelength, &rindex);
                    index4 = find_linear_v(h2oprofile_height, cheight*1000.0 + groundlevel, &rindex);
                    transmissionMoon4 += h2onorm*h2oprofile[index4]*airmass_heightMoon[k]*h2ocs[index][index2]*nprofile[index3]*(prheight - cheight)*1e5;
                }
            }

            // putting at 1.5 km (hack)
            if (k == 167) transmissionMoon5 += aerosoltau*airmass_heightMoon[k]*pow(wavelength/0.5, aerosolindex);


            minchi = fabs(log10(70.0) - log10(cheight));
            closestlayer = -1;
            for (long l = 0; l < layers; l++) {
                chi = fabs(log10(height[l]) - log10(cheight));
                if (chi < minchi) {
                    closestlayer = l;
                    minchi = chi;
                }
            }

            *(tau[5*(closestlayer + 1) + 0] + i) += transmission;
            *(tau[5*(closestlayer + 1) + 1] + i) += transmission2;
            *(tau[5*(closestlayer + 1) + 2] + i) += transmission3;
            *(tau[5*(closestlayer + 1) + 3] + i) += transmission4;
            *(tau[5*(closestlayer + 1) + 4] + i) += transmission5;
            *(tauMoon[5*(closestlayer + 1) + 0] + i) += transmissionMoon;
            *(tauMoon[5*(closestlayer + 1) + 1] + i) += transmissionMoon2;
            *(tauMoon[5*(closestlayer + 1) + 2] + i) += transmissionMoon3;
            *(tauMoon[5*(closestlayer + 1) + 3] + i) += transmissionMoon4;
            *(tauMoon[5*(closestlayer + 1) + 4] + i) += transmissionMoon5;
            prheight = cheight;

        }
    }

}
