///
/// @package phosim
/// @file instrumentfiles.cpp
/// @brief instrument  applications.
///
/// @brief Created by
/// @author Glenn Sembroski (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "instrumentfiles.h"

InstrumentFiles::InstrumentFiles() {
    //Nothing to do
}

InstrumentFiles::~InstrumentFiles() {
    //nothing to do
}

void InstrumentFiles::makeSurfaceMap(std::string opticsFile) {

    // Go through optics_0.txt file thats in the specified instrument directory
    // and create a map of surface names to surface index. Ignore "none" in
    // both key and surface type
    // Also get last device and last surface numbers

    readText opticsPars(opticsFile);

    std::string surfaceName;
    std::string surfaceType;
    _lastSurface = -1;
    int dev = 0;
    for (size_t t(0); t < opticsPars.getSize(); t++) {
        std::istringstream iss(opticsPars[t]);
        iss >> surfaceName;
        iss >> surfaceType;

        if (surfaceName != "none" && surfaceType != "none"){
            fSurfaceMapPos = fSurfaceMap.find(surfaceName);
            if (fSurfaceMapPos != fSurfaceMap.end()){
                fSurfaceMap[surfaceName] = dev;
                dev++;
            }
        }

        if (surfaceType != "none"){
            _lastSurface++;
        }
    }
    return;

}

void InstrumentFiles::makeTrackingFile(std::string trackingFileName, double vistime, double jittertime) {

    std::ofstream outputDataFile(trackingFileName.c_str());
    double starttime = -0.5*vistime;
    double tempf1 = 0;
    double tempf2 = 0;
    double tempf3 = 0;
    double tempf4 = 0;
    double tempf5 = 0;
    double tempf6 = 0;
    double tempf7 = 0;
    double tempf8 = 0;
    double tempf9 = 0;
    long prevj = 0;
    long currj = 0;
    double jitterele[10000];
    double jitterazi[10000];
    double jitterrot[10000];
    double jitterwind[10000];
    double windshake[10000];

    for (long j = 0; j < 10000; j++) {
        if (j == currj) {
            tempf1 = random.uniform()*jittertime;
            if (tempf1 < vistime/10000.0){
                tempf1 = vistime/10000.0;
            }
            prevj = currj;
            currj = static_cast<long>(currj + (tempf1/vistime)*10000.0);
            if (currj == prevj) {
                currj += 1;
            }
            tempf5 = tempf2;
            tempf6 = tempf3;
            tempf7 = tempf4;
            tempf8 = tempf9;
            tempf2 = random.normal()*sqrt((jittertime/2.0)/vistime) + tempf5;
            tempf3 = random.normal()*sqrt((jittertime/2.0)/vistime) + tempf6;
            tempf4 = random.normal()*sqrt((jittertime/2.0)/vistime) + tempf7;
            tempf9 = random.normal()*sqrt((jittertime/2.0)/vistime) + tempf8;
        }
        jitterele[j] = tempf2*(j - prevj)/(currj - prevj) +
            tempf5*(currj - j)/(currj - prevj);
        jitterazi[j] = tempf3*(j - prevj)/(currj - prevj) +
            tempf6*(currj - j)/(currj - prevj);
        jitterrot[j] = tempf4*(j - prevj)/(currj - prevj) +
            tempf7*(currj - j)/(currj - prevj);
        jitterwind[j] = tempf9*(j - prevj)/(currj - prevj) +
            tempf8*(currj - j)/(currj - prevj);
    }

    currj=0;
    tempf2=0.0;
    tempf5=0.0;
        for (long j = 0; j < 10000; j++) {
        if (j == currj) {
            tempf1 = random.normal();
            tempf1=exp(log(90.0)+tempf1*1.4+tempf1*tempf1*0.5+tempf1*tempf1*tempf1*0.8);
            if (tempf1 < vistime/10000.0){
                tempf1 = vistime/10000.0;
            }
            prevj = currj;
            currj = static_cast<long>(currj + (tempf1/vistime)*10000.0);
            if (currj == prevj) {
                currj += 1;
            }
            tempf5 = tempf2;
            tempf2 = random.normal()*sqrt(90.0/vistime) + tempf5;
        }
        windshake[j] = tempf2*(j - prevj)/(currj - prevj) +
            tempf5*(currj - j)/(currj - prevj);
     }

    double meanjitter = jitterwind[5000];
    for (long j = 0; j < 10000; j++) {
        jitterwind[j] = jitterwind[j] - meanjitter;
    }
    for (long j = 4999; j >= 0; j--) {
        jitterwind[j] = jitterwind[j] + jitterwind[j + 1];
    }
    for (long j = 5001; j < 10000; j++) {
        jitterwind[j] = jitterwind[j] + jitterwind[j - 1];
    }
    for (long j = 4999; j >= 0; j--) {
        jitterwind[j] = jitterwind[j]/(5000 - j + 1);
    }
    for (long j = 5001; j < 10000; j++) {
        jitterwind[j] = jitterwind[j]/(j - 5000 + 1);
    }

    for (long j = 0; j < 10000; j++) {
        double time = starttime + vistime/10000.0*j;
        outputDataFile << std::fixed << std::setprecision(6) << time << " " << jitterele[j]
                       << " " << jitterazi[j] << " " << jitterrot[j] << " " << jitterwind[j]
                       << " " << windshake[j] << std::endl;
    }
    outputDataFile.close();
    return;
}


void InstrumentFiles::readoutPars(readText& focalPlaneLayoutPars,
                                  readText& segmentationPars, 
                                  std::string readoutString, int camConfig) 
// From the focalplanelayout.txt file and the segmentation.txt files
// make up and write out the readout_*.pars
{
    //Another random init.
    random.setSeed32Fixed(1000);



    bool useGroup0 = false;
    bool useGroup1 = false;
    bool useGroup2 = false;

    if (camConfig & 1){
        useGroup0 = true;
    }
    if (camConfig & 2){
        useGroup1 = true;
    }
    if (camConfig & 4){
        useGroup2 = true;
    }

    // *****************************************************************
    // iterate through sets of chips (selected by the useGroup variables above)
    // Go through the focalPlaneLayout
    // *****************************************************************
    int numfpLKeys = focalPlaneLayoutPars.getSize();
    int numSegKeys = segmentationPars.getSize();
    for (int keyNum = 0; keyNum < numfpLKeys; keyNum++) {
        std::string line = focalPlaneLayoutPars[keyNum];
        std::istringstream issfpL(line); //Pick up the line
        std::string chipID;
        issfpL>>chipID;                  //Pickup the chip ID (first argument)
        // **********************************************
        // Look for the Group designation
        // **********************************************

        std::string::size_type idx = line.find("Group");
        bool hasWriteOut = false;
        if (idx != std::string::npos){
            std::string groupID = line.substr(idx + 5, 1);// gets character at end
                                                       // of Group
            if (( groupID == "0" && useGroup0 ) ||
                ( groupID == "1" && useGroup1 ) ||
                ( groupID == "2" && useGroup2 ) ) {
                hasWriteOut = true;
            }
        }

        // ********************************************************************
        // If we are writing out this chip (ie it has a known group number) find 
        // the chipID in the segmentation file (its so easy now!)
        // That line will have the number of amplifies for this chip 
        // that follow in the segmentation.txt file. They have keys 
        // like R00_S21_C06
        if (hasWriteOut){
            bool hasKey = false;
            for (int segNum = 0; segNum < numSegKeys; segNum++ ) {
                std::string lineSeg = segmentationPars[segNum];
                std::istringstream issSeg(lineSeg); //Pick up the line
                std::string segChipID;
                issSeg >> segChipID;           //Pickup the chip ID (first argument)

                if (segChipID == chipID) {
                    // We have the chip. Note that this is a chipID of the form 
                    //  R00_S21. There is no amplifier designation.
                    // Get the number of amplifiers for this chip. We make a seperate
                    // output file for each chip
                    hasKey = true;
                    int numAmplifiers = 0;
                    issSeg >> numAmplifiers;   //number ampifier lines (second argument)

                    //Setup Output chip file
                    std::string outputChipFileName = readoutString + chipID + ".pars";
                    std::ofstream outChipStream(outputChipFileName.c_str());

                    std::vector<std::string> amplifiers;
                    for (int j = 0; j <  numAmplifiers; j++ ) {
                        int ampIndx = segNum + j + 1;
                        //std::cout<<"ampIndx: "<<ampIndx<<std::endl;
                        std::string lineAmp = segmentationPars[ampIndx];
                        std::istringstream issAmp(lineAmp); //Pick up the line
                        std::string ampName;
                        issAmp>>ampName;
                        std::vector< double > tokens;
                        double value;
                        while(issAmp>>value){
                            tokens.push_back(value);
                        }
                        //int numTokens=tokens.size();
                        //std::cout<<"ampName,numTokens: "<<ampName<<" "<<numTokens
                        //         <<std::endl;

                        int serialread   = tokens.at(4);
                        int parallelread = tokens.at(5);

                        double mean1 = tokens.at(6);
                        double var1  = tokens.at(7);
                        double gain = mean1*(1 + var1*random.normalFixed()/100);

                        double mean2 = tokens.at(8);
                        double var2  = tokens.at(9);
                        double bias = mean2*(1 + var2*random.normalFixed()/100);

                        double mean3 = tokens.at(10);
                        double var3  = tokens.at(11);
                        double readnoise = mean3*(1 + var3*random.normalFixed()/100);

                        double mean4 = tokens.at(12);
                        double var4  = tokens.at(13);
                        double darkcurrent = mean4*(1 + var4*random.normalFixed()/100);

                        int parallelPrescan = tokens.at(14);
                        int serialOverscan = tokens.at(15);
                        int serialPrescan = tokens.at(16);
                        int parallelOverscan = tokens.at(17);
                        double hotpixel = tokens.at(18);
                        double hotcolumn = tokens.at(19);

                        // Write it all out
                        outChipStream << "serialread    " << j << " " << serialread << std::endl;
                        outChipStream << "parallelread  " << j << " " << parallelread << std::endl;
                        outChipStream << "gain          " << j << " " << gain << std::endl;
                        outChipStream << "bias          " << j << " " << bias << std::endl;
                        outChipStream << "readnoise     " << j << " " << readnoise << std::endl;
                        outChipStream << "darkcurrent   " << j << " " << darkcurrent << std::endl;
                        outChipStream << "parallelprescan   " << j << " " << parallelPrescan << std::endl;
                        outChipStream << "serialoverscan    " << j << " " << serialOverscan << std::endl;
                        outChipStream << "serialprescan     " << j << " " << serialPrescan << std::endl;
                        outChipStream << "paralleloverscan  " << j << " " << parallelOverscan << std::endl;
                        outChipStream << "hotpixelrate  " << j << " " << hotpixel << std::endl;
                        outChipStream << "hotcolumnrate " << j << " " << hotcolumn << std::endl;
                    }//End amlplifier loop
                    outChipStream.close();
                    break;
                } //end of found Chip if
            }   //end of chip loop

            if (!hasKey) {                            //segmentation.txt file!
                std::cout << "Error#2 in readoutPars. Cannot find key:" << chipID
                          << " in segmentation.txt  file" << std::endl;
                return;
            }
        }//end of chip write out if
    }//end segment file chip loop
    return;
}


void InstrumentFiles::focalPlanePars(readText& fpLPars, std::string outChipString, int camConfig, int perturbationMode, std::vector<int> pertSurf, std::vector<int> pertSecondSurf, int pertFlag)

// This function makes the chip_999999_R00_S22_C1.pars type files that have
// all the body and zenike values for the chip. That info comes from the
// focalplanelayout file.  We do watch the camConfig value to decide which
// chip files to makw (cut on Group setting).
{
    // ************************************************************************
    // Now go through the focalplanelayout.txt file making a chip file for
    // any chip in an acceptable Group
    // ***********************************************************************
    // We need to iterate thorough the lines (a line per chip)  from the
    // focalplanelayout.txt we just parsed and find those chips that match our
    // camConfig specified Group.
    // **********************************************************************
    // We use the bit settings in camConfig to specify the allowed groups
    // Ie. camConfig= 7 meas groups 0 and 1 and 2.
    // decode camcomfig
    // **********************************************************************
    bool hasGroup0 = false;
    bool hasGroup1 = false;
    bool hasGroup2 = false;

    if (camConfig & 1){
        hasGroup0 = true;
    }
    if (camConfig & 2){
        hasGroup1 = true;
    }
    if (camConfig & 4){
        hasGroup2 = true;
    }

    // Go through the fpLayoutPars
    int numKeys = fpLPars.getSize();
    for (int keyNum = 0; keyNum < numKeys; keyNum++){
        std::string linefpL = fpLPars[keyNum];
        std::istringstream issfpL(linefpL); //Pick up the line
        std::string chipID;
        issfpL >> chipID;           //Pickup the chip ID (first argument)

        // Look for the Group designation
        std::string::size_type idx = linefpL.find("Group");
        bool hasWriteOut = false;
        if (idx != std::string::npos) {
            std::string groupID = linefpL.substr(idx + 5, 1);
            //gets character at end of Group
            if ((groupID == "0" && hasGroup0) || (groupID == "1" && hasGroup1) ||
                (groupID == "2" && hasGroup2)) {
                hasWriteOut = true;
            }
        }

        if (hasWriteOut) {
            std::string outChipFileName = outChipString + chipID + ".pars";
            std::ofstream outChipFile(outChipFileName.c_str());

            // ********************************************************************
            // Write out the body commands for this chip
            // ********************************************************************
            std::string line = linefpL.substr(idx + 6);
            //focalplanlayout string after "Group"*
            std::istringstream iss(line);

            //Body values first
            for (int i = 0; i < 3; i++) {
                double bodyValue;
                iss >> bodyValue;
                outChipFile << "body " << (_lastSurface + 1)  << " " << i << " "
                            << std::fixed << std::setprecision(7)
                            << bodyValue*M_PI/180.0 << std::endl;
            }
            for (int i = 0; i < 2; i++) {
                double bodyValue;
                iss >> bodyValue;
                outChipFile << "body " << (_lastSurface + 1)  << " " <<  i + 3  << " "
                            << std::fixed << std::setprecision(7) << bodyValue
                            << std::endl;
            }

            double defocus;
            iss >> defocus;
            std::string pertType;
            iss >> pertType;
            double defocus2;
            iss >> defocus2;
            defocus2 /= 1000.0;
            if (perturbationMode == 2 || perturbationMode == 3) defocus2 += defocus;
            outChipFile << "body " << (_lastSurface + 1)  << " " <<  5  << " "
                        << std::fixed << std::setprecision(7) << defocus2
                        << std::endl;


            //zernike values next
            if (pertType == "zern") {
                for (int i = 0; i < 21; i++) {
                    double zernikeValue;
                    iss >> zernikeValue;
                    if (perturbationMode == 2 || perturbationMode == 3)
                    outChipFile << "izernike " << _lastSurface << " " << i << " "
                                << std::scientific << std::setprecision(6)
                                << zernikeValue/1000.0 << std::endl;
                }
            } else if (pertType == "chebyshev") {
                for (int i = 0; i < NCHEB; i++) {
                    double chebyshevValue;
                    iss >> chebyshevValue;
                    if (perturbationMode == 2 || perturbationMode == 3)
                    outChipFile << "ichebyshev " << _lastSurface << " " << i << " "
                                << std::scientific << std::setprecision(6)
                                << chebyshevValue/1000.0 << std::endl;
                }
            } else {
                std::cout << "Error: Unknown perturbation type " << pertType << std::endl;
                return;
            }

            // QE variation last
            double qeVar;
            iss >> qeVar;
            outChipFile << "qevariation " << std::fixed << std::setprecision(6)
                        << qeVar << std::endl;


            //fea files
            if (pertFlag==1) {
                for (size_t t(0); t < pertSurf.size(); t++){
                    outChipFile << "distortion " << pertSurf[t] << " " << pertSecondSurf[t] << std::endl;
                }
            
            }
            
        } //Group test
    } //chip loop




    
    return;
}
// **************************************************************************
