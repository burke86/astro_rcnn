///
/// @package phosim
/// @file main.cpp
/// @brief main for atmosphere program
///
/// @brief Created by
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <time.h>
#include "atmosphere/atmosphere.h"
#include "atmosphere/operator.cpp"

int main(void) {
    // Default values.
    std::string obshistid = "9999";
    double outerx = 500000.0;
    double pix = 256.0; // pixel size in cm
    int telconfig = 0;
    float constrainseeing = -1;

    double seeing;
    double constrainairglow = -1.0;
    char tempstring[4096];
    FILE *fptr, *faptr;
    char outputfilename[4096];
    double zenith, altitude;
    long monthnum;
    long seed = 0;
    double tai = 0.0;

    int atmospheremode = 1;
    int opticsonlymode = 0;
    int turbulencemode = 1;
    int opacitymode = 1;
    long atmosphericDispersion = 1;
    int skipAtmosphere = 0;
    int spaceMode = -1;

    int haveMonth = 0;
    int haveTai = 0;
    int haveSunAlt = 0;
    double sunalt = -90.0;

    int haveMoonRa = 0;
    double moonra = 0.0;
    int haveMoonDec = 0;
    double moondec = 0.0;
    int haveMoonAlt = 0;
    double moonalt = 0.0;
    int haveMoonPhase = 0;
    double moonphase = 0.0;

    int havePointingRA = 0;
    double pointingra = 0.0;
    int havePointingDec = 0;
    double pointingdec = 0.0;
    int haveAltitude = 0;
    int haveAzimuth = 0;
    double azimuth = 0.0;
    int haveDistMoon = 0;
    double dist2moon = 0.0;

    int haveRotTelPos = 0;
    int haveRotSkyPos = 0;
    double rotTelPos = 0.0;
    double rotSkyPos = 0.0;

    int haveFilter = 0;
    int filter = 0;
    int control = 0;

    int haveTemperature = 0;
    double temperature = 20.0;
    int haveTemperatureVariation = 0;
    double temperatureVariation = 0.0;
    int havePressure = 0;
    double pressure = 520.0;

    Atmosphere atmosphere;
    std::string mjdstring;
    int haveMjdString = 0;

    atmosphere.numlevel = 8;
    atmosphere.constrainclouds = -1.0;
    static int maxlevel = 100;
    std::vector<int> cloudscreen(maxlevel, 0);
    cloudscreen[1] = 1;
    cloudscreen[2] = 1;
    std::vector<float> outerscale(maxlevel, -1);

    // Set some default values.
    atmosphere.datadir = "../data";
    atmosphere.instrdir = "../data/lsst";
    std::string workdir(".");

    // Read parameters from stdin.
    readText pars(std::cin);
    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;
        readText::get(line, "workdir", workdir);
        readText::get(line, "instrdir", atmosphere.instrdir);
        readText::get(line, "datadir", atmosphere.datadir);
        readText::get(line, "numlevel", atmosphere.numlevel);
        readText::get(line, "telconfig", telconfig);
        readText::get(line, "obshistid", obshistid);
        readText::get(line, "outerx", outerx);
        readText::get(line, "cloudscreen", cloudscreen);
        readText::get(line, "outerscale", outerscale);
        readText::get(line, "obsseed", seed);
        if (readText::getKey(line, "monthnum", monthnum)) haveMonth=1;
        if (readText::getKey(line, "solaralt", sunalt)) haveSunAlt=1;
        if (readText::getKey(line, "moonra", moonra)) haveMoonRa=1;
        if (readText::getKey(line, "moondec", moondec)) haveMoonDec=1;
        if (readText::getKey(line, "moonalt", moonalt)) haveMoonAlt=1;
        if (readText::getKey(line, "phaseang", moonphase)) haveMoonPhase=1;
        if (readText::getKey(line, "tai", tai)) haveTai=1;
        if (readText::getKey(line, "pointingra", pointingra)) havePointingRA=1;
        if (readText::getKey(line, "pointingdec", pointingdec)) havePointingDec=1;
        if (readText::getKey(line, "azimuth", azimuth)) haveAzimuth=1;
        if (readText::getKey(line, "moondist", dist2moon)) haveDistMoon=1;
        if (readText::getKey(line, "tai", mjdstring)) haveMjdString=1;
        if (readText::getKey(line, "rotationangle", rotSkyPos)) haveRotSkyPos=1;
        if (readText::getKey(line, "spiderangle", rotTelPos)) haveRotTelPos=1;
        if (readText::getKey(line, "filter", filter)) haveFilter=1;
        if (readText::getKey(line, "temperature", temperature)) haveTemperature=1;
        if (readText::getKey(line, "tempvar", temperatureVariation)) haveTemperatureVariation=1;
        if (readText::getKey(line, "pressure", pressure)) havePressure=1;
        readText::get(line, "atmospheremode", atmospheremode);
        readText::get(line, "opticsonlymode", opticsonlymode);
        readText::get(line, "atmosphericdispersion", atmosphericDispersion);
        readText::get(line, "spacemode", spaceMode);
        readText::get(line, "control", control);

        // check to see if atmosphere is off
        if (keyName == "cleareverything" || opticsonlymode == 1 || atmospheremode == 0) {
             skipAtmosphere = 1;
        }
        if (keyName == "clearturbulence") turbulencemode = 0;
        if (keyName == "clearopacity") opacitymode = 0;
        if (turbulencemode == 0 && atmosphericDispersion == 0.0 && opacitymode == 0) {
            skipAtmosphere = 1;
        }
        if (control == 1) skipAtmosphere=0;
 
        if (readText::getKey(line, "constrainseeing", seeing)) constrainseeing = seeing/2.35482;
        if (readText::getKey(line, "altitude", altitude)) {
            haveAltitude = 1;
            zenith = 90 - altitude;
        }
        readText::get(line, "zenith", zenith);
        readText::getKey(line, "constrainclouds", atmosphere.constrainclouds);
        readText::getKey(line, "constrainairglow", constrainairglow);

}

    atmosphere.random.setSeed32(seed);
    atmosphere.random.unwind(10000);

    std::cout << "------------------------------------------------------------------------------------------" << std::endl
    << "Operator" << std::endl
    << "------------------------------------------------------------------------------------------" << std::endl;

     // Get site location
    double latitude = 0.0, longitude = 0.0;
    std:: string ss;
    ss = atmosphere.instrdir + "/location.txt";
    std::ifstream newStream(ss.c_str());
    if (newStream) {
        readText locationPars(atmosphere.instrdir + "/location.txt");
        for (size_t t(0); t < locationPars.getSize(); t++) {
            std::string line(locationPars[t]);
            readText::getKey(line, "latitude", latitude);
            readText::getKey(line, "longitude", longitude);
            readText::get(line, "groundlevel", atmosphere.groundlevel);
            if (spaceMode == -1) readText::getKey(line, "spacemode", spaceMode);
        }
    }
    if (spaceMode == -1) spaceMode = 0;

    if (spaceMode > 0) {
        latitude = 0;
        longitude = 0;
        skipAtmosphere = 1;
        // say spaceMode 1 = low eath orbit
        // 2 = LagrangePt. L1
        // 3=L2, (JWST) ...etc
    }

    // Get time
    if (haveTai == 0) {
    redo2:;
        tai = 51544.5 + atmosphere.random.uniform()*(36525);
        double lstt = localSiderealTime(tai, longitude);
        double tsunAlt, tsunAz, tsunRa, tsunDec;
        sunPos(tai,&tsunRa,&tsunDec);
        raDectoAltAz(lstt, latitude, tsunRa, tsunDec, &tsunAlt, &tsunAz);
        if ((tsunAlt > -18) && (spaceMode == 0)) goto redo2;
    }
    if (mjdstring == "now" || mjdstring == "tonight") {
        time_t now = time(0);
        struct tm tstruct;
        tstruct = *gmtime(&now);
        long curryr, currmo, currdy, currhr, currmn, currsc;
        curryr = tstruct.tm_year+1900;
        currmo = tstruct.tm_mon+1;
        currdy = tstruct.tm_mday;
        currhr = tstruct.tm_hour;
        currmn = tstruct.tm_min;
        currsc = tstruct.tm_sec;
        if ((currmo == 1) || (currmo ==2)) {
            curryr = curryr -1;
            currmo = currmo + 12;
        }
        long curra, currb;
        curra = floor ( curryr /100);
        currb= 2 - curra + floor(curra/4);
        double currjd = floor(365.25*(curryr + 4716)) + floor(30.6001*(currmo + 1)) + currdy + currb - 1524.5;
        double currmjd = currjd - 2400000.5;
        currmjd += currhr/24.0 + currmn/24.0/60.0 + currsc/24.0/60.0/60.0;
        tai = currmjd;
    redo:;
        if (mjdstring == "tonight") tai = currmjd + atmosphere.random.uniform();
        double lstt = localSiderealTime(tai, longitude);
        double tsunAlt, tsunAz, tsunRa, tsunDec;
        sunPos(tai,&tsunRa,&tsunDec);
        raDectoAltAz(lstt, latitude, tsunRa, tsunDec, &tsunAlt, &tsunAz);
        if ((mjdstring == "tonight") && (tsunAlt > -18) && (spaceMode == 0)) goto redo;
        if ((mjdstring == "now") && (tsunAlt > -18) && (spaceMode == 0)) {
            std::cout << "Warning:  This observation is in daytime or twilight." << std::endl;
        }
    }


    //Calendar Date
    long year, month, day;
    double frac;
    jdToCalendar(tai, &month, &day, &year, &frac);
    if (haveMonth == 1) {
        std::cout << "Warning Redundant Setting: Using the Month input as " << month << " instead of " << monthnum << "." << std::endl;
    } else {
        monthnum = month;
    }

    //Day of year
    double dayOfYear = calculateDayOfYear(month, day, year, frac);

    //Local Sidereal Time
    double lst = localSiderealTime(tai, longitude);

    //Sun position
    double sunAlt, sunAz, sunRa, sunDec;
    sunPos(tai,&sunRa,&sunDec);
    raDectoAltAz(lst, latitude, sunRa, sunDec, &sunAlt, &sunAz);
    if (haveSunAlt == 1) {
        std::cout << "Warning Redundant Setting: Using the Sun Altitude input as " << sunalt << " instead of " << sunAlt << "." << std::endl;
    } else {
        sunalt = sunAlt;
    }


    atmosphere.random.setSeed32(seed);
    atmosphere.random.unwind(10000);
    double newra=0.0, newdec=0.0;

    if ((havePointingRA == 0) && (havePointingDec == 0) && (haveAltitude == 0) && (haveAzimuth == 0)) {
        pointingra = atmosphere.random.uniform()*360.0;
        pointingdec = (acos(2.0*atmosphere.random.uniform()-1.0)-PI/2.0)/DEGREE;
        raDectoAltAz(lst, latitude, pointingra, pointingdec, &altitude, &azimuth);
    } else {
        if (haveAltitude == 0) altitude = (acos(0.5*atmosphere.random.uniform()-1.0)-PI/2.0)/DEGREE;
        if (haveAzimuth == 0) azimuth = atmosphere.random.uniform()*360.;
        altAztoRaDec(lst, latitude, altitude, azimuth, &newra, &newdec);
        if (havePointingRA == 0) {
            pointingra=newra;
        } else {
            if ((haveAltitude != 0) && (haveAzimuth !=0)) std::cout << "Warning Redundant Setting: Using the Right Ascension input as " << pointingra << " instead of " << newra << "." << std::endl;
        }
        if (havePointingDec == 0) {
            pointingdec=newdec;
        } else {
            if ((haveAltitude != 0) && (haveAzimuth !=0)) std::cout << "Warning Redundant Setting: Using the Declination input as " << pointingdec << " instead of " << newdec << "." << std::endl;
        }
    }


    //Moon position
    double moonAlt, moonAz, moonRa, moonDec, newMoonPhase;
    moonPos(tai,&moonRa,&moonDec);
    raDectoAltAz(lst, latitude, moonRa, moonDec, &moonAlt, &moonAz);
    newMoonPhase = moonPhase(sunRa,sunDec,moonRa,moonDec);
    double distMoon = distSphere(pointingra, pointingdec, moonRa, moonDec);
    if (haveMoonRa == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon RA input as " << moonra << " instead of " << moonRa << "." << std::endl;
    } else {
        moonra = moonRa;
    }
    if (haveMoonDec == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon Dec input as " << moondec << " instead of " << moonDec << "." << std::endl;
    } else {
        moondec = moonDec;
    }
    if (haveMoonAlt == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon Alt input as " << moonalt << " instead of " << moonAlt << "." << std::endl;
    } else {
        moonalt = moonAlt;
    }
    if (haveMoonPhase == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon Phase input as " << moonphase << " instead of " << newMoonPhase << "." << std::endl;
    } else {
        moonphase = newMoonPhase;
    }
    if (haveDistMoon == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon Distance input as " << dist2moon << " instead of " << distMoon << "." << std::endl;
    } else {
        dist2moon = distMoon;
    }

    double ha = lst*15.0 - pointingra;

    // rottelpos = rotator angle
    // rotskypos + rotator = parallactic angle
    double newRotSkyPos = 0.0;
    if (haveRotTelPos == 1) {
        newRotSkyPos=calculateRotSkyPos(ha, latitude, pointingdec, rotTelPos);
        if (haveRotSkyPos == 1) {
            std::cout << "Warning Redundant Setting: Using the Parallactic - Rotator input as " << rotSkyPos << " instead of " << newRotSkyPos << std::endl;
        } else {
            rotSkyPos = newRotSkyPos;
        }
    } else {
        if (haveRotSkyPos == 1) {
            rotTelPos=calculateRotTelPos(ha, latitude, pointingdec, rotSkyPos);
        } else {
            rotTelPos = atmosphere.random.uniform()*360.0;
            rotSkyPos=calculateRotSkyPos(ha, latitude, pointingdec, rotTelPos);
        }
    }

    // Filter
    int numberFilter=0;
    char opticsFile[4096];
    int goodFilter=1;
    while (goodFilter) {
        sprintf(opticsFile, "%s/optics_%d.txt", atmosphere.instrdir.c_str(), numberFilter);
        std::ifstream f(opticsFile);
        goodFilter = f.good();
        if (goodFilter) numberFilter++;
    }
    if (haveFilter == 1) {
        if ((filter < 0) || (filter > (numberFilter-1))) {
            std::cout << "Filter/Optics Configuration is out of range." << std::endl;
            exit(1);
        }
    } else {
        filter = floor(atmosphere.random.uniform()*numberFilter);
    }

    // Weather & Climate
    if (spaceMode <= 0) {
        double meanHigh;
        if (abs(latitude) < 18) {
            meanHigh = 27.0;
        } else {
            meanHigh = 27.0 - 0.75 * (abs(latitude) - 18.0);
        }
        double oceanDistance = 1000.0;
        char filename[4096];
        sprintf(filename, "%s/atmosphere/earth.txt", atmosphere.datadir.c_str());
        FILE *f = fopen(filename, "r");
        double mindist = 1e30;
        for (long i = 0; i < 23328 ; i++) {
            double x = 0.0, y = 0.0, o = 0.0, z = 0.0;
            fscanf(f,"%lf %lf %lf %lf\n",&x,&y,&o,&z);
            double phi, theta;
            phi=x/432.0*2*PI;
            theta=y/216.0*PI-PI/2.0;
            double distance = acos(sin(latitude*DEGREE)*sin(theta) + cos(latitude*DEGREE)*cos(theta)*cos(longitude*DEGREE-phi));
            if (distance < mindist) {
                mindist=distance;
                oceanDistance = o;
            }
        }
        fclose(f);

        double dayNightDifference = 7.0 + 1.5*pow(oceanDistance, 0.2);
        double janJulyDifference = 1.0 + 0.13*latitude*pow(oceanDistance, 0.2);
        double seasonalLag = 45.0 - 6.0*pow(oceanDistance, 0.2);
        double Tmax = meanHigh + 0.5*janJulyDifference*sin(2*PI*(dayOfYear - 77.0 + seasonalLag)/365.25);
        double Tmin = Tmax - dayNightDifference;
        double delta = 23.45*sin(360.0/365.0*(284.0+dayOfYear));
        double omega = 2.0/15.0*acos(-tan(latitude*DEGREE)*tan(delta*DEGREE))/DEGREE;
        double thermalMax = 12.0 + 2.2;
        double tSunset = 12.0 + omega/2.0;
        double hour = frac*24.0;
        double localTemperature;
        if ((hour >= thermalMax - omega/2.0) && (hour <= tSunset)) {
            localTemperature = Tmin + (Tmax - Tmin)*cos(PI/omega*(hour - thermalMax));
        } else {
            double k = omega/PI*atan(PI/omega*(tSunset-thermalMax));
            if (hour < thermalMax - omega/2.0) hour = hour + 24.0;
            localTemperature = Tmin + (Tmax - Tmin)*cos(PI/omega*(tSunset - thermalMax))*exp(-(hour-tSunset)/k);
        }
        double randomT = 1.0 + pow(oceanDistance, 0.2);
        localTemperature += randomT*atmosphere.random.normalCorrel(tai + 12000, 2.0);
        double lapseRate = 4.65;
        localTemperature -= lapseRate*atmosphere.groundlevel/1000.0;
        if (haveTemperature == 0) temperature = localTemperature;

        // previous hour temperature
        double prTemperature;
        double prHour = frac*24.0 - 1.0;
        if ((prHour >= thermalMax - omega/2.0) && (prHour <= tSunset)) {
            prTemperature = Tmin + (Tmax - Tmin)*cos(PI/omega*(prHour - thermalMax));
        } else {
            double k = omega/PI*atan(PI/omega*(tSunset-thermalMax));
            if (prHour < thermalMax - omega/2.0) prHour = prHour + 24.0;
            prTemperature = Tmin + (Tmax - Tmin)*cos(PI/omega*(tSunset - thermalMax))*exp(-(prHour-tSunset)/k);
        }
        prTemperature += randomT*atmosphere.random.normalCorrel(tai + 12000 - 1.0/24.0, 2.0);
        prTemperature -= lapseRate*atmosphere.groundlevel/1000.0;
        if (haveTemperatureVariation == 0) temperatureVariation = (localTemperature - prTemperature);

    } else {
        if (haveTemperature == 0) temperature = 20.0;
    }
    double waterVaporPressure = 0.0;
    double relativeHumidity = 0.60 + 0.20*atmosphere.random.normalCorrel(tai + 14000, 2.0) - 0.04*atmosphere.groundlevel/1000.0;
    if (relativeHumidity > 1.0) {
        if (temperature > 0) std::cout << "Warning:  It is raining." << std::endl;
        if (temperature <= 0) std::cout << "Warning:  It is snowing." << std::endl;
        relativeHumidity=1.0;
    }
    if (relativeHumidity < 0.0) relativeHumidity=0.0;
    double temperatureKelvin = temperature + 273.15;
    if (temperature > 0) {
        double Tc = 647.096;
        double theta = 1.0 - temperatureKelvin/Tc;
        double Pc=220640.0;
        double c1=-7.85951783;
        double c2=1.84408259;
        double c3=-11.7866497;
        double c4=22.6807411;
        double c5=-15.9618719;
        double c6=1.80122502;
        waterVaporPressure = Pc*100.0/101325.0*760.0*
            exp(Tc/temperatureKelvin*(c1*theta+c2*pow(theta,1.5)+c3*pow(theta,3.0)+
                                      c4*pow(theta,3.5)+c5*pow(theta,4.0)+
                                      c6*pow(theta,7.5)));
    } else {
        double Tn = 273.16;
        double theta = temperatureKelvin/Tn;
        double Pn = 6.11657;
        double a0 = -13.928169;
        double a1 = 34.707823;
        waterVaporPressure = Pn*100.0/101325.0*760.0*
            exp(a0*(1 - pow(theta, -1.5)) + a1*(1 - pow(theta, -1.25)));
    }
    waterVaporPressure=relativeHumidity*waterVaporPressure;

    if (havePressure == 0) {
        if (spaceMode <= 0) {
            double pressureScaleHeight = 7640.0;
            pressure = 760.0 * exp(-atmosphere.groundlevel/pressureScaleHeight);
            pressure = pressure * (1.0 + 0.005*atmosphere.random.normalCorrel(tai + 13000, 2.0));
        } else {
            pressure = 520.0;
        }
    }

    // Calculated inputs
    sprintf(outputfilename, "%s/obs_%s.pars", workdir.c_str(), obshistid.c_str());
    fptr = fopen(outputfilename, "a+");
    fprintf(fptr, "pointingra %lf\n", pointingra);
    fprintf(fptr, "pointingdec %lf\n", pointingdec);
    fprintf(fptr, "altitude %lf\n", altitude);
    fprintf(fptr, "azimuth %lf\n", azimuth);
    fprintf(fptr, "tai %lf\n", tai);
    fprintf(fptr, "moonra %lf\n", moonra);
    fprintf(fptr, "moondec %lf\n", moondec);
    fprintf(fptr, "phaseang %lf\n", moonphase);
    fprintf(fptr, "moondist %lf\n", dist2moon);
    fprintf(fptr, "moonalt %lf\n", moonalt);
    fprintf(fptr, "monthnum %ld\n", monthnum);
    fprintf(fptr, "solaralt %lf\n", sunalt);
    fprintf(fptr, "spiderangle %lf\n", rotTelPos);
    fprintf(fptr, "rotationangle %lf\n", rotSkyPos);
    fprintf(fptr, "filter %d\n", filter);
    fprintf(fptr, "dayofyear %lf\n", dayOfYear);
    fprintf(fptr, "temperature %lf\n", temperature);
    fprintf(fptr, "pressure %lf\n", pressure);
    fprintf(fptr, "tempvar %lf\n", temperatureVariation);
    fprintf(fptr, "waterpressure %lf\n", waterVaporPressure);
    fprintf(fptr, "hourangle %lf\n", ha);
    fprintf(fptr, "lst %lf\n", lst);

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;

    printf("[rightascension] Pointing Right Ascension (degrees)        %lf\n", pointingra);
    printf("[declination] Pointing Declination (degrees)               %lf\n", pointingdec);
    if (spaceMode == 0) {
        printf("[altitude] Pointing Altitude (degrees)                     %lf\n", altitude);
        printf("[azimuth] Pointing Azimuth (degrees)                       %lf\n", azimuth);
    }
    printf("[mjd] Modified Julian Date (days)                          %lf\n", tai);
    printf("[date] Date (month/day/year/day-fraction)                  %02ld/%02ld/%2ld/%lf\n", monthnum,day,year,frac);
    printf("[rottelpos] Rotator Angle (degrees)                        %lf\n", rotTelPos);
    printf("[moonra] Moon Right Ascension (degrees)                    %lf\n", moonra);
    printf("[moondec] Moon Declination (degrees)                       %lf\n", moondec);
    printf("[moonphase] Moon Phase (percent)                           %lf\n", moonphase);
    printf("[moondist] Distance to Moon (degrees)                      %lf\n", dist2moon);
    printf("[moonalt] Altitude of Moon (degrees)                       %lf\n", moonalt);
    printf("[sunalt] Altitude of Sun (degrees)                         %lf\n", sunalt);
    printf("[rotskypos] Parallactic - Rotator Angle (degrees)          %lf\n", rotSkyPos);
    printf("[pressure] Ambient Pressure (mmHg)                         %lf\n", pressure);
    printf("[waterpressure] Water Vapor Pressure (mmHg)                %lf\n", waterVaporPressure);
    printf("[temperature] Ambient Temperature (degrees C)              %lf\n", temperature);
    printf("[tempvar] Temperature Variation (degrees/hr)               %lf\n", temperatureVariation);
    double totalseeing=constrainseeing;
    if (skipAtmosphere == 0) {

        std::cout << "------------------------------------------------------------------------------------------" << std::endl
        << "Atmosphere Creator" << std::endl
        << "------------------------------------------------------------------------------------------" << std::endl;

        atmosphere.osests.reserve(atmosphere.numlevel);
        atmosphere.altitudes.reserve(atmosphere.numlevel);
        atmosphere.jests.reserve(atmosphere.numlevel);

        sprintf(outputfilename, "%s/atmosphere_%s.pars", workdir.c_str(), obshistid.c_str());

        atmosphere.createAtmosphere(monthnum, constrainseeing, outputfilename, cloudscreen, seed, tai, &totalseeing);

        // overwrite outer scales if they have been provided.
        for (int m(0); m < atmosphere.numlevel; m++) {
            if (outerscale[m] >= 0) atmosphere.osests[m] = outerscale[m];
        }

        faptr = fopen(outputfilename, "a+");
        for (int i = 0; i < atmosphere.numlevel; i++) {
            printf("Creating layer %d.\n", i);
            double outers = atmosphere.osests[i]*100.0;
            sprintf(tempstring, "%s/atmospherescreen_%s_%d", workdir.c_str(), obshistid.c_str(), i);
            atmosphere.turb2d(seed*10 + i, seeing, outerx, outers, zenith, 0.5, tempstring);
            fprintf(faptr, "atmospherefile %d atmospherescreen_%s_%d\n", i, obshistid.c_str(), i);
        }
        for (int i = 0; i < atmosphere.numlevel; i++) {
            if (cloudscreen[i]) {
                double height = (atmosphere.altitudes[i] - atmosphere.groundlevel)/1000.0;
                sprintf(tempstring, "%s/cloudscreen_%s_%d", workdir.c_str(), obshistid.c_str(), i);
                atmosphere.cloud(seed*10 + i, height, pix, tempstring);
                fprintf(faptr, "cloudfile %d cloudscreen_%s_%d\n", i, obshistid.c_str(), i);
            }
        }

        sprintf(tempstring, "%s/airglowscreen_%s", workdir.c_str(), obshistid.c_str());
        atmosphere.airglow(seed*10, tempstring);

        atmosphere.random.setSeed32(seed);
        atmosphere.random.unwind(10000);
        if (constrainairglow > 0.0) {
            fprintf(faptr, "airglowpintensity %.3f\n", constrainairglow);
            fprintf(faptr, "airglowcintensity %.3f\n", constrainairglow);
        } else {
            fprintf(faptr, "airglowpintensity %.3f\n",
                    (21.7 + 0.2*atmosphere.random.normalCorrel(tai, 500.0/(24.0*3600.0))
                     + 2.5*log10(cos(zenith*PI/180.0))
                     + (0.50 + 0.50*sin(2*PI*(tai - 54466.0)/(365.24*11.0)
                                        - PI/2.0))));
            fprintf(faptr, "airglowcintensity %.3f\n",
                    (21.7 + 0.3*atmosphere.random.normalCorrel(tai+10000, 500.0/(24.0*3600.0))
                     + 2.5*log10(cos(zenith*PI/180.0))));
        }
        fclose(faptr);

    }

    fprintf(fptr, "constrainseeing %lf\n", totalseeing);
    fclose(fptr);

    return(0);
}
