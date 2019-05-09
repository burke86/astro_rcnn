///
/// @package phosim
/// @file readtext.h
/// @brief class to extract parameters from text files
///
/// @brief Created by
/// @author En-Hsing Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "ancillary/readtext.h"


readText::readText(std::istream & inStream) {
    readStream(inStream);
}

readText::readText(const std::string & txtFile) {
    std::ifstream inStream(txtFile.c_str());
    if (!inStream) {
        std::cout << "Error reading " << txtFile << std::endl;
        exit(1);
    }
    readStream(inStream);
}

void readText::readStream(std::istream & inStream) {
    std::string line;
    while (getLine(inStream, line)) {
        m_data.push_back(line);
    }
}

const std::string& readText::operator[] (const size_t index) const {
    return m_data[index];
}

bool readText::goodLine(std::string & line) {
    //ignore empty lines, lines starting with #, and all chars following #
    std::string::size_type idx = line.find("#");
    if (idx == 0) return false;
    if (idx != std::string::npos) line = line.substr(0, idx - 1);
    if (line.empty()) return false;
    return true;
}

std::istream& readText::getLine(std::istream & inStream, std::string & line) {
    if (std::getline(inStream, line)) {
        while (!goodLine(line)) {
            if (!std::getline(inStream, line)) break;
        }
    }
    return inStream;
}

void readText::readCol(const std::string & txtFile, std::vector<double> & col1, std::vector<double> & col2) {
    std::ifstream inStream(txtFile.c_str());
    if (inStream.is_open()){
        std::string line;
        while (getLine(inStream, line)) {
            std::istringstream iss(line);
            double d1, d2;
            iss>>d1>>d2;
            col1.push_back(d1);
            col2.push_back(d2);
        }
    } else {
        std::cout << "Error reading " + txtFile << std::endl;
        exit(1);
    }
}

void readText::readColArr(const std::string & txtFile, double* & col1, double* & col2) {
    std::ifstream inStream(txtFile.c_str());
    if (inStream.is_open()){
        std::string line;
        size_t m(0);
        while (getLine(inStream, line)) m++;
        col1 = new double[m];
        col2 = new double[m];
        inStream.clear();
        inStream.seekg(0);
        size_t i(0);
        while (getLine(inStream, line)) {
            std::istringstream iss(line);
            iss>>col1[i]>>col2[i];
            i++;
        }
    } else {
        std::cout << "Error reading " + txtFile << std::endl;
        exit(1);
    }
}


void readText::readMultiCol(const std::string & txtFile, std::vector<double> & col1, std::vector<std::vector<double> > & col2) {
    std::ifstream inStream(txtFile.c_str());
    if (inStream.is_open()){
        std::string line;
        while (getLine(inStream, line)) {
            std::istringstream iss(line);
            double d1, v;
            std::vector<double> d2;
            iss>>d1;
            col1.push_back(d1);
            while (iss >> v) {
                d2.push_back(v);
            }
            col2.push_back(d2);
        }
    } else {
        std::cout << "Error reading " + txtFile << std::endl;
        exit(1);
    }
}

void readText::readMultiColArr(const std::string & txtFile, double* & col1, double** & col2) {

    std::ifstream inStream(txtFile.c_str());
    if (inStream.is_open()){
        std::string line;
        size_t m(0);
        size_t n(0);
        while (getLine(inStream, line)) {
            if (m == 0) {
                double v;
                std::istringstream iss(line);
                while (iss >> v) {
                    n++;
                }
            }
            m++;
        }
        col1 = new double[m];
        col2 = new double*[m];
        for (size_t i(0); i < m; i++) {
            col2[i] = new double[n];
        }
        inStream.clear();
        inStream.seekg(0);
        size_t i(0);
        while (getLine(inStream, line)) {
            std::istringstream iss(line);
            iss >> col1[i];
            for (size_t j(0); j < n; j++) {
                iss >> col2[i][j];
            }
            i++;
        }
    } else {
        std::cout << "Error reading " + txtFile << std::endl;
        exit(1);
    }
}

void readText::readSegmentation(const std::string & txtFile, const std::string & key, std::vector<std::string> & value) {
    std::ifstream inStream(txtFile.c_str());
    if (inStream.is_open()){
        std::string line;
        while (getLine(inStream, line)) {
            std::istringstream iss(line);
            std::string keyName;
            iss >> keyName;
            if (keyName == key) {
                int i(0);
                int numAmplifiers;
                iss >> numAmplifiers;
                value.resize(numAmplifiers);
                while (getLine(inStream, value[i])) {
                    i++;
                    if (i == numAmplifiers) break;
                }
                break;
            }
        }
    } else {
        std::cout << "Error reading " + txtFile << std::endl;
        exit(1);
    }

}
