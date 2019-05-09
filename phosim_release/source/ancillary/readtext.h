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

#ifndef READTEXT_H
#define READTEXT_H

#include <cstdlib>
#include <fstream>
#include <istream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

class readText {

 public:
    readText(std::istream & inStream);
    readText(const std::string & txtFile);
    void readStream(std::istream & inStream);
    const std::string& operator[] (size_t index) const;
    void add(const std::string & line) {
        m_data.push_back(line);
    }
    size_t getSize() {
        return m_data.size();
    }
    bool static goodLine(std::string & line);
    static std::istream& getLine(std::istream & inStream, std::string & line);
    static void readCol(const std::string & txtFile, std::vector<double> & col1, std::vector<double> & col2);
    static void readColArr(const std::string & txtFile, double* & col1, double* & col2);
    static void readMultiCol(const std::string & txtFile, std::vector<double> & col1, std::vector<std::vector<double> > & col2);
    static void readMultiColArr(const std::string & txtFile, double* & col1, double** & col2);

    template <class T> static bool getkey(const std::string & line, const T & key, std::string & value) {
        std::istringstream iss(line);
        T keyName;
        iss >> keyName;
        if (keyName == key) {
            std::getline(iss, value);
            return true;
        }
        return false;
    }
    template <class T> static bool getKey(const std::string & line, const std::string & key, T & value) {
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;
        if (keyName == key) {
            iss >> value;
            return true;
        }
        return false;
    }
    template <class T> static bool getKeyVector(const std::string & line, const std::string & key, std::vector<T> & value) {
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;
        if (keyName == key) {
            size_t i;
            iss >> i;
            if (i >= value.size()) value.resize(i + 1);
            iss >> value[i];
            return true;
        }
        return false;
    }
    template <class T> static bool getKeyMatrix(const std::string & line, const std::string & key, std::vector<std::vector<T> > & value) {
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;
        if (keyName == key) {
            int i, j;
            iss >> i >> j;
            iss >> value[i][j];
            return true;
        }
        return false;
    }
    template <class T> static void get(const std::string & line, const std::string & key, T & value) {
        (void) getKey<T>(line, key, value);
    };
    template <class T> static void get(const std::string & line, const std::string & key, std::vector<T> & value) {
        (void) getKeyVector<T>(line, key, value);
    };
    template <class T> static void get(const std::string & line, const std::string & key, std::vector<std::vector<T> > & value) {
        (void) getKeyMatrix<T>(line, key, value);
    };

    template <class T> static std::string get(const std::string & txtFile, const T & key) {
        std::ifstream inStream(txtFile.c_str());
        if (inStream.is_open()){
            std::string line, value("");
            while (getLine(inStream, line)) {
                if (getkey(line, key, value)) break;
            }
            return value;
        } else {
            std::cout << "Error reading " + txtFile << std::endl;
            exit(1);
        }
    }

    static void readSegmentation(const std::string & txtFile, const std::string & key, std::vector<std::string> & value); 


    void setLine(int index, std::string& line) {
        m_data.at(index) = line;
        return;
    }

    void eraseLine(int index) {
        m_data.erase(m_data.begin() + index);
        return;
    }




 private:
    std::vector<std::string> m_data;
};


#endif
