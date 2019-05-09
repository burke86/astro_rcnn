///
/// @package phosim
/// @file test_readtext.cpp
/// @brief Unit tests for readText class.
///
/// @brief Created by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <math.h>
#include <stdio.h>
#include "ancillary/readtext.h"
#include "validation/unittest.h"

void test_readText() {

    int my_integer;
    long my_long;
    float my_float;
    double my_double;
    std::vector<double> seeing(7);
    std::vector<std::vector<double> > body(13, std::vector<double>(6));

    readText pars(std::cin);

    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        readText::get(line, "seeing", seeing);
        readText::get(line, "body", body);
        readText::get(line, "my_integer", my_integer);
        readText::get(line, "my_long", my_long);
        readText::get(line, "my_float", my_float);
        readText::get(line, "my_double", my_double);
    }

    std::string chipid("R22_S11");
    double centerx, centery;
    std::istringstream iss(readText::get("test.txt", chipid));
    iss >> centerx >> centery;

    unitTestOutput(pars.getSize()/100.0, 0.94, "readText::getSize", "File length", "x 100");
    unitTestOutput(my_double, 3.14, "readText::get", "Get parameter", "none");
    unitTestOutput(seeing[6], 0.230562, "readText::get", "Get vector", "none");
    unitTestOutput(body[8][1], 3.16162571e+00, "readText::get", "Get matrix", "none");
    unitTestOutput(centery, 0.04, "readText::get", "Read focalplanelayout", "none");

    std::vector<double> col1, col2;
    readText::readCol("../../data/atmosphere/nprofile.txt", col1, col2);
    unitTestOutput(col2[9]/1e20, 0.23959746, "readText::readCol", "Read 2-column file","none");

    std::vector<double> col3;
    std::vector<std::vector<double> > col4;
    readText::readMultiCol("../../data/atmosphere/o2cs.txt", col3, col4);
    unitTestOutput(col4[9][1]*1e32, 0.813033, "readText::readMultiCol", "Read multi-column file", "none");
    unitTestOutput(col4.size()/1000000.0, 0.180001, "readText::readMultiCol", "Column size", "x 1000000");

    double *col5, *col6;
    readText::readColArr("../../data/atmosphere/nprofile.txt", col5, col6);
    unitTestOutput(col6[9]/1e20,0.23959746, "readText::readColArr", "Read 2-column file", "none");

    double *col7;
    double **col8;
    readText::readMultiColArr("../../data/atmosphere/o2cs.txt", col7, col8);
    unitTestOutput(col8[9][1]*1e32, 0.813033, "readText::readMultiColArr", "Read multi-column file", "none");



}

int main() {
    test_readText();
    return 0;
}
