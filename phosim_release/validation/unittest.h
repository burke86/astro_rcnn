#include <stdio.h>

inline void unitTestOutput(double result, double expected, const char *functionName, const char *testName, const char *unit) {

    if (expected != 0.0) {
        if (fabs(result-expected)/expected < 1e-6) {
            printf("%-32s|%-32s|%12lf|%12lf|%-16s|%s\n",
                   functionName,testName,result,expected,unit,"Pass");
        } else {
            printf("%-32s|%-32s|%12lf|%12lf|%-16s|%s\n",
                   functionName,testName,result,expected,unit,"Fail");
        }
    } else {
        if (fabs(result-expected)/1.0 < 1e-6) {
            printf("%-32s|%-32s|%12lf|%12lf|%-16s|%s\n",
                   functionName,testName,result,expected,unit,"Pass");
        } else {
            printf("%-32s|%-32s|%12lf|%12lf|%-16s|%s\n",
                   functionName,testName,result,expected,unit,"Fail");
        }
    }

}

inline void unitTestOutput(double result, double expected, double tolerance, const char *functionName, const char *testName, const char *unit) {

    if (expected != 0.0) {
        if (fabs(result-expected)/expected < tolerance) {
            printf("%-32s|%-32s|%12lf|%12lf|%-16s|%s\n",
                   functionName,testName,result,expected,unit,"Pass");
        } else {
            printf("%-32s|%-32s|%12lf|%12lf|%-16s|%s\n",
                   functionName,testName,result,expected,unit,"Fail");
        }
    } else {
        if (fabs(result-expected)/1.0 < tolerance) {
            printf("%-32s|%-32s|%12lf|%12lf|%-16s|%s\n",
                   functionName,testName,result,expected,unit,"Pass");
        } else {
            printf("%-32s|%-32s|%12lf|%12lf|%-16s|%s\n",
                   functionName,testName,result,expected,unit,"Fail");
        }
    }

}
