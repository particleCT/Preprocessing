#ifndef TVcorrection_h
#define TVcorrection_h
// To read the calibration constants from a file for the T,V energy detector
// corrections
// R.P. Johnson  5/22/2016 to encapculate the TV calibration code of Vladimir
// Bashkirov
// R.P. Johnson  10/2/2016 added a binlinear interpolation of the TV map
// R.P. Johnson  10/21/2016 added optional checks on date and run ranges

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "Util.h"

using namespace std;

class TVcorrection {

#define nStage 5
#define nPix 380

  float tPlaneShifts[4][4];  // planes sensors array for sensors shifts (stored
                             // in TVcorr.txt but no longer used)
  float TVmap[nStage][nPix]; // 5 TV maps with 1cm pixels, 10x38 cm2 area, total
                             // 10x38=380 pixels

public:
  float ped[5];   // Calibrated pedestals (the program will generally use values
                  // generated on the fly, however)
  float Eempt[6]; // ADC pedestals and energy depositions from Empty run

  TVcorrection(const char *TVfile, int year, int month, int day,
               int run) { // pass 0 for all of year, month, day, run to avoid
                          // checks on those values
    cout << "TVcorrection: opening TV correction file " << TVfile << endl;
    fstream TVcalfile(TVfile, ios_base::in); // Open text file with TV-corr etc. data
    if (!TVcalfile.is_open()) {
      perror("Error opening TV correction file");
      return;
    }
    for (int i = 0; i < nStage; ++i) {
      for (int j = 0; j < nPix; ++j) {
        TVcalfile >> TVmap[i][j];
      }
    }
    for (int i = 0; i < nStage; ++i) {
      TVcalfile >> ped[i];
      cout << ped[i] << " ";
    }
    for (int i = 0; i < 6; ++i) {
      TVcalfile >> Eempt[i];
      cout << Eempt[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; ++j)
        TVcalfile >> tPlaneShifts[i][j];
    }
    TVcalfile.close();
    if (ped[0] == 0 || ped[4] == 0 || Eempt[5] < 100.) {
      perror("TVcorrection: Error reading TV correction file.");
      exit(1);
    }
    string minDate = "";
    string maxDate = "";
    int minRun = 0;
    int maxRun = 0;
    Util util;
    string line;
    ifstream infile(TVfile);
    if (infile) {
      cout << endl << "Echoing comment lines from the TV calibration file " << TVfile << ": " << endl;
      while (getline(infile, line)) {
        size_t found = line.find_first_not_of(" ");
        if (line[found] == '#') {
          cout << line << endl;
          string key, value;
          if (util.getKeyValue(line, key, value)) {
            if (key == "minDate")
              minDate = value;
            if (key == "maxDate")
              maxDate = value;
            char *ePtr;
            if (key == "minRun")
              minRun = strtol(value.c_str(), &ePtr, 10);
            if (key == "maxRun")
              maxRun = strtol(value.c_str(), &ePtr, 10);
          }
        }
      }
    }
    infile.close();
    bool checked = false;
    if (minRun < maxRun && !(year == 0 && month == 0 && day == 0 && run == 0)) {
      if (minDate != "" && maxDate != "") {
        int yr1, mn1, dy1;
        if (util.parseDate(minDate, yr1, mn1, dy1)) {
          int yr2, mn2, dy2;
          if (util.parseDate(maxDate, yr2, mn2, dy2)) {
            int minD = (yr1 - 1) * 372 + mn1 * 31 + dy1;
            int maxD = (yr2 - 1) * 372 + mn2 * 31 + dy2;
            int D = (year - 1) * 372 + month * 31 + day;
            if (D <= maxD && D >= minD && run <= maxRun && run >= minRun) {
              cout << "TVcorrection: the run and date of the data match the "
                      "calibration file ranges.  Success!" << endl;
              checked = true;
            } else {
              cout << "TVcorrection: the calibration file date and run are not "
                      "appropriate for these data with year=" << year << " month= " << month << " day= " << day
                   << " run= " << run << endl;
              exit(-1);
            }
          }
        }
      }
    }
    cout << "TVcorrection: completed initialization of the TV corrections.\n";
    if (!checked) {
      cout << "********** WARNING!  No check has been made of the validity of "
              "the date and run ranges of the TVcorrection calibration "
              "constants! **************" << endl;
    }
    cout << "   Pedestals = " << ped[0] << " " << ped[1] << " " << ped[2] << " " << ped[3] << " " << ped[4] << endl;
  }

  float corrFactorInt(float TVmapr[nStage][nPix], int stage, float T, float V,
                      bool &inBounds) { // Can be called with an external TVmap
    inBounds = true;                    // No class members or variables are referenced by this function
    int tPix = floor(0.1 * (T + 190.));
    if (tPix < 0) {
      tPix = 0;
      inBounds = false;
    }
    if (tPix > 37) {
      tPix = 37;
      inBounds = false;
    }
    int vPix = floor(0.1 * (V + 50.));
    if (vPix < 0) {
      vPix = 0;
      inBounds = false;
    }
    if (vPix > 9) {
      vPix = 9;
      inBounds = false;
    }
    if (inBounds) {
      int tLow = floor(0.1 * (T + 185.));
      int tHigh = tLow + 1;
      int vLow = floor(0.1 * (V + 45.));
      int vHigh = vLow + 1;
      if (tLow < 37 && tHigh > 0 && vLow < 9 && vHigh > 0) { // Bilinear interpolation
        float f00 = TVmapr[stage][vLow + 10 * tLow];
        float f10 = TVmapr[stage][vLow + 10 * tHigh];
        float f01 = TVmapr[stage][vHigh + 10 * tLow];
        float f11 = TVmapr[stage][vHigh + 10 * tHigh];
        if (f00 > 0.9 || f10 > 0.9 || f01 > 0.9 || f11 > 0.9) { // Can't include any bad pixels in the interpolation
          return TVmapr[stage][vPix + 10 * tPix];
        } else {
          float T0 = -185.0 + 10.0 * tLow;
          float V0 = -45.0 + 10.0 * vLow;
          float VS = (V - V0) / 10.;
          float TS = (T - T0) / 10.;
          float result = f00 + (f10 - f00) * TS + (f01 - f00) * VS + (f11 + f00 - f10 - f01) * TS * VS;
          return result;
        }
      } else {
        return TVmapr[stage][vPix + 10 * tPix];
      }
    } else {
      return TVmapr[stage][vPix + 10 * tPix];
    }
  }

  float corrFactor(int stage, float T, float V, bool &inBounds) { // Return the calibration factor for
                                                                  // location T,V in the specified stage
    return corrFactorInt(TVmap, stage, T, V, inBounds);
  }

}; // end of TVcorrection class
#endif