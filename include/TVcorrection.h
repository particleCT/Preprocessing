#ifndef TVcorrection_h
#define TVcorrection_h
// To read the calibration constants from a file for the T,V energy detector corrections
// R.P. Johnson  5/22/2016 to encapculate the TV calibration code of Vladimir Bashkirov
// R.P. Johnson  10/2/2016 added a binlinear interpolation of the TV map
// R.P. Johnson  10/21/2016 added optional checks on date and run ranges

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "Util.h"

using namespace std;
#define nStage 5
#define nPix 380

class TVcorrection {
public:
  float tPlaneShifts[4][4]  = {{0}}; // planes sensors array for sensors shifts (stored in TVcorr.txt but no longer used)
  float TVmap[nStage][nPix] = {{0}}; // 5 TV maps with 1cm pixels, 10x38 cm2 area, total 10x38=380 pixels
  float ped[5];   // Calibrated pedestals (the program will generally use values  generated on the fly, however)
  float Eempt[6]; // ADC pedestals and energy depositions from Empty run

  // Constructor -- read from the TV calibration file
  TVcorrection(const char *TVfile, int year, int month, int day, int run) { // pass 0 for all of year, month, day, run to avoid checks on those values
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
    cout << "TVcorrection: completed initialization of the TV corrections.\n";
    cout << "   Pedestals = " << ped[0] << " " << ped[1] << " " << ped[2] << " " << ped[3] << " " << ped[4] << endl;
  }

  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  float corrFactor(int stage, float T, float V,bool &inBounds) {
    // Can be called with an external TVmap
    // No class members or variables are referenced by this function
    inBounds = true;                    
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
        float f00 = TVmap[stage][vLow + 10 * tLow];
        float f10 = TVmap[stage][vLow + 10 * tHigh];
        float f01 = TVmap[stage][vHigh + 10 * tLow];
        float f11 = TVmap[stage][vHigh + 10 * tHigh];
        if (f00 > 0.9 || f10 > 0.9 || f01 > 0.9 || f11 > 0.9) { // Can't include any bad pixels in the interpolation
          return TVmap[stage][vPix + 10 * tPix];
        } else {
          float T0 = -185.0 + 10.0 * tLow;
          float V0 = -45.0 + 10.0 * vLow;
          float VS = (V - V0) / 10.;
          float TS = (T - T0) / 10.;
          float result = f00 + (f10 - f00) * TS + (f01 - f00) * VS + (f11 + f00 - f10 - f01) * TS * VS;
          return result;
        }
      } else {
        return TVmap[stage][vPix + 10 * tPix];
      }
    } else {
      return TVmap[stage][vPix + 10 * tPix];
    }
  }


}; // end of TVcorrection class
#endif
