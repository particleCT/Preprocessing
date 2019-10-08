/* -----------------------------------------------------------------------------

Wepl.h        Class merged from Wepl1 and Wepl2 [C.E. Ordonez, Aug 2016]

New Data Member:

bool wedge_calibration    Flag to indicate whether to use wedge calibration or
        step calibration. The flag is set based on the length -- or the
        number of fields found in the WEPL calibration file. For now,
        there are 1700 parameters for wedge calibration, but only 445
        for step calibration.
Methods:

Wepl    Constructor, reads calibration data from file previously
        prepared by one of the following:
        PreProcessingWcalib: step calibration, scans before Aug 2016
        PreProcessingWcalibW: wedge calibration, scans starting Aug 2016
        The filenames are no longer restricted to Wcalib.txt and
        WcalibW.txt; nor the calibration files themselves need to be in
        the same output directory for the projection files.

float EtoWEPL    Converts energy (in Mev) deposited in 5 stages to WEPL (in mm)

----------------------------------------------------------------------------- */

// Class Wepl2 constructor Wepl2  which reads WcalibW.txt file
// with calibration data prepared by PreProcessingWcalibW.exe according
// to the new ("Wedge") calibration scheme.
// Function EtoWepl2 converts  energy deposited in stages  (MeV)
// to WEPL (in mm) using look-up tables.
// Returned WEPL equal to -1000   means that energy deposition
// pattern is inconsistent with Bragg curve for stages 1-4.
// Returned WEPL equal to 1000 means that energy deposition in a
// stage is well above maximum for a single proton event.
// WEPL = 2000 if enrgy depositions in all stages are below thresholds.
// Function SetEthresholds is used to set energy thresholds for 5 stages,
// Created by V.A.Bashkirov, last revision  Aug 2016
//
// Added plotting of calibration curves and all of the comment lines at the end
// of calibration files, R. Johnson
// Added automatic dEEFiltter, L. Volz (Nov, 2018)

#ifndef _WEPL_H_
#define _WEPL_H_
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "Util.h"

using namespace std;

class Wepl {

#define nEnrg 340

  bool wedge_calibration;
  string partType;
  bool dodEEFilter;
  float cut0;
  float cut1;
  float cut2;
  float cut3;
  float cut4;
  float cut5;
  float binFac;
  float dEEparams[4][6];

  // Parameters for stairs calibration
  int Np0, Np1, Np2, Np3, Np4, Npc0, Npc1, Npc2, Npc3;
  float e0[12], e1[12], e2[12], e3[12], e4[12];
  float wetd0[12], wetd1[12], wetd2[12], wetd3[12], wetd4[12];
  float ec0[12], ec1[12], ec2[12], ec3[12];
  float wetd0c[12], wetd1c[12], wetd2c[12], wetd3c[12];
  float a0[10][3], a1[10][3], a2[10][3], a3[10][3], a4[10][3];
  float ac0[10][3], ac1[10][3], ac2[10][3], ac3[10][3];

  // Parameters for wedge calibration
  float rvse0[nEnrg], rvse1[nEnrg], rvse2[nEnrg], rvse3[nEnrg], rvse4[nEnrg];

  // Common parameters
  float thr0, thr1, thr2, thr3,
      thr4;  // Stage thresholds to define when a stage has a real signal
  float RSP; // Known phantom RSP

public:
  // Explicit constructor
  Wepl(const char *WcalibFile, int year, int month, int day,
       int run, string partType, bool dodEEFilter, string outputDir);

  // Set energy thresholds in stages
  void SetEthresholds1(float, float, float, float, float);

  // Convert energy in stages to WEPL (mm)
  float EtoWEPL(float Estage[5]);
  float EtoWEPL1(float Estage[5]); // use step calibration
  float EtoWEPL2(float Estage[5]); // use wedge calibration

  //  private:
  //    std::fstream Wcalfile;
  // some util for splitting key values for the dEE
  vector<string> split(string str, char delimiter);
};

inline Wepl::Wepl(const char *WCalibFile, int year, int month,
                  int day, int run, string partType, bool dodEEFilter,
                  string outputDir) {
  this->partType = partType;
  this->dodEEFilter = dodEEFilter;

  float scale = 1.0;
  if (partType == "He") {
    cut0 = 87.;
    cut1 = 95.;
    cut2 = 100.;
    cut3 = 120.;
    cut4 = 300.;
    cut5 = 270.;
  } else {
    cut0 = 19.;
    cut1 = 21.;
    cut2 = 25.;
    cut3 = 35.;
    cut4 = 75.;
    cut5 = 77.;
  }
  float EnergyBinWidth = 0.25;
  if (partType == "He")
    EnergyBinWidth = 1.0;
  binFac = 1. / EnergyBinWidth;

  if (partType == "He") {
    // lower limits
    dEEparams[0][0] = 0.00085;
    dEEparams[0][1] = -0.5664;
    dEEparams[0][2] = 220;
    dEEparams[1][0] = 0.00085;
    dEEparams[1][1] = -0.5664;
    dEEparams[1][2] = 222;
    dEEparams[2][0] = 0.00085;
    dEEparams[2][1] = -0.5664;
    dEEparams[2][2] = 226;
    dEEparams[3][0] = 0.00045;
    dEEparams[3][1] = -0.42;
    dEEparams[3][2] = 260;
    // upper limits
    dEEparams[0][3] = 0.00085;
    dEEparams[0][4] = -0.5664;
    dEEparams[0][5] = 246;
    dEEparams[1][3] = 0.00085;
    dEEparams[1][4] = -0.5664;
    dEEparams[1][5] = 248;
    dEEparams[2][3] = 0.00085;
    dEEparams[2][4] = -0.5664;
    dEEparams[2][5] = 254;
    dEEparams[3][3] = 0.00045;
    dEEparams[3][4] = -0.42;
    dEEparams[3][5] = 237;

  } else {
    // lower limits
    dEEparams[0][0] = 0.0040656;
    dEEparams[0][1] = -0.719491;
    dEEparams[0][2] = 65.9463;
    dEEparams[1][0] = 0.00378842;
    dEEparams[1][1] = -0.701085;
    dEEparams[1][2] = 66.1436;
    dEEparams[2][0] = 0.00380467;
    dEEparams[2][1] = -0.71088;
    dEEparams[2][2] = 67.9228;
    dEEparams[3][0] = 0.0032026;
    dEEparams[3][1] = -0.663234;
    dEEparams[3][2] = 69.0353;
    // upper limits
    dEEparams[0][3] = 0.00365163;
    dEEparams[0][4] = -0.735325;
    dEEparams[0][5] = 76.3504;
    dEEparams[1][3] = 0.00384016;
    dEEparams[1][4] = -0.732287;
    dEEparams[1][5] = 76.5896;
    dEEparams[2][3] = 0.00385641;
    dEEparams[2][4] = -0.742083;
    dEEparams[2][5] = 78.3689;
    dEEparams[3][3] = 0.00372006;
    dEEparams[3][4] = -0.709805;
    dEEparams[3][5] = 79.5232;
  }

  cout << "Wepl::Wepl: Reading WEPL calibration file: " << WCalibFile
       << " for particle type " << partType << endl;
  FILE *fp = fopen(WCalibFile, "r");
  if (fp == NULL) {
    cerr << "Error opening W correction file " << WCalibFile << endl;
    exit(EXIT_FAILURE);
  }

  // Count number of parameters in file
  float xtemp;
  int nwords = 0;
  while (fscanf(fp, "%f", &xtemp) == 1)
    nwords++;

  // As of Aug 2016: stairs calibration uses 445 parameters
  //               wedge calibration uses 1700 parameters

  rewind(fp);
  if (nwords >= 5 * nEnrg) {
    wedge_calibration = true;
    cout << "Wepl::Wepl: the calibration phantom is assumed to be the wedge."
         << endl;
    string minDate = "";
    string maxDate = "";
    int minRun = 0;
    int maxRun = 0;
    Util util;
    RSP = 1.030; // 1.0300 is for real calibration phantoms
    for (int k = 0; k < nEnrg; ++k)
      fscanf(fp, "%f", &rvse0[k]);
    for (int k = 0; k < nEnrg; ++k)
      fscanf(fp, "%f", &rvse1[k]);
    for (int k = 0; k < nEnrg; ++k)
      fscanf(fp, "%f", &rvse2[k]);
    for (int k = 0; k < nEnrg; ++k)
      fscanf(fp, "%f", &rvse3[k]);
    for (int k = 0; k < nEnrg; ++k)
      fscanf(fp, "%f", &rvse4[k]);
    fclose(fp);

    string fileName =
        outputDir + "/CalWepl.gp"; // To plot the calibration curves
    FILE *oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile,
              "#*** This file is intended to be displayed by gnuplot.\n");
      fprintf(oFile, "#*** Either double click on the file (works in Windows "
                     "at least),\n");
      fprintf(oFile, "#*** or else start up gnuplot and use the load command "
                     "to display the plot.\n");
      fprintf(oFile, "#set term X11 persist\n");
      fprintf(oFile, "set xtics font 'Verdana,12'\n");
      fprintf(oFile, "set ytics font 'Verdana,12'\n");
      fprintf(oFile, "set title 'WEPL Calibration Curves' font 'Verdana,12'\n");
      fprintf(oFile, "set xlabel 'Energy' font 'Verdana,12'\n");
      fprintf(oFile, "set ylabel 'WEPL' font 'Verdana,12'\n");
      fprintf(oFile, "set xrange[0 : %7.4f]\n", EnergyBinWidth * nEnrg);
      fprintf(oFile, "set nokey\n");
      fprintf(oFile, "$stg0 << EOD\n");
      for (int k = 0; k < nEnrg; k++) {
        float x = (float(k) + 0.5) * EnergyBinWidth;
        fprintf(oFile, "%7.4f %7.4f\n", x, rvse0[k]);
      }
      fprintf(oFile, "EOD\n");
      fprintf(oFile, "$stg1 << EOD\n");
      for (int k = 0; k < nEnrg; k++) {
        float x = (float(k) + 0.5) * EnergyBinWidth;
        fprintf(oFile, "%7.4f %7.4f\n", x, rvse1[k]);
      }
      fprintf(oFile, "EOD\n");
      fprintf(oFile, "$stg2 << EOD\n");
      for (int k = 0; k < nEnrg; k++) {
        float x = (float(k) + 0.5) * EnergyBinWidth;
        fprintf(oFile, "%7.4f %7.4f\n", x, rvse2[k]);
      }
      fprintf(oFile, "EOD\n");
      fprintf(oFile, "$stg3 << EOD\n");
      for (int k = 0; k < nEnrg; k++) {
        float x = (float(k) + 0.5) * EnergyBinWidth;
        fprintf(oFile, "%7.4f %7.4f\n", x, rvse3[k]);
      }
      fprintf(oFile, "EOD\n");
      fprintf(oFile, "$stg4 << EOD\n");
      for (int k = 0; k < nEnrg; k++) {
        float x = (float(k) + 0.5) * EnergyBinWidth;
        fprintf(oFile, "%7.4f %7.4f\n", x, rvse4[k]);
      }
      fprintf(oFile, "EOD\n");
      fprintf(oFile, "plot $stg0 with lines lw 3, $stg1 with lines lw 3, $stg2 "
                     "with lines lw 3, $stg3 with lines lw 3, $stg4 with lines "
                     "lw 3\n");
      fprintf(oFile, "show label\n");
      fclose(oFile);
    }

    string line;
    ifstream infile(WCalibFile);
    if (infile) {
      cout << endl << "Echoing comment lines from the WEPL calibration file "
           << WCalibFile << ": " << endl;
      while (getline(infile, line)) {
        size_t found = line.find_first_not_of(" ");
        if (found == line.npos)
          continue;
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
            if (key == "particle") {
              if (value.compare(partType) != 0) {
                cout
                    << "********* WARNING: the WEPL calibration particle type '"
                    << value << "' does not match the beam particle type '"
                    << partType << "'" << endl;
              }
            }
            if (key == "EnergyBinWidth") {
              float enrgBW = stof(value);
              if (abs(enrgBW - EnergyBinWidth) > 0.0001) {
                cout << "*********** WARNING: the WEPL array energy bin width "
                        "in the calibration file does not match the expected "
                        "value." << endl;
                cout << "            Please check the particle type "
                        "specification." << endl;
              }
            }
            if (key == "Threshold0") {
              float f_value = stof(value);
              if (f_value != thr0) {
                cout << "Wepl.h: we are overwriting the threshold setting "
                     << thr0 << " for stage 0 by the value " << f_value
                     << " from the file " << WCalibFile << endl;
                thr0 = f_value;
              }
            }
            if (key == "Threshold1") {
              float f_value = stof(value);
              if (f_value != thr1) {
                cout << "Wepl.h: we are overwriting the threshold setting "
                     << thr1 << " for stage 1 by the value " << f_value
                     << " from the file " << WCalibFile << endl;
                thr1 = f_value;
              }
            }
            if (key == "Threshold2") {
              float f_value = stof(value);
              if (f_value != thr2) {
                cout << "Wepl.h: we are overwriting the threshold setting "
                     << thr2 << " for stage 2 by the value " << f_value
                     << " from the file " << WCalibFile << endl;
                thr2 = f_value;
              }
            }
            if (key == "Threshold3") {
              float f_value = stof(value);
              if (f_value != thr3) {
                cout << "Wepl.h: we are overwriting the threshold setting "
                     << thr3 << " for stage 3 by the value " << f_value
                     << " from the file " << WCalibFile << endl;
                thr3 = f_value;
              }
            }
            if (key == "Threshold4") {
              float f_value = stof(value);
              if (f_value != thr4) {
                cout << "Wepl.h: we are overwriting the threshold setting "
                     << thr4 << " for stage 4 by the value " << f_value
                     << " from the file " << WCalibFile << endl;
                thr4 = f_value;
              }
            }
            if (dodEEFilter && key == "dEE1") {
              vector<string> sep = split(value, ',');
              if (sep.size() == 6) {
                for (int i = 0; i++; i < 6)
                  dEEparams[0][i] = stof(sep.at(i));
              } else {
                cout << "Need 6 parameters for dEE filtering: the default "
                        "values will be used for stage 1" << endl;
              }
            }
            if (dodEEFilter && key == "dEE2") {
              vector<string> sep = split(value, ',');
              if (sep.size() == 6) {
                for (int i = 0; i++; i < 6)
                  dEEparams[1][i] = stof(sep.at(i));
              } else {
                cout << "Need 6 parameters for dEE filtering: the default "
                        "values will be used for stage 2" << endl;
              }
            }
            if (dodEEFilter && key == "dEE3") {
              vector<string> sep = split(value, ',');
              if (sep.size() == 6) {
                for (int i = 0; i++; i < 6)
                  dEEparams[2][i] = stof(sep.at(i));
              } else {
                cout << "Need 6 parameters for dEE filtering: the default "
                        "values will be used for stage 3" << endl;
              }
            }
            if (dodEEFilter && key == "dEE4") {
              vector<string> sep = split(value, ',');
              if (sep.size() == 6) {
                for (int i = 0; i++; i < 6)
                  dEEparams[3][i] = stof(sep.at(i));
              } else {
                cout << "Need 6 parameters for dEE filtering: the default "
                        "values will be used for stage 4" << endl;
              }
            }
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
              cout << "WEPL: the run and date of the data match the "
                      "calibration file ranges.  Success!" << endl;
              checked = true;
            } else {
              cout << "WEPL: the WEPL calibration file date and run are not "
                      "appropriate for these data with year=" << year
                   << " month= " << month << " day= " << day << " run= " << run
                   << endl;
              exit(-1);
            }
          }
        }
      }
    }
    cout << "Completed loading of the WEPL calibration constants." << endl;
    if (!checked) {
      cout << "********** WARNING!  No check has been made of the validity of "
              "the date and run ranges of the WEPL calibration constants! "
              "**************" << endl;
    }
  } else {
    wedge_calibration = false;
    RSP = 1.031;
    fscanf(fp, "%d %d %d %d %d %d %d %d %d", &Np0, &Np1, &Np2, &Np3, &Np4,
           &Npc0, &Npc1, &Npc2, &Npc3);
    for (int i = 0; i < Np0; ++i)
      fscanf(fp, "%f %f", &e0[i], &wetd0[i]);
    for (int i = 0; i < Np1; ++i)
      fscanf(fp, "%f %f", &e1[i], &wetd1[i]);
    for (int i = 0; i < Np2; ++i)
      fscanf(fp, "%f %f", &e2[i], &wetd2[i]);
    for (int i = 0; i < Np3; ++i)
      fscanf(fp, "%f %f", &e3[i], &wetd3[i]);
    for (int i = 0; i < Np4; ++i)
      fscanf(fp, "%f %f", &e4[i], &wetd4[i]);
    for (int i = 0; i < Npc0; ++i)
      fscanf(fp, "%f %f", &ec0[i], &wetd0c[i]);
    for (int i = 0; i < Npc1; ++i)
      fscanf(fp, "%f %f", &ec1[i], &wetd1c[i]);
    for (int i = 0; i < Npc2; ++i)
      fscanf(fp, "%f %f", &ec2[i], &wetd2c[i]);
    for (int i = 0; i < Npc3; ++i)
      fscanf(fp, "%f %f", &ec3[i], &wetd3c[i]);
    for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 10; ++k)
        fscanf(fp, "%f", &a0[k][i]);
      for (int k = 0; k < 10; ++k)
        fscanf(fp, "%f", &a1[k][i]);
      for (int k = 0; k < 10; ++k)
        fscanf(fp, "%f", &a2[k][i]);
      for (int k = 0; k < 10; ++k)
        fscanf(fp, "%f", &a3[k][i]);
      for (int k = 0; k < 10; ++k)
        fscanf(fp, "%f", &a4[k][i]);
      for (int k = 0; k < 10; ++k)
        fscanf(fp, "%f", &ac0[k][i]);
      for (int k = 0; k < 10; ++k)
        fscanf(fp, "%f", &ac1[k][i]);
      for (int k = 0; k < 10; ++k)
        fscanf(fp, "%f", &ac2[k][i]);
      for (int k = 0; k < 10; ++k)
        fscanf(fp, "%f", &ac3[k][i]);
    }
    cout << "From WEPL stairs calibration file: " << Np0 << " " << Np1 << " "
         << Np2 << " " << Np3 << " " << Np4 << " " << Npc0 << " " << Npc1 << " "
         << Npc2 << " " << Npc3 << endl;
    fclose(fp);
  }
  if (wedge_calibration) {
    cout << "Based on the number of parameters in the WEPL calibration file, "
            "the wedge calibration method is assumed." << endl;
  } else {
    cout << "Based on the number of parameters in the WEPL calibration file, "
            "the stairs calibration method is assumed." << endl;
  }
  cout << "The RSP of the WEPL calibration phantom is assumed to be " << RSP
       << endl;
}

inline float Wepl::EtoWEPL(float Estage[5]) {
  return (wedge_calibration ? EtoWEPL2(Estage) : EtoWEPL1(Estage));
}

inline float Wepl::EtoWEPL1(float Estage[5]) // Return calibrated WEPL from the
                                             // stairs calibration.
{
  float wet_wm, wet_pr, x;
  int k;

  if (Estage[4] > thr4) {
    if (Estage[4] > 80.)
      return 1000;
    // Check if energy deposition pattern is compatible to Bragg curve (within 5
    // sigma):
    if (Estage[3] < 35 || Estage[2] < 25 || Estage[1] < 21 || Estage[0] < 19)
      return -1000;
    k = Np4 - 2;
    x = Estage[4];
    for (int i = 1; i < Np4 / 2; i++) {
      if (x < e4[2 * i]) {
        k = 2 * i - 1;
        break;
      }
    }
    wet_wm = a4[k][0] * x * x + a4[k][1] * x + a4[k][2];
    if (x < 3 && (Estage[3] < 70.5 - x || Estage[3] > 81 - x))
      return -1000;
    //      if(Estage[3]>70 && x<1.5) goto St3;
    //      if(Estage[3]<60 && x<4) return -1000;
    if (x < .12) {
      if (Estage[3] > 70 && x < 3)
        goto St3;
      k = Npc3 - 2;
      x = Estage[3] + Estage[2] + Estage[1] + Estage[0];
      for (int i = 1; i < Npc3 / 2; i++) {
        if (x < ec3[2 * i]) {
          k = 2 * i - 1;
          break;
        }
      }
      wet_pr = ac3[k][0] * x * x + ac3[k][1] * x + ac3[k][2];
      if (fabs(wet_pr - wet_wm) > 10)
        return -1000;
      //       if(wet_pr>wet_wm+3) wet_wm=wet_pr;
      //    if(fabs(wet_wm-wet_pr)>3) wet_wm=.5*wet_wm+.5*wet_pr;
    }
    return wet_wm * RSP; // polystyrene to water equivalent
  }

  if (Estage[3] > thr3) {
    if (Estage[3] > 80.)
      return 1000;
    if (Estage[2] < 25 || Estage[1] < 21 || Estage[0] < 19)
      return -1000;

  St3:
    k = Np3 - 2;
    x = Estage[3];
    for (int i = 1; i < Np3 / 2; i++) {
      if (x < e3[2 * i]) {
        k = 2 * i - 1;
        break;
      }
    }
    wet_wm = a3[k][0] * x * x + a3[k][1] * x + a3[k][2];
    if (x < 3 && (Estage[2] < 70.0 - x || Estage[2] > 80 - x))
      return -1000;
    return wet_wm * RSP;
  }

  if (Estage[2] > thr2) {
    if (Estage[2] > 79.)
      return 1000;
    if (Estage[1] < 21 || Estage[0] < 19)
      return -1000;
    k = Np2 - 2;
    x = Estage[2];
    for (int i = 1; i < Np2 / 2; i++) {
      if (x < e2[2 * i]) {
        k = 2 * i - 1;
        break;
      }
    }
    wet_wm = a2[k][0] * x * x + a2[k][1] * x + a2[k][2];
    if (x < 3 && (Estage[1] < 68.0 - x || Estage[1] > 78 - x))
      return -1000;
    return wet_wm * RSP;
  }

  if (Estage[1] > thr1) {
    if (Estage[1] > 78.)
      return 1000;
    if (Estage[0] < 19)
      return -1000;
    k = Np1 - 2;
    x = Estage[1];
    for (int i = 1; i < Np1 / 2; i++) {
      if (x < e1[2 * i]) {
        k = 2 * i - 1;
        break;
      }
    }
    wet_wm = a1[k][0] * x * x + a1[k][1] * x + a1[k][2];
    if (x < 3 && (Estage[0] < 67.0 - x || Estage[0] > 77 - x))
      return -1000;
    return wet_wm * RSP;
  }

  if (Estage[0] > thr0) {
    if (Estage[0] > 80.)
      return 1000;
    k = Np0 - 2;
    x = Estage[0];
    for (int i = 1; i < Np0 / 2; i++) {
      if (x < e0[2 * i]) {
        k = 2 * i - 1;
        break;
      }
    }
    wet_wm = a0[k][0] * x * x + a0[k][1] * x + a0[k][2];
    return wet_wm * RSP;
  }

  return 2000;
}

inline float Wepl::EtoWEPL2(float Estage[5]) // Return calibrated WEPL from the
                                             // wedge calibration
{
  int k;
  float e0 = Estage[0]; // These are just to make the values visible in the
                        // Windows debugger
  float e1 = Estage[1];
  float e2 = Estage[2];
  float e3 = Estage[3];
  float e4 = Estage[4];

  //
  // dE-E parameterization and energy check to cut out fragments for helium
  //

  if (Estage[4] > thr4) {
    if (Estage[4] > cut4)
      return (1000 + (rvse4[k] * RSP)); // 1000
    if (Estage[3] < cut3 || Estage[2] < cut2 || Estage[1] < cut1 ||
        Estage[0] < cut0)
      return (2000 + (rvse4[k] * RSP)); //-1000;
    if (dodEEFilter) {
      if (Estage[3] < (dEEparams[3][0] * Estage[4] * Estage[4] +
                       dEEparams[3][1] * Estage[4] + dEEparams[3][2]) ||
          Estage[3] > (dEEparams[3][3] * Estage[4] * Estage[4] +
                       dEEparams[3][4] * Estage[4] + dEEparams[3][5]))
        return (5000 + (rvse4[k] * RSP)); //-1000
    }
    k = int(binFac * Estage[4]);
    return rvse4[k] * RSP; // polystyrene to water equivalent
  } else if (Estage[3] > thr3) {
    if (Estage[3] > cut5)
      return (1000 + (rvse3[k] * RSP)); // 1000
    if (Estage[2] < cut2 || Estage[1] < cut1 || Estage[0] < cut0)
      return (2000 + rvse3[k] * RSP); //-1000
    if (dodEEFilter) {
      if (Estage[2] < (dEEparams[2][0] * Estage[3] * Estage[3] +
                       dEEparams[2][1] * Estage[3] + dEEparams[2][2]) ||
          Estage[2] > (dEEparams[2][3] * Estage[3] * Estage[3] +
                       dEEparams[2][4] * Estage[3] + dEEparams[2][5]))
        return (5000 + (rvse3[k] * RSP)); //-1000
    }
    k = int(binFac * Estage[3]);
    return rvse3[k] * RSP; // polystyrene to water equivalent
  } else if (Estage[2] > thr2) {
    if (Estage[2] > cut5)
      return (1000 + (rvse2[k] * RSP)); // 1000
    if (Estage[1] < cut1 || Estage[0] < cut0)
      return (2000 + (rvse2[k] * RSP)); //-1000;
    if (dodEEFilter) {
      if (Estage[1] < (dEEparams[1][0] * Estage[2] * Estage[2] +
                       dEEparams[1][1] * Estage[2] + dEEparams[1][2]) ||
          Estage[1] > (dEEparams[1][3] * Estage[2] * Estage[2] +
                       dEEparams[1][4] * Estage[2] + dEEparams[1][5]))
        return (5000 + (rvse2[k] * RSP)); //-1000;
    }
    k = int(binFac * Estage[2]);
    return rvse2[k] * RSP; // polystyrene to water equivalent
  } else if (Estage[1] > thr1) {
    if (Estage[1] > cut5)
      return (1000 + (rvse1[k] * RSP)); // 1000
    if (Estage[0] < cut0)
      return (2000 + (rvse1[k] * RSP)); //-1000;
    if (dodEEFilter) {
      if (Estage[0] < (dEEparams[0][0] * Estage[1] * Estage[1] +
                       dEEparams[0][1] * Estage[1] + dEEparams[0][2]) ||
          Estage[0] > (dEEparams[0][3] * Estage[1] * Estage[1] +
                       dEEparams[0][4] * Estage[1] + dEEparams[0][5]))
        return (5000 + (rvse1[k] * RSP)); //-1000;
    }
    k = int(binFac * Estage[1]);
    return rvse1[k] * RSP; // polystyrene to water equivalent
  } else if (Estage[0] > thr0) {
    if (Estage[0] > cut5)
      return (1000 + (rvse0[k] * RSP)); // 1000;
    k = int(binFac * Estage[0]);
    return rvse0[k] * RSP; // polystyrene to water equivalent
  } else
    return 2000;
}

inline void Wepl::SetEthresholds1(float t0, float t1, float t2, float t3,
                                  float t4) {
  thr0 = t0;
  thr1 = t1;
  thr2 = t2;
  thr3 = t3;
  thr4 = t4;
  cout << "Wepl::SetEthresholds1: WEPL detector stage energy thresholds are "
          "set to " << thr0 << " " << thr1 << " " << thr2 << " " << thr3 << " "
       << thr4 << " MeV" << endl;
}

vector<string> Wepl::split(string str, char delimiter) {
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;

  while (getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }

  return internal;
}

#endif // #ifndef Wepl2_h
