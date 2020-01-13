/* -----------------------------------------------------------------------------
Wepl.h        Class merged from Wepl1 and Wepl2 [C.E. Ordonez, Aug 2016]
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
#ifndef _WEPL_H_
#define _WEPL_H_
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TH2D.h"
#include "Util.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"
using namespace std;
#define nEnrg 340
class Wepl {
 public:
  string partType;
  bool dodEEFilter;
  float cut0, cut1, cut2, cut3, cut4, cut5;
  float binFac;
  float dEEparams[4][6];

  // Parameters for wedge calibration
  float rvse0[nEnrg], rvse1[nEnrg], rvse2[nEnrg], rvse3[nEnrg], rvse4[nEnrg];

  // Common parameters
  float thr0, thr1, thr2, thr3, thr4; // Stage thresholds to define when a stage has a real signal
  float RSP;                          // Known phantom RSP

  TProfile* calWEPLprofile[5];
  TH2D* dEEhist_root[5];
  // Explicit constructor
  Wepl(const char *WcalibFile, int year, int month, int day, int run, string partType, bool dodEEFilter, string outputDir);
       
  // Explicit destructor
  ~Wepl();
  // Set energy thresholds in stages
  void SetEthresholds1(float, float, float, float, float);

  // Convert energy in stages to WEPL (mm)
  float EtoWEPL(float Estage[5]);
  // some util for splitting key values for the dEE
  vector<string> split(string str, char delimiter);
};

inline Wepl::~Wepl(){
  TFile* Weplroot_file = new TFile("wepl.root", "update");
  Weplroot_file->mkdir("dEE");
  Weplroot_file->cd("dEE");
  for(int i =1;i<5;i++) dEEhist_root[i]->Write("",TObject::kOverwrite);
  for(int i =1;i<5;i++) dEEhist_root[0]->Add(dEEhist_root[i]);
  dEEhist_root[0]->Write("dEEhist_tot",TObject::kOverwrite);

  Weplroot_file->cd();
  Weplroot_file->mkdir("calWEPL");
  Weplroot_file->cd("calWEPL");
  for(int i =0;i<5;i++) calWEPLprofile[i]->Write("",TObject::kOverwrite);
  Weplroot_file->Close();

}  

inline Wepl::Wepl(const char *WCalibFile, int year, int month, int day, int run, string partType, bool dodEEFilter,
                  string outputDir) {

  this->partType = partType;
  this->dodEEFilter = dodEEFilter;
  float scale = 1.0;
  if (partType == "He") {  // sanity cuts
    cut0 = 87. ; cut1 = 95.;  cut2 = 100.;
    cut3 = 120.; cut4 = 300.; cut5 = 270.;
  } else {
    cut0 = 19.; cut1 = 21.; cut2 = 25.;
    cut3 = 35.; cut4 = 75.; cut5 = 77.;
  }

  float EnergyBinWidth = 0.25;
  if (partType == "He") EnergyBinWidth = 1.0;
  binFac = 1. / EnergyBinWidth;
  if (partType == "He") {
    // lower limits
    dEEparams[0][0] = 0.00085; dEEparams[0][1] = -0.5664; dEEparams[0][2] = 220;
    dEEparams[1][0] = 0.00085; dEEparams[1][1] = -0.5664; dEEparams[1][2] = 222;
    dEEparams[2][0] = 0.00085; dEEparams[2][1] = -0.5664; dEEparams[2][2] = 226;
    dEEparams[3][0] = 0.00045; dEEparams[3][1] = -0.42;   dEEparams[3][2] = 237;
    
    // upper limits
    dEEparams[0][3] = 0.00085; dEEparams[0][4] = -0.5664; dEEparams[0][5] = 246;
    dEEparams[1][3] = 0.00085; dEEparams[1][4] = -0.5664; dEEparams[1][5] = 248;
    dEEparams[2][3] = 0.00085; dEEparams[2][4] = -0.5664; dEEparams[2][5] = 254;
    dEEparams[3][3] = 0.00045; dEEparams[3][4] = -0.42;   dEEparams[3][5] = 260;

  } else {
    // lower limits
    dEEparams[0][0] = 0.0040656; dEEparams[0][1] = -0.719491; dEEparams[0][2] = 65.9463;
    dEEparams[1][0] = 0.00378842;dEEparams[1][1] = -0.701085; dEEparams[1][2] = 66.1436;
    dEEparams[2][0] = 0.00380467;dEEparams[2][1] = -0.71088;  dEEparams[2][2] = 67.9228;
    dEEparams[3][0] = 0.0032026; dEEparams[3][1] = -0.663234; dEEparams[3][2] = 69.0353;
    // upper limits
    dEEparams[0][3] = 0.00365163;dEEparams[0][4] = -0.735325; dEEparams[0][5] = 76.3504;
    dEEparams[1][3] = 0.00384016;dEEparams[1][4] = -0.732287; dEEparams[1][5] = 76.5896;
    dEEparams[2][3] = 0.00385641;dEEparams[2][4] = -0.742083; dEEparams[2][5] = 78.3689;
    dEEparams[3][3] = 0.00372006;dEEparams[3][4] = -0.709805; dEEparams[3][5] = 79.5232;
  }

  for(int i =0; i<5; i++){
    dEEhist_root[i] = new TH2D(Form("dE-EStage_%d_bricks",i), Form("dE-E spectra for stage %d ", i), // Name-Title
			       nEnrg, 0., nEnrg * EnergyBinWidth, nEnrg, 0., nEnrg * EnergyBinWidth);} // X-Y binning 
      
  cout << "Wepl::Wepl: Reading WEPL calibration file: " << WCalibFile << " for particle type " << partType << endl;
  FILE *fp = fopen(WCalibFile, "r");
  if (fp == NULL) {
    cerr << "Error opening W correction file " << WCalibFile << endl;
    exit(EXIT_FAILURE);
  }
  // Count number of parameters in file
  float xtemp;
  int nwords = 0;
  int ret;
  while (fscanf(fp, "%f", &xtemp) == 1) nwords++;
  rewind(fp);
  for(int i =0;i<5;i++) calWEPLprofile[i]= new TProfile(Form("calWepl_%d",i), Form("calWEPl_%d",i), 500,0,340);
    string minDate = "";
    string maxDate = "";
    int minRun = 0;
    int maxRun = 0;
    Util util;
    RSP = 1.030; // 1.0300 is for real calibration phantoms
    for (int k = 0; k < nEnrg; ++k) ret = fscanf(fp, "%f", &rvse0[k]);
    for (int k = 0; k < nEnrg; ++k) ret = fscanf(fp, "%f", &rvse1[k]);
    for (int k = 0; k < nEnrg; ++k) ret = fscanf(fp, "%f", &rvse2[k]);
    for (int k = 0; k < nEnrg; ++k) ret = fscanf(fp, "%f", &rvse3[k]);
    for (int k = 0; k < nEnrg; ++k) ret = fscanf(fp, "%f", &rvse4[k]);      
    fclose(fp);
    for (int k = 0; k < nEnrg; k++) {
      float x = (float(k) + 0.5) * EnergyBinWidth;
      calWEPLprofile[0]->Fill(x,rvse0[k]);
      calWEPLprofile[1]->Fill(x,rvse1[k]);
      calWEPLprofile[2]->Fill(x,rvse2[k]);
      calWEPLprofile[3]->Fill(x,rvse3[k]);
      calWEPLprofile[4]->Fill(x,rvse4[k]);
    }

    string line;
    ifstream infile(WCalibFile);
    if (infile) {
      cout << endl << "Echoing comment lines from the WEPL calibration file " << WCalibFile << ": " << endl;
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
            if (key == "minRun") minRun = strtol(value.c_str(), &ePtr, 10);
            if (key == "maxRun") maxRun = strtol(value.c_str(), &ePtr, 10);
	    if (key == "particle") {
              if (value.compare(partType) != 0) {
                cout << "********* WARNING: the WEPL calibration particle type '" << value
                     << "' does not match the beam particle type '" << partType << "'" << endl;
              }
            }
            if (key == "EnergyBinWidth") {
              float enrgBW = stof(value);
              if (abs(enrgBW - EnergyBinWidth) > 0.0001) {
                cout << "*********** WARNING: the WEPL array energy bin width "
                        "in the calibration file does not match the expected value." << endl;
                cout << "            Please check the particle type specification." << endl;		  
              }
            }
            if (key == "Threshold0") {
              float f_value = stof(value);
              if (f_value != thr0) {
                cout << "Wepl.h: we are overwriting the threshold setting " << thr0 << " for stage 0 by the value "
                     << f_value << " from the file " << WCalibFile << endl;
                thr0 = f_value;
              }
            }
            if (key == "Threshold1") {
              float f_value = stof(value);
              if (f_value != thr1) {
                cout << "Wepl.h: we are overwriting the threshold setting " << thr1 << " for stage 1 by the value "
                     << f_value << " from the file " << WCalibFile << endl;
                thr1 = f_value;
              }
            }
            if (key == "Threshold2") {
              float f_value = stof(value);
              if (f_value != thr2) {
                cout << "Wepl.h: we are overwriting the threshold setting " << thr2 << " for stage 2 by the value "
                     << f_value << " from the file " << WCalibFile << endl;
                thr2 = f_value;
              }
            }
            if (key == "Threshold3") {
              float f_value = stof(value);
              if (f_value != thr3) {
                cout << "Wepl.h: we are overwriting the threshold setting " << thr3 << " for stage 3 by the value "
                     << f_value << " from the file " << WCalibFile << endl;
                thr3 = f_value;
              }
            }
            if (key == "Threshold4") {
              float f_value = stof(value);
              if (f_value != thr4) {
                cout << "Wepl.h: we are overwriting the threshold setting " << thr4 << " for stage 4 by the value "
                     << f_value << " from the file " << WCalibFile << endl;
                thr4 = f_value;
              }
            }
	    
            if (dodEEFilter && key == "dEE1") {
              vector<string> sep = split(value, ',');
              if (sep.size() == 6) {
                for (int i = 0; i < 6; i++) dEEparams[0][i] = stof(sep.at(i));
              } else {
                cout << "Need 6 parameters for dEE filtering: the default "
                        "values will be used for stage 1" << endl;
              }
            }
            if (dodEEFilter && key == "dEE2") {
              vector<string> sep = split(value, ',');
              if (sep.size() == 6) {
                for (int i = 0; i < 6; i++) dEEparams[1][i] = stof(sep.at(i));
              } else {
                cout << "Need 6 parameters for dEE filtering: the default "
                        "values will be used for stage 2" << endl;
              }
            }
            if (dodEEFilter && key == "dEE3") {
              vector<string> sep = split(value, ',');
              if (sep.size() == 6) {
                for (int i = 0; i<6; i++) dEEparams[2][i] = stof(sep.at(i));
              } else {
                cout << "Need 6 parameters for dEE filtering: the default "
                        "values will be used for stage 3" << endl;
              }
            }
            if (dodEEFilter && key == "dEE4") {
              vector<string> sep = split(value, ',');
              if (sep.size() == 6) {
                for (int i = 0; i<6;  i++) dEEparams[3][i] = stof(sep.at(i));
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
                      "appropriate for these data with year=" << year << " month= " << month << " day= " << day << " run= " << run << endl;
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
    cout << "The RSP of the WEPL calibration phantom is assumed to be " << RSP << endl;
}
inline float Wepl::EtoWEPL(float Estage[5]) // Return calibrated WEPL from the wedge calibration
{
  int k;
  float e0 = Estage[0]; // These are just to make the values visible in the Windows debugger
  float e1 = Estage[1];
  float e2 = Estage[2];
  float e3 = Estage[3];
  float e4 = Estage[4];
  
  // dE-E parameterization and energy check to cut out fragments for helium
  // if particle stop in Stage 4
  if (Estage[4] > thr4) { // 1 MeV
    dEEhist_root[4]->Fill(Estage[3], Estage[4]);
    if (Estage[4] > cut4) return (1000 + (rvse4[k] * RSP)); // 1000 // Max Trans Filter
    if (Estage[3] < cut3 || Estage[2] < cut2 || Estage[1] < cut1 || Estage[0] < cut0) return (2000 + (rvse4[k] * RSP)); //-1000;  Threshold filter

    if (dodEEFilter) { // dEEFilter
      if (Estage[3] < (dEEparams[3][0] * Estage[4] * Estage[4] + dEEparams[3][1] * Estage[4] + dEEparams[3][2]) ||
          Estage[3] > (dEEparams[3][3] * Estage[4] * Estage[4] + dEEparams[3][4] * Estage[4] + dEEparams[3][5]))
        return (5000 + (rvse4[k] * RSP)); //5000
    }
    k = int(binFac * Estage[4]);
    return rvse4[k] * RSP; // polystyrene to water equivalent // everything passed
  }
  // if particle stop in Stage 3
  else if (Estage[3] > thr3) {
    dEEhist_root[3]->Fill(Estage[2], Estage[3]);
    if (Estage[3] > cut5) return (1000 + (rvse3[k] * RSP)); // 1000 // Max Trans Filter
    if (Estage[2] < cut2 || Estage[1] < cut1 || Estage[0] < cut0) return (2000 + rvse3[k] * RSP); //-1000 Threshold Filter
      
    if (dodEEFilter) { // dEEFilter
      if (Estage[2] < (dEEparams[2][0] * Estage[3] * Estage[3] + dEEparams[2][1] * Estage[3] + dEEparams[2][2]) ||
          Estage[2] > (dEEparams[2][3] * Estage[3] * Estage[3] + dEEparams[2][4] * Estage[3] + dEEparams[2][5]))
        return (5000 + (rvse3[k] * RSP)); //-1000
    }
    k = int(binFac * Estage[3]);
    return rvse3[k] * RSP; // polystyrene to water equivalent
  }
  // if particle stop in Stage 2
  else if (Estage[2] > thr2) { // 1 MeV
    dEEhist_root[2]->Fill(Estage[1], Estage[2]);
    if (Estage[2] > cut5) return (1000 + (rvse2[k] * RSP)); // 1000 // Max Trans Filter
    if (Estage[1] < cut1 || Estage[0] < cut0) return (2000 + (rvse2[k] * RSP)); //-1000; Threshold Filter

    if (dodEEFilter) { // dEEFilter
      if (Estage[1] < (dEEparams[1][0] * Estage[2] * Estage[2] + dEEparams[1][1] * Estage[2] + dEEparams[1][2]) ||
          Estage[1] > (dEEparams[1][3] * Estage[2] * Estage[2] + dEEparams[1][4] * Estage[2] + dEEparams[1][5]))
        return (5000 + (rvse2[k] * RSP)); //-1000;
    }
    k = int(binFac * Estage[2]);
    return rvse2[k] * RSP; // polystyrene to water equivalent
  }
  // if particle stop in Stage 1
  else if (Estage[1] > thr1) {
    dEEhist_root[1]->Fill(Estage[0], Estage[1]);
    if (Estage[1] > cut5) return (1000 + (rvse1[k] * RSP)); // 1000 // 1000 // Max Trans Filter      
    if (Estage[0] < cut0) return (2000 + (rvse1[k] * RSP)); //-1000; Threshold Filter
      
    if (dodEEFilter) { // dEEFilter
      if (Estage[0] < (dEEparams[0][0] * Estage[1] * Estage[1] + dEEparams[0][1] * Estage[1] + dEEparams[0][2]) ||
          Estage[0] > (dEEparams[0][3] * Estage[1] * Estage[1] + dEEparams[0][4] * Estage[1] + dEEparams[0][5]))
        return (5000 + (rvse1[k] * RSP)); //-1000;
    }
    k = int(binFac * Estage[1]);
    return rvse1[k] * RSP; // polystyrene to water equivalent
  }
  // if particle stop in Stage 0
  else if (Estage[0] > thr0) {
    if (Estage[0] > cut5) return (1000 + (rvse0[k] * RSP)); // 1000; // 1000 // Max Trans Filter
    k = int(binFac * Estage[0]);
    return rvse0[k] * RSP; // polystyrene to water equivalent
  }
  else return 2000;   
}

inline void Wepl::SetEthresholds1(float t0, float t1, float t2, float t3, float t4) {
  thr0 = t0;
  thr1 = t1;
  thr2 = t2;
  thr3 = t3;
  thr4 = t4;
  cout << "Wepl::SetEthresholds1: WEPL detector stage energy thresholds are "
          "set to " << thr0 << " " << thr1 << " " << thr2 << " " << thr3 << " " << thr4 << " MeV" << endl;
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

#endif 
