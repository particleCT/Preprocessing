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
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "pCTcut.h"
#include "pCTconfig.h"
using namespace std;
#define nEnrg 340
class Wepl {
 public:
  float cut0, cut1, cut2, cut3, cut4, maxEnergy, minEnergy;
  float dEElow[4][3], dEEhigh[4][3];

  // Common parameters
  float RSP;                          // Known phantom RSP
  TFile* projectionROOT;
  TH2D* dEEhist_root[5];
  TGraphErrors* calWEPL[5];
  TGraphErrors* RangeVsEnergy;
  TGraphErrors* EnergyVsRange;
  pCTcut *theCuts;
  void WriteHist(TFile*);
  Wepl(TFile*);       
  ~Wepl();
  // Set energy thresholds in stages
  // Convert energy in stages to WEPL (mm)
  float EtoWEPL(float Estage[5], Int_t &, Int_t &, Int_t &, float &);
 private:
  pCTconfig* theConfig;
  
};


#endif 
