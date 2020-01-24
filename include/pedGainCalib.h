#ifndef pedGainCalib_h
#define pedGainCalib_h

// Class to encapsulate all of Vladimir Bashkirov pedestal and gain on-the-fly calibration
// R.P. Johnson   5/22/2016

#include <iostream>
#include <cstdio>
#include <string>
#include <ctime>

#include "TFile.h"
#include "TVcorrection.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "pCTraw.h"

using namespace std;

class pedGainCalib {
 public:
  pedGainCalib(TFile* root, int pedMin[5], float oldPed[5], float t1, float t2, float t3, float t4,string partType);
  ~pedGainCalib();  
  //Variables
  double Ped[5];
  TH1D* hPed[5];
  TH1D* hEnrg[5];
  TH1D* hEnrgTot;
  TProfile2D* hProfT;
  TH1D* hTedet;
  TFile* RootFile;
  float emtrng1, emtrng2, emtrng3, emtrng4;
  string partType;
  float GainFac[5]; // Gain correction factor derived here for each stage
  
  // functions
  void ClearHist();
  void FillPeds(pCTraw &rawEvt);
  void GetPeds(const char *inFileName);
  void FillGains(float Vedet, float Tedet, float Ene[5]);
  void GetGains(TVcorrection *TVcorr, const char *inFileName);

}; // end of class pedGainCalib
#endif
