#ifndef pedGainCalib_h
#define pedGainCalib_h

// Class to encapsulate all of Vladimir Bashkirov pedestal and gain on-the-fly
// calibration
// R.P. Johnson   5/22/2016

#include <iostream>
#include <cstdio>
#include <string>
#include <ctime>

#include "TFile.h"
#include "TVcorrection.h"
#include "TH1D.h"
#include "TProfile.h"
#include "pCTraw.h"

using namespace std;

class pedGainCalib {
 public:
  pedGainCalib(string Outputdir, int pdstlr[5], float oldPed[5], int thread, float t1, float t2, float t3, float t4,
               string partType);
  ~pedGainCalib();  

  //Variables
  double Ped[5];
  TH1D* hPed[5];
  TH1D* hEnrg[5];
  TH1D* hEnrgTot;
  TProfile* hProfT;
  TH1D* hTedet;
  FILE *oFile;
  string setTerm;
  char fileName[256];
  char fileName2[256];
  char fileName3[256];
  int threadNumber;
  float emtrng1, emtrng2, emtrng3, emtrng4;
  string partType;
  float corrFac[5]; // Gain correction factor derived here for each stage
  
  // functions
  void rawPh(pCTraw &rawEvt);
  void getPeds(TFile*, const char *inFileName,int run_number, int program_version, float proj_angle, int nKeep, string start_time);
  void weplEvt(float Vedet, float Tedet, float Ene[5]);
  void getGains(TVcorrection *TVcorr, TFile*, const char *inFileName, int run_number, int program_version, int proj_angle, int nKeep, string start_time);

}; // end of class pedGainCalib
#endif
