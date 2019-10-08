#ifndef pedGainCalib_h
#define pedGainCalib_h

// Class to encapsulate all of Vladimir Bashkirov pedestal and gain on-the-fly
// calibration
// R.P. Johnson   5/22/2016

#include <iostream>
#include <cstdio>
#include <string>
#include <ctime>

#include "pCTraw.h"
#include "TVcorrection.h"
#include "Histogram.h"

using namespace std;

class pedGainCalib {

  double Ped[5];
  Histogram *hPed[5];
  Histogram *hEnrg[5];
  Histogram *hEnrgTot;
  Histogram *hTedet;
  ProfilePlot *hProfT;
  FILE *oFile;
  string setTerm;
  char fileName[256];
  char fileName2[256];
  char fileName3[256];
  int threadNumber;
  float emtrng1, emtrng2, emtrng3, emtrng4;
  string partType;

public:
  float corrFac[5]; // Gain correction factor derived here for each stage

  pedGainCalib(string Outputdir, int pdstlr[5], float oldPed[5], int thread,
               float t1, float t2, float t3, float t4, string partType,
               string OsName);

  inline double newPed(int stage) { return Ped[stage]; }

  void rawPh(pCTraw &rawEvt);

  void getPeds(const char *inFileName, int run_number, int program_version,
               float proj_angle, int nKeep, string start_time);

  void weplEvt(float Vedet, float Tedet, float Ene[5]);

  void getGains(TVcorrection *TVcorr, const char *inFileName, int run_number,
                int program_version, int proj_angle, int nKeep,
                string start_time);

}; // end of class pedGainCalib
#endif