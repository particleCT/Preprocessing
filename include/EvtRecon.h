#ifndef _EVTRECON_H_
#define _EVTRECON_H_

// Unpack and reconstruct raw data, and store the results in an event list
#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>
#include "TkrHits.h"
#include "pCT_Tracking.h"
#include "pCTcut.h"
#include "pCTgeo.h"
#include "pedGainCalib.h"
#include "pCTcalib.h"
#include "pCTconfig.h"
#include "TFile.h"

using namespace std;

struct Event{
  float timestamp;
  float Thit[4];
  float Vhit[4];
  int ADC[5];
};

class pCTcalib;
class pCTcut;
class EvtRecon {
 public:
  EvtRecon();
  ~EvtRecon();
  //Functions
  void dumpEvt(Event);
  void ReadInputFile(pCTgeo* Geometry, TVcorrection *const TVcorr, string, pedGainCalib* Calibrate);
  string to_str(int i) { // To fix some stupid compiler problem on my linux box
    long long int j = i;
    return to_string(j);
  }  
  // Variable
  time_t start_time;
  struct tm *now;
  FILE *in_file;
  size_t file_size;
  char outBuff[92];
  int nBuffBytes;
  string evtFileName; // Temporary file for large runs
  FILE *evtFile;
  void writeTmp(Event &evt);
  bool gainAnalysis;
  int Nblk, n_debug, max_time; 
  bool useTmpFile;       
  vector<Event> evtList; // This list doesn't get used if a temporary file is employed instead
  int nEvents;
  float uhitV[4]; // u value at each V layer, assumed to be the same for all events
  float uhitT[4]; // u value at each T layer, assumed to be the same for all events
  int runNumber;
  string runStartTime;
  int study_date;
  float stage_angle;
  int program_version;
  float Peds[5];    // Energy detector pedestals measured from the processed data set
  float GainFac[5]; // Gain correction factors
 private:
  pCTconfig* theConfig;
  
};

#endif
