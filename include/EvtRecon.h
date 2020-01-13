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
#include "pctConfig.h"

using namespace std;

struct Event{
  float Thit[4];
  float Vhit[4];
  int ADC[5];
};
class pCTcalib;

class EvtRecon {

 public:

  EvtRecon(pctConfig conf);
  time_t start_time;
  struct tm *now;
  FILE *in_file;
  size_t file_size;
  char outBuff[92];
  int nBuffBytes;
  string evtFileName; // Temporary file for large runs
  FILE *evtFile;
  void writeTmp(Event &evt);
  void delTmpFile();
  bool gainAnalysis;
  int Nblk, n_debug, max_time;
  inline ~EvtRecon() { delTmpFile(); }

  pctConfig config;  
  //Functions   
  void ReadInputFile(pCTgeo* Geometry, TVcorrection *const TVcorr, string);
  void readTmp(Event &evt);
  void dumpTmp(Event evt);
  void reopenTmpFile();
  void rewindTmpFile();
  string to_str(int i) { // To fix some stupid compiler problem on my linux box
    long long int j = i;
    return to_string(j);
    }
  // Variable
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
  float CorFacs[5]; // Gain correction factors
};

#endif
