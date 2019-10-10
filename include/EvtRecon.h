#ifndef _EVTRECON_H_
#define _EVTRECON_H_

// Unpack and reconstruct raw data, and store the results in an event list

#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>

#include "TkrHits.h"
#include "pCTgeo.h"
#include "pCT_Tracking.h"
#include "pCTcut.h"
#include "pedGainCalib.h"

using namespace std;

struct Event {
  float Thit[4];
  float Vhit[4];
  int ADC[5];
};

class EvtRecon {

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
  int Nblk;
  string OsName;

public:
  inline ~EvtRecon() { delTmpFile(); }
  EvtRecon(pCTgeo &Geometry, TVcorrection *const TVcorr, string inputFileName, string OutputDir, int max_events,
           int max_time, int n_debug, int n_plot, int nBlocks, bool useTemp, bool doGains, string partType,
           int pdstlr[5], bool reCalibrate, std::string OsName);
  void readTmp(Event &evt);
  void dumpTmp(Event evt);
  void reopenTmpFile();
  void rewindTmpFile();
  string to_str(int i) { // To fix some stupid compiler problem on my linux box
    long long int j = i;
    return to_string(j);
  }

  bool useTmpFile;       // Set true to use temporary files for event storage instead
                         // of local memory (needed for huge runs or little computers)
  vector<Event> evtList; // This list doesn't get used if a temporary file is
                         // employed instead
  int nEvents;
  float uhitV[4]; // u value at each V layer, assumed to be the same for all
                  // events
  float uhitT[4]; // u value at each T layer, assumed to be the same for all
                  // events
  int runNumber;
  string runStartTime;
  int study_date;
  float stage_angle;
  int program_version;
  float Peds[5];    // Energy detector pedestals measured from the processed data
                    // set
  float CorFacs[5]; // Gain correction factors
};

#endif
