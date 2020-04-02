#ifndef _PREPROCESSING_H_
#define _PREPROCESSING_H_
// Phase-II scanner data preprocessing code, produces the input needed by image
// reconstruction.
// This code also serves for monitoring raw data during runs.
// R.P. Johnson  May 8, 2016
//
// The following classes are used:
// - Preprocessing.h and Preprocessing.cpp: the top level driving program,
// called either by main here or from within another program
// - pCTgeo.h: a geometry package to encapsulate all of the geometry constants
// - Wepl1.h: WEPL calibration class, from V. Bashkirov
// - pedGainCalib.cpp and pedGainCalib.h: on-the-fly adjustment of pedestals and
// gains
// - TVcorrection.h: class for TV corrections of energy detector data
// - pCTraw.cpp and pCTraw.h: raw data input and unpacking class, encapsulates
// the original code of Piersimoni et al.
// - TkrHits.cpp and TkrHits.h: class for calculation of tracker coordinates
// from strip clusters
// - pCT_Tracking.cpp and pCT_Tracking.h: pattern recognition for the tracking;
// includes printing and plotting
#include "pCTconfig.h"
#include "TVcorrection.h"
#include "pedGainCalib.h"
#include "pCTgeo.h"
#include "TkrHits.h"
#include "pCT_Tracking.h"
#include "pCTraw.h"
#include "pCTcut.h"
#include "TFile.h"
class Preprocessing { // Top level program from the pCT preprocessing task.
 public:
  Preprocessing();
  TFile* projectionROOT;
  TFile* pCTcalibRootFile;
  int ProcessFile(float, int, int);
  int ret;
  size_t file_size;
  float Version, beamEnergy, StgThr[5], initialAngle, proj_angle;
  int fileBins, analysisLevel, max_events, max_time, n_debug, n_plot;
  bool callUser, continuous_scan, eventOrder, energyOutput, timeStampOutput, eventIDOutput, dodEEFilter;
  double Uhit[4];
  std::string study_name, Outputdir;
  std::string WcalibFile;
  std::string TVcorrFile;
  TVcorrection* theTVcorr;
  FILE *in_file;
  char inFileName[256];
  time_t start_time;
  pCTcut* theCuts;
  int ADC[5];
  struct tm *now;
  static int findEvt(FILE *fp);
  void pCTevents(pCTgeo* Geometry, pCTraw rawEvt, pedGainCalib *Calibrate, double Uhit[]);
 private:
  pCTconfig* theConfig;
  
};
#endif
