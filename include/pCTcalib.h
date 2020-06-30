#ifndef _pCTcalib_H_
#define _pCTcalib_H_
// Phase-II scanner data preprocessing code.  This class, the first to be run,  carries out the TV calibration task.
// R.P. Johnson  September 18, 2016, adapted from Vladimir Bashkirov's calibration code.

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <ctime>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TSpectrum.h"
#include "pCTgeo.h"
#include "EvtRecon.h"
#include "pCTcut.h"
#include "TGraphErrors.h"
#include "pCTconfig.h"
using namespace std;
#define nStage 5
#define nBricks 5
#define nRange 260
#define nEnrg 340

//#define nRange 520
//#define nEnrg 680
class pCTcut;
class EvtRecon;
class pCTcalib {
 public:
  
  // Functions
  pCTcalib(string inputFileNameIn);
  ~pCTcalib();

  //pCTconfig config;
  int Wcalib();
  void procWEPLcal(TH2D*[nStage], TH2D*[nStage], TH2D*);

  // NIST PSTAR proton range in polystyrene vs E corrected for Birk's law with Kb=.02 cm/meV
  // Experimentally measured R vs E for our 5 stage detector:
  inline float Rend(float x) {   return 0.0057 * x * x + 0.2463 * x; }    // -0.366;
  inline float RendHe(float x) { return -0.0057 * x * x + 0.2463 * x; } // -0.366; }
  
  void plot2D(string fn, string T, string TX, string TY, int N, float X[], float Y[], float E[]);
  int TVmapper();
  int TVmapper_FlatBricks();
  void enrgDep();
  void writeCalibfile();
  bool getLineIntersection(float, float, float, float, float, float, float, float, float &, float &);
  bool EnrgCut(float [5], float, float, float, float);
  void FilldEE(TH2D* dEEhist[nStage]);
  // Other classes
  EvtRecon *theEvtRecon;
  pCTgeo *theGeometry;
  TVcorrection *theTVcorr;
  pCTcut *theCuts;
  pedGainCalib* theCalibration;

  //pedGainCalib *theCalibration;
  // Variables
  vector<string> calFileNames;
  string CalFile;
  float EnergyBinWidth, RangeBinWidth;
  float topW[2]; // min/max t value for range in top of wedge, relative to the center, both sides
  float brickW;  // half width of range in t for brick-only protons
  float emptyW;  // half width of range in t for selecting empty events

  float Est[nStage] = {0}; // energy in each stage
  float EnS; // energy sum in all stage

  float Rst[nStage][nEnrg] = {{0.}}; // 5 arrays to store range vs E
  float Sst[nStage][nEnrg] = {{0.}}; // Peak width for each energy (sigma)
  float TestArray[nEnrg] = {0.};
  float est[nEnrg] = { 0. };         // Corresponding E array in MeV/4

  float dEElow[5][3];
  float dEEhigh[5][3];
    
  TFile* pCTcalibRootFile = new TFile("pCTcalib.root", "recreate"); // General File for Recalibration
  TH1D *pxHistE[nStage][nPix];
  TH1D *pxHistADC[nStage][nPix];
  TH1D *stgHistE[nStage];
  TH1D *EsumH;
  //TH2D* TVcorrHist[5]; // TVcorr histogram for each stage
  //TH2D *stgE[nStage];

  time_t currentTime;
  struct tm *now;
  TProfile* ProfileE[nStage];
  TH2D* ProfileE_Tot;
  TH2D* ProfileE_Tot_f;
  float V[2], T[2], Ut[2], Uv[2], Uft[2], Ufv[2], Tf[2], Vf[2]; 
  // Here are a bunch of parameters used to extract the calibration
  float EG4stage[nStage]; // MC derived stage energies, used to calibrate to MeV (CDH setup)
  float TVnormalizeFactor;
  float Teststage[nStage];
  int k1[nStage];         // Lots of interpolation parameters for cleaning up calibration curves
  int j1[nStage], j2[nStage], j3[nStage], j4[nStage];
  int i1[nStage], i2[nStage], i3;
 private:
  pCTconfig* theConfig;
  
};
#endif
