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
#include "Histogram.h"
#include "pCTgeo.h"
#include "EvtRecon.h"
#include "pCTcut.h"
#include "TGraphErrors.h"
#include "pCTconfig.h"
using namespace std;

#define nRange 260
#define nEnrg 340
class pCTcut;
class EvtRecon;
class pCTcalib {
 public:
  
  // Functions
  pCTcalib(pCTconfig cfg, string inputFileNameIn);
  ~pCTcalib();

  pCTconfig config;
  int Wcalib();
  void procWEPLcal(Histogram2D *[nStage], TH2D*[nStage], TH2D*[nStage]);

  // NIST PSTAR proton range in polystyrene vs E corrected for Birk's law with Kb=.02 cm/meV
  // Experimentally measured R vs E for our 5 stage detector:
  inline float Rend(float x) {   return 0.0057 * x * x + 0.2463 * x; }    // -0.366;
  inline float RendHe(float x) { return -0.0057 * x * x + 0.2463 * x; } // -0.366; }
  
  void plot2D(string fn, string T, string TX, string TY, int N, float X[], float Y[], float E[]);
  int TVmapper();
  void enrgDep();
  void writeCalibfile();
  bool getLineIntersection(double, double, double, double, double, double, double, double, double &, double &);
  bool EnrgCut(float [5], float, float, float, float);

  // Other classes
  EvtRecon *theEvtRecon;
  pCTgeo *theGeometry;
  TVcorrection *theTVcorr;
  pCTcut *theCuts;
  //pedGainCalib *theCalibration;
  // Variables
  vector<string> calFileNames;
  string CalFile;
  float EnergyBinWidth, RangeBinWidth;
  double topW[2]; // min/max t value for range in top of wedge, relative to the center, both sides
  double brickW;  // half width of range in t for brick-only protons
  double emptyW;  // half width of range in t for selecting empty events

  float Est[nStage] = {0}; // energy in each stage
  float EnS; // energy sum in all stage

  float Rst[nStage][nEnrg] = {{0.}}; // 5 arrays to store range vs E
  float Sst[nStage][nEnrg] = {{0.}}; // Peak width for each energy (sigma)
  float est[nEnrg] = { 0. };         // Corresponding E array in MeV/4

  float dEElow[5][3];
  float dEEhigh[5][3];
    
  TFile* pCTcalibRootFile = new TFile("pCTcalib.root", "recreate"); // General File for Recalibration
  TH1D *pxHistE_root[nStage][nPixRoot];
  TH1D *pxHistADC_root[nStage][nPixRoot];
  TH1D *stgHistE_root[nStage];
  TH1D *EsumH;
  TH2D* TVcorrHist[5]; // TVcorr histogram for each stage
  TProfile2D *stgE[nStage];
  Histogram *pxHistADC[nStage][nPix];
  Histogram *stgHistE[nStage];
  Histogram *pxHistE[nStage][nPix];

  time_t currentTime;
  struct tm *now;

  double V[2], T[2], Ut[2], Uv[2];
  // Here are a bunch of parameters used to extract the calibration
  float EG4stage[nStage]; // MC derived stage energies, used to calibrate to MeV (CDH setup)
  float Teststage[nStage];
  int k1[nStage];         // Lots of interpolation parameters for cleaning up calibration curves
  int j1[nStage], j2[nStage], j3[nStage], j4[nStage];
  int i1[nStage], i2[nStage], i3;
};
#endif
