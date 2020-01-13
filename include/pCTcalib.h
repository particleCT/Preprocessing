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
#include "TF1.h"
#include "TProfile.h"
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
  inline float Rend(float x) {   return  0.0057 * x * x + 0.2463 * x; }    // -0.366;
  inline float RendHe(float x) { return -0.0057 * x * x + 0.2463 * x; } // -0.366; }
  
  void plot2D(string fn, string T, string TX, string TY, int N, float X[], float Y[], float E[]);
  int TVmapper();
  void enrgDep();
  void writeTVfile();
  bool getLineIntersection(double, double, double, double, double, double, double, double, double &, double &);
  bool EnrgCut(float [5], float, float, float, float);

  // Variables
  vector<string> calFileNames;
  string CalFile;
  float EnergyBinWidth, RangeBinWidth;
  double topW[2]; // min/max t value for range in top of wedge, relative to the center, both sides
  double brickW;  // half width of range in t for brick-only protons
  double emptyW;  // half width of range in t for selecting empty events

  // Other classes
  EvtRecon *procEvt;
  pCTgeo *Geometry;
  TVcorrection *TVcal;
  pCTcut *cuts;
  
  //Variables
  float TVmap[nStage][nPix]; // 5 TV maps with 1cm pixels, 10x38 cm2 area, total 10x38=380 pixels per stage
  float tPlaneShifts[4][4] = {{0}};
  float Est[nStage] = {0};
  float EnS;
  TFile* pCTcalibRootFile = new TFile("pCTcalib.root", "recreate"); // General File for Recalibration
  TH1D *pxHistE_root[nStage][nPix];
  TH1D *pxHistADC_root[nStage][nPix];
  TH1D *stgHistE_root[nStage];
  TH1D *EsumH;
  TProfile *stgEvsT[nStage];
  TProfile *stgEvsV[nStage];

  Histogram *pxHistADC[nStage][nPix];
  Histogram *stgHistE[nStage];
  Histogram *pxHistE[nStage][nPix];

  time_t currentTime;
  struct tm *now;

  double V[2], T[2], Ut[2], Uv[2];
  // Here are a bunch of parameters used to extract the calibration
  float EG4stage[nStage]; // MC derived stage energies, used to calibrate to MeV (CDH setup)
  int k1[nStage];         // Lots of interpolation parameters for cleaning up calibration curves
  int j1[nStage], j2[nStage], j3[nStage], j4[nStage];
  int i1[nStage], i2[nStage], i3;
};

class quadFit { // Fit a set of points to a quadratic function
  double AA, BB, CC;
  double sumX, sumX2, sumX3, sumX4;
  double sumY, sumXY, sumX2Y;
  int nPnt;

  // Gaussian elimination to solve a set of 3 linear equations
  void gauss(double A[3][4]) {
    int n = 3;

    for (int i = 0; i < n; i++) {
      // Search for maximum in this column
      double maxEl = abs(A[i][i]);
      int maxRow = i;
      for (int k = i + 1; k < n; k++) {
        if (abs(A[k][i]) > maxEl) {
          maxEl = abs(A[k][i]);
          maxRow = k;
        }
      }

      // Swap maximum row with current row (column by column)
      for (int k = i; k < n + 1; k++) {
        double tmp = A[maxRow][k];
        A[maxRow][k] = A[i][k];
        A[i][k] = tmp;
      }

      // Make all rows below this one 0 in current column
      for (int k = i + 1; k < n; k++) {
        double c = -A[k][i] / A[i][i];
        for (int j = i; j < n + 1; j++) {
          if (i == j) {
            A[k][j] = 0;
          } else {
            A[k][j] += c * A[i][j];
          }
        }
      }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    double x[3];
    for (int i = n - 1; i >= 0; i--) {
      x[i] = A[i][n] / A[i][i];
      for (int k = i - 1; k >= 0; k--) {
        A[k][n] -= A[k][i] * x[i];
      }
    }
    AA = x[0];
    BB = x[1];
    CC = x[2];
  }

public:
  quadFit() {
    sumX = 0.;
    sumX2 = 0.;
    sumX3 = 0.;
    sumX4 = 0.;
    sumY = 0.;
    sumXY = 0.;
    sumX2Y = 0.;
    nPnt = 0;
    AA = 0.;
    BB = 0.;
    CC = 0.;
  }
  void addPnt(float x, float y) {
    sumX += x;
    sumX2 += x * x;
    sumX3 += x * x * x;
    sumX4 += x * x * x * x;
    sumY += y;
    sumXY += x * y;
    sumX2Y += x * x * y;
    nPnt++;
  }
  int solve() {
    if (nPnt < 3)
      return -1;
    double A[3][4];
    A[0][0] = nPnt;
    A[0][1] = sumX;
    A[0][2] = sumX2;
    A[1][0] = sumX;
    A[1][1] = sumX2;
    A[1][2] = sumX3;
    A[2][0] = sumX2;
    A[2][1] = sumX3;
    A[2][2] = sumX4;
    A[0][3] = sumY;
    A[1][3] = sumXY;
    A[2][3] = sumX2Y;
    this->gauss(A);
    return 0;
  }

  int nPoints() { return nPnt; }

  inline double eval(double x) { return AA + (BB + CC * x) * x; }
};

class qSpline { // Simple quadratic spline extrapolation
  double A, B, C;

public:
  qSpline(float R[nEnrg], int i1, int i2, double x0) {
    double x1 = (double)i1 + 0.5;
    double x2 = (double)i2 + 0.5;
    double y1 = (R[i1 - 1] + R[i1] + R[i1 + 1]) / 3.0;
    double y2 = (R[i2 - 1] + R[i2] + R[i2 + 1]) / 3.0;
    A = (y2 * (x1 - x0) - y1 * (x2 - x0)) / ((x1 - x0) * (x2 * x2 - x0 * x0) - (x2 - x0) * (x1 * x1 - x0 * x0));
    B = (y1 - A * (x1 * x1 - x0 * x0)) / (x1 - x0);
    C = y1 - (A * x1 + B) * x1;
  }
  inline double eval(double x) { return (A * x + B) * x + C; }
};


#endif
