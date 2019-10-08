#ifndef UserAnalysis_h
#define UserAnalysis_h

// The UserAnalysis class is defined to provide entry points for the user's
// private analysis code

// This code can slow down the processing significantly, so it is called only if
// the -u option is given.

// This is just an example user routine, which corresponds to the monitoring
// code used to verify
// the scanner operation during runs.  You may delete it all and/or replace with
// whatever you want,
// because nothing here is needed in order to produce the output file for image
// reconstruction.
// Please put your private code for peaking at the data here, instead of
// cluttering up the
// public code needed by all to preprocess or calibrate data.

// Note that the UserAnalysis entry points are called for each event by only the
// mother thread, so if
// you run your preprocessing job in multiple threads (-t option) only a
// fraction of the events will
// get analyzed in UserAnalysis (but all will be processed and written out to
// the binary file).

using namespace std;
#include <iostream>
#include <vector>
#include <cstdio>

#include "Histogram.h"
#include "pCTraw.h"
#include "pCT_Tracking.h"
#include "TkrHits.h"
#include "LineFit.h"

// Gaussian elimination to solve a set of linear equations
vector<double> gauss(vector<vector<double> > A) {
  int n = A.size();

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
  vector<double> x(n);
  for (int i = n - 1; i >= 0; i--) {
    x[i] = A[i][n] / A[i][i];
    for (int k = i - 1; k >= 0; k--) {
      A[k][n] -= A[k][i] * x[i];
    }
  }
  return x;
}

class Stats { // This class is needed only for the monitoring code
  double mean;
  double stdev;

public:
  Stats(std::vector<double> &X) {
    double sum = 0.;
    double sum2 = 0.;
    int N = X.size();
    for (int i = 0; i < N; i++) {
      sum = sum + X[i];
      sum2 = sum2 + X[i] * X[i];
    }
    if (N > 0) {
      mean = sum / float(N);
      double var = sum2 / float(N) - mean * mean;
      stdev = sqrt(var);
    } else {
      mean = 0.;
      stdev = 0.;
    }
  }
  double getMean() { return mean; }
  double getStdev() { return stdev; }
};

class UserAnalysis {

  pCTgeo *uGeometry;
  int numEvents;
  std::string fname;
  std::string OsName;

  float timeFirst;
  double uBeamSpot, uPosE;
  unsigned char mask[5];

  int nEvntOcc;
  int tkrChCnt[num_tkr_fpga][max_chips][64];

  // Create histogram and profile plot pointers

  Histogram *hClSz;
  Histogram *timeStamp;
  Histogram *hTdiff;
  Histogram *hTkSlopeFV;
  Histogram *hTkSlopeFT;
  Histogram *hTkSlopeRV;
  Histogram *hTkSlopeRT;
  Histogram *hTkMissV;
  Histogram *hTkMissT;
  Histogram *hnTkV;
  Histogram *hnTkT;
  Histogram *hEcal;
  Histogram *hWEPL;
  Histogram *hWEPLzero;
  Histogram *hTrig;
  Histogram *deltaT;
  Histogram *deltaTmissV;
  Histogram *deltaTmissT;
  Histogram *deltaTstrip;
  Histogram *deltaTevt;
  Histogram *hVbeamSpot;
  Histogram *hTbeamSpot;
  Histogram *hMxLyr;
  Histogram *hVcalInt, *hTcalInt;
  Histogram *hEnrgTrg;
  Histogram *MxPh0, *MxPh1, *MxPh2, *MxPh3, *MxPh4;
  Histogram *MxSample0, *MxSample1;
  ProfilePlot2Dpeak *hRadioGraph;
  ProfilePlot2D *hRadioGraph2;
  ProfilePlot *hEvsTheta;
  Histogram *hScat;
  Histogram *hEstop[5];
  Histogram *hOTR;
  Histogram *hOTR0;
  Histogram *hStop4[4];
  Histogram *hStop3[3];
  Histogram *hStop2[2];
  Histogram *hStop1;
  Histogram *hTheta;
  Histogram *hErrInt;
  Histogram2D *hTedet;
  Histogram2D *hVedet;

  std::vector<Histogram> hPeds;
  std::vector<Histogram> nClustV;
  std::vector<Histogram> nClustT;
  std::vector<Histogram> PHedet;
  std::vector<Histogram> tHist;
  std::vector<Histogram> vHist;
  std::vector<Histogram> hResidV;
  std::vector<Histogram> hResidT;
  std::vector<Histogram> hSmpMx;
  std::vector<Histogram> hSamplePH;

  std::vector<ProfilePlot> PHvsV;
  std::vector<ProfilePlot> PHvsT;
  std::vector<ProfilePlot> EvsV;
  std::vector<ProfilePlot> EvsT;

  std::vector<Histogram> Vresidual;
  std::vector<Histogram> Tresidual;

  long long minTimeDiff;
  int analysisLevel;
  int runNumber;
  int n_wepl;
  int iEvt;
  int programVersion;
  std::string startTime;
  float stageAngle;
  unsigned long long maxTime, lastTime, lastError;
  int nTgrLyr01, nTgrLyr02, nTgrLyr1, nTgrLyr2;
  vector<int> hits;
  int nWarning;
  float mxChi2;
  float hitWindow;

  int lyrFit[3]; // List of layers to which the straight-line track is to be
                 // fit.
  int lyrA; // The layer being analyzed for efficiency.
  double xHit[4];
  double yHit[4];
  int nVtrialL[4];
  int nVsuccessL[4];
  int nVtrialR[4];
  int nVsuccessR[4];
  int nTtrial[4];
  int nTsuccess[4];

  int nTkEvV; // Number of events with a perfect V track
  int nTkEvT; // Number of events with a perfect T track
  int nMultTk;
  int nInside; // Number of tracks passing through the center of the phantom
  int nTkEvnts; // Number of events with a good track
  int nVecInsideFront; // Number of events with vectors pointing to the phantom
                       // center
  int nVecInsideBack; // Number of events with vectors pointing to the phantom
                      // center
  double mxSlopeFront[2], mxSlopeBack[2];
  string partType;

public:
  UserAnalysis(const char *fnameIn, pCTgeo &Geometry, string partType,
               int analysisLevel, std::string OsName) { // User analysis constructor
    cout << "Now constructing the user analysis instance." << endl;

    // Initialize the member variables

    this->analysisLevel = analysisLevel;
    this->partType = partType;
    this->OsName = OsName;
    fname = fnameIn;
    cout << "User analysis: name of the input data file: " << fnameIn << "\n";

    // Parameters for the efficiency analysis
    mxChi2 = 60.;
    hitWindow = 2.5;

    nEvntOcc = 0;
    for (int FPGA = 0; FPGA < num_tkr_fpga; FPGA++) {
      for (int chip = 0; chip < max_chips; chip++) {
        for (int ch = 0; ch < 64; ch++) {
          tkrChCnt[FPGA][chip][ch] = 0;
        }
      }
    }

    numEvents = 0;
    timeFirst = 1.E26;
    hits.reserve(50);
    nWarning = 0;

    minTimeDiff = 99999999999;
    runNumber = -999;
    programVersion = 0;
    startTime = "none";
    stageAngle = 0.;
    maxTime = 0;
    lastTime = 0;
    lastError = 0;
    nTgrLyr1 = 0;
    nTgrLyr01 = 0;
    nTgrLyr02 = 0;
    nTgrLyr2 = 0;
    uGeometry = &Geometry;

    nTkEvV = 0;
    nTkEvT = 0;
    nMultTk = 0;
    nInside = 0;
    nTkEvnts = 0;
    nVecInsideFront = 0;
    nVecInsideBack = 0;

    memset(nVtrialL, 0, sizeof(nVtrialL));
    memset(nVsuccessL, 0, sizeof(nVsuccessL));
    memset(nVtrialR, 0, sizeof(nVtrialR));
    memset(nVsuccessR, 0, sizeof(nVsuccessR));
    memset(nTtrial, 0, sizeof(nTtrial));
    memset(nTsuccess, 0, sizeof(nTsuccess));

    mask[0] = 0x01; // For unpacking the out of range indicator
    mask[1] = 0x02;
    mask[2] = 0x04;
    mask[3] = 0x08;
    mask[4] = 0x10;

    // Define the histograms and profile plots
    hClSz = new Histogram(64, 0., 1., "Tracker cluster size", "number strips",
                          "clusters");
    timeStamp =
        new Histogram(400, 0., 2., "Event time stamp", "time (s)", "events");
    hTdiff = new Histogram(400, -800., 4., "Time since previous event",
                           "time (10 ns)", "events");
    hErrInt =
        new Histogram(400, 0., 100., "Time since previous time stamp error",
                      "time (10 ns)", "errors");
    hnTkV = new Histogram(20, 0., 1., "Number of tracks in the V view",
                          "Tracks", "events");
    hnTkT = new Histogram(20, 0., 1., "Number of tracks in the T view",
                          "Tracks", "events");
    hTkSlopeFV = new Histogram(
        100, -0.2, .004, "Slope of front track vector in V", "slope", "tracks");
    hTkSlopeFT = new Histogram(
        100, -0.2, .004, "Slope of front track vector in T", "slope", "tracks");
    hTkSlopeRV = new Histogram(
        100, -0.5, .01, "Slope of rear track vector in V", "slope", "tracks");
    hTkSlopeRT = new Histogram(
        100, -0.5, .01, "Slope of rear track vector in T", "slope", "tracks");
    hTkMissV = new Histogram(100, -20., 0.4,
                             "Distance between front and rear V tracks at u=0",
                             "mm", "tracks");
    hTkMissT = new Histogram(100, -20., 0.4,
                             "Distance between front and rear T tracks at u=0",
                             "mm", "tracks");

    hTrig = new Histogram(10, 0., 1., "Trigger bits", "bit", "events");
    deltaT = new Histogram(100, 0., 10., "Time between triggers",
                           "delta-t (ns)", "triggers");
    deltaTmissV =
        new Histogram(100, 0., 10., "Time since last trigger, missed V hit",
                      "delta-t (ns)", "triggers");
    deltaTmissT =
        new Histogram(100, 0., 10., "Time since last trigger, missed T hit",
                      "delta-t (ns)", "triggers");
    deltaTstrip =
        new Histogram(100, 200., 20., "Time between hits on same strip",
                      "delta-t (ns)", "repeated hits");
    deltaTevt = new Histogram(100, 0., 20., "Time between events",
                              "delta-t (ns)", "triggers");
    hEnrgTrg = new Histogram(7, -1., 1., "Energy Detector Trigger Bits Set",
                             "Bit", "events");
    hScat = new Histogram(180, 0., 0.1,
                          "Angle between incoming and outgoing tracks",
                          "degrees", "tracks");
    hOTR0 = new Histogram(10, 0., 1., "Stage overflow indicators", "Stage",
                          "particles");
    hOTR = new Histogram(10, 0., 1.,
                         "Stage overflow indicators for stopping particles",
                         "Stage", "particles");

    hVbeamSpot =
        new Histogram(100, -50., 1., "Front vector V extrapolated to u=-3000",
                      "V0", "vectors");
    hTbeamSpot =
        new Histogram(200, -200., 2., "Front vector T extrapolated to u=-3000",
                      "T0", "vectors");

    hVcalInt = new Histogram(100, -50., 1.,
                             "Track intersection with energy detector in V",
                             "V (mm)", "Tracks");
    hTcalInt = new Histogram(200, -200., 2.,
                             "Track intersection with energy detector in T",
                             "T (mm)", "Tracks");

    MxPh0 = new Histogram(100, 0., 100.,
                          "Stage 0 maximum sample pulse height in 1-trk events",
                          "ph", "events");
    MxPh1 = new Histogram(100, 0., 100.,
                          "Stage 1 maximum sample pulse height in 1-trk events",
                          "ph", "events");
    MxPh2 = new Histogram(100, 0., 100.,
                          "Stage 2 maximum sample pulse height in 1-trk events",
                          "ph", "events");
    MxPh3 = new Histogram(100, 0., 100.,
                          "Stage 3 maximum sample pulse height in 1-trk events",
                          "ph", "events");
    MxPh4 = new Histogram(100, 0., 100.,
                          "Stage 4 maximum sample pulse height in 1-trk events",
                          "ph", "events");
    MxSample0 =
        new Histogram(16, 0, 1., "Stage 0 maximum sample location in pulse",
                      "sample", "events");
    MxSample1 =
        new Histogram(16, 0, 1., "Stage 1 maximum sample location in pulse",
                      "sample", "events");

    PHvsV.push_back(ProfilePlot(100, -50., 1.0,
                                "Energy Detector Stage 0 Pulse Height vs V  ",
                                "V (mm)", "ADC Sum"));
    PHvsV.push_back(ProfilePlot(100, -50., 1.0,
                                "Energy Detector Stage 1 Pulse Height vs V  ",
                                "V (mm)", "ADC Sum"));
    PHvsV.push_back(ProfilePlot(100, -50., 1.0,
                                "Energy Detector Stage 2 Pulse Height vs V  ",
                                "V (mm)", "ADC Sum"));
    PHvsV.push_back(ProfilePlot(100, -50., 1.0,
                                "Energy Detector Stage 3 Pulse Height vs V  ",
                                "V (mm)", "ADC Sum"));
    PHvsV.push_back(ProfilePlot(100, -50., 1.0,
                                "Energy Detector Stage 4 Pulse Height vs V  ",
                                "V (mm)", "ADC Sum"));

    PHvsT.push_back(ProfilePlot(100, -200., 4.0,
                                "Energy Detector Stage 0 Pulse Height vs T  ",
                                "T (mm)", "ADC Sum"));
    PHvsT.push_back(ProfilePlot(100, -200., 4.0,
                                "Energy Detector Stage 1 Pulse Height vs T  ",
                                "T (mm)", "ADC Sum"));
    PHvsT.push_back(ProfilePlot(100, -200., 4.0,
                                "Energy Detector Stage 2 Pulse Height vs T  ",
                                "T (mm)", "ADC Sum"));
    PHvsT.push_back(ProfilePlot(100, -200., 4.0,
                                "Energy Detector Stage 3 Pulse Height vs T  ",
                                "T (mm)", "ADC Sum"));
    PHvsT.push_back(ProfilePlot(100, -200., 4.0,
                                "Energy Detector Stage 4 Pulse Height vs T  ",
                                "T (mm)", "ADC Sum"));

    hSmpMx.push_back(Histogram(100, -1000., 100.,
                               "Energy Detector Stage 0 Maximum Sample PH",
                               "ADC", "protons"));
    hSmpMx.push_back(Histogram(100, -1000., 100.,
                               "Energy Detector Stage 1 Maximum Sample PH",
                               "ADC", "protons"));
    hSmpMx.push_back(Histogram(100, -1000., 100.,
                               "Energy Detector Stage 2 Maximum Sample PH",
                               "ADC", "protons"));
    hSmpMx.push_back(Histogram(100, -1000., 100.,
                               "Energy Detector Stage 3 Maximum Sample PH",
                               "ADC", "protons"));
    hSmpMx.push_back(Histogram(100, -1000., 100.,
                               "Energy Detector Stage 4 Maximum Sample PH",
                               "ADC", "protons"));
    hSmpMx.push_back(Histogram(100, -1000., 100.,
                               "Energy Detector Stage 5 Maximum Sample PH",
                               "ADC", "protons"));

    hSamplePH.push_back(Histogram(
        200, -50., 1., "Energy Detector Stage 0 Sample PH", "ADC", "protons"));
    hSamplePH.push_back(Histogram(
        200, -50., 1., "Energy Detector Stage 1 Sample PH", "ADC", "protons"));
    hSamplePH.push_back(Histogram(
        200, -50., 1., "Energy Detector Stage 2 Sample PH", "ADC", "protons"));
    hSamplePH.push_back(Histogram(
        200, -50., 1., "Energy Detector Stage 3 Sample PH", "ADC", "protons"));
    hSamplePH.push_back(Histogram(
        200, -50., 1., "Energy Detector Stage 4 Sample PH", "ADC", "protons"));
    hSamplePH.push_back(Histogram(
        200, -50., 1., "Energy Detector Stage 5 Sample PH", "ADC", "protons"));

    nClustT.push_back(Histogram(11, -1., 1., "Number of clusters in T layer 1",
                                "Clusters", "Events"));
    nClustT.push_back(Histogram(11, -1., 1., "Number of clusters in T layer 2",
                                "Clusters", "Events"));
    nClustT.push_back(Histogram(11, -1., 1., "Number of clusters in T layer 3",
                                "Clusters", "Events"));
    nClustT.push_back(Histogram(11, -1., 1., "Number of clusters in T layer 4",
                                "Clusters", "Events"));

    nClustV.push_back(Histogram(11, -1., 1., "Number of clusters in V layer 1",
                                "Clusters", "Events"));
    nClustV.push_back(Histogram(11, -1., 1., "Number of clusters in V layer 2",
                                "Clusters", "Events"));
    nClustV.push_back(Histogram(11, -1., 1., "Number of clusters in V layer 3",
                                "Clusters", "Events"));
    nClustV.push_back(Histogram(11, -1., 1., "Number of clusters in V layer 4",
                                "Clusters", "Events"));

    hPeds.push_back(Histogram(100, -300., 6., "Stage 1 FPGA pedestals",
                              "ADC counts", "Triggers"));
    hPeds.push_back(Histogram(100, -300., 6., "Stage 2 FPGA pedestals",
                              "ADC counts", "Triggers"));
    hPeds.push_back(Histogram(100, -300., 6., "Stage 3 FPGA pedestals",
                              "ADC counts", "Triggers"));
    hPeds.push_back(Histogram(100, -300., 6., "Stage 4 FPGA pedestals",
                              "ADC counts", "Triggers"));
    hPeds.push_back(Histogram(100, -300., 6., "Stage 5 FPGA pedestals",
                              "ADC counts", "Triggers"));

    PHedet.push_back(Histogram(100, -1000., 110., "Layer 1 Pulse Height",
                               "ADC sum", "Triggers"));
    PHedet.push_back(Histogram(100, -1000., 110., "Layer 2 Pulse Height",
                               "ADC sum", "Triggers"));
    PHedet.push_back(Histogram(100, -1000., 110., "Layer 3 Pulse Height",
                               "ADC sum", "Triggers"));
    PHedet.push_back(Histogram(100, -1000., 110., "Layer 4 Pulse Height",
                               "ADC sum", "Triggers"));
    PHedet.push_back(Histogram(100, -1000., 110., "Layer 5 Pulse Height",
                               "ADC sum", "Triggers"));

    tHist.push_back(Histogram(440, -200.64, .912,
                              "Tracker layer 1 t Distribution", "t (mm)",
                              "Tracks"));
    tHist.push_back(Histogram(440, -200.64, .912,
                              "Tracker layer 2 t Distribution", "t (mm)",
                              "Tracks"));
    tHist.push_back(Histogram(440, -200.64, .912,
                              "Tracker layer 3 t Distribution", "t (mm)",
                              "Tracks"));
    tHist.push_back(Histogram(440, -200.64, .912,
                              "Tracker layer 4 t Distribution", "t (mm)",
                              "Tracks"));

    vHist.push_back(Histogram(440, -50.16, .228,
                              "Tracker layer 1 v Distribution", "v (mm)",
                              "Tracks"));
    vHist.push_back(Histogram(440, -50.16, .228,
                              "Tracker layer 2 v Distribution", "v (mm)",
                              "Tracks"));
    vHist.push_back(Histogram(440, -50.16, .228,
                              "Tracker layer 3 v Distribution", "v (mm)",
                              "Tracks"));
    vHist.push_back(Histogram(440, -50.16, .228,
                              "Tracker layer 4 v Distribution", "v (mm)",
                              "Tracks"));

    hResidV.push_back(Histogram(200, -1.5, 0.015, "Tracker layer 1 V residuals",
                                "residual (mm)", "Tracks"));
    hResidV.push_back(Histogram(200, -1.5, 0.015, "Tracker layer 2 V residuals",
                                "residual (mm)", "Tracks"));
    hResidV.push_back(Histogram(200, -1.5, 0.015, "Tracker layer 3 V residuals",
                                "residual (mm)", "Tracks"));
    hResidV.push_back(Histogram(200, -1.5, 0.015, "Tracker layer 4 V residuals",
                                "residual (mm)", "Tracks"));

    hResidT.push_back(Histogram(200, -1.5, 0.015, "Tracker layer 1 T residuals",
                                "residual (mm)", "Tracks"));
    hResidT.push_back(Histogram(200, -1.5, 0.015, "Tracker layer 2 T residuals",
                                "residual (mm)", "Tracks"));
    hResidT.push_back(Histogram(200, -1.5, 0.015, "Tracker layer 3 T residuals",
                                "residual (mm)", "Tracks"));
    hResidT.push_back(Histogram(200, -1.5, 0.015, "Tracker layer 4 T residuals",
                                "residual (mm)", "Tracks"));

    Vresidual.push_back(Histogram(80, -4., 0.1, "Layer 1 V Unbiased Residuals",
                                  "mm", "Tracks"));
    Vresidual.push_back(Histogram(80, -4., 0.1, "Layer 2 V Unbiased Residuals",
                                  "mm", "Tracks"));
    Vresidual.push_back(Histogram(80, -4., 0.1, "Layer 3 V Unbiased Residuals",
                                  "mm", "Tracks"));
    Vresidual.push_back(Histogram(80, -4., 0.1, "Layer 4 V Unbiased Residuals",
                                  "mm", "Tracks"));

    Tresidual.push_back(Histogram(80, -4., 0.1, "Layer 1 T Unbiased Residuals",
                                  "mm", "Tracks"));
    Tresidual.push_back(Histogram(80, -4., 0.1, "Layer 2 T Unbiased Residuals",
                                  "mm", "Tracks"));
    Tresidual.push_back(Histogram(80, -4., 0.1, "Layer 3 T Unbiased Residuals",
                                  "mm", "Tracks"));
    Tresidual.push_back(Histogram(80, -4., 0.1, "Layer 4 T Unbiased Residuals",
                                  "mm", "Tracks"));

    hTheta = new Histogram(720, 0., 8., "Projection angle (stage angle)",
                           "Degrees", "events");
    hEvsTheta = new ProfilePlot(720, 0., 8., "Energy vs projection angle",
                                "Degrees", "MeV");

    if (analysisLevel >= 2) {
      hVedet = new Histogram2D(50, -50., 2.0, 50, -100., 8.0,
                               "Projection of back vector to edet in V",
                               "V (mm)", "WEPL (mm)", "Tracks");
      hTedet = new Histogram2D(50, -150., 6.0, 50, -100., 8.0,
                               "Projection of back vector to edet in T",
                               "V (mm)", "WEPL (mm)", "Tracks");
      hMxLyr =
          new Histogram(7, -1., 1., "Deepest stage hit", "stage", "events");

      EvsV.push_back(ProfilePlot(100, -50., 1.0,
                                 "Energy Detector Stage 0 Energy vs V  ",
                                 "V (mm)", "MeV"));
      EvsV.push_back(ProfilePlot(100, -50., 1.0,
                                 "Energy Detector Stage 1 Energy vs V  ",
                                 "V (mm)", "MeV"));
      EvsV.push_back(ProfilePlot(100, -50., 1.0,
                                 "Energy Detector Stage 2 Energy vs V  ",
                                 "V (mm)", "MeV"));
      EvsV.push_back(ProfilePlot(100, -50., 1.0,
                                 "Energy Detector Stage 3 Energy vs V  ",
                                 "V (mm)", "MeV"));
      EvsV.push_back(ProfilePlot(100, -50., 1.0,
                                 "Energy Detector Stage 4 Energyt vs V  ",
                                 "V (mm)", "MeV"));

      EvsT.push_back(ProfilePlot(100, -200., 4.0,
                                 "Energy Detector Stage 0 Energy vs T  ",
                                 "T (mm)", "MeV"));
      EvsT.push_back(ProfilePlot(100, -200., 4.0,
                                 "Energy Detector Stage 1 Energy vs T  ",
                                 "T (mm)", "MeV"));
      EvsT.push_back(ProfilePlot(100, -200., 4.0,
                                 "Energy Detector Stage 2 Energy vs T  ",
                                 "T (mm)", "MeV"));
      EvsT.push_back(ProfilePlot(100, -200., 4.0,
                                 "Energy Detector Stage 3 Energy vs T  ",
                                 "T (mm)", "MeV"));
      EvsT.push_back(ProfilePlot(100, -200., 4.0,
                                 "Energy Detector Stage 4 Energy vs T  ",
                                 "T (mm)", "MeV"));

      if (partType == "H") {
        hEcal = new Histogram(150, 0., 2., "Calibrated sum of stage energies",
                              "mm", "protons");
      } else {
        hEcal = new Histogram(150, 0., 8., "Calibrated sum of stage energies",
                              "mm", "He ions");
      }
      hWEPL = new Histogram(200, -100., 2., "Calibrated WEPL", "mm", "protons");
      hWEPLzero = new Histogram(
          200, -100., 2., "Calibrated WEPL for |t|>110 mm", "mm", "protons");
      hRadioGraph =
          new ProfilePlot2Dpeak(220, -90., 1.0, 80, -40., 1.0, 100, -50., 3.,
                                "Radiograph", "T (mm)", "V (mm)", "WEPL (mm)");
      hRadioGraph2 =
          new ProfilePlot2D(440, -90., 0.5, 160, -40., 0.5, "Radiograph",
                            "T (mm)", "V (mm)", "WEPL (mm)");
      if (partType == "H") {
        hEstop[0] = new Histogram(400, 15., 0.35,
                                  "Stopping energy distribution for stage 0",
                                  "Energy (MeV)", "proton events");
        hEstop[1] = new Histogram(400, 15., 0.35,
                                  "Stopping energy distribution for stage 1",
                                  "Energy (MeV)", "proton events");
        hEstop[2] = new Histogram(400, 15., 0.35,
                                  "Stopping energy distribution for stage 2",
                                  "Energy (MeV)", "proton events");
        hEstop[3] = new Histogram(400, 15., 0.35,
                                  "Stopping energy distribution for stage 3",
                                  "Energy (MeV)", "proton events");
        hEstop[4] = new Histogram(400, 15., 0.35,
                                  "Stopping energy distribution for stage 4",
                                  "Energy (MeV)", "proton events");
        hStop4[0] = new Histogram(400, 60., 0.35,
                                  "Stage 0 energy for stopping in stage 4",
                                  "Energy (MeV)", "proton Events");
        hStop4[1] = new Histogram(400, 60., 0.35,
                                  "Stage 1 energy for stopping in stage 4",
                                  "Energy (MeV)", "proton Events");
        hStop4[2] = new Histogram(400, 60., 0.35,
                                  "Stage 2 energy for stopping in stage 4",
                                  "Energy (MeV)", "proton Events");
        hStop4[3] = new Histogram(400, 60., 0.35,
                                  "Stage 3 energy for stopping in stage 4",
                                  "Energy (MeV)", "proton Events");
        hStop3[0] = new Histogram(400, 60., 0.35,
                                  "Stage 0 energy for stopping in stage 3",
                                  "Energy (MeV)", "proton Events");
        hStop3[1] = new Histogram(400, 60., 0.35,
                                  "Stage 1 energy for stopping in stage 3",
                                  "Energy (MeV)", "proton Events");
        hStop3[2] = new Histogram(400, 60., 0.35,
                                  "Stage 2 energy for stopping in stage 3",
                                  "Energy (MeV)", "proton Events");
        hStop2[0] = new Histogram(400, 60., 0.35,
                                  "Stage 0 energy for stopping in stage 2",
                                  "Energy (MeV)", "proton Events");
        hStop2[1] = new Histogram(400, 60., 0.35,
                                  "Stage 1 energy for stopping in stage 2",
                                  "Energy (MeV)", "proton Events");
        hStop1 = new Histogram(400, 60., 0.35,
                               "Stage 0 energy for stopping in stage 1",
                               "Energy (MeV)", "proton Events");
      } else {
        hEstop[0] = new Histogram(400, 60., 1.4,
                                  "Stopping energy distribution for stage 0",
                                  "Energy (MeV)", "He events");
        hEstop[1] = new Histogram(400, 60., 1.4,
                                  "Stopping energy distribution for stage 1",
                                  "Energy (MeV)", "He events");
        hEstop[2] = new Histogram(400, 60., 1.4,
                                  "Stopping energy distribution for stage 2",
                                  "Energy (MeV)", "He events");
        hEstop[3] = new Histogram(400, 60., 1.4,
                                  "Stopping energy distribution for stage 3",
                                  "Energy (MeV)", "He events");
        hEstop[4] = new Histogram(400, 60., 1.4,
                                  "Stopping energy distribution for stage 4",
                                  "Energy (MeV)", "He events");
        hStop4[0] = new Histogram(400, 60., 1.4,
                                  "Stage 0 energy for stopping in stage 4",
                                  "Energy (MeV)", "He Events");
        hStop4[1] = new Histogram(400, 60., 1.4,
                                  "Stage 1 energy for stopping in stage 4",
                                  "Energy (MeV)", "He Events");
        hStop4[2] = new Histogram(400, 60., 1.4,
                                  "Stage 2 energy for stopping in stage 4",
                                  "Energy (MeV)", "He Events");
        hStop4[3] = new Histogram(400, 60., 1.4,
                                  "Stage 3 energy for stopping in stage 4",
                                  "Energy (MeV)", "He Events");
        hStop3[0] = new Histogram(400, 60., 1.4,
                                  "Stage 0 energy for stopping in stage 3",
                                  "Energy (MeV)", "He Events");
        hStop3[1] = new Histogram(400, 60., 1.4,
                                  "Stage 1 energy for stopping in stage 3",
                                  "Energy (MeV)", "He Events");
        hStop3[2] = new Histogram(400, 60., 1.4,
                                  "Stage 2 energy for stopping in stage 3",
                                  "Energy (MeV)", "He Events");
        hStop2[0] = new Histogram(400, 60., 1.4,
                                  "Stage 0 energy for stopping in stage 2",
                                  "Energy (MeV)", "He Events");
        hStop2[1] = new Histogram(400, 60., 1.4,
                                  "Stage 1 energy for stopping in stage 2",
                                  "Energy (MeV)", "He Events");
        hStop1 = new Histogram(400, 60., 1.4,
                               "Stage 0 energy for stopping in stage 1",
                               "Energy (MeV)", "He Events");
      }
    }
  }; // End of user analysis constructor

  void initialize(const pCTraw &pCTEvent) { // User analysis initialization,
                                            // called after reading the run
                                            // header
    cout << "UserAnalysis run initialization. . .\n";
    n_wepl = 0;
    if (pCTEvent.run_number >= 0) {
      runNumber = pCTEvent.run_number;
      programVersion = pCTEvent.program_version;
      startTime = pCTEvent.start_time;
      stageAngle = pCTEvent.stage_angle;
    }

    uBeamSpot = -2000.; // Assumed location of the beam origin (just a guess)
    uPosE = 215.;       // Rough u position of the front of the energy detector
    mxSlopeFront[0] = 0.03;
    mxSlopeFront[1] = 0.09;
    mxSlopeBack[0] = 0.1;
    mxSlopeBack[1] = 0.15;

  }; // End of user analysis initialization

  void weplEvent(float theta, float WEPL, float energy[5], float Thit[4],
                 double Uhit[4], float Vhit[4],
                 unsigned char OTR) { // User event analysis called for each
                                      // WEPL reconstructed
    n_wepl++;
    float Etot = 0.;
    for (int stage = 0; stage < 5; stage++) {
      Etot += energy[stage];
    }

    int lstLyr = -1;
    for (int k = 0; k < 5; k++) {
      if (energy[k] < 1.0)
        break;
      lstLyr = k;
    }
    hMxLyr->entry((double)lstLyr);

    // Look at the energy detector T and V dependence
    double Tback[2], Vback[2];
    double Tfront[2], Vfront[2];
    for (int lyr = 2; lyr < 4; lyr++) {
      Tback[lyr - 2] = Thit[lyr];
      Vback[lyr - 2] = Vhit[lyr];
    }
    for (int lyr = 0; lyr < 2; lyr++) {
      Tfront[lyr] = Thit[lyr];
      Vfront[lyr] = Vhit[lyr];
    }
    float Vedet[5], Tedet[5];
    float Uphantom =
        -2.0 * uGeometry->getBrickThickness(); // rough U location of the
                                               // calibration phantom
    float Vphantom = uGeometry->extrap2D(Uhit, Vfront, Uphantom);
    float Tphantom = uGeometry->extrap2D(Uhit, Tfront, Uphantom);
    for (int lyr = 0; lyr < 5; ++lyr) {
      Vedet[lyr] =
          uGeometry->extrap2D(&Uhit[2], Vback, uGeometry->energyDetectorU(lyr));
      Tedet[lyr] =
          uGeometry->extrap2D(&Uhit[2], Tback, uGeometry->energyDetectorU(lyr));
      EvsV[lyr].entry(Vphantom, energy[lyr]);
      EvsT[lyr].entry(Tphantom, energy[lyr]);
    }
    hVedet->entry((double)Vedet[0], (double)WEPL);
    hTedet->entry((double)Tedet[0], (double)WEPL);
    hEcal->entry(Etot);
    hWEPL->entry(WEPL);
    if (Tedet[0] > 110. || Tedet[0] < -110.) {
      hEvsTheta->entry(theta, Etot);
      hWEPLzero->entry(WEPL);
    }

    unsigned char tst[5];
    for (int stage = 0; stage < 5; ++stage)
      tst[stage] = mask[stage] & OTR;
    if (energy[4] > 1.0) {
      hEstop[4]->entry(energy[4]);
      if (tst[4] != 0)
        hOTR->entry(4);
      for (int stage = 0; stage < 4; ++stage)
        hStop4[stage]->entry(energy[stage]);
    } else if (energy[3] > 1.0) {
      hEstop[3]->entry(energy[3]);
      if (tst[3] != 0)
        hOTR->entry(3);
      for (int stage = 0; stage < 3; ++stage)
        hStop3[stage]->entry(energy[stage]);
    } else if (energy[2] > 1.0) {
      hEstop[2]->entry(energy[2]);
      if (tst[2] != 0)
        hOTR->entry(2);
      for (int stage = 0; stage < 2; ++stage)
        hStop2[stage]->entry(energy[stage]);
    } else if (energy[1] > 1.0) {
      hEstop[1]->entry(energy[1]);
      if (tst[1] != 0)
        hOTR->entry(1);
      hStop1->entry(energy[0]);
    } else if (energy[0] > 1.0) {
      hEstop[0]->entry(energy[0]);
      if (tst[0] != 0)
        hOTR->entry(0);
    }

    double vin[3], vout[3];
    vin[0] = Vhit[1] - Vhit[0];
    vin[1] = Thit[1] - Thit[0];
    vin[2] = Uhit[1] - Uhit[0];
    vout[0] = Vhit[3] - Vhit[2];
    vout[1] = Thit[3] - Thit[2];
    vout[2] = Uhit[3] - Uhit[2];
    double vinNrm = sqrt(vin[0] * vin[0] + vin[1] * vin[1] + vin[2] * vin[2]);
    double voutNrm =
        sqrt(vout[0] * vout[0] + vout[1] * vout[1] + vout[2] * vout[2]);
    for (int i = 0; i < 3; ++i) {
      vin[i] = vin[i] / vinNrm;
      vout[i] = vout[i] / voutNrm;
    }
    double ct = vin[0] * vout[0] + vin[1] * vout[1] + vin[2] * vout[2];
    double thetaScat = acos(ct) * 57.295779513;
    hScat->entry(thetaScat);

    // Interpolate with cubic polynomials to find a position in the middle of
    // the phantom location, to make a simple radiograph
    std::vector<std::vector<double> > A(4, vector<double>(5, 0.));
    A[0][0] = 1.;
    A[1][0] = 1.;
    A[2][0] = 1.;
    A[3][0] = 1.;
    A[0][1] = Uhit[0];
    A[1][1] = Uhit[1];
    A[2][1] = Uhit[2];
    A[3][1] = Uhit[3];
    A[0][2] = 0.01 * pow(Uhit[0], 2);
    A[1][2] = 0.01 * pow(Uhit[1], 2);
    A[2][2] = 0.01 * pow(Uhit[2], 2);
    A[3][2] = 0.01 * pow(Uhit[3], 2);
    A[0][3] = 0.0001 * pow(Uhit[0], 3);
    A[1][3] = 0.0001 * pow(Uhit[1], 3);
    A[2][3] = 0.0001 * pow(Uhit[2], 3);
    A[3][3] = 0.0001 * pow(Uhit[3], 3);
    A[0][4] = Thit[0];
    A[1][4] = Thit[1];
    A[2][4] = Thit[2];
    A[3][4] = Thit[3];
    std::vector<double> bt = gauss(A);
    A[0][4] = Vhit[0];
    A[1][4] = Vhit[1];
    A[2][4] = Vhit[2];
    A[3][4] = Vhit[3];
    std::vector<double> bv = gauss(A);
    double WEPL2 = (double)WEPL;
    if (WEPL2 > -20. && WEPL2 < 300. && thetaScat < 2.) {
      hRadioGraph->entry(bt[0], bv[0], WEPL2);
      hRadioGraph2->entry(bt[0], bv[0], WEPL2);
    }
    if (n_wepl % 100000000 == 0) {
      cout << "UserAnalysis weplEvent " << n_wepl << ":  WEPL=" << WEPL
           << "  E= ";
      for (int i = 0; i < 5; i++)
        cout << energy[i] << " ";
      cout << endl << "  T= ";
      for (int i = 0; i < 4; i++)
        cout << Thit[i] << " ";
      cout << endl << "  U= ";
      for (int i = 0; i < 4; i++)
        cout << Uhit[i] << " ";
      cout << endl << "  V= ";
      for (int i = 0; i < 4; i++)
        cout << Vhit[i] << " ";
      cout << endl << "UserAnalysis::weplEvent: U,T polynomial coefficients:";
      for (int i = 0; i < 4; i++)
        cout << bt[i] << " ";
      cout << endl << "UserAnalysis::weplEvent: U,V polynomial coefficients:";
      for (int i = 0; i < 4; i++)
        cout << bv[i] << " ";
      cout << endl;
    }
  }

  bool rawEvent(const pCTraw &pCTEvent, const TkrHits &pCThits,
                pCT_Tracking &pCTtracks, pCTgeo &geometry,
                float theta) { // User event analysis called for each raw event
                               // read in

    bool cutEvent =
        false; // Set this true to eliminate this event from the processing

    numEvents++;
    if (numEvents % 100000 == 0)
      cout << "UserAnalysis rawEvent " << numEvents
           << " Time tag= " << pCTEvent.time_tag << "\n";
    if (pCTEvent.time_tag > maxTime)
      maxTime = pCTEvent.time_tag;

    float timeSeconds = float(pCTEvent.time_tag) * 1.0E-8;
    if (timeSeconds < timeFirst) {
      cout << "User rawEvent: Setting the earliest time stamp to "
           << timeSeconds << "\n";
      timeFirst = timeSeconds;
    }

    long long timeDiff;
    if (pCTEvent.time_tag < lastTime) {
      timeDiff = lastTime - pCTEvent.time_tag;
      timeDiff = -timeDiff;
      float errorInt = (float)(lastTime - lastError);
      hErrInt->entry(errorInt);
      lastError = lastTime;
      if (nWarning < 11) {
        cout << "Warning: TIME TAG DECREASE encountered!  Old=" << lastTime
             << " New=" << pCTEvent.time_tag;
        cout << "  Difference=" << timeDiff
             << "  Time since last trigger = " << pCTEvent.delta_t << endl;
        cout << "  Time interval from the last tag error = " << errorInt
             << "     Time tag=" << pCTEvent.time_tag << endl;
        nWarning++;
      }
    } else {
      timeDiff = pCTEvent.time_tag - lastTime;
    }
    if (timeDiff < minTimeDiff)
      minTimeDiff = timeDiff;
    hTdiff->entry((float)timeDiff);
    timeStamp->entry(timeSeconds - timeFirst);
    hTheta->entry(theta);

    lastTime = pCTEvent.time_tag;

    for (int bit = 0; bit < 6; bit++) {
      if (pCTEvent.trigger_bits[bit])
        hTrig->entry(bit);
    }

    int phSum[5];
    phSum[0] = pCTEvent.enrg_fpga[0].pulse_sum[0];
    phSum[1] = pCTEvent.enrg_fpga[0].pulse_sum[1];
    phSum[2] = pCTEvent.enrg_fpga[0].pulse_sum[2];
    phSum[3] = pCTEvent.enrg_fpga[1].pulse_sum[0];
    phSum[4] = pCTEvent.enrg_fpga[1].pulse_sum[1];
    for (int i = 0; i < 5; i++) {
      PHedet[i].entry(phSum[i]);
    }

    // Energy Detector Pedestals
    if (pCTEvent.enrg_fpga[0].peds_out) {
      hPeds[0].entry(pCTEvent.enrg_fpga[0].pedestal[0]);
      hPeds[1].entry(pCTEvent.enrg_fpga[0].pedestal[1]);
      hPeds[2].entry(pCTEvent.enrg_fpga[0].pedestal[2]);
      hPeds[3].entry(pCTEvent.enrg_fpga[1].pedestal[0]);
      hPeds[4].entry(pCTEvent.enrg_fpga[1].pedestal[1]);
    }

    // Energy Detector Overflows
    if (pCTEvent.enrg_fpga[0].OTR[0])
      hOTR0->entry(0);
    if (pCTEvent.enrg_fpga[0].OTR[1])
      hOTR0->entry(1);
    if (pCTEvent.enrg_fpga[0].OTR[2])
      hOTR0->entry(2);
    if (pCTEvent.enrg_fpga[1].OTR[0])
      hOTR0->entry(3);
    if (pCTEvent.enrg_fpga[1].OTR[1])
      hOTR0->entry(4);

    // Energy Detector Samples
    for (int fpga = 0; fpga < 2; fpga++) {
      int nSamp = pCTEvent.enrg_fpga[fpga].num_samples;
      if (nSamp > 0) {
        for (int chan = 0; chan < 3; chan++) {
          int smpMax = -999;
          for (int smp = 0; smp < 16; smp++) {
            int smpl = pCTEvent.enrg_fpga[fpga].sample[chan][smp];
            hSamplePH[3 * fpga + chan].entry(smpl);
            if (smpl > smpMax)
              smpMax = smpl;
          }
          hSmpMx[3 * fpga + chan].entry(smpMax);
        }
      }
    }

    deltaT->entry(10.0 * (double)pCTEvent.delta_t);

    // Tracker channel occupancy counts
    nEvntOcc++;
    for (int FPGA = 0; FPGA < num_tkr_fpga; FPGA++) {
      for (int chip = 0; chip < pCTEvent.tkr_fpga[FPGA].num_chips; chip++) {
        int address = pCTEvent.tkr_fpga[FPGA].chip[chip].address;
        int nhits = pCTEvent.tkr_fpga[FPGA].chip[chip].num_clusts;
        for (int hit = 0; hit < nhits; hit++) {
          int frst = pCTEvent.tkr_fpga[FPGA].chip[chip].cluster[hit].first;
          int len = pCTEvent.tkr_fpga[FPGA].chip[chip].cluster[hit].length + 1;
          hClSz->entry(len);
          for (int i = 0; i < len; i++) {
            tkrChCnt[FPGA][address][frst + i]++;
          }
        }
      }
    }

    bool oneTkr = true;
    for (int lyr = 0; lyr < 4; lyr++) {
      nClustV[lyr].entry(pCThits.Lyr[lyr].N[0]);
      nClustT[lyr].entry(pCThits.Lyr[lyr].N[1]);
      if (pCThits.Lyr[lyr].N[0] != 1)
        oneTkr = false;
      if (pCThits.Lyr[lyr].N[1] != 1)
        oneTkr = false;
    }

    // Count number of tracks passing through the phantom center
    if (pCTtracks.nTracks > 0) {
      double vCenter = (pCTtracks.frontPredV(pCTtracks.itkV, 0.) +
                        pCTtracks.backPredV(pCTtracks.itkV, 0.)) /
                       2.;
      double tCenter = (pCTtracks.frontPredT(pCTtracks.itkT, 0.) +
                        pCTtracks.backPredT(pCTtracks.itkT, 0.)) /
                       2.;
      if (abs(vCenter) < 5.0 && abs(tCenter) < 5.0)
        nInside++;
      nTkEvnts++;
    }

    // Count front/back vectors pointing to the phantom center
    bool box[2] = { false, false };
    for (int Idx = 0; Idx < 2; Idx++) {
      for (int i = 0; i < pCThits.Lyr[0].N[Idx]; i++) {
        for (int j = 0; j < pCThits.Lyr[1].N[Idx]; j++) {
          double slope =
              (pCThits.Lyr[1].Y[Idx].at(j) - pCThits.Lyr[0].Y[Idx].at(i)) /
              (pCThits.Lyr[1].U[Idx].at(j) - pCThits.Lyr[0].U[Idx].at(i));
          if (abs(slope) < mxSlopeFront[Idx]) {
            double Y = pCThits.Lyr[0].Y[Idx].at(i);
            double U = pCThits.Lyr[0].U[Idx].at(i);
            double intercept = Y - slope * U;
            if (abs(intercept) < 5.0)
              box[Idx] = true;
          }
        }
      }
    }
    if (box[0] && box[1])
      nVecInsideFront++;
    box[0] = false;
    box[1] = false;
    for (int Idx = 0; Idx < 2; Idx++) {
      for (int i = 0; i < pCThits.Lyr[2].N[Idx]; i++) {
        for (int j = 0; j < pCThits.Lyr[3].N[Idx]; j++) {
          double slope =
              (pCThits.Lyr[3].Y[Idx].at(j) - pCThits.Lyr[2].Y[Idx].at(i)) /
              (pCThits.Lyr[3].U[Idx].at(j) - pCThits.Lyr[2].U[Idx].at(i));
          if (abs(slope) < mxSlopeBack[Idx]) {
            double Y = pCThits.Lyr[2].Y[Idx].at(i);
            double U = pCThits.Lyr[2].U[Idx].at(i);
            double intercept = Y - slope * U;
            if (abs(intercept) < 5.0)
              box[Idx] = true;
          }
        }
      }
    }
    if (box[0] && box[1])
      nVecInsideBack++;

    // Analyze the cut variables used to define tracks

    int nVtracks = pCTtracks.VTracks.size();
    hnTkV->entry(nVtracks);
    int nTtracks = pCTtracks.TTracks.size();
    hnTkT->entry(nTtracks);
    if (pCTtracks.nTracks == 1) {
      double missV = pCTtracks.VTracks[pCTtracks.itkV].Miss;
      hTkMissV->entry(missV);
      double missT = pCTtracks.TTracks[pCTtracks.itkT].Miss;
      hTkMissT->entry(missT);
      double slopeRV = (pCTtracks.VTracks[pCTtracks.itkV].X[3] -
                        pCTtracks.VTracks[pCTtracks.itkV].X[2]) /
                       (pCTtracks.VTracks[pCTtracks.itkV].U[3] -
                        pCTtracks.VTracks[pCTtracks.itkV].U[2]);
      hTkSlopeRV->entry(slopeRV);
      double slopeRT = (pCTtracks.TTracks[pCTtracks.itkT].X[3] -
                        pCTtracks.TTracks[pCTtracks.itkT].X[2]) /
                       (pCTtracks.TTracks[pCTtracks.itkT].U[3] -
                        pCTtracks.TTracks[pCTtracks.itkT].U[2]);
      hTkSlopeRT->entry(slopeRT);
      double slopeFV = (pCTtracks.VTracks[pCTtracks.itkV].X[1] -
                        pCTtracks.VTracks[0].X[pCTtracks.itkV]) /
                       (pCTtracks.VTracks[pCTtracks.itkV].U[1] -
                        pCTtracks.VTracks[pCTtracks.itkV].U[0]);
      hTkSlopeFV->entry(slopeFV);
      double slopeFT = (pCTtracks.TTracks[pCTtracks.itkT].X[1] -
                        pCTtracks.TTracks[pCTtracks.itkT].X[0]) /
                       (pCTtracks.TTracks[pCTtracks.itkT].U[1] -
                        pCTtracks.TTracks[pCTtracks.itkT].U[0]);
      hTkSlopeFT->entry(slopeFT);

      hVbeamSpot->entry(pCTtracks.frontPredV(pCTtracks.itkV, uBeamSpot));
      hTbeamSpot->entry(pCTtracks.frontPredT(pCTtracks.itkT, uBeamSpot));

      // Measure and make histograms of tracker residuals (most useful in empty
      // runs, where the tracks should be more-or-less straight)
      double vHits[4], uvHits[4], tHits[4], utHits[4];
      for (int lyr = 0; lyr < 4; lyr++) {
        vHits[lyr] = pCTtracks.VTracks[pCTtracks.itkV].X[lyr];
        uvHits[lyr] = pCTtracks.VTracks[pCTtracks.itkV].U[lyr];
        tHits[lyr] = pCTtracks.TTracks[pCTtracks.itkT].X[lyr];
        utHits[lyr] = pCTtracks.TTracks[pCTtracks.itkT].U[lyr];
      }
      LineFit lineFitV(4, uvHits, vHits); // least-squares fit of the V-view
                                          // hits to a straight line
      LineFit lineFitT(4, utHits, tHits);
      if (lineFitV.getChi2() < 4.0 && lineFitT.getChi2() < 4.0) {
        for (int lyr = 0; lyr < 4; lyr++) {
          hResidV[lyr].entry(vHits[lyr] - lineFitV.eval(uvHits[lyr]));
          hResidT[lyr].entry(tHits[lyr] - lineFitT.eval(utHits[lyr]));
        }
      }

      // Trigger efficiency analysis
      double vCal = -999.;
      vCal = pCTtracks.backPredV(pCTtracks.itkV, geometry.energyDetectorU(0));
      hVcalInt->entry(vCal);
      double tCal = -999.;
      tCal = pCTtracks.backPredT(pCTtracks.itkT, geometry.energyDetectorU(0));
      hTcalInt->entry(tCal);
      if (abs(vCal) < 35. && abs(tCal) < 125.) {
        for (int i = 0; i < 6; i++)
          if (pCTEvent.trigger_bits[i])
            hEnrgTrg->entry(i);
        hEnrgTrg->entry(-1);
        if (pCTEvent.trigger_bits[1]) {
          nTgrLyr1++;
          if (pCTEvent.trigger_bits[0])
            nTgrLyr01++;
        }
        if (pCTEvent.trigger_bits[2]) {
          nTgrLyr2++;
          if (pCTEvent.trigger_bits[0])
            nTgrLyr02++;
        }
      }

      // Look at the energy detector T and V dependence
      if (oneTkr) { // Extra hard cut to select single-track events
        for (int lyr = 0; lyr < 5; ++lyr) {
          PHvsV[lyr].entry(vCal, phSum[lyr]);
          PHvsT[lyr].entry(tCal, phSum[lyr]);
        }
      }

      // Look at the raw digitizations in the case that there is a single track
      int valmax0 = -1000;
      int valmax1 = -1000;
      int valmax2 = -1000;
      int valmax3 = -1000;
      int valmax4 = -1000;
      int smpMax0 = -1;
      int smpMax1 = -1;
      if (oneTkr) {
        if (phSum[1] < 6000000) {
          //                if (phSum[1]<6000 && phSum[2]< 4500 && phSum[3] <
          // 8000 && phSum[4]<7000) {   // These cuts are relevant for protons
          // only
          for (int ismp = 0; ismp < 16; ismp++) {
            int val0 = pCTEvent.enrg_fpga[0].sample[0][ismp];
            int val1 = pCTEvent.enrg_fpga[0].sample[1][ismp];
            int val2 = pCTEvent.enrg_fpga[0].sample[2][ismp];
            int val3 = pCTEvent.enrg_fpga[1].sample[0][ismp];
            int val4 = pCTEvent.enrg_fpga[1].sample[1][ismp];
            if (val0 > valmax0) {
              valmax0 = val0;
              smpMax0 = ismp;
            }
            if (val1 > valmax1) {
              valmax1 = val1;
              smpMax1 = ismp;
            }
            if (val2 > valmax2) {
              valmax2 = val2;
            }
            if (val3 > valmax3) {
              valmax3 = val3;
            }
            if (val4 > valmax4) {
              valmax4 = val4;
            }
          }
        }
      }

      MxPh0->entry(valmax0);
      MxPh1->entry(valmax1);
      MxPh2->entry(valmax2);
      MxPh3->entry(valmax3);
      MxPh4->entry(valmax4);
      MxSample0->entry(smpMax0);
      MxSample1->entry(smpMax1);
    }

    // Count numbers of events with perfect V and T tracks
    int Nperfect = 0;
    for (int i = 0; i < pCTtracks.VTracks.size(); i++) {
      if (pCTtracks.VTracks.at(i).Q == 3)
        Nperfect++;
    }
    if (Nperfect > 0)
      nTkEvV++;
    Nperfect = 0;
    for (int i = 0; i < pCTtracks.TTracks.size(); i++) {
      if (pCTtracks.TTracks.at(i).Q == 3)
        Nperfect++;
    }
    if (Nperfect > 0)
      nTkEvT++;
    if (pCTtracks.VTracks.size() > 1 || pCTtracks.TTracks.size() > 1)
      nMultTk++;

    // Accumulate tracker hit distributions
    for (int lyr = 0; lyr < 4; lyr++) {
      for (int i = 0; i < pCThits.Lyr[lyr].N[0]; i++) {
        vHist[lyr].entry(pCThits.Lyr[lyr].Y[0].at(i));
      }
      for (int i = 0; i < pCThits.Lyr[lyr].N[1]; i++) {
        tHist[lyr].entry(pCThits.Lyr[lyr].Y[1].at(i));
      }
    }

    int nVtkrs = 0;
    int nTtkrs = 0;
    for (int i = 0; i < pCTtracks.VTracks.size(); i++) {
      if (pCTtracks.VTracks[i].Good) {
        nVtkrs++;
      }
    }
    for (int i = 0; i < pCTtracks.TTracks.size(); i++) {
      if (pCTtracks.TTracks[i].Good) {
        nTtkrs++;
      }
    }
    // Tracker hit-efficiency analysis
    for (lyrA = 0; lyrA < 4;
         lyrA++) { // Study all 4 layers successively. lyrA is the one being
                   // measured for efficiency
      int cnt = 0;
      for (int lyr = 0; lyr < 4; lyr++) {
        if (lyr != lyrA) {
          lyrFit[cnt] = lyr;
          cnt++;
        }
      }
      // Analyze the hit efficiency in the V layer specified by lyrA
      if (nTtkrs == 1 && pCThits.Lyr[0].N[1] == 1 && pCThits.Lyr[1].N[1] == 1 &&
          pCThits.Lyr[2].N[1] == 1 &&
          pCThits.Lyr[3].N[1] == 1) { // 1 track in T with 4 perfect hits
        if (pCThits.Lyr[lyrFit[0]].N[0] == 1 &&
            pCThits.Lyr[lyrFit[1]].N[0] == 1 &&
            pCThits.Lyr[lyrFit[2]].N[0] ==
                1) { // All other V layers hit, but just once
          if (phSum[0] > 2500 && phSum[0] < 8500 && phSum[4] > -3200 &&
              phSum[4] < 4800) { // Filter out multi-proton events
            for (int lyr = 0; lyr < 4; lyr++) {
              xHit[lyr] = pCThits.Lyr[lyr].U[1].at(0);
              yHit[lyr] = pCThits.Lyr[lyr].Y[1].at(0);
            }
            LineFit LineFitT(4, xHit, yHit);
            if (LineFitT.getChi2() < mxChi2) { // Good track in T that is within
                                               // bounds at the V layer of
                                               // interest
              double Textrap = LineFitT.eval(geometry.uV(lyrA));
              if (Textrap > -155. && Textrap < 155.) {
                for (int lyr = 0; lyr < 3; lyr++) {
                  xHit[lyr] = pCThits.Lyr[lyrFit[lyr]].U[0].at(0);
                  yHit[lyr] = pCThits.Lyr[lyrFit[lyr]].Y[0].at(0);
                }
                LineFit LineFitV(3, xHit, yHit); // Fit a track in V to the 3
                                                 // layers not being studied
                double vPos = LineFitV.eval(geometry.uV(lyrA));
                if (LineFitV.getChi2() < mxChi2 && vPos > -40. &&
                    vPos < 40.) { // Must be a good fit in V
                  if (Textrap > 0.)
                    nVtrialL[lyrA]++;
                  else
                    nVtrialR[lyrA]++; // Study the two V ladders separately
                  bool success = false;
                  for (int hit = 0; hit < pCThits.Lyr[lyrA].N[0];
                       hit++) { // Note: I don't know which ladder the hit was
                                // in!!!
                    double resid = pCThits.Lyr[lyrA].Y[0].at(hit) - vPos;
                    Vresidual[lyrA].entry(resid);
                    if (abs(resid) < hitWindow) { // Is there a V hit in the
                                                  // layer under study close to
                                                  // the track?
                      if (Textrap > 0.)
                        nVsuccessL[lyrA]++;
                      else
                        nVsuccessR[lyrA]++;
                      success = true;
                      break;
                    }
                  }
                  if (!success) { // We missed a hit on this layer. What does
                                  // the time since previous trigger look like?
                    deltaTmissV->entry(10.0 * (double)pCTEvent.delta_t);
                  }
                }
              }
            }
          }
        }
      }

      // Analyze the hit efficiency in the T layer specified by lyrA
      if (nVtkrs == 1 && pCThits.Lyr[0].N[0] == 1 && pCThits.Lyr[1].N[0] == 1 &&
          pCThits.Lyr[2].N[0] == 1 && pCThits.Lyr[3].N[0] == 1) {
        if (pCThits.Lyr[lyrFit[0]].N[1] == 1 &&
            pCThits.Lyr[lyrFit[1]].N[1] == 1 &&
            pCThits.Lyr[lyrFit[2]].N[1] == 1) {
          if (phSum[0] > 2500 && phSum[0] < 8500 && phSum[4] > -3200 &&
              phSum[4] < 4800) { // Select single proton events
            for (int lyr = 0; lyr < 4; lyr++) {
              xHit[lyr] = pCThits.Lyr[lyr].U[0].at(0);
              yHit[lyr] = pCThits.Lyr[lyr].Y[0].at(0);
            }
            LineFit LineFitV(4, xHit, yHit);
            if (LineFitV.getChi2() < mxChi2) {
              double Vextrap = LineFitV.eval(geometry.uT(lyrA));
              if (Vextrap > -40. && Vextrap < 40.) {
                for (int lyr = 0; lyr < 3; lyr++) {
                  xHit[lyr] = pCThits.Lyr[lyrFit[lyr]].U[1].at(0);
                  yHit[lyr] = pCThits.Lyr[lyrFit[lyr]].Y[1].at(0);
                }
                LineFit LineFitT(3, xHit, yHit);
                double tPos = LineFitT.eval(geometry.uT(lyrA));
                if (LineFitT.getChi2() < mxChi2 && tPos > -155. &&
                    tPos < 155.) {
                  nTtrial[lyrA]++;
                  bool success = false;
                  for (int hit = 0; hit < pCThits.Lyr[lyrA].N[1]; hit++) {
                    double resid = pCThits.Lyr[lyrA].Y[1].at(hit) - tPos;
                    Tresidual[lyrA].entry(resid);
                    if (abs(resid) < hitWindow) {
                      nTsuccess[lyrA]++;
                      success = true;
                      break;
                    }
                  }
                  if (!success) { // We missed a hit on this layer. What does
                                  // the time since previous trigger look like?
                    deltaTmissT->entry(10.0 * (double)pCTEvent.delta_t);
                  }
                }
              }
            }
          }
        }
      }
    }

    // Look at the minimum time between hits on a given strip
    deltaTevt->entry(10 * (int)timeDiff);
    vector<int>::iterator it;
    for (int fpga = 0; fpga < num_tkr_fpga; fpga++) {
      int nChips = pCTEvent.tkr_fpga[fpga].num_chips;
      for (int chip = 0; chip < nChips; chip++) {
        int chipNum = pCTEvent.tkr_fpga[fpga].chip[chip].address;
        int nClus = pCTEvent.tkr_fpga[fpga].chip[chip].num_clusts;
        for (int clust = 0; clust < nClus; clust++) {
          int frstStrip =
              pCTEvent.tkr_fpga[fpga].chip[chip].cluster[clust].first;
          int len = pCTEvent.tkr_fpga[fpga].chip[chip].cluster[clust].length;
          for (int strip = 0; strip < len; strip++) {
            int address = frstStrip + strip;
            int code = address + 100 * chip + 10000 * fpga;
            for (it = hits.begin(); it != hits.end(); it++) {
              if (*it == code) {
                deltaTstrip->entry(10 * (int)timeDiff);
              }
            }
          }
        }
      }
    }
    hits.clear();
    for (int fpga = 0; fpga < num_tkr_fpga; fpga++) {
      int nChips = pCTEvent.tkr_fpga[fpga].num_chips;
      for (int chip = 0; chip < nChips; chip++) {
        int chipNum = pCTEvent.tkr_fpga[fpga].chip[chip].address;
        int nClus = pCTEvent.tkr_fpga[fpga].chip[chip].num_clusts;
        for (int clust = 0; clust < nClus; clust++) {
          int frstStrip =
              pCTEvent.tkr_fpga[fpga].chip[chip].cluster[clust].first;
          int len = pCTEvent.tkr_fpga[fpga].chip[chip].cluster[clust].length;
          for (int strip = 0; strip < len; strip++) {
            int address = frstStrip + strip;
            int code = address + 100 * chip + 10000 * fpga;
            hits.push_back(code);
          }
        }
      }
    }

    return cutEvent;

  }; // End of user event analysis

  void summary(std::string outputDir) { // User analysis summary called after
                                        // completing the loop over all events.
                                        // Write out results and print plots.

    cout << "\n UserAnalysis summary: maximum time tag encountered = "
         << maxTime << "\n";
    cout << " UserAnalysis summary: Minimum time since the last event, in 10ns "
            "units (negative indicates time tag error) = " << minTimeDiff
         << "\n\n";
    cout << " UserAnalysis summary: number of events analyzed = " << numEvents
         << "\n\n";
    cout << " UserAnalysis summary: number of events with at least one good "
            "track = " << nTkEvnts << "\n";
    cout << " UserAnalysis summary: number of events with multiple tracks in V "
            "and/or T = " << nMultTk << "\n";
    cout << " UserAnalysis summary: number of tracks through the phantom "
            "center = " << nInside << "\n";
    cout << " UserAnalysis summary: number of events with a front vector "
            "pointing to the phantom center = " << nVecInsideFront << "\n\n";
    cout << " UserAnalysis summary: number of events with a back vector "
            "pointing to the phantom center = " << nVecInsideBack << "\n\n";
    cout << " UserAnalysis summary: number of events with at least one perfect "
            "V track (i.e. with 4 hits) = " << nTkEvV << "\n";
    cout << " UserAnalysis summary: number of events with at least one perfect "
            "T track (i.e. with 4 hits) = " << nTkEvT << "\n\n";

    double VeffL[4] = { 0., 0., 0., 0. };
    double VeffR[4] = { 0., 0., 0., 0. };
    double Teff[4] = { 0., 0., 0., 0. };
    cout << "Chi^2 cut for the efficiency analysis = " << mxChi2 << endl;
    cout << "Residual cut for the efficiency analysis = " << hitWindow << " mm"
         << endl;
    for (lyrA = 0; lyrA < 4; lyrA++) {
      cout << "Number of V-layer efficiency trials for negative T=   "
           << nVtrialL[lyrA] << "\n";
      cout << "Number of V-layer efficiency successes for negative T="
           << nVsuccessL[lyrA] << "\n";
      if (nVtrialL > 0) {
        VeffL[lyrA] = (double)nVsuccessL[lyrA] / (double)nVtrialL[lyrA];
        cout << "Efficiency of V layer " << lyrA + 1 << " for negative T is "
             << VeffL[lyrA] << "\n";
      }
      cout << "Number of V-layer efficiency trials for positive T=   "
           << nVtrialR[lyrA] << "\n";
      cout << "Number of V-layer efficiency successes for positive T="
           << nVsuccessR[lyrA] << "\n";
      if (nVtrialR > 0) {
        VeffR[lyrA] = (double)nVsuccessR[lyrA] / (double)nVtrialR[lyrA];
        cout << "Efficiency of V layer " << lyrA + 1 << " for positive T is "
             << VeffR[lyrA] << "\n";
      }
      cout << "Number of T-layer efficiency trials=   " << nTtrial[lyrA]
           << "\n";
      cout << "Number of T-layer efficiency successes=" << nTsuccess[lyrA]
           << "\n";
      if (nTtrial > 0) {
        Teff[lyrA] = (double)nTsuccess[lyrA] / (double)nTtrial[lyrA];
        cout << "Efficiency of T layer " << lyrA + 1 << " is " << Teff[lyrA]
             << "\n";
      }
    }

    // Occupancy summary
    if (nEvntOcc > 0) {
      cout << endl << "Tracker channel occupancy anomaly summary:" << endl;
      int nchan = 0;
      int nhits = 0;
      for (int FPGA = 0; FPGA < num_tkr_fpga; FPGA++) {
        for (int chip = 0; chip < max_chips; chip++) {
          for (int ch = 0; ch < 64; ch++) {
            nchan++;
            nhits += tkrChCnt[FPGA][chip][ch];
            if (nEvntOcc > 1000000) {
              if (tkrChCnt[FPGA][chip][ch] == 0) {
                cout << "  FPGA " << FPGA << " chip " << chip << " channel "
                     << ch << " has zero occupancy and may be dead" << endl;
              }
            }
          }
        }
      }
      float avgOcc = float(nhits) / float(nchan) / float(nEvntOcc);
      cout << "  Average occupancy = " << avgOcc << endl;
      for (int FPGA = 0; FPGA < num_tkr_fpga; FPGA++) {
        for (int chip = 0; chip < max_chips; chip++) {
          for (int ch = 0; ch < 64; ch++) {
            float occ = float(tkrChCnt[FPGA][chip][ch]) / float(nEvntOcc);
            float focc = occ / avgOcc;
            if (focc > 3.5) {
              cout << "    FPGA " << FPGA << " chip " << chip << " channel "
                   << ch << " has occupancy=" << occ << ",   " << focc << endl;
            }
          }
        }
      }
    }

    char strbuf3[100];
    time_t t = time(NULL);
    struct tm *now = localtime(&t);
    sprintf(strbuf3, " %d-%d-%d", now->tm_year + 1900, now->tm_mon + 1,
            now->tm_mday);
    std::string FNpD = fname + strbuf3;

    FILE *oFile;
    string fileName;

    string setTerm;
    if (OsName == "Windows")
      setTerm = "set terminal wxt size 1300, 600\n";
    else
      setTerm = "set terminal x11 size 1300, 600\n";

    //  Plot time-stamp histogram
    fileName = outputDir + "/timeStamp.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Time stamps for file %s' layout 3,1 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      timeStamp->plot(oFile);
      hTdiff->plot(oFile);
      hErrInt->plot(oFile);
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //  Plot the trigger bit histogram
    fileName = outputDir + "/TrgBits.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      hTrig->plot(oFile, false, "", FNpD);
      fclose(oFile);
    }

    //  Tracker cluster size histogram
    fileName = outputDir + "/ClustSize.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      hClSz->plot(oFile, false, "", FNpD);
      fclose(oFile);
    }

    //  Plot histogram of track position at energy detector
    fileName = outputDir + "/TkrEpos.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Track at E-det %s' layout 3,1 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hVcalInt->plot(oFile);
      hTcalInt->plot(oFile);
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //  Plot the histogram of trigger spacing
    fileName = outputDir + "/deltaT.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      deltaT->plot(oFile, true, "", FNpD);
      deltaT->print("deltaT.txt");
      fclose(oFile);
    }
    fileName = outputDir + "/deltaTstrip.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Event delta-t %s' layout 1,2 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      deltaTevt->plot(oFile, true, "", FNpD);
      deltaTstrip->plot(oFile, true, "", FNpD);
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }
    fileName = outputDir + "/deltaTmiss.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Missed hits %s' layout 1,2 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      deltaTmissV->plot(oFile, true, "", FNpD);
      deltaTmissT->plot(oFile, true, "", FNpD);
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //    Plot histogram of trigger bits set

    fileName = outputDir + "/triggerBits.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      hEnrgTrg->plot(oFile, false, "", FNpD);
      fileName = outputDir + "/triggerBits.txt";
      hEnrgTrg->print(fileName);
      fclose(oFile);
    }

    //  Plot the raw energy detector data
    fileName = outputDir + "/EnergyDetector.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Energy Detector %s' layout 3,2 "
                     "columnsfirst scale 1.0,1.0\n",
              fname.c_str());
      PHedet[0].plot(oFile);
      PHedet[1].plot(oFile);
      PHedet[2].plot(oFile);
      PHedet[3].plot(oFile);
      fprintf(oFile, "set label 4 'Current date and time: %d-%d-%d    %d:%d' "
                     "at screen 0.55, 0.28\n",
              now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour,
              now->tm_min);
      fprintf(oFile, "set label 5 'Run Number %d' at screen 0.55, 0.25 left\n",
              runNumber);
      fprintf(
          oFile,
          "set label 6 'FPGA Program Version %d' at screen 0.55, 0.22 left\n",
          programVersion);
      fprintf(oFile, "set label 7 'Stage Angle = %5.1f degrees' at screen "
                     "0.55, 0.19 left\n",
              stageAngle);
      fprintf(oFile,
              "set label 8 'Start Time = %s ' at screen 0.55, 0.16 left\n",
              startTime.c_str());
      fprintf(oFile,
              "set label 9 'Number of events = %d' at screen 0.55, 0.13 left\n",
              numEvents);
      PHedet[4].plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/OTR.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Energy Detector ADC overflows for "
                     "file %s' layout 2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hOTR0->plot(oFile);
      hOTR->plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fclose(oFile);
    }

    //    Plot energy detector maximum samples
    fileName = outputDir + "/SampleMax.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Energy Detector Samples from %s' "
                     "layout 3,2 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hSmpMx[0].plot(oFile);
      hSmpMx[1].plot(oFile);
      hSmpMx[2].plot(oFile);
      hSmpMx[3].plot(oFile);
      hSmpMx[4].plot(oFile);
      hSmpMx[5].plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    // Plot all energy detector samples
    fileName = outputDir + "/SamplePH.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Energy Detector Samples from %s' "
                     "layout 3,2 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hSamplePH[0].plot(oFile, true);
      hSamplePH[1].plot(oFile, true);
      hSamplePH[2].plot(oFile, true);
      hSamplePH[3].plot(oFile, true);
      hSamplePH[4].plot(oFile, true);
      hSamplePH[5].plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //    Plot histograms of the energy detector response

    fileName = outputDir + "/PH1Profile.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Stage 1 ADC Profiles %s' layout 2,1 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      PHvsV[0].plot(oFile);
      PHvsT[0].plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/PH2Profile.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Stage 2 ADC Profiles %s' layout 2,1 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      PHvsV[1].plot(oFile);
      PHvsT[1].plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/PH3Profile.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Stage 3 ADC Profiles %s' layout 2,1 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      PHvsV[2].plot(oFile);
      PHvsT[2].plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/PH4Profile.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Stage 4 ADC Profiles %s' layout 2,1 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      PHvsV[3].plot(oFile);
      PHvsT[3].plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/PH5Profile.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Stage 5 ADC Profiles %s' layout 2,1 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      PHvsV[4].plot(oFile);
      PHvsT[4].plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/MxPh0.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Maximum Sample Pulse Height %s' "
                     "layout 3,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      MxPh0->plot(oFile, true);
      MxPh1->plot(oFile, true);
      MxPh2->plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }
    fileName = outputDir + "/MxPh1.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Maximum Sample Pulse Height %s' "
                     "layout 2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      MxPh3->plot(oFile, true);
      MxPh4->plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/MxSample.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Maximum Sample Location %s' layout "
                     "2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      MxSample0->plot(oFile, true);
      MxSample1->plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //  Plot pedestal histograms
    fileName = outputDir + "/FPGApedestals.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Energy pedestals %s' layout 3,2 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hPeds[0].plot(oFile, true);
      hPeds[1].plot(oFile, true);
      hPeds[2].plot(oFile, true);
      hPeds[3].plot(oFile, true);
      hPeds[4].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //  Plot histograms of cluster frequency

    fileName = outputDir + "/Vclusters.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Tracker V-Layer Cluster Frequency "
                     "%s' layout 2,2 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      nClustV[0].plot(oFile, true);
      nClustV[1].plot(oFile, true);
      nClustV[2].plot(oFile, true);
      nClustV[3].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/Tclusters.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Tracker T-Layer Cluster Frequency "
                     "%s' layout 2,2 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      nClustT[0].plot(oFile, true);
      nClustT[1].plot(oFile, true);
      nClustT[2].plot(oFile, true);
      nClustT[3].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //  Front tracker vectors

    fileName = outputDir + "/BeamSpot.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Front Vectors  %s' layout 2,1 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hVbeamSpot->plot(oFile);
      hTbeamSpot->plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //  Plot the tracking histograms
    fileName = outputDir + "/nTracks.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Number of Tracks Found %s' layout "
                     "2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hnTkV->plot(oFile, true);
      hnTkT->plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }
    fileName = outputDir + "/tMiss.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Track miss distance at u=0 %s' "
                     "layout 2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hTkMissV->plot(oFile);
      hTkMissT->plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }
    fileName = outputDir + "/tFrontSlope.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Front track vector slope %s' layout "
                     "2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hTkSlopeFV->plot(oFile);
      hTkSlopeFT->plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }
    fileName = outputDir + "/tRearSlope.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Rear track vector slope %s' layout "
                     "2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hTkSlopeRV->plot(oFile);
      hTkSlopeRT->plot(oFile);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //    Plot histograms of tracker strip hits

    fileName = outputDir + "/Cassette1.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Tracking Detector First Cassette  "
                     "%s' layout 2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      vHist[0].plot(oFile, true);
      tHist[0].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/Cassette2.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Tracking Detector Second Cassette  "
                     "%s' layout 2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      vHist[1].plot(oFile, true);
      tHist[1].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/Cassette3.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Tracking Detector Third Cassette  "
                     "%s' layout 2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      vHist[2].plot(oFile, true);
      tHist[2].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fclose(oFile);
    }

    fileName = outputDir + "/Cassette4.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Tracking Detector Fourth Cassette  "
                     "%s' layout 2,1 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      vHist[3].plot(oFile, true);
      tHist[3].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fclose(oFile);
    }

    // Plot the tracker residual distributions
    fileName = outputDir + "/Vresiduals.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Tracker V Residuals from %s' layout "
                     "2,2 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hResidV[0].plot(oFile, true);
      hResidV[1].plot(oFile, true);
      hResidV[2].plot(oFile, true);
      hResidV[3].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fclose(oFile);
    }
    fileName = outputDir + "/Tresiduals.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Tracker T Residuals from %s' layout "
                     "2,2 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hResidT[0].plot(oFile, true);
      hResidT[1].plot(oFile, true);
      hResidT[2].plot(oFile, true);
      hResidT[3].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fclose(oFile);
    }

    //  Plot histograms of tracker V residuals and print efficiencies

    fileName = outputDir + "/Vefficiency.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "binwidth=0.1\n");
      fprintf(oFile, "set multiplot title 'Tracker V-Layer Residuals %s' "
                     "layout 2,2 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      fprintf(oFile, "set label 1 'T<0 Eff.=%6.4f' at graph 0.05, 0.9 left\n",
              VeffL[0]);
      fprintf(oFile, "set label 2 'T>0 Eff.=%6.4f' at graph 0.05, 0.8 left\n",
              VeffR[0]);
      Vresidual[0].plot(oFile, true);
      fprintf(oFile, "set label 1 'T<0 Eff.=%6.4f' at graph 0.05, 0.9 left\n",
              VeffL[1]);
      fprintf(oFile, "set label 2 'T>0 Eff.=%6.4f' at graph 0.05, 0.8 left\n",
              VeffR[1]);
      Vresidual[1].plot(oFile, true);
      fprintf(oFile, "set label 1 'T<0 Eff.=%6.4f' at graph 0.05, 0.9 left\n",
              VeffL[2]);
      fprintf(oFile, "set label 2 'T>0 Eff.=%6.4f' at graph 0.05, 0.8 left\n",
              VeffR[2]);
      Vresidual[2].plot(oFile, true);
      fprintf(oFile, "set label 1 'T<0 Eff.=%6.4f' at graph 0.05, 0.9 left\n",
              VeffL[3]);
      fprintf(oFile, "set label 2 'T>0 Eff.=%6.4f' at graph 0.05, 0.8 left\n",
              VeffR[3]);
      Vresidual[3].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fprintf(oFile, "unset label 1 \n");
      fprintf(oFile, "unset label 2 \n");
      fclose(oFile);
    }

    // Plot histograms of T residuals and print efficiencies

    fileName = outputDir + "/Tefficiency.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "binwidth=0.1\n");
      fprintf(oFile, "set multiplot title 'Tracker T-Layer Residuals %s' "
                     "layout 2,2 columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      fprintf(oFile, "set label 3 'Efficiency=%6.4f' at graph 0.05, 0.9 left\n",
              Teff[0]);
      Tresidual[0].plot(oFile, true);
      fprintf(oFile, "set label 3 'Efficiency=%6.4f' at graph 0.05, 0.9 left\n",
              Teff[1]);
      Tresidual[1].plot(oFile, true);
      fprintf(oFile, "set label 3 'Efficiency=%6.4f' at graph 0.05, 0.9 left\n",
              Teff[2]);
      Tresidual[2].plot(oFile, true);
      fprintf(oFile, "set label 3 'Efficiency=%6.4f' at graph 0.05, 0.9 left\n",
              Teff[3]);
      Tresidual[3].plot(oFile, true);
      fprintf(oFile, "unset multiplot\n");
      fprintf(oFile, "show label\n");
      fprintf(oFile, "unset label\n");
      fclose(oFile);
    }

    //  Plot histogram of the stage angle
    fileName = outputDir + "/stageAngle.gp";
    oFile = fopen(fileName.c_str(), "w");
    if (oFile != NULL) {
      fprintf(oFile, setTerm.c_str());
      fprintf(oFile, "set multiplot title 'Stage angle %s' layout 2,1 "
                     "columnsfirst scale 1.0,1.0\n",
              FNpD.c_str());
      hTheta->plot(oFile);
      hEvsTheta->plot(oFile);
      fclose(oFile);
    }

    if (analysisLevel >= 2) {

      // Plot of WEPL vs position
      fileName = outputDir + "/TvsWEPL.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        hTedet->plot(oFile);
        fclose(oFile);
      }
      fileName = outputDir + "/VvsWEPL.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        hVedet->plot(oFile);
        fclose(oFile);
      }

      //  Plot the energy histogram
      fileName = outputDir + "/Etotal.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        hEcal->plot(oFile, true, "", FNpD);
        fclose(oFile);
      }

      //  Plot histogram of the deepest stage hit in the energy detector

      fileName = outputDir + "/maxLayer.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        hMxLyr->plot(oFile, true, "", FNpD);
        fclose(oFile);
      }

      // Plot histograms for stopping particles
      fileName = outputDir + "/StoppingEnergies.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        time_t t = time(NULL);
        struct tm *now = localtime(&t);
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Energy Detector Stage Stopping "
                       "Energies for file %s' layout 3,2 columnsfirst scale "
                       "1.0,1.0\n",
                fname.c_str());
        hEstop[0]->plot(oFile);
        hEstop[1]->plot(oFile);
        hEstop[2]->plot(oFile);
        hEstop[3]->plot(oFile);
        fprintf(oFile, "set label 4 'Current date and time: %d-%d-%d    %d:%d' "
                       "at screen 0.55, 0.28\n",
                now->tm_year + 1900, now->tm_mon + 1, now->tm_mday,
                now->tm_hour, now->tm_min);
        fprintf(oFile,
                "set label 5 'Run Number %d' at screen 0.55, 0.25 left\n",
                runNumber);
        fprintf(
            oFile,
            "set label 6 'FPGA Program Version %d' at screen 0.55, 0.22 left\n",
            programVersion);
        fprintf(oFile, "set label 7 'Stage Angle = %5.1f degrees' at screen "
                       "0.55, 0.19 left\n",
                stageAngle);
        fprintf(oFile,
                "set label 8 'Start Time = %s ' at screen 0.55, 0.16 left\n",
                startTime.c_str());
        fprintf(
            oFile,
            "set label 9 'Number of events = %d' at screen 0.55, 0.13 left\n",
            numEvents);
        hEstop[4]->plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fprintf(oFile, "show label\n");
        fprintf(oFile, "unset label\n");
        fclose(oFile);
      }

      //  Plot the scattering angle histogram
      fileName = outputDir + "/scatt.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        hScat->plot(oFile, true, "", FNpD);
        fclose(oFile);
      }

      //  Plot the WEPL histograms
      fileName = outputDir + "/WEPL.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Calibrated WEPL %s' layout 2,1 "
                       "columnsfirst scale 1.0,1.0\n",
                FNpD.c_str());
        hWEPL->plot(oFile);
        hWEPLzero->plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fprintf(oFile, "show label\n");
        fprintf(oFile, "unset label\n");
        fclose(oFile);
      }

      //  Plot the radiograph
      fileName = outputDir + "/RadioGraph.txt";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        hRadioGraph->print(oFile);
        fclose(oFile);
      }
      delete hRadioGraph;
      fileName = outputDir + "/RadioGraph2.txt";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        hRadioGraph2->print(oFile);
        fclose(oFile);
      }
      delete hRadioGraph2;

      fileName = outputDir + "/Energy1Profile.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Stage 1 Energy Profiles %s' "
                       "layout 2,1 columnsfirst scale 1.0,1.0\n",
                FNpD.c_str());
        EvsV[0].plot(oFile);
        EvsT[0].plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fprintf(oFile, "show label\n");
        fprintf(oFile, "unset label\n");
        fclose(oFile);
      }

      fileName = outputDir + "/Energy2Profile.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Stage 2 Energy Profiles %s' "
                       "layout 2,1 columnsfirst scale 1.0,1.0\n",
                FNpD.c_str());
        EvsV[1].plot(oFile);
        EvsT[1].plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fprintf(oFile, "show label\n");
        fprintf(oFile, "unset label\n");
        fclose(oFile);
      }

      fileName = outputDir + "/Energy3Profile.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Stage 3 Energy Profiles %s' "
                       "layout 2,1 columnsfirst scale 1.0,1.0\n",
                FNpD.c_str());
        EvsV[2].plot(oFile);
        EvsT[2].plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fprintf(oFile, "show label\n");
        fprintf(oFile, "unset label\n");
        fclose(oFile);
      }

      fileName = outputDir + "/Energy4Profile.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Stage 4 Energy Profiles %s' "
                       "layout 2,1 columnsfirst scale 1.0,1.0\n",
                FNpD.c_str());
        EvsV[3].plot(oFile);
        EvsT[3].plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fprintf(oFile, "show label\n");
        fprintf(oFile, "unset label\n");
        fclose(oFile);
      }

      fileName = outputDir + "/Energy5Profile.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Stage 5 Energy Profiles %s' "
                       "layout 2,1 columnsfirst scale 1.0,1.0\n",
                FNpD.c_str());
        EvsV[4].plot(oFile);
        EvsT[4].plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fprintf(oFile, "show label\n");
        fprintf(oFile, "unset label\n");
        fclose(oFile);
      }

      fileName = outputDir + "/Stage4.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Energies for particles stopping "
                       "in stage 4 %s' layout 2,2 columnsfirst scale 1.0,1.0\n",
                FNpD.c_str());
        hStop4[0]->plot(oFile);
        hStop4[1]->plot(oFile);
        hStop4[2]->plot(oFile);
        hStop4[3]->plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fclose(oFile);
      }
      fileName = outputDir + "/Stage3.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Energies for particles stopping "
                       "in stage 3 %s' layout 3,1 columnsfirst scale 1.0,1.0\n",
                FNpD.c_str());
        hStop3[0]->plot(oFile);
        hStop3[1]->plot(oFile);
        hStop3[2]->plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fclose(oFile);
      }
      fileName = outputDir + "/Stage2.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        fprintf(oFile, setTerm.c_str());
        fprintf(oFile, "set multiplot title 'Energies for particles stopping "
                       "in stage 2 %s' layout 2,1 columnsfirst scale 1.0,1.0\n",
                FNpD.c_str());
        hStop2[0]->plot(oFile);
        hStop2[1]->plot(oFile);
        fprintf(oFile, "unset multiplot\n");
        fclose(oFile);
      }
      fileName = outputDir + "/Stage1.gp";
      oFile = fopen(fileName.c_str(), "w");
      if (oFile != NULL) {
        hStop1->plot(oFile);
        fclose(oFile);
      }
    }


    cout << "Number of events with the 2nd stage trigger bits fired="
         << nTgrLyr1 << endl;
    cout << "Number of events with the 3rd stage trigger bits fired="
         << nTgrLyr2 << endl;
    cout << "Number of events with both 1st and 2nd stage trigger bits fired="
         << nTgrLyr01 << endl;
    cout << "Number of events with both 1st and 3rd stage trigger bits fired="
         << nTgrLyr02 << endl;
    if (nTgrLyr1 > 0) {
      float trigEff = float(nTgrLyr01) / float(nTgrLyr1);
      cout << "Efficiency of the layer-0 trigger for protons firing 2nd stage="
           << trigEff << endl;
    }
    if (nTgrLyr2 > 0) {
      float trigEff = float(nTgrLyr02) / float(nTgrLyr2);
      cout << "Efficiency of the layer-0 trigger for protons firing 3rd stage="
           << trigEff << endl;
    }
    cout << "Done with the user analysis summary\n";

  }; // End of user analysis summary

}; // End UserAnalysis class definition
#endif
