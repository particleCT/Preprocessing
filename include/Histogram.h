#ifndef Histogram_h
#define Histogram_h
// Simple classes to make histograms and profile plots, using Gnuplot to do the
// plotting
// R.P. Johnson   5/22/2016

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

class Histogram {
  int N;
  double BW, B0;
  vector<double> counts;
  string T, XL, YL;
  double sumX, sumX2;
  int nEntry;
  vector<double> X, Y, Ex, Ey;

public:
  Histogram(int nBins, double bin0, double binWidth, string title, string xlabel, string ylabel);
  string Title() { return T; }

  inline void entry(double x) {
    int bin = floor((x - B0) / BW);
    if (bin >= 0 && bin < N)
      counts[bin] += 1.0;
    if (x >= B0 && x <= B0 + N * BW) {
      sumX += x;
      sumX2 += x * x;
      nEntry++;
    }
  }
  inline void entry(float xf) {
    double x = (double)xf;
    int bin = floor((x - B0) / BW);
    if (bin >= 0 && bin < N)
      counts[bin] += 1.0;
    if (x >= B0 && x <= B0 + N * BW) {
      sumX += x;
      sumX2 += x * x;
      nEntry++;
    }
  }
  inline void entry(int i) {
    double x = double(i);
    int bin = floor((x - B0) / BW);
    if (bin >= 0 && bin < N)
      counts[bin] += 1.0;
    if (x >= B0 && x <= B0 + N * BW) {
      sumX += x;
      sumX2 += x * x;
      nEntry++;
    }
  }

  void plot(FILE *oFile, bool stats = false, string choice = "", std::string txt = " ");
  void print(string fn);
  vector<double> getBins();
  vector<double> getContents();
  double mode();
  int cnts(float min, float max);
  int cnts();
  double mean(float min, float max);
  double mean();
  double max();
  int FWHMboundaries(float &xLow, float &xHigh);
  bool setContents(vector<double> stuff);
  int imode();
};

// Class to make 2d histograms, using Gnuplot for display
class Histogram2D {
  int NX, NY;
  double BXW, BX0;
  double BYW, BY0;
  vector<vector<double> > counts;
  string T, XL, YL, ZL;
  int nEntry;

public:
  Histogram2D(int nXbins, double Xbin0, double XbinWidth, int nYbins, double Ybin0, double YbinWidth, string title,
              string xlabel, string ylabel, string zlabel);
  void entry(double x, double y);
  void entry(int i, int j);
  void plot(FILE *oFile);
  void read(string fn);
  bool add(Histogram2D *otherHist, string Title);
  void normalizeColumns(int ymin, int ymax, double area, double cut);
  Histogram rowSlice(int i);
  Histogram yProj();
  Histogram yProj(int bMin, int bMax);
  string Title() { return T; }
  void stats();
};

// Class to make "profile plots" in Gnuplot, i.e. a histogram that displays the
// mean and error on the mean in each bin.
// To do this, when processing data the sum and sum-squared have to be
// accumulated in addition to the number of counts, for each bin.
class ProfilePlot {
  int N;
  double BW, B0;
  string T, XL, YL;
  vector<int> counts;
  vector<double> sumY, sumY2;
  vector<double> X, Y, Ex, Ey;
  int nEntries;

public:
  ProfilePlot(int nBins, double bin0, double binWidth, string title, string xlabel, string ylabel);

  inline void entry(double x, double y) {
    int bin = floor((x - B0) / BW);
    if (bin < 0)
      bin = 0;
    if (bin >= N)
      bin = N - 1;
    counts[bin]++;
    sumY[bin] += y;
    sumY2[bin] += y * y;
    nEntries++;
  }
  inline void entry(float xf, float yf) {
    double x = (double)xf;
    double y = (double)yf;
    int bin = floor((x - B0) / BW);
    if (bin < 0)
      bin = 0;
    if (bin >= N)
      bin = N - 1;
    counts[bin]++;
    sumY[bin] += y;
    sumY2[bin] += y * y;
    nEntries++;
  }

  void plot(FILE *oFile);
  void print(std::string fn);
};

class ProfilePlot2D {
  int Nx, Ny;
  double xBW, xB0;
  double yBW, yB0;
  string T, XL, YL, ZL;
  vector<vector<int> > counts;
  vector<vector<double> > sumZ, sumZ2;
  int nEntries;

public:
  ProfilePlot2D(int nxBins, double xbin0, double xbinWidth, int nyBins, double ybin0, double ybinWidth, string title,
                string xlabel, string ylabel, string zlabel);

  inline void entry(double x, double y, double z) {
    int xbin = floor((x - xB0) / xBW);
    if (xbin < 0)
      xbin = 0;
    if (xbin >= Nx)
      xbin = Nx - 1;
    int ybin = floor((y - yB0) / yBW);
    if (ybin < 0)
      ybin = 0;
    if (ybin >= Ny)
      ybin = Ny - 1;
    counts[xbin][ybin]++;
    sumZ[xbin][ybin] += z;
    sumZ2[xbin][ybin] += z * z;
    nEntries++;
  }
  void plot(FILE *oFile);
  void print(FILE *oFile);
};

class ProfilePlot2Dpeak {
  int Nx, Ny, Nz;
  double xBW, xB0;
  double yBW, yB0;
  double zBW, zB0;
  string T, XL, YL, ZL;
  vector<vector<Histogram *> > counts;
  int nEntries;

public:
  ProfilePlot2Dpeak(int nxBins, double xbin0, double xbinWidth, int nyBins, double ybin0, double ybinWidth, int nzBins,
                    double zbin0, double zbinWidth, string title, string xlabel, string ylabel, string zlabel);
  ~ProfilePlot2Dpeak();
  void entry(double x, double y, double z);
  void print(FILE *oFile);
};

#endif
