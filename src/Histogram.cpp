// Simple classes to make histograms and profile plots, using Gnuplot to do the
// plotting
// R.P. Johnson   5/22/2016

#include "Histogram.h"
#include "Util.h"

Histogram::Histogram(int nBins, double bin0, double binWidth, std::string title, std::string xlabel,
                     std::string ylabel) {
  N = nBins;
  BW = binWidth;
  B0 = bin0;
  T = title;
  XL = xlabel;
  YL = ylabel;
  sumX = 0.;
  sumX2 = 0.;
  nEntry = 0;
  counts.resize(nBins);
  for (int i = 0; i < nBins; ++i) {
    counts[i] = 0.0;
  }
}

void Histogram::plot(FILE *oFile, bool stats, std::string choice,
                     std::string txt) { // Create a Gnuplot file to display the
                                        // histogram
  X.resize(N);
  Y.resize(N);
  Ex.resize(N);
  Ey.resize(N);
  for (int i = 0; i < N; ++i) {
    X[i] = B0 + (((double)i) + 0.5) * BW;
    Y[i] = counts[i];
    Ex[i] = BW / 2.;
    Ey[i] = sqrt(counts[i]);
  }
  double mean = -999.;
  double rms = -999.;
  if (nEntry > 0) {
    mean = sumX / float(nEntry);
    rms = sqrt(sumX2 / nEntry - mean * mean);
  }
  fprintf(oFile, "#*** This file is intended to be displayed by gnuplot.\n");
  fprintf(oFile, "#*** Either double click on the file (works in Windows at least),\n");
  fprintf(oFile, "#*** or else start up gnuplot and use the load command to "
                 "display the plot.\n");
  // The following line will make the plot persist in linux when double clicking
  // on it, but it doesn't work in multiplot mode.
  // fprintf(oFile,"set term X11 persist\n");
  if (stats) {
    fprintf(oFile, "set label 777 'mean=%7.3f' at graph 0.67, 0.9 left font "
                   "'Verdana,12'\n",
            mean);
    fprintf(oFile, "set label 778 'rms=%8.3f' at graph 0.67, 0.85 left font "
                   "'Verdana,12'\n",
            rms);
    fprintf(oFile, "set label 779 'counts=%d' at graph 0.67, 0.80 left font "
                   "'Verdana,12'\n",
            nEntry);
  }
  fprintf(oFile, "set xtics font 'Verdana,12'\n");
  fprintf(oFile, "set ytics font 'Verdana,12'\n");
  fprintf(oFile, "set title '%s %s' font 'Verdana,12'\n", T.c_str(), txt.c_str());
  fprintf(oFile, "set xlabel '%s' font 'Verdana,12'\n", XL.c_str());
  fprintf(oFile, "set ylabel '%s' font 'Verdana,12'\n", YL.c_str());
  fprintf(oFile, "set xrange[%7.4f : %7.4f]\n", B0, B0 + N * BW);
  fprintf(oFile, "set nokey\n");
  if (choice == "errors") {
    fprintf(oFile, "plot '-' with xyerrorbars \n");
    for (int i = 0; i < N; i++) {
      fprintf(oFile, "%8.3e %8.3e %8.3e %8.3e\n", X[i], Y[i], Ex[i], Ey[i]);
    }
    fprintf(oFile, "e\n");
  } else {
    fprintf(oFile, "plot '-' with boxes\n");
    for (int i = 0; i < N; i++) {
      fprintf(oFile, "%8.3e %8.3e\n", X[i], Y[i]);
    }
    fprintf(oFile, "e\n");
  }
}

bool Histogram::setContents(vector<double> stuff) {
  if (stuff.size() != counts.size()) {
    cout << "Histogram::setContents, mismatch of vector sizes " << endl;
    return false;
  }
  for (int i = 0; i < counts.size(); ++i) {
    counts[i] = stuff[i];
    nEntry += stuff[i];
    float x = B0 + (((double)i) + 0.5) * BW;
    sumX += ((double)stuff[i]) * x;
    sumX2 += ((double)stuff[i]) * x * x;
  }
  return true;
};

void Histogram::print(std::string fn) { // Print the histogram contents to an
                                        // ASCII file
  std::cout << "printing histogram " << T.c_str() << " to file " << fn << "\n";
  std::cout << "Number of bins = " << N << "\n";

  double mean = -999.;
  double rms = -999.;
  if (nEntry > 0) {
    mean = sumX / float(nEntry);
    rms = sqrt(sumX2 / nEntry - mean * mean);
  }
  FILE *oFile = fopen(fn.c_str(), "w");
  if (oFile != NULL) {
    fprintf(oFile, "Printing Histogram %s: %s vs %s\n", T.c_str(), YL.c_str(), XL.c_str());
    fprintf(oFile, "  Number of entries=%d\n", nEntry);
    fprintf(oFile, "  Mean=%e\n", mean);
    fprintf(oFile, "  rms= %e\n", rms);
    for (int i = 0; i < N; ++i) {
      double XX, EEy;
      double YY;
      XX = B0 + (((double)i) + 0.5) * BW;
      YY = counts[i];
      EEy = sqrt(counts[i]);
      if (BW > 0.1) {
        //fprintf(oFile, "  %d  %8.3f  %9.3 +0 %8.3f\n", i, XX, YY, EEy);
      } else {
        fprintf(oFile, "  %d  %e  %e\n", i, XX, YY);
      }
    }
    fclose(oFile);
  }
};

std::vector<double> Histogram::getBins() { // Return an array of the bin locations
  std::vector<double> tmp;
  for (int i = 0; i < N; ++i) {
    tmp.push_back(B0 + (((double)i) + 0.5) * BW);
  }
  return tmp;
};

std::vector<double> Histogram::getContents() { // Return an array of the contents of the bins
  std::vector<double> tmp;
  for (int i = 0; i < N; ++i) {
    tmp.push_back(counts[i]);
  }
  return tmp;
};

double Histogram::max() {
  double maxC = -1.0;
  for (int i = 0; i < N; ++i) {
    if (counts[i] > maxC) {
      maxC = counts[i];
    }
  }
  return maxC;
}

double Histogram::mode() { // Return the mode of the distribution
  double max = -1.0;
  int Imax = 0;
  for (int i = 0; i < N; ++i) {
    if (counts[i] > max) {
      max = counts[i];
      Imax = i;
    }
  }
  return B0 + (((double)Imax) + 0.5) * BW;
}

int Histogram::imode() { // Return the bin number for the mode of the distribution
  double max = -1.0;
  int Imax = 0;
  for (int i = 0; i < N; ++i) {
    if (counts[i] > max) {
      max = counts[i];
      Imax = i;
    }
  }
  return Imax;
}

int Histogram::FWHMboundaries(float &xLow, float &xHigh) {
  float max = -1.0;
  int Imax = -1;
  // Find the max bin and binId
  for (int i = 1; i < N - 1; ++i) {
    if (counts[i] > max) {
      // Try to avoid spurious single-bin 'noise' peaks
      //if (counts[i - 1] < 0.05*counts[i] && counts[i +  1] < 0.05*counts[i]) continue;
      max = counts[i];
      Imax = i;
    }
  }

  //binId < 0 
  if (Imax < 0) return -1; // safety check

  // higher bin Id
  int iUp = -1;
  for (int i = Imax; i < N; ++i) { // higher bound
    if (counts[i] < max / 2.)
      {
	iUp = i;
	break;
      }
    if (i == N - 1 && Imax != N - 2)
      {
	iUp = i;
	break;
      }
  }
  if (iUp < 0) return -2; 
    

  int iDn = -1;
  for (int i = Imax; i >= 0; --i) { // lower bound
    if (counts[i] < max / 2.) {
      iDn = i;
      break;
    }
    if (i == 0 && Imax != 1) {
      iDn = i;
      break;
    }
  }

  
  if (iDn < 0) return -3;
    

  xLow = B0 + (((double)iDn) + 0.5) * BW;
  xHigh = B0 + (((double)iUp) + 0.5) * BW;

  return 0;
}

int Histogram::cnts(float min, float max) { // Return counts within a range
  double tCnts = 0;
  for (int i = 0; i < N; ++i) {
    double XX = B0 + (((double)i) + 0.5) * BW;
    if (XX >= min && XX <= max) {
      tCnts += counts[i];
    }
  }
  return (int)tCnts;
}

int Histogram::cnts() {
  double tCnts = 0;
  for (int i = 0; i < N; ++i) {
    tCnts += counts[i];
  }
  return (int)tCnts;
}

double Histogram::mean(float min, float max) { // Return the mean within a range
  double tCnts = 0;
  double W = 0.;
  for (int i = 0; i < N; ++i) {
    double XX = B0 + (((double)i) + 0.5) * BW;
    if (XX >= min && XX <= max) {
      tCnts += counts[i];
      W += XX * ((double)counts[i]);
    }
  }
  double avg = 0.;
  if (tCnts > 0.0)
    avg = W / (tCnts);
  return avg;
}

double Histogram::mean() { // Return the mean of the full distribution
  double avg = 0.;
  if (nEntry > 0)
    avg = sumX / float(nEntry);
  return avg;
}

Histogram2D::Histogram2D(int nXbins, double Xbin0, double XbinWidth, int nYbins, double Ybin0, double YbinWidth,
                         string title, string xlabel, string ylabel, string zlabel) {
  NX = nXbins;
  NY = nYbins;
  BXW = XbinWidth;
  BYW = YbinWidth;
  BX0 = Xbin0;
  BY0 = Ybin0;
  T = title;
  XL = xlabel;
  YL = ylabel;
  ZL = zlabel;
  nEntry = 0;
  counts.resize(nXbins, vector<double>(nYbins, 0));
};

// Method to normalize all columns to the same area
// This is constructed specifically to aid in the pCT WEPL calibration
void Histogram2D::normalizeColumns(int ymin, int ymax, double area, double cut) {
  cout << "Normalizing the columns of 2D histogram " << T << " to equal areas " << area << endl;
  cout << "      Columns with area less than " << cut << " will be left alone." << endl;
  cout << "      Columns without a clear broad peak will be left alone." << endl;
  double *peaks = new double[NX];
  for (int col = 0; col < NX; ++col) {
    double A = 0.;
    int NNY = 0;
    double max = 0.;
    for (int row = 0; row < NY; ++row) {
      if (row < ymin)
        continue;
      if (row > ymax)
        break;
      NNY++;
      A += counts[col][row];
      if (counts[col][row] > max)
        max = counts[col][row];
    }
    peaks[col] = max;
    double AAvg = A / (double)NNY;
    double Aabove = 0.;
    int nConsec = 0;
    bool goodRow = false;
    for (int row = 0; row < NY; ++row) { // Noise filter
      if (row < ymin)
        continue;
      if (row > ymax)
        break;
      if (counts[col][row] > 1.5 * AAvg) {
        ++nConsec;
        Aabove += counts[col][row];
        if (nConsec > 9) {
          if (Aabove / A > 0.25) {
            goodRow = true;
            break;
          }
        }
      } else {
        nConsec = 0;
        Aabove = 0.;
      }
    }
    if (A == 0. || A < cut || !goodRow) {
      cout << "   Column " << col << ", total counts = " << A << " max=" << max << " avg=" << AAvg
           << " nConsec=" << nConsec << " Aabove=" << Aabove << endl;
      continue;
    }
    double norm = area / A;
    cout << "   Column " << col << ", total counts = " << A << " max=" << max << " avg=" << AAvg
         << " nConsec=" << nConsec << " Aabove=" << Aabove << " norm=" << norm << endl;
    for (int row = 0; row < NY; ++row) {
      counts[col][row] *= norm;
    }
    peaks[col] *= norm;
    if (col > 1) { // look for columns that are too low in their peak value and adjust
      if (peaks[col - 1] < 0.9 * peaks[col] && peaks[col - 1] < 0.9 * peaks[col - 2]) {
        double adj = 0.5 * (peaks[col - 2] + peaks[col]) / peaks[col - 1];
        for (int row = 0; row < NY; ++row) {
          counts[col - 1][row] *= adj;
        }
      }
    }
  }
  delete[] peaks;
}

void Histogram2D::entry(double x, double y) {
  int binx = floor((x - BX0) / BXW);
  int biny = floor((y - BY0) / BYW);
  if (binx >= 0 && binx < NX) {
    if (biny >= 0 && biny < NY) {
      counts[binx][biny] += 1.0;
      nEntry++;
    }
  }
};
void Histogram2D::entry(int i, int j) {
  double x = float(i);
  double y = float(j);
  int binx = floor((x - BX0) / BXW);
  int biny = floor((x - BY0) / BYW);
  if (binx >= 0 && binx < NX) {
    if (biny >= 0 && biny < NY) {
      counts[binx][biny] += 1.0;
      nEntry++;
    }
  }
};

// Add the contents of another histogram to those of this histogram
bool Histogram2D::add(Histogram2D *otherHist, string Title) {
  if (NX != otherHist->NX || NY != otherHist->NY) {
    cout << "Histogram2D::add, bin numbers do not match, cannot add" << endl;
    return false;
  }
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NY; ++j) {
      counts[i][j] += otherHist->counts[i][j];
    }
  }
  nEntry += otherHist->nEntry;
  T = Title;

  return true;
};

void Histogram2D::stats() {
  cout << endl << endl << "Statistics for histogram " << T << endl;
  cout << "Number of entries = " << nEntry << endl;
  double xMean = 0.;
  double yMean = 0.;
  for (int i = 0; i < NX; ++i) {
    double x = BX0 + i * BXW;
    for (int j = 0; j < NY; ++j) {
      double y = BY0 + j * BXW;
      xMean += x * counts[i][j];
      yMean += y * counts[i][j];
    }
  }
  if (nEntry > 0) {
    xMean = xMean / float(nEntry);
    yMean = yMean / float(nEntry);
  }
  cout << "Average X value = " << xMean << endl;
  cout << "Average Y value = " << yMean << endl << endl << endl;
};

// Return a row (constant y) of the 2D histogram as a 1D histogram:
Histogram Histogram2D::rowSlice(int i) {
  string Title;
  long long int slice = i;
  Title = T + " slice " + to_string(slice);
  Histogram tmpHst(NX, BX0, BXW, Title, XL, "n");

  vector<double> tmp;
  tmp.resize(NX);
  for (int j = 0; j < NX; ++j) {
    tmp[j] = counts[j][i];
  }
  tmpHst.setContents(tmp);

  return tmpHst;
};

// Return a Y projection of the 2D histogram as a 1D histogram:
Histogram Histogram2D::yProj() {
  string Title;
  Title = T + " y projection";
  Histogram tmpHst(NY, BY0, BYW, Title, YL, "y");

  vector<double> tmp;
  tmp.resize(NY);
  for (int j = 0; j < NY; ++j) {
    double cnt = 0.;
    for (int i = 0; i < NX; ++i) {
      cnt += counts[i][j];
    }
    tmp[j] = cnt;
  }
  tmpHst.setContents(tmp);

  return tmpHst;
};

Histogram Histogram2D::yProj(int bMin, int bMax) {
  string Title;
  Title = T + " y projection";
  Histogram tmpHst(NY, BY0, BYW, Title, YL, "y");

  vector<double> tmp;
  tmp.resize(NY);
  for (int j = 0; j < NY; ++j) {
    double cnt = 0.;
    for (int i = 0; i < NX; ++i) {
      if (i < bMin)
        continue;
      if (i > bMax)
        continue;
      cnt += counts[i][j];
    }
    tmp[j] = cnt;
  }
  tmpHst.setContents(tmp);

  return tmpHst;
};

void Histogram2D::read(string fn) {
  ifstream infile(fn);
  if (infile) {
    string line;
    while (getline(infile, line)) {
      size_t found = line.find("splot");
      if (found != line.npos) {
        cout << "Histogram2D::read: line " << line << " found in file " << fn << endl;
        break;
      }
    }

    for (int i = 0; i < NX; ++i) {
      for (int j = 0; j < NY; ++j) {
        if (getline(infile, line)) {
          vector<string> tokens = Util::getTokens(line);
          if (tokens.size() >= 3) {
            float X = stof(tokens[0]);
            float Y = stof(tokens[1]);
            float Z = stof(tokens[2]);
            counts[i][j] = Z;
          }
        } else {
          cout << "Histogram2D::read: a line is missing in file " << fn << " at i=" << i << " j=" << j << endl;
        }
      }
    }
  } else {
    cout << "Histogram2D::read: cannot find or open the input file " << fn << endl;
  }
}

void Histogram2D::plot(FILE *oFile) {

  fprintf(oFile, "#*** This file is intended to be displayed by gnuplot.\n");
  fprintf(oFile, "#*** Either double click on the file (works in Windows at least),\n");
  fprintf(oFile, "#*** or else start up gnuplot and use the load command to "
                 "display the plot.\n");
  // The following line will make the plot persist in linux when double clicking
  // on it, but it doesn't work in multiplot mode.
  // fprintf(oFile,"set term X11 persist\n");
  fprintf(oFile, "set title '%s' \n", T.c_str());
  fprintf(oFile, "set xlabel '%s' \n", XL.c_str());
  fprintf(oFile, "set ylabel '%s' \n", YL.c_str());
  fprintf(oFile, "set zlabel '%s' \n", ZL.c_str());
  fprintf(oFile, "#*** Modify these ranges by hand for better viewing\n");
  fprintf(oFile, "set xrange[%7.4f : %7.4f]\n", BX0, BX0 + NX * BXW);
  fprintf(oFile, "set yrange[%7.4f : %7.4f]\n", BY0, BY0 + NY * BYW);
  fprintf(oFile, "set zrange[0.0 : ]\n");
  fprintf(oFile, "set nokey\n");
  fprintf(oFile, "#*** Hidden 3d makes it slow to plot and manipulate the "
                 "graph; not worth it usually.\n");
  fprintf(oFile, "#set hidden3d\n");
  fprintf(oFile, "splot '-' using 1:2:3 with lines  \n");
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NY; ++j) {
      float X = BX0 + (((double)i) + 0.5) * BXW;
      float Y = BY0 + (((double)j) + 0.5) * BYW;
      float Z = counts[i][j];
      fprintf(oFile, "%8.3e %8.3e %8.3e\n", X, Y, Z);
    }
  }
  fprintf(oFile, "e\n");
};

ProfilePlot2D::ProfilePlot2D(int nxBins, double xbin0, double xbinWidth, int nyBins, double ybin0, double ybinWidth,
                             string title, string xlabel, string ylabel, string zlabel) {
  Nx = nxBins;
  Ny = nyBins;
  xB0 = xbin0;
  yB0 = ybin0;
  xBW = xbinWidth;
  yBW = ybinWidth;
  counts.resize(nxBins);
  sumZ.resize(nxBins);
  sumZ2.resize(nxBins);
  for (int i = 0; i < nxBins; i++) {
    counts[i].resize(nyBins, 0.);
    sumZ[i].resize(nyBins, 0.);
    sumZ2[i].resize(nyBins, 0.);
  }
  T = title;
  XL = xlabel;
  YL = ylabel;
  ZL = zlabel;
}

void ProfilePlot2D::plot(FILE *oFile) {

  fprintf(oFile, "#*** This file is intended to be displayed by gnuplot.\n");
  fprintf(oFile, "#*** Either double click on the file (works in Windows at least),\n");
  fprintf(oFile, "#*** or else start up gnuplot and use the load command to "
                 "display the plot.\n");
  fprintf(oFile, "set title '%s' \n", T.c_str());
  fprintf(oFile, "set xlabel '%s' \n", XL.c_str());
  fprintf(oFile, "set ylabel '%s' \n", YL.c_str());
  fprintf(oFile, "set zlabel '%s' \n", ZL.c_str());
  fprintf(oFile, "set xrange[%7.4f : %7.4f]\n", xB0, xB0 + Nx * xBW);
  fprintf(oFile, "set yrange[%7.4f : %7.4f]\n", yB0, yB0 + Ny * yBW);
  fprintf(oFile, "#set zrange[0.0 : ]\n");
  fprintf(oFile, "set nokey\n");
  fprintf(oFile, "#*** Hidden 3d makes it slow to plot and manipulate the "
                 "graph; not worth it usually.\n");
  fprintf(oFile, "#set hidden3d\n");
  fprintf(oFile, "$map << EOD\n");
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      double X = xB0 + (((double)i) + 0.5) * xBW;
      double Y = yB0 + (((double)j) + 0.5) * yBW;
      double Z;
      if (counts[i][j] > 0)
        Z = sumZ[i][j] / (double)counts[i][j];
      else
        Z = 0.;
      // double Ez;
      // if (counts[i][j] > 1) Ez= sqrt((sumZ2[i][j]/(double)counts[i][j] -
      // Z*Z)/((double)(counts[i][j]-1))); else Ez = 0.;
      if (counts[i][j] > 4)
        fprintf(oFile, "%8.3e %8.3e %8.3e\n", X, Y, Z);
    }
    fprintf(oFile, " \n");
  }
  fprintf(oFile, "EOD\n");
  fprintf(oFile, "plot '$map' using 1:2:3 with image pixels\n");
}

void ProfilePlot2D::print(FILE *oFile) {
  fprintf(oFile, "%s  Number of entries=%d\n", T.c_str(), nEntries);
  for (int j = 0; j < Ny; ++j) {
    double Y = yB0 + (((double)j) + 0.5) * yBW;
    fprintf(oFile, "%8.3e ", Y);
  }
  fprintf(oFile, " \n");
  for (int i = 0; i < Nx; ++i) {
    double X = xB0 + (((double)i) + 0.5) * xBW;
    fprintf(oFile, "%8.3e ", X);
    for (int j = 0; j < Ny; ++j) {
      double Y = yB0 + (((double)j) + 0.5) * yBW;
      double Z;
      if (counts[i][j] > 0)
        Z = sumZ[i][j] / (double)counts[i][j];
      else
        Z = 0.;
      fprintf(oFile, "%8.3e ", Z);
    }
    fprintf(oFile, " \n");
  }
}

ProfilePlot::ProfilePlot(int nBins, double bin0, double binWidth, std::string title, std::string xlabel,
                         std::string ylabel) {
  N = nBins;
  BW = binWidth;
  B0 = bin0;
  T = title;
  XL = xlabel;
  YL = ylabel;
  counts.resize(nBins);
  sumY.resize(nBins);
  sumY2.resize(nBins);
  for (int i = 0; i < nBins; ++i) {
    counts[i] = 0;
    sumY[i] = 0.;
    sumY2[i] = 0.;
  }
}

void ProfilePlot::print(std::string fn) {
  std::cout << "printing profile plot " << T.c_str() << " to file " << fn << "\n";
  std::cout << "Number of bins = " << N << "\n";

  FILE *oFile = fopen(fn.c_str(), "w");
  if (oFile != NULL) {
    fprintf(oFile, "Printing profile plot %s: %s vs %s\n", T.c_str(), YL.c_str(), XL.c_str());
    fprintf(oFile, "  Number of entries=%d\n", nEntries);
    for (int i = 0; i < N; ++i) {
      double XX, EEy, YY;
      int C;
      XX = B0 + (((double)i) + 0.5) * BW;
      C = counts[i];
      if (counts[i] > 0)
        YY = sumY[i] / (double)counts[i];
      else
        YY = 0.;
      if (counts[i] > 1)
        EEy = sqrt((sumY2[i] / (double)counts[i] - YY * YY) / ((double)(counts[i] - 1)));
      else
        EEy = 0.;
      fprintf(oFile, "  %d  %8.3f %d  %8.3f+-%8.3f\n", i, XX, C, YY, EEy);
    }
    fclose(oFile);
  }
}

void ProfilePlot::plot(FILE *oFile) {
  X.resize(N);
  Y.resize(N);
  Ex.resize(N);
  Ey.resize(N);
  for (int i = 0; i < N; ++i) {
    if (counts[i] > 0) {
      X[i] = B0 + (((double)i) + 0.5) * BW;
      Y[i] = sumY[i] / (double)counts[i];
      Ex[i] = BW / 2.;
      if (counts[i] > 1)
        Ey[i] = sqrt((sumY2[i] / (double)counts[i] - Y[i] * Y[i]) / ((double)(counts[i] - 1)));
      else
        Ey[i] = 0.;
    } else {
      X[i] = B0 + (((double)i) + 0.5) * BW;
      Y[i] = 0.;
      Ex[i] = 0.;
      Ey[i] = 0.;
    }
  }
  fprintf(oFile, "#*** This file is intended to be displayed by gnuplot.\n");
  fprintf(oFile, "#*** Either double click on the file (works in Windows at least),\n");
  fprintf(oFile, "#*** or else start up gnuplot and use the load command to "
                 "display the plot.\n");
  fprintf(oFile, "set title '%s' \n", T.c_str());
  fprintf(oFile, "set xlabel '%s' \n", XL.c_str());
  fprintf(oFile, "set ylabel '%s' \n", YL.c_str());
  fprintf(oFile, "set xrange[%7.4f : %7.4f]\n", B0, B0 + N * BW);
  fprintf(oFile, "set nokey\n");
  fprintf(oFile, "plot '-' with xyerrorbars \n");
  for (int i = 0; i < N; i++) {
    fprintf(oFile, "%8.3e %8.3e %8.3e %8.3e\n", X[i], Y[i], Ex[i], Ey[i]);
  }
  fprintf(oFile, "e\n");
}

// This class used a histogram in each bin of the 2D space to find the peak of
// the distribution in the bin, instead of the average value used in the normal
// profile plot
// Consequently it uses a lot of memory but works well for the radiograph
ProfilePlot2Dpeak::ProfilePlot2Dpeak(int nxBins, double xbin0, double xbinWidth, int nyBins, double ybin0,
                                     double ybinWidth, int nzBins, double zbin0, double zbinWidth, string title,
                                     string xlabel, string ylabel, string zlabel) {
  Nx = nxBins;
  Ny = nyBins;
  Nz = nzBins;
  xB0 = xbin0;
  yB0 = ybin0;
  zB0 = zbin0;
  xBW = xbinWidth;
  yBW = ybinWidth;
  zBW = zbinWidth;
  counts.resize(nxBins);

  for (int i = 0; i < nxBins; i++) {
    counts[i].resize(nyBins);
  }
  T = title;
  XL = xlabel;
  YL = ylabel;
  ZL = zlabel;
  nEntries = 0;

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      string title = T + " bin " + to_string((long double)i) + "," + to_string((long double)j);
      counts[i][j] = new Histogram(Nz, zB0, zBW, T, ZL, "Counts");
    }
  }
}

ProfilePlot2Dpeak::~ProfilePlot2Dpeak() {
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      delete counts[i][j];
    }
  }
}

void ProfilePlot2Dpeak::entry(double x, double y, double z) {
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
  counts[xbin][ybin]->entry(z);
  nEntries++;
}
void ProfilePlot2Dpeak::print(FILE *oFile) {
  fprintf(oFile, "%s  Number of entries=%d\n", T.c_str(), nEntries);
  for (int j = 0; j < Ny; ++j) {
    double Y = yB0 + (((double)j) + 0.5) * yBW;
    fprintf(oFile, "%8.3e ", Y);
  }
  fprintf(oFile, " \n");
  for (int i = 0; i < Nx; ++i) {
    double X = xB0 + (((double)i) + 0.5) * xBW;
    fprintf(oFile, "%8.3e ", X);
    for (int j = 0; j < Ny; ++j) {
      double Y = yB0 + (((double)j) + 0.5) * yBW;
      double Z;
      float xLow, xHigh;
      int ret = counts[i][j]->FWHMboundaries(xLow, xHigh);
      if (ret == 0) {
        Z = counts[i][j]->mean(xLow, xHigh);
      } else {
        float mode = counts[i][j]->mode();
        Z = counts[i][j]->mean(mode - 2. * zBW, mode + 2 * zBW);
      }
      fprintf(oFile, "%8.3e ", Z);
      /*			if (i>=40 && i<=59 && j==40) {
                                      string fn = "hist_" + to_string((long long
         int)i) + "_" + to_string((long long int)j) + ".gp";
                                      FILE* dFile= fopen(fn.c_str(),"w");
                                      counts[i][j]->plot(dFile);
                                      float mode = counts[i][j]->mode();
                                      float mean = counts[i][j]->mean();
                                      cout << "ProfilePlot2Dpeak::print: hist "
         << i << "," << j <<"100,80, xLow=" << xLow << " xHigh=" << xHigh << "
         Z=" << Z << " mode=" << mode << " mean= " << mean << endl;
                              }  */
    }
    fprintf(oFile, " \n");
  }
}
