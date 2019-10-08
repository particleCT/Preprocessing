#ifndef LineFit_h
#define LineFit_h
// This least squares line fit is not used in the pCT event processing but is
// useful for analyzing empty events.
// R. Johnson   5/22/2016
#include <iostream>
#include <cstdio>
#include <cmath>
class LineFit { // Least squares fit to a straight line
  double slope;
  double intercept;
  double chi2;

public:
  LineFit(int nPnts, double x[], double y[]) {
    double xbar = 0., ybar = 0., x2bar = 0., xybar = 0.;
    for (int i = 0; i < nPnts; ++i) {
      xbar += x[i];
      ybar += y[i];
      x2bar += x[i] * x[i];
      xybar += x[i] * y[i];
    }
    xbar = xbar / nPnts;
    ybar = ybar / nPnts;
    x2bar = x2bar / nPnts;
    xybar = xybar / nPnts;
    slope = (xybar - xbar * ybar) / (x2bar - xbar * xbar);
    intercept = ybar - slope * xbar;
    chi2 = 0.;
    for (int i = 0; i < nPnts; ++i) {
      chi2 += pow((intercept + slope * x[i] - y[i]) / 0.155,
                  2); // Fudged to make the mean ~1 for 1 d.o.f.
    }
    //        std::cout << "LineFit: ";
    //        for (int lyr; lyr<nPnts; lyr++) {cout << x[lyr] << " " << y[lyr]
    // << " ";}
    //        std::cout << "  chi2=" << chi2 << "\n";
  }
  inline double eval(double x) { return (intercept + slope * x); }
  inline double getSlope() { return slope; }
  inline double getIntercept() { return intercept; }
  inline double getChi2() { return chi2; }
};
#endif