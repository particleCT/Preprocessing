#ifndef TkrHits_h
#define TkrHits_h
// A class to find and store all of the pCT tracker hits for an event, starting
// from raw clusters
// R.P. Johnson   5/22/2016
#include <iostream>
#include <vector>
#include "pCTraw.h"
#include "pCTgeo.h"

class TkrHits {
 public:
  TkrHits(pCTraw &pCTEvent, const pCTgeo* Geometry, bool print);
  ~TkrHits();
  ofstream TkrLogFile;
  struct LyrHits {
    int N[2];                       // Number of hits in each of the V and T views for a given layer
    std::vector<double> Y[2], U[2]; // 0=V and 1=T
    std::vector<int> F[2];          // Track number; -1 if not used
  };

  LyrHits Lyr[4]; // Hit list for each layer.  0,1= front tracker    2,3= back
                  // tracker
  //  Convert the raw tracker strip information into coordinates
  
  // Method to print the list of hits
  void dumpHits(int eventNumber);

}; // End of the TkrHits class
#endif
