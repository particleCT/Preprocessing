#ifndef _PCTCUT_H_
#define _PCTCUT_H_
// Class for selecting events to pass on to the image reconstruction
#include "TkrHits.h"
#include "pCT_Tracking.h"
#include "pCTconfig.h"
using namespace std;
class TkrHits;
class pCT_Tracking;
class pCTcut {
  int event_counter;
  int n1track; // Various counters for summarizing the number of events killed by cuts
  int nLT8hits;
  int nGoodXtra;
  int nBack;
  int mxXlyr;
  int mxTotXlyr;
  int mxXhits;
  int Thread;
  int minTkrs;
  int maxTkrs;

 public:
  int nKeep;
  // constructor
  pCTcut( pCTconfig cfg); // Class constructor called prior to the event loop
  pCTconfig config;
  // classes
  bool cutEvt(bool userKill, pCT_Tracking&, TkrHits&);
  bool EnrgCut(float [5], float, float, float, float );
    
  void summary();
};
#endif
