#ifndef _PCTCUT_H_
#define _PCTCUT_H_
// Class for selecting events to pass on to the image reconstruction

using namespace std;

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

  pCTcut(int thread, int cut1 = 3, int cut2 = 5, int cut3 = 8) { // Class constructor called prior to the event loop
    n1track = 0; // Various counters for summarizing the number of events killed by cuts
    nLT8hits = 0;
    nGoodXtra = 0;
    nBack = 0;
    nKeep = 0;
    event_counter = 0;
    Thread = thread;

    // Definitions of the event selection cuts:
    mxXlyr = cut1;    // Maximum number of V layers with extra, unused hits; and the same for T layers
    mxTotXlyr = cut2; // Maximum total number of layers (V and T) with extra hits
    mxXhits = cut3;   // Maximum number of extra, unused hits allowed in the tracker

    minTkrs = 1; // Mininum number of good tracks
    maxTkrs = 1; // Maximum number of good tracks
  }

  // Call for each raw event after the tracking is completed
  bool cutEvt(bool userKill, pCT_Tracking &pCTtracks, TkrHits &pCThits) { 
    event_counter++;
    bool good = false;
    if (pCTtracks.nTracks >= minTkrs && pCTtracks.nTracks <= maxTkrs) { // Exactly 1 V track and 1 T track.  No 2-track events allowed.
      n1track++;
      int nXtraHits = 0; // Number of unused hits in the tracker
      int nLyrXtraV = 0; // Number of layers with extra hits in V
      int nLyrXtraT = 0; // Number of layers with extra hits in T
      for (int lyr = 0; lyr < 4; lyr++) {
        bool extra = false;
        for (int i = 0; i < pCThits.Lyr[lyr].N[0]; i++) {
          int tk = pCThits.Lyr[lyr].F[0].at(i);
          if (tk < 0) {
            nXtraHits++;
            extra = true;
          } else if (!pCTtracks.VTracks[tk].Good) {
            nXtraHits++;
            extra = true;
          }
        }
        if (extra)
          nLyrXtraV++;
        extra = false;
        for (int i = 0; i < pCThits.Lyr[lyr].N[1]; i++) {
          int tk = pCThits.Lyr[lyr].F[1].at(i);
          if (tk < 0) {
            nXtraHits++;
            extra = true;
          } else {
            if (!pCTtracks.TTracks[tk].Good) {
              nXtraHits++;
              extra = true;
            }
          }
        }
        if (extra)
          nLyrXtraT++;
      }
      if (nXtraHits < mxXhits) {
        nLT8hits++;
        if (nLyrXtraV < mxXlyr && nLyrXtraT < mxXlyr && nLyrXtraV + nLyrXtraT < mxTotXlyr) {
          nGoodXtra++;
          if (!userKill) { // Here we have a 'good' event for image
                           // reconstruction
            nKeep++;
            good = true;
          }
        }
      }
    } else { // Here, just select on the number of hits in the rear tracker
      good = true;
      for (int lyr = 2; lyr < 4; lyr++) {
        if (pCThits.Lyr[lyr].N[0] == 0 || pCThits.Lyr[lyr].N[0] > -1)
          good = false;
        if (pCThits.Lyr[lyr].N[1] == 0 || pCThits.Lyr[lyr].N[0] > -1)
          good = false;
      }
      if (good) {
        nBack++;
        nKeep++;
      }
    }
    return good;
  }

  void summary() { // Summary of the processing up to the point of selecting
                   // events based on tracking
    cout << "pCTcut thread " << Thread << ": number of raw events processed = " << event_counter << endl;
    cout << "pCTcut thread " << Thread << ": number of events with exactly 1 track = " << n1track << endl;
    cout << "pCTcut thread " << Thread << ": number events with less than " << mxXhits << " unused hits = " << nLT8hits
         << endl;
    cout << "pCTcut thread " << Thread << ": number after requiring number of V or T layers to have " << mxXlyr
         << " unused hits \n";
    cout << "        and number of V plus T layers to have < " << mxTotXlyr << " extra hits= " << nGoodXtra << endl;
    cout << "pCTcut thread " << Thread << endl;
    cout << "pCTcut thread " << Thread << ": final number remaining after the private user cuts = " << nKeep << endl;
  }
};
#endif
