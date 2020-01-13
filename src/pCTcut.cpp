#include "pCTcut.h"

pCTcut::pCTcut(pCTconfig cfg): config(cfg)
{ // Class constructor called prior to the event loop
  n1track = 0; // Various counters for summarizing the number of events killed by cuts
  nLT8hits = 0;
  nGoodXtra = 0;
  nBack = 0;
  nKeep = 0;
  event_counter = 0;
  Thread = 0;

  // Definitions of the event selection cuts:
  mxXlyr = 3;    // Maximum number of V layers with extra, unused hits; and the same for T layers
  mxTotXlyr = 5; // Maximum total number of layers (V and T) with extra hits
  mxXhits = 8;   // Maximum number of extra, unused hits allowed in the tracker

  minTkrs = 1; // Mininum number of good tracks
  maxTkrs = 1; // Maximum number of good tracks
}


////////////////////////////////////////////////////////////////////
// Energy cuts
////////////////////////////////////////////////////////////////////
bool pCTcut::EnrgCut(float Estage[5], float Etot, float cut0, float cut1, float cut3) {

  bool dropEvent = false;
  // Cut on energy deposit not compatible with Bragg-Peak (cut0, cut1) adapted from Vladimir Bashkirov
  // For helium: Use the stages as dE-E detector and check if the particle's energy
  // deposit is compatible with the parameterized stage response
  // (Last modified: Lennart Volz, December 2017)

  if (Estage[4] >  config.item_float["thr4"]) { // Particule stop in stage 4
    if (Estage[3] < cut0 || Estage[2] < cut0 || Estage[1] < cut0 || Estage[0] < cut0 || Estage[4] > cut3) dropEvent = true;
    /*if (config.item_str["partType"] == "He") {
      if (Estage[3] < (0.00045 * Estage[4] * Estage[4] - 0.42 * Estage[4] + 237) ||
          Estage[3] > (0.00045 * Estage[4] * Estage[4] - 0.42 * Estage[4] + 260)) dropEvent = true;
          }*/

    // DELTA E-E FILTER DOES NOT SIGNIFICANTLY CHANGE CALIBRATION
    /*  if (config.item_str["partType"] == "H") {
        if (Estage[3]<(0.0034*Estage[4] * Estage[4] - 0.672145*Estage[4]
        + 69.) || Estage[3]>(0.00351528*Estage[4] * Estage[4] -
        0.703048*Estage[4] + 79.5048)) dropEvent = true;
        }*/
  }
  else if (Estage[3] > config.item_float["thr3"]) { // Particule stop in stage 3
    if (Estage[2] < cut0 || Estage[1] < cut0 || Estage[0] < cut0 || Estage[3] > cut3) dropEvent = true;
    /*if (config.item_str["partType"] == "He") {
      if (Estage[2] < (0.00085 * Estage[3] * Estage[3] - 0.5664 * Estage[3] + 226) ||
          Estage[2] > (0.00085 * Estage[3] * Estage[3] - 0.5664 * Estage[3] + 254)) dropEvent = true;
          }*/
    // DELTA    E-E FILTER DOES NOT SIGNIFICANTLY CHANGE CALIBRATION
    /*  if (config.item_str["partType"] == "H") {
        if (Estage[2]<(0.0038*Estage[3] * Estage[3] - 0.71*Estage[3] +
        67.92) || Estage[2]>(0.0039*Estage[3] * Estage[3] - 0.74*Estage[4] +
        78.37)) dropEvent = true;
        }*/
  }
  else if (Estage[2] > config.item_float["thr2"]) { // Particule stop in stage 2
    if (Estage[1] < cut0 || Estage[0] < cut0 || Estage[2] > cut3) dropEvent = true;
    /*if (config.item_str["partType"] == "He") {
      if (Estage[1] < (0.00085 * Estage[2] * Estage[2] - 0.5664 * Estage[2] + 220) ||
          Estage[1] > (0.00085 * Estage[2] * Estage[2] - 0.5664 * Estage[2] + 248))
        dropEvent = true; 
	}*/
    // DELTA    E-E FILTER DOES NOT SIGNIFICANTLY CHANGE CALIBRATION
    /*  if (config.item_str["partType"] == "H") {
        if (Estage[1]<(0.0040*Estage[2] * Estage[2] - 0.71*Estage[2] +
        66.17) || Estage[1]>(0.0036*Estage[2] * Estage[2] - 0.73*Estage[2] +
        76.57)) dropEvent = true;
        }*/
  }
  else if (Estage[1] > config.item_float["thr1"]) { // Particule stop in stage 1
    if (Estage[0] < cut1 || Estage[1] > cut3) dropEvent = true;
    /*if (config.item_str["partType"] == "He") {                                                                                                                           
      if (Estage[0] < (0.00085 * Estage[1] * Estage[1] - 0.5664 * Estage[1] + 220) ||
          Estage[0] > (0.00085 * Estage[1] * Estage[1] - 0.5664 * Estage[1] + 246))
        dropEvent = true;
        }*/
  }
  else if (Estage[0]> config.item_float["thr0"]){ // Particule stop in stage 0
    if(Estage[0] > cut3) dropEvent = true; //
  }

  // Maximal energy filter
  if (config.item_str["partType"] == "He") {
    if (Etot > 801.52) dropEvent = true; // Intitial Energy of helium at HIT = 200.36 MeV/u
  }
  else{ // Particle is H
    if(Etot > 200.) dropEvent = true; // Initial Energy of protons at HIT = 200.11 MeV
  }
  return dropEvent;
}


////////////////////////////////////////////////////////////////////
// Tracking cuts
////////////////////////////////////////////////////////////////////
// Call for each raw event after the tracking is completed
bool pCTcut::cutEvt(pCT_Tracking &pCTtracks, TkrHits &pCThits) {
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
	nKeep++;
	good = true;
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
void pCTcut::summary() { // Summary of the processing up to the point of selecting
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

