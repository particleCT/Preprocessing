#include "pCTcut.h"

pCTcut::pCTcut()
{ // Class constructor called prior to the event loop

  theConfig = pCTconfig::GetInstance();
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

//////////////////////////////////////////////////////////////////////
// dEE Parameters
//////////////////////////////////////////////////////////////////////
void pCTcut::dEEFilterParameters(TH2D* dEEhist, float dEElow[3], float dEEhigh[3], int stage){
  float EnergyBinWidth;
  if (theConfig->item_str["partType"] == "H") EnergyBinWidth = 0.5;
  else EnergyBinWidth = 1.0;

  /// Lennart Volz, November 2018 dE-E parameter evaluation:
  /*  int Estep[5][3];
  if(theConfig->item_str["partType"] == "H"){
    Estep[0][0] = 40; Estep[0][1] = 0;   Estep[0][2] = 0;
    Estep[1][0] = 40; Estep[1][1] = 120; Estep[1][2] = 200;
    Estep[2][0] = 40; Estep[2][1] = 100; Estep[2][2] = 160;
    Estep[3][0] = 40; Estep[3][1] = 80;  Estep[3][2] = 120;
    Estep[4][0] = 20; Estep[4][1] = 50;  Estep[4][2] = 70;
  }
  else{
    Estep[0][0] = 40; Estep[0][1] = 0;   Estep[0][2] = 0;
    Estep[1][0] = 40; Estep[1][1] = 120; Estep[1][2] = 230;
    Estep[2][0] = 40; Estep[2][1] = 100; Estep[2][2] = 230;
    Estep[3][0] = 40; Estep[3][1] = 100; Estep[3][2] = 230;
    Estep[4][0] = 40; Estep[4][1] = 100; Estep[4][2] = 180;
    }*/
  int Estep[3] = { 240, 120, 12 }; 
  int bin_cut  = 0;
  float E[3] = {Estep[0]*EnergyBinWidth, Estep[1]*EnergyBinWidth, Estep[2]*EnergyBinWidth};
  float xlow[3], xhigh[3], xmax, max;
  for (int j = 0; j < 3; j++) {
    TH1D* dEESlice = dEEhist->ProjectionX(Form("ProjX_Stage%d_Bin%d",stage,Estep[j]),Estep[j] ,Estep[j]+1);
    // Cheesy fix to remove the low energy spike
    if(theConfig->item_str["partType"] == "H") bin_cut = dEESlice->FindBin(40);
    else bin_cut = dEESlice->FindBin(60);
    dEESlice->GetXaxis()->SetRange(bin_cut,dEESlice->GetNbinsX());
    max   =  dEESlice->GetMaximum();
    xmax  =  dEESlice->GetBinCenter(dEESlice->GetMaximumBin());
    TF1 *f1 = new TF1("f1", "gaus", xmax-10, xmax+10);
    f1->SetParameter(0, max);
    f1->SetParameter(1, xmax);
    Int_t fitStatus = dEESlice->Fit(f1,"QR"); // Q means quiet, R means ue the same range as the histogram
    if(fitStatus==0) // Everything passed
      {
	xmax      = f1->GetParameter(1);// mean
	xhigh[j]  = xmax + f1->GetParameter(2); //sigma
	xlow[j]   = xmax - f1->GetParameter(2); //sigma
      }
    else{
      xmax        =  dEESlice->GetBinCenter(dEESlice->GetMaximumBin());
      xhigh[j]    =  dEESlice->GetBinCenter(dEESlice->FindLastBinAbove(dEESlice->GetMaximum()/2));
      xlow[j]     =  xmax - (xhigh[j] -xmax);
    }
    xlow[j]     = xlow[j] - (xhigh[j] - xlow[j]) * 0.7848;
    xhigh[j]    =   xhigh[j] + (xhigh[j] - xlow[j]) * 0.7848;
  }
  dEElow[0] = (E[0] * (xlow[2] - xlow[1]) + E[1] * (xlow[0] - xlow[2]) + E[2] * (xlow[1] - xlow[0])) /
              ((E[0] - E[1]) * (E[0] - E[2]) * (E[1] - E[2]));
  dEElow[1]  = (xlow[1] - xlow[0]) / (E[1] - E[0]) - dEElow[0] * (E[0] + E[1]);
  dEElow[2]  = xlow[0] - dEElow[0] * E[0] * E[0] - dEElow[1] * E[0];

  dEEhigh[0] = (E[0] * (xhigh[2] - xhigh[1]) + E[1] * (xhigh[0] - xhigh[2]) + E[2] * (xhigh[1] - xhigh[0])) /
               ((E[0] - E[1]) * (E[0] - E[2]) * (E[1] - E[2]));
  dEEhigh[1] = (xhigh[1] - xhigh[0]) / (E[1] - E[0]) - dEEhigh[0] * (E[0] + E[1]);
  dEEhigh[2] = xhigh[0] - dEEhigh[0] * E[0] * E[0] - dEEhigh[1] * E[0];
}
////////////////////////////////////////////////////////////////////
// dEE cuts
////////////////////////////////////////////////////////////////////
bool pCTcut::dEEFilter(float Elow, float Ehigh, float dEElow[3], float dEEhigh[3]){
  if (Elow < (dEElow[0]  * pow(Ehigh,2) + dEElow[1] * Ehigh + dEElow[2]) ||
      Elow > (dEEhigh[0] * pow(Ehigh,2) + dEEhigh[1] * Ehigh + dEEhigh[2])){
    return false;}
  else return true;
}
////////////////////////////////////////////////////////////////////
// Energy cuts
////////////////////////////////////////////////////////////////////
bool pCTcut::EnrgCut(float Estage[5], float Etot, float cut0, float cut1, float cut3) {
  bool dropEvent = false;
  if (Estage[4] >  theConfig->item_float["thr4"]) { // Particule stop in stage 4
    if (Estage[3] < cut0 || Estage[2] < cut0 || Estage[1] < cut0 || Estage[0] < cut0 || Estage[4] > cut3) dropEvent = true;
  }
  else if (Estage[3] > theConfig->item_float["thr3"]) { // Particule stop in stage 3
    if (Estage[2] < cut0 || Estage[1] < cut0 || Estage[0] < cut0 || Estage[3] > cut3) dropEvent = true;
  }
  else if (Estage[2] > theConfig->item_float["thr2"]) { // Particule stop in stage 2
    if (Estage[1] < cut0 || Estage[0] < cut0 || Estage[2] > cut3) dropEvent = true;
  }
  else if (Estage[1] > theConfig->item_float["thr1"]) { // Particule stop in stage 1
    if (Estage[0] < cut1 || Estage[1] > cut3) dropEvent = true;
  }
  else if (Estage[0]> theConfig->item_float["thr0"]){ // Particule stop in stage 0
    if(Estage[0] > cut3) dropEvent = true; //
  }

  // Maximal energy filter
  if (theConfig->item_str["partType"] == "He") {
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

