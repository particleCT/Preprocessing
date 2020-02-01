/* -----------------------------------------------------------------------------
Wepl.h        Class merged from Wepl1 and Wepl2 [C.E. Ordonez, Aug 2016]
Methods:
Wepl    Constructor, reads calibration data from file previously
        prepared by one of the following:
        PreProcessingWcalib: step calibration, scans before Aug 2016
        PreProcessingWcalibW: wedge calibration, scans starting Aug 2016
        The filenames are no longer restricted to Wcalib.txt and
        WcalibW.txt; nor the calibration files themselves need to be in
        the same output directory for the projection files.
float EtoWEPL    Converts energy (in Mev) deposited in 5 stages to WEPL (in mm)
----------------------------------------------------------------------------- */
#ifndef _WEPL_H_
#define _WEPL_H_
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TH2D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "pCTconfig.h"
using namespace std;
#define nEnrg 340
class Wepl {
 public:
  float cut0, cut1, cut2, cut3, cut4, cut5;
  float dEElow[4][3], dEEhigh[4][3];

  // Common parameters
  float thr[5];//0, thr1, thr2, thr3, thr4; // Stage thresholds to define when a stage has a real signal
  float RSP;                          // Known phantom RSP
  TFile* projectionROOT;
  pCTconfig config;
  TGraphErrors* calWEPL[5];
  TH2D* dEEhist_root[5];
  // Explicit constructor
  Wepl(pCTconfig, TFile*, TFile*);
       
  // Explicit destructor
  ~Wepl();
  // Set energy thresholds in stages
  void SetEthresholds1(float, float, float, float, float);

  // Convert energy in stages to WEPL (mm)
  float EtoWEPL(float Estage[5], Int_t &, Int_t &, Int_t &);
  // some util for splitting key values for the dEE
  vector<string> split(string str, char delimiter);
};

inline Wepl::~Wepl(){
  projectionROOT->mkdir("dEE");
  projectionROOT->cd("dEE");
  for(int i =1;i<5;i++) dEEhist_root[i]->Write("",TObject::kOverwrite);
  for(int i =1;i<5;i++) dEEhist_root[0]->Add(dEEhist_root[i]);
  dEEhist_root[0]->Write("dEEhist_tot",TObject::kOverwrite);
  
  projectionROOT->cd();
  projectionROOT->mkdir("calWEPL");
  projectionROOT->cd("calWEPL");
  for(int i =0;i<5;i++) calWEPL[i]->Write("",TObject::kOverwrite);

}  

inline Wepl::Wepl(pCTconfig cfg, TFile* calibFile, TFile* outputFile): config(cfg)
{
  this->projectionROOT = outputFile;
  TTree* header= (TTree*)calibFile->Get("header");
  float scale = 1.0;
  if (config.item_str["partType"] == "He") {  // sanity cuts
    cut0 = 87. ; cut1 = 95.;  cut2 = 100.;
    cut3 = 120.; cut4 = 300.; cut5 = 270.;
  } else {
    cut0 = 19.; cut1 = 21.; cut2 = 25.;
    cut3 = 35.; cut4 = 75.; cut5 = 77.;
  }
  float EnergyBinWidth = 0.25;
  if (config.item_str["partType"] == "He") EnergyBinWidth = 1.0;
  if (config.item_str["partType"] == "He") {
    // lower limits
    dEElow[0][0] = 0.00085; dEElow[0][1] = -0.5664; dEElow[0][2] = 220;
    dEElow[1][0] = 0.00085; dEElow[1][1] = -0.5664; dEElow[1][2] = 222;
    dEElow[2][0] = 0.00085; dEElow[2][1] = -0.5664; dEElow[2][2] = 226;
    dEElow[3][0] = 0.00045; dEElow[3][1] = -0.42;   dEElow[3][2] = 237;
    
    // upper limits
    dEEhigh[0][3] = 0.00085; dEEhigh[0][4] = -0.5664; dEEhigh[0][5] = 246;
    dEEhigh[1][3] = 0.00085; dEEhigh[1][4] = -0.5664; dEEhigh[1][5] = 248;
    dEEhigh[2][3] = 0.00085; dEEhigh[2][4] = -0.5664; dEEhigh[2][5] = 254;
    dEEhigh[3][3] = 0.00045; dEEhigh[3][4] = -0.42;   dEEhigh[3][5] = 260;

  } else {
    // lower limits
    dEElow[0][0] = 0.0040656; dEElow[0][1] = -0.719491; dEElow[0][2] = 65.9463;
    dEElow[1][0] = 0.00378842;dEElow[1][1] = -0.701085; dEElow[1][2] = 66.1436;
    dEElow[2][0] = 0.00380467;dEElow[2][1] = -0.71088;  dEElow[2][2] = 67.9228;
    dEElow[3][0] = 0.0032026; dEElow[3][1] = -0.663234; dEElow[3][2] = 69.0353;
    // upper limits
    dEEhigh[0][3] = 0.00365163;dEEhigh[0][4] = -0.735325; dEEhigh[0][5] = 76.3504;
    dEEhigh[1][3] = 0.00384016;dEEhigh[1][4] = -0.732287; dEEhigh[1][5] = 76.5896;
    dEEhigh[2][3] = 0.00385641;dEEhigh[2][4] = -0.742083; dEEhigh[2][5] = 78.3689;
    dEEhigh[3][3] = 0.00372006;dEEhigh[3][4] = -0.709805; dEEhigh[3][5] = 79.5232;
  }

  for(int i =0; i<5; i++){
    dEEhist_root[i] = new TH2D(Form("dE-EStage_%d_bricks",i), Form("dE-E spectra for stage %d ", i), // Name-Title
			       nEnrg, 0., nEnrg * EnergyBinWidth, nEnrg, 0., nEnrg * EnergyBinWidth);} // X-Y binning 
      
  for(int stage =0;stage<5;stage++){
    calWEPL[stage]= (TGraphErrors*)calibFile->Get(Form("RangeVsEnergy_Uncorr/RangeVsEnergy_Uncorr_%d",stage)); // Load the Range-Energy calibration
    header->SetBranchAddress(Form("thr_%d", stage),&thr[stage]); // Load the threshold
  }

  for(int stage =0;stage<4;stage++){ // 1 stage less due to the fact that the 1st-stage can't do dEE
    for(int i =0; i<3; i++){
      header->SetBranchAddress(Form("dEElow_Stage%d_%d", stage,i),&dEElow[stage][i]); // Load the dEE parameters
      header->SetBranchAddress(Form("dEEhigh_Stage%d_%d", stage,i),&dEEhigh[stage][i]);      
    }
  }
  header->GetEntry(0);
  RSP = 1.030; // 1.0300 is for real calibration phantoms
  cout << "The RSP of the WEPL calibration phantom is assumed to be " << RSP << endl;
}

inline float Wepl::EtoWEPL(float Estage[5], Int_t &MaxTrans, Int_t &Threshold, Int_t &dEE) // Return calibrated WEPL from the wedge calibration
{
  double WET;
  // dE-E parameterization and energy check to cut out fragments for helium if particle stop in Stage 4
  if (Estage[4] > thr[4]) { // 1 MeV
    dEEhist_root[4]->Fill(Estage[3], Estage[4]);
    WET = calWEPL[4]->Eval(Estage[4])*RSP;
    if( WET < 0 ) return -1000;  // Unit Test
    if (Estage[4] > cut4) MaxTrans = 0; // 1000 // Max Trans Filter
    if (Estage[3] < cut3 || Estage[2] < cut2 || Estage[1] < cut1 || Estage[0] < cut0) Threshold = 0; 
    if (config.item_int["dEEFilter"]) { // dEEFilter
      if (Estage[3] < (dEElow[3][0]  * pow(Estage[4],2) + dEElow[3][1] * Estage[4] + dEElow[3][2]) ||
          Estage[3] > (dEEhigh[3][3] * pow(Estage[4],2) + dEEhigh[3][4] * Estage[4] + dEEhigh[3][5])) dEE=0;
    }
    return WET; // polystyrene to water equivalent // everything passed
  }
  // if particle stop in Stage 3
  else if (Estage[3] > thr[3]) {
    dEEhist_root[3]->Fill(Estage[2], Estage[3]);
    WET = calWEPL[3]->Eval(Estage[3])*RSP;
    if( WET < 0 ) return -1000; // Unit Test
    if (Estage[3] > cut5) MaxTrans = 0;
    if (Estage[2] < cut2 || Estage[1] < cut1 || Estage[0] < cut0) Threshold = 0; 
    if (config.item_int["dEEFilter"]) { // dEEFilter
      if (Estage[2] < (dEElow[2][0]  * pow(Estage[3],2) + dEElow[2][1] * Estage[3] + dEElow[2][2]) ||
          Estage[2] > (dEEhigh[2][3] * pow(Estage[3],2) + dEEhigh[2][4] * Estage[3] + dEEhigh[2][5])) dEE = 0;        
    }
    return WET; // polystyrene to water equivalent
  }
  // if particle stop in Stage 2
  else if (Estage[2] > thr[2]) { // 1 MeV
    dEEhist_root[2]->Fill(Estage[1], Estage[2]);
    WET = calWEPL[2]->Eval(Estage[2])*RSP;
    if (Estage[2] > cut5) MaxTrans = 0;
    if (Estage[1] < cut1 || Estage[0] < cut0) Threshold = 0;
    if (config.item_int["dEEFilter"]) { // dEEFilter
      if (Estage[1] < (dEElow[1][0]  * pow(Estage[2],2) + dEElow[1][1] * Estage[2] + dEElow[1][2]) ||
          Estage[1] > (dEEhigh[1][3] * pow(Estage[2],2) + dEEhigh[1][4] * Estage[2] + dEEhigh[1][5])) dEE = 0;
    }
    return WET; // polystyrene to water equivalent
  }
  // if particle stop in Stage 1
  else if (Estage[1] > thr[1]) {
    dEEhist_root[1]->Fill(Estage[0], Estage[1]);
    WET = calWEPL[1]->Eval(Estage[1])*RSP;
    if( WET < 0 ) return -1000;  // Unit Test
    if (Estage[1] > cut5) MaxTrans = 0;
    if (Estage[0] < cut0) Threshold = 0;
    if (config.item_int["dEEFilter"]) { // dEEFilter
      if (Estage[0] < (dEElow[0][0]  * pow(Estage[1],2) + dEElow[0][1] * Estage[1] + dEElow[0][2]) ||
          Estage[0] > (dEEhigh[0][3] * pow(Estage[1],2) + dEEhigh[0][4] * Estage[1] + dEEhigh[0][5])) dEE = 0;
        
    }
    return WET; // polystyrene to water equivalent
  }
  // if particle stop in Stage 0
  else if (Estage[0] > thr[0]) {
    WET = calWEPL[0]->Eval(Estage[0])*RSP;
    if(WET>1000) cout<<"0 "<<WET<<endl;
    if( WET < 0 ) return -1000;  // Unit Test
    if (Estage[0] > cut5) MaxTrans = 0;
    return WET; // polystyrene to water equivalent
  }
  else return 2000;    // No Energy found in any stage above the thresold, returning big WEPL
}
inline void Wepl::SetEthresholds1(float t0, float t1, float t2, float t3, float t4) {
  thr[0] = t0; thr[1] = t1; thr[2] = t2;
  thr[3] = t3; thr[4] = t4;
  cout << "Wepl::SetEthresholds1: WEPL energy thresholds are set to " << thr[0] << " " << thr[1]<< " " << thr[2] << " " << thr[3] << " " << thr[4]<< " MeV" << endl;
}
#endif 
