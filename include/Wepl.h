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
#include "pCTcut.h"
#include "pCTconfig.h"
using namespace std;
#define nEnrg 340
class Wepl {
 public:
  float cut0, cut1, cut2, cut3, cut4, maxEnergy, minEnergy;
  float dEElow[4][3], dEEhigh[4][3];

  // Common parameters
  float thr[5];//0, thr1, thr2, thr3, thr4; // Stage thresholds to define when a stage has a real signal
  float RSP;                          // Known phantom RSP
  TFile* projectionROOT;
  pCTconfig config;
  TGraphErrors* calWEPL[5];
  TH2D* dEEhist_root[5];
  TGraphErrors* Test;
  TGraphErrors* Test2;
  pCTcut *theCuts;
  // Explicit constructor
  Wepl(pCTconfig, TFile*, TFile*);
       
  // Explicit destructor
  ~Wepl();
  // Set energy thresholds in stages
  void SetEthresholds(float, float, float, float, float);

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
  Test->Write("",TObject::kOverwrite);
  Test2->Write("",TObject::kOverwrite);
  for(int i =0;i<5;i++) calWEPL[i]->Write("",TObject::kOverwrite);
  
}  

inline Wepl::Wepl(pCTconfig cfg, TFile* calibFile, TFile* outputFile): config(cfg)
{
  this->projectionROOT = outputFile;
  TTree* header= (TTree*)calibFile->Get("header");
  float scale = 1.0;
  if (config.item_str["partType"] == "He") {  // sanity cuts
    cut0 = 87. ; cut1 = 95.;  cut2 = 100.;
    cut3 = 120.; cut4 = 300.;
    maxEnergy = 270;
    minEnergy = 1;
  } else {
    cut0 = 19.; cut1 = 21.; cut2 = 25.;
    cut3 = 35.; cut4 = 75.;
    maxEnergy = 77.;
    minEnergy = 1;
  }

  theCuts = new pCTcut(config);
  
  float EnergyBinWidth = 0.5;
  if (config.item_str["partType"] == "He") EnergyBinWidth = 1.0;
  for(int i =0; i<5; i++){
    dEEhist_root[i] = new TH2D(Form("dE-EStage_%d_bricks",i), Form("dE-E spectra for stage %d ", i), // Name-Title
			       nEnrg, 0., nEnrg * EnergyBinWidth, nEnrg, 0., nEnrg * EnergyBinWidth);
  } // X-Y binning 
      
  for(int stage =0;stage<5;stage++){
    calWEPL[stage]= (TGraphErrors*)calibFile->Get(Form("RangeVsEnergy/RangeVsEnergy_%d",stage)); // Load the Range-Energy calibration
    header->SetBranchAddress(Form("thr_%d", stage),&thr[stage]); // Load the threshold
  }
  //Test  = (TGraphErrors*)calibFile->Get("NewRange_YX");
  Test2 = (TGraphErrors*)calibFile->Get("NewRange_XY");
  Test = (TGraphErrors*)calibFile->Get("RangeVsEnergy/RangeVsEnergy_YX");
  for(int stage =0;stage<5;stage++){ // 1 stage less due to the fact that the 1st-stage can't do dEE
    for(int i =0; i<3; i++){
      header->SetBranchAddress(Form("dEElow_Stage%d_%d", stage,i),&dEElow[stage][i]); // Load the dEE parameters
      header->SetBranchAddress(Form("dEEhigh_Stage%d_%d", stage,i),&dEEhigh[stage][i]);

    }
  }
  header->GetEntry(0);
}

inline float Wepl::EtoWEPL(float Estage[5], Int_t &MaxTrans, Int_t &Threshold, Int_t &dEE) // Return calibrated WEPL from the wedge calibration
{
  double WET;
  double E_tot;
  E_tot = Estage[0] + Estage[1] + Estage[2] + Estage[3] + Estage[4];
  // dE-E parameterization and energy check to cut out fragments for helium if particle stop in Stage 4
  if (Estage[4] > thr[4]) { // 1 MeV
    dEEhist_root[4]->Fill(Estage[3], Estage[4]);
    WET = Test->Eval(E_tot);
    //WET = calWEPL[4]->Eval(Estage[4]);    
    if( WET < 0 ) return -1000;  // Unit Test
    if (Estage[4] > maxEnergy || Estage[4] < minEnergy) MaxTrans = 0;
    if (Estage[3] < cut3 || Estage[2] < cut2 || Estage[1] < cut1 || Estage[0] < cut0) Threshold = 0;
    if (config.item_int["dEEFilter"] && !theCuts->dEEFilter(Estage[3], Estage[4], dEElow[4], dEEhigh[4])) dEE = 0;
    return WET; // polystyrene to water equivalent // everything passed
  }
  // if particle stop in Stage 3
  else if (Estage[3] > thr[3]) {

    WET = Test->Eval(E_tot);
    dEEhist_root[3]->Fill(Estage[2], Estage[3]);
    //WET = calWEPL[3]->Eval(Estage[3]);
    if( WET < 0 ) return -1000; // Unit Test
    if (Estage[3] > maxEnergy || Estage[3] < minEnergy) MaxTrans = 0;
    if (Estage[2] < cut2 || Estage[1] < cut1 || Estage[0] < cut0) Threshold = 0; 
    if (config.item_int["dEEFilter"] && !theCuts->dEEFilter(Estage[2], Estage[3], dEElow[3], dEEhigh[3])) dEE = 0;
    return WET; // polystyrene to water equivalent
  }
  // if particle stop in Stage 2
  else if (Estage[2] > thr[2]) { // 1 MeV
    dEEhist_root[2]->Fill(Estage[1], Estage[2]);
    WET = Test->Eval(E_tot);
    //WET = calWEPL[2]->Eval(Estage[2]);
    if (Estage[2] > maxEnergy || Estage[2] < minEnergy) MaxTrans = 0;
    if (Estage[1] < cut1 || Estage[0] < cut0) Threshold = 0;
    if (config.item_int["dEEFilter"] && !theCuts->dEEFilter(Estage[1], Estage[2], dEElow[2], dEEhigh[2])) dEE = 0;
    return WET; // polystyrene to water equivalent
  }
  // if particle stop in Stage 1
  else if (Estage[1] > thr[1]) {
    dEEhist_root[1]->Fill(Estage[0], Estage[1]);
    WET = Test->Eval(E_tot);
    // WET = calWEPL[1]->Eval(Estage[1]);
    if( WET < 0 ) return -1000;  // Unit Test
    if (Estage[1] > maxEnergy || Estage[1] < minEnergy) MaxTrans = 0;
    if (Estage[0] < cut0) Threshold = 0;
    if (config.item_int["dEEFilter"] && !theCuts->dEEFilter(Estage[0], Estage[1], dEElow[1], dEEhigh[1])) dEE = 0;
    return WET; // polystyrene to water equivalent
  }
  // if particle stop in Stage 0
  else if (Estage[0] > thr[0]) {
    WET = Test->Eval(E_tot);
    //WET = calWEPL[0]->Eval(Estage[0]);
    if( WET < 0 ) return -1000;  // Unit Test
    if (Estage[0] > maxEnergy || Estage[0] < minEnergy) MaxTrans = 0;
    return WET; // polystyrene to water equivalent
  }
  else return 2000;    // No Energy found in any stage above the thresold, returning big WEPL
}
inline void Wepl::SetEthresholds(float t0, float t1, float t2, float t3, float t4) {
  thr[0] = t0; thr[1] = t1; thr[2] = t2; thr[3] = t3; thr[4] = t4;
  cout << "Wepl::SetEthresholds1: WEPL energy thresholds are set to " << thr[0] << " " << thr[1]<< " " << thr[2] << " " << thr[3] << " " << thr[4]<< " MeV" << endl;
}
#endif 
