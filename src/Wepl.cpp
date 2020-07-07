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
#include "Wepl.h"
#include "TTree.h"
Wepl::~Wepl(){}

Wepl::Wepl(TFile* calibFile){
  theConfig = pCTconfig::GetInstance();

  //Load from the calibration file
  TTree* header= (TTree*)calibFile->Get("header");

  // Load the Range-Energy calibration individual stage
  for(int stage =0;stage<5;stage++) calWEPL[stage]= (TGraphErrors*)calibFile->Get(Form("CalibrationCurves/RangeVsEnergy_%d",stage));
  // Combined all stages
  RangeVsEnergy = (TGraphErrors*)calibFile->Get("CalibrationCurves/RangeVsEnergy");
  EnergyVsRange = (TGraphErrors*)calibFile->Get("CalibrationCurves/EnergyVsRange");

  // dEE Filter parameters
  for(int stage =0;stage<5;stage++){ // 1 stage less due to the fact that the 1st-stage can't do dEE
    for(int i =0; i<3; i++){
      header->SetBranchAddress(Form("dEElow_Stage%d_%d", stage,i),&dEElow[stage][i]); // Load the dEE parameters
      header->SetBranchAddress(Form("dEEhigh_Stage%d_%d", stage,i),&dEEhigh[stage][i]);
    }
  }

  
  header->GetEntry(0);
  if (theConfig->item_str["partType"] == "He") {  // sanity cuts
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

  theCuts = pCTcut::GetInstance();
  float EnergyBinWidth = 0.5;
  if (theConfig->item_str["partType"] == "He") EnergyBinWidth = 1.0;
  for(int i =0; i<5; i++){
    dEEhist_root[i] = new TH2D(Form("dE-EStage_%d",i), Form("dE-E spectra for stage %d ", i), // Name-Title
			       nEnrg, 0., nEnrg * EnergyBinWidth, nEnrg, 0., nEnrg * EnergyBinWidth);
  }
}

float Wepl::EtoWEPL(float Estage[5], Int_t &MaxTrans, Int_t &Threshold, Int_t &dEE, float &E_tot) // Return calibrated WEPL from the wedge calibration
{
  float WET;

  if (Estage[4] > theConfig->item_float["thr4"] ) { // // if particle stop in Stage 4
    E_tot = Estage[0]+Estage[1]+Estage[2]+Estage[3]+ Estage[4];
    if(theConfig->item_int["CalibCurve"]==2) WET = RangeVsEnergy->Eval(E_tot);
    else if(theConfig->item_int["CalibCurve"]==1) WET = EnergyVsRange->Eval(E_tot);
    else WET = calWEPL[4]->Eval(Estage[4]);
    
    if (Estage[4] > maxEnergy || Estage[4] < minEnergy) MaxTrans = 0;
    if (Estage[3] < cut3 || Estage[2] < cut2 || Estage[1] < cut1 || Estage[0] < cut0) Threshold = 0;

    if (theConfig->item_int["dEEFilter"] && ! theCuts->dEEFilter(Estage[3], Estage[4], dEElow[4], dEEhigh[4])) dEE = 0;
    else dEEhist_root[4]->Fill(Estage[3], Estage[4]); // Passed dEE
    return WET; // polystyrene to water equivalent // everything passed
  }

  else if (Estage[3] > theConfig->item_float["thr3"]){ // if particle stop in Stage 3
    E_tot = Estage[0]+Estage[1]+Estage[2]+Estage[3];
    if(theConfig->item_int["CalibCurve"]==2) WET = RangeVsEnergy->Eval(E_tot);
    else if(theConfig->item_int["CalibCurve"]==1) WET = EnergyVsRange->Eval(E_tot);
    else WET = calWEPL[3]->Eval(Estage[3]);

    if (Estage[3] > maxEnergy || Estage[3] < minEnergy) MaxTrans = 0;
    if (Estage[2] < cut2 || Estage[1] < cut1 || Estage[0] < cut0) Threshold = 0;
    if (theConfig->item_int["dEEFilter"] && !theCuts->dEEFilter(Estage[2], Estage[3], dEElow[3], dEEhigh[3])) dEE = 0;
    else dEEhist_root[3]->Fill(Estage[2], Estage[3]);
    return WET; // polystyrene to water equivalent
  }
  else if (Estage[2] > theConfig->item_float["thr2"]) { // if particle stop in Stage 2
    E_tot = Estage[0]+Estage[1]+Estage[2];
    if(theConfig->item_int["CalibCurve"]==2) WET = RangeVsEnergy->Eval(E_tot);
    else if(theConfig->item_int["CalibCurve"]==1) WET = EnergyVsRange->Eval(E_tot);
    else WET = calWEPL[2]->Eval(Estage[2]);

    if (Estage[2] > maxEnergy || Estage[2] < minEnergy) MaxTrans = 0;
    if (Estage[1] < cut1 || Estage[0] < cut0) Threshold = 0;
    if (theConfig->item_int["dEEFilter"] && !theCuts->dEEFilter(Estage[1], Estage[2], dEElow[2], dEEhigh[2])) dEE = 0;
    else dEEhist_root[2]->Fill(Estage[1], Estage[2]);
    return WET; // polystyrene to water equivalent
  }

  else if (Estage[1] > theConfig->item_float["thr1"]) {   // if particle stop in Stage 1
    E_tot = Estage[0]+Estage[1];
    if(theConfig->item_int["CalibCurve"]==2) WET = RangeVsEnergy->Eval(E_tot);
    else if(theConfig->item_int["CalibCurve"]==1) WET = EnergyVsRange->Eval(E_tot);
    else WET = calWEPL[1]->Eval(Estage[1]);

    if (Estage[1] > maxEnergy || Estage[1] < minEnergy) MaxTrans = 0;
    if (Estage[0] < cut0) Threshold = 0;
    if (theConfig->item_int["dEEFilter"] && !theCuts->dEEFilter(Estage[0], Estage[1], dEElow[1], dEEhigh[1])) dEE = 0;
    else dEEhist_root[1]->Fill(Estage[0], Estage[1]);

    return WET; // polystyrene to water equivalent
  }
  else if (Estage[0] > theConfig->item_float["thr0"]) {   // if particle stop in Stage 0
    E_tot = Estage[0];
    if(theConfig->item_int["CalibCurve"]==2) WET = RangeVsEnergy->Eval(E_tot);
    if(theConfig->item_int["CalibCurve"]==1) WET = EnergyVsRange->Eval(E_tot);
    else WET = calWEPL[0]->Eval(Estage[0]);

    if( WET < 0 ) return -1000;  // Unit Test
    if (Estage[0] > maxEnergy || Estage[0] < minEnergy) MaxTrans = 0;
    return WET; // polystyrene to water equivalent
  }
  else return 2000;    // No Energy found in any stage above the thresold, returning big WEPL
}
void Wepl::WriteHist(TFile* projectionROOT){

  projectionROOT->cd();
  projectionROOT->mkdir("dEE");
  projectionROOT->cd("dEE");
  for(int i =1;i<5;i++) dEEhist_root[i]->Write("",TObject::kOverwrite);
  projectionROOT->cd();
  projectionROOT->mkdir("calWEPL");
  projectionROOT->cd("calWEPL");
  RangeVsEnergy->Write("",TObject::kOverwrite);
  EnergyVsRange->Write("",TObject::kOverwrite);
  for(int i =0;i<5;i++) calWEPL[i]->Write("",TObject::kOverwrite);

}
