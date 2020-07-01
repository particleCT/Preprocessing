#ifndef TVcorrection_h
#define TVcorrection_h
// To read the calibration constants from a file for the T,V energy detector corrections
// R.P. Johnson  5/22/2016 to encapculate the TV calibration code of Vladimir Bashkirov
// R.P. Johnson  10/2/2016 added a binlinear interpolation of the TV map
// R.P. Johnson  10/21/2016 added optional checks on date and run ranges

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "TH2D.h"
#include "Util.h"
#include "TTree.h"
#include "TFile.h"
using namespace std;
#define nStage 5
#define nPixX 34
#define nPixY 10  // under and overflow
#define nPix (nPixX+2)*(nPixY+2)

class TVcorrection {
public:
  float ped[5] = {0};   // Calibrated pedestals (the program will generally use values  generated on the fly, however)
  float Eempt[6] = {0}; // ADC pedestals and energy depositions from Empty run  
  TH2D* TVcorrHist[5]; // TVcorr histogram for each stage

  // Constructor -- read from the TV calibration file
  TVcorrection(TFile* calibFile, int calib) { 

    if(calib){// Calibration we initialize the data be empty
      for(int stage =0; stage<nStage; stage++) TVcorrHist[stage] = new TH2D(Form("TVcorrMap_%d", stage), "", nPixX, -190, 190, nPixY, -50, 50);
    }
    else{ // Preprocessing we fill from the calibration file
      TTree* header= (TTree*)calibFile->Get("header");
      for(int stage =0; stage<nStage; stage++){
	header->SetBranchAddress(Form("Est_%d", stage),&Eempt[stage]);
	header->SetBranchAddress(Form("peds_%d",stage),&ped[stage]);
	TVcorrHist[stage] = (TH2D*)calibFile->Get(Form("TVcorrProfile/TVcorrMap_%d",stage));
      }
      header->GetEntry(0);
    }
  }
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  float corrFactor(int stage, float T, float V, bool &inBounds) {

    inBounds = true;                    
    int tPix = TVcorrHist[stage]->GetXaxis()->FindBin(T);
    int vPix = TVcorrHist[stage]->GetYaxis()->FindBin(V);

    if (tPix < 0 || tPix > TVcorrHist[stage]->GetNbinsX()-1) {
      inBounds = false;
    }
    if (vPix < 0 || vPix > TVcorrHist[stage]->GetNbinsY()-1) {
      vPix = 0;
      inBounds = false;
    }
    if (inBounds){
      Int_t binx = TVcorrHist[stage]->GetXaxis()->FindFixBin(T);
      Int_t biny = TVcorrHist[stage]->GetYaxis()->FindFixBin(V);
      for(int i =binx-1; i<binx+2; i++){
	for(int j =biny-1; j<biny+2; j++){
	  if(TVcorrHist[stage]->GetBinContent(binx,biny) >0.9){ inBounds = false;
	    return TVcorrHist[stage]->GetBinContent( TVcorrHist[stage]->FindBin(T,V));
	  }
	  else return TVcorrHist[stage]->Interpolate(T,V); // bilinear interpolation
	}      
      }
    }
    else return 0.;
  } 
}; // end of TVcorrection class
#endif
