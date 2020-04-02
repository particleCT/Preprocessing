// Encapsulate all of Vladimir Bashkirov pedestal and gain on-the-fly
// calibration
// R.P. Johnson   5/22/2016

#include "pedGainCalib.h"
#include "TFile.h"
#include "TROOT.h"
#include "TFitResultPtr.h"
pedGainCalib::~pedGainCalib()
{
}

pedGainCalib::pedGainCalib(TFile* root, int pedMin[5], float oldPed[5], float t1, float t2, float t3, float t4):RootFile(root)
{  
  theConfig = pCTconfig::GetInstance();
  // Two ranges in t occupied by unobstructed (empty) protons
  emtrng1 = t1; // negative t side
  emtrng2 = t2;
  emtrng3 = t3; // positive t side
  emtrng4 = t4;
  for (int stage = 0; stage < 5; stage++) {
    Ped[stage] = oldPed[stage];
    GainFac[stage] = 1.0;

    cout << "pedGainCalib setting the default pedestal for stage " << stage << " to " << Ped[stage] << endl;
  }
  // Define histograms for pedestal calculation
  for (int i =0; i<5; i++) hPed[i] = new TH1D(Form("PedestalStage_%i",i), Form("Pedestal region for stage %i",i), 400, pedMin[i], pedMin[i] +400*5);
  for (int i =0; i<5; i++) hTotFil[i] = new TH1D(Form("FullADCStage_Filtered_%i",i), Form("Full ADC for stage %i",i), 400, pedMin[i], pedMin[i] +400*5);
  for (int i =0; i<5; i++) hTotUnFil[i] = new TH1D(Form("FullADCStage_Unfiltered_%i",i), Form("Full ADC for stage %i",i), 400, pedMin[i], pedMin[i] +400*50);


  // Define histograms for pedestal calculation
  for (int i =0; i<5; i++) hPed_In[i] = new TH1D(Form("PedestalStage_%i_In",i), Form("Pedestal region for stage %i",i), 400, pedMin[i], pedMin[i] +400*5);
  for (int i =0; i<5; i++) hPed_Out[i] = new TH1D(Form("PedestalStage_%i_Out",i), Form("Pedestal region for stage %i",i), 400, pedMin[i], pedMin[i] +400*5);

  // Define histograms for gain calibration
  if (theConfig->item_str["partType"] == "H") {
    for (int i =0; i<5; i++) hEnrg[i] = new TH1D(Form("EnergyDistribution_%i",i), Form("Energy Distribution for stage %i",i), 800, 15, 15 + 800*0.175);
    hEnrgTot = new TH1D("SumStageEnergies", "Sum of stage energies", 800, 0, 0 + 0.3*800);
  }
  else {
    for (int i =0; i<5; i++) hEnrg[i] = new TH1D(Form("EnergyDistribution_%i",i), Form("Energy Distribution for stage %i",i), 800, 15, 15 + 800*0.9);
    hEnrgTot = new TH1D("SumStageEnergies", "Sum of stage energies", 800, 0, 0 + 1.2*800);
  }
  // Profile plot to make sure that the phantom does not extend into the regions used for gain calibration
  hProfT = new TProfile2D("Stage0EnergyProfile", "Stage 0 energy profile", 100, -150, -150 + 100*3.0, 100, -50., -50 + 1.0*100);
  hTedet = new TH1D("T_Ions_GainRecalib", "T of ions used for gain recalibration", 100, -150, -150 + 100*3.0);
  RootFile->mkdir("Pedestals");
}

void pedGainCalib::FillPeds(pCTraw &rawEvt) { // Called for each raw event read in from the input data file
  // Accumulate histograms for pedestal analysis  
  hPed[0]->Fill(rawEvt.enrg_fpga[0].pulse_sum[0]);
  hPed[1]->Fill(rawEvt.enrg_fpga[0].pulse_sum[1]);
  hPed[2]->Fill(rawEvt.enrg_fpga[0].pulse_sum[2]);
  hPed[3]->Fill(rawEvt.enrg_fpga[1].pulse_sum[0]);
  hPed[4]->Fill(rawEvt.enrg_fpga[1].pulse_sum[1]);
  hTotUnFil[0]->Fill(rawEvt.enrg_fpga[0].pulse_sum[0]);
  hTotUnFil[1]->Fill(rawEvt.enrg_fpga[0].pulse_sum[1]);
  hTotUnFil[2]->Fill(rawEvt.enrg_fpga[0].pulse_sum[2]);
  hTotUnFil[3]->Fill(rawEvt.enrg_fpga[1].pulse_sum[0]);
  hTotUnFil[4]->Fill(rawEvt.enrg_fpga[1].pulse_sum[1]);
}
  
void pedGainCalib::ResetHist(){
  for(int i =0; i<5; i++){
    hTotFil[i]->Reset();
    hTotUnFil[i]->Reset();
    hPed[i]->Reset();
    hPed_Out[i]->Reset();
    hPed_In[i]->Reset();
    hEnrg[i]->Reset();
  }
  hEnrgTot->Reset();
  hTedet->Reset();
  hProfT->Reset();
}

void pedGainCalib::FillADC(int ADC[5]){ for(int i =0; i<5; i++) hTotFil[i]->Fill(ADC[i]);}

void pedGainCalib::WriteHist(){ // For analysis
  inFileName_s = theConfig->item_str["inputFileName"];
  inFileName_s = inFileName_s.substr(inFileName_s.find_last_of("\\/") +1,  inFileName_s.size());
  RootFile->mkdir(Form("Pedestals/%s",inFileName_s.c_str()));
  RootFile->cd(Form("Pedestals/%s",inFileName_s.c_str()));
  for(int i =0; i<5; i++) {
    hTotFil[i]->Write("", TObject::kOverwrite);
    hTotUnFil[i]->Write("", TObject::kOverwrite);
    hPed[i]->Write("", TObject::kOverwrite);
    hPed_Out[i]->Write("", TObject::kOverwrite);
    hPed_In[i]->Write("", TObject::kOverwrite);
    hEnrg[i]->Write("", TObject::kOverwrite);
    hTotFil[i]->Write("", TObject::kOverwrite);
  }
  hEnrgTot->Write("",TObject::kOverwrite);
  hTedet->Write("", TObject::kOverwrite);
  hProfT->Write("", TObject::kOverwrite);  
}

void pedGainCalib::GetPeds() {
  // Calculate the pedestal
  double std[5];
  for (int stage = 0; stage < 5; stage++) {
    float xLow, xHigh;
    float max, xmax, xpeak;
    max   = hPed[stage]->GetMaximum();
    xmax  = hPed[stage]->GetBinCenter(hPed[stage]->GetMaximumBin());
    TFitResultPtr r = hPed[stage]->Fit("gaus","QNS","",xmax-100,xmax+100); // Q means quiet, R is using TF1 range
    if(int(r)==0) // Everything passed
    {
      Ped[stage]  = r->Parameter(1);// mean
      std[stage]  = r->Parameter(2);// std
    }
    else{  // sanity
      Ped[stage] = 0.;
      cout << "pedGainCalib::getPeds ERROR: could not find the peak of the pedestal distribution for stage ************" << stage << endl;
    }
    cout<< inFileName_s <<" "<<hPed[stage]->Integral()<< " getPeds: measured pedestal for stage " << stage << " is " << Ped[stage] << " with a width of "<<std[stage]<<endl;
  }
}

void pedGainCalib::FillGains(float Vedet, float Tedet, float Ene[5], int phSum[5]) {// Called for each event in the temporary file of proton histories -- gain recalibration
  float Esum = 0.;
  //between t1 and t2 or between t3 and t4
  if ((Tedet > emtrng1 && Tedet < emtrng2) || (Tedet > emtrng3 && Tedet < emtrng4)) { // Analyze full-energy protons
                                                                                      // outside of the phantom region
    if (fabs(Vedet) < 40.) { // not outside the detector
      for (int stage = 0; stage < 5; stage++) {
        float stgEne = Ene[stage];
        hEnrg[stage]->Fill(stgEne);
        Esum += Ene[stage];
	hPed_Out[stage]->Fill(phSum[stage]);
	
      }
      hEnrgTot->Fill(Esum);
      hTedet->Fill(Tedet);
    }
  }

  else if(Tedet<-50){ // Analyze full energy particles inside the phantom region
    if (fabs(Vedet) < 40.) for(int stage =0; stage<5; stage++)hPed_In[stage]->Fill(phSum[stage]);
  }    
  
  if (Ene[0] > 10.) hProfT->Fill(Tedet, Vedet, Ene[0]);
}
void pedGainCalib::GetGains(TVcorrection *TVcorr) {
                            
  // Called prior to the final loop over protons histories to calculate WEPL

  // Calculate the gain
  double Peak[5];
  double std[5];
  for (int stage = 0; stage < 5; stage++) {
    float xLow, xHigh;
    float max, xmax, xpeak,min;
    max   =  hEnrg[stage]->GetMaximum();
    min   =  hEnrg[stage]->GetMinimum();
    xmax  =  hEnrg[stage]->GetBinCenter(hEnrg[stage]->GetMaximumBin());
    TFitResultPtr r =hEnrg[stage]->Fit("gaus","QNS","",xmax-10,xmax+10); // Q means quiet,
    if(int(r)==0) // Everything passed
      {
	Peak[stage]  = r->Parameter(1);// mean
	std[stage]   = r->Parameter(2);//std
      }
      else Peak[stage] = 0.0;
    cout << "getGains: measured peak location for stage " << stage << " is " << Peak[stage] << " and a width of "<<std[stage]<<endl;
    if (Peak[stage] > 0.01) GainFac[stage] = TVcorr->Eempt[stage] / Peak[stage];
    else {
      cout << "getGains ERROR: unable to find a peak position; leaving the gains unchanged. **********" << endl;
      GainFac[stage] = 1.0;
    }
  }

}
