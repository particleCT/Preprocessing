// Encapsulate all of Vladimir Bashkirov pedestal and gain on-the-fly
// calibration
// R.P. Johnson   5/22/2016

#include "pedGainCalib.h"
#include "TFile.h"
pedGainCalib::~pedGainCalib()
{
}
pedGainCalib::pedGainCalib(TFile* root, int pedMin[5], float oldPed[5], float t1, float t2, float t3, float t4, string partType):RootFile(root)
{  
  // Two ranges in t occupied by unobstructed (empty) protons
  emtrng1 = t1; // negative t side
  emtrng2 = t2;
  emtrng3 = t3; // positive t side
  emtrng4 = t4;
  this->partType = partType;
  for (int stage = 0; stage < 5; stage++) {
    Ped[stage] = oldPed[stage];
    GainFac[stage] = 1.0;
    cout << "pedGainCalib setting the default pedestal for stage " << stage << " to " << Ped[stage] << endl;
  }

  // Define histograms for pedestal calculation
  for (int i =0; i<5; i++) hPed[i] = new TH1D(Form("PedestalStage_%i",i), Form("Pedestal region for stage %i",i), 400, pedMin[0], pedMin[0] +400*5); 
  // Define histograms for gain calibration
  if (partType == "H") {
    for (int i =0; i<5; i++) hEnrg[i] = new TH1D(Form("EnergyDistribution_%i",i), Form("Energy Distribution for stage %i",i), 400, 15, 15 + 400*0.175);
    hEnrgTot = new TH1D("SumStageEnergies", "Sum of stage energies", 800, 0, 0 + 0.3*800);
  }
  else {
    for (int i =0; i<5; i++) hEnrg[i] = new TH1D(Form("EnergyDistribution_%i",i), Form("Energy Distribution for stage %i",i), 400, 60, 60 + 400*0.7);
    hEnrgTot = new TH1D("SumStageEnergies", "Sum of stage energies", 800, 0, 0 + 1.2*800);
  }
  // Profile plot to make sure that the phantom does not extend into the regions used for gain calibration
  hProfT = new TProfile2D("Stage0EnergyProfile", "Stage 0 energy profile", 100, -150, -150 + 100*3.0, 100, -50., -50 + 1.0*100);
  hTedet = new TH1D("T_Ions_GainRecalib", "T of ions used for gain recalibration", 100, -150, -150 + 100*3.0);
  RootFile->mkdir("Pedestals");
  
  
}

void pedGainCalib::ClearHist(){
  hProfT->Reset();
  hTedet->Reset();
  hEnrgTot->Reset();
  for (int i =0; i<5; i++){
    hEnrg[i]->Reset();
    hPed[i]->Reset();
  }
}

void pedGainCalib::FillPeds(pCTraw &rawEvt) { // Called for each raw event read in from the input data file
  // Accumulate histograms for pedestal analysis 
  hPed[0]->Fill(rawEvt.enrg_fpga[0].pulse_sum[0]);
  hPed[1]->Fill(rawEvt.enrg_fpga[0].pulse_sum[1]);
  hPed[2]->Fill(rawEvt.enrg_fpga[0].pulse_sum[2]);
  hPed[3]->Fill(rawEvt.enrg_fpga[1].pulse_sum[0]);
  hPed[4]->Fill(rawEvt.enrg_fpga[1].pulse_sum[1]);
}

void pedGainCalib::GetPeds(const char *inFileName, int run_number, int program_version, float proj_angle, int nKeep, string start_time) {
  // Called after completion of  the loop over all input raw data
  // Save plots of the histograms so that the pedestal region can be visualized  
  std::string inFileName_s = inFileName;
  inFileName_s = inFileName_s.substr(inFileName_s.find_last_of("\\/") +1,  inFileName_s.size());  
  RootFile->mkdir(Form("Pedestals/%s",inFileName_s.c_str()));
  RootFile->cd(Form("Pedestals/%s",inFileName_s.c_str()));
  // Calculate the pedestal
  for(int i =0; i<5; i++) hPed[i]->Write("", TObject::kOverwrite);
  for (int stage = 0; stage < 5; stage++) {
    float xLow, xHigh;
    //xLow  =  hPed[stage]->GetBinCenter(hPed[stage]->FindFirstBinAbove( hPed[stage]->GetMaximum()/2));
    //xHigh =  hPed[stage]->GetBinCenter(hPed[stage]->FindLastBinAbove(  hPed[stage]->GetMaximum()/2));
    //hPed[stage]->GetXaxis()->SetRange(xLow,xHigh);
    if ( hPed[stage]->GetEntries()>100 ) Ped[stage] = hPed[stage]->GetBinCenter(hPed[stage]->GetMaximumBin());
    else if(xLow > xHigh || xLow < 0 || xHigh < 0){  // sanity
      Ped[stage] = 0.;
      cout << "pedGainCalib::getPeds ERROR: could not find the peak of the pedestal distribution for stage ************" << stage << endl;
    }
    
    //hPed[stage]->GetXaxis()->SetRange();
    cout << inFileName_s <<" "<<hPed[stage]->Integral()<< " getPeds: measured pedestal for stage " << stage << " is " << Ped[stage] << endl;
  }
}

void pedGainCalib::FillGains(float Vedet, float Tedet, float Ene[5]) { // Called for each event in the temporary file of proton histories -- gain recalibration
  float Esum = 0.;
  //between t1 and t2 or between t3 and t4
  if ((Tedet > emtrng1 && Tedet < emtrng2) || (Tedet > emtrng3 && Tedet < emtrng4)) { // Analyze full-energy protons
                                                                                      // outside of the phantom region
    if (fabs(Vedet) < 40.) {
      for (int stage = 0; stage < 5; stage++) {
        float stgEne = Ene[stage];
        hEnrg[stage]->Fill(stgEne);
        Esum += Ene[stage];
      }
      hEnrgTot->Fill(Esum);
      hTedet->Fill(Tedet);
    }
  }
  if (Ene[0] > 10.) hProfT->Fill(Tedet, Vedet, Ene[0]);
}
void pedGainCalib::GetGains(TVcorrection *TVcorr, const char *inFileName, int run_number, int program_version, int proj_angle, int nKeep, string start_time) {
                            
  // Called prior to the final loop over protons histories to calculate WEPL
  // Save plots of the histograms so that the gains can be visualized
  // Create a folder with the filename
  // Write the files
  std::string inFileName_s = inFileName;
  inFileName_s = inFileName_s.substr(inFileName_s.find_last_of("\\/") +1,  inFileName_s.size());  
  RootFile->mkdir(Form("Pedestals/%s",inFileName_s.c_str()));
  RootFile->cd(Form("Pedestals/%s",inFileName_s.c_str()));
  for(int i =0; i<5; i++) hEnrg[i]->Write("", TObject::kOverwrite);
  hEnrgTot->Write("",TObject::kOverwrite);
  hTedet->Write("", TObject::kOverwrite);
  hProfT->Write("", TObject::kOverwrite);  

  // Calculate the gain
  double Peak[5];
  for (int stage = 0; stage < 5; stage++) {
    float xLow, xHigh;
    //xLow       =  hEnrg[stage]->GetBinCenter(hEnrg[stage]->FindFirstBinAbove( hEnrg[stage]->GetMaximum()/2));
    //xHigh      =  hEnrg[stage]->GetBinCenter(hEnrg[stage]->FindLastBinAbove(  hEnrg[stage]->GetMaximum()/2));
    //hEnrg[stage]->GetXaxis()->SetRange(xLow,xHigh);
    if ( hEnrg[stage]->GetEntries()>100 )  Peak[stage] = hEnrg[stage]->GetBinCenter(hEnrg[stage]->GetMaximumBin());
    else if(xLow > xHigh || xLow < 0 || xHigh < 0) Peak[stage] = 0.0;
    //hEnrg[stage]->GetXaxis()->SetRange();
    cout << "getGains: measured peak location for stage " << stage << " is " << Peak[stage] << endl;
    if (Peak[stage] > 0.01) GainFac[stage] = TVcorr->Eempt[stage] / Peak[stage];
    else {
      cout << "getGains ERROR: unable to find a peak position; leaving the gains unchanged. **********" << endl;
      GainFac[stage] = 1.0;
    }
  }

}
