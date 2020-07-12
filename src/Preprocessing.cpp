// Top-level driving routines for the pCT preprocessing task
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cstring>
#include <cstdint>
#include <cstddef>
#include <thread>
#include <cmath>
#include <ctime>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "Wepl.h"
#include "pCTcut.h"
#include "Preprocessing.h"
#include "BadEvent.h"
using namespace std;
Preprocessing::Preprocessing(){
  theConfig = pCTconfig::GetInstance();
  cout << "*********** Entering the driver program for pCT preprocessing **************" << endl;
  energyOutput = 0;
  timeStampOutput = 0;
  eventIDOutput = 0;
  start_time = time(NULL);
  now = localtime(&start_time);
  printf("Current local time and date: %s", asctime(now));

  if (theConfig->item_int["continuous"])  cout << "A continuous scan will be analyzed\n. The initial stage angle, at time 0, is " << initialAngle << " degrees." << endl;
  else {
    cout << "A fixed angle scan will be analyzed\n";
    if (theConfig->item_float["projection"] > -360.0)
      cout << "The assumed projection angle of " << theConfig->item_float["projection"] << " will override what comes from the data file " << endl;
  }
  cout << "The file will be split into  " << theConfig->item_int["bins"] << " sub-files " << endl;
  if (theConfig->item_int["max_events"] > 0) cout << "The preprocessing will halt after processing " << theConfig->item_int["max_events"] << " events\n";

  //Setting the Filename and creating the output file
  size_t found1 = theConfig->item_str["inputFileName"].find_last_of("/\\"); // for removing absolute path in file (works for both windows and linux)
  size_t found2 = theConfig->item_str["inputFileName"].find('.');
  string outputFile;
  if(found1!=string::npos) outputFile = theConfig->item_str["outputDir"] + theConfig->item_str["inputFileName"].substr(found1,found2-found1)+".root";
  else outputFile = theConfig->item_str["outputDir"] + "/" + theConfig->item_str["inputFileName"].substr(0,theConfig->item_str["intputFileName"].size()-4) + ".root";
   
  pCTcalibRootFile = new TFile(theConfig->item_str["calib"].c_str());
  theTVcorr        = new TVcorrection(pCTcalibRootFile, 0);  
  projectionROOT   = new TFile(outputFile.c_str(),"recreate");
  theCuts          = new pCTcut();// Initialize the code for event selection
  theEvtRecon      = new EvtRecon();
  theWEPL          = new Wepl(pCTcalibRootFile);
  theGeometry      = new pCTgeo(theConfig->item_float["wedgeoffset"]);

  float wedgeLimit = theGeometry->getTWedgeBreaks(4) + 25.0; // NO BRICK OFFSET //+5.
  float openRange  = 20.0;  
  /*float t1 = -100.; // These define two ranges for finding protons passing through zero phantom material, for gain calibration
  float t2 = -theConfig->item_float["sizeleft"]; 
  float t3 = theConfig->item_float["sizeright"];// + theConfig->item_float["wedgeoffset"];
  float t4 = 100.;*/
  float t1 = -150.; // These define two ranges for finding protons passing through zero phantom material, for gain calibration
  float t2 = -151; 
  float t3 = wedgeLimit;
  float t4 = wedgeLimit + openRange;
  float pedestals[5];
  for (int stage = 0; stage < 5; stage++) pedestals[stage] = theTVcorr->ped[stage];
  int pdstlr[5];
  for (int stage = 0; stage < nStage; stage++) pdstlr[stage] = theConfig->item_int[Form("pedrng%d",stage)];
  theCalibration = new pedGainCalib(projectionROOT, pdstlr, pedestals, t1, t2, t3, t4);
};
// ******************************* ******************************* *******************************
// end of the Preprocessing constructor
// ******************************* ******************************* *******************************
//*********** Driving program for pCT preprocessing **************
int Preprocessing::ProcessFile(float fileFraction, int numbTkrFPGA, int numbEdetFPGA) {
  cout << "Preprocessing.cpp: Entering the driver routine for pCT preprocessing. . ." << endl;
  cout << "The phantom is assumed to be less than " << theConfig->item_float["sizeright"] << " in radius, for recalibration of gains." << endl;
  cout << "The wedge phantom offset is assumed to be " << theConfig->item_float["wedgeoffset"] << endl;
  cout << "There are " << numbTkrFPGA << " tracker FPGAs and " << numbEdetFPGA << " energy detector FPGAs" << endl;
  cout << "The file fraction to use is " << fileFraction << endl;
  cout << "Reading the input raw data file " << theConfig->item_str["inputFileName"] << endl;
  in_file = fopen(theConfig->item_str["inputFileName"].c_str(), "rb");  
  if (in_file == NULL) {
    perror("Error opening the input raw data file.");
    exit(1);
  }

  fseek(in_file, 0L, SEEK_END);
  file_size = ftell(in_file);
  rewind(in_file);  
  cout << "Input raw data file size=" << file_size << endl;
  if (fileFraction > 1.0) fileFraction = 1.0;
  // Divide the file into pieces, one for each file
  float fractSize = fileFraction * static_cast<float>(file_size);
  size_t sizeToUse = static_cast<size_t>(fractSize);
  cout << "Preprocessing::ProcessFile, file_size=" << file_size << " sizeToUse=" << sizeToUse << endl;
  size_t fileSize = sizeToUse;
  rewind(in_file);
  // Create an instance of the class for parsing and storing the raw data from the input file
  pCTraw rawEvt(in_file, fileSize, 0, numbTkrFPGA, numbEdetFPGA); 
  rawEvt.readRunHeader(theConfig->item_str["inputFileName"].c_str()); // Look for the run header bits and parse them
  cout << "Preprocessing.cpp: The output directory is " << theConfig->item_str["outputDir"] << endl;
  // Check whether the specified stage angle agrees with what is in the data file
  if (!theConfig->item_int["continuous"]) {
    if (theConfig->item_float["projection"] > -360.0) {
      if (abs(theConfig->item_float["projection"] - rawEvt.stage_angle) / theConfig->item_float["projection"] > 0.001) {
        cout << "Preprocessing.cpp: The provided projection angle does not match the input file run header!\n";
        cout << "The provided projection angle = " << theConfig->item_float["projection"] << endl;
        cout << "The stage angle from the file = " << rawEvt.stage_angle << endl;
        cout << "We are overriding the value from the input file run header.\n";
      }
    } else {
      theConfig->item_float["projection"] = rawEvt.stage_angle;
      cout << "Preprocessing.cpp: We are setting the projection angle according to the input file value of " << theConfig->item_float["projection"] << endl;
      cout << "The input file in general should contain the true reading from the stage for non-continuous-scan runs.\n";
    }
  }
  int year, month, day;
  if (rawEvt.parseDate(year, month, day)) {
    cout << "Preprocessing.cpp: Parsing the run start time date from " << rawEvt.start_time << endl;
    cout << "    Year = " << year << endl << "    Month= " << month << endl << "    Day=   " << day << endl;
  } else cout << "Preprocessing.cpp: Was not able to parse the run date from " << rawEvt.start_time << endl;    
  /////////////////////////////////////////////////////////////////
  // Call the routine that reads the data and does the analysis.
  /////////////////////////////////////////////////////////////////
  theConfig->item_int["doGains"] = true;   // Equivalent to recalibrate
  theEvtRecon->ReadInputFile(theGeometry, theTVcorr, theConfig->item_str["inputFileName"], theCalibration);
  theCalibration->WriteHist();

  // Prepare the ROOT File header
  projectionROOT->cd();
  TTree* header;
  char magic_number[] = "PCTD";
  const char *PREPARED_BY = getenv("USER");
  if (PREPARED_BY == NULL) { // The getenv fails in Windows
    std::string prdby = "Dolittle";
    PREPARED_BY = prdby.c_str();
  }
  float versionNumber = Version;
  int version_id = 0;
  if (energyOutput) version_id += 10;
  if (timeStampOutput) version_id += 100;
  if (eventIDOutput) version_id += 1000;

  string data_source_string = theConfig->item_str["inputFileName"];
  string study_name_string = theConfig->item_str["study"];
  string prepared_by_string = string(PREPARED_BY);
  int current_time      = time(NULL);
  int recalibrate       = theConfig->item_int["recalibrate"];
  int study_name_size = theConfig->item_str["study"].size(); 
  int data_source_size  = theConfig->item_str["inputFileName"].size();
  int prepared_by_size  = strlen(PREPARED_BY);
  
  header = new TTree("header", "meta-data");
  header->Branch("beamEnergy",&beamEnergy,"beamEnergy/F");
  header->Branch("recalibrate",&recalibrate,"recalibrate/I");
  header->Branch("study_date",&rawEvt.study_date,"study_date/I");
  header->Branch("preprocess_date",&current_time,"preprocess_date/I");
  header->Branch("study_name_size",&study_name_size,"study_name_size/I");
  header->Branch("data_source_size",&data_source_size,"data_source_size/I");
  header->Branch("prepared_by_size",&prepared_by_size,"prepared_by_size/I");
  header->Branch("study_name",&study_name_string);
  header->Branch("data_source",&data_source_string);
  header->Branch("prepared_by",&prepared_by_string);
  for (int stage = 0; stage < nStage; ++stage) header->Branch(Form("Gain_%d",stage),&theCalibration->GainFac[stage]);
  header->Fill();
  // Prepare the root file phasespace
  TTree* phase;
  float Ene[5], Wet, theta;
  float Vhit[4], Thit[4];
  float x0,y0,z0;
  float x1,y1,z1;
  float px0,py0,pz0;
  float px1,py1,pz1;
  float E_tot;
  Int_t MaxEnergyTransFilter, ThresholdFilter, dEEFilter;
  unsigned int timeStamp;
  phase = new TTree("phase", "bin tree");
  phase->Branch("t", &Thit, "t[4]/F");  
  phase->Branch("v", &Vhit, "v[4]/F");
  phase->Branch("u", &Uhit, "u[4]/F");
  phase->Branch("wepl", &Wet, "wepl/F");
  phase->Branch("theta", &theta, "theta/F");
  phase->Branch("MaxEnergyTransFilter",&MaxEnergyTransFilter,"MaxEnergyTransFilter/I");
  phase->Branch("ThresholdFilter",&ThresholdFilter,"ThresholdFilter/I");
  phase->Branch("dEEFilter",&dEEFilter,"dEEFilter/I");

  if(!theConfig->item_int["CTOutput"]){
      phase->Branch("ADC", &ADC, "ADC[5]/I");
      phase->Branch("timeStamp", &timeStamp, "timeStamp/I");
      phase->Branch("E", &Ene, "E[5]/F");
      phase->Branch("E_tot", &E_tot, "E_tot/F");      
      
      phase->Branch("x0",&x0,"x0/F");
      phase->Branch("y0",&y0,"y0/F");
      phase->Branch("z0",&z0,"z0/F");
      phase->Branch("x1",&x1,"x1/F");
      phase->Branch("y1",&y1,"y1/F");
      phase->Branch("z1",&z1,"z1/F");
      phase->Branch("px0",&px0,"px0/F");
      phase->Branch("py0",&py0,"py0/F");
      phase->Branch("pz0",&pz0,"pz0/F");
      phase->Branch("px1",&px1,"px1/F");      
      phase->Branch("py1",&py1,"py1/F");
      phase->Branch("pz1",&pz1,"pz1/F");
  }
  int nBadWEPL = 0, nBadTimeStamp = 0, nMaxTrans = 0, nThreshold = 0, ndEEFilter = 0, nTot = 0;
  unsigned int timeStampOld = 0;
  long long timeStampOffset = 0;
  float V[2], T[2], Ut[2], Uv[2], Uft[2], Ufv[2], Tf[2], Vf[2];
  Uft[0] = theEvtRecon->uhitT[0]; Uft[1] = theEvtRecon->uhitT[1];
  Ufv[0] = theEvtRecon->uhitV[0]; Ufv[1] = theEvtRecon->uhitV[1];
  Ut[0]  = theEvtRecon->uhitT[2];  Ut[1] = theEvtRecon->uhitT[3];
  Uv[0]  = theEvtRecon->uhitV[2];  Uv[1] = theEvtRecon->uhitV[3];
  
  for (int EvtNum = 0; EvtNum < theEvtRecon->nEvents; ++EvtNum) {
    Event thisEvent;
    thisEvent = theEvtRecon->evtList[EvtNum];
    dEEFilter = 1;
    MaxEnergyTransFilter = 1;
    ThresholdFilter = 1;

    thisEvent = theEvtRecon->evtList[EvtNum];
    Tf[0] = thisEvent.Thit[0]; Tf[1] = thisEvent.Thit[1];// Front tracker coordinates
    Vf[0] = thisEvent.Vhit[0]; Vf[1] = thisEvent.Vhit[1];
    T[0]  = thisEvent.Thit[2];  T[1] = thisEvent.Thit[3];// Rear tracker coordinates
    V[0]  = thisEvent.Vhit[2];  V[1] = thisEvent.Vhit[3];        
    // Check for a flaky time stamp, and watch out for roll-over of the time stamp counter (after about 10 minutes)!
    // Before V65 of the event builder firmware there were frequent overflows of the time stamp buffer, causing decreasing values for short times    
    if (thisEvent.timestamp < timeStampOld) {
      nBadTimeStamp++;
      if (timeStampOld - timeStamp > 10000000.) {
	cout << "***** Preprocessing: the time stamp decreased a lot since the previous event; we will assume that the counter rolled  over.";
	cout << "  Previous time stamp = " << timeStampOld << "  Time stamp = " << timeStamp << "   Time stamp offset = " << timeStampOffset << endl;
	timeStampOffset = timeStampOffset + pow(2, 36);
      }
    }
    timeStampOld = thisEvent.timestamp;
    if (theConfig->item_int["continuous"]) {
      long long longTimeStamp = 16 * ((long long)thisEvent.timestamp);
      theta = ((float)(longTimeStamp + timeStampOffset)) * theGeometry->timeRes() * theGeometry->stageSpeed() + initialAngle;
    } 
    else  theta = proj_angle;    
    bool inBounds;
    int nGood = 0;
    for (int stage = 0; stage < 5; stage++) {
      float Tedet = theGeometry->extrap2D(Ut, T, theGeometry->energyDetectorU(stage));
      float Vedet = theGeometry->extrap2D(Uv, V, theGeometry->energyDetectorU(stage));
      inBounds = true;
      float TVCorrFactor = theTVcorr->corrFactor(stage, Tedet, Vedet, inBounds);
      Ene[stage] = theCalibration->GainFac[stage] * (thisEvent.ADC[stage] - theCalibration->Ped[stage]) * TVCorrFactor;
      if(inBounds) nGood++;
    }
    Wet = theWEPL->EtoWEPL(Ene, MaxEnergyTransFilter, ThresholdFilter, dEEFilter,E_tot); // Energy to WEPL conversion
    if(!dEEFilter) ++ndEEFilter;
    if(!ThresholdFilter) ++nThreshold;
    if(!MaxEnergyTransFilter) ++nMaxTrans;
    if (Wet < 0. || Wet > 999.) ++nBadWEPL;
    if(Wet > 0 && Wet < 260  && MaxEnergyTransFilter && ThresholdFilter && dEEFilter) theCalibration->FillADC(ADC);
    else ++nTot;
    x0   = Ufv[1]; y0 = Tf[1]; z0 = Vf[1];
    x1   = Uv[0] ; y1 = T[1];  z1 = V[1];    
    px0  = Ufv[1] -  Ufv[0];  py0  = Tf[1]   -  Tf[0]; pz0  = Vf[1]   -  Vf[0];
    px1  = Uv[1]  -  Uv[0];   py1  = T[1]    -  T[1]; pz1  =  V[1]    -  V[0];

    Uhit[0] = Ufv[0]; Uhit[1] = Ufv[1]; Uhit[2] = Uv[0]; Uhit[3] = Uv[1];
    Thit[0] = Tf[0];  Thit[1] = Tf[1];  Thit[2] = T[0];  Thit[3] = T[1];
    Vhit[0] = Vf[0];  Vhit[1] = Vf[1];  Vhit[2] = V[0];  Vhit[3] = V[1];    
    float Length_0 = sqrt(px0*px0 + py0*py0 + pz0*pz0);
    float Length_1 = sqrt(px1*px1 + py1*py1 + pz1*pz1);

    px0 /=  Length_0; py0 /= Length_0; pz0 /= Length_0;
    px1 /=  Length_1; py1 /= Length_1; pz1 /= Length_1;
    phase->Fill();    
    if (EvtNum % 1000000 == 0) cout << " Processing event " << EvtNum << ", time stamp=" << thisEvent.timestamp <<" angle=" << theta << endl;
  }
  theWEPL->WriteHist(projectionROOT);
  cout << "Preprocessing.cpp: The total number of events saved for output was " << theEvtRecon->nEvents << endl;
  cout << "                   The total number of events rejected with was " << nTot << endl;
  cout << "                   The number of events rejected with bad WEPL was " << nBadWEPL << endl;
  cout << "                   The number of events rejected by the Max Transfer Filter "<<nMaxTrans<<endl;
  cout << "                   The number of events rejected by the Threshold Filter "<<nThreshold<<endl;
  cout << "                   The number of events rejected by the dE-E Filter "<<ndEEFilter<<endl;
  cout << "                   The number of decreasing time stamps was " << nBadTimeStamp << endl;
  cout << "                   The accumulated time-stamp correction is, in 10ns units, " << timeStampOffset << endl << endl;
  cout << "Preprocessing.cpp: Write binary file for Output Filename : " << projectionROOT->GetName() << endl;
  projectionROOT->cd();
  header->Write("",TObject::kOverwrite);
  phase->Write("",TObject::kOverwrite);
  projectionROOT->Close();
  time_t end_time = time(NULL);
  now = localtime(&end_time);
  printf("Preprocessing.cpp: Local time and date at end of execution: %s", asctime(now));
  float seconds = difftime(end_time, start_time);
  cout << "Preprocessing.cpp: The total time lapse during execution was " << seconds << " seconds.\n";
  cout << "Preprocessing.cpp: pCT_Preprocessing is all done, including output " "of the projection data." << endl;
  delete theCalibration;
  delete theWEPL; 
  return 0;
}
