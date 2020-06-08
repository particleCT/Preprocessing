// Top-level driving routines for the pCT preprocessing task
// R.P. Johnson   September 15, 2016
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

  if (theConfig->item_int["continuous"]) {
    cout << "A continuous scan will be analyzed\n";
    cout << "The initial stage angle, at time 0, is " << initialAngle << " degrees." << endl;
  } else {
    cout << "A fixed angle scan will be analyzed\n";
    if (theConfig->item_float["projection"] > -360.0)
      cout << "The assumed projection angle of " << theConfig->item_float["projection"] << " will override what comes from the data file " << endl;
  }
  cout << "The file will be split into  " << theConfig->item_int["bins"] << " sub-files " << endl;

  if (theConfig->item_int["max_events"] > 0) cout << "The preprocessing will halt after processing " << theConfig->item_int["max_events"] << " events\n";


  pCTcalibRootFile = new TFile(theConfig->item_str["calib"].c_str());
  TString filename = Form("%s/%s.root",
			  theConfig->item_str["outputDir"].c_str(),
			  theConfig->item_str["inputFileName"].substr(7, theConfig->item_str["inputFileName"].size()-4).c_str());
			  
  projectionROOT = new TFile(filename,"recreate");
  theCuts = new pCTcut();// Initialize the code for event selection  
};
// ******************************* ******************************* *******************************
// end of the Preprocessing constructor
// ******************************* ******************************* *******************************
// Routine called to read the raw data, analyze it, write results to a temporary file, and analyze the WEPL calibration pedestals and gains.  
void Preprocessing::pCTevents(pCTgeo* Geometry, pCTraw rawEvt, pedGainCalib *Calibrate, double Uhit[]) {
  cout << "****************** Executing pCTevents *******************************" << endl;
  char outBuff[93]; // Buffer for writing or reading the temporary file
  int nBuffBytes = 8 * sizeof(float) + 7 * sizeof(int) + sizeof(unsigned char);
  cout << "PCTevents buffer size for file writing is " << nBuffBytes << " bytes." << endl;
  memset(outBuff, 0, nBuffBytes);
  string tempfile = theConfig->item_str["outputDir"] + "/extracted_data_0d.tmp";
  FILE *fptmp;
  cout << "Opening the temporary file " << tempfile << endl;
    fptmp = fopen(tempfile.c_str(), "wb");
    if (fptmp == NULL) {
      perror("Failed to open the temporary file");
      cout << "Program is terminating because the temporary file could not be "
              "opened." << endl;
      exit(1);
    }

  /////////////////////////////////////////////////////////////////////////////////////////
  // Loop over the events in the input data file and store the raw data in the
  // pCTraw class
  /////////////////////////////////////////////////////////////////////////////////////////

  int nErrDmp = 0;
  while (!rawEvt.stop_reading) { // event loop
    try {
      bool debug = rawEvt.event_counter < theConfig->item_int["n_debug"];
      bool Eureka = rawEvt.findEvtHdr(debug); // Search for the bits that indicate the beginning of an event
      if (!Eureka) {
        cout << "pCTevents event header not found after " << rawEvt.event_counter << " events.\n";
        break;
      }
      /////////////////////////////////////////////////
      // Call the method to unpack the raw event data
      /////////////////////////////////////////////////

      rawEvt.readOneEvent(debug);
      Calibrate->FillPeds(rawEvt); // Accumulating histograms for pedestal analysis
      bool daqErr = rawEvt.bad_fpga_address || rawEvt.tag_mismatch || rawEvt.CRC_error || rawEvt.chip_error || rawEvt.bad_strip_address;
      if (daqErr) nErrDmp++; 
      if (debug || (nErrDmp < 10 && daqErr)) {
        if (daqErr) {
          cout << " dumping event with DAQ error: Bad FPGA=" << rawEvt.bad_fpga_address << " Tag mismatch=" << rawEvt.tag_mismatch;
          cout << " CRC error=" << rawEvt.CRC_error << " Chip error=" << rawEvt.chip_error << " Bad strip=" << rawEvt.bad_strip_address << endl;
        }
        rawEvt.dumpEvt(); // Detailed print-out of the raw data
      }
      if (rawEvt.DAQ_error) { // Don't try to reconstruct error events
        // Decide whether to read another event
        cout << "Preprocessing skipping reconstruction at event count " << rawEvt.event_counter << " due to DAQ error\n";
        rawEvt.doWeStop(theConfig->item_int["max_events"], theConfig->item_int["max_time"]); // This will set the stop_reading flag inside the rawEvt instance
        if (rawEvt.stop_reading) cout << "Preprocessing stopping after " << rawEvt.event_counter << " events.\n";
        continue;
      }
      unsigned int timeStampOut = rawEvt.time_tag / 16; // Reduced precision to fit into 32 bits
      unsigned int eventIdOut = rawEvt.event_number;

      // EXTRACT THE TRACKER COORDINATE INFORMATION
      TkrHits pCThits(rawEvt, Geometry, debug);

      // PERFORM THE TRACKER PATTERN RECOGNITION
      pCT_Tracking pCTtracks(pCThits, Geometry);
      // WRITE OUT ONLY EVENTS THAT ARE SUITABLE FOR IMAGE RECONSTRUCTION
      if (theCuts->cutEvt(pCTtracks, pCThits)){
	  if (rawEvt.event_counter % 100000 == 0)
	    cout << "Event " << rawEvt.event_number << ", GoodEvtCount " << theCuts->nKeep << " timeTag " << rawEvt.time_tag << endl;
	  // Write the good events out to a temporary file
      float Vhit[4], Thit[4];
      if (pCTtracks.nTracks == 1) {
	for (int lyr = 0; lyr < 4; lyr++) {
	  Vhit[lyr] = pCTtracks.VTracks[pCTtracks.itkV].X[lyr];
	  Uhit[lyr] = pCTtracks.VTracks[pCTtracks.itkV].U[lyr];
	  if (lyr < 2) Thit[lyr] = pCTtracks.frontPredT(pCTtracks.itkT, Uhit[lyr]); // Extrapolate the T tracks to the U planes occupied by the V layers
	  else Thit[lyr] = pCTtracks.backPredT(pCTtracks.itkT, Uhit[lyr]);
	}
      } else { // Use just the first hit in each layer if there is no track
	double UhitT[4], ThitD[4];
	for (int lyr = 0; lyr < 4; lyr++) {
	  if (pCThits.Lyr[lyr].N[0] > 0) {
	    Vhit[lyr] = pCThits.Lyr[lyr].Y[0].at(0);
	    Uhit[lyr] = pCThits.Lyr[lyr].U[0].at(0);
	  } else {
	    Vhit[lyr] = 0.;
	    Uhit[lyr] = -999.;
	  }
	  if (pCThits.Lyr[lyr].N[1] > 0) {
	    ThitD[lyr] = pCThits.Lyr[lyr].Y[1].at(0);
	    UhitT[lyr] = pCThits.Lyr[lyr].U[1].at(0);
	  } else {
	    ThitD[lyr] = 0.;
	    UhitT[lyr] = -999.;
	  }
	}
	Thit[0] = ThitD[0];
	Thit[1] = ThitD[1];
	if (Uhit[2] > -900. && Uhit[3] > -900. && UhitT[2] > -900. && UhitT[3] > -900.) {// Complete rear tracker vector
	  Thit[2] = Geometry->extrap2D(&UhitT[2], &ThitD[2], Uhit[2]); // Displace the T measurements
	  Thit[3] = Geometry->extrap2D(&UhitT[2], &ThitD[2], Uhit[3]); // to same U as V measurements
	} else {
	  Thit[2] = ThitD[2];
	  Thit[3] = ThitD[3];
	}
      }
      
      unsigned char mask[5] = { 0x01, 0x02, 0x04, 0x08, 0x10 };
      unsigned char OTR = 0;
      if (rawEvt.enrg_fpga[0].OTR[0]) OTR = OTR | mask[0];
      if (rawEvt.enrg_fpga[0].OTR[1]) OTR = OTR | mask[1];
      if (rawEvt.enrg_fpga[0].OTR[2]) OTR = OTR | mask[2];
      if (rawEvt.enrg_fpga[1].OTR[0]) OTR = OTR | mask[3];
      if (rawEvt.enrg_fpga[1].OTR[1]) OTR = OTR | mask[4];
      
      ADC[0] = rawEvt.enrg_fpga[0].pulse_sum[0];
      ADC[1] = rawEvt.enrg_fpga[0].pulse_sum[1];
      ADC[2] = rawEvt.enrg_fpga[0].pulse_sum[2];
      ADC[3] = rawEvt.enrg_fpga[1].pulse_sum[0];
      ADC[4] = rawEvt.enrg_fpga[1].pulse_sum[1];
      
      memcpy(&outBuff[0], &timeStampOut, sizeof(int));
      memcpy(&outBuff[sizeof(int)], &eventIdOut, sizeof(int));
      memcpy(&outBuff[2 * sizeof(int)], Thit, 4 * sizeof(float));
      memcpy(&outBuff[2 * sizeof(int) + 4 * sizeof(float)], Vhit, 4 * sizeof(float));
      memcpy(&outBuff[2 * sizeof(int) + 8 * sizeof(float)], ADC, 5 * sizeof(int));
      memcpy(&outBuff[8 * sizeof(float) + 7 * sizeof(int)], &OTR, sizeof(unsigned char));
      fwrite(outBuff, sizeof(char), nBuffBytes, fptmp); // Writing the processed event data to a temporary file
	}

      // Decide whether to read another event
      rawEvt.doWeStop(theConfig->item_int["max_events"], theConfig->item_int["max_time"]); // This will set the
      // stop_reading flag inside the rawEvt instance
      if (rawEvt.stop_reading) cout << "Preprocessing stopping after " << rawEvt.event_counter << " events.\n";
    }
    catch (const BadEvent &badEvent) {
      // do nothing - go back and find next event header
      cout << badEvent.what() << endl;
      cout << "Preprocessing: bad event data input caught at event number " << rawEvt.event_counter << endl;
    }
  } // End of the loop over events

  cout << endl << " Number of events with DAQ errors reported = " << nErrDmp << endl;
  theCuts->summary(); // Summarize the event counts
  cout << " V-layer U coordinates = " << Uhit[0] << " " << Uhit[1] << " " << Uhit[2] << " " << Uhit[3] << " assumed the same for all events\n";
  fclose(fptmp);
  /////////////////////////////////////////////////////////
  // PEDESTAL ANALYSIS FOR THE WEPL DETECTOR
  /////////////////////////////////////////////////////////
  if (rawEvt.event_counter > 90000) {
    Calibrate->GetPeds();
    cout <<"We are updating the energy detector pedestal settings to the on-the-fly measurement values.\n";
    for (int stage = 0; stage < 5; stage++) {
      if (Calibrate->Ped[stage] == 0.) cout <<"Keeping the theTVcorr pedestal value " << theTVcorr->ped[stage] << " for stage " << stage << endl;
      else cout <<"Pedestal from calibration file = " << theTVcorr->ped[stage] << "Drift= " << Calibrate->Ped[stage] - theTVcorr->ped[stage] << endl;
    }
  }
  else cout<<"Not enough events to recalibrate the pedestals." << endl;
  for (int stage = 0; stage < 5; stage++) cout <<"The energy detector pedestal for stage " << stage << " is " << Calibrate->Ped[stage] << endl;

  /////////////////////////////////////////////////////////
  // Read PROCESSED DATA BACK FROM THE TEMPORARY FILE TO CALCULATE GAIN CORRECTION FACTORS
  /////////////////////////////////////////////////////////
  if (theConfig->item_int["recalibrate"]) {
    fptmp = fopen(tempfile.c_str(), "rb");
    if (fptmp == NULL) {
      cout <<"Failed to reopen the temporary file\n";
      exit(1);
    }
    fseek(fptmp, 0L, SEEK_END);
    size_t file_size = ftell(fptmp);
    rewind(fptmp);
    cout <<"Temporary data file size=" << file_size << endl;

    for (int EvtNum = 0; EvtNum < theCuts->nKeep; EvtNum++) {
      if (EvtNum % 100000 == 0) {
        cout << "Processing event " << EvtNum
             << " from the temp file for calibration." << endl;
      }
      int ret = fread(outBuff, sizeof(char), nBuffBytes, fptmp);
      float Vhit[4], Thit[4];
      //int ADC[5];
      memcpy(Thit, &outBuff[2 * sizeof(int)], 4 * sizeof(float));
      memcpy(Vhit, &outBuff[2 * sizeof(int) + 4 * sizeof(float)], 4 * sizeof(float));
      memcpy(ADC,  &outBuff[2 * sizeof(int) + 8 * sizeof(float)], 5 * sizeof(int));
      if (EvtNum <= theConfig->item_int["n_debug"]) {
        cout << "pCTevents reading back the temporary file for calibrations. . .\n";
        cout << "  Vhit     Uhit    Thit\n";
        for (int lyr = 0; lyr < 4; lyr++)  cout << "  " << Vhit[lyr] << "  " << Uhit[lyr] << "  " << Thit[lyr] << endl;
        cout << "  Stage pulse sums= ";
        for (int stage = 0; stage < 5; stage++) cout << ADC[stage] << "  ";
        cout << endl;
      }
      // Calculate the calibrated WEPL
      float Ene[5], Vedet[5], Tedet[5];
      double Tback[2], Vback[2];
      double Tfront[2], Vfront[2];
      for (int lyr = 2; lyr < 4; lyr++) {
        Tback[lyr - 2] = Thit[lyr];
        Vback[lyr - 2] = Vhit[lyr];
      }
      for (int lyr = 0; lyr < 2; lyr++) {
        Tfront[lyr] = Thit[lyr];
        Vfront[lyr] = Vhit[lyr];
      }

      float Tphantom = Geometry->extrap2D(Uhit, Tfront, 0.); // Center of the phantom
      int nGood = 0;
      for (int stage = 0; stage < 5; stage++) {
        // Extrapolate the rear track vector to the energy detector stage
        Vedet[stage] = Geometry->extrap2D(&Uhit[2], Vback, Geometry->energyDetectorU(stage));
        Tedet[stage] = Geometry->extrap2D(&Uhit[2], Tback, Geometry->energyDetectorU(stage));
        bool inBounds;
        //Ene[stage] = ((float)ADC[stage] - Calibrate->Ped[stage]) * theTVcorr->corrFactor(stage, Tedet[stage], Vedet[stage], inBounds);
	Ene[stage] = ((float)ADC[stage]) * theTVcorr->corrFactor(stage, Tedet[stage], Vedet[stage], inBounds);
	if (inBounds) nGood++;
      }
      if(nGood==5) Calibrate->FillGains(Vedet[0], Tphantom, Ene, ADC); // Accumulate histograms for gain recalibration 
    }
    Calibrate->GetGains(theTVcorr);
    Calibrate->WriteHist();
    fclose(fptmp);
    cout << "closed the temporary file " << tempfile << endl;
  }
  cout <<"Gain correction factors: ";
  for (int stage = 0; stage < 5; stage++) cout << Calibrate->GainFac[stage] << "  ";
  cout << endl;

  cout << "done with pCTevents.  Number of events kept=" << theCuts->nKeep << endl;
  return;
}

//*********** Driving program for pCT preprocessing **************
int Preprocessing::ProcessFile(float fileFraction, int numbTkrFPGA, int numbEdetFPGA) {
  cout << "Preprocessing.cpp: Entering the driver routine for pCT preprocessing. . ." << endl;
  cout << "The phantom is assumed to be less than " << theConfig->item_float["size"] << " in radius, for recalibration of gains." << endl;
  cout << "The wedge phantom offset is assumed to be " << theConfig->item_float["wedgeoffset"] << endl;
  cout << "There are " << numbTkrFPGA << " tracker FPGAs and " << numbEdetFPGA << " energy detector FPGAs" << endl;
  cout << "The file fraction to use is " << fileFraction << endl;

  //////////////////////////////////////////////////////////
  // Opening the input file
  //////////////////////////////////////////////////////////
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
  pCTgeo* Geometry = new pCTgeo(0.);   // Create a class instance with all of the geometry information

  cout << "Preprocessing.cpp: The output directory is " << theConfig->item_str["outputDir"] << endl;
  // Check whether the specified stage angle agrees with what is in the data file
  if (!theConfig->item_int["continuous"]) {
    if (theConfig->item_float["projection"] > -360.0) {
      if (abs(theConfig->item_float["projection"] - rawEvt.stage_angle) / theConfig->item_float["projection"] > 0.001) {
        cout << "Preprocessing.cpp: The provided projection angle does not "
	  "match the input file run header!\n";
        cout << "The provided projection angle = " << theConfig->item_float["projection"] << endl;
        cout << "The stage angle from the file = " << rawEvt.stage_angle << endl;
        cout << "We are overriding the value from the input file run header.\n";
      }
    } else {
      theConfig->item_float["projection"] = rawEvt.stage_angle;
      cout << "Preprocessing.cpp: We are setting the projection angle "
              "according to the input file value of " << theConfig->item_float["projection"] << endl;
      cout << "The input file in general should contain the true reading from "
              "the stage for non-continuous-scan runs.\n";
    }
  }

  int year, month, day;
  if (rawEvt.parseDate(year, month, day)) {
    cout << "Preprocessing.cpp: Parsing the run start time date from " << rawEvt.start_time << endl;
    cout << "    Year = " << year << endl;
    cout << "    Month= " << month << endl;
    cout << "    Day=   " << day << endl;
  } else cout << "Preprocessing.cpp: Was not able to parse the run date from " << rawEvt.start_time << endl;
    
  ///////////////////////////////////////////// ////////////
  // GET EXISTING WEPL CALIBRATION CONSTANTS  [CEO Aug 2016]
  /////////////////////////////////////////////////////////

  for(int i =0; i<5; i++)StgThr[i] = theConfig->item_float[Form("thr%d",i)];
  Wepl *WEPL = new Wepl(pCTcalibRootFile, projectionROOT);
  theTVcorr = new TVcorrection(pCTcalibRootFile, 0);
  WEPL->SetEthresholds(StgThr[0], StgThr[1], StgThr[2], StgThr[3], StgThr[4]);  
  // Create a vector of pointers to instances of the pedestal and gain calibration class
  float t1 = -100.; // These define two ranges for finding protons passing through zero phantom material, for gain calibration
  float t2 = -100.; // ****** Let's keep this to one side only, for now, to accommodate the wedge calibration runs with bricks
  float t3 = theConfig->item_float["size"];// + theConfig->item_float["wedgeoffset"];
  float t4 = 100.;
  float pedestals[5];
  for (int stage = 0; stage < 5; stage++) pedestals[stage] = theTVcorr->ped[stage];
  int pdstlr[5];
  for (int stage = 0; stage < nStage; stage++) pdstlr[stage] = theConfig->item_int[Form("pedrng%d",stage)];
  pedGainCalib* Calibrate = new pedGainCalib(projectionROOT, pdstlr, pedestals, t1, t2, t3, t4);
  /////////////////////////////////////////////////////////////////
  // Call the routine that reads the data and does the analysis.
  /////////////////////////////////////////////////////////////////
  // Uhit is filled and returned for use below, just to save space in the temporary file.
  pCTevents(Geometry, rawEvt, Calibrate, Uhit);

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
  for (int stage = 0; stage < nStage; ++stage) header->Branch(Form("Gain_%d",stage),&Calibrate->GainFac[stage]);
  header->Fill();

  // Prepare the root file phasespace
  TTree* phase;
  float Ene[5], Wet, theta;
  float Vhit[4], Thit[4];
  float x0,y0,z0;
  float x1,y1,z1;
  float px0,py0,pz0;
  float px1,py1,pz1;
  Int_t MaxEnergyTransFilter, ThresholdFilter, dEEFilter;
  unsigned int timeStamp;
  phase = new TTree("phase", "bin tree");
  phase->Branch("t", &Thit, "t[4]/F");  
  phase->Branch("v", &Vhit, "v[4]/F");
  phase->Branch("u", &Uhit, "u[4]/F");
  phase->Branch("ADC", &ADC, "ADC[5]/I");
  phase->Branch("timeStamp", &timeStamp, "timeStamp/I");
  phase->Branch("E", &Ene, "E[5]/F");
  phase->Branch("wepl", &Wet, "wepl/F");
  phase->Branch("theta", &theta, "theta/F");

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
  phase->Branch("MaxEnergyTransFilter",&MaxEnergyTransFilter,"MaxEnergyTransFilter/I");
  phase->Branch("ThresholdFilter",&ThresholdFilter,"ThresholdFilter/I");
  phase->Branch("dEEFilter",&dEEFilter,"dEEFilter/I");

  char outBuff[93]; // Make this extra large, just in case floats are not 4 bytes in size
  int nBuffBytes = 8 * sizeof(float) + 7 * sizeof(int) + sizeof(unsigned char);
  cout << "Preprocessing.cpp: Buffer size for writing the temporary file is " << nBuffBytes << " bytes\n";
  memset(outBuff, 0, nBuffBytes);

  FILE *fptmp;
  string tempfile;
  int nEvtot = 0;
  int nBadWEPL = 0, nBadTimeStamp = 0, nMaxTrans = 0, nThreshold = 0, ndEEFilter = 0, nTot = 0;
  unsigned int timeStampOld = 0;
  long long timeStampOffset = 0;
  
  tempfile = theConfig->item_str["outputDir"] + "/extracted_data_0d.tmp";
  fptmp = fopen(tempfile.c_str(), "rb");
  if (fptmp == NULL) {
    perror("Preprocessing.cpp: Failed to open the temporary file");
    exit(1);
  }
  for (int EvtNum = 0; EvtNum < theCuts->nKeep; EvtNum++) { // Analyze data from the temporary file event by event    
    ret = fread(outBuff, sizeof(char), nBuffBytes, fptmp);
    dEEFilter = 1;
    MaxEnergyTransFilter = 1;
    ThresholdFilter = 1;


    unsigned int event_id;
    unsigned char OTR;
    memcpy(&timeStamp, &outBuff[0], sizeof(int));
    memcpy(&event_id, &outBuff[sizeof(int)], sizeof(int));
    memcpy(Thit, &outBuff[2 * sizeof(int)], 4 * sizeof(float));
    memcpy(Vhit, &outBuff[2 * sizeof(int) + 4 * sizeof(float)], 4 * sizeof(float));
    memcpy(ADC, &outBuff[2 * sizeof(int) + 8 * sizeof(float)], 5 * sizeof(int));
    memcpy(&OTR, &outBuff[8 * sizeof(float) + 7 * sizeof(int)], sizeof(unsigned char));
    // Check for a flaky time stamp, and watch out for roll-over of the time stamp counter (after about 10 minutes)!
    // Before V65 of the event builder firmware there were frequent overflows
    // of the time stamp buffer, causing decreasing values for short times    
    if (timeStamp < timeStampOld) {
      nBadTimeStamp++;
      if (timeStampOld - timeStamp > 10000000.) {
	cout << "***** Preprocessing: the time stamp decreased a lot since the previous event; we will assume that the counter rolled  over.";
	cout << "  Previous time stamp = " << timeStampOld << "  Time stamp = " << timeStamp << "   Time stamp offset = " << timeStampOffset << endl;
	timeStampOffset = timeStampOffset + pow(2, 36);
      }
    }

    timeStampOld = timeStamp;
    if (theConfig->item_int["continuous"]) {
      long long longTimeStamp = 16 * ((long long)timeStamp);
      theta = ((double)(longTimeStamp + timeStampOffset)) * Geometry->timeRes() * Geometry->stageSpeed() + initialAngle;
    } 
    else  theta = proj_angle;    
    // Calculate the calibrated WEPL
    double Tback[2], Vback[2];
    for (int lyr = 2; lyr < 4; lyr++){
	Tback[lyr - 2] = Thit[lyr];
	Vback[lyr - 2] = Vhit[lyr];
    }

    bool inBounds;
    int nGood = 0;
    for (int stage = 0; stage < 5; stage++) {
      // TVcorrection and MeV conversion.  Extrapolate the rear track vector
      // to the energy detector location.
      float Vedet = Geometry->extrap2D(&Uhit[2], Vback, Geometry->energyDetectorU(stage));
      float Tedet = Geometry->extrap2D(&Uhit[2], Tback, Geometry->energyDetectorU(stage));
      float TVCorrFactor = theTVcorr->corrFactor(stage, Tedet, Vedet, inBounds);
      //Ene[stage] = Calibrate->GainFac[stage] * ((float)ADC[stage] - Calibrate->Ped[stage]) * TVCorrFactor;
      Ene[stage] = Calibrate->GainFac[stage] * ((float)ADC[stage]) * TVCorrFactor;
      if(inBounds) nGood++;

    }

    Wet = WEPL->EtoWEPL(Ene, MaxEnergyTransFilter, ThresholdFilter, dEEFilter); // Energy to WEPL conversion
    if(!dEEFilter) ++ndEEFilter;
    if(!ThresholdFilter) ++nThreshold;
    if(!MaxEnergyTransFilter) ++nMaxTrans;
    if (Wet < 0. || Wet > 999.) ++nBadWEPL;
    if(Wet > 0 && Wet < 260  && MaxEnergyTransFilter && ThresholdFilter && dEEFilter) Calibrate->FillADC(ADC);
    else ++nTot;
    
    x0   = Uhit[1]; y0 = Thit[1]; z0 = Vhit[1];
    x1   = Uhit[2]; y1 = Thit[2]; z1 = Vhit[2];
    
    px0  = Uhit[1] -  Uhit[0];  py0  = Thit[1]   -  Thit[0]; pz0  = Vhit[1]   -  Vhit[0];
    px1  = Uhit[3] -  Uhit[2];  py1  = Thit[3]   -  Thit[2]; pz1  = Vhit[3]   -  Vhit[2];

    //Normalise
    double Length_0 = sqrt(px0*px0 + py0*py0 + pz0*pz0);
    double Length_1 = sqrt(px1*px1 + py1*py1 + pz1*pz1);

    px0 /=  Length_0; py0 /= Length_0; pz0 /= Length_0;
    px1 /=  Length_1; py1 /= Length_1; pz1 /= Length_1;

    phase->Fill();    
    if (EvtNum % 100000 == 0) cout << " Processing event " << EvtNum << " from the temp file, time stamp=" << timeStamp <<" angle=" << theta << endl;
    ++nEvtot;
  }
  Calibrate->WriteHist();
  delete WEPL; // the class
  fclose(fptmp);
  string cmd  = "rm -f " + tempfile;
  cout << "Preprocessing.cpp: Try to delete the temporary file " << tempfile << "with command "<<cmd<<endl;
  ret = system(cmd.c_str());
  cout << "Preprocessing.cpp: The total number of events saved for output was " << nEvtot << endl;
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
  double seconds = difftime(end_time, start_time);
  cout << "Preprocessing.cpp: The total time lapse during execution was " << seconds << " seconds.\n";
  cout << "Preprocessing.cpp: pCT_Preprocessing is all done, including output " "of the projection data." << endl;
  delete Calibrate;
  return 0;
}
