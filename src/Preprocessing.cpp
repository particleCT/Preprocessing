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
#include "arg.h" // Command line options

#include "TTree.h"
#include "TFile.h"
// Effective Aug 2016: Use new WEPL calibration
#include "Wepl.h"
#include "pCTcut.h"
#include "Preprocessing.h"
#include "BadEvent.h"

using namespace std;

Preprocessing::Preprocessing(pCTconfig cfg): config(cfg){
  
  cout << "*********** Entering the driver program for pCT preprocessing **************" << endl;
  energyOutput = 0;
  timeStampOutput = 0;
  eventIDOutput = 0;
  start_time = time(NULL);
  now = localtime(&start_time);
  printf("Current local time and date: %s", asctime(now));

  if (config.item_int["continuous"]) {
    cout << "A continuous scan will be analyzed\n";
     cout << "The initial stage angle, at time 0, is " << initialAngle << " degrees." << endl;
  } else {
    cout << "A fixed angle scan will be analyzed\n";
    if (config.item_float["projection"] > -360.0)
      cout << "The assumed projection angle of " << config.item_float["projection"] << " will override what comes from the data file " << endl;
  }
  cout << "The file will be split into  " << fileBins << " sub-files " << endl;

  if (config.item_int["max_events"] > 0)
    cout << "The preprocessing will halt after processing " << config.item_int["max_events"] << " events\n";

  cout << "Reading the input raw data file " << config.item_str["inputFileName"] << endl;
  in_file = fopen(config.item_str["inputFileName"].c_str(), "rb");
  if (in_file == NULL) {
    perror("Error opening the input raw data file.");
    exit(1);
  }

  fseek(in_file, 0L, SEEK_END);
  file_size = ftell(in_file);
  rewind(in_file);
  cout << "Input raw data file size=" << file_size << endl;


};

// ******************************* ******************************* *******************************
// end of the Preprocessing constructor
// ******************************* ******************************* *******************************
void Preprocessing::WriteBinaryFile(bool timeStampOutput, bool energyOutput, bool eventIDOutput, float AngleNb,
				    const char OutputFilename[], const char DATA_SOURCE[], const char PHANTOM_NAME[],
                                     int study_date, int event_counter, double u[], float V0[], float V1[], float V2[],
                                     float V3[], float T0[], float T1[], float T2[], float T3[], float E1[], float E2[],
                                     float E3[], float E4[], float E5[], float WetBinary[], float ProjAngle[],
                                     unsigned int TimeStamp[], unsigned int EventIDs[])
{
  ofstream data_file;
  data_file.open(OutputFilename, ios::binary | ios::trunc);
  char magic_number[] = "PCTD";
  const char *PREPARED_BY = getenv("USER");
  if (PREPARED_BY == NULL) { // The getenv fails in Windows
    std::string prdby = "Dolittle";
    PREPARED_BY = prdby.c_str();
  }
  float versionNumber = Version;
  float projection_angle = (float)AngleNb; 

  int version_id = 0;
  if (energyOutput)
    version_id += 10;
  if (timeStampOutput)
    version_id += 100;
  if (eventIDOutput)
    version_id += 1000;
  int current_time      = time(NULL);
  int phantom_name_size = strlen(PHANTOM_NAME);
  int data_source_size  = strlen(DATA_SOURCE);
  int prepared_by_size  = strlen(PREPARED_BY);

  // Write headers:
  data_file.write(magic_number, 4);                                            // magic number identifier (note that it doesn't include null terminator '\0')
  data_file.write(reinterpret_cast<char *>(&version_id), sizeof(int));         // format version identifier
  data_file.write(reinterpret_cast<char *>(&event_counter), sizeof(int));      // number of events in file
  data_file.write(reinterpret_cast<char *>(&projection_angle), sizeof(float)); // projection angle
  data_file.write(reinterpret_cast<char *>(&beamEnergy), sizeof(float));    // beam energy
  data_file.write(reinterpret_cast<char *>(&study_date), sizeof(int));         // generation date
  data_file.write(reinterpret_cast<char *>(&current_time), sizeof(int));       // pre-process date
  data_file.write(reinterpret_cast<char *>(&phantom_name_size), sizeof(int));
  data_file.write(PHANTOM_NAME, phantom_name_size); // phantom name or description (string)
  data_file.write(reinterpret_cast<char *>(&data_source_size), sizeof(int));
  data_file.write(DATA_SOURCE, data_source_size); // data source (string)
  data_file.write(reinterpret_cast<char *>(&prepared_by_size), sizeof(int));
  data_file.write(PREPARED_BY, prepared_by_size); // prepared by (string)

  // Write event data:
  data_file.write(reinterpret_cast<char *>(T0), event_counter * sizeof(float));
  data_file.write(reinterpret_cast<char *>(T1), event_counter * sizeof(float));
  data_file.write(reinterpret_cast<char *>(T2), event_counter * sizeof(float));
  data_file.write(reinterpret_cast<char *>(T3), event_counter * sizeof(float));

  data_file.write(reinterpret_cast<char *>(V0), event_counter * sizeof(float));
  data_file.write(reinterpret_cast<char *>(V1), event_counter * sizeof(float));
  data_file.write(reinterpret_cast<char *>(V2), event_counter * sizeof(float));
  data_file.write(reinterpret_cast<char *>(V3), event_counter * sizeof(float));

  float *tmp = new float[event_counter];
  for (int n = 0; n < event_counter; n++) {
    tmp[n] = u[0];
  }
  data_file.write(reinterpret_cast<char *>(tmp), event_counter * sizeof(float));
  for (int n = 0; n < event_counter; n++) {
    tmp[n] = u[1];
  }
  data_file.write(reinterpret_cast<char *>(tmp), event_counter * sizeof(float));
  for (int n = 0; n < event_counter; n++) {
    tmp[n] = u[2];
  }
  data_file.write(reinterpret_cast<char *>(tmp), event_counter * sizeof(float));
  for (int n = 0; n < event_counter; n++) {
    tmp[n] = u[3];
  }
  data_file.write(reinterpret_cast<char *>(tmp), event_counter * sizeof(float));

  if (energyOutput) {
    cout << "WriteBinaryFile3: writing out the stage energy arrays" << endl;
    data_file.write(reinterpret_cast<char *>(E1), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char *>(E2), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char *>(E3), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char *>(E4), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char *>(E5), event_counter * sizeof(float));
  }

  data_file.write(reinterpret_cast<char *>(WetBinary), event_counter * sizeof(float));
  data_file.write(reinterpret_cast<char *>(ProjAngle), event_counter * sizeof(float));

  if (timeStampOutput) {
    cout << "WriteBinaryFile3: writing out the time stamps" << endl;
    data_file.write(reinterpret_cast<char *>(TimeStamp), event_counter * sizeof(unsigned int));
  }

  if (eventIDOutput) {
    cout << "WriteBinaryFile3: writing out the event ids" << endl;
    data_file.write(reinterpret_cast<char *>(EventIDs), event_counter * sizeof(int));
  }
  data_file.close();
  delete[] tmp;
};

// ******************************* ******************************* *******************************
// end of the Preprocessing : WriteBinaryToFile
// ******************************* ******************************* *******************************
void Preprocessing::WriteRootFile(bool timeStampOutput, bool energyOutput, bool eventIDOutput, int fileNb,
				  int study_date, int event_counter, double u[],
				  float V0[], float V1[], float V2[], float V3[],
				  float T0[], float T1[], float T2[], float T3[],
				  float E1[], float E2[], float E3[], float E4[], float E5[],
				  float WetBinary[], float ProjAngle[], unsigned int TimeStamp[], unsigned int EventIDs[]) 
{

  TString filename = Form("%s/%s_%d.root",
			  config.item_str["outputDir"].c_str(),
			  config.item_str["inputFileName"].substr(7, config.item_str["inputFileName"].size()-4).c_str(),
			  fileNb);
  TFile* projectionROOT = new TFile( filename,"update");

  TTree* header;
  TTree* phase;

  float t[4], v[4], E[5], wepl, angle;
  float x0,y0,z0;
  float x1,y1,z1;
  float px0,py0,pz0;
  float px1,py1,pz1;
  Int_t MaxEnergyTransFilter, ThresholdFilter, dEEFilter;
  
  
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

  string data_source_string = config.item_str["inputFileName"];
  string study_name_string = config.item_str["study"];
  string prepared_by_string = string(PREPARED_BY);
  
  int current_time      = time(NULL);
  int recalibrate       = config.item_int["recalibrate"];
  int study_name_size = config.item_str["study"].size(); 
  int data_source_size  = config.item_str["inputFileName"].size();
  int prepared_by_size  = strlen(PREPARED_BY);
  
  header = new TTree("header", "meta-data");
  header->Branch("beamEnergy",&beamEnergy,"beamEnergy/F");
  header->Branch("recalibrate",&recalibrate,"recalibrate/I");
  header->Branch("study_date",&study_date,"study_date/I");
  header->Branch("preprocess_date",&current_time,"preprocess_date/I");
  header->Branch("study_name_size",&study_name_size,"study_name_size/I");
  header->Branch("data_source_size",&data_source_size,"data_source_size/I");
  header->Branch("prepared_by_size",&prepared_by_size,"prepared_by_size/I");
  header->Branch("study_name",&study_name_string);
  header->Branch("data_source",&data_source_string);
  header->Branch("prepared_by",&prepared_by_string);
  header->Fill();
  

  phase = new TTree("phase", "bin tree");
  phase->Branch("t", &t, "t[4]/F");  // duplicate information -- uncomment for DROPTVS
  phase->Branch("v", &v, "v[4]/F");
  phase->Branch("u", &u, "u[4]/F");
  phase->Branch("E", &u, "E[4]/F");
  phase->Branch("wepl", &wepl, "wepl/F");
  phase->Branch("angle", &angle, "angle/F");

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

  for( int i = 0; i < event_counter; i++ ){
    
    dEEFilter = 1;
    MaxEnergyTransFilter = 1;
    ThresholdFilter = 1;

    
    t[0] = T0[i]  ; t[1] = T1[i]; t[2] = T2[i]; t[3] = T3[i];
    v[0] = V0[i]  ; v[1] = V1[i]; v[2] = V2[i]; v[3] = V3[i];
    //u should be already done
    x0   = u[1]; y0 = T1[i]; z0 = V1[i];
    x1   = u[2]; y1 = T2[i]; z1 = V2[i];
    
    px0  = u[1] -  u[0];  py0  = T1[i]   -  T0[i]; pz0  = V1[i]   -  V0[i];
    px1  = u[3] -  u[2];  py1  = T3[i]   -  T2[i]; pz1  = V3[i]   -  V2[i];

    //Normalise
    double Length_0 = sqrt(px0*px0 + py0*py0 + pz0*pz0);
    double Length_1 = sqrt(px1*px1 + py1*py1 + pz1*pz1);

    px0 /=  Length_0; py0 /= Length_0; pz0 /= Length_0;
    px1 /=  Length_1; py1 /= Length_1; pz1 /= Length_1;
    
    
    
    if (energyOutput) {
      cout << "WriteRootFile: writing out the stage energy arrays" << endl;
      E[0] = E1[i]; E[1] = E2[i]; E[2] = E3[i]; E[3] = E4[i]; E[4] = E5[i];
    }
    wepl = WetBinary[i];
    angle = ProjAngle[i];
    if(WetBinary[i] >= 5000){
      wepl = WetBinary[i]-5000;
      dEEFilter = 0;
      MaxEnergyTransFilter = 1;
      ThresholdFilter = 1;
    }

    else if(WetBinary[i] >= 2000 && WetBinary[i] < 5000){
      wepl = WetBinary[i]-2000;
      dEEFilter = 1;
      ThresholdFilter = 0;
      MaxEnergyTransFilter = 1;

    }

    else if(WetBinary[i] >= 1000 && WetBinary[i] < 2000){
      wepl = WetBinary[i]-1000;
      dEEFilter = 1;
      ThresholdFilter = 1;
      MaxEnergyTransFilter = 0;
    }
    else wepl = WetBinary[i];
    phase->Fill();
    
  }
  header->Write("",TObject::kOverwrite);
  phase->Write("",TObject::kOverwrite);
  projectionROOT->Close();
};

// ******************************* ******************************* *******************************
// end of the Preprocessing : WriteRootFile
// ******************************* ******************************* *******************************

// Routine called to read the raw data, analyze it, write results to a temporary
// file, and analyze the WEPL calibration pedestals and gains.  Ideally this would be a private member
// of the Preprocessing class, but then the compiler doesn't allow it to be passed to multiple threads.
void Preprocessing::pCTevents(pCTconfig config, pCTgeo* Geometry, pCTraw rawEvt, pedGainCalib *Calibrate,
               TVcorrection *const TVcorr, int &nKeep, double Uhit[]) {

  // Multiple instances of this program can execute in parallel threads.
  // Uhit is passed by reference so that the U coordinates can be returned by
  // thread 0, to save space by avoiding writing them in the temporary file.
  //      Therefore, each instance needs to be given separate memory locations
  // for Uhit, to avoid writing conflict.
  // All threads return a value nKeep for the number of events passing cuts and
  // written to the temporary file.
  // All threads write the output data into separate temporary files.
  // All threads get their own instance of pedGainCalib, so that each can return
  // the updated pedestals and gains.
  // The other objects are common to all instances and should not be modified.

  cout << "****************** Executing pCTevents *******************************" << endl;

  pCTcut cuts(config);//.threadNum, 3, 5, 8); // Initialize the code for event selection

  char outBuff[93]; // Buffer for writing or reading the temporary file
  int nBuffBytes = 8 * sizeof(float) + 7 * sizeof(int) + sizeof(unsigned char);
  cout << "PCTevents buffer size for file writing is " << nBuffBytes << " bytes."
       << endl;
  memset(outBuff, 0, nBuffBytes);

  string tempfile = config.item_str["outputDir"] + "/extracted_data_0d.tmp";
  FILE *fptmp;
  if (config.item_int["level"] != 0) {
    cout << "Opening the temporary file " << tempfile << endl;
    fptmp = fopen(tempfile.c_str(), "wb");
    if (fptmp == NULL) {
      perror("Failed to open the temporary file");
      cout << "Program is terminating because the temporary file could not be "
              "opened." << endl;
      exit(1);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  // Loop over the events in the input data file and store the raw data in the
  // pCTraw class
  /////////////////////////////////////////////////////////////////////////////////////////

  int nErrDmp = 0;
  while (!rawEvt.stop_reading) {
    try {
      bool debug = rawEvt.event_counter < config.item_int["n_debug"];
      bool Eureka = rawEvt.findEvtHdr(debug); // Search for the bits that indicate the beginning of an event
      if (!Eureka) {
        cout << "pCTevents event header not found after " << rawEvt.event_counter << " events.\n";
        break;
      }

      if (debug) cout << "Event beginning found " << rawEvt.event_counter << endl;

      /////////////////////////////////////////////////
      // Call the method to unpack the raw event data
      /////////////////////////////////////////////////

      rawEvt.readOneEvent(debug);
      Calibrate->rawPh(rawEvt); // Accumulating histograms for pedestal analysis

      bool daqErr = rawEvt.bad_fpga_address || rawEvt.tag_mismatch || rawEvt.CRC_error || rawEvt.chip_error ||
                    rawEvt.bad_strip_address;
      if (daqErr)
        nErrDmp++;
      if (debug || (nErrDmp < 10 && daqErr)) {
        if (daqErr) {
          cout << " dumping event with DAQ error: Bad FPGA=" << rawEvt.bad_fpga_address
               << " Tag mismatch=" << rawEvt.tag_mismatch;
          cout << " CRC error=" << rawEvt.CRC_error << " Chip error=" << rawEvt.chip_error
               << " Bad strip=" << rawEvt.bad_strip_address << endl;
        }
        rawEvt.dumpEvt(); // Detailed print-out of the raw data
      }
      if (rawEvt.DAQ_error) { // Don't try to reconstruct error events
        // Decide whether to read another event
        cout << "Preprocessing skipping reconstruction at event count "
             << rawEvt.event_counter << " due to DAQ error\n";
        rawEvt.doWeStop(config.item_int["max_events"], config.item_int["max_time"]); // This will set the stop_reading flag inside the rawEvt instance

        if (rawEvt.stop_reading) cout << "Preprocessing stopping after " << rawEvt.event_counter << " events.\n";
	
        continue;
      }
      unsigned int timeStampOut = rawEvt.time_tag / 16; // Reduced precision to fit into 32 bits
      unsigned int eventIdOut = rawEvt.event_number;

      // Calculate the stage angle in the case of a continuous scan, assuming a known constant rotation velocity
      float theta;
      if (config.item_int["continuous"]) {
        double theta2 = ((double)(rawEvt.time_tag)) * Geometry->timeRes() * Geometry->stageSpeed();
        theta = (float)theta2;
      } 
      else theta = config.item_float["projection"];
        
      
      /////////////////////////////////////////////////
      // EXTRACT THE TRACKER COORDINATE INFORMATION
      /////////////////////////////////////////////////

      TkrHits pCThits(rawEvt, Geometry, debug);
      if (debug)
        pCThits.dumpHits(rawEvt.event_number);

      /////////////////////////////////////////////////
      // PERFORM THE TRACKER PATTERN RECOGNITION
      /////////////////////////////////////////////////

      pCT_Tracking pCTtracks(pCThits, Geometry);
      if (debug)
        pCTtracks.dumpTracks(rawEvt.event_number);
      if (rawEvt.event_counter < config.item_int["n_plot"])
        pCTtracks.displayEvent(rawEvt.event_number, pCThits, config.item_str["outputDir"]);

      ////////////////////////////////////////////////////////////////////
      // WRITE OUT ONLY EVENTS THAT ARE SUITABLE FOR IMAGE RECONSTRUCTION
      ////////////////////////////////////////////////////////////////////

      if (config.item_int["level"] > 0 && cuts.cutEvt(pCTtracks, pCThits)) {
        if (rawEvt.event_counter % 100000 == 0)
          cout << "Event " << rawEvt.event_number << ", GoodEvtCount " << cuts.nKeep
               << " timeTag " << rawEvt.time_tag << ", theta  " << theta << endl;

        // Write the good events out to a temporary file
        float Vhit[4], Thit[4];
        if (pCTtracks.nTracks == 1) {
          for (int lyr = 0; lyr < 4; lyr++) {
            Vhit[lyr] = pCTtracks.VTracks[pCTtracks.itkV].X[lyr];
            Uhit[lyr] = pCTtracks.VTracks[pCTtracks.itkV].U[lyr];
            if (lyr < 2)
              Thit[lyr] = pCTtracks.frontPredT(pCTtracks.itkT, Uhit[lyr]); // Extrapolate the T tracks to the
                                                                           // U planes occupied by the V
                                                                           // layers
            else
              Thit[lyr] = pCTtracks.backPredT(pCTtracks.itkT, Uhit[lyr]);
            if (debug) cout << "Layer " << lyr << " T,U,V= " << Thit[lyr] << " " << Uhit[lyr] << " " << Vhit[lyr] << endl;
              
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
          if (Uhit[2] > -900. && Uhit[3] > -900. && UhitT[2] > -900. &&
              UhitT[3] > -900.) {                                       // Complete rear tracker vector
            Thit[2] = Geometry->extrap2D(&UhitT[2], &ThitD[2], Uhit[2]); // Displace the T measurements
                                                                        // to same U as V measurements
            Thit[3] = Geometry->extrap2D(&UhitT[2], &ThitD[2], Uhit[3]);
          } else {
            Thit[2] = ThitD[2];
            Thit[3] = ThitD[3];
          }
          if (debug) {
            for (int lyr = 0; lyr < 4; lyr++) {
              cout << "  Layer " << lyr << ": Uv=" << Uhit[lyr] << ", Ut=" << UhitT[lyr] << ", V=" << Vhit[lyr]
                   << ", Ti=" << ThitD[lyr] << ", Tf=" << Thit[lyr] << endl;
            }
          }
        }

        unsigned char mask[5] = { 0x01, 0x02, 0x04, 0x08, 0x10 };
        unsigned char OTR = 0;
        if (rawEvt.enrg_fpga[0].OTR[0])
          OTR = OTR | mask[0];
        if (rawEvt.enrg_fpga[0].OTR[1])
          OTR = OTR | mask[1];
        if (rawEvt.enrg_fpga[0].OTR[2])
          OTR = OTR | mask[2];
        if (rawEvt.enrg_fpga[1].OTR[0])
          OTR = OTR | mask[3];
        if (rawEvt.enrg_fpga[1].OTR[1])
          OTR = OTR | mask[4];
        int phSum[5];
        phSum[0] = rawEvt.enrg_fpga[0].pulse_sum[0];
        phSum[1] = rawEvt.enrg_fpga[0].pulse_sum[1];
        phSum[2] = rawEvt.enrg_fpga[0].pulse_sum[2];
        phSum[3] = rawEvt.enrg_fpga[1].pulse_sum[0];
        phSum[4] = rawEvt.enrg_fpga[1].pulse_sum[1];
        if (debug)
          cout << "Pulse heights= " << phSum[0] << " " << phSum[1] << " " << phSum[2] << " " << phSum[3] << " "
               << phSum[4] << endl;
        memcpy(&outBuff[0], &timeStampOut, sizeof(int));
        memcpy(&outBuff[sizeof(int)], &eventIdOut, sizeof(int));
        memcpy(&outBuff[2 * sizeof(int)], Thit, 4 * sizeof(float));
        memcpy(&outBuff[2 * sizeof(int) + 4 * sizeof(float)], Vhit, 4 * sizeof(float));
        memcpy(&outBuff[2 * sizeof(int) + 8 * sizeof(float)], phSum, 5 * sizeof(int));
        memcpy(&outBuff[8 * sizeof(float) + 7 * sizeof(int)], &OTR, sizeof(unsigned char));
        fwrite(outBuff, sizeof(char), nBuffBytes, fptmp); // Writing the processed event data to a temporary file
      }

      // Decide whether to read another event
      rawEvt.doWeStop(config.item_int["max_events"], config.item_int["max_time"]); // This will set the
      // stop_reading flag
      // inside the rawEvt
      // instance
      if (rawEvt.stop_reading) cout << "Preprocessing stopping after " << rawEvt.event_counter << " events.\n";
             
    }
    catch (const BadEvent &badEvent) {
      // do nothing - go back and find next event header
      cout << badEvent.what() << endl;
      cout << "Preprocessing: bad event data input caught at event number " << rawEvt.event_counter << endl;
    }
  } // End of the loop over events

  cout << endl << " Number of events with DAQ errors reported = " << nErrDmp << endl
       << endl;

  cuts.summary(); // Summarize the event counts

  cout << " V-layer U coordinates = " << Uhit[0] << " " << Uhit[1] << " " << Uhit[2]
       << " " << Uhit[3] << " assumed the same for all events\n";

  if (config.item_int["level"] == 0) {
    nKeep = 0;
    return;
  }
  fclose(fptmp);

  /////////////////////////////////////////////////////////
  // PEDESTAL ANALYSIS FOR THE WEPL DETECTOR
  /////////////////////////////////////////////////////////

  if (rawEvt.event_counter > 90000) {
    Calibrate->getPeds(pedCalibrateROOTFile,config.item_str["inputFileName"].c_str(), rawEvt.run_number, rawEvt.program_version, config.item_float["projection"],
                       cuts.nKeep, rawEvt.start_time);
    cout <<"We are updating the energy detector pedestal settings to the "
      "on-the-fly measurement values.\n";
    for (int stage = 0; stage < 5; stage++) {
      if (Calibrate->Ped[stage] == 0.) {
        cout <<"Keeping the TVcorr pedestal value " << TVcorr->ped[stage]
             << " for stage " << stage << endl;
      } else {
        cout <<"Pedestal from calibration file = " << TVcorr->ped[stage]
             << "     Drift= " << Calibrate->Ped[stage] - TVcorr->ped[stage] << endl;
      }
    }
  }
  else cout<<"Not enough events to recalibrate the pedestals." << endl;
  for (int stage = 0; stage < 5; stage++) {
    cout <<"The energy detector pedestal for stage " << stage << " is " << Calibrate->Ped[stage] << endl;
  }

  /////////////////////////////////////////////////////////
  // Read PROCESSED DATA BACK FROM THE TEMPORARY FILE
  // AND CALCULATE GAIN CORRECTION FACTORS
  /////////////////////////////////////////////////////////

  if (config.item_int["recalibrate"]) {
    fptmp = fopen(tempfile.c_str(), "rb");
    if (fptmp == NULL) {
      cout <<"Failed to reopen the temporary file\n";
      exit(1);
    }
    fseek(fptmp, 0L, SEEK_END);
    size_t file_size = ftell(fptmp);
    rewind(fptmp);
    cout <<"Temporary data file size=" << file_size << endl;

    for (int EvtNum = 0; EvtNum < cuts.nKeep; EvtNum++) {
      if (EvtNum % 100000 == 0) {
        cout << "Processing event " << EvtNum
             << " from the temp file for calibration." << endl;
      }
      int ret = fread(outBuff, sizeof(char), nBuffBytes, fptmp);
      float Vhit[4], Thit[4];
      int phSum[5];
      memcpy(Thit, &outBuff[2 * sizeof(int)], 4 * sizeof(float));
      memcpy(Vhit, &outBuff[2 * sizeof(int) + 4 * sizeof(float)], 4 * sizeof(float));
      memcpy(phSum, &outBuff[2 * sizeof(int) + 8 * sizeof(float)], 5 * sizeof(int));
      if (EvtNum <= config.item_int["n_debug"]) {
        cout << "pCTevents reading back the temporary file for calibrations. . .\n";
        cout << "  Vhit     Uhit    Thit\n";
        for (int lyr = 0; lyr < 4; lyr++) {
          cout << "  " << Vhit[lyr] << "  " << Uhit[lyr] << "  " << Thit[lyr] << endl;
        }

        cout << "  Stage pulse sums= ";
        for (int stage = 0; stage < 5; stage++)
          cout << phSum[stage] << "  ";
        cout << endl;
      }

      // Calculate the calibrated WEPL
      float Ene[5], Vedet[5], Tedet[5];
      double Tback[2], Vback[2];
      double Tfront[2];
      // double Vfront[2];
      for (int lyr = 2; lyr < 4; lyr++) {
        Tback[lyr - 2] = Thit[lyr];
        Vback[lyr - 2] = Vhit[lyr];
      }
      for (int lyr = 0; lyr < 2; lyr++) {
        Tfront[lyr] = Thit[lyr];
        // Vfront[lyr] = Vhit[lyr];
      }
      float Tphantom = Geometry->extrap2D(Uhit, Tfront, 0.); // Center of the stage
      // float Vphantom = Geometry->extrap2D(Uhit, Vfront, 0.);
      for (int stage = 0; stage < 5; stage++) {
        // Extrapolate the rear track vector to the energy detector stage
        Vedet[stage] = Geometry->extrap2D(&Uhit[2], Vback, Geometry->energyDetectorU(stage));
        Tedet[stage] = Geometry->extrap2D(&Uhit[2], Tback, Geometry->energyDetectorU(stage));
        bool inBounds;
        Ene[stage] = ((float)phSum[stage] - Calibrate->Ped[stage]) *
                     TVcorr->corrFactor(stage, Tedet[stage], Vedet[stage], inBounds);
      }

      Calibrate->weplEvt(Vedet[0], Tphantom, Ene); // Accumulate histograms for gain recalibration
    }
    Calibrate->getGains(TVcorr, pedCalibrateROOTFile, config.item_str["inputFileName"].c_str(), rawEvt.run_number, rawEvt.program_version, config.item_float["projection"],
                        cuts.nKeep, rawEvt.start_time);

    fclose(fptmp);
    cout << "closed the temporary file " << tempfile << endl;
  }
  cout <<"Gain correction factors: ";
  for (int stage = 0; stage < 5; stage++)
    cout << Calibrate->corrFac[stage] << "  ";
  cout << endl;

  nKeep = cuts.nKeep;
  cout << "done with pCTevents.  Number of events kept=" << nKeep << endl;

  return;
}

//*********** Driving program for pCT preprocessing **************
int Preprocessing::ProcessFile(float phantomSize, string partType, float wedgeOffset, float fileFraction,
                               int numbTkrFPGA, int numbEdetFPGA, std::string KillCh) {

  cout << "Preprocessing.cpp: Entering the driver routine for pCT preprocessing. . ." << endl;
  cout << "The phantom is assumed to be less than " << phantomSize << " in radius, for recalibration of gains." << endl;
  cout << "The wedge phantom offset is assumed to be " << wedgeOffset << endl;
  cout << "There are " << numbTkrFPGA << " tracker FPGAs and " << numbEdetFPGA << " energy detector FPGAs" << endl;
  cout << "The file with the list of channels to kill is " << KillCh << endl;
  cout << "The file fraction to use is " << fileFraction << endl;
  if (fileFraction > 1.0)     fileFraction = 1.0;

  // Divide the file into pieces, one for each file
  float fractSize = fileFraction * static_cast<float>(file_size);
  size_t sizeToUse = static_cast<size_t>(fractSize);
  cout << "Preprocessing::ProcessFile, file_size=" << file_size << " sizeToUse=" << sizeToUse << endl;
  size_t fileSize = sizeToUse;
  rewind(in_file);

  struct chKill {
    int FPGA;
    int chip;
    int channel;
  };
  vector<chKill> chToKill;

  // Get the list of tracker channels to kill
  if (KillCh != "null") {
    FILE *fkill = fopen(KillCh.c_str(), "r");
    if (fkill != NULL) {
      int FPGA, chip, ch;
      while (fscanf(fkill, "%d %d %d", &FPGA, &chip, &ch) == 3) {
        //cout << "Tracker FPGA " << FPGA << " chip " << chip << " channel " << ch << " will be killed." << endl;
        chKill entry;
        entry.FPGA = FPGA;
        entry.chip = chip;
        entry.channel = ch;
        chToKill.push_back(entry);
      }
    } else {
      cout << "**** Error: could not open the file " << KillCh << " giving the tracker channels to be killed." << endl;
    }
  }

  // Create an instance of the class for parsing and storing the raw data from
  // the input file
  pCTraw rawEvt(in_file, fileSize, 0, numbTkrFPGA, numbEdetFPGA); // For the mother thread only
  for (auto entry : chToKill)
    rawEvt.pCTkillStrip(entry.FPGA, entry.chip, entry.channel);

  rawEvt.readRunHeader(config.item_str["inputFileName"].c_str());    // Look for the run header bits and parse them
  pCTgeo* Geometry = new pCTgeo(0.);   // Create a class instance with all of the geometry information

    

  if (config.item_int["level"] < 0 || config.item_int["level"] > 2)
    config.item_int["level"] = 2;
  switch (config.item_int["level"]) {
  case 0:
    cout << "Preprocessing.cpp: Only raw and tracker data will be analyzed\n";
    break;
  case 1:
    cout << "Preprocessing.cpp: Raw data, tracker data, and WEPL data will be analyzed\n";
    break;
  case 2:
    cout << "Preprocessing.cpp: The full analysis will be executed and projection data output\n";
    break;
  }

  cout << "Preprocessing.cpp: The output directory is " << config.item_str["outputDir"] << endl;

  // Check whether the specified stage angle agrees with what is in the data
  // file
  if (!config.item_int["continuous"]) {
    if (config.item_float["projection"] > -360.0) {
      if (abs(config.item_float["projection"] - rawEvt.stage_angle) / config.item_float["projection"] > 0.001) {
        cout << "Preprocessing.cpp: The provided projection angle does not "
	  "match the input file run header!\n";
        cout << "The provided projection angle = " << config.item_float["projection"] << endl;
        cout << "The stage angle from the file = " << rawEvt.stage_angle << endl;
        cout << "We are overriding the value from the input file run header.\n";
      }
    } else {
      config.item_float["projection"] = rawEvt.stage_angle;
      cout << "Preprocessing.cpp: We are setting the projection angle "
              "according to the input file value of " << config.item_float["projection"] << endl;
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
  } else
    cout << "Preprocessing.cpp: Was not able to parse the run date from " << rawEvt.start_time << endl;

  ///////////////////////////////////////////// ////////////
  // GET EXISTING WEPL CALIBRATION CONSTANTS  [CEO Aug 2016]
  /////////////////////////////////////////////////////////

  if (WcalibFile.size() == 0) {
    WcalibFile = config.item_str["outputDir"] + "/WcalibW.txt";
    cout << "Preprocessing.cpp: Since no WEPL calibration filename and path "
            "was specified, we will look for it in " << WcalibFile << endl;
  }
  Wepl *WEPL = new Wepl(WcalibFile.c_str(), year, month, day, rawEvt.run_number, partType, dodEEFilter, config.item_str["outputDir"]);
  WEPL->SetEthresholds1(StgThr[0], StgThr[1], StgThr[2], StgThr[3], StgThr[4]);
  
  if (TVcorrFile.size() == 0) {
    TVcorrFile = config.item_str["outputDir"] + "/TVcorr.txt"; // Default TVcorr calibration file location.
    cout << "Preprocessing.cpp: Since no TV calibration filename and path was "
            "specified, we will look for it in " << TVcorrFile << endl;
  }
  TVcorrection *TVcorr = new TVcorrection(TVcorrFile.c_str(), year, month, day, rawEvt.run_number);

  // Create a vector of pointers to instances of the pedestal and gain calibration class
    float t1 = -150.; // These define two ranges for finding protons passing
                    // through zero phantom material, for gain calibration
  float t2 = -150.; // ****** Let's keep this to one side only, for now, to
                    // accommodate the wedge calibration runs with bricks
                    //    float t2= -phantomSize + wedgeOffset;
  float t3 = phantomSize + wedgeOffset;
  float t4 = 150.;
  float pedestals[5];
  for (int stage = 0; stage < 5; stage++) pedestals[stage] = TVcorr->ped[stage];


  int pdstlr[5];
  for (int stage = 0; stage < nStage; stage++) pdstlr[stage] = config.item_int[Form("pedrng%d",stage)];
  pedGainCalib* Calibrate = new pedGainCalib(config.item_str["outputDir"], pdstlr, pedestals, 0, t1, t2, t3, t4, partType); 

  /////////////////////////////////////////////////////////////////
  // Call the routine that reads the data and does the analysis.
  /////////////////////////////////////////////////////////////////

  int nKeep;
  // Uhit is filled and returned for use below, just to save space in the temporary file.
  pCTevents(config, Geometry, rawEvt, Calibrate, TVcorr, std::ref(nKeep), Uhit);
  if (config.item_int["level"] < 2) {
    cout << "Preprocessing.cpp: pCT_Preprocessing all done with the raw data "
            "monitoring task." << endl;
    return 0;
  }

  /////////////////////////////////////////////////////////
  // NOW READ PROCESSED DATA BACK FROM THE TEMPORARY FILE
  // AND COMPLETE THE WEPL RECONSTRUCTION
  /////////////////////////////////////////////////////////

  //        [CEO Jan 2016] Buffer to store event data for writing to a temporary file
  //        buff[0:3]          angle
  //        buff[4:19]         v_hitsV[0:3]    T locations of the T-projection hits
  //        buff[20:35]        t_hitsV[0:3]    V locations of the V-projection hits
  //        buff[36:55]        energy[0:4]     Raw ADC counts from WEPL detector
  //        buff[56]           OTR             ADC out of range indicators
  char outBuff[93]; // Make this extra large, just in case floats are not 4
                    // bytes in size
  int nBuffBytes = 8 * sizeof(float) + 7 * sizeof(int) + sizeof(unsigned char);
  cout << "Preprocessing.cpp: Buffer size for writing the temporary file is " << nBuffBytes << " bytes\n";
  memset(outBuff, 0, nBuffBytes);

  // Vectors of vectors to hold the proton histories pending writing out the
  // projection files
  vector<vector<float> > V0(fileBins);
  vector<vector<float> > V1(fileBins);
  vector<vector<float> > V2(fileBins);
  vector<vector<float> > V3(fileBins);
  vector<vector<float> > T0(fileBins);
  vector<vector<float> > T1(fileBins);
  vector<vector<float> > T2(fileBins);
  vector<vector<float> > T3(fileBins);
  int sizeE = 1;
  if (energyOutput) sizeE = fileBins;
    
  vector<vector<float> > E1(sizeE);
  vector<vector<float> > E2(sizeE);
  vector<vector<float> > E3(sizeE);
  vector<vector<float> > E4(sizeE);
  vector<vector<float> > E5(sizeE);
  int sizeT = 1;
  if (timeStampOutput)
    sizeT = fileBins;
  vector<vector<unsigned int> > TS(sizeT);
  int sizeI = 1;
  if (eventIDOutput)
    sizeI = fileBins;
  vector<vector<unsigned int> > EventIDs(sizeI);
  vector<vector<float> > WetBinary(fileBins);
  vector<vector<float> > ProjAngle(fileBins);
  int totEvt = 0;
  cout << "Preprocessing.cpp: the total number of events in the temporary "
          "files is " << totEvt << endl;
  totEvt = nKeep;
  int alloc = totEvt / fileBins;
  if (alloc > 10) { // Reserve memory space for the vectors -- this is bonkers, there is no guarantee that fileBins is
                    // the same as angle
    for (int bin = 0; bin < fileBins; ++bin) {
      V0[bin].reserve(alloc);
      V1[bin].reserve(alloc);
      V2[bin].reserve(alloc);
      V3[bin].reserve(alloc);
      T0[bin].reserve(alloc);
      T1[bin].reserve(alloc);
      T2[bin].reserve(alloc);
      T3[bin].reserve(alloc);
      if (energyOutput) {
        E1[bin].reserve(alloc);
        E2[bin].reserve(alloc);
        E3[bin].reserve(alloc);
        E4[bin].reserve(alloc);
        E5[bin].reserve(alloc);
      }
      if (timeStampOutput)
        TS[bin].reserve(alloc);
      if (eventIDOutput)
        EventIDs[bin].reserve(alloc);
      WetBinary[bin].reserve(alloc);
      ProjAngle[bin].reserve(alloc);
    }
  }

  FILE *fptmp;
  string tempfile;
  int nEvtot = 0;
  int nBadWEPL = 0;
  int nBadTimeStamp = 0;
  unsigned int timeStampOld = 0;
  long long timeStampOffset = 0;
  
  tempfile = config.item_str["outputDir"] + "/extracted_data_0d.tmp";
  fptmp = fopen(tempfile.c_str(), "rb");
  if (fptmp == NULL) {
    perror("Preprocessing.cpp: Failed to open the temporary file");
    exit(1);
  }

  
  for (int EvtNum = 0; EvtNum < nKeep; EvtNum++) { // Analyze data from the temporary file event by event
    
    ret = fread(outBuff, sizeof(char), nBuffBytes, fptmp);
    float Vhit[4], Thit[4];
    int phSum[5];
    unsigned int timeStamp;
    unsigned int event_id;
    unsigned char OTR;
    memcpy(&timeStamp, &outBuff[0], sizeof(int));
    memcpy(&event_id, &outBuff[sizeof(int)], sizeof(int));
    memcpy(Thit, &outBuff[2 * sizeof(int)], 4 * sizeof(float));
    memcpy(Vhit, &outBuff[2 * sizeof(int) + 4 * sizeof(float)], 4 * sizeof(float));
    memcpy(phSum, &outBuff[2 * sizeof(int) + 8 * sizeof(float)], 5 * sizeof(int));
    memcpy(&OTR, &outBuff[8 * sizeof(float) + 7 * sizeof(int)], sizeof(unsigned char));
    
    // Check for a flaky time stamp, and watch out for roll-over of the time
    // stamp counter (after about 10 minutes)!
    // Before V65 of the event builder firmware there were frequent overflows
    // of the time stamp buffer, causing decreasing values for short times
    
    if (timeStamp < timeStampOld) {
      nBadTimeStamp++;
      // cout << "ERROR: The time stamp decreased since the previous event:
      // tDiff = " << tDiff << endl;
      if (timeStampOld - timeStamp > 10000000.) {
	cout << "***** Preprocessing: the time stamp decreased a lot since "
	  "the previous event; we will assume that the counter rolled "
	  "over.";
	cout << "  Previous time stamp = " << timeStampOld << "  Time stamp = " << timeStamp
	     << "   Time stamp offset = " << timeStampOffset << endl;
	timeStampOffset = timeStampOffset + pow(2, 36);
      }
    }
    timeStampOld = timeStamp;
    float theta;
    if (config.item_int["continuous"]) {
      long long longTimeStamp = 16 * ((long long)timeStamp);
      double theta2 = ((double)(longTimeStamp + timeStampOffset)) * Geometry->timeRes() * Geometry->stageSpeed();
      theta = (float)theta2 + initialAngle;
    } 
    else  theta = proj_angle;
    bool debug = EvtNum < config.item_int["n_debug"];
    if (debug) {
      cout << EvtNum << " Preprocessing.cpp: File position for temp file 0 = " << ftell(fptmp) << endl;
      cout << EvtNum << " Reading back event data from the temporary file 0"  << endl;
      cout << "   Time resolution = " << Geometry->timeRes() << " stage speed = " << Geometry->stageSpeed() << endl;
      cout << "   Time Stamp = " << (float)timeStamp * 16.0 * Geometry->timeRes() << " seconds" << endl;
      cout << "   Event ID = " << event_id << endl;
      cout << "   Theta = " << theta << endl;
      cout << "  Vhit     Uhit    Thit\n";
      for (int lyr = 0; lyr < 4; lyr++) cout << "  " << Vhit[lyr] << "  " << Uhit[lyr] << "  " << Thit[lyr] << endl;
	
      cout << "  Stage pulse sums= ";
      for (int stage = 0; stage < 5; stage++) cout << phSum[stage] << "  ";
	
      cout << endl;
      if (OTR != 0)
          cout << "   At least one stage has a sample out of range " << endl;
    }
    
    // Calculate the calibrated WEPL
    double Tback[2], Vback[2];
    for (int lyr = 2; lyr < 4; lyr++)
      {
	Tback[lyr - 2] = Thit[lyr];
	Vback[lyr - 2] = Vhit[lyr];
      }
    float Ene[5];
    bool inBounds;
    for (int stage = 0; stage < 5; stage++) {
      // TVcorrection and MeV conversion.  Extrapolate the rear track vector
      // to the energy detector location.
      float Vedet = Geometry->extrap2D(&Uhit[2], Vback, Geometry->energyDetectorU(stage));
      float Tedet = Geometry->extrap2D(&Uhit[2], Tback, Geometry->energyDetectorU(stage));
      float TVCorrFactor = TVcorr->corrFactor(stage, Tedet, Vedet, inBounds);
      Ene[stage] = Calibrate->corrFac[stage] * ((float)phSum[stage] - Calibrate->Ped[stage]) * TVCorrFactor;      
    }
    
    float Wet;
    Wet = WEPL->EtoWEPL(Ene); // Energy to WEPL conversion
    if (Wet < -999. || Wet > 999.)
      ++nBadWEPL;
    if (debug) {
      cout << "  WEPL = " << Wet << endl;
      cout << "  Corrected stage energies= ";
      for (int stage = 0; stage < 5; stage++)
	cout << Ene[stage] << "  ";
      cout << endl;
    }
    // Store the data in vectors according to the angular bin (angle index k)
    // First eliminate events rejected by the WEPL analysis (WEPL = -1000 or
    // WEPL = 1000)
    if (EvtNum % 100000 == 0) {
      cout << " Processing event " << EvtNum
	   << " from the temp file, time stamp=" << timeStamp;
      cout << " angle=" << theta << endl;
    }
    if (config.item_int["level"] == 2) {
      ++nEvtot;
      int k = (EvtNum - EvtNum % alloc) / alloc; // Split the event per files
      if (k == fileBins) k--; // Sanity check
      
      V0[k].push_back(Vhit[0]);
      V1[k].push_back(Vhit[1]);
      V2[k].push_back(Vhit[2]);
      V3[k].push_back(Vhit[3]);
      
      T0[k].push_back(Thit[0]);
      T1[k].push_back(Thit[1]);
      T2[k].push_back(Thit[2]);
      T3[k].push_back(Thit[3]);

      ProjAngle[k].push_back(theta);
      WetBinary[k].push_back(Wet);
      if (energyOutput)
	{
	  E1[k].push_back(Ene[0]);
	  E2[k].push_back(Ene[1]);
	  E3[k].push_back(Ene[2]);
	  E4[k].push_back(Ene[3]);
	  E5[k].push_back(Ene[4]);
	}
      if (timeStampOutput)  TS[k].push_back(timeStamp);
      if (eventIDOutput) EventIDs[k].push_back(event_id);
    }
  }
  delete WEPL;
  fclose(fptmp);
  cout << "Preprocessing.cpp: Try to delete the temporary file " << tempfile << endl;
  string cmd;
  cmd = "rm -f " + tempfile;
  cout << "The command to delete the temporary file is " << cmd << endl;
  ret = system(cmd.c_str());

  cout << endl;
  cout << "Preprocessing.cpp: The total number of events saved for output was " << nEvtot << endl;
  cout << "                   The number of events rejected with bad WEPL was " << nBadWEPL << endl;
  cout << "                   The number of decreasing time stamps was " << nBadTimeStamp << endl;
  cout << "                   The accumulated time-stamp correction is, in 10ns units, " << timeStampOffset << endl
       << endl;

  if(config.item_int["level"] == 1) {
    cout << "Preprocess.cpp: The pCT_Preprocessing monitoring task is all done. . ." << endl;
    
    return 0;
  }

  // [CEO Jan 2016] Save the projection data for each angle bin (single bin in the case it is not a continuous scan)
  for (int k = 0; k < fileBins; k++) {
    float AngleNb;
    int Event_Counter = V0[k].size();
    //sprintf(OutputFilename, "%s/projection_%1d.bin", config.item_str["outputDir"].c_str(), k); // float
    char OutputFilename[512];
    //sprintf(OutputFilename, "%s/%s.root", config.item_str["outputDir"], config.item_str["inputFileName"].substr(7, config.item_str["inputFileName"].size()).c_str() );
    cout << "Preprocessing.cpp: Write binary file for file number " << k <<" with "<< Event_Counter << " histories."
         << " Output Filename : " << OutputFilename << endl;

    /*WriteBinaryFile(timeStampOutput, energyOutput, eventIDOutput, AngleNb, OutputFilename, config.item_str["inputFileName"],
                     study_name.c_str(), rawEvt.study_date, Event_Counter, Uhit, V0[k].data(), V1[k].data(),
                     V2[k].data(), V3[k].data(), T0[k].data(), T1[k].data(), T2[k].data(), T3[k].data(), E1[k].data(),
                     E2[k].data(), E3[k].data(), E4[k].data(), E5[k].data(), WetBinary[k].data(), ProjAngle[k].data(),
                     TS[k].data(), EventIDs[k].data());*/

    WriteRootFile(timeStampOutput, energyOutput, eventIDOutput, k, rawEvt.study_date, Event_Counter, Uhit,
		  V0[k].data(), V1[k].data(), V2[k].data(), V3[k].data(),
		  T0[k].data(), T1[k].data(), T2[k].data(), T3[k].data(),
		  E1[k].data(), E2[k].data(), E3[k].data(), E4[k].data(), E5[k].data(),
		  WetBinary[k].data(), ProjAngle[k].data(), TS[k].data(), EventIDs[k].data());
  }


  time_t end_time = time(NULL);
  now = localtime(&end_time);
  printf("Preprocessing.cpp: Local time and date at end of execution: %s", asctime(now));
  double seconds = difftime(end_time, start_time);
  cout << "Preprocessing.cpp: The total time lapse during execution was " << seconds << " seconds.\n";

  cout << "Preprocessing.cpp: pCT_Preprocessing is all done, including output "
          "of the projection data." << endl;

  delete Calibrate;
  return 0;
}
