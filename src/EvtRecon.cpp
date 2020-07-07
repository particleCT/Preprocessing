// Process the raw data to produce a simple track list; used in the calibration tasks
#include "EvtRecon.h"
#include "BadEvent.h"
#include "pCTcut.h"
// Constructors/Destructor
EvtRecon::EvtRecon(){
  theConfig = pCTconfig::GetInstance();
}
EvtRecon::~EvtRecon() {}

////////////////////////////////////////////////////////////////////////////////////
// Read a file and fill the different temporary files needed
////////////////////////////////////////////////////////////////////////////////////
void EvtRecon::ReadInputFile(pCTgeo* Geometry, TVcorrection *const TVcorr , string inputFileName, pedGainCalib* Calibrate){
  // Restart anew
  evtList.clear();

  cout << "*********** Entering EvtRecon " << theConfig->item_int["Nbricks"] << " for processing raw calibration data **************" << endl;
  if (theConfig->item_int["max_events"] > 0) cout << "The preprocessing will halt after processing " << theConfig->item_int["max_events"] << " events\n";

  // Read the input file
  nEvents = 0;  
  cout << "EvtRecon_" << theConfig->item_int["Nbricks"] << ": Reading the input raw data file " << inputFileName << endl;
  in_file = fopen(inputFileName.c_str(), "rb");
  if (in_file == NULL) {
    cout << "EvtRecon: error opening the input raw data file " << inputFileName << " for nBlocks= " << theConfig->item_int["Nbricks"] << endl;
    string msg = "EvtRecon: Error opening the input raw data file " + inputFileName + " for nBlocks=" + to_str(theConfig->item_int["Nbricks"]);
    cout << "EvtRecon is aborting the run" << endl;
    perror(msg.c_str());
    exit(1);
  }

  fseek(in_file, 0L, SEEK_END);
  file_size = ftell(in_file);
  rewind(in_file);
  
  cout << "EvtRecon_" << theConfig->item_int["Nbricks"] << ": Input raw data file size=" << file_size << endl;
  pCTcut cuts;
  // Create an instance of the class for parsing and storing the raw data from the input file
  pCTraw rawEvt(in_file, file_size, 0, num_tkr_fpga, num_enrg_fpga);

  rawEvt.readRunHeader(inputFileName.c_str()); // Look for the run header bits and parse them
  runNumber = rawEvt.run_number;
  runStartTime = rawEvt.start_time;
  study_date = rawEvt.study_date;
  stage_angle = rawEvt.stage_angle;
  program_version = rawEvt.program_version;

  Event tmpEvt;
  // Set the range in T that should be occupied by unobstructed protons. Only the +T side!
  // September 12, 2018 we started moving the edge of the bricks 2 cm past the wedge,
  // so this range had to be moved.  It should be okay also for earlier runs with the wedge phantom.

  // Here is the start of the event loop
  while (!rawEvt.stop_reading) {
    try {
      bool Eureka = rawEvt.findEvtHdr(false); // Search for the bits that indicate the beginning of an event
      if (!Eureka) {
        cout << "EvtRecon_" << theConfig->item_int["Nbricks"] << ": event header not found after " << rawEvt.event_counter << " events.\n";
        break;
      }
      if (rawEvt.event_counter % 1000000 == 0) cout << "EvtRecon_" << theConfig->item_int["Nbricks"] << ": Processing raw event " << rawEvt.event_counter << endl;
        
      rawEvt.readOneEvent(false);
      Calibrate->FillPeds(rawEvt); // Accumulate pedestal histograms

      TkrHits pCThits(rawEvt, Geometry, false); // Reconstruct the tracker hits from the raw strip data
      pCT_Tracking pCTtracks(pCThits, Geometry); // Track pattern recognition

      //Store the reconstructed track and raw energy information into the output list
      if (cuts.cutEvt(pCTtracks, pCThits)) {
        for (int lyr = 0; lyr < 4; lyr++) {
          tmpEvt.Thit[lyr] = pCTtracks.TTracks[pCTtracks.itkT].X[lyr];
          uhitT[lyr]       = pCTtracks.TTracks[pCTtracks.itkT].U[lyr];
          tmpEvt.Vhit[lyr] = pCTtracks.VTracks[pCTtracks.itkV].X[lyr];
          uhitV[lyr]       = pCTtracks.VTracks[pCTtracks.itkV].U[lyr];
        }
        tmpEvt.ADC[0] = rawEvt.enrg_fpga[0].pulse_sum[0];
        tmpEvt.ADC[1] = rawEvt.enrg_fpga[0].pulse_sum[1];
        tmpEvt.ADC[2] = rawEvt.enrg_fpga[0].pulse_sum[2];
        tmpEvt.ADC[3] = rawEvt.enrg_fpga[1].pulse_sum[0];
        tmpEvt.ADC[4] = rawEvt.enrg_fpga[1].pulse_sum[1];

	evtList.push_back(tmpEvt);
        nEvents++;
      }
      // Decide whether to read another event
      rawEvt.doWeStop(theConfig->item_int["max_events"], theConfig->item_int["max_time"]); // This will set the stop_reading flag inside the rawEvt instance
    }
    catch (const BadEvent &badEvent) {
      // do nothing - go back and find next event header
      cout << badEvent.what() << endl;
      cout << "EvtRecon " << theConfig->item_int["Nbricks"] << ": bad event data input caught at event number " << rawEvt.event_counter << endl;
    }
  } // End of the loop over events

  cuts.summary(); // Summarize the event counts
  Calibrate->GetPeds();
  
  for (int stage = 0; stage < 5; stage++) Peds[stage] = Calibrate->Ped[stage];  
  if (theConfig->item_int["doGains"]) {
    cout << "Starting the gain analysis for nBlocks= " << theConfig->item_int["Nbricks"] << endl;
    for (int EvtNum = 0; EvtNum < nEvents; EvtNum++) {
      if (EvtNum % 1000000 == 0) cout << "EvtRecon_" << theConfig->item_int["Nbricks"] << " gain analysis, processing event " << EvtNum << endl;	
      Event thisEvent;
      thisEvent = evtList[EvtNum];

      // Calculate the calibrated WEPL
      float Ene[5], Vedet[5], Tedet[5];
      float UV[4], UT[4], V[4], T[4];
      for (int i = 0; i < 4; ++i) {
        UV[i]  = uhitV[i];
        UT[i]  = uhitT[i];
        V[i]  =  thisEvent.Vhit[i];
        T[i]  =  thisEvent.Thit[i];
      }
      int nGood = 0;
      for (int stage = 0; stage < 5; stage++) {
        // Extrapolate the rear track vector to the energy detector stage
        Vedet[stage] = Geometry->extrap2D(&UV[2], &V[2], Geometry->energyDetectorU(stage));
        Tedet[stage] = Geometry->extrap2D(&UT[2], &T[2], Geometry->energyDetectorU(stage));

        bool inBounds;
        Ene[stage] = ((float)thisEvent.ADC[stage] - Calibrate->Ped[stage]) * TVcorr->corrFactor(stage, Tedet[stage], Vedet[stage], inBounds);
        if (inBounds) nGood++;

      }
      // set it at the entrance brick position
      float vPh = Geometry->extrap2D(UV, V, -76.2); 
      float tPh = Geometry->extrap2D(UT, T, -76.2);
      if (nGood == 5) Calibrate->FillGains(vPh, tPh, Ene, thisEvent.ADC);
    } // end of the event loop

    if (theConfig->item_int["recalibrate"]) {
      Calibrate->GetGains(TVcorr);

      cout << "EvtRecon nBricks=" << theConfig->item_int["Nbricks"] << ": Gain correction factors: ";
      for (int stage = 0; stage < 5; stage++) {
        cout << Calibrate->GainFac[stage] << "  ";
        GainFac[stage] = Calibrate->GainFac[stage];
      }
      cout << endl;
    }
    
    else {
      cout << "EvtRecon nBricks=" << theConfig->item_int["Nbricks"] << ": Recalibration turned off, Gain correction factors are ";
      for (int stage = 0; stage < 5; stage++) {
        GainFac[stage] = 1.0;
        cout << GainFac[stage] << "  ";
      }
      cout << endl;
    }
  }
  cout << "EvtRecon_" << theConfig->item_int["Nbricks"] << " is complete.  " << nEvents << " events were stored in the event list " << endl;

} 
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////             
void EvtRecon::dumpEvt(Event evt) {
  cout << "EvtRecon: event data from the temporary file: " << endl;
  for (int i = 0; i < 4; ++i) {
    cout << "  Layer " << i << " T=" << evt.Thit[i] << ",  V=" << evt.Vhit[i] << endl;
  }
  cout << "  Stage ADC values: ";
  for (int stage = 0; stage < 5; ++stage) {
    cout << evt.ADC[stage] << " ";
  }
  cout << endl;
}
