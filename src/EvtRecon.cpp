// Process the raw data to produce a simple track list; used in the calibration tasks
#include "EvtRecon.h"
#include "BadEvent.h"

// Constructors/Destructor
EvtRecon::EvtRecon(pCTconfig conf, TFile* root): config(conf), rootfile(root)
{}
EvtRecon::~EvtRecon() {}

////////////////////////////////////////////////////////////////////////////////////
// Read a file and fill the different temporary files needed
////////////////////////////////////////////////////////////////////////////////////
void EvtRecon::ReadInputFile(pCTgeo* Geometry, TVcorrection *const TVcorr , string inputFileName){
  
  cout << "*********** Entering EvtRecon " << config.item_int["Nbricks"] << " for processing raw calibration data **************" << endl;
  if (config.item_int["max_events"] > 0) cout << "The preprocessing will halt after processing " << config.item_int["max_events"] << " events\n";

  // Create the Temporary files.
  nEvents = 0;  
  cout << "EvtRecon_" << config.item_int["Nbricks"] << ": Reading the input raw data file " << inputFileName << endl;
  in_file = fopen(inputFileName.c_str(), "rb");
  if (in_file == NULL) {
    cout << "EvtRecon: error opening the input raw data file " << inputFileName << " for nBlocks= " << config.item_int["Nbricks"] << endl;
    string msg = "EvtRecon: Error opening the input raw data file " + inputFileName + " for nBlocks=" + to_str(config.item_int["Nbricks"]);
    cout << "EvtRecon is aborting the run" << endl;
    perror(msg.c_str());
    exit(1);
  }

  fseek(in_file, 0L, SEEK_END);
  file_size = ftell(in_file);
  rewind(in_file);
  
  cout << "EvtRecon_" << config.item_int["Nbricks"] << ": Input raw data file size=" << file_size << endl;
  pCTcut cuts(config);//.item_int["Nbricks"], 1, 2, 2); // Initialize the code for event selection
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

  float wedgeLimit = Geometry->getTWedgeBreaks(4) + 25.0; // NO BRICK OFFSET //+5.
  float openRange = 20.0;
  cout << "EvtRecon_ " << config.item_int["Nbricks"] << ": range in T used for recalibration is " << wedgeLimit << " to " << wedgeLimit + openRange << " mm\n";
  float pedestals[nStage];
  int pdstlr[5];

  for (int stage = 0; stage < nStage; stage++) pedestals[stage] = 0.;
  for (int stage = 0; stage < nStage; stage++) pdstlr[stage] = config.item_int[Form("pedrng%d",stage)];

  pedGainCalib Calibrate(config.item_str["outputDir"], pdstlr, pedestals, config.item_int["Nbricks"], -150., -151., wedgeLimit, wedgeLimit + openRange,
                         config.item_str["partType"]); // Specifies range in t where particles miss the wedge phantom


  // Here is the start of the event loop
  while (!rawEvt.stop_reading) {
    try {
      bool Eureka = rawEvt.findEvtHdr(false); // Search for the bits that indicate the beginning of an event
      if (!Eureka) {
        cout << "EvtRecon_" << config.item_int["Nbricks"] << ": event header not found after " << rawEvt.event_counter << " events.\n";
        break;
      }

      if (rawEvt.event_counter % 1000000 == 0) cout << "EvtRecon_" << config.item_int["Nbricks"] << ": Processing raw event " << rawEvt.event_counter << endl;
        
      bool debug = false;//rawEvt.event_counter < config.item_int["n_debug"];
      //if (debug) cout << "EvtRecon_" << config.item_int["Nbricks"] << ": Event beginning found for event count " << rawEvt.event_counter << endl;
        
      rawEvt.readOneEvent(debug);
      
      //if (debug) rawEvt.dumpEvt(); // Detailed print-out of the raw data

      Calibrate.rawPh(rawEvt); // Accumulate pedestal histograms

      TkrHits pCThits(rawEvt, Geometry, debug); // Reconstruct the tracker hits from the raw strip data
      //if (debug) pCThits.dumpHits(rawEvt.event_number);

      pCT_Tracking pCTtracks(pCThits, Geometry); // Track pattern recognition
      //if (debug) pCTtracks.dumpTracks(rawEvt.event_number);
      //if (rawEvt.event_counter < config.item_int["n_plot"]) pCTtracks.displayEvent(rawEvt.event_number, pCThits, config.item_str["outputDir"]);
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
      rawEvt.doWeStop(config.item_int["max_events"], config.item_int["max_time"]); // This will set the stop_reading flag inside the rawEvt instance
    }
    catch (const BadEvent &badEvent) {
      // do nothing - go back and find next event header
      cout << badEvent.what() << endl;
      cout << "EvtRecon " << config.item_int["Nbricks"] << ": bad event data input caught at event number " << rawEvt.event_counter
           << endl;
    }
  } // End of the loop over events

  cuts.summary(); // Summarize the event counts
  Calibrate.getPeds(rootfile, inputFileName.c_str(), runNumber, program_version, stage_angle, nEvents, runStartTime);
  
  for (int stage = 0; stage < 5; stage++) Peds[stage] = Calibrate.Ped[stage];  
  if (config.item_int["doGains"]) {
    cout << "Starting the gain analysis for nBlocks= " << config.item_int["Nbricks"] << endl;
    for (int EvtNum = 0; EvtNum < nEvents; EvtNum++) {
      if (EvtNum % 1000000 == 0) cout << "EvtRecon_" << config.item_int["Nbricks"] << " gain analysis, processing event " << EvtNum << endl;	
      Event thisEvent;
      thisEvent = evtList[EvtNum];

      // Calculate the calibrated WEPL
      float Ene[5], Vedet[5], Tedet[5];
      double UV[4], UT[4], V[4], T[4];
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
        Ene[stage] = ((float)thisEvent.ADC[stage] - Peds[stage]) * TVcorr->corrFactor(stage, Tedet[stage], Vedet[stage], inBounds);
        if (inBounds) nGood++;

      }
      float vPh = Geometry->extrap2D(UV, V, -76.2);
      float tPh = Geometry->extrap2D(UT, T, -76.2);
      if (nGood == 5) Calibrate.weplEvt(vPh, tPh, Ene);
    } // end of the event loop

    if (config.item_int["recalibrate"]) {
      Calibrate.getGains(TVcorr, rootfile, inputFileName.c_str(), rawEvt.run_number, rawEvt.program_version, 0., cuts.nKeep, rawEvt.start_time);

      cout << "EvtRecon nBricks=" << config.item_int["Nbricks"] << ": Gain correction factors: ";
      for (int stage = 0; stage < 5; stage++) {
        cout << Calibrate.corrFac[stage] << "  ";
        CorFacs[stage] = Calibrate.corrFac[stage];
      }
      cout << endl;
    }
    
    else {
      cout << "EvtRecon nBricks=" << config.item_int["Nbricks"] << ": Recalibration turned off, Gain correction factors are ";
      for (int stage = 0; stage < 5; stage++) {
        CorFacs[stage] = 1.0;
        cout << CorFacs[stage] << "  ";
      }
      cout << endl;
    }
  }
  cout << "EvtRecon_" << config.item_int["Nbricks"] << " is complete.  " << nEvents << " events were stored in the event list " << endl;
  evtList.clear();
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
