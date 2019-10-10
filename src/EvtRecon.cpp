// Process the raw data to produce a simple track list; used in the calibration
// tasks

#include "EvtRecon.h"
#include "BadEvent.h"

void EvtRecon::writeTmp(Event &evt) {
  memcpy(&outBuff[0], evt.Thit, 4 * sizeof(float));
  memcpy(&outBuff[4 * sizeof(float)], evt.Vhit, 4 * sizeof(float));
  memcpy(&outBuff[8 * sizeof(float)], evt.ADC, 5 * sizeof(int));
  fwrite(outBuff, sizeof(char), nBuffBytes, evtFile);
}
void EvtRecon::readTmp(Event &evt) {
  fread(outBuff, sizeof(char), nBuffBytes, evtFile);
  memcpy(evt.Thit, &outBuff[0], 4 * sizeof(float));
  memcpy(evt.Vhit, &outBuff[4 * sizeof(float)], 4 * sizeof(float));
  memcpy(evt.ADC, &outBuff[8 * sizeof(float)], 5 * sizeof(int));
}
void EvtRecon::dumpTmp(Event evt) {
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
void EvtRecon::reopenTmpFile() {
  evtFile = fopen(evtFileName.c_str(), "rb");
  if (evtFile == NULL) {
    cout << "EvtRecon_" << Nblk << ":  Failed to reopen the temporary file" << evtFileName << endl;
    exit(1);
  }
  fseek(evtFile, 0L, SEEK_END);
  size_t file_size = ftell(evtFile);
  rewind(evtFile);
  cout << "EvtRecon_" << Nblk << ": temporary data file " << evtFileName << " size=" << file_size << endl;
}
void EvtRecon::rewindTmpFile() {
  rewind(evtFile);
  fseek(evtFile, 0L, SEEK_END);
  size_t file_size = ftell(evtFile);
  rewind(evtFile);
  cout << "EvtRecon_" << Nblk << ": temporary data file " << evtFileName << " size=" << file_size << endl;
}
void EvtRecon::delTmpFile() {
  fclose(evtFile);
  string cmd;
  if (OsName == "Windows") {
    cmd = "del " + evtFileName;
  } else {
    cmd = "rm -f " + evtFileName;
  }
  system(cmd.c_str()); // Delete the temporary file
}

EvtRecon::EvtRecon(pCTgeo &Geometry, TVcorrection *const TVcorr, string inputFileName, string OutputDir, int max_events,
                   int max_time, int n_debug, int n_plot, int nBlocks, bool useTemp, bool doGains, string partType,
                   int pdstlr[5], bool reCalibrate, std::string OsName) {

  cout << "*********** Entering EvtRecon " << nBlocks << " for processing raw calibration data **************" << endl;
  start_time = time(NULL);
  now = localtime(&start_time);
  printf("Current local time and date: %s", asctime(now));

  if (max_events > 0)
    cout << "The preprocessing will halt after processing " << max_events << " events\n";
  cout << "The input data are assumed to be real, not simulated" << endl;

  this->OsName = OsName;

  Nblk = nBlocks;
  gainAnalysis = doGains;
  useTmpFile = useTemp;
  nBuffBytes = 8 * sizeof(float) + 5 * sizeof(int);
  memset(outBuff, 0, nBuffBytes);
  evtFileName = OutputDir + "/EvtTmpFile" + to_str(nBlocks) + ".dat";
  if (useTmpFile) {
    evtFile = fopen(evtFileName.c_str(), "wb");
    if (evtFile == NULL) {
      cout << "EvtRecon: failed to open the temporary file " << evtFileName << endl;
      exit(1);
    }
    cout << "EvtRecon_" << nBlocks << ": Opened file " << evtFileName << " for temporary storage of event data."
         << endl;
  }
  nEvents = 0;

  // ******************** Open the input data file
  // *******************************

  cout << "EvtRecon_" << nBlocks << ": Reading the input raw data file " << inputFileName << endl;

  in_file = fopen(inputFileName.c_str(), "rb");
  if (in_file == NULL) {
    cout << "EvtRecon: error opening the input raw data file " << inputFileName << " for nBlocks= " << nBlocks << endl;
    string msg = "EvtRecon: Error opening the input raw data file " + inputFileName + " for nBlocks=" + to_str(nBlocks);
    cout << "EvtRecon is aborting the run" << endl;
    perror(msg.c_str());
    exit(1);
  }

  fseek(in_file, 0L, SEEK_END);
  file_size = ftell(in_file);
  rewind(in_file);
  cout << "EvtRecon_" << nBlocks << ": Input raw data file size=" << file_size << endl;

  pCTcut cuts(nBlocks, 1, 2, 2); // Initialize the code for event selection

  // Create an instance of the class for parsing and storing the raw data from
  // the input file
  pCTraw rawEvt(in_file, file_size, 0, num_tkr_fpga, num_enrg_fpga);

  rawEvt.readRunHeader(inputFileName.c_str()); // Look for the run header bits and parse them
  runNumber = rawEvt.run_number;
  runStartTime = rawEvt.start_time;
  study_date = rawEvt.study_date;
  stage_angle = rawEvt.stage_angle;
  program_version = rawEvt.program_version;

  // Allocate memory for the output event data (for efficiency---it can get more
  // memory later if needed)
  int esL = ((int)file_size) / 50;
  unsigned int estLength = (esL < max_events) ? esL : max_events;
  if (!useTmpFile) {
    cout << "Estimated number of bytes needed for temporary storage = " << estLength << endl;
    evtList.reserve(estLength);
  }

  Event tmpEvt;

  // Set the range in T that should be occupied by unobstructed protons. Only
  // the +T side!
  // September 12, 2018 we started moving the edge of the bricks 2 cm past the
  // wedge,
  // so this range had to be moved.  It should be okay also for earlier runs
  // with the wedge phantom.
  float wedgeLimit = Geometry.getTWedgeBreaks(4) + 5.0; // NO BRICK OFFSET //+5.
  float openRange = 20.0;
  cout << "EvtRecon_ " << nBlocks << ": range in T used for recalibration is " << wedgeLimit << " to "
       << wedgeLimit + openRange << " mm\n";
  float pedestals[nStage];
  for (int stage = 0; stage < nStage; stage++)
    pedestals[stage] = 0.;
  pedGainCalib Calibrate(OutputDir, pdstlr, pedestals, nBlocks, -150., -151., wedgeLimit, wedgeLimit + openRange,
                         partType, OsName); // Specifies range in t where particles miss the wedge phantom

  // Here is the start of the event loop
  while (!rawEvt.stop_reading) {
    try {
      bool Eureka = rawEvt.findEvtHdr(false); // Search for the bits that indicate the beginning of an event
      if (!Eureka) {
        cout << "EvtRecon_" << nBlocks << ": event header not found after " << rawEvt.event_counter << " events.\n";
        break;
      }

      if (rawEvt.event_counter % 1000000 == 0)
        cout << "EvtRecon_" << nBlocks << ": Processing raw event " << rawEvt.event_counter << endl;
      bool debug = rawEvt.event_counter < n_debug;
      if (debug)
        cout << "EvtRecon_" << nBlocks << ": Event beginning found for event count " << rawEvt.event_counter << endl;

      rawEvt.readOneEvent(debug);
      if (debug)
        rawEvt.dumpEvt(); // Detailed print-out of the raw data

      Calibrate.rawPh(rawEvt); // Accumulate pedestal histograms

      TkrHits pCThits(rawEvt, Geometry, debug); // Reconstruct the tracker hits
                                                // from the raw strip data
      if (debug)
        pCThits.dumpHits(rawEvt.event_number);

      pCT_Tracking pCTtracks(pCThits, Geometry); // Track pattern recognition
      if (debug)
        pCTtracks.dumpTracks(rawEvt.event_number);
      if (rawEvt.event_counter < n_plot)
        pCTtracks.displayEvent(rawEvt.event_number, pCThits, OutputDir);

      // Store the reconstructed track and raw energy information into the
      // output list
      if (cuts.cutEvt(false, pCTtracks, pCThits)) {
        for (int lyr = 0; lyr < 4; lyr++) {
          tmpEvt.Thit[lyr] = pCTtracks.TTracks[pCTtracks.itkT].X[lyr];
          uhitT[lyr] = pCTtracks.TTracks[pCTtracks.itkT].U[lyr];
          tmpEvt.Vhit[lyr] = pCTtracks.VTracks[pCTtracks.itkV].X[lyr];
          uhitV[lyr] = pCTtracks.VTracks[pCTtracks.itkV].U[lyr];
        }
        tmpEvt.ADC[0] = rawEvt.enrg_fpga[0].pulse_sum[0];
        tmpEvt.ADC[1] = rawEvt.enrg_fpga[0].pulse_sum[1];
        tmpEvt.ADC[2] = rawEvt.enrg_fpga[0].pulse_sum[2];
        tmpEvt.ADC[3] = rawEvt.enrg_fpga[1].pulse_sum[0];
        tmpEvt.ADC[4] = rawEvt.enrg_fpga[1].pulse_sum[1];
        if (useTmpFile)
          writeTmp(tmpEvt);
        else
          evtList.push_back(tmpEvt);
        nEvents++;
      }

      // Decide whether to read another event
      rawEvt.doWeStop(max_events, max_time); // This will set the stop_reading
                                             // flag inside the rawEvt instance
    }
    catch (const BadEvent &badEvent) {
      // do nothing - go back and find next event header
      cout << badEvent.what() << endl;
      cout << "EvtRecon " << nBlocks << ": bad event data input caught at event number " << rawEvt.event_counter
           << endl;
    }
  } // End of the loop over events

  // evtList.shrink_to_fit();  // Exists only in c++11 onward

  cuts.summary(); // Summarize the event counts

  Calibrate.getPeds(inputFileName.c_str(), runNumber, program_version, stage_angle, nEvents, runStartTime);
  for (int stage = 0; stage < 5; stage++) {
    Peds[stage] = Calibrate.newPed(stage);
  }

  if (useTmpFile)
    fclose(evtFile);

  if (gainAnalysis) {
    cout << "Starting the gain analysis for nBlocks= " << nBlocks << endl;
    if (useTmpFile) {
      evtFile = fopen(evtFileName.c_str(), "rb");
      if (evtFile == NULL) {
        cout << "EvtRecon nBlocks=" << nBlocks << ":  Failed to reopen the temporary file\n";
        exit(1);
      }
      fseek(evtFile, 0L, SEEK_END);
      size_t file_size = ftell(evtFile);
      rewind(evtFile);
      cout << "EvtRecon nBlocks=" << nBlocks << ": Temporary data file size=" << file_size << endl;
    }
    for (int EvtNum = 0; EvtNum < nEvents; EvtNum++) {
      if (EvtNum % 1000000 == 0)
        cout << "EvtRecon_" << nBlocks << " gain analysis, processing event " << EvtNum << endl;
      Event thisEvent;
      if (useTmpFile)
        readTmp(thisEvent);
      else
        thisEvent = evtList[EvtNum];
      if (EvtNum <= n_debug) {
        cout << "EvtRecon gain analysis nBlocks=" << nBlocks
             << ", reading back the temporary file for gain calibrations. . .\n";
        cout << "EvtRecon__" << nBlocks << "  Vhit     Thit\n";
        for (int lyr = 0; lyr < 4; lyr++) {
          cout << "  " << thisEvent.Vhit[lyr] << "  " << thisEvent.Thit[lyr] << endl;
        }
        cout << "EvtRecon_" << nBlocks << "  Stage pulse sums= ";
        for (int stage = 0; stage < 5; stage++)
          cout << thisEvent.ADC[stage] << "  ";
        cout << endl;
      }

      // Calculate the calibrated WEPL
      float Ene[5], Vedet[5], Tedet[5];
      double UV[4], UT[4], V[4], T[4];
      for (int i = 0; i < 4; ++i) {
        UV[i] = uhitV[i];
        UT[i] = uhitT[i];
        V[i] = thisEvent.Vhit[i];
        T[i] = thisEvent.Thit[i];
      }
      int nGood = 0;
      for (int stage = 0; stage < 5; stage++) {
        // Extrapolate the rear track vector to the energy detector stage
        Vedet[stage] = Geometry.extrap2D(&UV[2], &V[2], Geometry.energyDetectorU(stage));
        Tedet[stage] = Geometry.extrap2D(&UT[2], &T[2], Geometry.energyDetectorU(stage));

        bool inBounds;
        Ene[stage] = ((float)thisEvent.ADC[stage] - Peds[stage]) *
                     TVcorr->corrFactor(stage, Tedet[stage], Vedet[stage], inBounds);
        if (inBounds)
          nGood++;
      }
      float vPh = Geometry.extrap2D(UV, V, -76.2);
      float tPh = Geometry.extrap2D(UT, T, -76.2);
      if (nGood == 5)
        Calibrate.weplEvt(vPh, tPh, Ene);
    }

    if (reCalibrate) {
      Calibrate.getGains(TVcorr, inputFileName.c_str(), rawEvt.run_number, rawEvt.program_version, 0., cuts.nKeep,
                         rawEvt.start_time);
      cout << "EvtRecon nBricks=" << nBlocks << ": Gain correction factors: ";
      for (int stage = 0; stage < 5; stage++) {
        cout << Calibrate.corrFac[stage] << "  ";
        CorFacs[stage] = Calibrate.corrFac[stage];
      }
      cout << endl;
    } else {
      cout << "EvtRecon nBricks=" << nBlocks << ": Recalibration turned off, Gain correction factors are ";
      for (int stage = 0; stage < 5; stage++) {
        CorFacs[stage] = 1.0;
        cout << CorFacs[stage] << "  ";
      }
      cout << endl;
    }

    fclose(evtFile);
  }

  cout << "EvtRecon_" << nBlocks << " is complete.  " << nEvents << " events were stored in the event list " << endl;

}; // end of the EvtRecon constructor