// Encapsulate all of Vladimir Bashkirov pedestal and gain on-the-fly
// calibration
// R.P. Johnson   5/22/2016

#include "pedGainCalib.h"

pedGainCalib::pedGainCalib(string Outputdir, int pdstlr[5], float oldPed[5],
                           int thread, float t1, float t2, float t3, float t4,
                           string partType, string OsName) {

  // Two ranges in t occupied by unobstructed (empty) protons
  emtrng1 = t1; // negative t side
  emtrng2 = t2;
  emtrng3 = t3; // positive t side
  emtrng4 = t4;
  this->partType = partType;

  if (OsName == "Windows")
    setTerm = "set terminal wxt size 1300, 600\n";
  else
    setTerm = "set terminal x11 size 1300, 600\n";

  for (int stage = 0; stage < 5; stage++) {
    Ped[stage] = oldPed[stage];
    corrFac[stage] = 1.0;
    cout << "pedGainCalib " << thread
         << ": setting the default pedestal for stage " << stage << " to "
         << Ped[stage] << endl;
  }

  // Define histograms for pedestal calculation
  hPed[0] = new Histogram(400, pdstlr[0], 5., "Pedestal region for stage 0",
                          "ADC", "events");
  hPed[1] = new Histogram(400, pdstlr[1], 5., "Pedestal region for stage 1",
                          "ADC", "events");
  hPed[2] = new Histogram(400, pdstlr[2], 5., "Pedestal region for stage 2",
                          "ADC", "events");
  hPed[3] = new Histogram(400, pdstlr[3], 5., "Pedestal region for stage 3",
                          "ADC", "events");
  hPed[4] = new Histogram(400, pdstlr[4], 5., "Pedestal region for stage 4",
                          "ADC", "events");

  // Define histograms for gain calibration
  if (partType == "H") {
    hEnrg[0] = new Histogram(400, 15., 0.175, "Energy distribution for stage 0",
                             "Energy (MeV)", "proton events");
    hEnrg[1] = new Histogram(400, 15., 0.175, "Energy distribution for stage 1",
                             "Energy (MeV)", "proton events");
    hEnrg[2] = new Histogram(400, 15., 0.175, "Energy distribution for stage 2",
                             "Energy (MeV)", "proton events");
    hEnrg[3] = new Histogram(400, 15., 0.175, "Energy distribution for stage 3",
                             "Energy (MeV)", "proton events");
    hEnrg[4] = new Histogram(400, 25., 0.175, "Energy distribution for stage 4",
                             "Energy (MeV)", "proton events");
    hEnrgTot = new Histogram(800, 0., 0.3, "Sum of stage energies",
                             "Energy (MeV)", "proton events");
  } else {
    hEnrg[0] = new Histogram(400, 60., 0.7, "Energy distribution for stage 0",
                             "Energy (MeV)", "He events");
    hEnrg[1] = new Histogram(400, 60., 0.7, "Energy distribution for stage 1",
                             "Energy (MeV)", "He events");
    hEnrg[2] = new Histogram(400, 60., 0.7, "Energy distribution for stage 2",
                             "Energy (MeV)", "He events");
    hEnrg[3] = new Histogram(400, 60., 0.7, "Energy distribution for stage 3",
                             "Energy (MeV)", "He events");
    hEnrg[4] = new Histogram(400, 60., 0.7, "Energy distribution for stage 4",
                             "Energy (MeV)", "He events");
    hEnrgTot = new Histogram(800, 0., 1.2, "Sum of stage energies",
                             "Energy (MeV)", "He events");
  }

  // Profile plot to make sure that the phantom does not extend into the regions
  // used for gain calibration
  hProfT = new ProfilePlot(100, -150., 3.0, "Stage 0 energy profile in T",
                           "T (mm)", "mean energy (MeV)");
  hTedet =
      new Histogram(100, -150., 3.0, "T of ions used for gain recalibration",
                    "T (mm)", "ions");

  sprintf(fileName, "%s/Pedestals_%d.gp", Outputdir.c_str(), thread);
  sprintf(fileName2, "%s/Energies_%d.gp", Outputdir.c_str(), thread);
  sprintf(fileName3, "%s/E_Profile_%d.gp", Outputdir.c_str(), thread);
  threadNumber = thread;
}

void pedGainCalib::rawPh(pCTraw &rawEvt) { // Called for each raw event read in
                                           // from the input data file
  // Accumulate histograms for pedestal analysis
  hPed[0]->entry(rawEvt.enrg_fpga[0].pulse_sum[0]);
  hPed[1]->entry(rawEvt.enrg_fpga[0].pulse_sum[1]);
  hPed[2]->entry(rawEvt.enrg_fpga[0].pulse_sum[2]);
  hPed[3]->entry(rawEvt.enrg_fpga[1].pulse_sum[0]);
  hPed[4]->entry(rawEvt.enrg_fpga[1].pulse_sum[1]);
}

void pedGainCalib::getPeds(const char *inFileName, int run_number,
                           int program_version, float proj_angle, int nKeep,
                           string start_time) { // Called after comopletion of
                                                // the loop over all input raw
                                                // data

  // Save plots of the histograms so that the pedestal region can be visualized
  oFile = fopen(fileName, "w");
  if (oFile != NULL) {
    time_t t = time(NULL);
    struct tm *now = localtime(&t);
    fprintf(oFile, setTerm.c_str());
    fprintf(oFile, "set multiplot title 'PedGainCalib::getPeds: Energy "
                   "Detector Pedestal Region for file %s' layout 3,2 "
                   "columnsfirst scale 1.0,1.0\n",
            inFileName);
    hPed[0]->plot(oFile, true);
    hPed[1]->plot(oFile, true);
    hPed[2]->plot(oFile, true);
    hPed[3]->plot(oFile, true);
    fprintf(oFile, "set label 4 'Current date and time: %d-%d-%d    %d:%d' at "
                   "screen 0.55, 0.28\n",
            now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour,
            now->tm_min);
    fprintf(oFile, "set label 5 'Run Number %d' at screen 0.55, 0.25 left\n",
            run_number);
    fprintf(oFile,
            "set label 6 'FPGA Program Version %d' at screen 0.55, 0.22 left\n",
            program_version);
    fprintf(
        oFile,
        "set label 7 'Stage Angle = %5.1f degrees' at screen 0.55, 0.19 left\n",
        proj_angle);
    fprintf(oFile, "set label 8 'Start Time = %s ' at screen 0.55, 0.16 left\n",
            start_time.c_str());
    fprintf(
        oFile,
        "set label 9 'Number of good events = %d' at screen 0.55, 0.13 left\n",
        nKeep);
    hPed[4]->plot(oFile, true);
    fprintf(oFile, "unset multiplot\n");
    fprintf(oFile, "show label\n");
    fprintf(oFile, "unset label\n");
    fclose(oFile);
  } else {
    cout << "Unable to open the histogram file " << fileName << endl;
  }

  // oFile = fopen("allpeds.txt","a");
  // fprintf(oFile,"Run ,%d, Pedestals=,",run_number);
  for (int stage = 0; stage < 5; stage++) {
    float xLow, xHigh;
    int ret = hPed[stage]->FWHMboundaries(xLow, xHigh);
    if (ret == 0 && hPed[stage]->cnts(xLow, xHigh) > 100) {
      Ped[stage] = hPed[stage]->mean(xLow, xHigh);
    } else {
      Ped[stage] = 0.;
      cout << "pedGainCalib::getPeds thread " << threadNumber
           << ": ERROR: could not find the peak of the pedestal distribution "
              "for stage ************" << stage << endl;
    }
    cout << "getPeds thread " << threadNumber
         << ": measured pedestal for stage " << stage << " is " << Ped[stage]
         << endl;
    // fprintf(oFile,"%10.2f,",Ped[stage]);
  }
  // fprintf(oFile,"\n");
  // fclose(oFile);
}

void pedGainCalib::weplEvt(float Vedet, float Tedet,
                           float Ene[5]) { // Called for each event in the
                                           // temporary file of proton histories
  float Esum = 0.;
  if ((Tedet > emtrng1 && Tedet < emtrng2) ||
      (Tedet > emtrng3 && Tedet < emtrng4)) { // Analyze full-energy protons
                                              // outside of the phantom region
    if (fabs(Vedet) < 40.) {
      for (int stage = 0; stage < 5; stage++) {
        float stgEne = Ene[stage];
        hEnrg[stage]->entry(stgEne);
        Esum += Ene[stage];
      }
      hEnrgTot->entry(Esum);
      hTedet->entry(Tedet);
    }
  }
  if (Ene[0] > 10.)
    hProfT->entry(Tedet, Ene[0]);
}

void pedGainCalib::getGains(TVcorrection *TVcorr, const char *inFileName,
                            int run_number, int program_version, int proj_angle,
                            int nKeep,
                            string start_time) { // Called prior to the final
                                                 // loop over protons histories
                                                 // to calculate WEPL
  // Save plots of the histograms so that the gains can be visualized
  oFile = fopen(fileName2, "w");
  if (oFile != NULL) {
    time_t t = time(NULL);
    struct tm *now = localtime(&t);
    fprintf(oFile, setTerm.c_str());
    fprintf(oFile, "set multiplot title 'PedGainCalib::getGains: Energy "
                   "Detector Stage Energies for file %s' layout 3,2 "
                   "columnsfirst scale 1.0,1.0\n",
            inFileName);
    hEnrg[0]->plot(oFile);
    hEnrg[1]->plot(oFile);
    hEnrg[2]->plot(oFile);
    hEnrg[3]->plot(oFile);
    hEnrg[4]->plot(oFile);
    hEnrgTot->plot(oFile);
    fprintf(oFile, "unset multiplot\n");
    fprintf(oFile, "show label\n");
    fprintf(oFile, "unset label\n");
    fclose(oFile);
  } else {
    cout << "Unable to open the histogram file " << fileName2 << endl;
  }

  oFile = fopen(fileName3, "w");
  if (oFile != NULL) {
    fprintf(oFile, setTerm.c_str());
    fprintf(oFile, "set multiplot title 'Energy Detector for file %s' layout "
                   "2,1 columnsfirst scale 1.0,1.0\n",
            inFileName);
    hTedet->plot(oFile);
    hProfT->plot(oFile);
    fprintf(oFile, "unset multiplot\n");
    fprintf(oFile, "show label\n");
    fprintf(oFile, "unset label\n");
    fclose(oFile);
  } else {
    cout << "getGains thread " << threadNumber
         << ": unable to open the profile plot file " << fileName3 << endl;
  }

  double EBragg = 6.0;
  double Peak[5];
  double rng[5] = { 0.5, 0.5, 0.5, 0.5, 1.0 };
  for (int stage = 0; stage < 5; stage++) {
    float xLow, xHigh;
    int ret = hEnrg[stage]->FWHMboundaries(xLow, xHigh);
    if (ret == 0 && hEnrg[stage]->cnts(xLow, xHigh) > 100) {
      Peak[stage] = hEnrg[stage]->mean(xLow, xHigh);
    } else {
      Peak[stage] = 0.0;
      //            double mode= hEnrg[stage]->mode();
      // // The location of the peak in the energy distribution
      //            Peak[stage]=
      // hEnrg[stage]->mean(mode-rng[stage]*EBragg,mode+rng[stage]*EBragg);  //
      // Mean of the distribution around the peak
      //            cout << "getGains thread " << threadNumber << ": Could not
      // find the FWHM points of the energy distribution for stage " << stage <<
      // endl;
      //            cout << "getGains thread " << threadNumber << ": mode= " <<
      // mode << " and peak=" << Peak[stage] << endl;
    }
    cout << "getGains thread " << threadNumber
         << ": measured peak location for stage " << stage << " is "
         << Peak[stage] << endl;
    if (Peak[stage] > 0.01) {
      corrFac[stage] = TVcorr->Eempt[stage] / Peak[stage];
    } else {
      cout << "getGains thread " << threadNumber
           << ", ERROR: unable to find a peak position; leaving the gains "
              "unchanged. **********" << endl;
      corrFac[stage] = 1.0;
    }
  }
}
