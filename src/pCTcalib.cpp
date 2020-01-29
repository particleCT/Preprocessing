// Code to process raw calibration data to derive the TV and WEPL calibration constants.
#include "pCTcalib.h"
pCTcalib::~pCTcalib() {pCTcalibRootFile->Close(); }

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
pCTcalib::pCTcalib(pCTconfig cfg, string inputFileName): config(cfg)
{
  cout << "\nEntering pCTcalib, the constructor for the pCT calibration "<<endl;    
  cout << "Filename for the list of calibration data files = " << inputFileName << endl;
  cout << "The output directory is " << cfg.item_str["outputDir"] << endl;
  cout << "The maximum number of events is " << cfg.item_int["max_events"] << endl;
  cout << "The maximum time is " << cfg.item_int["max_time"] << endl;
  cout << "The number of debug event printouts is " << cfg.item_int["n_debug"] << endl;
  cout << "The number of events to plot is " << cfg.item_int["n_plot"] << endl;
  cout << "The input WEPL calibration file is " << cfg.item_str["Wcalib"] << endl;
  cout << "The input TV correction calibration file is " << cfg.item_str["TVcorr"] << endl;
  cout << "Input date range = " << cfg.item_int["minDate"] << " to " << cfg.item_int["maxDate"] << endl;
  cout << "The input run range is " << cfg.item_int["minRun"] << " to " << cfg.item_int["maxRun"] << endl;
  cout << "The calibration wedge offset is assumed to be " << cfg.item_float["wedgeoffset"] << endl;
  cout << "The particle type is assumed to be " << cfg.item_str["partType"] << endl;
  cout << "Real-time pedestal and gain recalibration is done? " << cfg.item_int["recalibrate"] << endl;
  for (int stage = 0; stage < nStage; ++stage) {
    cout << "  For stage " << stage << " the range in which to look for pedestals begins at " << cfg.item_int[Form("pedrng%d",stage)] << " ADC counts " << endl;
  }
  CalFile = inputFileName;
  if (config.item_str["partType"] == "H") EnergyBinWidth = 0.25;
  else EnergyBinWidth = 1.0;
  RangeBinWidth = 1.0;

  // LOW E INTERPOLATION: NEVER USED DUE OF THRESHOLD OF 1MeV
  if (config.item_str["partType"] == "He"){ k1[0] = 0; k1[1] = 2; k1[2] = 2; k1[3] = 2; k1[4] = 3;}
  else{ k1[0] = 0; k1[1] = 5; k1[2] = 5; k1[3] = 5; k1[4] = 6;}
  
  if (config.item_str["partType"] == "He") {
    // EXTRAPOLATION AT HIGH ENERGIES (LOW STATISTICS THERE) 
    i1[0] = 230; i1[1] = 230; i1[2] = 247;
    i1[3] = 247; i1[4] = 250; 

    i2[0] = 242; i2[1] = 240; i2[2] = 255;
    i2[3] = 255; i2[4] = 265;
    i3    = 300;
  }
  else {
    // NO EXTRAPOLATION
    i1[0] = 280; i1[1] = 280; i1[2] = 288; i1[3] = 288; i1[4] = 200;
    i2[0] = 292; i2[1] = 292; i2[2] = 304;
    i2[3] = 304; i2[4] = 240; i3    = 330;   
  }
  if (config.item_str["partType"] == "H") {
    EG4stage[0] = 25.25; EG4stage[1] = 28.01; EG4stage[2] = 32.76;// MC derived stage energies, used to calibrate to MeV
    EG4stage[3] = 42.62; EG4stage[4] = 67.71;                     // for protons
  } else {
    EG4stage[0] = 100.; EG4stage[1] = 111.; EG4stage[2] = 129.;// MC derived stage energies, used to calibrate to MeV for He
    EG4stage[3] = 166.; EG4stage[4] = 279.;
  }
  cout << "pCTcalib: reading the list of raw data file names from " << CalFile << endl;
  cout << "pCTcalib: the list of raw data file names for calibration is as follows:" << endl;

  // Read the CalFile to find the calibration files
  int linecount = 0;
  string line;
  ifstream infile(CalFile);
  if (infile) {
    while (getline(infile, line)) {
      size_t found = line.find_first_not_of(" ");
      size_t end = line.find_first_of(" \n", found);
      string vfn = line.substr(found, end);
      FILE *cFile;
      cFile = fopen(vfn.c_str(), "r");
      if (cFile == NULL) {
	cout << "Cannot open the calibration raw data file '" << line << "' " << endl;
	exit(3);
      }
      cout << "Calibration raw data file name " << linecount << " is " << vfn << endl;
      calFileNames.push_back(line);
      fclose(cFile);
      linecount++;
      
    }
    if (linecount < 6) {
      cout << "Failed to find all 6 filename paths for the calibration; linecount=" << linecount << endl;
      return;
    }
  } else {
    cout << "Failed to open the list of calibration files in " << CalFile << endl;
    return;
  }
  
  infile.close();
  if (calFileNames.size() == 0) return;

  cout << "The beam particle type is " << config.item_str["partType"] << endl;
  float EG4tot = 0.;
  for (int i = 0; i < nStage; ++i) {
    cout << "MC derived stage energy for an empty event, for stage " << i << " is " << EG4stage[i] << endl;
    EG4tot += EG4stage[i];
  }
  cout << "MC derived total energy for an empty event is " << EG4tot << endl << endl;
  currentTime = time(NULL);
  now = localtime(&currentTime);  
  // Initialize class
  theGeometry    = new pCTgeo(config.item_float["wedgeoffset"]);                 
  theCuts        = new pCTcut(config);
  theTVcorr      = new TVcorrection(pCTcalibRootFile, 1);
  theEvtRecon    = new EvtRecon(config);


}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int pCTcalib::TVmapper() {
  if (calFileNames.size() == 0) {
    cout << "TVmaper: only " << calFileNames.size() << " calibration file names found.  Need at least 1 for the TV map."
         << endl;
    return -2;
  }
  config.item_int["Nbricks"] = 0; // This is for the empty one for the TV correction
  config.item_int["doGains"] = false;


  // Calibration
  float wedgeLimit = theGeometry->getTWedgeBreaks(4) + 25.0; // NO BRICK OFFSET //+5.
  float openRange = 20.0;
  float pedestals[nStage];
  int pedMin[5];
  for (int stage = 0; stage < nStage; stage++) pedestals[stage] = 0.;
  for (int stage = 0; stage < nStage; stage++){
    pedMin[stage] = config.item_int[Form("pedrng%d",stage)];
  }
  pedGainCalib* theCalibration = new pedGainCalib(pCTcalibRootFile, pedMin, pedestals,-150., -151., wedgeLimit, wedgeLimit + openRange, config.item_str["partType"]);  
  theEvtRecon->ReadInputFile(theGeometry, theTVcorr, calFileNames[0], theCalibration); // Pedestal are determined here

  // Position of the trackers
  Ut[0] = theEvtRecon->uhitT[2]; Ut[1] = theEvtRecon->uhitT[3];
  Uv[0] = theEvtRecon->uhitV[2]; Uv[1] = theEvtRecon->uhitV[3];

  cout << "pCTcalib: U values for T coordinates = " << Ut[0] << " and " << Ut[1] << endl;
  cout << "pCTcalib: U values for V coordinates = " << Uv[0] << " and " << Uv[1] << endl;  
  for (int stage = 0; stage < nStage; ++stage) cout << "Pedestal measured for stage " << stage << " is " << theEvtRecon->Peds[stage] << " ADC counts " << endl;
  
  // Define a signal histogram for each pixel
  float bin0[nStage] = { 1000., 1000., 1000., 1000., 1000. };
  for (int Stage = 0; Stage < nStage; Stage++) {
    for (int pix = 0; pix < nPixRoot ; pix++) { // 480 with the underflow overflow 
      string Title = "Stage " + to_string((long long int)Stage) + ", Pixel " + to_string((long long int)pix);
      if(config.item_str["partType"] == "H") pxHistADC_root[Stage][pix] = new TH1D(Title.c_str(), Title.c_str(), 2000, bin0[Stage], bin0[Stage] +4.*2000);
      else pxHistADC_root[Stage][pix] = new TH1D(Title.c_str(), Title.c_str(), 2000, bin0[Stage], bin0[Stage] +8.*2000);
    }
  }
  TH2D* hTVmap_root = new TH2D("Profile of tracks at energy detector","",150, -150, 150 , 80, -40, 40);
  cout << "TVmapper: starting event loop" << endl;
  for (int EvtNum = 0; EvtNum < theEvtRecon->nEvents; EvtNum++) {
    Event thisEvent;
    thisEvent = theEvtRecon->evtList[EvtNum];
    T[0] = thisEvent.Thit[2]; T[1] = thisEvent.Thit[3];
    V[0] = thisEvent.Vhit[2]; V[1] = thisEvent.Vhit[3];

    if (EvtNum % 1000000 == 0) {
      cout << "  Processing event " << EvtNum << endl;
      theEvtRecon->dumpEvt(thisEvent);
    }
    hTVmap_root->Fill(T[1], V[1]);
    for (int stage = 0; stage < nStage; stage++) {
      // Extrapolate the rear track vector to the energy detector stage
      double Tcorr = theGeometry->extrap2D(Ut, T, theGeometry->energyDetectorU(stage));
      double Vcorr = theGeometry->extrap2D(Uv, V, theGeometry->energyDetectorU(stage));
      // Outside of the detector
      if (fabs(Tcorr) > 189.0) continue;
      if (fabs(Vcorr) > 49.0)  continue;

      int global = theTVcorr->TVcorrHist[0]->FindBin(Tcorr,Vcorr);
      
      pxHistADC_root[stage][global]->Fill(thisEvent.ADC[stage] - theEvtRecon->Peds[stage]);
    }
  }
  pCTcalibRootFile->cd("");
  pCTcalibRootFile->mkdir("tracksProfile");
  pCTcalibRootFile->cd("tracksProfile");
  hTVmap_root->Write("");
  
  // Save plots of the histograms
  pCTcalibRootFile->mkdir("ADC_Stage");
  pCTcalibRootFile->cd("ADC_Stage");
  cout << "TVmapper: Starting evaluation of the TV correction factors." << endl;

  //Fill the TV corrs
  for (int pix = 0; pix < nPixRoot; pix++) {
    cout << "TVmapper pixel stage ADC values:  " << pix;
    for (int stage = 0; stage < nStage; stage++) {
      pxHistADC_root[stage][pix]->Write("");
      float xLow, xHigh;
      float xpeak, xmax, max;
      if(pxHistADC_root[stage][pix]->GetEntries()>10){
	int imax;
	imax  =  pxHistADC_root[stage][pix]->GetMaximumBin();
	max   =  pxHistADC_root[stage][pix]->GetMaximum();
	xmax  =  pxHistADC_root[stage][pix]->GetBinCenter(imax); 
	TF1 *f1 = new TF1("f1", "gaus", xmax-100, xmax+100);
	f1->SetParameter(0, max);
	f1->SetParameter(1, xmax);
	Int_t fitStatus = pxHistADC_root[stage][pix]->Fit(f1,"QR"); // Q means quiet, R is 
	if(fitStatus==0) // Everything passed
	  {
	    xpeak  = f1->GetParameter(1);// mean
	    xHigh  = xpeak + f1->GetParameter(2); //sigma
	    xLow   = xpeak - f1->GetParameter(2); //sigma
	    theTVcorr->TVcorrHist[stage]->SetBinContent(pix, EG4stage[stage] / xpeak);
	  }
	else theTVcorr->TVcorrHist[stage]->SetBinContent(pix, 0.999);
      }
      else theTVcorr->TVcorrHist[stage]->SetBinContent(pix, 0.999);
      cout << " stg" << stage << "=       "<< theTVcorr->TVcorrHist[stage]->GetBinContent(pix);
      delete pxHistADC_root[stage][pix]; // Free up the memory that was sucked up by the pixel histograms
   }
   cout << endl;
  }

  pCTcalibRootFile->cd("");
  pCTcalibRootFile->mkdir("TVcorrProfile");
  pCTcalibRootFile->cd("TVcorrProfile");  
  for(int stage =0; stage<nStage; stage++) theTVcorr->TVcorrHist[stage]->Write("", TObject::kOverwrite);
  cout << "TVmapper: Finished with mapping the TV calibration" << endl;

  delete (theCalibration);
  return 0;
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void pCTcalib::enrgDep() { // Calculate energy depositions in each stage and their sum using new TV maps, to check that calibration worked
  // Create Histogram  for each stage
  for (int stage = 0; stage < nStage; ++stage) {
    string Title = "pCTcalib::enrgDep Energy of stage " + to_string((long long int)stage);
    if (config.item_str["partType"] == "H"){
      stgHistE_root[stage]  = new TH1D(Title.c_str(), Title.c_str(),  400, 15, 15 + 400*0.175); 
      stgHistE[stage]       = new Histogram(400, 15., 0.175, Title, "E (MeV)", "protons");
    }
    else{
      stgHistE_root[stage]  = new TH1D(Title.c_str(), Title.c_str(),  400, 15, 15 + 400*0.9); 
      stgHistE[stage]       = new Histogram(400, 15., 0.9, Title, "E (MeV)", "He ions");
    }
    Title = "pCTcalib::enrgDep Energy of stage " + to_string((long long int)stage);
    stgE[stage] = new TProfile2D(Form("Profile Energy of stage %d",stage), "", 100, -150, -150+3.0*100, 100, -50., -50 + 1.0*100);
  }
  // Sum of stage energies
  if (config.item_str["partType"] == "H") EsumH = new TH1D("pCTcalib::enrgDep Sum of stage energies","pCTcalib::enrgDep Sum of stage energies", 900, 0, 0.25*900);
  else EsumH = new TH1D("pCTcalib::enrgDep Sum of stage energies", "pCTcalib::enrgDep Sum of stage energies", 900, 0, 1.0*900);

  // Create Histogram for each pixel
  for (int iPix = 0; iPix < nPix; iPix++) {
    for (int stage = 0; stage < nStage; ++stage) {
      string Title = "pCTcalib::enrgDep Energy of stage " + to_string((long long int)stage) + " for pixel " +
                     to_string((long long int)iPix);
      if (config.item_str["partType"] == "H"){
        pxHistE[stage][iPix] = new Histogram(200, 15., .35, Title, "E (MeV)", "protons");
	pxHistE_root[stage][iPix] = new TH1D(Title.c_str(), Title.c_str(), 200, 15., 15 +200*0.35);
      }
      else{
        pxHistE[stage][iPix] = new Histogram(200, 60., 1.4, Title, "E (MeV)", "He ions");
	pxHistE_root[stage][iPix] = new TH1D(Title.c_str(), Title.c_str(), 200, 60., 60 +200*1.4);
      }
    }
  }
  
  float mxEstage, mnEstage;
  if (config.item_str["partType"] == "H") { mnEstage = 15.; mxEstage = 100.; }
  else { mnEstage = 60.; mxEstage = 400.; }
  
  for (int EvtNum = 0; EvtNum < theEvtRecon->nEvents; EvtNum++) { // loop over events
    if (EvtNum % 1000000 == 0) cout << "pCTcalib::enrgDep, Processing event " << EvtNum << endl;
    Event thisEvent;
    thisEvent = theEvtRecon->evtList[EvtNum];
    T[0] = thisEvent.Thit[2]; T[1] = thisEvent.Thit[3];
    V[0] = thisEvent.Vhit[2]; V[1] = thisEvent.Vhit[3];

    float Esum = 0.;
    for (int stage = 0; stage < nStage; stage++) {
      // Extrapolate the rear track vector to the energy detector stage
      float Tcorr = theGeometry->extrap2D(Ut, T, theGeometry->energyDetectorU(stage));
      float Vcorr = theGeometry->extrap2D(Uv, V, theGeometry->energyDetectorU(stage));
      if (fabs(Tcorr) > 150.0) continue;
      if (fabs(Vcorr) > 45.0)  continue; 

      bool inBounds = true;
      float Ecorr = ((float)thisEvent.ADC[stage] - theEvtRecon->Peds[stage]) * theTVcorr->corrFactor(stage, Tcorr, Vcorr, inBounds);
      if (!inBounds) continue;

      int tPix = floor(0.1 * (Tcorr + 190.));
      int vPix = floor(0.1 * (Vcorr + 50.));
      int iPix = vPix + 10 * tPix;
      if (iPix >= nPix) continue;
      if (vPix > 9) {
        cout << "pCTcalib: iPix=" << iPix << " vPix=" << vPix << endl;
        perror("pCTcalib: invalid pixel");
        iPix = 0;
      }
      pxHistE[stage][iPix]->entry(Ecorr);
      pxHistE_root[stage][iPix]->Fill(Ecorr);
      stgHistE_root[stage]->Fill(Ecorr);
      stgHistE[stage]->entry(Ecorr);
      //if (Ecorr > mnEstage && Ecorr < mxEstage) {
      //if (fabs(Vcorr) < 35.0)  stgEvsT[stage]->Fill(Tcorr, Ecorr);
      //if (fabs(Tcorr) < 140.0) stgEvsV[stage]->Fill(Vcorr, Ecorr);
      if (fabs(Vcorr) < 35.0 && fabs(Tcorr) < 140.0)  stgE[stage]->Fill(Tcorr, Vcorr, Ecorr);
	//}
      Esum += Ecorr;
    }
    EsumH->Fill(Esum);
  }
  pCTcalibRootFile->mkdir("EnergyHistogram");  
  pCTcalibRootFile->cd("EnergyHistogram");
  for(int i =0; i<5; i++) stgHistE_root[i]->Write("", TObject::kOverwrite);
  EsumH->Write("TotalEnergy",TObject::kOverwrite);
  pCTcalibRootFile->cd();

  pCTcalibRootFile->mkdir("EnergyProfile");  
  pCTcalibRootFile->cd("EnergyProfile");
  for(int i =0; i<5; i++) stgE[i]->Write("", TObject::kOverwrite);
  pCTcalibRootFile->cd();  
  //
  for (int stage = 0; stage < nStage; ++stage) {
    float xLow, xHigh, Sadc;
    /*if(stgHistE_root[stage]->GetEntries()>500){
      Sadc = stgHistE_root[stage]->GetBinCenter(stgHistE_root[stage]->GetMaximumBin());
      }*/
    int ret = stgHistE[stage]->FWHMboundaries(xLow, xHigh);
    if (ret == 0) Sadc = stgHistE[stage]->mean(xLow, xHigh); 
    else {
      Sadc = 0.;
      cout << "pCTcalib::enrgDep, no data for stage " << stage << endl;
    }
    Est[stage] = Sadc;
    theTVcorr->ped[stage]   = theEvtRecon->Peds[stage];
    theTVcorr->Eempt[stage] = Est[stage];
    cout << "pCTcalib::enrgDep, Energy in stage " << stage << " = " << Est[stage] << "; dE=" << Est[stage] - EG4stage[stage] << " MeV" << endl;
    delete stgHistE[stage];
    //delete stgEvsT[stage];
    //delete stgEvsV[stage];
    delete stgE[stage];
  }
  cout << "pCTcalib::enrgDep, MC - experiment difference |dE| should be <0.1 MeV ! " << endl;
  EnS = EsumH->GetBinCenter(EsumH->GetMaximumBin());//mean(xLow, xHigh);
  cout << "pCTcalib::enrgDep, Energy sum = " << EnS << endl;
  delete EsumH;
  cout << "pCTcalib::enrgDep,  Done with the TV calibration task" << endl;
}
//////////////////////////////////////////////////////////////////////
// Wcalib function
//////////////////////////////////////////////////////////////////////
int pCTcalib::Wcalib(){
  cout << endl << "Entering pCTcalib::Wcalib to execute the WEPL calibration process" << endl;
  if (calFileNames.size() != 6) {
    cout << "Only " << calFileNames.size() << " calibration raw data file names found.  Need all 6 to include "
      "the WEPL calibration." << endl;
    return -1;
  }
  for (int nBricks = 0; nBricks < 5; nBricks++) cout << "The raw data file for " << nBricks << " bricks is " << calFileNames[nBricks] << endl;
  TH2D* dEEhist[nStage][5];
  TH2D* REhist_root[nStage][5];
  struct { Histogram2D *h[nStage];}  REhist[5];// Range Energy histogram
  for (int nBricks = 0; nBricks < 5; ++nBricks) {
    for (int stage = 0; stage < nStage; ++stage) {

      dEEhist[nBricks][stage] = new TH2D(Form("dE-EStage_%d_bricks_%d",stage,nBricks), Form("dE-E spectra for stage %d and %d bricks ", stage, nBricks), // Title
					      nEnrg, 0., nEnrg * EnergyBinWidth, nEnrg, 0., nEnrg * EnergyBinWidth);
					      
      float enrgB0 = 0.; // The analysis of energy slices further down assumes that the first energy bin starts at zero
      string Title = "WEPL calibration array for stage " + to_string((long long int)stage) + " and nBricks=" + to_string((long long int)nBricks);
      //REhist_root[nBricks][stage] = new TH2D(Form("WEPLcalib_stage%d_bricks%d",stage,nBricks), Form("WEPLcalib_stage%d_bricks%d",stage,nBricks),
      //nRange, 0, RangeBinWidth*nRange, nEnrg, enrgB0, enrgB0 + nEnrg*EnergyBinWidth);
      REhist_root[nBricks][stage] = new TH2D(Form("WEPLcalib_stage%d_bricks%d",stage,nBricks), Form("WEPLcalib_stage%d_bricks%d",stage,nBricks),
					     nRange, 0, RangeBinWidth*nRange, nEnrg, enrgB0, enrgB0 + nEnrg*EnergyBinWidth);
      REhist_root[nBricks][stage]->GetXaxis()->SetTitle("Range (mm)");
      REhist_root[nBricks][stage]->GetYaxis()->SetTitle("Energy (MeV)");      

      if (config.item_str["partType"] == "H") REhist[nBricks].h[stage] = new Histogram2D(nRange, 0., RangeBinWidth, nEnrg, enrgB0, EnergyBinWidth, Title,
											 "proton range (mm)", "energy (MeV)", "counts");
      else REhist[nBricks].h[stage] = new Histogram2D(nRange, 0., RangeBinWidth, nEnrg, enrgB0, EnergyBinWidth, Title,
						      "He range (mm)", "energy (MeV)", "counts");
    }
  }
  // Submit the event processing  -- fill the histograms
  for (int nBricks = 1; nBricks < 5; nBricks++) {
    config.item_str["inputFileName"] = calFileNames[nBricks + 1];
    config.item_int["Nbricks"] = nBricks;
    procWEPLcal(REhist[nBricks].h, REhist_root[nBricks] , dEEhist[nBricks]);
  }
  config.item_str["inputFileName"] = calFileNames[1];
  config.item_int["Nbricks"] = 0;
  procWEPLcal(REhist[0].h, REhist_root[0], dEEhist[0]);

  cout << "pCTcalib::Wcalib, begin adding together the range-energy tables from the runs with different numbers of bricks." << endl;
  // Combine the resulting maps into one map for each stage by simply adding them
  pCTcalibRootFile->mkdir("dEE");
  pCTcalibRootFile->mkdir("REhist");
  for (int stage = 0; stage < nStage; ++stage) {
    
    // Save the individual histograms
    pCTcalibRootFile->cd("dEE");
    for (int nBricks = 0; nBricks < 5; ++nBricks) dEEhist[nBricks][stage]->Write("",TObject::kOverwrite);
    pCTcalibRootFile->cd("REhist");
    for (int nBricks = 0; nBricks < 5; ++nBricks) REhist_root[nBricks][stage]->Write("",TObject::kOverwrite);
    
    //Add them up together
    for (int nBricks = 1; nBricks < 5; ++nBricks) {
      string Title = "Summed WEPL calibration array for stage " + to_string((long long int)stage);
      REhist[0].h[stage]->add(REhist[nBricks].h[stage],Title);
      delete REhist[nBricks].h[stage];
      dEEhist[0][stage]->Add(dEEhist[nBricks][stage]);
      REhist_root[0][stage]->Add(REhist_root[nBricks][stage]);
    }

    // Save the added histograms
    pCTcalibRootFile->cd("dEE");
    dEEhist[0][stage]->Write(Form("dEE_Tot_stage%d",stage),TObject::kOverwrite);
    pCTcalibRootFile->cd("REhist");
    REhist_root[0][stage]->Write(Form("RE_Tot_stage%d",stage),TObject::kOverwrite);
  }
  // Analyze energy slices of the summed maps, and plot the slices -- Roberts
  for (int stage = 0; stage < nStage; ++stage) {
    int nrg = 0;
    for (int hst = 0; hst < nEnrg / 6 + 1; ++hst) { // Group by 6 slices, to plot 6 histograms per page
      for (int j = 0; j < 6; ++j) {
        if (nrg >= nEnrg) break;
	float xLow, xHigh, Sadc;
	Histogram Slice = REhist[0].h[stage]->rowSlice(nrg);
	int ret = Slice.FWHMboundaries(xLow, xHigh);
        if (ret == 0) {
	Sadc = Slice.mean(xLow, xHigh);
        }
	else {
          if (stage == 4) {
            if (Slice.imode() == 0)  // max
              Sadc = 0.; // Up against the wall at zero range
            else {
              cout << ret << " Cannot find the FWHM bounds for stage " << stage << ", slice " << nrg << endl;
              float mode = Slice.mode();
              Sadc = Slice.mean(mode - 2. * RangeBinWidth, mode + 2 * RangeBinWidth);
            }
          }
	  else {
            cout << ret << " Cannot find the FWHM bounds for stage " << stage << ", slice " << nrg << endl;
            float mode = Slice.mode();
            Sadc = Slice.mean(mode - 2. * RangeBinWidth, mode + 2 * RangeBinWidth);
          }
        }
	// Fit the Range from the adc value
        Rst[stage][nrg] = Sadc;
        Sst[stage][nrg] = (xHigh - xLow) / 2.2;
        nrg++;
	}
      }
  }
  // Analyze energy slices of the summed maps, and plot the slices -- Charles
  //Define a function to find the peaks
  
  /*for (int stage = 0; stage < nStage; ++stage) {
    for (int nrg = 0; nrg < nEnrg; nrg++) { 
      float xLow, xHigh, Sadc,xpeak;
      TH1D* RESlice =  REhist_root[0][stage]->ProjectionX("ProjX_Stage",nrg,nrg+1);
      Int_t NEntries= RESlice->GetEntries();
      xpeak  =  RESlice->GetBinCenter(RESlice->GetMaximumBin()); // cheesy fix because of high level of low energy data
      xHigh  =  RESlice->GetBinCenter(RESlice->FindLastBinAbove(  RESlice->GetMaximum()/2));
      xLow   =  xpeak - (xHigh -xpeak);	
      RESlice->GetXaxis()->SetRange(xLow,xHigh);  
      Sadc   =  RESlice->GetMean();
      if(NEntries>1000){

	xpeak  =  RESlice->GetBinCenter(RESlice->GetMaximumBin()); // cheesy fix because of high level of low energy data
	xHigh  =  RESlice->GetBinCenter(RESlice->FindLastBinAbove(  RESlice->GetMaximum()/2));
	xLow   =  xpeak - (xHigh -xpeak);	
	Sadc   =  xpeak;

	TF1 *f1 = new TF1("f1", "gaus", xLow, xHigh);
	f1->SetParameter(0, xpeak);
	Int_t fitStatus = RESlice->Fit(f1,"Q"); // Q means quiet
	if(fitStatus==0) // Everything passed
	  { 
	    //TF1* g = RESlice->GetFunction("gaus");
	    Sadc  = f1->GetParameter(1);// mean
	    xHigh = Sadc + f1->GetParameter(2); //sigma
	    xLow  = Sadc - f1->GetParameter(2); //sigma
	    //delete g;
	  }
      }

      else if (NEntries>200){
	xpeak       =  RESlice->GetBinCenter(RESlice->GetMaximumBin()); // cheesy fix because of high level of low energy data
	xHigh       =  RESlice->GetBinCenter(RESlice->FindLastBinAbove(  RESlice->GetMaximum()/2));
	xLow        =  xpeak - (xHigh -xpeak);
	Sadc        =  xpeak;
	}

	//else{Sadc  = 0;	 xLow  = 0;  xHigh = 0;}
      delete RESlice;
      
      Rst[stage][nrg] = Sadc;
      Sst[stage][nrg] = (xHigh - xLow) / 2.2;
      
    }
   }*/


  // Plot the final calibration results before corrections
  float rngB[nEnrg];
  pCTcalibRootFile->cd();
  pCTcalibRootFile->mkdir("RangeVsEnergy_Uncorr");
  pCTcalibRootFile->cd("RangeVsEnergy_Uncorr");
  TGraphErrors* rngEnrg_unCorr[nStage];
  
  for (int stage = 0; stage < nStage; ++stage) {
    string Title = "Range vs Energy Uncorrected for Stage " + to_string((long long int)stage);
    rngEnrg_unCorr[stage] = new TGraphErrors(nEnrg);
    rngEnrg_unCorr[stage]->SetTitle(Title.c_str());
    rngEnrg_unCorr[stage]->SetName(Form("RangeVsEnergy_Uncorr_%d",stage));

    for (int nrg = 0; nrg < nEnrg; ++nrg){
      rngB[nrg] = float(nrg) * EnergyBinWidth;
      rngEnrg_unCorr[stage]->SetPoint(nrg, double(rngB[nrg]), double(Rst[stage][nrg]));
      rngEnrg_unCorr[stage]->SetPointError(nrg, 0.,   double(Sst[stage][nrg]));
    }
    rngEnrg_unCorr[stage]->GetXaxis()->SetTitle("Energy (MeV)");
    rngEnrg_unCorr[stage]->GetYaxis()->SetTitle("Range (mm)");      
    rngEnrg_unCorr[stage]->Write("", TObject::kOverwrite);
  }
  // Correction of the range-energy curve starts here
  cout << "Wcalib: interpolation parameters:" << endl;
  cout << "        k1 is for extrapolating below the threshold" << endl;
  cout << "        j1 through j4 are for interpolating through a kink region (if there is one)" << endl;
  cout << "        i1 and i2 are for interpolating between stages or, for the "
          "last stage, extrapolating to negative WET" << endl;
  cout << "        These parameters can be set by comments in the WEPL calibration file." << endl;
    
  for (int stage = 0; stage < nStage; stage++) {
    cout << "Wcalib: stage " << stage << "; k1=" << k1[stage] << endl;
    cout << "Wcalib: stage " << stage << "; j1=" << j1[stage] << endl;
    cout << "Wcalib: stage " << stage << "; j2=" << j2[stage] << endl;
    cout << "Wcalib: stage " << stage << "; j3=" << j3[stage] << endl;
    cout << "Wcalib: stage " << stage << "; j4=" << j4[stage] << endl;
    cout << "Wcalib: stage " << stage << "; i1=" << i1[stage] << endl;
    cout << "Wcalib: stage " << stage << "; i2=" << i2[stage] << endl;
  }
  
  // Correct stages 1-4 for 1 MeV threshold (no data for E<1MeV, use first
  // non-zero bin)
  for (int stage = 1; stage < nStage; stage++) {
    for (int k = 0; k < k1[stage]; k++) {
      Rst[stage][k] = Rst[stage][k1[stage]];
    }
  }

  /// Lennart Volz, November 2018
  /// dE-E parameter evaluation:
  int Estep[3] = { 240, 120, 12 }; // work in the same way for helium and proton, but could also be set manually for optimization
  float E[3] = { Estep[0] * EnergyBinWidth, Estep[1] * EnergyBinWidth, Estep[2] * EnergyBinWidth };

  FILE oFile;
  string fn;
  for (int stage = 1; stage < 5; stage++) { // stage 
    float xlow[3], xhigh[3];
    float xpeak, xmax, max;
    for (int j = 0; j < 3; j++) {
      TH1D* dEESlice = dEEhist[0][stage]->ProjectionX(Form("ProjX_Stage%d_Bin%d",stage,Estep[j]),Estep[j],Estep[j]+1);
      int imax;
      imax  =  dEESlice->GetMaximumBin();
      max   =  dEESlice->GetMaximum();
      xmax  =  dEESlice->GetBinCenter(imax);
      TF1 *f1 = new TF1("f1", "gaus", xmax-10, xmax+10);
      f1->SetParameter(0, max);
      f1->SetParameter(1, xmax);
      Int_t fitStatus = dEESlice->Fit(f1,"QR"); // Q means quiet, R is
      if(fitStatus==0) // Everything passed
	{
	  xpeak     = f1->GetParameter(1);// mean
	  xhigh[j]  = xpeak + f1->GetParameter(2); //sigma
	  xlow[j]   = xpeak - f1->GetParameter(2); //sigma
	}
      
      /*float energy_cut = 130; // MeV
      int bin_cut = dEESlice->FindBin(energy_cut);
      dEESlice->GetXaxis()->SetRange(bin_cut,dEESlice->GetNbinsX());
      float xpeak =  dEESlice->GetBinCenter(dEESlice->GetMaximumBin());     
      xhigh[j]    =  dEESlice->GetBinCenter(dEESlice->FindLastBinAbove(dEESlice->GetMaximum()/2));
      xlow[j]     =  xpeak - (xhigh[j] -xpeak); */
      // extend to 3 sigma region instead of 1 FWHM
      xlow[j]     = xlow[j] - (xhigh[j] - xlow[j]) * 0.7848;
      xhigh[j]    =   xhigh[j] + (xhigh[j] - xlow[j]) * 0.7848;          
      
    }    
    dEElow[stage][0] = (E[0] * (xlow[2] - xlow[1]) + E[1] * (xlow[0] - xlow[2]) + E[2] * (xlow[1] - xlow[0])) /
                       ((E[0] - E[1]) * (E[0] - E[2]) * (E[1] - E[2]));
    dEElow[stage][1] = (xlow[1] - xlow[0]) / (E[1] - E[0]) - dEElow[stage][0] * (E[0] + E[1]);
    dEElow[stage][2] = xlow[0] - dEElow[stage][0] * E[0] * E[0] - dEElow[stage][1] * E[0];
    dEEhigh[stage][0] = (E[0] * (xhigh[2] - xhigh[1]) + E[1] * (xhigh[0] - xhigh[2]) + E[2] * (xhigh[1] - xhigh[0])) /
                        ((E[0] - E[1]) * (E[0] - E[2]) * (E[1] - E[2]));
    dEEhigh[stage][1] = (xhigh[1] - xhigh[0]) / (E[1] - E[0]) - dEEhigh[stage][0] * (E[0] + E[1]);
    dEEhigh[stage][2] = xhigh[0] - dEEhigh[stage][0] * E[0] * E[0] - dEEhigh[stage][1] * E[0];
  }  
  return 0;
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void pCTcalib::writeCalibfile() {

  ////////////////////////////////////////////////////////////////
  // TV Calibration
  ////////////////////////////////////////////////////////////////
  pCTcalibRootFile->cd();
  TTree* header = new TTree("header","meta-data");
  header->Branch("particle",&config.item_str["partType"]);
  header->Branch("minDate",&config.item_str["minDate"]);
  header->Branch("maxDate",&config.item_str["maxDate"]);
  header->Branch("outputDir",&config.item_str["outputDir"]);
  header->Branch("minrun",&config.item_int["minrun"], "minrun/I");
  header->Branch("maxrun",&config.item_int["maxrun"], "maxrun/I" );
  for(int i=0; i<=5; i++) header->Branch(Form("calFileName_%d",i),&calFileNames[i]);
  for (int stage = 0; stage < nStage; ++stage) header->Branch(Form("peds_%d",stage),&theEvtRecon->Peds[stage]);
  for (int stage = 0; stage < nStage; ++stage) header->Branch(Form("Est_%d",stage),&Est[stage]);
  header->Branch("Esum",&EnS, "Esum/F");
  header->Fill();
  header->Write("", TObject::kOverwrite);

  ////////////////////////////////////////////////////////////////
  // W Calibration
  ////////////////////////////////////////////////////////////////
  /// Write the Rst0 through Rst4 range vs energy tables to the calibration text file
  string Wfilename = config.item_str["outputDir"] + "/" + config.item_str["Wcalib"];
  cout << "Opening output calibration file: " << Wfilename << endl;
  ofstream Wcalfile;
  Wcalfile.open(Wfilename, ios::out | ios::trunc); // text file R vs E tables
  if (!Wcalfile.is_open()) {
    cout << "Error opening W correction output file " << Wfilename << endl;
    cout << "Trying to open " << config.item_str["Wcalib"] << " instead. . ." << endl;
    Wcalfile.open(config.item_str["Wcalib"], ios::out | ios::trunc);
    if (Wcalfile.is_open())
      cout << "Successful file open!" << endl;
    else {
      cout << "Rats!  That didn't work either." << endl;
      exit(1);
    }
  }
  
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k) Wcalfile << Rst[0][k + i * 10] << " ";
    Wcalfile << endl;
  }
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k) Wcalfile << Rst[1][k + i * 10] << " ";
    Wcalfile << endl;
  }
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k) Wcalfile << Rst[2][k + i * 10] << " ";
    Wcalfile << endl;
  }
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k) Wcalfile << Rst[3][k + i * 10] << " ";
    Wcalfile << endl;
  }
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k) Wcalfile << Rst[4][k + i * 10] << " ";
    Wcalfile << endl;
  }
  
  Wcalfile << endl;
  Wcalfile << "# particle = " << config.item_str["partType"] << endl;
  Wcalfile << "# EnergyBinWidth = " << EnergyBinWidth << endl;
  Wcalfile << "# RangeBinWidth = " << RangeBinWidth << endl;
  Wcalfile << "# Valid date and run ranges for this calibration:" << endl;
  Wcalfile << "# minDate = " << config.item_str["mindate"] << endl;
  Wcalfile << "# maxDate = " << config.item_str["maxdate"] << endl;
  Wcalfile << "# minRun = " << config.item_int["minrun"] << endl;
  Wcalfile << "# maxRun = " << config.item_int["maxrun"] << endl;
  Wcalfile << "# Date and time of calibration run: " << now->tm_year + 1900 << "-" << now->tm_mon + 1 << "-"
	   << now->tm_mday << "    " << now->tm_hour << ":" << now->tm_min << endl;
  Wcalfile << "# File used for empty run for calibrating TV corrections was " << calFileNames[0] << endl;
  Wcalfile << "# The file to where the TV corrections were written was " << config.item_str["outputDir"] + "/TVcorr.txt" << endl;
  Wcalfile << "# File used for calibration with wedge and 0 bricks was " << calFileNames[1] << endl;
  Wcalfile << "# File used for calibration with wedge and 1 brick  was " << calFileNames[2] << endl;
  Wcalfile << "# File used for calibration with wedge and 2 bricks was " << calFileNames[3] << endl;
  Wcalfile << "# File used for calibration with wedge and 3 bricks was " << calFileNames[4] << endl;
  Wcalfile << "# File used for calibration with wedge and 4 bricks was " << calFileNames[5] << endl;
  Wcalfile << "# Threshold0 = " << config.item_float["thr0"] << endl;
  Wcalfile << "# Threshold1 = " << config.item_float["thr1"] << endl;
  Wcalfile << "# Threshold2 = " << config.item_float["thr2"] << endl;
  Wcalfile << "# Threshold3 = " << config.item_float["thr3"] << endl;
  Wcalfile << "# Threshold4 = " << config.item_float["thr4"] << endl;
  for (int stage = 1; stage < 5; stage++) {
    Wcalfile << "# dEE" + to_string(stage) + " = ";
    for (int i = 0; i < 3; i++)
      Wcalfile << dEElow[stage][i] << ",";
    for (int i = 0; i < 3; i++)
      Wcalfile << dEEhigh[stage][i] << ",";
    Wcalfile << endl;
  }
  Wcalfile.close();
};
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void pCTcalib::procWEPLcal(Histogram2D *REhist[nStage],TH2D* REhist_root[nStage], TH2D* dEEhist[nStage]) {
  // Routine passed to process the WEPL calibration.
  cout << "Entering procWEPLcal for Nbricks=" << config.item_int["Nbricks"] << ".  The list of histograms to fill is" << endl;
  for(int stage = 0; stage < nStage; ++stage) cout << config.item_int["Nbricks"] << " bricks, stage " << stage << ", title=  " << REhist[stage]->Title() << endl;
  theEvtRecon->config.item_int["doGains"] = true;
  theEvtRecon->config.item_int["Nbricks"] = config.item_int["Nbricks"];

  // Calibration
  float wedgeLimit = theGeometry->getTWedgeBreaks(4) + 25.0; // NO BRICK OFFSET //+5.
  float openRange = 20.0;
  float pedestals[nStage];
  int pedMin[5];
  for (int stage = 0; stage < nStage; stage++) pedestals[stage] = 0.;
  for (int stage = 0; stage < nStage; stage++){
    pedMin[stage] = config.item_int[Form("pedrng%d",stage)];
  }
  pedGainCalib* theCalibration = new pedGainCalib(pCTcalibRootFile, pedMin, pedestals,-150., -151., wedgeLimit, wedgeLimit + openRange, config.item_str["partType"]);
  //theCalibration->ClearHist();
  theEvtRecon->ReadInputFile(theGeometry, theTVcorr, config.item_str["inputFileName"], theCalibration);

  double Ut[2] , Uv[2] , T[2] , V[2];
  double Uft[2], Ufv[2], Tf[2], Vf[2];

  Uft[0] = theEvtRecon->uhitT[0]; Uft[1] = theEvtRecon->uhitT[1];
  Ufv[0] = theEvtRecon->uhitV[0]; Ufv[1] = theEvtRecon->uhitV[1];
  Ut[0]  = theEvtRecon->uhitT[2]; Ut[1]  = theEvtRecon->uhitT[3];
  Uv[0]  = theEvtRecon->uhitV[2]; Uv[1]  = theEvtRecon->uhitV[3];
  
  cout << "procWEPLcal: Tracker plane u locations:" << endl;
  cout << "Front tracker t: " << Uft[0] << " " << Uft[1] << endl;
  cout << "Front tracker v: " << Ufv[0] << " " << Ufv[1] << endl;
  cout << "Front tracker t: " << Ut[0] << " " << Ut[1] << endl;
  cout << "Front tracker t: " << Uv[0] << " " << Uv[1] << endl;

  //Define histograms
  Histogram *hStgEcorr[nStage];
  TH1D *hStgEcorr_root[nStage];
  for (int stage = 0; stage < nStage; ++stage) {
    cout << "procWEPLcal " << config.item_int["Nbricks"] << ": Pedestal measured for stage " << stage << " is " << theEvtRecon->Peds[stage] << " ADC counts " << endl;
    cout << "procWEPLcal " << config.item_int["Nbricks"] << ": Gain correction factor for stage " << stage << " is " << theEvtRecon->GainFac[stage] << endl;
    string Title = "nBricks= " + to_string((long long int)config.item_int["Nbricks"]) + " corrected energy for stage " +
      to_string((long long int)stage);
    if (config.item_str["partType"] == "H"){
      hStgEcorr[stage] = new Histogram(400, 15., 0.175, Title, "Energy (MeV)", "protons");
      hStgEcorr_root[stage] = new TH1D(Title.c_str(), Title.c_str(), 400, 15.,  15. +400*0.175);
    }
    else{
      hStgEcorr[stage] = new Histogram(400, 15., 0.7, Title, "Energy (MeV)", "He ions");
      hStgEcorr_root[stage] = new TH1D(Title.c_str(), Title.c_str(), 400, 15.,  15. +400*0.7);
    }
  }

  //Define more histograms
  string Title = "nBricks " + to_string((long long int)config.item_int["Nbricks"]) + " sum of corrected stage energies";
  Histogram *hTotEcorr;
  TH1D *hTotEcorr_root;
  if (config.item_str["partType"] == "H"){
    hTotEcorr = new Histogram(800, 0., 0.3, Title, "Energy (MeV)", "protons");
    hTotEcorr_root = new TH1D(Title.c_str(), Title.c_str(), 800, 0., 0.3*800);
    }
  else{
    hTotEcorr = new Histogram(840, 0., 1.0, Title, "Energy (MeV)", "He ions");
    hTotEcorr_root = new TH1D(Title.c_str(), Title.c_str(), 840, 0., 1.0*840);
  }

  Title = "nBricks " + to_string((long long int)config.item_int["Nbricks"]) + " calibration phantom thickness";
  TProfile* hLengthProf = new TProfile(Title.c_str(), Title.c_str(), 100, -150, -150 +3.0*100);

  // Wedge phantom geometry and locations,
  double Tw1 = theGeometry->getTWedgeBreaks(1);
  double Tw2 = theGeometry->getTWedgeBreaks(2);
  double Tw3 = theGeometry->getTWedgeBreaks(3);
  double Tw4 = theGeometry->getTWedgeBreaks(4);
  double tBrickEnd = theGeometry->getTWedgeBreaks(4);
  double BrickThickness = theGeometry->getBrickThickness(); // Brick thickness; also the wedge thickness

  double Ust[2]; // Polystyrene wedge phantom U coordinates, front and back
  Ust[0] = -2.5 * BrickThickness;
  Ust[1] = Ust[0] + BrickThickness;
  double Wbricks = BrickThickness * config.item_int["Nbricks"];
  double Uout = Ust[1] + Wbricks; // U coordinate of the last brick, downstream side
  cout << "procWEPLcal: The wedge phantom goes from u=" << Ust[0] << " to u=" << Ust[1] << endl;
  cout << "procWEPLcal: The " << config.item_int["Nbricks"] << " bricks go from u=" << Ust[1] << " to u=" << Uout << endl;
  cout << "procWEPLcal: Stage thresholds are " << config.item_float["thr0"] << " " << config.item_float["thr1"] << " " << config.item_float["thr2"] << " "
       << config.item_float["thr3"] << " " << config.item_float["thr4"] << endl;
  cout << "procWEPLcal: The wedge phantom break points are, in mm: " << Tw1 << ", " << Tw2 << ", " << Tw3 << ", " << Tw4 << endl;
  cout << "procWEPLcal: The calibration brick thickness is " << BrickThickness << " mm " << endl;

  float cut0, cut1, cut3; // energy cuts to reject events incompatible with expected Bragg peak.
  if (config.item_str["partType"] == "He"){ cut0 = 60.; cut1 = 80.; cut3 = 300;}
  else{  cut0 = 15.; cut1 = 20.; cut3 = 100.;}

  double maxV     = 45.;    // limit of range in V for calibration
  double maxT     = 170.;   // limit of range in T for calibration
  double deltaInt = 20.;    // intercept and slope to calculate maximum difference in V or T vs stage number
  double deltaSlope = 30.;
  //double brickC     = 6.5;  // Distance in t from wedge for center of range to take protons through bricks only
  float mxEstage;

  if (config.item_str["partType"] == "He") mxEstage = 400.;
  else  mxEstage = 100.;
  double diff = 0;
  int N = 0;
  for (int EvtNum = 0; EvtNum < theEvtRecon->nEvents; ++EvtNum) {
    Event thisEvent;
    thisEvent = theEvtRecon->evtList[EvtNum];

    Tf[0] = thisEvent.Thit[0]; Tf[1] = thisEvent.Thit[1];// Front tracker coordinates
    Vf[0] = thisEvent.Vhit[0]; Vf[1] = thisEvent.Vhit[1];
    T[0]  = thisEvent.Thit[2];  T[1] = thisEvent.Thit[3];// Rear tracker coordinates
    V[0]  = thisEvent.Vhit[2];  V[1] = thisEvent.Vhit[3];

    float eStage[nStage]; // energy in this stage
    float eTot = 0;
    bool good = true;
    for (int stage = 0; stage < nStage; ++stage) {
      // Extrapolate the rear track vector to the energy detector stage
      double Tcorr = theGeometry->extrap2D(Ut, T, theGeometry->energyDetectorU(stage));
      double Vcorr = theGeometry->extrap2D(Uv, V, theGeometry->energyDetectorU(stage));

      // Apply the pedestals and TV correction to the stage ADC values
      bool inBounds;
      eStage[stage] = theEvtRecon->GainFac[stage] * theTVcorr->corrFactor(stage, Tcorr, Vcorr, inBounds) * (thisEvent.ADC[stage] - theEvtRecon->Peds[stage]);
      if (!inBounds){// || eStage[stage] > mxEstage) {
        good = false;
        break;
      }
      if (abs(Vcorr) > 40.0 || abs(Tcorr) > 160.0) { // To avoid protons that scatter out of the device
        good = false;
        break;
      } 
      eTot += eStage[stage];
    }
    if (!good) continue;
    // Skip events not compatible with parameterized energy deposit
    // (Lennart Volz, December 2017)
    // bool dropEvt = pCTcut->EnrgCut(eStage, eTot, cut0, cut1,cut3);
    // ENERGY CUTS FOR REMOVING EVENT PILE-UP AND BAD EVENTS
    // if (dropEvt) continue;

    // calculate geometric WET "wt0" for the wedge phantom using tracker info.
    double Uin  = Ust[0];

    // First estimate of the u coordinate at entry into the phantom wedge
    double Tin  = theGeometry->extrap2D(Uft, Tf, Uin);
    double Vin  = theGeometry->extrap2D(Ufv, Vf, Uin);

    // First brick past the wedge
    double TinB = theGeometry->extrap2D(Uft, Tf, Ust[1]); 
    double VinB = theGeometry->extrap2D(Ufv, Vf, Ust[1]);
    
    if (fabs(Vin) > maxV || fabs(Tin) > maxT) continue; // if track is outside of detector in T/V - skip it

    //Back of bricks
    double Tout = theGeometry->extrap2D(Ut, T, Uout); 
    double Vout = theGeometry->extrap2D(Uv, V, Uout);
    if (fabs(Vout) > maxV || fabs(Tout) > maxT) continue;  // if track is outside of detector in T/V- skip it
        
    bool emptyEvt = false;
    int cs = 0;
    
    if ((Tin > Tw2 && Tin < Tw3) || (TinB < Tw1) || (TinB > Tw4))
      {// ONLY WEDGE SLOPE TO BE USED EVERYTHING ELSE CREATES ISSUES!
	if(Tin < Tw3) continue;
	if (TinB > tBrickEnd) emptyEvt = true;	
	cs = 1;
	continue;
      }

    else if (TinB >= Tw1 && Tin <= Tw2) { // in the -t wedge of the phantom; Uin
                                          // and Tin get overwritten here
      if (getLineIntersection(theEvtRecon->uhitT[0], thisEvent.Thit[0], theEvtRecon->uhitT[1], thisEvent.Thit[1],
                              Uin + BrickThickness, Tw1, Uin, Tw2, Uin, Tin)) {
        Vin = theGeometry->extrap2D(Ufv, Vf, Uin);
        cs = 4;
      } else
        continue;
    } else if (Tin >= Tw3 && TinB <= Tw4) { //&& Tout<tBrickEnd) {    // in the +t wedge of
      // the phantom; Uin and Tin get overwritten here
      // ALSO EXCLUDE PARTICLES THAT LEFT THE BRICKS
      if (getLineIntersection(theEvtRecon->uhitT[0], thisEvent.Thit[0], theEvtRecon->uhitT[1], thisEvent.Thit[1], Uin, Tw3,
                              Uin + BrickThickness, Tw4, Uin, Tin)) {
        Vin = theGeometry->extrap2D(Ufv, Vf, Uin);
        // if (config.item_int["Nbricks"] == 0 || (Tw4-TinB)>0.0) {  // stay away from the
        // large discontinuity at edges of bricks??
        cs = 5;
        //} else continue;
      } else
        continue;
    } else {
      cs = 2;
      continue;
    }
    // calculate polystyrene thickness
    double Length;
    if (emptyEvt) {
      Length = 0.0;
      hTotEcorr->entry(eTot);
      hTotEcorr_root->Fill(eTot);
      for (int stage = 0; stage < nStage; ++stage){
        hStgEcorr[stage]->entry(eStage[stage]);
	hStgEcorr_root[stage]->Fill(eStage[stage]);
      }
    }
    else {
      Vin = theGeometry->extrap2D(Ufv, Vf, Uin);
      // if (fabs(Vin-Vout) > deltaInt +
      // deltaSlope*sqrt(float(config.item_int["Nbricks"]))) continue;   // Too much scatter;
      // skip the event //NO SCATTERING FILTER NEEDED
      Length = sqrt( pow((Tin - Tout),2) + pow((Vin - Vout),2) + pow((Uin - Uout),2) ); // Length  crossed (converted to WET in Wepl.h)
    }
    hLengthProf->Fill(TinB,Length);
    // For each stage where the proton stops, increment the corresponding range-energy table cell content.
    // Energy cuts are now performed in seperate funcion (Last modified: Lenny, December 2017)
    
    if (eStage[4] > config.item_float["thr4"]) { // Particles end in stage 4
      if (cs > 3) {
        REhist[4]->entry(Length, eStage[4]);
	REhist_root[4]->Fill(Length, eStage[4]);
	dEEhist[4]->Fill(eStage[3], eStage[4]);
      }
    }
    
    else if (eStage[3] > config.item_float["thr3"]) { // Particles end in stage 3
      if (cs > 3) {
        REhist[3]->entry(Length, eStage[3]);
	REhist_root[3]->Fill(Length, eStage[3]);
	dEEhist[3]->Fill(eStage[2], eStage[3]);
      }
    }

    else if (eStage[2] > config.item_float["thr2"] ) { // Particles end in stage 2
      if (cs > 3) {
        REhist[2]->entry(Length, eStage[2]);
	REhist_root[2]->Fill(Length, eStage[2]);
	dEEhist[2]->Fill(eStage[1], eStage[2]);
      }
    }
    else if (eStage[1] > config.item_float["thr1"]) { // Particles end in stage 1
      if (cs > 3) {
        REhist[1]->entry(Length, eStage[1]);
	REhist_root[1]->Fill(Length, eStage[1]);
	dEEhist[1]->Fill(eStage[0], eStage[1]);
      }
    }
    else if (eStage[0] > config.item_float["thr0"]) { // Particles end in stage 0
	if (cs > 3){
	REhist_root[0]->Fill(Length, eStage[0]);
        REhist[0]->entry(Length, eStage[0]);
	}
    }
  } // End of the event loop

  pCTcalibRootFile->mkdir("runCorrEnrgs");
  pCTcalibRootFile->cd("runCorrEnrgs");
  for (int stage = 0; stage < nStage; ++stage) hStgEcorr_root[stage]->Write("",TObject::kOverwrite);
  hTotEcorr_root->Write("",TObject::kOverwrite);

  pCTcalibRootFile->cd();
  pCTcalibRootFile->mkdir("Phantom");  
  pCTcalibRootFile->cd("Phantom");
  hLengthProf->Write("",TObject::kOverwrite);
  pCTcalibRootFile->cd();
  // memory management
  delete (hTotEcorr);
  delete (theCalibration);
  for (int stage = 0; stage < nStage; ++stage) delete (hStgEcorr[stage]);
}

//////////////////////////////////////////////////////////////////////
// 2D equation to find line intersection
//////////////////////////////////////////////////////////////////////
bool pCTcalib::getLineIntersection(double p0u, double p0t,                   // Two points on the first line
                         double p1u, double p1t, double p2u, double p2t,     // Two points on the second line
                         double p3u, double p3t, double &ipu, double &ipt) { // Return the intersection point
  if (p0u == p1u || p2u == p3u)
    return false;
  double slope1 = (p1t - p0t) / (p1u - p0u); // Slope of the first line
  double slope2 = (p3t - p2t) / (p3u - p2u); // Slope of the second line
  if (slope1 == slope2)
    return false;                   // There is no intersection
  double int1 = p0t - slope1 * p0u; // Intercept of the first line
  double int2 = p2t - slope2 * p2u; // Intercept of the second line

  ipu = (int1 - int2) / (slope2 - slope1);
  ipt = (slope2 * int1 - slope1 * int2) / (slope2 - slope1);
  return true;
}
