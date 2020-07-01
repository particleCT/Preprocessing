// Code to process raw calibration data to derive the TV and WEPL calibration constants.
#include "pCTcalib.h"
pCTcalib::~pCTcalib() {pCTcalibRootFile->Close(); }
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
pCTcalib::pCTcalib(string inputFileName)
{

  theConfig = pCTconfig::GetInstance();
  
  cout << "\nEntering pCTcalib, the constructor for the pCT calibration "<<endl;    
  cout << "Filename for the list of calibration data files = " << inputFileName << endl;
  cout << "The output directory is " << theConfig->item_str["outputDir"] << endl;
  cout << "The maximum number of events is " << theConfig->item_int["max_events"] << endl;
  cout << "The maximum time is " << theConfig->item_int["max_time"] << endl;
  cout << "The number of debug event printouts is " << theConfig->item_int["n_debug"] << endl;
  cout << "The number of events to plot is " << theConfig->item_int["n_plot"] << endl;
  cout << "The input WEPL calibration file is " << theConfig->item_str["Wcalib"] << endl;
  cout << "The input TV correction calibration file is " << theConfig->item_str["TVcorr"] << endl;
  cout << "Input date range = " << theConfig->item_int["minDate"] << " to " << theConfig->item_int["maxDate"] << endl;
  cout << "The input run range is " << theConfig->item_int["minRun"] << " to " << theConfig->item_int["maxRun"] << endl;
  cout << "The calibration wedge offset is assumed to be " << theConfig->item_float["wedgeoffset"] << endl;
  cout << "The particle type is assumed to be " << theConfig->item_str["partType"] << endl;
  cout << "Real-time pedestal and gain recalibration is done? " << theConfig->item_int["recalibrate"] << endl;
  for (int stage = 0; stage < nStage; ++stage) {
    cout << "  For stage " << stage << " the range in which to look for pedestals begins at " << theConfig->item_int[Form("pedrng%d",stage)] << " ADC counts " << endl;
  }
  CalFile = inputFileName;
  if (theConfig->item_str["partType"] == "H") EnergyBinWidth = 0.5;
  else EnergyBinWidth = 1.0;
  RangeBinWidth = 1.0;

  if (theConfig->item_str["partType"] == "H") {
    EG4stage[0] = 25.25; EG4stage[1] = 28.01; EG4stage[2] = 32.76;// MC derived stage energies, used to calibrate to MeV
      EG4stage[3] = 42.62; EG4stage[4] = 67.71;                     // for protons*/
      TVnormalizeFactor = 50;
  } else {
    EG4stage[0] = 100.; EG4stage[1] = 111.; EG4stage[2] = 129.;// MC derived stage energies, used to calibrate to MeV for He
      EG4stage[3] = 166.; EG4stage[4] = 215.;
      TVnormalizeFactor = 200;
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

  cout << "The beam particle type is " << theConfig->item_str["partType"] << endl;
  float EG4tot = 0.;
  for (int i = 0; i < nStage; ++i) {
    cout << "MC derived stage energy for an empty event, for stage " << i << " is " << EG4stage[i] << endl;
    EG4tot += EG4stage[i];
  }
  cout << "MC derived total energy for an empty event is " << EG4tot << endl << endl;
  currentTime = time(NULL);
  now = localtime(&currentTime);  



  // Initialize class
  theGeometry    = new pCTgeo(theConfig->item_float["wedgeoffset"]);                 
  theCuts        = new pCTcut();
  theTVcorr      = new TVcorrection(pCTcalibRootFile, 1);
  theEvtRecon    = new EvtRecon();

  float wedgeLimit = theGeometry->getTWedgeBreaks(4) + 25.0; // NO BRICK OFFSET //+5.
  float openRange  = 20.0;
  float pedestals[nStage] = {0};
  int pedMin[nStage];  
  for (int stage = 0; stage < nStage; stage++) pedMin[stage]   = theConfig->item_int[Form("pedrng%d",stage)];
  ProfileE_Tot      = new TH2D("ProfileE_Tot","",300,-100,150, 1000, 0, 1000);//
  ProfileE_Tot_f    = new TH2D("ProfileE_Tot_f","",600,-100,150, 1000, 0, 1000);//
  theCalibration    = new pedGainCalib(pCTcalibRootFile, pedMin, pedestals,-150., -151., wedgeLimit, wedgeLimit + openRange);  

}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int pCTcalib::TVmapper() {
  if (calFileNames.size() == 0) {
    cout << "TVmaper: only " << calFileNames.size() << " calibration file names found.  Need at least 1 for the TV map." << endl;
    return -2;
  }
  theConfig->item_int["doGains"] = false; // This is for the empty one for the TV correction
  theConfig->item_int["Nbricks"] = 0;

  // Calibration
  //TH2D* dEETest[nStage];
  //for (int stage = 0; stage < nStage; stage++) //dEETest[stage] = new TH2D(Form("TestdEE_%d",stage),"", 400, 0, 300, 400, 0, 300);

  
  float bin0[nStage] = { 1000., 1000., 1000., 1000., 1000. };

  //---------------------------------------------------------------------------
  // Alignement
  //--------------------------------------------------------------------------- 
  TH1D* TalignementPos01 = new TH1D("TalignementPosT0->T1","",300, -20, -20 +300*0.1);
  TH1D* TalignementPos23 = new TH1D("TalignementPosT2->T3","",300, -20, -20 +300*0.1);

  //Forward Project
  TH1D* TalignementPos12 = new TH1D("TalignementPosT1->T2","",300, -20, -20 +300*0.1);
  TH1D* ValignementPos12 = new TH1D("ValignementPosV1->V2","",300, -20, -20 +300*0.1);
  TH1D* TalignementPos13 = new TH1D("TalignementPosT1->T3","",300, -20, -20 +300*0.1);
  TH1D* ValignementPos13 = new TH1D("ValignementPosV1->V3","",300, -20, -20 +300*0.1);
  
  // Backward project
  TH1D* TalignementPos20 = new TH1D("TalignementPosT2->T0","",300, -20, -20 +300*0.1);
  TH1D* ValignementPos20 = new TH1D("ValignementPosV2->V0","",300, -20, -20 +300*0.1);
  TH1D* TalignementPos21 = new TH1D("TalignementPosT2->T1","",300, -20, -20 +300*0.1);
  TH1D* ValignementPos21 = new TH1D("ValignementPosV2->V1","",300, -20, -20 +300*0.1);

  // Direction
  TH1D* TalignementDir01 = new TH1D("TalignementDirT0T1","",300, -.1, -.1 +200*0.001);
  TH1D* TalignementDir23 = new TH1D("TalignementDirT2T3","",300, -.1, -.1 +200*0.001);
  TH1D* TalignementDir12 = new TH1D("TalignementDirT1T2","",300, -.1, -.1 +200*0.001);
  TH1D* ValignementDir01 = new TH1D("ValignementDirV0V1","",300, -20, -20 +300*0.1);
  TH1D* ValignementDir23 = new TH1D("ValignementDirV2V3","",300, -20, -20 +300*0.1);
  
  // On trackers
  TH1D* TalignementPos0 = new TH1D("TalignementPosT0","",300, -200, -200 +300*1.0);
  TH1D* TalignementPos1 = new TH1D("TalignementPosT1","",300, -200, -200 +300*1.0);
  TH1D* TalignementPos2 = new TH1D("TalignementPosT2","",300, -200, -200 +300*1.0);
  TH1D* TalignementPos3 = new TH1D("TalignementPosT3","",300, -200, -200 +300*1.0);

  TH1D* ValignementPos0 = new TH1D("ValignementPosV0","",300, -200, -200 +300*1.0);
  TH1D* ValignementPos1 = new TH1D("ValignementPosV1","",300, -200, -200 +300*1.0);
  TH1D* ValignementPos2 = new TH1D("ValignementPosV2","",300, -200, -200 +300*1.0);
  TH1D* ValignementPos3 = new TH1D("ValignementPosV3","",300, -200, -200 +300*1.0);

  //---------------------------------------------------------------------------
  // Source Position in U
  //---------------------------------------------------------------------------
  // V source position
  TH1D* Vsource = new TH1D("Vsource","", 900, -10000.,-1000.); //Between 1 and 10 m in 1cm precision
 
  // T source position
  TH1D* Tsource = new TH1D("Tsource","", 900, -10000.,-1000.);


  cout<<"File:"<<calFileNames[0]<<endl;
  theConfig->item_str["inputFileName"] = calFileNames[0];
  //theCalibration.config = theConfig;
  theCalibration->ResetHist();
  theEvtRecon->ReadInputFile(theGeometry, theTVcorr,  theConfig->item_str["inputFileName"], theCalibration); // Pedestal are determined here
  theCalibration->WriteHist(); // For analysis

  // U Position of the trackers
  Uft[0] = theEvtRecon->uhitT[0]; Uft[1] = theEvtRecon->uhitT[1];
  Ufv[0] = theEvtRecon->uhitV[0]; Ufv[1] = theEvtRecon->uhitV[1];
  Ut[0]  = theEvtRecon->uhitT[2];  Ut[1] = theEvtRecon->uhitT[3];
  Uv[0]  = theEvtRecon->uhitV[2];  Uv[1] = theEvtRecon->uhitV[3];
    
  cout << "pCTcalib: U values for T coordinates = " << Ut[0] << " and " << Ut[1] << endl;
  cout << "pCTcalib: U values for V coordinates = " << Uv[0] << " and " << Uv[1] << endl;  
  for (int stage = 0; stage <nStage; stage++) cout << "Pedestal measured for stage " << stage << " is " << theEvtRecon->Peds[stage] << " ADC counts " << endl;
  Int_t binx,biny,binz;
  // Define a signal histogram for each pixel
  for (int stage = 0; stage < nStage; stage++) {
    for (int pix = 0; pix < nPix ; pix++) { // 480 with the underflow overflow
      theTVcorr->TVcorrHist[stage]->GetBinXYZ(pix, binx, biny, binz);
      string Title = "Stage " + to_string((long long int)stage) + " Global " + to_string((long long int)pix) +
	" BinX " +  to_string((long long int)binx) +  "BinY" + to_string((long long int)biny);
      if(theConfig->item_str["partType"] == "H") pxHistADC[stage][pix] = new TH1D(Title.c_str(),Title.c_str(), 2000, bin0[stage], bin0[stage] +4.*2000);
      else pxHistADC[stage][pix] = new TH1D(Title.c_str(), Title.c_str(), 2000, bin0[stage], bin0[stage] +8.*2000);
    }
  }
  
  TH2D* hTVmap = new TH2D("Profile of tracks at energy detector","",150, -150, 150 , 80, -40, 40);
  cout << "TVmapper: starting event loop" << endl;
  for (int EvtNum = 0; EvtNum < theEvtRecon->nEvents; EvtNum++) {
    Event thisEvent;
    thisEvent = theEvtRecon->evtList[EvtNum];
    // Front tracker T-V coordinates
    Tf[0] = thisEvent.Thit[0]; Tf[1] = thisEvent.Thit[1];
    Vf[0] = thisEvent.Vhit[0]; Vf[1] = thisEvent.Vhit[1];
    T[0] = thisEvent.Thit[2]; T[1] = thisEvent.Thit[3];
    V[0] = thisEvent.Vhit[2]; V[1] = thisEvent.Vhit[3];
    
    if (EvtNum % 1000000 == 0) {
      cout << "  Processing event " << EvtNum << endl;
      theEvtRecon->dumpEvt(thisEvent);
    }
    hTVmap->Fill(T[1], V[1]);
    // Extrapolate the rear track vector to the energy detector stage
    for (int stage = 0; stage < nStage; stage++) {
      
      float Tcorr = theGeometry->extrap2D(Ut, T, theGeometry->energyDetectorU(stage));
      float Vcorr = theGeometry->extrap2D(Uv, V, theGeometry->energyDetectorU(stage));
      // Outside of the detector
      if (fabs(Tcorr) > 150.0) continue;
      if (fabs(Vcorr) > 40.0)  continue;
      int global = theTVcorr->TVcorrHist[0]->FindBin(Tcorr,Vcorr);
      pxHistADC[stage][global]->Fill(thisEvent.ADC[stage] - theEvtRecon->Peds[stage]);
      //pxHistADC[stage][global]->Fill(thisEvent.ADC[stage]);
    }


    TalignementPos01->Fill(Tf[1]-Tf[0]);
    TalignementPos23->Fill(T[1]-T[0]);


    // Projected to the exit tracker
    float Tin12 = theGeometry->extrap2D(Uft,Tf, Ut[0]);
    float Vin12 = theGeometry->extrap2D(Ufv,Vf, Uv[0]);
    float Tin13 = theGeometry->extrap2D(Uft,Tf, Ut[1]);
    float Vin13 = theGeometry->extrap2D(Ufv,Vf, Uv[1]);
    TalignementPos12->Fill(Tin12 - T[0]);
    ValignementPos12->Fill(Vin12 - V[0]);
    TalignementPos13->Fill(Tin13 - T[1]);
    ValignementPos13->Fill(Vin13 - V[1]);
    
    // Projected to the front tracker
    float Tout20 = theGeometry->extrap2D(Ut,T, Uft[0]);
    float Vout20 = theGeometry->extrap2D(Uv,V, Ufv[0]);
    float Tout21 = theGeometry->extrap2D(Ut,T, Uft[1]);
    float Vout21 = theGeometry->extrap2D(Uv,V, Ufv[1]);
    TalignementPos20->Fill(Tout20 - Tf[0]);
    ValignementPos20->Fill(Vout20 - Vf[0]);
    TalignementPos21->Fill(Tout21 - Tf[1]);
    ValignementPos21->Fill(Vout21 - Vf[1]);
    
    // Direction
    TalignementDir01->Fill((Tf[1]-Tf[0])/(Uft[1]-Uft[0]));
    TalignementDir23->Fill((T[1]-T[0])/(Ut[1]-Ut[0]));
    TalignementDir12->Fill((T[2]-Tf[1])/(Ut[0]-Uft[1]));
    ValignementDir01->Fill(Vf[1]-Vf[0]);
    ValignementDir23->Fill(V[1]-V[0]);
    
    // On the tracker
    TalignementPos0->Fill(Tf[0]); TalignementPos1->Fill(Tf[1]);      
    TalignementPos2->Fill(T[0]);  TalignementPos3->Fill(T[1]);
    ValignementPos0->Fill(Vf[0]); ValignementPos1->Fill(Vf[1]);      
    ValignementPos2->Fill(V[0]);  ValignementPos3->Fill(V[1]);      

    //Source position alignment V
    double dirV = (Vf[1] - Vf[0])/(Ufv[1]-Ufv[0]); //50mm is the distance between tracker planes
    double dirT = (Tf[1] - Tf[0])/(Uft[1]-Uft[0]); 
    Vsource->Fill(-(Vf[0]/dirV));
    Tsource->Fill(-(Tf[0]/dirT));
          
  }
  pCTcalibRootFile->cd("");
  pCTcalibRootFile->mkdir("tracksProfile");
  pCTcalibRootFile->cd("tracksProfile");
  hTVmap->Write("");
  
  pCTcalibRootFile->cd();
  pCTcalibRootFile->mkdir("Alignement");  
  pCTcalibRootFile->cd("Alignement");

  TalignementPos01->Write("",TObject::kOverwrite);
  TalignementPos23->Write("",TObject::kOverwrite);
  
  TalignementPos12->Write("",TObject::kOverwrite);
  ValignementPos12->Write("",TObject::kOverwrite);
  TalignementPos13->Write("",TObject::kOverwrite);
  ValignementPos13->Write("",TObject::kOverwrite);
  
  TalignementPos20->Write("",TObject::kOverwrite);
  ValignementPos20->Write("",TObject::kOverwrite);
  TalignementPos21->Write("",TObject::kOverwrite);
  ValignementPos21->Write("",TObject::kOverwrite);
  
  TalignementDir01->Write("",TObject::kOverwrite);
  TalignementDir23->Write("",TObject::kOverwrite);
  TalignementDir12->Write("",TObject::kOverwrite);
  ValignementDir01->Write("",TObject::kOverwrite);
  ValignementDir23->Write("",TObject::kOverwrite);
  
  TalignementPos0->Write("",TObject::kOverwrite);
  TalignementPos1->Write("",TObject::kOverwrite);
  TalignementPos2->Write("",TObject::kOverwrite);
  TalignementPos3->Write("",TObject::kOverwrite);
  ValignementPos0->Write("",TObject::kOverwrite);
  ValignementPos1->Write("",TObject::kOverwrite);
  ValignementPos2->Write("",TObject::kOverwrite);
  ValignementPos3->Write("",TObject::kOverwrite);
  
  Vsource->Write("",TObject::kOverwrite);
  Tsource->Write("",TObject::kOverwrite);

  // Save plots of the histograms
  pCTcalibRootFile->mkdir("ADC_Stage");
  pCTcalibRootFile->cd("ADC_Stage");
  cout << "TVmapper: Starting evaluation of the TV correction factors." << endl;
  
  //Fill the TV corrs
  for (int pix = 0; pix < nPix; pix++) {
    cout << "TVmapper pixel stage ADC values:  " << pix;
    for (int stage = 0; stage < nStage; stage++) {
      pxHistADC[stage][pix]->Write("");
      float xLow, xHigh;
      float xpeak, xmax, max;
      if(pxHistADC[stage][pix]->GetEntries()>10){
	max   =  pxHistADC[stage][pix]->GetMaximum();
	xmax  =  pxHistADC[stage][pix]->GetBinCenter(pxHistADC[stage][pix]->GetMaximumBin()); 
	TF1 *f1 = new TF1("f1", "gaus", xmax-100, xmax+100);
	f1->SetParameter(0, max);
	f1->SetParameter(1, xmax);
	Int_t fitStatus = pxHistADC[stage][pix]->Fit(f1,"QR"); // Q means quiet, R is 
	if(fitStatus==0) // Everything passed
	  {
	    xpeak  = f1->GetParameter(1);// mean
	    xHigh  = xpeak + f1->GetParameter(2); //sigma
	    xLow   = xpeak - f1->GetParameter(2); //sigma
	    theTVcorr->TVcorrHist[stage]->SetBinContent(pix, EG4stage[stage]/ xpeak);
	  }
	else theTVcorr->TVcorrHist[stage]->SetBinContent(pix, 0.99);
      }
      else theTVcorr->TVcorrHist[stage]->SetBinContent(pix, 0.99);
      cout << " stg" << stage << "=       "<< theTVcorr->TVcorrHist[stage]->GetBinContent(pix);
      delete pxHistADC[stage][pix]; // Free up the memory that was sucked up by the pixel histograms
      
    }
    cout<<endl;
  }
  
  pCTcalibRootFile->cd("");
  pCTcalibRootFile->mkdir("TVcorrProfile");
  pCTcalibRootFile->cd("TVcorrProfile");
  
  for(int stage =0; stage<nStage; stage++) theTVcorr->TVcorrHist[stage]->Write("", TObject::kOverwrite);  
  cout << "TVmapper: Finished with mapping the TV calibration" << endl;
  return 0;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void pCTcalib::enrgDep() { // Calculate energy depositions in each stage and their sum using new TV maps, to check that calibration worked
  // Create Histogram  for each stage
  
  for (int stage = 0; stage < nStage; ++stage) {
    string Title = "pCTcalib::enrgDep Energy of stage " + to_string((long long int)stage);
    if (theConfig->item_str["partType"] == "H") stgHistE[stage]  = new TH1D(Title.c_str(), Title.c_str(),  800, 15, 15 + 800*0.175); 
    else stgHistE[stage]  = new TH1D(Title.c_str(), Title.c_str(),  800, 15, 15 + 800*0.9); 
  }
  // Sum of stage energies
  if (theConfig->item_str["partType"] == "H") EsumH = new TH1D("pCTcalib::enrgDep Sum of stage energies","pCTcalib::enrgDep Sum of stage energies", 900, 0, 0.25*900);
  else EsumH = new TH1D("pCTcalib::enrgDep Sum of stage energies", "pCTcalib::enrgDep Sum of stage energies", 900, 0, 1.0*900);
  
  // Create Histogram for each pixel
  for (int iPix = 0; iPix < nPix; iPix++) {
    for (int stage = 0; stage < nStage; ++stage) {
      string Title = "pCTcalib::enrgDep Energy of stage " + to_string((long long int)stage) + " for pixel " +
                     to_string((long long int)iPix);
      if (theConfig->item_str["partType"] == "H") pxHistE[stage][iPix] = new TH1D(Title.c_str(), Title.c_str(), 200, 15., 15 +200*0.35);
      else pxHistE[stage][iPix] = new TH1D(Title.c_str(), Title.c_str(), 200, 60., 60 +200*1.4);
    }
  }
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
      pxHistE[stage][iPix]->Fill(Ecorr);
      stgHistE[stage]->Fill(Ecorr);
      Esum += Ecorr;
    }
    EsumH->Fill(Esum);
  }
  pCTcalibRootFile->mkdir("EnergyHistogram");  
  pCTcalibRootFile->cd("EnergyHistogram");
  for(int i =0; i<5; i++) stgHistE[i]->Write("", TObject::kOverwrite);
  EsumH->Write("TotalEnergy",TObject::kOverwrite);
  pCTcalibRootFile->cd();

  for (int stage = 0; stage < nStage; ++stage) {
    float Sadc, xmax, max;
    if(stgHistE[stage]->GetEntries()>10){
      max   =  stgHistE[stage]->GetMaximum();
      xmax  =  stgHistE[stage]->GetBinCenter(stgHistE[stage]->GetMaximumBin());
      TF1 *f1 = new TF1("f1", "gaus", xmax-10, xmax+10);
      f1->SetParameter(0, max);
      f1->SetParameter(1, xmax);
      Int_t fitStatus = stgHistE[stage]->Fit(f1,"QR"); // Q means quiet, R is    
      if(fitStatus==0) Sadc = f1->GetParameter(1);// mean // Everything passed
      else { Sadc = 0.;
      cout << "pCTcalib::enrgDep, no data for stage " << stage << endl;
      }
    }
    else Sadc = 0;
    Est[stage] = Sadc;

    cout << "pCTcalib::enrgDep, Energy in stage " << stage << " = " << Est[stage] << "; dE=" << Est[stage] - EG4stage[stage] << " MeV" << endl;
    delete stgHistE[stage];
  }
  // Using the last empty run to calculate gain recalibration and pedestals
  for (int stage = 0; stage < nStage; ++stage) {
    theTVcorr->ped[stage]   = theEvtRecon->Peds[stage];
    theTVcorr->Eempt[stage] = Est[stage];
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
  if (calFileNames.size() < 6) {
    cout << "Only " << calFileNames.size() << " calibration raw data file names found.  Need all 6 to include "
      "the WEPL calibration." << endl;
    return -1;
  }
  for (int bricks = 0; bricks < nBricks; bricks++) cout << "The raw data file for " << bricks << " bricks is " << calFileNames[bricks+1] << endl; //bug
  TH2D* dEEhist[nStage];
  TH2D* REhist[nBricks][nStage];

  for (int stage = 0; stage < nStage; ++stage) {
    dEEhist[stage] = new TH2D(Form("dE-EStage_%d",stage), Form("dE-E spectra for stage %d", stage), // Title
			      nEnrg, 0., nEnrg * EnergyBinWidth, nEnrg, 0., nEnrg * EnergyBinWidth);
    for (int bricks = 0; bricks < nBricks; ++bricks) {
      float enrgB0 = 0.; // The analysis of energy slices further down assumes that the first energy bin starts at zero
      string Title = "WEPL calibration array for stage " + to_string((long long int)stage) + " and nBricks=" + to_string((long long int)bricks);
      REhist[bricks][stage] = new TH2D(Form("WEPLcalib_stage%d_bricks%d",stage,bricks), Form("WEPLcalib_stage%d_bricks%d",stage,bricks),
					     nRange, 0, RangeBinWidth*nRange, nEnrg, enrgB0, enrgB0 + nEnrg*EnergyBinWidth);
      REhist[bricks][stage]->GetXaxis()->SetTitle("Range (mm)");
      REhist[bricks][stage]->GetYaxis()->SetTitle("Energy (MeV)");      
    }
  }

  TH2D* REhist_Tot = new TH2D("RE_Tot","Range-Energy Full calibration",nRange,0,RangeBinWidth*nRange, nEnrg*5, 0, nEnrg * 5*EnergyBinWidth);  
  for (int bricks = 0; bricks < nBricks; bricks++) {
    theConfig->item_str["inputFileName"] = calFileNames[bricks+1];
    theConfig->item_int["Nbricks"] = bricks;
    FilldEE(dEEhist);
  }
  for (int stage = 0; stage < nStage; stage++) theCuts->CalculatedEEFilterParameters(dEEhist[stage],dEElow[stage],dEEhigh[stage],stage);    

  // Submit the event processing  -- fill the histograms
  for (int bricks = 0; bricks < nBricks; bricks++) {
    theConfig->item_str["inputFileName"] = calFileNames[bricks+1]; // Without flat bricks
    theConfig->item_int["Nbricks"] = bricks;
    procWEPLcal(REhist[bricks] , dEEhist, REhist_Tot);
  }
  cout << "pCTcalib::Wcalib, begin adding together the range-energy tables from the runs with different numbers of bricks." << endl;
  // Combine the resulting maps into one map for each stage by adding them for all bricks
  pCTcalibRootFile->mkdir("dEE");
  pCTcalibRootFile->mkdir("REhist");
  for (int stage = 0; stage < nStage; ++stage) {
    // Save the individual histograms
    pCTcalibRootFile->cd("dEE");
    dEEhist[stage]->Write("",TObject::kOverwrite);
    pCTcalibRootFile->cd("REhist");
    for (int bricks = 0; bricks < nBricks; ++bricks) REhist[bricks][stage]->Write("",TObject::kOverwrite);
    REhist_Tot->Write("",TObject::kOverwrite);

    //Add them up together
    for (int bricks = 1; bricks < nBricks; ++bricks) {
      string Title = "Summed WEPL calibration array for stage " + to_string((long long int)stage);
      REhist[0][stage]->Add(REhist[bricks][stage]);
      
    }
    // Save the added histograms
    pCTcalibRootFile->cd("REhist");
    REhist[0][stage]->Write(Form("RE_Tot_stage%d",stage),TObject::kOverwrite);
  }
  ProfileE_Tot->Write("",TObject::kOverwrite);
  ProfileE_Tot_f->Write("",TObject::kOverwrite);
  
  // Analyze energy slices of the summed maps, and plot the slices  stage by stage --Most likely range for an energy
  /*for (int stage = 0; stage < nStage; ++stage) {
    for (int nrg = 0; nrg < nEnrg; nrg++) { 
      float xLow, xHigh, Sadc,xmax, max;
      TH1D* RESlice =  REhist[0][stage]->ProjectionX("ProjX_Stage",nrg,nrg+1);
      Int_t NEntries= RESlice->GetEntries();
      if(NEntries>10){
      max   =  RESlice->GetMaximum();
      xmax  =  RESlice->GetBinCenter(RESlice->GetMaximumBin());
      TF1 *f1 = new TF1("f1", "gaus", xmax-10, xmax+10);
      f1->SetParameter(0, max);
      f1->SetParameter(1, xmax);
      Int_t fitStatus = RESlice->Fit(f1,"QR"); // Q means quiet, R is
      if(fitStatus==0) // Everything passed
	{
	  Sadc      = f1->GetParameter(1);// mean
	  xHigh     = Sadc + f1->GetParameter(2); //sigma
	  xLow      = Sadc - f1->GetParameter(2); //sigma
	}
      else if (NEntries>200){
	Sadc        =  RESlice->GetBinCenter(RESlice->GetMaximumBin()); // cheesy fix because of high level of low energy data
	xHigh       =  RESlice->GetBinCenter(RESlice->FindLastBinAbove(  RESlice->GetMaximum()/2));
	xLow        =  Sadc - (xHigh -xmax);
	}
      else{Sadc  = 0;	 xLow  = 0;  xHigh = 0;}
      }
      else{Sadc  = 0;	 xLow  = 0;  xHigh = 0;}
      delete RESlice;
      Rst[stage][nrg] = Sadc;
      Sst[stage][nrg] = (xHigh - xLow) / 2.2;      
    }
  }*/

  //Fit all stages by Range Projection -- Most likely energy for a given range
  for (int stage = 0; stage < nStage; ++stage) {
    for (int nrg = 0; nrg < nRange; nrg++) { 
      float xLow, xHigh, Sadc,xmax, max;
      TH1D* RESlice =  REhist[0][stage]->ProjectionY("ProjY_Stage",nrg,nrg+1);
      int bin_cut = RESlice->FindBin(5); // MeV
      RESlice->GetXaxis()->SetRange(bin_cut,RESlice->GetNbinsX());      
      Int_t NEntries= RESlice->GetEntries();
      if(NEntries>200){
      max   =  RESlice->GetMaximum();
      xmax  =  RESlice->GetBinCenter(RESlice->GetMaximumBin());
      TF1 *f1 = new TF1("f1", "gaus", xmax-10, xmax+10);
      f1->SetParameter(0, max);
      f1->SetParameter(1, xmax);
      Int_t fitStatus = RESlice->Fit(f1,"QR"); // Q means quiet, R is
      if(fitStatus==0) // Everything passed
	{
	  Sadc      = f1->GetParameter(1);// mean
	  xHigh     = Sadc + f1->GetParameter(2); //sigma
	  xLow      = Sadc - f1->GetParameter(2); //sigma
	}
      else{Sadc  = 0;	 xLow  = 0;  xHigh = 0;}
      }
      else{Sadc  = 0;	 xLow  = 0;  xHigh = 0;}
      delete RESlice;
      Rst[stage][nrg] = Sadc;
      Sst[stage][nrg] = (xHigh - xLow) / 2.2;      
    }
  }

  // Save the hist
  float rngB[nEnrg];
  TGraphErrors* rngEnrg[nStage];
  TGraphErrors* EnrgRng[nStage];
  pCTcalibRootFile->cd();
  pCTcalibRootFile->mkdir("CalibrationCurves");
  pCTcalibRootFile->cd("CalibrationCurves");

  for (int stage = 0; stage < nStage; ++stage) { // Normal direction
    string Title = "Range vs Energy Uncorrected for Stage " + to_string((long long int)stage);

    rngEnrg[stage] = new TGraphErrors(nEnrg);
    rngEnrg[stage]->SetTitle(Title.c_str());
    rngEnrg[stage]->SetName(Form("RangeVsEnergy_%d",stage));

    for (int nrg = 0; nrg < nEnrg; ++nrg){
      rngB[nrg] = float(nrg) * EnergyBinWidth;
      rngEnrg[stage]->SetPointError(nrg, 0.,   float(Sst[stage][nrg]));
    }
    rngEnrg[stage]->GetXaxis()->SetTitle("Energy (MeV)");
    rngEnrg[stage]->GetYaxis()->SetTitle("Range (mm)");      
    rngEnrg[stage]->Write("", TObject::kOverwrite);
  }
  
  // Define a region of validity for this brick
  float xmin = -13;
  float xmax = 42; // arbitrary;
  for (int stage = nStage-1 ; stage >= 0 ; --stage) { // X vs Y
    string Title = "Range vs Energy Uncorrected for Stage " + to_string((long long int)stage);
    rngEnrg[stage]    = new TGraphErrors();
    EnrgRng[stage] = new TGraphErrors();
    rngEnrg[stage]->SetTitle(Title.c_str());
    EnrgRng[stage]->SetTitle(Title.c_str());    
    rngEnrg[stage]->SetName(Form("RangeVsEnergy_%d",stage));
    EnrgRng[stage]->SetName(Form("EnergyVsRange_%d",stage));

    for (int nrg = 0; nrg < nRange; ++nrg){
      rngB[nrg] = float(nrg) * (RangeBinWidth) ;
      if(rngB[nrg] > xmin && rngB[nrg] < xmax) {
	rngB[nrg] += 0.5*RangeBinWidth;
	int N =  rngEnrg[stage]->GetN();
	EnrgRng[stage]->SetPoint(N, float(rngB[nrg]),float(Rst[stage][nrg]));	
	EnrgRng[stage]->SetPointError(N, 0.,float(Sst[stage][nrg]));
	rngEnrg[stage]->SetPoint(N, float(Rst[stage][nrg]), float(rngB[nrg]));
	rngEnrg[stage]->SetPointError(N, float(Sst[stage][nrg]), 0.);
      }
    }
    xmin += 50.8*1.03; // One brick shift
    xmax += 50.8*1.03;
    EnrgRng[stage]->GetXaxis()->SetTitle("Range (mm)");
    EnrgRng[stage]->GetYaxis()->SetTitle("Energy (MeV)");
    EnrgRng[stage]->Write("", TObject::kOverwrite);
    rngEnrg[stage]->GetXaxis()->SetTitle("Energy (MeV)");
    rngEnrg[stage]->GetYaxis()->SetTitle("Range (mm)");          
    rngEnrg[stage]->Write("", TObject::kOverwrite);    
  }

  // Combination of all stages
  TGraphErrors* RangeVsEnergy_Tot = new TGraphErrors(nEnrg);  
  TGraphErrors* EnergyVsRange_Tot = new TGraphErrors(nEnrg);
  RangeVsEnergy_Tot->SetName("RangeVsEnergy");
  EnergyVsRange_Tot->SetName("EnergyVsRange");
  
  for (int nrg = 0; nrg < nRange; nrg++) {
    float xLow, xHigh, Sadc,xmax, max;
    TH1D* RESlice =  REhist_Tot->ProjectionY("ProjY_Stage",nrg,nrg+1);
    int bin_cut = RESlice->FindBin(5); // MeV
    RESlice->GetXaxis()->SetRange(bin_cut,RESlice->GetNbinsX());
    Int_t NEntries= RESlice->GetEntries();
    if(NEntries>200){
      max   =  RESlice->GetMaximum();
      xmax  =  RESlice->GetBinCenter(RESlice->GetMaximumBin());
      TF1 *f1 = new TF1("f1", "gaus", xmax-10, xmax+10);
      f1->SetParameter(0, max);
      f1->SetParameter(1, xmax);
      Int_t fitStatus = RESlice->Fit(f1,"QR"); // Q means quiet, R is
      if(fitStatus==0) // Everything passed
	{
	  Sadc      = f1->GetParameter(1);// mean
	  xHigh     = Sadc + f1->GetParameter(2); //sigma
	  xLow      = Sadc - f1->GetParameter(2); //sigma
	}
      else{Sadc  = 0;    xLow  = 0;  xHigh = 0;}
    }
    else{Sadc  = 0;    xLow  = 0;  xHigh = 0;}
    delete RESlice;
    float rngB = REhist_Tot->GetXaxis()->GetBinCenter(nrg);
    float Sst  = (xHigh - xLow) / 2.2;
    
    int N =  RangeVsEnergy_Tot->GetN();
    RangeVsEnergy_Tot->SetPoint(N, rngB, Sadc);    
    RangeVsEnergy_Tot->SetPointError(N, 0.,Sst);
    EnergyVsRange_Tot->SetPoint(N, Sadc,rngB);    
    EnergyVsRange_Tot->SetPointError(N, Sst, 0.);
  }
  RangeVsEnergy_Tot->Write("", TObject::kOverwrite);
  EnergyVsRange_Tot->Write("", TObject::kOverwrite);
  return 0;
};    
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void pCTcalib::writeCalibfile() {

  // Dump all relevant info in the header
  pCTcalibRootFile->cd();
  TTree* header = new TTree("header","meta-data");
  header->Branch("particle",&theConfig->item_str["partType"]);
  header->Branch("minDate",&theConfig->item_str["minDate"]);
  header->Branch("maxDate",&theConfig->item_str["maxDate"]);
  header->Branch("minrun",&theConfig->item_int["minrun"], "minrun/I");
  header->Branch("maxrun",&theConfig->item_int["maxrun"], "maxrun/I" );
  header->Branch("outputDir",&theConfig->item_str["outputDir"]);
  for(int i=0; i<=5; i++) header->Branch(Form("calFileName_%d",i),&calFileNames[i]);
  for (int stage = 0; stage < nStage; ++stage) header->Branch(Form("peds_%d",stage),&theEvtRecon->Peds[stage]);
  for (int stage = 0; stage < nStage; ++stage) header->Branch(Form("Est_%d",stage),&Est[stage]);
  header->Branch("Esum",&EnS, "Esum/F");
  header->Branch("EnergyBinWidth",&EnergyBinWidth, "EnergyBinWidth/F");
  header->Branch("RangeBinWidth",&RangeBinWidth, "RangeBinWidth/F");
  header->Branch("Year",&now->tm_year + 1900, "Year/I");
  header->Branch("Month",&now->tm_mon + 1, "Month/I");
  header->Branch("Day",&now->tm_mday , "Day/I");
  header->Branch("Hour",&now->tm_hour , "Hour/I");
  header->Branch("Min",&now->tm_min , "Min/I");

  float MiddleWedge;
  header->Branch("MidlleWedge",&MiddleWedge , "MiddleWedge/F");
  float Tw2 = theGeometry->getTWedgeBreaks(2);
  float Tw3 = theGeometry->getTWedgeBreaks(3);
  
  MiddleWedge =  (Tw3 + Tw2)/2;

    
  for (int stage = 0; stage < nStage; ++stage) header->Branch(Form("thr_%d",stage),&theConfig->item_float[Form("thr%d",stage)],Form("thr_%d/F",stage));
  for (int stage = 0; stage < nStage; ++stage){
    for (int i= 0;  i < 3; ++i) header->Branch(Form("dEElow_Stage%d_%d",stage,i),&dEElow[stage][i]);
    for (int i= 0;  i < 3; ++i) header->Branch(Form("dEEhigh_Stage%d_%d",stage,i),&dEEhigh[stage][i]);
  }  
  header->Fill();
  header->Write("", TObject::kOverwrite);
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void pCTcalib::FilldEE(TH2D* dEEhist[nStage]) {
  cout<<"Fill dEE"<<endl;
  theConfig->item_int["doGains"] = true;
  //theConfig->item_int["Nbricks"] = theConfig->item_int["Nbricks"];

  theCalibration->ResetHist();
  theEvtRecon->ReadInputFile(theGeometry, theTVcorr, theConfig->item_str["inputFileName"], theCalibration);

  for (int EvtNum = 0; EvtNum < theEvtRecon->nEvents; ++EvtNum) {
    Event thisEvent;
    thisEvent = theEvtRecon->evtList[EvtNum];
    Tf[0] = thisEvent.Thit[0]; Tf[1] = thisEvent.Thit[1];// Front tracker coordinates
    Vf[0] = thisEvent.Vhit[0]; Vf[1] = thisEvent.Vhit[1];
    T[0]  = thisEvent.Thit[2];  T[1] = thisEvent.Thit[3];// Rear tracker coordinates
    V[0]  = thisEvent.Vhit[2];  V[1] = thisEvent.Vhit[3];
    float eStage[nStage];

    for (int stage = 0; stage < nStage; ++stage) {
      
      // Extrapolate the rear track vector to the energy detector stage
      float Tcorr = theGeometry->extrap2D(Ut, T, theGeometry->energyDetectorU(stage));
      float Vcorr = theGeometry->extrap2D(Uv, V, theGeometry->energyDetectorU(stage));
      bool inBounds;

      // Apply the pedestals and TV correction to the stage ADC values
      eStage[stage] = theEvtRecon->GainFac[stage] * theTVcorr->corrFactor(stage, Tcorr, Vcorr, inBounds) * (thisEvent.ADC[stage] - theEvtRecon->Peds[stage]);
    }

    if (eStage[4] > theConfig->item_float["thr4"])       dEEhist[4]->Fill(eStage[3], eStage[4]);
    else if (eStage[3] > theConfig->item_float["thr3"])  dEEhist[3]->Fill(eStage[2], eStage[3]);
    else if (eStage[2] > theConfig->item_float["thr2"] ) dEEhist[2]->Fill(eStage[1], eStage[2]);
    else if (eStage[1] > theConfig->item_float["thr1"])  dEEhist[1]->Fill(eStage[0], eStage[1]);
    else if (eStage[0] > theConfig->item_float["thr0"])  continue;
  }
}  
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void pCTcalib::procWEPLcal(TH2D* REhist[nStage], TH2D* dEEhist[nStage], TH2D* REhist_Tot) {


  // Routine passed to process the WEPL calibration.
  cout << "Entering procWEPLcal for Nbricks=" << theConfig->item_int["Nbricks"] << ".  The list of histograms to fill is" << endl;
  for(int stage = 0; stage < nStage; ++stage) cout << theConfig->item_int["Nbricks"] << " bricks, stage " << stage << ", title=  " << REhist[stage]->GetTitle() << endl;
  theConfig->item_int["doGains"] = true;
  theConfig->item_int["Nbricks"] = theConfig->item_int["Nbricks"];

  theCalibration->ResetHist();
  theEvtRecon->ReadInputFile(theGeometry, theTVcorr, theConfig->item_str["inputFileName"], theCalibration);
  theCalibration->WriteHist(); // For analysis
  
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
  TH1D *hStgEcorr[nStage];
  for (int stage = 0; stage < nStage; ++stage) {
    cout << "procWEPLcal " << theConfig->item_int["Nbricks"] << ": Pedestal measured for stage " << stage << " is " << theEvtRecon->Peds[stage] << " ADC counts " << endl;
    cout << "procWEPLcal " << theConfig->item_int["Nbricks"] << ": Gain correction factor for stage " << stage << " is " << theEvtRecon->GainFac[stage] << endl;
    string Title = "nBricks= " + to_string((long long int)theConfig->item_int["Nbricks"]) + " corrected energy for stage " +
      to_string((long long int)stage);
    if (theConfig->item_str["partType"] == "H") hStgEcorr[stage] = new TH1D(Title.c_str(), Title.c_str(), 400, 15.,  15. +400*0.175);
    else hStgEcorr[stage] = new TH1D(Title.c_str(), Title.c_str(), 400, 15.,  15. +400*0.7);
  }

  //Define more histograms
  string Title = "nBricks " + to_string((long long int)theConfig->item_int["Nbricks"]) + " sum of corrected stage energies";
  TH1D *hTotEcorr;
  if (theConfig->item_str["partType"] == "H") hTotEcorr = new TH1D(Title.c_str(), Title.c_str(), 800, 0., 0.3*800);
  else hTotEcorr = new TH1D(Title.c_str(), Title.c_str(), 840, 0., 1.0*840);
  Title = "nBricks " + to_string((long long int)theConfig->item_int["Nbricks"]) + " calibration phantom thickness";
  TProfile* hLengthProf = new TProfile(Title.c_str(), Title.c_str(), 300, -150, -150 +3.0*100);
  
  float cut0, cut1, cut2, cut3, cut4, cut5;
    
  // Wedge phantom geometry and locations in the T axis,
  float Tw1 = theGeometry->getTWedgeBreaks(1);
  float Tw2 = theGeometry->getTWedgeBreaks(2);
  float Tw3 = theGeometry->getTWedgeBreaks(3);
  float Tw4 = theGeometry->getTWedgeBreaks(4);
  float tBrickEnd = theGeometry->getTWedgeBreaks(4) + 10.0; // Bricks shifted but we take only the first ten cm//???

  // Define the U Position of the bricks now
  float BrickThickness = theGeometry->getBrickThickness(); // Brick thickness; also the wedge thickness
  float Ust[2]; // Polystyrene wedge phantom U coordinates, front and back
  Ust[0] = -2.5 * BrickThickness; // centered around the middle of the bricks, 4 bricks + 1 wedge
  Ust[1] = Ust[0] + BrickThickness; // this is the plane after the wedge
  float Wbricks = BrickThickness * theConfig->item_int["Nbricks"];
  float Uout = Ust[1] + Wbricks; // U coordinate of the last brick, downstream side

  cout << "procWEPLcal: The wedge phantom goes from u=" << Ust[0] << " to u=" << Ust[1] << " for a thickness of "<< Ust[1] - Ust[0]<<" or in WET "
       << 1.03*(Ust[1] - Ust[0])<<endl;
  cout << "procWEPLcal: The " << theConfig->item_int["Nbricks"] << " bricks go from u=" << Ust[1] << " to u=" << Uout << " for a thickness of "<< Uout - Ust[1]<<
    " or in WET "<< 1.03*(Uout - Ust[1])<<endl;
  cout << "The total thickness is "<<Uout - Ust[0]<<" or in WET "<< 1.03*(Uout -Ust[0]) <<endl;
  cout << "procWEPLcal: Stage thresholds are " << theConfig->item_float["thr0"] << " " << theConfig->item_float["thr1"] << " " << theConfig->item_float["thr2"] << " "
       << theConfig->item_float["thr3"] << " " << theConfig->item_float["thr4"] << endl;
  cout << "procWEPLcal: The wedge phantom break points are, in mm: " << Tw1 << ", " << Tw2 << ", " << Tw3 << ", " << Tw4 << endl;
  cout << "procWEPLcal: The calibration brick thickness is " << BrickThickness << " mm " << endl;


  float maxV       = 45.;    // limit of range in V for calibration
  float maxT       = 170.;   // limit of range in T for calibration
  // Start of the event loop

  for (int EvtNum = 0; EvtNum < theEvtRecon->nEvents; ++EvtNum) {
    Event thisEvent;
    thisEvent = theEvtRecon->evtList[EvtNum];

    //Tf[0] = thisEvent.Thit[0] - 1; Tf[1] = thisEvent.Thit[1] - 1;// Front tracker coordinates 1 mm shift
    Tf[0] = thisEvent.Thit[0]; Tf[1] = thisEvent.Thit[1];// Front tracker coordinates    
    Vf[0] = thisEvent.Vhit[0]; Vf[1] = thisEvent.Vhit[1];
    T[0]  = thisEvent.Thit[2];  T[1] = thisEvent.Thit[3];// Rear tracker coordinates
    V[0]  = thisEvent.Vhit[2];  V[1] = thisEvent.Vhit[3];

    float eStage[nStage]; // energy in this stage
    float eTot = 0;
    bool good = true;
    for (int stage = 0; stage < nStage; ++stage) {
      // Extrapolate the rear track vector to the energy detector stage
      float Tcorr = theGeometry->extrap2D(Ut, T, theGeometry->energyDetectorU(stage));
      float Vcorr = theGeometry->extrap2D(Uv, V, theGeometry->energyDetectorU(stage));
      
      // Apply the pedestals and TV correction to the stage ADC values
      bool inBounds;
      eStage[stage] = theEvtRecon->GainFac[stage] * theTVcorr->corrFactor(stage, Tcorr, Vcorr, inBounds) * (thisEvent.ADC[stage] - theEvtRecon->Peds[stage]);
      if (!inBounds){
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

    // calculate geometric WET "wt0" for the wedge phantom using tracker info.

    // First estimate of the u coordinate at entry into the phantom wedge
    float Tin  = theGeometry->extrap2D(Uft, Tf, Ust[0]);
    float Vin  = theGeometry->extrap2D(Ufv, Vf, Ust[0]);

    // First brick past the wedge
    float TinB = theGeometry->extrap2D(Uft, Tf, Ust[1]); 
    float VinB = theGeometry->extrap2D(Ufv, Vf, Ust[1]);

    if (fabs(Vin) > maxV || fabs(Tin) > maxT) continue; // if track is outside of detector in T/V - skip it

    //Back of bricks
    float Tout = theGeometry->extrap2D(Ut, T, Uout);//theGeometry->extrap2D(Uft, Tf, Uout); // theGeometry->extrap2D(Ut, T, Uout); 
    float Vout = theGeometry->extrap2D(Uv, V, Uout);//theGeometry->extrap2D(Ufv, Vf, Uout); // theGeometry->extrap2D(Uv, V, Uout);
    if (fabs(Vout) > maxV || fabs(Tout) > maxT) continue;  // if track is outside of detector in T/V- skip it

    float Uin    = Ust[0];
    bool emptyEvt = false;
    
    // Before the negative slope -- NOT ACCOUNTED
    if(TinB < Tw1) continue;

    // Negative Slope
    else if (TinB >= Tw1 && Tin <= Tw2) {
      if(!getLineIntersection(theEvtRecon->uhitT[0], thisEvent.Thit[0], theEvtRecon->uhitT[1], thisEvent.Thit[1],
			      Uin + BrickThickness, Tw1, Uin, Tw2, Uin, Tin)) continue;
    }
    // Flat part of the wedge - No need to change Vin and Tin
    else if(Tin > Tw2 && Tin < Tw3){
      Uin = Ust[0];
    }
    
    // Positive Slope
    else if (Tin >= Tw3 && TinB <= Tw4) {  // Verify the particle hasn't weirdly scattered
      if (!getLineIntersection(theEvtRecon->uhitT[0], thisEvent.Thit[0], theEvtRecon->uhitT[1], thisEvent.Thit[1],
			       Uin, Tw3 ,Uin + BrickThickness, Tw4, Uin, Tin)) continue;
    }
    // Past the positive slope 
    else if(TinB > Tw4 && TinB < tBrickEnd){
      Uin = Ust[1];
      Tin = TinB;
      Vin = VinB;
    }

    // Past the end of the brick
    else if(TinB > tBrickEnd){
      continue;
    }
    else continue;
    Vin  = theGeometry->extrap2D(Ufv, Vf, Uin);
    float Length = sqrt( pow((Tin - Tout),2) + pow((Vin - Vout),2) + pow((Uin - Uout),2) ); // Length  crossed (converted to WET in Wepl.h)
    if (emptyEvt) {
      hTotEcorr->Fill(eTot);
      for (int stage = 0; stage < nStage; ++stage){
	hStgEcorr[stage]->Fill(eStage[stage]);
      }
    }

    float RSP = 1.030; // Transformation to WET
    hLengthProf->Fill(TinB,Length*RSP);
    if(emptyEvt) continue;
    //Energy cuts
    if (theConfig->item_str["partType"] == "He") {  
      cut0 = 87. ; cut1 = 95.;  cut2 = 100.;
      cut3 = 120.; cut4 = 300.; cut5 = 270.;
    } else {
      cut0 = 19.; cut1 = 21.; cut2 = 25.;
      cut3 = 35.; cut4 = 75.; cut5 = 77.;
    }
    float E_tot = eStage[0]+eStage[1]+eStage[2]+eStage[3]+ eStage[4];
    ProfileE_Tot->Fill(T[0], E_tot);
    ProfileE_Tot_f->Fill(Tf[0], E_tot);
    REhist_Tot->Fill(Length*RSP,E_tot);

    if (eStage[4] > theConfig->item_float["thr4"]){
      float E_tot = eStage[0]+eStage[1]+eStage[2]+eStage[3]+ eStage[4];
      REhist[4]->Fill(Length*RSP, eStage[4]);
      //REhist_Tot->Fill(Length*RSP,E_tot);
    }
    else if (eStage[3] > theConfig->item_float["thr3"]){
      float E_tot = eStage[0] +eStage[1] +eStage[2] +eStage[3];
      REhist[3]->Fill(Length*RSP, eStage[3]);
      //REhist_Tot->Fill(Length*RSP,E_tot);
    }
    else if (eStage[2] > theConfig->item_float["thr2"]){
      float E_tot = eStage[0]+eStage[1]+eStage[2];
      REhist[2]->Fill(Length*RSP, eStage[2]);
      //REhist_Tot->Fill(Length*RSP,E_tot);
    }
    else if (eStage[1] > theConfig->item_float["thr1"]){
      float E_tot = eStage[0]+eStage[1];
      REhist[1]->Fill(Length*RSP, eStage[1]);
      //REhist_Tot->Fill(Length*RSP,E_tot);
    }
    else if (eStage[0] > theConfig->item_float["thr0"]){
      float E_tot = eStage[0];
      REhist[0]->Fill(Length*RSP, eStage[0]); // Particles end in stage 0
      //REhist_Tot->Fill(Length*RSP,E_tot);
    }    
  } // End of the event loop
  pCTcalibRootFile->mkdir("runCorrEnrgs");
  pCTcalibRootFile->cd("runCorrEnrgs");
  for (int stage = 0; stage < nStage; ++stage) hStgEcorr[stage]->Write("",TObject::kOverwrite);
  hTotEcorr->Write("",TObject::kOverwrite);

  pCTcalibRootFile->cd();
  pCTcalibRootFile->mkdir("Phantom");  
  pCTcalibRootFile->cd("Phantom");
  hLengthProf->Write("",TObject::kOverwrite);
  pCTcalibRootFile->cd();

  // memory management
  delete hTotEcorr;
  for (int stage = 0; stage < nStage; ++stage) delete (hStgEcorr[stage]);
}
///

///////////////////////////////////////////////////////////////////
// 2D equation to find line intersection
//////////////////////////////////////////////////////////////////////
bool pCTcalib::getLineIntersection(float p0u, float p0t, float p1u, float p1t,  // Two points on the first line
				   float p2u, float p2t, float p3u, float p3t,  // Two points on the second line
				   float &ipu, float &ipt) { // Return the intersection point
  if (p0u == p1u || p2u == p3u) return false;
    
  float slope1 = (p1t - p0t) / (p1u - p0u); // Slope of the first line
  float slope2 = (p3t - p2t) / (p3u - p2u); // Slope of the second line
  if (slope1 == slope2) return false;                   // There is no intersection
  float int1 = p0t - slope1 * p0u; // Intercept of the first line
  float int2 = p2t - slope2 * p2u; // Intercept of the second line
  
  ipu = (int1 - int2) / (slope2 - slope1); // intersection point U
  ipt = (slope2 * int1 - slope1 * int2) / (slope2 - slope1); // Intersection Point T
  return true;
}
