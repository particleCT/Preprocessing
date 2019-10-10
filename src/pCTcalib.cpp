// Code to process raw calibration data to derive the TV and WEPL calibration
// constants.
// R.P. Johnson  September 18, 2016, adapted from Vladimir Bashkirov's
// calibration code.
// Modified by Lennart Volz December 8, 2017:
//        Applied modifications: New energy cuts for helium fragmentation in
// bool EnergyCheck()
//                               Different values for the kink and high E
// inter-/extrapolation for He
//                               for smoother calibration curves

#include "pCTcalib.h"

bool getLineIntersection(double p0u, double p0t,                             // Two points on the first line
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

  /*  cout << "Point 1: u=" << p0u << " t=" << p0t << endl;
      cout << "Point 2: u=" << p1u << " t=" << p1t << endl;
      cout << "Point 3: u=" << p2u << " t=" << p2t << endl;
      cout << "Point 4: u=" << p3u << " t=" << p3t << endl;
      cout << "Line 1: slope=" << slope1 << " intercept=" << int1 << endl;
      cout << "Line 2: slope=" << slope2 << " intercept=" << int2 << endl;
      cout << "Intersection u=" << ipu << " t= " << ipt << endl;
      cout << "Line 1: t predicted= " << ipu*slope1 + int1 << endl;
      cout << "Line 2: t predicted= " << ipu*slope2 + int2 << endl;  */

  return true;
}

struct calStuff { // Structure used to facilitate passing information to the 5
                  // parallel calibration threads when they start
  string inputFileName;
  string Outputdir;
  int Nbricks;
  int max_events;
  int max_time;
  bool useTemp;
  string partType;
  float Thr[5];
  int pdstlr[5];
  string OsName;
  bool reCalibrate;
  double topW1;  // minimum t value for range in top of wedge, relative to the
                 // center, both sides
  double topW2;  // maximum t value for range in top of wedge, relative to the
                 // center, both sides
  double brickW; // half width of range in t for brick-only protons
  double emptyW; // half width of empty region in t to use for calibration
};

bool EnrgCut(calStuff stuff, float Estage[5], float Etot, float cut0, float cut1, float cut3) {

  bool dropEvent = false;

  //
  // Cut on energy deposit not compatible with Bragg-Peak (cut0, cut1) adapted
  // from Vladimir Bashkirov
  // For helium: Use the stages as dE-E detector and check if the particle's
  // energy
  // deposit is compatible with the parameterized stage response
  // (Last modified: Lennart Volz, December 2017)
  //

  if (Estage[4] > stuff.Thr[4]) {
    if (stuff.partType == "He") {
      if (Estage[3] < (0.00045 * Estage[4] * Estage[4] - 0.42 * Estage[4] + 237) ||
          Estage[3] > (0.00045 * Estage[4] * Estage[4] - 0.42 * Estage[4] + 260))
        dropEvent = true;
    }

    // DELTA E-E FILTER DOES NOT SIGNIFICANTLY CHANGE CALIBRATION
    /*	if (stuff.partType == "H") {
                if (Estage[3]<(0.0034*Estage[4] * Estage[4] - 0.672145*Estage[4]
       + 69.) || Estage[3]>(0.00351528*Estage[4] * Estage[4] -
       0.703048*Estage[4] + 79.5048)) dropEvent = true;
            }*/
    if (Estage[3] < cut0 || Estage[2] < cut0 || Estage[1] < cut0 || Estage[0] < cut0 || Estage[4] > cut3)
      dropEvent = true;
  } else if (Estage[3] > stuff.Thr[3]) {
    if (stuff.partType == "He") {
      if (Estage[2] < (0.00085 * Estage[3] * Estage[3] - 0.5664 * Estage[3] + 226) ||
          Estage[2] > (0.00085 * Estage[3] * Estage[3] - 0.5664 * Estage[3] + 254))
        dropEvent = true;
    }
    // DELTA	E-E FILTER DOES	NOT SIGNIFICANTLY CHANGE CALIBRATION
    /*	if (stuff.partType == "H") {
                if (Estage[2]<(0.0038*Estage[3] * Estage[3] - 0.71*Estage[3] +
       67.92) || Estage[2]>(0.0039*Estage[3] * Estage[3] - 0.74*Estage[4] +
       78.37)) dropEvent = true;
            }*/

    if (Estage[2] < cut0 || Estage[1] < cut0 || Estage[0] < cut0 || Estage[3] > cut3)
      dropEvent = true;
  } else if (Estage[2] > stuff.Thr[2]) {
    if (stuff.partType == "He") {
      if (Estage[1] < (0.00085 * Estage[2] * Estage[2] - 0.5664 * Estage[2] + 220) ||
          Estage[1] > (0.00085 * Estage[2] * Estage[2] - 0.5664 * Estage[2] + 248))
        dropEvent = true;
    }
    // DELTA	E-E FILTER DOES	NOT SIGNIFICANTLY CHANGE CALIBRATION
    /*	if (stuff.partType == "H") {
                if (Estage[1]<(0.0040*Estage[2] * Estage[2] - 0.71*Estage[2] +
       66.17) || Estage[1]>(0.0036*Estage[2] * Estage[2] - 0.73*Estage[2] +
       76.57)) dropEvent = true;
            }*/

    if (Estage[1] < cut0 || Estage[0] < cut0 || Estage[2] > cut3)
      dropEvent = true;
  } else if (Estage[1] > stuff.Thr[1]) {
    if (stuff.partType == "He") {
      if (Estage[0] < (0.00085 * Estage[1] * Estage[1] - 0.5664 * Estage[1] + 220) ||
          Estage[0] > (0.00085 * Estage[1] * Estage[1] - 0.5664 * Estage[1] + 246))
        dropEvent = true;
    }
    // DELTA	E-E FILTER DOES	NOT SIGNIFICANTLY CHANGE CALIBRATION
    /*	if (stuff.partType == "H") {
                if (Estage[0]<(0.0040*Estage[1] * Estage[1] - 0.72*Estage[1] +
       66.93) || Estage[0]>(0.0040*Estage[1] * Estage[1] - 0.75*Estage[1] +
       76.38)) dropEvent = true;
            }*/

    if (Estage[0] < cut1 || Estage[1] > cut3)
      dropEvent = true;
  } else if (Estage[0] > cut3)
    dropEvent = true;

  if (stuff.partType == "He") {
    if (Etot > 801.52)
      dropEvent = true; // Intitial Energy of helium at HIT = 200.36 MeV/u
  } else if (Etot > 200.)
    dropEvent = true; // Initial Energy of protons at HIT = 200.11 MeV

  return dropEvent;

  /***DEFAULT PARAMETERS FOR PROTON dE-E
      # dEE1 = 0.00395087,-0.715705,65.936,0.00400261,-0.746907,76.382,
  # dEE2 = 0.00404936,-0.709696,66.167,0.00363538,-0.72553,76.5712,
  # dEE3 = 0.00380467,-0.71088,67.9228,0.00385641,-0.742083,78.3689,
  # dEE4 = 0.00346354,-0.671845,69.0588,0.00351528,-0.703048,79.5048,
  */
}

void procWEPLcal(calStuff stuff, TVcorrection TVcal, Histogram2D *REhist[nStage], Histogram2D *dEEhist[nStage],
                 pCTgeo Geometry) {
  // Routine passed to multiple threads to process the WEPL calibration.
  // In principle it is a private method of the pCTcalib class, but the
  // compilers will not allow it
  // to be multithreaded if it is a class member.

  cout << "Entering procWEPLcal for Nbricks=" << stuff.Nbricks << ".  The list of histograms to fill is" << endl;
  for (int stage = 0; stage < nStage; ++stage) {
    cout << stuff.Nbricks << " bricks, stage " << stage << ", title=  " << REhist[stage]->Title() << endl;
  }

  TVcorrection *TVcorr = &TVcal;
  EvtRecon *procEvt =
      new EvtRecon(Geometry, TVcorr, stuff.inputFileName, stuff.Outputdir, stuff.max_events, stuff.max_time, 0, 0,
                   stuff.Nbricks, stuff.useTemp, true, stuff.partType, stuff.pdstlr, stuff.reCalibrate,
                   stuff.OsName); // Process the raw data file to produce the event list

  double Ut[2], Uv[2], T[2], V[2];
  double Uft[2], Ufv[2], Tf[2], Vf[2];

  Uft[0] = procEvt->uhitT[0];
  Uft[1] = procEvt->uhitT[1];
  Ufv[0] = procEvt->uhitV[0];
  Ufv[1] = procEvt->uhitV[1];
  Ut[0] = procEvt->uhitT[2];
  Ut[1] = procEvt->uhitT[3];
  Uv[0] = procEvt->uhitV[2];
  Uv[1] = procEvt->uhitV[3];
  cout << "procWEPLcal: Tracker plane u locations:" << endl;
  cout << "Front tracker t: " << Uft[0] << " " << Uft[1] << endl;
  cout << "Front tracker v: " << Ufv[0] << " " << Ufv[1] << endl;
  cout << "Front tracker t: " << Ut[0] << " " << Ut[1] << endl;
  cout << "Front tracker t: " << Uv[0] << " " << Uv[1] << endl;
  Histogram *hStgEcorr[nStage];
  for (int stage = 0; stage < nStage; ++stage) {
    cout << "procWEPLcal " << stuff.Nbricks << ": Pedestal measured for stage " << stage << " is "
         << procEvt->Peds[stage] << " ADC counts " << endl;
    cout << "procWEPLcal " << stuff.Nbricks << ": gain correction factor for stage " << stage << " is "
         << procEvt->CorFacs[stage] << endl;
    string Title = "nBricks= " + to_string((long long int)stuff.Nbricks) + " corrected energy for stage " +
                   to_string((long long int)stage);
    if (stuff.partType == "H")
      hStgEcorr[stage] = new Histogram(400, 15., 0.175, Title, "Energy (MeV)", "protons");
    else
      hStgEcorr[stage] = new Histogram(400, 15., 0.7, Title, "Energy (MeV)", "He ions");
  }
  string Title = "nBricks " + to_string((long long int)stuff.Nbricks) + " sum of corrected stage energies";
  Histogram *hTotEcorr;
  if (stuff.partType == "H")
    hTotEcorr = new Histogram(800, 0., 0.3, Title, "Energy (MeV)", "protons");
  else
    hTotEcorr = new Histogram(840, 0., 1.0, Title, "Energy (MeV)", "He ions");
  Title = "nBricks " + to_string((long long int)stuff.Nbricks) + " calibration phantom thickness";
  ProfilePlot *hwt0Prof = new ProfilePlot(100, -150, 3.0, Title, "T (mm)", "Thickness (mm)");

  // Wedge phantom geometry and locations,
  double Tw1 = Geometry.getTWedgeBreaks(1);
  double Tw2 = Geometry.getTWedgeBreaks(2);
  double Tw3 = Geometry.getTWedgeBreaks(3);
  double Tw4 = Geometry.getTWedgeBreaks(4);
  double tBrickEnd = Geometry.getTWedgeBreaks(4);
  double BrickThickness = Geometry.getBrickThickness(); // Brick thickness; also the wedge thickness

  double Ust[2]; // Polystyrene wedge phantom U coordinates, front and back
  Ust[0] = -2.5 * BrickThickness;
  Ust[1] = Ust[0] + BrickThickness;
  double Wbricks = BrickThickness * stuff.Nbricks;
  double Uout = Ust[1] + Wbricks; // U coordinate of the last brick, downstream side
  cout << "procWEPLcal: The wedge phantom goes from u=" << Ust[0] << " to u=" << Ust[1] << endl;
  cout << "procWEPLcal: The " << stuff.Nbricks << " bricks go from u=" << Ust[1] << " to u=" << Uout << endl;
  cout << "procWEPLcal: Stage thresholds are " << stuff.Thr[0] << " " << stuff.Thr[1] << " " << stuff.Thr[2] << " "
       << stuff.Thr[3] << " " << stuff.Thr[4] << endl;
  cout << "procWEPLcal: The wedge phantom break points are, in mm: " << Tw1 << ", " << Tw2 << ", " << Tw3 << ", " << Tw4
       << endl;
  cout << "procWEPLcal: The calibration brick thickness is " << BrickThickness << " mm " << endl;

  float cut0, cut1, cut3; // energy cuts to reject events incompatible with
                          // expected Bragg peak.

  if (stuff.partType == "He") {
    cut0 = 60.;
    cut1 = 80.;
    cut3 = 300;
  } else {
    cut0 = 15.;
    cut1 = 20.;
    cut3 = 100.;
  }

  double maxV = 45.;     // limit of range in V for calibration
  double maxT = 170.;    // limit of range in T for calibration
  double deltaInt = 20.; // intercept and slope to calculate maximum difference
                         // in V or T vs stage number
  double deltaSlope = 30.;
  double tEmpty = 30.;      // t location of center of range for selecting empty events
  double eMin4Empty = 40.0; // minimum energy in stage 4 for an empty event
  double brickC = 6.5;      // Distance in t from wedge for center of range to take
                            // protons through bricks only

  float mxEstage;
  if (stuff.partType == "He") {
    mxEstage = 400.;
  } else {
    mxEstage = 100.;
  }
  if (procEvt->useTmpFile)
    procEvt->reopenTmpFile();
  double diff = 0;
  int N = 0;
  for (int EvtNum = 0; EvtNum < procEvt->nEvents; ++EvtNum) {
    Event thisEvent;
    if (procEvt->useTmpFile)
      procEvt->readTmp(thisEvent);
    else
      thisEvent = procEvt->evtList[EvtNum];
    Tf[0] = thisEvent.Thit[0]; // Front tracker coordinates
    Tf[1] = thisEvent.Thit[1];
    Vf[0] = thisEvent.Vhit[0];
    Vf[1] = thisEvent.Vhit[1];
    T[0] = thisEvent.Thit[2]; // Rear tracker coordinates
    T[1] = thisEvent.Thit[3];
    V[0] = thisEvent.Vhit[2];
    V[1] = thisEvent.Vhit[3];

    float eStage[nStage];
    float eTot = 0;
    bool good = true;
    for (int stage = 0; stage < nStage; ++stage) {
      // Extrapolate the rear track vector to the energy detector stage
      double Tcorr = Geometry.extrap2D(Ut, T, Geometry.energyDetectorU(stage));
      double Vcorr = Geometry.extrap2D(Uv, V, Geometry.energyDetectorU(stage));

      // Apply the pedestals and TV correction to the stage ADC values
      bool inBounds;
      eStage[stage] = procEvt->CorFacs[stage] * TVcal.corrFactor(stage, Tcorr, Vcorr, inBounds) *
                      (thisEvent.ADC[stage] - procEvt->Peds[stage]);
      if (!inBounds || eStage[stage] > mxEstage) {
        good = false;
        break;
      }
      if (abs(Vcorr) > 40.0 || abs(Tcorr) > 160.0) {
        good = false;
        break;
      } // To avoid protons that scatter out of the device
      eTot += eStage[stage];
    }
    if (!good)
      continue;
    //
    // Skip events not compatible with parameterized energy deposit
    // (Lennart Volz, December 2017)
    //
    // bool dropEvt = EnrgCut(stuff, eStage, eTot, cut0, cut1,cut3); //ENERGY
    // CUTS FOR REMOVING EVENT PILE-UP AND BAD EVENTS
    // if (dropEvt) continue;

    // calculate geometric WET "wt0" for the wedge phantom using tracker info.
    double Uin = Ust[0];                          // First estimate of the u coordinate at entry into the phantom
    double Tin = Geometry.extrap2D(Uft, Tf, Uin); // extrapolate the front track to the phantom entrance
    double Vin = Geometry.extrap2D(Ufv, Vf, Uin);
    double TinB = Geometry.extrap2D(Uft, Tf, Ust[1]); // extrapolate to the front of the first brick
    double VinB = Geometry.extrap2D(Ufv, Vf, Ust[1]);
    if (fabs(Vin) > maxV || fabs(Tin) > maxT)
      continue; // if track is outside of FOW - skip it

    double Tout = Geometry.extrap2D(Ut, T, Uout); // extrapolate rear track to back of phantom
    double Vout = Geometry.extrap2D(Uv, V, Uout);

    if (fabs(Vout) > maxV ||
        fabs(Tout) > maxT /*|| fabs(TinB-Tout) > deltaInt + deltaSlope*sqrt(float(stuff.Nbricks))*/)
      continue; // skip scattered protons //NO SCATTERING FILTER NEEDED

    bool emptyEvt = false;
    int cs;

    if ((Tin > Tw2 && Tin < Tw3) || (TinB < Tw1) ||
        (TinB > Tw4)) { // ONLY WEDGE SLOPE TO BE USED EVERYTHING ELSE CREATES ISSUES!
      if (TinB > tBrickEnd)
        emptyEvt = true;
      cs = 1;
      continue;
    }

    /*************NOT USED IN CURRETN VERSION ANYMORE
            if (TinB > tBrickEnd) {                                       // the
    bricks were chopped off at or near the +t end of the wedge
                if (fabs(TinB - (Tw4 + tEmpty))<stuff.emptyW &&
    eStage[4]>eMin4Empty) cs=6;  // we'll use the straight-through protons in
    this range
                else cs=1;                                                //
    Offset increased 9/12/2018 so we can move bricks away from wedge end
                emptyEvt = true;                                          //
    missed all of the phantom plastic
            } else if (Tin >= Tw2 && Tin <= Tw3) {                        //
    proton goes through the wedge top
                float dC = abs(Tin - (Tw2+Tw3)/2.);
                if (dC>stuff.topW1 && dC<stuff.topW2) cs=8;               //
    Keep some to get data with full bricks, but avoid the center gap
                else cs=2;
            } else if (TinB <= Tw1) {
                Tin=TinB;
                Vin=VinB;
                Uin=Ust[1];
                if (stuff.Nbricks > 0 && abs(TinB - (Tw1 - brickC)) <
    stuff.brickW) cs=7;  // no wedges hit, just bricks; we'll use them in this
    range
                else cs=3;
            }
    FLAT BRICKS NOT USED ANYMORE***************/
    else if (TinB >= Tw1 && Tin <= Tw2) { // in the -t wedge of the phantom; Uin
                                          // and Tin get overwritten here
      if (getLineIntersection(procEvt->uhitT[0], thisEvent.Thit[0], procEvt->uhitT[1], thisEvent.Thit[1],
                              Uin + BrickThickness, Tw1, Uin, Tw2, Uin, Tin)) {
        Vin = Geometry.extrap2D(Ufv, Vf, Uin);
        cs = 4;
      } else
        continue;
    } else if (Tin >= Tw3 && TinB <= Tw4) { //&& Tout<tBrickEnd) {    // in the +t wedge of
      // the phantom; Uin and Tin get overwritten here
      // ALSO EXCLUDE PARTICLES THAT LEFT THE BRICKS
      if (getLineIntersection(procEvt->uhitT[0], thisEvent.Thit[0], procEvt->uhitT[1], thisEvent.Thit[1], Uin, Tw3,
                              Uin + BrickThickness, Tw4, Uin, Tin)) {
        Vin = Geometry.extrap2D(Ufv, Vf, Uin);
        // if (stuff.Nbricks == 0 || (Tw4-TinB)>0.0) {  // stay away from the
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
    double wt0;
    double wt0old;
    if (emptyEvt) {
      wt0 = 0.0;
      hTotEcorr->entry(eTot);
      for (int stage = 0; stage < nStage; ++stage)
        hStgEcorr[stage]->entry(eStage[stage]);
    } else {
      Vin = Geometry.extrap2D(Ufv, Vf, Uin);
      //            if (fabs(Vin-Vout) > deltaInt +
      // deltaSlope*sqrt(float(stuff.Nbricks))) continue;   // Too much scatter;
      // skip the event //NO SCATTERING FILTER NEEDED

      wt0 = sqrt((Tin - Tout) * (Tin - Tout) + (Vin - Vout) * (Vin - Vout) +
                 (Uin - Uout) * (Uin - Uout)); // Straight line path length to
                                               // get the physical thickness
                                               // crossed (is converted to WET
                                               // in Wepl.h)

      /*
      //POTENTIAL TO USE CSP CODE FOR MORE ACCURACY
                  wt0old = wt0;
                  //More accurate path estimation with CSP algorithm
                  float alpha1 = 1.01+0.43*pow((wt0*1.03)/261.0,2);//wt0*1.03:
      rough estimate of water equivalent thickness crossed; 261.0 range of 200
      MeV protons in Water, 1.03 RSP of calibration phantom
                  float alpha2 = 0.99-0.46*pow((wt0*1.03)/261.0,2);

                  float lambda1 = alpha1*wt0; // in Collins-Fekete et al. P_opt
      = alpha*|X2-X0|*P/||P||, despite the naming scheme wt0 is the physical
      thickness
                  float lambda2 = alpha2*wt0;
                  float t_mlp, v_mlp;

                  float distV = 217.3 -167.2;    // U positions of the 8 tracker
      boards in the Phase-II scanner
                  float distT = 211.4 - 161.4;

                  float dir_t_0 = Tf[1]-Tf[0];
                  float dir_v_0 = Vf[1]-Vf[0];
                  dir_t_0 = dir_t_0/sqrt(dir_t_0*dir_t_0 + dir_v_0*dir_v_0 +
      distT*distT); //Normalizing the direction vector
                  dir_v_0 = dir_v_0/sqrt(dir_t_0*dir_t_0 + dir_v_0*dir_v_0 +
      distT*distT);

                  float dir_t_2 = T[1]-T[0];
                  float dir_v_2 = V[1]-V[0]);

                  dir_t_2 = dir_t_2/sqrt(dir_t_2*dir_t_2 + dir_v_2*dir_v_2 +
      distT*distT); //Normalizing the direction vector
                  dir_v_2 = dir_v_2/sqrt(dir_t_2*dir_t_2 + dir_v_2*dir_v_2 +
      distT*distT);

                  float P0_t = dir_t_0*lambda1;
                  float P2_t = dir_t_2*lambda2;

                  float A_t = (2*Tin + P0_t - 2*Tout + P2_t);
                  float B_t = (3*Tout - 3*Tin - 2*P0_t - P2_t);

                  float P0_v = dir_v_0*lambda1;
                  float P2_v = dir_v_2*lambda2;

                  float A_v = (2*Vin + P0_v - 2*Vout + P2_v);
                  float B_v = (3*Vout - 3*Vin - 2*P0_v - P2_v);

                  float t_mlp_old=Tin;
                  float v_mlp_old=Vin;
                  float u_old=Uin;

                  float step = 0.1; //0.1 mm steps, could be reduced for speed
      (up to 250mm path length for calibration run with 4 bricks)
                  float t;
                  float u = Uin+step; //first step

                  // Step through until reach the exit depth
                  wt0 = 0;

                  while(u<=Uout){
                              // Inside the object so calculate MLP in water
                              //S = At^3+Bt^2 + Ct + D

                              t = (u-Uin)/(Uout-Uin);

                              t_mlp = t*(t*(t*A_t + B_t) + P0_t) + Tin;

                              v_mlp = t*(t*(t*A_v + B_v) + P0_v) + Vin;

                              wt0+=sqrt((t_mlp-t_mlp_old)*(t_mlp-t_mlp_old)+(v_mlp-v_mlp_old)*(v_mlp-v_mlp_old)+(u-u_old)*(u-u_old));

                              t_mlp_old=t_mlp; v_mlp_old=v_mlp; u_old=u;
                              u=u+step;
                  }
                  u = u-step;
                  sqrt((Tout-t_mlp_old)*(Tout-T_mlp_old)+(Vout-v_mlp_old)*(Vout-v_mlp_old)+(Uout-u)*(Uout-u));
      //don't forget to take also teh last step that was not computed due to the
      while loop break

      */
      // END CSP OPTIMIZED ESTIMATE
    }

    hwt0Prof->entry(TinB, wt0);
    //
    // For each stage where the proton stops, increment the corresponding
    // range-energy table cell content.
    // Energy cuts are now performed in seperate funcion (Last modified: Lenny,
    // December 2017)

    if (eStage[4] > stuff.Thr[4]) {
      if (cs > 3) {
        if (stuff.Nbricks >= 2)
          continue;
        // if(wt0>50) cout << TinB << " " << Tin << " " << Tout << endl;
        REhist[4]->entry(wt0, eStage[4]);
        dEEhist[4]->entry(eStage[3], eStage[4]);
      }
    } else if (eStage[3] > stuff.Thr[3]) {
      if (cs > 3) {
        REhist[3]->entry(wt0, eStage[3]);
        dEEhist[3]->entry(eStage[2], eStage[3]);
      }
    } else if (eStage[2] > stuff.Thr[2]) {
      if (cs > 3) {
        REhist[2]->entry(wt0, eStage[2]);
        dEEhist[2]->entry(eStage[1], eStage[2]);
      }
    } else if (eStage[1] > stuff.Thr[1]) {
      if (cs > 3) {
        REhist[1]->entry(wt0, eStage[1]);
        dEEhist[1]->entry(eStage[0], eStage[1]);
      }
    } else if (eStage[0] > stuff.Thr[0]) {
      if (cs > 3)
        REhist[0]->entry(wt0, eStage[0]);
    }

  } // End of the event loop

  string histFileName = stuff.Outputdir + "/runCorrEnrgs_" + to_string((long long int)stuff.Nbricks) + ".gp";
  FILE *oFile = fopen(histFileName.c_str(), "w");
  string setTerm;
  if (stuff.OsName == "Windows")
    setTerm = "set terminal wxt size 1300, 600\n";
  else
    setTerm = "set terminal x11 size 1300, 600\n";
  if (oFile != NULL) {
    fprintf(oFile, setTerm.c_str());
    fprintf(oFile, "set multiplot title 'procWEPLcal %d: Energies for file %s' "
                   "layout 3,2 columnsfirst scale 1.0,1.0\n",
            stuff.Nbricks, stuff.inputFileName.c_str());
    float xLow, xHigh;
    for (int stage = 0; stage < nStage; ++stage) {
      int ret = hStgEcorr[stage]->FWHMboundaries(xLow, xHigh);
      if (ret == 0) {
        float peak = hStgEcorr[stage]->mean(xLow, xHigh);
        fprintf(oFile, "#Energy peak for stage %d is at E=%9.3f\n", stage, peak);
      }
    }
    int ret = hTotEcorr->FWHMboundaries(xLow, xHigh);
    if (ret == 0) {
      float peak = hTotEcorr->mean(xLow, xHigh);
      fprintf(oFile, "#Energy peak for the sum is at E=%9.3f\n", peak);
    }
    for (int stage = 0; stage < nStage; ++stage) {
      hStgEcorr[stage]->plot(oFile);
    }
    hTotEcorr->plot(oFile);
    fprintf(oFile, "unset multiplot\n");
    fprintf(oFile, "show label\n");
    fprintf(oFile, "unset label\n");
    fclose(oFile);
  } else {
    cout << "Unable to open the histogram file " << histFileName << endl;
  }

  /*    for (int stage=0; stage<nStage; stage++) {
          histFileName = Outputdir + "/RngEnrg2D_" + to_string((long long
     int)Nbricks) + "_" + to_string((long long int)stage) +".gp";
          FILE* oFile = fopen(histFileName.c_str(),"w");
          if (oFile != NULL) {
              REhist[stage]->plot(oFile);
          } else {cout << "Unable to open the histogram file " << histFileName
     << endl;}
      }
  */
  histFileName = stuff.Outputdir + "/Phantom_" + to_string((long long int)stuff.Nbricks) + ".gp";
  oFile = fopen(histFileName.c_str(), "w");
  if (oFile != NULL) {
    hwt0Prof->plot(oFile);
  } else {
    cout << "Unable to open the profile plot file " << histFileName << endl;
  }
  hwt0Prof->print(stuff.Outputdir + "/Phantom_" + to_string((long long int)stuff.Nbricks) + ".txt");

  delete (procEvt);
  delete (hTotEcorr);
  for (int stage = 0; stage < nStage; ++stage)
    delete (hStgEcorr[stage]);
};

pCTcalib::pCTcalib(string inputFileName, string Outputdir, int max_events, int max_time, int n_debug, int n_plot,
                   bool useTemp, string version, string WcalibFile, string TVcorrFile, string minDate, string maxDate,
                   int minRun, int maxRun, double tWedgeOffset, string partType, int pdstlr[5], bool realTimeCal,
                   std::string OsName, bool redo, bool normalise) {

  cout << "\nEntering pCTcalib, the constructor for the pCT calibration "
          "program, version " << version << endl;
  if (redo)
    cout << "This run will only analyze the existing calibration histograms." << endl;
  cout << "Filename for the list of calibration data files = " << inputFileName << endl;
  cout << "The output directory is " << Outputdir << endl;
  cout << "The maximum number of events is " << max_events << endl;
  cout << "The maximum time is " << max_time << endl;
  cout << "The number of debug event printouts is " << n_debug << endl;
  cout << "The number of events to plot is " << n_plot << endl;
  if (useTemp)
    cout << "Temporary data will be stored in disk files " << endl;
  cout << "The input WEPL calibration file is " << WcalibFile << endl;
  cout << "The input TV correction calibration file is " << TVcorrFile << endl;
  cout << "Input date range = " << minDate << " to " << maxDate << endl;
  cout << "The input run range is " << minRun << " to " << maxRun << endl;
  cout << "The calibration wedge offset is assumed to be " << tWedgeOffset << endl;
  cout << "The particle type is assumed to be " << partType << endl;
  cout << "Real-time pedestal and gain recalibration is done? " << realTimeCal << endl;
  cout << "The operating system is " << OsName << endl;
  for (int stage = 0; stage < nStage; ++stage) {
    cout << "  For stage " << stage << " the range in which to look for pedestals begins at " << pdstlr[stage]
         << " ADC counts " << endl;
  }
  CalFile = inputFileName;
  this->normalise = normalise;
  this->Outputdir = Outputdir;
  this->max_events = max_events;
  this->max_time = max_time;
  this->n_debug = n_debug;
  this->n_plot = n_plot;
  useTempFile = useTemp;
  programVersion = version;
  this->WcalibFile = WcalibFile;
  this->TVcorrFile = TVcorrFile;
  this->minDate = minDate;
  this->maxDate = maxDate;
  this->minRun = minRun;
  this->maxRun = maxRun;
  this->partType = partType;
  this->OsName = OsName;
  this->realTimeCal = realTimeCal;
  this->redo = redo;
  if (partType == "H")
    EnergyBinWidth = 0.25;
  else
    EnergyBinWidth = 1.0;
  RangeBinWidth = 1.0;

  // Operating system dependent setting of the Gnuplot terminal type for
  // histogram plotting
  if (OsName == "Windows")
    setTerm = "set terminal wxt size 1300, 600\n";
  else
    setTerm = "set terminal x11 size 1300, 600\n";

  // Set default parameters for this code

  // LOW E INTERPOLATION: NEVER USED DUE OF THRESHOLD OF 1MeV
  if (partType == "He") {
    k1[0] = 0;
    k1[1] = 2;
    k1[2] = 2;
    k1[3] = 2;
    k1[4] = 3;
  } else {
    k1[0] = 0;
    k1[1] = 5;
    k1[2] = 5;
    k1[3] = 5;
    k1[4] = 6;
  }

  // KINK REGION INTERPOLATION: DISCONTINUITY IN STATISTICS WHEN ANOTHER BRICK
  // IS ADDED (PROTONS CAN SCATTER OUTSIDE OF BRICKS IN OLD STANDARD)
  if (partType == "He") { // Modified by Lennart, December 2017
    j1[0] = 180;
    j1[1] = 180;
    j1[2] = 180;
    j1[3] = 183;
    j2[0] = 192;
    j2[1] = 190;
    j2[2] = 192;
    j2[3] = 193;
    j3[0] = 208;
    j3[1] = 208;
    j3[2] = 210;
    j3[3] = 213;
    j4[0] = 216;
    j4[1] = 218;
    j4[2] = 218;
    j4[3] = 222;
  } else { //
    j1[0] = 222;
    j1[1] = 222;
    j1[2] = 222;
    j1[3] = 222;
    j2[0] = 230;
    j2[1] = 230;
    j2[2] = 230;
    j2[3] = 230;
    j3[0] = 264;
    j3[1] = 264;
    j3[2] = 268;
    j3[3] = 272;
    j4[0] = 272;
    j4[1] = 272;
    j4[2] = 274;
    j4[3] = 278;
  }

  // EXTRAPOLATION AT HIGH ENERGIES (LOW STATISTICS THERE - EXTRAPOLATION
  // NEEDED)
  if (partType == "He") {
    i1[4] = 250;
    i1[3] = 247;
    i1[2] = 247;
    i1[1] = 230;
    i1[0] = 230;
    i3 = 300;
    i2[4] = 265;
    i2[3] = 255;
    i2[2] = 255;
    i2[1] = 240;
    i2[0] = 242;
  } else {
    // NO EXTRAPOLATION
    i1[4] = 200;
    i1[3] = 288;
    i1[2] = 288;
    i1[1] = 280;
    i1[0] = 280;
    i3 = 330;
    i2[4] = 240;
    i2[3] = 304;
    i2[2] = 304;
    i2[1] = 292;
    i2[0] = 292;
  }

  if (partType == "H") {
    EG4stage[0] = 25.25; // MC derived stage energies, used to calibrate to MeV
                         // for protons
    EG4stage[1] = 28.01;
    EG4stage[2] = 32.76;
    EG4stage[3] = 42.62;
    EG4stage[4] = 67.71;
  } else {
    EG4stage[0] = 100.; // MC derived stage energies, used to calibrate to MeV for He
    EG4stage[1] = 111.;
    EG4stage[2] = 129.;
    EG4stage[3] = 166.;
    EG4stage[4] = 279.;
  }
  sG4stage[0] = .04; // From MC, expected peak half widths (not used. . .)
  sG4stage[1] = .04;
  sG4stage[2] = .04;
  sG4stage[3] = .05;
  sG4stage[4] = .07;
  topW[0] = 3.00; // minimum t value for range in top of wedge, relative to the
                  // center, both sides
  topW[1] = 3.75; // maximum t value for range in top of wedge, relative to the
                  // center, both sides
  brickW = 0.30;  // half width of range in t for brick-only protons
  emptyW = 2.0;   // half width of range in t for empty protons

  cout << "pCTcalib: reading the list of raw data file names from " << CalFile << endl;
  cout << "pCTcalib: the list of raw data file names for calibration is as "
          "follows:" << endl;
  int linecount = 0;
  string line;
  ifstream infile(CalFile);
  if (infile) {
    while (getline(infile, line)) {
      // cout << "From the file " << CalFile << " line " << linecount << ": " <<
      // line << endl;
      size_t found = line.find_first_not_of(" ");
      if (line[found] != '#') {
        size_t end = line.find_first_of(" \n", found);
        string vfn = line.substr(found, end);
        // cout << "Calibration file: '" << vfn << "'" << endl;
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
      } else {
        // The comment lines may contain values to override default calibration
        // parameters
        cout << line << endl;
        found = line.find_first_not_of("# ");
        line = line.substr(found);
        found = line.find_first_of(" =");
        string key = line.substr(0, found);
        found = line.find("=");
        if (found == line.npos) {
          cout << "No = sign found. There must be an = sign between the key "
                  "and value(s)" << endl;
          continue;
        }
        line = line.substr(found + 1);
        std::string::size_type sz = 0;
        int ival;
        float val;
        // cout << "The key is " << key << endl;
        if (key == "k1" || key == "K1") {
          for (int i = 0; i < nStage; i++) {
            line = line.substr(sz);
            try {
              ival = stoi(line, &sz);
            }
            catch (...) {
              break;
            }
            k1[i] = ival;
            cout << "   pCTcalib: setting calibration parameter k1[" << i << "]=" << k1[i] << endl;
          }
        } else if (key == "j1" || key == "J1") {
          for (int i = 0; i < nStage; i++) {
            line = line.substr(sz);
            try {
              ival = stoi(line, &sz);
            }
            catch (...) {
              break;
            }
            j1[i] = ival;
            cout << "   pCTcalib: setting calibration parameter j1[" << i << "]=" << j1[i] << endl;
          }
        } else if (key == "j2" || key == "J2") {
          for (int i = 0; i < nStage; i++) {
            line = line.substr(sz);
            try {
              ival = stoi(line, &sz);
            }
            catch (...) {
              break;
            }
            j2[i] = ival;
            cout << "   pCTcalib: setting calibration parameter j2[" << i << "]=" << j2[i] << endl;
          }
        } else if (key == "j3" || key == "J3") {
          for (int i = 0; i < nStage; i++) {
            line = line.substr(sz);
            try {
              ival = stoi(line, &sz);
            }
            catch (...) {
              break;
            }
            j3[i] = ival;
            cout << "   pCTcalib: setting calibration parameter j3[" << i << "]=" << j3[i] << endl;
          }
        } else if (key == "j4" || key == "J4") {
          for (int i = 0; i < nStage; i++) {
            line = line.substr(sz);
            try {
              ival = stoi(line, &sz);
            }
            catch (...) {
              break;
            }
            j4[i] = ival;
            cout << "   pCTcalib: setting calibration parameter j4[" << i << "]=" << j4[i] << endl;
          }
        } else if (key == "i1" || key == "I1") {
          for (int i = 0; i < nStage; i++) {
            line = line.substr(sz);
            try {
              ival = stoi(line, &sz);
            }
            catch (...) {
              break;
            }
            i1[i] = ival;
            cout << "   pCTcalib: setting calibration parameter i1[" << i << "]=" << i1[i] << endl;
          }
        } else if (key == "i2" || key == "I2") {
          for (int i = 0; i < nStage; i++) {
            line = line.substr(sz);
            try {
              ival = stoi(line, &sz);
            }
            catch (...) {
              break;
            }
            i2[i] = ival;
            cout << "   pCTcalib: setting calibration parameter i2[" << i << "]=" << i2[i] << endl;
          }
        } else if (key == "i3" || key == "I3") {
          line = line.substr(sz);
          i3 = stoi(line, &sz);
          cout << "   pCTcalib: setting calibration parameter i3=" << i3 << endl;
        } else if (key == "EG4stage" || key == "EG4STAGE") {
          for (int i = 0; i < nStage; i++) {
            line = line.substr(sz);
            try {
              val = stof(line, &sz);
            }
            catch (...) {
              break;
            }
            EG4stage[i] = val;
            cout << "   pCTcalib: setting calibration parameter EG4stage[" << i << "]=" << EG4stage[i] << endl;
          }
        } else if (key == "topW" || key == "TOPW" || key == "topw") {
          for (int i = 0; i < 2; i++) {
            line = line.substr(sz);
            try {
              val = stof(line, &sz);
            }
            catch (...) {
              break;
            }
            topW[i] = val;
            cout << "   pCTcalib: setting calibration parameter topW[" << i << "]=" << topW[i] << endl;
          }
        } else if (key == "brickW" || key == "BRICKW" || key == "brickw") {
          line = line.substr(sz);
          try {
            val = stof(line, &sz);
          }
          catch (...) {
            break;
          }
          brickW = val;
          cout << "   pCTcalib: setting calibration parameter brickW = " << brickW << endl;
        } else if (key == "emptyW" || key == "EMPTYW" || key == "emptyw") {
          line = line.substr(sz);
          try {
            val = stof(line, &sz);
          }
          catch (...) {
            break;
          }
          emptyW = val;
          cout << "   pCTcalib: setting calibration parameter emptyW = " << emptyW << endl;
        }
      }
    }
    if (linecount < 6) {
      cout << "Failed to find all 6 filename paths for the calibration; "
              "linecount=" << linecount << endl;
      return;
    }
  } else {
    cout << "Failed to open the list of calibration files in " << CalFile << endl;
    return;
  }
  infile.close();
  if (calFileNames.size() == 0)
    return;
  cout << "Set n_plotIn > 0 to produce lots of debug histograms." << endl;

  cout << endl;
  cout << "The beam particle type is " << partType << endl;
  float EG4tot = 0.;
  for (int i = 0; i < nStage; ++i) {
    cout << "MC derived stage energy for an empty event, for stage " << i << " is " << EG4stage[i] << endl;
    EG4tot += EG4stage[i];
  }
  cout << "MC derived total energy for an empty event is " << EG4tot << endl << endl;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      tPlaneShifts[i][j] = 0.0;
    }
  }

  for (int stage = 0; stage < nStage; ++stage) {
    Est[stage] = 0.0;
  }

  currentTime = time(NULL);
  now = localtime(&currentTime);

  Geometry = new pCTgeo(tWedgeOffset); // Create a class instance
                                       // with all of the geometry
                                       // information

  TVcal = new TVcorrection(TVcorrFile.c_str(), 0, 0, 0, 0); // The constants read in are not used in this
                                                            // round, but we need to pass the object to
                                                            // procEvt

  pCTgeo GeoInst = *Geometry;
  if (!redo) {
    procEvt = new EvtRecon(GeoInst, TVcal, calFileNames[0], Outputdir, max_events, max_time, n_debug, n_plot, 5,
                           useTemp, false, partType, pdstlr, realTimeCal,
                           OsName); // Process the raw data file to produce the event list
    Ut[0] = procEvt->uhitT[2];
    Ut[1] = procEvt->uhitT[3];
    Uv[0] = procEvt->uhitV[2];
    Uv[1] = procEvt->uhitV[3];
    cout << "pCTcalib: U values for T coordinates = " << Ut[0] << " and " << Ut[1] << endl;
    cout << "pCTcalib: U values for V coordinates = " << Uv[0] << " and " << Uv[1] << endl;
    for (int stage = 0; stage < nStage; ++stage) {
      cout << "Pedestal measured for stage " << stage << " is " << procEvt->Peds[stage] << " ADC counts " << endl;
    }
  }

} // end of the pCTcalib constructor

int pCTcalib::TVmapper() {
  if (redo) {
    cout << "pCTcalib::TVmapper, doing nothing in this run." << endl;
    return 0;
  }

  if (calFileNames.size() == 0) {
    cout << "TVmaper: only " << calFileNames.size() << " calibration file names found.  Need at least 1 for the TV map."
         << endl;
    return -2;
  }
  // Define a signal histogram for each pixel
  float bin0[nStage] = { 1000., 1000., 1000., 1000., 1000. };
  for (int Stage = 0; Stage < nStage; Stage++) {
    for (int pix = 0; pix < nPix; pix++) {
      int tPix = pix / 10;
      int vPix = pix % 10;
      string Title = "Stage " + to_string((long long int)Stage) + ", Pixel " + to_string((long long int)pix) +
                     "; tPix=" + to_string((long long int)tPix) + " vPix=" + to_string((long long int)vPix);
      if (partType == "H") {
        pxHistADC[Stage][pix] = new Histogram(2000, bin0[Stage], 4., Title, "ADC counts", "protons");
      } else {
        pxHistADC[Stage][pix] = new Histogram(2000, bin0[Stage], 8., Title, "ADC counts", "He ions");
      }
    }
  }
  Histogram2D hTVmap(150, -150., 2., 80, -40., 1., "V vs T distribution of tracks at energy detector", "T (mm)",
                     "V (mm)", "N");

  if (procEvt->useTmpFile)
    procEvt->reopenTmpFile();
  cout << "TVmapper: starting event loop" << endl;
  for (int EvtNum = 0; EvtNum < procEvt->nEvents; EvtNum++) {
    Event thisEvent;
    if (procEvt->useTmpFile)
      procEvt->readTmp(thisEvent);
    else
      thisEvent = procEvt->evtList[EvtNum];
    T[0] = thisEvent.Thit[2];
    V[0] = thisEvent.Vhit[2];
    T[1] = thisEvent.Thit[3];
    V[1] = thisEvent.Vhit[3];
    if (EvtNum % 1000000 == 0) {
      cout << "  Processing event " << EvtNum << endl;
      procEvt->dumpTmp(thisEvent);
    }

    for (int stage = 0; stage < nStage; stage++) {
      // Extrapolate the rear track vector to the energy detector stage
      double Tcorr = Geometry->extrap2D(Ut, T, Geometry->energyDetectorU(stage));
      double Vcorr = Geometry->extrap2D(Uv, V, Geometry->energyDetectorU(stage));
      if (fabs(Tcorr) > 189.0)
        continue;
      if (fabs(Vcorr) > 49.0)
        continue;

      int tPix = floor(0.1 * (Tcorr + 190.)); // calculate pixel indices in TV
                                              // correction maps
      int vPix = floor(0.1 * (Vcorr + 50.));
      int iPix = vPix + 10 * tPix;
      if (EvtNum % 1000000 == 0) {
        cout << "TVmapper: Stage " << stage << ", Tcorr=" << Tcorr << ", Vcorr=" << Vcorr << ", pixel number = ";
        cout << iPix << " U=" << Geometry->energyDetectorU(stage) << endl;
      }
      hTVmap.entry(Tcorr, Vcorr);
      if (iPix >= nPix || vPix > 9) {
        cout << "pCTcalib: iPix=" << iPix << " vPix=" << vPix << endl;
        perror("pCTcalib: invalid pixel");
        iPix = 0;
      }
      pxHistADC[stage][iPix]->entry(thisEvent.ADC[stage] - procEvt->Peds[stage]); // Histogram the raw
                                                                                  // ADC minus pedestal
                                                                                  // in each pixel of
                                                                                  // each stage
    }
  }

  // Save plots of some of the histograms as gnuplot files
  FILE *oFile;
  if (n_plot > 0) {
    for (int stage = 0; stage < nStage; stage++) {
      for (int hst = 0; hst < nPix / 9 + 1; ++hst) {
        string fileName =
            Outputdir + "/ADC_Stage_" + to_string((long long int)stage) + "_" + to_string((long long int)hst) + ".gp";
        oFile = fopen(fileName.c_str(), "w");
        if (oFile != NULL) {
          fprintf(oFile, setTerm.c_str());
          fprintf(oFile, "set multiplot title 'Stage %d Signal Region for file "
                         "%s' layout 3,3 columnsfirst scale 1.0,1.0\n",
                  stage, calFileNames[0].c_str());
          for (int j = 0; j < 9; j++) {
            int pix = hst * 9 + j;
            if (pix < nPix)
              pxHistADC[stage][pix]->plot(oFile);
          }
          fprintf(oFile, "unset multiplot\n");
          fprintf(oFile, "show label\n");
          fprintf(oFile, "unset label\n");
          fclose(oFile);
        } else {
          cout << "TVmapper: unable to open the ADC histogram file " << fileName << endl;
        }
      }
    }
  }
  hTVmap.stats();
  string fileName = Outputdir + "/track_VvsT.gp";
  oFile = fopen(fileName.c_str(), "w");
  if (oFile != NULL) {
    hTVmap.plot(oFile);
  }

  cout << "TVmapper: Starting evaluation of the TV correction factors." << endl;

  for (int pix = 0; pix < nPix; pix++) {
    cout << "TVmapper pixel stage ADC values:  " << pix;
    for (int stage = 0; stage < nStage; stage++) {
      float xLow, xHigh;
      int ret =
          pxHistADC[stage][pix]->FWHMboundaries(xLow, xHigh); // Find the full-width-half-max upper and lower bounds
      if (ret == 0) {
        float Sadc = pxHistADC[stage][pix]->mean(xLow, xHigh); // Evaluate the histogram mean between those bounds
        TVmap[stage][pix] = EG4stage[stage] / Sadc;            // Correction factor and conversion to MeV
        cout << "  stg" << stage << "=" << Sadc;
      } else {
        TVmap[stage][pix] = 0.999;
        cout << "  stg" << stage << "=       ";
      }
      delete pxHistADC[stage][pix]; // Free up the memory that was sucked up by
                                    // the pixel histograms
    }
    cout << endl;
  }

  cout << "TVmapper: Finished with mapping the TV calibration" << endl;
  return 0;
}

void pCTcalib::enrgDep() { // Calculate energy depositions in each stage and
                           // their sum using new TV maps, to check that
                           // calibration worked
  if (redo) {
    cout << "pCTcalib::enrgDep, doing nothing in this run." << endl;
    return;
  }
  for (int stage = 0; stage < nStage; ++stage) {
    string Title = "pCTcalib::enrgDep Energy of stage " + to_string((long long int)stage);
    if (partType == "H")
      stgHistE[stage] = new Histogram(400, 15., 0.175, Title, "E (MeV)", "protons");
    else
      stgHistE[stage] = new Histogram(400, 15., 0.7, Title, "E (MeV)", "He ions");
    Title = "pCTcalib::enrgDep Energy of stage " + to_string((long long int)stage) + " vs T";
    stgEvsT[stage] = new ProfilePlot(100, -150., 3.0, Title, "T (mm)", "E (MeV)");
    Title = "pCTcalib::enrgDep Energy of stage " + to_string((long long int)stage) + " vs V";
    stgEvsV[stage] = new ProfilePlot(100, -50., 1.0, Title, "V (mm)", "E (MeV)");
  }
  if (partType == "H")
    EsumH = new Histogram(900, 0., 0.25, "pCTcalib::enrgDep Sum of stage energies", "Etot (MeV)", "protons");
  else
    EsumH = new Histogram(900, 0., 1.0, "pCTcalib::enrgDep Sum of stage energies", "Etot (MeV)", "He ions");
  for (int iPix = 0; iPix < nPix; iPix++) {
    for (int stage = 0; stage < nStage; ++stage) {
      string Title = "pCTcalib::enrgDep Energy of stage " + to_string((long long int)stage) + " for pixel " +
                     to_string((long long int)iPix);
      if (partType == "H")
        pxHistE[stage][iPix] = new Histogram(200, 15., .35, Title, "E (MeV)", "protons");
      else
        pxHistE[stage][iPix] = new Histogram(200, 60., 1.4, Title, "E (MeV)", "He ions");
    }
  }

  float mxEstage, mnEstage;
  if (partType == "H") {
    mnEstage = 15.;
    mxEstage = 100.;
  } else {
    mnEstage = 60.;
    mxEstage = 400.;
  }
  if (procEvt->useTmpFile)
    procEvt->rewindTmpFile();
  for (int EvtNum = 0; EvtNum < procEvt->nEvents; EvtNum++) {
    if (EvtNum % 1000000 == 0)
      cout << "pCTcalib::enrgDep, Processing event " << EvtNum << endl;
    Event thisEvent;
    if (procEvt->useTmpFile)
      procEvt->readTmp(thisEvent);
    else
      thisEvent = procEvt->evtList[EvtNum];
    T[0] = thisEvent.Thit[2];
    V[0] = thisEvent.Vhit[2];
    T[1] = thisEvent.Thit[3];
    V[1] = thisEvent.Vhit[3];

    float Esum = 0.;
    for (int stage = 0; stage < nStage; stage++) {
      // Extrapolate the rear track vector to the energy detector stage
      float Tcorr = Geometry->extrap2D(Ut, T, Geometry->energyDetectorU(stage));
      float Vcorr = Geometry->extrap2D(Uv, V, Geometry->energyDetectorU(stage));
      if (fabs(Tcorr) > 150.0)
        goto nextEvt;
      if (fabs(Vcorr) > 45.0)
        goto nextEvt;

      bool inBounds;
      float Ecorr = ((float)thisEvent.ADC[stage] - procEvt->Peds[stage]) *
                    TVcal->corrFactorInt(TVmap, stage, Tcorr, Vcorr, inBounds);
      if (!inBounds)
        goto nextEvt;
      int tPix = floor(0.1 * (Tcorr + 190.));
      int vPix = floor(0.1 * (Vcorr + 50.));
      int iPix = vPix + 10 * tPix;
      if (iPix >= nPix)
        goto nextEvt;
      if (vPix > 9) {
        cout << "pCTcalib: iPix=" << iPix << " vPix=" << vPix << endl;
        perror("pCTcalib: invalid pixel");
        iPix = 0;
      }
      pxHistE[stage][iPix]->entry(Ecorr);
      stgHistE[stage]->entry(Ecorr);
      if (Ecorr > mnEstage && Ecorr < mxEstage) {
        if (fabs(Vcorr) < 35.0)
          stgEvsT[stage]->entry(Tcorr, Ecorr);
        if (fabs(Tcorr) < 140.0)
          stgEvsV[stage]->entry(Vcorr, Ecorr);
      }
      Esum = Esum + Ecorr;
    }
    EsumH->entry(Esum);
  nextEvt:
    ;
  }
  cout << "pCTcalib::enrgDep, finished the loop over " << procEvt->nEvents << " events." << endl;

  for (int iPix = 0; iPix < nPix; ++iPix) {
    cout << "Corrected stage energies for pixel " << iPix;
    for (int stage = 0; stage < nStage; ++stage) {
      float xLow, xHigh;
      int ret = pxHistE[stage][iPix]->FWHMboundaries(xLow, xHigh);
      if (ret == 0) {
        float Enrg = pxHistE[stage][iPix]->mean(xLow, xHigh);
        cout << "  stg" << stage << "=" << Enrg;
      } else {
        cout << "  stg" << stage << "=      ";
      }
      delete pxHistE[stage][iPix];
    }
    cout << endl;
  }

  string fileName = Outputdir + "/StageEnergy.gp";
  FILE *oFile = fopen(fileName.c_str(), "w");
  if (oFile != NULL) {
    fprintf(oFile, setTerm.c_str());
    fprintf(oFile, "set multiplot title 'Calibrated Energy Detector Response "
                   "for file %s' layout 3,2 columnsfirst scale 1.0,1.0\n",
            calFileNames[0].c_str());
    stgHistE[0]->plot(oFile);
    stgHistE[1]->plot(oFile);
    stgHistE[2]->plot(oFile);
    stgHistE[3]->plot(oFile);
    fprintf(oFile, "set label 4 'Current date and time: %d-%d-%d    %d:%d' at "
                   "screen 0.55, 0.28\n",
            now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min);
    fprintf(oFile, "set label 5 'Run Number %d' at screen 0.55, 0.25 left\n", procEvt->runNumber);
    fprintf(oFile, "set label 6 'FPGA Program Version %d' at screen 0.55, 0.22 left\n", procEvt->program_version);
    fprintf(oFile, "set label 7 'Stage Angle = %5.1f degrees' at screen 0.55, 0.19 left\n", procEvt->stage_angle);
    fprintf(oFile, "set label 8 'Start Time = %s ' at screen 0.55, 0.16 left\n", procEvt->runStartTime.c_str());
    fprintf(oFile, "set label 9 'Number of good events = %d' at screen 0.55, 0.13 left\n", procEvt->nEvents);
    stgHistE[4]->plot(oFile);
    fprintf(oFile, "unset multiplot\n");
    fprintf(oFile, "show label\n");
    fprintf(oFile, "unset label\n");
    fclose(oFile);
  } else {
    cout << "pCTcalib::enrgDep Unable to open the stage energy histogram file " << fileName << endl;
  }

  fileName = Outputdir + "/StageEvsT.gp";
  oFile = fopen(fileName.c_str(), "w");
  if (oFile != NULL) {
    fprintf(oFile, setTerm.c_str());
    fprintf(oFile, "set multiplot title 'Calibrated Energy Detector Response "
                   "vs T for file %s' layout 3,2 columnsfirst scale 1.0,1.0\n",
            calFileNames[0].c_str());
    stgEvsT[0]->plot(oFile);
    stgEvsT[1]->plot(oFile);
    stgEvsT[2]->plot(oFile);
    stgEvsT[3]->plot(oFile);
    fprintf(oFile, "set label 4 'Current date and time: %d-%d-%d    %d:%d' at "
                   "screen 0.55, 0.28\n",
            now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min);
    fprintf(oFile, "set label 5 'Run Number %d' at screen 0.55, 0.25 left\n", procEvt->runNumber);
    fprintf(oFile, "set label 6 'FPGA Program Version %d' at screen 0.55, 0.22 left\n", procEvt->program_version);
    fprintf(oFile, "set label 7 'Stage Angle = %5.1f degrees' at screen 0.55, 0.19 left\n", procEvt->stage_angle);
    fprintf(oFile, "set label 8 'Start Time = %s ' at screen 0.55, 0.16 left\n", procEvt->runStartTime.c_str());
    fprintf(oFile, "set label 9 'Number of good events = %d' at screen 0.55, 0.13 left\n", procEvt->nEvents);
    stgEvsT[4]->plot(oFile);
    fprintf(oFile, "unset multiplot\n");
    fprintf(oFile, "show label\n");
    fprintf(oFile, "unset label\n");
    fclose(oFile);
  } else {
    cout << "pCTcalib::enrgDep Unable to open the stage energy histogram file " << fileName << endl;
  }

  fileName = Outputdir + "/StageEvsV.gp";
  oFile = fopen(fileName.c_str(), "w");
  if (oFile != NULL) {
    fprintf(oFile, setTerm.c_str());
    fprintf(oFile, "set multiplot title 'Calibrated Energy Detector Response "
                   "vs V for file %s' layout 3,2 columnsfirst scale 1.0,1.0\n",
            calFileNames[0].c_str());
    stgEvsV[0]->plot(oFile);
    stgEvsV[1]->plot(oFile);
    stgEvsV[2]->plot(oFile);
    stgEvsV[3]->plot(oFile);
    fprintf(oFile, "set label 4 'Current date and time: %d-%d-%d    %d:%d' at "
                   "screen 0.55, 0.28\n",
            now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min);
    fprintf(oFile, "set label 5 'Run Number %d' at screen 0.55, 0.25 left\n", procEvt->runNumber);
    fprintf(oFile, "set label 6 'FPGA Program Version %d' at screen 0.55, 0.22 left\n", procEvt->program_version);
    fprintf(oFile, "set label 7 'Stage Angle = %5.1f degrees' at screen 0.55, 0.19 left\n", procEvt->stage_angle);
    fprintf(oFile, "set label 8 'Start Time = %s ' at screen 0.55, 0.16 left\n", procEvt->runStartTime.c_str());
    fprintf(oFile, "set label 9 'Number of good events = %d' at screen 0.55, 0.13 left\n", procEvt->nEvents);
    stgEvsV[4]->plot(oFile);
    fprintf(oFile, "unset multiplot\n");
    fprintf(oFile, "show label\n");
    fprintf(oFile, "unset label\n");
    fclose(oFile);
  } else {
    cout << "pCTcalib::enrgDep Unable to open the stage energy histogram file " << fileName << endl;
  }

  fileName = Outputdir + "/TotalEnergy.gp";
  oFile = fopen(fileName.c_str(), "w");
  if (oFile != NULL) {
    EsumH->plot(oFile);
  } else {
    cout << "pCTcalib::enrgDep Unable to open the total energy histogram file " << fileName << endl;
  }

  for (int stage = 0; stage < nStage; ++stage) {
    float xLow, xHigh;
    int ret = stgHistE[stage]->FWHMboundaries(xLow, xHigh);
    float Sadc;
    if (ret == 0) {
      Sadc = stgHistE[stage]->mean(xLow, xHigh);
    } else {
      Sadc = 0.;
      cout << "pCTcalib::enrgDep, no data for stage " << stage << endl;
    }
    Est[stage] = Sadc;
    cout << "pCTcalib::enrgDep, Energy in stage " << stage << " = " << Est[stage]
         << "; dE=" << Est[stage] - EG4stage[stage] << " MeV" << endl;

    delete stgHistE[stage];
    delete stgEvsT[stage];
    delete stgEvsV[stage];
  }
  cout << "pCTcalib::enrgDep, MC - experiment difference |dE| should be <0.1 "
          "MeV ! " << endl;

  float xLow, xHigh;
  int ret = EsumH->FWHMboundaries(xLow, xHigh);
  if (ret == 0) {
    EnS = EsumH->mean(xLow, xHigh);
    cout << "pCTcalib::enrgDep, Energy sum = " << EnS << endl;
  } else {
    EnS = 0.;
    cout << "pCTcalib::enrgDep Failed to find the peak of the energy "
            "distribution" << endl;
  }
  delete EsumH;
  cout << "pCTcalib::enrgDep,  Done with the TV calibration task" << endl;
}

void pCTcalib::writeTVfile() {
  if (redo) {
    cout << "pCTcalib::writeTVfile, doing nothing in this run." << endl;
    return;
  }
  cout << "pCTcalib::writeTVfile, the TV calibration file will be opened in " << Outputdir << endl;
  string TVfile = Outputdir + "/" + TVcorrFile;
  cout << "Opening TV correction file: " << TVfile << " for output." << endl;

  ofstream TVcalfile;
  TVcalfile.open(TVfile, ios::out | ios::trunc); // text file for TV-corr data
  if (!TVcalfile.is_open()) {
    cout << "Error opening TV correction output file " << TVfile << endl;
    cout << "Trying to open " << TVcorrFile << " instead. . ." << endl;
    TVcalfile.open(TVcorrFile, ios::out | ios::trunc);
    if (TVcalfile.is_open()) {
      cout << "Successful file open!" << endl;
      TVfile = TVcorrFile;
    } else {
      cout << "Rats!  That didn't work either." << endl;
      exit(1);
    }
  }
  TVcorrFile = TVfile; // Update this path to be sure that Wcalib uses the same file.

  for (int stage = 0; stage < nStage; ++stage) {
    for (int j = 0; j < nPix / 10; ++j) {
      for (int k = 0; k < 10; ++k)
        TVcalfile << TVmap[stage][10 * j + k] << " ";
      TVcalfile << endl;
    }
  }
  for (int stage = 0; stage < nStage; ++stage)
    TVcalfile << procEvt->Peds[stage] << " ";
  TVcalfile << endl;
  for (int stage = 0; stage < nStage; ++stage)
    TVcalfile << Est[stage] << " ";
  TVcalfile << EnS << endl;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; ++j)
      TVcalfile << tPlaneShifts[i][j] << " ";
    TVcalfile << endl;
  }
  TVcalfile << "# particle = " << partType << endl;
  TVcalfile << "# Valid date and run ranges for this calibration:" << endl;
  TVcalfile << "# minDate = " << minDate << endl;
  TVcalfile << "# maxDate = " << maxDate << endl;
  TVcalfile << "# minRun = " << minRun << endl;
  TVcalfile << "# maxRun = " << maxRun << endl;
  TVcalfile << "# Date and time of this calibration run: " << now->tm_year + 1900 << "-" << now->tm_mon + 1 << "-"
            << now->tm_mday << "    " << now->tm_hour << ":" << now->tm_min << endl;
  TVcalfile << "# This file was originally created in " << Outputdir
            << " directory. It contains TV correction maps evaluated" << endl;
  TVcalfile << "# using the 'empty' raw data file " << calFileNames[0] << endl;
  TVcalfile << "# The pre-processing code version number is [" << programVersion << "]" << endl;

  TVcalfile.close();

  cout << "writeTVfile: done with outputting the new TV calibration "
          "corrections in " << TVfile << endl;

  delete procEvt;
};

int pCTcalib::Wcalib(float Thresholds[nStage], int pdstlr[nStage], int n_threads) {

  cout << endl << "Entering pCTcalib::Wcalib to execute the WEPL calibration process" << endl;

  if (calFileNames.size() != 6) {
    cout << "Only " << calFileNames.size() << " calibration raw data file names found.  Need all 6 to include "
                                              "the WEPL calibration." << endl;
    return -1;
  }
  for (int nBricks = 0; nBricks < 5; nBricks++) {
    cout << "The raw data file for " << nBricks << " bricks is " << calFileNames[nBricks] << endl;
  }

  cout << "pCTcalib::Wcalib: we will read TV correction constants from " << TVcorrFile << endl;
  TVcorrection TVcal(TVcorrFile.c_str(), 0, 0, 0, 0);

  struct {
    Histogram2D *h[nStage];
  } REhist[5];
  struct {
    Histogram2D *h[nStage];
  } dEEhist[5];
  for (int nBricks = 0; nBricks < 5; ++nBricks) {
    for (int stage = 0; stage < nStage; ++stage) {
      float enrgB0 = 0.; // The analysis of energy slices further down assumes
                         // that the first energy bin starts at zero
      string Title = "WEPL calibration array for stage " + to_string((long long int)stage) + " and nBricks=" +
                     to_string((long long int)nBricks);
      string TitledEE = "dE-E spectra for stage" + to_string((long long int)stage) + " and nBricks=" +
                        to_string((long long int)nBricks);
      dEEhist[nBricks].h[stage] = new Histogram2D(nEnrg, 0., EnergyBinWidth, nEnrg, 0., EnergyBinWidth, TitledEE,
                                                  "Energy dE stage (MeV)", "Energy Bragg-peak stage (MeV)", "counts");

      if (partType == "H")
        REhist[nBricks].h[stage] = new Histogram2D(nRange, 0., RangeBinWidth, nEnrg, enrgB0, EnergyBinWidth, Title,
                                                   "proton range (mm)", "energy (MeV)", "counts");
      else
        REhist[nBricks].h[stage] = new Histogram2D(nRange, 0., RangeBinWidth, nEnrg, enrgB0, EnergyBinWidth, Title,
                                                   "He range (mm)", "energy (MeV)", "counts");
    }
  }
  calStuff stuff;
  stuff.Outputdir = Outputdir; // This structure is just used to reduce the
                               // number of arguments, which is limited for
                               // spawning threads, at least in Windows.
  stuff.max_events = max_events;
  stuff.max_time = max_time;
  stuff.partType = partType;
  stuff.OsName = OsName;
  stuff.topW1 = topW[0];
  stuff.topW2 = topW[1];
  stuff.brickW = brickW;
  stuff.emptyW = emptyW;
  stuff.reCalibrate = realTimeCal;
  for (int stage = 0; stage < nStage; stage++) {
    stuff.Thr[stage] = Thresholds[stage];
    stuff.pdstlr[stage] = pdstlr[stage];
  }
  stuff.useTemp = useTempFile;
  if (!redo) {
    // Submit the event processing in 5 either sequential or parallel threads
    pCTgeo GeoInst = *Geometry;
    vector<std::thread> thr; // Vector of thread identifiers

    for (int nBricks = 1; nBricks < 5; nBricks++) {
      stuff.inputFileName = calFileNames[nBricks + 1];
      stuff.Nbricks = nBricks;
      if (n_threads <= 1)
        procWEPLcal(stuff, TVcal, REhist[nBricks].h, dEEhist[nBricks].h, GeoInst);
      else
        thr.push_back(thread(procWEPLcal, stuff, TVcal, REhist[nBricks].h, dEEhist[nBricks].h, GeoInst));
    }
    stuff.inputFileName = calFileNames[1];
    stuff.Nbricks = 0;
    procWEPLcal(stuff, TVcal, REhist[0].h, dEEhist[0].h, GeoInst);
    if (n_threads > 1) {
      for (int n = 0; n < thr.size(); ++n) {
        thr[n].join(); // Waiting here for all the daughter threads to finish.
        cout << "pCTcalib::Wcalib, joined thread number " << n << endl;
      }
    }

    // Debug printout of some histogram projections
    if (n_plot > 0) {
      for (int nBricks = 0; nBricks < 5; ++nBricks) {
        string fn = Outputdir + "/Eproj_" + to_string((long long int)nBricks) + ".gp";
        FILE *oFile = fopen(fn.c_str(), "w");
        if (oFile != NULL) {
          fprintf(oFile, setTerm.c_str());
          fprintf(oFile, "set multiplot title 'pCTcalib::Wcalib: Energy "
                         "projections for %d bricks' layout 3,2 columnsfirst "
                         "scale 1.0,1.0\n",
                  nBricks);
          for (int stage = 0; stage < nStage; ++stage) {
            Histogram Eproj = REhist[nBricks].h[stage]->yProj();
            Eproj.plot(oFile);
          }
          fprintf(oFile, "unset multiplot\n");
          fprintf(oFile, "show label\n");
          fprintf(oFile, "unset label\n");
          fclose(oFile);
        }
      }
    }

    // Write out the individual 2-D histograms, so that they might be reanalyzed
    // later
    for (int stage = 0; stage < nStage; ++stage) {
      for (int nBricks = 0; nBricks < 5; ++nBricks) {
        string fn = Outputdir + "/calMap2D_" + to_string((long long int)stage) + "_B" +
                    to_string((long long int)nBricks) + ".gp";
        FILE *oFile = fopen(fn.c_str(), "w");
        if (oFile != NULL) {
          REhist[nBricks].h[stage]->plot(oFile);
          fclose(oFile);
        }
      }
    }
  } else {
    cout << "pCTcalib::Wcalib: redoing analysis of calibration histograms "
            "already generated." << endl;
    for (int stage = 0; stage < nStage; ++stage) {
      for (int nBricks = 0; nBricks < 5; ++nBricks) {
        string fn = Outputdir + "/calMap2D_" + to_string((long long int)stage) + "_B" +
                    to_string((long long int)nBricks) + ".gp";
        REhist[nBricks].h[stage]->read(fn);
      }
    }
  }

  cout << "pCTcalib::Wcalib, begin adding together the range-energy tables "
          "from the runs with different numbers of bricks." << endl;

  // Combine the resulting maps into one map for each stage by simply adding
  // them
  for (int stage = 0; stage < nStage; ++stage) {
    for (int nBricks = 1; nBricks < 5; ++nBricks) {
      string Title = "Summed WEPL calibration array for stage " + to_string((long long int)stage);
      cout << "pCTcalib::Wcalib, stage " << stage << " add in the histogram for " << nBricks << " bricks" << endl;
      REhist[0].h[stage]->add(REhist[nBricks].h[stage], Title);
      delete REhist[nBricks].h[stage];
      dEEhist[0].h[stage]->add(dEEhist[nBricks].h[stage], Title);
      delete dEEhist[nBricks].h[stage];
    }
    string fn = Outputdir + "/calMap2Dorig_" + to_string((long long int)stage) + ".gp";
    FILE *oFile = fopen(fn.c_str(), "w");
    if (oFile != NULL) {
      REhist[0].h[stage]->plot(oFile);
      fclose(oFile);
    } else
      cout << "Unable to open the stage energy histogram file " << fn << endl;

    // Normalize all the columns to have equal total counts (equal counts for
    // all ranges)
    // However, we have to take care not to amplify nearly empty columns (i.e.
    // noise)
    if (normalise) {    // DEFAULT IS NO
      int bMin = 0;     // 20;  // Only count the area between these bin limits in energy
      int bMax = nEnrg; //-20
      Histogram hYproj = REhist[0].h[stage]->yProj(bMin, bMax);
      double area = hYproj.max();
      cout << "pCTcalib::Wcalib: normalizing columns in the 2D WET vs E plot "
              "for stage " << stage << ", max area=" << area << endl;
      double threshold = area / 20.;
      REhist[0].h[stage]->normalizeColumns(bMin, bMax, area, threshold);
    }

    fn = Outputdir + "/calMap2D_" + to_string((long long int)stage) + ".gp";
    oFile = fopen(fn.c_str(), "w");
    if (oFile != NULL) {
      REhist[0].h[stage]->plot(oFile);
      fclose(oFile);
    } else
      cout << "Unable to open the stage energy histogram file " << fn << endl;
    fn = Outputdir + "/dEEMap2D_" + to_string((long long int)stage) + ".gp";
    oFile = fopen(fn.c_str(), "w");
    if (oFile != NULL) {
      dEEhist[0].h[stage]->plot(oFile);
      fclose(oFile);
    } else
      cout << "Unable to open the stage dE-E histogram file " << fn << endl;
  }

  float Rst[nStage][nEnrg] = { 0. }; // 5 arrays to store range vs E
  float Sst[nStage][nEnrg] = { 0. }; // Peak width for each energy (sigma)
  float est[nEnrg] = { 0. };         // Corresponding E array in MeV/4

  // Analyze energy slices of the summed maps, and plot the slices
  for (int stage = 0; stage < nStage; ++stage) {
    int nrg = 0;
    for (int hst = 0; hst < nEnrg / 6 + 1; ++hst) { // Group by 6 slices, to plot 6 histograms per page
      string fn =
          Outputdir + "/Slices_Stage" + to_string((long long int)stage) + "_" + to_string((long long int)hst) + ".gp";
      FILE *oFile;
      if (n_plot > 0) {
        oFile = fopen(fn.c_str(), "w");
        if (oFile != NULL) {
          fprintf(oFile, setTerm.c_str());
          fprintf(oFile, "set multiplot title 'Slices at constant energy for "
                         "stage %d' layout 3,2 columnsfirst scale 1.0,1.0\n",
                  stage);
        }
      }
      for (int j = 0; j < 6; ++j) {
        if (nrg >= nEnrg)
          break;
        Histogram Slice = REhist[0].h[stage]->rowSlice(nrg);
        float xLow, xHigh;
        int ret = Slice.FWHMboundaries(xLow, xHigh);
        float Sadc;
        if (ret == 0) {
          Sadc = Slice.mean(xLow, xHigh);
        } else {
          if (stage == 4) {
            if (Slice.imode() == 0)
              Sadc = 0.; // Up against the wall at zero range
            else {
              cout << ret << " Cannot find the FWHM bounds for stage " << stage << ", slice " << nrg << endl;
              float mode = Slice.mode();
              Sadc = Slice.mean(mode - 2. * RangeBinWidth, mode + 2 * RangeBinWidth);
            }
          } else {
            cout << ret << " Cannot find the FWHM bounds for stage " << stage << ", slice " << nrg << endl;
            float mode = Slice.mode();
            Sadc = Slice.mean(mode - 2. * RangeBinWidth, mode + 2 * RangeBinWidth);
          }
        }
        Rst[stage][nrg] = Sadc;
        Sst[stage][nrg] = (xHigh - xLow) / 2.2;
        if (oFile != NULL && n_plot > 0)
          Slice.plot(oFile);
        nrg++;
      }
      if (oFile != NULL && n_plot > 0) {
        fprintf(oFile, "unset multiplot\n");
        fclose(oFile);
      }
    }
  }

  // Plot the results before corrections
  if (n_plot > 0) {
    float rngB[nEnrg];
    for (int nrg = 0; nrg < nEnrg; ++nrg)
      rngB[nrg] = float(nrg) * EnergyBinWidth;
    for (int stage = 0; stage < nStage; ++stage) {
      string Title = "Range vs Energy for Stage " + to_string((long long int)stage);
      string fn = Outputdir + "/rngEnrgUncorr_" + to_string((long long int)stage) + ".gp";
      plot2D(fn, Title, "Energy (MeV)", "Range (mm)", nEnrg, rngB, Rst[stage], Sst[stage]);
    }
  }

  cout << "Wcalib: interpolation parameters:" << endl;
  cout << "        k1 is for extrapolating below the threshold" << endl;
  cout << "        j1 through j4 are for interpolating through a kink region "
          "(if there is one)" << endl;
  cout << "        i1 and i2 are for interpolating between stages or, for the "
          "last stage, extrapolating to negative WET" << endl;
  cout << "        These parameters can be set by comments in the WEPL "
          "calibration file." << endl;
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

  cout << "Wcalib: Interpolating the calibration curves through the kink region" << endl;
  for (int stage = 0; stage < 4; ++stage) {
    if (false) { //(j1[stage] > nEnrg) {
      cout << "Wcalib: skipping kink-region interpolation for stage " << stage << " (no kink)." << endl;
      continue;
    }
    quadFit qF;
    for (int i = j1[stage]; i < j2[stage]; i++) {
      qF.addPnt(float(i) + 0.5, Rst[stage][i]);
    }
    for (int i = j3[stage]; i < j4[stage]; i++) {
      qF.addPnt(float(i) + 0.5, Rst[stage][i]);
    }
    if (qF.solve() == 0) {
      for (int i = j2[stage]; i < j3[stage]; i++) {
        Rst[stage][i] = qF.eval(float(i) + 0.5);
      }
    } else
      cout << "***************WARNING***************: Could not fit kink "
              "region for stage " << stage << "!" << endl;
  }

  // approximate the curve for stage 4 near range R = 0  with a quadratic spline
  // to exclude R(E) distortion due to range R straggling in the MSS and allow
  // negative range values
  // for correct air RSP calculations in reconstruction.
  // Fill calibration curve bins above 60 MeV (for protons) using quadratic
  // spline
  // Also, fill in the high energy ends of the other curves with fitted
  // quadratic curves, although in general that
  // energy range won't be used.
  double x0 = 4. * TVcal.Eempt[4]; // Sadc / Sweight; // mean energy deposition
                                   // for air (x 4)

  qSpline Spl(&Rst[4][0], i1[4], i2[4], x0); // Quadratic spline extrapolation for the last stage
  for (int i = i2[4] + 1; i < nEnrg; i++) {
    Rst[4][i] = Spl.eval(float(i) + 0.5);
  }

  for (int stage = 0; stage < 5; stage++) { // Don't allow the extrapolations to curve back upwards
    int iMin = i2[stage];
    float rMin = Rst[stage][i2[stage]];
    for (int i = i2[stage]; i < nEnrg; i++) {
      if (Rst[stage][i] < rMin) {
        rMin = Rst[stage][i];
        iMin = i;
      }
    }
    for (int i = iMin; i < nEnrg; i++) {
      Rst[stage][i] = rMin;
    }
  }

  cout << "Corrected stage 4 air W=0. E, dE: " << 0.25 * x0 << "  " << Est[4] - x0 / 4 << endl;

  // correct stage 0 low energy part for lack of data due to high threshold (up
  // to 13 MeV Aug 2016,  18 MeV Oct 2016)
  // using experimentally measured R(E) dependance approximated with pol2
  // (quadratic) function Rend.
  // First, find where the calibration data start to look reasonable
  int kGood = 79;
  for (int k = 0; k < nEnrg - 7; k++) {
    float thisVal = Rst[0][k];
    if (thisVal == 0.)
      continue;
    float diff;
    for (int n = k + 1; n < k + 8; n++) {
      diff = (Rst[0][n] - thisVal) / thisVal;
      if (abs(diff) > 0.05)
        goto nextVal;
    }
    diff = (Rst[0][k + 1] - thisVal) / thisVal;
    if (abs(diff) > 0.005)
      continue;
    kGood = k;
    cout << "pCTcalib: first good calibration point for stage 0 is at k=" << kGood << endl;
    goto foundIt;
  nextVal:
    ;
  }
  cout << "pCTcalib: failed to find a good starting value for stage 0.  "
          "Setting kGood to " << kGood << endl;
foundIt:
  for (int k = 0; k < kGood; k++) {
    if (partType == "H")
      Rst[0][k] =
          Rst[0][kGood] + Rend((kGood - k) * EnergyBinWidth) - float(kGood - k) * (52.16 - 51.12) / float(kGood);
    else
      Rst[0][k] =
          Rst[0][kGood] + RendHe((kGood - k) * EnergyBinWidth) - float(kGood - k) * (52.39 - 51.12) / float(kGood);
  }

  cout << "Thickness of the stages 0 - 3 (should be 51.12 +/- 1mm) :" << endl;
  cout << Rst[0][0] - Rst[1][0] << " " << Rst[1][0] - Rst[2][0] << " " << Rst[2][0] - Rst[3][0] << " "
       << Rst[3][0] - Rst[4][0] << endl;

  // Plot the final calibration results
  float rngB[nEnrg];
  for (int nrg = 0; nrg < nEnrg; ++nrg)
    rngB[nrg] = float(nrg) * EnergyBinWidth;
  for (int stage = 0; stage < nStage; ++stage) {
    string Title = "Range vs Energy for Stage " + to_string((long long int)stage);
    string fn = Outputdir + "/rngEnrg_" + to_string((long long int)stage) + ".gp";
    plot2D(fn, Title, "Energy (MeV)", "Range (mm)", nEnrg, rngB, Rst[stage], Sst[stage]);
  }

  /// Lennart Volz, November 2018
  /// dE-E parameter evaluation:
  int Estep[3] = { 240, 120, 12 }; // work in the same way for helium and
                                   // proton, but could also be set manually for
                                   // optimization
  float E[3] = { Estep[0] * EnergyBinWidth, Estep[1] * EnergyBinWidth, Estep[2] * EnergyBinWidth };
  float dEElow[5][3];
  float dEEhigh[5][3];
  FILE oFile;
  string fn;
  for (int stage = 1; stage < 5; stage++) {
    float xlow[3], xhigh[3];
    for (int j = 0; j < 3; j++) {
      Histogram dEESlice =
          dEEhist[0].h[stage]->rowSlice(Estep[j]); // gets the dE stage distribution for a given energy E
      int ret = dEESlice.FWHMboundaries(xlow[j], xhigh[j]);
      xlow[j] = xlow[j] - (xhigh[j] - xlow[j]) * 0.7848;   // extend to 3 sigma region instead of 1 FWHM
      xhigh[j] = xhigh[j] + (xhigh[j] - xlow[j]) * 0.7848; // might be more optimal to use 2.5 sigma

      if (n_plot > 0) {
        string fn =
            Outputdir + "/dEESlice_" + to_string((long long int)stage) + "_" + to_string((long long int)j) + ".gp";
        FILE *oFile = fopen(fn.c_str(), "w");
        dEESlice.plot(oFile);
        fclose(oFile);
      }
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

  /// Write the Rst0 through Rst4 range vs energy tables to the calibration text
  /// file
  string Wfile = Outputdir + "/" + WcalibFile;
  cout << "Opening output calibration file: " << Wfile << endl;
  ofstream tWcalfile;
  tWcalfile.open(Wfile, ios::out | ios::trunc); // text file R vs E tables
  if (!tWcalfile.is_open()) {
    cout << "Error opening W correction output file " << Wfile << endl;
    cout << "Trying to open " << WcalibFile << " instead. . ." << endl;
    tWcalfile.open(WcalibFile, ios::out | ios::trunc);
    if (tWcalfile.is_open())
      cout << "Successful file open!" << endl;
    else {
      cout << "Rats!  That didn't work either." << endl;
      exit(1);
    }
  }
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k)
      tWcalfile << Rst[0][k + i * 10] << " ";
    tWcalfile << endl;
  }
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k)
      tWcalfile << Rst[1][k + i * 10] << " ";
    tWcalfile << endl;
  }
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k)
      tWcalfile << Rst[2][k + i * 10] << " ";
    tWcalfile << endl;
  }
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k)
      tWcalfile << Rst[3][k + i * 10] << " ";
    tWcalfile << endl;
  }
  for (int i = 0; i < nEnrg / 10; ++i) {
    for (int k = 0; k < 10; ++k)
      tWcalfile << Rst[4][k + i * 10] << " ";
    tWcalfile << endl;
  }
  tWcalfile << endl;
  tWcalfile << "# particle = " << partType << endl;
  tWcalfile << "# EnergyBinWidth = " << EnergyBinWidth << endl;
  tWcalfile << "# RangeBinWidth = " << RangeBinWidth << endl;
  tWcalfile << "# Valid date and run ranges for this calibration:" << endl;
  tWcalfile << "# minDate = " << minDate << endl;
  tWcalfile << "# maxDate = " << maxDate << endl;
  tWcalfile << "# minRun = " << minRun << endl;
  tWcalfile << "# maxRun = " << maxRun << endl;
  tWcalfile << "# Date and time of calibration run: " << now->tm_year + 1900 << "-" << now->tm_mon + 1 << "-"
            << now->tm_mday << "    " << now->tm_hour << ":" << now->tm_min << endl;
  tWcalfile << "# The pre-processing code version number is [" << programVersion << "]" << endl;
  tWcalfile << "# File used for empty run for calibrating TV corrections was " << calFileNames[0] << endl;
  tWcalfile << "# The file to where the TV corrections were written was " << Outputdir + "/TVcorr.txt" << endl;
  tWcalfile << "# File used for calibration with wedge and 0 bricks was " << calFileNames[1] << endl;
  tWcalfile << "# File used for calibration with wedge and 1 brick  was " << calFileNames[2] << endl;
  tWcalfile << "# File used for calibration with wedge and 2 bricks was " << calFileNames[3] << endl;
  tWcalfile << "# File used for calibration with wedge and 3 bricks was " << calFileNames[4] << endl;
  tWcalfile << "# File used for calibration with wedge and 4 bricks was " << calFileNames[5] << endl;
  tWcalfile << "# Threshold0 = " << stuff.Thr[0] << endl;
  tWcalfile << "# Threshold1 = " << stuff.Thr[1] << endl;
  tWcalfile << "# Threshold2 = " << stuff.Thr[2] << endl;
  tWcalfile << "# Threshold3 = " << stuff.Thr[3] << endl;
  tWcalfile << "# Threshold4 = " << stuff.Thr[4] << endl;
  for (int stage = 1; stage < 5; stage++) {
    tWcalfile << "# dEE" + to_string(stage) + " = ";
    for (int i = 0; i < 3; i++)
      tWcalfile << dEElow[stage][i] << ",";
    for (int i = 0; i < 3; i++)
      tWcalfile << dEEhigh[stage][i] << ",";
    tWcalfile << endl;
  }

  tWcalfile.close();

  return 0;
};

void pCTcalib::plot2D(string fn, string T, string TX, string TY, int N, float X[], float Y[], float E[]) {

  FILE *oFile = fopen(fn.c_str(), "w");
  if (oFile != NULL) {
    fprintf(oFile, "#*** This file is intended to be displayed by gnuplot.\n");
    fprintf(oFile, "#*** Either double click on the file (works in Windows at least),\n");
    fprintf(oFile, "#*** or else start up gnuplot and use the load command to "
                   "display the plot.\n");
    // fprintf(oFile,"set term X11 persist\n");
    fprintf(oFile, "set title '%s' \n", T.c_str());
    fprintf(oFile, "set xlabel '%s' \n", TX.c_str());
    fprintf(oFile, "set ylabel '%s' \n", TY.c_str());
    fprintf(oFile, "set xrange[%7.4f : %7.4f]\n", X[0], X[N - 1]);
    fprintf(oFile, "set nokey\n");
    fprintf(oFile, "plot '-' with xyerrorbars\n");
    for (int i = 0; i < N; i++) {
      fprintf(oFile, "%8.3e %8.3e %8.3e\n", X[i], Y[i], E[i]);
    }
    fprintf(oFile, "e\n");
    fclose(oFile);
  }
};
