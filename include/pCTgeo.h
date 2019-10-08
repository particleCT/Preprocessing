#ifndef pCTgeo_h
#define pCTgeo_h
// Class to contain all of the pCT scanner geometry information, plus other
// constants, needed for preprocessing
// R.P. Johnson 5/22/2016

using namespace std;

#include <vector>
#include <iostream>
#include <fstream>

class pCTgeo {
  int fpgaLyr[12]; // mapping from fpga number to layer

  double uPosV[4]; // u locations of each of the V layers
  int VBoard[4]; // hardware identifier for each V layer
  double Vpin[4]; // cassette alignment in V layers
  double firstStripV[7][2]; // offset of the first strip on each sensor of each
                            // hardware V board
  // relative to the alignment pin, as measured under a microscope by Forest

  double uPosT[4]; // u locations of each of the T layers
  double Tpin[4]; // cassette alignment in T layers
  double Tdir[4]; // T-board orientation
  int TBoard[4]; // hardware identifier for each T layer
  double firstStripT[7][4]; // offset of the first strip on each sensor of each
                            // hardware T board
  // relative to the alignment pin, as measured under a microscope by Forest

  double tGap[3][4]; // gap locations for each T layer

  double stripPitch;
  std::vector<double> BeamVtxV; // approximate origin of the proton beam (the
                                // lead foil in LLUMC runs)
  std::vector<double> BeamVtxT; // at CPC it is different in the V and T views
  double uEdet; // approximate u coordinate of the front face of the energy
                // detector

  double TwedgeBreaks[4]; // Break points in the wedge phantom geometry
  double BrickThickness;  // Brick thickness; also the wedge thickness

  double StageThickness; // Thickness in u of each energy detector stage

  double rotationSpeed; // Speed of the stage rotation in continuous scans, in
                        // degrees/s
  double timeStampRes;  // Hardware time-stamp resolution in seconds

public:
  pCTgeo(double tWedgeOffset = 0.) { // Most of this stuff needs to be kept in a database and read from there in this constructor
    ofstream geoLogFile;
    geoLogFile.open ("pCTGeometry.log");
    geoLogFile << "pCTgeo.h: loading geometry constants for the real Phase-II pCT Scanner." << endl;
    fpgaLyr[0] = 0; // Translation from FPGA to layer
    fpgaLyr[1] = 1;
    fpgaLyr[2] = 2;
    fpgaLyr[3] = 3;
    fpgaLyr[4] = 0;
    fpgaLyr[5] = 0;
    fpgaLyr[6] = 1;
    fpgaLyr[7] = 1;
    fpgaLyr[8] = 2;
    fpgaLyr[9] = 2;
    fpgaLyr[10] = 3;
    fpgaLyr[11] = 3;
    for (int i = 0; i < 12; i++)
      geoLogFile << "pCTgeo.h: FPGA number " << i << " is located in layer "
           << fpgaLyr[i] << endl;

    rotationSpeed = 6.0; // The stage velocity is 1 rpm or 6 deg/s
    timeStampRes = 1.0e-8;
    geoLogFile << "pCTgeo.h: stage rotation speed is assumed to equal "
         << rotationSpeed << " degrees/s" << endl;
    geoLogFile << "pCTgeo.h: time-stamp resolution is assumed to equal "
         << timeStampRes << " seconds" << endl;
    double beamZV = -3000; //-1850.;   // Average value for CPC, from Mark
                           //Pankuch, Sept 2018
    double beamZT = -3000; //-1850.; // NEW MEASURED DISTANCE GIVES WORSE
                           //RESULTS FOR WHATEVER REASON
    double beamX = 0.;
    double beamY = 0.;
    BeamVtxV.push_back(beamX);
    BeamVtxV.push_back(beamY);
    BeamVtxV.push_back(beamZV);
    BeamVtxT.push_back(beamX);
    BeamVtxT.push_back(beamY);
    BeamVtxT.push_back(beamZT);
    geoLogFile << "pCTgeo.h: the V-view beam origin is assumed to be at x=" << beamX
         << " mm, y=" << beamY << " mm, z=" << beamZV << "mm\n";
    geoLogFile << "pCTgeo.h: the T-view beam origin is assumed to be at x=" << beamX
         << " mm, y=" << beamY << " mm, z=" << beamZT << "mm\n";

    uEdet =
        216.9 + 40.0; // Approximate location of the energy detector entrance
    geoLogFile << "pCTgeo.h: the entrance to the energy detector is assumed to be "
            "located at u=" << uEdet << endl;
    
    uPosV[0] =-217.3; // U positions of the 8 tracker boards in the Phase-II scanner
    uPosV[1] = -167.2;
    uPosV[2] = 167.2;
    uPosV[3] = 217.3;
    uPosT[0] = -211.4;
    uPosT[1] = -161.4;
    uPosT[2] = 161.4;
    uPosT[3] = 211.4;

    stripPitch = 0.228;
    geoLogFile << "pCTgeo.h: Silicon-strip detector strip pitch = " << stripPitch
         << endl;
    geoLogFile << "pCTgeo.h: u positions of the V and T tracker boards for tracker "
            "layer 0 are " << uPosV[0] << " " << uPosT[0] << endl;
    geoLogFile << "pCTgeo.h: u positions of the V and T tracker boards for tracker "
            "layer 1 are " << uPosV[1] << " " << uPosT[1] << endl;
    geoLogFile << "pCTgeo.h: u positions of the V and T tracker boards for tracker "
            "layer 2 are " << uPosV[2] << " " << uPosT[2] << endl;
    geoLogFile << "pCTgeo.h: u positions of the V and T tracker boards for tracker "
            "layer 3 are " << uPosV[3] << " " << uPosT[3] << endl;

    VBoard[0] = 6; // Hardware configuration of V boards
    VBoard[1] = 4;
    VBoard[2] = 2;
    VBoard[3] = 1;

    Vpin[0] = 0.001; // V board alignment pin locations, including offsets derived from data analysis
    Vpin[1] = 0.115;
    Vpin[2] = 0.063;
    Vpin[3] = 0.040;

    for (int i = 0; i < 4; i++)
      geoLogFile << "pCTgeo.h: V layer " << i << " is board " << VBoard[i]
           << " with alignment pin at v=" << Vpin[i] << " mm\n";

    double fSV[7][2] = // V board internal alignment from optical surveys of
                       // physics boards
        { { -43.7193, -43.716 },
          { -43.727, -43.687 },
          { -43.682, -43.681 },
          { -43.702, -43.713 },
          { -43.7193, -43.716 },
          { -43.7193, -43.716 },
          { -43.686, -43.687 } };
    double fSV_MC[2] = { -43.717, -43.717 }; // V board alignment assumed in the
                                             // simulation

    Tpin[0] = 215.168; // T board alignment pin locations, including corrections derived from data analysis
    Tpin[1] = 211.373;
    Tpin[2] = -203.373;
    Tpin[3] = -207.168;

    Tdir[0] = -1.; // T board orientations (front and back trackers are
                   // reflected in u)
    Tdir[1] = -1.;
    Tdir[2] = 1.;
    Tdir[3] = 1;

    TBoard[0] = 5; // Hardware configuration of T boards
    TBoard[1] = 4;
    TBoard[2] = 1;
    TBoard[3] = 2;
    for (int i = 0; i < 4; i++)
      geoLogFile << "pCTgeo.h: T layer " << i << " is board " << TBoard[i]
           << " with alignment pin at t=" << Tpin[i] << " mm and direction "
           << Tdir[i] << endl;

    double fST[7][4] = { // T board internal alignment from optical surveys of
                         // the physical boards
      { -999, -999, -999, -999 }, { 38.60, 126.87, 215.15, 303.42 },
      { 38.48, 126.76, 215.04, 303.32 }, { 38.69, 126.95, 215.23, 303.57 },
      { 38.58, 126.85, 215.11, 303.37 }, { 38.62, 126.90, 215.16, 303.41 },
      { 38.58, 126.85, 215.11, 303.37 } };
    double fST_MC[4][4] = { // T board alignment assumed in the simulation, from
                            // Hartmut's ppt file
      { 38.556, 126.825, 215.111, 303.382 },
      { 38.63, 126.91, 215.187, 303.474 }, { 38.835, 127.104, 215.38, 303.719 },
      { 38.542, 126.812, 215.071, 303.331 }, };


      for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 4; j++) {
          if (j < 2)
            firstStripV[i][j] = fSV[i][j];
          firstStripT[i][j] = fST[i][j];
        }
      }

    for (int i = 0; i < 4; i++) {
      geoLogFile << "pCTgeo.h: V board layer " << i
           << " first strip locations per ladder:";
      for (int j = 0; j < 2; j++)
        geoLogFile << " " << firstStripV[VBoard[i]][j];
      geoLogFile << endl;
    }
    for (int i = 0; i < 4; i++) {
      geoLogFile << "pCTgeo.h: T board layer " << i
           << " first strip locations per wafer:";
      for (int j = 0; j < 4; j++)
        geoLogFile << " " << firstStripT[TBoard[i]][j];
      geoLogFile << endl;
    }

    // Find the T-layer gap locations, half way between the surrounding strips
    for (int lyr = 0; lyr < 4; lyr++) {
      int brd = TBoard[lyr];
      geoLogFile << "pCTgeo.h: t-layer gap locations for layer " << lyr << ":";
      for (int gap = 0; gap < 3; gap++) {
        double t1 = Tpin[lyr] + Tdir[lyr] * (firstStripT[brd][gap] +
                                             (6.0 * 64.0 - 1.0) * stripPitch);
        double t2 = Tpin[lyr] + Tdir[lyr] * firstStripT[brd][gap + 1];
        // geoLogFile << " t1=" << t1 << " t2=" << t2;
        tGap[gap][lyr] = 0.5 * (t1 + t2);
        geoLogFile << "  " << tGap[gap][lyr];
      }
      geoLogFile << std::endl;
    }

    // Geometry for the wedge phantom, all in mm (INCLUDING 0.1MM SHIFT AND
    // 0.1MM SLIT between wedges):
    TwedgeBreaks[0] = -104.50 + tWedgeOffset; // Start of the wedge slope
    TwedgeBreaks[1] =
        -4.75 + tWedgeOffset; // End of the slope, start of the flat
    TwedgeBreaks[2] =
        4.75 + tWedgeOffset; // End of the slope, start of the opposite slope
    TwedgeBreaks[3] = 104.50 + tWedgeOffset; // End of the opposite slope.
                                             // Should also be the end of the
                                             // bricks, if positioned correctly
    BrickThickness =
        50.8; // Brick thickness in U; also the wedge maximum thickness
    geoLogFile << "pCTgeo.h: calibration brick and wedge thickness = "
         << BrickThickness << endl;
    geoLogFile << "pCTgeo.h: calibration phantom wedge break locations in t are ";
    for (int i = 0; i < 4; i++)
      geoLogFile << TwedgeBreaks[i] << ", ";
    geoLogFile << endl;

    StageThickness = 51.;
    geoLogFile << "pCTgeo.h: thickness of each energy detector stage = "
         << StageThickness << " mm" << endl;
    geoLogFile.close();
  };

  inline double stageSpeed() { return rotationSpeed; }

  inline double timeRes() { return timeStampRes; }

  inline int getFPGAlayer(int FPGA) { return fpgaLyr[FPGA]; }

  inline double getTWedgeBreaks(int n) const { // Return the 4 calibration wedge
                                               // break points in t; n=1,2,3,4
                                               // in increasing t
    return TwedgeBreaks[n - 1];
  }

  inline double getBrickThickness() const { // Return the calibration brick and
                                            // wedge thickness
    return BrickThickness;
  }

  inline double energyDetectorU(int stage) const { // stage from 0 through 4
    return uEdet + stage * (StageThickness);
  }

  inline double tBoardGap(int gap, int lyr) const { // Method to return the
                                                    // t-layer gap locations
    return tGap[gap][lyr];
  }

  inline int
  Lyr(int FPGA) const { // Method to translate from FPGA to layer number
    return fpgaLyr[FPGA];
  };

  inline double
  uV(int Lyr) const { // Method to return the U coordinate of V boards
    return uPosV[Lyr];
  }

  inline double
  uT(int Lyr) const { // Method to return the U coordinate of T boards
    return uPosT[Lyr];
  }

  inline std::vector<double> BeamVertex(int view) const {
    if (view == 0)
      return BeamVtxV;
    else
      return BeamVtxT;
  }

  inline double extrap2D(double X[2], double Y[2], double Xnew) {
    double dX = X[1] - X[0];
    if (dX == 0.) {
      cout << "pCTgeo::extrap2D, division by zero; x values must be different."
           << endl;
      dX = 1.0e-23;
    }
    double slope = (Y[1] - Y[0]) / (dX);
    return Y[1] + (Xnew - X[1]) * slope;
  }

  // Method to translate strip number to coordinate for V layers
  inline double strip_position_V(int fpga, int chip, double strip) const {
    double global_strip;
    int side;
    if (chip < 6) {
      global_strip = (5 - chip) * 64.0 + strip;
      side = 0;
    } else {
      global_strip = (chip - 6) * 64.0 + (63.0 - strip);
      side = 1;
    }
    int brd = VBoard[fpga];
    // std::cout << "strip_position_V, fpga " << fpga << " chip " << chip << "
    // strip " << strip << std::endl;
    // std::cout << "   brd=" << brd << " 1st strip=" << firstStripV[brd][side]
    // << " glob=" << global_strip << " Vpin=" << Vpin[fpga] << std::endl;
    return firstStripV[brd][side] + Vpin[fpga] + global_strip * stripPitch;
  }

  // Method to translate strip number to coordinate for T layers
  inline double strip_position_T(int fpga, int chip, double strip) const {
    int sensor = 2 * (fpga % 2) + chip / 6;
    int chip_in_sensor = chip % 6;
    double global_strip = 64.0 * chip_in_sensor + (63.0 - strip);
    int lyr = fpgaLyr[fpga];
    int brd = TBoard[lyr];
    // std::cout << "strip_position_T, fpga " << fpga << " chip " << chip << "
    // strip " << strip << " layer " << lyr << std::endl;
    // std::cout << "   brd=" << brd << " 1st strip=" <<
    // firstStripT[brd][sensor] << " glob=" << global_strip << " Tpin=" <<
    // Tpin[lyr] << std::endl;
    return Tpin[lyr] +
           Tdir[lyr] * (firstStripT[brd][sensor] + global_strip * stripPitch);

  }

};
#endif
