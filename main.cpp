//*******************************************************************************************************
//
// Reorganization and rewrite of the Phase-II scanner data preprocessing code,
// starting from the
// code of P. Piersimoni and including the continuous scan facility added by
// C.E. Ordonez.  Also,
// the WEPL detector calibration code of V. Bashkirov is included.
// R.P. Johnson  May 8, 2016
// R.P. Johnson  October 4, 2016  Integrated all of the TV and WEPL calibration
// code of Vladimir Bashkirov.
// The objective was to package all of the different pieces of the program into
// C++ classes, to make
// the code easier to follow and maintain. This includes
// - Preprocessing.h and Preprocessing.cpp: the top level driving program,
// called either by main here or from within another program
// - pCTgeo.h: a geometry package to encapsulate all of the geometry constants
// - Wepl1.h: WEPL calibration class, from V. Bashkirov
// - pedGainCalib.cpp and pedGainCalib.h: on-the-fly adjustment of pedestals and
// gains
// - TVcorrection.h: class for TV corrections of energy detector data
// - pCTraw.cpp and pCTraw.h: raw data input and unpacking class, encapsulates
// the original code of Piersimoni et al.
// - TkrHits.cpp and TkrHits.h: class for calculation of tracker coordinates
// from strip clusters
// - pCT_Tracking.cpp and pCT_Tracking.h: pattern recognition for the tracking;
// includes printing and plotting
// - pCTcalib executes the calibration sequence when provided with the necessary
// 6 calibration run raw data files.
// - EvtRecon does the raw data event reconstruction in the case of calibration
// runs.
//*******************************************************************************************************

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
#include <signal.h>
#include "Preprocessing.h"
#include "pCTcalib.h"
#include "pCTconfig.h"
#include "Util.h"

using namespace std;
int main(int argc, char *argv[]) {
  // Entry point for preprocessing of data from the pCT phase-II scanner.
  // R.P. Johnson 5/22/2016
  // R.P. Johnson 10/2/2016  Integrated Vladimir's TV and WEPL calibration
  // algorithms into this framework.

  float Version = 4.0;

  std::string version = std::to_string(Version);

  for (int i = 0; i < argc; i++)
    cout << argv[i] << " ";
  cout << endl;

  const string configFile = argv[1];//"pCT_config.txt";
  pCTconfig cfg(configFile); // Create a class instance for parsing the configuration file

  int numbTkrFPGA = 12;
  cfg.addItem("nTracker", numbTkrFPGA);

  int numbEdetFPGA = 2;
  cfg.addItem("nEdet", numbEdetFPGA);

  int pdstlr[5];
  pdstlr[0] = -500;
  cfg.addItem("pedrng0", pdstlr[0]);

  pdstlr[1] = -500;
  cfg.addItem("pedrng1", pdstlr[1]);

  pdstlr[2] = -500;
  cfg.addItem("pedrng2", pdstlr[2]);

  pdstlr[3] = -500;
  cfg.addItem("pedrng3", pdstlr[3]);

  pdstlr[4] = -500;
  cfg.addItem("pedrng4", pdstlr[4]);

  float fileFraction = 1.0;
  cfg.addItem("fraction", fileFraction);

  string partType = "H";
  cfg.addItem("partType", partType);

  int fileBins = 1;
  cfg.addItem("bins", fileBins);

  int continuous_scan = 1;
  cfg.addItem("continuous", continuous_scan);

  string Outputdir = ".";
  cfg.addItem("outputDir", Outputdir);
  
  int max_events = 0;
  cfg.addItem("max_events", max_events);

  int max_time = 0;
  cfg.addItem("max_time", max_time);

  int n_debug = 1;
  cfg.addItem("n_debug", n_debug);

  int n_plot = 0;
  cfg.addItem("n_plot", n_plot);

  string logFile = "";
  cfg.addItem("log", logFile);

  float initialAngle = 0.0;
  cfg.addItem("angle", initialAngle);

  float beamEnergy = 200.; // Unless a value is supplied by the user, this angle will be taken from the input file
  cfg.addItem("energy", beamEnergy);

  float proj_angle = -999.; // Unless a value is supplied by the user, this angle will be taken from the input file
  cfg.addItem("projection", proj_angle);

  int reCalibrate = 1;
  cfg.addItem("recalibrate", reCalibrate);

  float phantomSizeRight = 100.;
  cfg.addItem("sizeright", phantomSizeRight);

  float phantomSizeLeft = 100.;
  cfg.addItem("sizeleft", phantomSizeLeft);

  int Calibrate = 0;
  cfg.addItem("calibrate", Calibrate);

  int Normalize = 0;
  cfg.addItem("normalize", Normalize);
  
  float wedgeOff = 0.0;
  cfg.addItem("wedgeoffset", wedgeOff);

  string minDate = "2030/01/01";
  cfg.addItem("minDate", minDate);

  string maxDate = "2000/01/01";
  cfg.addItem("maxDate", maxDate);

  int minRun = 999;
  cfg.addItem("minrun", minRun);
  int maxRun = -1;
  cfg.addItem("maxrun", maxRun);

  string study_name = "";
  cfg.addItem("study", study_name);

  string rootCalibFile = "pCTcalib.root";
  cfg.addItem("calib", rootCalibFile);

  int CalibCurve = 0; 
  cfg.addItem("CalibCurve", CalibCurve);

  float thr[5]; // Array of stage thresholds for WEPL analysis
  thr[0] = 1.0;
  cfg.addItem("thr0", thr[0]);

  thr[1] = 1.0;
  cfg.addItem("thr1", thr[1]);

  thr[2] = 1.0;
  cfg.addItem("thr2", thr[2]);

  thr[3] = 1.0;
  cfg.addItem("thr3", thr[3]);

  thr[4] = 1.0;
  cfg.addItem("thr4", thr[4]);

  int dodEEFilter = 1; // changed default to yes
  cfg.addItem("dEEFilter", dodEEFilter); // Also add the option to the list used for parsing the config file
 
  int maxFluence = 1000; //max fluence per mm^2 (full Nb of part = maxFluence x fieldSize)
  cfg.addItem("maxFluence", maxFluence);
 
  int MultiTrackReject = 0; 
  cfg.addItem("MultiTrackReject",MultiTrackReject);

  int CTOutput = 0;
  cfg.addItem("CTOutput",CTOutput);


  float TpinOff1 = -1.; // Offset of the Tpins 
  cfg.addItem("TpinOff1", TpinOff1); 
  float TpinOff2 = -1.; // Offset of the Tpins 
  cfg.addItem("TpinOff2", TpinOff2);
  float TpinOff3 = 1.218; // Offset of the Tpins 
  cfg.addItem("TpinOff3", TpinOff3);
  float TpinOff4 = 1.363; // Offset of the Tpins //1.218 + 0.145
  cfg.addItem("TpinOff4", TpinOff4);

  float VpinOff1 = 0; // Offset of the Tpins 
  cfg.addItem("VpinOff1", VpinOff1);
  float VpinOff2 = 0; // Offset of the Tpins 
  cfg.addItem("VpinOff2", VpinOff2);
  float VpinOff3 = 1.218; // Offset of the Tpins 
  cfg.addItem("VpinOff3", VpinOff3);
  float VpinOff4 = 1.363; // Offset of the Tpins //1.218 + 0.145
  cfg.addItem("VpinOff4", VpinOff4); 


  // Read the default configuration from the config file
  if (cfg.Configure() != 0) {
    cout << "Was not able to read a default configuration from " << configFile << endl;
    cout << "The hardwired default configuration will be used." << endl;
  }

  //////////////////////////////////////////////////////
  // Printing out some number
  //////////////////////////////////////////////////////
  if (numbTkrFPGA  != 12) cout << "Non-standard number of tracker FPGAs in the readout = " << numbTkrFPGA << endl;
  if (numbEdetFPGA != 2)  cout << "Non-standard number of energy detector FPGAs in the readout = " << numbEdetFPGA << endl;
    
  // Weird validation
  if (cfg.item_int["recalibrate"]) cout << "Energy detector stage gains will be recalibrated on the fly during processing." << endl;
  else cout << "Energy detector stage gains will NOT be recalibrated on the fly during processing." << endl;
  if (fileBins <= 0) {
    cout << "************ The number of files was specified to be 0 or negative. Resetting to equal 1 bin. **********" << endl;
    fileBins = 1; // Protects from crashing due to bad input
  }
  if (max_events > 0) {
    cout << "The maximum number of events to analyze is set to " << max_events << endl;
  } else cout << "No restriction is set on the maximum number of events to analyze." << endl;
  if (max_time > 0) {
    cout << "The maximum time stamp to analyze is set to " << max_time << " seconds." << endl;
  } else cout << "No restriction is set on the maximum time stamp to analyze." << endl;

  // Most of the floating point variables are specified double precision on the assumption
  // that this is going to execute anyway on a 64-bit machine.  An exception is the
  // temporary data file and the binary output, which are intended to use 4-byte floating point.

  cout << "Executing " << argv[0] << " version " << version << endl;
  cout << "  float is " << sizeof(float) << " bytes\n";
  cout << "  double is " << sizeof(double) << " bytes\n";
  cout << "  char is " << sizeof(char) << " bytes\n";
  cout << "  int is " << sizeof(int) << " bytes\n";
  cout << "  long long is " << sizeof(long long) << " bytes\n";
  cout << "  long is " << sizeof(long) << " bytes\n";

  // The user has to enter the full filename including path
  string CalFile = "CalFileList.txt"; // For calibration runs the input filenames are taken from here
  string inputFileName;

  // Get the list of required, position-sensitive arguments, in this case just the input filename
  //vector<string> requiredArgs = //parser.args();
  if (argc == 0) {
    if (!cfg.item_int["calibrate"]) {
      cout << "pCT_Preprocessing: no input raw data file was specified!\n";
      exit(1);
    }
  } else {
    inputFileName = argv[2];
    CalFile = argv[2];
  }
  if (cfg.item_int["calibrate"]) cout << "Calibration run.  The list of input files is from " << CalFile << endl;
  else cout << "Preprocessing run, the input file name is " << inputFileName << endl;
  cout << "The TV calibration file is " << cfg.item_str["TVcorr"] << endl;
  cout << "The WEPL calibration file is " << cfg.item_str["Wcalib"] << endl;

  if (continuous_scan) {
    cout << "The data are assumed to be from a continuous scan." << endl;
    if (initialAngle <= 0.0) {
      if (logFile == "" || logFile == "NULL" || logFile == "null" || logFile == "Null") {
        size_t found = inputFileName.find_last_of(".");
        if (found != inputFileName.npos) {
          logFile = inputFileName.substr(0, found) + ".log";
        }
      }
      FILE *ftst = fopen(logFile.c_str(), "r"); // Just to test whether the log file really exists. . .
      if (ftst != NULL) {
        fclose(ftst);
        Util util;
        initialAngle = util.getStartAngle(logFile);
        cout << "The initial stage angle " << initialAngle << " was calculated from the log file " << logFile << endl;
      } else {
        initialAngle = 0.0;
        cout << "Could not open the log file.  Setting the initial stage angle to zero." << endl;
      }
    } else {
      cout << "The initial stage angle was set by the user to " << initialAngle << " degrees." << endl;
    }
  } else
    cout << "The data are assumed to be from a single projection of a stepped scan." << endl;

  cout << "The number of events for debug printing is " << n_debug << endl;
  cout << "The number of events for which to plot the tracker hits and tracks is " << n_plot << endl;

  if (cfg.item_int["calibrate"]) cout << "Set the number of events to plot > 0 to get loads of debug histograms in calibration runs." << endl;
  cout << "Fraction of the input file to be analyzed is " << fileFraction << endl;
  cout << "The phantom size for preprocessing is assumed to be " << phantomSizeLeft << " mm in extent in -T and " << phantomSizeRight << " mm in +T." << endl;
    
  if (dodEEFilter) {
    if (!cfg.item_int["calibrate"])
      cout << "The dE-E filtering of nuclear interactions will be used before WEPL reconstruction" << endl;
    else if (!cfg.item_int["calibrate"] && cfg.item_str["partType"] == "He")
      cout << "WARNING: helium fragments will be included in the analysis!" << endl;
  } else
    cout << "No dE-E filtering of nuclear interactions will be used" << endl;
  ////////////////////////////////////////////////////
  // Calibration run
  ////////////////////////////////////////////////////
  if (cfg.item_int["calibrate"]) { // Calibration Run
    cout << "Running the pCT TV and WEPL calibration sequence" << endl;
    if (Normalize) cout << "Will normalize columns in the WET vs E plot to have equal area." << endl;

    // Few more options to be used later in the calibration
    int Nbricks = 1; 
    cfg.addItem("Nbricks", Nbricks);
    int doGains = 1; 
    cfg.addItem("doGains", doGains);
    
    pCTcalib calibProcessor(CalFile);
    if (calibProcessor.TVmapper() == 0) { // First the TVmapper
    //if (calibProcessor.TVmapper_FlatBricks() == 0) { // First the TVmapper
      calibProcessor.enrgDep(); // Verify the energy dependence      
      calibProcessor.Wcalib();
      calibProcessor.writeCalibfile();
    }
    else {
      cout << "The calibration run failed in calibProcessor.TVmapper; WEPL calibration will not be run." << endl;
    }
  }

  ////////////////////////////////////////////////////
  // Real run
  ////////////////////////////////////////////////////
  else {
    cfg.addItem("inputFileName", inputFileName);
    cout << "Executing a pCT data pre-processing run" << endl;
    // Here we call the complete preprocessing program
    Preprocessing pCTpreprocessor;
    int errorCode = pCTpreprocessor.ProcessFile(fileFraction, numbTkrFPGA, numbEdetFPGA);
    return errorCode;
  }
} // end of the main program
