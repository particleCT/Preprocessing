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
  cfg.addItem('q', "nTracker", numbTkrFPGA);

  int numbEdetFPGA = 2;
  cfg.addItem('Q', "nEdet", numbEdetFPGA);

  int pdstlr[5];
  pdstlr[0] = -500;
  cfg.addItem('5', "pedrng0", pdstlr[0]);

  pdstlr[1] = -500;
  cfg.addItem('6', "pedrng1", pdstlr[1]);

  pdstlr[2] = -500;
  cfg.addItem('7', "pedrng2", pdstlr[2]);

  pdstlr[3] = -500;
  cfg.addItem('8', "pedrng3", pdstlr[3]);

  pdstlr[4] = -500;
  cfg.addItem('9', "pedrng4", pdstlr[4]);

  float fileFraction = 1.0;
  cfg.addItem('f', "fraction", fileFraction);

  string partType = "H";
  cfg.addItem('S', "partType", partType);

  int fileBins = 1;
  cfg.addItem('b', "bins", fileBins);

  int continuous_scan = 1;
  cfg.addItem('c', "continuous", continuous_scan);

  string Outputdir = ".";
  cfg.addItem('o', "outputDir", Outputdir);
  
  int max_events = 0;
  cfg.addItem('n', "max_events", max_events);

  int max_time = 0;
  cfg.addItem('M', "max_time", max_time);

  int n_debug = 1;
  cfg.addItem('d', "n_debug", n_debug);

  int n_plot = 0;
  cfg.addItem('j', "n_plot", n_plot);

  string logFile = "";
  cfg.addItem('g', "log", logFile);

  float initialAngle = 0.0;
  cfg.addItem('a', "angle", initialAngle);

  float beamEnergy = 200.; // Unless a value is supplied by the user, this angle will be taken from the input file
  cfg.addItem('e', "energy", beamEnergy);

  float proj_angle = -999.; // Unless a value is supplied by the user, this angle will be taken from the input file
  cfg.addItem('p', "projection", proj_angle);

  int reCalibrate = 1;
  cfg.addItem('r', "recalibrate", reCalibrate);

  float phantomSize = 110.;
  cfg.addItem('S', "size", phantomSize);

  int Calibrate = 0;
  cfg.addItem('C', "calibrate", Calibrate);

  int Normalize = 0;
  cfg.addItem('L', "normalize", Normalize);
  
  float wedgeOff = 0.0;
  cfg.addItem('F', "wedgeoffset", wedgeOff);

  string minDate = "2030/01/01";
  cfg.addItem('x', "minDate", minDate);

  string maxDate = "2000/01/01";
  cfg.addItem('y', "maxDate", maxDate);

  int minRun = 999;
  cfg.addItem('w', "minrun", minRun);

  int maxRun = -1;
  cfg.addItem('z', "maxrun", maxRun);

  string study_name = "";
  cfg.addItem('s', "study", study_name);

  string WcalibFile = "Wcalib.txt";
  cfg.addItem('W', "Wcalib", WcalibFile);

  string TVcorrFile = "TVcorr.txt";
  cfg.addItem('T', "TVcorr", TVcorrFile);

  string rootCalibFile = "pCTcalib.root";
  cfg.addItem('R', "calib", rootCalibFile);
  
  float thr[5]; // Array of stage thresholds for WEPL analysis
  thr[0] = 1.0;
  cfg.addItem('0', "thr0", thr[0]);

  thr[1] = 1.0;
  cfg.addItem('1', "thr1", thr[1]);

  thr[2] = 1.0;
  cfg.addItem('2', "thr2", thr[2]);

  thr[3] = 1.0;
  cfg.addItem('3', "thr3", thr[3]);

  thr[4] = 1.0;
  cfg.addItem('4', "thr4", thr[4]);
  cout<<thr[2]<<endl;
  int dodEEFilter = 1; // changed default to yes
  cfg.addItem('e', "dEEFilter", dodEEFilter); // Also add the option to the list used for parsing the config file

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
  cout << "The phantom size for preprocessing is assumed to be " << phantomSize << " mm in radius." << endl;
    
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
    cfg.addItem('b', "Nbricks", Nbricks);
    int doGains = 1; 
    cfg.addItem('g', "doGains", doGains);
    
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
    cfg.addItem('i', "inputFileName", inputFileName);
    cout << "Executing a pCT data pre-processing run" << endl;
    // Here we call the complete preprocessing program
    Preprocessing pCTpreprocessor;
    int errorCode = pCTpreprocessor.ProcessFile(fileFraction, numbTkrFPGA, numbEdetFPGA);
    return errorCode;
  }
} // end of the main program
