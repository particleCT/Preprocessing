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
//
// The objective was to package all of the different pieces of the program into
// C++ classes, to make
// the code easier to follow and maintain. This includes
// - arg.h and arg.cpp: quasi-unix-style command-line parsing (Copyright (C)
// 2010 Chun-Chung Chen <cjj@u.washington.edu>)
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
// - Histogram.cpp and Histogram.h: simple facilities for making histograms and
// profile plots
//
// Compile it with g++ using -std=c++0x -lpthread -g -rdynamic (the -g and
// -rdynamic can be omitted if all is working well)
// For Linux, see the provided Makefile.
//
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
#include "arg.h" // Command line options
#include "Preprocessing.h"
#include "pCTcalib.h"
#include "pCTconfig.h"
#include "Util.h"

using namespace std;

// |\\  //|    //\    ||  |\\  || 
// ||\\//||   // \\   ||  ||\\ || 
// || \/ ||  //ZZZ\\  ||  || \\|| 
// ||    || //     \\ ||  ||  \\| 

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

  const string configFile = "pCT_config.txt";
  pCTconfig cfg(configFile); // Create a class instance for parsing the configuration file
  arg::Parser parser; // Create a class instance for parsing the command line (see the program arg.cc and header arg.h)

  parser.set_header(" \n ******** pCT preprocessing version " + version + " **********");
  parser.add_help("");
  parser.add_help(" There is one positional argument: file_name");
  parser.add_help("   The file name is either the raw data file or, for "
                  "calibration, the text file with a list of ");
  parser.add_help("   calibration raw data files (default 'CalFileList.txt').");
  parser.add_help(" Default values of all parameters may be set with a "
                  "pCT_config.txt file.");
  parser.add_help(" Command line options will override the pCT_config.txt settings.");
  parser.add_help("   The command line option syntax can take either a short "
                  "or long form, e.g.:");
  parser.add_help("          -n 3000000");
  parser.add_help("               or equivalently");
  parser.add_help("          --number=3000000");
  parser.add_help("   The parser for pCT_config.txt can recognize the short "
                  "and long names, but the syntax is always as follows:");
  parser.add_help("          n = 3000000");
  parser.add_help("               or equivalently");
  parser.add_help("          number = 3000000");
  parser.add_help("   Use the string NULL to represent a null or empty string.");
  parser.add_help(" Available options are:");

  // Define all of the command line options
  string KillCh = "null";
  parser.add_opt('k', "kill")
      .stow(KillCh)
      .help("Filename for the list of tracker channels to kill", "STRING")
      .show_default();
  cfg.addItem('k', "kill", KillCh);

  int numbTkrFPGA = 12;
  parser.add_opt('q', "nTracker")
      .stow(numbTkrFPGA)
      .help("number of tracker FPGAs in the readout", "INT")
      .show_default();
  cfg.addItem('q', "nTracker", numbTkrFPGA);

  int numbEdetFPGA = 2;
  parser.add_opt('Q', "nEdet")
      .stow(numbEdetFPGA)
      .help("number of energy detector FPGAs in the readout", "INT")
      .show_default();
  cfg.addItem('Q', "nEdet", numbEdetFPGA);

  int pdstlr[5];
  pdstlr[0] = -500;
  parser.add_opt('5', "pedrng0")
      .stow(pdstlr[0])
      .help("stage 0 start of pedestal region (ADC counts)", "INT")
      .show_default();
  cfg.addItem('5', "pedrng0", pdstlr[0]);

  pdstlr[1] = -500;
  parser.add_opt('6', "pedrng1")
      .stow(pdstlr[1])
      .help("stage 1 start of pedestal region (ADC counts)", "INT")
      .show_default();
  cfg.addItem('6', "pedrng1", pdstlr[1]);

  pdstlr[2] = -500;
  parser.add_opt('7', "pedrng2")
      .stow(pdstlr[2])
      .help("stage 2 start of pedestal region (ADC counts)", "INT")
      .show_default();
  cfg.addItem('7', "pedrng2", pdstlr[2]);

  pdstlr[3] = -500;
  parser.add_opt('8', "pedrng3")
      .stow(pdstlr[3])
      .help("stage 3 start of pedestal region (ADC counts)", "INT")
      .show_default();
  cfg.addItem('8', "pedrng3", pdstlr[3]);

  pdstlr[4] = -500;
  parser.add_opt('9', "pedrng4")
      .stow(pdstlr[4])
      .help("stage 4 start of pedestal region (ADC counts)", "INT")
      .show_default();
  cfg.addItem('9', "pedrng4", pdstlr[4]);

  float fileFraction = 1.0;
  parser.add_opt('f', "fraction").stow(fileFraction).help("Fraction of the input file to use", "FLOAT").show_default();
  cfg.addItem('f', "fraction", fileFraction);

  string partType = "H";
  parser.add_opt('S', "particle")
      .stow(partType)
      .help("Particle type, hydrogen (H) or helium (He)", "STRING")
      .show_default();
  cfg.addItem('S', "partType", partType);

  int fileBins = 1;
  parser.add_opt('b', "bins")
      .stow(fileBins)
      .help("Number of files to sub-divise the results into", "INT")
      .show_default();
  cfg.addItem('b', "bins", fileBins);

  int continuous_scan = 1;
  parser.add_opt('b', "bins")
    .stow(continuous_scan)
    .help("for a continuous scan, set to 1", "INT")
    .show_default();
  cfg.addItem('c', "continuous", continuous_scan);

  string Outputdir = ".";
  parser.add_opt('o', "outputDir")
    .stow(Outputdir)
    .help("set the output directory", "STRING")
    .show_default();
  cfg.addItem('o', "outputDir", Outputdir);
  
  int max_events = 0;
  parser.add_opt('n', "max_events")
      .stow(max_events)
      .help("set maximum number of events to read and analyze", "INT")
      .show_default();
  cfg.addItem('n', "max_events", max_events);

  int max_time = 0;
  parser.add_opt('M', "time")
      .stow(max_time)
      .help("set the maximum time stamp to read and analyze", "INT")
      .show_default();
  cfg.addItem('M', "max_time", max_time);

  int n_debug = 1;
  parser.add_opt('d', "debug")
    .stow(n_debug)
    .help("set the number of events for debug printout", "INT")
    .show_default();
  cfg.addItem('d', "n_debug", n_debug);

  int n_plot = 0;
  parser.add_opt('j', "plot")
    .stow(n_plot)
    .help("set the number of tracker events to plot", "INT")
    .show_default();
  cfg.addItem('j', "n_plot", n_plot);

  int analysisLevel = 2;
  parser.add_opt('l', "level")
    .stow(analysisLevel)
    .help("0= monitor raw data only, 1= monitor raw and WEPL, 2= output ","INT")
    .show_default();
  cfg.addItem('l', "level", analysisLevel);

  string logFile = "";
  parser.add_opt('g', "log")
      .stow(logFile)
      .help("path to the log file for the run; used for continuous scans to "
            "find the start angle",
            "STRING")
      .show_default();
  cfg.addItem('g', "log", logFile);

  float initialAngle = 0.0;
  parser.add_opt('a', "angle")
      .stow(initialAngle)
      .help("initial angle, at time 0, for a continuous scan (better to provide the log file)", "FLOAT")
      .show_default();
  cfg.addItem('a', "angle", initialAngle);

  float beamEnergy = 200.; // Unless a value is supplied by the user, this angle will be taken from the input file
  parser.add_opt('p', "energy")
      .stow(beamEnergy)
      .help("set the energy (Mev/u) to override the value from the input file", "FLOAT")
      .show_default();
  cfg.addItem('e', "energy", beamEnergy);

  float proj_angle = -999.; // Unless a value is supplied by the user, this angle will be taken from the input file
  parser.add_opt('p', "projection")
      .stow(proj_angle)
      .help("set the projection angle to override the value from the input file", "FLOAT")
      .show_default();
  cfg.addItem('p', "projection", proj_angle);

  int reCalibrate = 1;
  parser.add_opt('r', "recalibrate")
      .stow(reCalibrate)
      .help("Execute real-time energy detector gain corrections? yes or no", "INT")
      .show_default();
  cfg.addItem('r', "recalibrate", reCalibrate);

  float phantomSize = 110.;
  parser.add_opt('Z', "size")
      .stow(phantomSize)
      .help("Maximum radius of the phantom in mm, for real-time gain "
            "recalibration",
            "FLOAT")
      .show_default();
  cfg.addItem('S', "size", phantomSize);

  int Calibrate = 0;
  parser.add_opt('C', "calibrate")
      .stow(Calibrate)
      .help("Produce the TV and WEPL calibration constants from calibration "
            "data? yes or no",
            "INT")
      .show_default();
  cfg.addItem('C', "calibrate", Calibrate);

  int Normalize = 0;
  parser.add_opt('L', "normalize")
      .stow(Normalize)
      .help("Normalize the colums of the WET vs E plot all to be equal area? yes or no", "INT")
      .show_default();
  cfg.addItem('L', "normalize", Normalize);
  
  float wedgeOff = 0.0;
  parser.add_opt('F', "offset")
      .stow(wedgeOff)
      .help("Offset of the calibration wedge phantom from center (normally "
            "zero)",
            "FLOAT")
      .show_default();
  cfg.addItem('F', "wedgeoffset", wedgeOff);

  string minDate = "2030/01/01";
  parser.add_opt('x', "mindate")
      .stow(minDate)
      .help("Minimum valid date for the calibration, for calibration runs, to "
            "write in the output file",
            "STRING")
      .show_default();
  cfg.addItem('x', "mindate", minDate);

  string maxDate = "2000/01/01";
  parser.add_opt('y', "maxdate")
      .stow(maxDate)
      .help("Maximum valid date for the calibration, for calibration runs, to "
            "write in the output file",
            "STRING")
      .show_default();
  cfg.addItem('y', "maxdate", maxDate);

  int minRun = 999;
  parser.add_opt('w', "minrun")
      .stow(minRun)
      .help("Minimum valid run number for the calibration, for calibration "
            "runs, to write in the output file",
            "INT")
      .show_default();
  cfg.addItem('w', "minrun", minRun);

  int maxRun = -1;
  parser.add_opt('z', "maxrun")
      .stow(maxRun)
      .help("Maximum valid run number for the calibration, for calibration "
            "runs, to write in the output file",
            "INT")
      .show_default();
  cfg.addItem('z', "maxrun", maxRun);

  string study_name = "";
  parser.add_opt('s', "study")
      .stow(study_name)
      .help("set the study name to override the one derived from the filename", "STRING")
      .show_default();
  cfg.addItem('s', "study", study_name);

  string WcalibFile = "Wcalib.txt";
  parser.add_opt('W', "Wcalib")
      .stow(WcalibFile)
      .help("path and name of WEPL calibration file", "STRING")
      .show_default();
  cfg.addItem('W', "Wcalib", WcalibFile);

  string TVcorrFile = "TVcorr.txt";
  parser.add_opt('T', "TVcorr")
      .stow(TVcorrFile)
      .help("path and name of TV gain-correction calibration file", "STRING")
      .show_default();
  cfg.addItem('T', "TVcorr", TVcorrFile);

  float thr[5]; // Array of stage thresholds for WEPL analysis
  thr[0] = 1.0;
  parser.add_opt('0', "thr0").stow(thr[0]).help("stage 0 threshold (MeV)", "FLOAT").show_default();
  cfg.addItem('0', "thr0", thr[0]);

  thr[1] = 1.0;
  parser.add_opt('1', "thr1").stow(thr[1]).help("stage 1 threshold (MeV)", "FLOAT").show_default();
  cfg.addItem('1', "thr1", thr[1]);

  thr[2] = 1.0;
  parser.add_opt('2', "thr2").stow(thr[2]).help("stage 2 threshold (MeV)", "FLOAT").show_default();
  cfg.addItem('2', "thr2", thr[2]);

  thr[3] = 1.0;
  parser.add_opt('3', "thr3").stow(thr[3]).help("stage 3 threshold (MeV)", "FLOAT").show_default();
  cfg.addItem('3', "thr3", thr[3]);

  thr[4] = 1.0;
  parser.add_opt('4', "thr4").stow(thr[4]).help("stage 4 threshold (MeV)", "FLOAT").show_default();
  cfg.addItem('4', "thr4", thr[4]);

  int useTemp = 1;
  parser.add_opt('E', "useTemp")
    .stow(useTemp) // Add an option to the command line parser
    .help("Use a temporary file instead of local storage for the CALIBRATION process? yes or no", "INT")
    .show_default();
  cfg.addItem('E', "useTemp", useTemp); // Also add the option to the list used for parsing the config file

  int dodEEFilter = 1; // changed default to yes
  parser.add_opt('e', "dEEFilter")
      .stow(dodEEFilter) // Add an option to the command line parser
      .help("Use the dEE filter to filter out nuclear interactions/fragments? 1 or 0", "INT")
      .show_default();
  cfg.addItem('e', "dEEFilter", dodEEFilter); // Also add the option to the list used for parsing the config file

  parser.add_opt_help();
  parser.add_opt_version(version);

  // Read the default configuration from the config file
  if (cfg.Configure() != 0) {
    cout << "Was not able to read a default configuration from " << configFile << endl;
    cout << "The hardwired default configuration will be used." << endl;
  }

  // Parse the command options
  try { parser.parse(argc, argv); }
  catch (arg::Error e) {
    cout << " Error parsing command line: " << e.get_msg() << '\n';
    return 1;
  }

  // Printing out some number
  //////////////////////////////////////////////////////
  if (numbTkrFPGA  != 12) cout << "Non-standard number of tracker FPGAs in the readout = " << numbTkrFPGA << endl;
  if (numbEdetFPGA != 2)  cout << "Non-standard number of energy detector FPGAs in the readout = " << numbEdetFPGA << endl;
    
  // Weird validation
  if (reCalibrate) cout << "Energy detector stage gains will be recalibrated on the fly during processing." << endl;
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


  if (partType == "p" || partType == "proton" || partType == "hydrogen" || partType == "Hydrogen") partType = "H";
  if (partType == "helium" || partType == "alpha" || partType == "Helium") partType = "He";    
  if (partType == "H" || partType == "He") {
    cout << "The beam particle type assumed for analysis and calibration is " << cfg.item_str["partType"] << endl;
  } else {
    cout << "Unrecognized particle type '" << partType << "'" << endl;
    exit(2);
  }

  // The user has to enter the full filename including path
  string CalFile = "CalFileList.txt"; // For calibration runs the input filenames are taken from here
  string inputFileName;

  // Get the list of required, position-sensitive arguments, in this case just the input filename
  vector<string> requiredArgs = parser.args();
  if (requiredArgs.size() == 0) {
    if (!Calibrate) {
      cout << "pCT_Preprocessing: no input raw data file was specified!\n";
      exit(1);
    }
  } else {
    inputFileName = requiredArgs[0];
    CalFile = requiredArgs[0];
  }
  if (Calibrate) cout << "Calibration run.  The list of input files is from " << CalFile << endl;
  else cout << "Preprocessing run, the input file name is " << inputFileName << endl;
  cout << "The TV calibration file is " << TVcorrFile << endl;
  cout << "The WEPL calibration file is " << WcalibFile << endl;

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

  if (Calibrate) cout << "Set the number of events to plot > 0 to get loads of debug histograms in calibration runs." << endl;
  cout << "Fraction of the input file to be analyzed is " << fileFraction << endl;

  cout << "The phantom size for preprocessing is assumed to be " << phantomSize << " mm in radius." << endl;

  if (KillCh != "null")
    cout << "The list of tracker channels to kill will be taken from file " << KillCh << endl;
  if (dodEEFilter) {
    if (!Calibrate)
      cout << "The dE-E filtering of nuclear interactions will be used before WEPL reconstruction" << endl;
    else if (!Calibrate && partType == "He")
      cout << "WARNING: helium fragments will be included in the analysis!" << endl;
  } else
    cout << "No dE-E filtering of nuclear interactions will be used" << endl;

    cfg.addItem('g', "inputFileName", inputFileName);

  
  ////////////////////////////////////////////////////
  // Calibration run
  ////////////////////////////////////////////////////
  if (Calibrate) { // Calibration Run
    cout << "Running the pCT TV and WEPL calibration sequence" << endl;
    if (Normalize) cout << "Will normalize columns in the WET vs E plot to have equal area." << endl;

    // Few more options to be used later in the calibration
    int Nbricks = 1; 
    cfg.addItem('b', "NBricks", Nbricks); 
    int doGains = 1; 
    cfg.addItem('g', "doGains", doGains);
    
    pCTcalib calibProcessor(cfg, CalFile);
    if (calibProcessor.TVmapper() == 0) { // First the TVmapper
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
    cout << "Executing a pCT data pre-processing run" << endl;
    // Here we call the complete preprocessing program
    Preprocessing pCTpreprocessor(cfg);
    int errorCode = pCTpreprocessor.ProcessFile(phantomSize, partType, wedgeOff, fileFraction, numbTkrFPGA, numbEdetFPGA, KillCh);
    return errorCode;
  }

} // end of the main program
