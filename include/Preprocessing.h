#ifndef _PREPROCESSING_H_
#define _PREPROCESSING_H_
// Phase-II scanner data preprocessing code, produces the input needed by image
// reconstruction.
// This code also serves for monitoring raw data during runs.
// R.P. Johnson  May 8, 2016
//
// The following classes are used:
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
// - UserAnalysis.h: optional entry points for private user data analysis, to
// keep those hacks separate from public code
// - Histogram.cpp and Histogram.h: simple facilities for making histograms and
// profile plots
//
// - TDB: the geometry and calibration constants need to be accessed from a
// database with date key
//
// Compile it with g++ using -std=c++0x -lpthread -g -rdynamic (the -g and
// -rdynamic can be omitted if all is working well)
//
struct generalparam { // Variables needed by the pCTevents program, put together
                      // to reduce the size of the parameter list.
  // This is to enable passing of the raw data processing task "pCTevents" to multiple threads.
  std::string inFileName; // Raw data file name
  std::string Outputdir;  // User supplied location for output files and input of
                          // calibration files
  int threadNum;          // To identify which thread is printing out stuff
  int n_debug;            // Number of events for debug print-out
  int n_plot;             // Number of events for which to plot tracking information
  int max_events;         // Maximum number of events to process (for each thread)
  int max_time;
  bool continuous_scan; // True if a continuous scan is being executed, in which
                        // case the time-stamp is used to calculate the
                        // projection angle
  float proj_angle;     // Projection angle from the raw data file or else
                        // overridden by a value supplied by the user command line
  int analysisLevel;    // Analysis level supplied by the user command line
  bool callUser;        // Command line option to call, or not call, the UserAnalysis
                        // entry points. Only applies to thread 0 in any case.
  bool reCalibrate;     // Update the gain calibrations on the fly during
                        // processing.
  int pdstlr[5];        // Lower end of the pedestal histogram range
};

class Preprocessing { // Top level program from the pCT preprocessing task.
  generalparam config;
  size_t file_size;
  float Version;
  float beamEnergy;
  int n_threads;
  double Uhit[4];
  std::string study_name;
  std::string Outputdir;
  std::string WcalibFile;
  std::string TVcorrFile;
  std::string OsName;
  float StgThr[5];
  int fileBins;
  int analysisLevel;
  bool callUser;
  bool continuous_scan;
  float initialAngle;
  int max_events;
  int max_time;
  int n_debug;
  int n_plot;
  float proj_angle;
  FILE *in_file;
  char inFileName[256];
  time_t start_time;
  struct tm *now;
  float dTheta;
  bool eventOrder;
  bool energyOutput;
  bool timeStampOutput;
  bool eventIDOutput;
  bool dodEEFilter;

  static int findEvt(FILE *fp);
  void WriteBinaryFile3(bool timeStampOutput, bool energyOutput, bool eventIDOutput, float AngleNb,
                        const char OutputFilename[], const char DATA_SOURCE[], const char PHANTOM_NAME[],
                        int study_date, int event_counter, double u[], float V0[], float V1[], float V2[], float V3[],
                        float T0[], float T1[], float T2[], float T3[], float E1[], float E2[], float E3[], float E4[],
                        float E5[], float WetBinary[], float ProjAngle[], unsigned int TimeStamp[],
                        unsigned int EventIDs[]); // float AngleNb

public:
  Preprocessing(std::string inputFileName, std::string study_name, std::string Outputdir, std::string WcalibFile,
                std::string TVcorrFile, int n_threads, float StgThr[5], int angleBins, int analysisLevel, bool callUser,
                bool continuous_scan, float initialAngle, bool realTimeCal, int max_events, int max_time, int n_debug,
                int n_plot, float proj_angle, bool dodEEFilter, int pdstlr[5], std::string OsName, float BeamEnergy, float Version);

  int ProcessFile(float phantomSize, std::string partType, float wedgeOffset, float fileFraction, int numbTkrFPGA,
                  int numbEdetFPGA, std::string KillCh);
};
#endif
