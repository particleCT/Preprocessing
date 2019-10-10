#ifndef pCTraw_h
#define pCTraw_h
// Class for reading the raw data file. It was put together starting from
// Blake's and Pier's bit parsing functions, but this also includes
// a data structure to hold the raw data for an event in a convenient, organized
// manner, allowing a clean interface
// to other parts of the program and to the private user analysis.
// R. Johnson   5/21/2016
#include <iostream>
#include <cstdlib>
#include <cstddef>
#include <cstdio>
#include <string>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <vector>
#include <cmath>

class pCTraw {

  typedef unsigned char byte;
  std::size_t BITS_PER_BYTE;

  long long stream_position;
  unsigned long long current_bits; // holds the ULL representation of all the
                                   // bits read into memory
                                   // and placed in the correct order
  unsigned int queued_bits;        // # of bits currently represented by the ULL
                                   // variable "current_bits"
  // (i.e. # bits in queue waiting to be processed)
  unsigned long long extracted_bits; // bit containing information needed
  int required_bits;
  int reduceFlagCheck;

  FILE *in_file; // pointer to the input data file
  size_t file_size;
  long long evtStart;
  int threadNumber;
  int nDaqErr;

  // This union is used to reverse byte order
  union int_reversal {
    unsigned int int_val;
    struct {
      unsigned short int_hi;
      unsigned short int_lo;
    };
    char each_byte[4]; // unused
  };

public:
#define num_tkr_fpga 12
#define num_enrg_fpga 2
#define max_chips 12
#define max_clusts 10
#define mxPrnt 10000
#define max_cluster_size 9

  // Lists of tracker channels to suppress
  int killTotal;
  std::vector<int> killCh[num_tkr_fpga][max_chips];

  // Run header information extracted from the input file:
  int run_number;         // Filled only if run header found
  std::string start_time; // Filled only if run header found
  int study_date;         // Filled only if run header found (integer representation of
                          // the start_time)
  float stage_angle;      // Filled only if run header found
  int program_version;    // Filled only if run header found
  bool TimeTags;          // Assumed true if no run header is found
  // Information from a single event.  This is just the raw data organized for
  // convenient access
  bool DAQ_error;
  int event_number;            // Warning, this can repeat in a long run
  bool trigger_bits[6];        // Energy detector trigger bits (all 0 for
                               // program_version<64)
  int raw_length;              // Number of bits in the raw event
  unsigned long long time_tag; // Time tag
  int delta_t;                 // Time in 10ns units since previous trigger condition
  bool bad_fpga_address;
  bool tag_mismatch;
  bool CRC_error;
  bool chip_error;
  bool bad_strip_address;
  struct EnergyFPGA {
    bool peds_out; // 0=no pedestals written out, 1=pedestals written
    // Only relevant for reduced data
    bool error;
    bool OTR[3]; // out-of-range indicators per channel; calculate for
                 // non-reduced data
    int tag;
    bool type;         // 0=samples written out, 1=reduced in FPGA
    int num_samples;   // 0 if samples are not written out
    int num_channels;  // 2 or 3 for reduced data, always 3 for samples
    int sample[3][16]; // up to 16 samples for each channel; none for reduced
                       // data
    int pedestal[3];   // equals the first sample
    int pulse_sum[3];  // calculate from the samples for non-reduced data
  } enrg_fpga[num_enrg_fpga];
  struct TrackerFPGA {
    bool error;
    int tag;
    int num_chips;
    struct TrackerChip {
      bool cluster_overflow;
      bool error;
      bool parity_error;
      int address;
      int num_clusts;
      struct TrackerCluster {
        int length;          // Number of strips in cluster (1 to 64)
        int first;           // Address of the first strip in cluster (0 to 63)
      } cluster[max_clusts]; // Entries beyond num_clusts are unfilled
    } chip[max_chips];       // Entries beyond num_chips are unfilled
  } tkr_fpga[num_tkr_fpga];  // Index is the FPGA address

  int event_counter;
  bool stop_reading;
  std::string Months[12];

  pCTraw(FILE *in_file, size_t file_size, int thread, int numbTkrFPGA, int numbEdetFPGA);
  void pCTkillStrip(int FPGA, int chip, int channel);
  void dumpEvt();
  void readRunHeader(const char *inFileName);
  bool findEvtHdr(bool debug);
  void readOneEvent(bool debug);
  void doWeStop(int max_events, int max_time);
  bool parseDate(int &year, int &month, int &day);

private: // Blakes's and Piersimoni's bit parsing methods, plus more, are all
         // hidden here:
  int numbTkrFPGA;
  int numbEdetFPGA;

  bool findRunHdr();

  void read_append_data(FILE *in_file, unsigned long long &bit_container, unsigned int &num_bits, long long &origin,
                        size_t file_size, bool &stop_reading);

  unsigned long long extract_N_bits(FILE *in_file, long long stream_position, unsigned long long &bit_container,
                                    unsigned int &num_bits, unsigned int bits_2_extract, size_t file_size,
                                    bool &stop_reading);

  unsigned long long just_read(FILE *in_file, size_t file_size, bool &stop_reading, unsigned long long &bit_container,
                               unsigned int &num_bits, unsigned int bits_2_extract, long long &origin);

  int read_file_header(unsigned long long fileHeader_bits, const char *namefile);

  int read_run_number(unsigned long long runNumber_bits);

  std::string read_run_startTime(unsigned long long startTime_bits);

  inline bool read_statusBits(unsigned long long status_bits) {
    // if (status_bits==1) cout << "Time tags are included in the data stream."
    // << endl;
    // else cout << "Time tags are NOT included in the data stream." << endl;
    return status_bits;
  }

  inline int read_programVersion(unsigned long long programVersion_bits) {
    // cout << "The pCT scanner event builder firmware program version is " <<
    // programVersion_bits << endl;
    return programVersion_bits;
  }

  // Projection angle
  inline int read_projection_angle(unsigned long long projectionangle_bits) {
    // cout << "The stage angle is " << projectionangle_bits << endl;
    return projectionangle_bits;
  }

  bool read_BegOfEvent(unsigned long long BegOfEvent_bits, unsigned long long &bit_container,
                       unsigned long long temp_cont, unsigned int &num_bits, unsigned int temp_queu, int &event_counter,
                       bool debug);

  bool read_BegOfRun(unsigned long long BegOfRun_bits, unsigned long long &bit_container, unsigned long long temp_cont,
                     unsigned int &num_bits, unsigned int temp_queu);

  inline unsigned long long read_timeTag(unsigned long long timeTag_bits) {
    // cout << "Event time tag (in clock cycles): " << timeTag_bits << endl;
    return timeTag_bits;
  }

  // Time since previous trigger function
  inline int read_timeDelta(unsigned long long timeDelta_bits) {
    // cout << "Time since the previous trigger (in clock cycles): " <<
    // timeDelta_bits  << endl;
    return timeDelta_bits;
  }

  int read_eventHeader(unsigned long long eventHeader_bits, int &TrgBits);

  int read_FPGAsHeader(unsigned long long FPGAsHeader_bits, int FPGA_num, int &tag, int &err);

  int read_ASICHeader(unsigned long long ASICHeader_bits, int &ASIC_add, int &err, int &perr, int &ovrflw);

  int read_stripHeader(unsigned long long stripHeader_bits, int &strip_add);

  int read_EnFPGAs(unsigned long long EnFPGAs_bits, bool &ped_flag, bool &Entype_flag, unsigned int &nsamp,
                   unsigned int &OTR, int &frontend_tag);

}; // End of the pCTraw class
#endif
