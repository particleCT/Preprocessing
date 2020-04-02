// Routines for reading the raw data file. It was put together starting from
// Pier's bit parsing functions, but this also includes
// a data structure to hold the raw data for an event in a convenient, organized
// manner, allowing a clean interface
// to other parts of the program and to the private user analysis.
// R. Johnson   5/21/2016

#include "pCTraw.h"
#include "BadEvent.h"
#include <fstream>

std::ofstream pCTrawLogFile;
pCTraw::~pCTraw() {pCTrawLogFile.close();}
pCTraw::pCTraw(FILE *in_file, size_t file_size, int thread, int numbTkrFPGA, int numbEdetFPGA) { // Class constructor
  std::cout << "pCTraw constructor thread " << thread << ", file_size=" << file_size << std::endl;
  
  pCTrawLogFile.open("pCTRaw.log");  
  BITS_PER_BYTE = 8; // CHAR_BIT;
  stream_position = 0;
  current_bits = 0;
  queued_bits = 0;
  required_bits = 0;
  reduceFlagCheck = -999;
  killTotal = 0;
  if (numbTkrFPGA < 0 || numbTkrFPGA > num_tkr_fpga) this->numbTkrFPGA = num_tkr_fpga;
  else this->numbTkrFPGA = numbTkrFPGA;

  if (numbEdetFPGA < 0 || numbEdetFPGA > num_enrg_fpga) this->numbEdetFPGA = num_enrg_fpga;
    
  else this->numbEdetFPGA = numbEdetFPGA;
  this->in_file = in_file;
  this->file_size = file_size;
  if (this->in_file == NULL) {
    perror(" pCTraw: Input file is invalid.");
    exit(1);
  }
  event_counter = 0;
  stop_reading = false;
  threadNumber = thread;
  nDaqErr = 0;
  Months[0] = "Jan";
  Months[1] = "Feb";
  Months[2] = "Mar";
  Months[3] = "Apr";
  Months[4] = "May";
  Months[5] = "Jun";
  Months[6] = "Jul";
  Months[7] = "Aug";
  Months[8] = "Sep";
  Months[9] = "Oct";
  Months[10] = "Nov";
  Months[11] = "Dec";
};

// ******************************* ******************************* *******************************
// end of the pCTraw  constructor
// ******************************* ******************************* *******************************

void pCTraw::pCTkillStrip(int FPGA, int chip, int channel) { // Kill a (noisy) channel in the
                                                             // tracker readout
  killCh[FPGA][chip].push_back(channel);
  killTotal++;
  pCTrawLogFile << "pCTraw for thread " << threadNumber << " killing tracker channel " << channel << " in chip " << chip << " of FPGA " << FPGA << std::endl;
}

// ******************************* ******************************* *******************************
// end of the pCTraw  killstrip
// ******************************* ******************************* *******************************

bool pCTraw::parseDate(int &year, int &month, int &day) {

  std::string Date = start_time;
  // pCTrawLogFile << "pCTraw::parseDate: input date string is " << Date <<
  // std::endl;
  std::vector<std::string> tokens;
  size_t found = Date.find_first_not_of(" ");
  while (found != Date.npos) {
    Date = Date.substr(found);
    found = Date.find_first_of(" ");
    if (found == Date.npos)
      found = Date.size();
    std::string token = Date.substr(0, found);
    tokens.push_back(token);
    Date = Date.substr(found);
    found = Date.find_first_not_of(" ");
  }
  if (tokens.size() > 4) {
    month = 13;
    for (int m = 0; m < 12; m++) {
      std::string tokenSub = tokens[1].substr(0, 3);
      if (Months[m].compare(tokenSub) == 0) {
        month = m + 1;
        break;
      }
    }
    if (month > 12) {
      pCTrawLogFile << "pCTraw::parseDate: no match for the month in " << tokens[1] << std::endl;

      return false;
    }
    char *pEnd;
    day = strtol(tokens[2].c_str(), &pEnd, 10);
    if (day == 0) {
      pCTrawLogFile << "pCTraw::parseDate: no match for the day in " << tokens[2] << std::endl;
      return false;
    }
    year = strtol(tokens[4].c_str(), &pEnd, 10);
    if (year == 0) {
      pCTrawLogFile << "pCTraw::parseDate: no match for the year in " << tokens[4] << std::endl;
      return false;
    }
    return true;
  } else {
    pCTrawLogFile << "pCTraw::parseDate: not enough tokens found in date " << start_time << std::endl;
    return false;
  }
}

// Method to print out all of the raw data for a single event (except for
// individual energy samples)
void pCTraw::dumpEvt() {
  pCTrawLogFile << "Dump of event number " << event_number << ", event count " << event_counter << "\n";
  if (DAQ_error)
    pCTrawLogFile << "  ** There was a DAQ error present in the event\n";
  pCTrawLogFile << "  Time tag=" << time_tag << " Delta-t=" << delta_t << "  Raw length=" << raw_length << "\n";
  pCTrawLogFile << "  Error flags=" << bad_fpga_address << " " << tag_mismatch << " " << CRC_error << "\n";
  for (int i = 0; i < numbEdetFPGA; i++) {
    int nsamp = enrg_fpga[i].num_samples;
    int nchan = enrg_fpga[i].num_channels;
    pCTrawLogFile << "  Energy FPGA " << i << "  peds_out= " << enrg_fpga[i].peds_out << " error= " << enrg_fpga[i].error
              << " type= " << enrg_fpga[i].type << "\n";
    pCTrawLogFile << "    tag=" << enrg_fpga[i].tag << "  Number samples=" << nsamp << "  Number channels=" << nchan
              << "\n";
    for (int j = 0; j < nchan; j++) {
      pCTrawLogFile << "      " << j << "  ped=" << enrg_fpga[i].pedestal[j] << "  ph=" << enrg_fpga[i].pulse_sum[j]
                << "\n";
    }
  }
  for (int i = 0; i < numbTkrFPGA; i++) {
    int nchips = tkr_fpga[i].num_chips;
    pCTrawLogFile << "  Tracker FPGA " << i << "  Number chips=" << nchips << "  error=" << tkr_fpga[i].error
              << "  tag=" << tkr_fpga[i].tag << "\n";
    for (int j = 0; j < nchips; j++) {
      int nclus = tkr_fpga[i].chip[j].num_clusts;
      int address = tkr_fpga[i].chip[j].address;
      pCTrawLogFile << "    " << j << "  chip " << address << " # clust=" << nclus
                << " Overflow=" << tkr_fpga[i].chip[j].cluster_overflow << "  error=" << tkr_fpga[i].chip[j].error
                << "  parity error=" << tkr_fpga[i].chip[j].parity_error << "\n";
      for (int c = 0; c < nclus; c++) {
        pCTrawLogFile << "          cluster " << c << " 1st strip=" << tkr_fpga[i].chip[j].cluster[c].first
                  << "   # strips=" << tkr_fpga[i].chip[j].cluster[c].length + 1 << "\n";
      }
    }
  }

}; // End of event dump

// ******************************* ******************************* *******************************
// end of the pCTraw parseDate
// ******************************* ******************************* *******************************

void pCTraw::readRunHeader(const char *inFileName) { // This is called once after opening the input file.
  // This for loop extracts the first 8 bytes from the data file
  for (int j = 0; j < 2; j++) {
    read_append_data(in_file, current_bits, queued_bits, stream_position, file_size, stop_reading);
  }

  ///////////////////////////////////////
  // Reading the FILE header
  // The file header is made by 4 bytes identifing a pCT file, 32 bits are
  // needed
  ///////////////////////////////////////

  // Search for the run header start string
  if (!findRunHdr()) pCTrawLogFile << "Unable to find the run header start string in file " << inFileName << std::endl;

  if (read_file_header(extracted_bits, inFileName) == 0) {

    // reading Run Number: 24 bits
    required_bits = 24;
    extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
    run_number = read_run_number(extracted_bits);
    pCTrawLogFile << "Input data file run header found.  Run number = " << run_number << std::endl;

    // reading Run start time:32 bits
    required_bits = 32;
    extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
    start_time = read_run_startTime(extracted_bits);
    study_date = (int)(extracted_bits); // Needed in integer format by the final output file

    // reading status bits:8 bits
    required_bits  = 8;
    extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
    int flgMsk[2] = { 0x01, 0x02 };
    TimeTags = extracted_bits & flgMsk[0];
    pCTrawLogFile << "The input data file is assumed to be real data from the Phase-II pCT scanner.\n";

    // bool time_tags = read_statusBits(extracted_bits);
    if (TimeTags) pCTrawLogFile << "Time tags are included in the data stream.\n";
    else pCTrawLogFile << "From the run header, time tags are NOT included in the data stream!\n";

    // reading program version number:8 bits
    required_bits  = 8;
    extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
        
    program_version = read_programVersion(extracted_bits);
    pCTrawLogFile << "Event builder FPGA firmware version number = " << program_version << std::endl;

    // reading projection_angle:12 bits
    required_bits  = 12;
    extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
      
    stage_angle = read_projection_angle(extracted_bits) / 10.;
    pCTrawLogFile << "Projection angle from the run header = " << stage_angle << " degrees.\n";
  } else {
    pCTrawLogFile << "Failed to find the run header in the pCT raw data file " << inFileName << std::endl;
    pCTrawLogFile << "Will rewind the file and try to proceed without the run "
                 "header information. . ." << std::endl;
    pCTrawLogFile << "The data are assumed to be real and to include time tags." << std::endl;

    TimeTags = true;
    rewind(in_file);
  }
};
// ******************************* ******************************* *******************************
// end of the pCTraw readRunHeader
// ******************************* ******************************* *******************************

// Method to search for a run header
bool pCTraw::findRunHdr() {
  bool Eureka = false;
  required_bits = 24;
  while (!Eureka) {
    if (stream_position >= file_size)
      return Eureka;
    if (queued_bits < required_bits)
      read_append_data(in_file, current_bits, queued_bits, stream_position, file_size, stop_reading);
    unsigned long long temp_container = current_bits;
    unsigned int temp_queued = queued_bits;
    extracted_bits =
        extract_N_bits(in_file, stream_position, current_bits, queued_bits, required_bits, file_size, stop_reading);
    Eureka = read_BegOfRun(extracted_bits, current_bits, temp_container, queued_bits, temp_queued);
  }
  return Eureka;
};

// ******************************* ******************************* *******************************
// end of the pCTraw findRunHeader
// ******************************* ******************************* *******************************
// Method to search for the beginning of an event
bool pCTraw::findEvtHdr(bool debug) {
  // reading Beginning-of-Event identifier: 24 bits. the string is: 1111 0000 0100 0011 0101 0100 = "1pCT" in ASCII std= 15745876(int dec)
  bool Eureka = false;
  required_bits = 24;
  while (!Eureka) {
    if (stream_position >= file_size) return Eureka; // EOF
    if (queued_bits < required_bits) read_append_data(in_file, current_bits, queued_bits, stream_position, file_size, stop_reading);
      
    unsigned long long temp_container = current_bits;
    unsigned int temp_queued = queued_bits;
    extracted_bits = extract_N_bits(in_file, stream_position, current_bits, queued_bits, required_bits, file_size, stop_reading);
      
    //if (debug) pCTrawLogFile << std::hex << "   extracted_bits=" << extracted_bits << " current_bits=" << current_bits
    //<< " queued_bits=" << queued_bits << std::dec << std::endl;
    Eureka = read_BegOfEvent(extracted_bits, current_bits, temp_container, queued_bits, temp_queued, event_counter, debug);

    //if (debug) pCTrawLogFile << "   stream_position=" << stream_position << " file_size=" << file_size << std::endl;
      
  }
  evtStart = stream_position;
  DAQ_error = false;
  return Eureka;
};

// ******************************* ******************************* *******************************
// end of the pCTraw findEvtHeader
// ******************************* ******************************* *******************************
void pCTraw::readOneEvent(bool debug) { // This is called once for each event after the event header has been located.
  // reading the  event time tag: 36 bits
  if (TimeTags) {
    required_bits  = 36;
    extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
    time_tag       = read_timeTag(extracted_bits);
    //if (debug) pCTrawLogFile << "   Time tag = " << time_tag << std::endl;
  }else time_tag   = 0;
    
  // reading the time since previous trigger: 12 bits
  required_bits    = 12;
  extracted_bits   = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
  delta_t          = read_timeDelta(extracted_bits);
  //if (debug) pCTrawLogFile << "   Time since last trigger = " << delta_t << std::endl;

  // reading Event header: 24 bits
  required_bits    = 24;
  extracted_bits   = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
    
  int TrgBits;
  //if (debug) pCTrawLogFile << std::hex << "   Event header bits = " << extracted_bits << std::dec << std::endl;
  event_number = read_eventHeader(extracted_bits, TrgBits);
  //if (debug) pCTrawLogFile << "Event number: " << event_number << std::endl;
    
  int trgMsk[6] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20 };
  for (int i = 0; i < 6; i++)
    trigger_bits[i] = (trgMsk[i] & TrgBits) != 0; // Will be all zero for program_version<64

  /*if (debug) {
    pCTrawLogFile << "Entering pCTraw::readOneEvent, TimeTags=" << TimeTags << " Run Number=" << run_number
              << " Event number=" << event_number;
    pCTrawLogFile << " Program version=" << program_version << " Stage angle=" << stage_angle;
    pCTrawLogFile << " Start time=" << start_time << " Event count=" << event_counter << std::endl;
    }*/


  
  // reading 12 tracker FPGA headers with 12 bits each
  required_bits = 12;
  for (int FPGA_numb = 0; FPGA_numb < num_tkr_fpga; FPGA_numb++) { // All 12 headers will always be there
    size_t strip_counter = 0;      
    int tkr_tag;
    int tkr_err;
    extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
    int hit_chips  = read_FPGAsHeader(extracted_bits, FPGA_numb, tkr_tag, tkr_err);
    //if (debug) pCTrawLogFile << "    FPGA " << FPGA_numb << " Number of chips hit=" << hit_chips << std::endl;
      
    tkr_fpga[FPGA_numb].num_chips = hit_chips;
    tkr_fpga[FPGA_numb].error     = (tkr_err != 0);
    tkr_fpga[FPGA_numb].tag       = tkr_tag;
    if (hit_chips > 0) {
      required_bits = 12;
      for (int i = 0; i < hit_chips; i++) {
        extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
        int tkerr;
        int tkperr;
        int ovrflw;
        int ASIC_address;
        int clusters_num = read_ASICHeader(extracted_bits, ASIC_address, tkerr, tkperr, ovrflw);
        //if (debug) pCTrawLogFile << "      ASIC address " << ASIC_address << "  Number of clusters= " << clusters_num << std::endl;
          
        tkr_fpga[FPGA_numb].chip[i].address = ASIC_address;
        if (clusters_num <= max_clusts) tkr_fpga[FPGA_numb].chip[i].num_clusts = clusters_num;
        else tkr_fpga[FPGA_numb].chip[i].num_clusts = max_clusts;
          
        tkr_fpga[FPGA_numb].chip[i].error = (tkerr != 0);
        tkr_fpga[FPGA_numb].chip[i].parity_error = (tkperr != 0);
        tkr_fpga[FPGA_numb].chip[i].cluster_overflow = ovrflw;
        required_bits = 12;
        int jcl = 0;
        for (int icl = 0; icl < clusters_num; icl++) {
          extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
          int strip_address;
          int strip_num = read_stripHeader(extracted_bits, strip_address); // if strip number is 0 it means 1 strip is hit
	    
          //if (debug) pCTrawLogFile << " Cluster: # strips= " << strip_num + 1 << " First strip= " << strip_address << std::endl;
          if (icl < max_clusts) {
            bool good = (strip_num <= max_cluster_size);
            if (killTotal > 0 && good) {
              for (int s = 0; s <= strip_num; s++) {
                for (auto ch : killCh[FPGA_numb][ASIC_address]) {
                  if (ch == strip_address + s) {
                    good = false;
                    break;
                  }
                }
                if (!good)
                  break;
              }
            }
            if (good) {
              tkr_fpga[FPGA_numb].chip[i].cluster[jcl].length = strip_num;
              tkr_fpga[FPGA_numb].chip[i].cluster[jcl].first  = strip_address;
              jcl++;
            } else {
              tkr_fpga[FPGA_numb].chip[i].num_clusts--;
            }
          }
        }
      }
    }
  } // ending the Tracker FPGA reading

  // Read data from the energy detector boards
  int count_channel = 0;
  for (int i = 0; i < numbEdetFPGA; i++) {
    bool pedestal_flag;
    bool reduce_flag;
    unsigned int n_samples;
    required_bits  = 12;
    extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
      
    unsigned int OTR;
    int enrg_tag;
    int channels = read_EnFPGAs(extracted_bits, pedestal_flag, reduce_flag, n_samples, OTR, enrg_tag);
    if (reduceFlagCheck != -999) {
      if ((reduce_flag && reduceFlagCheck == 0) || (!reduce_flag && reduceFlagCheck == 1)) {
        pCTrawLogFile << "pCTraw: event found with an inconsistent flag setting "
                     "for data reduction in the FPGA." << std::endl;
        throw BadEvent("bad data reduction flag");
      }
    } else {
      reduceFlagCheck = reduce_flag;
    }
    enrg_fpga[i].peds_out = pedestal_flag;
    if (channels == 0) {
      if (nDaqErr < mxPrnt)
        pCTrawLogFile << "Input data error: Energy FPGA " << i << " number of channels = 0!!!"
                  << "\n";
      if (i == 0) channels = 3;
      else channels = 2;
      nDaqErr++;
      DAQ_error = true;
      throw BadEvent("Zero channels in FPGA");
    }
    enrg_fpga[i].num_channels = channels;
    enrg_fpga[i].type = reduce_flag;
    enrg_fpga[i].tag = enrg_tag;
    enrg_fpga[i].num_samples = n_samples;
    //if (debug)
    // pCTrawLogFile << "  Energy board " << i << "  Number of channels= " << channels << " Reduced=" << reduce_flag
    //            << " Pedestals=" << pedestal_flag << " #samples=" << n_samples << std::endl;
    if (reduce_flag) {
      enrg_fpga[i].OTR[0] = (OTR & 1);
      int OTR1 = ((OTR & 2) >> 1);
      enrg_fpga[i].OTR[1] = (OTR1 != 0);
      int OTR2 = ((OTR & 4) >> 2);
      enrg_fpga[i].OTR[2] = (OTR2 != 0);
      for (int j = 0; j < channels; j++) {
        if (pedestal_flag) {
          required_bits  = 8;
          extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
	    
          int value = (int8_t)extracted_bits;
          // pCTrawLogFile << "    " << j << " pedestal= " <<
          // (bitset<8>)extracted_bits << " v " << value;
          enrg_fpga[i].pedestal[j] = value;
        }
        required_bits  = 16;
        extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
        int value1 = (int16_t)extracted_bits;
        enrg_fpga[i].pulse_sum[j] = value1;
        count_channel++;
        if (pedestal_flag == 0 && i == 1) {
          required_bits  = 4;
          extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
        }
      }
    } else { // Data with samples written out
      int imax[3];
      for (int k = 0; k < channels; k++) {
        enrg_fpga[i].OTR[k] = false;
        imax[k] = 0;
      }
      for (int smp = 0; smp < n_samples; smp++) {
        required_bits  = 3;
        extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
	  
        //if (debug) pCTrawLogFile << "         Sample: " << smp << " " << extracted_bits << " ";
          
        required_bits = 15;
        for (int ch = 0; ch < channels; ch++) {
          extracted_bits = just_read(in_file, file_size, stop_reading, current_bits, queued_bits, required_bits, stream_position);
	    
          unsigned int raw_phsmp = extracted_bits & 0x3fff;
          int otr = (extracted_bits & 0x4000) >> 14;
          int sign = (extracted_bits & 0x2000) >> 13;
          if (otr)
            enrg_fpga[i].OTR[ch] = true;
          int phsmp;
          if (sign) {
            phsmp = (-1) * ((raw_phsmp ^ 0x3fff) + 1); // For negative numbers
          } else {
            phsmp = raw_phsmp;
          }
          enrg_fpga[i].sample[ch][smp] = phsmp;
          //if (debug) pCTrawLogFile << phsmp << "  ";
            
          if (smp == 0) {
            enrg_fpga[i].pedestal[ch] = phsmp;
            enrg_fpga[i].pulse_sum[ch] = phsmp;
          } else {
            if (smp < 9 && phsmp > enrg_fpga[i].sample[ch][imax[ch]]) {
              imax[ch] = smp; // Locate the peak
              enrg_fpga[i].pulse_sum[ch] = enrg_fpga[i].sample[ch][smp - 1] + phsmp; // Add the peak pulse plus the previous one
            }
          }
        }
        //if (debug) pCTrawLogFile << std::endl;
          
      }
      for (int ch = 0; ch < channels; ch++) { // 6-sample around-the-peak data reduction algorithm, as in
                                              // the event builder.
        int endsmp;
        if (imax[ch] + 4 < enrg_fpga[i].num_samples)
          endsmp = imax[ch] + 5;
        else
          endsmp = enrg_fpga[i].num_samples;
        for (int smp = imax[ch] + 1; smp < endsmp; smp++) {
          enrg_fpga[i].pulse_sum[ch] += enrg_fpga[i].sample[ch][smp];
        }
        //if (debug) pCTrawLogFile << "    Channel " << ch << " sample sum=" << enrg_fpga[i].pulse_sum[ch] << std::endl;
          
      }
    }
  }
  raw_length = stream_position - evtStart;
  };
// ******************************* ******************************* *******************************
// end of the pCTraw readOneEvent
// ******************************* ******************************* *******************************
void pCTraw::doWeStop(int max_events, int max_time) { // Check whether we are done or need to
                                                      // read another event
  if (stream_position > file_size) {
    stop_reading = true;
    pCTrawLogFile << "pCTraw: stopping reading for thread " << threadNumber << " because at end of file\n";
  }
  if (max_events > 0) {
    if (event_counter >= max_events) {
      stop_reading = true;
      pCTrawLogFile << "pCTraw: stopping reading for thread " << threadNumber
                << " because the maximum number of events was reached.\n";
    }
  }
  if (max_time > 0) {
    if (time_tag * 1.0E-8 >= max_time) {
      stop_reading = true;
      pCTrawLogFile << "pCTraw: stopping reading for thread " << threadNumber
                << " because the maximum time stamp was reached.\n";
    }
  }
};

// ******************************* ******************************* *******************************
// end of the pCTraw doWeStop
// ******************************* ******************************* *******************************

// [CEO Jan 2016] Simpler code to reverse order of 4-byte integer
//          Reverse byte order: [0123] --> [3210]
uint32_t reverse_int_bytes(unsigned int x) {
  unsigned int result;
  unsigned char n1[4], n2[4];
  memcpy(n1, &x, sizeof(unsigned int));
  n2[0] = n1[3];
  n2[1] = n1[2];
  n2[2] = n1[1];
  n2[3] = n1[0];
  memcpy(&result, n2, sizeof(unsigned int));
  return result;
}

// ******************************* ******************************* *******************************
// end of the pCTraw reverse int bytes
// ******************************* ******************************* *******************************

void pCTraw::read_append_data(FILE *in_file, unsigned long long &bit_container, unsigned int &num_bits,
                              long long &origin, size_t file_size, bool &stop_reading) {
  unsigned int buffer; // 4 bytes
  unsigned int shift_by = BITS_PER_BYTE * sizeof(buffer); // shift is equivalent to the number of
  // bits in a byte (8) multiplied by the number of bytes composing the buffer (4)= 32 bits in total
  bit_container <<= shift_by; // Shift current_bits over by 4-bytes so another
                              // 4-bytes can be appended to the end

  int ret = fread(&buffer, sizeof(buffer), 1, in_file); // Read next 4-bytes from the current reading position inside the file
  origin = origin + sizeof(buffer);

  auto temp = reverse_int_bytes(buffer); // Reverse the order of the bytes previously read into the buffer
  bit_container |= temp;                 // Add these bits to the end of current_bits into the 4-bytes slot
  // opened up by shifting the existing bits by 4-bytes
  num_bits += shift_by; // Having just read and appended 4-bytes = 32-bits to current_bits,
  // add 32 to the count of # bits represented by current_bits
  if (origin > file_size) {
    stop_reading = 1;
    pCTrawLogFile << "read_append_data:  stopping file reading because we've "
                 "arrived at the end of the file.\n";
  }
}

// ******************************* ******************************* *******************************
// end of the pCTraw read_append_data
// ******************************* ******************************* *******************************
unsigned long long pCTraw::extract_N_bits(FILE *in_file, long long stream_position, unsigned long long &bit_container,
                                          unsigned int &num_bits, unsigned int bits_2_extract, size_t file_size,
                                          bool &stop_reading) {
  unsigned long long extracted_bits = bit_container;
  if (num_bits < bits_2_extract)
    read_append_data(in_file, bit_container, num_bits, stream_position, file_size, stop_reading);
  unsigned int shift_size = num_bits - bits_2_extract;
  extracted_bits >>= shift_size;
  bit_container = bit_container - (extracted_bits << shift_size);
  num_bits -= bits_2_extract;
  return extracted_bits;
}
// ******************************* ******************************* *******************************
// ******************************* ******************************* *******************************
// This function groups together the read and extract functions, and makes the
// check on the queued bits
unsigned long long pCTraw::just_read(FILE *in_file, size_t file_size, bool &stop_reading,
                                     unsigned long long &bit_container, unsigned int &num_bits,
                                     unsigned int bits_2_extract, long long &origin) {
  if (origin > file_size) {
    stop_reading = 1;
    pCTrawLogFile << "just_read thread " << threadNumber << ":  stopping reading because we're at the end of the input "
                                                        "raw data file\n";
    return 0;
  }
  else if (num_bits < bits_2_extract) read_append_data(in_file, bit_container, num_bits, origin, file_size, stop_reading);    
  unsigned long long required_bts = extract_N_bits(in_file, origin, bit_container, num_bits, bits_2_extract, file_size, stop_reading);
  return required_bts;
}


// ******************************* ******************************* *******************************
// ******************************* ******************************* *******************************
// File header function
int pCTraw::read_file_header(unsigned long long fileHeader_bits, const char *namefile) {
  // the number 13784398 corresponds to the bit string: 1101 0010 0101 0101 0100 1110 = "1RUN" in std ASCII (in a fancy
  // way)
  unsigned long long mask = pow(2, 24) - 1;
  if ((fileHeader_bits & mask) == 13784398)
    pCTrawLogFile << "File run header found for the pCT raw data file: " << namefile << std::endl;
  else {
    pCTrawLogFile << "No run header found for the pCT raw data file " << namefile << std::endl;
    perror("Did not find the expected header of a pCT file.\n");
    return -1;
  }
  return 0;
}
// ******************************* ******************************* *******************************
// ******************************* ******************************* *******************************
// Run number function
int pCTraw::read_run_number(unsigned long long runNumber_bits) {
  pCTrawLogFile << "Run Number: " << runNumber_bits << std::endl;
  return runNumber_bits;
}

// ******************************* ******************************* *******************************
// ******************************* ******************************* *******************************
// Start time function
std::string pCTraw::read_run_startTime(unsigned long long startTime_bits) {
  time_t RealTime = startTime_bits;
  std::string Tout = ctime(&RealTime);
  std::string whitespaces(" \t\f\v\n\r");
  std::size_t found = Tout.find_last_not_of(whitespaces);
  if (found != std::string::npos)
    Tout.erase(found + 1);
  else
    Tout.clear();
  pCTrawLogFile << "The run start time is " << Tout << "." << std::endl;
  return Tout;
}

// ******************************* ******************************* *******************************
// ******************************* ******************************* *******************************
// File header function
bool pCTraw::read_BegOfEvent(unsigned long long BegOfEvent_bits, unsigned long long &bit_container,
                             unsigned long long temp_cont, unsigned int &num_bits, unsigned int temp_queu,
                             int &event_counter, bool debug) {
  // the string is: 1111 0000 0100 0011 0101 0100= "1pCT" in ASCII std= 15745876
  bool Found = 0;
  unsigned long long mask = pow(2, 24) - 1;
  if ((BegOfEvent_bits & mask) == 15745876) {
    event_counter++;
    Found = 1;
  } else {
    bit_container = temp_cont;
    num_bits = temp_queu - 4; // I can do 4-bit steps, because the BoE usually
                              // starts at the byte or half byte beginning
  }
  //if (debug)
  // pCTrawLogFile << std::hex << "   read_BegOfEvent: bits=" << (BegOfEvent_bits & mask)
  //           << " bit_container=" << bit_container << std::dec << " num_bits=" << num_bits << std::endl;
  return Found;
}
// ******************************* ******************************* *******************************
// ******************************* ******************************* *******************************
bool pCTraw::read_BegOfRun(unsigned long long BegOfRun_bits, unsigned long long &bit_container,
                           unsigned long long temp_cont, unsigned int &num_bits, unsigned int temp_queu) {
  // the string is: 1101 0010 0101 0101 0100 1110 = d2554e
  bool Found = 0;
  unsigned long long mask = pow(2, 24) - 1;
  if ((BegOfRun_bits & mask) == 0xd2554e) {
    Found = 1;
  } else {
    bit_container = temp_cont;
    num_bits = temp_queu - 4; // I can do 4-bit steps, because the BoR usually
                              // starts at the byte or half byte beginning
  }
  return Found;
}

// Event header function
// the event header is made of 6 pieces (pz):
//    -start bits : 10 (11+22zeros)
//    -4 error flags (1+21,20,19,18, zeros)
//    -18 bits: event counter.(18 ones)
// I use a mask for each piece to extract the desired bits and I shift  right
// the extracted bits to have an integer
int pCTraw::read_eventHeader(unsigned long long eventHeader_bits, int &TrgBits) {
  unsigned long long eventHeader_mask[7] = { 0xC00000, 0x200000, 0x100000, 0x080000, 0x040000, 0x03F000, 0x000FFF };
  if (program_version < 64) {
    eventHeader_mask[5] = 0;
    eventHeader_mask[6] = 0x03FFFF;
  }
  unsigned long long eventHeader_pz[7] = {};
  for (int l = 0; l < 7; l++)
    eventHeader_pz[l] = eventHeader_bits & eventHeader_mask[l];

  if ((eventHeader_pz[0] >>= 22) != 2)
    pCTrawLogFile << "*********** DAQ Error flag: Event header starting bits are "
                 "not 10! ***" << std::endl;
  if ((eventHeader_pz[1] >>= 21) == 1) {
    if (nDaqErr < mxPrnt)
      pCTrawLogFile << "*********** DAQ Error flag: !Incorrect FPGA address received! ***" << std::endl;
    bad_fpga_address = true;
    nDaqErr++;
    DAQ_error = true;
    throw(BadEvent("incorrect FPGA address"));
  } else
    bad_fpga_address = false;
  if ((eventHeader_pz[2] >>= 20) == 1) {
    if (nDaqErr < mxPrnt)
      pCTrawLogFile << "*********** DAQ Error flag: !Tag mismatch error! ***" << std::endl;
    tag_mismatch = true;
    nDaqErr++;
    DAQ_error = true;
    throw(BadEvent("tag mismatch"));
  } else
    tag_mismatch = false;
  if ((eventHeader_pz[3] >>= 19) == 1) {
    if (nDaqErr < mxPrnt)
      pCTrawLogFile << "*********** DAQ Error flag: !CRC error! ***" << std::endl;
    CRC_error = true;
    nDaqErr++;
    throw(BadEvent("CRC error"));
  } else
    CRC_error = false;
  if ((eventHeader_pz[4] >>= 18) == 1) {
    if (nDaqErr < mxPrnt)
      pCTrawLogFile << "*********** DAQ Error flag: !Chip error! ***" << std::endl;
    chip_error = true;
    nDaqErr++;
    DAQ_error = true;
    throw(BadEvent("chip error"));
  } else
    chip_error = false;
  bad_strip_address = false;

  TrgBits = (eventHeader_pz[5] >>= 12);
  // if ((eventHeader_pz[0] >>= 22 )==2 && (eventHeader_pz[5]%1000000)==0)
  // pCTrawLogFile << "Starting event n " << eventHeader_pz[5] << " " <<
  // event_counter << std::endl;
  return eventHeader_pz[6];
}

// FPGAs header function (12 bits)
// FPGAs header is made of 6 pieces (pz):
//    -4 bits : FPGA adress 0-11, 1111+8zeros
//    -4 error flags (1+7,6,5,4, zeros){ignored}
//    -4 bits: cips number.(1111)
// I use a mask for each piece to extract the desired bits and I shift right the
// extracted bits to have an integer
int pCTraw::read_FPGAsHeader(unsigned long long FPGAsHeader_bits, int FPGA_num, int &tag, int &err) {
  unsigned long long FPGAsHeader_mask[4] = { 0xF00, 0xE0, 0x10, 0xF };
  // here I am ignoring some pieces
  unsigned long long FPGAsHeader_pz[4] = {};
  for (int l = 0; l < 4; l++) {
    FPGAsHeader_pz[l] = FPGAsHeader_bits & FPGAsHeader_mask[l];
  }
  tag = FPGAsHeader_pz[1] >>= 7;
  err = FPGAsHeader_pz[2] >>= 4;
  int numbChips = FPGAsHeader_pz[3];
  if (numbChips > max_chips) {
    pCTrawLogFile << "DAQ Error: number of chips " << numbChips << " exceeds the hardware limit " << max_chips << "\n";
    DAQ_error = true;
    numbChips = max_chips;
    throw(BadEvent("too many chips"));
  }
  if (FPGA_num == (FPGAsHeader_pz[0] >>= 8)) {
    //        pCTrawLogFile<<"reading FPGA "<<FPGAsHeader_pz[0]<<"..."<<std::endl;
    return numbChips;
  } else {
    if (nDaqErr < mxPrnt)
      pCTrawLogFile << "DAQ Error: FPGA Address mismatch: " << FPGA_num << " vs " << FPGAsHeader_pz[0] << "\n";
    return numbChips;
    nDaqErr++;
    DAQ_error = true;
    throw(BadEvent("bad FPGA address"));
  }
}

// ASICs header is made of 6 pieces (pz):
//    -2 bits (overflow+unused):  11+10zeros=0xC00
//    -4 bits: number of clusters (111+6zeros)=0x3C0
//    -2bits:Chip error and parity (ignored)=0x30
//    -4 bits: chips number.(1111)=0xF
// I use a mask for each piece to extract the desired bits and I shift  right
// the extracted bits to have an integer
int pCTraw::read_ASICHeader(unsigned long long ASICHeader_bits, int &ASIC_add, int &err, int &perr, int &ovrflw) {
  unsigned long long ASICHeader_mask[5] = { 0xC00, 0x3C0, 0x20, 0x10, 0xF }; // representaion HEX to do a mask
  unsigned long long ASICHeader_pz[5] = {};                                  // 1100 0000 0000
  for (int i = 0; i < 5; i++) {                                              // 0011 1100 0000
    ASICHeader_pz[i] = ASICHeader_bits & ASICHeader_mask[i];                 // 0000 0010 0000
  }                                                                          // 0000 0001 0000
  err = ASICHeader_pz[2] >>= 5;                                              // 0000 0000 1111
  perr = ASICHeader_pz[3] >>= 4;
  ASIC_add = ASICHeader_pz[4];
  if (ASIC_add < 0 || ASIC_add > 11) {
    if (nDaqErr < mxPrnt)
      pCTrawLogFile << "**** DAQ error, pCTraw::read_ASICHeader, ASIC address = " << ASIC_add << " is out of range."
                << std::endl;
    nDaqErr++;
    DAQ_error = true;
    ASIC_add = 0; // Not correct, but at least maybe better than out of range!
    throw(BadEvent("bad ASIC address"));
  }
  ovrflw = ASICHeader_pz[0] >>= 11;
  int nClust = ASICHeader_pz[1] >>= 6;
  if (nClust > max_clusts) {
    if (nDaqErr < mxPrnt)
      pCTrawLogFile << "**** DAQ error, pCTraw::read_ASICHeader, number of clusters = " << nClust << " is greater than "
                << max_clusts << std::endl;
    nDaqErr++;
    DAQ_error = true;
    throw(BadEvent("too many clusters"));
  }
  return (nClust);
}

// strip header is made of 2 pieces (pz):
//    -6 bits :  111111+6zeros = 0xFC0
//    -6 bits :  111111 = 0x3F
// I use a mask for each piece to extract the desired bits and I shift  right
// the extracted bits to have an integer
int pCTraw::read_stripHeader(unsigned long long stripHeader_bits, int &strip_add) {
  unsigned long long stripHeader_mask[2] = { 0xFC0, 0x3F }; // representaion HEX
                                                            // to do a mask
  unsigned long long stripHeader_pz[2] = {};
  for (int l = 0; l < 2; l++)
    stripHeader_pz[l] = stripHeader_bits & stripHeader_mask[l];
  strip_add = stripHeader_pz[1];
  int strip_num = (stripHeader_pz[0] >>= 6);
  if (strip_add + strip_num > 63) {
    if (nDaqErr < mxPrnt) {
      pCTrawLogFile << "**** DAQ Error in read_stripHeader****   First strip=" << strip_add
                << "  Number of strips minus 1=" << strip_num << "\n";
    }
    strip_num = strip_num + 63 - (strip_add + strip_num);
    if (nDaqErr < mxPrnt) {
      pCTrawLogFile << "                                         Resetting the "
                   "number of strips minus 1 to " << strip_num << "\n";
    }
    nDaqErr++;
    DAQ_error = true;
    bad_strip_address = true;
  }
  return strip_num;
}

// Energy FPGAs header is made of 6 pieces (pz):
//    -1 bit: pedestal_flag( 1= ped, 0=no ped) 1+12zeros=0x800
//    -4 bits:  front-end tag 1111+7zeros=0x780
//    -1 bit: data type ( 0= samples, 1= reduced) 1+6zeros= 0x40
//    -1 bit: error flag 1+5zeros =0x20
//    -2 bits: number of channel per fpga (3 or 2) 11+000= 0x18
//    -3 bits: OTR( not used if reduced is active)111=7
// I use a mask for each piece to extract the desired bits and I shift right
// the extracted bits to have an integer
int pCTraw::read_EnFPGAs(unsigned long long EnFPGAs_bits, bool &ped_flag, bool &Entype_flag, unsigned int &nsamp,
                         unsigned int &OTR, int &frontend_tag) {
  unsigned long long EnFPGAs_mask[8] = { 0x800, 0x780, 0x40, 0x20, 0x18, 7, 0x1f, 0xf80 }; // representation HEX to do a
                                                                                           // mask
  unsigned long long EnFPGAs_pz[8] = {};
  // if (EnFPGAs_bits == 1280) { pCTrawLogFile << "***** pCTraw::read_EnFPGAs, Error flag encountered!\n";}
  for (int l = 0; l < 7; l++)
    EnFPGAs_pz[l] = EnFPGAs_bits & EnFPGAs_mask[l];

  Entype_flag = (EnFPGAs_pz[2] >>= 6);
  ped_flag = (EnFPGAs_pz[0] >>= 11);
  bool EnFPGAs_err = (EnFPGAs_pz[3] >>= 5);
  int channels_num;
  if (Entype_flag) {
    channels_num = (EnFPGAs_pz[4] >>= 3);
    frontend_tag = (EnFPGAs_pz[1] >>= 7);
    nsamp = 0;
  } else {
    channels_num = 3;
    nsamp = EnFPGAs_pz[6];
    EnFPGAs_pz[7] = EnFPGAs_bits & EnFPGAs_mask[7];
    frontend_tag = (EnFPGAs_pz[7] >>= 7);
  }
  OTR = EnFPGAs_pz[5];
  if (EnFPGAs_err == 1 || (channels_num != 2 && channels_num != 3) || nsamp > 16) {
    if (nDaqErr < mxPrnt) {
      pCTrawLogFile << "***** pCTraw::read_EnFPGAs, Error flag encountered!\n";
      pCTrawLogFile << "      channels_num=" << channels_num << "  OTR=" << OTR << "  nsamp=" << nsamp
                << "  frontend_tag=" << frontend_tag;
      pCTrawLogFile << "  Entype_flag=" << Entype_flag << "  ped_flag=" << ped_flag << std::endl;
    }
    nDaqErr++;
    DAQ_error = true;
  }
  return channels_num;
}
