#ifndef _UTIL_H_
#define _UTIL_H_

#include <cstdio>
#include <cstdlib>

using namespace std;

class Util { // Just of list of useful functions packaged together in this
             // class.  The constructor is the default.

public:
  // Parse a string of the form "key = value" into the key and value
  static bool getKeyValue(string line, string &key, string &value) {
    size_t found = line.find_first_not_of(" #");
    if (found != line.npos) {
      line = line.substr(found);
      found = line.find_first_of(" =");
      if (found != line.npos) {
        key = line.substr(0, found);
        found = line.find_first_of("=");
        if (found != line.npos) {
          line = line.substr(found + 1);
          found = line.find_first_not_of(" ");
          if (found != line.npos) {
            line = line.substr(found);
            found = line.find_first_of(" ");
            if (found == line.npos)
              found = line.size();
            value = line.substr(0, found);
            return true;
          }
        }
      }
    }
    return false;
  }

  // Turn a date string into an integer to be used for algebraic comparisons
  static bool parseDate(string line, int &year, int &month, int &day) {
    size_t found = line.find_first_not_of(" ");
    if (found != line.npos) {
      line = line.substr(found);
      found = line.find_first_of("/");
      if (found != line.npos) {
        string yearStr = line.substr(0, found);
        line = line.substr(found + 1);
        found = line.find_first_of("/");
        if (found != line.npos) {
          string monthStr = line.substr(0, found);
          line = line.substr(found + 1);
          found = line.find_first_of(" ");
          if (found == line.npos)
            found = line.size();
          string dayStr = line.substr(0, found);
          char *ePtr;
          year = strtol(yearStr.c_str(), &ePtr, 10);
          month = strtol(monthStr.c_str(), &ePtr, 10);
          day = strtol(dayStr.c_str(), &ePtr, 10);
          return true;
        }
      }
    }
    return false;
  }

  // Separate a string into a list of tokens delimited by spaces
  static vector<string> getTokens(string line) {
    vector<string> temp;
    size_t found = line.find_first_not_of(" ");
    while (found != line.npos) {
      line = line.substr(found);
      found = line.find_first_of(" ");
      if (found == line.npos)
        found = line.size();
      temp.push_back(line.substr(0, found));
      line = line.substr(found);
      found = line.find_first_not_of(" ");
    }
    return temp;
  }

  // Change H:M:S to decimal
  static double expandTime(string Time) {
    char *ePtr;
    if (Time.size() < 15) {
      cout << "Util::expandTime:  input string " << Time << " is not long enough." << endl;
      return 0.0;
    }
    if (Time.c_str()[2] != ':' || Time.c_str()[5] != ':' || Time.c_str()[8] != '.') {
      cout << "Util::expandTime:  input string " << Time << " does not have the expected format." << endl;
      return 0.0;
    }
    int Hour = strtol(Time.substr(0, 2).c_str(), &ePtr, 10);
    int Minutes = strtol(Time.substr(3, 4).c_str(), &ePtr, 10);
    int Seconds = strtol(Time.substr(6, 7).c_str(), &ePtr, 10);
    int muSec = strtol(Time.substr(9, 14).c_str(), &ePtr, 10);
    long long All = (Hour * 3600 + Minutes * 60 + Seconds);
    double retVal = (double)All + (double)muSec / 1000000.;
    return retVal;
  }

  // Function to scan the run log file for continuous scans to try to extract
  // the angle of the stage at the time of the start of run
  static float getStartAngle(string logFile) {
    char *ePtr;
    ifstream infile(logFile);
    if (infile) {
      string line;
      double startTime = 0.0;
      while (getline(infile, line)) {
        size_t found = line.find("Run start time is");
        if (found == line.npos)
          continue;
        vector<string> tokens = getTokens(line);
        if (tokens.size() < 2)
          continue;
        startTime = expandTime(tokens.at(1));
        cout << "The start time of the run is " << startTime << " from the log file " << endl;
        break;
      }
      if (startTime <= 0.0) {
        cout << "Util::getStartAngle: unable to find the run start time." << endl;
        infile.close();
        return 0.0;
      }
      while (getline(infile, line)) {
        size_t found = line.find("Current stage rotation:");
        if (found == line.npos)
          continue;
        vector<string> tokens = getTokens(line);
        if (tokens.size() < 9)
          continue;
        double endTime = expandTime(tokens.at(1));
        float stageAngle = atof(tokens.at(8).c_str());
        cout << "From the log file the stage angle at time " << endTime << " is " << stageAngle << endl;
        infile.close();
        float offset = stageAngle - (endTime - startTime) * 6.0;
        cout << "The stage angle at the start of run was " << offset << endl;
        return offset; // Rotation rate of 6 degrees per second.
      }
      infile.close();
      cout << "Util::getStartAngle: unable to find the stage angle." << endl;
      return 0.0;
    } else {
      cout << "Util::getStartAngle: unable to open the log file " << logFile << endl;
      return 0.0;
    }
  }
};

#endif