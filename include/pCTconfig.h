#ifndef PCTCONFIG_h
#define PCTCONFIG_h
// Utility class to set the pCT_Preprocessing command-line defaults from a text
// file
// The text file format is a set of lines of format:
// key = value
// where key and value are text strings that can be interpreted as string,
// integer, or float
// Comment lines beginning with # can be embedded in the file
// R.P. Johnson  10/8/2016

#include "Util.h"
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <fstream> 
using namespace std;

class pCTconfig {
public:
  struct cfgItm {
    string key;
    string longKey;
    string type;
    int *valInt;
    float *valFloat;
    string *valString;
  };
  vector<cfgItm> itemList;
  string configFileName;
  map<string, string> item_str;
  map<string, int> item_int;
  map<string, float> item_float;

  pCTconfig(string fileName) {
    configFileName = fileName;
    return;
  }
  
  void addItem(char key, string longKey, string &value) { // Call this to add a key to the list for string
    cfgItm tmpItem;

    tmpItem.key = key;
    tmpItem.longKey = longKey;
    tmpItem.type = "STRING";
    tmpItem.valString = &value; // Save a pointer to the item value, so that it can be altered later.
    pair<map<string,string>::iterator,bool> ret;
    ret = item_str.insert(pair<string,string>(longKey, value));
    if(!ret.second) item_str[longKey] = value;
    itemList.push_back(tmpItem);
    
  }
  void addItem(char key, string longKey, int &value) { // For the case that the value is integer
    cfgItm tmpItem;
    tmpItem.key = key;
    tmpItem.longKey = longKey;
    tmpItem.type = "INT";
    tmpItem.valInt = &value;
    pair<map<string,int>::iterator,bool> ret;
    ret = item_int.insert(pair<string,int>(longKey, value));
    if(!ret.second) item_int[longKey] = value;
    itemList.push_back(tmpItem);
  }
  void addItem(char key, string longKey, float &value) { // For the case that the value is float
    cfgItm tmpItem;
    tmpItem.key = key;
    tmpItem.longKey = longKey;
    tmpItem.type = "FLOAT";
    tmpItem.valFloat = &value;
    pair<map<string,float>::iterator,bool> ret;
    ret = item_float.insert(pair<string,float>(longKey, value));
    if(!ret.second) item_float[longKey] = value;
    itemList.push_back(tmpItem);
    
  }
  int Configure() { // Read the config file and try to match keys with those in the list that has been built with addItem.
    Util U;
    string line;
    ifstream infile(configFileName);
    int linecount = 0;
    cout << "pCTconfig::Configure: setting option defaults from file " << configFileName << ":" << endl;
    if (infile) {
      while (getline(infile, line)) {
        if (line == "") continue; // Skip blank lines
        size_t found = line.find_first_not_of(" ");
        if (line[found] == '#') continue; // Skip comment lines
        line = line.substr(found);
        string key;
        string value;
        U.getKeyValue(line, key, value);
        for(int i = 0; i < itemList.size(); ++i) {
          if(key.compare(itemList[i].key) == 0 || key.compare(itemList[i].longKey) == 0) {
            if (itemList[i].type == "STRING") {
	      pair<map<string,string>::iterator,bool> ret;
	      ret = item_str.insert(pair<string,string>(key, value));
	      if(!ret.second) item_str[key] = value;
              break;
            } else if (itemList[i].type == "INT") {
	      pair<map<string,int>::iterator,bool> ret;
	      ret = item_int.insert(pair<string,int>(key, stoi(value)));
	      if(!ret.second) item_int[key] = stoi(value);
              break;
	    } else if (itemList[i].type == "FLOAT") {
	      pair<map<string,float>::iterator,bool> ret;
	      ret = item_float.insert(pair<string,float>(key, stof(value)));
	      if(!ret.second) item_float[key] = stof(value);
              break;
            }
	    else cout << "pCTconfig::Configure: no match found for key " << key << endl;
          }
	}
        linecount++;
      }
      return 0;
    }
    return 1;
  }
};
#endif
