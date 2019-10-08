/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   BadEvent.h
 * Author: swr
 *
 * Created on January 14, 2019, 11:35 AM
 */

#ifndef BADEVENT_H
#define BADEVENT_H
#include <exception>

class BadEvent : public std::exception {
  std::string why;
  char *strbuf = (char *)malloc(sizeof(char) * 100);

public:
  BadEvent(std::string why) { this->why = why; }
  const char *what() const throw() {
    sprintf(strbuf, "Bad event input exception thrown because %s", why.c_str());
    return strbuf;
  }
};

#endif /* BADEVENT_H */
