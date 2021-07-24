#ifndef _imadreader_h
#define _imadreader_h

#include <stdint.h>
#include <stdlib.h>
#include "tdcrb.hh"

#include <vector>
#include <map>
#include <string>
#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
class imadreader : public rbProcessor
  {
  public:
    imadreader();
    virtual void init(uint32_t run=0);
    virtual void end(uint32_t run=0);
    virtual  void processEvent(rbEvent* e);
    virtual  void processRunHeader(std::vector<uint32_t> header);
    virtual void loadParameters(Json::Value params);
  private:
    uint32_t _run,_event,_totalSize;
    int32_t _fdOut;
    bool _started,_dummy;
  };

#endif
