#ifndef _zdaq_tdcrb_h
#define _zdaq_tdcrb_h

#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include "jsonGeo.hh"
#include "DCHistogramHandler.hh"
#include "zmBuffer.hh"
#ifndef FEBCMS
#define TDC_TRIGGER_CHANNEL 24
#else
#define TDC_TRIGGER_CHANNEL 0
#endif
#include "rbevent.hh"


class tdcrb
{
public:
  tdcrb(std::string dire="/tmp");
  void registerProcessor(std::string name);
  void Read();
  void end();
  void open(std::string name);
  void close();
  void read();
   

  void geometry(std::string name);
  inline jsonGeo* getGeometry(){ return _geo;}
  uint32_t totalSize();
  uint32_t eventNumber();
  uint32_t runNumber();
  void addRun(uint32_t r,std::string name) { _files.push_back(std::pair<uint32_t,std::string>(r,name));}
  void addFiles();
  void setRun(int r){_run=r;}
  void setOutFileId(int32_t fid){_fdOut=fid;}

  void clearDataSet(){_files.clear();}
  void findDataSet(std::string dir,uint32_t run);
  void stop(){_started=false;}
  void processRawEvent(uint64_t idx);
  void clearShm();
  void monitor();
  void pull(std::string name,zdaq::buffer* buf,std::string sourcedir);
  uint32_t numberOfDataSource();

  void setNFirst(uint32_t f) {_nrfirst=f;}
  void setNMax(uint32_t f) {_nrmax=f;}
private:
  std::vector<std::pair<uint32_t,std::string> > _files;
  uint64_t _bxId,_bxId0;
  uint32_t _gtc;
  double _t,_t0,_tspill;
  std::string _directory;
  uint32_t _run,_event,_totalSize,_nread,_nrmax,_nrfirst;
  uint32_t _nevt,_ntrigger,_nfound,_nbside;
  int32_t _fdIn,_fdOut;
  bool _started;
  unsigned char _buf[32*1024*1024];
  uint32_t _idx;
  jsonGeo* _geo;
  double _readoutTime,_readoutTotalTime;

  DCHistogramHandler* _rh;
  uint32_t _runType,_dacSet,_vthSet,_mezzanine,_difId;
  
  std::map<uint64_t,std::vector<zdaq::buffer*> > _eventMap;
  rbEvent _theEvent;
  std::vector<rbProcessor* > _processors;
};
#endif
