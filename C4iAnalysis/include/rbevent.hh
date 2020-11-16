#ifndef _RBEVENT_HH_
#define _RBEVENT_HH_
#include "SdhcalDifAccess.hh"
#include <json/json.h>
#define MAXDIF 256
#define MAXFRAME 128
#define FSIZE 20
using namespace sdhcal;
class rbEvent
{
public:
  rbEvent(){_runType=0;_vthSet=0;}
  inline void init( uint32_t run, uint32_t event,uint32_t gtc,uint64_t abcid)
  {  _run=run;_event=event;_gtc=gtc;_abcid=abcid; memset(_frameCount,0,MAXDIF*sizeof(uint32_t));}
  inline void setCalibrationInfos(uint32_t runtype,uint32_t vth){_runType=runtype;_vthSet=vth;}
  inline uint32_t frameCount(uint32_t d) {return _frameCount[d];}
  inline void setFrameCount(uint32_t d,uint32_t n){_frameCount[d]=n;} 
  inline uint8_t*  frameBuffer(){return _fBuf;}
  inline uint32_t run(){return _run;}
  inline uint32_t event(){return _event;}
  inline uint32_t gtc(){return _gtc;}
  inline uint64_t abcid(){return _abcid;}
  inline uint32_t iPtr(uint32_t d,uint32_t f){ return (d*MAXFRAME+f)*FSIZE;}
  inline uint32_t header(uint32_t iptr){return _fBuf[iptr];}
  inline uint32_t bcid(uint32_t iptr){ return DIFUnpacker::getFrameBCID(&_fBuf[iptr]);}
  inline bool  padLevel(uint32_t iptr,uint32_t ipad,uint32_t ilevel){return DIFUnpacker::getFrameLevel(&_fBuf[iptr],ipad,ilevel);}
  inline bool  pad0(uint32_t iptr,uint32_t ipad) {return padLevel(iptr,ipad,0);}
  inline bool  pad1(uint32_t iptr,uint32_t ipad) {return padLevel(iptr,ipad,1);}
  inline int32_t runType(){return _runType;}
  inline int32_t seuil(){return _vthSet;}

private:
  uint32_t _run,_event,_gtc,_runType,_vthSet;
  uint64_t _abcid;
  uint32_t _frameCount[MAXDIF];
  uint8_t _fBuf[MAXDIF*MAXFRAME*FSIZE];
};


class rbProcessor
{
public:
  virtual void init(uint32_t run=0)=0;
  virtual void end(uint32_t run=0)=0;
  virtual  void processEvent(rbEvent* e)=0;
  virtual  void processRunHeader(std::vector<uint32_t> header)=0;
  virtual  void loadParameters(Json::Value params)=0;
};
#endif
