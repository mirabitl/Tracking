#ifndef _RBEVENT_HH_
#define _RBEVENT_HH_
#include "LirocAccess.hh"
#include <json/json.h>
#define MAXDIF 256
#define MAXFRAME 100000
#define FSIZE 4
using namespace liroc;
class rbEvent
{
public:
  rbEvent(): _difType(false) {_runType=0;_vthSet=0;}
  inline void init( uint32_t run, uint32_t event,uint32_t gtc,uint64_t abcid)
  {  _run=run;_event=event;_gtc=gtc;_abcid=abcid; memset(_frameCount,0,MAXDIF*sizeof(uint32_t));

    if (_difType)
      {
	for (uint32_t i=0;i<theLirocPtrList_.size();i++) 
	  {   
	  
	    delete theLirocPtrList_[i];
	  }
	theLirocPtrList_.clear();

      }
  }
  inline void setCalibrationInfos(uint32_t runtype,uint32_t vth){_runType=runtype;_vthSet=vth;}
  inline void setDifType(bool t) {_difType=t;}
  inline bool isDifType(){return _difType;}
  inline uint32_t frameCount(uint32_t d) {return _frameCount[d];}
  inline void setFrameCount(uint32_t d,uint32_t n){_frameCount[d]=n;} 
  inline uint8_t*  frameBuffer(){return _fBuf;}
  inline uint32_t run(){return _run;}
  inline uint32_t event(){return _event;}
  inline uint32_t gtc(){return _gtc;}
  inline uint64_t abcid(){return _abcid;}
  inline uint32_t iPtr(uint32_t d,uint32_t f){ return (d*MAXFRAME+f)*FSIZE;}
  inline uint32_t header(uint32_t iptr){return _fBuf[iptr];}
  inline uint8_t* frame(uint32_t iptr){return &_fBuf[iptr];}
  inline uint32_t bcid(uint32_t iptr){ return liroc::Unpacker::getFrameBCID(&_fBuf[iptr]);}
  inline uint32_t  channel(uint32_t iptr){return liroc::Unpacker::getFrameChannel(&_fBuf[iptr]);}
  inline  std::map<uint32_t,std::vector<std::pair<liroc::Ptr*,uint32_t> > > &tFrame(){ return _tframe;}
  inline std::map<uint32_t,std::bitset<64> > &tCount(){return _tcount;}
  inline std::vector<liroc::Ptr*> &difList(){return theLirocPtrList_;}
  inline int32_t seuil(){return _vthSet;}

private:
  uint32_t _run,_event,_gtc,_runType,_vthSet;
  uint64_t _abcid;
  uint32_t _frameCount[MAXDIF];
  uint8_t _fBuf[MAXDIF*MAXFRAME*FSIZE];
  bool _difType;
  std::vector<liroc::Ptr*> theLirocPtrList_;
  std::map<uint32_t,std::vector<std::pair<liroc::Ptr*,uint32_t> > > _tframe;
  std::map<uint32_t,std::bitset<64> > _tcount;

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
