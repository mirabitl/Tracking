#ifndef _TDCCHANNEL_HH
#define _TDCCHANNEL_HH
#define TDC_COARSE_TIME 2.5
#define TDC_FINE_TIME 0.009765625
#define COARSEMAX 0XFFFFFF
#include <stdint.h>
#include <stdio.h>
#include "TdcMapping.hh"
#include "jsonGeo.hh"
using namespace std;
namespace lydaq {

#define LASTCHAN 32

  class TdcChannel
{
public:
  TdcChannel() :_fr(NULL),_used(false) {;}
  TdcChannel(uint8_t*  b,uint8_t feb=0) :_fr(b),_used(false),_feb(feb),_0c(0),_0f(0),_0t(0.0) {;}
  inline uint8_t channel() {
    if (falling())
      return (_fr[0]-0x80)-LASTCHAN;
    else
      return _fr[0];}
  inline bool falling() {return  ((_fr[0]&0X80)>>7)==1;}
  inline uint16_t pr() {return  TDC2PR[channel()];}
  inline uint16_t lemo() {return  strip()+12*side();}
  inline uint16_t side() {return  SIDE[channel()];}
  inline uint16_t strip() {return  70+STRIP[channel()];}
  inline uint16_t lemo(jsonFebInfo& f) {return  strip(f)+12*side(f);}
  inline uint16_t side(jsonFebInfo& f) {return  f.tdc2side[channel()];}
  //  inline uint16_t strip(jsonFebInfo& f) {return  70+f.tdc2strip[channel()];}
  inline uint16_t strip(jsonFebInfo& f) {return  f.tdc2strip[channel()];}
  // Old style
  // inline double pedSubTime(jsonFebInfo& f) {return  tdcTime()-f.timeped[channel()];}
  inline double pedSubTime(jsonFebInfo& f) {
    if (side(f)==0)
      return  tdcTime()-f.st0[detectorStrip(f)];
    else
      return  tdcTime()-f.st1[detectorStrip(f)];

  }
  inline uint16_t feb(){return _feb;}


  inline uint16_t detectorStrip(uint32_t feb) {return  strip()+FEB2STRIP[feb];}
  inline uint16_t detectorStrip(jsonFebInfo& f) {return  strip(f)+f.stripShift;}

  inline uint8_t length(){return 6;}
  inline uint64_t coarse() const {return ((uint64_t)_fr[4])|((uint64_t)_fr[3]<<8)|((uint64_t)_fr[2]<<16)|((uint64_t)_fr[1]<<24);}
  inline uint8_t fine() const {return _fr[5];}

  inline uint32_t bcid(){return (uint32_t) (coarse()*TDC_COARSE_TIME/200);}
  inline  double rawTime(uint64_t c, uint8_t f) const { return (c+f/256.0)*TDC_COARSE_TIME;}

  inline  double tdcTime() const { double rt=rawTime(coarse(),fine());
    return (_0t<=rt)?rt-_0t:rt+COARSEMAX*TDC_COARSE_TIME-_0t;}
  inline uint8_t* frame(){ return _fr;}
  inline bool used(){return _used;}
  inline void setUsed(bool t){_used=t;}
  bool operator<(const TdcChannel &ipaddr){
    if( coarse() < ipaddr.coarse())
    return true;
  else
    return false;
    
}
  void dump()
  {
    for (int i=0;i<6;i++)
      printf("%.2x ",_fr[i]);
    printf("%d %ld %d %f \n",channel(),coarse(),fine(),tdcTime());
    printf("\n");
  }
  void setZero(uint64_t c, uint8_t f){_0c=c;_0f=f;_0t=rawTime(c,f);}
private:
  uint8_t* _fr;
  bool _used;
  uint8_t _feb;
  uint8_t _0f;
  uint64_t _0c;
  double _0t;
};


};
#endif
