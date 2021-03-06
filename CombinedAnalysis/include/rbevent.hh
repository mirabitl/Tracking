#ifndef _RBEVENT_HH_
#define _RBEVENT_HH_
#include "SdhcalPmrAccess.hh"
#include "TdcChannel.hh"
#include <json/json.h>
#include <algorithm>
#define MAXDIF 256
#define MAXFRAME 128
#define FSIZE 20
using namespace sdhcal;
#define LPCB 160.
#define VPCB (160. / 8.7)
#define VEM888 15.58
#define VFR4 14.62

namespace Lmana
{
  class TdcStrip
  {
  public:
    TdcStrip() : _ch(0), _dif(0), _str(0), _t0(0), _t1(0), _shift(0) { ; }
    TdcStrip(uint16_t dif, uint16_t st, double t0, double t1, double shift = 0) : _dif(dif), _str(st), _t0(t0), _t1(t1), _shift(shift), _ch(1) { ; }
    TdcStrip(uint16_t ch, uint16_t dif, uint16_t st, double t0, double t1, double shift = 0) : _ch(ch), _dif(dif), _str(st), _t0(t0), _t1(t1), _shift(shift) { ; }
    inline uint16_t strip() const { return _str; }
    inline uint16_t chamber() const { return _ch; }
    inline uint16_t dif() const { return _dif; }
    inline double t0() const { return _t0; }
    inline double t1() const { return _t1; }
    inline double TM() const { return (_t1+_t0)/2.; }
    inline double shift() const { return _shift; }
    //inline double ypos() const {return (_t0-_t1-_shift)/1.;}
    // Avant inline double ypos() const {return _shift+80.0+(160.-(_t1-_t0)*18.39)/2.0;}
    inline double ypos() const {
      float Ls=146.3;
      float lc[33]={0.,
		    3.0,3.0,3.0,3.0,
		    3.0,3.99,4.98,5.97,
		    6.96,7.95,8.94,9.93,
		    10.92,11.91,12.90,13.89,
		    13.89,12.90,11.91,10.92,
		    9.93,8.94,7.95,6.96,
		    5.97,4.98,3.99,3.0,
		    3.0,3.0,3.0,3.0};
      float lr0=35.8+6.9+146.2+8.6+2.5+22;
      float dlr[33]={0.,
		     9.0,8.4,7.8,7.2,6.6,6.0,5.4,4.8,
		     4.2,3.6,3.0,2.4,1.8,1.2,0.6,0.,
		     0.,0.6,1.2,1.8,2.4,3.0,3.6,4.2,
		     4.8,5.4,6.0,6.6,7.2,7.8,8.4,9.0,
		     };
      float Lr=lr0-dlr[_str];
      float Lc=lc[_str];
      float yloc=0;
      if (_str<=16)
	yloc= _shift+(Ls+Lr-Lc-VEM888*(_t1-_t0))/2.0;
      else
	yloc= _shift+(Ls+Lr-Lc-VFR4*(_t1-_t0))/2.0;

      yloc=yloc/cos(abs(_str-16)/96.*3.1415926535897932384626433/9);

      // 1786 yloc = yloc -2.18658 -(10.487-0.159423*yloc);
      // 1818 yloc = yloc -5.6153200 -(10.487-0.159423*yloc);
      return yloc;
	//cos(abs(_str-16)*20./96*3.14159265359/180.)*yloc;
      //return _shift + 160. - (_t1 - _t0) * 18.39 / 2.0;
    }
    inline double xpos() const
    {
      return _str * 0.96-(_str-16)*0.08;
    }
    inline double X(){return xpos();}
    inline double Y(){return ypos();}
  private:
    uint16_t _dif, _str, _ch;
    double _t0, _t1, _shift;
  };
  // DTA 2.5 DTY 1.5 puis 5. 5.
// was 4
#define DTA 3.0
// was 5.5
#define DTY 8.5
  class TdcCluster
  {
  public:
    TdcCluster() : _x(0), _y(0), _t0(0), _t1(0) { _strips.clear(); }
    void Print()
    {
      printf("X %f Y %f Size %ld \n", _x, _y, _strips.size());
      for (auto x : _strips)
      {
        printf("\t %d %f %f \n", x.strip(), x.xpos(), x.ypos());
      }
    }
    bool isAdjacent(TdcStrip &s, float step = 3)
    {

      for (auto x : _strips)
      {
        //if (abs(x.xpos()-s.xpos())<step && abs(x.ypos()-s.ypos())<2)
        if (x.chamber() != s.chamber())
          continue;
        float dta = DTA;
        if (x.dif() != s.dif())
          dta = 3 * DTA;
        if (abs(x.xpos() - s.xpos()) < step && abs((x.t0() + x.t1()) / 2 - (s.t0() + s.t1()) / 2) < dta && abs(x.ypos() - s.ypos()) < DTY)
        {
          return true;
        }
      }
      return false;
    }
    void addStrip(TdcStrip &s)
    {
      _strips.push_back(s);
      this->calcpos();
    }
    void calcpos()
    {
        _x = 0;
        _y = 0, _t0 = 0, _t1 = 0;
	_bp=0;
        for (int i = 0; i < _strips.size() ; i++)
        {
	  _bp|=(1<<_strips[i].strip());
          _x += _strips[i].xpos();
          _y += _strips[i].ypos();
          _t0 += _strips[i].t0();
          _t1 += _strips[i].t1();
        }
        _x /= (_strips.size() );
        _y /= (_strips.size());
        _t0 /= (_strips.size());
        _t1 /= (_strips.size());
        return;
      std::sort(_strips.begin(), _strips.end(), [](const TdcStrip &lhs, const TdcStrip &rhs) {
        return lhs.strip() < rhs.strip();
      });

      if (_strips.size() == 1)
      {
        _x = _strips[0].xpos();
        _y = _strips[0].ypos();
        _t0 = _strips[0].t0();
        _t1 = _strips[0].t1();
      }
      if (_strips.size() == 2)
      {
        _x = (_strips[0].xpos() + _strips[1].xpos()) / 2.;
        _y = (_strips[0].ypos() + _strips[1].ypos()) / 2.;
        _t0 = (_strips[0].t0() + _strips[1].t0()) / 2.;
        _t1 = (_strips[0].t1() + _strips[1].t1()) / 2.;
      }
      if (_strips.size() == 3)
      {
        _x = _strips[1].xpos();
        _y = _strips[1].ypos();
        _t0 = _strips[1].t0();
        _t1 = _strips[1].t1();
      }
      if (_strips.size() == 4)
      {
        _x = (_strips[2].xpos() + _strips[1].xpos()) / 2.;
        _y = (_strips[2].ypos() + _strips[1].ypos()) / 2.;
        _t0 = (_strips[2].t0() + _strips[1].t0()) / 2.;
        _t1 = (_strips[2].t1() + _strips[1].t1()) / 2.;
      }
      if (_strips.size() >= 5 && _strips.size() <= 12)
      {
        _x = 0;
        _y = 0, _t0 = 0, _t1 = 0;
        for (int i = 2; i < _strips.size() - 2; i++)
        {
          _x += _strips[i].xpos();
          _y += _strips[i].ypos();
          _t0 += _strips[i].t0();
          _t1 += _strips[i].t1();
        }
        _x /= (_strips.size() - 4);
        _y /= (_strips.size() - 4);
        _t0 /= (_strips.size() - 4);
        _t1 /= (_strips.size() - 4);
      }

      
    }
    inline double X() { return _x; }
    inline double Y() { return _y; }
    inline double T0() { return _t0; }
    inline double T1() { return _t1; }
    inline double TM() { return (_t0 + _t1) / 2.; }
    uint32_t size() { return _strips.size(); }
    Lmana::TdcStrip &strip(int n) { return _strips[n]; }
    inline uint16_t chamber() const { return (_strips.size() > 0) ? _strips[0].chamber() : 0; }
    inline uint16_t dif() const { return (_strips.size() > 0) ? _strips[0].dif() : 0; }
    inline uint64_t bp(){return _bp;}
  private:
    double _x, _y, _t0, _t1;
    std::vector<Lmana::TdcStrip> _strips;
    uint64_t _bp;
  };

class HR2Pad : public ROOT::Math::XYZPoint
{
public:
  HR2Pad(double x,double y,double z,uint32_t plan,uint32_t dif,uint32_t i,uint32_t j)
  {
    this->SetXYZ(x,y,z);
    valid_=true;_used=false;
    _plan=plan;
    _dif=dif;
    _I=i;
    _J=j;
    
  }
  HR2Pad(uint32_t plan,uint32_t dif,uint32_t i,uint32_t j){_plan=plan; _dif=dif;
    _I=i;
    _J=j;}
  ~HR2Pad(){;}
  inline double const dX(){return 0.3;}
  inline double const dY(){return 0.3;}
  inline uint32_t  plan(){return _plan;}
  inline uint32_t  dif(){return _dif;}
  inline uint32_t  I(){return _I;}
  inline uint32_t  J(){return _J;}

  inline bool isUsed() const {return _used;}
  inline void setUse(bool t ){_used=t;}
  bool operator < (const Lmana::HR2Pad& str) const
  {
    return (Z() < str.Z());
  }
  void setValidity(bool t){valid_=t;}
  bool isValid() const {return valid_;}
  void setPlan(uint32_t p){_plan=p;}
protected:
  bool valid_,_used;
  uint32_t _plan,_dif,_I,_J;
};

  
    class HR2Cluster
  {
  public:
    HR2Cluster() : _x(0), _y(0) { _pads.clear(); }
    void Print()
    {
      printf("X %f Y %f Size %ld \n", _x, _y, _pads.size());
      for (auto x : _pads)
      {
        printf("\t %d %d %d %f %f \n", x.dif(),x.I(),x.J(), x.X(), x.Y());
      }
    }
    bool isAdjacent(Lmana::HR2Pad &s, float step = 3)
    {

      for (auto x : _pads)
      {
        //if (abs(x.xpos()-s.xpos())<step && abs(x.ypos()-s.ypos())<2)
        if (x.dif() != s.dif())
          continue;
        if (abs(x.X() - s.X()) < step && abs(x.Y() - s.Y()) < step) 
        {
          return true;
        }
      }
      return false;
    }
    void addPad(HR2Pad &s)
    {
      _pads.push_back(s);
      this->calcpos();
    }
    void calcpos()
    {
        _x = 0;
        _y = 0;

        for (int i = 0; i < _pads.size() ; i++)
        {
          _x += _pads[i].X();
          _y += _pads[i].Y();
        }
        _x /= (_pads.size() );
        _y /= (_pads.size());
        return;
      
    }
    inline double X() { return _x; }
    inline double Y() { return _y; }
    inline double Z() { return  (_pads.size() > 0) ? _pads[0].Z() : 0;  }
    uint32_t size() { return _pads.size(); }
    Lmana::HR2Pad &pad(int n) { return _pads[n]; }
    inline uint32_t chamber()  { return (_pads.size() > 0) ? _pads[0].plan() : 0; }
    inline uint32_t dif()  { return (_pads.size() > 0) ? _pads[0].dif() : 0; }

  private:
    double _x, _y;
    std::vector<Lmana::HR2Pad> _pads;
  };

}; // namespace Lmana

class rbEvent
{
public:
  rbEvent() : _difType(false) { ; }
  inline void init(uint32_t run, uint32_t event, uint32_t gtc, uint64_t abcid)
  {
    _run = run;
    _event = event;
    _gtc = gtc;
    _abcid = abcid;
    memset(_frameCount, 0, MAXDIF * sizeof(uint32_t));

    if (_difType)
    {
      for (uint32_t i = 0; i < thePMRPtrList_.size(); i++)
      {

        delete thePMRPtrList_[i];
      }
      thePMRPtrList_.clear();
    }
  }
  inline void setDifType(bool t) { _difType = t; }
  inline bool isDifType() { return _difType; }
  inline uint32_t frameCount(uint32_t d) { return _frameCount[d]; }
  inline void setFrameCount(uint32_t d, uint32_t n) { _frameCount[d] = n; }
  inline uint8_t *frameBuffer() { return _fBuf; }
  inline uint32_t run() { return _run; }
  inline uint32_t event() { return _event; }
  inline uint32_t gtc() { return _gtc; }
  inline uint64_t abcid() { return _abcid; }
  inline uint32_t iPtr(uint32_t d, uint32_t f) { return (d * MAXFRAME + f) * FSIZE; }
  inline uint32_t header(uint32_t iptr) { return _fBuf[iptr]; }
  inline uint32_t bcid(uint32_t iptr) { return PMRUnpacker::getFrameBCID(&_fBuf[iptr]); }
  inline bool padLevel(uint32_t iptr, uint32_t ipad, uint32_t ilevel) { return PMRUnpacker::getFrameLevel(&_fBuf[iptr], ipad, ilevel); }
  inline bool pad0(uint32_t iptr, uint32_t ipad) { return padLevel(iptr, ipad, 0); }
  inline bool pad1(uint32_t iptr, uint32_t ipad) { return padLevel(iptr, ipad, 1); }
  inline std::map<uint32_t, std::vector<std::pair<sdhcal::PMRPtr *, uint32_t>>> &tFrame() { return _tframe; }
  inline std::map<uint32_t, std::bitset<64>> &tCount() { return _tcount; }
  inline std::vector<sdhcal::PMRPtr *> &difList() { return thePMRPtrList_; }
  inline std::vector<lydaq::TdcChannel> &tdcChannels() { return _vAll; }
  inline void setNoiseRun(bool n){_noiseRun=n;}
  inline bool noiseRun(){return _noiseRun;}
private:
  uint32_t _run, _event, _gtc;
  uint64_t _abcid;
  uint32_t _frameCount[MAXDIF];
  uint8_t _fBuf[MAXDIF * MAXFRAME * FSIZE];
  bool _difType;
  std::vector<sdhcal::PMRPtr *> thePMRPtrList_;
  std::map<uint32_t, std::vector<std::pair<sdhcal::PMRPtr *, uint32_t>>> _tframe;
  std::map<uint32_t, std::bitset<64>> _tcount;

  std::vector<lydaq::TdcChannel> _vAll;
  bool _noiseRun;
};

class rbProcessor
{
public:
  virtual void init(uint32_t run = 0) = 0;
  virtual void end(uint32_t run = 0) = 0;
  virtual void processEvent(rbEvent *e) = 0;
  virtual void processRunHeader(std::vector<uint32_t> header) = 0;
  virtual void loadParameters(Json::Value params) = 0;
};
#endif
