#ifndef _tricotreader_h
#define _tricotreader_h

#include <stdint.h>
#include <stdlib.h>
#include "rbevent.hh"
#include "DCHistogramHandler.hh"
#include <Math/PositionVector3D.h>
#include <Math/Point3Dfwd.h>
#include <Math/Vector3Dfwd.h>
#include <Math/DisplacementVector3D.h>



#include <vector>
#include <map>
#include <string>
#include <bitset>
#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include "recoTrack.hh"
#include "TTree.h"
#include "TFile.h"

//#define TSTEP 0.432
#define TSTEP 0.4988
#define XSTEP 0.43197347
#define T30 0.5773502691896257

//1.732050808
typedef struct {
  double xm[3];
  double lambda[3];
  double l0[3];
  double l1[3];
  double l2[3];
  double fp;
  double lp;
  double fx,lx;
  double fy,ly;
} ShowerParams;

class tricotreader : public rbProcessor
  {
  public:
    tricotreader();
    virtual void info();
    virtual void init(uint32_t run=0);
    virtual void end(uint32_t run=0);
    virtual  void processEvent(rbEvent* e);
    virtual  void processRunHeader(std::vector<uint32_t> header);
    virtual void loadParameters(Json::Value params);
    void scurveAnalysis(rbEvent* e);
    void fillTimeMap(rbEvent *e);
    void buildPlaneHits(rbEvent *e, std::vector<uint32_t> &hits);

  private:
    uint32_t _run,_event,_totalSize,_gtc;
    int32_t _fdOut;
    bool _started,_dummy;
    DCHistogramHandler* _rh;
    std::map<uint32_t,std::vector<uint32_t> > _timeMap;
    uint32_t _maxTime;
    std::vector<recoPoint> _vPoints;
    Json::Value _jparams;
    std::map<uint32_t,Json::Value> _plinfo;
    std::bitset<16> _hplanes;
    uint64_t _bsplanes;

  };

#endif
