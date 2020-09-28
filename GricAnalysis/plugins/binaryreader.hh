#ifndef _binaryreader_h
#define _binaryreader_h

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

class binaryreader : public rbProcessor
  {
  public:
    binaryreader();
    virtual void init(uint32_t run=0);
    virtual void end(uint32_t run=0);
    virtual  void processEvent(rbEvent* e);
    virtual  void processRunHeader(std::vector<uint32_t> header);
    virtual void loadParameters(Json::Value params);
    void fillTimeMap(rbEvent* e);
    void buildTracks();
    void fillTracks();
    void kickSearch();
    void drawHits();
    void buildPosition(rbEvent* e,uint32_t plane,uint32_t ddmax=0,bool all=false);
    void buildPlaneHits(rbEvent* e,uint32_t plane,std::vector<uint32_t>& hits);
    int TPrincipalComponents(double result[21],float zmin,float zmax);
    void createTrees(std::string s);
    void closeTrees();
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
    // Tracking
    recoTrack top_tk;
    recoTrack bot_tk;
    recoTrack a_tk;
    
    recoTrack vstop[10000];
    uint32_t nstop,lastprocessed;
    // TTree
    // Event gtc timestamp
    uint32_t _bxdif;
    // Bit set 128 strip *3 directions *8 plans
    uint64_t _bs[48];
    // Tracks
    uint8_t _t_h,_b_h,_a_h;
    double _t_x[3],_t_v[3],_b_x[3],_b_v[3],_a_x[3],_a_v[3];
    double _t_c2,_b_c2,_a_c2;
    double _cos_th,_th,_xcross,_ycross,_zcross,_dist;
    double _rd3,_probd3;
    TTree* tEvents_;
    TFile* treeFile_;

  };

#endif
