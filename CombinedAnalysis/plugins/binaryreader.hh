#ifndef _binaryreader_h
#define _binaryreader_h

#include <stdint.h>
#include <stdlib.h>
#include "rbevent.hh"
#include "DCHistogramHandler.hh"

#include "jsonGeo.hh"
#include "recoTrack.hh"
#include <Math/PositionVector3D.h>
#include <Math/Point3Dfwd.h>
#include <Math/Vector3Dfwd.h>
#include <Math/DisplacementVector3D.h>


#include <vector>
#include <map>
#include <string>
#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
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
    void processCoincidence(rbEvent* e,uint32_t ibc);
    int TPrincipalComponents(double result[21],float zmin,float zmax);
    bool stripStudy(std::vector<lydaq::TdcChannel>& vChannel,std::string subdir);
  private:
    uint32_t _run,_event,_totalSize;
    int32_t _fdOut;
    bool _started,_dummy;
    DCHistogramHandler* _rh;
    Json::Value _geoRoot;
    jsonGeo* _geo;
    std::vector<recoPoint> _vPoints;


    recoTrack top_tk;
    std::vector<Lmana::TdcStrip> _strips;
    std::vector<Lmana::TdcCluster> _clusters;

    ROOT::Math::XYZPoint _pex;
  };

#endif
