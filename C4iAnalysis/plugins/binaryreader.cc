#include "binaryreader.hh"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <fcntl.h>
#include <iostream>
#include <sstream>

#include <TCanvas.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TFitter.h>
#include <TF1.h>
#include <TPluginManager.h>
#include <stdint.h>
#include <math.h>
#include "TPolyLine3D.h"
#include "TVirtualPad.h"
#include "recoTrack.hh"
#include "HoughLocal.hh"
#include "TPrincipal.h"

static TCanvas *TCPlot = NULL;
static TCanvas *TCHits = NULL;
static TCanvas *TCShower = NULL;
static TCanvas *TCEdge = NULL;
static TCanvas *TCHT = NULL;
static TCanvas *TCCluster = NULL;
static float z[16] = {0, 54., 0, 0, 76., 10., 0, 0, 32, 0, 0, 0, 0, 0, 0, 0};
static int32_t togric[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//using namespace zdaq;

// TTree info

typedef struct
{
  uint32_t run;
  uint32_t event;
  uint32_t gtc;
  uint32_t idx;
  uint64_t bcid;
  uint64_t bsplan;
  float t_ax, t_bx, t_ay, t_by, t_chi2, t_pchi2;
  float b_ax, b_bx, b_ay, b_by, b_chi2, b_pchi2;
  float ctheta, dist, xcross, ycross, zcross;
  float ax, bx, ay, by, chi2, pchi2;
} event_t;

event_t theEvent;

binaryreader::binaryreader() : _run(0), _started(false), _fdOut(-1), _totalSize(0), _event(0),tEvents_(NULL) {}
void binaryreader::init(uint32_t run)
{
  _run = run;
  _event = 0;
  _started = true;
  _rh = DCHistogramHandler::instance();
  nstop = 0;
}
void binaryreader::loadParameters(Json::Value params)
{
  _jparams = params;
  if (params.isMember("general"))
    std::cout << "DIRECTORY " << params["general"]["directory"].asString() << std::endl;
  //getchar();
  if (params.isMember("setup"))
    {
      const Json::Value &b = params["setup"];

      for (Json::ValueConstIterator it = b.begin(); it != b.end(); ++it)
	{
	  const Json::Value &ch = *it;
	  uint32_t pl = ch["id"].asUInt();
	  uint32_t plane = ch["num"].asUInt();
	  std::pair<uint32_t, Json::Value> p(pl, ch);
	  z[plane] = ch["z"].asFloat();
	  togric[plane]=pl;
	  _plinfo.insert(p);
	  printf("%s GRIC %d Plane  %d is at %f cm from floor \n", _plinfo[pl]["name"].asString().c_str(), pl, plane, z[plane]);
	}
      //getchar();
    }
}
void binaryreader::end(uint32_t run)
{
  //this->closeTrees();
  _started = false;
}

void binaryreader::processRunHeader(std::vector<uint32_t> header)
{
}
void binaryreader::fillTimeMap(rbEvent *e)
{
  _timeMap.clear();
  std::stringstream sraw;
  sraw << "/gric/";
  TH1 *hfc = _rh->GetTH1(sraw.str() + "FrameCount");
  TH1 *hftm = _rh->GetTH1(sraw.str() + "MaxTime");


  if (hfc == NULL)
    {

      hfc = _rh->BookTH1(sraw.str() + "FrameCount", 255, 0., 255.);
      hftm = _rh->BookTH1(sraw.str() + "MaxTime", 65536, 0., 2.);
    }

  // Find Maximum time and ramFull
  _maxTime = 0;
  bool _ramFull = false;
  for (int id = 0; id < MAXDIF; id++)
    {
      hfc->Fill(id * 1., e->frameCount(id));
      _ramFull = _ramFull || (e->frameCount(id) > 126);
      for (int j = 0; j < e->frameCount(id); j++)
	{
	  uint32_t idx = e->iPtr(id, j);

	  if (e->bcid(idx) > _maxTime)
	    _maxTime = e->bcid(idx);
	}
    }
  hftm->Fill(_maxTime*2E-7);
  // Fill the Map
   int32_t bcidmin=1;
  for (int id = 0; id < MAXDIF; id++)
    if (e->frameCount(id))
      {
	 #define RAWHIST
  #ifdef RAWHIST
  std::stringstream sraw1;
  sraw1 << "/gric/ASIC" << std::hex << id << std::dec << "/";
  
  TH1 *hp1 = _rh->GetTH1(sraw1.str() + "Pad1");
  TH1 *hft = _rh->GetTH1(sraw1.str() + "FrameTime");
  TH1 *hfc = _rh->GetTH1(sraw1.str() + "FrameCount");
  if (hp1 == NULL)
    {
      hp1 = _rh->BookTH1(sraw1.str() + "Pad1", 64, 0., 64.);
      hft = _rh->BookTH1(sraw1.str() + "FrameTime", 65536, 0., 2.);
      hfc = _rh->BookTH1(sraw1.str() + "FrameCount", 65, 0., 65.);
    }


  #endif
 
	uint16_t igplane = (id >> 4) & 0xF;
	uint16_t plane = _plinfo[igplane]["num"].asUInt();
			
	for (int j = 0; j < e->frameCount(id); j++)
	  {
	    uint32_t idx = e->iPtr(id, j);
	    int32_t dd = _maxTime - e->bcid(idx);
	    // Cleaning 
	    // Remove frame at beginning of windows
	    if (e->bcid(idx) < bcidmin)
	      continue;
	    //if (e->frameCount(id)==127 && (dd<10)) continue;
	    // Remove full ASIC fired
	    std::bitset<64> bs;
	    bs.reset();
	    for (int k = 0; k < 64; k++)
	      if (e->pad0(idx, k) || e->pad1(idx, k))
		{
		  bs.set(k);
		  #ifdef RAWHIST
		  hp1->Fill(k * 1.);
		  #endif
		}
	    #ifdef RAWHIST
	    hfc->Fill(bs.count() * 1.);
	    hft->Fill((_maxTime - e->bcid(idx)) * 2E-7);

	    #endif
	    if (bs.count() > 60)
	      continue;
				
				
	    bool found = false;
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itm = _timeMap.find(e->bcid(idx));
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itmm = _timeMap.find(e->bcid(idx) - 1);
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itmp = _timeMap.find(e->bcid(idx) + 1);
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itmmm = _timeMap.find(e->bcid(idx) - 2);
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itmpp = _timeMap.find(e->bcid(idx) + 2);
	    found = (itm != _timeMap.end()) || (itmp != _timeMap.end()) || (itmm != _timeMap.end()) || (itmpp != _timeMap.end()) || (itmmm != _timeMap.end());
	    found =(itm!=_timeMap.end())||(itmp!=_timeMap.end())||(itmm!=_timeMap.end());
	    if (!found)
	      {
		std::vector<uint32_t> v;
		v.push_back(idx);
		std::pair<uint32_t, std::vector<uint32_t>> p(e->bcid(idx), v);
		_timeMap.insert(p);
	      }
	    else
	      {
		if (itm != _timeMap.end())
		  itm->second.push_back(idx);
		if (itmm != _timeMap.end())
		  itmm->second.push_back(idx);
		if (itmp != _timeMap.end())
		  itmp->second.push_back(idx);
		/*
		if (itmmm != _timeMap.end())
		  itmmm->second.push_back(idx);
		if (itmpp != _timeMap.end())
		  itmpp->second.push_back(idx);
		*/
	      }

	  }
      }

  // Now select time slice with at least NPLANESMIN hit
  int32_t nplanesmin = _jparams["general"]["minplan"].asUInt();
  std::bitset<16> planes;
  planes.reset();
  uint8_t u[16], v[16], w[16];
  for (auto it=_timeMap.begin();it!=_timeMap.end();)
    {
      auto x=(*it);
      _bxdif = x.first;
      int32_t dd = _maxTime - x.first;
      // On veut 2 directions sur nplanesmin
      bool candidate=false;
      if (x.second.size() > nplanesmin * 2 && dd>20)
	{
	      
	  planes.reset();
	  memset(u, 0, 16);
	  memset(v, 0, 16);
	  memset(w, 0, 16);
	  for (auto y : x.second)
	    {
	      uint32_t ig = y / MAXFRAME / FSIZE;
	      uint32_t igplane = (ig >> 4) & 0xF;
	      uint32_t plane = _plinfo[igplane]["num"].asUInt();
	      u[plane] = u[plane] || ((ig & 0xF) == 2 || (ig & 0xF) == 4);
	      v[plane] = v[plane] || ((ig & 0xF) == 3 || (ig & 0xF) == 5);
	      w[plane] = w[plane] || ((ig & 0xF) == 6 || (ig & 0xF) == 7);
	      planes.set(plane);
	      //std::cout<<std::hex<<ig<<std::dec<<" IGPL"<<" "<<igplane<<"("<<plane<<")"<<std::endl;
	    }
	  planes.reset();
	  for (int ip = 0; ip < 16; ip++)
	    {
	      if (u[ip] && v[ip])
		planes.set(ip);
	      if (u[ip] && w[ip])
		planes.set(ip);
	      if (w[ip] && v[ip])
		planes.set(ip);
	    }
	  uint32_t nplanes = planes.count();
	  candidate= (nplanes >= nplanesmin);
	}
      if (candidate)
	{
	  it++;
	}
      else
	{
	  it->second.clear();
	  _timeMap.erase(it++);
	}
    }
  // Now Build the planes hits
  if (_timeMap.size()==0) return;
  bool display = _jparams["general"]["display"].asUInt() == 1;
  for (auto x:_timeMap)
    {
      // Build Points list
      _bxdif=x.first;
      _vPoints.clear();
      _hplanes.reset();
      for (int i = 1; i <= 8; i++)
	buildPlaneHits(e, i, x.second);
      if (_hplanes.count()<nplanesmin) continue;
      this->fillTracks();
      uint64_t ltopbs=(_hplanes.to_ulong()&30);
      uint64_t lbotbs=(_hplanes.to_ulong()&480);
      _bsplanes= _hplanes.to_ulong();
      std::bitset<16> topbs(ltopbs);
      std::bitset<16> botbs(lbotbs);
      if (topbs.count()<2) continue;
      if (botbs.count()<2) continue;
      //std::cout<<_run<<" "<<_event<<" Candidate " <<x.first<<" Pt "<<_vPoints.size()<<" pattern "<<_hplanes<<" "<<topbs <<" "<<botbs<<std::endl;
      this->buildTracks();
      this->kickSearch();
      if (display) this->drawHits();	  
    }
  //getchar();

}

void binaryreader::buildTracks()
{
  // Build a top segments
  top_tk.clear();
  bot_tk.clear();
  float px2cut = 0.01, thcut = 1.73, dcut = 7.;
  ShowerParams isha;
  double ax = 0, bx = 0, ay = 0, by = 0;
  float zmin = z[4] - 1, zmax = z[1] + 1;
  int ier = TPrincipalComponents((double *)&isha, zmin, zmax);
  double *sx = isha.xm;
  double *sv = isha.l2;
  double z0 = sx[2];
  double x0 = sx[0];
  double y0 = sx[1];
  double x1 = sx[0] + sv[0];
  double y1 = sx[1] + sv[1];
  double z1 = sx[2] + sv[2];
  //double ax,ay,bx,by;
  if (ier == 0)
    {
      ax = (x1 - x0) / (z1 - z0);
      bx = x1 - ax * z1;
      ay = (y1 - y0) / (z1 - z0);
      by = y1 - ay * z1;
    }
  else
    return;
  

  _t_h = 0;
  _b_h = 0;
  _a_h = 0;
  top_tk.setDir(ax, ay, 1.);
  top_tk.setOrig(bx, by, 0);
  //printf(" ax %f bx %f ay %f by %f \n",ax,ay,bx,by);
  for (auto x = _vPoints.begin(); x != _vPoints.end(); x++)
    {
      recoPoint &p = (*x);
      if (p.Z() < zmin)
	continue;
      if (p.Z() > zmax)
	continue;
      
      if (top_tk.distance(&p) < dcut || ier != 0)
	{
	  //	  std::cout<<x->X()<<" "<<x->Y()<<" "<<x->Z()<<" TOPP_TK"<<x->plan()<<std::endl;

	  top_tk.addPoint(&p);
	  _t_h |= (1 << p.plan() - 1);
	}
    }
  // Ask at least 2 points
  if (top_tk.size() < 2)
    return;
  top_tk.regression();
  top_tk.calculateChi2();
  
  if (top_tk.pchi2() < px2cut &&top_tk.size()>2)
    {
      //fprintf(stderr,"top Chi2 cut \n");
      top_tk.clear();
      return;
    }
  if (abs(top_tk.dir().X()) > thcut)
    {
      //fprintf(stderr,"top X thcut \n");
      top_tk.clear();
      return;
    }


  if (abs(top_tk.dir().Y()) > thcut)
    {
      //fprintf(stderr,"top  Y thcut \n");
      top_tk.clear();
      return;
    }


  
  _t_c2 = top_tk.pchi2();
  _t_x[0] = top_tk.orig().X();
  _t_x[1] = top_tk.orig().Y();
  _t_x[2] = top_tk.orig().Z();
  _t_v[0] = top_tk.dir().X();
  _t_v[1] = top_tk.dir().Y();
  _t_v[2] = top_tk.dir().Z();

  // Bottom track
  zmin = z[8] - 1, zmax = z[5] + 1;
  memset(&isha, 0, sizeof(ShowerParams));
  ier = TPrincipalComponents((double *)&isha, zmin, zmax);
  sx = isha.xm;
  sv = isha.l2;
  z0 = sx[2];
  x0 = sx[0];
  y0 = sx[1];
  x1 = sx[0] + sv[0];
  y1 = sx[1] + sv[1];
  z1 = sx[2] + sv[2];
  //double ax,ay,bx,by;
  if (ier == 0)
    {
      ax = (x1 - x0) / (z1 - z0);
      bx = x1 - ax * z1;
      ay = (y1 - y0) / (z1 - z0);
      by = y1 - ay * z1;
    }
  else
    {
      fprintf(stderr,"bot No shower \n");
      return;
    }
  bot_tk.clear();
  bot_tk.setDir(ax, ay, 1.);
  bot_tk.setOrig(bx, by, 0);

  for (auto x = _vPoints.begin(); x != _vPoints.end(); x++)
    {
      recoPoint &p = (*x);
      if (p.Z() < zmin)
	continue;
      if (p.Z() > zmax)
	continue;
      
      if (bot_tk.distance(&p) < dcut || ier != 0)
	{
	  //std::cout<<x->X()<<" "<<x->Y()<<" "<<x->Z()<<" BOT_TK"<<x->plan()<<std::endl;
	  bot_tk.addPoint(&p);
	  _b_h |= (1 << p.plan() - 1);
	}
    }
  if (bot_tk.size() >=2)
    {
      bot_tk.regression();
      bot_tk.calculateChi2();
      _b_c2 = bot_tk.pchi2();
      _b_x[0] = bot_tk.orig().X();
      _b_x[1] = bot_tk.orig().Y();
      _b_x[2] = bot_tk.orig().Z();
      _b_v[0] = bot_tk.dir().X();
      _b_v[1] = bot_tk.dir().Y();
      _b_v[2] = bot_tk.dir().Z();
      //this->fillTracks();

    }
  else
    {
      fprintf(stderr,"bot Size too small \n");
      bot_tk.clear();
      return;
    }
  return;

}
void binaryreader::fillTracks()
{
  a_tk.clear();
  float dcut=7.0;
  ShowerParams isha;
  double ax = 0, bx = 0, ay = 0, by = 0;
  float zmin = z[8] - 1, zmax = z[1] + 1;
  int ier = TPrincipalComponents((double *)&isha, zmin, zmax);
  double *sx = isha.xm;
  double *sv = isha.l2;
  double z0 = sx[2];
  double x0 = sx[0];
  double y0 = sx[1];
  double x1 = sx[0] + sv[0];
  double y1 = sx[1] + sv[1];
  double z1 = sx[2] + sv[2];
  //double ax,ay,bx,by;
  if (ier == 0)
    {
      ax = (x1 - x0) / (z1 - z0);
      bx = x1 - ax * z1;
      ay = (y1 - y0) / (z1 - z0);
      by = y1 - ay * z1;
    }
  else
    return;
  

  _a_h = 0;
  a_tk.setDir(ax, ay, 1.);
  a_tk.setOrig(bx, by, 0);
  //printf("A_TK  ax %f bx %f ay %f by %f \n",ax,ay,bx,by);
  for (auto x = _vPoints.begin(); x != _vPoints.end(); x++)
    {
      recoPoint &p = (*x);
      if (p.Z() < zmin)
	continue;
      if (p.Z() > zmax)
	continue;
      
      if (a_tk.distance(&p) < dcut || ier != 0)
	{
	  //std::cout<<x->X()<<" "<<x->Y()<<" "<<x->Z()<<" A_TK"<<x->plan()<<std::endl;
	  //fflush(stdout);

	  a_tk.addPoint(&p);
	  _a_h |= (1 << p.plan() - 1);
	}
    }
  if (a_tk.size()<=2)
    {
      a_tk.clear();
      return;
    }
  //std::cout<<a_tk<<std::endl<<"============================="<<std::endl;;
  if ((_a_h & 1)&& (_a_h &(1<<7)))
    {
      recoTrack ta;ta.clear();
      for (auto x = a_tk.points().begin(); x != a_tk.points().end(); x++)
	{
	  recoPoint *p = (*x);
	  if (p->plan() !=1 && p->plan()!=8)
	    continue;
	  //std::cout<<p->X()<<" "<<p->Y()<<" "<<p->Z()<<" TP_TK"<<p->plan()<<std::endl;
	  //fflush(stdout);

	  ta.addPoint(p);
	}
      ta.regression();
      for (int ip = 1; ip <= 8; ip++)
	{
	  if (z[ip] == 0)
	    continue;
	  ROOT::Math::XYZPoint pe = ta.extrapolate(z[ip]);
	  for (auto x = _vPoints.begin(); x != _vPoints.end(); x++)
	    {
	      if (x->plan() != ip)
		continue;
	      std::stringstream splane;
	      splane << "/gric/ALG" << ip << "_G"<<togric[ip]<<"/";

	      TH1 *hdistx = _rh->GetTH1(splane.str() + "DX");
	      TH1 *hdisty = _rh->GetTH1(splane.str() + "DY");
	      if (hdistx == NULL)
		{
		  hdistx = _rh->BookTH1(splane.str() + "DX", 100, -10., 10.);
		  hdisty = _rh->BookTH1(splane.str() + "DY", 100, -10., 10.);
		}

	      ROOT::Math::XYZPoint &p = (*x);
	      ROOT::Math::XYZVector d = p - pe;
	      float dx = d.X();
	      float dy = d.Y();
	  //std::cout<<p1.X()<<" "<<p1.Y()<<" "<<p1.Z()<<std::endl;
	      double dist = sqrt(d.Mag2());
	  //std::cout<<dist<<" "<<hdistx->GetEntries()<<std::endl;
	      hdistx->Fill(dx);
	      //std::cout<<dx<<" "<<dy<<" "<<hdistx->GetEntries()<<std::endl;
	      hdisty->Fill(dy);
	}

	}

    }
  for (int ip = 1; ip <= 8; ip++)
    {
      if (z[ip] == 0)
	continue;
      uint16_t bs=0;
      recoTrack tp;
      tp.clear();
      // for (auto x = a_tk.points().begin(); x != a_tk.points().end(); x++)
      // 	{
      // 	  recoPoint *p = (*x);
      // 	  //	  std::cout<<p->X()<<" "<<p->Y()<<" "<<p->Z()<<" A_TK"<<p->plan()<<std::endl;
      // 	}
      for (auto x = a_tk.points().begin(); x != a_tk.points().end(); x++)
	{
	  recoPoint *p = (*x);
	  if (p->plan() == ip)
	    continue;
	  //std::cout<<p->X()<<" "<<p->Y()<<" "<<p->Z()<<" TP_TK"<<p->plan()<<std::endl;
	  //fflush(stdout);

	  tp.addPoint(p);
	}
      if (tp.size() < 3)
	continue;
      tp.regression();
      tp.calculateChi2();
      

      if (tp.pchi2() < 0.01)
	continue;
      /*
      if (ip == 1 && !(tp.planUsed(2) && tp.planUsed(3)))
	continue;
      if (ip == 1 && !(tp.planUsed(2) && tp.planUsed(8)))
	continue;

      if (ip == 8 && !(tp.planUsed(6) && tp.planUsed(7)))
	continue;
      if (ip == 8 && !(tp.planUsed(6) && tp.planUsed(1)))
	continue;
      if (ip > 1 && ip < 8 && !(tp.planUsed(ip - 1) && tp.planUsed(ip + 1)))
	continue;
      if (ip > 1 && ip < 8 && !(tp.planUsed(1) && tp.planUsed(8)))
	continue;
      */
      if (tp.size()<4)
	continue;
      ROOT::Math::XYZPoint pe = tp.extrapolate(z[ip]);

      if (pe.X() < -27.648 + 3)
	continue;
      if (pe.X() > 27.648 - 3)
	continue;
      if (pe.X() < 0)
	{
	  if (pe.Y() < -1. * T30 * pe.X() + 3)
	    continue;
	  if (pe.Y() > T30 * pe.X() + 55.296 - 3)
	    continue;
	}
      else
	{
	  if (pe.Y() < T30 * pe.X() + 3)
	    continue;
	  if (pe.Y() > -1. * T30 * pe.X() + 55.296 - 3)
	    continue;
	}
      //std::cout<<"PLAN "<<ip<<" "<<tp<<std::endl;
      //getchar();


      std::stringstream splane;
      splane << "/gric/EFF" << ip << "_G"<<togric[ip]<<"/";
      TH2 *hext = _rh->GetTH2(splane.str() + "XYext");
      TH2 *hfop = _rh->GetTH2(splane.str() + "XYfound");

      TH1 *hdistx = _rh->GetTH1(splane.str() + "DX");
      TH1 *hdisty = _rh->GetTH1(splane.str() + "DY");
      TH1 *hax = _rh->GetTH1(splane.str() + "AX");
      TH1 *hay = _rh->GetTH1(splane.str() + "AY");
      TH1 *hpc2 = _rh->GetTH1(splane.str() + "PChi2");
      if (hdistx == NULL)
	{
	  hpc2 = _rh->BookTH1(splane.str() + "PChi2", 1000, 0., 1.);
	  hax = _rh->BookTH1(splane.str() + "AX", 1000, -1., 1.);
	  hay = _rh->BookTH1(splane.str() + "AY", 1000, -1., 1.);
	  hdistx = _rh->BookTH1(splane.str() + "DX", 100, -10., 10.);
	  hdisty = _rh->BookTH1(splane.str() + "DY", 100, -10., 10.);
	  hext = _rh->BookTH2(splane.str() + "XYext", 32, -30., 30., 32, 0., 60.);
	  hfop = _rh->BookTH2(splane.str() + "XYfound", 32, -30., 30., 32, 0., 60.);
	}
      hpc2->Fill(tp.pchi2());
      hax->Fill(tp.dir().X());
      hay->Fill(tp.dir().Y());

      //if ( abs(tp.dir().X())>thcut || abs(tp.dir().Y())>thcut) continue;
      hext->Fill(pe.X(), pe.Y());

      for (auto x = _vPoints.begin(); x != _vPoints.end(); x++)
	{
	  if (x->plan() != ip)
	    continue;
	  ROOT::Math::XYZPoint &p = (*x);
	  ROOT::Math::XYZVector d = p - pe;
	  float dx = d.X();
	  float dy = d.Y();
	  //std::cout<<p1.X()<<" "<<p1.Y()<<" "<<p1.Z()<<std::endl;
	  double dist = sqrt(d.Mag2());
	  //std::cout<<dist<<" "<<hdistx->GetEntries()<<std::endl;
	  hdistx->Fill(dx);
	  //std::cout<<dx<<" "<<dy<<" "<<hdistx->GetEntries()<<std::endl;
	  hdisty->Fill(dy);
	  if (dist < 5.)
	    {
	      hfop->Fill(pe.X(), pe.Y());
	      break;
	    }
	}
    }

}
void binaryreader::kickSearch()
{
  if (tEvents_==NULL)
    {
      std::stringstream ss;
      ss<<"/tmp/tree"<<_run<<"_"<<_gtc<<".root";
      this->createTrees(ss.str());			    
    }
  if (top_tk.size()<2) return;
  if (bot_tk.size()<2) return;
  float thetacut = _jparams["general"]["thetacut"].asFloat();
  float pcut = _jparams["general"]["pcut"].asFloat();

  float cpacut = _jparams["general"]["cpacut"].asFloat();

  double sc = top_tk.dir().Dot(bot_tk.dir()) / sqrt(top_tk.dir().Mag2()) / sqrt(bot_tk.dir().Mag2());
  double th = acos(sc) * 180 / M_PI;
  //std::cout << "ICI COSTH " << sc << " DEG " << th << std::endl;

  _cos_th = sc;
  _th = th;
  if (_cos_th < 0.5)
    return;

  double adist;
  ROOT::Math::XYZPoint p1, p2;
  top_tk.cap(bot_tk, adist, p1, p2);

  //std::cout << p1.X() << ":" << p1.Y() << ":" << p1.Z() << std::endl;
  //std::cout << p2.X() << ":" << p2.Y() << ":" << p2.Z() << std::endl;
  float zcross = (p1.Z() + p2.Z()) / 2.;
  //printf("%x %x %x \n",hcpa,hcpai,hcpao);
  //std::cout << "LA DISTANCE " << adist << " Z CROSS " << zcross << std::endl;
  if (zcross < z[5] + 1 || zcross > z[4] - 1)
    return;
  fflush(stdout);

  _xcross = (p1.X() + p2.X()) / 2.;
  _ycross = (p1.Y() + p2.Y()) / 2.;
  _zcross = zcross;
  _dist = adist;

  ROOT::Math::XYZPoint po3 = bot_tk.extrapolate(z[5]);
  ROOT::Math::XYZPoint pi3 = top_tk.extrapolate(z[5]);
  ROOT::Math::XYZVector d3 = po3 - pi3;

  //std::cout<<p1.X()<<" "<<p1.Y()<<" "<<p1.Z()<<std::endl;
  double rd3 = sqrt(d3.Mag2());
  double prob = 1.;

  if (rd3 < 800)
    {
      prob = 1. - erf(rd3 / sqrt(0.5) / 2.);
      if (prob < 1E-20)
	prob = 1E-20;
      //printf("Distance %f %f %f %f \n", rd3, prob, th, thetacut);
    }
  _rd3=rd3;
  _probd3=prob;
  if (tEvents_!=NULL)
    {
      treeFile_->cd();
                
      tEvents_->Fill();
    }


}
void binaryreader::drawHits()
{
  
  TH2 *hzx = _rh->GetTH2("ZX");
  TH2 *hzy = _rh->GetTH2("ZY");
  if (hzx == NULL)
    {
      hzx = _rh->BookTH2("ZX", 100., 0., 400., 120., -40., 40.);
      hzy = _rh->BookTH2("ZY", 100., 0., 400., 180., 0., 90.);
    }
  hzx->Reset();
  hzy->Reset();

  for (auto x:_vPoints)
    {
      hzx->Fill(x.Z(),x.X());
      hzy->Fill(x.Z(),x.Y());
      fprintf(stderr,"%d %d %d %f %f %f \n",_run,_event,_bxdif,x.Z(),x.X(),x.Y());
    }
  if (TCHits == NULL)
    {
      TCHits = new TCanvas("TCHits", "tChits1", 900, 900);
      TCHits->Modified();
      TCHits->Draw();
      TCHits->Divide(1, 2);
    }
  TCHits->cd(1);
  hzx->SetMarkerStyle(25);
  hzx->SetMarkerColor(kRed);
  hzx->Draw("P");
  if (top_tk.size()>1) top_tk.linex()->Draw("SAME");
  if (bot_tk.size()>1) bot_tk.linex()->Draw("SAME");
  TCHits->Modified();
  TCHits->Draw();
  TCHits->Update();
  TCHits->cd(2);
  hzy->SetMarkerStyle(22);
  hzy->SetMarkerColor(kGreen);
  hzy->Draw("P");
  if (top_tk.size()>1) top_tk.liney()->Draw("SAME");
  if (bot_tk.size()>1) bot_tk.liney()->Draw("SAME");
  TCHits->Modified();
  TCHits->Draw();
  TCHits->Update();
  ::usleep(100);
  TCHits->Update();
  getchar();
    

  
}
void binaryreader::processEvent(rbEvent *e)
{
  uint8_t u[16], v[16], w[16];
  if (!_started)
    return;
  if (e->gtc()%100==1)
    printf("BR => %d %d %d %d \n",e->run(),e->event(),e->gtc(),e->seuil());
  _event = e->gtc();
  _run = e->run();
  _gtc=e->gtc();
  
  if (e->seuil()!=0)
    this->scurveAnalysis(e);

  this->fillTimeMap(e);
  return;

  _timeMap.clear();
  _maxTime = 0;
  std::stringstream sraw;
  sraw << "/gric/";
  TH1 *hfc = _rh->GetTH1(sraw.str() + "FrameCount");
  TH1 *hst = _rh->GetTH1("/gric/Stat");
  TH1 *hft = _rh->GetTH1(sraw.str() + "FrameTime");
  TH1 *hftm = _rh->GetTH1(sraw.str() + "MaxTime");
  TH1 *hfts = _rh->GetTH1(sraw.str() + "FrameTimeSelected");

  if (hfc == NULL)
    {
      hst = _rh->BookTH1("/gric/Stat", 255, 0., 255.);
      hfc = _rh->BookTH1(sraw.str() + "FrameCount", 255, 0., 255.);
      hft = _rh->BookTH1(sraw.str() + "FrameTime", 65536, 0., 2.);
      hftm = _rh->BookTH1(sraw.str() + "MaxTime", 65536, 0., 2.);

      hfts = _rh->BookTH1(sraw.str() + "FrameTimeSelected", 10000, 0., 10000.);
    }
  hst->Fill(1.);
  bool ramf = false;
  for (int id = 0; id < MAXDIF; id++)
    {
      ramf = ramf || (e->frameCount(id) > 126);
      for (int j = 0; j < e->frameCount(id); j++)
	{
	  uint32_t idx = e->iPtr(id, j);

	  if (e->bcid(idx) > _maxTime)
	    _maxTime = e->bcid(idx);
	}
    }

  bool coinc = false;

  bool trigger = true;

  int32_t nplanesmin = _jparams["general"]["minplan"].asUInt();
  float pcut = 0.1;
  float thetacut = 1.6;
  float cpacut = 2.0;
  if (_jparams.isMember("general"))
    {
      trigger = _jparams["general"]["trigger"].asUInt() == 1;
      coinc = _jparams["general"]["coinc"].asUInt() == 1;
      pcut = _jparams["general"]["pcut"].asFloat();
      thetacut = _jparams["general"]["thetacut"].asFloat();
      cpacut = _jparams["general"]["cpacut"].asFloat();
    }
  hftm->Fill(_maxTime * 2E-7);
  if (_maxTime < 100)
    return;
  //if (rf) return;
  bool rf = false;
  std::bitset<16> plh;
  plh.reset();
  for (int id = 0; id < MAXDIF; id++)
    if (e->frameCount(id))
      {
	uint16_t igplane = (id >> 4) & 0xF;
	uint16_t plane = _plinfo[igplane]["num"].asUInt();
	plh.set(plane, 1);
	hfc->Fill(id * 1., e->frameCount(id));
	//printf("GRIC %x %d frames \n",id,e->frameCount(id));
	std::stringstream sraw1;
	sraw1 << "/gric/ASIC" << std::hex << id << std::dec << "/";

	TH1 *hp1 = _rh->GetTH1(sraw1.str() + "Pad1");
	TH1 *hftg = _rh->GetTH1(sraw1.str() + "FrameTime");
	TH1 *hfc = _rh->GetTH1(sraw1.str() + "FrameCount");
	if (hp1 == NULL)
	  {
	    hp1 = _rh->BookTH1(sraw1.str() + "Pad1", 64, 0., 64.);
	    hftg = _rh->BookTH1(sraw1.str() + "FrameTime", 65536, 0., 2.);
	    hfc = _rh->BookTH1(sraw1.str() + "FrameCount", 65, 0., 65.);
	  }

	for (int j = 0; j < e->frameCount(id); j++)
	  {
	    uint32_t idx = e->iPtr(id, j);
	    int32_t dd = _maxTime - e->bcid(idx);

	    if (e->bcid(idx) < 50)
	      continue;
	    //if (e->frameCount(id)==127 && (dd<10)) continue;
	    std::bitset<64> bs, bs0, bs1;
	    bs.reset();
	    bs1.reset();
	    bs0.reset();
	    for (int k = 0; k < 64; k++)
	      {
		if (e->pad0(idx, k))
		  bs0.set(k);
		if (e->pad1(idx, k))
		  bs1.set(k);
		if (e->pad0(idx, k) || e->pad1(idx, k))
		  {
		    bs.set(k);
		    //if (dd<100)
		    // hp1->Fill(k*1.);
		  }
	      }
	    hfc->Fill(bs.count() * 1.);
	    if (bs.count() > 60)
	      continue;
	    for (int k = 0; k < 64; k++)
	      {
		if (e->pad0(idx, k))
		  bs0.set(k);
		if (e->pad1(idx, k))
		  bs1.set(k);
		if (e->pad0(idx, k) || e->pad1(idx, k))
		  {
		    bs.set(k);
		    //if (dd<100)
		    hp1->Fill(k * 1.);
		  }
	      }
	    hft->Fill((_maxTime - e->bcid(idx)) * 2E-7);
	    hftg->Fill((_maxTime - e->bcid(idx)) * 2E-7);
	    bool found = false;
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itm = _timeMap.find(e->bcid(idx));
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itmm = _timeMap.find(e->bcid(idx) - 1);
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itmp = _timeMap.find(e->bcid(idx) + 1);
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itmmm = _timeMap.find(e->bcid(idx) - 2);
	    std::map<uint32_t, std::vector<uint32_t>>::iterator itmpp = _timeMap.find(e->bcid(idx) + 2);
	    found = (itm != _timeMap.end()) || (itmp != _timeMap.end()) || (itmm != _timeMap.end()) || (itmpp != _timeMap.end()) || (itmmm != _timeMap.end());
	    //found =(itm!=_timeMap.end())||(itmp!=_timeMap.end())||(itmm!=_timeMap.end());
	    if (!found)
	      {
		std::vector<uint32_t> v;
		v.push_back(idx);
		std::pair<uint32_t, std::vector<uint32_t>> p(e->bcid(idx), v);
		_timeMap.insert(p);
	      }
	    else
	      {
		if (itm != _timeMap.end())
		  itm->second.push_back(idx);
		if (itmm != _timeMap.end())
		  itmm->second.push_back(idx);
		if (itmp != _timeMap.end())
		  itmp->second.push_back(idx);
		if (itmmm != _timeMap.end())
		  itmmm->second.push_back(idx);
		if (itmpp != _timeMap.end())
		  itmpp->second.push_back(idx);
	      }

	    rf = rf | (itm->second.size() >= nplanesmin * 2);
	    /*	    
		    printf("%x frame bcid %.9d %8.5f %s %d\n",id,e->bcid(idx),e->bcid(idx)*2E-7,bs.to_string().c_str(),bs.count());
		    printf("%x frame bcid %.9d %8.5f %s %d\n",id,e->bcid(idx),e->bcid(idx)*2E-7,bs0.to_string().c_str(),bs0.count());
		    printf("%x frame bcid %.9d %8.5f %s %d\n ----------\n",id,e->bcid(idx),e->bcid(idx)*2E-7,bs1.to_string().c_str(),bs1.count());
	    */
	  }
      }
  if (ramf)
    hst->Fill(2.);
  else
    hst->Fill(3.);
  if (ramf && e->gtc() % 100 == 0)
    for (int ip = 1; ip < 9; ip++)
      //if (plh[ip]==1)
      buildPosition(e, ip, 0, true);
  if (rf && trigger)
    return;
  if (!rf && coinc)
    return;
  //std::cout<<"Selected-->"<<_timeMap.size()<<std::endl;
  //  getchar();
  if (rf)
    {
      std::stringstream sres;
      sres.clear();
      std::bitset<16> planes;
      planes.reset();

      if (trigger)
	nplanesmin = 0;
      if (trigger)
	for (int ip = 1; ip < 9; ip++)
	  if (plh[ip] == 1)
	    buildPosition(e, ip, 6, false);

      //buildPosition(e,7,6,false);
      for (auto x : _timeMap)
	{
	  _bxdif = x.first;
	  // On veut 2 directions sur nplanesmin
	  if (x.second.size() < nplanesmin * 2)
	    continue;
	  int32_t dd = _maxTime - x.first;
	  if (dd > 20 and trigger)
	    continue;
	  if (dd < 20 && coinc)
	    continue;
	  sres.str("");
	  sres << "\t" << x.first << " : " << x.second.size() << " => " << _maxTime << " " << dd << std::hex;
	  planes.reset();
	  memset(u, 0, 16);
	  memset(v, 0, 16);
	  memset(w, 0, 16);
	  for (auto y : x.second)
	    {
	      uint32_t ig = y / MAXFRAME / FSIZE;
	      uint32_t igplane = (ig >> 4) & 0xF;
	      uint32_t plane = _plinfo[igplane]["num"].asUInt();
	      u[plane] = u[plane] || ((ig & 0xF) == 2 || (ig & 0xF) == 4);
	      v[plane] = v[plane] || ((ig & 0xF) == 3 || (ig & 0xF) == 5);
	      w[plane] = w[plane] || ((ig & 0xF) == 6 || (ig & 0xF) == 7);
	      planes.set(plane);
	      //std::cout<<std::hex<<ig<<std::dec<<" IGPL"<<" "<<igplane<<"("<<plane<<")"<<std::endl;
	    }
	  planes.reset();
	  for (int ip = 0; ip < 16; ip++)
	    {
	      if (u[ip] && v[ip])
		planes.set(ip);
	      if (u[ip] && w[ip])
		planes.set(ip);
	      if (w[ip] && v[ip])
		planes.set(ip);
	    }
	  uint32_t nplanes = planes.count();
	  if (nplanes < nplanesmin)
	    continue;
	  sres << std::dec << " Plans " << x.second.size() << "->" << nplanes << ": " << planes << std::endl;
	  // On veut 2 des 4 premiers plans touches
	  bool goodmu = (planes[1] && planes[2] && planes[3] && planes[4]);
	  goodmu = goodmu || (planes[1] && planes[2] && planes[3]);
	  goodmu = goodmu || (planes[1] && planes[2] && planes[4]);
	  goodmu = goodmu || (planes[4] && planes[2] && planes[3]);
	  for (int ik = 1; ik <= 4; ik++)
	    for (int jk = ik + 1; jk <= 4; jk++)
	      goodmu = goodmu || (planes[ik] && planes[jk]);
	  if (!goodmu)
	    continue;
	  if (!coinc)
	    continue;

	  TH2 *hzx = _rh->GetTH2("/gric/ZX");
	  TH2 *hzy = _rh->GetTH2("/gric/ZY");
	  if (hzx == NULL)
	    {
	      hzx = _rh->BookTH2("/gric/ZX", 100., 0., 400., 120., -40., 40.);
	      hzy = _rh->BookTH2("/gric/ZY", 100., 0., 400., 180., 0., 90.);
	    }
	  hzx->Reset();
	  hzy->Reset();

	  // Build Points list
	  _vPoints.clear();
	  _hplanes.reset();
	  for (int i = 1; i <= 8; i++)
	    if (planes[i] == 1)
	      buildPlaneHits(e, i, x.second);

	  // Track segments
	  // High part
	  bool goodhmu = (_hplanes[1] && _hplanes[2] && _hplanes[3] && _hplanes[4]);
	  goodhmu = goodhmu || (_hplanes[1] && _hplanes[2] && _hplanes[3]);
	  goodhmu = goodhmu || (_hplanes[1] && _hplanes[2] && _hplanes[4]);
	  goodhmu = goodhmu || (_hplanes[1] && _hplanes[3] && _hplanes[4]);
	  goodhmu = goodhmu || (_hplanes[4] && _hplanes[2] && _hplanes[3]);

	  for (int ik = 1; ik <= 4; ik++)
	    for (int jk = ik + 1; jk <= 4; jk++)
	      goodhmu = goodhmu || (_hplanes[ik] && _hplanes[jk]);
	  if (!goodhmu)
	    continue;
	  // Low part
	  bool goodhkick = (_hplanes[5] && _hplanes[6] && _hplanes[7] && _hplanes[8]);
	  goodhkick = goodhkick || (_hplanes[5] && _hplanes[6] && _hplanes[7]);
	  goodhkick = goodhkick || (_hplanes[5] && _hplanes[6] && _hplanes[8]);
	  goodhkick = goodhkick || (_hplanes[8] && _hplanes[6] && _hplanes[7]);
	  goodhkick = goodhkick || (_hplanes[8] && _hplanes[5] && _hplanes[7]);

	  for (int ik = 5; ik <= 8; ik++)
	    for (int jk = ik + 1; jk <= 8; jk++)
	      goodhkick = goodhkick || (_hplanes[ik] && _hplanes[jk]);

	  // Stopped
	  bool stopmu = !(_hplanes[5] || _hplanes[6] || _hplanes[7] || _hplanes[8]);

	  // Allways ask for a top segment
	  if (!goodhmu)
	    continue;
	  if (!goodhkick)
	    continue;
	  // Build a top segments
	  float px2cut = 0.01, thcut = 0.6, dcut = 3.;
	  ShowerParams isha;
	  double ax = 0, bx = 0, ay = 0, by = 0;
	  float zmin = z[4] - 1, zmax = z[1] + 1;
	  int ier = TPrincipalComponents((double *)&isha, zmin, zmax);
	  double *sx = isha.xm;
	  double *sv = isha.l2;
	  double z0 = sx[2];
	  double x0 = sx[0];
	  double y0 = sx[1];
	  double x1 = sx[0] + sv[0];
	  double y1 = sx[1] + sv[1];
	  double z1 = sx[2] + sv[2];
	  //double ax,ay,bx,by;
	  if (ier == 0)
	    {
	      ax = (x1 - x0) / (z1 - z0);
	      bx = x1 - ax * z1;
	      ay = (y1 - y0) / (z1 - z0);
	      by = y1 - ay * z1;
	    }

	  top_tk.clear();
	  _t_h = 0;
	  _b_h = 0;
	  _a_h = 0;
	  top_tk.setDir(ax, ay, 1.);
	  top_tk.setOrig(bx, by, 0);
	  //printf(" ax %f bx %f ay %f by %f \n",ax,ay,bx,by);
	  for (auto x = _vPoints.begin(); x != _vPoints.end(); x++)
	    {
	      recoPoint &p = (*x);
	      if (p.Z() < zmin)
		continue;
	      if (p.Z() > zmax)
		continue;

	      if (top_tk.distance(&p) < dcut || ier != 0)
		{
		  top_tk.addPoint(&p);
		  _t_h |= (1 << p.plan() - 1);
		}
	    }
	  // Ask at least 3 points
	  if (top_tk.size() < 2)
	    continue;
	  top_tk.regression();
	  top_tk.calculateChi2();

	  if (top_tk.pchi2() < px2cut)
	    continue;
	  if (abs(top_tk.dir().X()) > thcut)
	    continue;
	  if (abs(top_tk.dir().Y()) > thcut)
	    continue;

	  _t_c2 = top_tk.pchi2();
	  _t_x[0] = top_tk.orig().X();
	  _t_x[1] = top_tk.orig().Y();
	  _t_x[2] = top_tk.orig().Z();
	  _t_v[0] = top_tk.dir().X();
	  _t_v[1] = top_tk.dir().Y();
	  _t_v[2] = top_tk.dir().Z();

	  if (goodhkick || true)
	    {
	      zmin = z[8] - 1, zmax = z[5] + 1;
	      memset(&isha, 0, sizeof(ShowerParams));
	      ier = TPrincipalComponents((double *)&isha, zmin, zmax);
	      sx = isha.xm;
	      sv = isha.l2;
	      z0 = sx[2];
	      x0 = sx[0];
	      y0 = sx[1];
	      x1 = sx[0] + sv[0];
	      y1 = sx[1] + sv[1];
	      z1 = sx[2] + sv[2];
	      //double ax,ay,bx,by;
	      if (ier == 0)
		{
		  ax = (x1 - x0) / (z1 - z0);
		  bx = x1 - ax * z1;
		  ay = (y1 - y0) / (z1 - z0);
		  by = y1 - ay * z1;
		}
	      bot_tk.clear();
	      bot_tk.setDir(ax, ay, 1.);
	      bot_tk.setOrig(bx, by, 0);
	      //printf("BOT ax %f bx %f ay %f by %f \n",ax,ay,bx,by);
	      //fflush(stdout);
	      for (auto x = _vPoints.begin(); x != _vPoints.end(); x++)
		{
		  recoPoint &p = (*x);
		  if (p.Z() < zmin)
		    continue;
		  if (p.Z() > zmax)
		    continue;

		  if (bot_tk.distance(&p) < dcut || ier != 0)
		    {
		      bot_tk.addPoint(&p);
		      _b_h |= (1 << p.plan() - 1);
		    }
		}
	      if (bot_tk.size() > 2)
		{
		  bot_tk.regression();
		  bot_tk.calculateChi2();
		}
	      //stopmu=true;
	      // Ask at least 3 points
	      ROOT::Math::XYZPoint pe = top_tk.extrapolate((z[5] + z[4]) / 2.);
	      if (bot_tk.size() >= 2 && bot_tk.pchi2() > px2cut)
		{
		  //stopmu=false;
		  bot_tk.regression();
		  bot_tk.calculateChi2();

		  _b_c2 = bot_tk.pchi2();
		  _b_x[0] = bot_tk.orig().X();
		  _b_x[1] = bot_tk.orig().Y();
		  _b_x[2] = bot_tk.orig().Z();
		  _b_v[0] = bot_tk.dir().X();
		  _b_v[1] = bot_tk.dir().Y();
		  _b_v[2] = bot_tk.dir().Z();

		  std::stringstream splane;
		  splane << "/gric/SEG/";
		  TH2 *hextb = _rh->GetTH2(splane.str() + "XYextBig");
		  TH3 *h3b = _rh->GetTH3(splane.str() + "XYZBig");
		  TH2 *hzxb = _rh->GetTH2(splane.str() + "ZXextBig");
		  TH2 *hzyb = _rh->GetTH2(splane.str() + "ZYextBig");
		  TH2 *hextc = _rh->GetTH2(splane.str() + "XYextCut");

		  TH2 *hextl = _rh->GetTH2(splane.str() + "XYextLow");

		  TH2 *hcpathet = _rh->GetTH2(splane.str() + "CPAvsTheta");
		  TH2 *hcpatheto = _rh->GetTH2(splane.str() + "CPAvsThetaOut");
		  TH1 *hthet = _rh->GetTH1(splane.str() + "Theta");
		  TH1 *hzbig = _rh->GetTH1(splane.str() + "zBig");
		  TH1 *hzlow = _rh->GetTH1(splane.str() + "zLow");
		  TH1 *hthi = _rh->GetTH1(splane.str() + "ThetaIn");
		  TH1 *htho = _rh->GetTH1(splane.str() + "ThetaOut");
		  TH1 *hcpa = _rh->GetTH1(splane.str() + "CPA");
		  TH1 *hcpai = _rh->GetTH1(splane.str() + "CPAI");
		  TH1 *hcpao = _rh->GetTH1(splane.str() + "CPAO");
		  TH1 *hrd3 = _rh->GetTH1(splane.str() + "Distance");
		  TH1 *hprob = _rh->GetTH1(splane.str() + "Prob");
		  TH1 *hprobc = _rh->GetTH1(splane.str() + "ProbCut");
		  if (hextb == NULL)
		    {
		      hthet = _rh->BookTH1(splane.str() + "Theta", 360, 0., 90.);
		      hcpathet = _rh->BookTH2(splane.str() + "CPAvsTheta", 80, 0., 20., 50, 0., 12.5);
		      hcpatheto = _rh->BookTH2(splane.str() + "CPAvsThetaOut", 80, 0., 20., 50, 0., 12.5);
		      hthi = _rh->BookTH1(splane.str() + "ThetaIn", 360, 0., 90.);
		      hzbig = _rh->BookTH1(splane.str() + "zBig", 300, 0., 300.);
		      hzlow = _rh->BookTH1(splane.str() + "zLow", 300, 0., 300.);
		      htho = _rh->BookTH1(splane.str() + "ThetaOut", 360, 0., 90.);
		      hcpa = _rh->BookTH1(splane.str() + "CPA", 200, 0., 50.);
		      hrd3 = _rh->BookTH1(splane.str() + "Distance", 200, 0., 10.);
		      hprob = _rh->BookTH1(splane.str() + "Prob", 200, 0., 1.);
		      hprobc = _rh->BookTH1(splane.str() + "ProbCut", 200, 0., 1.);
		      hcpai = _rh->BookTH1(splane.str() + "CPAI", 200, 0., 50.);
		      hcpao = _rh->BookTH1(splane.str() + "CPAO", 200, 0., 50.);
		      hextb = _rh->BookTH2(splane.str() + "XYextBig", 24, -30., 30., 24, 0., 60.);
		      h3b = _rh->BookTH3(splane.str() + "XYZBig", 16, -30., 30., 16, 0., 60., 60, 0., 240.);
		      hzxb = _rh->BookTH2(splane.str() + "ZXextBig", 50, z[5], z[4], 32, -30., 30.);
		      hzyb = _rh->BookTH2(splane.str() + "ZYextBig", 50, z[5], z[4], 32, 0., 60.);
		      hextc = _rh->BookTH2(splane.str() + "XYextCut", 32, -30., 30., 32, 0., 60.);
		      hextl = _rh->BookTH2(splane.str() + "XYextLow", 32, -30., 30., 32, 0., 60.);
		    }

		  double sc = top_tk.dir().Dot(bot_tk.dir()) / sqrt(top_tk.dir().Mag2()) / sqrt(bot_tk.dir().Mag2());
		  double th = acos(sc) * 180 / M_PI;
		  std::cout << "ICI COSTH " << sc << " DEG " << th << std::endl;

		  _cos_th = sc;
		  _th = th;
		  if (th > 15)
		    continue;

		  double adist;
		  ROOT::Math::XYZPoint p1, p2;
		  top_tk.cap(bot_tk, adist, p1, p2);

		  std::cout << p1.X() << ":" << p1.Y() << ":" << p1.Z() << std::endl;
		  std::cout << p2.X() << ":" << p2.Y() << ":" << p2.Z() << std::endl;
		  float zcross = (p1.Z() + p2.Z()) / 2.;
		  //printf("%x %x %x \n",hcpa,hcpai,hcpao);
		  std::cout << "LA DISTANCE " << adist << " Z CROSS " << zcross << std::endl;
		  if (zcross < z[5] + 1 || zcross > z[4] - 1)
		    continue;
		  fflush(stdout);

		  _xcross = (p1.X() + p2.X()) / 2.;
		  _ycross = (p1.Y() + p2.Y()) / 2.;
		  _zcross = zcross;
		  _dist = adist;

		  ROOT::Math::XYZPoint po3 = bot_tk.extrapolate(z[4]);
		  ROOT::Math::XYZPoint pi3 = top_tk.extrapolate(z[4]);
		  ROOT::Math::XYZVector d3 = po3 - pi3;

		  //std::cout<<p1.X()<<" "<<p1.Y()<<" "<<p1.Z()<<std::endl;
		  double rd3 = sqrt(d3.Mag2());
		  double prob = 1.;
		  hrd3->Fill(rd3);
		  if (rd3 < 8)
		    {
		      prob = 1. - erf(rd3 / sqrt(2.) / 2.);
		      if (prob < 1E-20)
			prob = 1E-20;
		      printf("Distance %f %f %f %f \n", rd3, prob, th, thetacut);
		      if (th > thetacut)
			hprobc->Fill(prob);
		      else
			hprob->Fill(prob);
		    }

		  float fdist = adist;
		  hcpa->Fill(fdist);

		  if (prob < pcut)
		    hthi->Fill(th);
		  else
		    htho->Fill(th);
		  if (zcross < z[4] - 1 && zcross > z[5] + 1)
		    {

		      hcpathet->Fill(th, fdist);
		    }
		  else
		    {
		      hcpatheto->Fill(th, fdist);
		    }
		  if (th > thetacut && prob < pcut && th < 15 && fdist < cpacut)
		    {
		      hcpai->Fill(fdist);
		      hzbig->Fill(zcross);
		      hextc->Fill((p1.X() + p2.X()) / 2., (p1.Y() + p2.Y()) / 2.);
		    }
		  else
		    {
		      hcpao->Fill(fdist);
		      hzlow->Fill(zcross);
		    }
		  hthet->Fill(th);

		  if (th > thetacut && fdist < cpacut && prob < pcut && zcross < z[4] - 1 && zcross > z[5] + 1)
		    {
		      hextb->Fill((p1.X() + p2.X()) / 2., (p1.Y() + p2.Y()) / 2.);
		      h3b->Fill((p1.X() + p2.X()) / 2., (p1.Y() + p2.Y()) / 2., _zcross);
		      if (th > 1.5 * thetacut)
			{
			  hzxb->Fill(zcross, (p1.X() + p2.X()) / 2.);
			  hzyb->Fill(zcross, (p1.Y() + p2.Y()) / 2.);
			}
		    }
		  else
		    hextl->Fill(pe.X(), pe.Y());
		}
	      else if (stopmu)
		{

		  //stop mu
		  std::stringstream splane;
		  splane << "/gric/SEG/";
		  TH2 *hstop = _rh->GetTH2(splane.str() + "XYStop");
		  TH2 *hcomb = _rh->GetTH2(splane.str() + "XYComb");
		  TH2 *hzxcomb = _rh->GetTH2(splane.str() + "ZXComb");
		  TH2 *hzycomb = _rh->GetTH2(splane.str() + "ZYComb");
		  if (hstop == NULL)
		    {
		      hstop = _rh->BookTH2(splane.str() + "XYStop", 32, -30., 30., 32, 0., 64.);
		      hcomb = _rh->BookTH2(splane.str() + "XYComb", 32, -30., 30., 32, 0., 64.);
		      hzxcomb = _rh->BookTH2(splane.str() + "ZXComb", 50, z[5], z[4], 32, -30., 30.);
		      hzycomb = _rh->BookTH2(splane.str() + "ZYComb", 50, z[5], z[4], 32, 0., 64.);
		    }

		  hstop->Fill(pe.X(), pe.Y());
		  pe = top_tk.extrapolate(z[8]);
		  bool tkok = top_tk.size() >= 3;
		  if (pe.X() < -27.648 + 3)
		    tkok = false;
		  if (pe.X() > 27.648 - 3)
		    tkok = false;
		  if (pe.X() < 0)
		    {
		      if (pe.Y() < -1. * T30 * pe.X() + 3)
			tkok = false;
		      if (pe.Y() > T30 * pe.X() + 55.296 - 3)
			tkok = false;
		    }
		  else
		    {
		      if (pe.Y() < T30 * pe.X() + 3)
			tkok = false;
		      if (pe.Y() > -1. * T30 * pe.X() + 55.296 - 3)
			tkok = false;
		    }
		  if (tkok)
		    {
		      vstop[nstop] = top_tk;
		      nstop++;
		    }
		  if (_event % 10 == 0)
		    {
		      TH1 *hftm = _rh->GetTH1("/gric/MaxTime");
		      if (hftm != NULL)
			printf("============================> Event %d Stopped %d  %f %f \n", _event, nstop, hftm->GetEntries(), hftm->GetEntries() * hftm->GetMean());

		      fflush(stdout);
		    }
		  if (nstop % 20 == 0 && nstop > 0)
		    {
		      if (nstop > 1000)
			continue;
		      if (nstop == lastprocessed)
			continue;
		      lastprocessed = nstop;
		      hcomb->Reset();
		      for (int i = 0; i < nstop; i++)
			for (int j = i + 1; j < nstop; j++)
			  {
			    double adist;
			    ROOT::Math::XYZPoint p1, p2;
			    vstop[i].cap(vstop[j], adist, p1, p2);
			    //std::cout<<adist<<std::endl;
			    ///std::cout<<p1.X()<<":"<<p1.Y()<<":"<<p1.Z()<<std::endl;
			    //std::cout<<p2.X()<<":"<<p2.Y()<<":"<<p2.Z()<<std::endl;
			    double zcross = (p1.Z() + p2.Z()) / 2.;
			    //printf("%x %x %x \n",hcpa,hcpai,hcpao);
			    //fflush(stdout);

			    double xcross = (p1.X() + p2.X()) / 2.;
			    double ycross = (p1.Y() + p2.Y()) / 2.;
			    if (zcross > z[5] && zcross < z[4] && adist > 3. && adist < 13.)
			      {
				hcomb->Fill(xcross, ycross);
				hzxcomb->Fill(zcross, xcross);
				hzycomb->Fill(zcross, ycross);
			      }
			  }

		      std::cout << "Combinatories " << hcomb->GetEntries() << std::endl;
		    }
		}
	    }

	  // Fill Histo

	  for (int ip = 1; ip <= 8; ip++)
	    {
	      if (z[ip] == 0)
		continue;
	      recoTrack tp;
	      tp.clear();
	      for (auto x = top_tk.points().begin(); x != top_tk.points().end(); x++)
		{
		  recoPoint *p = (*x);
		  if (p->plan() == ip)
		    continue;
		  tp.addPoint(p);
		}
	      for (auto x = bot_tk.points().begin(); x != bot_tk.points().end(); x++)
		{
		  recoPoint *p = (*x);
		  if (p->plan() == ip)
		    continue;
		  tp.addPoint(p);
		}
	      if (tp.size() < 3)
		continue;
	      tp.regression();
	      tp.calculateChi2();
	      if (tp.pchi2() < px2cut)
		continue;

	      if (ip == 1 && !(tp.planUsed(2) && tp.planUsed(3)))
		continue;
	      if (ip == 1 && !(tp.planUsed(2) && tp.planUsed(8)))
		continue;

	      if (ip == 8 && !(tp.planUsed(6) && tp.planUsed(7)))
		continue;
	      if (ip == 8 && !(tp.planUsed(6) && tp.planUsed(1)))
		continue;
	      if (ip > 1 && ip < 8 && !(tp.planUsed(ip - 1) && tp.planUsed(ip + 1)))
		continue;
	      if (ip > 1 && ip < 8 && !(tp.planUsed(1) && tp.planUsed(8)))
		continue;
	      ROOT::Math::XYZPoint pe = tp.extrapolate(z[ip]);

	      if (pe.X() < -27.648 + 3)
		continue;
	      if (pe.X() > 27.648 - 3)
		continue;
	      if (pe.X() < 0)
		{
		  if (pe.Y() < -1. * T30 * pe.X() + 3)
		    continue;
		  if (pe.Y() > T30 * pe.X() + 55.296 - 3)
		    continue;
		}
	      else
		{
		  if (pe.Y() < T30 * pe.X() + 3)
		    continue;
		  if (pe.Y() > -1. * T30 * pe.X() + 55.296 - 3)
		    continue;
		}
	      std::stringstream splane;
	      splane << "/gric/PLANE" << ip << "/";
	      TH2 *hext = _rh->GetTH2(splane.str() + "XYext");
	      TH2 *hfop = _rh->GetTH2(splane.str() + "XYfound");

	      TH1 *hdistx = _rh->GetTH1(splane.str() + "DX");
	      TH1 *hdisty = _rh->GetTH1(splane.str() + "DY");
	      TH1 *hax = _rh->GetTH1(splane.str() + "AX");
	      TH1 *hay = _rh->GetTH1(splane.str() + "AY");
	      TH1 *hpc2 = _rh->GetTH1(splane.str() + "PChi2");
	      if (hdistx == NULL)
		{
		  hpc2 = _rh->BookTH1(splane.str() + "PChi2", 1000, 0., 1.);
		  hax = _rh->BookTH1(splane.str() + "AX", 1000, -1., 1.);
		  hay = _rh->BookTH1(splane.str() + "AY", 1000, -1., 1.);
		  hdistx = _rh->BookTH1(splane.str() + "DX", 100, -10., 10.);
		  hdisty = _rh->BookTH1(splane.str() + "DY", 100, -10., 10.);
		  hext = _rh->BookTH2(splane.str() + "XYext", 32, -30., 30., 32, 0., 60.);
		  hfop = _rh->BookTH2(splane.str() + "XYfound", 32, -30., 30., 32, 0., 60.);
		}
	      hpc2->Fill(tp.pchi2());
	      hax->Fill(tp.dir().X());
	      hay->Fill(tp.dir().Y());

	      //if ( abs(tp.dir().X())>thcut || abs(tp.dir().Y())>thcut) continue;
	      hext->Fill(pe.X(), pe.Y());

	      for (auto x = _vPoints.begin(); x != _vPoints.end(); x++)
		{
		  if (x->plan() != ip)
		    continue;
		  ROOT::Math::XYZPoint &p = (*x);
		  ROOT::Math::XYZVector d = p - pe;
		  float dx = d.X();
		  float dy = d.Y();
		  //std::cout<<p1.X()<<" "<<p1.Y()<<" "<<p1.Z()<<std::endl;
		  double dist = sqrt(d.Mag2());
		  //std::cout<<dist<<" "<<hdistx->GetEntries()<<std::endl;
		  hdistx->Fill(dx);
		  //std::cout<<dx<<" "<<dy<<" "<<hdistx->GetEntries()<<std::endl;
		  hdisty->Fill(dy);
		  if (dist < 5.)
		    {
		      hfop->Fill(pe.X(), pe.Y());
		      break;
		    }
		}
	    }

	  bool display = _jparams["general"]["display"].asUInt() == 1;
	  if (display)
	    {
	      if (TCHits == NULL)
		{
		  TCHits = new TCanvas("TCHits", "tChits1", 900, 900);
		  TCHits->Modified();
		  TCHits->Draw();
		  TCHits->Divide(1, 2);
		}
	      TCHits->cd(1);
	      hzx->SetMarkerStyle(25);
	      hzx->SetMarkerColor(kRed);
	      hzx->Draw("P");
	      top_tk.linex()->Draw("SAME");

	      std::cout << top_tk << std::endl;
	      std::cout << top_tk.pchi2() << std::endl;
	      std::cout << "KICK ?:" << goodhkick << " HPLANE" << _hplanes << std::endl;

	      if (bot_tk.size() > 1)
		{
		  bot_tk.linex()->Draw("SAME");

		  std::cout << bot_tk << std::endl;
		  std::cout << bot_tk.pchi2() << std::endl;

		  double sc = top_tk.dir().Dot(bot_tk.dir()) / sqrt(top_tk.dir().Mag2()) / sqrt(bot_tk.dir().Mag2());
		  std::cout << "COSTH " << sc << " DEG " << acos(sc) * 180 / M_PI << std::endl;
		}
	      TCHits->Modified();
	      TCHits->Draw();
	      TCHits->Update();
	      TCHits->cd(2);
	      hzy->SetMarkerStyle(22);
	      hzy->SetMarkerColor(kGreen);
	      hzy->Draw("P");
	      top_tk.liney()->Draw("SAME");

	      if (bot_tk.size() > 1)
		{
		  bot_tk.liney()->Draw("SAME");
		}

	      TCHits->Modified();
	      TCHits->Draw();
	      TCHits->Update();
	      sleep(1);
	      TCHits->Update();
	      //std::stringstream ss("");
	      //ss<<"/tmp/Display_"<<evt_->getRunNumber()<<"_"<<evt_->getEventNumber()<<"_"<<vgood[is]<<".png";
	      //TCHits->SaveAs(ss.str().c_str());

	      //std::cout<<"HITS :"<<_hplanes.count()<<" "<<_hplanes<<std::endl;

	      // Hough Transform
#undef HTFIT
#ifdef HTFIT
	      HoughLocal htl(0, M_PI, -360., 360, 16, 64);
	      htl.clear();
	      htl.fill(&_vPoints);
	      std::vector<std::pair<uint32_t, uint32_t>> hbins;
	      hbins.clear();
	      htl.findMaxima(hbins, nplanesmin - 1);
	      /*
		for (std::vector < std::pair<uint32_t,uint32_t> >::iterator ihb=hbins.begin();ihb<hbins.end();ihb++)
		{
		printf("Bin %d %d => %d \n",(*ihb).first,(*ihb).second,htl.getValue((*ihb).first,(*ihb).second));
		}
	      */
	      htl.draw(_rh);

	      HoughLocal htly(0, M_PI, -350., 350, 16, 64);
	      htly.clear();
	      htly.fill(&_vPoints, 1);
	      std::vector<std::pair<uint32_t, uint32_t>> hbinsy;
	      hbinsy.clear();
	      htly.findMaxima(hbinsy, nplanesmin - 1);
	      /*
		for (std::vector < std::pair<uint32_t,uint32_t> >::iterator ihb=hbinsy.begin();ihb<hbinsy.end();ihb++)
		{
		printf("Bin %d %d => %d \n",(*ihb).first,(*ihb).second,htly.getValue((*ihb).first,(*ihb).second));
		}
	      */
	      htly.draw(_rh);

	      std::cout << _vPoints.size() << " -> " << hbins.size() << " " << hbinsy.size() << std::endl;
	      for (auto x = _vPoints.begin(); x != _vPoints.end(); x++)
		{
		  ROOT::Math::XYZPoint &p = (*x);

		  std::cout << p.X() << ":" << p.Y() << ":" << p.Z() << std::endl;
		}

	      if (hbins.size() == 0)
		continue;
	      for (std::vector<std::pair<uint32_t, uint32_t>>::iterator ihb = hbins.begin(); ihb != hbins.end(); ihb++)
		{
		  uint32_t ith = (*ihb).first;
		  uint32_t ir = (*ihb).second;
		  float theta = (htl.getTheta(ith) + htl.getTheta(ith + 1)) / 2;
		  float r = (htl.getR(ir) + htl.getR(ir + 1)) / 2;
		  float a = -1. / tan(theta);
		  float b = r / sin(theta);
		  ax = a;
		  bx = b;
		}
	      if (hbinsy.size() == 0)
		continue;
	      for (std::vector<std::pair<uint32_t, uint32_t>>::iterator ihb = hbinsy.begin(); ihb != hbinsy.end(); ihb++)
		{
		  uint32_t ith = (*ihb).first;
		  uint32_t ir = (*ihb).second;
		  float theta = (htly.getTheta(ith) + htly.getTheta(ith + 1)) / 2;
		  float r = (htly.getR(ir) + htly.getR(ir + 1)) / 2;
		  float a = -1. / tan(theta);
		  float b = r / sin(theta);
		  ay = a;
		  by = b;
		}

#endif
	      getchar();
	    }
	  if (true)
	    continue;

	  for (auto y : x.second)
	    {
	      uint32_t ig = y / MAXFRAME / FSIZE;
	      uint16_t plane = (ig >> 4) & 0xF;
	      bool muv = u[plane] && v[plane];
	      bool muw = u[plane] && w[plane];
	      bool mvw = v[plane] && w[plane];
	      if (!(muv || muw || mvw))
		continue;

	      uint8_t cpos = (ig & 0xF);
	      bool udir = ((ig & 0xF) == 2 || (ig & 0xF) == 4);
	      bool vdir = ((ig & 0xF) == 3 || (ig & 0xF) == 5);
	      bool wdir = ((ig & 0xF) == 6 || (ig & 0xF) == 7);
	      int sshift = 0;
	      if (cpos == 4 || cpos == 5 || cpos == 7)
		sshift = 64;
	      //if (!measure) continue;
	      std::stringstream sraw1;
	      std::stringstream splane;

	      sraw1 << "/gric/ASIC" << std::hex << ig << std::dec << "/";
	      splane << "/gric/PLANE" << plane << "/";

	      TH1 *hp1 = _rh->GetTH1(sraw1.str() + "Pad1Sel");
	      TH1 *hsu = _rh->GetTH1(splane.str() + "Ustrip");
	      TH1 *hsv = _rh->GetTH1(splane.str() + "Vstrip");
	      TH1 *hsw = _rh->GetTH1(splane.str() + "Wstrip");
	      TH1 *hfcs = _rh->GetTH1(sraw1.str() + "FrameCountSel");
	      TH1 *hftse = _rh->GetTH1(sraw1.str() + "FrameTimeSel");

	      if (hp1 == NULL)
		{
		  hp1 = _rh->BookTH1(sraw1.str() + "Pad1Sel", 64, 0., 64.);
		  hfcs = _rh->BookTH1(sraw1.str() + "FrameCountSel", 65, 0., 65.);
		  hftse = _rh->BookTH1(sraw1.str() + "FrameTimeSel", 500, 0., 500.);
		  hsu = _rh->BookTH1(splane.str() + "Ustrip", 128, 0., 128.);
		  hsv = _rh->BookTH1(splane.str() + "Vstrip", 128, 0., 128.);
		  hsw = _rh->BookTH1(splane.str() + "Wstrip", 128, 0., 128.);
		}
	      hftse->Fill(_maxTime - e->bcid(y) * 1.);
	      uint32_t idx = y;
	      std::bitset<64> bs;
	      bs.reset();
	      for (int k = 0; k < 64; k++)
		if (e->pad0(idx, k) || e->pad1(idx, k))
		  {
		    bs.set(k);
		    hp1->Fill(k * 1.);
		    if (udir)
		      hsu->Fill(k + sshift * 1.);
		    if (vdir)
		      hsv->Fill(k + sshift * 1.);
		    if (wdir)
		      hsw->Fill(k + sshift * 1.);
		  }
	      hfcs->Fill(bs.count() * 1.);
	      //fprintf(stderr,"%x frame bcid %.9d %8.5f %s %d\n",ig,e->bcid(idx),e->bcid(idx)*2E-7,bs.to_string().c_str(),bs.count());
	      // getchar();
	    }
	}
      //if (trigger&&!ramf) getchar();
    }

  return;
}
void binaryreader::buildPlaneHits(rbEvent *e, uint32_t plane, std::vector<uint32_t> &hits)

{

  memset(&_bs[6 * plane], 0, 6 * sizeof(uint64_t));
  TH2 *hzx = _rh->GetTH2("/gric/ZX");
  TH2 *hzy = _rh->GetTH2("/gric/ZY");
  if (hzx == NULL)
    {
      hzx = _rh->BookTH2("/gric/ZX", 100., 0., 100., 120., -60., 60.);
      hzy = _rh->BookTH2("/gric/ZY", 100., 0., 100., 120., -60., 60.);
    }
  std::stringstream splane;

  splane << "/gric/PLANE" << plane << "/";
  std::bitset<128> ub;
  ub.reset();
  std::bitset<128> vb;
  vb.reset();
  std::bitset<128> wb;
  wb.reset();
  std::bitset<128> ub2;
  ub2.reset();
  std::bitset<128> vb2;
  vb2.reset();
  std::bitset<128> wb2;
  wb2.reset();
  std::bitset<128> ub3;
  ub3.reset();
  std::bitset<128> vb3;
  vb3.reset();
  std::bitset<128> wb3;
  wb3.reset();
  _hplanes.reset(plane);
  //  std::cout<<"Hit Plane before "<<plane<<" BS "<<_hplanes<<std::endl;
  uint32_t igplane = 0;
  for (auto y : hits)
    {
      // Select the plane
      uint32_t ig = y / MAXFRAME / FSIZE;
      uint16_t igpl = (ig >> 4) & 0xF;
      uint16_t pl = _plinfo[igpl]["num"].asUInt();

      if (pl != plane)
	continue;
      igplane = igpl;
      uint16_t dir = (ig & 0xF);
      bool udir = (dir == 2 || dir == 4);
      bool vdir = (dir == 3 || dir == 5);
      bool wdir = (dir == 6 || dir == 7);
      int sshift = 0;
      if (dir == 4 || dir == 5 || dir == 7)
	sshift = 64;
      for (int k = 0; k < 64; k++)
	{
	  if (e->pad0(y, k) || e->pad1(y, k))
	    {
	      if (udir)
		ub.set(k + sshift);
	      if (vdir)
		vb.set(k + sshift);
	      if (wdir)
		wb.set(k + sshift);
	    }
	  if (e->pad0(y, k) && !e->pad1(y, k))
	    {
	      if (udir)
		ub2.set(k + sshift);
	      if (vdir)
		vb2.set(k + sshift);
	      if (wdir)
		wb2.set(k + sshift);
	    }
	  if (e->pad0(y, k) && e->pad1(y, k))
	    {
	      if (udir)
		ub3.set(k + sshift);
	      if (vdir)
		vb3.set(k + sshift);
	      if (wdir)
		wb3.set(k + sshift);
	    }
	}
    }
  float x_u = -1, x_v = -1, x_w = -1;
  float s_u = -1, s_v = -1, s_w = -1;

  std::vector<float> vx_u;
  vx_u.clear();
  std::vector<float> vx_v;
  vx_v.clear();
  std::vector<float> vx_w;
  vx_w.clear();
  std::vector<float> sx_u;
  std::vector<float> sx_v;
  std::vector<float> sx_w;

#define MAXHITONLY
#ifdef MAXHITONLY
  if (ub2.count() != 0)
    ub = ub2;
  if (vb2.count() != 0)
    vb = vb2;
  if (wb2.count() != 0)
    wb = wb2;

  if (ub3.count() != 0)
    ub = ub3;
  if (vb3.count() != 0)
    vb = vb3;
  if (wb3.count() != 0)
    wb = wb3;
  /*  */
#endif
  if (ub.count() > 0)
    {
      int f = 1000, l = -1;
      std::vector<int> first, last;
      for (int i = 0; i < 128; i++)
	{
	  if (ub[i] && f > 128)
	    f = i;
	  if (i < 64 && ub[i])
	    _bs[6 * (plane - 1)] |= (1 << i);
	  if (i >= 64 && ub[i])
	    _bs[6 * (plane - 1) + 1] |= (1 << (i - 64));
	  if (ub[i] == 0 && f < 128 && l < 0)
	    {
	      l = i - 1;
	      first.push_back(f);
	      last.push_back(l);
	      f = 1000;
	      l = -1;
	    }
	}
      if (f < 128 && l < 0)
	{
	  first.push_back(f);
	  last.push_back(127);
	}

      int idx = -1;
      //printf("vx u size %d \n",first.size());
      for (int i = 0; i < first.size(); i++)
	{
	  vx_u.push_back((first[i] + last[i]) / 2.);
	  sx_u.push_back(last[i] - first[i] + 1);
	  if ((last[i] - first[i] + 1) > s_u)
	    {
	      s_u = (last[i] - first[i] + 1);
	      idx = i;
	    }
	}
      //std::cout<<"X IDX "<<idx<<std::endl;
      x_u = (first[idx] + last[idx]) / 2.;

      //std::cout<<" U :"<<f<<":"<<l<<":"<<x_u<<"=> "<<ub<<std::endl;
      //std::cout<<" U2 :"<<f<<":"<<l<<":"<<x_u<<"=> "<<ub2<<std::endl;
      TH1 *hnu = _rh->GetTH1(splane.str() + "UCount");

      if (hnu == NULL)
	{
	  hnu = _rh->BookTH1(splane.str() + "UCount", 128, 0., 128.);
	}
      hnu->Fill(s_u);
    }
  if (vb.count() > 0)
    {
      int f = 1000, l = -1;
      std::vector<int> first, last;
      for (int i = 0; i < 128; i++)
	{
	  if (vb[i] && f > 128)
	    f = i;
	  if (i < 64 && vb[i])
	    _bs[6 * (plane - 1) + 2] |= (1 << i);
	  if (i >= 64 && vb[i])
	    _bs[6 * (plane - 1) + 3] |= (1 << (i - 64));

	  if (vb[i] == 0 && f < 128 && l < 0)
	    {
	      l = i - 1;
	      first.push_back(f);
	      last.push_back(l);
	      f = 1000;
	      l = -1;
	    }
	}
      if (f < 128 && l < 0)
	{
	  first.push_back(f);
	  last.push_back(127);
	}

      int idx = -1;
      //printf("vx v size %d \n",first.size());
      for (int i = 0; i < first.size(); i++)
	{
	  vx_v.push_back((first[i] + last[i]) / 2.);
	  sx_v.push_back(last[i] - first[i] + 1);

	  if ((last[i] - first[i] + 1) > s_v)
	    {
	      s_v = (last[i] - first[i] + 1);
	      idx = i;
	    }
	}
      //std::cout<<"V IDX "<<first.size()<<idx<<std::endl;
      x_v = (first[idx] + last[idx]) / 2.;
      //std::cout<<" V :"<<f<<":"<<l<<":"<<x_v<<"=> "<<vb<<std::endl;
      //std::cout<<" V2 :"<<f<<":"<<l<<":"<<x_v<<"=> "<<vb2<<std::endl;
      TH1 *hnv = _rh->GetTH1(splane.str() + "VCount");

      if (hnv == NULL)
	{
	  hnv = _rh->BookTH1(splane.str() + "VCount", 128, 0., 128.);
	}
      hnv->Fill(s_v);
    }
  if (wb.count() > 0)
    {
      int f = 1000, l = -1;
      std::vector<int> first, last;
      for (int i = 0; i < 128; i++)
	{
	  if (wb[i] && f > 128)
	    f = i;
	  if (i < 64 && wb[i])
	    _bs[6 * (plane - 1) + 4] |= (1 << i);
	  if (i >= 64 && wb[i])
	    _bs[6 * (plane - 1) + 5] |= (1 << (i - 64));

	  if (wb[i] == 0 && f < 128 && l < 0)
	    {
	      l = i - 1;
	      first.push_back(f);
	      last.push_back(l);
	      f = 1000;
	      l = -1;
	    }
	}
      if (f < 128 && l < 0)
	{
	  first.push_back(f);
	  last.push_back(127);
	}

      // printf("vx w size %d \n",first.size());
      int idx = -1;
      for (int i = 0; i < first.size(); i++)
	{
	  vx_w.push_back((first[i] + last[i]) / 2.);
	  sx_w.push_back(last[i] - first[i] + 1);

	  if ((last[i] - first[i] + 1) > s_w)
	    {
	      s_w = (last[i] - first[i] + 1);
	      idx = i;
	    }
	}
      //std::cout<<"W IDX "<<idx<<std::endl;
      x_w = (first[idx] + last[idx]) / 2.;
      //std::cout<<" W :"<<f<<":"<<l<<":"<<x_w<<"=> "<<wb<<std::endl;
      TH1 *hnw = _rh->GetTH1(splane.str() + "WCount");

      if (hnw == NULL)
	{
	  hnw = _rh->BookTH1(splane.str() + "WCount", 128, 0., 128.);
	}
      hnw->Fill(s_w);
    }
  if (ub.count() > 0 || vb.count() > 0 || wb.count() > 0)
    {
      //std::cout<<ub.count()<<":"<<vb.count()<<":"<<wb.count()<<std::endl;
      //std::cout<<"COUCOU"<<_maxTime-x.first<<std::endl;
      //getchar();
    }
  float X_uv, X_uw, X_vw, X_uvw;
  float Y_uv, Y_uw, Y_vw, Y_uvw;
  bool buv = false, buw = false, bvw = false, buvw = false;
  TH2 *huvwo = _rh->GetTH2(splane.str() + "UVWOr");

  if (huvwo == NULL)
    {

      huvwo = _rh->BookTH2(splane.str() + "UVWOr", 100, -30., 30., 100, 0., 60.);
    }

  if (x_u >= 0 && x_v >= 0)
    {
      float a_v = -1. * T30, b_v = (128 - x_v) * TSTEP;
      float X = (x_u - 64) * XSTEP;
      float Y = a_v * X + b_v;
      X_uv = X;
      Y_uv = Y;
      TH2 *huv = _rh->GetTH2(splane.str() + "UV");
      TH2 *huvw = _rh->GetTH2(splane.str() + "UVW");

      if (huv == NULL)
	{
	  huv = _rh->BookTH2(splane.str() + "UV", 100, -30., 30., 100, 0., 60.);
	  huvw = _rh->BookTH2(splane.str() + "UVW", 100, -30., 30., 100, 0., 60.);
	}
      huv->Fill(X, Y);
      //	  hst->Fill(11.);
      if (x_w > 0)
	{
	  buvw = true;
	  X_uvw = X;
	  Y_uvw = Y;
	  huvw->Fill(X, Y);

	  TH1 *hn3u = _rh->GetTH1(splane.str() + "U3Count");
	  TH1 *hn3v = _rh->GetTH1(splane.str() + "V3Count");
	  TH1 *hn3w = _rh->GetTH1(splane.str() + "W3Count");

	  if (hn3w == NULL)
	    {
	      hn3w = _rh->BookTH1(splane.str() + "W3Count", 128, 0., 128.);
	      hn3v = _rh->BookTH1(splane.str() + "V3Count", 128, 0., 128.);
	      hn3u = _rh->BookTH1(splane.str() + "U3Count", 128, 0., 128.);
	    }
	  hn3u->Fill(s_u);
	  hn3v->Fill(s_v);
	  hn3w->Fill(s_w);
	}
      buv = true;
    }
  if (x_u >= 0 && x_w >= 0)
    {
      float a_w = 1. * T30, b_w = x_w * TSTEP;
      float X = (x_u - 64) * XSTEP;
      float Y = a_w * X + b_w;

      X_uw = X;
      Y_uw = Y;
      //std::cout<<x_u<<":"<<x_w<<"=>"<<X<<","<<Y<<std::endl;
      //getchar();

      TH2 *huw = _rh->GetTH2(splane.str() + "UW");

      if (huw == NULL)
	{
	  huw = _rh->BookTH2(splane.str() + "UW", 100, -30., 30., 100, 0., 60.);
	}
      huw->Fill(X, Y);

      buw = true;
    }
  if (x_w >= 0 && x_v >= 0)
    {
      float a_w = T30, a_v = -1. * T30, b_w = x_w * TSTEP, b_v = (128 - x_v) * TSTEP;
      float X = (b_v - b_w) / (a_w - a_v);
      float Y = a_w * X + b_w;

      X_vw = X;
      Y_vw = Y;
      TH2 *hvw = _rh->GetTH2(splane.str() + "VW");

      if (hvw == NULL)
	{
	  hvw = _rh->BookTH2(splane.str() + "VW", 100, -30., 30., 100, 0., 60.);
	}
      hvw->Fill(X, Y);

      bvw = true;
    }
  float X_g, Y_g;
  if (buvw)
    {
      X_g = (X_uv + X_uw + X_vw) / 3.;
      Y_g = (Y_uv + Y_uw + Y_vw) / 3.;

      for (auto xu : vx_u)
	for (auto xv : vx_v)
	  for (auto xw : vx_w)
	    {

	      float a_v = -1. * T30, b_v = (128 - xv) * TSTEP;
	      float Xuv = (xu - 64) * XSTEP;
	      float Yuv = a_v * Xuv + b_v;
	      float a_w = 1. * T30, b_w = xw * TSTEP;
	      float Xuw = (x_u - 64) * XSTEP;
	      float Yuw = a_w * Xuw + b_w;
	      a_w = T30;
	      a_v = -1. * T30, b_w = xw * TSTEP;
	      b_v = (128 - xv) * TSTEP;
	      float Xvw = (b_v - b_w) / (a_w - a_v);
	      float Yvw = a_w * Xvw + b_w;
	      float Xg = (Xuv + Xuw + Xvw) / 3., Yg = (Yuv + Yuw + Yvw) / 3.;

	      double a = sqrt((Xuv - Xuw) * (Xuv - Xuw) + (Yuv - Yuw) * (Yuv - Yuw));
	      double b = sqrt((Xuv - Xvw) * (Xuv - Xvw) + (Yuv - Yvw) * (Yuv - Yvw));
	      double c = sqrt((Xvw - Xuw) * (Xvw - Xuw) + (Yvw - Yuw) * (Yvw - Yuw));
	      double p = (a + b + c) / 2.;
	      double Sg = sqrt(p * (p - a) * (p - b) * (p - c));
	      if (Sg > 15)
		continue;
	      _hplanes.set(plane);
	      //printf("%d %f %f %f \n",plane,Xg,Yg,Sg);
	      recoPoint pg(Xg - _plinfo[igplane]["x"].asFloat(), Yg - _plinfo[igplane]["y"].asFloat(), z[plane], plane);
	      //std::cout<<" UVW "<<plane<<" "<<Xg-_plinfo[igplane]["x"].asFloat()<<" "<<Yg-_plinfo[igplane]["y"].asFloat()<<" "<<z[plane]<<std::endl;
	      _vPoints.push_back(pg);
	      hzx->Fill(z[plane], Xg);
	      hzy->Fill(z[plane], Yg);

	      //getchar();
	    }
    }
  else if (buv)
    {
      for (auto xu : vx_u)
	for (auto xv : vx_v)
	  {
	    _hplanes.set(plane);
	    float a_v = -1. * T30, b_v = (128 - xv) * TSTEP;
	    float Xuv = (xu - 64) * XSTEP;
	    float Yuv = a_v * Xuv + b_v;
	    //recoPoint pg(Xuv,Yuv,z[plane]);
	    recoPoint pg(Xuv - _plinfo[igplane]["x"].asFloat(), Yuv - _plinfo[igplane]["y"].asFloat(), z[plane], plane);
	    //std::cout<<" UV "<<plane<<" "<<Xuv-_plinfo[igplane]["x"].asFloat()<<" "<<Yuv-_plinfo[igplane]["y"].asFloat()<<" "<<z[plane]<<std::endl;
	    _vPoints.push_back(pg);
	    hzx->Fill(z[plane], Xuv);
	    hzy->Fill(z[plane], Yuv);

	    //getchar();
	  }

      X_g = X_uv;
      Y_g = Y_uv;
    }
  else if (buw)
    {

      X_g = X_uw;
      Y_g = Y_uw;
      for (auto xu : vx_u)
	for (auto xw : vx_w)
	  {
	    _hplanes.set(plane);
	    float a_w = 1. * T30, b_w = xw * TSTEP;
	    float Xuw = (x_u - 64) * XSTEP;
	    float Yuw = a_w * Xuw + b_w;
	    //printf("%d %f %f %f \n",plane,Xg,Yg,Sg);
	    //recoPoint pg(Xuw,Yuw,z[plane]);
	    recoPoint pg(Xuw - _plinfo[igplane]["x"].asFloat(), Yuw - _plinfo[igplane]["y"].asFloat(), z[plane], plane);
	    //std::cout<<" UW "<<plane<<" "<<Xuw-_plinfo[igplane]["x"].asFloat()<<" "<<Yuw-_plinfo[igplane]["y"].asFloat()<<" "<<z[plane]<<std::endl;
	    _vPoints.push_back(pg);
	    hzx->Fill(z[plane], Xuw);
	    hzy->Fill(z[plane], Yuw);

	    //getchar();
	  }
    }
  else if (bvw)
    {

      for (auto xv : vx_v)
	for (auto xw : vx_w)
	  {
	    _hplanes.set(plane);
	    float a_w = T30, a_v = -1. * T30, b_w = xw * TSTEP, b_v = (128 - xv) * TSTEP;
	    float Xvw = (b_v - b_w) / (a_w - a_v);
	    float Yvw = a_w * Xvw + b_w;
	    //printf("%d %f %f %f \n",plane,Xg,Yg,Sg);
	    //recoPoint pg(Xvw,Yvw,z[plane]);
	    recoPoint pg(Xvw - _plinfo[igplane]["x"].asFloat(), Yvw - _plinfo[igplane]["y"].asFloat(), z[plane], plane);
	    //std::cout<<" VW "<<plane<<" "<<Xvw-_plinfo[igplane]["x"].asFloat()<<" "<<Yvw-_plinfo[igplane]["y"].asFloat()<<" "<<z[plane]<<std::endl;
	    _vPoints.push_back(pg);
	    hzx->Fill(z[plane], Xvw);
	    hzy->Fill(z[plane], Yvw);

	    //getchar();
	  }
    }

  if (buv || buw || bvw)
    {
      huvwo->Fill(X_g, Y_g);
      //hzx->Fill(z[plane],X_g);
      //hzy->Fill(z[plane],Y_g);
      //ROOT::Math::XYZPoint p(X_g,Y_g,z[plane]);
      // Test      _vPoints.push_back(p);
    }
  //std::cout<<"Hit Plane after "<<plane<<" BS "<<_hplanes<<std::endl;
}
void binaryreader::buildPosition(rbEvent *e, uint32_t plane, uint32_t ddmax, bool all)
{
  TH1 *hst = _rh->GetTH1("/gric/Stat");
  std::stringstream splane;

  std::string suf;
  if (all)
    suf = "ALL";
  else
    suf = "TRIG";
  splane << "/gric/" << suf << "/" << _plinfo[plane]["name"].asString() << "_" << plane << "/";

  for (auto x : _timeMap)
    {
      // First cut on time
      if (ddmax > 0 && (_maxTime - x.first) > ddmax)
	continue;
      if (ddmax = 0 && (_maxTime - x.first) < 30)
	continue;
      if (x.second.size() < 2)
	continue; // At least 2 direction
      std::bitset<128> ub;
      ub.reset();
      std::bitset<128> vb;
      vb.reset();
      std::bitset<128> wb;
      wb.reset();
      for (auto y : x.second)
	{
	  // Select the plane
	  uint32_t ig = y / MAXFRAME / FSIZE;
	  uint16_t pl = (ig >> 4) & 0xF;

	  if (pl != plane)
	    continue;
	  uint16_t dir = (ig & 0xF);
	  bool udir = (dir == 2 || dir == 4);
	  bool vdir = (dir == 3 || dir == 5);
	  bool wdir = (dir == 6 || dir == 7);
	  int sshift = 0;
	  if (dir == 4 || dir == 5 || dir == 7)
	    sshift = 64;
	  for (int k = 0; k < 64; k++)
	    if (e->pad0(y, k) || e->pad1(y, k))
	      {
		if (udir)
		  ub.set(k + sshift);
		if (vdir)
		  vb.set(k + sshift);
		if (wdir)
		  wb.set(k + sshift);
	      }
	}
      float x_u = -1, x_v = -1, x_w = -1;
      float s_u = -1, s_v = -1, s_w = -1;
      if (ub.count() > 0)
	{
	  int f = 1000, l = -1;
	  std::vector<int> first, last;
	  for (int i = 0; i < 128; i++)
	    {
	      if (ub[i] && f > 128)
		f = i;
	      if (ub[i] == 0 && f < 128 && l < 0)
		{
		  l = i - 1;
		  first.push_back(f);
		  last.push_back(l);
		  f = 1000;
		  l = -1;
		}
	    }
	  if (f < 128 && l < 0)
	    {
	      first.push_back(f);
	      last.push_back(127);
	    }

	  int idx = -1;
	  for (int i = 0; i < first.size(); i++)
	    if ((last[i] - first[i] + 1) > s_u)
	      {
		s_u = (last[i] - first[i] + 1);
		idx = i;
	      }

	  //std::cout<<"X IDX "<<idx<<std::endl;
	  x_u = (first[idx] + last[idx]) / 2.;

	  //std::cout<<" U :"<<f<<":"<<l<<":"<<x_u<<"=> "<<ub<<std::endl;
	  TH1 *hnu = _rh->GetTH1(splane.str() + "UCount");

	  if (hnu == NULL)
	    {
	      hnu = _rh->BookTH1(splane.str() + "UCount", 128, 0., 128.);
	    }
	  hnu->Fill(s_u);
	}
      if (vb.count() > 0)
	{
	  int f = 1000, l = -1;
	  std::vector<int> first, last;
	  for (int i = 0; i < 128; i++)
	    {
	      if (vb[i] && f > 128)
		f = i;
	      if (vb[i] == 0 && f < 128 && l < 0)
		{
		  l = i - 1;
		  first.push_back(f);
		  last.push_back(l);
		  f = 1000;
		  l = -1;
		}
	    }
	  if (f < 128 && l < 0)
	    {
	      first.push_back(f);
	      last.push_back(127);
	    }

	  int idx = -1;
	  for (int i = 0; i < first.size(); i++)
	    if ((last[i] - first[i] + 1) > s_v)
	      {
		s_v = (last[i] - first[i] + 1);
		idx = i;
	      }
	  //std::cout<<"V IDX "<<first.size()<<idx<<std::endl;
	  x_v = (first[idx] + last[idx]) / 2.;
	  //std::cout<<" V :"<<f<<":"<<l<<":"<<x_v<<"=> "<<vb<<std::endl;
	  TH1 *hnv = _rh->GetTH1(splane.str() + "VCount");

	  if (hnv == NULL)
	    {
	      hnv = _rh->BookTH1(splane.str() + "VCount", 128, 0., 128.);
	    }
	  hnv->Fill(s_v);
	}
      if (wb.count() > 0)
	{
	  int f = 1000, l = -1;
	  std::vector<int> first, last;
	  for (int i = 0; i < 128; i++)
	    {
	      if (wb[i] && f > 128)
		f = i;
	      if (wb[i] == 0 && f < 128 && l < 0)
		{
		  l = i - 1;
		  first.push_back(f);
		  last.push_back(l);
		  f = 1000;
		  l = -1;
		}
	    }
	  if (f < 128 && l < 0)
	    {
	      first.push_back(f);
	      last.push_back(127);
	    }

	  int idx = -1;
	  for (int i = 0; i < first.size(); i++)
	    if ((last[i] - first[i] + 1) > s_w)
	      {
		s_w = (last[i] - first[i] + 1);
		idx = i;
	      }
	  //std::cout<<"W IDX "<<idx<<std::endl;
	  x_w = (first[idx] + last[idx]) / 2.;
	  //std::cout<<" W :"<<f<<":"<<l<<":"<<x_w<<"=> "<<wb<<std::endl;
	  TH1 *hnw = _rh->GetTH1(splane.str() + "WCount");

	  if (hnw == NULL)
	    {
	      hnw = _rh->BookTH1(splane.str() + "WCount", 128, 0., 128.);
	    }
	  hnw->Fill(s_w);
	}
      if (ub.count() > 0 || vb.count() > 0 || wb.count() > 0)
	{
	  //std::cout<<ub.count()<<":"<<vb.count()<<":"<<wb.count()<<std::endl;
	  //std::cout<<"COUCOU"<<_maxTime-x.first<<std::endl;
	  //getchar();
	}
      float X_uv, X_uw, X_vw, X_uvw;
      float Y_uv, Y_uw, Y_vw, Y_uvw;
      bool buv = false, buw = false, bvw = false, buvw = false;
      TH2 *huvwo = _rh->GetTH2(splane.str() + "UVWOr");

      if (huvwo == NULL)
	{

	  huvwo = _rh->BookTH2(splane.str() + "UVWOr", 100, -60., 60., 100, 0., 120.);
	}

      if (x_u >= 0 && x_v >= 0)
	{
	  float a_v = -1. * T30, b_v = (128 - x_v) * TSTEP;
	  float X = (x_u - 64) * XSTEP;
	  float Y = a_v * X + b_v;
	  X_uv = X;
	  Y_uv = Y;
	  TH2 *huv = _rh->GetTH2(splane.str() + "UV");
	  TH2 *huvw = _rh->GetTH2(splane.str() + "UVW");

	  if (huv == NULL)
	    {
	      huv = _rh->BookTH2(splane.str() + "UV", 100, -60., 60., 100, 0., 120.);
	      huvw = _rh->BookTH2(splane.str() + "UVW", 100, -60., 60., 100, 0., 120.);
	    }
	  huv->Fill(X, Y);
	  if (!all)
	    hst->Fill(11.);
	  if (x_w > 0)
	    {
	      buvw = true;
	      X_uvw = X;
	      Y_uvw = Y;
	      huvw->Fill(X, Y);
	      if (!all)
		hst->Fill(21.);

	      TH1 *hn3u = _rh->GetTH1(splane.str() + "U3Count");
	      TH1 *hn3v = _rh->GetTH1(splane.str() + "V3Count");
	      TH1 *hn3w = _rh->GetTH1(splane.str() + "W3Count");

	      if (hn3w == NULL)
		{
		  hn3w = _rh->BookTH1(splane.str() + "W3Count", 128, 0., 128.);
		  hn3v = _rh->BookTH1(splane.str() + "V3Count", 128, 0., 128.);
		  hn3u = _rh->BookTH1(splane.str() + "U3Count", 128, 0., 128.);
		}
	      hn3u->Fill(s_u);
	      hn3v->Fill(s_v);
	      hn3w->Fill(s_w);
	    }
	  buv = true;
	}
      if (x_u >= 0 && x_w >= 0)
	{
	  float a_w = 1. * T30, b_w = x_w * TSTEP;
	  float X = (x_u - 64) * XSTEP;
	  float Y = a_w * X + b_w;

	  X_uw = X;
	  Y_uw = Y;
	  //std::cout<<x_u<<":"<<x_w<<"=>"<<X<<","<<Y<<std::endl;
	  //getchar();

	  TH2 *huw = _rh->GetTH2(splane.str() + "UW");

	  if (huw == NULL)
	    {
	      huw = _rh->BookTH2(splane.str() + "UW", 100, -60., 60., 100, 0., 120.);
	    }
	  huw->Fill(X, Y);
	  if (!all)
	    hst->Fill(12.);
	  buw = true;
	}
      if (x_w >= 0 && x_v >= 0)
	{
	  float a_w = T30, a_v = -1. * T30, b_w = x_w * TSTEP, b_v = (128 - x_v) * TSTEP;
	  float X = (b_v - b_w) / (a_w - a_v);
	  float Y = a_w * X + b_w;

	  X_vw = X;
	  Y_vw = Y;
	  TH2 *hvw = _rh->GetTH2(splane.str() + "VW");

	  if (hvw == NULL)
	    {
	      hvw = _rh->BookTH2(splane.str() + "VW", 100, -60., 60., 100, 0., 120.);
	    }
	  hvw->Fill(X, Y);
	  if (!all)
	    hst->Fill(13.);
	  bvw = true;
	}
      if (buv || buw || bvw)
	if (!all)
	  hst->Fill(22.);
      if (buvw)
	huvwo->Fill(X_uvw, Y_uvw);
      else if (buv)
	huvwo->Fill(X_uv, Y_uv);
      else if (buw)
	huvwo->Fill(X_uw, Y_uw);
      else if (bvw)
	huvwo->Fill(X_vw, Y_vw);
    }
}

int32_t binaryreader::TPrincipalComponents(double result[21], float zmin, float zmax)
{
  double resultc[21];
  uint32_t nh = 0;
  double xb = 0, yb = 0, zb = 0;
  double wt = 0.;

  double fp = DBL_MAX;
  double lp = -DBL_MAX;
  double fx = DBL_MAX;
  double lx = -DBL_MAX;
  double fy = DBL_MAX;
  double ly = -DBL_MAX;
  TPrincipal tp(3, "D");
  double xp[3];
  memset(result, 0, 21 * sizeof(double));
  //INFO_PRINT("%d vector size\n",v.size());
  for (auto it = _vPoints.begin(); it != _vPoints.end(); it++)
    {
      ROOT::Math::XYZPoint &iht = (*it);
      if (iht.Z() < zmin)
	continue;
      if (iht.Z() > zmax)
	continue;
      //INFO_PRINT("%x %d %d \n",iht,iht.I(),iht.J());
      //INFO_PRINT("%f %f \n",iht.x(),iht.y());
      //INFO_PRINT("%f %f \n",iht.X(),iht.Y());
      double w = 1.;
      xb += iht.X() * w;
      yb += iht.Y() * w;
      zb += iht.Z() * w;
      wt += w;
      nh++;
      xp[0] = iht.X();
      xp[1] = iht.Y();
      xp[2] = iht.Z();
      //printf("XP  %f %f %f \n",xp[0],xp[1],xp[2]);
      tp.AddRow(xp);
    }

  if (nh < 2)
    return -1;
  tp.MakePrincipals();
  // store barycenter
  const TVectorD *fvb = tp.GetMeanValues();
  // printf("barycentre %f %f %f \n \t %f %f %f \n", xb/wt,yb/wt,zb/wt,(*fvb)[0],(*fvb)[1],(*fvb)[2]);
  result[0] = (*fvb)[0];
  result[1] = (*fvb)[1];
  result[2] = (*fvb)[2];

  const TVectorD *fva = tp.GetEigenValues();
  result[3] = (*fva)[2];
  result[4] = (*fva)[1];
  result[5] = (*fva)[0];

  //fva->Print();
  // printf("eigen results %g %g %g \n",(*fva)[0],(*fva)[1],(*fva)[2]);

  //tp.Print("MSEV");

  // store principal axis
  const TMatrixD *fvv = tp.GetEigenVectors();
  //Matrix<double,3,3> vv=eigensolver.eigenvectors();
  result[6] = (*fvv)(0, 2);
  result[7] = (*fvv)(1, 2);
  result[8] = (*fvv)(2, 2);
  // printf("eigen vector results %g %g %g \n",result[7],result[7],result[8]);
  result[9] = (*fvv)(0, 1);

  result[10] = (*fvv)(1, 1);
  result[11] = (*fvv)(2, 1);
  // printf("eigen vector results %g %g %g \n",result[9],result[10],result[11]);
  result[12] = (*fvv)(0, 0);
  result[13] = (*fvv)(1, 0);
  result[14] = (*fvv)(2, 0);
  // printf("eigen vector results %g %g %g \n",result[12],result[13],result[14]);
  // Store First and last Z
  result[15] = fp;
  result[16] = lp;
  result[17] = fx;
  result[18] = lx;
  result[19] = fy;
  result[20] = ly;
  //INFO_PRINT("=11\n");
  //getchar();
  /*
    if (TCPC==NULL)
    {
    TCPC=new TCanvas("TCPC","tChits1",1300,600);
    TCPC->Modified();
    TCPC->Draw();

    }
    TCPC->cd(0);
    tp.Test();
    TCPC->Modified();
    TCPC->Draw();
    TCPC->Update();
  */
  return 0;
}

void binaryreader::createTrees(std::string s)
{

  treeFile_ = new TFile(s.c_str(), "recreate");
  treeFile_->cd();

  tEvents_ = new TTree("events", "Events");
  tEvents_->SetAutoSave(50000000);
  
  tEvents_->Branch("bcid", &_bxdif, "bcid/l");
  tEvents_->Branch("run", &_run, "run/i");
  tEvents_->Branch("event", &_event, "event/i ");
  tEvents_->Branch("bsplan", &_bsplanes, "bsplan/l");
  tEvents_->Branch("top_hit", &_t_h,"top_hit/b");
  tEvents_->Branch("bot_hit", &_b_h,"bot_hit/b");
  tEvents_->Branch("all_hit", &_a_h,"bot_hit/b");
  tEvents_->Branch("t_x", &_t_x, "t_x[3]/D");
  tEvents_->Branch("t_v", &_t_v, "t_v[3]/D");
  tEvents_->Branch("b_x", &_t_x, "b_x[3]/D");
  tEvents_->Branch("b_v", &_t_v, "b_v[3]/D");
  tEvents_->Branch("a_x", &_t_x, "a_x[3]/D");
  tEvents_->Branch("a_v", &_t_v, "a_v[3]/D");
  tEvents_->Branch("t_c2", &_t_c2, "t_c2/D");
  tEvents_->Branch("b_c2", &_b_c2, "b_c2/D");
  tEvents_->Branch("dist", &_dist, "dist/D");
  tEvents_->Branch("cos_th", &_cos_th, "cos_th/D");
  tEvents_->Branch("xcross", &_xcross, "xcross/D");
  tEvents_->Branch("ycross", &_ycross, "ycross/D");
  tEvents_->Branch("zcross", &_zcross, "zcross/D");
  tEvents_->Branch("rd3", &_rd3, "rd3/D");
  tEvents_->Branch("probd3", &_probd3, "probd3/D");
  

  std::cout << " create Trees" << std::endl;
}
void binaryreader::closeTrees()
{
  if (tEvents_!=0)
    {
  treeFile_->cd();
  tEvents_->Write();
  treeFile_->ls();
  treeFile_->Close();
    }
}
void binaryreader::scurveAnalysis(rbEvent *e)
{

  if (e->seuil()==0) return;
  //std::cout<<"Event "<<_event<<" GTC"<<_gtc<<" Vth set "<<e->seuil()<<std::endl;
  //fflush(stdout);
  for (int id = 0; id < MAXDIF; id++)
    if (e->frameCount(id))
      {
	std::stringstream sraw1;
	std::stringstream sraw2;
	std::stringstream sraw3;
	sraw1 << "/gric/B01SCURVE" << std::hex << id << std::dec << "/";
	sraw2 << "/gric/B10SCURVE" << std::hex << id << std::dec << "/";
	sraw3 << "/gric/B11SCURVE" << std::hex << id << std::dec << "/";
	
	TH1 *hp1 = _rh->GetTH1(sraw1.str() + "Padc1");
	if (hp1 == NULL)
	  {
	    for (int i=0;i<64;i++)
	      {
		std::stringstream srpc("");
		srpc<<sraw1.str()<<"Padc"<<i;
		TH1* hpc = _rh->BookTH1(srpc.str(),1024,0.,1024);
		srpc.str( std::string() );srpc.clear();
		srpc<<sraw2.str()<<"Padc"<<i;
		TH1* hpc2 = _rh->BookTH1(srpc.str(),1024,0.,1024);
		srpc.str( std::string() );srpc.clear();
		srpc<<sraw3.str()<<"Padc"<<i;
		TH1* hpc3 = _rh->BookTH1(srpc.str(),1024,0.,1024);
	      }
	  }

	for (int j = 0; j < e->frameCount(id); j++)
	  {
	    uint32_t idx = e->iPtr(id, j);

	    if (e->bcid(idx) < 4)
	      continue;
	    for (int k = 0; k < 64; k++)
	      {
		if (e->pad0(idx, k) && !e->pad1(idx,k) )
		  {
		    std::stringstream srpc("");
		    srpc<<sraw1.str()<<"Padc"<<k;

		    TH1* hpc= _rh->GetTH1(srpc.str());
		    //std::cout<<srpc.str()<<std::endl;
		    hpc->Fill(e->seuil()*1.);

		  }
		if (e->pad1(idx, k) && !e->pad0(idx,k))
		  {
		    std::stringstream srpc("");
		    srpc<<sraw2.str()<<"Padc"<<k;

		    TH1* hpc= _rh->GetTH1(srpc.str());
		    //std::cout<<srpc.str()<<std::endl;
		    hpc->Fill(e->seuil()*1.);

		  }
		if (e->pad0(idx, k) && e->pad1(idx, k))
		  {
		    std::stringstream srpc("");
		    srpc<<sraw3.str()<<"Padc"<<k;

		    TH1* hpc= _rh->GetTH1(srpc.str());
		    //std::cout<<srpc.str()<<std::endl;
		    hpc->Fill(e->seuil()*1.);

		  }
	      }
	  }
      }

}

extern "C"
{
  // loadDHCALAnalyzer function creates new LowPassDHCALAnalyzer object and returns it.
  rbProcessor *loadProcessor(void)
  {
    return (new binaryreader);
  }
  // The deleteDHCALAnalyzer function deletes the LowPassDHCALAnalyzer that is passed
  // to it.  This isn't a very safe function, since there's no
  // way to ensure that the object provided is indeed a LowPassDHCALAnalyzer.
  void deleteProcessor(rbProcessor *obj)
  {
    delete obj;
  }
}
