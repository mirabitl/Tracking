#include "tricotreader.hh"
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

#include <stdint.h>
#include <math.h>

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
void tricotreader::info(){std::cout<<"JE SUIS UN TRICOTREADER"<<std::endl;}
tricotreader::tricotreader() : _run(0), _started(false), _fdOut(-1), _totalSize(0), _event(0) {}
void tricotreader::init(uint32_t run)
{
  _run = run;
  _event = 0;
  _started = true;
  _rh = DCHistogramHandler::instance();

}
void tricotreader::loadParameters(Json::Value params)
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
	  //z[plane] = ch["z"].asFloat();
	  //togric[plane]=pl;
	  _plinfo.insert(p);
	  printf("%s GRIC %d Plane  %d is at 0 cm from floor \n", _plinfo[pl]["name"].asString().c_str(), pl, plane);
	}
      //getchar();
    }
}
void tricotreader::end(uint32_t run)
{
  //this->closeTrees();
  _started = false;
}

void tricotreader::processRunHeader(std::vector<uint32_t> header)
{
}
void tricotreader::processEvent(rbEvent *e)
{
  uint8_t u[16], v[16], w[16];
  if (!_started)
    return;
  if (e->gtc()%100!=100)
    printf("BR => %d %d %d %d \n",e->run(),e->event(),e->gtc(),e->seuil());
  _event = e->gtc();
  _run = e->run();
  _gtc=e->gtc();
  
  if (e->seuil()!=0)
    this->scurveAnalysis(e);

  this->fillTimeMap(e);

  return;
}
void tricotreader::fillTimeMap(rbEvent *e)
{
  _timeMap.clear();
  std::stringstream sraw;
  sraw << "/gric/";
  TH1 *hfc = _rh->GetTH1(sraw.str() + "FrameCount");
  TH1 *hftm = _rh->GetTH1(sraw.str() + "MaxTime");
  TH1 *hfound = _rh->GetTH1(sraw.str() + "NFrameInTime");


  if (hfc == NULL)
    {

      hfc = _rh->BookTH1(sraw.str() + "FrameCount", 255, 0., 255.);
      hftm = _rh->BookTH1(sraw.str() + "MaxTime", 65536, 0., 2.);
      hfound = _rh->BookTH1(sraw.str() + "NFrameInTime", 30, -0.1, 29.9);
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
  TH1 *hp1s = _rh->GetTH1(sraw1.str() + "Pad1Selected");
  TH1 *hp1 = _rh->GetTH1(sraw1.str() + "Pad1");
  TH1 *hft = _rh->GetTH1(sraw1.str() + "FrameTime");
  TH1 *htt = _rh->GetTH1(sraw1.str() + "Time2Trigger");
  TH1 *hfc = _rh->GetTH1(sraw1.str() + "FrameCount");
  if (hp1 == NULL)
    {
      hp1 = _rh->BookTH1(sraw1.str() + "Pad1", 64, 0., 64.);
      hp1s = _rh->BookTH1(sraw1.str() + "Pad1Selected", 64, 0., 64.);
      hft = _rh->BookTH1(sraw1.str() + "FrameTime", 65536, 0., 2.);
      htt = _rh->BookTH1(sraw1.str() + "Time2Trigger", 1000, 0., 1000.);
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
	    bool selected=false;
	    if (e->trgBcid(id)!=0)
	      {
		htt->Fill(e->trgBcid(id)-e->bcid(idx));
		selected=(e->trgBcid(id)-e->bcid(idx))<15;
	      }
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
		  if (selected)  hp1s->Fill(k * 1.);
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
		/*
		if (itmmm != _timeMap.end())
		  itmmm->second.push_back(idx);
		if (itmpp != _timeMap.end())
		  itmpp->second.push_back(idx);
		*/
	      }

	  }
      }
  // Ask at least 2 frames selected in time

  for (auto it=_timeMap.begin();it!=_timeMap.end();)
    {
      auto x=(*it);
      uint32_t nh=0;
      std::bitset<3> dirs(0);
      for (auto y : x.second)
	{
	      // ID
	      uint32_t id = y / MAXFRAME / FSIZE;
	      if ((e->trgBcid(id)-e->bcid(y))<15)
		{
		  uint32_t dir=((id>>4)&0xF)-1;
		  dirs.set(dir);
		  nh++;
		}
	}

      
      if (nh!=0)
	{
	std::cout<<"Map size "<<x.first<<" nh "<<nh<<" dirs "<<dirs<<std::endl;
	}
      if (dirs.count()>1)
	{
	  for (auto y : x.second)
	    {
	      // ID
	      uint32_t id = y / MAXFRAME / FSIZE;
	      uint32_t dir=((id>>4)&0xF);
	      uint32_t dif=id&0xF;
	      std::stringstream sraw1;

	      sraw1 << "/gric/DIR"<<dir<< std::dec << "/";
	      TH1 *hp1r = _rh->GetTH1(sraw1.str() + "Pad1Reco");
	      if (hp1r==NULL)
		{
		  hp1r = _rh->BookTH1(sraw1.str() + "Pad1Reco", 192, 0., 192.);
 
		}
	      for (int k = 0; k < 64; k++)
		if (e->pad0(y, k) || e->pad1(y, k))
		{
		  hp1r->Fill(k * 1.+(dif-1)*64);

		}
	    }
	  this->buildPlaneHits(e, x.second);
	  //
	  it++;
	}
      else
	{
	  it->second.clear();
	  _timeMap.erase(it++);
	}
    }

  hfound->Fill(_timeMap.size()*1.);
  //if (_timeMap.size()>1)
  //  getchar();
}
void tricotreader::buildPlaneHits(rbEvent *e, std::vector<uint32_t> &hits)

{


  std::stringstream splane;

  splane << "/gric/PLANE/";
  std::bitset<384> ub;
  ub.reset();
  std::bitset<384> vb;
  vb.reset();
  std::bitset<384> wb;
  wb.reset();
  std::bitset<384> ub2;
  ub2.reset();
  std::bitset<384> vb2;
  vb2.reset();
  std::bitset<384> wb2;
  wb2.reset();
  std::bitset<384> ub3;
  ub3.reset();
  std::bitset<384> vb3;
  vb3.reset();
  std::bitset<384> wb3;
  wb3.reset();

  //  std::cout<<"Hit Plane before "<<plane<<" BS "<<_hplanes<<std::endl;
  uint32_t igplane = 0;
  for (auto y : hits)
    {
      // Select the plane
      uint32_t ig = y / MAXFRAME / FSIZE;
      uint16_t igpl = (ig >> 4) & 0xF;
      uint32_t idif=(ig&0xF)-1;
      igplane = igpl;
      uint16_t dir = (igpl & 0xF);
      bool udir = (dir == 1);
      bool vdir = (dir == 3);
      bool wdir = (dir == 2);
      int sshift = idif*64;
      // if (dir == 3)
      //  	sshift = (2-idif)*64;
      for (int k = 0; k < 64; k++)
	{
	  if (e->pad0(y, k) || e->pad1(y, k))
	    {
	      if (udir)
		ub.set(k + sshift);
	      if (vdir)
		vb.set(63-k + sshift);
	      if (wdir)
		wb.set(k + sshift);
	    }
	  if (e->pad0(y, k) && !e->pad1(y, k))
	    {
	      if (udir)
		ub2.set(k + sshift);
	      if (vdir)
		vb2.set(63-k + sshift);
	      if (wdir)
		wb2.set(k + sshift);
	    }
	  if (e->pad0(y, k) && e->pad1(y, k))
	    {
	      if (udir)
		ub3.set(k + sshift);
	      if (vdir)
		vb3.set(63-k + sshift);
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
      for (int i = 0; i < 384; i++)
	{
	  if (ub[i] && f > 384)
	    f = i;
	  if (ub[i] == 0 && f < 384 && l < 0)
	    {
	      l = i - 1;
	      first.push_back(f);
	      last.push_back(l);
	      f = 1000;
	      l = -1;
	    }
	}
      if (f < 384 && l < 0)
	{
	  first.push_back(f);
	  last.push_back(255);
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
	  hnu = _rh->BookTH1(splane.str() + "UCount", 64, 0., 64.);
	}
      hnu->Fill(s_u);
    }
  if (vb.count() > 0)
    {
      int f = 1000, l = -1;
      std::vector<int> first, last;
      for (int i = 0; i < 384; i++)
	{
	  if (vb[i] && f > 384)
	    f = i;

	  if (vb[i] == 0 && f < 384 && l < 0)
	    {
	      l = i - 1;
	      first.push_back(f);
	      last.push_back(l);
	      f = 1000;
	      l = -1;
	    }
	}
      if (f < 384 && l < 0)
	{
	  first.push_back(f);
	  last.push_back(255);
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
	  hnv = _rh->BookTH1(splane.str() + "VCount", 64, 0., 64.);
	}
      hnv->Fill(s_v);
    }
  if (wb.count() > 0)
    {
      int f = 1000, l = -1;
      std::vector<int> first, last;
      for (int i = 0; i < 384; i++)
	{
	  if (wb[i] && f > 384)
	    f = i;

	  if (wb[i] == 0 && f < 384 && l < 0)
	    {
	      l = i - 1;
	      first.push_back(f);
	      last.push_back(l);
	      f = 1000;
	      l = -1;
	    }
	}
      if (f < 384 && l < 0)
	{
	  first.push_back(f);
	  last.push_back(255);
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
	  hnw = _rh->BookTH1(splane.str() + "WCount", 64, 0., 64.);
	}
      hnw->Fill(s_w);
    }
  if (ub.count() > 0 || vb.count() > 0 || wb.count() > 0)
    {
      // std::cout<<ub.count()<<":"<<vb.count()<<":"<<wb.count()<<std::endl;
      // std::cout<<"COUCOU"<<x_u<<" "<<x_v<<" "<<x_w<<std::endl;
      // getchar();
    }
  float X_uv, X_uw, X_vw, X_uvw;
  float Y_uv, Y_uw, Y_vw, Y_uvw;
  bool buv = false, buw = false, bvw = false, buvw = false;
  TH2 *huvwo = _rh->GetTH2(splane.str() + "UVWOr");

  if (huvwo == NULL)
    {

      huvwo = _rh->BookTH2(splane.str() + "UVWOr", 400, -200., 200., 400,0.,200.);
    }

  if (x_u >= 0 && x_v >= 0)
    {
      float a_v = -1. * T30, b_v = 48. + x_v * TSTEP;
      float X = (x_u) * XSTEP;
      float Y = a_v * X + b_v;
      X_uv = X;
      Y_uv = Y;
      TH2 *huv = _rh->GetTH2(splane.str() + "UV");
      TH2 *huvw = _rh->GetTH2(splane.str() + "UVW");

      if (huv == NULL)
	{
	  huv = _rh->BookTH2(splane.str() + "UV",  400, -200., 200., 400,0.,200.);
	  huvw = _rh->BookTH2(splane.str() + "UVW", 400, -200., 200., 400,0.,200.);
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
	      hn3w = _rh->BookTH1(splane.str() + "W3Count", 64, 0., 64.);
	      hn3v = _rh->BookTH1(splane.str() + "V3Count", 64, 0., 64.);
	      hn3u = _rh->BookTH1(splane.str() + "U3Count", 64, 0., 64.);
	    }
	  hn3u->Fill(s_u);
	  hn3v->Fill(s_v);
	  hn3w->Fill(s_w);
	}
      buv = true;
    }
  if (x_u >= 0 && x_w >= 0)
    {
      float a_w = 1. * T30, b_w = 16. +x_w * TSTEP;
      float X = (x_u) * XSTEP;
      float Y = a_w * X + b_w;

      X_uw = X;
      Y_uw = Y;
      //std::cout<<x_u<<":"<<x_w<<"=>"<<X<<","<<Y<<std::endl;
      //getchar();

      TH2 *huw = _rh->GetTH2(splane.str() + "UW");

      if (huw == NULL)
	{
	  huw = _rh->BookTH2(splane.str() + "UW", 400, -200., 200., 400,0.,200.);
	}
      huw->Fill(X, Y);

      buw = true;
    }
  if (x_w >= 0 && x_v >= 0)
    {
      float a_w = 1.* T30, a_v = -1. * T30, b_v = 48. + x_v * TSTEP, b_w = 16. +x_w * TSTEP;
      float X = (b_v - b_w) / (a_w - a_v);
      float Y = a_v * X + b_v;

      X_vw = X;
      Y_vw = Y;
      TH2 *hvw = _rh->GetTH2(splane.str() + "VW");

      if (hvw == NULL)
	{
	  hvw = _rh->BookTH2(splane.str() + "VW",  400, -200., 200., 400,0.,200.);
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

	      float a_v = -1. * T30, b_v = (384 - xv) * TSTEP;
	      float Xuv = (xu - 64) * XSTEP;
	      float Yuv = a_v * Xuv + b_v;
	      float a_w = 1. * T30, b_w = xw * TSTEP;
	      float Xuw = (x_u - 64) * XSTEP;
	      float Yuw = a_w * Xuw + b_w;
	      a_w = T30;
	      a_v = -1. * T30, b_w = xw * TSTEP;
	      b_v = (384 - xv) * TSTEP;
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
	    }
    }
  else if (buv)
    {
      for (auto xu : vx_u)
	for (auto xv : vx_v)
	  {

	    float a_v = -1. * T30, b_v = (384 - xv) * TSTEP;
	    float Xuv = (xu - 64) * XSTEP;
	    float Yuv = a_v * Xuv + b_v;
	    //recoPoint pg(Xuv,Yuv,z[plane]);
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

	    float a_w = 1. * T30, b_w = xw * TSTEP;
	    float Xuw = (x_u - 64) * XSTEP;
	    float Yuw = a_w * Xuw + b_w;

	    //getchar();
	  }
    }
  else if (bvw)
    {

      for (auto xv : vx_v)
	for (auto xw : vx_w)
	  {

	    float a_w = T30, a_v = -1. * T30, b_w = xw * TSTEP, b_v = (384 - xv) * TSTEP;
	    float Xvw = (b_v - b_w) / (a_w - a_v);
	    float Yvw = a_w * Xvw + b_w;

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

void tricotreader::scurveAnalysis(rbEvent *e)
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
	uint64_t lb0=0,lb1=0,lb2=0;
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
		    lb0 |=(1<<k);

		  }
		if (e->pad1(idx, k) && !e->pad0(idx,k))
		  {
		    std::stringstream srpc("");
		    srpc<<sraw2.str()<<"Padc"<<k;

		    TH1* hpc= _rh->GetTH1(srpc.str());
		    //std::cout<<srpc.str()<<std::endl;
		    hpc->Fill(e->seuil()*1.);
		    lb1 |=(1<<k);
		  }
		if (e->pad0(idx, k) && e->pad1(idx, k))
		  {
		    std::stringstream srpc("");
		    srpc<<sraw3.str()<<"Padc"<<k;

		    TH1* hpc= _rh->GetTH1(srpc.str());
		    //std::cout<<srpc.str()<<std::endl;
		    hpc->Fill(e->seuil()*1.);
		    lb2 |=(1<<k);
		  }
	      }
	  }
	for (int i=0;i<64;i++)
	  {
	    if (lb0==0 &&lb1==0 && lb2==0) continue;
	    if ((lb0>>i)&1==1)
	      {
		std::stringstream srpc("");
		srpc<<sraw1.str()<<"Pade"<<i;
		TH1* hpe= _rh->GetTH1(srpc.str());
		if (hpe==0)
		   hpe = _rh->BookTH1(srpc.str(),1024,0.,1024);
		hpe->Fill(e->seuil()*1.);
	      }
	    if ((lb1>>i)&1==1)
	      {
		std::stringstream srpc("");
		srpc<<sraw2.str()<<"Pade"<<i;
		TH1* hpe= _rh->GetTH1(srpc.str());
		if (hpe==0)
		   hpe = _rh->BookTH1(srpc.str(),1024,0.,1024);
		hpe->Fill(e->seuil()*1.);
	      }
	    if ((lb2>>i)&1==1)
	      {
		std::stringstream srpc("");
		srpc<<sraw3.str()<<"Pade"<<i;
		TH1* hpe= _rh->GetTH1(srpc.str());
		if (hpe==0)
		   hpe = _rh->BookTH1(srpc.str(),1024,0.,1024);
		hpe->Fill(e->seuil()*1.);
	      }

	  }
	  

      }

}

extern "C"
{
  // loadDHCALAnalyzer function creates new LowPassDHCALAnalyzer object and returns it.
  rbProcessor *loadProcessor(void)
  {
    return (new tricotreader);
  }
  // The deleteDHCALAnalyzer function deletes the LowPassDHCALAnalyzer that is passed
  // to it.  This isn't a very safe function, since there's no
  // way to ensure that the object provided is indeed a LowPassDHCALAnalyzer.
  void deleteProcessor(rbProcessor *obj)
  {
    delete obj;
  }
}
