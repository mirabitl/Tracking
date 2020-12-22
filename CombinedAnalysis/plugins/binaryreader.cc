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
#include "TPrincipal.h"
#include "HoughLocal.hh"

using namespace lydaq;
using namespace Lmana;
binaryreader::binaryreader() : _run(0),_started(false),_fdOut(-1),_totalSize(0),_event(0),tEvents_(NULL) {}
void binaryreader::init(uint32_t run)
{
  _run=run; 
  _event=0;
  _started=true;
  _rh=DCHistogramHandler::instance();
}
void binaryreader::loadParameters(Json::Value params)
{
  
  if (params.isMember("general"))
    std::cout<<"DIRECTORY "<<params["general"]["directory"].asString()<<std::endl;
  _geoRoot=params;
  _geo=new jsonGeo();
  _geo->init(_geoRoot);
  //getchar();
}
void binaryreader::end(uint32_t run)
{
  this->closeTrees();

  _started=false;

}

void binaryreader::processRunHeader(std::vector<uint32_t> header)
{
}
static TCanvas* TCHits=NULL;
void binaryreader::processCoincidence(rbEvent* e,uint32_t ibc)
{
  memset(&_fevt,0,sizeof(struct FullEventTree));
  _fevt.bc=ibc;
  _fevt.run=e->run();
  _fevt.gtc=e->gtc();
  _fevt.event=e->event();

  if (tEvents_==NULL)
    {
      std::stringstream ss;
      if (_geoRoot["general"]["noise"].asUInt()==0)
	ss<<_geo->general()["directory"].asString()<<"/tree"<<e->run()<<"_"<<e->gtc()<<".root";
      else
	ss<<_geo->general()["directory"].asString()<<"/Noisetree"<<e->run()<<"_"<<e->gtc()<<".root";
      this->createTrees(ss.str());			    
    }
  
  if (ibc<30) return;
  TH2* hzx=_rh->GetTH2("/gric/ZX");
  TH2* hzy=_rh->GetTH2("/gric/ZY");
  TH2* hxy=_rh->GetTH2("/gric/XY");
  TH2* hxyt=_rh->GetTH2("/gric/XYT");
  TH2* hxyf=_rh->GetTH2("/gric/XYF");
  TH2* hxyf14=_rh->GetTH2("/gric/XYF14");
  TH2* hxyf15=_rh->GetTH2("/gric/XYF15");
  TH1* hcount=_rh->GetTH1("/gric/Count");
  TH1* hdt=_rh->GetTH1("/gric/dt");
  TH1* hch=_rh->GetTH1("/gric/tdcchannel");
  TH1* hinti=_rh->GetTH1("/gric/InTime");
  TH1* htkchi2=_rh->GetTH1("/gric/tkchi2");
  TH1* htkax=_rh->GetTH1("/gric/tkax");
  TH1* htkay=_rh->GetTH1("/gric/tkay");
  if (hzx==NULL)
    {
      htkchi2=_rh->BookTH1("/gric/tkchi2",300,0.,1.01);
      htkax=_rh->BookTH1("/gric/tkax",300,-2.,2.);
      htkay=_rh->BookTH1("/gric/tkay",300,-2.,2.);
      hcount=_rh->BookTH1("/gric/Count",30,0.1,30.1);
      hdt=_rh->BookTH1("/gric/dt",500,-250.,250.);
      hch=_rh->BookTH1("/gric/tdcchannel",70,0.,70.);
      hinti=_rh->BookTH1("/gric/InTime",70,-0.1,69.9);
      hzx=_rh->BookTH2("/gric/ZX",400,0.,200.,100.,0.,100.);
      hzy=_rh->BookTH2("/gric/ZY",400,0.,200.,80.,0.,80.);
      hxy=_rh->BookTH2("/gric/XY",100,0.,50.,100.,0.,100.);
      hxyt=_rh->BookTH2("/gric/XYT",100,0.,50.,100.,0.,100.);
      hxyf=_rh->BookTH2("/gric/XYF",100,0.,50.,100.,0.,100.);
      hxyf14=_rh->BookTH2("/gric/XYF14",100,0.,50.,100.,0.,100.);
      hxyf15=_rh->BookTH2("/gric/XYF15",100,0.,50.,100.,0.,100.);
	    }
  hzx->Reset();
  hzy->Reset();
  hcount->Fill(1.);
  _vPoints.clear();
  _vPads.clear();
  _vHRCl.clear();
  hcount->Fill(1.);
  for (auto it=e->tFrame()[ibc].begin();it!=e->tFrame()[ibc].end();it++)
    {
      //std::cout<<it->first->getID()<<" "<<it->first->getFrameTimeToTrigger(it->second)<<std::endl;
      uint32_t ifra=it->second;
      sdhcal::PMRPtr* d=it->first;
      uint32_t chid=_geo->difInfo(d->getID()).chamber;

      std::stringstream splane;

      
      splane<<"/pmr/PLANE"<<chid<<"/";
      TH1* hfcs=_rh->GetTH1(splane.str()+"FrameCountSel");
      //      TH1* hax=_rh->GetTH1(splane.str()+"ax");
      if (hfcs==NULL)
	{
	  hfcs=_rh->BookTH1(splane.str()+"FrameCountSel",65,0.,65.0);
	}
      hcount->Fill(chid+1.);
      std::bitset<64> ph;
      ph.reset();
      for (int ipad=0;ipad<64;ipad++)
	{
	  
	  if (d->getFrameLevel(ifra,ipad,0) || d->getFrameLevel(ifra,ipad,1))
	    {
	      ph.set(ipad,1);
	      //recoPoint pt(chid);
	      HR2Pad pt(chid,d->getID(),d->getFrameAsicHeader(ifra),ipad);
	      //pt.setPlan(chid);
	      _geo->convert(d->getID(),d->getFrameAsicHeader(ifra),ipad,&pt);
	      pt.SetZ(_geo->chamberInfo(chid).z0);
	      _vPads.push_back(pt);
	      _fevt.pad_dif[ _fevt.npad]=d->getID();
	      _fevt.pad_asic[ _fevt.npad]=d->getFrameAsicHeader(ifra);
	      _fevt.pad_channel[ _fevt.npad]=ipad;
	      _fevt.npad++;
	      if (_geoRoot["general"]["display"].asUInt()!=0)
	      //printf("Chamber %d %f \n",chid,_geo->chamberInfo(chid).z0);
		printf("Point %d %d %d %f %f %f \n",d->getID(),d->getFrameAsicHeader(ifra),ipad,pt.X(),pt.Y(),pt.Z());
	    }
	}
      hfcs->Fill(ph.count()*1.);
    }
  // Now build clusters
  for (auto it=_vPads.begin();it!=_vPads.end();it++)
	{
	  bool found=false;
	  for (auto ic=_vHRCl.begin();ic!=_vHRCl.end();ic++)
	    {
	      if (ic->isAdjacent((*it),3))
		{
		  ic->addPad((*it));
		  found=true;
		  break;
		}
	    }
	  if(!found)
	    {
	      Lmana::HR2Cluster c;
	      c.addPad((*it));
	      _vHRCl.push_back(c);
	    }
	}

      // Merge adjacent cluster
      bool merged=false;
      //printf("HRCL size %ld \n",_vHRCl.size());

      for (auto it=_vHRCl.begin();it!=_vHRCl.end();it++)
	{
	  for (auto jt=it+1;jt!=_vHRCl.end();)
	    {
	      bool adj=false;
	      for (int i=0;i<jt->size();i++)
		{
		  if (it->isAdjacent(jt->pad(i),3))
		    {adj=true;break;}
		}
	      if (adj)
		{
		  merged=true;
		  printf("Merigng cluster \n");
		  for (int i=0;i<jt->size();i++)
		    {
		      it->addPad(jt->pad(i));
		    }
		  _vHRCl.erase(jt);
		}
	      else
		++jt;
	    }
	}
      printf("HR2Cluster  size after %ld \n",_vHRCl.size());
      //getchar();
      for (auto x:_vHRCl)
	{
	  recoPoint pt(x.chamber());
	  float dx0=_geo->chamberInfo(x.chamber()).x0;
	  float dy0=_geo->chamberInfo(x.chamber()).y0;
	  pt.SetXYZ(x.X()+dx0,x.Y()+dy0,x.Z());
	  _vPoints.push_back(pt);
	  hzx->Fill(x.Z(),x.X());
	  hzy->Fill(x.Z(),x.Y());
	  _fevt.tel_x[_fevt.ntel]=x.X();
	  _fevt.tel_y[_fevt.ntel]=x.Y();
	  _fevt.tel_z[_fevt.ntel]=x.Z();
	  _fevt.ntel++;
	}
  ShowerParams isha;
  double ax=0,bx=0,ay=0,by=0;
  float zmin=0,zmax=200;
  int ier=TPrincipalComponents((double*) &isha,zmin,zmax);
  double* sx=isha.xm;
  double* sv=isha.l2;
  double z0=sx[2];
  double x0=sx[0];
  double y0=sx[1];
  double x1=sx[0]+sv[0];
  double y1=sx[1]+sv[1];
  double z1=sx[2]+sv[2];
  //double ax,ay,bx,by;
  if (ier==0)
    {
      ax=(x1-x0)/(z1-z0);
      bx=x1-ax*z1;
      ay=(y1-y0)/(z1-z0);
      by=y1-ay*z1;
    }
  
  top_tk.clear();

  top_tk.setDir(ax,ay,1.);
  top_tk.setOrig(bx,by,0);
  //printf("Shower ax %f ay %f bx %f by %f \n",ax,ay,bx,by);
  for (auto x=_vPoints.begin();x!=_vPoints.end();x++)
    {
      recoPoint &p=(*x);
      
      
      if (top_tk.distance(&p)<5.|| ier!=0)
	{
	  top_tk.addPoint(&p);
	}
    }
  // Ask at least 3 points
  if (top_tk.size()>=3) {
    top_tk.regression();
    top_tk.calculateChi2();
  }
  else
    return;
  if (top_tk.orig().X()>45) return;
  htkchi2->Fill(top_tk.pchi2());
  htkax->Fill(top_tk.dir().X());
  htkay->Fill(top_tk.dir().Y());
  if (top_tk.pchi2()<0.4) return;
  if (abs(top_tk.dir().X())>0.5) return;
  if (abs(top_tk.dir().Y())>0.5) return;
  //if (top_tk.plans()!=30) return;
  hcount->Fill(20.);
  uint32_t cnt[20];
  memset(cnt,0,20*sizeof(uint32_t));
  t12_tk.clear();
  for (auto x=_vPoints.begin();x!=_vPoints.end();x++)
    {
      hcount->Fill(22+(*x).plan());
      cnt[x->plan()]++;
      recoPoint &p=(*x);
      if (p.plan()<3)
	t12_tk.addPoint(&p);
    }
    t12_tk.regression();
    t12_tk.calculateChi2();
    if (cnt[3]>3 || cnt[4]>3) return;

    for (auto x=_vPoints.begin();x!=_vPoints.end();x++)
    {
      if (x->plan()<3) continue;
      ROOT::Math::XYZPoint pe=top_tk.extrapolate(x->Z());
      std::stringstream splane;

      
      splane<<"/pmr/ALIGN"<<x->plan()<<"/";
      TH1* hdx=_rh->GetTH1(splane.str()+"dx");
      TH1* hdy=_rh->GetTH1(splane.str()+"dy");
      if (hdx==NULL)
	{
	  hdx=_rh->BookTH1(splane.str()+"dx",500,-10.,10.0);
	  hdy=_rh->BookTH1(splane.str()+"dy",500,-10.,10.0);
	}
      hdx->Fill(pe.X()-x->X());
      hdy->Fill(pe.Y()-x->Y());
    }
    
  // printf(" TK ax %f ay %f bx %f by %f \n",
  // 	   top_tk.dir().X(),
  // 	   top_tk.dir().Y(),
  // 	   top_tk.orig().X(),
  // 	 top_tk.orig().Y());
    _fevt.tk_x[0]=top_tk.orig().X();
    _fevt.tk_x[1]=top_tk.orig().Y();
    _fevt.tk_x[2]=top_tk.orig().Z();
    _fevt.tk_v[0]=top_tk.dir().X();
    _fevt.tk_v[1]=top_tk.dir().Y();
    _fevt.tk_v[2]=top_tk.dir().Z();
    _fevt.tk_pchi2=top_tk.pchi2();
    _fevt.tk_plans=top_tk.plans();
  ROOT::Math::XYZPoint p=top_tk.extrapolate(82);
  _pex.SetXYZ(33-p.Y(),88-p.X(),p.Z()); //was 50-p.X() (100-12 connecteur strip)

  _fevt.pex_x[0]=_pex.X();
  _fevt.pex_x[1]=_pex.Y();
  _fevt.pex_x[2]=_pex.Z();
  hxy->Fill(_pex.X(),_pex.Y());
  hcount->Fill(9.);
  // printf("NCH %ld \n",e->tdcChannels().size());
  // getchar();
  uint32_t ninti=0;
  std::vector<lydaq::TdcChannel> vChannel; vChannel.clear();
  for (auto x:e->tdcChannels())
    {
      uint32_t coarsedif=ibc*80;
     
      while (coarsedif>((1<<24)-1))
	coarsedif -=((1<<24)-1);
      double ddt=(coarsedif-(x.tdcTime()/2.5))*2.5/200.;
      //ddt=ddt*2.5/200.;
       bool timeselected =
	 ((_geoRoot["general"]["noise"].asUInt()==0)&&(ddt<-3 && ddt>-8))||
	 ((_geoRoot["general"]["noise"].asUInt()==1)&&(ddt<-13 && ddt>-113));
	 
      // Noise if (ddt<-53 && ddt>-153)
	 if (timeselected)
	{
	  vChannel.push_back(x);
	  //printf("=======> %d %d %d %d \n",ibc,int(ddt),x.feb(),x.channel());
	  hch->Fill((x.feb()-14.)+x.channel());
	  _fevt.f_feb[ninti]=x.feb();
	  _fevt.f_channel[ninti]=x.channel();
	  _fevt.f_coarse[ninti]=x.coarse();
	  _fevt.f_fine[ninti]=x.fine();

	  ninti++;
	  _fevt.ninti=ninti;
	}
      hdt->Fill(ddt);

    }

  hinti->Fill(ninti*1.);
  if (ninti>0) {
    printf("BCID %d found in NCH %ld channels \n",ibc,e->tdcChannels().size());
    hxyt->Fill(_pex.X(),_pex.Y());

      
  }
  _run=e->run();

  bool cfound=this->stripStudy(vChannel,"FEB");
  if (cfound)
    {
      _fevt.found_feb=_selfeb;
      hxyf->Fill(_pex.X(),_pex.Y());
      if (_selfeb==14)
	hxyf14->Fill(_pex.X(),_pex.Y());
      if (_selfeb==15)
	hxyf15->Fill(_pex.X(),_pex.Y());
    }
  if (tEvents_!=NULL)
    {
      treeFile_->cd();
                
      tEvents_->Fill();
    }

  //getchar();
  if (_geoRoot["general"]["display"].asUInt()==0) return;
  if (TCHits==NULL)
    {
      TCHits=new TCanvas("TCHits","tChits1",900,900);
      TCHits->Modified();
      TCHits->Draw();
      TCHits->Divide(1,2);
    }
  TCHits->cd(1);
  hzx->SetMarkerStyle(25);
  hzx->SetMarkerColor(kRed);
  hzx->Draw("P");
  if (top_tk.size()>1) top_tk.linex()->Draw("SAME");
  TCHits->Modified();
  TCHits->Draw();
  TCHits->Update();
  TCHits->cd(2);
  hzy->SetMarkerStyle(22);
  hzy->SetMarkerColor(kGreen);
  hzy->Draw("P");
  if (top_tk.size()>1) top_tk.liney()->Draw("SAME");
  TCHits->Modified();
  TCHits->Draw();
  TCHits->Update();
  usleep(10000);
  TCHits->Update();

  getchar();
}
void binaryreader::processEvent(rbEvent* e)
{
  if (_event==0)
    {
      _geo->fillFebs(e->run());
      _geo->fillAlign(e->run());
      
      for (int i=0;i<255;i++)
	{
	  if (_geo->feb(i).isEmpty()) continue;
	  printf("FEB %d %d %d \n",i,_geo->feb(i).id,_geo->feb(i).tdc2strip[15]);
	  

	}
      

    }

  _event=e->gtc();
  
  uint8_t u[16],v[16],w[16];
  if (!_started) return;
  //printf("BR => %d %d %d \n",e->run(),e->event(),e->gtc());


  

  std::map<uint32_t,std::vector<uint32_t> > tm;tm.clear();
  uint32_t maxt=0;
  std::stringstream sraw;
  sraw<<"/pmr/";
  TH1* hfc=_rh->GetTH1(sraw.str()+"FrameCount");
  TH1* hft=_rh->GetTH1(sraw.str()+"FrameTime");
  TH1* hftm=_rh->GetTH1(sraw.str()+"MaxTime");
  TH1* hfts=_rh->GetTH1(sraw.str()+"FrameTimeSelected");
  TH1* hcount=_rh->GetTH1(sraw.str()+"Count");



  if (hfc==NULL)
    {
      hcount=_rh->BookTH1(sraw.str()+"Count",10,0.1,10.1);
      hfc=_rh->BookTH1(sraw.str()+"FrameCount",255,0.,255.);
      hft=_rh->BookTH1(sraw.str()+"FrameTime",65536 ,0.,2.);
      hftm=_rh->BookTH1(sraw.str()+"MaxTime",65536,0.,2.);
      
      hfts=_rh->BookTH1(sraw.str()+"FrameTimeSelected",10000,0.,10000.);
    }

  hcount->Fill(1.);
  bool f3pl=false;
  for (auto it=e->tCount().begin();it!=e->tCount().end();it++)
    {
      //std::cout<<"BC "<<it->first<<" Count "<<it->second<<std::endl;
    if ( it->second.count()>2)
      {
	//std::cout<<"BC "<<it->first<<" Count "<<it->second<<std::endl;
      f3pl=true;
      //  getchar();
      hcount->Fill(it->second.count()*1.);
      auto itm=e->tFrame().find(it->first);
      if (itm!=e->tFrame().end())
	{
	std::cout<<"BCM "<<itm->first<<" Count "<<itm->second.size()<<std::endl;
      	processCoincidence(e,it->first);
	}
      }
    }

  if (true) return;

  
  bool ramf=false;
  for (int id=0;id<MAXDIF;id++)
    {
      ramf=ramf||(e->frameCount(id)>126);
      if (e->frameCount(id))
      printf("ID %d => Count %d \n",id,e->frameCount(id));
      for (int j=0;j<e->frameCount(id);j++)
	{
	  uint32_t idx=e->iPtr(id,j);
	
	  if (e->bcid(idx)>maxt) maxt=e->bcid(idx);
	}
    }
  bool coinc=true;
  
  bool trigger=false;
  hftm->Fill(maxt*2E-7);
  printf(" Max time %d \n", maxt);
  if (maxt<100) return;
  //if (rf) return;
  bool rf=false;
  for (int id=0;id<MAXDIF;id++)
    if (e->frameCount(id))
      {
	uint16_t plane=(id>>4)&0xF;
	hfc->Fill(id*1.,e->frameCount(id));
	printf("PMR %x %d frames \n",id,e->frameCount(id));
	std::stringstream sraw1;
	sraw1<<"/pmr/ASIC"<<std::hex<<id<<std::dec<<"/";

	TH1* hp1=_rh->GetTH1(sraw1.str()+"Pad1");
	TH1* hftg=_rh->GetTH1(sraw1.str()+"FrameTime");
	TH1* hfc=_rh->GetTH1(sraw1.str()+"FrameCount");
	if (hp1==NULL)
	  {
	    hp1=_rh->BookTH1(sraw1.str()+"Pad1",64,0.,64.);
	    hftg=_rh->BookTH1(sraw1.str()+"FrameTime",65536,0.,2.);
	    hfc=_rh->BookTH1(sraw1.str()+"FrameCount",65,0.,65.);
	  
	  }

	for (int j=0;j<e->frameCount(id);j++)
	  {
	    uint32_t idx=e->iPtr(id,j);
	    int32_t dd=maxt-e->bcid(idx);

	    if (e->bcid(idx)<50) continue;
	    //if (e->frameCount(id)==127 && (dd<10)) continue;
	    std::bitset<64> bs,bs0,bs1;bs.reset();bs1.reset();bs0.reset();
	    for (int k=0;k<64;k++)
	      {
		if (e->pad0(idx,k)) bs0.set(k);
		if (e->pad1(idx,k)) bs1.set(k);
		  if (e->pad0(idx,k)||e->pad1(idx,k)) {bs.set(k);
		    //if (dd<100)
		    // hp1->Fill(k*1.);
		  }
		}
	    hfc->Fill(bs.count()*1.);
	    if (bs.count()>60) continue;
	    for (int k=0;k<64;k++)
	      {
		if (e->pad0(idx,k)) bs0.set(k);
		if (e->pad1(idx,k)) bs1.set(k);
		  if (e->pad0(idx,k)||e->pad1(idx,k)) {bs.set(k);
		    //if (dd<100)
		     hp1->Fill(k*1.);}
		}
	    hft->Fill((maxt-e->bcid(idx))*2E-7);
	    hftg->Fill((maxt-e->bcid(idx))*2E-7);
	    bool found=false;
	    std::map<uint32_t,std::vector<uint32_t> >::iterator itm=tm.find(e->bcid(idx));
	    std::map<uint32_t,std::vector<uint32_t> >::iterator itmm=tm.find(e->bcid(idx)-1);
	    std::map<uint32_t,std::vector<uint32_t> >::iterator itmp=tm.find(e->bcid(idx)+1);
	    found =(itm!=tm.end())||(itmp!=tm.end())||(itmm!=tm.end());
	    if (!found)
	      {
		std::vector<uint32_t> v;
		v.push_back(idx);
		std::pair<uint32_t,std::vector<uint32_t> > p(e->bcid(idx),v);
		tm.insert(p);

	      }
	    else
	      {
		if (itm!=tm.end()) itm->second.push_back(idx);
		if (itmm!=tm.end()) itmm->second.push_back(idx);
		if (itmp!=tm.end()) itmp->second.push_back(idx);
	      }

	    rf=rf|(itm->second.size()>=3);
	    /*	    
	    printf("%x frame bcid %.9d %8.5f %s %d\n",id,e->bcid(idx),e->bcid(idx)*2E-7,bs.to_string().c_str(),bs.count());
	    printf("%x frame bcid %.9d %8.5f %s %d\n",id,e->bcid(idx),e->bcid(idx)*2E-7,bs0.to_string().c_str(),bs0.count());
	    printf("%x frame bcid %.9d %8.5f %s %d\n ----------\n",id,e->bcid(idx),e->bcid(idx)*2E-7,bs1.to_string().c_str(),bs1.count());
	    */	    
	    if (j>5)
	      {
		//printf("....\n");
		//break;
	      }
	  }
      }

  if (ramf && trigger) return;
  if (!ramf && coinc) return;
  //std::cout<<"Selected-->"<<tm.size()<<std::endl;
  //  getchar();
  if (rf || true)
    {
      std::stringstream sres;
      sres.clear();
      std::bitset<16> planes;planes.reset();
      int32_t nplanesmin=2;
      if (trigger)  nplanesmin=0;
      for (auto x:tm)
	{
	  if (x.second.size()<3) continue;
	  int32_t dd=maxt-x.first;
	  if (dd>20 and trigger) continue;
	  if (dd<20 && coinc)  continue;
	  sres.str("");
	  sres<<"\t"<<x.first<<" : "<<x.second.size()<<" => "<<maxt<<" "<<dd<<std::hex;
	  planes.reset();
	  memset(u,0,16);
	  memset(v,0,16);
	  memset(w,0,16);
	  for (auto y:x.second)
	    {
	      uint32_t ig=y/MAXFRAME/FSIZE;
	      uint32_t plane=(ig>>4)&0xF;
	      u[plane] = u[plane] || ( (ig&0xF)==2 || (ig&0xF)==4);
	      v[plane] = v[plane] || ( (ig&0xF)==3 || (ig&0xF)==5);
	      w[plane] = w[plane] || ( (ig&0xF)==6 || (ig&0xF)==7);
	      planes.set(plane,1);
	      sres<<" "<<ig<<"("<<plane<<")";
	    }
	  uint32_t nplanes=planes.count();
	  sres<<std::dec<<" Plans "<<x.second.size()<<"->"<<nplanes<<": "<<planes<<std::endl;
	  //if (nplanes==1) continue;

	  
	  if (nplanes>nplanesmin)
	    {
	  if (nplanes>nplanesmin)
	    std::cout<<x.first<<"=>SELECTEDS "<<planes.count()<<" "<<dd<<" "<<x.first<<" " <<maxt<<"\t RESP "<<sres.str()<<std::endl;
	  //getchar();
	  for (auto y:x.second)
	    {
	      uint32_t ig=y/MAXFRAME/FSIZE;
	      uint16_t plane=(ig>>4)&0xF;
	      bool muv=u[plane]&&v[plane];
	      bool muw=u[plane]&&w[plane];
	      bool mvw=v[plane]&&w[plane];
	      if (!(muv||muw||mvw)) continue;
	      uint8_t cpos=(ig&0xF);
	      bool udir= ( (ig&0xF)==2 || (ig&0xF)==4);
	      bool vdir= ( (ig&0xF)==3 || (ig&0xF)==5);
	      bool wdir= ( (ig&0xF)==6 || (ig&0xF)==7);
	      int sshift=0;
	      if (cpos==4 || cpos==5 || cpos==7) sshift=64;
	      //if (!measure) continue;	      
	      std::stringstream sraw1;
	      std::stringstream splane;

	      sraw1<<"/pmr/ASIC"<<std::hex<<ig<<std::dec<<"/";
	      splane<<"/pmr/PLANE"<<plane<<"/";
	      
	      TH1* hp1=_rh->GetTH1(sraw1.str()+"Pad1Sel");
	      TH1* hsu=_rh->GetTH1(splane.str()+"Ustrip");
	      TH1* hsv=_rh->GetTH1(splane.str()+"Vstrip");
	      TH1* hsw=_rh->GetTH1(splane.str()+"Wstrip");
	      TH1* hfcs=_rh->GetTH1(sraw1.str()+"FrameCountSel");
	      TH1* hftse=_rh->GetTH1(sraw1.str()+"FrameTimeSel");

	      if (hp1==NULL)
		{
		  hp1=_rh->BookTH1(sraw1.str()+"Pad1Sel",64,0.,64.);
		  hfcs=_rh->BookTH1(sraw1.str()+"FrameCountSel",65,0.,65.);
		  hftse=_rh->BookTH1(sraw1.str()+"FrameTimeSel",500,0.,500.);
		  hsu=_rh->BookTH1(splane.str()+"Ustrip",128,0.,128.);
		  hsv=_rh->BookTH1(splane.str()+"Vstrip",128,0.,128.);
		  hsw=_rh->BookTH1(splane.str()+"Wstrip",128,0.,128.);

		}
	      hftse->Fill(maxt-e->bcid(y)*1.);
	      uint32_t idx=y;
	      std::bitset<64> bs;bs.reset();
	      for (int k=0;k<64;k++)
		if (e->pad0(idx,k)||e->pad1(idx,k)) {
		  bs.set(k);hp1->Fill(k*1.);
		  if (udir) hsu->Fill(k+sshift*1.);
		  if (vdir) hsv->Fill(k+sshift*1.);
		  if (wdir) hsw->Fill(k+sshift*1.);
		}
	      hfcs->Fill(bs.count()*1.);
	      //fprintf(stderr,"%x frame bcid %.9d %8.5f %s %d\n",ig,e->bcid(idx),e->bcid(idx)*2E-7,bs.to_string().c_str(),bs.count());
	      // getchar();
	    }
    
	}
	}
    }
  
  return;
  // Selected
  bool sel=false;
  for (int id=0;id<MAXDIF;id++)
    if (e->frameCount(id))
      {
	std::stringstream sraw1;
	sraw1<<"/pmr/ASIC"<<std::hex<<id<<std::dec<<"/";
	
	TH1* hp1=_rh->GetTH1(sraw1.str()+"Pad1Sel");
	if (hp1==NULL)
	  {
	    hp1=_rh->BookTH1(sraw1.str()+"Pad1Sel",64,0.,64.);
	  }
	for (int j=0;j<e->frameCount(id);j++)
	  {
	    uint32_t idx=e->iPtr(id,j);


	    auto itm=tm.find(e->bcid(idx));
	    if (itm==tm.end())
	      {
		std::cout<<"OOOPS"<<std::endl;
		continue;
	      }
	    if (itm->second.size()<3) continue;
	    int32_t dd=maxt-itm->first;
	    //if (dd>20) continue;
	    std::bitset<64> bs;bs.reset();
	    for (int k=0;k<64;k++)
	      if (e->pad0(idx,k)||e->pad1(idx,k)) {bs.set(k);hp1->Fill(k*1.);}
	    if (bs.count()>50) continue;
	    sel=true;
	    //if (maxt-itm->first>50) continue;
	    //if (maxt-itm->first>5000) continue;
	    hfts->Fill((maxt-e->bcid(idx))*1.);
	 
	    //if (bs.count()>20) continue;
	    //printf("%x frame bcid %d %f %s %d\n",id,e->bcid(idx),e->bcid(idx)*2E-7,bs.to_string().c_str(),bs.count());
	    //Xgetchar();
	  }
      }
  //if (sel) getchar();
  
  
  //getchar();
}
int32_t binaryreader::TPrincipalComponents(double result[21],float zmin,float zmax)
{
	double resultc[21];
	uint32_t nh=0;
	double xb=0,yb=0,zb=0;
	double wt=0.;

	double fp=DBL_MAX;
	double lp=-DBL_MAX;
	double fx=DBL_MAX;
	double lx=-DBL_MAX;
	double fy=DBL_MAX;
	double ly=-DBL_MAX;
	TPrincipal tp(3,"D");
	double xp[3];
	memset(result,0,21*sizeof(double));
	//INFO_PRINT("%d vector size\n",v.size());
	for (auto it=_vPoints.begin();it!=_vPoints.end();it++)
	{
	  ROOT::Math::XYZPoint& iht=(*it);
	  if (iht.Z()<zmin) continue;
	  if (iht.Z()>zmax) continue;
		//INFO_PRINT("%x %d %d \n",iht,iht.I(),iht.J());
		//INFO_PRINT("%f %f \n",iht.x(),iht.y());
		//INFO_PRINT("%f %f \n",iht.X(),iht.Y());
		double w=1.;
		xb+=iht.X()*w;
		yb+=iht.Y()*w;
		zb+=iht.Z()*w;
		wt+=w;
		nh++;
		xp[0]=iht.X();
		xp[1]=iht.Y();
		xp[2]=iht.Z();
		//printf("XP  %f %f %f \n",xp[0],xp[1],xp[2]);
		tp.AddRow(xp);
	}

	if (nh<2) return -1;
	tp.MakePrincipals();
	// store barycenter
	const TVectorD* fvb=tp.GetMeanValues();
	// printf("barycentre %f %f %f \n \t %f %f %f \n", xb/wt,yb/wt,zb/wt,(*fvb)[0],(*fvb)[1],(*fvb)[2]);
	result[0]=(*fvb)[0];
	result[1]=(*fvb)[1];
	result[2]=(*fvb)[2];

	const TVectorD* fva=tp.GetEigenValues();
	result[3]=(*fva)[2];
	result[4]=(*fva)[1];
	result[5]=(*fva)[0];

	//fva->Print();
	// printf("eigen results %g %g %g \n",(*fva)[0],(*fva)[1],(*fva)[2]);

	//tp.Print("MSEV");

	// store principal axis
	const TMatrixD* fvv=tp.GetEigenVectors();
	//Matrix<double,3,3> vv=eigensolver.eigenvectors();
	result[6]=(*fvv)(0,2);
	result[7]=(*fvv)(1,2);
	result[8]=(*fvv)(2,2);
	// printf("eigen vector results %g %g %g \n",result[7],result[7],result[8]);
	result[9]=(*fvv)(0,1);

	result[10]=(*fvv)(1,1);
	result[11]=(*fvv)(2,1);
	// printf("eigen vector results %g %g %g \n",result[9],result[10],result[11]);
	result[12]=(*fvv)(0,0);
	result[13]=(*fvv)(1,0);
	result[14]=(*fvv)(2,0);
	// printf("eigen vector results %g %g %g \n",result[12],result[13],result[14]);
	// Store First and last Z
	result[15]=fp;
	result[16]=lp;
	result[17]=fx;
	result[18]=lx;
	result[19]=fy;
	result[20]=ly;
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


bool binaryreader::stripStudy(std::vector<lydaq::TdcChannel>& vChannel,std::string subdir)
{
  float ch1_dt[128];

  bool clusterFound=false;
  memset(ch1_dt,0,128*sizeof(float));

  uint32_t triggerChannel=1;
  float dtmin=-540,dtmax=-585;
  bool noisy=false,_display=false;
  uint32_t chamber=1;
  //dtmin+=100;dtmax+=100;
  _strips.clear();
  std::vector<TdcChannel*> c_strip[128];
  for (int i=0;i<128;i++) c_strip[i].clear();
  float maxtime=0,mttime=1E12;uint32_t nch=0,ntrg=0;
  std::bitset<48> stfeb(0);
   for (auto x=vChannel.begin();x!=vChannel.end();x++)
    {
      if (x->channel()==1) continue;
      if (x->pedSubTime(_geo->feb(x->feb()))< mttime)
	mttime= x->pedSubTime(_geo->feb(x->feb()));
    }
  for (auto x=vChannel.begin();x!=vChannel.end();x++)
    {
      // Drop Trigger Channel
      if (x->channel()==1) continue;
      if (x->pedSubTime(_geo->feb(x->feb()))>(40+mttime)) continue;
      std::stringstream sraw;
      sraw<<"/"<<subdir<<"/Chamber"<<chamber<<"/Raw/";
      TH1* hchan=_rh->GetTH1(sraw.str()+"Channels");
      TH1* hstrips=_rh->GetTH1(sraw.str()+"Strips");
      if (hchan==NULL)
	{
	  hchan=_rh->BookTH1(sraw.str()+"Channels",48*16,0.,48.*16);
	  hstrips=_rh->BookTH1(sraw.str()+"Strips",96,0,96);
	}
      hchan->Fill(x->feb()*48+x->channel());
      hstrips->Fill( x->side(_geo->feb(x->feb()))*48+x->detectorStrip(_geo->feb(x->feb())));
	    
	  // printf("%f %f \n",dtmin,dtmax);
	  // getchar();
	  // Book and fill time to trigger
	  //dtm[x->feb()][ 
      c_strip[x->detectorStrip(_geo->feb(x->feb()))].push_back(&(*x));

      //jsonFebInfo f=_geo->feb(x->feb());
      //for (int ic=0;ic<49;ic++) std::cout<<f.tdc2strip[ic]<<" ";
      printf("FEB %d Channel %d Strip %d Shift %d Side %d %f raw %f\n",x->feb(),x->channel(),x->detectorStrip(_geo->feb(x->feb())),_geo->feb(x->feb()).stripShift,x->side(_geo->feb(x->feb())),x->pedSubTime(_geo->feb(x->feb())),x->pedSubTime(_geo->feb(x->feb())));
	      //if (x->feb()==11)
      



	  
      nch++;
    }
  //if (nch) getchar();
  std::stringstream srcc;
  srcc<<"/"<<subdir<<"/Chamber"<<chamber<<"/";
      
  TH1* hfr=_rh->GetTH1(srcc.str()+"FebCount");
  TH1* hfrs=_rh->GetTH1(srcc.str()+"FebCountSel");
  TH2* hcor=_rh->GetTH2(srcc.str()+"CCOR");

      
  if (hfr==NULL)
    {
      
      hfr=_rh->BookTH1(srcc.str()+"FebCount",49,0.,49.);
      hfrs=_rh->BookTH1(srcc.str()+"FebCountSel",49,0.,49.);
      hcor=_rh->BookTH2(srcc.str()+"CCOR",100,0,100,100,0,100);

    }
  hfr->Fill(48.);
  hfrs->Fill(48.);
  for (int i=0;i<48;i++)
    if (stfeb[i]>0) hfr->Fill(i*1.);




  for (int i=0;i<128;i++)
    {
      if (c_strip[i].size()==0) continue;
      for (auto x:c_strip[i])
	{

	  if (x->side(_geo->feb(x->feb()))==0) continue;
	  for (int j=0;j<128;j++)
	    {
	      if (c_strip[j].size()==0) continue;
	      for (auto y:c_strip[j])
		{
		  if (y->side(_geo->feb(y->feb()))==1) continue;
		  if (x->channel()!=y->channel())
		    hcor->Fill((x->feb()-14)*48+x->channel()+0.2,(y->feb()-14)*48+y->channel()+0.2);
		}
	    }
	}

    }
      maxtime=maxtime*1E-9;

      //fprintf(stderr," Maxtime %d %f %d %f \n",chamber,maxtime,nch,nch/maxtime/6500);
      //getchar();
      bool dostop=false;int nstrip=0;
      uint16_t febc[48];
      memset(febc,0,48*sizeof(uint16_t));
      std::bitset<64> stb(0);
      std::bitset<64> stb0(0);
      std::bitset<64> stb1(0);
#ifdef ALGOSTANDARD
      for (int i=0;i<128;i++)
	{
	  if (c_strip[i].size()>0)
	    {
	      fprintf(stderr,"Chamber %d Strip %d # %d \n",chamber,i,c_strip[i].size());
	      nstrip++;
	      if (i<48)
		stb.set(i,1);
	    }

	  for (auto x:c_strip[i])
	    {
	      if (x->side(_geo->feb(x->feb()))==0 && i<48)
		stb0.set(i,1);
	      else
		if (i<48)
		  stb1.set(i,1);
	    }
	
	  if (c_strip[i].size()>2)
	    {
	      dostop=true;
	      fprintf(stderr," too many hits \n");
	      for (auto x:c_strip[i])
		{
		  //  fprintf(stderr,"\t %d %d %d %f %d \n",x->feb(),x->channel(), x->side(_geo->feb(x->feb())),x->tdcTime(),x->detectorStrip(_geo->feb(x->feb())));
		}
	      //getchar();
	    }
	  if (c_strip[i].size()==2)
	    {
	      double t0=-1,t1=-1;
	      for (auto x:c_strip[i])
		{

		  //fprintf(stderr,"\t %d %d %f %f \n",x->channel(), x->side(_geo->feb(x->feb())),x->tdcTime(),x->tdcTime());
		  double dt=_geo->feb(x->feb()).dtc[x->channel()];


		  if (t0<0 &&  x->side(_geo->feb(x->feb()))==0)
		    {
		      t0=x->pedSubTime(_geo->feb(x->feb()))-dt;

		      printf("T0 %d %d %d %d %f %f dt=%f \n",x->feb(),x->channel(),x->coarse(),x->fine(),x->tdcTime(),t0,dt);
		    }
		  if (t1<0 &&  x->side(_geo->feb(x->feb()))==1)
		    {
		      t1=x->pedSubTime(_geo->feb(x->feb()))-dt;
		      printf("T1 %d %d %d %d %f %f \n",x->feb(),x->channel(),x->coarse(),x->fine(),x->tdcTime(),t1);
		    }
		  if(t0>0 && t1>0 )
		    {
		      //if (abs(t1-t0)<5) continue;
		      febc[x->feb()]++;
		      std::cout<<x->feb()<<" FEBC "<< febc[x->feb()]<<std::endl;
		      if (_geo->feb(x->feb()).polarity==-1)
			{
			  double tt=t1;
			  t1=t0;
			  t0=tt;
			}
		    
		      //Lmana::TdcStrip ts(_geo->feb(x->feb()).chamber,x->feb(),x->detectorStrip(_geo->feb(x->feb())),t0,t1,_geo->feb(x->feb()).timePedestal[x->detectorStrip( _geo->feb(x->feb()))]);
		      if (chamber==1)
			{
			  Lmana::TdcStrip ts(_geo->feb(x->feb()).chamber,x->feb(),x->detectorStrip(_geo->feb(x->feb())),t0,t1,_geo->feb(x->feb()).timePedestal[x->detectorStrip( _geo->feb(x->feb()))]);
					     //ch1_dt[x->detectorStrip( _geo->feb(x->feb()))+1]);
			    _strips.push_back(ts);
			}
		      // else
		      // 	{
		      // 	  Lmana::TdcStrip ts(_geo->feb(x->feb()).chamber,x->feb(),x->detectorStrip(_geo->feb(x->feb())),t0,t1,ch2_dt[x->detectorStrip( _geo->feb(x->feb()))+1]);
		      // 	  _strips.push_back(ts);
		      // 	}

		    }
		}



	    
	    }
	}
#else
      for (int i=0;i<128;i++)
	{
	  if (c_strip[i].size()>0)
	    {
	      fprintf(stderr,"Chamber %d Strip %d # %ld \n",chamber,i,c_strip[i].size());
	      nstrip++;
	      if (i<48)
		stb.set(i,1);
	    }
	  else
	    continue;
	  
	  for (int ic=0;ic<c_strip[i].size();ic++)
	    {
	      lydaq::TdcChannel* x =c_strip[i][ic];
	      if (x->side(_geo->feb(x->feb()))!=1) continue;
	      double t1=x->pedSubTime(_geo->feb(x->feb()));
	      //fprintf(stderr,"Side 1 FEB %d Channel %d Strip %d Shift %d Side %d %f raw %f\n",x->feb(),x->channel(),x->detectorStrip(_geo->feb(x->feb())),_geo->feb(x->feb()).stripShift,x->side(_geo->feb(x->feb())),x->pedSubTime(_geo->feb(x->feb())),t0);
	      for (int jc=0;jc<c_strip[i].size();jc++)
		{
		  lydaq::TdcChannel* y =c_strip[i][jc];
		  if (y->side(_geo->feb(y->feb()))!=0) continue;
		  if (y->used()) continue;
		  double t0=y->pedSubTime(_geo->feb(y->feb()));
		  //  fprintf(stderr,"Side 0 FEB %d Channel %d Strip %d Shift %d Side %d %f raw %f\n",y->feb(),y->channel(),y->detectorStrip(_geo->feb(y->feb())),_geo->feb(y->feb()).stripShift,y->side(_geo->feb(y->feb())),y->pedSubTime(_geo->feb(y->feb())),t1);
		  //fprintf(stderr,"Strip  t0 %f  \n",t0);
		  //fprintf(stderr,"Strip  t1 %f  \n",t1); 

		  Lmana::TdcStrip ts(_geo->feb(x->feb()).chamber,x->feb(),x->detectorStrip(_geo->feb(x->feb())),t0,t1,_geo->feb(x->feb()).timePedestal[x->detectorStrip( _geo->feb(x->feb()))]);
		  //if (i==1008)
		    fprintf(stderr,"Strip %d c1 %d c0 %d  t1 %f t0 %f ypos %f expected %f \n",
			    i,x->channel(),y->channel(),t1,t0,ts.ypos(),_pex.Y()); 


		  std::stringstream srcs;
		  srcs<<"/Align/strip"<<i;
		  //std:cout<<srcs.str()<<std::endl;
		  TH1* hdts=_rh->GetTH1(srcs.str());
		  if (hdts==NULL)
		    hdts=_rh->BookTH1(srcs.str(),300,-300,300.);
		  std::stringstream srcd;
		  srcd<<"/Align/dt"<<i;
		  //std:cout<<srcs.str()<<std::endl;
		  TH1* hddt=_rh->GetTH1(srcd.str());
		  if (hddt==NULL)
		    hddt=_rh->BookTH1(srcd.str(),300,-10.,40.);

		  if (abs(_pex.X()-ts.xpos())<3.)
		    {
		      hdts->Fill(_pex.Y()-ts.ypos());
		      hddt->Fill(t1-t0);
		    }
		
		  if (ts.ypos()<-100 || ts.ypos()>350.) continue;
		  x->setUsed(true);
		  y->setUsed(true);
		   _strips.push_back(ts);
		}
	    }
	  
	  if (i==1008)
	    getchar();


	}
	
#endif

      for (int i=0;i<48;i++)
	if (febc[i]>0) hfrs->Fill(i*1.);

      if (dostop) return clusterFound;
      if (stb.count()>48) return clusterFound;
      noisy=(stb.count()>48);
      //for (int i=0;i<48;i++)
      // if (febc[i]>=10) return true;
      //std::cout<<stb<<std::endl;
      /////////////////////////////////////////////////////////// STRIP STUDIES ////////////////////////////////////////////////////////////////
      std::stringstream s_src;
      s_src<<"/"<<subdir<<"/Chamber"<<chamber<<"/Strips/";
		  

      TH2* hs_xy=_rh->GetTH2(s_src.str()+"XY");

      TH2* hs_xym=_rh->GetTH2(s_src.str()+"XYMore");
      TH2* hs_xym20=_rh->GetTH2(s_src.str()+"XYMore20");
      TH2* hs_xym50=_rh->GetTH2(s_src.str()+"XYMore50");
      TH2* hs_xysel=_rh->GetTH2(s_src.str()+"XYSel");
      TH2* hs_xymin=_rh->GetTH2(s_src.str()+"XYMin");
      TH2* hs_dxdy=_rh->GetTH2(s_src.str()+"DXDY");
      TH2* hs_dxdy20=_rh->GetTH2(s_src.str()+"DXDY20");
      TH2* hs_dxdy50=_rh->GetTH2(s_src.str()+"DXDY50");
      TH2* hs_dxvsx=_rh->GetTH2(s_src.str()+"dxvsx");
      TH2* hs_dyvsx=_rh->GetTH2(s_src.str()+"dyvsx");
      TH2* hs_dxvsy=_rh->GetTH2(s_src.str()+"dxvsy");
      TH2* hs_dyvsy=_rh->GetTH2(s_src.str()+"dyvsy");
      TH1* hs_abst=_rh->GetTH1(s_src.str()+"TOA");
      TH1* hs_dtmore=_rh->GetTH1(s_src.str()+"AbsoluteTimeDelta");
      TH1* hs_dxmore=_rh->GetTH1(s_src.str()+"Xtosel");
      TH1* hs_dxsel=_rh->GetTH1(s_src.str()+"Xtoex");
      TH1* hs_dthrem8=_rh->GetTH1(s_src.str()+"ADTHREM888");
      TH1* hs_dtlrem8=_rh->GetTH1(s_src.str()+"ADTLREM888");
      TH1* hs_dthrfr4=_rh->GetTH1(s_src.str()+"ADTHRFR4");
      TH1* hs_dtlrfr4=_rh->GetTH1(s_src.str()+"ADTLRFR4");
      TH1* hs_t1t0150=_rh->GetTH1(s_src.str()+"T1T0_150");
      TH1* hs_dist=_rh->GetTH1(s_src.str()+"Dist2Ext");

      if (hs_xy==NULL)
	{
	      
	  hs_xy=_rh->BookTH2(s_src.str()+"XY",48,0.,48.,420,-60.,250.);
	  hs_xym=_rh->BookTH2(s_src.str()+"XYMore",48,0.,48.,420,-60.,250.);
	  hs_xym20=_rh->BookTH2(s_src.str()+"XYMore20",48,0.,48.,420,-60.,250.);
	  hs_xym50=_rh->BookTH2(s_src.str()+"XYMore50",48,0.,48.,420,-60.,250.);
	  hs_xysel=_rh->BookTH2(s_src.str()+"XYSel",48,0.,48.,420,-60.,250.);
	  hs_xymin=_rh->BookTH2(s_src.str()+"XYMin",48,0.,48.,420,-60.,250.);
	  hs_dxdy=_rh->BookTH2(s_src.str()+"DXDY",100,-10.,10.,600,-20.,220.);
	  hs_dxdy20=_rh->BookTH2(s_src.str()+"DXDY20",100,-10.,10.,600,-20.,220.);
	  hs_dxdy50=_rh->BookTH2(s_src.str()+"DXDY50",100,-10.,10.,600,-20.,220.);
	  
	  hs_dxvsx=_rh->BookTH2(s_src.str()+"dxvsx",48,0.,48.,100,-10.,10.);
	  hs_dyvsx=_rh->BookTH2(s_src.str()+"dyvsx",48,0.,48.,100,-20.,20.);
	  hs_dxvsy=_rh->BookTH2(s_src.str()+"dxvsy",105,-60.,250.,100,-10.,10.);
	  hs_dyvsy=_rh->BookTH2(s_src.str()+"dyvsy",105,-60.,250.,100,-20.,20.);
	  hs_abst=_rh->BookTH1(s_src.str()+"TOA",100,-570.,-600.);
	  hs_dtmore=_rh->BookTH1(s_src.str()+"AbsoluteTimeDelta",20000,-2000.,2000.);
	  hs_dxmore=_rh->BookTH1(s_src.str()+"Xtosel",200,-30.,30.);
	  hs_dxsel=_rh->BookTH1(s_src.str()+"Xtoex",200,-30.,30.);
	  hs_dthrem8=_rh->BookTH1(s_src.str()+"ADTHREM888",20000,-2000.,2000.);
	  hs_dtlrem8=_rh->BookTH1(s_src.str()+"ADTLREM888",20000,-2000.,2000.);
	  hs_dthrfr4=_rh->BookTH1(s_src.str()+"ADTHRFR4",20000,-2000.,2000.);
	  hs_dtlrfr4=_rh->BookTH1(s_src.str()+"ADTLRFR4",20000,-2000.,2000.);
	  hs_t1t0150=_rh->BookTH1(s_src.str()+"T1T0_150",2000,-5.,100.);

	  hs_dist=_rh->BookTH1(s_src.str()+"Dist2Ext",2000,0.,200.);
	}

      double s_distmin=99999999.;
      double s_tmsel=0,s_xmsel=0,s_ymsel=0;
      double t0_sel=0,t1_sel=0;

      for (auto x:_strips)
	{
	  hs_abst->Fill(x.TM());
	  hs_xy->Fill(x.X(),x.Y());
	  double dist=sqrt((_pex.X()-x.X())*(_pex.X()-x.X())+
			   (_pex.Y()-x.Y())*(_pex.Y()-x.Y()));
	  hs_dist->Fill(dist);
	  if (dist<s_distmin)
	    {
	      s_distmin=dist;
	      s_tmsel=x.TM();
	      t0_sel=x.t0();
	      t1_sel=x.t1();
	      s_xmsel=x.X();
	      s_ymsel=x.Y();
	    }
	}
      
      for (auto x:_strips)
	{
	  if (x.TM()==s_tmsel) 
	    {
	       hs_xymin->Fill(x.X(),x.Y());
	       hs_dxdy->Fill(x.X()-_pex.X(),x.Y()-_pex.Y());
	       if (x.Y()>125 &&x.Y()<160)
		 hs_dxsel->Fill(x.X()-_pex.X());
	      if (s_distmin<8)
		{
	      //hposcma->Fill(x.X(),x.Y());
		  hs_xysel->Fill(x.X(),x.Y());

		  hs_dxvsx->Fill(_pex.X(),x.X()-_pex.X());
		  hs_dyvsx->Fill(_pex.X(),x.Y()-_pex.Y());
		  hs_dxvsy->Fill(_pex.Y(),x.X()-_pex.X());
		  hs_dyvsy->Fill(_pex.Y(),x.Y()-_pex.Y());
		}
	     
	    }
	  else
	    {
	      hs_xym->Fill(x.X(),x.Y());
	      

	       double dist=sqrt((_pex.X()-x.X())*(_pex.X()-x.X())+
				(_pex.Y()-x.Y())*(_pex.Y()-x.Y()));
	       if (dist>15. && x.Y()>125 &&x.Y()<160)
		 {
		   printf("End of strip issue %d \n",_event);
		   // getchar();
		   hs_dtmore->Fill(x.TM()-s_tmsel);
		   hs_dxmore->Fill(x.X()-s_xmsel);
		   hs_t1t0150->Fill(x.t1()-x.t0());
		 }

		 if (dist>15. && x.Y()>170)
		 {
		   printf("Xtalk issue %d \n",_event);
		   // getchar();
		 
		 }


	       if (x.strip()<=16 && dist>15. && x.Y()>125 &&x.Y()<160)
		{
		  hs_dthrem8->Fill(x.t0()-t0_sel);
		  hs_dtlrem8->Fill(x.t1()-t1_sel);
		}
	      if (x.strip()>=17 && dist>15. && x.Y()>125 &&x.Y()<160)
		{
		  hs_dthrfr4->Fill(x.t0()-t0_sel);
		  hs_dtlrfr4->Fill(x.t1()-t1_sel);
		}
		
	      if (abs(x.TM()-s_tmsel)<20)
		{
		  hs_xym20->Fill(x.X(),x.Y());
		  hs_dxdy20->Fill(x.X()-_pex.X(),x.Y()-_pex.Y());
		}
	      else
		{
		  hs_xym50->Fill(x.X(),x.Y());
		  hs_dxdy50->Fill(x.X()-_pex.X(),x.Y()-_pex.Y());
		}

	    }
	}
      
      /////////////////////////////////////////////////////////////////////////// CLUSTER STUDIES //////////////////////////////////////////////////////////////////////////////
      std::vector<Lmana::TdcCluster> vclus;
      vclus.clear();
      float step=2.;
      //if (chamber==1) step=4.;
      for (auto it=_strips.begin();it!=_strips.end();it++)
	{
	  //fprintf(stderr,"%d %d %f %f %f %f \n",it->chamber(),it->strip(),it->xpos(),it->ypos(),it->t0(),it->t1());
	  if (it->chamber()!=chamber) continue;
	  //	      if (it->ypos()<-10 || it->ypos()>-0.2) continue;
	  bool found=false;
	  for (auto ic=vclus.begin();ic!=vclus.end();ic++)
	    {
	      if (ic->isAdjacent((*it),step))
		{
		  ic->addStrip((*it));
		  found=true;
		  break;
		}
	    }
	  if(!found)
	    {
	      Lmana::TdcCluster c;
	      c.addStrip((*it));
	      vclus.push_back(c);
	    }
	}

      // Merge adjacent cluster
      bool merged=false;
      // printf("vclus size %ld \n",vclus.size());

      for (auto it=vclus.begin();it!=vclus.end();it++)
	{
	  for (auto jt=it+1;jt!=vclus.end();)
	    {
	      bool adj=false;
	      for (int i=0;i<jt->size();i++)
		{
		  if (it->isAdjacent(jt->strip(i),step))
		    {adj=true;break;}
		}
	      if (adj)
		{
		  merged=true;
		  printf("Merigng cluster \n");
		  for (int i=0;i<jt->size();i++)
		    {
		      it->addStrip(jt->strip(i));
		    }
		  vclus.erase(jt);
		}
	      else
		++jt;
	    }
	}
      printf("vclus size after %ld \n",vclus.size());
      //getchar();
      //if (merged)
      //	getchar();
      std::stringstream src;
      src<<"/"<<subdir<<"/Chamber"<<chamber<<"/ClusterNew/";
		  
      TH2* hposs=_rh->GetTH2(src.str()+"XYStrip");
      TH2* hposc=_rh->GetTH2(src.str()+"XY");
      TH2* hposc1=_rh->GetTH2(src.str()+"XY1");
      TH2* hposcm=_rh->GetTH2(src.str()+"XYMore");
      TH2* hposcm20=_rh->GetTH2(src.str()+"XYMore20");
      TH2* hposcm50=_rh->GetTH2(src.str()+"XYMore50");
      TH2* hposcma=_rh->GetTH2(src.str()+"XYMax");
      TH2* hposok=_rh->GetTH2(src.str()+"XYSel");
      TH2* hposdma=_rh->GetTH2(src.str()+"DIST");
      TH2* hposdma20=_rh->GetTH2(src.str()+"DIST20");
      TH2* hposdma50=_rh->GetTH2(src.str()+"DIST50");
      TH2* hpdxvsx=_rh->GetTH2(src.str()+"dxvsx");
      TH2* hpdyvsx=_rh->GetTH2(src.str()+"dyvsx");
      TH2* hpdxvsy=_rh->GetTH2(src.str()+"dxvsy");
      TH2* hpdyvsy=_rh->GetTH2(src.str()+"dyvsy");
      TH2* hposx=_rh->GetTH2(src.str()+"XYX");
      TH1* hncl=_rh->GetTH1(src.str()+"Clusters");
      TH1* hmulc=_rh->GetTH1(src.str()+"ClusterSize");
      TH1* hmulc1=_rh->GetTH1(src.str()+"ClusterSize1");
      TH1* hns=_rh->GetTH1(src.str()+"nstrip");
      TH1* hns0=_rh->GetTH1(src.str()+"nstrip0");
      TH1* hns1=_rh->GetTH1(src.str()+"nstrip1");
      TH1* hns2=_rh->GetTH1(src.str()+"nstrip2");
      TH1* htoa=_rh->GetTH1(src.str()+"TOA");
      TH1* hdtr0=_rh->GetTH1(src.str()+"DT0");
      TH1* hdtr1=_rh->GetTH1(src.str()+"DT1");
      TH1* hdtmore=_rh->GetTH1(src.str()+"AbsoluteTimeDelta");

      if (hposc==NULL)
	{
	      
	  hposc=_rh->BookTH2(src.str()+"XY",48,0.,48.,420,-60.,250.);
	  hposs=_rh->BookTH2(src.str()+"XYStrip",48,0.,48.,420,-60.,250.);
	  hposc1=_rh->BookTH2(src.str()+"XY1",48,0.,48.,420,-60.,250.);
	  hposcm=_rh->BookTH2(src.str()+"XYMore",48,0.,48.,420,-60.,250.);
	  hposcm20=_rh->BookTH2(src.str()+"XYMore20",48,0.,48.,420,-60.,250.);
	  hposcm50=_rh->BookTH2(src.str()+"XYMore50",48,0.,48.,420,-60.,250.);
	  hposcma=_rh->BookTH2(src.str()+"XYMax",48,0.,48.,420,-60.,250.);
	  hposok=_rh->BookTH2(src.str()+"XYSel",48,0.,48.,180,-20.,170.);
	  hposdma=_rh->BookTH2(src.str()+"DIST",100,-10.,10.,600,-20.,220.);
	  hposdma20=_rh->BookTH2(src.str()+"DIST20",100,-10.,10.,600,-20.,220.);
	  hposdma50=_rh->BookTH2(src.str()+"DIST50",100,-10.,10.,600,-20.,220.);
	  hpdxvsx=_rh->BookTH2(src.str()+"dxvsx",48,0.,48.,100,-10.,10.);
	  hpdyvsx=_rh->BookTH2(src.str()+"dyvsx",48,0.,48.,100,-20.,20.);
	  hpdxvsy=_rh->BookTH2(src.str()+"dxvsy",105,-60.,250.,100,-10.,10.);
	  hpdyvsy=_rh->BookTH2(src.str()+"dyvsy",105,-60.,250.,100,-20.,20.);
	  hposx=_rh->BookTH2(src.str()+"XYX",48,0.,48.,420,-60.,250.);
	  hncl=_rh->BookTH1(src.str()+"Clusters",32,0.,32.);
	  hmulc=_rh->BookTH1(src.str()+"ClusterSize",32,0.,32.);
	  hmulc1=_rh->BookTH1(src.str()+"ClusterSize1",32,0.,32.);
	  hns=_rh->BookTH1(src.str()+"nstrip",48,0.,48.);
	  hns0=_rh->BookTH1(src.str()+"nstrip0",48,0.,48.);
	  hns1=_rh->BookTH1(src.str()+"nstrip1",48,0.,48.);
	  hns2=_rh->BookTH1(src.str()+"nstrip2",48,0.,48.);
	  htoa=_rh->BookTH1(src.str()+"TOA",100,-570.,-600.);
	  hdtr0=_rh->BookTH1(src.str()+"DT0",400,-20.,20.);
	  hdtr1=_rh->BookTH1(src.str()+"DT1",20000,-2000.,2000.);
	  hdtmore=_rh->BookTH1(src.str()+"AbsoluteTimeDelta",20000,-2000.,2000.);
	}
      hns->Fill(nstrip);
      hns0->Fill(stb0.count());
      hns1->Fill(stb1.count());
      hns2->Fill(_strips.size());
      printf(" ===> %d  %ld strips , Number of clusters %ld \n",chamber,_strips.size(),vclus.size());
      //if (vclus.size()>0)
      hncl->Fill(vclus.size()*1.);
      uint32_t maxs=0;
      double distmin=99999999.;
      double tmsel=0;
      for (auto x:vclus)
	{
	  // if (vclus.size()>1)
	  //   {
	  if (_display)
	    {
	      fprintf(stderr,"\t %f %f %d \n",x.X(),x.Y(),x.size());
	      for (int i=0;i<x.size();i++)
		{
		  fprintf(stderr,"\t \t  %5.1f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f \n",
		  x.strip(i).xpos(),x.strip(i).ypos(),x.strip(i).shift(),x.strip(i).t0(),x.strip(i).t1(),(x.strip(i).t0()+x.strip(i).t1())/2.,0.0);
		}
	    }
	  //   }
	  if (vclus.size()==1)
	    for (int i=0;i<x.size();i++)
	      {
		htoa->Fill((x.strip(i).t0()+x.strip(i).t1())/2.);
		//hdtr0->Fill(x.strip(i).t0());
		//hdtr1->Fill(x.strip(i).t1());
	      }

	  if (x.size()>30) continue;
	  hposc->Fill(x.X(),x.Y());
	  if (vclus.size()==1)
	    {hposc1->Fill(x.X(),x.Y());

	    }
	  
	  
	  hmulc->Fill(x.size()*1.);
	  //printf("%d %f %f \n",x.size(),x.X(),x.Y());
	  double dist=sqrt((_pex.X()-x.X())*(_pex.X()-x.X())+
			   (_pex.Y()-x.Y())*(_pex.Y()-x.Y()));
	  if (dist<distmin)
	    {
	      distmin=dist;
	      tmsel=x.TM();
	    }
      
	  if (x.size()>maxs) maxs=x.size();
	}
      //if (vclus.size()>1) getchar();
      int nc=0;
      // float dy[48]={0.,
      // 		    63.6,-18.1,-13.1,27.1,-39.5,22.2,4.9,39.5,12.3,8.6,29.7,22.4,
      // 		    -16.1,-30.3,0,-42.8,65.4,4.3,-1.5,56.0,-48.1,30.3,39.0,-20.9,
      // 		    -12.0,36.1,25.8,16.8,5.1,13.9,32.1,-31.11,
      // 		    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      float dy[33]={0.,
      		    62.8,-6.9,-10.9,49.3,-28.8,26.5,14.1,-2.3,
		    5.9,29.6,-3.2,13.7,23.1,10.3,0,0,
		    43.3,4.7,-0.9,56.5,-47.9,54.2,40.0,-20.1,
		    11.6,36.4,25.9,16.7,5.1,14.0,31.9,-30.6};
      memset(dy,0,33*sizeof(float));
      //clusterFound=vclus.size()>0;
      for (auto x:vclus)
	{
	  nc++;
	  if (x.size()>16) continue;
	  float L=160.;
	  float v=160./8.7;
	  float xl=(L-x.Y()*v)/2.+L/2.;
	  

	  hposx->Fill(x.X(),x.Y());
	  fprintf(stderr,"Xext %f X %f Yex %f Y %f \n",_pex.X(),x.X(),_pex.Y(),x.Y());
	  x.Print();
	  if (abs(x.X()-_pex.X())>6) continue;
	  //getchar();

	  if (x.TM()==tmsel)
	    {
	      //hposcma->Fill(x.X(),x.Y());
	      hposcma->Fill(x.X(),x.Y()+dy[int(x.X())]);
	      hposdma->Fill(x.X()-_pex.X(),x.Y()-_pex.Y());
	      hpdxvsx->Fill(_pex.X(),x.X()-_pex.X());
	      hpdyvsx->Fill(_pex.X(),x.Y()-_pex.Y());
	      hpdxvsy->Fill(_pex.Y(),x.X()-_pex.X());
	      hpdyvsy->Fill(_pex.Y(),x.Y()-_pex.Y());
	      
	      clusterFound|=abs(x.X()-_pex.X())<6&& abs(x.Y()+dy[int(x.X())]-_pex.Y())<9.5 ;
	      if (abs(x.X()-_pex.X())<6&& abs(x.Y()+dy[int(x.X())]-_pex.Y())<9.5)
		{
		  hposok->Fill(x.X(),x.Y());
		  hmulc1->Fill(x.size()*1.);
		}
	      hdtr0->Fill(x.X()-_pex.X());
	      _selfeb=x.strip(0).dif();
	      for (int i=0;i<x.size();i++)
		{
		  hposs->Fill(x.strip(i).xpos(),x.strip(i).ypos());
		  std::stringstream srcs;
		  srcs<<"/Align/Sstrip"<<int(x.strip(i).xpos());
		  //std:cout<<srcs.str()<<std::endl;
		  TH1* hdts=_rh->GetTH1(srcs.str());
		  if (hdts==NULL)
		    hdts=_rh->BookTH1(srcs.str(),200,-20.,20.);

		  hdts->Fill(_pex.Y()-x.strip(i).ypos()-dy[int(x.X())]);
		}
	      break;
	    }
	}
      int ncp=0;
      for (auto x:vclus)
	{
	  ncp++;
	  if (x.size()>16) continue;

	  //if (ncp==nc) continue;
	  if (x.TM()==tmsel) continue;
	  hposcm->Fill(x.X(),x.Y());
	  hdtmore->Fill(x.TM()-tmsel);
	  if (abs(x.TM()-tmsel)<20)
	    {
	      hposcm20->Fill(x.X(),x.Y());
	      hposdma20->Fill(x.X()-_pex.X(),x.Y()-_pex.Y());
	    }
	  else
	    {
	      hposcm50->Fill(x.X(),x.Y());
	      hposdma50->Fill(x.X()-_pex.X(),x.Y()-_pex.Y());
	    }
	}
	//if (nstrip>0&& !clusterFound)
	//	getchar();
      // if (_display)
      // 	{
      // 	  this->drawHits(chamber);
      // 	  if (chamber==2 ) getchar();
      // 	}


  return clusterFound;
}

void binaryreader::createTrees(std::string s)
{

  treeFile_ = new TFile(s.c_str(), "recreate");
  treeFile_->cd();

  tEvents_ = new TTree("events", "Events");
  tEvents_->SetAutoSave(50000000);
  
  tEvents_->Branch("bcid", &_fevt.bc, "bcid/l");
  tEvents_->Branch("run", &_fevt.run, "run/i");
  tEvents_->Branch("event", &_fevt.event, "event/i");
  tEvents_->Branch("gtc", &_fevt.event, "gtc/i");
  tEvents_->Branch("npad", &_fevt.npad, "npad/s");
  tEvents_->Branch("pad_dif", &_fevt.pad_dif, "pad_dif[npad]/s");
  tEvents_->Branch("pad_asic", &_fevt.pad_asic, "pad_asic[npad]/b");
  tEvents_->Branch("pad_channel", &_fevt.pad_channel, "pad_channel[npad]/b");
  tEvents_->Branch("ntel", &_fevt.ntel, "ntel/s");
  tEvents_->Branch("tel_x", &_fevt.tel_x, "tel_x[ntel]/F");
  tEvents_->Branch("tel_y", &_fevt.tel_y, "tel_y[ntel]/F");
  tEvents_->Branch("tel_z", &_fevt.tel_z, "tel_z[ntel]/F");
  tEvents_->Branch("tk_x", &_fevt.tk_x, "tk_x[3]/F");
  tEvents_->Branch("tk_v", &_fevt.tk_x, "tk_v[3]/F");
  tEvents_->Branch("tk_pchi2", &_fevt.tk_pchi2, "tk_pchi2/F");
  tEvents_->Branch("tk_plans", &_fevt.tk_plans, "tk_plans/l");
  tEvents_->Branch("pex_x", &_fevt.pex_x, "pex_x[3]/F");

  tEvents_->Branch("found_feb", &_fevt.found_feb, "found_feb/b");
  tEvents_->Branch("ninti", &_fevt.ninti, "ninti/s");
  tEvents_->Branch("f_feb", &_fevt.f_feb, "f_feb[ninti]/b");
  tEvents_->Branch("f_channel", &_fevt.f_channel, "f_channel[ninti]/b");
  tEvents_->Branch("f_fine", &_fevt.f_fine, "f_fine[ninti]/b");
  tEvents_->Branch("f_coarse", &_fevt.f_coarse, "f_coarse[ninti]/l");



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


extern "C" 
{
  // loadDHCALAnalyzer function creates new LowPassDHCALAnalyzer object and returns it.  
  rbProcessor* loadProcessor(void)
  {
    return (new binaryreader);
  }
  // The deleteDHCALAnalyzer function deletes the LowPassDHCALAnalyzer that is passed 
  // to it.  This isn't a very safe function, since there's no 
  // way to ensure that the object provided is indeed a LowPassDHCALAnalyzer.
  void deleteProcessor(rbProcessor* obj)
  {
    delete obj;
  }
}
