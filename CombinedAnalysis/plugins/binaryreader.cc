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

//using namespace zdaq;
binaryreader::binaryreader() : _run(0),_started(false),_fdOut(-1),_totalSize(0),_event(0) {}
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
  _started=false;

}

void binaryreader::processRunHeader(std::vector<uint32_t> header)
{
}
static TCanvas* TCHits=NULL;
void binaryreader::processCoincidence(rbEvent* e,uint32_t ibc)
{
  TH2* hzx=_rh->GetTH2("/gric/ZX");
  TH2* hzy=_rh->GetTH2("/gric/ZY");
  TH2* hxy=_rh->GetTH2("/gric/XY");
  TH1* hcount=_rh->GetTH1("/gric/Count");
  if (hzx==NULL)
    {
      hcount=_rh->BookTH1("/gric/Count",10,0.1,10.1);
      hzx=_rh->BookTH2("/gric/ZX",400,0.,200.,100.,0.,100.);
      hzy=_rh->BookTH2("/gric/ZY",400,0.,200.,80.,0.,80.);
      hxy=_rh->BookTH2("/gric/XY",100,0.,100.,80.,0.,80.);
	    }
  hzx->Reset();
  hzy->Reset();
  hcount->Fill(1.);
  _vPoints.clear();
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
	      recoPoint pt(chid);
	      //pt.setPlan(chid);
	      _geo->convert(d->getID(),d->getFrameAsicHeader(ifra),ipad,&pt);
	      pt.SetZ(_geo->chamberInfo(chid).z0);
	      _vPoints.push_back(pt);
	      if (_geoRoot["general"]["display"].asUInt()!=0)
	      //printf("Chamber %d %f \n",chid,_geo->chamberInfo(chid).z0);
		printf("Point %d %d %d %f %f %f \n",d->getID(),d->getFrameAsicHeader(ifra),ipad,pt.X(),pt.Y(),pt.Z());
	      hzx->Fill(pt.Z(),pt.X());
	      hzy->Fill(pt.Z(),pt.Y());
	    }
	}
      hfcs->Fill(ph.count()*1.);
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
  printf("Shower ax %f ay %f bx %f by %f \n",ax,ay,bx,by);
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
  if (top_tk.orig().X()>47) return;
  printf(" TK ax %f ay %f bx %f by %f \n",
	   top_tk.dir().X(),
	   top_tk.dir().Y(),
	   top_tk.orig().X(),
	 top_tk.orig().Y());

  ROOT::Math::XYZPoint p=top_tk.extrapolate(85);
  hxy->Fill(p.X(),p.Y());
  hcount->Fill(9.);
  printf("NCH %d \n",e->tdcChannels().size());
  // getchar();
  for (auto x:e->tdcChannels())
    {
      if ((ibc-x.bcid())<200 &&(ibc-x.bcid())>-200 )
	printf("%d %d %d \n",ibc,x.bcid(),ibc-x.bcid());
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
  TCHits->Modified();
  TCHits->Draw();
  TCHits->Update();
  TCHits->cd(2);
  hzy->SetMarkerStyle(22);
  hzy->SetMarkerColor(kGreen);
  hzy->Draw("P");
  TCHits->Modified();
  TCHits->Draw();
  TCHits->Update();
  usleep(10000);
  TCHits->Update();

  getchar();
}
void binaryreader::processEvent(rbEvent* e)
{
  uint8_t u[16],v[16],w[16];
  if (!_started) return;
  printf("BR => %d %d %d \n",e->run(),e->event(),e->gtc());

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
      std::cout<<"BC "<<it->first<<" Count "<<it->second<<std::endl;
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
  printf(" Max time %f \n", maxt);
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
