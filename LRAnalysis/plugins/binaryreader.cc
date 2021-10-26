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
binaryreader::binaryreader() : _run(0),_started(true),_fdOut(-1),_totalSize(0),_event(0) {}
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
void binaryreader::processEvent(rbEvent* e)
{
  uint8_t u[16],v[16],w[16];
  //if (!_started) return;
  printf("BR => %d %d %d %d \n",e->run(),e->event(),e->gtc(),e->seuil());
  if (e->seuil()!=0)
    {this->scurveAnalysis(e);return;}

}

void binaryreader::scurveAnalysis(rbEvent *e)
{

  if (e->seuil()==0) return;
  std::cout<<"Event "<<_event<<" GTC"<<e->gtc()<<" Vth set "<<e->seuil()<<std::endl;
  //fflush(stdout);
  for (int id = 0; id < MAXDIF; id++)
    {
      //std::cout<<"Event "<<_event<<" GTC"<<e->gtc()<<" ID "<<id<<" " <<e->frameCount(id)<<std::endl;
      if (e->frameCount(id))
	{
	
	  std::stringstream sraw1;
	  sraw1 << "/LR/SCURVE" << std::hex << id << std::dec << "/";
	
	  TH1 *hp1 = _rh->GetTH1(sraw1.str() + "Padc1");
	  if (hp1 == NULL)
	    {
	      for (int i=0;i<64;i++)
		{
		  std::stringstream srpc("");
		  srpc<<sraw1.str()<<"Padc"<<i;
		  TH1* hpc = _rh->BookTH1(srpc.str(),1024,0.,1024);
		
		}
	    }
	  uint64_t lb0=0,lb1=0,lb2=0,lbcid=0;
	  for (int j = 0; j < e->frameCount(id); j++)
	    {
	      uint32_t idx = e->iPtr(id, j);
	    
	      if (e->bcid(idx) < 0)
		continue;
	      if ((e->bcid(idx)-lbcid)<3)
		continue;
	      lbcid=e->bcid(idx);
	      uint32_t k=e->channel(idx);
	      if (e->seuil()==556 && k==10)
		{
		  fprintf(stderr,"iptr %d bcid %d \n",idx,e->bcid(idx));
		  getchar();
		}

	      std::stringstream srpc("");
	      srpc<<sraw1.str()<<"Padc"<<k;
	    
	      TH1* hpc= _rh->GetTH1(srpc.str());
	      //std::cout<<srpc.str()<<std::endl;
	      hpc->Fill(e->seuil()*1.);
	      lb0 |=(1<<k);


	    }

	}

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
