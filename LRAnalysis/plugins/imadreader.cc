#include "imadreader.hh"
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
using namespace zdaq;
imadreader::imadreader() : _run(0),_started(false),_fdOut(-1),_totalSize(0),_event(0) {}
void imadreader::init(uint32_t run)
{
  _run=run; 
  _event=0;
  _started=true;
}
void imadreader::loadParameters(Json::Value params)
{
  /*
  if (params.isMember("directory"))
    _directory=params["directory"].asString();
   if (params.isMember("dummy"))
     _dummy=(params["dummy"].asUInt()!=0);
  */
}
void imadreader::end(uint32_t run)
{
  _started=false;

}

void imadreader::processRunHeader(std::vector<uint32_t> header)
{
}
void imadreader::processEvent(rbEvent* e)
{

  if (!_started) return;
  printf("BR => %d %d %d \n",e->run(),e->event(),e->gtc());
  for (int id=0;id<MAXDIF;id++)
    if (e->frameCount(id))
      {
	printf("GRIC %x %d frames \n",id,e->frameCount(id));
	for (int j=0;j<e->frameCount(id);j++)
	  {
	    uint32_t idx=e->iPtr(id,j);
	    printf("frame bcid %d %f \n",e->bcid(idx),e->bcid(idx)*2E-7);
	    for (int k=0;k<64;k++)
	      if (e->pad1(idx,k)) printf("pad %d touches\n",k);
	  }
      }
}
extern "C" 
{
    // loadDHCALAnalyzer function creates new LowPassDHCALAnalyzer object and returns it.  
  rbProcessor* loadProcessor(void)
    {
      return (new imadreader);
    }
    // The deleteDHCALAnalyzer function deletes the LowPassDHCALAnalyzer that is passed 
    // to it.  This isn't a very safe function, since there's no 
    // way to ensure that the object provided is indeed a LowPassDHCALAnalyzer.
  void deleteProcessor(rbProcessor* obj)
    {
      delete obj;
    }
}
