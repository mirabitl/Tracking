#include <dlfcn.h>
#include "tdcrb.hh"
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
#include <map>
#include <bitset>
#include "TCanvas.h"
#include "TdcMapping.hh"
#include <fstream>
#include <dirent.h>
#include <fnmatch.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "shmwriterProcessor.hh"
#include <arpa/inet.h>
#include "LirocAccess.hh"







using namespace zdaq;
tdcrb::tdcrb(std::string dire) : _directory(dire),_run(0),_started(false),_fdIn(-1),_totalSize(0),_event(0),_geo(NULL),_t0(2E50),_t(0),_tspill(0)
			       ,_readoutTotalTime(0),_runType(0),_dacSet(0),_fdOut(-1),_bxId0(0),_nrmax(5000000),_nrfirst(0)
{_rh=DCHistogramHandler::instance();

}

void  tdcrb::registerProcessor(std::string name)
{
  std::stringstream s;
  s<<"lib"<<name<<".so";
  void* library = dlopen(s.str().c_str(), RTLD_NOW);


  //LOG4CXX_INFO(_logZdaq," Error "<<dlerror()<<" Library open address "<<std::hex<<library<<std::dec);
  // Get the loadFilter function, for loading objects
  rbProcessor* (*create)();
  create = (rbProcessor* (*)())dlsym(library, "loadProcessor");
  //LOG4CXX_INFO(_logZdaq," Error "<<dlerror()<<" file "<<s.str()<<" loads to processor address "<<std::hex<<create<<std::dec);
  //printf("%s %x \n",dlerror(),(unsigned int) create);
  // printf("%s lods to %x \n",s.str().c_str(),(unsigned int) create); 
  //void (*destroy)(Filter*);
  // destroy = (void (*)(Filter*))dlsym(library, "deleteFilter");
  // Get a new filter object
  rbProcessor* a=(rbProcessor*) create();
  _processors.push_back(a);
}

void tdcrb::findDataSet(std::string dirp,uint32_t runask)
{
  std::stringstream spat;
  //int runask=atol(argv[1]);
  spat<<"SMM*"<<runask<<"*.dat";
  //spat<<"SMM*"<<argv[1]<<"*.dat";
  struct dirent **namelist;
  int n;
  std::cout<<"Pattern "<<spat.str()<<std::endl;
  //std::string dirp="/data/srv02/RAID6/Dome0718";

  //dirp=".";
  n = scandir(dirp.c_str(), &namelist, NULL, alphasort);
  if (n < 0)
    perror("scandir");
  else {
    while (n--) {

      if (fnmatch(spat.str().c_str(), namelist[n]->d_name, 0)==0)
	{
	  printf("%s %d \n", namelist[n]->d_name,fnmatch(spat.str().c_str(), namelist[n]->d_name, 0));
	  printf("found\n");
	  std::stringstream sf;
	  sf<<dirp<<"/"<< namelist[n]->d_name;
	  this->addRun(runask,sf.str());
	}
      free(namelist[n]);
    }
    free(namelist);
  }

}
void tdcrb::geometry(std::string name)
{
  _geo=new jsonGeo(name);
  //_analyzer->setGeometry(_geo);
  for (auto x:_processors)
    {
      x->loadParameters(_geo->root());
    }
}

void tdcrb::pull(std::string name,zdaq::buffer* buf,std::string sourcedir)
{
  std::stringstream sc,sd;
  sc.str(std::string());
  sd.str(std::string());
  sc<<sourcedir<<"/closed/"<<name;
  sd<<sourcedir<<"/"<<name;
  int fd=::open(sd.str().c_str(),O_RDONLY);
  if (fd<0) 
    {
      printf("%s  Cannot open file %s : return code %d \n",__PRETTY_FUNCTION__,sd.str().c_str(),fd);
      //LOG4CXX_FATAL(_logShm," Cannot open shm file "<<fname);
      return ;
    }
  int size_buf=::read(fd,buf->ptr(),0x20000);
  buf->setPayloadSize(size_buf-(3*sizeof(uint32_t)+sizeof(uint64_t)));
  //printf("%d bytes read %x %d \n",size_buf,cbuf[0],cbuf[1]);
  ::close(fd);
  ::unlink(sc.str().c_str());
  ::unlink(sd.str().c_str());
}

uint32_t tdcrb::numberOfDataSource() {return 1;}



void tdcrb::processRawEvent(uint64_t idx)
{
}
void tdcrb::clearShm()
{
  std::vector<std::string> vnames;
  shmwriterProcessor::ls("/dev/shm/monitor",vnames);
  std::stringstream sc,sd;
  sc.str(std::string());
  sd.str(std::string());
  for (auto name:vnames)
    {
      sc<<"/dev/shm/monitor/closed/"<<name;
      sd<<"/dev/shm/monitor/"<<name;
      ::unlink(sc.str().c_str());
      ::unlink(sd.str().c_str());

    }
}
void tdcrb::monitor()
{
  _started=true;
  _geo->fillFebs(_run);
  _geo->fillAlign(_run);

  while (_started)
    {
      std::vector<std::string> vnames;
      shmwriterProcessor::ls("/dev/shm/monitor",vnames);
      for (auto x:vnames)
	{std::cout<<x<<std::endl;
	  //EUDAQ_WARN("Find file "+x);
	  //continue;
	  zdaq::buffer* b=new zdaq::buffer(0x80000);
	  this->pull(x,b,"/dev/shm/monitor");
	  if (b->detectorId()==255)
	    {
	      uint32_t* buf=(uint32_t*) b->payload();
	      printf("NEW RUN %d \n",b->dataSourceId());
	      _run=b->dataSourceId();


	      for (int i=0;i<b->payloadSize()/4;i++)
		{
		  printf("%d ",buf[i]);
		}

	      _runType=buf[0];
	      if (_runType==1)
		_dacSet=buf[1];
	      if (_runType==2)
		_vthSet=buf[1];
	      printf("\n Run type %d DAC set %d VTH set %d \n",_runType,_dacSet,_vthSet);
	      //getchar();
	      //	      _analyzer->jEvent()["runtype"]=_runType;
	      continue;
	    }
	  uint64_t idx_storage=b->eventId(); // usually abcid
	  std::map<uint64_t,std::vector<zdaq::buffer*> >::iterator it_gtc=_eventMap.find(idx_storage);
	  if (it_gtc!=_eventMap.end())
	    it_gtc->second.push_back(b);
	  else
	    {
	      std::vector<zdaq::buffer*> v;
	      v.clear();
	      v.push_back(b);
          
	      std::pair<uint64_t,std::vector<zdaq::buffer*> > p(idx_storage,v);
	      _eventMap.insert(p);
	      it_gtc=_eventMap.find(idx_storage);
	    }
	  if (it_gtc->second.size()==this->numberOfDataSource())
	    {
	      if (it_gtc->first%100==0)
		printf("GTC %lu %lu  %d\n",it_gtc->first,it_gtc->second.size(),this->numberOfDataSource());
	      this->processRawEvent(idx_storage);
	    }
	}
      usleep(50000);	
  
    }
  
}

void tdcrb::open(std::string filename)
{
  if (_geo==NULL)
    {
      std::cout<<"Please speicfy a geometry"<<std::endl;
      exit(0);
    }
  _fdIn= ::open(filename.c_str(), O_RDONLY | O_NONBLOCK,S_IRWXU);
  if (_fdIn<0)
    {
      perror("Ici No way to store to file :");
      //std::cout<<" No way to store to file"<<std::endl;
      return;
    }  
  _event=0;
  _started=true;
}
void tdcrb::close()
{
  _started=false;
  ::sleep(1);
  if (_fdIn>0)
    {
      ::close(_fdIn);
      _fdIn=-1;
    }


}
uint32_t tdcrb::totalSize(){return _totalSize;}
uint32_t tdcrb::eventNumber(){return _event;}
uint32_t tdcrb::runNumber(){return _run;}
void tdcrb::Read()
{
  //  if (_analyzer == NULL) _analyzer= new lmana::TdcAnalyzer(_rh);
  for (auto p:_processors)
    p->init();
  int nfile=0;
  _nread=0;
  for (std::vector<std::pair<uint32_t,std::string> >::iterator it=_files.begin();it!=_files.end();it++)
    {
      std::cout<<"NEW File "<<it->first<<" "<<it->second<<std::endl;
      _run=it->first;
      std::stringstream sff;
      sff<<"sudo chmod o+r "<<it->second;
      system(sff.str().c_str());
      this->open(it->second);
      if (_fdOut<=0)
	this->read();
      //else
      //this->streamout(4);
      
      this->close();
      // std::stringstream sroot;
      // sroot<<"/tmp/histo"<<_run<<"_"<<nfile<<".root";
      //  _rh->writeHistograms(sroot.str());

      //  std::cout<<"Writing histos "<<sroot.str()<<std::endl;
      //getchar();
      nfile++;
    }
}

void tdcrb::read()
{
  
  zdaq::buffer b(0x100000);
  int last=-1;
  uint64_t _eventChannel[4096*8];
  //  std::vector<lydaq::TdcChannel> _vAll;
  uint32_t _eventChannels;
  _geo->fillFebs(_run);
  _geo->fillAlign(_run);
  uint32_t nid=0;
  while (_started)
    {
      if (!_started) return;
      uint32_t theNumberOfDIF=0;
      // To be implemented
      if (_fdIn>0)
	{
	  _idx=0;

	  
	  int ier=::read(_fdIn,&_event,sizeof(uint32_t));
	  _nread++;
	  if (last == _event) nid++;
	  if (ier<0 || nid>2)
	    {
	      printf("Cannot read Event anymore %d %d %d \n ",ier,last,_event);return;
	    }
	  //if (last == _event) continue;
	   if (_nread>(_nrmax+_nrfirst)) return;
	  //else
	  last=_event;
	  if (_event%100==0)
	    INFO_PRINTF("Event read %d \n",_event);
      
	  ier=::read(_fdIn,&theNumberOfDIF,sizeof(uint32_t));
	  if (ier<0)
	    {
	      printf("Cannot read anymore number of DIF %d \n ",ier);return;
	    }
	  else
	    if (_event%100==0)
	      INFO_PRINTF("================> Event %d Number of DIF found %d \n",_event,theNumberOfDIF);
	  INFO_PRINTF("================> Event %d Number of DIF found %d \n",_event,theNumberOfDIF);
	  uint32_t difFound[256];
	  memset(difFound,0,256*sizeof(uint32_t));
	  uint32_t trigFound[256];
	  memset(trigFound,0,256*sizeof(uint32_t));
	  //	  _analyzer->clear();
	  
	  memset(_eventChannel,0,4096*8*sizeof(uint64_t));
	  _eventChannels=0;
	  bool _initialised=false;
	  for (uint32_t idif=0;idif<theNumberOfDIF;idif++) 
	    {
	      uint32_t tbcid=0;
	      //DEBUG_PRINTF("\t writing %d bytes",idata[SHM_BUFFER_SIZE]);
	      //(*iv)->compress();
	      uint32_t bsize=0;
	      // _totalSize+=bsize;
	      ier=::read(_fdIn,&bsize,sizeof(uint32_t));
	      if (ier<0)
		{
		  printf("Cannot read anymore  DIF Size %d \n ",ier);return;
		}
	      else
		if (_event%100==0)
		  DEBUG_PRINTF("\t DIF size %d \n",bsize);
	  
	      ier=::read(_fdIn,b.ptr(),bsize);
	      if (ier<0)
		{
		  printf("Cannot read anymore Read data %d \n ",ier);return;
		}
	      b.setPayloadSize(bsize-(3*sizeof(uint32_t)+sizeof(uint64_t)));
	      b.uncompress();
	      memcpy(&_buf[_idx], b.payload(),b.payloadSize());

	      b.setDetectorId(b.detectorId()&0xFF);
	      INFO_PRINTF("\t \t %d %x %d %ld %d %d %d\n",b.detectorId()&0XFF,b.dataSourceId(),b.eventId(),b.bxId(),b.payloadSize(),bsize,_idx);
	      
	      _bxId=b.bxId();
	      if (_bxId0==0) _bxId0=_bxId;
	      uint32_t _detId=b.detectorId()&0xFF;
	      //DEBUG_PRINTF("DETID %d \n",_detId);
	      // getchar();
	      if (_detId==255)
		{
		  uint32_t* buf=(uint32_t*) b.payload();
		  printf("NEW RUN %d \n",_event);
		  _run=_event;


		  for (int i=0;i<b.payloadSize()/4;i++)
		    {
		      printf("%d ",buf[i]);
		    }
		  _difId=b.dataSourceId();
		  _runType=buf[0];
		  if (_runType==1)
		    _dacSet=buf[1];
		  if (_runType==2)
		    _vthSet=buf[1];
		  printf("\n Run type %d DAC set %d VTH set %d \n",_runType,_dacSet,_vthSet);
		  // getchar();

		}
	      if (_detId==170)
		{
		  _theEvent.setDifType(true);
		  //uint32_t* ibuf=(uint32_t*) b.payload();

		  //uint64_t* lbuf=(uint64_t*) b.payload();
		  uint8_t* bb=&_buf[_idx];
		  uint16_t* sbuf=(uint16_t*)&bb[1];
		  uint32_t* ibuf=(uint32_t*) &bb[1];
		  uint64_t* lbuf=(uint64_t*) &bb[1];
		  //		    itemp[0]=_event;
		  //
		  //  itemp[1]=_lastGTC;
		  //  ltemp[1]=_lastABCID;
		  //  itemp[4]= _event;
		  //  itemp[5]=_adr;
		  //  itemp[6]=length;
		  // printf("%d %ld %d %x %d \n",ibuf[1],lbuf[1],ibuf[4],ibuf[5],ibuf[6]);
		  /**
		     int cnt=0,cntmax=20;
		     for (int i=0;i<b.payloadSize();i++)
		     {
		     printf("%.2x ",bb[i]);
		     cnt++;
		     if (cnt==cntmax)
		     {
		     uint8_t ub1=0,ub2=0,ub3=0;
		     for (int ii=0;ii<8;ii++)
		     {
		     if (bb[i-19] & (1<<ii))
		     ub1|=(1<<(7-ii));
		     }
			  
		     printf("=> %d %d\n",bb[i-19],ub1);
		     cnt=0;
		     if (cntmax==20) cntmax=20;
		     }
		     }
		  
		     getchar();
		  */
		  uint32_t idstart=liroc::Unpacker::getStartOfLIROC(bb,b.payloadSize(),0);
		  fprintf(stdout,"\n Size %d IDSTART  %d %d  %d %llu => %d %d Dif %d \n",b.payloadSize(),idstart,liroc::Unpacker::getID(bb,idstart),liroc::Unpacker::getGTC(bb,idstart),liroc::Unpacker::getAbsoluteBCID(bb,idstart),liroc::Unpacker::getLastTriggerBCID(bb,idstart),liroc::Unpacker::getBCID(bb,idstart),liroc::Unpacker::getLastTriggerBCID(bb,idstart)-liroc::Unpacker::getBCID(bb,idstart));
		  //getchar();
		  if (!_initialised)
		    {
		      _theEvent.init(_run,_event,liroc::Unpacker::getGTC(bb,idstart),liroc::Unpacker::getAbsoluteBCID(bb,idstart));
		      _initialised=true;
		    }
		  
		  _theEvent.setCalibrationInfos(_runType,_vthSet);
		  _theEvent.setFrameCount(liroc::Unpacker::getID(bb,idstart),(b.payloadSize()-_idx-20)/DU_FRAME_SIZE);

		  fprintf(stderr,"%d vthset %d \n",_runType,_vthSet);
		  liroc::Ptr* d= new liroc::Ptr(&_buf[_idx],b.payloadSize());
		  _theEvent.difList().push_back(d);
		  _idx+=b.payloadSize();
		}
     
	    }
	  bool hadprint=false;
	  for (auto itd= _theEvent.difList().begin();itd!= _theEvent.difList().end();itd++)
	    {
	      _theEvent.setFrameCount((*itd)->getID(),(*itd)->getNumberOfFrames());
	      std::cout<<(*itd)->getID()<<"  frames "<<(*itd)->getNumberOfFrames()<<std::endl;
	      for (int j=0;j<(*itd)->getNumberOfFrames();j++)
		{
		  //uint64_t* bp=( uint64_t*)&(*itd)->getFramePtr(j)[4];
		  uint32_t* bc=( uint32_t*)&(*itd)->getFramePtr(j)[0];
		  uint32_t channel=(*itd)->getFrameChannel(j);
		  memcpy(_theEvent.frame(_theEvent.iPtr((*itd)->getID(),j)),(*itd)->getFramePtr(j),FSIZE);
		  bool doprint=(_theEvent.bcid(_theEvent.iPtr((*itd)->getID(),j))<100000) &&( channel==3 || channel ==2);
		  hadprint=hadprint||doprint;
		  //doprint=(channel==10 && _vthSet==554);
		  doprint=true;
		  //std::cout<<channel<<std::endl;
		  if (doprint){
		  for ( int ii=0;ii<4;ii++)
		    fprintf(stderr,"%.2x ",(*itd)->getFramePtr(j)[ii]);
		    fprintf(stderr,"\t ");
		    
		  std::cout<<j<<" "<<_theEvent.iPtr((*itd)->getID(),j)<<" "<<std::hex<<(*bc)<<std::dec<<std::endl;
		  }

		 if (doprint){
		std::cout<<" BCID "<<_theEvent.bcid(_theEvent.iPtr((*itd)->getID(),j));
		std::cout<<" Channel "<<(*itd)->getFrameChannel(j)<<std::endl;
			  
		 }
		}
		 if (hadprint)
		   getchar();

	    }
	  printf(" %d %d TCOUNT %d TMAP %d \n",_run,_event,_theEvent.frameCount(1),_theEvent.tFrame().size());

	  //if (oneselected)
	    for (auto p:_processors)
	      p->processEvent(&_theEvent);

	  //_analyzer->fullAnalysis(_vAll);

	}

    }
}
void tdcrb::end()
{
  for (auto p:_processors)
    p->end();

  //  _analyzer->end();
  std::stringstream sr;
  //sr<<_geo->general()["directory"].asString()<<"/histo"<<_run<<"_0.root";
  sr<<"/tmp/histo"<<_run<<"_0.root";
  _rh->writeHistograms(sr.str());


  
}
