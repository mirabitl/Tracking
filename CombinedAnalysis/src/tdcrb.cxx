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
#include "SdhcalPmrAccess.hh"







using namespace zdaq;
using namespace sdhcal;
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
  std::map<uint32_t,uint64_t> mbx;
  mbx.clear();
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
	  if (ier<0 || last==_event)
	    {
	      printf("Cannot read Event anymore %d %d %d \n ",ier,last,_event);return;
	    }
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
	  _theEvent.tdcChannels().clear();
	  bool _initialised=false;
	  float treal=0;
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

	      auto ib=mbx.find(b.dataSourceId());
	      if (ib==mbx.end())
		{
		  std::pair<uint32_t,uint64_t> p(b.dataSourceId(),0);
		  mbx.insert(p);
		}
	      
	      b.setDetectorId(b.detectorId()&0xFF);
	      INFO_PRINTF("\t \t %d %x %d %lx %d %d %d delta %ld %ld \n",b.detectorId()&0XFF,b.dataSourceId(),b.eventId(),b.bxId(),b.payloadSize(),bsize,_idx, mbx[b.dataSourceId()], b.bxId()-mbx[b.dataSourceId()]);
	      treal=(b.bxId()-mbx[b.dataSourceId()])*2E-7;
	      _bxId=b.bxId();
	      mbx[b.dataSourceId()]=_bxId;
	      //continue;
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
	      if (_detId==150)
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
		  
		     int cnt=0,cntmax=20;
		     for (int i=0;i<30;i++)
		     {
		     printf("%.2x ",bb[i]);
		     }
		     printf("\n");
		     /**
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
		  uint32_t idstart=sdhcal::PMRUnpacker::getStartOfPMR(bb,b.payloadSize(),0);
		  //fprintf(stdout,"\n IDSTART %d %d  %d %llu => %d %d Dif %d \n",idstart,sdhcal::PMRUnpacker::getID(bb,idstart),sdhcal::PMRUnpacker::getGTC(bb,idstart),sdhcal::PMRUnpacker::getAbsoluteBCID(bb,idstart),sdhcal::PMRUnpacker::getLastTriggerBCID(bb,idstart),sdhcal::PMRUnpacker::getBCID(bb,idstart),sdhcal::PMRUnpacker::getLastTriggerBCID(bb,idstart)-sdhcal::PMRUnpacker::getBCID(bb,idstart));
		  //getchar();
		  if (!_initialised)
		    {
		      _theEvent.init(_run,_event,sdhcal::PMRUnpacker::getGTC(bb,idstart),sdhcal::PMRUnpacker::getAbsoluteBCID(bb,idstart));
		      _initialised=true;
		    }
		  

		  sdhcal::PMRPtr* d= new sdhcal::PMRPtr(&_buf[_idx],b.payloadSize());
		  _theEvent.difList().push_back(d);
		  _idx+=b.payloadSize();
		}


	      if (_detId==130)
		{
		  uint8_t* bb=&_buf[_idx];
		  uint32_t* ibuf=(uint32_t*) b.payload();
		  
		  for (int i=0;i<7;i++)
		    {
		      printf("%d ",ibuf[i]);
		    }
		  uint32_t nch=ibuf[6];
		  printf("\n channels -> %d \n",nch);
		  _mezzanine=ibuf[4];
		  _difId=(ibuf[5]>>24)&0xFF;
		  _gtc=ibuf[1];
		  
		  if (!_initialised)
		    {
		      _theEvent.init(_run,_event,_bxId,_gtc);
		      _initialised=true;
		    }
		  INFO_PRINTF("\t \t \t %d %d GTC %d NCH %d \n",_mezzanine,_difId,_gtc,nch);

		  
		  if (ibuf[6]>=0)
		    {
		      uint8_t* cbuf=( uint8_t*)&ibuf[7];
		      bool tfound=false;
		      for (int i=0;i<nch;i++)
			{
#define DUMPCHANSN	
#ifdef DUMPCHANS			  
			  for (int j=0;j<6;j++)
			     INFO_PRINTF("\t %.2x ",cbuf[i*6+j]);
			   INFO_PRINTF("\n");
#endif
			  memcpy(&_eventChannel[_eventChannels],&cbuf[6*i],6*sizeof(uint8_t));
			  lydaq::TdcChannel ca((uint8_t*) &_eventChannel[_eventChannels],_difId&0xFF);
			  //c.dump();
			 
			  _eventChannels++;
			  if (ca.channel()==0)
			    {
			      _0coarse[_difId&0xFF]=ca.coarse();
			      _0fine[_difId&0xFF]=ca.fine();
			    }
			  //_mezMap[_difId].push_back(c);
			  ca.setZero(_0coarse[_difId&0xFF], _0fine[_difId&0xFF]);
			  //  ca.dump();
			  //std::cout<<_difId<<" "<<(int) _0coarse[_difId&0xFF]<<" "<< (int) _0fine[_difId&0xFF]<<std::endl;
			  _theEvent.tdcChannels().push_back(ca);
			  
			}
		      //if (nch>0) getchar();
#ifdef DUMPCHANS
		      if (nch>0) getchar();
#endif
		    }

		}

	    }

	  /// Timing
	  bool oneselected=false;
	  if (_theEvent.isDifType())
	    {
	      TH1* hacqtim= _rh->GetTH1("/BR/AcquistionTime");
	      TH1* hrealtim= _rh->GetTH1("/BR/RealTime");
	      TH1* hmt= _rh->GetTH1("/BR/MaxTime");
	      TH1* hc= _rh->GetTH1("/BR/Count");
	      TH1* hcs= _rh->GetTH1("/BR/PadCountSelected");
	      TH1* hfs= _rh->GetTH1("/BR/FrameCountSelected");
	      TH2* hxy= _rh->GetTH2("/BR/XY");
	      TH2* hxys= _rh->GetTH2("/BR/XYSelected");
	      
	      if (hacqtim==NULL)
		{
		  hacqtim=_rh->BookTH1("/BR/AcquistionTime",10000.,0.,10000.);
		  hrealtim=_rh->BookTH1("/BR/RealTime",10000.,0.,4.);
		  hmt=_rh->BookTH1("/BR/MaxTime",10000.,0.,5.);
		  hc=_rh->BookTH1("/BR/Count",100.,-0.1,99.9);
		  hcs=_rh->BookTH1("/BR/PadCountSelected",100.,-0.1,99.9);
		  hfs=_rh->BookTH1("/BR/FrameCountSelected",100.,-0.1,99.9);
		  hxy=_rh->BookTH2("/BR/XY",64,0.1,64.1,35,0.1,35.1);
		  hxys=_rh->BookTH2("/BR/XYSelected",48,0.1,48.1,24,0.1,24.1);
		}
	      hrealtim->Fill(treal);
	      hc->Fill(1.);
	      _theEvent.tFrame().clear();
	      _theEvent.tCount().clear();
	      std::map<uint32_t,std::vector<std::pair<sdhcal::PMRPtr*,uint32_t> > >::iterator ifm=_theEvent.tFrame().end();
	      std::map<uint32_t,std::bitset<64> >::iterator im=_theEvent.tCount().end();
	      ROOT::Math::XYZPoint pt;
	      float mt=0;

	      for (std::vector<sdhcal::PMRPtr*>::iterator it = _theEvent.difList().begin();it!=_theEvent.difList().end();it++) // Loop on DIF
		{
		  sdhcal::PMRPtr* d = (*it);
		  uint32_t chid= _geo->difInfo(d->getID()).chamber;
		  fprintf(stderr,"CHAMBER %d %d \n",d->getID(),chid);
		  // LMTest      uint32_t bc = rint(f->getBunchCrossingTime()/DCBufferReader::getDAQ_BC_Period());
		  //uint32_t chid = 1;
		  uint32_t window=2;
    
		  // Loop on Frames
		  bool trigged=false;
		  uint32_t nps=0,nfs=0;

		  for (uint32_t ifra=0;ifra<d->getNumberOfFrames();ifra++)
		    {
		      int32_t bc=d->getFrameTimeToTrigger(ifra);
		      bc=d->getFrameBCID(ifra);
		      //printf(" Frame BCID ch %d dif %d asic %d framebc %d difbc %d delta %d \n",chid,d->getID(),d->getFrameAsicHeader(ifra),d->getFrameBCID(ifra),d->getBCID(),bc);
		      if ((bc*2E-7)>mt) mt=bc*2E-7;
		      fflush(stdout);
		      hacqtim->Fill(bc);
		      std::stringstream sas,sdif,sch;
		      sas.str("");
		      sdif.str("");
		      sch.str("");
		      sch<<"/BR/CH_"<<chid<<"/";
		      sas<<"/BR/DIF_"<<d->getID()<<"/A_"<<d->getFrameAsicHeader(ifra)<<"/";
		      sdif<<"/BR/DIF_"<<d->getID()<<"/";
		      TH1* has=_rh->GetTH1(sdif.str()+"Asics");
		      TH1* hch=_rh->GetTH1(sas.str()+"Pads");
		      TH1* hchs=_rh->GetTH1(sas.str()+"PadsSelected");
		      TH1* hacs=_rh->GetTH1(sch.str()+"AcquisitionTime");
		      TH2* hchxy=_rh->GetTH2(sch.str()+"XY");
		      
		      if (hch==NULL)
			{
			  has=_rh->BookTH1(sdif.str()+"Asics",48,0.,48.);
			  hch=_rh->BookTH1(sas.str()+"Pads",64,0.,64.);
			  hchs=_rh->BookTH1(sas.str()+"PadsSelected",64,0.,64.);
			  hchxy=_rh->BookTH2(sch.str()+"XY",64,0.1,64.1,35,0.1,35.1);
			  hacs=_rh->BookTH1(sch.str()+"AcquisitionTime",10000.,0.,10000.);
			}
		      bool selected=(bc>62 && bc<70);
		      std::bitset<64> ph;
		      ph.reset();
		      for (int ipad=0;ipad<64;ipad++)
			{
			  
			  if (d->getFrameLevel(ifra,ipad,0) || d->getFrameLevel(ifra,ipad,1))
			    {
			      ph.set(ipad,1);
			      _geo->convert(d->getID(),d->getFrameAsicHeader(ifra),ipad,&pt);
			      //printf("Point %d %d %d %f %f \n",d->getID(),d->getFrameAsicHeader(ifra),ipad,pt.X(),pt.Y());
			      hxy->Fill(pt.X(),pt.Y());
			      hchxy->Fill(pt.X(),pt.Y());
			      hch->Fill(ipad*1.);
			      has->Fill(d->getFrameAsicHeader(ifra)*1.);
			      hacs->Fill(bc);
			      if (selected)
				{
				  hxys->Fill(pt.X(),pt.Y());


				  hchs->Fill(ipad*1.);
				  nps++;
				}
			    }
			}
		      if (selected) {
			if (nps>0)
			  {
			    printf(" Frame BCID %d %d %d %d %d \n",chid,d->getID(),d->getFrameAsicHeader(ifra),d->getFrameBCID(ifra),bc);
			    std::cout<<ph<<std::endl;
			    fflush(stdout);
			    //getchar();
			  }
			nfs++;
		      }
		      //printf(" Frame BCID %d %d %d %d \n",d->getID(),d->getFrameAsicHeader(ifra),d->getFrameBCID(ifra),bc);
		      trigged=trigged||selected;
		      oneselected=oneselected||selected;
		      bool found=false;
		      for (int dt=-2;dt<=2;dt++)
			{
			  uint32_t tc=bc+dt;
			  //if (selected) printf("%d %d %d \n",bc,dt,tc);
			  im=_theEvent.tCount().find(tc);
	  

			  if (im!=_theEvent.tCount().end())
			    {
			      im->second.set(chid,1);

			      //if (selected) std::cout<<tc<<"CH "<<chid<<" "<<im->second<<std::endl;
			      ifm=_theEvent.tFrame().find(tc);
			      ifm->second.push_back(std::pair<sdhcal::PMRPtr*,uint32_t>(d,ifra));
			      found=true;
			      // std::cout<<bc<<" add "<<_theEvent.tFrame()[tc].size()<<std::endl;
			      break;
			    }
			}
		      
		      if (found) continue;
		      std::bitset<64> v(0);
		      v.set(chid,1);
		      std::pair<uint32_t,std::bitset<64> > p(bc,v);
		      //if (selected) std::cout<<bc<<"CH0 "<<chid<<" "<<p.second<<std::endl;
		      _theEvent.tCount().insert(p);
		      std::vector<std::pair<sdhcal::PMRPtr*,uint32_t> > vf;
		      vf.clear();
		      vf.push_back(std::pair<sdhcal::PMRPtr*,uint32_t>(d,ifra));
		      std::pair<uint32_t,std::vector<std::pair<sdhcal::PMRPtr*,uint32_t> > > pf(bc,vf);
		      _theEvent.tFrame().insert(pf);
		      //std::cout<<bc<<" create "<<_theEvent.tFrame()[bc].size()<<std::endl;		      
		    }
		  

		  //getchar();
		  if (trigged) hc->Fill(2.);
		  if (nps>0) hcs->Fill(nps*1.0);
		  if (nfs>0) hfs->Fill(nfs*1.0);
		}
	      hmt->Fill(mt);
	      

	    }
	  printf(" TCOUNT %ld TMAP %ld \n",_theEvent.tCount().size(),_theEvent.tFrame().size());
	  bool coinc=false;
	  for (auto x:_theEvent.tCount())
	    {

	      if (x.second.count()>2)
		{
		  fprintf(stderr," bcid %d cnt %ld \n",x.first,x.second.count());
		  auto tf=_theEvent.tFrame()[x.first];
		 for (auto it=tf.begin();it!=tf.end();it++)
		   {
		     std::cout<<it->first->getID()<<" "<<it->first->getFrameTimeToTrigger(it->second)<<std::endl;
		       }
		coinc=true;
		}
	    }
	  //getchar();
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
