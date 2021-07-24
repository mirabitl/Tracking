#ifndef _SDHCALDIFACCESS_HH
#define _SDHCALDIFACCESS_HH
#include <string>
#include <vector>
#include <stdint.h>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <bitset>
#include <string.h>


#define DU_DATA_FORMAT_VERSION 13

#define DU_START_OF_DIF      0xB0
#define DU_START_OF_DIF_TEMP 0xBB
#define DU_END_OF_DIF        0xA0
#define DU_START_OF_LINES    0xC4
#define DU_END_OF_LINES      0xD4

#define DU_START_OF_FRAME    0xB4
#define DU_END_OF_FRAME      0xA3


#define DU_ID_SHIFT      1
#define DU_DTC_SHIFT     2
#define DU_GTC_SHIFT    10
#define DU_ABCID_SHIFT  14
#define DU_BCID_SHIFT   20
#define DU_LINES_SHIFT  23
#define DU_TASU1_SHIFT  24
#define DU_TASU2_SHIFT  28
#define DU_TDIF_SHIFT   32

#define DU_FRAME_ASIC_HEADER_SHIFT 0
#define DU_FRAME_BCID_SHIFT        1
#define DU_FRAME_DATA_SHIFT        4
#define DU_FRAME_SIZE              20

#define MAX_NUMBER_OF_HIT 4000000
namespace sdhcal
{
  class DIFUnpacker {
  public:
    static unsigned long long GrayToBin(unsigned long long n);
    static uint32_t getStartOfDIF(unsigned char* cbuf,uint32_t size_buf,uint32_t start=92);
    static uint32_t getID(unsigned char* cb,uint32_t idx=0);
    static uint32_t getDTC(unsigned char* cb,uint32_t idx=0);
    static uint32_t getGTC(unsigned char* cb,uint32_t idx=0);
    static unsigned long long getAbsoluteBCID(unsigned char* cb,uint32_t idx=0);
    static uint32_t getBCID(unsigned char* cb,uint32_t idx=0);
    static uint32_t getLines(unsigned char* cb,uint32_t idx=0);
    static bool hasLine(uint32_t line,unsigned char* cb,uint32_t idx=0);
    static uint32_t getTASU1(unsigned char* cb,uint32_t idx=0);
    static uint32_t getTASU2(unsigned char* cb,uint32_t idx=0);
    static uint32_t getTDIF(unsigned char* cb,uint32_t idx=0);
    static bool hasTemperature(unsigned char* cb,uint32_t idx=0);
    static bool hasAnalogReadout(unsigned char* cb,uint32_t idx=0);

    static uint32_t getFrameAsicHeader(unsigned char* framePtr);
    static uint32_t getFrameBCID(unsigned char* framePtr);
    static uint8_t getFrameData(unsigned char* framePtr,uint32_t ip);

    static bool getFramePAD(unsigned char* framePtr,uint32_t ip);
    static bool getFrameLevel(unsigned char* framePtr,uint32_t ip,uint32_t level);
  
    static uint32_t getAnalogPtr(std::vector<unsigned char*> &vLines,unsigned char* cb,uint32_t idx=0) throw (std::string);
    static uint32_t getFramePtr(std::vector<unsigned char*> &vFrame,std::vector<unsigned char*> &vLines,uint32_t max_size,unsigned char* cb,uint32_t idx=0) throw (std::string);
    static void dumpFrameOld(unsigned char* buf);
    static unsigned long swap_bytes(unsigned int n,unsigned char* buf);
  };

  class DIFPtr
  {
  public:
    DIFPtr(unsigned char* p,uint32_t max_size) : theSize_(max_size),theDIF_(p)
    {
      theFrames_.clear();theLines_.clear();
      try
	{
	  sdhcal::DIFUnpacker::getFramePtr(theFrames_,theLines_,theSize_,theDIF_);
	}
      catch (std::string e)
	{
	  std::cout<<"DIF "<<getID()<<" T ? "<<hasTemperature()<<" " <<e<<std::endl;
	}
    }
    inline unsigned char* getPtr(){return theDIF_;}
    inline std::vector<unsigned char*>& getFramesVector(){return theFrames_;}
    inline std::vector<unsigned char*>& getLinesVector(){return theLines_;}
    inline  uint32_t getID(){return sdhcal::DIFUnpacker::getID(theDIF_);}
    inline  uint32_t getDTC(){return sdhcal::DIFUnpacker::getDTC(theDIF_);}
    inline  uint32_t getGTC(){return sdhcal::DIFUnpacker::getGTC(theDIF_);}
    inline  unsigned long long getAbsoluteBCID(){return sdhcal::DIFUnpacker::getAbsoluteBCID(theDIF_);}
    inline  uint32_t getBCID(){return sdhcal::DIFUnpacker::getBCID(theDIF_);}
    inline  uint32_t getLines(){return sdhcal::DIFUnpacker::getLines(theDIF_);}
    inline  bool hasLine(uint32_t line){return sdhcal::DIFUnpacker::hasLine(line,theDIF_);}
    inline  uint32_t getTASU1(){return sdhcal::DIFUnpacker::getTASU1(theDIF_);}
    inline  uint32_t getTASU2(){return sdhcal::DIFUnpacker::getTASU2(theDIF_);}
    inline  uint32_t getTDIF(){return sdhcal::DIFUnpacker::getTDIF(theDIF_);}
    inline  float getTemperatureDIF(){return 0.508*getTDIF()-9.659;}
    inline  float getTemperatureASU1(){return (getTASU1()>>3)*0.0625;}
    inline  float getTemperatureASU2(){return (getTASU2()>>3)*0.0625;}
    inline  bool hasTemperature(){return sdhcal::DIFUnpacker::hasTemperature(theDIF_);}
    inline  bool hasAnalogReadout(){return sdhcal::DIFUnpacker::hasAnalogReadout(theDIF_);}
    inline uint32_t getNumberOfFrames(){return theFrames_.size();}
    inline unsigned char* getFramePtr(uint32_t i){return theFrames_[i];}
    inline uint32_t getFrameAsicHeader(uint32_t i){return sdhcal::DIFUnpacker::getFrameAsicHeader(theFrames_[i]);}
    inline uint32_t getFrameBCID(uint32_t i){return sdhcal::DIFUnpacker::getFrameBCID(theFrames_[i]);}
    inline uint32_t getFrameTimeToTrigger(uint32_t i){return getBCID()-getFrameBCID(i);}
    inline bool getFrameLevel(uint32_t i,uint32_t ipad,uint32_t ilevel){return sdhcal::DIFUnpacker::getFrameLevel(theFrames_[i],ipad,ilevel);}
    inline uint32_t getFrameData(uint32_t i,uint32_t iword){return sdhcal::DIFUnpacker::getFrameData(theFrames_[i],iword);}
    void dumpDIFInfo()
    {
      printf("DIF %d DTC %d GTC %d ABCID %lld BCID %d Lines %d Temperature %d \n",
	     getID(),
	     getDTC(),
	     getGTC(),
	     getAbsoluteBCID(),
	     getBCID(),
	     getLines(),
	     hasTemperature());

      if (hasTemperature())
	printf("T: ASU1 %d %f ASU2 %d %f DIF %d  %f \n",getTASU1(),getTemperatureASU1(),getTASU2(),getTemperatureASU2(),getTDIF(),getTemperatureDIF());
      printf("Found %d Lines and %d Frames \n",(int) theLines_.size(),(int) theFrames_.size());
    }

  private:
    uint32_t theSize_;
    unsigned char* theDIF_;
    std::vector<unsigned char*> theFrames_;
    std::vector<unsigned char*> theLines_;
  };

  
};
#endif
