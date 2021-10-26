#pragma once
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
#define DU_END_OF_DIF        0xA0


#define DU_START_OF_FRAME    0xB4
#define DU_END_OF_FRAME      0xA3

#define LIROC_ID_SHIFT 1
#define LIROC_NBASIC_SHIFT 2
#define LIROC_FORMAT_SHIFT 3
#define LIROC_GTC_SHIFT 4
#define LIROC_ABCID_SHIFT 7
#define LIROC_BCID_SHIFT 13
#define LIROC_LTRG_SHIFT 16
#define LIROC_HEADER_SHIFT 20

#define DU_ID_SHIFT     LIROC_ID_SHIFT
#define DU_GTC_SHIFT    LIROC_GTC_SHIFT
#define DU_ABCID_SHIFT  LIROC_ABCID_SHIFT
#define DU_BCID_SHIFT   LIROC_BCID_SHIFT
#define DU_LTRG_SHIFT   LIROC_LTRG_SHIFT

#define DU_FRAME_CHANNEL_SHIFT 3
#define DU_FRAME_BCID_SHIFT        0
#define DU_FRAME_SIZE              4

#define MAX_NUMBER_OF_HIT 4000000
namespace liroc
{
  class Unpacker {
  public:
    static uint32_t getStartOfLIROC(unsigned char* cbuf,uint32_t size_buf,uint32_t start=92);
    static uint32_t getID(unsigned char* cb,uint32_t idx=0);
    static uint32_t getGTC(unsigned char* cb,uint32_t idx=0);
    static unsigned long long getAbsoluteBCID(unsigned char* cb,uint32_t idx=0);
    static uint32_t getBCID(unsigned char* cb,uint32_t idx=0);
    static uint32_t getLastTriggerBCID(unsigned char* cb,uint32_t idx=0);
    static uint32_t getFrameChannel(unsigned char* framePtr);
    static uint32_t getFrameBCID(unsigned char* framePtr);
//static uint8_t getFrameData(unsigned char* framePtr,uint32_t ip);



  
    static uint32_t getFramePtr(std::vector<unsigned char*> &vFrame,uint32_t max_size,unsigned char* cb,uint32_t idx=0) throw (std::string);
    static void dumpFrameOld(unsigned char* buf);

    static unsigned long swap_bytes(unsigned int n,unsigned char* buf);
  };

  class Ptr
  {
  public:
    Ptr(unsigned char* p,uint32_t max_size) : theSize_(max_size),theLIROC_(p)
    {
      theFrames_.clear();
      try
	{
	  liroc::Unpacker::getFramePtr(theFrames_,theSize_,theLIROC_);
	}
      catch (std::string e)
	{
	  std::cout<<"LIROC "<<getID()<<e<<std::endl;
	}
    }
    inline unsigned char* getPtr(){return theLIROC_;}
    inline std::vector<unsigned char*>& getFramesVector(){return theFrames_;}
    inline  uint32_t getID(){return liroc::Unpacker::getID(theLIROC_);}
    inline  uint32_t getGTC(){return liroc::Unpacker::getGTC(theLIROC_);}
    inline  unsigned long long getAbsoluteBCID(){return liroc::Unpacker::getAbsoluteBCID(theLIROC_);}
    inline  uint32_t getBCID(){return liroc::Unpacker::getBCID(theLIROC_);}
    inline uint32_t getNumberOfFrames(){return theFrames_.size();}
    inline unsigned char* getFramePtr(uint32_t i){return theFrames_[i];}
    inline uint32_t getFrameChannel(uint32_t i){return liroc::Unpacker::getFrameChannel(theFrames_[i]);}
    inline uint32_t getFrameBCID(uint32_t i){return liroc::Unpacker::getFrameBCID(theFrames_[i]);}
					//inline uint32_t getFrameData(uint32_t i,uint32_t iword){return liroc::Unpacker::getFrameData(theFrames_[i],iword);}
    void dumpInfo()
    {
      printf("LIROC %d  GTC %d ABCID %lld BCID %d \n",
	     getID(),
	     getGTC(),
	     getAbsoluteBCID(),
	     getBCID());
      printf("Found  %d Frames \n",(int) theFrames_.size());
    }

  private:
    uint32_t theSize_;
    unsigned char* theLIROC_;
    std::vector<unsigned char*> theFrames_;

  };

  
};

