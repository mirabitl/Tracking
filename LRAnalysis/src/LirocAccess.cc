#include "LirocAccess.hh"

#include <sstream>
#include <iostream>
#include <stdio.h>
#include <bitset>
#include <string.h>

using namespace std;


uint32_t liroc::Unpacker::getStartOfLIROC(unsigned char* cbuf,uint32_t size_buf,uint32_t start)
{
  uint32_t id0=0;

  //for (int i=0;i<128;i++) {printf("%02x ",cbuf[i]); if ((i+1)%32==0) 	printf("\n");};	  printf("\n");
  
  
  for (uint32_t i=start;i<size_buf;i++)
    {
      if (cbuf[i]!=DU_START_OF_DIF) continue;
      id0=i;
      //if (cbuf[id0+DU_ID_SHIFT]>0xFF) continue; 
      break;
    }
  //printf("DIF %d Found at %d \n",cbuf[id0+1],id0);
  return id0;
}


uint32_t liroc::Unpacker::getID(unsigned char* cb,uint32_t idx)
{
  return cb[idx+DU_ID_SHIFT];
}


uint32_t liroc::Unpacker::getGTC(unsigned char* cb,uint32_t idx)
{
  return (cb[idx+DU_GTC_SHIFT]<<16)+(cb[idx+DU_GTC_SHIFT+1]<<8)+cb[idx+DU_GTC_SHIFT+2];
}

unsigned long long liroc::Unpacker::getAbsoluteBCID(unsigned char* cb,uint32_t idx)
{
  unsigned long long Shift=16777216ULL;//to shift the value from the 24 first bits
  unsigned long long LBC= ((cb[idx+DU_ABCID_SHIFT]<<8) | (cb[idx+DU_ABCID_SHIFT+1]) )*Shift+( (cb[idx+DU_ABCID_SHIFT+2]<<24) | (cb[idx+DU_ABCID_SHIFT+3]<<16)| (cb[idx+DU_ABCID_SHIFT+4]<<8) | (cb[idx+DU_ABCID_SHIFT+5]));
  return LBC;
}

uint32_t liroc::Unpacker::getBCID(unsigned char* cb,uint32_t idx)
{
  return (cb[idx+DU_BCID_SHIFT]<<16)+(cb[idx+DU_BCID_SHIFT+1]<<8)+cb[idx+DU_BCID_SHIFT+2];
}
uint32_t liroc::Unpacker::getLastTriggerBCID(unsigned char* cb,uint32_t idx)
{
  return (cb[idx+DU_LTRG_SHIFT+1]<<16)+(cb[idx+DU_LTRG_SHIFT+2]<<8)+cb[idx+DU_LTRG_SHIFT+3];
}





uint32_t liroc::Unpacker::getFrameChannel(unsigned char* framePtr)
{
  return (framePtr[DU_FRAME_CHANNEL_SHIFT]>>2)&0x3F;
}
uint32_t liroc::Unpacker::getFrameBCID(unsigned char* framePtr)
{

  // for (int i=0;i<4;i++)
  //   fprintf(stderr,"%.2x ",framePtr[i]);
  // fprintf(stderr,"\n ");
  
  uint32_t igray=(framePtr[DU_FRAME_BCID_SHIFT+2]<<16)+(framePtr[DU_FRAME_BCID_SHIFT+1]<<8)+framePtr[DU_FRAME_BCID_SHIFT];
  return igray;
}




uint32_t liroc::Unpacker::getFramePtr(std::vector<unsigned char*> &vFrame,uint32_t max_size,unsigned char* cb,uint32_t idx) throw (std::string)
{
  uint32_t fshift=idx+LIROC_HEADER_SHIFT;


  do
    {
      // printf("fshift %d and %d \n",fshift,max_size);

      /// Pas de b4 if (cb[fshift]==DU_START_OF_FRAME) fshift++;
      /// Pas de  A3 if (cb[fshift]==DU_END_OF_FRAME) {fshift++;continue;}
      //uint32_t header =liroc::Unpacker::getFrameAsicHeader(&cb[fshift]);
      /// Ilogique if (header == DU_END_OF_FRAME) return (fshift+2);

      //// INVERT HEADER NOt OK
      /**
      uint8_t headerb=0;
      for (int i=0;i<8;i++)
	{
	  if (header & (1<<i))
	    headerb|=(1<<(7-i));
	}
      // fprintf(stderr,"Header found %x %x Shift %d \n",header,headerb,fshift);
      setFrameAsicHeader(&cb[fshift],headerb);

      
      header=headerb;
      */
      // A VOIR if (cb[fshift]==DU_END_OF_DIF) return fshift;
      // if (header<1 || (header>48 && header!=129)) // correction d'antoine pour la BIF
      // 	{
      // 	  std::stringstream s("");
      // 	  s<<header<<" Header problem "<<fshift<<std::endl;
      // 	  throw  s.str();
      // 	  return fshift;

      // 	}
      vFrame.push_back(&cb[fshift]);fshift+=DU_FRAME_SIZE;
      if (fshift>max_size) 
	{
	  printf("fshift %d exceed %d \n",fshift,max_size);
	  return fshift;
	}
      //printf("%x \n",cb[fshift]);
      /// Pas de A3 if (cb[fshift]==DU_END_OF_FRAME) fshift++;

    } while (1);
}
unsigned long liroc::Unpacker::swap_bytes(unsigned int n,unsigned char* buf)
{
  unsigned char Swapped[4];
  memset(Swapped,0,4);
  for (unsigned int i=0;i<n;i++)
    Swapped[i] = buf[n-1-i];
  

  unsigned long *temp =(unsigned long*) &Swapped;

  return (*temp);
}
void liroc::Unpacker::dumpFrameOld(unsigned char* buf)
{
  bool PAD[128];
  bool l0[64];
  bool l1[64];
  unsigned short un = 1;

  for(int ip= 0; ip<128; ip++){PAD[ip]=0;} //init PADs
  uint32_t idx1=4;
  for(int ik=0;ik<4;ik++)
    {
      
      unsigned long PadEtat= swap_bytes(4,&buf[idx1]);
      idx1+=4;
      
      for(int e=0;e<32;e++)
	{	
	  PAD[((3-ik)*32)+(31-e)]=PadEtat & un; //binary operation
	  PadEtat=PadEtat>>1;	//décalage des bit de 1
	}
    }
  // fill bool arrays
  for(int p=0; p<64;p++)
    {
      l0[p]=(bool)PAD[(2*p)]; //_Lev0 (PAD paire)
      l1[p]=(bool)PAD[(2*p)+1]; //_Lev1 (PAD impaires)
      
    }
  std::bitset<64> bs0(0);
  std::bitset<64> bs1(0);
  for (uint32_t ip=0;ip<64;ip++) {bs0.set(ip,l0[ip]);bs1.set(ip,l1[ip]);}

  std::cout<<"\t \t"<<bs0<<std::endl;
  std::cout<<"\t \t"<<bs1<<std::endl;      
}
