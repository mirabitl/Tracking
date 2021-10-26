#pragma once
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstdint>
#include <arpa/inet.h>

#define RUN_HEADER_START 0x5b
#define RUN_HEADER_END 0x5d
#define MAX_RUN_HEADER_SIZE 200
#define READOUT_START 0x28
#define READOUT_END 0x29

class FEB_v2_data_reader
{
public:
  FEB_v2_data_reader() {}
  ~FEB_v2_data_reader() { data_file.close(); free(buffer);}
  void close(){data_file.close();}
  bool openDataFile(const char* file_name);
  const uint32_t* runHeader() {return run_header;} // might put decoding run header later
  //run header decoding
  uint32_t get_fc7_general_register_value() {return run_header[0];}
  unsigned int get_acquisition_mode() {return (run_header[0]>>6)&0x7;}
  unsigned int get_nb_frame_acquisition() {return run_header[1];} //for mode TDC_Count
  unsigned int get_time_window() {return run_header[2];} //for mode time_window
  unsigned int get_circular_buffer_size() {return run_header[3];} //for mode circular_buffer trigger internal or external
  unsigned int get_FC7_transmission_register_value() {return run_header[4];}
  unsigned int get_N_orbits() {return (run_header[4]>>21)&0xff;}
  float get_sensor1_temperature() {return ((run_header[5]>>18)&0x1ff)/2.0;}
  float get_sensor2_temperature() {return ((run_header[5]>>9)&0x1ff)/2.0;}
  float get_sensor3_temperature() {return ((run_header[5])&0x1ff)/2.0;}
  float get_sensor4_temperature() {return ((run_header[6]>>9)&0x1ff)/2.0;}
  float get_sensor5_temperature() {return ((run_header[6])&0x1ff)/2.0;}
  uint32_t get_fpga0_strip_TDCchannel_enable_mask() {return run_header[7];}
  uint32_t get_fpga1_strip_TDCchannel_enable_mask() {return run_header[8];}
  uint32_t get_fpga2_strip_TDCchannel_enable_mask() {return run_header[9];}
  unsigned int get_fpga0_BC0_TDC_enable() {return ((run_header[10]>>4)&1);}
  unsigned int get_fpga1_BC0_TDC_enable() {return ((run_header[10]>>2)&1);}
  unsigned int get_fpga2_BC0_TDC_enable() {return ((run_header[10])&1);}
  unsigned int get_fpga0_Resync_TDC_enable() {return ((run_header[10]>>5)&1);}
  unsigned int get_fpga1_Resync_TDC_enable() {return ((run_header[10]>>3)&1);}
  unsigned int get_fpga2_Resync_TDC_enable() {return ((run_header[10]>>1)&1);}
  unsigned int get_fpga0_Global_TDC_enable() {return ((run_header[10]>>6)&1);}
  unsigned int get_fpga1_Global_TDC_enable() {return ((run_header[10]>>7)&1);}
  unsigned int get_fpga2_Global_TDC_enable() {return ((run_header[10]>>8)&1);}
  //
  unsigned int runNumber() {return run_number;}
  unsigned int eventNumber() {return event_number;}
  bool nextReadout();
  bool nextGBTFrame();
  bool nextTDCFrame();
  bool checkGBTFramePadding();
  uint32_t BC0id() { uint32_t *p= (uint32_t*)(current_GBT_frame+14); return ntohl(*p);} 
  char validFrameStatus() {return ( *(current_GBT_frame+19)) & 0x7;}
  unsigned int frameIsValid() { unsigned int shift=2-current_TDC_frame; return (validFrameStatus()>>shift)&1;}
  unsigned int fpga() {return (unsigned int) (((*currentTDCframe())&0xff)>>6);}
  unsigned int TDC_channel() {return (unsigned int) ((*currentTDCframe())&0x3f);}
  unsigned int TDC_value() { uint32_t* pval=(uint32_t*) currentTDCframe(); unsigned int val = ntohl(*pval); return val&0xffffff;}
private:
  char* currentTDCframe() {return current_GBT_frame+(20+4*current_TDC_frame);}
  void revertBuffer(uint32_t* begin, uint32_t*end);
  std::ifstream data_file;
  uint32_t run_number;
  uint32_t run_header[MAX_RUN_HEADER_SIZE];
  uint32_t event_number;
  unsigned int buffer_size=0;
  char* buffer=nullptr;
  char* buffer_end=nullptr;
  char* current_GBT_frame=nullptr;
  unsigned int current_TDC_frame=0;
  bool newfile=true;
  void debugPrint(); 
};

void FEB_v2_data_reader::debugPrint()
{
  std::cout << "---> DEBUG : " << buffer_size << " buffer start " << (void*)buffer << " buffer end " <<  (void*)buffer_end << " GBT frame " << (void*)current_GBT_frame << " TDC frame " <<  current_TDC_frame << std::endl;
}

bool FEB_v2_data_reader::openDataFile(const char* file_name)
{
  data_file.close();
  data_file.open(file_name);
  newfile=true;
  if (not data_file.good()) { std::cerr << "Can't open file " << file_name << std::endl; return false;}
  char key;
  data_file.read(&key,1);
  if (key == READOUT_START) {data_file.unget(); return true;}
  if (key != RUN_HEADER_START) { std::cerr << "File " << file_name << " doesn't start with run header key nor event readout key" << std::endl; return false;}
  data_file.read((char*) &run_number, 4);
  std::cout << "run Number is " << run_number << std::endl;
  uint32_t nint=0;
  data_file.read((char*) &nint, 4);
  //std::cout << "There are "<< nint << " words to read" << std::endl;
  data_file.read((char*) run_header, nint*4);
  //for (uint32_t i=0; i<nint; ++i) std::cout << run_header[i] << " // ";
  //std::cout << std::endl;
  data_file.read(&key,1); 
  if (key != RUN_HEADER_END) { std::cerr << "File " << file_name << " corrupted : missing run header end key" << std::endl; return false;}
  return true;
}

void FEB_v2_data_reader::revertBuffer(uint32_t* begin, uint32_t*end)
{
  for (uint32_t *it=begin; it<end; ++it)
    (*it)=ntohl(*it); //big-endian to little-endian conversion
}

bool FEB_v2_data_reader::nextReadout()
{
  newfile=false;
  char key=0;
  data_file.read(&key,1);
  if (data_file.eof())
    {
      std::cout << "End of file reached" << std::endl;
      return false;
    }
  if (key != READOUT_START)
    {
      std::cerr << " Readout start key not found " << std::endl;
      return false;
    }
  data_file.read((char*) &event_number,4);
  //std::cout << "event number is " << event_number << std::endl;
  uint32_t nWords=0;
  data_file.read((char*) &nWords,4);
  //std::cout << "word to read = " << nWords << std::endl;
  unsigned int sizeToRead = nWords*4;
  if (sizeToRead> buffer_size)
    {
      free(buffer);
      buffer = (char*) malloc(sizeToRead);
      buffer_size = sizeToRead;
    }
  data_file.read(buffer,sizeToRead);
  buffer_end=buffer+sizeToRead;
  revertBuffer((uint32_t*) buffer, (uint32_t*) buffer_end); 
  current_GBT_frame=buffer;
  data_file.read(&key,1);
  if (key != READOUT_END)
    {
      std::cerr << " Readout end key not found , got " << key << "("<<int(key)<<")" << std::endl;
      return false;
    }  
  return true;
}

bool FEB_v2_data_reader::nextGBTFrame()
{
  if (newfile) return nextReadout();
  current_GBT_frame=current_GBT_frame+32; //a readout is 256 bits
  if (current_GBT_frame>=buffer_end) return nextReadout();
  return true;
}


bool FEB_v2_data_reader::checkGBTFramePadding()
{
  char* p=current_GBT_frame;
  bool OK=true;
  for (int i=0; i<14; ++i)
    {
      OK = (OK and ((*p) ==0));
      ++p;
    }
  return OK;
}

bool FEB_v2_data_reader::nextTDCFrame()
{
  //debugPrint();
  ++current_TDC_frame;
  if (newfile or current_TDC_frame==3)
    {
      bool OK=nextGBTFrame();
      if (not OK) return false;
      current_TDC_frame=0;
    }
  if ( frameIsValid() )
    return true;
  else
    return nextTDCFrame();
}
