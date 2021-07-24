#pragma once
#include "FEB_v2_reader.cpp"
#include <utility>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>

typedef unsigned int TDC_value_type;
typedef unsigned int channel_type;
typedef unsigned int fpga_type;
typedef uint32_t BC0_id_type;

typedef std::pair<channel_type,TDC_value_type> TDC_data_type;
typedef std::vector<TDC_data_type> one_FPGA_data_type;
typedef std::array<one_FPGA_data_type,3> three_FPGA_data_type;
typedef std::map<BC0_id_type,three_FPGA_data_type> FEB_frame_data_type;

bool decode(const char* file_name, FEB_frame_data_type& FEB_data)
{
  FEB_v2_data_reader f;
  bool OK=f.openDataFile(file_name);
  if (not OK) return false;
  while (f.nextTDCFrame())
    {
      FEB_data[f.BC0id()][f.fpga()].push_back(TDC_data_type(f.TDC_channel(),f.TDC_value()));
    }
  for (FEB_frame_data_type::iterator it=FEB_data.begin(); it!=FEB_data.end(); ++it)
    for (three_FPGA_data_type::iterator it2=it->second.begin(); it2 != it->second.end(); ++it2)
      std::sort(it2->begin(),it2->end());
  return true;
}

void print(const FEB_frame_data_type& FEB_data, std::ostream& out=std::cout)
{
  for (FEB_frame_data_type::const_iterator it=FEB_data.begin(); it!=FEB_data.end(); ++it)
    {
      out << "BC0_id " << it->first << std::endl;
      for (int i=0; i< it->second.size(); ++i)
	{
	  out << "#fpga_" << i << std::endl;
	  const one_FPGA_data_type& afpga=(it->second)[i];
	  for (one_FPGA_data_type::const_iterator itc=afpga.begin(); itc !=afpga.end(); ++itc)
	    out << itc->first <<"," << itc->second << std::endl;
	}
    }
}

bool read(const char* filename, FEB_frame_data_type& FEB_data)
{
  std::ifstream input(filename);
  if (not input.good()) return false;
  char c;
  BC0_id_type BC0id;
  fpga_type fpga;
  channel_type chan;
  TDC_value_type tdcval;
  input >> c;
  while (input.good())
    {
      if (c=='B')
	{
	  //std::cout << " got B" << std::endl;
	  input >> c; if (c != 'C') {std::cout << "error C expected got " << c << std::endl; return false;}  
	  input >> c; if (c != '0') {std::cout << "error 0 expected got " << c << std::endl; return false;}  
	  input >> c; if (c != '_') {std::cout << "error _ expected got " << c << std::endl; return false;}  
	  input >> c; if (c != 'i') {std::cout << "error i expected got " << c << std::endl; return false;}  
	  input >> c; if (c != 'd') {std::cout << "error d expected got " << c << std::endl; return false;}
	  input >> BC0id;
	  std::cout << "got BC0_id = " << BC0id << std::endl;
	}
      else if (c=='#')
	{
	  //std::cout << " got #" << std::endl;
	  input >> c; if (c != 'f') {std::cout << "error f expected got " << c << std::endl; return false;}  
	  input >> c; if (c != 'p') {std::cout << "error p expected got " << c << std::endl; return false;}  
	  input >> c; if (c != 'g') {std::cout << "error g expected got " << c << std::endl; return false;}  
	  input >> c; if (c != 'a') {std::cout << "error a expected got " << c << std::endl; return false;}  
	  input >> c; if (c != '_') {std::cout << "error _ expected got " << c << std::endl; return false;}
	  input >> fpga;
	  //std::cout  << "got fpga = " << fpga << std::endl;
	}
      else if ( (c >= '0') && (c <= '9')) 
	{
	  input.unget();
	  input >> chan;
	  //std::cout << " got channel=" << chan << " ";
	  input >> c; if (c != ',') {std::cout << "error , expected got " << c << std::endl; return false;}
	  input >> tdcval;
	  //std::cout << tdcval << "// ";
	  FEB_data[BC0id][fpga].push_back(TDC_data_type(chan,tdcval));				  
	}
      else
	{
	  std::cout << " got unexpected character " << c << std::endl;
	  return false;
	}
      input >> c;
    }
  for (FEB_frame_data_type::iterator it=FEB_data.begin(); it!=FEB_data.end(); ++it)
    for (three_FPGA_data_type::iterator it2=it->second.begin(); it2 != it->second.end(); ++it2)
      std::sort(it2->begin(),it2->end());
  return true;
}
