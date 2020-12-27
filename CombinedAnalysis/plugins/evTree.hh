#ifndef _evtree_h
#define _evtree_h
struct FullEventTree {
  uint32_t bc;
  uint32_t run;
  uint32_t gtc;
  uint32_t event;
  uint16_t npad;
  uint16_t pad_dif[1024];
  uint8_t pad_asic[1024];
  uint8_t pad_channel[1024];
  uint16_t ntel;
  float tel_x[1024];
  float tel_y[1024];
  float tel_z[1024];
  float tk_x[3];
  float tk_v[3];
  float tk_pchi2;
  uint64_t tk_plans;
  float pex_x[3];
  uint16_t ninti;
  uint8_t f_feb[2048];
  uint8_t f_channel[2048];
  uint64_t f_coarse[2048];
  uint8_t f_fine[2048];
  uint8_t found_feb;
};
#endif
