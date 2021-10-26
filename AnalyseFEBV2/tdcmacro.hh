#pragma once
#define m_encode(t,fpga,ch,side,strip,prc) ((uint64_t) t +( (uint64_t) (fpga&0x3)<<32)+( (uint64_t) (ch&0x3F)<<34)+( (uint64_t) (side&0x1)<<40)+( (uint64_t) (strip&0x3F)<<41)+( (uint64_t) (prc&0x1F)<<47))

#define m_traw(m) (m&0XFFFFFFFF)
#define m_fpga(m) ((m>>32)&0X3)
#define m_channel(m) ((m>>34)&0X3F)
#define m_side(m) ((m>>40)&0X1)
#define m_strip(m) ((m>>41)&0X3F)
#define m_pr_channel(m) ((m>>47)&0X1F)
#define m_tcor(t,tbc0) (t-tbc0)*2.5/256

// Side 0 LR 1 HR
#define c_side(ch) ((ch>15)?ch%2:(ch+1)%2)
#define c_local_strip(ch) ((ch>15)?((ch-16)/2):16-(ch/2+1))
#define c_strip(fpga,ch) (fpga*16+c_local_strip(ch))
#define c_petiroc(ch) ((ch>15)?(ch-16)*2:ch*2)

