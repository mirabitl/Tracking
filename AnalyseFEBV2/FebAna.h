//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 21 14:16:02 2021 by ROOT version 6.18/04
// from TTree evt/a Tree with SDHCAL frame storage
// found on file: /data/SMM_210721_134227_1248.root
//////////////////////////////////////////////////////////

#ifndef FebAna_h
#define FebAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "DCHistogramHandler.hh"
#include "tdcmacro.hh"
#include "stdint.h"
#define CALIB_RUN2216
#undef CALIBRATION
// Header file for the classes stored in the TTree if any.
/**
 * 
 **/
class TdcChannel
{
public:
  TdcChannel() : _fr(0), _used(false), _pedMap(NULL) { ; }
  TdcChannel(ULong64_t b, uint8_t feb = 0, std::map<uint64_t, double> *pedmap = NULL) : _fr(b), _used(false), _feb(feb), _0t(0.0), _pedMap(pedmap) { ; }
  inline uint16_t channel() { return m_channel(_fr); }
  inline uint16_t fpga() { return m_fpga(_fr); }
  inline uint16_t pr() { return c_petiroc(channel()); }
  inline uint16_t side() { return c_side(channel()); }
  inline uint16_t strip() { return c_strip(fpga(), channel()); }
  // inline double pedSubTime(jsonFebInfo& f) {return  tdcTime()-f.timeped[channel()];}
  inline double pedestal(uint16_t fp, uint16_t ch)
  {
    if (_pedMap)
      return (*_pedMap)[m_encode(0, fp, ch, c_side(ch), c_strip(fp, ch), c_petiroc(ch))];
    else
      return 0.0;
  }
  inline double pedSubTime() { return tdcTime() - pedestal(fpga(), channel()); }
  inline uint16_t feb() { return _feb; }

  inline uint16_t detectorStrip(uint32_t feb) { return strip() + feb * 48; }

  inline uint8_t length() { return 6; }
  inline uint32_t bcid() { return (uint32_t)(tdcTime() / 200); }
  inline double rawTime() const { return m_traw(_fr); }

  inline double tdcTime() const { return m_tcor(rawTime(), 0); }
  inline uint64_t frame() { return _fr; }
  inline bool used() { return _used; }
  inline void setUsed(bool t) { _used = t; }
  bool operator<(const TdcChannel &ipaddr)
  {
    if (tdcTime() < ipaddr.tdcTime())
      return true;
    else
      return false;
  }
  void dump()
  {
    printf("%d %d %d %f \n", channel(), fpga(), strip(), tdcTime());
  }
  void setZero(double t) { _0t = t; }

private:
  ULong64_t _fr;
  bool _used;
  uint8_t _feb;
  double _0t;
  std::map<uint64_t, double> *_pedMap;
};
#define VEM888 15.58
#define VSTRIP 15.58

class TdcStrip
{
public:
  TdcStrip() : _ch(0), _dif(0), _str(0), _t0(0), _t1(0), _shift(0), _xloc(0), _yloc(0) { ; }
  TdcStrip(uint16_t dif, uint16_t st, double t0, double t1, double shift = -93.48) : _dif(dif), _str(st), _t0(t0), _t1(t1), _shift(shift), _ch(1), _xloc(0), _yloc(0) { getPositions(); }
  TdcStrip(uint16_t ch, uint16_t dif, uint16_t st, double t0, double t1, double shift = -93.48) : _ch(ch), _dif(dif), _str(st), _t0(t0), _t1(t1), _shift(shift), _xloc(0), _yloc(0) { getPositions(); }
  inline uint16_t strip() const { return _str; }
  inline uint16_t chamber() const { return _ch; }
  inline uint16_t dif() const { return _dif; }
  inline double t0() const { return _t0; }
  inline double t1() const { return _t1; }
  inline double TM() const { return (_t1 + _t0) / 2.; }
  inline double shift() const { return _shift; }

  inline void getPositions()
  {
    // Strip in array of XTOF/IVAN
    uint32_t strip_nb = _str; //47-_str;
#define RE41
#undef RE31
#ifdef RE31
    // Xtof to connector (Top of strip to connector)
    double clc[48] = {90.543564, 80.410461, 70.289078, 60.173553, 50.069733, 40.903408, 36.713425, 32.523438, 31.166557, 36.356541, 40.546543, 49.202293, 59.317810, 69.433334, 79.548859, 89.664398, 56.088192, 66.203720, 76.319244, 86.428909, 96.550293, 106.665825, 116.775482, 126.896873, 137.006546, 147.122055, 157.243454, 167.358978, 177.474487, 187.590027, 197.699677, 207.821075, 175.118668, 185.379196, 195.481796, 205.597305, 215.848694, 225.834229, 235.949753, 246.065292, 256.180817, 266.302185, 276.671875, 286.533234, 301.558472, 314.650421, 328.493317, 343.113739};
    // Xtof Return (Low radius strip to connector)
    double clr[48] = {2313.407715, 2303.935059, 2294.468262, 2285.000977, 2275.528320, 2266.067383, 2256.606445, 2247.133789, 2238.666504, 2228.200195, 2218.733398, 2209.266113, 2199.787598, 2190.332520, 2180.877441, 2171.387207, 2120.251953, 2110.785156, 2101.318359, 2091.851074, 2082.384277, 2072.917480, 2063.444824, 2053.983887, 2044.516846, 2035.049927, 2025.583008, 2016.116333, 2006.649414, 1997.182495, 1987.721558, 1978.243042, 1926.101929, 1916.629150, 1907.162231, 1897.701416, 1888.234497, 1878.773438, 1869.300659, 1855.762939, 1844.198853, 1832.120117, 1819.977539, 1807.770386, 1795.510132, 1783.159668, 1770.754272, 1757.014038};
    // Xtof strip length
    double cls[48] = {1468.002441, 1468.022339, 1468.062012, 1468.121460, 1468.200928, 1468.300049, 1468.419067, 1468.557861, 1467.717041, 1468.895020, 1469.093262, 1469.311279, 1469.549194, 1469.806763, 1470.084229, 1470.381348, 1470.698242, 1471.034790, 1471.391235, 1471.767212, 1472.162964, 1472.578247, 1473.013306, 1473.468018, 1473.942383, 1474.436279, 1474.949707, 1475.482788, 1476.035400, 1476.607544, 1477.199097, 1477.810181, 1477.132324, 1479.090576, 1479.760010, 1480.448730, 1480.165649, 1481.884155, 1482.630737, 1481.655518, 1477.538940, 1472.481323, 1466.406738, 1462.278442, 1452.549316, 1436.809814, 1420.182617, 1403.590942};
#endif
#ifdef RE41
 double clc[48]={
		 90.894295,80.668900,70.553459,60.438061,50.316757,40.429695,36.239754,0.890015,
		 31.640133,34.830097,39.020042,44.321323,57.034809,68.150246,77.265640,79.041321,
		 56.276508,64.931931,75.047348,85.162766,95.278206,105.393593,115.509033,125.624458,135.739899,145.855316,155.970718,166.086151,143.920959,187.317017,141.053879,207.547836,0.880005,184.098709,194.214096,204.329529,214.444931,224.560394,234.675842,244.791199,10.737075,265.022095,275.137512,285.252930,300.251160,313.273132,327.085388,340.790131};
double clr[48]={2105.610596,2095.231201,2084.851562,2074.472168,2064.092773,2053.713379,2043.333496,2064.113770,
		2021.574585,2012.194946,2001.815796,1994.034180,1981.056641,1969.677124,1960.297363,1958.257812,
		1895.858643,1886.479004,1876.099609,1865.719971,1855.340576,1844.960938,1834.581787,1824.202026,1813.822510,1803.442993,1793.063354,1782.683838,1804.585083,1760.925049,1807.088135,1740.165894,1861.209717,1677.726807,1667.347534,1656.968018,1646.588623,1636.208740,1625.829346,1615.449707,1845.542236,1588.633545,1575.151978,1561.595825,1547.964844,1534.258301,1520.645142,1505.342773};
double cls[48]={1227.423828,1227.126099,1220.051514,1220.100952,1220.166870,1220.249268,1220.129150,1220.048218,
		1220.830566,1220.743408,1220.908203,1220.385437,1220.453369,1221.862793,1221.177856,1221.770844,1221.770844,1222.520996,1222.817017,1223.129395,1223.458008,1223.803101,1224.164551,1224.542236,1224.936279,1225.346436,1225.773071,1226.215820,1225.540161,1227.914307,1226.424438,1228.963013,1228.963013,1229.212769,1229.768799,1230.340942,1230.929077,1231.533203,1232.153442,1232.789673,1234.478271,1225.268066,1220.836914,1214.699097,1203.984131,1187.263916,1179.133911,1151.987427};
#endif 

    // Calculate position inside the strip
    /* 
       t1 is T High radius
       t0 is T Low radius
       
       lc: distance High radius of strip to connector
       ls : strip length
       lr : distance Low radius of strip to connector
       z : distance along the strip from High radius

       t1= (lc+z)/V
       t0= (ls-z+lr)/V

       => V(t1-t0)=lc+2z -ls -lr

       => z= (lr+ls-lc+V(t1-t0))/2

       Passage mm->cm 
    */

    _Lr = clr[strip_nb] / 10.;
    _Ls = cls[strip_nb] / 10.;
    _Lc = clc[strip_nb] / 10.;

    _z_in_strip = (_Ls + _Lr - _Lc + VEM888 * (_t1 - _t0)) / 2.0;

    //Strips poligon
    /*
PT1     PT0
	+++++
	 \   \
	  \   \
	   \   \
     +++++
	 PB1     PB0

        WE calculate:

         - P Top ( (xt1+xt0)/2,(yt1+yt0)/2)
         - P Bot  ( (xb1+xb0)/2,(yb1+yb0)/2)

         Then Cos(theta) = (Xbot-Xtop)/Lstrip  Sin(theta)= (Ybot-Ytop)/Lstrip

        and 
          X = Xtop+cos(theta)*z
	  Y = Ytop+sin(theta)*z
      */
#ifdef RE31

    double XT0[48] = {0.75, 12.135525, 23.521051, 34.906576, 46.292101, 57.677626, 69.063152, 80.448677, 91.834202, 103.219728, 114.605253, 125.990778, 137.376303, 148.761829, 160.147354, 171.532879, 182.918405, 194.30393, 205.689455, 217.07498, 228.460506, 239.846031, 251.231556, 262.617082, 274.002607, 285.388132, 296.773657, 308.159183, 319.544708, 330.930233, 342.315758, 353.701284, 365.086809, 376.472334, 387.85786, 399.243385, 410.62891, 422.014435, 433.399961, 444.785486, 456.171011, 467.556537, 478.942062, 490.327587, 501.713112, 511.520357, 520.962511, 530.335895};
    double YT0[48] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.54039499999999, 20.82492, 32.027257};

    double XB0[48] = {0.75, 6.74246299999999, 12.734926, 18.727389, 24.719852, 30.712315, 36.704778, 42.697241, 48.689704, 54.682167, 60.67463, 66.667093, 72.659557, 78.65202, 84.644483, 90.636946, 96.629409, 102.621872, 108.614335, 114.606798, 120.599261, 126.591724, 132.584187, 138.57665, 144.569113, 150.561576, 156.554039, 162.546502, 168.538965, 174.531428, 180.523891, 186.516354, 192.508817, 198.50128, 204.493744, 210.486207, 216.47867, 222.471133, 228.463596, 234.456059, 241.025566, 247.902391, 254.827969, 261.802821, 268.827475, 275.902467, 283.028338, 290.20564};
    double YB0[48] = {1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1463.075864, 1457.304243, 1451.491705, 1445.637811, 1439.742119, 1433.80418, 1427.823538, 1421.799731};

    double XB1[48] = {6.24246299999999, 12.234926, 18.227389, 24.219852, 30.212315, 36.204778, 42.197241, 48.189704, 54.182167, 60.17463, 66.167093, 72.159557, 78.15202, 84.144483, 90.136946, 96.129409, 102.121872, 108.114335, 114.106798, 120.099261, 126.091724, 132.084187, 138.07665, 144.069113, 150.061576, 156.054039, 162.046502, 168.038965, 174.031428, 180.023891, 186.016354, 192.008817, 198.00128, 203.993744, 209.986207, 215.97867, 221.971133, 227.963596, 233.956059, 236.35, 247.329979, 254.253528, 261.226337, 268.248933, 275.321852, 282.445635, 289.620835, 296.848009};

    double YB1[48] = {1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1467, 1457.78466, 1451.973824, 1446.121646, 1440.227681, 1434.291482, 1428.312592, 1422.29055, 1416.224886};

    double XT1[48] = {11.635525, 23.021051, 34.406576, 45.792101, 57.177626, 68.563152, 79.948677, 91.334202, 102.719728, 114.105253, 125.490778, 136.876303, 148.261829, 159.647354, 171.032879, 182.418405, 193.80393, 205.189455, 216.57498, 227.960506, 239.346031, 250.731556, 262.117082, 273.502607, 284.888132, 296.273657, 307.659183, 319.044708, 330.430233, 341.815758, 353.201284, 364.586809, 375.972334, 387.35786, 398.743385, 410.12891, 421.514435, 432.899961, 444.285486, 455.671011, 467.056537, 478.442062, 489.827587, 501.213112, 503.537578, 520.546573, 529.921472, 539.22834};
    double YT1[48] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20.327824, 31.531972, 42.654814};
#endif
#ifdef RE41

double XT0[48]={ 0.00,12.00,23.00,34.00,46.00,57.00,69.00,80.00,91.00,103.00,114.00,125.00,137.00,148.00,160.00,171.00,182.00,194.00,205.00,217.00,228.00,239.00,251.00,262.00,274.00,285.00,296.00,308.00,319.00,330.00,342.00,353.00,365.00,376.00,387.00,399.00,410.00,422.00,433.00,444.00,456.00,467.00,478.00,490.00,501.00,511.00,520.00,530.00};
double YT0[48]={ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,10.00,21.00,32.00};
double XB0[48]={ 0.00, 7.00,14.00,21.00,28.00,35.00,42.00,49.00,55.00,62.00,69.00,76.00,83.00,90.00,97.00,104.00,111.00,118.00,125.00,131.00,138.00,145.00,152.00,159.00,166.00,173.00,180.00,187.00,194.00,200.00,207.00,214.00,221.00,228.00,235.00,242.00,249.00,256.00,263.00,270.00,276.00,284.00,292.00,300.00,308.00,316.00,324.00,333.00};
double YB0[48]={1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1216.00,1209.00,1202.00,1195.00,1189.00,1182.00,1175.00};
double XT1[48]={11.00,23.00,34.00,45.00,57.00,68.00,79.00,91.00,102.00,114.00,125.00,136.00,148.00,159.00,171.00,182.00,193.00,205.00,216.00,227.00,239.00,250.00,262.00,273.00,284.00,296.00,307.00,319.00,330.00,341.00,353.00,364.00,375.00,387.00,398.00,410.00,421.00,432.00,444.00,455.00,467.00,478.00,489.00,501.00,503.00,520.00,529.00,539.00};
double YT1[48]={ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,21.00,32.00,43.00};
double XB1[48]={ 7.00,14.00,20.00,27.00,34.00,41.00,48.00,55.00,62.00,69.00,76.00,83.00,90.00,96.00,103.00,110.00,117.00,124.00,131.00,138.00,145.00,152.00,159.00,165.00,172.00,179.00,186.00,193.00,200.00,207.00,214.00,221.00,228.00,235.00,241.00,248.00,255.00,262.00,269.00,276.00,279.00,291.00,299.00,307.00,316.00,324.00,332.00,340.00};
double YB1[48]={1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1220.00,1209.00,1203.00,1196.00,1189.00,1182.00,1175.00,1168.00};



    

#endif
    double xt = (XT0[strip_nb] + XT1[strip_nb]) / 20.;
    double yt = (YT0[strip_nb] + YT1[strip_nb]) / 20.;
    double xb = (XB0[strip_nb] + XB1[strip_nb]) / 20.;
    double yb = (YB0[strip_nb] + YB1[strip_nb]) / 20.;

    double lp = sqrt((xb - xt) * (xb - xt) + (yb - yt) * (yb - yt));
    double ct = (xb - xt) / lp;
    double st = (yb - yt) / lp;

    _xloc = xt + _z_in_strip * ct;
    _yloc = yt + _z_in_strip * st;
  }

  inline void calculate_re31_left()
  {
    // Ivan Strip
    double ils[48] = {1467.009913, 1467.03965, 1467.089211, 1467.158594, 1467.247795, 1467.356812, 1467.356812, 1467.48564, 1467.802739, 1467.802739, 1467.991004, 1468.426773, 1468.674331, 1468.674379, 1469.22871, 1469.535567, 1469.862091, 1470.208325, 1470.574255, 1470.959865, 1471.365142, 1471.790067, 1472.234709, 1472.698886, 1473.182658, 1473.686007, 1474.208912, 1474.751353, 1475.313309, 1475.894866, 1476.495786, 1477.116151, 1477.755937, 1478.415119, 1479.093798, 1479.791695, 1480.508908, 1481.245409, 1482.001169, 1480.6476825, 1476.5305035, 1471.4726985, 1466.38505, 1461.268888, 1451.541249, 1435.797934, 1419.1703525, 1402.5796075};
    // Ivan to connector
    double ilc[48] = {89.234025, 79.118885, 69.003025, 58.8744665, 48.772025, 38.8797705, 34.689425, 30.527362, 29.090515, 33.27841, 37.4704695, 45.3659465, 55.483534, 65.5993955, 75.714525, 85.830384, 53.268025, 63.383882, 73.4332645, 83.6148815, 93.730025, 103.846025, 113.961525, 124.077025, 134.192525, 144.308025, 154.423525, 164.539025, 174.6548805, 184.770025, 194.8858805, 205.001025, 172.438939, 182.557883, 192.66788, 202.783025, 212.89888, 223.014025, 233.12988, 243.245025, 253.361025, 263.482383, 273.592025, 283.707525, 297.651346, 311.827659, 325.665812, 339.3975755};
    // Ivan return line
    double ilr[48] = {2311.5805375, 2302.115654, 2292.6482265, 2283.181437, 2273.714415, 2264.247734, 2257.3496125, 2245.3141605, 2235.846714, 2226.3802165, 2216.9127225, 2207.446365, 2197.978995, 2188.512669, 2179.0455665, 2169.573028, 2117.4320435, 2107.9653625, 2098.4981145, 2089.0313245, 2079.564658, 2070.0974425, 2060.630839, 2051.1578355, 2041.697214, 2032.2301775, 2022.7633945, 2013.2954325, 2003.829547, 1994.356909, 1984.895486, 1975.4288055, 1923.281931, 1913.8154285, 1904.348002, 1894.8815075, 1885.4141275, 1875.9478055, 1866.480487, 1851.9954605, 1841.383921, 1829.306444, 1817.1639285, 1804.956784, 1792.68389, 1780.345062, 1767.862138, 1754.195906};
    // Ivan total
    double ilt[48] = {3867.8244755, 3848.274189, 3828.7404625, 3809.2144975, 3789.734235, 3770.4843165, 3759.3958495, 3743.3271625, 3732.739968, 3727.4613655, 3722.374196, 3721.2390845, 3722.13686, 3722.7864435, 3723.9888015, 3724.938979, 3640.5621595, 3641.5575695, 3642.505634, 3643.606071, 3644.659825, 3645.7335345, 3646.827073, 3647.9337465, 3649.072397, 3650.2242095, 3651.3958315, 3652.5858105, 3653.7977365, 3655.0218, 3656.2771525, 3657.5459815, 3573.476807, 3574.7884305, 3576.10968, 3577.4562275, 3578.8219155, 3580.2072395, 3581.611536, 3575.888168, 3571.2754495, 3564.2615255, 3557.1410035, 3549.933197, 3541.876485, 3527.970655, 3512.6983025, 3496.173089};
    //_z_in_strip=VSTRIP*(_t1-_t0)/2+VSTRIP/VEM888/2.*(_Lr-_Lc)+_Ls/2.;
  }
  //inline double ypos() const {return (_t0-_t1-_shift)/1.;}
  // Avant inline double ypos() const {return _shift+80.0+(160.-(_t1-_t0)*18.39)/2.0;}
  inline double ypos() const
  {
    return _yloc;
  }
  inline double xpos() const
  {
    return _xloc;
  }
  inline double X() { return xpos(); }
  inline double Y() { return ypos(); }
  inline double Lr() { return _Lr; }
  inline double Ls() { return _Ls; }
  inline double Lc() { return _Lc; }
  inline double Zs() { return _z_in_strip; }

private:
  uint16_t _dif, _str, _ch;
  double _t0, _t1, _shift, _xloc, _yloc, _Lr, _Ls, _Lc, _z_in_strip;
};
#define MAX_FRAMES 65535
class FebAna
{
public:
  TTree *fChain;  //!pointer to the analyzed TTree or TChain
  Int_t fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  UInt_t run;
  UInt_t event;
  UInt_t bc0;
  UInt_t nframe;
  ULong64_t frame[MAX_FRAMES]; //[nframe]

  // List of branches
  TBranch *b_run;    //!
  TBranch *b_event;  //!
  TBranch *b_bc0;    //!
  TBranch *b_nframe; //!
  TBranch *b_frame;  //!

  FebAna(TTree *tree = 0) { ; }
  FebAna(std::string fn);
  virtual ~FebAna();
  virtual Int_t Cut(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);
  virtual void Loop();
  virtual Bool_t Notify();
  virtual void Show(Long64_t entry = -1);
  DCHistogramHandler *_rh;
  TH1 *GetTH1(std::string name) { return _rh->GetTH1(name); }
  TH2 *GetTH2(std::string name) { return _rh->GetTH2(name); }
  void writeHistograms(std::string n) { _rh->writeHistograms(n); }
  std::map<uint64_t, double> *pedMap() { return &_pedMap; }
  std::vector<TdcChannel *> &channels() { return _channels; }

  std::vector<TdcChannel *> _channels;
  std::map<uint64_t, double> _pedMap;
};

#endif

#ifdef FebAna_cxx

FebAna::FebAna(std::string fn) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  TTree *tree = NULL;
  if (tree == 0)
  {
    TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(fn.c_str());
    if (!f || !f->IsOpen())
    {
      f = new TFile(fn.c_str());
    }

    f->GetObject("evt", tree);
  }

  Init(tree);

  // Create the histogram handler
  _rh = DCHistogramHandler::instance();

  // Store the time pedestal for the 3 FPGAs


#ifdef CALIB_RUN1216
  double dt0[32] = {0.00, -0.37, 1.06, 5.59, 0.59, 3.83, 0.15, 2.84, 2.19, 0.95, 2.52, 3.89, 9.14, 2.62, 4.05, 7.11, 3.48, 8.20, 7.89, 2.97, 7.50, 4.08, 4.79, 3.10, 4.83, 2.84, 8.87, 7.16, -0.01, 1.85, 3.49, 0.57};
  double dt1[32] = {0.00, -0.47, 1.11, 5.64, 0.57, 3.80, 0.13, 2.86, 2.23, 0.82, 2.61, 4.06, 9.23, 2.59, 4.10, 7.06, 3.32, 8.13, 8.62, 2.60, 7.67, 3.87, 4.64, 3.05, 4.64, 2.58, 8.72, 7.11, -0.15, 1.62, 3.45, 0.49};
  double dt2[32] = {0.00, -0.53, 1.22, 5.74, 0.67, 4.06, 0.19, 2.88, 2.30, 1.08, 2.63, 4.15, 9.32, 2.78, 4.32, 7.33, 3.81, 8.44, 8.03, 3.14, 7.87, 4.32, 5.20, 3.42, 4.99, 3.03, 9.22, 7.59, 0.19, 2.05, 3.76, 0.91};
  double d0 = 12.39, d1 = 12.24, d2 = 12.5;
  for (int i = 0; i < 32; i++)
  {
    _pedMap.insert(std::pair<uint64_t, double>(m_encode(0, 0, i, c_side(i), c_strip(0, i), c_petiroc(i)), dt0[i] + d0));
    _pedMap.insert(std::pair<uint64_t, double>(m_encode(0, 1, i, c_side(i), c_strip(1, i), c_petiroc(i)), dt1[i] + d1));
    _pedMap.insert(std::pair<uint64_t, double>(m_encode(0, 2, i, c_side(i), c_strip(2, i), c_petiroc(i)), dt2[i] + d2));
  }
#endif
#ifdef CALIB_RUN2216 
  double dt0[32]={ 0.00,-0.35, 1.03, 5.61, 0.71, 4.03, 0.09, 2.86, 2.27, 0.95, 2.58, 4.06, 9.37, 2.64, 4.23, 7.06, 3.70, 8.31, 7.72, 2.81, 7.77, 3.96, 5.01, 3.22, 4.72, 2.76, 9.02, 7.25, 0.13, 1.74, 3.43, 0.41};
double dt1[32]={ 0.00,-0.47, 0.93, 5.81, 0.62, 3.83,-0.02, 2.82, 2.16, 0.89, 2.40, 3.88, 9.21, 2.48, 3.97, 6.83, 3.43, 8.30, 7.76, 2.75, 7.72, 4.02, 4.72, 3.04, 4.70, 2.58, 8.90, 7.24,-0.06, 1.66, 3.26, 0.10};
double dt2[32]={ 0.00,-0.24, 1.18, 5.72, 0.59, 4.02, 0.11, 2.87, 2.28, 1.18, 2.71, 4.07, 9.40, 2.76, 4.25, 7.16, 3.60, 8.42, 7.84, 2.89, 7.69, 4.00, 4.92, 3.24, 4.65, 2.88, 8.81, 7.33, 0.12, 1.75, 3.27, 0.42};
  double d0 = 12.44, d1 = 12.50, d2 = 12.54;
  for (int i = 0; i < 32; i++)
  {
    _pedMap.insert(std::pair<uint64_t, double>(m_encode(0, 0, i, c_side(i), c_strip(0, i), c_petiroc(i)), dt0[i] + d0));
    _pedMap.insert(std::pair<uint64_t, double>(m_encode(0, 1, i, c_side(i), c_strip(1, i), c_petiroc(i)), dt1[i] + d1));
    _pedMap.insert(std::pair<uint64_t, double>(m_encode(0, 2, i, c_side(i), c_strip(2, i), c_petiroc(i)), dt2[i] + d2));
  }
#endif
}

FebAna::~FebAna()
{
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
}

Int_t FebAna::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain)
    return 0;
  return fChain->GetEntry(entry);
}
Long64_t FebAna::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (fChain->GetTreeNumber() != fCurrent)
  {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void FebAna::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree)
    return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("bc0", &bc0, &b_bc0);
  fChain->SetBranchAddress("nframe", &nframe, &b_nframe);
  fChain->SetBranchAddress("frame", frame, &b_frame);
  Notify();
}

Bool_t FebAna::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void FebAna::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain)
    return;
  fChain->Show(entry);
}
Int_t FebAna::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef FebAna_cxx
