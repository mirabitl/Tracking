#ifndef _RECOTRACK_HH

#define _RECOTRACK_HH
#include <limits.h>
#include <stdint.h>
#include <vector>
#include <TLine.h>


#include <Math/PositionVector3D.h>
#include <Math/Point3Dfwd.h>
#include <Math/Vector3Dfwd.h>
#include <Math/DisplacementVector3D.h>
#include <ostream>
#include <bitset>



class recoPoint : public ROOT::Math::XYZPoint
{
public:
  recoPoint(double x,double y,double z,uint32_t plan)
  {
    this->SetXYZ(x,y,z);
    valid_=true;_used=false;
    _plan=plan;
    
  }
  recoPoint(uint32_t plan){_plan=plan;}
  ~recoPoint(){;}
  inline double const dX(){return 0.8;}
  inline double const dY(){return 0.8;}
  inline uint32_t const plan(){return _plan;}
  inline bool isUsed() const {return _used;}
  inline void setUse(bool t ){_used=t;}
  bool operator < (const recoPoint& str) const
  {
    return (Z() < str.Z());
  }
  void setValidity(bool t){valid_=t;}
  bool isValid() const {return valid_;}
  void setPlan(uint32_t p){_plan=p;}
protected:
  bool valid_,_used;
  uint32_t _plan;
};


class recoTrack
{
public:
  recoTrack();
  ~recoTrack(){;}
  recoTrack(const recoTrack& t)
  {
    for (uint32_t i=0;i<t.size();i++)
      _points.push_back(t.at(i));
    _orig=t.orig();
    _dir=t.dir();
    _chi2=t.chi2();
    _pchi2=t.pchi2();
    zmin_=t.zmin();
    zmax_=t.zmax();
  }



  TLine* linex();
  TLine* liney();
  void clear();
  //double cap(recoTrack& t);
  void cap(recoTrack& t,double &d,  ROOT::Math::XYZPoint  &p1, ROOT::Math::XYZPoint &p2);
  void clean(float cut);
  void addPoint(recoPoint* p);
  void remove(recoPoint* p);
  double distance(ROOT::Math::XYZPoint* p);
  ROOT::Math::XYZPoint extrapolate(double z);
  void regression();
  inline bool isValid(){return _valid;}
  inline void setValid(bool t){_valid=t;}
  inline uint32_t size() const {return _points.size();}
   ROOT::Math::XYZPoint orig()  const {return _orig;}
  ROOT::Math::XYZVector dir() const {return _dir;}
  void setOrig(float x,float y,float z){_orig.SetXYZ(x,y,z);}
  void setDir(float x,float y,float z){_dir.SetXYZ(x,y,z);}
  recoPoint* at(uint32_t i) const {return _points[i];}
  std::vector<recoPoint*>& points() { return _points;}
  void Dump( std::ostream &os=std::cout );
  void setChi2(double chi2);
  void calculateChi2();
  inline double chi2() const {return _chi2;}
  inline double pchi2() const {return _pchi2;}
  inline float zmin() const { return zmin_;}
  inline float zmax() const { return zmax_;}
  friend std::ostream &operator<<( std::ostream &os, 
                                       const recoTrack &obj)
  {
    // write obj to stream
    os<<" X:"<<obj.orig().X();
    os<<" Y:"<<obj.orig().Y();
    os<<" Z:"<<obj.orig().Z();
    os<<" AX:"<<obj.dir().X();
    os<<" AY:"<<obj.dir().Y()<<std::endl;
    os<<" N points : "<<obj.size()<<std::endl;
    return os;
  }
  bool planUsed(uint32_t p){return _plh[p]==1;}
  uint64_t plans(){return _plh.to_ulong();}
private:
  bool _valid;
  float zmin_,zmax_;
  std::vector<recoPoint*> _points;
  ROOT::Math::XYZVector _dir;
  ROOT::Math::XYZPoint  _orig;
  double _chi2,_pchi2;
  std::bitset<16> _plh;
};
#endif
