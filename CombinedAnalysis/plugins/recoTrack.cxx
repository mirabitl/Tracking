#include "recoTrack.hh"
#include "TMath.h"
recoTrack::recoTrack() :_valid(false)
{
  clear();
}




void recoTrack::clear()
{
  _points.clear();
  _orig.SetXYZ(0,0,0);
  _dir.SetXYZ(0,0,0);
  _plh.reset();
}

void recoTrack::addPoint(recoPoint* p)
{
  _points.push_back(p);
  _valid=(_points.size()>=2);

  float z0=999999.,zm=-999999.;
  for (auto ip=_points.begin();ip!=_points.end();ip++)
    {
      if ((*ip)->Z()<z0) z0=(*ip)->Z();
      if ((*ip)->Z()>zm) zm=(*ip)->Z();
    }
  _valid=(zm-z0)>1.;
  _plh.set(p->plan());
  this->regression();
}

void recoTrack::cap(recoTrack& t,double &dist, ROOT::Math::XYZPoint  &p1, ROOT::Math::XYZPoint &p2)
{
  // Closets time approach

  
  ROOT::Math::XYZVector dv = dir() - t.dir();
  float cpatime=0;
  float    dv2 = dv.Mag2();
  if (dv2 > 1E-6)      // the  tracks are almost parallel
    {
      ROOT::Math::XYZVector  w0 = orig() - t.orig();
      cpatime = -1*w0.Dot(dv) / dv2;
    }

  p1=orig()+cpatime*dir();
   p2=t.orig()+cpatime*t.dir();

  ROOT::Math::XYZVector d=p1-p2;

  //std::cout<<p1.X()<<" "<<p1.Y()<<" "<<p1.Z()<<std::endl;
  dist=sqrt(d.Mag2());

  
}
double recoTrack::distance(ROOT::Math::XYZPoint* p)
{
  if (_dir.Mag2()<1E-3)
    return 999999.;
  ROOT::Math::XYZVector ba=(*p)-_orig;
  ROOT::Math::XYZVector cross=_dir.Cross(ba);
  double m2=cross.Mag2()/_dir.Mag2();
  m2=sqrt(m2);
  return m2;
}
ROOT::Math::XYZPoint recoTrack::extrapolate(double z)
{
  ROOT::Math::XYZPoint e;
  e.SetXYZ(_dir.X()*z+_orig.X(),_dir.Y()*z+_orig.Y(),z);
  return e;
}
void recoTrack::regression()
{
  unsigned int n = _points.size();

  if (!_valid) return;
  double zbar=0;
  double xbar=0;
  double ybar=0;
  double z2bar =0;
  double zxbar=0;
  double zybar=0;
  zmin_=1000000;
  zmax_=-100000;
  double wxt=0;
  double wyt=0;
  _plh.reset();
  for (std::vector<recoPoint*>::iterator ip=_points.begin();ip!=_points.end();ip++)
    {
      if ((*ip)->Z()<zmin_) zmin_=(*ip)->Z();
      if ((*ip)->Z()>zmax_) zmax_=(*ip)->Z();
      _plh.set((*ip)->plan());
      zbar+=(*ip)->Z();
      z2bar+=(*ip)->Z()*(*ip)->Z();
      zxbar+=(*ip)->Z()*(*ip)->X();
      zybar+=(*ip)->Z()*(*ip)->Y();
      xbar+=(*ip)->X();
      ybar+=(*ip)->Y();
    }
  zbar /=n;
  z2bar /=n;
  zxbar /=n;
  zybar /=n;
  ybar /=n;
  xbar /=n;
  double s2z = z2bar-zbar*zbar;
  double szx = zxbar-zbar*xbar;
  double szy = zybar-zbar*ybar;
  double ax_ = szx/s2z;double bx_=xbar -ax_*zbar;
  double ay_ = szy/s2z;double by_=ybar -ay_*zbar;
  _orig.SetXYZ(bx_,by_,0);
  _dir.SetXYZ(ax_,ay_,1.0);

 
}
void recoTrack::Dump(std::ostream &os) 
{
    for (std::vector<recoPoint*>::iterator ip=_points.begin();ip!=_points.end();ip++)
    {
      os<<"\t ptr:"<<(*ip)<<" ("<<(*ip)->X()<<":"<<(*ip)->Y()<<":"<<(*ip)->Z()<<")"<<" ==>"<<distance((*ip))<<std::endl;
    }

}

TLine* recoTrack::linex()
{
  return new TLine(zmin_,_orig.X()+_dir.X()*zmin_,zmax_,_orig.X()+_dir.X()*zmax_);
}
TLine* recoTrack::liney()
{
  return new TLine(zmin_,_orig.Y()+_dir.Y()*zmin_,zmax_,_orig.Y()+_dir.Y()*zmax_);
}
void recoTrack::clean(float cut)
{
  for (std::vector<recoPoint*>::iterator ip=_points.begin();ip!=_points.end();)
    {
      if (distance((*ip))>cut)
	_points.erase(ip);
      else
	ip++;
    }
  regression();		

}
void recoTrack::remove(recoPoint* p)
{
  for (std::vector<recoPoint*>::iterator ip=_points.begin();ip!=_points.end();)
    {
      if ((*ip)==p)
	_points.erase(ip);
      else
	ip++;
    }
  regression();		

}
void recoTrack::setChi2(double chi2)
{
  _chi2=chi2;
  _pchi2=TMath::Prob(_chi2,_points.size()*2-4);
}

void recoTrack::calculateChi2()
{
  double chi2=0,paderr=100./96./sqrt(12.);
  for (std::vector<recoPoint*>::iterator ip=this->points().begin();ip!=this->points().end();ip++)
    {
      double cont=this->distance((*ip));
      double err=0;
      bool found=false;
      
      double errx=0.6;
      double erry=0.6;
      err=sqrt(errx*errx+erry*erry);
      chi2+=cont*cont/err/err;
	  //nc++;


    }

  //std::cout<<"chi2 "<<chi2<<" ndf"<<this->points().size()*2-4<<" "<<TMath::Prob(chi2,this->points().size()*2-4)<<std::endl;
  //getchar();
  this->setChi2(chi2);

}

