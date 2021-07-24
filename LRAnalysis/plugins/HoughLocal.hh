#ifndef _HOUGHLOCAL_H

#define _HOUGHLOCAL_H
#include <vector>
#include<stdint.h>
#include "DCHistogramHandler.hh"
#include <Math/PositionVector3D.h>
#include <Math/Point3Dfwd.h>
#include <Math/Vector3Dfwd.h>
#include <Math/DisplacementVector3D.h>
class HoughLocal
{
public:
	HoughLocal(float thmin,float thmax,float rmin,float rmax,uint32_t nbintheta=8,uint32_t nbinr=8);
	~HoughLocal();
	void fill(std::vector<ROOT::Math::XYZPoint> *vp,uint32_t dir=0);
	void clear();
	float getTheta(uint32_t i);
	float getR(uint32_t i);
	uint16_t getValue(uint32_t i,uint32_t j);
	void findMaxima(std::vector< std::pair<uint32_t,uint32_t> >& maxval,uint32_t cut=5);

	void draw(DCHistogramHandler* h,std::vector< std::pair<uint32_t,uint32_t> > *maxval=NULL);
private:
	std::vector<ROOT::Math::XYZPoint> *thePoints_;
	double theSin_[512];
	double theCos_[512];

	uint16_t theHoughImage_[512][512];
	float theThetaMin_;
	float theThetaMax_;
	float theRMin_,theRMax_;
	float theThetaBin_,theRBin_;
	uint32_t theNbinTheta_;
	uint32_t theNbinR_;
	uint16_t theVoteMax_;
};



#endif
