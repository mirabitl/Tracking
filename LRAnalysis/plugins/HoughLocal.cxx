#include "HoughLocal.hh"

#include <TLine.h>
HoughLocal::HoughLocal(float thmin,float thmax,float rmin,float rmax,uint32_t nbintheta,uint32_t nbinr) : theThetaMin_(thmin),theThetaMax_(thmax),theRMin_(rmin),theRMax_(rmax),theNbinTheta_(nbintheta),theNbinR_(nbinr)
{
  this->clear();
  theThetaBin_=(theThetaMax_-theThetaMin_)*1./theNbinTheta_;
  theRBin_=(theRMax_-theRMin_)*1./theNbinR_;
	
  for (uint32_t i=0;i<theNbinTheta_;i++)
    {
      double theta=theThetaMin_+(i+0.5)*theThetaBin_;
      theSin_[i]=sin(theta);
      theCos_[i]=cos(theta);

    }
}
HoughLocal::~HoughLocal(){;}
void HoughLocal::fill(std::vector<ROOT::Math::XYZPoint> *vp,uint32_t dir)
{
  thePoints_=vp;
  theVoteMax_=0;
  for (std::vector<ROOT::Math::XYZPoint>::iterator ip=thePoints_->begin();ip!=thePoints_->end();ip++)
    {
      float z=(*ip).Z();
      float x=(*ip).X();
      if (dir!=0)  x=(*ip).Y();
      for (uint32_t i=0;i<theNbinTheta_;i++)
	{
	  float r = z*theCos_[i]+x*theSin_[i];
	  if (r>theRMax_ || r<theRMin_) continue;
	  uint32_t ir=int(round((r-theRMin_)/theRBin_));
	  theHoughImage_[i][ir]+=1;
	  if (theHoughImage_[i][ir]>theVoteMax_) theVoteMax_=theHoughImage_[i][ir];
	}
    }
  //printf("The Vote maximum is %d \n",theVoteMax_);
}
void HoughLocal::clear() 
{
  memset(theHoughImage_,0,512*512*sizeof(uint16_t));

}
float HoughLocal::getTheta(uint32_t i) {return theThetaMin_+i*theThetaBin_;}
float HoughLocal::getR(uint32_t i) {return theRMin_+i*theRBin_;}
uint16_t HoughLocal::getValue(uint32_t i,uint32_t j){return theHoughImage_[i][j];}

void HoughLocal::findMaxima(std::vector< std::pair<uint32_t,uint32_t> >& maxval,uint32_t cut)
{
  if (theVoteMax_ <= cut) return;
  double theVoteBin_=theVoteMax_*1./(theNbinTheta_*theNbinR_);
  //printf("the Vote bin %f \n",theVoteBin_);
  uint16_t theVotes_[theNbinTheta_*theNbinR_];
  memset(theVotes_,0,theNbinTheta_*theNbinR_*sizeof(uint16_t));
  for (uint32_t ith=0;ith<theNbinTheta_;ith++)
    for (uint32_t ir=0;ir<theNbinR_;ir++)
      for (uint32_t iv=0;iv<theNbinTheta_*theNbinR_;iv++) if (theHoughImage_[ith][ir]>iv*theVoteBin_) theVotes_[iv]++;
	
  uint32_t ivmin=0;
  for (uint32_t iv=theNbinTheta_*theNbinR_-1;iv>=0;iv--) if (theVotes_[iv]>cut) {ivmin=iv+1;break;}
  if (ivmin == theNbinTheta_*theNbinR_) ivmin--;
  float vcut=ivmin*theVoteBin_;
  //printf("The vote cut is %f \n",vcut);
  for (uint32_t ith=0;ith<theNbinTheta_;ith++)
    for (uint32_t ir=0;ir<theNbinR_;ir++)
      if (theHoughImage_[ith][ir]>vcut)
	{
	  std::pair <uint32_t,uint32_t> p(ith,ir);
	  maxval.push_back(p);
	}
	
}
static TCanvas* CanvasHough=NULL;
void HoughLocal::draw(DCHistogramHandler* h,std::vector< std::pair<uint32_t,uint32_t> > *maxval)
{
  if (CanvasHough==NULL)
    {
      CanvasHough=new TCanvas("CanvasHough","hough",800,900);
      CanvasHough->Modified();
      CanvasHough->Draw();
      CanvasHough->Divide(1,2);
    }
  CanvasHough->cd();

  TH2F* hhtx = (TH2F*) h->GetTH2("HoughLocal");
  TH2F* hx = (TH2F*) h->GetTH2("LocalImage");
  if (hhtx==NULL)
    {
      hhtx =(TH2F*)h->BookTH2("HoughTransform",theNbinTheta_,theThetaMin_,theThetaMax_,theNbinR_,theRMin_,theRMax_);
      hx =(TH2F*)h->BookTH2("HOriginalX",200,0.,2.8*60.,200,0.,100.);
    }
  else
    {
      hhtx->Reset();
      hx->Reset();
    }
  for (std::vector<ROOT::Math::XYZPoint>::iterator ip=thePoints_->begin();ip!=thePoints_->end();ip++)
    {
      hx->Fill((*ip).Z(),(*ip).X());
		
    }
		
	
  for (uint32_t i=0;i<theNbinTheta_;i++)
    for (uint32_t j=0;j<theNbinR_;j++)
      {
	hhtx->SetBinContent(i+1,j+1,theHoughImage_[i][j]*1.);

		
      }
	
  CanvasHough->cd(1);
  hhtx->Draw("COLZ");
  //hw->Draw();
  CanvasHough->cd(2);
  hx->SetFillColor(3);
  hx->Draw("BOX");
  std::vector<TLine*> vline;vline.clear();
  if (maxval!=NULL)
    {
      for (std::vector < std::pair<uint32_t,uint32_t> >::iterator ihb=maxval->begin();ihb<maxval->end();ihb++)
	{
	  uint32_t ith=(*ihb).first;
	  uint32_t ir=(*ihb).second;
	  float theta = (this->getTheta(ith)+this->getTheta(ith+1))/2;
	  float r = (this->getR(ir)+this->getR(ir+1))/2;
	  float a=-1./tan(theta);
	  float b=r/sin(theta);
	  TLine* l =new TLine(0.,b,2.8*60,a*2.8*60+b);
	  l->SetLineColor(2);
	  l->Draw("SAME");
	  vline.push_back(l);
	}
    }
  CanvasHough->Modified();
  CanvasHough->Draw();
  CanvasHough->Update();

  char c;c=getchar();putchar(c); if (c=='.') exit(0);
  for (std::vector<TLine*>::iterator il=vline.begin();il!=vline.end();il++) delete (*il);

}
