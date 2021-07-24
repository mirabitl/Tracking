#define FebAna_cxx
#include "FebAna.h"
#include "tdcmacro.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void FebAna::Loop()
{
//   In a ROOT session, you can do:
//      root> .L FebAna.C
//      root> FebAna t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   uint64_t t_bc0[3],trig[3];
   //TH1F* htime=new TH1F("htime","Time du 24",400,15560.,15580.);

   bool debug=false;
   uint32_t lastEvt=0,lastbc0=0,delbc0=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (event%10000==1 || (event-lastEvt)>1)
	{std::cout<<"run "<<run<<"/"<<event<<"/"<<bc0*89E-6<<" Frames "<<nframe<<" entires "<<jentry<<"/"<<nentries<<" last Evt "<<delbc0<<std::flush<<std::endl;
	  if ( (event-lastEvt)>1)
	    getchar();
	}
      if (lastEvt!=event)
	delbc0=bc0-lastbc0;
      lastEvt=event;
      lastbc0=bc0;

      if (nframe>80)
	{
	  fprintf(stderr,"More than 20 frames %d\n",nframe);
	  debug=true;
	}
      else
	debug=false;
      for (uint32_t i=0;i<nframe;i++)
	{
	  if (debug)
	    std::cout<<"\t "<<i<<"/"<<m_fpga(frame[i])<<" "<<m_channel(frame[i])<<" "<<m_traw(frame[i])<<std::endl;
	  if (m_channel(frame[i])==32)
	    t_bc0[m_fpga(frame[i])]=m_traw(frame[i])-9139200;
	  if (m_channel(frame[i])==33)
	    trig[m_fpga(frame[i])]=m_traw(frame[i]);
	}
      //std::cout<<t_bc0[0]<<" "<<t_bc0[1]<<" "<<t_bc0[2]<<std::endl;
      for (uint32_t i=0;i<nframe;i++)
	{
	  if (m_channel(frame[i])==32) continue;
	  //std::cout<<"\t "<<i<<" "<<m_fpga(frame[i])<<" "<<m_channel(frame[i])<<" "<<m_traw(frame[i])<<" corrected "<<m_tcor(m_traw(frame[i]),0)<<std::endl;
	  std::stringstream ss;
	  ss.str("");
	  ss.clear(); 
	  uint32_t cha=m_fpga(frame[i])*48+m_channel(frame[i]);
	  uint32_t ch=m_channel(frame[i]);
	  ss<<"/HitStudy/time"<<cha;
	  TH1* htime= _rh->GetTH1(ss.str());
	  float tc=m_tcor(m_traw(frame[i]),0);
	  if (htime==NULL)
	   {
	    std::cout<<"Booking "<<ss.str()<<std::endl;
	    htime=_rh->BookTH1(ss.str(),5000,tc-50,tc+50);

	  }
	  if (!debug) 
	    htime->Fill(tc);
	  ss.str("");
	  ss.clear(); 
	  ss<<"/HitStudy/Rate/Count"<<m_fpga(frame[i]);
	  TH1* hcount= _rh->GetTH1(ss.str());
	  if (hcount==NULL)
	   {
	    std::cout<<"Booking "<<ss.str()<<std::endl;
	    hcount=_rh->BookTH1(ss.str(),40,-0.1,39.9);
	  }
	  hcount->Fill(ch);
	  ss.str("");
	  ss.clear(); 
	  ss<<"/HitStudy/Rate/Top/Count"<<m_fpga(frame[i]);
	  TH1* hctop= _rh->GetTH1(ss.str());
	  if (hctop==NULL)
	   {
	    std::cout<<"Booking "<<ss.str()<<std::endl;
	    hctop=_rh->BookTH1(ss.str(),40,-0.1,39.9);
	  }
	  if (ch<16 || ch>31)
	    hctop->Fill(ch);
	  ss.str("");
	  ss.clear(); 
	  ss<<"/HitStudy/Rate/Bottom/Count"<<m_fpga(frame[i]);
	  TH1* hcbot= _rh->GetTH1(ss.str());
	  if (hcbot==NULL)
	   {
	    std::cout<<"Booking "<<ss.str()<<std::endl;
	    hcbot=_rh->BookTH1(ss.str(),40,-0.1,39.9);
	  }
	  if (ch>15)
	    hcbot->Fill(ch);
	  
	  
	}
      //if (debug) getchar();
	
      // if (Cut(ientry) < 0) continue;
   }
}
