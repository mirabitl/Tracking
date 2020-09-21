#include <TROOT.h>
#include <TCanvas.h>
#include <TH3.h>
#include <TF1.h>
#include <TList.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>

double transparency(double* x, double* param)
{
	//This function is required to suppress
	//boxes for empty bins - make them transparent.
	if(x)
	{
		if ( *x <= 0 )
			return 0.0 ;
		else if ( *x > param[0] )
			return 0.0 ;
		else if ( *x <= param[0] )
			return 1.0 - 0.95*TMath::TanH( 1.7*(param[0] - *x) ) ;
		else
			return 1.0 ;
	}

	return 0.0 ;
}



void D3V()
{


  
      //  gSystem->Unlink("testAnim.gif") ; // delete old file

	gStyle->SetCanvasPreferGL(kTRUE) ;
	gStyle->SetPalette(1);


	TFile *_file0 = TFile::Open("track1111B.root");

	TCanvas* c1 = new TCanvas("c1","c1",800,200,800,800) ;

    TH3D* const histo = new TH3D ("hhxyz","hhxyz",35,35,35,35,-35,35,30,130,210)
 ;
    TNtuple *Smiths= (TNtuple*)_file0->Get("Smiths");
    
    Smiths->Draw("zc:xc:yc>>histo","zc>140&&zc<230&&costh<.9995&&distdca<5","");

    	histo->Draw("colz") ;
	c1->Modified();
	c1->Update();
	
	getchar();
	TList * const lof = histo->GetListOfFunctions() ;
	if (!lof)
	{
		std::cout << "List of functions is null" << std::endl ;
		delete histo ;
		return ;
	}

	TF1* func = new TF1("TransferFunction", transparency, 0.0, 100.0, 1) ;
	lof->Add( func ) ;

	for ( int i = 0 ; i < 10 ; i++ )
    {
    	double time = (i+1)/10.0 ;

    	  std::cout << "time : " << time << std::endl ;
		func->SetParameter(0,time) ;
		histo->Draw("colz") ;
		// histo->Draw("glcol") ;
		   getchar();
		
		   //		c1->Modified();
		   // c1->Update() ;

		  
    }

    
}
