//
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "FEB_v2_reader.cpp"

#include "FEB_v2_data.cpp"
#include "tdcmacro.hh"


#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/dir.h>  
#include <sys/param.h>  
#include <stdio.h>  
#include <stdlib.h>  
#include <unistd.h>  
#include <string.h>
#include <stdint.h>
#include <fcntl.h>
#include <iostream>
#include <sstream>
#include <map>
#include <bitset>
#include "TCanvas.h"
#include <dirent.h>
#include <fnmatch.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>




#define MAXFRAME 65535
typedef struct {
  UInt_t run;
  UInt_t event;
  UInt_t bc0;
  UInt_t    nframe;
  ULong64_t   frame[MAXFRAME];
} bin23Event_t;
TFile* _fdOut=NULL;
TTree* _tree=NULL;
bin23Event_t _gEvent;
void closeTree()
{
  printf("bin2writer::closeTree \n");
  if (_fdOut!=NULL)
    {
      printf("bin2writer Tree header is writen \n");
      if (_tree)
	_tree->Write();
      printf("bin2writer File is closed \n");
      _fdOut->Close();
      //delete _fdOut;
      _fdOut=NULL;
    }
  if (_tree!=NULL)
    {
      //delete _tree;
      _tree=NULL;
    }
}

void createTree(std::string dirp,std::string pref="Run")
{
  std::string _directory=dirp;
  if (_fdOut!=NULL)
    {
      if (_tree)
	_tree->Write();
      _fdOut->Close();
      printf("On detruit fdOut \n");fflush(stdout);
      delete _fdOut;
      _fdOut=NULL;
    }
  if (_tree!=NULL)
    {
      printf("On detruit le ttree \n");fflush(stdout);
      //delete _tree;
      _tree=NULL;
    }
  std::stringstream filename("");    
  char dateStr [64];
  
  time_t tm= time(NULL);
  strftime(dateStr,20,"SMM_%d%m%y_%H%M%S",localtime(&tm));
  filename<<_directory<<"/"<<pref<<"_"<<_gEvent.run<<".root";
  
  _fdOut = new TFile(filename.str().c_str(),"recreate");
  _tree= new TTree("evt","a Tree with SDHCAL frame storage");
  _tree->SetAutoSave(0);
  _tree->Branch("run",&_gEvent.run,"run/i");
  _tree->Branch("event",&_gEvent.event,"event/i");
  _tree->Branch("bc0",&_gEvent.bc0,"bc0/i");
  _tree->Branch("nframe",&_gEvent.nframe,"nframe/i");
  _tree->Branch("frame",_gEvent.frame,"frame[nframe]/l");

}

int main(int argc, char *argv[])
{


  std::string dirp="./rawdata";
  std::string diro="./";
  std::string pats="HV_3_SN_13";
  int32_t runask=0,nmax=5000000;
  char c;
  while ( (c = getopt(argc, argv, "O:r:d:s:m:p:hvN")) != -1 ) {
    fprintf(stderr,"%c read\n",c);
     
    switch ( c ) {
    case 'N':
      /* On stocke le pointeur vers l'argument
       * (car sa valeur pourrait être écrasée
       * au prochain tour de boucle). */
      fprintf(stderr, "Noise run \n");
      break;

    case 'O':
      /* Chaque fois que l'option -v est utilisée,
       * on augmente le degré de verbosité. */
	      
      diro.assign(optarg);
      fprintf(stderr, "Output directory %s \n",optarg);
      break;
    case 'd':
      /* Chaque fois que l'option -v est utilisée,
       * on augmente le degré de verbosité. */
	      
      dirp.assign(optarg);
      fprintf(stderr, "Directory %s \n",optarg);
      break;
    case 's':
      /* Chaque fois que l'option -v est utilisée,
       * on augmente le degré de verbosité. */
	      
      pats.assign(optarg);
      fprintf(stderr, "pattern %s \n",optarg);
      break;
    case 'p':
      /* Chaque fois que l'option -v est utilisée,
       * on augmente le degré de verbosité. */
	      
      //prn.assign(optarg);
      fprintf(stderr, "Processor %s \n",optarg);
      break;


    case 'r':
      /* Ici, on convertit la valeur de optarg
       * en entier. Il conviendrait de gérer les
       * cas où la conversion est impossible (si
       * l'utilisateur a fait n'importe quoi...). */
      runask = atoi(optarg);
      fprintf(stderr, "Run %s \n",optarg);

      break;

    case 'h':
      /* On affiche l'aide et on termine. */
      fprintf(stderr, "Usage: tdcr [options] \n"
	      "Options:\n"
	      "  -v         Augmente la verbosité.\n"
	      "  -V n       Définit le niveau de verbosité à n.\n"
	      "  -h         Affiche ce message d'aide.\n"
	      "  -f FICHIER Écrit dans FICHIER plutôt que sur la sortie standard.\n"
	      );
      return 0;

    case 'v':
      /* On stocke le pointeur vers l'argument
       * (car sa valeur pourrait être écrasée
       * au prochain tour de boucle). */
      fprintf(stderr, "On serait verbose Usage: tdcr [options] \n");
      break;
    case 'm':
      /* Ici, on convertit la valeur de optarg
       * en entier. Il conviendrait de gérer les
       * cas où la conversion est impossible (si
       * l'utilisateur a fait n'importe quoi...). */
      nmax = atoi(optarg);
      fprintf(stderr, "Max evt %s \n",optarg);

      break;
    case '?':
      /* getopt renvoie (par défaut) '?' en cas
       * d'erreur, si une option non acceptée
       * est utilisée ou si une option attendant
       * un argument n'est pas suivi de son argument. */
      fprintf(stderr, "Ligne de commande incorrecte.\n");
      return 1;
    }
  }







  
  //bs.geometry("gifpp_geom.json");
  std::stringstream spat;
  //int runask=atol(argv[1]);
  spat<<"*"<<pats<<"*R"<<runask<<"_SMM*.dat";
  //spat<<"SMM*"<<argv[1]<<"*.dat";
  std::vector<std::string> vfile;
  vfile.clear();
  struct dirent **namelist;
  int n;
  std::cout<<"Pattern "<<spat.str()<<std::endl;
  //std::string dirp="/data/srv02/RAID6/Dome0718";

  //dirp=".";
  n = scandir(dirp.c_str(), &namelist, NULL, alphasort);
  if (n < 0)
    perror("scandir");
  else {
    while (n--) {

      if (fnmatch(spat.str().c_str(), namelist[n]->d_name, 0)==0)
	{
	  printf("%s %d \n", namelist[n]->d_name,fnmatch(spat.str().c_str(), namelist[n]->d_name, 0));
	  printf("found\n");
	  std::stringstream sf;
	  sf<<dirp<<"/"<< namelist[n]->d_name;
	  vfile.push_back(sf.str());
	}
      free(namelist[n]);
    }
    free(namelist);
  }

  





  std::sort(vfile.begin(), vfile.end());

  FEB_v2_data_reader f;
  uint32_t evn=1,shift=0;
  int lastBc0=0,evtCounter=0,lastEvt=0;
  for (auto x:vfile)
    {
      std::cout<<"Opening "<<x<<std::endl;
      //getchar();
      bool OK=f.openDataFile(x.c_str());
      if (not OK) continue;
      int count=0;
      while (f.nextTDCFrame())
	{
	  if (not f.checkGBTFramePadding()) std::cout << "Arrghh" << std::endl;



	  if (f.BC0id()!=lastBc0)
	    {
	  
	      if (lastBc0!=0)
		{
		  //Fill ntuple
		  if (_gEvent.event%1000==1)
		    std::cout<<"Filling "<<_gEvent.run <<"/"<<_gEvent.event<<" BC0 "<<_gEvent.bc0<<" Nframes "<<_gEvent.nframe<<std::endl;


		  if ((f.BC0id()-lastBc0)>1000)
		    {
		      fprintf(stderr,"\t \t ===> IPBUS error %d %d %f \n",f.BC0id(),lastBc0,_gEvent.bc0*1./_gEvent.event);
		      std::flush(std::cout);
		      //getchar();
		    }
		  _tree->Fill();
	       	  evn++;
		}
	      else
		{
		  _gEvent.run=f.runNumber();

		  createTree(diro,pats);
		}
	      if (f.eventNumber()!=lastEvt)
		{
		  evtCounter++;
		  lastEvt=f.eventNumber();
		}
	      _gEvent.run=f.runNumber();
	      _gEvent.event=evtCounter;
	  
	      _gEvent.bc0=f.BC0id();
	      _gEvent.nframe=0;
	      lastBc0=f.BC0id();
	    }
	  //std::cout << "Event R "<<f.eventNumber()<< " BC0 " << f.BC0id() << " fpga " << f.fpga() << " channel " << f.TDC_channel() << " timestamp " << f.TDC_value()  << std::endl;
	  _gEvent.frame[_gEvent.nframe]=m_encode(f.TDC_value(),f.fpga(),f.TDC_channel(),c_side(f.TDC_channel()),c_strip(f.fpga(),f.TDC_channel()),c_petiroc(f.TDC_channel()));
	  //uint64_t m=_gEvent.frame[_gEvent.nframe];
	  //std::cout << "Event T "<<_gEvent.event<< " BC0 " << _gEvent.bc0 << " fpga " << m_fpga(m) << " channel " << m_channel(m) << " timestamp " << m_traw(m)<< std::endl;
      
	  _gEvent.nframe++;
	  ++count;
	}
      shift+=_gEvent.event+1;
      f.close();
      std::cout << "file"<<x<<" contains "  << count << " frame" << std::endl;
      //getchar();
      count=0;
    }

  closeTree();
  /*
    std::ofstream fout("test_out.txt");
    FEB_frame_data_type fd;
    decode (argv[1],fd);
    print(fd,fout);
    if (argc>2)
    {      
    FEB_frame_data_type falban;
    read(argv[2],falban);
    std::ofstream foutb("test_out_check.txt");
    print(falban,foutb);
    }
  */
  return 0;
}
