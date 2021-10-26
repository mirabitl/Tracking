#define FebAna_cxx
#include "FebAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <bitset>
#define TTRIGBEAMDUMP -909
#define TTRIG_NOSOURCE -1388
#define TTRIG -1374
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
	if (fChain == 0)
		return;
	Long64_t nentries = fChain->GetEntriesFast();
	// Limit the analysis to 50000  events
	if (nentries > 450000)
		nentries = 450000;
	Long64_t nbytes = 0, nb = 0;

	uint64_t t_bc0[3], ntrig[3];
	double trig[3];

	bool debug = false;
	bool info=false;
	uint32_t lastEvt = 0, lastbc0 = 100000000, delbc0 = 0, max_strip_count = 0, max_strip10_count;
	bool newevt = false;
	bool dumpinfo = false;

	//
	//Loop on entry
	// Each entry contains data of an orbit
	//
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		//
		// Clear channel array
		//
		for (auto x = channels().begin(); x != channels().end(); x++)
			delete (*x);
		_channels.clear();
		// Trigger channel array and counter
		TdcChannel *trig_ch[3];
		for (int i = 0; i < 3; i++)
		{
			trig_ch[i] = 0;
			ntrig[i] = 0;
		}
		// Load the netry
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		//
		// Dump Evnts bc0 number of frame evry 10000 entries
		//
		if (event % 10000 == 1 || (event - lastEvt) > 1)
		{
			std::cout << "run " << run << "/" << event << "/" << bc0 * 89E-6 << " Frames " << nframe << " entires " << jentry << "/" << nentries << " last Evt " << delbc0 << std::flush << std::endl;
			if ((event - lastEvt) > 1)
				getchar();
		}
		// Just to know the time to read a full event (norbit)
		if (lastEvt != event)
		{
			newevt = true;
			delbc0 = bc0 - lastbc0;
			if (dumpinfo)
				getchar();
		}
		else
			newevt = false;
		//std::cout << event << " " << lastEvt << " " << newevt << std::endl;
		lastEvt = event;
		lastbc0 = bc0;
		// Debug if too many frames
		if (nframe > MAX_FRAMES/2)
		{
			fprintf(stderr, "More than %d frames %d\n", MAX_FRAMES/2,nframe);
			debug = true;
		}
		else
			debug = false;
		// Initialise strips channels vector
		std::vector<TdcChannel *> c_strip[48];
		bool badevt = false;

		// Loop on frames and find trigger time if any in this orbit
		for (uint32_t i = 0; i < nframe; i++)
		{
			// Store channels in the vector
			if (m_channel(frame[i]) < 32)
				_channels.push_back(new TdcChannel(frame[i], 0, pedMap()));
			if (debug)
				std::cout << "\t " << i << "/" << m_fpga(frame[i]) << " " << m_channel(frame[i]) << " " << m_traw(frame[i]) << std::endl;

			// Store BC0 time
			if (m_channel(frame[i]) == 32)
				t_bc0[m_fpga(frame[i])] = m_traw(frame[i]) - 9139200;

			// Trigger channel
			if (m_channel(frame[i]) == 33)
			{
				// Bad event if more than 1 trigger
				if (trig_ch[m_fpga(frame[i])] != 0)
				{
					badevt = true;
				}
				// Store the trigger time, the trigger channel and count 
				trig[m_fpga(frame[i])] = m_tcor(m_traw(frame[i]), 0); //m_traw(frame[i]);
				trig_ch[m_fpga(frame[i])] = new TdcChannel(frame[i]);
				ntrig[m_fpga(frame[i])]++;
			}
		}
		
		/// Histogram of trigger time and event counters
		TH1 *htrig = _rh->GetTH1("Time_trigger");
		TH1 *hcounter = _rh->GetTH1("Counter");
		if (htrig == NULL)
		{

			htrig = _rh->BookTH1("Time_trigger", 90000, 0., 90000.);
			hcounter = _rh->BookTH1("Counter", 50, 0.1, 50.1);
		}
		for (int it = 0; it < 3; it++)
			if (trig_ch[it] != NULL)
			{
				htrig->Fill(trig[it]);
				// Remove trigger 2 micros @ beginning and end of orbits
				#ifndef CALIBRATION
				badevt = badevt || (trig[it] < 2000) || (trig[it] > 87000);
				#endif
				//printf("%d  %f\n",
			}
		// Only One trigger in the orbit
		badevt = badevt || !(ntrig[0] == 1 && ntrig[1] == 1 && ntrig[2] == 1);
		if (info)
		  std::cout << bc0 << " " << ntrig[0] << " " << ntrig[1] << " " << ntrig[2] << std::endl;
		if (debug) badevt=false;
		// Make analysis only for good orbits
		if (!badevt)
		{
			// Good orbits count
			hcounter->Fill(1.);
			// Time calibration histogram
			TProfile *hcalib = (TProfile *)_rh->GetTH1("/Calibration/Distance2Trigger");
			if (hcalib == NULL)
				hcalib = _rh->BookProfile("/Calibration/Distance2Trigger", 400, -1.1, 98.9, -50., 50.);

			// Channel hit study
			for (uint32_t i = 0; i < nframe; i++)
			{
				if (m_channel(frame[i]) == 32)
					continue;
				//std::cout<<"\t "<<i<<" "<<m_fpga(frame[i])<<" "<<m_channel(frame[i])<<" "<<m_traw(frame[i])<<" corrected "<<m_tcor(m_traw(frame[i]),0)<<std::endl;
				std::stringstream ss;
				ss.str("");
				ss.clear();
				uint32_t cha = m_fpga(frame[i]) * 48 + m_channel(frame[i]);
				uint32_t ch = m_channel(frame[i]);
				uint32_t strip = m_strip(frame[i]);
				uint32_t side = m_side(frame[i]);
				ss << "/HitStudy/chan" << cha;
				TH1 *htime = _rh->GetTH1(ss.str() + "/rawtime");
				TH1 *httrig = _rh->GetTH1(ss.str() + "/dist_trigger");
				TH1 *hctrig = _rh->GetTH1(ss.str() + "/cal_trigger");
				float tc = m_tcor(m_traw(frame[i]), 0);

				if (htime == NULL)
				{
					//std::cout << "Booking " << ss.str() << std::endl;
					htime = _rh->BookTH1(ss.str() + "/rawtime", 5000, tc - 50, tc + 50);
					httrig = _rh->BookTH1(ss.str() + "/dist_trigger", 3000, -1500, 1500);
					hctrig = _rh->BookTH1(ss.str() + "/cal_trigger", 2500, -10, 40);
				}

				// Selected events are +- 40 ns around TTRIG

				float tmin = TTRIG - 40;
				float tmax = TTRIG + 40;
				bool selected = (tc - trig[m_fpga(frame[i])]) > tmin && (tc - trig[m_fpga(frame[i])]) < tmax;

				if (!debug)
				{
					htime->Fill(tc);
					httrig->Fill(tc - trig[m_fpga(frame[i])]);
					hctrig->Fill(tc - trig[m_fpga(frame[i])]);
					if (ch < 32)
						hcalib->Fill(ch + 32 * m_fpga(frame[i]), (tc - trig[m_fpga(frame[i])]));
				}

				ss.str("");
				ss.clear();
				ss << "/HitStudy/Rate/Count" << m_fpga(frame[i]);
				TH1 *hcount = _rh->GetTH1(ss.str());
				if (hcount == NULL)
				{
					std::cout << "Booking " << ss.str() << std::endl;
					hcount = _rh->BookTH1(ss.str(), 40, -0.1, 39.9);
				}
				hcount->Fill(ch);
				ss.str("");
				ss.clear();
				ss << "/HitStudy/Rate/Top/Count" << m_fpga(frame[i]);
				TH1 *hctop = _rh->GetTH1(ss.str());
				if (hctop == NULL)
				{
					std::cout << "Booking " << ss.str() << std::endl;
					hctop = _rh->BookTH1(ss.str(), 40, -0.1, 39.9);
				}
				if (side == 1 || ch > 31)
					hctop->Fill(ch);
				ss.str("");
				ss.clear();
				ss << "/HitStudy/Rate/Bottom/Count" << m_fpga(frame[i]);
				TH1 *hcbot = _rh->GetTH1(ss.str());
				if (hcbot == NULL)
				{
					std::cout << "Booking " << ss.str() << std::endl;
					hcbot = _rh->BookTH1(ss.str(), 40, -0.1, 39.9);
				}
				if (side == 0 || ch > 31)
					hcbot->Fill(ch);

				ss.str("");
				ss.clear();
				ss << "/HitStudy/Rate/Bottom/Strip";
				TH1 *hsbot = _rh->GetTH1(ss.str());
				if (hsbot == NULL)
					hsbot = _rh->BookTH1(ss.str(), 50, -0.1, 49.9);
				ss.str("");
				ss.clear();
				ss << "/HitStudy/Rate/Bottom/SelectedStrip";
				TH1 *hssbot = _rh->GetTH1(ss.str());
				if (hssbot == NULL)
					hssbot = _rh->BookTH1(ss.str(), 50, -0.1, 49.9);

				if (side == 0)
				{
					hsbot->Fill(strip);
					if (selected)
						hssbot->Fill(strip);
				}
				ss.str("");
				ss.clear();
				ss << "/HitStudy/Rate/Top/Strip";
				TH1 *hstop = _rh->GetTH1(ss.str());
				if (hstop == NULL)
					hstop = _rh->BookTH1(ss.str(), 50, -0.1, 49.9);
				ss.str("");
				ss.clear();
				ss << "/HitStudy/Rate/Top/SelectedStrip";
				TH1 *hsstop = _rh->GetTH1(ss.str());
				if (hsstop == NULL)
					hsstop = _rh->BookTH1(ss.str(), 50, -0.1, 49.9);

				if (side == 1)
				{
					hstop->Fill(strip);
					if (selected)
						hsstop->Fill(strip);
				}
			}
			// Delay studies and stip arry filling
			//printf("On passe 1\n");
			float tmin = TTRIG - 40;
			float tmax = TTRIG + 40;

			TProfile *hcalibc = (TProfile *)_rh->GetTH1("/Calibration/Corrected");
			if (hcalibc == NULL)
				hcalibc = _rh->BookProfile("/Calibration/Corrected", 400, -1.1, 98.9, -50., 50.);

			for (auto x = channels().begin(); x != channels().end(); x++)
			{
				TdcChannel *pch = (*x);
				std::stringstream ss;
				ss.str("");
				ss.clear();
				ss << "/Delay/FPGA" << pch->fpga();

				TH2 *hdist = _rh->GetTH2(ss.str() + "/dT");
				TH2 *hdist0 = _rh->GetTH2(ss.str() + "/dT0");
				TH2 *hdist1 = _rh->GetTH2(ss.str() + "/dT1");
				if (hdist == NULL)
				  {
				    hdist = _rh->BookTH2(ss.str() + "/dT", 3000, -2500, 500, 32, 0, 32);
				    hdist0 = _rh->BookTH2(ss.str() + "/dT0", 3000, -2500, 500, 32, 0, 32);
				    hdist1 = _rh->BookTH2(ss.str() + "/dT1", 3000, -2500, 500, 32, 0, 32);
				  }
				if (trig_ch[pch->fpga()] != 0)
				{
					hcalibc->Fill(32 * pch->fpga() + pch->channel(), pch->pedSubTime() - trig_ch[pch->fpga()]->tdcTime());
					double deltat = pch->pedSubTime() - trig_ch[pch->fpga()]->tdcTime();
					hdist->Fill(deltat, pch->channel());
					if (pch->side()==0)
					  hdist0->Fill(deltat, pch->channel());
					else
					  hdist1->Fill(deltat, pch->channel());
					// Flag selected channel
					if (deltat > tmin && deltat < tmax)
						pch->setUsed(true);
				}

				if (pch->used())
				{
					ss.str("");
					ss.clear();
					ss << "/Strip/Global";
					if (pch->side() == 1)
					{
						TH1 *hs = _rh->GetTH1(ss.str() + "/HRSelected");
						if (hs == NULL)
							hs = _rh->BookTH1(ss.str() + "/HRSelected", 50, -0.1, 49.9);
						hs->Fill(pch->strip());
					}
					else
					{
						TH1 *hs = _rh->GetTH1(ss.str() + "/LRSelected");
						if (hs == NULL)
							hs = _rh->BookTH1(ss.str() + "/LRSelected", 50, -0.1, 49.9);
						hs->Fill(pch->strip());
					}
					// fill strip vector of channels with selected channels
					c_strip[pch->strip()].push_back(pch);
				}
			}


			// Stip Study
			TH2 *hstrippos = _rh->GetTH2("/Strip/Global/TvsStrip");

			TH2 *hstripxy = _rh->GetTH2("/Strip/Global/XY");
			TH2 *hstripsy = _rh->GetTH2("/Strip/Global/YvsStrip");
			TH2 *hstriplty = _rh->GetTH2("/Strip/Global/LTY");
			TH2 *hstripltz = _rh->GetTH2("/Strip/Global/LTZ");
			TH1 *hstripcount = _rh->GetTH1("/Strip/Global/multiplicity");
			TH1 *hstrip10count = _rh->GetTH1("/Strip/InChamber/multiplicity");
			TH2 *hstrip10xy = _rh->GetTH2("/Strip/InChamber/XY");
			if (hstrippos == NULL)
			{
				hstrippos = _rh->BookTH2("/Strip/Global/TvsStrip", 50, -0.1, 49.9, 600, -30., 30.);
				hstripxy = _rh->BookTH2("/Strip/Global/XY", 100, 0., 80., 500, -50., 200.);
				hstrip10xy = _rh->BookTH2("/Strip/InChamber/XY", 100, 0., 80., 500, -50., 200.);
				hstripsy = _rh->BookTH2("/Strip/Global/YvsStrip", 50, -0.1, 49.9, 500, -50., 200.);
				hstriplty = _rh->BookTH2("/Strip/Global/LTY", 1000, 0., 400., 500, -50., 200.);
				hstripltz = _rh->BookTH2("/Strip/Global/LTZ", 1000, 0., 400., 500, -50., 200.);
				hstripcount = _rh->BookTH1("/Strip/Global/multiplicity", 50, -1., 49.);
				hstrip10count = _rh->BookTH1("/Strip/InChamber/multiplicity", 50, -1., 49.);
			}

			// Loop on all strips vector
			std::bitset<48> bs(0);
			std::bitset<48> bs10(0);
			for (int i = 0; i < 48; i++)
				// Ask only strips with 2 measurements
				if (c_strip[i].size() == 2)
				{
					// the strip is selected
					bs.set(i);
					float ts[2];
					TdcChannel *chh=NULL, *chl=NULL;
					for (int j = 0; j < 2; j++)
					{
						ts[c_strip[i][j]->side()] = c_strip[i][j]->pedSubTime();
						if (c_strip[i][j]->side() == 1)
							chh = c_strip[i][j];
						else
							chl = c_strip[i][j];
					}
					if (chh==NULL || chl==NULL) continue;
					// In Chamber condition: TL-TH > 10 ns
					bool inCh=(ts[0] - ts[1]) > 7 && chl->strip()>10 &&chl->strip()<30;
					if (inCh)
						bs10.set(i);
					// BUild a strip and fill histos
					TdcStrip a(0, chl->fpga(), chl->strip(), ts[0], ts[1], 0);
					if (dumpinfo)
						printf("Strip %d (%d,%d) HR %5.2f  (%d,%d)LR %5.2f delta(LR-HR) %5.2f X %5.2f Y %5.2f \n", i,
							   chh->fpga(), chh->channel(), ts[1],
							   chl->fpga(), chl->channel(), ts[0], ts[0] - ts[1], a.xpos(), a.ypos());

					hstrippos->Fill(i * 1., ts[0] - ts[1]);
					hstripxy->Fill(a.xpos(), a.ypos());
					if (inCh) hstrip10xy->Fill(a.xpos(), a.ypos());
					hstripsy->Fill(i*1., a.ypos());
					hstriplty->Fill(a.Lr(), a.ypos());
					hstripltz->Fill(a.Lr(), a.Zs());
				}
			//printf("%d %d \n",i,c_strip[i].size());
			#ifndef CALIBRATION
			if (info)
			  {
			    std::cout << run << " " << event << " " << bc0 << " " << bs << " " << bs.count() << std::endl;
			    std::cout << run << " " << event << " " << bc0 << " " << bs10 << " " << bs10.count() << " " << newevt << std::endl;
			  }
			#endif
			// At least one strip
			if (bs.count() > 0)
			{
				hcounter->Fill(2.);
				hstripcount->Fill(bs.count());
			}
			// At least one strip in chamber acceptance
			if (bs10.count() > 0)
			{
				hcounter->Fill(3.);
				hstrip10count->Fill(bs10.count());
			}
		}
		// clear memory
		for (int i = 0; i < 3; i++)
			if (trig_ch[i] != 0)
				delete trig_ch[i];

		//if (debug) getchar();

		// if (Cut(ientry) < 0) continue;
	}
}
#ifdef AFAIRE
void FebAna::fillTimePedestal(std::vector<lydaq::TdcChannel *> c_strip[])
{

	for (int is = 1; is < 32; is++)
		if (((bp >> is) & 1 == 1) && ((bp >> (is + 1)) & 1 == 1))
		{
			// Strip is et is+1 hits
			float t0 = 0, t0p = 0, t1 = 0, t1p = 0;
			for (int j = 0; j < c_strip[is].size(); j++)
			{
				TdcChannel *x = c_strip[is][j];
				if (x->side(_geo->feb(x->feb())) == 0)
					t0 = x->pedSubTime(_geo->feb(x->feb()));
				else
					t1 = x->pedSubTime(_geo->feb(x->feb()));
			}
			for (int j = 0; j < c_strip[is + 1].size(); j++)
			{
				TdcChannel *x = c_strip[is + 1][j];
				if (x->side(_geo->feb(x->feb())) == 0)
					t0p = x->pedSubTime(_geo->feb(x->feb()));
				else
					t1p = x->pedSubTime(_geo->feb(x->feb()));
			}

			std::stringstream srcs;
			srcs << "/Align/Pedestal" << is + 1;
			//std:cout<<srcs.str()<<std::endl;
			TH1 *hdt0 = _rh->GetTH1(srcs.str() + "HR");
			TH1 *hdt1 = _rh->GetTH1(srcs.str() + "LR");
			if (hdt0 == NULL)
			{
				hdt0 = _rh->BookTH1(srcs.str() + "HR", 90, -30., 30.);
				hdt1 = _rh->BookTH1(srcs.str() + "LR", 90, -30., 30.);
			}

			hdt0->Fill(t0p - t0);
			hdt1->Fill(t1p - t1);
		}
}
#endif
