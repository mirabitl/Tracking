import sqlite3
import ROOT
import ROOT.TF1
import json
import math
import time
from array import array
import copy
import numpy as np
import os
class prettyfloat(float):
    def __repr__(self):
        return "%0.2f" % self


defped=[31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31]
#defped=[30, 30, 34, 29, 32, 32, 37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 33, 27, 31, 29, 25, 35, 33, 9]
#defped=[39,38,42,37,41,40,46,0,43,48,37,49,0,0,0,0,0,0,0,0,39,39,36,35,41,35,40,34,32,44,43,34]

defped=[38,38,41,39,41,40,45,0,42,47,37,47,0,0,0,0,0,0,0,0,38,38,35,36,42,36,40,36,33,43,42,35]
def calP(p0,alt):
  return p0*(1-0.0065*alt/288.15)**5.255

def calV(V,P,T):
  #print 1-(0.2+0.8*P/990.*293./T)
  return  V/(0.2+0.8*P/990.*293./T)

def calValice(V,P,T):
  print (1-(1*P/990.*293./T))
  return  V/(1.0*P/990.*293./T)

def calApp(V,P,T):
  print (1-(0.2+0.8*P/990.*293./T))
  return  V*(0.2+0.8*P/990.*293./T)

def createTree(scan,idx):
    strcmd="./produceTree 1000 HV_%d_SN_%d" % (idx,scan)
    os.system(strcmd)
    fmn="/tmp/proc%d_%d.C" % (idx,scan)
    fout=open(fmn,"w");
    fout.write("{");
    fout.write('FebAna mon("/data/beamdump/root_trees/HV_%d_SN_%d_1000.root");' % (idx,scan));
    fout.write("mon.Loop();");
    fout.write('mon.writeHistograms("Results/ResHV_%d_SN_%d_1000.root");' % (idx,scan))
    fout.write("exit(0);}");
    fout.close()
    os.system("root -l -b %s" % fmn)
    
def getInfo(scan,hvr,reprocess=False,sub="./Results",debug=False,comment="",source=0):
    res_all={}
    # Chamber surface
    SCH=6022.0
    if (debug):
        c=ROOT.TCanvas("irpc","Alignement Studies",545,342)
        c.cd()
    else:
        c=None
    for i in range(len(hvr)):
        if (reprocess):
            createTree(scan,i+1)
        f=ROOT.TFile(sub+"/ResHV_%d_SN_%d_1000.root" % (i+1,scan))
        
        hcnt=f.Get("Counter")
        if (hcnt ==None):
            continue;
        ham=f.Get("/Strip/Global/multiplicity")
        hcm=f.Get("/Strip/InChamber/multiplicity")
        ntrg=hcnt.GetBinContent(1)
        if (ntrg<100):
            continue

        none=hcnt.GetBinContent(2)
        ninc=hcnt.GetBinContent(3)
        mul_all=int(ham.GetMean()*10)/10.
        mul_inc=int(hcm.GetMean()*10)/10.
        EffMax=100*none/ntrg
        EffInc=100*ninc/ntrg
        trg0=[]
        bkg0=[]
        trg1=[]
        bkg1=[]
        fopt="Q"
        if (debug):
            fopt=""
        for fpga in range(3):
            if (EffMax<20):
                continue
            hdt0=f.Get("/Delay/FPGA%d/dT0" %fpga)
            #print(hdt0.GetEntries())
            hdt0x=hdt0.ProjectionX()
            hdt0x.Rebin(4)

            xm=-100000.
            maxh=0
            hst=hdt0x
            for ib in range(1,hst.GetNbinsX()):
                if (hst.GetBinContent(ib)>maxh):
                    maxh=hst.GetBinContent(ib)
                    xm=hst.GetXaxis().GetBinCenter(ib)
            scfit=ROOT.TF1("scfit","[3]+[2]*TMath::Gaus(x,[0],[1])",xm-160,xm-160)
            scfit.SetParameter(0,xm)
            scfit.SetParameter(1,10)
            
            hdt0x.Fit("scfit",fopt,"",xm-160,xm+160);
            dt0mean=scfit.GetParameter(0)
            dt0res=scfit.GetParameter(1)
            dt0bkg=(scfit.GetParameter(3)/SCH/ntrg/4E-9/mul_all/EffMax*100)
            #print(" back ",scfit.GetParameter(3),SCH,ntrg,mul_all,dt0bkg)
            if (debug):
                hdt0x.Draw()
                c.Update()
                print("Hi " + input("Say something: "))
            hdt1=f.Get("/Delay/FPGA%d/dT1" %fpga)
            hdt1x=hdt1.ProjectionX()
            hdt1x.Rebin(6)
            xm=hdt1x.GetMean()
            xm=-100000.
            maxh=0
            hst=hdt1x
            for ib in range(1,hst.GetNbinsX()):
                if (hst.GetBinContent(ib)>maxh):
                    maxh=hst.GetBinContent(ib)
                    xm=hst.GetXaxis().GetBinCenter(ib)
            scfit=ROOT.TF1("scfit","[3]+[2]*TMath::Gaus(x,[0],[1])",xm-160,xm-160)
            scfit.SetParameter(0,xm)
            scfit.SetParameter(1,10)
            
            hdt1x.Fit("scfit",fopt,"",xm-160,xm+160);
            dt1mean=scfit.GetParameter(0)
            dt1res=scfit.GetParameter(1)
            dt1bkg=scfit.GetParameter(3)
            dt1bkg=(scfit.GetParameter(3)/SCH/ntrg/4E-9/mul_all/EffMax*100)
            if (debug):
                hdt1x.Draw()
                c.Update()
                print("Hi " + input("Say something: "))

            trg0.append(int(dt0mean*10)/10.)
            bkg0.append(int(dt0bkg*10)/10.)
            trg1.append(int(dt1mean*10)/10.)
            bkg1.append(int(dt1bkg*10)/10.)

        res={}
        res["scan"]=scan
        res["id"]=i+1
        res["hv"]=hvr[i]
        res["orbits"]=int(ntrg)
        res["OneInWindow"]=int(none)
        res["OneInChamber"]=int(ninc)
        res["EffWindow"]=int(EffMax*10)/10.
        res["EffChamber"]=int(EffInc*10)/10.
        res["MultiplicityWindow"]=mul_all
        res["MultiplicityChamber"]=mul_inc
        res["dt0"]=trg0
        res["dt1"]=trg1
        res["rate0"]=bkg0
        res["rate1"]=bkg1
        res["source"]=source
        res["comment"]=comment
        
        print(scan,i+1,hvr[i],int(ntrg),int(none),int(ninc),int(EffMax*10)/10.,int(EffInc*10)/10.,mul_all,mul_inc,trg0,trg1,bkg0,bkg1,"SOURCE ",source,comment)
        if (debug):
            print(res)
        res_all["HV_%d_SN_%d" % (i+1,scan)]=res
        f.Close()
    return res_all
        
def getdy(run,hn="strip",sub="/tmp/"):
  c=ROOT.TCanvas("irpc","Alignement Studies",545,342)
  f82=ROOT.TFile(sub+"histo%d_0.root" % run);
  f82.cd("/Align")
  gStyle.SetOptFit();
  gStyle.SetOptStat(0);
  res=[]
  dres=[]
  for i in range(48):
    res.append(0)
    dres.append(0)
  for ist in range(1,33):
    hst=f82.Get("/Align/%s%d" % (hn,ist))
    
    if (hst==None):
      continue
    if (hst.GetEntries()<100):
      hst.Rebin(2)
    hst.SetTitle("Y_{ext} - Y_{strip} Strip %d" % ist)
    hst.GetXaxis().SetTitle("#Delta Y (cm)")
    hst.GetYaxis().SetTitle("Number of event")
    xm=-100000.
    maxh=0
    for ib in range(1,hst.GetNbinsX()):
      if (hst.GetBinContent(ib)>maxh):
        maxh=hst.GetBinContent(ib)
        xm=hst.GetXaxis().GetBinCenter(ib)
    scfit=TF1("scfit","gaus",xm-10,xm-10)
    hst.Fit("scfit","Q","",xm-20,xm+20);
    dtmean=scfit.GetParameter(1)
    dtres=scfit.GetParameter(2)
    print(ist,hst.GetEntries(),hst.GetMean(),xm,dtmean,dtres)
    res[ist]=dtmean
    dres[ist]=dtres
    c.cd()
    hst.Draw()
    c.Modified()
    c.Update()
    c.SaveAs("IRPC_dy_strip%d.png" % ist)
  val=raw_input()
  y=np.array(res)
  dy=np.array(dres)
  #np.set_printoptions(precision=1)
  np.set_printoptions(formatter={'float': '{: 0.2f}'.format})
  print(y)
  print(dy)

def getratio(run,chamber=1,sub="/tmp/",rebin=1,ncut=30):
  c=ROOT.TCanvas("irpc","IRPC Studies",545,842)
  c.Divide(2,4)
  f82=ROOT.TFile(sub+"histo%d_0.root" % run);
  f82.cd("/FEB/Chamber%d/Raw" % chamber)
  hst=f82.Get("/FEB/Chamber%d/Raw/Strips" % chamber)
  c.cd(1)
  hst.Draw()
  hrat=ROOT.TH1F("hrat","Ratio LR/HR ",50,0,50.)
  for i in range(1,33):
    x0=hst.GetBinContent(i+1)
    x1=hst.GetBinContent(i+1+48)
    r=0
    if (x0>x1):
      r=x1/x0
      #print x0,x1,r
    hrat.SetBinContent(i+1,r)
  c.cd(2)
  hrat.Draw()
  c.Modified()
  c.Update()
  #val=raw_input()
  f82.cd("/gric")
  hxy=f82.Get("/gric/XY")
  hxy.SetTitle(" 4 points Track extrapolation to the iRPC")
  hxyf=f82.Get("/gric/XYF")
  hxyf14=f82.Get("/gric/XYF14")
  hxyf15=f82.Get("/gric/XYF15")
  hxyf.SetTitle(" 4 points Track extrapolation to the iRPC when cluster found nearby (5cm,20cm)")
  hxyt=f82.Get("/gric/XYT")
  result=[]
  efft=hxyt.GetEntries()/hxy.GetEntries()
  result.append(run)
  result.append(hxy.GetEntries())
  result.append(hxyt.GetEntries())
  result.append(efft*100)
  effo=hxyf.GetEntries()/hxy.GetEntries()
  effo14=hxyf14.GetEntries()/hxy.GetEntries()
  effo15=hxyf15.GetEntries()/hxy.GetEntries()
  result.append(hxyf.GetEntries())
  result.append(effo*100)
  result.append(hxyf14.GetEntries())
  result.append(effo14*100)
  result.append(hxyf15.GetEntries())
  result.append(effo15*100)
  #print hxy.GetEntries(),efft*100,effo*100,effo14*100,effo15*100
  c.cd(3)
  hxy.Draw("COLZ")
  c.cd(4)
  hxyf.Draw("COLZ")
  rbx=rebin
  rby=2
  hxy.Rebin2D(rbx,rby)
  hxyf.Rebin2D(rbx,rby)
  hxyf14.Rebin2D(rbx,rby)
  hxyf15.Rebin2D(rbx,rby)
  hxyt.Rebin2D(rbx,rby)
  hxy.SetAxisRange(0.,35.,"X")
  hxy.SetAxisRange(0.,60.,"Y")
  hxyf.SetAxisRange(0.,35.,"X")
  hxyf.SetAxisRange(0.,60.,"Y")
  hxyf14.SetAxisRange(0.,35.,"X")
  hxyf14.SetAxisRange(0.,60.,"Y")
  hxyf15.SetAxisRange(0.,35.,"X")
  hxyf15.SetAxisRange(0.,60.,"Y")
  hxyEff=hxyf.Clone("hxyEff")
  hxyEff.Divide(hxy)
  hxyEff.SetTitle("Local Efficiency map")
  hxyEff14=hxyf14.Clone("hxyEff14")
  hxyEff14.Divide(hxy)
  hxyEff14.SetTitle("Local Efficiency map FR4")
  hxyEff15=hxyf15.Clone("hxyEff15")
  hxyEff15.Divide(hxy)
  hxyEff15.SetTitle("Local Efficiency map EM888")
  heff=ROOT.TH1F("heff","Local Efficiency Ntk ext>15",110,0,1.1)
  dmin=max(0,effo*0.7)
  dmax=min(1.1,effo*1.6)
  nb=int((dmax-dmin)/0.005)
  heff1=ROOT.TH1F("heff1","Local Efficiency EM888",nb,dmin,dmax)
  heff32=ROOT.TH1F("heff32","Local Efficiency FR4",nb,dmin,dmax)

  for i in range(1,hxy.GetNbinsX()):
    for j in range(1,hxy.GetNbinsY()):
      if (hxy.GetBinContent(i,j)>ncut):
        heff.Fill(hxyEff.GetBinContent(i,j))
        xp=hxy.GetXaxis().GetBinCenter(i)
        yp=hxy.GetYaxis().GetBinCenter(j)
# Jusqu'au run 1721
        #if (xp<=13):
        #    heff1.Fill(hxyEff14.GetBinContent(i,j))
        #if (xp>=19):
        #    heff32.Fill(hxyEff15.GetBinContent(i,j))
# A partir du 1722 (14 et 15 inverted)
        if (xp>=7 and xp<=15 and yp<46 and yp>3):
            heff1.Fill(hxyEff15.GetBinContent(i,j))
            print (hxy.GetBinContent(i,j), hxyf15.GetBinContent(i,j),hxyEff15.GetBinContent(i,j))
        if (xp>=17 and xp<=26 and yp<46 and yp>3):
            heff32.Fill(hxyEff14.GetBinContent(i,j))

  print (heff1.GetMean()*100,heff32.GetMean()*100)

  deff1=(heff1.GetRMS()/math.sqrt(heff1.GetEntries()))*100
  deff32=(heff32.GetRMS()/math.sqrt(heff32.GetEntries()))*100
  result.append(heff1.GetEntries())
  result.append(heff1.GetMean()*100)
  result.append(deff1)
  result.append(heff32.GetEntries())
  result.append(heff32.GetMean()*100)
  result.append(deff32)
  y=np.array(result)
  #np.set_printoptions(precision=1)
  np.set_printoptions(formatter={'float': '{: 0.1f}'.format})
  y.tofile("result%d.csv" % run,sep='|',format='%7.1f')
  print (y)
  
  c.cd(5)
  hxyEff15.Draw("COLZ")
  c.Modified()
  c.Update()
  #val=raw_input()
  c.cd(6)
  hxyEff14.Draw("COLZ")
  c.Modified()
  c.Update()
  #val=raw_input()
  c.cd(7)
  heff1.Draw()
  c.Modified()
  c.Update()
  #val=raw_input()
  c.cd(8)
  heff32.Draw()
  c.Modified()
  c.Update()
  #
  c.SaveAs("Analyse%d.pdf" % run)
  val=raw_input()
def drawtdc(run,sub="./"):
    f82=ROOT.TFile(sub+"res%d.root" % run);
    f82.cd()
    hp=f82.Get("/Calibration/Distance2Trigger")
    fp={}
    fp[0]=[]
    fp[1]=[]
    fp[2]=[]
    i=5
    for k in range(0,10000):
        k=i-5
        ch=k/4

        fpga=ch//32
        lch=ch%32
        #print(i,ch,fpga,lch,hp.GetBinContent(i));
        fp[fpga].append(hp.GetBinContent(i))
        i=i+4
        if (ch==95):
            break
    fpc=copy.deepcopy(fp)
    for i in range(0,32):
        #print(i)
        fpc[0][i]=round((fp[0][i]-fp[0][0])*100)/100
        fpc[1][i]=round((fp[1][i]-fp[1][0])*100)/100
        fpc[2][i]=round((fp[2][i]-fp[2][0])*100)/100
    #print(fpc)
    
    s0="double dt0[32]={"
    for i in range(0,31):
        s0=s0+"%5.2f," % fpc[0][i]
    s0=s0+"%5.2f};"% fpc[0][31]
    s1="double dt1[32]={"
    for i in range(0,31):
        s1=s1+"%5.2f," % fpc[1][i]
    s1=s1+"%5.2f};"% fpc[1][31]
    s2="double dt2[32]={"
    for i in range(0,31):
        s2=s2+"%5.2f," % fpc[2][i]
    s2=s2+"%5.2f};"% fpc[2][31]
    print(s0)
    print(s1)
    print(s2)
def getdt(run,chamber,feb,sub=""):
  c=ROOT.TCanvas()
  f82=ROOT.TFile(sub+"histo%d_0.root" % run);
  r=[]
  for i in range(0,49):
    r.append(0)
  for i in range(1,49):
    #print "/run%d/Chamber%d/FEB/%d/Side0/channel%d" % (run,chamber,feb,i)
    f82.cd("/run%d/InTime/Chamber%d/FEB/%d/Side0/" % (run,chamber,feb))
    hch=f82.Get("/run%d/InTime/Chamber%d/FEB/%d/Side0/channel%d" % (run,chamber,feb,i))
    #print hch
    if (hch!=None):
      #print i,hch.GetEntries(),hch.GetMean();
      if (hch.GetEntries()>25):
          if (hch.GetEntries()<100):
              hch.Rebin(2)
          if (hch.GetEntries()<50):
              hch.Rebin(2)
          nmax=0
          imax=0
          for ib in range(1,hch.GetNbinsX()):
              if (hch.GetBinContent(ib)>nmax):
                  nmax=hch.GetBinContent(ib)
                  imax=ib

          vmin=hch.GetBinCenter(imax)-7;
          vmax=hch.GetBinCenter(imax)+7;
          print(vmin,vmax,imax,nmax)
          scfit=TF1("scfit","gaus",vmin,vmax)
          hch.Fit("scfit","Q","",vmin,vmax);
          dtmean=scfit.GetParameter(1)
          dtres=scfit.GetParameter(2)
          print (i,hch.GetEntries(),hch.GetMean(),dtmean,dtres)
          c.cd()
          hch.Draw()
          c.Modified()
          c.Update()
          val=raw_input()
          r[i]=dtmean
          if (abs(dtmean-hch.GetMean())>hch.GetRMS()):
              r[i]=hch.GetMean()
      continue
    
    f82.cd("/run%d/InTime/Chamber%d/FEB/%d/Side1/" % (run,chamber,feb))
    hch1=f82.Get("/run%d/InTime/Chamber%d/FEB/%d/Side1/channel%d" % (run,chamber,feb,i))
    #print hch1
    if (hch1!=None):
      #print i,hch1.GetEntries(),hch1.GetMean();
      if (hch1.GetEntries()>25):
          if (hch1.GetEntries()<100):
              hch1.Rebin(2)
          if (hch1.GetEntries()<50):
              hch1.Rebin(2)

          nmax=0
          imax=0
          for ib in range(1,hch1.GetNbinsX()):
              if (hch1.GetBinContent(ib)>nmax):
                  nmax=hch1.GetBinContent(ib)
                  imax=ib

          vmin=hch1.GetBinCenter(imax)-7;
          vmax=hch1.GetBinCenter(imax)+7;

          scfit=TF1("scfit","gaus",vmin,vmax)
          hch1.Fit("scfit","Q","",vmin,vmax);
          dtmean=scfit.GetParameter(1)
          dtres=scfit.GetParameter(2)
          print(i,hch1.GetEntries(),hch1.GetMean(),dtmean,dtres)
          c.cd()
          hch1.Draw()
          c.Modified()
          c.Update()
          val=raw_input()

          r[i]=dtmean
          #hch1.GetMean()
          if (abs(dtmean-hch1.GetMean())>hch1.GetRMS()):
              r[i]=hch1.GetMean()

      continue

  #print r
  hfeb0=f82.Get("/run%d/InTime/Chamber%d/FEB/%d/Side0/DTall" % (run,chamber,feb))
  hfeb1=f82.Get("/run%d/InTime/Chamber%d/FEB/%d/Side1/DTall" % (run,chamber,feb))
  sfit=TF1("sfit","gaus",-150.,-100.)
  hfeb0.Fit("sfit","Q","",-150.,-100.);
  c.cd()
  hfeb0.Draw()
  c.Modified()
  c.Update()
  #val=raw_input()
  print ("dt0: %5.1f" % (sfit.GetParameter(1)))
  hfeb1.Fit("sfit","Q","",-150.,-100.);
  c.cd()
  hfeb1.Draw()
  c.Modified()
  c.Update()
  #val=raw_input()
  print("dt1: %5.1f"  % sfit.GetParameter(1))
  
  r = map(prettyfloat, r)
  print ('"dtc":',r)
def fitProfile(run,sel=92):
  f82=ROOT.TFile("./histo%d_0.root" % run);
  #f82.cd("/run%d/TDC%d/LmAnalysis/Timing" % (run,tdc));
  f82.cd("/run%d/ChamberAll" % (run));
  c1=ROOT.TCanvas();
  gStyle.SetOptFit();
  pos0=[]
  pmean0=[]
  pos1=[]
  pmean1=[]
  hall=f82.Get("/run%d/ChamberAll/XY" % (run));
  for i in range(sel,sel+1):
    pos0.append(0)
    pmean0.append(0)
    #hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/Timing/hdtpos%d" % (run,tdc,i+71));
    hstrip=hall.ProjectionY("strip%d" % (i),(i),i+1)
    if (hstrip.GetEntries()<100):
      continue
    scfit=TF1("scfit","gaus",hstrip.GetMean()-2.*hstrip.GetRMS(),hstrip.GetMean()+2.*hstrip.GetRMS())
    hstrip.Fit("scfit","Q","",hstrip.GetMean()-2.*hstrip.GetRMS(),hstrip.GetMean()+2.*hstrip.GetRMS())
    dtmean=scfit.GetParameter(1)
    dtres=scfit.GetParameter(2)
    print (run,i,dtmean,dtres,dtmean*80./8.3,dtres*80/8.3)
    c1.Update()
    #val = raw_input()
    time.sleep(2)


def fitdif(run):
  f82=ROOT.TFile("./histo%d_0.root" % run);
  #f82.cd("/run%d/TDC%d/LmAnalysis/Timing" % (run,tdc));
  f82.cd("/run%d/ChamberDif" % (run));
  c1=ROOT.TCanvas();
  gStyle.SetOptFit();
  pos0=[]
  pmean0=[]
  pos1=[]
  pmean1=[]
  for i in range(48):
    pos0.append(0)
    pmean0.append(0)
    #hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/Timing/hdtpos%d" % (run,tdc,i+71));
    hstrip=f82.Get("/run%d/Timing/OneStrip/Side0/hdtr_%d" % (run,i+71));
    if (hstrip == None):
      continue
    if (hstrip.GetEntries()<20):
      continue
    hstrip.Rebin(2)
    
    hstrip.Draw()
    c1.Update()

    print ("Enter min max")
    #hmin = float(raw_input())
    #hmax = float(raw_input())
    hmin=hstrip.GetMean()-5.*hstrip.GetRMS()
    hmax=hstrip.GetMean()+5.*hstrip.GetRMS()
    print (hmin,hmax)
    #scfit=TF1("scfit","gaus",hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    #hstrip.GetXaxis().SetRangeUser(hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    
    scfit=TF1("scfit","gaus",hmin,hmax)
    hstrip.GetXaxis().SetRangeUser(hmin,hmax)
    hstrip.Fit("scfit","Q");
    dtmean=scfit.GetParameter(1)
    dtres=scfit.GetParameter(2)
    hstrip.Draw()
    c1.Update()
    #c1.SaveAs("Run%d_Strip_pos.png" % (run));

    val = raw_input()
    pmean0[i]=hstrip.GetMean()
    if (dtres<hstrip.GetRMS()):
      pos0[i]=dtmean
    else:
      pos0[i]=hstrip.GetMean()
      
  print (pos0)
  print(pmean0)
  for i in range(48):
    pos1.append(0)
    pmean1.append(0)
    #hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/Timing/hdtpos%d" % (run,tdc,i+71));
    hstrip=f82.Get("/run%d/Timing/OneStrip/Side1/hdtr_%d" % (run,i+71));
    if (hstrip == None):
      continue
    if (hstrip.GetEntries()<20):
      continue
    hstrip.Rebin(2)
    
    hstrip.Draw()
    c1.Update()

    print("Enter min max")
    #hmin = float(raw_input())
    #hmax = float(raw_input())
    hmin=hstrip.GetMean()-5.*hstrip.GetRMS()
    hmax=hstrip.GetMean()+5.*hstrip.GetRMS()
    print (hmin,hmax)
    #scfit=TF1("scfit","gaus",hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    #hstrip.GetXaxis().SetRangeUser(hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    
    scfit=TF1("scfit","gaus",hmin,hmax)
    hstrip.GetXaxis().SetRangeUser(hmin,hmax)
    hstrip.Fit("scfit","Q");
    dtmean=scfit.GetParameter(1)
    dtres=scfit.GetParameter(2)
    hstrip.Draw()
    c1.Update()
    #c1.SaveAs("Run%d_Strip_pos.png" % (run));

    val = raw_input()
    pmean1[i]=hstrip.GetMean()
    if (dtres<hstrip.GetRMS()):
      pos1[i]=dtmean
    else:
      pos1[i]=hstrip.GetMean()
      
  print (pos0)
  print (pmean0)





  
  dt=0
  for i in range(48):
    print( "fe1_2tr[%d]=%5.3f;" % (i+71,pos0[i]-pos1[i]))

  #print pos



