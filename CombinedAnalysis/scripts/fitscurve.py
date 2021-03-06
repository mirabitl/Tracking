import sqlite3
from ROOT import *
import json
import math
import time
from array import array
import numpy as np
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
  print 1-(1*P/990.*293./T)
  return  V/(1.0*P/990.*293./T)

def calApp(V,P,T):
  print 1-(0.2+0.8*P/990.*293./T)
  return  V*(0.2+0.8*P/990.*293./T)

def getdy(run,hn="strip",sub="/tmp/"):
  c=TCanvas("irpc","Alignement Studies",545,342)
  f82=TFile(sub+"histo%d_0.root" % run);
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
    hst.Fit("scfit","","",xm-20,xm+20);
    dtmean=scfit.GetParameter(1)
    dtres=scfit.GetParameter(2)
    print ist,hst.GetEntries(),hst.GetMean(),xm,dtmean,dtres
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
  print y
  print dy

def getratio(run,chamber=1,sub="/tmp/",rebin=1,ncut=30):
  c=TCanvas("irpc","IRPC Studies",545,842)
  c.Divide(2,4)
  f82=TFile(sub+"histo%d_0.root" % run);
  f82.cd("/FEB/Chamber%d/Raw" % chamber)
  hst=f82.Get("/FEB/Chamber%d/Raw/Strips" % chamber)
  c.cd(1)
  hst.Draw()
  hrat=TH1F("hrat","Ratio LR/HR ",50,0,50.)
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
  heff=TH1F("heff","Local Efficiency Ntk ext>15",110,0,1.1)
  dmin=max(0,effo*0.7)
  dmax=min(1.1,effo*1.6)
  nb=int((dmax-dmin)/0.005)
  heff1=TH1F("heff1","Local Efficiency EM888",nb,dmin,dmax)
  heff32=TH1F("heff32","Local Efficiency FR4",nb,dmin,dmax)

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
            print hxy.GetBinContent(i,j), hxyf15.GetBinContent(i,j),hxyEff15.GetBinContent(i,j)
        if (xp>=17 and xp<=26 and yp<46 and yp>3):
            heff32.Fill(hxyEff14.GetBinContent(i,j))

  print heff1.GetMean()*100,heff32.GetMean()*100

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
  print y
  
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
def drawtdc(run,sub="Histos/InTime/"):
    c=TCanvas()
    hmean=TH1F("hmean","Summary of mean distance ",1000,-20.,20.)
    hrms=TH1F("hrms","Summary of dispersion ",200,0,0.04)
    f82=TFile(sub+"histo%d_0.root" % run);
    f82.cd("/run%d/Raw/DT/" % run )
    for i in range(64):
        for j in range(i+1,64):
            hch=f82.Get("/run%d/Raw/DT/ch_%d_%d" % (run,i,j))
            if (hch==None):
                continue
            m=hch.GetMean()
            r=hch.GetRMS()
            hmean.Fill(m)
            hrms.Fill(r)
    c.cd()
    hmean.Draw()
    c.SaveAs("Dispersion_channel_to_channel_%d.png" % run)
    v=raw_input()
    
    hrms.Draw()
    c.Update()
    c.SaveAs("RMS_channel_to_channel_%d.png" % run)
    v=raw_input()
def getdt(run,chamber,feb,sub=""):
  c=TCanvas()
  f82=TFile(sub+"histo%d_0.root" % run);
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
          print vmin,vmax,imax,nmax
          scfit=TF1("scfit","gaus",vmin,vmax)
          hch.Fit("scfit","Q","",vmin,vmax);
          dtmean=scfit.GetParameter(1)
          dtres=scfit.GetParameter(2)
          print i,hch.GetEntries(),hch.GetMean(),dtmean,dtres
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
          print i,hch1.GetEntries(),hch1.GetMean(),dtmean,dtres
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
  print "dt0: %5.1f" % (sfit.GetParameter(1))
  hfeb1.Fit("sfit","Q","",-150.,-100.);
  c.cd()
  hfeb1.Draw()
  c.Modified()
  c.Update()
  #val=raw_input()
  print "dt1: %5.1f"  % sfit.GetParameter(1)
  
  r = map(prettyfloat, r)
  print '"dtc":',r
def fitProfile(run,sel=92):
  f82=TFile("./histo%d_0.root" % run);
  #f82.cd("/run%d/TDC%d/LmAnalysis/Timing" % (run,tdc));
  f82.cd("/run%d/ChamberAll" % (run));
  c1=TCanvas();
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
    print run,i,dtmean,dtres,dtmean*80./8.3,dtres*80/8.3
    c1.Update()
    #val = raw_input()
    time.sleep(2)


def fitdif(run):
  f82=TFile("./histo%d_0.root" % run);
  #f82.cd("/run%d/TDC%d/LmAnalysis/Timing" % (run,tdc));
  f82.cd("/run%d/ChamberDif" % (run));
  c1=TCanvas();
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

    print "Enter min max"
    #hmin = float(raw_input())
    #hmax = float(raw_input())
    hmin=hstrip.GetMean()-5.*hstrip.GetRMS()
    hmax=hstrip.GetMean()+5.*hstrip.GetRMS()
    print hmin,hmax
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
      
  print pos0
  print pmean0
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

    print "Enter min max"
    #hmin = float(raw_input())
    #hmax = float(raw_input())
    hmin=hstrip.GetMean()-5.*hstrip.GetRMS()
    hmax=hstrip.GetMean()+5.*hstrip.GetRMS()
    print hmin,hmax
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
      
  print pos0
  print pmean0





  
  dt=0
  for i in range(48):
    print "fe1_2tr[%d]=%5.3f;" % (i+71,pos0[i]-pos1[i])

  #print pos

def fitpos(run,tdc):
  f82=TFile("./histo%d_0.root" % run);
  #f82.cd("/run%d/TDC%d/LmAnalysis/Timing" % (run,tdc));
  f82.cd("/run%d/Chamber%d/Timing" % (run,tdc));
  c1=TCanvas();
  gStyle.SetOptFit();
  pos=[]
  pmean=[]
  for i in range(48):
    pos.append(0)
    pmean.append(0)
    #hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/Timing/hdtpos%d" % (run,tdc,i+71));
    hstrip=f82.Get("/run%d/Chamber%d/Timing/All/hdtpos%d" % (run,tdc,i+71));
    if (hstrip == None):
      continue
    hstrip.Rebin(2)
    
    hstrip.Draw()
    c1.Update()

    print "Enter min max"
    #hmin = float(raw_input())
    #hmax = float(raw_input())
    hmin=-20.
    hmax=20.
    print hmin,hmax
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
    pmean[i]=hstrip.GetMean()
    if (dtres>1. and  dtres<4):
      pos[i]=dtmean
    else:
      pos[i]=hstrip.GetMean()
  print pos
  print pmean
  
  for i in range(12):
    print "fe1_2tr[%d]=%5.3f;" % (i,pos[i])
  #print pos

def calcefn(run,chamber,hv=0,dirp="."):
  f82=TFile("%s/histo%d_0.root" % (dirp,run));
  f82.cd("/run%d/Chamber%d/" % (run,chamber));
  c1=TCanvas();
  gStyle.SetOptFit();

  hns=f82.Get("/run%d/Chamber%d/Efficiency" % (run,chamber));
  #hns2=f82.Get("/run%d/TDC%d/LmAnalysis/NStrips2" % (run,tdc));
  hstrip=f82.Get("/run%d/Chamber%d/NSTRIP" % (run,chamber));
  hrate=f82.Get("/run%d/Chamber%d/Rate" % (run,chamber));

  hclusters=f82.Get("/run%d/Chamber%d/ClusterNew/Clusters" % (run,chamber));
  hfr=f82.Get("/run%d/Chamber%d/FebCount" % (run,chamber));
  hfrs=f82.Get("/run%d/Chamber%d/FebCountSel" % (run,chamber));
  hclusterm=f82.Get("/run%d/Chamber%d/ClusterNew/ClusterSize1" % (run,chamber));
  #hstrip.Rebin(2)nCont(I)
  febrate=0
  febs=(0,6500./4,7500./4)
  #print hfr
  febratesel=0
  if (hfr!=None):
    nevt_fr=hfr.GetBinContent(25);
    nmax=0
    nfb=0
    for i in range(0,24):
      #print i,hfr.GetBinContent(i),nevt_fr
      if (hfr.GetBinContent(i)>0):
        nmax=nmax+hfr.GetBinContent(i)
        nfb=nfb+1
    #print nevt_fr,nfb,nmax
    if (nevt_fr>0 and nfb>0):
      febrate=nmax*1./nfb/(nevt_fr*20E-9)
  if (hfrs!=None):
    nevt_frs=hfrs.GetBinContent(25);
    nmaxs=0
    nfbs=0
    for i in range(0,24):
      #print i,hfr.GetBinContent(i),nevt_fr
      if (hfrs.GetBinContent(i)>0):
        nmaxs=nmaxs+hfrs.GetBinContent(i)
        nfbs=nfbs+1
    #print nevt_frs,nfbs,nmaxs
    if (nevt_frs>0 and nfbs>0):
      febratesel=nmaxs*1./nfbs/(nevt_frs*20E-9)


  csize=0.1
  effc=0.0
  deffc=0.0
  ncev=0
  ncl=0
  if (hclusterm!=None):
    csize=hclusterm.GetMean()
    ncev=hclusters.GetEntries()
    nc=ncev-hclusters.GetBinContent(1)
    ncl=0
    nevcl=0
    for i in range(2,16):
        x=i-1.;
        y=hclusters.GetBinContent(i)
        nevcl=nevcl+y
        ncl=ncl+x*y
    if (nevcl>0):
        ncl=ncl*1./nevcl
    effc=nc*1./ncev
    deffc=math.sqrt(effc*(1-effc)/ncev)
  ntrg=hns.GetBinContent(3)
  nall=hns.GetBinContent(4)
  nxy=hns.GetBinContent(5)

  hstrip.GetXaxis().SetRangeUser(0.5,9.5)
  hrate.GetXaxis().SetRangeUser(0.1,60000.)
  mul=hstrip.GetMean()

  
  #print ntrg,hns.GetBinContent(1),hns2.GetBinContent(1)
  eff=nall*1./ntrg
  deff=math.sqrt(eff*(1-eff)/ntrg)
  effp=nxy*1./ntrg
  deffp=math.sqrt(effp*(1-effp)/ntrg)
  
  print "|%d|%d|%7.1f|%d|%d|%d|%5.2f|%5.2f|%5.2f|%5.2f|%5.1f|%7.1f|%d|%5.2f|%5.2f|%5.2f|%5.2f|%5.1f|%5.2f|%5.1f" % (run,chamber,hv,int(ntrg),int(nall),int(nxy),eff*100,deff*100,effp*100,deffp*100,mul,hrate.GetMean(),ncev,effc*100,deffc*100,csize,ncl,febrate/febs[chamber],-febrate*10E-7,febratesel/febs[chamber])
  #hstrip.Draw()
  #c1.Update()
  #c1.SaveAs("Run%d_Strip_pos.png" % (run));

  #val = raw_input()
  r=(run,chamber,hv,int(ntrg),int(nall),int(nxy),eff*100,deff*100,effp*100,deffp*100,mul,hrate.GetMean(),ncev,effc*100,deffc*100,csize,ncl,febrate/febs[chamber],-febrate*10E-7,febratesel/febs[chamber])
  #r = map(prettyfloat, r)
  return r


def extractEfficiency(run,chamber,hv=0,dirp="."):
  f82=TFile("%s/histo%d_0.root" % (dirp,run));
  f82.cd("/run%d/InTime/Chamber%d/" % (run,chamber));
  #c1=TCanvas();
  #gStyle.SetOptFit();

  hclusters=f82.Get("/run%d/InTime/Chamber%d/ClusterNew/Clusters" % (run,chamber));
  hclustero=f82.Get("/run%d/OffTime/Chamber%d/ClusterNew/Clusters" % (run,chamber));

  hfr=f82.Get("/run%d/InTime/Chamber%d/FebCount" % (run,chamber));
  hfrs=f82.Get("/run%d/OffTime/Chamber%d/FebCountSel" % (run,chamber));
  hclusterm=f82.Get("/run%d/InTime/Chamber%d/ClusterNew/ClusterSize1" % (run,chamber));
  #hstrip.Rebin(2)nCont(I)
  febrate=0
  febs=(0,6500./4,7500./4)
  #print hfr
  febratesel=0
  if (hfr!=None):
    nevt_fr=hfr.GetBinContent(25);
    nmax=0
    nfb=0
    for i in range(0,24):
      #print i,hfr.GetBinContent(i),nevt_fr
      if (hfr.GetBinContent(i)>0):
        nmax=nmax+hfr.GetBinContent(i)
        nfb=nfb+1
    #print nevt_fr,nfb,nmax
    if (nevt_fr>0 and nfb>0):
      febrate=nmax*1./nfb/(nevt_fr*20E-9)
  if (hfrs!=None):
    nevt_frs=hfrs.GetBinContent(25);
    nmaxs=0
    nfbs=0
    for i in range(0,24):
      #print i,hfr.GetBinContent(i),nevt_fr
      if (hfrs.GetBinContent(i)>0):
        nmaxs=nmaxs+hfrs.GetBinContent(i)
        nfbs=nfbs+1
    #print nevt_frs,nfbs,nmaxs
    if (nevt_frs>0 and nfbs>0):
      febratesel=nmaxs*1./nfbs/(nevt_frs*20E-9)


  csize=0.1
  effc=0.0
  deffc=0.0
  effco=0.0
  deffco=0.0
  ncevo=0
  ncev=0
  ncl=0
  if (hclusterm!=None):
    csize=hclusterm.GetMean()
    ncev=hclusters.GetEntries()
    nc=ncev-hclusters.GetBinContent(1)
    ncl=0
    nevcl=0
    for i in range(2,32):
        x=i-1.;
        y=hclusters.GetBinContent(i)
        nevcl=nevcl+y
        ncl=ncl+x*y
    if (nevcl>0):
        ncl=ncl*1./nevcl
    effc=nc*1./ncev
    deffc=math.sqrt(effc*(1-effc)/ncev)
    ncevo=hclustero.GetEntries()
    nco=ncevo-hclustero.GetBinContent(1)
    effco=nco*1./ncevo
    if (ncevo>0 and effco>0):
        #print ncevo,effco
        deffco=math.sqrt(effco*(1-effco)/ncevo)

  
  print "|%d|%d|%7.1f|%d|%5.2f|%5.2f|%5.2f|%5.2f|%5.1f|%5.2f|%5.1f|%5.2f|%5.2f|" % (run,chamber,hv,ncev,effc*100,deffc*100,csize,ncl,febrate/febs[chamber],-febrate*10E-7,febratesel/febs[chamber],effco*100,100*(effc-effco)/(1.-effco))
  #hstrip.Draw()
  #c1.Update()
  #c1.SaveAs("Run%d_Strip_pos.png" % (run));

  #val = raw_input()
  r=(run,chamber,hv,ncev,effc*100,deffc*100,csize,ncl,febrate/febs[chamber],-febrate*10E-7,febratesel/febs[chamber],effco*100,100*(effc-effco)/(1.-effco))
  #r = map(prettyfloat, r)
  return r



def calceff(run,tdc,strip=71):
  f82=TFile("./histo%d_0.root" % run);
  f82.cd("/run%d/TDC%d/LmAnalysis" % (run,tdc));
  c1=TCanvas();
  gStyle.SetOptFit();

  heff=f82.Get("/run%d/TDC%d/LmAnalysis/Efficiency" % (run,tdc));
  hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/hdts%d" % (run,tdc,strip));
  hstrip.Rebin(2)
  scfit=TF1("scfit","gaus",hstrip.GetMean()-6.*hstrip.GetRMS(),hstrip.GetMean()+6.*hstrip.GetRMS())
  hstrip.GetXaxis().SetRangeUser(hstrip.GetMean()-6.*hstrip.GetRMS(),hstrip.GetMean()+6.*hstrip.GetRMS())
  hstrip.Fit("scfit","Q");
  dtmean=scfit.GetParameter(1)
  dtres=scfit.GetParameter(2)

  ntrg=heff.GetBinContent(2)
  n1=heff.GetBinContent(3)
  n2=heff.GetBinContent(4)
  eff=n1/ntrg
  deff=math.sqrt(eff*(1-eff)/ntrg)
  effp=n2/ntrg
  deffp=math.sqrt(effp*(1-effp)/ntrg)
  
  print "|%d|150|%d|%d|%d|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.3f|" % (run,int(ntrg),int(n1),int(n2),eff*100,deff*100,effp*100,deffp*100,dtmean,dtres)
  hstrip.Draw()
  c1.Update()
  c1.SaveAs("Run%d_Strip%d_pos.png" % (run,strip));

  val = raw_input()

def fitped(run,tdc,vthmin,vthmax,asic=1,ncha=24,rising=True,old=defped):
  rb=1
  fi=1
  la=ncha+1
  #if (asic==2):
  #    fi=ncha+1
  #    la=2*ncha+1
  asicmap={}
  for i in {1,2}:
      asicmap[i]=[]
      for j in range(49):
          asicmap[i].append(0)

  asicmap[1][1]=30
  asicmap[1][2]=28
  asicmap[1][3]=26
  asicmap[1][4]=24
  asicmap[1][5]=23
  asicmap[1][6]=22
  asicmap[1][7]=21
  asicmap[1][8]=20
  asicmap[1][9]=19
  asicmap[1][10]=18
  asicmap[1][11]=17
  asicmap[1][12]=16
  asicmap[1][13]=15
  asicmap[1][14]=14
  asicmap[1][15]=13
  asicmap[1][16]=12
  asicmap[1][17]=11
  asicmap[1][18]=10
  asicmap[1][19]=9
  asicmap[1][20]=8
  asicmap[1][21]=7
  asicmap[1][22]=6
  asicmap[1][23]=5
  asicmap[1][24]=4
  asicmap[2][25]=30
  asicmap[2][26]=28
  asicmap[2][27]=26
  asicmap[2][28]=24
  asicmap[2][29]=23
  asicmap[2][30]=22
  asicmap[2][31]=21
  asicmap[2][32]=20
  asicmap[2][33]=19
  asicmap[2][34]=18
  asicmap[2][35]=17
  asicmap[2][36]=16
  asicmap[2][37]=15
  asicmap[2][38]=14
  asicmap[2][39]=13
  asicmap[2][40]=12
  asicmap[2][41]=11
  asicmap[2][42]=10
  asicmap[2][43]=9
  asicmap[2][44]=8
  asicmap[2][45]=7
  asicmap[2][46]=6
  asicmap[2][47]=5
  asicmap[2][48]=4


  asicmap={}
  for i in {1,2}:
    asicmap[i]=[]
    for j in range(50):
      asicmap[i].append(0)
  f=open("/opt/TdcAnalysis/feb_mapping.json")
  s=json.loads(f.read())
  prh=s["v1_56"]["FlexTop"]["High"]["PR"]
  prl=s["v1_56"]["FlexTop"]["Low"]["PR"]
  tdch0=s["v1_56"]["FlexTop"]["High"]["TDC"][0]
  tdch1=s["v1_56"]["FlexTop"]["High"]["TDC"][1]
  tdcl0=s["v1_56"]["FlexTop"]["Low"]["TDC"][0]
  tdcl1=s["v1_56"]["FlexTop"]["Low"]["TDC"][1]
  for i in range(12):
      asicmap[1][tdch0[i]]=prh[i]
      asicmap[1][tdcl0[i]]=prl[i]
      asicmap[2][tdcl1[i]]=prl[i]
      asicmap[2][tdch1[i]]=prh[i]
  print asicmap
  ped=[]
  for i in range(32):
    ped.append(0)
  f82=TFile("Histos/InTime/histo%d_0.root" % run);
  f82.cd("/run%d/TDC%d" % (run,tdc));
  c1=TCanvas();
  #c2=TCanvas("c2","Test",1400,900);
  #c2.cd()
  #c2.Divide(6,4)
  #c2.Draw()
  #c2.Update()
  #val = raw_input()
  #c2.Draw()
  fout=open("summary_pedestal_%d_tdc%d.txt" % (run,tdc),"w");
  fout.write("+--+-----+-----+-----+ \n");
  gStyle.SetOptFit();
  hmean=TH1F("hmean","Summary %d %d " %(run,tdc),vthmax-vthmin+1,vthmin,vthmax)
  hnoise=TH1F("hnoise","Summary noise %d %d " %(run,tdc),100,0.,30.)
  hpmean=TH1F("hpmean","Summary %d %d " %(run,tdc),2*ncha,0.,2.*ncha);
  hpnoise=TH1F("hpnoise","Summary noise %d %d " %(run,tdc),2*ncha,0.,2.*ncha);
  scfit=TF1("scfit","[0]*TMath::Erfc((x-[1])/[2])",vthmin+1,vthmax);
  
  for ip in range(la-1,fi,-1):
      #c2.cd()
      if (asicmap[asic][ip]==0):
          continue;
      hs=None
      if (rising):
          hs=f82.Get("/run%d/TDC%d/vthc%d" % (run,tdc,ip));
      else:
          hs=f82.Get("/run%d/TDC%d/vthd%d" % (run,tdc,ip));
      if (hs==None):
          continue;
      if (hs.GetEntries()==0):
        continue
      print ip,fi,la," found"
      hs.Scale(1./2700.);
      nmax=0
      for i in range(1,hs.GetNbinsX()):
        if (hs.GetBinContent(i)==0):
              if (hs.GetBinContent(i-1)!=0 and hs.GetBinContent(i+1)!=0):
                  hs.SetBinContent(i,(hs.GetBinContent(i-1)+hs.GetBinContent(i+1))/2.)
        else:
          if (hs.GetBinContent(i)>nmax):
              nmax=hs.GetBinContent(i)

      hs.GetXaxis().SetRangeUser(vthmin-1,vthmax);
      icolor= ip%4 +1
      istyle= ip/4+1
      hs.SetLineColor(icolor)
      hs.SetLineStyle(istyle)
      hs.SetLineWidth(2)
      c1.cd()
      c1.Draw()
      



      if (ip==0):
        hs.Draw()
      else:
        hs.Draw("SAME")
  c1.Update()
  c1.SaveAs("Run%d_AllStrip%d_%d.root" % (run,tdc,asic));
  c1.SaveAs("Run%d_AllStrip%d_%d.png" % (run,tdc,asic));

  val = raw_input()
      
  for ip in range(fi,la):
      #c2.cd()
      if (asicmap[asic][ip]==0):
          continue;

      hs=None
      if (rising):
          hs=f82.Get("/run%d/TDC%d/vthc%d" % (run,tdc,ip));
      else:
          hs=f82.Get("/run%d/TDC%d/vthd%d" % (run,tdc,ip));
      if (hs==None):
          continue;
      if (hs.GetEntries()==0):
        continue
      #hs.Scale(1./2700.);
      hder=TH1F("hder%d" % ip,"derivative",900/rb,0.,900.)	
      hs.Rebin(rb)
      nmax=0
      for i in range(1,hs.GetNbinsX()):
          if (hs.GetBinContent(i)==0):
              if (hs.GetBinContent(i-1)!=0 and hs.GetBinContent(i+1)!=0):
                  hs.SetBinContent(i,(hs.GetBinContent(i-1)+hs.GetBinContent(i+1))/2.)
          else:
            if (hs.GetBinContent(i)>nmax):
              nmax=hs.GetBinContent(i)
      for i in range(1,hs.GetNbinsX()):
        if (hs.GetBinContent(i)-hs.GetBinContent(i+1)>-10):
          hder.SetBinContent(i,hs.GetBinContent(i)-hs.GetBinContent(i+1))
      hder.Rebin(2)
      hder.GetXaxis().SetRangeUser(vthmin-1,vthmax);
      scfit.SetParameter(0,nmax/2.);
      scfit.SetParameter(1,hder.GetMean());
      scfit.SetParameter(2,hder.GetRMS());



      hs.GetXaxis().SetRangeUser(vthmin-1,vthmax);
      hs.Fit("scfit","Q","",vthmin+2,vthmax);
      #hs.GetXaxis().SetRangeUser(vthmin-1,scfit.GetParameter(1)+60);
      #gPad.SetLogy();
      rped=scfit.GetParameter(1)
      c1.cd()
      c1.Draw()
      

      hder.Draw()
      c1.Update()
      val1 = raw_input()

      print "heho ",val1,rped,hder.GetMean()
      rped=hder.GetMean()
      if (len(val1)>0):
          rped=float(val1)
      hs.Draw()
      

      c1.cd()
      c1.Draw()
      c1.Update()

      fout.write("|%2d|%5.1f|%5.1f|%5.2f| \n" % (ip,scfit.GetParameter(0),rped,scfit.GetParameter(2)));
      ipr=0
      if (ip%2==1):
        ipr=ip/2
      else:
        ipr=31-ip/2
      firmwaret=[31,29,27,25,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6]
      #firmwaco=[31,29,27,25,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4]
      firmware2=[24,5,3,1,0,2,4,6,7,8,9,10,26,28,30,31,29,27,25,23,22,21,20,19]
      
      firmwareta1=[21,20,23,22,25,24,27,26,29,28,31,30,1,0,3,2,5,4,7,6,10,8,15,12]
      firmwareta2=[21,20,23,22,25,24,27,26,29,28,31,30,1,0,3,2,5,4,7,6,10,8,14,12]
      firmwaretb=[30,26,28,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,30,26,28,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4]

      if (asic==1):
          firmware=firmwareta1
      else:
          firmware=firmwareta2
      firmware=firmwaretb
      #if (ip>0):
      #  ipr=firmware[ip-fi-1]
      #else:
      #  ipr=0
      print ip,fi,ip-fi
      ipr=firmware[ip-fi]
      ipr=asicmap[asic][ip]
      ped[ipr]=rped
      print ip,ipr,rped,scfit.GetParameter(2)
      hmean.Fill(rped)
      hnoise.Fill(scfit.GetParameter(2))
      hpmean.SetBinContent(ip+1,rped);
      hpnoise.SetBinContent(ip+1,scfit.GetParameter(2))
      #c1.SaveAs("Run%d_Strip%d.root" % (run,ip));
      val = raw_input()

      #hder.Draw()
      
      #c1.Update()
      #val = raw_input()
  c1.cd()
  hmean.Draw()
  hpmean.GetYaxis().SetRangeUser(vthmin,vthmax)
  hpmean.Draw()
  c1.Update()
  c1.SaveAs("Summary_%d_TDC%d.png" % (run,tdc));
  val = raw_input()
  hnoise.Draw()
  c1.Update()
  val = raw_input()
  hpnoise.Draw()
  c1.Update()
  val = raw_input()
  c1.Update()
  c1.SaveAs("Summary_Noise_%d_TDC%d.png" % (run,tdc));

  fout.write("+--+-----+-----+-----+ \n");
  fout.close()
  print ped
  val = raw_input()
  med=5550.0
  for i in range(32):
    if (ped[i]==0):
      continue;
    if (ped[i] < med):
      #print med,ped[i]
      med=ped[i]

  med=med+5
  med=480
  print "Alignment to :",med
  dac=ped
  for i in range(32):
    if (ped[i]==0):
      continue;
    old[i]=0
    dac[i]=int(round(old[i]+(med-ped[i])*1./2.97))
  print "cor%d_%d=" % (tdc,asic),dac
  return dac

def calcped(oldpr,ped,median):
  print "Alignment to :",median
  dac=[]
  for i in range(32):
    dac.append(0)

  for i in range(32):
    if (ped[i]<=0):
      continue;
    dac[i]=int(round(oldpr[i]+(median-ped[i])*1./2.97))
  print dac
  return dac
import os
def process(runs,proc=True):
  if (proc):  
    for run in runs:
      os.system("./bin/tdcr -ggifpp_geom.json -r%d " % run)

  v=6800
  for run in runs:
    calcefn(run,1,v)
    v=v+200

  v=6800
  for run in runs:
    calcefn(run,2,v)
    v=v+200

def proclist(first,last,proc=True,vf=6700,step=100):
  if (proc):  
    for run in range(first,last+1):
      os.system("./bin/tdcr -ggifpp_geom.json -r%d " % run)

  v=vf
  for run in range(first,last+1):
    calcefn(run,1,v)
    v=v+step

  v=vf
  for run in range(first,last+1):
    calcefn(run,2,v)
    v=v+step
def processDCS(fdb,webdcs,proc=True,diro="."):
    conn = sqlite3.connect(fdb)
    conn.text_factory = str
    curs = conn.cursor()
    curs.execute("SELECT RUN,HV FROM runs WHERE DCS=%d" % webdcs)
    v=curs.fetchall()
    #print v[0]
    #print v[0][0]
    #return
    if (proc):  
        for x in v:
            os.system("./bin/tdcr -ggifpp_geom.json -r%d " % x[0])
            
    res=[]        
    for x in v:
        res.append(extractEfficiency(x[0],1,x[1],dirp=diro))
    for x in v:
        res.append(extractEfficiency(x[0],2,x[1],dirp=diro))
    return res

def storeResults(fdbi,fdbo,webdcs,histod=".",runtype=1):
    conn = sqlite3.connect(fdbo)
    conn.text_factory = str

    sql_create_res_table = """CREATE TABLE IF NOT EXISTS corana (
                                    ID integer PRIMARY KEY,
                                    RUN INTEGER NOT NULL,
                                    CHAMBER INTEGER NOT NULL,
                                    TYPE INTEGER NOT NULL,
                                    HV REAL,
                                    NCEVT REAL,
                                    EFFC REAL,
                                    DEFFC REAL,
                                    CSIZE REAL,
                                    NCLUS REAL,
                                    DAQFEBRATE REAL,
                                    DAQEFFLOSS REAL,
                                    XYFEBRATE REAL,
                                    EFFBACK REAL,
                                    EFFCOR  REAL

                                );"""

#febrate/febs[chamber],-febrate*10E-7,febratesel/febs[chamber])
    try:
        c = conn.cursor()
        c.execute(sql_create_res_table)
        conn.commit()
    except Error as e:
        print(e)






    
    curs = conn.cursor()
    res=processDCS(fdbi,webdcs,False,diro=histod)
    #run,chamber,hv,int(ntrg),int(nall),int(nxy),eff*100,deff*100,effp*100,deffp*100,mul,hrate.GetMean(),ncev,effc*100,deffc*100,csize,ncl]
    for r in res:
        #print len(r)
        sql_ins='''INSERT INTO corana(run,chamber,type,hv,ncevt,effc,deffc,csize,nclus,daqfebrate,daqeffloss,xyfebrate,effback,effcor) VALUES(%d,%d,%d,%d,%5.1f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f)''' % (r[0],r[1],runtype,r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9],r[10],r[11],r[12])
        print sql_ins
        curs.execute(sql_ins)
        conn.commit()
    
def buildTGraph(title,vx,dvx,vy,dvy,tx,ty):
    r_ti, r_p = array( 'd' ), array( 'd' )
    er_ti, er_p = array( 'd' ), array( 'd' )
    for x in vx:
      r_ti.append(x)
    for y in vy:
      r_p.append(y)
    for x in dvx:
      er_ti.append(x)
    for y in dvy:
      er_p.append(y)

    gr = TGraphErrors( len(vx), r_ti, r_p,er_ti,er_p )
    gr.SetLineColor( 1 )
    gr.SetLineWidth( 1 )
    gr.SetMarkerColor( 2 )
    gr.SetMarkerStyle( 21 )
    gr.SetMarkerSize( 0.4 )
    gr.SetTitle(title)
    gr.GetXaxis().SetTitle(tx )
    gr.GetYaxis().SetTitle( ty )
    return gr

def buildTGraph1(title,vx,vy,tx,ty):
    r_ti, r_p = array( 'd' ), array( 'd' )
    er_ti, er_p = array( 'd' ), array( 'd' )
    for x in vx:
      r_ti.append(x)
    for y in vy:
      r_p.append(y)
   
    
    gr = TGraph( len(vx), r_ti, r_p)
    gr.SetLineColor( 1 )
    gr.SetLineWidth( 1 )
    gr.SetMarkerColor( 2 )
    gr.SetMarkerStyle( 21 )
    gr.SetMarkerSize( 0.4 )
    gr.SetTitle(title)
    gr.GetXaxis().SetTitle(tx )
    gr.GetYaxis().SetTitle( ty )
    return gr
def  drawDCS(fdbi,webdcs,chamber,c=None):
    if (c==None):
        c=TCanvas()
    jdict={}
    jdict['settings']={}
    gStyle.SetOptFit(1)
    c.Clear()
    #c.Divide(2,4)
    conn = sqlite3.connect(fdbi)
    conn.text_factory = str
    curs = conn.cursor()
    sql_dcs="select ATT,DEAD,TRET,TCOAX,TRIGGER,ACTIVE,START from webdcs WHERE dcs=%d" % webdcs
    curs.execute(sql_dcs)
    v=curs.fetchall()
    att=-1
    dthr=-1
    dead=-1
    trig="UNKNOwN"
    if (len(v)<1):
        return
    jdict['settings']['dcs']=webdcs
    jdict['settings']['att']=v[0][0]
    jdict['settings']['dead']=v[0][1]
    jdict['settings']['tret']=v[0][2]
    jdict['settings']['tcoax']=v[0][3]
    jdict['settings']['trigger']=v[0][4]
    jdict['settings']['active']=v[0][5]
    jdict['settings']['start']=v[0][6]

    dirout="./results/dcs/%d_ATT%3.1f_DT%d_THR%d_%s/chamber%d/" % (webdcs,v[0][0],v[0][1],v[0][2]-500,v[0][4].replace('"','').replace(' ','_'),chamber)
    os.system("mkdir -p %s" % dirout)
    fout=open(dirout+"summary%d_ATT%3.1f_DT%d_THR%d_%s.txt" % (webdcs,v[0][0],v[0][1],v[0][2]-500,v[0][4].replace('"','').replace(' ','_')),"w")
    for x in v:
        fout.write("Cuts: %5.2f %d %d %d %s %d %s\n" % x)
        att=x[0]
        dead=x[1]
        dthr=x[2]-500
        trig=x[4]
    sql_query=" select EFFCOR,DEFFC,(SELECT HV FROM runs WHERE runs.RUN=corana.RUN),DAQFEBRATE,DAQEFFLOSS,RUN,NCEVT,EFFC,CSIZE,NCLUS,XYFEBRATE,EFFBACK  from corana WHERE RUN IN (SELECT RUN FROM runs WHERE DCS=%d) AND CHAMBER=%d" % (webdcs,chamber)
    curs.execute(sql_query)
    vo=curs.fetchall()
    if (len(vo)<1):
        return;
    hv=[]
    eff=[]
    dhv=[]
    deff=[]
    febrate=[]
    dfebrate=[]
    effloss=[]
    csize=[]
    nclus=[]
    effback=[]
    jdict['runs']={}
    jdict['runs']['effcor']=[]
    jdict['runs']['hv']=[]
    jdict['runs']['febrate']=[]
    jdict['runs']['effloss']=[]
    jdict['runs']['runid']=[]
    jdict['runs']['nevt']=[]
    jdict['runs']['effclu']=[]
    jdict['runs']['csize']=[]
    jdict['runs']['nclus']=[]
    jdict['runs']['xyrate']=[]
    jdict['runs']['effback']=[]

    
    for x in vo:
        if (x[1]==0):
            continue
        fout.write("Results: %5.2f %5.2f %5.2f %5.2f %5.2f %d %d %5.2f %5.2f %5.2f %5.2f %5.2f \n" % x)
        hv.append(x[2])
        dhv.append(10.)
        eff.append(x[0])
        deff.append(x[1])
        febrate.append(x[3])
        effloss.append(x[4])
        csize.append(x[8])
        nclus.append(x[9])
        effback.append(x[11])
        jdict['runs']['effcor'].append(x[0])
        jdict['runs']['hv'].append(x[2])
        jdict['runs']['febrate'].append(x[3])
        jdict['runs']['effloss'].append(x[4])
        jdict['runs']['runid'].append(x[5])
        jdict['runs']['nevt'].append(x[6])
        jdict['runs']['effclu'].append(x[7])
        jdict['runs']['csize'].append(x[8])
        jdict['runs']['nclus'].append(x[9])
        jdict['runs']['xyrate'].append(x[10])
        jdict['runs']['effback'].append(x[11])
    stitle="DCS%d_TRG%s_ATT%3.1f_THR%d_DT%d_CH%d" % (webdcs,trig,att,dthr,dead,chamber)
    gr=buildTGraph("effi",hv,dhv,eff,deff,"HV effective (V)","efficiency (%)")

    func = TF1("func", "([0]/(1+ TMath::Exp(-[1]*(x-[2]))))", 6500,8200)
    func.SetParameters(90, 9.E-3, 7000)
    print 100, 9.E-3, 7100
    gr.Fit(func,"","",6700,8200)
    hv95=func.GetX(func.GetParameter(0)*0.95)
    hv99=func.GetX(func.GetParameter(0)*0.99)
    print "HV95",hv95,hv95+150,hv99,hv99-hv95
    wp=hv99
    jdict['fit']={}
    jdict['fit']['Efficiency']=func.GetParameter(0)
    jdict['fit']['Slope']=func.GetParameter(1)
    jdict['fit']['HV50']=func.GetParameter(2)
    jdict['fit']['HV95']=hv95
    jdict['fit']['HV99']=wp

    fout.write("FIT results: %f %f %f %f %f \n" % (func.GetParameter(0),func.GetParameter(1),func.GetParameter(2),hv95,wp))
    title=stitle+"_HV95_%4.0f_WP_%4.0f" % (hv95,wp)
    gr.SetTitle(title)
    gStyle.SetStatX(0.85)
    gStyle.SetStatY(0.7)
    c.cd(1)
    gr.Draw("AP")
    c.Update()
    #val=raw_input()
    c.SaveAs(dirout+"%s.png" % title)

    tgr=[]
    tgr.append(buildTGraph1('FEB Rate',hv,febrate,'HV eff (V)','Rate (Hz/cm^2)'))
    tgr.append(buildTGraph1('Dead Time Loss',hv,effloss,'HV eff (V)','Dead time loss (%)'))
    tgr.append(buildTGraph1('Cluster Size',hv,csize,'HV eff (V)','Cluster size'))
    tgr.append(buildTGraph1('Cluster Number',hv,nclus,'HV eff (V)','Clusters'))
    tgr.append(buildTGraph1('Background Efficiency',hv,effback,'HV eff (V)','Efficiency (%)'))
    #val=raw_input()
    # grb=buildTGraph("background",hv,dhv,febrate,dfebrate,"HV effective (V)","FEB rate (Hz/cm^2)")
    # title=stitle
    # grb.SetTitle(title)
    # gStyle.SetStatX(0.85)
    # gStyle.SetStatY(0.7)
    # c.cd(2)
    # grb.Draw("AP")
    # c.Update()
    # c.SaveAs(dirout+"%s.png" % title)
    #val=raw_input()
    sql_query="select * from RESULTS WHERE  CTIME>=(SELECT CFIRST FROM webdcs WHERE DCS=%d) AND CTIME<=(SELECT CLAST FROM webdcs WHERE DCS=%d)-200 AND HARDWARE='BMP'" % (webdcs,webdcs)
    curs.execute(sql_query)
    v=curs.fetchall()
    bmp=[]
    for x in v:
        bmp.append([x[3],json.loads(x[4].decode('latin-1').encode("utf-8"))])
    #print bmp
    PB=0
    TB=0
    NB=0
    for x in bmp:
      PB=PB+x[1]['pressure']
      TB=TB+x[1]['temperature']
      NB=NB+1
    if (NB>0):
      PB=PB/NB
      TB=TB/NB+273.15
    fout.write("Pressure: %d %5.2f %5.2f \n" % (NB,PB,TB))
    jdict['BMP']={"P":PB,"T":TB}

    #val=raw_input()
    sql_query="select * from RESULTS WHERE  CTIME>=(SELECT CFIRST FROM webdcs WHERE DCS=%d) AND CTIME<=(SELECT CLAST FROM webdcs WHERE DCS=%d)-200 AND HARDWARE='SY1527'" % (webdcs,webdcs)
    curs.execute(sql_query)
    v=curs.fetchall()
    
    a=[]
    for x in v:
        a.append([x[3],json.loads(x[4].decode('latin-1').encode("utf-8"))])
    #print a
    #val=raw_input()
    if (chamber==1):
      chan=[4,5]
    else:
      chan=[1,2] 
    idx=0
    #tgr=[]
    jdict['HV']={}
    for ch in chan:
          x_t=[]
          y_vs=[]
          z_vm=[]
          w_im=[]
          g_g=[]
          chname=""
          first=0
          for x in a:
            if (x[1]['channels'][ch]['rampup']==0):
              continue
            if (float(x[1]['channels'][ch]['vset']) < 5500.0):
              continue
            #print float(x[1]['channels'][ch]['vset'])
            if (first==0):
              first=x[0]
            x_t.append(x[0]-first)
            y_vs.append(x[1]['channels'][ch]['vset'])
            z_vm.append(x[1]['channels'][ch]['vout'])
            w_im.append(x[1]['channels'][ch]['iout'])
            g_g.append( x[1]['channels'][ch]['iout']/ x[1]['channels'][ch]['vout'])
            #print ch,x
            chname=x[1]['channels'][ch]['name']
          if (len(x_t)<1):
            continue
          dy_vs=[]
          if (len(y_vs)>1):
            for i in range(0,len(y_vs)):
              if (i==0 and i!=len(y_vs)-1):
                dy_vs.append(y_vs[i+1]-y_vs[i])
                continue
              if (i==len(y_vs)-1 and i!=0):
                dy_vs.append(y_vs[i]-y_vs[i-1])
                continue
              dy_vs.append((y_vs[i+1]-y_vs[i-1])/2.)
          else:
            dy_vs.append(0)
          vset=[]
          vmon=[]
          veff=[]
          imon=[]
          nval=0
          vmon_sum=0
          vset_sum=0
          imon_sum=0
          for i in range(0,len(y_vs)):
            if (abs(dy_vs[i])>10):
              if (nval>0):
                vmon.append(vmon_sum/nval)
                vset.append(vset_sum/nval)
                imon.append(imon_sum/nval)
                if (NB>0):
                  veff.append(calV(vmon_sum/nval,PB,TB))
                else:
                  veff.append(vmon_sum/nval)
              nval=0
              vmon_sum=0
              vset_sum=0
              imon_sum=0
              continue
            nval=nval+1
            vmon_sum=vmon_sum+z_vm[i]
            vset_sum=vset_sum+y_vs[i]
            imon_sum=imon_sum+w_im[i]
          if (nval>0):
            vmon.append(vmon_sum/nval)
            vset.append(vset_sum/nval)
            imon.append(imon_sum/nval)
            if (NB>0):
              veff.append(calV(vmon_sum/nval,PB,TB))
            else:
              veff.append(vmon_sum/nval)
          #print chname,"VSET",vset
          #print chname,"VMON",vmon
          #print chname,"VEFF",veff
          #print chname,"IMON",imon
          fout.write(chname+"\n")
          jdict['HV'][chname]={}
          jdict['HV'][chname]['vset']=vset
          jdict['HV'][chname]['vmon']=vmon
          jdict['HV'][chname]['imon']=imon
          jdict['HV'][chname]['veff']=veff

          for ip in range(0,len(vmon)):
            fout.write("%5.2f %5.2f %5.2f %5.2f \n" % (vset[ip],vmon[ip],veff[ip],imon[ip]))
          #tgr.append(buildTGraph1('V set vs t  %s' % chname,x_t,y_vs,'t(s)','V set (V)'))
          tgr.append(buildTGraph1('V Mon vs t  %s' % chname,x_t,z_vm,'t(s)','V mon (V)'))
          tgr.append(buildTGraph1('I Mon vs t  %s' % chname,x_t,w_im,'t(s)','I mon ([m]A)'))
          tgr.append(buildTGraph1('I Mon vs V eff %s' % chname,veff,imon,'V eff(V)','I mon ([m]A)'))
    icd=3
    for x in tgr:
      c.cd(icd)
      icd=icd+1
      x.Draw("AP")
      c.Update()
      c.SaveAs(dirout+"%s.png" % x.GetTitle().replace(" ","_"))
      #val=raw_input()      
    # resume
    chans=[]
    if (chamber == 2):
      chans=["COAX-BOT","COAX-TOP"]
    else:
      chans=["RETURN-BOT","RETURN-TOP"]
    jdict['AWP']={}
    jdict['AWP']['hv']=wp
    ihv=-1
    for i in range(0,len(jdict['runs']['hv'])-1):
      if (wp>jdict['runs']['hv'][i] and wp<=jdict['runs']['hv'][i+1]):
        ihv=i
        break
    if (ihv!=-1):
      jdict['AWP']['plateau']=jdict['fit']['Efficiency']
      jdict['AWP']['febrate']=approx(wp,jdict['runs']['hv'][ihv],jdict['runs']['hv'][ihv+1],jdict['runs']['febrate'][ihv],jdict['runs']['febrate'][ihv+1])
      jdict['AWP']['effloss']=approx(wp,jdict['runs']['hv'][ihv],jdict['runs']['hv'][ihv+1],jdict['runs']['effloss'][ihv],jdict['runs']['effloss'][ihv+1])
      jdict['AWP']['effcor']=approx(wp,jdict['runs']['hv'][ihv],jdict['runs']['hv'][ihv+1],jdict['runs']['effcor'][ihv],jdict['runs']['effcor'][ihv+1])
      jdict['AWP']['csize']=approx(wp,jdict['runs']['hv'][ihv],jdict['runs']['hv'][ihv+1],jdict['runs']['csize'][ihv],jdict['runs']['csize'][ihv+1])
      jdict['AWP']['nclus']=approx(wp,jdict['runs']['hv'][ihv],jdict['runs']['hv'][ihv+1],jdict['runs']['nclus'][ihv],jdict['runs']['nclus'][ihv+1])
    for x in chans:
      ihv=-1
      if (x in jdict['HV']):
        for i in range(0,len(jdict['HV'][x]['veff'])-1):
          if (wp>jdict['HV'][x]['veff'][i] and wp<=jdict['HV'][x]['veff'][i+1]):
            ihv=i
            break
      if (ihv!=-1):
         jdict['AWP'][x]=approx(wp,jdict['HV'][x]['veff'][ihv],jdict['HV'][x]['veff'][ihv+1],jdict['HV'][x]['imon'][ihv],jdict['HV'][x]['imon'][ihv+1])
    itot=0
    surf=0
    if (chamber==1):
      surf=13000
      if ('RETURN-TOP' in jdict['AWP']):
        itot=itot+jdict['AWP']['RETURN-TOP']
      if ('RETURN-BOT' in jdict['AWP']):
        itot=itot+jdict['AWP']['RETURN-BOT']
    if (chamber==2):
      surf=15000
      if ('COAX-TOP' in jdict['AWP']):
        itot=itot+jdict['AWP']['COAX-TOP']
      if ('COAX-BOT' in jdict['AWP']):
        itot=itot+jdict['AWP']['COAX-BOT']
    if ('febrate' in jdict['AWP']):
      if (jdict['AWP']['febrate']!=0):
        jdict['AWP']['ITOT']=itot
        jdict['AWP']['QSEEN']=itot/jdict['AWP']['febrate']/surf*1E6
        fout.write("%d|%5.1f|%5.0f|%5.2f|%5.2f|%5.2f|%5.2f|%5.1f|%5.1f|%5.1f|%5.1f\n" % 
        (webdcs,
        jdict['settings']['att'],
        jdict['AWP']['hv'],
        jdict['AWP']['plateau'],
        jdict['AWP']['effcor'],
        jdict['AWP']['effloss'],
        jdict['AWP']['plateau']-jdict['AWP']['effloss'],
        jdict['AWP']['febrate'],
        jdict['AWP']['ITOT'],
        jdict['AWP']['QSEEN'],jdict['AWP']['csize']))
    with open(dirout+'summary.json', 'w') as outfile:
      dd= json.dumps(jdict,sort_keys=True, indent=2,separators=(',', ': '))
      outfile.write(dd)
      outfile.close()
    fout.close()
    return jdict


def processAllDCS(fdb,dbo,diro=".",proc=True,store=True,draw=False,canvas=None,first=0,last=100000):
    fout=open('resume.txt','w')
    conn = sqlite3.connect(fdb)
    conn.text_factory = str
    curs = conn.cursor()
    curs.execute("select distinct(dcs) from webdcs where dcs>=%d and dcs<=%d" % (first,last))
    v=curs.fetchall()
    for x in v:
        if (proc):
            processDCS(fdb,x[0],proc,diro)
        if (store):
          storeResults(fdb,dbo,x[0],diro)
        if (draw):
          jdict=drawDCS(fdb,x[0],1,canvas)
          if ('QSEEN' in jdict['AWP']):
            fout.write("%d|%d|%s|%5.1f|%5.1f|%5.0f|%5.2f|%5.2f|%5.2f|%5.2f|%5.1f|%5.1f|%5.1f\n" % 
            (jdict['settings']['dcs'],1,
            jdict['settings']['trigger'],
            jdict['AWP']['febrate'],
            jdict['settings']['att'],
            jdict['AWP']['hv'],
            jdict['AWP']['plateau'],
            jdict['AWP']['effcor'],
            jdict['AWP']['effloss'],
            jdict['AWP']['plateau']-jdict['AWP']['effloss'],
            jdict['AWP']['ITOT'],
            jdict['AWP']['QSEEN'],jdict['AWP']['csize']))
          jdict=drawDCS(fdb,x[0],2,canvas)
          if ('QSEEN' in jdict['AWP']):
            fout.write("%d|%d|%s|%5.1f|%5.1f|%5.0f|%5.2f|%5.2f|%5.2f|%5.2f|%5.1f|%5.1f|%5.1f\n" % 
            (jdict['settings']['dcs'],2,
            jdict['settings']['trigger'],
            jdict['AWP']['febrate'],
            jdict['settings']['att'],
            jdict['AWP']['hv'],
            jdict['AWP']['plateau'],
            jdict['AWP']['effcor'],
            jdict['AWP']['effloss'],
            jdict['AWP']['plateau']-jdict['AWP']['effloss'],
            jdict['AWP']['ITOT'],
            jdict['AWP']['QSEEN'],jdict['AWP']['csize']))
    fout.close()
def approx(v,v1,v2,a1,a2):
  return a1+(v-v1)*(a2-a1)*1./(v2-v1)

def pltres(run,direc="/data/NAS/EM888/Current"):

    _file0 = TFile("%s/histo%d_0.root" % (direc,run));
    c=TCanvas("c","Strip studies %d" % run ,545,345);
    _file0.cd("/FEB/Chamber1/Strips/");
    
    XY = _file0.Get("/FEB/Chamber1/Strips/XY");									       
    XY.SetTitle("Strip hits position in the time window (7400 V Thr 525)");
    XY.GetXaxis().SetRangeUser(0.,35.);
    XY.GetYaxis().SetTitle("Y_{Strip} (cm)");
    XY.GetXaxis().SetTitle("X_{Strip} (cm)");
    XY.Draw();
    c.SaveAs("Run%d_StripXYAll.png" % run);
    XYMin = _file0.Get("/FEB/Chamber1/Strips/XYMin");									       
    XYMin.SetTitle("Strip hits position nearest to the extrapolation (7400 V Thr 525)");
    XYMin.GetXaxis().SetRangeUser(0.,35.);
    XYMin.GetYaxis().SetTitle("Y_{Strip} (cm)");
    XYMin.GetXaxis().SetTitle("X_{Strip} (cm)");
    XYMin.Draw();
    c.SaveAs("Run%d_StripXYMin.png" % run);
    XYSel = _file0.Get("/FEB/Chamber1/Strips/XYSel");
    XYSel.SetTitle("Strip hit associated (7400 V Thr 525)");
    XYSel.GetXaxis().SetRangeUser(0.,35.);
    XYSel.GetYaxis().SetTitle("Y_{Strip} (cm)");
    XYSel.GetXaxis().SetTitle("X_{Strip} (cm)");
    XYSel.Draw();
    c.SaveAs("Run%d_StripXYSelected.png" % run);

    DXDY = _file0.Get("/FEB/Chamber1/Strips/DXDY");
    DXDY.SetTitle("Distance of nearest hit to the extrapolation (7400 V Thr 525)");
    DXDY.GetYaxis().SetTitle("#delta Y_{Strip} (cm)");
    DXDY.GetXaxis().SetTitle("#delta X_{Strip} (cm)");
    DXDY.Draw();
    c.SaveAs("Run%d_StripDXDYSelected.png" % run);

    AbsoluteTimeDelta= _file0.Get("/FEB/Chamber1/Strips/AbsoluteTimeDelta");
    AbsoluteTimeDelta.SetTitle("Other hits, Absolute Time distance to the selected hit time (7400 V Thr 525)");
    AbsoluteTimeDelta.GetXaxis().SetTitle("#delta T_{Absolute} (ns)");
    AbsoluteTimeDelta.GetYaxis().SetTitle("Number of hits");
    AbsoluteTimeDelta.GetXaxis().SetRangeUser(-300.,300.);
    AbsoluteTimeDelta.Draw();
    c.SaveAs("Run%d_StripAbsoluteDeltaT.png" % run);
    AbsoluteTimeDelta.GetXaxis().SetRangeUser(-30.,130.);
    AbsoluteTimeDelta.Draw();
    c.SaveAs("Run%d_StripAbsoluteDeltaTZoomed.png" % run);

    _file0.cd("/FEB/Chamber1/ClusterNew/");
  
    XYC = _file0.Get("/FEB/Chamber1/ClusterNew/XY");									       
    XYC.SetTitle("Cluster position in the time window (7400 V Thr 525)");
    XYC.GetXaxis().SetRangeUser(0.,35.);
    XYC.GetYaxis().SetTitle("Y_{Strip} (cm)");
    XYC.GetXaxis().SetTitle("X_{Strip} (cm)");
    XYC.Draw();
    c.SaveAs("Run%d_ClusterXYAll.png" % run);
    
    XYMinC = _file0.Get("/FEB/Chamber1/ClusterNew/XYMax");									       
    XYMinC.SetTitle("Cluster position nearest to the extrapolation (7400 V Thr 525)");
    XYMinC.GetXaxis().SetRangeUser(0.,35.);
    XYMinC.GetYaxis().SetTitle("Y_{Strip} (cm)");
    XYMinC.GetXaxis().SetTitle("X_{Strip} (cm)");
    XYMinC.Draw();
    c.SaveAs("Run%d_ClusterXYMin.png" % run);
    
    XYSelC = _file0.Get("/FEB/Chamber1/ClusterNew/XYSel");
    XYSelC.SetTitle("Cluster associated (7400 V Thr 525)");
    XYSelC.GetXaxis().SetRangeUser(0.,35.);
    XYSelC.GetYaxis().SetTitle("Y_{Strip} (cm)");
    XYSelC.GetXaxis().SetTitle("X_{Strip} (cm)");
    XYSelC.Draw();
    c.SaveAs("Run%d_ClusterXYSelected.png" % run);

    DXDYC = _file0.Get("/FEB/Chamber1/ClusterNew/DIST");
    DXDYC.SetTitle("Distance of nearest cluster to the extrapolation (7400 V Thr 525)");
    DXDYC.GetYaxis().SetTitle("#delta Y_{Strip} (cm)");
    DXDYC.GetXaxis().SetTitle("#delta X_{Strip} (cm)");
    DXDYC.Draw();
    c.SaveAs("Run%d_ClusterDXDYSelected.png" % run);

    hnclus= _file0.Get("/FEB/Chamber1/ClusterNew/Clusters");
    hnclus.SetTitle("Number of clusters (7400 V Thr 525)");
    hnclus.GetXaxis().SetTitle("Number of clusters reconstructed in the event");
    hnclus.GetYaxis().SetTitle("Number of events");
    hnclus.Draw();
    c.SaveAs("Run%d_ClusterNclus.png" % run);
    
    cluss= _file0.Get("/FEB/Chamber1/ClusterNew/ClusterSize");
    cluss.SetTitle("All clusters size (7400 V Thr 525)");
    cluss.GetXaxis().SetTitle("Number of strip hits in cluster");
    cluss.GetYaxis().SetTitle("Number of events");
    cluss.Draw();
    c.SaveAs("Run%d_ClusterClusterSize.png" % run);
    
    cluss1= _file0.Get("/FEB/Chamber1/ClusterNew/ClusterSize1");
    cluss1.SetTitle("Associated cluster  size (7400 V Thr 525)");
    cluss1.GetXaxis().SetTitle("Number of strip hits in cluster");
    cluss1.GetYaxis().SetTitle("Number of events");
    cluss1.Draw();
    c.SaveAs("Run%d_ClusterClusterSizeAssociated.png" % run);
 



