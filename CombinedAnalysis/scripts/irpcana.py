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

class analyse:
    def __init__(self,run,jsfile,direc="/data/NAS/EM888/Current"):
        self.run = run
        self.direc = direc
        self.f0 = TFile(direc+"/histo%d_0.root" % run)
        self.f1 = TFile(direc+"/Noisehisto%d_0.root" % run)
        self.P=0
        self.P0=0
        self.HVapp=0
        self.HVeff=0
        self.T0=0
        self.T=0
        self.Threshold=0
        self.jsfile=jsfile
        f=open(jsfile)
        self.params=json.loads(f.read())
        
    def setConditions(self,hvapp,thr,p=0,t=0,p0=0,t0=0):
        self.HVapp=hvapp
        self.Threshold=thr
        if (p!=0):
            self.P=p
        if (t!=0):
            self.T=t
        if (t0!=0):
            self.T0=t0
            self.T=t0+273.15+2
        if (p0!=0):
            self.P0=p0
            self.P=self.calP(p0,170)
        if (self.P!=0 and self.T!=0):
            self.HVeff=self.calV(self.HVapp,self.P,self.T)
        self.comment="R%d (%d/%d) VTH %d(%.1f fC)" % (self.run,self.HVapp,self.HVeff,self.Threshold,(thr-480)*2.4)
    def calP(self,p0,alt):
        return p0*(1-0.0065*alt/288.15)**5.255

    def calV(self,V,P,T):
        return  V/(0.2+0.8*P/990.*293./T)

    def calValice(self,V,P,T):
        print 1-(1*P/990.*293./T)
        return  V/(1.0*P/990.*293./T)

    def calApp(self,V,P,T):
        print 1-(0.2+0.8*P/990.*293./T)
        return  V*(0.2+0.8*P/990.*293./T)

    def getdT(self):

        c=TCanvas("irpc","Alignement Studies",545,342)
        self.f0.cd("/Align")

        gStyle.SetOptStat(0);
        dT0=[]
        dT1=[]
        T0=[]
        T1=[]
        for i in range(48):
            dT0.append(0)
            T0.append(0)
            dT1.append(0)
            T1.append(0)
        for ist in range(2,33):
            hst0=self.f0.Get("/Align/Pedestal%dHR" % (ist))
            hst1=self.f0.Get("/Align/Pedestal%dLR" % (ist))
    
            if (hst0!=None):
                dT0[ist]=hst0.GetMean()
            if (hst1!=None):
                dT1[ist]=hst1.GetMean()
        print dT0
        print dT1
        t0=0
        t1=0
        for ist in range(1,48):
            T0[ist]=t0
            t0=t0+dT0[ist]
            T1[ist]=t1
            t1=t1+dT1[ist]
        print T0
        print T1
        ST0=[]
        ST1=[]
        ST0.append(0)
        ST1.append(0)
        
        str_res='"st0":[0.'
        for i in range(2,48):
            str_res=str_res+",%.2f" % T0[i]
            ST0.append(T0[i])
        str_res=str_res+"],\n"
        str_res+='"st1":[0.'
        for i in range(2,48):
            str_res=str_res+",%.2f" % T1[i]
            ST1.append(T1[i])
        str_res=str_res+"],\n"
        
        print str_res
        fout=open("DT_%d.json" % (self.run),"w");
        fout.write(str_res);
        fout.close()

        self.params["febs"][0]["st0"]=ST0
        self.params["febs"][1]["st0"]=ST0
        self.params["febs"][0]["st1"]=ST1
        self.params["febs"][1]["st1"]=ST1
        fo=open('DTalign_%d.json' % self.run, 'w')
        json.dump(self.params,fo)
        fo.close()
    def getdy(self,hn="strip"):

        c=TCanvas("irpc","Alignement Studies",545,342)
        self.f0.cd("/Align")
        gStyle.SetOptFit();
        gStyle.SetOptStat(0);
        res=[]
        dres=[]
        for i in range(48):
            res.append(0)
            dres.append(0)
        for ist in range(1,33):
            hst=self.f0.Get("/Align/%s%d" % (hn,ist))
    
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
            c.SaveAs("IRPC%d_dy_%s%d.png" % (self.run,hn,ist))
            
        val=raw_input()
        str_res='"delta":[0.'
        for i in range(1,48):
            str_res=str_res+",%.2f" % res[i]
        str_res=str_res+"],\n"
        print str_res
        fout=open("align_%d.json" % (self.run),"w");
        fout.write(str_res);
        fout.close()
        y=np.array(res)
        dy=np.array(dres)
        #np.set_printoptions(precision=1)
        np.set_printoptions(formatter={'float': '{: 0.2f}'.format})
        print y
        print dy
        self.params["febs"][0]["delta"]=res
        self.params["febs"][1]["delta"]=res
        fo=open('DYalign_%d.json' % self.run, 'w')
        json.dump(self.params,fo)
        fo.close()


    def getratio(self,chamber=1,rebin=1,ncut=30):
        c=TCanvas("irpc","IRPC Studies",545,842)
        c.Divide(2,4)
        self.f0.cd("/FEB/Chamber%d/Raw" % chamber)
        hst=self.f0.Get("/FEB/Chamber%d/Raw/Strips" % chamber)
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
        self.f0.cd("/gric")
        hxy=self.f0.Get("/gric/XY")
        hxy.SetTitle(" 4 points Track extrapolation to the iRPC")
        hxyf=self.f0.Get("/gric/XYF")
        hxyf14=self.f0.Get("/gric/XYF14")
        hxyf15=self.f0.Get("/gric/XYF15")
        hxyf.SetTitle(" 4 points Track extrapolation to the iRPC when cluster found nearby (5cm,20cm)")
        hxyt=self.f0.Get("/gric/XYT")
        result=[]
        efft=hxyt.GetEntries()/hxy.GetEntries()
        result.append(self.run)
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
        chxy=hxy.Clone("chxy")
        chxyt=hxyt.Clone("chxyt")
        chxyf=hxyf.Clone("chxyf")
        chxyf14=hxyf14.Clone("chxyf14")
        chxyf15=hxyf15.Clone("chxyf15")
        
        chxy.Rebin2D(rbx,rby)
        chxyf.Rebin2D(rbx,rby)
        chxyf14.Rebin2D(rbx,rby)
        chxyf15.Rebin2D(rbx,rby)
        chxyt.Rebin2D(rbx,rby)
        chxy.SetAxisRange(0.,35.,"X")
        chxy.SetAxisRange(20.,120.,"Y")
        chxyf.SetAxisRange(0.,35.,"X")
        chxyf.SetAxisRange(20.,120.,"Y")
        chxyf14.SetAxisRange(0.,35.,"X")
        chxyf14.SetAxisRange(20.,120.,"Y")
        chxyf15.SetAxisRange(0.,35.,"X")
        chxyf15.SetAxisRange(20.,120.,"Y")
        chxyEff=chxyf.Clone("chxyEff")
        chxyEff.Divide(chxy)
        chxyEff.SetTitle("Local Efficiency map")
        chxyEff14=chxyf14.Clone("chxyEff14")
        chxyEff14.Divide(chxy)
        chxyEff14.SetTitle("Local Efficiency map FR4")
        chxyEff15=chxyf15.Clone("chxyEff15")
        chxyEff15.Divide(chxy)
        chxyEff15.SetTitle("Local Efficiency map EM888")
        heff=TH1F("heff","Local Efficiency Ntk ext>15",110,0,1.1)
        dmin=max(0,effo*0.7)
        dmax=min(1.1,effo*1.6)
        nb=int((dmax-dmin)/0.005)
        heff1=TH1F("heff1","Local Efficiency EM888",nb,dmin,dmax)
        heff32=TH1F("heff32","Local Efficiency FR4",nb,dmin,dmax)

        for i in range(1,chxy.GetNbinsX()):
            for j in range(1,chxy.GetNbinsY()):
                if (chxy.GetBinContent(i,j)>ncut):
                    heff.Fill(chxyEff.GetBinContent(i,j))
                    xp=chxy.GetXaxis().GetBinCenter(i)
                    yp=chxy.GetYaxis().GetBinCenter(j)
                    # Jusqu'au run 1721
                    #if (xp<=13):
                    #    heff1.Fill(chxyEff14.GetBinContent(i,j))
                    #if (xp>=19):
                    #    heff32.Fill(chxyEff15.GetBinContent(i,j))
                    # A partir du 1722 (14 et 15 inverted)
                    if (xp>=8 and xp<=15 and yp<84 and yp>41):
                        heff1.Fill(chxyEff15.GetBinContent(i,j))
                        print chxy.GetBinContent(i,j), chxyf15.GetBinContent(i,j),chxyEff15.GetBinContent(i,j)
                    if (xp>=17 and xp<=26 and yp<84 and yp>41):
                        heff32.Fill(chxyEff14.GetBinContent(i,j))

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
        y.tofile("result%d.csv" % self.run,sep='|',format='%7.1f')
        print y
  
        c.cd(5)
        chxyEff15.Draw("COLZ")
        c.Modified()
        c.Update()
        #val=raw_input()
        c.cd(6)
        chxyEff14.Draw("COLZ")
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
        c.SaveAs("Analyse%d.pdf" % self.run)
        val=raw_input()

    def pltres(self):
        c=TCanvas("c","IRPC studies %d" % self.run ,545,345);
        self.f0.cd("/FEB/Chamber1/Strips/");
        
        XY = self.f0.Get("/FEB/Chamber1/Strips/XY");									       
        XY.SetTitle("Strip hits position in the time window %s" % self.comment);
        XY.GetXaxis().SetRangeUser(0.,35.);
        XY.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XY.GetXaxis().SetTitle("X_{Strip} (cm)");
        XY.Draw();
        c.SaveAs("Run%d_StripXYAll.png" % self.run);
        XYMin = self.f0.Get("/FEB/Chamber1/Strips/XYMin");									       
        XYMin.SetTitle("Strip hits position nearest to the extrapolation %s" % self.comment);
        XYMin.GetXaxis().SetRangeUser(0.,35.);
        XYMin.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYMin.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYMin.Draw();
        c.SaveAs("Run%d_StripXYMin.png" % self.run);
        XYSel = self.f0.Get("/FEB/Chamber1/Strips/XYSel");
        XYSel.SetTitle("Strip hit associated %s" % self.comment);
        XYSel.GetXaxis().SetRangeUser(0.,35.);
        XYSel.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYSel.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYSel.Draw();
        c.SaveAs("Run%d_StripXYSelected.png" % self.run);
        
        DXDY = self.f0.Get("/FEB/Chamber1/Strips/DXDY");
        DXDY.SetTitle("Distance of nearest hit to the extrapolation %s" % self.comment);
        DXDY.GetYaxis().SetTitle("#delta Y_{Strip} (cm)");
        DXDY.GetXaxis().SetTitle("#delta X_{Strip} (cm)");
        DXDY.Draw();
        c.SaveAs("Run%d_StripDXDYSelected.png" % self.run);
        
        AbsoluteTimeDelta= self.f0.Get("/FEB/Chamber1/Strips/AbsoluteTimeDelta");
        AbsoluteTimeDelta.SetTitle("Other hits, Absolute Time distance to the selected hit time %s" % self.comment);
        AbsoluteTimeDelta.GetXaxis().SetTitle("#delta T_{Absolute} (ns)");
        AbsoluteTimeDelta.GetYaxis().SetTitle("Number of hits");
        AbsoluteTimeDelta.GetXaxis().SetRangeUser(-300.,300.);
        AbsoluteTimeDelta.Draw();
        c.SaveAs("Run%d_StripAbsoluteDeltaT.png" % self.run);
        AbsoluteTimeDelta.GetXaxis().SetRangeUser(-30.,130.);
        AbsoluteTimeDelta.Draw();
        c.SaveAs("Run%d_StripAbsoluteDeltaTZoomed.png" % self.run);

        self.f0.cd("/FEB/Chamber1/ClusterNew/");
        
        XYC = self.f0.Get("/FEB/Chamber1/ClusterNew/XY");									       
        XYC.SetTitle("Cluster position in the time window %s" % self.comment);
        XYC.GetXaxis().SetRangeUser(0.,35.);
        XYC.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYC.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYC.Draw();
        c.SaveAs("Run%d_ClusterXYAll.png" % self.run);
        
        XYMinC = self.f0.Get("/FEB/Chamber1/ClusterNew/XYMax");									       
        XYMinC.SetTitle("Cluster position nearest to the extrapolation %s" % self.comment);
        XYMinC.GetXaxis().SetRangeUser(0.,35.);
        XYMinC.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYMinC.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYMinC.Draw();
        c.SaveAs("Run%d_ClusterXYMin.png" % self.run);
        
        XYSelC = self.f0.Get("/FEB/Chamber1/ClusterNew/XYSel");
        XYSelC.SetTitle("Cluster associated %s" % self.comment);
        XYSelC.GetXaxis().SetRangeUser(0.,35.);
        XYSelC.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYSelC.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYSelC.Draw();
        c.SaveAs("Run%d_ClusterXYSelected.png" % self.run);
        
        DXDYC = self.f0.Get("/FEB/Chamber1/ClusterNew/DIST");
        DXDYC.SetTitle("Distance of nearest cluster to the extrapolation %s" % self.comment);
        DXDYC.GetYaxis().SetTitle("#delta Y_{Strip} (cm)");
        DXDYC.GetXaxis().SetTitle("#delta X_{Strip} (cm)");
        DXDYC.Draw();
        c.SaveAs("Run%d_ClusterDXDYSelected.png" % self.run);
        
        hnclus= self.f0.Get("/FEB/Chamber1/ClusterNew/Clusters");
        hnclus.SetTitle("Number of clusters %s" % self.comment);
        hnclus.GetXaxis().SetTitle("Number of clusters reconstructed in the event");
        hnclus.GetYaxis().SetTitle("Number of events");
        hnclus.Draw();
        c.SaveAs("Run%d_ClusterNclus.png" % self.run);
    
        cluss= self.f0.Get("/FEB/Chamber1/ClusterNew/ClusterSize");
        cluss.SetTitle("All clusters size %s" % self.comment);
        cluss.GetXaxis().SetTitle("Number of strip hits in cluster");
        cluss.GetYaxis().SetTitle("Number of events");
        cluss.Draw();
        c.SaveAs("Run%d_ClusterClusterSize.png" % self.run);
    
        cluss1= self.f0.Get("/FEB/Chamber1/ClusterNew/ClusterSize1");
        cluss1.SetTitle("Associated cluster  size %s" % self.comment);
        cluss1.GetXaxis().SetTitle("Number of strip hits in cluster");
        cluss1.GetYaxis().SetTitle("Number of events");
        cluss1.Draw();
        c.SaveAs("Run%d_ClusterClusterSizeAssociated.png" % self.run);
 



