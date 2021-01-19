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
        try:
            self.f1 = TFile(direc+"/Noisehisto%d_0.root" % run)
        except:
            self.f1=None
        #self.f1=None
        self.P=0
        self.P0=0
        self.HVapp=0
        self.HVeff=0
        self.T0=0
        self.T=0
        self.fC=0
        self.Threshold=0
        self.jsfile=jsfile
        f=open(jsfile)
        self.params=json.loads(f.read())
        fr=open("etc/runsum.json")
        self.runparam=json.loads(fr.read())
        for x,y in self.runparam.items():
            if (int(x)==run):
                self.HVapp=y["HVApplied"]
                self.Threshold=y["Threshold"]
                if (y["TDome"]==0):
                    self.T0=y["TMeteo"]
                    self.T=self.T0+273.15+3.0
                else:
                    self.T=y["TDome"]
                if (y["PDome"]==0):
                    self.P0=y["PMeteo"]
                    self.P= self.calP(self.P0,170)
                else:
                    self.P=y["PDome"]
                self.HVeff=self.calV(self.HVapp,self.P,self.T)
                self.fC=(self.Threshold-480)*2.4
                self.LR=self.Threshold
                if ("LR" in y):
                    self.LR=self.Threshold+y["LR"]
                self.lrfC=(self.LR-480)*2.4
                self.comment="R%d (%d/%d) VTH %d(%.1f/%.1f fC)" % (self.run,self.HVapp,self.HVeff,self.Threshold,self.fC,self.lrfC)

                #print self.HVapp,self.HVeff,self.Threshold,self.fC,self.P,self.T
                #print self.comment
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
                xm=-100000.
                maxh=0
                for ib in range(1,hst0.GetNbinsX()):
                    #print ib,hst0.GetBinContent(ib),maxh,hst0.GetXaxis().GetBinCenter(ib)
                    if (hst0.GetBinContent(ib)>maxh):
                        maxh=hst0.GetBinContent(ib)
                        xm=hst0.GetXaxis().GetBinCenter(ib)
                        #print "found"
                hst0.GetXaxis().SetRangeUser(xm-5.,xm+5.);
                dT0[ist]=hst0.GetMean()
                print xm,xm-5.,xm+5.,hst0.GetMean()
                c.cd()
                hst0.Draw()
                c.Modified()
                c.Update()
                val=raw_input()
            if (hst1!=None):
                xm=-100000.
                maxh=0
                for ib in range(1,hst1.GetNbinsX()):
                    if (hst1.GetBinContent(ib)>maxh):
                        maxh=hst1.GetBinContent(ib)
                        xm=hst1.GetXaxis().GetBinCenter(ib)
                hst1.GetXaxis().SetRangeUser(xm-5.,xm+5.);
                dT1[ist]=hst1.GetMean()
                print xm,hst1.GetMean()
                c.cd()
                hst1.Draw()
                c.Modified()
                c.Update()
                val=raw_input()
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
            ST0.append(round(T0[i],2))
        str_res=str_res+"],\n"
        str_res+='"st1":[0.'
        for i in range(2,48):
            str_res=str_res+",%.2f" % T1[i]
            ST1.append(round(T1[i],2))
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
            
            if (hst.GetEntries()<400):
                hst.Rebin(2)
                
            hst.SetTitle("Y_{ext} - Y_{strip} Strip %d %s" % (ist,self.comment))
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
            res[ist]=round(dtmean,2)
            dres[ist]=round(dtres,3)
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
        c.Divide(2,3)
        #self.f0.cd("/FEB/Chamber%d/Raw" % chamber)
        #hst=self.f0.Get("/FEB/Chamber%d/Raw/Strips" % chamber)
        #c.cd(1)
        #hst.Draw()
        #hrat=TH1F("hrat","Ratio LR/HR ",50,0,50.)
        #for i in range(1,33):
        #    x0=hst.GetBinContent(i+1)
        #    x1=hst.GetBinContent(i+1+48)
        #    r=0
        #    if (x0>x1):
        #        r=x1/x0
        #    #print x0,x1,r
        #    hrat.SetBinContent(i+1,r)


        #c.cd(2)
        #hrat.Draw()
        #c.Modified()
        #c.Update()
        #val=raw_input()
        noiseff=0
        if (self.f1!=None):
             self.f1.cd("/gric")
             hxyn=self.f1.Get("/gric/XYT")
             noiseff=hxyn.GetEntries()
        self.f0.cd("/gric")
        hxy=self.f0.Get("/gric/XY")
        hxy.SetTitle(" 4 points Track extrapolation to the iRPC %s" % self.comment)
        hxya=self.f0.Get("/gric/XYA")
        hxya14=self.f0.Get("/gric/XYA14")
        hxya15=self.f0.Get("/gric/XYA15")

        hxyf=self.f0.Get("/gric/XYF")
        hxyf14=self.f0.Get("/gric/XYF14")
        hxyf15=self.f0.Get("/gric/XYF15")
        hxyf.SetTitle(" 4 points Track extrapolation to the iRPC when cluster selected %s" % self.comment)
        hxyt=self.f0.Get("/gric/XYT")
        result={}
        efft=hxyt.GetEntries()/hxy.GetEntries()
        effn=noiseff/hxy.GetEntries()/2E-5/3300./efft
        result["run"]=int(self.run)
        result["threshold"]=int(self.Threshold)
        result["hrq"]=round(self.fC,1)
        result["lrq"]=round(self.lrfC,1)
        result["hvapp"]=round(self.HVapp,2)
        result["pressure"]=round(self.P,1)
        result["temperature"]=round(self.T,1)
        result["hveff"]=round(self.HVeff,2)
        result["ntk"]=int(hxy.GetEntries())
        result["nintime"]=int(hxyt.GetEntries())
        result["efftime"]=round(efft*100,2)
        result["noiserate"]=round(effn,2)
        effo=hxyf.GetEntries()/hxy.GetEntries()
        effo14=hxyf14.GetEntries()/hxy.GetEntries()
        effo15=hxyf15.GetEntries()/hxy.GetEntries()
        result["nfound"]=int(hxyf.GetEntries())
        result["effound"]=round(effo*100,2)
        result["nfound14"]=int(hxyf14.GetEntries())
        result["efffound14"]=round(effo14*100,2)
        result["nfound15"]=int(hxyf15.GetEntries())
        result["efffound15"]=round(effo15*100,2)
        #print hxy.GetEntries(),efft*100,effo*100,effo14*100,effo15*100
        c.cd(1)
        hxy.Draw("COLZ")
        c.cd(2)
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


        chxya=hxya.Clone("chxya")
        chxya14=hxya14.Clone("chxya14")
        chxya15=hxya15.Clone("chxya15")
        

        chxya.Rebin2D(rbx,rby)
        chxya14.Rebin2D(rbx,rby)
        chxya15.Rebin2D(rbx,rby)

        
        chxyt.Rebin2D(rbx,rby)
        chxy.SetAxisRange(0.,35.,"X")
        chxy.SetAxisRange(20.,120.,"Y")
        chxyf.SetAxisRange(0.,35.,"X")
        chxyf.SetAxisRange(20.,120.,"Y")
        chxyf14.SetAxisRange(0.,35.,"X")
        chxyf14.SetAxisRange(20.,120.,"Y")
        chxyf15.SetAxisRange(0.,35.,"X")
        chxyf15.SetAxisRange(20.,120.,"Y")
        chxya.SetAxisRange(0.,35.,"X")
        chxya.SetAxisRange(20.,120.,"Y")
        chxya14.SetAxisRange(0.,35.,"X")
        chxya14.SetAxisRange(20.,120.,"Y")
        chxya15.SetAxisRange(0.,35.,"X")
        chxya15.SetAxisRange(20.,120.,"Y")

        chxyEff=chxyf.Clone("chxyEff")
        chxyEff.Divide(chxy)
        chxyEff.SetTitle("Local Efficiency map %s " % self.comment)
        chxyEff14=chxyf14.Clone("chxyEff14")
        chxyEff14.Divide(chxy)
        chxyEff14.SetTitle("Local Efficiency map FR4 %s " % self.comment)
        chxyEff15=chxyf15.Clone("chxyEff15")
        chxyEff15.Divide(chxy)
        chxyEff15.SetTitle("Local Efficiency map EM888 %s " % self.comment)
        heff=TH1F("heff","Local Efficiency Ntk ext>15",110,0,1.1)
        dmin=max(0,effo*0.1)
        dmax=min(1.1,effo*1.9)
        nb=int((dmax-dmin)/0.005)
        heff1=TH1F("heff1","Local Efficiency EM888 %s " % self.comment,nb,dmin,dmax)
        heff32=TH1F("heff32","Local Efficiency FR4 %s " % self.comment,nb,dmin,dmax)

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
                    if (xp>=5 and xp<=15 and yp<84 and yp>41):
                        heff1.Fill(chxyEff15.GetBinContent(i,j))
                        #print chxy.GetBinContent(i,j), chxyf15.GetBinContent(i,j),chxyEff15.GetBinContent(i,j)
                    if (xp>=17 and xp<=26 and yp<84 and yp>41):
                        heff32.Fill(chxyEff14.GetBinContent(i,j))

        #print heff1.GetMean()*100,heff32.GetMean()*100

        deff1=(heff1.GetRMS()/math.sqrt(heff1.GetEntries()))*100
        deff32=(heff32.GetRMS()/math.sqrt(heff32.GetEntries()))*100
        result["nblem888"]=int(heff1.GetEntries())
        result["eflem888"]=round(heff1.GetMean()*100,2)
        result["deflem888"]=round(deff1,2)
        result["nblfr4"]=int(heff32.GetEntries())
        result["eflfr4"]=round(heff32.GetMean()*100,2)
        result["deflfr4"]=round(deff32,2)


        effao=hxya.GetEntries()/hxy.GetEntries()
        effao14=hxya14.GetEntries()/hxy.GetEntries()
        effao15=hxya15.GetEntries()/hxy.GetEntries()

        result["nand"]=int(hxya.GetEntries())
        result["effand"]=round(effao*100,2)
        result["nand14"]=int(hxya14.GetEntries())
        result["effand14"]=round(effao14*100,2)
        result["nand15"]=int(hxya15.GetEntries())
        result["effand15"]=round(effao15*100,2)

        
        chxyaEff=chxya.Clone("chxyaEff")
        chxyaEff.Divide(chxy)
        chxyaEff14=chxya14.Clone("chxyaEff14")
        chxyaEff14.Divide(chxy)
        chxyaEff14.SetTitle("Local Efficiency map FR4 %s " % self.comment)
        chxyaEff15=chxya15.Clone("chxyaEff15")
        chxyaEff15.Divide(chxy)
        chxyaEff15.SetTitle("Local Efficiency map EM888 %s " % self.comment)
        heffa=TH1F("heffa","Local Efficiency Ntk ext>15",110,0,1.1)
        dmin=max(0,effao*0.1)
        dmax=min(1.1,effao*1.9)
        nb=int((dmax-dmin)/0.005)
        heffa1=TH1F("heffa1","Local Efficiency EM888 %s " % self.comment,nb,dmin,dmax)
        heffa32=TH1F("heffa32","Local Efficiency FR4 %s " % self.comment,nb,dmin,dmax)

        for i in range(1,chxy.GetNbinsX()):
            for j in range(1,chxy.GetNbinsY()):
                if (chxy.GetBinContent(i,j)>ncut):
                    heffa.Fill(chxyaEff.GetBinContent(i,j))
                    xp=chxy.GetXaxis().GetBinCenter(i)
                    yp=chxy.GetYaxis().GetBinCenter(j)
                    if (xp>=5 and xp<=15 and yp<84 and yp>41):
                        heffa1.Fill(chxyaEff15.GetBinContent(i,j))
                        #print chxy.GetBinContent(i,j), chxyf15.GetBinContent(i,j),chxyEff15.GetBinContent(i,j)
                    if (xp>=17 and xp<=26 and yp<84 and yp>41):
                        heffa32.Fill(chxyaEff14.GetBinContent(i,j))

        #print heff1.GetMean()*100,heff32.GetMean()*100

        deffa1=(heffa1.GetRMS()/math.sqrt(heffa1.GetEntries()))*100
        deffa32=(heffa32.GetRMS()/math.sqrt(heffa32.GetEntries()))*100
        result["nbaem888"]=int(heffa1.GetEntries())
        result["efaem888"]=round(heffa1.GetMean()*100,2)
        result["defaem888"]=round(deffa1,2)
        result["nbafr4"]=int(heffa32.GetEntries())
        result["efafr4"]=round(heffa32.GetMean()*100,2)
        result["defafr4"]=round(deffa32,2)





        
        print result
       
        c.cd(3)
        chxyEff15.Draw("COLZ")
        c.Modified()
        c.Update()
        #val=raw_input()
        c.cd(4)
        chxyEff14.Draw("COLZ")
        c.Modified()
        c.Update()
        #val=raw_input()
        c.cd(5)
        heff1.Draw()
        c.Modified()
        c.Update()
        #val=raw_input()
        c.cd(6)
        heff32.Draw()
        c.Modified()
        c.Update()
        #
        #c.SaveAs("Analyse%d.pdf" % self.run)
        #val=raw_input()
        return result
    def pltres(self):
        c=TCanvas("c","IRPC studies %d" % self.run ,545,345);
        c.GetFrame().SetBorderSize( 12 )
        gStyle.SetOptStat(0)
        self.f0.cd("/FEB/Chamber1/Strips/");
        
        XY = self.f0.Get("/FEB/Chamber1/Strips/XY");									       
        XY.SetTitle("Strip hits position in the time window %s" % self.comment);
        XY.GetXaxis().SetRangeUser(0.,35.);
        XY.GetYaxis().SetRangeUser(0.,200.);
        XY.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XY.GetXaxis().SetTitle("X_{Strip} (cm)");
        XY.Draw();
        t8 =TText(8.,120.,"EM888");
        t8.SetTextAlign(22);
        t8.SetTextColor(2);
        t8.SetTextFont(43);
        t8.SetTextSize(12);
        t8.SetTextAngle(90);
        t8.Draw("SAME");
        t4 =TText(22.,120.,"FR4");
        t4.SetTextAlign(22);
        t4.SetTextColor(3);
        t4.SetTextFont(43);
        t4.SetTextSize(12);
        t4.SetTextAngle(90);
        t4.Draw("SAME");
        l=TLine(0.,150.0,35.,150.);
        l.SetLineColor(4);
        l.SetLineWidth(2);
        l.Draw("SAME");
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_StripXYAll.png" % self.run);
        
        XYMin = self.f0.Get("/FEB/Chamber1/Strips/XYMin");									       
        XYMin.SetTitle("Strip hits position nearest to the extrapolation %s" % self.comment);
        XYMin.GetXaxis().SetRangeUser(0.,35.);
        XYMin.GetYaxis().SetRangeUser(0.,200.);
        XYMin.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYMin.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYMin.Draw();
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_StripXYMin.png" % self.run);
        
        XYSel = self.f0.Get("/FEB/Chamber1/Strips/XYSel");
        XYSel.SetTitle("Strip hit associated %s" % self.comment);
        XYSel.GetXaxis().SetRangeUser(0.,35.);
        XYSel.GetYaxis().SetRangeUser(0.,200.);
        XYSel.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYSel.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYSel.Draw();
        c.Modified()
        c.Update()
        val=raw_input()

        c.SaveAs("Run%d_StripXYSelected.png" % self.run);
        
        DXDY = self.f0.Get("/FEB/Chamber1/Strips/DXDY");
        DXDY.SetTitle("Distance of nearest hit to the extrapolation %s" % self.comment);
        DXDY.GetYaxis().SetTitle("#delta Y_{Strip} (cm)");
        DXDY.GetXaxis().SetTitle("#delta X_{Strip} (cm)");
        DXDY.Draw();
        c.Modified()
        c.Update()
        val=raw_input()

        c.SaveAs("Run%d_StripDXDYSelected.png" % self.run);
        
        self.f0.cd("/FEB/Chamber1/ClusterNew/");
        
        XYC = self.f0.Get("/FEB/Chamber1/ClusterNew/XY");									       
        XYC.SetTitle("Cluster position in the time window %s" % self.comment);
        XYC.GetXaxis().SetRangeUser(0.,35.);
        XYC.GetYaxis().SetRangeUser(0.,200.);
        XYC.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYC.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYC.Draw();
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_ClusterXYAll.png" % self.run);
        
        XYMinC = self.f0.Get("/FEB/Chamber1/ClusterNew/XYMax");									       
        XYMinC.SetTitle("Cluster position nearest to the extrapolation %s" % self.comment);
        XYMinC.GetXaxis().SetRangeUser(0.,35.);
        XYMinC.GetYaxis().SetRangeUser(0.,200.);
        XYMinC.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYMinC.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYMinC.Draw();
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_ClusterXYMin.png" % self.run);
        
        XYSelC = self.f0.Get("/FEB/Chamber1/ClusterNew/XYSel");
        XYSelC.SetTitle("Cluster associated %s" % self.comment);
        XYSelC.GetXaxis().SetRangeUser(0.,35.);
        XYSelC.GetYaxis().SetRangeUser(0.,200.);
        XYSelC.GetYaxis().SetTitle("Y_{Strip} (cm)");
        XYSelC.GetXaxis().SetTitle("X_{Strip} (cm)");
        XYSelC.Draw();
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_ClusterXYSelected.png" % self.run);
        
        DXDYC = self.f0.Get("/FEB/Chamber1/ClusterNew/DIST");
        DXDYC.SetTitle("Distance of nearest cluster to the extrapolation %s" % self.comment);
        DXDYC.GetYaxis().SetTitle("#delta Y_{Strip} (cm)");
        DXDYC.GetXaxis().SetTitle("#delta X_{Strip} (cm)");
        DXDYC.Draw();
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_ClusterDXDYSelected.png" % self.run);
        
        hnclus= self.f0.Get("/FEB/Chamber1/ClusterNew/Clusters");
        hnclus.SetTitle("Number of clusters %s" % self.comment);
        hnclus.GetXaxis().SetTitle("Number of clusters reconstructed in the event");
        hnclus.GetYaxis().SetTitle("Number of events");
        hnclus.Draw();
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_ClusterNclus.png" % self.run);
    
        cluss= self.f0.Get("/FEB/Chamber1/ClusterNew/ClusterSize");
        cluss.SetTitle("All clusters size %s" % self.comment);
        cluss.GetXaxis().SetTitle("Number of strip hits in cluster");
        cluss.GetYaxis().SetTitle("Number of events");
        cluss.Draw();
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_ClusterClusterSize.png" % self.run);
    
        cluss1= self.f0.Get("/FEB/Chamber1/ClusterNew/ClusterSize1");
        cluss1.SetTitle("Associated cluster  size %s" % self.comment);
        cluss1.GetXaxis().SetTitle("Number of strip hits in cluster");
        cluss1.GetYaxis().SetTitle("Number of events");
        cluss1.Draw();
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_ClusterClusterSizeAssociated.png" % self.run);


        c.Clear()
        c.Divide(2,2)

        c.cd(1)
        hrsel=self.f0.Get("/FEB/Clustering1D/HRCount")
        hrsel.GetXaxis().SetTitle("Number of strip")
        hrsel.GetYaxis().SetTitle("Number of events")
        lrsel=self.f0.Get("/FEB/Clustering1D/LRCount")
        lrsel.GetXaxis().SetTitle("Number of strip")
        lrsel.GetYaxis().SetTitle("Number of events")
        hrmis=self.f0.Get("/Missed/Clustering1D/HRCount")
        hrmis.GetXaxis().SetTitle("Number of strip")
        hrmis.GetYaxis().SetTitle("Number of events")
        lrmis=self.f0.Get("/Missed/Clustering1D/LRCount")
        lrmis.GetXaxis().SetTitle("Number of strip")
        lrmis.GetYaxis().SetTitle("Number of events")
        gStyle.SetOptStat(110)
        gStyle.SetOptTitle(0);
        hrsel.Draw()
        c.cd(3)
        hrmis.Draw()
        c.cd(2)
        lrsel.Draw()
        c.cd(4)
        lrmis.Draw()
        c.Modified()
        c.Update()
        val=raw_input()
        c.SaveAs("Run%d_MissedPattern.png" % self.run);
        

    def setTDRStyle(self):

        tdrStyle = TStyle("tdrStyle","Style for P-TDR");

        # For the canvas:
        tdrStyle.SetCanvasBorderMode(0);
        tdrStyle.SetCanvasColor(kWhite);
        tdrStyle.SetCanvasDefH(600); #Height of canvas
        tdrStyle.SetCanvasDefW(600); #Width of canvas
        tdrStyle.SetCanvasDefX(0);   #POsition on screen
        tdrStyle.SetCanvasDefY(0);

        # For the Pad:
        tdrStyle.SetPadBorderMode(0);
        # tdrStyle.SetPadBorderSize(Width_t size = 1);
        tdrStyle.SetPadColor(kWhite);
        tdrStyle.SetPadGridX(False);
        tdrStyle.SetPadGridY(False);
        tdrStyle.SetGridColor(0);
        tdrStyle.SetGridStyle(3);
        tdrStyle.SetGridWidth(1);

        # For the frame:
        tdrStyle.SetFrameBorderMode(0);
        tdrStyle.SetFrameBorderSize(1);
        tdrStyle.SetFrameFillColor(0);
        tdrStyle.SetFrameFillStyle(0);
        tdrStyle.SetFrameLineColor(1);
        tdrStyle.SetFrameLineStyle(1);
        tdrStyle.SetFrameLineWidth(1);
        
        # For the histo:
        # tdrStyle.SetHistFillColor(1);
        # tdrStyle.SetHistFillStyle(0);
        tdrStyle.SetHistLineColor(1);
        tdrStyle.SetHistLineStyle(0);
        tdrStyle.SetHistLineWidth(1);
        # tdrStyle.SetLegoInnerR(Float_t rad = 0.5);
        # tdrStyle.SetNumberContours(Int_t number = 20);

        tdrStyle.SetEndErrorSize(2);
        # tdrStyle.SetErrorMarker(20);
        #tdrStyle.SetErrorX(0.);
        
        tdrStyle.SetMarkerStyle(20);
  
        #For the fit/function:
        tdrStyle.SetOptFit(1);
        tdrStyle.SetFitFormat("5.4g");
        tdrStyle.SetFuncColor(2);
        tdrStyle.SetFuncStyle(1);
        tdrStyle.SetFuncWidth(1);
    
        #For the date:
        tdrStyle.SetOptDate(0);
        # tdrStyle.SetDateX(Float_t x = 0.01);
        # tdrStyle.SetDateY(Float_t y = 0.01);
        
        # For the statistics box:
        tdrStyle.SetOptFile(0);
        tdrStyle.SetOptStat(0); # To display the mean and RMS:   SetOptStat("mr");
        tdrStyle.SetStatColor(kWhite);
        tdrStyle.SetStatFont(42);
        tdrStyle.SetStatFontSize(0.025);
        tdrStyle.SetStatTextColor(1);
        tdrStyle.SetStatFormat("6.4g");
        tdrStyle.SetStatBorderSize(1);
        tdrStyle.SetStatH(0.1);
        tdrStyle.SetStatW(0.15);
        # tdrStyle.SetStatStyle(Style_t style = 1001);
        # tdrStyle.SetStatX(Float_t x = 0);
        # tdrStyle.SetStatY(Float_t y = 0);

        # Margins:
        tdrStyle.SetPadTopMargin(0.10);
        tdrStyle.SetPadBottomMargin(0.13);
        tdrStyle.SetPadLeftMargin(0.16);
        tdrStyle.SetPadRightMargin(0.02);

        # For the Global title:

        tdrStyle.SetOptTitle(1);
        tdrStyle.SetTitleFont(42);
        tdrStyle.SetTitleColor(1);
        tdrStyle.SetTitleTextColor(1);
        tdrStyle.SetTitleFillColor(10);
        tdrStyle.SetTitleFontSize(0.05);
        tdrStyle.SetTitleH(0); # Set the height of the title box
        tdrStyle.SetTitleW(0); # Set the width of the title box
        tdrStyle.SetTitleX(0); # Set the position of the title box
        tdrStyle.SetTitleY(0.985); # Set the position of the title box
        tdrStyle.SetTitleStyle(1001);
        tdrStyle.SetTitleBorderSize(2);

        # For the axis titles:
        
        tdrStyle.SetTitleColor(1, "XYZ");
        tdrStyle.SetTitleFont(42, "XYZ");
        tdrStyle.SetTitleSize(0.06, "XYZ");
        # tdrStyle.SetTitleXSize(Float_t size = 0.02); # Another way to set the size?
        # tdrStyle.SetTitleYSize(Float_t size = 0.02);
        tdrStyle.SetTitleXOffset(0.9);
        tdrStyle.SetTitleYOffset(1.25);
        # tdrStyle.SetTitleOffset(1.1, "Y"); # Another way to set the Offset

        # For the axis labels:

        tdrStyle.SetLabelColor(1, "XYZ");
        tdrStyle.SetLabelFont(42, "XYZ");
        tdrStyle.SetLabelOffset(0.007, "XYZ");
        tdrStyle.SetLabelSize(0.05, "XYZ");

        # For the axis:

        tdrStyle.SetAxisColor(1, "XYZ");
        tdrStyle.SetStripDecimals(True);
        tdrStyle.SetTickLength(0.03, "XYZ");
        tdrStyle.SetNdivisions(510, "XYZ");
        tdrStyle.SetPadTickX(1);  # To get tick marks on the opposite side of the frame
        tdrStyle.SetPadTickY(1);
        
        # Change for log plots:
        tdrStyle.SetOptLogx(0);
        tdrStyle.SetOptLogy(0);
        tdrStyle.SetOptLogz(0);

        # Postscript options:
        tdrStyle.SetPaperSize(20.,20.);
        # tdrStyle.SetLineScalePS(Float_t scale = 3);
        # tdrStyle.SetLineStyleString(Int_t i, const char* text);
        # tdrStyle.SetHeaderPS(const char* header);
        # tdrStyle.SetTitlePS(const char* pstitle);
        
        # tdrStyle.SetBarOffset(Float_t baroff = 0.5);
        # tdrStyle.SetBarWidth(Float_t barwidth = 0.5);
        # tdrStyle.SetPaintTextFormat(const char* format = "g");
        # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
        # tdrStyle.SetTimeOffset(Double_t toffset);
        # tdrStyle.SetHistMinimumZero(kTRUE);
        
        tdrStyle.SetHatchesLineWidth(5);
        tdrStyle.SetHatchesSpacing(0.05);

        tdrStyle.cd();

