#!/usr/bin/python2
from ROOT import *
import json
import math
import time
from array import array

class prettyfloat(float):
    def __repr__(self):
        return "%0.2f" % self

class chana:
    def __init__(self,fname):
        self.root_file=TFile(fname)
        self.c=TCanvas()
        self.c.Divide(2,2)
        self.c.Draw()

        
    def calP(self,p0,alt):
        return p0*(1-0.0065*alt/288.15)**5.255

    def calV(self,V,P,T):
        #print 1-(0.2+0.8*P/990.*293./T)
        return  V/(0.2+0.8*P/990.*293./T)

    def calValice(self,V,P,T):
        print 1-(1*P/990.*293./T)
        return  V/(1.0*P/990.*293./T)

    def calApp(self,V,P,T):
        print 1-(0.2+0.8*P/990.*293./T)
        return  V*(0.2+0.8*P/990.*293./T)

    def approx(self,v,v1,v2,a1,a2):
        return a1+(v-v1)*(a2-a1)*1./(v2-v1)

    def parse(self,fjson):
        f=open(fjson)
        self.par=json.loads(f.read())
        #print self.par
        f.close()
    def draweff(self,fjson,plane):
        self.parse(fjson)

        for x in self.par["setup"]:
            if (x["num"]!=plane):
                continue
            self.root_file.cd("/gric/EFF%d_G%d/" % (x["num"],x["id"]))
            hext=self.root_file.Get("/gric/EFF%d_G%d/XYext" % (x["num"],x["id"]))
            hfound=self.root_file.Get("/gric/EFF%d_G%d/XYfound" % (x["num"],x["id"]))
            n_ext=hext.GetEntries()
            n_found=hfound.GetEntries()
            hext.Rebin2D(2,2)
            hfound.Rebin2D(2,2)
            self.c.cd(1)
            hext.Draw("COLZ")
            self.c.cd(2)
            hfound.Draw("COLZ")
            self.c.Update()
            #b=raw_input()
            self.c.cd(3)
            self.heff2=hfound
            self.heff2.Divide(hext)
            self.heff2.Draw("COLZ")
            self.heff=TH1F("eff","Efficiency",202,0.,1.05)
            for i in range(1,17):
                for j in range(1,16):
                    t=hext.GetBinContent(i,j)
                    if (t>10):
                        self.heff.Fill(hfound.GetBinContent(i,j))
            self.c.cd(4)
            self.heff.Draw()
            print "Efficiency : %d %d %d %f %f \n" % (plane,n_ext,n_found,n_found*100./n_ext,self.heff.GetMean()*100)
            self.c.Update()
            #b=raw_input()
