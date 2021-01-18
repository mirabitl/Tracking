from __future__ import print_function
from ROOT import TCanvas, TGraph,TGraphErrors,TLegend
from ROOT import gROOT
from math import sin,sqrt
from array import array
import json
c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )

#c1.SetFillColor( 42 )
c1.SetGrid()
fout=open("etc/SummaryPR1.json")
result=json.loads(fout.read())
print(result)

#runs=[1786,1787,1788,1789,1790,1791,1792]
#runs=[1795,1801,1802,1803,1804,1805]
runs=[1792,1793,1794,1795,1797,1798,1799,1800]
n = 0
x, y4, y5,ya4,ya5,yt,nr = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
dx, dy4, dy5,dya4,dya5,dyt = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
hv=0
leff=100
thrc=""
for i in runs:
    r=result["%d" %i]
    #x.append(r["hrq"])
    x.append(r["hveff"])
    dx.append(1E-3*r["hveff"])
    yt.append(r["efftime"] )
    eff =r["efftime"]/100.
    nt=r["nintime"]
    deff=sqrt(eff*(1-eff)/nt)
    dyt.append(deff)
    y4.append(r["eflfr4"] )
    dy4.append(r["deflfr4"] )
    y5.append(r["eflem888"] )
    dy5.append(r["deflem888"] )
    ya4.append(r["efafr4"] )
    dya4.append(r["defafr4"] )
    ya5.append(r["efaem888"] )
    dya5.append(r["defaem888"] )
    hv=r["hveff"]
    thrc="(H%.0f/L%.0f) fC" % (r["hrq"],r["lrq"])
    print(' i %i %f %f ' % (i,x[n],y5[n]))
    if (y4[n]<leff):
        leff=y4[n]
    if (y5[n]<leff):
        leff=y5[n]
    n=n+1
 
gr4 = TGraphErrors( n, x, y4,dx,dy4 )

gr4.SetMarkerColor( 1 )
gr4.SetMarkerStyle( 20 )
gr4.SetTitle( 'HV scan Threshold=%s' % thrc )
#gr4.GetXaxis().SetTitle( 'Threshold (fC)' )
gr4.GetXaxis().SetTitle( 'HV_{eff} (V)' )
gr4.GetYaxis().SetTitle( 'Efficiency (%)' )


gr4.GetYaxis().SetRangeUser(leff-10,103.)
gr4.Draw( 'AP' )
gr5 = TGraphErrors( n, x, y5,dx,dy5 )
gr5.SetMarkerColor( 2 )
gr5.SetMarkerStyle( 22 )
gr5.Draw( 'PSAME' )


gra4 = TGraphErrors( n, x, ya4,dx,dya4 )
gra4.SetMarkerColor( 3 )
gra4.SetMarkerStyle( 20 )
gra4.Draw( 'PSAME' )
gra5 = TGraphErrors( n, x, ya5,dx,dya5 )
gra5.SetMarkerColor( 4 )
gra5.SetMarkerStyle( 22 )
gra5.Draw( 'PSAME' )

grt = TGraphErrors( n, x, yt,dx,dyt )
grt.SetMarkerColor( 6 )
grt.SetMarkerStyle( 23 )
grt.Draw( 'PSAME' )

legend =TLegend()
#legend.SetHeader("Header","C")
legend.AddEntry(gra5,"EM888 HR & LR","p");
legend.AddEntry(gr5,"EM888 d<10 cm","p");
legend.AddEntry(gra4,"FR4 HR & LR ","p");
legend.AddEntry(gr4,"FR4 d<10 cm","p");
legend.AddEntry(grt,"Full chamber HR || LR","p");

legend.Draw();

 # TCanvas.Update() draws the frame, after which one can change it
c1.Update()
#c1.GetFrame().SetFillColor( 21 )
c1.GetFrame().SetBorderSize( 12 )
c1.Modified()
c1.Update()
val=raw_input()
c1.SaveAs("%s.png" % val )
