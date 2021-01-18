from __future__ import print_function
from ROOT import TCanvas, TGraph,TLegend
from ROOT import gROOT
from math import sin
from array import array
import json

c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )

#c1.SetFillColor( 42 )
c1.SetGrid()
fout=open("etc/SummaryPR2.json")
result=json.loads(fout.read())
print(result)

#runs=[1786,1787,1788,1789,1790,1791,1792]
#runs=[1795,1801,1802,1803,1804,1805]
runs=[1819,1820,1821,1822,1823]
#runs=[1824,1825,1826,1827,1828]
n = 0
x, y4, y5,ya4,ya5 = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
hv=0
thrc=""
leff=100
for i in runs:
    r=result["%d" %i]
    #x.append(r["hrq"])
    x.append(r["hveff"])
    y4.append(r["eflfr4"] )
    y5.append(r["eflem888"] )
    ya4.append(r["efafr4"] )
    ya5.append(r["efaem888"] )
    hv=r["hveff"]
    thrc="Cut %.1f/%.1f" % (r["hrq"],r["lrq"])
    print(' i %i %f %f ' % (i,x[n],y5[n]))
    if (y4[n]<leff):
        leff=y4[n]
    if (y5[n]<leff):
        leff=y5[n]
    n=n+1
 
gr4 = TGraph( n, x, y4 )
gr4.SetLineColor( 2 )
gr4.SetLineWidth( 4 )
gr4.SetMarkerColor( 1 )
gr4.SetMarkerStyle( 20 )
#gr4.SetTitle( 'Threshold scan HV=%f' % hv )
gr4.SetTitle( 'HV scan  Threshold=%s' % thrc )
gr4.GetXaxis().SetTitle( 'HV (V)' )
gr4.GetYaxis().SetTitle( 'Efficiency (%)' )
gr4.GetYaxis().SetRangeUser(leff-5,103.)
gr4.Draw( 'AP' )
gr5 = TGraph( n, x, y5 )
gr5.SetMarkerColor( 2 )
gr5.SetMarkerStyle( 22 )
gr5.SetTitle( 'a simple graph' )
gr5.GetXaxis().SetTitle( 'X5 title' )
gr5.GetYaxis().SetTitle( 'Y5 title' )
gr5.Draw( 'PSAME' )


gra4 = TGraph( n, x, ya4 )
gra4.SetMarkerColor( 3 )
gra4.SetMarkerStyle( 20 )
gra4.Draw( 'PSAME' )
gra5 = TGraph( n, x, ya5 )
gra5.SetMarkerColor( 4 )
gra5.SetMarkerStyle( 22 )
gra5.Draw( 'PSAME' )

legend =TLegend()
#legend.SetHeader("Header","C")
legend.AddEntry(gra5,"EM888 HR & LR","p");
legend.AddEntry(gr5,"EM888 d<10 cm","p");
legend.AddEntry(gra4,"FR4 HR & LR ","p");
legend.AddEntry(gr4,"FR4 d<10 cm","p");

legend.Draw();
 # TCanvas.Update() draws the frame, after which one can change it
c1.Update()
#c1.GetFrame().SetFillColor( 21 )
c1.GetFrame().SetBorderSize( 12 )
c1.Modified()
c1.Update()
val=raw_input()
#print val
c1.SaveAs("%s.png" % val)
