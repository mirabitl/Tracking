from __future__ import print_function
from ROOT import TCanvas, TGraph
from ROOT import gROOT
from math import sin
from array import array
import json

c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )

c1.SetFillColor( 42 )
c1.SetGrid()
fout=open("SummaryPR1.json")
result=json.loads(fout.read())
print(result)

#runs=[1786,1787,1788,1789,1790,1791,1792]
#runs=[1795,1801,1802,1803,1804,1805]
runs=[1792,1793,1794,1795,1797,1798,1799,1800]
n = 0
x, y4, y5,ya4,ya5 = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
hv=0 
for i in runs:
    r=result["%d" %i]
    #x.append(r["hrq"])
    x.append(r["hveff"])
    y4.append(r["eflfr4"] )
    y5.append(r["eflem888"] )
    ya4.append(r["efafr4"] )
    ya5.append(r["efaem888"] )
    #hv=r["hveff"]
    hv=r["hrq"]
    print(' i %i %f %f ' % (i,x[n],y5[n]))
    n=n+1
 
gr4 = TGraph( n, x, y4 )
gr4.SetLineColor( 2 )
gr4.SetLineWidth( 4 )
gr4.SetMarkerColor( 4 )
gr4.SetMarkerStyle( 20 )
gr4.SetTitle( 'Threshold scan HV=%f' % hv )
gr4.GetXaxis().SetTitle( 'X title' )
gr4.GetYaxis().SetTitle( 'Y title' )
gr4.GetYaxis().SetRangeUser(20.,103.)
gr4.Draw( 'AP' )
gr5 = TGraph( n, x, y5 )
gr5.SetMarkerColor( 5 )
gr5.SetMarkerStyle( 22 )
gr5.SetTitle( 'a simple graph' )
gr5.GetXaxis().SetTitle( 'X5 title' )
gr5.GetYaxis().SetTitle( 'Y5 title' )
gr5.Draw( 'PSAME' )


gra4 = TGraph( n, x, ya4 )
gra4.SetMarkerColor( 6 )
gra4.SetMarkerStyle( 20 )
gra4.Draw( 'PSAME' )
gra5 = TGraph( n, x, ya5 )
gra5.SetMarkerColor( 7 )
gra5.SetMarkerStyle( 22 )
gra5.Draw( 'PSAME' )



 # TCanvas.Update() draws the frame, after which one can change it
c1.Update()
c1.GetFrame().SetFillColor( 21 )
c1.GetFrame().SetBorderSize( 12 )
c1.Modified()
c1.Update()
val=raw_input()
