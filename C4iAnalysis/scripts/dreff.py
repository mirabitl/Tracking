from ROOT import *
import accessHisto as ah


def analyze(run,seq=0):

  f2=TFile("/tmp/histo%d_%d.root" % (run,seq) )
  f2.cd("/gric")
  print "%.2f %.6f %.2f \n" % (MaxTime.GetEntries(),MaxTime.GetMean(),MaxTime.GetEntries()*MaxTime.GetMean())
  #print f2
  a=[]
  rb=4
  a.append(f2)
  for ip in range(1,9):
    hext1=ah.getth2("/gric/PLANE%d/XYext" % ip)
    hext1.Rebin2D(rb,rb)
    hfo1=ah.getth2("/gric/PLANE%d/XYfound" %ip )
    hfo1.Rebin2D(rb,rb)
    hfo1.Divide(hext1)
    hsum1=TH1F("sum%d" % ip ,"sum%d" % ip,110,0.,1.1)
    for i in range(hext1.GetNbinsX()):
      for j in range(hext1.GetNbinsY()):
        if (hext1.GetBinContent(i+1,j+1)>15):
          hsum1.Fill(hfo1.GetBinContent(i+1,j+1))

    hsum1.Draw()
    a.append(hext1)
    a.append(hfo1)
    a.append(hsum1)

 
  return a

