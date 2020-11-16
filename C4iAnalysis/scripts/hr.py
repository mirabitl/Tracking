from ROOT import *
import accessHisto as ah

class hr:
  def __init__(self,run,seq):
      """
      Handle all application definition and p  arameters , It controls the acquisition via the FDAQ application and the Slow control via the FSLOW application
      """
      self.f=TFile("/tmp/histo%d_%d.root" % (run,seq) )
      self.canvas=TCanvas("Run_%d" % run)



  def draw2(self,name):
    self.f.cd()
    hd=ah.getth2(name)
    self.canvas.cd()
    hd.Draw("COLZ")
    return hd
  def draw1(self,name):
    self.f.cd()
    hd=ah.getth1(name)
    self.canvas.cd()
    hd.Draw("COLZ")
    return hd
