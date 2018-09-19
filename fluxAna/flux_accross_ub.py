import sys
import ROOT
from ROOT import TFile, TTree, TH1D, TCanvas, gROOT, TPad, TGaxis, TColor, TLegend
#import matplotlib

file1_path = "/uboone/data/users/guzowski/numi_flux/fgd.root"
file1 = TFile(file1_path, 'READ')
if(file1.IsOpen()):
    print 'File ', file1_path, ' is open'
if(file1.IsOpen() == False):
    quit()

# these are the corrected gSimple flux histograms
flux_nue_000 = file1.Get("nueFluxHisto000")
flux_nuebar_000 = file1.Get("anueFluxHisto000")
flux_nue_555 = file1.Get("nueFluxHisto555")
flux_nuebar_555 = file1.Get("anueFluxHisto555")
flux_nue_999 = file1.Get("nueFluxHisto999")
flux_nuebar_999 = file1.Get("anueFluxHisto999")

c1 = TCanvas("c1", "c1", 800, 600)
c1.cd()

flux_nue_555.SetLineColor(46)
flux_nue_999.SetLineColor(32)

flux_nue_000.Draw("hist")
flux_nue_555.Draw("hist same")
flux_nue_999.Draw("hist same")

c1.Print("plots/nue_flux_inside_detector_position.pdf")

c2 = TCanvas("c2", "c2", 800, 600)
c2.cd()

flux_nuebar_555.SetLineColor(46)
flux_nuebar_999.SetLineColor(32)

flux_nuebar_000.Draw("hist")
flux_nuebar_555.Draw("hist same")
flux_nuebar_999.Draw("hist same")

c2.Print("plots/nuebar_flux_inside_detector_position.pdf")

wait = raw_input("PRESS ENTER TO CONTINUE.")
