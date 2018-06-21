import sys
import ROOT
from ROOT import TFile, TTree, TH1D, TCanvas, gROOT
import matplotlib

file1_path = "ppfx_ME000Z200i_NuMI_at_MicroBooNE.root"
file1 = TFile(file1_path, 'READ')
if(file1.IsOpen()):
    print 'File ', file1_path, ' is open'
if(file1.IsOpen() == False):
    quit()
dir1 = file1.results

# these are the central value histograms from the dk2nu flux version
dk2nu_nue = dir1.Get("hcv_nue")
dk2nu_nuebar = dir1.Get("hcv_nuebar")

file2_path = "../arxiv/NuMIFlux.root"
file2 = TFile(file2_path, 'READ')
if(file2.IsOpen()):
    print 'File ', file2_path, ' is open'
if(file2.IsOpen() == False):
    quit()

# these are the flux histograms coming from the gsimple files
simple_nue = file2.Get("nueFluxHisto")
simple_nuebar = file2.Get("anueFluxHisto")

# gsimple flux is in /cm2 and dk2nu is in /m2
# converting to /cm2
dk2nu_nue.Scale(1. / (100. * 100.))
dk2nu_nuebar.Scale(1. / (100. * 100.))

# gsimple is per 6e20 POT and dk2nu is per 1e6
# converting dk2nu to per 6e20
dk2nu_nue.Scale(6. * pow(10, 14))
dk2nu_nuebar.Scale(6. * pow(10, 14))

dk2nu_nue.GetYaxis().SetTitle("#nu / cm2 / 6e20 POT / GeV")
dk2nu_nuebar.GetYaxis().SetTitle("#nu / cm2 / 6e20 POT / GeV")

c1 = TCanvas("c1", "c1", 800, 600)
c1.cd()

dk2nu_nue.Draw("hist")
dk2nu_nuebar.Draw("hist same")

c2 = TCanvas("c2", "c2", 800, 600)
c2.cd()

simple_nue.Draw("hist")
simple_nuebar.Draw("hist same")

# we also want their ratios

# first change the binning
# print dk2nu_nue.GetNbinsX(), dk2nu_nuebar.GetNbinsX()
# print '-----------------------------------------'
# print simple_nue.GetNbinsX(), simple_nuebar.GetNbinsX()

rebin_nue = TH1D("rebin_nue", "rebin_nue", 400, 0, 20)
rebin_nuebar = TH1D("rebin_nuebar", "rebin_nuebar", 400, 0, 20)
ratio_nue = TH1D("ratio_nue", "ratio_nue", 400, 0, 20)
ratio_nuebar = TH1D("ratio_nuebar", "ratio_nuebar", 400, 0, 20)

dk2nu_nue_bins = dk2nu_nue.GetNbinsX()
dk2nu_nuebar_bins = dk2nu_nuebar.GetNbinsX()

for i in range(dk2nu_nue_bins):
    # skip the underflow bin
    if(i == 0):
        continue
    bin_content = dk2nu_nue.GetBinContent(i)
    bin_center = dk2nu_nue.GetBinCenter(i)
    bin_width = dk2nu_nue.GetBinWidth(i)
    # print bin_content, bin_center
    rebin_nue.Fill(bin_center, bin_content * bin_width)
    ratio_nue.Fill(bin_center, bin_content * bin_width)

for i in range(dk2nu_nuebar_bins):
    # skip the underflow bin
    if(i == 0):
        continue
    bin_content = dk2nu_nuebar.GetBinContent(i)
    bin_center = dk2nu_nuebar.GetBinCenter(i)
    bin_width = dk2nu_nuebar.GetBinWidth(i)
    # print bin_content, bin_center
    rebin_nuebar.Fill(bin_center, bin_content * bin_width)
    ratio_nuebar.Fill(bin_center, bin_content * bin_width)

c3a = TCanvas("c3a", "c3a", 800, 600)
c3a.cd()
rebin_nue.GetYaxis().SetTitle("#nu / cm2 / 6e20 POT")
rebin_nue.Draw("hist")

c3b = TCanvas("c3b", "c3b", 800, 600)
c3b.cd()
rebin_nuebar.GetYaxis().SetTitle("#nu / cm2 / 6e20 POT")
rebin_nuebar.Draw("hist")

c4 = TCanvas("c4", "c4", 800, 600)
c4.cd()
ratio_nue.Divide(simple_nue)
ratio_nue.GetYaxis().SetTitle("dk2nu : gsimple Fluxes (#nu e)")
ratio_nue.Draw("ep")

c5 = TCanvas("c5", "c5", 800, 600)
c5.cd()
ratio_nuebar.Divide(simple_nuebar)
ratio_nuebar.GetYaxis().SetTitle("dk2nu : gsimple Fluxes (#nu e-bar)")
ratio_nuebar.Draw("ep")

wait = raw_input("PRESS ENTER TO CONTINUE.")
