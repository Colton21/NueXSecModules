import sys
import ROOT
from ROOT import TFile, TTree, TH1D, TCanvas, gROOT, TPad, TGaxis, TColor, TLegend
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
dk2nu_numu = dir1.Get("hcv_numu")
dk2nu_numubar = dir1.Get("hcv_numubar")

file2_path = "../arxiv/NuMIFlux.root"
file2 = TFile(file2_path, 'READ')
if(file2.IsOpen()):
    print 'File ', file2_path, ' is open'
if(file2.IsOpen() == False):
    quit()

# these are the flux histograms coming from the gsimple files
simple_nue = file2.Get("nueFluxHisto")
simple_nuebar = file2.Get("anueFluxHisto")
simple_numu = file2.Get("numuFluxHisto")
simple_numubar = file2.Get("anumuFluxHisto")

# gsimple flux is in /cm2 and dk2nu is in /m2
# converting to /cm2
dk2nu_nue.Scale(1. / (100. * 100.))
dk2nu_nuebar.Scale(1. / (100. * 100.))
dk2nu_numu.Scale(1. / (100. * 100.))
dk2nu_numubar.Scale(1. / (100. * 100.))

# gsimple is per 6e20 POT and dk2nu is per 1e6
# converting dk2nu to per 6e20
dk2nu_nue.Scale(6. * pow(10, 14))
dk2nu_nuebar.Scale(6. * pow(10, 14))
dk2nu_numu.Scale(6. * pow(10, 14))
dk2nu_numubar.Scale(6. * pow(10, 14))

dk2nu_nue.GetYaxis().SetTitle("#nu / cm2 / 6e20 POT / GeV")
dk2nu_nuebar.GetYaxis().SetTitle("#nu / cm2 / 6e20 POT / GeV")
dk2nu_numu.GetYaxis().SetTitle("#nu / cm2 / 6e20 POT / GeV")
dk2nu_numubar.GetYaxis().SetTitle("#nu / cm2 / 6e20 POT / GeV")

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

rebin_numu = TH1D("rebin_nue", "rebin_nue", 400, 0, 20)
rebin_numubar = TH1D("rebin_nuebar", "rebin_nuebar", 400, 0, 20)
ratio_numu = TH1D("ratio_nue", "ratio_nue", 400, 0, 20)
ratio_numubar = TH1D("ratio_nuebar", "ratio_nuebar", 400, 0, 20)


dk2nu_nue_bins = dk2nu_nue.GetNbinsX()
dk2nu_nuebar_bins = dk2nu_nuebar.GetNbinsX()
dk2nu_numu_bins = dk2nu_numu.GetNbinsX()
dk2nu_numubar_bins = dk2nu_numubar.GetNbinsX()

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

for i in range(dk2nu_numu_bins):
    # skip the underflow bin
    if(i == 0):
        continue
    bin_content = dk2nu_numu.GetBinContent(i)
    bin_center = dk2nu_numu.GetBinCenter(i)
    bin_width = dk2nu_numu.GetBinWidth(i)
    # print bin_content, bin_center
    rebin_numu.Fill(bin_center, bin_content * bin_width)
    ratio_numu.Fill(bin_center, bin_content * bin_width)

for i in range(dk2nu_nuebar_bins):
    # skip the underflow bin
    if(i == 0):
        continue
    bin_content = dk2nu_numubar.GetBinContent(i)
    bin_center = dk2nu_numubar.GetBinCenter(i)
    bin_width = dk2nu_numubar.GetBinWidth(i)
    # print bin_content, bin_center
    rebin_numubar.Fill(bin_center, bin_content * bin_width)
    ratio_numubar.Fill(bin_center, bin_content * bin_width)

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
ratio_nue.GetYaxis().SetTitle("dk2nu : gsimple (#nu e)")
ratio_nue.Draw("ep")

c5 = TCanvas("c5", "c5", 800, 600)
c5.cd()
ratio_nuebar.Divide(simple_nuebar)
ratio_nuebar.GetYaxis().SetTitle("dk2nu : gsimple (#nu e-bar)")
ratio_nuebar.Draw("ep")

c4b = TCanvas("c4b", "c4b", 800, 600)
c4b.cd()
ratio_numu.Divide(simple_numu)
ratio_numu.GetYaxis().SetTitle("dk2nu : gsimple (#nu #mu)")
ratio_numu.Draw("ep")

c5b = TCanvas("c5b", "c5b", 800, 600)
c5b.cd()
ratio_numubar.Divide(simple_numubar)
ratio_numubar.GetYaxis().SetTitle("dk2nu : gsimple (#nu #mu-bar)")
ratio_numubar.Draw("ep")

#-----------------------------------------------------------------------

c6 = TCanvas("c6", "c6", 800, 800)
c6.cd()

pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad1.SetGridx()         # Vertical grid
pad1.Draw()             # Draw the upper pad: pad1
pad1.cd()               # pad1 becomes the current pad
rebin_nue.SetStats(0)   # No statistics on upper plot
simple_nue.SetStats(0)
rebin_nue.SetLineColor(38)
simple_nue.SetLineColor(46)
rebin_nue.SetLineWidth(6)
simple_nue.SetLineWidth(6)
simple_nue.GetXaxis().SetRangeUser(0, 4)
rebin_nue.GetXaxis().SetRangeUser(0, 4)

simple_nue.Draw("hist")
rebin_nue.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
rebin_nue.GetYaxis().SetLabelSize(0.)
axis = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis.SetLabelSize(15)
axis.Draw()

leg_nue = TLegend(0.75, 0.95, 0.75, 0.95)
leg_nue.AddEntry(simple_nue, "nue gSimple Flux",  "l")
leg_nue.AddEntry(rebin_nue,  "nue dk2nu Flux",    "l")
leg_nue.Draw()

c6.cd()
pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
pad2.SetTopMargin(0.1)
pad2.SetBottomMargin(0.1)
pad2.SetGridx()  # vertical grid
pad2.Draw()
pad2.cd()  # pad2 becomes the current pad

ratio_nue.SetLineColor(12)
ratio_nue.SetMinimum(0.0)  # Define Y ..
ratio_nue.SetMaximum(3.0)  # .. range
ratio_nue.GetXaxis().SetRangeUser(0, 4)
ratio_nue.Sumw2()
ratio_nue.SetStats(0)  # No statistics on lower plot
ratio_nue.SetMarkerStyle(21)
ratio_nue.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
rebin_nue.GetYaxis().SetTitleSize(20)
rebin_nue.GetYaxis().SetTitleFont(43)
rebin_nue.GetYaxis().SetTitleOffset(1.55)
ratio_nue.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_nue.GetYaxis().SetNdivisions(505)
ratio_nue.GetYaxis().SetTitleSize(14)
ratio_nue.GetYaxis().SetTitleFont(43)
ratio_nue.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_nue.GetYaxis().SetLabelFont(43)
ratio_nue.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_nue.GetXaxis().SetTitleSize(20)
ratio_nue.GetXaxis().SetTitleFont(43)
ratio_nue.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_nue.GetXaxis().SetLabelFont(43)
ratio_nue.GetXaxis().SetLabelSize(15)

#------------------------------------------------------------------

c7 = TCanvas("c7", "c7", 800, 800)
c7.cd()

pad3 = TPad("pad3", "pad3", 0, 0.3, 1, 1.0)
pad3.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad3.SetGridx()         # Vertical grid
pad3.Draw()             # Draw the upper pad: pad1
pad3.cd()               # pad1 becomes the current pad
rebin_nuebar.SetStats(0)   # No statistics on upper plot
simple_nuebar.SetStats(0)
rebin_nuebar.SetLineColor(38)
simple_nuebar.SetLineColor(46)
rebin_nuebar.SetLineWidth(6)
simple_nuebar.SetLineWidth(6)
simple_nuebar.GetXaxis().SetRangeUser(0, 4)
rebin_nuebar.GetXaxis().SetRangeUser(0, 4)

simple_nuebar.Draw("hist")
rebin_nuebar.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
rebin_nuebar.GetYaxis().SetLabelSize(0.)
axis2 = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis2.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis2.SetLabelSize(15)
axis2.Draw()

leg_nuebar = TLegend(0.75, 0.95, 0.75, 0.95)
leg_nuebar.AddEntry(simple_nuebar, "nuebar gSimple Flux",  "l")
leg_nuebar.AddEntry(rebin_nuebar,  "nuebar dk2nu Flux",    "l")
leg_nuebar.Draw()

c7.cd()
pad4 = TPad("pad4", "pad4", 0, 0.05, 1, 0.3)
pad4.SetTopMargin(0.1)
pad4.SetBottomMargin(0.1)
pad4.SetGridx()  # vertical grid
pad4.Draw()
pad4.cd()  # pad2 becomes the current pad

ratio_nuebar.SetLineColor(12)
ratio_nuebar.SetMinimum(0.0)  # Define Y ..
ratio_nuebar.SetMaximum(3.0)  # .. range
ratio_nuebar.GetXaxis().SetRangeUser(0, 4)
ratio_nuebar.Sumw2()
ratio_nuebar.SetStats(0)  # No statistics on lower plot
ratio_nuebar.SetMarkerStyle(21)
ratio_nuebar.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
rebin_nuebar.GetYaxis().SetTitleSize(20)
rebin_nuebar.GetYaxis().SetTitleFont(43)
rebin_nuebar.GetYaxis().SetTitleOffset(1.55)
ratio_nuebar.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_nuebar.GetYaxis().SetNdivisions(505)
ratio_nuebar.GetYaxis().SetTitleSize(14)
ratio_nuebar.GetYaxis().SetTitleFont(43)
ratio_nuebar.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_nuebar.GetYaxis().SetLabelFont(43)
ratio_nuebar.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_nuebar.GetXaxis().SetTitleSize(20)
ratio_nuebar.GetXaxis().SetTitleFont(43)
ratio_nuebar.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_nuebar.GetXaxis().SetLabelFont(43)
ratio_nuebar.GetXaxis().SetLabelSize(15)

# numu------------------------------------------------------------------

c8 = TCanvas("c8", "c8", 800, 800)
c8.cd()

pad5 = TPad("pad5", "pad5", 0, 0.3, 1, 1.0)
pad5.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad5.SetGridx()         # Vertical grid
pad5.Draw()             # Draw the upper pad: pad1
pad5.cd()               # pad1 becomes the current pad
rebin_numu.SetStats(0)   # No statistics on upper plot
simple_numu.SetStats(0)
rebin_numu.SetLineColor(38)
simple_numu.SetLineColor(46)
rebin_numu.SetLineWidth(6)
simple_numu.SetLineWidth(6)
simple_numu.GetXaxis().SetRangeUser(0, 4)
rebin_numu.GetXaxis().SetRangeUser(0, 4)

simple_numu.Draw("hist")
rebin_numu.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
rebin_numu.GetYaxis().SetLabelSize(0.)
axis3 = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis3.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis3.SetLabelSize(15)
axis3.Draw()

leg_numu = TLegend(0.75, 0.95, 0.75, 0.95)
leg_numu.AddEntry(simple_numu, "numu gSimple Flux",  "l")
leg_numu.AddEntry(rebin_numu,  "numu dk2nu Flux",    "l")
leg_numu.Draw()

c8.cd()
pad6 = TPad("pad6", "pad6", 0, 0.05, 1, 0.3)
pad6.SetTopMargin(0.1)
pad6.SetBottomMargin(0.1)
pad6.SetGridx()  # vertical grid
pad6.Draw()
pad6.cd()  # pad2 becomes the current pad

ratio_numu.SetLineColor(12)
ratio_numu.SetMinimum(0.0)  # Define Y ..
ratio_numu.SetMaximum(3.0)  # .. range
ratio_numu.GetXaxis().SetRangeUser(0, 4)
ratio_numu.Sumw2()
ratio_numu.SetStats(0)  # No statistics on lower plot
ratio_numu.SetMarkerStyle(21)
ratio_numu.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
rebin_numu.GetYaxis().SetTitleSize(20)
rebin_numu.GetYaxis().SetTitleFont(43)
rebin_numu.GetYaxis().SetTitleOffset(1.55)
ratio_numu.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_numu.GetYaxis().SetNdivisions(505)
ratio_numu.GetYaxis().SetTitleSize(14)
ratio_numu.GetYaxis().SetTitleFont(43)
ratio_numu.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_numu.GetYaxis().SetLabelFont(43)
ratio_numu.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_numu.GetXaxis().SetTitleSize(20)
ratio_numu.GetXaxis().SetTitleFont(43)
ratio_numu.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_numu.GetXaxis().SetLabelFont(43)
ratio_numu.GetXaxis().SetLabelSize(15)


# numu bar------------------------------------------------------------------

c9 = TCanvas("c9", "c9", 800, 800)
c9.cd()

pad7 = TPad("pad7", "pad7", 0, 0.3, 1, 1.0)
pad7.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad7.SetGridx()         # Vertical grid
pad7.Draw()             # Draw the upper pad: pad1
pad7.cd()               # pad1 becomes the current pad
rebin_numubar.SetStats(0)   # No statistics on upper plot
simple_numubar.SetStats(0)
rebin_numubar.SetLineColor(38)
simple_numubar.SetLineColor(46)
rebin_numubar.SetLineWidth(6)
simple_numubar.SetLineWidth(6)
simple_numubar.GetXaxis().SetRangeUser(0, 4)
rebin_numubar.GetXaxis().SetRangeUser(0, 4)

simple_numubar.Draw("hist")
rebin_numubar.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
rebin_numubar.GetYaxis().SetLabelSize(0.)
axis4 = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis4.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis4.SetLabelSize(15)
axis4.Draw()

leg_numubar = TLegend(0.75, 0.95, 0.75, 0.95)
leg_numubar.AddEntry(simple_numubar, "numubar gSimple Flux",  "l")
leg_numubar.AddEntry(rebin_numubar,  "numubar dk2nu Flux",    "l")
leg_numubar.Draw()

c9.cd()
pad8 = TPad("pad8", "pad8", 0, 0.05, 1, 0.3)
pad8.SetTopMargin(0.1)
pad8.SetBottomMargin(0.1)
pad8.SetGridx()  # vertical grid
pad8.Draw()
pad8.cd()  # pad2 becomes the current pad

ratio_numubar.SetLineColor(12)
ratio_numubar.SetMinimum(0.0)  # Define Y ..
ratio_numubar.SetMaximum(3.0)  # .. range
ratio_numubar.GetXaxis().SetRangeUser(0, 4)
ratio_numubar.Sumw2()
ratio_numubar.SetStats(0)  # No statistics on lower plot
ratio_numubar.SetMarkerStyle(21)
ratio_numubar.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
rebin_numubar.GetYaxis().SetTitleSize(20)
rebin_numubar.GetYaxis().SetTitleFont(43)
rebin_numubar.GetYaxis().SetTitleOffset(1.55)
ratio_numubar.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_numubar.GetYaxis().SetNdivisions(505)
ratio_numubar.GetYaxis().SetTitleSize(14)
ratio_numubar.GetYaxis().SetTitleFont(43)
ratio_numubar.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_numubar.GetYaxis().SetLabelFont(43)
ratio_numubar.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_numubar.GetXaxis().SetTitleSize(20)
ratio_numubar.GetXaxis().SetTitleFont(43)
ratio_numubar.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_numubar.GetXaxis().SetLabelFont(43)
ratio_numubar.GetXaxis().SetLabelSize(15)

wait = raw_input("PRESS ENTER TO CONTINUE.")
