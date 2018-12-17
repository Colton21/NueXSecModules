import sys
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TCanvas, gROOT, TPad, TGaxis, TColor, TLegend
import matplotlib

file1_path = "../arxiv/NuMIFlux_map_corrected.root"
file1 = TFile(file1_path, 'READ')
if(file1.IsOpen()):
    print 'File ', file1_path, ' is open'
if(file1.IsOpen() == False):
    quit()

# these are the corrected gSimple flux histograms
corrected_nue           = file1.Get("nueFluxHisto")
corrected_nuebar        = file1.Get("anueFluxHisto")
corrected_numu          = file1.Get("numuFluxHisto")
corrected_numubar       = file1.Get("anumuFluxHisto")
corrected_theta_nue     = file1.Get("nueFluxHistoTheta")
corrected_phi_nue       = file1.Get("nueFluxHistoPhi")
corrected_theta_nuebar  = file1.Get("anueFluxHistoTheta")
corrected_phi_nuebar    = file1.Get("anueFluxHistoPhi")
corrected_theta_numu    = file1.Get("numuFluxHistoTheta")
corrected_phi_numu      = file1.Get("numuFluxHistoPhi")
corrected_theta_numubar = file1.Get("anumuFluxHistoTheta")
corrected_phi_numubar   = file1.Get("anumuFluxHistoPhi")

#Also grabbing the 2D histograms

corrected_eng_theta_nue     = file1.Get("nueFluxHistoEngTheta")
corrected_eng_theta_nuebar  = file1.Get("anueFluxHistoEngTheta")
corrected_eng_theta_numu    = file1.Get("numuFluxHistoEngTheta")
corrected_eng_theta_numubar = file1.Get("anumuFluxHistoEngTheta")

corrected_eng_phi_nue     = file1.Get("nueFluxHistoEngPhi")
corrected_eng_phi_nuebar  = file1.Get("anueFluxHistoEngPhi")
corrected_eng_phi_numu    = file1.Get("numuFluxHistoEngPhi")
corrected_eng_phi_numubar = file1.Get("anumuFluxHistoEngPhi")

corrected_theta_phi_nue     = file1.Get("nueFluxHistoThetaPhi")
corrected_theta_phi_nuebar  = file1.Get("anueFluxHistoThetaPhi")
corrected_theta_phi_numu    = file1.Get("numuFluxHistoThetaPhi")
corrected_theta_phi_numubar = file1.Get("anumuFluxHistoThetaPhi")

file2_path = "../arxiv/NuMIFlux_map_old.root"
file2 = TFile(file2_path, 'READ')
if(file2.IsOpen()):
    print 'File ', file2_path, ' is open'
if(file2.IsOpen() == False):
    quit()

# these are the flux histograms coming from the gsimple files
simple_nue           = file2.Get("nueFluxHisto")
simple_nuebar        = file2.Get("anueFluxHisto")
simple_numu          = file2.Get("numuFluxHisto")
simple_numubar       = file2.Get("anumuFluxHisto")
simple_theta_nue     = file2.Get("nueFluxHistoTheta")
simple_phi_nue       = file2.Get("nueFluxHistoPhi")
simple_theta_nuebar  = file2.Get("anueFluxHistoTheta")
simple_phi_nuebar    = file2.Get("anueFluxHistoPhi")
simple_theta_numu    = file2.Get("numuFluxHistoTheta")
simple_phi_numu      = file2.Get("numuFluxHistoPhi")
simple_theta_numubar = file2.Get("anumuFluxHistoTheta")
simple_phi_numubar   = file2.Get("anumuFluxHistoPhi")

#2D histograms

simple_eng_theta_nue     = file1.Get("nueFluxHistoEngTheta")
simple_eng_theta_nuebar  = file1.Get("anueFluxHistoEngTheta")
simple_eng_theta_numu    = file1.Get("numuFluxHistoEngTheta")
simple_eng_theta_numubar = file1.Get("anumuFluxHistoEngTheta")

simple_eng_phi_nue     = file1.Get("nueFluxHistoEngPhi")
simple_eng_phi_nuebar  = file1.Get("anueFluxHistoEngPhi")
simple_eng_phi_numu    = file1.Get("numuFluxHistoEngPhi")
simple_eng_phi_numubar = file1.Get("anumuFluxHistoEngPhi")

simple_theta_phi_nue     = file1.Get("nueFluxHistoThetaPhi")
simple_theta_phi_nuebar  = file1.Get("anueFluxHistoThetaPhi")
simple_theta_phi_numu    = file1.Get("numuFluxHistoThetaPhi")
simple_theta_phi_numubar = file1.Get("anumuFluxHistoThetaPhi")

c1 = TCanvas("c1", "c1", 800, 600)
c1.cd()

corrected_nue.Draw("hist")
corrected_nuebar.Draw("hist same")

c2 = TCanvas("c2", "c2", 800, 600)
c2.cd()

simple_nue.Draw("hist")
simple_nuebar.Draw("hist same")

# we  want their ratios
ratio_nue = TH1D("ratio_nue", "ratio_nue", 400, 0, 20)
ratio_nuebar = TH1D("ratio_nuebar", "ratio_nuebar", 400, 0, 20)
ratio_numu = TH1D("ratio_nue", "ratio_nue", 400, 0, 20)
ratio_numubar = TH1D("ratio_nuebar", "ratio_nuebar", 400, 0, 20)

ratio_theta_nue     = TH1D("ratio_theta_nue", "ratio_theta_nue", 400, 0, 180)
ratio_phi_nue       = TH1D("ratio_phi_nue", "ratio_phi_nue", 400, -180, 180)
ratio_theta_nuebar  = TH1D("ratio_theta_nuebar", "ratio_theta_nuebar", 400, 0, 180)
ratio_phi_nuebar    = TH1D("ratio_phi_nuebar", "ratio_phi_nuebar", 400, -180, 180)
ratio_theta_numu    = TH1D("ratio_theta_numu", "ratio_theta_numu", 400, 0, 180)
ratio_phi_numu      = TH1D("ratio_phi_numu", "ratio_phi_numu", 400, -180, 180)
ratio_theta_numubar = TH1D("ratio_theta_numubar", "ratio_theta_numubar", 400, 0, 180)
ratio_phi_numubar   = TH1D("ratio_phi_numubar", "ratio_phi_numubar", 400, -180, 180)

ratio_eng_theta_nue     = TH2D("ratio_eng_theta_nue",      "ratio_eng_theta_nue",     400, 0, 20, 400, 0, 180)
ratio_eng_theta_nuebar  = TH2D("ratio_eng_theta_nuebar",   "ratio_eng_theta_nuebar",  400, 0, 20, 400, 0, 180)
ratio_eng_theta_numu    = TH2D("ratio_eng_theta_numu",     "ratio_eng_theta_numu",    400, 0, 20, 400, 0, 180)
ratio_eng_theta_numubar = TH2D("ratio_eng_theta_numubar",  "ratio_eng_theta_numubar", 400, 0, 20, 400, 0, 180)
ratio_eng_theta_total   = TH2D("ratio_eng_theta_total",    "ratio_eng_theta_total",   400, 0, 20, 400, 0, 180)

ratio_eng_phi_nue     = TH2D("ratio_eng_phi_nue",     "ratio_eng_phi_nue",     400, 0, 20, 400, -180, 180)
ratio_eng_phi_nuebar  = TH2D("ratio_eng_phi_nuebar",  "ratio_eng_phi_nuebar",  400, 0, 20, 400, -180, 180)
ratio_eng_phi_numu    = TH2D("ratio_eng_phi_numu",    "ratio_eng_phi_numu",    400, 0, 20, 400, -180, 180)
ratio_eng_phi_numubar = TH2D("ratio_eng_phi_numubar", "ratio_eng_phi_numubar", 400, 0, 20, 400, -180, 180)
ratio_eng_phi_total   = TH2D("ratio_eng_phi_total",   "ratio_eng_phi_total",   400, 0, 20, 400, -180, 180)

ratio_theta_phi_nue     = TH2D("ratio_theta_phi_nue",     "ratio_theta_phi_nue",     400, 0, 180, 400, -180, 180)
ratio_theta_phi_nuebar  = TH2D("ratio_theta_phi_nuebar",  "ratio_theta_phi_nuebar",  400, 0, 180, 400, -180, 180)
ratio_theta_phi_numu    = TH2D("ratio_theta_phi_numu",    "ratio_theta_phi_numu",    400, 0, 180, 400, -180, 180)
ratio_theta_phi_numubar = TH2D("ratio_theta_phi_numubar", "ratio_theta_phi_numubar", 400, 0, 180, 400, -180, 180)
ratio_theta_phi_total   = TH2D("ratio_theta_phi_total",   "ratio_theta_phi_total",   400, 0, 180, 400, -180, 180)

ratio_nue = corrected_nue.Clone()
ratio_nuebar = corrected_nuebar.Clone()
ratio_numu = corrected_numu.Clone()
ratio_numubar = corrected_numubar.Clone()

ratio_theta_nue     = corrected_theta_nue.Clone()
ratio_phi_nue       = corrected_phi_nue.Clone()
ratio_theta_nuebar  = corrected_theta_nuebar.Clone()
ratio_phi_nuebar    = corrected_phi_nuebar.Clone()
ratio_theta_numu    = corrected_theta_numu.Clone()
ratio_phi_numu      = corrected_phi_numu.Clone()
ratio_theta_numubar = corrected_theta_numubar.Clone()
ratio_phi_numubar   = corrected_phi_numubar.Clone()

ratio_eng_theta_nue     = corrected_eng_theta_nue.Clone()
ratio_eng_theta_nuebar  = corrected_eng_theta_nuebar.Clone()
ratio_eng_theta_numu    = corrected_eng_theta_numu.Clone()
ratio_eng_theta_numubar = corrected_eng_theta_numubar.Clone()

ratio_eng_phi_nue     = corrected_eng_phi_nue.Clone()
ratio_eng_phi_nuebar  = corrected_eng_phi_nuebar.Clone()
ratio_eng_phi_numu    = corrected_eng_phi_numu.Clone()
ratio_eng_phi_numubar = corrected_eng_phi_numubar.Clone()

ratio_theta_phi_nue     = corrected_theta_phi_nue.Clone()
ratio_theta_phi_nuebar  = corrected_theta_phi_nuebar.Clone()
ratio_theta_phi_numu    = corrected_theta_phi_numu.Clone()
ratio_theta_phi_numubar = corrected_theta_phi_numubar.Clone()


c4 = TCanvas("c4", "c4", 800, 600)
c4.cd()
ratio_nue.Divide(simple_nue)
ratio_nue.GetYaxis().SetTitle("corrected : std (#nu e)")
ratio_nue.Draw("ep")

c5 = TCanvas("c5", "c5", 800, 600)
c5.cd()
ratio_nuebar.Divide(simple_nuebar)
ratio_nuebar.GetYaxis().SetTitle("corrected : std (#nu e-bar)")
ratio_nuebar.Draw("ep")

c4b = TCanvas("c4b", "c4b", 800, 600)
c4b.cd()
ratio_numu.Divide(simple_numu)
ratio_numu.GetYaxis().SetTitle("corrected : std (#nu #mu)")
ratio_numu.Draw("ep")

c5b = TCanvas("c5b", "c5b", 800, 600)
c5b.cd()
ratio_numubar.Divide(simple_numubar)
ratio_numubar.GetYaxis().SetTitle("corrected : std (#nu #mu-bar)")
ratio_numubar.Draw("ep")

c4c = TCanvas("c4c", "c4c", 800, 600)
c4c.cd()
ratio_theta_nue.Divide(simple_theta_nue)
ratio_theta_nue.GetYaxis().SetTitle("corrected : std (#nu_{e})")
ratio_theta_nue.Draw("ep")

c4d = TCanvas("c4d", "c4d", 800, 600)
c4d.cd()
ratio_phi_nue.Divide(simple_phi_nue)
ratio_phi_nue.GetYaxis().SetTitle("corrected : std (#nu}_{e})")
ratio_phi_nue.Draw("ep")

c4e = TCanvas("c4e", "c4e", 800, 600)
c4e.cd()
ratio_theta_nuebar.Divide(simple_theta_nuebar)
ratio_theta_nuebar.GetYaxis().SetTitle("corrected : std (#bar{#nu}_{e})")
ratio_theta_nuebar.Draw("ep")

c4f = TCanvas("c4f", "c4f", 800, 600)
c4f.cd()
ratio_phi_nuebar.Divide(simple_phi_nuebar)
ratio_phi_nuebar.GetYaxis().SetTitle("corrected : std (#bar{#nu}_{e})")
ratio_phi_nuebar.Draw("ep")

c5c = TCanvas("c5c", "c5c", 800, 600)
c5c.cd()
ratio_theta_numu.Divide(simple_theta_numu)
ratio_theta_numu.GetYaxis().SetTitle("corrected : std (#nu_{#mu})")
ratio_theta_numu.Draw("ep")

c5d = TCanvas("c5d", "c5d", 800, 600)
c5d.cd()
ratio_phi_numu.Divide(simple_phi_numu)
ratio_phi_numu.GetYaxis().SetTitle("corrected : std (#nu}_{#mu})")
ratio_phi_numu.Draw("ep")

c5e = TCanvas("c5e", "c5e", 800, 600)
c5e.cd()
ratio_theta_numubar.Divide(simple_theta_numubar)
ratio_theta_numubar.GetYaxis().SetTitle("corrected : std (#bar{#nu}_{#mu})")
ratio_theta_numubar.Draw("ep")

c5f = TCanvas("c5f", "c5f", 800, 600)
c5f.cd()
ratio_phi_numubar.Divide(simple_phi_numubar)
ratio_phi_numubar.GetYaxis().SetTitle("corrected : std (#bar{#nu}_{#mu})")
ratio_phi_numubar.Draw("ep")

#2D ratio plots start here
c5g1 = TCanvas("c5g1", "c5g1", 800, 600)
c5g1.cd()
ratio_eng_theta_nue.Divide(simple_eng_theta_nue)
ratio_eng_theta_nue.SetTitle("corrected : std (#nu_{e})")
ratio_eng_theta_nue.GetXaxis().SetTitle("True Neutrino Energy [GeV]")
ratio_eng_theta_nue.GetYaxis().SetTitle("True Neutrino Theta [Degrees]"))
ratio_eng_theta_nue.Draw("colz")

c5g2 = TCanvas("c5g2", "c5g2", 800, 600)
c5g2.cd()
ratio_eng_phi_nue.Divide(simple_eng_phi_nue)
ratio_eng_phi_nue.SetTitle("corrected : std (#nu_{e})")
ratio_eng_phi_nue.GetXaxis().SetTitle("True Neutrino Energy[GeV]")
ratio_eng_phi_nue.GetYaxis().SetTitle("True Neutrino Phi [Degrees]"))
ratio_eng_phi_nue.Draw("colz")

c5g3 = TCanvas("c5g3", "c5g3", 800, 600)
c5g3.cd()
ratio_theta_phi_nue.Divide(simple_theta_phi_nue)
ratio_theta_phi_nue.SetTitle("corrected : std (#nu_{e})")
ratio_theta_phi_nue.GetXaxis().SetTitle("True Neutrino Theta [Degrees]")
ratio_theta_phi_nue.GetYaxis().SetTitle("True Neutrino Phi [Degrees]")
ratio_theta_phi_nue.Draw("colz")

c5h1 = TCanvas("c5h1", "c5h1", 800, 600)
c5h1.cd()
ratio_eng_theta_nuebar.Divide(simple_eng_theta_nuebar)
ratio_eng_theta_nuebar.SetTitle("corrected : std (#bar{#nu}_{e})")
ratio_eng_theta_nuebar.GetXaxis().SetTitle("True Neutrino Energy [GeV]")
ratio_eng_theta_nuebar.GetYaxis().SetTitle("True Neutrino Theta [Degrees]")
ratio_eng_theta_nuebar.Draw("colz")

c5h2 = TCanvas("c5h2", "c5h2", 800, 600)
c5h2.cd()
ratio_eng_phi_nuebar.Divide(simple_eng_phi_nuebar)
ratio_eng_phi_nuebar.SetTitle("corrected : std (#bar{#nu}_{e})")
ratio_eng_phi_nuebar.GetXaxis().SetTitle("True Neutrino Energy [GeV]")
ratio_eng_phi_nuebar.GetYaxis().SetTitle("True Neutrino Phi [Degrees]")
ratio_eng_phi_nuebar.Draw("colz")

c5h3 = TCanvas("c5h3", "c5h3", 800, 600)
c5h3.cd()
ratio_theta_phi_nuebar.Divide(simple_theta_phi_nuebar)
ratio_theta_phi_nuebar.SetTitle("corrected : std (#bar{#nu}_{e})")
ratio_theta_phi_nuebar.GetXaxis().SetTitle("True Neutrino Theta [Degrees]")
ratio_theta_phi_nuebar.GetYaxis().SetTitle("True Neutrino Phi [Degrees]")
ratio_theta_phi_nuebar.Draw("colz")

c5i1 = TCanvas("c5i1", "c5i1", 800, 600)
c5i1.cd()
ratio_eng_theta_numu.Divide(simple_eng_theta_numu)
ratio_eng_theta_numu.SetTitle("corrected : std (#nu_{#mu})")
ratio_eng_theta_numu.GetXaxis().SetTitle("True Neutrino Energy [GeV]")
ratio_eng_theta_numu.GetYaxis().SetTitle("True Neutrino Theta [Degrees]"))
ratio_eng_theta_numu.Draw("colz")

c5i2 = TCanvas("c5i2", "c5i2", 800, 600)
c5i2.cd()
ratio_eng_phi_numu.Divide(simple_eng_phi_numu)
ratio_eng_phi_numu.SetTitle("corrected : std (#nu_{#mu})")
ratio_eng_phi_numu.GetXaxis().SetTitle("True Neutrino Energy [GeV]")
ratio_eng_phi_numu.GetYaxis().SetTitle("True Neutrino Phi [Degrees]"))
ratio_eng_phi_numu.Draw("colz")

c5i3 = TCanvas("c5i3", "c5i3", 800, 600)
c5i3.cd()
ratio_theta_phi_numu.Divide(simple_theta_phi_numu)
ratio_theta_phi_numu.SetTitle("corrected : std (#nu_{#mu})")
ratio_theta_phi_numu.GetXaxis().SetTitle("True Neutrino Theta [Degrees]")
ratio_theta_phi_numu.GetYaxis().SetTitle()"True Neutrino Phi [Degrees]")
ratio_theta_phi_numu.Draw("colz")

c5j1 = TCanvas("c5j1", "c5j1", 800, 600)
c5j1.cd()
ratio_eng_theta_numubar.Divide(simple_eng_theta_numubar)
ratio_eng_theta_numubar.SetTitle("corrected : std (#bar{#nu}_{#mu})")
ratio_eng_theta_numubar.GetXaxis().SetTitle("True Neutrino Energy [GeV]")
ratio_eng_theta_numubar.GetYaxis().SetTitle("True Neutrino Theta [Degrees]")
ratio_eng_theta_numubar.Draw("colz")

c5j2 = TCanvas("c5j2", "c5j2", 800, 600)
c5j2.cd()
ratio_eng_phi_numubar.Divide(simple_eng_phi_numubar)
ratio_eng_phi_numubar.SetTitle("corrected : std (#bar{#nu}_{#mu})")
ratio_eng_phi_numubar.GetXaxis().SetTitle("True Neutrino Energy [GeV]")
ratio_eng_phi_numubar.GetYaxis().SetTitle("True Neutrino Phi [Degrees]")
ratio_eng_phi_numubar.Draw("colz")

c5j3 = TCanvas("c5j3", "c5j3", 800, 600)
c5j3.cd()
ratio_theta_phi_numubar.Divide(simple_theta_phi_numubar)
ratio_theta_phi_numubar.SetTitle("corrected : std (#bar{#nu}_{#mu})")
ratio_theta_phi_numubar.GetXaxis().SetTitle("True Neutrino Theta [Degrees]")
ratio_theta_phi_numubar.GetYaxis().SetTitle("True Neutrino Phi [Degrees]")
ratio_theta_phi_numubar.Draw("colz")

#2D ratio totals

c5k1 = TCanvas("c5k1", "c5k1", 800, 600)
c5k1.cd()
ratio_eng_theta_total.Add(ratio_eng_theta_nue,     1)
ratio_eng_theta_total.Add(ratio_eng_theta_nuebar,  1)
ratio_eng_theta_total.Add(ratio_eng_theta_numu,    1)
ratio_eng_theta_total.Add(ratio_eng_theta_numubar, 1)
ratio_eng_theta_total.SetTitle("corrected : std (total)")
ratio_eng_theta_total.GetXaxis().SetTitle("True Neutrino Energy [GeV]")
ratio_eng_theta_total.GetYaxis().SetTitle("True Neutrino Theta [Degrees]")
ratio_eng_theta_total.Draw("colz")

ratio_eng_phi_total.Add(ratio_eng_phi_nue,     1)
ratio_eng_phi_total.Add(ratio_eng_phi_nuebar,  1)
ratio_eng_phi_total.Add(ratio_eng_phi_numu,    1)
ratio_eng_phi_total.Add(ratio_eng_phi_numubar, 1)
ratio_eng_phi_total.SetTitle("corrected : std (total)")
ratio_eng_phi_total.GetXaxis().SetTitle("True Neutrino Energy [GeV]")
ratio_eng_phi_total.GetYaxis().SetTitle("True Neutirno Phi [Degrees]")
ratio_eng_phi_total.Draw("colz")

ratio_theta_phi_total.Add(ratio_theta_phi_nue,     1)
ratio_theta_phi_total.Add(ratio_theta_phi_nuebar,  1)
ratio_theta_phi_total.Add(ratio_theta_phi_numu,    1)
ratio_theta_phi_total.Add(ratio_theta_phi_numubar, 1)
ratio_theta_phi_total.SetTitle("corrected : std (total)")
ratio_theta_phi_total.GetXaxis().SetTitle("True Neutrino Theta [Degrees]")
ratio_theta_phi_total.GetYaxis().SetTitle("True Neutrino Phi [Degrees]")
ratio_theta_phi_total.Draw("colz")

#-----------------------------------------------------------------------

c6 = TCanvas("c6", "c6", 800, 800)
c6.cd()

pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad1.SetGridx()         # Vertical grid
pad1.Draw()             # Draw the upper pad: pad1
pad1.cd()               # pad1 becomes the current pad
corrected_nue.SetStats(0)   # No statistics on upper plot
simple_nue.SetStats(0)
corrected_nue.SetLineColor(38)
simple_nue.SetLineColor(46)
corrected_nue.SetLineWidth(4)
simple_nue.SetLineWidth(4)
simple_nue.GetXaxis().SetRangeUser(0, 4)
corrected_nue.GetXaxis().SetRangeUser(0, 4)

simple_nue.Draw("hist")
corrected_nue.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_nue.GetYaxis().SetLabelSize(0.)
axis = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis.SetLabelSize(15)
axis.Draw()

leg_nue = TLegend(0.75, 0.95, 0.75, 0.95)
leg_nue.AddEntry(simple_nue, "nue std. Flux",  "l")
leg_nue.AddEntry(corrected_nue,  "nue cor. Flux",    "l")
leg_nue.Draw()

c6.cd()
pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
pad2.SetTopMargin(0.1)
pad2.SetBottomMargin(0.1)
pad2.SetGridx()  # vertical grid
pad2.Draw()
pad2.cd()  # pad2 becomes the current pad

ratio_nue.SetLineColor(12)
ratio_nue.SetMinimum(0.5)  # Define Y ..
ratio_nue.SetMaximum(1.5)  # .. range
ratio_nue.GetXaxis().SetRangeUser(0, 4)
ratio_nue.Sumw2()
ratio_nue.SetStats(0)  # No statistics on lower plot
ratio_nue.SetMarkerStyle(3)
ratio_nue.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_nue.GetYaxis().SetTitleSize(20)
corrected_nue.GetYaxis().SetTitleFont(43)
corrected_nue.GetYaxis().SetTitleOffset(1.55)
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
corrected_nuebar.SetStats(0)   # No statistics on upper plot
simple_nuebar.SetStats(0)
corrected_nuebar.SetLineColor(38)
simple_nuebar.SetLineColor(46)
corrected_nuebar.SetLineWidth(4)
simple_nuebar.SetLineWidth(4)
simple_nuebar.GetXaxis().SetRangeUser(0, 4)
corrected_nuebar.GetXaxis().SetRangeUser(0, 4)

simple_nuebar.Draw("hist")
corrected_nuebar.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_nuebar.GetYaxis().SetLabelSize(0.)
axis2 = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis2.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis2.SetLabelSize(15)
axis2.Draw()

leg_nuebar = TLegend(0.75, 0.95, 0.75, 0.95)
leg_nuebar.AddEntry(simple_nuebar, "nuebar std. Flux",  "l")
leg_nuebar.AddEntry(corrected_nuebar,  "nuebar cor. Flux",    "l")
leg_nuebar.Draw()

c7.cd()
pad4 = TPad("pad4", "pad4", 0, 0.05, 1, 0.3)
pad4.SetTopMargin(0.1)
pad4.SetBottomMargin(0.1)
pad4.SetGridx()  # vertical grid
pad4.Draw()
pad4.cd()  # pad2 becomes the current pad

ratio_nuebar.SetLineColor(12)
ratio_nuebar.SetMinimum(0.5)  # Define Y ..
ratio_nuebar.SetMaximum(1.5)  # .. range
ratio_nuebar.GetXaxis().SetRangeUser(0, 4)
ratio_nuebar.Sumw2()
ratio_nuebar.SetStats(0)  # No statistics on lower plot
ratio_nuebar.SetMarkerStyle(3)
ratio_nuebar.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_nuebar.GetYaxis().SetTitleSize(20)
corrected_nuebar.GetYaxis().SetTitleFont(43)
corrected_nuebar.GetYaxis().SetTitleOffset(1.55)
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
corrected_numu.SetStats(0)   # No statistics on upper plot
simple_numu.SetStats(0)
corrected_numu.SetLineColor(38)
simple_numu.SetLineColor(46)
corrected_numu.SetLineWidth(4)
simple_numu.SetLineWidth(4)
simple_numu.GetXaxis().SetRangeUser(0, 4)
corrected_numu.GetXaxis().SetRangeUser(0, 4)

simple_numu.Draw("hist")
corrected_numu.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_numu.GetYaxis().SetLabelSize(0.)
axis3 = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis3.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis3.SetLabelSize(15)
axis3.Draw()

leg_numu = TLegend(0.75, 0.95, 0.75, 0.95)
leg_numu.AddEntry(simple_numu, "numu std. Flux",  "l")
leg_numu.AddEntry(corrected_numu,  "numu cor. Flux",    "l")
leg_numu.Draw()

c8.cd()
pad6 = TPad("pad6", "pad6", 0, 0.05, 1, 0.3)
pad6.SetTopMargin(0.1)
pad6.SetBottomMargin(0.1)
pad6.SetGridx()  # vertical grid
pad6.Draw()
pad6.cd()  # pad2 becomes the current pad

ratio_numu.SetLineColor(12)
ratio_numu.SetMinimum(0.5)  # Define Y ..
ratio_numu.SetMaximum(1.5)  # .. range
ratio_numu.GetXaxis().SetRangeUser(0, 4)
ratio_numu.Sumw2()
ratio_numu.SetStats(0)  # No statistics on lower plot
ratio_numu.SetMarkerStyle(3)
ratio_numu.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_numu.GetYaxis().SetTitleSize(20)
corrected_numu.GetYaxis().SetTitleFont(43)
corrected_numu.GetYaxis().SetTitleOffset(1.55)
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
corrected_numubar.SetStats(0)   # No statistics on upper plot
simple_numubar.SetStats(0)
corrected_numubar.SetLineColor(38)
simple_numubar.SetLineColor(46)
corrected_numubar.SetLineWidth(4)
simple_numubar.SetLineWidth(4)
simple_numubar.GetXaxis().SetRangeUser(0, 4)
corrected_numubar.GetXaxis().SetRangeUser(0, 4)

simple_numubar.Draw("hist")
corrected_numubar.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_numubar.GetYaxis().SetLabelSize(0.)
axis4 = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis4.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis4.SetLabelSize(15)
axis4.Draw()

leg_numubar = TLegend(0.75, 0.95, 0.75, 0.95)
leg_numubar.AddEntry(simple_numubar, "numubar std. Flux",  "l")
leg_numubar.AddEntry(corrected_numubar,  "numubar cor. Flux",    "l")
leg_numubar.Draw()

c9.cd()
pad8 = TPad("pad8", "pad8", 0, 0.05, 1, 0.3)
pad8.SetTopMargin(0.1)
pad8.SetBottomMargin(0.1)
pad8.SetGridx()  # vertical grid
pad8.Draw()
pad8.cd()  # pad2 becomes the current pad

ratio_numubar.SetLineColor(12)
ratio_numubar.SetMinimum(0.5)  # Define Y ..
ratio_numubar.SetMaximum(1.5)  # .. range
ratio_numubar.GetXaxis().SetRangeUser(0, 4)
ratio_numubar.Sumw2()
ratio_numubar.SetStats(0)  # No statistics on lower plot
ratio_numubar.SetMarkerStyle(3)
ratio_numubar.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_numubar.GetYaxis().SetTitleSize(20)
corrected_numubar.GetYaxis().SetTitleFont(43)
corrected_numubar.GetYaxis().SetTitleOffset(1.55)
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

#----------------------------------------------------
# Also replicate the overlay plots for theta and phi
#----------------------------------------------------

#-----------------------------------------------------------------------

c6b = TCanvas("c6b", "c6b", 800, 800)
c6b.cd()

pad1b = TPad("pad1b", "pad1b", 0, 0.3, 1, 1.0)
pad1b.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad1b.SetGridx()         # Vertical grid
pad1b.Draw()             # Draw the upper pad: pad1
pad1b.cd()               # pad1 becomes the current pad
corrected_theta_nue.SetStats(0)   # No statistics on upper plot
simple_theta_nue.SetStats(0)
corrected_theta_nue.SetLineColor(38)
simple_theta_nue.SetLineColor(46)
corrected_theta_nue.SetLineWidth(4)
simple_theta_nue.SetLineWidth(4)
simple_theta_nue.GetXaxis().SetRangeUser(0, 180)
corrected_theta_nue.GetXaxis().SetRangeUser(0, 180)

simple_theta_nue.Draw("hist")
corrected_theta_nue.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_theta_nue.GetYaxis().SetLabelSize(0.)
axisb = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axisb.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axisb.SetLabelSize(15)
axisb.Draw()

leg_nue_b = TLegend(0.75, 0.95, 0.75, 0.95)
leg_nue_b.AddEntry(simple_theta_nue, "nue std. Flux",  "l")
leg_nue_b.AddEntry(corrected_theta_nue,  "nue cor. Flux",    "l")
leg_nue_b.Draw()

c6b.cd()
pad2b = TPad("pad2b", "pad2b", 0, 0.05, 1, 0.3)
pad2b.SetTopMargin(0.1)
pad2b.SetBottomMargin(0.1)
pad2b.SetGridx()  # vertical grid
pad2b.Draw()
pad2b.cd()  # pad2 becomes the current pad

ratio_theta_nue.SetLineColor(12)
ratio_theta_nue.SetMinimum(0.5)  # Define Y ..
ratio_theta_nue.SetMaximum(1.5)  # .. range
ratio_theta_nue.GetXaxis().SetRangeUser(0, 180)
ratio_theta_nue.Sumw2()
ratio_theta_nue.SetStats(0)  # No statistics on lower plot
ratio_theta_nue.SetMarkerStyle(3)
ratio_theta_nue.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_theta_nue.GetYaxis().SetTitleSize(20)
corrected_theta_nue.GetYaxis().SetTitleFont(43)
corrected_theta_nue.GetYaxis().SetTitleOffset(1.55)
ratio_theta_nue.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_theta_nue.GetYaxis().SetNdivisions(505)
ratio_theta_nue.GetYaxis().SetTitleSize(14)
ratio_theta_nue.GetYaxis().SetTitleFont(43)
ratio_theta_nue.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_theta_nue.GetYaxis().SetLabelFont(43)
ratio_theta_nue.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_theta_nue.GetXaxis().SetTitleSize(20)
ratio_theta_nue.GetXaxis().SetTitleFont(43)
ratio_theta_nue.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_theta_nue.GetXaxis().SetLabelFont(43)
ratio_theta_nue.GetXaxis().SetLabelSize(15)

#------------------------------------------------------------------

c7b = TCanvas("c7b", "c7b", 800, 800)
c7b.cd()

pad3b = TPad("pad3b", "pad3b", 0, 0.3, 1, 1.0)
pad3b.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad3b.SetGridx()         # Vertical grid
pad3b.Draw()             # Draw the upper pad: pad1
pad3b.cd()               # pad1 becomes the current pad
corrected_theta_nuebar.SetStats(0)   # No statistics on upper plot
simple_theta_nuebar.SetStats(0)
corrected_theta_nuebar.SetLineColor(38)
simple_theta_nuebar.SetLineColor(46)
corrected_theta_nuebar.SetLineWidth(4)
simple_theta_nuebar.SetLineWidth(4)
simple_theta_nuebar.GetXaxis().SetRangeUser(0, 180)
corrected_theta_nuebar.GetXaxis().SetRangeUser(0, 180)

simple_theta_nuebar.Draw("hist")
corrected_theta_nuebar.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_theta_nuebar.GetYaxis().SetLabelSize(0.)
axis2b = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis2b.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis2b.SetLabelSize(15)
axis2b.Draw()

leg_nuebar_b = TLegend(0.75, 0.95, 0.75, 0.95)
leg_nuebar_b.AddEntry(simple_theta_nuebar, "nuebar std. Flux",  "l")
leg_nuebar_b.AddEntry(corrected_theta_nuebar,  "nuebar cor. Flux",    "l")
leg_nuebar_b.Draw()

c7b.cd()
pad4b = TPad("pad4b", "pad4b", 0, 0.05, 1, 0.3)
pad4b.SetTopMargin(0.1)
pad4b.SetBottomMargin(0.1)
pad4b.SetGridx()  # vertical grid
pad4b.Draw()
pad4b.cd()  # pad2 becomes the current pad

ratio_theta_nuebar.SetLineColor(12)
ratio_theta_nuebar.SetMinimum(0.5)  # Define Y ..
ratio_theta_nuebar.SetMaximum(1.5)  # .. range
ratio_theta_nuebar.GetXaxis().SetRangeUser(0, 180)
ratio_theta_nuebar.Sumw2()
ratio_theta_nuebar.SetStats(0)  # No statistics on lower plot
ratio_theta_nuebar.SetMarkerStyle(3)
ratio_theta_nuebar.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_theta_nuebar.GetYaxis().SetTitleSize(20)
corrected_theta_nuebar.GetYaxis().SetTitleFont(43)
corrected_theta_nuebar.GetYaxis().SetTitleOffset(1.55)
ratio_theta_nuebar.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_theta_nuebar.GetYaxis().SetNdivisions(505)
ratio_theta_nuebar.GetYaxis().SetTitleSize(14)
ratio_theta_nuebar.GetYaxis().SetTitleFont(43)
ratio_theta_nuebar.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_theta_nuebar.GetYaxis().SetLabelFont(43)
ratio_theta_nuebar.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_theta_nuebar.GetXaxis().SetTitleSize(20)
ratio_theta_nuebar.GetXaxis().SetTitleFont(43)
ratio_theta_nuebar.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_theta_nuebar.GetXaxis().SetLabelFont(43)
ratio_theta_nuebar.GetXaxis().SetLabelSize(15)

# numu------------------------------------------------------------------

c8b = TCanvas("c8b", "c8b", 800, 800)
c8b.cd()

pad5b = TPad("pad5b", "pad5b", 0, 0.3, 1, 1.0)
pad5b.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad5b.SetGridx()         # Vertical grid
pad5b.Draw()             # Draw the upper pad: pad1
pad5b.cd()               # pad1 becomes the current pad
corrected_theta_numu.SetStats(0)   # No statistics on upper plot
simple_theta_numu.SetStats(0)
corrected_theta_numu.SetLineColor(38)
simple_theta_numu.SetLineColor(46)
corrected_theta_numu.SetLineWidth(4)
simple_theta_numu.SetLineWidth(4)
simple_theta_numu.GetXaxis().SetRangeUser(0, 180)
corrected_theta_numu.GetXaxis().SetRangeUser(0, 180)

simple_theta_numu.Draw("hist")
corrected_theta_numu.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_theta_numu.GetYaxis().SetLabelSize(0.)
axis3b = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis3b.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis3b.SetLabelSize(15)
axis3b.Draw()

leg_numu_b = TLegend(0.75, 0.95, 0.75, 0.95)
leg_numu_b.AddEntry(simple_theta_numu, "numu std. Flux",  "l")
leg_numu_b.AddEntry(corrected_theta_numu,  "numu cor. Flux",    "l")
leg_numu_b.Draw()

c8b.cd()
pad6b = TPad("pad6b", "pad6b", 0, 0.05, 1, 0.3)
pad6b.SetTopMargin(0.1)
pad6b.SetBottomMargin(0.1)
pad6b.SetGridx()  # vertical grid
pad6b.Draw()
pad6b.cd()  # pad2 becomes the current pad

ratio_theta_numu.SetLineColor(12)
ratio_theta_numu.SetMinimum(0.5)  # Define Y ..
ratio_theta_numu.SetMaximum(1.5)  # .. range
ratio_theta_numu.GetXaxis().SetRangeUser(0, 180)
ratio_theta_numu.Sumw2()
ratio_theta_numu.SetStats(0)  # No statistics on lower plot
ratio_theta_numu.SetMarkerStyle(3)
ratio_theta_numu.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_theta_numu.GetYaxis().SetTitleSize(20)
corrected_theta_numu.GetYaxis().SetTitleFont(43)
corrected_theta_numu.GetYaxis().SetTitleOffset(1.55)
ratio_theta_numu.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_theta_numu.GetYaxis().SetNdivisions(505)
ratio_theta_numu.GetYaxis().SetTitleSize(14)
ratio_theta_numu.GetYaxis().SetTitleFont(43)
ratio_theta_numu.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_theta_numu.GetYaxis().SetLabelFont(43)
ratio_theta_numu.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_theta_numu.GetXaxis().SetTitleSize(20)
ratio_theta_numu.GetXaxis().SetTitleFont(43)
ratio_theta_numu.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_theta_numu.GetXaxis().SetLabelFont(43)
ratio_theta_numu.GetXaxis().SetLabelSize(15)


# numu bar------------------------------------------------------------------

c9b = TCanvas("c9b", "c9b", 800, 800)
c9b.cd()

pad7b = TPad("pad7b", "pad7b", 0, 0.3, 1, 1.0)
pad7b.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad7b.SetGridx()         # Vertical grid
pad7b.Draw()             # Draw the upper pad: pad1
pad7b.cd()               # pad1 becomes the current pad
corrected_theta_numubar.SetStats(0)   # No statistics on upper plot
simple_theta_numubar.SetStats(0)
corrected_theta_numubar.SetLineColor(38)
simple_theta_numubar.SetLineColor(46)
corrected_theta_numubar.SetLineWidth(4)
simple_theta_numubar.SetLineWidth(4)
simple_theta_numubar.GetXaxis().SetRangeUser(0, 180)
corrected_theta_numubar.GetXaxis().SetRangeUser(0, 180)

simple_theta_numubar.Draw("hist")
corrected_theta_numubar.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_theta_numubar.GetYaxis().SetLabelSize(0.)
axis4b = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis4b.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis4b.SetLabelSize(15)
axis4b.Draw()

leg_numubar_b = TLegend(0.75, 0.95, 0.75, 0.95)
leg_numubar_b.AddEntry(simple_theta_numubar, "numubar std. Flux",  "l")
leg_numubar_b.AddEntry(corrected_theta_numubar,  "numubar cor. Flux",    "l")
leg_numubar_b.Draw()

c9b.cd()
pad8b = TPad("pad8b", "pad8b", 0, 0.05, 1, 0.3)
pad8b.SetTopMargin(0.1)
pad8b.SetBottomMargin(0.1)
pad8b.SetGridx()  # vertical grid
pad8b.Draw()
pad8b.cd()  # pad2 becomes the current pad

ratio_theta_numubar.SetLineColor(12)
ratio_theta_numubar.SetMinimum(0.5)  # Define Y ..
ratio_theta_numubar.SetMaximum(1.5)  # .. range
ratio_theta_numubar.GetXaxis().SetRangeUser(0, 180)
ratio_theta_numubar.Sumw2()
ratio_theta_numubar.SetStats(0)  # No statistics on lower plot
ratio_theta_numubar.SetMarkerStyle(3)
ratio_theta_numubar.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_theta_numubar.GetYaxis().SetTitleSize(20)
corrected_theta_numubar.GetYaxis().SetTitleFont(43)
corrected_theta_numubar.GetYaxis().SetTitleOffset(1.55)
ratio_theta_numubar.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_theta_numubar.GetYaxis().SetNdivisions(505)
ratio_theta_numubar.GetYaxis().SetTitleSize(14)
ratio_theta_numubar.GetYaxis().SetTitleFont(43)
ratio_theta_numubar.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_theta_numubar.GetYaxis().SetLabelFont(43)
ratio_theta_numubar.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_theta_numubar.GetXaxis().SetTitleSize(20)
ratio_theta_numubar.GetXaxis().SetTitleFont(43)
ratio_theta_numubar.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_theta_numubar.GetXaxis().SetLabelFont(43)
ratio_theta_numubar.GetXaxis().SetLabelSize(15)


#-----------------------------------------------------------------------

c6c = TCanvas("c6c", "c6c", 800, 800)
c6c.cd()

pad1c = TPad("pad1c", "pad1c", 0, 0.3, 1, 1.0)
pad1c.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad1c.SetGridx()         # Vertical grid
pad1c.Draw()             # Draw the upper pad: pad1
pad1c.cd()               # pad1 becomes the current pad
corrected_phi_nue.SetStats(0)   # No statistics on upper plot
simple_phi_nue.SetStats(0)
corrected_phi_nue.SetLineColor(38)
simple_phi_nue.SetLineColor(46)
corrected_phi_nue.SetLineWidth(4)
simple_phi_nue.SetLineWidth(4)
simple_phi_nue.GetXaxis().SetRangeUser(-180, 180)
corrected_phi_nue.GetXaxis().SetRangeUser(-180, 180)

simple_phi_nue.Draw("hist")
corrected_phi_nue.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_phi_nue.GetYaxis().SetLabelSize(0.)
axisc = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axisc.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axisc.SetLabelSize(15)
axisc.Draw()

leg_nue_c = TLegend(0.75, 0.95, 0.75, 0.95)
leg_nue_c.AddEntry(simple_phi_nue, "nue std. Flux",  "l")
leg_nue_c.AddEntry(corrected_phi_nue,  "nue cor. Flux",    "l")
leg_nue_c.Draw()

c6c.cd()
pad2c = TPad("pad2c", "pad2c", 0, 0.05, 1, 0.3)
pad2c.SetTopMargin(0.1)
pad2c.SetBottomMargin(0.1)
pad2c.SetGridx()  # vertical grid
pad2c.Draw()
pad2c.cd()  # pad2 becomes the current pad

ratio_phi_nue.SetLineColor(12)
ratio_phi_nue.SetMinimum(0.5)  # Define Y ..
ratio_phi_nue.SetMaximum(1.5)  # .. range
ratio_phi_nue.GetXaxis().SetRangeUser(-180, 180)
ratio_phi_nue.Sumw2()
ratio_phi_nue.SetStats(0)  # No statistics on lower plot
ratio_phi_nue.SetMarkerStyle(3)
ratio_phi_nue.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_phi_nue.GetYaxis().SetTitleSize(20)
corrected_phi_nue.GetYaxis().SetTitleFont(43)
corrected_phi_nue.GetYaxis().SetTitleOffset(1.55)
ratio_phi_nue.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_phi_nue.GetYaxis().SetNdivisions(505)
ratio_phi_nue.GetYaxis().SetTitleSize(14)
ratio_phi_nue.GetYaxis().SetTitleFont(43)
ratio_phi_nue.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_phi_nue.GetYaxis().SetLabelFont(43)
ratio_phi_nue.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_phi_nue.GetXaxis().SetTitleSize(20)
ratio_phi_nue.GetXaxis().SetTitleFont(43)
ratio_phi_nue.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_phi_nue.GetXaxis().SetLabelFont(43)
ratio_phi_nue.GetXaxis().SetLabelSize(15)

#------------------------------------------------------------------

c7c = TCanvas("c7c", "c7c", 800, 800)
c7c.cd()

pad3c = TPad("pad3c", "pad3c", 0, 0.3, 1, 1.0)
pad3c.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad3c.SetGridx()         # Vertical grid
pad3c.Draw()             # Draw the upper pad: pad1
pad3c.cd()               # pad1 becomes the current pad
corrected_phi_nuebar.SetStats(0)   # No statistics on upper plot
simple_phi_nuebar.SetStats(0)
corrected_phi_nuebar.SetLineColor(38)
simple_phi_nuebar.SetLineColor(46)
corrected_phi_nuebar.SetLineWidth(4)
simple_phi_nuebar.SetLineWidth(4)
simple_phi_nuebar.GetXaxis().SetRangeUser(-180, 180)
corrected_phi_nuebar.GetXaxis().SetRangeUser(-180, 180)

simple_phi_nuebar.Draw("hist")
corrected_phi_nuebar.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_phi_nuebar.GetYaxis().SetLabelSize(0.)
axis2c = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis2c.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis2c.SetLabelSize(15)
axis2c.Draw()

leg_nuebar_c = TLegend(0.75, 0.95, 0.75, 0.95)
leg_nuebar_c.AddEntry(simple_phi_nuebar, "nuebar std. Flux",  "l")
leg_nuebar_c.AddEntry(corrected_phi_nuebar,  "nuebar cor. Flux",    "l")
leg_nuebar_c.Draw()

c7c.cd()
pad4c = TPad("pad4c", "pad4c", 0, 0.05, 1, 0.3)
pad4c.SetTopMargin(0.1)
pad4c.SetBottomMargin(0.1)
pad4c.SetGridx()  # vertical grid
pad4c.Draw()
pad4c.cd()  # pad2 becomes the current pad

ratio_phi_nuebar.SetLineColor(12)
ratio_phi_nuebar.SetMinimum(0.5)  # Define Y ..
ratio_phi_nuebar.SetMaximum(1.5)  # .. range
ratio_phi_nuebar.GetXaxis().SetRangeUser(-180, 180)
ratio_phi_nuebar.Sumw2()
ratio_phi_nuebar.SetStats(0)  # No statistics on lower plot
ratio_phi_nuebar.SetMarkerStyle(3)
ratio_phi_nuebar.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_phi_nuebar.GetYaxis().SetTitleSize(20)
corrected_phi_nuebar.GetYaxis().SetTitleFont(43)
corrected_phi_nuebar.GetYaxis().SetTitleOffset(1.55)
ratio_phi_nuebar.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_phi_nuebar.GetYaxis().SetNdivisions(505)
ratio_phi_nuebar.GetYaxis().SetTitleSize(14)
ratio_phi_nuebar.GetYaxis().SetTitleFont(43)
ratio_phi_nuebar.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_phi_nuebar.GetYaxis().SetLabelFont(43)
ratio_phi_nuebar.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_phi_nuebar.GetXaxis().SetTitleSize(20)
ratio_phi_nuebar.GetXaxis().SetTitleFont(43)
ratio_phi_nuebar.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_phi_nuebar.GetXaxis().SetLabelFont(43)
ratio_phi_nuebar.GetXaxis().SetLabelSize(15)

# numu------------------------------------------------------------------

c8c = TCanvas("c8c", "c8c", 800, 800)
c8c.cd()

pad5c = TPad("pad5c", "pad5c", 0, 0.3, 1, 1.0)
pad5c.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad5c.SetGridx()         # Vertical grid
pad5c.Draw()             # Draw the upper pad: pad1
pad5c.cd()               # pad1 becomes the current pad
corrected_phi_numu.SetStats(0)   # No statistics on upper plot
simple_phi_numu.SetStats(0)
corrected_phi_numu.SetLineColor(38)
simple_phi_numu.SetLineColor(46)
corrected_phi_numu.SetLineWidth(4)
simple_phi_numu.SetLineWidth(4)
simple_phi_numu.GetXaxis().SetRangeUser(-180, 180)
corrected_phi_numu.GetXaxis().SetRangeUser(-180, 180)

simple_phi_numu.Draw("hist")
corrected_phi_numu.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_phi_numu.GetYaxis().SetLabelSize(0.)
axis3c = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis3c.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis3c.SetLabelSize(15)
axis3c.Draw()

leg_numu_c = TLegend(0.75, 0.95, 0.75, 0.95)
leg_numu_c.AddEntry(simple_phi_numu, "numu std. Flux",  "l")
leg_numu_c.AddEntry(corrected_phi_numu,  "numu cor. Flux",    "l")
leg_numu_c.Draw()

c8c.cd()
pad6c = TPad("pad6c", "pad6c", 0, 0.05, 1, 0.3)
pad6c.SetTopMargin(0.1)
pad6c.SetBottomMargin(0.1)
pad6c.SetGridx()  # vertical grid
pad6c.Draw()
pad6c.cd()  # pad2 becomes the current pad

ratio_phi_numu.SetLineColor(12)
ratio_phi_numu.SetMinimum(0.5)  # Define Y ..
ratio_phi_numu.SetMaximum(1.5)  # .. range
ratio_phi_numu.GetXaxis().SetRangeUser(-180, 180)
ratio_phi_numu.Sumw2()
ratio_phi_numu.SetStats(0)  # No statistics on lower plot
ratio_phi_numu.SetMarkerStyle(3)
ratio_phi_numu.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_phi_numu.GetYaxis().SetTitleSize(20)
corrected_phi_numu.GetYaxis().SetTitleFont(43)
corrected_phi_numu.GetYaxis().SetTitleOffset(1.55)
ratio_phi_numu.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_phi_numu.GetYaxis().SetNdivisions(505)
ratio_phi_numu.GetYaxis().SetTitleSize(14)
ratio_phi_numu.GetYaxis().SetTitleFont(43)
ratio_phi_numu.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_phi_numu.GetYaxis().SetLabelFont(43)
ratio_phi_numu.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_phi_numu.GetXaxis().SetTitleSize(20)
ratio_phi_numu.GetXaxis().SetTitleFont(43)
ratio_phi_numu.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_phi_numu.GetXaxis().SetLabelFont(43)
ratio_phi_numu.GetXaxis().SetLabelSize(15)


# numu bar------------------------------------------------------------------

c9b = TCanvas("c9b", "c9b", 800, 800)
c9b.cd()

pad7c = TPad("pad7c", "pad7c", 0, 0.3, 1, 1.0)
pad7c.SetBottomMargin(0.1)  # Upper and lower plot are joined
pad7c.SetGridx()         # Vertical grid
pad7c.Draw()             # Draw the upper pad: pad1
pad7c.cd()               # pad1 becomes the current pad
corrected_phi_numubar.SetStats(0)   # No statistics on upper plot
simple_phi_numubar.SetStats(0)
corrected_phi_numubar.SetLineColor(38)
simple_phi_numubar.SetLineColor(46)
corrected_phi_numubar.SetLineWidth(4)
simple_phi_numubar.SetLineWidth(4)
simple_phi_numubar.GetXaxis().SetRangeUser(0, 180)
corrected_phi_numubar.GetXaxis().SetRangeUser(0, 180)

simple_phi_numubar.Draw("hist")
corrected_phi_numubar.Draw("hist same")


# Do not draw the Y axis label on the upper plot and redraw a small
# axis instead, in order to avoid the first label (0) to be clipped.
corrected_phi_numubar.GetYaxis().SetLabelSize(0.)
axis4c = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis4c.SetLabelFont(43)  # Absolute font size in pixel (precision 3)
axis4c.SetLabelSize(15)
axis4c.Draw()

leg_numubar_c = TLegend(0.75, 0.95, 0.75, 0.95)
leg_numubar_c.AddEntry(simple_phi_numubar, "numubar std. Flux",  "l")
leg_numubar_c.AddEntry(corrected_phi_numubar,  "numubar cor. Flux",    "l")
leg_numubar_c.Draw()

c9c.cd()
pad8c = TPad("pad8c", "pad8c", 0, 0.05, 1, 0.3)
pad8c.SetTopMargin(0.1)
pad8c.SetBottomMargin(0.1)
pad8c.SetGridx()  # vertical grid
pad8c.Draw()
pad8c.cd()  # pad2 becomes the current pad

ratio_phi_numubar.SetLineColor(12)
ratio_phi_numubar.SetMinimum(0.5)  # Define Y ..
ratio_phi_numubar.SetMaximum(1.5)  # .. range
ratio_phi_numubar.GetXaxis().SetRangeUser(-180, 180)
ratio_phi_numubar.Sumw2()
ratio_phi_numubar.SetStats(0)  # No statistics on lower plot
ratio_phi_numubar.SetMarkerStyle(3)
ratio_phi_numubar.Draw("hist p")  # Draw the ratio plot

# Y axis h1 plot settings
corrected_phi_numubar.GetYaxis().SetTitleSize(20)
corrected_phi_numubar.GetYaxis().SetTitleFont(43)
corrected_phi_numubar.GetYaxis().SetTitleOffset(1.55)
ratio_phi_numubar.SetTitle("")  # Remove the ratio title

# Y axis ratio plot settings
ratio_phi_numubar.GetYaxis().SetNdivisions(505)
ratio_phi_numubar.GetYaxis().SetTitleSize(14)
ratio_phi_numubar.GetYaxis().SetTitleFont(43)
ratio_phi_numubar.GetYaxis().SetTitleOffset(1.55)
# Absolute font size in pixel (precision 3)
ratio_phi_numubar.GetYaxis().SetLabelFont(43)
ratio_phi_numubar.GetYaxis().SetLabelSize(15)

# X axis ratio plot settings
ratio_phi_numubar.GetXaxis().SetTitleSize(20)
ratio_phi_numubar.GetXaxis().SetTitleFont(43)
ratio_phi_numubar.GetXaxis().SetTitleOffset(4.)
# Absolute font size in pixel (precision 3)
ratio_phi_numubar.GetXaxis().SetLabelFont(43)
ratio_phi_numubar.GetXaxis().SetLabelSize(15)

wait = raw_input("PRESS ENTER TO CONTINUE.")
