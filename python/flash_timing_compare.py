import sys
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TCanvas, gROOT, TPad, TGaxis, TColor, TLegend, TLine


def median(lst):
    n = len(lst)
    if n < 1:
        return None
    if n % 2 == 1:
        return sorted(lst)[n // 2]
    else:
        return sum(sorted(lst)[n // 2 - 1:n // 2 + 1]) / 2.0

print "**************************************************"
print "Expecting Samples 1-4 of EXT and then On-Beam Data"
print "**************************************************"

file1_path = sys.argv[1]
file1 = TFile(file1_path, 'READ')
if(file1.IsOpen()):
    print 'File ', file1_path, ' is open'
if(file1.IsOpen() == False):
    quit()
dir1 = file1.AnalyzeTPCO
opt_tree1 = dir1.Get("optical_tree")
num_events1 = opt_tree1.GetEntries()
print "Number of Events: ", num_events1


file2_path = sys.argv[2]
file2 = TFile(file2_path, 'READ')
if(file2.IsOpen()):
    print 'File ', file2_path, ' is open'
if(file2.IsOpen() == False):
    quit()
dir2 = file2.AnalyzeTPCO
opt_tree2 = dir2.Get("optical_tree")
num_events2 = opt_tree2.GetEntries()
print "Number of Events: ", num_events2


file3_path = sys.argv[3]
file3 = TFile(file3_path, 'READ')
if(file3.IsOpen()):
    print 'File ', file3_path, ' is open'
if(file3.IsOpen() == False):
    quit()
dir3 = file3.AnalyzeTPCO
opt_tree3 = dir3.Get("optical_tree")
num_events3 = opt_tree3.GetEntries()
print "Number of Events: ", num_events3


file4_path = sys.argv[4]
file4 = TFile(file4_path, 'READ')
if(file4.IsOpen()):
    print 'File ', file4_path, ' is open'
if(file4.IsOpen() == False):
    quit()
dir4 = file4.AnalyzeTPCO
opt_tree4 = dir4.Get("optical_tree")
num_events4 = opt_tree4.GetEntries()
print "Number of Events: ", num_events4

h_timing_1 = TH1D("h_timing_1", "h_timing_1", 60, 0, 20)
h_timing_2 = TH1D("h_timing_2", "h_timing_2", 60, 0, 20)
h_timing_3 = TH1D("h_timing_3", "h_timing_3", 60, 0, 20)
h_timing_4 = TH1D("h_timing_4", "h_timing_4", 60, 0, 20)

h_timing_1_2 = TH1D("h_timing_1_2", "h_timing_1_2", 60, 0, 20)
h_timing_2_2 = TH1D("h_timing_2_2", "h_timing_2_2", 60, 0, 20)
h_timing_3_2 = TH1D("h_timing_3_2", "h_timing_3_2", 60, 0, 20)
h_timing_4_2 = TH1D("h_timing_4_2", "h_timing_4_2", 60, 0, 20)

h_pe_1 = TH1D("h_pe_1", "h_pe_1", 80, 0, 2000)
h_pe_2 = TH1D("h_pe_2", "h_pe_2", 80, 0, 2000)
h_pe_3 = TH1D("h_pe_3", "h_pe_3", 80, 0, 2000)
h_pe_4 = TH1D("h_pe_4", "h_pe_4", 80, 0, 2000)

h_pe_1_2 = TH1D("h_pe_1_2", "h_pe_1_2", 80, 0, 2000)
h_pe_2_2 = TH1D("h_pe_2_2", "h_pe_2_2", 80, 0, 2000)
h_pe_3_2 = TH1D("h_pe_3_2", "h_pe_3_2", 80, 0, 2000)
h_pe_4_2 = TH1D("h_pe_4_2", "h_pe_4_2", 80, 0, 2000)

h_timing_pe_1 = TH2D("h_timing_pe_1", "h_timing_pe_1", 80, 0, 20, 80, 0, 2000)
h_timing_pe_2 = TH2D("h_timing_pe_2", "h_timing_pe_2", 80, 0, 20, 80, 0, 2000)
h_timing_pe_3 = TH2D("h_timing_pe_3", "h_timing_pe_3", 80, 0, 20, 80, 0, 2000)
h_timing_pe_4 = TH2D("h_timing_pe_4", "h_timing_pe_4", 80, 0, 20, 80, 0, 2000)

h_timing_pe_1_2 = TH2D(
    "h_timing_pe_1_2", "h_timing_pe_1_2", 80, 0, 20, 80, 0, 2000)
h_timing_pe_2_2 = TH2D(
    "h_timing_pe_2_2", "h_timing_pe_2_2", 80, 0, 20, 80, 0, 2000)
h_timing_pe_3_2 = TH2D(
    "h_timing_pe_3_2", "h_timing_pe_3_2", 80, 0, 20, 80, 0, 2000)
h_timing_pe_4_2 = TH2D(
    "h_timing_pe_4_2", "h_timing_pe_4_2", 80, 0, 20, 80, 0, 2000)

# differences in trigger board configurations
ext_off_set_factor = 0.343
#ext_off_set_factor = 0.406

for num in range(num_events1):
    opt_tree1.GetEntry(num)
    flash_time = opt_tree1.OpFlashTime
    flash_pe = opt_tree1.OpFlashPE
    h_timing_1.Fill(flash_time - ext_off_set_factor)
    h_timing_1_2.Fill(flash_time - ext_off_set_factor)
    h_pe_1.Fill(flash_pe)
    h_pe_1_2.Fill(flash_pe)
    h_timing_pe_1.Fill(flash_time, flash_pe)
    h_timing_pe_1_2.Fill(flash_time, flash_pe)

for num in range(num_events2):
    opt_tree2.GetEntry(num)
    flash_time = opt_tree2.OpFlashTime
    flash_pe = opt_tree2.OpFlashPE
    h_timing_2.Fill(flash_time - ext_off_set_factor)
    h_timing_2_2.Fill(flash_time - ext_off_set_factor)
    h_pe_2.Fill(flash_pe)
    h_pe_2_2.Fill(flash_pe)
    h_timing_pe_2.Fill(flash_time, flash_pe)
    h_timing_pe_2_2.Fill(flash_time, flash_pe)

for num in range(num_events3):
    opt_tree3.GetEntry(num)
    flash_time = opt_tree3.OpFlashTime
    flash_pe = opt_tree3.OpFlashPE
    h_timing_3.Fill(flash_time - ext_off_set_factor)
    h_timing_3_2.Fill(flash_time - ext_off_set_factor)
    h_pe_3.Fill(flash_pe)
    h_pe_3_2.Fill(flash_pe)
    h_timing_pe_3.Fill(flash_time, flash_pe)
    h_timing_pe_3_2.Fill(flash_time, flash_pe)

for num in range(num_events4):
    opt_tree4.GetEntry(num)
    flash_time = opt_tree4.OpFlashTime
    flash_pe = opt_tree4.OpFlashPE
    h_timing_4.Fill(flash_time - ext_off_set_factor)
    h_timing_4_2.Fill(flash_time - ext_off_set_factor)
    h_pe_4.Fill(flash_pe)
    h_pe_4_2.Fill(flash_pe)
    h_timing_pe_4.Fill(flash_time, flash_pe)
    h_timing_pe_4_2.Fill(flash_time, flash_pe)

intime_scale_factor_1 = 6180310. / 1141794.9
intime_scale_factor_2 = 6180310. / 1619476.0
intime_scale_factor_3 = 6180310. / 1937816.9
intime_scale_factor_4 = 6180310. / 588091.2
intime_scale_factor_234 = 6180310. / (1619476.0 + 1937816.9 + 588091.2)
intime_scale_factor_total = 6180310. / \
    (1141794.9 + 1619476.0 + 1937816.9 + 588091.2)

print "Scale Factors: "
print '\t', "1: ", intime_scale_factor_1
print '\t', "2: ", intime_scale_factor_2
print '\t', "3: ", intime_scale_factor_3
print '\t', "4: ", intime_scale_factor_4
print '\t', "Total: ", intime_scale_factor_total

line = TLine(0, 1, 20, 1)
line.SetLineColor(2)

intime_scale_factor_12 = intime_scale_factor_1 / intime_scale_factor_2
intime_scale_factor_13 = intime_scale_factor_1 / intime_scale_factor_3
intime_scale_factor_14 = intime_scale_factor_1 / intime_scale_factor_4
intime_scale_factor_23 = intime_scale_factor_2 / intime_scale_factor_3
intime_scale_factor_24 = intime_scale_factor_2 / intime_scale_factor_4
intime_scale_factor_34 = intime_scale_factor_3 / intime_scale_factor_4

h_timing_1.Sumw2()
h_timing_2.Sumw2()
h_timing_3.Sumw2()
h_timing_4.Sumw2()
h_timing_1_2.Sumw2()
h_timing_2_2.Sumw2()
h_timing_3_2.Sumw2()
h_timing_4_2.Sumw2()

h_divide_12 = h_timing_1.Clone()
h_divide_13 = h_timing_1.Clone()
h_divide_14 = h_timing_1.Clone()
h_divide_23 = h_timing_2.Clone()
h_divide_24 = h_timing_2.Clone()
h_divide_34 = h_timing_3.Clone()

h_divide_12.Divide(h_timing_2)
h_divide_12.Scale(intime_scale_factor_12)

h_divide_13.Divide(h_timing_3)
h_divide_13.Scale(intime_scale_factor_13)

h_divide_14.Divide(h_timing_4)
h_divide_14.Scale(intime_scale_factor_14)

h_divide_23.Divide(h_timing_3)
h_divide_23.Scale(intime_scale_factor_23)

h_divide_24.Divide(h_timing_4)
h_divide_24.Scale(intime_scale_factor_24)

h_divide_34.Divide(h_timing_4)
h_divide_34.Scale(intime_scale_factor_34)


c1a = TCanvas()
c1a.cd()
h_timing_1.Draw("e")
c1a.Print("plots/ext_timing_1.pdf")

c1b = TCanvas()
c1b.cd()
h_timing_2.Draw("e")
c1b.Print("plots/ext_timing_2.pdf")

c1c = TCanvas()
c1c.cd()
h_timing_3.Draw("e")
c1c.Print("plots/ext_timing_3.pdf")

c1d = TCanvas()
c1d.cd()
h_timing_4.Draw("e")
c1d.Print("plots/ext_timing_4.pdf")

c2a = TCanvas()
c2a.cd()
h_divide_12.Draw("e")
line.Draw("same")
c2a.Print("plots/ext_timing_divide_12.pdf")

c2b = TCanvas()
c2b.cd()
h_divide_13.Draw("e")
line.Draw("same")
c2b.Print("plots/ext_timing_divide_13.pdf")

c2c = TCanvas()
c2c.cd()
h_divide_14.Draw("e")
line.Draw("same")
c2c.Print("plots/ext_timing_divide_14.pdf")

c2di = TCanvas()
c2di.cd()
h_divide_23.Draw("e")
line.Draw("same")
c2di.Print("plots/ext_timing_divide_23.pdf")

c2dii = TCanvas()
c2dii.cd()
h_divide_24.Draw("e")
line.Draw("same")
c2dii.Print("plots/ext_timing_divide_24.pdf")

c2diii = TCanvas()
c2diii.cd()
h_divide_34.Draw("e")
line.Draw("same")
c2diii.Print("plots/ext_timing_divide_34.pdf")

c1e = TCanvas()
c1e.cd()
h_timing_1.Scale(intime_scale_factor_1)
h_timing_2.Scale(intime_scale_factor_2)
h_timing_3.Scale(intime_scale_factor_3)
h_timing_4.Scale(intime_scale_factor_4)
h_timing_1.SetLineColor(1)
h_timing_1.Draw("e hist")
h_timing_2.SetLineColor(46)
h_timing_2.Draw("e hist same")
h_timing_3.SetLineColor(30)
h_timing_3.Draw("e hist same")
h_timing_4.SetLineColor(9)
h_timing_4.Draw("e hist same")
c1e.Print("plots/ext_timing_overlay.pdf")

c2d = TCanvas()
c2d.cd()
h_pe_1.Draw()
c2d.Print("plots/ext_pe_1.pdf")

c2e = TCanvas()
c2e.cd()
h_pe_2.Draw()
c2e.Print("plots/ext_pe_2.pdf")

c2f = TCanvas()
c2f.cd()
h_pe_3.Draw()
c2f.Print("plots/ext_pe_3.pdf")

c2g = TCanvas()
c2g.cd()
h_pe_4.Draw()
c2g.Print("plots/ext_pe_4.pdf")

c2h = TCanvas()
c2h.cd()
h_pe_1.SetLineColor(1)
h_pe_1.Scale(intime_scale_factor_1)
h_pe_1.Draw()
h_pe_2.SetLineColor(46)
h_pe_2.Scale(intime_scale_factor_2)
h_pe_2.Draw("same")
h_pe_3.SetLineColor(30)
h_pe_3.Scale(intime_scale_factor_3)
h_pe_3.Draw("same")
h_pe_4.SetLineColor(9)
h_pe_4.Scale(intime_scale_factor_4)
h_pe_4.Draw("same")
c2h.Print("plots/ext_pe_overlay.pdf")

c2i = TCanvas()
c2i.cd()
h_timing_pe_1.GetXaxis().SetTitle("Flash Time [#mus]")
h_timing_pe_1.GetYaxis().SetTitle("Flash PE")
h_timing_pe_1.Draw("colz")
c2i.Print("plots/ext_timing_pe_1.pdf")

c2j = TCanvas()
c2j.cd()
h_timing_pe_2.GetXaxis().SetTitle("Flash Time [#mus]")
h_timing_pe_2.GetYaxis().SetTitle("Flash PE")
h_timing_pe_2.Draw("colz")
c2j.Print("plots/ext_timing_pe_2.pdf")

c2k = TCanvas()
c2k.cd()
h_timing_pe_3.GetXaxis().SetTitle("Flash Time [#mus]")
h_timing_pe_3.GetYaxis().SetTitle("Flash PE")
h_timing_pe_3.Draw("colz")
c2k.Print("plots/ext_timing_pe_3.pdf")

c2l = TCanvas()
c2l.cd()
h_timing_pe_4.GetXaxis().SetTitle("Flash Time [#mus]")
h_timing_pe_4.GetYaxis().SetTitle("Flash PE")
h_timing_pe_4.Draw("colz")
c2l.Print("plots/ext_timing_pe_4.pdf")

# c2m = TCanvas()
# c2m.cd()

# now let's bring in the On-Beam so we can calculate the average factor
# we should be scaling-up by

file5_path = sys.argv[5]
file5 = TFile(file5_path, 'READ')
if(file5.IsOpen()):
    print 'File ', file5_path, ' is open'
if(file5.IsOpen() == False):
    quit()
dir5 = file5.AnalyzeTPCO
opt_tree5 = dir5.Get("optical_tree")
num_events5 = opt_tree5.GetEntries()
print "Number of Events (On-Beam): ", num_events5

h_timing_5 = TH1D("h_timing_5", "h_timing_5", 60, 0, 20)
h_pe_5 = TH1D("h_pe_5", "h_pe_5", 80, 0, 2000)
h_timing_pe_5 = TH2D("h_timing_pe_5", "h_timing_pe_5", 80, 0, 20, 80, 0, 2000)

h_timing_ext_total = TH1D("h_timing_ext_total",
                          "h_timing_ext_total", 60, 0, 20)
h_timing_ext_234 = TH1D("h_timing_ext_234", "h_timing_ext_234", 60, 0, 20)

h_pe_ext_total = TH1D("h_pe_ext_total", "h_pe_ext_total", 80, 0, 2000)
h_timing_pe_ext_total = TH2D("h_timing_pe_ext_total",
                             "h_timing_pe_ext_total", 80, 0, 20, 80, 0, 2000)

h_timing_ext_234.Add(h_timing_2_2, 1)
h_timing_ext_234.Add(h_timing_3_2, 1)
h_timing_ext_234.Add(h_timing_4_2, 1)
h_timing_ext_234.Scale(intime_scale_factor_234)

h_timing_ext_total.Add(h_timing_1_2, 1)
h_timing_ext_total.Add(h_timing_2_2, 1)
h_timing_ext_total.Add(h_timing_3_2, 1)
h_timing_ext_total.Add(h_timing_4_2, 1)
h_timing_ext_total.Scale(intime_scale_factor_total)

h_pe_ext_total.Add(h_pe_1_2, 1)
h_pe_ext_total.Add(h_pe_2_2, 1)
h_pe_ext_total.Add(h_pe_3_2, 1)
h_pe_ext_total.Add(h_pe_4_2, 1)
h_pe_ext_total.Scale(intime_scale_factor_total)

h_timing_pe_ext_total.Add(h_timing_pe_1_2, 1)
h_timing_pe_ext_total.Add(h_timing_pe_2_2, 1)
h_timing_pe_ext_total.Add(h_timing_pe_3_2, 1)
h_timing_pe_ext_total.Add(h_timing_pe_4_2, 1)
h_timing_pe_ext_total.Scale(intime_scale_factor_total)

for num in range(num_events5):
    opt_tree5.GetEntry(num)
    flash_time = opt_tree5.OpFlashTime
    flash_pe = opt_tree5.OpFlashPE
    h_timing_5.Fill(flash_time)
    h_pe_5.Fill(flash_pe)
    h_timing_pe_5.Fill(flash_time, flash_pe)


h_on_minus_off_1 = h_timing_5.Clone()
h_on_minus_off_2 = h_timing_5.Clone()
h_on_minus_off_3 = h_timing_5.Clone()
h_on_minus_off_4 = h_timing_5.Clone()
h_on_minus_off_total = h_timing_5.Clone()

h_on_minus_off_1.Add(h_timing_1, -1)
h_on_minus_off_2.Add(h_timing_2, -1)
h_on_minus_off_3.Add(h_timing_3, -1)
h_on_minus_off_4.Add(h_timing_4, -1)
h_on_minus_off_total.Add(h_timing_ext_total, -1)


h_divide_on_off = h_timing_5.Clone()
h_divide_on_off.Divide(h_timing_ext_total)

h_divide_on_off_234 = h_timing_5.Clone()
h_divide_on_off_234.Divide(h_timing_ext_234)

h_divide_on_off_1 = h_timing_5.Clone()
h_divide_on_off_1.Divide(h_timing_1)

h_divide_on_off_2 = h_timing_5.Clone()
h_divide_on_off_2.Divide(h_timing_2)

h_divide_on_off_3 = h_timing_5.Clone()
h_divide_on_off_3.Divide(h_timing_3)

h_divide_on_off_4 = h_timing_5.Clone()
h_divide_on_off_4.Divide(h_timing_4)

h_divide_timing_pe_on_off = h_timing_pe_5.Clone()
h_divide_timing_pe_on_off.Divide(h_timing_pe_ext_total)

# now let's see if we can get some normalisation factors
bin_val_v_1 = []
bin_val_v_2 = []
for bin_num in range(h_divide_on_off.GetNbinsX() - 1):
    if(bin_num == 0):
        continue
    bin_val = h_divide_on_off.GetBinContent(bin_num)
    bin_center = h_divide_on_off.GetBinCenter(bin_num)
    if(bin_center < 6):
        bin_val_v_1.append(bin_val)
    if(bin_center > 15):
        bin_val_v_2.append(bin_val)

median_pre_beam = median(bin_val_v_1)
median_post_beam = median(bin_val_v_2)

print "Median scaling before beam window: ", median_pre_beam
print "Median scaling after  beam window: ", median_post_beam

c3a = TCanvas()
c3a.cd()
h_timing_5.Draw("e")
c3a.Print("plots/on_beam_timing.pdf")

c3b = TCanvas()
c3b.cd()
h_timing_5.Draw("e")
h_timing_ext_total.SetLineColor(46)
h_timing_ext_total.Draw("e same")
c3b.Print("plots/on_beam_ext_overlay.pdf")

c3c = TCanvas()
c3c.cd()
h_divide_on_off.Draw("e")
line.Draw("same")
c3c.Print("plots/on_beam_ext_divide.pdf")

c3d = TCanvas()
c3d.cd()
h_pe_5.Draw()
c3d.Print("plots/on_beam_pe.pdf")

c3e = TCanvas()
c3e.cd()
h_pe_5.Draw()
h_pe_ext_total.SetLineColor(46)
h_pe_ext_total.Draw("same")
c3e.Print("plots/on_beam_ext_pe_overlay.pdf")

c3f = TCanvas()
c3f.cd()
h_timing_pe_5.GetXaxis().SetTitle("Flash Time [#mus]")
h_timing_pe_5.GetYaxis().SetTitle("Flash PE")
h_timing_pe_5.Draw("colz")
c3f.Print("plots/on_beam_timing_pe.pdf")

c3g = TCanvas()
c3g.cd()
h_divide_timing_pe_on_off.GetXaxis().SetTitle("Flash Time [#mus]")
h_divide_timing_pe_on_off.GetYaxis().SetTitle("Flash PE")
c3g.SetLogz()
h_divide_timing_pe_on_off.Draw("colz")
c3g.Print("plots/on_beam_ext_timing_pe_divide.pdf")

c3h = TCanvas()
c3h.cd()
h_on_minus_off_1.Draw()
line.Draw("same")
c3h.Print("plots/on_minus_off_1.pdf")

c3i = TCanvas()
c3i.cd()
h_on_minus_off_2.Draw()
line.Draw("same")
c3i.Print("plots/on_minus_off_2.pdf")

c3j = TCanvas()
c3j.cd()
h_on_minus_off_3.Draw()
line.Draw("same")
c3j.Print("plots/on_minus_off_3.pdf")

c3k = TCanvas()
c3k.cd()
h_on_minus_off_4.Draw()
line.Draw("same")
c3k.Print("plots/on_minus_off_4.pdf")

c3l = TCanvas()
c3l.cd()
h_on_minus_off_total.Draw()
line.Draw("same")
c3l.Print("plots/on_minus_off_total.pdf")

c3m = TCanvas()
c3m.cd()
h_on_minus_off_1.SetLineColor(1)
h_on_minus_off_1.Draw()
h_on_minus_off_2.SetLineColor(46)
h_on_minus_off_2.Draw("same")
h_on_minus_off_3.SetLineColor(30)
h_on_minus_off_3.Draw("same")
h_on_minus_off_4.SetLineColor(9)
h_on_minus_off_4.Draw("same")
line.Draw("same")
c3m.Print("plots/on_minus_off_overlay.pdf")

c3n = TCanvas()
c3n.cd()
h_divide_on_off_234.Draw("e")
line.Draw("same")
c3n.Print("plots/on_off_divide_234.pdf")

c3o = TCanvas()
c3o.cd()
h_divide_on_off_234.SetLineColor(1)
h_divide_on_off_234.Draw("e")
h_divide_on_off.SetLineColor(30)
h_divide_on_off.Draw("e same")
line.Draw("same")
c3o.Print("plots/on_off_divide_overlay.pdf")

c3p = TCanvas()
c3p.cd()
h_divide_on_off_1.SetLineColor(1)
h_divide_on_off_1.Draw("e")
h_divide_on_off_2.SetLineColor(2)
h_divide_on_off_2.Draw("e same")
h_divide_on_off_3.SetLineColor(30)
h_divide_on_off_3.Draw("e same")
h_divide_on_off_4.SetLineColor(9)
h_divide_on_off_4.Draw("e same")
line.Draw("same")
leg_divide = TLegend()
leg_divide.AddEntry(h_divide_on_off_1, "Sample 1", "l")
leg_divide.AddEntry(h_divide_on_off_2, "Sample 2", "l")
leg_divide.AddEntry(h_divide_on_off_3, "Sample 3", "l")
leg_divide.AddEntry(h_divide_on_off_4, "Sample 4", "l")
leg_divide.Draw()
c3p.Print("plots/on_off_divide_full_overlay.pdf")
