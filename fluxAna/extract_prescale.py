import sys
import ROOT
from ROOT import TFile, TTree, TH1D, TCanvas, gROOT, TPad, TGaxis, TColor, TLegend
#import matplotlib

#file1_path = "/uboone/data/users/chill2/prescalecheck_1.root"
file1_path = sys.argv[1]
file1 = TFile(file1_path, 'READ')
if(file1.IsOpen()):
    print 'File ', file1_path, ' is open'
if(file1.IsOpen() == False):
    quit()
dir1 = file1.PrescaleCheck
prescale_tree = dir1.Get("prescale_tree")
num_events = prescale_tree.GetEntries()
print "Number of Events: ", num_events

h_prescale = TH1D("h_prescale", "h_prescale", 2700, 4900, 7600)
h_prescale_val = TH1D("h_prescale_val", "h_prescale_val", 200, 0, 110)
h_unique_run = TH1D("h_unique_run", "h_unique_run", 2700, 4900, 7600)

h_prescale_val_bnb = TH1D("h_prescale_val_bnb", "h_prescale_val_bnb", 200, 0, 110)
h_prescale_val_unbiased = TH1D("h_prescale_val_unbiased", "h_prescale_val_unbiased", 200, 0, 110)

h_pass_numi = TH1D("h_pass_numi", "h_pass_numi", 2, 0, 2)
h_pass_bnb = TH1D("h_pass_bnb", "h_pass_bnb", 2, 0, 2)
h_pass_unbiased = TH1D("h_pass_unbiased", "h_pass_unbiased", 2, 0, 2)
h_pass_combo = TH1D("h_pass_combo", "h_pass_combo", 9, 0, 9)

prev_run_number = 0
num_unique_runs = 0

prescale_factor_v = []
run_number_v = []

pair_v = []

for event in range(num_events):
  prescale_tree.GetEntry(event)
  prescale_factor = prescale_tree.PrescaleNuMI
  prescale_factor_bnb = prescale_tree.PrescaleBNB
  prescale_factor_unbiased = prescale_tree.PrescaleUnbiased
  run_number = prescale_tree.Run

  pass_numi = prescale_tree.PassedNuMI
  pass_bnb = prescale_tree.PassedBNB
  pass_unbiased = prescale_tree.PassedUnbiased

  prescale_factor_v.append(prescale_factor)
  run_number_v.append(run_number)
  pair_v.append((run_number, prescale_factor, prescale_factor_bnb, prescale_factor_unbiased, pass_numi, pass_bnb, pass_unbiased))

pair_v.sort(key=lambda tup: tup[0])

unique_prescale = []

for pair in pair_v:
  run_number = pair[0]
  prescale_factor = pair[1]
  prescale_factor_bnb = pair[2]
  prescale_factor_unbiased = pair[3]

  pass_numi = pair[4]
  pass_bnb = pair[5]
  pass_unbiased = pair[6]

  #print run_number, prescale_factor

  if(run_number != prev_run_number):
    print run_number
    h_prescale.Fill(run_number, prescale_factor * 100.)

    h_prescale_val.Fill(prescale_factor * 100.) 
    h_prescale_val_bnb.Fill(prescale_factor_bnb * 100.)
    h_prescale_val_unbiased.Fill(prescale_factor_unbiased * 100.)

    h_pass_numi.Fill(pass_numi)
    h_pass_bnb.Fill(pass_bnb)
    h_pass_unbiased.Fill(pass_unbiased)

    if(not pass_numi and not pass_bnb and not pass_unbiased):
      h_pass_combo.Fill(0)
    if(not pass_numi and pass_bnb and pass_unbiased):
      h_pass_combo.Fill(1)
    if(not pass_numi and not pass_bnb and pass_unbiased):
      h_pass_combo.Fill(2)
    if(not pass_numi and pass_bnb and not pass_unbiased):
      h_pass_combo.Fill(3)
    if(pass_numi):
      if(pass_bnb):
        if(pass_unbiased):
          h_pass_combo.Fill(4)
        if(not pass_unbiased):
          h_pass_combo.Fill(5)
      if(not pass_bnb):
        if(pass_unbiased):
          h_pass_combo.Fill(6)
        if(not pass_unbiased):
          h_pass_combo.Fill(7)

    h_unique_run.Fill(run_number)
    num_unique_runs = num_unique_runs + 1

    unique_prescale.append(prescale_factor * 100.)

  prev_run_number = run_number  


print "Average Prescale Factor: ", sum(unique_prescale)/float(len(unique_prescale))

print "Number of Unique Runs: ", num_unique_runs

c1 = TCanvas()
c1.cd()
h_prescale.GetXaxis().SetTitle("Run Number")
h_prescale.GetYaxis().SetTitle("Prescale Factor")
h_prescale.Draw("hist")
c1.Print("prescale_plots/prescale_factor.pdf")

c2 = TCanvas()
c2.cd()
h_prescale_val.GetXaxis().SetTitle("Prescale Value [%]")
h_prescale_val.Draw("hist")
c2.Print("prescale_plots/prescale_value.pdf")

c3 = TCanvas()
c3.cd()
h_unique_run.GetXaxis().SetTitle("Run Number")
h_unique_run.Draw()
c3.Print("prescale_plots/unique_run_number.pdf")

c4 = TCanvas()
c4.cd()
h_prescale_val_bnb.GetXaxis().SetTitle("Prescale Value [%]")
h_prescale_val_bnb.Draw()
c4.Print("prescale_plots/prescale_value_bnb.pdf")

c5 = TCanvas()
c5.cd()
h_prescale_val_unbiased.GetXaxis().SetTitle("Prescale Value [%]")
h_prescale_val_unbiased.Draw()
c5.Print("prescale_plots/prescale_value_unbiased.pdf")

c6a = TCanvas()
c6a.cd()
h_pass_numi.Draw()
c6a.Print("prescale_plots/pass_numi.pdf")

c6b = TCanvas()
c6b.cd()
h_pass_bnb.Draw()
c6b.Print("prescale_plots/pass_bnb.pdf")

c6c = TCanvas()
c6c.cd()
h_pass_unbiased.Draw()
c6c.Print("prescale_plots/pass_unbiased.pdf")

c6d = TCanvas()
c6d.cd()
h_pass_combo.Draw()
c6d.Print("prescale_plots/pass_combo.pdf")

wait = raw_input("PRESS ENTER TO CONTINUE.")
