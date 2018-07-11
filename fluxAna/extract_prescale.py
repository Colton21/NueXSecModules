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

prev_run_number = 0
num_unique_runs = 0

prescale_factor_v = []
run_number_v = []

pair_v = []

for event in range(num_events):
  prescale_tree.GetEntry(event)
  prescale_factor = prescale_tree.PrescaleNuMI
  run_number = prescale_tree.Run

  prescale_factor_v.append(prescale_factor)
  run_number_v.append(run_number)
  pair_v.append((run_number, prescale_factor))

pair_v.sort(key=lambda tup: tup[0])

unique_prescale = []

for pair in pair_v:
  run_number = pair[0]
  prescale_factor = pair[1]

  #print run_number, prescale_factor

  if(run_number != prev_run_number):
    print run_number
    h_prescale.Fill(run_number, prescale_factor * 100.)
    h_prescale_val.Fill(prescale_factor * 100.) 
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
c1.Print("prescale_factor.pdf")

c2 = TCanvas()
c2.cd()
h_prescale_val.GetXaxis().SetTitle("Prescale Value [%]")
h_prescale_val.Draw("hist")
c2.Print("prescale_value.pdf")

c3 = TCanvas()
c3.cd()
h_unique_run.GetXaxis().SetTitle("Run Number")
h_unique_run.Draw()
c3.Print("unique_run_number.pdf")

wait = raw_input("PRESS ENTER TO CONTINUE.")
