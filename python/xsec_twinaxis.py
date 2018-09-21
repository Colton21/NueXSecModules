import sys
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TCanvas, gROOT, TPad, TGaxis, TColor, TLegend, TLine
import matplotlib.pyplot as plt
import numpy as np

file1_path = sys.argv[1]
file1 = TFile(file1_path, 'READ')
if(file1.IsOpen()):
    print 'File ', file1_path, ' is open'
if(file1.IsOpen() == False):
    quit()
nue_flux = file1.Get("nueFluxHisto")
anue_flux = file1.Get("anueFluxHisto")

nue_flux_list = []
anue_flux_list = []
bin_flux_list = []
a_bin_flux_list = []

nue_flux_bins = nue_flux.GetNbinsX()
anue_flux_bins = anue_flux.GetNbinsX()
average_num = 0
#average_den = (nue_flux.GetEntries() + anue_flux.GetEntries())
average_den = (nue_flux.Integral() + anue_flux.Integral())

for bin in range(nue_flux_bins):
    nue_flux_list.append(nue_flux.GetBinContent(bin))
    bin_flux_list.append(nue_flux.GetBinCenter(bin))
    average_num = average_num + \
        (nue_flux.GetBinContent(bin) * nue_flux.GetBinCenter(bin))

for bin in range(anue_flux_bins):
    anue_flux_list.append(anue_flux.GetBinContent(bin))
    a_bin_flux_list.append(anue_flux.GetBinCenter(bin))
    average_num = average_num + \
        (anue_flux.GetBinContent(bin) * anue_flux.GetBinCenter(bin))

average = float(average_num) / float(average_den)
print 'Average: ', average

file2_path = sys.argv[2]
file2 = TFile(file2_path, 'READ')
if(file2.IsOpen()):
    print 'File ', file2_path, ' is open'
if(file2.IsOpen() == False):
    quit()
dir2 = file2.nu_e_Ar40
nue_xsec = dir2.Get("tot_cc")
dir3 = file2.nu_e_bar_Ar40
anue_xsec = dir3.Get("tot_cc")

# Create buffers
x_buff_nue = nue_xsec.GetX()
y_buff_nue = nue_xsec.GetY()
N_nue = nue_xsec.GetN()
x_buff_nue.SetSize(N_nue)
y_buff_nue.SetSize(N_nue)
# Create arrays from buffers, copy to prevent data loss
x_arr_nue = np.array(x_buff_nue, copy=True)
y_arr_nue = np.array(y_buff_nue, copy=True)

# Create buffers
x_buff_anue = anue_xsec.GetX()
y_buff_anue = anue_xsec.GetY()
N_anue = anue_xsec.GetN()
x_buff_anue.SetSize(N_anue)
y_buff_anue.SetSize(N_anue)
# Create arrays from buffers, copy to prevent data loss
x_arr_anue = np.array(x_buff_anue, copy=True)
y_arr_anue = np.array(y_buff_anue, copy=True)

y_arr_nue = np.multiply(y_arr_nue, (1.0e-38))
y_arr_anue = np.multiply(y_arr_anue, (1.0e-38))

y_arr_nue = np.multiply(y_arr_nue, (1. / 40.))
y_arr_anue = np.multiply(y_arr_anue, (1. / 40.))

fig, ax1 = plt.subplots()
line1a = ax1.plot(bin_flux_list, nue_flux_list, 'darkcyan',
                  linewidth=2.0, label=r'NuMI $\nu_{e}$ Flux')
ax1.fill_between(bin_flux_list, nue_flux_list, 0,
                 facecolor='darkgray', alpha=0.2)
line1b = ax1.plot(a_bin_flux_list, anue_flux_list, 'royalblue',
                  linewidth=2.0, label=r'NuMI $\bar{\nu}_{e}$ Flux')
ax1.fill_between(a_bin_flux_list, anue_flux_list,
                 0, facecolor='darkgray', alpha=0.2)
ax1.set_xlabel('Neutrino Energy [GeV]')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r'$\nu_{e} | \bar{\nu_{e}}$ / cm$^{2}$ / 6e20 POT')
# ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
line2a = ax2.plot(x_arr_nue, y_arr_nue, 'darkcyan',
                  linewidth=2.0, label=r'GENIE $\nu_{e}$ Spline', linestyle=':')
line2b = ax2.plot(x_arr_anue, y_arr_anue, 'royalblue',
                  linewidth=2.0, label=r'GENIE $\bar{\nu}_{e}$ Spline', linestyle=':')
ax2.set_ylabel(r'$\nu_{e}^{(-)}$ CC Cross Section [cm$^{2}$]')
# ax2.tick_params('y', colors='r')

ax3 = ax2
genie_xsec_point = 4.83114e-39
# for now we're taking the Stat from data
genie_xsec_point_stat_err = 0.703e-39
genie_xsec_point_sys_err = 0.25 * genie_xsec_point
genie_xsec_point_total_err = np.sqrt(
    genie_xsec_point_stat_err**2 + genie_xsec_point_sys_err**2)
genie_xsec_energy = average

genie_xsec_point_x_err = genie_xsec_energy * 0.2

line3a = ax3.plot(genie_xsec_energy, genie_xsec_point,
                  'salmon', label=r'MC Xsec')
line3b = ax3.errorbar(genie_xsec_energy, genie_xsec_point, yerr=[genie_xsec_point_total_err], xerr=[
    genie_xsec_point_x_err], fmt='--o', color='salmon', ecolor='salmon', linewidth=1.0, markersize=6.0, label='xsec_point')

lines = line1a + line1b + line3a + line2a + line2b
labels = [l.get_label() for l in lines]
# ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)
ax3.legend(lines, labels, loc=9, ncol=2)

ax2.set_ylim(0, 35e-39)
plt.xlim(0.0, 3.0)

plt.show()
