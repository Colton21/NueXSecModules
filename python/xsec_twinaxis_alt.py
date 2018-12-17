import sys
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TCanvas, gROOT, TPad, TGaxis, TColor, TLegend, TLine, TGraph
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
# average_den = (nue_flux.GetEntries() + anue_flux.GetEntries())
average_den = (nue_flux.Integral() + anue_flux.Integral())

summed_flux = TH1D("summed_flux", "summed_flux", 400, 0, 20)
summed_flux.Add(nue_flux, 1)
summed_flux.Add(anue_flux, 1)

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

average_bin = 0
for bin in range(nue_flux_bins):
    if(summed_flux.GetBinCenter(bin) <= average):
        continue
    # here we should get the bin with the average in it
    average_bin = bin
    break

max_sum = 0
max_bin_val = 0
summed_flux_integral = summed_flux.Integral()
print summed_flux_integral
for bin in range(average_bin, nue_flux_bins):
    max_sum = max_sum + summed_flux.GetBinContent(bin)
    if(max_sum / summed_flux_integral >= 0.34):
        max_bin_val = summed_flux.GetBinCenter(bin)
        break

min_sum = 0
min_bin_val = 0
for bin in range(average_bin, 0, -1):
    print "Bin: ", bin
    min_sum = min_sum + summed_flux.GetBinContent(bin)
    if(min_sum / summed_flux_integral >= 0.34):
        min_bin_val = summed_flux.GetBinCenter(bin)
        break

print "MaxBinVal: ", max_bin_val
print "MinBinVal: ", min_bin_val

#
# res = 100
# elow = nue_flux.GetXaxis().GetBinLowEdge(1)
# ehigh = nue_flux.GetXaxis().GetBinLowEdge(nue_flux_bins + 1)

# fine = TH1D("fine", "fine", nue_flux_bins * res, elow, ehigh)
# temp = TGraph()

# i = 0
# for bin in range(nue_flux_bins):
#    E = nue_flux.GetXaxis().GetBinCenter(i + 1)
#    C = nue_flux.GetBinContent(i + 1)
#    W = nue_flux.GetXaxis().GetBinWidth(i + 1)

#    if (W != 0.0):
#        temp.SetPoint(temp.GetN(), E, C / W)
#    i = i + 1

# j = 0
# for bin in range(fine.GetNbinsX()):
#    E = fine.GetXaxis().GetBinCenter(j + 1)
#    W = fine.GetBinWidth(j + 1)

#    fine.SetBinContent(j + 1, temp.Eval(E, 0, "S") * W)
#    j = j + 1

# fine.Scale(nue_flux.Integral(1, nue_flux.GetNbinsX() + 1) /
#           fine.Integral(1, fine.GetNbinsX() + 1))

# print "Interpolation Difference = ", fine.Integral(1, fine.GetNbinsX() +
# 1), "/", nue_flux.Integral(1, nue_flux.GetNbinsX() + 1)

# fine_flux_list = []
# fine_flux_bins = []
# for bin in range(fine.GetNbinsX()):
#    fine_flux_list.append(nue_flux.GetBinContent(bin))
#    fine_flux_bins.append(nue_flux.GetBinCenter(bin))

##############################################################################

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

##############################################################################

fig, ax1 = plt.subplots()
# line0 = ax1.plot(fine_flux_bins, fine_flux_list, 'black', linewidth=2.0)
line1a = ax1.plot(bin_flux_list, nue_flux_list, 'mediumblue',
                  linewidth=2.0, label=r'NuMI $\nu_{e}$ Flux')
ax1.fill_between(bin_flux_list, nue_flux_list, 0,
                 facecolor='darkgray', alpha=0.2)
line1b = ax1.plot(a_bin_flux_list, anue_flux_list, 'darkgreen',
                  linewidth=2.0, label=r'NuMI $\bar{\nu}_{e}$ Flux')
ax1.fill_between(a_bin_flux_list, anue_flux_list,
                 0, facecolor='darkgray', alpha=0.2)
ax1.set_xlabel('Neutrino Energy [GeV]')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r'$\nu_{e}/\bar{\nu_{e}}$  / cm$^{2}$ / 6e20 POT')
# ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
line2a = ax2.plot(x_arr_nue, y_arr_nue, 'mediumblue',
                  linewidth=2.4, label=r'GENIE $\nu_{e}$ Cross Section', linestyle=':')
line2b = ax2.plot(x_arr_anue, y_arr_anue, 'darkgreen',
                  linewidth=2.4, label=r'GENIE $\bar{\nu}_{e}$ Cross Section', linestyle=':')
ax2.set_ylabel(r'$\nu_{e}/\bar{\nu}_{e}$ CC Cross Section [cm$^{2}$]')
# ax2.tick_params('y', colors='r')

ax3 = ax2
genie_xsec_point = 4.83114e-39
# for now we're taking the Stat from data
genie_xsec_point_stat_err = 0.703e-39
genie_xsec_point_sys_err = 0.30 * genie_xsec_point
genie_xsec_point_total_err = np.sqrt(
    genie_xsec_point_stat_err**2 + genie_xsec_point_sys_err**2)
genie_xsec_energy = average
asymmetric_x_err = np.array([[average - min_bin_val, max_bin_val - average]]).T
symmetric_y_err = np.array(
    [[genie_xsec_point_total_err, genie_xsec_point_total_err]]).T

genie_xsec_point_x_err = genie_xsec_energy * 0.2

# line3a = ax3.plot(genie_xsec_energy, genie_xsec_point,
#                  'black', label = r'MC Xsec')
# line3b = ax3.errorbar(genie_xsec_energy, genie_xsec_point,
# xerr=asymmetric_x_err, yerr=symmetric_y_err, fmt='--o',
# color = 'black', ecolor = 'black', linewidth = 1.0, markersize = 6.0,
# label = 'xsec_point')

lines = line1a + line1b + line2a + line2b
labels = [l.get_label() for l in lines]
# ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)
ax3.legend(lines, labels, loc=9, ncol=2)

ax2.set_ylim(0, 35e-39)
plt.xlim(0.0, 3.0)

plt.show()

# figure 2
fig2, ax_2 = plt.subplots()

# these are to make the legend have a fill, rather than line
dummy_hist_1a = ax_2.hist(
    (-1, -2), 1, fill=True, label=r'Expected Total Bound', color='darkgray', alpha=0.2)
# ax.legend(markerscale=2,fontsize=20)
dummy_hist_1b = ax_2.hist(
    (-1, -2), 1, fill=True, label=r'Expected Stats Bound', color='maroon', alpha=0.2)

genie_xsec_point_nue = (6.34569e-39, 6.34569e-39)
genie_xsec_point_nue_bar = (2.24685e-39, 2.24685e-39)
genie_xsec_point_nue_nue_bar = (5.15e-39, 5.15e-39)

x_edge = (-0.9, 1.1)

# point_1a = ax_2.plot(0.5, genie_xsec_point, 'black',
#                     label=r'MC $\nu_{e}$ + $\bar{\nu}_{e}$')
# point_1b = ax_2.errorbar(0.5, genie_xsec_point, yerr=symmetric_y_err,
# fmt='--o', color='black', ecolor='black', linewidth=1.0, markersize=6.0)

lower_1a = genie_xsec_point + (genie_xsec_point * 0.14)
upper_1a = genie_xsec_point - (genie_xsec_point * 0.14)

lower_1a_list = (lower_1a, lower_1a)
upper_1a_list = (upper_1a, upper_1a)

quadrature = np.sqrt(0.14**2 + 0.30**2)

lower_1b = genie_xsec_point + (genie_xsec_point * quadrature)
upper_1b = genie_xsec_point - (genie_xsec_point * quadrature)

lower_1b_list = (lower_1b, lower_1b)
upper_1b_list = (upper_1b, upper_1b)

band_1a_lower = ax_2.plot(x_edge, lower_1a_list,
                          color='maroon')
band_1a_upper = ax_2.plot(x_edge, upper_1a_list, color='maroon')

band_1b_lower = ax_2.plot(x_edge, lower_1b_list,
                          color='black')
band_1b_upper = ax_2.plot(x_edge, upper_1b_list, color='black')

ax_2.fill_between(x_edge, upper_1a_list, lower_1a_list, hatch='//',
                  facecolor='maroon', alpha=0.2)

ax_2.fill_between(x_edge, upper_1b_list, lower_1b_list,
                  facecolor='darkgray', alpha=0.2)

genie_xsec_point_list = (genie_xsec_point, genie_xsec_point)

line_2a = ax_2.plot(x_edge, genie_xsec_point_nue,
                    color='mediumblue', linewidth=1.5, label=r'Genie $\nu_{e}$ Cross Section', linestyle='--')
line_2b = ax_2.plot(x_edge, genie_xsec_point_nue_bar,
                    color='darkgreen', linewidth=1.5, label=r'Genie $\bar{\nu}_{e}$ Cross Section', linestyle='--')
line_2c = ax_2.plot(x_edge, genie_xsec_point_list,
                    color='#AA66CC', linewidth=1.5, label=r'Genie $\nu_{e} + \bar{\nu}_{e}$ Cross Section', linestyle='--')


ax_2.set_ylim(0, 10e-39)
ax_2.set_xlim(0, 1)
ax_2.set_ylabel(r'$\nu_{e}/\bar{\nu}_{e}$ CC Cross Section [cm$^{2}$]')
ax_2.xaxis.label.set_size(0)
ax_2.xaxis.set_ticklabels([])
ax_2.xaxis.set_visible(False)

#lines_2 = band_1a_lower + band_1b_lower + line_2a + line_2b + line_2c
lines_2 = (dummy_hist_1a, dummy_hist_1b, line_2a, line_2b, line_2c)
#labels_2 = [l.get_label() for l in lines_2]
# ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)
labels_list = (r'Expected Total Bound', r'Expected Stats Bound',
               r'Genie $\nu_{e}$ Cross Section', r'Genie $\bar{\nu}_{e}$ Cross Section', r'Genie $\nu_{e} + \bar{\nu}_{e}$ Cross Section')
ax_2.legend(loc=9, ncol=2)

plt.show()

##############################################################################

fig3, ax3 = plt.subplots()
# line0 = ax1.plot(fine_flux_bins, fine_flux_list, 'black', linewidth=2.0)
line3a = ax3.plot(bin_flux_list, nue_flux_list, 'mediumblue',
                  linewidth=2.0, label=r'NuMI $\nu_{e}$ Flux')
ax3.fill_between(bin_flux_list, nue_flux_list, 0,
                 facecolor='darkgray', alpha=0.2)
line3b = ax3.plot(a_bin_flux_list, anue_flux_list, 'darkgreen',
                  linewidth=2.0, label=r'NuMI $\bar{\nu}_{e}$ Flux')
ax3.fill_between(a_bin_flux_list, anue_flux_list,
                 0, facecolor='darkgray', alpha=0.2)
ax3.set_xlabel('Neutrino Energy [GeV]')
# make the y-axis label, ticks and tick labels match the line color.
ax3.set_ylabel(r'$\nu_{e}/\bar{\nu_{e}}$  / cm$^{2}$ / 6e20 POT')
# ax1.tick_params('y', colors='b')

ax4 = ax3.twinx()
line4a = ax4.plot(x_arr_nue, y_arr_nue, 'mediumblue',
                  linewidth=2.4, label=r'GENIE $\nu_{e}$ Cross Section', linestyle=':')
line4b = ax4.plot(x_arr_anue, y_arr_anue, 'darkgreen',
                  linewidth=2.4, label=r'GENIE $\bar{\nu}_{e}$ Cross Section', linestyle=':')
ax4.set_ylabel(r'$\nu_{e}/\bar{\nu}_{e}$ CC Cross Section [cm$^{2}$]')
# ax2.tick_params('y', colors='r')

ax5 = ax4
#genie_xsec_point = 4.83114e-39
# for now we're taking the Stat from data
#genie_xsec_point_stat_err = 0.703e-39
#genie_xsec_point_sys_err = 0.25 * genie_xsec_point
# genie_xsec_point_total_err = np.sqrt(
#    genie_xsec_point_stat_err**2 + genie_xsec_point_sys_err**2)
#genie_xsec_energy = average
#asymmetric_x_err = np.array([[average - min_bin_val, max_bin_val - average]]).T
# symmetric_y_err = np.array(
#    [[genie_xsec_point_total_err, genie_xsec_point_total_err]]).T

#genie_xsec_point_x_err = genie_xsec_energy * 0.2

line5a = ax4.plot(genie_xsec_energy, genie_xsec_point,
                  'black', label=r'$\nu_e$ + $\bar{\nu}_{e}$ Data Cross Section', linewidth=2.0)
line5b = ax4.errorbar(genie_xsec_energy, genie_xsec_point,
                      xerr=asymmetric_x_err, yerr=symmetric_y_err, fmt='--o',
                      color='black', ecolor='black', linewidth=2.0, markersize=6.0, capsize=5.0, markeredgewidth=2,
                      label='xsec_point')

symmetric_y_err_2 = np.array(
    [[genie_xsec_point_stat_err, genie_xsec_point_stat_err]]).T
line6a = ax4.errorbar(genie_xsec_energy, genie_xsec_point,
                      xerr=asymmetric_x_err, yerr=symmetric_y_err_2, fmt='--o', color='black', ecolor='black', linewidth=2.0, markersize=6.0,
                      capsize=5.0, markeredgewidth=2, label='xsec_point_2', linestyle='--')

lines = line3a + line3b + line5a + line4a + line4b
labels = [l.get_label() for l in lines]
# ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)
ax5.legend(lines, labels, loc=9, ncol=2)

ax4.set_ylim(0, 35e-39)
plt.xlim(0.0, 3.0)

plt.show()
