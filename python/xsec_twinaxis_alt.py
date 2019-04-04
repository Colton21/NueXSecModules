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

# down sample the flux to match with the xsec - xsec is 1000 bins from 0 to 125 GeV
# flux is 400 bins from 0 to 20 GeV, need to get down to 160 bins
nue_flux_rebin = TH1D("nue_flux_rebin", "nue_flux_rebin", 160, 0, 20)
anue_flux_rebin = TH1D("anue_flux_rebin", "anue_flux_rebin", 160, 0, 20)

nue_flux_list = []
anue_flux_list = []
bin_flux_list = []
a_bin_flux_list = []

nue_flux_list_rebin = []
anue_flux_list_rebin = []
bin_flux_list_rebin = []
a_bin_flux_list_rebin = []

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
    nue_flux_rebin.Fill(nue_flux.GetBinCenter(bin),
                        nue_flux.GetBinContent(bin))
    average_num = average_num + \
        (nue_flux.GetBinContent(bin) * nue_flux.GetBinCenter(bin))

print 'Number of Bins (Nue Flux): ', nue_flux_bins, '\n'
print 'Number of Bins (Rebin Nue): ', nue_flux_rebin.GetNbinsX(), '\n'

for bin in range(nue_flux_rebin.GetNbinsX()):
    nue_flux_list_rebin.append(nue_flux_rebin.GetBinContent(bin))
    bin_flux_list_rebin.append(nue_flux_rebin.GetBinCenter(bin))
    # print nue_flux_rebin.GetBinCenter(bin), ' ',
    # nue_flux_rebin.GetBinContent(bin), '\n'

for bin in range(anue_flux_bins):
    anue_flux_list.append(anue_flux.GetBinContent(bin))
    a_bin_flux_list.append(anue_flux.GetBinCenter(bin))
    anue_flux_rebin.Fill(anue_flux.GetBinCenter(bin),
                         anue_flux.GetBinContent(bin))
    average_num = average_num + \
        (anue_flux.GetBinContent(bin) * anue_flux.GetBinCenter(bin))

for bin in range(anue_flux_rebin.GetNbinsX()):
    anue_flux_list_rebin.append(anue_flux_rebin.GetBinContent(bin))
    a_bin_flux_list_rebin.append(anue_flux_rebin.GetBinCenter(bin))

print 'Number of Bins (NueBar Flux): ', anue_flux_bins, '\n'

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
    # print "Bin: ", bin
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

# bonus figures 0 and 01 and 001
#############################################################################
# flux times xsec
# old arrays - not as good as prod_nue_flux_xsec_list
prod_nue_flux_list = []
prod_anue_flux_list = []

# these lists are for forwarding to see how low energy flux increase
# changes results
prod_nue_flux_up_list = []
prod_anue_flux_up_list = []
prod_nue_flux_xsec_up_list = []
prod_anue_flux_xsec_up_list = []
low_eng_scale = 1.0  # 100% in this case

prod_nue_flux_xsec_list = []
prod_anue_flux_xsec_list = []

plotting_prod_nue_flux_xsec_list = []
plotting_prod_anue_flux_xsec_list = []

# better way of doing it - more bins and eval
j = 0
for val in bin_flux_list:
    xsec_value = nue_xsec.Eval(val) * 1.0e-38 / 40.
    flux_value = nue_flux_list[j]
    # if(xsec_value < 0):
    if(j < 6):
        print "Xsec (nue): ", xsec_value, " at: ", val
        xsec_value = 0.0
    if(val <= 0.5):
        prod_nue_flux_xsec_up_list.append(
            (flux_value + flux_value * low_eng_scale) * xsec_value)
    if(val > 0.5):
        prod_nue_flux_xsec_up_list.append(xsec_value * flux_value)
    prod_nue_flux_xsec_list.append(xsec_value * flux_value)
    plotting_prod_nue_flux_xsec_list.append(
        xsec_value * flux_value / 0.005 * 3.50191e+31)
    j = j + 1

j = 0
for val in bin_flux_list:
    xsec_value = anue_xsec.Eval(val) * 1.0e-38 / 40.
    flux_value = anue_flux_list[j]
    if(xsec_value < 0):
        # print "Xsec (anue): ", xsec_value, " at: ", val
        xsec_value = 0.0
    if(val <= 0.5):
        prod_anue_flux_xsec_up_list.append(
            (flux_value + flux_value * low_eng_scale) * xsec_value)
    if(val > 0.5):
        prod_anue_flux_xsec_up_list.append(xsec_value * flux_value)
    prod_anue_flux_xsec_list.append(xsec_value * flux_value)
    plotting_prod_anue_flux_xsec_list.append(
        xsec_value * flux_value / 0.005 * 3.50191e+31)
    j = j + 1


print "Integral: ", sum(plotting_prod_nue_flux_xsec_list) * 0.005

fig0, ax0 = plt.subplots()

line0a = ax0.plot(bin_flux_list,
                  plotting_prod_nue_flux_xsec_list, 'goldenrod', linewidth=2.0, label=r'NuMI $\nu_{e}$ Flux $\times$ Cross Section')
line0b = ax0.plot(a_bin_flux_list,
                  plotting_prod_anue_flux_xsec_list, 'slategray', linewidth=2.0, label=r'NuMI $\bar{\nu}_{e}$ Flux $\times$ Cross Section')

ax0.fill_between(bin_flux_list, plotting_prod_nue_flux_xsec_list, 0,
                 facecolor='darkgray', alpha=0.2)
ax0.fill_between(a_bin_flux_list, plotting_prod_anue_flux_xsec_list,
                 0, facecolor='darkgray', alpha=0.2)


ax0.set_xlabel('Neutrino Energy [GeV]', fontsize=18)
ax0.set_ylabel(r'$\nu_{e}/\bar{\nu_{e}}$ / GeV / 6e20 POT', fontsize=18)

ax0b = ax0.twinx()

line00a = ax0b.plot(x_arr_nue, y_arr_nue, 'goldenrod',
                    linewidth=2.4, label=r'GENIE $\nu_{e}$ Cross Section', linestyle=':')
line00b = ax0b.plot(x_arr_anue, y_arr_anue, 'slategray',
                    linewidth=2.4, label=r'GENIE $\bar{\nu}_{e}$ Cross Section', linestyle=':')
ax0b.set_ylabel(
    r'$\nu_{e}/\bar{\nu}_{e}$ CC Cross Section [cm$^{2}$]', fontsize=18)

lines0 = line0a + line0b + line00a + line00b
labels0 = [l.get_label() for l in lines0]

ax0.legend(lines0, labels0, loc=9, ncol=2)

ax0b.set_ylim(0, 35e-39)
# ax0.set_ylim(0, 1.0e9)
plt.xlim(0.0, 4.0)
plt.show()

#############################################################################
# figure01

# manually grab values for the efficiency, this is probably easiest,
# and we're not overly concerned with the exact values, just estiimates

# this is for after all selection cuts are applied
# binning is not uniform in energy, so giving binning scheme here
energy_efficiency_bin_list = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 4.0]
energy_efficiency_list = [0.005, 0.055, 0.085, 0.11, 0.105, 0.075]

# now bring the flux * xsec into this binning scheme
nue_rate_efficiency = [0] * 6
anue_rate_efficiency = [0] * 6

nue_rate_efficiency_up = [0] * 6
anue_rate_efficiency_up = [0] * 6

for n in range(6):
    bin_max = energy_efficiency_bin_list[n + 1]
    bin_min = energy_efficiency_bin_list[n]
    bin_width = bin_max - bin_min

    av_prod_nue_flux_xsec_list = []
    av_prod_anue_flux_xsec_list = []
    av_prod_nue_flux_xsec_up_list = []
    av_prod_anue_flux_xsec_up_list = []
    for m in range(len(bin_flux_list)):
        if(bin_flux_list[m] > energy_efficiency_bin_list[n + 1]):
            break
        if(bin_flux_list[m] >= energy_efficiency_bin_list[n] and bin_flux_list[m] < energy_efficiency_bin_list[n + 1]):

            av_prod_nue_flux_xsec_list.append(prod_nue_flux_xsec_list[m])
            av_prod_anue_flux_xsec_list.append(prod_anue_flux_xsec_list[m])
            av_prod_nue_flux_xsec_up_list.append(prod_nue_flux_xsec_up_list[m])
            av_prod_anue_flux_xsec_up_list.append(
                prod_anue_flux_xsec_up_list[m])

    av_prod_nue_flux_xsec = sum(av_prod_nue_flux_xsec_list)
    av_prod_anue_flux_xsec = sum(av_prod_anue_flux_xsec_list)
    av_prod_nue_flux_xsec_up = sum(av_prod_nue_flux_xsec_up_list)
    av_prod_anue_flux_xsec_up = sum(av_prod_anue_flux_xsec_up_list)

    eff = energy_efficiency_list[n]

    nue_rate_efficiency[n] = av_prod_nue_flux_xsec * eff / bin_width
    anue_rate_efficiency[n] = av_prod_anue_flux_xsec * eff / bin_width
    nue_rate_efficiency_up[n] = av_prod_nue_flux_xsec_up * eff / bin_width
    anue_rate_efficiency_up[n] = av_prod_anue_flux_xsec_up * eff / bin_width

# if we're multiplying by efficiency, also include num targets
num_targets = 3.50191e+31
nue_rate_efficiency = np.multiply(nue_rate_efficiency, num_targets)
anue_rate_efficiency = np.multiply(anue_rate_efficiency, num_targets)

nue_rate_efficiency_up = np.multiply(nue_rate_efficiency_up, num_targets)
anue_rate_efficiency_up = np.multiply(anue_rate_efficiency_up, num_targets)

#########
fig, ax00 = plt.subplots()

energy_efficiency_bin_position = [0.125, 0.375, 0.625, 0.875, 1.25, 2.75]

# we want individual points for each value with errors, not a line
efficiency_y_err = [0.015 / energy_efficiency_list[0], 0.02 / energy_efficiency_list[1], 0.02 / energy_efficiency_list[2],
                    0.02 / energy_efficiency_list[3], 0.02 / energy_efficiency_list[4], 0.015 / energy_efficiency_list[5]]
efficiency_x_err = [0.125, 0.125, 0.125, 0.125, 0.25, 1.25]

efficiency_nue_y_err = [0] * 6

efficiency_nue_y_err[0] = efficiency_y_err[0] * nue_rate_efficiency[0]
efficiency_nue_y_err[1] = efficiency_y_err[1] * nue_rate_efficiency[1]
efficiency_nue_y_err[2] = efficiency_y_err[2] * nue_rate_efficiency[2]
efficiency_nue_y_err[3] = efficiency_y_err[3] * nue_rate_efficiency[3]
efficiency_nue_y_err[4] = efficiency_y_err[4] * nue_rate_efficiency[4]
efficiency_nue_y_err[5] = efficiency_y_err[5] * nue_rate_efficiency[5]

efficiency_anue_y_err = [0] * 6

efficiency_anue_y_err[0] = efficiency_y_err[0] * anue_rate_efficiency[0]
efficiency_anue_y_err[1] = efficiency_y_err[1] * anue_rate_efficiency[1]
efficiency_anue_y_err[2] = efficiency_y_err[2] * anue_rate_efficiency[2]
efficiency_anue_y_err[3] = efficiency_y_err[3] * anue_rate_efficiency[3]
efficiency_anue_y_err[4] = efficiency_y_err[4] * anue_rate_efficiency[4]
efficiency_anue_y_err[5] = efficiency_y_err[5] * anue_rate_efficiency[5]

line00a = ax00.errorbar(energy_efficiency_bin_position, nue_rate_efficiency,
                        xerr=efficiency_x_err, yerr=efficiency_nue_y_err, fmt='o', color='goldenrod', ecolor='goldenrod', linewidth=2.0, markersize=6.0,
                        capsize=5.0, markeredgewidth=2, label='xsec_point_2', markeredgecolor='goldenrod')


line00b = ax00.errorbar(energy_efficiency_bin_position, anue_rate_efficiency,
                        xerr=efficiency_x_err, yerr=efficiency_anue_y_err, fmt='o', color='slategray', ecolor='slategray', linewidth=2.0, markersize=6.0,
                        capsize=5.0, markeredgewidth=2, label='xsec_point_2', markeredgecolor='slategray')

ax00.legend((line00a, line00b),
            (r'$\nu_{e}$ Selection Rate', r'$\bar{\nu}_{e}$ Selection Rate'), numpoints=1, fontsize=16.0)

ax00.set_xlabel('Neutrino Energy [GeV]')
ax00.set_ylabel(r'$\nu_{e}/\bar{\nu_{e}}$ Selected / GeV / 6e20 POT')
plt.xlim(0.0, 4.1)
ax00.set_ylim(0, 250)
plt.show()

#############################################################################
# figure001 - as 01, but with fluctuated up low Energy flux

efficiency_nue_up_y_err = [0] * 6

efficiency_nue_up_y_err[0] = efficiency_y_err[0] * nue_rate_efficiency_up[0]
efficiency_nue_up_y_err[1] = efficiency_y_err[1] * nue_rate_efficiency_up[1]
efficiency_nue_up_y_err[2] = efficiency_y_err[2] * nue_rate_efficiency_up[2]
efficiency_nue_up_y_err[3] = efficiency_y_err[3] * nue_rate_efficiency_up[3]
efficiency_nue_up_y_err[4] = efficiency_y_err[4] * nue_rate_efficiency_up[4]
efficiency_nue_up_y_err[5] = efficiency_y_err[5] * nue_rate_efficiency_up[5]

efficiency_anue_up_y_err = [0] * 6

efficiency_anue_up_y_err[0] = efficiency_y_err[0] * anue_rate_efficiency_up[0]
efficiency_anue_up_y_err[1] = efficiency_y_err[1] * anue_rate_efficiency_up[1]
efficiency_anue_up_y_err[2] = efficiency_y_err[2] * anue_rate_efficiency_up[2]
efficiency_anue_up_y_err[3] = efficiency_y_err[3] * anue_rate_efficiency_up[3]
efficiency_anue_up_y_err[4] = efficiency_y_err[4] * anue_rate_efficiency_up[4]
efficiency_anue_up_y_err[5] = efficiency_y_err[5] * anue_rate_efficiency_up[5]

fig, ax000 = plt.subplots()

line000a = ax000.errorbar(energy_efficiency_bin_position, nue_rate_efficiency,
                          xerr=efficiency_x_err, yerr=efficiency_nue_y_err, fmt='o', color='goldenrod', ecolor='goldenrod', linewidth=2.0, markersize=6.0,
                          capsize=5.0, markeredgewidth=2, label='xsec_point_2', markeredgecolor='goldenrod')

line000b = ax000.errorbar(energy_efficiency_bin_position, anue_rate_efficiency,
                          xerr=efficiency_x_err, yerr=efficiency_anue_y_err, fmt='o', color='slategray', ecolor='slategray', linewidth=2.0, markersize=6.0,
                          capsize=5.0, markeredgewidth=2, label='xsec_point_2', markeredgecolor='slategray')


line000c = ax000.errorbar(energy_efficiency_bin_position, nue_rate_efficiency_up,
                          xerr=efficiency_x_err, yerr=efficiency_nue_up_y_err, fmt='o', color='darkgoldenrod', ecolor='darkgoldenrod', linewidth=2.0, markersize=6.0,
                          capsize=5.0, markeredgewidth=2, label='xsec_point_2', markeredgecolor='darkgoldenrod')

line000d = ax000.errorbar(energy_efficiency_bin_position, anue_rate_efficiency_up,
                          xerr=efficiency_x_err, yerr=efficiency_anue_up_y_err, fmt='o', color='darkslategray', ecolor='darkslategray', linewidth=2.0, markersize=6.0,
                          capsize=5.0, markeredgewidth=2, label='xsec_point_2', markeredgecolor='darkslategray')


ax000.legend((line000a, line000b, line000c, line000d),
             (r'$\nu_{e}$ Selection Rate', r'$\bar{\nu}_{e}$ Selection Rate', r'$\nu_{e}$ Selection Rate Enhanced Low E', r'$\bar{\nu}_{e}$ Selection Rate Enhanced Low E'), numpoints=1, fontsize=16.0)

ax000.set_xlabel('Neutrino Energy [GeV]')
ax000.set_ylabel(r'$\nu_{e}/\bar{\nu_{e}}$ Selected / GeV / 6e20 POT')
plt.xlim(0.0, 4.1)
ax000.set_ylim(0, 250)
plt.show()

##############################################################################
# figure1

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
ax1.set_xlabel('Neutrino Energy [GeV]', fontsize=18)
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r'$\nu_{e}/\bar{\nu_{e}}$  / cm$^{2}$ / 6e20 POT', fontsize=18)
# ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
line2a = ax2.plot(x_arr_nue, y_arr_nue, 'mediumblue',
                  linewidth=2.4, label=r'GENIE $\nu_{e}$ Cross Section', linestyle=':')
line2b = ax2.plot(x_arr_anue, y_arr_anue, 'darkgreen',
                  linewidth=2.4, label=r'GENIE $\bar{\nu}_{e}$ Cross Section', linestyle=':')
ax2.set_ylabel(
    r'$\nu_{e}/\bar{\nu}_{e}$ CC Cross Section [cm$^{2}$]', fontsize=18)
# ax2.tick_params('y', colors='r')

ax3 = ax2
genie_xsec_point = 4.83114e-39
data_xsec_point = 4.67e-39
# genie_xsec_point = 4.93656e-39
# for now we're taking the Stat from data
genie_xsec_point_stat_err = 0.057e-39
data_xsec_point_stat_err = 1.01e-39
genie_xsec_point_sys_err = 0.0
data_xsec_point_sys_err = 0.0
genie_xsec_point_total_err = np.sqrt(
    genie_xsec_point_stat_err**2 + genie_xsec_point_sys_err**2)
data_xsec_point_total_err = np.sqrt(
    data_xsec_point_stat_err**2 + data_xsec_point_sys_err**2)
genie_xsec_energy = average
data_xsec_energy = average
asymmetric_x_err = np.array([[average - min_bin_val, max_bin_val - average]]).T
genie_symmetric_y_err = np.array(
    [[genie_xsec_point_total_err, genie_xsec_point_total_err]]).T
data_symmetric_y_err = np.array(
    [[data_xsec_point_total_err, data_xsec_point_total_err]]).T

genie_xsec_point_x_err = genie_xsec_energy * 0.2

# line3a = ax3.plot(genie_xsec_energy, genie_xsec_point,
#                   '#AA66CC', label=r'MC \sigma_{\nu_{e} + \bar{\nu}_{e}}')
# line3b = ax3.errorbar(genie_xsec_energy, genie_xsec_point,
#                       xerr=asymmetric_x_err, yerr=genie_symmetric_y_err, fmt='--o',
#                       color='#AA66CC', ecolor='#AA66CC', linewidth=1.0, markersize=6.0,
#                       label='mc_xsec_point')

line3c = ax3.plot(data_xsec_energy, data_xsec_point,
                  'black', label=r'Data $\sigma_{\nu_{e} + \bar{\nu}_{e}}$')
line3d = ax3.errorbar(data_xsec_energy, data_xsec_point,
                      xerr=asymmetric_x_err, yerr=data_symmetric_y_err, fmt='--o',
                      color='black', ecolor='black', linewidth=1.0, markersize=6.0,
                      label='data_xsec_point')

lines = line1a + line1b + line3c + line2a + line2b
labels = [l.get_label() for l in lines]
# ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)
ax3.legend(lines, labels, loc=9, ncol=2, prop={'size': 16})

ax2.set_ylim(0, 30e-39)
plt.xlim(0.0, 3.0)
plt.show()

###################
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
genie_xsec_point_nue_nue_bar = (4.83e-39, 4.83e-39)

x_edge = (-0.9, 1.1)

# point_1a = ax_2.plot(0.5, genie_xsec_point, 'black',
#                     label=r'MC $\nu_{e}$ + $\bar{\nu}_{e}$')
# point_1b = ax_2.errorbar(0.5, genie_xsec_point, yerr=symmetric_y_err,
# fmt='--o', color='black', ecolor='black', linewidth=1.0, markersize=6.0)

data_xsec_point_list = (4.67e-39, 4.67e-39)
data_xsec_point = 4.67e-39

# lower_1a = genie_xsec_point + (genie_xsec_point * 0.14)
# upper_1a = genie_xsec_point - (genie_xsec_point * 0.14)

# stats uncertainty
data_stat_err = 0.2167
data_sys_err = 0.30
lower_1a = data_xsec_point + (data_xsec_point * data_stat_err)
upper_1a = data_xsec_point - (data_xsec_point * data_stat_err)

lower_1a_list = (lower_1a, lower_1a)
upper_1a_list = (upper_1a, upper_1a)

quadrature = np.sqrt(data_stat_err**2 + data_sys_err**2)

lower_1b = data_xsec_point + (data_xsec_point * quadrature)
upper_1b = data_xsec_point - (data_xsec_point * quadrature)

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

# lines_2 = band_1a_lower + band_1b_lower + line_2a + line_2b + line_2c
lines_2 = (dummy_hist_1a, dummy_hist_1b, line_2a, line_2b, line_2c)
# labels_2 = [l.get_label() for l in lines_2]
# ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)
labels_list = (r'Expected Total Bound', r'Expected Stats Bound',
               r'Genie $\nu_{e}$ Cross Section', r'Genie $\bar{\nu}_{e}$ Cross Section', r'Genie $\nu_{e} + \bar{\nu}_{e}$ Cross Section')
ax_2.legend(loc=9, ncol=2)

plt.show()

##############################################
##############################################
# figure 2b
fig2b, ax_2b = plt.subplots()

# these are to make the legend have a fill, rather than line
dummy_hist_1a = ax_2b.hist(
    (-1, -2), 1, fill=True, label=r'Total Bound', color='darkgray', alpha=0.2)
# ax.legend(markerscale=2,fontsize=20)
dummy_hist_1b = ax_2b.hist(
    (-1, -2), 1, fill=True, label=r'Stats Bound', color='maroon', alpha=0.2)

band_1a_lower = ax_2b.plot(x_edge, lower_1a_list,
                           color='maroon')
band_1a_upper = ax_2b.plot(x_edge, upper_1a_list, color='maroon')

band_1b_lower = ax_2b.plot(x_edge, lower_1b_list,
                           color='black')
band_1b_upper = ax_2b.plot(x_edge, upper_1b_list, color='black')

ax_2b.fill_between(x_edge, upper_1a_list, lower_1a_list, hatch='//',
                   facecolor='maroon', alpha=0.2)

ax_2b.fill_between(x_edge, upper_1b_list, lower_1b_list,
                   facecolor='darkgray', alpha=0.2)

line_2a = ax_2b.plot(x_edge, genie_xsec_point_nue,
                     color='mediumblue', linewidth=1.5, label=r'Genie $\nu_{e}$ Cross Section', linestyle='--')
line_2b = ax_2b.plot(x_edge, genie_xsec_point_nue_bar,
                     color='darkgreen', linewidth=1.5, label=r'Genie $\bar{\nu}_{e}$ Cross Section', linestyle='--')
line_2c = ax_2b.plot(x_edge, genie_xsec_point_list,
                     color='cadetblue', linewidth=1.5, label=r'Genie $\nu_{e} + \bar{\nu}_{e}$ Cross Section', linestyle='--')

line_2d = ax_2b.plot(x_edge, data_xsec_point_list, color='black', linewidth=1.8,
                     label=r'Data $\nu_{e} + \bar{\nu}_{e}$ Cross Section')

ax_2b.set_ylim(0, 10e-39)
ax_2b.set_xlim(0, 1)
ax_2b.set_ylabel(
    r'$\nu_{e}/\bar{\nu}_{e}$ CC Cross Section [cm$^{2}$]', fontsize=18)
ax_2b.xaxis.label.set_size(0)
ax_2b.xaxis.set_ticklabels([])
ax_2b.xaxis.set_visible(False)

# lines_2 = band_1a_lower + band_1b_lower + line_2a + line_2b + line_2c
# lines_list = (dummy_hist_1a, dummy_hist_1b, line_2a, line_2b, line_2c)
# labels_2 = [l.get_label() for l in lines_2]
# ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)
# labels_list = (r'Expected Total Bound', r'Expected Stats Bound',
# r'Genie $\nu_{e}$ Cross Section', r'Genie $\bar{\nu}_{e}$ Cross
# Section', r'Genie $\nu_{e} + \bar{\nu}_{e}$ Cross Section')
ax_2b.legend(loc=9, ncol=2)

plt.show()


##############################################################################

# fig3, ax3 = plt.subplots()
# # line0 = ax1.plot(fine_flux_bins, fine_flux_list, 'black', linewidth=2.0)
# line3a = ax3.plot(bin_flux_list, nue_flux_list, 'mediumblue',
#                   linewidth=2.0, label=r'NuMI $\nu_{e}$ Flux')
# ax3.fill_between(bin_flux_list, nue_flux_list, 0,
#                  facecolor='darkgray', alpha=0.2)
# line3b = ax3.plot(a_bin_flux_list, anue_flux_list, 'darkgreen',
#                   linewidth=2.0, label=r'NuMI $\bar{\nu}_{e}$ Flux')
# ax3.fill_between(a_bin_flux_list, anue_flux_list,
#                  0, facecolor='darkgray', alpha=0.2)
# ax3.set_xlabel('Neutrino Energy [GeV]')
# # make the y-axis label, ticks and tick labels match the line color.
# ax3.set_ylabel(r'$\nu_{e}/\bar{\nu_{e}}$  / cm$^{2}$ / 6e20 POT')
# # ax1.tick_params('y', colors='b')
#
# ax4 = ax3.twinx()
# line4a = ax4.plot(x_arr_nue, y_arr_nue, 'mediumblue',
#                   linewidth=2.4, label=r'GENIE $\nu_{e}$ Cross Section', linestyle=':')
# line4b = ax4.plot(x_arr_anue, y_arr_anue, 'darkgreen',
#                   linewidth=2.4, label=r'GENIE $\bar{\nu}_{e}$ Cross Section', linestyle=':')
# ax4.set_ylabel(r'$\nu_{e}/\bar{\nu}_{e}$ CC Cross Section [cm$^{2}$]')
# # ax2.tick_params('y', colors='r')
#
# ax5 = ax4
# # genie_xsec_point = 4.83114e-39
# # for now we're taking the Stat from data
# # genie_xsec_point_stat_err = 0.703e-39
# # genie_xsec_point_sys_err = 0.25 * genie_xsec_point
# # genie_xsec_point_total_err = np.sqrt(
# #    genie_xsec_point_stat_err**2 + genie_xsec_point_sys_err**2)
# # genie_xsec_energy = average
# # asymmetric_x_err = np.array([[average - min_bin_val, max_bin_val - average]]).T
# # symmetric_y_err = np.array(
# #    [[genie_xsec_point_total_err, genie_xsec_point_total_err]]).T
#
# # genie_xsec_point_x_err = genie_xsec_energy * 0.2
#
# line5a = ax4.plot(genie_xsec_energy, genie_xsec_point,
#                   'black', label=r'$\nu_e$ + $\bar{\nu}_{e}$ MC Cross Section', linewidth=2.0)
# line5b = ax4.errorbar(genie_xsec_energy, genie_xsec_point,
#                       xerr=asymmetric_x_err, yerr=genie_symmetric_y_err, fmt='--o',
#                       color='black', ecolor='black', linewidth=2.0, markersize=6.0, capsize=5.0, markeredgewidth=2,
#                       label='xsec_point')
#
# symmetric_y_err_2 = np.array(
#     [[genie_xsec_point_stat_err, genie_xsec_point_stat_err]]).T
# line6a = ax4.errorbar(genie_xsec_energy, genie_xsec_point,
#                       xerr=asymmetric_x_err, yerr=symmetric_y_err_2, fmt='--o', color='black', ecolor='black', linewidth=2.0, markersize=6.0,
#                       capsize=5.0, markeredgewidth=2, label='xsec_point_2', linestyle='--')
#
# lines = line3a + line3b + line5a + line4a + line4b
# labels = [l.get_label() for l in lines]
# # ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
# #           ncol=2, mode="expand", borderaxespad=0.)
# ax5.legend(lines, labels, loc=9, ncol=2)
#
# ax4.set_ylim(0, 35e-39)
# plt.xlim(0.0, 3.0)
#
# plt.show()
