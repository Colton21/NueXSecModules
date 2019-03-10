import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# use this module to quickly calculate some errors


stats_error = 0.
detector_error = 0.

# errors are as a %
dataSCE = 3.38
LArG4BugFix = 1.38
DLup = 3.97
DLdown = 2.56
DTup = 1.14
DTdown = 0.83
noiseAmpUp = 1.22
noiseAmpDown = 2.04
withDIC = 6.25
squeezeResp = 2.81
stretchResp = 1.83
altDeadChannels = 0.354
deadSaturatedChannels = 5.90
TPCVisUp = 0.413
BirksRecomb = 3.99
UpPEnoise = 2.75
DownPEnoise = 1.61

# calculate the detector error
DL = 0.
DT = 0.
noiseAmp = 0.
PEnoise = 0.

if(DLup > DLdown):
    DL = DLup
if(DLup < DLdown):
    DL = DLdown

if(DTup > DTdown):
    DT = DTup
if(DTup < DTdown):
    DT = DTdown

if(noiseAmpUp > noiseAmpDown):
    noiseAmp = noiseAmpUp
if(noiseAmpUp < noiseAmpDown):
    noiseAmp = noiseAmpDown

if(UpPEnoise > DownPEnoise):
    PEnoise = UpPEnoise
if(UpPEnoise < DownPEnoise):
    PEnoise = DownPEnoise

detector_error = np.sqrt(dataSCE**2 + LArG4BugFix**2 + DL**2 + DT**2 +
                         noiseAmp**2 + withDIC**2 + squeezeResp**2 + stretchResp**2 +
                         altDeadChannels**2 + deadSaturatedChannels**2 + TPCVisUp**2 +
                         BirksRecomb**2 + PEnoise**2)

print 'Detector Uncertainty: ', detector_error, '%'

detector_error_list = (dataSCE, LArG4BugFix, DL, DT, noiseAmp, PEnoise, withDIC, squeezeResp,
                       stretchResp, altDeadChannels, deadSaturatedChannels, TPCVisUp, BirksRecomb, round(detector_error, 2))

bin_list = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
labels = ("Data-driven SCE", "LArG4 Bug Fix", "Longitudinal Diffusion", "Transverse Diffusion", "Wire Noise",
          "PE Noise", "Induced Charge", "Squeeze Resp.", "Stretch Resp.", "Dead Channels",
          "Saturated Channels", "TPC Vis. Up", "Birks Recombination", "Total")

fig0, ax0 = plt.subplots()

line0a = ax0.bar(bin_list, detector_error_list,
                 color='darkslategray', width=0.9)

for i, v in enumerate(detector_error_list):
    ax0.text(i + 1.2, v / detector_error_list[i] + 0.005, detector_error_list[i],
             fontsize=16, color='salmon')

plt.xticks(np.arange(1.5, 15.5, step=1.0), labels,
           rotation=90, horizontalalignment='center')


# ax0.fill_between(bin_flux_list, plotting_prod_nue_flux_xsec_list, 0,
# facecolor='darkgray', alpha=0.2)

ax0.set_xlabel('Detector Variation', fontsize=18, horizontalalignment='center')
ax0.set_ylabel(
    r'$\sigma_{var}$ - $\sigma_{CV}$ / $\sigma_{CV}$ [%]', fontsize=18)

ax0.set_ylim(0, 14.0)
#plt.xlim(0.0, 4.0)
# ax0.align_xlabels()  # same as fig.align_xlabels(); fig.align_ylabels()

plt.show()

# detector error - efficiency
e_detector_cv = 7.57
e_dataSCE = abs(7.32 - e_detector_cv) * 100.0 / e_detector_cv
e_LArG4BugFix = abs(7.47 - e_detector_cv) * 100.0 / e_detector_cv
e_DLup = abs(7.27 - e_detector_cv) * 100.0 / e_detector_cv
e_DLdown = abs(7.77 - e_detector_cv) * 100.0 / e_detector_cv
e_DTup = abs(7.48 - e_detector_cv) * 100.0 / e_detector_cv
e_DTdown = abs(7.63 - e_detector_cv) * 100.0 / e_detector_cv
e_noiseAmpUp = abs(7.66 - e_detector_cv) * 100.0 / e_detector_cv
e_noiseAmpDown = abs(7.77 - e_detector_cv) * 100.0 / e_detector_cv
e_withDIC = abs(7.11 - e_detector_cv) * 100.0 / e_detector_cv
e_squeezeResp = abs(7.35 - e_detector_cv) * 100.0 / e_detector_cv
e_stretchResp = abs(7.71 - e_detector_cv) * 100.0 / e_detector_cv
e_altDeadChannels = abs(7.59 - e_detector_cv) * 100.0 / e_detector_cv
e_deadSaturatedChannels = abs(7.14 - e_detector_cv) * 100.0 / e_detector_cv
e_TPCVisUp = abs(7.54 - e_detector_cv) * 100.0 / e_detector_cv
e_BirksRecomb = abs(7.27 - e_detector_cv) * 100.0 / e_detector_cv
e_UpPEnoise = abs(7.36 - e_detector_cv) * 100.0 / e_detector_cv
e_DownPEnoise = abs(7.44 - e_detector_cv) * 100.0 / e_detector_cv

# calculate the detector error
e_DL = 0.
e_DT = 0.
e_noiseAmp = 0.
e_PEnoise = 0.

if(e_DLup > e_DLdown):
    e_DL = e_DLup
if(e_DLup < e_DLdown):
    e_DL = e_DLdown

if(e_DTup > e_DTdown):
    e_DT = e_DTup
if(e_DTup < e_DTdown):
    e_DT = e_DTdown

if(e_noiseAmpUp > e_noiseAmpDown):
    e_noiseAmp = e_noiseAmpUp
if(e_noiseAmpUp < e_noiseAmpDown):
    e_noiseAmp = e_noiseAmpDown

if(e_UpPEnoise > e_DownPEnoise):
    e_PEnoise = e_UpPEnoise
if(e_UpPEnoise < e_DownPEnoise):
    e_PEnoise = e_DownPEnoise

e_detector_error = np.sqrt(e_dataSCE**2 + e_LArG4BugFix**2 + e_DL**2 + e_DT**2 +
                           e_noiseAmp**2 + e_withDIC**2 + e_squeezeResp**2 + e_stretchResp**2 +
                           e_altDeadChannels**2 + e_deadSaturatedChannels**2 + e_TPCVisUp**2 +
                           e_BirksRecomb**2 + e_PEnoise**2)

print 'Detector Uncertainty (efficiency): ', e_detector_error, '%'

e_detector_error_list = (round(e_dataSCE, 2), round(e_LArG4BugFix, 2), round(e_DL, 2), round(e_DT, 2), round(e_noiseAmp, 2),
                         round(e_PEnoise, 2), round(
                             e_withDIC, 2), round(e_squeezeResp, 2),
                         round(e_stretchResp, 2), round(e_altDeadChannels, 2), round(
                             e_deadSaturatedChannels, 2),
                         round(e_TPCVisUp, 2), round(e_BirksRecomb, 2), round(e_detector_error, 2))

fig1, ax1 = plt.subplots()

line1a = ax1.bar(bin_list, e_detector_error_list,
                 color='darkslategray', width=0.9)

for i, v in enumerate(e_detector_error_list):
    ax1.text(i + 1.2, v / e_detector_error_list[i] + 0.005, e_detector_error_list[i],
             fontsize=16, color='salmon')

plt.xticks(np.arange(1.5, 15.5, step=1.0), labels,
           rotation=90, horizontalalignment='center')


# ax0.fill_between(bin_flux_list, plotting_prod_nue_flux_xsec_list, 0,
# facecolor='darkgray', alpha=0.2)

ax1.set_xlabel('Detector Variation', fontsize=18, horizontalalignment='center')
ax1.set_ylabel(
    r'$\epsilon_{var}$ - $\epsilon_{CV}$ / $\epsilon_{CV}$ [%]', fontsize=18)

ax1.set_ylim(0, 14)
#plt.xlim(0.0, 4.0)
# ax0.align_xlabels()  # same as fig.align_xlabels(); fig.align_ylabels()

plt.show()


##########
# statistical uncertainty

n_mc_bkg = 329.
n_mc_dirt = 30.
n_ext = 81.
n_data = 203.

# statistical_error = np.sqrt((np.sqrt(n_data) / n_data)**2 + (np.sqrt(n_mc_dirt) /
# n_mc_dirt)**2 + (np.sqrt(n_ext) / n_ext)**2 + (np.sqrt(n_mc_bkg) /
# n_mc_bkg)**2)

statistical_error = np.sqrt(
    n_data + n_ext + n_mc_bkg * (0.1301)**2 + n_mc_dirt * (0.16411)**2)

print 'Number of events error for N-B: ', statistical_error

statistical_error = statistical_error * 100. / 73.026  # from N - B
print 'Statistical Error: ', statistical_error, '%'

# efficiency error

n_mc_signal = 728500.
efficiency = (619. / 7103.)

efficiency_error = (1. / np.sqrt(n_mc_signal)) * \
    np.sqrt(efficiency * (1 - efficiency))
efficiency_error = efficiency_error * 100.
print 'Efficiency Error: ', efficiency_error, '%'

#############################################################################
