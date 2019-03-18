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
           rotation=90, horizontalalignment='center', fontsize=18)


# ax0.fill_between(bin_flux_list, plotting_prod_nue_flux_xsec_list, 0,
# facecolor='darkgray', alpha=0.2)

ax0.set_xlabel('Detector Variation', fontsize=20, horizontalalignment='center')
ax0.set_ylabel(
    r'$\sigma_{var}$ - $\sigma_{CV}$ / $\sigma_{CV}$ [%]', fontsize=20)

ax0.set_ylim(0, 14.0)
# plt.xlim(0.0, 4.0)
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
           rotation=90, horizontalalignment='center', fontsize=18)


# ax0.fill_between(bin_flux_list, plotting_prod_nue_flux_xsec_list, 0,
# facecolor='darkgray', alpha=0.2)

ax1.set_xlabel('Detector Variation', fontsize=20, horizontalalignment='center')
ax1.set_ylabel(
    r'$\epsilon_{var}$ - $\epsilon_{CV}$ / $\epsilon_{CV}$ [%]', fontsize=20)

ax1.set_ylim(0, 14)
# plt.xlim(0.0, 4.0)
# ax0.align_xlabels()  # same as fig.align_xlabels(); fig.align_ylabels()

plt.show()


##########
# statistical uncertainty

n_mc_bkg = 356.
n_mc_dirt = 30.
n_ext = 81.
n_data = 214.

# statistical_error = np.sqrt((np.sqrt(n_data) / n_data)**2 + (np.sqrt(n_mc_dirt) /
# n_mc_dirt)**2 + (np.sqrt(n_ext) / n_ext)**2 + (np.sqrt(n_mc_bkg) /
# n_mc_bkg)**2)

statistical_error = np.sqrt(
    n_data + n_ext * (1.0154)**2 + n_mc_bkg * (0.1301)**2 + n_mc_dirt * (0.16411)**2)

print 'N_data ', n_data, ' sqrt(N_data) ', np.sqrt(n_data), ' sqrt(N_data) / N_data ', np.sqrt(n_data) / n_data
print 'N_ext ', n_ext, ' sqrt(N_ext) ', np.sqrt(n_ext), ' sqrt(N_ext) / N_ext ', np.sqrt(n_ext) / n_ext
print 'N_dirt ', n_mc_dirt, ' sqrt(N_dirt) ', np.sqrt(n_mc_dirt), ' sqrt(N_dirt) / N_dirt ', np.sqrt(n_mc_dirt) / n_mc_dirt
print 'N_bkg ', n_mc_bkg, ' sqrt(N_bkg) ', np.sqrt(n_mc_bkg), ' sqrt(N_mc_bkg) / N_mc_bkg ', np.sqrt(n_mc_bkg) / n_mc_bkg
print 'Total Sats Error (alt): ', np.sqrt((np.sqrt(n_data) / n_data)**2 +
                                          (np.sqrt(n_ext) / n_ext)**2 +
                                          (np.sqrt(n_mc_dirt) / n_mc_dirt)**2 +
                                          (np.sqrt(n_mc_bkg) / n_mc_bkg)**2)
print '--------------------------'
print 'sum ', np.sqrt(n_data + n_ext)
print 'N_mc_bkg ', n_mc_bkg, ' and N_mc_bkg * scaling^2: ', n_mc_bkg * (0.1301)**2
print 'N_mc_dirt ', n_mc_dirt, ' and N_mc_dirt * scaling^2: ', n_mc_dirt * (0.16411)**2

print 'Number of events error for N-B: ', statistical_error

statistical_error = statistical_error * 100. / 80.51  # from N - B
print 'Statistical Error: ', statistical_error, '%'

# efficiency error

n_mc_signal = 728500.
efficiency = (619. / 7103.)

efficiency_error = (1. / np.sqrt(n_mc_signal)) * \
    np.sqrt(efficiency * (1 - efficiency))
efficiency_error = efficiency_error * 100.
print 'Efficiency Error: ', efficiency_error, '%'

#############################################################################
# flux uncertainty - beamline variations
# these are a series of uni-sims

beamline_error = 0.

nominal = 0.0

# errors are as a %
beam_spot_1_1mm = (0.0 - nominal) / nominal
beam_spot_1_5mm = (0.0 - nominal) / nominal
horns_water_0mm = (0.0 - nominal) / nominal
horns_water_2mm = (0.0 - nominal) / nominal
old_horn = (0.0 - nominal) / nominal

horn_plus_2ka = (0.0 - nominal) / nominal
horn_minus_2ka = (0.0 - nominal) / nominal

horn1_x_plus_3mm = (0.0 - nominal) / nominal
horn1_x_minus_3mm = (0.0 - nominal) / nominal

horn1_y_plus_3mm = (0.0 - nominal) / nominal
horn1_y_minus_3mm = (0.0 - nominal) / nominal

horn2_x_plus_3mm = (0.0 - nominal) / nominal
horn2_x_minus_3mm = (0.0 - nominal) / nominal

horn2_y_plus_3mm = (0.0 - nominal) / nominal
horn2_y_minus_3mm = (0.0 - nominal) / nominal

beam_shift_x_plus_1mm = (0.0 - nominal) / nominal
beam_shift_x_minus_1mm = (0.0 - nominal) / nominal


# calculate the detector error
horn_2ka = 0.
horn1_x_3mm = 0.
horn1_y_3mm = 0.
horn2_x_3mm = 0.
horn2_y_3mm = 0.
beam_shift_x_1mm = 0.

if(horn_plus_2ka > horn_minus_2ka):
    horn_2ka = horn_plus_2ka
if(horn_plus_2ka < horn_minus_2ka):
    horn_2ka = horn_minus_2ka

if(horn1_x_plus_3mm > horn1_x_minus_3mm):
    horn1_x_3mm = horn1_x_plus_3mm
if(horn1_x_plus_3mm < horn1_x_minus_3mm):
    horn1_x_3mm = horn1_x_minus_3mm

if(horn1_y_plus_3mm > horn1_y_minus_3mm):
    horn1_y_3mm = horn1_y_plus_3mm
if(horn1_y_plus_3mm < horn1_y_minus_3mm):
    horn1_y_3mm = horn1_y_minus_3mm

if(horn2_x_plus_3mm > horn2_x_minus_3mm):
    horn2_x_3mm = horn2_x_plus_3mm
if(horn2_x_plus_3mm < horn2_x_minus_3mm):
    horn2_x_3mm = horn2_x_minus_3mm

if(horn2_y_plus_3mm > horn2_y_minus_3mm):
    horn2_y_3mm = horn2_y_plus_3mm
if(horn2_y_plus_3mm < horn2_y_minus_3mm):
    horn2_y_3mm = horn2_y_minus_3mm

if(beam_shift_x_plus_1mm > beam_shift_x_minus_1mm):
    beam_shift_x_1mm = beam_shift_x_plus_1mm
if(beam_shift_x_plus_1mm < beam_shift_x_minus_1mm):
    beam_shift_x_1mm = beam_shift_x_minus_1mm

beamline_error = np.sqrt(beam_spot_1_1mm**2 + beam_spot_1_5mm**2 + horns_water_0mm**2 + horns_water_2mm**2 +
                         old_horn**2 + horn_2ka**2 + horn1_x_3mm**2 + horn1_y_3mm**2 +
                         horn2_x_3mm**2 + horn2_y_3mm**2 + beam_shift_x_1mm**2)

print 'Beamline Variation Uncertainty: ', beamline_error, '%'

beamline_error_list = (round(beam_spot_1_1mm, 2), round(beam_spot_1_5mm, 2), round(horns_water_0mm, 2), round(horns_water_2mm, 2),
                       round(old_horn, 2),
                       round(horn_2ka, 2),
                       round(horn1_x_3mm, 2),
                       round(horn1_y_3mm, 2),
                       round(horn2_x_3mm, 2),
                       round(horn2_y_3mm, 2),
                       round(beam_shift_x_1mm, 2),
                       round(beamline_error, 2))

beamline_bin_list = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
beamline_labels = ("Beam Spot 1mm", "Beam Spot 5mm", "Horn Water 0mm", "Horn Water 2mm", "Old Horn",
                   "Horn 2kA", "Horn1 X 3mm", "Horn1 Y 3mm", "Horn2 X 3mm", "Horn2 Y 3mm",
                   "Beam Shift X 1mm", "Total")

fig2, ax2 = plt.subplots()

line2a = ax2.bar(beamline_bin_list, beamline_error_list,
                 color='darkgray', width=0.9)

for i, v in enumerate(beamline_error_list):
    ax2.text(i + 1.2, v / beamline_error_list[i] + 0.005, beamline_error_list[i],
             fontsize=16, color='salmon')

plt.xticks(np.arange(1.5, 12.5, step=1.0), labels,
           rotation=90, horizontalalignment='center', fontsize=18)


# ax0.fill_between(bin_flux_list, plotting_prod_nue_flux_xsec_list, 0,
# facecolor='darkgray', alpha=0.2)

ax0.set_xlabel('Beamline Variation', fontsize=20, horizontalalignment='center')
ax0.set_ylabel(
    r'$\sigma_{var}$ - $\sigma_{CV}$ / $\sigma_{CV}$ [%]', fontsize=20)

ax0.set_ylim(0, 20.0)
# plt.xlim(0.0, 4.0)
# ax0.align_xlabels()  # same as fig.align_xlabels(); fig.align_ylabels()

plt.show()
