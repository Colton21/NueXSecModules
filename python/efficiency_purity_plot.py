import matplotlib.pyplot as plt

file_string = "../scripts/selection_output.txt"

file = open(file_string, "r")
print file_string, " is open!"

cut_list = []
num_list = []
eff_list = []
pur_list = []
eff_pur_list = []

# should be size 15
# this is purity without Off-Beam

eff_list_short = (0.5489, 0.38, 0.315, 0.261, 0.132, 0.0903)
pur_list_short = (0.00736, 0.03569, 0.0418, 0.1514, 0.2928, 0.3973)
num_list_edit = (1, 2, 3, 4, 5, 6)

pur_nu_only_list = (0.1225, 0.1612, 0.1943, 0.3253, 0.5142, 0.6978)

cut_list_edit = (r'Simple Pre-Cuts (1)', r'Flash-Matching (2)', r'Reco. Vertex Quality (3)',
                 r'Shower Hit Threshold (4)', r'Electron-like Shower (5)', r'Finalisation Cuts (6)')

# pur_nu_only_list = (0.00329, 0.00338, 0.0244, 0.0427, 0.1257, 0.1451,
# 0.1488, 0.272, 0.2927, 0.3375, 0.4986, 0.5765, 0.6091, 0.6259, 0.6530)

# cut_list_edit = (r'In Time', r'PE Threshold', r'Reco $\nu_{e}$', r'In Fid. Vol.', r'Vertex to Flash', r'Shower to $\nu$', r'Track to $\nu$', r'Hit Threshold',
#                 r'YPlane Hit Threshold', r'Opening Angle', r'dE/dx', r'Secondary Shower Dist.', r'Hit/Length Ratio',
#                 r'Track Len / Shower Len Ratio', r'Track Contained')

num = 1
for line in file:
    line_info = line.split(",")
    efficiency = float(line_info[1])
    purity = float(line_info[2])
    cut_list.append(line_info[0])
    eff_list.append(efficiency)
    pur_list.append(purity)
    num_list.append(num)
    eff_pur_list.append(efficiency * purity)
    num = num + 1

file.close()

fig, ax1 = plt.subplots()
line1a = ax1.plot(num_list_edit, eff_list_short, 'darkcyan',
                  linewidth=2.0, label=r'Efficiency')
line1b = ax1.plot(num_list_edit, pur_list_short, 'royalblue',
                  linewidth=2.0, label=r'Purity')

ax1.set_xlabel('Selection Cut')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r'Selection Performance')
# ax1.tick_params('y', colors='b')

#ax2 = ax1.twinx()
# line1c = ax2.plot(num_list, eff_pur_list, 'blueviolet',
#                  linewidth=2.0, label=r'Efficiency x Purity')
#ax2.set_ylabel(r'Efficiency x Purity')
# line2a = ax1.plot(num_list, pur_list, 'royalblue',
#                  linewidth=2.0, label=r'Purity')

line1d = ax1.plot(num_list_edit, pur_nu_only_list, 'blueviolet',
                  linewidth=2.0, label=r'Purity (Beam Only)')

lines = line1a + line1b + line1d
labels = [l.get_label() for l in lines]
# ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)
ax1.legend(lines, labels, loc=9, ncol=3)

ax1.set_ylim(0, 1.0)
#ax2.set_ylim(0, 0.1)
#plt.xlim(0.0, 3.0)
#plt.xticks(num_list, cut_list)
# for tick in ax1.get_xticklabels():
#    tick.set_rotation(45)

ax1.set_xticks(num_list_edit)
ax1.set_xticklabels(cut_list_edit, rotation=40, ha='right')

plt.grid(True)
plt.show()
