import matplotlib.pyplot as plt

file_string = "../scripts/selection_output.txt"

file = open(file_string, "r")
print file_string, " is open!"

cut_list = []
num_list = []
eff_list = []
pur_list = []
eff_pur_list = []

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
line1a = ax1.plot(num_list, eff_list, 'darkcyan',
                  linewidth=2.0, label=r'Efficiency')
line1b = ax1.plot(num_list, pur_list, 'royalblue',
                  linewidth=2.0, label=r'Purity')

ax1.set_xlabel('Selection Cut')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r'Selection Performance')
# ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
line1c = ax2.plot(num_list, eff_pur_list, 'blueviolet',
                  linewidth=2.0, label=r'Efficiency x Purity')
ax2.set_ylabel(r'Efficiency x Purity')
# line2a = ax2.plot(num_list, pur_list, 'royalblue',
#                  linewidth=2.0, label=r'Purity')

lines = line1a + line1b + line1c
labels = [l.get_label() for l in lines]
# ax3.legend(lines, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)
ax1.legend(lines, labels, loc=9, ncol=3)

ax1.set_ylim(0, 1.0)
ax2.set_ylim(0, 0.1)
#plt.xlim(0.0, 3.0)
#plt.xticks(num_list, cut_list)
# for tick in ax1.get_xticklabels():
#    tick.set_rotation(45)

ax1.set_xticks(num_list)
ax1.set_xticklabels(cut_list, rotation=40, ha='right')

plt.show()
