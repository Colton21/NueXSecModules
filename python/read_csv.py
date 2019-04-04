import numpy as np

#f = open("../scripts/neutrino_in_tpc_list.txt", "r")
f = open("../scripts/selected_events_for_truth_script.txt", "r")

print f.readline()

strange_list_nue = []
strange_list_nuebar = []

for line in f:
    split = line.split(',')
    energy_str = split[6]
    id_str = split[2]
    x_dir_str = split[3]
    y_dir_str = split[4]
    z_dir_str = split[5]
    energy = float(energy_str[1:])
    id = float(id_str[1:])
    x_dir = float(x_dir_str[1:])
    y_dir = float(y_dir_str[1:])
    z_dir = float(z_dir_str[1:])
    if(energy >= 3.125 and energy <= 3.130):
        if(id == 1):
            strange_list_nue.append([x_dir, y_dir, z_dir, energy])
        if(id == 5):
            strange_list_nuebar.append([x_dir, y_dir, z_dir, energy])

print 'Nue'
print strange_list_nue
print 'Nuebar'
print strange_list_nuebar

#f2 = open("../scripts/selected_events_for_truth_script.txt", "r")
# print f2.readline()
