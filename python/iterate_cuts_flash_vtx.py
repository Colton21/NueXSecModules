# '**********************************'
# 'Script to generate a custom cut configuration file' and run selection
# '***********************************'

# ************************************
# HOW TO USE: execute this script as a normal python script
# so: python run_selection.py MC_file EXT_file Data_file
# Including any file is optional
# ************************************

# This version of the python script runs the slimmed version of the
# selection by default --> this means no histograms will be generated

# writes variable as string to file


def write_to_file(file, var):
    file.write(str(var))
    file.write('\n')

from subprocess import call
import sys
import matplotlib
import matplotlib.pyplot as plt

# need clean output file for every run
call(["rm", "-rf", "selection_output.txt"])


# Now for the code which calls the c++ code multiple times with a varied
# input parameter

num_cuts = 21
run = 0

loop_cut = []

for cut in range(num_cuts):

    f = open('config.txt', 'w')

    # write the value configurations here
    # if you're looking for the default parameters, and some more explanation,
    # see main.h

    # configurable values:

    # fiducial volume values
    x1 = 40.0
    x2 = 40.0
    y1 = 40.0
    y2 = 40.0
    z1 = 40.0
    z2 = 40.0

    write_to_file(f, x1)
    write_to_file(f, x2)
    write_to_file(f, y1)
    write_to_file(f, y2)
    write_to_file(f, z1)
    write_to_file(f, z2)

    # optical cuts
    flash_pe_threshold = 50.0
    flash_time_start = 5.0
    flash_time_end = 16.0

    write_to_file(f, flash_pe_threshold)
    write_to_file(f, flash_time_start)
    write_to_file(f, flash_time_end)

    # 2D distance between nu-vtx and largest flash centre
    tolerance = 80.0

    #***************************************************
    # tolerance is a value we'd like to investigate:
    tolerance = 60.0 + (run * 2)
    # threshold - scan between 60.0 - 100.0 PE
    loop_cut.append(tolerance)

    write_to_file(f, tolerance)

    # reco track/shower to nu-vtx distance
    shwr_nue_tolerance = 4.0
    trk_nue_tolerance = 4.0

    write_to_file(f, shwr_nue_tolerance)
    write_to_file(f, trk_nue_tolerance)

    # leading shower hit thresholds
    shwr_hit_threshold = 200.0
    shwr_hit_threshold_collection = 80.0

    write_to_file(f, shwr_hit_threshold)
    write_to_file(f, shwr_hit_threshold_collection)

    # leading shower opening angle
    tolerance_open_angle_min = 2.0
    tolerance_open_angle_max = 15.0

    write_to_file(f, tolerance_open_angle_min)
    write_to_file(f, tolerance_open_angle_max)

    # leading shower dE/dx
    tolerance_dedx_min = 1.4
    tolerance_dedx_max = 3.0

    write_to_file(f, tolerance_dedx_min)
    write_to_file(f, tolerance_dedx_max)

    # maximum distance secondary shower from nu-vtx
    # 22 cm is roughly two radiation lenghts, which should still keep Pi0
    # events
    dist_tolerance = 22.0

    write_to_file(f, dist_tolerance)

    # ratio of hits / length (hits / cm)
    pfp_hits_length_tolerance = 3.0

    write_to_file(f, pfp_hits_length_tolerance)

    # ratio of longest track / leading shower lengths
    ratio_tolerance = 1.0

    write_to_file(f, ratio_tolerance)

    # not listed here as there's no configuration - also require that all
    # tracks are contained!

    # close input file
    f.close()

    # script name + MC name
    if(len(sys.argv) == 2):
        call(["../scripts/main.exe", sys.argv[1], "-c", "config.txt", "--slim"])

    # script name + MC name + EXT name
    if(len(sys.argv) == 3):
        call(["../scripts/main.exe", sys.argv[1],
              sys.argv[2], "-c", "config.txt", "--slim"])

    # script name + MC name + EXT name
    if(len(sys.argv) == 4):
        call(["../scripts/main.exe", sys.argv[1],
              sys.argv[2], sys.argv[3], "-c", "config.txt", "--slim"])

    run = run + 1
# END for loop

# After the function runs - run the rest of the python script

# we only need read access to this file
purity = []
efficiency = []

purity_final = []
efficiency_final = []

import csv
with open('selection_output.txt', 'rb') as out_file:
    reader = csv.reader(out_file)
    for row in reader:
        if(row[0] == "Vtx-to-Flash"):
            efficiency.append(float(row[1]) * 100.0)
            purity.append(float(row[2]) * 100.0)
        if(row[0] == "Track Containment"):
            efficiency_final.append(float(row[1]) * 100.0)
            purity_final.append(float(row[2]) * 100.0)


figure_1 = plt.figure()
plt.plot(loop_cut, efficiency_final, 'ko')
f1_axes = figure_1.add_subplot(111)
f1_axes.set_title("Efficiency vs Vtx-to-Flash")
f1_axes.set_autoscaley_on(False)
f1_axes.set_ylim([7, 11])
f1_axes.set_xlabel("Vtx-to-Flash [cm]")
f1_axes.set_xlim([58, 102])
f1_axes.set_ylabel("Efficiency")
figure_1.savefig("plots/final_efficiency_flash_vtx.pdf")

figure_2 = plt.figure()
plt.plot(loop_cut, purity_final, 'bo')
f2_axes = figure_2.add_subplot(111)
f2_axes.set_title("Purity vs Vtx-to-Flash")
f2_axes.set_autoscaley_on(False)
f2_axes.set_ylim([50, 65])
f2_axes.set_xlabel("Vtx-to-Flash [cm]")
f2_axes.set_xlim([58, 102])
f2_axes.set_ylabel("Purity")
figure_2.savefig("plots/final_purity_flash_vtx.pdf")


# then delete the file, the c++ code appends values so needs to be fresh
# every time this script is run
#call(["rm", "-rf", "selection_output.txt"])
# this is done at the start now
