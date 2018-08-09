# '**********************************'
# 'Script to generate a custom cut configuration file' and run selection
# '***********************************'

# ************************************
# HOW TO USE: execute this script as a normal python script
# so: python run_selection.py MC_file EXT_file Data_file
# Including any file is optional
# ************************************

# writes variable as string to file


def write_to_file(file, var):
    file.write(str(var))
    file.write('\n')

from subprocess import call
import sys

f = open('config.txt', 'w')
f2 = open('photon_config.txt', 'w')

# write the value configurations here
# if you're looking for the default parameters, and some more explanation,
# see main.h

# configurable values:

# fiducial volume values
x1 = 20.0
x2 = 20.0
y1 = 20.0
y2 = 20.0
z1 = 20.0
z2 = 20.0

write_to_file(f, x1)
write_to_file(f, x2)
write_to_file(f, y1)
write_to_file(f, y2)
write_to_file(f, z1)
write_to_file(f, z2)

write_to_file(f2, x1)
write_to_file(f2, x2)
write_to_file(f2, y1)
write_to_file(f2, y2)
write_to_file(f2, z1)
write_to_file(f2, z2)

# optical cuts
flash_pe_threshold = 50.0
flash_time_start = 5.8
flash_time_end = 15.5

write_to_file(f, flash_pe_threshold)
write_to_file(f, flash_time_start)
write_to_file(f, flash_time_end)

write_to_file(f2, flash_pe_threshold)
write_to_file(f2, flash_time_start)
write_to_file(f2, flash_time_end)

# 2D distance between nu-vtx and largest flash centre
tolerance = 80.0

write_to_file(f, tolerance)
write_to_file(f2, tolerance)

# reco track/shower to nu-vtx distance
shwr_nue_tolerance = 4.0
trk_nue_tolerance = 4.0

write_to_file(f, shwr_nue_tolerance)
write_to_file(f, trk_nue_tolerance)

write_to_file(f2, shwr_nue_tolerance)
write_to_file(f2, trk_nue_tolerance)

# leading shower hit thresholds
shwr_hit_threshold = 200.0
shwr_hit_threshold_collection = 80.0

write_to_file(f, shwr_hit_threshold)
write_to_file(f, shwr_hit_threshold_collection)

write_to_file(f2, shwr_hit_threshold)
write_to_file(f2, shwr_hit_threshold_collection)

# leading shower opening angle
tolerance_open_angle_min = 2.0
tolerance_open_angle_max = 15.0

write_to_file(f, tolerance_open_angle_min)
write_to_file(f, tolerance_open_angle_max)

write_to_file(f2, tolerance_open_angle_min)
write_to_file(f2, tolerance_open_angle_max)

# leading shower dE/dx
tolerance_dedx_min = 1.4
tolerance_dedx_max = 3.0

tolerance_photon_dedx_min = 3.5
tolerance_photon_dedx_max = 5.0

write_to_file(f, tolerance_dedx_min)
write_to_file(f, tolerance_dedx_max)

write_to_file(f2, tolerance_photon_dedx_min)
write_to_file(f2, tolerance_photon_dedx_max)

# maximum distance secondary shower from nu-vtx
# 22 cm is roughly two radiation lenghts, which should still keep Pi0 events
dist_tolerance = 22.0

write_to_file(f, dist_tolerance)
write_to_file(f2, dist_tolerance)

# ratio of hits / length (hits / cm)
pfp_hits_length_tolerance = 3.0

write_to_file(f, pfp_hits_length_tolerance)
write_to_file(f2, pfp_hits_length_tolerance)

# ratio of longest track / leading shower lengths
ratio_tolerance = 1.0

write_to_file(f, ratio_tolerance)
write_to_file(f2, ratio_tolerance)

# not listed here as there's no configuration - also require that all
# tracks are contained!

# close input file
f.close()
f2.close()

print "*** *** ***"
print " This script will run the selection twice - once as standard, once for the dE/dx cut on photons "
print "*** *** ***"

# script name + MC name
if(len(sys.argv) == 2):
    call(["../scripts/main.exe", sys.argv[1], "-c",
          "config.txt", "-f", "../scripts/plots/"])
    call(["../scripts/main.exe", sys.argv[1], "-c",
          "photon_config.txt", "-f", "../scripts/bkgplots/"])

# script name + MC name + EXT name
if(len(sys.argv) == 3):
    call(["../scripts/main.exe", sys.argv[1], sys.argv[2],
          "-c", "config.txt", "-f", "../scripts/plots/"])
    call(["../scripts/main.exe", sys.argv[1], sys.argv[2],
          "-c", "photon_config.txt", "-f", "../scripts/bkgplots/"])

# script name + MC name + EXT name
if(len(sys.argv) == 4):
    call(["../scripts/main.exe", sys.argv[1],
          sys.argv[2], sys.argv[3], "-c", "config.txt", "-f", "../scripts/plots/"])
    call(["../scripts/main.exe", sys.argv[1],
          sys.argv[2], sys.argv[3], "-c", "photon_config.txt", "-f", "../scripts/bkgplots/"])
