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
# 22 cm is roughly two radiation lenghts, which should still keep Pi0 events
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
    call(["../scripts/main.exe", sys.argv[1], "-c", "config.txt"])

# script name + MC name + EXT name
if(len(sys.argv) == 3):
    call(["../scripts/main.exe", sys.argv[1], sys.argv[2], "-c", "config.txt"])

# script name + MC name + EXT name
if(len(sys.argv) == 4):
    call(["../scripts/main.exe", sys.argv[1],
          sys.argv[2], sys.argv[3], "-c", "config.txt"])
