#ifndef MAIN_h
#define MAIN_h

#include "selection.h"
#include "utility.h"

//*******************
// Cut Values
//*******************

//fiducial volume
// const double _x1 = 0;
// const double _x2 = 0;
// const double _y1 = 0;
// const double _y2 = 0;
// const double _z1 = 0;
// const double _z2 = 0;

const double _x1 = 40;
const double _x2 = 40;
const double _y1 = 40;
const double _y2 = 40;
const double _z1 = 40;
const double _z2 = 40;

//in time flash
const int flash_pe_threshold = 50;
const double flash_time_start = 5;
const double flash_time_end = 16;

//vertex to flash
//standard 100 cm
const double tolerance = 80;//cm

//distance between pfp shower and nue object
//standard 50 cm
//Roberto runs with 4cm
const double shwr_nue_tolerance = 4;//cm
const double trk_nue_tolerance = 4;//cm

//hit threshold for at least one shower
//standard 120 hits
//higher cuts 200 hits
const double shwr_hit_threshold = 200;//hits

//hit threshold for at least one shower on collection plane
//standard 50 hits
const double shwr_hit_threshold_collection = 80;//hits

//tolerance for leading shower open angle
//standard 20 degrees
//Roberto uses 15 degrees
//const std::vector<double> tolerance_open_angle {2, 15};//degrees
const double tolerance_open_angle_min = 2;
const double tolerance_open_angle_max = 15;

//tolerance for dedx of leading shower
//Roberto uses: 1.4 - 3 MeV / cm
//standard is: 0, 3.5
const double tolerance_dedx_min = 1.4;
const double tolerance_dedx_max = 3;

//tolerance for distance from the reco nue vtx for TPCO w/ >3 showers
const double dist_tolerance = 22; //cm
//22 cm is something like ~2 radiation lengths - I expect that TPCO w/ >3 showers
//should have most of the activity around the nucleus right? or at least within the 22 cm

//tolerance for hits/length - these should be a property of a shower if it's true
const double pfp_hits_length_tolerance = 3; //hits/cm

//tolerance for longest track length / leading shower length
const double ratio_tolerance = 1;

#endif
