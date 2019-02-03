#ifndef MAIN_h
#define MAIN_h

#include "selection.h"
#include "selection_slim.h"
#include "variation_output.h"
#include "utility.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

/*

   This is the main header for the seleciton functions

   details about the ordering of the selection cuts can be most easily found in
   selection_slim.cxx.

   The slim and non-slim selection simply have the difference of creating plots (histograms)

   The main controls mostly the file i/o and handling inputs from the python scrips
   which can help with conveniently running the selection.

   To investigate exactly how each cut functions, look at selection_cuts.h/.cxx

   selection_functions / selection_functions_data operate primarily as classes to hold
   functions used in filling histograms

   You can always run the selection "manually" from the main.exe created in the makefile,
   just be careful of formatting the arguments (look in main.cxx). The scripts in
   ../python/ handle some of the arguments for you.

 */

//this is the default place
char * file_locate_prefix = (char*)"../scripts/plots/";
std::string file_locate_prefix_2;

//These values are the defaults and will run if no parameter list is given
//The parameter list is generated when running any of the python scripts

//*******************
// Cut Values
//*******************

//fiducial volume
const double _x1 = 20;
const double _x2 = 20;
const double _y1 = 20;
const double _y2 = 20;
const double _z1 = 20;
const double _z2 = 20;

// const double _x1 = 40;
// const double _x2 = 40;
// const double _y1 = 40;
// const double _y2 = 40;
// const double _z1 = 40;
// const double _z2 = 40;

//in time flash
const int flash_pe_threshold = 50;
//const double flash_time_start = 5.8;
//const double flash_time_end = 15.5;
const double flash_time_start = 5.5;
const double flash_time_end = 16.0;

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
//const double pfp_hits_length_tolerance = 5; //hits/cm

//tolerance for longest track length / leading shower length
const double ratio_tolerance = 1;

const bool detector_variations = false;

#endif
