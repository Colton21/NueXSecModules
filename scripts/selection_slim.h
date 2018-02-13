#ifndef SELECTION_h
#define SELECTION_h

#include "selection_functions.h"
#include "selection_cuts.h"

#include "../xsecAna/LinkDef.h"

namespace xsecSelection {


//const double POT = 4.05982e+19;      //POT - all NuMI + cosmics
const double POT = 1.23206e+20; //POT - all NuMI + cosmics, bigger sample
//const double POT = 2.90469e+21;    //POT - nue + cosmics
const double scaling = 1.52938e-11;  //nues / POT / cm^2
const double genie_xsec = 5.05191e-39; //cm^2

//*******************
// Cut Values
//*******************

//fiducial volume
const double _x1 = 15;
const double _x2 = 15;
const double _y1 = 30;
const double _y2 = 15;
const double _z1 = 15;
const double _z2 = 15;

// const double _x1 = 0;
// const double _x2 = 0;
// const double _y1 = 0;
// const double _y2 = 0;
// const double _z1 = 0;
// const double _z2 = 0;

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

//hit threshold for showers
//standard 50 hits
const double shwr_hit_threshold = 200;//hits

//hit threshold for at least one shower on collection plane
//standard 50 hits
const double shwr_hit_threshold_collection = 80;//hits

//tolerance for leading shower open angle
//standard 20 degrees
//Roberto uses 15 degrees
const std::vector<double> tolerance_open_angle {2, 15};//degrees

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

//********************
//********************

int total_mc_entries_inFV =0;

int mc_nue_cc_counter = 0;
int mc_nue_nc_counter = 0;
int mc_numu_cc_counter = 0;
int mc_numu_nc_counter = 0;
int mc_nue_cc_counter_bar = 0;
int mc_numu_cc_counter_bar = 0;
int mc_nue_nc_counter_bar = 0;
int mc_numu_nc_counter_bar = 0;
double mc_nu_energy = 0;
double mc_nu_momentum = 0;
int mc_nu_id = -1;
double mc_nu_vtx_x = -999;
double mc_nu_vtx_y = -999;
double mc_nu_vtx_z = -999;

double mc_nu_dir_x = -999;
double mc_nu_dir_y = -999;
double mc_nu_dir_z = -999;
double mc_ele_dir_x = -999;
double mc_ele_dir_y = -999;
double mc_ele_dir_z = -999;
double mc_ele_energy = 0;
double mc_ele_momentum = 0;
bool has_pi0 = false;
double mc_nu_time = -1;

int mc_nu_num_particles = 0;
int mc_nu_num_charged_particles = 0;

int run_sum = 0;
int out_of_time_sum = 0;
int low_pe_sum = 0;

std::vector<int> tabulated_origins;

}//end namespace

#endif
