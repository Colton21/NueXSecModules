#ifndef SELECTION_h
#define SELECTION_h

#include "selection_functions.h"

#include "../xsecAna/LinkDef.h"

namespace xsecSelection {


const double POT = 4.05982e+19;      //POT - all NuMI + cosmics
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
const double shwr_hit_threshold = 100;//hits

//tolerance for leading shower open angle
//standard 20 degrees
//Roberto uses 15 degrees
const double tolerance_open_angle = 15;//degrees

//tolerance for dedx of leading shower
//Roberto uses: 1.4 - 3 MeV / cm
//standard is:
const double tolerance_dedx_min = 0;
const double tolerance_dedx_max = 10000;

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

int mc_nu_num_particles = 0;
int mc_nu_num_charged_particles = 0;

int run_sum = 0;

int reco_nue_counter = 0;
int reco_nue_counter_nue_cc = 0;
int reco_nue_counter_nue_cc_mixed = 0;
int reco_nue_counter_nue_cc_out_fv = 0;
int reco_nue_counter_cosmic = 0;
int reco_nue_counter_nue_nc = 0;
int reco_nue_counter_numu_cc = 0;
int reco_nue_counter_numu_cc_mixed = 0;
int reco_nue_counter_numu_nc = 0;
int reco_nue_counter_unmatched = 0;
int reco_nue_counter_other_mixed = 0;
int in_fv_counter = 0;
int in_fv_counter_nue_cc = 0;
int in_fv_counter_nue_cc_mixed = 0;
int in_fv_counter_nue_cc_out_fv = 0;
int in_fv_counter_cosmic = 0;
int in_fv_counter_nue_nc = 0;
int in_fv_counter_numu_cc = 0;
int in_fv_counter_numu_cc_mixed = 0;
int in_fv_counter_numu_nc = 0;
int in_fv_counter_unmatched = 0;
int in_fv_counter_other_mixed = 0;
int vtx_flash_counter = 0;
int vtx_flash_counter_nue_cc = 0;
int vtx_flash_counter_nue_cc_mixed = 0;
int vtx_flash_counter_nue_cc_out_fv = 0;
int vtx_flash_counter_cosmic = 0;
int vtx_flash_counter_nue_nc = 0;
int vtx_flash_counter_numu_cc = 0;
int vtx_flash_counter_numu_cc_mixed = 0;
int vtx_flash_counter_numu_nc = 0;
int vtx_flash_counter_unmatched = 0;
int vtx_flash_counter_other_mixed = 0;
int shwr_tpco_counter = 0;
int shwr_tpco_counter_nue_cc = 0;
int shwr_tpco_counter_nue_cc_mixed = 0;
int shwr_tpco_counter_nue_cc_out_fv = 0;
int shwr_tpco_counter_cosmic = 0;
int shwr_tpco_counter_nue_nc = 0;
int shwr_tpco_counter_numu_cc = 0;
int shwr_tpco_counter_numu_cc_mixed = 0;
int shwr_tpco_counter_numu_nc = 0;
int shwr_tpco_counter_unmatched = 0;
int shwr_tpco_counter_other_mixed = 0;
int trk_tpco_counter = 0;
int trk_tpco_counter_nue_cc = 0;
int trk_tpco_counter_nue_cc_mixed = 0;
int trk_tpco_counter_nue_cc_out_fv = 0;
int trk_tpco_counter_cosmic = 0;
int trk_tpco_counter_nue_nc = 0;
int trk_tpco_counter_numu_cc = 0;
int trk_tpco_counter_numu_cc_mixed = 0;
int trk_tpco_counter_numu_nc = 0;
int trk_tpco_counter_unmatched = 0;
int trk_tpco_counter_other_mixed = 0;
int hit_threshold_counter = 0;
int hit_threshold_counter_nue_cc = 0;
int hit_threshold_counter_nue_cc_mixed = 0;
int hit_threshold_counter_nue_cc_out_fv = 0;
int hit_threshold_counter_cosmic = 0;
int hit_threshold_counter_nue_nc = 0;
int hit_threshold_counter_numu_cc = 0;
int hit_threshold_counter_numu_cc_mixed = 0;
int hit_threshold_counter_numu_nc = 0;
int hit_threshold_counter_unmatched = 0;
int hit_threshold_counter_other_mixed = 0;
int open_angle_counter = 0;
int open_angle_counter_nue_cc = 0;
int open_angle_counter_nue_cc_mixed = 0;
int open_angle_counter_nue_cc_out_fv = 0;
int open_angle_counter_cosmic = 0;
int open_angle_counter_nue_nc = 0;
int open_angle_counter_numu_cc = 0;
int open_angle_counter_numu_cc_mixed = 0;
int open_angle_counter_numu_nc = 0;
int open_angle_counter_unmatched = 0;
int open_angle_counter_other_mixed = 0;
std::vector<int> tabulated_origins;


}//end namespace

#endif
