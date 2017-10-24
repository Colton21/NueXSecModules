#ifndef SELECTION_h
#define SELECTION_h

#include "selection_functions.h"

#include "../xsecAna/LinkDef.h"

namespace xsecSelection {


//const double POT = 4.05982e+19;      //POT - all NuMI + cosmics
const double POT = 2.90469e+21;    //POT - nue + cosmics
const double scaling = 1.52938e-11;  //nues / POT / cm^2
const double genie_xsec = 5.05191e-39; //cm^2

//*******************
// Cut Values
//*******************

//fiducial volume
const double _x1 = 0;
const double _x2 = 0;
const double _y1 = 0;
const double _y2 = 0;
const double _z1 = 0;
const double _z2 = 0;

//in time flash
const int flash_pe_threshold = 50;
const double flash_time_start = 5;
const double flash_time_end = 16;

//vertex to flash
const double tolerance = 100;//cm

//distance between pfp shower and nue object
const double shwr_nue_tolerance = 50;//cm

//hit threshold for showers
const double shwr_hit_threshold = 50;//hits

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
int mc_nu_id = -1;
double mc_nu_vtx_x = -999;
double mc_nu_vtx_y = -999;
double mc_nu_vtx_z = -999;

int run_sum = 0;

int reco_nue_counter = 0;
int reco_nue_counter_nue_cc = 0;
int reco_nue_counter_nue_cc_mixed = 0;
int reco_nue_counter_cosmic = 0;
int reco_nue_counter_nue_nc = 0;
int reco_nue_counter_numu = 0;
int reco_nue_counter_unmatched = 0;
int reco_nue_counter_other_mixed = 0;
int in_fv_counter = 0;
int in_fv_counter_nue_cc = 0;
int in_fv_counter_nue_cc_mixed = 0;
int in_fv_counter_cosmic = 0;
int in_fv_counter_nue_nc = 0;
int in_fv_counter_numu = 0;
int in_fv_counter_unmatched = 0;
int in_fv_counter_other_mixed = 0;
int vtx_flash_counter = 0;
int vtx_flash_counter_nue_cc = 0;
int vtx_flash_counter_nue_cc_mixed = 0;
int vtx_flash_counter_cosmic = 0;
int vtx_flash_counter_nue_nc = 0;
int vtx_flash_counter_numu = 0;
int vtx_flash_counter_unmatched = 0;
int vtx_flash_counter_other_mixed = 0;
int shwr_tpco_counter = 0;
int shwr_tpco_counter_nue_cc = 0;
int shwr_tpco_counter_nue_cc_mixed = 0;
int shwr_tpco_counter_cosmic = 0;
int shwr_tpco_counter_nue_nc = 0;
int shwr_tpco_counter_numu = 0;
int shwr_tpco_counter_unmatched = 0;
int shwr_tpco_counter_other_mixed = 0;
int hit_threshold_counter = 0;
int hit_threshold_counter_nue_cc = 0;
int hit_threshold_counter_nue_cc_mixed = 0;
int hit_threshold_counter_cosmic = 0;
int hit_threshold_counter_nue_nc = 0;
int hit_threshold_counter_numu = 0;
int hit_threshold_counter_unmatched = 0;
int hit_threshold_counter_other_mixed = 0;
std::vector<int> tabulated_origins;


}//end namespace

#endif
