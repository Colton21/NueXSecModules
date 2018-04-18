#ifndef SELECTION_h
#define SELECTION_h

#include "selection_functions.h"
#include "selection_cuts.h"
#include "selection_functions_data.h"
#include "histogram_functions.h"

#include "../xsecAna/LinkDef.h"

#include <fstream>

namespace xsecSelection {


const double POT = 1.82027e21;

//const double POT = 4.05982e+19;      //POT - all NuMI + cosmics
//const double POT = 1.23206e+20; //POT - all NuMI + cosmics, bigger sample
//const double POT = 2.90469e+21;    //POT - nue + cosmics
const double scaling_nue = 1.52938e-11;  //nues / POT / cm^2
const double scaling_nue_bar = 7.77111e-12; //anues / POT / cm^2
const double scaling = (scaling_nue + scaling_nue_bar);

//since I have both nue and nue_bar as signal definition need to adjust for this
const double genie_xsec_nue = 5.63067e-39; //cm^2
const double genie_xsec_nue_bar = 2.0893e-39; //cm^2
const double genie_xsec = genie_xsec_nue + genie_xsec_nue_bar;

/*
   3e13 POT / spills for NuMI -> 6.0675667e7 triggers for MC, or 5.9966667e5 triggers when MC scaled to data POT
   2.571102 Million EXT spills
   scale EXT by: 0.23323333 when MC scaled to data
   otherwise scale EXT by: intime_scale_factor / data_scale_factor (23.599090)

   For new EXT: EXT 352180
   But looking just at samdef = 2571102 ??? why ???
   intime_scale_factor = 0.59966667e6 / 352180
   intime_scale_factor =
 */
//const double intime_scale_factor = 23.5991 * (0.0098831492);
//const double intime_scale_factor = 23.5991 * (0.0364451);
//const double intime_scale_factor = 23.5991 * (0.0556456);

/*
   scale via POT - 1.23e20 MC / 2.189e19 POT for full april sample
   upscale Data by: 5.62
   but if MC is downscaled by EXT difference, then we do 5.62 / 1.5946470
   = 3.52
 */

//new data set has : 1.799e19 POT
//EA9CNT: 472210
//const double data_scale_factor = 0.0098831492;

//tor101_wcut for whole samdef = 6.634e+19 POT
//EA9CNT_wcut = 1678059
const double intime_scale_factor = 0.655210;
const double data_scale_factor = 0.0364451;

//tor101_wcut for summing 1,2,3,4 samdef manually = 1.0129e20
//const double data_scale_factor = 0.0556456;


const double theta_translation = 29.36 * (3.1415/180);
const double phi_translation = 8.121 * (3.1415/180);
// const double theta_translation = 0.0;
// const double phi_translation = 0.0;

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

//TEfficiency histograms
TH1D * h_nue_eng_eff_den           = new TH1D("h_nue_eng_eff_den", "h_nue_eng_eff_den", 6, 0, 4);
TH1D * h_nue_eng_eff_num           = new TH1D("h_nue_eng_eff_num", "h_nue_eng_eff_num", 6, 0, 4);
TH1D * h_ele_eng_eff_den           = new TH1D("h_ele_eng_eff_den", "h_nue_ele_eff_den", 6, 0, 4);
TH1D * h_ele_eng_eff_num           = new TH1D("h_ele_eng_eff_num", "h_nue_ele_eff_num", 6, 0, 4);
TH1D * h_ele_eng_eff_num_pre_cuts  = new TH1D("h_ele_eng_eff_num_pre_cuts", "h_ele_eng_eff_num_pre_cuts", 6, 0, 4);
TH1D * h_nue_vtx_x_eff_den         = new TH1D("h_nue_vtx_x_eff_den", "h_nue_vtx_x_eff_den", 6, 0, 256.35);
TH1D * h_nue_vtx_x_eff_num         = new TH1D("h_nue_vtx_x_eff_num", "h_nue_vtx_x_eff_num", 6, 0, 256.35);
TH1D * h_nue_vtx_y_eff_den         = new TH1D("h_nue_vtx_y_eff_den", "h_nue_vtx_y_eff_den", 6, -116.5, 116.5);
TH1D * h_nue_vtx_y_eff_num         = new TH1D("h_nue_vtx_y_eff_num", "h_nue_vtx_y_eff_num", 6, -116.5, 116.5);
TH1D * h_nue_vtx_z_eff_den         = new TH1D("h_nue_vtx_z_eff_den", "h_nue_vtx_z_eff_den", 6, 0, 1036.8);
TH1D * h_nue_vtx_z_eff_num         = new TH1D("h_nue_vtx_z_eff_num", "h_nue_vtx_z_eff_num", 6, 0, 1036.8);
TH1D * h_nue_dir_x_eff_den         = new TH1D("h_nue_dir_x_eff_den", "h_nue_dir_x_eff_den", 10, -1, 1);
TH1D * h_nue_dir_x_eff_num         = new TH1D("h_nue_dir_x_eff_num", "h_nue_dir_x_eff_num", 10, -1, 1);
TH1D * h_nue_dir_y_eff_den         = new TH1D("h_nue_dir_y_eff_den", "h_nue_dir_y_eff_den", 10, -1, 1);
TH1D * h_nue_dir_y_eff_num         = new TH1D("h_nue_dir_y_eff_num", "h_nue_dir_y_eff_num", 10, -1, 1);
TH1D * h_nue_dir_z_eff_den         = new TH1D("h_nue_dir_z_eff_den", "h_nue_dir_z_eff_den", 10, -1, 1);
TH1D * h_nue_dir_z_eff_num         = new TH1D("h_nue_dir_z_eff_num", "h_nue_dir_z_eff_num", 10, -1, 1);
TH1D * h_ele_dir_x_eff_den         = new TH1D("h_ele_dir_x_eff_den", "h_ele_dir_x_eff_den", 10, -1, 1);
TH1D * h_ele_dir_x_eff_num         = new TH1D("h_ele_dir_x_eff_num", "h_ele_dir_x_eff_num", 10, -1, 1);
TH1D * h_ele_dir_y_eff_den         = new TH1D("h_ele_dir_y_eff_den", "h_ele_dir_y_eff_den", 10, -1, 1);
TH1D * h_ele_dir_y_eff_num         = new TH1D("h_ele_dir_y_eff_num", "h_ele_dir_y_eff_num", 10, -1, 1);
TH1D * h_ele_dir_z_eff_den         = new TH1D("h_ele_dir_z_eff_den", "h_ele_dir_z_eff_den", 10, -1, 1);
TH1D * h_ele_dir_z_eff_num         = new TH1D("h_ele_dir_z_eff_num", "h_ele_dir_z_eff_num", 10, -1, 1);
TH1D * h_nue_num_part_eff_den      = new TH1D("h_nue_num_part_eff_den", "h_nue_num_part_eff_den", 10, 0, 20);
TH1D * h_nue_num_part_eff_num      = new TH1D("h_nue_num_part_eff_num", "h_nue_num_part_eff_num", 10, 0, 20);
TH1D * h_nue_num_chrg_part_eff_den = new TH1D("h_nue_num_chrg_part_eff_den", "h_nue_num_chrg_part_eff_den", 10, 0, 20);
TH1D * h_nue_num_chrg_part_eff_num = new TH1D("h_nue_num_chrg_part_eff_num", "h_nue_num_chrg_part_eff_num", 10, 0, 20);
TH1D * h_nue_cos_theta_eff_den     = new TH1D("h_nue_cos_theta_eff_den", "h_nue_cos_theta_eff_den", 6, -1, 1);
TH1D * h_nue_cos_theta_eff_num     = new TH1D("h_nue_cos_theta_eff_num", "h_nue_cos_theta_eff_num", 6, -1, 1);
TH1D * h_nue_phi_eff_den           = new TH1D("h_nue_phi_eff_den", "h_nue_phi_eff_den", 10, -180, 180);
TH1D * h_nue_phi_eff_num           = new TH1D("h_nue_phi_eff_num", "h_nue_phi_eff_num", 10, -180, 180);
TH1D * h_ele_cos_theta_eff_den     = new TH1D("h_ele_cos_theta_eff_den", "h_ele_cos_theta_eff_den", 6, -1, 1);
TH1D * h_ele_cos_theta_eff_num     = new TH1D("h_ele_cos_theta_eff_num", "h_ele_cos_theta_eff_num", 6, -1, 1);
TH1D * h_ele_cos_theta_eff_num_pre_cuts = new TH1D ("h_ele_cos_theta_eff_num_pre_cuts", "h_ele_cos_theta_eff_num_pre_cuts", 6, -1, 1);
TH1D * h_ele_phi_eff_den           = new TH1D("h_ele_phi_eff_den", "h_ele_phi_eff_den", 10, -180, 180);
TH1D * h_ele_phi_eff_num           = new TH1D("h_ele_phi_eff_num", "h_ele_phi_eff_num", 10, -180, 180);
TH1D * h_ele_theta_eff_num         = new TH1D("h_ele_theta_eff_num", "h_ele_theta_eff_num", 10, 0, 180);
TH1D * h_ele_theta_eff_den         = new TH1D("h_ele_theta_eff_den", "h_ele_theta_eff_den", 10, 0, 180);
//

TH1D * h_flash_time        = new TH1D ("h_flash_time",        "h_flash_time",        20, 0, 20);
TH1D * h_flash_time_intime = new TH1D ("h_flash_time_intime", "h_flash_time_intime", 20, 0, 20);
TH1D * h_flash_time_data   = new TH1D ("h_flash_time_data",   "h_flash_time_data",   20, 0, 20);

//
TH2I * h_tracks_showers         = new TH2I("h_tracks_showers", "h_tracks_showers", 8, 0, 8, 8, 0, 8);
TH2I * h_tracks_showers_cosmic  = new TH2I("h_tracks_showers_cosmic", "h_tracks_showers_cosmic", 8, 0, 8, 8, 0, 8);
TH2I * h_tracks_showers_numu    = new TH2I("h_tracks_showers_numu", "h_tracks_showers_numu", 8, 0, 8, 8, 0, 8);

TH1D * h_leading_shower_open_angle_nue_cc        = new TH1D("h_leading_shower_open_angle_nue_cc",        "h_leading_shower_open_angle_nue_cc",        25, 0, 50);
TH1D * h_leading_shower_open_angle_nue_cc_mixed  = new TH1D("h_leading_shower_open_angle_nue_cc_mixed",  "h_leading_shower_open_angle_nue_cc_mixed",  25, 0, 50);
TH1D * h_leading_shower_open_angle_nue_cc_out_fv = new TH1D("h_leading_shower_open_angle_nue_cc_out_fv", "h_leading_shower_open_angle_nue_cc_out_fv", 25, 0, 50);
TH1D * h_leading_shower_open_angle_numu_cc       = new TH1D("h_leading_shower_open_angle_numu_cc",       "h_leading_shower_open_angle_numu_cc",       25, 0, 50);
TH1D * h_leading_shower_open_angle_nc            = new TH1D("h_leading_shower_open_angle_nc",            "h_leading_shower_open_angle_nc",            25, 0, 50);
TH1D * h_leading_shower_open_angle_cosmic        = new TH1D("h_leading_shower_open_angle_cosmic",        "h_leading_shower_open_angle_cosmic",        25, 0, 50);
TH1D * h_leading_shower_open_angle_nc_pi0        = new TH1D("h_leading_shower_open_angle_nc_pi0",        "h_leading_shower_open_angle_nc_pi0",        25, 0, 50);
TH1D * h_leading_shower_open_angle_numu_cc_mixed = new TH1D("h_leading_shower_open_angle_numu_cc_mixed", "h_leading_shower_open_angle_numu_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_other_mixed   = new TH1D("h_leading_shower_open_angle_other_mixed",   "h_leading_shower_open_angle_other_mixed",   25, 0, 50);
TH1D * h_leading_shower_open_angle_unmatched     = new TH1D("h_leading_shower_open_angle_unmatched",     "h_leading_shower_open_angle_unmatched",     25, 0, 50);
TH1D * h_leading_shower_open_angle_intime        = new TH1D("h_leading_shower_open_angle_intime",        "h_leading_shower_open_angle_intime",        25, 0, 50);
TH1D * h_leading_shower_open_angle_data          = new TH1D("h_leading_shower_open_angle_data",          "h_leading_shower_open_angle_data",          25, 0, 50);

TH1D * h_leading_shower_open_angle_nue_cc_after        = new TH1D("h_leading_shower_open_angle_nue_cc_after",        "h_leading_shower_open_angle_nue_cc_after",        25, 0, 50);
TH1D * h_leading_shower_open_angle_nue_cc_mixed_after  = new TH1D("h_leading_shower_open_angle_nue_cc_mixed_after",  "h_leading_shower_open_angle_nue_cc_mixed_after",  25, 0, 50);
TH1D * h_leading_shower_open_angle_nue_cc_out_fv_after = new TH1D("h_leading_shower_open_angle_nue_cc_out_fv_after", "h_leading_shower_open_angle_nue_cc_out_fv_after", 25, 0, 50);
TH1D * h_leading_shower_open_angle_numu_cc_after       = new TH1D("h_leading_shower_open_angle_numu_cc_after",       "h_leading_shower_open_angle_numu_cc_after",       25, 0, 50);
TH1D * h_leading_shower_open_angle_nc_after            = new TH1D("h_leading_shower_open_angle_nc_after",            "h_leading_shower_open_angle_nc_after",            25, 0, 50);
TH1D * h_leading_shower_open_angle_cosmic_after        = new TH1D("h_leading_shower_open_angle_cosmic_after",        "h_leading_shower_open_angle_cosmic_after",        25, 0, 50);
TH1D * h_leading_shower_open_angle_nc_pi0_after        = new TH1D("h_leading_shower_open_angle_nc_pi0_after",        "h_leading_shower_open_angle_nc_pi0_after",        25, 0, 50);
TH1D * h_leading_shower_open_angle_numu_cc_mixed_after = new TH1D("h_leading_shower_open_angle_numu_cc_mixed_after", "h_leading_shower_open_angle_numu_cc_mixed_after", 25, 0, 50);
TH1D * h_leading_shower_open_angle_other_mixed_after   = new TH1D("h_leading_shower_open_angle_other_mixed_after",   "h_leading_shower_open_angle_other_mixed_after",   25, 0, 50);
TH1D * h_leading_shower_open_angle_unmatched_after     = new TH1D("h_leading_shower_open_angle_unmatched_after",     "h_leading_shower_open_angle_unmatched_after",     25, 0, 50);
TH1D * h_leading_shower_open_angle_intime_after        = new TH1D("h_leading_shower_open_angle_intime_after",        "h_leading_shower_open_angle_intime_after",        25, 0, 50);
TH1D * h_leading_shower_open_angle_data_after          = new TH1D("h_leading_shower_open_angle_data_after",          "h_leading_shower_open_angle_data_after",          25, 0, 50);

TH1D * h_leading_shower_open_angle_1_nue_cc        = new TH1D("h_leading_shower_open_angle_1_nue_cc",        "h_leading_shower_open_angle_1_nue_cc",        25, 0, 50);
TH1D * h_leading_shower_open_angle_1_nue_cc_mixed  = new TH1D("h_leading_shower_open_angle_1_nue_cc_mixed",  "h_leading_shower_open_angle_1_nue_cc_mixed",  25, 0, 50);
TH1D * h_leading_shower_open_angle_1_nue_cc_out_fv = new TH1D("h_leading_shower_open_angle_1_nue_cc_out_fv", "h_leading_shower_oepn_angle_1_nue_cc_out_fv", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_numu_cc       = new TH1D("h_leading_shower_open_angle_1_numu_cc",       "h_leading_shower_open_angle_1_numu_cc",       25, 0, 50);
TH1D * h_leading_shower_open_angle_1_nc            = new TH1D("h_leading_shower_open_angle_1_nc",            "h_leading_shower_open_angle_1_nc",            25, 0, 50);
TH1D * h_leading_shower_open_angle_1_cosmic        = new TH1D("h_leading_shower_open_angle_1_cosmic",        "h_leading_shower_open_angle_1_cosmic",        25, 0, 50);
TH1D * h_leading_shower_open_angle_1_nc_pi0        = new TH1D("h_leading_shower_open_angle_1_nc_pi0",        "h_leading_shower_open_angle_1_nc_pi0",        25, 0, 50);
TH1D * h_leading_shower_open_angle_1_numu_cc_mixed = new TH1D("h_leading_shower_open_angle_1_numu_cc_mixed", "h_leading_shower_open_angle_1_numu_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_other_mixed   = new TH1D("h_leading_shower_open_angle_1_other_mixed",   "h_leading_shower_open_angle_1_other_mixed",   25, 0, 50);
TH1D * h_leading_shower_open_angle_1_unmatched     = new TH1D("h_leading_shower_open_angle_1_unmatched",     "h_leading_shower_open_angle_1_unmatched",     25, 0, 50);
TH1D * h_leading_shower_open_angle_1_intime        = new TH1D("h_leading_shower_open_angle_1_intime",        "h_leading_shower_open_angle_1_intime",        25, 0, 50);
TH1D * h_leading_shower_open_angle_1_data          = new TH1D("h_leading_shower_open_angle_1_data",          "h_leading_shower_open_angle_1_data",          25, 0, 50);

TH1D * h_leading_shower_open_angle_2plus_nue_cc        = new TH1D("h_leading_shower_open_angle_2plus_nue_cc",        "h_leading_shower_open_angle_2plus_nue_cc",        25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_nue_cc_mixed  = new TH1D("h_leading_shower_open_angle_2plus_nue_cc_mixed",  "h_leading_shower_open_angle_2plus_nue_cc_mixed",  25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_nue_cc_out_fv = new TH1D("h_leading_shower_open_angle_2plus_nue_cc_out_fv", "h_leading_shower_open_angle_2plus_nue_cc_out_fv", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_numu_cc       = new TH1D("h_leading_shower_open_angle_2plus_numu_cc",       "h_leading_shower_open_angle_2plus_numu_cc",       25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_nc            = new TH1D("h_leading_shower_open_angle_2plus_nc",            "h_leading_shower_open_angle_2plus_nc",            25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_cosmic        = new TH1D("h_leading_shower_open_angle_2plus_cosmic",        "h_leading_shower_open_angle_2plus_cosmic",        25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_nc_pi0        = new TH1D("h_leading_shower_open_angle_2plus_nc_pi0",        "h_leading_shower_open_angle_2plus_nc_pi0",        25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_numu_cc_mixed = new TH1D("h_leading_shower_open_angle_2plus_numu_cc_mixed", "h_leading_shower_open_angle_2plus_numu_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_other_mixed   = new TH1D("h_leading_shower_open_angle_2plus_other_mixed",   "h_leading_shower_open_angle_2plus_other_mixed",   25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_unmatched     = new TH1D("h_leading_shower_open_angle_2plus_unmatched",     "h_leading_shower_open_angle_2plus_unmatched",     25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_intime        = new TH1D("h_leading_shower_open_angle_2plus_intime",        "h_leading_shower_open_angle_2plus_intime",        25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_data          = new TH1D("h_leading_shower_open_angle_2plus_data",          "h_leading_shower_open_angle_2plus_data",          25, 0, 50);

TH1D * h_trk_vtx_dist_nue_cc        = new TH1D("h_trk_vtx_dist_nue_cc",        "h_trk_vtx_dist_nue_cc",        25, 0, 20);
TH1D * h_trk_vtx_dist_nue_cc_mixed  = new TH1D("h_trk_vtx_dist_nue_cc_mixed",  "h_trk_vtx_dist_nue_cc_mixed",  25, 0, 20);
TH1D * h_trk_vtx_dist_nue_cc_out_fv = new TH1D("h_trk_vtx_dist_nue_cc_out_fv", "h_trk_vtx_dist_nue_cc_out_fv", 25, 0, 20);
TH1D * h_trk_vtx_dist_numu_cc       = new TH1D("h_trk_vtx_dist_numu_cc",       "h_trk_vtx_dist_numu_cc",       25, 0, 20);
TH1D * h_trk_vtx_dist_nc            = new TH1D("h_trk_vtx_dist_nc",            "h_trk_vtx_dist_nc",            25, 0, 20);
TH1D * h_trk_vtx_dist_cosmic        = new TH1D("h_trk_vtx_dist_cosmic",        "h_trk_vtx_dist_cosmic",        25, 0, 20);
TH1D * h_trk_vtx_dist_nc_pi0        = new TH1D("h_trk_vtx_dist_nc_pi0",        "h_trk_vtx_dist_nc_pi0",        25, 0, 20);
TH1D * h_trk_vtx_dist_numu_cc_mixed = new TH1D("h_trk_vtx_dist_numu_cc_mixed", "h_trk_vtx_dist_numu_cc_mixed", 25, 0, 20);
TH1D * h_trk_vtx_dist_other_mixed   = new TH1D("h_trk_vtx_dist_other_mixed",   "h_trk_vtx_dist_other_mixed",   25, 0, 20);
TH1D * h_trk_vtx_dist_unmatched     = new TH1D("h_trk_vtx_dist_unmatched",     "h_trk_vtx_dist_unmatched",     25, 0, 20);
TH1D * h_trk_vtx_dist_intime        = new TH1D("h_trk_vtx_dist_intime",        "h_trk_vtx_dist_intime",        25, 0, 20);
TH1D * h_trk_vtx_dist_data          = new TH1D("h_trk_vtx_dist_data",          "h_trk_vtx_dist_data",          25, 0, 20);

TH1D * h_trk_vtx_dist_nue_cc_after        = new TH1D("h_trk_vtx_dist_nue_cc_after",        "h_trk_vtx_dist_nue_cc_after",        25, 0, 20);
TH1D * h_trk_vtx_dist_nue_cc_mixed_after  = new TH1D("h_trk_vtx_dist_nue_cc_mixed_after",  "h_trk_vtx_dist_nue_cc_mixed_after",  25, 0, 20);
TH1D * h_trk_vtx_dist_nue_cc_out_fv_after = new TH1D("h_trk_vtx_dist_nue_cc_out_fv_after", "h_trk_vtx_dist_nue_cc_out_fv_after", 25, 0, 20);
TH1D * h_trk_vtx_dist_numu_cc_after       = new TH1D("h_trk_vtx_dist_numu_cc_after",       "h_trk_vtx_dist_numu_cc_after",       25, 0, 20);
TH1D * h_trk_vtx_dist_nc_after            = new TH1D("h_trk_vtx_dist_nc_after",            "h_trk_vtx_dist_nc_after",            25, 0, 20);
TH1D * h_trk_vtx_dist_cosmic_after        = new TH1D("h_trk_vtx_dist_cosmic_after",        "h_trk_vtx_dist_cosmic_after",        25, 0, 20);
TH1D * h_trk_vtx_dist_nc_pi0_after        = new TH1D("h_trk_vtx_dist_nc_pi0_after",        "h_trk_vtx_dist_nc_pi0_after",        25, 0, 20);
TH1D * h_trk_vtx_dist_numu_cc_mixed_after = new TH1D("h_trk_vtx_dist_numu_cc_mixed_after", "h_trk_vtx_dist_numu_cc_mixed_after", 25, 0, 20);
TH1D * h_trk_vtx_dist_other_mixed_after   = new TH1D("h_trk_vtx_dist_other_mixed_after",   "h_trk_vtx_dist_other_mixed_after",   25, 0, 20);
TH1D * h_trk_vtx_dist_unmatched_after     = new TH1D("h_trk_vtx_dist_unmatched_after",     "h_trk_vtx_dist_unmatched_after",     25, 0, 20);
TH1D * h_trk_vtx_dist_intime_after        = new TH1D("h_trk_vtx_dist_intime_after",        "h_trk_vtx_dist_intime_after",        25, 0, 20);
TH1D * h_trk_vtx_dist_data_after          = new TH1D("h_trk_vtx_dist_data_after",          "h_trk_vtx_dist_data_after",          25, 0, 20);

TH2D * h_pfp_track_shower_nue_cc_qe     = new TH2D("h_pfp_track_shower_nue_cc_qe",     "h_pfp_track_shower_nue_cc_qe", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_out_fv = new TH2D("h_pfp_track_shower_nue_cc_out_fv", "h_pfp_track_shower_nue_cc_out_fv", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_res    = new TH2D("h_pfp_track_shower_nue_cc_res",    "h_pfp_track_shower_nue_cc_res", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_dis    = new TH2D("h_pfp_track_shower_nue_cc_dis",    "h_pfp_track_shower_nue_cc_dis", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_coh    = new TH2D("h_pfp_track_shower_nue_cc_coh",    "h_pfp_track_shower_nue_cc_coh", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_mec    = new TH2D("h_pfp_track_shower_nue_cc_mec",    "h_pfp_track_shower_nue_cc_mec", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nc            = new TH2D("h_pfp_track_shower_nc",            "h_pfp_track_shower_nc", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_qe    = new TH2D("h_pfp_track_shower_numu_cc_qe",    "h_pfp_track_shower_numu_cc_qe", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_res   = new TH2D("h_pfp_track_shower_numu_cc_res",   "h_pfp_track_shower_numu_cc_res", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_dis   = new TH2D("h_pfp_track_shower_numu_cc_dis",   "h_pfp_track_shower_numu_cc_dis", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_coh   = new TH2D("h_pfp_track_shower_numu_cc_coh",   "h_pfp_track_shower_numu_cc_coh", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_mec   = new TH2D("h_pfp_track_shower_numu_cc_mec",   "h_pfp_track_shower_numu_cc_mec", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nc_pi0        = new TH2D("h_pfp_track_shower_nc_pi0",        "h_pfp_track_shower_nc_pi0", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_mixed  = new TH2D("h_pfp_track_shower_nue_cc_mixed",  "h_pfp_track_shower_nue_cc_mixed", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_mixed = new TH2D("h_pfp_track_shower_numu_cc_mixed", "h_pfp_track_shower_numu_cc_mixed", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_cosmic        = new TH2D("h_pfp_track_shower_cosmic",        "h_pfp_track_shower_cosmic", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_other_mixed   = new TH2D("h_pfp_track_shower_other_mixed",   "h_pfp_track_shower_other_mixed", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_unmatched     = new TH2D("h_pfp_track_shower_unmatched",     "h_pfp_track_shower_unmatched", 10, 0, 10, 10, 0, 10);

TH2D * h_leading_shower_mc_pdg_nue_cc_qe     = new TH2D("h_leading_shower_mc_pdg_nue_cc_qe",     "h_leading_shower_mc_pdg_nue_cc_qe",     3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_out_fv = new TH2D("h_leading_shower_mc_pdg_nue_cc_out_fv", "h_leading_shower_mc_pdg_nue_cc_out_fv", 3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_res    = new TH2D("h_leading_shower_mc_pdg_nue_cc_res",    "h_leading_shower_mc_pdg_nue_cc_res",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_dis    = new TH2D("h_leading_shower_mc_pdg_nue_cc_dis",    "h_leading_shower_mc_pdg_nue_cc_dis",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_coh    = new TH2D("h_leading_shower_mc_pdg_nue_cc_coh",    "h_leading_shower_mc_pdg_nue_cc_coh",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_mec    = new TH2D("h_leading_shower_mc_pdg_nue_cc_mec",    "h_leading_shower_mc_pdg_nue_cc_mec",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nc            = new TH2D("h_leading_shower_mc_pdg_nc",            "h_leading_shower_mc_pdg_nc",            3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_qe    = new TH2D("h_leading_shower_mc_pdg_numu_cc_qe",    "h_leading_shower_mc_pdg_numu_cc_qe",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_res   = new TH2D("h_leading_shower_mc_pdg_numu_cc_res",   "h_leading_shower_mc_pdg_numu_cc_res",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_dis   = new TH2D("h_leading_shower_mc_pdg_numu_cc_dis",   "h_leading_shower_mc_pdg_numu_cc_dis",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_coh   = new TH2D("h_leading_shower_mc_pdg_numu_cc_coh",   "h_leading_shower_mc_pdg_numu_cc_coh",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_mec   = new TH2D("h_leading_shower_mc_pdg_numu_cc_mec",   "h_leading_shower_mc_pdg_numu_cc_mec",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nc_pi0        = new TH2D("h_leading_shower_mc_pdg_nc_pi0",        "h_leading_shower_mc_pdg_nc_pi0",        3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_mixed  = new TH2D("h_leading_shower_mc_pdg_nue_cc_mixed",  "h_leading_shower_mc_pdg_nue_cc_mixed",  3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_mixed = new TH2D("h_leading_shower_mc_pdg_numu_cc_mixed", "h_leading_shower_mc_pdg_numu_cc_mixed", 3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_cosmic        = new TH2D("h_leading_shower_mc_pdg_cosmic",        "h_leading_shower_mc_pdg_cosmic",        3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_other_mixed   = new TH2D("h_leading_shower_mc_pdg_other_mixed",   "h_leading_shower_mc_pdg_other_mixed",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_unmatched     = new TH2D("h_leading_shower_mc_pdg_unmatched",     "h_leading_shower_mc_pdg_unmatched",     3, 0, 3, 10, 0, 10);

TH1D * h_pfp_track_nue_cc_qe     = new TH1D("h_pfp_track_nue_cc_qe",     "h_pfp_track_nue_cc_qe",     10, 0, 10);
TH1D * h_pfp_track_nue_cc_out_fv = new TH1D("h_pfp_track_nue_cc_out_fv", "h_pfp_track_nue_cc_out_fv", 10, 0, 10);
TH1D * h_pfp_track_nue_cc_res    = new TH1D("h_pfp_track_nue_cc_res",    "h_pfp_track_nue_cc_res",    10, 0, 10);
TH1D * h_pfp_track_nue_cc_dis    = new TH1D("h_pfp_track_nue_cc_dis",    "h_pfp_track_nue_cc_dis",    10, 0, 10);
TH1D * h_pfp_track_nue_cc_coh    = new TH1D("h_pfp_track_nue_cc_coh",    "h_pfp_track_nue_cc_coh",    10, 0, 10);
TH1D * h_pfp_track_nue_cc_mec    = new TH1D("h_pfp_track_nue_cc_mec",    "h_pfp_track_nue_cc_mec",    10, 0, 10);
TH1D * h_pfp_track_nc            = new TH1D("h_pfp_track_nc",            "h_pfp_track_nc",            10, 0, 10);
TH1D * h_pfp_track_numu_cc_qe    = new TH1D("h_pfp_track_numu_cc_qe",    "h_pfp_track_numu_cc_qe",    10, 0, 10);
TH1D * h_pfp_track_numu_cc_res   = new TH1D("h_pfp_track_numu_cc_res",   "h_pfp_track_numu_cc_res",   10, 0, 10);
TH1D * h_pfp_track_numu_cc_dis   = new TH1D("h_pfp_track_numu_cc_dis",   "h_pfp_track_numu_cc_dis",   10, 0, 10);
TH1D * h_pfp_track_numu_cc_coh   = new TH1D("h_pfp_track_numu_cc_coh",   "h_pfp_track_numu_cc_coh",   10, 0, 10);
TH1D * h_pfp_track_numu_cc_mec   = new TH1D("h_pfp_track_numu_cc_mec",   "h_pfp_track_numu_cc_mec",   10, 0, 10);
TH1D * h_pfp_track_nc_pi0        = new TH1D("h_pfp_track_nc_pi0",        "h_pfp_track_nc_pi0",        10, 0, 10);
TH1D * h_pfp_track_nue_cc_mixed  = new TH1D("h_pfp_track_nue_cc_mixed",  "h_pfp_track_nue_cc_mixed",  10, 0, 10);
TH1D * h_pfp_track_numu_cc_mixed = new TH1D("h_pfp_track_numu_cc_mixed", "h_pfp_track_numu_cc_mixed", 10, 0, 10);
TH1D * h_pfp_track_cosmic        = new TH1D("h_pfp_track_cosmic",        "h_pfp_track_cosmic",        10, 0, 10);
TH1D * h_pfp_track_other_mixed   = new TH1D("h_pfp_track_other_mixed",   "h_pfp_track_other_mixed",   10, 0, 10);
TH1D * h_pfp_track_unmatched     = new TH1D("h_pfp_track_unmatched",     "h_pfp_track_unmatched",     10, 0, 10);

TH1D * h_pfp_shower_nue_cc_qe     = new TH1D("h_pfp_shower_nue_cc_qe",     "h_pfp_shower_nue_cc_qe",     10, 0, 10);
TH1D * h_pfp_shower_nue_cc_out_fv = new TH1D("h_pfp_shower_nue_cc_out_fv", "h_pfp_shower_nue_cc_out_fv", 10, 0, 10);
TH1D * h_pfp_shower_nue_cc_res    = new TH1D("h_pfp_shower_nue_cc_res",    "h_pfp_shower_nue_cc_res",    10, 0, 10);
TH1D * h_pfp_shower_nue_cc_dis    = new TH1D("h_pfp_shower_nue_cc_dis",    "h_pfp_shower_nue_cc_dis",    10, 0, 10);
TH1D * h_pfp_shower_nue_cc_coh    = new TH1D("h_pfp_shower_nue_cc_coh",    "h_pfp_shower_nue_cc_coh",    10, 0, 10);
TH1D * h_pfp_shower_nue_cc_mec    = new TH1D("h_pfp_shower_nue_cc_mec",    "h_pfp_shower_nue_cc_mec",    10, 0, 10);
TH1D * h_pfp_shower_nc            = new TH1D("h_pfp_shower_nc",            "h_pfp_shower_nc",            10, 0, 10);
TH1D * h_pfp_shower_numu_cc_qe    = new TH1D("h_pfp_shower_numu_cc_qe",    "h_pfp_shower_numu_cc_qe",    10, 0, 10);
TH1D * h_pfp_shower_numu_cc_res   = new TH1D("h_pfp_shower_numu_cc_res",   "h_pfp_shower_numu_cc_res",   10, 0, 10);
TH1D * h_pfp_shower_numu_cc_dis   = new TH1D("h_pfp_shower_numu_cc_dis",   "h_pfp_shower_numu_cc_dis",   10, 0, 10);
TH1D * h_pfp_shower_numu_cc_coh   = new TH1D("h_pfp_shower_numu_cc_coh",   "h_pfp_shower_numu_cc_coh",   10, 0, 10);
TH1D * h_pfp_shower_numu_cc_mec   = new TH1D("h_pfp_shower_numu_cc_mec",   "h_pfp_shower_numu_cc_mec",   10, 0, 10);
TH1D * h_pfp_shower_nc_pi0        = new TH1D("h_pfp_shower_nc_pi0",        "h_pfp_shower_nc_pi0",        10, 0, 10);
TH1D * h_pfp_shower_nue_cc_mixed  = new TH1D("h_pfp_shower_nue_cc_mixed",  "h_pfp_shower_nue_cc_mixed",  10, 0, 10);
TH1D * h_pfp_shower_numu_cc_mixed = new TH1D("h_pfp_shower_numu_cc_mixed", "h_pfp_shower_numu_cc_mixed", 10, 0, 10);
TH1D * h_pfp_shower_cosmic        = new TH1D("h_pfp_shower_cosmic",        "h_pfp_shower_cosmic",        10, 0, 10);
TH1D * h_pfp_shower_other_mixed   = new TH1D("h_pfp_shower_other_mixed",   "h_pfp_shower_other_mixed",   10, 0, 10);
TH1D * h_pfp_shower_unmatched     = new TH1D("h_pfp_shower_unmatched",     "h_pfp_shower_unmatched",     10, 0, 10);

TH1D * h_pfp_shower_open_angle_nue_cc_qe     = new TH1D("h_pfp_shower_open_angle_nue_cc_qe",     "h_pfp_shower_open_angle_nue_cc_qe",     10, 0, 10);
TH1D * h_pfp_shower_open_angle_nue_cc_out_fv = new TH1D("h_pfp_shower_open_angle_nue_cc_out_fv", "h_pfp_shower_open_angle_nue_cc_out_fv", 10, 0, 10);
TH1D * h_pfp_shower_open_angle_nue_cc_res    = new TH1D("h_pfp_shower_open_angle_nue_cc_res",    "h_pfp_shower_open_angle_nue_cc_res",    10, 0, 10);
TH1D * h_pfp_shower_open_angle_nue_cc_dis    = new TH1D("h_pfp_shower_open_angle_nue_cc_dis",    "h_pfp_shower_open_angle_nue_cc_dis",    10, 0, 10);
TH1D * h_pfp_shower_open_angle_nue_cc_coh    = new TH1D("h_pfp_shower_open_angle_nue_cc_coh",    "h_pfp_shower_open_angle_nue_cc_coh",    10, 0, 10);
TH1D * h_pfp_shower_open_angle_nue_cc_mec    = new TH1D("h_pfp_shower_open_angle_nue_cc_mec",    "h_pfp_shower_open_angle_nue_cc_mec",    10, 0, 10);
TH1D * h_pfp_shower_open_angle_nc            = new TH1D("h_pfp_shower_open_angle_nc",            "h_pfp_shower_open_angle_nc",            10, 0, 10);
TH1D * h_pfp_shower_open_angle_numu_cc_qe    = new TH1D("h_pfp_shower_open_angle_numu_cc_qe",    "h_pfp_shower_open_angle_numu_cc_qe",    10, 0, 10);
TH1D * h_pfp_shower_open_angle_numu_cc_res   = new TH1D("h_pfp_shower_open_angle_numu_cc_res",   "h_pfp_shower_open_angle_numu_cc_res",   10, 0, 10);
TH1D * h_pfp_shower_open_angle_numu_cc_dis   = new TH1D("h_pfp_shower_open_angle_numu_cc_dis",   "h_pfp_shower_open_angle_numu_cc_dis",   10, 0, 10);
TH1D * h_pfp_shower_open_angle_numu_cc_coh   = new TH1D("h_pfp_shower_open_angle_numu_cc_coh",   "h_pfp_shower_open_angle_numu_cc_coh",   10, 0, 10);
TH1D * h_pfp_shower_open_angle_numu_cc_mec   = new TH1D("h_pfp_shower_open_angle_numu_cc_mec",   "h_pfp_shower_open_angle_numu_cc_mec",   10, 0, 10);
TH1D * h_pfp_shower_open_angle_nc_pi0        = new TH1D("h_pfp_shower_open_angle_nc_pi0",        "h_pfp_shower_open_angle_nc_pi0",        10, 0, 10);
TH1D * h_pfp_shower_open_angle_nue_cc_mixed  = new TH1D("h_pfp_shower_open_angle_nue_cc_mixed",  "h_pfp_shower_open_angle_nue_cc_mixed",  10, 0, 10);
TH1D * h_pfp_shower_open_angle_numu_cc_mixed = new TH1D("h_pfp_shower_open_angle_numu_cc_mixed", "h_pfp_shower_open_angle_numu_cc_mixed", 10, 0, 10);
TH1D * h_pfp_shower_open_angle_cosmic        = new TH1D("h_pfp_shower_open_angle_cosmic",        "h_pfp_shower_open_angle_cosmic",        10, 0, 10);
TH1D * h_pfp_shower_open_angle_other_mixed   = new TH1D("h_pfp_shower_open_angle_other_mixed",   "h_pfp_shower_open_angle_other_mixed",   10, 0, 10);
TH1D * h_pfp_shower_open_angle_unmatched     = new TH1D("h_pfp_shower_open_angle_unmatched",     "h_pfp_shower_open_angle_unmatched",     10, 0, 10);
TH1D * h_pfp_shower_open_angle_intime        = new TH1D("h_pfp_shower_open_angle_intime",        "h_pfp_shower_open_angle_intime",        10, 0, 10);
TH1D * h_pfp_shower_open_angle_data          = new TH1D("h_pfp_shower_open_angle_data",          "h_pfp_shower_open_angle_data",          10, 0, 10);

TH1D * h_pfp_shower_dedx_nue_cc_qe     = new TH1D("h_pfp_shower_dedx_nue_cc_qe",     "h_pfp_dedx_nue_cc_qe",     10, 0, 10);
TH1D * h_pfp_shower_dedx_nue_cc_out_fv = new TH1D("h_pfp_shower_dedx_nue_cc_out_fv", "h_pfp_dedx_nue_cc_out_fv", 10, 0, 10);
TH1D * h_pfp_shower_dedx_nue_cc_res    = new TH1D("h_pfp_shower_dedx_nue_cc_res",    "h_pfp_dedx_nue_cc_res",    10, 0, 10);
TH1D * h_pfp_shower_dedx_nue_cc_dis    = new TH1D("h_pfp_shower_dedx_nue_cc_dis",    "h_pfp_dedx_nue_cc_dis",    10, 0, 10);
TH1D * h_pfp_shower_dedx_nue_cc_coh    = new TH1D("h_pfp_shower_dedx_nue_cc_coh",    "h_pfp_dedx_nue_cc_coh",    10, 0, 10);
TH1D * h_pfp_shower_dedx_nue_cc_mec    = new TH1D("h_pfp_shower_dedx_nue_cc_mec",    "h_pfp_dedx_nue_cc_mec",    10, 0, 10);
TH1D * h_pfp_shower_dedx_nc            = new TH1D("h_pfp_shower_dedx_nc",            "h_pfp_dedx_nc",            10, 0, 10);
TH1D * h_pfp_shower_dedx_numu_cc_qe    = new TH1D("h_pfp_shower_dedx_numu_cc_qe",    "h_pfp_dedx_numu_cc_qe",    10, 0, 10);
TH1D * h_pfp_shower_dedx_numu_cc_res   = new TH1D("h_pfp_shower_dedx_numu_cc_res",   "h_pfp_dedx_numu_cc_res",   10, 0, 10);
TH1D * h_pfp_shower_dedx_numu_cc_dis   = new TH1D("h_pfp_shower_dedx_numu_cc_dis",   "h_pfp_dedx_numu_cc_dis",   10, 0, 10);
TH1D * h_pfp_shower_dedx_numu_cc_coh   = new TH1D("h_pfp_shower_dedx_numu_cc_coh",   "h_pfp_dedx_numu_cc_coh",   10, 0, 10);
TH1D * h_pfp_shower_dedx_numu_cc_mec   = new TH1D("h_pfp_shower_dedx_numu_cc_mec",   "h_pfp_dedx_numu_cc_mec",   10, 0, 10);
TH1D * h_pfp_shower_dedx_nc_pi0        = new TH1D("h_pfp_shower_dedx_nc_pi0",        "h_pfp_dedx_nc_pi0",        10, 0, 10);
TH1D * h_pfp_shower_dedx_nue_cc_mixed  = new TH1D("h_pfp_shower_dedx_nue_cc_mixed",  "h_pfp_dedx_nue_cc_mixed",  10, 0, 10);
TH1D * h_pfp_shower_dedx_numu_cc_mixed = new TH1D("h_pfp_shower_dedx_numu_cc_mixed", "h_pfp_dedx_numu_cc_mixed", 10, 0, 10);
TH1D * h_pfp_shower_dedx_cosmic        = new TH1D("h_pfp_shower_dedx_cosmic",        "h_pfp_dedx_cosmic",        10, 0, 10);
TH1D * h_pfp_shower_dedx_other_mixed   = new TH1D("h_pfp_shower_dedx_other_mixed",   "h_pfp_dedx_other_mixed",   10, 0, 10);
TH1D * h_pfp_shower_dedx_unmatched     = new TH1D("h_pfp_shower_dedx_unmatched",     "h_pfp_dedx_unmatched",     10, 0, 10);
TH1D * h_pfp_shower_dedx_intime        = new TH1D("h_pfp_shower_dedx_intime",        "h_pfp_dedx_intime",        10, 0, 10);
TH1D * h_pfp_shower_dedx_data          = new TH1D("h_pfp_shower_dedx_data",          "h_pfp_dedx_data",          10, 0, 10);

TH2D * h_pfp_track_shower_nue_cc_qe_last     = new TH2D("h_pfp_track_shower_nue_cc_qe_last",     "h_pfp_track_shower_nue_cc_qe_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_out_fv_last = new TH2D("h_pfp_track_shower_nue_cc_out_fv_last", "h_pfp_track_shower_nue_cc_out_fv_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_res_last    = new TH2D("h_pfp_track_shower_nue_cc_res_last",    "h_pfp_track_shower_nue_cc_res_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_dis_last    = new TH2D("h_pfp_track_shower_nue_cc_dis_last",    "h_pfp_track_shower_nue_cc_dis_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_coh_last    = new TH2D("h_pfp_track_shower_nue_cc_coh_last",    "h_pfp_track_shower_nue_cc_coh_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_mec_last    = new TH2D("h_pfp_track_shower_nue_cc_mec_last",    "h_pfp_track_shower_nue_cc_mec_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nc_last            = new TH2D("h_pfp_track_shower_nc_last",            "h_pfp_track_shower_nc_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_qe_last    = new TH2D("h_pfp_track_shower_numu_cc_qe_last",    "h_pfp_track_shower_numu_cc_qe_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_res_last   = new TH2D("h_pfp_track_shower_numu_cc_res_last",   "h_pfp_track_shower_numu_cc_res_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_dis_last   = new TH2D("h_pfp_track_shower_numu_cc_dis_last",   "h_pfp_track_shower_numu_cc_dis_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_coh_last   = new TH2D("h_pfp_track_shower_numu_cc_coh_last",   "h_pfp_track_shower_numu_cc_coh_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_mec_last   = new TH2D("h_pfp_track_shower_numu_cc_mec_last",   "h_pfp_track_shower_numu_cc_mec_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nc_pi0_last        = new TH2D("h_pfp_track_shower_nc_pi0_last",        "h_pfp_track_shower_nc_pi0_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_mixed_last  = new TH2D("h_pfp_track_shower_nue_cc_mixed_last",  "h_pfp_track_shower_nue_cc_mixed_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_mixed_last = new TH2D("h_pfp_track_shower_numu_cc_mixed_last", "h_pfp_track_shower_numu_cc_mixed_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_cosmic_last        = new TH2D("h_pfp_track_shower_cosmic_last",        "h_pfp_track_shower_cosmic_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_other_mixed_last   = new TH2D("h_pfp_track_shower_other_mixed_last",   "h_pfp_track_shower_other_mixed_last", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_unmatched_last     = new TH2D("h_pfp_track_shower_unmatched_last",     "h_pfp_track_shower_unmatched_last", 10, 0, 10, 10, 0, 10);

TH2D * h_leading_shower_mc_pdg_nue_cc_qe_last     = new TH2D("h_leading_shower_mc_pdg_nue_cc_qe_last",     "h_leading_shower_mc_pdg_nue_cc_qe_last",     3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_out_fv_last = new TH2D("h_leading_shower_mc_pdg_nue_cc_out_fv_last", "h_leading_shower_mc_pdg_nue_cc_out_fv_last", 3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_res_last    = new TH2D("h_leading_shower_mc_pdg_nue_cc_res_last",    "h_leading_shower_mc_pdg_nue_cc_res_last",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_dis_last    = new TH2D("h_leading_shower_mc_pdg_nue_cc_dis_last",    "h_leading_shower_mc_pdg_nue_cc_dis_last",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_coh_last    = new TH2D("h_leading_shower_mc_pdg_nue_cc_coh_last",    "h_leading_shower_mc_pdg_nue_cc_coh_last",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_mec_last    = new TH2D("h_leading_shower_mc_pdg_nue_cc_mec_last",    "h_leading_shower_mc_pdg_nue_cc_mec_last",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nc_last            = new TH2D("h_leading_shower_mc_pdg_nc_last",            "h_leading_shower_mc_pdg_nc_last",            3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_qe_last    = new TH2D("h_leading_shower_mc_pdg_numu_cc_qe_last",    "h_leading_shower_mc_pdg_numu_cc_qe_last",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_res_last   = new TH2D("h_leading_shower_mc_pdg_numu_cc_res_last",   "h_leading_shower_mc_pdg_numu_cc_res_last",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_dis_last   = new TH2D("h_leading_shower_mc_pdg_numu_cc_dis_last",   "h_leading_shower_mc_pdg_numu_cc_dis_last",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_coh_last   = new TH2D("h_leading_shower_mc_pdg_numu_cc_coh_last",   "h_leading_shower_mc_pdg_numu_cc_coh_last",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_mec_last   = new TH2D("h_leading_shower_mc_pdg_numu_cc_mec_last",   "h_leading_shower_mc_pdg_numu_cc_mec_last",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nc_pi0_last        = new TH2D("h_leading_shower_mc_pdg_nc_pi0_last",        "h_leading_shower_mc_pdg_nc_pi0_last",        3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_nue_cc_mixed_last  = new TH2D("h_leading_shower_mc_pdg_nue_cc_mixed_last",  "h_leading_shower_mc_pdg_nue_cc_mixed_last",  3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_mixed_last = new TH2D("h_leading_shower_mc_pdg_numu_cc_mixed_last", "h_leading_shower_mc_pdg_numu_cc_mixed_last", 3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_cosmic_last        = new TH2D("h_leading_shower_mc_pdg_cosmic_last",        "h_leading_shower_mc_pdg_cosmic_last",        3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_other_mixed_last   = new TH2D("h_leading_shower_mc_pdg_other_mixed_last",   "h_leading_shower_mc_pdg_other_mixed_last",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_unmatched_last     = new TH2D("h_leading_shower_mc_pdg_unmatched_last",     "h_leading_shower_mc_pdg_unmatched_last",     3, 0, 3, 10, 0, 10);

TH1D * h_pfp_track_nue_cc_qe_last     = new TH1D("h_pfp_track_nue_cc_qe_last",     "h_pfp_track_nue_cc_qe_last",     10, 0, 10);
TH1D * h_pfp_track_nue_cc_out_fv_last = new TH1D("h_pfp_track_nue_cc_out_fv_last", "h_pfp_track_nue_cc_out_fv_last", 10, 0, 10);
TH1D * h_pfp_track_nue_cc_res_last    = new TH1D("h_pfp_track_nue_cc_res_last",    "h_pfp_track_nue_cc_res_last",    10, 0, 10);
TH1D * h_pfp_track_nue_cc_dis_last    = new TH1D("h_pfp_track_nue_cc_dis_last",    "h_pfp_track_nue_cc_dis_last",    10, 0, 10);
TH1D * h_pfp_track_nue_cc_coh_last    = new TH1D("h_pfp_track_nue_cc_coh_last",    "h_pfp_track_nue_cc_coh_last",    10, 0, 10);
TH1D * h_pfp_track_nue_cc_mec_last    = new TH1D("h_pfp_track_nue_cc_mec_last",    "h_pfp_track_nue_cc_mec_last",    10, 0, 10);
TH1D * h_pfp_track_nc_last            = new TH1D("h_pfp_track_nc_last",            "h_pfp_track_nc_last",            10, 0, 10);
TH1D * h_pfp_track_numu_cc_qe_last    = new TH1D("h_pfp_track_numu_cc_qe_last",    "h_pfp_track_numu_cc_qe_last",    10, 0, 10);
TH1D * h_pfp_track_numu_cc_res_last   = new TH1D("h_pfp_track_numu_cc_res_last",   "h_pfp_track_numu_cc_res_last",   10, 0, 10);
TH1D * h_pfp_track_numu_cc_dis_last   = new TH1D("h_pfp_track_numu_cc_dis_last",   "h_pfp_track_numu_cc_dis_last",   10, 0, 10);
TH1D * h_pfp_track_numu_cc_coh_last   = new TH1D("h_pfp_track_numu_cc_coh_last",   "h_pfp_track_numu_cc_coh_last",   10, 0, 10);
TH1D * h_pfp_track_numu_cc_mec_last   = new TH1D("h_pfp_track_numu_cc_mec_last",   "h_pfp_track_numu_cc_mec_last",   10, 0, 10);
TH1D * h_pfp_track_nc_pi0_last        = new TH1D("h_pfp_track_nc_pi0_last",        "h_pfp_track_nc_pi0_last",        10, 0, 10);
TH1D * h_pfp_track_nue_cc_mixed_last  = new TH1D("h_pfp_track_nue_cc_mixed_last",  "h_pfp_track_nue_cc_mixed_last",  10, 0, 10);
TH1D * h_pfp_track_numu_cc_mixed_last = new TH1D("h_pfp_track_numu_cc_mixed_last", "h_pfp_track_numu_cc_mixed_last", 10, 0, 10);
TH1D * h_pfp_track_cosmic_last        = new TH1D("h_pfp_track_cosmic_last",        "h_pfp_track_cosmic_last",        10, 0, 10);
TH1D * h_pfp_track_other_mixed_last   = new TH1D("h_pfp_track_other_mixed_last",   "h_pfp_track_other_mixed_last",   10, 0, 10);
TH1D * h_pfp_track_unmatched_last     = new TH1D("h_pfp_track_unmatched_last",     "h_pfp_track_unmatched_last",     10, 0, 10);

TH1D * h_pfp_shower_nue_cc_qe_last     = new TH1D("h_pfp_shower_nue_cc_qe_last",     "h_pfp_shower_nue_cc_qe_last",     10, 0, 10);
TH1D * h_pfp_shower_nue_cc_out_fv_last = new TH1D("h_pfp_shower_nue_cc_out_fv_last", "h_pfp_shower_nue_cc_out_fv_last", 10, 0, 10);
TH1D * h_pfp_shower_nue_cc_res_last    = new TH1D("h_pfp_shower_nue_cc_res_last",    "h_pfp_shower_nue_cc_res_last",    10, 0, 10);
TH1D * h_pfp_shower_nue_cc_dis_last    = new TH1D("h_pfp_shower_nue_cc_dis_last",    "h_pfp_shower_nue_cc_dis_last",    10, 0, 10);
TH1D * h_pfp_shower_nue_cc_coh_last    = new TH1D("h_pfp_shower_nue_cc_coh_last",    "h_pfp_shower_nue_cc_coh_last",    10, 0, 10);
TH1D * h_pfp_shower_nue_cc_mec_last    = new TH1D("h_pfp_shower_nue_cc_mec_last",    "h_pfp_shower_nue_cc_mec_last",    10, 0, 10);
TH1D * h_pfp_shower_nc_last            = new TH1D("h_pfp_shower_nc_last",            "h_pfp_shower_nc_last",            10, 0, 10);
TH1D * h_pfp_shower_numu_cc_qe_last    = new TH1D("h_pfp_shower_numu_cc_qe_last",    "h_pfp_shower_numu_cc_qe_last",    10, 0, 10);
TH1D * h_pfp_shower_numu_cc_res_last   = new TH1D("h_pfp_shower_numu_cc_res_last",   "h_pfp_shower_numu_cc_res_last",   10, 0, 10);
TH1D * h_pfp_shower_numu_cc_dis_last   = new TH1D("h_pfp_shower_numu_cc_dis_last",   "h_pfp_shower_numu_cc_dis_last",   10, 0, 10);
TH1D * h_pfp_shower_numu_cc_coh_last   = new TH1D("h_pfp_shower_numu_cc_coh_last",   "h_pfp_shower_numu_cc_coh_last",   10, 0, 10);
TH1D * h_pfp_shower_numu_cc_mec_last   = new TH1D("h_pfp_shower_numu_cc_mec_last",   "h_pfp_shower_numu_cc_mec_last",   10, 0, 10);
TH1D * h_pfp_shower_nc_pi0_last        = new TH1D("h_pfp_shower_nc_pi0_last",        "h_pfp_shower_nc_pi0_last",        10, 0, 10);
TH1D * h_pfp_shower_nue_cc_mixed_last  = new TH1D("h_pfp_shower_nue_cc_mixed_last",  "h_pfp_shower_nue_cc_mixed_last",  10, 0, 10);
TH1D * h_pfp_shower_numu_cc_mixed_last = new TH1D("h_pfp_shower_numu_cc_mixed_last", "h_pfp_shower_numu_cc_mixed_last", 10, 0, 10);
TH1D * h_pfp_shower_cosmic_last        = new TH1D("h_pfp_shower_cosmic_last",        "h_pfp_shower_cosmic_last",        10, 0, 10);
TH1D * h_pfp_shower_other_mixed_last   = new TH1D("h_pfp_shower_other_mixed_last",   "h_pfp_shower_other_mixed_last",   10, 0, 10);
TH1D * h_pfp_shower_unmatched_last     = new TH1D("h_pfp_shower_unmatched_last",     "h_pfp_shower_unmatched_last",     10, 0, 10);


TH1D * h_vtx_flash_nue_cc        = new TH1D("h_vtx_flash_nue_cc",         "h_vtx_flash_nue_cc",        40, 0, 200);
TH1D * h_vtx_flash_nue_cc_mixed  = new TH1D("h_vtx_flash_nue_cc_mixed",   "h_vtx_flash_nue_cc_mixed",  40, 0, 200);
TH1D * h_vtx_flash_nue_cc_out_fv = new TH1D("h_vtx_flash_nue_cc_out_fv",  "h_vtx_flash_nue_cc_out_fv", 40, 0, 200);
TH1D * h_vtx_flash_numu_cc       = new TH1D("h_vtx_flash_numu_cc",        "h_vtx_flash_numu_cc_mixed", 40, 0, 200);
TH1D * h_vtx_flash_nc_pi0        = new TH1D("h_vtx_flash_nc_pi0",         "h_vtx_flash_nc_pi0",        40, 0, 200);
TH1D * h_vtx_flash_cosmic        = new TH1D("h_vtx_flash_cosmic",         "h_vtx_flash_cosmic",        40, 0, 200);
TH1D * h_vtx_flash_nc            = new TH1D("h_vtx_flash_nc",             "h_vtx_flash_nc",            40, 0, 200);
TH1D * h_vtx_flash_numu_cc_mixed = new TH1D("h_vtx_flash_numu_cc_mixed",  "h_vtx_flash_numu_cc_mixed", 40, 0, 200);
TH1D * h_vtx_flash_other_mixed   = new TH1D("h_vtx_flash_other_mixed",    "h_vtx_flash_other_mixed",   40, 0, 200);
TH1D * h_vtx_flash_unmatched     = new TH1D("h_vtx_flash_unmatched",      "h_vtx_flash_unmatched",     40, 0, 200);
TH1D * h_vtx_flash_intime        = new TH1D("h_vtx_flash_intime",         "h_vtx_flash_intime",        40, 0, 200);
TH1D * h_vtx_flash_data          = new TH1D("h_vtx_flash_data",           "h_vtx_flash_data",          40, 0, 200);

TH1D * h_vtx_flash_nue_cc_after        = new TH1D("h_vtx_flash_nue_cc_after",         "h_vtx_flash_nue_cc_after",        40, 0, 200);
TH1D * h_vtx_flash_nue_cc_mixed_after  = new TH1D("h_vtx_flash_nue_cc_mixed_after",   "h_vtx_flash_nue_cc_mixed_after",  40, 0, 200);
TH1D * h_vtx_flash_nue_cc_out_fv_after = new TH1D("h_vtx_flash_nue_cc_out_fv_after",  "h_vtx_flash_nue_cc_out_fv_after", 40, 0, 200);
TH1D * h_vtx_flash_numu_cc_after       = new TH1D("h_vtx_flash_numu_cc_after",        "h_vtx_flash_numu_cc_mixed_after", 40, 0, 200);
TH1D * h_vtx_flash_nc_pi0_after        = new TH1D("h_vtx_flash_nc_pi0_after",         "h_vtx_flash_nc_pi0_after",        40, 0, 200);
TH1D * h_vtx_flash_cosmic_after        = new TH1D("h_vtx_flash_cosmic_after",         "h_vtx_flash_cosmic_after",        40, 0, 200);
TH1D * h_vtx_flash_nc_after            = new TH1D("h_vtx_flash_nc_after",             "h_vtx_flash_nc_after",            40, 0, 200);
TH1D * h_vtx_flash_numu_cc_mixed_after = new TH1D("h_vtx_flash_numu_cc_mixed_after",  "h_vtx_flash_numu_cc_mixed_after", 40, 0, 200);
TH1D * h_vtx_flash_other_mixed_after   = new TH1D("h_vtx_flash_other_mixed_after",    "h_vtx_flash_other_mixed_after",   40, 0, 200);
TH1D * h_vtx_flash_unmatched_after     = new TH1D("h_vtx_flash_unmatched_after",      "h_vtx_flash_unmatched_after",     40, 0, 200);
TH1D * h_vtx_flash_intime_after        = new TH1D("h_vtx_flash_intime_after",         "h_vtx_flash_intime_after",        40, 0, 200);
TH1D * h_vtx_flash_data_after          = new TH1D("h_vtx_flash_data_after",           "h_vtx_flash_data_after",          40, 0, 200);

TH1D * h_shwr_vtx_dist_nue_cc        = new TH1D("h_shwr_vtx_dist_nue_cc",         "h_shwr_vtx_dist_nue_cc",        20, 0, 20);
TH1D * h_shwr_vtx_dist_nue_cc_mixed  = new TH1D("h_shwr_vtx_dist_nue_cc_mixed",   "h_shwr_vtx_dist_nue_cc_mixed",  20, 0, 20);
TH1D * h_shwr_vtx_dist_nue_cc_out_fv = new TH1D("h_shwr_vtx_dist_nue_cc_out_fv",  "h_shwr_vtx_dist_nue_cc_out_fv", 20, 0, 20);
TH1D * h_shwr_vtx_dist_numu_cc       = new TH1D("h_shwr_vtx_dist_numu_cc",        "h_shwr_vtx_dist_numu_cc",       20, 0, 20);
TH1D * h_shwr_vtx_dist_nc_pi0        = new TH1D("h_shwr_vtx_dist_nc_pi0",         "h_shwr_vtx_dist_nc_pi0",        20, 0, 20);
TH1D * h_shwr_vtx_dist_cosmic        = new TH1D("h_shwr_vtx_dist_cosmic",         "h_shwr_vtx_dist_cosmic",        20, 0, 20);
TH1D * h_shwr_vtx_dist_nc            = new TH1D("h_shwr_vtx_dist_nc",             "h_shwr_vtx_dist_nc",            20, 0, 20);
TH1D * h_shwr_vtx_dist_numu_cc_mixed = new TH1D("h_shwr_vtx_dist_numu_cc_mixed",  "h_shwr_vtx_dist_numu_cc_mixed", 20, 0, 20);
TH1D * h_shwr_vtx_dist_other_mixed   = new TH1D("h_shwr_vtx_dist_other_mixed",    "h_shwr_vtx_dist_other_mixed",   20, 0, 20);
TH1D * h_shwr_vtx_dist_unmatched     = new TH1D("h_shwr_vtx_dist_unmatched",      "h_shwr_vtx_dist_unmatched",     20, 0, 20);
TH1D * h_shwr_vtx_dist_intime        = new TH1D("h_shwr_vtx_dist_intime",         "h_shwr_vtx_dist_intime",        20, 0, 20);
TH1D * h_shwr_vtx_dist_data          = new TH1D("h_shwr_vtx_dist_data",           "h_shwr_vtx_dist_data",          20, 0, 20);

TH1D * h_shwr_vtx_dist_nue_cc_after        = new TH1D("h_shwr_vtx_dist_nue_cc_after",         "h_shwr_vtx_dist_nue_cc_after",        20, 0, 20);
TH1D * h_shwr_vtx_dist_nue_cc_mixed_after  = new TH1D("h_shwr_vtx_dist_nue_cc_mixed_after",   "h_shwr_vtx_dist_nue_cc_mixed_after",  20, 0, 20);
TH1D * h_shwr_vtx_dist_nue_cc_out_fv_after = new TH1D("h_shwr_vtx_dist_nue_cc_out_fv_after",  "h_shwr_vtx_dist_nue_cc_out_fv_after", 20, 0, 20);
TH1D * h_shwr_vtx_dist_numu_cc_after       = new TH1D("h_shwr_vtx_dist_numu_cc_after",        "h_shwr_vtx_dist_numu_cc_after",       20, 0, 20);
TH1D * h_shwr_vtx_dist_nc_pi0_after        = new TH1D("h_shwr_vtx_dist_nc_pi0_after",         "h_shwr_vtx_dist_nc_pi0_after",        20, 0, 20);
TH1D * h_shwr_vtx_dist_cosmic_after        = new TH1D("h_shwr_vtx_dist_cosmic_after",         "h_shwr_vtx_dist_cosmic_after",        20, 0, 20);
TH1D * h_shwr_vtx_dist_nc_after            = new TH1D("h_shwr_vtx_dist_nc_after",             "h_shwr_vtx_dist_nc_after",            20, 0, 20);
TH1D * h_shwr_vtx_dist_numu_cc_mixed_after = new TH1D("h_shwr_vtx_dist_numu_cc_mixed_after",  "h_shwr_vtx_dist_numu_cc_mixed_after", 20, 0, 20);
TH1D * h_shwr_vtx_dist_other_mixed_after   = new TH1D("h_shwr_vtx_dist_other_mixed_after",    "h_shwr_vtx_dist_other_mixed_after",   20, 0, 20);
TH1D * h_shwr_vtx_dist_unmatched_after     = new TH1D("h_shwr_vtx_dist_unmatched_after",      "h_shwr_vtx_dist_unmatched_after",     20, 0, 20);
TH1D * h_shwr_vtx_dist_intime_after        = new TH1D("h_shwr_vtx_dist_intime_after",         "h_shwr_vtx_dist_intime_after",        20, 0, 20);
TH1D * h_shwr_vtx_dist_data_after          = new TH1D("h_shwr_vtx_dist_data_after",           "h_shwr_vtx_dist_data_after",          20, 0, 20);

TH1D * h_dedx_cuts_nue_cc        = new TH1D("h_dedx_cuts_nue_cc",         "h_dedx_cuts_nue_cc",        20, 0, 10);
TH1D * h_dedx_cuts_nue_cc_mixed  = new TH1D("h_dedx_cuts_nue_cc_mixed",   "h_dedx_cuts_nue_cc_mixed",  20, 0, 10);
TH1D * h_dedx_cuts_nue_cc_out_fv = new TH1D("h_dedx_cuts_nue_cc_out_fv",  "h_dedx_cuts_nue_cc_out_fv", 20, 0, 10);
TH1D * h_dedx_cuts_numu_cc       = new TH1D("h_dedx_cuts_numu_cc",        "h_dedx_cuts_numu_cc",       20, 0, 10);
TH1D * h_dedx_cuts_nc_pi0        = new TH1D("h_dedx_cuts_nc_pi0",         "h_dedx_cuts_nc_pi0",        20, 0, 10);
TH1D * h_dedx_cuts_cosmic        = new TH1D("h_dedx_cuts_cosmic",         "h_dedx_cuts_cosmic",        20, 0, 10);
TH1D * h_dedx_cuts_nc            = new TH1D("h_dedx_cuts_nc",             "h_dedx_cuts_nc",            20, 0, 10);
TH1D * h_dedx_cuts_numu_cc_mixed = new TH1D("h_dedx_cuts_numu_cc_mixed",  "h_dedx_cuts_numu_cc_mixed", 20, 0, 10);
TH1D * h_dedx_cuts_other_mixed   = new TH1D("h_dedx_cuts_other_mixed",    "h_dedx_cuts_other_mixed",   20, 0, 10);
TH1D * h_dedx_cuts_unmatched     = new TH1D("h_dedx_cuts_unmatched",      "h_dedx_cuts_unmatched",     20, 0, 10);
TH1D * h_dedx_cuts_intime        = new TH1D("h_dedx_cuts_intime",         "h_dedx_cuts_intime",        20, 0, 10);
TH1D * h_dedx_cuts_data          = new TH1D("h_dedx_cuts_data",           "h_dedx_cuts_data",          20, 0, 10);

TH1D * h_dedx_cuts_nue_cc_after        = new TH1D("h_dedx_cuts_nue_cc_after",         "h_dedx_cuts_nue_cc_after",        20, 0, 10);
TH1D * h_dedx_cuts_nue_cc_mixed_after  = new TH1D("h_dedx_cuts_nue_cc_mixed_after",   "h_dedx_cuts_nue_cc_mixed_after",  20, 0, 10);
TH1D * h_dedx_cuts_nue_cc_out_fv_after = new TH1D("h_dedx_cuts_nue_cc_out_fv_after",  "h_dedx_cuts_nue_cc_out_fv_after", 20, 0, 10);
TH1D * h_dedx_cuts_numu_cc_after       = new TH1D("h_dedx_cuts_numu_cc_after",        "h_dedx_cuts_numu_cc_after",       20, 0, 10);
TH1D * h_dedx_cuts_nc_pi0_after        = new TH1D("h_dedx_cuts_nc_pi0_after",         "h_dedx_cuts_nc_pi0_after",        20, 0, 10);
TH1D * h_dedx_cuts_cosmic_after        = new TH1D("h_dedx_cuts_cosmic_after",         "h_dedx_cuts_cosmic_after",        20, 0, 10);
TH1D * h_dedx_cuts_nc_after            = new TH1D("h_dedx_cuts_nc_after",             "h_dedx_cuts_nc_after",            20, 0, 10);
TH1D * h_dedx_cuts_numu_cc_mixed_after = new TH1D("h_dedx_cuts_numu_cc_mixed_after",  "h_dedx_cuts_numu_cc_mixed_after", 20, 0, 10);
TH1D * h_dedx_cuts_other_mixed_after   = new TH1D("h_dedx_cuts_other_mixed_after",    "h_dedx_cuts_other_mixed_after",   20, 0, 10);
TH1D * h_dedx_cuts_unmatched_after     = new TH1D("h_dedx_cuts_unmatched_after",      "h_dedx_cuts_unmatched_after",     20, 0, 10);
TH1D * h_dedx_cuts_intime_after        = new TH1D("h_dedx_cuts_intime_after",         "h_dedx_cuts_intime_after",        20, 0, 10);
TH1D * h_dedx_cuts_data_after          = new TH1D("h_dedx_cuts_data_after",           "h_dedx_cuts_data_after",          20, 0, 10);

TH2D * h_shwr_hits_nu_eng_zoom  = new TH2D ("h_shwr_hits_nu_eng_zoom",  "h_shwr_hits_nu_eng_zoom", 20, 0, 4, 20, 0, 500);
TH2D * h_shwr_hits_ele_eng_zoom = new TH2D ("h_shwr_hits_ele_eng_zoom", "h_shwr_hits_ele_eng_zoom", 20, 0, 2, 20, 0, 500);
TH2D * h_shwr_hits_nu_eng       = new TH2D ("h_shwr_hits_nu_eng",       "h_shwr_hits_nu_eng", 20, 0, 4, 20, 0, 5000);
TH2D * h_shwr_hits_ele_eng      = new TH2D ("h_shwr_hits_ele_eng",      "h_shwr_hits_ele_eng", 20, 0, 2, 20, 0, 5000);

TH2D * h_shwr_hits_nu_eng_zoom_last  = new TH2D ("h_shwr_hits_nu_eng_zoom_last",  "h_shwr_hits_nu_eng_zoom_last", 20, 0, 4, 20, 0, 500);
TH2D * h_shwr_hits_ele_eng_zoom_last = new TH2D ("h_shwr_hits_ele_eng_zoom_last", "h_shwr_hits_ele_eng_zoom_last", 20, 0, 2, 20, 0, 500);
TH2D * h_shwr_hits_nu_eng_last       = new TH2D ("h_shwr_hits_nu_eng_last",       "h_shwr_hits_nu_eng_last", 20, 0, 4, 20, 0, 5000);
TH2D * h_shwr_hits_ele_eng_last      = new TH2D ("h_shwr_hits_ele_eng_last",      "h_shwr_hits_ele_eng_last", 20, 0, 2, 20, 0, 5000);


TH1D * h_selected_nu_energy_no_cut         = new TH1D ("h_selected_nu_energy_no_cut",         "h_selected_nu_energy_no_cut", 20, 0, 5);
TH1D * h_selected_ele_energy_no_cut        = new TH1D ("h_selected_ele_energy_no_cut",        "h_selected_ele_energy_no_cut", 20, 0, 5);
TH1D * h_selected_nu_energy_reco_nue       = new TH1D ("h_selected_nu_energy_reco_nue",       "h_selected_nu_energy_reco_nue", 20, 0, 5);
TH1D * h_selected_ele_energy_reco_nue      = new TH1D ("h_selected_ele_energy_reco_nue",      "h_selected_ele_energy_reco_nue", 20, 0, 5);
TH1D * h_selected_nu_energy_in_fv          = new TH1D ("h_selected_nu_energy_in_fv",          "h_selected_nu_energy_in_fv", 20, 0, 5);
TH1D * h_selected_ele_energy_in_fv         = new TH1D ("h_selected_ele_energy_in_fv",         "h_selected_ele_energy_in_fv", 20, 0, 5);
TH1D * h_selected_nu_energy_vtx_flash      = new TH1D ("h_selected_nu_energy_vtx_flash",      "h_selected_nu_energy_vtx_flash", 20, 0, 5);
TH1D * h_selected_ele_energy_vtx_flash     = new TH1D ("h_selected_ele_energy_vtx_flash",     "h_selected_ele_energy_vtx_flash", 20, 0, 5);
TH1D * h_selected_nu_energy_shwr_vtx       = new TH1D ("h_selected_nu_energy_shwr_vtx",       "h_selected_nu_energy_shwr_vtx", 20, 0, 5);
TH1D * h_selected_ele_energy_shwr_vtx      = new TH1D ("h_selected_ele_energy_shwr_vtx",      "h_selected_ele_energy_shwr_vtx", 20, 0, 5);
TH1D * h_selected_nu_energy_trk_vtx        = new TH1D ("h_selected_nu_energy_trk_vtx",        "h_selected_nu_energy_trk_vtx", 20, 0, 5);
TH1D * h_selected_ele_energy_trk_vtx       = new TH1D ("h_selected_ele_energy_trk_vtx",       "h_selected_ele_energy_trk_vtx", 20, 0, 5);
TH1D * h_selected_nu_energy_hit_threshold  = new TH1D ("h_selected_nu_energy_hit_threshold",  "h_selected_nu_energy_nu_threshold", 20, 0, 5);
TH1D * h_selected_ele_energy_hit_threshold = new TH1D ("h_selected_ele_energy_hit_threshold", "h_selected_ele_energy_hit_threshold", 20, 0, 5);
TH1D * h_selected_nu_energy_open_angle     = new TH1D ("h_selected_nu_energy_open_angle",     "h_selected_nu_energy_open_angle", 20, 0, 5);
TH1D * h_selected_ele_energy_open_angle    = new TH1D ("h_selected_ele_energy_open_angle",    "h_selected_ele_energy_open_angle", 20, 0, 5);
TH1D * h_selected_nu_energy_dedx           = new TH1D ("h_selected_nu_energy_dedx",           "h_selected_nu_energy_dedx", 20, 0, 5);
TH1D * h_selected_ele_energy_dedx          = new TH1D ("h_selected_ele_energy_dedx",          "h_selected_ele_energy_dedx", 20, 0, 5);

TH1D * h_charge_share_nue_cc_mixed   = new TH1D("h_charge_share_nue_cc_mixed", "h_charge_share_nue_cc_mixed", 10, 0, 1);

TH1D * h_flash_t0_diff = new TH1D ("h_flash_t0_diff", "h_flash_t0_diff", 40, -10, 10);


TH2D * h_dedx_open_angle_nue_cc        = new TH2D("h_dedx_open_angle_nue_cc",         "h_dedx_open_angle_nue_cc",        20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_nue_cc_mixed  = new TH2D("h_dedx_open_angle_nue_cc_mixed",   "h_dedx_open_angle_nue_cc_mixed",  20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_nue_cc_out_fv = new TH2D("h_dedx_open_angle_nue_cc_out_fv",  "h_dedx_open_angle_nue_cc_out_fv", 20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_numu_cc       = new TH2D("h_dedx_open_angle_numu_cc",        "h_dedx_open_angle_numu_cc",       20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_nc_pi0        = new TH2D("h_dedx_open_angle_nc_pi0",         "h_dedx_open_angle_nc_pi0",        20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_cosmic        = new TH2D("h_dedx_open_angle_cosmic",         "h_dedx_open_angle_cosmic",        20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_nc            = new TH2D("h_dedx_open_angle_nc",             "h_dedx_open_angle_nc",            20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_numu_cc_mixed = new TH2D("h_dedx_open_angle_numu_cc_mixed",  "h_dedx_open_angle_numu_cc_mixed", 20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_other_mixed   = new TH2D("h_dedx_open_angle_other_mixed",    "h_dedx_open_angle_other_mixed",   20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_unmatched     = new TH2D("h_dedx_open_angle_unmatched",      "h_dedx_open_angle_unmatched",     20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_intime        = new TH2D("h_dedx_open_angle_intime",         "h_dedx_open_angle_intime",        20, 0, 10, 20, 0, 45);
TH2D * h_dedx_open_angle_data          = new TH2D("h_dedx_open_angle_data",           "h_dedx_open_angle_data",          20, 0, 10, 20, 0, 45);

TH2D * h_shwr_len_hits_nue_cc          = new TH2D("h_shwr_len_hits_nue_cc",        "h_shwr_len_hits_nue_cc",          20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_nue_cc_out_fv   = new TH2D("h_shwr_len_hits_nue_cc_out_fv", "h_shwr_len_hits_nue_cc_out_fv",   20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_nue_cc_mixed    = new TH2D("h_shwr_len_hits_nue_cc_mixed",  "h_shwr_len_hits_nue_cc_mixed",    20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_numu_cc         = new TH2D("h_shwr_len_hits_numu_cc",       "h_shwr_len_hits_numu_cc",         20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_numu_cc_mixed   = new TH2D("h_shwr_len_hits_numu_cc_mixed", "h_shwr_len_hits_numu_cc_mixed",   20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_nc              = new TH2D("h_shwr_len_hits_nc",            "h_shwr_len_hits_nc",              20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_nc_pi0          = new TH2D("h_shwr_len_hits_nc_pi0",        "h_shwr_len_hits_nc_pi0",          20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_cosmic          = new TH2D("h_shwr_len_hits_cosmic",        "h_shwr_len_hits_cosmic",          20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_other_mixed     = new TH2D("h_shwr_len_hits_other_mixed",   "h_shwr_len_hits_other_mixed",     20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_unmatched       = new TH2D("h_shwr_len_hits_unmatched",     "h_shwr_len_hits_unmatched",       20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_intime          = new TH2D("h_shwr_len_hits_intime",        "h_shwr_len_hits_intime",          20, 0, 50, 20, 0, 500);
TH2D * h_shwr_len_hits_data            = new TH2D("h_shwr_len_hits_data",          "h_shwr_len_hits_data",            20, 0, 50, 20, 0, 500);

TH1D * h_second_shwr_dist_nue_cc         = new TH1D ("h_second_shwr_dist_nue_cc",         "h_second_shwr_dist_nue_cc",         20, 0, 160);
TH1D * h_second_shwr_dist_nue_cc_out_fv  = new TH1D ("h_second_shwr_dist_nue_cc_out_fv",  "h_second_shwr_dist_nue_cc_out_fv",  20, 0, 160);
TH1D * h_second_shwr_dist_nue_cc_mixed   = new TH1D ("h_second_shwr_dist_nue_cc_mixed",   "h_second_shwr_dist_nue_cc_mixed",   20, 0, 160);
TH1D * h_second_shwr_dist_numu_cc        = new TH1D ("h_second_shwr_dist_numu_cc",        "h_second_shwr_dist_numu_cc",        20, 0, 160);
TH1D * h_second_shwr_dist_numu_cc_mixed  = new TH1D ("h_second_shwr_dist_numu_cc_mixed",  "h_second_shwr_dist_numu_cc_mixed",  20, 0, 160);
TH1D * h_second_shwr_dist_nc             = new TH1D ("h_second_shwr_dist_nc",             "h_second_shwr_dist_nc",             20, 0, 160);
TH1D * h_second_shwr_dist_nc_pi0         = new TH1D ("h_second_shwr_dist_nc_pi0",         "h_second_shwr_dist_nc_pi0",         20, 0, 160);
TH1D * h_second_shwr_dist_cosmic         = new TH1D ("h_second_shwr_dist_cosmic",         "h_second_shwr_dist_cosmic",         20, 0, 160);
TH1D * h_second_shwr_dist_other_mixed    = new TH1D ("h_second_shwr_dist_other_mixed",    "h_second_shwr_dist_other_mixed",    20, 0, 160);
TH1D * h_second_shwr_dist_unmatched      = new TH1D ("h_second_shwr_dist_unmatched",      "h_second_shwr_dist_unmatched",      20, 0, 160);
TH1D * h_second_shwr_dist_intime         = new TH1D ("h_second_shwr_dist_intime",         "h_second_shwr_dist_intime",         20, 0, 160);
TH1D * h_second_shwr_dist_data           = new TH1D ("h_second_shwr_dist_data",           "h_second_shwr_dist_data",           20, 0, 160);

TH1D * h_second_shwr_dist_nue_cc_after         = new TH1D ("h_second_shwr_dist_nue_cc_after",         "h_second_shwr_dist_nue_cc_after",         20, 0, 160);
TH1D * h_second_shwr_dist_nue_cc_out_fv_after  = new TH1D ("h_second_shwr_dist_nue_cc_out_fv_after",  "h_second_shwr_dist_nue_cc_out_fv_after",  20, 0, 160);
TH1D * h_second_shwr_dist_nue_cc_mixed_after   = new TH1D ("h_second_shwr_dist_nue_cc_mixed_after",   "h_second_shwr_dist_nue_cc_mixed_after",   20, 0, 160);
TH1D * h_second_shwr_dist_numu_cc_after        = new TH1D ("h_second_shwr_dist_numu_cc_after",        "h_second_shwr_dist_numu_cc_after",        20, 0, 160);
TH1D * h_second_shwr_dist_numu_cc_mixed_after  = new TH1D ("h_second_shwr_dist_numu_cc_mixed_after",  "h_second_shwr_dist_numu_cc_mixed_after",  20, 0, 160);
TH1D * h_second_shwr_dist_nc_after             = new TH1D ("h_second_shwr_dist_nc_after",             "h_second_shwr_dist_nc_after",             20, 0, 160);
TH1D * h_second_shwr_dist_nc_pi0_after         = new TH1D ("h_second_shwr_dist_nc_pi0_after",         "h_second_shwr_dist_nc_pi0_after",         20, 0, 160);
TH1D * h_second_shwr_dist_cosmic_after         = new TH1D ("h_second_shwr_dist_cosmic_after",         "h_second_shwr_dist_cosmic_after",         20, 0, 160);
TH1D * h_second_shwr_dist_other_mixed_after    = new TH1D ("h_second_shwr_dist_other_mixed_after",    "h_second_shwr_dist_other_mixed_after",    20, 0, 160);
TH1D * h_second_shwr_dist_unmatched_after      = new TH1D ("h_second_shwr_dist_unmatched_after",      "h_second_shwr_dist_unmatched_after",      20, 0, 160);
TH1D * h_second_shwr_dist_intime_after         = new TH1D ("h_second_shwr_dist_intime_after",         "h_second_shwr_dist_intime_after",         20, 0, 160);
TH1D * h_second_shwr_dist_data_after           = new TH1D ("h_second_shwr_dist_data_after",           "h_second_shwr_dist_data_after",           20, 0, 160);

TH1D * h_hit_length_ratio_nue_cc         = new TH1D ("h_hit_length_ratio_nue_cc",         "h_hit_length_ratio_nue_cc",         20, 0, 20);
TH1D * h_hit_length_ratio_nue_cc_out_fv  = new TH1D ("h_hit_length_ratio_nue_cc_out_fv",  "h_hit_length_ratio_nue_cc_out_fv",  20, 0, 20);
TH1D * h_hit_length_ratio_nue_cc_mixed   = new TH1D ("h_hit_length_ratio_nue_cc_mixed",   "h_hit_length_ratio_nue_cc_mixed",   20, 0, 20);
TH1D * h_hit_length_ratio_numu_cc        = new TH1D ("h_hit_length_ratio_numu_cc",        "h_hit_length_ratio_numu_cc",        20, 0, 20);
TH1D * h_hit_length_ratio_numu_cc_mixed  = new TH1D ("h_hit_length_ratio_numu_cc_mixed",  "h_hit_length_ratio_numu_cc_mixed",  20, 0, 20);
TH1D * h_hit_length_ratio_nc             = new TH1D ("h_hit_length_ratio_nc",             "h_hit_length_ratio_nc",             20, 0, 20);
TH1D * h_hit_length_ratio_nc_pi0         = new TH1D ("h_hit_length_ratio_nc_pi0",         "h_hit_length_ratio_nc_pi0",         20, 0, 20);
TH1D * h_hit_length_ratio_cosmic         = new TH1D ("h_hit_length_ratio_cosmic",         "h_hit_length_ratio_cosmic",         20, 0, 20);
TH1D * h_hit_length_ratio_other_mixed    = new TH1D ("h_hit_length_ratio_other_mixed",    "h_hit_length_ratio_other_mixed",    20, 0, 20);
TH1D * h_hit_length_ratio_unmatched      = new TH1D ("h_hit_length_ratio_unmatched",      "h_hit_length_ratio_unmatched",      20, 0, 20);
TH1D * h_hit_length_ratio_intime         = new TH1D ("h_hit_length_ratio_intime",         "h_hit_length_ratio_intime",         20, 0, 20);
TH1D * h_hit_length_ratio_data           = new TH1D ("h_hit_length_ratio_data",           "h_hit_length_ratio_data",           20, 0, 20);

TH1D * h_hit_length_ratio_nue_cc_after         = new TH1D ("h_hit_length_ratio_nue_cc_after",         "h_hit_length_ratio_nue_cc_after",         20, 0, 20);
TH1D * h_hit_length_ratio_nue_cc_out_fv_after  = new TH1D ("h_hit_length_ratio_nue_cc_out_fv_after",  "h_hit_length_ratio_nue_cc_out_fv_after",  20, 0, 20);
TH1D * h_hit_length_ratio_nue_cc_mixed_after   = new TH1D ("h_hit_length_ratio_nue_cc_mixed_after",   "h_hit_length_ratio_nue_cc_mixed_after",   20, 0, 20);
TH1D * h_hit_length_ratio_numu_cc_after        = new TH1D ("h_hit_length_ratio_numu_cc_after",        "h_hit_length_ratio_numu_cc_after",        20, 0, 20);
TH1D * h_hit_length_ratio_numu_cc_mixed_after  = new TH1D ("h_hit_length_ratio_numu_cc_mixed_after",  "h_hit_length_ratio_numu_cc_mixed_after",  20, 0, 20);
TH1D * h_hit_length_ratio_nc_after             = new TH1D ("h_hit_length_ratio_nc_after",             "h_hit_length_ratio_nc_after",             20, 0, 20);
TH1D * h_hit_length_ratio_nc_pi0_after         = new TH1D ("h_hit_length_ratio_nc_pi0_after",         "h_hit_length_ratio_nc_pi0_after",         20, 0, 20);
TH1D * h_hit_length_ratio_cosmic_after         = new TH1D ("h_hit_length_ratio_cosmic_after",         "h_hit_length_ratio_cosmic_after",         20, 0, 20);
TH1D * h_hit_length_ratio_other_mixed_after    = new TH1D ("h_hit_length_ratio_other_mixed_after",    "h_hit_length_ratio_other_mixed_after",    20, 0, 20);
TH1D * h_hit_length_ratio_unmatched_after      = new TH1D ("h_hit_length_ratio_unmatched_after",      "h_hit_length_ratio_unmatched_after",      20, 0, 20);
TH1D * h_hit_length_ratio_intime_after         = new TH1D ("h_hit_length_ratio_intime_after",         "h_hit_length_ratio_intime_after",         20, 0, 20);
TH1D * h_hit_length_ratio_data_after           = new TH1D ("h_hit_length_ratio_data_after",           "h_hit_length_ratio_data_after",           20, 0, 20);

TH1D * h_trk_length_nue_cc         = new TH1D ("h_trk_length_nue_cc",         "h_trk_length_nue_cc",         20, 0, 150);
TH1D * h_trk_length_nue_cc_out_fv  = new TH1D ("h_trk_length_nue_cc_out_fv",  "h_trk_length_nue_cc_out_fv",  20, 0, 150);
TH1D * h_trk_length_nue_cc_mixed   = new TH1D ("h_trk_length_nue_cc_mixed",   "h_trk_length_nue_cc_mixed",   20, 0, 150);
TH1D * h_trk_length_numu_cc        = new TH1D ("h_trk_length_numu_cc",        "h_trk_length_numu_cc",        20, 0, 150);
TH1D * h_trk_length_numu_cc_mixed  = new TH1D ("h_trk_length_numu_cc_mixed",  "h_trk_length_numu_cc_mixed",  20, 0, 150);
TH1D * h_trk_length_nc             = new TH1D ("h_trk_length_nc",             "h_trk_length_nc",             20, 0, 150);
TH1D * h_trk_length_nc_pi0         = new TH1D ("h_trk_length_nc_pi0",         "h_trk_length_nc_pi0",         20, 0, 150);
TH1D * h_trk_length_cosmic         = new TH1D ("h_trk_length_cosmic",         "h_trk_length_cosmic",         20, 0, 150);
TH1D * h_trk_length_other_mixed    = new TH1D ("h_trk_length_other_mixed",    "h_trk_length_other_mixed",    20, 0, 150);
TH1D * h_trk_length_unmatched      = new TH1D ("h_trk_length_unmatched",      "h_trk_length_unmatched",      20, 0, 150);
TH1D * h_trk_length_intime         = new TH1D ("h_trk_length_intime",         "h_trk_length_intime",         20, 0, 150);
TH1D * h_trk_length_data           = new TH1D ("h_trk_length_data",           "h_trk_length_data",           20, 0, 150);

TH1D * h_longest_trk_length_nue_cc         = new TH1D ("h_longest_trk_length_nue_cc",         "h_longest_trk_length_nue_cc",         20, 0, 150);
TH1D * h_longest_trk_length_nue_cc_out_fv  = new TH1D ("h_longest_trk_length_nue_cc_out_fv",  "h_longest_trk_length_nue_cc_out_fv",  20, 0, 150);
TH1D * h_longest_trk_length_nue_cc_mixed   = new TH1D ("h_longest_trk_length_nue_cc_mixed",   "h_longest_trk_length_nue_cc_mixed",   20, 0, 150);
TH1D * h_longest_trk_length_numu_cc        = new TH1D ("h_longest_trk_length_numu_cc",        "h_longest_trk_length_numu_cc",        20, 0, 150);
TH1D * h_longest_trk_length_numu_cc_mixed  = new TH1D ("h_longest_trk_length_numu_cc_mixed",  "h_longest_trk_length_numu_cc_mixed",  20, 0, 150);
TH1D * h_longest_trk_length_nc             = new TH1D ("h_longest_trk_length_nc",             "h_longest_trk_length_nc",             20, 0, 150);
TH1D * h_longest_trk_length_nc_pi0         = new TH1D ("h_longest_trk_length_nc_pi0",         "h_longest_trk_length_nc_pi0",         20, 0, 150);
TH1D * h_longest_trk_length_cosmic         = new TH1D ("h_longest_trk_length_cosmic",         "h_longest_trk_length_cosmic",         20, 0, 150);
TH1D * h_longest_trk_length_other_mixed    = new TH1D ("h_longest_trk_length_other_mixed",    "h_longest_trk_length_other_mixed",    20, 0, 150);
TH1D * h_longest_trk_length_unmatched      = new TH1D ("h_longest_trk_length_unmatched",      "h_longest_trk_length_unmatched",      20, 0, 150);
TH1D * h_longest_trk_length_intime         = new TH1D ("h_longest_trk_length_intime",         "h_longest_trk_length_intime",         20, 0, 150);
TH1D * h_longest_trk_length_data           = new TH1D ("h_longest_trk_length_data",           "h_longest_trk_length_data",           20, 0, 150);

TH1D * h_shwr_length_nue_cc         = new TH1D ("h_shwr_length_nue_cc",         "h_shwr_length_nue_cc",         30, 0, 300);
TH1D * h_shwr_length_nue_cc_out_fv  = new TH1D ("h_shwr_length_nue_cc_out_fv",  "h_shwr_length_nue_cc_out_fv",  30, 0, 300);
TH1D * h_shwr_length_nue_cc_mixed   = new TH1D ("h_shwr_length_nue_cc_mixed",   "h_shwr_length_nue_cc_mixed",   30, 0, 300);
TH1D * h_shwr_length_numu_cc        = new TH1D ("h_shwr_length_numu_cc",        "h_shwr_length_numu_cc",        30, 0, 300);
TH1D * h_shwr_length_numu_cc_mixed  = new TH1D ("h_shwr_length_numu_cc_mixed",  "h_shwr_length_numu_cc_mixed",  30, 0, 300);
TH1D * h_shwr_length_nc             = new TH1D ("h_shwr_length_nc",             "h_shwr_length_nc",             30, 0, 300);
TH1D * h_shwr_length_nc_pi0         = new TH1D ("h_shwr_length_nc_pi0",         "h_shwr_length_nc_pi0",         30, 0, 300);
TH1D * h_shwr_length_cosmic         = new TH1D ("h_shwr_length_cosmic",         "h_shwr_length_cosmic",         30, 0, 300);
TH1D * h_shwr_length_other_mixed    = new TH1D ("h_shwr_length_other_mixed",    "h_shwr_length_other_mixed",    30, 0, 300);
TH1D * h_shwr_length_unmatched      = new TH1D ("h_shwr_length_unmatched",      "h_shwr_length_unmatched",      30, 0, 300);
TH1D * h_shwr_length_intime         = new TH1D ("h_shwr_length_intime",         "h_shwr_length_intime",         30, 0, 300);
TH1D * h_shwr_length_data           = new TH1D ("h_shwr_length_data",           "h_shwr_length_data",           30, 0, 300);

TH1D * h_longest_shwr_length_nue_cc         = new TH1D ("h_longest_shwr_length_nue_cc",         "h_longest_shwr_length_nue_cc",         30, 0, 300);
TH1D * h_longest_shwr_length_nue_cc_out_fv  = new TH1D ("h_longest_shwr_length_nue_cc_out_fv",  "h_longest_shwr_length_nue_cc_out_fv",  30, 0, 300);
TH1D * h_longest_shwr_length_nue_cc_mixed   = new TH1D ("h_longest_shwr_length_nue_cc_mixed",   "h_longest_shwr_length_nue_cc_mixed",   30, 0, 300);
TH1D * h_longest_shwr_length_numu_cc        = new TH1D ("h_longest_shwr_length_numu_cc",        "h_longest_shwr_length_numu_cc",        30, 0, 300);
TH1D * h_longest_shwr_length_numu_cc_mixed  = new TH1D ("h_longest_shwr_length_numu_cc_mixed",  "h_longest_shwr_length_numu_cc_mixed",  30, 0, 300);
TH1D * h_longest_shwr_length_nc             = new TH1D ("h_longest_shwr_length_nc",             "h_longest_shwr_length_nc",             30, 0, 300);
TH1D * h_longest_shwr_length_nc_pi0         = new TH1D ("h_longest_shwr_length_nc_pi0",         "h_longest_shwr_length_nc_pi0",         30, 0, 300);
TH1D * h_longest_shwr_length_cosmic         = new TH1D ("h_longest_shwr_length_cosmic",         "h_longest_shwr_length_cosmic",         30, 0, 300);
TH1D * h_longest_shwr_length_other_mixed    = new TH1D ("h_longest_shwr_length_other_mixed",    "h_longest_shwr_length_other_mixed",    30, 0, 300);
TH1D * h_longest_shwr_length_unmatched      = new TH1D ("h_longest_shwr_length_unmatched",      "h_longest_shwr_length_unmatched",      30, 0, 300);
TH1D * h_longest_shwr_length_intime         = new TH1D ("h_longest_shwr_length_intime",         "h_longest_shwr_length_intime",         30, 0, 300);
TH1D * h_longest_shwr_length_data           = new TH1D ("h_longest_shwr_length_data",           "h_longest_shwr_length_data",           30, 0, 300);

TH1D * h_leading_shwr_length_nue_cc         = new TH1D ("h_leading_shwr_length_nue_cc",         "h_leading_shwr_length_nue_cc",         30, 0, 300);
TH1D * h_leading_shwr_length_nue_cc_out_fv  = new TH1D ("h_leading_shwr_length_nue_cc_out_fv",  "h_leading_shwr_length_nue_cc_out_fv",  30, 0, 300);
TH1D * h_leading_shwr_length_nue_cc_mixed   = new TH1D ("h_leading_shwr_length_nue_cc_mixed",   "h_leading_shwr_length_nue_cc_mixed",   30, 0, 300);
TH1D * h_leading_shwr_length_numu_cc        = new TH1D ("h_leading_shwr_length_numu_cc",        "h_leading_shwr_length_numu_cc",        30, 0, 300);
TH1D * h_leading_shwr_length_numu_cc_mixed  = new TH1D ("h_leading_shwr_length_numu_cc_mixed",  "h_leading_shwr_length_numu_cc_mixed",  30, 0, 300);
TH1D * h_leading_shwr_length_nc             = new TH1D ("h_leading_shwr_length_nc",             "h_leading_shwr_length_nc",             30, 0, 300);
TH1D * h_leading_shwr_length_nc_pi0         = new TH1D ("h_leading_shwr_length_nc_pi0",         "h_leading_shwr_length_nc_pi0",         30, 0, 300);
TH1D * h_leading_shwr_length_cosmic         = new TH1D ("h_leading_shwr_length_cosmic",         "h_leading_shwr_length_cosmic",         30, 0, 300);
TH1D * h_leading_shwr_length_other_mixed    = new TH1D ("h_leading_shwr_length_other_mixed",    "h_leading_shwr_length_other_mixed",    30, 0, 300);
TH1D * h_leading_shwr_length_unmatched      = new TH1D ("h_leading_shwr_length_unmatched",      "h_leading_shwr_length_unmatched",      30, 0, 300);
TH1D * h_leading_shwr_length_intime         = new TH1D ("h_leading_shwr_length_intime",         "h_leading_shwr_length_intime",         30, 0, 300);
TH1D * h_leading_shwr_length_data           = new TH1D ("h_leading_shwr_length_data",           "h_leading_shwr_length_data",           30, 0, 300);

TH1D * h_leading_shwr_length_nue_cc_after         = new TH1D ("h_leading_shwr_length_nue_cc_after",         "h_leading_shwr_length_nue_cc_after",         30, 0, 300);
TH1D * h_leading_shwr_length_nue_cc_out_fv_after  = new TH1D ("h_leading_shwr_length_nue_cc_out_fv_after",  "h_leading_shwr_length_nue_cc_out_fv_after",  30, 0, 300);
TH1D * h_leading_shwr_length_nue_cc_mixed_after   = new TH1D ("h_leading_shwr_length_nue_cc_mixed_after",   "h_leading_shwr_length_nue_cc_mixed_after",   30, 0, 300);
TH1D * h_leading_shwr_length_numu_cc_after        = new TH1D ("h_leading_shwr_length_numu_cc_after",        "h_leading_shwr_length_numu_cc_after",        30, 0, 300);
TH1D * h_leading_shwr_length_numu_cc_mixed_after  = new TH1D ("h_leading_shwr_length_numu_cc_mixed_after",  "h_leading_shwr_length_numu_cc_mixed_after",  30, 0, 300);
TH1D * h_leading_shwr_length_nc_after             = new TH1D ("h_leading_shwr_length_nc_after",             "h_leading_shwr_length_nc_after",             30, 0, 300);
TH1D * h_leading_shwr_length_nc_pi0_after         = new TH1D ("h_leading_shwr_length_nc_pi0_after",         "h_leading_shwr_length_nc_pi0_after",         30, 0, 300);
TH1D * h_leading_shwr_length_cosmic_after         = new TH1D ("h_leading_shwr_length_cosmic_after",         "h_leading_shwr_length_cosmic_after",         30, 0, 300);
TH1D * h_leading_shwr_length_other_mixed_after    = new TH1D ("h_leading_shwr_length_other_mixed_after",    "h_leading_shwr_length_other_mixed_after",    30, 0, 300);
TH1D * h_leading_shwr_length_unmatched_after      = new TH1D ("h_leading_shwr_length_unmatched_after",      "h_leading_shwr_length_unmatched_after",      30, 0, 300);
TH1D * h_leading_shwr_length_intime_after         = new TH1D ("h_leading_shwr_length_intime_after",         "h_leading_shwr_length_intime_after",         30, 0, 300);
TH1D * h_leading_shwr_length_data_after           = new TH1D ("h_leading_shwr_length_data_after",           "h_leading_shwr_length_data_after",           30, 0, 300);

TH1D * h_leading_shwr_trk_length_nue_cc         = new TH1D ("h_leading_shwr_trk_length_nue_cc",         "h_leading_shwr_trk_length_nue_cc",         20, 0, 3);
TH1D * h_leading_shwr_trk_length_nue_cc_out_fv  = new TH1D ("h_leading_shwr_trk_length_nue_cc_out_fv",  "h_leading_shwr_trk_length_nue_cc_out_fv",  20, 0, 3);
TH1D * h_leading_shwr_trk_length_nue_cc_mixed   = new TH1D ("h_leading_shwr_trk_length_nue_cc_mixed",   "h_leading_shwr_trk_length_nue_cc_mixed",   20, 0, 3);
TH1D * h_leading_shwr_trk_length_numu_cc        = new TH1D ("h_leading_shwr_trk_length_numu_cc",        "h_leading_shwr_trk_length_numu_cc",        20, 0, 3);
TH1D * h_leading_shwr_trk_length_numu_cc_mixed  = new TH1D ("h_leading_shwr_trk_length_numu_cc_mixed",  "h_leading_shwr_trk_length_numu_cc_mixed",  20, 0, 3);
TH1D * h_leading_shwr_trk_length_nc             = new TH1D ("h_leading_shwr_trk_length_nc",             "h_leading_shwr_trk_length_nc",             20, 0, 3);
TH1D * h_leading_shwr_trk_length_nc_pi0         = new TH1D ("h_leading_shwr_trk_length_nc_pi0",         "h_leading_shwr_trk_length_nc_pi0",         20, 0, 3);
TH1D * h_leading_shwr_trk_length_cosmic         = new TH1D ("h_leading_shwr_trk_length_cosmic",         "h_leading_shwr_trk_length_cosmic",         20, 0, 3);
TH1D * h_leading_shwr_trk_length_other_mixed    = new TH1D ("h_leading_shwr_trk_length_other_mixed",    "h_leading_shwr_trk_length_other_mixed",    20, 0, 3);
TH1D * h_leading_shwr_trk_length_unmatched      = new TH1D ("h_leading_shwr_trk_length_unmatched",      "h_leading_shwr_trk_length_unmatched",      20, 0, 3);
TH1D * h_leading_shwr_trk_length_intime         = new TH1D ("h_leading_shwr_trk_length_intime",         "h_leading_shwr_trk_length_intime",         20, 0, 3);
TH1D * h_leading_shwr_trk_length_data           = new TH1D ("h_leading_shwr_trk_length_data",           "h_leading_shwr_trk_length_data",           20, 0, 3);

TH1D * h_leading_shwr_trk_length_nue_cc_after         = new TH1D ("h_leading_shwr_trk_length_nue_cc_after",         "h_leading_shwr_trk_length_nue_cc_after",         20, 0, 3);
TH1D * h_leading_shwr_trk_length_nue_cc_out_fv_after  = new TH1D ("h_leading_shwr_trk_length_nue_cc_out_fv_after",  "h_leading_shwr_trk_length_nue_cc_out_fv_after",  20, 0, 3);
TH1D * h_leading_shwr_trk_length_nue_cc_mixed_after   = new TH1D ("h_leading_shwr_trk_length_nue_cc_mixed_after",   "h_leading_shwr_trk_length_nue_cc_mixed_after",   20, 0, 3);
TH1D * h_leading_shwr_trk_length_numu_cc_after        = new TH1D ("h_leading_shwr_trk_length_numu_cc_after",        "h_leading_shwr_trk_length_numu_cc_after",        20, 0, 3);
TH1D * h_leading_shwr_trk_length_numu_cc_mixed_after  = new TH1D ("h_leading_shwr_trk_length_numu_cc_mixed_after",  "h_leading_shwr_trk_length_numu_cc_mixed_after",  20, 0, 3);
TH1D * h_leading_shwr_trk_length_nc_after             = new TH1D ("h_leading_shwr_trk_length_nc_after",             "h_leading_shwr_trk_length_nc_after",             20, 0, 3);
TH1D * h_leading_shwr_trk_length_nc_pi0_after         = new TH1D ("h_leading_shwr_trk_length_nc_pi0_after",         "h_leading_shwr_trk_length_nc_pi0_after",         20, 0, 3);
TH1D * h_leading_shwr_trk_length_cosmic_after         = new TH1D ("h_leading_shwr_trk_length_cosmic_after",         "h_leading_shwr_trk_length_cosmic_after",         20, 0, 3);
TH1D * h_leading_shwr_trk_length_other_mixed_after    = new TH1D ("h_leading_shwr_trk_length_other_mixed_after",    "h_leading_shwr_trk_length_other_mixed_after",    20, 0, 3);
TH1D * h_leading_shwr_trk_length_unmatched_after      = new TH1D ("h_leading_shwr_trk_length_unmatched_after",      "h_leading_shwr_trk_length_unmatched_after",      20, 0, 3);
TH1D * h_leading_shwr_trk_length_intime_after         = new TH1D ("h_leading_shwr_trk_length_intime_after",         "h_leading_shwr_trk_length_intime_after",         20, 0, 3);
TH1D * h_leading_shwr_trk_length_data_after           = new TH1D ("h_leading_shwr_trk_length_data_after",           "h_leading_shwr_trk_length_data_after",           20, 0, 3);

TH1D * h_longest_shwr_trk_length_nue_cc         = new TH1D ("h_longest_shwr_trk_length_nue_cc",         "h_longest_shwr_trk_length_nue_cc",         20, 0, 3);
TH1D * h_longest_shwr_trk_length_nue_cc_out_fv  = new TH1D ("h_longest_shwr_trk_length_nue_cc_out_fv",  "h_longest_shwr_trk_length_nue_cc_out_fv",  20, 0, 3);
TH1D * h_longest_shwr_trk_length_nue_cc_mixed   = new TH1D ("h_longest_shwr_trk_length_nue_cc_mixed",   "h_longest_shwr_trk_length_nue_cc_mixed",   20, 0, 3);
TH1D * h_longest_shwr_trk_length_numu_cc        = new TH1D ("h_longest_shwr_trk_length_numu_cc",        "h_longest_shwr_trk_length_numu_cc",        20, 0, 3);
TH1D * h_longest_shwr_trk_length_numu_cc_mixed  = new TH1D ("h_longest_shwr_trk_length_numu_cc_mixed",  "h_longest_shwr_trk_length_numu_cc_mixed",  20, 0, 3);
TH1D * h_longest_shwr_trk_length_nc             = new TH1D ("h_longest_shwr_trk_length_nc",             "h_longest_shwr_trk_length_nc",             20, 0, 3);
TH1D * h_longest_shwr_trk_length_nc_pi0         = new TH1D ("h_longest_shwr_trk_length_nc_pi0",         "h_longest_shwr_trk_length_nc_pi0",         20, 0, 3);
TH1D * h_longest_shwr_trk_length_cosmic         = new TH1D ("h_longest_shwr_trk_length_cosmic",         "h_longest_shwr_trk_length_cosmic",         20, 0, 3);
TH1D * h_longest_shwr_trk_length_other_mixed    = new TH1D ("h_longest_shwr_trk_length_other_mixed",    "h_longest_shwr_trk_length_other_mixed",    20, 0, 3);
TH1D * h_longest_shwr_trk_length_unmatched      = new TH1D ("h_longest_shwr_trk_length_unmatched",      "h_longest_shwr_trk_length_unmatched",      20, 0, 3);
TH1D * h_longest_shwr_trk_length_intime         = new TH1D ("h_longest_shwr_trk_length_intime",         "h_longest_shwr_trk_length_intime",         20, 0, 3);
TH1D * h_longest_shwr_trk_length_data           = new TH1D ("h_longest_shwr_trk_length_data",           "h_longest_shwr_trk_length_data",           20, 0, 3);

TH2D * h_collection_total_hits_track_nue_cc         = new TH2D ("h_collection_total_hits_track_nue_cc",         "h_collection_total_hits_track_nue_cc",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_nue_cc_out_fv  = new TH2D ("h_collection_total_hits_track_nue_cc_out_fv",  "h_collection_total_hits_track_nue_cc_out_fv",  20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_nue_cc_mixed   = new TH2D ("h_collection_total_hits_track_nue_cc_mixed",   "h_collection_total_hits_track_nue_cc_mixed",   20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_numu_cc        = new TH2D ("h_collection_total_hits_track_numu_cc",        "h_collection_total_hits_track_numu_cc",        20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_numu_cc_mixed  = new TH2D ("h_collection_total_hits_track_numu_cc_mixed",  "h_collection_total_hits_track_numu_cc_mixed",  20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_nc             = new TH2D ("h_collection_total_hits_track_nc",             "h_collection_total_hits_track_nc",             20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_nc_pi0         = new TH2D ("h_collection_total_hits_track_nc_pi0",         "h_collection_total_hits_track_nc_pi0",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_cosmic         = new TH2D ("h_collection_total_hits_track_cosmic",         "h_collection_total_hits_track_cosmic",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_other_mixed    = new TH2D ("h_collection_total_hits_track_other_mixed",    "h_collection_total_hits_track_other_mixed",    20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_unmatched      = new TH2D ("h_collection_total_hits_track_unmatched",      "h_collection_total_hits_track_unmatched",      20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_intime         = new TH2D ("h_collection_total_hits_track_intime",         "h_collection_total_hits_track_intime",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_track_data           = new TH2D ("h_collection_total_hits_track_data",           "h_collection_total_hits_track_data",           20, 0, 200, 20, 0, 400);

TH2D * h_collection_total_hits_shower_nue_cc         = new TH2D ("h_collection_total_hits_shower_nue_cc",         "h_collection_total_hits_shower_nue_cc",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_nue_cc_out_fv  = new TH2D ("h_collection_total_hits_shower_nue_cc_out_fv",  "h_collection_total_hits_shower_nue_cc_out_fv",  20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_nue_cc_mixed   = new TH2D ("h_collection_total_hits_shower_nue_cc_mixed",   "h_collection_total_hits_shower_nue_cc_mixed",   20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_numu_cc        = new TH2D ("h_collection_total_hits_shower_numu_cc",        "h_collection_total_hits_shower_numu_cc",        20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_numu_cc_mixed  = new TH2D ("h_collection_total_hits_shower_numu_cc_mixed",  "h_collection_total_hits_shower_numu_cc_mixed",  20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_nc             = new TH2D ("h_collection_total_hits_shower_nc",             "h_collection_total_hits_shower_nc",             20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_nc_pi0         = new TH2D ("h_collection_total_hits_shower_nc_pi0",         "h_collection_total_hits_shower_nc_pi0",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_cosmic         = new TH2D ("h_collection_total_hits_shower_cosmic",         "h_collection_total_hits_shower_cosmic",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_other_mixed    = new TH2D ("h_collection_total_hits_shower_other_mixed",    "h_collection_total_hits_shower_other_mixed",    20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_unmatched      = new TH2D ("h_collection_total_hits_shower_unmatched",      "h_collection_total_hits_shower_unmatched",      20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_intime         = new TH2D ("h_collection_total_hits_shower_intime",         "h_collection_total_hits_shower_intime",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_shower_data           = new TH2D ("h_collection_total_hits_shower_data",           "h_collection_total_hits_shower_data",           20, 0, 200, 20, 0, 400);

TH2D * h_collection_total_hits_leading_shower_nue_cc         = new TH2D ("h_collection_total_hits_leading_shower_nue_cc",         "h_collection_total_hits_leading_shower_nue_cc",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_nue_cc_out_fv  = new TH2D ("h_collection_total_hits_leading_shower_nue_cc_out_fv",  "h_collection_total_hits_leading_shower_nue_cc_out_fv",  20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_nue_cc_mixed   = new TH2D ("h_collection_total_hits_leading_shower_nue_cc_mixed",   "h_collection_total_hits_leading_shower_nue_cc_mixed",   20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_numu_cc        = new TH2D ("h_collection_total_hits_leading_shower_numu_cc",        "h_collection_total_hits_leading_shower_numu_cc",        20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_numu_cc_mixed  = new TH2D ("h_collection_total_hits_leading_shower_numu_cc_mixed",  "h_collection_total_hits_leading_shower_numu_cc_mixed",  20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_nc             = new TH2D ("h_collection_total_hits_leading_shower_nc",             "h_collection_total_hits_leading_shower_nc",             20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_nc_pi0         = new TH2D ("h_collection_total_hits_leading_shower_nc_pi0",         "h_collection_total_hits_leading_shower_nc_pi0",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_cosmic         = new TH2D ("h_collection_total_hits_leading_shower_cosmic",         "h_collection_total_hits_leading_shower_cosmic",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_other_mixed    = new TH2D ("h_collection_total_hits_leading_shower_other_mixed",    "h_collection_total_hits_leading_shower_other_mixed",    20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_unmatched      = new TH2D ("h_collection_total_hits_leading_shower_unmatched",      "h_collection_total_hits_leading_shower_unmatched",      20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_intime         = new TH2D ("h_collection_total_hits_leading_shower_intime",         "h_collection_total_hits_leading_shower_intime",         20, 0, 200, 20, 0, 400);
TH2D * h_collection_total_hits_leading_shower_data           = new TH2D ("h_collection_total_hits_leading_shower_data",           "h_collection_total_hits_leading_shower_data",           20, 0, 200, 20, 0, 400);

TH1D * h_collection_hits_track_nue_cc         = new TH1D ("h_collection_hits_track_nue_cc",         "h_collection_hits_track_nue_cc",         20, 0, 250);
TH1D * h_collection_hits_track_nue_cc_out_fv  = new TH1D ("h_collection_hits_track_nue_cc_out_fv",  "h_collection_hits_track_nue_cc_out_fv",  20, 0, 250);
TH1D * h_collection_hits_track_nue_cc_mixed   = new TH1D ("h_collection_hits_track_nue_cc_mixed",   "h_collection_hits_track_nue_cc_mixed",   20, 0, 250);
TH1D * h_collection_hits_track_numu_cc        = new TH1D ("h_collection_hits_track_numu_cc",        "h_collection_hits_track_numu_cc",        20, 0, 250);
TH1D * h_collection_hits_track_numu_cc_mixed  = new TH1D ("h_collection_hits_track_numu_cc_mixed",  "h_collection_hits_track_numu_cc_mixed",  20, 0, 250);
TH1D * h_collection_hits_track_nc             = new TH1D ("h_collection_hits_track_nc",             "h_collection_hits_track_nc",             20, 0, 250);
TH1D * h_collection_hits_track_nc_pi0         = new TH1D ("h_collection_hits_track_nc_pi0",         "h_collection_hits_track_nc_pi0",         20, 0, 250);
TH1D * h_collection_hits_track_cosmic         = new TH1D ("h_collection_hits_track_cosmic",         "h_collection_hits_track_cosmic",         20, 0, 250);
TH1D * h_collection_hits_track_other_mixed    = new TH1D ("h_collection_hits_track_other_mixed",    "h_collection_hits_track_other_mixed",    20, 0, 250);
TH1D * h_collection_hits_track_unmatched      = new TH1D ("h_collection_hits_track_unmatched",      "h_collection_hits_track_unmatched",      20, 0, 250);
TH1D * h_collection_hits_track_intime         = new TH1D ("h_collection_hits_track_intime",         "h_collection_hits_track_intime",         20, 0, 250);
TH1D * h_collection_hits_track_data           = new TH1D ("h_collection_hits_track_data",           "h_collection_hits_track_data",           20, 0, 250);

TH1D * h_collection_hits_track_nue_cc_after         = new TH1D ("h_collection_hits_track_nue_cc_after",         "h_collection_hits_track_nue_cc_after",         20, 0, 250);
TH1D * h_collection_hits_track_nue_cc_out_fv_after  = new TH1D ("h_collection_hits_track_nue_cc_out_fv_after",  "h_collection_hits_track_nue_cc_out_fv_after",  20, 0, 250);
TH1D * h_collection_hits_track_nue_cc_mixed_after   = new TH1D ("h_collection_hits_track_nue_cc_mixed_after",   "h_collection_hits_track_nue_cc_mixed_after",   20, 0, 250);
TH1D * h_collection_hits_track_numu_cc_after        = new TH1D ("h_collection_hits_track_numu_cc_after",        "h_collection_hits_track_numu_cc_after",        20, 0, 250);
TH1D * h_collection_hits_track_numu_cc_mixed_after  = new TH1D ("h_collection_hits_track_numu_cc_mixed_after",  "h_collection_hits_track_numu_cc_mixed_after",  20, 0, 250);
TH1D * h_collection_hits_track_nc_after             = new TH1D ("h_collection_hits_track_nc_after",             "h_collection_hits_track_nc_after",             20, 0, 250);
TH1D * h_collection_hits_track_nc_pi0_after         = new TH1D ("h_collection_hits_track_nc_pi0_after",         "h_collection_hits_track_nc_pi0_after",         20, 0, 250);
TH1D * h_collection_hits_track_cosmic_after         = new TH1D ("h_collection_hits_track_cosmic_after",         "h_collection_hits_track_cosmic_after",         20, 0, 250);
TH1D * h_collection_hits_track_other_mixed_after    = new TH1D ("h_collection_hits_track_other_mixed_after",    "h_collection_hits_track_other_mixed_after",    20, 0, 250);
TH1D * h_collection_hits_track_unmatched_after      = new TH1D ("h_collection_hits_track_unmatched_after",      "h_collection_hits_track_unmatched_after",      20, 0, 250);
TH1D * h_collection_hits_track_intime_after         = new TH1D ("h_collection_hits_track_intime_after",         "h_collection_hits_track_intime_after",         20, 0, 250);
TH1D * h_collection_hits_track_data_after           = new TH1D ("h_collection_hits_track_data_after",           "h_collection_hits_track_data_after",           20, 0, 250);

TH1D * h_collection_hits_shower_nue_cc         = new TH1D ("h_collection_hits_shower_nue_cc",         "h_collection_hits_shower_nue_cc",         20, 0, 250);
TH1D * h_collection_hits_shower_nue_cc_out_fv  = new TH1D ("h_collection_hits_shower_nue_cc_out_fv",  "h_collection_hits_shower_nue_cc_out_fv",  20, 0, 250);
TH1D * h_collection_hits_shower_nue_cc_mixed   = new TH1D ("h_collection_hits_shower_nue_cc_mixed",   "h_collection_hits_shower_nue_cc_mixed",   20, 0, 250);
TH1D * h_collection_hits_shower_numu_cc        = new TH1D ("h_collection_hits_shower_numu_cc",        "h_collection_hits_shower_numu_cc",        20, 0, 250);
TH1D * h_collection_hits_shower_numu_cc_mixed  = new TH1D ("h_collection_hits_shower_numu_cc_mixed",  "h_collection_hits_shower_numu_cc_mixed",  20, 0, 250);
TH1D * h_collection_hits_shower_nc             = new TH1D ("h_collection_hits_shower_nc",             "h_collection_hits_shower_nc",             20, 0, 250);
TH1D * h_collection_hits_shower_nc_pi0         = new TH1D ("h_collection_hits_shower_nc_pi0",         "h_collection_hits_shower_nc_pi0",         20, 0, 250);
TH1D * h_collection_hits_shower_cosmic         = new TH1D ("h_collection_hits_shower_cosmic",         "h_collection_hits_shower_cosmic",         20, 0, 250);
TH1D * h_collection_hits_shower_other_mixed    = new TH1D ("h_collection_hits_shower_other_mixed",    "h_collection_hits_shower_other_mixed",    20, 0, 250);
TH1D * h_collection_hits_shower_unmatched      = new TH1D ("h_collection_hits_shower_unmatched",      "h_collection_hits_shower_unmatched",      20, 0, 250);
TH1D * h_collection_hits_shower_intime         = new TH1D ("h_collection_hits_shower_intime",         "h_collection_hits_shower_intime",         20, 0, 250);
TH1D * h_collection_hits_shower_data           = new TH1D ("h_collection_hits_shower_data",           "h_collection_hits_shower_data",           20, 0, 250);

TH1D * h_collection_hits_shower_nue_cc_after         = new TH1D ("h_collection_hits_shower_nue_cc_after",         "h_collection_hits_shower_nue_cc_after",         20, 0, 250);
TH1D * h_collection_hits_shower_nue_cc_out_fv_after  = new TH1D ("h_collection_hits_shower_nue_cc_out_fv_after",  "h_collection_hits_shower_nue_cc_out_fv_after",  20, 0, 250);
TH1D * h_collection_hits_shower_nue_cc_mixed_after   = new TH1D ("h_collection_hits_shower_nue_cc_mixed_after",   "h_collection_hits_shower_nue_cc_mixed_after",   20, 0, 250);
TH1D * h_collection_hits_shower_numu_cc_after        = new TH1D ("h_collection_hits_shower_numu_cc_after",        "h_collection_hits_shower_numu_cc_after",        20, 0, 250);
TH1D * h_collection_hits_shower_numu_cc_mixed_after  = new TH1D ("h_collection_hits_shower_numu_cc_mixed_after",  "h_collection_hits_shower_numu_cc_mixed_after",  20, 0, 250);
TH1D * h_collection_hits_shower_nc_after             = new TH1D ("h_collection_hits_shower_nc_after",             "h_collection_hits_shower_nc_after",             20, 0, 250);
TH1D * h_collection_hits_shower_nc_pi0_after         = new TH1D ("h_collection_hits_shower_nc_pi0_after",         "h_collection_hits_shower_nc_pi0_after",         20, 0, 250);
TH1D * h_collection_hits_shower_cosmic_after         = new TH1D ("h_collection_hits_shower_cosmic_after",         "h_collection_hits_shower_cosmic_after",         20, 0, 250);
TH1D * h_collection_hits_shower_other_mixed_after    = new TH1D ("h_collection_hits_shower_other_mixed_after",    "h_collection_hits_shower_other_mixed_after",    20, 0, 250);
TH1D * h_collection_hits_shower_unmatched_after      = new TH1D ("h_collection_hits_shower_unmatched_after",      "h_collection_hits_shower_unmatched_after",      20, 0, 250);
TH1D * h_collection_hits_shower_intime_after         = new TH1D ("h_collection_hits_shower_intime_after",         "h_collection_hits_shower_intime_after",         20, 0, 250);
TH1D * h_collection_hits_shower_data_after           = new TH1D ("h_collection_hits_shower_data_after",           "h_collection_hits_shower_data_after",           20, 0, 250);

TH1D * h_collection_hits_leading_shower_nue_cc         = new TH1D ("h_collection_hits_leading_shower_nue_cc",         "h_collection_hits_leading_shower_nue_cc",         30, 0, 300);
TH1D * h_collection_hits_leading_shower_nue_cc_out_fv  = new TH1D ("h_collection_hits_leading_shower_nue_cc_out_fv",  "h_collection_hits_leading_shower_nue_cc_out_fv",  30, 0, 300);
TH1D * h_collection_hits_leading_shower_nue_cc_mixed   = new TH1D ("h_collection_hits_leading_shower_nue_cc_mixed",   "h_collection_hits_leading_shower_nue_cc_mixed",   30, 0, 300);
TH1D * h_collection_hits_leading_shower_numu_cc        = new TH1D ("h_collection_hits_leading_shower_numu_cc",        "h_collection_hits_leading_shower_numu_cc",        30, 0, 300);
TH1D * h_collection_hits_leading_shower_numu_cc_mixed  = new TH1D ("h_collection_hits_leading_shower_numu_cc_mixed",  "h_collection_hits_leading_shower_numu_cc_mixed",  30, 0, 300);
TH1D * h_collection_hits_leading_shower_nc             = new TH1D ("h_collection_hits_leading_shower_nc",             "h_collection_hits_leading_shower_nc",             30, 0, 300);
TH1D * h_collection_hits_leading_shower_nc_pi0         = new TH1D ("h_collection_hits_leading_shower_nc_pi0",         "h_collection_hits_leading_shower_nc_pi0",         30, 0, 300);
TH1D * h_collection_hits_leading_shower_cosmic         = new TH1D ("h_collection_hits_leading_shower_cosmic",         "h_collection_hits_leading_shower_cosmic",         30, 0, 300);
TH1D * h_collection_hits_leading_shower_other_mixed    = new TH1D ("h_collection_hits_leading_shower_other_mixed",    "h_collection_hits_leading_shower_other_mixed",    30, 0, 300);
TH1D * h_collection_hits_leading_shower_unmatched      = new TH1D ("h_collection_hits_leading_shower_unmatched",      "h_collection_hits_leading_shower_unmatched",      30, 0, 300);
TH1D * h_collection_hits_leading_shower_intime         = new TH1D ("h_collection_hits_leading_shower_intime",         "h_collection_hits_leading_shower_intime",         30, 0, 300);
TH1D * h_collection_hits_leading_shower_data           = new TH1D ("h_collection_hits_leading_shower_data",           "h_collection_hits_leading_shower_data",           30, 0, 300);

TH1D * h_collection_hits_leading_shower_nue_cc_after         = new TH1D ("h_collection_hits_leading_shower_nue_cc_after",         "h_collection_hits_leading_shower_nue_cc_after",         30, 0, 300);
TH1D * h_collection_hits_leading_shower_nue_cc_out_fv_after  = new TH1D ("h_collection_hits_leading_shower_nue_cc_out_fv_after",  "h_collection_hits_leading_shower_nue_cc_out_fv_after",  30, 0, 300);
TH1D * h_collection_hits_leading_shower_nue_cc_mixed_after   = new TH1D ("h_collection_hits_leading_shower_nue_cc_mixed_after",   "h_collection_hits_leading_shower_nue_cc_mixed_after",   30, 0, 300);
TH1D * h_collection_hits_leading_shower_numu_cc_after        = new TH1D ("h_collection_hits_leading_shower_numu_cc_after",        "h_collection_hits_leading_shower_numu_cc_after",        30, 0, 300);
TH1D * h_collection_hits_leading_shower_numu_cc_mixed_after  = new TH1D ("h_collection_hits_leading_shower_numu_cc_mixed_after",  "h_collection_hits_leading_shower_numu_cc_mixed_after",  30, 0, 300);
TH1D * h_collection_hits_leading_shower_nc_after             = new TH1D ("h_collection_hits_leading_shower_nc_after",             "h_collection_hits_leading_shower_nc_after",             30, 0, 300);
TH1D * h_collection_hits_leading_shower_nc_pi0_after         = new TH1D ("h_collection_hits_leading_shower_nc_pi0_after",         "h_collection_hits_leading_shower_nc_pi0_after",         30, 0, 300);
TH1D * h_collection_hits_leading_shower_cosmic_after         = new TH1D ("h_collection_hits_leading_shower_cosmic_after",         "h_collection_hits_leading_shower_cosmic_after",         30, 0, 300);
TH1D * h_collection_hits_leading_shower_other_mixed_after    = new TH1D ("h_collection_hits_leading_shower_other_mixed_after",    "h_collection_hits_leading_shower_other_mixed_after",    30, 0, 300);
TH1D * h_collection_hits_leading_shower_unmatched_after      = new TH1D ("h_collection_hits_leading_shower_unmatched_after",      "h_collection_hits_leading_shower_unmatched_after",      30, 0, 300);
TH1D * h_collection_hits_leading_shower_intime_after         = new TH1D ("h_collection_hits_leading_shower_intime_after",         "h_collection_hits_leading_shower_intime_after",         30, 0, 300);
TH1D * h_collection_hits_leading_shower_data_after           = new TH1D ("h_collection_hits_leading_shower_data_after",           "h_collection_hits_leading_shower_data_after",           30, 0, 300);

TH1D * h_total_hits_leading_shower_nue_cc         = new TH1D ("h_total_hits_leading_shower_nue_cc",         "h_total_hits_leading_shower_nue_cc",         30, 0, 600);
TH1D * h_total_hits_leading_shower_nue_cc_out_fv  = new TH1D ("h_total_hits_leading_shower_nue_cc_out_fv",  "h_total_hits_leading_shower_nue_cc_out_fv",  30, 0, 600);
TH1D * h_total_hits_leading_shower_nue_cc_mixed   = new TH1D ("h_total_hits_leading_shower_nue_cc_mixed",   "h_total_hits_leading_shower_nue_cc_mixed",   30, 0, 600);
TH1D * h_total_hits_leading_shower_numu_cc        = new TH1D ("h_total_hits_leading_shower_numu_cc",        "h_total_hits_leading_shower_numu_cc",        30, 0, 600);
TH1D * h_total_hits_leading_shower_numu_cc_mixed  = new TH1D ("h_total_hits_leading_shower_numu_cc_mixed",  "h_total_hits_leading_shower_numu_cc_mixed",  30, 0, 600);
TH1D * h_total_hits_leading_shower_nc             = new TH1D ("h_total_hits_leading_shower_nc",             "h_total_hits_leading_shower_nc",             30, 0, 600);
TH1D * h_total_hits_leading_shower_nc_pi0         = new TH1D ("h_total_hits_leading_shower_nc_pi0",         "h_total_hits_leading_shower_nc_pi0",         30, 0, 600);
TH1D * h_total_hits_leading_shower_cosmic         = new TH1D ("h_total_hits_leading_shower_cosmic",         "h_total_hits_leading_shower_cosmic",         30, 0, 600);
TH1D * h_total_hits_leading_shower_other_mixed    = new TH1D ("h_total_hits_leading_shower_other_mixed",    "h_total_hits_leading_shower_other_mixed",    30, 0, 600);
TH1D * h_total_hits_leading_shower_unmatched      = new TH1D ("h_total_hits_leading_shower_unmatched",      "h_total_hits_leading_shower_unmatched",      30, 0, 600);
TH1D * h_total_hits_leading_shower_intime         = new TH1D ("h_total_hits_leading_shower_intime",         "h_total_hits_leading_shower_intime",         30, 0, 600);
TH1D * h_total_hits_leading_shower_data           = new TH1D ("h_total_hits_leading_shower_data",           "h_total_hits_leading_shower_data",           30, 0, 600);

TH1D * h_total_hits_leading_shower_nue_cc_after         = new TH1D ("h_total_hits_leading_shower_nue_cc_after",         "h_total_hits_leading_shower_nue_cc_after",         30, 0, 600);
TH1D * h_total_hits_leading_shower_nue_cc_out_fv_after  = new TH1D ("h_total_hits_leading_shower_nue_cc_out_fv_after",  "h_total_hits_leading_shower_nue_cc_out_fv_after",  30, 0, 600);
TH1D * h_total_hits_leading_shower_nue_cc_mixed_after   = new TH1D ("h_total_hits_leading_shower_nue_cc_mixed_after",   "h_total_hits_leading_shower_nue_cc_mixed_after",   30, 0, 600);
TH1D * h_total_hits_leading_shower_numu_cc_after        = new TH1D ("h_total_hits_leading_shower_numu_cc_after",        "h_total_hits_leading_shower_numu_cc_after",        30, 0, 600);
TH1D * h_total_hits_leading_shower_numu_cc_mixed_after  = new TH1D ("h_total_hits_leading_shower_numu_cc_mixed_after",  "h_total_hits_leading_shower_numu_cc_mixed_after",  30, 0, 600);
TH1D * h_total_hits_leading_shower_nc_after             = new TH1D ("h_total_hits_leading_shower_nc_after",             "h_total_hits_leading_shower_nc_after",             30, 0, 600);
TH1D * h_total_hits_leading_shower_nc_pi0_after         = new TH1D ("h_total_hits_leading_shower_nc_pi0_after",         "h_total_hits_leading_shower_nc_pi0_after",         30, 0, 600);
TH1D * h_total_hits_leading_shower_cosmic_after         = new TH1D ("h_total_hits_leading_shower_cosmic_after",         "h_total_hits_leading_shower_cosmic_after",         30, 0, 600);
TH1D * h_total_hits_leading_shower_other_mixed_after    = new TH1D ("h_total_hits_leading_shower_other_mixed_after",    "h_total_hits_leading_shower_other_mixed_after",    30, 0, 600);
TH1D * h_total_hits_leading_shower_unmatched_after      = new TH1D ("h_total_hits_leading_shower_unmatched_after",      "h_total_hits_leading_shower_unmatched_after",      30, 0, 600);
TH1D * h_total_hits_leading_shower_intime_after         = new TH1D ("h_total_hits_leading_shower_intime_after",         "h_total_hits_leading_shower_intime_after",         30, 0, 600);
TH1D * h_total_hits_leading_shower_data_after           = new TH1D ("h_total_hits_leading_shower_data_after",           "h_total_hits_leading_shower_data_after",           30, 0, 600);

TH1D * h_pre_cut_collection_hits_track_nue_cc         = new TH1D ("h_pre_cut_collection_hits_track_nue_cc",         "h_pre_cut_collection_hits_track_nue_cc",         20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_nue_cc_out_fv  = new TH1D ("h_pre_cut_collection_hits_track_nue_cc_out_fv",  "h_pre_cut_collection_hits_track_nue_cc_out_fv",  20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_nue_cc_mixed   = new TH1D ("h_pre_cut_collection_hits_track_nue_cc_mixed",   "h_pre_cut_collection_hits_track_nue_cc_mixed",   20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_numu_cc        = new TH1D ("h_pre_cut_collection_hits_track_numu_cc",        "h_pre_cut_collection_hits_track_numu_cc",        20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_numu_cc_mixed  = new TH1D ("h_pre_cut_collection_hits_track_numu_cc_mixed",  "h_pre_cut_collection_hits_track_numu_cc_mixed",  20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_nc             = new TH1D ("h_pre_cut_collection_hits_track_nc",             "h_pre_cut_collection_hits_track_nc",             20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_nc_pi0         = new TH1D ("h_pre_cut_collection_hits_track_nc_pi0",         "h_pre_cut_collection_hits_track_nc_pi0",         20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_cosmic         = new TH1D ("h_pre_cut_collection_hits_track_cosmic",         "h_pre_cut_collection_hits_track_cosmic",         20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_other_mixed    = new TH1D ("h_pre_cut_collection_hits_track_other_mixed",    "h_pre_cut_collection_hits_track_other_mixed",    20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_unmatched      = new TH1D ("h_pre_cut_collection_hits_track_unmatched",      "h_pre_cut_collection_hits_track_unmatched",      20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_intime         = new TH1D ("h_pre_cut_collection_hits_track_intime",         "h_pre_cut_collection_hits_track_intime",         20, 0, 250);
TH1D * h_pre_cut_collection_hits_track_data           = new TH1D ("h_pre_cut_collection_hits_track_data",           "h_pre_cut_collection_hits_track_data",           20, 0, 250);

TH1D * h_pre_cut_collection_hits_shower_nue_cc         = new TH1D ("h_pre_cut_collection_hits_shower_nue_cc",         "h_pre_cut_collection_hits_shower_nue_cc",         20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_nue_cc_out_fv  = new TH1D ("h_pre_cut_collection_hits_shower_nue_cc_out_fv",  "h_pre_cut_collection_hits_shower_nue_cc_out_fv",  20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_nue_cc_mixed   = new TH1D ("h_pre_cut_collection_hits_shower_nue_cc_mixed",   "h_pre_cut_collection_hits_shower_nue_cc_mixed",   20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_numu_cc        = new TH1D ("h_pre_cut_collection_hits_shower_numu_cc",        "h_pre_cut_collection_hits_shower_numu_cc",        20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_numu_cc_mixed  = new TH1D ("h_pre_cut_collection_hits_shower_numu_cc_mixed",  "h_pre_cut_collection_hits_shower_numu_cc_mixed",  20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_nc             = new TH1D ("h_pre_cut_collection_hits_shower_nc",             "h_pre_cut_collection_hits_shower_nc",             20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_nc_pi0         = new TH1D ("h_pre_cut_collection_hits_shower_nc_pi0",         "h_pre_cut_collection_hits_shower_nc_pi0",         20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_cosmic         = new TH1D ("h_pre_cut_collection_hits_shower_cosmic",         "h_pre_cut_collection_hits_shower_cosmic",         20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_other_mixed    = new TH1D ("h_pre_cut_collection_hits_shower_other_mixed",    "h_pre_cut_collection_hits_shower_other_mixed",    20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_unmatched      = new TH1D ("h_pre_cut_collection_hits_shower_unmatched",      "h_pre_cut_collection_hits_shower_unmatched",      20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_intime         = new TH1D ("h_pre_cut_collection_hits_shower_intime",         "h_pre_cut_collection_hits_shower_intime",         20, 0, 250);
TH1D * h_pre_cut_collection_hits_shower_data           = new TH1D ("h_pre_cut_collection_hits_shower_data",           "h_pre_cut_collection_hits_shower_data",           20, 0, 250);

TH1D * h_pre_cut_collection_hits_leading_shower_nue_cc         = new TH1D ("h_pre_cut_collection_hits_leading_shower_nue_cc",         "h_pre_cut_collection_hits_leading_shower_nue_cc",         30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv  = new TH1D ("h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv",  "h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv",  30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_nue_cc_mixed   = new TH1D ("h_pre_cut_collection_hits_leading_shower_nue_cc_mixed",   "h_pre_cut_collection_hits_leading_shower_nue_cc_mixed",   30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_numu_cc        = new TH1D ("h_pre_cut_collection_hits_leading_shower_numu_cc",        "h_pre_cut_collection_hits_leading_shower_numu_cc",        30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_numu_cc_mixed  = new TH1D ("h_pre_cut_collection_hits_leading_shower_numu_cc_mixed",  "h_pre_cut_collection_hits_leading_shower_numu_cc_mixed",  30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_nc             = new TH1D ("h_pre_cut_collection_hits_leading_shower_nc",             "h_pre_cut_collection_hits_leading_shower_nc",             30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_nc_pi0         = new TH1D ("h_pre_cut_collection_hits_leading_shower_nc_pi0",         "h_pre_cut_collection_hits_leading_shower_nc_pi0",         30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_cosmic         = new TH1D ("h_pre_cut_collection_hits_leading_shower_cosmic",         "h_pre_cut_collection_hits_leading_shower_cosmic",         30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_other_mixed    = new TH1D ("h_pre_cut_collection_hits_leading_shower_other_mixed",    "h_pre_cut_collection_hits_leading_shower_other_mixed",    30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_unmatched      = new TH1D ("h_pre_cut_collection_hits_leading_shower_unmatched",      "h_pre_cut_collection_hits_leading_shower_unmatched",      30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_intime         = new TH1D ("h_pre_cut_collection_hits_leading_shower_intime",         "h_pre_cut_collection_hits_leading_shower_intime",         30, 0, 300);
TH1D * h_pre_cut_collection_hits_leading_shower_data           = new TH1D ("h_pre_cut_collection_hits_leading_shower_data",           "h_pre_cut_collection_hits_leading_shower_data",           30, 0, 300);

TH1D * h_pre_cut_total_hits_leading_shower_nue_cc         = new TH1D ("h_pre_cut_total_hits_leading_shower_nue_cc",         "h_pre_cut_total_hits_leading_shower_nue_cc",         30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_nue_cc_out_fv  = new TH1D ("h_pre_cut_total_hits_leading_shower_nue_cc_out_fv",  "h_pre_cut_total_hits_leading_shower_nue_cc_out_fv",  30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_nue_cc_mixed   = new TH1D ("h_pre_cut_total_hits_leading_shower_nue_cc_mixed",   "h_pre_cut_total_hits_leading_shower_nue_cc_mixed",   30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_numu_cc        = new TH1D ("h_pre_cut_total_hits_leading_shower_numu_cc",        "h_pre_cut_total_hits_leading_shower_numu_cc",        30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_numu_cc_mixed  = new TH1D ("h_pre_cut_total_hits_leading_shower_numu_cc_mixed",  "h_pre_cut_total_hits_leading_shower_numu_cc_mixed",  30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_nc             = new TH1D ("h_pre_cut_total_hits_leading_shower_nc",             "h_pre_cut_total_hits_leading_shower_nc",             30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_nc_pi0         = new TH1D ("h_pre_cut_total_hits_leading_shower_nc_pi0",         "h_pre_cut_total_hits_leading_shower_nc_pi0",         30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_cosmic         = new TH1D ("h_pre_cut_total_hits_leading_shower_cosmic",         "h_pre_cut_total_hits_leading_shower_cosmic",         30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_other_mixed    = new TH1D ("h_pre_cut_total_hits_leading_shower_other_mixed",    "h_pre_cut_total_hits_leading_shower_other_mixed",    30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_unmatched      = new TH1D ("h_pre_cut_total_hits_leading_shower_unmatched",      "h_pre_cut_total_hits_leading_shower_unmatched",      30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_intime         = new TH1D ("h_pre_cut_total_hits_leading_shower_intime",         "h_pre_cut_total_hits_leading_shower_intime",         30, 0, 600);
TH1D * h_pre_cut_total_hits_leading_shower_data           = new TH1D ("h_pre_cut_total_hits_leading_shower_data",           "h_pre_cut_total_hits_leading_shower_data",           30, 0, 600);

TH1D * h_failure_reason_nue_cc         = new TH1D ("h_failure_reason_nue_cc",         "h_failure_reason_nue_cc",         24, 0, 12);
TH1D * h_failure_reason_nue_cc_out_fv  = new TH1D ("h_failure_reason_nue_cc_out_fv",  "h_failure_reason_nue_cc_out_fv",  24, 0, 12);
TH1D * h_failure_reason_nue_cc_mixed   = new TH1D ("h_failure_reason_nue_cc_mixed",   "h_failure_reason_nue_cc_mixed",   24, 0, 12);
TH1D * h_failure_reason_numu_cc        = new TH1D ("h_failure_reason_numu_cc",        "h_failure_reason_numu_cc",        24, 0, 12);
TH1D * h_failure_reason_numu_cc_mixed  = new TH1D ("h_failure_reason_numu_cc_mixed",  "h_failure_reason_numu_cc_mixed",  24, 0, 12);
TH1D * h_failure_reason_nc             = new TH1D ("h_failure_reason_nc",             "h_failure_reason_nc",             24, 0, 12);
TH1D * h_failure_reason_nc_pi0         = new TH1D ("h_failure_reason_nc_pi0",         "h_failure_reason_nc_pi0",         24, 0, 12);
TH1D * h_failure_reason_cosmic         = new TH1D ("h_failure_reason_cosmic",         "h_failure_reason_cosmic",         24, 0, 12);
TH1D * h_failure_reason_other_mixed    = new TH1D ("h_failure_reason_other_mixed",    "h_failure_reason_other_mixed",    24, 0, 12);
TH1D * h_failure_reason_unmatched      = new TH1D ("h_failure_reason_unmatched",      "h_failure_reason_unmatched",      24, 0, 12);

TH1D * h_post_cuts_num_showers_purity_qe     = new TH1D("h_post_cuts_num_showers_purity_qe",             "h_post_cuts_num_showers_purity_qe",           4, 1, 5);
TH1D * h_post_cuts_num_showers_purity_res    = new TH1D("h_post_cuts_num_showers_purity_res",            "h_post_cuts_num_showers_purity_res",          4, 1, 5);
TH1D * h_post_cuts_num_showers_purity_dis    = new TH1D("h_post_cuts_num_showers_purity_dis",            "h_post_cuts_num_showers_purity_dis",          4, 1, 5);
TH1D * h_post_cuts_num_showers_purity_coh    = new TH1D("h_post_cuts_num_showers_purity_coh",            "h_post_cuts_num_showers_purity_coh",          4, 1, 5);
TH1D * h_post_cuts_num_showers_purity_mec    = new TH1D("h_post_cuts_num_showers_purity_mec",            "h_post_cuts_num_showers_purity_mec",          4, 1, 5);

TH1D * h_post_open_angle_cuts_num_showers_purity_qe   = new TH1D("h_post_open_angle_cuts_num_showers_purity_qe",  "h_post_open_angle_cuts_num_showers_purity", 4, 1, 5);
TH1D * h_post_open_angle_cuts_num_showers_purity_res  = new TH1D("h_post_open_angle_cuts_num_showers_purity_res", "h_post_open_angle_cuts_num_showers_purity", 4, 1, 5);
TH1D * h_post_open_angle_cuts_num_showers_purity_dis  = new TH1D("h_post_open_angle_cuts_num_showers_purity_dis", "h_post_open_angle_cuts_num_showers_purity", 4, 1, 5);
TH1D * h_post_open_angle_cuts_num_showers_purity_coh  = new TH1D("h_post_open_angle_cuts_num_showers_purity_coh", "h_post_open_angle_cuts_num_showers_purity", 4, 1, 5);
TH1D * h_post_open_angle_cuts_num_showers_purity_mec  = new TH1D("h_post_open_angle_cuts_num_showers_purity_mec", "h_post_open_angle_cuts_num_showers_purity", 4, 1, 5);

TH2D * h_ele_eng_total_hits       = new TH2D ("h_ele_eng_total_hits",      "h_ele_eng_total_hits",      20, 0, 600, 20, 0, 3);
TH2D * h_ele_eng_colleciton_hits  = new TH2D ("h_ele_eng_colleciton_hits", "h_ele_eng_colleciton_hits", 20, 0, 300, 20, 0, 3);
TH2D * h_nu_eng_total_hits        = new TH2D ("h_nu_eng_total_hits",       "h_nu_eng_total_hits",       20, 0, 600, 20, 0, 3);
TH2D * h_nu_eng_collection_hits   = new TH2D ("h_nu_eng_collection_hits",  "h_nu_eng_collection_hits",  20, 0, 300, 20, 0, 3);

TH1D * h_ele_cos_theta_nue_cc         = new TH1D ("h_ele_cos_theta_nue_cc",         "h_ele_cos_theta_nue_cc",         20, -1, 1);
TH1D * h_ele_cos_theta_nue_cc_out_fv  = new TH1D ("h_ele_cos_theta_nue_cc_out_fv",  "h_ele_cos_theta_nue_cc_out_fv",  20, -1, 1);
TH1D * h_ele_cos_theta_nue_cc_mixed   = new TH1D ("h_ele_cos_theta_nue_cc_mixed",   "h_ele_cos_theta_nue_cc_mixed",   20, -1, 1);
TH1D * h_ele_cos_theta_numu_cc        = new TH1D ("h_ele_cos_theta_numu_cc",        "h_ele_cos_theta_numu_cc",        20, -1, 1);
TH1D * h_ele_cos_theta_numu_cc_mixed  = new TH1D ("h_ele_cos_theta_numu_cc_mixed",  "h_ele_cos_theta_numu_cc_mixed",  20, -1, 1);
TH1D * h_ele_cos_theta_nc             = new TH1D ("h_ele_cos_theta_nc",             "h_ele_cos_theta_nc",             20, -1, 1);
TH1D * h_ele_cos_theta_nc_pi0         = new TH1D ("h_ele_cos_theta_nc_pi0",         "h_ele_cos_theta_nc_pi0",         20, -1, 1);
TH1D * h_ele_cos_theta_cosmic         = new TH1D ("h_ele_cos_theta_cosmic",         "h_ele_cos_theta_cosmic",         20, -1, 1);
TH1D * h_ele_cos_theta_other_mixed    = new TH1D ("h_ele_cos_theta_other_mixed",    "h_ele_cos_theta_other_mixed",    20, -1, 1);
TH1D * h_ele_cos_theta_unmatched      = new TH1D ("h_ele_cos_theta_unmatched",      "h_ele_cos_theta_unmatched",      20, -1, 1);
TH1D * h_ele_cos_theta_intime         = new TH1D ("h_ele_cos_theta_intime",         "h_ele_cos_theta_intime",         20, -1, 1);
TH1D * h_ele_cos_theta_data           = new TH1D ("h_ele_cos_theta_data",           "h_ele_cos_theta_data",           20, -1, 1);

TH1D * h_ele_cos_theta_last_nue_cc         = new TH1D ("h_ele_cos_theta_last_nue_cc",         "h_ele_cos_theta_last_nue_cc",         20, -1, 1);
TH1D * h_ele_cos_theta_last_nue_cc_out_fv  = new TH1D ("h_ele_cos_theta_last_nue_cc_out_fv",  "h_ele_cos_theta_last_nue_cc_out_fv",  20, -1, 1);
TH1D * h_ele_cos_theta_last_nue_cc_mixed   = new TH1D ("h_ele_cos_theta_last_nue_cc_mixed",   "h_ele_cos_theta_last_nue_cc_mixed",   20, -1, 1);
TH1D * h_ele_cos_theta_last_numu_cc        = new TH1D ("h_ele_cos_theta_last_numu_cc",        "h_ele_cos_theta_last_numu_cc",        20, -1, 1);
TH1D * h_ele_cos_theta_last_numu_cc_mixed  = new TH1D ("h_ele_cos_theta_last_numu_cc_mixed",  "h_ele_cos_theta_last_numu_cc_mixed",  20, -1, 1);
TH1D * h_ele_cos_theta_last_nc             = new TH1D ("h_ele_cos_theta_last_nc",             "h_ele_cos_theta_last_nc",             20, -1, 1);
TH1D * h_ele_cos_theta_last_nc_pi0         = new TH1D ("h_ele_cos_theta_last_nc_pi0",         "h_ele_cos_theta_last_nc_pi0",         20, -1, 1);
TH1D * h_ele_cos_theta_last_cosmic         = new TH1D ("h_ele_cos_theta_last_cosmic",         "h_ele_cos_theta_last_cosmic",         20, -1, 1);
TH1D * h_ele_cos_theta_last_other_mixed    = new TH1D ("h_ele_cos_theta_last_other_mixed",    "h_ele_cos_theta_last_other_mixed",    20, -1, 1);
TH1D * h_ele_cos_theta_last_unmatched      = new TH1D ("h_ele_cos_theta_last_unmatched",      "h_ele_cos_theta_last_unmatched",      20, -1, 1);
TH1D * h_ele_cos_theta_last_intime         = new TH1D ("h_ele_cos_theta_last_intime",         "h_ele_cos_theta_last_intime",         20, -1, 1);
TH1D * h_ele_cos_theta_last_data           = new TH1D ("h_ele_cos_theta_last_data",           "h_ele_cos_theta_last_data",           20, -1, 1);

TH1D * h_ele_cos_theta_last_trans_nue_cc         = new TH1D ("h_ele_cos_theta_last_trans_nue_cc",         "h_ele_cos_theta_last_trans_nue_cc",         20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_nue_cc_out_fv  = new TH1D ("h_ele_cos_theta_last_trans_nue_cc_out_fv",  "h_ele_cos_theta_last_trans_nue_cc_out_fv",  20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_nue_cc_mixed   = new TH1D ("h_ele_cos_theta_last_trans_nue_cc_mixed",   "h_ele_cos_theta_last_trans_nue_cc_mixed",   20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_numu_cc        = new TH1D ("h_ele_cos_theta_last_trans_numu_cc",        "h_ele_cos_theta_last_trans_numu_cc",        20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_numu_cc_mixed  = new TH1D ("h_ele_cos_theta_last_trans_numu_cc_mixed",  "h_ele_cos_theta_last_trans_numu_cc_mixed",  20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_nc             = new TH1D ("h_ele_cos_theta_last_trans_nc",             "h_ele_cos_theta_last_trans_nc",             20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_nc_pi0         = new TH1D ("h_ele_cos_theta_last_trans_nc_pi0",         "h_ele_cos_theta_last_trans_nc_pi0",         20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_cosmic         = new TH1D ("h_ele_cos_theta_last_trans_cosmic",         "h_ele_cos_theta_last_trans_cosmic",         20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_other_mixed    = new TH1D ("h_ele_cos_theta_last_trans_other_mixed",    "h_ele_cos_theta_last_trans_other_mixed",    20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_unmatched      = new TH1D ("h_ele_cos_theta_last_trans_unmatched",      "h_ele_cos_theta_last_trans_unmatched",      20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_intime         = new TH1D ("h_ele_cos_theta_last_trans_intime",         "h_ele_cos_theta_last_trans_intime",         20, -1, 1);
TH1D * h_ele_cos_theta_last_trans_data           = new TH1D ("h_ele_cos_theta_last_trans_data",           "h_ele_cos_theta_last_trans_data",           20, -1, 1);

TH1D * h_ele_pfp_momentum_nue_cc         = new TH1D ("h_ele_pfp_momentum_nue_cc",         "h_ele_pfp_momentum_nue_cc",         10, 0, 2);
TH1D * h_ele_pfp_momentum_nue_cc_out_fv  = new TH1D ("h_ele_pfp_momentum_nue_cc_out_fv",  "h_ele_pfp_momentum_nue_cc_out_fv",  10, 0, 2);
TH1D * h_ele_pfp_momentum_nue_cc_mixed   = new TH1D ("h_ele_pfp_momentum_nue_cc_mixed",   "h_ele_pfp_momentum_nue_cc_mixed",   10, 0, 2);
TH1D * h_ele_pfp_momentum_numu_cc        = new TH1D ("h_ele_pfp_momentum_numu_cc",        "h_ele_pfp_momentum_numu_cc",        10, 0, 2);
TH1D * h_ele_pfp_momentum_numu_cc_mixed  = new TH1D ("h_ele_pfp_momentum_numu_cc_mixed",  "h_ele_pfp_momentum_numu_cc_mixed",  10, 0, 2);
TH1D * h_ele_pfp_momentum_nc             = new TH1D ("h_ele_pfp_momentum_nc",             "h_ele_pfp_momentum_nc",             10, 0, 2);
TH1D * h_ele_pfp_momentum_nc_pi0         = new TH1D ("h_ele_pfp_momentum_nc_pi0",         "h_ele_pfp_momentum_nc_pi0",         10, 0, 2);
TH1D * h_ele_pfp_momentum_cosmic         = new TH1D ("h_ele_pfp_momentum_cosmic",         "h_ele_pfp_momentum_cosmic",         10, 0, 2);
TH1D * h_ele_pfp_momentum_other_mixed    = new TH1D ("h_ele_pfp_momentum_other_mixed",    "h_ele_pfp_momentum_other_mixed",    10, 0, 2);
TH1D * h_ele_pfp_momentum_unmatched      = new TH1D ("h_ele_pfp_momentum_unmatched",      "h_ele_pfp_momentum_unmatched",      10, 0, 2);
TH1D * h_ele_pfp_momentum_intime         = new TH1D ("h_ele_pfp_momentum_intime",         "h_ele_pfp_momentum_intime",         10, 0, 2);
TH1D * h_ele_pfp_momentum_data           = new TH1D ("h_ele_pfp_momentum_data",           "h_ele_pfp_momentum_data",           10, 0, 2);

TH1D * h_ele_pfp_theta_nue_cc         = new TH1D ("h_ele_pfp_theta_nue_cc",         "h_ele_pfp_theta_nue_cc",         20, 0, 180);
TH1D * h_ele_pfp_theta_nue_cc_out_fv  = new TH1D ("h_ele_pfp_theta_nue_cc_out_fv",  "h_ele_pfp_theta_nue_cc_out_fv",  20, 0, 180);
TH1D * h_ele_pfp_theta_nue_cc_mixed   = new TH1D ("h_ele_pfp_theta_nue_cc_mixed",   "h_ele_pfp_theta_nue_cc_mixed",   20, 0, 180);
TH1D * h_ele_pfp_theta_numu_cc        = new TH1D ("h_ele_pfp_theta_numu_cc",        "h_ele_pfp_theta_numu_cc",        20, 0, 180);
TH1D * h_ele_pfp_theta_numu_cc_mixed  = new TH1D ("h_ele_pfp_theta_numu_cc_mixed",  "h_ele_pfp_theta_numu_cc_mixed",  20, 0, 180);
TH1D * h_ele_pfp_theta_nc             = new TH1D ("h_ele_pfp_theta_nc",             "h_ele_pfp_theta_nc",             20, 0, 180);
TH1D * h_ele_pfp_theta_nc_pi0         = new TH1D ("h_ele_pfp_theta_nc_pi0",         "h_ele_pfp_theta_nc_pi0",         20, 0, 180);
TH1D * h_ele_pfp_theta_cosmic         = new TH1D ("h_ele_pfp_theta_cosmic",         "h_ele_pfp_theta_cosmic",         20, 0, 180);
TH1D * h_ele_pfp_theta_other_mixed    = new TH1D ("h_ele_pfp_theta_other_mixed",    "h_ele_pfp_theta_other_mixed",    20, 0, 180);
TH1D * h_ele_pfp_theta_unmatched      = new TH1D ("h_ele_pfp_theta_unmatched",      "h_ele_pfp_theta_unmatched",      20, 0, 180);
TH1D * h_ele_pfp_theta_intime         = new TH1D ("h_ele_pfp_theta_intime",         "h_ele_pfp_theta_intime",         20, 0, 180);
TH1D * h_ele_pfp_theta_data           = new TH1D ("h_ele_pfp_theta_data",           "h_ele_pfp_theta_data",           20, 0, 180);

TH1D * h_ele_pfp_theta_after_nue_cc         = new TH1D ("h_ele_pfp_theta_after_nue_cc",         "h_ele_pfp_theta_after_nue_cc",         20, 0, 180);
TH1D * h_ele_pfp_theta_after_nue_cc_out_fv  = new TH1D ("h_ele_pfp_theta_after_nue_cc_out_fv",  "h_ele_pfp_theta_after_nue_cc_out_fv",  20, 0, 180);
TH1D * h_ele_pfp_theta_after_nue_cc_mixed   = new TH1D ("h_ele_pfp_theta_after_nue_cc_mixed",   "h_ele_pfp_theta_after_nue_cc_mixed",   20, 0, 180);
TH1D * h_ele_pfp_theta_after_numu_cc        = new TH1D ("h_ele_pfp_theta_after_numu_cc",        "h_ele_pfp_theta_after_numu_cc",        20, 0, 180);
TH1D * h_ele_pfp_theta_after_numu_cc_mixed  = new TH1D ("h_ele_pfp_theta_after_numu_cc_mixed",  "h_ele_pfp_theta_after_numu_cc_mixed",  20, 0, 180);
TH1D * h_ele_pfp_theta_after_nc             = new TH1D ("h_ele_pfp_theta_after_nc",             "h_ele_pfp_theta_after_nc",             20, 0, 180);
TH1D * h_ele_pfp_theta_after_nc_pi0         = new TH1D ("h_ele_pfp_theta_after_nc_pi0",         "h_ele_pfp_theta_after_nc_pi0",         20, 0, 180);
TH1D * h_ele_pfp_theta_after_cosmic         = new TH1D ("h_ele_pfp_theta_after_cosmic",         "h_ele_pfp_theta_after_cosmic",         20, 0, 180);
TH1D * h_ele_pfp_theta_after_other_mixed    = new TH1D ("h_ele_pfp_theta_after_other_mixed",    "h_ele_pfp_theta_after_other_mixed",    20, 0, 180);
TH1D * h_ele_pfp_theta_after_unmatched      = new TH1D ("h_ele_pfp_theta_after_unmatched",      "h_ele_pfp_theta_after_unmatched",      20, 0, 180);
TH1D * h_ele_pfp_theta_after_intime         = new TH1D ("h_ele_pfp_theta_after_intime",         "h_ele_pfp_theta_after_intime",         20, 0, 180);
TH1D * h_ele_pfp_theta_after_data           = new TH1D ("h_ele_pfp_theta_after_data",           "h_ele_pfp_theta_after_data",           20, 0, 180);

TH1D * h_ele_pfp_theta_last_nue_cc         = new TH1D ("h_ele_pfp_theta_last_nue_cc",         "h_ele_pfp_theta_last_nue_cc",         20, 0, 180);
TH1D * h_ele_pfp_theta_last_nue_cc_out_fv  = new TH1D ("h_ele_pfp_theta_last_nue_cc_out_fv",  "h_ele_pfp_theta_last_nue_cc_out_fv",  20, 0, 180);
TH1D * h_ele_pfp_theta_last_nue_cc_mixed   = new TH1D ("h_ele_pfp_theta_last_nue_cc_mixed",   "h_ele_pfp_theta_last_nue_cc_mixed",   20, 0, 180);
TH1D * h_ele_pfp_theta_last_numu_cc        = new TH1D ("h_ele_pfp_theta_last_numu_cc",        "h_ele_pfp_theta_last_numu_cc",        20, 0, 180);
TH1D * h_ele_pfp_theta_last_numu_cc_mixed  = new TH1D ("h_ele_pfp_theta_last_numu_cc_mixed",  "h_ele_pfp_theta_last_numu_cc_mixed",  20, 0, 180);
TH1D * h_ele_pfp_theta_last_nc             = new TH1D ("h_ele_pfp_theta_last_nc",             "h_ele_pfp_theta_last_nc",             20, 0, 180);
TH1D * h_ele_pfp_theta_last_nc_pi0         = new TH1D ("h_ele_pfp_theta_last_nc_pi0",         "h_ele_pfp_theta_last_nc_pi0",         20, 0, 180);
TH1D * h_ele_pfp_theta_last_cosmic         = new TH1D ("h_ele_pfp_theta_last_cosmic",         "h_ele_pfp_theta_last_cosmic",         20, 0, 180);
TH1D * h_ele_pfp_theta_last_other_mixed    = new TH1D ("h_ele_pfp_theta_last_other_mixed",    "h_ele_pfp_theta_last_other_mixed",    20, 0, 180);
TH1D * h_ele_pfp_theta_last_unmatched      = new TH1D ("h_ele_pfp_theta_last_unmatched",      "h_ele_pfp_theta_last_unmatched",      20, 0, 180);
TH1D * h_ele_pfp_theta_last_intime         = new TH1D ("h_ele_pfp_theta_last_intime",         "h_ele_pfp_theta_last_intime",         20, 0, 180);
TH1D * h_ele_pfp_theta_last_data           = new TH1D ("h_ele_pfp_theta_last_data",           "h_ele_pfp_theta_last_data",           20, 0, 180);

TH1D * h_ele_pfp_phi_nue_cc         = new TH1D ("h_ele_pfp_phi_nue_cc",         "h_ele_pfp_phi_nue_cc",         20, -180, 180);
TH1D * h_ele_pfp_phi_nue_cc_out_fv  = new TH1D ("h_ele_pfp_phi_nue_cc_out_fv",  "h_ele_pfp_phi_nue_cc_out_fv",  20, -180, 180);
TH1D * h_ele_pfp_phi_nue_cc_mixed   = new TH1D ("h_ele_pfp_phi_nue_cc_mixed",   "h_ele_pfp_phi_nue_cc_mixed",   20, -180, 180);
TH1D * h_ele_pfp_phi_numu_cc        = new TH1D ("h_ele_pfp_phi_numu_cc",        "h_ele_pfp_phi_numu_cc",        20, -180, 180);
TH1D * h_ele_pfp_phi_numu_cc_mixed  = new TH1D ("h_ele_pfp_phi_numu_cc_mixed",  "h_ele_pfp_phi_numu_cc_mixed",  20, -180, 180);
TH1D * h_ele_pfp_phi_nc             = new TH1D ("h_ele_pfp_phi_nc",             "h_ele_pfp_phi_nc",             20, -180, 180);
TH1D * h_ele_pfp_phi_nc_pi0         = new TH1D ("h_ele_pfp_phi_nc_pi0",         "h_ele_pfp_phi_nc_pi0",         20, -180, 180);
TH1D * h_ele_pfp_phi_cosmic         = new TH1D ("h_ele_pfp_phi_cosmic",         "h_ele_pfp_phi_cosmic",         20, -180, 180);
TH1D * h_ele_pfp_phi_other_mixed    = new TH1D ("h_ele_pfp_phi_other_mixed",    "h_ele_pfp_phi_other_mixed",    20, -180, 180);
TH1D * h_ele_pfp_phi_unmatched      = new TH1D ("h_ele_pfp_phi_unmatched",      "h_ele_pfp_phi_unmatched",      20, -180, 180);
TH1D * h_ele_pfp_phi_intime         = new TH1D ("h_ele_pfp_phi_intime",         "h_ele_pfp_phi_intime",         20, -180, 180);
TH1D * h_ele_pfp_phi_data           = new TH1D ("h_ele_pfp_phi_data",           "h_ele_pfp_phi_data",           20, -180, 180);

TH1D * h_ele_pfp_phi_after_nue_cc         = new TH1D ("h_ele_pfp_phi_after_nue_cc",         "h_ele_pfp_phi_after_nue_cc",         20, -180, 180);
TH1D * h_ele_pfp_phi_after_nue_cc_out_fv  = new TH1D ("h_ele_pfp_phi_after_nue_cc_out_fv",  "h_ele_pfp_phi_after_nue_cc_out_fv",  20, -180, 180);
TH1D * h_ele_pfp_phi_after_nue_cc_mixed   = new TH1D ("h_ele_pfp_phi_after_nue_cc_mixed",   "h_ele_pfp_phi_after_nue_cc_mixed",   20, -180, 180);
TH1D * h_ele_pfp_phi_after_numu_cc        = new TH1D ("h_ele_pfp_phi_after_numu_cc",        "h_ele_pfp_phi_after_numu_cc",        20, -180, 180);
TH1D * h_ele_pfp_phi_after_numu_cc_mixed  = new TH1D ("h_ele_pfp_phi_after_numu_cc_mixed",  "h_ele_pfp_phi_after_numu_cc_mixed",  20, -180, 180);
TH1D * h_ele_pfp_phi_after_nc             = new TH1D ("h_ele_pfp_phi_after_nc",             "h_ele_pfp_phi_after_nc",             20, -180, 180);
TH1D * h_ele_pfp_phi_after_nc_pi0         = new TH1D ("h_ele_pfp_phi_after_nc_pi0",         "h_ele_pfp_phi_after_nc_pi0",         20, -180, 180);
TH1D * h_ele_pfp_phi_after_cosmic         = new TH1D ("h_ele_pfp_phi_after_cosmic",         "h_ele_pfp_phi_after_cosmic",         20, -180, 180);
TH1D * h_ele_pfp_phi_after_other_mixed    = new TH1D ("h_ele_pfp_phi_after_other_mixed",    "h_ele_pfp_phi_after_other_mixed",    20, -180, 180);
TH1D * h_ele_pfp_phi_after_unmatched      = new TH1D ("h_ele_pfp_phi_after_unmatched",      "h_ele_pfp_phi_after_unmatched",      20, -180, 180);
TH1D * h_ele_pfp_phi_after_intime         = new TH1D ("h_ele_pfp_phi_after_intime",         "h_ele_pfp_phi_after_intime",         20, -180, 180);
TH1D * h_ele_pfp_phi_after_data           = new TH1D ("h_ele_pfp_phi_after_data",           "h_ele_pfp_phi_after_data",           20, -180, 180);

TH1D * h_ele_pfp_phi_last_nue_cc         = new TH1D ("h_ele_pfp_phi_last_nue_cc",         "h_ele_pfp_phi_last_nue_cc",         20, -180, 180);
TH1D * h_ele_pfp_phi_last_nue_cc_out_fv  = new TH1D ("h_ele_pfp_phi_last_nue_cc_out_fv",  "h_ele_pfp_phi_last_nue_cc_out_fv",  20, -180, 180);
TH1D * h_ele_pfp_phi_last_nue_cc_mixed   = new TH1D ("h_ele_pfp_phi_last_nue_cc_mixed",   "h_ele_pfp_phi_last_nue_cc_mixed",   20, -180, 180);
TH1D * h_ele_pfp_phi_last_numu_cc        = new TH1D ("h_ele_pfp_phi_last_numu_cc",        "h_ele_pfp_phi_last_numu_cc",        20, -180, 180);
TH1D * h_ele_pfp_phi_last_numu_cc_mixed  = new TH1D ("h_ele_pfp_phi_last_numu_cc_mixed",  "h_ele_pfp_phi_last_numu_cc_mixed",  20, -180, 180);
TH1D * h_ele_pfp_phi_last_nc             = new TH1D ("h_ele_pfp_phi_last_nc",             "h_ele_pfp_phi_last_nc",             20, -180, 180);
TH1D * h_ele_pfp_phi_last_nc_pi0         = new TH1D ("h_ele_pfp_phi_last_nc_pi0",         "h_ele_pfp_phi_last_nc_pi0",         20, -180, 180);
TH1D * h_ele_pfp_phi_last_cosmic         = new TH1D ("h_ele_pfp_phi_last_cosmic",         "h_ele_pfp_phi_last_cosmic",         20, -180, 180);
TH1D * h_ele_pfp_phi_last_other_mixed    = new TH1D ("h_ele_pfp_phi_last_other_mixed",    "h_ele_pfp_phi_last_other_mixed",    20, -180, 180);
TH1D * h_ele_pfp_phi_last_unmatched      = new TH1D ("h_ele_pfp_phi_last_unmatched",      "h_ele_pfp_phi_last_unmatched",      20, -180, 180);
TH1D * h_ele_pfp_phi_last_intime         = new TH1D ("h_ele_pfp_phi_last_intime",         "h_ele_pfp_phi_last_intime",         20, -180, 180);
TH1D * h_ele_pfp_phi_last_data           = new TH1D ("h_ele_pfp_phi_last_data",           "h_ele_pfp_phi_last_data",           20, -180, 180);

TH1D * h_leading_shwr_length_1shwr_nue_cc         = new TH1D ("h_leading_shwr_length_1shwr_nue_cc",         "h_leading_shwr_length_1shwr_nue_cc",         20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_nue_cc_out_fv  = new TH1D ("h_leading_shwr_length_1shwr_nue_cc_out_fv",  "h_leading_shwr_length_1shwr_nue_cc_out_fv",  20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_nue_cc_mixed   = new TH1D ("h_leading_shwr_length_1shwr_nue_cc_mixed",   "h_leading_shwr_length_1shwr_nue_cc_mixed",   20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_numu_cc        = new TH1D ("h_leading_shwr_length_1shwr_numu_cc",        "h_leading_shwr_length_1shwr_numu_cc",        20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_numu_cc_mixed  = new TH1D ("h_leading_shwr_length_1shwr_numu_cc_mixed",  "h_leading_shwr_length_1shwr_numu_cc_mixed",  20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_nc             = new TH1D ("h_leading_shwr_length_1shwr_nc",             "h_leading_shwr_length_1shwr_nc",             20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_nc_pi0         = new TH1D ("h_leading_shwr_length_1shwr_nc_pi0",         "h_leading_shwr_length_1shwr_nc_pi0",         20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_cosmic         = new TH1D ("h_leading_shwr_length_1shwr_cosmic",         "h_leading_shwr_length_1shwr_cosmic",         20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_other_mixed    = new TH1D ("h_leading_shwr_length_1shwr_other_mixed",    "h_leading_shwr_length_1shwr_other_mixed",    20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_unmatched      = new TH1D ("h_leading_shwr_length_1shwr_unmatched",      "h_leading_shwr_length_1shwr_unmatched",      20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_intime         = new TH1D ("h_leading_shwr_length_1shwr_intime",         "h_leading_shwr_length_1shwr_intime",         20, 0, 300);
TH1D * h_leading_shwr_length_1shwr_data           = new TH1D ("h_leading_shwr_length_1shwr_data",           "h_leading_shwr_length_1shwr_data",           20, 0, 300);

TH1D * h_leading_shwr_length_2shwr_nue_cc         = new TH1D ("h_leading_shwr_length_2shwr_nue_cc",         "h_leading_shwr_length_2shwr_nue_cc",         20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_nue_cc_out_fv  = new TH1D ("h_leading_shwr_length_2shwr_nue_cc_out_fv",  "h_leading_shwr_length_2shwr_nue_cc_out_fv",  20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_nue_cc_mixed   = new TH1D ("h_leading_shwr_length_2shwr_nue_cc_mixed",   "h_leading_shwr_length_2shwr_nue_cc_mixed",   20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_numu_cc        = new TH1D ("h_leading_shwr_length_2shwr_numu_cc",        "h_leading_shwr_length_2shwr_numu_cc",        20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_numu_cc_mixed  = new TH1D ("h_leading_shwr_length_2shwr_numu_cc_mixed",  "h_leading_shwr_length_2shwr_numu_cc_mixed",  20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_nc             = new TH1D ("h_leading_shwr_length_2shwr_nc",             "h_leading_shwr_length_2shwr_nc",             20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_nc_pi0         = new TH1D ("h_leading_shwr_length_2shwr_nc_pi0",         "h_leading_shwr_length_2shwr_nc_pi0",         20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_cosmic         = new TH1D ("h_leading_shwr_length_2shwr_cosmic",         "h_leading_shwr_length_2shwr_cosmic",         20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_other_mixed    = new TH1D ("h_leading_shwr_length_2shwr_other_mixed",    "h_leading_shwr_length_2shwr_other_mixed",    20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_unmatched      = new TH1D ("h_leading_shwr_length_2shwr_unmatched",      "h_leading_shwr_length_2shwr_unmatched",      20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_intime         = new TH1D ("h_leading_shwr_length_2shwr_intime",         "h_leading_shwr_length_2shwr_intime",         20, 0, 240);
TH1D * h_leading_shwr_length_2shwr_data           = new TH1D ("h_leading_shwr_length_2shwr_data",           "h_leading_shwr_length_2shwr_data",           20, 0, 240);

TH1D * h_leading_shwr_hits_1shwr_nue_cc         = new TH1D ("h_leading_shwr_hits_1shwr_nue_cc",         "h_leading_shwr_hits_1shwr_nue_cc",         20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_nue_cc_out_fv  = new TH1D ("h_leading_shwr_hits_1shwr_nue_cc_out_fv",  "h_leading_shwr_hits_1shwr_nue_cc_out_fv",  20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_nue_cc_mixed   = new TH1D ("h_leading_shwr_hits_1shwr_nue_cc_mixed",   "h_leading_shwr_hits_1shwr_nue_cc_mixed",   20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_numu_cc        = new TH1D ("h_leading_shwr_hits_1shwr_numu_cc",        "h_leading_shwr_hits_1shwr_numu_cc",        20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_numu_cc_mixed  = new TH1D ("h_leading_shwr_hits_1shwr_numu_cc_mixed",  "h_leading_shwr_hits_1shwr_numu_cc_mixed",  20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_nc             = new TH1D ("h_leading_shwr_hits_1shwr_nc",             "h_leading_shwr_hits_1shwr_nc",             20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_nc_pi0         = new TH1D ("h_leading_shwr_hits_1shwr_nc_pi0",         "h_leading_shwr_hits_1shwr_nc_pi0",         20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_cosmic         = new TH1D ("h_leading_shwr_hits_1shwr_cosmic",         "h_leading_shwr_hits_1shwr_cosmic",         20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_other_mixed    = new TH1D ("h_leading_shwr_hits_1shwr_other_mixed",    "h_leading_shwr_hits_1shwr_other_mixed",    20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_unmatched      = new TH1D ("h_leading_shwr_hits_1shwr_unmatched",      "h_leading_shwr_hits_1shwr_unmatched",      20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_intime         = new TH1D ("h_leading_shwr_hits_1shwr_intime",         "h_leading_shwr_hits_1shwr_intime",         20, 0, 1500);
TH1D * h_leading_shwr_hits_1shwr_data           = new TH1D ("h_leading_shwr_hits_1shwr_data",           "h_leading_shwr_hits_1shwr_data",           20, 0, 1500);

TH1D * h_leading_shwr_hits_2shwr_nue_cc         = new TH1D ("h_leading_shwr_hits_2shwr_nue_cc",         "h_leading_shwr_hits_2shwr_nue_cc",         20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_nue_cc_out_fv  = new TH1D ("h_leading_shwr_hits_2shwr_nue_cc_out_fv",  "h_leading_shwr_hits_2shwr_nue_cc_out_fv",  20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_nue_cc_mixed   = new TH1D ("h_leading_shwr_hits_2shwr_nue_cc_mixed",   "h_leading_shwr_hits_2shwr_nue_cc_mixed",   20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_numu_cc        = new TH1D ("h_leading_shwr_hits_2shwr_numu_cc",        "h_leading_shwr_hits_2shwr_numu_cc",        20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_numu_cc_mixed  = new TH1D ("h_leading_shwr_hits_2shwr_numu_cc_mixed",  "h_leading_shwr_hits_2shwr_numu_cc_mixed",  20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_nc             = new TH1D ("h_leading_shwr_hits_2shwr_nc",             "h_leading_shwr_hits_2shwr_nc",             20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_nc_pi0         = new TH1D ("h_leading_shwr_hits_2shwr_nc_pi0",         "h_leading_shwr_hits_2shwr_nc_pi0",         20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_cosmic         = new TH1D ("h_leading_shwr_hits_2shwr_cosmic",         "h_leading_shwr_hits_2shwr_cosmic",         20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_other_mixed    = new TH1D ("h_leading_shwr_hits_2shwr_other_mixed",    "h_leading_shwr_hits_2shwr_other_mixed",    20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_unmatched      = new TH1D ("h_leading_shwr_hits_2shwr_unmatched",      "h_leading_shwr_hits_2shwr_unmatched",      20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_intime         = new TH1D ("h_leading_shwr_hits_2shwr_intime",         "h_leading_shwr_hits_2shwr_intime",         20, 0, 1500);
TH1D * h_leading_shwr_hits_2shwr_data           = new TH1D ("h_leading_shwr_hits_2shwr_data",           "h_leading_shwr_hits_2shwr_data",           20, 0, 1500);

TH1D * h_ele_pfp_x_nue_cc         = new TH1D ("h_ele_pfp_x_nue_cc",         "h_ele_pfp_x_nue_cc",         20, -10, 270);
TH1D * h_ele_pfp_x_nue_cc_out_fv  = new TH1D ("h_ele_pfp_x_nue_cc_out_fv",  "h_ele_pfp_x_nue_cc_out_fv",  20, -10, 270);
TH1D * h_ele_pfp_x_nue_cc_mixed   = new TH1D ("h_ele_pfp_x_nue_cc_mixed",   "h_ele_pfp_x_nue_cc_mixed",   20, -10, 270);
TH1D * h_ele_pfp_x_numu_cc        = new TH1D ("h_ele_pfp_x_numu_cc",        "h_ele_pfp_x_numu_cc",        20, -10, 270);
TH1D * h_ele_pfp_x_numu_cc_mixed  = new TH1D ("h_ele_pfp_x_numu_cc_mixed",  "h_ele_pfp_x_numu_cc_mixed",  20, -10, 270);
TH1D * h_ele_pfp_x_nc             = new TH1D ("h_ele_pfp_x_nc",             "h_ele_pfp_x_nc",             20, -10, 270);
TH1D * h_ele_pfp_x_nc_pi0         = new TH1D ("h_ele_pfp_x_nc_pi0",         "h_ele_pfp_x_nc_pi0",         20, -10, 270);
TH1D * h_ele_pfp_x_cosmic         = new TH1D ("h_ele_pfp_x_cosmic",         "h_ele_pfp_x_cosmic",         20, -10, 270);
TH1D * h_ele_pfp_x_other_mixed    = new TH1D ("h_ele_pfp_x_other_mixed",    "h_ele_pfp_x_other_mixed",    20, -10, 270);
TH1D * h_ele_pfp_x_unmatched      = new TH1D ("h_ele_pfp_x_unmatched",      "h_ele_pfp_x_unmatched",      20, -10, 270);
TH1D * h_ele_pfp_x_intime         = new TH1D ("h_ele_pfp_x_intime",         "h_ele_pfp_x_intime",         20, -10, 270);
TH1D * h_ele_pfp_x_data           = new TH1D ("h_ele_pfp_x_data",           "h_ele_pfp_x_data",           20, -10, 270);

TH1D * h_ele_pfp_y_nue_cc         = new TH1D ("h_ele_pfp_y_nue_cc",         "h_ele_pfp_y_nue_cc",         20, -120, 120);
TH1D * h_ele_pfp_y_nue_cc_out_fv  = new TH1D ("h_ele_pfp_y_nue_cc_out_fv",  "h_ele_pfp_y_nue_cc_out_fv",  20, -120, 120);
TH1D * h_ele_pfp_y_nue_cc_mixed   = new TH1D ("h_ele_pfp_y_nue_cc_mixed",   "h_ele_pfp_y_nue_cc_mixed",   20, -120, 120);
TH1D * h_ele_pfp_y_numu_cc        = new TH1D ("h_ele_pfp_y_numu_cc",        "h_ele_pfp_y_numu_cc",        20, -120, 120);
TH1D * h_ele_pfp_y_numu_cc_mixed  = new TH1D ("h_ele_pfp_y_numu_cc_mixed",  "h_ele_pfp_y_numu_cc_mixed",  20, -120, 120);
TH1D * h_ele_pfp_y_nc             = new TH1D ("h_ele_pfp_y_nc",             "h_ele_pfp_y_nc",             20, -120, 120);
TH1D * h_ele_pfp_y_nc_pi0         = new TH1D ("h_ele_pfp_y_nc_pi0",         "h_ele_pfp_y_nc_pi0",         20, -120, 120);
TH1D * h_ele_pfp_y_cosmic         = new TH1D ("h_ele_pfp_y_cosmic",         "h_ele_pfp_y_cosmic",         20, -120, 120);
TH1D * h_ele_pfp_y_other_mixed    = new TH1D ("h_ele_pfp_y_other_mixed",    "h_ele_pfp_y_other_mixed",    20, -120, 120);
TH1D * h_ele_pfp_y_unmatched      = new TH1D ("h_ele_pfp_y_unmatched",      "h_ele_pfp_y_unmatched",      20, -120, 120);
TH1D * h_ele_pfp_y_intime         = new TH1D ("h_ele_pfp_y_intime",         "h_ele_pfp_y_intime",         20, -120, 120);
TH1D * h_ele_pfp_y_data           = new TH1D ("h_ele_pfp_y_data",           "h_ele_pfp_y_data",           20, -120, 120);

TH1D * h_ele_pfp_z_nue_cc         = new TH1D ("h_ele_pfp_z_nue_cc",         "h_ele_pfp_z_nue_cc",         40, -40, 1050);
TH1D * h_ele_pfp_z_nue_cc_out_fv  = new TH1D ("h_ele_pfp_z_nue_cc_out_fv",  "h_ele_pfp_z_nue_cc_out_fv",  40, -40, 1050);
TH1D * h_ele_pfp_z_nue_cc_mixed   = new TH1D ("h_ele_pfp_z_nue_cc_mixed",   "h_ele_pfp_z_nue_cc_mixed",   40, -40, 1050);
TH1D * h_ele_pfp_z_numu_cc        = new TH1D ("h_ele_pfp_z_numu_cc",        "h_ele_pfp_z_numu_cc",        40, -40, 1050);
TH1D * h_ele_pfp_z_numu_cc_mixed  = new TH1D ("h_ele_pfp_z_numu_cc_mixed",  "h_ele_pfp_z_numu_cc_mixed",  40, -40, 1050);
TH1D * h_ele_pfp_z_nc             = new TH1D ("h_ele_pfp_z_nc",             "h_ele_pfp_z_nc",             40, -40, 1050);
TH1D * h_ele_pfp_z_nc_pi0         = new TH1D ("h_ele_pfp_z_nc_pi0",         "h_ele_pfp_z_nc_pi0",         40, -40, 1050);
TH1D * h_ele_pfp_z_cosmic         = new TH1D ("h_ele_pfp_z_cosmic",         "h_ele_pfp_z_cosmic",         40, -40, 1050);
TH1D * h_ele_pfp_z_other_mixed    = new TH1D ("h_ele_pfp_z_other_mixed",    "h_ele_pfp_z_other_mixed",    40, -40, 1050);
TH1D * h_ele_pfp_z_unmatched      = new TH1D ("h_ele_pfp_z_unmatched",      "h_ele_pfp_z_unmatched",      40, -40, 1050);
TH1D * h_ele_pfp_z_intime         = new TH1D ("h_ele_pfp_z_intime",         "h_ele_pfp_z_intime",         40, -40, 1050);
TH1D * h_ele_pfp_z_data           = new TH1D ("h_ele_pfp_z_data",           "h_ele_pfp_z_data",           40, -40, 1050);

TH2D * h_post_cuts_num_tracks_showers_purity_qe     = new TH2D ("h_post_cuts_num_tracks_showers_purity_qe",    "h_post_cuts_num_tracks_showers_purity_qe",    3, 1, 4, 2, 0, 2);
TH2D * h_post_cuts_num_tracks_showers_purity_res    = new TH2D ("h_post_cuts_num_tracks_showers_purity_res",   "h_post_cuts_num_tracks_showers_purity_res",   3, 1, 4, 2, 0, 2);
TH2D * h_post_cuts_num_tracks_showers_purity_dis    = new TH2D ("h_post_cuts_num_tracks_showers_purity_dis",   "h_post_cuts_num_tracks_showers_purity_dis",   3, 1, 4, 2, 0, 2);
TH2D * h_post_cuts_num_tracks_showers_purity_coh    = new TH2D ("h_post_cuts_num_tracks_showers_purity_coh",   "h_post_cuts_num_tracks_showers_purity_coh",   3, 1, 4, 2, 0, 2);
TH2D * h_post_cuts_num_tracks_showers_purity_mec    = new TH2D ("h_post_cuts_num_tracks_showers_purity_mec",   "h_post_cuts_num_tracks_showers_purity_mec",   3, 1, 4, 2, 0, 2);
TH2D * h_post_cuts_num_tracks_showers_purity_total  = new TH2D ("h_post_cuts_num_tracks_showers_purity_total", "h_post_cuts_num_tracks_showers_purity_total", 3, 1, 4, 2, 0, 2);
TH2D * h_post_cuts_num_tracks_showers_total_total   = new TH2D ("h_post_cuts_num_tracks_showers_total_total",  "h_post_cuts_num_tracks_showers_total_total",  3, 1, 4, 2, 0, 2);

TH2D * h_ele_theta_phi_nue_cc         = new TH2D ("h_ele_theta_phi_nue_cc",         "h_ele_theta_phi_nue_cc",         14, -180, 180, 14, 0, 180);
TH2D * h_ele_theta_phi_nue_cc_out_fv  = new TH2D ("h_ele_theta_phi_nue_cc_out_fv",  "h_ele_theta_phi_nue_cc_out_fv",  14, -180, 180, 14, 0, 180);
TH2D * h_ele_theta_phi_nue_cc_mixed   = new TH2D ("h_ele_theta_phi_nue_cc_mixed",   "h_ele_theta_phi_nue_cc_mixed",   14, -180, 180, 14, 0, 180);
TH2D * h_ele_theta_phi_numu_cc        = new TH2D ("h_ele_theta_phi_numu_cc",        "h_ele_theta_phi_numu_cc",        14, -180, 180, 14, 0, 180);
TH2D * h_ele_theta_phi_numu_cc_mixed  = new TH2D ("h_ele_theta_phi_numu_cc_mixed",  "h_ele_theta_phi_numu_cc_mixed",  14, -180, 180, 14, 0, 180);
TH2D * h_ele_theta_phi_nc             = new TH2D ("h_ele_theta_phi_nc",             "h_ele_theta_phi_nc",             14, -180, 180, 14, 0, 180);
TH2D * h_ele_theta_phi_nc_pi0         = new TH2D ("h_ele_theta_phi_nc_pi0",         "h_ele_theta_phi_nc_pi0",         14, -180, 180, 14, 0, 180);
TH2D * h_ele_theta_phi_cosmic         = new TH2D ("h_ele_theta_phi_cosmic",         "h_ele_theta_phi_cosmic",         14, -180, 180, 14, 0, 180);
TH2D * h_ele_theta_phi_other_mixed    = new TH2D ("h_ele_theta_phi_other_mixed",    "h_ele_theta_phi_other_mixed",    14, -180, 180, 14, 0, 180);
TH2D * h_ele_theta_phi_unmatched      = new TH2D ("h_ele_theta_phi_unmatched",      "h_ele_theta_phi_unmatched",      14, -180, 180, 14, 0, 180);

TH2D * h_ele_eng_costheta_nue_cc         = new TH2D ("h_ele_eng_costheta_nue_cc",         "h_ele_eng_costheta_nue_cc",         10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_nue_cc_out_fv  = new TH2D ("h_ele_eng_costheta_nue_cc_out_fv",  "h_ele_eng_costheta_nue_cc_out_fv",  10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_nue_cc_mixed   = new TH2D ("h_ele_eng_costheta_nue_cc_mixed",   "h_ele_eng_costheta_nue_cc_mixed",   10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_numu_cc        = new TH2D ("h_ele_eng_costheta_numu_cc",        "h_ele_eng_costheta_numu_cc",        10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_numu_cc_mixed  = new TH2D ("h_ele_eng_costheta_numu_cc_mixed",  "h_ele_eng_costheta_numu_cc_mixed",  10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_nc             = new TH2D ("h_ele_eng_costheta_nc",             "h_ele_eng_costheta_nc",             10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_nc_pi0         = new TH2D ("h_ele_eng_costheta_nc_pi0",         "h_ele_eng_costheta_nc_pi0",         10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_cosmic         = new TH2D ("h_ele_eng_costheta_cosmic",         "h_ele_eng_costheta_cosmic",         10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_other_mixed    = new TH2D ("h_ele_eng_costheta_other_mixed",    "h_ele_eng_costheta_other_mixed",    10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_unmatched      = new TH2D ("h_ele_eng_costheta_unmatched",      "h_ele_eng_costheta_unmatched",      10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_intime         = new TH2D ("h_ele_eng_costheta_intime",         "h_ele_eng_costheta_intime",         10, 0, 2, 10, -1, 1);
TH2D * h_ele_eng_costheta_data           = new TH2D ("h_ele_eng_costheta_data",           "h_ele_eng_costheta_data",           10, 0, 2, 10, -1, 1);

TH1D * h_nue_true_theta  = new TH1D ("h_nue_true_theta", "h_nue_true_theta", 14,    0, 180);
TH1D * h_nue_true_phi    = new TH1D ("h_nue_true_phi", "h_nue_true_phi",   14, -180, 180);
TH2D * h_nue_true_theta_phi = new TH2D ("h_nue_true_theta_phi", "h_nue_true_theta_phi", 14, -180, 180, 14, 0, 180);

TH2D * h_true_reco_ele_momentum = new TH2D ("h_true_reco_ele_momentum", "h_true_reco_ele_momentum", 10, 0, 3, 10, 0, 3);
TH2D * h_true_reco_ele_costheta = new TH2D ("h_true_reco_ele_costheta", "h_true_reco_ele_costheta", 10, -1, 1, 10, -1, 1);
TH1D * h_true_num_e = new TH1D ("h_true_num_e", "h_true_num_e", 5, 0, 5);
TH2D * h_true_reco_ele_momentum_pre = new TH2D ("h_true_reco_ele_momentum_pre", "h_true_reco_ele_momentum_pre", 10, 0, 3, 10, 0, 3);
TH2D * h_true_reco_ele_costheta_pre = new TH2D ("h_true_reco_ele_costheta_pre", "h_true_reco_ele_costheta_pre", 10, -1, 1, 10, -1, 1);
TH1D * h_true_num_e_pre = new TH1D ("h_true_num_e_pre", "h_true_num_e_pre", 5, 0, 5);

TH1D * h_ele_eng_for_nue_cc         = new TH1D ("h_ele_eng_for_nue_cc",         "h_ele_eng_for_nue_cc",         10, 0, 2);
TH1D * h_ele_eng_for_nue_cc_out_fv  = new TH1D ("h_ele_eng_for_nue_cc_out_fv",  "h_ele_eng_for_nue_cc_out_fv",  10, 0, 2);
TH1D * h_ele_eng_for_nue_cc_mixed   = new TH1D ("h_ele_eng_for_nue_cc_mixed",   "h_ele_eng_for_nue_cc_mixed",   10, 0, 2);
TH1D * h_ele_eng_for_numu_cc        = new TH1D ("h_ele_eng_for_numu_cc",        "h_ele_eng_for_numu_cc",        10, 0, 2);
TH1D * h_ele_eng_for_numu_cc_mixed  = new TH1D ("h_ele_eng_for_numu_cc_mixed",  "h_ele_eng_for_numu_cc_mixed",  10, 0, 2);
TH1D * h_ele_eng_for_nc             = new TH1D ("h_ele_eng_for_nc",             "h_ele_eng_for_nc",             10, 0, 2);
TH1D * h_ele_eng_for_nc_pi0         = new TH1D ("h_ele_eng_for_nc_pi0",         "h_ele_eng_for_nc_pi0",         10, 0, 2);
TH1D * h_ele_eng_for_cosmic         = new TH1D ("h_ele_eng_for_cosmic",         "h_ele_eng_for_cosmic",         10, 0, 2);
TH1D * h_ele_eng_for_other_mixed    = new TH1D ("h_ele_eng_for_other_mixed",    "h_ele_eng_for_other_mixed",    10, 0, 2);
TH1D * h_ele_eng_for_unmatched      = new TH1D ("h_ele_eng_for_unmatched",      "h_ele_eng_for_unmatched",      10, 0, 2);
TH1D * h_ele_eng_for_intime         = new TH1D ("h_ele_eng_for_intime",         "h_ele_eng_for_intime",         10, 0, 2);
TH1D * h_ele_eng_for_data           = new TH1D ("h_ele_eng_for_data",           "h_ele_eng_for_data",           10, 0, 2);
TH1D * h_ele_eng_for_trans_nue_cc         = new TH1D ("h_ele_eng_for_trans_nue_cc",         "h_ele_eng_for_trans_nue_cc",         10, 0, 2);
TH1D * h_ele_eng_for_trans_nue_cc_out_fv  = new TH1D ("h_ele_eng_for_trans_nue_cc_out_fv",  "h_ele_eng_for_trans_nue_cc_out_fv",  10, 0, 2);
TH1D * h_ele_eng_for_trans_nue_cc_mixed   = new TH1D ("h_ele_eng_for_trans_nue_cc_mixed",   "h_ele_eng_for_trans_nue_cc_mixed",   10, 0, 2);
TH1D * h_ele_eng_for_trans_numu_cc        = new TH1D ("h_ele_eng_for_trans_numu_cc",        "h_ele_eng_for_trans_numu_cc",        10, 0, 2);
TH1D * h_ele_eng_for_trans_numu_cc_mixed  = new TH1D ("h_ele_eng_for_trans_numu_cc_mixed",  "h_ele_eng_for_trans_numu_cc_mixed",  10, 0, 2);
TH1D * h_ele_eng_for_trans_nc             = new TH1D ("h_ele_eng_for_trans_nc",             "h_ele_eng_for_trans_nc",             10, 0, 2);
TH1D * h_ele_eng_for_trans_nc_pi0         = new TH1D ("h_ele_eng_for_trans_nc_pi0",         "h_ele_eng_for_trans_nc_pi0",         10, 0, 2);
TH1D * h_ele_eng_for_trans_cosmic         = new TH1D ("h_ele_eng_for_trans_cosmic",         "h_ele_eng_for_trans_cosmic",         10, 0, 2);
TH1D * h_ele_eng_for_trans_other_mixed    = new TH1D ("h_ele_eng_for_trans_other_mixed",    "h_ele_eng_for_trans_other_mixed",    10, 0, 2);
TH1D * h_ele_eng_for_trans_unmatched      = new TH1D ("h_ele_eng_for_trans_unmatched",      "h_ele_eng_for_trans_unmatched",      10, 0, 2);
TH1D * h_ele_eng_for_trans_intime         = new TH1D ("h_ele_eng_for_trans_intime",         "h_ele_eng_for_trans_intime",         10, 0, 2);
TH1D * h_ele_eng_for_trans_data           = new TH1D ("h_ele_eng_for_trans_data",           "h_ele_eng_for_trans_data",           10, 0, 2);

TH1D * h_ele_eng_mid_nue_cc         = new TH1D ("h_ele_eng_mid_nue_cc",         "h_ele_eng_mid_nue_cc",         10, 0, 2);
TH1D * h_ele_eng_mid_nue_cc_out_fv  = new TH1D ("h_ele_eng_mid_nue_cc_out_fv",  "h_ele_eng_mid_nue_cc_out_fv",  10, 0, 2);
TH1D * h_ele_eng_mid_nue_cc_mixed   = new TH1D ("h_ele_eng_mid_nue_cc_mixed",   "h_ele_eng_mid_nue_cc_mixed",   10, 0, 2);
TH1D * h_ele_eng_mid_numu_cc        = new TH1D ("h_ele_eng_mid_numu_cc",        "h_ele_eng_mid_numu_cc",        10, 0, 2);
TH1D * h_ele_eng_mid_numu_cc_mixed  = new TH1D ("h_ele_eng_mid_numu_cc_mixed",  "h_ele_eng_mid_numu_cc_mixed",  10, 0, 2);
TH1D * h_ele_eng_mid_nc             = new TH1D ("h_ele_eng_mid_nc",             "h_ele_eng_mid_nc",             10, 0, 2);
TH1D * h_ele_eng_mid_nc_pi0         = new TH1D ("h_ele_eng_mid_nc_pi0",         "h_ele_eng_mid_nc_pi0",         10, 0, 2);
TH1D * h_ele_eng_mid_cosmic         = new TH1D ("h_ele_eng_mid_cosmic",         "h_ele_eng_mid_cosmic",         10, 0, 2);
TH1D * h_ele_eng_mid_other_mixed    = new TH1D ("h_ele_eng_mid_other_mixed",    "h_ele_eng_mid_other_mixed",    10, 0, 2);
TH1D * h_ele_eng_mid_unmatched      = new TH1D ("h_ele_eng_mid_unmatched",      "h_ele_eng_mid_unmatched",      10, 0, 2);
TH1D * h_ele_eng_mid_intime         = new TH1D ("h_ele_eng_mid_intime",         "h_ele_eng_mid_intime",         10, 0, 2);
TH1D * h_ele_eng_mid_data           = new TH1D ("h_ele_eng_mid_data",           "h_ele_eng_mid_data",           10, 0, 2);
TH1D * h_ele_eng_mid_trans_nue_cc         = new TH1D ("h_ele_eng_mid_trans_nue_cc",         "h_ele_eng_mid_trans_nue_cc",         10, 0, 2);
TH1D * h_ele_eng_mid_trans_nue_cc_out_fv  = new TH1D ("h_ele_eng_mid_trans_nue_cc_out_fv",  "h_ele_eng_mid_trans_nue_cc_out_fv",  10, 0, 2);
TH1D * h_ele_eng_mid_trans_nue_cc_mixed   = new TH1D ("h_ele_eng_mid_trans_nue_cc_mixed",   "h_ele_eng_mid_trans_nue_cc_mixed",   10, 0, 2);
TH1D * h_ele_eng_mid_trans_numu_cc        = new TH1D ("h_ele_eng_mid_trans_numu_cc",        "h_ele_eng_mid_trans_numu_cc",        10, 0, 2);
TH1D * h_ele_eng_mid_trans_numu_cc_mixed  = new TH1D ("h_ele_eng_mid_trans_numu_cc_mixed",  "h_ele_eng_mid_trans_numu_cc_mixed",  10, 0, 2);
TH1D * h_ele_eng_mid_trans_nc             = new TH1D ("h_ele_eng_mid_trans_nc",             "h_ele_eng_mid_trans_nc",             10, 0, 2);
TH1D * h_ele_eng_mid_trans_nc_pi0         = new TH1D ("h_ele_eng_mid_trans_nc_pi0",         "h_ele_eng_mid_trans_nc_pi0",         10, 0, 2);
TH1D * h_ele_eng_mid_trans_cosmic         = new TH1D ("h_ele_eng_mid_trans_cosmic",         "h_ele_eng_mid_trans_cosmic",         10, 0, 2);
TH1D * h_ele_eng_mid_trans_other_mixed    = new TH1D ("h_ele_eng_mid_trans_other_mixed",    "h_ele_eng_mid_trans_other_mixed",    10, 0, 2);
TH1D * h_ele_eng_mid_trans_unmatched      = new TH1D ("h_ele_eng_mid_trans_unmatched",      "h_ele_eng_mid_trans_unmatched",      10, 0, 2);
TH1D * h_ele_eng_mid_trans_intime         = new TH1D ("h_ele_eng_mid_trans_intime",         "h_ele_eng_mid_trans_intime",         10, 0, 2);
TH1D * h_ele_eng_mid_trans_data           = new TH1D ("h_ele_eng_mid_trans_data",           "h_ele_eng_mid_trans_data",           10, 0, 2);

TH1D * h_ele_eng_back_nue_cc         = new TH1D ("h_ele_eng_back_nue_cc",         "h_ele_eng_back_nue_cc",         10, 0, 2);
TH1D * h_ele_eng_back_nue_cc_out_fv  = new TH1D ("h_ele_eng_back_nue_cc_out_fv",  "h_ele_eng_back_nue_cc_out_fv",  10, 0, 2);
TH1D * h_ele_eng_back_nue_cc_mixed   = new TH1D ("h_ele_eng_back_nue_cc_mixed",   "h_ele_eng_back_nue_cc_mixed",   10, 0, 2);
TH1D * h_ele_eng_back_numu_cc        = new TH1D ("h_ele_eng_back_numu_cc",        "h_ele_eng_back_numu_cc",        10, 0, 2);
TH1D * h_ele_eng_back_numu_cc_mixed  = new TH1D ("h_ele_eng_back_numu_cc_mixed",  "h_ele_eng_back_numu_cc_mixed",  10, 0, 2);
TH1D * h_ele_eng_back_nc             = new TH1D ("h_ele_eng_back_nc",             "h_ele_eng_back_nc",             10, 0, 2);
TH1D * h_ele_eng_back_nc_pi0         = new TH1D ("h_ele_eng_back_nc_pi0",         "h_ele_eng_back_nc_pi0",         10, 0, 2);
TH1D * h_ele_eng_back_cosmic         = new TH1D ("h_ele_eng_back_cosmic",         "h_ele_eng_back_cosmic",         10, 0, 2);
TH1D * h_ele_eng_back_other_mixed    = new TH1D ("h_ele_eng_back_other_mixed",    "h_ele_eng_back_other_mixed",    10, 0, 2);
TH1D * h_ele_eng_back_unmatched      = new TH1D ("h_ele_eng_back_unmatched",      "h_ele_eng_back_unmatched",      10, 0, 2);
TH1D * h_ele_eng_back_intime         = new TH1D ("h_ele_eng_back_intime",         "h_ele_eng_back_intime",         10, 0, 2);
TH1D * h_ele_eng_back_data           = new TH1D ("h_ele_eng_back_data",           "h_ele_eng_back_data",           10, 0, 2);
TH1D * h_ele_eng_back_trans_nue_cc         = new TH1D ("h_ele_eng_back_trans_nue_cc",         "h_ele_eng_back_trans_nue_cc",         10, 0, 2);
TH1D * h_ele_eng_back_trans_nue_cc_out_fv  = new TH1D ("h_ele_eng_back_trans_nue_cc_out_fv",  "h_ele_eng_back_trans_nue_cc_out_fv",  10, 0, 2);
TH1D * h_ele_eng_back_trans_nue_cc_mixed   = new TH1D ("h_ele_eng_back_trans_nue_cc_mixed",   "h_ele_eng_back_trans_nue_cc_mixed",   10, 0, 2);
TH1D * h_ele_eng_back_trans_numu_cc        = new TH1D ("h_ele_eng_back_trans_numu_cc",        "h_ele_eng_back_trans_numu_cc",        10, 0, 2);
TH1D * h_ele_eng_back_trans_numu_cc_mixed  = new TH1D ("h_ele_eng_back_trans_numu_cc_mixed",  "h_ele_eng_back_trans_numu_cc_mixed",  10, 0, 2);
TH1D * h_ele_eng_back_trans_nc             = new TH1D ("h_ele_eng_back_trans_nc",             "h_ele_eng_back_trans_nc",             10, 0, 2);
TH1D * h_ele_eng_back_trans_nc_pi0         = new TH1D ("h_ele_eng_back_trans_nc_pi0",         "h_ele_eng_back_trans_nc_pi0",         10, 0, 2);
TH1D * h_ele_eng_back_trans_cosmic         = new TH1D ("h_ele_eng_back_trans_cosmic",         "h_ele_eng_back_trans_cosmic",         10, 0, 2);
TH1D * h_ele_eng_back_trans_other_mixed    = new TH1D ("h_ele_eng_back_trans_other_mixed",    "h_ele_eng_back_trans_other_mixed",    10, 0, 2);
TH1D * h_ele_eng_back_trans_unmatched      = new TH1D ("h_ele_eng_back_trans_unmatched",      "h_ele_eng_back_trans_unmatched",      10, 0, 2);
TH1D * h_ele_eng_back_trans_intime         = new TH1D ("h_ele_eng_back_trans_intime",         "h_ele_eng_back_trans_intime",         10, 0, 2);
TH1D * h_ele_eng_back_trans_data           = new TH1D ("h_ele_eng_back_trans_data",           "h_ele_eng_back_trans_data",           10, 0, 2);

TH2D * h_mc_vtx_xy_nue_cc     = new TH2D ("h_mc_vtx_xy_nue_cc",     "h_mc_vtx_xy_nue_cc",     20,    0,  260, 20, -120, 120 );
TH2D * h_mc_vtx_xz_nue_cc     = new TH2D ("h_mc_vtx_xz_nue_cc",     "h_mc_vtx_xz_nue_cc",     20,    0,  260, 20,    0, 1060);
TH2D * h_mc_vtx_yz_nue_cc     = new TH2D ("h_mc_vtx_yz_nue_cc",     "h_mc_vtx_yz_nue_cc",     20,    0, 1060, 20, -120, 120 );
TH2D * h_reco_vtx_xy_nue_cc   = new TH2D ("h_reco_vtx_xy_nue_cc",   "h_reco_vtx_xy_nue_cc",   20,    0,  260, 20, -120, 120 );
TH2D * h_reco_vtx_xz_nue_cc   = new TH2D ("h_reco_vtx_xz_nue_cc",   "h_reco_vtx_xz_nue_cc",   20,    0,  260, 20,    0, 1060);
TH2D * h_reco_vtx_yz_nue_cc   = new TH2D ("h_reco_vtx_yz_nue_cc",   "h_reco_vtx_yz_nue_cc",   20,    0, 1060, 20, -120, 120 );
TH2D * h_mc_reco_vtx_x_nue_cc = new TH2D ("h_mc_reco_vtx_x_nue_cc", "h_mc_reco_vtx_x_nue_cc", 20,    0,  260, 20,    0, 260 );
TH2D * h_mc_reco_vtx_y_nue_cc = new TH2D ("h_mc_reco_vtx_y_nue_cc", "h_mc_reco_vtx_y_nue_cc", 20, -120,  120, 20, -120, 120 );
TH2D * h_mc_reco_vtx_z_nue_cc = new TH2D ("h_mc_reco_vtx_z_nue_cc", "h_mc_reco_vtx_z_nue_cc", 20,    0, 1060, 20,    0, 1060);
TH2D * h_mc_vtx_xy_nue_cc_out_fv     = new TH2D ("h_mc_vtx_xy_nue_cc_out_fv",     "h_mc_vtx_xy_nue_cc_out_fv",     20,    0,  260, 20, -120, 120 );
TH2D * h_mc_vtx_xz_nue_cc_out_fv     = new TH2D ("h_mc_vtx_xz_nue_cc_out_fv",     "h_mc_vtx_xz_nue_cc_out_fv",     20,    0,  260, 20,    0, 1060);
TH2D * h_mc_vtx_yz_nue_cc_out_fv     = new TH2D ("h_mc_vtx_yz_nue_cc_out_fv",     "h_mc_vtx_yz_nue_cc_out_fv",     20,    0, 1060, 20, -120, 120 );
TH2D * h_reco_vtx_xy_nue_cc_out_fv   = new TH2D ("h_reco_vtx_xy_nue_cc_out_fv",   "h_reco_vtx_xy_nue_cc_out_fv",   20,    0,  260, 20, -120, 120 );
TH2D * h_reco_vtx_xz_nue_cc_out_fv   = new TH2D ("h_reco_vtx_xz_nue_cc_out_fv",   "h_reco_vtx_xz_nue_cc_out_fv",   20,    0,  260, 20,    0, 1060);
TH2D * h_reco_vtx_yz_nue_cc_out_fv   = new TH2D ("h_reco_vtx_yz_nue_cc_out_fv",   "h_reco_vtx_yz_nue_cc_out_fv",   20,    0, 1060, 20, -120, 120 );
TH2D * h_mc_reco_vtx_x_nue_cc_out_fv = new TH2D ("h_mc_reco_vtx_x_nue_cc_out_fv", "h_mc_reco_vtx_x_nue_cc_out_fv", 20,    0,  260, 20,    0, 260 );
TH2D * h_mc_reco_vtx_y_nue_cc_out_fv = new TH2D ("h_mc_reco_vtx_y_nue_cc_out_fv", "h_mc_reco_vtx_y_nue_cc_out_fv", 20, -120,  120, 20, -120, 120 );
TH2D * h_mc_reco_vtx_z_nue_cc_out_fv = new TH2D ("h_mc_reco_vtx_z_nue_cc_out_fv", "h_mc_reco_vtx_z_nue_cc_out_fv", 20,    0, 1060, 20,    0, 1060);

TH2D * h_post_cuts_num_tracks_showers_signal_total  = new TH2D ("h_post_cuts_num_tracks_showers_signal_total", "h_post_cuts_num_tracks_showers_signal_total", 3, 1, 4, 2, 0, 2);
TH2D * h_post_cuts_num_tracks_showers_bkg_total     = new TH2D ("h_post_cuts_num_tracks_showers_bkg_total",    "h_post_cuts_num_tracks_showers_bkg_total",    3, 1, 4, 2, 0, 2);


TH1D * h_low_true_momentum = new TH1D  ("h_low_true_momentum",  "h_low_true_momentum",  20, 0, 4);
TH1D * h_med_true_momentum = new TH1D  ("h_med_true_momentum",  "h_med_true_momentum",  20, 0, 4);
TH1D * h_high_true_momentum = new TH1D ("h_high_true_momentum", "h_high_true_momentum", 20, 0, 4);

}//end namespace

#endif
