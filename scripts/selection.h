#ifndef SELECTION_h
#define SELECTION_h

#include "selection_functions.h"
#include "selection_cuts.h"
#include "selection_functions_data.h"
#include "histogram_functions.h"

#include "../xsecAna/LinkDef.h"

#include <fstream>

namespace xsecSelection {

class selection {

private:

//newest values (July 2018):
// Run1 On-Beam Data:
// 3.22250e20 POT (tor101_wcut)
// 8.287403e6 EA9CNT_wcut triggers
//these values come from the sam definitions

// MC : 1.820227e21 POT
//this value comes from my scripts

// EXT: 18.732136e6 EXT triggers
// this value comes from the sam definition

//This means that POT Scaling of On-Beam : MC = 5.648
// And EXT On-Beam : EXT = 0.442416

//const double POT = 1.82027e+21;
//const double data_scale_factor = (1 / 5.648) * 1.1; //ie scale MC down by factor
//const double intime_scale_factor = 0.442416 * 5.7 * 0.9; //ie scale EXT down by factor
//const double intime_scale_factor = 0.56940408 * 5.7 * 1.1; //for first two datasets,
//factor 5.7 comes from prescale on EXT data!
//extra factor is a temporary measure for counting triggers and POT in data

//as of applying DQM and removing aNuMI running
const double POT = 1.82027e+21;
//scale MC down by this factor (i.e. should be less than 1)
// 2.369e+20 / 1.82027e+21 = 0.13
const double data_scale_factor = 0.13;

//scale EXT down by this factor
//(26749613 (ext) * 0.2064 (prescale)) = 5,521,120.1232 - per Run weighted
// 6180310 (data) / 5521120.1232 = 1.12

// per event weighted prescale is 0.51
//(26749613 (ext) * 0.51 (prescale)) = 13642302.63
// 6180310 / 13642302.63 = 0.453

//according to Zarko's new script - 5,287,179 EXT (prescale folded in already)
// 6180310 / 5287179 = 1.16

const double intime_scale_factor = 1.16;



//these are for the flux calculations
const double scaling_nue = 1.52938e-11;        //nues  / POT / cm^2
const double scaling_nue_bar = 7.77111e-12;    //anues / POT / cm^2
const double scaling = scaling_nue + scaling_nue_bar;

//these values come from GENIE file - no MEC!
// const double genie_xsec_nue = 5.63067e-39;   //cm^2
// const double genie_xsec_nue_bar = 2.0893e-39;   //cm^2
//const double genie_xsec = genie_xsec_nue + genie_xsec_nue_bar;
//this values come from GENIE file and have MEC!
const double genie_xsec_nue = 6.34569e-39; //cm2
const double genie_xsec_nue_bar = 2.24685e-39; //cm2


// older values used
// const double POT = 1.82027e21;
//
// //const double POT = 4.05982e+19;      //POT - all NuMI + cosmics
// //const double POT = 1.23206e+20; //POT - all NuMI + cosmics, bigger sample
// //const double POT = 2.90469e+21;    //POT - nue + cosmics
// const double scaling_nue = 1.52938e-11;    //nues / POT / cm^2
// const double scaling_nue_bar = 7.77111e-12;   //anues / POT / cm^2
// const double scaling = (scaling_nue + scaling_nue_bar);
//
// //since I have both nue and nue_bar as signal definition need to adjust for this
// const double genie_xsec_nue = 5.63067e-39;   //cm^2
// const double genie_xsec_nue_bar = 2.0893e-39;   //cm^2
// const double genie_xsec = genie_xsec_nue + genie_xsec_nue_bar;
//
// /*
//    3e13 POT / spills for NuMI -> 6.0675667e7 triggers for MC, or 5.9966667e5 triggers when MC scaled to data POT
//    2.571102 Million EXT spills
//    scale EXT by: 0.23323333 when MC scaled to data
//    otherwise scale EXT by: intime_scale_factor / data_scale_factor (23.599090)
//
//    For new EXT: EXT 352180
//    But looking just at samdef = 2571102 ??? why ???
//    intime_scale_factor = 0.59966667e6 / 352180
//    intime_scale_factor =
//  */
// //const double intime_scale_factor = 23.5991 * (0.0098831492);
// //const double intime_scale_factor = 23.5991 * (0.0364451);
// //const double intime_scale_factor = 23.5991 * (0.0556456);
//
// /*
//    scale via POT - 1.23e20 MC / 2.189e19 POT for full april sample
//    upscale Data by: 5.62
//    but if MC is downscaled by EXT difference, then we do 5.62 / 1.5946470
//    = 3.52
//  */
//
// //new data set has : 1.799e19 POT
// //EA9CNT: 472210
// //const double data_scale_factor = 0.0098831492;
//
// //tor101_wcut for whole samdef = 6.634e+19 POT
// //EA9CNT_wcut = 1678059
// const double intime_scale_factor = 0.655210;
// const double data_scale_factor = 0.0364451;
//
// //tor101_wcut for summing 1,2,3,4 samdef manually = 1.0129e20
// //const double data_scale_factor = 0.0556456;


const double theta_translation = 29.36 * (3.1415/180);
const double phi_translation = 8.121 * (3.1415/180);
// const double theta_translation = 0.0;
// const double phi_translation = 0.0;

double _x1;
double _x2;
double _y1;
double _y2;
double _z1;
double _z2;
double flash_pe_threshold;
double flash_time_start;
double flash_time_end;
double tolerance;
double shwr_nue_tolerance;
double trk_nue_tolerance;
double shwr_hit_threshold;
double shwr_hit_threshold_collection;
double tolerance_open_angle_min;
double tolerance_open_angle_max;
double tolerance_dedx_min;
double tolerance_dedx_max;
double dist_tolerance;
double pfp_hits_length_tolerance;
double ratio_tolerance;

public:

selection() = default;

void make_selection(
        const char * _file1,
        const char * _file2,
        const char * _file3,
        const std::vector<double> _config,
        std::vector<std::tuple<double, double, std::string> > * results_v
        );

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
TH1D * h_nue_eng_eff_den           = new TH1D("h_nue_eng_eff_den", "h_nue_eng_eff_den", 40, 0, 4);
TH1D * h_nue_eng_eff_num           = new TH1D("h_nue_eng_eff_num", "h_nue_eng_eff_num", 40, 0, 4);
TH1D * h_ele_eng_eff_den           = new TH1D("h_ele_eng_eff_den", "h_ele_eng_eff_den", 40, 0, 4);
TH1D * h_ele_eng_eff_num           = new TH1D("h_ele_eng_eff_num", "h_ele_eng_eff_num", 40, 0, 4);
TH1D * h_ele_eng_eff_num_pre_cuts  = new TH1D("h_ele_eng_eff_num_pre_cuts", "h_ele_eng_eff_num_pre_cuts", 40, 0, 4);
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

TH1D * h_flash_time        = new TH1D ("h_flash_time",        "h_flash_time",        160, 0, 20);
TH1D * h_flash_time_intime = new TH1D ("h_flash_time_intime", "h_flash_time_intime", 160, 0, 20);
TH1D * h_flash_time_data   = new TH1D ("h_flash_time_data",   "h_flash_time_data",   160, 0, 20);
//TH1D * h_flash_time_data_first_half   = new TH1D ("h_flash_time_data_first_half",    "h_flash_time_data_first_half",   180, 0, 20);
//TH1D * h_flash_time_data_second_half  = new TH1D ("h_flash_time_data_second_half",   "h_flash_time_data_second_half",  180, 0, 20);

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

TH1D * h_vtx_flash_upstream_nue_cc        = new TH1D("h_vtx_flash_upstream_nue_cc",         "h_vtx_flash_upstream_nue_cc",        40, 0, 200);
TH1D * h_vtx_flash_upstream_nue_cc_mixed  = new TH1D("h_vtx_flash_upstream_nue_cc_mixed",   "h_vtx_flash_upstream_nue_cc_mixed",  40, 0, 200);
TH1D * h_vtx_flash_upstream_nue_cc_out_fv = new TH1D("h_vtx_flash_upstream_nue_cc_out_fv",  "h_vtx_flash_upstream_nue_cc_out_fv", 40, 0, 200);
TH1D * h_vtx_flash_upstream_numu_cc       = new TH1D("h_vtx_flash_upstream_numu_cc",        "h_vtx_flash_upstream_numu_cc_mixed", 40, 0, 200);
TH1D * h_vtx_flash_upstream_nc_pi0        = new TH1D("h_vtx_flash_upstream_nc_pi0",         "h_vtx_flash_upstream_nc_pi0",        40, 0, 200);
TH1D * h_vtx_flash_upstream_cosmic        = new TH1D("h_vtx_flash_upstream_cosmic",         "h_vtx_flash_upstream_cosmic",        40, 0, 200);
TH1D * h_vtx_flash_upstream_nc            = new TH1D("h_vtx_flash_upstream_nc",             "h_vtx_flash_upstream_nc",            40, 0, 200);
TH1D * h_vtx_flash_upstream_numu_cc_mixed = new TH1D("h_vtx_flash_upstream_numu_cc_mixed",  "h_vtx_flash_upstream_numu_cc_mixed", 40, 0, 200);
TH1D * h_vtx_flash_upstream_other_mixed   = new TH1D("h_vtx_flash_upstream_other_mixed",    "h_vtx_flash_upstream_other_mixed",   40, 0, 200);
TH1D * h_vtx_flash_upstream_unmatched     = new TH1D("h_vtx_flash_upstream_unmatched",      "h_vtx_flash_upstream_unmatched",     40, 0, 200);
TH1D * h_vtx_flash_upstream_intime        = new TH1D("h_vtx_flash_upstream_intime",         "h_vtx_flash_upstream_intime",        40, 0, 200);
TH1D * h_vtx_flash_upstream_data          = new TH1D("h_vtx_flash_upstream_data",           "h_vtx_flash_upstream_data",          40, 0, 200);

TH1D * h_vtx_flash_downstream_nue_cc        = new TH1D("h_vtx_flash_downstream_nue_cc",         "h_vtx_flash_downstream_nue_cc",        40, 0, 200);
TH1D * h_vtx_flash_downstream_nue_cc_mixed  = new TH1D("h_vtx_flash_downstream_nue_cc_mixed",   "h_vtx_flash_downstream_nue_cc_mixed",  40, 0, 200);
TH1D * h_vtx_flash_downstream_nue_cc_out_fv = new TH1D("h_vtx_flash_downstream_nue_cc_out_fv",  "h_vtx_flash_downstream_nue_cc_out_fv", 40, 0, 200);
TH1D * h_vtx_flash_downstream_numu_cc       = new TH1D("h_vtx_flash_downstream_numu_cc",        "h_vtx_flash_downstream_numu_cc_mixed", 40, 0, 200);
TH1D * h_vtx_flash_downstream_nc_pi0        = new TH1D("h_vtx_flash_downstream_nc_pi0",         "h_vtx_flash_downstream_nc_pi0",        40, 0, 200);
TH1D * h_vtx_flash_downstream_cosmic        = new TH1D("h_vtx_flash_downstream_cosmic",         "h_vtx_flash_downstream_cosmic",        40, 0, 200);
TH1D * h_vtx_flash_downstream_nc            = new TH1D("h_vtx_flash_downstream_nc",             "h_vtx_flash_downstream_nc",            40, 0, 200);
TH1D * h_vtx_flash_downstream_numu_cc_mixed = new TH1D("h_vtx_flash_downstream_numu_cc_mixed",  "h_vtx_flash_downstream_numu_cc_mixed", 40, 0, 200);
TH1D * h_vtx_flash_downstream_other_mixed   = new TH1D("h_vtx_flash_downstream_other_mixed",    "h_vtx_flash_downstream_other_mixed",   40, 0, 200);
TH1D * h_vtx_flash_downstream_unmatched     = new TH1D("h_vtx_flash_downstream_unmatched",      "h_vtx_flash_downstream_unmatched",     40, 0, 200);
TH1D * h_vtx_flash_downstream_intime        = new TH1D("h_vtx_flash_downstream_intime",         "h_vtx_flash_downstream_intime",        40, 0, 200);
TH1D * h_vtx_flash_downstream_data          = new TH1D("h_vtx_flash_downstream_data",           "h_vtx_flash_downstream_data",          40, 0, 200);

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

TH1D * h_leading_momentum_nue_cut_electron       = new TH1D("h_leading_momentum_nue_cut_electron",      "h_leading_momentum_nue_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_nue_cut_photon         = new TH1D("h_leading_momentum_nue_cut_photon",        "h_leading_momentum_nue_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_nue_cut_proton         = new TH1D("h_leading_momentum_nue_cut_proton",        "h_leading_momentum_nue_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_nue_cut_pion           = new TH1D("h_leading_momentum_nue_cut_pion",          "h_leading_momentum_nue_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_nue_cut_muon           = new TH1D("h_leading_momentum_nue_cut_muon",          "h_leading_momentum_nue_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_nue_cut_kaon           = new TH1D("h_leading_momentum_nue_cut_kaon",          "h_leading_momentum_nue_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_nue_cut_neutron        = new TH1D("h_leading_momentum_nue_cut_neutron",       "h_leading_momentum_nue_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_nue_cut_mc_unmatched   = new TH1D("h_leading_momentum_nue_cut_mc_unmatched",  "h_leading_momentum_nue_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_nue_cut_ext_unmatched  = new TH1D("h_leading_momentum_nue_cut_ext_unmatched", "h_leading_momentum_nue_cut_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_fv_cut_electron       = new TH1D("h_leading_momentum_fv_cut_electron",      "h_leading_momentum_fv_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_fv_cut_photon         = new TH1D("h_leading_momentum_fv_cut_photon",        "h_leading_momentum_fv_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_fv_cut_proton         = new TH1D("h_leading_momentum_fv_cut_proton",        "h_leading_momentum_fv_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_fv_cut_pion           = new TH1D("h_leading_momentum_fv_cut_pion",          "h_leading_momentum_fv_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_fv_cut_muon           = new TH1D("h_leading_momentum_fv_cut_muon",          "h_leading_momentum_fv_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_fv_cut_kaon           = new TH1D("h_leading_momentum_fv_cut_kaon",          "h_leading_momentum_fv_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_fv_cut_neutron        = new TH1D("h_leading_momentum_fv_cut_neutron",       "h_leading_momentum_fv_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_fv_cut_mc_unmatched   = new TH1D("h_leading_momentum_fv_cut_mc_unmatched",  "h_leading_momentum_fv_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_fv_cut_ext_unmatched  = new TH1D("h_leading_momentum_fv_cut_ext_unmatched", "h_leading_momentum_fv_cut_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_flash_vtx_electron        = new TH1D("h_leading_momentum_flash_vtx_electron",      "h_leading_momentum_flash_vtx_electron",      20, 0, 4);
TH1D * h_leading_momentum_flash_vtx_photon          = new TH1D("h_leading_momentum_flash_vtx_photon",        "h_leading_momentum_flash_vtx_photon",        20, 0, 4);
TH1D * h_leading_momentum_flash_vtx_proton          = new TH1D("h_leading_momentum_flash_vtx_proton",        "h_leading_momentum_flash_vtx_proton",        20, 0, 4);
TH1D * h_leading_momentum_flash_vtx_pion            = new TH1D("h_leading_momentum_flash_vtx_pion",          "h_leading_momentum_flash_vtx_pion",          20, 0, 4);
TH1D * h_leading_momentum_flash_vtx_muon            = new TH1D("h_leading_momentum_flash_vtx_muon",          "h_leading_momentum_flash_vtx_muon",          20, 0, 4);
TH1D * h_leading_momentum_flash_vtx_kaon            = new TH1D("h_leading_momentum_flash_vtx_kaon",          "h_leading_momentum_flash_vtx_kaon",          20, 0, 4);
TH1D * h_leading_momentum_flash_vtx_neutron         = new TH1D("h_leading_momentum_flash_vtx_neutron",       "h_leading_momentum_flash_vtx_neutron",       20, 0, 4);
TH1D * h_leading_momentum_flash_vtx_mc_unmatched    = new TH1D("h_leading_momentum_flash_vtx_mc_unmatched",  "h_leading_momentum_flash_vtx_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_flash_vtx_ext_unmatched   = new TH1D("h_leading_momentum_flash_vtx_ext_unmatched", "h_leading_momentum_flash_vtx_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_shwr_vtx_electron       = new TH1D("h_leading_momentum_shwr_vtx_electron",      "h_leading_momentum_shwr_vtx_electron",      20, 0, 4);
TH1D * h_leading_momentum_shwr_vtx_photon         = new TH1D("h_leading_momentum_shwr_vtx_photon",        "h_leading_momentum_shwr_vtx_photon",        20, 0, 4);
TH1D * h_leading_momentum_shwr_vtx_proton         = new TH1D("h_leading_momentum_shwr_vtx_proton",        "h_leading_momentum_shwr_vtx_proton",        20, 0, 4);
TH1D * h_leading_momentum_shwr_vtx_pion           = new TH1D("h_leading_momentum_shwr_vtx_pion",          "h_leading_momentum_shwr_vtx_pion",          20, 0, 4);
TH1D * h_leading_momentum_shwr_vtx_muon           = new TH1D("h_leading_momentum_shwr_vtx_muon",          "h_leading_momentum_shwr_vtx_muon",          20, 0, 4);
TH1D * h_leading_momentum_shwr_vtx_kaon           = new TH1D("h_leading_momentum_shwr_vtx_kaon",          "h_leading_momentum_shwr_vtx_kaon",          20, 0, 4);
TH1D * h_leading_momentum_shwr_vtx_neutron        = new TH1D("h_leading_momentum_shwr_vtx_neutron",       "h_leading_momentum_shwr_vtx_neutron",       20, 0, 4);
TH1D * h_leading_momentum_shwr_vtx_mc_unmatched   = new TH1D("h_leading_momentum_shwr_vtx_mc_unmatched",  "h_leading_momentum_shwr_vtx_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_shwr_vtx_ext_unmatched  = new TH1D("h_leading_momentum_shwr_vtx_ext_unmatched", "h_leading_momentum_shwr_vtx_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_trk_vtx_electron       = new TH1D("h_leading_momentum_trk_vtx_electron",      "h_leading_momentum_trk_vtx_electron",      20, 0, 4);
TH1D * h_leading_momentum_trk_vtx_photon         = new TH1D("h_leading_momentum_trk_vtx_photon",        "h_leading_momentum_trk_vtx_photon",        20, 0, 4);
TH1D * h_leading_momentum_trk_vtx_proton         = new TH1D("h_leading_momentum_trk_vtx_proton",        "h_leading_momentum_trk_vtx_proton",        20, 0, 4);
TH1D * h_leading_momentum_trk_vtx_pion           = new TH1D("h_leading_momentum_trk_vtx_pion",          "h_leading_momentum_trk_vtx_pion",          20, 0, 4);
TH1D * h_leading_momentum_trk_vtx_muon           = new TH1D("h_leading_momentum_trk_vtx_muon",          "h_leading_momentum_trk_vtx_muon",          20, 0, 4);
TH1D * h_leading_momentum_trk_vtx_kaon           = new TH1D("h_leading_momentum_trk_vtx_kaon",          "h_leading_momentum_trk_vtx_kaon",          20, 0, 4);
TH1D * h_leading_momentum_trk_vtx_neutron        = new TH1D("h_leading_momentum_trk_vtx_neutron",       "h_leading_momentum_trk_vtx_neutron",       20, 0, 4);
TH1D * h_leading_momentum_trk_vtx_mc_unmatched   = new TH1D("h_leading_momentum_trk_vtx_mc_unmatched",  "h_leading_momentum_trk_vtx_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_trk_vtx_ext_unmatched  = new TH1D("h_leading_momentum_trk_vtx_ext_unmatched", "h_leading_momentum_trk_vtx_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_hit_cut_electron       = new TH1D("h_leading_momentum_hit_cut_electron",      "h_leading_momentum_hit_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_hit_cut_photon         = new TH1D("h_leading_momentum_hit_cut_photon",        "h_leading_momentum_hit_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_hit_cut_proton         = new TH1D("h_leading_momentum_hit_cut_proton",        "h_leading_momentum_hit_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_hit_cut_pion           = new TH1D("h_leading_momentum_hit_cut_pion",          "h_leading_momentum_hit_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_hit_cut_muon           = new TH1D("h_leading_momentum_hit_cut_muon",          "h_leading_momentum_hit_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_hit_cut_kaon           = new TH1D("h_leading_momentum_hit_cut_kaon",          "h_leading_momentum_hit_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_hit_cut_neutron        = new TH1D("h_leading_momentum_hit_cut_neutron",       "h_leading_momentum_hit_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_hit_cut_mc_unmatched   = new TH1D("h_leading_momentum_hit_cut_mc_unmatched",  "h_leading_momentum_hit_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_hit_cut_ext_unmatched  = new TH1D("h_leading_momentum_hit_cut_ext_unmatched", "h_leading_momentum_hit_cut_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_yhit_cut_electron       = new TH1D("h_leading_momentum_yhit_cut_electron",      "h_leading_momentum_yhit_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_yhit_cut_photon         = new TH1D("h_leading_momentum_yhit_cut_photon",        "h_leading_momentum_yhit_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_yhit_cut_proton         = new TH1D("h_leading_momentum_yhit_cut_proton",        "h_leading_momentum_yhit_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_yhit_cut_pion           = new TH1D("h_leading_momentum_yhit_cut_pion",          "h_leading_momentum_yhit_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_yhit_cut_muon           = new TH1D("h_leading_momentum_yhit_cut_muon",          "h_leading_momentum_yhit_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_yhit_cut_kaon           = new TH1D("h_leading_momentum_yhit_cut_kaon",          "h_leading_momentum_yhit_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_yhit_cut_neutron        = new TH1D("h_leading_momentum_yhit_cut_neutron",       "h_leading_momentum_yhit_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_yhit_cut_mc_unmatched   = new TH1D("h_leading_momentum_yhit_cut_mc_unmatched",  "h_leading_momentum_yhit_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_yhit_cut_ext_unmatched  = new TH1D("h_leading_momentum_yhit_cut_ext_unmatched", "h_leading_momentum_yhit_cut_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_open_angle_cut_electron       = new TH1D("h_leading_momentum_open_angle_cut_electron",      "h_leading_momentum_open_angle_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_open_angle_cut_photon         = new TH1D("h_leading_momentum_open_angle_cut_photon",        "h_leading_momentum_open_angle_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_open_angle_cut_proton         = new TH1D("h_leading_momentum_open_angle_cut_proton",        "h_leading_momentum_open_angle_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_open_angle_cut_pion           = new TH1D("h_leading_momentum_open_angle_cut_pion",          "h_leading_momentum_open_angle_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_open_angle_cut_muon           = new TH1D("h_leading_momentum_open_angle_cut_muon",          "h_leading_momentum_open_angle_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_open_angle_cut_kaon           = new TH1D("h_leading_momentum_open_angle_cut_kaon",          "h_leading_momentum_open_angle_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_open_angle_cut_neutron        = new TH1D("h_leading_momentum_open_angle_cut_neutron",       "h_leading_momentum_open_angle_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_open_angle_cut_mc_unmatched   = new TH1D("h_leading_momentum_open_angle_cut_mc_unmatched",  "h_leading_momentum_open_angle_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_open_angle_cut_ext_unmatched  = new TH1D("h_leading_momentum_open_angle_cut_ext_unmatched", "h_leading_momentum_open_angle_cut_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_dedx_cut_electron       = new TH1D("h_leading_momentum_dedx_cut_electron",      "h_leading_momentum_dedx_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_dedx_cut_photon         = new TH1D("h_leading_momentum_dedx_cut_photon",        "h_leading_momentum_dedx_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_dedx_cut_proton         = new TH1D("h_leading_momentum_dedx_cut_proton",        "h_leading_momentum_dedx_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_dedx_cut_pion           = new TH1D("h_leading_momentum_dedx_cut_pion",          "h_leading_momentum_dedx_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_dedx_cut_muon           = new TH1D("h_leading_momentum_dedx_cut_muon",          "h_leading_momentum_dedx_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_dedx_cut_kaon           = new TH1D("h_leading_momentum_dedx_cut_kaon",          "h_leading_momentum_dedx_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_dedx_cut_neutron        = new TH1D("h_leading_momentum_dedx_cut_neutron",       "h_leading_momentum_dedx_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_dedx_cut_mc_unmatched   = new TH1D("h_leading_momentum_dedx_cut_mc_unmatched",  "h_leading_momentum_dedx_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_dedx_cut_ext_unmatched  = new TH1D("h_leading_momentum_dedx_cut_ext_unmatched", "h_leading_momentum_dedx_cut_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_2shwr_cut_electron       = new TH1D("h_leading_momentum_2shwr_cut_electron",      "h_leading_momentum_2shwr_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_2shwr_cut_photon         = new TH1D("h_leading_momentum_2shwr_cut_photon",        "h_leading_momentum_2shwr_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_2shwr_cut_proton         = new TH1D("h_leading_momentum_2shwr_cut_proton",        "h_leading_momentum_2shwr_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_2shwr_cut_pion           = new TH1D("h_leading_momentum_2shwr_cut_pion",          "h_leading_momentum_2shwr_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_2shwr_cut_muon           = new TH1D("h_leading_momentum_2shwr_cut_muon",          "h_leading_momentum_2shwr_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_2shwr_cut_kaon           = new TH1D("h_leading_momentum_2shwr_cut_kaon",          "h_leading_momentum_2shwr_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_2shwr_cut_neutron        = new TH1D("h_leading_momentum_2shwr_cut_neutron",       "h_leading_momentum_2shwr_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_2shwr_cut_mc_unmatched   = new TH1D("h_leading_momentum_2shwr_cut_mc_unmatched",  "h_leading_momentum_2shwr_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_2shwr_cut_ext_unmatched  = new TH1D("h_leading_momentum_2shwr_cut_ext_unmatched", "h_leading_momentum_2shwr_cut_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_hit_length_cut_electron       = new TH1D("h_leading_momentum_hit_length_cut_electron",      "h_leading_momentum_hit_length_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_hit_length_cut_photon         = new TH1D("h_leading_momentum_hit_length_cut_photon",        "h_leading_momentum_hit_length_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_hit_length_cut_proton         = new TH1D("h_leading_momentum_hit_length_cut_proton",        "h_leading_momentum_hit_length_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_hit_length_cut_pion           = new TH1D("h_leading_momentum_hit_length_cut_pion",          "h_leading_momentum_hit_length_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_hit_length_cut_muon           = new TH1D("h_leading_momentum_hit_length_cut_muon",          "h_leading_momentum_hit_length_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_hit_length_cut_kaon           = new TH1D("h_leading_momentum_hit_length_cut_kaon",          "h_leading_momentum_hit_length_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_hit_length_cut_neutron        = new TH1D("h_leading_momentum_hit_length_cut_neutron",       "h_leading_momentum_hit_length_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_hit_length_cut_mc_unmatched   = new TH1D("h_leading_momentum_hit_length_cut_mc_unmatched",  "h_leading_momentum_hit_length_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_hit_length_cut_ext_unmatched  = new TH1D("h_leading_momentum_hit_length_cut_ext_unmatched", "h_leading_momentum_hit_length_cut_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_length_ratio_cut_electron       = new TH1D("h_leading_momentum_length_ratio_cut_electron",      "h_leading_momentum_length_ratio_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_length_ratio_cut_photon         = new TH1D("h_leading_momentum_length_ratio_cut_photon",        "h_leading_momentum_length_ratio_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_length_ratio_cut_proton         = new TH1D("h_leading_momentum_length_ratio_cut_proton",        "h_leading_momentum_length_ratio_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_length_ratio_cut_pion           = new TH1D("h_leading_momentum_length_ratio_cut_pion",          "h_leading_momentum_length_ratio_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_length_ratio_cut_muon           = new TH1D("h_leading_momentum_length_ratio_cut_muon",          "h_leading_momentum_length_ratio_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_length_ratio_cut_kaon           = new TH1D("h_leading_momentum_length_ratio_cut_kaon",          "h_leading_momentum_length_ratio_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_length_ratio_cut_neutron        = new TH1D("h_leading_momentum_length_ratio_cut_neutron",       "h_leading_momentum_length_ratio_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_length_ratio_cut_mc_unmatched   = new TH1D("h_leading_momentum_length_ratio_cut_mc_unmatched",  "h_leading_momentum_length_ratio_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_length_ratio_cut_ext_unmatched  = new TH1D("h_leading_momentum_length_ratio_cut_ext_unmatched", "h_leading_momentum_length_ratio_cut_ext_unmatched", 20, 0, 4);

TH1D * h_leading_momentum_containment_cut_electron       = new TH1D("h_leading_momentum_containment_cut_electron",      "h_leading_momentum_containment_cut_electron",      20, 0, 4);
TH1D * h_leading_momentum_containment_cut_photon         = new TH1D("h_leading_momentum_containment_cut_photon",        "h_leading_momentum_containment_cut_photon",        20, 0, 4);
TH1D * h_leading_momentum_containment_cut_proton         = new TH1D("h_leading_momentum_containment_cut_proton",        "h_leading_momentum_containment_cut_proton",        20, 0, 4);
TH1D * h_leading_momentum_containment_cut_pion           = new TH1D("h_leading_momentum_containment_cut_pion",          "h_leading_momentum_containment_cut_pion",          20, 0, 4);
TH1D * h_leading_momentum_containment_cut_muon           = new TH1D("h_leading_momentum_containment_cut_muon",          "h_leading_momentum_containment_cut_muon",          20, 0, 4);
TH1D * h_leading_momentum_containment_cut_kaon           = new TH1D("h_leading_momentum_containment_cut_kaon",          "h_leading_momentum_containment_cut_kaon",          20, 0, 4);
TH1D * h_leading_momentum_containment_cut_neutron        = new TH1D("h_leading_momentum_containment_cut_neutron",       "h_leading_momentum_containment_cut_neutron",       20, 0, 4);
TH1D * h_leading_momentum_containment_cut_mc_unmatched   = new TH1D("h_leading_momentum_containment_cut_mc_unmatched",  "h_leading_momentum_containment_cut_mc_unmatched",  20, 0, 4);
TH1D * h_leading_momentum_containment_cut_ext_unmatched  = new TH1D("h_leading_momentum_containment_cut_ext_unmatched", "h_leading_momentum_containment_cut_ext_unmatched", 20, 0, 4);

TH1D * h_dedx_cuts_electron       = new TH1D("h_dedx_cuts_electron",      "h_dedx_cuts_electron",      20, 0, 10);
TH1D * h_dedx_cuts_photon         = new TH1D("h_dedx_cuts_photon",        "h_dedx_cuts_photon",        20, 0, 10);
TH1D * h_dedx_cuts_proton         = new TH1D("h_dedx_cuts_proton",        "h_dedx_cuts_proton",        20, 0, 10);
TH1D * h_dedx_cuts_pion           = new TH1D("h_dedx_cuts_pion",          "h_dedx_cuts_pion",          20, 0, 10);
TH1D * h_dedx_cuts_muon           = new TH1D("h_dedx_cuts_muon",          "h_dedx_cuts_muon",          20, 0, 10);
TH1D * h_dedx_cuts_kaon           = new TH1D("h_dedx_cuts_kaon",          "h_dedx_cuts_kaon",          20, 0, 10);
TH1D * h_dedx_cuts_neutron        = new TH1D("h_dedx_cuts_neutron",       "h_dedx_cuts_neutron",       20, 0, 10);
TH1D * h_dedx_cuts_mc_unmatched   = new TH1D("h_dedx_cuts_mc_unmatched",  "h_dedx_cuts_mc_unmatched",  20, 0, 10);
TH1D * h_dedx_cuts_ext_unmatched  = new TH1D("h_dedx_cuts_ext_unmatched", "h_dedx_cuts_ext_unmatched", 20, 0, 10);

TH1D * h_dedx_cuts_last_electron       = new TH1D("h_dedx_cuts_last_electron",      "h_dedx_cuts_last_electron",      20, 0, 10);
TH1D * h_dedx_cuts_last_photon         = new TH1D("h_dedx_cuts_last_photon",        "h_dedx_cuts_last_photon",        20, 0, 10);
TH1D * h_dedx_cuts_last_proton         = new TH1D("h_dedx_cuts_last_proton",        "h_dedx_cuts_last_proton",        20, 0, 10);
TH1D * h_dedx_cuts_last_pion           = new TH1D("h_dedx_cuts_last_pion",          "h_dedx_cuts_last_pion",          20, 0, 10);
TH1D * h_dedx_cuts_last_muon           = new TH1D("h_dedx_cuts_last_muon",          "h_dedx_cuts_last_muon",          20, 0, 10);
TH1D * h_dedx_cuts_last_kaon           = new TH1D("h_dedx_cuts_last_kaon",          "h_dedx_cuts_last_kaon",          20, 0, 10);
TH1D * h_dedx_cuts_last_neutron        = new TH1D("h_dedx_cuts_last_neutron",       "h_dedx_cuts_last_neutron",       20, 0, 10);
TH1D * h_dedx_cuts_last_mc_unmatched   = new TH1D("h_dedx_cuts_last_mc_unmatched",  "h_dedx_cuts_last_mc_unmatched",  20, 0, 10);
TH1D * h_dedx_cuts_last_ext_unmatched  = new TH1D("h_dedx_cuts_last_ext_unmatched", "h_dedx_cuts_last_ext_unmatched", 20, 0, 10);

TH2D * h_dedx_cuts_hits_electron       = new TH2D("h_dedx_cuts_hits_electron",     "h_dedx_cuts_hits_electron",     20, 0, 10, 20, 0, 1600);
TH2D * h_dedx_cuts_hits_photon         = new TH2D("h_dedx_cuts_hits_photon",       "h_dedx_cuts_hits_photon",       20, 0, 10, 20, 0, 1600);
TH2D * h_dedx_cuts_hits_proton         = new TH2D("h_dedx_cuts_hits_proton",       "h_dedx_cuts_hits_proton",       20, 0, 10, 20, 0, 1600);
TH2D * h_dedx_cuts_hits_pion           = new TH2D("h_dedx_cuts_hits_pion",         "h_dedx_cuts_hits_pion",         20, 0, 10, 20, 0, 1600);
TH2D * h_dedx_cuts_hits_muon           = new TH2D("h_dedx_cuts_hits_muon",         "h_dedx_cuts_hits_muon",         20, 0, 10, 20, 0, 1600);
TH2D * h_dedx_cuts_hits_kaon           = new TH2D("h_dedx_cuts_hits_kaon",         "h_dedx_cuts_hits_kaon",         20, 0, 10, 20, 0, 1600);
TH2D * h_dedx_cuts_hits_neutron        = new TH2D("h_dedx_cuts_hits_neutron",      "h_dedx_cuts_hits_neutron",      20, 0, 10, 20, 0, 1600);
TH2D * h_dedx_cuts_hits_mc_unmatched   = new TH2D("h_dedx_cuts_hits_mc_unmatched", "h_dedx_cuts_hits_mc_unmatched", 20, 0, 10, 20, 0, 1600);

TH2D * h_dedx_cuts_collection_hits_electron       = new TH2D("h_dedx_cuts_collection_hits_electron",     "h_dedx_cuts_collection_hits_electron",     20, 0, 10, 20, 0, 1000);
TH2D * h_dedx_cuts_collection_hits_photon         = new TH2D("h_dedx_cuts_collection_hits_photon",       "h_dedx_cuts_collection_hits_photon",       20, 0, 10, 20, 0, 1000);
TH2D * h_dedx_cuts_collection_hits_proton         = new TH2D("h_dedx_cuts_collection_hits_proton",       "h_dedx_cuts_collection_hits_proton",       20, 0, 10, 20, 0, 1000);
TH2D * h_dedx_cuts_collection_hits_pion           = new TH2D("h_dedx_cuts_collection_hits_pion",         "h_dedx_cuts_collection_hits_pion",         20, 0, 10, 20, 0, 1000);
TH2D * h_dedx_cuts_collection_hits_muon           = new TH2D("h_dedx_cuts_collection_hits_muon",         "h_dedx_cuts_collection_hits_muon",         20, 0, 10, 20, 0, 1000);
TH2D * h_dedx_cuts_collection_hits_kaon           = new TH2D("h_dedx_cuts_collection_hits_kaon",         "h_dedx_cuts_collection_hits_kaon",         20, 0, 10, 20, 0, 1000);
TH2D * h_dedx_cuts_collection_hits_neutron        = new TH2D("h_dedx_cuts_collection_hits_neutron",      "h_dedx_cuts_collection_hits_neutron",      20, 0, 10, 20, 0, 1000);
TH2D * h_dedx_cuts_collection_hits_mc_unmatched   = new TH2D("h_dedx_cuts_collection_hits_mc_unmatched", "h_dedx_cuts_collection_hits_mc_unmatched", 20, 0, 10, 20, 0, 1000);


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

TH1D * h_ele_cos_theta_last_nue_cc         = new TH1D ("h_ele_cos_theta_last_nue_cc",         "h_ele_cos_theta_last_nue_cc",         16, -1, 1);
TH1D * h_ele_cos_theta_last_nue_cc_out_fv  = new TH1D ("h_ele_cos_theta_last_nue_cc_out_fv",  "h_ele_cos_theta_last_nue_cc_out_fv",  16, -1, 1);
TH1D * h_ele_cos_theta_last_nue_cc_mixed   = new TH1D ("h_ele_cos_theta_last_nue_cc_mixed",   "h_ele_cos_theta_last_nue_cc_mixed",   16, -1, 1);
TH1D * h_ele_cos_theta_last_numu_cc        = new TH1D ("h_ele_cos_theta_last_numu_cc",        "h_ele_cos_theta_last_numu_cc",        16, -1, 1);
TH1D * h_ele_cos_theta_last_numu_cc_mixed  = new TH1D ("h_ele_cos_theta_last_numu_cc_mixed",  "h_ele_cos_theta_last_numu_cc_mixed",  16, -1, 1);
TH1D * h_ele_cos_theta_last_nc             = new TH1D ("h_ele_cos_theta_last_nc",             "h_ele_cos_theta_last_nc",             16, -1, 1);
TH1D * h_ele_cos_theta_last_nc_pi0         = new TH1D ("h_ele_cos_theta_last_nc_pi0",         "h_ele_cos_theta_last_nc_pi0",         16, -1, 1);
TH1D * h_ele_cos_theta_last_cosmic         = new TH1D ("h_ele_cos_theta_last_cosmic",         "h_ele_cos_theta_last_cosmic",         16, -1, 1);
TH1D * h_ele_cos_theta_last_other_mixed    = new TH1D ("h_ele_cos_theta_last_other_mixed",    "h_ele_cos_theta_last_other_mixed",    16, -1, 1);
TH1D * h_ele_cos_theta_last_unmatched      = new TH1D ("h_ele_cos_theta_last_unmatched",      "h_ele_cos_theta_last_unmatched",      16, -1, 1);
TH1D * h_ele_cos_theta_last_intime         = new TH1D ("h_ele_cos_theta_last_intime",         "h_ele_cos_theta_last_intime",         16, -1, 1);
TH1D * h_ele_cos_theta_last_data           = new TH1D ("h_ele_cos_theta_last_data",           "h_ele_cos_theta_last_data",           16, -1, 1);

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

TH1D * h_ele_pfp_momentum_no_track_nue_cc         = new TH1D ("h_ele_pfp_momentum_no_track_nue_cc",         "h_ele_pfp_momentum_no_track_nue_cc",         10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_nue_cc_out_fv  = new TH1D ("h_ele_pfp_momentum_no_track_nue_cc_out_fv",  "h_ele_pfp_momentum_no_track_nue_cc_out_fv",  10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_nue_cc_mixed   = new TH1D ("h_ele_pfp_momentum_no_track_nue_cc_mixed",   "h_ele_pfp_momentum_no_track_nue_cc_mixed",   10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_numu_cc        = new TH1D ("h_ele_pfp_momentum_no_track_numu_cc",        "h_ele_pfp_momentum_no_track_numu_cc",        10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_numu_cc_mixed  = new TH1D ("h_ele_pfp_momentum_no_track_numu_cc_mixed",  "h_ele_pfp_momentum_no_track_numu_cc_mixed",  10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_nc             = new TH1D ("h_ele_pfp_momentum_no_track_nc",             "h_ele_pfp_momentum_no_track_nc",             10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_nc_pi0         = new TH1D ("h_ele_pfp_momentum_no_track_nc_pi0",         "h_ele_pfp_momentum_no_track_nc_pi0",         10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_cosmic         = new TH1D ("h_ele_pfp_momentum_no_track_cosmic",         "h_ele_pfp_momentum_no_track_cosmic",         10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_other_mixed    = new TH1D ("h_ele_pfp_momentum_no_track_other_mixed",    "h_ele_pfp_momentum_no_track_other_mixed",    10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_unmatched      = new TH1D ("h_ele_pfp_momentum_no_track_unmatched",      "h_ele_pfp_momentum_no_track_unmatched",      10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_intime         = new TH1D ("h_ele_pfp_momentum_no_track_intime",         "h_ele_pfp_momentum_no_track_intime",         10, 0, 2);
TH1D * h_ele_pfp_momentum_no_track_data           = new TH1D ("h_ele_pfp_momentum_no_track_data",           "h_ele_pfp_momentum_no_track_data",           10, 0, 2);

TH1D * h_ele_pfp_momentum_has_track_nue_cc         = new TH1D ("h_ele_pfp_momentum_has_track_nue_cc",         "h_ele_pfp_momentum_has_track_nue_cc",         10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_nue_cc_out_fv  = new TH1D ("h_ele_pfp_momentum_has_track_nue_cc_out_fv",  "h_ele_pfp_momentum_has_track_nue_cc_out_fv",  10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_nue_cc_mixed   = new TH1D ("h_ele_pfp_momentum_has_track_nue_cc_mixed",   "h_ele_pfp_momentum_has_track_nue_cc_mixed",   10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_numu_cc        = new TH1D ("h_ele_pfp_momentum_has_track_numu_cc",        "h_ele_pfp_momentum_has_track_numu_cc",        10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_numu_cc_mixed  = new TH1D ("h_ele_pfp_momentum_has_track_numu_cc_mixed",  "h_ele_pfp_momentum_has_track_numu_cc_mixed",  10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_nc             = new TH1D ("h_ele_pfp_momentum_has_track_nc",             "h_ele_pfp_momentum_has_track_nc",             10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_nc_pi0         = new TH1D ("h_ele_pfp_momentum_has_track_nc_pi0",         "h_ele_pfp_momentum_has_track_nc_pi0",         10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_cosmic         = new TH1D ("h_ele_pfp_momentum_has_track_cosmic",         "h_ele_pfp_momentum_has_track_cosmic",         10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_other_mixed    = new TH1D ("h_ele_pfp_momentum_has_track_other_mixed",    "h_ele_pfp_momentum_has_track_other_mixed",    10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_unmatched      = new TH1D ("h_ele_pfp_momentum_has_track_unmatched",      "h_ele_pfp_momentum_has_track_unmatched",      10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_intime         = new TH1D ("h_ele_pfp_momentum_has_track_intime",         "h_ele_pfp_momentum_has_track_intime",         10, 0, 2);
TH1D * h_ele_pfp_momentum_has_track_data           = new TH1D ("h_ele_pfp_momentum_has_track_data",           "h_ele_pfp_momentum_has_track_data",           10, 0, 2);

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

TH1D * h_any_pfp_x_nue_cc         = new TH1D ("h_any_pfp_x_nue_cc",         "h_any_pfp_x_nue_cc",         20, -10, 270);
TH1D * h_any_pfp_x_nue_cc_out_fv  = new TH1D ("h_any_pfp_x_nue_cc_out_fv",  "h_any_pfp_x_nue_cc_out_fv",  20, -10, 270);
TH1D * h_any_pfp_x_nue_cc_mixed   = new TH1D ("h_any_pfp_x_nue_cc_mixed",   "h_any_pfp_x_nue_cc_mixed",   20, -10, 270);
TH1D * h_any_pfp_x_numu_cc        = new TH1D ("h_any_pfp_x_numu_cc",        "h_any_pfp_x_numu_cc",        20, -10, 270);
TH1D * h_any_pfp_x_numu_cc_mixed  = new TH1D ("h_any_pfp_x_numu_cc_mixed",  "h_any_pfp_x_numu_cc_mixed",  20, -10, 270);
TH1D * h_any_pfp_x_nc             = new TH1D ("h_any_pfp_x_nc",             "h_any_pfp_x_nc",             20, -10, 270);
TH1D * h_any_pfp_x_nc_pi0         = new TH1D ("h_any_pfp_x_nc_pi0",         "h_any_pfp_x_nc_pi0",         20, -10, 270);
TH1D * h_any_pfp_x_cosmic         = new TH1D ("h_any_pfp_x_cosmic",         "h_any_pfp_x_cosmic",         20, -10, 270);
TH1D * h_any_pfp_x_other_mixed    = new TH1D ("h_any_pfp_x_other_mixed",    "h_any_pfp_x_other_mixed",    20, -10, 270);
TH1D * h_any_pfp_x_unmatched      = new TH1D ("h_any_pfp_x_unmatched",      "h_any_pfp_x_unmatched",      20, -10, 270);
TH1D * h_any_pfp_x_intime         = new TH1D ("h_any_pfp_x_intime",         "h_any_pfp_x_intime",         20, -10, 270);
TH1D * h_any_pfp_x_data           = new TH1D ("h_any_pfp_x_data",           "h_any_pfp_x_data",           20, -10, 270);

TH1D * h_any_pfp_y_nue_cc         = new TH1D ("h_any_pfp_y_nue_cc",         "h_any_pfp_y_nue_cc",         20, -120, 120);
TH1D * h_any_pfp_y_nue_cc_out_fv  = new TH1D ("h_any_pfp_y_nue_cc_out_fv",  "h_any_pfp_y_nue_cc_out_fv",  20, -120, 120);
TH1D * h_any_pfp_y_nue_cc_mixed   = new TH1D ("h_any_pfp_y_nue_cc_mixed",   "h_any_pfp_y_nue_cc_mixed",   20, -120, 120);
TH1D * h_any_pfp_y_numu_cc        = new TH1D ("h_any_pfp_y_numu_cc",        "h_any_pfp_y_numu_cc",        20, -120, 120);
TH1D * h_any_pfp_y_numu_cc_mixed  = new TH1D ("h_any_pfp_y_numu_cc_mixed",  "h_any_pfp_y_numu_cc_mixed",  20, -120, 120);
TH1D * h_any_pfp_y_nc             = new TH1D ("h_any_pfp_y_nc",             "h_any_pfp_y_nc",             20, -120, 120);
TH1D * h_any_pfp_y_nc_pi0         = new TH1D ("h_any_pfp_y_nc_pi0",         "h_any_pfp_y_nc_pi0",         20, -120, 120);
TH1D * h_any_pfp_y_cosmic         = new TH1D ("h_any_pfp_y_cosmic",         "h_any_pfp_y_cosmic",         20, -120, 120);
TH1D * h_any_pfp_y_other_mixed    = new TH1D ("h_any_pfp_y_other_mixed",    "h_any_pfp_y_other_mixed",    20, -120, 120);
TH1D * h_any_pfp_y_unmatched      = new TH1D ("h_any_pfp_y_unmatched",      "h_any_pfp_y_unmatched",      20, -120, 120);
TH1D * h_any_pfp_y_intime         = new TH1D ("h_any_pfp_y_intime",         "h_any_pfp_y_intime",         20, -120, 120);
TH1D * h_any_pfp_y_data           = new TH1D ("h_any_pfp_y_data",           "h_any_pfp_y_data",           20, -120, 120);

TH1D * h_any_pfp_z_nue_cc         = new TH1D ("h_any_pfp_z_nue_cc",         "h_any_pfp_z_nue_cc",         40, -40, 1050);
TH1D * h_any_pfp_z_nue_cc_out_fv  = new TH1D ("h_any_pfp_z_nue_cc_out_fv",  "h_any_pfp_z_nue_cc_out_fv",  40, -40, 1050);
TH1D * h_any_pfp_z_nue_cc_mixed   = new TH1D ("h_any_pfp_z_nue_cc_mixed",   "h_any_pfp_z_nue_cc_mixed",   40, -40, 1050);
TH1D * h_any_pfp_z_numu_cc        = new TH1D ("h_any_pfp_z_numu_cc",        "h_any_pfp_z_numu_cc",        40, -40, 1050);
TH1D * h_any_pfp_z_numu_cc_mixed  = new TH1D ("h_any_pfp_z_numu_cc_mixed",  "h_any_pfp_z_numu_cc_mixed",  40, -40, 1050);
TH1D * h_any_pfp_z_nc             = new TH1D ("h_any_pfp_z_nc",             "h_any_pfp_z_nc",             40, -40, 1050);
TH1D * h_any_pfp_z_nc_pi0         = new TH1D ("h_any_pfp_z_nc_pi0",         "h_any_pfp_z_nc_pi0",         40, -40, 1050);
TH1D * h_any_pfp_z_cosmic         = new TH1D ("h_any_pfp_z_cosmic",         "h_any_pfp_z_cosmic",         40, -40, 1050);
TH1D * h_any_pfp_z_other_mixed    = new TH1D ("h_any_pfp_z_other_mixed",    "h_any_pfp_z_other_mixed",    40, -40, 1050);
TH1D * h_any_pfp_z_unmatched      = new TH1D ("h_any_pfp_z_unmatched",      "h_any_pfp_z_unmatched",      40, -40, 1050);
TH1D * h_any_pfp_z_intime         = new TH1D ("h_any_pfp_z_intime",         "h_any_pfp_z_intime",         40, -40, 1050);
TH1D * h_any_pfp_z_data           = new TH1D ("h_any_pfp_z_data",           "h_any_pfp_z_data",           40, -40, 1050);

TH1D * h_any_pfp_x_last_nue_cc         = new TH1D ("h_any_pfp_x_last_nue_cc",         "h_any_pfp_x_last_nue_cc",         12, -10, 270);
TH1D * h_any_pfp_x_last_nue_cc_out_fv  = new TH1D ("h_any_pfp_x_last_nue_cc_out_fv",  "h_any_pfp_x_last_nue_cc_out_fv",  12, -10, 270);
TH1D * h_any_pfp_x_last_nue_cc_mixed   = new TH1D ("h_any_pfp_x_last_nue_cc_mixed",   "h_any_pfp_x_last_nue_cc_mixed",   12, -10, 270);
TH1D * h_any_pfp_x_last_numu_cc        = new TH1D ("h_any_pfp_x_last_numu_cc",        "h_any_pfp_x_last_numu_cc",        12, -10, 270);
TH1D * h_any_pfp_x_last_numu_cc_mixed  = new TH1D ("h_any_pfp_x_last_numu_cc_mixed",  "h_any_pfp_x_last_numu_cc_mixed",  12, -10, 270);
TH1D * h_any_pfp_x_last_nc             = new TH1D ("h_any_pfp_x_last_nc",             "h_any_pfp_x_last_nc",             12, -10, 270);
TH1D * h_any_pfp_x_last_nc_pi0         = new TH1D ("h_any_pfp_x_last_nc_pi0",         "h_any_pfp_x_last_nc_pi0",         12, -10, 270);
TH1D * h_any_pfp_x_last_cosmic         = new TH1D ("h_any_pfp_x_last_cosmic",         "h_any_pfp_x_last_cosmic",         12, -10, 270);
TH1D * h_any_pfp_x_last_other_mixed    = new TH1D ("h_any_pfp_x_last_other_mixed",    "h_any_pfp_x_last_other_mixed",    12, -10, 270);
TH1D * h_any_pfp_x_last_unmatched      = new TH1D ("h_any_pfp_x_last_unmatched",      "h_any_pfp_x_last_unmatched",      12, -10, 270);
TH1D * h_any_pfp_x_last_intime         = new TH1D ("h_any_pfp_x_last_intime",         "h_any_pfp_x_last_intime",         12, -10, 270);
TH1D * h_any_pfp_x_last_data           = new TH1D ("h_any_pfp_x_last_data",           "h_any_pfp_x_last_data",           12, -10, 270);

TH1D * h_any_pfp_y_last_nue_cc         = new TH1D ("h_any_pfp_y_last_nue_cc",         "h_any_pfp_y_last_nue_cc",         12, -120, 120);
TH1D * h_any_pfp_y_last_nue_cc_out_fv  = new TH1D ("h_any_pfp_y_last_nue_cc_out_fv",  "h_any_pfp_y_last_nue_cc_out_fv",  12, -120, 120);
TH1D * h_any_pfp_y_last_nue_cc_mixed   = new TH1D ("h_any_pfp_y_last_nue_cc_mixed",   "h_any_pfp_y_last_nue_cc_mixed",   12, -120, 120);
TH1D * h_any_pfp_y_last_numu_cc        = new TH1D ("h_any_pfp_y_last_numu_cc",        "h_any_pfp_y_last_numu_cc",        12, -120, 120);
TH1D * h_any_pfp_y_last_numu_cc_mixed  = new TH1D ("h_any_pfp_y_last_numu_cc_mixed",  "h_any_pfp_y_last_numu_cc_mixed",  12, -120, 120);
TH1D * h_any_pfp_y_last_nc             = new TH1D ("h_any_pfp_y_last_nc",             "h_any_pfp_y_last_nc",             12, -120, 120);
TH1D * h_any_pfp_y_last_nc_pi0         = new TH1D ("h_any_pfp_y_last_nc_pi0",         "h_any_pfp_y_last_nc_pi0",         12, -120, 120);
TH1D * h_any_pfp_y_last_cosmic         = new TH1D ("h_any_pfp_y_last_cosmic",         "h_any_pfp_y_last_cosmic",         12, -120, 120);
TH1D * h_any_pfp_y_last_other_mixed    = new TH1D ("h_any_pfp_y_last_other_mixed",    "h_any_pfp_y_last_other_mixed",    12, -120, 120);
TH1D * h_any_pfp_y_last_unmatched      = new TH1D ("h_any_pfp_y_last_unmatched",      "h_any_pfp_y_last_unmatched",      12, -120, 120);
TH1D * h_any_pfp_y_last_intime         = new TH1D ("h_any_pfp_y_last_intime",         "h_any_pfp_y_last_intime",         12, -120, 120);
TH1D * h_any_pfp_y_last_data           = new TH1D ("h_any_pfp_y_last_data",           "h_any_pfp_y_last_data",           12, -120, 120);

TH1D * h_any_pfp_z_last_nue_cc         = new TH1D ("h_any_pfp_z_last_nue_cc",         "h_any_pfp_z_last_nue_cc",         12, -40, 1050);
TH1D * h_any_pfp_z_last_nue_cc_out_fv  = new TH1D ("h_any_pfp_z_last_nue_cc_out_fv",  "h_any_pfp_z_last_nue_cc_out_fv",  12, -40, 1050);
TH1D * h_any_pfp_z_last_nue_cc_mixed   = new TH1D ("h_any_pfp_z_last_nue_cc_mixed",   "h_any_pfp_z_last_nue_cc_mixed",   12, -40, 1050);
TH1D * h_any_pfp_z_last_numu_cc        = new TH1D ("h_any_pfp_z_last_numu_cc",        "h_any_pfp_z_last_numu_cc",        12, -40, 1050);
TH1D * h_any_pfp_z_last_numu_cc_mixed  = new TH1D ("h_any_pfp_z_last_numu_cc_mixed",  "h_any_pfp_z_last_numu_cc_mixed",  12, -40, 1050);
TH1D * h_any_pfp_z_last_nc             = new TH1D ("h_any_pfp_z_last_nc",             "h_any_pfp_z_last_nc",             12, -40, 1050);
TH1D * h_any_pfp_z_last_nc_pi0         = new TH1D ("h_any_pfp_z_last_nc_pi0",         "h_any_pfp_z_last_nc_pi0",         12, -40, 1050);
TH1D * h_any_pfp_z_last_cosmic         = new TH1D ("h_any_pfp_z_last_cosmic",         "h_any_pfp_z_last_cosmic",         12, -40, 1050);
TH1D * h_any_pfp_z_last_other_mixed    = new TH1D ("h_any_pfp_z_last_other_mixed",    "h_any_pfp_z_last_other_mixed",    12, -40, 1050);
TH1D * h_any_pfp_z_last_unmatched      = new TH1D ("h_any_pfp_z_last_unmatched",      "h_any_pfp_z_last_unmatched",      12, -40, 1050);
TH1D * h_any_pfp_z_last_intime         = new TH1D ("h_any_pfp_z_last_intime",         "h_any_pfp_z_last_intime",         12, -40, 1050);
TH1D * h_any_pfp_z_last_data           = new TH1D ("h_any_pfp_z_last_data",           "h_any_pfp_z_last_data",           12, -40, 1050);

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
TH2D * h_ele_theta_phi_intime         = new TH2D ("h_ele_theta_phi_intime",         "h_ele_theta_phi_intime",         14, -180, 180, 14, 0, 180);

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
TH2D * h_nue_true_energy_theta = new TH2D("h_nue_true_energy_theta", "h_nue_true_energy_theta", 20, 0, 10, 14, 0, 180);
TH2D * h_all_true_energy_theta = new TH2D("h_all_true_energy_theta", "h_all_true_energy_theta", 20, 0, 10, 14, 0, 180);
TH2D * h_nue_true_energy_phi   = new TH2D("h_nue_true_energy_phi",   "h_nue_true_energy_phi",   20, 0, 10, 14, -180, 180);
TH2D * h_ele_true_energy_theta = new TH2D("h_ele_true_energy_theta", "h_ele_true_energy_theta", 20, 0, 10, 14, 0, 180);
TH2D * h_ele_true_energy_phi   = new TH2D("h_ele_true_energy_phi",   "h_ele_true_energy_phi",   20, 0, 10, 14, -180, 180);

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

TH2 * h_pfp_zy_vtx_nue_cc = new TH2D ("h_pfp_zy_vtx_nue_cc", "h_pfp_zy_vtx_nue_cc", 40, 0, 1038, 40, -117, 117);
TH2 * h_pfp_zy_vtx_all    = new TH2D ("h_pfp_zy_vtx_all",    "h_pfp_zy_vtx_all",    40, 0, 1038, 40, -117, 117);
TH2 * h_pfp_zy_vtx_ext    = new TH2D ("h_pfp_zy_vtx_ext",    "h_pfp_zy_vtx_ext",    40, 0, 1038, 40, -117, 117);
TH2 * h_pfp_zy_vtx_data   = new TH2D ("h_pfp_zy_vtx_data",   "h_pfp_zy_vtx_data",   40, 0, 1038, 40, -117, 117);

TH2D * h_post_cuts_num_tracks_showers_signal_total  = new TH2D ("h_post_cuts_num_tracks_showers_signal_total", "h_post_cuts_num_tracks_showers_signal_total", 3, 1, 4, 2, 0, 2);
TH2D * h_post_cuts_num_tracks_showers_bkg_total     = new TH2D ("h_post_cuts_num_tracks_showers_bkg_total",    "h_post_cuts_num_tracks_showers_bkg_total",    3, 1, 4, 2, 0, 2);

TH1D * h_track_containment_nue_cc          = new TH1D ("h_track_containment_nue_cc",        "h_track_containment_nue_cc",        2, 0, 2);
TH1D * h_track_containment_nue_cc_out_fv   = new TH1D ("h_track_containment_nue_cc_out_fv", "h_track_containment_nue_cc_out_fv", 2, 0, 2);
TH1D * h_track_containment_nue_cc_mixed    = new TH1D ("h_track_containment_nue_cc_mixed",  "h_track_containment_nue_cc_mixed",  2, 0, 2);
TH1D * h_track_containment_numu_cc         = new TH1D ("h_track_containment_numu_cc",       "h_track_containment_numu_cc",       2, 0, 2);
TH1D * h_track_containment_numu_cc_mixed   = new TH1D ("h_track_containment_numu_cc_mixed", "h_track_containment_numu_cc_mixed", 2, 0, 2);
TH1D * h_track_containment_nc              = new TH1D ("h_track_containment_nc",            "h_track_containment_nc",            2, 0, 2);
TH1D * h_track_containment_nc_pi0          = new TH1D ("h_track_containment_nc_pi0",        "h_track_containment_nc_pi0",        2, 0, 2);
TH1D * h_track_containment_cosmic          = new TH1D ("h_track_containment_cosmic",        "h_track_containment_cosmic",        2, 0, 2);
TH1D * h_track_containment_other_mixed     = new TH1D ("h_track_containment_other_mixed",   "h_track_containment_other_mixed",   2, 0, 2);
TH1D * h_track_containment_unmatched       = new TH1D ("h_track_containment_unmatched",     "h_track_containment_unmatched",     2, 0, 2);
TH1D * h_track_containment_intime          = new TH1D ("h_track_containment_intime",        "h_track_containment_intime",        2, 0, 2);
TH1D * h_track_containment_data            = new TH1D ("h_track_containment_data",          "h_track_containment_data",          2, 0, 2);

TH1D * h_low_true_momentum = new TH1D  ("h_low_true_momentum",  "h_low_true_momentum",  20, 0, 4);
TH1D * h_med_true_momentum = new TH1D  ("h_med_true_momentum",  "h_med_true_momentum",  20, 0, 4);
TH1D * h_high_true_momentum = new TH1D ("h_high_true_momentum", "h_high_true_momentum", 20, 0, 4);

TH2D * h_dedx_collection_angle_nue_cc         = new TH2D ("h_dedx_collection_angle_nue_cc",         "h_dedx_collection_nue_cc",         20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_nue_cc_out_fv  = new TH2D ("h_dedx_collection_angle_nue_cc_out_fv",  "h_dedx_collection_nue_cc_out_fv",  20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_nue_cc_mixed   = new TH2D ("h_dedx_collection_angle_nue_cc_mixed",   "h_dedx_collection_nue_cc_mixed",   20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_numu_cc        = new TH2D ("h_dedx_collection_angle_numu_cc",        "h_dedx_collection_numu_cc",        20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_numu_cc_mixed  = new TH2D ("h_dedx_collection_angle_numu_cc_mixed",  "h_dedx_collection_numu_cc_mixed",  20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_nc             = new TH2D ("h_dedx_collection_angle_nc",             "h_dedx_collection_nc",             20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_nc_pi0         = new TH2D ("h_dedx_collection_angle_nc_pi0",         "h_dedx_collection_nc_pi0",         20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_cosmic         = new TH2D ("h_dedx_collection_angle_cosmic",         "h_dedx_collection_cosmic",         20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_other_mixed    = new TH2D ("h_dedx_collection_angle_other_mixed",    "h_dedx_collection_other_mixed",    20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_unmatched      = new TH2D ("h_dedx_collection_angle_unmatched",      "h_dedx_collection_unmatched",      20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_intime         = new TH2D ("h_dedx_collection_angle_intime",         "h_dedx_collection_intime",         20, 0, 10, 20, -3.5, 3.5);
TH2D * h_dedx_collection_angle_data           = new TH2D ("h_dedx_collection_angle_data",           "h_dedx_collection_data",           20, 0, 10, 20, -3.5, 3.5);

TH2D * h_dedx_theta_nue_cc         = new TH2D ("h_dedx_theta_nue_cc",         "h_dedx_theta_nue_cc",         10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_nue_cc_out_fv  = new TH2D ("h_dedx_theta_nue_cc_out_fv",  "h_dedx_theta_nue_cc_out_fv",  10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_nue_cc_mixed   = new TH2D ("h_dedx_theta_nue_cc_mixed",   "h_dedx_theta_nue_cc_mixed",   10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_numu_cc        = new TH2D ("h_dedx_theta_numu_cc",        "h_dedx_theta_numu_cc",        10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_numu_cc_mixed  = new TH2D ("h_dedx_theta_numu_cc_mixed",  "h_dedx_theta_numu_cc_mixed",  10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_nc             = new TH2D ("h_dedx_theta_nc",             "h_dedx_theta_nc",             10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_nc_pi0         = new TH2D ("h_dedx_theta_nc_pi0",         "h_dedx_theta_nc_pi0",         10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_cosmic         = new TH2D ("h_dedx_theta_cosmic",         "h_dedx_theta_cosmic",         10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_other_mixed    = new TH2D ("h_dedx_theta_other_mixed",    "h_dedx_theta_other_mixed",    10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_unmatched      = new TH2D ("h_dedx_theta_unmatched",      "h_dedx_theta_unmatched",      10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_intime         = new TH2D ("h_dedx_theta_intime",         "h_dedx_theta_intime",         10, 0, 5, 10, 0, 180);
TH2D * h_dedx_theta_data           = new TH2D ("h_dedx_theta_data",           "h_dedx_theta_data",           10, 0, 5, 10, 0, 180);

TH2D * h_dedx_theta_pre_cuts_nue_cc         = new TH2D ("h_dedx_theta_pre_cuts_nue_cc",         "h_dedx_theta_pre_cuts_nue_cc",         20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_nue_cc_out_fv  = new TH2D ("h_dedx_theta_pre_cuts_nue_cc_out_fv",  "h_dedx_theta_pre_cuts_nue_cc_out_fv",  20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_nue_cc_mixed   = new TH2D ("h_dedx_theta_pre_cuts_nue_cc_mixed",   "h_dedx_theta_pre_cuts_nue_cc_mixed",   20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_numu_cc        = new TH2D ("h_dedx_theta_pre_cuts_numu_cc",        "h_dedx_theta_pre_cuts_numu_cc",        20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_numu_cc_mixed  = new TH2D ("h_dedx_theta_pre_cuts_numu_cc_mixed",  "h_dedx_theta_pre_cuts_numu_cc_mixed",  20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_nc             = new TH2D ("h_dedx_theta_pre_cuts_nc",             "h_dedx_theta_pre_cuts_nc",             20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_nc_pi0         = new TH2D ("h_dedx_theta_pre_cuts_nc_pi0",         "h_dedx_theta_pre_cuts_nc_pi0",         20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_cosmic         = new TH2D ("h_dedx_theta_pre_cuts_cosmic",         "h_dedx_theta_pre_cuts_cosmic",         20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_other_mixed    = new TH2D ("h_dedx_theta_pre_cuts_other_mixed",    "h_dedx_theta_pre_cuts_other_mixed",    20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_unmatched      = new TH2D ("h_dedx_theta_pre_cuts_unmatched",      "h_dedx_theta_pre_cuts_unmatched",      20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_intime         = new TH2D ("h_dedx_theta_pre_cuts_intime",         "h_dedx_theta_pre_cuts_intime",         20, 0, 5, 20, 0, 180);
TH2D * h_dedx_theta_pre_cuts_data           = new TH2D ("h_dedx_theta_pre_cuts_data",           "h_dedx_theta_pre_cuts_data",           20, 0, 5, 20, 0, 180);

TH1D * h_mc_ele_e_0  = new TH1D ("h_mc_ele_e_0",  "h_mc_ele_e_0",  20, 0, 5);
TH1D * h_mc_ele_e_1  = new TH1D ("h_mc_ele_e_1",  "h_mc_ele_e_1",  20, 0, 5);
TH1D * h_mc_ele_e_2  = new TH1D ("h_mc_ele_e_2",  "h_mc_ele_e_2",  20, 0, 5);
TH1D * h_mc_ele_e_3  = new TH1D ("h_mc_ele_e_3",  "h_mc_ele_e_3",  20, 0, 5);
TH1D * h_mc_ele_e_4  = new TH1D ("h_mc_ele_e_4",  "h_mc_ele_e_4",  20, 0, 5);
TH1D * h_mc_ele_e_5  = new TH1D ("h_mc_ele_e_5",  "h_mc_ele_e_5",  20, 0, 5);
TH1D * h_mc_ele_e_6  = new TH1D ("h_mc_ele_e_6",  "h_mc_ele_e_6",  20, 0, 5);
TH1D * h_mc_ele_e_7  = new TH1D ("h_mc_ele_e_7",  "h_mc_ele_e_7",  20, 0, 5);
TH1D * h_mc_ele_e_8  = new TH1D ("h_mc_ele_e_8",  "h_mc_ele_e_8",  20, 0, 5);
TH1D * h_mc_ele_e_9  = new TH1D ("h_mc_ele_e_9",  "h_mc_ele_e_9",  20, 0, 5);
TH1D * h_mc_ele_e_10 = new TH1D ("h_mc_ele_e_10", "h_mc_ele_e_10", 20, 0, 5);
TH1D * h_mc_ele_e_11 = new TH1D ("h_mc_ele_e_11", "h_mc_ele_e_11", 20, 0, 5);
TH1D * h_mc_ele_e_12 = new TH1D ("h_mc_ele_e_12", "h_mc_ele_e_12", 20, 0, 5);
TH1D * h_mc_ele_e_13 = new TH1D ("h_mc_ele_e_13", "h_mc_ele_e_13", 20, 0, 5);

TH1D * h_reco_ele_e_0  = new TH1D ("h_reco_ele_e_0",  "h_reco_ele_e_0",  20, 0, 5);
TH1D * h_reco_ele_e_1  = new TH1D ("h_reco_ele_e_1",  "h_reco_ele_e_1",  20, 0, 5);
TH1D * h_reco_ele_e_2  = new TH1D ("h_reco_ele_e_2",  "h_reco_ele_e_2",  20, 0, 5);
TH1D * h_reco_ele_e_3  = new TH1D ("h_reco_ele_e_3",  "h_reco_ele_e_3",  20, 0, 5);
TH1D * h_reco_ele_e_4  = new TH1D ("h_reco_ele_e_4",  "h_reco_ele_e_4",  20, 0, 5);
TH1D * h_reco_ele_e_5  = new TH1D ("h_reco_ele_e_5",  "h_reco_ele_e_5",  20, 0, 5);
TH1D * h_reco_ele_e_6  = new TH1D ("h_reco_ele_e_6",  "h_reco_ele_e_6",  20, 0, 5);
TH1D * h_reco_ele_e_7  = new TH1D ("h_reco_ele_e_7",  "h_reco_ele_e_7",  20, 0, 5);
TH1D * h_reco_ele_e_8  = new TH1D ("h_reco_ele_e_8",  "h_reco_ele_e_8",  20, 0, 5);
TH1D * h_reco_ele_e_9  = new TH1D ("h_reco_ele_e_9",  "h_reco_ele_e_9",  20, 0, 5);
TH1D * h_reco_ele_e_10 = new TH1D ("h_reco_ele_e_10", "h_reco_ele_e_10", 20, 0, 5);
TH1D * h_reco_ele_e_11 = new TH1D ("h_reco_ele_e_11", "h_reco_ele_e_11", 20, 0, 5);
TH1D * h_reco_ele_e_12 = new TH1D ("h_reco_ele_e_12", "h_reco_ele_e_12", 20, 0, 5);
TH1D * h_reco_ele_e_13 = new TH1D ("h_reco_ele_e_13", "h_reco_ele_e_13", 20, 0, 5);

TH2D * h_mc_reco_ele_e_0  = new TH2D ("h_mc_reco_ele_e_0",  "h_mc_reco_ele_e_0",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_1  = new TH2D ("h_mc_reco_ele_e_1",  "h_mc_reco_ele_e_1",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_2  = new TH2D ("h_mc_reco_ele_e_2",  "h_mc_reco_ele_e_2",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_3  = new TH2D ("h_mc_reco_ele_e_3",  "h_mc_reco_ele_e_3",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_4  = new TH2D ("h_mc_reco_ele_e_4",  "h_mc_reco_ele_e_4",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_5  = new TH2D ("h_mc_reco_ele_e_5",  "h_mc_reco_ele_e_5",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_6  = new TH2D ("h_mc_reco_ele_e_6",  "h_mc_reco_ele_e_6",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_7  = new TH2D ("h_mc_reco_ele_e_7",  "h_mc_reco_ele_e_7",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_8  = new TH2D ("h_mc_reco_ele_e_8",  "h_mc_reco_ele_e_8",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_9  = new TH2D ("h_mc_reco_ele_e_9",  "h_mc_reco_ele_e_9",  20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_10 = new TH2D ("h_mc_reco_ele_e_10", "h_mc_reco_ele_e_10", 20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_11 = new TH2D ("h_mc_reco_ele_e_11", "h_mc_reco_ele_e_11", 20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_12 = new TH2D ("h_mc_reco_ele_e_12", "h_mc_reco_ele_e_12", 20, 0, 5, 20, 0, 5);
TH2D * h_mc_reco_ele_e_13 = new TH2D ("h_mc_reco_ele_e_13", "h_mc_reco_ele_e_13", 20, 0, 5, 20, 0, 5);

//
/*
   TH1D * h_mc_ele_e_1  = new TH1D ("h_mc_ele_e_1",  "h_mc_ele_e_1",  20, 0, 5);
   TH1D * h_mc_ele_e_2  = new TH1D ("h_mc_ele_e_2",  "h_mc_ele_e_2",  20, 0, 5);
   TH1D * h_mc_ele_e_3  = new TH1D ("h_mc_ele_e_3",  "h_mc_ele_e_3",  20, 0, 5);
   TH1D * h_mc_ele_e_4  = new TH1D ("h_mc_ele_e_4",  "h_mc_ele_e_4",  20, 0, 5);
   TH1D * h_mc_ele_e_5  = new TH1D ("h_mc_ele_e_5",  "h_mc_ele_e_5",  20, 0, 5);
   TH1D * h_mc_ele_e_6  = new TH1D ("h_mc_ele_e_6",  "h_mc_ele_e_6",  20, 0, 5);
   TH1D * h_mc_ele_e_7  = new TH1D ("h_mc_ele_e_7",  "h_mc_ele_e_7",  20, 0, 5);
   TH1D * h_mc_ele_e_8  = new TH1D ("h_mc_ele_e_8",  "h_mc_ele_e_8",  20, 0, 5);
   TH1D * h_mc_ele_e_9  = new TH1D ("h_mc_ele_e_9",  "h_mc_ele_e_9",  20, 0, 5);
   TH1D * h_mc_ele_e_10 = new TH1D ("h_mc_ele_e_10", "h_mc_ele_e_10", 20, 0, 5);
   TH1D * h_mc_ele_e_11 = new TH1D ("h_mc_ele_e_11", "h_mc_ele_e_11", 20, 0, 5);
   TH1D * h_mc_ele_e_12 = new TH1D ("h_mc_ele_e_12", "h_mc_ele_e_12", 20, 0, 5);
   TH1D * h_mc_ele_e_13 = new TH1D ("h_mc_ele_e_13", "h_mc_ele_e_13", 20, 0, 5);
 */

TH1D * h_ele_momentum_no_cut_nue_cc         = new TH1D ("h_ele_momentum_no_cut_nue_cc",          "h_ele_momentum_no_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_no_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_no_cut_nue_cc_out_fv",   "h_ele_momentum_no_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_no_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_no_cut_nue_cc_mixed",    "h_ele_mometnum_no_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_no_cut_numu_cc        = new TH1D ("h_ele_momentum_no_cut_numu_cc",         "h_ele_momentum_no_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_no_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_no_cut_numu_cc_mixed",   "h_ele_momentum_no_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_no_cut_nc             = new TH1D ("h_ele_momentum_no_cut_nc",              "h_ele_momentum_no_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_no_cut_nc_pi0         = new TH1D ("h_ele_momentum_no_cut_nc_pi0",          "h_ele_momentum_no_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_no_cut_cosmic         = new TH1D ("h_ele_momentum_no_cut_cosmic",          "h_ele_momentum_no_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_no_cut_other_mixed    = new TH1D ("h_ele_momentum_no_cut_other_mixed",     "h_ele_momentum_no_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_no_cut_unmatched      = new TH1D ("h_ele_momentum_no_cut_unmatched",       "h_ele_momentum_no_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_no_cut_intime         = new TH1D ("h_ele_mometnum_no_cut_intime",          "h_ele_momentum_no_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_no_cut_data           = new TH1D ("h_ele_momentum_no_cut_data",            "h_ele_momentum_no_cut_data",          20, 0, 4);


TH1D * h_ele_momentum_nue_cut_nue_cc         = new TH1D ("h_ele_momentum_nue_cut_nue_cc",          "h_ele_momentum_nue_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_nue_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_nue_cut_nue_cc_out_fv",   "h_ele_momentum_nue_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_nue_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_nue_cut_nue_cc_mixed",    "h_ele_mometnum_nue_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_nue_cut_numu_cc        = new TH1D ("h_ele_momentum_nue_cut_numu_cc",         "h_ele_momentum_nue_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_nue_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_nue_cut_numu_cc_mixed",   "h_ele_momentum_nue_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_nue_cut_nc             = new TH1D ("h_ele_momentum_nue_cut_nc",              "h_ele_momentum_nue_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_nue_cut_nc_pi0         = new TH1D ("h_ele_momentum_nue_cut_nc_pi0",          "h_ele_momentum_nue_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_nue_cut_cosmic         = new TH1D ("h_ele_momentum_nue_cut_cosmic",          "h_ele_momentum_nue_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_nue_cut_other_mixed    = new TH1D ("h_ele_momentum_nue_cut_other_mixed",     "h_ele_momentum_nue_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_nue_cut_unmatched      = new TH1D ("h_ele_momentum_nue_cut_unmatched",       "h_ele_momentum_nue_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_nue_cut_intime         = new TH1D ("h_ele_mometnum_nue_cut_intime",          "h_ele_momentum_nue_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_nue_cut_data           = new TH1D ("h_ele_momentum_nue_cut_data",            "h_ele_momentum_nue_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_fv_cut_nue_cc         = new TH1D ("h_ele_momentum_fv_cut_nue_cc",          "h_ele_momentum_fv_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_fv_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_fv_cut_nue_cc_out_fv",   "h_ele_momentum_fv_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_fv_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_fv_cut_nue_cc_mixed",    "h_ele_mometnum_fv_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_fv_cut_numu_cc        = new TH1D ("h_ele_momentum_fv_cut_numu_cc",         "h_ele_momentum_fv_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_fv_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_fv_cut_numu_cc_mixed",   "h_ele_momentum_fv_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_fv_cut_nc             = new TH1D ("h_ele_momentum_fv_cut_nc",              "h_ele_momentum_fv_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_fv_cut_nc_pi0         = new TH1D ("h_ele_momentum_fv_cut_nc_pi0",          "h_ele_momentum_fv_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_fv_cut_cosmic         = new TH1D ("h_ele_momentum_fv_cut_cosmic",          "h_ele_momentum_fv_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_fv_cut_other_mixed    = new TH1D ("h_ele_momentum_fv_cut_other_mixed",     "h_ele_momentum_fv_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_fv_cut_unmatched      = new TH1D ("h_ele_momentum_fv_cut_unmatched",       "h_ele_momentum_fv_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_fv_cut_intime         = new TH1D ("h_ele_mometnum_fv_cut_intime",          "h_ele_momentum_fv_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_fv_cut_data           = new TH1D ("h_ele_momentum_fv_cut_data",            "h_ele_momentum_fv_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_flash_vtx_cut_nue_cc         = new TH1D ("h_ele_momentum_flash_vtx_cut_nue_cc",          "h_ele_momentum_flash_vtx_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_flash_vtx_cut_nue_cc_out_fv",   "h_ele_momentum_flash_vtx_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_flash_vtx_cut_nue_cc_mixed",    "h_ele_mometnum_flash_vtx_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_numu_cc        = new TH1D ("h_ele_momentum_flash_vtx_cut_numu_cc",         "h_ele_momentum_flash_vtx_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_flash_vtx_cut_numu_cc_mixed",   "h_ele_momentum_flash_vtx_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_nc             = new TH1D ("h_ele_momentum_flash_vtx_cut_nc",              "h_ele_momentum_flash_vtx_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_nc_pi0         = new TH1D ("h_ele_momentum_flash_vtx_cut_nc_pi0",          "h_ele_momentum_flash_vtx_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_cosmic         = new TH1D ("h_ele_momentum_flash_vtx_cut_cosmic",          "h_ele_momentum_flash_vtx_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_other_mixed    = new TH1D ("h_ele_momentum_flash_vtx_cut_other_mixed",     "h_ele_momentum_flash_vtx_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_unmatched      = new TH1D ("h_ele_momentum_flash_vtx_cut_unmatched",       "h_ele_momentum_flash_vtx_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_intime         = new TH1D ("h_ele_mometnum_flash_vtx_cut_intime",          "h_ele_momentum_flash_vtx_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_flash_vtx_cut_data           = new TH1D ("h_ele_momentum_flash_vtx_cut_data",            "h_ele_momentum_flash_vtx_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_shwr_vtx_cut_nue_cc         = new TH1D ("h_ele_momentum_shwr_vtx_cut_nue_cc",          "h_ele_momentum_shwr_vtx_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_shwr_vtx_cut_nue_cc_out_fv",   "h_ele_momentum_shwr_vtx_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_shwr_vtx_cut_nue_cc_mixed",    "h_ele_mometnum_shwr_vtx_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_numu_cc        = new TH1D ("h_ele_momentum_shwr_vtx_cut_numu_cc",         "h_ele_momentum_shwr_vtx_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_shwr_vtx_cut_numu_cc_mixed",   "h_ele_momentum_shwr_vtx_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_nc             = new TH1D ("h_ele_momentum_shwr_vtx_cut_nc",              "h_ele_momentum_shwr_vtx_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_nc_pi0         = new TH1D ("h_ele_momentum_shwr_vtx_cut_nc_pi0",          "h_ele_momentum_shwr_vtx_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_cosmic         = new TH1D ("h_ele_momentum_shwr_vtx_cut_cosmic",          "h_ele_momentum_shwr_vtx_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_other_mixed    = new TH1D ("h_ele_momentum_shwr_vtx_cut_other_mixed",     "h_ele_momentum_shwr_vtx_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_unmatched      = new TH1D ("h_ele_momentum_shwr_vtx_cut_unmatched",       "h_ele_momentum_shwr_vtx_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_intime         = new TH1D ("h_ele_mometnum_shwr_vtx_cut_intime",          "h_ele_momentum_shwr_vtx_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_shwr_vtx_cut_data           = new TH1D ("h_ele_momentum_shwr_vtx_cut_data",            "h_ele_momentum_shwr_vtx_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_trk_vtx_cut_nue_cc         = new TH1D ("h_ele_momentum_trk_vtx_cut_nue_cc",          "h_ele_momentum_trk_vtx_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_trk_vtx_cut_nue_cc_out_fv",   "h_ele_momentum_trk_vtx_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_trk_vtx_cut_nue_cc_mixed",    "h_ele_mometnum_trk_vtx_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_numu_cc        = new TH1D ("h_ele_momentum_trk_vtx_cut_numu_cc",         "h_ele_momentum_trk_vtx_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_trk_vtx_cut_numu_cc_mixed",   "h_ele_momentum_trk_vtx_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_nc             = new TH1D ("h_ele_momentum_trk_vtx_cut_nc",              "h_ele_momentum_trk_vtx_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_nc_pi0         = new TH1D ("h_ele_momentum_trk_vtx_cut_nc_pi0",          "h_ele_momentum_trk_vtx_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_cosmic         = new TH1D ("h_ele_momentum_trk_vtx_cut_cosmic",          "h_ele_momentum_trk_vtx_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_other_mixed    = new TH1D ("h_ele_momentum_trk_vtx_cut_other_mixed",     "h_ele_momentum_trk_vtx_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_unmatched      = new TH1D ("h_ele_momentum_trk_vtx_cut_unmatched",       "h_ele_momentum_trk_vtx_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_intime         = new TH1D ("h_ele_mometnum_trk_vtx_cut_intime",          "h_ele_momentum_trk_vtx_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_trk_vtx_cut_data           = new TH1D ("h_ele_momentum_trk_vtx_cut_data",            "h_ele_momentum_trk_vtx_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_hit_cut_nue_cc         = new TH1D ("h_ele_momentum_hit_cut_nue_cc",          "h_ele_momentum_hit_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_hit_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_hit_cut_nue_cc_out_fv",   "h_ele_momentum_hit_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_hit_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_hit_cut_nue_cc_mixed",    "h_ele_mometnum_hit_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_hit_cut_numu_cc        = new TH1D ("h_ele_momentum_hit_cut_numu_cc",         "h_ele_momentum_hit_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_hit_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_hit_cut_numu_cc_mixed",   "h_ele_momentum_hit_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_hit_cut_nc             = new TH1D ("h_ele_momentum_hit_cut_nc",              "h_ele_momentum_hit_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_hit_cut_nc_pi0         = new TH1D ("h_ele_momentum_hit_cut_nc_pi0",          "h_ele_momentum_hit_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_hit_cut_cosmic         = new TH1D ("h_ele_momentum_hit_cut_cosmic",          "h_ele_momentum_hit_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_hit_cut_other_mixed    = new TH1D ("h_ele_momentum_hit_cut_other_mixed",     "h_ele_momentum_hit_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_hit_cut_unmatched      = new TH1D ("h_ele_momentum_hit_cut_unmatched",       "h_ele_momentum_hit_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_hit_cut_intime         = new TH1D ("h_ele_mometnum_hit_cut_intime",          "h_ele_momentum_hit_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_hit_cut_data           = new TH1D ("h_ele_momentum_hit_cut_data",            "h_ele_momentum_hit_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_yhit_cut_nue_cc         = new TH1D ("h_ele_momentum_yhit_cut_nue_cc",          "h_ele_momentum_yhit_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_yhit_cut_nue_cc_out_fv",   "h_ele_momentum_yhit_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_yhit_cut_nue_cc_mixed",    "h_ele_mometnum_yhit_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_numu_cc        = new TH1D ("h_ele_momentum_yhit_cut_numu_cc",         "h_ele_momentum_yhit_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_yhit_cut_numu_cc_mixed",   "h_ele_momentum_yhit_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_nc             = new TH1D ("h_ele_momentum_yhit_cut_nc",              "h_ele_momentum_yhit_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_nc_pi0         = new TH1D ("h_ele_momentum_yhit_cut_nc_pi0",          "h_ele_momentum_yhit_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_cosmic         = new TH1D ("h_ele_momentum_yhit_cut_cosmic",          "h_ele_momentum_yhit_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_other_mixed    = new TH1D ("h_ele_momentum_yhit_cut_other_mixed",     "h_ele_momentum_yhit_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_unmatched      = new TH1D ("h_ele_momentum_yhit_cut_unmatched",       "h_ele_momentum_yhit_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_intime         = new TH1D ("h_ele_mometnum_yhit_cut_intime",          "h_ele_momentum_yhit_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_yhit_cut_data           = new TH1D ("h_ele_momentum_yhit_cut_data",            "h_ele_momentum_yhit_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_open_angle_cut_nue_cc         = new TH1D ("h_ele_momentum_open_angle_cut_nue_cc",          "h_ele_momentum_open_angle_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_open_angle_cut_nue_cc_out_fv",   "h_ele_momentum_open_angle_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_open_angle_cut_nue_cc_mixed",    "h_ele_mometnum_open_angle_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_numu_cc        = new TH1D ("h_ele_momentum_open_angle_cut_numu_cc",         "h_ele_momentum_open_angle_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_open_angle_cut_numu_cc_mixed",   "h_ele_momentum_open_angle_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_nc             = new TH1D ("h_ele_momentum_open_angle_cut_nc",              "h_ele_momentum_open_angle_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_nc_pi0         = new TH1D ("h_ele_momentum_open_angle_cut_nc_pi0",          "h_ele_momentum_open_angle_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_cosmic         = new TH1D ("h_ele_momentum_open_angle_cut_cosmic",          "h_ele_momentum_open_angle_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_other_mixed    = new TH1D ("h_ele_momentum_open_angle_cut_other_mixed",     "h_ele_momentum_open_angle_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_unmatched      = new TH1D ("h_ele_momentum_open_angle_cut_unmatched",       "h_ele_momentum_open_angle_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_intime         = new TH1D ("h_ele_mometnum_open_angle_cut_intime",          "h_ele_momentum_open_angle_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_open_angle_cut_data           = new TH1D ("h_ele_momentum_open_angle_cut_data",            "h_ele_momentum_open_angle_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_dedx_cut_nue_cc         = new TH1D ("h_ele_momentum_dedx_cut_nue_cc",          "h_ele_momentum_dedx_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_dedx_cut_nue_cc_out_fv",   "h_ele_momentum_dedx_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_dedx_cut_nue_cc_mixed",    "h_ele_mometnum_dedx_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_numu_cc        = new TH1D ("h_ele_momentum_dedx_cut_numu_cc",         "h_ele_momentum_dedx_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_dedx_cut_numu_cc_mixed",   "h_ele_momentum_dedx_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_nc             = new TH1D ("h_ele_momentum_dedx_cut_nc",              "h_ele_momentum_dedx_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_nc_pi0         = new TH1D ("h_ele_momentum_dedx_cut_nc_pi0",          "h_ele_momentum_dedx_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_cosmic         = new TH1D ("h_ele_momentum_dedx_cut_cosmic",          "h_ele_momentum_dedx_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_other_mixed    = new TH1D ("h_ele_momentum_dedx_cut_other_mixed",     "h_ele_momentum_dedx_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_unmatched      = new TH1D ("h_ele_momentum_dedx_cut_unmatched",       "h_ele_momentum_dedx_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_intime         = new TH1D ("h_ele_mometnum_dedx_cut_intime",          "h_ele_momentum_dedx_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_dedx_cut_data           = new TH1D ("h_ele_momentum_dedx_cut_data",            "h_ele_momentum_dedx_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_2shwr_cut_nue_cc         = new TH1D ("h_ele_momentum_2shwr_cut_nue_cc",          "h_ele_momentum_2shwr_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_2shwr_cut_nue_cc_out_fv",   "h_ele_momentum_2shwr_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_2shwr_cut_nue_cc_mixed",    "h_ele_mometnum_2shwr_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_numu_cc        = new TH1D ("h_ele_momentum_2shwr_cut_numu_cc",         "h_ele_momentum_2shwr_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_2shwr_cut_numu_cc_mixed",   "h_ele_momentum_2shwr_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_nc             = new TH1D ("h_ele_momentum_2shwr_cut_nc",              "h_ele_momentum_2shwr_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_nc_pi0         = new TH1D ("h_ele_momentum_2shwr_cut_nc_pi0",          "h_ele_momentum_2shwr_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_cosmic         = new TH1D ("h_ele_momentum_2shwr_cut_cosmic",          "h_ele_momentum_2shwr_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_other_mixed    = new TH1D ("h_ele_momentum_2shwr_cut_other_mixed",     "h_ele_momentum_2shwr_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_unmatched      = new TH1D ("h_ele_momentum_2shwr_cut_unmatched",       "h_ele_momentum_2shwr_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_intime         = new TH1D ("h_ele_mometnum_2shwr_cut_intime",          "h_ele_momentum_2shwr_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_2shwr_cut_data           = new TH1D ("h_ele_momentum_2shwr_cut_data",            "h_ele_momentum_2shwr_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_hit_length_cut_nue_cc         = new TH1D ("h_ele_momentum_hit_length_cut_nue_cc",          "h_ele_momentum_hit_length_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_hit_length_cut_nue_cc_out_fv",   "h_ele_momentum_hit_length_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_hit_length_cut_nue_cc_mixed",    "h_ele_mometnum_hit_length_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_numu_cc        = new TH1D ("h_ele_momentum_hit_length_cut_numu_cc",         "h_ele_momentum_hit_length_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_hit_length_cut_numu_cc_mixed",   "h_ele_momentum_hit_length_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_nc             = new TH1D ("h_ele_momentum_hit_length_cut_nc",              "h_ele_momentum_hit_length_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_nc_pi0         = new TH1D ("h_ele_momentum_hit_length_cut_nc_pi0",          "h_ele_momentum_hit_length_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_cosmic         = new TH1D ("h_ele_momentum_hit_length_cut_cosmic",          "h_ele_momentum_hit_length_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_other_mixed    = new TH1D ("h_ele_momentum_hit_length_cut_other_mixed",     "h_ele_momentum_hit_length_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_unmatched      = new TH1D ("h_ele_momentum_hit_length_cut_unmatched",       "h_ele_momentum_hit_length_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_intime         = new TH1D ("h_ele_mometnum_hit_length_cut_intime",          "h_ele_momentum_hit_length_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_hit_length_cut_data           = new TH1D ("h_ele_momentum_hit_length_cut_data",            "h_ele_momentum_hit_length_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_length_ratio_cut_nue_cc         = new TH1D ("h_ele_momentum_length_ratio_cut_nue_cc",          "h_ele_momentum_length_ratio_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_length_ratio_cut_nue_cc_out_fv",   "h_ele_momentum_length_ratio_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_length_ratio_cut_nue_cc_mixed",    "h_ele_mometnum_length_ratio_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_numu_cc        = new TH1D ("h_ele_momentum_length_ratio_cut_numu_cc",         "h_ele_momentum_length_ratio_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_length_ratio_cut_numu_cc_mixed",   "h_ele_momentum_length_ratio_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_nc             = new TH1D ("h_ele_momentum_length_ratio_cut_nc",              "h_ele_momentum_length_ratio_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_nc_pi0         = new TH1D ("h_ele_momentum_length_ratio_cut_nc_pi0",          "h_ele_momentum_length_ratio_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_cosmic         = new TH1D ("h_ele_momentum_length_ratio_cut_cosmic",          "h_ele_momentum_length_ratio_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_other_mixed    = new TH1D ("h_ele_momentum_length_ratio_cut_other_mixed",     "h_ele_momentum_length_ratio_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_unmatched      = new TH1D ("h_ele_momentum_length_ratio_cut_unmatched",       "h_ele_momentum_length_ratio_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_intime         = new TH1D ("h_ele_mometnum_length_ratio_cut_intime",          "h_ele_momentum_length_ratio_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_length_ratio_cut_data           = new TH1D ("h_ele_momentum_length_ratio_cut_data",            "h_ele_momentum_length_ratio_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_containment_cut_nue_cc         = new TH1D ("h_ele_momentum_containment_cut_nue_cc",          "h_ele_momentum_containment_cut_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_containment_cut_nue_cc_out_fv  = new TH1D ("h_ele_momentum_containment_cut_nue_cc_out_fv",   "h_ele_momentum_containment_cut_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_containment_cut_nue_cc_mixed   = new TH1D ("h_ele_momentum_containment_cut_nue_cc_mixed",    "h_ele_mometnum_containment_cut_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_containment_cut_numu_cc        = new TH1D ("h_ele_momentum_containment_cut_numu_cc",         "h_ele_momentum_containment_cut_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_containment_cut_numu_cc_mixed  = new TH1D ("h_ele_momentum_containment_cut_numu_cc_mixed",   "h_ele_momentum_containment_cut_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_containment_cut_nc             = new TH1D ("h_ele_momentum_containment_cut_nc",              "h_ele_momentum_containment_cut_nc",            20, 0, 4);
TH1D * h_ele_momentum_containment_cut_nc_pi0         = new TH1D ("h_ele_momentum_containment_cut_nc_pi0",          "h_ele_momentum_containment_cut_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_containment_cut_cosmic         = new TH1D ("h_ele_momentum_containment_cut_cosmic",          "h_ele_momentum_containment_cut_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_containment_cut_other_mixed    = new TH1D ("h_ele_momentum_containment_cut_other_mixed",     "h_ele_momentum_containment_cut_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_containment_cut_unmatched      = new TH1D ("h_ele_momentum_containment_cut_unmatched",       "h_ele_momentum_containment_cut_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_containment_cut_intime         = new TH1D ("h_ele_mometnum_containment_cut_intime",          "h_ele_momentum_containment_cut_intime",        20, 0, 4);
TH1D * h_ele_momentum_containment_cut_data           = new TH1D ("h_ele_momentum_containment_cut_data",            "h_ele_momentum_containment_cut_data",          20, 0, 4);

TH1D * h_ele_momentum_slice_1_nue_cc         = new TH1D ("h_ele_momentum_slice_1_nue_cc",          "h_ele_momentum_slice_1_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_slice_1_nue_cc_out_fv  = new TH1D ("h_ele_momentum_slice_1_nue_cc_out_fv",   "h_ele_momentum_slice_1_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_slice_1_nue_cc_mixed   = new TH1D ("h_ele_momentum_slice_1_nue_cc_mixed",    "h_ele_mometnum_slice_1_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_slice_1_numu_cc        = new TH1D ("h_ele_momentum_slice_1_numu_cc",         "h_ele_momentum_slice_1_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_slice_1_numu_cc_mixed  = new TH1D ("h_ele_momentum_slice_1_numu_cc_mixed",   "h_ele_momentum_slice_1_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_slice_1_nc             = new TH1D ("h_ele_momentum_slice_1_nc",              "h_ele_momentum_slice_1_nc",            20, 0, 4);
TH1D * h_ele_momentum_slice_1_nc_pi0         = new TH1D ("h_ele_momentum_slice_1_nc_pi0",          "h_ele_momentum_slice_1_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_slice_1_cosmic         = new TH1D ("h_ele_momentum_slice_1_cosmic",          "h_ele_momentum_slice_1_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_slice_1_other_mixed    = new TH1D ("h_ele_momentum_slice_1_other_mixed",     "h_ele_momentum_slice_1_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_slice_1_unmatched      = new TH1D ("h_ele_momentum_slice_1_unmatched",       "h_ele_momentum_slice_1_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_slice_1_intime         = new TH1D ("h_ele_mometnum_slice_1_intime",          "h_ele_momentum_slice_1_intime",        20, 0, 4);
TH1D * h_ele_momentum_slice_1_data           = new TH1D ("h_ele_momentum_slice_1_data",            "h_ele_momentum_slice_1_data",          20, 0, 4);

TH1D * h_ele_momentum_slice_2_nue_cc         = new TH1D ("h_ele_momentum_slice_2_nue_cc",          "h_ele_momentum_slice_2_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_slice_2_nue_cc_out_fv  = new TH1D ("h_ele_momentum_slice_2_nue_cc_out_fv",   "h_ele_momentum_slice_2_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_slice_2_nue_cc_mixed   = new TH1D ("h_ele_momentum_slice_2_nue_cc_mixed",    "h_ele_mometnum_slice_2_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_slice_2_numu_cc        = new TH1D ("h_ele_momentum_slice_2_numu_cc",         "h_ele_momentum_slice_2_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_slice_2_numu_cc_mixed  = new TH1D ("h_ele_momentum_slice_2_numu_cc_mixed",   "h_ele_momentum_slice_2_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_slice_2_nc             = new TH1D ("h_ele_momentum_slice_2_nc",              "h_ele_momentum_slice_2_nc",            20, 0, 4);
TH1D * h_ele_momentum_slice_2_nc_pi0         = new TH1D ("h_ele_momentum_slice_2_nc_pi0",          "h_ele_momentum_slice_2_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_slice_2_cosmic         = new TH1D ("h_ele_momentum_slice_2_cosmic",          "h_ele_momentum_slice_2_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_slice_2_other_mixed    = new TH1D ("h_ele_momentum_slice_2_other_mixed",     "h_ele_momentum_slice_2_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_slice_2_unmatched      = new TH1D ("h_ele_momentum_slice_2_unmatched",       "h_ele_momentum_slice_2_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_slice_2_intime         = new TH1D ("h_ele_mometnum_slice_2_intime",          "h_ele_momentum_slice_2_intime",        20, 0, 4);
TH1D * h_ele_momentum_slice_2_data           = new TH1D ("h_ele_momentum_slice_2_data",            "h_ele_momentum_slice_2_data",          20, 0, 4);

TH1D * h_ele_momentum_slice_3_nue_cc         = new TH1D ("h_ele_momentum_slice_3_nue_cc",          "h_ele_momentum_slice_3_nue_cc",        20, 0, 4);
TH1D * h_ele_momentum_slice_3_nue_cc_out_fv  = new TH1D ("h_ele_momentum_slice_3_nue_cc_out_fv",   "h_ele_momentum_slice_3_nue_cc_out_fv", 20, 0, 4);
TH1D * h_ele_momentum_slice_3_nue_cc_mixed   = new TH1D ("h_ele_momentum_slice_3_nue_cc_mixed",    "h_ele_mometnum_slice_3_nue_cc_mixed",  20, 0, 4);
TH1D * h_ele_momentum_slice_3_numu_cc        = new TH1D ("h_ele_momentum_slice_3_numu_cc",         "h_ele_momentum_slice_3_numu_cc",       20, 0, 4);
TH1D * h_ele_momentum_slice_3_numu_cc_mixed  = new TH1D ("h_ele_momentum_slice_3_numu_cc_mixed",   "h_ele_momentum_slice_3_numu_cc_mixed", 20, 0, 4);
TH1D * h_ele_momentum_slice_3_nc             = new TH1D ("h_ele_momentum_slice_3_nc",              "h_ele_momentum_slice_3_nc",            20, 0, 4);
TH1D * h_ele_momentum_slice_3_nc_pi0         = new TH1D ("h_ele_momentum_slice_3_nc_pi0",          "h_ele_momentum_slice_3_nc_pi0",        20, 0, 4);
TH1D * h_ele_momentum_slice_3_cosmic         = new TH1D ("h_ele_momentum_slice_3_cosmic",          "h_ele_momentum_slice_3_cosmic",        20, 0, 4);
TH1D * h_ele_momentum_slice_3_other_mixed    = new TH1D ("h_ele_momentum_slice_3_other_mixed",     "h_ele_momentum_slice_3_other_mixed",   20, 0, 4);
TH1D * h_ele_momentum_slice_3_unmatched      = new TH1D ("h_ele_momentum_slice_3_unmatched",       "h_ele_momentum_slice_3_unmatched",     20, 0, 4);
TH1D * h_ele_momentum_slice_3_intime         = new TH1D ("h_ele_mometnum_slice_3_intime",          "h_ele_momentum_slice_3_intime",        20, 0, 4);
TH1D * h_ele_momentum_slice_3_data           = new TH1D ("h_ele_momentum_slice_3_data",            "h_ele_momentum_slice_3_data",          20, 0, 4);

TH1D * h_dedx_slice_1_nue_cc         = new TH1D ("h_dedx_slice_1_nue_cc",          "h_dedx_slice_1_nue_cc",        20, 0, 4);
TH1D * h_dedx_slice_1_nue_cc_out_fv  = new TH1D ("h_dedx_slice_1_nue_cc_out_fv",   "h_dedx_slice_1_nue_cc_out_fv", 20, 0, 4);
TH1D * h_dedx_slice_1_nue_cc_mixed   = new TH1D ("h_dedx_slice_1_nue_cc_mixed",    "h_dedx_slice_1_nue_cc_mixed",  20, 0, 4);
TH1D * h_dedx_slice_1_numu_cc        = new TH1D ("h_dedx_slice_1_numu_cc",         "h_dedx_slice_1_numu_cc",       20, 0, 4);
TH1D * h_dedx_slice_1_numu_cc_mixed  = new TH1D ("h_dedx_slice_1_numu_cc_mixed",   "h_dedx_slice_1_numu_cc_mixed", 20, 0, 4);
TH1D * h_dedx_slice_1_nc             = new TH1D ("h_dedx_slice_1_nc",              "h_dedx_slice_1_nc",            20, 0, 4);
TH1D * h_dedx_slice_1_nc_pi0         = new TH1D ("h_dedx_slice_1_nc_pi0",          "h_dedx_slice_1_nc_pi0",        20, 0, 4);
TH1D * h_dedx_slice_1_cosmic         = new TH1D ("h_dedx_slice_1_cosmic",          "h_dedx_slice_1_cosmic",        20, 0, 4);
TH1D * h_dedx_slice_1_other_mixed    = new TH1D ("h_dedx_slice_1_other_mixed",     "h_dedx_slice_1_other_mixed",   20, 0, 4);
TH1D * h_dedx_slice_1_unmatched      = new TH1D ("h_dedx_slice_1_unmatched",       "h_dedx_slice_1_unmatched",     20, 0, 4);
TH1D * h_dedx_slice_1_intime         = new TH1D ("h_dedx_slice_1_intime",          "h_dedx_slice_1_intime",        20, 0, 4);
TH1D * h_dedx_slice_1_data           = new TH1D ("h_dedx_slice_1_data",            "h_dedx_slice_1_data",          20, 0, 4);

TH1D * h_dedx_slice_2_nue_cc         = new TH1D ("h_dedx_slice_2_nue_cc",          "h_dedx_slice_2_nue_cc",        20, 0, 4);
TH1D * h_dedx_slice_2_nue_cc_out_fv  = new TH1D ("h_dedx_slice_2_nue_cc_out_fv",   "h_dedx_slice_2_nue_cc_out_fv", 20, 0, 4);
TH1D * h_dedx_slice_2_nue_cc_mixed   = new TH1D ("h_dedx_slice_2_nue_cc_mixed",    "h_dedx_slice_2_nue_cc_mixed",  20, 0, 4);
TH1D * h_dedx_slice_2_numu_cc        = new TH1D ("h_dedx_slice_2_numu_cc",         "h_dedx_slice_2_numu_cc",       20, 0, 4);
TH1D * h_dedx_slice_2_numu_cc_mixed  = new TH1D ("h_dedx_slice_2_numu_cc_mixed",   "h_dedx_slice_2_numu_cc_mixed", 20, 0, 4);
TH1D * h_dedx_slice_2_nc             = new TH1D ("h_dedx_slice_2_nc",              "h_dedx_slice_2_nc",            20, 0, 4);
TH1D * h_dedx_slice_2_nc_pi0         = new TH1D ("h_dedx_slice_2_nc_pi0",          "h_dedx_slice_2_nc_pi0",        20, 0, 4);
TH1D * h_dedx_slice_2_cosmic         = new TH1D ("h_dedx_slice_2_cosmic",          "h_dedx_slice_2_cosmic",        20, 0, 4);
TH1D * h_dedx_slice_2_other_mixed    = new TH1D ("h_dedx_slice_2_other_mixed",     "h_dedx_slice_2_other_mixed",   20, 0, 4);
TH1D * h_dedx_slice_2_unmatched      = new TH1D ("h_dedx_slice_2_unmatched",       "h_dedx_slice_2_unmatched",     20, 0, 4);
TH1D * h_dedx_slice_2_intime         = new TH1D ("h_dedx_slice_2_intime",          "h_dedx_slice_2_intime",        20, 0, 4);
TH1D * h_dedx_slice_2_data           = new TH1D ("h_dedx_slice_2_data",            "h_dedx_slice_2_data",          20, 0, 4);

TH1D * h_dedx_slice_3_nue_cc         = new TH1D ("h_dedx_slice_3_nue_cc",          "h_dedx_slice_3_nue_cc",        20, 0, 4);
TH1D * h_dedx_slice_3_nue_cc_out_fv  = new TH1D ("h_dedx_slice_3_nue_cc_out_fv",   "h_dedx_slice_3_nue_cc_out_fv", 20, 0, 4);
TH1D * h_dedx_slice_3_nue_cc_mixed   = new TH1D ("h_dedx_slice_3_nue_cc_mixed",    "h_dedx_slice_3_nue_cc_mixed",  20, 0, 4);
TH1D * h_dedx_slice_3_numu_cc        = new TH1D ("h_dedx_slice_3_numu_cc",         "h_dedx_slice_3_numu_cc",       20, 0, 4);
TH1D * h_dedx_slice_3_numu_cc_mixed  = new TH1D ("h_dedx_slice_3_numu_cc_mixed",   "h_dedx_slice_3_numu_cc_mixed", 20, 0, 4);
TH1D * h_dedx_slice_3_nc             = new TH1D ("h_dedx_slice_3_nc",              "h_dedx_slice_3_nc",            20, 0, 4);
TH1D * h_dedx_slice_3_nc_pi0         = new TH1D ("h_dedx_slice_3_nc_pi0",          "h_dedx_slice_3_nc_pi0",        20, 0, 4);
TH1D * h_dedx_slice_3_cosmic         = new TH1D ("h_dedx_slice_3_cosmic",          "h_dedx_slice_3_cosmic",        20, 0, 4);
TH1D * h_dedx_slice_3_other_mixed    = new TH1D ("h_dedx_slice_3_other_mixed",     "h_dedx_slice_3_other_mixed",   20, 0, 4);
TH1D * h_dedx_slice_3_unmatched      = new TH1D ("h_dedx_slice_3_unmatched",       "h_dedx_slice_3_unmatched",     20, 0, 4);
TH1D * h_dedx_slice_3_intime         = new TH1D ("h_dedx_slice_3_intime",          "h_dedx_slice_3_intime",        20, 0, 4);
TH1D * h_dedx_slice_3_data           = new TH1D ("h_dedx_slice_3_data",            "h_dedx_slice_3_data",          20, 0, 4);

}; //end class

}//end namespace

#endif
