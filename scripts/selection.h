#ifndef SELECTION_h
#define SELECTION_h

#include "selection_functions.h"
#include "selection_cuts.h"
#include "histogram_functions.h"

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

//hit threshold for at least one shower
//standard 100 hits
const double shwr_hit_threshold = 120;//hits

//hit threshold for at least one shower on collection plane
//standard 50 hits
const double shwr_hit_threshold_collection = 50;//hits

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
TH1D * h_nue_eng_eff_den           = new TH1D("h_nue_eng_eff_den", "h_nue_eng_eff_den", 8, 0, 4);
TH1D * h_nue_eng_eff_num           = new TH1D("h_nue_eng_eff_num", "h_nue_eng_eff_num", 8, 0, 4);
TH1D * h_ele_eng_eff_den           = new TH1D("h_ele_eng_eff_den", "h_nue_ele_eff_den", 8, 0, 4);
TH1D * h_ele_eng_eff_num           = new TH1D("h_ele_eng_eff_num", "h_nue_ele_eff_num", 8, 0, 4);
TH1D * h_nue_vtx_x_eff_den         = new TH1D("h_nue_vtx_x_eff_den", "h_nue_vtx_x_eff_den", 10, 0, 256.35);
TH1D * h_nue_vtx_x_eff_num         = new TH1D("h_nue_vtx_x_eff_num", "h_nue_vtx_x_eff_num", 10, 0, 256.35);
TH1D * h_nue_vtx_y_eff_den         = new TH1D("h_nue_vtx_y_eff_den", "h_nue_vtx_y_eff_den", 10, -116.5, 116.5);
TH1D * h_nue_vtx_y_eff_num         = new TH1D("h_nue_vtx_y_eff_num", "h_nue_vtx_y_eff_num", 10, -116.5, 116.5);
TH1D * h_nue_vtx_z_eff_den         = new TH1D("h_nue_vtx_z_eff_den", "h_nue_vtx_z_eff_den", 10, 0, 1036.8);
TH1D * h_nue_vtx_z_eff_num         = new TH1D("h_nue_vtx_z_eff_num", "h_nue_vtx_z_eff_num", 10, 0, 1036.8);
TH1D * h_nue_dir_x_eff_den         = new TH1D("h_nue_dir_x_eff_den", "h_nue_dir_x_eff_den", 20, -1, 1);
TH1D * h_nue_dir_x_eff_num         = new TH1D("h_nue_dir_x_eff_num", "h_nue_dir_x_eff_num", 20, -1, 1);
TH1D * h_nue_dir_y_eff_den         = new TH1D("h_nue_dir_y_eff_den", "h_nue_dir_y_eff_den", 20, -1, 1);
TH1D * h_nue_dir_y_eff_num         = new TH1D("h_nue_dir_y_eff_num", "h_nue_dir_y_eff_num", 20, -1, 1);
TH1D * h_nue_dir_z_eff_den         = new TH1D("h_nue_dir_z_eff_den", "h_nue_dir_z_eff_den", 20, -1, 1);
TH1D * h_nue_dir_z_eff_num         = new TH1D("h_nue_dir_z_eff_num", "h_nue_dir_z_eff_num", 20, -1, 1);
TH1D * h_ele_dir_x_eff_den         = new TH1D("h_ele_dir_x_eff_den", "h_ele_dir_x_eff_den", 20, -1, 1);
TH1D * h_ele_dir_x_eff_num         = new TH1D("h_ele_dir_x_eff_num", "h_ele_dir_x_eff_num", 20, -1, 1);
TH1D * h_ele_dir_y_eff_den         = new TH1D("h_ele_dir_y_eff_den", "h_ele_dir_y_eff_den", 20, -1, 1);
TH1D * h_ele_dir_y_eff_num         = new TH1D("h_ele_dir_y_eff_num", "h_ele_dir_y_eff_num", 20, -1, 1);
TH1D * h_ele_dir_z_eff_den         = new TH1D("h_ele_dir_z_eff_den", "h_ele_dir_z_eff_den", 20, -1, 1);
TH1D * h_ele_dir_z_eff_num         = new TH1D("h_ele_dir_z_eff_num", "h_ele_dir_z_eff_num", 20, -1, 1);
TH1D * h_nue_num_part_eff_den      = new TH1D("h_nue_num_part_eff_den", "h_nue_num_part_eff_den", 20, 0, 20);
TH1D * h_nue_num_part_eff_num      = new TH1D("h_nue_num_part_eff_num", "h_nue_num_part_eff_num", 20, 0, 20);
TH1D * h_nue_num_chrg_part_eff_den = new TH1D("h_nue_num_chrg_part_eff_den", "h_nue_num_chrg_part_eff_den", 20, 0, 20);
TH1D * h_nue_num_chrg_part_eff_num = new TH1D("h_nue_num_chrg_part_eff_num", "h_nue_num_chrg_part_eff_num", 20, 0, 20);
TH1D * h_nue_cos_theta_eff_den     = new TH1D("h_nue_cos_theta_eff_den", "h_nue_cos_theta_eff_den", 10, -1, 1);
TH1D * h_nue_cos_theta_eff_num     = new TH1D("h_nue_cos_theta_eff_num", "h_nue_cos_theta_eff_num", 10, -1, 1);
TH1D * h_nue_phi_eff_den           = new TH1D("h_nue_phi_eff_den", "h_nue_phi_eff_den", 60, -180, 180);
TH1D * h_nue_phi_eff_num           = new TH1D("h_nue_phi_eff_num", "h_nue_phi_eff_num", 60, -180, 180);
TH1D * h_ele_cos_theta_eff_den     = new TH1D("h_ele_cos_theta_eff_den", "h_ele_cos_theta_eff_den", 10, -1, 1);
TH1D * h_ele_cos_theta_eff_num     = new TH1D("h_ele_cos_theta_eff_num", "h_ele_cos_theta_eff_num", 10, -1, 1);
TH1D * h_ele_phi_eff_den           = new TH1D("h_ele_phi_eff_den", "h_ele_phi_eff_den", 60, -180, 180);
TH1D * h_ele_phi_eff_num           = new TH1D("h_ele_phi_eff_num", "h_ele_phi_eff_num", 60, -180, 180);

TH2I * h_tracks_showers         = new TH2I("h_tracks_showers", "h_tracks_showers", 8, 0, 8, 8, 0, 8);
TH2I * h_tracks_showers_cosmic  = new TH2I("h_tracks_showers_cosmic", "h_tracks_showers_cosmic", 8, 0, 8, 8, 0, 8);
TH2I * h_tracks_showers_numu    = new TH2I("h_tracks_showers_numu", "h_tracks_showers_numu", 8, 0, 8, 8, 0, 8);

TH1D * h_leading_shower_open_angle_nue_cc        = new TH1D("h_leading_shower_open_angle_nue_cc",        "h_leading_shower_open_angle_nue_cc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_nue_cc_mixed  = new TH1D("h_leading_shower_open_angle_nue_cc_mixed",  "h_leading_shower_open_angle_nue_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_numu_cc       = new TH1D("h_leading_shower_open_angle_numu_cc",       "h_leading_shower_open_angle_numu_cc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_nc            = new TH1D("h_leading_shower_open_angle_nc",            "h_leading_shower_open_angle_nc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_cosmic        = new TH1D("h_leading_shower_open_angle_cosmic",        "h_leading_shower_open_angle_cosmic", 25, 0, 50);
TH1D * h_leading_shower_open_angle_nc_pi0        = new TH1D("h_leading_shower_open_angle_nc_pi0",        "h_leading_shower_open_angle_nc_pi0", 25, 0, 50);
TH1D * h_leading_shower_open_angle_numu_cc_mixed = new TH1D("h_leading_shower_open_angle_numu_cc_mixed", "h_leading_shower_open_angle_numu_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_other_mixed   = new TH1D("h_leading_shower_open_angle_other_mixed",   "h_leading_shower_open_angle_other_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_unmatched     = new TH1D("h_leading_shower_open_angle_unmatched",     "h_leading_shower_open_angle_unmatched", 25, 0, 50);

TH1D * h_leading_shower_open_angle_1_nue_cc        = new TH1D("h_leading_shower_open_angle_1_nue_cc",        "h_leading_shower_open_angle_1_nue_cc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_nue_cc_mixed  = new TH1D("h_leading_shower_open_angle_1_nue_cc_mixed",  "h_leading_shower_open_angle_1_nue_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_numu_cc       = new TH1D("h_leading_shower_open_angle_1_numu_cc",       "h_leading_shower_open_angle_1_numu_cc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_nc            = new TH1D("h_leading_shower_open_angle_1_nc",            "h_leading_shower_open_angle_1_nc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_cosmic        = new TH1D("h_leading_shower_open_angle_1_cosmic",        "h_leading_shower_open_angle_1_cosmic", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_nc_pi0        = new TH1D("h_leading_shower_open_angle_1_nc_pi0",        "h_leading_shower_open_angle_1_nc_pi0", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_numu_cc_mixed = new TH1D("h_leading_shower_open_angle_1_numu_cc_mixed", "h_leading_shower_open_angle_1_numu_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_other_mixed   = new TH1D("h_leading_shower_open_angle_1_other_mixed",   "h_leading_shower_open_angle_1_other_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_1_unmatched     = new TH1D("h_leading_shower_open_angle_1_unmatched",     "h_leading_shower_open_angle_1_unmatched", 25, 0, 50);

TH1D * h_leading_shower_open_angle_2plus_nue_cc        = new TH1D("h_leading_shower_open_angle_2plus_nue_cc",        "h_leading_shower_open_angle_2plus_nue_cc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_nue_cc_mixed  = new TH1D("h_leading_shower_open_angle_2plus_nue_cc_mixed",  "h_leading_shower_open_angle_2plus_nue_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_numu_cc       = new TH1D("h_leading_shower_open_angle_2plus_numu_cc",       "h_leading_shower_open_angle_2plus_numu_cc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_nc            = new TH1D("h_leading_shower_open_angle_2plus_nc",            "h_leading_shower_open_angle_2plus_nc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_cosmic        = new TH1D("h_leading_shower_open_angle_2plus_cosmic",        "h_leading_shower_open_angle_2plus_cosmic", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_nc_pi0        = new TH1D("h_leading_shower_open_angle_2plus_nc_pi0",        "h_leading_shower_open_angle_2plus_nc_pi0", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_numu_cc_mixed = new TH1D("h_leading_shower_open_angle_2plus_numu_cc_mixed", "h_leading_shower_open_angle_2plus_numu_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_other_mixed   = new TH1D("h_leading_shower_open_angle_2plus_other_mixed",   "h_leading_shower_open_angle_2plus_other_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_2plus_unmatched     = new TH1D("h_leading_shower_open_angle_2plus_unmatched",     "h_leading_shower_open_angle_2plus_unmatched", 25, 0, 50);

TH1D * h_trk_vtx_dist_nue_cc        = new TH1D("h_trk_vtx_dist_nue_cc", "h_trk_vtx_dist_nue_cc", 25, 0, 20);
TH1D * h_trk_vtx_dist_nue_cc_mixed  = new TH1D("h_trk_vtx_dist_nue_cc_mixed", "h_trk_vtx_dist_nue_cc_mixed", 25, 0, 20);
TH1D * h_trk_vtx_dist_numu_cc       = new TH1D("h_trk_vtx_dist_numu_cc", "h_trk_vtx_dist_numu_cc", 25, 0, 20);
TH1D * h_trk_vtx_dist_nc            = new TH1D("h_trk_vtx_dist_nc", "h_trk_vtx_dist_nc", 25, 0, 20);
TH1D * h_trk_vtx_dist_cosmic        = new TH1D("h_trk_vtx_dist_cosmic", "h_trk_vtx_dist_cosmic", 25, 0, 20);
TH1D * h_trk_vtx_dist_nc_pi0        = new TH1D("h_trk_vtx_dist_nc_pi0", "h_trk_vtx_dist_nc_pi0", 25, 0, 20);
TH1D * h_trk_vtx_dist_numu_cc_mixed = new TH1D("h_trk_vtx_dist_numu_cc_mixed", "h_trk_vtx_dist_numu_cc_mixed", 25, 0, 20);
TH1D * h_trk_vtx_dist_other_mixed   = new TH1D("h_trk_vtx_dist_other_mixed", "h_trk_vtx_dist_other_mixed", 25, 0, 20);
TH1D * h_trk_vtx_dist_unmatched     = new TH1D("h_trk_vtx_dist_unmatched", "h_trk_vtx_dist_unmatched", 25, 0, 20);

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
TH1D * h_vtx_flash_numu_cc       = new TH1D("h_vtx_flash_numu_cc",        "h_vtx_flash_numu_cc_mixed", 40, 0, 200);
TH1D * h_vtx_flash_nc_pi0        = new TH1D("h_vtx_flash_nc_pi0",         "h_vtx_flash_nc_pi0",        40, 0, 200);
TH1D * h_vtx_flash_cosmic        = new TH1D("h_vtx_flash_cosmic",         "h_vtx_flash_cosmic",        40, 0, 200);
TH1D * h_vtx_flash_nc            = new TH1D("h_vtx_flash_nc",             "h_vtx_flash_nc",            40, 0, 200);
TH1D * h_vtx_flash_numu_cc_mixed = new TH1D("h_vtx_flash_numu_cc_mixed",  "h_vtx_flash_numu_cc_mixed", 40, 0, 200);
TH1D * h_vtx_flash_other_mixed   = new TH1D("h_vtx_flash_other_mixed",    "h_vtx_flash_other_mixed",   40, 0, 200);
TH1D * h_vtx_flash_unmatched     = new TH1D("h_vtx_flash_unmatched",      "h_vtx_flash_unmatched",     40, 0, 200);

TH1D * h_shwr_vtx_dist_nue_cc        = new TH1D("h_shwr_vtx_dist_nue_cc",         "h_shwr_vtx_dist_nue_cc",        20, 0, 20);
TH1D * h_shwr_vtx_dist_nue_cc_mixed  = new TH1D("h_shwr_vtx_dist_nue_cc_mixed",   "h_shwr_vtx_dist_nue_cc_mixed",  20, 0, 20);
TH1D * h_shwr_vtx_dist_numu_cc       = new TH1D("h_shwr_vtx_dist_numu_cc",        "h_shwr_vtx_dist_numu_cc",       20, 0, 20);
TH1D * h_shwr_vtx_dist_nc_pi0        = new TH1D("h_shwr_vtx_dist_nc_pi0",         "h_shwr_vtx_dist_nc_pi0",        20, 0, 20);
TH1D * h_shwr_vtx_dist_cosmic        = new TH1D("h_shwr_vtx_dist_cosmic",         "h_shwr_vtx_dist_cosmic",        20, 0, 20);
TH1D * h_shwr_vtx_dist_nc            = new TH1D("h_shwr_vtx_dist_nc",             "h_shwr_vtx_dist_nc",            20, 0, 20);
TH1D * h_shwr_vtx_dist_numu_cc_mixed = new TH1D("h_shwr_vtx_dist_numu_cc_mixed",  "h_shwr_vtx_dist_numu_cc_mixed", 20, 0, 20);
TH1D * h_shwr_vtx_dist_other_mixed   = new TH1D("h_shwr_vtx_dist_other_mixed",    "h_shwr_vtx_dist_other_mixed",   20, 0, 20);
TH1D * h_shwr_vtx_dist_unmatched     = new TH1D("h_shwr_vtx_dist_unmatched",      "h_shwr_vtx_dist_unmatched",     20, 0, 20);

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

TH1D * h_failure_reason_nue_cc         = new TH1D ("h_failure_reason_nue_cc",         "h_failure_reason_nue_cc",         22, 0, 11);
TH1D * h_failure_reason_nue_cc_out_fv  = new TH1D ("h_failure_reason_nue_cc_out_fv",  "h_failure_reason_nue_cc_out_fv",  22, 0, 11);
TH1D * h_failure_reason_nue_cc_mixed   = new TH1D ("h_failure_reason_nue_cc_mixed",   "h_failure_reason_nue_cc_mixed",   22, 0, 11);
TH1D * h_failure_reason_numu_cc        = new TH1D ("h_failure_reason_numu_cc",        "h_failure_reason_numu_cc",        22, 0, 11);
TH1D * h_failure_reason_numu_cc_mixed  = new TH1D ("h_failure_reason_numu_cc_mixed",  "h_failure_reason_numu_cc_mixed",  22, 0, 11);
TH1D * h_failure_reason_nc             = new TH1D ("h_failure_reason_nc",             "h_failure_reason_nc",             22, 0, 11);
TH1D * h_failure_reason_nc_pi0         = new TH1D ("h_failure_reason_nc_pi0",         "h_failure_reason_nc_pi0",         22, 0, 11);
TH1D * h_failure_reason_cosmic         = new TH1D ("h_failure_reason_cosmic",         "h_failure_reason_cosmic",         22, 0, 11);
TH1D * h_failure_reason_other_mixed    = new TH1D ("h_failure_reason_other_mixed",    "h_failure_reason_other_mixed",    22, 0, 11);
TH1D * h_failure_reason_unmatched      = new TH1D ("h_failure_reason_unmatched",      "h_failure_reason_unmatched",      22, 0, 11);

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


}//end namespace

#endif
