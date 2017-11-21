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
//standard is: 0, 3.5
const double tolerance_dedx_min = 1.4;
const double tolerance_dedx_max = 3;

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
int reco_nue_counter_nue_cc_qe  = 0;
int reco_nue_counter_nue_cc_res = 0;
int reco_nue_counter_nue_cc_dis = 0;
int reco_nue_counter_nue_cc_coh = 0;
int reco_nue_counter_nue_cc_mec = 0;
int reco_nue_counter_numu_cc_qe  = 0;
int reco_nue_counter_numu_cc_res = 0;
int reco_nue_counter_numu_cc_dis = 0;
int reco_nue_counter_numu_cc_coh = 0;
int reco_nue_counter_numu_cc_mec = 0;
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
int in_fv_counter_nue_cc_qe  = 0;
int in_fv_counter_nue_cc_res = 0;
int in_fv_counter_nue_cc_dis = 0;
int in_fv_counter_nue_cc_coh = 0;
int in_fv_counter_nue_cc_mec = 0;
int in_fv_counter_numu_cc_qe  = 0;
int in_fv_counter_numu_cc_res = 0;
int in_fv_counter_numu_cc_dis = 0;
int in_fv_counter_numu_cc_coh = 0;
int in_fv_counter_numu_cc_mec = 0;
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
int vtx_flash_counter_nue_cc_qe  = 0;
int vtx_flash_counter_nue_cc_res = 0;
int vtx_flash_counter_nue_cc_dis = 0;
int vtx_flash_counter_nue_cc_coh = 0;
int vtx_flash_counter_nue_cc_mec = 0;
int vtx_flash_counter_numu_cc_qe  = 0;
int vtx_flash_counter_numu_cc_res = 0;
int vtx_flash_counter_numu_cc_dis = 0;
int vtx_flash_counter_numu_cc_coh = 0;
int vtx_flash_counter_numu_cc_mec = 0;
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
int shwr_tpco_counter_nue_cc_qe  = 0;
int shwr_tpco_counter_nue_cc_res = 0;
int shwr_tpco_counter_nue_cc_dis = 0;
int shwr_tpco_counter_nue_cc_coh = 0;
int shwr_tpco_counter_nue_cc_mec = 0;
int shwr_tpco_counter_numu_cc_qe  = 0;
int shwr_tpco_counter_numu_cc_res = 0;
int shwr_tpco_counter_numu_cc_dis = 0;
int shwr_tpco_counter_numu_cc_coh = 0;
int shwr_tpco_counter_numu_cc_mec = 0;
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
int trk_tpco_counter_nue_cc_qe  = 0;
int trk_tpco_counter_nue_cc_res = 0;
int trk_tpco_counter_nue_cc_dis = 0;
int trk_tpco_counter_nue_cc_coh = 0;
int trk_tpco_counter_nue_cc_mec = 0;
int trk_tpco_counter_numu_cc_qe  = 0;
int trk_tpco_counter_numu_cc_res = 0;
int trk_tpco_counter_numu_cc_dis = 0;
int trk_tpco_counter_numu_cc_coh = 0;
int trk_tpco_counter_numu_cc_mec = 0;
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
int hit_threshold_counter_nue_cc_qe  = 0;
int hit_threshold_counter_nue_cc_res = 0;
int hit_threshold_counter_nue_cc_dis = 0;
int hit_threshold_counter_nue_cc_coh = 0;
int hit_threshold_counter_nue_cc_mec = 0;
int hit_threshold_counter_numu_cc_qe  = 0;
int hit_threshold_counter_numu_cc_res = 0;
int hit_threshold_counter_numu_cc_dis = 0;
int hit_threshold_counter_numu_cc_coh = 0;
int hit_threshold_counter_numu_cc_mec = 0;
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
int open_angle_counter_nue_cc_qe  = 0;
int open_angle_counter_nue_cc_res = 0;
int open_angle_counter_nue_cc_dis = 0;
int open_angle_counter_nue_cc_coh = 0;
int open_angle_counter_nue_cc_mec = 0;
int open_angle_counter_numu_cc_qe  = 0;
int open_angle_counter_numu_cc_res = 0;
int open_angle_counter_numu_cc_dis = 0;
int open_angle_counter_numu_cc_coh = 0;
int open_angle_counter_numu_cc_mec = 0;
int dedx_counter = 0;
int dedx_counter_nue_cc = 0;
int dedx_counter_nue_cc_mixed = 0;
int dedx_counter_nue_cc_out_fv = 0;
int dedx_counter_cosmic = 0;
int dedx_counter_nue_nc = 0;
int dedx_counter_numu_cc = 0;
int dedx_counter_numu_cc_mixed = 0;
int dedx_counter_numu_nc = 0;
int dedx_counter_unmatched = 0;
int dedx_counter_other_mixed = 0;
int dedx_counter_nue_cc_qe  = 0;
int dedx_counter_nue_cc_res = 0;
int dedx_counter_nue_cc_dis = 0;
int dedx_counter_nue_cc_coh = 0;
int dedx_counter_nue_cc_mec = 0;
int dedx_counter_numu_cc_qe  = 0;
int dedx_counter_numu_cc_res = 0;
int dedx_counter_numu_cc_dis = 0;
int dedx_counter_numu_cc_coh = 0;
int dedx_counter_numu_cc_mec = 0;
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

TH1D * h_leading_shower_open_angle_nue_cc        = new TH1D("h_leading_shower_open_angle_nue_cc", "h_leading_shower_open_angle_nue_cc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_nue_cc_mixed  = new TH1D("h_leading_shower_open_angle_nue_cc_mixed", "h_leading_shower_open_angle_nue_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_numu_cc       = new TH1D("h_leading_shower_open_angle_numu_cc", "h_leading_shower_open_angle_numu_cc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_numu_nc       = new TH1D("h_leading_shower_open_angle_numu_nc", "h_leading_shower_open_angle_numu_nc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_cosmic        = new TH1D("h_leading_shower_open_angle_cosmic", "h_leading_shower_open_angle_cosmic", 25, 0, 50);
TH1D * h_leading_shower_open_angle_nue_nc        = new TH1D("h_leading_shower_open_angle_nue_nc", "h_leading_shower_open_angle_nue_nc", 25, 0, 50);
TH1D * h_leading_shower_open_angle_numu_cc_mixed = new TH1D("h_leading_shower_open_angle_numu_cc_mixed", "h_leading_shower_open_angle_numu_cc_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_other_mixed   = new TH1D("h_leading_shower_open_angle_other_mixed", "h_leading_shower_open_angle_other_mixed", 25, 0, 50);
TH1D * h_leading_shower_open_angle_unmatched     = new TH1D("h_leading_shower_open_angle_unmatched", "h_leading_shower_open_angle_unmatched", 25, 0, 50);

TH1D * h_trk_vtx_dist_nue_cc        = new TH1D("h_trk_vtx_dist_nue_cc", "h_trk_vtx_dist_nue_cc", 25, 0, 20);
TH1D * h_trk_vtx_dist_nue_cc_mixed  = new TH1D("h_trk_vtx_dist_nue_cc_mixed", "h_trk_vtx_dist_nue_cc_mixed", 25, 0, 20);
TH1D * h_trk_vtx_dist_numu_cc       = new TH1D("h_trk_vtx_dist_numu_cc", "h_trk_vtx_dist_numu_cc", 25, 0, 20);
TH1D * h_trk_vtx_dist_numu_nc       = new TH1D("h_trk_vtx_dist_numu_nc", "h_trk_vtx_dist_numu_nc", 25, 0, 20);
TH1D * h_trk_vtx_dist_cosmic        = new TH1D("h_trk_vtx_dist_cosmic", "h_trk_vtx_dist_cosmic", 25, 0, 20);
TH1D * h_trk_vtx_dist_nue_nc        = new TH1D("h_trk_vtx_dist_nue_nc", "h_trk_vtx_dist_nue_nc", 25, 0, 20);
TH1D * h_trk_vtx_dist_numu_cc_mixed = new TH1D("h_trk_vtx_dist_numu_cc_mixed", "h_trk_vtx_dist_numu_cc_mixed", 25, 0, 20);
TH1D * h_trk_vtx_dist_other_mixed   = new TH1D("h_trk_vtx_dist_other_mixed", "h_trk_vtx_dist_other_mixed", 25, 0, 20);
TH1D * h_trk_vtx_dist_unmatched     = new TH1D("h_trk_vtx_dist_unmatched", "h_trk_vtx_dist_unmatched", 25, 0, 20);

TH2D * h_pfp_track_shower_nue_cc_qe     = new TH2D("h_pfp_track_shower_nue_cc_qe",     "h_pfp_track_shower_nue_cc_qe", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_out_fv = new TH2D("h_pfp_track_shower_nue_cc_out_fv", "h_pfp_track_shower_nue_cc_out_fv", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_res    = new TH2D("h_pfp_track_shower_nue_cc_res",    "h_pfp_track_shower_nue_cc_res", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_dis    = new TH2D("h_pfp_track_shower_nue_cc_dis",    "h_pfp_track_shower_nue_cc_dis", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_coh    = new TH2D("h_pfp_track_shower_nue_cc_coh",    "h_pfp_track_shower_nue_cc_coh", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_cc_mec    = new TH2D("h_pfp_track_shower_nue_cc_mec",    "h_pfp_track_shower_nue_cc_mec", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_nue_nc        = new TH2D("h_pfp_track_shower_nue_nc",        "h_pfp_track_shower_nue_nc", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_qe    = new TH2D("h_pfp_track_shower_numu_cc_qe",    "h_pfp_track_shower_numu_cc_qe", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_res   = new TH2D("h_pfp_track_shower_numu_cc_res",   "h_pfp_track_shower_numu_cc_res", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_dis   = new TH2D("h_pfp_track_shower_numu_cc_dis",   "h_pfp_track_shower_numu_cc_dis", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_coh   = new TH2D("h_pfp_track_shower_numu_cc_coh",   "h_pfp_track_shower_numu_cc_coh", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_cc_mec   = new TH2D("h_pfp_track_shower_numu_cc_mec",   "h_pfp_track_shower_numu_cc_mec", 10, 0, 10, 10, 0, 10);
TH2D * h_pfp_track_shower_numu_nc       = new TH2D("h_pfp_track_shower_numu_nc",       "h_pfp_track_shower_numu_nc", 10, 0, 10, 10, 0, 10);
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
TH2D * h_leading_shower_mc_pdg_nue_nc        = new TH2D("h_leading_shower_mc_pdg_nue_nc",        "h_leading_shower_mc_pdg_nue_nc",        3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_qe    = new TH2D("h_leading_shower_mc_pdg_numu_cc_qe",    "h_leading_shower_mc_pdg_numu_cc_qe",    3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_res   = new TH2D("h_leading_shower_mc_pdg_numu_cc_res",   "h_leading_shower_mc_pdg_numu_cc_res",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_dis   = new TH2D("h_leading_shower_mc_pdg_numu_cc_dis",   "h_leading_shower_mc_pdg_numu_cc_dis",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_coh   = new TH2D("h_leading_shower_mc_pdg_numu_cc_coh",   "h_leading_shower_mc_pdg_numu_cc_coh",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_cc_mec   = new TH2D("h_leading_shower_mc_pdg_numu_cc_mec",   "h_leading_shower_mc_pdg_numu_cc_mec",   3, 0, 3, 10, 0, 10);
TH2D * h_leading_shower_mc_pdg_numu_nc       = new TH2D("h_leading_shower_mc_pdg_numu_nc",       "h_leading_shower_mc_pdg_numu_nc",       3, 0, 3, 10, 0, 10);
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
TH1D * h_pfp_track_nue_nc        = new TH1D("h_pfp_track_nue_nc",        "h_pfp_track_nue_nc",        10, 0, 10);
TH1D * h_pfp_track_numu_cc_qe    = new TH1D("h_pfp_track_numu_cc_qe",    "h_pfp_track_numu_cc_qe",    10, 0, 10);
TH1D * h_pfp_track_numu_cc_res   = new TH1D("h_pfp_track_numu_cc_res",   "h_pfp_track_numu_cc_res",   10, 0, 10);
TH1D * h_pfp_track_numu_cc_dis   = new TH1D("h_pfp_track_numu_cc_dis",   "h_pfp_track_numu_cc_dis",   10, 0, 10);
TH1D * h_pfp_track_numu_cc_coh   = new TH1D("h_pfp_track_numu_cc_coh",   "h_pfp_track_numu_cc_coh",   10, 0, 10);
TH1D * h_pfp_track_numu_cc_mec   = new TH1D("h_pfp_track_numu_cc_mec",   "h_pfp_track_numu_cc_mec",   10, 0, 10);
TH1D * h_pfp_track_numu_nc       = new TH1D("h_pfp_track_numu_nc",       "h_pfp_track_numu_nc",       10, 0, 10);
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
TH1D * h_pfp_shower_nue_nc        = new TH1D("h_pfp_shower_nue_nc",        "h_pfp_shower_nue_nc",        10, 0, 10);
TH1D * h_pfp_shower_numu_cc_qe    = new TH1D("h_pfp_shower_numu_cc_qe",    "h_pfp_shower_numu_cc_qe",    10, 0, 10);
TH1D * h_pfp_shower_numu_cc_res   = new TH1D("h_pfp_shower_numu_cc_res",   "h_pfp_shower_numu_cc_res",   10, 0, 10);
TH1D * h_pfp_shower_numu_cc_dis   = new TH1D("h_pfp_shower_numu_cc_dis",   "h_pfp_shower_numu_cc_dis",   10, 0, 10);
TH1D * h_pfp_shower_numu_cc_coh   = new TH1D("h_pfp_shower_numu_cc_coh",   "h_pfp_shower_numu_cc_coh",   10, 0, 10);
TH1D * h_pfp_shower_numu_cc_mec   = new TH1D("h_pfp_shower_numu_cc_mec",   "h_pfp_shower_numu_cc_mec",   10, 0, 10);
TH1D * h_pfp_shower_numu_nc       = new TH1D("h_pfp_shower_numu_nc",       "h_pfp_shower_numu_nc",       10, 0, 10);
TH1D * h_pfp_shower_nue_cc_mixed  = new TH1D("h_pfp_shower_nue_cc_mixed",  "h_pfp_shower_nue_cc_mixed",  10, 0, 10);
TH1D * h_pfp_shower_numu_cc_mixed = new TH1D("h_pfp_shower_numu_cc_mixed", "h_pfp_shower_numu_cc_mixed", 10, 0, 10);
TH1D * h_pfp_shower_cosmic        = new TH1D("h_pfp_shower_cosmic",        "h_pfp_shower_cosmic",        10, 0, 10);
TH1D * h_pfp_shower_other_mixed   = new TH1D("h_pfp_shower_other_mixed",   "h_pfp_shower_other_mixed",   10, 0, 10);
TH1D * h_pfp_shower_unmatched     = new TH1D("h_pfp_shower_unmatched",     "h_pfp_shower_unmatched",     10, 0, 10);

}//end namespace

#endif
