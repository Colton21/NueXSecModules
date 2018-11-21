#ifndef SELECTION_SLIM_h
#define SELECTION_SLIM_h

#include "selection_functions.h"
#include "selection_cuts.h"
#include "selection_functions_data.h"
#include "histogram_functions.h"

#include "../xsecAna/LinkDef.h"

#include <fstream>

namespace xsecSelection {

class selection_slim {

private:

//values are old - need to be updated - circa Nov 2018
const double POT = 1.82027e+21;
const double data_scale_factor = 1 / 5.648;   //ie scale MC down by factor
//const double intime_scale_factor = 0.442416; //ie scale EXT down by factor
const double intime_scale_factor = 0.56940408;   //for first two datasets

//dirt scaling
const double dirt_scale_factor = 0;

//these are for the flux calculations
const double scaling_nue = 1.52938e-11;          //nues  / POT / cm^2
const double scaling_nue_bar = 7.77111e-12;      //anues / POT / cm^2
const double scaling = scaling_nue + scaling_nue_bar;

const double genie_xsec_nue = 5.63067e-39;     //cm^2
const double genie_xsec_nue_bar = 2.0893e-39;     //cm^2
//const double genie_xsec = genie_xsec_nue + genie_xsec_nue_bar;

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

selection_slim() = default;

void make_selection_slim(
        const char * _file1,
        const char * _file2,
        const char * _file3,
        const char * _file4,
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

}; //end class

}//end namespace

#endif
