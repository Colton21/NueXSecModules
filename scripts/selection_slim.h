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

const double POT = 1.82027e21;

//const double POT = 4.05982e+19;      //POT - all NuMI + cosmics
//const double POT = 1.23206e+20; //POT - all NuMI + cosmics, bigger sample
//const double POT = 2.90469e+21;    //POT - nue + cosmics
const double scaling_nue = 1.52938e-11;    //nues / POT / cm^2
const double scaling_nue_bar = 7.77111e-12;   //anues / POT / cm^2
const double scaling = (scaling_nue + scaling_nue_bar);

//since I have both nue and nue_bar as signal definition need to adjust for this
const double genie_xsec_nue = 5.63067e-39;   //cm^2
const double genie_xsec_nue_bar = 2.0893e-39;   //cm^2
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
