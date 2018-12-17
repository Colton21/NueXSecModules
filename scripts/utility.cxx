#include "utility.h"

namespace utility {

std::vector<double> configure_cuts(
        double _x1,
        double _x2,
        double _y1,
        double _y2,
        double _z1,
        double _z2,
        double flash_pe_threshold,
        double flash_time_start,
        double flash_time_end,
        double tolerance,
        double shwr_nue_tolerance,
        double trk_nue_tolerance,
        double shwr_hit_threshold,
        double shwr_hit_threshold_collection,
        double tolerance_open_angle_min,
        double tolerance_open_angle_max,
        double tolerance_dedx_min,
        double tolerance_dedx_max,
        double dist_tolerance,
        double pfp_hits_length_tolerance,
        double ratio_tolerance,
        bool do_variations
        )
{
	std::vector<double> config;
	config.resize(22,0);

	config[0] = _x1;
	config[1] = _x2;
	config[2] = _y1;
	config[3] = _y2;
	config[4] = _z1;
	config[5] = _z2;
	config[6] = flash_pe_threshold;
	config[7] = flash_time_start;
	config[8] = flash_time_end;
	config[9] = tolerance;
	config[10] = shwr_nue_tolerance;
	config[11] = trk_nue_tolerance;
	config[12] = shwr_hit_threshold;
	config[13] = shwr_hit_threshold_collection;
	config[14] = tolerance_open_angle_min;
	config[15] = tolerance_open_angle_max;
	config[16] = tolerance_dedx_min;
	config[17] = tolerance_dedx_max;
	config[18] = dist_tolerance;
	config[19] = pfp_hits_length_tolerance;
	config[20] = ratio_tolerance;
	config[21] = do_variations;

	return config;
}//end config function
}//end namespace utlity
