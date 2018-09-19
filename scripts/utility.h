#ifndef UTILITY_h
#define UTILITY_h

#include <vector>

namespace utility {

std::vector<double> configure_cuts ( double _x1,
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
                                     double ratio_tolerance
                                     );
}//end namespace

#endif
