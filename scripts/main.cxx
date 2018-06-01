#include "main.h"

int main(int argc, char *argv[]){

	const char * file1 = argv[1];
	const char * file2 = argv[2];
	const char * file3 = argv[3];

	std::cout << "INPUT FORMAT: MC_FILE INTIME_FILE DATA_FILE" << std::endl;

	std::vector<double> config = utility::configure_cuts(
	        _x1,
	        _x2,
	        _y1,
	        _y2,
	        _z1,
	        _z2,
	        flash_pe_threshold,
	        flash_time_start,
	        flash_time_end,
	        tolerance,
	        shwr_nue_tolerance,
	        trk_nue_tolerance,
	        shwr_hit_threshold,
	        shwr_hit_threshold_collection,
	        tolerance_open_angle_min,
	        tolerance_open_angle_max,
	        tolerance_dedx_min,
	        tolerance_dedx_max,
	        dist_tolerance,
	        pfp_hits_length_tolerance,
	        ratio_tolerance
	        );

	xsecSelection::selection _selection_instance;
	if(argc < 2 )  { std::cout << "Please inclue the input file path" << std::endl; exit(1); }

	if(argc != 3 && argc != 4)
	{
		std::cout << "Running without in-time cosmics " << std::endl;
		std::cout << "Running without data" << std::endl;
		_selection_instance.xsecSelection::selection::make_selection(file1, "empty", "empty", config);
		return 0;
	}
	if(argc != 4 )
	{
		std::cout << "Running without in-time data " << std::endl;
		_selection_instance.xsecSelection::selection::make_selection(file1, file2, "empty", config);
		return 0;
	}

	_selection_instance.xsecSelection::selection::make_selection(file1, file2, file3, config);
	return 0;
}
