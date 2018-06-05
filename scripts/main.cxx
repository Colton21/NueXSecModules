#include "main.h"

int main(int argc, char *argv[]){

	const char * file1 = argv[1];//mc
	const char * file2 = argv[2];//ext
	const char * file3 = argv[3];//on-beam data
	const char * input_config_file_name;

	bool using_default_config = true;

	//start after name of .exe
	for(int i =1; i < argc; i++)
	{
		auto const arg = argv[i];
		std::cout << arg << std::endl;
		if(strcmp(arg, "-c") == 0)
		{
			using_default_config = false;
			input_config_file_name = argv[i+1];
			break;
		}
	}

	std::cout << "INPUT FORMAT: MC_FILE INTIME_FILE DATA_FILE" << std::endl;
	if(argc < 2 )  { std::cout << " \n Please inclue the input file path \n " << std::endl; exit(1); }

	std::vector<double> config;
	std::vector<double> input_config;

	xsecSelection::selection _selection_instance;

	std::vector<double> default_config = utility::configure_cuts(
	        _x1, _x2, _y1, _y2, _z1, _z2,
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

	std::ifstream input_config_file;
	//check if string for input file is empty -
	//if not empty try to open file
	if(using_default_config == false)
	{
		input_config_file.open(input_config_file_name);
		if(!input_config_file.is_open())
		{
			std::cout << "*** \t File did not open! \t ***" << std::endl;
			using_default_config = true;
		}
		if(input_config_file.is_open())
		{
			std::cout << "*** \t File Opened \t *** " << std::endl;
			std::string line;
			while (!input_config_file.eof())
			{
				std::getline (input_config_file, line);
				if(!line.empty())
				{
					input_config.push_back(std::stof(line));
				}
			}
			input_config_file.close();
			std::cout << "*** \t File Closed \t *** " << std::endl;
		}
	}

	if(using_default_config == true)  {config = default_config; }
	if(using_default_config == false) {config = input_config; }

	std::vector<std::tuple<double, double, std::string> > * results_v = new std::vector<std::tuple<double, double, std::string> >;

	if(using_default_config == true)
	{
		std::cout << "\n --- Using Default Configuration --- \n" << std::endl;
		if(argc == 2)
		{
			std::cout << "Running without in-time cosmics " << std::endl;
			std::cout << "Running without data" << std::endl;
			_selection_instance.xsecSelection::selection::make_selection(file1, "empty", "empty", config, results_v);
			return 0;
		}
		if(argc == 3)
		{
			std::cout << "Running without data " << std::endl;
			_selection_instance.xsecSelection::selection::make_selection(file1, file2, "empty", config, results_v);
			return 0;
		}
		if(argc == 4)
		{
			std::cout << "Running with MC, EXT, and Data" << std::endl;
			_selection_instance.xsecSelection::selection::make_selection(file1, file2, file3, config, results_v);
			return 0;
		}
	}

	if(using_default_config == false)
	{
		std::cout << "\n --- Using Input Configuration --- \n" << std::endl;
		if(argc == 4)
		{
			std::cout << "Running without in-time cosmics " << std::endl;
			std::cout << "Running without data" << std::endl;
			_selection_instance.xsecSelection::selection::make_selection(file1, "empty", "empty", config, results_v);
			return 0;
		}
		if(argc == 5)
		{
			std::cout << "Running without data " << std::endl;
			_selection_instance.xsecSelection::selection::make_selection(file1, file2, "empty", config, results_v);
			return 0;
		}
		if(argc == 6)
		{
			std::cout << "Running with MC, EXT, and Data" << std::endl;
			_selection_instance.xsecSelection::selection::make_selection(file1, file2, file3, config, results_v);
			for(auto const results : * results_v) {std::cout << std::get<0>(results) << ", "
				                                         << std::get<1>(results) << ", "
				                                         << std::get<2>(results) << std::endl; }
			return 0;
		}
	}

	std::cout << "Returned without running..." << std::endl;
	return 0;
}
