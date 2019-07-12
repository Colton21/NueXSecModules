#include "main.h"

int main(int argc, char *argv[]){

	const char * input_config_file_name;
	const char * input_file_locate_prefix;

	bool using_default_config = true;
	bool using_slim_version = false;
	bool using_dynamic_file_prefix = false;

	bool monte_carlo_file = false;
	bool cosmic_file = false;
	bool data_file = false;
	bool dirt_file = false;
	bool variation_file = false;

	char * monte_carlo_file_path = (char *)"empty";
	char * cosmic_file_path = (char *)"empty";
	char * data_file_path = (char *)"empty";
	char * dirt_file_path = (char *)"empty";
	char * variation_file_path = (char *)"empty";

	bool variation_mode;
	bool variation_mode_bkg = false;

	//start after name of .exe
	for(int i =1; i < argc; i++)
	{
		auto const arg = argv[i];
		//std::cout << arg << std::endl; //this is for debugging
		if(strcmp(arg, "--slim") == 0)
		{
			using_slim_version = true;
			std::cout << " *** \t Running with Slimmed Selection \t *** " << std::endl;
		}
		if(strcmp(arg, "-c") == 0)
		{
			using_default_config = false;
			input_config_file_name = argv[i+1];
			//break; //is the break necessary?
		}
		if(strcmp(arg, "-f") == 0)
		{
			using_dynamic_file_prefix = true;
			file_locate_prefix = argv[i+1];
		}
		//checking file paths
		if(strcmp(arg, "--mc") == 0)
		{
			monte_carlo_file = true;
			monte_carlo_file_path = argv[i+1];
		}
		if(strcmp(arg, "--cosmic") == 0)
		{
			cosmic_file = true;
			cosmic_file_path = argv[i+1];
		}
		if(strcmp(arg, "--data") == 0)
		{
			data_file = true;
			data_file_path = argv[i+1];
		}
		if(strcmp(arg, "--dirt") == 0)
		{
			dirt_file = true;
			dirt_file_path = argv[i+1];
		}
		if(strcmp(arg, "--var") == 0)
		{
			variation_file = true;
			variation_file_path = argv[i+1];
			std::string variation_type = variation_file_path;
			//std::string delimiter_1 = "var_";
			//std::string delimiter_2 = ".root";
			//std::string token = variation_type.substr(0, variation_type.find(delimiter_2));
			//file_locate_prefix_2 = token;
		}
		if(strcmp(arg, "--var_mode") == 0)
		{
			variation_mode = true;
		}
		if(strcmp(arg, "--var_mode_bkg") == 0)
		{
			variation_mode_bkg = true;
		}
	}
	if(argc < 2 )  { std::cout << " \n Please inclue the input file path \n " << std::endl; exit(1); }

	//variation mode
	if(variation_mode == true)
	{
		variation_output _var_mode_instance;
		// ./main.exe --var_mode file_name <"same">/""
		_var_mode_instance.run_var(argv[2], argv[3]);
		return 0;
	}
	

	std::vector<double> config;
	std::vector<double> input_config;

	xsecSelection::selection _selection_instance;
	xsecSelection::selection_slim _selection_slim_instance;

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
	        ratio_tolerance,
	        detector_variations
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

	if(using_default_config == true)
	{
		config = default_config;
		std::cout << "(Default) Parameter Config from header" << std::endl;
	}
	if(using_default_config == false)
	{
		config = input_config;
		std::cout << "Paramerter Config from input file" << std::endl;
	}

	std::vector<std::tuple<double, double, std::string> > * results_v = new std::vector<std::tuple<double, double, std::string> >;

	// double _argc = argc;
	// if(using_slim_version == true)    {_argc = _argc - 1; }//this is to account for the "--slim"
	// if(using_default_config == false) {_argc = _argc - 2; }//this is to account for the "-c" and "config_file"
	// if(using_dynamic_file_prefix == true) { _argc = _argc - 2; } // this is to account for the "-f" and "path"


	//variation mode for backgrounds
	if(variation_mode_bkg == true)
	{
		variation_output_bkg _var_mode_instance_bkg;
		// ./main.exe --var_mode_bkg file_name <"same">/""
		_var_mode_instance_bkg.run_var(argv[2], argv[3], config);
		return 0;
	}


	//default state
	if(using_slim_version == false)
	{
		if(variation_file == false)
		{
			_selection_instance.xsecSelection::selection::make_selection(monte_carlo_file_path, cosmic_file_path, data_file_path, dirt_file_path,
			                                                             variation_file_path, config, results_v, file_locate_prefix);
		}
		if(variation_file == true)
		{
			_selection_instance.xsecSelection::selection::make_selection(monte_carlo_file_path, cosmic_file_path, data_file_path, dirt_file_path,
			                                                             variation_file_path, config, results_v, file_locate_prefix);
		}
	}
	//slim selection
	if(using_slim_version == true)
	{
		_selection_slim_instance.xsecSelection::selection_slim::make_selection_slim(monte_carlo_file_path, cosmic_file_path, data_file_path, dirt_file_path,
		                                                                            variation_file_path, config, results_v);
	}
	//write results from selection to output file
	//python script does most of the file managing, since script needs to be contained
	std::ofstream output_file;
	const char * output_file_name = "selection_output.txt"; //I'll expand this to also come from the python script
	//output_file.open(output_file_name, std::ios_base::app); //Only append to this file
	//we want this to be only the current run, for now
	output_file.open(output_file_name, std::ios_base::trunc);
	if(!output_file.is_open())
	{
		std::cout << "*** \t Output File did not open! \t ***" << std::endl;
	}
	if(output_file.is_open())
	{
		std::cout << "*** \t Output File Opened \t *** " << std::endl;

		for(auto const results : * results_v)
		{
			output_file << std::get<2>(results) << ", "
			            << std::get<0>(results) << ", "
			            << std::get<1>(results) << "\n";
		}
		input_config_file.close();
		std::cout << "*** \t Output File Closed \t *** " << std::endl;
	}

	std::cout << "*** \t Exiting C++ Code... \t *** " << std::endl;
	return 0;
}
