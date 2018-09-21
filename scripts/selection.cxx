#include "selection.h"

namespace xsecSelection {
void selection::make_selection( const char * _file1,
                                const char * _file2,
                                const char * _file3,
                                const std::vector<double> _config,
                                std::vector<std::tuple<double, double, std::string> > * results_v,
                                const char * file_locate_prefix
                                )
{

	std::cout << "File Path: " << _file1 << std::endl;
	const bool _verbose = false;
	const bool _post_cuts_verbose = false;
	gErrorIgnoreLevel = kWarning;

	//first we need to open the root file
	TFile * f = new TFile(_file1);
	if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; exit(1); }
	TTree * mytree = (TTree*)f->Get("AnalyzeTPCO/tree");
	TTree * optree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");
	TTree * mctree = (TTree*)f->Get("AnalyzeTPCO/mcparticle_tree");
	TTree * mctruth_counter_tree = (TTree*)f->Get("AnalyzeTPCO/mctruth_counter_tree");

	std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
	mytree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

	selection_functions _functions_instance;
	selection_functions_data _data_functions_instance;
	selection_cuts _cuts_instance;

	std::cout << "Running With: " << POT << " POT " << std::endl;
	const double flux_nue = POT * scaling_nue;
	const double flux_nue_bar = POT * scaling_nue_bar;
	const double flux = flux_nue + flux_nue_bar;


	std::vector<double> selected_energy_vector;

	mctruth_counter_tree->SetBranchAddress("mc_nue_cc_counter",      &mc_nue_cc_counter);
	mctruth_counter_tree->SetBranchAddress("mc_nue_nc_counter",      &mc_nue_nc_counter);
	mctruth_counter_tree->SetBranchAddress("mc_numu_cc_counter",     &mc_numu_cc_counter);
	mctruth_counter_tree->SetBranchAddress("mc_numu_nc_counter",     &mc_numu_nc_counter);
	mctruth_counter_tree->SetBranchAddress("mc_nue_cc_counter_bar",  &mc_nue_cc_counter_bar);
	mctruth_counter_tree->SetBranchAddress("mc_numu_cc_counter_bar", &mc_numu_cc_counter_bar);
	mctruth_counter_tree->SetBranchAddress("mc_nue_nc_counter_bar",  &mc_nue_nc_counter_bar);
	mctruth_counter_tree->SetBranchAddress("mc_numu_nc_counter_bar", &mc_numu_nc_counter_bar);
	mctruth_counter_tree->SetBranchAddress("fMCNuEnergy", &mc_nu_energy);
	mctruth_counter_tree->SetBranchAddress("fMCNuMomentum", &mc_nu_momentum);
	mctruth_counter_tree->SetBranchAddress("fMCNuID", &mc_nu_id);
	mctruth_counter_tree->SetBranchAddress("fMCNuVtxX", &mc_nu_vtx_x);
	mctruth_counter_tree->SetBranchAddress("fMCNuVtxY", &mc_nu_vtx_y);
	mctruth_counter_tree->SetBranchAddress("fMCNuVtxZ", &mc_nu_vtx_z);
	mctruth_counter_tree->SetBranchAddress("fMCNuDirX", &mc_nu_dir_x);
	mctruth_counter_tree->SetBranchAddress("fMCNuDirY", &mc_nu_dir_y);
	mctruth_counter_tree->SetBranchAddress("fMCNuDirZ", &mc_nu_dir_z);
	mctruth_counter_tree->SetBranchAddress("fMCNumParticles", &mc_nu_num_particles);
	mctruth_counter_tree->SetBranchAddress("fMCNumChargedParticles", &mc_nu_num_charged_particles);
	mctruth_counter_tree->SetBranchAddress("fMCEleDirX", &mc_ele_dir_x);
	mctruth_counter_tree->SetBranchAddress("fMCEleDirY", &mc_ele_dir_y);
	mctruth_counter_tree->SetBranchAddress("fMCEleDirZ", &mc_ele_dir_z);
	mctruth_counter_tree->SetBranchAddress("fMCEleEnergy", &mc_ele_energy);
	mctruth_counter_tree->SetBranchAddress("fMCEleMomentum", &mc_ele_momentum);
	mctruth_counter_tree->SetBranchAddress("has_pi0", &has_pi0);
	mctruth_counter_tree->SetBranchAddress("fMCNuTime", &mc_nu_time);

	//configure the externally configurable cut parameters
	std::cout << "\n --- Configuring Parameters --- \n" << std::endl;
	_x1 = _config[0];
	_x2 = _config[1];
	_y1 = _config[2];
	_y2 = _config[3];
	_z1 = _config[4];
	_z2 = _config[5];
	flash_pe_threshold = _config[6];
	flash_time_start = _config[7];
	flash_time_end = _config[8];
	tolerance = _config[9];
	shwr_nue_tolerance = _config[10];
	trk_nue_tolerance = _config[11];
	shwr_hit_threshold = _config[12];
	shwr_hit_threshold_collection = _config[13];
	tolerance_open_angle_min = _config[14];
	tolerance_open_angle_max = _config[15];
	tolerance_dedx_min = _config[16];
	tolerance_dedx_max = _config[17];
	dist_tolerance = _config[18];
	pfp_hits_length_tolerance = _config[19];
	ratio_tolerance = _config[20];
	const std::vector<double> tolerance_open_angle {tolerance_open_angle_min, tolerance_open_angle_max};


	std::vector<double> fv_boundary_v = {_x1, _x2, _y1, _y2, _z1, _z2};

	const int total_mc_entries = mctruth_counter_tree->GetEntries();
	std::cout << "Total MC Entries: " << total_mc_entries << std::endl;

	int _mc_nue_cc_counter = 0;
	int _mc_nue_cc_counter_bar = 0;
	int _mc_numu_cc_counter = 0;
	int _mc_numu_cc_counter_bar = 0;
	int _mc_nue_nc_counter = 0;
	int _mc_nue_nc_counter_bar = 0;
	int _mc_numu_nc_counter = 0;
	int _mc_numu_nc_counter_bar = 0;
	std::vector<bool> true_in_tpc_v;
	true_in_tpc_v.resize(total_mc_entries, false);
	int total_mc_entries_inFV_nue = 0;
	int total_mc_entries_inFV_nue_bar = 0;
	int total_mc_entries_inFV_numu_cc = 0;
	int total_mc_entries_inFV_nue_nc = 0;
	int total_mc_entries_inFV_numu_nc = 0;
	int total_mc_entries_inFV_numu_cc_bar = 0;
	int total_mc_entries_inFV_nue_nc_bar = 0;
	int total_mc_entries_inFV_numu_nc_bar = 0;
	for(int i = 0; i < total_mc_entries; i++)
	{
		mctruth_counter_tree->GetEntry(i);
		if(mc_nu_id == 1) {_mc_nue_cc_counter++; }
		if(mc_nu_id == 2) {_mc_numu_cc_counter++; }
		if(mc_nu_id == 3) {_mc_nue_nc_counter++; }
		if(mc_nu_id == 4) {_mc_numu_nc_counter++; }
		if(mc_nu_id == 5) {_mc_nue_cc_counter_bar++; }
		if(mc_nu_id == 6) {_mc_numu_cc_counter_bar++; }
		if(mc_nu_id == 7) {_mc_nue_nc_counter_bar++; }
		if(mc_nu_id == 8) {_mc_numu_nc_counter_bar++; }
		const bool true_in_tpc = _cuts_instance.selection_cuts::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v);
		true_in_tpc_v.at(i) = true_in_tpc;
		if(true_in_tpc == true && (mc_nu_id == 1)) {total_mc_entries_inFV_nue++; }
		if(true_in_tpc == true && (mc_nu_id == 5)) {total_mc_entries_inFV_nue_bar++; }
		if(true_in_tpc == true && (mc_nu_id == 2)) {total_mc_entries_inFV_numu_cc++; }
		if(true_in_tpc == true && (mc_nu_id == 3)) {total_mc_entries_inFV_nue_nc++; }
		if(true_in_tpc == true && (mc_nu_id == 4)) {total_mc_entries_inFV_numu_nc++; }
		if(true_in_tpc == true && (mc_nu_id == 6)) {total_mc_entries_inFV_numu_cc_bar++; }
		if(true_in_tpc == true && (mc_nu_id == 7)) {total_mc_entries_inFV_nue_nc_bar++; }
		if(true_in_tpc == true && (mc_nu_id == 8)) {total_mc_entries_inFV_numu_nc_bar++; }
	}
	int total_mc_entries_inFV = total_mc_entries_inFV_nue + total_mc_entries_inFV_nue_bar;

	std::cout << "MC Nue CC Counter      --- " << _mc_nue_cc_counter << std::endl;
	std::cout << "MC Nue NC Counter      --- " << _mc_nue_nc_counter << std::endl;
	std::cout << "MC Numu CC Counter     --- " << _mc_numu_cc_counter << std::endl;
	std::cout << "MC Numu NC Counter     --- " << _mc_numu_nc_counter << std::endl;
	std::cout << "MC Nue CC Counter Bar  --- " << _mc_nue_cc_counter_bar << std::endl;
	std::cout << "MC Nue NC Counter Bar  --- " << _mc_nue_nc_counter_bar << std::endl;
	std::cout << "MC Numu CC Counter Bar --- " << _mc_numu_cc_counter_bar << std::endl;
	std::cout << "MC Numu NC Counter Bar --- " << _mc_numu_nc_counter_bar << std::endl;

	double xyz_near_mc = 0;
	double xyz_near_ext = 0;
	double xyz_near_data = 0;

	double xyz_far_mc = 0;
	double xyz_far_ext = 0;
	double xyz_far_data = 0;

	//*****************
	//***** DATA *****
	//*****************

	std::vector<int> * data_in_time_counter_v = new std::vector<int>;
	data_in_time_counter_v->resize(24, 0);
	std::vector<int> * data_pe_counter_v = new std::vector<int>;
	data_pe_counter_v->resize(24, 0);
	std::vector<int> * data_reco_nue_counter_v = new std::vector<int>;
	data_reco_nue_counter_v->resize(24, 0);
	std::vector<int> * data_in_fv_counter_v = new std::vector<int>;
	data_in_fv_counter_v->resize(24, 0);
	std::vector<int> * data_vtx_flash_counter_v = new std::vector<int>;
	data_vtx_flash_counter_v->resize(24, 0);
	std::vector<int> * data_shwr_tpco_counter_v = new std::vector<int>;
	data_shwr_tpco_counter_v->resize(24, 0);
	std::vector<int> * data_trk_tpco_counter_v = new std::vector<int>;
	data_trk_tpco_counter_v->resize(24, 0);
	std::vector<int> * data_hit_threshold_counter_v = new std::vector<int>;
	data_hit_threshold_counter_v->resize(24, 0);
	std::vector<int> * data_open_angle_counter_v = new std::vector<int>;
	data_open_angle_counter_v->resize(24, 0);
	std::vector<int> * data_dedx_counter_v = new std::vector<int>;
	data_dedx_counter_v->resize(24, 0);
	std::vector<int> * data_secondary_shower_counter_v = new std::vector<int>;
	data_secondary_shower_counter_v->resize(24, 0);
	std::vector<int> * data_hit_lengthRatio_counter_v = new std::vector<int>;
	data_hit_lengthRatio_counter_v->resize(24, 0);
	std::vector<int> * data_hit_threshold_collection_counter_v = new std::vector<int>;
	data_hit_threshold_collection_counter_v->resize(24, 0);
	std::vector<int> * data_trk_len_shwr_len_ratio_counter_v = new std::vector<int>;
	data_trk_len_shwr_len_ratio_counter_v->resize(24, 0);
	std::vector<int> * data_track_containment_counter_v = new std::vector<int>;
	data_track_containment_counter_v->resize(24, 0);

	std::vector<std::pair<double, int> > * data_flash_time = new std::vector<std::pair<double, int> >;

	std::vector<TH1 * > * h_ele_pfp_xyz_data = new std::vector<TH1 * >;
	h_ele_pfp_xyz_data->push_back(h_ele_pfp_x_data);
	h_ele_pfp_xyz_data->push_back(h_ele_pfp_y_data);
	h_ele_pfp_xyz_data->push_back(h_ele_pfp_z_data);

	std::vector<TH1 * > * h_any_pfp_xyz_data = new std::vector<TH1 * >;
	h_any_pfp_xyz_data->push_back(h_any_pfp_x_data);
	h_any_pfp_xyz_data->push_back(h_any_pfp_y_data);
	h_any_pfp_xyz_data->push_back(h_any_pfp_z_data);

	std::vector<TH1 * > * h_any_pfp_xyz_last_data = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_data->push_back(h_any_pfp_x_last_data);
	h_any_pfp_xyz_last_data->push_back(h_any_pfp_y_last_data);
	h_any_pfp_xyz_last_data->push_back(h_any_pfp_z_last_data);

	std::vector<int> * tabulated_origins_data = new std::vector<int>;
	tabulated_origins_data->resize(24, 0);

	std::vector<std::tuple<int, int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v_data
	        = new std::vector<std::tuple<int, int, int, double, double, double, std::string, std::string, int, int, double> >;

	std::vector<std::tuple<int, int, int> > duplicate_test;

	//check if a 3rd input parameter was given - if not skip the loop for data
	if(strcmp(_file3, "empty") != 0)
	{
		std::cout << "File Path: " << _file3 << std::endl;
		TFile * data_f = new TFile(_file3);
		if(!data_f->IsOpen()) {std::cout << "Could not open file!" << std::endl; exit(1); }
		TTree * data_tree   = (TTree*)data_f->Get("AnalyzeTPCO/tree");
		TTree * data_optree = (TTree*)data_f->Get("AnalyzeTPCO/optical_tree");

		std::vector<xsecAna::TPCObjectContainer> * data_tpc_object_container_v = nullptr;
		data_tree->SetBranchAddress("TpcObjectContainerV", &data_tpc_object_container_v);

		//need a text file with the run, subrun output
		std::ofstream run_subrun_file;
		run_subrun_file.open("run_subrun_list_data.txt");
		int data_run = 0;
		int data_subrun = 0;
		int last_data_run = 0;
		int last_data_subrun = 0;

		std::cout << "=====================" << std::endl;
		std::cout << "======== Data =======" << std::endl;
		std::cout << "=====================" << std::endl;

		const int data_total_entries = data_tree->GetEntries();
		std::cout << "Total Events: " << data_total_entries << std::endl;

		std::vector<int> * data_passed_runs = new std::vector<int>;
		//passed runs is filled with 0, 1, or 2
		//0 = not in time
		//1 = passed - in time and PE threshold
		// 2 = in-time, but not enough PE -- this counts against my efficiency
		data_passed_runs->resize(data_total_entries);

		_cuts_instance.selection_cuts::loop_flashes(data_f, data_optree, flash_pe_threshold, flash_time_start, flash_time_end,
		                                            data_passed_runs, data_flash_time, 0);
		for(auto const run : * data_passed_runs)
		{
			if(run == 1) {run_sum++; }
			if(run == 0) {out_of_time_sum++; }
			if(run == 2) {low_pe_sum++; }
		}
		std::cout << " -------------------------------------- " << std::endl;
		std::cout << "Passed Runs Vector Size: " << data_passed_runs->size() << std::endl;
		std::cout << "Number Events In-Time & > 50 PE: " << run_sum << std::endl;
		std::cout << "Number Events Not In-Time      : " << out_of_time_sum << std::endl;
		std::cout << "Number Events In-Time & < 50 PE: " << low_pe_sum << std::endl;
		std::cout << " -------------------------------------- " << std::endl;

		//get vector with largest flashes y,z positions
		std::vector< std::vector< double> > * data_largest_flash_v_v = new std::vector < std::vector < double > >;
		_cuts_instance.selection_cuts::SetXYflashVector(data_f, data_optree, data_largest_flash_v_v, flash_time_start, flash_time_end, flash_pe_threshold);
		std::cout << "[Data] Largest Flash Vector Size: " << data_largest_flash_v_v->size() << std::endl;

		for(int event = 0; event < data_total_entries; event++)
		{
			if(_verbose)
			{
				std::cout << "----------------------" << std::endl;
				std::cout << "[DATA EVENT NUMBER] \t " << event << std::endl;
				std::cout << "----------------------" << std::endl;
			}
			data_tree->GetEntry(event);
			//***********************************************************
			//this is where the in-time optical cut actually takes effect
			//***********************************************************
			if(data_passed_runs->at(event) == 0)
			{
				if(_verbose) std::cout << "[Failed In-Time Cut]" << std::endl;
				continue;
			}//false

			//writing the run and subrun values to a text file -
			//this can be used as a cross-check for POT counting
			for(auto const tpc_obj : * data_tpc_object_container_v)
			{
				data_run    = tpc_obj.RunNumber();
				data_subrun = tpc_obj.SubRunNumber();
				if(data_run != last_data_run && data_subrun != last_data_subrun)
				{
					run_subrun_file << data_run << " " << data_subrun << "\n";
					break;
				}
			}
			last_data_run = data_run;
			last_data_subrun = data_subrun;

			for(auto const tpc_obj : * data_tpc_object_container_v)
			{
				duplicate_test.push_back(std::make_tuple(tpc_obj.RunNumber(), tpc_obj.SubRunNumber(), tpc_obj.EventNumber()));
				break;
			}

			//YZ Position of largest flash
			std::vector < double > largest_flash_v = data_largest_flash_v_v->at(event);
			//control for poorly reco flashes
			if(largest_flash_v.at(1) != 0) {h_flash_z_data->Fill(largest_flash_v.at(1)); }

			//List of TPC Objects which pass the cuts
			std::vector<std::pair<int, std::string> > * passed_tpco_data = new std::vector<std::pair<int, std::string> >;
			passed_tpco_data->resize(data_tpc_object_container_v->size());
			std::vector<std::pair<int, std::string> > * dummy_passed_tpco_data = new std::vector<std::pair<int, std::string> >;
			dummy_passed_tpco_data->resize(data_tpc_object_container_v->size());
			//set initial state of objects
			for(int i = 0; i < passed_tpco_data->size(); i++)
			{
				passed_tpco_data->at(i).first = 1;
				passed_tpco_data->at(i).second = "Passed";
				dummy_passed_tpco_data->at(i).first = 1;
				dummy_passed_tpco_data->at(i).second = "Passed";
			}
			//***********************************************************
			//this is where the in-time optical cut again takes effect
			//***********************************************************
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_in_time_counter_v);

			_data_functions_instance.selection_functions_data::XYZPositionData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                   h_any_pfp_xyz_data, h_pfp_zy_vtx_data, xyz_near_data, xyz_far_data);

			//PE threshold cut
			if(data_passed_runs->at(event) == 2)
			{
				if(_verbose) std::cout << "[Passed In-Time Cut] [Failed PE Threshold] " << std::endl;
				continue;
			}

			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_pe_counter_v);

			//****************************
			// ****** reco nue cut *******
			//****************************
			_cuts_instance.selection_cuts::HasNue(data_tpc_object_container_v, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_reco_nue_counter_v);

			//** Testing leading shower length vs hits **//
			_data_functions_instance.selection_functions_data::ShowerLengthvsHitsData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_shwr_len_hits_data);
			_data_functions_instance.selection_functions_data::XYZPositionData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_xyz_data);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_nue_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_nue_cut_data,
			                                                                         h_multiplicity_track_nue_cut_data);

			// _data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, dummy_passed_tpco_data,
			//                                                                                _verbose, h_ele_momentum_no_cut_data);

			//************************
			//******** in fv cut *****
			//************************
			_cuts_instance.selection_cuts::fiducial_volume_cut(data_tpc_object_container_v, fv_boundary_v, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_in_fv_counter_v);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_fv_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_fv_cut_data,
			                                                                         h_multiplicity_track_fv_cut_data);

			//*****************************
			//**** vertex to flash cut ****
			//*****************************
			_data_functions_instance.selection_functions_data::PostCutsVtxFlashData(largest_flash_v, data_tpc_object_container_v, passed_tpco_data,
			                                                                        h_vtx_flash_data);
			_data_functions_instance.selection_functions_data::PostCutsVtxFlashUpstreamData(largest_flash_v, data_tpc_object_container_v, passed_tpco_data,
			                                                                                h_vtx_flash_upstream_data);
			_data_functions_instance.selection_functions_data::PostCutsVtxFlashDownstreamData(largest_flash_v, data_tpc_object_container_v, passed_tpco_data,
			                                                                                  h_vtx_flash_downstream_data);

			_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, data_tpc_object_container_v, tolerance, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_vtx_flash_counter_v);

			_data_functions_instance.selection_functions_data::PostCutsVtxFlashData(largest_flash_v, data_tpc_object_container_v, passed_tpco_data,
			                                                                        h_vtx_flash_data_after);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_flash_vtx_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_flash_vtx_cut_data,
			                                                                         h_multiplicity_track_flash_vtx_cut_data);

			//******************************************************
			//*** distance between pfp shower and nue object cut ***
			//******************************************************
			_data_functions_instance.selection_functions_data::PostCutsShwrVtxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_shwr_vtx_dist_data);

			_cuts_instance.selection_cuts::VtxNuDistance(data_tpc_object_container_v, shwr_nue_tolerance, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_shwr_tpco_counter_v);

			_data_functions_instance.selection_functions_data::PostCutsShwrVtxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_shwr_vtx_dist_data_after);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_shwr_vtx_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_shwr_vtx_cut_data,
			                                                                         h_multiplicity_track_shwr_vtx_cut_data);

			//******************************************************
			// **** distance between pfp track and nue object cut **
			//******************************************************
			_data_functions_instance.selection_functions_data::PostCutTrkVtxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_trk_vtx_dist_data);
			_cuts_instance.selection_cuts::VtxTrackNuDistance(data_tpc_object_container_v, trk_nue_tolerance, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_trk_tpco_counter_v);
			_data_functions_instance.selection_functions_data::PostCutTrkVtxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_trk_vtx_dist_data_after);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_trk_vtx_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_trk_vtx_cut_data,
			                                                                         h_multiplicity_track_trk_vtx_cut_data);

			//****************************************************
			// ******** hit threshold for showers cut *************
			//******************************************************
			_data_functions_instance.selection_functions_data::HitsPlots1DData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                   h_pre_cut_collection_hits_track_data,
			                                                                   h_pre_cut_collection_hits_shower_data,
			                                                                   h_pre_cut_collection_hits_leading_shower_data,
			                                                                   h_pre_cut_total_hits_leading_shower_data);

			_cuts_instance.selection_cuts::HitThreshold(data_tpc_object_container_v, shwr_hit_threshold, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_hit_threshold_counter_v);

			_data_functions_instance.selection_functions_data::dEdxVsOpenAngleData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_dedx_open_angle_data);
			_data_functions_instance.selection_functions_data::LeadingCosThetaData(data_tpc_object_container_v, passed_tpco_data, 0, 0, _verbose, h_ele_cos_theta_data);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_hit_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_hit_cut_data,
			                                                                         h_multiplicity_track_hit_cut_data);

			//***************************************//
			//*** Collection Plane Hits Threshold ***//
			//***************************************//
			_data_functions_instance.selection_functions_data::LeadingPhiData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                  h_ele_pfp_phi_data);
			_data_functions_instance.selection_functions_data::LeadingThetaData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                    h_ele_pfp_theta_data);

			_data_functions_instance.selection_functions_data::HitsPlots1DData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                   h_collection_hits_track_data,
			                                                                   h_collection_hits_shower_data,
			                                                                   h_collection_hits_leading_shower_data,
			                                                                   h_total_hits_leading_shower_data);

			_cuts_instance.selection_cuts::HitThresholdCollection(data_tpc_object_container_v, shwr_hit_threshold_collection, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_hit_threshold_collection_counter_v);

			_data_functions_instance.selection_functions_data::LeadingPhiData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_phi_after_data);
			_data_functions_instance.selection_functions_data::LeadingThetaData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_theta_after_data);

			_data_functions_instance.selection_functions_data::HitsPlots1DData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                   h_collection_hits_track_data_after,
			                                                                   h_collection_hits_shower_data_after,
			                                                                   h_collection_hits_leading_shower_data_after,
			                                                                   h_total_hits_leading_shower_data_after);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_yhit_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_yhit_cut_data,
			                                                                         h_multiplicity_track_yhit_cut_data);


			//*****************************************************
			//****** open angle cut for the leading shower ********
			//******************************************************
			_data_functions_instance.selection_functions_data::NumShowersOpenAngleData(data_tpc_object_container_v, passed_tpco_data, h_pfp_shower_open_angle_data);
			_data_functions_instance.selection_functions_data::PostCutOpenAngleData(data_tpc_object_container_v, passed_tpco_data,
			                                                                        _verbose, h_leading_shower_open_angle_data);
			_data_functions_instance.selection_functions_data::PostCutOpenAngle1ShowerData(tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_leading_shower_open_angle_1_data);
			_data_functions_instance.selection_functions_data::PostCutOpenAngle2PlusShowerData(data_tpc_object_container_v, passed_tpco_data,
			                                                                                   _verbose, h_leading_shower_open_angle_2plus_data);
			_cuts_instance.selection_cuts::OpenAngleCut(data_tpc_object_container_v, passed_tpco_data, tolerance_open_angle, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_open_angle_counter_v);

			_data_functions_instance.selection_functions_data::PostCutOpenAngleData(data_tpc_object_container_v, passed_tpco_data,
			                                                                        _verbose, h_leading_shower_open_angle_data_after);
			_data_functions_instance.selection_functions_data::NumShowersOpenAngleData(data_tpc_object_container_v, passed_tpco_data, h_pfp_shower_dedx_data);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_open_angle_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_open_angle_cut_data,
			                                                                         h_multiplicity_track_open_angle_cut_data);

			//*****************************************************
			//*********** dEdx cut for the leading shower *********
			//******************************************************
			_data_functions_instance.selection_functions_data::PostCutsdEdxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_dedx_cuts_data);

			_data_functions_instance.selection_functions_data::PostCutsdEdxAltScaleData(data_tpc_object_container_v, passed_tpco_data,
			                                                                            _verbose, 0.98, h_dedx_cuts_scale_1_data);
			_data_functions_instance.selection_functions_data::PostCutsdEdxAltScaleData(data_tpc_object_container_v, passed_tpco_data,
			                                                                            _verbose, 0.97, h_dedx_cuts_scale_2_data);
			_data_functions_instance.selection_functions_data::PostCutsdEdxAltScaleData(data_tpc_object_container_v, passed_tpco_data,
			                                                                            _verbose, 0.95, h_dedx_cuts_scale_3_data);

			_data_functions_instance.selection_functions_data::dEdxThetaData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_dedx_theta_pre_cuts_data);
			_data_functions_instance.selection_functions_data::dedxThetaSliceData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                      h_dedx_slice_1_data,
			                                                                      h_dedx_slice_2_data,
			                                                                      h_dedx_slice_3_data);

			_data_functions_instance.selection_functions_data::dedxThetaSliceData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                      h_dedx_slice_1_zoom_data,
			                                                                      h_dedx_slice_2_zoom_data,
			                                                                      h_dedx_slice_3_zoom_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_pre_dedx_data, h_multiplicity_track_pre_dedx_data);

			_cuts_instance.selection_cuts::dEdxCut(data_tpc_object_container_v, passed_tpco_data, tolerance_dedx_min, tolerance_dedx_max, _verbose, false);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_dedx_counter_v);

			_data_functions_instance.selection_functions_data::PostCutsdEdxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_dedx_cuts_data_after);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_dedx_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_dedx_cut_data,
			                                                                         h_multiplicity_track_dedx_cut_data);

			//***************************************************************************
			// ******* Secondary Showers Distance Cut *****************
			//***************************************************************************
			_data_functions_instance.selection_functions_data::SecondaryShowersDistData(data_tpc_object_container_v, passed_tpco_data,
			                                                                            _verbose, h_second_shwr_dist_data);
			_cuts_instance.selection_cuts::SecondaryShowersDistCut(data_tpc_object_container_v, passed_tpco_data, _verbose, dist_tolerance);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_secondary_shower_counter_v);
			_data_functions_instance.selection_functions_data::SecondaryShowersDistData(data_tpc_object_container_v, passed_tpco_data,
			                                                                            _verbose, h_second_shwr_dist_data_after);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_2shwr_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_2shwr_cut_data,
			                                                                         h_multiplicity_track_2shwr_cut_data);

			//******************************************************************************
			// ********** Hit Length Ratio Cut *************
			//******************************************************************************
			_data_functions_instance.selection_functions_data::HitLengthRatioData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_hit_length_ratio_data);
			_cuts_instance.selection_cuts::HitLengthRatioCut(data_tpc_object_container_v, passed_tpco_data, _verbose, pfp_hits_length_tolerance);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_hit_lengthRatio_counter_v);

			_data_functions_instance.selection_functions_data::HitLengthRatioData(data_tpc_object_container_v, passed_tpco_data,
			                                                                      _verbose, h_hit_length_ratio_data_after);
			_data_functions_instance.selection_functions_data::PlaneHitsComparisonTrackData(data_tpc_object_container_v, passed_tpco_data,
			                                                                                _verbose, h_collection_total_hits_track_data);
			_data_functions_instance.selection_functions_data::PlaneHitsComparisonShowerData(data_tpc_object_container_v, passed_tpco_data,
			                                                                                 _verbose, h_collection_total_hits_shower_data);
			_data_functions_instance.selection_functions_data::PlaneHitsComparisonLeadingShowerData(data_tpc_object_container_v, passed_tpco_data,
			                                                                                        _verbose, h_collection_total_hits_leading_shower_data);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_hit_length_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_hit_length_cut_data,
			                                                                         h_multiplicity_track_hit_length_cut_data);

			//******************************************************************************
			//*** cut for longest track / leading shower ratio *** //
			//******************************************************************************
			_data_functions_instance.selection_functions_data::LeadingShowerLengthData(data_tpc_object_container_v, passed_tpco_data,
			                                                                           _verbose, h_leading_shwr_length_data);
			_data_functions_instance.selection_functions_data::LeadingShowerTrackLengthsData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                                 h_leading_shwr_trk_length_data);

			_cuts_instance.selection_cuts::LongestTrackLeadingShowerCut(data_tpc_object_container_v, passed_tpco_data, _verbose, ratio_tolerance);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_trk_len_shwr_len_ratio_counter_v);

			_data_functions_instance.selection_functions_data::LeadingShowerLengthData(data_tpc_object_container_v, passed_tpco_data,
			                                                                           _verbose, h_leading_shwr_length_data_after);
			_data_functions_instance.selection_functions_data::LeadingShowerTrackLengthsData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                                 h_leading_shwr_trk_length_data_after);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_length_ratio_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_length_ratio_cut_data,
			                                                                         h_multiplicity_track_length_ratio_cut_data);

			//***************************************************************
			//*** contained track cut *** //
			//**************************************************************
			_data_functions_instance.selection_functions_data::IsContainedPlotData(data_tpc_object_container_v, passed_tpco_data,
			                                                                       _verbose, fv_boundary_v, h_track_containment_data);

			_cuts_instance.selection_cuts::ContainedTracksCut(data_tpc_object_container_v, passed_tpco_data, _verbose, fv_boundary_v, true);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_track_containment_counter_v);

			_data_functions_instance.selection_functions_data::PostCutsLeadingMomentumData(data_tpc_object_container_v, passed_tpco_data,
			                                                                               _verbose, h_ele_momentum_containment_cut_data);

			_data_functions_instance.selection_functions_data::EventMultiplicityData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                         h_multiplicity_shower_containment_cut_data,
			                                                                         h_multiplicity_track_containment_cut_data);

			//*********** Data *********
			_functions_instance.selection_functions::FillPostCutVector(data_tpc_object_container_v, passed_tpco_data, post_cuts_v_data);
			//*************************************
			// ******** End Selection Cuts! *******
			//*************************************
			_data_functions_instance.selection_functions_data::LeadingMomentumData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_momentum_data);
			_data_functions_instance.selection_functions_data::LeadingMomentumTrackTopologyData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                                    h_ele_pfp_momentum_no_track_data, h_ele_pfp_momentum_has_track_data);

			_data_functions_instance.selection_functions_data::LeadingPhiTrackTopologyData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                               h_ele_pfp_phi_no_track_data, h_ele_pfp_phi_has_track_data);

			_data_functions_instance.selection_functions_data::LeadingThetaTrackTopologyData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                                 h_ele_pfp_theta_no_track_data, h_ele_pfp_theta_has_track_data);

			_data_functions_instance.selection_functions_data::LeadingPhiData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_phi_last_data);
			_data_functions_instance.selection_functions_data::LeadingThetaData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_theta_last_data);
			_data_functions_instance.selection_functions_data::LeadingCosThetaData(data_tpc_object_container_v, passed_tpco_data, theta_translation, phi_translation,
			                                                                       _verbose, h_ele_cos_theta_last_trans_data);
			_data_functions_instance.selection_functions_data::LeadingCosThetaData(data_tpc_object_container_v, passed_tpco_data, 0, 0,
			                                                                       _verbose, h_ele_cos_theta_last_data);

			_data_functions_instance.selection_functions_data::XYZPositionData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_any_pfp_xyz_last_data);

			_data_functions_instance.selection_functions_data::dEdxThetaData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_dedx_theta_data);

			_data_functions_instance.selection_functions_data::LeadingMomentumThetaSliceData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                                 h_ele_momentum_slice_1_data,
			                                                                                 h_ele_momentum_slice_2_data,
			                                                                                 h_ele_momentum_slice_3_data);

			_data_functions_instance.selection_functions_data::EnergyCosThetaData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_eng_costheta_data);
			_data_functions_instance.selection_functions_data::EnergyCosThetaSlicesData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                            0, 0,
			                                                                            h_ele_eng_for_data,
			                                                                            h_ele_eng_mid_data,
			                                                                            h_ele_eng_back_data);
			_data_functions_instance.selection_functions_data::EnergyCosThetaSlicesData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                            theta_translation, phi_translation,
			                                                                            h_ele_eng_for_trans_data,
			                                                                            h_ele_eng_mid_trans_data,
			                                                                            h_ele_eng_back_trans_data);

			_data_functions_instance.selection_functions_data::LeadingKinematicsShowerTopologyData(data_tpc_object_container_v, passed_tpco_data, _verbose,
			                                                                                       h_ele_pfp_momentum_1shwr_data, h_ele_pfp_momentum_2shwr_data,
			                                                                                       h_ele_pfp_theta_1shwr_data, h_ele_pfp_theta_2shwr_data,
			                                                                                       h_ele_pfp_phi_1shwr_data, h_ele_pfp_phi_2shwr_data);

			//delete at the very end!
			delete passed_tpco_data;
		}
	}//end data running
	 //****************************
	 //*** END Data Calculation ***
	 //****************************

//***********************************
	//***********************************
	//*** In-time Cosmics Calculation ***
	//***********************************
//***********************************

	std::vector<int> * intime_in_time_counter_v = new std::vector<int>;
	intime_in_time_counter_v->resize(24, 0);
	std::vector<int> * intime_pe_counter_v = new std::vector<int>;
	intime_pe_counter_v->resize(24, 0);
	std::vector<int> * intime_reco_nue_counter_v = new std::vector<int>;
	intime_reco_nue_counter_v->resize(24, 0);
	std::vector<int> * intime_in_fv_counter_v = new std::vector<int>;
	intime_in_fv_counter_v->resize(24, 0);
	std::vector<int> * intime_vtx_flash_counter_v = new std::vector<int>;
	intime_vtx_flash_counter_v->resize(24, 0);
	std::vector<int> * intime_shwr_tpco_counter_v = new std::vector<int>;
	intime_shwr_tpco_counter_v->resize(24, 0);
	std::vector<int> * intime_trk_tpco_counter_v = new std::vector<int>;
	intime_trk_tpco_counter_v->resize(24, 0);
	std::vector<int> * intime_hit_threshold_counter_v = new std::vector<int>;
	intime_hit_threshold_counter_v->resize(24, 0);
	std::vector<int> * intime_open_angle_counter_v = new std::vector<int>;
	intime_open_angle_counter_v->resize(24, 0);
	std::vector<int> * intime_dedx_counter_v = new std::vector<int>;
	intime_dedx_counter_v->resize(24, 0);
	std::vector<int> * intime_secondary_shower_counter_v = new std::vector<int>;
	intime_secondary_shower_counter_v->resize(24, 0);
	std::vector<int> * intime_hit_lengthRatio_counter_v = new std::vector<int>;
	intime_hit_lengthRatio_counter_v->resize(24, 0);
	std::vector<int> * intime_hit_threshold_collection_counter_v = new std::vector<int>;
	intime_hit_threshold_collection_counter_v->resize(24, 0);
	std::vector<int> * intime_trk_len_shwr_len_ratio_counter_v = new std::vector<int>;
	intime_trk_len_shwr_len_ratio_counter_v->resize(24, 0);
	std::vector<int> * intime_track_containment_counter_v = new std::vector<int>;
	intime_track_containment_counter_v->resize(24, 0);

	std::vector<std::pair<double, int> > * intime_flash_time = new std::vector<std::pair<double, int> >;

	std::vector<TH1 * > * h_ele_pfp_xyz_intime = new std::vector<TH1 * >;
	h_ele_pfp_xyz_intime->push_back(h_ele_pfp_x_intime);
	h_ele_pfp_xyz_intime->push_back(h_ele_pfp_y_intime);
	h_ele_pfp_xyz_intime->push_back(h_ele_pfp_z_intime);

	std::vector<TH1 * > * h_any_pfp_xyz_intime = new std::vector<TH1 * >;
	h_any_pfp_xyz_intime->push_back(h_any_pfp_x_intime);
	h_any_pfp_xyz_intime->push_back(h_any_pfp_y_intime);
	h_any_pfp_xyz_intime->push_back(h_any_pfp_z_intime);

	std::vector<TH1 * > * h_any_pfp_xyz_last_intime = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_intime->push_back(h_any_pfp_x_last_intime);
	h_any_pfp_xyz_last_intime->push_back(h_any_pfp_y_last_intime);
	h_any_pfp_xyz_last_intime->push_back(h_any_pfp_z_last_intime);

	std::vector<int> * tabulated_origins = new std::vector<int>;
	tabulated_origins->resize(24, 0);
	std::vector<int> * tabulated_origins_intime = new std::vector<int>;
	tabulated_origins_intime->resize(24, 0);

	std::vector<std::tuple<int, int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v
	        = new std::vector<std::tuple<int, int, int, double, double, double, std::string, std::string, int, int, double> >;

	if(strcmp(_file2, "empty") != 0)
	{
		std::cout << "File Path: " << _file2 << std::endl;
		TFile * intime_f = new TFile(_file2);
		if(!intime_f->IsOpen()) {std::cout << "Could not open file!" << std::endl; exit(1); }
		TTree * intime_tree   = (TTree*)intime_f->Get("AnalyzeTPCO/tree");
		TTree * intime_optree = (TTree*)intime_f->Get("AnalyzeTPCO/optical_tree");

		std::vector<xsecAna::TPCObjectContainer> * intime_tpc_object_container_v = new std::vector<xsecAna::TPCObjectContainer>;
		intime_tree->SetBranchAddress("TpcObjectContainerV", &intime_tpc_object_container_v);

		std::ofstream run_subrun_file;
		run_subrun_file.open("run_subrun_list_intime.txt");
		int data_run = 0;
		int data_subrun = 0;
		int last_data_run = 0;
		int last_data_subrun = 0;

		std::cout << "=====================" << std::endl;
		std::cout << "== In-Time Cosmics ==" << std::endl;
		std::cout << "=====================" << std::endl;

		const int in_time_total_entries = intime_tree->GetEntries();
		std::cout << "Total Events: " << in_time_total_entries << std::endl;

		std::vector<int> * intime_passed_runs = new std::vector<int>;
		//passed runs is filled with 0, 1, or 2
		//0 = not in time
		//1 = passed - in time and PE threshold
		// 2 = in-time, but not enough PE -- this counts against my efficiency
		intime_passed_runs->resize(in_time_total_entries);

		_cuts_instance.selection_cuts::loop_flashes(intime_f, intime_optree, flash_pe_threshold,
		                                            flash_time_start, flash_time_end, intime_passed_runs, intime_flash_time, 1);
		for(auto const run : * intime_passed_runs)
		{
			if(run == 1) {run_sum++; }
			if(run == 0) {out_of_time_sum++; }
			if(run == 2) {low_pe_sum++; }
		}
		std::cout << " -------------------------------------- " << std::endl;
		std::cout << "Passed Runs Vector Size: " << intime_passed_runs->size() << std::endl;
		std::cout << "Number Events In-Time & > 50 PE: " << run_sum << std::endl;
		std::cout << "Number Events Not In-Time      : " << out_of_time_sum << std::endl;
		std::cout << "Number Events In-Time & < 50 PE: " << low_pe_sum << std::endl;
		std::cout << " -------------------------------------- " << std::endl;

		//get vector with largest flashes y,z positions
		std::vector< std::vector< double> > * intime_largest_flash_v_v = new std::vector < std::vector < double > >;
		_cuts_instance.selection_cuts::SetXYflashVector(intime_f, intime_optree, intime_largest_flash_v_v, flash_time_start, flash_time_end, flash_pe_threshold);
		std::cout << "Largest Flash Vector Size: " << intime_largest_flash_v_v->size() << std::endl;

		for(int event = 0; event < in_time_total_entries; event++)
		{
			if(_verbose)
			{
				std::cout << "----------------------" << std::endl;
				std::cout << "[IN-TIME EVENT NUMBER] \t " << event << std::endl;
				std::cout << "----------------------" << std::endl;
			}
			intime_tree->GetEntry(event);

			//writing the run and subrun values to a text file -
			//this can be used as a cross-check for POT counting

			for(auto const tpc_obj : * intime_tpc_object_container_v)
			{
				data_run    = tpc_obj.RunNumber();
				data_subrun = tpc_obj.SubRunNumber();
				if(data_run != last_data_run && data_subrun != last_data_subrun)
				{
					run_subrun_file << data_run << " " << data_subrun << "\n";
					break;
				}
			}
			last_data_run = data_run;
			last_data_subrun = data_subrun;

			//***********************************************************
			//this is where the in-time optical cut actually takes effect
			//***********************************************************
			if(_verbose) {std::cout << "In-Time Optical Cut (1)" << std::endl; }
			if(intime_passed_runs->at(event) == 0)
			{
				if(_verbose) std::cout << "[Failed In-Time Cut]" << std::endl;
				continue;
			}//false

			//YZ Position of largest flash
			std::vector < double > largest_flash_v = intime_largest_flash_v_v->at(event);
			//control for poorly reco flashes
			if(largest_flash_v.at(1) != 0) {h_flash_z_intime->Fill(largest_flash_v.at(1)); }

			//List of TPC Objects which pass the cuts
			std::vector<std::pair<int, std::string> > * passed_tpco_intime = new std::vector<std::pair<int, std::string> >;
			passed_tpco_intime->resize(intime_tpc_object_container_v->size());
			std::vector<std::pair<int, std::string> > * dummy_passed_tpco_intime = new std::vector<std::pair<int, std::string> >;
			dummy_passed_tpco_intime->resize(intime_tpc_object_container_v->size());
			//set initial state of objects
			for(int i = 0; i < passed_tpco_intime->size(); i++)
			{
				passed_tpco_intime->at(i).first = 1;
				passed_tpco_intime->at(i).second = "Passed";
				dummy_passed_tpco_intime->at(i).first = 1;
				dummy_passed_tpco_intime->at(i).second = "Passed";
			}
			//***********************************************************
			//this is where the in-time optical cut again takes effect
			//***********************************************************
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_in_time_counter_v);

			_functions_instance.selection_functions::XYZPositionInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                           h_any_pfp_xyz_intime, h_pfp_zy_vtx_ext, xyz_near_ext, xyz_far_ext);

			//PE threshold cut
			if(_verbose) {std::cout << "In-Time Optical Cut (2)" << std::endl; }
			if(intime_passed_runs->at(event) == 2)
			{
				if(_verbose) std::cout << "[Passed In-Time Cut] [Failed PE Threshold] " << std::endl;
				continue;
			}

			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_pe_counter_v);

			//****************************
			// ****** reco nue cut *******
			//****************************
			if(_verbose) {std::cout << "In-Time Reco Nue Cut" << std::endl; }
			_cuts_instance.selection_cuts::HasNue(intime_tpc_object_container_v, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_reco_nue_counter_v);

			//** Testing leading shower length vs hits **//
			_functions_instance.selection_functions::ShowerLengthvsHitsInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_shwr_len_hits_intime);
			_functions_instance.selection_functions::XYZPositionInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_ele_pfp_xyz_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_nue_cut_intime);
			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, dummy_passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_no_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_nue_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_nue_cut_intime,
			                                                                 h_multiplicity_track_nue_cut_intime);

			//************************
			//******** in fv cut *****
			//************************
			if(_verbose) {std::cout << "In-Time FV Cut" << std::endl; }
			_cuts_instance.selection_cuts::fiducial_volume_cut(intime_tpc_object_container_v, fv_boundary_v, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_in_fv_counter_v);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_fv_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_fv_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_fv_cut_intime,
			                                                                 h_multiplicity_track_fv_cut_intime);

			//*****************************
			//**** vertex to flash cut ****
			//*****************************
			if(_verbose) {std::cout << "In-Time Vertex-Flash Cut" << std::endl; }
			_functions_instance.selection_functions::PostCutsVtxFlashInTime(largest_flash_v, intime_tpc_object_container_v, passed_tpco_intime,
			                                                                _verbose, h_vtx_flash_intime);
			_functions_instance.selection_functions::PostCutsVtxFlashUpstreamInTime(largest_flash_v, intime_tpc_object_container_v, passed_tpco_intime,
			                                                                        _verbose, h_vtx_flash_upstream_intime);
			_functions_instance.selection_functions::PostCutsVtxFlashDownstreamInTime(largest_flash_v, intime_tpc_object_container_v, passed_tpco_intime,
			                                                                          _verbose, h_vtx_flash_downstream_intime);

			_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, intime_tpc_object_container_v, tolerance, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_vtx_flash_counter_v);

			_functions_instance.selection_functions::PostCutsVtxFlashInTime(largest_flash_v, intime_tpc_object_container_v, passed_tpco_intime,
			                                                                _verbose, h_vtx_flash_intime_after);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_flash_vtx_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_flash_vtx_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_flash_vtx_cut_intime,
			                                                                 h_multiplicity_track_flash_vtx_cut_intime);
			//******************************************************
			//*** distance between pfp shower and nue object cut ***
			//******************************************************
			if(_verbose) {std::cout << "In-Time Shower-Vtx Cut" << std::endl; }
			_functions_instance.selection_functions::PostCutsShwrVtxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_shwr_vtx_dist_intime);

			_cuts_instance.selection_cuts::VtxNuDistance(intime_tpc_object_container_v, shwr_nue_tolerance, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_shwr_tpco_counter_v);

			_functions_instance.selection_functions::PostCutsShwrVtxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_shwr_vtx_dist_intime_after);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_shwr_vtx_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_shwr_vtx_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_shwr_vtx_cut_intime,
			                                                                 h_multiplicity_track_shwr_vtx_cut_intime);

			//******************************************************
			// **** distance between pfp track and nue object cut **
			//******************************************************
			if(_verbose) {std::cout << "In-Time Track-Vtx Cut" << std::endl; }
			_functions_instance.selection_functions::PostCutTrkVtxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_trk_vtx_dist_intime);

			_cuts_instance.selection_cuts::VtxTrackNuDistance(intime_tpc_object_container_v, trk_nue_tolerance, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_trk_tpco_counter_v);

			_functions_instance.selection_functions::PostCutTrkVtxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_trk_vtx_dist_intime_after);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_trk_vtx_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_trk_vtx_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_trk_vtx_cut_intime,
			                                                                 h_multiplicity_track_trk_vtx_cut_intime);

			//****************************************************
			// ******** hit threshold for showers cut *************
			//******************************************************
			_functions_instance.selection_functions::HitsPlots1DInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                           h_pre_cut_collection_hits_track_intime,
			                                                           h_pre_cut_collection_hits_shower_intime,
			                                                           h_pre_cut_collection_hits_leading_shower_intime,
			                                                           h_pre_cut_total_hits_leading_shower_intime);

			_cuts_instance.selection_cuts::HitThreshold(intime_tpc_object_container_v, shwr_hit_threshold, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_hit_threshold_counter_v);

			_functions_instance.selection_functions::dEdxVsOpenAngleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_dedx_open_angle_intime);
			_functions_instance.selection_functions::LeadingCosThetaInTime(intime_tpc_object_container_v, passed_tpco_intime, 0, 0,
			                                                               _verbose, h_ele_cos_theta_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_hit_cut_intime);


			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_hit_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_hit_cut_intime,
			                                                                 h_multiplicity_track_hit_cut_intime);

			//***************************************//
			//*** Collection Plane Hits Threshold ***//
			//***************************************//
			_functions_instance.selection_functions::LeadingPhiInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                          h_ele_pfp_phi_intime);

			_functions_instance.selection_functions::LeadingThetaInTime(intime_tpc_object_container_v, passed_tpco_intime, 0, 0, _verbose,
			                                                            h_ele_pfp_theta_intime);

			_functions_instance.selection_functions::HitsPlots1DInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                           h_collection_hits_track_intime,
			                                                           h_collection_hits_shower_intime,
			                                                           h_collection_hits_leading_shower_intime,
			                                                           h_total_hits_leading_shower_intime);

			_cuts_instance.selection_cuts::HitThresholdCollection(intime_tpc_object_container_v, shwr_hit_threshold_collection, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_hit_threshold_collection_counter_v);

			_functions_instance.selection_functions::LeadingPhiInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_ele_pfp_phi_after_intime);

			_functions_instance.selection_functions::LeadingThetaInTime(intime_tpc_object_container_v, passed_tpco_intime, 0, 0, _verbose, h_ele_pfp_theta_after_intime);

			_functions_instance.selection_functions::HitsPlots1DInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                           h_collection_hits_track_intime_after,
			                                                           h_collection_hits_shower_intime_after,
			                                                           h_collection_hits_leading_shower_intime_after,
			                                                           h_total_hits_leading_shower_intime_after);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_yhit_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_yhit_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_yhit_cut_intime,
			                                                                 h_multiplicity_track_yhit_cut_intime);

			//*****************************************************
			//****** open angle cut for the leading shower ********
			//******************************************************
			_functions_instance.selection_functions::NumShowersOpenAngleInTime(intime_tpc_object_container_v, passed_tpco_intime, h_pfp_shower_open_angle_intime);
			_functions_instance.selection_functions::PostCutOpenAngleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                _verbose, h_leading_shower_open_angle_intime);
			_functions_instance.selection_functions::PostCutOpenAngle1ShowerInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                       _verbose, h_leading_shower_open_angle_1_intime);
			_functions_instance.selection_functions::PostCutOpenAngle2PlusShowerInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                           _verbose, h_leading_shower_open_angle_2plus_intime);
			_cuts_instance.selection_cuts::OpenAngleCut(intime_tpc_object_container_v, passed_tpco_intime, tolerance_open_angle, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_open_angle_counter_v);

			_functions_instance.selection_functions::NumShowersOpenAngleInTime(intime_tpc_object_container_v, passed_tpco_intime, h_pfp_shower_dedx_intime);

			_functions_instance.selection_functions::PostCutOpenAngleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                _verbose, h_leading_shower_open_angle_intime_after);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_open_angle_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_open_angle_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_open_angle_cut_intime,
			                                                                 h_multiplicity_track_open_angle_cut_intime);

			//*****************************************************
			//*********** dEdx cut for the leading shower *********
			//******************************************************
			_functions_instance.selection_functions::PostCutsdEdxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_dedx_cuts_intime);

			_functions_instance.selection_functions::PostCutsdEdxAltScaleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                    _verbose, 0.98, h_dedx_cuts_scale_1_intime);
			_functions_instance.selection_functions::PostCutsdEdxAltScaleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                    _verbose, 0.97, h_dedx_cuts_scale_2_intime);
			_functions_instance.selection_functions::PostCutsdEdxAltScaleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                    _verbose, 0.95, h_dedx_cuts_scale_3_intime);

			_functions_instance.selection_functions::dEdxCollectionAngleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                   _verbose, h_dedx_collection_angle_intime);
			_functions_instance.selection_functions::dEdxThetaInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                         _verbose, h_dedx_theta_pre_cuts_intime);

			_functions_instance.selection_functions::PostCutsdEdxTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                        _verbose, h_dedx_cuts_ext_unmatched);

			_functions_instance.selection_functions::PostCutsdedxThetaSliceInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                      h_dedx_slice_1_intime,
			                                                                      h_dedx_slice_2_intime,
			                                                                      h_dedx_slice_3_intime);

			_functions_instance.selection_functions::PostCutsdedxThetaSliceInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                      h_dedx_slice_1_zoom_intime,
			                                                                      h_dedx_slice_2_zoom_intime,
			                                                                      h_dedx_slice_3_zoom_intime);

			_cuts_instance.selection_cuts::dEdxCut(intime_tpc_object_container_v, passed_tpco_intime, tolerance_dedx_min, tolerance_dedx_max, _verbose, true);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_dedx_counter_v);

			_functions_instance.selection_functions::PostCutsdEdxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_dedx_cuts_intime_after);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_dedx_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_dedx_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_dedx_cut_intime,
			                                                                 h_multiplicity_track_dedx_cut_intime);

			//***************************************************************************
			// ******* Secondary Showers Distance Cut *****************
			//***************************************************************************
			_functions_instance.selection_functions::SecondaryShowersDistInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                    _verbose, h_second_shwr_dist_intime);
			_cuts_instance.selection_cuts::SecondaryShowersDistCut(intime_tpc_object_container_v, passed_tpco_intime, _verbose, dist_tolerance);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_secondary_shower_counter_v);

			_functions_instance.selection_functions::SecondaryShowersDistInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                    _verbose, h_second_shwr_dist_intime_after);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_2shwr_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_2shwr_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_2shwr_cut_intime,
			                                                                 h_multiplicity_track_2shwr_cut_intime);

			//******************************************************************************
			// ********** Hit Length Ratio Cut *************
			//******************************************************************************
			_functions_instance.selection_functions::HitLengthRatioInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_hit_length_ratio_intime);

			_cuts_instance.selection_cuts::HitLengthRatioCut(intime_tpc_object_container_v, passed_tpco_intime, _verbose, pfp_hits_length_tolerance);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_hit_lengthRatio_counter_v);

			_functions_instance.selection_functions::PlaneHitsComparisonTrackInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                        _verbose, h_collection_total_hits_track_intime);
			_functions_instance.selection_functions::PlaneHitsComparisonShowerInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                         _verbose, h_collection_total_hits_shower_intime);
			_functions_instance.selection_functions::PlaneHitsComparisonLeadingShowerInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                                _verbose, h_collection_total_hits_leading_shower_intime);

			_functions_instance.selection_functions::HitLengthRatioInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                              _verbose, h_hit_length_ratio_intime_after);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_hit_length_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_hit_length_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_hit_length_cut_intime,
			                                                                 h_multiplicity_track_hit_length_cut_intime);

			//******************************************************************************
			//*** cut for longest track / leading shower ratio *** //
			//******************************************************************************
			_functions_instance.selection_functions::LeadingShowerLengthInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                   _verbose, h_leading_shwr_length_intime);

			_functions_instance.selection_functions::LeadingShowerTrackLengthsInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                         h_leading_shwr_trk_length_intime);

			_cuts_instance.selection_cuts::LongestTrackLeadingShowerCut(intime_tpc_object_container_v, passed_tpco_intime, _verbose, ratio_tolerance);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_trk_len_shwr_len_ratio_counter_v);

			_functions_instance.selection_functions::LeadingShowerLengthInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                   _verbose, h_leading_shwr_length_intime_after);

			_functions_instance.selection_functions::LeadingShowerTrackLengthsInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                         h_leading_shwr_trk_length_intime_after);

			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_length_ratio_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_length_ratio_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_length_ratio_cut_intime,
			                                                                 h_multiplicity_track_length_ratio_cut_intime);

			//******************************************************************************
			//*** contained track cut *** //
			//******************************************************************************
			_functions_instance.selection_functions::IsContainedPlotInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                               _verbose, fv_boundary_v, h_track_containment_intime);

			_cuts_instance.selection_cuts::ContainedTracksCut(intime_tpc_object_container_v, passed_tpco_intime, _verbose, fv_boundary_v, true);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_track_containment_counter_v);


			_functions_instance.selection_functions::PostCutsLeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_momentum_containment_cut_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                   h_leading_momentum_containment_cut_ext_unmatched);

			_functions_instance.selection_functions::EventMultiplicityInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                 h_multiplicity_shower_containment_cut_intime,
			                                                                 h_multiplicity_track_containment_cut_intime);

			//*********** In-time Cosmics *********
			//*************************************
			// ******** End Selection Cuts! *******
			//*************************************
			_functions_instance.selection_functions::LeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_ele_pfp_momentum_intime);
			_functions_instance.selection_functions::LeadingMomentumTrackTopologyInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                            h_ele_pfp_momentum_no_track_intime, h_ele_pfp_momentum_has_track_intime);

			_functions_instance.selection_functions::LeadingPhiTrackTopologyInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                       h_ele_pfp_phi_no_track_intime, h_ele_pfp_phi_has_track_intime);

			_functions_instance.selection_functions::LeadingThetaTrackTopologyInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                         h_ele_pfp_theta_no_track_intime, h_ele_pfp_theta_has_track_intime);

			_functions_instance.selection_functions::PostCutsdEdxTrueParticleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                        _verbose, h_dedx_cuts_last_ext_unmatched);

			_functions_instance.selection_functions::LeadingPhiInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_ele_pfp_phi_last_intime);
			_functions_instance.selection_functions::LeadingThetaInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                            theta_translation, phi_translation, _verbose, h_ele_pfp_theta_last_intime);
			_functions_instance.selection_functions::LeadingCosThetaInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                               0, 0, _verbose, h_ele_cos_theta_last_intime);
			_functions_instance.selection_functions::LeadingCosThetaInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                               theta_translation, phi_translation, _verbose, h_ele_cos_theta_last_trans_intime);
			_functions_instance.selection_functions::EnergyCosThetaInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_ele_eng_costheta_intime);
			_functions_instance.selection_functions::EnergyCosThetaSlicesInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                    0, 0,
			                                                                    h_ele_eng_for_intime,
			                                                                    h_ele_eng_mid_intime,
			                                                                    h_ele_eng_back_intime);
			_functions_instance.selection_functions::EnergyCosThetaSlicesInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                    theta_translation, phi_translation,
			                                                                    h_ele_eng_for_trans_intime,
			                                                                    h_ele_eng_mid_trans_intime,
			                                                                    h_ele_eng_back_trans_intime);

			_functions_instance.selection_functions::PostCutsLeadingMomentumThetaSliceInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                                 h_ele_momentum_slice_1_intime,
			                                                                                 h_ele_momentum_slice_2_intime,
			                                                                                 h_ele_momentum_slice_3_intime);

			_functions_instance.selection_functions::dEdxThetaInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                         _verbose, h_dedx_theta_intime);

			_functions_instance.selection_functions::XYZPositionInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_any_pfp_xyz_last_intime);

			_functions_instance.selection_functions::FillPostCutVector(intime_tpc_object_container_v, passed_tpco_intime, post_cuts_v);
			_functions_instance.selection_functions::LeadingThetaPhiInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_ele_theta_phi_intime);

			_functions_instance.selection_functions::LeadingKinematicsShowerTopologyInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose,
			                                                                               h_ele_pfp_momentum_1shwr_intime, h_ele_pfp_momentum_2shwr_intime,
			                                                                               h_ele_pfp_theta_1shwr_intime, h_ele_pfp_theta_2shwr_intime,
			                                                                               h_ele_pfp_phi_1shwr_intime, h_ele_pfp_phi_2shwr_intime);

			//delete at the very end!
			delete passed_tpco_intime;
		}
	}//end in-time cosmic running
//*********************************************************
//*********************************************************
	std::vector<int> * in_time_counter_v = new std::vector<int>;
	in_time_counter_v->resize(24, 0);
	std::vector<int> * pe_counter_v = new std::vector<int>;
	pe_counter_v->resize(24, 0);
	std::vector<int> * reco_nue_counter_v = new std::vector<int>;
	reco_nue_counter_v->resize(24, 0);
	std::vector<int> * in_fv_counter_v = new std::vector<int>;
	in_fv_counter_v->resize(24, 0);
	std::vector<int> * vtx_flash_counter_v = new std::vector<int>;
	vtx_flash_counter_v->resize(24, 0);
	std::vector<int> * shwr_tpco_counter_v = new std::vector<int>;
	shwr_tpco_counter_v->resize(24, 0);
	std::vector<int> * trk_tpco_counter_v = new std::vector<int>;
	trk_tpco_counter_v->resize(24, 0);
	std::vector<int> * hit_threshold_counter_v = new std::vector<int>;
	hit_threshold_counter_v->resize(24, 0);
	std::vector<int> * open_angle_counter_v = new std::vector<int>;
	open_angle_counter_v->resize(24, 0);
	std::vector<int> * dedx_counter_v = new std::vector<int>;
	dedx_counter_v->resize(24, 0);
	std::vector<int> * secondary_shower_counter_v = new std::vector<int>;
	secondary_shower_counter_v->resize(24, 0);
	std::vector<int> * hit_lengthRatio_counter_v = new std::vector<int>;
	hit_lengthRatio_counter_v->resize(24, 0);
	std::vector<int> * hit_threshold_collection_counter_v = new std::vector<int>;
	hit_threshold_collection_counter_v->resize(24, 0);
	std::vector<int> * trk_len_shwr_len_ratio_counter_v = new std::vector<int>;
	trk_len_shwr_len_ratio_counter_v->resize(24, 0);
	std::vector<int> * track_containment_counter_v = new std::vector<int>;
	track_containment_counter_v->resize(24, 0);

	std::vector<std::pair<double, int> > * flash_time = new std::vector<std::pair<double, int> >;

	std::vector<int> * has_track = new std::vector<int>;
	has_track->resize(2, 0);
	std::vector<int> * no_track = new std::vector<int>;
	no_track->resize(2, 0);
	std::vector<int> * _1_shwr = new std::vector<int>;
	_1_shwr->resize(2, 0);
	std::vector<int> * _2_shwr = new std::vector<int>;
	_2_shwr->resize(2, 0);
	std::vector<int> * _3_shwr = new std::vector<int>;
	_3_shwr->resize(2, 0);
	std::vector<int> * _4_shwr = new std::vector<int>;
	_4_shwr->resize(2, 0);

	std::vector<TH1 * > * h_ele_pfp_xyz_nue_cc         = new std::vector<TH1 * >;
	h_ele_pfp_xyz_nue_cc->push_back(h_ele_pfp_x_nue_cc);
	h_ele_pfp_xyz_nue_cc->push_back(h_ele_pfp_y_nue_cc);
	h_ele_pfp_xyz_nue_cc->push_back(h_ele_pfp_z_nue_cc);
	std::vector<TH1 * > * h_ele_pfp_xyz_nue_cc_out_fv  = new std::vector<TH1 * >;
	h_ele_pfp_xyz_nue_cc_out_fv->push_back(h_ele_pfp_x_nue_cc_out_fv);
	h_ele_pfp_xyz_nue_cc_out_fv->push_back(h_ele_pfp_y_nue_cc_out_fv);
	h_ele_pfp_xyz_nue_cc_out_fv->push_back(h_ele_pfp_z_nue_cc_out_fv);
	std::vector<TH1 * > * h_ele_pfp_xyz_nue_cc_mixed   = new std::vector<TH1 * >;
	h_ele_pfp_xyz_nue_cc_mixed->push_back(h_ele_pfp_x_nue_cc_mixed);
	h_ele_pfp_xyz_nue_cc_mixed->push_back(h_ele_pfp_y_nue_cc_mixed);
	h_ele_pfp_xyz_nue_cc_mixed->push_back(h_ele_pfp_z_nue_cc_mixed);
	std::vector<TH1 * > * h_ele_pfp_xyz_numu_cc        = new std::vector<TH1 * >;
	h_ele_pfp_xyz_numu_cc->push_back(h_ele_pfp_x_numu_cc);
	h_ele_pfp_xyz_numu_cc->push_back(h_ele_pfp_y_numu_cc);
	h_ele_pfp_xyz_numu_cc->push_back(h_ele_pfp_z_numu_cc);
	std::vector<TH1 * > * h_ele_pfp_xyz_numu_cc_mixed  = new std::vector<TH1 * >;
	h_ele_pfp_xyz_numu_cc_mixed->push_back(h_ele_pfp_x_numu_cc_mixed);
	h_ele_pfp_xyz_numu_cc_mixed->push_back(h_ele_pfp_y_numu_cc_mixed);
	h_ele_pfp_xyz_numu_cc_mixed->push_back(h_ele_pfp_z_numu_cc_mixed);
	std::vector<TH1 * > * h_ele_pfp_xyz_nc             = new std::vector<TH1 * >;
	h_ele_pfp_xyz_nc->push_back(h_ele_pfp_x_nc);
	h_ele_pfp_xyz_nc->push_back(h_ele_pfp_y_nc);
	h_ele_pfp_xyz_nc->push_back(h_ele_pfp_z_nc);
	std::vector<TH1 * > * h_ele_pfp_xyz_nc_pi0         = new std::vector<TH1 * >;
	h_ele_pfp_xyz_nc_pi0->push_back(h_ele_pfp_x_nc_pi0);
	h_ele_pfp_xyz_nc_pi0->push_back(h_ele_pfp_y_nc_pi0);
	h_ele_pfp_xyz_nc_pi0->push_back(h_ele_pfp_z_nc_pi0);
	std::vector<TH1 * > * h_ele_pfp_xyz_cosmic         = new std::vector<TH1 * >;
	h_ele_pfp_xyz_cosmic->push_back(h_ele_pfp_x_cosmic);
	h_ele_pfp_xyz_cosmic->push_back(h_ele_pfp_y_cosmic);
	h_ele_pfp_xyz_cosmic->push_back(h_ele_pfp_z_cosmic);
	std::vector<TH1 * > * h_ele_pfp_xyz_other_mixed    = new std::vector<TH1 * >;
	h_ele_pfp_xyz_other_mixed->push_back(h_ele_pfp_x_other_mixed);
	h_ele_pfp_xyz_other_mixed->push_back(h_ele_pfp_y_other_mixed);
	h_ele_pfp_xyz_other_mixed->push_back(h_ele_pfp_z_other_mixed);
	std::vector<TH1 * > * h_ele_pfp_xyz_unmatched      = new std::vector<TH1 * >;
	h_ele_pfp_xyz_unmatched->push_back(h_ele_pfp_x_unmatched);
	h_ele_pfp_xyz_unmatched->push_back(h_ele_pfp_y_unmatched);
	h_ele_pfp_xyz_unmatched->push_back(h_ele_pfp_z_unmatched);

	std::vector<TH1 * > * h_any_pfp_xyz_nue_cc         = new std::vector<TH1 * >;
	h_any_pfp_xyz_nue_cc->push_back(h_any_pfp_x_nue_cc);
	h_any_pfp_xyz_nue_cc->push_back(h_any_pfp_y_nue_cc);
	h_any_pfp_xyz_nue_cc->push_back(h_any_pfp_z_nue_cc);
	std::vector<TH1 * > * h_any_pfp_xyz_nue_cc_out_fv  = new std::vector<TH1 * >;
	h_any_pfp_xyz_nue_cc_out_fv->push_back(h_any_pfp_x_nue_cc_out_fv);
	h_any_pfp_xyz_nue_cc_out_fv->push_back(h_any_pfp_y_nue_cc_out_fv);
	h_any_pfp_xyz_nue_cc_out_fv->push_back(h_any_pfp_z_nue_cc_out_fv);
	std::vector<TH1 * > * h_any_pfp_xyz_nue_cc_mixed   = new std::vector<TH1 * >;
	h_any_pfp_xyz_nue_cc_mixed->push_back(h_any_pfp_x_nue_cc_mixed);
	h_any_pfp_xyz_nue_cc_mixed->push_back(h_any_pfp_y_nue_cc_mixed);
	h_any_pfp_xyz_nue_cc_mixed->push_back(h_any_pfp_z_nue_cc_mixed);
	std::vector<TH1 * > * h_any_pfp_xyz_numu_cc        = new std::vector<TH1 * >;
	h_any_pfp_xyz_numu_cc->push_back(h_any_pfp_x_numu_cc);
	h_any_pfp_xyz_numu_cc->push_back(h_any_pfp_y_numu_cc);
	h_any_pfp_xyz_numu_cc->push_back(h_any_pfp_z_numu_cc);
	std::vector<TH1 * > * h_any_pfp_xyz_numu_cc_mixed  = new std::vector<TH1 * >;
	h_any_pfp_xyz_numu_cc_mixed->push_back(h_any_pfp_x_numu_cc_mixed);
	h_any_pfp_xyz_numu_cc_mixed->push_back(h_any_pfp_y_numu_cc_mixed);
	h_any_pfp_xyz_numu_cc_mixed->push_back(h_any_pfp_z_numu_cc_mixed);
	std::vector<TH1 * > * h_any_pfp_xyz_nc             = new std::vector<TH1 * >;
	h_any_pfp_xyz_nc->push_back(h_any_pfp_x_nc);
	h_any_pfp_xyz_nc->push_back(h_any_pfp_y_nc);
	h_any_pfp_xyz_nc->push_back(h_any_pfp_z_nc);
	std::vector<TH1 * > * h_any_pfp_xyz_nc_pi0         = new std::vector<TH1 * >;
	h_any_pfp_xyz_nc_pi0->push_back(h_any_pfp_x_nc_pi0);
	h_any_pfp_xyz_nc_pi0->push_back(h_any_pfp_y_nc_pi0);
	h_any_pfp_xyz_nc_pi0->push_back(h_any_pfp_z_nc_pi0);
	std::vector<TH1 * > * h_any_pfp_xyz_cosmic         = new std::vector<TH1 * >;
	h_any_pfp_xyz_cosmic->push_back(h_any_pfp_x_cosmic);
	h_any_pfp_xyz_cosmic->push_back(h_any_pfp_y_cosmic);
	h_any_pfp_xyz_cosmic->push_back(h_any_pfp_z_cosmic);
	std::vector<TH1 * > * h_any_pfp_xyz_other_mixed    = new std::vector<TH1 * >;
	h_any_pfp_xyz_other_mixed->push_back(h_any_pfp_x_other_mixed);
	h_any_pfp_xyz_other_mixed->push_back(h_any_pfp_y_other_mixed);
	h_any_pfp_xyz_other_mixed->push_back(h_any_pfp_z_other_mixed);
	std::vector<TH1 * > * h_any_pfp_xyz_unmatched      = new std::vector<TH1 * >;
	h_any_pfp_xyz_unmatched->push_back(h_any_pfp_x_unmatched);
	h_any_pfp_xyz_unmatched->push_back(h_any_pfp_y_unmatched);
	h_any_pfp_xyz_unmatched->push_back(h_any_pfp_z_unmatched);

	std::vector<TH1 * > * h_any_pfp_xyz_last_nue_cc         = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_nue_cc->push_back(h_any_pfp_x_last_nue_cc);
	h_any_pfp_xyz_last_nue_cc->push_back(h_any_pfp_y_last_nue_cc);
	h_any_pfp_xyz_last_nue_cc->push_back(h_any_pfp_z_last_nue_cc);
	std::vector<TH1 * > * h_any_pfp_xyz_last_nue_cc_out_fv  = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_nue_cc_out_fv->push_back(h_any_pfp_x_last_nue_cc_out_fv);
	h_any_pfp_xyz_last_nue_cc_out_fv->push_back(h_any_pfp_y_last_nue_cc_out_fv);
	h_any_pfp_xyz_last_nue_cc_out_fv->push_back(h_any_pfp_z_last_nue_cc_out_fv);
	std::vector<TH1 * > * h_any_pfp_xyz_last_nue_cc_mixed   = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_nue_cc_mixed->push_back(h_any_pfp_x_last_nue_cc_mixed);
	h_any_pfp_xyz_last_nue_cc_mixed->push_back(h_any_pfp_y_last_nue_cc_mixed);
	h_any_pfp_xyz_last_nue_cc_mixed->push_back(h_any_pfp_z_last_nue_cc_mixed);
	std::vector<TH1 * > * h_any_pfp_xyz_last_numu_cc        = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_numu_cc->push_back(h_any_pfp_x_last_numu_cc);
	h_any_pfp_xyz_last_numu_cc->push_back(h_any_pfp_y_last_numu_cc);
	h_any_pfp_xyz_last_numu_cc->push_back(h_any_pfp_z_last_numu_cc);
	std::vector<TH1 * > * h_any_pfp_xyz_last_numu_cc_mixed  = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_numu_cc_mixed->push_back(h_any_pfp_x_last_numu_cc_mixed);
	h_any_pfp_xyz_last_numu_cc_mixed->push_back(h_any_pfp_y_last_numu_cc_mixed);
	h_any_pfp_xyz_last_numu_cc_mixed->push_back(h_any_pfp_z_last_numu_cc_mixed);
	std::vector<TH1 * > * h_any_pfp_xyz_last_nc             = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_nc->push_back(h_any_pfp_x_last_nc);
	h_any_pfp_xyz_last_nc->push_back(h_any_pfp_y_last_nc);
	h_any_pfp_xyz_last_nc->push_back(h_any_pfp_z_last_nc);
	std::vector<TH1 * > * h_any_pfp_xyz_last_nc_pi0         = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_nc_pi0->push_back(h_any_pfp_x_last_nc_pi0);
	h_any_pfp_xyz_last_nc_pi0->push_back(h_any_pfp_y_last_nc_pi0);
	h_any_pfp_xyz_last_nc_pi0->push_back(h_any_pfp_z_last_nc_pi0);
	std::vector<TH1 * > * h_any_pfp_xyz_last_cosmic         = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_cosmic->push_back(h_any_pfp_x_last_cosmic);
	h_any_pfp_xyz_last_cosmic->push_back(h_any_pfp_y_last_cosmic);
	h_any_pfp_xyz_last_cosmic->push_back(h_any_pfp_z_last_cosmic);
	std::vector<TH1 * > * h_any_pfp_xyz_last_other_mixed    = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_other_mixed->push_back(h_any_pfp_x_last_other_mixed);
	h_any_pfp_xyz_last_other_mixed->push_back(h_any_pfp_y_last_other_mixed);
	h_any_pfp_xyz_last_other_mixed->push_back(h_any_pfp_z_last_other_mixed);
	std::vector<TH1 * > * h_any_pfp_xyz_last_unmatched      = new std::vector<TH1 * >;
	h_any_pfp_xyz_last_unmatched->push_back(h_any_pfp_x_last_unmatched);
	h_any_pfp_xyz_last_unmatched->push_back(h_any_pfp_y_last_unmatched);
	h_any_pfp_xyz_last_unmatched->push_back(h_any_pfp_z_last_unmatched);

	//Event, Run, VtxX, VtxY, VtxZ, pass/fail reason
	// std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v
	//         = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> >;
	std::vector<std::tuple<int, int, int, double, double, double, std::string, std::string, int, int, double> > * post_open_angle_cuts_v
	        = new std::vector<std::tuple<int, int, int, double, double, double, std::string, std::string, int, int, double> >;

	std::cout << "========================" << std::endl;
	std::cout << "== Begin MC Selection ==" << std::endl;
	std::cout << "========================" << std::endl;

	const int total_entries = mytree->GetEntries();
	std::cout << "Total Events     : " << total_entries << std::endl;
	std::cout << "Total Events (MC): " << mctruth_counter_tree->GetEntries() << std::endl;

	std::vector<int> * passed_runs = new std::vector<int>;
	//passed runs is filled with 0, 1, or 2
	//0 = not in time
	//1 = passed - in time and PE threshold
	// 2 = in-time, but not enough PE -- this counts against my efficiency
	passed_runs->resize(total_entries);

	_cuts_instance.selection_cuts::loop_flashes(f, optree, flash_pe_threshold, flash_time_start, flash_time_end, passed_runs, flash_time, 2);
	run_sum = 0;
	out_of_time_sum = 0;
	low_pe_sum = 0;
	for(auto const run : * passed_runs)
	{
		if(run == 1) {run_sum++; }
		if(run == 0) {out_of_time_sum++; }
		if(run == 2) {low_pe_sum++; }
	}
	std::cout << " -------------------------------------- " << std::endl;
	std::cout << "Passed Runs Vector Size: " << passed_runs->size() << std::endl;
	std::cout << "Number Events In-Time & >= 50 PE: " << run_sum << std::endl;
	std::cout << "Number Events Not In-Time       : " << out_of_time_sum << std::endl;
	std::cout << "Number Events In-Time & <= 50 PE: " << low_pe_sum << std::endl;
	std::cout << " -------------------------------------- " << std::endl;

	//get vector with largest flashes y,z positions
	std::vector< std::vector< double> > * largest_flash_v_v = new std::vector < std::vector < double > >;
	_cuts_instance.selection_cuts::SetXYflashVector(f, optree, largest_flash_v_v, flash_time_start, flash_time_end, flash_pe_threshold);

	int test_mc_nue_cc_counter = 0;
	int test_mc_numu_cc_counter = 0;
	int test_mc_nue_nc_counter = 0;
	int test_mc_numu_nc_counter = 0;
	int test_mc_nue_cc_counter_bar = 0;
	int test_mc_numu_cc_counter_bar = 0;
	int test_mc_nue_nc_counter_bar = 0;
	int test_mc_numu_nc_counter_bar = 0;

	int infv_test_mc_nue_cc_counter = 0;
	int infv_test_mc_numu_cc_counter = 0;
	int infv_test_mc_nue_nc_counter = 0;
	int infv_test_mc_numu_nc_counter = 0;
	int infv_test_mc_nue_cc_counter_bar = 0;
	int infv_test_mc_numu_cc_counter_bar = 0;
	int infv_test_mc_nue_nc_counter_bar = 0;
	int infv_test_mc_numu_nc_counter_bar = 0;

	//**********************************
	//now let's do the cuts
	//*********************************
	for(int event = 0; event < total_entries; event++)
	{
		if(_verbose)
		{
			std::cout << "----------------------" << std::endl;
			std::cout << "[EVENT NUMBER] \t " << event << std::endl;
			std::cout << "----------------------" << std::endl;
		}
		mytree->GetEntry(event);
		mctruth_counter_tree->GetEntry(event);

		//********************************
		//before Any cuts!!!
		//********************************

		//check if nue interaction has true vtx in TPC
		//std::cout << tpc_object_container_v->size() << std::endl;
		// bool true_in_tpc = false;
		// for (int i = 0; i < tpc_object_container_v->size(); i++)
		// {
		//      auto const tpc_object_container = tpc_object_container_v->at(i);
		//      double _mc_nu_vtx_x = 0.;
		//      double _mc_nu_vtx_y = 0.;
		//      double _mc_nu_vtx_z = 0.;
		//      _mc_nu_vtx_x = tpc_object_container.mcVtxX();
		//      _mc_nu_vtx_y = tpc_object_container.mcVtxY();
		//      _mc_nu_vtx_z = tpc_object_container.mcVtxZ();
		//
		//      true_in_tpc = _cuts_instance.selection_cuts::in_fv(_mc_nu_vtx_x, _mc_nu_vtx_y, _mc_nu_vtx_z, fv_boundary_v);
		//      if(true_in_tpc == true) {break; }
		// }
		const bool true_in_tpc = true_in_tpc_v.at(event);

		if(mc_nu_id == 1) {test_mc_nue_cc_counter++; }
		if(mc_nu_id == 2) {test_mc_numu_cc_counter++; }
		if(mc_nu_id == 3) {test_mc_nue_nc_counter++; }
		if(mc_nu_id == 4) {test_mc_numu_nc_counter++; }
		if(mc_nu_id == 5) {test_mc_nue_cc_counter_bar++; }
		if(mc_nu_id == 6) {test_mc_numu_cc_counter_bar++; }
		if(mc_nu_id == 7) {test_mc_nue_nc_counter_bar++; }
		if(mc_nu_id == 8) {test_mc_numu_nc_counter_bar++; }

		if(mc_nu_id == 1 && true_in_tpc == true) {infv_test_mc_nue_cc_counter++; }
		if(mc_nu_id == 2 && true_in_tpc == true) {infv_test_mc_numu_cc_counter++; }
		if(mc_nu_id == 3 && true_in_tpc == true) {infv_test_mc_nue_nc_counter++; }
		if(mc_nu_id == 4 && true_in_tpc == true) {infv_test_mc_numu_nc_counter++; }
		if(mc_nu_id == 5 && true_in_tpc == true) {infv_test_mc_nue_cc_counter_bar++; }
		if(mc_nu_id == 6 && true_in_tpc == true) {infv_test_mc_numu_cc_counter_bar++; }
		if(mc_nu_id == 7 && true_in_tpc == true) {infv_test_mc_nue_nc_counter_bar++; }
		if(mc_nu_id == 8 && true_in_tpc == true) {infv_test_mc_numu_nc_counter_bar++; }


		//now we apply the classifier to all TPC Objects in this event
		std::vector<std::pair<std::string, int> > * tpco_classifier_v = new std::vector<std::pair<std::string, int> >;
		_functions_instance.selection_functions::FillTPCOClassV(tpc_object_container_v, true_in_tpc, has_pi0, tpco_classifier_v);

		//YZ Position of largest flash
		std::vector < double > largest_flash_v = largest_flash_v_v->at(event);
		//control for poorly reco flashes
		if(largest_flash_v.at(1) != 0) {h_flash_z_mc->Fill(largest_flash_v.at(1)); }

		//List of TPC Objects which pass the cuts
		std::vector<std::pair<int, std::string> > * passed_tpco = new std::vector<std::pair<int, std::string> >;
		std::vector<std::pair<int, std::string> > * dummy_passed_tpco = new std::vector<std::pair<int, std::string> >;
		passed_tpco->resize(tpc_object_container_v->size());
		dummy_passed_tpco->resize(tpc_object_container_v->size());
		for(int i = 0; i < passed_tpco->size(); i++)
		{
			passed_tpco->at(i).first = 1;
			passed_tpco->at(i).second = "Passed";
			dummy_passed_tpco->at(i).first = 1;
			dummy_passed_tpco->at(i).second = "Passed";
		}

		double mc_cos_theta = -999;
		if(mc_nu_momentum != 0) {mc_cos_theta = mc_nu_dir_z; }
		const double mc_phi       = atan2(mc_nu_dir_y, mc_nu_dir_x);
		double mc_ele_cos_theta = -999;
		double mc_ele_theta = -999;
		if(mc_ele_momentum != 0)
		{
			mc_ele_cos_theta = mc_ele_dir_z;
			mc_ele_theta = acos(mc_ele_dir_z) * (180/3.1415);
		}
		const double mc_ele_phi = atan2(mc_ele_dir_y, mc_ele_dir_x);
		if(true_in_tpc == true)
		{
			h_all_true_energy_theta->Fill(mc_nu_energy, acos(mc_cos_theta) * (180 / 3.1415));
			if(mc_nu_id == 1 || mc_nu_id == 5)
			{
				h_nue_true_energy_theta->Fill(mc_nu_energy, acos(mc_cos_theta) * (180 / 3.1415));
				h_nue_true_energy_phi->Fill(mc_nu_energy, mc_phi * (180 / 3.1415));
				h_ele_true_energy_theta->Fill(mc_ele_energy, acos(mc_ele_cos_theta) * (180 / 3.1415));
				h_ele_true_energy_phi->Fill(mc_ele_energy, mc_ele_phi * (180 / 3.1415));
			}
		}
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			h_nue_eng_eff_den->Fill(mc_nu_energy);
			h_ele_eng_eff_den->Fill(mc_ele_energy);
			h_nue_vtx_x_eff_den->Fill(mc_nu_vtx_x);
			h_nue_vtx_y_eff_den->Fill(mc_nu_vtx_y);
			h_nue_vtx_z_eff_den->Fill(mc_nu_vtx_z);
			h_nue_dir_x_eff_den->Fill(mc_nu_dir_x);
			h_nue_dir_y_eff_den->Fill(mc_nu_dir_y);
			h_nue_dir_z_eff_den->Fill(mc_nu_dir_z);
			h_ele_dir_x_eff_den->Fill(mc_ele_dir_x);
			h_ele_dir_y_eff_den->Fill(mc_ele_dir_y);
			h_ele_dir_z_eff_den->Fill(mc_ele_dir_z);
			h_ele_theta_eff_den->Fill(mc_ele_theta);
			h_nue_num_part_eff_den->Fill(mc_nu_num_particles);
			h_nue_num_chrg_part_eff_den->Fill(mc_nu_num_charged_particles);
			h_nue_cos_theta_eff_den->Fill(mc_cos_theta);
			h_nue_phi_eff_den->Fill(mc_phi * (180 / 3.1415));
			h_ele_cos_theta_eff_den->Fill(mc_ele_cos_theta);
			h_ele_phi_eff_den->Fill(mc_ele_phi * (180 / 3.1415));
			h_nue_true_theta->Fill( acos(mc_cos_theta) * (180 / 3.1415));
			h_nue_true_phi->Fill(mc_phi * (180 / 3.1415));
			h_nue_true_theta_phi->Fill(mc_phi * (180 / 3.1415), acos(mc_cos_theta) * (180 / 3.1415));
		}
		//********************************
		//begin cuts!!!
		//********************************

		//***********************************************************
		//this is where the in-time optical cut actually takes effect
		//***********************************************************
		if(passed_runs->at(event) == 0)
		{
			if(_verbose) std::cout << "[Failed In-Time Cut]" << std::endl;
			continue;
		}//false

		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, in_time_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_no_cut, h_selected_ele_energy_no_cut);

		_functions_instance.selection_functions::XYZPosition(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                     mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                     h_any_pfp_xyz_nue_cc,
		                                                     h_any_pfp_xyz_nue_cc_out_fv,
		                                                     h_any_pfp_xyz_nue_cc_mixed,
		                                                     h_any_pfp_xyz_numu_cc,
		                                                     h_any_pfp_xyz_numu_cc_mixed,
		                                                     h_any_pfp_xyz_nc,
		                                                     h_any_pfp_xyz_nc_pi0,
		                                                     h_any_pfp_xyz_cosmic,
		                                                     h_any_pfp_xyz_other_mixed,
		                                                     h_any_pfp_xyz_unmatched,
		                                                     h_pfp_zy_vtx_nue_cc,
		                                                     h_pfp_zy_vtx_all, xyz_near_mc, xyz_far_mc);

		//***********************************************************
		//this is where the pe optical cut takes effect
		//***********************************************************
		//PE threshold cut
		if(passed_runs->at(event) == 2)
		{
			if(_verbose) std::cout << "[Passed In-Time Cut] [Failed PE Threshold] " << std::endl;
			continue;
		}
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, pe_counter_v);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_num->Fill(mc_ele_energy); }
		//****************************
		// ****** reco nue cut *******
		//****************************
		_cuts_instance.selection_cuts::HasNue(tpc_object_container_v, passed_tpco, _verbose);

		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, reco_nue_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_reco_nue, h_selected_ele_energy_reco_nue);
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_1, h_reco_ele_e_1, h_mc_reco_ele_e_1);
			// _functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, dummy_passed_tpco, tpco_classifier_v, mc_ele_energy,
			//                                                             h_mc_ele_e_0, h_reco_ele_e_0, h_mc_reco_ele_e_0);
		}
		if(mc_nu_id == 1 || mc_nu_id == 5) {
			if(true_in_tpc == true)
			{
				int pass = 0;
				bool case_1 = false;
				bool case_2 = false;
				for(auto const passed_tpco_pair : * passed_tpco) {pass += passed_tpco_pair.first; }
				if(pass > 0) {case_1 = true; } //lets more pass?
				if(tabulated_origins->at(0) >= 1) {case_2 = true; } //lets less pass, not reflected in printed values
				// if(case_1 && !case_2)
				// {
				//      std::cout << "mc_nu_id: " << mc_nu_id << std::endl;
				//      int num_tpc_obj = 0;
				//      for(const auto tpco_classifier : * tpco_classifier_v)
				//      {
				//              auto const tpc_obj = tpc_object_container_v->at(num_tpc_obj);
				//
				//              if(tpco_classifier.first != "cosmic")
				//              {
				//                      std::cout << '\t' << "classifier_id: " << tpco_classifier.first;
				//                      std::cout << ", Mode: " << tpc_obj.Mode() << ", PFP_PDG: " << tpc_obj.PFParticlePdgCode();
				//                      std::cout << ", CCNC: " << tpc_obj.CCNC() << ", HasMCPi0: " << tpc_obj.HasMCPi0() << std::endl;
				//                      const int n_pfp = tpc_obj.NumPFParticles();
				//                      for(int j = 0; j < n_pfp; j++)
				//                      {
				//                              auto const part = tpc_obj.GetParticle(j);
				//                              //const int n_pfp_hits = part.NumPFPHits();
				//                              const int mc_parent_pdg = part.MCParentPdg();
				//                              const int pfp_pdg = part.PFParticlePdgCode();
				//                              std::cout << '\t' << '\t' << "MC Parent PDG: " << mc_parent_pdg << ", PFP PDG: " << pfp_pdg << std::endl;
				//                      }
				//              }
				//              num_tpc_obj++;
				//      }
				//      num_debug_events++;
				// }
			}
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_nue_cut_nue_cc,
		                                                                 h_ele_momentum_nue_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_nue_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_nue_cut_numu_cc,
		                                                                 h_ele_momentum_nue_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_nue_cut_nc,
		                                                                 h_ele_momentum_nue_cut_nc_pi0,
		                                                                 h_ele_momentum_nue_cut_cosmic,
		                                                                 h_ele_momentum_nue_cut_other_mixed,
		                                                                 h_ele_momentum_nue_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_nue_cut_electron,
		                                                                             h_leading_momentum_nue_cut_photon,
		                                                                             h_leading_momentum_nue_cut_proton,
		                                                                             h_leading_momentum_nue_cut_pion,
		                                                                             h_leading_momentum_nue_cut_muon,
		                                                                             h_leading_momentum_nue_cut_kaon,
		                                                                             h_leading_momentum_nue_cut_neutron,
		                                                                             h_leading_momentum_nue_cut_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_nue_cut_nue_cc,
		                                                           h_multiplicity_shower_nue_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_nue_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_nue_cut_numu_cc,
		                                                           h_multiplicity_shower_nue_cut_nc,
		                                                           h_multiplicity_shower_nue_cut_nc_pi0,
		                                                           h_multiplicity_shower_nue_cut_cosmic,
		                                                           h_multiplicity_shower_nue_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_nue_cut_other_mixed,
		                                                           h_multiplicity_shower_nue_cut_unmatched,
		                                                           h_multiplicity_track_nue_cut_nue_cc,
		                                                           h_multiplicity_track_nue_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_nue_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_nue_cut_numu_cc,
		                                                           h_multiplicity_track_nue_cut_nc,
		                                                           h_multiplicity_track_nue_cut_nc_pi0,
		                                                           h_multiplicity_track_nue_cut_cosmic,
		                                                           h_multiplicity_track_nue_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_nue_cut_other_mixed,
		                                                           h_multiplicity_track_nue_cut_unmatched);


		//** Testing flash vs neutrino interaction for origin **
		_functions_instance.selection_functions::FlashTot0(largest_flash_v, mc_nu_time, mc_nu_id, tabulated_origins,
		                                                   fv_boundary_v, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,h_flash_t0_diff);

		//** Testing leading shower length vs hits **//
		_functions_instance.selection_functions::ShowerLengthvsHits(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                            h_shwr_len_hits_nue_cc, h_shwr_len_hits_nue_cc_out_fv,
		                                                            h_shwr_len_hits_nue_cc_mixed, h_shwr_len_hits_numu_cc,
		                                                            h_shwr_len_hits_numu_cc_mixed, h_shwr_len_hits_nc,
		                                                            h_shwr_len_hits_nc_pi0, h_shwr_len_hits_cosmic,
		                                                            h_shwr_len_hits_other_mixed, h_shwr_len_hits_unmatched);

		//pre most cuts hits
		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins->at(0) >= 1 && true_in_tpc == true)
		{
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
			                                                             mc_nu_energy, mc_ele_energy, h_shwr_hits_nu_eng, h_shwr_hits_ele_eng,
			                                                             h_shwr_collection_hits_nu_eng, h_shwr_collection_hits_ele_eng);
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
			                                                             mc_nu_energy, mc_ele_energy, h_shwr_hits_nu_eng_zoom, h_shwr_hits_ele_eng_zoom,
			                                                             h_shwr_collection_hits_nu_eng_zoom, h_shwr_collection_hits_ele_eng_zoom);

			_functions_instance.selection_functions::TrueRecoEle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
			                                                     mc_ele_momentum, mc_ele_cos_theta,
			                                                     h_true_reco_ele_momentum_pre, h_true_reco_ele_costheta_pre, h_true_num_e_pre);
		}

		_functions_instance.selection_functions::TopologyPlots1(tpc_object_container_v, passed_tpco, tpco_classifier_v,
		                                                        h_pfp_track_shower_nue_cc_qe, h_pfp_track_shower_nue_cc_out_fv,
		                                                        h_pfp_track_shower_nue_cc_res, h_pfp_track_shower_nue_cc_dis,
		                                                        h_pfp_track_shower_nue_cc_coh, h_pfp_track_shower_nue_cc_mec,
		                                                        h_pfp_track_shower_nc, h_pfp_track_shower_numu_cc_qe,
		                                                        h_pfp_track_shower_numu_cc_res, h_pfp_track_shower_numu_cc_dis,
		                                                        h_pfp_track_shower_numu_cc_coh, h_pfp_track_shower_numu_cc_mec,
		                                                        h_pfp_track_shower_nc_pi0, h_pfp_track_shower_nue_cc_mixed,
		                                                        h_pfp_track_shower_numu_cc_mixed, h_pfp_track_shower_cosmic,
		                                                        h_pfp_track_shower_other_mixed, h_pfp_track_shower_unmatched,
		                                                        h_leading_shower_mc_pdg_nue_cc_qe, h_leading_shower_mc_pdg_nue_cc_out_fv,
		                                                        h_leading_shower_mc_pdg_nue_cc_res, h_leading_shower_mc_pdg_nue_cc_dis,
		                                                        h_leading_shower_mc_pdg_nue_cc_coh, h_leading_shower_mc_pdg_nue_cc_mec,
		                                                        h_leading_shower_mc_pdg_nc, h_leading_shower_mc_pdg_numu_cc_qe,
		                                                        h_leading_shower_mc_pdg_numu_cc_res, h_leading_shower_mc_pdg_numu_cc_dis,
		                                                        h_leading_shower_mc_pdg_numu_cc_coh, h_leading_shower_mc_pdg_numu_cc_mec,
		                                                        h_leading_shower_mc_pdg_nc_pi0, h_leading_shower_mc_pdg_nue_cc_mixed,
		                                                        h_leading_shower_mc_pdg_numu_cc_mixed, h_leading_shower_mc_pdg_cosmic,
		                                                        h_leading_shower_mc_pdg_other_mixed, h_leading_shower_mc_pdg_unmatched,
		                                                        h_pfp_track_nue_cc_qe, h_pfp_track_nue_cc_out_fv,
		                                                        h_pfp_track_nue_cc_res, h_pfp_track_nue_cc_dis,
		                                                        h_pfp_track_nue_cc_coh, h_pfp_track_nue_cc_mec,
		                                                        h_pfp_track_nc, h_pfp_track_numu_cc_qe,
		                                                        h_pfp_track_numu_cc_res, h_pfp_track_numu_cc_dis,
		                                                        h_pfp_track_numu_cc_coh, h_pfp_track_numu_cc_mec,
		                                                        h_pfp_track_nc_pi0, h_pfp_track_nue_cc_mixed,
		                                                        h_pfp_track_numu_cc_mixed, h_pfp_track_cosmic,
		                                                        h_pfp_track_other_mixed, h_pfp_track_unmatched,
		                                                        h_pfp_shower_nue_cc_qe, h_pfp_shower_nue_cc_out_fv,
		                                                        h_pfp_shower_nue_cc_res, h_pfp_shower_nue_cc_dis,
		                                                        h_pfp_shower_nue_cc_coh, h_pfp_shower_nue_cc_mec,
		                                                        h_pfp_shower_nc, h_pfp_shower_numu_cc_qe,
		                                                        h_pfp_shower_numu_cc_res, h_pfp_shower_numu_cc_dis,
		                                                        h_pfp_shower_numu_cc_coh, h_pfp_shower_numu_cc_mec,
		                                                        h_pfp_shower_nc_pi0, h_pfp_shower_nue_cc_mixed,
		                                                        h_pfp_shower_numu_cc_mixed, h_pfp_shower_cosmic,
		                                                        h_pfp_shower_other_mixed, h_pfp_shower_unmatched);

		_functions_instance.selection_functions::XYZPosition(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                     mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                     h_ele_pfp_xyz_nue_cc,
		                                                     h_ele_pfp_xyz_nue_cc_out_fv,
		                                                     h_ele_pfp_xyz_nue_cc_mixed,
		                                                     h_ele_pfp_xyz_numu_cc,
		                                                     h_ele_pfp_xyz_numu_cc_mixed,
		                                                     h_ele_pfp_xyz_nc,
		                                                     h_ele_pfp_xyz_nc_pi0,
		                                                     h_ele_pfp_xyz_cosmic,
		                                                     h_ele_pfp_xyz_other_mixed,
		                                                     h_ele_pfp_xyz_unmatched,
		                                                     h_mc_vtx_xy_nue_cc,
		                                                     h_mc_vtx_xz_nue_cc,
		                                                     h_mc_vtx_yz_nue_cc,
		                                                     h_reco_vtx_xy_nue_cc,
		                                                     h_reco_vtx_xz_nue_cc,
		                                                     h_reco_vtx_yz_nue_cc,
		                                                     h_mc_vtx_xy_nue_cc_out_fv,
		                                                     h_mc_vtx_xz_nue_cc_out_fv,
		                                                     h_mc_vtx_yz_nue_cc_out_fv,
		                                                     h_reco_vtx_xy_nue_cc_out_fv,
		                                                     h_reco_vtx_xz_nue_cc_out_fv,
		                                                     h_reco_vtx_yz_nue_cc_out_fv,
		                                                     h_mc_reco_vtx_x_nue_cc,
		                                                     h_mc_reco_vtx_y_nue_cc,
		                                                     h_mc_reco_vtx_z_nue_cc,
		                                                     h_mc_reco_vtx_x_nue_cc_out_fv,
		                                                     h_mc_reco_vtx_y_nue_cc_out_fv,
		                                                     h_mc_reco_vtx_z_nue_cc_out_fv);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_reco_nue->Fill(mc_ele_energy); }
		//************************
		//******** in fv cut *****
		//************************
		_cuts_instance.selection_cuts::fiducial_volume_cut(tpc_object_container_v, fv_boundary_v, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, in_fv_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_in_fv, h_selected_ele_energy_in_fv);
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_2, h_reco_ele_e_2, h_mc_reco_ele_e_2);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_fv_cut_nue_cc,
		                                                                 h_ele_momentum_fv_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_fv_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_fv_cut_numu_cc,
		                                                                 h_ele_momentum_fv_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_fv_cut_nc,
		                                                                 h_ele_momentum_fv_cut_nc_pi0,
		                                                                 h_ele_momentum_fv_cut_cosmic,
		                                                                 h_ele_momentum_fv_cut_other_mixed,
		                                                                 h_ele_momentum_fv_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_fv_cut_electron,
		                                                                             h_leading_momentum_fv_cut_photon,
		                                                                             h_leading_momentum_fv_cut_proton,
		                                                                             h_leading_momentum_fv_cut_pion,
		                                                                             h_leading_momentum_fv_cut_muon,
		                                                                             h_leading_momentum_fv_cut_kaon,
		                                                                             h_leading_momentum_fv_cut_neutron,
		                                                                             h_leading_momentum_fv_cut_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_fv_cut_nue_cc,
		                                                           h_multiplicity_shower_fv_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_fv_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_fv_cut_numu_cc,
		                                                           h_multiplicity_shower_fv_cut_nc,
		                                                           h_multiplicity_shower_fv_cut_nc_pi0,
		                                                           h_multiplicity_shower_fv_cut_cosmic,
		                                                           h_multiplicity_shower_fv_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_fv_cut_other_mixed,
		                                                           h_multiplicity_shower_fv_cut_unmatched,
		                                                           h_multiplicity_track_fv_cut_nue_cc,
		                                                           h_multiplicity_track_fv_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_fv_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_fv_cut_numu_cc,
		                                                           h_multiplicity_track_fv_cut_nc,
		                                                           h_multiplicity_track_fv_cut_nc_pi0,
		                                                           h_multiplicity_track_fv_cut_cosmic,
		                                                           h_multiplicity_track_fv_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_fv_cut_other_mixed,
		                                                           h_multiplicity_track_fv_cut_unmatched);

		//we also want to look at the cos(theta) and energy efficiency before we make selection cuts
		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins->at(0) >= 1 && true_in_tpc == true)
		{
			h_ele_eng_eff_num_pre_cuts->Fill(mc_ele_energy);
			h_ele_cos_theta_eff_num_pre_cuts->Fill(mc_ele_cos_theta);
		}
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_in_fv->Fill(mc_ele_energy); }
		//*****************************
		//**** vertex to flash cut ****
		//*****************************
		_functions_instance.selection_functions::PostCutsVtxFlash(largest_flash_v, tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                          h_vtx_flash_nue_cc, h_vtx_flash_nue_cc_out_fv, h_vtx_flash_nue_cc_mixed,
		                                                          h_vtx_flash_numu_cc, h_vtx_flash_nc,
		                                                          h_vtx_flash_cosmic, h_vtx_flash_nc_pi0,
		                                                          h_vtx_flash_numu_cc_mixed, h_vtx_flash_other_mixed,
		                                                          h_vtx_flash_unmatched);
		_functions_instance.selection_functions::PostCutsVtxFlashUpstream(largest_flash_v, tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                  h_vtx_flash_upstream_nue_cc,        h_vtx_flash_upstream_nue_cc_out_fv,
		                                                                  h_vtx_flash_upstream_nue_cc_mixed,
		                                                                  h_vtx_flash_upstream_numu_cc,       h_vtx_flash_upstream_nc,
		                                                                  h_vtx_flash_upstream_cosmic,        h_vtx_flash_upstream_nc_pi0,
		                                                                  h_vtx_flash_upstream_numu_cc_mixed, h_vtx_flash_upstream_other_mixed,
		                                                                  h_vtx_flash_upstream_unmatched);
		_functions_instance.selection_functions::PostCutsVtxFlashDownstream(largest_flash_v, tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                    h_vtx_flash_downstream_nue_cc,        h_vtx_flash_downstream_nue_cc_out_fv,
		                                                                    h_vtx_flash_downstream_nue_cc_mixed,
		                                                                    h_vtx_flash_downstream_numu_cc,       h_vtx_flash_downstream_nc,
		                                                                    h_vtx_flash_downstream_cosmic,        h_vtx_flash_downstream_nc_pi0,
		                                                                    h_vtx_flash_downstream_numu_cc_mixed, h_vtx_flash_downstream_other_mixed,
		                                                                    h_vtx_flash_downstream_unmatched);

		_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, tpc_object_container_v, tolerance, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, vtx_flash_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_vtx_flash, h_selected_ele_energy_vtx_flash);
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_3, h_reco_ele_e_3, h_mc_reco_ele_e_3);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_flash_vtx_cut_nue_cc,
		                                                                 h_ele_momentum_flash_vtx_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_flash_vtx_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_flash_vtx_cut_numu_cc,
		                                                                 h_ele_momentum_flash_vtx_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_flash_vtx_cut_nc,
		                                                                 h_ele_momentum_flash_vtx_cut_nc_pi0,
		                                                                 h_ele_momentum_flash_vtx_cut_cosmic,
		                                                                 h_ele_momentum_flash_vtx_cut_other_mixed,
		                                                                 h_ele_momentum_flash_vtx_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_flash_vtx_electron,
		                                                                             h_leading_momentum_flash_vtx_photon,
		                                                                             h_leading_momentum_flash_vtx_proton,
		                                                                             h_leading_momentum_flash_vtx_pion,
		                                                                             h_leading_momentum_flash_vtx_muon,
		                                                                             h_leading_momentum_flash_vtx_kaon,
		                                                                             h_leading_momentum_flash_vtx_neutron,
		                                                                             h_leading_momentum_flash_vtx_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_flash_vtx_cut_nue_cc,
		                                                           h_multiplicity_shower_flash_vtx_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_flash_vtx_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_flash_vtx_cut_numu_cc,
		                                                           h_multiplicity_shower_flash_vtx_cut_nc,
		                                                           h_multiplicity_shower_flash_vtx_cut_nc_pi0,
		                                                           h_multiplicity_shower_flash_vtx_cut_cosmic,
		                                                           h_multiplicity_shower_flash_vtx_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_flash_vtx_cut_other_mixed,
		                                                           h_multiplicity_shower_flash_vtx_cut_unmatched,
		                                                           h_multiplicity_track_flash_vtx_cut_nue_cc,
		                                                           h_multiplicity_track_flash_vtx_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_flash_vtx_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_flash_vtx_cut_numu_cc,
		                                                           h_multiplicity_track_flash_vtx_cut_nc,
		                                                           h_multiplicity_track_flash_vtx_cut_nc_pi0,
		                                                           h_multiplicity_track_flash_vtx_cut_cosmic,
		                                                           h_multiplicity_track_flash_vtx_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_flash_vtx_cut_other_mixed,
		                                                           h_multiplicity_track_flash_vtx_cut_unmatched);

		_functions_instance.selection_functions::PostCutsVtxFlash(largest_flash_v, tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                          h_vtx_flash_nue_cc_after, h_vtx_flash_nue_cc_out_fv_after, h_vtx_flash_nue_cc_mixed_after,
		                                                          h_vtx_flash_numu_cc_after, h_vtx_flash_nc_after,
		                                                          h_vtx_flash_cosmic_after, h_vtx_flash_nc_pi0_after,
		                                                          h_vtx_flash_numu_cc_mixed_after, h_vtx_flash_other_mixed_after,
		                                                          h_vtx_flash_unmatched_after);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_vtx_flash->Fill(mc_ele_energy); }
		//******************************************************
		//*** distance between pfp shower and nue object cut ***
		//******************************************************
		_functions_instance.selection_functions::PostCutsShwrVtx(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                         h_shwr_vtx_dist_nue_cc,
		                                                         h_shwr_vtx_dist_nue_cc_out_fv,
		                                                         h_shwr_vtx_dist_nue_cc_mixed,
		                                                         h_shwr_vtx_dist_numu_cc,
		                                                         h_shwr_vtx_dist_nc,
		                                                         h_shwr_vtx_dist_cosmic,
		                                                         h_shwr_vtx_dist_nc_pi0,
		                                                         h_shwr_vtx_dist_numu_cc_mixed,
		                                                         h_shwr_vtx_dist_other_mixed,
		                                                         h_shwr_vtx_dist_unmatched);

		_cuts_instance.selection_cuts::VtxNuDistance(tpc_object_container_v, shwr_nue_tolerance, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, shwr_tpco_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins,  mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_shwr_vtx, h_selected_ele_energy_shwr_vtx);
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_4, h_reco_ele_e_4, h_mc_reco_ele_e_4);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_shwr_vtx_cut_nue_cc,
		                                                                 h_ele_momentum_shwr_vtx_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_shwr_vtx_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_shwr_vtx_cut_numu_cc,
		                                                                 h_ele_momentum_shwr_vtx_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_shwr_vtx_cut_nc,
		                                                                 h_ele_momentum_shwr_vtx_cut_nc_pi0,
		                                                                 h_ele_momentum_shwr_vtx_cut_cosmic,
		                                                                 h_ele_momentum_shwr_vtx_cut_other_mixed,
		                                                                 h_ele_momentum_shwr_vtx_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_shwr_vtx_electron,
		                                                                             h_leading_momentum_shwr_vtx_photon,
		                                                                             h_leading_momentum_shwr_vtx_proton,
		                                                                             h_leading_momentum_shwr_vtx_pion,
		                                                                             h_leading_momentum_shwr_vtx_muon,
		                                                                             h_leading_momentum_shwr_vtx_kaon,
		                                                                             h_leading_momentum_shwr_vtx_neutron,
		                                                                             h_leading_momentum_shwr_vtx_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_shwr_vtx_cut_nue_cc,
		                                                           h_multiplicity_shower_shwr_vtx_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_shwr_vtx_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_shwr_vtx_cut_numu_cc,
		                                                           h_multiplicity_shower_shwr_vtx_cut_nc,
		                                                           h_multiplicity_shower_shwr_vtx_cut_nc_pi0,
		                                                           h_multiplicity_shower_shwr_vtx_cut_cosmic,
		                                                           h_multiplicity_shower_shwr_vtx_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_shwr_vtx_cut_other_mixed,
		                                                           h_multiplicity_shower_shwr_vtx_cut_unmatched,
		                                                           h_multiplicity_track_shwr_vtx_cut_nue_cc,
		                                                           h_multiplicity_track_shwr_vtx_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_shwr_vtx_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_shwr_vtx_cut_numu_cc,
		                                                           h_multiplicity_track_shwr_vtx_cut_nc,
		                                                           h_multiplicity_track_shwr_vtx_cut_nc_pi0,
		                                                           h_multiplicity_track_shwr_vtx_cut_cosmic,
		                                                           h_multiplicity_track_shwr_vtx_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_shwr_vtx_cut_other_mixed,
		                                                           h_multiplicity_track_shwr_vtx_cut_unmatched);

		_functions_instance.selection_functions::PostCutsShwrVtx(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                         h_shwr_vtx_dist_nue_cc_after,
		                                                         h_shwr_vtx_dist_nue_cc_out_fv_after,
		                                                         h_shwr_vtx_dist_nue_cc_mixed_after,
		                                                         h_shwr_vtx_dist_numu_cc_after,
		                                                         h_shwr_vtx_dist_nc_after,
		                                                         h_shwr_vtx_dist_cosmic_after,
		                                                         h_shwr_vtx_dist_nc_pi0_after,
		                                                         h_shwr_vtx_dist_numu_cc_mixed_after,
		                                                         h_shwr_vtx_dist_other_mixed_after,
		                                                         h_shwr_vtx_dist_unmatched_after);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_shwr_vtx->Fill(mc_ele_energy); }
		//******************************************************
		// **** distance between pfp track and nue object cut **
		//******************************************************
		_functions_instance.selection_functions::PostCutTrkVtx(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                       h_trk_vtx_dist_nue_cc, h_trk_vtx_dist_nue_cc_out_fv,
		                                                       h_trk_vtx_dist_nue_cc_mixed,
		                                                       h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_nc,
		                                                       h_trk_vtx_dist_cosmic, h_trk_vtx_dist_nc_pi0,
		                                                       h_trk_vtx_dist_numu_cc_mixed, h_trk_vtx_dist_other_mixed,
		                                                       h_trk_vtx_dist_unmatched);
		_cuts_instance.selection_cuts::VtxTrackNuDistance(tpc_object_container_v, trk_nue_tolerance, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, trk_tpco_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_trk_vtx, h_selected_ele_energy_trk_vtx);
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_5, h_reco_ele_e_5, h_mc_reco_ele_e_5);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_trk_vtx_cut_nue_cc,
		                                                                 h_ele_momentum_trk_vtx_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_trk_vtx_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_trk_vtx_cut_numu_cc,
		                                                                 h_ele_momentum_trk_vtx_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_trk_vtx_cut_nc,
		                                                                 h_ele_momentum_trk_vtx_cut_nc_pi0,
		                                                                 h_ele_momentum_trk_vtx_cut_cosmic,
		                                                                 h_ele_momentum_trk_vtx_cut_other_mixed,
		                                                                 h_ele_momentum_trk_vtx_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_trk_vtx_electron,
		                                                                             h_leading_momentum_trk_vtx_photon,
		                                                                             h_leading_momentum_trk_vtx_proton,
		                                                                             h_leading_momentum_trk_vtx_pion,
		                                                                             h_leading_momentum_trk_vtx_muon,
		                                                                             h_leading_momentum_trk_vtx_kaon,
		                                                                             h_leading_momentum_trk_vtx_neutron,
		                                                                             h_leading_momentum_trk_vtx_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_trk_vtx_cut_nue_cc,
		                                                           h_multiplicity_shower_trk_vtx_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_trk_vtx_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_trk_vtx_cut_numu_cc,
		                                                           h_multiplicity_shower_trk_vtx_cut_nc,
		                                                           h_multiplicity_shower_trk_vtx_cut_nc_pi0,
		                                                           h_multiplicity_shower_trk_vtx_cut_cosmic,
		                                                           h_multiplicity_shower_trk_vtx_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_trk_vtx_cut_other_mixed,
		                                                           h_multiplicity_shower_trk_vtx_cut_unmatched,
		                                                           h_multiplicity_track_trk_vtx_cut_nue_cc,
		                                                           h_multiplicity_track_trk_vtx_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_trk_vtx_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_trk_vtx_cut_numu_cc,
		                                                           h_multiplicity_track_trk_vtx_cut_nc,
		                                                           h_multiplicity_track_trk_vtx_cut_nc_pi0,
		                                                           h_multiplicity_track_trk_vtx_cut_cosmic,
		                                                           h_multiplicity_track_trk_vtx_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_trk_vtx_cut_other_mixed,
		                                                           h_multiplicity_track_trk_vtx_cut_unmatched);

		_functions_instance.selection_functions::PostCutTrkVtx(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                       h_trk_vtx_dist_nue_cc_after, h_trk_vtx_dist_nue_cc_out_fv_after,
		                                                       h_trk_vtx_dist_nue_cc_mixed_after,
		                                                       h_trk_vtx_dist_numu_cc_after, h_trk_vtx_dist_nc_after,
		                                                       h_trk_vtx_dist_cosmic_after, h_trk_vtx_dist_nc_pi0_after,
		                                                       h_trk_vtx_dist_numu_cc_mixed_after, h_trk_vtx_dist_other_mixed_after,
		                                                       h_trk_vtx_dist_unmatched_after);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_trk_vtx->Fill(mc_ele_energy); }
		//****************************************************
		// ******** hit threshold for showers cut *************
		//******************************************************
		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins->at(0) >= 1 && true_in_tpc == true)
		{
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose,
			                                                             tpco_classifier_v, mc_nu_energy, mc_ele_energy,
			                                                             h_shwr_hits_nu_eng_last, h_shwr_hits_ele_eng_last,
			                                                             h_shwr_collection_hits_nu_eng_last, h_shwr_collection_hits_ele_eng_last);
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose,
			                                                             tpco_classifier_v, mc_nu_energy, mc_ele_energy,
			                                                             h_shwr_hits_nu_eng_zoom_last, h_shwr_hits_ele_eng_zoom_last,
			                                                             h_shwr_collection_hits_nu_eng_zoom_last, h_shwr_collection_hits_ele_eng_zoom_last);
		}

		_functions_instance.selection_functions::HitsPlots1D(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                     h_pre_cut_collection_hits_track_nue_cc,
		                                                     h_pre_cut_collection_hits_track_nue_cc_out_fv,
		                                                     h_pre_cut_collection_hits_track_nue_cc_mixed,
		                                                     h_pre_cut_collection_hits_track_numu_cc,
		                                                     h_pre_cut_collection_hits_track_numu_cc_mixed,
		                                                     h_pre_cut_collection_hits_track_nc,
		                                                     h_pre_cut_collection_hits_track_nc_pi0,
		                                                     h_pre_cut_collection_hits_track_cosmic,
		                                                     h_pre_cut_collection_hits_track_other_mixed,
		                                                     h_pre_cut_collection_hits_track_unmatched,
		                                                     h_pre_cut_collection_hits_shower_nue_cc,
		                                                     h_pre_cut_collection_hits_shower_nue_cc_out_fv,
		                                                     h_pre_cut_collection_hits_shower_nue_cc_mixed,
		                                                     h_pre_cut_collection_hits_shower_numu_cc,
		                                                     h_pre_cut_collection_hits_shower_numu_cc_mixed,
		                                                     h_pre_cut_collection_hits_shower_nc,
		                                                     h_pre_cut_collection_hits_shower_nc_pi0,
		                                                     h_pre_cut_collection_hits_shower_cosmic,
		                                                     h_pre_cut_collection_hits_shower_other_mixed,
		                                                     h_pre_cut_collection_hits_shower_unmatched,
		                                                     h_pre_cut_collection_hits_leading_shower_nue_cc,
		                                                     h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv,
		                                                     h_pre_cut_collection_hits_leading_shower_nue_cc_mixed,
		                                                     h_pre_cut_collection_hits_leading_shower_numu_cc,
		                                                     h_pre_cut_collection_hits_leading_shower_numu_cc_mixed,
		                                                     h_pre_cut_collection_hits_leading_shower_nc,
		                                                     h_pre_cut_collection_hits_leading_shower_nc_pi0,
		                                                     h_pre_cut_collection_hits_leading_shower_cosmic,
		                                                     h_pre_cut_collection_hits_leading_shower_other_mixed,
		                                                     h_pre_cut_collection_hits_leading_shower_unmatched,
		                                                     h_pre_cut_total_hits_leading_shower_nue_cc,
		                                                     h_pre_cut_total_hits_leading_shower_nue_cc_out_fv,
		                                                     h_pre_cut_total_hits_leading_shower_nue_cc_mixed,
		                                                     h_pre_cut_total_hits_leading_shower_numu_cc,
		                                                     h_pre_cut_total_hits_leading_shower_numu_cc_mixed,
		                                                     h_pre_cut_total_hits_leading_shower_nc,
		                                                     h_pre_cut_total_hits_leading_shower_nc_pi0,
		                                                     h_pre_cut_total_hits_leading_shower_cosmic,
		                                                     h_pre_cut_total_hits_leading_shower_other_mixed,
		                                                     h_pre_cut_total_hits_leading_shower_unmatched);

		_cuts_instance.selection_cuts::HitThreshold(tpc_object_container_v, shwr_hit_threshold, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_threshold_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_hit_threshold, h_selected_ele_energy_hit_threshold);
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_6, h_reco_ele_e_6, h_mc_reco_ele_e_6);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_hit_cut_nue_cc,
		                                                                 h_ele_momentum_hit_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_hit_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_hit_cut_numu_cc,
		                                                                 h_ele_momentum_hit_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_hit_cut_nc,
		                                                                 h_ele_momentum_hit_cut_nc_pi0,
		                                                                 h_ele_momentum_hit_cut_cosmic,
		                                                                 h_ele_momentum_hit_cut_other_mixed,
		                                                                 h_ele_momentum_hit_cut_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_hit_cut_nue_cc,
		                                                           h_multiplicity_shower_hit_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_hit_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_hit_cut_numu_cc,
		                                                           h_multiplicity_shower_hit_cut_nc,
		                                                           h_multiplicity_shower_hit_cut_nc_pi0,
		                                                           h_multiplicity_shower_hit_cut_cosmic,
		                                                           h_multiplicity_shower_hit_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_hit_cut_other_mixed,
		                                                           h_multiplicity_shower_hit_cut_unmatched,
		                                                           h_multiplicity_track_hit_cut_nue_cc,
		                                                           h_multiplicity_track_hit_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_hit_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_hit_cut_numu_cc,
		                                                           h_multiplicity_track_hit_cut_nc,
		                                                           h_multiplicity_track_hit_cut_nc_pi0,
		                                                           h_multiplicity_track_hit_cut_cosmic,
		                                                           h_multiplicity_track_hit_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_hit_cut_other_mixed,
		                                                           h_multiplicity_track_hit_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_hit_cut_electron,
		                                                                             h_leading_momentum_hit_cut_photon,
		                                                                             h_leading_momentum_hit_cut_proton,
		                                                                             h_leading_momentum_hit_cut_pion,
		                                                                             h_leading_momentum_hit_cut_muon,
		                                                                             h_leading_momentum_hit_cut_kaon,
		                                                                             h_leading_momentum_hit_cut_neutron,
		                                                                             h_leading_momentum_hit_cut_mc_unmatched);

		_functions_instance.selection_functions::PlaneHitsComparisonTrack(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                  h_collection_total_hits_track_nue_cc,
		                                                                  h_collection_total_hits_track_nue_cc_out_fv,
		                                                                  h_collection_total_hits_track_nue_cc_mixed,
		                                                                  h_collection_total_hits_track_numu_cc,
		                                                                  h_collection_total_hits_track_numu_cc_mixed,
		                                                                  h_collection_total_hits_track_nc,
		                                                                  h_collection_total_hits_track_nc_pi0,
		                                                                  h_collection_total_hits_track_cosmic,
		                                                                  h_collection_total_hits_track_other_mixed,
		                                                                  h_collection_total_hits_track_unmatched);
		_functions_instance.selection_functions::PlaneHitsComparisonShower(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                   h_collection_total_hits_shower_nue_cc,
		                                                                   h_collection_total_hits_shower_nue_cc_out_fv,
		                                                                   h_collection_total_hits_shower_nue_cc_mixed,
		                                                                   h_collection_total_hits_shower_numu_cc,
		                                                                   h_collection_total_hits_shower_numu_cc_mixed,
		                                                                   h_collection_total_hits_shower_nc,
		                                                                   h_collection_total_hits_shower_nc_pi0,
		                                                                   h_collection_total_hits_shower_cosmic,
		                                                                   h_collection_total_hits_shower_other_mixed,
		                                                                   h_collection_total_hits_shower_unmatched);
		_functions_instance.selection_functions::PlaneHitsComparisonLeadingShower(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                          h_collection_total_hits_leading_shower_nue_cc,
		                                                                          h_collection_total_hits_leading_shower_nue_cc_out_fv,
		                                                                          h_collection_total_hits_leading_shower_nue_cc_mixed,
		                                                                          h_collection_total_hits_leading_shower_numu_cc,
		                                                                          h_collection_total_hits_leading_shower_numu_cc_mixed,
		                                                                          h_collection_total_hits_leading_shower_nc,
		                                                                          h_collection_total_hits_leading_shower_nc_pi0,
		                                                                          h_collection_total_hits_leading_shower_cosmic,
		                                                                          h_collection_total_hits_leading_shower_other_mixed,
		                                                                          h_collection_total_hits_leading_shower_unmatched);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_hit->Fill(mc_ele_energy); }

		//***************************************//
		//*** Collection Plane Hits Threshold ***//
		//***************************************//
		_functions_instance.selection_functions::LeadingPhi(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                    h_ele_pfp_phi_nue_cc,
		                                                    h_ele_pfp_phi_nue_cc_out_fv,
		                                                    h_ele_pfp_phi_nue_cc_mixed,
		                                                    h_ele_pfp_phi_numu_cc,
		                                                    h_ele_pfp_phi_numu_cc_mixed,
		                                                    h_ele_pfp_phi_nc,
		                                                    h_ele_pfp_phi_nc_pi0,
		                                                    h_ele_pfp_phi_cosmic,
		                                                    h_ele_pfp_phi_other_mixed,
		                                                    h_ele_pfp_phi_unmatched);

		_functions_instance.selection_functions::LeadingTheta(tpc_object_container_v, passed_tpco, 0, 0, _verbose, tpco_classifier_v,
		                                                      h_ele_pfp_theta_nue_cc,
		                                                      h_ele_pfp_theta_nue_cc_out_fv,
		                                                      h_ele_pfp_theta_nue_cc_mixed,
		                                                      h_ele_pfp_theta_numu_cc,
		                                                      h_ele_pfp_theta_numu_cc_mixed,
		                                                      h_ele_pfp_theta_nc,
		                                                      h_ele_pfp_theta_nc_pi0,
		                                                      h_ele_pfp_theta_cosmic,
		                                                      h_ele_pfp_theta_other_mixed,
		                                                      h_ele_pfp_theta_unmatched);

		_functions_instance.selection_functions::HitsPlots1D(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                     h_collection_hits_track_nue_cc,
		                                                     h_collection_hits_track_nue_cc_out_fv,
		                                                     h_collection_hits_track_nue_cc_mixed,
		                                                     h_collection_hits_track_numu_cc,
		                                                     h_collection_hits_track_numu_cc_mixed,
		                                                     h_collection_hits_track_nc,
		                                                     h_collection_hits_track_nc_pi0,
		                                                     h_collection_hits_track_cosmic,
		                                                     h_collection_hits_track_other_mixed,
		                                                     h_collection_hits_track_unmatched,
		                                                     h_collection_hits_shower_nue_cc,
		                                                     h_collection_hits_shower_nue_cc_out_fv,
		                                                     h_collection_hits_shower_nue_cc_mixed,
		                                                     h_collection_hits_shower_numu_cc,
		                                                     h_collection_hits_shower_numu_cc_mixed,
		                                                     h_collection_hits_shower_nc,
		                                                     h_collection_hits_shower_nc_pi0,
		                                                     h_collection_hits_shower_cosmic,
		                                                     h_collection_hits_shower_other_mixed,
		                                                     h_collection_hits_shower_unmatched,
		                                                     h_collection_hits_leading_shower_nue_cc,
		                                                     h_collection_hits_leading_shower_nue_cc_out_fv,
		                                                     h_collection_hits_leading_shower_nue_cc_mixed,
		                                                     h_collection_hits_leading_shower_numu_cc,
		                                                     h_collection_hits_leading_shower_numu_cc_mixed,
		                                                     h_collection_hits_leading_shower_nc,
		                                                     h_collection_hits_leading_shower_nc_pi0,
		                                                     h_collection_hits_leading_shower_cosmic,
		                                                     h_collection_hits_leading_shower_other_mixed,
		                                                     h_collection_hits_leading_shower_unmatched,
		                                                     h_total_hits_leading_shower_nue_cc,
		                                                     h_total_hits_leading_shower_nue_cc_out_fv,
		                                                     h_total_hits_leading_shower_nue_cc_mixed,
		                                                     h_total_hits_leading_shower_numu_cc,
		                                                     h_total_hits_leading_shower_numu_cc_mixed,
		                                                     h_total_hits_leading_shower_nc,
		                                                     h_total_hits_leading_shower_nc_pi0,
		                                                     h_total_hits_leading_shower_cosmic,
		                                                     h_total_hits_leading_shower_other_mixed,
		                                                     h_total_hits_leading_shower_unmatched);

		_cuts_instance.selection_cuts::HitThresholdCollection(tpc_object_container_v, shwr_hit_threshold_collection, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_threshold_collection_counter_v);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_7, h_reco_ele_e_7, h_mc_reco_ele_e_7);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_yhit_cut_nue_cc,
		                                                                 h_ele_momentum_yhit_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_yhit_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_yhit_cut_numu_cc,
		                                                                 h_ele_momentum_yhit_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_yhit_cut_nc,
		                                                                 h_ele_momentum_yhit_cut_nc_pi0,
		                                                                 h_ele_momentum_yhit_cut_cosmic,
		                                                                 h_ele_momentum_yhit_cut_other_mixed,
		                                                                 h_ele_momentum_yhit_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_yhit_cut_electron,
		                                                                             h_leading_momentum_yhit_cut_photon,
		                                                                             h_leading_momentum_yhit_cut_proton,
		                                                                             h_leading_momentum_yhit_cut_pion,
		                                                                             h_leading_momentum_yhit_cut_muon,
		                                                                             h_leading_momentum_yhit_cut_kaon,
		                                                                             h_leading_momentum_yhit_cut_neutron,
		                                                                             h_leading_momentum_yhit_cut_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_yhit_cut_nue_cc,
		                                                           h_multiplicity_shower_yhit_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_yhit_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_yhit_cut_numu_cc,
		                                                           h_multiplicity_shower_yhit_cut_nc,
		                                                           h_multiplicity_shower_yhit_cut_nc_pi0,
		                                                           h_multiplicity_shower_yhit_cut_cosmic,
		                                                           h_multiplicity_shower_yhit_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_yhit_cut_other_mixed,
		                                                           h_multiplicity_shower_yhit_cut_unmatched,
		                                                           h_multiplicity_track_yhit_cut_nue_cc,
		                                                           h_multiplicity_track_yhit_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_yhit_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_yhit_cut_numu_cc,
		                                                           h_multiplicity_track_yhit_cut_nc,
		                                                           h_multiplicity_track_yhit_cut_nc_pi0,
		                                                           h_multiplicity_track_yhit_cut_cosmic,
		                                                           h_multiplicity_track_yhit_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_yhit_cut_other_mixed,
		                                                           h_multiplicity_track_yhit_cut_unmatched);

		_functions_instance.selection_functions::HitsPlots1D(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                     h_collection_hits_track_nue_cc_after,
		                                                     h_collection_hits_track_nue_cc_out_fv_after,
		                                                     h_collection_hits_track_nue_cc_mixed_after,
		                                                     h_collection_hits_track_numu_cc_after,
		                                                     h_collection_hits_track_numu_cc_mixed_after,
		                                                     h_collection_hits_track_nc_after,
		                                                     h_collection_hits_track_nc_pi0_after,
		                                                     h_collection_hits_track_cosmic_after,
		                                                     h_collection_hits_track_other_mixed_after,
		                                                     h_collection_hits_track_unmatched_after,
		                                                     h_collection_hits_shower_nue_cc_after,
		                                                     h_collection_hits_shower_nue_cc_out_fv_after,
		                                                     h_collection_hits_shower_nue_cc_mixed_after,
		                                                     h_collection_hits_shower_numu_cc_after,
		                                                     h_collection_hits_shower_numu_cc_mixed_after,
		                                                     h_collection_hits_shower_nc_after,
		                                                     h_collection_hits_shower_nc_pi0_after,
		                                                     h_collection_hits_shower_cosmic_after,
		                                                     h_collection_hits_shower_other_mixed_after,
		                                                     h_collection_hits_shower_unmatched_after,
		                                                     h_collection_hits_leading_shower_nue_cc_after,
		                                                     h_collection_hits_leading_shower_nue_cc_out_fv_after,
		                                                     h_collection_hits_leading_shower_nue_cc_mixed_after,
		                                                     h_collection_hits_leading_shower_numu_cc_after,
		                                                     h_collection_hits_leading_shower_numu_cc_mixed_after,
		                                                     h_collection_hits_leading_shower_nc_after,
		                                                     h_collection_hits_leading_shower_nc_pi0_after,
		                                                     h_collection_hits_leading_shower_cosmic_after,
		                                                     h_collection_hits_leading_shower_other_mixed_after,
		                                                     h_collection_hits_leading_shower_unmatched_after,
		                                                     h_total_hits_leading_shower_nue_cc_after,
		                                                     h_total_hits_leading_shower_nue_cc_out_fv_after,
		                                                     h_total_hits_leading_shower_nue_cc_mixed_after,
		                                                     h_total_hits_leading_shower_numu_cc_after,
		                                                     h_total_hits_leading_shower_numu_cc_mixed_after,
		                                                     h_total_hits_leading_shower_nc_after,
		                                                     h_total_hits_leading_shower_nc_pi0_after,
		                                                     h_total_hits_leading_shower_cosmic_after,
		                                                     h_total_hits_leading_shower_other_mixed_after,
		                                                     h_total_hits_leading_shower_unmatched_after);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_yhit->Fill(mc_ele_energy); }
		//*****************************************************
		//****** open angle cut for the leading shower ********
		//******************************************************
		_functions_instance.selection_functions::dEdxVsOpenAngle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                         h_dedx_open_angle_nue_cc, h_dedx_open_angle_nue_cc_out_fv,
		                                                         h_dedx_open_angle_nue_cc_mixed, h_dedx_open_angle_numu_cc,
		                                                         h_dedx_open_angle_numu_cc_mixed, h_dedx_open_angle_nc,
		                                                         h_dedx_open_angle_nc_pi0, h_dedx_open_angle_cosmic,
		                                                         h_dedx_open_angle_other_mixed, h_dedx_open_angle_unmatched);

		_functions_instance.selection_functions::LeadingCosTheta(tpc_object_container_v, passed_tpco, 0, 0, _verbose, tpco_classifier_v,
		                                                         h_ele_cos_theta_nue_cc,
		                                                         h_ele_cos_theta_nue_cc_out_fv,
		                                                         h_ele_cos_theta_nue_cc_mixed,
		                                                         h_ele_cos_theta_numu_cc,
		                                                         h_ele_cos_theta_numu_cc_mixed,
		                                                         h_ele_cos_theta_nc,
		                                                         h_ele_cos_theta_nc_pi0,
		                                                         h_ele_cos_theta_cosmic,
		                                                         h_ele_cos_theta_other_mixed,
		                                                         h_ele_cos_theta_unmatched);

		_functions_instance.selection_functions::FillPostCutVector(tpc_object_container_v, passed_tpco, tpco_classifier_v, post_open_angle_cuts_v);
		_functions_instance.selection_functions::NumShowersOpenAngle(tpc_object_container_v, passed_tpco, tpco_classifier_v,
		                                                             h_pfp_shower_open_angle_nue_cc_qe,
		                                                             h_pfp_shower_open_angle_nue_cc_out_fv,
		                                                             h_pfp_shower_open_angle_nue_cc_res,
		                                                             h_pfp_shower_open_angle_nue_cc_dis,
		                                                             h_pfp_shower_open_angle_nue_cc_coh,
		                                                             h_pfp_shower_open_angle_nue_cc_mec,
		                                                             h_pfp_shower_open_angle_nc,
		                                                             h_pfp_shower_open_angle_numu_cc_qe,
		                                                             h_pfp_shower_open_angle_numu_cc_res,
		                                                             h_pfp_shower_open_angle_numu_cc_dis,
		                                                             h_pfp_shower_open_angle_numu_cc_coh,
		                                                             h_pfp_shower_open_angle_numu_cc_mec,
		                                                             h_pfp_shower_open_angle_nc_pi0,
		                                                             h_pfp_shower_open_angle_nue_cc_mixed,
		                                                             h_pfp_shower_open_angle_numu_cc_mixed,
		                                                             h_pfp_shower_open_angle_cosmic,
		                                                             h_pfp_shower_open_angle_other_mixed,
		                                                             h_pfp_shower_open_angle_unmatched
		                                                             );
		_functions_instance.selection_functions::PostCutOpenAngle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                          h_leading_shower_open_angle_nue_cc, h_leading_shower_open_angle_nue_cc_out_fv,
		                                                          h_leading_shower_open_angle_nue_cc_mixed,
		                                                          h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_nc,
		                                                          h_leading_shower_open_angle_cosmic, h_leading_shower_open_angle_nc_pi0,
		                                                          h_leading_shower_open_angle_numu_cc_mixed, h_leading_shower_open_angle_other_mixed,
		                                                          h_leading_shower_open_angle_unmatched);
		_functions_instance.selection_functions::PostCutOpenAngle1Shower(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_leading_shower_open_angle_1_nue_cc, h_leading_shower_open_angle_1_nue_cc_out_fv,
		                                                                 h_leading_shower_open_angle_1_nue_cc_mixed,
		                                                                 h_leading_shower_open_angle_1_numu_cc, h_leading_shower_open_angle_1_nc,
		                                                                 h_leading_shower_open_angle_1_cosmic, h_leading_shower_open_angle_1_nc_pi0,
		                                                                 h_leading_shower_open_angle_1_numu_cc_mixed, h_leading_shower_open_angle_1_other_mixed,
		                                                                 h_leading_shower_open_angle_1_unmatched);
		_functions_instance.selection_functions::PostCutOpenAngle2PlusShower(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                     h_leading_shower_open_angle_2plus_nue_cc,
		                                                                     h_leading_shower_open_angle_2plus_nue_cc_out_fv,
		                                                                     h_leading_shower_open_angle_2plus_nue_cc_mixed,
		                                                                     h_leading_shower_open_angle_2plus_numu_cc, h_leading_shower_open_angle_2plus_nc,
		                                                                     h_leading_shower_open_angle_2plus_cosmic, h_leading_shower_open_angle_2plus_nc_pi0,
		                                                                     h_leading_shower_open_angle_2plus_numu_cc_mixed,
		                                                                     h_leading_shower_open_angle_2plus_other_mixed,
		                                                                     h_leading_shower_open_angle_2plus_unmatched);
		_cuts_instance.selection_cuts::OpenAngleCut(tpc_object_container_v, passed_tpco, tolerance_open_angle, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, open_angle_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_open_angle, h_selected_ele_energy_open_angle);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_8, h_reco_ele_e_8, h_mc_reco_ele_e_8);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_open_angle_cut_nue_cc,
		                                                                 h_ele_momentum_open_angle_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_open_angle_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_open_angle_cut_numu_cc,
		                                                                 h_ele_momentum_open_angle_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_open_angle_cut_nc,
		                                                                 h_ele_momentum_open_angle_cut_nc_pi0,
		                                                                 h_ele_momentum_open_angle_cut_cosmic,
		                                                                 h_ele_momentum_open_angle_cut_other_mixed,
		                                                                 h_ele_momentum_open_angle_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_open_angle_cut_electron,
		                                                                             h_leading_momentum_open_angle_cut_photon,
		                                                                             h_leading_momentum_open_angle_cut_proton,
		                                                                             h_leading_momentum_open_angle_cut_pion,
		                                                                             h_leading_momentum_open_angle_cut_muon,
		                                                                             h_leading_momentum_open_angle_cut_kaon,
		                                                                             h_leading_momentum_open_angle_cut_neutron,
		                                                                             h_leading_momentum_open_angle_cut_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_open_angle_cut_nue_cc,
		                                                           h_multiplicity_shower_open_angle_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_open_angle_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_open_angle_cut_numu_cc,
		                                                           h_multiplicity_shower_open_angle_cut_nc,
		                                                           h_multiplicity_shower_open_angle_cut_nc_pi0,
		                                                           h_multiplicity_shower_open_angle_cut_cosmic,
		                                                           h_multiplicity_shower_open_angle_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_open_angle_cut_other_mixed,
		                                                           h_multiplicity_shower_open_angle_cut_unmatched,
		                                                           h_multiplicity_track_open_angle_cut_nue_cc,
		                                                           h_multiplicity_track_open_angle_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_open_angle_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_open_angle_cut_numu_cc,
		                                                           h_multiplicity_track_open_angle_cut_nc,
		                                                           h_multiplicity_track_open_angle_cut_nc_pi0,
		                                                           h_multiplicity_track_open_angle_cut_cosmic,
		                                                           h_multiplicity_track_open_angle_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_open_angle_cut_other_mixed,
		                                                           h_multiplicity_track_open_angle_cut_unmatched);

		_functions_instance.selection_functions::NumShowersOpenAngle(tpc_object_container_v, passed_tpco, tpco_classifier_v,
		                                                             h_pfp_shower_dedx_nue_cc_qe,
		                                                             h_pfp_shower_dedx_nue_cc_out_fv,
		                                                             h_pfp_shower_dedx_nue_cc_res,
		                                                             h_pfp_shower_dedx_nue_cc_dis,
		                                                             h_pfp_shower_dedx_nue_cc_coh,
		                                                             h_pfp_shower_dedx_nue_cc_mec,
		                                                             h_pfp_shower_dedx_nc,
		                                                             h_pfp_shower_dedx_numu_cc_qe,
		                                                             h_pfp_shower_dedx_numu_cc_res,
		                                                             h_pfp_shower_dedx_numu_cc_dis,
		                                                             h_pfp_shower_dedx_numu_cc_coh,
		                                                             h_pfp_shower_dedx_numu_cc_mec,
		                                                             h_pfp_shower_dedx_nc_pi0,
		                                                             h_pfp_shower_dedx_nue_cc_mixed,
		                                                             h_pfp_shower_dedx_numu_cc_mixed,
		                                                             h_pfp_shower_dedx_cosmic,
		                                                             h_pfp_shower_dedx_other_mixed,
		                                                             h_pfp_shower_dedx_unmatched);

		_functions_instance.selection_functions::PostCutOpenAngle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                          h_leading_shower_open_angle_nue_cc_after, h_leading_shower_open_angle_nue_cc_out_fv_after,
		                                                          h_leading_shower_open_angle_nue_cc_mixed_after,
		                                                          h_leading_shower_open_angle_numu_cc_after, h_leading_shower_open_angle_nc_after,
		                                                          h_leading_shower_open_angle_cosmic_after, h_leading_shower_open_angle_nc_pi0_after,
		                                                          h_leading_shower_open_angle_numu_cc_mixed_after, h_leading_shower_open_angle_other_mixed_after,
		                                                          h_leading_shower_open_angle_unmatched_after);


		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_open_angle->Fill(mc_ele_energy); }
		//*****************************************************
		//*********** dEdx cut for the leading shower *********
		//******************************************************
		_functions_instance.selection_functions::PostCutsdEdx(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                      h_dedx_cuts_nue_cc,        h_dedx_cuts_nue_cc_mixed,
		                                                      h_dedx_cuts_nue_cc_out_fv,
		                                                      h_dedx_cuts_numu_cc,       h_dedx_cuts_nc,
		                                                      h_dedx_cuts_cosmic,        h_dedx_cuts_nc_pi0,
		                                                      h_dedx_cuts_numu_cc_mixed, h_dedx_cuts_other_mixed,
		                                                      h_dedx_cuts_unmatched);

		_functions_instance.selection_functions::dEdxCollectionAngle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                             h_dedx_collection_angle_nue_cc,        h_dedx_collection_angle_nue_cc_mixed,
		                                                             h_dedx_collection_angle_nue_cc_out_fv,
		                                                             h_dedx_collection_angle_numu_cc,       h_dedx_collection_angle_nc,
		                                                             h_dedx_collection_angle_cosmic,        h_dedx_collection_angle_nc_pi0,
		                                                             h_dedx_collection_angle_numu_cc_mixed, h_dedx_collection_angle_other_mixed,
		                                                             h_dedx_collection_angle_unmatched);

		_functions_instance.selection_functions::PostCutsdEdxTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                  h_dedx_cuts_electron, h_dedx_cuts_photon, h_dedx_cuts_proton, h_dedx_cuts_pion,
		                                                                  h_dedx_cuts_muon, h_dedx_cuts_kaon, h_dedx_cuts_neutron, h_dedx_cuts_mc_unmatched);
		_functions_instance.selection_functions::PostCutsdEdxHitsTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                      h_dedx_cuts_hits_electron, h_dedx_cuts_hits_photon, h_dedx_cuts_hits_proton,
		                                                                      h_dedx_cuts_hits_pion, h_dedx_cuts_hits_muon, h_dedx_cuts_hits_kaon,
		                                                                      h_dedx_cuts_hits_neutron, h_dedx_cuts_hits_mc_unmatched,
		                                                                      h_dedx_cuts_collection_hits_electron, h_dedx_cuts_collection_hits_photon,
		                                                                      h_dedx_cuts_collection_hits_proton, h_dedx_cuts_collection_hits_pion,
		                                                                      h_dedx_cuts_collection_hits_muon, h_dedx_cuts_collection_hits_kaon,
		                                                                      h_dedx_cuts_collection_hits_neutron, h_dedx_cuts_collection_hits_mc_unmatched);

		_functions_instance.selection_functions::dEdxTheta(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                   h_dedx_theta_pre_cuts_nue_cc,        h_dedx_theta_pre_cuts_nue_cc_mixed,
		                                                   h_dedx_theta_pre_cuts_nue_cc_out_fv,
		                                                   h_dedx_theta_pre_cuts_numu_cc,       h_dedx_theta_pre_cuts_nc,
		                                                   h_dedx_theta_pre_cuts_cosmic,        h_dedx_theta_pre_cuts_nc_pi0,
		                                                   h_dedx_theta_pre_cuts_numu_cc_mixed, h_dedx_theta_pre_cuts_other_mixed,
		                                                   h_dedx_theta_pre_cuts_unmatched);

		_functions_instance.selection_functions::PostCutsdedxThetaSlice(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                h_dedx_slice_1_nue_cc,
		                                                                h_dedx_slice_1_nue_cc_out_fv,
		                                                                h_dedx_slice_1_nue_cc_mixed,
		                                                                h_dedx_slice_1_numu_cc,
		                                                                h_dedx_slice_1_numu_cc_mixed,
		                                                                h_dedx_slice_1_nc,
		                                                                h_dedx_slice_1_nc_pi0,
		                                                                h_dedx_slice_1_cosmic,
		                                                                h_dedx_slice_1_other_mixed,
		                                                                h_dedx_slice_1_unmatched,
		                                                                h_dedx_slice_2_nue_cc,
		                                                                h_dedx_slice_2_nue_cc_out_fv,
		                                                                h_dedx_slice_2_nue_cc_mixed,
		                                                                h_dedx_slice_2_numu_cc,
		                                                                h_dedx_slice_2_numu_cc_mixed,
		                                                                h_dedx_slice_2_nc,
		                                                                h_dedx_slice_2_nc_pi0,
		                                                                h_dedx_slice_2_cosmic,
		                                                                h_dedx_slice_2_other_mixed,
		                                                                h_dedx_slice_2_unmatched,
		                                                                h_dedx_slice_3_nue_cc,
		                                                                h_dedx_slice_3_nue_cc_out_fv,
		                                                                h_dedx_slice_3_nue_cc_mixed,
		                                                                h_dedx_slice_3_numu_cc,
		                                                                h_dedx_slice_3_numu_cc_mixed,
		                                                                h_dedx_slice_3_nc,
		                                                                h_dedx_slice_3_nc_pi0,
		                                                                h_dedx_slice_3_cosmic,
		                                                                h_dedx_slice_3_other_mixed,
		                                                                h_dedx_slice_3_unmatched);

		_functions_instance.selection_functions::PostCutsdedxThetaSlice(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                h_dedx_slice_1_zoom_nue_cc,
		                                                                h_dedx_slice_1_zoom_nue_cc_out_fv,
		                                                                h_dedx_slice_1_zoom_nue_cc_mixed,
		                                                                h_dedx_slice_1_zoom_numu_cc,
		                                                                h_dedx_slice_1_zoom_numu_cc_mixed,
		                                                                h_dedx_slice_1_zoom_nc,
		                                                                h_dedx_slice_1_zoom_nc_pi0,
		                                                                h_dedx_slice_1_zoom_cosmic,
		                                                                h_dedx_slice_1_zoom_other_mixed,
		                                                                h_dedx_slice_1_zoom_unmatched,
		                                                                h_dedx_slice_2_zoom_nue_cc,
		                                                                h_dedx_slice_2_zoom_nue_cc_out_fv,
		                                                                h_dedx_slice_2_zoom_nue_cc_mixed,
		                                                                h_dedx_slice_2_zoom_numu_cc,
		                                                                h_dedx_slice_2_zoom_numu_cc_mixed,
		                                                                h_dedx_slice_2_zoom_nc,
		                                                                h_dedx_slice_2_zoom_nc_pi0,
		                                                                h_dedx_slice_2_zoom_cosmic,
		                                                                h_dedx_slice_2_zoom_other_mixed,
		                                                                h_dedx_slice_2_zoom_unmatched,
		                                                                h_dedx_slice_3_zoom_nue_cc,
		                                                                h_dedx_slice_3_zoom_nue_cc_out_fv,
		                                                                h_dedx_slice_3_zoom_nue_cc_mixed,
		                                                                h_dedx_slice_3_zoom_numu_cc,
		                                                                h_dedx_slice_3_zoom_numu_cc_mixed,
		                                                                h_dedx_slice_3_zoom_nc,
		                                                                h_dedx_slice_3_zoom_nc_pi0,
		                                                                h_dedx_slice_3_zoom_cosmic,
		                                                                h_dedx_slice_3_zoom_other_mixed,
		                                                                h_dedx_slice_3_zoom_unmatched);

		_cuts_instance.selection_cuts::dEdxCut(tpc_object_container_v, passed_tpco, tolerance_dedx_min, tolerance_dedx_max, _verbose, false);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, dedx_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_dedx, h_selected_ele_energy_dedx);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_9, h_reco_ele_e_9, h_mc_reco_ele_e_9);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_dedx_cut_nue_cc,
		                                                                 h_ele_momentum_dedx_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_dedx_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_dedx_cut_numu_cc,
		                                                                 h_ele_momentum_dedx_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_dedx_cut_nc,
		                                                                 h_ele_momentum_dedx_cut_nc_pi0,
		                                                                 h_ele_momentum_dedx_cut_cosmic,
		                                                                 h_ele_momentum_dedx_cut_other_mixed,
		                                                                 h_ele_momentum_dedx_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_dedx_cut_electron,
		                                                                             h_leading_momentum_dedx_cut_photon,
		                                                                             h_leading_momentum_dedx_cut_proton,
		                                                                             h_leading_momentum_dedx_cut_pion,
		                                                                             h_leading_momentum_dedx_cut_muon,
		                                                                             h_leading_momentum_dedx_cut_kaon,
		                                                                             h_leading_momentum_dedx_cut_neutron,
		                                                                             h_leading_momentum_dedx_cut_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_dedx_cut_nue_cc,
		                                                           h_multiplicity_shower_dedx_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_dedx_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_dedx_cut_numu_cc,
		                                                           h_multiplicity_shower_dedx_cut_nc,
		                                                           h_multiplicity_shower_dedx_cut_nc_pi0,
		                                                           h_multiplicity_shower_dedx_cut_cosmic,
		                                                           h_multiplicity_shower_dedx_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_dedx_cut_other_mixed,
		                                                           h_multiplicity_shower_dedx_cut_unmatched,
		                                                           h_multiplicity_track_dedx_cut_nue_cc,
		                                                           h_multiplicity_track_dedx_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_dedx_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_dedx_cut_numu_cc,
		                                                           h_multiplicity_track_dedx_cut_nc,
		                                                           h_multiplicity_track_dedx_cut_nc_pi0,
		                                                           h_multiplicity_track_dedx_cut_cosmic,
		                                                           h_multiplicity_track_dedx_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_dedx_cut_other_mixed,
		                                                           h_multiplicity_track_dedx_cut_unmatched);

		_functions_instance.selection_functions::PostCutsdEdx(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                      h_dedx_cuts_nue_cc_after, h_dedx_cuts_nue_cc_mixed_after,
		                                                      h_dedx_cuts_nue_cc_out_fv_after,
		                                                      h_dedx_cuts_numu_cc_after, h_dedx_cuts_nc_after,
		                                                      h_dedx_cuts_cosmic_after, h_dedx_cuts_nc_pi0_after,
		                                                      h_dedx_cuts_numu_cc_mixed_after, h_dedx_cuts_other_mixed_after,
		                                                      h_dedx_cuts_unmatched_after);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_dedx->Fill(mc_ele_energy); }
		//***************************************************************************
		// ******* Secondary Showers Distance Cut *****************
		//***************************************************************************
		_functions_instance.selection_functions::SecondaryShowersDist(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                              h_second_shwr_dist_nue_cc, h_second_shwr_dist_nue_cc_out_fv,
		                                                              h_second_shwr_dist_nue_cc_mixed, h_second_shwr_dist_numu_cc,
		                                                              h_second_shwr_dist_numu_cc_mixed, h_second_shwr_dist_nc,
		                                                              h_second_shwr_dist_nc_pi0, h_second_shwr_dist_cosmic,
		                                                              h_second_shwr_dist_other_mixed, h_second_shwr_dist_unmatched);

		_cuts_instance.selection_cuts::SecondaryShowersDistCut(tpc_object_container_v, passed_tpco, _verbose, dist_tolerance);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, secondary_shower_counter_v);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_10, h_reco_ele_e_10, h_mc_reco_ele_e_10);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_2shwr_cut_nue_cc,
		                                                                 h_ele_momentum_2shwr_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_2shwr_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_2shwr_cut_numu_cc,
		                                                                 h_ele_momentum_2shwr_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_2shwr_cut_nc,
		                                                                 h_ele_momentum_2shwr_cut_nc_pi0,
		                                                                 h_ele_momentum_2shwr_cut_cosmic,
		                                                                 h_ele_momentum_2shwr_cut_other_mixed,
		                                                                 h_ele_momentum_2shwr_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_2shwr_cut_electron,
		                                                                             h_leading_momentum_2shwr_cut_photon,
		                                                                             h_leading_momentum_2shwr_cut_proton,
		                                                                             h_leading_momentum_2shwr_cut_pion,
		                                                                             h_leading_momentum_2shwr_cut_muon,
		                                                                             h_leading_momentum_2shwr_cut_kaon,
		                                                                             h_leading_momentum_2shwr_cut_neutron,
		                                                                             h_leading_momentum_2shwr_cut_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_2shwr_cut_nue_cc,
		                                                           h_multiplicity_shower_2shwr_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_2shwr_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_2shwr_cut_numu_cc,
		                                                           h_multiplicity_shower_2shwr_cut_nc,
		                                                           h_multiplicity_shower_2shwr_cut_nc_pi0,
		                                                           h_multiplicity_shower_2shwr_cut_cosmic,
		                                                           h_multiplicity_shower_2shwr_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_2shwr_cut_other_mixed,
		                                                           h_multiplicity_shower_2shwr_cut_unmatched,
		                                                           h_multiplicity_track_2shwr_cut_nue_cc,
		                                                           h_multiplicity_track_2shwr_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_2shwr_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_2shwr_cut_numu_cc,
		                                                           h_multiplicity_track_2shwr_cut_nc,
		                                                           h_multiplicity_track_2shwr_cut_nc_pi0,
		                                                           h_multiplicity_track_2shwr_cut_cosmic,
		                                                           h_multiplicity_track_2shwr_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_2shwr_cut_other_mixed,
		                                                           h_multiplicity_track_2shwr_cut_unmatched);

		_functions_instance.selection_functions::SecondaryShowersDist(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                              h_second_shwr_dist_nue_cc_after, h_second_shwr_dist_nue_cc_out_fv_after,
		                                                              h_second_shwr_dist_nue_cc_mixed_after, h_second_shwr_dist_numu_cc_after,
		                                                              h_second_shwr_dist_numu_cc_mixed_after, h_second_shwr_dist_nc_after,
		                                                              h_second_shwr_dist_nc_pi0_after, h_second_shwr_dist_cosmic_after,
		                                                              h_second_shwr_dist_other_mixed_after, h_second_shwr_dist_unmatched_after);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_2shwr->Fill(mc_ele_energy); }
		//******************************************************************************
		// ********** Hit Length Ratio Cut *************
		//******************************************************************************
		_functions_instance.selection_functions::HitLengthRatio(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                        h_hit_length_ratio_nue_cc,
		                                                        h_hit_length_ratio_nue_cc_out_fv,
		                                                        h_hit_length_ratio_nue_cc_mixed,
		                                                        h_hit_length_ratio_numu_cc,
		                                                        h_hit_length_ratio_numu_cc_mixed,
		                                                        h_hit_length_ratio_nc,
		                                                        h_hit_length_ratio_nc_pi0,
		                                                        h_hit_length_ratio_cosmic,
		                                                        h_hit_length_ratio_other_mixed,
		                                                        h_hit_length_ratio_unmatched);

		_cuts_instance.selection_cuts::HitLengthRatioCut(tpc_object_container_v, passed_tpco, _verbose, pfp_hits_length_tolerance);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_lengthRatio_counter_v);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_11, h_reco_ele_e_11, h_mc_reco_ele_e_11);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_hit_length_cut_nue_cc,
		                                                                 h_ele_momentum_hit_length_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_hit_length_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_hit_length_cut_numu_cc,
		                                                                 h_ele_momentum_hit_length_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_hit_length_cut_nc,
		                                                                 h_ele_momentum_hit_length_cut_nc_pi0,
		                                                                 h_ele_momentum_hit_length_cut_cosmic,
		                                                                 h_ele_momentum_hit_length_cut_other_mixed,
		                                                                 h_ele_momentum_hit_length_cut_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_hit_length_cut_nue_cc,
		                                                           h_multiplicity_shower_hit_length_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_hit_length_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_hit_length_cut_numu_cc,
		                                                           h_multiplicity_shower_hit_length_cut_nc,
		                                                           h_multiplicity_shower_hit_length_cut_nc_pi0,
		                                                           h_multiplicity_shower_hit_length_cut_cosmic,
		                                                           h_multiplicity_shower_hit_length_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_hit_length_cut_other_mixed,
		                                                           h_multiplicity_shower_hit_length_cut_unmatched,
		                                                           h_multiplicity_track_hit_length_cut_nue_cc,
		                                                           h_multiplicity_track_hit_length_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_hit_length_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_hit_length_cut_numu_cc,
		                                                           h_multiplicity_track_hit_length_cut_nc,
		                                                           h_multiplicity_track_hit_length_cut_nc_pi0,
		                                                           h_multiplicity_track_hit_length_cut_cosmic,
		                                                           h_multiplicity_track_hit_length_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_hit_length_cut_other_mixed,
		                                                           h_multiplicity_track_hit_length_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_hit_length_cut_electron,
		                                                                             h_leading_momentum_hit_length_cut_photon,
		                                                                             h_leading_momentum_hit_length_cut_proton,
		                                                                             h_leading_momentum_hit_length_cut_pion,
		                                                                             h_leading_momentum_hit_length_cut_muon,
		                                                                             h_leading_momentum_hit_length_cut_kaon,
		                                                                             h_leading_momentum_hit_length_cut_neutron,
		                                                                             h_leading_momentum_hit_length_cut_mc_unmatched);

		_functions_instance.selection_functions::HitLengthRatio(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                        h_hit_length_ratio_nue_cc_after,
		                                                        h_hit_length_ratio_nue_cc_out_fv_after,
		                                                        h_hit_length_ratio_nue_cc_mixed_after,
		                                                        h_hit_length_ratio_numu_cc_after,
		                                                        h_hit_length_ratio_numu_cc_mixed_after,
		                                                        h_hit_length_ratio_nc_after,
		                                                        h_hit_length_ratio_nc_pi0_after,
		                                                        h_hit_length_ratio_cosmic_after,
		                                                        h_hit_length_ratio_other_mixed_after,
		                                                        h_hit_length_ratio_unmatched_after);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_hit_len->Fill(mc_ele_energy); }
		//******************************************************************************
		//*** cut for longest track / leading shower ratio *** //
		//******************************************************************************
		_functions_instance.selection_functions::LeadingShowerLength(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                             h_leading_shwr_length_nue_cc,
		                                                             h_leading_shwr_length_nue_cc_out_fv,
		                                                             h_leading_shwr_length_nue_cc_mixed,
		                                                             h_leading_shwr_length_numu_cc,
		                                                             h_leading_shwr_length_numu_cc_mixed,
		                                                             h_leading_shwr_length_nc,
		                                                             h_leading_shwr_length_nc_pi0,
		                                                             h_leading_shwr_length_cosmic,
		                                                             h_leading_shwr_length_other_mixed,
		                                                             h_leading_shwr_length_unmatched);

		_functions_instance.selection_functions::LeadingShowerTrackLengths(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                   h_leading_shwr_trk_length_nue_cc,
		                                                                   h_leading_shwr_trk_length_nue_cc_out_fv,
		                                                                   h_leading_shwr_trk_length_nue_cc_mixed,
		                                                                   h_leading_shwr_trk_length_numu_cc,
		                                                                   h_leading_shwr_trk_length_numu_cc_mixed,
		                                                                   h_leading_shwr_trk_length_nc,
		                                                                   h_leading_shwr_trk_length_nc_pi0,
		                                                                   h_leading_shwr_trk_length_cosmic,
		                                                                   h_leading_shwr_trk_length_other_mixed,
		                                                                   h_leading_shwr_trk_length_unmatched);

		_cuts_instance.selection_cuts::LongestTrackLeadingShowerCut(tpc_object_container_v, passed_tpco, _verbose, ratio_tolerance);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, trk_len_shwr_len_ratio_counter_v);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_12, h_reco_ele_e_12, h_mc_reco_ele_e_12);
		}
		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_length_ratio_cut_nue_cc,
		                                                                 h_ele_momentum_length_ratio_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_length_ratio_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_length_ratio_cut_numu_cc,
		                                                                 h_ele_momentum_length_ratio_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_length_ratio_cut_nc,
		                                                                 h_ele_momentum_length_ratio_cut_nc_pi0,
		                                                                 h_ele_momentum_length_ratio_cut_cosmic,
		                                                                 h_ele_momentum_length_ratio_cut_other_mixed,
		                                                                 h_ele_momentum_length_ratio_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_length_ratio_cut_electron,
		                                                                             h_leading_momentum_length_ratio_cut_photon,
		                                                                             h_leading_momentum_length_ratio_cut_proton,
		                                                                             h_leading_momentum_length_ratio_cut_pion,
		                                                                             h_leading_momentum_length_ratio_cut_muon,
		                                                                             h_leading_momentum_length_ratio_cut_kaon,
		                                                                             h_leading_momentum_length_ratio_cut_neutron,
		                                                                             h_leading_momentum_length_ratio_cut_mc_unmatched);

		_functions_instance.selection_functions::LeadingShowerLength(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                             h_leading_shwr_length_nue_cc_after,
		                                                             h_leading_shwr_length_nue_cc_out_fv_after,
		                                                             h_leading_shwr_length_nue_cc_mixed_after,
		                                                             h_leading_shwr_length_numu_cc_after,
		                                                             h_leading_shwr_length_numu_cc_mixed_after,
		                                                             h_leading_shwr_length_nc_after,
		                                                             h_leading_shwr_length_nc_pi0_after,
		                                                             h_leading_shwr_length_cosmic_after,
		                                                             h_leading_shwr_length_other_mixed_after,
		                                                             h_leading_shwr_length_unmatched_after);

		_functions_instance.selection_functions::LeadingShowerTrackLengths(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                   h_leading_shwr_trk_length_nue_cc_after,
		                                                                   h_leading_shwr_trk_length_nue_cc_out_fv_after,
		                                                                   h_leading_shwr_trk_length_nue_cc_mixed_after,
		                                                                   h_leading_shwr_trk_length_numu_cc_after,
		                                                                   h_leading_shwr_trk_length_numu_cc_mixed_after,
		                                                                   h_leading_shwr_trk_length_nc_after,
		                                                                   h_leading_shwr_trk_length_nc_pi0_after,
		                                                                   h_leading_shwr_trk_length_cosmic_after,
		                                                                   h_leading_shwr_trk_length_other_mixed_after,
		                                                                   h_leading_shwr_trk_length_unmatched_after);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_length_ratio_cut_nue_cc,
		                                                           h_multiplicity_shower_length_ratio_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_length_ratio_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_length_ratio_cut_numu_cc,
		                                                           h_multiplicity_shower_length_ratio_cut_nc,
		                                                           h_multiplicity_shower_length_ratio_cut_nc_pi0,
		                                                           h_multiplicity_shower_length_ratio_cut_cosmic,
		                                                           h_multiplicity_shower_length_ratio_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_length_ratio_cut_other_mixed,
		                                                           h_multiplicity_shower_length_ratio_cut_unmatched,
		                                                           h_multiplicity_track_length_ratio_cut_nue_cc,
		                                                           h_multiplicity_track_length_ratio_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_length_ratio_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_length_ratio_cut_numu_cc,
		                                                           h_multiplicity_track_length_ratio_cut_nc,
		                                                           h_multiplicity_track_length_ratio_cut_nc_pi0,
		                                                           h_multiplicity_track_length_ratio_cut_cosmic,
		                                                           h_multiplicity_track_length_ratio_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_length_ratio_cut_other_mixed,
		                                                           h_multiplicity_track_length_ratio_cut_unmatched);


		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_trk_shwr->Fill(mc_ele_energy); }
		//******************************************************************************
		//*** track containment -- all tracks need to be contained *** //
		//******************************************************************************
		_functions_instance.selection_functions::IsContainedPlot(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v, fv_boundary_v,
		                                                         h_track_containment_nue_cc,
		                                                         h_track_containment_nue_cc_out_fv,
		                                                         h_track_containment_nue_cc_mixed,
		                                                         h_track_containment_numu_cc,
		                                                         h_track_containment_numu_cc_mixed,
		                                                         h_track_containment_nc,
		                                                         h_track_containment_nc_pi0,
		                                                         h_track_containment_cosmic,
		                                                         h_track_containment_other_mixed,
		                                                         h_track_containment_unmatched);

		_cuts_instance.selection_cuts::ContainedTracksCut(tpc_object_container_v, passed_tpco, _verbose, fv_boundary_v, true);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, track_containment_counter_v);

		_functions_instance.selection_functions::PostCutsLeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_momentum_containment_cut_nue_cc,
		                                                                 h_ele_momentum_containment_cut_nue_cc_out_fv,
		                                                                 h_ele_momentum_containment_cut_nue_cc_mixed,
		                                                                 h_ele_momentum_containment_cut_numu_cc,
		                                                                 h_ele_momentum_containment_cut_numu_cc_mixed,
		                                                                 h_ele_momentum_containment_cut_nc,
		                                                                 h_ele_momentum_containment_cut_nc_pi0,
		                                                                 h_ele_momentum_containment_cut_cosmic,
		                                                                 h_ele_momentum_containment_cut_other_mixed,
		                                                                 h_ele_momentum_containment_cut_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                             h_leading_momentum_containment_cut_electron,
		                                                                             h_leading_momentum_containment_cut_photon,
		                                                                             h_leading_momentum_containment_cut_proton,
		                                                                             h_leading_momentum_containment_cut_pion,
		                                                                             h_leading_momentum_containment_cut_muon,
		                                                                             h_leading_momentum_containment_cut_kaon,
		                                                                             h_leading_momentum_containment_cut_neutron,
		                                                                             h_leading_momentum_containment_cut_mc_unmatched);

		_functions_instance.selection_functions::EventMultiplicity(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_multiplicity_shower_containment_cut_nue_cc,
		                                                           h_multiplicity_shower_containment_cut_nue_cc_mixed,
		                                                           h_multiplicity_shower_containment_cut_nue_cc_out_fv,
		                                                           h_multiplicity_shower_containment_cut_numu_cc,
		                                                           h_multiplicity_shower_containment_cut_nc,
		                                                           h_multiplicity_shower_containment_cut_nc_pi0,
		                                                           h_multiplicity_shower_containment_cut_cosmic,
		                                                           h_multiplicity_shower_containment_cut_numu_cc_mixed,
		                                                           h_multiplicity_shower_containment_cut_other_mixed,
		                                                           h_multiplicity_shower_containment_cut_unmatched,
		                                                           h_multiplicity_track_containment_cut_nue_cc,
		                                                           h_multiplicity_track_containment_cut_nue_cc_mixed,
		                                                           h_multiplicity_track_containment_cut_nue_cc_out_fv,
		                                                           h_multiplicity_track_containment_cut_numu_cc,
		                                                           h_multiplicity_track_containment_cut_nc,
		                                                           h_multiplicity_track_containment_cut_nc_pi0,
		                                                           h_multiplicity_track_containment_cut_cosmic,
		                                                           h_multiplicity_track_containment_cut_numu_cc_mixed,
		                                                           h_multiplicity_track_containment_cut_other_mixed,
		                                                           h_multiplicity_track_containment_cut_unmatched);

		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true && tabulated_origins->at(0) >= 1) {h_ele_eng_eff_contain->Fill(mc_ele_energy); }
		//*************************************
		// ******** End Selection Cuts! ******
		//*************************************
		_functions_instance.selection_functions::LeadingPhi(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                    h_ele_pfp_phi_after_nue_cc,
		                                                    h_ele_pfp_phi_after_nue_cc_out_fv,
		                                                    h_ele_pfp_phi_after_nue_cc_mixed,
		                                                    h_ele_pfp_phi_after_numu_cc,
		                                                    h_ele_pfp_phi_after_numu_cc_mixed,
		                                                    h_ele_pfp_phi_after_nc,
		                                                    h_ele_pfp_phi_after_nc_pi0,
		                                                    h_ele_pfp_phi_after_cosmic,
		                                                    h_ele_pfp_phi_after_other_mixed,
		                                                    h_ele_pfp_phi_after_unmatched);

		_functions_instance.selection_functions::LeadingTheta(tpc_object_container_v, passed_tpco, 0, 0, _verbose, tpco_classifier_v,
		                                                      h_ele_pfp_theta_after_nue_cc,
		                                                      h_ele_pfp_theta_after_nue_cc_out_fv,
		                                                      h_ele_pfp_theta_after_nue_cc_mixed,
		                                                      h_ele_pfp_theta_after_numu_cc,
		                                                      h_ele_pfp_theta_after_numu_cc_mixed,
		                                                      h_ele_pfp_theta_after_nc,
		                                                      h_ele_pfp_theta_after_nc_pi0,
		                                                      h_ele_pfp_theta_after_cosmic,
		                                                      h_ele_pfp_theta_after_other_mixed,
		                                                      h_ele_pfp_theta_after_unmatched);

		_functions_instance.selection_functions::XYZPosition(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                     mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                     h_any_pfp_xyz_last_nue_cc,
		                                                     h_any_pfp_xyz_last_nue_cc_out_fv,
		                                                     h_any_pfp_xyz_last_nue_cc_mixed,
		                                                     h_any_pfp_xyz_last_numu_cc,
		                                                     h_any_pfp_xyz_last_numu_cc_mixed,
		                                                     h_any_pfp_xyz_last_nc,
		                                                     h_any_pfp_xyz_last_nc_pi0,
		                                                     h_any_pfp_xyz_last_cosmic,
		                                                     h_any_pfp_xyz_last_other_mixed,
		                                                     h_any_pfp_xyz_last_unmatched);

		_functions_instance.selection_functions::dEdxTheta(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                   h_dedx_theta_nue_cc,        h_dedx_theta_nue_cc_mixed,
		                                                   h_dedx_theta_nue_cc_out_fv,
		                                                   h_dedx_theta_numu_cc,       h_dedx_theta_nc,
		                                                   h_dedx_theta_cosmic,        h_dedx_theta_nc_pi0,
		                                                   h_dedx_theta_numu_cc_mixed, h_dedx_theta_other_mixed,
		                                                   h_dedx_theta_unmatched);

		_functions_instance.selection_functions::PostCutsdEdxTrueParticle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                  h_dedx_cuts_last_electron, h_dedx_cuts_last_photon,
		                                                                  h_dedx_cuts_last_proton,   h_dedx_cuts_last_pion,
		                                                                  h_dedx_cuts_last_muon,     h_dedx_cuts_last_kaon,
		                                                                  h_dedx_cuts_last_neutron,  h_dedx_cuts_last_mc_unmatched);

		_functions_instance.selection_functions::FailureReason(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                       h_failure_reason_nue_cc, h_failure_reason_nue_cc_out_fv,
		                                                       h_failure_reason_nue_cc_mixed, h_failure_reason_numu_cc,
		                                                       h_failure_reason_numu_cc_mixed, h_failure_reason_nc,
		                                                       h_failure_reason_nc_pi0, h_failure_reason_cosmic,
		                                                       h_failure_reason_other_mixed, h_failure_reason_unmatched);

		//these are for the tefficiency plots, post all cuts
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			// int pass = 0;
			// for(auto const passed_tpco_pair : * passed_tpco) {pass += passed_tpco_pair.first; }
			// if(pass > 0)
			_functions_instance.selection_functions::FillTrueRecoEnergy(tpc_object_container_v, passed_tpco, tpco_classifier_v, mc_ele_energy,
			                                                            h_mc_ele_e_13, h_reco_ele_e_13, h_mc_reco_ele_e_13);
			if(tabulated_origins->at(0) >= 1)
			{
				selected_energy_vector.push_back(mc_nu_energy);
				h_nue_eng_eff_num->Fill(mc_nu_energy);
				//h_ele_eng_eff_num->Fill(mc_ele_energy);
				h_nue_vtx_x_eff_num->Fill(mc_nu_vtx_x);
				h_nue_vtx_y_eff_num->Fill(mc_nu_vtx_y);
				h_nue_vtx_z_eff_num->Fill(mc_nu_vtx_z);
				h_nue_dir_x_eff_num->Fill(mc_nu_dir_x);
				h_nue_dir_y_eff_num->Fill(mc_nu_dir_y);
				h_nue_dir_z_eff_num->Fill(mc_nu_dir_z);
				h_ele_dir_x_eff_num->Fill(mc_ele_dir_x);
				h_ele_dir_y_eff_num->Fill(mc_ele_dir_y);
				h_ele_dir_z_eff_num->Fill(mc_ele_dir_z);
				h_ele_theta_eff_num->Fill(mc_ele_theta);
				h_nue_num_part_eff_num->Fill(mc_nu_num_particles);
				h_nue_num_chrg_part_eff_num->Fill(mc_nu_num_charged_particles);
				h_nue_cos_theta_eff_num->Fill(mc_cos_theta);
				h_nue_phi_eff_num->Fill(mc_phi);
				h_ele_cos_theta_eff_num->Fill(mc_ele_cos_theta);
				h_ele_phi_eff_num->Fill(mc_ele_phi * (180/3.1415));
				_functions_instance.selection_functions::EnergyHits(tpc_object_container_v, passed_tpco, _verbose,
				                                                    tpco_classifier_v, mc_nu_energy, mc_ele_energy,
				                                                    h_ele_eng_total_hits, h_ele_eng_colleciton_hits, h_nu_eng_total_hits, h_nu_eng_collection_hits);
				_functions_instance.selection_functions::TrueRecoEle(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v, mc_ele_momentum, mc_ele_cos_theta,
				                                                     h_true_reco_ele_momentum, h_true_reco_ele_costheta, h_true_num_e);
				_functions_instance.selection_functions::TrueEleResolution(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
				                                                           mc_ele_momentum, mc_ele_phi  * (180/3.1415), mc_ele_theta,
				                                                           mc_ele_dir_x, mc_ele_dir_y, mc_ele_dir_z,
				                                                           h_ele_resolution_momentum, h_ele_resolution_phi, h_ele_resolution_theta,
				                                                           h_ele_resolution_dot_prod, h_ele_resolution_momentum_dot_prod,
				                                                           h_ele_resolution_momentum_dot_prod_zoom_y);
			}
		}

		_functions_instance.selection_functions::TopologyPlots2(tpc_object_container_v, passed_tpco, tpco_classifier_v,
		                                                        h_pfp_track_shower_nue_cc_qe_last, h_pfp_track_shower_nue_cc_out_fv_last,
		                                                        h_pfp_track_shower_nue_cc_res_last, h_pfp_track_shower_nue_cc_dis_last,
		                                                        h_pfp_track_shower_nue_cc_coh_last, h_pfp_track_shower_nue_cc_mec_last,
		                                                        h_pfp_track_shower_nc_last, h_pfp_track_shower_numu_cc_qe_last,
		                                                        h_pfp_track_shower_numu_cc_res_last, h_pfp_track_shower_numu_cc_dis_last,
		                                                        h_pfp_track_shower_numu_cc_coh_last, h_pfp_track_shower_numu_cc_mec_last,
		                                                        h_pfp_track_shower_nc_pi0_last, h_pfp_track_shower_nue_cc_mixed_last,
		                                                        h_pfp_track_shower_numu_cc_mixed_last, h_pfp_track_shower_cosmic_last,
		                                                        h_pfp_track_shower_other_mixed_last, h_pfp_track_shower_unmatched_last,
		                                                        h_leading_shower_mc_pdg_nue_cc_qe_last, h_leading_shower_mc_pdg_nue_cc_out_fv_last,
		                                                        h_leading_shower_mc_pdg_nue_cc_res_last, h_leading_shower_mc_pdg_nue_cc_dis_last,
		                                                        h_leading_shower_mc_pdg_nue_cc_coh_last, h_leading_shower_mc_pdg_nue_cc_mec_last,
		                                                        h_leading_shower_mc_pdg_nc_last, h_leading_shower_mc_pdg_numu_cc_qe_last,
		                                                        h_leading_shower_mc_pdg_numu_cc_res_last, h_leading_shower_mc_pdg_numu_cc_dis_last,
		                                                        h_leading_shower_mc_pdg_numu_cc_coh_last, h_leading_shower_mc_pdg_numu_cc_mec_last,
		                                                        h_leading_shower_mc_pdg_nc_pi0_last, h_leading_shower_mc_pdg_nue_cc_mixed_last,
		                                                        h_leading_shower_mc_pdg_numu_cc_mixed_last, h_leading_shower_mc_pdg_cosmic_last,
		                                                        h_leading_shower_mc_pdg_other_mixed_last, h_leading_shower_mc_pdg_unmatched_last,
		                                                        h_pfp_track_nue_cc_qe_last, h_pfp_track_nue_cc_out_fv_last,
		                                                        h_pfp_track_nue_cc_res_last, h_pfp_track_nue_cc_dis_last,
		                                                        h_pfp_track_nue_cc_coh_last, h_pfp_track_nue_cc_mec_last,
		                                                        h_pfp_track_nc_last, h_pfp_track_numu_cc_qe_last,
		                                                        h_pfp_track_numu_cc_res_last, h_pfp_track_numu_cc_dis_last,
		                                                        h_pfp_track_numu_cc_coh_last, h_pfp_track_numu_cc_mec_last,
		                                                        h_pfp_track_nc_pi0_last, h_pfp_track_nue_cc_mixed_last,
		                                                        h_pfp_track_numu_cc_mixed_last, h_pfp_track_cosmic_last,
		                                                        h_pfp_track_other_mixed_last, h_pfp_track_unmatched_last,
		                                                        h_pfp_shower_nue_cc_qe_last, h_pfp_shower_nue_cc_out_fv_last,
		                                                        h_pfp_shower_nue_cc_res_last, h_pfp_shower_nue_cc_dis_last,
		                                                        h_pfp_shower_nue_cc_coh_last, h_pfp_shower_nue_cc_mec_last,
		                                                        h_pfp_shower_nc_last, h_pfp_shower_numu_cc_qe_last,
		                                                        h_pfp_shower_numu_cc_res_last, h_pfp_shower_numu_cc_dis_last,
		                                                        h_pfp_shower_numu_cc_coh_last, h_pfp_shower_numu_cc_mec_last,
		                                                        h_pfp_shower_nc_pi0_last, h_pfp_shower_nue_cc_mixed_last,
		                                                        h_pfp_shower_numu_cc_mixed_last, h_pfp_shower_cosmic_last,
		                                                        h_pfp_shower_other_mixed_last, h_pfp_shower_unmatched_last);

		_functions_instance.selection_functions::LeadingKinematicsShowerTopology(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                         h_ele_pfp_momentum_1shwr_nue_cc,
		                                                                         h_ele_pfp_momentum_1shwr_nue_cc_out_fv,
		                                                                         h_ele_pfp_momentum_1shwr_nue_cc_mixed,
		                                                                         h_ele_pfp_momentum_1shwr_numu_cc,
		                                                                         h_ele_pfp_momentum_1shwr_numu_cc_mixed,
		                                                                         h_ele_pfp_momentum_1shwr_nc,
		                                                                         h_ele_pfp_momentum_1shwr_nc_pi0,
		                                                                         h_ele_pfp_momentum_1shwr_cosmic,
		                                                                         h_ele_pfp_momentum_1shwr_other_mixed,
		                                                                         h_ele_pfp_momentum_1shwr_unmatched,
		                                                                         h_ele_pfp_momentum_2shwr_nue_cc,
		                                                                         h_ele_pfp_momentum_2shwr_nue_cc_out_fv,
		                                                                         h_ele_pfp_momentum_2shwr_nue_cc_mixed,
		                                                                         h_ele_pfp_momentum_2shwr_numu_cc,
		                                                                         h_ele_pfp_momentum_2shwr_numu_cc_mixed,
		                                                                         h_ele_pfp_momentum_2shwr_nc,
		                                                                         h_ele_pfp_momentum_2shwr_nc_pi0,
		                                                                         h_ele_pfp_momentum_2shwr_cosmic,
		                                                                         h_ele_pfp_momentum_2shwr_other_mixed,
		                                                                         h_ele_pfp_momentum_2shwr_unmatched,
		                                                                         h_ele_pfp_theta_1shwr_nue_cc,
		                                                                         h_ele_pfp_theta_1shwr_nue_cc_out_fv,
		                                                                         h_ele_pfp_theta_1shwr_nue_cc_mixed,
		                                                                         h_ele_pfp_theta_1shwr_numu_cc,
		                                                                         h_ele_pfp_theta_1shwr_numu_cc_mixed,
		                                                                         h_ele_pfp_theta_1shwr_nc,
		                                                                         h_ele_pfp_theta_1shwr_nc_pi0,
		                                                                         h_ele_pfp_theta_1shwr_cosmic,
		                                                                         h_ele_pfp_theta_1shwr_other_mixed,
		                                                                         h_ele_pfp_theta_1shwr_unmatched,
		                                                                         h_ele_pfp_theta_2shwr_nue_cc,
		                                                                         h_ele_pfp_theta_2shwr_nue_cc_out_fv,
		                                                                         h_ele_pfp_theta_2shwr_nue_cc_mixed,
		                                                                         h_ele_pfp_theta_2shwr_numu_cc,
		                                                                         h_ele_pfp_theta_2shwr_numu_cc_mixed,
		                                                                         h_ele_pfp_theta_2shwr_nc,
		                                                                         h_ele_pfp_theta_2shwr_nc_pi0,
		                                                                         h_ele_pfp_theta_2shwr_cosmic,
		                                                                         h_ele_pfp_theta_2shwr_other_mixed,
		                                                                         h_ele_pfp_theta_2shwr_unmatched,
		                                                                         h_ele_pfp_phi_1shwr_nue_cc,
		                                                                         h_ele_pfp_phi_1shwr_nue_cc_out_fv,
		                                                                         h_ele_pfp_phi_1shwr_nue_cc_mixed,
		                                                                         h_ele_pfp_phi_1shwr_numu_cc,
		                                                                         h_ele_pfp_phi_1shwr_numu_cc_mixed,
		                                                                         h_ele_pfp_phi_1shwr_nc,
		                                                                         h_ele_pfp_phi_1shwr_nc_pi0,
		                                                                         h_ele_pfp_phi_1shwr_cosmic,
		                                                                         h_ele_pfp_phi_1shwr_other_mixed,
		                                                                         h_ele_pfp_phi_1shwr_unmatched,
		                                                                         h_ele_pfp_phi_2shwr_nue_cc,
		                                                                         h_ele_pfp_phi_2shwr_nue_cc_out_fv,
		                                                                         h_ele_pfp_phi_2shwr_nue_cc_mixed,
		                                                                         h_ele_pfp_phi_2shwr_numu_cc,
		                                                                         h_ele_pfp_phi_2shwr_numu_cc_mixed,
		                                                                         h_ele_pfp_phi_2shwr_nc,
		                                                                         h_ele_pfp_phi_2shwr_nc_pi0,
		                                                                         h_ele_pfp_phi_2shwr_cosmic,
		                                                                         h_ele_pfp_phi_2shwr_other_mixed,
		                                                                         h_ele_pfp_phi_2shwr_unmatched);

		_functions_instance.selection_functions::TopologyEfficiency(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                            no_track, has_track, _1_shwr, _2_shwr, _3_shwr, _4_shwr);

		_functions_instance.selection_functions::ChargeShare(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                     h_charge_share_nue_cc_mixed);

		_functions_instance.selection_functions::FillPostCutVector(tpc_object_container_v, passed_tpco, tpco_classifier_v, post_cuts_v);

		_functions_instance.selection_functions::TrackLength(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                     h_trk_length_nue_cc,
		                                                     h_trk_length_nue_cc_out_fv,
		                                                     h_trk_length_nue_cc_mixed,
		                                                     h_trk_length_numu_cc,
		                                                     h_trk_length_numu_cc_mixed,
		                                                     h_trk_length_nc,
		                                                     h_trk_length_nc_pi0,
		                                                     h_trk_length_cosmic,
		                                                     h_trk_length_other_mixed,
		                                                     h_trk_length_unmatched);
		_functions_instance.selection_functions::LongestTrackLength(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                            h_longest_trk_length_nue_cc,
		                                                            h_longest_trk_length_nue_cc_out_fv,
		                                                            h_longest_trk_length_nue_cc_mixed,
		                                                            h_longest_trk_length_numu_cc,
		                                                            h_longest_trk_length_numu_cc_mixed,
		                                                            h_longest_trk_length_nc,
		                                                            h_longest_trk_length_nc_pi0,
		                                                            h_longest_trk_length_cosmic,
		                                                            h_longest_trk_length_other_mixed,
		                                                            h_longest_trk_length_unmatched);

		_functions_instance.selection_functions::ShowerLength(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                      h_shwr_length_nue_cc,
		                                                      h_shwr_length_nue_cc_out_fv,
		                                                      h_shwr_length_nue_cc_mixed,
		                                                      h_shwr_length_numu_cc,
		                                                      h_shwr_length_numu_cc_mixed,
		                                                      h_shwr_length_nc,
		                                                      h_shwr_length_nc_pi0,
		                                                      h_shwr_length_cosmic,
		                                                      h_shwr_length_other_mixed,
		                                                      h_shwr_length_unmatched);
		_functions_instance.selection_functions::LongestShowerLength(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                             h_longest_shwr_length_nue_cc,
		                                                             h_longest_shwr_length_nue_cc_out_fv,
		                                                             h_longest_shwr_length_nue_cc_mixed,
		                                                             h_longest_shwr_length_numu_cc,
		                                                             h_longest_shwr_length_numu_cc_mixed,
		                                                             h_longest_shwr_length_nc,
		                                                             h_longest_shwr_length_nc_pi0,
		                                                             h_longest_shwr_length_cosmic,
		                                                             h_longest_shwr_length_other_mixed,
		                                                             h_longest_shwr_length_unmatched);

		_functions_instance.selection_functions::LongestShowerTrackLengths(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                   h_longest_shwr_trk_length_nue_cc,
		                                                                   h_longest_shwr_trk_length_nue_cc_out_fv,
		                                                                   h_longest_shwr_trk_length_nue_cc_mixed,
		                                                                   h_longest_shwr_trk_length_numu_cc,
		                                                                   h_longest_shwr_trk_length_numu_cc_mixed,
		                                                                   h_longest_shwr_trk_length_nc,
		                                                                   h_longest_shwr_trk_length_nc_pi0,
		                                                                   h_longest_shwr_trk_length_cosmic,
		                                                                   h_longest_shwr_trk_length_other_mixed,
		                                                                   h_longest_shwr_trk_length_unmatched);
		_functions_instance.selection_functions::LeadingCosTheta(tpc_object_container_v, passed_tpco, 0, 0, _verbose, tpco_classifier_v,
		                                                         h_ele_cos_theta_last_nue_cc,
		                                                         h_ele_cos_theta_last_nue_cc_out_fv,
		                                                         h_ele_cos_theta_last_nue_cc_mixed,
		                                                         h_ele_cos_theta_last_numu_cc,
		                                                         h_ele_cos_theta_last_numu_cc_mixed,
		                                                         h_ele_cos_theta_last_nc,
		                                                         h_ele_cos_theta_last_nc_pi0,
		                                                         h_ele_cos_theta_last_cosmic,
		                                                         h_ele_cos_theta_last_other_mixed,
		                                                         h_ele_cos_theta_last_unmatched);
		_functions_instance.selection_functions::LeadingCosTheta(tpc_object_container_v, passed_tpco, theta_translation, phi_translation, _verbose, tpco_classifier_v,
		                                                         h_ele_cos_theta_last_trans_nue_cc,
		                                                         h_ele_cos_theta_last_trans_nue_cc_out_fv,
		                                                         h_ele_cos_theta_last_trans_nue_cc_mixed,
		                                                         h_ele_cos_theta_last_trans_numu_cc,
		                                                         h_ele_cos_theta_last_trans_numu_cc_mixed,
		                                                         h_ele_cos_theta_last_trans_nc,
		                                                         h_ele_cos_theta_last_trans_nc_pi0,
		                                                         h_ele_cos_theta_last_trans_cosmic,
		                                                         h_ele_cos_theta_last_trans_other_mixed,
		                                                         h_ele_cos_theta_last_trans_unmatched);

		_functions_instance.selection_functions::LeadingMomentum(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                         h_ele_pfp_momentum_nue_cc,
		                                                         h_ele_pfp_momentum_nue_cc_out_fv,
		                                                         h_ele_pfp_momentum_nue_cc_mixed,
		                                                         h_ele_pfp_momentum_numu_cc,
		                                                         h_ele_pfp_momentum_numu_cc_mixed,
		                                                         h_ele_pfp_momentum_nc,
		                                                         h_ele_pfp_momentum_nc_pi0,
		                                                         h_ele_pfp_momentum_cosmic,
		                                                         h_ele_pfp_momentum_other_mixed,
		                                                         h_ele_pfp_momentum_unmatched);

		_functions_instance.selection_functions::LeadingMomentumTrackTopology(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                      h_ele_pfp_momentum_no_track_nue_cc,
		                                                                      h_ele_pfp_momentum_no_track_nue_cc_out_fv,
		                                                                      h_ele_pfp_momentum_no_track_nue_cc_mixed,
		                                                                      h_ele_pfp_momentum_no_track_numu_cc,
		                                                                      h_ele_pfp_momentum_no_track_numu_cc_mixed,
		                                                                      h_ele_pfp_momentum_no_track_nc,
		                                                                      h_ele_pfp_momentum_no_track_nc_pi0,
		                                                                      h_ele_pfp_momentum_no_track_cosmic,
		                                                                      h_ele_pfp_momentum_no_track_other_mixed,
		                                                                      h_ele_pfp_momentum_no_track_unmatched,
		                                                                      h_ele_pfp_momentum_has_track_nue_cc,
		                                                                      h_ele_pfp_momentum_has_track_nue_cc_out_fv,
		                                                                      h_ele_pfp_momentum_has_track_nue_cc_mixed,
		                                                                      h_ele_pfp_momentum_has_track_numu_cc,
		                                                                      h_ele_pfp_momentum_has_track_numu_cc_mixed,
		                                                                      h_ele_pfp_momentum_has_track_nc,
		                                                                      h_ele_pfp_momentum_has_track_nc_pi0,
		                                                                      h_ele_pfp_momentum_has_track_cosmic,
		                                                                      h_ele_pfp_momentum_has_track_other_mixed,
		                                                                      h_ele_pfp_momentum_has_track_unmatched);

		_functions_instance.selection_functions::LeadingPhi(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                    h_ele_pfp_phi_last_nue_cc,
		                                                    h_ele_pfp_phi_last_nue_cc_out_fv,
		                                                    h_ele_pfp_phi_last_nue_cc_mixed,
		                                                    h_ele_pfp_phi_last_numu_cc,
		                                                    h_ele_pfp_phi_last_numu_cc_mixed,
		                                                    h_ele_pfp_phi_last_nc,
		                                                    h_ele_pfp_phi_last_nc_pi0,
		                                                    h_ele_pfp_phi_last_cosmic,
		                                                    h_ele_pfp_phi_last_other_mixed,
		                                                    h_ele_pfp_phi_last_unmatched);

		_functions_instance.selection_functions::LeadingPhiTrackTopology(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                 h_ele_pfp_phi_no_track_nue_cc,
		                                                                 h_ele_pfp_phi_no_track_nue_cc_out_fv,
		                                                                 h_ele_pfp_phi_no_track_nue_cc_mixed,
		                                                                 h_ele_pfp_phi_no_track_numu_cc,
		                                                                 h_ele_pfp_phi_no_track_numu_cc_mixed,
		                                                                 h_ele_pfp_phi_no_track_nc,
		                                                                 h_ele_pfp_phi_no_track_nc_pi0,
		                                                                 h_ele_pfp_phi_no_track_cosmic,
		                                                                 h_ele_pfp_phi_no_track_other_mixed,
		                                                                 h_ele_pfp_phi_no_track_unmatched,
		                                                                 h_ele_pfp_phi_has_track_nue_cc,
		                                                                 h_ele_pfp_phi_has_track_nue_cc_out_fv,
		                                                                 h_ele_pfp_phi_has_track_nue_cc_mixed,
		                                                                 h_ele_pfp_phi_has_track_numu_cc,
		                                                                 h_ele_pfp_phi_has_track_numu_cc_mixed,
		                                                                 h_ele_pfp_phi_has_track_nc,
		                                                                 h_ele_pfp_phi_has_track_nc_pi0,
		                                                                 h_ele_pfp_phi_has_track_cosmic,
		                                                                 h_ele_pfp_phi_has_track_other_mixed,
		                                                                 h_ele_pfp_phi_has_track_unmatched);

		_functions_instance.selection_functions::LeadingTheta(tpc_object_container_v, passed_tpco, theta_translation, phi_translation,
		                                                      _verbose, tpco_classifier_v,
		                                                      h_ele_pfp_theta_last_nue_cc,
		                                                      h_ele_pfp_theta_last_nue_cc_out_fv,
		                                                      h_ele_pfp_theta_last_nue_cc_mixed,
		                                                      h_ele_pfp_theta_last_numu_cc,
		                                                      h_ele_pfp_theta_last_numu_cc_mixed,
		                                                      h_ele_pfp_theta_last_nc,
		                                                      h_ele_pfp_theta_last_nc_pi0,
		                                                      h_ele_pfp_theta_last_cosmic,
		                                                      h_ele_pfp_theta_last_other_mixed,
		                                                      h_ele_pfp_theta_last_unmatched);

		_functions_instance.selection_functions::LeadingThetaTrackTopology(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                   h_ele_pfp_theta_no_track_nue_cc,
		                                                                   h_ele_pfp_theta_no_track_nue_cc_out_fv,
		                                                                   h_ele_pfp_theta_no_track_nue_cc_mixed,
		                                                                   h_ele_pfp_theta_no_track_numu_cc,
		                                                                   h_ele_pfp_theta_no_track_numu_cc_mixed,
		                                                                   h_ele_pfp_theta_no_track_nc,
		                                                                   h_ele_pfp_theta_no_track_nc_pi0,
		                                                                   h_ele_pfp_theta_no_track_cosmic,
		                                                                   h_ele_pfp_theta_no_track_other_mixed,
		                                                                   h_ele_pfp_theta_no_track_unmatched,
		                                                                   h_ele_pfp_theta_has_track_nue_cc,
		                                                                   h_ele_pfp_theta_has_track_nue_cc_out_fv,
		                                                                   h_ele_pfp_theta_has_track_nue_cc_mixed,
		                                                                   h_ele_pfp_theta_has_track_numu_cc,
		                                                                   h_ele_pfp_theta_has_track_numu_cc_mixed,
		                                                                   h_ele_pfp_theta_has_track_nc,
		                                                                   h_ele_pfp_theta_has_track_nc_pi0,
		                                                                   h_ele_pfp_theta_has_track_cosmic,
		                                                                   h_ele_pfp_theta_has_track_other_mixed,
		                                                                   h_ele_pfp_theta_has_track_unmatched);

		_functions_instance.selection_functions::Leading1Shwr2Shwr(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                           h_leading_shwr_length_1shwr_nue_cc,
		                                                           h_leading_shwr_length_1shwr_nue_cc_out_fv,
		                                                           h_leading_shwr_length_1shwr_nue_cc_mixed,
		                                                           h_leading_shwr_length_1shwr_numu_cc,
		                                                           h_leading_shwr_length_1shwr_numu_cc_mixed,
		                                                           h_leading_shwr_length_1shwr_nc,
		                                                           h_leading_shwr_length_1shwr_nc_pi0,
		                                                           h_leading_shwr_length_1shwr_cosmic,
		                                                           h_leading_shwr_length_1shwr_other_mixed,
		                                                           h_leading_shwr_length_1shwr_unmatched,
		                                                           h_leading_shwr_length_2shwr_nue_cc,
		                                                           h_leading_shwr_length_2shwr_nue_cc_out_fv,
		                                                           h_leading_shwr_length_2shwr_nue_cc_mixed,
		                                                           h_leading_shwr_length_2shwr_numu_cc,
		                                                           h_leading_shwr_length_2shwr_numu_cc_mixed,
		                                                           h_leading_shwr_length_2shwr_nc,
		                                                           h_leading_shwr_length_2shwr_nc_pi0,
		                                                           h_leading_shwr_length_2shwr_cosmic,
		                                                           h_leading_shwr_length_2shwr_other_mixed,
		                                                           h_leading_shwr_length_2shwr_unmatched,
		                                                           h_leading_shwr_hits_1shwr_nue_cc,
		                                                           h_leading_shwr_hits_1shwr_nue_cc_out_fv,
		                                                           h_leading_shwr_hits_1shwr_nue_cc_mixed,
		                                                           h_leading_shwr_hits_1shwr_numu_cc,
		                                                           h_leading_shwr_hits_1shwr_numu_cc_mixed,
		                                                           h_leading_shwr_hits_1shwr_nc,
		                                                           h_leading_shwr_hits_1shwr_nc_pi0,
		                                                           h_leading_shwr_hits_1shwr_cosmic,
		                                                           h_leading_shwr_hits_1shwr_other_mixed,
		                                                           h_leading_shwr_hits_1shwr_unmatched,
		                                                           h_leading_shwr_hits_2shwr_nue_cc,
		                                                           h_leading_shwr_hits_2shwr_nue_cc_out_fv,
		                                                           h_leading_shwr_hits_2shwr_nue_cc_mixed,
		                                                           h_leading_shwr_hits_2shwr_numu_cc,
		                                                           h_leading_shwr_hits_2shwr_numu_cc_mixed,
		                                                           h_leading_shwr_hits_2shwr_nc,
		                                                           h_leading_shwr_hits_2shwr_nc_pi0,
		                                                           h_leading_shwr_hits_2shwr_cosmic,
		                                                           h_leading_shwr_hits_2shwr_other_mixed,
		                                                           h_leading_shwr_hits_2shwr_unmatched);

		_functions_instance.selection_functions::LeadingThetaPhi(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                         h_ele_theta_phi_nue_cc,
		                                                         h_ele_theta_phi_nue_cc_out_fv,
		                                                         h_ele_theta_phi_nue_cc_mixed,
		                                                         h_ele_theta_phi_numu_cc,
		                                                         h_ele_theta_phi_numu_cc_mixed,
		                                                         h_ele_theta_phi_nc,
		                                                         h_ele_theta_phi_nc_pi0,
		                                                         h_ele_theta_phi_cosmic,
		                                                         h_ele_theta_phi_other_mixed,
		                                                         h_ele_theta_phi_unmatched);

		_functions_instance.selection_functions::EnergyCosTheta(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                        h_ele_eng_costheta_nue_cc,
		                                                        h_ele_eng_costheta_nue_cc_out_fv,
		                                                        h_ele_eng_costheta_nue_cc_mixed,
		                                                        h_ele_eng_costheta_numu_cc,
		                                                        h_ele_eng_costheta_numu_cc_mixed,
		                                                        h_ele_eng_costheta_nc,
		                                                        h_ele_eng_costheta_nc_pi0,
		                                                        h_ele_eng_costheta_cosmic,
		                                                        h_ele_eng_costheta_other_mixed,
		                                                        h_ele_eng_costheta_unmatched);

		_functions_instance.selection_functions::PostCutsLeadingMomentumThetaSlice(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                           h_ele_momentum_slice_1_nue_cc,
		                                                                           h_ele_momentum_slice_1_nue_cc_out_fv,
		                                                                           h_ele_momentum_slice_1_nue_cc_mixed,
		                                                                           h_ele_momentum_slice_1_numu_cc,
		                                                                           h_ele_momentum_slice_1_numu_cc_mixed,
		                                                                           h_ele_momentum_slice_1_nc,
		                                                                           h_ele_momentum_slice_1_nc_pi0,
		                                                                           h_ele_momentum_slice_1_cosmic,
		                                                                           h_ele_momentum_slice_1_other_mixed,
		                                                                           h_ele_momentum_slice_1_unmatched,
		                                                                           h_ele_momentum_slice_2_nue_cc,
		                                                                           h_ele_momentum_slice_2_nue_cc_out_fv,
		                                                                           h_ele_momentum_slice_2_nue_cc_mixed,
		                                                                           h_ele_momentum_slice_2_numu_cc,
		                                                                           h_ele_momentum_slice_2_numu_cc_mixed,
		                                                                           h_ele_momentum_slice_2_nc,
		                                                                           h_ele_momentum_slice_2_nc_pi0,
		                                                                           h_ele_momentum_slice_2_cosmic,
		                                                                           h_ele_momentum_slice_2_other_mixed,
		                                                                           h_ele_momentum_slice_2_unmatched,
		                                                                           h_ele_momentum_slice_3_nue_cc,
		                                                                           h_ele_momentum_slice_3_nue_cc_out_fv,
		                                                                           h_ele_momentum_slice_3_nue_cc_mixed,
		                                                                           h_ele_momentum_slice_3_numu_cc,
		                                                                           h_ele_momentum_slice_3_numu_cc_mixed,
		                                                                           h_ele_momentum_slice_3_nc,
		                                                                           h_ele_momentum_slice_3_nc_pi0,
		                                                                           h_ele_momentum_slice_3_cosmic,
		                                                                           h_ele_momentum_slice_3_other_mixed,
		                                                                           h_ele_momentum_slice_3_unmatched);

		// _functions_instance.selection_functions::PostCutsdedxThetaSlice(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		//                                                                 h_dedx_slice_1_nue_cc,
		//                                                                 h_dedx_slice_1_nue_cc_out_fv,
		//                                                                 h_dedx_slice_1_nue_cc_mixed,
		//                                                                 h_dedx_slice_1_numu_cc,
		//                                                                 h_dedx_slice_1_numu_cc_mixed,
		//                                                                 h_dedx_slice_1_nc,
		//                                                                 h_dedx_slice_1_nc_pi0,
		//                                                                 h_dedx_slice_1_cosmic,
		//                                                                 h_dedx_slice_1_other_mixed,
		//                                                                 h_dedx_slice_1_unmatched,
		//                                                                 h_dedx_slice_2_nue_cc,
		//                                                                 h_dedx_slice_2_nue_cc_out_fv,
		//                                                                 h_dedx_slice_2_nue_cc_mixed,
		//                                                                 h_dedx_slice_2_numu_cc,
		//                                                                 h_dedx_slice_2_numu_cc_mixed,
		//                                                                 h_dedx_slice_2_nc,
		//                                                                 h_dedx_slice_2_nc_pi0,
		//                                                                 h_dedx_slice_2_cosmic,
		//                                                                 h_dedx_slice_2_other_mixed,
		//                                                                 h_dedx_slice_2_unmatched,
		//                                                                 h_dedx_slice_3_nue_cc,
		//                                                                 h_dedx_slice_3_nue_cc_out_fv,
		//                                                                 h_dedx_slice_3_nue_cc_mixed,
		//                                                                 h_dedx_slice_3_numu_cc,
		//                                                                 h_dedx_slice_3_numu_cc_mixed,
		//                                                                 h_dedx_slice_3_nc,
		//                                                                 h_dedx_slice_3_nc_pi0,
		//                                                                 h_dedx_slice_3_cosmic,
		//                                                                 h_dedx_slice_3_other_mixed,
		//                                                                 h_dedx_slice_3_unmatched);

		_functions_instance.selection_functions::EnergyCosThetaSlices(tpc_object_container_v, passed_tpco, _verbose,
		                                                              0, 0, tpco_classifier_v,
		                                                              h_ele_eng_for_nue_cc,
		                                                              h_ele_eng_for_nue_cc_out_fv,
		                                                              h_ele_eng_for_nue_cc_mixed,
		                                                              h_ele_eng_for_numu_cc,
		                                                              h_ele_eng_for_numu_cc_mixed,
		                                                              h_ele_eng_for_nc,
		                                                              h_ele_eng_for_nc_pi0,
		                                                              h_ele_eng_for_cosmic,
		                                                              h_ele_eng_for_other_mixed,
		                                                              h_ele_eng_for_unmatched,
		                                                              h_ele_eng_mid_nue_cc,
		                                                              h_ele_eng_mid_nue_cc_out_fv,
		                                                              h_ele_eng_mid_nue_cc_mixed,
		                                                              h_ele_eng_mid_numu_cc,
		                                                              h_ele_eng_mid_numu_cc_mixed,
		                                                              h_ele_eng_mid_nc,
		                                                              h_ele_eng_mid_nc_pi0,
		                                                              h_ele_eng_mid_cosmic,
		                                                              h_ele_eng_mid_other_mixed,
		                                                              h_ele_eng_mid_unmatched,
		                                                              h_ele_eng_back_nue_cc,
		                                                              h_ele_eng_back_nue_cc_out_fv,
		                                                              h_ele_eng_back_nue_cc_mixed,
		                                                              h_ele_eng_back_numu_cc,
		                                                              h_ele_eng_back_numu_cc_mixed,
		                                                              h_ele_eng_back_nc,
		                                                              h_ele_eng_back_nc_pi0,
		                                                              h_ele_eng_back_cosmic,
		                                                              h_ele_eng_back_other_mixed,
		                                                              h_ele_eng_back_unmatched);
		_functions_instance.selection_functions::EnergyCosThetaSlices(tpc_object_container_v, passed_tpco, _verbose,
		                                                              theta_translation, phi_translation, tpco_classifier_v,
		                                                              h_ele_eng_for_trans_nue_cc,
		                                                              h_ele_eng_for_trans_nue_cc_out_fv,
		                                                              h_ele_eng_for_trans_nue_cc_mixed,
		                                                              h_ele_eng_for_trans_numu_cc,
		                                                              h_ele_eng_for_trans_numu_cc_mixed,
		                                                              h_ele_eng_for_trans_nc,
		                                                              h_ele_eng_for_trans_nc_pi0,
		                                                              h_ele_eng_for_trans_cosmic,
		                                                              h_ele_eng_for_trans_other_mixed,
		                                                              h_ele_eng_for_trans_unmatched,
		                                                              h_ele_eng_mid_trans_nue_cc,
		                                                              h_ele_eng_mid_trans_nue_cc_out_fv,
		                                                              h_ele_eng_mid_trans_nue_cc_mixed,
		                                                              h_ele_eng_mid_trans_numu_cc,
		                                                              h_ele_eng_mid_trans_numu_cc_mixed,
		                                                              h_ele_eng_mid_trans_nc,
		                                                              h_ele_eng_mid_trans_nc_pi0,
		                                                              h_ele_eng_mid_trans_cosmic,
		                                                              h_ele_eng_mid_trans_other_mixed,
		                                                              h_ele_eng_mid_trans_unmatched,
		                                                              h_ele_eng_back_trans_nue_cc,
		                                                              h_ele_eng_back_trans_nue_cc_out_fv,
		                                                              h_ele_eng_back_trans_nue_cc_mixed,
		                                                              h_ele_eng_back_trans_numu_cc,
		                                                              h_ele_eng_back_trans_numu_cc_mixed,
		                                                              h_ele_eng_back_trans_nc,
		                                                              h_ele_eng_back_trans_nc_pi0,
		                                                              h_ele_eng_back_trans_cosmic,
		                                                              h_ele_eng_back_trans_other_mixed,
		                                                              h_ele_eng_back_trans_unmatched);

		_functions_instance.selection_functions::DifferentialEnergySlices(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                                  h_low_true_momentum, h_med_true_momentum, h_high_true_momentum);

	}//end event loop
	std::cout << "Debugging Events: " << num_debug_events << std::endl;

	std::cout << "------------------ " << std::endl;
	std::cout << " MC Nue          : " << total_mc_entries_inFV_nue << std::endl;
	std::cout << " MC NueBar       : " << total_mc_entries_inFV_nue_bar << std::endl;
	std::cout << " MC Entries in FV: " << total_mc_entries_inFV << std::endl;
	std::cout << total_mc_entries_inFV_numu_cc << std::endl;
	std::cout << total_mc_entries_inFV_nue_nc << std::endl;
	std::cout << total_mc_entries_inFV_numu_nc << std::endl;
	std::cout << total_mc_entries_inFV_numu_cc_bar << std::endl;
	std::cout << total_mc_entries_inFV_nue_nc_bar << std::endl;
	std::cout << total_mc_entries_inFV_numu_nc_bar << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "  End Selection    " << std::endl;
	std::cout << "------------------ " << std::endl;


	std::cout << test_mc_nue_cc_counter << std::endl;
	std::cout << test_mc_numu_cc_counter << std::endl;
	std::cout << test_mc_nue_nc_counter << std::endl;
	std::cout << test_mc_numu_nc_counter << std::endl;
	std::cout << test_mc_nue_cc_counter_bar << std::endl;
	std::cout << test_mc_numu_cc_counter_bar << std::endl;
	std::cout << test_mc_nue_nc_counter_bar << std::endl;
	std::cout << test_mc_numu_nc_counter_bar << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	std::cout << infv_test_mc_nue_cc_counter << std::endl;
	std::cout << infv_test_mc_numu_cc_counter << std::endl;
	std::cout << infv_test_mc_nue_nc_counter << std::endl;
	std::cout << infv_test_mc_numu_nc_counter << std::endl;
	std::cout << infv_test_mc_nue_cc_counter_bar << std::endl;
	std::cout << infv_test_mc_numu_cc_counter_bar << std::endl;
	std::cout << infv_test_mc_nue_nc_counter_bar << std::endl;
	std::cout << infv_test_mc_numu_nc_counter_bar << std::endl;

	//duplicate test here
	//first need to sort by the run number
	// std::sort(begin(duplicate_test), end(duplicate_test), [](std::tuple<int, int, int> const &t1, std::tuple<int,int,int> const &t2) {
	//              return std::get<0>(t1) < std::get<0>(t2);
	//      });
	//
	// std::tuple<int, int, int> last_tuple (0,0,0);
	// for(auto const tuple : duplicate_test)
	// {
	//      if(std::get<0>(tuple) == std::get<0>(last_tuple))
	//      {
	//              if(std::get<1>(tuple) == std::get<1>(last_tuple))
	//              {
	//                      if(std::get<2>(tuple) == std::get<2>(last_tuple))
	//                      {
	//                              std::cout << "Duplicate Events at Num: " << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ", " << std::get<2>(tuple) << std::endl;
	//                      }
	//              }
	//      }
	//      last_tuple = tuple;
	// }

	//we also want some metrics to print at the end
	//*************************************************************************************************************************
	//*************************************************************************************************************************
	selection_functions::PrintInfo( total_mc_entries_inFV, in_time_counter_v, intime_in_time_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "In Time");
	selection_functions_data::PrintInfoData(1 * data_in_time_counter_v->at(0),                                  "In Time");
	selection_functions::PrintInfo( total_mc_entries_inFV, pe_counter_v, intime_pe_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "PE Threshold");
	selection_functions_data::PrintInfoData(1 * data_pe_counter_v->at(0),                                       "PE Threshold");
	selection_functions::PrintInfo( total_mc_entries_inFV, reco_nue_counter_v, intime_reco_nue_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "Reco Nue");
	selection_functions_data::PrintInfoData(1 * data_reco_nue_counter_v->at(0),                                 "Reco Nue");
	selection_functions::PrintInfo( total_mc_entries_inFV, in_fv_counter_v, intime_in_fv_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "In FV");
	selection_functions_data::PrintInfoData(1 * data_in_fv_counter_v->at(0),                                    "In FV");
	selection_functions::PrintInfo( total_mc_entries_inFV, vtx_flash_counter_v, intime_vtx_flash_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "Vtx-to-Flash");
	selection_functions_data::PrintInfoData(1 * data_vtx_flash_counter_v->at(0),                                "Vtx-to-Flash");
	selection_functions::PrintInfo( total_mc_entries_inFV, shwr_tpco_counter_v, intime_shwr_tpco_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "Shower-to-TPCO");
	selection_functions_data::PrintInfoData(1 * data_shwr_tpco_counter_v->at(0),                                "Shower-to-TPCO");
	selection_functions::PrintInfo( total_mc_entries_inFV, trk_tpco_counter_v, intime_trk_tpco_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "Track-to-TPCO");
	selection_functions_data::PrintInfoData(1 * data_trk_tpco_counter_v->at(0),                                 "Track-to-TPCO");
	selection_functions::PrintInfo( total_mc_entries_inFV, hit_threshold_counter_v, intime_hit_threshold_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "Hit Threshold");
	selection_functions_data::PrintInfoData(1 * data_hit_threshold_counter_v->at(0),                            "Hit Threshold");
	selection_functions::PrintInfo( total_mc_entries_inFV, hit_threshold_collection_counter_v, intime_hit_threshold_collection_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                      "YPlane Hit Threshold");
	selection_functions_data::PrintInfoData(1 * data_hit_threshold_collection_counter_v->at(0),                  "YPlane hit Threshold");
	selection_functions::PrintInfo( total_mc_entries_inFV, open_angle_counter_v, intime_open_angle_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "Open Angle");
	selection_functions_data::PrintInfoData(1 * data_open_angle_counter_v->at(0),                               "Open Angle");
	selection_functions::PrintInfo( total_mc_entries_inFV, dedx_counter_v, intime_dedx_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     " dE / dx ");
	selection_functions_data::PrintInfoData(1 * data_dedx_counter_v->at(0),                                     " dE / dx ");
	selection_functions::PrintInfo( total_mc_entries_inFV, secondary_shower_counter_v, intime_secondary_shower_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     ">1 Shower TPCO Dist");
	selection_functions_data::PrintInfoData(1 * data_secondary_shower_counter_v->at(0),                         ">1 Shower TPCO Dist");
	selection_functions::PrintInfo( total_mc_entries_inFV, hit_lengthRatio_counter_v, intime_hit_lengthRatio_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "Hit Length Ratio");
	selection_functions_data::PrintInfoData(1 * data_hit_lengthRatio_counter_v->at(0),                          "Hit Length Ratio");
	selection_functions::PrintInfo( total_mc_entries_inFV, trk_len_shwr_len_ratio_counter_v, intime_trk_len_shwr_len_ratio_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "TrkLen/ShwrLen Ratio");
	selection_functions_data::PrintInfoData(1 * data_trk_len_shwr_len_ratio_counter_v->at(0),                   "TrkLen/ShwrLen Ratio");
	selection_functions::PrintInfo( total_mc_entries_inFV, track_containment_counter_v, intime_track_containment_counter_v->at(0),
	                                intime_scale_factor, data_scale_factor,                                     "Track Containment");
	selection_functions_data::PrintInfoData(1 * data_track_containment_counter_v->at(0),                        "Track Containment");

	//***********************************************************************************************************************
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, in_time_counter_v, intime_in_time_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "In Time", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, pe_counter_v, intime_pe_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "PE Threshold", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, reco_nue_counter_v, intime_reco_nue_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Reco Nue", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, in_fv_counter_v, intime_in_fv_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "In FV", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, vtx_flash_counter_v, intime_vtx_flash_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Vtx-to-Flash", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, shwr_tpco_counter_v, intime_shwr_tpco_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Shwr-to-TPCO", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, trk_tpco_counter_v, intime_trk_tpco_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Trk-to-TPCO", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, hit_threshold_counter_v, intime_hit_threshold_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Hit Threshold", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, hit_threshold_collection_counter_v, intime_hit_threshold_collection_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "YPlane Hit Threshold", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, open_angle_counter_v, intime_open_angle_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Open Angle", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, dedx_counter_v, intime_dedx_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "dE / dx", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, secondary_shower_counter_v, intime_secondary_shower_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, ">1 Shower-Dist", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, hit_lengthRatio_counter_v, intime_hit_lengthRatio_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Hit-Len-Ratio", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, trk_len_shwr_len_ratio_counter_v, intime_trk_len_shwr_len_ratio_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "TrkLen/ShwrLen Ratio", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, track_containment_counter_v, intime_track_containment_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Track Contained", results_v);

	selection_functions::PrintTopologyPurity(no_track, has_track, _1_shwr, _2_shwr, _3_shwr, _4_shwr);
	//*************************************************************************************************************************
	//*************************************************************************************************************************
	selection_functions::XSecWork(track_containment_counter_v->at(7),
	                              track_containment_counter_v->at(22), track_containment_counter_v->at(23), track_containment_counter_v->at(1),
	                              track_containment_counter_v->at(9), track_containment_counter_v->at(2),
	                              track_containment_counter_v->at(3), track_containment_counter_v->at(4),
	                              track_containment_counter_v->at(11), track_containment_counter_v->at(10), track_containment_counter_v->at(5),
	                              track_containment_counter_v->at(6), intime_track_containment_counter_v->at(0),
	                              intime_scale_factor, data_track_containment_counter_v->at(0), data_scale_factor,
	                              fv_boundary_v, flux_nue, flux_nue_bar, selected_energy_vector, genie_xsec_nue, genie_xsec_nue_bar,
	                              total_mc_entries_inFV_nue, total_mc_entries_inFV_nue_bar);
	//*************************************************************************************************************************
	//*************************************************************************************************************************

	if(_post_cuts_verbose == true)
	{
		std::cout << "Print Post Cuts MC: " << std::endl;
		_functions_instance.selection_functions::PrintPostCutVector(post_cuts_v,      _post_cuts_verbose);
		//filling a file with all of the run subrun and event numbers
	}
	if(_post_cuts_verbose == true)
	{
		std::cout << "Print Post Cuts On-Beam Data: " << std::endl;
		_functions_instance.selection_functions::PrintPostCutVector(post_cuts_v_data, _post_cuts_verbose);
	}
	std::ofstream selected_run_subrun_event_file;
	selected_run_subrun_event_file.open("selected_run_subrun_event_list.txt");
	for(auto const post_cuts : * post_cuts_v)
	{
		const int event_num = std::get<0>(post_cuts);
		const int run_num = std::get<1>(post_cuts);
		const int sub_run_num = std::get<2>(post_cuts);
		const std::string tpco_classification = std::get<7>(post_cuts);
		selected_run_subrun_event_file << tpco_classification << " " << run_num << " " << sub_run_num << " " << event_num << "\n";
	}
	selected_run_subrun_event_file.close();
	_functions_instance.selection_functions::PostCutVectorPlots(post_cuts_v, _post_cuts_verbose,
	                                                            h_post_cuts_num_showers_purity_qe,
	                                                            h_post_cuts_num_showers_purity_res,
	                                                            h_post_cuts_num_showers_purity_dis,
	                                                            h_post_cuts_num_showers_purity_coh,
	                                                            h_post_cuts_num_showers_purity_mec);
	_functions_instance.selection_functions::PostCutVectorPlots(post_open_angle_cuts_v, _post_cuts_verbose,
	                                                            h_post_open_angle_cuts_num_showers_purity_qe,
	                                                            h_post_open_angle_cuts_num_showers_purity_res,
	                                                            h_post_open_angle_cuts_num_showers_purity_dis,
	                                                            h_post_open_angle_cuts_num_showers_purity_coh,
	                                                            h_post_open_angle_cuts_num_showers_purity_mec);

	_functions_instance.selection_functions::PostCutVector2DPlots(post_cuts_v, _post_cuts_verbose, intime_scale_factor, data_scale_factor,
	                                                              total_mc_entries_inFV,
	                                                              h_post_cuts_num_tracks_showers_purity_qe,
	                                                              h_post_cuts_num_tracks_showers_purity_res,
	                                                              h_post_cuts_num_tracks_showers_purity_dis,
	                                                              h_post_cuts_num_tracks_showers_purity_coh,
	                                                              h_post_cuts_num_tracks_showers_purity_mec,
	                                                              h_post_cuts_num_tracks_showers_purity_total,
	                                                              h_post_cuts_num_tracks_showers_signal_total,
	                                                              h_post_cuts_num_tracks_showers_bkg_total,
	                                                              h_post_cuts_num_tracks_showers_total_total);

	std::cout << "*************************" << std::endl;
	std::cout << "Dirt Test? (Scaled to On Beam)" << std::endl;
	std::cout << "Near MC XYZ   : " << xyz_near_mc * data_scale_factor << std::endl;
	std::cout << "Near EXT XYZ  : " << xyz_near_ext * intime_scale_factor << std::endl;
	std::cout << "Total (MC+EXT): " << (xyz_near_mc + xyz_near_ext) << " | (" << xyz_near_mc << " , " << xyz_near_ext << ")" << std::endl;
	std::cout << "Near Data XYZ : " << xyz_near_data << std::endl;
	std::cout << "=========================" << std::endl;
	std::cout << "Far MC XYZ    : "  << xyz_far_mc * data_scale_factor << std::endl;
	std::cout << "Far EXT XYZ   : "  << xyz_far_ext * intime_scale_factor << std::endl;
	std::cout << "Total (MC+EXT): "  << (xyz_far_mc + xyz_far_ext) << " | (" << xyz_far_mc << " , " << xyz_far_ext << ")" << std::endl;
	std::cout << "Far Data XYZ  : "  << xyz_far_data << std::endl;
	std::cout << "*************************" << std::endl;


//********************//
//**** Histograms ****//
//*******************//
	for(auto const flash_timing : * flash_time)        {h_flash_time->Fill(flash_timing.first);        }
	for(auto const flash_timing : * data_flash_time)   {h_flash_time_data->Fill(flash_timing.first);   }
	for(auto const flash_timing : * intime_flash_time) {h_flash_time_intime->Fill(flash_timing.first); }
	histogram_functions::Plot1DHistogram(h_flash_time,        "Flash Time [#mus]", Form("%s%s", file_locate_prefix, "flash_time.pdf"));
	histogram_functions::Plot1DHistogram(h_flash_time_intime, "Flash Time [#mus]", Form("%s%s", file_locate_prefix, "flash_time_intime.pdf"));
	histogram_functions::Plot1DHistogram(h_flash_time_data,   "Flash Time [#mus]", Form("%s%s", file_locate_prefix, "flash_time_data.pdf"));

	histogram_functions::PlotFlashInfo(h_flash_z_mc, h_flash_z_intime, h_flash_z_mc, intime_scale_factor, data_scale_factor, "Largest Flash Z [cm]",
	                                   Form("%s%s", file_locate_prefix, "flash_z_data.pdf"));

	histogram_functions::TimingHistograms(h_flash_time, h_flash_time_intime, h_flash_time_data, data_scale_factor, intime_scale_factor,
	                                      "Flash Time [#mus]", Form("%s%s", file_locate_prefix, "flash_time_data_subtraction_mc.pdf"));
	histogram_functions::TimingHistogramsOverlay(data_flash_time, h_flash_time_intime, h_flash_time_data, intime_scale_factor, "Flash Time [#mus]",
	                                             Form("%s%s", file_locate_prefix, "flash_time_data_overlay.pdf"),
	                                             Form("%s%s", file_locate_prefix, "flash_time_data_subtraction.pdf"));

	histogram_functions::Plot1DHistogram (h_nue_eng_eff_num, "True Neutrino Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_true_neutrino_energy.pdf"));
	histogram_functions::Plot1DHistogram (h_nue_eng_eff_den, "True Neutrino Energy [GeV]", Form("%s%s", file_locate_prefix, "all_true_neutrino_energy.pdf"));
	histogram_functions::Plot1DHistogram (h_nue_num_part_eff_den, "True Particle Multiplicity",
	                                      Form("%s%s", file_locate_prefix, "all_true_neutrino_num_particles.pdf"));
	histogram_functions::Plot1DHistogram (h_nue_num_part_eff_num, "Selected True Particle Multiplicity",
	                                      Form("%s%s", file_locate_prefix, "selected_true_neutrino_num_particles.pdf"));
	histogram_functions::Plot1DHistogram (h_nue_num_chrg_part_eff_den, "True Charged Particle Multiplicity",
	                                      Form("%s%s", file_locate_prefix, "all_true_neutrino_num_charged_particles.pdf"));
	histogram_functions::Plot1DHistogram (h_nue_num_chrg_part_eff_num, "Selected Charged True Particle Multiplicity",
	                                      Form("%s%s", file_locate_prefix, "selected_true_neutrino_num_charged_particles.pdf"));
	histogram_functions::Plot1DHistogram(h_nue_true_theta, "True Neutrino Theta [Degrees]", Form("%s%s", file_locate_prefix, "true_nue_theta.pdf"));
	histogram_functions::Plot1DHistogram(h_nue_true_phi,   "True Neutrino Phi [Degrees]",   Form("%s%s", file_locate_prefix, "true_nue_phi.pdf"));
	histogram_functions::Plot2DHistogramLogZ(h_nue_true_theta_phi, "True Nue CC", "Phi [Degrees]", "Theta [Degrees]",
	                                         Form("%s%s", file_locate_prefix, "true_nue_theta_phi.pdf"), "colz");
	histogram_functions::Plot2DHistogram(h_all_true_energy_theta, "True All Neutrino Flavours", "True Neutrino Energy [GeV]", "True Neutrino Theta [Degrees]",
	                                     Form("%s%s", file_locate_prefix, "true_all_energy_theta.pdf"));
	histogram_functions::Plot2DHistogram(h_nue_true_energy_theta, "True Nue/Nue-bar CC", "True Neutrino Energy [GeV]", "True Neutrino Theta [Degrees]",
	                                     Form("%s%s", file_locate_prefix, "true_nue_energy_theta.pdf"));
	histogram_functions::Plot2DHistogram(h_nue_true_energy_phi, "True Nue/Nue-bar CC", "True Neutrino Energy [GeV]", "True Neutrino Phi [Degrees]",
	                                     Form("%s%s", file_locate_prefix, "true_nue_energy_phi.pdf"));
	histogram_functions::Plot2DHistogram(h_ele_true_energy_theta, "True Electron", "True Electron Energy [GeV]", "True Electron Theta [Degrees]",
	                                     Form("%s%s", file_locate_prefix, "true_ele_energy_theta.pdf"));
	histogram_functions::Plot2DHistogram(h_ele_true_energy_phi, "True Electron", "True Electron Energy [GeV]", "True Electron Phi [Degrees]",
	                                     Form("%s%s", file_locate_prefix, "true_ele_energy_phi.pdf"));

	histogram_functions::Plot1DHistogram(h_ele_phi_eff_num, "Selected True Electron Phi [Degrees]",
	                                     Form("%s%s", file_locate_prefix, "selected_true_electron_phi.pdf"));
	histogram_functions::Plot1DHistogram(h_ele_phi_eff_den, "True Electron Phi [Degrees]",
	                                     Form("%s%s", file_locate_prefix, "true_electron_phi.pdf"));
	histogram_functions::Plot1DHistogram(h_ele_cos_theta_eff_num, "Selected True Electron Cos(#theta)",
	                                     Form("%s%s", file_locate_prefix, "selected_true_electron_cos_theta.pdf"));
	histogram_functions::Plot1DHistogram(h_ele_cos_theta_eff_den, "True Electron Cos(#theta)",
	                                     Form("%s%s", file_locate_prefix, "true_electron_cos_theta.pdf"));
	histogram_functions::Plot1DHistogram(h_ele_theta_eff_num, "Selected True Electron Theta [Degrees]",
	                                     Form("%s%s", file_locate_prefix, "selected_true_electron_theta.pdf"));
	histogram_functions::Plot1DHistogram(h_ele_theta_eff_den, "True Electron Theta [Degrees]",
	                                     Form("%s%s", file_locate_prefix, "true_electron_theta.pdf"));

	histogram_functions::PlotTEfficiency (h_nue_eng_eff_num, h_nue_eng_eff_den,
	                                      ";True Neutrino Energy [GeV];Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_energy_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_eng_eff_num, h_nue_eng_eff_den, true,
	                                      ";True Neutrino Energy [GeV];Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_energy_efficiency_rebin.pdf"));
	histogram_functions::PlotTEfficiencyOverlay(h_ele_eng_eff_num,
	                                            h_ele_eng_eff_reco_nue,
	                                            h_ele_eng_eff_in_fv,
	                                            h_ele_eng_eff_vtx_flash,
	                                            h_ele_eng_eff_shwr_vtx,
	                                            h_ele_eng_eff_trk_vtx,
	                                            h_ele_eng_eff_hit,
	                                            h_ele_eng_eff_yhit,
	                                            h_ele_eng_eff_open_angle,
	                                            h_ele_eng_eff_dedx,
	                                            h_ele_eng_eff_2shwr,
	                                            h_ele_eng_eff_hit_len,
	                                            h_ele_eng_eff_trk_shwr,
	                                            h_ele_eng_eff_contain,
	                                            h_ele_eng_eff_den, true, ";True Electron Energy [GeV]",
	                                            Form("%s%s", file_locate_prefix, "signal_selection_ele_energy_efficiency_stack_rebin.pdf"));

	histogram_functions::PlotTEfficiency (h_nue_vtx_x_eff_num, h_nue_vtx_x_eff_den,
	                                      ";True Neutrino Vtx X [cm];Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_vtx_x_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_vtx_y_eff_num, h_nue_vtx_y_eff_den,
	                                      ";True Neutrino Vtx Y [cm];Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_vtx_y_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_vtx_z_eff_num, h_nue_vtx_z_eff_den,
	                                      ";True Neutrino Vtx Z [cm];Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_vtx_z_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_dir_x_eff_num, h_nue_dir_x_eff_den,
	                                      ";True Neutrino Dir X;Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_dir_x_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_dir_y_eff_num, h_nue_dir_y_eff_den,
	                                      ";True Neutrino Dir Y;Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_dir_y_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_dir_z_eff_num, h_nue_dir_z_eff_den,
	                                      ";True Neutrino Dir Z;Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_dir_z_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_num_part_eff_num, h_nue_num_part_eff_den,
	                                      ";True Particle Multiplicity;Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_num_particles_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_num_chrg_part_eff_num, h_nue_num_chrg_part_eff_den, ";True Charged Particle Multiplicity;Efficiency",
	                                      Form("%s%s", file_locate_prefix, "signal_selection_nu_num_charged_particles_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_cos_theta_eff_num, h_nue_cos_theta_eff_den,
	                                      ";True Neutrino Cos(#theta);Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_nu_cos_theta_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_nue_phi_eff_num, h_nue_phi_eff_den, ";True Neutrino Phi;Efficiency",
	                                      Form("%s%s", file_locate_prefix, "signal_selection_nu_phi_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_ele_cos_theta_eff_num, h_ele_cos_theta_eff_den, ";True Nue Electron Cos(#theta);Efficiency",
	                                      Form("%s%s", file_locate_prefix, "signal_selection_ele_cos_theta_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_ele_phi_eff_num, h_ele_phi_eff_den, ";True Nue Electron Phi;Efficiency",
	                                      Form("%s%s", file_locate_prefix, "signal_selection_ele_phi_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_ele_dir_x_eff_num, h_ele_dir_x_eff_den,
	                                      ";True Nue Electron Dir X;Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_ele_dir_x_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_ele_dir_y_eff_num, h_ele_dir_y_eff_den,
	                                      ";True Nue Electron Dir Y;Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_ele_dir_y_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_ele_dir_z_eff_num, h_ele_dir_z_eff_den,
	                                      ";True Nue Electron Dir Z;Efficiency", Form("%s%s", file_locate_prefix, "signal_selection_ele_dir_z_efficiency.pdf"));
	// histogram_functions::PlotTEfficiency (h_post_cuts_num_tracks_showers_signal_total, h_post_cuts_num_tracks_showers_bkg_total,
	//                                       ";Num Reco Showers; Num Reco Tracks",
	//Form("%s%s", file_locate_prefix, "post_cuts_num_tracks_showers_eff_purity_total.pdf"));
	histogram_functions::PlotTEfficiency(h_ele_cos_theta_eff_num_pre_cuts, h_ele_cos_theta_eff_den, ";True Electron Cos(#theta) Pre-Cuts;Efficiency",
	                                     Form("%s%s", file_locate_prefix, "signal_selection_ele_cos_theta_pre_cuts_efficiency.pdf"));

	histogram_functions::PlotTEfficiency(h_ele_eng_eff_num_pre_cuts, h_ele_eng_eff_den,
	                                     ";True Electron Energy Pre-Cuts;Efficiency",
	                                     Form("%s%s", file_locate_prefix, "signal_selection_ele_eng_pre_cuts_efficiency.pdf"));
	histogram_functions::PlotTEfficiency (h_ele_eng_eff_contain, h_ele_eng_eff_den, ";True Electron Energy [GeV];Efficiency",
	                                      Form("%s%s", file_locate_prefix, "signal_selection_ele_energy_efficiency.pdf"));
	//these histograms are the same, but are being rebinned
	histogram_functions::PlotTEfficiency(h_ele_eng_eff_num_pre_cuts, h_ele_eng_eff_den, true,
	                                     ";True Electron Energy Pre-Cuts [GeV];Efficiency",
	                                     Form("%s%s", file_locate_prefix, "signal_selection_ele_eng_pre_cuts_efficiency_rebin.pdf"));
	histogram_functions::PlotTEfficiency (h_ele_eng_eff_contain, h_ele_eng_eff_den, true,
	                                      ";True Electron Energy [GeV];Efficiency",
	                                      Form("%s%s", file_locate_prefix, "signal_selection_ele_energy_efficiency_rebin.pdf"));


	histogram_functions::Plot2DHistogram (h_tracks_showers, "Post Cuts - Showers/Tracks per Candidate Nue TPC Object",
	                                      "Reco Tracks", "Reco Showers", Form("%s%s", file_locate_prefix, "post_cuts_showers_tracks.pdf"));
	histogram_functions::Plot2DHistogram (h_tracks_showers, "Post Cuts - Showers/Tracks per Candidate Nue TPC Object",
	                                      "Reco Tracks", "Reco Showers", Form("%s%s", file_locate_prefix, "post_cuts_showers_tracks_cosmic.pdf"));
	histogram_functions::Plot2DHistogram (h_tracks_showers, "Post Cuts - Showers/Tracks per Candidate Nue TPC Object",
	                                      "Reco Tracks", "Reco Showers", Form("%s%s", file_locate_prefix, "post_cuts_showers_tracks_numu.pdf"));

	histogram_functions::PlotSimpleStack (h_leading_shower_open_angle_nue_cc,  h_leading_shower_open_angle_nue_cc_mixed,
	                                      h_leading_shower_open_angle_nue_cc_out_fv,
	                                      h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_numu_cc_mixed,
	                                      h_leading_shower_open_angle_cosmic,  h_leading_shower_open_angle_nc,
	                                      h_leading_shower_open_angle_nc_pi0,  h_leading_shower_open_angle_other_mixed,
	                                      h_leading_shower_open_angle_unmatched, "",
	                                      "Shower Opening Angle [Degrees]", "", Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_leading_shower_open_angle_nue_cc,  h_leading_shower_open_angle_nue_cc_mixed,
	                                            h_leading_shower_open_angle_nue_cc_out_fv,
	                                            h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_numu_cc_mixed,
	                                            h_leading_shower_open_angle_cosmic,  h_leading_shower_open_angle_nc,
	                                            h_leading_shower_open_angle_nc_pi0,  h_leading_shower_open_angle_other_mixed,
	                                            h_leading_shower_open_angle_unmatched, h_leading_shower_open_angle_intime, intime_scale_factor, data_scale_factor,
	                                            "", "Shower Opening Angle [Degrees]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_leading_shower_open_angle_nue_cc,  h_leading_shower_open_angle_nue_cc_mixed,
	                                          h_leading_shower_open_angle_nue_cc_out_fv,
	                                          h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_numu_cc_mixed,
	                                          h_leading_shower_open_angle_cosmic,  h_leading_shower_open_angle_nc,
	                                          h_leading_shower_open_angle_nc_pi0,  h_leading_shower_open_angle_other_mixed,
	                                          h_leading_shower_open_angle_unmatched, h_leading_shower_open_angle_intime, intime_scale_factor,
	                                          h_leading_shower_open_angle_data, data_scale_factor,
	                                          "", "Shower Opening Angle [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_leading_shower_open_angle_nue_cc,  h_leading_shower_open_angle_nue_cc_mixed,
		                                          h_leading_shower_open_angle_nue_cc_out_fv,
		                                          h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_numu_cc_mixed,
		                                          h_leading_shower_open_angle_cosmic,  h_leading_shower_open_angle_nc,
		                                          h_leading_shower_open_angle_nc_pi0,  h_leading_shower_open_angle_other_mixed,
		                                          h_leading_shower_open_angle_unmatched, h_leading_shower_open_angle_intime, scaled_intime_scale_factor,
		                                          h_leading_shower_open_angle_data, data_scale_factor,
		                                          "", "Shower Opening Angle [Degrees] (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle_scaled_data.pdf"));
	}
	histogram_functions::PlotSimpleStackData (h_leading_shower_open_angle_nue_cc_after,  h_leading_shower_open_angle_nue_cc_mixed_after,
	                                          h_leading_shower_open_angle_nue_cc_out_fv_after,
	                                          h_leading_shower_open_angle_numu_cc_after, h_leading_shower_open_angle_numu_cc_mixed_after,
	                                          h_leading_shower_open_angle_cosmic_after,  h_leading_shower_open_angle_nc_after,
	                                          h_leading_shower_open_angle_nc_pi0_after,  h_leading_shower_open_angle_other_mixed_after,
	                                          h_leading_shower_open_angle_unmatched_after, h_leading_shower_open_angle_intime_after, intime_scale_factor,
	                                          h_leading_shower_open_angle_data_after, data_scale_factor,
	                                          "", "Shower Opening Angle [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle_data_after.pdf"));

	histogram_functions::PlotSimpleStack (h_leading_shower_open_angle_1_nue_cc,  h_leading_shower_open_angle_1_nue_cc_mixed,
	                                      h_leading_shower_open_angle_1_nue_cc_out_fv,
	                                      h_leading_shower_open_angle_1_numu_cc, h_leading_shower_open_angle_1_numu_cc_mixed,
	                                      h_leading_shower_open_angle_1_cosmic,  h_leading_shower_open_angle_1_nc,
	                                      h_leading_shower_open_angle_1_nc_pi0,  h_leading_shower_open_angle_1_other_mixed,
	                                      h_leading_shower_open_angle_1_unmatched, "",
	                                      "Shower Opening Angle [Degrees]", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle_1_shower.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_leading_shower_open_angle_1_nue_cc,  h_leading_shower_open_angle_1_nue_cc_mixed,
	                                            h_leading_shower_open_angle_1_nue_cc_out_fv,
	                                            h_leading_shower_open_angle_1_numu_cc, h_leading_shower_open_angle_1_numu_cc_mixed,
	                                            h_leading_shower_open_angle_1_cosmic,  h_leading_shower_open_angle_1_nc,
	                                            h_leading_shower_open_angle_1_nc_pi0,  h_leading_shower_open_angle_1_other_mixed,
	                                            h_leading_shower_open_angle_1_unmatched, h_leading_shower_open_angle_1_intime,
	                                            intime_scale_factor, data_scale_factor,
	                                            "", "Shower Opening Angle [Degrees]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle_1_shower_intime.pdf"));

	histogram_functions::PlotSimpleStack (h_leading_shower_open_angle_2plus_nue_cc,  h_leading_shower_open_angle_2plus_nue_cc_mixed,
	                                      h_leading_shower_open_angle_2plus_nue_cc_out_fv,
	                                      h_leading_shower_open_angle_2plus_numu_cc, h_leading_shower_open_angle_2plus_numu_cc_mixed,
	                                      h_leading_shower_open_angle_2plus_cosmic,  h_leading_shower_open_angle_2plus_nc,
	                                      h_leading_shower_open_angle_2plus_nc_pi0,  h_leading_shower_open_angle_2plus_other_mixed,
	                                      h_leading_shower_open_angle_2plus_unmatched, "",
	                                      "Shower Opening Angle [Degrees]", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle_2plus_showers.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_leading_shower_open_angle_2plus_nue_cc,  h_leading_shower_open_angle_2plus_nue_cc_mixed,
	                                            h_leading_shower_open_angle_2plus_nue_cc_out_fv,
	                                            h_leading_shower_open_angle_2plus_numu_cc, h_leading_shower_open_angle_2plus_numu_cc_mixed,
	                                            h_leading_shower_open_angle_2plus_cosmic,  h_leading_shower_open_angle_2plus_nc,
	                                            h_leading_shower_open_angle_2plus_nc_pi0,  h_leading_shower_open_angle_2plus_other_mixed,
	                                            h_leading_shower_open_angle_2plus_unmatched, h_leading_shower_open_angle_2plus_intime,
	                                            intime_scale_factor, data_scale_factor, "", "Shower Opening Angle [Degrees]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle_2plus_showers_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_leading_shower_open_angle_2plus_nue_cc,  h_leading_shower_open_angle_2plus_nue_cc_mixed,
	                                          h_leading_shower_open_angle_2plus_nue_cc_out_fv,
	                                          h_leading_shower_open_angle_2plus_numu_cc, h_leading_shower_open_angle_2plus_numu_cc_mixed,
	                                          h_leading_shower_open_angle_2plus_cosmic,  h_leading_shower_open_angle_2plus_nc,
	                                          h_leading_shower_open_angle_2plus_nc_pi0,  h_leading_shower_open_angle_2plus_other_mixed,
	                                          h_leading_shower_open_angle_2plus_unmatched, h_leading_shower_open_angle_2plus_intime,
	                                          intime_scale_factor, h_leading_shower_open_angle_2plus_data, data_scale_factor,
	                                          "", "Shower Opening Angle [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_shower_open_angle_2plus_showers_data.pdf"));

	histogram_functions::PlotSimpleStack (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                      h_dedx_cuts_nue_cc_out_fv,
	                                      h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                      h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                      h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                      h_dedx_cuts_unmatched, "",
	                                      "Collection Plane dE/dx [MeV/cm]", "", Form("%s%s", file_locate_prefix, "post_cuts_dedx_cuts.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                            h_dedx_cuts_nue_cc_out_fv,
	                                            h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                            h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                            h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                            h_dedx_cuts_unmatched, h_dedx_cuts_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Collection Plane dE/dx [MeV/cm]", "", Form("%s%s", file_locate_prefix, "post_cuts_dedx_cuts_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                          h_dedx_cuts_nue_cc_out_fv,
	                                          h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                          h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                          h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                          h_dedx_cuts_unmatched, h_dedx_cuts_intime, intime_scale_factor,
	                                          h_dedx_cuts_data, data_scale_factor, "",
	                                          "Collection Plane dE/dx [MeV/cm]", "", Form("%s%s", file_locate_prefix, "post_cuts_dedx_cuts_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
		                                          h_dedx_cuts_nue_cc_out_fv,
		                                          h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
		                                          h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
		                                          h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
		                                          h_dedx_cuts_unmatched, h_dedx_cuts_intime, scaled_intime_scale_factor,
		                                          h_dedx_cuts_data, data_scale_factor, "",
		                                          "Collection Plane dE/dx [MeV/cm] (Scaled)", "", Form("%s%s", file_locate_prefix, "post_cuts_dedx_cuts_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStackData (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                          h_dedx_cuts_nue_cc_out_fv,
	                                          h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                          h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                          h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                          h_dedx_cuts_unmatched, h_dedx_cuts_scale_1_intime, intime_scale_factor,
	                                          h_dedx_cuts_scale_1_data, data_scale_factor, "",
	                                          "Collection Plane dE/dx [MeV/cm]", "", Form("%s%s", file_locate_prefix, "post_cuts_dedx_cuts_scale_1_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                          h_dedx_cuts_nue_cc_out_fv,
	                                          h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                          h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                          h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                          h_dedx_cuts_unmatched, h_dedx_cuts_scale_2_intime, intime_scale_factor,
	                                          h_dedx_cuts_scale_2_data, data_scale_factor, "",
	                                          "Collection Plane dE/dx [MeV/cm]", "", Form("%s%s", file_locate_prefix, "post_cuts_dedx_cuts_scale_2_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                          h_dedx_cuts_nue_cc_out_fv,
	                                          h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                          h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                          h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                          h_dedx_cuts_unmatched, h_dedx_cuts_scale_3_intime, intime_scale_factor,
	                                          h_dedx_cuts_scale_3_data, data_scale_factor, "",
	                                          "Collection Plane dE/dx [MeV/cm]", "", Form("%s%s", file_locate_prefix, "post_cuts_dedx_cuts_scale_3_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                          h_dedx_cuts_nue_cc_out_fv,
	                                          h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                          h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                          h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                          h_dedx_cuts_unmatched, h_dedx_cuts_scale_3_intime, scaled_intime_scale_factor,
	                                          h_dedx_cuts_scale_3_data, data_scale_factor, "",
	                                          "Collection Plane dE/dx [MeV/cm (Scaled)]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_dedx_cuts_scale_3_scaled_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_dedx_cuts_nue_cc_after,  h_dedx_cuts_nue_cc_mixed_after,
	                                          h_dedx_cuts_nue_cc_out_fv_after,
	                                          h_dedx_cuts_numu_cc_after, h_dedx_cuts_numu_cc_mixed_after,
	                                          h_dedx_cuts_cosmic_after,  h_dedx_cuts_nc_after,
	                                          h_dedx_cuts_nc_pi0_after,  h_dedx_cuts_other_mixed_after,
	                                          h_dedx_cuts_unmatched_after, h_dedx_cuts_intime_after, intime_scale_factor,
	                                          h_dedx_cuts_data_after, data_scale_factor, "",
	                                          "Collection Plane dE/dx [MeV/cm]", "", Form("%s%s", file_locate_prefix, "post_cuts_dedx_cuts_data_after.pdf"));


	histogram_functions::PlotSimpleStack (h_vtx_flash_nue_cc,  h_vtx_flash_nue_cc_mixed,
	                                      h_vtx_flash_nue_cc_out_fv,
	                                      h_vtx_flash_numu_cc, h_vtx_flash_numu_cc_mixed,
	                                      h_vtx_flash_cosmic,  h_vtx_flash_nc,
	                                      h_vtx_flash_nc_pi0,  h_vtx_flash_other_mixed,
	                                      h_vtx_flash_unmatched, "",
	                                      "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_vtx_to_flash_distance.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_vtx_flash_nue_cc,  h_vtx_flash_nue_cc_mixed,
	                                            h_vtx_flash_nue_cc_out_fv,
	                                            h_vtx_flash_numu_cc, h_vtx_flash_numu_cc_mixed,
	                                            h_vtx_flash_cosmic,  h_vtx_flash_nc,
	                                            h_vtx_flash_nc_pi0,  h_vtx_flash_other_mixed,
	                                            h_vtx_flash_unmatched, h_vtx_flash_intime, intime_scale_factor, data_scale_factor, "",
	                                            "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_vtx_to_flash_distance_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_vtx_flash_nue_cc,  h_vtx_flash_nue_cc_mixed,
	                                          h_vtx_flash_nue_cc_out_fv,
	                                          h_vtx_flash_numu_cc, h_vtx_flash_numu_cc_mixed,
	                                          h_vtx_flash_cosmic,  h_vtx_flash_nc,
	                                          h_vtx_flash_nc_pi0,  h_vtx_flash_other_mixed,
	                                          h_vtx_flash_unmatched, h_vtx_flash_intime, intime_scale_factor,
	                                          h_vtx_flash_data, data_scale_factor, "",
	                                          "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_vtx_to_flash_distance_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_vtx_flash_nue_cc,  h_vtx_flash_nue_cc_mixed,
		                                          h_vtx_flash_nue_cc_out_fv,
		                                          h_vtx_flash_numu_cc, h_vtx_flash_numu_cc_mixed,
		                                          h_vtx_flash_cosmic,  h_vtx_flash_nc,
		                                          h_vtx_flash_nc_pi0,  h_vtx_flash_other_mixed,
		                                          h_vtx_flash_unmatched, h_vtx_flash_intime, scaled_intime_scale_factor,
		                                          h_vtx_flash_data, data_scale_factor, "",
		                                          "2D Distance From Largest Flash to Reco Nu Vtx [cm] (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_cuts_vtx_to_flash_distance_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStackData (h_vtx_flash_nue_cc,  h_vtx_flash_nue_cc_mixed,
	                                          h_vtx_flash_nue_cc_out_fv,
	                                          h_vtx_flash_numu_cc, h_vtx_flash_numu_cc_mixed,
	                                          h_vtx_flash_cosmic,  h_vtx_flash_nc,
	                                          h_vtx_flash_nc_pi0,  h_vtx_flash_other_mixed,
	                                          h_vtx_flash_unmatched, h_vtx_flash_intime, intime_scale_factor,
	                                          h_vtx_flash_data, data_scale_factor,
	                                          0.73, 0.98, 0.98, 0.50, true, "",
	                                          "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_vtx_to_flash_distance_data_logy.pdf"));

	histogram_functions::PlotSimpleStackData (h_vtx_flash_upstream_nue_cc,  h_vtx_flash_upstream_nue_cc_mixed,
	                                          h_vtx_flash_upstream_nue_cc_out_fv,
	                                          h_vtx_flash_upstream_numu_cc, h_vtx_flash_upstream_numu_cc_mixed,
	                                          h_vtx_flash_upstream_cosmic,  h_vtx_flash_upstream_nc,
	                                          h_vtx_flash_upstream_nc_pi0,  h_vtx_flash_upstream_other_mixed,
	                                          h_vtx_flash_upstream_unmatched, h_vtx_flash_upstream_intime, intime_scale_factor,
	                                          h_vtx_flash_upstream_data, data_scale_factor, "",
	                                          "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_vtx_to_flash_distance_upstream_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_vtx_flash_downstream_nue_cc,  h_vtx_flash_downstream_nue_cc_mixed,
	                                          h_vtx_flash_downstream_nue_cc_out_fv,
	                                          h_vtx_flash_downstream_numu_cc, h_vtx_flash_downstream_numu_cc_mixed,
	                                          h_vtx_flash_downstream_cosmic,  h_vtx_flash_downstream_nc,
	                                          h_vtx_flash_downstream_nc_pi0,  h_vtx_flash_downstream_other_mixed,
	                                          h_vtx_flash_downstream_unmatched, h_vtx_flash_downstream_intime, intime_scale_factor,
	                                          h_vtx_flash_downstream_data, data_scale_factor, "",
	                                          "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_vtx_to_flash_distance_downstream_data.pdf"));


	histogram_functions::PlotSimpleStackData (h_vtx_flash_nue_cc_after,  h_vtx_flash_nue_cc_mixed_after,
	                                          h_vtx_flash_nue_cc_out_fv_after,
	                                          h_vtx_flash_numu_cc_after, h_vtx_flash_numu_cc_mixed_after,
	                                          h_vtx_flash_cosmic_after,  h_vtx_flash_nc_after,
	                                          h_vtx_flash_nc_pi0_after,  h_vtx_flash_other_mixed_after,
	                                          h_vtx_flash_unmatched_after, h_vtx_flash_intime_after, intime_scale_factor,
	                                          h_vtx_flash_data_after, data_scale_factor, "",
	                                          "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_vtx_to_flash_distance_data_after.pdf"));

	histogram_functions::PlotSimpleStack (h_trk_vtx_dist_nue_cc,  h_trk_vtx_dist_nue_cc_mixed,
	                                      h_trk_vtx_dist_nue_cc_out_fv,
	                                      h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed,
	                                      h_trk_vtx_dist_cosmic,  h_trk_vtx_dist_nc,
	                                      h_trk_vtx_dist_nc_pi0,  h_trk_vtx_dist_other_mixed,
	                                      h_trk_vtx_dist_unmatched, "",
	                                      "Track to Nue Candidate Vertex Distance [cm]", "", Form("%s%s", file_locate_prefix, "post_cuts_track_to_vtx.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_trk_vtx_dist_nue_cc,  h_trk_vtx_dist_nue_cc_mixed,
	                                            h_trk_vtx_dist_nue_cc_out_fv,
	                                            h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed,
	                                            h_trk_vtx_dist_cosmic,  h_trk_vtx_dist_nc,
	                                            h_trk_vtx_dist_nc_pi0,  h_trk_vtx_dist_other_mixed,
	                                            h_trk_vtx_dist_unmatched, h_trk_vtx_dist_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Track to Nue Candidate Vertex Distance [cm]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_track_to_vtx_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_trk_vtx_dist_nue_cc,  h_trk_vtx_dist_nue_cc_mixed,
	                                          h_trk_vtx_dist_nue_cc_out_fv,
	                                          h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed,
	                                          h_trk_vtx_dist_cosmic,  h_trk_vtx_dist_nc,
	                                          h_trk_vtx_dist_nc_pi0,  h_trk_vtx_dist_other_mixed,
	                                          h_trk_vtx_dist_unmatched, h_trk_vtx_dist_intime, intime_scale_factor,
	                                          h_trk_vtx_dist_data, data_scale_factor, "",
	                                          "Track to Nue Candidate Vertex Distance [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_track_to_vtx_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_trk_vtx_dist_nue_cc,  h_trk_vtx_dist_nue_cc_mixed,
	                                          h_trk_vtx_dist_nue_cc_out_fv,
	                                          h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed,
	                                          h_trk_vtx_dist_cosmic,  h_trk_vtx_dist_nc,
	                                          h_trk_vtx_dist_nc_pi0,  h_trk_vtx_dist_other_mixed,
	                                          h_trk_vtx_dist_unmatched, h_trk_vtx_dist_intime, intime_scale_factor,
	                                          h_trk_vtx_dist_data, data_scale_factor,
	                                          0.73, 0.98, 0.98, 0.50, true, "",
	                                          "Track to Nue Candidate Vertex Distance [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_track_to_vtx_data_logy.pdf"));

	int tracks_num_signals = 0;
	int tracks_num_total = 0;
	for(int i = 0; i < h_trk_vtx_dist_nue_cc->GetNbinsX() + 1; i++)
	{
		tracks_num_signals += h_trk_vtx_dist_nue_cc->GetBinContent(i) * data_scale_factor;
		tracks_num_total += ((h_trk_vtx_dist_nue_cc->GetBinContent(i) + h_trk_vtx_dist_nue_cc_mixed->GetBinContent(i) +
		                      h_trk_vtx_dist_nue_cc_out_fv->GetBinContent(i) + h_trk_vtx_dist_numu_cc->GetBinContent(i) +
		                      h_trk_vtx_dist_numu_cc_mixed->GetBinContent(i) +
		                      h_trk_vtx_dist_cosmic->GetBinContent(i) + h_trk_vtx_dist_nc->GetBinContent(i) + h_trk_vtx_dist_nc_pi0->GetBinContent(i) +
		                      h_trk_vtx_dist_other_mixed->GetBinContent(i) + h_trk_vtx_dist_unmatched->GetBinContent(i)) * data_scale_factor +
		                     (h_trk_vtx_dist_intime->GetBinContent(i) * intime_scale_factor));
	}
	std::cout << "--- At this point in the selection we have: " << tracks_num_signals << " signal events with >= 1 Reco Track" << std::endl;
	std::cout << "--- And a total of : " << tracks_num_total <<  " events with >= 1 Reco Track" << std::endl;

	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_trk_vtx_dist_nue_cc,  h_trk_vtx_dist_nue_cc_mixed,
		                                          h_trk_vtx_dist_nue_cc_out_fv,
		                                          h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed,
		                                          h_trk_vtx_dist_cosmic,  h_trk_vtx_dist_nc,
		                                          h_trk_vtx_dist_nc_pi0,  h_trk_vtx_dist_other_mixed,
		                                          h_trk_vtx_dist_unmatched, h_trk_vtx_dist_intime, scaled_intime_scale_factor,
		                                          h_trk_vtx_dist_data, data_scale_factor,
		                                          0.73, 0.98, 0.98, 0.50, true, "",
		                                          "Track to Nue Candidate Vertex Distance [cm] (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_cuts_track_to_vtx_data_scaled_logy.pdf"));
	}

	histogram_functions::PlotSimpleStackData (h_trk_vtx_dist_nue_cc_after,  h_trk_vtx_dist_nue_cc_mixed_after,
	                                          h_trk_vtx_dist_nue_cc_out_fv_after,
	                                          h_trk_vtx_dist_numu_cc_after, h_trk_vtx_dist_numu_cc_mixed_after,
	                                          h_trk_vtx_dist_cosmic_after,  h_trk_vtx_dist_nc_after,
	                                          h_trk_vtx_dist_nc_pi0_after,  h_trk_vtx_dist_other_mixed_after,
	                                          h_trk_vtx_dist_unmatched_after, h_trk_vtx_dist_intime_after, intime_scale_factor,
	                                          h_trk_vtx_dist_data_after, data_scale_factor, "",
	                                          "Track to Nue Candidate Vertex Distance [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_track_to_vtx_data_after.pdf"));

	histogram_functions::PlotSimpleStack (h_shwr_vtx_dist_nue_cc,  h_shwr_vtx_dist_nue_cc_mixed,
	                                      h_shwr_vtx_dist_nue_cc_out_fv,
	                                      h_shwr_vtx_dist_numu_cc, h_shwr_vtx_dist_numu_cc_mixed,
	                                      h_shwr_vtx_dist_cosmic,  h_shwr_vtx_dist_nc,
	                                      h_shwr_vtx_dist_nc_pi0,  h_shwr_vtx_dist_other_mixed,
	                                      h_shwr_vtx_dist_unmatched, "",
	                                      "Shower to Nue Candidate Vertex Distance [cm]", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_shower_to_vtx.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_shwr_vtx_dist_nue_cc,  h_shwr_vtx_dist_nue_cc_mixed,
	                                            h_shwr_vtx_dist_nue_cc_out_fv,
	                                            h_shwr_vtx_dist_numu_cc, h_shwr_vtx_dist_numu_cc_mixed,
	                                            h_shwr_vtx_dist_cosmic,  h_shwr_vtx_dist_nc,
	                                            h_shwr_vtx_dist_nc_pi0,  h_shwr_vtx_dist_other_mixed,
	                                            h_shwr_vtx_dist_unmatched, h_shwr_vtx_dist_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Shower to Nue Candidate Vertex Distance [cm]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_shower_to_vtx_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_shwr_vtx_dist_nue_cc,  h_shwr_vtx_dist_nue_cc_mixed,
	                                          h_shwr_vtx_dist_nue_cc_out_fv,
	                                          h_shwr_vtx_dist_numu_cc, h_shwr_vtx_dist_numu_cc_mixed,
	                                          h_shwr_vtx_dist_cosmic,  h_shwr_vtx_dist_nc,
	                                          h_shwr_vtx_dist_nc_pi0,  h_shwr_vtx_dist_other_mixed,
	                                          h_shwr_vtx_dist_unmatched, h_shwr_vtx_dist_intime, intime_scale_factor,
	                                          h_shwr_vtx_dist_data, data_scale_factor, "",
	                                          "Shower to Nue Candidate Vertex Distance [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_shower_to_vtx_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_shwr_vtx_dist_nue_cc,  h_shwr_vtx_dist_nue_cc_mixed,
	                                          h_shwr_vtx_dist_nue_cc_out_fv,
	                                          h_shwr_vtx_dist_numu_cc, h_shwr_vtx_dist_numu_cc_mixed,
	                                          h_shwr_vtx_dist_cosmic,  h_shwr_vtx_dist_nc,
	                                          h_shwr_vtx_dist_nc_pi0,  h_shwr_vtx_dist_other_mixed,
	                                          h_shwr_vtx_dist_unmatched, h_shwr_vtx_dist_intime, intime_scale_factor,
	                                          h_shwr_vtx_dist_data, data_scale_factor,
	                                          0.73, 0.98, 0.98, 0.50, true, "",
	                                          "Shower to Nue Candidate Vertex Distance [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_shower_to_vtx_data_logy.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_shwr_vtx_dist_nue_cc,  h_shwr_vtx_dist_nue_cc_mixed,
		                                          h_shwr_vtx_dist_nue_cc_out_fv,
		                                          h_shwr_vtx_dist_numu_cc, h_shwr_vtx_dist_numu_cc_mixed,
		                                          h_shwr_vtx_dist_cosmic,  h_shwr_vtx_dist_nc,
		                                          h_shwr_vtx_dist_nc_pi0,  h_shwr_vtx_dist_other_mixed,
		                                          h_shwr_vtx_dist_unmatched, h_shwr_vtx_dist_intime, scaled_intime_scale_factor,
		                                          h_shwr_vtx_dist_data, data_scale_factor,
		                                          0.73, 0.98, 0.98, 0.50, true, "",
		                                          "Shower to Nue Candidate Vertex Distance [cm] (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_cuts_shower_to_vtx_data_scaled_logy.pdf"));
	}

	int total_total = 0;
	double bin_val_total = 0;
	for(int i = 0; i < h_shwr_vtx_dist_nue_cc->GetNbinsX()+2; i++)
	{
		const double bin_val = h_shwr_vtx_dist_nue_cc->GetBinContent(i);
		std::cout << "Distance: " << h_shwr_vtx_dist_nue_cc->GetBinLowEdge(i) << ", Num Signal: " <<  bin_val << std::endl;
		if(h_shwr_vtx_dist_nue_cc->GetBinLowEdge(i) > 4)
		{
			bin_val_total += bin_val;
		}
		total_total += bin_val;
	}
	std::cout << "Num Signals: " << total_total << std::endl;
	std::cout << "Removed: " << bin_val_total << std::endl;

	histogram_functions::PlotSimpleStackData (h_shwr_vtx_dist_nue_cc_after,  h_shwr_vtx_dist_nue_cc_mixed_after,
	                                          h_shwr_vtx_dist_nue_cc_out_fv_after,
	                                          h_shwr_vtx_dist_numu_cc_after, h_shwr_vtx_dist_numu_cc_mixed_after,
	                                          h_shwr_vtx_dist_cosmic_after,  h_shwr_vtx_dist_nc_after,
	                                          h_shwr_vtx_dist_nc_pi0_after,  h_shwr_vtx_dist_other_mixed_after,
	                                          h_shwr_vtx_dist_unmatched_after, h_shwr_vtx_dist_intime_after, intime_scale_factor,
	                                          h_shwr_vtx_dist_data_after, data_scale_factor, "",
	                                          "Shower to Nue Candidate Vertex Distance [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_shower_to_vtx_data_after.pdf"));

	histogram_functions::PlotDetailStack(h_pfp_track_nue_cc_qe,
	                                     h_pfp_track_nue_cc_out_fv,
	                                     h_pfp_track_nue_cc_res,
	                                     h_pfp_track_nue_cc_dis,
	                                     h_pfp_track_nue_cc_coh,
	                                     h_pfp_track_nue_cc_mec,
	                                     h_pfp_track_nue_cc_mixed,
	                                     h_pfp_track_numu_cc_qe,
	                                     h_pfp_track_numu_cc_res,
	                                     h_pfp_track_numu_cc_dis,
	                                     h_pfp_track_numu_cc_coh,
	                                     h_pfp_track_numu_cc_mec,
	                                     h_pfp_track_numu_cc_mixed,
	                                     h_pfp_track_cosmic,
	                                     h_pfp_track_nc,
	                                     h_pfp_track_nc_pi0,
	                                     h_pfp_track_other_mixed,
	                                     h_pfp_track_unmatched,
	                                     "Reconstructed Tracks in Candidate Neutrino Object", "",
	                                     Form("%s%s", file_locate_prefix, "selected_pfp_track_stack.pdf"));

	histogram_functions::PlotDetailStack(h_pfp_shower_nue_cc_qe,
	                                     h_pfp_shower_nue_cc_out_fv,
	                                     h_pfp_shower_nue_cc_res,
	                                     h_pfp_shower_nue_cc_dis,
	                                     h_pfp_shower_nue_cc_coh,
	                                     h_pfp_shower_nue_cc_mec,
	                                     h_pfp_shower_nue_cc_mixed,
	                                     h_pfp_shower_numu_cc_qe,
	                                     h_pfp_shower_numu_cc_res,
	                                     h_pfp_shower_numu_cc_dis,
	                                     h_pfp_shower_numu_cc_coh,
	                                     h_pfp_shower_numu_cc_mec,
	                                     h_pfp_shower_numu_cc_mixed,
	                                     h_pfp_shower_cosmic,
	                                     h_pfp_shower_nc,
	                                     h_pfp_shower_nc_pi0,
	                                     h_pfp_shower_other_mixed,
	                                     h_pfp_shower_unmatched,
	                                     "Reconstructed Showers in Candidate Neutrino Object", "",
	                                     Form("%s%s", file_locate_prefix, "selected_pfp_shower_stack.pdf"));

	histogram_functions::PlotDetailStack(h_pfp_shower_open_angle_nue_cc_qe,
	                                     h_pfp_shower_open_angle_nue_cc_out_fv,
	                                     h_pfp_shower_open_angle_nue_cc_res,
	                                     h_pfp_shower_open_angle_nue_cc_dis,
	                                     h_pfp_shower_open_angle_nue_cc_coh,
	                                     h_pfp_shower_open_angle_nue_cc_mec,
	                                     h_pfp_shower_open_angle_nue_cc_mixed,
	                                     h_pfp_shower_open_angle_numu_cc_qe,
	                                     h_pfp_shower_open_angle_numu_cc_res,
	                                     h_pfp_shower_open_angle_numu_cc_dis,
	                                     h_pfp_shower_open_angle_numu_cc_coh,
	                                     h_pfp_shower_open_angle_numu_cc_mec,
	                                     h_pfp_shower_open_angle_numu_cc_mixed,
	                                     h_pfp_shower_open_angle_cosmic,
	                                     h_pfp_shower_open_angle_nc,
	                                     h_pfp_shower_open_angle_nc_pi0,
	                                     h_pfp_shower_open_angle_other_mixed,
	                                     h_pfp_shower_open_angle_unmatched,
	                                     "Reconstructed Showers in Candidate Neutrino Object", "",
	                                     Form("%s%s", file_locate_prefix, "selected_pfp_shower_open_angle_stack.pdf"));

	histogram_functions::PlotDetailStack(h_pfp_shower_dedx_nue_cc_qe,
	                                     h_pfp_shower_dedx_nue_cc_out_fv,
	                                     h_pfp_shower_dedx_nue_cc_res,
	                                     h_pfp_shower_dedx_nue_cc_dis,
	                                     h_pfp_shower_dedx_nue_cc_coh,
	                                     h_pfp_shower_dedx_nue_cc_mec,
	                                     h_pfp_shower_dedx_nue_cc_mixed,
	                                     h_pfp_shower_dedx_numu_cc_qe,
	                                     h_pfp_shower_dedx_numu_cc_res,
	                                     h_pfp_shower_dedx_numu_cc_dis,
	                                     h_pfp_shower_dedx_numu_cc_coh,
	                                     h_pfp_shower_dedx_numu_cc_mec,
	                                     h_pfp_shower_dedx_numu_cc_mixed,
	                                     h_pfp_shower_dedx_cosmic,
	                                     h_pfp_shower_dedx_nc,
	                                     h_pfp_shower_dedx_nc_pi0,
	                                     h_pfp_shower_dedx_other_mixed,
	                                     h_pfp_shower_dedx_unmatched,
	                                     "Reconstructed Showers in Candidate Neutrino Object", "",
	                                     Form("%s%s", file_locate_prefix, "selected_pfp_shower_dedx_stack.pdf"));

	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_qe, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_qe.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_out_fv, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_res, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_res.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_dis, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_dis.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_coh, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_coh.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mec, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_mec.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mixed, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_qe, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_qe.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_res, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_res.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_dis, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_dis.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_coh, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_coh.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mec, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_mec.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mixed, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nc.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc_pi0, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nc_pi0.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_cosmic, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_cosmic.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_other_mixed, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_other_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_unmatched, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_unmatched.pdf"));


	histogram_functions::PlotDetailStack(h_pfp_track_nue_cc_qe_last,
	                                     h_pfp_track_nue_cc_out_fv_last,
	                                     h_pfp_track_nue_cc_res_last,
	                                     h_pfp_track_nue_cc_dis_last,
	                                     h_pfp_track_nue_cc_coh_last,
	                                     h_pfp_track_nue_cc_mec_last,
	                                     h_pfp_track_nue_cc_mixed_last,
	                                     h_pfp_track_numu_cc_qe_last,
	                                     h_pfp_track_numu_cc_res_last,
	                                     h_pfp_track_numu_cc_dis_last,
	                                     h_pfp_track_numu_cc_coh_last,
	                                     h_pfp_track_numu_cc_mec_last,
	                                     h_pfp_track_numu_cc_mixed_last,
	                                     h_pfp_track_cosmic_last,
	                                     h_pfp_track_nc_last,
	                                     h_pfp_track_nc_pi0_last,
	                                     h_pfp_track_other_mixed_last,
	                                     h_pfp_track_unmatched_last,
	                                     "Reconstructed Tracks in Candidate Neutrino Object", "",
	                                     Form("%s%s", file_locate_prefix, "selected_pfp_track_stack_last.pdf"));

	histogram_functions::PlotDetailStack(h_pfp_shower_nue_cc_qe_last,
	                                     h_pfp_shower_nue_cc_out_fv_last,
	                                     h_pfp_shower_nue_cc_res_last,
	                                     h_pfp_shower_nue_cc_dis_last,
	                                     h_pfp_shower_nue_cc_coh_last,
	                                     h_pfp_shower_nue_cc_mec_last,
	                                     h_pfp_shower_nue_cc_mixed_last,
	                                     h_pfp_shower_numu_cc_qe_last,
	                                     h_pfp_shower_numu_cc_res_last,
	                                     h_pfp_shower_numu_cc_dis_last,
	                                     h_pfp_shower_numu_cc_coh_last,
	                                     h_pfp_shower_numu_cc_mec_last,
	                                     h_pfp_shower_numu_cc_mixed_last,
	                                     h_pfp_shower_cosmic_last,
	                                     h_pfp_shower_nc_last,
	                                     h_pfp_shower_nc_pi0_last,
	                                     h_pfp_shower_other_mixed_last,
	                                     h_pfp_shower_unmatched_last,
	                                     "Reconstructed Showers in Candidate Neutrino Object", "",
	                                     Form("%s%s", file_locate_prefix, "selected_pfp_shower_stack_last.pdf"));

	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_qe_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_qe_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_out_fv_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_out_fv_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_res_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_res_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_dis_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_dis_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_coh_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_coh_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mec_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_mec_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mixed_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nue_cc_mixed_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_qe_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_qe_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_res_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_res_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_dis_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_dis_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_coh_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_coh_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mec_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_mec_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mixed_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_numu_cc_mixed_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nc_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc_pi0_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_nc_pi0_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_cosmic_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_cosmic_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_other_mixed_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_other_mixed_last.pdf"));
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_unmatched_last, "", "PFP Tracks", "PFP Showers",
	                                      Form("%s%s", file_locate_prefix, "selected_pfp_track_shower_unmatched_last.pdf"));

	const char * str_origin[3] = {"kBeamNeutrino", "kCosmicRay", "kUnknown"};
	for (int i=1; i<= 3; i++)
	{
		h_leading_shower_mc_pdg_nue_cc_qe->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_nue_cc_out_fv->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_nue_cc_res->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_nue_cc_coh->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_nue_cc_dis->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_nue_cc_mec->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_numu_cc_qe->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_numu_cc_res->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_numu_cc_dis->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_numu_cc_coh->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_numu_cc_mec->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_nc->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_nc_pi0->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_nue_cc_mixed->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_numu_cc_mixed->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_other_mixed->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_unmatched->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_cosmic->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
	}
	const char * str_mc_particle[10] = {"Electron", "Positron", "Muon", "Mu+", "Photon",
		                            "Pion", "Proton", "Neutron", "Kaon", "No Match"};
	for (int i = 1; i <= 10; i++)
	{
		h_leading_shower_mc_pdg_nue_cc_qe->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_nue_cc_out_fv->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_nue_cc_res->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_nue_cc_coh->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_nue_cc_dis->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_nue_cc_mec->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_numu_cc_qe->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_numu_cc_res->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_numu_cc_dis->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_numu_cc_coh->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_numu_cc_mec->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_nc->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_nc_pi0->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_nue_cc_mixed->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_numu_cc_mixed->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_other_mixed->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_unmatched->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_cosmic->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
	}

	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_qe, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_nue_cc_qe.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_out_fv, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_res, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_nue_cc_res.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_dis, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_nue_cc_dis.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_coh, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_nue_cc_coh.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_mec, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_nue_cc_mec.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_mixed, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_nue_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_qe, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_numu_cc_qe.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_res, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_numu_cc_res.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_dis, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_numu_cc_dis.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_coh, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_numu_cc_coh.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_mec, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_numu_cc_mec.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_mixed, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_numu_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nc, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_nc.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nc_pi0, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_nc_pi0.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_cosmic, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_cosmic.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_other_mixed, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdf_other_mixed.pdf"));
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_unmatched, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle",
	                                            Form("%s%s", file_locate_prefix, "selected_leading_shower_mc_pdg_unmatched.pdf"));

	histogram_functions::Plot2DHistogram (h_shwr_hits_nu_eng, "", "True Neutrino Energy [GeV]", "Signal Electron Shower Hits",
	                                      Form("%s%s", file_locate_prefix, "shwr_hits_nu_eng.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_hits_ele_eng, "", "True Electron Energy [GeV]", "Signal Electron Shower Hits",
	                                      Form("%s%s", file_locate_prefix, "shwr_hits_ele_eng.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_hits_nu_eng_zoom, "", "True Neutrino Energy [GeV]", "Signal Electron Shower Hits",
	                                      Form("%s%s", file_locate_prefix, "shwr_hits_nu_eng_zoom.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_hits_ele_eng_zoom, "", "True Electron Energy [GeV]", "Signal Electron Shower Hits",
	                                      Form("%s%s", file_locate_prefix, "shwr_hits_ele_eng_zoom.pdf"));

	histogram_functions::Plot2DHistogram (h_shwr_hits_nu_eng_last, "", "True Neutrino Energy [GeV]", "Signal Electron Shower Hits",
	                                      Form("%s%s", file_locate_prefix, "shwr_hits_nu_eng_last.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_hits_ele_eng_last, "", "True Electron Energy [GeV]", "Signal Electron Shower Hits",
	                                      Form("%s%s", file_locate_prefix, "shwr_hits_ele_eng_last.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_hits_nu_eng_zoom_last, "", "True Neutrino Energy [GeV]", "Signal Electron Shower Hits",
	                                      Form("%s%s", file_locate_prefix, "shwr_hits_nu_eng_zoom_last.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_hits_ele_eng_zoom_last, "", "True Electron Energy [GeV]", "Signal Electron Shower Hits",
	                                      Form("%s%s", file_locate_prefix, "shwr_hits_ele_eng_zoom_last.pdf"));

	histogram_functions::Plot2DHistogram (h_shwr_collection_hits_nu_eng, "", "True Neutrino Energy [GeV]",
	                                      "Signal Electron Shower Hits (Collection)",
	                                      Form("%s%s", file_locate_prefix, "shwr_collection_hits_nu_eng.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_collection_hits_ele_eng, "", "True Electron Energy [GeV]",
	                                      "Signal Electron Shower Hits (Collection)",
	                                      Form("%s%s", file_locate_prefix, "shwr_collection_hits_ele_eng.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_collection_hits_nu_eng_zoom, "", "True Neutrino Energy [GeV]",
	                                      "Signal Electron Shower Hits (Collection)",
	                                      Form("%s%s", file_locate_prefix, "shwr_collection_hits_nu_eng_zoom.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_collection_hits_ele_eng_zoom, "", "True Electron Energy [GeV]",
	                                      "Signal Electron Shower Hits (Collection)",
	                                      Form("%s%s", file_locate_prefix, "shwr_collection_hits_ele_eng_zoom.pdf"));

	histogram_functions::Plot2DHistogram (h_shwr_collection_hits_nu_eng_last, "", "True Neutrino Energy [GeV]",
	                                      "Signal Electron Shower Hits (Collection)",
	                                      Form("%s%s", file_locate_prefix, "shwr_collection_hits_nu_eng_last.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_collection_hits_ele_eng_last, "", "True Electron Energy [GeV]",
	                                      "Signal Electron Shower Hits (Collection)",
	                                      Form("%s%s", file_locate_prefix, "shwr_collection_hits_ele_eng_last.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_collection_hits_nu_eng_zoom_last, "", "True Neutrino Energy [GeV]",
	                                      "Signal Electron Shower Hits (Collection)",
	                                      Form("%s%s", file_locate_prefix, "shwr_collection_hits_nu_eng_zoom_last.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_collection_hits_ele_eng_zoom_last, "", "True Electron Energy [GeV]",
	                                      "Signal Electron Shower Hits (Collection)",
	                                      Form("%s%s", file_locate_prefix, "shwr_collection_hits_ele_eng_zoom_last.pdf"));

	h_pfp_shower_nue_cc_qe_last->SetFillColor(30);
	h_pfp_shower_nue_cc_out_fv_last->SetFillColor(45);
	h_pfp_shower_nue_cc_res_last->SetFillColor(31);
	h_pfp_shower_nue_cc_dis_last->SetFillColor(32);
	h_pfp_shower_nue_cc_coh_last->SetFillColor(33);
	h_pfp_shower_nue_cc_mec_last->SetFillColor(34);
	h_pfp_shower_nc_last->SetFillColor(46);
	h_pfp_shower_numu_cc_qe_last->SetFillColor(28);
	h_pfp_shower_numu_cc_res_last->SetFillColor(27);
	h_pfp_shower_numu_cc_dis_last->SetFillColor(26);
	h_pfp_shower_numu_cc_coh_last->SetFillColor(23);
	h_pfp_shower_numu_cc_mec_last->SetFillColor(22);
	h_pfp_shower_nc_pi0_last->SetFillColor(36);
	h_pfp_shower_nue_cc_mixed_last->SetFillColor(38);
	h_pfp_shower_numu_cc_mixed_last->SetFillColor(25);
	h_pfp_shower_cosmic_last->SetFillColor(39);
	h_pfp_shower_other_mixed_last->SetFillColor(42);
	h_pfp_shower_unmatched_last->SetFillColor(12);


	histogram_functions::PostHistogramOverlay(h_selected_nu_energy_no_cut,        h_selected_nu_energy_reco_nue,
	                                          h_selected_nu_energy_in_fv,         h_selected_nu_energy_vtx_flash,
	                                          h_selected_nu_energy_shwr_vtx,      h_selected_nu_energy_trk_vtx,
	                                          h_selected_nu_energy_hit_threshold, h_selected_nu_energy_open_angle,
	                                          h_selected_nu_energy_dedx,
	                                          "True Signal Neutrino Energy [GeV]", "", Form("%s%s", file_locate_prefix, "sequential_nu_energy.pdf"));
	histogram_functions::PostHistogramOverlay(h_selected_ele_energy_no_cut,       h_selected_ele_energy_reco_nue,
	                                          h_selected_ele_energy_in_fv,        h_selected_ele_energy_vtx_flash,
	                                          h_selected_ele_energy_shwr_vtx,     h_selected_ele_energy_trk_vtx,
	                                          h_selected_ele_energy_hit_threshold,h_selected_ele_energy_open_angle,
	                                          h_selected_ele_energy_dedx,
	                                          "True Signal Electron Energy [GeV]", "", Form("%s%s", file_locate_prefix, "sequential_ele_energy.pdf"));


	histogram_functions::Plot1DHistogram (h_charge_share_nue_cc_mixed, "Neutrino Charge Fraction - Selected Nue CC Mixed",
	                                      Form("%s%s", file_locate_prefix, "charge_fraction_nue_cc_mixed.pdf"));
	histogram_functions::Plot1DHistogram (h_flash_t0_diff, "Largest Flash Time - True Neutrino Interaction Time [#mus]",
	                                      Form("%s%s", file_locate_prefix, "flash_t0_diff.pdf"));

	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nue_cc, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nue_cc_out_fv, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nue_cc_mixed, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_nue_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nue_cc, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_numu_cc, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_numu_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_numu_cc_mixed, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_numu_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nc, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_nc.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nc_pi0, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_nc_pi0.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_cosmic, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_cosmic.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_other_mixed, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_other_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_unmatched, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", Form("%s%s", file_locate_prefix, "dedx_open_angle_unmatched.pdf"));

	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nue_cc, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nue_cc_out_fv, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nue_cc_mixed, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_nue_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_numu_cc, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_numu_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_numu_cc_mixed, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_numu_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nc, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_nc.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nc_pi0, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_nc_pi0.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_cosmic, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_cosmic.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_other_mixed, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_other_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_unmatched, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", Form("%s%s", file_locate_prefix, "shwr_len_hits_unmatched.pdf"));

	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nue_cc, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_nue_cc.pdf"));
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nue_cc_out_fv, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_nue_cc_out_fv.pdf"));
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nue_cc_mixed, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_nue_cc_mixed.pdf"));
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_numu_cc, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_numu_cc.pdf"));
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_numu_cc_mixed, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_numu_cc_mixed.pdf"));
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nc, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_nc.pdf"));
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nc_pi0, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_nc_pi0.pdf"));
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_cosmic, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_cosmic.pdf"));
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_other_mixed, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_other_mixed.pdf"));
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_unmatched, "(TPCO w/ > 1 showers) Shower-Vtx Distance [cm]",
	                                      Form("%s%s", file_locate_prefix, "second_shwr_dist_unmatched.pdf"));

	histogram_functions::PlotSimpleStack (h_second_shwr_dist_nue_cc,  h_second_shwr_dist_nue_cc_mixed,
	                                      h_second_shwr_dist_nue_cc_out_fv,
	                                      h_second_shwr_dist_numu_cc, h_second_shwr_dist_numu_cc_mixed,
	                                      h_second_shwr_dist_cosmic,  h_second_shwr_dist_nc,
	                                      h_second_shwr_dist_nc_pi0,  h_second_shwr_dist_other_mixed,
	                                      h_second_shwr_dist_unmatched, "",
	                                      "(TPCO > 1 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "",
	                                      Form("%s%s", file_locate_prefix, "post_second_shwr_dist.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_second_shwr_dist_nue_cc,  h_second_shwr_dist_nue_cc_mixed,
	                                            h_second_shwr_dist_nue_cc_out_fv,
	                                            h_second_shwr_dist_numu_cc, h_second_shwr_dist_numu_cc_mixed,
	                                            h_second_shwr_dist_cosmic,  h_second_shwr_dist_nc,
	                                            h_second_shwr_dist_nc_pi0,  h_second_shwr_dist_other_mixed,
	                                            h_second_shwr_dist_unmatched, h_second_shwr_dist_intime, intime_scale_factor, data_scale_factor, "",
	                                            "(TPCO > 1 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "",
	                                            Form("%s%s", file_locate_prefix, "post_second_shwr_dist_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_second_shwr_dist_nue_cc,  h_second_shwr_dist_nue_cc_mixed,
	                                          h_second_shwr_dist_nue_cc_out_fv,
	                                          h_second_shwr_dist_numu_cc, h_second_shwr_dist_numu_cc_mixed,
	                                          h_second_shwr_dist_cosmic,  h_second_shwr_dist_nc,
	                                          h_second_shwr_dist_nc_pi0,  h_second_shwr_dist_other_mixed,
	                                          h_second_shwr_dist_unmatched, h_second_shwr_dist_intime, intime_scale_factor,
	                                          h_second_shwr_dist_data, data_scale_factor, "",
	                                          "(TPCO > 1 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_second_shwr_dist_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_second_shwr_dist_nue_cc,  h_second_shwr_dist_nue_cc_mixed,
		                                          h_second_shwr_dist_nue_cc_out_fv,
		                                          h_second_shwr_dist_numu_cc, h_second_shwr_dist_numu_cc_mixed,
		                                          h_second_shwr_dist_cosmic,  h_second_shwr_dist_nc,
		                                          h_second_shwr_dist_nc_pi0,  h_second_shwr_dist_other_mixed,
		                                          h_second_shwr_dist_unmatched, h_second_shwr_dist_intime, scaled_intime_scale_factor,
		                                          h_second_shwr_dist_data, data_scale_factor, "",
		                                          "(TPCO > 1 Reco Shower) Secondary Shwr-Vtx Distance [cm] (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_second_shwr_dist_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStackData (h_second_shwr_dist_nue_cc,  h_second_shwr_dist_nue_cc_mixed,
	                                          h_second_shwr_dist_nue_cc_out_fv,
	                                          h_second_shwr_dist_numu_cc, h_second_shwr_dist_numu_cc_mixed,
	                                          h_second_shwr_dist_cosmic,  h_second_shwr_dist_nc,
	                                          h_second_shwr_dist_nc_pi0,  h_second_shwr_dist_other_mixed,
	                                          h_second_shwr_dist_unmatched, h_second_shwr_dist_intime, intime_scale_factor,
	                                          h_second_shwr_dist_data, data_scale_factor,
	                                          0.73, 0.98, 0.98, 0.50, true, "",
	                                          "(TPCO > 1 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_second_shwr_dist_data_logy.pdf"));

	histogram_functions::PlotSimpleStackData (h_second_shwr_dist_nue_cc_after,  h_second_shwr_dist_nue_cc_mixed_after,
	                                          h_second_shwr_dist_nue_cc_out_fv_after,
	                                          h_second_shwr_dist_numu_cc_after, h_second_shwr_dist_numu_cc_mixed_after,
	                                          h_second_shwr_dist_cosmic_after,  h_second_shwr_dist_nc_after,
	                                          h_second_shwr_dist_nc_pi0_after,  h_second_shwr_dist_other_mixed_after,
	                                          h_second_shwr_dist_unmatched_after, h_second_shwr_dist_intime_after, intime_scale_factor,
	                                          h_second_shwr_dist_data_after, data_scale_factor, "",
	                                          "(TPCO > 1 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_second_shwr_dist_data_after.pdf"));

	histogram_functions::PlotSimpleStack (h_hit_length_ratio_nue_cc,  h_hit_length_ratio_nue_cc_mixed,
	                                      h_hit_length_ratio_nue_cc_out_fv,
	                                      h_hit_length_ratio_numu_cc, h_hit_length_ratio_numu_cc_mixed,
	                                      h_hit_length_ratio_cosmic,  h_hit_length_ratio_nc,
	                                      h_hit_length_ratio_nc_pi0,  h_hit_length_ratio_other_mixed,
	                                      h_hit_length_ratio_unmatched, "",
	                                      "Leading Shower (Hits / Length) [cm^-1]", "", Form("%s%s", file_locate_prefix, "post_hit_length_ratio.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_hit_length_ratio_nue_cc,  h_hit_length_ratio_nue_cc_mixed,
	                                            h_hit_length_ratio_nue_cc_out_fv,
	                                            h_hit_length_ratio_numu_cc, h_hit_length_ratio_numu_cc_mixed,
	                                            h_hit_length_ratio_cosmic,  h_hit_length_ratio_nc,
	                                            h_hit_length_ratio_nc_pi0,  h_hit_length_ratio_other_mixed,
	                                            h_hit_length_ratio_unmatched, h_hit_length_ratio_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower (Hits / Length) [cm^-1]", "",
	                                            Form("%s%s", file_locate_prefix, "post_hit_length_ratio_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_hit_length_ratio_nue_cc,  h_hit_length_ratio_nue_cc_mixed,
	                                          h_hit_length_ratio_nue_cc_out_fv,
	                                          h_hit_length_ratio_numu_cc, h_hit_length_ratio_numu_cc_mixed,
	                                          h_hit_length_ratio_cosmic,  h_hit_length_ratio_nc,
	                                          h_hit_length_ratio_nc_pi0,  h_hit_length_ratio_other_mixed,
	                                          h_hit_length_ratio_unmatched, h_hit_length_ratio_intime, intime_scale_factor,
	                                          h_hit_length_ratio_data, data_scale_factor, "",
	                                          "Leading Shower (Hits / Length) [cm^-1]", "",
	                                          Form("%s%s", file_locate_prefix, "post_hit_length_ratio_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_hit_length_ratio_nue_cc,  h_hit_length_ratio_nue_cc_mixed,
		                                          h_hit_length_ratio_nue_cc_out_fv,
		                                          h_hit_length_ratio_numu_cc, h_hit_length_ratio_numu_cc_mixed,
		                                          h_hit_length_ratio_cosmic,  h_hit_length_ratio_nc,
		                                          h_hit_length_ratio_nc_pi0,  h_hit_length_ratio_other_mixed,
		                                          h_hit_length_ratio_unmatched, h_hit_length_ratio_intime, scaled_intime_scale_factor,
		                                          h_hit_length_ratio_data, data_scale_factor, "",
		                                          "Leading Shower (Hits / Length) [cm^-1] (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_hit_length_ratio_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStackData (h_hit_length_ratio_nue_cc_after,  h_hit_length_ratio_nue_cc_mixed_after,
	                                          h_hit_length_ratio_nue_cc_out_fv_after,
	                                          h_hit_length_ratio_numu_cc_after, h_hit_length_ratio_numu_cc_mixed_after,
	                                          h_hit_length_ratio_cosmic_after,  h_hit_length_ratio_nc_after,
	                                          h_hit_length_ratio_nc_pi0_after,  h_hit_length_ratio_other_mixed_after,
	                                          h_hit_length_ratio_unmatched_after, h_hit_length_ratio_intime_after, intime_scale_factor,
	                                          h_hit_length_ratio_data_after, data_scale_factor, "",
	                                          "Leading Shower (Hits / Length) [cm^-1]", "",
	                                          Form("%s%s", file_locate_prefix, "post_hit_length_ratio_data_after.pdf"));

	histogram_functions::PlotSimpleStack (h_trk_length_nue_cc,  h_trk_length_nue_cc_mixed,
	                                      h_trk_length_nue_cc_out_fv,
	                                      h_trk_length_numu_cc, h_trk_length_numu_cc_mixed,
	                                      h_trk_length_cosmic,  h_trk_length_nc,
	                                      h_trk_length_nc_pi0,  h_trk_length_other_mixed,
	                                      h_trk_length_unmatched, "",
	                                      "All Track Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_all_track_lengths.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_trk_length_nue_cc,  h_trk_length_nue_cc_mixed,
	                                            h_trk_length_nue_cc_out_fv,
	                                            h_trk_length_numu_cc, h_trk_length_numu_cc_mixed,
	                                            h_trk_length_cosmic,  h_trk_length_nc,
	                                            h_trk_length_nc_pi0,  h_trk_length_other_mixed,
	                                            h_trk_length_unmatched, h_trk_length_intime, intime_scale_factor, data_scale_factor, "",
	                                            "All Track Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_all_track_lengths_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_trk_length_nue_cc,  h_trk_length_nue_cc_mixed,
	                                          h_trk_length_nue_cc_out_fv,
	                                          h_trk_length_numu_cc, h_trk_length_numu_cc_mixed,
	                                          h_trk_length_cosmic,  h_trk_length_nc,
	                                          h_trk_length_nc_pi0,  h_trk_length_other_mixed,
	                                          h_trk_length_unmatched, h_trk_length_intime, intime_scale_factor,
	                                          h_trk_length_data, data_scale_factor, "",
	                                          "All Track Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_all_track_lengths_data.pdf"));

	histogram_functions::PlotSimpleStack (h_longest_trk_length_nue_cc,  h_longest_trk_length_nue_cc_mixed,
	                                      h_longest_trk_length_nue_cc_out_fv,
	                                      h_longest_trk_length_numu_cc, h_longest_trk_length_numu_cc_mixed,
	                                      h_longest_trk_length_cosmic,  h_longest_trk_length_nc,
	                                      h_longest_trk_length_nc_pi0,  h_longest_trk_length_other_mixed,
	                                      h_longest_trk_length_unmatched, "",
	                                      "Longest Track Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_longest_track_lengths.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_longest_trk_length_nue_cc,  h_longest_trk_length_nue_cc_mixed,
	                                            h_longest_trk_length_nue_cc_out_fv,
	                                            h_longest_trk_length_numu_cc, h_longest_trk_length_numu_cc_mixed,
	                                            h_longest_trk_length_cosmic,  h_longest_trk_length_nc,
	                                            h_longest_trk_length_nc_pi0,  h_longest_trk_length_other_mixed,
	                                            h_longest_trk_length_unmatched, h_longest_trk_length_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Longest Track Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_longest_track_lengths_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_longest_trk_length_nue_cc,  h_longest_trk_length_nue_cc_mixed,
	                                          h_longest_trk_length_nue_cc_out_fv,
	                                          h_longest_trk_length_numu_cc, h_longest_trk_length_numu_cc_mixed,
	                                          h_longest_trk_length_cosmic,  h_longest_trk_length_nc,
	                                          h_longest_trk_length_nc_pi0,  h_longest_trk_length_other_mixed,
	                                          h_longest_trk_length_unmatched, h_longest_trk_length_intime, intime_scale_factor,
	                                          h_longest_trk_length_data, data_scale_factor, "",
	                                          "Longest Track Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_longest_track_lengths_data.pdf"));

	histogram_functions::PlotSimpleStack (h_shwr_length_nue_cc,  h_shwr_length_nue_cc_mixed,
	                                      h_shwr_length_nue_cc_out_fv,
	                                      h_shwr_length_numu_cc, h_shwr_length_numu_cc_mixed,
	                                      h_shwr_length_cosmic,  h_shwr_length_nc,
	                                      h_shwr_length_nc_pi0,  h_shwr_length_other_mixed,
	                                      h_shwr_length_unmatched, "",
	                                      "All Shower Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_all_shower_lengths.pdf"));
	histogram_functions::PlotSimpleStack (h_longest_shwr_length_nue_cc,  h_longest_shwr_length_nue_cc_mixed,
	                                      h_longest_shwr_length_nue_cc_out_fv,
	                                      h_longest_shwr_length_numu_cc, h_longest_shwr_length_numu_cc_mixed,
	                                      h_longest_shwr_length_cosmic,  h_longest_shwr_length_nc,
	                                      h_longest_shwr_length_nc_pi0,  h_longest_shwr_length_other_mixed,
	                                      h_longest_shwr_length_unmatched, "",
	                                      "Longest Shower Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_longest_shower_lengths.pdf"));
	histogram_functions::PlotSimpleStack (h_leading_shwr_length_nue_cc,  h_leading_shwr_length_nue_cc_mixed,
	                                      h_leading_shwr_length_nue_cc_out_fv,
	                                      h_leading_shwr_length_numu_cc, h_leading_shwr_length_numu_cc_mixed,
	                                      h_leading_shwr_length_cosmic,  h_leading_shwr_length_nc,
	                                      h_leading_shwr_length_nc_pi0,  h_leading_shwr_length_other_mixed,
	                                      h_leading_shwr_length_unmatched, "",
	                                      "Leading Shower Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_leading_shower_lengths.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_shwr_length_nue_cc,  h_shwr_length_nue_cc_mixed,
	                                            h_shwr_length_nue_cc_out_fv,
	                                            h_shwr_length_numu_cc, h_shwr_length_numu_cc_mixed,
	                                            h_shwr_length_cosmic,  h_shwr_length_nc,
	                                            h_shwr_length_nc_pi0,  h_shwr_length_other_mixed,
	                                            h_shwr_length_unmatched, h_shwr_length_intime, intime_scale_factor, data_scale_factor, "",
	                                            "All Shower Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_all_shower_lengths_intime.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_longest_shwr_length_nue_cc,  h_longest_shwr_length_nue_cc_mixed,
	                                            h_longest_shwr_length_nue_cc_out_fv,
	                                            h_longest_shwr_length_numu_cc, h_longest_shwr_length_numu_cc_mixed,
	                                            h_longest_shwr_length_cosmic,  h_longest_shwr_length_nc,
	                                            h_longest_shwr_length_nc_pi0,  h_longest_shwr_length_other_mixed,
	                                            h_longest_shwr_length_unmatched, h_longest_shwr_length_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Longest Shower Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_longest_shower_lengths_intime.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_length_nue_cc,  h_leading_shwr_length_nue_cc_mixed,
	                                            h_leading_shwr_length_nue_cc_out_fv,
	                                            h_leading_shwr_length_numu_cc, h_leading_shwr_length_numu_cc_mixed,
	                                            h_leading_shwr_length_cosmic,  h_leading_shwr_length_nc,
	                                            h_leading_shwr_length_nc_pi0,  h_leading_shwr_length_other_mixed,
	                                            h_leading_shwr_length_unmatched, h_leading_shwr_length_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_leading_shower_lengths_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_shwr_length_nue_cc,  h_shwr_length_nue_cc_mixed,
	                                          h_shwr_length_nue_cc_out_fv,
	                                          h_shwr_length_numu_cc, h_shwr_length_numu_cc_mixed,
	                                          h_shwr_length_cosmic,  h_shwr_length_nc,
	                                          h_shwr_length_nc_pi0,  h_shwr_length_other_mixed,
	                                          h_shwr_length_unmatched, h_shwr_length_intime, intime_scale_factor,
	                                          h_shwr_length_data, data_scale_factor, "",
	                                          "All Shower Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_all_shower_lengths_data.pdf"));
	histogram_functions::PlotSimpleStackData (h_longest_shwr_length_nue_cc,  h_longest_shwr_length_nue_cc_mixed,
	                                          h_longest_shwr_length_nue_cc_out_fv,
	                                          h_longest_shwr_length_numu_cc, h_longest_shwr_length_numu_cc_mixed,
	                                          h_longest_shwr_length_cosmic,  h_longest_shwr_length_nc,
	                                          h_longest_shwr_length_nc_pi0,  h_longest_shwr_length_other_mixed,
	                                          h_longest_shwr_length_unmatched, h_longest_shwr_length_intime, intime_scale_factor,
	                                          h_longest_shwr_length_data, data_scale_factor, "",
	                                          "Longest Shower Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_longest_shower_lengths_data.pdf"));
	histogram_functions::PlotSimpleStackData (h_leading_shwr_length_nue_cc,  h_leading_shwr_length_nue_cc_mixed,
	                                          h_leading_shwr_length_nue_cc_out_fv,
	                                          h_leading_shwr_length_numu_cc, h_leading_shwr_length_numu_cc_mixed,
	                                          h_leading_shwr_length_cosmic,  h_leading_shwr_length_nc,
	                                          h_leading_shwr_length_nc_pi0,  h_leading_shwr_length_other_mixed,
	                                          h_leading_shwr_length_unmatched, h_leading_shwr_length_intime, intime_scale_factor,
	                                          h_leading_shwr_length_data, data_scale_factor, "",
	                                          "Leading Shower Lengths [cm]", "", Form("%s%s", file_locate_prefix, "post_leading_shower_lengths_data.pdf"));

	histogram_functions::PlotSimpleStack (h_leading_shwr_trk_length_nue_cc,  h_leading_shwr_trk_length_nue_cc_mixed,
	                                      h_leading_shwr_trk_length_nue_cc_out_fv,
	                                      h_leading_shwr_trk_length_numu_cc, h_leading_shwr_trk_length_numu_cc_mixed,
	                                      h_leading_shwr_trk_length_cosmic,  h_leading_shwr_trk_length_nc,
	                                      h_leading_shwr_trk_length_nc_pi0,  h_leading_shwr_trk_length_other_mixed,
	                                      h_leading_shwr_trk_length_unmatched, "",
	                                      "Longest Track / Leading Shower Lengths", "", Form("%s%s", file_locate_prefix, "post_leading_shower_trk_lengths.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_trk_length_nue_cc,  h_leading_shwr_trk_length_nue_cc_mixed,
	                                            h_leading_shwr_trk_length_nue_cc_out_fv,
	                                            h_leading_shwr_trk_length_numu_cc, h_leading_shwr_trk_length_numu_cc_mixed,
	                                            h_leading_shwr_trk_length_cosmic,  h_leading_shwr_trk_length_nc,
	                                            h_leading_shwr_trk_length_nc_pi0,  h_leading_shwr_trk_length_other_mixed,
	                                            h_leading_shwr_trk_length_unmatched, h_leading_shwr_trk_length_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Longest Track / Leading Shower Lengths", "",
	                                            Form("%s%s", file_locate_prefix, "post_leading_shower_trk_lengths_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_leading_shwr_trk_length_nue_cc,  h_leading_shwr_trk_length_nue_cc_mixed,
	                                          h_leading_shwr_trk_length_nue_cc_out_fv,
	                                          h_leading_shwr_trk_length_numu_cc, h_leading_shwr_trk_length_numu_cc_mixed,
	                                          h_leading_shwr_trk_length_cosmic,  h_leading_shwr_trk_length_nc,
	                                          h_leading_shwr_trk_length_nc_pi0,  h_leading_shwr_trk_length_other_mixed,
	                                          h_leading_shwr_trk_length_unmatched, h_leading_shwr_trk_length_intime, intime_scale_factor,
	                                          h_leading_shwr_trk_length_data, data_scale_factor, "",
	                                          "Longest Track / Leading Shower Lengths", "",
	                                          Form("%s%s", file_locate_prefix, "post_leading_shower_trk_lengths_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_leading_shwr_trk_length_nue_cc,  h_leading_shwr_trk_length_nue_cc_mixed,
		                                          h_leading_shwr_trk_length_nue_cc_out_fv,
		                                          h_leading_shwr_trk_length_numu_cc, h_leading_shwr_trk_length_numu_cc_mixed,
		                                          h_leading_shwr_trk_length_cosmic,  h_leading_shwr_trk_length_nc,
		                                          h_leading_shwr_trk_length_nc_pi0,  h_leading_shwr_trk_length_other_mixed,
		                                          h_leading_shwr_trk_length_unmatched, h_leading_shwr_trk_length_intime, scaled_intime_scale_factor,
		                                          h_leading_shwr_trk_length_data, data_scale_factor, "",
		                                          "Longest Track / Leading Shower Lengths (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_leading_shower_trk_lengths_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStackData (h_leading_shwr_trk_length_nue_cc_after,  h_leading_shwr_trk_length_nue_cc_mixed_after,
	                                          h_leading_shwr_trk_length_nue_cc_out_fv_after,
	                                          h_leading_shwr_trk_length_numu_cc_after, h_leading_shwr_trk_length_numu_cc_mixed_after,
	                                          h_leading_shwr_trk_length_cosmic_after,  h_leading_shwr_trk_length_nc_after,
	                                          h_leading_shwr_trk_length_nc_pi0_after,  h_leading_shwr_trk_length_other_mixed_after,
	                                          h_leading_shwr_trk_length_unmatched_after, h_leading_shwr_trk_length_intime_after, intime_scale_factor,
	                                          h_leading_shwr_trk_length_data_after, data_scale_factor, "",
	                                          "Longest Track / Leading Shower Lengths", "",
	                                          Form("%s%s", file_locate_prefix, "post_leading_shower_trk_lengths_data_after.pdf"));

	histogram_functions::PlotSimpleStack (h_longest_shwr_trk_length_nue_cc,  h_longest_shwr_trk_length_nue_cc_mixed,
	                                      h_longest_shwr_trk_length_nue_cc_out_fv,
	                                      h_longest_shwr_trk_length_numu_cc, h_longest_shwr_trk_length_numu_cc_mixed,
	                                      h_longest_shwr_trk_length_cosmic,  h_longest_shwr_trk_length_nc,
	                                      h_longest_shwr_trk_length_nc_pi0,  h_longest_shwr_trk_length_other_mixed,
	                                      h_longest_shwr_trk_length_unmatched, "",
	                                      "Longest Track / Longest Shower Lengths", "",
	                                      Form("%s%s", file_locate_prefix, "post_longest_shower_trk_lengths.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_longest_shwr_trk_length_nue_cc,  h_longest_shwr_trk_length_nue_cc_mixed,
	                                            h_longest_shwr_trk_length_nue_cc_out_fv,
	                                            h_longest_shwr_trk_length_numu_cc, h_longest_shwr_trk_length_numu_cc_mixed,
	                                            h_longest_shwr_trk_length_cosmic,  h_longest_shwr_trk_length_nc,
	                                            h_longest_shwr_trk_length_nc_pi0,  h_longest_shwr_trk_length_other_mixed,
	                                            h_longest_shwr_trk_length_unmatched, h_longest_shwr_trk_length_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Longest Track / Longest Shower Lengths", "",
	                                            Form("%s%s", file_locate_prefix, "post_longest_shower_trk_lengths_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_longest_shwr_trk_length_nue_cc,  h_longest_shwr_trk_length_nue_cc_mixed,
	                                          h_longest_shwr_trk_length_nue_cc_out_fv,
	                                          h_longest_shwr_trk_length_numu_cc, h_longest_shwr_trk_length_numu_cc_mixed,
	                                          h_longest_shwr_trk_length_cosmic,  h_longest_shwr_trk_length_nc,
	                                          h_longest_shwr_trk_length_nc_pi0,  h_longest_shwr_trk_length_other_mixed,
	                                          h_longest_shwr_trk_length_unmatched, h_longest_shwr_trk_length_intime, intime_scale_factor,
	                                          h_longest_shwr_trk_length_data, data_scale_factor, "",
	                                          "Longest Track / Longest Shower Lengths", "",
	                                          Form("%s%s", file_locate_prefix, "post_longest_shower_trk_lengths_data.pdf"));

	// histogram_functions::Plot1DHistogram (h_post_cuts_num_showers_purity, "Number of Reconstructed Showers per TPCO", "post_cuts_num_showers_purity.pdf"));
	// histogram_functions::Plot1DHistogram (h_post_open_angle_cuts_num_showers_purity, "Number of Reconstructed Showers per TPCO - After Open Angle Cut",
	//                                       "post_open_angle_cuts_num_showers_purity.pdf"));

	histogram_functions::PurityStack(h_post_cuts_num_showers_purity_qe, h_post_cuts_num_showers_purity_res, h_post_cuts_num_showers_purity_dis,
	                                 h_post_cuts_num_showers_purity_coh, h_post_cuts_num_showers_purity_mec,
	                                 "Num Reco Showers per TPCO", Form("%s%s", file_locate_prefix, "post_cuts_num_showers_purity.pdf"));
	histogram_functions::PurityStack(h_post_open_angle_cuts_num_showers_purity_qe, h_post_open_angle_cuts_num_showers_purity_res,
	                                 h_post_open_angle_cuts_num_showers_purity_dis, h_post_open_angle_cuts_num_showers_purity_coh,
	                                 h_post_open_angle_cuts_num_showers_purity_mec,
	                                 "Num Reco Showers per TPCO - After Open Angle Cut",
	                                 Form("%s%s", file_locate_prefix, "post_open_angle_cuts_num_showers_purity.pdf"));

	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nue_cc, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nue_cc_out_fv, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nue_cc_mixed, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_nue_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_numu_cc, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_numu_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_numu_cc_mixed, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_numu_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nc, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_nc.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nc_pi0, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_nc_pi0.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_cosmic, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_cosmic.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_other_mixed, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_other_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_unmatched, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_track_unmatched.pdf"));

	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nue_cc, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nue_cc_out_fv, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nue_cc_mixed, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_nue_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_numu_cc, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_numu_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_numu_cc_mixed, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_numu_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nc, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_nc.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nc_pi0, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_nc_pi0.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_cosmic, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_cosmic.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_other_mixed, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_other_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_unmatched, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_shower_unmatched.pdf"));

	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nue_cc, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nue_cc_out_fv, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nue_cc_mixed, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_nue_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_numu_cc, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_numu_cc.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_numu_cc_mixed, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_numu_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nc, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_nc.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nc_pi0, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_nc_pi0.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_cosmic, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_cosmic.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_other_mixed, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_other_mixed.pdf"));
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_unmatched, "", "Hits - Collection Plane", "Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_total_hits_leading_shower_unmatched.pdf"));

	histogram_functions::PlotSimpleStack (h_collection_hits_track_nue_cc,  h_collection_hits_track_nue_cc_mixed,
	                                      h_collection_hits_track_nue_cc_out_fv,
	                                      h_collection_hits_track_numu_cc, h_collection_hits_track_numu_cc_mixed,
	                                      h_collection_hits_track_cosmic,  h_collection_hits_track_nc,
	                                      h_collection_hits_track_nc_pi0,  h_collection_hits_track_other_mixed,
	                                      h_collection_hits_track_unmatched, "",
	                                      "Tracks Hits - Collection Plane", "", Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_tracks.pdf"));
	histogram_functions::PlotSimpleStack (h_collection_hits_shower_nue_cc,  h_collection_hits_shower_nue_cc_mixed,
	                                      h_collection_hits_shower_nue_cc_out_fv,
	                                      h_collection_hits_shower_numu_cc, h_collection_hits_shower_numu_cc_mixed,
	                                      h_collection_hits_shower_cosmic,  h_collection_hits_shower_nc,
	                                      h_collection_hits_shower_nc_pi0,  h_collection_hits_shower_other_mixed,
	                                      h_collection_hits_shower_unmatched, "",
	                                      "Showers Hits - Collection Plane", "", Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_showers.pdf"));
	histogram_functions::PlotSimpleStack (h_collection_hits_leading_shower_nue_cc,  h_collection_hits_leading_shower_nue_cc_mixed,
	                                      h_collection_hits_leading_shower_nue_cc_out_fv,
	                                      h_collection_hits_leading_shower_numu_cc, h_collection_hits_leading_shower_numu_cc_mixed,
	                                      h_collection_hits_leading_shower_cosmic,  h_collection_hits_leading_shower_nc,
	                                      h_collection_hits_leading_shower_nc_pi0,  h_collection_hits_leading_shower_other_mixed,
	                                      h_collection_hits_leading_shower_unmatched, "",
	                                      "Leading Shower Hits - Collection Plane", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_leading_shower.pdf"));
	histogram_functions::PlotSimpleStack (h_total_hits_leading_shower_nue_cc,  h_total_hits_leading_shower_nue_cc_mixed,
	                                      h_total_hits_leading_shower_nue_cc_out_fv,
	                                      h_total_hits_leading_shower_numu_cc, h_total_hits_leading_shower_numu_cc_mixed,
	                                      h_total_hits_leading_shower_cosmic,  h_total_hits_leading_shower_nc,
	                                      h_total_hits_leading_shower_nc_pi0,  h_total_hits_leading_shower_other_mixed,
	                                      h_total_hits_leading_shower_unmatched, "",
	                                      "Leading Shower Hits - All Planes", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_total_hits_leading_shower.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_collection_hits_track_nue_cc,  h_collection_hits_track_nue_cc_mixed,
	                                            h_collection_hits_track_nue_cc_out_fv,
	                                            h_collection_hits_track_numu_cc, h_collection_hits_track_numu_cc_mixed,
	                                            h_collection_hits_track_cosmic,  h_collection_hits_track_nc,
	                                            h_collection_hits_track_nc_pi0,  h_collection_hits_track_other_mixed,
	                                            h_collection_hits_track_unmatched, h_collection_hits_track_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Tracks Hits - Collection Plane", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_tracks_intime.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_collection_hits_shower_nue_cc,  h_collection_hits_shower_nue_cc_mixed,
	                                            h_collection_hits_shower_nue_cc_out_fv,
	                                            h_collection_hits_shower_numu_cc, h_collection_hits_shower_numu_cc_mixed,
	                                            h_collection_hits_shower_cosmic,  h_collection_hits_shower_nc,
	                                            h_collection_hits_shower_nc_pi0,  h_collection_hits_shower_other_mixed,
	                                            h_collection_hits_shower_unmatched, h_collection_hits_shower_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Showers Hits - Collection Plane", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_showers_intime.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_collection_hits_leading_shower_nue_cc,  h_collection_hits_leading_shower_nue_cc_mixed,
	                                            h_collection_hits_leading_shower_nue_cc_out_fv,
	                                            h_collection_hits_leading_shower_numu_cc, h_collection_hits_leading_shower_numu_cc_mixed,
	                                            h_collection_hits_leading_shower_cosmic,  h_collection_hits_leading_shower_nc,
	                                            h_collection_hits_leading_shower_nc_pi0,  h_collection_hits_leading_shower_other_mixed,
	                                            h_collection_hits_leading_shower_unmatched, h_collection_hits_leading_shower_intime,
	                                            intime_scale_factor, data_scale_factor,
	                                            "", "Leading Shower Hits - Collection Plane", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_leading_shower_intime.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_total_hits_leading_shower_nue_cc,  h_total_hits_leading_shower_nue_cc_mixed,
	                                            h_total_hits_leading_shower_nue_cc_out_fv,
	                                            h_total_hits_leading_shower_numu_cc, h_total_hits_leading_shower_numu_cc_mixed,
	                                            h_total_hits_leading_shower_cosmic,  h_total_hits_leading_shower_nc,
	                                            h_total_hits_leading_shower_nc_pi0,  h_total_hits_leading_shower_other_mixed,
	                                            h_total_hits_leading_shower_unmatched, h_total_hits_leading_shower_intime,
	                                            intime_scale_factor, data_scale_factor,
	                                            "", "Leading Shower Hits - All Planes", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_total_hits_leading_shower_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_collection_hits_track_nue_cc,  h_collection_hits_track_nue_cc_mixed,
	                                          h_collection_hits_track_nue_cc_out_fv,
	                                          h_collection_hits_track_numu_cc, h_collection_hits_track_numu_cc_mixed,
	                                          h_collection_hits_track_cosmic,  h_collection_hits_track_nc,
	                                          h_collection_hits_track_nc_pi0,  h_collection_hits_track_other_mixed,
	                                          h_collection_hits_track_unmatched, h_collection_hits_track_intime, intime_scale_factor,
	                                          h_collection_hits_track_data, data_scale_factor, "",
	                                          "Tracks Hits - Collection Plane", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_tracks_data.pdf"));
	histogram_functions::PlotSimpleStackData (h_collection_hits_shower_nue_cc,  h_collection_hits_shower_nue_cc_mixed,
	                                          h_collection_hits_shower_nue_cc_out_fv,
	                                          h_collection_hits_shower_numu_cc, h_collection_hits_shower_numu_cc_mixed,
	                                          h_collection_hits_shower_cosmic,  h_collection_hits_shower_nc,
	                                          h_collection_hits_shower_nc_pi0,  h_collection_hits_shower_other_mixed,
	                                          h_collection_hits_shower_unmatched, h_collection_hits_shower_intime, intime_scale_factor,
	                                          h_collection_hits_shower_data, data_scale_factor, "",
	                                          "Showers Hits - Collection Plane", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_showers_data.pdf"));
	histogram_functions::PlotSimpleStackData (h_collection_hits_leading_shower_nue_cc,  h_collection_hits_leading_shower_nue_cc_mixed,
	                                          h_collection_hits_leading_shower_nue_cc_out_fv,
	                                          h_collection_hits_leading_shower_numu_cc, h_collection_hits_leading_shower_numu_cc_mixed,
	                                          h_collection_hits_leading_shower_cosmic,  h_collection_hits_leading_shower_nc,
	                                          h_collection_hits_leading_shower_nc_pi0,  h_collection_hits_leading_shower_other_mixed,
	                                          h_collection_hits_leading_shower_unmatched,
	                                          h_collection_hits_leading_shower_intime, intime_scale_factor,
	                                          h_collection_hits_leading_shower_data, data_scale_factor,
	                                          "", "Leading Shower Hits - Collection Plane", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_leading_shower_data.pdf"));
	histogram_functions::PlotSimpleStackData (h_total_hits_leading_shower_nue_cc,  h_total_hits_leading_shower_nue_cc_mixed,
	                                          h_total_hits_leading_shower_nue_cc_out_fv,
	                                          h_total_hits_leading_shower_numu_cc, h_total_hits_leading_shower_numu_cc_mixed,
	                                          h_total_hits_leading_shower_cosmic,  h_total_hits_leading_shower_nc,
	                                          h_total_hits_leading_shower_nc_pi0,  h_total_hits_leading_shower_other_mixed,
	                                          h_total_hits_leading_shower_unmatched, h_total_hits_leading_shower_intime, intime_scale_factor,
	                                          h_total_hits_leading_shower_data, data_scale_factor,
	                                          "", "Leading Shower Hits - All Planes", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_total_hits_leading_shower_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_collection_hits_leading_shower_nue_cc,  h_collection_hits_leading_shower_nue_cc_mixed,
		                                          h_collection_hits_leading_shower_nue_cc_out_fv,
		                                          h_collection_hits_leading_shower_numu_cc, h_collection_hits_leading_shower_numu_cc_mixed,
		                                          h_collection_hits_leading_shower_cosmic,  h_collection_hits_leading_shower_nc,
		                                          h_collection_hits_leading_shower_nc_pi0,  h_collection_hits_leading_shower_other_mixed,
		                                          h_collection_hits_leading_shower_unmatched,
		                                          h_collection_hits_leading_shower_intime, scaled_intime_scale_factor,
		                                          h_collection_hits_leading_shower_data, data_scale_factor,
		                                          "", "Leading Shower Hits - Collection Plane (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_leading_shower_scaled_data.pdf"));
		histogram_functions::PlotSimpleStackData (h_total_hits_leading_shower_nue_cc,  h_total_hits_leading_shower_nue_cc_mixed,
		                                          h_total_hits_leading_shower_nue_cc_out_fv,
		                                          h_total_hits_leading_shower_numu_cc, h_total_hits_leading_shower_numu_cc_mixed,
		                                          h_total_hits_leading_shower_cosmic,  h_total_hits_leading_shower_nc,
		                                          h_total_hits_leading_shower_nc_pi0,  h_total_hits_leading_shower_other_mixed,
		                                          h_total_hits_leading_shower_unmatched, h_total_hits_leading_shower_intime, scaled_intime_scale_factor,
		                                          h_total_hits_leading_shower_data, data_scale_factor,
		                                          "", "Leading Shower Hits - All Planes", "",
		                                          Form("%s%s", file_locate_prefix, "post_cuts_total_hits_leading_shower_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStackData (h_collection_hits_track_nue_cc_after,  h_collection_hits_track_nue_cc_mixed_after,
	                                          h_collection_hits_track_nue_cc_out_fv_after,
	                                          h_collection_hits_track_numu_cc_after, h_collection_hits_track_numu_cc_mixed_after,
	                                          h_collection_hits_track_cosmic_after,  h_collection_hits_track_nc_after,
	                                          h_collection_hits_track_nc_pi0_after,  h_collection_hits_track_other_mixed_after,
	                                          h_collection_hits_track_unmatched_after, h_collection_hits_track_intime_after, intime_scale_factor,
	                                          h_collection_hits_track_data_after, data_scale_factor, "",
	                                          "Tracks Hits - Collection Plane", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_tracks_data_after.pdf"));
	histogram_functions::PlotSimpleStackData (h_collection_hits_shower_nue_cc_after,  h_collection_hits_shower_nue_cc_mixed_after,
	                                          h_collection_hits_shower_nue_cc_out_fv_after,
	                                          h_collection_hits_shower_numu_cc_after, h_collection_hits_shower_numu_cc_mixed_after,
	                                          h_collection_hits_shower_cosmic_after,  h_collection_hits_shower_nc_after,
	                                          h_collection_hits_shower_nc_pi0_after,  h_collection_hits_shower_other_mixed_after,
	                                          h_collection_hits_shower_unmatched_after, h_collection_hits_shower_intime_after, intime_scale_factor,
	                                          h_collection_hits_shower_data_after, data_scale_factor, "",
	                                          "Showers Hits - Collection Plane", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_showers_data_after.pdf"));
	histogram_functions::PlotSimpleStackData (h_collection_hits_leading_shower_nue_cc_after,  h_collection_hits_leading_shower_nue_cc_mixed_after,
	                                          h_collection_hits_leading_shower_nue_cc_out_fv_after,
	                                          h_collection_hits_leading_shower_numu_cc_after, h_collection_hits_leading_shower_numu_cc_mixed_after,
	                                          h_collection_hits_leading_shower_cosmic_after,  h_collection_hits_leading_shower_nc_after,
	                                          h_collection_hits_leading_shower_nc_pi0_after,  h_collection_hits_leading_shower_other_mixed_after,
	                                          h_collection_hits_leading_shower_unmatched_after,
	                                          h_collection_hits_leading_shower_intime_after, intime_scale_factor,
	                                          h_collection_hits_leading_shower_data_after, data_scale_factor,
	                                          "", "Leading Shower Hits - Collection Plane", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_collection_hits_leading_shower_data_after.pdf"));
	histogram_functions::PlotSimpleStackData (h_total_hits_leading_shower_nue_cc_after,  h_total_hits_leading_shower_nue_cc_mixed_after,
	                                          h_total_hits_leading_shower_nue_cc_out_fv_after,
	                                          h_total_hits_leading_shower_numu_cc_after, h_total_hits_leading_shower_numu_cc_mixed_after,
	                                          h_total_hits_leading_shower_cosmic_after,  h_total_hits_leading_shower_nc_after,
	                                          h_total_hits_leading_shower_nc_pi0_after,  h_total_hits_leading_shower_other_mixed_after,
	                                          h_total_hits_leading_shower_unmatched_after, h_total_hits_leading_shower_intime_after, intime_scale_factor,
	                                          h_total_hits_leading_shower_data_after, data_scale_factor,
	                                          "", "Leading Shower Hits - All Planes", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_total_hits_leading_shower_data_after.pdf"));

	histogram_functions::PlotSimpleStack (h_pre_cut_collection_hits_track_nue_cc,  h_pre_cut_collection_hits_track_nue_cc_mixed,
	                                      h_pre_cut_collection_hits_track_nue_cc_out_fv,
	                                      h_pre_cut_collection_hits_track_numu_cc, h_pre_cut_collection_hits_track_numu_cc_mixed,
	                                      h_pre_cut_collection_hits_track_cosmic,  h_pre_cut_collection_hits_track_nc,
	                                      h_pre_cut_collection_hits_track_nc_pi0,  h_pre_cut_collection_hits_track_other_mixed,
	                                      h_pre_cut_collection_hits_track_unmatched, "",
	                                      "Tracks Hits - Collection Plane", "",
	                                      Form("%s%s", file_locate_prefix, "pre_hit_cut_collection_hits_tracks.pdf"));
	histogram_functions::PlotSimpleStack (h_pre_cut_collection_hits_shower_nue_cc,  h_pre_cut_collection_hits_shower_nue_cc_mixed,
	                                      h_pre_cut_collection_hits_shower_nue_cc_out_fv,
	                                      h_pre_cut_collection_hits_shower_numu_cc, h_pre_cut_collection_hits_shower_numu_cc_mixed,
	                                      h_pre_cut_collection_hits_shower_cosmic,  h_pre_cut_collection_hits_shower_nc,
	                                      h_pre_cut_collection_hits_shower_nc_pi0,  h_pre_cut_collection_hits_shower_other_mixed,
	                                      h_pre_cut_collection_hits_shower_unmatched, "",
	                                      "Showers Hits - Collection Plane", "",
	                                      Form("%s%s", file_locate_prefix, "pre_hit_cut_collection_hits_showers.pdf"));
	histogram_functions::PlotSimpleStack (h_pre_cut_collection_hits_leading_shower_nue_cc,  h_pre_cut_collection_hits_leading_shower_nue_cc_mixed,
	                                      h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv,
	                                      h_pre_cut_collection_hits_leading_shower_numu_cc, h_pre_cut_collection_hits_leading_shower_numu_cc_mixed,
	                                      h_pre_cut_collection_hits_leading_shower_cosmic,  h_pre_cut_collection_hits_leading_shower_nc,
	                                      h_pre_cut_collection_hits_leading_shower_nc_pi0,  h_pre_cut_collection_hits_leading_shower_other_mixed,
	                                      h_pre_cut_collection_hits_leading_shower_unmatched, "",
	                                      "Leading Shower Hits - Collection Plane", "",
	                                      Form("%s%s", file_locate_prefix, "pre_hit_cut_collection_hits_leading_shower.pdf"));
	histogram_functions::PlotSimpleStack (h_pre_cut_total_hits_leading_shower_nue_cc,  h_pre_cut_total_hits_leading_shower_nue_cc_mixed,
	                                      h_pre_cut_total_hits_leading_shower_nue_cc_out_fv,
	                                      h_pre_cut_total_hits_leading_shower_numu_cc, h_pre_cut_total_hits_leading_shower_numu_cc_mixed,
	                                      h_pre_cut_total_hits_leading_shower_cosmic,  h_pre_cut_total_hits_leading_shower_nc,
	                                      h_pre_cut_total_hits_leading_shower_nc_pi0,  h_pre_cut_total_hits_leading_shower_other_mixed,
	                                      h_pre_cut_total_hits_leading_shower_unmatched, "",
	                                      "Leading Shower Hits - All Planes", "",
	                                      Form("%s%s", file_locate_prefix, "pre_hit_cut_total_hits_leading_shower.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_pre_cut_collection_hits_track_nue_cc,  h_pre_cut_collection_hits_track_nue_cc_mixed,
	                                            h_pre_cut_collection_hits_track_nue_cc_out_fv,
	                                            h_pre_cut_collection_hits_track_numu_cc, h_pre_cut_collection_hits_track_numu_cc_mixed,
	                                            h_pre_cut_collection_hits_track_cosmic,  h_pre_cut_collection_hits_track_nc,
	                                            h_pre_cut_collection_hits_track_nc_pi0,  h_pre_cut_collection_hits_track_other_mixed,
	                                            h_pre_cut_collection_hits_track_unmatched, h_pre_cut_collection_hits_track_intime,
	                                            intime_scale_factor, data_scale_factor,
	                                            "", "Tracks Hits - Collection Plane", "",
	                                            Form("%s%s", file_locate_prefix, "pre_hit_cut_collection_hits_tracks_intime.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_pre_cut_collection_hits_shower_nue_cc,  h_pre_cut_collection_hits_shower_nue_cc_mixed,
	                                            h_pre_cut_collection_hits_shower_nue_cc_out_fv,
	                                            h_pre_cut_collection_hits_shower_numu_cc, h_pre_cut_collection_hits_shower_numu_cc_mixed,
	                                            h_pre_cut_collection_hits_shower_cosmic,  h_pre_cut_collection_hits_shower_nc,
	                                            h_pre_cut_collection_hits_shower_nc_pi0,  h_pre_cut_collection_hits_shower_other_mixed,
	                                            h_pre_cut_collection_hits_shower_unmatched, h_pre_cut_collection_hits_shower_intime,
	                                            intime_scale_factor, data_scale_factor,
	                                            "", "Showers Hits - Collection Plane", "",
	                                            Form("%s%s", file_locate_prefix, "pre_hit_cut_collection_hits_showers_intime.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_pre_cut_collection_hits_leading_shower_nue_cc,  h_pre_cut_collection_hits_leading_shower_nue_cc_mixed,
	                                            h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv,
	                                            h_pre_cut_collection_hits_leading_shower_numu_cc, h_pre_cut_collection_hits_leading_shower_numu_cc_mixed,
	                                            h_pre_cut_collection_hits_leading_shower_cosmic,  h_pre_cut_collection_hits_leading_shower_nc,
	                                            h_pre_cut_collection_hits_leading_shower_nc_pi0,  h_pre_cut_collection_hits_leading_shower_other_mixed,
	                                            h_pre_cut_collection_hits_leading_shower_unmatched, h_pre_cut_collection_hits_leading_shower_intime,
	                                            intime_scale_factor, data_scale_factor,
	                                            "", "Leading Shower Hits - Collection Plane", "",
	                                            Form("%s%s", file_locate_prefix, "pre_hit_cut_collection_hits_leading_shower_intime.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_pre_cut_total_hits_leading_shower_nue_cc,  h_pre_cut_total_hits_leading_shower_nue_cc_mixed,
	                                            h_pre_cut_total_hits_leading_shower_nue_cc_out_fv,
	                                            h_pre_cut_total_hits_leading_shower_numu_cc, h_pre_cut_total_hits_leading_shower_numu_cc_mixed,
	                                            h_pre_cut_total_hits_leading_shower_cosmic,  h_pre_cut_total_hits_leading_shower_nc,
	                                            h_pre_cut_total_hits_leading_shower_nc_pi0,  h_pre_cut_total_hits_leading_shower_other_mixed,
	                                            h_pre_cut_total_hits_leading_shower_unmatched, h_pre_cut_total_hits_leading_shower_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Hits - All Planes", "",
	                                            Form("%s%s", file_locate_prefix, "pre_hit_cut_total_hits_leading_shower_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_pre_cut_collection_hits_track_nue_cc,  h_pre_cut_collection_hits_track_nue_cc_mixed,
	                                          h_pre_cut_collection_hits_track_nue_cc_out_fv,
	                                          h_pre_cut_collection_hits_track_numu_cc, h_pre_cut_collection_hits_track_numu_cc_mixed,
	                                          h_pre_cut_collection_hits_track_cosmic,  h_pre_cut_collection_hits_track_nc,
	                                          h_pre_cut_collection_hits_track_nc_pi0,  h_pre_cut_collection_hits_track_other_mixed,
	                                          h_pre_cut_collection_hits_track_unmatched,
	                                          h_pre_cut_collection_hits_track_intime, intime_scale_factor,
	                                          h_pre_cut_collection_hits_track_data, data_scale_factor,
	                                          "", "Tracks Hits - Collection Plane", "",
	                                          Form("%s%s", file_locate_prefix, "pre_hit_cut_collection_hits_tracks_data.pdf"));
	histogram_functions::PlotSimpleStackData (h_pre_cut_collection_hits_shower_nue_cc,  h_pre_cut_collection_hits_shower_nue_cc_mixed,
	                                          h_pre_cut_collection_hits_shower_nue_cc_out_fv,
	                                          h_pre_cut_collection_hits_shower_numu_cc, h_pre_cut_collection_hits_shower_numu_cc_mixed,
	                                          h_pre_cut_collection_hits_shower_cosmic,  h_pre_cut_collection_hits_shower_nc,
	                                          h_pre_cut_collection_hits_shower_nc_pi0,  h_pre_cut_collection_hits_shower_other_mixed,
	                                          h_pre_cut_collection_hits_shower_unmatched,
	                                          h_pre_cut_collection_hits_shower_intime, intime_scale_factor,
	                                          h_pre_cut_collection_hits_shower_data, data_scale_factor,
	                                          "", "Showers Hits - Collection Plane", "",
	                                          Form("%s%s", file_locate_prefix, "pre_hit_cut_collection_hits_showers_data.pdf"));
	histogram_functions::PlotSimpleStackData (h_pre_cut_collection_hits_leading_shower_nue_cc,  h_pre_cut_collection_hits_leading_shower_nue_cc_mixed,
	                                          h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv,
	                                          h_pre_cut_collection_hits_leading_shower_numu_cc, h_pre_cut_collection_hits_leading_shower_numu_cc_mixed,
	                                          h_pre_cut_collection_hits_leading_shower_cosmic,  h_pre_cut_collection_hits_leading_shower_nc,
	                                          h_pre_cut_collection_hits_leading_shower_nc_pi0,  h_pre_cut_collection_hits_leading_shower_other_mixed,
	                                          h_pre_cut_collection_hits_leading_shower_unmatched,
	                                          h_pre_cut_collection_hits_leading_shower_intime, intime_scale_factor,
	                                          h_pre_cut_collection_hits_leading_shower_data, data_scale_factor,
	                                          "", "Leading Shower Hits - Collection Plane", "",
	                                          Form("%s%s", file_locate_prefix, "pre_hit_cut_collection_hits_leading_shower_data.pdf"));
	histogram_functions::PlotSimpleStackData (h_pre_cut_total_hits_leading_shower_nue_cc,  h_pre_cut_total_hits_leading_shower_nue_cc_mixed,
	                                          h_pre_cut_total_hits_leading_shower_nue_cc_out_fv,
	                                          h_pre_cut_total_hits_leading_shower_numu_cc, h_pre_cut_total_hits_leading_shower_numu_cc_mixed,
	                                          h_pre_cut_total_hits_leading_shower_cosmic,  h_pre_cut_total_hits_leading_shower_nc,
	                                          h_pre_cut_total_hits_leading_shower_nc_pi0,  h_pre_cut_total_hits_leading_shower_other_mixed,
	                                          h_pre_cut_total_hits_leading_shower_unmatched,
	                                          h_pre_cut_total_hits_leading_shower_intime, intime_scale_factor,
	                                          h_pre_cut_total_hits_leading_shower_data, data_scale_factor,
	                                          "", "Leading Shower Hits - All Planes", "",
	                                          Form("%s%s", file_locate_prefix, "pre_hit_cut_total_hits_leading_shower_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_pre_cut_total_hits_leading_shower_nue_cc,  h_pre_cut_total_hits_leading_shower_nue_cc_mixed,
		                                          h_pre_cut_total_hits_leading_shower_nue_cc_out_fv,
		                                          h_pre_cut_total_hits_leading_shower_numu_cc, h_pre_cut_total_hits_leading_shower_numu_cc_mixed,
		                                          h_pre_cut_total_hits_leading_shower_cosmic,  h_pre_cut_total_hits_leading_shower_nc,
		                                          h_pre_cut_total_hits_leading_shower_nc_pi0,  h_pre_cut_total_hits_leading_shower_other_mixed,
		                                          h_pre_cut_total_hits_leading_shower_unmatched,
		                                          h_pre_cut_total_hits_leading_shower_intime, scaled_intime_scale_factor,
		                                          h_pre_cut_total_hits_leading_shower_data, data_scale_factor,
		                                          "", "Leading Shower Hits - All Planes (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "pre_hit_cut_total_hits_leading_shower_scaled_data.pdf"));
	}

	int total_total_hits = 0;
	double bin_val_total_hits = 0;
	for(int i = 0; i < h_pre_cut_total_hits_leading_shower_nue_cc->GetNbinsX()+2; i++)
	{
		const double bin_val = h_pre_cut_total_hits_leading_shower_nue_cc->GetBinContent(i);
		std::cout << "Hits: " << h_pre_cut_total_hits_leading_shower_nue_cc->GetBinLowEdge(i) << ", Num Signal: " <<  bin_val << std::endl;
		if(h_pre_cut_total_hits_leading_shower_nue_cc->GetBinLowEdge(i) < 200)
		{
			bin_val_total_hits += bin_val;
		}
		total_total_hits += bin_val;
	}
	std::cout << "Num Signals: " << total_total_hits << std::endl;
	std::cout << "Removed: " << bin_val_total_hits << std::endl;

	histogram_functions::PlotSimpleStackData (h_pre_cut_total_hits_leading_shower_nue_cc,  h_pre_cut_total_hits_leading_shower_nue_cc_mixed,
	                                          h_pre_cut_total_hits_leading_shower_nue_cc_out_fv,
	                                          h_pre_cut_total_hits_leading_shower_numu_cc, h_pre_cut_total_hits_leading_shower_numu_cc_mixed,
	                                          h_pre_cut_total_hits_leading_shower_cosmic,  h_pre_cut_total_hits_leading_shower_nc,
	                                          h_pre_cut_total_hits_leading_shower_nc_pi0,  h_pre_cut_total_hits_leading_shower_other_mixed,
	                                          h_pre_cut_total_hits_leading_shower_unmatched,
	                                          h_pre_cut_total_hits_leading_shower_intime,intime_scale_factor,
	                                          h_pre_cut_total_hits_leading_shower_data, data_scale_factor,
	                                          0.73, 0.98, 0.98, 0.50, true,
	                                          "", "Leading Shower Hits - All Planes", "",
	                                          Form("%s%s", file_locate_prefix, "pre_hit_cut_total_hits_leading_shower_data_logy.pdf"));

	histogram_functions::Plot1DHistogram (h_pre_cut_total_hits_leading_shower_nue_cc, "Signal Events Leading Shower Hits - All Planes",
	                                      Form("%s%s", file_locate_prefix, "pre_hit_cut_total_hits_leading_shower_signal_only_data.pdf"));

	histogram_functions::Plot2DHistogram (h_ele_eng_total_hits, "", "Hits - All Planes", "True Electron Energy [GeV]",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_ele_eng_total_hits.pdf"));
	histogram_functions::Plot2DHistogram (h_ele_eng_colleciton_hits, "", "Hits - Collection Plane", "True Electron Energy [GeV]",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_ele_eng_collection_hits.pdf"));
	histogram_functions::Plot2DHistogram (h_nu_eng_total_hits, "", "Hits - All Planes", "True Neutrino Energy [GeV]",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_nu_eng_total_hits.pdf"));
	histogram_functions::Plot2DHistogram (h_nu_eng_collection_hits, "", "Hits - Collection Planes", "True Neutrino Energy [GeV]",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_nu_eng_colleciton_hits.pdf"));
	histogram_functions::PlotSimpleStack (h_ele_cos_theta_last_nue_cc,  h_ele_cos_theta_last_nue_cc_mixed,
	                                      h_ele_cos_theta_last_nue_cc_out_fv,
	                                      h_ele_cos_theta_last_numu_cc, h_ele_cos_theta_last_numu_cc_mixed,
	                                      h_ele_cos_theta_last_cosmic,  h_ele_cos_theta_last_nc,
	                                      h_ele_cos_theta_last_nc_pi0,  h_ele_cos_theta_last_other_mixed,
	                                      h_ele_cos_theta_last_unmatched, 0.15, 0.40, 0.98, 0.50, "",
	                                      "Leading Shower Cos(#theta)", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_cos_theta_last.pdf"));
	histogram_functions::PlotSimpleStack (h_ele_cos_theta_nue_cc,  h_ele_cos_theta_nue_cc_mixed,
	                                      h_ele_cos_theta_nue_cc_out_fv,
	                                      h_ele_cos_theta_numu_cc, h_ele_cos_theta_numu_cc_mixed,
	                                      h_ele_cos_theta_cosmic,  h_ele_cos_theta_nc,
	                                      h_ele_cos_theta_nc_pi0,  h_ele_cos_theta_other_mixed,
	                                      h_ele_cos_theta_unmatched, 0.15, 0.40, 0.98, 0.50, "",
	                                      "(Before Open Angle Cut) Leading Shower Cos(#theta)", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_cos_theta.pdf"));

	histogram_functions::PlotSimpleStackInTime (h_ele_cos_theta_last_nue_cc,  h_ele_cos_theta_last_nue_cc_mixed,
	                                            h_ele_cos_theta_last_nue_cc_out_fv,
	                                            h_ele_cos_theta_last_numu_cc, h_ele_cos_theta_last_numu_cc_mixed,
	                                            h_ele_cos_theta_last_cosmic,  h_ele_cos_theta_last_nc,
	                                            h_ele_cos_theta_last_nc_pi0,  h_ele_cos_theta_last_other_mixed,
	                                            h_ele_cos_theta_last_unmatched, h_ele_cos_theta_last_intime,
	                                            intime_scale_factor, data_scale_factor,
	                                            0.15, 0.40, 0.98, 0.50, "",
	                                            "Leading Shower Cos(#theta)", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_cos_theta_last_intime.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_ele_cos_theta_nue_cc,  h_ele_cos_theta_nue_cc_mixed,
	                                            h_ele_cos_theta_nue_cc_out_fv,
	                                            h_ele_cos_theta_numu_cc, h_ele_cos_theta_numu_cc_mixed,
	                                            h_ele_cos_theta_cosmic,  h_ele_cos_theta_nc,
	                                            h_ele_cos_theta_nc_pi0,  h_ele_cos_theta_other_mixed,
	                                            h_ele_cos_theta_unmatched, h_ele_cos_theta_intime,
	                                            intime_scale_factor, data_scale_factor,
	                                            0.15, 0.40, 0.98, 0.50, "",
	                                            "(Before Open Angle Cut) Leading Shower Cos(#theta)", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_cos_theta_intime.pdf"));

	histogram_functions::PlotSimpleStackData (h_ele_cos_theta_last_nue_cc,  h_ele_cos_theta_last_nue_cc_mixed,
	                                          h_ele_cos_theta_last_nue_cc_out_fv,
	                                          h_ele_cos_theta_last_numu_cc, h_ele_cos_theta_last_numu_cc_mixed,
	                                          h_ele_cos_theta_last_cosmic,  h_ele_cos_theta_last_nc,
	                                          h_ele_cos_theta_last_nc_pi0,  h_ele_cos_theta_last_other_mixed,
	                                          h_ele_cos_theta_last_unmatched, h_ele_cos_theta_last_intime, intime_scale_factor,
	                                          h_ele_cos_theta_last_data, data_scale_factor,
	                                          0.15, 0.40, 0.98, 0.50, false, "",
	                                          "Leading Shower Cos(#theta)", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_cos_theta_last_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_ele_cos_theta_last_nue_cc,  h_ele_cos_theta_last_nue_cc_mixed,
		                                          h_ele_cos_theta_last_nue_cc_out_fv,
		                                          h_ele_cos_theta_last_numu_cc, h_ele_cos_theta_last_numu_cc_mixed,
		                                          h_ele_cos_theta_last_cosmic,  h_ele_cos_theta_last_nc,
		                                          h_ele_cos_theta_last_nc_pi0,  h_ele_cos_theta_last_other_mixed,
		                                          h_ele_cos_theta_last_unmatched, h_ele_cos_theta_last_intime, scaled_intime_scale_factor,
		                                          h_ele_cos_theta_last_data, data_scale_factor,
		                                          0.15, 0.40, 0.98, 0.50, false, "",
		                                          "Leading Shower Cos(#theta) (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_cos_theta_last_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStackData (h_ele_cos_theta_last_trans_nue_cc,    h_ele_cos_theta_last_trans_nue_cc_mixed,
	                                          h_ele_cos_theta_last_trans_nue_cc_out_fv,
	                                          h_ele_cos_theta_last_trans_numu_cc,   h_ele_cos_theta_last_trans_numu_cc_mixed,
	                                          h_ele_cos_theta_last_trans_cosmic,    h_ele_cos_theta_last_trans_nc,
	                                          h_ele_cos_theta_last_trans_nc_pi0,    h_ele_cos_theta_last_trans_other_mixed,
	                                          h_ele_cos_theta_last_trans_unmatched, h_ele_cos_theta_last_trans_intime, intime_scale_factor,
	                                          h_ele_cos_theta_last_trans_data, data_scale_factor,
	                                          0.15, 0.40, 0.98, 0.50, false, "",
	                                          "(NuMI Transform) Leading Shower Cos(#theta)", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_cos_theta_last_trans_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_ele_cos_theta_nue_cc,  h_ele_cos_theta_nue_cc_mixed,
	                                          h_ele_cos_theta_nue_cc_out_fv,
	                                          h_ele_cos_theta_numu_cc, h_ele_cos_theta_numu_cc_mixed,
	                                          h_ele_cos_theta_cosmic,  h_ele_cos_theta_nc,
	                                          h_ele_cos_theta_nc_pi0,  h_ele_cos_theta_other_mixed,
	                                          h_ele_cos_theta_unmatched, h_ele_cos_theta_intime, intime_scale_factor,
	                                          h_ele_cos_theta_data, data_scale_factor,
	                                          0.15, 0.40, 0.98, 0.50, false, "",
	                                          "(Before Open Angle Cut) Leading Shower Cos(#theta)", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_cos_theta_data.pdf"));

	histogram_functions::PlotSimpleStack (h_ele_pfp_momentum_nue_cc,  h_ele_pfp_momentum_nue_cc_mixed,
	                                      h_ele_pfp_momentum_nue_cc_out_fv,
	                                      h_ele_pfp_momentum_numu_cc, h_ele_pfp_momentum_numu_cc_mixed,
	                                      h_ele_pfp_momentum_cosmic,  h_ele_pfp_momentum_nc,
	                                      h_ele_pfp_momentum_nc_pi0,  h_ele_pfp_momentum_other_mixed,
	                                      h_ele_pfp_momentum_unmatched, "",
	                                      "Leading Shower Momentum [GeV]", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_momentum.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_momentum_nue_cc,  h_ele_pfp_momentum_nue_cc_mixed,
	                                            h_ele_pfp_momentum_nue_cc_out_fv,
	                                            h_ele_pfp_momentum_numu_cc, h_ele_pfp_momentum_numu_cc_mixed,
	                                            h_ele_pfp_momentum_cosmic,  h_ele_pfp_momentum_nc,
	                                            h_ele_pfp_momentum_nc_pi0,  h_ele_pfp_momentum_other_mixed,
	                                            h_ele_pfp_momentum_unmatched, h_ele_pfp_momentum_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Momentum [GeV]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_momentum_intime.pdf"));
	histogram_functions::PlotSimpleStackDataMomentumRebin (h_ele_pfp_momentum_nue_cc,  h_ele_pfp_momentum_nue_cc_mixed,
	                                                       h_ele_pfp_momentum_nue_cc_out_fv,
	                                                       h_ele_pfp_momentum_numu_cc, h_ele_pfp_momentum_numu_cc_mixed,
	                                                       h_ele_pfp_momentum_cosmic,  h_ele_pfp_momentum_nc,
	                                                       h_ele_pfp_momentum_nc_pi0,  h_ele_pfp_momentum_other_mixed,
	                                                       h_ele_pfp_momentum_unmatched, h_ele_pfp_momentum_intime, intime_scale_factor,
	                                                       h_ele_pfp_momentum_data, data_scale_factor,
	                                                       "", "Leading Shower Momentum [GeV]", "",
	                                                       Form("%s%s", file_locate_prefix, "post_cuts_leading_momentum_data.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin (h_ele_pfp_momentum_no_track_nue_cc,  h_ele_pfp_momentum_no_track_nue_cc_mixed,
	                                                       h_ele_pfp_momentum_no_track_nue_cc_out_fv,
	                                                       h_ele_pfp_momentum_no_track_numu_cc, h_ele_pfp_momentum_no_track_numu_cc_mixed,
	                                                       h_ele_pfp_momentum_no_track_cosmic,  h_ele_pfp_momentum_no_track_nc,
	                                                       h_ele_pfp_momentum_no_track_nc_pi0,  h_ele_pfp_momentum_no_track_other_mixed,
	                                                       h_ele_pfp_momentum_no_track_unmatched, h_ele_pfp_momentum_no_track_intime, intime_scale_factor,
	                                                       h_ele_pfp_momentum_no_track_data, data_scale_factor,
	                                                       "", "Leading Shower Momentum [GeV]", "",
	                                                       Form("%s%s", file_locate_prefix, "post_cuts_leading_momentum_no_track_data.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin (h_ele_pfp_momentum_has_track_nue_cc,  h_ele_pfp_momentum_has_track_nue_cc_mixed,
	                                                       h_ele_pfp_momentum_has_track_nue_cc_out_fv,
	                                                       h_ele_pfp_momentum_has_track_numu_cc, h_ele_pfp_momentum_has_track_numu_cc_mixed,
	                                                       h_ele_pfp_momentum_has_track_cosmic,  h_ele_pfp_momentum_has_track_nc,
	                                                       h_ele_pfp_momentum_has_track_nc_pi0,  h_ele_pfp_momentum_has_track_other_mixed,
	                                                       h_ele_pfp_momentum_has_track_unmatched, h_ele_pfp_momentum_has_track_intime, intime_scale_factor,
	                                                       h_ele_pfp_momentum_has_track_data, data_scale_factor,
	                                                       "", "Leading Shower Momentum [GeV]", "",
	                                                       Form("%s%s", file_locate_prefix, "post_cuts_leading_momentum_has_track_data.pdf"));

	int tracks_num_signals_2 = 0;
	int tracks_num_total_2 = 0;
	for(int i = 0; i < h_ele_pfp_momentum_has_track_nue_cc->GetNbinsX() + 1; i++)
	{
		tracks_num_signals += h_ele_pfp_momentum_has_track_nue_cc->GetBinContent(i) * data_scale_factor;
		tracks_num_total += ((h_ele_pfp_momentum_has_track_nue_cc->GetBinContent(i) + h_ele_pfp_momentum_has_track_nue_cc_mixed->GetBinContent(i) +
		                      h_ele_pfp_momentum_has_track_nue_cc_out_fv->GetBinContent(i) + h_ele_pfp_momentum_has_track_numu_cc->GetBinContent(i) +
		                      h_ele_pfp_momentum_has_track_numu_cc_mixed->GetBinContent(i) +
		                      h_ele_pfp_momentum_has_track_cosmic->GetBinContent(i) + h_ele_pfp_momentum_has_track_nc->GetBinContent(i) +
		                      h_ele_pfp_momentum_has_track_nc_pi0->GetBinContent(i) +
		                      h_ele_pfp_momentum_has_track_other_mixed->GetBinContent(i) +
		                      h_ele_pfp_momentum_has_track_unmatched->GetBinContent(i)) * data_scale_factor +
		                     (h_ele_pfp_momentum_has_track_intime->GetBinContent(i) * intime_scale_factor));
	}
	std::cout << "--- At the end of the selection we have: " << tracks_num_signals_2 << " signal events with >= 1 Reco Track" << std::endl;
	std::cout << "--- And a total of : " << tracks_num_total_2 <<  " events with >= 1 Reco Track" << std::endl;

	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackDataMomentumRebin (h_ele_pfp_momentum_nue_cc,  h_ele_pfp_momentum_nue_cc_mixed,
		                                                       h_ele_pfp_momentum_nue_cc_out_fv,
		                                                       h_ele_pfp_momentum_numu_cc, h_ele_pfp_momentum_numu_cc_mixed,
		                                                       h_ele_pfp_momentum_cosmic,  h_ele_pfp_momentum_nc,
		                                                       h_ele_pfp_momentum_nc_pi0,  h_ele_pfp_momentum_other_mixed,
		                                                       h_ele_pfp_momentum_unmatched, h_ele_pfp_momentum_intime, scaled_intime_scale_factor,
		                                                       h_ele_pfp_momentum_data, data_scale_factor,
		                                                       "", "Leading Shower Momentum [GeV] (Scaled)", "",
		                                                       Form("%s%s", file_locate_prefix, "post_cuts_leading_momentum_scaled_data.pdf"));

		histogram_functions::PlotSimpleStackDataMomentumRebin (h_ele_pfp_momentum_no_track_nue_cc,  h_ele_pfp_momentum_no_track_nue_cc_mixed,
		                                                       h_ele_pfp_momentum_no_track_nue_cc_out_fv,
		                                                       h_ele_pfp_momentum_no_track_numu_cc, h_ele_pfp_momentum_no_track_numu_cc_mixed,
		                                                       h_ele_pfp_momentum_no_track_cosmic,  h_ele_pfp_momentum_no_track_nc,
		                                                       h_ele_pfp_momentum_no_track_nc_pi0,  h_ele_pfp_momentum_no_track_other_mixed,
		                                                       h_ele_pfp_momentum_no_track_unmatched, h_ele_pfp_momentum_no_track_intime, scaled_intime_scale_factor,
		                                                       h_ele_pfp_momentum_no_track_data, data_scale_factor,
		                                                       "", "Leading Shower Momentum [GeV] (Scaled)", "",
		                                                       Form("%s%s", file_locate_prefix, "post_cuts_leading_momentum_no_track_scaled_data.pdf"));

		histogram_functions::PlotSimpleStackDataMomentumRebin (h_ele_pfp_momentum_has_track_nue_cc,  h_ele_pfp_momentum_has_track_nue_cc_mixed,
		                                                       h_ele_pfp_momentum_has_track_nue_cc_out_fv,
		                                                       h_ele_pfp_momentum_has_track_numu_cc, h_ele_pfp_momentum_has_track_numu_cc_mixed,
		                                                       h_ele_pfp_momentum_has_track_cosmic,  h_ele_pfp_momentum_has_track_nc,
		                                                       h_ele_pfp_momentum_has_track_nc_pi0,  h_ele_pfp_momentum_has_track_other_mixed,
		                                                       h_ele_pfp_momentum_has_track_unmatched, h_ele_pfp_momentum_has_track_intime, scaled_intime_scale_factor,
		                                                       h_ele_pfp_momentum_has_track_data, data_scale_factor,
		                                                       "", "Leading Shower Momentum [GeV] (Scaled)", "",
		                                                       Form("%s%s", file_locate_prefix, "post_cuts_leading_momentum_has_track_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStackData (h_ele_pfp_phi_no_track_nue_cc,  h_ele_pfp_phi_no_track_nue_cc_mixed,
	                                          h_ele_pfp_phi_no_track_nue_cc_out_fv,
	                                          h_ele_pfp_phi_no_track_numu_cc, h_ele_pfp_phi_no_track_numu_cc_mixed,
	                                          h_ele_pfp_phi_no_track_cosmic,  h_ele_pfp_phi_no_track_nc,
	                                          h_ele_pfp_phi_no_track_nc_pi0,  h_ele_pfp_phi_no_track_other_mixed,
	                                          h_ele_pfp_phi_no_track_unmatched, h_ele_pfp_phi_no_track_intime, intime_scale_factor,
	                                          h_ele_pfp_phi_no_track_data, data_scale_factor,
	                                          0.73, 0.98, 1.0, 0.60, false, 1.5,
	                                          "", "Leading Shower Phi [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_phi_no_track_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_ele_pfp_phi_has_track_nue_cc,  h_ele_pfp_phi_has_track_nue_cc_mixed,
	                                          h_ele_pfp_phi_has_track_nue_cc_out_fv,
	                                          h_ele_pfp_phi_has_track_numu_cc, h_ele_pfp_phi_has_track_numu_cc_mixed,
	                                          h_ele_pfp_phi_has_track_cosmic,  h_ele_pfp_phi_has_track_nc,
	                                          h_ele_pfp_phi_has_track_nc_pi0,  h_ele_pfp_phi_has_track_other_mixed,
	                                          h_ele_pfp_phi_has_track_unmatched, h_ele_pfp_phi_has_track_intime, intime_scale_factor,
	                                          h_ele_pfp_phi_has_track_data, data_scale_factor,
	                                          0.73, 0.98, 1.0, 0.60, false, 1.5,
	                                          "", "Leading Shower Phi [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_phi_has_track_data.pdf"));
//
	histogram_functions::PlotSimpleStackData (h_ele_pfp_theta_no_track_nue_cc,  h_ele_pfp_theta_no_track_nue_cc_mixed,
	                                          h_ele_pfp_theta_no_track_nue_cc_out_fv,
	                                          h_ele_pfp_theta_no_track_numu_cc, h_ele_pfp_theta_no_track_numu_cc_mixed,
	                                          h_ele_pfp_theta_no_track_cosmic,  h_ele_pfp_theta_no_track_nc,
	                                          h_ele_pfp_theta_no_track_nc_pi0,  h_ele_pfp_theta_no_track_other_mixed,
	                                          h_ele_pfp_theta_no_track_unmatched, h_ele_pfp_theta_no_track_intime, intime_scale_factor,
	                                          h_ele_pfp_theta_no_track_data, data_scale_factor,
	                                          "", "Leading Shower Theta [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_theta_no_track_data.pdf"));

	histogram_functions::PlotSimpleStackData (h_ele_pfp_theta_has_track_nue_cc,  h_ele_pfp_theta_has_track_nue_cc_mixed,
	                                          h_ele_pfp_theta_has_track_nue_cc_out_fv,
	                                          h_ele_pfp_theta_has_track_numu_cc, h_ele_pfp_theta_has_track_numu_cc_mixed,
	                                          h_ele_pfp_theta_has_track_cosmic,  h_ele_pfp_theta_has_track_nc,
	                                          h_ele_pfp_theta_has_track_nc_pi0,  h_ele_pfp_theta_has_track_other_mixed,
	                                          h_ele_pfp_theta_has_track_unmatched, h_ele_pfp_theta_has_track_intime, intime_scale_factor,
	                                          h_ele_pfp_theta_has_track_data, data_scale_factor,
	                                          "", "Leading Shower Theta [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_theta_has_track_data.pdf"));


	histogram_functions::PlotSimpleStack (h_ele_pfp_theta_nue_cc,  h_ele_pfp_theta_nue_cc_mixed,
	                                      h_ele_pfp_theta_nue_cc_out_fv,
	                                      h_ele_pfp_theta_numu_cc, h_ele_pfp_theta_numu_cc_mixed,
	                                      h_ele_pfp_theta_cosmic,  h_ele_pfp_theta_nc,
	                                      h_ele_pfp_theta_nc_pi0,  h_ele_pfp_theta_other_mixed,
	                                      h_ele_pfp_theta_unmatched, "",
	                                      "Leading Shower Theta [Degrees]", "",
	                                      Form("%s%s", file_locate_prefix, "pre_collection_cut_leading_theta.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_theta_nue_cc,  h_ele_pfp_theta_nue_cc_mixed,
	                                            h_ele_pfp_theta_nue_cc_out_fv,
	                                            h_ele_pfp_theta_numu_cc, h_ele_pfp_theta_numu_cc_mixed,
	                                            h_ele_pfp_theta_cosmic,  h_ele_pfp_theta_nc,
	                                            h_ele_pfp_theta_nc_pi0,  h_ele_pfp_theta_other_mixed,
	                                            h_ele_pfp_theta_unmatched, h_ele_pfp_theta_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Theta [Degrees]", "",
	                                            Form("%s%s", file_locate_prefix, "pre_collection_cut_leading_theta_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_ele_pfp_theta_nue_cc,  h_ele_pfp_theta_nue_cc_mixed,
	                                          h_ele_pfp_theta_nue_cc_out_fv,
	                                          h_ele_pfp_theta_numu_cc, h_ele_pfp_theta_numu_cc_mixed,
	                                          h_ele_pfp_theta_cosmic,  h_ele_pfp_theta_nc,
	                                          h_ele_pfp_theta_nc_pi0,  h_ele_pfp_theta_other_mixed,
	                                          h_ele_pfp_theta_unmatched, h_ele_pfp_theta_intime, intime_scale_factor,
	                                          h_ele_pfp_theta_data, data_scale_factor,
	                                          "", "Leading Shower Theta [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "pre_collection_cut_leading_theta_data.pdf"));

	histogram_functions::PlotSimpleStack (h_ele_pfp_theta_after_nue_cc,  h_ele_pfp_theta_after_nue_cc_mixed,
	                                      h_ele_pfp_theta_after_nue_cc_out_fv,
	                                      h_ele_pfp_theta_after_numu_cc, h_ele_pfp_theta_after_numu_cc_mixed,
	                                      h_ele_pfp_theta_after_cosmic,  h_ele_pfp_theta_after_nc,
	                                      h_ele_pfp_theta_after_nc_pi0,  h_ele_pfp_theta_after_other_mixed,
	                                      h_ele_pfp_theta_after_unmatched, "",
	                                      "Leading Shower Theta [Degrees]", "",
	                                      Form("%s%s", file_locate_prefix, "post_collection_cut_leading_theta.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_theta_after_nue_cc,  h_ele_pfp_theta_after_nue_cc_mixed,
	                                            h_ele_pfp_theta_after_nue_cc_out_fv,
	                                            h_ele_pfp_theta_after_numu_cc, h_ele_pfp_theta_after_numu_cc_mixed,
	                                            h_ele_pfp_theta_after_cosmic,  h_ele_pfp_theta_after_nc,
	                                            h_ele_pfp_theta_after_nc_pi0,  h_ele_pfp_theta_after_other_mixed,
	                                            h_ele_pfp_theta_after_unmatched, h_ele_pfp_theta_after_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Theta [Degrees]", "",
	                                            Form("%s%s", file_locate_prefix, "post_collection_cut_leading_theta_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_ele_pfp_theta_after_nue_cc,  h_ele_pfp_theta_after_nue_cc_mixed,
	                                          h_ele_pfp_theta_after_nue_cc_out_fv,
	                                          h_ele_pfp_theta_after_numu_cc, h_ele_pfp_theta_after_numu_cc_mixed,
	                                          h_ele_pfp_theta_after_cosmic,  h_ele_pfp_theta_after_nc,
	                                          h_ele_pfp_theta_after_nc_pi0,  h_ele_pfp_theta_after_other_mixed,
	                                          h_ele_pfp_theta_after_unmatched, h_ele_pfp_theta_after_intime, intime_scale_factor,
	                                          h_ele_pfp_theta_after_data, data_scale_factor,
	                                          "", "Leading Shower Theta [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_collection_cut_leading_theta_data.pdf"));

	histogram_functions::PlotSimpleStack (h_ele_pfp_theta_last_nue_cc,  h_ele_pfp_theta_last_nue_cc_mixed,
	                                      h_ele_pfp_theta_last_nue_cc_out_fv,
	                                      h_ele_pfp_theta_last_numu_cc, h_ele_pfp_theta_last_numu_cc_mixed,
	                                      h_ele_pfp_theta_last_cosmic,  h_ele_pfp_theta_last_nc,
	                                      h_ele_pfp_theta_last_nc_pi0,  h_ele_pfp_theta_last_other_mixed,
	                                      h_ele_pfp_theta_last_unmatched, "",
	                                      "Leading Shower Theta [Degrees]", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_theta.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_theta_last_nue_cc,  h_ele_pfp_theta_last_nue_cc_mixed,
	                                            h_ele_pfp_theta_last_nue_cc_out_fv,
	                                            h_ele_pfp_theta_last_numu_cc, h_ele_pfp_theta_last_numu_cc_mixed,
	                                            h_ele_pfp_theta_last_cosmic,  h_ele_pfp_theta_last_nc,
	                                            h_ele_pfp_theta_last_nc_pi0,  h_ele_pfp_theta_last_other_mixed,
	                                            h_ele_pfp_theta_last_unmatched, h_ele_pfp_theta_last_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Theta [Degrees]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_theta_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_ele_pfp_theta_last_nue_cc,  h_ele_pfp_theta_last_nue_cc_mixed,
	                                          h_ele_pfp_theta_last_nue_cc_out_fv,
	                                          h_ele_pfp_theta_last_numu_cc, h_ele_pfp_theta_last_numu_cc_mixed,
	                                          h_ele_pfp_theta_last_cosmic,  h_ele_pfp_theta_last_nc,
	                                          h_ele_pfp_theta_last_nc_pi0,  h_ele_pfp_theta_last_other_mixed,
	                                          h_ele_pfp_theta_last_unmatched, h_ele_pfp_theta_last_intime, intime_scale_factor,
	                                          h_ele_pfp_theta_last_data, data_scale_factor,
	                                          "", "Leading Shower Theta [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_theta_last_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_ele_pfp_theta_last_nue_cc,  h_ele_pfp_theta_last_nue_cc_mixed,
		                                          h_ele_pfp_theta_last_nue_cc_out_fv,
		                                          h_ele_pfp_theta_last_numu_cc, h_ele_pfp_theta_last_numu_cc_mixed,
		                                          h_ele_pfp_theta_last_cosmic,  h_ele_pfp_theta_last_nc,
		                                          h_ele_pfp_theta_last_nc_pi0,  h_ele_pfp_theta_last_other_mixed,
		                                          h_ele_pfp_theta_last_unmatched, h_ele_pfp_theta_last_intime, scaled_intime_scale_factor,
		                                          h_ele_pfp_theta_last_data, data_scale_factor,
		                                          "", "Leading Shower Theta [Degrees] (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_theta_last_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStack (h_ele_pfp_phi_nue_cc,  h_ele_pfp_phi_nue_cc_mixed,
	                                      h_ele_pfp_phi_nue_cc_out_fv,
	                                      h_ele_pfp_phi_numu_cc, h_ele_pfp_phi_numu_cc_mixed,
	                                      h_ele_pfp_phi_cosmic,  h_ele_pfp_phi_nc,
	                                      h_ele_pfp_phi_nc_pi0,  h_ele_pfp_phi_other_mixed,
	                                      h_ele_pfp_phi_unmatched, "",
	                                      "Leading Shower Phi [Degrees]", "",
	                                      Form("%s%s", file_locate_prefix, "pre_collection_cut_leading_phi.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_phi_nue_cc,  h_ele_pfp_phi_nue_cc_mixed,
	                                            h_ele_pfp_phi_nue_cc_out_fv,
	                                            h_ele_pfp_phi_numu_cc, h_ele_pfp_phi_numu_cc_mixed,
	                                            h_ele_pfp_phi_cosmic,  h_ele_pfp_phi_nc,
	                                            h_ele_pfp_phi_nc_pi0,  h_ele_pfp_phi_other_mixed,
	                                            h_ele_pfp_phi_unmatched, h_ele_pfp_phi_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Phi [Degrees]", "",
	                                            Form("%s%s", file_locate_prefix, "pre_collection_cut_leading_phi_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_ele_pfp_phi_nue_cc,  h_ele_pfp_phi_nue_cc_mixed,
	                                          h_ele_pfp_phi_nue_cc_out_fv,
	                                          h_ele_pfp_phi_numu_cc, h_ele_pfp_phi_numu_cc_mixed,
	                                          h_ele_pfp_phi_cosmic,  h_ele_pfp_phi_nc,
	                                          h_ele_pfp_phi_nc_pi0,  h_ele_pfp_phi_other_mixed,
	                                          h_ele_pfp_phi_unmatched, h_ele_pfp_phi_intime, intime_scale_factor,
	                                          h_ele_pfp_phi_data, data_scale_factor,
	                                          0.73, 0.98, 1.0, 0.60, false, 1.5,
	                                          "", "Leading Shower Phi [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "pre_collection_cut_leading_phi_data.pdf"));

	histogram_functions::PlotSimpleStack (h_ele_pfp_phi_after_nue_cc,  h_ele_pfp_phi_after_nue_cc_mixed,
	                                      h_ele_pfp_phi_after_nue_cc_out_fv,
	                                      h_ele_pfp_phi_after_numu_cc, h_ele_pfp_phi_after_numu_cc_mixed,
	                                      h_ele_pfp_phi_after_cosmic,  h_ele_pfp_phi_after_nc,
	                                      h_ele_pfp_phi_after_nc_pi0,  h_ele_pfp_phi_after_other_mixed,
	                                      h_ele_pfp_phi_after_unmatched, "",
	                                      "Leading Shower Phi [Degrees]", "",
	                                      Form("%s%s", file_locate_prefix, "post_collection_cut_leading_phi.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_phi_after_nue_cc,  h_ele_pfp_phi_after_nue_cc_mixed,
	                                            h_ele_pfp_phi_after_nue_cc_out_fv,
	                                            h_ele_pfp_phi_after_numu_cc, h_ele_pfp_phi_after_numu_cc_mixed,
	                                            h_ele_pfp_phi_after_cosmic,  h_ele_pfp_phi_after_nc,
	                                            h_ele_pfp_phi_after_nc_pi0,  h_ele_pfp_phi_after_other_mixed,
	                                            h_ele_pfp_phi_after_unmatched, h_ele_pfp_phi_after_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Phi [Degrees]", "",
	                                            Form("%s%s", file_locate_prefix, "post_collection_cut_leading_phi_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_ele_pfp_phi_after_nue_cc,  h_ele_pfp_phi_after_nue_cc_mixed,
	                                          h_ele_pfp_phi_after_nue_cc_out_fv,
	                                          h_ele_pfp_phi_after_numu_cc, h_ele_pfp_phi_after_numu_cc_mixed,
	                                          h_ele_pfp_phi_after_cosmic,  h_ele_pfp_phi_after_nc,
	                                          h_ele_pfp_phi_after_nc_pi0,  h_ele_pfp_phi_after_other_mixed,
	                                          h_ele_pfp_phi_after_unmatched, h_ele_pfp_phi_after_intime, intime_scale_factor,
	                                          h_ele_pfp_phi_after_data, data_scale_factor,
	                                          0.73, 0.98, 1.0, 0.60, false, 1.5,
	                                          "", "Leading Shower Phi [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_collection_cut_leading_phi_data.pdf"));

	histogram_functions::PlotSimpleStack (h_ele_pfp_phi_last_nue_cc,  h_ele_pfp_phi_last_nue_cc_mixed,
	                                      h_ele_pfp_phi_last_nue_cc_out_fv,
	                                      h_ele_pfp_phi_last_numu_cc, h_ele_pfp_phi_last_numu_cc_mixed,
	                                      h_ele_pfp_phi_last_cosmic,  h_ele_pfp_phi_last_nc,
	                                      h_ele_pfp_phi_last_nc_pi0,  h_ele_pfp_phi_last_other_mixed,
	                                      h_ele_pfp_phi_last_unmatched, "",
	                                      "Leading Shower Phi [Degrees]", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_phi.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_phi_last_nue_cc,  h_ele_pfp_phi_last_nue_cc_mixed,
	                                            h_ele_pfp_phi_last_nue_cc_out_fv,
	                                            h_ele_pfp_phi_last_numu_cc, h_ele_pfp_phi_last_numu_cc_mixed,
	                                            h_ele_pfp_phi_last_cosmic,  h_ele_pfp_phi_last_nc,
	                                            h_ele_pfp_phi_last_nc_pi0,  h_ele_pfp_phi_last_other_mixed,
	                                            h_ele_pfp_phi_last_unmatched, h_ele_pfp_phi_last_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Phi [Degrees]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_phi_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_ele_pfp_phi_last_nue_cc,  h_ele_pfp_phi_last_nue_cc_mixed,
	                                          h_ele_pfp_phi_last_nue_cc_out_fv,
	                                          h_ele_pfp_phi_last_numu_cc, h_ele_pfp_phi_last_numu_cc_mixed,
	                                          h_ele_pfp_phi_last_cosmic,  h_ele_pfp_phi_last_nc,
	                                          h_ele_pfp_phi_last_nc_pi0,  h_ele_pfp_phi_last_other_mixed,
	                                          h_ele_pfp_phi_last_unmatched, h_ele_pfp_phi_last_intime, intime_scale_factor,
	                                          h_ele_pfp_phi_last_data, data_scale_factor,
	                                          0.73, 0.98, 1.0, 0.60, false, 1.5,
	                                          "", "Leading Shower Phi [Degrees]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_phi_last_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_ele_pfp_phi_last_nue_cc,  h_ele_pfp_phi_last_nue_cc_mixed,
		                                          h_ele_pfp_phi_last_nue_cc_out_fv,
		                                          h_ele_pfp_phi_last_numu_cc, h_ele_pfp_phi_last_numu_cc_mixed,
		                                          h_ele_pfp_phi_last_cosmic,  h_ele_pfp_phi_last_nc,
		                                          h_ele_pfp_phi_last_nc_pi0,  h_ele_pfp_phi_last_other_mixed,
		                                          h_ele_pfp_phi_last_unmatched, h_ele_pfp_phi_last_intime, scaled_intime_scale_factor,
		                                          h_ele_pfp_phi_last_data, data_scale_factor,
		                                          0.73, 0.98, 1.0, 0.60, false, 1.5,
		                                          "", "Leading Shower Phi [Degrees] (Scaled)", "",
		                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_phi_last_scaled_data.pdf"));
	}

	histogram_functions::PlotSimpleStack (h_leading_shwr_length_1shwr_nue_cc,  h_leading_shwr_length_1shwr_nue_cc_mixed,
	                                      h_leading_shwr_length_1shwr_nue_cc_out_fv,
	                                      h_leading_shwr_length_1shwr_numu_cc, h_leading_shwr_length_1shwr_numu_cc_mixed,
	                                      h_leading_shwr_length_1shwr_cosmic,  h_leading_shwr_length_1shwr_nc,
	                                      h_leading_shwr_length_1shwr_nc_pi0,  h_leading_shwr_length_1shwr_other_mixed,
	                                      h_leading_shwr_length_1shwr_unmatched, "",
	                                      "Leading Shower Length (1 Shower Events) [cm]", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_length_1shwr.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_length_1shwr_nue_cc,  h_leading_shwr_length_1shwr_nue_cc_mixed,
	                                            h_leading_shwr_length_1shwr_nue_cc_out_fv,
	                                            h_leading_shwr_length_1shwr_numu_cc, h_leading_shwr_length_1shwr_numu_cc_mixed,
	                                            h_leading_shwr_length_1shwr_cosmic,  h_leading_shwr_length_1shwr_nc,
	                                            h_leading_shwr_length_1shwr_nc_pi0,  h_leading_shwr_length_1shwr_other_mixed,
	                                            h_leading_shwr_length_1shwr_unmatched, h_leading_shwr_length_1shwr_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Length (1 Shower Events) [cm]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_length_1shwr_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_leading_shwr_length_1shwr_nue_cc,  h_leading_shwr_length_1shwr_nue_cc_mixed,
	                                          h_leading_shwr_length_1shwr_nue_cc_out_fv,
	                                          h_leading_shwr_length_1shwr_numu_cc, h_leading_shwr_length_1shwr_numu_cc_mixed,
	                                          h_leading_shwr_length_1shwr_cosmic,  h_leading_shwr_length_1shwr_nc,
	                                          h_leading_shwr_length_1shwr_nc_pi0,  h_leading_shwr_length_1shwr_other_mixed,
	                                          h_leading_shwr_length_1shwr_unmatched, h_leading_shwr_length_1shwr_intime, intime_scale_factor,
	                                          h_leading_shwr_length_1shwr_data, data_scale_factor,
	                                          "", "Leading Shower Length (1 Shower Events) [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_length_1shwr_data.pdf"));

	histogram_functions::PlotSimpleStack (h_leading_shwr_length_2shwr_nue_cc,  h_leading_shwr_length_2shwr_nue_cc_mixed,
	                                      h_leading_shwr_length_2shwr_nue_cc_out_fv,
	                                      h_leading_shwr_length_2shwr_numu_cc, h_leading_shwr_length_2shwr_numu_cc_mixed,
	                                      h_leading_shwr_length_2shwr_cosmic,  h_leading_shwr_length_2shwr_nc,
	                                      h_leading_shwr_length_2shwr_nc_pi0,  h_leading_shwr_length_2shwr_other_mixed,
	                                      h_leading_shwr_length_2shwr_unmatched, "",
	                                      "Leading Shower Length (2+ Shower Events) [cm]", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_length_2shwr.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_length_2shwr_nue_cc,  h_leading_shwr_length_2shwr_nue_cc_mixed,
	                                            h_leading_shwr_length_2shwr_nue_cc_out_fv,
	                                            h_leading_shwr_length_2shwr_numu_cc, h_leading_shwr_length_2shwr_numu_cc_mixed,
	                                            h_leading_shwr_length_2shwr_cosmic,  h_leading_shwr_length_2shwr_nc,
	                                            h_leading_shwr_length_2shwr_nc_pi0,  h_leading_shwr_length_2shwr_other_mixed,
	                                            h_leading_shwr_length_2shwr_unmatched, h_leading_shwr_length_2shwr_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Length (2+ Shower Events) [cm]", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_length_2shwr_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_leading_shwr_length_2shwr_nue_cc,  h_leading_shwr_length_2shwr_nue_cc_mixed,
	                                          h_leading_shwr_length_2shwr_nue_cc_out_fv,
	                                          h_leading_shwr_length_2shwr_numu_cc, h_leading_shwr_length_2shwr_numu_cc_mixed,
	                                          h_leading_shwr_length_2shwr_cosmic,  h_leading_shwr_length_2shwr_nc,
	                                          h_leading_shwr_length_2shwr_nc_pi0,  h_leading_shwr_length_2shwr_other_mixed,
	                                          h_leading_shwr_length_2shwr_unmatched, h_leading_shwr_length_2shwr_intime, intime_scale_factor,
	                                          h_leading_shwr_length_2shwr_data, data_scale_factor,
	                                          "", "Leading Shower Length (2+ Shower Events) [cm]", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_length_2shwr_data.pdf"));

	histogram_functions::PlotSimpleStack (h_leading_shwr_hits_1shwr_nue_cc,  h_leading_shwr_hits_1shwr_nue_cc_mixed,
	                                      h_leading_shwr_hits_1shwr_nue_cc_out_fv,
	                                      h_leading_shwr_hits_1shwr_numu_cc, h_leading_shwr_hits_1shwr_numu_cc_mixed,
	                                      h_leading_shwr_hits_1shwr_cosmic,  h_leading_shwr_hits_1shwr_nc,
	                                      h_leading_shwr_hits_1shwr_nc_pi0,  h_leading_shwr_hits_1shwr_other_mixed,
	                                      h_leading_shwr_hits_1shwr_unmatched, "",
	                                      "Leading Shower Hits (1 Shower Events)", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_hits_1shwr.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_hits_1shwr_nue_cc,  h_leading_shwr_hits_1shwr_nue_cc_mixed,
	                                            h_leading_shwr_hits_1shwr_nue_cc_out_fv,
	                                            h_leading_shwr_hits_1shwr_numu_cc, h_leading_shwr_hits_1shwr_numu_cc_mixed,
	                                            h_leading_shwr_hits_1shwr_cosmic,  h_leading_shwr_hits_1shwr_nc,
	                                            h_leading_shwr_hits_1shwr_nc_pi0,  h_leading_shwr_hits_1shwr_other_mixed,
	                                            h_leading_shwr_hits_1shwr_unmatched, h_leading_shwr_hits_1shwr_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Hits (1 Shower Events)", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_hits_1shwr_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_leading_shwr_hits_1shwr_nue_cc,  h_leading_shwr_hits_1shwr_nue_cc_mixed,
	                                          h_leading_shwr_hits_1shwr_nue_cc_out_fv,
	                                          h_leading_shwr_hits_1shwr_numu_cc, h_leading_shwr_hits_1shwr_numu_cc_mixed,
	                                          h_leading_shwr_hits_1shwr_cosmic,  h_leading_shwr_hits_1shwr_nc,
	                                          h_leading_shwr_hits_1shwr_nc_pi0,  h_leading_shwr_hits_1shwr_other_mixed,
	                                          h_leading_shwr_hits_1shwr_unmatched, h_leading_shwr_hits_1shwr_intime, intime_scale_factor,
	                                          h_leading_shwr_hits_1shwr_data, data_scale_factor,
	                                          "", "Leading Shower Hits (1 Shower Events)", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_hits_1shwr_data.pdf"));

	histogram_functions::PlotSimpleStack (h_leading_shwr_hits_2shwr_nue_cc,  h_leading_shwr_hits_2shwr_nue_cc_mixed,
	                                      h_leading_shwr_hits_2shwr_nue_cc_out_fv,
	                                      h_leading_shwr_hits_2shwr_numu_cc, h_leading_shwr_hits_2shwr_numu_cc_mixed,
	                                      h_leading_shwr_hits_2shwr_cosmic,  h_leading_shwr_hits_2shwr_nc,
	                                      h_leading_shwr_hits_2shwr_nc_pi0,  h_leading_shwr_hits_2shwr_other_mixed,
	                                      h_leading_shwr_hits_2shwr_unmatched, "",
	                                      "Leading Shower Hits (2+ Shower Events)", "",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_leading_hits_2shwr.pdf"));
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_hits_2shwr_nue_cc,  h_leading_shwr_hits_2shwr_nue_cc_mixed,
	                                            h_leading_shwr_hits_2shwr_nue_cc_out_fv,
	                                            h_leading_shwr_hits_2shwr_numu_cc, h_leading_shwr_hits_2shwr_numu_cc_mixed,
	                                            h_leading_shwr_hits_2shwr_cosmic,  h_leading_shwr_hits_2shwr_nc,
	                                            h_leading_shwr_hits_2shwr_nc_pi0,  h_leading_shwr_hits_2shwr_other_mixed,
	                                            h_leading_shwr_hits_2shwr_unmatched, h_leading_shwr_hits_2shwr_intime,
	                                            intime_scale_factor, data_scale_factor, "",
	                                            "Leading Shower Hits (2+ Shower Events)", "",
	                                            Form("%s%s", file_locate_prefix, "post_cuts_leading_hits_2shwr_intime.pdf"));
	histogram_functions::PlotSimpleStackData (h_leading_shwr_hits_2shwr_nue_cc,  h_leading_shwr_hits_2shwr_nue_cc_mixed,
	                                          h_leading_shwr_hits_2shwr_nue_cc_out_fv,
	                                          h_leading_shwr_hits_2shwr_numu_cc, h_leading_shwr_hits_2shwr_numu_cc_mixed,
	                                          h_leading_shwr_hits_2shwr_cosmic,  h_leading_shwr_hits_2shwr_nc,
	                                          h_leading_shwr_hits_2shwr_nc_pi0,  h_leading_shwr_hits_2shwr_other_mixed,
	                                          h_leading_shwr_hits_2shwr_unmatched, h_leading_shwr_hits_2shwr_intime, intime_scale_factor,
	                                          h_leading_shwr_hits_2shwr_data, data_scale_factor,
	                                          "", "Leading Shower Hits (2+ Shower Events)", "",
	                                          Form("%s%s", file_locate_prefix, "post_cuts_leading_hits_2shwr_data.pdf"));

	histogram_functions::PlotSimpleStackData(h_ele_pfp_x_nue_cc, h_ele_pfp_x_nue_cc_mixed,
	                                         h_ele_pfp_x_nue_cc_out_fv, h_ele_pfp_x_numu_cc, h_ele_pfp_x_numu_cc_mixed,
	                                         h_ele_pfp_x_cosmic, h_ele_pfp_x_nc, h_ele_pfp_x_nc_pi0, h_ele_pfp_x_other_mixed,
	                                         h_ele_pfp_x_unmatched, h_ele_pfp_x_intime, intime_scale_factor,
	                                         h_ele_pfp_x_data, data_scale_factor,
	                                         "", "Reco Vertex X [cm]", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_x_data.pdf"));
	histogram_functions::PlotSimpleStackData(h_ele_pfp_y_nue_cc, h_ele_pfp_y_nue_cc_mixed,
	                                         h_ele_pfp_y_nue_cc_out_fv, h_ele_pfp_y_numu_cc, h_ele_pfp_y_numu_cc_mixed,
	                                         h_ele_pfp_y_cosmic, h_ele_pfp_y_nc, h_ele_pfp_y_nc_pi0, h_ele_pfp_y_other_mixed,
	                                         h_ele_pfp_y_unmatched, h_ele_pfp_y_intime, intime_scale_factor,
	                                         h_ele_pfp_y_data, data_scale_factor,
	                                         "", "Reco Vertex Y [cm]", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_y_data.pdf"));
	histogram_functions::PlotSimpleStackData(h_ele_pfp_z_nue_cc, h_ele_pfp_z_nue_cc_mixed,
	                                         h_ele_pfp_z_nue_cc_out_fv, h_ele_pfp_z_numu_cc, h_ele_pfp_z_numu_cc_mixed,
	                                         h_ele_pfp_z_cosmic, h_ele_pfp_z_nc, h_ele_pfp_z_nc_pi0, h_ele_pfp_z_other_mixed,
	                                         h_ele_pfp_z_unmatched, h_ele_pfp_z_intime, intime_scale_factor,
	                                         h_ele_pfp_z_data, data_scale_factor,
	                                         "", "Reco Vertex Z [cm]", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_z_data.pdf"));

	histogram_functions::PlotSimpleStackData(h_any_pfp_x_nue_cc,        h_any_pfp_x_nue_cc_mixed,
	                                         h_any_pfp_x_nue_cc_out_fv, h_any_pfp_x_numu_cc,      h_any_pfp_x_numu_cc_mixed,
	                                         h_any_pfp_x_cosmic,        h_any_pfp_x_nc,           h_any_pfp_x_nc_pi0, h_any_pfp_x_other_mixed,
	                                         h_any_pfp_x_unmatched,     h_any_pfp_x_intime, intime_scale_factor,
	                                         h_any_pfp_x_data, data_scale_factor,
	                                         "", "Reco Vertex X [cm]", "",
	                                         Form("%s%s", file_locate_prefix, "pre_cuts_leading_pfp_x_data.pdf"));
	histogram_functions::PlotSimpleStackData(h_any_pfp_y_nue_cc,        h_any_pfp_y_nue_cc_mixed,
	                                         h_any_pfp_y_nue_cc_out_fv, h_any_pfp_y_numu_cc,      h_any_pfp_y_numu_cc_mixed,
	                                         h_any_pfp_y_cosmic,        h_any_pfp_y_nc,           h_any_pfp_y_nc_pi0, h_any_pfp_y_other_mixed,
	                                         h_any_pfp_y_unmatched,     h_any_pfp_y_intime, intime_scale_factor,
	                                         h_any_pfp_y_data, data_scale_factor,
	                                         "", "Reco Vertex Y [cm]", "",
	                                         Form("%s%s", file_locate_prefix, "pre_cuts_leading_pfp_y_data.pdf"));
	histogram_functions::PlotSimpleStackData(h_any_pfp_z_nue_cc,        h_any_pfp_z_nue_cc_mixed,
	                                         h_any_pfp_z_nue_cc_out_fv, h_any_pfp_z_numu_cc,      h_any_pfp_z_numu_cc_mixed,
	                                         h_any_pfp_z_cosmic,        h_any_pfp_z_nc,           h_any_pfp_z_nc_pi0, h_any_pfp_z_other_mixed,
	                                         h_any_pfp_z_unmatched,     h_any_pfp_z_intime, intime_scale_factor,
	                                         h_any_pfp_z_data, data_scale_factor,
	                                         "", "Reco Vertex Z [cm]", "",
	                                         Form("%s%s", file_locate_prefix, "pre_cuts_leading_pfp_z_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData(h_ele_pfp_x_nue_cc, h_ele_pfp_x_nue_cc_mixed,
		                                         h_ele_pfp_x_nue_cc_out_fv, h_ele_pfp_x_numu_cc, h_ele_pfp_x_numu_cc_mixed,
		                                         h_ele_pfp_x_cosmic, h_ele_pfp_x_nc, h_ele_pfp_x_nc_pi0, h_ele_pfp_x_other_mixed,
		                                         h_ele_pfp_x_unmatched, h_ele_pfp_x_intime, scaled_intime_scale_factor,
		                                         h_ele_pfp_x_data, data_scale_factor,
		                                         "", "Reco Vertex X [cm] (Scaled)", "",
		                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_x_scaled_data.pdf"));
		histogram_functions::PlotSimpleStackData(h_ele_pfp_y_nue_cc, h_ele_pfp_y_nue_cc_mixed,
		                                         h_ele_pfp_y_nue_cc_out_fv, h_ele_pfp_y_numu_cc, h_ele_pfp_y_numu_cc_mixed,
		                                         h_ele_pfp_y_cosmic, h_ele_pfp_y_nc, h_ele_pfp_y_nc_pi0, h_ele_pfp_y_other_mixed,
		                                         h_ele_pfp_y_unmatched, h_ele_pfp_y_intime, scaled_intime_scale_factor,
		                                         h_ele_pfp_y_data, data_scale_factor,
		                                         "", "Reco Vertex Y [cm] (Scaled)", "",
		                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_y_scaled_data.pdf"));
		histogram_functions::PlotSimpleStackData(h_ele_pfp_z_nue_cc, h_ele_pfp_z_nue_cc_mixed,
		                                         h_ele_pfp_z_nue_cc_out_fv, h_ele_pfp_z_numu_cc, h_ele_pfp_z_numu_cc_mixed,
		                                         h_ele_pfp_z_cosmic, h_ele_pfp_z_nc, h_ele_pfp_z_nc_pi0, h_ele_pfp_z_other_mixed,
		                                         h_ele_pfp_z_unmatched, h_ele_pfp_z_intime, scaled_intime_scale_factor,
		                                         h_ele_pfp_z_data, data_scale_factor,
		                                         "", "Reco Vertex Z [cm] (Scaled)", "",
		                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_z_scaled_data.pdf"));

		histogram_functions::PlotSimpleStackData(h_any_pfp_x_nue_cc,        h_any_pfp_x_nue_cc_mixed,
		                                         h_any_pfp_x_nue_cc_out_fv, h_any_pfp_x_numu_cc,      h_any_pfp_x_numu_cc_mixed,
		                                         h_any_pfp_x_cosmic,        h_any_pfp_x_nc,           h_any_pfp_x_nc_pi0, h_any_pfp_x_other_mixed,
		                                         h_any_pfp_x_unmatched,     h_any_pfp_x_intime, scaled_intime_scale_factor,
		                                         h_any_pfp_x_data, data_scale_factor,
		                                         "", "Reco Vertex X [cm] (Scaled)", "",
		                                         Form("%s%s", file_locate_prefix, "pre_cuts_leading_pfp_x_scaled_data.pdf"));
		histogram_functions::PlotSimpleStackData(h_any_pfp_y_nue_cc,        h_any_pfp_y_nue_cc_mixed,
		                                         h_any_pfp_y_nue_cc_out_fv, h_any_pfp_y_numu_cc,      h_any_pfp_y_numu_cc_mixed,
		                                         h_any_pfp_y_cosmic,        h_any_pfp_y_nc,           h_any_pfp_y_nc_pi0, h_any_pfp_y_other_mixed,
		                                         h_any_pfp_y_unmatched,     h_any_pfp_y_intime, scaled_intime_scale_factor,
		                                         h_any_pfp_y_data, data_scale_factor,
		                                         "", "Reco Vertex Y [cm] (Scaled)", "",
		                                         Form("%s%s", file_locate_prefix, "pre_cuts_leading_pfp_y_scaled_data.pdf"));
		histogram_functions::PlotSimpleStackData(h_any_pfp_z_nue_cc,        h_any_pfp_z_nue_cc_mixed,
		                                         h_any_pfp_z_nue_cc_out_fv, h_any_pfp_z_numu_cc,      h_any_pfp_z_numu_cc_mixed,
		                                         h_any_pfp_z_cosmic,        h_any_pfp_z_nc,           h_any_pfp_z_nc_pi0, h_any_pfp_z_other_mixed,
		                                         h_any_pfp_z_unmatched,     h_any_pfp_z_intime, scaled_intime_scale_factor,
		                                         h_any_pfp_z_data, data_scale_factor,
		                                         "", "Reco Vertex Z [cm] (Scaled)", "",
		                                         Form("%s%s", file_locate_prefix, "pre_cuts_leading_pfp_z_scaled_data.pdf"));
	}


	histogram_functions::PlotSimpleStackData(h_any_pfp_x_last_nue_cc,        h_any_pfp_x_last_nue_cc_mixed,
	                                         h_any_pfp_x_last_nue_cc_out_fv, h_any_pfp_x_last_numu_cc,      h_any_pfp_x_last_numu_cc_mixed,
	                                         h_any_pfp_x_last_cosmic,        h_any_pfp_x_last_nc,           h_any_pfp_x_last_nc_pi0, h_any_pfp_x_last_other_mixed,
	                                         h_any_pfp_x_last_unmatched,     h_any_pfp_x_last_intime, intime_scale_factor,
	                                         h_any_pfp_x_last_data, data_scale_factor,
	                                         "", "Reco Vertex X [cm]", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_x_last_data.pdf"));
	histogram_functions::PlotSimpleStackData(h_any_pfp_y_last_nue_cc,        h_any_pfp_y_last_nue_cc_mixed,
	                                         h_any_pfp_y_last_nue_cc_out_fv, h_any_pfp_y_last_numu_cc,      h_any_pfp_y_last_numu_cc_mixed,
	                                         h_any_pfp_y_last_cosmic,        h_any_pfp_y_last_nc,           h_any_pfp_y_last_nc_pi0, h_any_pfp_y_last_other_mixed,
	                                         h_any_pfp_y_last_unmatched,     h_any_pfp_y_last_intime, intime_scale_factor,
	                                         h_any_pfp_y_last_data, data_scale_factor,
	                                         "", "Reco Vertex Y [cm]", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_y_last_data.pdf"));
	histogram_functions::PlotSimpleStackData(h_any_pfp_z_last_nue_cc,        h_any_pfp_z_last_nue_cc_mixed,
	                                         h_any_pfp_z_last_nue_cc_out_fv, h_any_pfp_z_last_numu_cc,      h_any_pfp_z_last_numu_cc_mixed,
	                                         h_any_pfp_z_last_cosmic,        h_any_pfp_z_last_nc,           h_any_pfp_z_last_nc_pi0, h_any_pfp_z_last_other_mixed,
	                                         h_any_pfp_z_last_unmatched,     h_any_pfp_z_last_intime, intime_scale_factor,
	                                         h_any_pfp_z_last_data, data_scale_factor,
	                                         "", "Reco Vertex Z [cm]", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_z_last_data.pdf"));


	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_qe, "Post Cuts - Showers/Tracks Purity - QE",
	                                      "Reco Showers", "Reco Tracks",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_showers_tracks_purity_qe.pdf"), "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_res, "Post Cuts - Showers/Tracks Purity - Res",
	                                      "Reco Showers", "Reco Tracks",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_showers_tracks_purity_res.pdf"), "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_dis, "Post Cuts - Showers/Tracks Purity - DIS",
	                                      "Reco Showers", "Reco Tracks",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_showers_tracks_purity_dis.pdf"), "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_coh, "Post Cuts - Showers/Tracks Purity - Coh",
	                                      "Reco Showers", "Reco Tracks",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_showers_tracks_purity_coh.pdf"), "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_mec, "Post Cuts - Showers/Tracks Purity - MEC",
	                                      "Reco Showers", "Reco Tracks",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_showers_tracks_purity_mec.pdf"), "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_total, "Post Cuts - Showers/Tracks Purity - Total",
	                                      "Reco Showers", "Reco Tracks",
	                                      Form("%s%s", file_locate_prefix, "post_cuts_showers_tracks_purity_total.pdf"), "colz text", 3, 2);

	histogram_functions::Plot2DHistogram (h_ele_theta_phi_intime, "Post Cuts - Theta Phi - In Time",
	                                      "Phi [Degrees]", "Theta [Degrees]", Form("%s%s", file_locate_prefix, "post_cuts_theta_phi_intime.pdf"));

	histogram_functions::OverlayScatter(h_ele_theta_phi_nue_cc, h_ele_theta_phi_nue_cc_mixed, h_ele_theta_phi_nue_cc_out_fv, h_ele_theta_phi_numu_cc,
	                                    h_ele_theta_phi_numu_cc_mixed, h_ele_theta_phi_cosmic, h_ele_theta_phi_nc,
	                                    h_ele_theta_phi_nc_pi0, h_ele_theta_phi_other_mixed, h_ele_theta_phi_unmatched,
	                                    0.15, 0.35, 0.65, 0.90, "", "Phi [Degrees]", "Theta [Degrees]",
	                                    Form("%s%s", file_locate_prefix, "post_cuts_theta_phi_scatter.pdf"));

	histogram_functions::Plot2DHistogram(h_ele_eng_costheta_nue_cc, "", "Reco Electron Momentum [GeV]", "Reco Electron Cos(#theta)",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_eng_costheta_nue_cc.pdf"));
	histogram_functions::LegoStackData(h_ele_eng_costheta_nue_cc, h_ele_eng_costheta_nue_cc_mixed, h_ele_eng_costheta_nue_cc_out_fv,
	                                   h_ele_eng_costheta_numu_cc,  h_ele_eng_costheta_cosmic, h_ele_eng_costheta_nc,
	                                   h_ele_eng_costheta_nc_pi0, h_ele_eng_costheta_other_mixed,
	                                   h_ele_eng_costheta_unmatched, h_ele_eng_costheta_intime, intime_scale_factor,
	                                   h_ele_eng_costheta_data, data_scale_factor, 0.75, 0.95, 0.70, 0.95,
	                                   "", "Reco Electron Momentum [GeV]", "Reco Electron Cos(#theta)",
	                                   Form("%s%s", file_locate_prefix, "post_cuts_leading_pfp_eng_costheta_data.pdf"));

	histogram_functions::Plot2DHistogram(h_true_reco_ele_momentum, "Selected True Electrons", "True Electron Momentum [GeV]", "Reco Electron Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_ele_true_reco_momentum.pdf"));
	histogram_functions::Plot2DHistogram(h_true_reco_ele_costheta, "Selected True Electrons", "True Electron Cos(#theta) [GeV]", "Reco Electron Cos(#theta) [GeV]",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_ele_true_reco_costheta.pdf"));
	histogram_functions::Plot1DHistogram(h_true_num_e, "Number of Selected True Electrons Per Event",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_num_true_ele.pdf"));

	histogram_functions::Plot2DHistogram(h_true_reco_ele_momentum_pre, "Pre Selection True Electrons",
	                                     "True Electron Momentum [GeV]", "Reco Electron Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_ele_true_reco_momentum_pre.pdf"));
	histogram_functions::Plot2DHistogram(h_true_reco_ele_costheta_pre, "Pre Selection True Electrons",
	                                     "True Electron Cos(#theta) [GeV]", "Reco Electron Cos(#theta) [GeV]",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_ele_true_reco_costheta_pre.pdf"));
	histogram_functions::Plot1DHistogram(h_true_num_e_pre, "Number of Pre Selection True Electrons Per Event",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_num_true_ele_pre.pdf"));

	histogram_functions::PlotSimpleStackData(h_ele_eng_for_nue_cc,        h_ele_eng_for_nue_cc_mixed,
	                                         h_ele_eng_for_nue_cc_out_fv, h_ele_eng_for_numu_cc,   h_ele_eng_for_numu_cc_mixed,
	                                         h_ele_eng_for_cosmic,        h_ele_eng_for_nc,        h_ele_eng_for_nc_pi0,
	                                         h_ele_eng_for_other_mixed,   h_ele_eng_for_unmatched, h_ele_eng_for_intime,
	                                         intime_scale_factor,         h_ele_eng_for_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (Cos(#theta) >= 0.5)", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_eng_for_data.pdf"));

	histogram_functions::PlotSimpleStackData(h_ele_eng_mid_nue_cc,        h_ele_eng_mid_nue_cc_mixed,
	                                         h_ele_eng_mid_nue_cc_out_fv, h_ele_eng_mid_numu_cc,   h_ele_eng_mid_numu_cc_mixed,
	                                         h_ele_eng_mid_cosmic,        h_ele_eng_mid_nc,        h_ele_eng_mid_nc_pi0,
	                                         h_ele_eng_mid_other_mixed,   h_ele_eng_mid_unmatched, h_ele_eng_mid_intime,
	                                         intime_scale_factor,         h_ele_eng_mid_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (0.5 > Cos(#theta) > -0.5)", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_eng_mid_data.pdf"));

	histogram_functions::PlotSimpleStackData(h_ele_eng_back_nue_cc,        h_ele_eng_back_nue_cc_mixed,
	                                         h_ele_eng_back_nue_cc_out_fv, h_ele_eng_back_numu_cc,   h_ele_eng_back_numu_cc_mixed,
	                                         h_ele_eng_back_cosmic,        h_ele_eng_back_nc,        h_ele_eng_back_nc_pi0,
	                                         h_ele_eng_back_other_mixed,   h_ele_eng_back_unmatched, h_ele_eng_back_intime,
	                                         intime_scale_factor,          h_ele_eng_back_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (Cos(#theta) <= -0.5)", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_eng_back_data.pdf"));

	histogram_functions::PlotSimpleStackData(h_ele_eng_for_trans_nue_cc,        h_ele_eng_for_trans_nue_cc_mixed,
	                                         h_ele_eng_for_trans_nue_cc_out_fv, h_ele_eng_for_trans_numu_cc,   h_ele_eng_for_trans_numu_cc_mixed,
	                                         h_ele_eng_for_trans_cosmic,        h_ele_eng_for_trans_nc,        h_ele_eng_for_trans_nc_pi0,
	                                         h_ele_eng_for_trans_other_mixed,   h_ele_eng_for_trans_unmatched, h_ele_eng_for_trans_intime,
	                                         intime_scale_factor,               h_ele_eng_for_trans_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (Cos(#theta_t) >= 0.5)", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_eng_for_trans_data.pdf"));

	histogram_functions::PlotSimpleStackData(h_ele_eng_mid_trans_nue_cc,        h_ele_eng_mid_trans_nue_cc_mixed,
	                                         h_ele_eng_mid_trans_nue_cc_out_fv, h_ele_eng_mid_trans_numu_cc,   h_ele_eng_mid_trans_numu_cc_mixed,
	                                         h_ele_eng_mid_trans_cosmic,        h_ele_eng_mid_trans_nc,        h_ele_eng_mid_trans_nc_pi0,
	                                         h_ele_eng_mid_trans_other_mixed,   h_ele_eng_mid_trans_unmatched, h_ele_eng_mid_trans_intime,
	                                         intime_scale_factor,         h_ele_eng_mid_trans_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (0.5 > Cos(#theta_t) > -0.5)", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_eng_mid_trans_data.pdf"));

	histogram_functions::PlotSimpleStackData(h_ele_eng_back_trans_nue_cc,        h_ele_eng_back_trans_nue_cc_mixed,
	                                         h_ele_eng_back_trans_nue_cc_out_fv, h_ele_eng_back_trans_numu_cc,   h_ele_eng_back_trans_numu_cc_mixed,
	                                         h_ele_eng_back_trans_cosmic,        h_ele_eng_back_trans_nc,        h_ele_eng_back_trans_nc_pi0,
	                                         h_ele_eng_back_trans_other_mixed,   h_ele_eng_back_trans_unmatched, h_ele_eng_back_trans_intime,
	                                         intime_scale_factor,          h_ele_eng_back_trans_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (Cos(#theta_t) <= -0.5)", "",
	                                         Form("%s%s", file_locate_prefix, "post_cuts_leading_eng_back_trans_data.pdf"));

	histogram_functions::Plot2DHistogram(h_mc_vtx_xy_nue_cc,     "", "True Signal Nue Vtx X [cm]", "True Signal Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_vtx_xy_nue_cc.pdf"    ));
	histogram_functions::Plot2DHistogram(h_mc_vtx_xz_nue_cc,     "", "True Signal Nue Vtx X [cm]", "True Signal Nue Vtx Z [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_vtx_xz_nue_cc.pdf"    ));
	histogram_functions::Plot2DHistogram(h_mc_vtx_yz_nue_cc,     "", "True Signal Nue Vtx Z [cm]", "True Signal Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_vtx_yz_nue_cc.pdf"    ));
	histogram_functions::Plot2DHistogram(h_reco_vtx_xy_nue_cc,   "", "Reco Signal Nue Vtx X [cm]", "Reco Signal Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "reco_vtx_xy_nue_cc.pdf"    ));
	histogram_functions::Plot2DHistogram(h_reco_vtx_xz_nue_cc,   "", "Reco Signal Nue Vtx X [cm]", "Reco Signal Nue Vtx Z [cm]",
	                                     Form("%s%s", file_locate_prefix, "reco_vtx_xz_nue_cc.pdf"    ));
	histogram_functions::Plot2DHistogram(h_reco_vtx_yz_nue_cc,   "", "Reco Signal Nue Vtx Z [cm]", "Reco Signal Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "reco vtx_yz_nue_cc.pdf"    ));
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_x_nue_cc, "", "True Signal Nue Vtx X [cm]", "Reco Signal Nue Vtx X [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_reco_vtx_x_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_y_nue_cc, "", "True Signal Nue Vtx Y [cm]", "Reco Signal Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_reco_vtx_y_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_z_nue_cc, "", "True Signal Nue Vtx Z [cm]", "Reco Signal Nue Vtx Z [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_reco_vtx_z_nue_cc.pdf"));

	histogram_functions::Plot2DHistogram(h_mc_vtx_xy_nue_cc_out_fv, "", "True OutFV Nue Vtx X [cm]", "True OutFV Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_vtx_xy_nue_cc_out_fv.pdf"    ));
	histogram_functions::Plot2DHistogram(h_mc_vtx_xz_nue_cc_out_fv, "", "True OutFV Nue Vtx X [cm]", "True OutFV Nue Vtx Z [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_vtx_xz_nue_cc_out_fv.pdf"    ));
	histogram_functions::Plot2DHistogram(h_mc_vtx_yz_nue_cc_out_fv, "", "True OutFV Nue Vtx Z [cm]", "True OutFV Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_vtx_yz_nue_cc_out_fv.pdf"    ));
	histogram_functions::Plot2DHistogram(h_reco_vtx_xy_nue_cc_out_fv, "", "Reco OutFV Nue Vtx X [cm]", "Reco OutFV Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "reco_vtx_xy_nue_cc_out_fv.pdf"    ));
	histogram_functions::Plot2DHistogram(h_reco_vtx_xz_nue_cc_out_fv, "", "Reco OutFV Nue Vtx X [cm]", "Reco OutFV Nue Vtx Z [cm]",
	                                     Form("%s%s", file_locate_prefix, "reco_vtx_xz_nue_cc_out_fv.pdf"    ));
	histogram_functions::Plot2DHistogram(h_reco_vtx_yz_nue_cc_out_fv, "", "Reco OutFV Nue Vtx Z [cm]", "Reco OutFV Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "reco vtx_yz_nue_cc_out_fv.pdf"    ));
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_x_nue_cc_out_fv, "", "True OutFV Nue Vtx X [cm]", "Reco OutFV Nue Vtx X [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_reco_vtx_x_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_y_nue_cc_out_fv, "", "True OutFV Nue Vtx Y [cm]", "Reco OutFV Nue Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_reco_vtx_y_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_z_nue_cc_out_fv, "", "True OutFV Nue Vtx Z [cm]", "Reco OutFV Nue Vtx Z [cm]",
	                                     Form("%s%s", file_locate_prefix, "true_reco_vtx_z_nue_cc_out_fv.pdf"));

	histogram_functions::Plot2DHistogram(h_dedx_cuts_hits_electron,     "", "Electron - dE/dx [MeV/cm]",   "Total Hits",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_dedx_hits_electron.pdf"    ));
	histogram_functions::Plot2DHistogram(h_dedx_cuts_hits_proton,       "", "Proton - dE/dx [MeV/cm]",     "Total Hits",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_dedx_hits_proton.pdf"      ));
	histogram_functions::Plot2DHistogram(h_dedx_cuts_hits_photon,       "", "Photon - dE/dx [MeV/cm]",     "Total Hits",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_dedx_hits_photon.pdf"      ));
	histogram_functions::Plot2DHistogram(h_dedx_cuts_hits_pion,         "", "Pion - dE/dx [MeV/cm]",       "Total Hits",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_dedx_hits_pion.pdf"        ));
	histogram_functions::Plot2DHistogram(h_dedx_cuts_hits_kaon,         "", "Kaon - dE/dx [MeV/cm]",       "Total Hits",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_dedx_hits_kaon.pdf"        ));
	histogram_functions::Plot2DHistogram(h_dedx_cuts_hits_muon,         "", "Muon - dE/dx [MeV/cm]",       "Total Hits",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_dedx_hits_muon.pdf"        ));
	histogram_functions::Plot2DHistogram(h_dedx_cuts_hits_neutron,      "", "Neutron - dE/dx [MeV/cm]",    "Total Hits",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_dedx_hits_neutron.pdf"     ));
	histogram_functions::Plot2DHistogram(h_dedx_cuts_hits_mc_unmatched, "", "Unmatched - dE/dx [MeV/cm]",  "Total Hits",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_dedx_hits_mc_unmatched.pdf"));

	histogram_functions::Plot2DHistogram(h_dedx_cuts_collection_hits_electron,     "", "Electron - dE/dx [MeV/cm]",
	                                     "Collection Hits", Form("%s%s", file_locate_prefix, "post_cuts_dedx_collection_hits_electron.pdf")    );
	histogram_functions::Plot2DHistogram(h_dedx_cuts_collection_hits_proton,       "", "Proton - dE/dx [MeV/cm]",
	                                     "Collection Hits", Form("%s%s", file_locate_prefix, "post_cuts_dedx_collection_hits_proton.pdf")      );
	histogram_functions::Plot2DHistogram(h_dedx_cuts_collection_hits_photon,       "", "Photon - dE/dx [MeV/cm]",
	                                     "Collection Hits", Form("%s%s", file_locate_prefix, "post_cuts_dedx_collection_hits_photon.pdf")      );
	histogram_functions::Plot2DHistogram(h_dedx_cuts_collection_hits_pion,         "", "Pion - dE/dx [MeV/cm]",
	                                     "Collection Hits", Form("%s%s", file_locate_prefix, "post_cuts_dedx_collection_hits_pion.pdf")        );
	histogram_functions::Plot2DHistogram(h_dedx_cuts_collection_hits_kaon,         "", "Kaon - dE/dx [MeV/cm]",
	                                     "Collection Hits", Form("%s%s", file_locate_prefix, "post_cuts_dedx_collection_hits_kaon.pdf")        );
	histogram_functions::Plot2DHistogram(h_dedx_cuts_collection_hits_muon,         "", "Muon - dE/dx [MeV/cm]",
	                                     "Collection Hits", Form("%s%s", file_locate_prefix, "post_cuts_dedx_collection_hits_muon.pdf")        );
	histogram_functions::Plot2DHistogram(h_dedx_cuts_collection_hits_neutron,      "", "Neutron - dE/dx [MeV/cm]",
	                                     "Collection Hits", Form("%s%s", file_locate_prefix, "post_cuts_dedx_collection_hits_neutron.pdf")     );
	histogram_functions::Plot2DHistogram(h_dedx_cuts_collection_hits_mc_unmatched, "", "Unmatched - dE/dx [MeV/cm]",
	                                     "Collection Hits", Form("%s%s", file_locate_prefix, "post_cuts_dedx_collection_hits_mc_unmatched.pdf"));

	histogram_functions::PlotSimpleStackParticle(h_dedx_cuts_electron, h_dedx_cuts_proton, h_dedx_cuts_photon, h_dedx_cuts_pion,
	                                             h_dedx_cuts_kaon, h_dedx_cuts_muon, h_dedx_cuts_neutron, h_dedx_cuts_mc_unmatched,
	                                             h_dedx_cuts_ext_unmatched, intime_scale_factor, data_scale_factor, 0.73, 0.98, 0.98, 0.50,
	                                             "True Particle Type - Reco dE/dx", "Leading Shower dE/dx [MeV/cm]", "",
	                                             Form("%s%s", file_locate_prefix, "post_cuts_dedx_true_particle.pdf"));
	histogram_functions::PlotSimpleStackParticle(h_dedx_cuts_last_electron, h_dedx_cuts_last_proton, h_dedx_cuts_last_photon, h_dedx_cuts_last_pion,
	                                             h_dedx_cuts_last_kaon, h_dedx_cuts_last_muon, h_dedx_cuts_last_neutron, h_dedx_cuts_last_mc_unmatched,
	                                             h_dedx_cuts_last_ext_unmatched, intime_scale_factor, data_scale_factor, 0.73, 0.98, 0.98, 0.50,
	                                             "True Particle Type - Reco dE/dx", "Leading Shower dE/dx [MeV/cm]",
	                                             "", Form("%s%s", file_locate_prefix, "post_cuts_dedx_true_particle_last.pdf"));


	histogram_functions::PlotSimpleStackInTime (h_track_containment_nue_cc,
	                                            h_track_containment_nue_cc_out_fv,
	                                            h_track_containment_nue_cc_mixed,
	                                            h_track_containment_numu_cc,
	                                            h_track_containment_numu_cc_mixed,
	                                            h_track_containment_nc,
	                                            h_track_containment_nc_pi0,
	                                            h_track_containment_cosmic,
	                                            h_track_containment_other_mixed,
	                                            h_track_containment_unmatched,
	                                            h_track_containment_intime, intime_scale_factor, data_scale_factor, "",
	                                            "Track Containment", "", Form("%s%s", file_locate_prefix, "track_containment.pdf"));

	histogram_functions::PlotSimpleStackData (h_track_containment_nue_cc,
	                                          h_track_containment_nue_cc_out_fv,
	                                          h_track_containment_nue_cc_mixed,
	                                          h_track_containment_numu_cc,
	                                          h_track_containment_numu_cc_mixed,
	                                          h_track_containment_nc,
	                                          h_track_containment_nc_pi0,
	                                          h_track_containment_cosmic,
	                                          h_track_containment_other_mixed,
	                                          h_track_containment_unmatched,
	                                          h_track_containment_intime, intime_scale_factor,
	                                          h_track_containment_data, data_scale_factor,
	                                          0.15, 0.40, 0.98, 0.50, false, " ",
	                                          "Track Containment", "", Form("%s%s", file_locate_prefix, "track_containment_data.pdf"));
	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData (h_track_containment_nue_cc,
		                                          h_track_containment_nue_cc_out_fv,
		                                          h_track_containment_nue_cc_mixed,
		                                          h_track_containment_numu_cc,
		                                          h_track_containment_numu_cc_mixed,
		                                          h_track_containment_nc,
		                                          h_track_containment_nc_pi0,
		                                          h_track_containment_cosmic,
		                                          h_track_containment_other_mixed,
		                                          h_track_containment_unmatched,
		                                          h_track_containment_intime, scaled_intime_scale_factor,
		                                          h_track_containment_data, data_scale_factor, "",
		                                          "Track Containment (Scaled)", "", Form("%s%s", file_locate_prefix, "track_containment_scaled_data.pdf"));
	}

	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_nue_cc, "", "Nue CC - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_nue_cc_out_fv, "", "Nue CC Out FV - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_nue_cc_out_fv.pdf"));
	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_nue_cc_mixed, "", "Nue CC Mixed - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_nue_cc_mixed.pdf"));
	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_numu_cc, "", "Numu CC - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_numu_cc.pdf"));
	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_nc, "", "NC - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_nc.pdf"));
	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_nc_pi0, "", "NC Pi 0 - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_nc_pi0.pdf"));
	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_cosmic, "", "Cosmic - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_cosmic.pdf"));
	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_other_mixed, "", "NC Mixed - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_nc_mixed.pdf"));
	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_unmatched, "", "Unmatched - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_unmatched.pdf"));
	histogram_functions::Plot2DHistogram(h_dedx_collection_angle_intime, "", "EXT In-Time - dE/dx [MeV/cm]", "Shower Angle to Collection Plane",
	                                     Form("%s%s", file_locate_prefix, "dedx_collection_angle_intime.pdf"));


	histogram_functions::Plot1DHistogram(h_mc_ele_e_1,  "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_1.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_2,  "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_2.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_3,  "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_3.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_4,  "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_4.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_5,  "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_5.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_6,  "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_6.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_7,  "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_7.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_8,  "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_8.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_9,  "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_9.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_10, "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_10.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_11, "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_11.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_12, "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_12.pdf"));
	histogram_functions::Plot1DHistogram(h_mc_ele_e_13, "True Selected Electron Energy [GeV]", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_13.pdf"));

	histogram_functions::Plot1DHistogram(h_reco_ele_e_1,  "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_1.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_2,  "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_2.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_3,  "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_3.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_4,  "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_4.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_5,  "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_5.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_6,  "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_6.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_7,  "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_7.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_8,  "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_8.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_9,  "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_9.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_10, "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_10.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_11, "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_11.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_12, "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_12.pdf"));
	histogram_functions::Plot1DHistogram(h_reco_ele_e_13, "Reco Selected Leading Shower Momentum [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_13.pdf"));

	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_1, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_1.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_2, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_2.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_3, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_3.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_4, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_4.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_5, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_5.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_6, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_6.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_7, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_7.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_8, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_8.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_9, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_9.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_10, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_10.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_11, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_11.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_12, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_12.pdf"));
	histogram_functions::Plot2DHistogram(h_mc_reco_ele_e_13, "", "True Electron Energy [GeV]", "Reco Leading Shower Energy [GeV]",
	                                     Form("%s%s", file_locate_prefix, "selected_mc_reco_ele_eng_13.pdf"));


	histogram_functions::EnergyOverlay(
	        h_mc_ele_e_0,
	        h_mc_ele_e_1,
	        h_mc_ele_e_2,
	        h_mc_ele_e_3,
	        h_mc_ele_e_4,
	        h_mc_ele_e_5,
	        h_mc_ele_e_6,
	        h_mc_ele_e_7,
	        h_mc_ele_e_8,
	        h_mc_ele_e_9,
	        h_mc_ele_e_10,
	        h_mc_ele_e_11,
	        h_mc_ele_e_12,
	        h_mc_ele_e_13,
	        "True Selected Electron Energy [GeV]",
	        "", Form("%s%s", file_locate_prefix, "selected_mc_ele_eng_overlay.pdf"));

	histogram_functions::EnergyOverlay(
	        h_reco_ele_e_0,
	        h_reco_ele_e_1,
	        h_reco_ele_e_2,
	        h_reco_ele_e_3,
	        h_reco_ele_e_4,
	        h_reco_ele_e_5,
	        h_reco_ele_e_6,
	        h_reco_ele_e_7,
	        h_reco_ele_e_8,
	        h_reco_ele_e_9,
	        h_reco_ele_e_10,
	        h_reco_ele_e_11,
	        h_reco_ele_e_12,
	        h_reco_ele_e_13,
	        "Reco Selected Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "selected_reco_ele_eng_overlay.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(h_ele_momentum_no_cut_nue_cc,        h_ele_momentum_no_cut_nue_cc_mixed,
	                                                      h_ele_momentum_no_cut_nue_cc_out_fv, h_ele_momentum_no_cut_numu_cc,   h_ele_momentum_no_cut_numu_cc_mixed,
	                                                      h_ele_momentum_no_cut_cosmic,        h_ele_momentum_no_cut_nc,        h_ele_momentum_no_cut_nc_pi0,
	                                                      h_ele_momentum_no_cut_other_mixed,   h_ele_momentum_no_cut_unmatched, h_ele_momentum_no_cut_intime,
	                                                      intime_scale_factor,                 h_ele_momentum_no_cut_data,      data_scale_factor,
	                                                      "No Cut", "Leading Shower Momentum [GeV]", "",
	                                                      Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_no_cut_data.pdf"));


	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_nue_cut_nue_cc,
	        h_ele_momentum_nue_cut_nue_cc_mixed,
	        h_ele_momentum_nue_cut_nue_cc_out_fv,
	        h_ele_momentum_nue_cut_numu_cc,
	        h_ele_momentum_nue_cut_numu_cc_mixed,
	        h_ele_momentum_nue_cut_cosmic,
	        h_ele_momentum_nue_cut_nc,
	        h_ele_momentum_nue_cut_nc_pi0,
	        h_ele_momentum_nue_cut_other_mixed,
	        h_ele_momentum_nue_cut_unmatched,
	        h_ele_momentum_nue_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_nue_cut_data,
	        data_scale_factor,
	        "Nue Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_nue_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_nue_cut_electron,
	        h_leading_momentum_nue_cut_proton,
	        h_leading_momentum_nue_cut_photon,
	        h_leading_momentum_nue_cut_pion,
	        h_leading_momentum_nue_cut_kaon,
	        h_leading_momentum_nue_cut_muon,
	        h_leading_momentum_nue_cut_neutron,
	        h_leading_momentum_nue_cut_mc_unmatched,
	        h_leading_momentum_nue_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_nue_cut.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_fv_cut_nue_cc,
	        h_ele_momentum_fv_cut_nue_cc_mixed,
	        h_ele_momentum_fv_cut_nue_cc_out_fv,
	        h_ele_momentum_fv_cut_numu_cc,
	        h_ele_momentum_fv_cut_numu_cc_mixed,
	        h_ele_momentum_fv_cut_cosmic,
	        h_ele_momentum_fv_cut_nc,
	        h_ele_momentum_fv_cut_nc_pi0,
	        h_ele_momentum_fv_cut_other_mixed,
	        h_ele_momentum_fv_cut_unmatched,
	        h_ele_momentum_fv_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_fv_cut_data,
	        data_scale_factor,
	        "Fiducial Volume Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_fv_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_fv_cut_electron,
	        h_leading_momentum_fv_cut_proton,
	        h_leading_momentum_fv_cut_photon,
	        h_leading_momentum_fv_cut_pion,
	        h_leading_momentum_fv_cut_kaon,
	        h_leading_momentum_fv_cut_muon,
	        h_leading_momentum_fv_cut_neutron,
	        h_leading_momentum_fv_cut_mc_unmatched,
	        h_leading_momentum_fv_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_fv_cut.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_flash_vtx_cut_nue_cc,
	        h_ele_momentum_flash_vtx_cut_nue_cc_mixed,
	        h_ele_momentum_flash_vtx_cut_nue_cc_out_fv,
	        h_ele_momentum_flash_vtx_cut_numu_cc,
	        h_ele_momentum_flash_vtx_cut_numu_cc_mixed,
	        h_ele_momentum_flash_vtx_cut_cosmic,
	        h_ele_momentum_flash_vtx_cut_nc,
	        h_ele_momentum_flash_vtx_cut_nc_pi0,
	        h_ele_momentum_flash_vtx_cut_other_mixed,
	        h_ele_momentum_flash_vtx_cut_unmatched,
	        h_ele_momentum_flash_vtx_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_flash_vtx_cut_data,
	        data_scale_factor,
	        "Flash Vertex Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_flash_vtx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_flash_vtx_electron,
	        h_leading_momentum_flash_vtx_proton,
	        h_leading_momentum_flash_vtx_photon,
	        h_leading_momentum_flash_vtx_pion,
	        h_leading_momentum_flash_vtx_kaon,
	        h_leading_momentum_flash_vtx_muon,
	        h_leading_momentum_flash_vtx_neutron,
	        h_leading_momentum_flash_vtx_mc_unmatched,
	        h_leading_momentum_flash_vtx_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_flash_vtx.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_shwr_vtx_cut_nue_cc,
	        h_ele_momentum_shwr_vtx_cut_nue_cc_mixed,
	        h_ele_momentum_shwr_vtx_cut_nue_cc_out_fv,
	        h_ele_momentum_shwr_vtx_cut_numu_cc,
	        h_ele_momentum_shwr_vtx_cut_numu_cc_mixed,
	        h_ele_momentum_shwr_vtx_cut_cosmic,
	        h_ele_momentum_shwr_vtx_cut_nc,
	        h_ele_momentum_shwr_vtx_cut_nc_pi0,
	        h_ele_momentum_shwr_vtx_cut_other_mixed,
	        h_ele_momentum_shwr_vtx_cut_unmatched,
	        h_ele_momentum_shwr_vtx_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_shwr_vtx_cut_data,
	        data_scale_factor,
	        "Shower Vertex Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_shwr_vtx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_shwr_vtx_electron,
	        h_leading_momentum_shwr_vtx_proton,
	        h_leading_momentum_shwr_vtx_photon,
	        h_leading_momentum_shwr_vtx_pion,
	        h_leading_momentum_shwr_vtx_kaon,
	        h_leading_momentum_shwr_vtx_muon,
	        h_leading_momentum_shwr_vtx_neutron,
	        h_leading_momentum_shwr_vtx_mc_unmatched,
	        h_leading_momentum_shwr_vtx_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_shwr_vtx.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_trk_vtx_cut_nue_cc,
	        h_ele_momentum_trk_vtx_cut_nue_cc_mixed,
	        h_ele_momentum_trk_vtx_cut_nue_cc_out_fv,
	        h_ele_momentum_trk_vtx_cut_numu_cc,
	        h_ele_momentum_trk_vtx_cut_numu_cc_mixed,
	        h_ele_momentum_trk_vtx_cut_cosmic,
	        h_ele_momentum_trk_vtx_cut_nc,
	        h_ele_momentum_trk_vtx_cut_nc_pi0,
	        h_ele_momentum_trk_vtx_cut_other_mixed,
	        h_ele_momentum_trk_vtx_cut_unmatched,
	        h_ele_momentum_trk_vtx_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_trk_vtx_cut_data,
	        data_scale_factor,
	        "Track Vertex Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_trk_vtx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_trk_vtx_electron,
	        h_leading_momentum_trk_vtx_proton,
	        h_leading_momentum_trk_vtx_photon,
	        h_leading_momentum_trk_vtx_pion,
	        h_leading_momentum_trk_vtx_kaon,
	        h_leading_momentum_trk_vtx_muon,
	        h_leading_momentum_trk_vtx_neutron,
	        h_leading_momentum_trk_vtx_mc_unmatched,
	        h_leading_momentum_trk_vtx_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_trk_vtx.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_hit_cut_nue_cc,
	        h_ele_momentum_hit_cut_nue_cc_mixed,
	        h_ele_momentum_hit_cut_nue_cc_out_fv,
	        h_ele_momentum_hit_cut_numu_cc,
	        h_ele_momentum_hit_cut_numu_cc_mixed,
	        h_ele_momentum_hit_cut_cosmic,
	        h_ele_momentum_hit_cut_nc,
	        h_ele_momentum_hit_cut_nc_pi0,
	        h_ele_momentum_hit_cut_other_mixed,
	        h_ele_momentum_hit_cut_unmatched,
	        h_ele_momentum_hit_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_hit_cut_data,
	        data_scale_factor,
	        "Total Hit Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_hit_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_hit_cut_electron,
	        h_leading_momentum_hit_cut_proton,
	        h_leading_momentum_hit_cut_photon,
	        h_leading_momentum_hit_cut_pion,
	        h_leading_momentum_hit_cut_kaon,
	        h_leading_momentum_hit_cut_muon,
	        h_leading_momentum_hit_cut_neutron,
	        h_leading_momentum_hit_cut_mc_unmatched,
	        h_leading_momentum_hit_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_hit_cut.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_yhit_cut_nue_cc,
	        h_ele_momentum_yhit_cut_nue_cc_mixed,
	        h_ele_momentum_yhit_cut_nue_cc_out_fv,
	        h_ele_momentum_yhit_cut_numu_cc,
	        h_ele_momentum_yhit_cut_numu_cc_mixed,
	        h_ele_momentum_yhit_cut_cosmic,
	        h_ele_momentum_yhit_cut_nc,
	        h_ele_momentum_yhit_cut_nc_pi0,
	        h_ele_momentum_yhit_cut_other_mixed,
	        h_ele_momentum_yhit_cut_unmatched,
	        h_ele_momentum_yhit_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_yhit_cut_data,
	        data_scale_factor,
	        "Collection Hit Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_yhit_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_yhit_cut_electron,
	        h_leading_momentum_yhit_cut_proton,
	        h_leading_momentum_yhit_cut_photon,
	        h_leading_momentum_yhit_cut_pion,
	        h_leading_momentum_yhit_cut_kaon,
	        h_leading_momentum_yhit_cut_muon,
	        h_leading_momentum_yhit_cut_neutron,
	        h_leading_momentum_yhit_cut_mc_unmatched,
	        h_leading_momentum_yhit_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_yhit_cut.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_open_angle_cut_nue_cc,
	        h_ele_momentum_open_angle_cut_nue_cc_mixed,
	        h_ele_momentum_open_angle_cut_nue_cc_out_fv,
	        h_ele_momentum_open_angle_cut_numu_cc,
	        h_ele_momentum_open_angle_cut_numu_cc_mixed,
	        h_ele_momentum_open_angle_cut_cosmic,
	        h_ele_momentum_open_angle_cut_nc,
	        h_ele_momentum_open_angle_cut_nc_pi0,
	        h_ele_momentum_open_angle_cut_other_mixed,
	        h_ele_momentum_open_angle_cut_unmatched,
	        h_ele_momentum_open_angle_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_open_angle_cut_data,
	        data_scale_factor,
	        "Open Angle Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_open_angle_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_open_angle_cut_electron,
	        h_leading_momentum_open_angle_cut_proton,
	        h_leading_momentum_open_angle_cut_photon,
	        h_leading_momentum_open_angle_cut_pion,
	        h_leading_momentum_open_angle_cut_kaon,
	        h_leading_momentum_open_angle_cut_muon,
	        h_leading_momentum_open_angle_cut_neutron,
	        h_leading_momentum_open_angle_cut_mc_unmatched,
	        h_leading_momentum_open_angle_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_open_angle_cut.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_dedx_cut_nue_cc,
	        h_ele_momentum_dedx_cut_nue_cc_mixed,
	        h_ele_momentum_dedx_cut_nue_cc_out_fv,
	        h_ele_momentum_dedx_cut_numu_cc,
	        h_ele_momentum_dedx_cut_numu_cc_mixed,
	        h_ele_momentum_dedx_cut_cosmic,
	        h_ele_momentum_dedx_cut_nc,
	        h_ele_momentum_dedx_cut_nc_pi0,
	        h_ele_momentum_dedx_cut_other_mixed,
	        h_ele_momentum_dedx_cut_unmatched,
	        h_ele_momentum_dedx_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_dedx_cut_data,
	        data_scale_factor,
	        "dE/dx Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_dedx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_dedx_cut_electron,
	        h_leading_momentum_dedx_cut_proton,
	        h_leading_momentum_dedx_cut_photon,
	        h_leading_momentum_dedx_cut_pion,
	        h_leading_momentum_dedx_cut_kaon,
	        h_leading_momentum_dedx_cut_muon,
	        h_leading_momentum_dedx_cut_neutron,
	        h_leading_momentum_dedx_cut_mc_unmatched,
	        h_leading_momentum_dedx_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_dedx_cut.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_2shwr_cut_nue_cc,
	        h_ele_momentum_2shwr_cut_nue_cc_mixed,
	        h_ele_momentum_2shwr_cut_nue_cc_out_fv,
	        h_ele_momentum_2shwr_cut_numu_cc,
	        h_ele_momentum_2shwr_cut_numu_cc_mixed,
	        h_ele_momentum_2shwr_cut_cosmic,
	        h_ele_momentum_2shwr_cut_nc,
	        h_ele_momentum_2shwr_cut_nc_pi0,
	        h_ele_momentum_2shwr_cut_other_mixed,
	        h_ele_momentum_2shwr_cut_unmatched,
	        h_ele_momentum_2shwr_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_2shwr_cut_data,
	        data_scale_factor,
	        "Secondary Shower Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_2shwr_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_2shwr_cut_electron,
	        h_leading_momentum_2shwr_cut_proton,
	        h_leading_momentum_2shwr_cut_photon,
	        h_leading_momentum_2shwr_cut_pion,
	        h_leading_momentum_2shwr_cut_kaon,
	        h_leading_momentum_2shwr_cut_muon,
	        h_leading_momentum_2shwr_cut_neutron,
	        h_leading_momentum_2shwr_cut_mc_unmatched,
	        h_leading_momentum_2shwr_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_2shwr_cut.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_hit_length_cut_nue_cc,
	        h_ele_momentum_hit_length_cut_nue_cc_mixed,
	        h_ele_momentum_hit_length_cut_nue_cc_out_fv,
	        h_ele_momentum_hit_length_cut_numu_cc,
	        h_ele_momentum_hit_length_cut_numu_cc_mixed,
	        h_ele_momentum_hit_length_cut_cosmic,
	        h_ele_momentum_hit_length_cut_nc,
	        h_ele_momentum_hit_length_cut_nc_pi0,
	        h_ele_momentum_hit_length_cut_other_mixed,
	        h_ele_momentum_hit_length_cut_unmatched,
	        h_ele_momentum_hit_length_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_hit_length_cut_data,
	        data_scale_factor,
	        "Hit/Length Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_hit_length_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_hit_length_cut_electron,
	        h_leading_momentum_hit_length_cut_proton,
	        h_leading_momentum_hit_length_cut_photon,
	        h_leading_momentum_hit_length_cut_pion,
	        h_leading_momentum_hit_length_cut_kaon,
	        h_leading_momentum_hit_length_cut_muon,
	        h_leading_momentum_hit_length_cut_neutron,
	        h_leading_momentum_hit_length_cut_mc_unmatched,
	        h_leading_momentum_hit_length_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_hit_length_cut.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_length_ratio_cut_nue_cc,
	        h_ele_momentum_length_ratio_cut_nue_cc_mixed,
	        h_ele_momentum_length_ratio_cut_nue_cc_out_fv,
	        h_ele_momentum_length_ratio_cut_numu_cc,
	        h_ele_momentum_length_ratio_cut_numu_cc_mixed,
	        h_ele_momentum_length_ratio_cut_cosmic,
	        h_ele_momentum_length_ratio_cut_nc,
	        h_ele_momentum_length_ratio_cut_nc_pi0,
	        h_ele_momentum_length_ratio_cut_other_mixed,
	        h_ele_momentum_length_ratio_cut_unmatched,
	        h_ele_momentum_length_ratio_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_length_ratio_cut_data,
	        data_scale_factor,
	        "Trk Length / Shower Length Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_length_ratio_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_length_ratio_cut_electron,
	        h_leading_momentum_length_ratio_cut_proton,
	        h_leading_momentum_length_ratio_cut_photon,
	        h_leading_momentum_length_ratio_cut_pion,
	        h_leading_momentum_length_ratio_cut_kaon,
	        h_leading_momentum_length_ratio_cut_muon,
	        h_leading_momentum_length_ratio_cut_neutron,
	        h_leading_momentum_length_ratio_cut_mc_unmatched,
	        h_leading_momentum_length_ratio_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_length_ratio_cut.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_containment_cut_nue_cc,
	        h_ele_momentum_containment_cut_nue_cc_mixed,
	        h_ele_momentum_containment_cut_nue_cc_out_fv,
	        h_ele_momentum_containment_cut_numu_cc,
	        h_ele_momentum_containment_cut_numu_cc_mixed,
	        h_ele_momentum_containment_cut_cosmic,
	        h_ele_momentum_containment_cut_nc,
	        h_ele_momentum_containment_cut_nc_pi0,
	        h_ele_momentum_containment_cut_other_mixed,
	        h_ele_momentum_containment_cut_unmatched,
	        h_ele_momentum_containment_cut_intime,
	        intime_scale_factor,
	        h_ele_momentum_containment_cut_data,
	        data_scale_factor,
	        "Track Containment Cut", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_containment_cut_data.pdf"));

	histogram_functions::PlotSimpleStackParticle(
	        h_leading_momentum_containment_cut_electron,
	        h_leading_momentum_containment_cut_proton,
	        h_leading_momentum_containment_cut_photon,
	        h_leading_momentum_containment_cut_pion,
	        h_leading_momentum_containment_cut_kaon,
	        h_leading_momentum_containment_cut_muon,
	        h_leading_momentum_containment_cut_neutron,
	        h_leading_momentum_containment_cut_mc_unmatched,
	        h_leading_momentum_containment_cut_ext_unmatched,
	        intime_scale_factor, data_scale_factor,
	        0.73, 0.98, 0.98, 0.50,
	        "",
	        "(True Particle) Leading Shower Momentum [GeV]",
	        "", Form("%s%s", file_locate_prefix, "leading_momentum_true_particle_containment_cut.pdf"));

	histogram_functions::PlotdEdxTheta(
	        h_dedx_theta_nue_cc,
	        h_dedx_theta_nue_cc_mixed,
	        h_dedx_theta_nue_cc_out_fv,
	        h_dedx_theta_numu_cc,
	        h_dedx_theta_numu_cc_mixed,
	        h_dedx_theta_cosmic,
	        h_dedx_theta_nc,
	        h_dedx_theta_nc_pi0,
	        h_dedx_theta_other_mixed,
	        h_dedx_theta_unmatched,
	        h_dedx_theta_intime,
	        intime_scale_factor,
	        h_dedx_theta_data,
	        data_scale_factor,
	        "", "Leading Shower dE/dx [MeV/cm]", "Leading Shower Theta [Degrees]",
	        Form("%s%s", file_locate_prefix, "post_cuts_dEdx_theta_mc_ext.pdf"),
	        Form("%s%s", file_locate_prefix, "post_cuts_dEdx_theta_data.pdf"),
	        Form("%s%s", file_locate_prefix, "post_cuts_dEdx_theta_diff.pdf"));

	histogram_functions::PlotdEdxTheta(
	        h_dedx_theta_pre_cuts_nue_cc,
	        h_dedx_theta_pre_cuts_nue_cc_mixed,
	        h_dedx_theta_pre_cuts_nue_cc_out_fv,
	        h_dedx_theta_pre_cuts_numu_cc,
	        h_dedx_theta_pre_cuts_numu_cc_mixed,
	        h_dedx_theta_pre_cuts_cosmic,
	        h_dedx_theta_pre_cuts_nc,
	        h_dedx_theta_pre_cuts_nc_pi0,
	        h_dedx_theta_pre_cuts_other_mixed,
	        h_dedx_theta_pre_cuts_unmatched,
	        h_dedx_theta_pre_cuts_intime,
	        intime_scale_factor,
	        h_dedx_theta_pre_cuts_data,
	        data_scale_factor,
	        "", "Leading Shower dE/dx [MeV/cm]", "Leading Shower Theta [Degrees]",
	        Form("%s%s", file_locate_prefix, "post_cuts_dEdx_theta_pre_cuts_mc_ext.pdf"),
	        Form("%s%s", file_locate_prefix, "post_cuts_dEdx_theta_pre_cuts_data.pdf"),
	        Form("%s%s", file_locate_prefix, "post_cuts_dEdx_theta_pre_cuts_diff.pdf"));

	histogram_functions::Plot2DHistogram(h_pfp_zy_vtx_nue_cc, "Signal", "Reco Pfp Vtx Z [cm]", "Reco pfp Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "pfp_zy_vtx_nue_cc.pdf"));
	histogram_functions::Plot2DHistogram(h_pfp_zy_vtx_all, "All", "Reco Pfp Vtx Z [cm]", "Reco pfp Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "pfp_zy_vtx_all.pdf"));
	histogram_functions::Plot2DHistogram(h_pfp_zy_vtx_ext, "EXT", "Reco Pfp Vtx Z [cm]", "Reco pfp Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "pfp_zy_vtx_ext.pdf"));
	histogram_functions::Plot2DHistogram(h_pfp_zy_vtx_data, "Data", "Reco Pfp Vtx Z [cm]", "Reco pfp Vtx Y [cm]",
	                                     Form("%s%s", file_locate_prefix, "pfp_zy_vtx_data.pdf"));

	histogram_functions::Plot2DHistogramNormZ(h_pfp_zy_vtx_all, h_pfp_zy_vtx_data, "All", "Data", "Reco Pfp Vtx Z [cm]", "Reco pfp Vtx Y [cm]",
	                                          Form("%s%s", file_locate_prefix, "pfp_zy_vtx_all_norm_z.pdf"),
	                                          Form("%s%s", file_locate_prefix, "pfp_zy_vtx_data_norm_z.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_nue_cut_nue_cc,
	        h_multiplicity_shower_nue_cut_nue_cc_mixed,
	        h_multiplicity_shower_nue_cut_nue_cc_out_fv,
	        h_multiplicity_shower_nue_cut_numu_cc,
	        h_multiplicity_shower_nue_cut_numu_cc_mixed,
	        h_multiplicity_shower_nue_cut_cosmic,
	        h_multiplicity_shower_nue_cut_nc,
	        h_multiplicity_shower_nue_cut_nc_pi0,
	        h_multiplicity_shower_nue_cut_other_mixed,
	        h_multiplicity_shower_nue_cut_unmatched,
	        h_multiplicity_shower_nue_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_nue_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_nue_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_nue_cut_nue_cc,
	        h_multiplicity_track_nue_cut_nue_cc_mixed,
	        h_multiplicity_track_nue_cut_nue_cc_out_fv,
	        h_multiplicity_track_nue_cut_numu_cc,
	        h_multiplicity_track_nue_cut_numu_cc_mixed,
	        h_multiplicity_track_nue_cut_cosmic,
	        h_multiplicity_track_nue_cut_nc,
	        h_multiplicity_track_nue_cut_nc_pi0,
	        h_multiplicity_track_nue_cut_other_mixed,
	        h_multiplicity_track_nue_cut_unmatched,
	        h_multiplicity_track_nue_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_nue_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_nue_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_fv_cut_nue_cc,
	        h_multiplicity_shower_fv_cut_nue_cc_mixed,
	        h_multiplicity_shower_fv_cut_nue_cc_out_fv,
	        h_multiplicity_shower_fv_cut_numu_cc,
	        h_multiplicity_shower_fv_cut_numu_cc_mixed,
	        h_multiplicity_shower_fv_cut_cosmic,
	        h_multiplicity_shower_fv_cut_nc,
	        h_multiplicity_shower_fv_cut_nc_pi0,
	        h_multiplicity_shower_fv_cut_other_mixed,
	        h_multiplicity_shower_fv_cut_unmatched,
	        h_multiplicity_shower_fv_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_fv_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_fv_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_fv_cut_nue_cc,
	        h_multiplicity_track_fv_cut_nue_cc_mixed,
	        h_multiplicity_track_fv_cut_nue_cc_out_fv,
	        h_multiplicity_track_fv_cut_numu_cc,
	        h_multiplicity_track_fv_cut_numu_cc_mixed,
	        h_multiplicity_track_fv_cut_cosmic,
	        h_multiplicity_track_fv_cut_nc,
	        h_multiplicity_track_fv_cut_nc_pi0,
	        h_multiplicity_track_fv_cut_other_mixed,
	        h_multiplicity_track_fv_cut_unmatched,
	        h_multiplicity_track_fv_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_fv_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_fv_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_flash_vtx_cut_nue_cc,
	        h_multiplicity_shower_flash_vtx_cut_nue_cc_mixed,
	        h_multiplicity_shower_flash_vtx_cut_nue_cc_out_fv,
	        h_multiplicity_shower_flash_vtx_cut_numu_cc,
	        h_multiplicity_shower_flash_vtx_cut_numu_cc_mixed,
	        h_multiplicity_shower_flash_vtx_cut_cosmic,
	        h_multiplicity_shower_flash_vtx_cut_nc,
	        h_multiplicity_shower_flash_vtx_cut_nc_pi0,
	        h_multiplicity_shower_flash_vtx_cut_other_mixed,
	        h_multiplicity_shower_flash_vtx_cut_unmatched,
	        h_multiplicity_shower_flash_vtx_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_flash_vtx_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_flash_vtx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_flash_vtx_cut_nue_cc,
	        h_multiplicity_track_flash_vtx_cut_nue_cc_mixed,
	        h_multiplicity_track_flash_vtx_cut_nue_cc_out_fv,
	        h_multiplicity_track_flash_vtx_cut_numu_cc,
	        h_multiplicity_track_flash_vtx_cut_numu_cc_mixed,
	        h_multiplicity_track_flash_vtx_cut_cosmic,
	        h_multiplicity_track_flash_vtx_cut_nc,
	        h_multiplicity_track_flash_vtx_cut_nc_pi0,
	        h_multiplicity_track_flash_vtx_cut_other_mixed,
	        h_multiplicity_track_flash_vtx_cut_unmatched,
	        h_multiplicity_track_flash_vtx_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_flash_vtx_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_flash_vtx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_shwr_vtx_cut_nue_cc,
	        h_multiplicity_shower_shwr_vtx_cut_nue_cc_mixed,
	        h_multiplicity_shower_shwr_vtx_cut_nue_cc_out_fv,
	        h_multiplicity_shower_shwr_vtx_cut_numu_cc,
	        h_multiplicity_shower_shwr_vtx_cut_numu_cc_mixed,
	        h_multiplicity_shower_shwr_vtx_cut_cosmic,
	        h_multiplicity_shower_shwr_vtx_cut_nc,
	        h_multiplicity_shower_shwr_vtx_cut_nc_pi0,
	        h_multiplicity_shower_shwr_vtx_cut_other_mixed,
	        h_multiplicity_shower_shwr_vtx_cut_unmatched,
	        h_multiplicity_shower_shwr_vtx_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_shwr_vtx_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_shwr_vtx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_shwr_vtx_cut_nue_cc,
	        h_multiplicity_track_shwr_vtx_cut_nue_cc_mixed,
	        h_multiplicity_track_shwr_vtx_cut_nue_cc_out_fv,
	        h_multiplicity_track_shwr_vtx_cut_numu_cc,
	        h_multiplicity_track_shwr_vtx_cut_numu_cc_mixed,
	        h_multiplicity_track_shwr_vtx_cut_cosmic,
	        h_multiplicity_track_shwr_vtx_cut_nc,
	        h_multiplicity_track_shwr_vtx_cut_nc_pi0,
	        h_multiplicity_track_shwr_vtx_cut_other_mixed,
	        h_multiplicity_track_shwr_vtx_cut_unmatched,
	        h_multiplicity_track_shwr_vtx_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_shwr_vtx_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_shwr_vtx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_trk_vtx_cut_nue_cc,
	        h_multiplicity_shower_trk_vtx_cut_nue_cc_mixed,
	        h_multiplicity_shower_trk_vtx_cut_nue_cc_out_fv,
	        h_multiplicity_shower_trk_vtx_cut_numu_cc,
	        h_multiplicity_shower_trk_vtx_cut_numu_cc_mixed,
	        h_multiplicity_shower_trk_vtx_cut_cosmic,
	        h_multiplicity_shower_trk_vtx_cut_nc,
	        h_multiplicity_shower_trk_vtx_cut_nc_pi0,
	        h_multiplicity_shower_trk_vtx_cut_other_mixed,
	        h_multiplicity_shower_trk_vtx_cut_unmatched,
	        h_multiplicity_shower_trk_vtx_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_trk_vtx_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_trk_vtx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_trk_vtx_cut_nue_cc,
	        h_multiplicity_track_trk_vtx_cut_nue_cc_mixed,
	        h_multiplicity_track_trk_vtx_cut_nue_cc_out_fv,
	        h_multiplicity_track_trk_vtx_cut_numu_cc,
	        h_multiplicity_track_trk_vtx_cut_numu_cc_mixed,
	        h_multiplicity_track_trk_vtx_cut_cosmic,
	        h_multiplicity_track_trk_vtx_cut_nc,
	        h_multiplicity_track_trk_vtx_cut_nc_pi0,
	        h_multiplicity_track_trk_vtx_cut_other_mixed,
	        h_multiplicity_track_trk_vtx_cut_unmatched,
	        h_multiplicity_track_trk_vtx_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_trk_vtx_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_trk_vtx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_hit_cut_nue_cc,
	        h_multiplicity_shower_hit_cut_nue_cc_mixed,
	        h_multiplicity_shower_hit_cut_nue_cc_out_fv,
	        h_multiplicity_shower_hit_cut_numu_cc,
	        h_multiplicity_shower_hit_cut_numu_cc_mixed,
	        h_multiplicity_shower_hit_cut_cosmic,
	        h_multiplicity_shower_hit_cut_nc,
	        h_multiplicity_shower_hit_cut_nc_pi0,
	        h_multiplicity_shower_hit_cut_other_mixed,
	        h_multiplicity_shower_hit_cut_unmatched,
	        h_multiplicity_shower_hit_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_hit_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_hit_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_hit_cut_nue_cc,
	        h_multiplicity_track_hit_cut_nue_cc_mixed,
	        h_multiplicity_track_hit_cut_nue_cc_out_fv,
	        h_multiplicity_track_hit_cut_numu_cc,
	        h_multiplicity_track_hit_cut_numu_cc_mixed,
	        h_multiplicity_track_hit_cut_cosmic,
	        h_multiplicity_track_hit_cut_nc,
	        h_multiplicity_track_hit_cut_nc_pi0,
	        h_multiplicity_track_hit_cut_other_mixed,
	        h_multiplicity_track_hit_cut_unmatched,
	        h_multiplicity_track_hit_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_hit_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_hit_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_yhit_cut_nue_cc,
	        h_multiplicity_shower_yhit_cut_nue_cc_mixed,
	        h_multiplicity_shower_yhit_cut_nue_cc_out_fv,
	        h_multiplicity_shower_yhit_cut_numu_cc,
	        h_multiplicity_shower_yhit_cut_numu_cc_mixed,
	        h_multiplicity_shower_yhit_cut_cosmic,
	        h_multiplicity_shower_yhit_cut_nc,
	        h_multiplicity_shower_yhit_cut_nc_pi0,
	        h_multiplicity_shower_yhit_cut_other_mixed,
	        h_multiplicity_shower_yhit_cut_unmatched,
	        h_multiplicity_shower_yhit_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_yhit_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_yhit_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_yhit_cut_nue_cc,
	        h_multiplicity_track_yhit_cut_nue_cc_mixed,
	        h_multiplicity_track_yhit_cut_nue_cc_out_fv,
	        h_multiplicity_track_yhit_cut_numu_cc,
	        h_multiplicity_track_yhit_cut_numu_cc_mixed,
	        h_multiplicity_track_yhit_cut_cosmic,
	        h_multiplicity_track_yhit_cut_nc,
	        h_multiplicity_track_yhit_cut_nc_pi0,
	        h_multiplicity_track_yhit_cut_other_mixed,
	        h_multiplicity_track_yhit_cut_unmatched,
	        h_multiplicity_track_yhit_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_yhit_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_yhit_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_open_angle_cut_nue_cc,
	        h_multiplicity_shower_open_angle_cut_nue_cc_mixed,
	        h_multiplicity_shower_open_angle_cut_nue_cc_out_fv,
	        h_multiplicity_shower_open_angle_cut_numu_cc,
	        h_multiplicity_shower_open_angle_cut_numu_cc_mixed,
	        h_multiplicity_shower_open_angle_cut_cosmic,
	        h_multiplicity_shower_open_angle_cut_nc,
	        h_multiplicity_shower_open_angle_cut_nc_pi0,
	        h_multiplicity_shower_open_angle_cut_other_mixed,
	        h_multiplicity_shower_open_angle_cut_unmatched,
	        h_multiplicity_shower_open_angle_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_open_angle_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_open_angle_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_open_angle_cut_nue_cc,
	        h_multiplicity_track_open_angle_cut_nue_cc_mixed,
	        h_multiplicity_track_open_angle_cut_nue_cc_out_fv,
	        h_multiplicity_track_open_angle_cut_numu_cc,
	        h_multiplicity_track_open_angle_cut_numu_cc_mixed,
	        h_multiplicity_track_open_angle_cut_cosmic,
	        h_multiplicity_track_open_angle_cut_nc,
	        h_multiplicity_track_open_angle_cut_nc_pi0,
	        h_multiplicity_track_open_angle_cut_other_mixed,
	        h_multiplicity_track_open_angle_cut_unmatched,
	        h_multiplicity_track_open_angle_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_open_angle_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_open_angle_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_dedx_cut_nue_cc,
	        h_multiplicity_shower_dedx_cut_nue_cc_mixed,
	        h_multiplicity_shower_dedx_cut_nue_cc_out_fv,
	        h_multiplicity_shower_dedx_cut_numu_cc,
	        h_multiplicity_shower_dedx_cut_numu_cc_mixed,
	        h_multiplicity_shower_dedx_cut_cosmic,
	        h_multiplicity_shower_dedx_cut_nc,
	        h_multiplicity_shower_dedx_cut_nc_pi0,
	        h_multiplicity_shower_dedx_cut_other_mixed,
	        h_multiplicity_shower_dedx_cut_unmatched,
	        h_multiplicity_shower_dedx_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_dedx_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_dedx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_dedx_cut_nue_cc,
	        h_multiplicity_track_dedx_cut_nue_cc_mixed,
	        h_multiplicity_track_dedx_cut_nue_cc_out_fv,
	        h_multiplicity_track_dedx_cut_numu_cc,
	        h_multiplicity_track_dedx_cut_numu_cc_mixed,
	        h_multiplicity_track_dedx_cut_cosmic,
	        h_multiplicity_track_dedx_cut_nc,
	        h_multiplicity_track_dedx_cut_nc_pi0,
	        h_multiplicity_track_dedx_cut_other_mixed,
	        h_multiplicity_track_dedx_cut_unmatched,
	        h_multiplicity_track_dedx_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_dedx_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_dedx_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_2shwr_cut_nue_cc,
	        h_multiplicity_shower_2shwr_cut_nue_cc_mixed,
	        h_multiplicity_shower_2shwr_cut_nue_cc_out_fv,
	        h_multiplicity_shower_2shwr_cut_numu_cc,
	        h_multiplicity_shower_2shwr_cut_numu_cc_mixed,
	        h_multiplicity_shower_2shwr_cut_cosmic,
	        h_multiplicity_shower_2shwr_cut_nc,
	        h_multiplicity_shower_2shwr_cut_nc_pi0,
	        h_multiplicity_shower_2shwr_cut_other_mixed,
	        h_multiplicity_shower_2shwr_cut_unmatched,
	        h_multiplicity_shower_2shwr_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_2shwr_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_2shwr_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_2shwr_cut_nue_cc,
	        h_multiplicity_track_2shwr_cut_nue_cc_mixed,
	        h_multiplicity_track_2shwr_cut_nue_cc_out_fv,
	        h_multiplicity_track_2shwr_cut_numu_cc,
	        h_multiplicity_track_2shwr_cut_numu_cc_mixed,
	        h_multiplicity_track_2shwr_cut_cosmic,
	        h_multiplicity_track_2shwr_cut_nc,
	        h_multiplicity_track_2shwr_cut_nc_pi0,
	        h_multiplicity_track_2shwr_cut_other_mixed,
	        h_multiplicity_track_2shwr_cut_unmatched,
	        h_multiplicity_track_2shwr_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_2shwr_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_2shwr_cut_data.pdf"));


	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_hit_length_cut_nue_cc,
	        h_multiplicity_shower_hit_length_cut_nue_cc_mixed,
	        h_multiplicity_shower_hit_length_cut_nue_cc_out_fv,
	        h_multiplicity_shower_hit_length_cut_numu_cc,
	        h_multiplicity_shower_hit_length_cut_numu_cc_mixed,
	        h_multiplicity_shower_hit_length_cut_cosmic,
	        h_multiplicity_shower_hit_length_cut_nc,
	        h_multiplicity_shower_hit_length_cut_nc_pi0,
	        h_multiplicity_shower_hit_length_cut_other_mixed,
	        h_multiplicity_shower_hit_length_cut_unmatched,
	        h_multiplicity_shower_hit_length_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_hit_length_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_hit_length_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_hit_length_cut_nue_cc,
	        h_multiplicity_track_hit_length_cut_nue_cc_mixed,
	        h_multiplicity_track_hit_length_cut_nue_cc_out_fv,
	        h_multiplicity_track_hit_length_cut_numu_cc,
	        h_multiplicity_track_hit_length_cut_numu_cc_mixed,
	        h_multiplicity_track_hit_length_cut_cosmic,
	        h_multiplicity_track_hit_length_cut_nc,
	        h_multiplicity_track_hit_length_cut_nc_pi0,
	        h_multiplicity_track_hit_length_cut_other_mixed,
	        h_multiplicity_track_hit_length_cut_unmatched,
	        h_multiplicity_track_hit_length_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_hit_length_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_hit_length_cut_data.pdf"));


	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_length_ratio_cut_nue_cc,
	        h_multiplicity_shower_length_ratio_cut_nue_cc_mixed,
	        h_multiplicity_shower_length_ratio_cut_nue_cc_out_fv,
	        h_multiplicity_shower_length_ratio_cut_numu_cc,
	        h_multiplicity_shower_length_ratio_cut_numu_cc_mixed,
	        h_multiplicity_shower_length_ratio_cut_cosmic,
	        h_multiplicity_shower_length_ratio_cut_nc,
	        h_multiplicity_shower_length_ratio_cut_nc_pi0,
	        h_multiplicity_shower_length_ratio_cut_other_mixed,
	        h_multiplicity_shower_length_ratio_cut_unmatched,
	        h_multiplicity_shower_length_ratio_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_length_ratio_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_length_ratio_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_length_ratio_cut_nue_cc,
	        h_multiplicity_track_length_ratio_cut_nue_cc_mixed,
	        h_multiplicity_track_length_ratio_cut_nue_cc_out_fv,
	        h_multiplicity_track_length_ratio_cut_numu_cc,
	        h_multiplicity_track_length_ratio_cut_numu_cc_mixed,
	        h_multiplicity_track_length_ratio_cut_cosmic,
	        h_multiplicity_track_length_ratio_cut_nc,
	        h_multiplicity_track_length_ratio_cut_nc_pi0,
	        h_multiplicity_track_length_ratio_cut_other_mixed,
	        h_multiplicity_track_length_ratio_cut_unmatched,
	        h_multiplicity_track_length_ratio_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_length_ratio_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_length_ratio_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_shower_containment_cut_nue_cc,
	        h_multiplicity_shower_containment_cut_nue_cc_mixed,
	        h_multiplicity_shower_containment_cut_nue_cc_out_fv,
	        h_multiplicity_shower_containment_cut_numu_cc,
	        h_multiplicity_shower_containment_cut_numu_cc_mixed,
	        h_multiplicity_shower_containment_cut_cosmic,
	        h_multiplicity_shower_containment_cut_nc,
	        h_multiplicity_shower_containment_cut_nc_pi0,
	        h_multiplicity_shower_containment_cut_other_mixed,
	        h_multiplicity_shower_containment_cut_unmatched,
	        h_multiplicity_shower_containment_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_shower_containment_cut_data,
	        data_scale_factor,
	        "", "Selected Shower Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_shower_multiplicity_containment_cut_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_multiplicity_track_containment_cut_nue_cc,
	        h_multiplicity_track_containment_cut_nue_cc_mixed,
	        h_multiplicity_track_containment_cut_nue_cc_out_fv,
	        h_multiplicity_track_containment_cut_numu_cc,
	        h_multiplicity_track_containment_cut_numu_cc_mixed,
	        h_multiplicity_track_containment_cut_cosmic,
	        h_multiplicity_track_containment_cut_nc,
	        h_multiplicity_track_containment_cut_nc_pi0,
	        h_multiplicity_track_containment_cut_other_mixed,
	        h_multiplicity_track_containment_cut_unmatched,
	        h_multiplicity_track_containment_cut_intime,
	        intime_scale_factor,
	        h_multiplicity_track_containment_cut_data,
	        data_scale_factor,
	        "", "Selected Track Multiplicity", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_track_multiplicity_containment_cut_data.pdf"));

	histogram_functions::Plot1DHistogramGausFit(h_low_true_momentum, "Selected Electron Momentum (True) [GeV]",
	                                            Form("%s%s", file_locate_prefix, "true_electron_momentum_low.pdf"));
	histogram_functions::Plot1DHistogramGausFit(h_med_true_momentum, "Selected Electron Momentum (True) [GeV]",
	                                            Form("%s%s", file_locate_prefix, "true_electron_momentum_med.pdf"));
	histogram_functions::Plot1DHistogramGausFit(h_high_true_momentum, "Selected Electron Momentum (True) [GeV]",
	                                            Form("%s%s", file_locate_prefix, "true_electron_momentum_high.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_slice_1_nue_cc,
	        h_ele_momentum_slice_1_nue_cc_mixed,
	        h_ele_momentum_slice_1_nue_cc_out_fv,
	        h_ele_momentum_slice_1_numu_cc,
	        h_ele_momentum_slice_1_numu_cc_mixed,
	        h_ele_momentum_slice_1_cosmic,
	        h_ele_momentum_slice_1_nc,
	        h_ele_momentum_slice_1_nc_pi0,
	        h_ele_momentum_slice_1_other_mixed,
	        h_ele_momentum_slice_1_unmatched,
	        h_ele_momentum_slice_1_intime,
	        intime_scale_factor,
	        h_ele_momentum_slice_1_data,
	        data_scale_factor,
	        "Theta Slice (0 - 40)", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_theta_slice_1_data.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_slice_2_nue_cc,
	        h_ele_momentum_slice_2_nue_cc_mixed,
	        h_ele_momentum_slice_2_nue_cc_out_fv,
	        h_ele_momentum_slice_2_numu_cc,
	        h_ele_momentum_slice_2_numu_cc_mixed,
	        h_ele_momentum_slice_2_cosmic,
	        h_ele_momentum_slice_2_nc,
	        h_ele_momentum_slice_2_nc_pi0,
	        h_ele_momentum_slice_2_other_mixed,
	        h_ele_momentum_slice_2_unmatched,
	        h_ele_momentum_slice_2_intime,
	        intime_scale_factor,
	        h_ele_momentum_slice_2_data,
	        data_scale_factor,
	        "Theta Slice (40 - 90)", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_theta_slice_2_data.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_momentum_slice_3_nue_cc,
	        h_ele_momentum_slice_3_nue_cc_mixed,
	        h_ele_momentum_slice_3_nue_cc_out_fv,
	        h_ele_momentum_slice_3_numu_cc,
	        h_ele_momentum_slice_3_numu_cc_mixed,
	        h_ele_momentum_slice_3_cosmic,
	        h_ele_momentum_slice_3_nc,
	        h_ele_momentum_slice_3_nc_pi0,
	        h_ele_momentum_slice_3_other_mixed,
	        h_ele_momentum_slice_3_unmatched,
	        h_ele_momentum_slice_3_intime,
	        intime_scale_factor,
	        h_ele_momentum_slice_3_data,
	        data_scale_factor,
	        "Theta Slice (90 - 180)", "Leading Shower Momentum [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_ele_momentum_theta_slice_3_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_dedx_slice_1_nue_cc,
	        h_dedx_slice_1_nue_cc_mixed,
	        h_dedx_slice_1_nue_cc_out_fv,
	        h_dedx_slice_1_numu_cc,
	        h_dedx_slice_1_numu_cc_mixed,
	        h_dedx_slice_1_cosmic,
	        h_dedx_slice_1_nc,
	        h_dedx_slice_1_nc_pi0,
	        h_dedx_slice_1_other_mixed,
	        h_dedx_slice_1_unmatched,
	        h_dedx_slice_1_intime,
	        intime_scale_factor,
	        h_dedx_slice_1_data,
	        data_scale_factor,
	        "Theta Slice (0 - 60)", "Leading Shower dE/dx [MeV/cm]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_1_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_dedx_slice_1_zoom_nue_cc,
	        h_dedx_slice_1_zoom_nue_cc_mixed,
	        h_dedx_slice_1_zoom_nue_cc_out_fv,
	        h_dedx_slice_1_zoom_numu_cc,
	        h_dedx_slice_1_zoom_numu_cc_mixed,
	        h_dedx_slice_1_zoom_cosmic,
	        h_dedx_slice_1_zoom_nc,
	        h_dedx_slice_1_zoom_nc_pi0,
	        h_dedx_slice_1_zoom_other_mixed,
	        h_dedx_slice_1_zoom_unmatched,
	        h_dedx_slice_1_zoom_intime,
	        intime_scale_factor,
	        h_dedx_slice_1_zoom_data,
	        data_scale_factor,
	        "Theta Slice (0 - 60)", "Leading Shower dE/dx [MeV/cm]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_1_zoom_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_dedx_slice_2_nue_cc,
	        h_dedx_slice_2_nue_cc_mixed,
	        h_dedx_slice_2_nue_cc_out_fv,
	        h_dedx_slice_2_numu_cc,
	        h_dedx_slice_2_numu_cc_mixed,
	        h_dedx_slice_2_cosmic,
	        h_dedx_slice_2_nc,
	        h_dedx_slice_2_nc_pi0,
	        h_dedx_slice_2_other_mixed,
	        h_dedx_slice_2_unmatched,
	        h_dedx_slice_2_intime,
	        intime_scale_factor,
	        h_dedx_slice_2_data,
	        data_scale_factor,
	        "Theta Slice (60 - 120)", "Leading Shower dE/dx [MeV/cm]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_2_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_dedx_slice_2_zoom_nue_cc,
	        h_dedx_slice_2_zoom_nue_cc_mixed,
	        h_dedx_slice_2_zoom_nue_cc_out_fv,
	        h_dedx_slice_2_zoom_numu_cc,
	        h_dedx_slice_2_zoom_numu_cc_mixed,
	        h_dedx_slice_2_zoom_cosmic,
	        h_dedx_slice_2_zoom_nc,
	        h_dedx_slice_2_zoom_nc_pi0,
	        h_dedx_slice_2_zoom_other_mixed,
	        h_dedx_slice_2_zoom_unmatched,
	        h_dedx_slice_2_zoom_intime,
	        intime_scale_factor,
	        h_dedx_slice_2_zoom_data,
	        data_scale_factor,
	        "Theta Slice (60 - 120)", "Leading Shower dE/dx [MeV/cm]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_2_zoom_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_dedx_slice_3_nue_cc,
	        h_dedx_slice_3_nue_cc_mixed,
	        h_dedx_slice_3_nue_cc_out_fv,
	        h_dedx_slice_3_numu_cc,
	        h_dedx_slice_3_numu_cc_mixed,
	        h_dedx_slice_3_cosmic,
	        h_dedx_slice_3_nc,
	        h_dedx_slice_3_nc_pi0,
	        h_dedx_slice_3_other_mixed,
	        h_dedx_slice_3_unmatched,
	        h_dedx_slice_3_intime,
	        intime_scale_factor,
	        h_dedx_slice_3_data,
	        data_scale_factor,
	        "Theta Slice (120 - 180)", "Leading Shower dE/dx [MeV/cm]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_3_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_dedx_slice_3_zoom_nue_cc,
	        h_dedx_slice_3_zoom_nue_cc_mixed,
	        h_dedx_slice_3_zoom_nue_cc_out_fv,
	        h_dedx_slice_3_zoom_numu_cc,
	        h_dedx_slice_3_zoom_numu_cc_mixed,
	        h_dedx_slice_3_zoom_cosmic,
	        h_dedx_slice_3_zoom_nc,
	        h_dedx_slice_3_zoom_nc_pi0,
	        h_dedx_slice_3_zoom_other_mixed,
	        h_dedx_slice_3_zoom_unmatched,
	        h_dedx_slice_3_zoom_intime,
	        intime_scale_factor,
	        h_dedx_slice_3_zoom_data,
	        data_scale_factor,
	        "Theta Slice (120 - 180)", "Leading Shower dE/dx [MeV/cm]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_3_zoom_data.pdf"));

	if(use_alt_scaling) {
		histogram_functions::PlotSimpleStackData(
		        h_dedx_slice_1_nue_cc,
		        h_dedx_slice_1_nue_cc_mixed,
		        h_dedx_slice_1_nue_cc_out_fv,
		        h_dedx_slice_1_numu_cc,
		        h_dedx_slice_1_numu_cc_mixed,
		        h_dedx_slice_1_cosmic,
		        h_dedx_slice_1_nc,
		        h_dedx_slice_1_nc_pi0,
		        h_dedx_slice_1_other_mixed,
		        h_dedx_slice_1_unmatched,
		        h_dedx_slice_1_intime,
		        scaled_intime_scale_factor,
		        h_dedx_slice_1_data,
		        data_scale_factor,
		        "Theta Slice (0 - 60)", "Leading Shower dE/dx [MeV/cm] (Scaled)", "",
		        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_1_scaled_data.pdf"));

		histogram_functions::PlotSimpleStackData(
		        h_dedx_slice_1_zoom_nue_cc,
		        h_dedx_slice_1_zoom_nue_cc_mixed,
		        h_dedx_slice_1_zoom_nue_cc_out_fv,
		        h_dedx_slice_1_zoom_numu_cc,
		        h_dedx_slice_1_zoom_numu_cc_mixed,
		        h_dedx_slice_1_zoom_cosmic,
		        h_dedx_slice_1_zoom_nc,
		        h_dedx_slice_1_zoom_nc_pi0,
		        h_dedx_slice_1_zoom_other_mixed,
		        h_dedx_slice_1_zoom_unmatched,
		        h_dedx_slice_1_zoom_intime,
		        scaled_intime_scale_factor,
		        h_dedx_slice_1_zoom_data,
		        data_scale_factor,
		        "Theta Slice (0 - 60)", "Leading Shower dE/dx [MeV/cm] (Scaled)", "",
		        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_1_zoom_scaled_data.pdf"));

		histogram_functions::PlotSimpleStackData(
		        h_dedx_slice_2_nue_cc,
		        h_dedx_slice_2_nue_cc_mixed,
		        h_dedx_slice_2_nue_cc_out_fv,
		        h_dedx_slice_2_numu_cc,
		        h_dedx_slice_2_numu_cc_mixed,
		        h_dedx_slice_2_cosmic,
		        h_dedx_slice_2_nc,
		        h_dedx_slice_2_nc_pi0,
		        h_dedx_slice_2_other_mixed,
		        h_dedx_slice_2_unmatched,
		        h_dedx_slice_2_intime,
		        scaled_intime_scale_factor,
		        h_dedx_slice_2_data,
		        data_scale_factor,
		        "Theta Slice (60 - 120)", "Leading Shower dE/dx [MeV/cm] (Scaled)", "",
		        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_2_scaled_data.pdf"));

		histogram_functions::PlotSimpleStackData(
		        h_dedx_slice_2_zoom_nue_cc,
		        h_dedx_slice_2_zoom_nue_cc_mixed,
		        h_dedx_slice_2_zoom_nue_cc_out_fv,
		        h_dedx_slice_2_zoom_numu_cc,
		        h_dedx_slice_2_zoom_numu_cc_mixed,
		        h_dedx_slice_2_zoom_cosmic,
		        h_dedx_slice_2_zoom_nc,
		        h_dedx_slice_2_zoom_nc_pi0,
		        h_dedx_slice_2_zoom_other_mixed,
		        h_dedx_slice_2_zoom_unmatched,
		        h_dedx_slice_2_zoom_intime,
		        scaled_intime_scale_factor,
		        h_dedx_slice_2_zoom_data,
		        data_scale_factor,
		        "Theta Slice (60 - 120)", "Leading Shower dE/dx [MeV/cm] (Scaled)", "",
		        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_2_zoom_scaled_data.pdf"));

		histogram_functions::PlotSimpleStackData(
		        h_dedx_slice_3_nue_cc,
		        h_dedx_slice_3_nue_cc_mixed,
		        h_dedx_slice_3_nue_cc_out_fv,
		        h_dedx_slice_3_numu_cc,
		        h_dedx_slice_3_numu_cc_mixed,
		        h_dedx_slice_3_cosmic,
		        h_dedx_slice_3_nc,
		        h_dedx_slice_3_nc_pi0,
		        h_dedx_slice_3_other_mixed,
		        h_dedx_slice_3_unmatched,
		        h_dedx_slice_3_intime,
		        scaled_intime_scale_factor,
		        h_dedx_slice_3_data,
		        data_scale_factor,
		        "Theta Slice (120 - 180)", "Leading Shower dE/dx [MeV/cm] (Scaled)", "",
		        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_3_scaled_data.pdf"));

		histogram_functions::PlotSimpleStackData(
		        h_dedx_slice_3_zoom_nue_cc,
		        h_dedx_slice_3_zoom_nue_cc_mixed,
		        h_dedx_slice_3_zoom_nue_cc_out_fv,
		        h_dedx_slice_3_zoom_numu_cc,
		        h_dedx_slice_3_zoom_numu_cc_mixed,
		        h_dedx_slice_3_zoom_cosmic,
		        h_dedx_slice_3_zoom_nc,
		        h_dedx_slice_3_zoom_nc_pi0,
		        h_dedx_slice_3_zoom_other_mixed,
		        h_dedx_slice_3_zoom_unmatched,
		        h_dedx_slice_3_zoom_intime,
		        scaled_intime_scale_factor,
		        h_dedx_slice_3_zoom_data,
		        data_scale_factor,
		        "Theta Slice (120 - 180)", "Leading Shower dE/dx [MeV/cm] (Scaled)", "",
		        Form("%s%s", file_locate_prefix, "post_cuts_dedx_theta_slice_3_zoom_scaled_data.pdf"));
	}

	histogram_functions::Plot1DHistogram(h_ele_resolution_momentum, "Momentum Resolution (True - Reco) / True",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_resolution_momentum.pdf"));
	histogram_functions::Plot1DHistogram(h_ele_resolution_phi, "Phi Resolution (True - Reco) / True",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_resolution_phi.pdf"));
	histogram_functions::Plot1DHistogram(h_ele_resolution_theta, "Theta Resolution (True - Reco) / True",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_resolution_theta.pdf"));
	histogram_functions::Plot1DHistogram(h_ele_resolution_dot_prod, "True.Reco Direction",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_resolution_dot_prod.pdf"));
	histogram_functions::Plot2DHistogram(h_ele_resolution_momentum_dot_prod, " ", "True Electron Momentum [GeV]", "True.Reco Shower Direction",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_resolution_momentum_dot_prod.pdf"));
	histogram_functions::Plot2DHistogram(h_ele_resolution_momentum_dot_prod_zoom_y, " ", "True Electron Momentum [GeV]", "True.Reco Shower Direction",
	                                     Form("%s%s", file_locate_prefix, "post_cuts_resolution_momentum_dot_prod_zoom_y.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_pfp_momentum_1shwr_nue_cc,
	        h_ele_pfp_momentum_1shwr_nue_cc_mixed,
	        h_ele_pfp_momentum_1shwr_nue_cc_out_fv,
	        h_ele_pfp_momentum_1shwr_numu_cc,
	        h_ele_pfp_momentum_1shwr_numu_cc_mixed,
	        h_ele_pfp_momentum_1shwr_cosmic,
	        h_ele_pfp_momentum_1shwr_nc,
	        h_ele_pfp_momentum_1shwr_nc_pi0,
	        h_ele_pfp_momentum_1shwr_other_mixed,
	        h_ele_pfp_momentum_1shwr_unmatched,
	        h_ele_pfp_momentum_1shwr_intime,
	        intime_scale_factor,
	        h_ele_pfp_momentum_1shwr_data,
	        data_scale_factor,
	        "", "Leading Shower Momentum (TPCO w/ 1 Shower) [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_pfp_momentum_1shwr_data.pdf"));

	histogram_functions::PlotSimpleStackDataMomentumRebin(
	        h_ele_pfp_momentum_2shwr_nue_cc,
	        h_ele_pfp_momentum_2shwr_nue_cc_mixed,
	        h_ele_pfp_momentum_2shwr_nue_cc_out_fv,
	        h_ele_pfp_momentum_2shwr_numu_cc,
	        h_ele_pfp_momentum_2shwr_numu_cc_mixed,
	        h_ele_pfp_momentum_2shwr_cosmic,
	        h_ele_pfp_momentum_2shwr_nc,
	        h_ele_pfp_momentum_2shwr_nc_pi0,
	        h_ele_pfp_momentum_2shwr_other_mixed,
	        h_ele_pfp_momentum_2shwr_unmatched,
	        h_ele_pfp_momentum_2shwr_intime,
	        intime_scale_factor,
	        h_ele_pfp_momentum_2shwr_data,
	        data_scale_factor,
	        "", "Leading Shower Momentum (TPCO w/ 2+ Showers) [GeV]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_pfp_momentum_2shwr_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_ele_pfp_theta_1shwr_nue_cc,
	        h_ele_pfp_theta_1shwr_nue_cc_mixed,
	        h_ele_pfp_theta_1shwr_nue_cc_out_fv,
	        h_ele_pfp_theta_1shwr_numu_cc,
	        h_ele_pfp_theta_1shwr_numu_cc_mixed,
	        h_ele_pfp_theta_1shwr_cosmic,
	        h_ele_pfp_theta_1shwr_nc,
	        h_ele_pfp_theta_1shwr_nc_pi0,
	        h_ele_pfp_theta_1shwr_other_mixed,
	        h_ele_pfp_theta_1shwr_unmatched,
	        h_ele_pfp_theta_1shwr_intime,
	        intime_scale_factor,
	        h_ele_pfp_theta_1shwr_data,
	        data_scale_factor,
	        "", "Leading Shower Theta (TPCO w/ 1 Shower) [Degrees]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_pfp_theta_1shwr_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_ele_pfp_theta_2shwr_nue_cc,
	        h_ele_pfp_theta_2shwr_nue_cc_mixed,
	        h_ele_pfp_theta_2shwr_nue_cc_out_fv,
	        h_ele_pfp_theta_2shwr_numu_cc,
	        h_ele_pfp_theta_2shwr_numu_cc_mixed,
	        h_ele_pfp_theta_2shwr_cosmic,
	        h_ele_pfp_theta_2shwr_nc,
	        h_ele_pfp_theta_2shwr_nc_pi0,
	        h_ele_pfp_theta_2shwr_other_mixed,
	        h_ele_pfp_theta_2shwr_unmatched,
	        h_ele_pfp_theta_2shwr_intime,
	        intime_scale_factor,
	        h_ele_pfp_theta_2shwr_data,
	        data_scale_factor,
	        "", "Leading Shower Theta (TPCO w/ 2+ Showers) [Degrees]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_pfp_theta_2shwr_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_ele_pfp_phi_1shwr_nue_cc,
	        h_ele_pfp_phi_1shwr_nue_cc_mixed,
	        h_ele_pfp_phi_1shwr_nue_cc_out_fv,
	        h_ele_pfp_phi_1shwr_numu_cc,
	        h_ele_pfp_phi_1shwr_numu_cc_mixed,
	        h_ele_pfp_phi_1shwr_cosmic,
	        h_ele_pfp_phi_1shwr_nc,
	        h_ele_pfp_phi_1shwr_nc_pi0,
	        h_ele_pfp_phi_1shwr_other_mixed,
	        h_ele_pfp_phi_1shwr_unmatched,
	        h_ele_pfp_phi_1shwr_intime,
	        intime_scale_factor,
	        h_ele_pfp_phi_1shwr_data,
	        data_scale_factor,
	        "", "Leading Shower Phi (TPCO w/ 1 Shower) [Degrees]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_pfp_phi_1shwr_data.pdf"));

	histogram_functions::PlotSimpleStackData(
	        h_ele_pfp_phi_2shwr_nue_cc,
	        h_ele_pfp_phi_2shwr_nue_cc_mixed,
	        h_ele_pfp_phi_2shwr_nue_cc_out_fv,
	        h_ele_pfp_phi_2shwr_numu_cc,
	        h_ele_pfp_phi_2shwr_numu_cc_mixed,
	        h_ele_pfp_phi_2shwr_cosmic,
	        h_ele_pfp_phi_2shwr_nc,
	        h_ele_pfp_phi_2shwr_nc_pi0,
	        h_ele_pfp_phi_2shwr_other_mixed,
	        h_ele_pfp_phi_2shwr_unmatched,
	        h_ele_pfp_phi_2shwr_intime,
	        intime_scale_factor,
	        h_ele_pfp_phi_2shwr_data,
	        data_scale_factor,
	        "", "Leading Shower Phi (TPCO w/ 2+ Showers) [Degrees]", "",
	        Form("%s%s", file_locate_prefix, "post_cuts_pfp_phi_2shwr_data.pdf"));

	TCanvas * failure_reason_stack_c1 = new TCanvas();
	failure_reason_stack_c1->cd();
	THStack * h_failure_reason_stack = new THStack();
	h_failure_reason_nue_cc->SetStats(kFALSE);
	h_failure_reason_nue_cc_mixed->SetStats(kFALSE);
	h_failure_reason_numu_cc->SetStats(kFALSE);
	h_failure_reason_nc_pi0->SetStats(kFALSE);
	h_failure_reason_cosmic->SetStats(kFALSE);
	h_failure_reason_nc->SetStats(kFALSE);
	h_failure_reason_numu_cc_mixed->SetStats(kFALSE);
	h_failure_reason_other_mixed->SetStats(kFALSE);
	h_failure_reason_unmatched->SetStats(kFALSE);
	h_failure_reason_nue_cc->SetFillColor(30);
	h_failure_reason_nue_cc_mixed->SetFillColor(38);
	h_failure_reason_numu_cc->SetFillColor(28);
	h_failure_reason_nc_pi0->SetFillColor(36);
	h_failure_reason_cosmic->SetFillColor(1);
	h_failure_reason_nc->SetFillColor(46);
	h_failure_reason_numu_cc_mixed->SetFillColor(20);
	h_failure_reason_other_mixed->SetFillColor(42);
	h_failure_reason_unmatched->SetFillColor(12);
	h_failure_reason_stack->Draw();

	std::vector<double> remaining_nue_cc;
	std::vector<double> remaining_nue_cc_mixed;
	std::vector<double> remaining_numu_cc;
	std::vector<double> remaining_numu_cc_mixed;
	std::vector<double> remaining_cosmic;
	std::vector<double> remaining_nc;
	std::vector<double> remaining_nc_pi0;
	std::vector<double> remaining_other_mixed;
	std::vector<double> remaining_unmatched;

	double bin_content_nue_cc_prev = pe_counter_v->at(0);
	double bin_content_nue_cc_mixed_prev = pe_counter_v->at(1);
	double bin_content_numu_cc_prev = pe_counter_v->at(4);
	double bin_content_nc_pi0_prev = pe_counter_v->at(10);
	double bin_content_cosmic_prev = pe_counter_v->at(2);
	double bin_content_nc_prev = pe_counter_v->at(3);
	double bin_content_numu_cc_mixed_prev = pe_counter_v->at(11);
	double bin_content_other_mixed_prev = pe_counter_v->at(6);
	double bin_content_unmatched_prev = pe_counter_v->at(5);

	for (int i=1; i<= 24; i++)
	{
		double bin_content_nue_cc                   = h_failure_reason_nue_cc->GetBinContent(i);
		double bin_content_nue_cc_mixed             = h_failure_reason_nue_cc_mixed->GetBinContent(i);
		double bin_content_numu_cc                  = h_failure_reason_numu_cc->GetBinContent(i);
		double bin_content_nc_pi0                   = h_failure_reason_nc_pi0->GetBinContent(i);
		double bin_content_cosmic                   = h_failure_reason_cosmic->GetBinContent(i);
		double bin_content_nc                       = h_failure_reason_nc->GetBinContent(i);
		double bin_content_numu_cc_mixed            = h_failure_reason_numu_cc_mixed->GetBinContent(i);
		double bin_content_other_mixed              = h_failure_reason_other_mixed->GetBinContent(i);
		double bin_content_unmatched                = h_failure_reason_unmatched->GetBinContent(i);
		if(i % 2 != 0)
		{
			const double ratio_nue_cc        = double(bin_content_nue_cc) / double(bin_content_nue_cc_prev);
			const double ratio_nue_cc_mixed  = double(bin_content_nue_cc_mixed) / double(bin_content_nue_cc_mixed_prev);
			const double ratio_numu_cc       = double(bin_content_numu_cc) / double(bin_content_numu_cc_prev);
			const double ratio_numu_cc_mixed = double(bin_content_numu_cc_mixed) / double(bin_content_numu_cc_mixed_prev);
			const double ratio_cosmic        = double(bin_content_cosmic) / double(bin_content_cosmic_prev);
			const double ratio_nc            = double(bin_content_nc) / double(bin_content_nc_prev);
			const double ratio_nc_pi0        = double(bin_content_nc_pi0) / double(bin_content_nc_pi0_prev);
			const double ratio_other_mixed   = double(bin_content_other_mixed) / double(bin_content_other_mixed_prev);
			const double ratio_unmatched     = double(bin_content_unmatched) / double(bin_content_unmatched_prev);

			bin_content_nue_cc_prev              -= bin_content_nue_cc;
			bin_content_nue_cc_mixed_prev        -= bin_content_nue_cc_mixed;
			bin_content_numu_cc_prev             -= bin_content_numu_cc;
			bin_content_nc_pi0_prev              -= bin_content_nc_pi0;
			bin_content_cosmic_prev              -= bin_content_cosmic;
			bin_content_nc_prev                  -= bin_content_nc;
			bin_content_numu_cc_mixed_prev       -= bin_content_numu_cc_mixed;
			bin_content_other_mixed_prev         -= bin_content_other_mixed;
			bin_content_unmatched_prev           -= bin_content_unmatched;

			remaining_nue_cc.push_back(ratio_nue_cc);
			remaining_nue_cc_mixed.push_back(ratio_nue_cc_mixed);
			remaining_numu_cc.push_back(ratio_numu_cc);
			remaining_numu_cc_mixed.push_back(ratio_numu_cc_mixed);
			remaining_cosmic.push_back(ratio_cosmic);
			remaining_nc.push_back(ratio_nc);
			remaining_nc_pi0.push_back(ratio_nc_pi0);
			remaining_other_mixed.push_back(ratio_other_mixed);
			remaining_unmatched.push_back(ratio_unmatched);
		}
	}
	for (int i = 0; i < remaining_nue_cc.size(); i++)
	{
		h_failure_reason_nue_cc->SetBinContent((i*2)+1,         remaining_nue_cc.at(i));
		h_failure_reason_nue_cc_mixed->SetBinContent((i*2)+1,   remaining_nue_cc_mixed.at(i));
		h_failure_reason_numu_cc->SetBinContent((i*2)+1,        remaining_numu_cc.at(i));
		h_failure_reason_nc_pi0->SetBinContent((i*2)+1,         remaining_nc_pi0.at(i));
		h_failure_reason_cosmic->SetBinContent((i*2)+1,         remaining_cosmic.at(i));
		h_failure_reason_nc->SetBinContent((i*2)+1,             remaining_nc.at(i));
		h_failure_reason_numu_cc_mixed->SetBinContent((i*2)+1,  remaining_numu_cc_mixed.at(i));
		h_failure_reason_other_mixed->SetBinContent((i*2)+1,    remaining_other_mixed.at(i));
		h_failure_reason_unmatched->SetBinContent((i*2)+1,      remaining_unmatched.at(i));
	}

	h_failure_reason_stack->Add(h_failure_reason_nue_cc);
	h_failure_reason_stack->Add(h_failure_reason_nue_cc_mixed);
	h_failure_reason_stack->Add(h_failure_reason_cosmic);
	h_failure_reason_stack->Add(h_failure_reason_numu_cc);
	h_failure_reason_stack->Add(h_failure_reason_numu_cc_mixed);
	h_failure_reason_stack->Add(h_failure_reason_nc);
	h_failure_reason_stack->Add(h_failure_reason_nc_pi0);
	h_failure_reason_stack->Add(h_failure_reason_other_mixed);
	h_failure_reason_stack->Add(h_failure_reason_unmatched);
	h_failure_reason_stack->Draw();
	h_failure_reason_stack->GetXaxis()->SetTitle("Failure Reason");
	h_failure_reason_stack->GetYaxis()->SetTitle("Rejected / (Rejected + Selected)");
	const char * str_cut[24] = {"HasNue", " ", "InFV", " ", "FlshDist", " ", "ShwrVtx", " ", "TrkVtx", " ", "Hits", " ",
		                    "OpenAngle", " ", "dEdx", " ", "2ndDist", " ", "HitLength", " ", "YHits", " ", "T/S-Ratio", " "};
	for (int i=1; i<= 24; i++)
	{
		h_failure_reason_stack->GetXaxis()->SetBinLabel(i,str_cut[(i-1)]);
	}

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack_failure_reason = new TLegend(0.865,0.75,0.99,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack_failure_reason->AddEntry(h_failure_reason_nue_cc,          "Nue CC", "f");
	leg_stack_failure_reason->AddEntry(h_failure_reason_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack_failure_reason->AddEntry(h_failure_reason_cosmic,          "Cosmic", "f");
	leg_stack_failure_reason->AddEntry(h_failure_reason_numu_cc,         "Numu CC", "f");
	leg_stack_failure_reason->AddEntry(h_failure_reason_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack_failure_reason->AddEntry(h_failure_reason_nc,              "NC", "f");
	leg_stack_failure_reason->AddEntry(h_failure_reason_nc_pi0,          "NC Pi0", "f");
	leg_stack_failure_reason->AddEntry(h_failure_reason_other_mixed,     "Other Mixed", "f");
	leg_stack_failure_reason->AddEntry(h_failure_reason_unmatched,       "Unmatched", "f");
	leg_stack_failure_reason->Draw();
	failure_reason_stack_c1->Print(Form("%s%s", file_locate_prefix, "failure_reason_stack.pdf"));


	std::cout << " --- End Cross Section Calculation --- " << std::endl;

	if(f->IsOpen()) {f->Close(); }
	//for some reason, the histogram is not being deleted upon function exit
	//trying to close root file instead
	//delete h_nue_eng_eff_den;
	//delete h_nue_eng_eff_num;

}//end selection
}//end namespace

#ifndef __ROOTCLING__

#endif
