#include "selection_slim.h"

namespace xsecSelection {
void selection_slim::make_selection_slim( const char * _file1,
                                          const char * _file2,
                                          const char * _file3,
                                          const std::vector<double> _config,
                                          std::vector<std::tuple<double, double, std::string> > * results_v
                                          )
{

	std::cout << "File Path: " << _file1 << std::endl;
	const bool _verbose = false;
	const bool _post_cuts_verbose = false;
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
	const double flux = POT * scaling;

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
	int total_mc_entries_inFV = 0;
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
		if(true_in_tpc == true && (mc_nu_id == 1 || mc_nu_id == 5)) {total_mc_entries_inFV++; }
	}

	std::cout << "MC Nue CC Counter      --- " << _mc_nue_cc_counter << std::endl;
	std::cout << "MC Nue NC Counter      --- " << _mc_nue_nc_counter << std::endl;
	std::cout << "MC Numu CC Counter     --- " << _mc_numu_cc_counter << std::endl;
	std::cout << "MC Numu NC Counter     --- " << _mc_numu_nc_counter << std::endl;
	std::cout << "MC Nue CC Counter Bar  --- " << _mc_nue_cc_counter_bar << std::endl;
	std::cout << "MC Nue NC Counter Bar  --- " << _mc_nue_nc_counter_bar << std::endl;
	std::cout << "MC Numu CC Counter Bar --- " << _mc_numu_cc_counter_bar << std::endl;
	std::cout << "MC Numu NC Counter Bar --- " << _mc_numu_nc_counter_bar << std::endl;

	//*****************
	//***** DATA *****
	//*****************

	std::vector<int> * data_in_time_counter_v = new std::vector<int>;
	data_in_time_counter_v->resize(22, 0);
	std::vector<int> * data_pe_counter_v = new std::vector<int>;
	data_pe_counter_v->resize(22, 0);
	std::vector<int> * data_reco_nue_counter_v = new std::vector<int>;
	data_reco_nue_counter_v->resize(22, 0);
	std::vector<int> * data_in_fv_counter_v = new std::vector<int>;
	data_in_fv_counter_v->resize(22, 0);
	std::vector<int> * data_vtx_flash_counter_v = new std::vector<int>;
	data_vtx_flash_counter_v->resize(22, 0);
	std::vector<int> * data_shwr_tpco_counter_v = new std::vector<int>;
	data_shwr_tpco_counter_v->resize(22, 0);
	std::vector<int> * data_trk_tpco_counter_v = new std::vector<int>;
	data_trk_tpco_counter_v->resize(22, 0);
	std::vector<int> * data_hit_threshold_counter_v = new std::vector<int>;
	data_hit_threshold_counter_v->resize(22, 0);
	std::vector<int> * data_open_angle_counter_v = new std::vector<int>;
	data_open_angle_counter_v->resize(22, 0);
	std::vector<int> * data_dedx_counter_v = new std::vector<int>;
	data_dedx_counter_v->resize(22, 0);
	std::vector<int> * data_secondary_shower_counter_v = new std::vector<int>;
	data_secondary_shower_counter_v->resize(22, 0);
	std::vector<int> * data_hit_lengthRatio_counter_v = new std::vector<int>;
	data_hit_lengthRatio_counter_v->resize(22, 0);
	std::vector<int> * data_hit_threshold_collection_counter_v = new std::vector<int>;
	data_hit_threshold_collection_counter_v->resize(22, 0);
	std::vector<int> * data_trk_len_shwr_len_ratio_counter_v = new std::vector<int>;
	data_trk_len_shwr_len_ratio_counter_v->resize(22, 0);
	std::vector<int> * data_track_containment_counter_v = new std::vector<int>;
	data_track_containment_counter_v->resize(22, 0);

	std::vector<double> * data_flash_time = new std::vector<double>;

	std::vector<int> * tabulated_origins_data = new std::vector<int>;
	tabulated_origins_data->resize(22, 0);

	std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v_data
	        = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> >;

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

		_cuts_instance.selection_cuts::loop_flashes(data_f, data_optree, flash_pe_threshold, flash_time_start, flash_time_end, data_passed_runs, data_flash_time);
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
			}        //false

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

			//XY Position of largest flash
			std::vector < double > largest_flash_v = data_largest_flash_v_v->at(event);

			//List of TPC Objects which pass the cuts
			std::vector<std::pair<int, std::string> > * passed_tpco_data = new std::vector<std::pair<int, std::string> >;
			passed_tpco_data->resize(data_tpc_object_container_v->size());
			//set initial state of objects
			for(int i = 0; i < passed_tpco_data->size(); i++)
			{
				passed_tpco_data->at(i).first = 1;
				passed_tpco_data->at(i).second = "Passed";
			}
			//***********************************************************
			//this is where the in-time optical cut again takes effect
			//***********************************************************
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_in_time_counter_v);

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

			//************************
			//******** in fv cut *****
			//************************

			_cuts_instance.selection_cuts::fiducial_volume_cut(data_tpc_object_container_v, fv_boundary_v, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_in_fv_counter_v);

			//*****************************
			//**** vertex to flash cut ****
			//*****************************

			_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, data_tpc_object_container_v, tolerance, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_vtx_flash_counter_v);

			//******************************************************
			//*** distance between pfp shower and nue object cut ***
			//******************************************************

			_cuts_instance.selection_cuts::VtxNuDistance(data_tpc_object_container_v, shwr_nue_tolerance, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_shwr_tpco_counter_v);

			//******************************************************
			// **** distance between pfp track and nue object cut **
			//******************************************************

			_cuts_instance.selection_cuts::VtxTrackNuDistance(data_tpc_object_container_v, trk_nue_tolerance, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_trk_tpco_counter_v);

			//****************************************************
			// ******** hit threshold for showers cut *************
			//******************************************************

			_cuts_instance.selection_cuts::HitThreshold(data_tpc_object_container_v, shwr_hit_threshold, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_hit_threshold_counter_v);

			//***************************************//
			//*** Collection Plane Hits Threshold ***//
			//***************************************//

			_cuts_instance.selection_cuts::HitThresholdCollection(data_tpc_object_container_v, shwr_hit_threshold_collection, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_hit_threshold_collection_counter_v);

			//*****************************************************
			//****** open angle cut for the leading shower ********
			//******************************************************

			_cuts_instance.selection_cuts::OpenAngleCut(data_tpc_object_container_v, passed_tpco_data, tolerance_open_angle, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_open_angle_counter_v);

			//*****************************************************
			//*********** dEdx cut for the leading shower *********
			//******************************************************

			_cuts_instance.selection_cuts::dEdxCut(data_tpc_object_container_v, passed_tpco_data, tolerance_dedx_min, tolerance_dedx_max, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_dedx_counter_v);

			//***************************************************************************
			// ******* Secondary Showers Distance Cut *****************
			//***************************************************************************

			_cuts_instance.selection_cuts::SecondaryShowersDistCut(data_tpc_object_container_v, passed_tpco_data, _verbose, dist_tolerance);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_secondary_shower_counter_v);

			//******************************************************************************
			// ********** Hit Length Ratio Cut *************
			//******************************************************************************

			_cuts_instance.selection_cuts::HitLengthRatioCut(data_tpc_object_container_v, passed_tpco_data, _verbose, pfp_hits_length_tolerance);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_hit_lengthRatio_counter_v);

			//******************************************************************************
			//*** cut for longest track / leading shower ratio *** //
			//******************************************************************************

			_cuts_instance.selection_cuts::LongestTrackLeadingShowerCut(data_tpc_object_container_v, passed_tpco_data, _verbose, ratio_tolerance);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_trk_len_shwr_len_ratio_counter_v);

			//***************************************************************
			//*** contained track cut *** //
			//**************************************************************

			_cuts_instance.selection_cuts::ContainedTracksCut(data_tpc_object_container_v, passed_tpco_data, _verbose, fv_boundary_v, true);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_track_containment_counter_v);

			//*********** Data *********
			_functions_instance.selection_functions::FillPostCutVector(data_tpc_object_container_v, passed_tpco_data, post_cuts_v_data);
			//*************************************
			// ******** End Selection Cuts! *******
			//*************************************

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
	intime_in_time_counter_v->resize(22, 0);
	std::vector<int> * intime_pe_counter_v = new std::vector<int>;
	intime_pe_counter_v->resize(22, 0);
	std::vector<int> * intime_reco_nue_counter_v = new std::vector<int>;
	intime_reco_nue_counter_v->resize(22, 0);
	std::vector<int> * intime_in_fv_counter_v = new std::vector<int>;
	intime_in_fv_counter_v->resize(22, 0);
	std::vector<int> * intime_vtx_flash_counter_v = new std::vector<int>;
	intime_vtx_flash_counter_v->resize(22, 0);
	std::vector<int> * intime_shwr_tpco_counter_v = new std::vector<int>;
	intime_shwr_tpco_counter_v->resize(22, 0);
	std::vector<int> * intime_trk_tpco_counter_v = new std::vector<int>;
	intime_trk_tpco_counter_v->resize(22, 0);
	std::vector<int> * intime_hit_threshold_counter_v = new std::vector<int>;
	intime_hit_threshold_counter_v->resize(22, 0);
	std::vector<int> * intime_open_angle_counter_v = new std::vector<int>;
	intime_open_angle_counter_v->resize(22, 0);
	std::vector<int> * intime_dedx_counter_v = new std::vector<int>;
	intime_dedx_counter_v->resize(22, 0);
	std::vector<int> * intime_secondary_shower_counter_v = new std::vector<int>;
	intime_secondary_shower_counter_v->resize(22, 0);
	std::vector<int> * intime_hit_lengthRatio_counter_v = new std::vector<int>;
	intime_hit_lengthRatio_counter_v->resize(22, 0);
	std::vector<int> * intime_hit_threshold_collection_counter_v = new std::vector<int>;
	intime_hit_threshold_collection_counter_v->resize(22, 0);
	std::vector<int> * intime_trk_len_shwr_len_ratio_counter_v = new std::vector<int>;
	intime_trk_len_shwr_len_ratio_counter_v->resize(22, 0);
	std::vector<int> * intime_track_containment_counter_v = new std::vector<int>;
	intime_track_containment_counter_v->resize(22, 0);

	std::vector<double> * intime_flash_time = new std::vector<double>;

	std::vector<int> * tabulated_origins = new std::vector<int>;
	tabulated_origins->resize(22, 0);
	std::vector<int> * tabulated_origins_intime = new std::vector<int>;
	tabulated_origins_intime->resize(22, 0);

	std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v
	        = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> >;

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
		                                            flash_time_start, flash_time_end, intime_passed_runs, intime_flash_time);
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

			//XY Position of largest flash
			std::vector < double > largest_flash_v = intime_largest_flash_v_v->at(event);

			//List of TPC Objects which pass the cuts
			std::vector<std::pair<int, std::string> > * passed_tpco_intime = new std::vector<std::pair<int, std::string> >;
			passed_tpco_intime->resize(intime_tpc_object_container_v->size());
			//set initial state of objects
			for(int i = 0; i < passed_tpco_intime->size(); i++)
			{
				passed_tpco_intime->at(i).first = 1;
				passed_tpco_intime->at(i).second = "Passed";
			}
			//***********************************************************
			//this is where the in-time optical cut again takes effect
			//***********************************************************
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_in_time_counter_v);

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

			_cuts_instance.selection_cuts::HasNue(intime_tpc_object_container_v, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_reco_nue_counter_v);

			//************************
			//******** in fv cut *****
			//************************

			_cuts_instance.selection_cuts::fiducial_volume_cut(intime_tpc_object_container_v, fv_boundary_v, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_in_fv_counter_v);

			//*****************************
			//**** vertex to flash cut ****
			//*****************************

			_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, intime_tpc_object_container_v, tolerance, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_vtx_flash_counter_v);

			//******************************************************
			//*** distance between pfp shower and nue object cut ***
			//******************************************************

			_cuts_instance.selection_cuts::VtxNuDistance(intime_tpc_object_container_v, shwr_nue_tolerance, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_shwr_tpco_counter_v);

			//******************************************************
			// **** distance between pfp track and nue object cut **
			//******************************************************

			_cuts_instance.selection_cuts::VtxTrackNuDistance(intime_tpc_object_container_v, trk_nue_tolerance, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_trk_tpco_counter_v);

			//****************************************************
			// ******** hit threshold for showers cut *************
			//******************************************************

			_cuts_instance.selection_cuts::HitThreshold(intime_tpc_object_container_v, shwr_hit_threshold, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_hit_threshold_counter_v);

			//***************************************//
			//*** Collection Plane Hits Threshold ***//
			//***************************************//

			_cuts_instance.selection_cuts::HitThresholdCollection(intime_tpc_object_container_v, shwr_hit_threshold_collection, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_hit_threshold_collection_counter_v);

			//*****************************************************
			//****** open angle cut for the leading shower ********
			//******************************************************

			_cuts_instance.selection_cuts::OpenAngleCut(intime_tpc_object_container_v, passed_tpco_intime, tolerance_open_angle, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_open_angle_counter_v);

			//*****************************************************
			//*********** dEdx cut for the leading shower *********
			//******************************************************

			_cuts_instance.selection_cuts::dEdxCut(intime_tpc_object_container_v, passed_tpco_intime, tolerance_dedx_min, tolerance_dedx_max, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_dedx_counter_v);

			//***************************************************************************
			// ******* Secondary Showers Distance Cut *****************
			//***************************************************************************

			_cuts_instance.selection_cuts::SecondaryShowersDistCut(intime_tpc_object_container_v, passed_tpco_intime, _verbose, dist_tolerance);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_secondary_shower_counter_v);

			//******************************************************************************
			// ********** Hit Length Ratio Cut *************
			//******************************************************************************

			_cuts_instance.selection_cuts::HitLengthRatioCut(intime_tpc_object_container_v, passed_tpco_intime, _verbose, pfp_hits_length_tolerance);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_hit_lengthRatio_counter_v);

			//******************************************************************************
			//*** cut for longest track / leading shower ratio *** //
			//******************************************************************************

			_cuts_instance.selection_cuts::LongestTrackLeadingShowerCut(intime_tpc_object_container_v, passed_tpco_intime, _verbose, ratio_tolerance);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_trk_len_shwr_len_ratio_counter_v);

			//******************************************************************************
			//*** contained track cut *** //
			//******************************************************************************

			_cuts_instance.selection_cuts::ContainedTracksCut(intime_tpc_object_container_v, passed_tpco_intime, _verbose, fv_boundary_v, true);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_track_containment_counter_v);

			//*********** In-time Cosmics *********
			//*************************************
			// ******** End Selection Cuts! *******
			//*************************************

			_functions_instance.selection_functions::FillPostCutVector(intime_tpc_object_container_v, passed_tpco_intime, post_cuts_v);
			//delete at the very end!
			delete passed_tpco_intime;
		}
	}//end in-time cosmic running
//*********************************************************
//*********************************************************
	std::vector<int> * in_time_counter_v = new std::vector<int>;
	in_time_counter_v->resize(22, 0);
	std::vector<int> * pe_counter_v = new std::vector<int>;
	pe_counter_v->resize(22, 0);
	std::vector<int> * reco_nue_counter_v = new std::vector<int>;
	reco_nue_counter_v->resize(22, 0);
	std::vector<int> * in_fv_counter_v = new std::vector<int>;
	in_fv_counter_v->resize(22, 0);
	std::vector<int> * vtx_flash_counter_v = new std::vector<int>;
	vtx_flash_counter_v->resize(22, 0);
	std::vector<int> * shwr_tpco_counter_v = new std::vector<int>;
	shwr_tpco_counter_v->resize(22, 0);
	std::vector<int> * trk_tpco_counter_v = new std::vector<int>;
	trk_tpco_counter_v->resize(22, 0);
	std::vector<int> * hit_threshold_counter_v = new std::vector<int>;
	hit_threshold_counter_v->resize(22, 0);
	std::vector<int> * open_angle_counter_v = new std::vector<int>;
	open_angle_counter_v->resize(22, 0);
	std::vector<int> * dedx_counter_v = new std::vector<int>;
	dedx_counter_v->resize(22, 0);
	std::vector<int> * secondary_shower_counter_v = new std::vector<int>;
	secondary_shower_counter_v->resize(22, 0);
	std::vector<int> * hit_lengthRatio_counter_v = new std::vector<int>;
	hit_lengthRatio_counter_v->resize(22, 0);
	std::vector<int> * hit_threshold_collection_counter_v = new std::vector<int>;
	hit_threshold_collection_counter_v->resize(22, 0);
	std::vector<int> * trk_len_shwr_len_ratio_counter_v = new std::vector<int>;
	trk_len_shwr_len_ratio_counter_v->resize(22, 0);
	std::vector<int> * track_containment_counter_v = new std::vector<int>;
	track_containment_counter_v->resize(22, 0);

	std::vector<double> * flash_time = new std::vector<double>;

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

	//Event, Run, VtxX, VtxY, VtxZ, pass/fail reason
	// std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v
	//         = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> >;
	std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_open_angle_cuts_v
	        = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> >;

	std::cout << "=====================" << std::endl;
	std::cout << "== Begin Selection ==" << std::endl;
	std::cout << "=====================" << std::endl;

	const int total_entries = mytree->GetEntries();
	std::cout << "Total Events     : " << total_entries << std::endl;
	std::cout << "Total Events (MC): " << mctruth_counter_tree->GetEntries() << std::endl;

	std::vector<int> * passed_runs = new std::vector<int>;
	//passed runs is filled with 0, 1, or 2
	//0 = not in time
	//1 = passed - in time and PE threshold
	// 2 = in-time, but not enough PE -- this counts against my efficiency
	passed_runs->resize(total_entries);

	//let's do the in-time cut as the very first thing
	std::cout << "=====================" << std::endl;
	std::cout << "==== In Time Cut ====" << std::endl;
	std::cout << "=====================" << std::endl;

	_cuts_instance.selection_cuts::loop_flashes(f, optree, flash_pe_threshold, flash_time_start, flash_time_end, passed_runs, flash_time);
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

	//**********************************
	//now let's do the TPCO related cuts
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
		//***********************************************************
		//this is where the in-time optical cut actually takes effect
		//***********************************************************
		if(passed_runs->at(event) == 0)
		{
			if(_verbose) std::cout << "[Failed In-Time Cut]" << std::endl;
			continue;
		}//false

		//check if nue interaction has true vtx in TPC
		//std::cout << tpc_object_container_v->size() << std::endl;
		bool true_in_tpc = false;
		for (int i = 0; i < tpc_object_container_v->size(); i++)
		{
			auto const tpc_object_container = tpc_object_container_v->at(i);
			double _mc_nu_vtx_x = 0.;
			double _mc_nu_vtx_y = 0.;
			double _mc_nu_vtx_z = 0.;
			_mc_nu_vtx_x = tpc_object_container.mcVtxX();
			_mc_nu_vtx_y = tpc_object_container.mcVtxY();
			_mc_nu_vtx_z = tpc_object_container.mcVtxZ();

			true_in_tpc = _cuts_instance.selection_cuts::in_fv(_mc_nu_vtx_x, _mc_nu_vtx_y, _mc_nu_vtx_z, fv_boundary_v);
			if(true_in_tpc == true) {break; }
		}
		//const bool true_in_tpc = true_in_tpc_v.at(event);

		//now we apply the classifier to all TPC Objects in this event
		std::vector<std::pair<std::string, int> > * tpco_classifier_v = new std::vector<std::pair<std::string, int> >;
		_functions_instance.selection_functions::FillTPCOClassV(tpc_object_container_v, true_in_tpc, has_pi0, tpco_classifier_v);

		//XY Position of largest flash
		std::vector < double > largest_flash_v = largest_flash_v_v->at(event);

		//List of TPC Objects which pass the cuts
		std::vector<std::pair<int, std::string> > * passed_tpco = new std::vector<std::pair<int, std::string> >;
		passed_tpco->resize(tpc_object_container_v->size());
		for(int i = 0; i < passed_tpco->size(); i++)
		{
			passed_tpco->at(i).first = 1;
			passed_tpco->at(i).second = "Passed";
		}

		//***********************************************************
		//this is where the in-time optical cut again takes effect
		//***********************************************************
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, in_time_counter_v);

		//PE threshold cut
		if(passed_runs->at(event) == 2)
		{
			if(_verbose) std::cout << "[Passed In-Time Cut] [Failed PE Threshold] " << std::endl;
			continue;
		}
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, pe_counter_v);

		//****************************
		// ****** reco nue cut *******
		//****************************

		_cuts_instance.selection_cuts::HasNue(tpc_object_container_v, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, reco_nue_counter_v);

		//************************
		//******** in fv cut *****
		//************************
		_cuts_instance.selection_cuts::fiducial_volume_cut(tpc_object_container_v, fv_boundary_v, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, in_fv_counter_v);

		//*****************************
		//**** vertex to flash cut ****
		//*****************************

		_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, tpc_object_container_v, tolerance, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, vtx_flash_counter_v);

		//******************************************************
		//*** distance between pfp shower and nue object cut ***
		//******************************************************

		_cuts_instance.selection_cuts::VtxNuDistance(tpc_object_container_v, shwr_nue_tolerance, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, shwr_tpco_counter_v);

		//******************************************************
		// **** distance between pfp track and nue object cut **
		//******************************************************

		_cuts_instance.selection_cuts::VtxTrackNuDistance(tpc_object_container_v, trk_nue_tolerance, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, trk_tpco_counter_v);

		//****************************************************
		// ******** hit threshold for showers cut *************
		//******************************************************

		_cuts_instance.selection_cuts::HitThreshold(tpc_object_container_v, shwr_hit_threshold, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_threshold_counter_v);

		//***************************************//
		//*** Collection Plane Hits Threshold ***//
		//***************************************//

		_cuts_instance.selection_cuts::HitThresholdCollection(tpc_object_container_v, shwr_hit_threshold_collection, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_threshold_collection_counter_v);

		//*****************************************************
		//****** open angle cut for the leading shower ********
		//******************************************************

		_cuts_instance.selection_cuts::OpenAngleCut(tpc_object_container_v, passed_tpco, tolerance_open_angle, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, open_angle_counter_v);

		//*****************************************************
		//*********** dEdx cut for the leading shower *********
		//******************************************************

		_cuts_instance.selection_cuts::dEdxCut(tpc_object_container_v, passed_tpco, tolerance_dedx_min, tolerance_dedx_max, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, dedx_counter_v);

		//***************************************************************************
		// ******* Secondary Showers Distance Cut *****************
		//***************************************************************************

		_cuts_instance.selection_cuts::SecondaryShowersDistCut(tpc_object_container_v, passed_tpco, _verbose, dist_tolerance);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, secondary_shower_counter_v);

		//******************************************************************************
		// ********** Hit Length Ratio Cut *************
		//******************************************************************************

		_cuts_instance.selection_cuts::HitLengthRatioCut(tpc_object_container_v, passed_tpco, _verbose, pfp_hits_length_tolerance);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_lengthRatio_counter_v);

		//******************************************************************************
		//*** cut for longest track / leading shower ratio *** //
		//******************************************************************************

		_cuts_instance.selection_cuts::LongestTrackLeadingShowerCut(tpc_object_container_v, passed_tpco, _verbose, ratio_tolerance);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, trk_len_shwr_len_ratio_counter_v);

		//******************************************************************************
		//*** track containment -- all tracks need to be contained *** //
		//******************************************************************************

		_cuts_instance.selection_cuts::ContainedTracksCut(tpc_object_container_v, passed_tpco, _verbose, fv_boundary_v, true);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, track_containment_counter_v);

		//*************************************
		// ******** End Selection Cuts! ******
		//*************************************

		//these are for the tefficiency plots, post all cuts
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			int pass = 0;
			for(auto const passed_tpco_pair : * passed_tpco) {pass += passed_tpco_pair.first; }
			if(pass > 0)
			{
				selected_energy_vector.push_back(mc_nu_energy);
			}
		}

		_functions_instance.selection_functions::TopologyEfficiency(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                            no_track, has_track, _1_shwr, _2_shwr, _3_shwr, _4_shwr);

		_functions_instance.selection_functions::FillPostCutVector(tpc_object_container_v, passed_tpco, tpco_classifier_v, post_cuts_v);

	}//end event loop

	std::cout << "------------------ " << std::endl;
	std::cout << " MC Entries in FV: " << total_mc_entries_inFV << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "  End Selection    " << std::endl;
	std::cout << "------------------ " << std::endl;

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
	                                            intime_scale_factor, data_scale_factor, ">1 Shower TPCO Dist", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, hit_lengthRatio_counter_v, intime_hit_lengthRatio_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Hit Length Ratio", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, trk_len_shwr_len_ratio_counter_v, intime_trk_len_shwr_len_ratio_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "TrkLen/ShwrLen Rati", results_v);
	selection_functions::ExportEfficiencyPurity(total_mc_entries_inFV, track_containment_counter_v, intime_track_containment_counter_v->at(0),
	                                            intime_scale_factor, data_scale_factor, "Track Containment", results_v);

	selection_functions::PrintTopologyPurity(no_track, has_track, _1_shwr, _2_shwr, _3_shwr, _4_shwr);
	//*************************************************************************************************************************
	//*************************************************************************************************************************
	selection_functions::XSecWork(trk_len_shwr_len_ratio_counter_v->at(7), trk_len_shwr_len_ratio_counter_v->at(0), trk_len_shwr_len_ratio_counter_v->at(1),
	                              trk_len_shwr_len_ratio_counter_v->at(9), trk_len_shwr_len_ratio_counter_v->at(2),
	                              trk_len_shwr_len_ratio_counter_v->at(3), trk_len_shwr_len_ratio_counter_v->at(4),
	                              trk_len_shwr_len_ratio_counter_v->at(11), trk_len_shwr_len_ratio_counter_v->at(10), trk_len_shwr_len_ratio_counter_v->at(5),
	                              trk_len_shwr_len_ratio_counter_v->at(6), intime_trk_len_shwr_len_ratio_counter_v->at(0),
	                              intime_scale_factor, data_trk_len_shwr_len_ratio_counter_v->at(0), data_scale_factor,
	                              fv_boundary_v, flux, selected_energy_vector, genie_xsec, total_mc_entries_inFV);
	//*************************************************************************************************************************
	//*************************************************************************************************************************

	if(_post_cuts_verbose == true) {_functions_instance.selection_functions::PrintPostCutVector(post_cuts_v, _post_cuts_verbose); }

	gErrorIgnoreLevel = kWarning;

	std::cout << " --- End Cross Section Calculation --- " << std::endl;

	if(f->IsOpen()) {f->Close(); }
	//for some reason, the histogram is not being deleted upon function exit
	//trying to close root file instead

}//end selection
}//end namespace

#ifndef __ROOTCLING__

#endif
