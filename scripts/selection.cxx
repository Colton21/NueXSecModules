#include "selection.h"

namespace xsecSelection {
int selection( const char * _file1, const char * _file2, const char * _file3){

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


	//this is the old method - can be used as a cross check now
	mctruth_counter_tree->GetEntry(total_mc_entries-1);

	std::cout << "MC Nue CC Counter      : " << mc_nue_cc_counter << std::endl;
	std::cout << "MC Nue NC Counter      : " << mc_nue_nc_counter << std::endl;
	std::cout << "MC Numu CC Counter     : " << mc_numu_cc_counter << std::endl;
	std::cout << "MC Numu NC Counter     : " << mc_numu_nc_counter << std::endl;
	std::cout << "MC Nue CC Counter Bar  : " << mc_nue_cc_counter_bar << std::endl;
	std::cout << "MC Nue NC Counter Bar  : " << mc_nue_nc_counter_bar << std::endl;
	std::cout << "MC Numu CC Counter Bar : " << mc_numu_cc_counter_bar << std::endl;
	std::cout << "MC Numu NC Counter Bar : " << mc_numu_nc_counter_bar << std::endl;

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

	std::vector<TH1 * > * h_ele_pfp_xyz_data = new std::vector<TH1 * >;
	h_ele_pfp_xyz_data->push_back(h_ele_pfp_x_data);
	h_ele_pfp_xyz_data->push_back(h_ele_pfp_y_data);
	h_ele_pfp_xyz_data->push_back(h_ele_pfp_z_data);

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
		run_subrun_file.open("run_subrun_list.txt");
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

		_cuts_instance.selection_cuts::loop_flashes(data_f, data_optree, flash_pe_threshold, flash_time_start, flash_time_end, data_passed_runs);
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

			//** Testing leading shower length vs hits **//
			_data_functions_instance.selection_functions_data::ShowerLengthvsHitsData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_shwr_len_hits_data);

			//************************
			//******** in fv cut *****
			//************************
			_cuts_instance.selection_cuts::fiducial_volume_cut(data_tpc_object_container_v, fv_boundary_v, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_in_fv_counter_v);

			//*****************************
			//**** vertex to flash cut ****
			//*****************************
			_data_functions_instance.selection_functions_data::PostCutsVtxFlashData(largest_flash_v, data_tpc_object_container_v, passed_tpco_data,
			                                                                        h_vtx_flash_data);

			_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, data_tpc_object_container_v, tolerance, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_vtx_flash_counter_v);

			_data_functions_instance.selection_functions_data::PostCutsVtxFlashData(largest_flash_v, data_tpc_object_container_v, passed_tpco_data,
			                                                                        h_vtx_flash_data_after);
			//******************************************************
			//*** distance between pfp shower and nue object cut ***
			//******************************************************
			_data_functions_instance.selection_functions_data::PostCutsShwrVtxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_shwr_vtx_dist_data);

			_cuts_instance.selection_cuts::VtxNuDistance(data_tpc_object_container_v, shwr_nue_tolerance, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_shwr_tpco_counter_v);

			_data_functions_instance.selection_functions_data::PostCutsShwrVtxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_shwr_vtx_dist_data_after);

			//******************************************************
			// **** distance between pfp track and nue object cut **
			//******************************************************
			_data_functions_instance.selection_functions_data::PostCutTrkVtxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_trk_vtx_dist_data);
			_cuts_instance.selection_cuts::VtxTrackNuDistance(data_tpc_object_container_v, trk_nue_tolerance, passed_tpco_data, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_trk_tpco_counter_v);
			_data_functions_instance.selection_functions_data::PostCutTrkVtxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_trk_vtx_dist_data_after);

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
			//*****************************************************
			//*********** dEdx cut for the leading shower *********
			//******************************************************
			_data_functions_instance.selection_functions_data::PostCutsdEdxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_dedx_cuts_data);

			_cuts_instance.selection_cuts::dEdxCut(data_tpc_object_container_v, passed_tpco_data, tolerance_dedx_min, tolerance_dedx_max, _verbose);
			_data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco_data, tabulated_origins_data);
			_functions_instance.selection_functions::TotalOrigins(tabulated_origins_data, data_dedx_counter_v);

			_data_functions_instance.selection_functions_data::PostCutsdEdxData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_dedx_cuts_data_after);
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
			                                                                                 h_leading_shwr_trk_length_data);
			//*********** Data *********
			_functions_instance.selection_functions::FillPostCutVector(data_tpc_object_container_v, passed_tpco_data, post_cuts_v_data);
			//*************************************
			// ******** End Selection Cuts! *******
			//*************************************
			_data_functions_instance.selection_functions_data::LeadingMomentumData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_momentum_data);
			_data_functions_instance.selection_functions_data::LeadingPhiData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_phi_last_data);
			_data_functions_instance.selection_functions_data::LeadingThetaData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_theta_last_data);
			_data_functions_instance.selection_functions_data::LeadingCosThetaData(data_tpc_object_container_v, passed_tpco_data, theta_translation, phi_translation,
			                                                                       _verbose, h_ele_cos_theta_last_trans_data);
			_data_functions_instance.selection_functions_data::LeadingCosThetaData(data_tpc_object_container_v, passed_tpco_data, 0, 0,
			                                                                       _verbose, h_ele_cos_theta_last_data);
			_data_functions_instance.selection_functions_data::XYZPositionData(data_tpc_object_container_v, passed_tpco_data, _verbose, h_ele_pfp_xyz_data);
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

			//delete at the very end!
			delete passed_tpco_data;
		}
	}//end data running
	 //****************************
	 //*** END Data Calculation ***
	 //****************************

	//***********************************
	//*** In-time Cosmics Calculation ***
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

	std::vector<TH1 * > * h_ele_pfp_xyz_intime = new std::vector<TH1 * >;
	h_ele_pfp_xyz_intime->push_back(h_ele_pfp_x_intime);
	h_ele_pfp_xyz_intime->push_back(h_ele_pfp_y_intime);
	h_ele_pfp_xyz_intime->push_back(h_ele_pfp_z_intime);

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

		std::vector<xsecAna::TPCObjectContainer> * intime_tpc_object_container_v = nullptr;
		intime_tree->SetBranchAddress("TpcObjectContainerV", &intime_tpc_object_container_v);

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

		_cuts_instance.selection_cuts::loop_flashes(intime_f, intime_optree, flash_pe_threshold, flash_time_start, flash_time_end, intime_passed_runs);
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
				std::cout << "[EVENT NUMBER] \t " << event << std::endl;
				std::cout << "----------------------" << std::endl;
			}
			intime_tree->GetEntry(event);
			//***********************************************************
			//this is where the in-time optical cut actually takes effect
			//***********************************************************
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

			//** Testing leading shower length vs hits **//
			_functions_instance.selection_functions::ShowerLengthvsHitsInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_shwr_len_hits_intime);

			//************************
			//******** in fv cut *****
			//************************
			_cuts_instance.selection_cuts::fiducial_volume_cut(intime_tpc_object_container_v, fv_boundary_v, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_in_fv_counter_v);

			//*****************************
			//**** vertex to flash cut ****
			//*****************************
			_functions_instance.selection_functions::PostCutsVtxFlashInTime(largest_flash_v, intime_tpc_object_container_v, passed_tpco_intime,
			                                                                _verbose, h_vtx_flash_intime);

			_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, intime_tpc_object_container_v, tolerance, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_vtx_flash_counter_v);

			_functions_instance.selection_functions::PostCutsVtxFlashInTime(largest_flash_v, intime_tpc_object_container_v, passed_tpco_intime,
			                                                                _verbose, h_vtx_flash_intime_after);
			//******************************************************
			//*** distance between pfp shower and nue object cut ***
			//******************************************************
			_functions_instance.selection_functions::PostCutsShwrVtxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_shwr_vtx_dist_intime);

			_cuts_instance.selection_cuts::VtxNuDistance(intime_tpc_object_container_v, shwr_nue_tolerance, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_shwr_tpco_counter_v);

			_functions_instance.selection_functions::PostCutsShwrVtxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_shwr_vtx_dist_intime_after);
			//******************************************************
			// **** distance between pfp track and nue object cut **
			//******************************************************
			_functions_instance.selection_functions::PostCutTrkVtxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_trk_vtx_dist_intime);

			_cuts_instance.selection_cuts::VtxTrackNuDistance(intime_tpc_object_container_v, trk_nue_tolerance, passed_tpco_intime, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_trk_tpco_counter_v);

			_functions_instance.selection_functions::PostCutTrkVtxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_trk_vtx_dist_intime_after);
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
			//*****************************************************
			//****** open angle cut for the leading shower ********
			//******************************************************
			_functions_instance.selection_functions::NumShowersOpenAngleInTime(intime_tpc_object_container_v, passed_tpco_intime, h_pfp_shower_open_angle_intime);
			_functions_instance.selection_functions::PostCutOpenAngleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                _verbose, h_leading_shower_open_angle_intime);
			_functions_instance.selection_functions::PostCutOpenAngle1ShowerInTime(tpc_object_container_v, passed_tpco_intime,
			                                                                       _verbose, h_leading_shower_open_angle_1_intime);
			_functions_instance.selection_functions::PostCutOpenAngle2PlusShowerInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                           _verbose, h_leading_shower_open_angle_2plus_intime);
			_cuts_instance.selection_cuts::OpenAngleCut(intime_tpc_object_container_v, passed_tpco_intime, tolerance_open_angle, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_open_angle_counter_v);

			_functions_instance.selection_functions::NumShowersOpenAngleInTime(intime_tpc_object_container_v, passed_tpco_intime, h_pfp_shower_dedx_intime);

			_functions_instance.selection_functions::PostCutOpenAngleInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                                _verbose, h_leading_shower_open_angle_intime_after);
			//*****************************************************
			//*********** dEdx cut for the leading shower *********
			//******************************************************
			_functions_instance.selection_functions::PostCutsdEdxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_dedx_cuts_intime);

			_cuts_instance.selection_cuts::dEdxCut(intime_tpc_object_container_v, passed_tpco_intime, tolerance_dedx_min, tolerance_dedx_max, _verbose);
			_functions_instance.selection_functions::TabulateOriginsInTime(intime_tpc_object_container_v, passed_tpco_intime, tabulated_origins_intime);
			_functions_instance.selection_functions::TotalOriginsInTime(tabulated_origins_intime, intime_dedx_counter_v);

			_functions_instance.selection_functions::PostCutsdEdxInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_dedx_cuts_intime_after);

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
			//*********** In-time Cosmics *********
			_functions_instance.selection_functions::FillPostCutVector(intime_tpc_object_container_v, passed_tpco_intime, post_cuts_v);
			//*************************************
			// ******** End Selection Cuts! *******
			//*************************************
			_functions_instance.selection_functions::LeadingMomentumInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_ele_pfp_momentum_intime);
			_functions_instance.selection_functions::LeadingPhiInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_ele_pfp_phi_last_intime);
			_functions_instance.selection_functions::LeadingThetaInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                            theta_translation, phi_translation, _verbose, h_ele_pfp_theta_last_intime);
			_functions_instance.selection_functions::LeadingCosThetaInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                               0, 0, _verbose, h_ele_cos_theta_last_intime);
			_functions_instance.selection_functions::LeadingCosThetaInTime(intime_tpc_object_container_v, passed_tpco_intime,
			                                                               theta_translation, phi_translation, _verbose, h_ele_cos_theta_last_trans_intime);
			_functions_instance.selection_functions::XYZPositionInTime(intime_tpc_object_container_v, passed_tpco_intime, _verbose, h_ele_pfp_xyz_intime);
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

	_cuts_instance.selection_cuts::loop_flashes(f, optree, flash_pe_threshold, flash_time_start, flash_time_end, passed_runs);
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
		//const bool true_in_tpc = _cuts_instance.selection_cuts::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v);
		const bool true_in_tpc = true_in_tpc_v.at(event);

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
		const double mc_ele_phi       = atan2(mc_ele_dir_y, mc_ele_dir_x);
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
		//***********************************************************
		//this is where the in-time optical cut again takes effect
		//***********************************************************
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, in_time_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_no_cut, h_selected_ele_energy_no_cut);
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
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_reco_nue, h_selected_ele_energy_reco_nue);

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
			                                                             mc_nu_energy, mc_ele_energy, h_shwr_hits_nu_eng, h_shwr_hits_ele_eng);
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
			                                                             mc_nu_energy, mc_ele_energy, h_shwr_hits_nu_eng_zoom, h_shwr_hits_ele_eng_zoom);

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

		//we also want to look at the cos(theta) and energy efficiency before we make selection cuts
		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins->at(0) >= 1 && true_in_tpc == true)
		{
			h_ele_eng_eff_num_pre_cuts->Fill(mc_ele_energy);
			h_ele_cos_theta_eff_num_pre_cuts->Fill(mc_ele_cos_theta);
		}

		//*****************************
		//**** vertex to flash cut ****
		//*****************************
		_functions_instance.selection_functions::PostCutsVtxFlash(largest_flash_v, tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                          h_vtx_flash_nue_cc, h_vtx_flash_nue_cc_out_fv, h_vtx_flash_nue_cc_mixed,
		                                                          h_vtx_flash_numu_cc, h_vtx_flash_nc,
		                                                          h_vtx_flash_cosmic, h_vtx_flash_nc_pi0,
		                                                          h_vtx_flash_numu_cc_mixed, h_vtx_flash_other_mixed,
		                                                          h_vtx_flash_unmatched);

		_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, tpc_object_container_v, tolerance, passed_tpco, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, vtx_flash_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_vtx_flash, h_selected_ele_energy_vtx_flash);

		_functions_instance.selection_functions::PostCutsVtxFlash(largest_flash_v, tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                          h_vtx_flash_nue_cc_after, h_vtx_flash_nue_cc_out_fv_after, h_vtx_flash_nue_cc_mixed_after,
		                                                          h_vtx_flash_numu_cc_after, h_vtx_flash_nc_after,
		                                                          h_vtx_flash_cosmic_after, h_vtx_flash_nc_pi0_after,
		                                                          h_vtx_flash_numu_cc_mixed_after, h_vtx_flash_other_mixed_after,
		                                                          h_vtx_flash_unmatched_after);
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

		_functions_instance.selection_functions::PostCutTrkVtx(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                       h_trk_vtx_dist_nue_cc_after, h_trk_vtx_dist_nue_cc_out_fv_after,
		                                                       h_trk_vtx_dist_nue_cc_mixed_after,
		                                                       h_trk_vtx_dist_numu_cc_after, h_trk_vtx_dist_nc_after,
		                                                       h_trk_vtx_dist_cosmic_after, h_trk_vtx_dist_nc_pi0_after,
		                                                       h_trk_vtx_dist_numu_cc_mixed_after, h_trk_vtx_dist_other_mixed_after,
		                                                       h_trk_vtx_dist_unmatched_after);
		//****************************************************
		// ******** hit threshold for showers cut *************
		//******************************************************
		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins->at(0) >= 1 && true_in_tpc == true)
		{
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose,
			                                                             tpco_classifier_v, mc_nu_energy, mc_ele_energy,
			                                                             h_shwr_hits_nu_eng_last, h_shwr_hits_ele_eng_last);
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose,
			                                                             tpco_classifier_v, mc_nu_energy, mc_ele_energy,
			                                                             h_shwr_hits_nu_eng_zoom_last, h_shwr_hits_ele_eng_zoom_last);
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

		//*****************************************************
		//*********** dEdx cut for the leading shower *********
		//******************************************************
		_functions_instance.selection_functions::PostCutsdEdx(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                      h_dedx_cuts_nue_cc, h_dedx_cuts_nue_cc_mixed,
		                                                      h_dedx_cuts_nue_cc_out_fv,
		                                                      h_dedx_cuts_numu_cc, h_dedx_cuts_nc,
		                                                      h_dedx_cuts_cosmic, h_dedx_cuts_nc_pi0,
		                                                      h_dedx_cuts_numu_cc_mixed, h_dedx_cuts_other_mixed,
		                                                      h_dedx_cuts_unmatched);

		_cuts_instance.selection_cuts::dEdxCut(tpc_object_container_v, passed_tpco, tolerance_dedx_min, tolerance_dedx_max, _verbose);
		_functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, tabulated_origins, tpco_classifier_v);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, dedx_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   fv_boundary_v,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_dedx, h_selected_ele_energy_dedx);

		_functions_instance.selection_functions::PostCutsdEdx(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                      h_dedx_cuts_nue_cc_after, h_dedx_cuts_nue_cc_mixed_after,
		                                                      h_dedx_cuts_nue_cc_out_fv_after,
		                                                      h_dedx_cuts_numu_cc_after, h_dedx_cuts_nc_after,
		                                                      h_dedx_cuts_cosmic_after, h_dedx_cuts_nc_pi0_after,
		                                                      h_dedx_cuts_numu_cc_mixed_after, h_dedx_cuts_other_mixed_after,
		                                                      h_dedx_cuts_unmatched_after);

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

		_functions_instance.selection_functions::SecondaryShowersDist(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                              h_second_shwr_dist_nue_cc_after, h_second_shwr_dist_nue_cc_out_fv_after,
		                                                              h_second_shwr_dist_nue_cc_mixed_after, h_second_shwr_dist_numu_cc_after,
		                                                              h_second_shwr_dist_numu_cc_mixed_after, h_second_shwr_dist_nc_after,
		                                                              h_second_shwr_dist_nc_pi0_after, h_second_shwr_dist_cosmic_after,
		                                                              h_second_shwr_dist_other_mixed_after, h_second_shwr_dist_unmatched_after);

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

		_functions_instance.selection_functions::FailureReason(tpc_object_container_v, passed_tpco, _verbose, tpco_classifier_v,
		                                                       h_failure_reason_nue_cc, h_failure_reason_nue_cc_out_fv,
		                                                       h_failure_reason_nue_cc_mixed, h_failure_reason_numu_cc,
		                                                       h_failure_reason_numu_cc_mixed, h_failure_reason_nc,
		                                                       h_failure_reason_nc_pi0, h_failure_reason_cosmic,
		                                                       h_failure_reason_other_mixed, h_failure_reason_unmatched);

		//these are for the tefficiency plots, post all cuts
		if((mc_nu_id == 1 || mc_nu_id == 5) && true_in_tpc == true)
		{
			int pass = 0;
			for(auto const passed_tpco_pair : * passed_tpco) {pass += passed_tpco_pair.first; }
			if(pass > 0)
			{
				selected_energy_vector.push_back(mc_nu_energy);
				h_nue_eng_eff_num->Fill(mc_nu_energy);
				h_ele_eng_eff_num->Fill(mc_ele_energy);
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
		_functions_instance.selection_functions::LeadingCosTheta(tpc_object_container_v, passed_tpco, 0, 0, _verbose, tpco_classifier_v,
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
	selection_functions::PrintInfo( total_mc_entries_inFV, in_time_counter_v,
	                                intime_scale_factor * intime_in_time_counter_v->at(0),                      "In Time");
	selection_functions_data::PrintInfoData(data_scale_factor * data_in_time_counter_v->at(0),                  "In Time");
	selection_functions::PrintInfo( total_mc_entries_inFV, pe_counter_v,
	                                intime_scale_factor * intime_pe_counter_v->at(0),                           "PE Threshold");
	selection_functions_data::PrintInfoData(data_scale_factor * data_pe_counter_v->at(0),                       "PE Threshold");
	selection_functions::PrintInfo( total_mc_entries_inFV, reco_nue_counter_v,
	                                intime_scale_factor * intime_reco_nue_counter_v->at(0),                     "Reco Nue");
	selection_functions_data::PrintInfoData(data_scale_factor * data_reco_nue_counter_v->at(0),                 "Reco Nue");
	selection_functions::PrintInfo( total_mc_entries_inFV, in_fv_counter_v,
	                                intime_scale_factor * intime_in_fv_counter_v->at(0),                        "In FV");
	selection_functions_data::PrintInfoData(data_scale_factor * data_in_fv_counter_v->at(0),                    "In FV");
	selection_functions::PrintInfo( total_mc_entries_inFV, vtx_flash_counter_v,
	                                intime_scale_factor * intime_vtx_flash_counter_v->at(0),                    "Vtx-to-Flash");
	selection_functions_data::PrintInfoData(data_scale_factor * data_vtx_flash_counter_v->at(0),                "Vtx-to-Flash");
	selection_functions::PrintInfo( total_mc_entries_inFV, shwr_tpco_counter_v,
	                                intime_scale_factor * intime_shwr_tpco_counter_v->at(0),                    "Shower-to-TPCO");
	selection_functions_data::PrintInfoData(data_scale_factor * data_shwr_tpco_counter_v->at(0),                "Shower-to-TPCO");
	selection_functions::PrintInfo( total_mc_entries_inFV, trk_tpco_counter_v,
	                                intime_scale_factor * intime_trk_tpco_counter_v->at(0),                     "Track-to-TPCO");
	selection_functions_data::PrintInfoData(data_scale_factor * data_trk_tpco_counter_v->at(0),                 "Track-to-TPCO");
	selection_functions::PrintInfo( total_mc_entries_inFV, hit_threshold_counter_v,
	                                intime_scale_factor * intime_hit_threshold_counter_v->at(0),                "Hit Threshold");
	selection_functions_data::PrintInfoData(data_scale_factor * data_hit_threshold_counter_v->at(0),            "Hit Threshold");
	selection_functions::PrintInfo( total_mc_entries_inFV, hit_threshold_collection_counter_v,
	                                intime_scale_factor * intime_hit_threshold_collection_counter_v->at(0),     "WPlane Hit Threshold");
	selection_functions_data::PrintInfoData(data_scale_factor * data_hit_threshold_collection_counter_v->at(0), "WPlane hit Threshold");
	selection_functions::PrintInfo( total_mc_entries_inFV, open_angle_counter_v,
	                                intime_scale_factor * intime_open_angle_counter_v->at(0),                   "Open Angle");
	selection_functions_data::PrintInfoData(data_scale_factor * data_open_angle_counter_v->at(0),               "Open Angle");
	selection_functions::PrintInfo( total_mc_entries_inFV, dedx_counter_v,
	                                intime_scale_factor * intime_dedx_counter_v->at(0),                         " dE / dx ");
	selection_functions_data::PrintInfoData(data_scale_factor * data_dedx_counter_v->at(0),                     " dE / dx ");
	selection_functions::PrintInfo( total_mc_entries_inFV, secondary_shower_counter_v,
	                                intime_scale_factor * intime_secondary_shower_counter_v->at(0),             ">1 Shower TPCO Dist");
	selection_functions_data::PrintInfoData(data_scale_factor * data_secondary_shower_counter_v->at(0),         ">1 Shower TPCO Dist");
	selection_functions::PrintInfo( total_mc_entries_inFV, hit_lengthRatio_counter_v,
	                                intime_scale_factor * intime_hit_lengthRatio_counter_v->at(0),              "Hit Length Ratio");
	selection_functions_data::PrintInfoData(data_scale_factor * data_hit_lengthRatio_counter_v->at(0),          "Hit Length Ratio");
	selection_functions::PrintInfo( total_mc_entries_inFV, trk_len_shwr_len_ratio_counter_v,
	                                intime_scale_factor * intime_trk_len_shwr_len_ratio_counter_v->at(0),       "TrkLen/ShwrLen Ratio");
	selection_functions_data::PrintInfoData(data_scale_factor * data_trk_len_shwr_len_ratio_counter_v->at(0),   "TrkLen/ShwrLen Ratio");

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

	_functions_instance.selection_functions::PostCutVector2DPlots(post_cuts_v, _post_cuts_verbose, intime_scale_factor,
	                                                              h_post_cuts_num_tracks_showers_purity_qe,
	                                                              h_post_cuts_num_tracks_showers_purity_res,
	                                                              h_post_cuts_num_tracks_showers_purity_dis,
	                                                              h_post_cuts_num_tracks_showers_purity_coh,
	                                                              h_post_cuts_num_tracks_showers_purity_mec,
	                                                              h_post_cuts_num_tracks_showers_purity_total,
	                                                              h_post_cuts_num_tracks_showers_signal_total,
	                                                              h_post_cuts_num_tracks_showers_bkg_total,
	                                                              h_post_cuts_num_tracks_showers_total_total);

//********************//
//**** Histograms ****//
//*******************//

	//histogram_functions::1DHistogram (TH1 histogram, const std::string x_axis_name, const std::string * print_name)

	histogram_functions::Plot1DHistogram (h_nue_eng_eff_num, "True Neutrino Energy [GeV]", "selected_true_neutrino_energy.pdf");
	histogram_functions::Plot1DHistogram (h_nue_eng_eff_den, "True Neutrino Energy [GeV]", "all_true_neutrino_energy.pdf");
	histogram_functions::Plot1DHistogram (h_nue_num_part_eff_den, "True Particle Multiplicity", "all_true_neutrino_num_particles.pdf");
	histogram_functions::Plot1DHistogram (h_nue_num_part_eff_num, "Selected True Particle Multiplicity", "selected_true_neutrino_num_particles.pdf");
	histogram_functions::Plot1DHistogram (h_nue_num_chrg_part_eff_den, "True Charged Particle Multiplicity", "all_true_neutrino_num_charged_particles.pdf");
	histogram_functions::Plot1DHistogram (h_nue_num_chrg_part_eff_num, "Selected Charged True Particle Multiplicity",
	                                      "selected_true_neutrino_num_charged_particles.pdf");
	histogram_functions::Plot1DHistogram(h_nue_true_theta, "True Neutrino Theta [Degrees]", "true_nue_theta.pdf");
	histogram_functions::Plot1DHistogram(h_nue_true_phi,   "True Neutrino Phi [Degrees]",   "true_nue_phi.pdf");
	histogram_functions::Plot2DHistogram(h_nue_true_theta_phi, "True Nue CC", "Phi [Degrees]", "Theta [Degrees]", "true_nue_theta_phi.pdf");

	histogram_functions::Plot1DHistogram(h_ele_phi_eff_num, "Selected True Electron Phi [Degrees]", "selected_true_electron_phi.pdf");
	histogram_functions::Plot1DHistogram(h_ele_phi_eff_den, "True Electron Phi [Degrees]", "true_electron_phi.pdf");
	histogram_functions::Plot1DHistogram(h_ele_cos_theta_eff_num, "Selected True Electron Cos(#theta)", "selected_true_electron_cos_theta.pdf");
	histogram_functions::Plot1DHistogram(h_ele_cos_theta_eff_den, "True Electron Cos(#theta)", "true_electron_cos_theta.pdf");
	histogram_functions::Plot1DHistogram(h_ele_theta_eff_num, "Selected True Electron Theta [Degrees]", "selected_true_electron_theta.pdf");
	histogram_functions::Plot1DHistogram(h_ele_theta_eff_den, "True Electron Theta [Degrees]", "true_electron_theta.pdf");

	histogram_functions::PlotTEfficiency (h_nue_eng_eff_num, h_nue_eng_eff_den,
	                                      ";True Neutrino Energy [GeV];Efficiency", "signal_selection_nu_energy_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_vtx_x_eff_num, h_nue_vtx_x_eff_den,
	                                      ";True Neutrino Vtx X [cm];Efficiency", "signal_selection_nu_vtx_x_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_vtx_y_eff_num, h_nue_vtx_y_eff_den,
	                                      ";True Neutrino Vtx Y [cm];Efficiency", "signal_selection_nu_vtx_y_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_vtx_z_eff_num, h_nue_vtx_z_eff_den,
	                                      ";True Neutrino Vtx Z [cm];Efficiency", "signal_selection_nu_vtx_z_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_dir_x_eff_num, h_nue_dir_x_eff_den,
	                                      ";True Neutrino Dir X;Efficiency", "signal_selection_nu_dir_x_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_dir_y_eff_num, h_nue_dir_y_eff_den,
	                                      ";True Neutrino Dir Y;Efficiency", "signal_selection_nu_dir_y_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_dir_z_eff_num, h_nue_dir_z_eff_den,
	                                      ";True Neutrino Dir Z;Efficiency", "signal_selection_nu_dir_z_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_num_part_eff_num, h_nue_num_part_eff_den,
	                                      ";True Particle Multiplicity;Efficiency", "signal_selection_nu_num_particles_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_num_chrg_part_eff_num, h_nue_num_chrg_part_eff_den,
	                                      ";True Charged Particle Multiplicity;Efficiency", "signal_selection_nu_num_charged_particles_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_cos_theta_eff_num, h_nue_cos_theta_eff_den,
	                                      ";True Neutrino Cos(#theta);Efficiency", "signal_selection_nu_cos_theta_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_nue_phi_eff_num, h_nue_phi_eff_den, ";True Neutrino Phi;Efficiency", "signal_selection_nu_phi_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_ele_cos_theta_eff_num, h_ele_cos_theta_eff_den,
	                                      ";True Nue Electron Cos(#theta);Efficiency", "signal_selection_ele_cos_theta_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_ele_phi_eff_num, h_ele_phi_eff_den, ";True Nue Electron Phi;Efficiency", "signal_selection_ele_phi_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_ele_dir_x_eff_num, h_ele_dir_x_eff_den,
	                                      ";True Nue Electron Dir X;Efficiency", "signal_selection_ele_dir_x_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_ele_dir_y_eff_num, h_ele_dir_y_eff_den,
	                                      ";True Nue Electron Dir Y;Efficiency", "signal_selection_ele_dir_y_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_ele_dir_z_eff_num, h_ele_dir_z_eff_den,
	                                      ";True Nue Electron Dir Z;Efficiency", "signal_selection_ele_dir_z_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_ele_eng_eff_num, h_ele_eng_eff_den, ";True Nue Electron Energy;Efficiency",
	                                      "signal_selection_ele_energy_efficiency.pdf");
	histogram_functions::PlotTEfficiency (h_post_cuts_num_tracks_showers_signal_total, h_post_cuts_num_tracks_showers_bkg_total,
	                                      ";Num Reco Showers; Num Reco Tracks", "post_cuts_num_tracks_showers_eff_purity_total.pdf");
	histogram_functions::PlotTEfficiency(h_ele_cos_theta_eff_num_pre_cuts, h_ele_cos_theta_eff_den,
	                                     ";True Electron Cos(#theta);Efficiency", "signal_selection_ele_cos_theta_pre_cuts_efficiency.pdf");
	histogram_functions::PlotTEfficiency(h_ele_eng_eff_num_pre_cuts, h_ele_eng_eff_den,
	                                     ";True Electron Energy;Efficiency", "signal_selection_ele_eng_pre_cuts_efficiency.pdf");

	histogram_functions::Plot2DHistogram (h_tracks_showers, "Post Cuts - Showers/Tracks per Candidate Nue TPC Object",
	                                      "Reco Tracks", "Reco Showers", "post_cuts_showers_tracks.pdf");
	histogram_functions::Plot2DHistogram (h_tracks_showers, "Post Cuts - Showers/Tracks per Candidate Nue TPC Object",
	                                      "Reco Tracks", "Reco Showers", "post_cuts_showers_tracks_cosmic.pdf");
	histogram_functions::Plot2DHistogram (h_tracks_showers, "Post Cuts - Showers/Tracks per Candidate Nue TPC Object",
	                                      "Reco Tracks", "Reco Showers", "post_cuts_showers_tracks_numu.pdf");

	histogram_functions::PlotSimpleStack (h_leading_shower_open_angle_nue_cc,  h_leading_shower_open_angle_nue_cc_mixed,
	                                      h_leading_shower_open_angle_nue_cc_out_fv,
	                                      h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_numu_cc_mixed,
	                                      h_leading_shower_open_angle_cosmic,  h_leading_shower_open_angle_nc,
	                                      h_leading_shower_open_angle_nc_pi0,  h_leading_shower_open_angle_other_mixed,
	                                      h_leading_shower_open_angle_unmatched, "",
	                                      "Shower Opening Angle [Degrees]", "", "post_cuts_leading_shower_open_angle.pdf");
	histogram_functions::PlotSimpleStackInTime (h_leading_shower_open_angle_nue_cc,  h_leading_shower_open_angle_nue_cc_mixed,
	                                            h_leading_shower_open_angle_nue_cc_out_fv,
	                                            h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_numu_cc_mixed,
	                                            h_leading_shower_open_angle_cosmic,  h_leading_shower_open_angle_nc,
	                                            h_leading_shower_open_angle_nc_pi0,  h_leading_shower_open_angle_other_mixed,
	                                            h_leading_shower_open_angle_unmatched, h_leading_shower_open_angle_intime, intime_scale_factor,
	                                            "", "Shower Opening Angle [Degrees]", "", "post_cuts_leading_shower_open_angle_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_leading_shower_open_angle_nue_cc,  h_leading_shower_open_angle_nue_cc_mixed,
	                                          h_leading_shower_open_angle_nue_cc_out_fv,
	                                          h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_numu_cc_mixed,
	                                          h_leading_shower_open_angle_cosmic,  h_leading_shower_open_angle_nc,
	                                          h_leading_shower_open_angle_nc_pi0,  h_leading_shower_open_angle_other_mixed,
	                                          h_leading_shower_open_angle_unmatched, h_leading_shower_open_angle_intime, intime_scale_factor,
	                                          h_leading_shower_open_angle_data, data_scale_factor,
	                                          "", "Shower Opening Angle [Degrees]", "", "post_cuts_leading_shower_open_angle_data.pdf");
	histogram_functions::PlotSimpleStackData (h_leading_shower_open_angle_nue_cc_after,  h_leading_shower_open_angle_nue_cc_mixed_after,
	                                          h_leading_shower_open_angle_nue_cc_out_fv_after,
	                                          h_leading_shower_open_angle_numu_cc_after, h_leading_shower_open_angle_numu_cc_mixed_after,
	                                          h_leading_shower_open_angle_cosmic_after,  h_leading_shower_open_angle_nc_after,
	                                          h_leading_shower_open_angle_nc_pi0_after,  h_leading_shower_open_angle_other_mixed_after,
	                                          h_leading_shower_open_angle_unmatched_after, h_leading_shower_open_angle_intime_after, intime_scale_factor,
	                                          h_leading_shower_open_angle_data_after, data_scale_factor,
	                                          "", "Shower Opening Angle [Degrees]", "", "post_cuts_leading_shower_open_angle_data_after.pdf");

	histogram_functions::PlotSimpleStack (h_leading_shower_open_angle_1_nue_cc,  h_leading_shower_open_angle_1_nue_cc_mixed,
	                                      h_leading_shower_open_angle_1_nue_cc_out_fv,
	                                      h_leading_shower_open_angle_1_numu_cc, h_leading_shower_open_angle_1_numu_cc_mixed,
	                                      h_leading_shower_open_angle_1_cosmic,  h_leading_shower_open_angle_1_nc,
	                                      h_leading_shower_open_angle_1_nc_pi0,  h_leading_shower_open_angle_1_other_mixed,
	                                      h_leading_shower_open_angle_1_unmatched, "",
	                                      "Shower Opening Angle [Degrees]", "", "post_cuts_leading_shower_open_angle_1_shower.pdf");

	histogram_functions::PlotSimpleStackInTime (h_leading_shower_open_angle_1_nue_cc,  h_leading_shower_open_angle_1_nue_cc_mixed,
	                                            h_leading_shower_open_angle_1_nue_cc_out_fv,
	                                            h_leading_shower_open_angle_1_numu_cc, h_leading_shower_open_angle_1_numu_cc_mixed,
	                                            h_leading_shower_open_angle_1_cosmic,  h_leading_shower_open_angle_1_nc,
	                                            h_leading_shower_open_angle_1_nc_pi0,  h_leading_shower_open_angle_1_other_mixed,
	                                            h_leading_shower_open_angle_1_unmatched, h_leading_shower_open_angle_1_intime, intime_scale_factor,
	                                            "", "Shower Opening Angle [Degrees]", "", "post_cuts_leading_shower_open_angle_1_shower_intime.pdf");

	histogram_functions::PlotSimpleStack (h_leading_shower_open_angle_2plus_nue_cc,  h_leading_shower_open_angle_2plus_nue_cc_mixed,
	                                      h_leading_shower_open_angle_2plus_nue_cc_out_fv,
	                                      h_leading_shower_open_angle_2plus_numu_cc, h_leading_shower_open_angle_2plus_numu_cc_mixed,
	                                      h_leading_shower_open_angle_2plus_cosmic,  h_leading_shower_open_angle_2plus_nc,
	                                      h_leading_shower_open_angle_2plus_nc_pi0,  h_leading_shower_open_angle_2plus_other_mixed,
	                                      h_leading_shower_open_angle_2plus_unmatched, "",
	                                      "Shower Opening Angle [Degrees]", "", "post_cuts_leading_shower_open_angle_2plus_showers.pdf");

	histogram_functions::PlotSimpleStackInTime (h_leading_shower_open_angle_2plus_nue_cc,  h_leading_shower_open_angle_2plus_nue_cc_mixed,
	                                            h_leading_shower_open_angle_2plus_nue_cc_out_fv,
	                                            h_leading_shower_open_angle_2plus_numu_cc, h_leading_shower_open_angle_2plus_numu_cc_mixed,
	                                            h_leading_shower_open_angle_2plus_cosmic,  h_leading_shower_open_angle_2plus_nc,
	                                            h_leading_shower_open_angle_2plus_nc_pi0,  h_leading_shower_open_angle_2plus_other_mixed,
	                                            h_leading_shower_open_angle_2plus_unmatched, h_leading_shower_open_angle_2plus_intime,
	                                            intime_scale_factor, "", "Shower Opening Angle [Degrees]", "",
	                                            "post_cuts_leading_shower_open_angle_2plus_showers_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_leading_shower_open_angle_2plus_nue_cc,  h_leading_shower_open_angle_2plus_nue_cc_mixed,
	                                          h_leading_shower_open_angle_2plus_nue_cc_out_fv,
	                                          h_leading_shower_open_angle_2plus_numu_cc, h_leading_shower_open_angle_2plus_numu_cc_mixed,
	                                          h_leading_shower_open_angle_2plus_cosmic,  h_leading_shower_open_angle_2plus_nc,
	                                          h_leading_shower_open_angle_2plus_nc_pi0,  h_leading_shower_open_angle_2plus_other_mixed,
	                                          h_leading_shower_open_angle_2plus_unmatched, h_leading_shower_open_angle_2plus_intime,
	                                          intime_scale_factor, h_leading_shower_open_angle_2plus_data, data_scale_factor,
	                                          "", "Shower Opening Angle [Degrees]", "",
	                                          "post_cuts_leading_shower_open_angle_2plus_showers_data.pdf");

	histogram_functions::PlotSimpleStack (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                      h_dedx_cuts_nue_cc_out_fv,
	                                      h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                      h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                      h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                      h_dedx_cuts_unmatched, "",
	                                      "Collection Plane dE/dx [MeV/cm]", "", "post_cuts_dedx_cuts.pdf");

	histogram_functions::PlotSimpleStackInTime (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                            h_dedx_cuts_nue_cc_out_fv,
	                                            h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                            h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                            h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                            h_dedx_cuts_unmatched, h_dedx_cuts_intime, intime_scale_factor, "",
	                                            "Collection Plane dE/dx [MeV/cm]", "", "post_cuts_dedx_cuts_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_dedx_cuts_nue_cc,  h_dedx_cuts_nue_cc_mixed,
	                                          h_dedx_cuts_nue_cc_out_fv,
	                                          h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                          h_dedx_cuts_cosmic,  h_dedx_cuts_nc,
	                                          h_dedx_cuts_nc_pi0,  h_dedx_cuts_other_mixed,
	                                          h_dedx_cuts_unmatched, h_dedx_cuts_intime, intime_scale_factor,
	                                          h_dedx_cuts_data, data_scale_factor, "",
	                                          "Collection Plane dE/dx [MeV/cm]", "", "post_cuts_dedx_cuts_data.pdf");
	histogram_functions::PlotSimpleStackData (h_dedx_cuts_nue_cc_after,  h_dedx_cuts_nue_cc_mixed_after,
	                                          h_dedx_cuts_nue_cc_out_fv_after,
	                                          h_dedx_cuts_numu_cc_after, h_dedx_cuts_numu_cc_mixed_after,
	                                          h_dedx_cuts_cosmic_after,  h_dedx_cuts_nc_after,
	                                          h_dedx_cuts_nc_pi0_after,  h_dedx_cuts_other_mixed_after,
	                                          h_dedx_cuts_unmatched_after, h_dedx_cuts_intime_after, intime_scale_factor,
	                                          h_dedx_cuts_data_after, data_scale_factor, "",
	                                          "Collection Plane dE/dx [MeV/cm]", "", "post_cuts_dedx_cuts_data_after.pdf");

	histogram_functions::PlotSimpleStack (h_vtx_flash_nue_cc,  h_vtx_flash_nue_cc_mixed,
	                                      h_vtx_flash_nue_cc_out_fv,
	                                      h_vtx_flash_numu_cc, h_vtx_flash_numu_cc_mixed,
	                                      h_vtx_flash_cosmic,  h_vtx_flash_nc,
	                                      h_vtx_flash_nc_pi0,  h_vtx_flash_other_mixed,
	                                      h_vtx_flash_unmatched, "",
	                                      "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "", "post_cuts_vtx_to_flash_distance.pdf");

	histogram_functions::PlotSimpleStackInTime (h_vtx_flash_nue_cc,  h_vtx_flash_nue_cc_mixed,
	                                            h_vtx_flash_nue_cc_out_fv,
	                                            h_vtx_flash_numu_cc, h_vtx_flash_numu_cc_mixed,
	                                            h_vtx_flash_cosmic,  h_vtx_flash_nc,
	                                            h_vtx_flash_nc_pi0,  h_vtx_flash_other_mixed,
	                                            h_vtx_flash_unmatched, h_vtx_flash_intime, intime_scale_factor, "",
	                                            "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "", "post_cuts_vtx_to_flash_distance_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_vtx_flash_nue_cc,  h_vtx_flash_nue_cc_mixed,
	                                          h_vtx_flash_nue_cc_out_fv,
	                                          h_vtx_flash_numu_cc, h_vtx_flash_numu_cc_mixed,
	                                          h_vtx_flash_cosmic,  h_vtx_flash_nc,
	                                          h_vtx_flash_nc_pi0,  h_vtx_flash_other_mixed,
	                                          h_vtx_flash_unmatched, h_vtx_flash_intime, intime_scale_factor,
	                                          h_vtx_flash_data, data_scale_factor, "",
	                                          "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "", "post_cuts_vtx_to_flash_distance_data.pdf");

	histogram_functions::PlotSimpleStackData (h_vtx_flash_nue_cc_after,  h_vtx_flash_nue_cc_mixed_after,
	                                          h_vtx_flash_nue_cc_out_fv_after,
	                                          h_vtx_flash_numu_cc_after, h_vtx_flash_numu_cc_mixed_after,
	                                          h_vtx_flash_cosmic_after,  h_vtx_flash_nc_after,
	                                          h_vtx_flash_nc_pi0_after,  h_vtx_flash_other_mixed_after,
	                                          h_vtx_flash_unmatched_after, h_vtx_flash_intime_after, intime_scale_factor,
	                                          h_vtx_flash_data_after, data_scale_factor, "",
	                                          "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "", "post_cuts_vtx_to_flash_distance_data_after.pdf");

	histogram_functions::PlotSimpleStack (h_trk_vtx_dist_nue_cc,  h_trk_vtx_dist_nue_cc_mixed,
	                                      h_trk_vtx_dist_nue_cc_out_fv,
	                                      h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed,
	                                      h_trk_vtx_dist_cosmic,  h_trk_vtx_dist_nc,
	                                      h_trk_vtx_dist_nc_pi0,  h_trk_vtx_dist_other_mixed,
	                                      h_trk_vtx_dist_unmatched, "",
	                                      "Track to Nue Candidate Vertex Distance [cm]", "", "post_cuts_track_to_vtx.pdf");

	histogram_functions::PlotSimpleStackInTime (h_trk_vtx_dist_nue_cc,  h_trk_vtx_dist_nue_cc_mixed,
	                                            h_trk_vtx_dist_nue_cc_out_fv,
	                                            h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed,
	                                            h_trk_vtx_dist_cosmic,  h_trk_vtx_dist_nc,
	                                            h_trk_vtx_dist_nc_pi0,  h_trk_vtx_dist_other_mixed,
	                                            h_trk_vtx_dist_unmatched, h_trk_vtx_dist_intime, intime_scale_factor, "",
	                                            "Track to Nue Candidate Vertex Distance [cm]", "", "post_cuts_track_to_vtx_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_trk_vtx_dist_nue_cc,  h_trk_vtx_dist_nue_cc_mixed,
	                                          h_trk_vtx_dist_nue_cc_out_fv,
	                                          h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed,
	                                          h_trk_vtx_dist_cosmic,  h_trk_vtx_dist_nc,
	                                          h_trk_vtx_dist_nc_pi0,  h_trk_vtx_dist_other_mixed,
	                                          h_trk_vtx_dist_unmatched, h_trk_vtx_dist_intime, intime_scale_factor,
	                                          h_trk_vtx_dist_data, data_scale_factor, "",
	                                          "Track to Nue Candidate Vertex Distance [cm]", "", "post_cuts_track_to_vtx_data.pdf");

	histogram_functions::PlotSimpleStackData (h_trk_vtx_dist_nue_cc_after,  h_trk_vtx_dist_nue_cc_mixed_after,
	                                          h_trk_vtx_dist_nue_cc_out_fv_after,
	                                          h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed_after,
	                                          h_trk_vtx_dist_cosmic_after,  h_trk_vtx_dist_nc_after,
	                                          h_trk_vtx_dist_nc_pi0_after,  h_trk_vtx_dist_other_mixed_after,
	                                          h_trk_vtx_dist_unmatched_after, h_trk_vtx_dist_intime_after, intime_scale_factor,
	                                          h_trk_vtx_dist_data_after, data_scale_factor, "",
	                                          "Track to Nue Candidate Vertex Distance [cm]", "", "post_cuts_track_to_vtx_data_after.pdf");

	histogram_functions::PlotSimpleStack (h_shwr_vtx_dist_nue_cc,  h_shwr_vtx_dist_nue_cc_mixed,
	                                      h_shwr_vtx_dist_nue_cc_out_fv,
	                                      h_shwr_vtx_dist_numu_cc, h_shwr_vtx_dist_numu_cc_mixed,
	                                      h_shwr_vtx_dist_cosmic,  h_shwr_vtx_dist_nc,
	                                      h_shwr_vtx_dist_nc_pi0,  h_shwr_vtx_dist_other_mixed,
	                                      h_shwr_vtx_dist_unmatched, "",
	                                      "Shower to Nue Candidate Vertex Distance [cm]", "", "post_cuts_shower_to_vtx.pdf");

	histogram_functions::PlotSimpleStackInTime (h_shwr_vtx_dist_nue_cc,  h_shwr_vtx_dist_nue_cc_mixed,
	                                            h_shwr_vtx_dist_nue_cc_out_fv,
	                                            h_shwr_vtx_dist_numu_cc, h_shwr_vtx_dist_numu_cc_mixed,
	                                            h_shwr_vtx_dist_cosmic,  h_shwr_vtx_dist_nc,
	                                            h_shwr_vtx_dist_nc_pi0,  h_shwr_vtx_dist_other_mixed,
	                                            h_shwr_vtx_dist_unmatched, h_shwr_vtx_dist_intime, intime_scale_factor, "",
	                                            "Shower to Nue Candidate Vertex Distance [cm]", "", "post_cuts_shower_to_vtx_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_shwr_vtx_dist_nue_cc,  h_shwr_vtx_dist_nue_cc_mixed,
	                                          h_shwr_vtx_dist_nue_cc_out_fv,
	                                          h_shwr_vtx_dist_numu_cc, h_shwr_vtx_dist_numu_cc_mixed,
	                                          h_shwr_vtx_dist_cosmic,  h_shwr_vtx_dist_nc,
	                                          h_shwr_vtx_dist_nc_pi0,  h_shwr_vtx_dist_other_mixed,
	                                          h_shwr_vtx_dist_unmatched, h_shwr_vtx_dist_intime, intime_scale_factor,
	                                          h_shwr_vtx_dist_data, data_scale_factor, "",
	                                          "Shower to Nue Candidate Vertex Distance [cm]", "", "post_cuts_shower_to_vtx_data.pdf");

	histogram_functions::PlotSimpleStackData (h_shwr_vtx_dist_nue_cc_after,  h_shwr_vtx_dist_nue_cc_mixed_after,
	                                          h_shwr_vtx_dist_nue_cc_out_fv_after,
	                                          h_shwr_vtx_dist_numu_cc_after, h_shwr_vtx_dist_numu_cc_mixed_after,
	                                          h_shwr_vtx_dist_cosmic_after,  h_shwr_vtx_dist_nc_after,
	                                          h_shwr_vtx_dist_nc_pi0_after,  h_shwr_vtx_dist_other_mixed_after,
	                                          h_shwr_vtx_dist_unmatched_after, h_shwr_vtx_dist_intime_after, intime_scale_factor,
	                                          h_shwr_vtx_dist_data_after, data_scale_factor, "",
	                                          "Shower to Nue Candidate Vertex Distance [cm]", "", "post_cuts_shower_to_vtx_data_after.pdf");

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
	                                     "Reconstructed Tracks in Candidate Neutrino Object", "", "selected_pfp_track_stack.pdf");

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
	                                     "Reconstructed Showers in Candidate Neutrino Object", "", "selected_pfp_shower_stack.pdf");

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
	                                     "Reconstructed Showers in Candidate Neutrino Object", "", "selected_pfp_shower_open_angle_stack.pdf");

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
	                                     "Reconstructed Showers in Candidate Neutrino Object", "", "selected_pfp_shower_dedx_stack.pdf");

	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_qe, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_qe.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_out_fv, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_res, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_res.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_dis, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_dis.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_coh, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_coh.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mec, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_mec.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mixed, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_qe, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_qe.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_res, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_res.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_dis, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_dis.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_coh, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_coh.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mec, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_mec.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mixed, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nc.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc_pi0, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nc_pi0.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_cosmic, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_cosmic.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_other_mixed, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_other_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_unmatched, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_unmatched.pdf");


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
	                                     "Reconstructed Tracks in Candidate Neutrino Object", "", "selected_pfp_track_stack_last.pdf");

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
	                                     "Reconstructed Showers in Candidate Neutrino Object", "", "selected_pfp_shower_stack_last.pdf");

	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_qe_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_qe_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_out_fv_last, "", "PFP Tracks", "PFP Showers",
	                                      "selected_pfp_track_shower_nue_cc_out_fv_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_res_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_res_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_dis_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_dis_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_coh_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_coh_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mec_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_mec_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mixed_last, "", "PFP Tracks", "PFP Showers",
	                                      "selected_pfp_track_shower_nue_cc_mixed_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_qe_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_qe_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_res_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_res_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_dis_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_dis_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_coh_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_coh_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mec_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_mec_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mixed_last, "", "PFP Tracks", "PFP Showers",
	                                      "selected_pfp_track_shower_numu_cc_mixed_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nc_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc_pi0_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nc_pi0_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_cosmic_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_cosmic_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_other_mixed_last, "", "PFP Tracks", "PFP Showers",
	                                      "selected_pfp_track_shower_other_mixed_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_unmatched_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_unmatched_last.pdf");

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
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_nue_cc_qe.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_out_fv, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_res, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_nue_cc_res.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_dis, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_nue_cc_dis.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_coh, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_nue_cc_coh.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_mec, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_nue_cc_mec.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nue_cc_mixed, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_nue_cc_mixed.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_qe, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_numu_cc_qe.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_res, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_numu_cc_res.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_dis, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_numu_cc_dis.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_coh, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_numu_cc_coh.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_mec, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_numu_cc_mec.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_numu_cc_mixed, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_numu_cc_mixed.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nc, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_nc.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_nc_pi0, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_nc_pi0.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_cosmic, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_cosmic.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_other_mixed, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdf_other_mixed.pdf");
	histogram_functions::Plot2DHistogramOffSet (h_leading_shower_mc_pdg_unmatched, 0.002, 1.35, "",
	                                            "Leading Shower Origin", "Leading Shower True Particle", "selected_leading_shower_mc_pdg_unmatched.pdf");

	histogram_functions::Plot2DHistogram (h_shwr_hits_nu_eng, "", "True Neutrino Energy [GeV]", "Signal Electron Shower Hits", "shwr_hits_nu_eng.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_hits_ele_eng, "", "True Electron Energy [GeV]", "Signal Electron Shower Hits", "shwr_hits_ele_eng.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_hits_nu_eng_zoom, "", "True Neutrino Energy [GeV]", "Signal Electron Shower Hits", "shwr_hits_nu_eng_zoom.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_hits_ele_eng_zoom, "", "True Electron Energy [GeV]", "Signal Electron Shower Hits", "shwr_hits_ele_eng_zoom.pdf");

	histogram_functions::Plot2DHistogram (h_shwr_hits_nu_eng_last, "", "True Neutrino Energy [GeV]", "Signal Electron Shower Hits", "shwr_hits_nu_eng_last.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_hits_ele_eng_last, "", "True Electron Energy [GeV]", "Signal Electron Shower Hits", "shwr_hits_ele_eng_last.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_hits_nu_eng_zoom_last, "", "True Neutrino Energy [GeV]",
	                                      "Signal Electron Shower Hits", "shwr_hits_nu_eng_zoom_last.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_hits_ele_eng_zoom_last, "", "True Electron Energy [GeV]",
	                                      "Signal Electron Shower Hits", "shwr_hits_ele_eng_zoom_last.pdf");

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
	                                          "True Signal Neutrino Energy [GeV]", "", "sequential_nu_energy.pdf");
	histogram_functions::PostHistogramOverlay(h_selected_ele_energy_no_cut,       h_selected_ele_energy_reco_nue,
	                                          h_selected_ele_energy_in_fv,        h_selected_ele_energy_vtx_flash,
	                                          h_selected_ele_energy_shwr_vtx,     h_selected_ele_energy_trk_vtx,
	                                          h_selected_ele_energy_hit_threshold,h_selected_ele_energy_open_angle,
	                                          h_selected_ele_energy_dedx,
	                                          "True Signal Electron Energy [GeV]", "", "sequential_ele_energy.pdf");


	histogram_functions::Plot1DHistogram (h_charge_share_nue_cc_mixed, "Neutrino Charge Fraction - Selected Nue CC Mixed", "charge_fraction_nue_cc_mixed.pdf");
	histogram_functions::Plot1DHistogram (h_flash_t0_diff, "Largest Flash Time - True Neutrino Interaction Time [us]", "flash_t0_diff.pdf");

	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nue_cc, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_nue_cc.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nue_cc_out_fv, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nue_cc_mixed, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_nue_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nue_cc, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_nue_cc.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_numu_cc, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_numu_cc.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_numu_cc_mixed, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_numu_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nc, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_nc.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_nc_pi0, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_nc_pi0.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_cosmic, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_cosmic.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_other_mixed, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_other_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_dedx_open_angle_unmatched, "", "Leading Shower dEdx [MeV/cm]",
	                                      "Leading Shower Open Angle [Degrees]", "dedx_open_angle_unmatched.pdf");

	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nue_cc, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_nue_cc.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nue_cc_out_fv, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nue_cc_mixed, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_nue_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_numu_cc, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_numu_cc.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_numu_cc_mixed, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_numu_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nc, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_nc.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nc_pi0, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_nc_pi0.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_cosmic, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_cosmic.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_other_mixed, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_other_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_unmatched, "", "Leading Shower Length [cm]",
	                                      "Leading Shower Hits", "shwr_len_hits_unmatched.pdf");

	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nue_cc, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_nue_cc.pdf");
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nue_cc_out_fv, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_nue_cc_out_fv.pdf");
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nue_cc_mixed, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_nue_cc_mixed.pdf");
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_numu_cc, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_numu_cc.pdf");
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_numu_cc_mixed, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_numu_cc_mixed.pdf");
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nc, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_nc.pdf");
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_nc_pi0, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_nc_pi0.pdf");
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_cosmic, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_cosmic.pdf");
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_other_mixed, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_other_mixed.pdf");
	histogram_functions::Plot1DHistogram (h_second_shwr_dist_unmatched, "(TPCO w/ > 3 showers) Shower-Vtx Distance [cm]", "second_shwr_dist_unmatched.pdf");

	histogram_functions::PlotSimpleStack (h_second_shwr_dist_nue_cc,  h_second_shwr_dist_nue_cc_mixed,
	                                      h_second_shwr_dist_nue_cc_out_fv,
	                                      h_second_shwr_dist_numu_cc, h_second_shwr_dist_numu_cc_mixed,
	                                      h_second_shwr_dist_cosmic,  h_second_shwr_dist_nc,
	                                      h_second_shwr_dist_nc_pi0,  h_second_shwr_dist_other_mixed,
	                                      h_second_shwr_dist_unmatched, "",
	                                      "(TPCO > 1 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "", "post_second_shwr_dist.pdf");

	histogram_functions::PlotSimpleStackInTime (h_second_shwr_dist_nue_cc,  h_second_shwr_dist_nue_cc_mixed,
	                                            h_second_shwr_dist_nue_cc_out_fv,
	                                            h_second_shwr_dist_numu_cc, h_second_shwr_dist_numu_cc_mixed,
	                                            h_second_shwr_dist_cosmic,  h_second_shwr_dist_nc,
	                                            h_second_shwr_dist_nc_pi0,  h_second_shwr_dist_other_mixed,
	                                            h_second_shwr_dist_unmatched, h_second_shwr_dist_intime, intime_scale_factor, "",
	                                            "(TPCO > 1 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "", "post_second_shwr_dist_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_second_shwr_dist_nue_cc,  h_second_shwr_dist_nue_cc_mixed,
	                                          h_second_shwr_dist_nue_cc_out_fv,
	                                          h_second_shwr_dist_numu_cc, h_second_shwr_dist_numu_cc_mixed,
	                                          h_second_shwr_dist_cosmic,  h_second_shwr_dist_nc,
	                                          h_second_shwr_dist_nc_pi0,  h_second_shwr_dist_other_mixed,
	                                          h_second_shwr_dist_unmatched, h_second_shwr_dist_intime, intime_scale_factor,
	                                          h_second_shwr_dist_data, data_scale_factor, "",
	                                          "(TPCO > 1 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "", "post_second_shwr_dist_data.pdf");

	histogram_functions::PlotSimpleStackData (h_second_shwr_dist_nue_cc_after,  h_second_shwr_dist_nue_cc_mixed_after,
	                                          h_second_shwr_dist_nue_cc_out_fv_after,
	                                          h_second_shwr_dist_numu_cc_after, h_second_shwr_dist_numu_cc_mixed_after,
	                                          h_second_shwr_dist_cosmic_after,  h_second_shwr_dist_nc_after,
	                                          h_second_shwr_dist_nc_pi0_after,  h_second_shwr_dist_other_mixed_after,
	                                          h_second_shwr_dist_unmatched_after, h_second_shwr_dist_intime_after, intime_scale_factor,
	                                          h_second_shwr_dist_data_after, data_scale_factor, "",
	                                          "(TPCO > 3 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "", "post_second_shwr_dist_data_after.pdf");

	histogram_functions::PlotSimpleStack (h_hit_length_ratio_nue_cc,  h_hit_length_ratio_nue_cc_mixed,
	                                      h_hit_length_ratio_nue_cc_out_fv,
	                                      h_hit_length_ratio_numu_cc, h_hit_length_ratio_numu_cc_mixed,
	                                      h_hit_length_ratio_cosmic,  h_hit_length_ratio_nc,
	                                      h_hit_length_ratio_nc_pi0,  h_hit_length_ratio_other_mixed,
	                                      h_hit_length_ratio_unmatched, "",
	                                      "Leading Shower (Hits / Length) [cm^-1]", "", "post_hit_length_ratio.pdf");

	histogram_functions::PlotSimpleStackInTime (h_hit_length_ratio_nue_cc,  h_hit_length_ratio_nue_cc_mixed,
	                                            h_hit_length_ratio_nue_cc_out_fv,
	                                            h_hit_length_ratio_numu_cc, h_hit_length_ratio_numu_cc_mixed,
	                                            h_hit_length_ratio_cosmic,  h_hit_length_ratio_nc,
	                                            h_hit_length_ratio_nc_pi0,  h_hit_length_ratio_other_mixed,
	                                            h_hit_length_ratio_unmatched, h_hit_length_ratio_intime, intime_scale_factor, "",
	                                            "Leading Shower (Hits / Length) [cm^-1]", "", "post_hit_length_ratio_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_hit_length_ratio_nue_cc,  h_hit_length_ratio_nue_cc_mixed,
	                                          h_hit_length_ratio_nue_cc_out_fv,
	                                          h_hit_length_ratio_numu_cc, h_hit_length_ratio_numu_cc_mixed,
	                                          h_hit_length_ratio_cosmic,  h_hit_length_ratio_nc,
	                                          h_hit_length_ratio_nc_pi0,  h_hit_length_ratio_other_mixed,
	                                          h_hit_length_ratio_unmatched, h_hit_length_ratio_intime, intime_scale_factor,
	                                          h_hit_length_ratio_data, data_scale_factor, "",
	                                          "Leading Shower (Hits / Length) [cm^-1]", "", "post_hit_length_ratio_data.pdf");

	histogram_functions::PlotSimpleStackData (h_hit_length_ratio_nue_cc_after,  h_hit_length_ratio_nue_cc_mixed_after,
	                                          h_hit_length_ratio_nue_cc_out_fv_after,
	                                          h_hit_length_ratio_numu_cc_after, h_hit_length_ratio_numu_cc_mixed_after,
	                                          h_hit_length_ratio_cosmic_after,  h_hit_length_ratio_nc_after,
	                                          h_hit_length_ratio_nc_pi0_after,  h_hit_length_ratio_other_mixed_after,
	                                          h_hit_length_ratio_unmatched_after, h_hit_length_ratio_intime_after, intime_scale_factor,
	                                          h_hit_length_ratio_data_after, data_scale_factor, "",
	                                          "Leading Shower (Hits / Length) [cm^-1]", "", "post_hit_length_ratio_data_after.pdf");

	histogram_functions::PlotSimpleStack (h_trk_length_nue_cc,  h_trk_length_nue_cc_mixed,
	                                      h_trk_length_nue_cc_out_fv,
	                                      h_trk_length_numu_cc, h_trk_length_numu_cc_mixed,
	                                      h_trk_length_cosmic,  h_trk_length_nc,
	                                      h_trk_length_nc_pi0,  h_trk_length_other_mixed,
	                                      h_trk_length_unmatched, "",
	                                      "All Track Lengths [cm]", "", "post_all_track_lengths.pdf");

	histogram_functions::PlotSimpleStackInTime (h_trk_length_nue_cc,  h_trk_length_nue_cc_mixed,
	                                            h_trk_length_nue_cc_out_fv,
	                                            h_trk_length_numu_cc, h_trk_length_numu_cc_mixed,
	                                            h_trk_length_cosmic,  h_trk_length_nc,
	                                            h_trk_length_nc_pi0,  h_trk_length_other_mixed,
	                                            h_trk_length_unmatched, h_trk_length_intime, intime_scale_factor, "",
	                                            "All Track Lengths [cm]", "", "post_all_track_lengths_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_trk_length_nue_cc,  h_trk_length_nue_cc_mixed,
	                                          h_trk_length_nue_cc_out_fv,
	                                          h_trk_length_numu_cc, h_trk_length_numu_cc_mixed,
	                                          h_trk_length_cosmic,  h_trk_length_nc,
	                                          h_trk_length_nc_pi0,  h_trk_length_other_mixed,
	                                          h_trk_length_unmatched, h_trk_length_intime, intime_scale_factor,
	                                          h_trk_length_data, data_scale_factor, "",
	                                          "All Track Lengths [cm]", "", "post_all_track_lengths_data.pdf");

	histogram_functions::PlotSimpleStack (h_longest_trk_length_nue_cc,  h_longest_trk_length_nue_cc_mixed,
	                                      h_longest_trk_length_nue_cc_out_fv,
	                                      h_longest_trk_length_numu_cc, h_longest_trk_length_numu_cc_mixed,
	                                      h_longest_trk_length_cosmic,  h_longest_trk_length_nc,
	                                      h_longest_trk_length_nc_pi0,  h_longest_trk_length_other_mixed,
	                                      h_longest_trk_length_unmatched, "",
	                                      "Longest Track Lengths [cm]", "", "post_longest_track_lengths.pdf");

	histogram_functions::PlotSimpleStackInTime (h_longest_trk_length_nue_cc,  h_longest_trk_length_nue_cc_mixed,
	                                            h_longest_trk_length_nue_cc_out_fv,
	                                            h_longest_trk_length_numu_cc, h_longest_trk_length_numu_cc_mixed,
	                                            h_longest_trk_length_cosmic,  h_longest_trk_length_nc,
	                                            h_longest_trk_length_nc_pi0,  h_longest_trk_length_other_mixed,
	                                            h_longest_trk_length_unmatched, h_longest_trk_length_intime, intime_scale_factor, "",
	                                            "Longest Track Lengths [cm]", "", "post_longest_track_lengths_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_longest_trk_length_nue_cc,  h_longest_trk_length_nue_cc_mixed,
	                                          h_longest_trk_length_nue_cc_out_fv,
	                                          h_longest_trk_length_numu_cc, h_longest_trk_length_numu_cc_mixed,
	                                          h_longest_trk_length_cosmic,  h_longest_trk_length_nc,
	                                          h_longest_trk_length_nc_pi0,  h_longest_trk_length_other_mixed,
	                                          h_longest_trk_length_unmatched, h_longest_trk_length_intime, intime_scale_factor,
	                                          h_longest_trk_length_data, data_scale_factor, "",
	                                          "Longest Track Lengths [cm]", "", "post_longest_track_lengths_data.pdf");

	histogram_functions::PlotSimpleStack (h_shwr_length_nue_cc,  h_shwr_length_nue_cc_mixed,
	                                      h_shwr_length_nue_cc_out_fv,
	                                      h_shwr_length_numu_cc, h_shwr_length_numu_cc_mixed,
	                                      h_shwr_length_cosmic,  h_shwr_length_nc,
	                                      h_shwr_length_nc_pi0,  h_shwr_length_other_mixed,
	                                      h_shwr_length_unmatched, "",
	                                      "All Shower Lengths [cm]", "", "post_all_shower_lengths.pdf");
	histogram_functions::PlotSimpleStack (h_longest_shwr_length_nue_cc,  h_longest_shwr_length_nue_cc_mixed,
	                                      h_longest_shwr_length_nue_cc_out_fv,
	                                      h_longest_shwr_length_numu_cc, h_longest_shwr_length_numu_cc_mixed,
	                                      h_longest_shwr_length_cosmic,  h_longest_shwr_length_nc,
	                                      h_longest_shwr_length_nc_pi0,  h_longest_shwr_length_other_mixed,
	                                      h_longest_shwr_length_unmatched, "",
	                                      "Longest Shower Lengths [cm]", "", "post_longest_shower_lengths.pdf");
	histogram_functions::PlotSimpleStack (h_leading_shwr_length_nue_cc,  h_leading_shwr_length_nue_cc_mixed,
	                                      h_leading_shwr_length_nue_cc_out_fv,
	                                      h_leading_shwr_length_numu_cc, h_leading_shwr_length_numu_cc_mixed,
	                                      h_leading_shwr_length_cosmic,  h_leading_shwr_length_nc,
	                                      h_leading_shwr_length_nc_pi0,  h_leading_shwr_length_other_mixed,
	                                      h_leading_shwr_length_unmatched, "",
	                                      "Leading Shower Lengths [cm]", "", "post_leading_shower_lengths.pdf");

	histogram_functions::PlotSimpleStackInTime (h_shwr_length_nue_cc,  h_shwr_length_nue_cc_mixed,
	                                            h_shwr_length_nue_cc_out_fv,
	                                            h_shwr_length_numu_cc, h_shwr_length_numu_cc_mixed,
	                                            h_shwr_length_cosmic,  h_shwr_length_nc,
	                                            h_shwr_length_nc_pi0,  h_shwr_length_other_mixed,
	                                            h_shwr_length_unmatched, h_shwr_length_intime, intime_scale_factor, "",
	                                            "All Shower Lengths [cm]", "", "post_all_shower_lengths_intime.pdf");
	histogram_functions::PlotSimpleStackInTime (h_longest_shwr_length_nue_cc,  h_longest_shwr_length_nue_cc_mixed,
	                                            h_longest_shwr_length_nue_cc_out_fv,
	                                            h_longest_shwr_length_numu_cc, h_longest_shwr_length_numu_cc_mixed,
	                                            h_longest_shwr_length_cosmic,  h_longest_shwr_length_nc,
	                                            h_longest_shwr_length_nc_pi0,  h_longest_shwr_length_other_mixed,
	                                            h_longest_shwr_length_unmatched, h_longest_shwr_length_intime, intime_scale_factor, "",
	                                            "Longest Shower Lengths [cm]", "", "post_longest_shower_lengths_intime.pdf");
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_length_nue_cc,  h_leading_shwr_length_nue_cc_mixed,
	                                            h_leading_shwr_length_nue_cc_out_fv,
	                                            h_leading_shwr_length_numu_cc, h_leading_shwr_length_numu_cc_mixed,
	                                            h_leading_shwr_length_cosmic,  h_leading_shwr_length_nc,
	                                            h_leading_shwr_length_nc_pi0,  h_leading_shwr_length_other_mixed,
	                                            h_leading_shwr_length_unmatched, h_leading_shwr_length_intime, intime_scale_factor, "",
	                                            "Leading Shower Lengths [cm]", "", "post_leading_shower_lengths_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_shwr_length_nue_cc,  h_shwr_length_nue_cc_mixed,
	                                          h_shwr_length_nue_cc_out_fv,
	                                          h_shwr_length_numu_cc, h_shwr_length_numu_cc_mixed,
	                                          h_shwr_length_cosmic,  h_shwr_length_nc,
	                                          h_shwr_length_nc_pi0,  h_shwr_length_other_mixed,
	                                          h_shwr_length_unmatched, h_shwr_length_intime, intime_scale_factor,
	                                          h_shwr_length_data, data_scale_factor, "",
	                                          "All Shower Lengths [cm]", "", "post_all_shower_lengths_data.pdf");
	histogram_functions::PlotSimpleStackData (h_longest_shwr_length_nue_cc,  h_longest_shwr_length_nue_cc_mixed,
	                                          h_longest_shwr_length_nue_cc_out_fv,
	                                          h_longest_shwr_length_numu_cc, h_longest_shwr_length_numu_cc_mixed,
	                                          h_longest_shwr_length_cosmic,  h_longest_shwr_length_nc,
	                                          h_longest_shwr_length_nc_pi0,  h_longest_shwr_length_other_mixed,
	                                          h_longest_shwr_length_unmatched, h_longest_shwr_length_intime, intime_scale_factor,
	                                          h_longest_shwr_length_data, data_scale_factor, "",
	                                          "Longest Shower Lengths [cm]", "", "post_longest_shower_lengths_data.pdf");
	histogram_functions::PlotSimpleStackData (h_leading_shwr_length_nue_cc,  h_leading_shwr_length_nue_cc_mixed,
	                                          h_leading_shwr_length_nue_cc_out_fv,
	                                          h_leading_shwr_length_numu_cc, h_leading_shwr_length_numu_cc_mixed,
	                                          h_leading_shwr_length_cosmic,  h_leading_shwr_length_nc,
	                                          h_leading_shwr_length_nc_pi0,  h_leading_shwr_length_other_mixed,
	                                          h_leading_shwr_length_unmatched, h_leading_shwr_length_intime, intime_scale_factor,
	                                          h_leading_shwr_length_data, data_scale_factor, "",
	                                          "Leading Shower Lengths [cm]", "", "post_leading_shower_lengths_data.pdf");

	histogram_functions::PlotSimpleStack (h_leading_shwr_trk_length_nue_cc,  h_leading_shwr_trk_length_nue_cc_mixed,
	                                      h_leading_shwr_trk_length_nue_cc_out_fv,
	                                      h_leading_shwr_trk_length_numu_cc, h_leading_shwr_trk_length_numu_cc_mixed,
	                                      h_leading_shwr_trk_length_cosmic,  h_leading_shwr_trk_length_nc,
	                                      h_leading_shwr_trk_length_nc_pi0,  h_leading_shwr_trk_length_other_mixed,
	                                      h_leading_shwr_trk_length_unmatched, "",
	                                      "Longest Track / Leading Shower Lengths", "", "post_leading_shower_trk_lengths.pdf");

	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_trk_length_nue_cc,  h_leading_shwr_trk_length_nue_cc_mixed,
	                                            h_leading_shwr_trk_length_nue_cc_out_fv,
	                                            h_leading_shwr_trk_length_numu_cc, h_leading_shwr_trk_length_numu_cc_mixed,
	                                            h_leading_shwr_trk_length_cosmic,  h_leading_shwr_trk_length_nc,
	                                            h_leading_shwr_trk_length_nc_pi0,  h_leading_shwr_trk_length_other_mixed,
	                                            h_leading_shwr_trk_length_unmatched, h_leading_shwr_trk_length_intime, intime_scale_factor, "",
	                                            "Longest Track / Leading Shower Lengths", "", "post_leading_shower_trk_lengths_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_leading_shwr_trk_length_nue_cc,  h_leading_shwr_trk_length_nue_cc_mixed,
	                                          h_leading_shwr_trk_length_nue_cc_out_fv,
	                                          h_leading_shwr_trk_length_numu_cc, h_leading_shwr_trk_length_numu_cc_mixed,
	                                          h_leading_shwr_trk_length_cosmic,  h_leading_shwr_trk_length_nc,
	                                          h_leading_shwr_trk_length_nc_pi0,  h_leading_shwr_trk_length_other_mixed,
	                                          h_leading_shwr_trk_length_unmatched, h_leading_shwr_trk_length_intime, intime_scale_factor,
	                                          h_leading_shwr_trk_length_data, data_scale_factor, "",
	                                          "Longest Track / Leading Shower Lengths", "", "post_leading_shower_trk_lengths_data.pdf");

	histogram_functions::PlotSimpleStackData (h_leading_shwr_trk_length_nue_cc_after,  h_leading_shwr_trk_length_nue_cc_mixed_after,
	                                          h_leading_shwr_trk_length_nue_cc_out_fv_after,
	                                          h_leading_shwr_trk_length_numu_cc_after, h_leading_shwr_trk_length_numu_cc_mixed_after,
	                                          h_leading_shwr_trk_length_cosmic_after,  h_leading_shwr_trk_length_nc_after,
	                                          h_leading_shwr_trk_length_nc_pi0_after,  h_leading_shwr_trk_length_other_mixed_after,
	                                          h_leading_shwr_trk_length_unmatched_after, h_leading_shwr_trk_length_intime_after, intime_scale_factor,
	                                          h_leading_shwr_trk_length_data_after, data_scale_factor, "",
	                                          "Longest Track / Leading Shower Lengths", "", "post_leading_shower_trk_lengths_data_after.pdf");

	histogram_functions::PlotSimpleStack (h_longest_shwr_trk_length_nue_cc,  h_longest_shwr_trk_length_nue_cc_mixed,
	                                      h_longest_shwr_trk_length_nue_cc_out_fv,
	                                      h_longest_shwr_trk_length_numu_cc, h_longest_shwr_trk_length_numu_cc_mixed,
	                                      h_longest_shwr_trk_length_cosmic,  h_longest_shwr_trk_length_nc,
	                                      h_longest_shwr_trk_length_nc_pi0,  h_longest_shwr_trk_length_other_mixed,
	                                      h_longest_shwr_trk_length_unmatched, "",
	                                      "Longest Track / Longest Shower Lengths", "", "post_longest_shower_trk_lengths.pdf");

	histogram_functions::PlotSimpleStackInTime (h_longest_shwr_trk_length_nue_cc,  h_longest_shwr_trk_length_nue_cc_mixed,
	                                            h_longest_shwr_trk_length_nue_cc_out_fv,
	                                            h_longest_shwr_trk_length_numu_cc, h_longest_shwr_trk_length_numu_cc_mixed,
	                                            h_longest_shwr_trk_length_cosmic,  h_longest_shwr_trk_length_nc,
	                                            h_longest_shwr_trk_length_nc_pi0,  h_longest_shwr_trk_length_other_mixed,
	                                            h_longest_shwr_trk_length_unmatched, h_longest_shwr_trk_length_intime, intime_scale_factor, "",
	                                            "Longest Track / Longest Shower Lengths", "", "post_longest_shower_trk_lengths_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_longest_shwr_trk_length_nue_cc,  h_longest_shwr_trk_length_nue_cc_mixed,
	                                          h_longest_shwr_trk_length_nue_cc_out_fv,
	                                          h_longest_shwr_trk_length_numu_cc, h_longest_shwr_trk_length_numu_cc_mixed,
	                                          h_longest_shwr_trk_length_cosmic,  h_longest_shwr_trk_length_nc,
	                                          h_longest_shwr_trk_length_nc_pi0,  h_longest_shwr_trk_length_other_mixed,
	                                          h_longest_shwr_trk_length_unmatched, h_longest_shwr_trk_length_intime, intime_scale_factor,
	                                          h_longest_shwr_trk_length_data, data_scale_factor, "",
	                                          "Longest Track / Longest Shower Lengths", "", "post_longest_shower_trk_lengths_data.pdf");

	// histogram_functions::Plot1DHistogram (h_post_cuts_num_showers_purity, "Number of Reconstructed Showers per TPCO", "post_cuts_num_showers_purity.pdf");
	// histogram_functions::Plot1DHistogram (h_post_open_angle_cuts_num_showers_purity, "Number of Reconstructed Showers per TPCO - After Open Angle Cut",
	//                                       "post_open_angle_cuts_num_showers_purity.pdf");

	histogram_functions::PurityStack(h_post_cuts_num_showers_purity_qe, h_post_cuts_num_showers_purity_res, h_post_cuts_num_showers_purity_dis,
	                                 h_post_cuts_num_showers_purity_coh, h_post_cuts_num_showers_purity_mec,
	                                 "Num Reco Showers per TPCO", "post_cuts_num_showers_purity.pdf");
	histogram_functions::PurityStack(h_post_open_angle_cuts_num_showers_purity_qe, h_post_open_angle_cuts_num_showers_purity_res,
	                                 h_post_open_angle_cuts_num_showers_purity_dis, h_post_open_angle_cuts_num_showers_purity_coh,
	                                 h_post_open_angle_cuts_num_showers_purity_mec,
	                                 "Num Reco Showers per TPCO - After Open Angle Cut", "post_open_angle_cuts_num_showers_purity.pdf");

	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nue_cc, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_nue_cc.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nue_cc_out_fv, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nue_cc_mixed, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_nue_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_numu_cc, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_numu_cc.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_numu_cc_mixed, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_numu_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nc, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_nc.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_nc_pi0, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_nc_pi0.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_cosmic, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_cosmic.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_other_mixed, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_other_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_track_unmatched, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_track_unmatched.pdf");

	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nue_cc, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_nue_cc.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nue_cc_out_fv, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nue_cc_mixed, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_nue_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_numu_cc, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_numu_cc.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_numu_cc_mixed, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_numu_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nc, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_nc.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_nc_pi0, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_nc_pi0.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_cosmic, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_cosmic.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_other_mixed, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_other_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_shower_unmatched, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_shower_unmatched.pdf");

	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nue_cc, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_nue_cc.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nue_cc_out_fv, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nue_cc_mixed, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_nue_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_numu_cc, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_numu_cc.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_numu_cc_mixed, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_numu_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nc, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_nc.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_nc_pi0, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_nc_pi0.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_cosmic, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_cosmic.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_other_mixed, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_other_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_collection_total_hits_leading_shower_unmatched, "",
	                                      "Hits - Collection Plane", "Hits - All Planes", "post_cuts_collection_total_hits_leading_shower_unmatched.pdf");

	histogram_functions::PlotSimpleStack (h_collection_hits_track_nue_cc,  h_collection_hits_track_nue_cc_mixed,
	                                      h_collection_hits_track_nue_cc_out_fv,
	                                      h_collection_hits_track_numu_cc, h_collection_hits_track_numu_cc_mixed,
	                                      h_collection_hits_track_cosmic,  h_collection_hits_track_nc,
	                                      h_collection_hits_track_nc_pi0,  h_collection_hits_track_other_mixed,
	                                      h_collection_hits_track_unmatched, "",
	                                      "Tracks Hits - Collection Plane", "", "post_cuts_collection_hits_tracks.pdf");
	histogram_functions::PlotSimpleStack (h_collection_hits_shower_nue_cc,  h_collection_hits_shower_nue_cc_mixed,
	                                      h_collection_hits_shower_nue_cc_out_fv,
	                                      h_collection_hits_shower_numu_cc, h_collection_hits_shower_numu_cc_mixed,
	                                      h_collection_hits_shower_cosmic,  h_collection_hits_shower_nc,
	                                      h_collection_hits_shower_nc_pi0,  h_collection_hits_shower_other_mixed,
	                                      h_collection_hits_shower_unmatched, "",
	                                      "Showers Hits - Collection Plane", "", "post_cuts_collection_hits_showers.pdf");
	histogram_functions::PlotSimpleStack (h_collection_hits_leading_shower_nue_cc,  h_collection_hits_leading_shower_nue_cc_mixed,
	                                      h_collection_hits_leading_shower_nue_cc_out_fv,
	                                      h_collection_hits_leading_shower_numu_cc, h_collection_hits_leading_shower_numu_cc_mixed,
	                                      h_collection_hits_leading_shower_cosmic,  h_collection_hits_leading_shower_nc,
	                                      h_collection_hits_leading_shower_nc_pi0,  h_collection_hits_leading_shower_other_mixed,
	                                      h_collection_hits_leading_shower_unmatched, "",
	                                      "Leading Shower Hits - Collection Plane", "", "post_cuts_collection_hits_leading_shower.pdf");
	histogram_functions::PlotSimpleStack (h_total_hits_leading_shower_nue_cc,  h_total_hits_leading_shower_nue_cc_mixed,
	                                      h_total_hits_leading_shower_nue_cc_out_fv,
	                                      h_total_hits_leading_shower_numu_cc, h_total_hits_leading_shower_numu_cc_mixed,
	                                      h_total_hits_leading_shower_cosmic,  h_total_hits_leading_shower_nc,
	                                      h_total_hits_leading_shower_nc_pi0,  h_total_hits_leading_shower_other_mixed,
	                                      h_total_hits_leading_shower_unmatched, "",
	                                      "Leading Shower Hits - All Planes", "", "post_cuts_total_hits_leading_shower.pdf");

	histogram_functions::PlotSimpleStackInTime (h_collection_hits_track_nue_cc,  h_collection_hits_track_nue_cc_mixed,
	                                            h_collection_hits_track_nue_cc_out_fv,
	                                            h_collection_hits_track_numu_cc, h_collection_hits_track_numu_cc_mixed,
	                                            h_collection_hits_track_cosmic,  h_collection_hits_track_nc,
	                                            h_collection_hits_track_nc_pi0,  h_collection_hits_track_other_mixed,
	                                            h_collection_hits_track_unmatched, h_collection_hits_track_intime, intime_scale_factor, "",
	                                            "Tracks Hits - Collection Plane", "", "post_cuts_collection_hits_tracks_intime.pdf");
	histogram_functions::PlotSimpleStackInTime (h_collection_hits_shower_nue_cc,  h_collection_hits_shower_nue_cc_mixed,
	                                            h_collection_hits_shower_nue_cc_out_fv,
	                                            h_collection_hits_shower_numu_cc, h_collection_hits_shower_numu_cc_mixed,
	                                            h_collection_hits_shower_cosmic,  h_collection_hits_shower_nc,
	                                            h_collection_hits_shower_nc_pi0,  h_collection_hits_shower_other_mixed,
	                                            h_collection_hits_shower_unmatched, h_collection_hits_shower_intime, intime_scale_factor, "",
	                                            "Showers Hits - Collection Plane", "", "post_cuts_collection_hits_showers_intime.pdf");
	histogram_functions::PlotSimpleStackInTime (h_collection_hits_leading_shower_nue_cc,  h_collection_hits_leading_shower_nue_cc_mixed,
	                                            h_collection_hits_leading_shower_nue_cc_out_fv,
	                                            h_collection_hits_leading_shower_numu_cc, h_collection_hits_leading_shower_numu_cc_mixed,
	                                            h_collection_hits_leading_shower_cosmic,  h_collection_hits_leading_shower_nc,
	                                            h_collection_hits_leading_shower_nc_pi0,  h_collection_hits_leading_shower_other_mixed,
	                                            h_collection_hits_leading_shower_unmatched, h_collection_hits_leading_shower_intime, intime_scale_factor,
	                                            "", "Leading Shower Hits - Collection Plane", "", "post_cuts_collection_hits_leading_shower_intime.pdf");
	histogram_functions::PlotSimpleStackInTime (h_total_hits_leading_shower_nue_cc,  h_total_hits_leading_shower_nue_cc_mixed,
	                                            h_total_hits_leading_shower_nue_cc_out_fv,
	                                            h_total_hits_leading_shower_numu_cc, h_total_hits_leading_shower_numu_cc_mixed,
	                                            h_total_hits_leading_shower_cosmic,  h_total_hits_leading_shower_nc,
	                                            h_total_hits_leading_shower_nc_pi0,  h_total_hits_leading_shower_other_mixed,
	                                            h_total_hits_leading_shower_unmatched, h_total_hits_leading_shower_intime, intime_scale_factor,
	                                            "", "Leading Shower Hits - All Planes", "", "post_cuts_total_hits_leading_shower_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_collection_hits_track_nue_cc,  h_collection_hits_track_nue_cc_mixed,
	                                          h_collection_hits_track_nue_cc_out_fv,
	                                          h_collection_hits_track_numu_cc, h_collection_hits_track_numu_cc_mixed,
	                                          h_collection_hits_track_cosmic,  h_collection_hits_track_nc,
	                                          h_collection_hits_track_nc_pi0,  h_collection_hits_track_other_mixed,
	                                          h_collection_hits_track_unmatched, h_collection_hits_track_intime, intime_scale_factor,
	                                          h_collection_hits_track_data, data_scale_factor, "",
	                                          "Tracks Hits - Collection Plane", "", "post_cuts_collection_hits_tracks_data.pdf");
	histogram_functions::PlotSimpleStackData (h_collection_hits_shower_nue_cc,  h_collection_hits_shower_nue_cc_mixed,
	                                          h_collection_hits_shower_nue_cc_out_fv,
	                                          h_collection_hits_shower_numu_cc, h_collection_hits_shower_numu_cc_mixed,
	                                          h_collection_hits_shower_cosmic,  h_collection_hits_shower_nc,
	                                          h_collection_hits_shower_nc_pi0,  h_collection_hits_shower_other_mixed,
	                                          h_collection_hits_shower_unmatched, h_collection_hits_shower_intime, intime_scale_factor,
	                                          h_collection_hits_shower_data, data_scale_factor, "",
	                                          "Showers Hits - Collection Plane", "", "post_cuts_collection_hits_showers_data.pdf");
	histogram_functions::PlotSimpleStackData (h_collection_hits_leading_shower_nue_cc,  h_collection_hits_leading_shower_nue_cc_mixed,
	                                          h_collection_hits_leading_shower_nue_cc_out_fv,
	                                          h_collection_hits_leading_shower_numu_cc, h_collection_hits_leading_shower_numu_cc_mixed,
	                                          h_collection_hits_leading_shower_cosmic,  h_collection_hits_leading_shower_nc,
	                                          h_collection_hits_leading_shower_nc_pi0,  h_collection_hits_leading_shower_other_mixed,
	                                          h_collection_hits_leading_shower_unmatched,
	                                          h_collection_hits_leading_shower_intime, intime_scale_factor,
	                                          h_collection_hits_leading_shower_data, data_scale_factor,
	                                          "", "Leading Shower Hits - Collection Plane", "", "post_cuts_collection_hits_leading_shower_data.pdf");
	histogram_functions::PlotSimpleStackData (h_total_hits_leading_shower_nue_cc,  h_total_hits_leading_shower_nue_cc_mixed,
	                                          h_total_hits_leading_shower_nue_cc_out_fv,
	                                          h_total_hits_leading_shower_numu_cc, h_total_hits_leading_shower_numu_cc_mixed,
	                                          h_total_hits_leading_shower_cosmic,  h_total_hits_leading_shower_nc,
	                                          h_total_hits_leading_shower_nc_pi0,  h_total_hits_leading_shower_other_mixed,
	                                          h_total_hits_leading_shower_unmatched, h_total_hits_leading_shower_intime, intime_scale_factor,
	                                          h_total_hits_leading_shower_data, data_scale_factor,
	                                          "", "Leading Shower Hits - All Planes", "", "post_cuts_total_hits_leading_shower_data.pdf");

	histogram_functions::PlotSimpleStackData (h_collection_hits_track_nue_cc_after,  h_collection_hits_track_nue_cc_mixed_after,
	                                          h_collection_hits_track_nue_cc_out_fv_after,
	                                          h_collection_hits_track_numu_cc_after, h_collection_hits_track_numu_cc_mixed_after,
	                                          h_collection_hits_track_cosmic_after,  h_collection_hits_track_nc_after,
	                                          h_collection_hits_track_nc_pi0_after,  h_collection_hits_track_other_mixed_after,
	                                          h_collection_hits_track_unmatched_after, h_collection_hits_track_intime_after, intime_scale_factor,
	                                          h_collection_hits_track_data_after, data_scale_factor, "",
	                                          "Tracks Hits - Collection Plane", "", "post_cuts_collection_hits_tracks_data_after.pdf");
	histogram_functions::PlotSimpleStackData (h_collection_hits_shower_nue_cc_after,  h_collection_hits_shower_nue_cc_mixed_after,
	                                          h_collection_hits_shower_nue_cc_out_fv_after,
	                                          h_collection_hits_shower_numu_cc_after, h_collection_hits_shower_numu_cc_mixed_after,
	                                          h_collection_hits_shower_cosmic_after,  h_collection_hits_shower_nc_after,
	                                          h_collection_hits_shower_nc_pi0_after,  h_collection_hits_shower_other_mixed_after,
	                                          h_collection_hits_shower_unmatched_after, h_collection_hits_shower_intime_after, intime_scale_factor,
	                                          h_collection_hits_shower_data_after, data_scale_factor, "",
	                                          "Showers Hits - Collection Plane", "", "post_cuts_collection_hits_showers_data_after.pdf");
	histogram_functions::PlotSimpleStackData (h_collection_hits_leading_shower_nue_cc_after,  h_collection_hits_leading_shower_nue_cc_mixed_after,
	                                          h_collection_hits_leading_shower_nue_cc_out_fv_after,
	                                          h_collection_hits_leading_shower_numu_cc_after, h_collection_hits_leading_shower_numu_cc_mixed_after,
	                                          h_collection_hits_leading_shower_cosmic_after,  h_collection_hits_leading_shower_nc_after,
	                                          h_collection_hits_leading_shower_nc_pi0_after,  h_collection_hits_leading_shower_other_mixed_after,
	                                          h_collection_hits_leading_shower_unmatched_after,
	                                          h_collection_hits_leading_shower_intime_after, intime_scale_factor,
	                                          h_collection_hits_leading_shower_data_after, data_scale_factor,
	                                          "", "Leading Shower Hits - Collection Plane", "", "post_cuts_collection_hits_leading_shower_data_after.pdf");
	histogram_functions::PlotSimpleStackData (h_total_hits_leading_shower_nue_cc_after,  h_total_hits_leading_shower_nue_cc_mixed_after,
	                                          h_total_hits_leading_shower_nue_cc_out_fv_after,
	                                          h_total_hits_leading_shower_numu_cc_after, h_total_hits_leading_shower_numu_cc_mixed_after,
	                                          h_total_hits_leading_shower_cosmic_after,  h_total_hits_leading_shower_nc_after,
	                                          h_total_hits_leading_shower_nc_pi0_after,  h_total_hits_leading_shower_other_mixed_after,
	                                          h_total_hits_leading_shower_unmatched_after, h_total_hits_leading_shower_intime_after, intime_scale_factor,
	                                          h_total_hits_leading_shower_data_after, data_scale_factor,
	                                          "", "Leading Shower Hits - All Planes", "", "post_cuts_total_hits_leading_shower_data_after.pdf");

	histogram_functions::PlotSimpleStack (h_pre_cut_collection_hits_track_nue_cc,  h_pre_cut_collection_hits_track_nue_cc_mixed,
	                                      h_pre_cut_collection_hits_track_nue_cc_out_fv,
	                                      h_pre_cut_collection_hits_track_numu_cc, h_pre_cut_collection_hits_track_numu_cc_mixed,
	                                      h_pre_cut_collection_hits_track_cosmic,  h_pre_cut_collection_hits_track_nc,
	                                      h_pre_cut_collection_hits_track_nc_pi0,  h_pre_cut_collection_hits_track_other_mixed,
	                                      h_pre_cut_collection_hits_track_unmatched, "",
	                                      "Tracks Hits - Collection Plane", "", "pre_hit_cut_collection_hits_tracks.pdf");
	histogram_functions::PlotSimpleStack (h_pre_cut_collection_hits_shower_nue_cc,  h_pre_cut_collection_hits_shower_nue_cc_mixed,
	                                      h_pre_cut_collection_hits_shower_nue_cc_out_fv,
	                                      h_pre_cut_collection_hits_shower_numu_cc, h_pre_cut_collection_hits_shower_numu_cc_mixed,
	                                      h_pre_cut_collection_hits_shower_cosmic,  h_pre_cut_collection_hits_shower_nc,
	                                      h_pre_cut_collection_hits_shower_nc_pi0,  h_pre_cut_collection_hits_shower_other_mixed,
	                                      h_pre_cut_collection_hits_shower_unmatched, "",
	                                      "Showers Hits - Collection Plane", "", "pre_hit_cut_collection_hits_showers.pdf");
	histogram_functions::PlotSimpleStack (h_pre_cut_collection_hits_leading_shower_nue_cc,  h_pre_cut_collection_hits_leading_shower_nue_cc_mixed,
	                                      h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv,
	                                      h_pre_cut_collection_hits_leading_shower_numu_cc, h_pre_cut_collection_hits_leading_shower_numu_cc_mixed,
	                                      h_pre_cut_collection_hits_leading_shower_cosmic,  h_pre_cut_collection_hits_leading_shower_nc,
	                                      h_pre_cut_collection_hits_leading_shower_nc_pi0,  h_pre_cut_collection_hits_leading_shower_other_mixed,
	                                      h_pre_cut_collection_hits_leading_shower_unmatched, "",
	                                      "Leading Shower Hits - Collection Plane", "", "pre_hit_cut_collection_hits_leading_shower.pdf");
	histogram_functions::PlotSimpleStack (h_pre_cut_total_hits_leading_shower_nue_cc,  h_pre_cut_total_hits_leading_shower_nue_cc_mixed,
	                                      h_pre_cut_total_hits_leading_shower_nue_cc_out_fv,
	                                      h_pre_cut_total_hits_leading_shower_numu_cc, h_pre_cut_total_hits_leading_shower_numu_cc_mixed,
	                                      h_pre_cut_total_hits_leading_shower_cosmic,  h_pre_cut_total_hits_leading_shower_nc,
	                                      h_pre_cut_total_hits_leading_shower_nc_pi0,  h_pre_cut_total_hits_leading_shower_other_mixed,
	                                      h_pre_cut_total_hits_leading_shower_unmatched, "",
	                                      "Leading Shower Hits - All Planes", "", "pre_hit_cut_total_hits_leading_shower.pdf");

	histogram_functions::PlotSimpleStackInTime (h_pre_cut_collection_hits_track_nue_cc,  h_pre_cut_collection_hits_track_nue_cc_mixed,
	                                            h_pre_cut_collection_hits_track_nue_cc_out_fv,
	                                            h_pre_cut_collection_hits_track_numu_cc, h_pre_cut_collection_hits_track_numu_cc_mixed,
	                                            h_pre_cut_collection_hits_track_cosmic,  h_pre_cut_collection_hits_track_nc,
	                                            h_pre_cut_collection_hits_track_nc_pi0,  h_pre_cut_collection_hits_track_other_mixed,
	                                            h_pre_cut_collection_hits_track_unmatched, h_pre_cut_collection_hits_track_intime, intime_scale_factor,
	                                            "", "Tracks Hits - Collection Plane", "", "pre_hit_cut_collection_hits_tracks_intime.pdf");
	histogram_functions::PlotSimpleStackInTime (h_pre_cut_collection_hits_shower_nue_cc,  h_pre_cut_collection_hits_shower_nue_cc_mixed,
	                                            h_pre_cut_collection_hits_shower_nue_cc_out_fv,
	                                            h_pre_cut_collection_hits_shower_numu_cc, h_pre_cut_collection_hits_shower_numu_cc_mixed,
	                                            h_pre_cut_collection_hits_shower_cosmic,  h_pre_cut_collection_hits_shower_nc,
	                                            h_pre_cut_collection_hits_shower_nc_pi0,  h_pre_cut_collection_hits_shower_other_mixed,
	                                            h_pre_cut_collection_hits_shower_unmatched, h_pre_cut_collection_hits_shower_intime, intime_scale_factor,
	                                            "", "Showers Hits - Collection Plane", "", "pre_hit_cut_collection_hits_showers_intime.pdf");
	histogram_functions::PlotSimpleStackInTime (h_pre_cut_collection_hits_leading_shower_nue_cc,  h_pre_cut_collection_hits_leading_shower_nue_cc_mixed,
	                                            h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv,
	                                            h_pre_cut_collection_hits_leading_shower_numu_cc, h_pre_cut_collection_hits_leading_shower_numu_cc_mixed,
	                                            h_pre_cut_collection_hits_leading_shower_cosmic,  h_pre_cut_collection_hits_leading_shower_nc,
	                                            h_pre_cut_collection_hits_leading_shower_nc_pi0,  h_pre_cut_collection_hits_leading_shower_other_mixed,
	                                            h_pre_cut_collection_hits_leading_shower_unmatched, h_pre_cut_collection_hits_leading_shower_intime,
	                                            intime_scale_factor,
	                                            "", "Leading Shower Hits - Collection Plane", "", "pre_hit_cut_collection_hits_leading_shower_intime.pdf");
	histogram_functions::PlotSimpleStackInTime (h_pre_cut_total_hits_leading_shower_nue_cc,  h_pre_cut_total_hits_leading_shower_nue_cc_mixed,
	                                            h_pre_cut_total_hits_leading_shower_nue_cc_out_fv,
	                                            h_pre_cut_total_hits_leading_shower_numu_cc, h_pre_cut_total_hits_leading_shower_numu_cc_mixed,
	                                            h_pre_cut_total_hits_leading_shower_cosmic,  h_pre_cut_total_hits_leading_shower_nc,
	                                            h_pre_cut_total_hits_leading_shower_nc_pi0,  h_pre_cut_total_hits_leading_shower_other_mixed,
	                                            h_pre_cut_total_hits_leading_shower_unmatched, h_pre_cut_total_hits_leading_shower_intime,
	                                            intime_scale_factor, "",
	                                            "Leading Shower Hits - All Planes", "", "pre_hit_cut_total_hits_leading_shower_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_pre_cut_collection_hits_track_nue_cc,  h_pre_cut_collection_hits_track_nue_cc_mixed,
	                                          h_pre_cut_collection_hits_track_nue_cc_out_fv,
	                                          h_pre_cut_collection_hits_track_numu_cc, h_pre_cut_collection_hits_track_numu_cc_mixed,
	                                          h_pre_cut_collection_hits_track_cosmic,  h_pre_cut_collection_hits_track_nc,
	                                          h_pre_cut_collection_hits_track_nc_pi0,  h_pre_cut_collection_hits_track_other_mixed,
	                                          h_pre_cut_collection_hits_track_unmatched,
	                                          h_pre_cut_collection_hits_track_intime, intime_scale_factor,
	                                          h_pre_cut_collection_hits_track_data, data_scale_factor,
	                                          "", "Tracks Hits - Collection Plane", "", "pre_hit_cut_collection_hits_tracks_data.pdf");
	histogram_functions::PlotSimpleStackData (h_pre_cut_collection_hits_shower_nue_cc,  h_pre_cut_collection_hits_shower_nue_cc_mixed,
	                                          h_pre_cut_collection_hits_shower_nue_cc_out_fv,
	                                          h_pre_cut_collection_hits_shower_numu_cc, h_pre_cut_collection_hits_shower_numu_cc_mixed,
	                                          h_pre_cut_collection_hits_shower_cosmic,  h_pre_cut_collection_hits_shower_nc,
	                                          h_pre_cut_collection_hits_shower_nc_pi0,  h_pre_cut_collection_hits_shower_other_mixed,
	                                          h_pre_cut_collection_hits_shower_unmatched,
	                                          h_pre_cut_collection_hits_shower_intime, intime_scale_factor,
	                                          h_pre_cut_collection_hits_shower_data, data_scale_factor,
	                                          "", "Showers Hits - Collection Plane", "", "pre_hit_cut_collection_hits_showers_data.pdf");
	histogram_functions::PlotSimpleStackData (h_pre_cut_collection_hits_leading_shower_nue_cc,  h_pre_cut_collection_hits_leading_shower_nue_cc_mixed,
	                                          h_pre_cut_collection_hits_leading_shower_nue_cc_out_fv,
	                                          h_pre_cut_collection_hits_leading_shower_numu_cc, h_pre_cut_collection_hits_leading_shower_numu_cc_mixed,
	                                          h_pre_cut_collection_hits_leading_shower_cosmic,  h_pre_cut_collection_hits_leading_shower_nc,
	                                          h_pre_cut_collection_hits_leading_shower_nc_pi0,  h_pre_cut_collection_hits_leading_shower_other_mixed,
	                                          h_pre_cut_collection_hits_leading_shower_unmatched,
	                                          h_pre_cut_collection_hits_leading_shower_intime,intime_scale_factor,
	                                          h_pre_cut_collection_hits_leading_shower_data, data_scale_factor,
	                                          "", "Leading Shower Hits - Collection Plane", "", "pre_hit_cut_collection_hits_leading_shower_data.pdf");
	histogram_functions::PlotSimpleStackData (h_pre_cut_total_hits_leading_shower_nue_cc,  h_pre_cut_total_hits_leading_shower_nue_cc_mixed,
	                                          h_pre_cut_total_hits_leading_shower_nue_cc_out_fv,
	                                          h_pre_cut_total_hits_leading_shower_numu_cc, h_pre_cut_total_hits_leading_shower_numu_cc_mixed,
	                                          h_pre_cut_total_hits_leading_shower_cosmic,  h_pre_cut_total_hits_leading_shower_nc,
	                                          h_pre_cut_total_hits_leading_shower_nc_pi0,  h_pre_cut_total_hits_leading_shower_other_mixed,
	                                          h_pre_cut_total_hits_leading_shower_unmatched,
	                                          h_pre_cut_total_hits_leading_shower_intime,intime_scale_factor,
	                                          h_pre_cut_total_hits_leading_shower_data, data_scale_factor,
	                                          "", "Leading Shower Hits - All Planes", "", "pre_hit_cut_total_hits_leading_shower_data.pdf");

	histogram_functions::Plot2DHistogram (h_ele_eng_total_hits, "", "Hits - All Planes", "True Electron Energy [GeV]", "post_cuts_ele_eng_total_hits.pdf");
	histogram_functions::Plot2DHistogram (h_ele_eng_colleciton_hits, "", "Hits - Collection Plane",
	                                      "True Electron Energy [GeV]", "post_cuts_ele_eng_collection_hits.pdf");
	histogram_functions::Plot2DHistogram (h_nu_eng_total_hits, "", "Hits - All Planes", "True Neutrino Energy [GeV]", "post_cuts_nu_eng_total_hits.pdf");
	histogram_functions::Plot2DHistogram (h_nu_eng_collection_hits, "", "Hits - Collection Planes",
	                                      "True Neutrino Energy [GeV]", "post_cuts_nu_eng_colleciton_hits.pdf");
	histogram_functions::PlotSimpleStack (h_ele_cos_theta_last_nue_cc,  h_ele_cos_theta_last_nue_cc_mixed,
	                                      h_ele_cos_theta_last_nue_cc_out_fv,
	                                      h_ele_cos_theta_last_numu_cc, h_ele_cos_theta_last_numu_cc_mixed,
	                                      h_ele_cos_theta_last_cosmic,  h_ele_cos_theta_last_nc,
	                                      h_ele_cos_theta_last_nc_pi0,  h_ele_cos_theta_last_other_mixed,
	                                      h_ele_cos_theta_last_unmatched, 0.15, 0.35, 0.70, 0.95, "",
	                                      "Leading Shower Cos(#theta)", "", "post_cuts_leading_cos_theta_last.pdf");
	histogram_functions::PlotSimpleStack (h_ele_cos_theta_nue_cc,  h_ele_cos_theta_nue_cc_mixed,
	                                      h_ele_cos_theta_nue_cc_out_fv,
	                                      h_ele_cos_theta_numu_cc, h_ele_cos_theta_numu_cc_mixed,
	                                      h_ele_cos_theta_cosmic,  h_ele_cos_theta_nc,
	                                      h_ele_cos_theta_nc_pi0,  h_ele_cos_theta_other_mixed,
	                                      h_ele_cos_theta_unmatched, 0.15, 0.35, 0.70, 0.95, "",
	                                      "Leading Shower Cos(#theta)", "", "post_cuts_leading_cos_theta.pdf");

	histogram_functions::PlotSimpleStackInTime (h_ele_cos_theta_last_nue_cc,  h_ele_cos_theta_last_nue_cc_mixed,
	                                            h_ele_cos_theta_last_nue_cc_out_fv,
	                                            h_ele_cos_theta_last_numu_cc, h_ele_cos_theta_last_numu_cc_mixed,
	                                            h_ele_cos_theta_last_cosmic,  h_ele_cos_theta_last_nc,
	                                            h_ele_cos_theta_last_nc_pi0,  h_ele_cos_theta_last_other_mixed,
	                                            h_ele_cos_theta_last_unmatched, h_ele_cos_theta_last_intime, intime_scale_factor,
	                                            0.15, 0.35, 0.70, 0.95, "",
	                                            "Leading Shower Cos(#theta)", "", "post_cuts_leading_cos_theta_last_intime.pdf");
	histogram_functions::PlotSimpleStackInTime (h_ele_cos_theta_nue_cc,  h_ele_cos_theta_nue_cc_mixed,
	                                            h_ele_cos_theta_last_nue_cc_out_fv,
	                                            h_ele_cos_theta_numu_cc, h_ele_cos_theta_numu_cc_mixed,
	                                            h_ele_cos_theta_cosmic,  h_ele_cos_theta_nc,
	                                            h_ele_cos_theta_nc_pi0,  h_ele_cos_theta_other_mixed,
	                                            h_ele_cos_theta_unmatched, h_ele_cos_theta_intime, intime_scale_factor,
	                                            0.15, 0.35, 0.70, 0.95, "",
	                                            "Leading Shower Cos(#theta)", "", "post_cuts_leading_cos_theta_intime.pdf");

	histogram_functions::PlotSimpleStackData (h_ele_cos_theta_last_nue_cc,  h_ele_cos_theta_last_nue_cc_mixed,
	                                          h_ele_cos_theta_last_nue_cc_out_fv,
	                                          h_ele_cos_theta_last_numu_cc, h_ele_cos_theta_last_numu_cc_mixed,
	                                          h_ele_cos_theta_last_cosmic,  h_ele_cos_theta_last_nc,
	                                          h_ele_cos_theta_last_nc_pi0,  h_ele_cos_theta_last_other_mixed,
	                                          h_ele_cos_theta_last_unmatched, h_ele_cos_theta_last_intime, intime_scale_factor,
	                                          h_ele_cos_theta_last_data, data_scale_factor,
	                                          0.15, 0.35, 0.70, 0.95, "",
	                                          "Leading Shower Cos(#theta)", "", "post_cuts_leading_cos_theta_last_data.pdf");

	histogram_functions::PlotSimpleStackData (h_ele_cos_theta_last_trans_nue_cc,    h_ele_cos_theta_last_trans_nue_cc_mixed,
	                                          h_ele_cos_theta_last_trans_nue_cc_out_fv,
	                                          h_ele_cos_theta_last_trans_numu_cc,   h_ele_cos_theta_last_trans_numu_cc_mixed,
	                                          h_ele_cos_theta_last_trans_cosmic,    h_ele_cos_theta_last_trans_nc,
	                                          h_ele_cos_theta_last_trans_nc_pi0,    h_ele_cos_theta_last_trans_other_mixed,
	                                          h_ele_cos_theta_last_trans_unmatched, h_ele_cos_theta_last_trans_intime, intime_scale_factor,
	                                          h_ele_cos_theta_last_trans_data, data_scale_factor,
	                                          0.15, 0.35, 0.70, 0.95, "",
	                                          "Leading Shower Cos(#theta)", "", "post_cuts_leading_cos_theta_last_trans_data.pdf");

	histogram_functions::PlotSimpleStackData (h_ele_cos_theta_nue_cc,  h_ele_cos_theta_nue_cc_mixed,
	                                          h_ele_cos_theta_nue_cc_out_fv,
	                                          h_ele_cos_theta_numu_cc, h_ele_cos_theta_numu_cc_mixed,
	                                          h_ele_cos_theta_cosmic,  h_ele_cos_theta_nc,
	                                          h_ele_cos_theta_nc_pi0,  h_ele_cos_theta_other_mixed,
	                                          h_ele_cos_theta_unmatched, h_ele_cos_theta_intime, intime_scale_factor,
	                                          h_ele_cos_theta_data, data_scale_factor,
	                                          0.15, 0.35, 0.70, 0.95, "",
	                                          "Leading Shower Cos(#theta)", "", "post_cuts_leading_cos_theta_data.pdf");

	histogram_functions::PlotSimpleStack (h_ele_pfp_momentum_nue_cc,  h_ele_pfp_momentum_nue_cc_mixed,
	                                      h_ele_pfp_momentum_nue_cc_out_fv,
	                                      h_ele_pfp_momentum_numu_cc, h_ele_pfp_momentum_numu_cc_mixed,
	                                      h_ele_pfp_momentum_cosmic,  h_ele_pfp_momentum_nc,
	                                      h_ele_pfp_momentum_nc_pi0,  h_ele_pfp_momentum_other_mixed,
	                                      h_ele_pfp_momentum_unmatched, "",
	                                      "Leading Shower Momentum [GeV]", "", "post_cuts_leading_momentum.pdf");
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_momentum_nue_cc,  h_ele_pfp_momentum_nue_cc_mixed,
	                                            h_ele_pfp_momentum_nue_cc_out_fv,
	                                            h_ele_pfp_momentum_numu_cc, h_ele_pfp_momentum_numu_cc_mixed,
	                                            h_ele_pfp_momentum_cosmic,  h_ele_pfp_momentum_nc,
	                                            h_ele_pfp_momentum_nc_pi0,  h_ele_pfp_momentum_other_mixed,
	                                            h_ele_pfp_momentum_unmatched, h_ele_pfp_momentum_intime, intime_scale_factor, "",
	                                            "Leading Shower Momentum [GeV]", "", "post_cuts_leading_momentum_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_ele_pfp_momentum_nue_cc,  h_ele_pfp_momentum_nue_cc_mixed,
	                                          h_ele_pfp_momentum_nue_cc_out_fv,
	                                          h_ele_pfp_momentum_numu_cc, h_ele_pfp_momentum_numu_cc_mixed,
	                                          h_ele_pfp_momentum_cosmic,  h_ele_pfp_momentum_nc,
	                                          h_ele_pfp_momentum_nc_pi0,  h_ele_pfp_momentum_other_mixed,
	                                          h_ele_pfp_momentum_unmatched, h_ele_pfp_momentum_intime, intime_scale_factor,
	                                          h_ele_pfp_momentum_data, data_scale_factor,
	                                          "", "Leading Shower Momentum [GeV]", "", "post_cuts_leading_momentum_data.pdf");

	histogram_functions::PlotSimpleStack (h_ele_pfp_theta_nue_cc,  h_ele_pfp_theta_nue_cc_mixed,
	                                      h_ele_pfp_theta_nue_cc_out_fv,
	                                      h_ele_pfp_theta_numu_cc, h_ele_pfp_theta_numu_cc_mixed,
	                                      h_ele_pfp_theta_cosmic,  h_ele_pfp_theta_nc,
	                                      h_ele_pfp_theta_nc_pi0,  h_ele_pfp_theta_other_mixed,
	                                      h_ele_pfp_theta_unmatched, "",
	                                      "Leading Shower Theta [Degrees]", "", "pre_collection_cut_leading_theta.pdf");
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_theta_nue_cc,  h_ele_pfp_theta_nue_cc_mixed,
	                                            h_ele_pfp_theta_nue_cc_out_fv,
	                                            h_ele_pfp_theta_numu_cc, h_ele_pfp_theta_numu_cc_mixed,
	                                            h_ele_pfp_theta_cosmic,  h_ele_pfp_theta_nc,
	                                            h_ele_pfp_theta_nc_pi0,  h_ele_pfp_theta_other_mixed,
	                                            h_ele_pfp_theta_unmatched, h_ele_pfp_theta_intime, intime_scale_factor, "",
	                                            "Leading Shower Theta [Degrees]", "", "pre_collection_cut_leading_theta_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_ele_pfp_theta_nue_cc,  h_ele_pfp_theta_nue_cc_mixed,
	                                          h_ele_pfp_theta_nue_cc_out_fv,
	                                          h_ele_pfp_theta_numu_cc, h_ele_pfp_theta_numu_cc_mixed,
	                                          h_ele_pfp_theta_cosmic,  h_ele_pfp_theta_nc,
	                                          h_ele_pfp_theta_nc_pi0,  h_ele_pfp_theta_other_mixed,
	                                          h_ele_pfp_theta_unmatched, h_ele_pfp_theta_intime, intime_scale_factor,
	                                          h_ele_pfp_theta_data, data_scale_factor,
	                                          "", "Leading Shower Theta [Degrees]", "", "pre_collection_cut_leading_theta_data.pdf");

	histogram_functions::PlotSimpleStack (h_ele_pfp_theta_after_nue_cc,  h_ele_pfp_theta_after_nue_cc_mixed,
	                                      h_ele_pfp_theta_after_nue_cc_out_fv,
	                                      h_ele_pfp_theta_after_numu_cc, h_ele_pfp_theta_after_numu_cc_mixed,
	                                      h_ele_pfp_theta_after_cosmic,  h_ele_pfp_theta_after_nc,
	                                      h_ele_pfp_theta_after_nc_pi0,  h_ele_pfp_theta_after_other_mixed,
	                                      h_ele_pfp_theta_after_unmatched, "",
	                                      "Leading Shower Theta [Degrees]", "", "post_collection_cut_leading_theta.pdf");
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_theta_after_nue_cc,  h_ele_pfp_theta_after_nue_cc_mixed,
	                                            h_ele_pfp_theta_after_nue_cc_out_fv,
	                                            h_ele_pfp_theta_after_numu_cc, h_ele_pfp_theta_after_numu_cc_mixed,
	                                            h_ele_pfp_theta_after_cosmic,  h_ele_pfp_theta_after_nc,
	                                            h_ele_pfp_theta_after_nc_pi0,  h_ele_pfp_theta_after_other_mixed,
	                                            h_ele_pfp_theta_after_unmatched, h_ele_pfp_theta_after_intime, intime_scale_factor, "",
	                                            "Leading Shower Theta [Degrees]", "", "post_collection_cut_leading_theta_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_ele_pfp_theta_after_nue_cc,  h_ele_pfp_theta_after_nue_cc_mixed,
	                                          h_ele_pfp_theta_after_nue_cc_out_fv,
	                                          h_ele_pfp_theta_after_numu_cc, h_ele_pfp_theta_after_numu_cc_mixed,
	                                          h_ele_pfp_theta_after_cosmic,  h_ele_pfp_theta_after_nc,
	                                          h_ele_pfp_theta_after_nc_pi0,  h_ele_pfp_theta_after_other_mixed,
	                                          h_ele_pfp_theta_after_unmatched, h_ele_pfp_theta_after_intime, intime_scale_factor,
	                                          h_ele_pfp_theta_after_data, data_scale_factor,
	                                          "", "Leading Shower Theta [Degrees]", "", "post_collection_cut_leading_theta_data.pdf");

	histogram_functions::PlotSimpleStack (h_ele_pfp_theta_last_nue_cc,  h_ele_pfp_theta_last_nue_cc_mixed,
	                                      h_ele_pfp_theta_last_nue_cc_out_fv,
	                                      h_ele_pfp_theta_last_numu_cc, h_ele_pfp_theta_last_numu_cc_mixed,
	                                      h_ele_pfp_theta_last_cosmic,  h_ele_pfp_theta_last_nc,
	                                      h_ele_pfp_theta_last_nc_pi0,  h_ele_pfp_theta_last_other_mixed,
	                                      h_ele_pfp_theta_last_unmatched, "",
	                                      "Leading Shower Theta [Degrees]", "", "post_cuts_leading_theta.pdf");
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_theta_last_nue_cc,  h_ele_pfp_theta_last_nue_cc_mixed,
	                                            h_ele_pfp_theta_last_nue_cc_out_fv,
	                                            h_ele_pfp_theta_last_numu_cc, h_ele_pfp_theta_last_numu_cc_mixed,
	                                            h_ele_pfp_theta_last_cosmic,  h_ele_pfp_theta_last_nc,
	                                            h_ele_pfp_theta_last_nc_pi0,  h_ele_pfp_theta_last_other_mixed,
	                                            h_ele_pfp_theta_last_unmatched, h_ele_pfp_theta_last_intime, intime_scale_factor, "",
	                                            "Leading Shower Theta [Degrees]", "", "post_cuts_leading_theta_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_ele_pfp_theta_last_nue_cc,  h_ele_pfp_theta_last_nue_cc_mixed,
	                                          h_ele_pfp_theta_last_nue_cc_out_fv,
	                                          h_ele_pfp_theta_last_numu_cc, h_ele_pfp_theta_last_numu_cc_mixed,
	                                          h_ele_pfp_theta_last_cosmic,  h_ele_pfp_theta_last_nc,
	                                          h_ele_pfp_theta_last_nc_pi0,  h_ele_pfp_theta_last_other_mixed,
	                                          h_ele_pfp_theta_last_unmatched, h_ele_pfp_theta_last_intime, intime_scale_factor,
	                                          h_ele_pfp_theta_last_data, data_scale_factor,
	                                          "", "Leading Shower Theta [Degrees]", "", "post_cuts_leading_theta_last_data.pdf");

	histogram_functions::PlotSimpleStack (h_ele_pfp_phi_nue_cc,  h_ele_pfp_phi_nue_cc_mixed,
	                                      h_ele_pfp_phi_nue_cc_out_fv,
	                                      h_ele_pfp_phi_numu_cc, h_ele_pfp_phi_numu_cc_mixed,
	                                      h_ele_pfp_phi_cosmic,  h_ele_pfp_phi_nc,
	                                      h_ele_pfp_phi_nc_pi0,  h_ele_pfp_phi_other_mixed,
	                                      h_ele_pfp_phi_unmatched, "",
	                                      "Leading Shower Phi [Degrees]", "", "pre_collection_cut_leading_phi.pdf");
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_phi_nue_cc,  h_ele_pfp_phi_nue_cc_mixed,
	                                            h_ele_pfp_phi_nue_cc_out_fv,
	                                            h_ele_pfp_phi_numu_cc, h_ele_pfp_phi_numu_cc_mixed,
	                                            h_ele_pfp_phi_cosmic,  h_ele_pfp_phi_nc,
	                                            h_ele_pfp_phi_nc_pi0,  h_ele_pfp_phi_other_mixed,
	                                            h_ele_pfp_phi_unmatched, h_ele_pfp_phi_intime, intime_scale_factor, "",
	                                            "Leading Shower Phi [Degrees]", "", "pre_collection_cut_leading_phi_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_ele_pfp_phi_nue_cc,  h_ele_pfp_phi_nue_cc_mixed,
	                                          h_ele_pfp_phi_nue_cc_out_fv,
	                                          h_ele_pfp_phi_numu_cc, h_ele_pfp_phi_numu_cc_mixed,
	                                          h_ele_pfp_phi_cosmic,  h_ele_pfp_phi_nc,
	                                          h_ele_pfp_phi_nc_pi0,  h_ele_pfp_phi_other_mixed,
	                                          h_ele_pfp_phi_unmatched, h_ele_pfp_phi_intime, intime_scale_factor,
	                                          h_ele_pfp_phi_data, data_scale_factor,
	                                          "", "Leading Shower Phi [Degrees]", "", "pre_collection_cut_leading_phi_data.pdf");

	histogram_functions::PlotSimpleStack (h_ele_pfp_phi_after_nue_cc,  h_ele_pfp_phi_after_nue_cc_mixed,
	                                      h_ele_pfp_phi_after_nue_cc_out_fv,
	                                      h_ele_pfp_phi_after_numu_cc, h_ele_pfp_phi_after_numu_cc_mixed,
	                                      h_ele_pfp_phi_after_cosmic,  h_ele_pfp_phi_after_nc,
	                                      h_ele_pfp_phi_after_nc_pi0,  h_ele_pfp_phi_after_other_mixed,
	                                      h_ele_pfp_phi_after_unmatched, "",
	                                      "Leading Shower Phi [Degrees]", "", "post_collection_cut_leading_phi.pdf");
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_phi_after_nue_cc,  h_ele_pfp_phi_after_nue_cc_mixed,
	                                            h_ele_pfp_phi_after_nue_cc_out_fv,
	                                            h_ele_pfp_phi_after_numu_cc, h_ele_pfp_phi_after_numu_cc_mixed,
	                                            h_ele_pfp_phi_after_cosmic,  h_ele_pfp_phi_after_nc,
	                                            h_ele_pfp_phi_after_nc_pi0,  h_ele_pfp_phi_after_other_mixed,
	                                            h_ele_pfp_phi_after_unmatched, h_ele_pfp_phi_after_intime, intime_scale_factor, "",
	                                            "Leading Shower Phi [Degrees]", "", "post_collection_cut_leading_phi_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_ele_pfp_phi_after_nue_cc,  h_ele_pfp_phi_after_nue_cc_mixed,
	                                          h_ele_pfp_phi_after_nue_cc_out_fv,
	                                          h_ele_pfp_phi_after_numu_cc, h_ele_pfp_phi_after_numu_cc_mixed,
	                                          h_ele_pfp_phi_after_cosmic,  h_ele_pfp_phi_after_nc,
	                                          h_ele_pfp_phi_after_nc_pi0,  h_ele_pfp_phi_after_other_mixed,
	                                          h_ele_pfp_phi_after_unmatched, h_ele_pfp_phi_after_intime, intime_scale_factor,
	                                          h_ele_pfp_phi_after_data, data_scale_factor,
	                                          "", "Leading Shower Phi [Degrees]", "", "post_collection_cut_leading_phi_data.pdf");

	histogram_functions::PlotSimpleStack (h_ele_pfp_phi_last_nue_cc,  h_ele_pfp_phi_last_nue_cc_mixed,
	                                      h_ele_pfp_phi_last_nue_cc_out_fv,
	                                      h_ele_pfp_phi_last_numu_cc, h_ele_pfp_phi_last_numu_cc_mixed,
	                                      h_ele_pfp_phi_last_cosmic,  h_ele_pfp_phi_last_nc,
	                                      h_ele_pfp_phi_last_nc_pi0,  h_ele_pfp_phi_last_other_mixed,
	                                      h_ele_pfp_phi_last_unmatched, "",
	                                      "Leading Shower Phi [Degrees]", "", "post_cuts_leading_phi.pdf");
	histogram_functions::PlotSimpleStackInTime (h_ele_pfp_phi_last_nue_cc,  h_ele_pfp_phi_last_nue_cc_mixed,
	                                            h_ele_pfp_phi_last_nue_cc_out_fv,
	                                            h_ele_pfp_phi_last_numu_cc, h_ele_pfp_phi_last_numu_cc_mixed,
	                                            h_ele_pfp_phi_last_cosmic,  h_ele_pfp_phi_last_nc,
	                                            h_ele_pfp_phi_last_nc_pi0,  h_ele_pfp_phi_last_other_mixed,
	                                            h_ele_pfp_phi_last_unmatched, h_ele_pfp_phi_last_intime, intime_scale_factor, "",
	                                            "Leading Shower Phi [Degrees]", "", "post_cuts_leading_phi_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_ele_pfp_phi_last_nue_cc,  h_ele_pfp_phi_last_nue_cc_mixed,
	                                          h_ele_pfp_phi_last_nue_cc_out_fv,
	                                          h_ele_pfp_phi_last_numu_cc, h_ele_pfp_phi_last_numu_cc_mixed,
	                                          h_ele_pfp_phi_last_cosmic,  h_ele_pfp_phi_last_nc,
	                                          h_ele_pfp_phi_last_nc_pi0,  h_ele_pfp_phi_last_other_mixed,
	                                          h_ele_pfp_phi_last_unmatched, h_ele_pfp_phi_last_intime, intime_scale_factor,
	                                          h_ele_pfp_phi_last_data, data_scale_factor,
	                                          "", "Leading Shower Phi [Degrees]", "", "post_cuts_leading_phi_last_data.pdf");

	histogram_functions::PlotSimpleStack (h_leading_shwr_length_1shwr_nue_cc,  h_leading_shwr_length_1shwr_nue_cc_mixed,
	                                      h_leading_shwr_length_1shwr_nue_cc_out_fv,
	                                      h_leading_shwr_length_1shwr_numu_cc, h_leading_shwr_length_1shwr_numu_cc_mixed,
	                                      h_leading_shwr_length_1shwr_cosmic,  h_leading_shwr_length_1shwr_nc,
	                                      h_leading_shwr_length_1shwr_nc_pi0,  h_leading_shwr_length_1shwr_other_mixed,
	                                      h_leading_shwr_length_1shwr_unmatched, "",
	                                      "Leading Shower Length (1 Shower Events) [cm]", "", "post_cuts_leading_length_1shwr.pdf");
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_length_1shwr_nue_cc,  h_leading_shwr_length_1shwr_nue_cc_mixed,
	                                            h_leading_shwr_length_1shwr_nue_cc_out_fv,
	                                            h_leading_shwr_length_1shwr_numu_cc, h_leading_shwr_length_1shwr_numu_cc_mixed,
	                                            h_leading_shwr_length_1shwr_cosmic,  h_leading_shwr_length_1shwr_nc,
	                                            h_leading_shwr_length_1shwr_nc_pi0,  h_leading_shwr_length_1shwr_other_mixed,
	                                            h_leading_shwr_length_1shwr_unmatched, h_leading_shwr_length_1shwr_intime, intime_scale_factor, "",
	                                            "Leading Shower Length (1 Shower Events) [cm]", "", "post_cuts_leading_length_1shwr_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_leading_shwr_length_1shwr_nue_cc,  h_leading_shwr_length_1shwr_nue_cc_mixed,
	                                          h_leading_shwr_length_1shwr_nue_cc_out_fv,
	                                          h_leading_shwr_length_1shwr_numu_cc, h_leading_shwr_length_1shwr_numu_cc_mixed,
	                                          h_leading_shwr_length_1shwr_cosmic,  h_leading_shwr_length_1shwr_nc,
	                                          h_leading_shwr_length_1shwr_nc_pi0,  h_leading_shwr_length_1shwr_other_mixed,
	                                          h_leading_shwr_length_1shwr_unmatched, h_leading_shwr_length_1shwr_intime, intime_scale_factor,
	                                          h_leading_shwr_length_1shwr_data, data_scale_factor,
	                                          "", "Leading Shower Length (1 Shower Events) [cm]", "", "post_cuts_leading_length_1shwr_data.pdf");

	histogram_functions::PlotSimpleStack (h_leading_shwr_length_2shwr_nue_cc,  h_leading_shwr_length_2shwr_nue_cc_mixed,
	                                      h_leading_shwr_length_2shwr_nue_cc_out_fv,
	                                      h_leading_shwr_length_2shwr_numu_cc, h_leading_shwr_length_2shwr_numu_cc_mixed,
	                                      h_leading_shwr_length_2shwr_cosmic,  h_leading_shwr_length_2shwr_nc,
	                                      h_leading_shwr_length_2shwr_nc_pi0,  h_leading_shwr_length_2shwr_other_mixed,
	                                      h_leading_shwr_length_2shwr_unmatched, "",
	                                      "Leading Shower Length (2+ Shower Events) [cm]", "", "post_cuts_leading_length_2shwr.pdf");
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_length_2shwr_nue_cc,  h_leading_shwr_length_2shwr_nue_cc_mixed,
	                                            h_leading_shwr_length_2shwr_nue_cc_out_fv,
	                                            h_leading_shwr_length_2shwr_numu_cc, h_leading_shwr_length_2shwr_numu_cc_mixed,
	                                            h_leading_shwr_length_2shwr_cosmic,  h_leading_shwr_length_2shwr_nc,
	                                            h_leading_shwr_length_2shwr_nc_pi0,  h_leading_shwr_length_2shwr_other_mixed,
	                                            h_leading_shwr_length_2shwr_unmatched, h_leading_shwr_length_2shwr_intime, intime_scale_factor, "",
	                                            "Leading Shower Length (2+ Shower Events) [cm]", "", "post_cuts_leading_length_2shwr_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_leading_shwr_length_2shwr_nue_cc,  h_leading_shwr_length_2shwr_nue_cc_mixed,
	                                          h_leading_shwr_length_2shwr_nue_cc_out_fv,
	                                          h_leading_shwr_length_2shwr_numu_cc, h_leading_shwr_length_2shwr_numu_cc_mixed,
	                                          h_leading_shwr_length_2shwr_cosmic,  h_leading_shwr_length_2shwr_nc,
	                                          h_leading_shwr_length_2shwr_nc_pi0,  h_leading_shwr_length_2shwr_other_mixed,
	                                          h_leading_shwr_length_2shwr_unmatched, h_leading_shwr_length_2shwr_intime, intime_scale_factor,
	                                          h_leading_shwr_length_2shwr_data, data_scale_factor,
	                                          "", "Leading Shower Length (2+ Shower Events) [cm]", "", "post_cuts_leading_length_2shwr_data.pdf");

	histogram_functions::PlotSimpleStack (h_leading_shwr_hits_1shwr_nue_cc,  h_leading_shwr_hits_1shwr_nue_cc_mixed,
	                                      h_leading_shwr_hits_1shwr_nue_cc_out_fv,
	                                      h_leading_shwr_hits_1shwr_numu_cc, h_leading_shwr_hits_1shwr_numu_cc_mixed,
	                                      h_leading_shwr_hits_1shwr_cosmic,  h_leading_shwr_hits_1shwr_nc,
	                                      h_leading_shwr_hits_1shwr_nc_pi0,  h_leading_shwr_hits_1shwr_other_mixed,
	                                      h_leading_shwr_hits_1shwr_unmatched, "",
	                                      "Leading Shower Hits (1 Shower Events)", "", "post_cuts_leading_hits_1shwr.pdf");
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_hits_1shwr_nue_cc,  h_leading_shwr_hits_1shwr_nue_cc_mixed,
	                                            h_leading_shwr_hits_1shwr_nue_cc_out_fv,
	                                            h_leading_shwr_hits_1shwr_numu_cc, h_leading_shwr_hits_1shwr_numu_cc_mixed,
	                                            h_leading_shwr_hits_1shwr_cosmic,  h_leading_shwr_hits_1shwr_nc,
	                                            h_leading_shwr_hits_1shwr_nc_pi0,  h_leading_shwr_hits_1shwr_other_mixed,
	                                            h_leading_shwr_hits_1shwr_unmatched, h_leading_shwr_hits_1shwr_intime, intime_scale_factor, "",
	                                            "Leading Shower Hits (1 Shower Events)", "", "post_cuts_leading_hits_1shwr_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_leading_shwr_hits_1shwr_nue_cc,  h_leading_shwr_hits_1shwr_nue_cc_mixed,
	                                          h_leading_shwr_hits_1shwr_nue_cc_out_fv,
	                                          h_leading_shwr_hits_1shwr_numu_cc, h_leading_shwr_hits_1shwr_numu_cc_mixed,
	                                          h_leading_shwr_hits_1shwr_cosmic,  h_leading_shwr_hits_1shwr_nc,
	                                          h_leading_shwr_hits_1shwr_nc_pi0,  h_leading_shwr_hits_1shwr_other_mixed,
	                                          h_leading_shwr_hits_1shwr_unmatched, h_leading_shwr_hits_1shwr_intime, intime_scale_factor,
	                                          h_leading_shwr_hits_1shwr_data, data_scale_factor,
	                                          "", "Leading Shower Hits (1 Shower Events)", "", "post_cuts_leading_hits_1shwr_data.pdf");

	histogram_functions::PlotSimpleStack (h_leading_shwr_hits_2shwr_nue_cc,  h_leading_shwr_hits_2shwr_nue_cc_mixed,
	                                      h_leading_shwr_hits_2shwr_nue_cc_out_fv,
	                                      h_leading_shwr_hits_2shwr_numu_cc, h_leading_shwr_hits_2shwr_numu_cc_mixed,
	                                      h_leading_shwr_hits_2shwr_cosmic,  h_leading_shwr_hits_2shwr_nc,
	                                      h_leading_shwr_hits_2shwr_nc_pi0,  h_leading_shwr_hits_2shwr_other_mixed,
	                                      h_leading_shwr_hits_2shwr_unmatched, "",
	                                      "Leading Shower Hits (2+ Shower Events)", "", "post_cuts_leading_hits_2shwr.pdf");
	histogram_functions::PlotSimpleStackInTime (h_leading_shwr_hits_2shwr_nue_cc,  h_leading_shwr_hits_2shwr_nue_cc_mixed,
	                                            h_leading_shwr_hits_2shwr_nue_cc_out_fv,
	                                            h_leading_shwr_hits_2shwr_numu_cc, h_leading_shwr_hits_2shwr_numu_cc_mixed,
	                                            h_leading_shwr_hits_2shwr_cosmic,  h_leading_shwr_hits_2shwr_nc,
	                                            h_leading_shwr_hits_2shwr_nc_pi0,  h_leading_shwr_hits_2shwr_other_mixed,
	                                            h_leading_shwr_hits_2shwr_unmatched, h_leading_shwr_hits_2shwr_intime, intime_scale_factor, "",
	                                            "Leading Shower Hits (2+ Shower Events)", "", "post_cuts_leading_hits_2shwr_intime.pdf");
	histogram_functions::PlotSimpleStackData (h_leading_shwr_hits_2shwr_nue_cc,  h_leading_shwr_hits_2shwr_nue_cc_mixed,
	                                          h_leading_shwr_hits_2shwr_nue_cc_out_fv,
	                                          h_leading_shwr_hits_2shwr_numu_cc, h_leading_shwr_hits_2shwr_numu_cc_mixed,
	                                          h_leading_shwr_hits_2shwr_cosmic,  h_leading_shwr_hits_2shwr_nc,
	                                          h_leading_shwr_hits_2shwr_nc_pi0,  h_leading_shwr_hits_2shwr_other_mixed,
	                                          h_leading_shwr_hits_2shwr_unmatched, h_leading_shwr_hits_2shwr_intime, intime_scale_factor,
	                                          h_leading_shwr_hits_2shwr_data, data_scale_factor,
	                                          "", "Leading Shower Hits (2+ Shower Events)", "", "post_cuts_leading_hits_2shwr_data.pdf");

	histogram_functions::PlotSimpleStackData(h_ele_pfp_x_nue_cc, h_ele_pfp_x_nue_cc_mixed,
	                                         h_ele_pfp_x_nue_cc_out_fv, h_ele_pfp_x_numu_cc, h_ele_pfp_x_numu_cc_mixed,
	                                         h_ele_pfp_x_cosmic, h_ele_pfp_x_nc, h_ele_pfp_x_nc_pi0, h_ele_pfp_x_other_mixed,
	                                         h_ele_pfp_x_unmatched, h_ele_pfp_x_intime, intime_scale_factor,
	                                         h_ele_pfp_x_data, data_scale_factor,
	                                         "", "Reco Vertex X [cm]", "", "post_cuts_leading_pfp_x_data.pdf");
	histogram_functions::PlotSimpleStackData(h_ele_pfp_y_nue_cc, h_ele_pfp_y_nue_cc_mixed,
	                                         h_ele_pfp_y_nue_cc_out_fv, h_ele_pfp_y_numu_cc, h_ele_pfp_y_numu_cc_mixed,
	                                         h_ele_pfp_y_cosmic, h_ele_pfp_y_nc, h_ele_pfp_y_nc_pi0, h_ele_pfp_y_other_mixed,
	                                         h_ele_pfp_y_unmatched, h_ele_pfp_y_intime, intime_scale_factor,
	                                         h_ele_pfp_y_data, data_scale_factor,
	                                         "", "Reco Vertex Y [cm]", "", "post_cuts_leading_pfp_y_data.pdf");
	histogram_functions::PlotSimpleStackData(h_ele_pfp_z_nue_cc, h_ele_pfp_z_nue_cc_mixed,
	                                         h_ele_pfp_z_nue_cc_out_fv, h_ele_pfp_z_numu_cc, h_ele_pfp_z_numu_cc_mixed,
	                                         h_ele_pfp_z_cosmic, h_ele_pfp_z_nc, h_ele_pfp_z_nc_pi0, h_ele_pfp_z_other_mixed,
	                                         h_ele_pfp_z_unmatched, h_ele_pfp_z_intime, intime_scale_factor,
	                                         h_ele_pfp_z_data, data_scale_factor,
	                                         "", "Reco Vertex Z [cm]", "", "post_cuts_leading_pfp_z_data.pdf");


	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_qe, "Post Cuts - Showers/Tracks Purity - QE",
	                                      "Reco Showers", "Reco Tracks", "post_cuts_showers_tracks_purity_qe.pdf", "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_res, "Post Cuts - Showers/Tracks Purity - Res",
	                                      "Reco Showers", "Reco Tracks", "post_cuts_showers_tracks_purity_res.pdf", "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_dis, "Post Cuts - Showers/Tracks Purity - DIS",
	                                      "Reco Showers", "Reco Tracks", "post_cuts_showers_tracks_purity_dis.pdf", "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_coh, "Post Cuts - Showers/Tracks Purity - Coh",
	                                      "Reco Showers", "Reco Tracks", "post_cuts_showers_tracks_purity_coh.pdf", "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_mec, "Post Cuts - Showers/Tracks Purity - MEC",
	                                      "Reco Showers", "Reco Tracks", "post_cuts_showers_tracks_purity_mec.pdf", "colz text", 3, 2);
	histogram_functions::Plot2DHistogram (h_post_cuts_num_tracks_showers_purity_total, "Post Cuts - Showers/Tracks Purity - Total",
	                                      "Reco Showers", "Reco Tracks", "post_cuts_showers_tracks_purity_total.pdf", "colz text", 3, 2);

	histogram_functions::OverlayScatter(h_ele_theta_phi_nue_cc, h_ele_theta_phi_nue_cc_mixed, h_ele_theta_phi_nue_cc_out_fv, h_ele_theta_phi_numu_cc,
	                                    h_ele_theta_phi_numu_cc_mixed, h_ele_theta_phi_cosmic, h_ele_theta_phi_nc,
	                                    h_ele_theta_phi_nc_pi0, h_ele_theta_phi_other_mixed, h_ele_theta_phi_unmatched,
	                                    0.15, 0.35, 0.65, 0.90, "", "Phi [Degrees]", "Theta [Degrees]", "post_cuts_theta_phi_scatter.pdf");

	histogram_functions::Plot2DHistogram(h_ele_eng_costheta_nue_cc, "", "Reco Electron Momentum [GeV]", "Reco Electron Cos(#theta)",
	                                     "post_cuts_leading_pfp_eng_costheta_nue_cc.pdf");
	histogram_functions::LegoStackData(h_ele_eng_costheta_nue_cc, h_ele_eng_costheta_nue_cc_mixed, h_ele_eng_costheta_numu_cc, h_ele_eng_costheta_numu_cc_mixed,
	                                   h_ele_eng_costheta_cosmic, h_ele_eng_costheta_nc, h_ele_eng_costheta_nc_pi0, h_ele_eng_costheta_other_mixed,
	                                   h_ele_eng_costheta_unmatched, h_ele_eng_costheta_intime, intime_scale_factor,
	                                   h_ele_eng_costheta_data, data_scale_factor, 0.75, 0.95, 0.70, 0.95,
	                                   "", "Reco Electron Momentum [GeV]", "Reco Electron Cos(#theta)", "post_cuts_leading_pfp_eng_costheta_data.pdf");

	histogram_functions::Plot2DHistogram(h_true_reco_ele_momentum, "Selected True Electrons", "Reco Electron Momentum [GeV]", "True Electron Momentum [GeV]",
	                                     "post_cuts_ele_true_reco_momentum.pdf");
	histogram_functions::Plot2DHistogram(h_true_reco_ele_costheta, "Selected True Electrons", "Reco Electron Cos(#theta) [GeV]", "True Electron Cos(#theta) [GeV]",
	                                     "post_cuts_ele_true_reco_costheta.pdf");
	histogram_functions::Plot1DHistogram(h_true_num_e, "Number of Selected True Electrons Per Event", "post_cuts_num_true_ele.pdf");

	histogram_functions::Plot2DHistogram(h_true_reco_ele_momentum_pre, "Pre Selection True Electrons",
	                                     "Reco Electron Momentum [GeV]", "True Electron Momentum [GeV]",
	                                     "post_cuts_ele_true_reco_momentum_pre.pdf");
	histogram_functions::Plot2DHistogram(h_true_reco_ele_costheta_pre, "Pre Selection True Electrons",
	                                     "Reco Electron Cos(#theta) [GeV]", "True Electron Cos(#theta) [GeV]",
	                                     "post_cuts_ele_true_reco_costheta_pre.pdf");
	histogram_functions::Plot1DHistogram(h_true_num_e_pre, "Number of Pre Selection True Electrons Per Event", "post_cuts_num_true_ele_pre.pdf");

	histogram_functions::PlotSimpleStackData(h_ele_eng_for_nue_cc,        h_ele_eng_for_nue_cc_mixed,
	                                         h_ele_eng_for_nue_cc_out_fv, h_ele_eng_for_numu_cc,   h_ele_eng_for_numu_cc_mixed,
	                                         h_ele_eng_for_cosmic,        h_ele_eng_for_nc,        h_ele_eng_for_nc_pi0,
	                                         h_ele_eng_for_other_mixed,   h_ele_eng_for_unmatched, h_ele_eng_for_intime,
	                                         intime_scale_factor,         h_ele_eng_for_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (Cos(#theta) >= 0.5)", "", "post_cuts_leading_eng_for_data.pdf");

	histogram_functions::PlotSimpleStackData(h_ele_eng_mid_nue_cc,        h_ele_eng_mid_nue_cc_mixed,
	                                         h_ele_eng_mid_nue_cc_out_fv, h_ele_eng_mid_numu_cc,   h_ele_eng_mid_numu_cc_mixed,
	                                         h_ele_eng_mid_cosmic,        h_ele_eng_mid_nc,        h_ele_eng_mid_nc_pi0,
	                                         h_ele_eng_mid_other_mixed,   h_ele_eng_mid_unmatched, h_ele_eng_mid_intime,
	                                         intime_scale_factor,         h_ele_eng_mid_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (0.5 > Cos(#theta) > -0.5)", "", "post_cuts_leading_eng_mid_data.pdf");

	histogram_functions::PlotSimpleStackData(h_ele_eng_back_nue_cc,        h_ele_eng_back_nue_cc_mixed,
	                                         h_ele_eng_back_nue_cc_out_fv, h_ele_eng_back_numu_cc,   h_ele_eng_back_numu_cc_mixed,
	                                         h_ele_eng_back_cosmic,        h_ele_eng_back_nc,        h_ele_eng_back_nc_pi0,
	                                         h_ele_eng_back_other_mixed,   h_ele_eng_back_unmatched, h_ele_eng_back_intime,
	                                         intime_scale_factor,          h_ele_eng_back_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (Cos(#theta) <= -0.5)", "", "post_cuts_leading_eng_back_data.pdf");

	histogram_functions::PlotSimpleStackData(h_ele_eng_for_trans_nue_cc,        h_ele_eng_for_trans_nue_cc_mixed,
	                                         h_ele_eng_for_trans_nue_cc_out_fv, h_ele_eng_for_trans_numu_cc,   h_ele_eng_for_trans_numu_cc_mixed,
	                                         h_ele_eng_for_trans_cosmic,        h_ele_eng_for_trans_nc,        h_ele_eng_for_trans_nc_pi0,
	                                         h_ele_eng_for_trans_other_mixed,   h_ele_eng_for_trans_unmatched, h_ele_eng_for_trans_intime,
	                                         intime_scale_factor,               h_ele_eng_for_trans_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (Cos(#theta_t) >= 0.5)", "", "post_cuts_leading_eng_for_trans_data.pdf");

	histogram_functions::PlotSimpleStackData(h_ele_eng_mid_trans_nue_cc,        h_ele_eng_mid_trans_nue_cc_mixed,
	                                         h_ele_eng_mid_trans_nue_cc_out_fv, h_ele_eng_mid_trans_numu_cc,   h_ele_eng_mid_trans_numu_cc_mixed,
	                                         h_ele_eng_mid_trans_cosmic,        h_ele_eng_mid_trans_nc,        h_ele_eng_mid_trans_nc_pi0,
	                                         h_ele_eng_mid_trans_other_mixed,   h_ele_eng_mid_trans_unmatched, h_ele_eng_mid_trans_intime,
	                                         intime_scale_factor,         h_ele_eng_mid_trans_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (0.5 > Cos(#theta_t) > -0.5)", "", "post_cuts_leading_eng_mid_trans_data.pdf");

	histogram_functions::PlotSimpleStackData(h_ele_eng_back_trans_nue_cc,        h_ele_eng_back_trans_nue_cc_mixed,
	                                         h_ele_eng_back_trans_nue_cc_out_fv, h_ele_eng_back_trans_numu_cc,   h_ele_eng_back_trans_numu_cc_mixed,
	                                         h_ele_eng_back_trans_cosmic,        h_ele_eng_back_trans_nc,        h_ele_eng_back_trans_nc_pi0,
	                                         h_ele_eng_back_trans_other_mixed,   h_ele_eng_back_trans_unmatched, h_ele_eng_back_trans_intime,
	                                         intime_scale_factor,          h_ele_eng_back_trans_data,      data_scale_factor,
	                                         "", "Reco Electron Energy [GeV] (Cos(#theta_t) <= -0.5)", "", "post_cuts_leading_eng_back_trans_data.pdf");

	histogram_functions::Plot2DHistogram(h_mc_vtx_xy_nue_cc,     "", "True Signal Nue Vtx X [cm]", "True Signal Nue Vtx Y [cm]", "true_vtx_xy_nue_cc.pdf"    );
	histogram_functions::Plot2DHistogram(h_mc_vtx_xz_nue_cc,     "", "True Signal Nue Vtx X [cm]", "True Signal Nue Vtx Z [cm]", "true_vtx_xz_nue_cc.pdf"    );
	histogram_functions::Plot2DHistogram(h_mc_vtx_yz_nue_cc,     "", "True Signal Nue Vtx Z [cm]", "True Signal Nue Vtx Y [cm]", "true_vtx_yz_nue_cc.pdf"    );
	histogram_functions::Plot2DHistogram(h_reco_vtx_xy_nue_cc,   "", "Reco Signal Nue Vtx X [cm]", "Reco Signal Nue Vtx Y [cm]", "reco_vtx_xy_nue_cc.pdf"    );
	histogram_functions::Plot2DHistogram(h_reco_vtx_xz_nue_cc,   "", "Reco Signal Nue Vtx X [cm]", "Reco Signal Nue Vtx Z [cm]", "reco_vtx_xz_nue_cc.pdf"    );
	histogram_functions::Plot2DHistogram(h_reco_vtx_yz_nue_cc,   "", "Reco Signal Nue Vtx Z [cm]", "Reco Signal Nue Vtx Y [cm]", "reco vtx_yz_nue_cc.pdf"    );
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_x_nue_cc, "", "True Signal Nue Vtx X [cm]", "Reco Signal Nue Vtx X [cm]", "true_reco_vtx_x_nue_cc.pdf");
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_y_nue_cc, "", "True Signal Nue Vtx Y [cm]", "Reco Signal Nue Vtx Y [cm]", "true_reco_vtx_y_nue_cc.pdf");
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_z_nue_cc, "", "True Signal Nue Vtx Z [cm]", "Reco Signal Nue Vtx Z [cm]", "true_reco_vtx_z_nue_cc.pdf");

	histogram_functions::Plot2DHistogram(h_mc_vtx_xy_nue_cc_out_fv,
	                                     "", "True OutFV Nue Vtx X [cm]", "True OutFV Nue Vtx Y [cm]", "true_vtx_xy_nue_cc_out_fv.pdf"    );
	histogram_functions::Plot2DHistogram(h_mc_vtx_xz_nue_cc_out_fv,
	                                     "", "True OutFV Nue Vtx X [cm]", "True OutFV Nue Vtx Z [cm]", "true_vtx_xz_nue_cc_out_fv.pdf"    );
	histogram_functions::Plot2DHistogram(h_mc_vtx_yz_nue_cc_out_fv,
	                                     "", "True OutFV Nue Vtx Z [cm]", "True OutFV Nue Vtx Y [cm]", "true_vtx_yz_nue_cc_out_fv.pdf"    );
	histogram_functions::Plot2DHistogram(h_reco_vtx_xy_nue_cc_out_fv,
	                                     "", "Reco OutFV Nue Vtx X [cm]", "Reco OutFV Nue Vtx Y [cm]", "reco_vtx_xy_nue_cc_out_fv.pdf"    );
	histogram_functions::Plot2DHistogram(h_reco_vtx_xz_nue_cc_out_fv,
	                                     "", "Reco OutFV Nue Vtx X [cm]", "Reco OutFV Nue Vtx Z [cm]", "reco_vtx_xz_nue_cc_out_fv.pdf"    );
	histogram_functions::Plot2DHistogram(h_reco_vtx_yz_nue_cc_out_fv,
	                                     "", "Reco OutFV Nue Vtx Z [cm]", "Reco OutFV Nue Vtx Y [cm]", "reco vtx_yz_nue_cc_out_fv.pdf"    );
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_x_nue_cc_out_fv,
	                                     "", "True OutFV Nue Vtx X [cm]", "Reco OutFV Nue Vtx X [cm]", "true_reco_vtx_x_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_y_nue_cc_out_fv,
	                                     "", "True OutFV Nue Vtx Y [cm]", "Reco OutFV Nue Vtx Y [cm]", "true_reco_vtx_y_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogram(h_mc_reco_vtx_z_nue_cc_out_fv,
	                                     "", "True OutFV Nue Vtx Z [cm]", "Reco OutFV Nue Vtx Z [cm]", "true_reco_vtx_z_nue_cc_out_fv.pdf");


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
	failure_reason_stack_c1->Print("failure_reason_stack.pdf");


	std::cout << " --- End Cross Section Calculation --- " << std::endl;

	if(f->IsOpen()) {f->Close(); }
	//for some reason, the histogram is not being deleted upon function exit
	//trying to close root file instead
	//delete h_nue_eng_eff_den;
	//delete h_nue_eng_eff_num;

	return 0;
}//end selection
}//end namespace

#ifndef __ROOTCLING__

int main(int argc, char *argv[]){

	const char * file1 = argv[1];
	const char * file2 = argv[2];
	const char * file3 = argv[3];

	std::cout << "INPUT FORMAT: MC_FILE INTIME_FILE DATA_FILE" << std::endl;

	if(argc != 3 && argc != 4)
	{
		std::cout << "Running without in-time cosmics " << std::endl;
		std::cout << "Running without data" << std::endl;
		return xsecSelection::selection(file1, "empty", "empty");
	}
	if(argc != 4 )
	{
		std::cout << "Running without in-time data " << std::endl;
		return xsecSelection::selection(file1, file2, "empty");
	}
	if(argc < 2 )  { std::cout << "Please inclue the input file path" << std::endl; exit(1); }
	return xsecSelection::selection(file1, file2, file3);
}

#endif
