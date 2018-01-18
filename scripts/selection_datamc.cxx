#include "selection_datamc.h"

namespace xsecSelection {
int selection( const char * _file1, const char * _file2){

	std::cout << "File Path: " << _file1 << std::endl;
	std::cout << "File Path: " << _file2 << std::endl;
	const bool _verbose = false;
	const bool _post_cuts_verbose = false;
	const bool data_post_cuts_verbose = false;
	//first we need to open the root file
	TFile * f = new TFile(_file1);
	if(!f->IsOpen()) {std::cout << "Could not open file: " << _file1 << "!" << std::endl; exit(1); }
	TTree * mytree = (TTree*)f->Get("AnalyzeTPCO/tree");
	TTree * optree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");
	TTree * mctree = (TTree*)f->Get("AnalyzeTPCO/mcparticle_tree");
	TTree * mctruth_counter_tree = (TTree*)f->Get("AnalyzeTPCO/mctruth_counter_tree");

	std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
	mytree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

	TFile * f_data = new TFile(_file2);
	if(!f_data->IsOpen()) {std::cout << "Could not open: " << _file2 << "!" << std::endl; exit(1); }
	TTree * data_tree   = (TTree*)f_data->Get("AnalyzeTPCO/tree");
	TTree * data_optree = (TTree*)f_data->Get("AnalyzeTPCO/optical_tree");

	std::vector<xsecAna::TPCObjectContainer> * data_tpc_object_container_v = nullptr;
	data_tree->SetBranchAddress("TpcObjectContainerV", &data_tpc_object_container_v);

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


	const int total_mc_entries = mctruth_counter_tree->GetEntries();
	std::cout << "Total MC Entries: " << total_mc_entries << std::endl;
	mctruth_counter_tree->GetEntry(total_mc_entries-1);

	std::cout << "MC Nue CC Counter      : " << mc_nue_cc_counter << std::endl;
	std::cout << "MC Nue NC Counter      : " << mc_nue_nc_counter << std::endl;
	std::cout << "MC Numu CC Counter     : " << mc_numu_cc_counter << std::endl;
	std::cout << "MC Numu NC Counter     : " << mc_numu_nc_counter << std::endl;
	std::cout << "MC Nue CC Counter Bar  : " << mc_nue_cc_counter_bar << std::endl;
	std::cout << "MC Nue NC Counter Bar  : " << mc_nue_nc_counter_bar << std::endl;
	std::cout << "MC Numu CC Counter Bar : " << mc_numu_cc_counter_bar << std::endl;
	std::cout << "MC Numu NC Counter Bar : " << mc_numu_nc_counter_bar << std::endl;
	std::cout << "----------------------------------" << std::endl;

	std::cout << "=====================" << std::endl;
	std::cout << "== Begin Selection ==" << std::endl;
	std::cout << "=====================" << std::endl;

	std::cout << "======================" << std::endl;
	std::cout << "== First We Do Data ==" << std::endl;
	std::cout << "======================" << std::endl;

	std::vector<int> * data_in_time_counter_v          = new std::vector<int>;
	data_in_time_counter_v->resize(22, 0);
	std::vector<int> * data_pe_counter_v               = new std::vector<int>;
	data_pe_counter_v->resize(22, 0);
	std::vector<int> * data_reco_nue_counter_v         = new std::vector<int>;
	data_reco_nue_counter_v->resize(22, 0);
	std::vector<int> * data_in_fv_counter_v            = new std::vector<int>;
	data_in_fv_counter_v->resize(22, 0);
	std::vector<int> * data_vtx_flash_counter_v        = new std::vector<int>;
	data_vtx_flash_counter_v->resize(22, 0);
	std::vector<int> * data_shwr_tpco_counter_v        = new std::vector<int>;
	data_shwr_tpco_counter_v->resize(22, 0);
	std::vector<int> * data_trk_tpco_counter_v         = new std::vector<int>;
	data_trk_tpco_counter_v->resize(22, 0);
	std::vector<int> * data_hit_threshold_counter_v    = new std::vector<int>;
	data_hit_threshold_counter_v->resize(22, 0);
	std::vector<int> * data_open_angle_counter_v       = new std::vector<int>;
	data_open_angle_counter_v->resize(22, 0);
	std::vector<int> * data_dedx_counter_v             = new std::vector<int>;
	data_dedx_counter_v->resize(22, 0);
	std::vector<int> * data_secondary_shower_counter_v = new std::vector<int>;
	data_secondary_shower_counter_v->resize(22, 0);
	std::vector<int> * data_hit_lengthRatio_counter_v  = new std::vector<int>;
	data_hit_lengthRatio_counter_v->resize(22, 0);

	std::vector<int> * data_has_track = new std::vector<int>;
	data_has_track->resize(2, 0);
	std::vector<int> * data_no_track = new std::vector<int>;
	data_no_track->resize(2, 0);

	//Event, Run, VtxX, VtxY, VtxZ, pass/fail reason
	std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int> > * data_post_cuts_v
	        = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int> >;

	//need a text file with the run, subrun output
	std::ofstream run_subrun_file;
	run_subrun_file.open("run_subrun_list.txt");


	const int data_total_entries = data_tree->GetEntries();
	std::cout << "Data Total Events: " << data_total_entries << std::endl;

	std::vector<int> * data_passed_runs = new std::vector<int>;
	//passed runs is filled with 0, 1, or 2
	//0 = not in time
	//1 = passed - in time and PE threshold
	// 2 = in-time, but not enough PE -- this counts against my efficiency
	data_passed_runs->resize(data_total_entries);

	_cuts_instance.selection_cuts::loop_flashes(f_data, data_optree, flash_pe_threshold, flash_time_start,
	                                            flash_time_end, data_passed_runs);
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
	_cuts_instance.selection_cuts::SetXYflashVector(f_data, data_optree, data_largest_flash_v_v, flash_time_start, flash_time_end);
	std::cout << "-Data- Largest Flash Vector Size: " << data_largest_flash_v_v->size() << std::endl;

	//**********************************
	//now let's do the TPCO related cuts
	//*********************************
	for(int event = 0; event < data_total_entries; event++)
	{
		if(_verbose)
		{
			std::cout << "----------------------" << std::endl;
			std::cout << "[EVENT NUMBER] \t " << event << std::endl;
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
		std::vector<std::string> * tpco_origin_v = new std::vector<std::string>;
		_cuts_instance.selection_cuts::GetOrigins(data_tpc_object_container_v, tpco_origin_v);
		//**********************************************************************
		//Get Run and Subrun*
		//run and subrun number should be the same for all tpc obj in one event
		//**********************************************************************
		//some events have no tpc objects...
		if(data_tpc_object_container_v->size() != 0)
		{
			const int data_run    = data_tpc_object_container_v->at(0).RunNumber();
			const int data_subrun = data_tpc_object_container_v->at(0).SubRunNumber();
			run_subrun_file << data_run << data_subrun << "\n";
		}
		//XY Position of largest flash
		std::vector < double > largest_flash_v = data_largest_flash_v_v->at(event);
		//List of TPC Objects which pass the cuts
		std::vector<std::pair<int, std::string> > * passed_tpco = new std::vector<std::pair<int, std::string> >;
		passed_tpco->resize(data_tpc_object_container_v->size());
		for(int i = 0; i < passed_tpco->size(); i++)
		{
			passed_tpco->at(i).first = 1;
			passed_tpco->at(i).second = "Passed";
		}
		//***********************************************************
		//this is where the in-time optical cut again takes effect
		//***********************************************************
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_in_time_counter_v);

		//PE threshold cut
		if(data_passed_runs->at(event) == 2)
		{
			if(_verbose) std::cout << "[Passed In-Time Cut] [Failed PE Threshold] " << std::endl;
			continue;
		}
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_pe_counter_v);

		//****************************
		// ****** reco nue cut *******
		//****************************
		_cuts_instance.selection_cuts::HasNue(data_tpc_object_container_v, passed_tpco, _verbose);
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_reco_nue_counter_v);

		//** Testing leading shower length vs hits **//
		_data_functions_instance.selection_functions_data::ShowerLengthvsHitsData(data_tpc_object_container_v, passed_tpco, _verbose,
		                                                                          h_shwr_len_hits_data);

		_data_functions_instance.selection_functions_data::TopologyPlots1Data(data_tpc_object_container_v, passed_tpco,
		                                                                      h_pfp_track_shower_data,
		                                                                      h_pfp_track_data,
		                                                                      h_pfp_shower_data);
		//************************
		//******** in fv cut *****
		//************************
		_cuts_instance.selection_cuts::fiducial_volume_cut(data_tpc_object_container_v, _x1, _x2, _y1, _y2, _z1, _z2, passed_tpco, _verbose);
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_in_fv_counter_v);

		//*****************************
		//**** vertex to flash cut ****
		//*****************************
		_data_functions_instance.selection_functions_data::PostCutsVtxFlashData(largest_flash_v, data_tpc_object_container_v, passed_tpco, h_vtx_flash_data);
		_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, data_tpc_object_container_v, tolerance, passed_tpco, _verbose);
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_vtx_flash_counter_v);

		//******************************************************
		//*** distance between pfp shower and nue object cut ***
		//******************************************************
		_data_functions_instance.selection_functions_data::PostCutsShwrVtxData(data_tpc_object_container_v, passed_tpco, _verbose, h_shwr_vtx_dist_data);

		_cuts_instance.selection_cuts::VtxNuDistance(data_tpc_object_container_v, shwr_nue_tolerance, passed_tpco, _verbose);
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_shwr_tpco_counter_v);

		//******************************************************
		// **** distance between pfp track and nue object cut **
		//******************************************************
		_data_functions_instance.selection_functions_data::PostCutTrkVtxData(data_tpc_object_container_v, passed_tpco, _verbose, h_trk_vtx_dist_data);
		_cuts_instance.selection_cuts::VtxTrackNuDistance(data_tpc_object_container_v, trk_nue_tolerance, passed_tpco, _verbose);
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_trk_tpco_counter_v);

		//****************************************************
		// ******** hit threshold for showers cut *************
		//******************************************************
		_cuts_instance.selection_cuts::HitThreshold(data_tpc_object_container_v, shwr_hit_threshold, passed_tpco, _verbose);
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_hit_threshold_counter_v);

		_data_functions_instance.selection_functions_data::dEdxVsOpenAngleData(data_tpc_object_container_v, passed_tpco, _verbose, h_dedx_open_angle_data);
		//*****************************************************
		//****** open angle cut for the leading shower ********
		//******************************************************
		_data_functions_instance.selection_functions_data::PostCutOpenAngleData(data_tpc_object_container_v, passed_tpco, _verbose, h_leading_shower_open_angle_data);
		_cuts_instance.selection_cuts::OpenAngleCut(data_tpc_object_container_v, passed_tpco, tolerance_open_angle, _verbose);
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_open_angle_counter_v);

		//*****************************************************
		//*********** dEdx cut for the leading shower *********
		//******************************************************
		_data_functions_instance.selection_functions_data::PostCutsdEdxData(data_tpc_object_container_v, passed_tpco, _verbose, h_dedx_cuts_data);
		_cuts_instance.selection_cuts::dEdxCut(data_tpc_object_container_v, passed_tpco, tolerance_dedx_min, tolerance_dedx_max, _verbose);
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_dedx_counter_v);

		//*************************************
		// ******** End Selection Cuts! ******
		//*************************************
		_data_functions_instance.selection_functions_data::TopologyPlots2Data(data_tpc_object_container_v, passed_tpco,
		                                                                      h_pfp_track_shower_data_last,
		                                                                      h_pfp_track_data_last,
		                                                                      h_pfp_shower_data_last);

		_data_functions_instance.selection_functions_data::TopologyEfficiencyData(data_tpc_object_container_v, passed_tpco, _verbose,
		                                                                          data_no_track, data_has_track);

		_data_functions_instance.selection_functions_data::FillPostCutVectorData(data_tpc_object_container_v, passed_tpco, data_post_cuts_v);
		_data_functions_instance.selection_functions_data::SecondaryShowersDistData(data_tpc_object_container_v, passed_tpco, _verbose, h_second_shwr_dist_data);
		_data_functions_instance.selection_functions_data::HitLengthRatioData(data_tpc_object_container_v, passed_tpco, _verbose, h_hit_length_ratio_data);

//***************************************************************************
		_cuts_instance.selection_cuts::SecondaryShowersDistCut(data_tpc_object_container_v, passed_tpco, _verbose, dist_tolerance);
		//if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_secondary_shower_counter_v);

		_cuts_instance.selection_cuts::HitLengthRatioCut(data_tpc_object_container_v, passed_tpco, _verbose, pfp_hits_length_tolerance);
		//if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _data_functions_instance.selection_functions_data::TabulateOriginsData(data_tpc_object_container_v, passed_tpco);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, data_hit_lengthRatio_counter_v);

		_data_functions_instance.selection_functions_data::FailureReasonData(data_tpc_object_container_v, passed_tpco, h_failure_reason_data);

	}//end event loop
	run_subrun_file.close();

	//we also want some metrics to print at the end
	_data_functions_instance.selection_functions_data::PrintInfoData( data_in_time_counter_v, "In Time");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_pe_counter_v, "PE Threshold");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_reco_nue_counter_v, "Reco Nue");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_in_fv_counter_v, "In FV");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_vtx_flash_counter_v, "Vtx-to-Flash");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_shwr_tpco_counter_v, "Shower-to-TPCO");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_trk_tpco_counter_v, "Track-to-TPCO");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_hit_threshold_counter_v,"Hit Threshold");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_open_angle_counter_v, "Open Angle");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_dedx_counter_v, " dE / dx ");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_secondary_shower_counter_v, ">3 Shower TPCO Dist");
	_data_functions_instance.selection_functions_data::PrintInfoData( data_hit_lengthRatio_counter_v, "Hit Length Ratio");

	if(data_post_cuts_verbose == true) {_functions_instance.selection_functions::PrintPostCutVector(data_post_cuts_v, data_post_cuts_verbose); }

	//*********************************************
	std::cout << "======================" << std::endl;
	std::cout << "==== Now We Do MC ====" << std::endl;
	std::cout << "======================" << std::endl;

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

	std::vector<int> * has_track = new std::vector<int>;
	has_track->resize(2, 0);
	std::vector<int> * no_track = new std::vector<int>;
	no_track->resize(2, 0);

	const int total_entries = mytree->GetEntries();
	std::cout << "Total Events: " << total_entries << std::endl;

	//Event, Run, VtxX, VtxY, VtxZ, pass/fail reason
	std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int> > * post_cuts_v
	        = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int> >;

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

	_cuts_instance.selection_cuts::loop_flashes(f, optree, flash_pe_threshold, flash_time_start,
	                                            flash_time_end, passed_runs);
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
	std::cout << "Number Events In-Time & > 50 PE: " << run_sum << std::endl;
	std::cout << "Number Events Not In-Time      : " << out_of_time_sum << std::endl;
	std::cout << "Number Events In-Time & < 50 PE: " << low_pe_sum << std::endl;
	std::cout << " -------------------------------------- " << std::endl;

	//get vector with largest flashes y,z positions
	std::vector< std::vector< double> > * largest_flash_v_v = new std::vector < std::vector < double > >;
	_cuts_instance.selection_cuts::SetXYflashVector(f, optree, largest_flash_v_v, flash_time_start, flash_time_end);
	std::cout << "Largest Flash Vector Size: " << largest_flash_v_v->size() << std::endl;

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

		double mc_cos_theta = -999;
		if(mc_nu_momentum != 0) {mc_cos_theta = mc_nu_dir_z / mc_nu_momentum; }
		const double mc_phi       = atan2(mc_nu_dir_y, mc_nu_dir_x);
		double mc_ele_cos_theta = -999;
		if(mc_ele_momentum != 0) {mc_ele_cos_theta = mc_ele_dir_z / mc_ele_momentum; }
		const double mc_ele_phi       = atan2(mc_ele_dir_y, mc_ele_dir_x);
		if(mc_nu_id == 1 || mc_nu_id == 5)
		//if this event is a true nue CC interaction and is inside the FV
		//also include nue_cc_bar as in the tpco classification I use the nue-bar as well
		{
			if(_cuts_instance.selection_cuts::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, _x1, _x2, _y1, _y2, _z1, _z2) == true)
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
				h_nue_num_part_eff_den->Fill(mc_nu_num_particles);
				h_nue_num_chrg_part_eff_den->Fill(mc_nu_num_charged_particles);
				h_nue_cos_theta_eff_den->Fill(mc_cos_theta);
				h_nue_phi_eff_den->Fill(mc_phi * (180 / 3.1415));
				h_ele_cos_theta_eff_den->Fill(mc_ele_cos_theta);
				h_ele_phi_eff_den->Fill(mc_ele_phi * (180 / 3.1415));
				//std::cout << "MC Phi: " << mc_phi * (180 / 3.1415) << std::endl;
				total_mc_entries_inFV++;
			}
		}
		std::vector<std::string> *tpco_origin_v = new std::vector<std::string>;
		_cuts_instance.selection_cuts::GetOrigins(tpc_object_container_v, tpco_origin_v);

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
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, in_time_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   _x1, _x2, _y1, _y2, _z1, _z2,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_no_cut, h_selected_ele_energy_no_cut);
		//PE threshold cut
		if(passed_runs->at(event) == 2)
		{
			if(_verbose) std::cout << "[Passed In-Time Cut] [Failed PE Threshold] " << std::endl;
			continue;
		}
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, pe_counter_v);

		//****************************
		// ****** reco nue cut *******
		//****************************
		_cuts_instance.selection_cuts::HasNue(tpc_object_container_v, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, reco_nue_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   _x1, _x2, _y1, _y2, _z1, _z2,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_reco_nue, h_selected_ele_energy_reco_nue);
		//** Testing flash vs neutrino interaction for origin **
		_functions_instance.selection_functions::FlashTot0(largest_flash_v, mc_nu_time, mc_nu_id, tabulated_origins,
		                                                   _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,h_flash_t0_diff);
		//** Testing leading shower length vs hits **//
		_functions_instance.selection_functions::ShowerLengthvsHits(tpc_object_container_v, passed_tpco, _verbose, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2,
		                                                            mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                            h_shwr_len_hits_nue_cc, h_shwr_len_hits_nue_cc_out_fv,
		                                                            h_shwr_len_hits_nue_cc_mixed, h_shwr_len_hits_numu_cc,
		                                                            h_shwr_len_hits_numu_cc_mixed, h_shwr_len_hits_nc,
		                                                            h_shwr_len_hits_nc_pi0, h_shwr_len_hits_cosmic,
		                                                            h_shwr_len_hits_other_mixed, h_shwr_len_hits_unmatched);

		//pre most cuts hits
		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins.at(0) == 1)
		{
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
			                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, mc_nu_energy, mc_ele_energy,
			                                                             h_shwr_hits_nu_eng, h_shwr_hits_ele_eng);
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
			                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, mc_nu_energy, mc_ele_energy,
			                                                             h_shwr_hits_nu_eng_zoom, h_shwr_hits_ele_eng_zoom);
		}
		_functions_instance.selection_functions::TopologyPlots1(tpc_object_container_v, passed_tpco, has_pi0,
		                                                        _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
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

		// _functions_instance.selection_functions::SecondaryShowersDist(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		//                                                               _x1, _x2, _y1, _y2, _z1, _z2,
		//                                                               mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		//                                                               h_second_shwr_dist_nue_cc, h_second_shwr_dist_nue_cc_out_fv,
		//                                                               h_second_shwr_dist_nue_cc_mixed, h_second_shwr_dist_numu_cc,
		//                                                               h_second_shwr_dist_numu_cc_mixed, h_second_shwr_dist_nc,
		//                                                               h_second_shwr_dist_nc_pi0, h_second_shwr_dist_cosmic,
		//                                                               h_second_shwr_dist_other_mixed, h_second_shwr_dist_unmatched);
		//************************
		//******** in fv cut *****
		//************************
		_cuts_instance.selection_cuts::fiducial_volume_cut(tpc_object_container_v, _x1, _x2, _y1, _y2, _z1, _z2, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, in_fv_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   _x1, _x2, _y1, _y2, _z1, _z2,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_in_fv, h_selected_ele_energy_in_fv);

		//*****************************
		//**** vertex to flash cut ****
		//*****************************
		_functions_instance.selection_functions::PostCutsVtxFlash(largest_flash_v, tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                          _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                          h_vtx_flash_nue_cc, h_vtx_flash_nue_cc_mixed,
		                                                          h_vtx_flash_numu_cc, h_vtx_flash_nc,
		                                                          h_vtx_flash_cosmic, h_vtx_flash_nc_pi0,
		                                                          h_vtx_flash_numu_cc_mixed, h_vtx_flash_other_mixed,
		                                                          h_vtx_flash_unmatched);

		_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, tpc_object_container_v, tolerance, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, vtx_flash_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   _x1, _x2, _y1, _y2, _z1, _z2,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_vtx_flash, h_selected_ele_energy_vtx_flash);

		//******************************************************
		//*** distance between pfp shower and nue object cut ***
		//******************************************************
		_functions_instance.selection_functions::PostCutsShwrVtx(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                         _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                         h_shwr_vtx_dist_nue_cc,
		                                                         h_shwr_vtx_dist_nue_cc_mixed,
		                                                         h_shwr_vtx_dist_numu_cc,
		                                                         h_shwr_vtx_dist_nc,
		                                                         h_shwr_vtx_dist_cosmic,
		                                                         h_shwr_vtx_dist_nc_pi0,
		                                                         h_shwr_vtx_dist_numu_cc_mixed,
		                                                         h_shwr_vtx_dist_other_mixed,
		                                                         h_shwr_vtx_dist_unmatched     );

		_cuts_instance.selection_cuts::VtxNuDistance(tpc_object_container_v, shwr_nue_tolerance, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, shwr_tpco_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   _x1, _x2, _y1, _y2, _z1, _z2,
		                                                                   tabulated_origins,  mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_shwr_vtx, h_selected_ele_energy_shwr_vtx);


		//******************************************************
		// **** distance between pfp track and nue object cut **
		//******************************************************
		_functions_instance.selection_functions::PostCutTrkVtx(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                       _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                       h_trk_vtx_dist_nue_cc, h_trk_vtx_dist_nue_cc_mixed,
		                                                       h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_nc,
		                                                       h_trk_vtx_dist_cosmic, h_trk_vtx_dist_nc_pi0,
		                                                       h_trk_vtx_dist_numu_cc_mixed, h_trk_vtx_dist_other_mixed,
		                                                       h_trk_vtx_dist_unmatched);
		_cuts_instance.selection_cuts::VtxTrackNuDistance(tpc_object_container_v, trk_nue_tolerance, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, trk_tpco_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   _x1, _x2, _y1, _y2, _z1, _z2,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_trk_vtx, h_selected_ele_energy_trk_vtx);

		//****************************************************
		// ******** hit threshold for showers cut *************
		//******************************************************
		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins.at(0) == 1)
		{
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
			                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, mc_nu_energy, mc_ele_energy,
			                                                             h_shwr_hits_nu_eng_last, h_shwr_hits_ele_eng_last);
			_functions_instance.selection_functions::PostCutHitThreshold(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
			                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, mc_nu_energy, mc_ele_energy,
			                                                             h_shwr_hits_nu_eng_zoom_last, h_shwr_hits_ele_eng_zoom_last);
		}

		_cuts_instance.selection_cuts::HitThreshold(tpc_object_container_v, shwr_hit_threshold, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_threshold_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   _x1, _x2, _y1, _y2, _z1, _z2,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_hit_threshold, h_selected_ele_energy_hit_threshold);

		_functions_instance.selection_functions::dEdxVsOpenAngle(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                         _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                         h_dedx_open_angle_nue_cc, h_dedx_open_angle_nue_cc_out_fv,
		                                                         h_dedx_open_angle_nue_cc_mixed, h_dedx_open_angle_numu_cc,
		                                                         h_dedx_open_angle_numu_cc_mixed, h_dedx_open_angle_nc,
		                                                         h_dedx_open_angle_nc_pi0, h_dedx_open_angle_cosmic,
		                                                         h_dedx_open_angle_other_mixed, h_dedx_open_angle_unmatched);
		//*****************************************************
		//****** open angle cut for the leading shower ********
		//******************************************************
		_functions_instance.selection_functions::PostCutOpenAngle(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                          _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                          h_leading_shower_open_angle_nue_cc, h_leading_shower_open_angle_nue_cc_mixed,
		                                                          h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_nc,
		                                                          h_leading_shower_open_angle_cosmic, h_leading_shower_open_angle_nc_pi0,
		                                                          h_leading_shower_open_angle_numu_cc_mixed, h_leading_shower_open_angle_other_mixed,
		                                                          h_leading_shower_open_angle_unmatched);
		_cuts_instance.selection_cuts::OpenAngleCut(tpc_object_container_v, passed_tpco, tolerance_open_angle, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, open_angle_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   _x1, _x2, _y1, _y2, _z1, _z2,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_open_angle, h_selected_ele_energy_open_angle);

		//*****************************************************
		//*********** dEdx cut for the leading shower *********
		//******************************************************
		_functions_instance.selection_functions::PostCutsdEdx(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                      _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                      h_dedx_cuts_nue_cc, h_dedx_cuts_nue_cc_mixed,
		                                                      h_dedx_cuts_nue_cc_out_fv,
		                                                      h_dedx_cuts_numu_cc, h_dedx_cuts_nc,
		                                                      h_dedx_cuts_cosmic, h_dedx_cuts_nc_pi0,
		                                                      h_dedx_cuts_numu_cc_mixed, h_dedx_cuts_other_mixed,
		                                                      h_dedx_cuts_unmatched     );

		_cuts_instance.selection_cuts::dEdxCut(tpc_object_container_v, passed_tpco, tolerance_dedx_min, tolerance_dedx_max, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, dedx_counter_v);
		_functions_instance.selection_functions::SequentialTrueEnergyPlots(mc_nu_id, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                                   _x1, _x2, _y1, _y2, _z1, _z2,
		                                                                   tabulated_origins, mc_nu_energy, mc_ele_energy,
		                                                                   h_selected_nu_energy_dedx, h_selected_ele_energy_dedx);
		//*************************************
		// ******** End Selection Cuts! ******
		//*************************************

		//these are for the tefficiency plots, post all cuts
		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins.at(0) == 1)
		{
			if(_cuts_instance.selection_cuts::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
			                                        _x1, _x2, _y1,
			                                        _y2, _z1, _z2) == true)
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
				h_nue_num_part_eff_num->Fill(mc_nu_num_particles);
				h_nue_num_chrg_part_eff_num->Fill(mc_nu_num_charged_particles);
				h_nue_cos_theta_eff_num->Fill(mc_cos_theta);
				h_nue_phi_eff_num->Fill(mc_phi);
				h_ele_cos_theta_eff_num->Fill(mc_ele_cos_theta);
				h_ele_phi_eff_num->Fill(mc_ele_phi);
			}
		}
		_functions_instance.selection_functions::TopologyPlots2(tpc_object_container_v, passed_tpco, has_pi0,
		                                                        _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
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

		_functions_instance.selection_functions::TopologyEfficiency(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                            _x1, _x2, _y1, _y2, _z1, _z2,
		                                                            mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                            no_track, has_track);

		_functions_instance.selection_functions::ChargeShare(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                     _x1, _x2, _y1, _y2, _z1, _z2,
		                                                     mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                     h_charge_share_nue_cc_mixed);

		_functions_instance.selection_functions::FillPostCutVector(tpc_object_container_v, passed_tpco, has_pi0,
		                                                           _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, post_cuts_v);

		_functions_instance.selection_functions::SecondaryShowersDist(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                              _x1, _x2, _y1, _y2, _z1, _z2,
		                                                              mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                              h_second_shwr_dist_nue_cc, h_second_shwr_dist_nue_cc_out_fv,
		                                                              h_second_shwr_dist_nue_cc_mixed, h_second_shwr_dist_numu_cc,
		                                                              h_second_shwr_dist_numu_cc_mixed, h_second_shwr_dist_nc,
		                                                              h_second_shwr_dist_nc_pi0, h_second_shwr_dist_cosmic,
		                                                              h_second_shwr_dist_other_mixed, h_second_shwr_dist_unmatched);

		_functions_instance.selection_functions::HitLengthRatio(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                        _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
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

//***************************************************************************
		_cuts_instance.selection_cuts::SecondaryShowersDistCut(tpc_object_container_v, passed_tpco, _verbose, dist_tolerance);
		//if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, secondary_shower_counter_v);




		_cuts_instance.selection_cuts::HitLengthRatioCut(tpc_object_container_v, passed_tpco, _verbose, pfp_hits_length_tolerance);
		//if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_lengthRatio_counter_v);

		_functions_instance.selection_functions::FailureReason(tpc_object_container_v, passed_tpco, has_pi0, _verbose, _x1, _x2, _y1, _y2, _z1, _z2,
		                                                       mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                       h_failure_reason_nue_cc, h_failure_reason_nue_cc_out_fv,
		                                                       h_failure_reason_nue_cc_mixed, h_failure_reason_numu_cc,
		                                                       h_failure_reason_numu_cc_mixed, h_failure_reason_nc,
		                                                       h_failure_reason_nc_pi0, h_failure_reason_cosmic,
		                                                       h_failure_reason_other_mixed, h_failure_reason_unmatched);

	}//end event loop

	std::cout << "------------------ " << std::endl;
	std::cout << " MC Entries in FV: " << total_mc_entries_inFV << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "  End Selection    " << std::endl;
	std::cout << "------------------ " << std::endl;

	//change mc_nue_cc_counter to total_mc_entries_inFV once files are ready!

	//we also want some metrics to print at the end
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, in_time_counter_v, "In Time");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, pe_counter_v, "PE Threshold");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, reco_nue_counter_v, "Reco Nue");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, in_fv_counter_v, "In FV");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, vtx_flash_counter_v, "Vtx-to-Flash");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, shwr_tpco_counter_v, "Shower-to-TPCO");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, trk_tpco_counter_v, "Track-to-TPCO");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, hit_threshold_counter_v,"Hit Threshold");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, open_angle_counter_v, "Open Angle");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, dedx_counter_v, " dE / dx ");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, secondary_shower_counter_v, ">3 Shower TPCO Dist");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, hit_lengthRatio_counter_v, "Hit Length Ratio");

	std::cout << "---------------------" << std::endl;
	std::cout << "No Track Signal: " << no_track->at(0) << std::endl;
	std::cout << "No Track Bkg   : " << no_track->at(1) << std::endl;
	std::cout << "Purity         : " << double(no_track->at(0)) / double(no_track->at(0) + no_track->at(1)) << std::endl;
	std::cout << " ******************* " << std::endl;
	std::cout << "1+ Track Signal: " << has_track->at(0) << std::endl;
	std::cout << "1+ Track Bkg   : " << has_track->at(1) << std::endl;
	std::cout << "Purity         : " << double(has_track->at(0)) / double(has_track->at(0) + has_track->at(1)) << std::endl;
	std::cout << "---------------------" << std::endl;

	std::vector<double> * xsec_cc = new std::vector<double>;
	const double final_counter                = dedx_counter_v->at(7);
	const double final_counter_nue_cc         = dedx_counter_v->at(0);
	const double final_counter_nue_cc_mixed   = dedx_counter_v->at(1);
	const double final_counter_nue_cc_out_fv  = dedx_counter_v->at(9);
	const double final_counter_cosmic         = dedx_counter_v->at(2);
	const double final_counter_nc         = dedx_counter_v->at(3);
	const double final_counter_numu_cc        = dedx_counter_v->at(4);
	const double final_counter_numu_cc_mixed  = dedx_counter_v->at(11);
	const double final_counter_nc_pi0        = dedx_counter_v->at(10);
	const double final_counter_unmatched      = dedx_counter_v->at(5);
	const double final_counter_other_mixed    = dedx_counter_v->at(6);
	const int n_total = final_counter;
	const int n_bkg = (final_counter_nue_cc_mixed + final_counter_nue_cc_out_fv + final_counter_cosmic + final_counter_nc
	                   + final_counter_numu_cc + final_counter_numu_cc_mixed + final_counter_nc_pi0 + final_counter_unmatched + final_counter_other_mixed);
	const double efficiency = final_counter_nue_cc / double(total_mc_entries_inFV);
	_functions_instance.selection_functions::calcXSec(_x1, _x2, _y1, _y2, _z1, _z2,
	                                                  n_total, n_bkg, flux,
	                                                  efficiency, xsec_cc);

	std::cout << "-------------------------" << std::endl;
	std::cout << " Cross Section Results:  " << std::endl;
	std::cout << " " << xsec_cc->at(0) << " +/- (stats) "
	          << xsec_cc->at(1) << " +/- (sys) "
	          << xsec_cc->at(2) << std::endl;
	std::cout << "-------------------------" << std::endl;
	xsec_cc->clear();
	_functions_instance.selection_functions::calcXSec(_x1, _x2, _y1, _y2, _z1, _z2,
	                                                  total_mc_entries_inFV, 0, flux,
	                                                  1, xsec_cc);
	std::cout << "-------------------------" << std::endl;
	std::cout << " Cross Section Results (Truth):  " << std::endl;
	std::cout << " " << xsec_cc->at(0) << " +/- (stats) "
	          << xsec_cc->at(1) << " +/- (sys) "
	          << xsec_cc->at(2) << std::endl;
	std::cout << "-------------------------" << std::endl;
	std::cout << "-------------------------" << std::endl;
	std::cout << " Genie value of Flux " << '\n' <<
	        " Integrated Xsec:    " << genie_xsec << std::endl;

	if(_post_cuts_verbose == true) {_functions_instance.selection_functions::PrintPostCutVector(post_cuts_v, _post_cuts_verbose); }

//********************//
//**** Histograms ****//
//*******************//

	TCanvas * test_c1 = new TCanvas();
	test_c1->cd();
	h_nue_eng_eff_num->GetXaxis()->SetTitle("Selected True Neutrino Energy [GeV]");
	h_nue_eng_eff_num->Draw();
	test_c1->Print("selected_true_neutrino_energy.pdf");
	TCanvas * test_c2 = new TCanvas();
	test_c2->cd();
	h_nue_eng_eff_den->GetXaxis()->SetTitle("True Neutrino Energy [GeV]");
	h_nue_eng_eff_den->GetYaxis()->SetTitle("Events");
	h_nue_eng_eff_den->Draw();
	test_c2->Print("all_true_neutrino_energy.pdf");
	TCanvas * num_part_c1 = new TCanvas();
	num_part_c1->cd();
	h_nue_num_part_eff_den->GetXaxis()->SetTitle("True Particle Multiplicity");
	h_nue_num_part_eff_den->Draw();
	num_part_c1->Print("all_true_neutrino_num_particles.pdf");
	TCanvas * num_part_c2 = new TCanvas();
	num_part_c2->cd();
	h_nue_num_part_eff_num->GetXaxis()->SetTitle("Selected True Particle Multiplicity");
	h_nue_num_part_eff_num->Draw();
	num_part_c2->Print("selected_true_neutrino_num_particles.pdf");
	TCanvas * num_part_c3 = new TCanvas();
	num_part_c3->cd();
	h_nue_num_chrg_part_eff_den->GetXaxis()->SetTitle("True Charged Particle Multiplicity");
	h_nue_num_chrg_part_eff_den->Draw();
	num_part_c3->Print("all_true_neutrino_num_charged_particles.pdf");
	TCanvas * num_part_c4 = new TCanvas();
	num_part_c4->cd();
	h_nue_num_chrg_part_eff_num->GetXaxis()->SetTitle("Selected Charged True Particle Multiplicity");
	h_nue_num_chrg_part_eff_num->Draw();
	num_part_c4->Print("selected_true_neutrino_num_charged_particles.pdf");


	TCanvas * efficiency_c1 = new TCanvas();
	efficiency_c1->cd();
	TEfficiency * eng_eff = new TEfficiency(*h_nue_eng_eff_num, *h_nue_eng_eff_den);
	eng_eff->SetTitle(";True Neutrino Energy [GeV];Efficiency");
	eng_eff->SetLineColor(kGreen+3);
	eng_eff->SetMarkerColor(kGreen+3);
	eng_eff->SetMarkerStyle(20);
	eng_eff->SetMarkerSize(0.5);
	eng_eff->Draw("AP");
	efficiency_c1->Print("signal_selection_nu_energy_efficiency.pdf");

	TCanvas * efficiency_c2 = new TCanvas();
	efficiency_c2->cd();
	TEfficiency * vtx_x_eff = new TEfficiency(*h_nue_vtx_x_eff_num, *h_nue_vtx_x_eff_den);
	vtx_x_eff->SetTitle(";True Neutrino Vtx X [cm];Efficiency");
	vtx_x_eff->SetLineColor(kGreen+3);
	vtx_x_eff->SetMarkerColor(kGreen+3);
	vtx_x_eff->SetMarkerStyle(20);
	vtx_x_eff->SetMarkerSize(0.5);
	vtx_x_eff->Draw("AP");
	efficiency_c2->Print("signal_selection_nu_vtx_x_efficiency.pdf");

	TCanvas * efficiency_c3 = new TCanvas();
	efficiency_c3->cd();
	TEfficiency * vtx_y_eff = new TEfficiency(*h_nue_vtx_y_eff_num, *h_nue_vtx_y_eff_den);
	vtx_y_eff->SetTitle(";True Neutrino Vtx Y [cm];Efficiency");
	vtx_y_eff->SetLineColor(kGreen+3);
	vtx_y_eff->SetMarkerColor(kGreen+3);
	vtx_y_eff->SetMarkerStyle(20);
	vtx_y_eff->SetMarkerSize(0.5);
	vtx_y_eff->Draw("AP");
	efficiency_c3->Print("signal_selection_nu_vtx_y_efficiency.pdf");

	TCanvas * efficiency_c4 = new TCanvas();
	efficiency_c4->cd();
	TEfficiency * vtx_z_eff = new TEfficiency(*h_nue_vtx_z_eff_num, *h_nue_vtx_z_eff_den);
	vtx_z_eff->SetTitle(";True Neutrino Vtx Z [cm];Efficiency");
	vtx_z_eff->SetLineColor(kGreen+3);
	vtx_z_eff->SetMarkerColor(kGreen+3);
	vtx_z_eff->SetMarkerStyle(20);
	vtx_z_eff->SetMarkerSize(0.5);
	vtx_z_eff->Draw("AP");
	efficiency_c4->Print("signal_selection_nu_vtx_z_efficiency.pdf");

	TCanvas * efficiency_c5 = new TCanvas();
	efficiency_c5->cd();
	TEfficiency * dir_x_eff = new TEfficiency(*h_nue_dir_x_eff_num, *h_nue_dir_x_eff_den);
	dir_x_eff->SetTitle(";True Neutrino Dir X;Efficiency");
	dir_x_eff->SetLineColor(kGreen+3);
	dir_x_eff->SetMarkerColor(kGreen+3);
	dir_x_eff->SetMarkerStyle(20);
	dir_x_eff->SetMarkerSize(0.5);
	dir_x_eff->Draw("AP");
	efficiency_c5->Print("signal_selection_nu_dir_x_efficiency.pdf");
	TCanvas * efficiency_c6 = new TCanvas();
	efficiency_c6->cd();
	TEfficiency * dir_y_eff = new TEfficiency(*h_nue_dir_y_eff_num, *h_nue_dir_y_eff_den);
	dir_y_eff->SetTitle(";True Neutrino Dir Y;Efficiency");
	dir_y_eff->SetLineColor(kGreen+3);
	dir_y_eff->SetMarkerColor(kGreen+3);
	dir_y_eff->SetMarkerStyle(20);
	dir_y_eff->SetMarkerSize(0.5);
	dir_y_eff->Draw("AP");
	efficiency_c6->Print("signal_selection_nu_dir_y_efficiency.pdf");
	TCanvas * efficiency_c7 = new TCanvas();
	efficiency_c7->cd();
	TEfficiency * dir_z_eff = new TEfficiency(*h_nue_dir_z_eff_num, *h_nue_dir_z_eff_den);
	dir_z_eff->SetTitle(";True Neutrino Dir Z;Efficiency");
	dir_z_eff->SetLineColor(kGreen+3);
	dir_z_eff->SetMarkerColor(kGreen+3);
	dir_z_eff->SetMarkerStyle(20);
	dir_z_eff->SetMarkerSize(0.5);
	dir_z_eff->Draw("AP");
	efficiency_c7->Print("signal_selection_nu_dir_z_efficiency.pdf");
	TCanvas * efficiency_c8 = new TCanvas();
	efficiency_c8->cd();
	TEfficiency * num_particles_eff = new TEfficiency(*h_nue_num_part_eff_num, *h_nue_num_part_eff_den);
	num_particles_eff->SetTitle(";True Particle Multiplicity;Efficiency");
	num_particles_eff->SetLineColor(kGreen+3);
	num_particles_eff->SetMarkerColor(kGreen+3);
	num_particles_eff->SetMarkerStyle(20);
	num_particles_eff->SetMarkerSize(0.5);
	num_particles_eff->Draw("AP");
	efficiency_c8->Print("signal_selection_nu_num_particles_efficiency.pdf");
	TCanvas * efficiency_c9 = new TCanvas();
	efficiency_c9->cd();
	TEfficiency * num_charged_particles_eff = new TEfficiency(*h_nue_num_chrg_part_eff_num, *h_nue_num_chrg_part_eff_den);
	num_charged_particles_eff->SetTitle(";True Charged Particle Multiplicity;Efficiency");
	num_charged_particles_eff->SetLineColor(kGreen+3);
	num_charged_particles_eff->SetMarkerColor(kGreen+3);
	num_charged_particles_eff->SetMarkerStyle(20);
	num_charged_particles_eff->SetMarkerSize(0.5);
	num_charged_particles_eff->Draw("AP");
	efficiency_c9->Print("signal_selection_nu_num_charged_particles_efficiency.pdf");
	TCanvas * efficiency_c10 = new TCanvas();
	efficiency_c10->cd();
	TEfficiency * cos_theta_eff = new TEfficiency(*h_nue_cos_theta_eff_num, *h_nue_cos_theta_eff_den);
	cos_theta_eff->SetTitle(";True Neutrino Cos(#theta);Efficiency");
	cos_theta_eff->SetLineColor(kGreen+3);
	cos_theta_eff->SetMarkerColor(kGreen+3);
	cos_theta_eff->SetMarkerStyle(20);
	cos_theta_eff->SetMarkerSize(0.5);
	cos_theta_eff->Draw("AP");
	efficiency_c10->Print("signal_selection_nu_cos_theta_efficiency.pdf");
	TCanvas * efficiency_c11 = new TCanvas();
	efficiency_c11->cd();
	TEfficiency * phi_eff = new TEfficiency(*h_nue_phi_eff_num, *h_nue_phi_eff_den);
	phi_eff->SetTitle(";True Neutrino Phi;Efficiency");
	phi_eff->SetLineColor(kGreen+3);
	phi_eff->SetMarkerColor(kGreen+3);
	phi_eff->SetMarkerStyle(20);
	phi_eff->SetMarkerSize(0.5);
	phi_eff->Draw("AP");
	efficiency_c11->Print("signal_selection_nu_phi_efficiency.pdf");
	TCanvas * efficiency_c12 = new TCanvas();
	efficiency_c12->cd();
	TEfficiency * cos_ele_theta_eff = new TEfficiency(*h_ele_cos_theta_eff_num, *h_ele_cos_theta_eff_den);
	cos_ele_theta_eff->SetTitle(";True Nue Electron Cos(#theta);Efficiency");
	cos_ele_theta_eff->SetLineColor(kGreen+3);
	cos_ele_theta_eff->SetMarkerColor(kGreen+3);
	cos_ele_theta_eff->SetMarkerStyle(20);
	cos_ele_theta_eff->SetMarkerSize(0.5);
	cos_ele_theta_eff->Draw("AP");
	efficiency_c12->Print("signal_selection_ele_cos_theta_efficiency.pdf");
	TCanvas * efficiency_c13 = new TCanvas();
	efficiency_c13->cd();
	TEfficiency * ele_phi_eff = new TEfficiency(*h_ele_phi_eff_num, *h_ele_phi_eff_den);
	ele_phi_eff->SetTitle(";True Nue Electron Phi;Efficiency");
	ele_phi_eff->SetLineColor(kGreen+3);
	ele_phi_eff->SetMarkerColor(kGreen+3);
	ele_phi_eff->SetMarkerStyle(20);
	ele_phi_eff->SetMarkerSize(0.5);
	ele_phi_eff->Draw("AP");
	efficiency_c13->Print("signal_selection_ele_phi_efficiency.pdf");
	TCanvas * efficiency_c14 = new TCanvas();
	efficiency_c14->cd();
	TEfficiency * ele_dir_x_eff = new TEfficiency(*h_ele_dir_x_eff_num, *h_ele_dir_x_eff_den);
	ele_dir_x_eff->SetTitle(";True Nue Electron Dir X;Efficiency");
	ele_dir_x_eff->SetLineColor(kGreen+3);
	ele_dir_x_eff->SetMarkerColor(kGreen+3);
	ele_dir_x_eff->SetMarkerStyle(20);
	ele_dir_x_eff->SetMarkerSize(0.5);
	ele_dir_x_eff->Draw("AP");
	efficiency_c14->Print("signal_selection_ele_dir_x_efficiency.pdf");
	TCanvas * efficiency_c15 = new TCanvas();
	efficiency_c15->cd();
	TEfficiency * ele_dir_y_eff = new TEfficiency(*h_ele_dir_y_eff_num, *h_ele_dir_y_eff_den);
	ele_dir_y_eff->SetTitle(";True Nue Electron Dir Y;Efficiency");
	ele_dir_y_eff->SetLineColor(kGreen+3);
	ele_dir_y_eff->SetMarkerColor(kGreen+3);
	ele_dir_y_eff->SetMarkerStyle(20);
	ele_dir_y_eff->SetMarkerSize(0.5);
	ele_dir_y_eff->Draw("AP");
	efficiency_c15->Print("signal_selection_ele_dir_y_efficiency.pdf");
	TCanvas * efficiency_c16 = new TCanvas();
	efficiency_c16->cd();
	TEfficiency * ele_dir_z_eff = new TEfficiency(*h_ele_dir_z_eff_num, *h_ele_dir_z_eff_den);
	ele_dir_z_eff->SetTitle(";True Nue Electron Dir Z;Efficiency");
	ele_dir_z_eff->SetLineColor(kGreen+3);
	ele_dir_z_eff->SetMarkerColor(kGreen+3);
	ele_dir_z_eff->SetMarkerStyle(20);
	ele_dir_z_eff->SetMarkerSize(0.5);
	ele_dir_z_eff->Draw("AP");
	efficiency_c16->Print("signal_selection_ele_dir_z_efficiency.pdf");
	TCanvas * efficiency_c17 = new TCanvas();
	efficiency_c17->cd();
	TEfficiency * ele_eng_eff = new TEfficiency(*h_ele_eng_eff_num, *h_ele_eng_eff_den);
	ele_eng_eff->SetTitle(";True Nue Electron Energy [GeV];Efficiency");
	ele_eng_eff->SetLineColor(kGreen+3);
	ele_eng_eff->SetMarkerColor(kGreen+3);
	ele_eng_eff->SetMarkerStyle(20);
	ele_eng_eff->SetMarkerSize(0.5);
	ele_eng_eff->Draw("AP");
	efficiency_c17->Print("signal_selection_ele_energy_efficiency.pdf");



	double all_energy = 0;
	for(auto const energy : selected_energy_vector) {all_energy += energy; }
	const double average_true_energy = all_energy / selected_energy_vector.size();
	_functions_instance.selection_functions::xsec_plot(_verbose, genie_xsec, xsec_cc->at(0), average_true_energy, xsec_cc->at(1));

	TCanvas * post_cuts_c1 = new TCanvas();
	post_cuts_c1->cd();
	h_tracks_showers->GetYaxis()->SetTitle("Reco Showers");
	h_tracks_showers->GetXaxis()->SetTitle("Reco Tracks");
	h_tracks_showers->SetTitle("Post Cuts - Showers/Tracks per Candidate Nue TPC Object");
	h_tracks_showers->SetStats(kFALSE);
	h_tracks_showers->Draw("colz");
	post_cuts_c1->Print("post_cuts_showers_tracks.pdf");

	TCanvas * post_cuts_c2 = new TCanvas();
	post_cuts_c2->cd();
	h_tracks_showers_cosmic->GetYaxis()->SetTitle("Reco Showers");
	h_tracks_showers_cosmic->GetXaxis()->SetTitle("Reco Tracks");
	h_tracks_showers_cosmic->SetTitle("Post Cuts - Showers/Tracks per Candidate Nue TPC Object");
	h_tracks_showers_cosmic->SetStats(kFALSE);
	h_tracks_showers_cosmic->Draw("colz");
	post_cuts_c2->Print("post_cuts_showers_tracks_cosmic.pdf");

	TCanvas * post_cuts_c3 = new TCanvas();
	post_cuts_c3->cd();
	h_tracks_showers_numu->GetYaxis()->SetTitle("Reco Showers");
	h_tracks_showers_numu->GetXaxis()->SetTitle("Reco Tracks");
	h_tracks_showers_numu->SetTitle("Post Cuts - Showers/Tracks per Candidate Nue TPC Object");
	h_tracks_showers_numu->SetStats(kFALSE);
	h_tracks_showers_numu->Draw("colz");
	post_cuts_c3->Print("post_cuts_showers_tracks_numu.pdf");

	TCanvas * open_angle_stack_c1 = new TCanvas();
	open_angle_stack_c1->cd();
	THStack * open_angle_stack = new THStack();
	h_leading_shower_open_angle_nue_cc->SetStats(kFALSE);
	h_leading_shower_open_angle_nue_cc_mixed->SetStats(kFALSE);
	h_leading_shower_open_angle_numu_cc->SetStats(kFALSE);
	h_leading_shower_open_angle_nc_pi0->SetStats(kFALSE);
	h_leading_shower_open_angle_cosmic->SetStats(kFALSE);
	h_leading_shower_open_angle_nc->SetStats(kFALSE);
	h_leading_shower_open_angle_numu_cc_mixed->SetStats(kFALSE);
	h_leading_shower_open_angle_other_mixed->SetStats(kFALSE);
	h_leading_shower_open_angle_unmatched->SetStats(kFALSE);
	h_leading_shower_open_angle_nue_cc->SetFillColor(30);
	h_leading_shower_open_angle_nue_cc_mixed->SetFillColor(38);
	h_leading_shower_open_angle_numu_cc->SetFillColor(28);
	h_leading_shower_open_angle_nc_pi0->SetFillColor(36);
	h_leading_shower_open_angle_cosmic->SetFillColor(1);
	h_leading_shower_open_angle_nc->SetFillColor(46);
	h_leading_shower_open_angle_numu_cc_mixed->SetFillColor(20);
	h_leading_shower_open_angle_other_mixed->SetFillColor(42);
	h_leading_shower_open_angle_unmatched->SetFillColor(12);
	open_angle_stack->Add(h_leading_shower_open_angle_nue_cc);
	open_angle_stack->Add(h_leading_shower_open_angle_nue_cc_mixed);
	open_angle_stack->Add(h_leading_shower_open_angle_cosmic);
	open_angle_stack->Add(h_leading_shower_open_angle_numu_cc);
	open_angle_stack->Add(h_leading_shower_open_angle_numu_cc_mixed);
	open_angle_stack->Add(h_leading_shower_open_angle_nc);
	open_angle_stack->Add(h_leading_shower_open_angle_nc_pi0);
	open_angle_stack->Add(h_leading_shower_open_angle_other_mixed);
	open_angle_stack->Add(h_leading_shower_open_angle_unmatched);
	open_angle_stack->Draw();
	open_angle_stack->GetXaxis()->SetTitle("Shower Opening Angle [Degrees]");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_leading_shower_open_angle_nue_cc,          "Nue CC", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_cosmic,          "Cosmic", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_numu_cc,         "Numu CC", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_nc,              "NC", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_nc_pi0,          "NC Pi0", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_other_mixed,     "Other Mixed", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_unmatched,       "Unmatched", "f");
	leg_stack->Draw();
	open_angle_stack_c1->Print("post_cuts_leading_shower_open_angle.pdf");

	double open_angle_integral_nue_cc = h_leading_shower_open_angle_nue_cc->Integral();
	double open_angle_integral_nue_cc_mixed = h_leading_shower_open_angle_nue_cc_mixed->Integral();
	double open_angle_integral_cosmic = h_leading_shower_open_angle_cosmic->Integral();
	double open_angle_integral_numu_cc = h_leading_shower_open_angle_numu_cc->Integral();
	double open_angle_integral_numu_cc_mixed = h_leading_shower_open_angle_numu_cc_mixed->Integral();
	double open_angle_integral_nc = h_leading_shower_open_angle_nc->Integral();
	double open_angle_integral_nc_pi0 = h_leading_shower_open_angle_nc_pi0->Integral();
	double open_angle_integral_other_mixed = h_leading_shower_open_angle_other_mixed->Integral();
	double open_angle_integral_unmatched = h_leading_shower_open_angle_unmatched->Integral();
	double open_angle_integral = open_angle_integral_nue_cc
	                             + open_angle_integral_nue_cc_mixed
	                             + open_angle_integral_cosmic
	                             + open_angle_integral_numu_cc
	                             + open_angle_integral_numu_cc_mixed
	                             + open_angle_integral_nc
	                             + open_angle_integral_nc_pi0
	                             + open_angle_integral_other_mixed
	                             + open_angle_integral_unmatched;
	//h_leading_shower_open_angle_data->Scale(1.29916); //- relative scaling (on-beam - off-beam) and POT
	//h_leading_shower_open_angle_data->Scale(5.619); //- POT scaling
	h_leading_shower_open_angle_data->Scale((open_angle_integral / h_leading_shower_open_angle_data->Integral()) * data_mc_scale_factor);
	h_leading_shower_open_angle_data->Draw("P E1 same");
	open_angle_stack_c1->Print("data_mc_post_cuts_leading_shower_open_angle.pdf");

	TCanvas * dedx_cuts_c1 = new TCanvas();
	dedx_cuts_c1->cd();
	THStack * dedx_cuts_stack = new THStack();
	h_dedx_cuts_nue_cc->SetStats(kFALSE);
	h_dedx_cuts_nue_cc_mixed->SetStats(kFALSE);
	h_dedx_cuts_numu_cc->SetStats(kFALSE);
	h_dedx_cuts_nc_pi0->SetStats(kFALSE);
	h_dedx_cuts_cosmic->SetStats(kFALSE);
	h_dedx_cuts_nc->SetStats(kFALSE);
	h_dedx_cuts_numu_cc_mixed->SetStats(kFALSE);
	h_dedx_cuts_other_mixed->SetStats(kFALSE);
	h_dedx_cuts_unmatched->SetStats(kFALSE);
	h_dedx_cuts_nue_cc->SetFillColor(30);
	h_dedx_cuts_nue_cc_mixed->SetFillColor(38);
	h_dedx_cuts_numu_cc->SetFillColor(28);
	h_dedx_cuts_nc_pi0->SetFillColor(36);
	h_dedx_cuts_cosmic->SetFillColor(1);
	h_dedx_cuts_nc->SetFillColor(46);
	h_dedx_cuts_numu_cc_mixed->SetFillColor(25);
	h_dedx_cuts_other_mixed->SetFillColor(42);
	h_dedx_cuts_unmatched->SetFillColor(12);
	dedx_cuts_stack->Add(h_dedx_cuts_nue_cc);
	dedx_cuts_stack->Add(h_dedx_cuts_nue_cc_mixed);
	dedx_cuts_stack->Add(h_dedx_cuts_cosmic);
	dedx_cuts_stack->Add(h_dedx_cuts_numu_cc);
	dedx_cuts_stack->Add(h_dedx_cuts_numu_cc_mixed);
	dedx_cuts_stack->Add(h_dedx_cuts_nc);
	dedx_cuts_stack->Add(h_dedx_cuts_nc_pi0);
	dedx_cuts_stack->Add(h_dedx_cuts_other_mixed);
	dedx_cuts_stack->Add(h_dedx_cuts_unmatched);
	dedx_cuts_stack->Draw();
	dedx_cuts_stack->GetXaxis()->SetTitle("Collection Plane dE/dx [MeV/cm]");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack_dedx = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack_dedx->AddEntry(h_dedx_cuts_nue_cc,          "Nue CC", "f");
	leg_stack_dedx->AddEntry(h_dedx_cuts_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack_dedx->AddEntry(h_dedx_cuts_cosmic,          "Cosmic", "f");
	leg_stack_dedx->AddEntry(h_dedx_cuts_numu_cc,         "Numu CC", "f");
	leg_stack_dedx->AddEntry(h_dedx_cuts_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack_dedx->AddEntry(h_dedx_cuts_nc,              "NC", "f");
	leg_stack_dedx->AddEntry(h_dedx_cuts_nc_pi0,          "NC Pi0", "f");
	leg_stack_dedx->AddEntry(h_dedx_cuts_other_mixed,     "Other Mixed", "f");
	leg_stack_dedx->AddEntry(h_dedx_cuts_unmatched,       "Unmatched", "f");
	leg_stack_dedx->Draw();
	dedx_cuts_c1->Print("post_cuts_dedx_cuts.pdf");

	double dedx_integral_nue_cc         = h_dedx_cuts_nue_cc->Integral();
	double dedx_integral_nue_cc_mixed   = h_dedx_cuts_nue_cc_mixed->Integral();
	double dedx_integral_cosmic         = h_dedx_cuts_cosmic->Integral();
	double dedx_integral_numu_cc        = h_dedx_cuts_numu_cc->Integral();
	double dedx_integral_numu_cc_mixed  = h_dedx_cuts_numu_cc_mixed->Integral();
	double dedx_integral_nc             = h_dedx_cuts_nc->Integral();
	double dedx_integral_nc_pi0         = h_dedx_cuts_nc_pi0->Integral();
	double dedx_integral_other_mixed    = h_dedx_cuts_other_mixed->Integral();
	double dedx_integral_unmatched      = h_dedx_cuts_unmatched->Integral();
	double dedx_integral = dedx_integral_nue_cc
	                       + dedx_integral_nue_cc_mixed
	                       + dedx_integral_cosmic
	                       + dedx_integral_numu_cc
	                       + dedx_integral_numu_cc_mixed
	                       + dedx_integral_nc
	                       + dedx_integral_nc_pi0
	                       + dedx_integral_other_mixed
	                       + dedx_integral_unmatched;

	//h_dedx_cuts_data->Scale(1.2991626);
	//h_dedx_cuts_data->Scale(5.619);
	h_dedx_cuts_data->Scale((dedx_integral / h_dedx_cuts_data->Integral()) * data_mc_scale_factor);
	h_dedx_cuts_data->Draw("P E1 same");
	dedx_cuts_c1->Print("data_mc_post_cuts_dedx_cuts.pdf");

	TCanvas * vtx_to_flash_c1 = new TCanvas();
	vtx_to_flash_c1->cd();
	THStack * vtx_to_flash_stack = new THStack();
	h_vtx_flash_nue_cc->SetStats(kFALSE);
	h_vtx_flash_nue_cc_mixed->SetStats(kFALSE);
	h_vtx_flash_numu_cc->SetStats(kFALSE);
	h_vtx_flash_nc_pi0->SetStats(kFALSE);
	h_vtx_flash_cosmic->SetStats(kFALSE);
	h_vtx_flash_nc->SetStats(kFALSE);
	h_vtx_flash_numu_cc_mixed->SetStats(kFALSE);
	h_vtx_flash_other_mixed->SetStats(kFALSE);
	h_vtx_flash_unmatched->SetStats(kFALSE);
	h_vtx_flash_nue_cc->SetFillColor(30);
	h_vtx_flash_nue_cc_mixed->SetFillColor(38);
	h_vtx_flash_numu_cc->SetFillColor(28);
	h_vtx_flash_nc_pi0->SetFillColor(36);
	h_vtx_flash_cosmic->SetFillColor(1);
	h_vtx_flash_nc->SetFillColor(46);
	h_vtx_flash_numu_cc_mixed->SetFillColor(25);
	h_vtx_flash_other_mixed->SetFillColor(42);
	h_vtx_flash_unmatched->SetFillColor(12);
	vtx_to_flash_stack->Add(h_vtx_flash_nue_cc);
	vtx_to_flash_stack->Add(h_vtx_flash_nue_cc_mixed);
	vtx_to_flash_stack->Add(h_vtx_flash_cosmic);
	vtx_to_flash_stack->Add(h_vtx_flash_numu_cc);
	vtx_to_flash_stack->Add(h_vtx_flash_numu_cc_mixed);
	vtx_to_flash_stack->Add(h_vtx_flash_nc);
	vtx_to_flash_stack->Add(h_vtx_flash_nc_pi0);
	vtx_to_flash_stack->Add(h_vtx_flash_other_mixed);
	vtx_to_flash_stack->Add(h_vtx_flash_unmatched);
	vtx_to_flash_stack->Draw();
	vtx_to_flash_stack->GetXaxis()->SetTitle("2D Distance From Largest Flash to Reco Nu Vtx [cm]");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack_flash = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack_flash->AddEntry(h_vtx_flash_nue_cc,          "Nue CC", "f");
	leg_stack_flash->AddEntry(h_vtx_flash_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack_flash->AddEntry(h_vtx_flash_cosmic,          "Cosmic", "f");
	leg_stack_flash->AddEntry(h_vtx_flash_numu_cc,         "Numu CC", "f");
	leg_stack_flash->AddEntry(h_vtx_flash_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack_flash->AddEntry(h_vtx_flash_nc,              "NC", "f");
	leg_stack_flash->AddEntry(h_vtx_flash_nc_pi0,          "NC Pi0", "f");
	leg_stack_flash->AddEntry(h_vtx_flash_other_mixed,     "Other Mixed", "f");
	leg_stack_flash->AddEntry(h_vtx_flash_unmatched,       "Unmatched", "f");
	leg_stack_flash->Draw();
	vtx_to_flash_c1->Print("post_cuts_vtx_to_flash_distance.pdf");

	double vtx_flash_integral_nue_cc         = h_vtx_flash_nue_cc->Integral();
	double vtx_flash_integral_nue_cc_mixed   = h_vtx_flash_nue_cc_mixed->Integral();
	double vtx_flash_integral_cosmic         = h_vtx_flash_cosmic->Integral();
	double vtx_flash_integral_numu_cc        = h_vtx_flash_numu_cc->Integral();
	double vtx_flash_integral_numu_cc_mixed  = h_vtx_flash_numu_cc_mixed->Integral();
	double vtx_flash_integral_nc             = h_vtx_flash_nc->Integral();
	double vtx_flash_integral_nc_pi0         = h_vtx_flash_nc_pi0->Integral();
	double vtx_flash_integral_other_mixed    = h_vtx_flash_other_mixed->Integral();
	double vtx_flash_integral_unmatched      = h_vtx_flash_unmatched->Integral();
	double vtx_flash_integral = vtx_flash_integral_nue_cc
	                            + vtx_flash_integral_nue_cc_mixed
	                            + vtx_flash_integral_cosmic
	                            + vtx_flash_integral_numu_cc
	                            + vtx_flash_integral_numu_cc_mixed
	                            + vtx_flash_integral_nc
	                            + vtx_flash_integral_nc_pi0
	                            + vtx_flash_integral_other_mixed
	                            + vtx_flash_integral_unmatched;

	//h_vtx_flash_data->Scale(1.2991626);
	//h_vtx_flash_data->Scale(5.619);
	h_vtx_flash_data->Scale((vtx_flash_integral / h_vtx_flash_data->Integral()) * data_mc_scale_factor);
	h_vtx_flash_data->Draw("P E1 same");
	vtx_to_flash_c1->Print("data_mc_vtx_to_flash_distance.pdf");

	TCanvas * trk_vtx_dist_stack_c1 = new TCanvas();
	trk_vtx_dist_stack_c1->cd();
	THStack * trk_vtx_dist_stack = new THStack();
	h_trk_vtx_dist_nue_cc->SetStats(kFALSE);
	h_trk_vtx_dist_nue_cc_mixed->SetStats(kFALSE);
	h_trk_vtx_dist_numu_cc->SetStats(kFALSE);
	h_trk_vtx_dist_nc_pi0->SetStats(kFALSE);
	h_trk_vtx_dist_cosmic->SetStats(kFALSE);
	h_trk_vtx_dist_nc->SetStats(kFALSE);
	h_trk_vtx_dist_numu_cc_mixed->SetStats(kFALSE);
	h_trk_vtx_dist_other_mixed->SetStats(kFALSE);
	h_trk_vtx_dist_unmatched->SetStats(kFALSE);
	h_trk_vtx_dist_nue_cc->SetFillColor(30);
	h_trk_vtx_dist_nue_cc_mixed->SetFillColor(38);
	h_trk_vtx_dist_numu_cc->SetFillColor(28);
	h_trk_vtx_dist_nc_pi0->SetFillColor(36);
	h_trk_vtx_dist_cosmic->SetFillColor(1);
	h_trk_vtx_dist_nc->SetFillColor(46);
	h_trk_vtx_dist_numu_cc_mixed->SetFillColor(25);
	h_trk_vtx_dist_other_mixed->SetFillColor(42);
	h_trk_vtx_dist_unmatched->SetFillColor(12);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_nue_cc);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_nue_cc_mixed);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_cosmic);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_numu_cc);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_numu_cc_mixed);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_nc);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_nc_pi0);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_other_mixed);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_unmatched);
	trk_vtx_dist_stack->Draw();
	trk_vtx_dist_stack->GetXaxis()->SetTitle("Track to Nue Candidate Vertex Distance [cm]");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack2 = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack2->AddEntry(h_trk_vtx_dist_nue_cc,          "Nue CC", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_cosmic,          "Cosmic", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_numu_cc,         "Numu CC", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_nc,              "NC", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_nc_pi0,          "NC Pi0", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_other_mixed,     "Other Mixed", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_unmatched,       "Unmatched", "f");
	leg_stack2->Draw();
	trk_vtx_dist_stack_c1->Print("post_cuts_track_to_vtx.pdf");

	double trk_vtx_integral_nue_cc          = h_trk_vtx_dist_nue_cc->Integral();
	double trk_vtx_integral_nue_cc_mixed    = h_trk_vtx_dist_nue_cc_mixed->Integral();
	double trk_vtx_integral_cosmic          = h_trk_vtx_dist_cosmic->Integral();
	double trk_vtx_integral_numu_cc         = h_trk_vtx_dist_numu_cc->Integral();
	double trk_vtx_integral_numu_cc_mixed   = h_trk_vtx_dist_numu_cc_mixed->Integral();
	double trk_vtx_integral_nc              = h_trk_vtx_dist_nc->Integral();
	double trk_vtx_integral_nc_pi0          = h_trk_vtx_dist_nc_pi0->Integral();
	double trk_vtx_integral_other_mixed     = h_trk_vtx_dist_other_mixed->Integral();
	double trk_vtx_integral_unmatched       = h_trk_vtx_dist_unmatched->Integral();
	double trk_vtx_integral = trk_vtx_integral_nue_cc
	                          + trk_vtx_integral_nue_cc_mixed
	                          + trk_vtx_integral_cosmic
	                          + trk_vtx_integral_numu_cc
	                          + trk_vtx_integral_numu_cc_mixed
	                          + trk_vtx_integral_nc
	                          + trk_vtx_integral_nc_pi0
	                          + trk_vtx_integral_other_mixed
	                          + trk_vtx_integral_unmatched;

	//h_trk_vtx_dist_data->Scale(1.2991626);
	//h_trk_vtx_dist_data->Scale(5.619);
	h_trk_vtx_dist_data->Scale((trk_vtx_integral / h_trk_vtx_dist_data->Integral()) * data_mc_scale_factor);
	h_trk_vtx_dist_data->Draw("P E1 same");
	trk_vtx_dist_stack_c1->Print("data_mc_post_cuts_track_to_vtx.pdf");

	TCanvas * shwr_vtx_dist_stack_c1 = new TCanvas();
	shwr_vtx_dist_stack_c1->cd();
	THStack * shwr_vtx_dist_stack = new THStack();
	h_shwr_vtx_dist_nue_cc->SetStats(kFALSE);
	h_shwr_vtx_dist_nue_cc_mixed->SetStats(kFALSE);
	h_shwr_vtx_dist_numu_cc->SetStats(kFALSE);
	h_shwr_vtx_dist_nc_pi0->SetStats(kFALSE);
	h_shwr_vtx_dist_cosmic->SetStats(kFALSE);
	h_shwr_vtx_dist_nc->SetStats(kFALSE);
	h_shwr_vtx_dist_numu_cc_mixed->SetStats(kFALSE);
	h_shwr_vtx_dist_other_mixed->SetStats(kFALSE);
	h_shwr_vtx_dist_unmatched->SetStats(kFALSE);
	h_shwr_vtx_dist_nue_cc->SetFillColor(30);
	h_shwr_vtx_dist_nue_cc_mixed->SetFillColor(38);
	h_shwr_vtx_dist_numu_cc->SetFillColor(28);
	h_shwr_vtx_dist_nc_pi0->SetFillColor(36);
	h_shwr_vtx_dist_cosmic->SetFillColor(1);
	h_shwr_vtx_dist_nc->SetFillColor(46);
	h_shwr_vtx_dist_numu_cc_mixed->SetFillColor(25);
	h_shwr_vtx_dist_other_mixed->SetFillColor(42);
	h_shwr_vtx_dist_unmatched->SetFillColor(12);
	shwr_vtx_dist_stack->Add(h_shwr_vtx_dist_nue_cc);
	shwr_vtx_dist_stack->Add(h_shwr_vtx_dist_nue_cc_mixed);
	shwr_vtx_dist_stack->Add(h_shwr_vtx_dist_cosmic);
	shwr_vtx_dist_stack->Add(h_shwr_vtx_dist_numu_cc);
	shwr_vtx_dist_stack->Add(h_shwr_vtx_dist_numu_cc_mixed);
	shwr_vtx_dist_stack->Add(h_shwr_vtx_dist_nc);
	shwr_vtx_dist_stack->Add(h_shwr_vtx_dist_nc_pi0);
	shwr_vtx_dist_stack->Add(h_shwr_vtx_dist_other_mixed);
	shwr_vtx_dist_stack->Add(h_shwr_vtx_dist_unmatched);
	shwr_vtx_dist_stack->Draw();
	shwr_vtx_dist_stack->GetXaxis()->SetTitle("Shower to Nue Candidate Vertex Distance [cm]");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack_shwr = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack_shwr->AddEntry(h_shwr_vtx_dist_nue_cc,          "Nue CC", "f");
	leg_stack_shwr->AddEntry(h_shwr_vtx_dist_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack_shwr->AddEntry(h_shwr_vtx_dist_cosmic,          "Cosmic", "f");
	leg_stack_shwr->AddEntry(h_shwr_vtx_dist_numu_cc,         "Numu CC", "f");
	leg_stack_shwr->AddEntry(h_shwr_vtx_dist_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack_shwr->AddEntry(h_shwr_vtx_dist_nc,              "NC", "f");
	leg_stack_shwr->AddEntry(h_shwr_vtx_dist_nc_pi0,          "NC Pi0", "f");
	leg_stack_shwr->AddEntry(h_shwr_vtx_dist_other_mixed,     "Other Mixed", "f");
	leg_stack_shwr->AddEntry(h_shwr_vtx_dist_unmatched,       "Unmatched", "f");
	leg_stack_shwr->Draw();
	shwr_vtx_dist_stack_c1->Print("post_cuts_shower_to_vtx.pdf");

	double shwr_vtx_integral_nue_cc         = h_shwr_vtx_dist_nue_cc->Integral();
	double shwr_vtx_integral_nue_cc_mixed   = h_shwr_vtx_dist_nue_cc_mixed->Integral();
	double shwr_vtx_integral_cosmic         = h_shwr_vtx_dist_cosmic->Integral();
	double shwr_vtx_integral_numu_cc        = h_shwr_vtx_dist_numu_cc->Integral();
	double shwr_vtx_integral_numu_cc_mixed  = h_shwr_vtx_dist_numu_cc_mixed->Integral();
	double shwr_vtx_integral_nc             = h_shwr_vtx_dist_nc->Integral();
	double shwr_vtx_integral_nc_pi0         = h_shwr_vtx_dist_nc_pi0->Integral();
	double shwr_vtx_integral_other_mixed    = h_shwr_vtx_dist_other_mixed->Integral();
	double shwr_vtx_integral_unmatched      = h_shwr_vtx_dist_unmatched->Integral();
	double shwr_vtx_integral = shwr_vtx_integral_nue_cc
	                           + shwr_vtx_integral_nue_cc_mixed
	                           + shwr_vtx_integral_cosmic
	                           + shwr_vtx_integral_numu_cc
	                           + shwr_vtx_integral_numu_cc_mixed
	                           + shwr_vtx_integral_nc
	                           + shwr_vtx_integral_nc_pi0
	                           + shwr_vtx_integral_other_mixed
	                           + shwr_vtx_integral_unmatched;

	//h_shwr_vtx_dist_data->Scale(1.2991626);
	//h_shwr_vtx_dist_data->Scale(5.619);
	h_shwr_vtx_dist_data->Scale((shwr_vtx_integral / h_shwr_vtx_dist_data->Integral()) * data_mc_scale_factor);
	h_shwr_vtx_dist_data->Draw("P E1 same");
	shwr_vtx_dist_stack_c1->Print("data_mc_post_cuts_shower_to_vtx.pdf");

	TCanvas * track_stack_c1 = new TCanvas();
	track_stack_c1->cd();
	THStack * h_track_stack = new THStack();
	h_pfp_track_nue_cc_qe->SetStats(kFALSE);
	h_pfp_track_nue_cc_out_fv->SetStats(kFALSE);
	h_pfp_track_nue_cc_res->SetStats(kFALSE);
	h_pfp_track_nue_cc_dis->SetStats(kFALSE);
	h_pfp_track_nue_cc_coh->SetStats(kFALSE);
	h_pfp_track_nue_cc_mec->SetStats(kFALSE);
	h_pfp_track_nc->SetStats(kFALSE);
	h_pfp_track_numu_cc_qe->SetStats(kFALSE);
	h_pfp_track_numu_cc_res->SetStats(kFALSE);
	h_pfp_track_numu_cc_dis->SetStats(kFALSE);
	h_pfp_track_numu_cc_coh->SetStats(kFALSE);
	h_pfp_track_numu_cc_mec->SetStats(kFALSE);
	h_pfp_track_nc_pi0->SetStats(kFALSE);
	h_pfp_track_nue_cc_mixed->SetStats(kFALSE);
	h_pfp_track_numu_cc_mixed->SetStats(kFALSE);
	h_pfp_track_cosmic->SetStats(kFALSE);
	h_pfp_track_other_mixed->SetStats(kFALSE);
	h_pfp_track_unmatched->SetStats(kFALSE);
	h_pfp_track_nue_cc_qe->SetFillColor(30);
	h_pfp_track_nue_cc_out_fv->SetFillColor(45);
	h_pfp_track_nue_cc_res->SetFillColor(31);
	h_pfp_track_nue_cc_dis->SetFillColor(32);
	h_pfp_track_nue_cc_coh->SetFillColor(33);
	h_pfp_track_nue_cc_mec->SetFillColor(34);
	h_pfp_track_nc->SetFillColor(46);
	h_pfp_track_numu_cc_qe->SetFillColor(28);
	h_pfp_track_numu_cc_res->SetFillColor(27);
	h_pfp_track_numu_cc_dis->SetFillColor(26);
	h_pfp_track_numu_cc_coh->SetFillColor(23);
	h_pfp_track_numu_cc_mec->SetFillColor(22);
	h_pfp_track_nc_pi0->SetFillColor(36);
	h_pfp_track_nue_cc_mixed->SetFillColor(38);
	h_pfp_track_numu_cc_mixed->SetFillColor(25);
	h_pfp_track_cosmic->SetFillColor(39);
	h_pfp_track_other_mixed->SetFillColor(42);
	h_pfp_track_unmatched->SetFillColor(12);
	h_track_stack->Add(h_pfp_track_nue_cc_qe    );
	h_track_stack->Add(h_pfp_track_nue_cc_out_fv);
	h_track_stack->Add(h_pfp_track_nue_cc_res   );
	h_track_stack->Add(h_pfp_track_nue_cc_dis   );
	h_track_stack->Add(h_pfp_track_nue_cc_coh   );
	h_track_stack->Add(h_pfp_track_nue_cc_mec   );
	h_track_stack->Add(h_pfp_track_nc       );
	h_track_stack->Add(h_pfp_track_numu_cc_qe   );
	h_track_stack->Add(h_pfp_track_numu_cc_res  );
	h_track_stack->Add(h_pfp_track_numu_cc_dis  );
	h_track_stack->Add(h_pfp_track_numu_cc_coh  );
	h_track_stack->Add(h_pfp_track_numu_cc_mec  );
	h_track_stack->Add(h_pfp_track_nc_pi0      );
	h_track_stack->Add(h_pfp_track_nue_cc_mixed );
	h_track_stack->Add(h_pfp_track_numu_cc_mixed);
	h_track_stack->Add(h_pfp_track_cosmic       );
	h_track_stack->Add(h_pfp_track_other_mixed  );
	h_track_stack->Add(h_pfp_track_unmatched    );
	h_track_stack->Draw();
	h_track_stack->GetXaxis()->SetTitle("Reconstructed Tracks in Candidate Neutrino Object");
	TLegend * leg_track_stack_l1 = new TLegend(0.75, 0.50, 0.95, 0.95);
	leg_track_stack_l1->AddEntry(h_pfp_track_nue_cc_qe,      "Nue CC QE", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_nue_cc_out_fv,  "Nue CC Out FV", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_nue_cc_res,     "Nue CC Res", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_nue_cc_dis,     "Nue CC DIS", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_nue_cc_coh,     "Nue CC Coh", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_nue_cc_mec,     "Nue CC MEC", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_nc,         "NC", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_qe,     "Numu CC QE", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_res,    "Numu CC Res", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_dis,    "Numu CC DIS", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_coh,    "Numu CC Coh", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_mec,    "Numu CC MEC", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_nc_pi0,        "NC Pi0", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_nue_cc_mixed,   "Nue CC Mixed", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_mixed,  "Numu CC Mixed", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_cosmic,         "Cosmic", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_other_mixed,    "Other Mixed", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_unmatched,      "Unmatched", "f");
	leg_track_stack_l1->Draw();
	track_stack_c1->Print("selected_pfp_track_stack.pdf");

	//h_pfp_track_data->Scale(1.2991626);
	h_pfp_track_data->Scale(5.619);
	h_pfp_track_data->Draw("P E1 same");
	track_stack_c1->Print("data_mc_selected_pfp_track_stack.pdf");

	TCanvas * shower_stack_c1 = new TCanvas();
	shower_stack_c1->cd();
	THStack * h_shower_stack = new THStack();
	h_pfp_shower_nue_cc_qe->SetStats(kFALSE);
	h_pfp_shower_nue_cc_out_fv->SetStats(kFALSE);
	h_pfp_shower_nue_cc_res->SetStats(kFALSE);
	h_pfp_shower_nue_cc_dis->SetStats(kFALSE);
	h_pfp_shower_nue_cc_coh->SetStats(kFALSE);
	h_pfp_shower_nue_cc_mec->SetStats(kFALSE);
	h_pfp_shower_nc->SetStats(kFALSE);
	h_pfp_shower_numu_cc_qe->SetStats(kFALSE);
	h_pfp_shower_numu_cc_res->SetStats(kFALSE);
	h_pfp_shower_numu_cc_dis->SetStats(kFALSE);
	h_pfp_shower_numu_cc_coh->SetStats(kFALSE);
	h_pfp_shower_numu_cc_mec->SetStats(kFALSE);
	h_pfp_shower_nc_pi0->SetStats(kFALSE);
	h_pfp_shower_nue_cc_mixed->SetStats(kFALSE);
	h_pfp_shower_numu_cc_mixed->SetStats(kFALSE);
	h_pfp_shower_cosmic->SetStats(kFALSE);
	h_pfp_shower_other_mixed->SetStats(kFALSE);
	h_pfp_shower_unmatched->SetStats(kFALSE);
	h_pfp_shower_nue_cc_qe->SetFillColor(30);
	h_pfp_shower_nue_cc_out_fv->SetFillColor(45);
	h_pfp_shower_nue_cc_res->SetFillColor(31);
	h_pfp_shower_nue_cc_dis->SetFillColor(32);
	h_pfp_shower_nue_cc_coh->SetFillColor(33);
	h_pfp_shower_nue_cc_mec->SetFillColor(34);
	h_pfp_shower_nc->SetFillColor(46);
	h_pfp_shower_numu_cc_qe->SetFillColor(28);
	h_pfp_shower_numu_cc_res->SetFillColor(27);
	h_pfp_shower_numu_cc_dis->SetFillColor(26);
	h_pfp_shower_numu_cc_coh->SetFillColor(23);
	h_pfp_shower_numu_cc_mec->SetFillColor(22);
	h_pfp_shower_nc_pi0->SetFillColor(36);
	h_pfp_shower_nue_cc_mixed->SetFillColor(38);
	h_pfp_shower_numu_cc_mixed->SetFillColor(25);
	h_pfp_shower_cosmic->SetFillColor(39);
	h_pfp_shower_other_mixed->SetFillColor(42);
	h_pfp_shower_unmatched->SetFillColor(12);
	h_shower_stack->Add(h_pfp_shower_nue_cc_qe    );
	h_shower_stack->Add(h_pfp_shower_nue_cc_out_fv);
	h_shower_stack->Add(h_pfp_shower_nue_cc_res   );
	h_shower_stack->Add(h_pfp_shower_nue_cc_dis   );
	h_shower_stack->Add(h_pfp_shower_nue_cc_coh   );
	h_shower_stack->Add(h_pfp_shower_nue_cc_mec   );
	h_shower_stack->Add(h_pfp_shower_nc       );
	h_shower_stack->Add(h_pfp_shower_numu_cc_qe   );
	h_shower_stack->Add(h_pfp_shower_numu_cc_res  );
	h_shower_stack->Add(h_pfp_shower_numu_cc_dis  );
	h_shower_stack->Add(h_pfp_shower_numu_cc_coh  );
	h_shower_stack->Add(h_pfp_shower_numu_cc_mec  );
	h_shower_stack->Add(h_pfp_shower_nc_pi0      );
	h_shower_stack->Add(h_pfp_shower_nue_cc_mixed );
	h_shower_stack->Add(h_pfp_shower_numu_cc_mixed);
	h_shower_stack->Add(h_pfp_shower_cosmic       );
	h_shower_stack->Add(h_pfp_shower_other_mixed  );
	h_shower_stack->Add(h_pfp_shower_unmatched    );
	h_shower_stack->Draw();
	h_shower_stack->GetXaxis()->SetTitle("Reconstructed Showers in Candidate Neutrino Object");
	TLegend * leg_shower_stack_l1 = new TLegend(0.75, 0.50, 0.95, 0.95);
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nue_cc_qe,      "Nue CC QE", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nue_cc_out_fv,  "Nue CC Out FV", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nue_cc_res,     "Nue CC Res", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nue_cc_dis,     "Nue CC DIS", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nue_cc_coh,     "Nue CC Coh", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nue_cc_mec,     "Nue CC MEC", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nc,         "NC", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_qe,     "Numu CC QE", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_res,    "Numu CC Res", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_dis,    "Numu CC DIS", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_coh,    "Numu CC Coh", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_mec,    "Numu CC MEC", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nc_pi0,        "NC Pi0", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nue_cc_mixed,   "Nue CC Mixed", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_mixed,  "Numu CC Mixed", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_cosmic,         "Cosmic", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_other_mixed,    "Other Mixed", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_unmatched,      "Unmatched", "f");
	leg_shower_stack_l1->Draw();
	shower_stack_c1->Print("selected_pfp_shower_stack.pdf");

	//h_pfp_shower_data->Scale(1.2991626);
	h_pfp_shower_data->Scale(5.619);
	h_pfp_shower_data->Draw("P E1 same");
	shower_stack_c1->Print("data_mc_selected_pfp_shower_stack.pdf");

	TCanvas * track_shower_c1 = new TCanvas();
	track_shower_c1->cd();
	h_pfp_track_shower_nue_cc_qe->Draw("colz");
	h_pfp_track_shower_nue_cc_qe->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_qe->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c1->Print("selected_pfp_track_shower_nue_cc_qe.pdf");
	TCanvas * track_shower_c2 = new TCanvas();
	track_shower_c2->cd();
	h_pfp_track_shower_nue_cc_out_fv->Draw("colz");
	h_pfp_track_shower_nue_cc_out_fv->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_out_fv->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c2->Print("selected_pfp_track_shower_nue_cc_out_fv.pdf");
	TCanvas * track_shower_c3 = new TCanvas();
	track_shower_c3->cd();
	h_pfp_track_shower_nue_cc_res->Draw("colz");
	h_pfp_track_shower_nue_cc_res->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_res->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c3->Print("selected_pfp_track_shower_nue_cc_res.pdf");
	TCanvas * track_shower_c4 = new TCanvas();
	track_shower_c4->cd();
	h_pfp_track_shower_nue_cc_dis->Draw("colz");
	h_pfp_track_shower_nue_cc_dis->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_dis->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c4->Print("selected_pfp_track_shower_nue_cc_dis.pdf");
	TCanvas * track_shower_c5 = new TCanvas();
	track_shower_c5->cd();
	h_pfp_track_shower_nue_cc_coh->Draw("colz");
	h_pfp_track_shower_nue_cc_coh->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_coh->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c5->Print("selected_pfp_track_shower_nue_cc_coh.pdf");
	TCanvas * track_shower_c6 = new TCanvas();
	track_shower_c6->cd();
	h_pfp_track_shower_nue_cc_mec->Draw("colz");
	h_pfp_track_shower_nue_cc_mec->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_mec->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c6->Print("selected_pfp_track_shower_nue_cc_mec.pdf");
	TCanvas * track_shower_c7 = new TCanvas();
	track_shower_c7->cd();
	h_pfp_track_shower_nc->Draw("colz");
	h_pfp_track_shower_nc->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nc->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c7->Print("selected_pfp_track_shower_nc.pdf");
	TCanvas * track_shower_c8 = new TCanvas();
	track_shower_c8->cd();
	h_pfp_track_shower_numu_cc_qe->Draw("colz");
	h_pfp_track_shower_numu_cc_qe->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_qe->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c8->Print("selected_pfp_track_shower_numu_cc_qe.pdf");
	TCanvas * track_shower_c9 = new TCanvas();
	track_shower_c9->cd();
	h_pfp_track_shower_numu_cc_res->Draw("colz");
	h_pfp_track_shower_numu_cc_res->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_res->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c9->Print("selected_pfp_track_shower_numu_cc_res.pdf");
	TCanvas * track_shower_c10 = new TCanvas();
	track_shower_c10->cd();
	h_pfp_track_shower_numu_cc_dis->Draw("colz");
	h_pfp_track_shower_numu_cc_dis->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_dis->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c10->Print("selected_pfp_track_shower_numu_cc_dis.pdf");
	TCanvas * track_shower_c11 = new TCanvas();
	track_shower_c11->cd();
	h_pfp_track_shower_numu_cc_coh->Draw("colz");
	h_pfp_track_shower_numu_cc_coh->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_coh->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c11->Print("selected_pfp_track_shower_numu_cc_coh.pdf");
	TCanvas * track_shower_c12 = new TCanvas();
	track_shower_c12->cd();
	h_pfp_track_shower_numu_cc_mec->Draw("colz");
	h_pfp_track_shower_numu_cc_mec->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_mec->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c12->Print("selected_pfp_track_shower_numu_cc_mec.pdf");
	TCanvas * track_shower_c13 = new TCanvas();
	track_shower_c13->cd();
	h_pfp_track_shower_nc_pi0->Draw("colz");
	h_pfp_track_shower_nc_pi0->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nc_pi0->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c13->Print("selected_pfp_track_shower_nc_pi0.pdf");
	TCanvas * track_shower_c14 = new TCanvas();
	track_shower_c14->cd();
	h_pfp_track_shower_nue_cc_mixed->Draw("colz");
	h_pfp_track_shower_nue_cc_mixed->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_mixed->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c14->Print("selected_pfp_track_shower_nue_cc_mixed.pdf");
	TCanvas * track_shower_c15 = new TCanvas();
	track_shower_c15->cd();
	h_pfp_track_shower_numu_cc_mixed->Draw("colz");
	h_pfp_track_shower_numu_cc_mixed->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_mixed->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c15->Print("selected_pfp_track_shower_numu_cc_mixed.pdf");
	TCanvas * track_shower_c16 = new TCanvas();
	track_shower_c16->cd();
	h_pfp_track_shower_cosmic->Draw("colz");
	h_pfp_track_shower_cosmic->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_cosmic->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c16->Print("selected_pfp_track_shower_cosmic.pdf");
	TCanvas * track_shower_c17 = new TCanvas();
	track_shower_c17->cd();
	h_pfp_track_shower_other_mixed->Draw("colz");
	h_pfp_track_shower_other_mixed->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_other_mixed->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c17->Print("selected_pfp_track_shower_other_mixed.pdf");
	TCanvas * track_shower_c18 = new TCanvas();
	track_shower_c18->cd();
	h_pfp_track_shower_unmatched->Draw("colz");
	h_pfp_track_shower_unmatched->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_unmatched->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c18->Print("selected_pfp_track_shower_unmatched.pdf");
	TCanvas * track_shower_c19 = new TCanvas();
	track_shower_c19->cd();
	h_pfp_track_shower_data->Draw("colz");
	h_pfp_track_shower_data->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_data->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c19->Print("selected_pfp_track_shower_data.pdf");

	TCanvas * track_stack_c2 = new TCanvas();
	track_stack_c2->cd();
	THStack * h_track_stack_last = new THStack();
	h_pfp_track_nue_cc_qe_last->SetStats(kFALSE);
	h_pfp_track_nue_cc_out_fv_last->SetStats(kFALSE);
	h_pfp_track_nue_cc_res_last->SetStats(kFALSE);
	h_pfp_track_nue_cc_dis_last->SetStats(kFALSE);
	h_pfp_track_nue_cc_coh_last->SetStats(kFALSE);
	h_pfp_track_nue_cc_mec_last->SetStats(kFALSE);
	h_pfp_track_nc_last->SetStats(kFALSE);
	h_pfp_track_numu_cc_qe_last->SetStats(kFALSE);
	h_pfp_track_numu_cc_res_last->SetStats(kFALSE);
	h_pfp_track_numu_cc_dis_last->SetStats(kFALSE);
	h_pfp_track_numu_cc_coh_last->SetStats(kFALSE);
	h_pfp_track_numu_cc_mec_last->SetStats(kFALSE);
	h_pfp_track_nc_pi0_last->SetStats(kFALSE);
	h_pfp_track_nue_cc_mixed_last->SetStats(kFALSE);
	h_pfp_track_numu_cc_mixed_last->SetStats(kFALSE);
	h_pfp_track_cosmic_last->SetStats(kFALSE);
	h_pfp_track_other_mixed_last->SetStats(kFALSE);
	h_pfp_track_unmatched_last->SetStats(kFALSE);
	h_pfp_track_nue_cc_qe_last->SetFillColor(30);
	h_pfp_track_nue_cc_out_fv_last->SetFillColor(45);
	h_pfp_track_nue_cc_res_last->SetFillColor(31);
	h_pfp_track_nue_cc_dis_last->SetFillColor(32);
	h_pfp_track_nue_cc_coh_last->SetFillColor(33);
	h_pfp_track_nue_cc_mec_last->SetFillColor(34);
	h_pfp_track_nc_last->SetFillColor(46);
	h_pfp_track_numu_cc_qe_last->SetFillColor(28);
	h_pfp_track_numu_cc_res_last->SetFillColor(27);
	h_pfp_track_numu_cc_dis_last->SetFillColor(26);
	h_pfp_track_numu_cc_coh_last->SetFillColor(23);
	h_pfp_track_numu_cc_mec_last->SetFillColor(22);
	h_pfp_track_nc_pi0_last->SetFillColor(36);
	h_pfp_track_nue_cc_mixed_last->SetFillColor(38);
	h_pfp_track_numu_cc_mixed_last->SetFillColor(25);
	h_pfp_track_cosmic_last->SetFillColor(39);
	h_pfp_track_other_mixed_last->SetFillColor(42);
	h_pfp_track_unmatched_last->SetFillColor(12);
	h_track_stack_last->Add(h_pfp_track_nue_cc_qe_last    );
	h_track_stack_last->Add(h_pfp_track_nue_cc_out_fv_last);
	h_track_stack_last->Add(h_pfp_track_nue_cc_res_last   );
	h_track_stack_last->Add(h_pfp_track_nue_cc_dis_last   );
	h_track_stack_last->Add(h_pfp_track_nue_cc_coh_last   );
	h_track_stack_last->Add(h_pfp_track_nue_cc_mec_last   );
	h_track_stack_last->Add(h_pfp_track_nc_last       );
	h_track_stack_last->Add(h_pfp_track_numu_cc_qe_last   );
	h_track_stack_last->Add(h_pfp_track_numu_cc_res_last  );
	h_track_stack_last->Add(h_pfp_track_numu_cc_dis_last  );
	h_track_stack_last->Add(h_pfp_track_numu_cc_coh_last  );
	h_track_stack_last->Add(h_pfp_track_numu_cc_mec_last  );
	h_track_stack_last->Add(h_pfp_track_nc_pi0_last      );
	h_track_stack_last->Add(h_pfp_track_nue_cc_mixed_last );
	h_track_stack_last->Add(h_pfp_track_numu_cc_mixed_last);
	h_track_stack_last->Add(h_pfp_track_cosmic_last       );
	h_track_stack_last->Add(h_pfp_track_other_mixed_last  );
	h_track_stack_last->Add(h_pfp_track_unmatched_last    );
	h_track_stack_last->Draw();
	h_track_stack_last->GetXaxis()->SetTitle("Reconstructed Tracks in Candidate Neutrino Object");
	TLegend * leg_track_stack_l2 = new TLegend(0.75, 0.50, 0.95, 0.95);
	leg_track_stack_l2->AddEntry(h_pfp_track_nue_cc_qe_last,      "Nue CC QE", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_nue_cc_out_fv_last,  "Nue CC Out FV", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_nue_cc_res_last,     "Nue CC Res", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_nue_cc_dis_last,     "Nue CC DIS", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_nue_cc_coh_last,     "Nue CC Coh", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_nue_cc_mec_last,     "Nue CC MEC", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_nc_last,         "NC", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_numu_cc_qe_last,     "Numu CC QE", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_numu_cc_res_last,    "Numu CC Res", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_numu_cc_dis_last,    "Numu CC DIS", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_numu_cc_coh_last,    "Numu CC Coh", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_numu_cc_mec_last,    "Numu CC MEC", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_nc_pi0_last,        "NC Pi0", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_nue_cc_mixed_last,   "Nue CC Mixed", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_numu_cc_mixed_last,  "Numu CC Mixed", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_cosmic_last,         "Cosmic", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_other_mixed_last,    "Other Mixed", "f");
	leg_track_stack_l2->AddEntry(h_pfp_track_unmatched_last,      "Unmatched", "f");
	leg_track_stack_l2->Draw();
	track_stack_c2->Print("selected_pfp_track_stack_last.pdf");

	//h_pfp_track_data_last->Scale(1.2991626);
	h_pfp_track_data_last->Scale(5.619);
	h_pfp_track_data_last->Draw("P E1 same");
	track_stack_c2->Print("data_mc_selected_pfp_track_stack_last.pdf");

	TCanvas * shower_stack_c2 = new TCanvas();
	shower_stack_c2->cd();
	THStack * h_shower_stack_last = new THStack();
	h_pfp_shower_nue_cc_qe_last->SetStats(kFALSE);
	h_pfp_shower_nue_cc_out_fv_last->SetStats(kFALSE);
	h_pfp_shower_nue_cc_res_last->SetStats(kFALSE);
	h_pfp_shower_nue_cc_dis_last->SetStats(kFALSE);
	h_pfp_shower_nue_cc_coh_last->SetStats(kFALSE);
	h_pfp_shower_nue_cc_mec_last->SetStats(kFALSE);
	h_pfp_shower_nc_last->SetStats(kFALSE);
	h_pfp_shower_numu_cc_qe_last->SetStats(kFALSE);
	h_pfp_shower_numu_cc_res_last->SetStats(kFALSE);
	h_pfp_shower_numu_cc_dis_last->SetStats(kFALSE);
	h_pfp_shower_numu_cc_coh_last->SetStats(kFALSE);
	h_pfp_shower_numu_cc_mec_last->SetStats(kFALSE);
	h_pfp_shower_nc_pi0_last->SetStats(kFALSE);
	h_pfp_shower_nue_cc_mixed_last->SetStats(kFALSE);
	h_pfp_shower_numu_cc_mixed_last->SetStats(kFALSE);
	h_pfp_shower_cosmic_last->SetStats(kFALSE);
	h_pfp_shower_other_mixed_last->SetStats(kFALSE);
	h_pfp_shower_unmatched_last->SetStats(kFALSE);
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
	h_shower_stack_last->Add(h_pfp_shower_nue_cc_qe_last    );
	h_shower_stack_last->Add(h_pfp_shower_nue_cc_out_fv_last);
	h_shower_stack_last->Add(h_pfp_shower_nue_cc_res_last   );
	h_shower_stack_last->Add(h_pfp_shower_nue_cc_dis_last   );
	h_shower_stack_last->Add(h_pfp_shower_nue_cc_coh_last   );
	h_shower_stack_last->Add(h_pfp_shower_nue_cc_mec_last   );
	h_shower_stack_last->Add(h_pfp_shower_nc_last       );
	h_shower_stack_last->Add(h_pfp_shower_numu_cc_qe_last   );
	h_shower_stack_last->Add(h_pfp_shower_numu_cc_res_last  );
	h_shower_stack_last->Add(h_pfp_shower_numu_cc_dis_last  );
	h_shower_stack_last->Add(h_pfp_shower_numu_cc_coh_last  );
	h_shower_stack_last->Add(h_pfp_shower_numu_cc_mec_last  );
	h_shower_stack_last->Add(h_pfp_shower_nc_pi0_last      );
	h_shower_stack_last->Add(h_pfp_shower_nue_cc_mixed_last );
	h_shower_stack_last->Add(h_pfp_shower_numu_cc_mixed_last);
	h_shower_stack_last->Add(h_pfp_shower_cosmic_last       );
	h_shower_stack_last->Add(h_pfp_shower_other_mixed_last  );
	h_shower_stack_last->Add(h_pfp_shower_unmatched_last    );
	h_shower_stack_last->Draw();
	h_shower_stack->GetXaxis()->SetTitle("Reconstructed Showers in Candidate Neutrino Object");
	TLegend * leg_shower_stack_l2 = new TLegend(0.75, 0.50, 0.95, 0.95);
	leg_shower_stack_l2->AddEntry(h_pfp_shower_nue_cc_qe_last,      "Nue CC QE", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_nue_cc_out_fv_last,  "Nue CC Out FV", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_nue_cc_res_last,     "Nue CC Res", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_nue_cc_dis_last,     "Nue CC DIS", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_nue_cc_coh_last,     "Nue CC Coh", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_nue_cc_mec_last,     "Nue CC MEC", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_nc_last,         "NC", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_numu_cc_qe_last,     "Numu CC QE", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_numu_cc_res_last,    "Numu CC Res", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_numu_cc_dis_last,    "Numu CC DIS", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_numu_cc_coh_last,    "Numu CC Coh", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_numu_cc_mec_last,    "Numu CC MEC", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_nc_pi0_last,        "NC Pi0", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_nue_cc_mixed_last,   "Nue CC Mixed", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_numu_cc_mixed_last,  "Numu CC Mixed", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_cosmic_last,         "Cosmic", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_other_mixed_last,    "Other Mixed", "f");
	leg_shower_stack_l2->AddEntry(h_pfp_shower_unmatched_last,      "Unmatched", "f");
	leg_shower_stack_l2->Draw();
	shower_stack_c2->Print("selected_pfp_shower_stack_last.pdf");

	//h_pfp_shower_data_last->Scale(1.2991626);
	h_pfp_shower_data_last->Scale(5.619);
	h_pfp_shower_data_last->Draw("P E1 same");
	shower_stack_c2->Print("data_mc_selected_pfp_shower_stack_last.pdf");

	TCanvas * track_shower_c1_last = new TCanvas();
	track_shower_c1_last->cd();
	h_pfp_track_shower_nue_cc_qe_last->Draw("colz");
	h_pfp_track_shower_nue_cc_qe_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_qe_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c1_last->Print("selected_pfp_track_shower_nue_cc_qe_last.pdf");
	TCanvas * track_shower_c2_last = new TCanvas();
	track_shower_c2_last->cd();
	h_pfp_track_shower_nue_cc_out_fv_last->Draw("colz");
	h_pfp_track_shower_nue_cc_out_fv_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_out_fv_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c2_last->Print("selected_pfp_track_shower_nue_cc_out_fv_last.pdf");
	TCanvas * track_shower_c3_last = new TCanvas();
	track_shower_c3_last->cd();
	h_pfp_track_shower_nue_cc_res_last->Draw("colz");
	h_pfp_track_shower_nue_cc_res_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_res_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c3_last->Print("selected_pfp_track_shower_nue_cc_res_last.pdf");
	TCanvas * track_shower_c4_last = new TCanvas();
	track_shower_c4_last->cd();
	h_pfp_track_shower_nue_cc_dis_last->Draw("colz");
	h_pfp_track_shower_nue_cc_dis_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_dis_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c4_last->Print("selected_pfp_track_shower_nue_cc_dis_last.pdf");
	TCanvas * track_shower_c5_last = new TCanvas();
	track_shower_c5_last->cd();
	h_pfp_track_shower_nue_cc_coh_last->Draw("colz");
	h_pfp_track_shower_nue_cc_coh_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_coh_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c5_last->Print("selected_pfp_track_shower_nue_cc_coh_last.pdf");
	TCanvas * track_shower_c6_last = new TCanvas();
	track_shower_c6_last->cd();
	h_pfp_track_shower_nue_cc_mec_last->Draw("colz");
	h_pfp_track_shower_nue_cc_mec_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_mec_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c6_last->Print("selected_pfp_track_shower_nue_cc_mec_last.pdf");
	TCanvas * track_shower_c7_last = new TCanvas();
	track_shower_c7_last->cd();
	h_pfp_track_shower_nc_last->Draw("colz");
	h_pfp_track_shower_nc_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nc_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c7_last->Print("selected_pfp_track_shower_nc_last.pdf");
	TCanvas * track_shower_c8_last = new TCanvas();
	track_shower_c8_last->cd();
	h_pfp_track_shower_numu_cc_qe_last->Draw("colz");
	h_pfp_track_shower_numu_cc_qe_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_qe_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c8_last->Print("selected_pfp_track_shower_numu_cc_qe_last.pdf");
	TCanvas * track_shower_c9_last = new TCanvas();
	track_shower_c9_last->cd();
	h_pfp_track_shower_numu_cc_res_last->Draw("colz");
	h_pfp_track_shower_numu_cc_res_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_res_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c9_last->Print("selected_pfp_track_shower_numu_cc_res_last.pdf");
	TCanvas * track_shower_c10_last = new TCanvas();
	track_shower_c10_last->cd();
	h_pfp_track_shower_numu_cc_dis_last->Draw("colz");
	h_pfp_track_shower_numu_cc_dis_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_dis_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c10_last->Print("selected_pfp_track_shower_numu_cc_dis_last.pdf");
	TCanvas * track_shower_c11_last = new TCanvas();
	track_shower_c11_last->cd();
	h_pfp_track_shower_numu_cc_coh_last->Draw("colz");
	h_pfp_track_shower_numu_cc_coh_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_coh_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c11_last->Print("selected_pfp_track_shower_numu_cc_coh_last.pdf");
	TCanvas * track_shower_c12_last = new TCanvas();
	track_shower_c12_last->cd();
	h_pfp_track_shower_numu_cc_mec_last->Draw("colz");
	h_pfp_track_shower_numu_cc_mec_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_mec_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c12_last->Print("selected_pfp_track_shower_numu_cc_mec_last.pdf");
	TCanvas * track_shower_c13_last = new TCanvas();
	track_shower_c13_last->cd();
	h_pfp_track_shower_nc_pi0_last->Draw("colz");
	h_pfp_track_shower_nc_pi0_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nc_pi0_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c13_last->Print("selected_pfp_track_shower_nc_pi0_last.pdf");
	TCanvas * track_shower_c14_last = new TCanvas();
	track_shower_c14_last->cd();
	h_pfp_track_shower_nue_cc_mixed_last->Draw("colz");
	h_pfp_track_shower_nue_cc_mixed_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_cc_mixed_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c14_last->Print("selected_pfp_track_shower_nue_cc_mixed_last.pdf");
	TCanvas * track_shower_c15_last = new TCanvas();
	track_shower_c15_last->cd();
	h_pfp_track_shower_numu_cc_mixed_last->Draw("colz");
	h_pfp_track_shower_numu_cc_mixed_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_cc_mixed_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c15_last->Print("selected_pfp_track_shower_numu_cc_mixed_last.pdf");
	TCanvas * track_shower_c16_last = new TCanvas();
	track_shower_c16_last->cd();
	h_pfp_track_shower_cosmic_last->Draw("colz");
	h_pfp_track_shower_cosmic_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_cosmic_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c16_last->Print("selected_pfp_track_shower_cosmic_last.pdf");
	TCanvas * track_shower_c17_last = new TCanvas();
	track_shower_c17_last->cd();
	h_pfp_track_shower_other_mixed_last->Draw("colz");
	h_pfp_track_shower_other_mixed_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_other_mixed_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c17_last->Print("selected_pfp_track_shower_other_mixed_last.pdf");
	TCanvas * track_shower_c18_last = new TCanvas();
	track_shower_c18_last->cd();
	h_pfp_track_shower_unmatched_last->Draw("colz");
	h_pfp_track_shower_unmatched_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_unmatched_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c18_last->Print("selected_pfp_track_shower_unmatched_last.pdf");
	TCanvas * track_shower_c19_last = new TCanvas();
	track_shower_c19_last->cd();
	h_pfp_track_shower_data_last->Draw("colz");
	h_pfp_track_shower_data_last->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_data_last->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c19_last->Print("selected_pfp_track_shower_data_last.pdf");
	//********************************************************************

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
	TCanvas * leading_c1 = new TCanvas();
	leading_c1->cd();
	h_leading_shower_mc_pdg_nue_cc_qe->Draw("colz");
	h_leading_shower_mc_pdg_nue_cc_qe->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nue_cc_qe->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nue_cc_qe->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nue_cc_qe->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c1->Print("selected_leading_shower_mc_pdg_nue_cc_qe.pdf");
	TCanvas * leading_c2 = new TCanvas();
	leading_c2->cd();
	h_leading_shower_mc_pdg_nue_cc_out_fv->Draw("colz");
	h_leading_shower_mc_pdg_nue_cc_out_fv->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nue_cc_out_fv->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nue_cc_out_fv->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nue_cc_out_fv->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c2->Print("selected_leading_shower_mc_pdg_nue_cc_out_fv.pdf");
	TCanvas * leading_c3 = new TCanvas();
	leading_c3->cd();
	h_leading_shower_mc_pdg_nue_cc_res->Draw("colz");
	h_leading_shower_mc_pdg_nue_cc_res->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nue_cc_res->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nue_cc_res->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nue_cc_res->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c3->Print("selected_leading_shower_mc_pdg_nue_cc_res.pdf");
	TCanvas * leading_c4 = new TCanvas();
	leading_c4->cd();
	h_leading_shower_mc_pdg_nue_cc_dis->Draw("colz");
	h_leading_shower_mc_pdg_nue_cc_dis->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nue_cc_dis->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nue_cc_dis->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nue_cc_dis->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c4->Print("selected_leading_shower_mc_pdg_nue_cc_dis.pdf");
	TCanvas * leading_c5 = new TCanvas();
	leading_c5->cd();
	h_leading_shower_mc_pdg_nue_cc_coh->Draw("colz");
	h_leading_shower_mc_pdg_nue_cc_coh->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nue_cc_coh->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nue_cc_coh->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nue_cc_coh->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c5->Print("selected_leading_shower_mc_pdg_nue_cc_coh.pdf");
	TCanvas * leading_c6 = new TCanvas();
	leading_c6->cd();
	h_leading_shower_mc_pdg_nue_cc_mec->Draw("colz");
	h_leading_shower_mc_pdg_nue_cc_mec->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nue_cc_mec->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nue_cc_mec->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nue_cc_mec->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c6->Print("selected_leading_shower_mc_pdg_nue_cc_mec.pdf");
	TCanvas * leading_c7 = new TCanvas();
	leading_c7->cd();
	h_leading_shower_mc_pdg_nc->Draw("colz");
	h_leading_shower_mc_pdg_nc->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nc->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nc->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nc->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c7->Print("selected_leading_shower_mc_pdg_nc.pdf");
	TCanvas * leading_c8 = new TCanvas();
	leading_c8->cd();
	h_leading_shower_mc_pdg_numu_cc_qe->Draw("colz");
	h_leading_shower_mc_pdg_numu_cc_qe->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_numu_cc_qe->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_numu_cc_qe->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_numu_cc_qe->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c8->Print("selected_leading_shower_mc_pdg_numu_cc_qe.pdf");
	TCanvas * leading_c9 = new TCanvas();
	leading_c9->cd();
	h_leading_shower_mc_pdg_numu_cc_res->Draw("colz");
	h_leading_shower_mc_pdg_numu_cc_res->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_numu_cc_res->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_numu_cc_res->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_numu_cc_res->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c9->Print("selected_leading_shower_mc_pdg_numu_cc_res.pdf");
	TCanvas * leading_c10 = new TCanvas();
	leading_c10->cd();
	h_leading_shower_mc_pdg_numu_cc_dis->Draw("colz");
	h_leading_shower_mc_pdg_numu_cc_dis->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_numu_cc_dis->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_numu_cc_dis->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_numu_cc_dis->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c10->Print("selected_leading_shower_mc_pdg_numu_cc_dis.pdf");
	TCanvas * leading_c11 = new TCanvas();
	leading_c11->cd();
	h_leading_shower_mc_pdg_numu_cc_coh->Draw("colz");
	h_leading_shower_mc_pdg_numu_cc_coh->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_numu_cc_coh->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_numu_cc_coh->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_numu_cc_coh->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c11->Print("selected_leading_shower_mc_pdg_numu_cc_coh.pdf");
	TCanvas * leading_c12 = new TCanvas();
	leading_c12->cd();
	h_leading_shower_mc_pdg_numu_cc_mec->Draw("colz");
	h_leading_shower_mc_pdg_numu_cc_mec->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_numu_cc_mec->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_numu_cc_mec->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_numu_cc_mec->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c12->Print("selected_leading_shower_mc_pdg_numu_cc_mec.pdf");
	TCanvas * leading_c13 = new TCanvas();
	leading_c13->cd();
	h_leading_shower_mc_pdg_nc_pi0->Draw("colz");
	h_leading_shower_mc_pdg_nc_pi0->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nc_pi0->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nc_pi0->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nc_pi0->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c13->Print("selected_leading_shower_mc_pdg_nc_pi0.pdf");
	TCanvas * leading_c14 = new TCanvas();
	leading_c14->cd();
	h_leading_shower_mc_pdg_nue_cc_mixed->Draw("colz");
	h_leading_shower_mc_pdg_nue_cc_mixed->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nue_cc_mixed->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nue_cc_mixed->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nue_cc_mixed->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c14->Print("selected_leading_shower_mc_pdg_nue_cc_mixed.pdf");
	TCanvas * leading_c15 = new TCanvas();
	leading_c15->cd();
	h_leading_shower_mc_pdg_numu_cc_mixed->Draw("colz");
	h_leading_shower_mc_pdg_numu_cc_mixed->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_numu_cc_mixed->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_numu_cc_mixed->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_numu_cc_mixed->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c15->Print("selected_leading_shower_mc_pdg_numu_cc_mixed.pdf");
	TCanvas * leading_c16 = new TCanvas();
	leading_c16->cd();
	h_leading_shower_mc_pdg_cosmic->Draw("colz");
	h_leading_shower_mc_pdg_cosmic->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_cosmic->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_cosmic->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_cosmic->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c16->Print("selected_leading_shower_mc_pdg_cosmic.pdf");
	TCanvas * leading_c17 = new TCanvas();
	leading_c17->cd();
	h_leading_shower_mc_pdg_other_mixed->Draw("colz");
	h_leading_shower_mc_pdg_other_mixed->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_other_mixed->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_other_mixed->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_other_mixed->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c17->Print("selected_leading_shower_mc_pdg_other_mixed.pdf");
	TCanvas * leading_c18 = new TCanvas();
	leading_c18->cd();
	h_leading_shower_mc_pdg_unmatched->Draw("colz");
	h_leading_shower_mc_pdg_unmatched->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_unmatched->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_unmatched->GetXaxis()->SetTitle("Leading Shower Origin");
	h_leading_shower_mc_pdg_unmatched->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c18->Print("selected_leading_shower_mc_pdg_unmatched.pdf");


	TCanvas * hits_eng_c1 = new TCanvas();
	hits_eng_c1->cd();
	h_shwr_hits_nu_eng->Draw("colz");
	h_shwr_hits_nu_eng->GetXaxis()->SetTitle("True Neutrino Energy [GeV]");
	h_shwr_hits_nu_eng->GetYaxis()->SetTitle("Signal Electron Shower Hits");
	hits_eng_c1->Print("shwr_hits_nu_eng.pdf");

	TCanvas * hits_eng_c2 = new TCanvas();
	hits_eng_c2->cd();
	h_shwr_hits_ele_eng->Draw("colz");
	h_shwr_hits_ele_eng->GetXaxis()->SetTitle("True Electron Energy [GeV]");
	h_shwr_hits_ele_eng->GetYaxis()->SetTitle("Signal Electron Shower Hits");
	hits_eng_c2->Print("shwr_hits_ele_eng.pdf");

	TCanvas * hits_eng_c3 = new TCanvas();
	hits_eng_c3->cd();
	h_shwr_hits_nu_eng_zoom->Draw("colz");
	h_shwr_hits_nu_eng_zoom->GetXaxis()->SetTitle("True Neutrino Energy [GeV]");
	h_shwr_hits_nu_eng_zoom->GetYaxis()->SetTitle("Signal Electron Shower Hits");
	hits_eng_c3->Print("shwr_hits_nu_eng_zoom.pdf");

	TCanvas * hits_eng_c4 = new TCanvas();
	hits_eng_c4->cd();
	h_shwr_hits_ele_eng_zoom->Draw("colz");
	h_shwr_hits_ele_eng_zoom->GetXaxis()->SetTitle("True Electron Energy [GeV]");
	h_shwr_hits_ele_eng_zoom->GetYaxis()->SetTitle("Signal Electron Shower Hits");
	hits_eng_c4->Print("shwr_hits_ele_eng_zoom.pdf");

	TCanvas * hits_eng_c1_last = new TCanvas();
	hits_eng_c1_last->cd();
	h_shwr_hits_nu_eng_last->Draw("colz");
	h_shwr_hits_nu_eng_last->GetXaxis()->SetTitle("True Selected Neutrino Energy [GeV]");
	h_shwr_hits_nu_eng_last->GetYaxis()->SetTitle("Selected Signal Electron Shower Hits");
	hits_eng_c1_last->Print("shwr_hits_nu_eng_last.pdf");

	TCanvas * hits_eng_c2_last = new TCanvas();
	hits_eng_c2_last->cd();
	h_shwr_hits_ele_eng_last->Draw("colz");
	h_shwr_hits_ele_eng_last->GetXaxis()->SetTitle("True Selected Electron Energy [GeV]");
	h_shwr_hits_ele_eng_last->GetYaxis()->SetTitle("Selected Signal Electron Shower Hits");
	hits_eng_c2_last->Print("shwr_hits_ele_eng_last.pdf");

	TCanvas * hits_eng_c3_last = new TCanvas();
	hits_eng_c3_last->cd();
	h_shwr_hits_nu_eng_zoom_last->Draw("colz");
	h_shwr_hits_nu_eng_zoom_last->GetXaxis()->SetTitle("True Selected Neutrino Energy [GeV]");
	h_shwr_hits_nu_eng_zoom_last->GetYaxis()->SetTitle("Selected Signal Electron Shower Hits");
	hits_eng_c3_last->Print("shwr_hits_nu_eng_zoom_last.pdf");

	TCanvas * hits_eng_c4_last = new TCanvas();
	hits_eng_c4_last->cd();
	h_shwr_hits_ele_eng_zoom_last->Draw("colz");
	h_shwr_hits_ele_eng_zoom_last->GetXaxis()->SetTitle("True Selected Electron Energy [GeV]");
	h_shwr_hits_ele_eng_zoom_last->GetYaxis()->SetTitle("Selected Signal Electron Shower Hits");
	hits_eng_c4_last->Print("shwr_hits_ele_eng_zoom_last.pdf");

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


	TCanvas * sequential_nu_energy_c1 = new TCanvas();
	sequential_nu_energy_c1->cd();
	h_selected_nu_energy_no_cut->SetStats(kFALSE);
	h_selected_nu_energy_no_cut->SetFillColor(29);
	h_selected_nu_energy_no_cut->GetXaxis()->SetTitle("True Signal Neutrino Energy [GeV]");
	const double nu_energy_no_cut_integral = h_selected_nu_energy_no_cut->Integral();
	h_selected_nu_energy_no_cut->Scale(1./nu_energy_no_cut_integral);
	h_selected_nu_energy_no_cut->Draw("hist");
	h_selected_nu_energy_reco_nue->SetFillColor(30);
	h_selected_nu_energy_reco_nue->SetStats(kFALSE);
	//h_selected_nu_energy_reco_nue->GetXaxis()->SetTitle("True Signal Neutrino Energy [GeV]");
	h_selected_nu_energy_reco_nue->Scale(1./nu_energy_no_cut_integral);
	h_selected_nu_energy_reco_nue->Draw("hist same");
	h_selected_nu_energy_in_fv->SetFillColor(45);
	h_selected_nu_energy_in_fv->SetStats(kFALSE);
	h_selected_nu_energy_in_fv->Scale(1./nu_energy_no_cut_integral);
	h_selected_nu_energy_in_fv->Draw("hist same");
	h_selected_nu_energy_vtx_flash->SetFillColor(28);
	h_selected_nu_energy_vtx_flash->SetStats(kFALSE);
	h_selected_nu_energy_vtx_flash->Scale(1./nu_energy_no_cut_integral);
	h_selected_nu_energy_vtx_flash->Draw("hist same");
	h_selected_nu_energy_shwr_vtx->SetFillColor(26);
	h_selected_nu_energy_shwr_vtx->SetStats(kFALSE);
	h_selected_nu_energy_shwr_vtx->Scale(1./nu_energy_no_cut_integral);
	h_selected_nu_energy_shwr_vtx->Draw("hist same");
	h_selected_nu_energy_trk_vtx->SetFillColor(36);
	h_selected_nu_energy_trk_vtx->SetStats(kFALSE);
	h_selected_nu_energy_trk_vtx->Scale(1./nu_energy_no_cut_integral);
	h_selected_nu_energy_trk_vtx->Draw("hist same");
	h_selected_nu_energy_hit_threshold->SetFillColor(39);
	h_selected_nu_energy_hit_threshold->SetStats(kFALSE);
	h_selected_nu_energy_hit_threshold->Scale(1./nu_energy_no_cut_integral);
	h_selected_nu_energy_hit_threshold->Draw("hist same");
	h_selected_nu_energy_open_angle->SetFillColor(42);
	h_selected_nu_energy_open_angle->SetStats(kFALSE);
	h_selected_nu_energy_open_angle->Scale(1./nu_energy_no_cut_integral);
	h_selected_nu_energy_open_angle->Draw("hist same");
	h_selected_nu_energy_dedx->SetFillColor(12);
	h_selected_nu_energy_dedx->SetStats(kFALSE);
	h_selected_nu_energy_dedx->Scale(1./nu_energy_no_cut_integral);
	h_selected_nu_energy_dedx->Draw("hist same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_sequential1 = new TLegend(0.65,0.65,0.85,0.85);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_sequential1->AddEntry(h_selected_nu_energy_no_cut,         "No Cuts",  "f");
	leg_sequential1->AddEntry(h_selected_nu_energy_reco_nue,       "Reco Nue", "f");
	leg_sequential1->AddEntry(h_selected_nu_energy_in_fv,          "Fiducial Volume", "f");
	leg_sequential1->AddEntry(h_selected_nu_energy_vtx_flash,      "Vtx To Flash", "f");
	leg_sequential1->AddEntry(h_selected_nu_energy_shwr_vtx,       "Shower Vtx", "f");
	leg_sequential1->AddEntry(h_selected_nu_energy_trk_vtx,        "Track Vtx", "f");
	leg_sequential1->AddEntry(h_selected_nu_energy_hit_threshold,  "Hits", "f");
	leg_sequential1->AddEntry(h_selected_nu_energy_open_angle,     "Open Angle", "f");
	leg_sequential1->AddEntry(h_selected_nu_energy_dedx,           "dE/dx", "f");
	leg_sequential1->Draw();
	sequential_nu_energy_c1->Print("sequential_nu_energy.pdf");

	TCanvas * sequential_ele_energy_c1 = new TCanvas();
	sequential_ele_energy_c1->cd();
	h_selected_ele_energy_no_cut->SetStats(kFALSE);
	h_selected_ele_energy_no_cut->SetFillColor(29);
	h_selected_ele_energy_no_cut->GetXaxis()->SetTitle("True Selected Electron Energy [GeV]");
	const double ele_energy_no_cut_integral = h_selected_ele_energy_no_cut->Integral();
	h_selected_ele_energy_no_cut->Scale(1./ele_energy_no_cut_integral);
	h_selected_ele_energy_no_cut->Draw("hist");
	h_selected_ele_energy_reco_nue->SetFillColor(30);
	h_selected_ele_energy_reco_nue->SetStats(kFALSE);
	//h_selected_ele_energy_reco_nue->GetXaxis()->SetTitle("True Signal Electron Energy [GeV]");
	h_selected_ele_energy_reco_nue->Scale(1./ele_energy_no_cut_integral);
	h_selected_ele_energy_reco_nue->Draw("hist same");
	h_selected_ele_energy_in_fv->SetFillColor(45);
	h_selected_ele_energy_in_fv->SetStats(kFALSE);
	h_selected_ele_energy_in_fv->Scale(1./ele_energy_no_cut_integral);
	h_selected_ele_energy_in_fv->Draw("hist same");
	h_selected_ele_energy_vtx_flash->SetFillColor(28);
	h_selected_ele_energy_vtx_flash->SetStats(kFALSE);
	h_selected_ele_energy_vtx_flash->Scale(1./ele_energy_no_cut_integral);
	h_selected_ele_energy_vtx_flash->Draw("hist same");
	h_selected_ele_energy_shwr_vtx->SetFillColor(26);
	h_selected_ele_energy_shwr_vtx->SetStats(kFALSE);
	h_selected_ele_energy_shwr_vtx->Scale(1./ele_energy_no_cut_integral);
	h_selected_ele_energy_shwr_vtx->Draw("hist same");
	h_selected_ele_energy_trk_vtx->SetFillColor(36);
	h_selected_ele_energy_trk_vtx->SetStats(kFALSE);
	h_selected_ele_energy_trk_vtx->Scale(1./ele_energy_no_cut_integral);
	h_selected_ele_energy_trk_vtx->Draw("hist same");
	h_selected_ele_energy_hit_threshold->SetFillColor(39);
	h_selected_ele_energy_hit_threshold->SetStats(kFALSE);
	h_selected_ele_energy_hit_threshold->Scale(1./ele_energy_no_cut_integral);
	h_selected_ele_energy_hit_threshold->Draw("hist same");
	h_selected_ele_energy_open_angle->SetFillColor(42);
	h_selected_ele_energy_open_angle->SetStats(kFALSE);
	h_selected_ele_energy_open_angle->Scale(1./ele_energy_no_cut_integral);
	h_selected_ele_energy_open_angle->Draw("hist same");
	h_selected_ele_energy_dedx->SetFillColor(12);
	h_selected_ele_energy_dedx->SetStats(kFALSE);
	h_selected_ele_energy_dedx->Scale(1./ele_energy_no_cut_integral);
	h_selected_ele_energy_dedx->Draw("hist same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_sequential2 = new TLegend(0.65,0.60,0.85,0.85);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_sequential2->AddEntry(h_selected_ele_energy_no_cut,         "No Cuts",  "f");
	leg_sequential2->AddEntry(h_selected_ele_energy_reco_nue,       "Reco Nue", "f");
	leg_sequential2->AddEntry(h_selected_ele_energy_in_fv,          "Fiducial Volume", "f");
	leg_sequential2->AddEntry(h_selected_ele_energy_vtx_flash,      "Vtx To Flash", "f");
	leg_sequential2->AddEntry(h_selected_ele_energy_shwr_vtx,       "Shower Vtx", "f");
	leg_sequential2->AddEntry(h_selected_ele_energy_trk_vtx,        "Track Vtx", "f");
	leg_sequential2->AddEntry(h_selected_ele_energy_hit_threshold,  "Hits", "f");
	leg_sequential2->AddEntry(h_selected_ele_energy_open_angle,     "Open Angle", "f");
	leg_sequential2->AddEntry(h_selected_ele_energy_dedx,           "dE/dx", "f");
	leg_sequential2->Draw();
	sequential_ele_energy_c1->Print("sequential_ele_energy.pdf");


	TCanvas * charge_share_nue_cc_mixed = new TCanvas();
	charge_share_nue_cc_mixed->cd();
	h_charge_share_nue_cc_mixed->GetXaxis()->SetTitle("Neutrino Charge Fraction - Selected Nue CC Mixed");
	h_charge_share_nue_cc_mixed->Draw();
	charge_share_nue_cc_mixed->Print("charge_fraction_nue_cc_mixed.pdf");

	TCanvas * flash_t0_diff_c1 = new TCanvas();
	flash_t0_diff_c1->cd();
	h_flash_t0_diff->GetXaxis()->SetTitle("Largest Flash Time - True Neutrino Interaction Time [us]");
	flash_t0_diff_c1->SetLogy();
	h_flash_t0_diff->Draw();
	flash_t0_diff_c1->Print("flash_t0_diff.pdf");


	TCanvas * dedx_open_angle_c1 = new TCanvas();
	dedx_open_angle_c1->cd();
	h_dedx_open_angle_nue_cc->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_nue_cc->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_nue_cc->Draw("colz");
	dedx_open_angle_c1->Print("dedx_open_angle_nue_cc.pdf");

	TCanvas * dedx_open_angle_c2 = new TCanvas();
	dedx_open_angle_c2->cd();
	h_dedx_open_angle_nue_cc_out_fv->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_nue_cc_out_fv->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_nue_cc_out_fv->Draw("colz");
	dedx_open_angle_c2->Print("dedx_open_angle_nue_cc_out_fv.pdf");

	TCanvas * dedx_open_angle_c3 = new TCanvas();
	dedx_open_angle_c3->cd();
	h_dedx_open_angle_nue_cc_mixed->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_nue_cc_mixed->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_nue_cc_mixed->Draw("colz");
	dedx_open_angle_c3->Print("dedx_open_angle_nue_cc_mixed.pdf");

	TCanvas * dedx_open_angle_c4 = new TCanvas();
	dedx_open_angle_c4->cd();
	h_dedx_open_angle_numu_cc->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_numu_cc->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_numu_cc->Draw("colz");
	dedx_open_angle_c4->Print("dedx_open_angle_numu_cc.pdf");

	TCanvas * dedx_open_angle_c5 = new TCanvas();
	dedx_open_angle_c5->cd();
	h_dedx_open_angle_numu_cc_mixed->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_numu_cc_mixed->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_numu_cc_mixed->Draw("colz");
	dedx_open_angle_c5->Print("dedx_open_angle_numu_cc_mixed.pdf");

	TCanvas * dedx_open_angle_c6 = new TCanvas();
	dedx_open_angle_c6->cd();
	h_dedx_open_angle_nc->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_nc->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_nc->Draw("colz");
	dedx_open_angle_c6->Print("dedx_open_angle_nc.pdf");

	TCanvas * dedx_open_angle_c7 = new TCanvas();
	dedx_open_angle_c7->cd();
	h_dedx_open_angle_nc_pi0->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_nc_pi0->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_nc_pi0->Draw("colz");
	dedx_open_angle_c7->Print("dedx_open_angle_nc_pi0.pdf");

	TCanvas * dedx_open_angle_c8 = new TCanvas();
	dedx_open_angle_c8->cd();
	h_dedx_open_angle_cosmic->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_cosmic->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_cosmic->Draw("colz");
	dedx_open_angle_c8->Print("dedx_open_angle_cosmic.pdf");

	TCanvas * dedx_open_angle_c9 = new TCanvas();
	dedx_open_angle_c9->cd();
	h_dedx_open_angle_other_mixed->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_other_mixed->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_other_mixed->Draw("colz");
	dedx_open_angle_c9->Print("dedx_open_angle_other_mixed.pdf");

	TCanvas * dedx_open_angle_c10 = new TCanvas();
	dedx_open_angle_c10->cd();
	h_dedx_open_angle_unmatched->GetXaxis()->SetTitle("Leading Shower dEdx [MeV/cm]");
	h_dedx_open_angle_unmatched->GetYaxis()->SetTitle("Leading Shower Open Angle [Degrees]");
	h_dedx_open_angle_unmatched->Draw("colz");
	dedx_open_angle_c10->Print("dedx_open_angle_unmatched.pdf");


	TCanvas * shwr_len_hits_c1 = new TCanvas();
	shwr_len_hits_c1->cd();
	h_shwr_len_hits_nue_cc->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_nue_cc->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_nue_cc->Draw("colz");
	shwr_len_hits_c1->Print("shwr_len_hits_nue_cc.pdf");

	TCanvas * shwr_len_hits_c2 = new TCanvas();
	shwr_len_hits_c2->cd();
	h_shwr_len_hits_nue_cc_out_fv->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_nue_cc_out_fv->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_nue_cc_out_fv->Draw("colz");
	shwr_len_hits_c2->Print("shwr_len_hits_nue_cc_out_fv.pdf");

	TCanvas * shwr_len_hits_c3 = new TCanvas();
	shwr_len_hits_c3->cd();
	h_shwr_len_hits_nue_cc_mixed->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_nue_cc_mixed->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_nue_cc_mixed->Draw("colz");
	shwr_len_hits_c3->Print("shwr_len_hits_nue_cc_mixed.pdf");

	TCanvas * shwr_len_hits_c4 = new TCanvas();
	shwr_len_hits_c4->cd();
	h_shwr_len_hits_numu_cc->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_numu_cc->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_numu_cc->Draw("colz");
	shwr_len_hits_c4->Print("shwr_len_hits_numu_cc.pdf");

	TCanvas * shwr_len_hits_c5 = new TCanvas();
	shwr_len_hits_c5->cd();
	h_shwr_len_hits_numu_cc_mixed->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_numu_cc_mixed->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_numu_cc_mixed->Draw("colz");
	shwr_len_hits_c5->Print("shwr_len_hits_numu_cc_mixed.pdf");

	TCanvas * shwr_len_hits_c6 = new TCanvas();
	shwr_len_hits_c6->cd();
	h_shwr_len_hits_nc->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_nc->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_nc->Draw("colz");
	shwr_len_hits_c6->Print("shwr_len_hits_nc.pdf");

	TCanvas * shwr_len_hits_c7 = new TCanvas();
	shwr_len_hits_c7->cd();
	h_shwr_len_hits_nc_pi0->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_nc_pi0->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_nc_pi0->Draw("colz");
	shwr_len_hits_c7->Print("shwr_len_hits_nc_pi0.pdf");

	TCanvas * shwr_len_hits_c8 = new TCanvas();
	shwr_len_hits_c8->cd();
	h_shwr_len_hits_cosmic->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_cosmic->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_cosmic->Draw("colz");
	shwr_len_hits_c8->Print("shwr_len_hits_cosmic.pdf");

	TCanvas * shwr_len_hits_c9 = new TCanvas();
	shwr_len_hits_c9->cd();
	h_shwr_len_hits_other_mixed->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_other_mixed->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_other_mixed->Draw("colz");
	shwr_len_hits_c9->Print("shwr_len_hits_other_mixed.pdf");

	TCanvas * shwr_len_hits_c10 = new TCanvas();
	shwr_len_hits_c10->cd();
	h_shwr_len_hits_unmatched->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_unmatched->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_unmatched->Draw("colz");
	shwr_len_hits_c10->Print("shwr_len_hits_unmatched.pdf");

	TCanvas * shwr_len_hits_c11 = new TCanvas();
	shwr_len_hits_c11->cd();
	h_shwr_len_hits_data->GetXaxis()->SetTitle("Leading Shower Length [cm]");
	h_shwr_len_hits_data->GetYaxis()->SetTitle("Leading Shower Hits");
	h_shwr_len_hits_data->Draw("colz");
	shwr_len_hits_c11->Print("shwr_len_hits_data.pdf");


	TCanvas * second_shwr_dist_c1 = new TCanvas();
	second_shwr_dist_c1->cd();
	h_second_shwr_dist_nue_cc->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_nue_cc->Draw();
	second_shwr_dist_c1->Print("second_shwr_dist_nue_cc.pdf");

	TCanvas * second_shwr_dist_c2 = new TCanvas();
	second_shwr_dist_c2->cd();
	h_second_shwr_dist_nue_cc_out_fv->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_nue_cc_out_fv->Draw();
	second_shwr_dist_c2->Print("second_shwr_dist_nue_cc_out_fv.pdf");

	TCanvas * second_shwr_dist_c3 = new TCanvas();
	second_shwr_dist_c3->cd();
	h_second_shwr_dist_nue_cc_mixed->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_nue_cc_mixed->Draw();
	second_shwr_dist_c3->Print("second_shwr_dist_nue_cc_mixed.pdf");

	TCanvas * second_shwr_dist_c4 = new TCanvas();
	second_shwr_dist_c4->cd();
	h_second_shwr_dist_numu_cc->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_numu_cc->Draw();
	second_shwr_dist_c4->Print("second_shwr_dist_numu_cc.pdf");

	TCanvas * second_shwr_dist_c5 = new TCanvas();
	second_shwr_dist_c5->cd();
	h_second_shwr_dist_numu_cc_mixed->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_numu_cc_mixed->Draw();
	second_shwr_dist_c5->Print("second_shwr_dist_numu_cc_mixed.pdf");

	TCanvas * second_shwr_dist_c6 = new TCanvas();
	second_shwr_dist_c6->cd();
	h_second_shwr_dist_nc->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_nc->Draw();
	second_shwr_dist_c6->Print("second_shwr_dist_nc.pdf");

	TCanvas * second_shwr_dist_c7 = new TCanvas();
	second_shwr_dist_c7->cd();
	h_second_shwr_dist_nc_pi0->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_nc_pi0->Draw();
	second_shwr_dist_c7->Print("second_shwr_dist_nc_pi0.pdf");

	TCanvas * second_shwr_dist_c8 = new TCanvas();
	second_shwr_dist_c8->cd();
	h_second_shwr_dist_cosmic->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_cosmic->Draw();
	second_shwr_dist_c8->Print("second_shwr_dist_cosmic.pdf");

	TCanvas * second_shwr_dist_c9 = new TCanvas();
	second_shwr_dist_c9->cd();
	h_second_shwr_dist_other_mixed->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_other_mixed->Draw();
	second_shwr_dist_c9->Print("second_shwr_dist_other_mixed.pdf");

	TCanvas * second_shwr_dist_c10 = new TCanvas();
	second_shwr_dist_c10->cd();
	h_second_shwr_dist_unmatched->GetXaxis()->SetTitle("(TPCO w/ > 3 Showers) Shower-Vtx Distance [cm]");
	h_second_shwr_dist_unmatched->Draw();
	second_shwr_dist_c10->Print("second_shwr_dist_unmatched.pdf");

	TCanvas * second_shwr_dist_stack_c1 = new TCanvas();
	second_shwr_dist_stack_c1->cd();
	THStack * h_second_shwr_dist_stack = new THStack();
	h_second_shwr_dist_nue_cc->SetStats(kFALSE);
	h_second_shwr_dist_nue_cc_mixed->SetStats(kFALSE);
	h_second_shwr_dist_numu_cc->SetStats(kFALSE);
	h_second_shwr_dist_nc_pi0->SetStats(kFALSE);
	h_second_shwr_dist_cosmic->SetStats(kFALSE);
	h_second_shwr_dist_nc->SetStats(kFALSE);
	h_second_shwr_dist_numu_cc_mixed->SetStats(kFALSE);
	h_second_shwr_dist_other_mixed->SetStats(kFALSE);
	h_second_shwr_dist_unmatched->SetStats(kFALSE);
	h_second_shwr_dist_nue_cc->SetFillColor(30);
	h_second_shwr_dist_nue_cc_mixed->SetFillColor(38);
	h_second_shwr_dist_numu_cc->SetFillColor(28);
	h_second_shwr_dist_nc_pi0->SetFillColor(36);
	h_second_shwr_dist_cosmic->SetFillColor(1);
	h_second_shwr_dist_nc->SetFillColor(46);
	h_second_shwr_dist_numu_cc_mixed->SetFillColor(20);
	h_second_shwr_dist_other_mixed->SetFillColor(42);
	h_second_shwr_dist_unmatched->SetFillColor(12);
	h_second_shwr_dist_stack->Add(h_second_shwr_dist_nue_cc);
	h_second_shwr_dist_stack->Add(h_second_shwr_dist_nue_cc_mixed);
	h_second_shwr_dist_stack->Add(h_second_shwr_dist_cosmic);
	h_second_shwr_dist_stack->Add(h_second_shwr_dist_numu_cc);
	h_second_shwr_dist_stack->Add(h_second_shwr_dist_numu_cc_mixed);
	h_second_shwr_dist_stack->Add(h_second_shwr_dist_nc);
	h_second_shwr_dist_stack->Add(h_second_shwr_dist_nc_pi0);
	h_second_shwr_dist_stack->Add(h_second_shwr_dist_other_mixed);
	h_second_shwr_dist_stack->Add(h_second_shwr_dist_unmatched);
	h_second_shwr_dist_stack->Draw();
	h_second_shwr_dist_stack->GetXaxis()->SetTitle("(TPCO > 3 Reco Shower) Secondary Shwr-Vtx Distance [cm]");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack_second = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack_second->AddEntry(h_second_shwr_dist_nue_cc,          "Nue CC", "f");
	leg_stack_second->AddEntry(h_second_shwr_dist_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack_second->AddEntry(h_second_shwr_dist_cosmic,          "Cosmic", "f");
	leg_stack_second->AddEntry(h_second_shwr_dist_numu_cc,         "Numu CC", "f");
	leg_stack_second->AddEntry(h_second_shwr_dist_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack_second->AddEntry(h_second_shwr_dist_nc,              "NC", "f");
	leg_stack_second->AddEntry(h_second_shwr_dist_nc_pi0,          "NC Pi0", "f");
	leg_stack_second->AddEntry(h_second_shwr_dist_other_mixed,     "Other Mixed", "f");
	leg_stack_second->AddEntry(h_second_shwr_dist_unmatched,       "Unmatched", "f");
	leg_stack_second->Draw();
	second_shwr_dist_stack_c1->Print("post_second_shwr_dist.pdf");

	double second_shwr_dist_integral_nue_cc         = h_second_shwr_dist_nue_cc->Integral();
	double second_shwr_dist_integral_nue_cc_mixed   = h_second_shwr_dist_nue_cc_mixed->Integral();
	double second_shwr_dist_integral_cosmic         = h_second_shwr_dist_cosmic->Integral();
	double second_shwr_dist_integral_numu_cc        = h_second_shwr_dist_numu_cc->Integral();
	double second_shwr_dist_integral_numu_cc_mixed  = h_second_shwr_dist_numu_cc_mixed->Integral();
	double second_shwr_dist_integral_nc             = h_second_shwr_dist_nc->Integral();
	double second_shwr_dist_integral_nc_pi0         = h_second_shwr_dist_nc_pi0->Integral();
	double second_shwr_dist_integral_other_mixed    = h_second_shwr_dist_other_mixed->Integral();
	double second_shwr_dist_integral_unmatched      = h_second_shwr_dist_unmatched->Integral();
	double second_shwr_dist_integral = second_shwr_dist_integral_nue_cc
	                                   + second_shwr_dist_integral_nue_cc_mixed
	                                   + second_shwr_dist_integral_cosmic
	                                   + second_shwr_dist_integral_numu_cc
	                                   + second_shwr_dist_integral_numu_cc_mixed
	                                   + second_shwr_dist_integral_nc
	                                   + second_shwr_dist_integral_nc_pi0
	                                   + second_shwr_dist_integral_other_mixed
	                                   + second_shwr_dist_integral_unmatched;

	//h_second_shwr_dist_data->Scale(1.2991626);
	//h_second_shwr_dist_data->Scale(5.619);
	h_second_shwr_dist_data->Scale((second_shwr_dist_integral / h_second_shwr_dist_data->Integral()) * data_mc_scale_factor);
	h_second_shwr_dist_data->Draw("P E1 same");
	second_shwr_dist_stack_c1->Print("data_mc_post_second_shwr_dist.pdf");

	TCanvas * hit_length_ratio_stack_c1 = new TCanvas();
	hit_length_ratio_stack_c1->cd();
	THStack * h_hit_length_ratio_stack = new THStack();
	h_hit_length_ratio_nue_cc->SetStats(kFALSE);
	h_hit_length_ratio_nue_cc_mixed->SetStats(kFALSE);
	h_hit_length_ratio_numu_cc->SetStats(kFALSE);
	h_hit_length_ratio_nc_pi0->SetStats(kFALSE);
	h_hit_length_ratio_cosmic->SetStats(kFALSE);
	h_hit_length_ratio_nc->SetStats(kFALSE);
	h_hit_length_ratio_numu_cc_mixed->SetStats(kFALSE);
	h_hit_length_ratio_other_mixed->SetStats(kFALSE);
	h_hit_length_ratio_unmatched->SetStats(kFALSE);
	h_hit_length_ratio_nue_cc->SetFillColor(30);
	h_hit_length_ratio_nue_cc_mixed->SetFillColor(38);
	h_hit_length_ratio_numu_cc->SetFillColor(28);
	h_hit_length_ratio_nc_pi0->SetFillColor(36);
	h_hit_length_ratio_cosmic->SetFillColor(1);
	h_hit_length_ratio_nc->SetFillColor(46);
	h_hit_length_ratio_numu_cc_mixed->SetFillColor(20);
	h_hit_length_ratio_other_mixed->SetFillColor(42);
	h_hit_length_ratio_unmatched->SetFillColor(12);
	h_hit_length_ratio_stack->Add(h_hit_length_ratio_nue_cc);
	h_hit_length_ratio_stack->Add(h_hit_length_ratio_nue_cc_mixed);
	h_hit_length_ratio_stack->Add(h_hit_length_ratio_cosmic);
	h_hit_length_ratio_stack->Add(h_hit_length_ratio_numu_cc);
	h_hit_length_ratio_stack->Add(h_hit_length_ratio_numu_cc_mixed);
	h_hit_length_ratio_stack->Add(h_hit_length_ratio_nc);
	h_hit_length_ratio_stack->Add(h_hit_length_ratio_nc_pi0);
	h_hit_length_ratio_stack->Add(h_hit_length_ratio_other_mixed);
	h_hit_length_ratio_stack->Add(h_hit_length_ratio_unmatched);
	h_hit_length_ratio_stack->Draw();
	h_hit_length_ratio_stack->GetXaxis()->SetTitle("Leading Shower (Hits / Length) [cm^-1]");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack_hit_length_ratio = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack_hit_length_ratio->AddEntry(h_hit_length_ratio_nue_cc,          "Nue CC", "f");
	leg_stack_hit_length_ratio->AddEntry(h_hit_length_ratio_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack_hit_length_ratio->AddEntry(h_hit_length_ratio_cosmic,          "Cosmic", "f");
	leg_stack_hit_length_ratio->AddEntry(h_hit_length_ratio_numu_cc,         "Numu CC", "f");
	leg_stack_hit_length_ratio->AddEntry(h_hit_length_ratio_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack_hit_length_ratio->AddEntry(h_hit_length_ratio_nc,              "NC", "f");
	leg_stack_hit_length_ratio->AddEntry(h_hit_length_ratio_nc_pi0,          "NC Pi0", "f");
	leg_stack_hit_length_ratio->AddEntry(h_hit_length_ratio_other_mixed,     "Other Mixed", "f");
	leg_stack_hit_length_ratio->AddEntry(h_hit_length_ratio_unmatched,       "Unmatched", "f");
	leg_stack_hit_length_ratio->Draw();
	hit_length_ratio_stack_c1->Print("post_hit_length_ratio.pdf");

	double hit_length_ratio_integral_nue_cc         = h_hit_length_ratio_nue_cc->Integral();
	double hit_length_ratio_integral_nue_cc_mixed   = h_hit_length_ratio_nue_cc_mixed->Integral();
	double hit_length_ratio_integral_cosmic         = h_hit_length_ratio_cosmic->Integral();
	double hit_length_ratio_integral_numu_cc        = h_hit_length_ratio_numu_cc->Integral();
	double hit_length_ratio_integral_numu_cc_mixed  = h_hit_length_ratio_numu_cc_mixed->Integral();
	double hit_length_ratio_integral_nc             = h_hit_length_ratio_nc->Integral();
	double hit_length_ratio_integral_nc_pi0         = h_hit_length_ratio_nc_pi0->Integral();
	double hit_length_ratio_integral_other_mixed    = h_hit_length_ratio_other_mixed->Integral();
	double hit_length_ratio_integral_unmatched      = h_hit_length_ratio_unmatched->Integral();
	double hit_length_ratio_integral = hit_length_ratio_integral_nue_cc
	                                   + hit_length_ratio_integral_nue_cc_mixed
	                                   + hit_length_ratio_integral_cosmic
	                                   + hit_length_ratio_integral_numu_cc
	                                   + hit_length_ratio_integral_numu_cc_mixed
	                                   + hit_length_ratio_integral_nc
	                                   + hit_length_ratio_integral_nc_pi0
	                                   + hit_length_ratio_integral_other_mixed
	                                   + hit_length_ratio_integral_unmatched;

	//h_hit_length_ratio_data->Scale(1.2991626);
	//h_hit_length_ratio_data->Scale(5.619);
	h_hit_length_ratio_data->Scale((hit_length_ratio_integral / h_hit_length_ratio_data->Integral()) * data_mc_scale_factor);
	h_hit_length_ratio_data->Draw("P E1 same");
	hit_length_ratio_stack_c1->Print("data_mc_post_hit_length_ratio.pdf");

	TCanvas * failure_reason_stack_c1 = new TCanvas();
	failure_reason_stack_c1->cd();
	THStack * h_failure_reason_stack = new THStack();
	h_failure_reason_data->SetStats(kFALSE);
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

	std::vector<double> remaining_data;
	std::vector<double> remaining_nue_cc;
	std::vector<double> remaining_nue_cc_mixed;
	std::vector<double> remaining_numu_cc;
	std::vector<double> remaining_numu_cc_mixed;
	std::vector<double> remaining_cosmic;
	std::vector<double> remaining_nc;
	std::vector<double> remaining_nc_pi0;
	std::vector<double> remaining_other_mixed;
	std::vector<double> remaining_unmatched;

	double bin_content_data_prev          = data_pe_counter_v->at(5);
	double bin_content_nue_cc_prev        = pe_counter_v->at(0);
	double bin_content_nue_cc_mixed_prev  = pe_counter_v->at(1);
	double bin_content_numu_cc_prev       = pe_counter_v->at(4);
	double bin_content_nc_pi0_prev        = pe_counter_v->at(10);
	double bin_content_cosmic_prev        = pe_counter_v->at(2);
	double bin_content_nc_prev            = pe_counter_v->at(3);
	double bin_content_numu_cc_mixed_prev = pe_counter_v->at(11);
	double bin_content_other_mixed_prev   = pe_counter_v->at(6);
	double bin_content_unmatched_prev     = pe_counter_v->at(5);

	for (int i=1; i<= 22; i++)
	{
		double bin_content_data                     = h_failure_reason_data->GetBinContent(i);
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
			const double ratio_data          = double(bin_content_data)   / double(bin_content_data_prev);
			const double ratio_nue_cc        = double(bin_content_nue_cc) / double(bin_content_nue_cc_prev);
			const double ratio_nue_cc_mixed  = double(bin_content_nue_cc_mixed) / double(bin_content_nue_cc_mixed_prev);
			const double ratio_numu_cc       = double(bin_content_numu_cc) / double(bin_content_numu_cc_prev);
			const double ratio_numu_cc_mixed = double(bin_content_numu_cc_mixed) / double(bin_content_numu_cc_mixed_prev);
			const double ratio_cosmic        = double(bin_content_cosmic) / double(bin_content_cosmic_prev);
			const double ratio_nc            = double(bin_content_nc) / double(bin_content_nc_prev);
			const double ratio_nc_pi0        = double(bin_content_nc_pi0) / double(bin_content_nc_pi0_prev);
			const double ratio_other_mixed   = double(bin_content_other_mixed) / double(bin_content_other_mixed_prev);
			const double ratio_unmatched     = double(bin_content_unmatched) / double(bin_content_unmatched_prev);

			bin_content_data_prev                -= bin_content_data;
			bin_content_nue_cc_prev              -= bin_content_nue_cc;
			bin_content_nue_cc_mixed_prev        -= bin_content_nue_cc_mixed;
			bin_content_numu_cc_prev             -= bin_content_numu_cc;
			bin_content_nc_pi0_prev              -= bin_content_nc_pi0;
			bin_content_cosmic_prev              -= bin_content_cosmic;
			bin_content_nc_prev                  -= bin_content_nc;
			bin_content_numu_cc_mixed_prev       -= bin_content_numu_cc_mixed;
			bin_content_other_mixed_prev         -= bin_content_other_mixed;
			bin_content_unmatched_prev           -= bin_content_unmatched;

			remaining_data.push_back(ratio_data);
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
		h_failure_reason_data->SetBinContent((i*2)+1,           remaining_data.at(i));
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

	h_failure_reason_stack->Add(h_failure_reason_data);
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
	const char * str_cut[22] = {"HasNue", " ", "InFV", " ", "FlshDist", " ", "ShwrVtx", " ", "TrkVtx", " ", "Hits", " ",
		                    "OpenAngle", " ", "dEdx", " ", "2ndDist", " ", "HitLength", " ", "Passed", " "};
	for (int i=1; i<= 22; i++)
	{
		h_failure_reason_stack->GetXaxis()->SetBinLabel(i,str_cut[(i-1)]);
	}

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack_failure_reason = new TLegend(0.865,0.75,0.99,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack_failure_reason->AddEntry(h_failure_reason_data,          "Data", "f");
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
	//if(argc == 2 ) { std::cout << "Please Include Both the MC and Data File" << std::endl; exit(1); }
	if(argc != 3 ) { std::cout << "Please inclue the input file path for MC and Data File" << std::endl; exit(1); }
	const char * file1 = argv[1];
	const char * file2 = argv[2];

	return xsecSelection::selection(file1, file2);
}

#endif
