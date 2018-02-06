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
	std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * data_post_cuts_v
	        = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> >;

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

	int data_run;
	int data_subrun;
	int last_data_run = 0;
	int last_data_subrun = 0;

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
		 //std::vector<std::string> * tpco_origin_v = new std::vector<std::string>;
		 //_cuts_instance.selection_cuts::GetOrigins(data_tpc_object_container_v, tpco_origin_v);
		 //**********************************************************************
		 //Get Run and Subrun*
		 //run and subrun number should be the same for all tpc obj in one event
		 //**********************************************************************
		 //some events have no tpc objects...
		for(auto const tpc_obj : * tpc_object_container_v)
		{
			data_run    = tpc_obj.RunNumber();
			data_subrun = tpc_obj.SubRunNumber();
			std::cout << data_run << " " << data_subrun << std::endl;
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
	std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v
	        = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> >;

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

		// _functions_instance.selection_functions::TopologyEfficiency(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		//                                                             _x1, _x2, _y1, _y2, _z1, _z2,
		//                                                             mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		//                                                             no_track, has_track);

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

	histogram_functions::Plot1DHistogram (h_nue_eng_eff_num, "True Neutrino Energy [GeV]", "selected_true_neutrino_energy.pdf");
	histogram_functions::Plot1DHistogram (h_nue_eng_eff_den, "True Neutrino Energy [GeV]", "all_true_neutrino_energy.pdf");
	histogram_functions::Plot1DHistogram (h_nue_num_part_eff_den, "True Particle Multiplicity", "all_true_neutrino_num_particles.pdf");
	histogram_functions::Plot1DHistogram (h_nue_num_part_eff_num, "Selected True Particle Multiplicity", "selected_true_neutrino_num_particles.pdf");
	histogram_functions::Plot1DHistogram (h_nue_num_chrg_part_eff_den, "True Charged Particle Multiplicity", "all_true_neutrino_num_charged_particles.pdf");
	histogram_functions::Plot1DHistogram (h_nue_num_chrg_part_eff_num, "Selected Charged True Particle Multiplicity",
	                                      "selected_true_neutrino_num_charged_particles.pdf");

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
	histogram_functions::PlotTEfficiency (h_ele_eng_eff_num, h_ele_eng_eff_den, ";True Nue Electron Energy;Efficiency", "signal_selection_ele_energy_efficiency.pdf");



	double all_energy = 0;
	for(auto const energy : selected_energy_vector) {all_energy += energy; }
	const double average_true_energy = all_energy / selected_energy_vector.size();
	_functions_instance.selection_functions::xsec_plot(_verbose, genie_xsec, xsec_cc->at(0), average_true_energy, xsec_cc->at(1));

	histogram_functions::Plot2DHistogram (h_tracks_showers, "Post Cuts - Showers/Tracks per Candidate Nue TPC Object",
	                                      "Reco Tracks", "Reco Showers", "post_cuts_showers_tracks.pdf");
	histogram_functions::Plot2DHistogram (h_tracks_showers, "Post Cuts - Showers/Tracks per Candidate Nue TPC Object",
	                                      "Reco Tracks", "Reco Showers", "post_cuts_showers_tracks_cosmic.pdf");
	histogram_functions::Plot2DHistogram (h_tracks_showers, "Post Cuts - Showers/Tracks per Candidate Nue TPC Object",
	                                      "Reco Tracks", "Reco Showers", "post_cuts_showers_tracks_numu.pdf");

	histogram_functions::PlotDataMC(h_leading_shower_open_angle_nue_cc, h_leading_shower_open_angle_nue_cc_mixed,
	                                h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_numu_cc_mixed,
	                                h_leading_shower_open_angle_cosmic, h_leading_shower_open_angle_nc,
	                                h_leading_shower_open_angle_nc_pi0, h_leading_shower_open_angle_other_mixed,
	                                h_leading_shower_open_angle_unmatched, h_leading_shower_open_angle_data, data_mc_scale_factor,
	                                "Shower Opening Angle [Degrees]", "",
	                                "post_cuts_leading_shower_open_angle.pdf", "data_mc_post_cuts_leading_shower_open_angle.pdf");

	histogram_functions::PlotDataMC(h_dedx_cuts_nue_cc, h_dedx_cuts_nue_cc_mixed,
	                                h_dedx_cuts_numu_cc, h_dedx_cuts_numu_cc_mixed,
	                                h_dedx_cuts_cosmic, h_dedx_cuts_nc,
	                                h_dedx_cuts_nc_pi0, h_dedx_cuts_other_mixed,
	                                h_dedx_cuts_unmatched, h_dedx_cuts_data, data_mc_scale_factor,
	                                "Collection Plane dE/dx [MeV/cm]", "",
	                                "post_cuts_dedx_cuts.pdf", "data_mc_post_cuts_dedx_cuts.pdf");

	histogram_functions::PlotDataMC(h_vtx_flash_nue_cc, h_vtx_flash_nue_cc_mixed,
	                                h_vtx_flash_numu_cc, h_vtx_flash_numu_cc_mixed,
	                                h_vtx_flash_cosmic, h_vtx_flash_nc,
	                                h_vtx_flash_nc_pi0, h_vtx_flash_other_mixed,
	                                h_vtx_flash_unmatched, h_vtx_flash_data, data_mc_scale_factor,
	                                "2D Distance From Largest Flash to Reco Nu Vtx [cm]", "",
	                                "post_cuts_vtx_to_flash_distance.pdf", "data_mc_post_cuts_vtx_to_flash_distance.pdf");

	histogram_functions::PlotDataMC(h_trk_vtx_dist_nue_cc, h_trk_vtx_dist_nue_cc_mixed,
	                                h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_cc_mixed,
	                                h_trk_vtx_dist_cosmic, h_trk_vtx_dist_nc,
	                                h_trk_vtx_dist_nc_pi0, h_trk_vtx_dist_other_mixed,
	                                h_trk_vtx_dist_unmatched, h_trk_vtx_dist_data, data_mc_scale_factor,
	                                "Track to Nue Candidate Vertex Distance [cm]", "",
	                                "post_cuts_track_to_vtx.pdf", "data_mc_post_cuts_track_to_vtx.pdf");

	histogram_functions::PlotDataMC(h_shwr_vtx_dist_nue_cc, h_shwr_vtx_dist_nue_cc_mixed,
	                                h_shwr_vtx_dist_numu_cc, h_shwr_vtx_dist_numu_cc_mixed,
	                                h_shwr_vtx_dist_cosmic, h_shwr_vtx_dist_nc,
	                                h_shwr_vtx_dist_nc_pi0, h_shwr_vtx_dist_other_mixed,
	                                h_shwr_vtx_dist_unmatched, h_shwr_vtx_dist_data, data_mc_scale_factor,
	                                "Shower to Nue Candidate Vertex Distance [cm]", "",
	                                "post_cuts_shower_to_vtx.pdf", "data_mc_post_cuts_shower_to_vtx.pdf");


	histogram_functions::PlotDetailDataMCStack(h_pfp_track_nue_cc_qe,
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
	                                           h_pfp_track_data,
	                                           data_mc_scale_factor,
	                                           "Reconstructed Tracks in Candidate Neutrino Object", "",
	                                           "selected_pfp_track_stack.pdf", "data_mc_selected_pfp_track_stack.pdf");

	histogram_functions::PlotDetailDataMCStack(h_pfp_shower_nue_cc_qe,
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
	                                           h_pfp_shower_data,
	                                           data_mc_scale_factor,
	                                           "Reconstructed Showers in Candidate Neutrino Object", "",
	                                           "selected_pfp_shower_stack.pdf", "data_mc_selected_pfp_shower_stack.pdf");

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
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_data, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_data.pdf");

	histogram_functions::PlotDetailDataMCStack(h_pfp_track_nue_cc_qe_last,
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
	                                           h_pfp_track_data_last,
	                                           data_mc_scale_factor,
	                                           "Reconstructed Tracks in Candidate Neutrino Object", "",
	                                           "selected_pfp_track_stack_last.pdf", "data_mc_selected_pfp_track_stack_last.pdf");

	histogram_functions::PlotDetailDataMCStack(h_pfp_shower_nue_cc_qe_last,
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
	                                           h_pfp_shower_data_last,
	                                           data_mc_scale_factor,
	                                           "Reconstructed Showers in Candidate Neutrino Object", "",
	                                           "selected_pfp_shower_stack_last.pdf", "data_mc_selected_pfp_shower_stack_last.pdf");

	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_qe_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_qe_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_out_fv_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_out_fv_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_res_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_res_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_dis_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_dis_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_coh_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_coh_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mec_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_mec_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nue_cc_mixed_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nue_cc_mixed_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_qe_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_qe_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_res_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_res_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_dis_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_dis_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_coh_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_coh_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mec_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_mec_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_numu_cc_mixed_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_numu_cc_mixed_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nc_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_nc_pi0_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_nc_pi0_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_cosmic_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_cosmic_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_other_mixed_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_other_mixed_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_unmatched_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_unmatched_last.pdf");
	histogram_functions::Plot2DHistogram (h_pfp_track_shower_data_last, "", "PFP Tracks", "PFP Showers", "selected_pfp_track_shower_data_last.pdf");
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


	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nue_cc, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_nue_cc.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nue_cc_out_fv, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_nue_cc_out_fv.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nue_cc_mixed, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_nue_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_numu_cc, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_numu_cc.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_numu_cc_mixed, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_numu_cc_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nc, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_nc.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_nc_pi0, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_nc_pi0.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_cosmic, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_cosmic.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_other_mixed, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_other_mixed.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_unmatched, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_unmatched.pdf");
	histogram_functions::Plot2DHistogram (h_shwr_len_hits_data, "", "Leading Shower Length [cm]", "Leading Shower Hits", "shwr_len_hits_data.pdf");

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


	histogram_functions::PlotDataMC(h_second_shwr_dist_nue_cc,  h_second_shwr_dist_nue_cc_mixed,
	                                h_second_shwr_dist_numu_cc, h_second_shwr_dist_numu_cc_mixed,
	                                h_second_shwr_dist_cosmic,  h_second_shwr_dist_nc,
	                                h_second_shwr_dist_nc_pi0,  h_second_shwr_dist_other_mixed,
	                                h_second_shwr_dist_unmatched, h_second_shwr_dist_data, data_mc_scale_factor,
	                                "(TPCO > 3 Reco Shower) Secondary Shwr-Vtx Distance [cm]", "",
	                                "post_second_shwr_dist.pdf", "data_mc_post_second_shwr_dist.pdf");

	histogram_functions::PlotDataMC(h_hit_length_ratio_nue_cc,  h_hit_length_ratio_nue_cc_mixed,
	                                h_hit_length_ratio_numu_cc, h_hit_length_ratio_numu_cc_mixed,
	                                h_hit_length_ratio_cosmic,  h_hit_length_ratio_nc,
	                                h_hit_length_ratio_nc_pi0,  h_hit_length_ratio_other_mixed,
	                                h_hit_length_ratio_unmatched, h_hit_length_ratio_data, data_mc_scale_factor,
	                                "Leading Shower (Hits / Length) [cm^-1]", "",
	                                "post_hit_length_ratio.pdf", "data_mc_post_hit_length_ratio.pdf");

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
