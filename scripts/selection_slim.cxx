#include "selection_slim.h"

namespace xsecSelection {
int selection( const char * _file1){

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

	//Event, Run, VtxX, VtxY, VtxZ, pass/fail reason
	std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v
	        = new std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> >;

	std::cout << "=====================" << std::endl;
	std::cout << "== Begin Selection ==" << std::endl;
	std::cout << "=====================" << std::endl;

	const int total_entries = mytree->GetEntries();
	std::cout << "Total Events: " << total_entries << std::endl;

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

		if(mc_nu_id == 1 || mc_nu_id == 5)
		//if this event is a true nue CC interaction and is inside the FV
		//also include nue_cc_bar as in the tpco classification I use the nue-bar as well
		{
			if(_cuts_instance.selection_cuts::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, _x1, _x2, _y1, _y2, _z1, _z2) == true) { total_mc_entries_inFV++; }
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

		//************************
		//******** in fv cut *****
		//************************
		_cuts_instance.selection_cuts::fiducial_volume_cut(tpc_object_container_v, _x1, _x2, _y1, _y2, _z1, _z2, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, in_fv_counter_v);

		//*****************************
		//**** vertex to flash cut ****
		//*****************************
		_cuts_instance.selection_cuts::flashRecoVtxDist(largest_flash_v, tpc_object_container_v, tolerance, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, vtx_flash_counter_v);

		//******************************************************
		//*** distance between pfp shower and nue object cut ***
		//******************************************************
		_cuts_instance.selection_cuts::VtxNuDistance(tpc_object_container_v, shwr_nue_tolerance, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, shwr_tpco_counter_v);

		//******************************************************
		// **** distance between pfp track and nue object cut **
		//******************************************************
		_cuts_instance.selection_cuts::VtxTrackNuDistance(tpc_object_container_v, trk_nue_tolerance, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, trk_tpco_counter_v);

		//****************************************************
		// ******** hit threshold for showers cut *************
		//******************************************************
		_cuts_instance.selection_cuts::HitThreshold(tpc_object_container_v, shwr_hit_threshold, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_threshold_counter_v);
		//*****************************************************
		//****** open angle cut for the leading shower ********
		//******************************************************
		_cuts_instance.selection_cuts::OpenAngleCut(tpc_object_container_v, passed_tpco, tolerance_open_angle, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, open_angle_counter_v);

		//*****************************************************
		//*********** dEdx cut for the leading shower *********
		//******************************************************
		_cuts_instance.selection_cuts::dEdxCut(tpc_object_container_v, passed_tpco, tolerance_dedx_min, tolerance_dedx_max, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, dedx_counter_v);
		//*************************************
		// ******** End Selection Cuts! ******
		//*************************************
		_functions_instance.selection_functions::TopologyEfficiency(tpc_object_container_v, passed_tpco, _verbose, has_pi0,
		                                                            _x1, _x2, _y1, _y2, _z1, _z2,
		                                                            mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                            no_track, has_track, _1_shwr, _2_shwr, _3_shwr, _4_shwr);

		_functions_instance.selection_functions::FillPostCutVector(tpc_object_container_v, passed_tpco, has_pi0,
		                                                           _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, post_cuts_v);

//***************************************************************************
		_cuts_instance.selection_cuts::SecondaryShowersDistCut(tpc_object_container_v, passed_tpco, _verbose, dist_tolerance);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, secondary_shower_counter_v);

		_cuts_instance.selection_cuts::HitLengthRatioCut(tpc_object_container_v, passed_tpco, _verbose, pfp_hits_length_tolerance);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_lengthRatio_counter_v);

		_cuts_instance.selection_cuts::HitThresholdCollection(tpc_object_container_v, shwr_hit_threshold_collection, passed_tpco, _verbose);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_threshold_collection_counter_v);

		//*** cut for longest track / leading shower ratio *** //
		_cuts_instance.selection_cuts::LongestTrackLeadingShowerCut(tpc_object_container_v, passed_tpco, _verbose, ratio_tolerance);
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco, has_pi0,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, trk_len_shwr_len_ratio_counter_v);


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
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, hit_threshold_collection_counter_v, "WPlane Hit Threshold");
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV, trk_len_shwr_len_ratio_counter_v,   "TrkLen/ShwrLen Ratio");

	std::cout << "---------------------" << std::endl;
	std::cout << "No Track Signal: " << no_track->at(0) << std::endl;
	std::cout << "No Track Bkg   : " << no_track->at(1) << std::endl;
	std::cout << "Purity         : " << double(no_track->at(0)) / double(no_track->at(0) + no_track->at(1)) << std::endl;
	std::cout << " ******************* " << std::endl;
	std::cout << "1+ Track Signal: " << has_track->at(0) << std::endl;
	std::cout << "1+ Track Bkg   : " << has_track->at(1) << std::endl;
	std::cout << "Purity         : " << double(has_track->at(0)) / double(has_track->at(0) + has_track->at(1)) << std::endl;
	std::cout << "---------------------" << std::endl;
	std::cout << "---------------------" << std::endl;
	std::cout << "1 Shower Signal : " << _1_shwr->at(0) << std::endl;
	std::cout << "1 Shower Bkg    : " << _1_shwr->at(1) << std::endl;
	std::cout << "Purity          : " << double(_1_shwr->at(0)) / double(_1_shwr->at(0) + _1_shwr->at(1)) << std::endl;
	std::cout << " ******************* " << std::endl;
	std::cout << "2 Shower Signal : " << _2_shwr->at(0) << std::endl;
	std::cout << "2 Shower Bkg    : " << _2_shwr->at(1) << std::endl;
	std::cout << "Purity          : " << double(_2_shwr->at(0)) / double(_2_shwr->at(0) + _2_shwr->at(1)) << std::endl;
	std::cout << " ******************* " << std::endl;
	std::cout << "3 Shower Signal : " << _3_shwr->at(0) << std::endl;
	std::cout << "3 Shower Bkg    : " << _3_shwr->at(1) << std::endl;
	std::cout << "Purity          : " << double(_3_shwr->at(0)) / double(_3_shwr->at(0) + _3_shwr->at(1)) << std::endl;
	std::cout << " ******************* " << std::endl;
	std::cout << "4+ Shower Signal: " << _4_shwr->at(0) << std::endl;
	std::cout << "4+ Shower Bkg   : " << _4_shwr->at(1) << std::endl;
	std::cout << "Purity          : " << double(_4_shwr->at(0)) / double(_4_shwr->at(0) + _4_shwr->at(1)) << std::endl;
	std::cout << "---------------------" << std::endl;
	std::cout << "---------------------" << std::endl;

	std::vector<double> * xsec_cc = new std::vector<double>;
	const double final_counter                = dedx_counter_v->at(7);
	const double final_counter_nue_cc         = dedx_counter_v->at(0);
	const double final_counter_nue_cc_mixed   = dedx_counter_v->at(1);
	const double final_counter_nue_cc_out_fv  = dedx_counter_v->at(9);
	const double final_counter_cosmic         = dedx_counter_v->at(2);
	const double final_counter_nc             = dedx_counter_v->at(3);
	const double final_counter_numu_cc        = dedx_counter_v->at(4);
	const double final_counter_numu_cc_mixed  = dedx_counter_v->at(11);
	const double final_counter_nc_pi0         = dedx_counter_v->at(10);
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
	if(argc != 2 ) { std::cout << "Please inclue the input file path" << std::endl; exit(1); }
	const char * file1 = argv[1];


	return xsecSelection::selection(file1);
}

#endif
