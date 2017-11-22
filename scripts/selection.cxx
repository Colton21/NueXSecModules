#include "selection.h"

namespace xsecSelection {
int selection( const char * _file1){

	std::cout << "File Path: " << _file1 << std::endl;
	const bool _verbose = false;
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

	std::cout << "=====================" << std::endl;
	std::cout << "== Begin Selection ==" << std::endl;
	std::cout << "=====================" << std::endl;

	std::vector<int> * passed_runs = new std::vector<int>;
	const int total_entries = mytree->GetEntries();
	std::cout << "Total Events: " << total_entries << std::endl;

	//let's do the in-time cut as the very first thing
	std::cout << "=====================" << std::endl;
	std::cout << "==== In Time Cut ====" << std::endl;
	std::cout << "=====================" << std::endl;

	_functions_instance.selection_functions::loop_flashes(f, optree, flash_pe_threshold, flash_time_start,
	                                                      flash_time_end, passed_runs);

	//get vector with largest flashes y,z positions
	std::vector< std::vector< double> > * largest_flash_v_v = new std::vector < std::vector < double > >;
	_functions_instance.selection_functions::SetXYflashVector(f, optree, largest_flash_v_v);
	std::cout << "Largest Flash Vector Size: " << largest_flash_v_v->size() << std::endl;

	for(auto const run : * passed_runs) {run_sum = run_sum + run; }
	std::cout << "Passed Runs Vector Size: " << passed_runs->size() << std::endl;
	std::cout << "Number Passed Events: " << run_sum << std::endl;

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

		mctruth_counter_tree->GetEntry(event);
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
			if(_functions_instance.selection_functions::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
			                                                  _x1, _x2, _y1,
			                                                  _y2, _z1, _z2) == true)
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

		mytree->GetEntry(event);
		//this is where the in-time optical cut actually takes effect
		if(passed_runs->at(event) == 0)
		{
			if(_verbose) std::cout << "[Failed In-Time Cut]" << std::endl;
			continue;
		}//false

		std::vector<std::string> *tpco_origin_v = new std::vector<std::string>;
		_functions_instance.selection_functions::GetOrigins(tpc_object_container_v, tpco_origin_v);

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

		//** start the cuts here **

		//reco nue cut
		_functions_instance.selection_functions::HasNue(tpc_object_container_v, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, reco_nue_counter_v);

		//this doesn't normally sit here!
		_functions_instance.selection_functions::TopologyPlots(tpc_object_container_v, passed_tpco,
		                                                       _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		                                                       h_pfp_track_shower_nue_cc_qe,
		                                                       h_pfp_track_shower_nue_cc_out_fv,
		                                                       h_pfp_track_shower_nue_cc_res,
		                                                       h_pfp_track_shower_nue_cc_dis,
		                                                       h_pfp_track_shower_nue_cc_coh,
		                                                       h_pfp_track_shower_nue_cc_mec,
		                                                       h_pfp_track_shower_nue_nc,
		                                                       h_pfp_track_shower_numu_cc_qe,
		                                                       h_pfp_track_shower_numu_cc_res,
		                                                       h_pfp_track_shower_numu_cc_dis,
		                                                       h_pfp_track_shower_numu_cc_coh,
		                                                       h_pfp_track_shower_numu_cc_mec,
		                                                       h_pfp_track_shower_numu_nc,
		                                                       h_pfp_track_shower_nue_cc_mixed,
		                                                       h_pfp_track_shower_numu_cc_mixed,
		                                                       h_pfp_track_shower_cosmic,
		                                                       h_pfp_track_shower_other_mixed,
		                                                       h_pfp_track_shower_unmatched,
		                                                       h_leading_shower_mc_pdg_nue_cc_qe,
		                                                       h_leading_shower_mc_pdg_nue_cc_out_fv,
		                                                       h_leading_shower_mc_pdg_nue_cc_res,
		                                                       h_leading_shower_mc_pdg_nue_cc_dis,
		                                                       h_leading_shower_mc_pdg_nue_cc_coh,
		                                                       h_leading_shower_mc_pdg_nue_cc_mec,
		                                                       h_leading_shower_mc_pdg_nue_nc,
		                                                       h_leading_shower_mc_pdg_numu_cc_qe,
		                                                       h_leading_shower_mc_pdg_numu_cc_res,
		                                                       h_leading_shower_mc_pdg_numu_cc_dis,
		                                                       h_leading_shower_mc_pdg_numu_cc_coh,
		                                                       h_leading_shower_mc_pdg_numu_cc_mec,
		                                                       h_leading_shower_mc_pdg_numu_nc,
		                                                       h_leading_shower_mc_pdg_nue_cc_mixed,
		                                                       h_leading_shower_mc_pdg_numu_cc_mixed,
		                                                       h_leading_shower_mc_pdg_cosmic,
		                                                       h_leading_shower_mc_pdg_other_mixed,
		                                                       h_leading_shower_mc_pdg_unmatched,
		                                                       h_pfp_track_nue_cc_qe,
		                                                       h_pfp_track_nue_cc_out_fv,
		                                                       h_pfp_track_nue_cc_res,
		                                                       h_pfp_track_nue_cc_dis,
		                                                       h_pfp_track_nue_cc_coh,
		                                                       h_pfp_track_nue_cc_mec,
		                                                       h_pfp_track_nue_nc,
		                                                       h_pfp_track_numu_cc_qe,
		                                                       h_pfp_track_numu_cc_res,
		                                                       h_pfp_track_numu_cc_dis,
		                                                       h_pfp_track_numu_cc_coh,
		                                                       h_pfp_track_numu_cc_mec,
		                                                       h_pfp_track_numu_nc,
		                                                       h_pfp_track_nue_cc_mixed,
		                                                       h_pfp_track_numu_cc_mixed,
		                                                       h_pfp_track_cosmic,
		                                                       h_pfp_track_other_mixed,
		                                                       h_pfp_track_unmatched,
		                                                       h_pfp_shower_nue_cc_qe,
		                                                       h_pfp_shower_nue_cc_out_fv,
		                                                       h_pfp_shower_nue_cc_res,
		                                                       h_pfp_shower_nue_cc_dis,
		                                                       h_pfp_shower_nue_cc_coh,
		                                                       h_pfp_shower_nue_cc_mec,
		                                                       h_pfp_shower_nue_nc,
		                                                       h_pfp_shower_numu_cc_qe,
		                                                       h_pfp_shower_numu_cc_res,
		                                                       h_pfp_shower_numu_cc_dis,
		                                                       h_pfp_shower_numu_cc_coh,
		                                                       h_pfp_shower_numu_cc_mec,
		                                                       h_pfp_shower_numu_nc,
		                                                       h_pfp_shower_nue_cc_mixed,
		                                                       h_pfp_shower_numu_cc_mixed,
		                                                       h_pfp_shower_cosmic,
		                                                       h_pfp_shower_other_mixed,
		                                                       h_pfp_shower_unmatched);


		//in fv cut
		_functions_instance.selection_functions::fiducial_volume_cut(tpc_object_container_v, _x1, _x2, _y1, _y2, _z1, _z2, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, in_fv_counter_v);

		//vertex to flash cut
		_functions_instance.selection_functions::flashRecoVtxDist(largest_flash_v, tpc_object_container_v,
		                                                          tolerance, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, vtx_flash_counter_v);

		//distance between pfp shower and nue object cut
		_functions_instance.selection_functions::VtxNuDistance(tpc_object_container_v, shwr_nue_tolerance, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, shwr_tpco_counter_v);

		_functions_instance.selection_functions::VtxTrackNuDistance(tpc_object_container_v, trk_nue_tolerance, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, trk_tpco_counter_v);

		//hit threshold for showers cut
		_functions_instance.selection_functions::HitThreshold(tpc_object_container_v, shwr_hit_threshold, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, hit_threshold_counter_v);

		//open angle cut for the leading shower
		_functions_instance.selection_functions::OpenAngleCut(tpc_object_container_v, passed_tpco, tolerance_open_angle, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, open_angle_counter_v);

		//dEdx cut for the leading shower
		_functions_instance.selection_functions::dEdxCut(tpc_object_container_v, passed_tpco, tolerance_dedx_min, tolerance_dedx_max, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		_functions_instance.selection_functions::TotalOrigins(tabulated_origins, dedx_counter_v);

		//these are for the tefficiency plots, post all cuts
		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins.at(0) == 1)
		{
			if(_functions_instance.selection_functions::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
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
		_functions_instance.selection_functions::PostCutPlots(tpc_object_container_v, passed_tpco, _verbose,
		                                                      h_tracks_showers, h_tracks_showers_cosmic, h_tracks_showers_numu,
		                                                      h_leading_shower_open_angle_nue_cc, h_leading_shower_open_angle_nue_cc_mixed,
		                                                      h_leading_shower_open_angle_numu_cc, h_leading_shower_open_angle_numu_nc,
		                                                      h_leading_shower_open_angle_cosmic, h_leading_shower_open_angle_nue_nc,
		                                                      h_leading_shower_open_angle_numu_cc_mixed, h_leading_shower_open_angle_other_mixed,
		                                                      h_leading_shower_open_angle_unmatched,
		                                                      h_trk_vtx_dist_nue_cc, h_trk_vtx_dist_nue_cc_mixed,
		                                                      h_trk_vtx_dist_numu_cc, h_trk_vtx_dist_numu_nc,
		                                                      h_trk_vtx_dist_cosmic, h_trk_vtx_dist_nue_nc,
		                                                      h_trk_vtx_dist_numu_cc_mixed, h_trk_vtx_dist_other_mixed,
		                                                      h_trk_vtx_dist_unmatched);

		// _functions_instance.selection_functions::TopologyPlots(tpc_object_container_v, passed_tpco,
		//                                                        _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
		//                                                        h_pfp_track_shower_nue_cc_qe,
		//                                                        h_pfp_track_shower_nue_cc_out_fv,
		//                                                        h_pfp_track_shower_nue_cc_res,
		//                                                        h_pfp_track_shower_nue_cc_dis,
		//                                                        h_pfp_track_shower_nue_cc_coh,
		//                                                        h_pfp_track_shower_nue_cc_mec,
		//                                                        h_pfp_track_shower_nue_nc,
		//                                                        h_pfp_track_shower_numu_cc_qe,
		//                                                        h_pfp_track_shower_numu_cc_res,
		//                                                        h_pfp_track_shower_numu_cc_dis,
		//                                                        h_pfp_track_shower_numu_cc_coh,
		//                                                        h_pfp_track_shower_numu_cc_mec,
		//                                                        h_pfp_track_shower_numu_nc,
		//                                                        h_pfp_track_shower_nue_cc_mixed,
		//                                                        h_pfp_track_shower_numu_cc_mixed,
		//                                                        h_pfp_track_shower_cosmic,
		//                                                        h_pfp_track_shower_other_mixed,
		//                                                        h_pfp_track_shower_unmatched,
		//                                                        h_leading_shower_mc_pdg_nue_cc_qe,
		//                                                        h_leading_shower_mc_pdg_nue_cc_out_fv,
		//                                                        h_leading_shower_mc_pdg_nue_cc_res,
		//                                                        h_leading_shower_mc_pdg_nue_cc_dis,
		//                                                        h_leading_shower_mc_pdg_nue_cc_coh,
		//                                                        h_leading_shower_mc_pdg_nue_cc_mec,
		//                                                        h_leading_shower_mc_pdg_nue_nc,
		//                                                        h_leading_shower_mc_pdg_numu_cc_qe,
		//                                                        h_leading_shower_mc_pdg_numu_cc_res,
		//                                                        h_leading_shower_mc_pdg_numu_cc_dis,
		//                                                        h_leading_shower_mc_pdg_numu_cc_coh,
		//                                                        h_leading_shower_mc_pdg_numu_cc_mec,
		//                                                        h_leading_shower_mc_pdg_numu_nc,
		//                                                        h_leading_shower_mc_pdg_nue_cc_mixed,
		//                                                        h_leading_shower_mc_pdg_numu_cc_mixed,
		//                                                        h_leading_shower_mc_pdg_cosmic,
		//                                                        h_leading_shower_mc_pdg_other_mixed,
		//                                                        h_leading_shower_mc_pdg_unmatched
		//                                                        h_pfp_track_nue_cc_qe,
		//                                                        h_pfp_track_nue_cc_out_fv,
		//                                                        h_pfp_track_nue_cc_res,
		//                                                        h_pfp_track_nue_cc_dis,
		//                                                        h_pfp_track_nue_cc_coh,
		//                                                        h_pfp_track_nue_cc_mec,
		//                                                        h_pfp_track_nue_nc,
		//                                                        h_pfp_track_numu_cc_qe,
		//                                                        h_pfp_track_numu_cc_res,
		//                                                        h_pfp_track_numu_cc_dis,
		//                                                        h_pfp_track_numu_cc_coh,
		//                                                        h_pfp_track_numu_cc_mec,
		//                                                        h_pfp_track_numu_nc,
		//                                                        h_pfp_track_nue_cc_mixed,
		//                                                        h_pfp_track_numu_cc_mixed,
		//                                                        h_pfp_track_cosmic,
		//                                                        h_pfp_track_other_mixed,
		//                                                        h_pfp_track_unmatched,
		//                                                        h_pfp_shower_nue_cc_qe,
		//                                                        h_pfp_shower_nue_cc_out_fv,
		//                                                        h_pfp_shower_nue_cc_res,
		//                                                        h_pfp_shower_nue_cc_dis,
		//                                                        h_pfp_shower_nue_cc_coh,
		//                                                        h_pfp_shower_nue_cc_mec,
		//                                                        h_pfp_shower_nue_nc,
		//                                                        h_pfp_shower_numu_cc_qe,
		//                                                        h_pfp_shower_numu_cc_res,
		//                                                        h_pfp_shower_numu_cc_dis,
		//                                                        h_pfp_shower_numu_cc_coh,
		//                                                        h_pfp_shower_numu_cc_mec,
		//                                                        h_pfp_shower_numu_nc,
		//                                                        h_pfp_shower_nue_cc_mixed,
		//                                                        h_pfp_shower_numu_cc_mixed,
		//                                                        h_pfp_shower_cosmic,
		//                                                        h_pfp_shower_other_mixed,
		//                                                        h_pfp_shower_unmatched);

	}//end event loop

	std::cout << "------------------ " << std::endl;
	std::cout << " MC Entries in FV: " << total_mc_entries_inFV << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "  End Selection    " << std::endl;
	std::cout << "------------------ " << std::endl;

	//change mc_nue_cc_counter to total_mc_entries_inFV once files are ready!

	//we also want some metrics to print at the end
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    reco_nue_counter_v,
	                                                    "Reco Nue"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    in_fv_counter_v,
	                                                    "In FV"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    vtx_flash_counter_v,
	                                                    "Vtx-to-Flash"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    shwr_tpco_counter_v,
	                                                    "Shower-to-TPCO"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    trk_tpco_counter_v,
	                                                    "Track-to-TPCO"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    hit_threshold_counter_v,
	                                                    "Hit Threshold"
	                                                    );

	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    open_angle_counter_v,
	                                                    "Open Angle"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    dedx_counter_v,
	                                                    " dE / dx "
	                                                    );

	std::vector<double> * xsec_cc = new std::vector<double>;
	const double final_counter                = dedx_counter_v->at(7);
	const double final_counter_nue_cc         = dedx_counter_v->at(0);
	const double final_counter_nue_cc_mixed   = dedx_counter_v->at(1);
	const double final_counter_nue_cc_out_fv  = dedx_counter_v->at(9);
	const double final_counter_cosmic         = dedx_counter_v->at(2);
	const double final_counter_nue_nc         = dedx_counter_v->at(3);
	const double final_counter_numu_cc        = dedx_counter_v->at(4);
	const double final_counter_numu_cc_mixed  = dedx_counter_v->at(11);
	const double final_counter_numu_nc        = dedx_counter_v->at(10);
	const double final_counter_unmatched      = dedx_counter_v->at(5);
	const double final_counter_other_mixed    = dedx_counter_v->at(6);
	const int n_total = final_counter;
	const int n_bkg = (final_counter_nue_cc_mixed + final_counter_nue_cc_out_fv + final_counter_cosmic + final_counter_nue_nc
	                   + final_counter_numu_cc + final_counter_numu_cc_mixed + final_counter_numu_nc + final_counter_unmatched + final_counter_other_mixed);
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
	h_leading_shower_open_angle_numu_nc->SetStats(kFALSE);
	h_leading_shower_open_angle_cosmic->SetStats(kFALSE);
	h_leading_shower_open_angle_nue_nc->SetStats(kFALSE);
	h_leading_shower_open_angle_numu_cc_mixed->SetStats(kFALSE);
	h_leading_shower_open_angle_other_mixed->SetStats(kFALSE);
	h_leading_shower_open_angle_unmatched->SetStats(kFALSE);
	h_leading_shower_open_angle_nue_cc->SetFillColor(30);
	h_leading_shower_open_angle_nue_cc_mixed->SetFillColor(38);
	h_leading_shower_open_angle_numu_cc->SetFillColor(28);
	h_leading_shower_open_angle_numu_nc->SetFillColor(36);
	h_leading_shower_open_angle_cosmic->SetFillColor(39);
	h_leading_shower_open_angle_nue_nc->SetFillColor(46);
	h_leading_shower_open_angle_numu_cc_mixed->SetFillColor(25);
	h_leading_shower_open_angle_other_mixed->SetFillColor(42);
	h_leading_shower_open_angle_unmatched->SetFillColor(12);
	open_angle_stack->Add(h_leading_shower_open_angle_unmatched);
	open_angle_stack->Add(h_leading_shower_open_angle_other_mixed);
	open_angle_stack->Add(h_leading_shower_open_angle_numu_nc);
	open_angle_stack->Add(h_leading_shower_open_angle_nue_nc);
	open_angle_stack->Add(h_leading_shower_open_angle_numu_cc_mixed);
	open_angle_stack->Add(h_leading_shower_open_angle_numu_cc);
	open_angle_stack->Add(h_leading_shower_open_angle_cosmic);
	open_angle_stack->Add(h_leading_shower_open_angle_nue_cc_mixed);
	open_angle_stack->Add(h_leading_shower_open_angle_nue_cc);
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
	leg_stack->AddEntry(h_leading_shower_open_angle_nue_nc,          "Nue NC", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_numu_nc,         "Numu NC", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_other_mixed,     "Other Mixed", "f");
	leg_stack->AddEntry(h_leading_shower_open_angle_unmatched,       "Unmatched", "f");
	leg_stack->Draw();
	open_angle_stack_c1->Print("post_cuts_leading_shower_open_angle.pdf");

	TCanvas * trk_vtx_dist_stack_c1 = new TCanvas();
	trk_vtx_dist_stack_c1->cd();
	THStack * trk_vtx_dist_stack = new THStack();
	h_trk_vtx_dist_nue_cc->SetStats(kFALSE);
	h_trk_vtx_dist_nue_cc_mixed->SetStats(kFALSE);
	h_trk_vtx_dist_numu_cc->SetStats(kFALSE);
	h_trk_vtx_dist_numu_nc->SetStats(kFALSE);
	h_trk_vtx_dist_cosmic->SetStats(kFALSE);
	h_trk_vtx_dist_nue_nc->SetStats(kFALSE);
	h_trk_vtx_dist_numu_cc_mixed->SetStats(kFALSE);
	h_trk_vtx_dist_other_mixed->SetStats(kFALSE);
	h_trk_vtx_dist_unmatched->SetStats(kFALSE);
	h_trk_vtx_dist_nue_cc->SetFillColor(30);
	h_trk_vtx_dist_nue_cc_mixed->SetFillColor(38);
	h_trk_vtx_dist_numu_cc->SetFillColor(28);
	h_trk_vtx_dist_numu_nc->SetFillColor(36);
	h_trk_vtx_dist_cosmic->SetFillColor(39);
	h_trk_vtx_dist_nue_nc->SetFillColor(46);
	h_trk_vtx_dist_numu_cc_mixed->SetFillColor(25);
	h_trk_vtx_dist_other_mixed->SetFillColor(42);
	h_trk_vtx_dist_unmatched->SetFillColor(12);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_unmatched);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_other_mixed);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_numu_nc);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_nue_nc);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_numu_cc_mixed);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_numu_cc);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_cosmic);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_nue_cc_mixed);
	trk_vtx_dist_stack->Add(h_trk_vtx_dist_nue_cc);
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
	leg_stack2->AddEntry(h_trk_vtx_dist_nue_nc,          "Nue NC", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_numu_nc,         "Numu NC", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_other_mixed,     "Other Mixed", "f");
	leg_stack2->AddEntry(h_trk_vtx_dist_unmatched,       "Unmatched", "f");
	leg_stack2->Draw();
	trk_vtx_dist_stack_c1->Print("post_cuts_track_to_vtx.pdf");

	TCanvas * track_stack_c1 = new TCanvas();
	track_stack_c1->cd();
	THStack * h_track_stack = new THStack();
	h_pfp_track_nue_cc_qe->SetStats(kFALSE);
	h_pfp_track_nue_cc_out_fv->SetStats(kFALSE);
	h_pfp_track_nue_cc_res->SetStats(kFALSE);
	h_pfp_track_nue_cc_dis->SetStats(kFALSE);
	h_pfp_track_nue_cc_coh->SetStats(kFALSE);
	h_pfp_track_nue_cc_mec->SetStats(kFALSE);
	h_pfp_track_nue_nc->SetStats(kFALSE);
	h_pfp_track_numu_cc_qe->SetStats(kFALSE);
	h_pfp_track_numu_cc_res->SetStats(kFALSE);
	h_pfp_track_numu_cc_dis->SetStats(kFALSE);
	h_pfp_track_numu_cc_coh->SetStats(kFALSE);
	h_pfp_track_numu_cc_mec->SetStats(kFALSE);
	h_pfp_track_numu_nc->SetStats(kFALSE);
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
	h_pfp_track_nue_nc->SetFillColor(46);
	h_pfp_track_numu_cc_qe->SetFillColor(28);
	h_pfp_track_numu_cc_res->SetFillColor(27);
	h_pfp_track_numu_cc_dis->SetFillColor(26);
	h_pfp_track_numu_cc_coh->SetFillColor(23);
	h_pfp_track_numu_cc_mec->SetFillColor(22);
	h_pfp_track_numu_nc->SetFillColor(36);
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
	h_track_stack->Add(h_pfp_track_nue_nc       );
	h_track_stack->Add(h_pfp_track_numu_cc_qe   );
	h_track_stack->Add(h_pfp_track_numu_cc_res  );
	h_track_stack->Add(h_pfp_track_numu_cc_dis  );
	h_track_stack->Add(h_pfp_track_numu_cc_coh  );
	h_track_stack->Add(h_pfp_track_numu_cc_mec  );
	h_track_stack->Add(h_pfp_track_numu_nc      );
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
	leg_track_stack_l1->AddEntry(h_pfp_track_nue_nc,         "Nue NC", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_qe,     "Numu CC QE", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_res,    "Numu CC Res", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_dis,    "Numu CC DIS", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_coh,    "Numu CC Coh", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_mec,    "Numu CC MEC", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_nc,        "Numu NC", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_nue_cc_mixed,   "Nue CC Mixed", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_numu_cc_mixed,  "Numu CC Mixed", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_cosmic,         "Cosmic", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_other_mixed,    "Other Mixed", "f");
	leg_track_stack_l1->AddEntry(h_pfp_track_unmatched,      "Unmatched", "f");
	leg_track_stack_l1->Draw();
	track_stack_c1->Print("selected_pfp_track_stack.pdf");

	TCanvas * shower_stack_c1 = new TCanvas();
	shower_stack_c1->cd();
	THStack * h_shower_stack = new THStack();
	h_pfp_shower_nue_cc_qe->SetStats(kFALSE);
	h_pfp_shower_nue_cc_out_fv->SetStats(kFALSE);
	h_pfp_shower_nue_cc_res->SetStats(kFALSE);
	h_pfp_shower_nue_cc_dis->SetStats(kFALSE);
	h_pfp_shower_nue_cc_coh->SetStats(kFALSE);
	h_pfp_shower_nue_cc_mec->SetStats(kFALSE);
	h_pfp_shower_nue_nc->SetStats(kFALSE);
	h_pfp_shower_numu_cc_qe->SetStats(kFALSE);
	h_pfp_shower_numu_cc_res->SetStats(kFALSE);
	h_pfp_shower_numu_cc_dis->SetStats(kFALSE);
	h_pfp_shower_numu_cc_coh->SetStats(kFALSE);
	h_pfp_shower_numu_cc_mec->SetStats(kFALSE);
	h_pfp_shower_numu_nc->SetStats(kFALSE);
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
	h_pfp_shower_nue_nc->SetFillColor(46);
	h_pfp_shower_numu_cc_qe->SetFillColor(28);
	h_pfp_shower_numu_cc_res->SetFillColor(27);
	h_pfp_shower_numu_cc_dis->SetFillColor(26);
	h_pfp_shower_numu_cc_coh->SetFillColor(23);
	h_pfp_shower_numu_cc_mec->SetFillColor(22);
	h_pfp_shower_numu_nc->SetFillColor(36);
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
	h_shower_stack->Add(h_pfp_shower_nue_nc       );
	h_shower_stack->Add(h_pfp_shower_numu_cc_qe   );
	h_shower_stack->Add(h_pfp_shower_numu_cc_res  );
	h_shower_stack->Add(h_pfp_shower_numu_cc_dis  );
	h_shower_stack->Add(h_pfp_shower_numu_cc_coh  );
	h_shower_stack->Add(h_pfp_shower_numu_cc_mec  );
	h_shower_stack->Add(h_pfp_shower_numu_nc      );
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
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nue_nc,         "Nue NC", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_qe,     "Numu CC QE", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_res,    "Numu CC Res", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_dis,    "Numu CC DIS", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_coh,    "Numu CC Coh", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_mec,    "Numu CC MEC", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_nc,        "Numu NC", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_nue_cc_mixed,   "Nue CC Mixed", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_numu_cc_mixed,  "Numu CC Mixed", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_cosmic,         "Cosmic", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_other_mixed,    "Other Mixed", "f");
	leg_shower_stack_l1->AddEntry(h_pfp_shower_unmatched,      "Unmatched", "f");
	leg_shower_stack_l1->Draw();
	shower_stack_c1->Print("selected_pfp_shower_stack.pdf");

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
	h_pfp_track_shower_nue_nc->Draw("colz");
	h_pfp_track_shower_nue_nc->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_nue_nc->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c7->Print("selected_pfp_track_shower_nue_nc.pdf");
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
	h_pfp_track_shower_numu_nc->Draw("colz");
	h_pfp_track_shower_numu_nc->GetXaxis()->SetTitle("PFP Tracks ");
	h_pfp_track_shower_numu_nc->GetYaxis()->SetTitle("PFP Showers");
	track_shower_c13->Print("selected_pfp_track_shower_numu_nc.pdf");
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
		h_leading_shower_mc_pdg_nue_nc->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_leading_shower_mc_pdg_numu_nc->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
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
		h_leading_shower_mc_pdg_nue_nc->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_numu_nc->GetYaxis()->SetBinLabel(i,str_mc_particle[i-1]);
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
	h_leading_shower_mc_pdg_nue_nc->Draw("colz");
	h_leading_shower_mc_pdg_nue_nc->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_nue_nc->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_nue_nc->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_nue_nc->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c7->Print("selected_leading_shower_mc_pdg_nue_nc.pdf");
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
	h_leading_shower_mc_pdg_numu_nc->Draw("colz");
	h_leading_shower_mc_pdg_numu_nc->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_mc_pdg_numu_nc->GetYaxis()->SetTitleOffset(1.35);
	h_leading_shower_mc_pdg_numu_nc->GetXaxis()->SetTitle("Leading Shower Origin ");
	h_leading_shower_mc_pdg_numu_nc->GetYaxis()->SetTitle("Leading Shower True Particle");
	leading_c13->Print("selected_leading_shower_mc_pdg_numu_nc.pdf");
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
