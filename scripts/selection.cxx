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

	//TEfficiency histograms
	TH1D * h_nue_eng_eff_den           = new TH1D("h_nue_eng_eff_den", "h_nue_eng_eff_den", 8, 0, 4);
	TH1D * h_nue_eng_eff_num           = new TH1D("h_nue_eng_eff_num", "h_nue_eng_eff_num", 8, 0, 4);
	TH1D * h_nue_vtx_x_eff_den         = new TH1D("h_nue_vtx_x_eff_den", "h_nue_vtx_x_eff_den", 10, 0, 256.35);
	TH1D * h_nue_vtx_x_eff_num         = new TH1D("h_nue_vtx_x_eff_num", "h_nue_vtx_x_eff_num", 10, 0, 256.35);
	TH1D * h_nue_vtx_y_eff_den         = new TH1D("h_nue_vtx_y_eff_den", "h_nue_vtx_y_eff_den", 10, -116.5, 116.5);
	TH1D * h_nue_vtx_y_eff_num         = new TH1D("h_nue_vtx_y_eff_num", "h_nue_vtx_y_eff_num", 10, -116.5, 116.5);
	TH1D * h_nue_vtx_z_eff_den         = new TH1D("h_nue_vtx_z_eff_den", "h_nue_vtx_z_eff_den", 10, 0, 1036.8);
	TH1D * h_nue_vtx_z_eff_num         = new TH1D("h_nue_vtx_z_eff_num", "h_nue_vtx_z_eff_num", 10, 0, 1036.8);
	TH1D * h_nue_dir_x_eff_den         = new TH1D("h_nue_dir_x_eff_den", "h_nue_dir_x_eff_den", 10, -1, 1);
	TH1D * h_nue_dir_x_eff_num         = new TH1D("h_nue_dir_x_eff_num", "h_nue_dir_x_eff_num", 10, -1, 1);
	TH1D * h_nue_dir_y_eff_den         = new TH1D("h_nue_dir_y_eff_den", "h_nue_dir_y_eff_den", 10, -1, 1);
	TH1D * h_nue_dir_y_eff_num         = new TH1D("h_nue_dir_y_eff_num", "h_nue_dir_y_eff_num", 10, -1, 1);
	TH1D * h_nue_dir_z_eff_den         = new TH1D("h_nue_dir_z_eff_den", "h_nue_dir_z_eff_den", 10, -1, 1);
	TH1D * h_nue_dir_z_eff_num         = new TH1D("h_nue_dir_z_eff_num", "h_nue_dir_z_eff_num", 10, -1, 1);
	TH1D * h_nue_num_part_eff_den      = new TH1D("h_nue_num_part_eff_den", "h_nue_num_part_eff_den", 10, 0, 20);
	TH1D * h_nue_num_part_eff_num      = new TH1D("h_nue_num_part_eff_num", "h_nue_num_part_eff_num", 10, 0, 20);
	TH1D * h_nue_num_chrg_part_eff_den = new TH1D("h_nue_num_chrg_part_eff_den", "h_nue_num_chrg_part_eff_den", 10, 0, 20);
	TH1D * h_nue_num_chrg_part_eff_num = new TH1D("h_nue_num_chrg_part_eff_num", "h_nue_num_chrg_part_eff_num", 10, 0, 20);
	TH1D * h_nue_cos_theta_eff_den     = new TH1D("h_nue_cos_theta_eff_den", "h_nue_cos_theta_eff_den", 10, -1, 1);
	TH1D * h_nue_cos_theta_eff_num     = new TH1D("h_nue_cos_theta_eff_num", "h_nue_cos_theta_eff_num", 10, -1, 1);
	TH1D * h_nue_phi_eff_den           = new TH1D("h_nue_phi_eff_den", "h_nue_phi_eff_den", 10, -180, 180);
	TH1D * h_nue_phi_eff_num           = new TH1D("h_nue_phi_eff_num", "h_nue_phi_eff_num", 10, -180, 180);

	TH2I * h_tracks_showers         = new TH2I("h_tracks_showers", "h_tracks_showers", 8, 0, 8, 8, 0, 8);
	TH2I * h_tracks_showers_cosmic  = new TH2I("h_tracks_showers_cosmic", "h_tracks_showers_cosmic", 8, 0, 8, 8, 0, 8);
	TH2I * h_tracks_showers_numu    = new TH2I("h_tracks_showers_numu", "h_tracks_showers_numu", 8, 0, 8, 8, 0, 8);

	TH1D * h_leading_shower_open_angle_nue_cc        = new TH1D("h_leading_shower_open_angle_nue_cc", "h_leading_shower_open_angle_nue_cc", 25, 0, 50);
	TH1D * h_leading_shower_open_angle_nue_cc_mixed  = new TH1D("h_leading_shower_open_angle_nue_cc_mixed", "h_leading_shower_open_angle_nue_cc_mixed", 25, 0, 50);
	TH1D * h_leading_shower_open_angle_numu_cc       = new TH1D("h_leading_shower_open_angle_numu_cc", "h_leading_shower_open_angle_numu_cc", 25, 0, 50);
	TH1D * h_leading_shower_open_angle_numu_nc       = new TH1D("h_leading_shower_open_angle_numu_nc", "h_leading_shower_open_angle_numu_nc", 25, 0, 50);
	TH1D * h_leading_shower_open_angle_cosmic        = new TH1D("h_leading_shower_open_angle_cosmic", "h_leading_shower_open_angle_cosmic", 25, 0, 50);
	TH1D * h_leading_shower_open_angle_nue_nc        = new TH1D("h_leading_shower_open_angle_nue_nc", "h_leading_shower_open_angle_nue_nc", 25, 0, 50);
	TH1D * h_leading_shower_open_angle_numu_cc_mixed = new TH1D("h_leading_shower_open_angle_numu_cc_mixed", "h_leading_shower_open_angle_numu_cc_mixed", 25, 0, 50);
	TH1D * h_leading_shower_open_angle_other_mixed   = new TH1D("h_leading_shower_open_angle_other_mixed", "h_leading_shower_open_angle_other_mixed", 25, 0, 50);
	TH1D * h_leading_shower_open_angle_unmatched     = new TH1D("h_leading_shower_open_angle_unmatched", "h_leading_shower_open_angle_unmatched", 25, 0, 50);

	TH1D * h_trk_vtx_dist_nue_cc        = new TH1D("h_trk_vtx_dist_nue_cc", "h_trk_vtx_dist_nue_cc", 25, 0, 20);
	TH1D * h_trk_vtx_dist_nue_cc_mixed  = new TH1D("h_trk_vtx_dist_nue_cc_mixed", "h_trk_vtx_dist_nue_cc_mixed", 25, 0, 20);
	TH1D * h_trk_vtx_dist_numu_cc       = new TH1D("h_trk_vtx_dist_numu_cc", "h_trk_vtx_dist_numu_cc", 25, 0, 20);
	TH1D * h_trk_vtx_dist_numu_nc       = new TH1D("h_trk_vtx_dist_numu_nc", "h_trk_vtx_dist_numu_nc", 25, 0, 20);
	TH1D * h_trk_vtx_dist_cosmic        = new TH1D("h_trk_vtx_dist_cosmic", "h_trk_vtx_dist_cosmic", 25, 0, 20);
	TH1D * h_trk_vtx_dist_nue_nc        = new TH1D("h_trk_vtx_dist_nue_nc", "h_trk_vtx_dist_nue_nc", 25, 0, 20);
	TH1D * h_trk_vtx_dist_numu_cc_mixed = new TH1D("h_trk_vtx_dist_numu_cc_mixed", "h_trk_vtx_dist_numu_cc_mixed", 25, 0, 20);
	TH1D * h_trk_vtx_dist_other_mixed   = new TH1D("h_trk_vtx_dist_other_mixed", "h_trk_vtx_dist_other_mixed", 25, 0, 20);
	TH1D * h_trk_vtx_dist_unmatched     = new TH1D("h_trk_vtx_dist_unmatched", "h_trk_vtx_dist_unmatched", 25, 0, 20);

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
	mctruth_counter_tree->SetBranchAddress("fMCNuEnegy", &mc_nu_energy);
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
		if(mc_nu_id == 1 || mc_nu_id == 5)
		//if this event is a true nue CC interaction and is inside the FV
		//also include nue_cc_bar as in the tpco classification I use the nue-bar as well
		{
			if(_functions_instance.selection_functions::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
			                                                  _x1, _x2, _y1,
			                                                  _y2, _z1, _z2) == true)
			{
				h_nue_eng_eff_den->Fill(mc_nu_energy);
				h_nue_vtx_x_eff_den->Fill(mc_nu_vtx_x);
				h_nue_vtx_y_eff_den->Fill(mc_nu_vtx_y);
				h_nue_vtx_z_eff_den->Fill(mc_nu_vtx_z);
				h_nue_dir_x_eff_den->Fill(mc_nu_dir_x);
				h_nue_dir_y_eff_den->Fill(mc_nu_dir_y);
				h_nue_dir_z_eff_den->Fill(mc_nu_dir_z);
				h_nue_num_part_eff_den->Fill(mc_nu_num_particles);
				h_nue_num_chrg_part_eff_den->Fill(mc_nu_num_charged_particles);
				h_nue_cos_theta_eff_den->Fill(mc_cos_theta);
				h_nue_phi_eff_den->Fill(mc_phi);
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
		std::vector<int> * passed_tpco = new std::vector<int>;
		passed_tpco->resize(tpc_object_container_v->size(), 1);

		//** start the cuts here **

		//reco nue cut
		_functions_instance.selection_functions::HasNue(tpc_object_container_v, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		reco_nue_counter_nue_cc        = reco_nue_counter_nue_cc + tabulated_origins.at(0);
		reco_nue_counter_nue_cc_mixed  = reco_nue_counter_nue_cc_mixed + tabulated_origins.at(1);
		reco_nue_counter_nue_cc_out_fv = reco_nue_counter_nue_cc_out_fv + tabulated_origins.at(9);
		reco_nue_counter_cosmic        = reco_nue_counter_cosmic + tabulated_origins.at(2);
		reco_nue_counter_nue_nc        = reco_nue_counter_nue_nc + tabulated_origins.at(3);
		reco_nue_counter_numu_cc       = reco_nue_counter_numu_cc + tabulated_origins.at(4);
		reco_nue_counter_numu_cc_mixed = reco_nue_counter_numu_cc_mixed + tabulated_origins.at(11);
		reco_nue_counter_numu_nc       = reco_nue_counter_numu_nc + tabulated_origins.at(10);
		reco_nue_counter_unmatched     = reco_nue_counter_unmatched + tabulated_origins.at(5);
		reco_nue_counter_other_mixed   = reco_nue_counter_other_mixed + tabulated_origins.at(6);
		reco_nue_counter = reco_nue_counter + tabulated_origins.at(7);

		//in fv cut
		_functions_instance.selection_functions::fiducial_volume_cut(tpc_object_container_v, _x1, _x2, _y1, _y2, _z1, _z2, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		in_fv_counter_nue_cc        = in_fv_counter_nue_cc + tabulated_origins.at(0);
		in_fv_counter_nue_cc_mixed  = in_fv_counter_nue_cc_mixed + tabulated_origins.at(1);
		in_fv_counter_nue_cc_out_fv = in_fv_counter_nue_cc_out_fv + tabulated_origins.at(9);
		in_fv_counter_cosmic        = in_fv_counter_cosmic + tabulated_origins.at(2);
		in_fv_counter_nue_nc        = in_fv_counter_nue_nc + tabulated_origins.at(3);
		in_fv_counter_numu_cc       = in_fv_counter_numu_cc + tabulated_origins.at(4);
		in_fv_counter_numu_cc_mixed = in_fv_counter_numu_cc_mixed + tabulated_origins.at(11);
		in_fv_counter_numu_nc       = in_fv_counter_numu_nc + tabulated_origins.at(10);
		in_fv_counter_unmatched     = in_fv_counter_unmatched + tabulated_origins.at(5);
		in_fv_counter_other_mixed   = in_fv_counter_other_mixed + tabulated_origins.at(6);
		in_fv_counter = in_fv_counter + tabulated_origins.at(7);

		//vertex to flash cut
		_functions_instance.selection_functions::flashRecoVtxDist(largest_flash_v, tpc_object_container_v,
		                                                          tolerance, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		vtx_flash_counter_nue_cc        = vtx_flash_counter_nue_cc + tabulated_origins.at(0);
		vtx_flash_counter_nue_cc_mixed  = vtx_flash_counter_nue_cc_mixed + tabulated_origins.at(1);
		vtx_flash_counter_nue_cc_out_fv = vtx_flash_counter_nue_cc_out_fv + tabulated_origins.at(9);
		vtx_flash_counter_cosmic        = vtx_flash_counter_cosmic + tabulated_origins.at(2);
		vtx_flash_counter_nue_nc        = vtx_flash_counter_nue_nc + tabulated_origins.at(3);
		vtx_flash_counter_numu_cc       = vtx_flash_counter_numu_cc + tabulated_origins.at(4);
		vtx_flash_counter_numu_cc_mixed = vtx_flash_counter_numu_cc_mixed + tabulated_origins.at(11);
		vtx_flash_counter_numu_nc       = vtx_flash_counter_numu_nc + tabulated_origins.at(10);
		vtx_flash_counter_unmatched     = vtx_flash_counter_unmatched + tabulated_origins.at(5);
		vtx_flash_counter_other_mixed   = vtx_flash_counter_other_mixed + tabulated_origins.at(6);
		vtx_flash_counter = vtx_flash_counter + tabulated_origins.at(7);

		//distance between pfp shower and nue object cut
		_functions_instance.selection_functions::VtxNuDistance(tpc_object_container_v, shwr_nue_tolerance, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		shwr_tpco_counter_nue_cc        = shwr_tpco_counter_nue_cc + tabulated_origins.at(0);
		shwr_tpco_counter_nue_cc_mixed  = shwr_tpco_counter_nue_cc_mixed + tabulated_origins.at(1);
		shwr_tpco_counter_nue_cc_out_fv = shwr_tpco_counter_nue_cc_out_fv + tabulated_origins.at(9);
		shwr_tpco_counter_cosmic        = shwr_tpco_counter_cosmic + tabulated_origins.at(2);
		shwr_tpco_counter_nue_nc        = shwr_tpco_counter_nue_nc + tabulated_origins.at(3);
		shwr_tpco_counter_numu_cc       = shwr_tpco_counter_numu_cc + tabulated_origins.at(4);
		shwr_tpco_counter_numu_cc_mixed = shwr_tpco_counter_numu_cc_mixed + tabulated_origins.at(11);
		shwr_tpco_counter_numu_nc       = shwr_tpco_counter_numu_nc + tabulated_origins.at(10);
		shwr_tpco_counter_unmatched     = shwr_tpco_counter_unmatched + tabulated_origins.at(5);
		shwr_tpco_counter_other_mixed   = shwr_tpco_counter_other_mixed + tabulated_origins.at(6);
		shwr_tpco_counter = shwr_tpco_counter + tabulated_origins.at(7);

		_functions_instance.selection_functions::VtxTrackNuDistance(tpc_object_container_v, trk_nue_tolerance, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		trk_tpco_counter_nue_cc        = trk_tpco_counter_nue_cc + tabulated_origins.at(0);
		trk_tpco_counter_nue_cc_mixed  = trk_tpco_counter_nue_cc_mixed + tabulated_origins.at(1);
		trk_tpco_counter_nue_cc_out_fv = trk_tpco_counter_nue_cc_out_fv + tabulated_origins.at(9);
		trk_tpco_counter_cosmic        = trk_tpco_counter_cosmic + tabulated_origins.at(2);
		trk_tpco_counter_nue_nc        = trk_tpco_counter_nue_nc + tabulated_origins.at(3);
		trk_tpco_counter_numu_cc       = trk_tpco_counter_numu_cc + tabulated_origins.at(4);
		trk_tpco_counter_numu_cc_mixed = trk_tpco_counter_numu_cc_mixed + tabulated_origins.at(11);
		trk_tpco_counter_numu_nc       = trk_tpco_counter_numu_nc + tabulated_origins.at(10);
		trk_tpco_counter_unmatched     = trk_tpco_counter_unmatched + tabulated_origins.at(5);
		trk_tpco_counter_other_mixed   = trk_tpco_counter_other_mixed + tabulated_origins.at(6);
		trk_tpco_counter               = trk_tpco_counter + tabulated_origins.at(7);

		//hit threshold for showers cut
		_functions_instance.selection_functions::HitThreshold(tpc_object_container_v, shwr_hit_threshold, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		hit_threshold_counter_nue_cc        = hit_threshold_counter_nue_cc + tabulated_origins.at(0);
		hit_threshold_counter_nue_cc_mixed  = hit_threshold_counter_nue_cc_mixed + tabulated_origins.at(1);
		hit_threshold_counter_nue_cc_out_fv = hit_threshold_counter_nue_cc_out_fv + tabulated_origins.at(9);
		hit_threshold_counter_cosmic        = hit_threshold_counter_cosmic + tabulated_origins.at(2);
		hit_threshold_counter_nue_nc        = hit_threshold_counter_nue_nc + tabulated_origins.at(3);
		hit_threshold_counter_numu_cc       = hit_threshold_counter_numu_cc + tabulated_origins.at(4);
		hit_threshold_counter_numu_cc_mixed = hit_threshold_counter_numu_cc_mixed + tabulated_origins.at(11);
		hit_threshold_counter_numu_nc       = hit_threshold_counter_numu_nc + tabulated_origins.at(10);
		hit_threshold_counter_unmatched     = hit_threshold_counter_unmatched + tabulated_origins.at(5);
		hit_threshold_counter_other_mixed   = hit_threshold_counter_other_mixed + tabulated_origins.at(6);
		hit_threshold_counter = hit_threshold_counter + tabulated_origins.at(7);

		//open angle cut for the leading shower
		_functions_instance.selection_functions::OpenAngleCut(tpc_object_container_v, passed_tpco, tolerance_open_angle, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco,
		                                                                             _x1, _x2, _y1, _y2, _z1, _z2, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z);
		open_angle_counter_nue_cc        = open_angle_counter_nue_cc + tabulated_origins.at(0);
		open_angle_counter_nue_cc_mixed  = open_angle_counter_nue_cc_mixed + tabulated_origins.at(1);
		open_angle_counter_nue_cc_out_fv = open_angle_counter_nue_cc_out_fv + tabulated_origins.at(9);
		open_angle_counter_cosmic        = open_angle_counter_cosmic + tabulated_origins.at(2);
		open_angle_counter_nue_nc        = open_angle_counter_nue_nc + tabulated_origins.at(3);
		open_angle_counter_numu_cc       = open_angle_counter_numu_cc + tabulated_origins.at(4);
		open_angle_counter_numu_cc_mixed = open_angle_counter_numu_cc_mixed + tabulated_origins.at(11);
		open_angle_counter_numu_nc       = open_angle_counter_numu_nc + tabulated_origins.at(10);
		open_angle_counter_unmatched     = open_angle_counter_unmatched + tabulated_origins.at(5);
		open_angle_counter_other_mixed   = open_angle_counter_other_mixed + tabulated_origins.at(6);
		open_angle_counter               = open_angle_counter + tabulated_origins.at(7);


		if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins.at(0) == 1)
		{
			if(_functions_instance.selection_functions::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
			                                                  _x1, _x2, _y1,
			                                                  _y2, _z1, _z2) == true)
			{
				selected_energy_vector.push_back(mc_nu_energy);
				h_nue_eng_eff_num->Fill(mc_nu_energy);
				h_nue_vtx_x_eff_num->Fill(mc_nu_vtx_x);
				h_nue_vtx_y_eff_num->Fill(mc_nu_vtx_y);
				h_nue_vtx_z_eff_num->Fill(mc_nu_vtx_z);
				h_nue_dir_x_eff_num->Fill(mc_nu_dir_x);
				h_nue_dir_y_eff_num->Fill(mc_nu_dir_y);
				h_nue_dir_z_eff_num->Fill(mc_nu_dir_z);
				h_nue_num_part_eff_num->Fill(mc_nu_num_particles);
				h_nue_num_chrg_part_eff_num->Fill(mc_nu_num_charged_particles);
				h_nue_cos_theta_eff_num->Fill(mc_cos_theta);
				h_nue_phi_eff_num->Fill(mc_phi);
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
	                                                    reco_nue_counter,
	                                                    reco_nue_counter_nue_cc,
	                                                    reco_nue_counter_nue_cc_mixed,
	                                                    reco_nue_counter_nue_cc_out_fv,
	                                                    reco_nue_counter_cosmic,
	                                                    reco_nue_counter_nue_nc,
	                                                    reco_nue_counter_numu_cc,
	                                                    reco_nue_counter_numu_cc_mixed,
	                                                    reco_nue_counter_numu_nc,
	                                                    reco_nue_counter_unmatched,
	                                                    reco_nue_counter_other_mixed,
	                                                    "Reco Nue"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    in_fv_counter,
	                                                    in_fv_counter_nue_cc,
	                                                    in_fv_counter_nue_cc_mixed,
	                                                    in_fv_counter_nue_cc_out_fv,
	                                                    in_fv_counter_cosmic,
	                                                    in_fv_counter_nue_nc,
	                                                    in_fv_counter_numu_cc,
	                                                    in_fv_counter_numu_cc_mixed,
	                                                    in_fv_counter_numu_nc,
	                                                    in_fv_counter_unmatched,
	                                                    in_fv_counter_other_mixed,
	                                                    "In FV"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    vtx_flash_counter,
	                                                    vtx_flash_counter_nue_cc,
	                                                    vtx_flash_counter_nue_cc_mixed,
	                                                    vtx_flash_counter_nue_cc_out_fv,
	                                                    vtx_flash_counter_cosmic,
	                                                    vtx_flash_counter_nue_nc,
	                                                    vtx_flash_counter_numu_cc,
	                                                    vtx_flash_counter_numu_cc_mixed,
	                                                    vtx_flash_counter_numu_nc,
	                                                    vtx_flash_counter_unmatched,
	                                                    vtx_flash_counter_other_mixed,
	                                                    "Vtx-to-Flash"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    shwr_tpco_counter,
	                                                    shwr_tpco_counter_nue_cc,
	                                                    shwr_tpco_counter_nue_cc_mixed,
	                                                    shwr_tpco_counter_nue_cc_out_fv,
	                                                    shwr_tpco_counter_cosmic,
	                                                    shwr_tpco_counter_nue_nc,
	                                                    shwr_tpco_counter_numu_cc,
	                                                    shwr_tpco_counter_numu_cc_mixed,
	                                                    shwr_tpco_counter_numu_nc,
	                                                    shwr_tpco_counter_unmatched,
	                                                    shwr_tpco_counter_other_mixed,
	                                                    "Shower-to-TPCO"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    trk_tpco_counter,
	                                                    trk_tpco_counter_nue_cc,
	                                                    trk_tpco_counter_nue_cc_mixed,
	                                                    trk_tpco_counter_nue_cc_out_fv,
	                                                    trk_tpco_counter_cosmic,
	                                                    trk_tpco_counter_nue_nc,
	                                                    trk_tpco_counter_numu_cc,
	                                                    trk_tpco_counter_numu_cc_mixed,
	                                                    trk_tpco_counter_numu_nc,
	                                                    trk_tpco_counter_unmatched,
	                                                    trk_tpco_counter_other_mixed,
	                                                    "Track-to-TPCO"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    hit_threshold_counter,
	                                                    hit_threshold_counter_nue_cc,
	                                                    hit_threshold_counter_nue_cc_mixed,
	                                                    hit_threshold_counter_nue_cc_out_fv,
	                                                    hit_threshold_counter_cosmic,
	                                                    hit_threshold_counter_nue_nc,
	                                                    hit_threshold_counter_numu_cc,
	                                                    hit_threshold_counter_numu_cc_mixed,
	                                                    hit_threshold_counter_numu_nc,
	                                                    hit_threshold_counter_unmatched,
	                                                    hit_threshold_counter_other_mixed,
	                                                    "Hit Threshold"
	                                                    );

	_functions_instance.selection_functions::PrintInfo( total_mc_entries_inFV,
	                                                    open_angle_counter,
	                                                    open_angle_counter_nue_cc,
	                                                    open_angle_counter_nue_cc_mixed,
	                                                    open_angle_counter_nue_cc_out_fv,
	                                                    open_angle_counter_cosmic,
	                                                    open_angle_counter_nue_nc,
	                                                    open_angle_counter_numu_cc,
	                                                    open_angle_counter_numu_cc_mixed,
	                                                    open_angle_counter_numu_nc,
	                                                    open_angle_counter_unmatched,
	                                                    open_angle_counter_other_mixed,
	                                                    "Open Angle"
	                                                    );

	std::vector<double> * xsec_cc = new std::vector<double>;
	const double final_counter                = open_angle_counter;
	const double final_counter_nue_cc         = open_angle_counter_nue_cc;
	const double final_counter_nue_cc_mixed   = open_angle_counter_nue_cc_mixed;
	const double final_counter_nue_cc_out_fv  = open_angle_counter_nue_cc_out_fv;
	const double final_counter_cosmic         = open_angle_counter_cosmic;
	const double final_counter_nue_nc         = open_angle_counter_nue_nc;
	const double final_counter_numu_cc        = open_angle_counter_numu_cc;
	const double final_counter_numu_cc_mixed  = open_angle_counter_numu_cc_mixed;
	const double final_counter_numu_nc        = open_angle_counter_numu_nc;
	const double final_counter_unmatched      = open_angle_counter_unmatched;
	const double final_counter_other_mixed    = open_angle_counter_other_mixed;
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
	h_nue_eng_eff_num->Draw();
	test_c1->Print("selected_true_neutrino_energy.pdf");
	TCanvas * test_c2 = new TCanvas();
	test_c2->cd();
	h_nue_eng_eff_den->Draw();
	test_c2->Print("all_true_neutrino_energy.pdf");


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
	TEfficiency * dir_y_eff = new TEfficiency(*h_nue_dir_x_eff_num, *h_nue_dir_x_eff_den);
	dir_y_eff->SetTitle(";True Neutrino Dir Y;Efficiency");
	dir_y_eff->SetLineColor(kGreen+3);
	dir_y_eff->SetMarkerColor(kGreen+3);
	dir_y_eff->SetMarkerStyle(20);
	dir_y_eff->SetMarkerSize(0.5);
	dir_y_eff->Draw("AP");
	efficiency_c6->Print("signal_selection_nu_dir_y_efficiency.pdf");
	TCanvas * efficiency_c7 = new TCanvas();
	efficiency_c7->cd();
	TEfficiency * dir_z_eff = new TEfficiency(*h_nue_dir_x_eff_num, *h_nue_dir_x_eff_den);
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
	open_angle_stack->GetXaxis()->SetTitle("Shower Opening Angle [Degrees]");
	open_angle_stack->Draw();

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
	trk_vtx_dist_stack->GetXaxis()->SetTitle("Track to Nue Candidate Vertex Distance [cm]");
	trk_vtx_dist_stack->Draw();

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
