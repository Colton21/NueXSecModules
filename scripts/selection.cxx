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

	TH1D * h_nue_eng_eff_den = new TH1D("h_nue_eng_eff_den", "h_nue_eng_eff_den", 12, 0, 6);
	TH1D * h_nue_eng_eff_num = new TH1D("h_nue_eng_eff_num", "h_nue_eng_eff_num", 12, 0, 6);
	TH2I * h_tracks_showers  = new TH2I("h_tracks_showers", "h_tracks_showers", 10, 0, 10, 10, 0, 10);

	std::cout << "Running With: " << POT << " POT " << std::endl;
	const double flux = POT * scaling;

	mctruth_counter_tree->SetBranchAddress("mc_nue_cc_counter",      &mc_nue_cc_counter);
	mctruth_counter_tree->SetBranchAddress("mc_nue_nc_counter",      &mc_nue_nc_counter);
	mctruth_counter_tree->SetBranchAddress("mc_numu_cc_counter",     &mc_numu_cc_counter);
	mctruth_counter_tree->SetBranchAddress("mc_numu_nc_counter",     &mc_numu_nc_counter);
	mctruth_counter_tree->SetBranchAddress("mc_nue_cc_counter_bar",  &mc_nue_cc_counter_bar);
	mctruth_counter_tree->SetBranchAddress("mc_numu_cc_counter_bar", &mc_numu_cc_counter_bar);
	mctruth_counter_tree->SetBranchAddress("mc_nue_nc_counter_bar",  &mc_nue_nc_counter_bar);
	mctruth_counter_tree->SetBranchAddress("mc_numu_nc_counter_bar", &mc_numu_nc_counter_bar);
	mctruth_counter_tree->SetBranchAddress("fMCNuEnegy", &mc_nu_energy);
	mctruth_counter_tree->SetBranchAddress("fMCNuID", &mc_nu_id);
	mctruth_counter_tree->SetBranchAddress("fMCNuVtxX", &mc_nu_vtx_x);
	mctruth_counter_tree->SetBranchAddress("fMCNuVtxY", &mc_nu_vtx_y);
	mctruth_counter_tree->SetBranchAddress("fMCNuVtxZ", &mc_nu_vtx_z);

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
		if(mc_nu_id == 1) //if this event is a true nue CC interaction
		{
			if(_functions_instance.selection_functions::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
			                                                  _x1, _x2, _y1,
			                                                  _y2, _z1, _z2) == true)
			{
				h_nue_eng_eff_den->Fill(mc_nu_energy);
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
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco);
		reco_nue_counter_nue_cc       = reco_nue_counter_nue_cc + tabulated_origins.at(0);
		reco_nue_counter_nue_cc_mixed = reco_nue_counter_nue_cc_mixed + tabulated_origins.at(1);
		reco_nue_counter_cosmic       = reco_nue_counter_cosmic + tabulated_origins.at(2);
		reco_nue_counter_nue_nc       = reco_nue_counter_nue_nc + tabulated_origins.at(3);
		reco_nue_counter_numu         = reco_nue_counter_numu + tabulated_origins.at(4);
		reco_nue_counter_unmatched    = reco_nue_counter_unmatched + tabulated_origins.at(5);
		reco_nue_counter_other_mixed  = reco_nue_counter_other_mixed + tabulated_origins.at(6);
		reco_nue_counter = reco_nue_counter + tabulated_origins.at(7);

		//in fv cut
		_functions_instance.selection_functions::fiducial_volume_cut(tpc_object_container_v, _x1, _x2, _y1, _y2, _z1, _z2, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco);
		in_fv_counter_nue_cc       = in_fv_counter_nue_cc + tabulated_origins.at(0);
		in_fv_counter_nue_cc_mixed = in_fv_counter_nue_cc_mixed + tabulated_origins.at(1);
		in_fv_counter_cosmic       = in_fv_counter_cosmic + tabulated_origins.at(2);
		in_fv_counter_nue_nc       = in_fv_counter_nue_nc + tabulated_origins.at(3);
		in_fv_counter_numu         = in_fv_counter_numu + tabulated_origins.at(4);
		in_fv_counter_unmatched    = in_fv_counter_unmatched + tabulated_origins.at(5);
		in_fv_counter_other_mixed  = in_fv_counter_other_mixed + tabulated_origins.at(6);
		in_fv_counter = in_fv_counter + tabulated_origins.at(7);

		//vertex to flash cut
		_functions_instance.selection_functions::flashRecoVtxDist(largest_flash_v, tpc_object_container_v,
		                                                          tolerance, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco);
		vtx_flash_counter_nue_cc       = vtx_flash_counter_nue_cc + tabulated_origins.at(0);
		vtx_flash_counter_nue_cc_mixed = vtx_flash_counter_nue_cc_mixed + tabulated_origins.at(1);
		vtx_flash_counter_cosmic       = vtx_flash_counter_cosmic + tabulated_origins.at(2);
		vtx_flash_counter_nue_nc       = vtx_flash_counter_nue_nc + tabulated_origins.at(3);
		vtx_flash_counter_numu         = vtx_flash_counter_numu + tabulated_origins.at(4);
		vtx_flash_counter_unmatched    = vtx_flash_counter_unmatched + tabulated_origins.at(5);
		vtx_flash_counter_other_mixed  = vtx_flash_counter_other_mixed + tabulated_origins.at(6);
		vtx_flash_counter = vtx_flash_counter + tabulated_origins.at(7);

		//distance between pfp shower and nue object cut
		_functions_instance.selection_functions::VtxNuDistance(tpc_object_container_v, shwr_nue_tolerance, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco);
		shwr_tpco_counter_nue_cc       = shwr_tpco_counter_nue_cc + tabulated_origins.at(0);
		shwr_tpco_counter_nue_cc_mixed = shwr_tpco_counter_nue_cc_mixed + tabulated_origins.at(1);
		shwr_tpco_counter_cosmic       = shwr_tpco_counter_cosmic + tabulated_origins.at(2);
		shwr_tpco_counter_nue_nc       = shwr_tpco_counter_nue_nc + tabulated_origins.at(3);
		shwr_tpco_counter_numu         = shwr_tpco_counter_numu + tabulated_origins.at(4);
		shwr_tpco_counter_unmatched    = shwr_tpco_counter_unmatched + tabulated_origins.at(5);
		shwr_tpco_counter_other_mixed  = shwr_tpco_counter_other_mixed + tabulated_origins.at(6);
		shwr_tpco_counter = shwr_tpco_counter + tabulated_origins.at(7);

		//hit threshold for showers cut
		_functions_instance.selection_functions::HitThreshold(tpc_object_container_v, shwr_hit_threshold, passed_tpco, _verbose);
		if(_functions_instance.selection_functions::ValidTPCObjects(passed_tpco) == false) {continue; }
		tabulated_origins = _functions_instance.selection_functions::TabulateOrigins(tpc_object_container_v, passed_tpco);
		hit_threshold_counter_nue_cc       = hit_threshold_counter_nue_cc + tabulated_origins.at(0);
		hit_threshold_counter_nue_cc_mixed = hit_threshold_counter_nue_cc_mixed + tabulated_origins.at(1);
		hit_threshold_counter_cosmic       = hit_threshold_counter_cosmic + tabulated_origins.at(2);
		hit_threshold_counter_nue_nc       = hit_threshold_counter_nue_nc + tabulated_origins.at(3);
		hit_threshold_counter_numu         = hit_threshold_counter_numu + tabulated_origins.at(4);
		hit_threshold_counter_unmatched    = hit_threshold_counter_unmatched + tabulated_origins.at(5);
		hit_threshold_counter_other_mixed  = hit_threshold_counter_other_mixed + tabulated_origins.at(6);
		hit_threshold_counter = hit_threshold_counter + tabulated_origins.at(7);

		if(mc_nu_id == 1 && tabulated_origins.at(0) == 1) {h_nue_eng_eff_num->Fill(mc_nu_energy); }
		_functions_instance.selection_functions::PostCutPlots(tpc_object_container_v, passed_tpco, _verbose, h_tracks_showers);

	}//end event loop

	std::cout << "------------------ " << std::endl;
	std::cout << " MC Entries in FV: " << total_mc_entries_inFV << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "------------------ " << std::endl;
	std::cout << "  End Selection    " << std::endl;
	std::cout << "------------------ " << std::endl;

	//change mc_nue_cc_counter to total_mc_entries_inFV once files are ready!

	//we also want some metrics to print at the end
	_functions_instance.selection_functions::PrintInfo( mc_nue_cc_counter,
	                                                    reco_nue_counter,
	                                                    reco_nue_counter_nue_cc,
	                                                    reco_nue_counter_nue_cc_mixed,
	                                                    reco_nue_counter_cosmic,
	                                                    reco_nue_counter_nue_nc,
	                                                    reco_nue_counter_numu,
	                                                    reco_nue_counter_unmatched,
	                                                    reco_nue_counter_other_mixed,
	                                                    "Reco Nue"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( mc_nue_cc_counter,
	                                                    in_fv_counter,
	                                                    in_fv_counter_nue_cc,
	                                                    in_fv_counter_nue_cc_mixed,
	                                                    in_fv_counter_cosmic,
	                                                    in_fv_counter_nue_nc,
	                                                    in_fv_counter_numu,
	                                                    in_fv_counter_unmatched,
	                                                    in_fv_counter_other_mixed,
	                                                    "In FV"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( mc_nue_cc_counter,
	                                                    vtx_flash_counter,
	                                                    vtx_flash_counter_nue_cc,
	                                                    vtx_flash_counter_nue_cc_mixed,
	                                                    vtx_flash_counter_cosmic,
	                                                    vtx_flash_counter_nue_nc,
	                                                    vtx_flash_counter_numu,
	                                                    vtx_flash_counter_unmatched,
	                                                    vtx_flash_counter_other_mixed,
	                                                    "Vtx-to-Flash"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( mc_nue_cc_counter,
	                                                    shwr_tpco_counter,
	                                                    shwr_tpco_counter_nue_cc,
	                                                    shwr_tpco_counter_nue_cc_mixed,
	                                                    shwr_tpco_counter_cosmic,
	                                                    shwr_tpco_counter_nue_nc,
	                                                    shwr_tpco_counter_numu,
	                                                    shwr_tpco_counter_unmatched,
	                                                    shwr_tpco_counter_other_mixed,
	                                                    "Shower-to-TPCO"
	                                                    );
	_functions_instance.selection_functions::PrintInfo( mc_nue_cc_counter,
	                                                    hit_threshold_counter,
	                                                    hit_threshold_counter_nue_cc,
	                                                    hit_threshold_counter_nue_cc_mixed,
	                                                    hit_threshold_counter_cosmic,
	                                                    hit_threshold_counter_nue_nc,
	                                                    hit_threshold_counter_numu,
	                                                    hit_threshold_counter_unmatched,
	                                                    hit_threshold_counter_other_mixed,
	                                                    "Hit Threshold"
	                                                    );

	std::vector<double> * xsec_cc = new std::vector<double>;
	const double final_counter = hit_threshold_counter;
	const double final_counter_nue_cc = hit_threshold_counter_nue_cc;
	const double final_counter_nue_cc_mixed = hit_threshold_counter_nue_cc_mixed;
	const double final_counter_cosmic = hit_threshold_counter_cosmic;
	const double final_counter_nue_nc = hit_threshold_counter_nue_nc;
	const double final_counter_numu = hit_threshold_counter_numu;
	const double final_counter_unmatched = hit_threshold_counter_unmatched;
	const double final_counter_other_mixed = hit_threshold_counter_other_mixed;
	const int n_total = final_counter;
	const int n_bkg = (final_counter_nue_cc_mixed + final_counter_cosmic + final_counter_nue_nc
	                   + final_counter_numu + final_counter_unmatched + final_counter_other_mixed);
	const double efficiency = final_counter_nue_cc / double(mc_nue_cc_counter);
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
	                                                  mc_nue_cc_counter, 0, flux,
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


	TCanvas * efficency_c1 = new TCanvas();
	efficency_c1->cd();
	TEfficiency * eng_eff = new TEfficiency(*h_nue_eng_eff_num, *h_nue_eng_eff_den);
	eng_eff->SetTitle(";True Neutrino Energy [GeV];Efficiency");
	eng_eff->SetLineColor(kGreen+3);
	eng_eff->SetMarkerColor(kGreen+3);
	eng_eff->SetMarkerStyle(20);
	eng_eff->SetMarkerSize(0.5);
	eng_eff->Draw("AP");
	efficency_c1->Print("signal_selection_nu_energy_efficiency.pdf");

	_functions_instance.selection_functions::xsec_plot(_verbose, genie_xsec, xsec_cc->at(1));

	TCanvas * post_cuts_c1 = new TCanvas();
	post_cuts_c1->cd();
	h_tracks_showers->GetYaxis()->SetTitle("Reco Showers");
	h_tracks_showers->GetXaxis()->SetTitle("Reco Tracks");
	h_tracks_showers->SetTitle("Post Cuts - Showers/Tracks per Candidate Nue TPC Object");
	h_tracks_showers->SetStats(kFALSE);
	h_tracks_showers->Draw("colz");
	post_cuts_c1->Print("post_cuts_showers_tracks.pdf");


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
