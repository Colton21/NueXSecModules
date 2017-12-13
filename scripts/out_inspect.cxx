#include "../xsecAna/TpcObjectContainer.h"
#include "../xsecAna/ParticleContainer.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"

#include <iostream>
#include <vector>

#include "../xsecAna/LinkDef.h"

int out_inspect(const char * _file1)
{
	//const char * _file1 = "../nue_xsec_extraction.root";
	//const char * _file1 = "../cosmic_extraction.root";
	std::cout << "File Path: " << _file1 << std::endl;
	//first we need to open the root file
	TFile * f = new TFile(_file1);
	if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return 1; }
	TTree * mytree = (TTree*)f->Get("AnalyzeTPCO/tree");
	TTree * optree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

	std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
	mytree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;
	double fOpFlashWidthY = 0;
	double fOpFlashWidthZ = 0;
	double fOpFlashCenterY = 0;
	double fOpFlashCenterZ = 0;
	optree->SetBranchAddress("Run", &fRun);
	optree->SetBranchAddress("Event", &fEvent);
	optree->SetBranchAddress("OpFlashPE",        &fOpFlashPE);
	optree->SetBranchAddress("OpFlashTime",      &fOpFlashTime);
	optree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY);
	optree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ);
	optree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY);
	optree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ);

	std::vector< int > reco_nue_v_pfp_pdg;
	std::vector< int > reco_nue_v_mc_pdg;
	std::vector< std::string > reco_nue_v_origin;
	std::vector< int > reco_nue_v_pfp_hits;
	std::vector < std::vector < double > > reco_nue_shower_vtx_v_v;

	TFile * histo_file = new TFile("histo_file.root", "RECREATE");
	const int total_entries = mytree->GetEntries();

	TH1I * h_tpco_origin_unknown  = new TH1I("h_tpco_origin_unknown", "h_tpco_origin_unknown", 10, 0, 10);
	TH1I * h_tpco_origin_cosmic   = new TH1I("h_tpco_origin_cosmic", "h_tpco_origin_cosmic", 10, 0, 10);
	TH1I * h_tpco_origin_neutrino = new TH1I("h_tpco_origin_neutrino", "h_tpco_origin_neutrino", 10, 0, 10);
	TH2I * h_tpco_origin_neutrino_cosmic = new TH2I("h_tpco_origin_neutrino_cosmic", "h_tpco_origin_neutrino_cosmic", 10, 0, 10, 10, 0, 10);
	TH2I * h_nue_daughter_track_shower = new TH2I("h_nue_daughter_track_shower",
	                                              "h_nue_daughter_track_shower", 5, 0, 5, 5, 0, 5);
	TH2I * h_leading_shower_track_shower_beam      = new TH2I("h_leading_shower_track_shower_beam",
	                                                          "h_leading_shower_track_shower_beam", 5, 0, 5, 5, 0, 5);
	TH2I * h_leading_shower_track_shower_cosmic    = new TH2I("h_leading_shower_track_shower_cosmic",
	                                                          "h_leading_shower_track_shower_cosmic", 5, 0, 5, 5, 0, 5);
	TH2I * h_leading_shower_track_shower_unknown    = new TH2I("h_leading_shower_track_shower_unknown",
	                                                           "h_leading_shower_track_shower_unknown", 5, 0, 5, 5, 0, 5);

	TH1I * h_tracks_tpco_cat_nue_cc       = new TH1I("h_tracks_tpco_cat_nue_cc", "h_tracks_tpco_cat_nue_cc", 5, 0, 5);
	TH1I * h_tracks_tpco_cat_nue_cc_mixed = new TH1I("h_tracks_tpco_cat_nue_cc", "h_tracks_tpco_cat_nue_cc_mixed", 5, 0, 5);
	TH1I * h_tracks_tpco_cat_cosmic       = new TH1I("h_tracks_tpco_cat_nue_cc", "h_tracks_tpco_cat_cosmic", 5, 0, 5);
	TH1I * h_tracks_tpco_cat_nue_nc       = new TH1I("h_tracks_tpco_cat_nue_cc", "h_tracks_tpco_cat_nue_nc", 5, 0, 5);
	TH1I * h_tracks_tpco_cat_numu         = new TH1I("h_tracks_tpco_cat_nue_cc", "h_tracks_tpco_cat_numu", 5, 0, 5);
	TH1I * h_tracks_tpco_cat_other_mixed  = new TH1I("h_tracks_tpco_cat_nue_cc", "h_tracks_tpco_cat_other_mixed", 5, 0, 5);
	TH1I * h_tracks_tpco_cat_unmatched    = new TH1I("h_tracks_tpco_cat_nue_cc", "h_tracks_tpco_cat_unmatched", 5, 0, 5);

	TH1I * h_showers_tpco_cat_nue_cc       = new TH1I("h_showers_tpco_cat_nue_cc", "h_showers_tpco_cat_nue_cc", 5, 0, 5);
	TH1I * h_showers_tpco_cat_nue_cc_mixed = new TH1I("h_showers_tpco_cat_nue_cc", "h_showers_tpco_cat_nue_cc_mixed", 5, 0, 5);
	TH1I * h_showers_tpco_cat_cosmic       = new TH1I("h_showers_tpco_cat_nue_cc", "h_showers_tpco_cat_cosmic", 5, 0, 5);
	TH1I * h_showers_tpco_cat_nue_nc       = new TH1I("h_showers_tpco_cat_nue_cc", "h_showers_tpco_cat_nue_nc", 5, 0, 5);
	TH1I * h_showers_tpco_cat_numu         = new TH1I("h_showers_tpco_cat_nue_cc", "h_showers_tpco_cat_numu", 5, 0, 5);
	TH1I * h_showers_tpco_cat_other_mixed  = new TH1I("h_showers_tpco_cat_nue_cc", "h_showers_tpco_cat_other_mixed", 5, 0, 5);
	TH1I * h_showers_tpco_cat_unmatched    = new TH1I("h_showers_tpco_cat_nue_cc", "h_showers_tpco_cat_unmatched", 5, 0, 5);

	// h_tpco_origin_unknown->SetStats(kFALSE);
	// h_tpco_origin_cosmic->SetStats(kFALSE);
	// h_tpco_origin_neutrino->SetStats(kFALSE);
	h_tpco_origin_neutrino_cosmic->SetStats(kFALSE);

	TH1I * h_leading_shower_origin = new TH1I("h_leading_shower_origin", "h_leading_shower_origin", 3, 0, 3);
	TH1I * h_leading_shower_mc_pdg = new TH1I("h_leading_shower_mc_pdg", "h_leading_shower_mc_pdg", 10, 0, 10);
	TH2I * h_leading_shower_origin_mc_pdg = new TH2I("h_leading_shower_origin_mc_pdg", "h_leading_shower_origin_mc_pdg", 3, 0, 3, 10, 0, 10);
	TH2I * h_leading_shower_origin_hits = new TH2I("h_leading_shower_origin_hits", "h_leading_shower_origin_hits", 3, 0, 3, 100, 0, 3000);
	TH2I * h_leading_shower_origin_tpco_cat = new TH2I("h_leading_shower_origin_tpco_cat", "h_leading_shower_origin_tpco_cat", 3, 0, 3, 7, 0, 7);

	TH2D * h_leading_shower_mc_pdg_pfp_hits_nue_cc        = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_nue_cc",
	                                                                  "h_leading_shower_mc_pdg_pfp_hits_nue_cc", 20, 0, 3000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed  = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed",
	                                                                  "h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed", 20, 0, 3000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_cosmic        = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_cosmic",
	                                                                  "h_leading_shower_mc_pdg_pfp_hits_cosmic", 20, 0, 3000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_nue_nc        = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_nue_nc",
	                                                                  "h_leading_shower_mc_pdg_pfp_hits_nue_nc", 20, 0, 3000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_numu          = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_numu",
	                                                                  "h_leading_shower_mc_pdg_pfp_hits_numu", 20, 0, 3000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_other_mixed   = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_other_mixed",
	                                                                  "h_leading_shower_mc_pdg_pfp_hits_other_mixed", 20, 0, 3000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_unmatched     = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_unmatched",
	                                                                  "h_leading_shower_mc_pdg_pfp_hits_unmatched", 20, 0, 3000, 10, 0, 10);

	TH2D * h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom        = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom",
	                                                                       "h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom", 20, 0, 1000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom  = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom",
	                                                                       "h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom", 20, 0, 1000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom        = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom",
	                                                                       "h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom", 20, 0, 1000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom        = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom",
	                                                                       "h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom", 20, 0, 1000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_numu_zoom          = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_numu_zoom",
	                                                                       "h_leading_shower_mc_pdg_pfp_hits_numu_zoom", 20, 0, 1000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom   = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom",
	                                                                       "h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom", 20, 0, 1000, 10, 0, 10);
	TH2D * h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom     = new TH2D ("h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom",
	                                                                       "h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom", 20, 0, 1000, 10, 0, 10);

	TH1D * h_shwr_to_vtx_beam_electron = new TH1D("h_shwr_to_vtx_beam_electron", "h_shwr_to_vtx_beam_electron", 40, 0, 80);
	TH1D * h_shwr_to_vtx_beam_photon = new TH1D("h_shwr_to_vtx_beam_photon", "h_shwr_to_vtx_beam_photon", 40, 0, 80);
	TH1D * h_shwr_to_vtx_cosmic_electron = new TH1D("h_shwr_to_vtx_cosmic_electron", "h_shwr_to_vtx_cosmic_electron", 40, 0, 80);
	TH1D * h_shwr_to_vtx_cosmic_photon = new TH1D("h_shwr_to_vtx_cosmic_photon", "h_shwr_to_vtx_cosmic_photon", 40, 0, 80);
	TH1D * h_shwr_to_vtx_unmatched = new TH1D("h_shwr_to_vtx_unmatched", "h_shwr_to_vtx_unmatched", 40, 0, 80);

	TH2D * h_shwr_to_vtx_beam_electron_tracks     = new TH2D("h_shwr_to_vtx_beam_electron_tracks", "h_shwr_to_vtx_beam_electron_tracks", 40, 0, 160, 8, 0, 8);
	TH2D * h_shwr_to_vtx_beam_photon_tracks       = new TH2D("h_shwr_to_vtx_beam_photon_tracks", "h_shwr_to_vtx_beam_photon_tracks", 40, 0, 160, 8, 0, 8);
	TH2D * h_shwr_to_vtx_cosmic_electron_tracks   = new TH2D("h_shwr_to_vtx_cosmic_electron_tracks", "h_shwr_to_vtx_cosmic_electron_tracks", 40, 0, 160, 8, 0, 8);
	TH2D * h_shwr_to_vtx_cosmic_photon_tracks     = new TH2D("h_shwr_to_vtx_cosmic_photon_tracks", "h_shwr_to_vtx_cosmic_photon_tracks", 40, 0, 160, 8, 0, 8);
	TH2D * h_shwr_to_vtx_unmatched_tracks         = new TH2D("h_shwr_to_vtx_unmatched_tracks", "h_shwr_to_vtx_unmatched_tracks", 40, 0, 160, 8, 0, 8);

	TH2D * h_shwr_to_vtx_beam_electron_showers     = new TH2D("h_shwr_to_vtx_beam_electron_showers", "h_shwr_to_vtx_beam_electron_showers", 40, 0, 160, 8, 0, 8);
	TH2D * h_shwr_to_vtx_beam_photon_showers       = new TH2D("h_shwr_to_vtx_beam_photon_showers", "h_shwr_to_vtx_beam_photon_showers", 40, 0, 160, 8, 0, 8);
	TH2D * h_shwr_to_vtx_cosmic_electron_showers   = new TH2D("h_shwr_to_vtx_cosmic_electron_showers", "h_shwr_to_vtx_cosmic_electron_showers", 40, 0, 160, 8, 0, 8);
	TH2D * h_shwr_to_vtx_cosmic_photon_showers     = new TH2D("h_shwr_to_vtx_cosmic_photon_showers", "h_shwr_to_vtx_cosmic_photon_showers", 40, 0, 160, 8, 0, 8);
	TH2D * h_shwr_to_vtx_unmatched_showers         = new TH2D("h_shwr_to_vtx_unmatched_showers", "h_shwr_to_vtx_unmatched_showers", 40, 0, 160, 8, 0, 8);

	TH1D * h_shwr_open_angle_nue_cc        = new TH1D("h_shwr_open_angle_nue_cc", "h_shwr_open_angle_nue_cc", 20, 0, 60);
	TH1D * h_shwr_open_angle_nue_cc_mixed  = new TH1D("h_shwr_open_angle_nue_cc_mixed", "h_shwr_open_angle_nue_cc_mixed", 20, 0, 60);
	TH1D * h_shwr_open_angle_numu_cc       = new TH1D("h_shwr_open_angle_numu_cc", "h_shwr_open_angle_numu_cc", 20, 0, 60);
	TH1D * h_shwr_open_angle_numu_cc_mixed = new TH1D("h_shwr_open_angle_numu_cc_mixed", "h_shwr_open_angle_numu_cc_mixed", 20, 0, 60);
	TH1D * h_shwr_open_angle_numu_nc       = new TH1D("h_shwr_open_angle_numu_nc", "h_shwr_open_angle_numu_nc", 20, 0, 60);
	TH1D * h_shwr_open_angle_cosmic        = new TH1D("h_shwr_open_angle_cosmic", "h_shwr_open_angle_cosmic", 20, 0, 60);

	//leading shower dEdx plots
	TH1D * h_shwr_dedx_nue_cc          = new TH1D("h_shwr_dedx_nue_cc", "h_shwr_dedx_nue_cc", 20, 0, 10);
	TH1D * h_shwr_dedx_nue_cc_mixed    = new TH1D("h_shwr_dedx_nue_cc_mixed", "h_shwr_dedx_nue_cc_mixed", 20, 0, 10);
	TH1D * h_shwr_dedx_nue_nc          = new TH1D("h_shwr_dedx_nue_nc", "h_shwr_dedx_nue_nc", 20, 0, 10);
	TH1D * h_shwr_dedx_numu_cc         = new TH1D("h_shwr_dedx_numu_cc", "h_shwr_dedx_numu_cc", 20, 0, 10);
	TH1D * h_shwr_dedx_numu_cc_mixed   = new TH1D("h_shwr_dedx_numu_cc_mixed", "h_shwr_dedx_numu_cc_mixed", 20, 0, 10);
	TH1D * h_shwr_dedx_numu_nc         = new TH1D("h_shwr_dedx_numu_nc", "h_shwr_dedx_numu_nc", 20, 0, 10);
	TH1D * h_shwr_dedx_cosmic          = new TH1D("h_shwr_dedx_cosmic", "h_shwr_dedx_cosmic", 20, 0, 10);
	TH1D * h_shwr_dedx_other_mixed     = new TH1D("h_shwr_dedx_other_mixed", "h_shwr_dedx_other_mixed", 20, 0, 10);
	TH1D * h_shwr_dedx_unmatched       = new TH1D("h_shwr_dedx_unmatched", "h_shwr_dedx_unmatched", 20, 0, 10);

	//leading shower dEdx vs cos(theta)
	TH2D * h_shwr_dedx_theta_nue_cc          = new TH2D("h_shwr_dedx_theta_nue_cc", "h_shwr_dedx_theta_nue_cc", 20, 0, 10, 10, -1, 1);
	TH2D * h_shwr_dedx_theta_nue_cc_mixed    = new TH2D("h_shwr_dedx_theta_nue_cc_mixed", "h_shwr_dedx_theta_nue_cc_mixed", 20, 0, 10, 10, -1, 1);
	TH2D * h_shwr_dedx_theta_nue_nc          = new TH2D("h_shwr_dedx_theta_nue_nc", "h_shwr_dedx_theta_nue_nc", 20, 0, 10, 10, -1, 1);
	TH2D * h_shwr_dedx_theta_numu_cc         = new TH2D("h_shwr_dedx_theta_numu_cc", "h_shwr_dedx_theta_numu_cc", 20, 0, 10, 10, -1, 1);
	TH2D * h_shwr_dedx_theta_numu_cc_mixed   = new TH2D("h_shwr_dedx_theta_numu_cc_mixed", "h_shwr_dedx_theta_numu_cc_mixed", 20, 0, 10, 10, -1, 1);
	TH2D * h_shwr_dedx_theta_numu_nc         = new TH2D("h_shwr_dedx_theta_numu_nc", "h_shwr_dedx_theta_numu_nc", 20, 0, 10, 10, -1, 1);
	TH2D * h_shwr_dedx_theta_cosmic          = new TH2D("h_shwr_dedx_theta_cosmic", "h_shwr_dedx_theta_cosmic", 20, 0, 10, 10, -1, 1);
	TH2D * h_shwr_dedx_theta_other_mixed     = new TH2D("h_shwr_dedx_theta_other_mixed", "h_shwr_dedx_theta_other_mixed", 20, 0, 10, 10, -1, 1);
	TH2D * h_shwr_dedx_theta_unmatched       = new TH2D("h_shwr_dedx_theta_unmatched", "h_shwr_dedx_theta_unmatched", 20, 0, 10, 10, -1,1);

	//THStack *hs = new THStack(h_title,"");
	int nue_cc = 0;
	int nue_cc_mixed = 0;
	int cosmic = 0;
	int nue_cc_out_fv = 0;
	int nue_nc        = 0;
	int numu_cc       = 0;
	int numu_cc_mixed = 0;
	int numu_nc       = 0;
	int unmatched     = 0;
	int other_mixed   = 0;
	int total = 0;

	for(int i = 0; i < total_entries; i++)
	{
		mytree->GetEntry(i);
		const int n_tpc_obj = tpc_object_container_v->size();
		std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
		for(auto const tpc_obj : *tpc_object_container_v)
		{
			const int run_number = tpc_obj.RunNumber();
			const int event_number = tpc_obj.EventNumber();
			const int n_pfp = tpc_obj.NumPFParticles();
			const int n_pfp_nu = tpc_obj.NumPFPNeutrinos();
			//const int n_mc_hits = tpc_obj.NumMCHits();
			const int n_pfp_hits = tpc_obj.NumPFPHits();
			const int is_cc = tpc_obj.CCNC();//0 for cc, 1 for nc
			const int mode = tpc_obj.Mode();
			const int tpco_pfp_pdg_code = tpc_obj.PFParticlePdgCode();

			//we only care about the nue tpco!
			if(tpco_pfp_pdg_code != 12 && tpco_pfp_pdg_code != -12) {continue; }

			const std::string tpc_obj_origin = tpc_obj.Origin();
			std::vector < double > tpc_obj_pfp_vtx;
			if(!tpc_obj_pfp_vtx.empty()) {tpc_obj_pfp_vtx.clear(); }
			tpc_obj_pfp_vtx.push_back(tpc_obj.pfpVtxX());
			tpc_obj_pfp_vtx.push_back(tpc_obj.pfpVtxY());
			tpc_obj_pfp_vtx.push_back(tpc_obj.pfpVtxZ());
			std::vector < double > tpc_obj_mc_vtx;
			if(!tpc_obj_mc_vtx.empty()) {tpc_obj_mc_vtx.clear(); }
			tpc_obj_mc_vtx.push_back(tpc_obj.mcVtxX());
			tpc_obj_mc_vtx.push_back(tpc_obj.mcVtxY());
			tpc_obj_mc_vtx.push_back(tpc_obj.mcVtxZ());

			int origin_cosmic = 0;
			int origin_beam = 0;
			int origin_unknown = 0;
			int n_showers = 0;
			int n_tracks = 0;

			std::cout << " \t Number of PFParticles: " << n_pfp << "\t Number of PFP Neutrinos: " << n_pfp_nu << std::endl;
			std::cout << " \t Run: " << run_number << " , Event: " << event_number << std::endl;
			std::cout << " \t Origin : " << tpc_obj_origin << std::endl;
			std::cout << " \t PDG Code: - Reco " << tpco_pfp_pdg_code << std::endl;
			std::cout << " \t Reco Hits: " << n_pfp_hits << std::endl;
			std::cout << " \t Interaction Mode: " << mode << " CC/NC: " << is_cc << std::endl;
			std::cout << " \t Vertex - Reco : " << tpc_obj_pfp_vtx.at(0) << ", " << tpc_obj_pfp_vtx.at(1) << ", " << tpc_obj_pfp_vtx.at(2) << std::endl;
			std::cout << " \t Vertex - MC   : " << tpc_obj_mc_vtx.at(0)  << ", " << tpc_obj_mc_vtx.at(1)  << ", " << tpc_obj_mc_vtx.at(2)  << std::endl;

			xsecAna::ParticleContainer leading_shower;
			int leading_index = 0;
			int most_hits = 0;
			int part_nue_cc        = 0;
			int part_cosmic        = 0;
			int part_nue_nc        = 0;
			int part_numu_cc       = 0;
			int part_numu_nc       = 0;
			int part_unmatched     = 0;
			std::string tpco_cat = "null";
			for(int j = 0; j < n_pfp; j++)
			{
				auto const part = tpc_obj.GetParticle(j);
				const int pfp_pdg = part.PFParticlePdgCode();
				if(pfp_pdg == 13) {n_tracks++; }
				if(pfp_pdg == 11) {n_showers++; }
			}


			for(int j = 0; j < n_pfp; j++)
			{
				auto const part = tpc_obj.GetParticle(j);
				const int pfp_pdg = part.PFParticlePdgCode();
				const int mc_pdg = part.MCPdgCode();
				const int mc_parent_pdg = part.MCParentPdg();
				const int pfp_parent_pdg = part.PFParticleParentPdgCode();
				const std::string origin = part.Origin();
				const int pfp_mode = part.Mode();
				const int pfp_CCNC = part.CCNC();
				std::vector < double > pfp_vtx;
				if(!pfp_vtx.empty()) {pfp_vtx.clear(); }
				pfp_vtx.push_back(part.pfpVtxX());
				pfp_vtx.push_back(part.pfpVtxY());
				pfp_vtx.push_back(part.pfpVtxZ());
				std::vector < double > mc_vtx;
				if(!mc_vtx.empty()) {mc_vtx.clear(); }
				mc_vtx.push_back(part.mcVtxX());
				mc_vtx.push_back(part.mcVtxY());
				mc_vtx.push_back(part.mcVtxZ());
				std::vector < double > pfp_dir;
				if(!pfp_dir.empty()) {pfp_dir.clear(); }
				pfp_dir.push_back(part.pfpDirX());
				pfp_dir.push_back(part.pfpDirY());
				pfp_dir.push_back(part.pfpDirZ());
				std::vector < double > mc_dir;
				if(!mc_dir.empty()) {mc_dir.clear(); }
				mc_dir.push_back(part.mcDirX());
				mc_dir.push_back(part.mcDirY());
				mc_dir.push_back(part.mcDirZ());
				const double mc_length = part.mcLength();
				const double pfp_length = part.pfpLength();
				const double mc_energy = part.mcEnergy();
				const double mc_momentum = part.mcMomentum();
				const double pfp_momentum = part.pfpMomentum();
				const int n_pfp_hits = part.NumPFPHits();
				const int n_mc_hits = part.NumMCHits();
				const double pfp_open_angle = part.pfpOpenAngle();
				std::vector < std::vector < double > > pfp_cluster_dqdx = part.PfpClusterdQdx();
				std::vector<double> pfp_dedx = part.PfpdEdx();
				//if(mc_parent_pdg != 14) {continue; }

				std::cout << " \t \t ----------------------------------------------------" << std::endl;
				std::cout << " \t \t Particle PDG Codes:        - Reco " << pfp_pdg << " - True " << mc_pdg << std::endl;
				std::cout << " \t \t Parent Particle PDG Codes: - Reco " << pfp_parent_pdg << " - True " << mc_parent_pdg << std::endl;
				std::cout << " \t \t Origin: " << origin << std::endl;
				std::cout << " \t \t Mode: " << pfp_mode        << "  CC/NC: " << pfp_CCNC << std::endl;
				std::cout << " \t \t Vertex - Reco : " << pfp_vtx.at(0) << ", " << pfp_vtx.at(1) << ", " << pfp_vtx.at(2) << std::endl;
				std::cout << " \t \t Vertex - MC   : " << mc_vtx.at(0)  << ", " << mc_vtx.at(1)  << ", " << mc_vtx.at(2)  << std::endl;
				std::cout << " \t \t Direction - Reco : " << pfp_dir.at(0) << ", " << pfp_dir.at(1) << ", " << pfp_dir.at(2) << std::endl;
				std::cout << " \t \t Direction - MC   : " << mc_dir.at(0)  << ", " << mc_dir.at(1)  << ", " << mc_dir.at(2)  << std::endl;
				std::cout << " \t \t Length - MC : " << mc_length << " - Reco : " << pfp_length << std::endl;
				std::cout << " \t \t Momentum - MC : " << mc_momentum << " - Reco : " << pfp_momentum << std::endl;
				std::cout << " \t \t Hits - Reco : " << n_pfp_hits << " - MC: " << "NOT SET" << std::endl;//n_mc_hits << std::endl;
				std::cout << " \t \t Open Angle: " << pfp_open_angle << std::endl;
				std::cout << " \t \t ----------------------------------------------------" << std::endl;

				// for(auto const cluster_dqdx : pfp_cluster_dqdx)
				// {
				//      std::cout << " \t \t Cluster Plane 0: " << cluster_dqdx.at(0) << std::endl;
				//      std::cout << " \t \t Cluster Plane 1: " << cluster_dqdx.at(1) << std::endl;
				//      std::cout << " \t \t Cluster Plane 2: " << cluster_dqdx.at(2) << std::endl;
				// }
				std::cout << " \t \t dEdx Vector Size: " << pfp_dedx.size() << std::endl;
				std::cout << " \t \t dEdx Plane 0: " << pfp_dedx.at(0) << std::endl;
				std::cout << " \t \t dEdx Plane 1: " << pfp_dedx.at(1) << std::endl;
				std::cout << " \t \t dEdx Plane 2: " << pfp_dedx.at(2) << std::endl;

				if(n_pfp_hits > most_hits) {leading_index = j; most_hits = n_pfp_hits; }
				//categorise the tpco based on the origins of the particles it contains
				if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_cc++; }
				if(part.CCNC() == 1 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_nc++; }
				if(part.Origin() == "kBeamNeutrino" && part.CCNC() == 0 && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu_cc++; }
				if(part.Origin() == "kBeamNeutrino" && part.CCNC() == 1 && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu_nc++; }
				if(part.Origin() == "kCosmicRay") { part_cosmic++; }
				if(part.Origin() == "kUnknown")   { part_unmatched++; }

				//this should only let daughters of pfp nue's through
				//we also don't want the neutrino pfps!
				if((pfp_parent_pdg == 12 || pfp_parent_pdg == -12) && (pfp_pdg != 12 && pfp_pdg != 14))
				{
					reco_nue_v_origin.push_back(origin);
					if(origin == "kUnknown") {origin_unknown++; }
					if(origin == "kBeamNeutrino") {origin_beam++; }
					if(origin == "kCosmicRay") {origin_cosmic++; }
					reco_nue_v_pfp_pdg.push_back(pfp_pdg);
					reco_nue_v_mc_pdg.push_back(mc_pdg);
					reco_nue_v_pfp_hits.push_back(n_pfp_hits);
					const double shower_to_vtx_distance = (pow((pfp_vtx.at(0) - tpc_obj_pfp_vtx.at(0)), 2) +
					                                       pow((pfp_vtx.at(1) - tpc_obj_pfp_vtx.at(1)), 2) +
					                                       pow((pfp_vtx.at(2) - tpc_obj_pfp_vtx.at(2)), 2));
					if(pfp_pdg == 11)
					{

						std::vector < double > reco_nue_shower_vtx_v;
						reco_nue_shower_vtx_v.push_back(pfp_vtx.at(0));
						reco_nue_shower_vtx_v.push_back(pfp_vtx.at(1));
						reco_nue_shower_vtx_v.push_back(pfp_vtx.at(2));
						reco_nue_shower_vtx_v.push_back(event_number);
						reco_nue_shower_vtx_v.push_back(mc_pdg);
						if(origin == "kBeamNeutrino") reco_nue_shower_vtx_v.push_back(0);
						if(origin == "kCosmicRay")    reco_nue_shower_vtx_v.push_back(1);
						if(origin == "kUnknown")      reco_nue_shower_vtx_v.push_back(2);
						reco_nue_shower_vtx_v_v.push_back(reco_nue_shower_vtx_v);
						//n_showers++;
						if(mc_pdg == 11)//electron
						{
							if(origin == "kBeamNeutrino")
							{
								h_shwr_to_vtx_beam_electron->Fill(shower_to_vtx_distance);
								h_shwr_to_vtx_beam_electron_tracks->Fill(shower_to_vtx_distance, n_tracks);
								h_shwr_to_vtx_beam_electron_showers->Fill(shower_to_vtx_distance, n_showers);
							}
							if(origin == "kCosmicRay")
							{
								h_shwr_to_vtx_cosmic_electron->Fill(shower_to_vtx_distance);
								h_shwr_to_vtx_cosmic_electron_tracks->Fill(shower_to_vtx_distance, n_tracks);
								h_shwr_to_vtx_cosmic_electron_showers->Fill(shower_to_vtx_distance, n_showers);
							}
						}
						if(mc_pdg == 22)//photon
						{
							if(origin == "kBeamNeutrino")
							{
								h_shwr_to_vtx_beam_photon->Fill(shower_to_vtx_distance);
								h_shwr_to_vtx_beam_photon_tracks->Fill(shower_to_vtx_distance, n_tracks);
								h_shwr_to_vtx_beam_photon_showers->Fill(shower_to_vtx_distance, n_showers);
							}
							if(origin == "kCosmicRay")
							{
								h_shwr_to_vtx_cosmic_photon->Fill(shower_to_vtx_distance);
								h_shwr_to_vtx_cosmic_photon_tracks->Fill(shower_to_vtx_distance, n_tracks);
								h_shwr_to_vtx_cosmic_photon_showers->Fill(shower_to_vtx_distance, n_showers);
							}
						}
						if(mc_pdg == 0)//unmatched
						{
							h_shwr_to_vtx_unmatched->Fill(shower_to_vtx_distance);
							h_shwr_to_vtx_unmatched_tracks->Fill(shower_to_vtx_distance, n_tracks);
							h_shwr_to_vtx_unmatched_showers->Fill(shower_to_vtx_distance, n_showers);
						}
					}
				}
			}//end looping particles


			leading_shower = tpc_obj.GetParticle(leading_index);
			const std::string leading_origin = leading_shower.Origin();
			const int leading_mc_pdg = leading_shower.MCPdgCode();
			const int leading_hits = leading_shower.NumPFPHits();
			const double leading_pfp_open_angle = leading_shower.pfpOpenAngle();
			const std::vector< double > leading_dedx = leading_shower.PfpdEdx();
			const double leading_cos_theta = cos(leading_shower.pfpTheta());

			if(part_cosmic > 0)
			{
				if(part_nue_cc  > 0                       && tpco_cat == "null")    {nue_cc_mixed++;  tpco_cat = "nue_cc_mixed"; }
				if(part_numu_cc > 0                       && tpco_cat == "null" )   {numu_cc_mixed++; tpco_cat = "numu_cc_mixed"; }
				if((part_nue_nc  > 0 || part_numu_nc > 0) && tpco_cat == "null")    {other_mixed++;   tpco_cat = "other_mixed"; }
				if(tpco_cat == "null") {cosmic++; tpco_cat = "cosmic"; }
			}
			if(part_cosmic == 0)
			{
				if(part_nue_cc    > 0 && tpco_cat == "null") {nue_cc++;       tpco_cat = "nue_cc"; }
				if(part_nue_nc    > 0 && tpco_cat == "null") {nue_nc++;       tpco_cat = "nue_nc"; }
				if(part_numu_cc   > 0 && tpco_cat == "null") {numu_cc++;      tpco_cat = "numu_cc"; }
				if(part_numu_nc   > 0 && tpco_cat == "null") {numu_nc++;      tpco_cat = "numu_nc"; }
				if(part_unmatched > 0 && tpco_cat == "null") {unmatched++;    tpco_cat = "unmatched"; }
			}

			if(tpco_cat == "nue_cc" )       {h_shwr_open_angle_nue_cc->Fill        ( leading_pfp_open_angle * (180 / 3.1415)); }
			if(tpco_cat == "nue_cc_mixed")  {h_shwr_open_angle_nue_cc_mixed->Fill  ( leading_pfp_open_angle * (180 / 3.1415)); }
			if(tpco_cat == "numu_cc")       {h_shwr_open_angle_numu_cc->Fill       ( leading_pfp_open_angle * (180 / 3.1415)); }
			if(tpco_cat == "numu_cc_mixed") {h_shwr_open_angle_numu_cc_mixed->Fill ( leading_pfp_open_angle * (180 / 3.1415)); }
			if(tpco_cat == "numu_nc")       {h_shwr_open_angle_numu_nc->Fill       ( leading_pfp_open_angle * (180 / 3.1415)); }
			if(tpco_cat == "cosmic")        {h_shwr_open_angle_cosmic->Fill        ( leading_pfp_open_angle * (180 / 3.1415)); }

			if(tpco_cat == "nue_cc")
			{
				h_tracks_tpco_cat_nue_cc->Fill(n_tracks);
				h_showers_tpco_cat_nue_cc->Fill(n_showers);
				h_shwr_dedx_nue_cc->Fill(leading_dedx.at(2));//collection plane only
				h_shwr_dedx_theta_nue_cc->Fill(leading_dedx.at(2), leading_cos_theta);
				if(leading_mc_pdg == 11)                            {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 0.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 0.0); }
				if(leading_mc_pdg == -11)                           {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 1.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 1.0); }
				if(leading_mc_pdg == 13)                            {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 2.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 2.0); }
				if(leading_mc_pdg == -13)                           {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 3.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 3.0); }
				if(leading_mc_pdg == 22)                            {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 4.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 5.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 5.0); }
				if(leading_mc_pdg == 2212)                          {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 6.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 6.0); }
				if(leading_mc_pdg == 2112)                          {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 7.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 8.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_mc_pdg_pfp_hits_nue_cc->Fill(leading_hits, 9.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Fill(leading_hits, 9.0); }

			}
			if(tpco_cat == "nue_cc_mixed")
			{
				h_tracks_tpco_cat_nue_cc_mixed->Fill(n_tracks);
				h_showers_tpco_cat_nue_cc_mixed->Fill(n_showers);
				h_shwr_dedx_nue_cc_mixed->Fill(leading_dedx.at(2));
				h_shwr_dedx_theta_nue_cc_mixed->Fill(leading_dedx.at(2), leading_cos_theta);
				if(leading_mc_pdg == 11)                            {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 0.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 0.0); }
				if(leading_mc_pdg == -11)                           {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 1.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 1.0); }
				if(leading_mc_pdg == 13)                            {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 2.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 2.0); }
				if(leading_mc_pdg == -13)                           {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 3.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 3.0); }
				if(leading_mc_pdg == 22)                            {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 4.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 5.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 5.0); }
				if(leading_mc_pdg == 2212)                          {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 6.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 6.0); }
				if(leading_mc_pdg == 2112)                          {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 7.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 8.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Fill(leading_hits, 9.0); h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Fill(leading_hits, 9.0); }
			}
			if(tpco_cat == "cosmic")
			{
				h_tracks_tpco_cat_cosmic->Fill(n_tracks);
				h_showers_tpco_cat_cosmic->Fill(n_showers);
				h_shwr_dedx_cosmic->Fill(leading_dedx.at(2));
				h_shwr_dedx_theta_cosmic->Fill(leading_dedx.at(2), leading_cos_theta);
				if(leading_mc_pdg == 11)                            {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 0.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 0.0); }
				if(leading_mc_pdg == -11)                           {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 1.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 1.0); }
				if(leading_mc_pdg == 13)                            {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 2.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 2.0); }
				if(leading_mc_pdg == -13)                           {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 3.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 3.0); }
				if(leading_mc_pdg == 22)                            {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 4.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 5.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 5.0); }
				if(leading_mc_pdg == 2212)                          {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 6.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 6.0); }
				if(leading_mc_pdg == 2112)                          {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 7.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 8.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_mc_pdg_pfp_hits_cosmic->Fill(leading_hits, 9.0); h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Fill(leading_hits, 9.0); }
			}
			if(tpco_cat == "nue_nc")
			{
				h_tracks_tpco_cat_nue_nc->Fill(n_tracks);
				h_showers_tpco_cat_nue_nc->Fill(n_showers);
				h_shwr_dedx_nue_nc->Fill(leading_dedx.at(2));
				h_shwr_dedx_theta_nue_nc->Fill(leading_dedx.at(2), leading_cos_theta);
				if(leading_mc_pdg == 11)                            {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 0.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 0.0); }
				if(leading_mc_pdg == -11)                           {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 1.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 1.0); }
				if(leading_mc_pdg == 13)                            {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 2.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 2.0); }
				if(leading_mc_pdg == -13)                           {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 3.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 3.0); }
				if(leading_mc_pdg == 22)                            {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 4.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 5.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 5.0); }
				if(leading_mc_pdg == 2212)                          {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 6.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 6.0); }
				if(leading_mc_pdg == 2112)                          {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 7.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 8.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_mc_pdg_pfp_hits_nue_nc->Fill(leading_hits, 9.0); h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Fill(leading_hits, 9.0); }
			}
			if(tpco_cat == "numu_cc")
			{
				h_tracks_tpco_cat_numu->Fill(n_tracks);
				h_showers_tpco_cat_numu->Fill(n_showers);
				h_shwr_dedx_numu_cc->Fill(leading_dedx.at(2));
				h_shwr_dedx_theta_numu_cc->Fill(leading_dedx.at(2), leading_cos_theta);
				if(leading_mc_pdg == 11)                            {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 0.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 0.0); }
				if(leading_mc_pdg == -11)                           {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 1.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 1.0); }
				if(leading_mc_pdg == 13)                            {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 2.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 2.0); }
				if(leading_mc_pdg == -13)                           {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 3.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 3.0); }
				if(leading_mc_pdg == 22)                            {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 4.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 5.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 5.0); }
				if(leading_mc_pdg == 2212)                          {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 6.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 6.0); }
				if(leading_mc_pdg == 2112)                          {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 7.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 8.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_mc_pdg_pfp_hits_numu->Fill(leading_hits, 9.0); h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Fill(leading_hits, 9.0); }
			}
			if(tpco_cat == "numu_cc_mixed")
			{
				h_shwr_dedx_numu_cc_mixed->Fill(leading_dedx.at(2));
				h_shwr_dedx_theta_numu_cc_mixed->Fill(leading_dedx.at(2), leading_cos_theta);
			}
			if(tpco_cat == "numu_nc")
			{
				h_shwr_dedx_numu_nc->Fill(leading_dedx.at(2));
				h_shwr_dedx_theta_numu_nc->Fill(leading_dedx.at(2), leading_cos_theta);
			}
			if(tpco_cat == "other_mixed")
			{
				h_tracks_tpco_cat_other_mixed->Fill(n_tracks);
				h_showers_tpco_cat_other_mixed->Fill(n_showers);
				h_shwr_dedx_other_mixed->Fill(leading_dedx.at(2));
				h_shwr_dedx_theta_other_mixed->Fill(leading_dedx.at(2), leading_cos_theta);
				if(leading_mc_pdg == 11)                            {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 0.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 0.0); }
				if(leading_mc_pdg == -11)                           {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 1.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 1.0); }
				if(leading_mc_pdg == 13)                            {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 2.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 2.0); }
				if(leading_mc_pdg == -13)                           {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 3.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 3.0); }
				if(leading_mc_pdg == 22)                            {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 4.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 5.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 5.0); }
				if(leading_mc_pdg == 2212)                          {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 6.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 6.0); }
				if(leading_mc_pdg == 2112)                          {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 7.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 8.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_mc_pdg_pfp_hits_other_mixed->Fill(leading_hits, 9.0); h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Fill(leading_hits, 9.0); }
			}
			if(tpco_cat == "unmatched")
			{
				h_tracks_tpco_cat_unmatched->Fill(n_tracks);
				h_showers_tpco_cat_unmatched->Fill(n_showers);
				h_shwr_dedx_unmatched->Fill(leading_dedx.at(2));
				h_shwr_dedx_theta_unmatched->Fill(leading_dedx.at(2), leading_cos_theta);
				if(leading_mc_pdg == 11)                            {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 0.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 0.0); }
				if(leading_mc_pdg == -11)                           {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 1.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 1.0); }
				if(leading_mc_pdg == 13)                            {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 2.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 2.0); }
				if(leading_mc_pdg == -13)                           {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 3.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 3.0); }
				if(leading_mc_pdg == 22)                            {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 4.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 5.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 5.0); }
				if(leading_mc_pdg == 2212)                          {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 6.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 6.0); }
				if(leading_mc_pdg == 2112)                          {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 7.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 8.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_mc_pdg_pfp_hits_unmatched->Fill(leading_hits, 9.0); h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Fill(leading_hits, 9.0); }
			}

			if(leading_origin == "kBeamNeutrino")
			{
				h_leading_shower_origin->Fill(0);
				h_leading_shower_origin_hits->Fill(0.0, leading_hits);
				h_leading_shower_track_shower_beam->Fill(n_tracks, n_showers);
				if(leading_mc_pdg == 11)                         {h_leading_shower_origin_mc_pdg->Fill(0.0, 0.0);  h_leading_shower_mc_pdg->Fill(0.0); }
				if(leading_mc_pdg == -11)                        {h_leading_shower_origin_mc_pdg->Fill(0.0, 1.0);  h_leading_shower_mc_pdg->Fill(1.0); }
				if(leading_mc_pdg == 13)                         {h_leading_shower_origin_mc_pdg->Fill(0.0, 2.0);  h_leading_shower_mc_pdg->Fill(2.0); }
				if(leading_mc_pdg == -13)                        {h_leading_shower_origin_mc_pdg->Fill(0.0, 3.0);  h_leading_shower_mc_pdg->Fill(3.0); }
				if(leading_mc_pdg == 22)                         {h_leading_shower_origin_mc_pdg->Fill(0.0, 4.0);  h_leading_shower_mc_pdg->Fill(4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_origin_mc_pdg->Fill(0.0, 5.0);  h_leading_shower_mc_pdg->Fill(5.0); }
				if(leading_mc_pdg == 2212)                       {h_leading_shower_origin_mc_pdg->Fill(0.0, 6.0);  h_leading_shower_mc_pdg->Fill(6.0); }
				if(leading_mc_pdg == 2112)                       {h_leading_shower_origin_mc_pdg->Fill(0.0, 7.0);  h_leading_shower_mc_pdg->Fill(7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_origin_mc_pdg->Fill(0.0, 8.0); h_leading_shower_mc_pdg->Fill(8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_origin_mc_pdg->Fill(0.0, 9.0); h_leading_shower_mc_pdg->Fill(9.0); }
				if(tpco_cat == "nue_cc")       {h_leading_shower_origin_tpco_cat->Fill(0.0, 0.0); }
				if(tpco_cat == "nue_cc_mixed") {h_leading_shower_origin_tpco_cat->Fill(0.0, 1.0); }
				if(tpco_cat == "cosmic")       {h_leading_shower_origin_tpco_cat->Fill(0.0, 2.0); }
				if(tpco_cat == "nue_nc")       {h_leading_shower_origin_tpco_cat->Fill(0.0, 3.0); }
				if(tpco_cat == "numu_cc")      {h_leading_shower_origin_tpco_cat->Fill(0.0, 4.0); }
				if(tpco_cat == "other_mixed")  {h_leading_shower_origin_tpco_cat->Fill(0.0, 5.0); }
				if(tpco_cat == "unmatched")    {h_leading_shower_origin_tpco_cat->Fill(0.0, 6.0); }
			}
			if(leading_origin == "kCosmicRay")
			{
				h_leading_shower_origin->Fill(1);
				h_leading_shower_origin_hits->Fill(1.0, leading_hits);
				h_leading_shower_track_shower_cosmic->Fill(n_tracks, n_showers);
				if(leading_mc_pdg == 11)                         {h_leading_shower_origin_mc_pdg->Fill(1.0, 0.0);  h_leading_shower_mc_pdg->Fill(0.0); }
				if(leading_mc_pdg == -11)                        {h_leading_shower_origin_mc_pdg->Fill(1.0, 1.0);  h_leading_shower_mc_pdg->Fill(1.0); }
				if(leading_mc_pdg == 13)                         {h_leading_shower_origin_mc_pdg->Fill(1.0, 2.0);  h_leading_shower_mc_pdg->Fill(2.0); }
				if(leading_mc_pdg == -13)                        {h_leading_shower_origin_mc_pdg->Fill(1.0, 3.0);  h_leading_shower_mc_pdg->Fill(3.0); }
				if(leading_mc_pdg == 22)                         {h_leading_shower_origin_mc_pdg->Fill(1.0, 4.0);  h_leading_shower_mc_pdg->Fill(4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_origin_mc_pdg->Fill(1.0, 5.0);  h_leading_shower_mc_pdg->Fill(5.0); }
				if(leading_mc_pdg == 2212)                       {h_leading_shower_origin_mc_pdg->Fill(1.0, 6.0);  h_leading_shower_mc_pdg->Fill(6.0); }
				if(leading_mc_pdg == 2112)                       {h_leading_shower_origin_mc_pdg->Fill(1.0, 7.0);  h_leading_shower_mc_pdg->Fill(7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_origin_mc_pdg->Fill(1.0, 8.0); h_leading_shower_mc_pdg->Fill(8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_origin_mc_pdg->Fill(1.0, 9.0); h_leading_shower_mc_pdg->Fill(9.0); }
				if(tpco_cat == "nue_cc")       {h_leading_shower_origin_tpco_cat->Fill(1.0, 0.0); }
				if(tpco_cat == "nue_cc_mixed") {h_leading_shower_origin_tpco_cat->Fill(1.0, 1.0); }
				if(tpco_cat == "cosmic")       {h_leading_shower_origin_tpco_cat->Fill(1.0, 2.0); }
				if(tpco_cat == "nue_nc")       {h_leading_shower_origin_tpco_cat->Fill(1.0, 3.0); }
				if(tpco_cat == "numu_cc")      {h_leading_shower_origin_tpco_cat->Fill(1.0, 4.0); }
				if(tpco_cat == "other_mixed")  {h_leading_shower_origin_tpco_cat->Fill(1.0, 5.0); }
				if(tpco_cat == "unmatched")    {h_leading_shower_origin_tpco_cat->Fill(1.0, 6.0); }
			}
			if(leading_origin == "kUnknown")
			{
				h_leading_shower_origin->Fill(2);
				h_leading_shower_origin_hits->Fill(2.0, leading_hits);
				h_leading_shower_track_shower_unknown->Fill(n_tracks, n_showers);
				if(leading_mc_pdg == 11)                         {h_leading_shower_origin_mc_pdg->Fill(2.0, 0.0);  h_leading_shower_mc_pdg->Fill(0.0); }
				if(leading_mc_pdg == -11)                        {h_leading_shower_origin_mc_pdg->Fill(2.0, 1.0);  h_leading_shower_mc_pdg->Fill(1.0); }
				if(leading_mc_pdg == 13)                         {h_leading_shower_origin_mc_pdg->Fill(2.0, 2.0);  h_leading_shower_mc_pdg->Fill(2.0); }
				if(leading_mc_pdg == -13)                        {h_leading_shower_origin_mc_pdg->Fill(2.0, 3.0);  h_leading_shower_mc_pdg->Fill(3.0); }
				if(leading_mc_pdg == 22)                         {h_leading_shower_origin_mc_pdg->Fill(2.0, 4.0);  h_leading_shower_mc_pdg->Fill(4.0); }
				if(leading_mc_pdg == 211 || leading_mc_pdg == -211) {h_leading_shower_origin_mc_pdg->Fill(2.0, 5.0);  h_leading_shower_mc_pdg->Fill(5.0); }
				if(leading_mc_pdg == 2212)                       {h_leading_shower_origin_mc_pdg->Fill(2.0, 6.0);  h_leading_shower_mc_pdg->Fill(6.0); }
				if(leading_mc_pdg == 2112)                       {h_leading_shower_origin_mc_pdg->Fill(2.0, 7.0);  h_leading_shower_mc_pdg->Fill(7.0); }
				if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
				   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
				   leading_mc_pdg == -321)                       {h_leading_shower_origin_mc_pdg->Fill(2.0, 8.0); h_leading_shower_mc_pdg->Fill(8.0); }
				if(leading_mc_pdg == 0)                          {h_leading_shower_origin_mc_pdg->Fill(2.0, 9.0); h_leading_shower_mc_pdg->Fill(9.0); }
				if(tpco_cat == "nue_cc")       {h_leading_shower_origin_tpco_cat->Fill(2.0, 0.0); }
				if(tpco_cat == "nue_cc_mixed") {h_leading_shower_origin_tpco_cat->Fill(2.0, 1.0); }
				if(tpco_cat == "cosmic")       {h_leading_shower_origin_tpco_cat->Fill(2.0, 2.0); }
				if(tpco_cat == "nue_nc")       {h_leading_shower_origin_tpco_cat->Fill(2.0, 3.0); }
				if(tpco_cat == "numu_cc")      {h_leading_shower_origin_tpco_cat->Fill(2.0, 4.0); }
				if(tpco_cat == "other_mixed")  {h_leading_shower_origin_tpco_cat->Fill(2.0, 5.0); }
				if(tpco_cat == "unmatched")    {h_leading_shower_origin_tpco_cat->Fill(2.0, 6.0); }
			}

			std::cout << "\n " << std::endl;
			h_tpco_origin_unknown->Fill(origin_unknown);
			h_tpco_origin_cosmic->Fill(origin_cosmic);
			h_tpco_origin_neutrino->Fill(origin_beam);
			h_tpco_origin_neutrino_cosmic->Fill(origin_beam, origin_cosmic);
			h_nue_daughter_track_shower->Fill(n_tracks, n_showers);

		}//end looping tpc objects
	}//end looping events
	 //*******************
	TCanvas * origin_c1 = new TCanvas();
	origin_c1->cd();
	h_tpco_origin_unknown->GetYaxis()->SetTitleOffset(1.3);
	h_tpco_origin_unknown->GetXaxis()->SetTitle("NPFP Origin = Unknown");
	h_tpco_origin_unknown->GetYaxis()->SetTitle("NPFP per TPC Objects reco as Nue");
	h_tpco_origin_unknown->Draw();
	origin_c1->Print("tpco_origin_unknown.pdf");
	TCanvas * origin_c2 = new TCanvas();
	origin_c2->cd();
	h_tpco_origin_cosmic->GetYaxis()->SetTitleOffset(1.3);
	h_tpco_origin_cosmic->GetXaxis()->SetTitle("NPFP Origin = Cosmic");
	h_tpco_origin_cosmic->GetYaxis()->SetTitle("NPFP per TPC Objects reco as Nue");
	h_tpco_origin_cosmic->Draw();
	origin_c2->Print("tpco_origin_cosmic.pdf");
	TCanvas * origin_c3 = new TCanvas();
	origin_c3->cd();
	h_tpco_origin_neutrino->GetYaxis()->SetTitleOffset(1.3);
	h_tpco_origin_neutrino->GetXaxis()->SetTitle("NPFP Origin = Neutrino");
	h_tpco_origin_neutrino->GetYaxis()->SetTitle("NPFP per TPC Objects reco as Nue");
	h_tpco_origin_neutrino->Draw();
	origin_c3->Print("tpco_origin_neutrino.pdf");
	TCanvas * origin_c4 = new TCanvas();
	origin_c4->cd();
	h_tpco_origin_neutrino_cosmic->GetXaxis()->SetTitle("NPFP Origin = Neutrino");
	h_tpco_origin_neutrino_cosmic->GetYaxis()->SetTitle("NPFP Origin = Cosmic");
	h_tpco_origin_neutrino_cosmic->Draw("colz");
	origin_c4->Print("tpco_origin_neutrino_cosmic.pdf");
	TCanvas * origin_c5 = new TCanvas();
	origin_c5->cd();
	h_nue_daughter_track_shower->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_nue_daughter_track_shower->GetYaxis()->SetTitle("N Showers per TPC Object");
	h_nue_daughter_track_shower->SetStats(kFALSE);
	h_nue_daughter_track_shower->Draw("colz");
	origin_c5->Print("nue_daughter_track_shower.pdf");
	TCanvas * origin_c6 = new TCanvas();
	origin_c6->cd();
	h_leading_shower_track_shower_beam->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_leading_shower_track_shower_beam->GetYaxis()->SetTitle("N Showers per TPC Object");
	h_leading_shower_track_shower_beam->SetStats(kFALSE);
	h_leading_shower_track_shower_beam->Draw("colz");
	origin_c6->Print("leading_shower_track_shower_beam.pdf");
	TCanvas * origin_c7 = new TCanvas();
	origin_c7->cd();
	h_leading_shower_track_shower_cosmic->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_leading_shower_track_shower_cosmic->GetYaxis()->SetTitle("N Showers per TPC Object");
	h_leading_shower_track_shower_cosmic->SetStats(kFALSE);
	h_leading_shower_track_shower_cosmic->Draw("colz");
	origin_c7->Print("leading_shower_track_shower_cosmic.pdf");
	TCanvas * origin_c8 = new TCanvas();
	origin_c8->cd();
	h_leading_shower_track_shower_unknown->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_leading_shower_track_shower_unknown->GetYaxis()->SetTitle("N Showers per TPC Object");
	h_leading_shower_track_shower_unknown->SetStats(kFALSE);
	h_leading_shower_track_shower_unknown->Draw("colz");
	origin_c8->Print("leading_shower_track_shower_unknown.pdf");

	TCanvas * track_c1 = new TCanvas();
	track_c1->cd();
	h_tracks_tpco_cat_nue_cc->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_tracks_tpco_cat_nue_cc->Draw();
	track_c1->Print("tpc_object_category_tracks_nue_cc.pdf");
	TCanvas * track_c2 = new TCanvas();
	track_c2->cd();
	h_tracks_tpco_cat_nue_cc_mixed->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_tracks_tpco_cat_nue_cc_mixed->Draw();
	track_c2->Print("tpc_object_category_tracks_nue_cc_mixed.pdf");
	TCanvas * track_c3 = new TCanvas();
	track_c3->cd();
	h_tracks_tpco_cat_cosmic->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_tracks_tpco_cat_cosmic->Draw();
	track_c3->Print("tpc_object_category_tracks_cosmic.pdf");
	TCanvas * track_c4 = new TCanvas();
	track_c4->cd();
	h_tracks_tpco_cat_nue_nc->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_tracks_tpco_cat_nue_nc->Draw();
	track_c4->Print("tpc_object_category_tracks_nue_nc.pdf");
	TCanvas * track_c5 = new TCanvas();
	track_c5->cd();
	h_tracks_tpco_cat_numu->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_tracks_tpco_cat_numu->Draw();
	track_c5->Print("tpc_object_category_tracks_numu.pdf");
	TCanvas * track_c6 = new TCanvas();
	track_c6->cd();
	h_tracks_tpco_cat_other_mixed->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_tracks_tpco_cat_other_mixed->Draw();
	track_c6->Print("tpc_object_category_tracks_other_mixed.pdf");
	TCanvas * track_c7 = new TCanvas();
	track_c7->cd();
	h_tracks_tpco_cat_unmatched->GetXaxis()->SetTitle("N Tracks per TPC Object");
	h_tracks_tpco_cat_unmatched->Draw();
	track_c7->Print("tpc_object_category_tracks_unmatched.pdf");

	TCanvas * track_c8 = new TCanvas();
	THStack * hstack_tracks_tpco_cat = new THStack("hstack_tracks_tpco_cat","hstack_tracks_tpco_cat");
	hstack_tracks_tpco_cat->Add(h_tracks_tpco_cat_nue_cc);
	hstack_tracks_tpco_cat->Add(h_tracks_tpco_cat_nue_cc_mixed);
	hstack_tracks_tpco_cat->Add(h_tracks_tpco_cat_cosmic);
	hstack_tracks_tpco_cat->Add(h_tracks_tpco_cat_nue_nc);
	hstack_tracks_tpco_cat->Add(h_tracks_tpco_cat_numu);
	hstack_tracks_tpco_cat->Add(h_tracks_tpco_cat_other_mixed);
	hstack_tracks_tpco_cat->Add(h_tracks_tpco_cat_unmatched);

	h_tracks_tpco_cat_nue_cc->SetFillColor(30);
	h_tracks_tpco_cat_nue_cc_mixed->SetFillColor(38);
	h_tracks_tpco_cat_cosmic->SetFillColor(39);
	h_tracks_tpco_cat_nue_nc->SetFillColor(46);
	h_tracks_tpco_cat_numu->SetFillColor(28);
	h_tracks_tpco_cat_other_mixed->SetFillColor(42);
	h_tracks_tpco_cat_unmatched->SetFillColor(12);

	h_tracks_tpco_cat_nue_cc->SetStats(kFALSE);
	h_tracks_tpco_cat_nue_cc_mixed->SetStats(kFALSE);
	h_tracks_tpco_cat_cosmic->SetStats(kFALSE);
	h_tracks_tpco_cat_nue_nc->SetStats(kFALSE);
	h_tracks_tpco_cat_numu->SetStats(kFALSE);
	h_tracks_tpco_cat_other_mixed->SetStats(kFALSE);
	h_tracks_tpco_cat_unmatched->SetStats(kFALSE);

	hstack_tracks_tpco_cat->Draw();

	hstack_tracks_tpco_cat->GetYaxis()->SetTitleOffset(1.2);
	hstack_tracks_tpco_cat->GetXaxis()->SetTitleOffset(0.9);
	//hstack_tracks_tpco_cat->GetXaxis()->SetTitleSize(axis_title_size);
	//hstack_tracks_tpco_cat->GetYaxis()->SetTitleSize(axis_title_size);
	hstack_tracks_tpco_cat->GetXaxis()->SetTitle("N Tracks per Nue Candidate");
	hstack_tracks_tpco_cat->GetYaxis()->SetTitle("Counts");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_track = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_track->AddEntry(h_tracks_tpco_cat_nue_cc,         "Nue CC", "f");
	leg_track->AddEntry(h_tracks_tpco_cat_nue_cc_mixed,   "Nue CC Mixed", "f");
	leg_track->AddEntry(h_tracks_tpco_cat_cosmic,         "Cosmic", "f");
	leg_track->AddEntry(h_tracks_tpco_cat_nue_nc,         "Nue NC", "f");
	leg_track->AddEntry(h_tracks_tpco_cat_numu,           "Numu", "f");
	leg_track->AddEntry(h_tracks_tpco_cat_other_mixed,    "Other Mixed", "f");
	leg_track->AddEntry(h_tracks_tpco_cat_unmatched,      "Unmatched", "f");
	leg_track->Draw();

	track_c8->Print("stack_tracks_tpco_cat.pdf");

	TCanvas * shower_c1 = new TCanvas();
	shower_c1->cd();
	h_showers_tpco_cat_nue_cc->GetXaxis()->SetTitle("N Showers per TPC Object");
	h_showers_tpco_cat_nue_cc->Draw();
	shower_c1->Print("tpc_object_category_showers_nue_cc.pdf");
	TCanvas * shower_c2 = new TCanvas();
	shower_c2->cd();
	h_showers_tpco_cat_nue_cc_mixed->GetXaxis()->SetTitle("N Showers per TPC Object");
	h_showers_tpco_cat_nue_cc_mixed->Draw();
	shower_c2->Print("tpc_object_category_showers_nue_cc_mixed.pdf");
	TCanvas * shower_c3 = new TCanvas();
	shower_c3->cd();
	h_showers_tpco_cat_cosmic->GetXaxis()->SetTitle("N Showers per TPC Object");
	h_showers_tpco_cat_cosmic->Draw();
	shower_c3->Print("tpc_object_category_showers_cosmic.pdf");
	TCanvas * shower_c4 = new TCanvas();
	shower_c4->cd();
	h_showers_tpco_cat_nue_nc->GetXaxis()->SetTitle("N Showers per TPC Object");
	h_showers_tpco_cat_nue_nc->Draw();
	shower_c4->Print("tpc_object_category_showers_nue_nc.pdf");
	TCanvas * shower_c5 = new TCanvas();
	shower_c5->cd();
	h_showers_tpco_cat_numu->GetXaxis()->SetTitle("N Showers per TPC Object");
	h_showers_tpco_cat_numu->Draw();
	shower_c5->Print("tpc_object_category_showers_numu.pdf");
	TCanvas * shower_c6 = new TCanvas();
	shower_c6->cd();
	h_showers_tpco_cat_other_mixed->GetXaxis()->SetTitle("N Showers per TPC Object");
	h_showers_tpco_cat_other_mixed->Draw();
	shower_c6->Print("tpc_object_category_showers_other_mixed.pdf");
	TCanvas * shower_c7 = new TCanvas();
	shower_c7->cd();
	h_showers_tpco_cat_unmatched->GetXaxis()->SetTitle("N Showers per TPC Object");
	h_showers_tpco_cat_unmatched->Draw();
	shower_c7->Print("tpc_object_category_showers_unmatched.pdf");

	TCanvas * shower_c8 = new TCanvas();
	THStack * hstack_showers_tpco_cat = new THStack("hstack_showers_tpco_cat","hstack_showers_tpco_cat");
	hstack_showers_tpco_cat->Add(h_showers_tpco_cat_nue_cc);
	hstack_showers_tpco_cat->Add(h_showers_tpco_cat_nue_cc_mixed);
	hstack_showers_tpco_cat->Add(h_showers_tpco_cat_cosmic);
	hstack_showers_tpco_cat->Add(h_showers_tpco_cat_nue_nc);
	hstack_showers_tpco_cat->Add(h_showers_tpco_cat_numu);
	hstack_showers_tpco_cat->Add(h_showers_tpco_cat_other_mixed);
	hstack_showers_tpco_cat->Add(h_showers_tpco_cat_unmatched);

	h_showers_tpco_cat_nue_cc->SetFillColor(30);
	h_showers_tpco_cat_nue_cc_mixed->SetFillColor(38);
	h_showers_tpco_cat_cosmic->SetFillColor(39);
	h_showers_tpco_cat_nue_nc->SetFillColor(46);
	h_showers_tpco_cat_numu->SetFillColor(28);
	h_showers_tpco_cat_other_mixed->SetFillColor(42);
	h_showers_tpco_cat_unmatched->SetFillColor(12);

	h_showers_tpco_cat_nue_cc->SetStats(kFALSE);
	h_showers_tpco_cat_nue_cc_mixed->SetStats(kFALSE);
	h_showers_tpco_cat_cosmic->SetStats(kFALSE);
	h_showers_tpco_cat_nue_nc->SetStats(kFALSE);
	h_showers_tpco_cat_numu->SetStats(kFALSE);
	h_showers_tpco_cat_other_mixed->SetStats(kFALSE);
	h_showers_tpco_cat_unmatched->SetStats(kFALSE);

	hstack_showers_tpco_cat->Draw();

	hstack_showers_tpco_cat->GetYaxis()->SetTitleOffset(1.2);
	hstack_showers_tpco_cat->GetXaxis()->SetTitleOffset(0.9);
	//hstack_showers_tpco_cat->GetXaxis()->SetTitleSize(axis_title_size);
	//hstack_showers_tpco_cat->GetYaxis()->SetTitleSize(axis_title_size);
	hstack_showers_tpco_cat->GetXaxis()->SetTitle("N Showers per Nue Candidate");
	hstack_showers_tpco_cat->GetYaxis()->SetTitle("Counts");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_shower = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_shower->AddEntry(h_showers_tpco_cat_nue_cc,         "Nue CC", "f");
	leg_shower->AddEntry(h_showers_tpco_cat_nue_cc_mixed,   "Nue CC Mixed", "f");
	leg_shower->AddEntry(h_showers_tpco_cat_cosmic,         "Cosmic", "f");
	leg_shower->AddEntry(h_showers_tpco_cat_nue_nc,         "Nue NC", "f");
	leg_shower->AddEntry(h_showers_tpco_cat_numu,           "Numu", "f");
	leg_shower->AddEntry(h_showers_tpco_cat_other_mixed,    "Other Mixed", "f");
	leg_shower->AddEntry(h_showers_tpco_cat_unmatched,      "Unmatched", "f");
	leg_shower->Draw();

	shower_c8->Print("stack_showers_tpco_cat.pdf");


	TCanvas * shwr_to_vtx_c1 = new TCanvas();
	shwr_to_vtx_c1->cd();
	h_shwr_to_vtx_beam_electron->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_beam_electron->Draw();
	shwr_to_vtx_c1->Print("shwr_to_vtx_beam_electron.pdf");
	TCanvas * shwr_to_vtx_c2 = new TCanvas();
	shwr_to_vtx_c2->cd();
	h_shwr_to_vtx_beam_photon->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_beam_photon->Draw();
	shwr_to_vtx_c2->Print("shwr_to_vtx_beam_photon.pdf");
	TCanvas * shwr_to_vtx_c3 = new TCanvas();
	shwr_to_vtx_c3->cd();
	h_shwr_to_vtx_cosmic_electron->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_cosmic_electron->Draw();
	shwr_to_vtx_c3->Print("shwr_to_vtx_cosmic_electron.pdf");
	TCanvas * shwr_to_vtx_c4 = new TCanvas();
	shwr_to_vtx_c4->cd();
	h_shwr_to_vtx_cosmic_photon->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_cosmic_photon->Draw();
	shwr_to_vtx_c4->Print("shwr_to_vtx_cosmic_photon.pdf");
	TCanvas * shwr_to_vtx_c5 = new TCanvas();
	shwr_to_vtx_c5->cd();
	h_shwr_to_vtx_unmatched->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_unmatched->Draw();
	shwr_to_vtx_c5->Print("shwr_to_vtx_unmatched.pdf");
	TCanvas * shwr_to_vtx_overlay_c1 = new TCanvas();
	shwr_to_vtx_overlay_c1->cd();
	h_shwr_to_vtx_beam_electron->SetLineColor(30);
	h_shwr_to_vtx_beam_photon->SetLineColor(38);
	h_shwr_to_vtx_cosmic_electron->SetLineColor(39);
	h_shwr_to_vtx_cosmic_photon->SetLineColor(46);
	h_shwr_to_vtx_unmatched->SetLineColor(28);
	h_shwr_to_vtx_beam_electron->SetStats(kFALSE);
	h_shwr_to_vtx_beam_photon->SetStats(kFALSE);
	h_shwr_to_vtx_cosmic_electron->SetStats(kFALSE);
	h_shwr_to_vtx_cosmic_photon->SetStats(kFALSE);
	h_shwr_to_vtx_unmatched->SetStats(kFALSE);
	h_shwr_to_vtx_beam_electron->Draw();
	h_shwr_to_vtx_beam_photon->Draw("SAME");
	h_shwr_to_vtx_cosmic_electron->Draw("SAME");
	h_shwr_to_vtx_cosmic_photon->Draw("SAME");
	h_shwr_to_vtx_unmatched->Draw("SAME");
	TLegend * leg_shwr_to_vtx = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_shwr_to_vtx->AddEntry(h_shwr_to_vtx_beam_electron,         "Beam Electron", "f");
	leg_shwr_to_vtx->AddEntry(h_shwr_to_vtx_beam_photon,           "Beam Photon", "f");
	leg_shwr_to_vtx->AddEntry(h_shwr_to_vtx_cosmic_electron,       "Cosmic Electron", "f");
	leg_shwr_to_vtx->AddEntry(h_shwr_to_vtx_cosmic_photon,         "Cosmic Photon", "f");
	leg_shwr_to_vtx->AddEntry(h_shwr_to_vtx_unmatched,             "Unmatched", "f");
	leg_shwr_to_vtx->Draw();
	shwr_to_vtx_overlay_c1->Print("shwr_to_vtx_overlay.pdf");

	TCanvas * shwr_dedx_c1 = new TCanvas();
	shwr_dedx_c1->cd();
	h_shwr_dedx_nue_cc->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_nue_cc->Draw();
	shwr_dedx_c1->Print("shwr_dedx_nue_cc.pdf");
	TCanvas * shwr_dedx_c2 = new TCanvas();
	shwr_dedx_c2->cd();
	h_shwr_dedx_nue_cc_mixed->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_nue_cc_mixed->Draw();
	shwr_dedx_c2->Print("shwr_dedx_nue_cc_mixed.pdf");
	TCanvas * shwr_dedx_c3 = new TCanvas();
	shwr_dedx_c3->cd();
	h_shwr_dedx_nue_nc->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_nue_nc->Draw();
	shwr_dedx_c3->Print("shwr_dedx_nue_nc.pdf");
	TCanvas * shwr_dedx_c4 = new TCanvas();
	shwr_dedx_c4->cd();
	h_shwr_dedx_numu_cc->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_numu_cc->Draw();
	shwr_dedx_c4->Print("shwr_dedx_numu_cc.pdf");
	TCanvas * shwr_dedx_c5 = new TCanvas();
	shwr_dedx_c5->cd();
	h_shwr_dedx_numu_cc_mixed->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_numu_cc_mixed->Draw();
	shwr_dedx_c5->Print("shwr_dedx_numu_cc_mixed.pdf");
	TCanvas * shwr_dedx_c6 = new TCanvas();
	shwr_dedx_c6->cd();
	h_shwr_dedx_numu_nc->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_numu_nc->Draw();
	shwr_dedx_c6->Print("shwr_dedx_numu_nc.pdf");
	TCanvas * shwr_dedx_c7 = new TCanvas();
	shwr_dedx_c7->cd();
	h_shwr_dedx_cosmic->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_cosmic->Draw();
	shwr_dedx_c7->Print("shwr_dedx_cosmic.pdf");
	TCanvas * shwr_dedx_c8 = new TCanvas();
	shwr_dedx_c8->cd();
	h_shwr_dedx_other_mixed->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_other_mixed->Draw();
	shwr_dedx_c8->Print("shwr_dedx_other_mixed.pdf");
	TCanvas * shwr_dedx_c9 = new TCanvas();
	shwr_dedx_c9->cd();
	h_shwr_dedx_unmatched->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_unmatched->Draw();
	shwr_dedx_c9->Print("shwr_dedx_unmatched.pdf");

	TCanvas * shwr_dedx_theta_c1 = new TCanvas();
	shwr_dedx_theta_c1->cd();
	h_shwr_dedx_theta_nue_cc->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_theta_nue_cc->GetYaxis()->SetTitle("Cos(#theta)");
	h_shwr_dedx_theta_nue_cc->Draw("colz");
	shwr_dedx_theta_c1->Print("shwr_dedx_theta_nue_cc.pdf");
	TCanvas * shwr_dedx_theta_c2 = new TCanvas();
	shwr_dedx_theta_c2->cd();
	h_shwr_dedx_theta_nue_cc_mixed->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_theta_nue_cc_mixed->GetYaxis()->SetTitle("Cos(#theta)");
	h_shwr_dedx_theta_nue_cc_mixed->Draw("colz");
	shwr_dedx_theta_c2->Print("shwr_dedx_theta_nue_cc_mixed.pdf");
	TCanvas * shwr_dedx_theta_c3 = new TCanvas();
	shwr_dedx_theta_c3->cd();
	h_shwr_dedx_theta_nue_nc->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_theta_nue_nc->GetYaxis()->SetTitle("Cos(#theta)");
	h_shwr_dedx_theta_nue_nc->Draw("colz");
	shwr_dedx_theta_c3->Print("shwr_dedx_theta_nue_nc.pdf");
	TCanvas * shwr_dedx_theta_c4 = new TCanvas();
	shwr_dedx_theta_c4->cd();
	h_shwr_dedx_theta_numu_cc->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_theta_numu_cc->GetYaxis()->SetTitle("Cos(#theta)");
	h_shwr_dedx_theta_numu_cc->Draw("colz");
	shwr_dedx_theta_c4->Print("shwr_dedx_theta_numu_cc.pdf");
	TCanvas * shwr_dedx_theta_c5 = new TCanvas();
	shwr_dedx_theta_c5->cd();
	h_shwr_dedx_theta_numu_cc_mixed->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_theta_numu_cc_mixed->GetYaxis()->SetTitle("Cos(#theta)");
	h_shwr_dedx_theta_numu_cc_mixed->Draw("colz");
	shwr_dedx_theta_c5->Print("shwr_dedx_theta_numu_cc_mixed.pdf");
	TCanvas * shwr_dedx_theta_c6 = new TCanvas();
	shwr_dedx_theta_c6->cd();
	h_shwr_dedx_theta_numu_nc->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_theta_numu_nc->GetYaxis()->SetTitle("Cos(#theta)");
	h_shwr_dedx_theta_numu_nc->Draw("colz");
	shwr_dedx_theta_c6->Print("shwr_dedx_theta_numu_nc.pdf");
	TCanvas * shwr_dedx_theta_c7 = new TCanvas();
	shwr_dedx_theta_c7->cd();
	h_shwr_dedx_theta_cosmic->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_theta_cosmic->GetYaxis()->SetTitle("Cos(#theta)");
	h_shwr_dedx_theta_cosmic->Draw("colz");
	shwr_dedx_theta_c7->Print("shwr_dedx_theta_cosmic.pdf");
	TCanvas * shwr_dedx_theta_c8 = new TCanvas();
	shwr_dedx_theta_c8->cd();
	h_shwr_dedx_theta_other_mixed->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_theta_other_mixed->GetYaxis()->SetTitle("Cos(#theta)");
	h_shwr_dedx_theta_other_mixed->Draw("colz");
	shwr_dedx_theta_c8->Print("shwr_dedx_theta_other_mixed.pdf");
	TCanvas * shwr_dedx_theta_c9 = new TCanvas();
	shwr_dedx_theta_c9->cd();
	h_shwr_dedx_theta_unmatched->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
	h_shwr_dedx_theta_unmatched->GetYaxis()->SetTitle("Cos(#theta)");
	h_shwr_dedx_theta_unmatched->Draw("colz");
	shwr_dedx_theta_c9->Print("shwr_dedx_theta_unmatched.pdf");


	TCanvas * shwr_dedx_stack_c1 = new TCanvas();
	shwr_dedx_stack_c1->cd();
	h_shwr_dedx_nue_cc->SetFillColor(30);
	h_shwr_dedx_nue_cc_mixed->SetFillColor(38);
	h_shwr_dedx_nue_nc->SetFillColor(28);
	h_shwr_dedx_numu_cc->SetFillColor(36);
	h_shwr_dedx_numu_cc_mixed->SetFillColor(39);
	h_shwr_dedx_numu_nc->SetFillColor(46);
	h_shwr_dedx_cosmic->SetFillColor(25);
	h_shwr_dedx_other_mixed->SetFillColor(42);
	h_shwr_dedx_unmatched->SetFillColor(12);
	h_shwr_dedx_nue_cc->SetStats(kFALSE);
	h_shwr_dedx_nue_cc_mixed->SetStats(kFALSE);
	h_shwr_dedx_nue_nc->SetStats(kFALSE);
	h_shwr_dedx_numu_cc->SetStats(kFALSE);
	h_shwr_dedx_numu_cc_mixed->SetStats(kFALSE);
	h_shwr_dedx_numu_nc->SetStats(kFALSE);
	h_shwr_dedx_cosmic->SetStats(kFALSE);
	h_shwr_dedx_other_mixed->SetStats(kFALSE);
	h_shwr_dedx_unmatched->SetStats(kFALSE);
	THStack * shwr_dedx_stack = new THStack();
	shwr_dedx_stack->Add(h_shwr_dedx_nue_cc);
	shwr_dedx_stack->Add(h_shwr_dedx_nue_cc_mixed);
	shwr_dedx_stack->Add(h_shwr_dedx_nue_nc);
	shwr_dedx_stack->Add(h_shwr_dedx_numu_cc);
	shwr_dedx_stack->Add(h_shwr_dedx_numu_cc_mixed);
	shwr_dedx_stack->Add(h_shwr_dedx_numu_nc);
	shwr_dedx_stack->Add(h_shwr_dedx_cosmic);
	shwr_dedx_stack->Add(h_shwr_dedx_other_mixed);
	shwr_dedx_stack->Add(h_shwr_dedx_unmatched);
	shwr_dedx_stack->Draw();
	shwr_dedx_stack->GetXaxis()->SetTitle("dE/dx [MeV / cm]");
	TLegend * leg_shwr_dedx = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_shwr_dedx->AddEntry(h_shwr_dedx_nue_cc,             "Nue CC", "f");
	leg_shwr_dedx->AddEntry(h_shwr_dedx_nue_cc_mixed,       "Nue CC Mixed", "f");
	leg_shwr_dedx->AddEntry(h_shwr_dedx_nue_nc,             "Nue NC", "f");
	leg_shwr_dedx->AddEntry(h_shwr_dedx_numu_cc,            "Numu CC", "f");
	leg_shwr_dedx->AddEntry(h_shwr_dedx_numu_cc_mixed,      "Numu CC Mixed", "f");
	leg_shwr_dedx->AddEntry(h_shwr_dedx_numu_nc,            "Numu NC", "f");
	leg_shwr_dedx->AddEntry(h_shwr_dedx_cosmic,             "Cosmic", "f");
	leg_shwr_dedx->AddEntry(h_shwr_dedx_other_mixed,        "Other Mixed", "f");
	leg_shwr_dedx->AddEntry(h_shwr_dedx_unmatched,          "Unmatched", "f");
	leg_shwr_dedx->Draw();
	shwr_dedx_stack_c1->Print("shwr_dedx_stack.pdf");


	TCanvas * shwr_to_vtx_tracks_c1 = new TCanvas();
	shwr_to_vtx_tracks_c1->cd();
	h_shwr_to_vtx_beam_electron_tracks->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_beam_electron_tracks->GetYaxis()->SetTitle("Number of Reco Tracks");
	h_shwr_to_vtx_beam_electron_tracks->Draw("colz");
	shwr_to_vtx_tracks_c1->Print("shwr_to_vtx_beam_electron_tracks.pdf");
	TCanvas * shwr_to_vtx_tracks_c2 = new TCanvas();
	shwr_to_vtx_tracks_c2->cd();
	h_shwr_to_vtx_beam_photon_tracks->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_beam_photon_tracks->GetYaxis()->SetTitle("Number of Reco Tracks");
	h_shwr_to_vtx_beam_photon_tracks->Draw("colz");
	shwr_to_vtx_tracks_c2->Print("shwr_to_vtx_beam_photon_tracks.pdf");
	TCanvas * shwr_to_vtx_tracks_c3 = new TCanvas();
	shwr_to_vtx_tracks_c3->cd();
	h_shwr_to_vtx_cosmic_electron_tracks->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_cosmic_electron_tracks->GetYaxis()->SetTitle("Number of Reco Tracks");
	h_shwr_to_vtx_cosmic_electron_tracks->Draw("colz");
	shwr_to_vtx_tracks_c3->Print("shwr_to_vtx_cosmic_electron_tracks.pdf");
	TCanvas * shwr_to_vtx_tracks_c4 = new TCanvas();
	shwr_to_vtx_tracks_c4->cd();
	h_shwr_to_vtx_cosmic_photon_tracks->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_cosmic_photon_tracks->GetYaxis()->SetTitle("Number of Reco Tracks");
	h_shwr_to_vtx_cosmic_photon_tracks->Draw("colz");
	shwr_to_vtx_tracks_c4->Print("shwr_to_vtx_cosmic_photon_tracks.pdf");
	TCanvas * shwr_to_vtx_tracks_c5 = new TCanvas();
	shwr_to_vtx_tracks_c5->cd();
	h_shwr_to_vtx_unmatched_tracks->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_unmatched_tracks->GetYaxis()->SetTitle("Number of Reco Tracks");
	h_shwr_to_vtx_unmatched_tracks->Draw("colz");
	shwr_to_vtx_tracks_c5->Print("shwr_to_vtx_unmatched_tracks.pdf");

	TCanvas * shwr_to_vtx_showers_c1 = new TCanvas();
	shwr_to_vtx_showers_c1->cd();
	h_shwr_to_vtx_beam_electron_showers->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_beam_electron_showers->GetYaxis()->SetTitle("Number of Reco Showers");
	h_shwr_to_vtx_beam_electron_showers->Draw("colz");
	shwr_to_vtx_showers_c1->Print("shwr_to_vtx_beam_electron_showers.pdf");
	TCanvas * shwr_to_vtx_showers_c2 = new TCanvas();
	shwr_to_vtx_showers_c2->cd();
	h_shwr_to_vtx_beam_photon_showers->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_beam_photon_showers->GetYaxis()->SetTitle("Number of Reco Showers");
	h_shwr_to_vtx_beam_photon_showers->Draw("colz");
	shwr_to_vtx_showers_c2->Print("shwr_to_vtx_beam_photon_showers.pdf");
	TCanvas * shwr_to_vtx_showers_c3 = new TCanvas();
	shwr_to_vtx_showers_c3->cd();
	h_shwr_to_vtx_cosmic_electron_showers->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_cosmic_electron_showers->GetYaxis()->SetTitle("Number of Reco Showers");
	h_shwr_to_vtx_cosmic_electron_showers->Draw("colz");
	shwr_to_vtx_showers_c3->Print("shwr_to_vtx_cosmic_electron_showers.pdf");
	TCanvas * shwr_to_vtx_showers_c4 = new TCanvas();
	shwr_to_vtx_showers_c4->cd();
	h_shwr_to_vtx_cosmic_photon_showers->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_cosmic_photon_showers->GetYaxis()->SetTitle("Number of Reco Showers");
	h_shwr_to_vtx_cosmic_photon_showers->Draw("colz");
	shwr_to_vtx_showers_c4->Print("shwr_to_vtx_cosmic_photon_showers.pdf");
	TCanvas * shwr_to_vtx_showers_c5 = new TCanvas();
	shwr_to_vtx_showers_c5->cd();
	h_shwr_to_vtx_unmatched_showers->GetXaxis()->SetTitle("PFP Shower Vtx to Neutrino Candidate Vertex [cm]");
	h_shwr_to_vtx_unmatched_showers->GetYaxis()->SetTitle("Number of Reco Showers");
	h_shwr_to_vtx_unmatched_showers->Draw("colz");
	shwr_to_vtx_showers_c5->Print("shwr_to_vtx_unmatched_showers.pdf");

	TCanvas * leading_open_angle_c1 = new TCanvas();
	leading_open_angle_c1->cd();
	h_shwr_open_angle_nue_cc->GetXaxis()->SetTitle("Leading Shower Open Angle");
	h_shwr_open_angle_nue_cc->Draw();
	leading_open_angle_c1->Print("leading_shower_open_angle_nue_cc.pdf");
	TCanvas * leading_open_angle_c2 = new TCanvas();
	leading_open_angle_c2->cd();
	h_shwr_open_angle_nue_cc_mixed->GetXaxis()->SetTitle("Leading Shower Open Angle");
	h_shwr_open_angle_nue_cc_mixed->Draw();
	leading_open_angle_c2->Print("leading_shower_open_angle_nue_cc_mixed.pdf");
	TCanvas * leading_open_angle_c3 = new TCanvas();
	leading_open_angle_c3->cd();
	h_shwr_open_angle_numu_cc->GetXaxis()->SetTitle("Leading Shower Open Angle");
	h_shwr_open_angle_numu_cc->Draw();
	leading_open_angle_c3->Print("leading_shower_open_angle_numu_cc.pdf");
	TCanvas * leading_open_angle_c4 = new TCanvas();
	leading_open_angle_c4->cd();
	h_shwr_open_angle_numu_cc_mixed->GetXaxis()->SetTitle("Leading Shower Open Angle");
	h_shwr_open_angle_numu_cc->Draw();
	leading_open_angle_c4->Print("leading_shower_open_angle_numu_cc.pdf");
	TCanvas * leading_open_angle_c5 = new TCanvas();
	leading_open_angle_c5->cd();
	h_shwr_open_angle_numu_nc->GetXaxis()->SetTitle("Leading Shower Open Angle");
	h_shwr_open_angle_numu_nc->Draw();
	leading_open_angle_c5->Print("leading_shower_open_angle_numu_nc.pdf");
	TCanvas * leading_open_angle_c6 = new TCanvas();
	leading_open_angle_c6->cd();
	h_shwr_open_angle_cosmic->GetXaxis()->SetTitle("Leading Shower Open Angle");
	h_shwr_open_angle_cosmic->Draw();
	leading_open_angle_c6->Print("leading_shower_open_angle_cosmic.pdf");
	TCanvas * leading_open_angle_stack = new TCanvas();
	leading_open_angle_stack->cd();
	h_shwr_open_angle_nue_cc->SetStats(kFALSE);
	h_shwr_open_angle_nue_cc_mixed->SetStats(kFALSE);
	h_shwr_open_angle_numu_cc->SetStats(kFALSE);
	h_shwr_open_angle_numu_cc_mixed->SetStats(kFALSE);
	h_shwr_open_angle_numu_nc->SetStats(kFALSE);
	h_shwr_open_angle_cosmic->SetStats(kFALSE);
	h_shwr_open_angle_nue_cc->SetFillColor(30);
	h_shwr_open_angle_nue_cc_mixed->SetFillColor(38);
	h_shwr_open_angle_numu_cc->SetFillColor(28);
	h_shwr_open_angle_numu_cc_mixed->SetFillColor(42);
	h_shwr_open_angle_numu_nc->SetFillColor(35);
	h_shwr_open_angle_cosmic->SetFillColor(39);
	THStack * open_angle_stack = new THStack();
	open_angle_stack->Add(h_shwr_open_angle_nue_cc);
	open_angle_stack->Add(h_shwr_open_angle_nue_cc_mixed);
	open_angle_stack->Add(h_shwr_open_angle_numu_cc);
	open_angle_stack->Add(h_shwr_open_angle_numu_cc_mixed);
	open_angle_stack->Add(h_shwr_open_angle_numu_nc);
	open_angle_stack->Add(h_shwr_open_angle_cosmic);
	open_angle_stack->Draw();

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_open_angle = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_open_angle->AddEntry(h_shwr_open_angle_nue_cc,         "Nue CC", "f");
	leg_open_angle->AddEntry(h_shwr_open_angle_nue_cc_mixed,   "Nue CC Mixed", "f");
	leg_open_angle->AddEntry(h_shwr_open_angle_cosmic,         "Cosmic", "f");
	leg_open_angle->AddEntry(h_shwr_open_angle_numu_cc,        "Numu CC", "f");
	leg_open_angle->AddEntry(h_shwr_open_angle_numu_cc_mixed,  "Numu Mixed", "f");
	leg_open_angle->AddEntry(h_shwr_open_angle_numu_nc,        "Numu NC", "f");
	leg_open_angle->Draw();
	leading_open_angle_stack->Print("leading_shower_open_angle_stack.pdf");

	TH1D * h_nue_daughter_origin = new TH1D("h_nue_daughter_origin", "h_nue_daughter_origin", 3, 0, 3);
	TH1D * h_nue_daughter_pfp_pdg = new TH1D("h_nue_daughter_pfp_pdg", "h_nue_daughter_pfp_pdg", 3, 0, 3);
	TH1D * h_nue_daughter_mc_pdg = new TH1D("h_nue_daughter_mc_pdg", "h_nue_daughter_mc_pdg", 10, 0, 10);
	TH1D * h_nue_daughter_pfp_hits = new TH1D("h_nue_daughter_pfp_hits", "h_nue_daughter_pfp_hits", 20, 0, 3000);
	TH2I * h_nue_daughter_mc_pfp_pdg = new TH2I("h_nue_daughter_mc_pfp_pdg", "h_nue_daughter_mc_pfp_pdg", 10, 0, 10, 2, 0, 2);
	TH2I * h_nue_daughter_shower_mc_pdg_pfp_hits = new TH2I ("h_nue_daughter_shower_mc_pdg_pfp_hits",
	                                                         "h_nue_daughter_shower_mc_pdg_pfp_hits", 20, 0, 3000, 10, 0, 10);
	TH2D * h_nue_daughter_track_mc_pdg_pfp_hits = new TH2D ("h_nue_daughter_track_mc_pdg_pfp_hits",
	                                                        "h_nue_daughter_track_mc_pdg_pfp_hits", 20, 0, 3000, 10, 0, 10);
	TH2I * h_nue_daughter_origin_mc_pdg = new TH2I ("h_nue_daughter_origin_mc_pdg", "h_nue_daughter_origin_mc_pdg", 3, 0, 3, 10, 0, 10);
	TH2I * h_nue_daughter_origin_mc_pdg_shwr = new TH2I ("h_nue_daughter_origin_mc_pdg_shwr", "h_nue_daughter_origin_mc_pdg_shwr", 3, 0, 3, 10, 0, 10);
	TH2I * h_nue_daughter_origin_mc_pdg_trk  = new TH2I ("h_nue_daughter_origin_mc_pdg_trk", "h_nue_daughter_origin_mc_pdg_trk", 3, 0, 3, 10, 0, 10);


	//here I modify the names of the axis labels
	const char * str_origin[3] = {"kBeamNeutrino", "kCosmicRay", "kUnknown"};
	for (int i=1; i<= 3; i++)
	{
		h_nue_daughter_origin->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_nue_daughter_origin_mc_pdg->GetXaxis()->SetBinLabel(i, str_origin[i-1]);
		h_nue_daughter_origin_mc_pdg_shwr->GetXaxis()->SetBinLabel(i, str_origin[i-1]);
		h_nue_daughter_origin_mc_pdg_trk->GetXaxis()->SetBinLabel(i, str_origin[i-1]);
		h_leading_shower_origin->GetXaxis()->SetBinLabel(i, str_origin[i-1]);
		h_leading_shower_origin_mc_pdg->GetXaxis()->SetBinLabel(i, str_origin[i-1]);
		h_leading_shower_origin_hits->GetXaxis()->SetBinLabel(i, str_origin[i-1]);
		h_leading_shower_origin_tpco_cat->GetXaxis()->SetBinLabel(i, str_origin[i-1]);
	}
	const char * str_pfp_particles[3] = {"Shower", "Track", "Other"};
	for (int i = 1; i <= 3; i++)
	{
		h_nue_daughter_pfp_pdg->GetXaxis()->SetBinLabel(i, str_pfp_particles[i-1]);
		h_nue_daughter_mc_pfp_pdg->GetYaxis()->SetBinLabel(i, str_pfp_particles[i-1]);
	}
	const char * str_mc_particle[12] = {"Electron", "Positron", "Muon", "Mu+", "Photon",
		                            "Pion", "Proton", "Neutron", "Kaon", "No Match"};
	for (int i = 1; i <= 10; i++)
	{
		h_nue_daughter_mc_pdg->GetXaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_nue_daughter_mc_pfp_pdg->GetXaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_nue_daughter_shower_mc_pdg_pfp_hits->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_nue_daughter_track_mc_pdg_pfp_hits->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_nue_daughter_origin_mc_pdg->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_nue_daughter_origin_mc_pdg_shwr->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_nue_daughter_origin_mc_pdg_trk->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg->GetXaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_origin_mc_pdg->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_nue_cc->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_cosmic->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_nue_nc->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_numu->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_other_mixed->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_unmatched->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_numu_zoom->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
		h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->GetYaxis()->SetBinLabel(i, str_mc_particle[i-1]);
	}
	const char * str_tpco_cat[7] = {"Nue CC", "Nue CC Mixed", "Cosmic", "Nue NC",
		                        "Numu", "Other Mixed", "Unmatched"};
	for (int i = 1; i <=7; i++)
	{
		h_leading_shower_origin_tpco_cat->GetYaxis()->SetBinLabel(i, str_tpco_cat[i-1]);
	}

	//--------------------
	//loop for plotting...
	for(int i = 0; i < reco_nue_v_origin.size(); i++)
	{
		auto const this_origin = reco_nue_v_origin.at(i);
		const int this_mc_pdg  = reco_nue_v_mc_pdg.at(i);
		const int this_pfp_pdg = reco_nue_v_pfp_pdg.at(i);

		if(this_origin == "kBeamNeutrino")
		{
			h_nue_daughter_origin->Fill(0);
			if(this_mc_pdg == 11)                         {h_nue_daughter_origin_mc_pdg->Fill(0.0, 0.0);  }
			if(this_mc_pdg == -11)                        {h_nue_daughter_origin_mc_pdg->Fill(0.0, 1.0);  }
			if(this_mc_pdg == 13)                         {h_nue_daughter_origin_mc_pdg->Fill(0.0, 2.0);  }
			if(this_mc_pdg == -13)                        {h_nue_daughter_origin_mc_pdg->Fill(0.0, 3.0);  }
			if(this_mc_pdg == 22)                         {h_nue_daughter_origin_mc_pdg->Fill(0.0, 4.0);  }
			if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_origin_mc_pdg->Fill(0.0, 5.0);  }
			if(this_mc_pdg == 2212)                       {h_nue_daughter_origin_mc_pdg->Fill(0.0, 6.0);  }
			if(this_mc_pdg == 2112)                       {h_nue_daughter_origin_mc_pdg->Fill(0.0, 7.0);  }
			if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
			   this_mc_pdg == 311 || this_mc_pdg == 321 ||
			   this_mc_pdg == -321)                       {h_nue_daughter_origin_mc_pdg->Fill(0.0, 8.0); }
			if(this_mc_pdg == 0)                          {h_nue_daughter_origin_mc_pdg->Fill(0.0, 9.0); }
		}
		if(this_origin == "kCosmicRay")
		{
			h_nue_daughter_origin->Fill(1);
			if(this_mc_pdg == 11)                         {h_nue_daughter_origin_mc_pdg->Fill(1.0, 0.0);  }
			if(this_mc_pdg == -11)                        {h_nue_daughter_origin_mc_pdg->Fill(1.0, 1.0);  }
			if(this_mc_pdg == 13)                         {h_nue_daughter_origin_mc_pdg->Fill(1.0, 2.0);  }
			if(this_mc_pdg == -13)                        {h_nue_daughter_origin_mc_pdg->Fill(1.0, 3.0);  }
			if(this_mc_pdg == 22)                         {h_nue_daughter_origin_mc_pdg->Fill(1.0, 4.0);  }
			if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_origin_mc_pdg->Fill(1.0, 5.0);  }
			if(this_mc_pdg == 2212)                       {h_nue_daughter_origin_mc_pdg->Fill(1.0, 6.0);  }
			if(this_mc_pdg == 2112)                       {h_nue_daughter_origin_mc_pdg->Fill(1.0, 7.0);  }
			if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
			   this_mc_pdg == 311 || this_mc_pdg == 321 ||
			   this_mc_pdg == -321)                       {h_nue_daughter_origin_mc_pdg->Fill(1.0, 8.0); }
			if(this_mc_pdg == 0)                          {h_nue_daughter_origin_mc_pdg->Fill(1.0, 9.0); }
		}
		if(this_origin == "kUnknown")
		{
			h_nue_daughter_origin->Fill(2);
			if(this_mc_pdg == 11)                         {h_nue_daughter_origin_mc_pdg->Fill(2.0, 0.0);  }
			if(this_mc_pdg == -11)                        {h_nue_daughter_origin_mc_pdg->Fill(2.0, 1.0);  }
			if(this_mc_pdg == 13)                         {h_nue_daughter_origin_mc_pdg->Fill(2.0, 2.0);  }
			if(this_mc_pdg == -13)                        {h_nue_daughter_origin_mc_pdg->Fill(2.0, 3.0);  }
			if(this_mc_pdg == 22)                         {h_nue_daughter_origin_mc_pdg->Fill(2.0, 4.0);  }
			if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_origin_mc_pdg->Fill(2.0, 5.0);  }
			if(this_mc_pdg == 2212)                       {h_nue_daughter_origin_mc_pdg->Fill(2.0, 6.0);  }
			if(this_mc_pdg == 2112)                       {h_nue_daughter_origin_mc_pdg->Fill(2.0, 7.0);  }
			if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
			   this_mc_pdg == 311 || this_mc_pdg == 321 ||
			   this_mc_pdg == -321)                       {h_nue_daughter_origin_mc_pdg->Fill(2.0, 8.0); }
			if(this_mc_pdg == 0)                          {h_nue_daughter_origin_mc_pdg->Fill(2.0, 9.0); }
		}
		if(this_pfp_pdg == 11)
		{
			if(this_origin == "kBeamNeutrino")
			{
				if(this_mc_pdg == 11)                         {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 0.0);  }
				if(this_mc_pdg == -11)                        {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 1.0);  }
				if(this_mc_pdg == 13)                         {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 2.0);  }
				if(this_mc_pdg == -13)                        {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 3.0);  }
				if(this_mc_pdg == 22)                         {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 4.0);  }
				if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 5.0);  }
				if(this_mc_pdg == 2212)                       {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 6.0);  }
				if(this_mc_pdg == 2112)                       {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 7.0);  }
				if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
				   this_mc_pdg == 311 || this_mc_pdg == 321 ||
				   this_mc_pdg == -321)                       {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 8.0); }
				if(this_mc_pdg == 0)                          {h_nue_daughter_origin_mc_pdg_shwr->Fill(0.0, 9.0); }
			}
			if(this_origin == "kCosmicRay")
			{
				if(this_mc_pdg == 11)                         {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 0.0);  }
				if(this_mc_pdg == -11)                        {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 1.0);  }
				if(this_mc_pdg == 13)                         {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 2.0);  }
				if(this_mc_pdg == -13)                        {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 3.0);  }
				if(this_mc_pdg == 22)                         {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 4.0);  }
				if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 5.0);  }
				if(this_mc_pdg == 2212)                       {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 6.0);  }
				if(this_mc_pdg == 2112)                       {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 7.0);  }
				if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
				   this_mc_pdg == 311 || this_mc_pdg == 321 ||
				   this_mc_pdg == -321)                       {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 8.0); }
				if(this_mc_pdg == 0)                          {h_nue_daughter_origin_mc_pdg_shwr->Fill(1.0, 9.0); }
			}
			if(this_origin == "kUnknown")
			{
				if(this_mc_pdg == 11)                         {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 0.0);  }
				if(this_mc_pdg == -11)                        {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 1.0);  }
				if(this_mc_pdg == 13)                         {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 2.0);  }
				if(this_mc_pdg == -13)                        {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 3.0);  }
				if(this_mc_pdg == 22)                         {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 4.0);  }
				if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 5.0);  }
				if(this_mc_pdg == 2212)                       {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 6.0);  }
				if(this_mc_pdg == 2112)                       {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 7.0);  }
				if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
				   this_mc_pdg == 311 || this_mc_pdg == 321 ||
				   this_mc_pdg == -321)                       {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 8.0); }
				if(this_mc_pdg == 0)                          {h_nue_daughter_origin_mc_pdg_shwr->Fill(2.0, 9.0); }
			}
		}
		if(this_pfp_pdg == 13)
		{
			if(this_origin == "kBeamNeutrino")
			{
				if(this_mc_pdg == 11)                         {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 0.0);  }
				if(this_mc_pdg == -11)                        {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 1.0);  }
				if(this_mc_pdg == 13)                         {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 2.0);  }
				if(this_mc_pdg == -13)                        {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 3.0);  }
				if(this_mc_pdg == 22)                         {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 4.0);  }
				if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 5.0);  }
				if(this_mc_pdg == 2212)                       {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 6.0);  }
				if(this_mc_pdg == 2112)                       {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 7.0);  }
				if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
				   this_mc_pdg == 311 || this_mc_pdg == 321 ||
				   this_mc_pdg == -321)                       {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 8.0); }
				if(this_mc_pdg == 0)                          {h_nue_daughter_origin_mc_pdg_trk->Fill(0.0, 9.0); }
			}
			if(this_origin == "kCosmicRay")
			{
				if(this_mc_pdg == 11)                         {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 0.0);  }
				if(this_mc_pdg == -11)                        {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 1.0);  }
				if(this_mc_pdg == 13)                         {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 2.0);  }
				if(this_mc_pdg == -13)                        {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 3.0);  }
				if(this_mc_pdg == 22)                         {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 4.0);  }
				if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 5.0);  }
				if(this_mc_pdg == 2212)                       {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 6.0);  }
				if(this_mc_pdg == 2112)                       {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 7.0);  }
				if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
				   this_mc_pdg == 311 || this_mc_pdg == 321 ||
				   this_mc_pdg == -321)                       {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 8.0); }
				if(this_mc_pdg == 0)                          {h_nue_daughter_origin_mc_pdg_trk->Fill(1.0, 9.0); }
			}
			if(this_origin == "kUnknown")
			{
				if(this_mc_pdg == 11)                         {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 0.0);  }
				if(this_mc_pdg == -11)                        {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 1.0);  }
				if(this_mc_pdg == 13)                         {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 2.0);  }
				if(this_mc_pdg == -13)                        {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 3.0);  }
				if(this_mc_pdg == 22)                         {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 4.0);  }
				if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 5.0);  }
				if(this_mc_pdg == 2212)                       {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 6.0);  }
				if(this_mc_pdg == 2112)                       {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 7.0);  }
				if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
				   this_mc_pdg == 311 || this_mc_pdg == 321 ||
				   this_mc_pdg == -321)                       {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 8.0); }
				if(this_mc_pdg == 0)                          {h_nue_daughter_origin_mc_pdg_trk->Fill(2.0, 9.0); }
			}
		}
		//this next line should be 0 as pfp neutrinos are removed already
		if(this_pfp_pdg != 11 && this_pfp_pdg != 13) {h_nue_daughter_pfp_pdg->Fill(2);  }

		if(this_mc_pdg == 11)                         {h_nue_daughter_mc_pdg->Fill(0);  }
		if(this_mc_pdg == -11)                        {h_nue_daughter_mc_pdg->Fill(1);  }
		if(this_mc_pdg == 13)                         {h_nue_daughter_mc_pdg->Fill(2);  }
		if(this_mc_pdg == -13)                        {h_nue_daughter_mc_pdg->Fill(3);  }
		if(this_mc_pdg == 22)                         {h_nue_daughter_mc_pdg->Fill(4);  }
		if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_mc_pdg->Fill(5);  }
		if(this_mc_pdg == 2212)                       {h_nue_daughter_mc_pdg->Fill(6);  }
		if(this_mc_pdg == 2112)                       {h_nue_daughter_mc_pdg->Fill(7);  }
		if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
		   this_mc_pdg == 311 || this_mc_pdg == 321 ||
		   this_mc_pdg == -321)                       {h_nue_daughter_mc_pdg->Fill(8);  }
		if(this_mc_pdg == 0)                          {h_nue_daughter_mc_pdg->Fill(9); }
		//else{h_nue_daughter_mc_pdg->Fill(10); }
		const int this_pfp_hits = reco_nue_v_pfp_hits.at(i);
		h_nue_daughter_pfp_hits->Fill(this_pfp_hits);

		if(this_pfp_pdg == 11 && this_mc_pdg == 11)   {h_nue_daughter_mc_pfp_pdg->Fill(0.0, 0.0); }
		if(this_pfp_pdg == 11 && this_mc_pdg == -11)  {h_nue_daughter_mc_pfp_pdg->Fill(1.0, 0.0); }
		if(this_pfp_pdg == 11 && this_mc_pdg == 13)   {h_nue_daughter_mc_pfp_pdg->Fill(2.0, 0.0); }
		if(this_pfp_pdg == 11 && this_mc_pdg == -13)  {h_nue_daughter_mc_pfp_pdg->Fill(3.0, 0.0); }
		if(this_pfp_pdg == 11 && this_mc_pdg == 22)   {h_nue_daughter_mc_pfp_pdg->Fill(4.0, 0.0); }
		if(this_pfp_pdg == 11 && (this_mc_pdg == 211 || this_mc_pdg == -211)) {h_nue_daughter_mc_pfp_pdg->Fill(5.0, 0.0); }
		if(this_pfp_pdg == 11 && this_mc_pdg == 2212) {h_nue_daughter_mc_pfp_pdg->Fill(6.0, 0.0); }
		if(this_pfp_pdg == 11 && this_mc_pdg == 2112) {h_nue_daughter_mc_pfp_pdg->Fill(7.0, 0.0); }
		if(this_pfp_pdg == 11 && (this_mc_pdg == 130 ||
		                          this_mc_pdg == 310 || this_mc_pdg == 311 ||
		                          this_mc_pdg == 321 || this_mc_pdg == -321)) {h_nue_daughter_mc_pfp_pdg->Fill(8.0, 0.0); }
		if(this_pfp_pdg == 11 && this_mc_pdg == 0)    {h_nue_daughter_mc_pfp_pdg->Fill(9.0, 0.0); }

		if(this_pfp_pdg == 13 && this_mc_pdg == 11)   {h_nue_daughter_mc_pfp_pdg->Fill(0.0, 1.0); }
		if(this_pfp_pdg == 13 && this_mc_pdg == -11)  {h_nue_daughter_mc_pfp_pdg->Fill(1.0, 1.0); }
		if(this_pfp_pdg == 13 && this_mc_pdg == 13)   {h_nue_daughter_mc_pfp_pdg->Fill(2.0, 1.0); }
		if(this_pfp_pdg == 13 && this_mc_pdg == -13)  {h_nue_daughter_mc_pfp_pdg->Fill(3.0, 1.0); }
		if(this_pfp_pdg == 13 && this_mc_pdg == 22)   {h_nue_daughter_mc_pfp_pdg->Fill(4.0, 1.0); }
		if(this_pfp_pdg == 13 && (this_mc_pdg == 211 || this_mc_pdg == -211)) {h_nue_daughter_mc_pfp_pdg->Fill(5.0, 1.0); }
		if(this_pfp_pdg == 13 && this_mc_pdg == 2212) {h_nue_daughter_mc_pfp_pdg->Fill(6.0, 1.0); }
		if(this_pfp_pdg == 13 && this_mc_pdg == 2112) {h_nue_daughter_mc_pfp_pdg->Fill(7.0, 1.0); }
		if(this_pfp_pdg == 11 && (this_mc_pdg == 130 ||
		                          this_mc_pdg == 310 || this_mc_pdg == 311 ||
		                          this_mc_pdg == 321 || this_mc_pdg == -321)) {h_nue_daughter_mc_pfp_pdg->Fill(8.0, 1.0); }
		if(this_pfp_pdg == 13 && this_mc_pdg == 0)    {h_nue_daughter_mc_pfp_pdg->Fill(9.0, 1.0); }

		if(this_pfp_pdg == 11)
		{
			if(this_mc_pdg == 11)                         {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 0.0); }
			if(this_mc_pdg == -11)                        {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 1.0); }
			if(this_mc_pdg == 13)                         {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 2.0); }
			if(this_mc_pdg == -13)                        {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 3.0); }
			if(this_mc_pdg == 22)                         {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 4.0); }
			if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 5.0); }
			if(this_mc_pdg == 2212)                       {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 6.0); }
			if(this_mc_pdg == 2112)                       {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 7.0); }
			if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
			   this_mc_pdg == 311 || this_mc_pdg == 321 ||
			   this_mc_pdg == -321)                       {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 8.0);  }
			if(this_mc_pdg == 0)                          {h_nue_daughter_shower_mc_pdg_pfp_hits->Fill(this_pfp_hits, 9.0); }
		}
		if(this_pfp_pdg == 13)
		{
			if(this_mc_pdg == 11)                         {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 0.0); }
			if(this_mc_pdg == -11)                        {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 1.0); }
			if(this_mc_pdg == 13)                         {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 2.0); }
			if(this_mc_pdg == -13)                        {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 3.0); }
			if(this_mc_pdg == 22)                         {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 4.0); }
			if(this_mc_pdg == 211 || this_mc_pdg == -211) {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 5.0); }
			if(this_mc_pdg == 2212)                       {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 6.0); }
			if(this_mc_pdg == 2112)                       {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 7.0); }
			if(this_mc_pdg == 130 || this_mc_pdg == 310 ||
			   this_mc_pdg == 311 || this_mc_pdg == 321 ||
			   this_mc_pdg == -321)                       {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 8.0);  }
			if(this_mc_pdg == 0)                          {h_nue_daughter_track_mc_pdg_pfp_hits->Fill(this_pfp_hits, 9.0); }
		}

	}
	TCanvas * tpco_c1 = new TCanvas();
	tpco_c1->cd();
	h_nue_daughter_origin->GetYaxis()->SetTitleOffset(1.3);
	h_nue_daughter_origin->GetXaxis()->SetTitle("Reco Nue Daughter Origin");
	h_nue_daughter_origin->GetYaxis()->SetTitle("PFPs");
	h_nue_daughter_origin->Draw();
	tpco_c1->Print("nue_daughter_origin.pdf");
	TCanvas * tpco_c2 = new TCanvas();
	tpco_c2->cd();
	h_nue_daughter_pfp_pdg->GetYaxis()->SetTitleOffset(1.3);
	h_nue_daughter_pfp_pdg->GetXaxis()->SetTitle("Reco Nue Daughter PF Particle");
	h_nue_daughter_pfp_pdg->GetYaxis()->SetTitle("PFPs");
	h_nue_daughter_pfp_pdg->Draw();
	tpco_c2->Print("nue_daughter_pfp_pdg.pdf");
	TCanvas * tpco_c3 = new TCanvas();
	tpco_c3->cd();
	h_nue_daughter_mc_pdg->GetYaxis()->SetTitleOffset(1.3);
	h_nue_daughter_mc_pdg->GetXaxis()->SetTitleOffset(1.2);
	h_nue_daughter_mc_pdg->GetXaxis()->SetTitle("Reco Nue Daughter MC Particle");
	h_nue_daughter_mc_pdg->GetYaxis()->SetTitle("PFPs");
	h_nue_daughter_mc_pdg->Draw();
	tpco_c3->Print("nue_daughter_mc_pdg.pdf");
	TCanvas * tpco_c4 = new TCanvas();
	tpco_c4->cd();
	h_nue_daughter_pfp_hits->GetYaxis()->SetTitleOffset(1.3);
	h_nue_daughter_pfp_hits->GetXaxis()->SetTitle("Reco Nue Daughter PFP Hits");
	h_nue_daughter_pfp_hits->GetYaxis()->SetTitle("PFPs");
	h_nue_daughter_pfp_hits->Draw();
	tpco_c4->Print("nue_daughter_pfp_hits.pdf");
	TCanvas * tpco_c5 = new TCanvas();
	tpco_c5->cd();
	tpco_c5->SetLogz();
	h_nue_daughter_mc_pfp_pdg->GetYaxis()->SetTitleOffset(1.3);
	h_nue_daughter_mc_pfp_pdg->GetXaxis()->SetTitle("Reco Nue Daughter MC Particle");
	h_nue_daughter_mc_pfp_pdg->GetYaxis()->SetTitle("Reco Nue Daughter PF Particle");
	h_nue_daughter_mc_pfp_pdg->SetStats(kFALSE);
	h_nue_daughter_mc_pfp_pdg->Draw("colz");
	tpco_c5->Print("nue_daughter_mc_pfp_pdg.pdf");
	TCanvas * tpco_c6 = new TCanvas();
	tpco_c6->cd();
	tpco_c6->SetLogz();
	h_nue_daughter_shower_mc_pdg_pfp_hits->GetYaxis()->SetLabelOffset(0.002);
	h_nue_daughter_shower_mc_pdg_pfp_hits->GetYaxis()->SetTitleOffset(1.35);
	h_nue_daughter_shower_mc_pdg_pfp_hits->GetXaxis()->SetTitle("Reco Nue Daughter - Shower PFP Hits");
	h_nue_daughter_shower_mc_pdg_pfp_hits->GetYaxis()->SetTitle("Reco Nue Daughter - Shower MC Particle");
	h_nue_daughter_shower_mc_pdg_pfp_hits->SetStats(kFALSE);
	h_nue_daughter_shower_mc_pdg_pfp_hits->Draw("colz");
	tpco_c6->Print("nue_daughter_shower_mc_pdg_pfp_hits.pdf");
	TCanvas * tpco_c7 = new TCanvas();
	tpco_c7->cd();
	tpco_c7->SetLogz();
	h_nue_daughter_track_mc_pdg_pfp_hits->GetYaxis()->SetLabelOffset(0.002);
	h_nue_daughter_track_mc_pdg_pfp_hits->GetYaxis()->SetTitleOffset(1.35);
	h_nue_daughter_track_mc_pdg_pfp_hits->GetXaxis()->SetTitle("Reco Nue Daughter - Track PFP Hits");
	h_nue_daughter_track_mc_pdg_pfp_hits->GetYaxis()->SetTitle("Reco Nue Daughter - Track MC Particle");
	h_nue_daughter_track_mc_pdg_pfp_hits->SetStats(kFALSE);
	h_nue_daughter_track_mc_pdg_pfp_hits->Draw("colz");
	tpco_c7->Print("nue_daughter_track_mc_pdg_pfp_hits.pdf");
	TCanvas * tpco_c8 = new TCanvas();
	tpco_c8->cd();
	tpco_c8->SetLogz();
	h_nue_daughter_origin_mc_pdg->GetYaxis()->SetLabelOffset(0.002);
	h_nue_daughter_origin_mc_pdg->GetYaxis()->SetTitleOffset(1.35);
	h_nue_daughter_origin_mc_pdg->GetYaxis()->SetTitle("Reco Nue Daughter MC Particle");
	h_nue_daughter_origin_mc_pdg->GetXaxis()->SetTitle("Reco Nue Daughter Origin");
	h_nue_daughter_origin_mc_pdg->SetStats(kFALSE);
	h_nue_daughter_origin_mc_pdg->Draw("colz");
	tpco_c8->Print("nue_daughter_origin_mc_pdg.pdf");
	TCanvas * tpco_c9 = new TCanvas();
	tpco_c9->cd();
	tpco_c9->SetLogz();
	h_nue_daughter_origin_mc_pdg_shwr->GetYaxis()->SetLabelOffset(0.002);
	h_nue_daughter_origin_mc_pdg_shwr->GetYaxis()->SetTitleOffset(1.35);
	h_nue_daughter_origin_mc_pdg_shwr->GetYaxis()->SetTitle("Reco Nue Daughter - Shower MC Particle");
	h_nue_daughter_origin_mc_pdg_shwr->GetXaxis()->SetTitle("Reco Nue Daughter - Shower Origin");
	h_nue_daughter_origin_mc_pdg_shwr->SetStats(kFALSE);
	h_nue_daughter_origin_mc_pdg_shwr->Draw("colz");
	tpco_c9->Print("nue_daughter_origin_mc_pdg_shwr.pdf");
	TCanvas * tpco_c10 = new TCanvas();
	tpco_c10->cd();
	tpco_c10->SetLogz();
	h_nue_daughter_origin_mc_pdg_trk->GetYaxis()->SetLabelOffset(0.002);
	h_nue_daughter_origin_mc_pdg_trk->GetYaxis()->SetTitleOffset(1.35);
	h_nue_daughter_origin_mc_pdg_trk->GetYaxis()->SetTitle("Reco Nue Daughter - Track MC Particle");
	h_nue_daughter_origin_mc_pdg_trk->GetXaxis()->SetTitle("Reco Nue Daughter - Track Origin");
	h_nue_daughter_origin_mc_pdg_trk->SetStats(kFALSE);
	h_nue_daughter_origin_mc_pdg_trk->Draw("colz");
	tpco_c10->Print("nue_daughter_origin_mc_pdg_trk.pdf");


	TH1D * h_opt_time = new TH1D("h_opt_time", "h_opt_time", 50, 0, 20);
	h_opt_time->GetXaxis()->SetTitle("Time [us]");
	h_opt_time->GetYaxis()->SetTitle("Flashes");
	TCanvas * opt_c1 = new TCanvas();
	TH1D * h_opt_pe = new TH1D("h_opt_pe", "h_opt_pe", 50, 0, 500);
	h_opt_pe->GetYaxis()->SetTitleOffset(1.3);
	h_opt_pe->GetXaxis()->SetTitle("Photoelectrons");
	h_opt_pe->GetYaxis()->SetTitle("Flashes");
	TCanvas * opt_c2 = new TCanvas();
	TH2D * h_opt_time_pe = new TH2D("h_opt_time_pe", "h_opt_time_pe", 20, 0, 20, 20, 0, 200);
	h_opt_time_pe->GetYaxis()->SetTitleOffset(1.3);
	h_opt_time_pe->GetXaxis()->SetTitle("Time [us]");
	h_opt_time_pe->GetYaxis()->SetTitle("Photoelectrons");
	h_opt_time_pe->SetStats(kFALSE);
	TCanvas * opt_c3 = new TCanvas();


	TCanvas * ls_c1 = new TCanvas();
	ls_c1->cd();
	h_leading_shower_mc_pdg->GetXaxis()->SetTitleOffset(1.2);
	h_leading_shower_mc_pdg->GetXaxis()->SetTitle("Nue Candidate Leading Shower True PDG");
	h_leading_shower_mc_pdg->Draw();
	ls_c1->Print("leading_shower_mc_pdg.pdf");
	TCanvas * ls_c2 = new TCanvas();
	ls_c2->cd();
	h_leading_shower_origin->GetXaxis()->SetTitle("Nue Candidate Leading Shower Origin");
	h_leading_shower_origin->Draw();
	ls_c2->Print("leading_shower_origin.pdf");
	TCanvas * ls_c3 = new TCanvas();
	ls_c3->cd();
	h_leading_shower_origin_mc_pdg->GetXaxis()->SetTitle("Nue Candidate Leading Shower Origin");
	h_leading_shower_origin_mc_pdg->GetYaxis()->SetTitleOffset(1.3);
	h_leading_shower_origin_mc_pdg->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_origin_mc_pdg->GetYaxis()->SetTitle("Nue Candidate Leading Shower True PDG");
	h_leading_shower_origin_mc_pdg->SetStats(kFALSE);
	h_leading_shower_origin_mc_pdg->Draw("colz");
	ls_c3->Print("leading_shower_origin_mc_pdg.pdf");
	TCanvas * ls_c4 = new TCanvas();
	ls_c4->cd();
	h_leading_shower_origin_hits->GetXaxis()->SetTitle("Nue Candidate Leading Shower Orgin");
	h_leading_shower_origin_hits->GetYaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_origin_hits->SetStats(kFALSE);
	h_leading_shower_origin_hits->Draw("colz");
	ls_c4->Print("leading_shower_origin_hits.pdf");
	TCanvas * ls_c5 = new TCanvas();
	ls_c5->cd();
	h_leading_shower_origin_tpco_cat->GetXaxis()->SetTitle("Nue Candidate Leading Shower Origin");
	h_leading_shower_origin_tpco_cat->GetYaxis()->SetTitleOffset(1.6);
	h_leading_shower_origin_tpco_cat->GetYaxis()->SetLabelOffset(0.002);
	h_leading_shower_origin_tpco_cat->GetYaxis()->SetTitle("TPC Object classification");
	h_leading_shower_origin_tpco_cat->SetStats(kFALSE);
	h_leading_shower_origin_tpco_cat->Draw("colz");
	ls_c5->Print("leading_shower_origin_tpco_cat.pdf");


	TCanvas * ls_c6 = new TCanvas();
	ls_c6->cd();
	h_leading_shower_mc_pdg_pfp_hits_nue_cc->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_nue_cc->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_nue_cc->Draw("colz");
	ls_c6->Print("leading_shower_mc_pdg_pfp_hits_nue_cc.pdf");
	TCanvas * ls_c7 = new TCanvas();
	ls_c7->cd();
	h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed->Draw("colz");
	ls_c7->Print("leading_shower_mc_pdg_pfp_hits_nue_cc_mixed.pdf");
	TCanvas * ls_c8 = new TCanvas();
	ls_c8->cd();
	h_leading_shower_mc_pdg_pfp_hits_cosmic->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_cosmic->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_cosmic->Draw("colz");
	ls_c8->Print("leading_shower_mc_pdg_pfp_hits_cosmic.pdf");
	TCanvas * ls_c9 = new TCanvas();
	ls_c9->cd();
	h_leading_shower_mc_pdg_pfp_hits_nue_nc->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_nue_nc->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_nue_nc->Draw("colz");
	ls_c9->Print("leading_shower_mc_pdg_pfp_hits_nue_nc.pdf");
	TCanvas * ls_c10 = new TCanvas();
	ls_c10->cd();
	h_leading_shower_mc_pdg_pfp_hits_numu->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_numu->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_numu->Draw("colz");
	ls_c10->Print("leading_shower_mc_pdg_pfp_hits_numu.pdf");
	TCanvas * ls_c11 = new TCanvas();
	ls_c11->cd();
	h_leading_shower_mc_pdg_pfp_hits_other_mixed->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_other_mixed->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_other_mixed->Draw("colz");
	ls_c11->Print("leading_shower_mc_pdg_pfp_hits_other_mixed.pdf");
	TCanvas * ls_c12 = new TCanvas();
	ls_c12->cd();
	h_leading_shower_mc_pdg_pfp_hits_unmatched->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_unmatched->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_unmatched->Draw("colz");
	ls_c12->Print("leading_shower_mc_pdg_pfp_hits_unmatched.pdf");

	TCanvas * ls_c6_zoom = new TCanvas();
	ls_c6_zoom->cd();
	h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_nue_cc_zoom->Draw("colz");
	ls_c6_zoom->Print("leading_shower_mc_pdg_pfp_hits_nue_cc_zoom.pdf");
	TCanvas * ls_c7_zoom = new TCanvas();
	ls_c7_zoom->cd();
	h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom->Draw("colz");
	ls_c7_zoom->Print("leading_shower_mc_pdg_pfp_hits_nue_cc_mixed_zoom.pdf");
	TCanvas * ls_c8_zoom = new TCanvas();
	ls_c8_zoom->cd();
	h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_cosmic_zoom->Draw("colz");
	ls_c8_zoom->Print("leading_shower_mc_pdg_pfp_hits_cosmic_zoom.pdf");
	TCanvas * ls_c9_zoom = new TCanvas();
	ls_c9_zoom->cd();
	h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_nue_nc_zoom->Draw("colz");
	ls_c9_zoom->Print("leading_shower_mc_pdg_pfp_hits_nue_nc_zoom.pdf");
	TCanvas * ls_c10_zoom = new TCanvas();
	ls_c10_zoom->cd();
	h_leading_shower_mc_pdg_pfp_hits_numu_zoom->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_numu_zoom->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_numu_zoom->Draw("colz");
	ls_c10_zoom->Print("leading_shower_mc_pdg_pfp_hits_numu_zoom.pdf");
	TCanvas * ls_c11_zoom = new TCanvas();
	ls_c11_zoom->cd();
	h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_other_mixed_zoom->Draw("colz");
	ls_c11_zoom->Print("leading_shower_mc_pdg_pfp_hits_other_mixed_zoom.pdf");
	TCanvas * ls_c12_zoom = new TCanvas();
	ls_c12_zoom->cd();
	h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->GetXaxis()->SetTitle("Reconstructed Hits");
	h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->GetYaxis()->SetTitle("Leading Shower True PDG");
	h_leading_shower_mc_pdg_pfp_hits_unmatched_zoom->Draw("colz");
	ls_c12_zoom->Print("leading_shower_mc_pdg_pfp_hits_unmatched_zoom.pdf");


	double largest_flash = 0;
	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;
	std::vector < std::vector < double > > largest_flash_v_v;
	const double flash_time_min = 5;
	const double flash_time_max = 16;

	const int opt_entries = optree->GetEntries();
	for(int i = 0; i < opt_entries; i++)
	{
		optree->GetEntry(i);
		// std::cout << "[Optical Info] \t Event:   " << fEvent << std::endl;
		// std::cout << "[Optical Info] \t Run:     " << fRun << std::endl;
		// std::cout << "[Optical Info] \t PE:      " << fOpFlashPE << std::endl;
		// std::cout << "[Optical Info] \t Time:    " << fOpFlashTime << std::endl;
		// std::cout << "[Optical Info] \t WidthY:  " << fOpFlashWidthY << std::endl;
		// std::cout << "[Optical Info] \t WidthZ:  " << fOpFlashWidthZ << std::endl;
		// std::cout << "[Optical Info] \t CenterY: " << fOpFlashCenterY << std::endl;
		// std::cout << "[Optical Info] \t CenterZ: " << fOpFlashCenterZ << std::endl;
		if(fOpFlashTime <= flash_time_max && fOpFlashTime >= flash_time_min) {h_opt_pe->Fill(fOpFlashPE); }
		h_opt_time->Fill(fOpFlashTime);
		h_opt_time_pe->Fill(fOpFlashTime, fOpFlashPE);
		current_run = fRun;
		current_event = fEvent;
		if(current_event != last_event) {largest_flash = 0; }
		double this_flash = fOpFlashPE;
		std::vector < double > largest_flash_v;                //contains the y,z for largest flash


		if(this_flash <= largest_flash)
		{
			last_event = current_event;
			last_run = current_run;
			continue;
		}
		if(current_event == last_event && current_run == last_run)
		{
			if(this_flash > largest_flash)
			{
				largest_flash = this_flash;
				largest_flash_v_v.pop_back();
			}
		}
		last_event = current_event;
		last_run = current_run;
		largest_flash = this_flash;
		largest_flash_v.push_back(fOpFlashCenterY);
		largest_flash_v.push_back(fOpFlashCenterZ);
		largest_flash_v.push_back(current_event);
		largest_flash_v.push_back(fOpFlashWidthZ);
		largest_flash_v_v.push_back(largest_flash_v);
		largest_flash_v.clear();
	}


	TH1D * h_flash_dist_beam_electron   = new TH1D("h_flash_dist_beam_electron", "h_flash_dist_beam_electron", 40, 0, 200);
	TH1D * h_flash_dist_beam_photon     = new TH1D("h_flash_dist_beam_photon", "h_flash_dist_beam_photon", 40, 0, 200);
	TH1D * h_flash_dist_cosmic_electron = new TH1D("h_flash_dist_cosmic_electron", "h_flash_dist_cosmic_electron", 40, 0, 200);
	TH1D * h_flash_dist_cosmic_photon   = new TH1D("h_flash_dist_cosmic_photon", "h_flash_dist_cosmic_photon", 40, 0, 200);
	TH1D * h_flash_dist_unknown         = new TH1D("h_flash_dist_unknown", "h_flash_dist_unknown", 40, 0, 200);

	TH1D * h_flash_dist_width_beam_electron   = new TH1D("h_flash_dist_width_beam_electron",   "h_flash_dist_width_beam_electron", 40, -160, 200);
	TH1D * h_flash_dist_width_beam_photon     = new TH1D("h_flash_dist_width_beam_photon",     "h_flash_dist_width_beam_photon", 40, -160, 200);
	TH1D * h_flash_dist_width_cosmic_electron = new TH1D("h_flash_dist_width_cosmic_electron", "h_flash_dist_width_cosmic_electron", 40, -160, 200);
	TH1D * h_flash_dist_width_cosmic_photon   = new TH1D("h_flash_dist_width_cosmic_photon",   "h_flash_dist_width_cosmic_photon", 40, -160, 200);
	TH1D * h_flash_dist_width_unknown         = new TH1D("h_flash_dist_width_unknown",         "h_flash_dist_width_unknown", 40, -160, 200);

	std::cout << "Comparing the Largest Flash to the Shower Vertex" << std::endl;
	for(auto const this_flash_v : largest_flash_v_v)
	{
		const int this_event = this_flash_v.at(2);
		const double flash_vtx_y = this_flash_v.at(0);
		const double flash_vtx_z = this_flash_v.at(1);
		const double flash_width_z = this_flash_v.at(3);
		for(auto const this_shower_v : reco_nue_shower_vtx_v_v)
		{
			if(this_event != this_shower_v.at(3)) {continue; }
			const double shwr_vtx_y = this_shower_v.at(1);
			const double shwr_vtx_z = this_shower_v.at(2);
			const double distance = sqrt(pow((shwr_vtx_y - flash_vtx_y), 2) + pow((shwr_vtx_z - flash_vtx_z), 2) );
			const double distance_width = distance - flash_width_z;
			if(this_shower_v.at(5) == 0)
			{
				if(this_shower_v.at(4) == 11) {h_flash_dist_beam_electron->Fill(distance); h_flash_dist_width_beam_electron->Fill(distance_width); }
				if(this_shower_v.at(4) == 22) {h_flash_dist_beam_photon->Fill(distance); h_flash_dist_width_beam_photon->Fill(distance_width); }
			}
			if(this_shower_v.at(5) == 1)
			{
				if(this_shower_v.at(4) == 11) {h_flash_dist_cosmic_electron->Fill(distance); h_flash_dist_width_cosmic_electron->Fill(distance_width); }
				if(this_shower_v.at(4) == 22) {h_flash_dist_cosmic_photon->Fill(distance); h_flash_dist_width_cosmic_photon->Fill(distance_width); }
			}
			if(this_shower_v.at(5) == 2) {h_flash_dist_unknown->Fill(distance); h_flash_dist_width_unknown->Fill(distance_width); }
		}
	}
	TCanvas * flash_dist_c1 = new TCanvas();
	flash_dist_c1->cd();
	h_flash_dist_beam_electron->GetXaxis()->SetTitle("Largest Flash to Shower Vertex [cm]");
	h_flash_dist_beam_electron->Draw();
	flash_dist_c1->Print("flash_dist_beam_electron.pdf");
	TCanvas * flash_dist_c2 = new TCanvas();
	flash_dist_c2->cd();
	h_flash_dist_beam_photon->GetXaxis()->SetTitle("Largest Flash to Shower Vertex [cm]");
	h_flash_dist_beam_photon->Draw();
	flash_dist_c2->Print("flash_dist_beam_photon.pdf");
	TCanvas * flash_dist_c3 = new TCanvas();
	flash_dist_c3->cd();
	h_flash_dist_cosmic_electron->GetXaxis()->SetTitle("Largest Flash to Shower Vertex [cm]");
	h_flash_dist_cosmic_electron->Draw();
	flash_dist_c3->Print("flash_dist_cosmic_electron.pdf");
	TCanvas * flash_dist_c4 = new TCanvas();
	flash_dist_c4->cd();
	h_flash_dist_cosmic_photon->GetXaxis()->SetTitle("Largest Flash to Shower Vertex [cm]");
	h_flash_dist_cosmic_photon->Draw();
	flash_dist_c4->Print("flash_dist_cosmic_photon.pdf");
	TCanvas * flash_dist_c5 = new TCanvas();
	flash_dist_c5->cd();
	h_flash_dist_unknown->GetXaxis()->SetTitle("Largest Flash to Shower Vertex [cm]");
	h_flash_dist_unknown->Draw();
	flash_dist_c5->Print("flash_dist_unknown.pdf");

	TCanvas * flash_dist_width_c1 = new TCanvas();
	flash_dist_width_c1->cd();
	h_flash_dist_width_beam_electron->GetXaxis()->SetTitle("Largest Flash to Shower Vertex - Flash Width Z [cm]");
	h_flash_dist_width_beam_electron->Draw();
	flash_dist_width_c1->Print("flash_dist_width_beam_electron.pdf");
	TCanvas * flash_dist_width_c2 = new TCanvas();
	flash_dist_width_c2->cd();
	h_flash_dist_width_beam_photon->GetXaxis()->SetTitle("Largest Flash to Shower Vertex - Flash Width Z [cm]");
	h_flash_dist_width_beam_photon->Draw();
	flash_dist_width_c2->Print("flash_dist_width_beam_photon.pdf");
	TCanvas * flash_dist_width_c3 = new TCanvas();
	flash_dist_width_c3->cd();
	h_flash_dist_width_cosmic_electron->GetXaxis()->SetTitle("Largest Flash to Shower Vertex - Flash Width Z [cm]");
	h_flash_dist_width_cosmic_electron->Draw();
	flash_dist_width_c3->Print("flash_dist_width_cosmic_electron.pdf");
	TCanvas * flash_dist_width_c4 = new TCanvas();
	flash_dist_width_c4->cd();
	h_flash_dist_width_cosmic_photon->GetXaxis()->SetTitle("Largest Flash to Shower Vertex - Flash Width Z [cm]");
	h_flash_dist_width_cosmic_photon->Draw();
	flash_dist_width_c4->Print("flash_dist_width_cosmic_photon.pdf");
	TCanvas * flash_dist_width_c5 = new TCanvas();
	flash_dist_width_c5->cd();
	h_flash_dist_width_unknown->GetXaxis()->SetTitle("Largest Flash to Shower Vertex - Flash Width Z [cm]");
	h_flash_dist_width_unknown->Draw();
	flash_dist_width_c5->Print("flash_dist_width_unknown.pdf");

	opt_c1->cd();
	h_opt_time->Draw();
	opt_c1->Print("opt_time.pdf");
	opt_c2->cd();
	h_opt_pe->Draw();
	opt_c2->Print("opt_pe.pdf");
	opt_c3->cd();
	h_opt_time_pe->Draw("colz");
	opt_c3->Print("opt_time_pe.pdf");

	if(f->IsOpen()) {f->Close(); }

	return 0;
}//end out_inspect

#ifndef __ROOTCLING__


int main(int argc, char *argv[]){

	argc = 2;
	const char * file1 = argv[1];
	return out_inspect(file1);
}

#endif
