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

#include <iostream>
#include <vector>

#include "../xsecAna/LinkDef.h"

int out_inspect()
{

	const char * _file1 = "../nue_xsec_extraction.root";
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

	TFile * histo_file = new TFile("histo_file.root", "RECREATE");
	const int total_entries = mytree->GetEntries();

	TH1I * h_tpco_origin_unknown  = new TH1I("h_tpco_origin_unknown", "h_tpco_origin_unknown", 10, 0, 10);
	TH1I * h_tpco_origin_cosmic   = new TH1I("h_tpco_origin_cosmic", "h_tpco_origin_cosmic", 10, 0, 10);
	TH1I * h_tpco_origin_neutrino = new TH1I("h_tpco_origin_neutrino", "h_tpco_origin_neutrino", 10, 0, 10);
	TH2I * h_tpco_origin_neutrino_cosmic = new TH2I("h_tpco_origin_neutrino_cosmic", "h_tpco_origin_neutrino_cosmic", 10, 0, 10, 10, 0, 10);
	TH2I * h_nue_daughter_track_shower = new TH2I("h_nue_daughter_track_shower",
	                                              "h_nue_daughter_track_shower", 5, 0, 5, 5, 0, 5);

	// h_tpco_origin_unknown->SetStats(kFALSE);
	// h_tpco_origin_cosmic->SetStats(kFALSE);
	// h_tpco_origin_neutrino->SetStats(kFALSE);
	h_tpco_origin_neutrino_cosmic->SetStats(kFALSE);


	//THStack *hs = new THStack(h_title,"");

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
			if(tpco_pfp_pdg_code != 12) {continue; }

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

				//this should only let daughters of pfp nue's through
				//we also don't want the neutrino pfps!
				if(pfp_parent_pdg == 12 && (pfp_pdg != 12 && pfp_pdg != 14))
				{
					reco_nue_v_origin.push_back(origin);
					if(origin == "kUnknown") {origin_unknown++; }
					if(origin == "kBeamNeutrino") {origin_beam++; }
					if(origin == "kCosmicRay") {origin_cosmic++; }
					reco_nue_v_pfp_pdg.push_back(pfp_pdg);
					reco_nue_v_mc_pdg.push_back(mc_pdg);
					reco_nue_v_pfp_hits.push_back(n_pfp_hits);
					if(pfp_pdg == 11) {n_showers++; }
					if(pfp_pdg == 13) {n_tracks++; }
				}
			}//end looping particles
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

	//here I modify the names of the axis labels
	const char * str_origin[3] = {"kBeamNeutrino", "kCosmicRay", "kUnknown"};
	for (int i=1; i<= 3; i++)
	{
		h_nue_daughter_origin->GetXaxis()->SetBinLabel(i,str_origin[i-1]);
		h_nue_daughter_origin_mc_pdg->GetXaxis()->SetBinLabel(i, str_origin[i-1]);
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
	}

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
		if(this_pfp_pdg == 11)                       {h_nue_daughter_pfp_pdg->Fill(0);  }
		if(this_pfp_pdg == 13)                       {h_nue_daughter_pfp_pdg->Fill(1);  }
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
	h_nue_daughter_mc_pdg->GetYaxis()->SetTitleOffset(1.3);
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


	TH1D * h_opt_time = new TH1D("h_opt_time", "h_opt_time", 50, 0, 20);
	h_opt_time->GetXaxis()->SetTitle("Time [us]");
	h_opt_time->GetYaxis()->SetTitle("Flashes");
	TCanvas * opt_c1 = new TCanvas();
	TH1D * h_opt_pe = new TH1D("h_opt_pe", "h_opt_pe", 50, 0, 15000);
	h_opt_pe->GetYaxis()->SetTitleOffset(1.3);
	h_opt_pe->GetXaxis()->SetTitle("Photoelectrons");
	h_opt_pe->GetYaxis()->SetTitle("Flashes");
	TCanvas * opt_c2 = new TCanvas();
	TH2D * h_opt_time_pe = new TH2D("h_opt_time_pe", "h_opt_time_pe", 20, 0, 20, 20, 0, 15000);
	h_opt_time_pe->GetYaxis()->SetTitleOffset(1.3);
	h_opt_time_pe->GetXaxis()->SetTitle("Time [us]");
	h_opt_time_pe->GetYaxis()->SetTitle("Photoelectrons");
	TCanvas * opt_c3 = new TCanvas();


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
		h_opt_pe->Fill(fOpFlashPE);
		h_opt_time->Fill(fOpFlashTime);
		h_opt_time_pe->Fill(fOpFlashTime, fOpFlashPE);
	}
	opt_c1->cd();
	h_opt_time->Draw();
	opt_c1->Print("opt_time.pdf");
	opt_c2->cd();
	h_opt_pe->Draw();
	opt_c2->Print("opt_pe.pdf");
	opt_c3->cd();
	h_opt_time_pe->Draw("colz");
	opt_c3->Print("opt_time_pe.pdf");


	return 0;
}//end out_inspect

#ifndef __ROOTCLING__


int main(){

	return out_inspect();
}

#endif
