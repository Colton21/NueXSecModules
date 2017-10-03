#include "../xsecAna/TpcObjectContainer.h"
#include "../xsecAna/ParticleContainer.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TInterpreter.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

#include "../xsecAna/LinkDef.h"

int out_inspect()
{

	const char * _file1 = "../nue_xsec_extraction.root";
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

	const int total_entries = mytree->GetEntries();

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
			const int n_mc_hits = tpc_obj.NumMCHits();
			const int n_pfp_hits = tpc_obj.NumPFPHits();
			const int is_cc = tpc_obj.IsCC();
			const int mode = tpc_obj.Mode();
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

			std::cout << " \t Number of PFParticles: " << n_pfp << "\t Number of PFP Neutrinos: " << n_pfp_nu << std::endl;
			std::cout << " \t Run: " << run_number << " , Event: " << event_number << std::endl;
			std::cout << " \t Origin : " << tpc_obj_origin << std::endl;
			std::cout << " \t Reco Hits: " << n_pfp_hits << "  MC Hits: " << "NOT SET" << std::endl;//n_mc_hits << std::endl;
			std::cout << " \t Interaction Mode: " << mode << " IsCC: " << is_cc << std::endl;
			std::cout << " \t Vertex - Reco : " << tpc_obj_pfp_vtx.at(0) << ", " << tpc_obj_pfp_vtx.at(1) << ", " << tpc_obj_pfp_vtx.at(2) << std::endl;
			std::cout << " \t Vertex - MC   : " << tpc_obj_mc_vtx.at(0)  << ", " << tpc_obj_mc_vtx.at(1)  << ", " << tpc_obj_mc_vtx.at(2)  << std::endl;

			for(int j = 0; j < n_pfp; j++)
			{
				auto const part = tpc_obj.GetParticle(j);
				const int pfp_pdg = part.PFParticlePdgCode();
				const int mc_pdg = part.MCPdgCode();
				//const int mc_parent_pdg = part.MCParentPdg();//not set in module
				const int pfp_parent_pdg = part.PFParticleParentPdgCode();
				const std::string origin = part.Origin();
				const int pfp_mode = part.Mode();
				const int pfp_isCC = part.IsCC();
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
				std::cout << " \t \t Particle PDG Codes: - Reco " << pfp_pdg << " - True " << mc_pdg << " - Reco Parent " << pfp_parent_pdg << std::endl;
				std::cout << " \t \t Origin: " << origin << std::endl;
				std::cout << " \t \t Mode: " << pfp_mode        << "  Is CC: " << pfp_isCC << std::endl;
				std::cout << " \t \t Vertex - Reco : " << pfp_vtx.at(0) << ", " << pfp_vtx.at(1) << ", " << pfp_vtx.at(2) << std::endl;
				std::cout << " \t \t Vertex - MC   : " << mc_vtx.at(0)  << ", " << mc_vtx.at(1)  << ", " << mc_vtx.at(2)  << std::endl;
				std::cout << " \t \t Direction - Reco : " << pfp_dir.at(0) << ", " << pfp_dir.at(1) << ", " << pfp_dir.at(2) << std::endl;
				std::cout << " \t \t Direction - MC   : " << mc_dir.at(0)  << ", " << mc_dir.at(1)  << ", " << mc_dir.at(2)  << std::endl;
				std::cout << " \t \t Length - MC : " << mc_length << " - Reco : " << pfp_length << std::endl;
				std::cout << " \t \t Momentum - MC : " << mc_momentum << " - Reco : " << pfp_momentum << std::endl;
				std::cout << " \t \t Hits - Reco : " << n_pfp_hits << " - MC: " << "NOT SET" << std::endl;//n_mc_hits << std::endl;
				std::cout << " \t \t Open Angle: " << pfp_open_angle << std::endl;
				std::cout << " \t \t ----------------------------------------------------" << std::endl;
			}//end looping particles
			std::cout << "\n " << std::endl;

		}//end looping tpc objects
	}//end looping events

	const int opt_entries = optree->GetEntries();
	for(int i = 0; i < opt_entries; i++)
	{
		optree->GetEntry(i);
		std::cout << "[Optical Info] \t Event:   " << fEvent << std::endl;
		std::cout << "[Optical Info] \t Run:     " << fRun << std::endl;
		std::cout << "[Optical Info] \t PE:      " << fOpFlashPE << std::endl;
		std::cout << "[Optical Info] \t Time:    " << fOpFlashTime << std::endl;
		std::cout << "[Optical Info] \t WidthY:  " << fOpFlashWidthY << std::endl;
		std::cout << "[Optical Info] \t WidthZ:  " << fOpFlashWidthZ << std::endl;
		std::cout << "[Optical Info] \t CenterY: " << fOpFlashCenterY << std::endl;
		std::cout << "[Optical Info] \t CenterZ: " << fOpFlashCenterZ << std::endl;
	}

	return 0;
}//end out_inspect

#ifndef __ROOTCLING__


int main(){

	return out_inspect();
}

#endif
