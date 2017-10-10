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

#ifndef __ROOTCLING__

//***************************************************************************
//***************************************************************************

int flash_in_time(int flash_pe, int flash_pe_threshold, double flash_time, double flash_start, double flash_end)
{
	if(flash_time >= flash_start && flash_time <= flash_end)
	{
		if(flash_pe >= flash_pe_threshold)
		{
			return 1;//true
		}
	}
	return 0;//false
}
//***************************************************************************
//***************************************************************************

void loop_flashes(TFile * f, TTree * optical_tree, int flash_pe_threshold, double flash_time_start,
                  double flash_time_end, std::vector<int> * _passed_runs)
{

	optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;

	optical_tree->SetBranchAddress("Run", &fRun);
	optical_tree->SetBranchAddress("Event", &fEvent);
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE);
	optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime);


	const int optical_entries = optical_tree->GetEntries();

	//int in_time = 0;//false
	/// I want passed_runs to be the same size as the number of TPCObject_v
	/// This will give the vector a 1:1 mapping to the events

	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;

	std::cout << "Optical Entries: " << optical_entries << std::endl;

	for(int i = 0; i < optical_entries; i++)
	{
		optical_tree->GetEntry(i);

		//this function here is meant to construct a vector mapping to which
		//events successfully pass this cut
		current_run = fRun;
		current_event = fEvent;

		int in_time = 0;//false
		if(_passed_runs->size() != 0)
		{
			if(current_event == last_event && current_run == last_run)
			{
				if(_passed_runs->back() == 0)//false
				{
					_passed_runs->pop_back();
				}
				if(_passed_runs->back() == 1)//true
				{
					last_event = current_event;
					last_run = current_run;
					continue;
				}
			}
		}
		last_event = current_event;
		last_run = current_run;

		in_time = flash_in_time(fOpFlashPE, flash_pe_threshold, fOpFlashTime, flash_time_start, flash_time_end);
		_passed_runs->push_back(in_time);

	}//end loop optical entries
	 //for(auto const run : * _passed_runs) {std::cout << "In Time: " << run << std::endl; }

}//end optical loop function

//***************************************************************************
//***************************************************************************

bool in_fv(double x, double y, double z,
           double _x1 =0, double _x2 =0, double _y1 =0,
           double _y2 =0, double _z1 =0, double _z2 =0 )
{
	double x1 = _x1;
	double x2 = _x2;
	double y1 = _y1;
	double y2 = _y2;
	double z1 = _z1;
	double z2 = _z2;

	const double det_x1 = 0;
	const double det_x2 = 256.35;
	const double det_y1 = -116.5;
	const double det_y2 = 116.5;
	const double det_z1 = 0;
	const double det_z2 = 1036.8;

	if(x <= det_x1 + x1 || x >= det_x2 - x2) {return false; }
	if(y <= det_y1 + y1 || y >= det_y2 - y2) {return false; }
	if(z <= det_z1 + z1 || z >= det_z2 - z2) {return false; }
	else{return true; }
}

//***************************************************************************
//***************************************************************************

void fiducial_volume_cut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                         std::vector<int> * passed_tpco)
{
	int n_tpc_obj = tpc_object_container_v->size();
	//std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
	if(passed_tpco->size() != n_tpc_obj) {std::cout << "Passed TPCO Vector Size != nTPCO!" << std::endl; }
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		//this vertex (tpcobject container) is the vertex of the pfp neutrino
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const bool InFV = in_fv(tpc_vtx_x, tpc_vtx_y, tpc_vtx_z, _x1, _x2, _y1, _y2, _z1, _z2);
		if(InFV == 1)//true
		{
			passed_tpco->at(i) = 1;
			std::cout << " \t " << i << "[Fid Volume Cut] \t Passed" << std::endl;
		}
		if(InFV == 0)//false
		{
			passed_tpco->at(i) = 0;
			//std::cout << " \t " << i << " TPC Object Outside FV" << std::endl;
		}
	}//end tpc object loop
}

//***************************************************************************
//***************************************************************************

bool opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance)
{
	const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
	if(distance <= tolerance) {return true; }
	return false;
}

//***************************************************************************
//***************************************************************************

void SetXYflashVector(TFile * f, TTree * optical_tree, std::vector< std::vector< double> > * largest_flash_v_v)
{
	optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");
	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashWidthY = 0;
	double fOpFlashWidthZ = 0;
	double fOpFlashCenterY = 0;
	double fOpFlashCenterZ = 0;

	optical_tree->SetBranchAddress("Run", &fRun);
	optical_tree->SetBranchAddress("Event", &fEvent);
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE);
	optical_tree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY);
	optical_tree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ);
	optical_tree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY);
	optical_tree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ);

	const int optical_entries = optical_tree->GetEntries();
	double largest_flash = 0;
	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;
	//since I only care about the largest flash in each event,
	//size(optical_entries)-->size(total_entries)
	for(int i = 0; i < optical_entries; i++)
	{
		optical_tree->GetEntry(i);
		current_run = fRun;
		current_event = fEvent;
		if(current_event != last_event) {largest_flash = 0; }
		double this_flash = fOpFlashPE;
		std::vector < double > largest_flash_v;//contains the y,z for largest flash
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
				largest_flash_v_v->pop_back();
			}
		}
		last_event = current_event;
		last_run = current_run;
		if(this_flash > largest_flash) {largest_flash = this_flash; }
		largest_flash_v.push_back(fOpFlashCenterY);
		largest_flash_v.push_back(fOpFlashCenterZ);
		largest_flash_v_v->push_back(largest_flash_v);
		largest_flash_v.clear();
	}
}

//***************************************************************************
//***************************************************************************

void flashRecoVtxDist(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                      double tolerance, std::vector<int> * passed_tpco)
{
	int n_tpc_obj = tpc_object_container_v->size();
	//std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const bool is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), tolerance);
		if(is_close == 1)//true
		{
			passed_tpco->at(i) = 1;
			std::cout << " \t " << i << "[Vertex-To-Flash] \t Passed " << std::endl;
		}
		if(is_close == 0)//false
		{
			passed_tpco->at(i) = 0;
			//std::cout << " \t " << i << " TPC Object Vtx far from Flash" << std::endl;
		}
	}        //end tpc object loop
}//end flashRecoVtxDist

//***************************************************************************
//***************************************************************************

bool shwr_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z, double tolerance)
{
	const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );
	if(distance <= tolerance) {return true; }
	return false;
}

//***************************************************************************
//***************************************************************************

//this function wants to remove particles too far from the reconstructed neutrino vertex
void VtxNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                   double tolerance, std::vector<int> * passed_tpco)
{
	int n_tpc_obj = tpc_object_container_v->size();
	//std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		bool close_shower = false;
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			//check if at least one shower is within the tolerance to tpco vtx
			//if no shower within tolerance, discard tpco (i.e. it's not a cc nue)
			if(pfp_pdg == 11)        //if it's a shower
			{
				const double pfp_vtx_x = part.pfpVtxX();
				const double pfp_vtx_y = part.pfpVtxY();
				const double pfp_vtx_z = part.pfpVtxZ();
				close_shower = shwr_vtx_distance(tpc_vtx_x, tpc_vtx_y, tpc_vtx_z,
				                                 pfp_vtx_x, pfp_vtx_y, pfp_vtx_z, tolerance);
				if(close_shower == true)
				{
					std::cout << " \t " << i << "[Shower-To-Nue] \t Passed " << std::endl;
					passed_tpco->at(i) = 1;
					break;
				}
			}
		} //end loop pfps
		if(close_shower == 0)//false
		{
			//std::cout << " \t " << i << " TPC Object Vtx far from All Shower Vtx" << std::endl;
			passed_tpco->at(i) = 0;
		}
	}//end tpc object loop
}

//***************************************************************************
//***************************************************************************

//this function wants to remove particles too far from the reconstructed neutrino vertex
void HitThreshold(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                  double threshold, std::vector<int> * passed_tpco)
{
	int n_tpc_obj = tpc_object_container_v->size();
	//std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		bool over_threshold = false;
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_tpco_hits = tpc_obj.NumPFPHits();
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			//check if at least one shower is within the tolerance to tpco vtx
			//if no shower within tolerance, discard tpco (i.e. it's not a cc nue)
			if(pfp_pdg == 11)        //if it's a shower
			{
				const int n_pfp_hits = part.NumPFPHits();
				if(n_pfp_hits >= threshold) {over_threshold = true; }
				if(over_threshold == true)
				{
					std::cout << " \t " << i << "[Hit Threshold] \t Passed " << std::endl;
					passed_tpco->at(i) = 1;
					break;
				}
			}
		} //end loop pfps
		if(over_threshold == 0)//false
		{
			//std::cout << " \t " << i << " No Showers over threshold" << std::endl;
			passed_tpco->at(i) = 0;
		}
	}//end tpc object loop
}

//***************************************************************************
//***************************************************************************
//this gives a list of all of the origins of the tpc objects
void GetOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::string> * tpco_origin_v)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		auto const tpc_obj = tpc_object_container_v->at(i);
		const std::string tpc_obj_origin = tpc_obj.Origin();
		tpco_origin_v->push_back(tpc_obj_origin);
	}
}

//***************************************************************************
//***************************************************************************
//this function simply checks if the tpc object is a nue
void HasNue(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();

		bool has_nue = false;
		for(int j = 0; j <n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			if(pfp_pdg == 12)
			{
				has_nue = true;
				break;
			}
		}//end loop pfparticles
		if(has_nue == false)
		{
			passed_tpco->at(i) = 0;
			//std::cout << " \t " << i << " No Nue for TPC Object " << std::endl;
		}
		if(has_nue == true)
		{
			passed_tpco->at(i) = 1;
			std::cout << " \t " << i << "[Reco Nue Cut] \t Passed" << std::endl;
		}
	}
}

//***************************************************************************
//***************************************************************************
//this function will help us check how each cut is performing
bool ValidTPCObjects(std::vector<int> * passed_tpco)
{
	const int n_tpco = passed_tpco->size();
	int pass_sum = 0;
	for(int const tpco : * passed_tpco)
	{
		pass_sum = pass_sum + tpco;
	}
	//0 --> failed tpco, 1 --> passed tpco
	if(pass_sum == 0) {return false; }
	return true;
}
//***************************************************************************
//***************************************************************************

std::vector<int> TabulateOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco)
{
	int nue_cc = 0;
	int cosmic = 0;
	int nue_nc = 0;
	int numu   = 0;
	std::vector<int> tabulated_origins;
	tabulated_origins.resize(4);

	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const std::string tpc_obj_origin = tpc_obj.Origin();
		std::cout << "\t " << i << " " << tpc_obj_origin << std::endl;
		//const int tpc_obj_pdg = tpc_obj.MCParticlePdgCode();
		if(tpc_obj_origin == "kCosmicRay") {cosmic++; }
		//if(tpc_obj_origin == "kBeamNeutrino") {nue_cc++; } //this is for the time being
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" ) {}
		}
	}
	tabulated_origins.at(0) = nue_cc;
	tabulated_origins.at(1) = cosmic;
	tabulated_origins.at(2) = nue_nc;
	tabulated_origins.at(3) = numu;
	return tabulated_origins;
}

//***************************************************************************
//***************************************************************************
//modify this so it takes a string of the cut name so I only pass it a few variable at a time,
//then I can call this function several times later at the bottom
void PrintInfo(int mc_nue_counter,
               int counter,
               int counter_nue_cc,
               int counter_cosmic,
               std::string cut_name)
{
	std::cout << " <" << cut_name << "> " << std::endl;
	std::cout << "Total Reco Nue    : " << counter << std::endl;
	std::cout << " Number of Nue CC : " << counter_nue_cc << std::endl;
	// std::cout << " Number of Nue NC : " << << std::endl;
	// std::cout << " Number of Numu   : " << << std::endl;
	std::cout << " Number of Cosmic : " << counter_cosmic << std::endl;
	// std::cout << " kBeamNeutrino    : " << << std::endl;
	// std::cout << " kCosmicRay       : " << << std::endl;
	// std::cout << " kMixed/Unknown   : " << << std::endl;
	std::cout << "------------------------" << std::endl;
	const double efficiency = counter_nue_cc / mc_nue_counter;
	const double purity = counter_nue_cc / counter;
	std::cout << " Efficiency       : " << efficiency << std::endl;
	std::cout << " Purity           : " << purity << std::endl;
	std::cout << "------------------------" << std::endl;
	std::cout << "------------------------" << std::endl;
}

//***************************************************************************
//***************************************************************************

int selection(){

	const char * _file1 = "../nue_xsec_extraction.root";
	//const char * _file1 = "../cosmic_extraction.root";
	std::cout << "File Path: " << _file1 << std::endl;
	//first we need to open the root file
	TFile * f = new TFile(_file1);
	if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return 1; }
	TTree * mytree = (TTree*)f->Get("AnalyzeTPCO/tree");
	TTree * optree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");
	TTree * mctree = (TTree*)f->Get("AnalyzeTPCO/mcparticle_tree");

	std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
	mytree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

	int fMC_PDG = 0;
	int fMCOrigin = -1;
	int fMCMother = 0;
	mctree->SetBranchAddress("MC_PDG", &fMC_PDG);
	mctree->SetBranchAddress("MC_Origin", &fMCOrigin);
	mctree->SetBranchAddress("MC_Mother", &fMCMother);
	const int total_mc_entires = mctree->GetEntries();
	std::cout << "Total MC Entries: " << total_mc_entires << std::endl;
	int mc_nue_counter = 0;
	for(int i = 0; i < total_mc_entires; i++)
	{
		mctree->GetEntry(i);
		//for now we'll just count the nue cc interactions - primary, beam electrons
		if(fMC_PDG == 11 && fMCMother == 0 && fMCOrigin == 0) {mc_nue_counter++; }
	}
	std::cout << "MC Nue Counter: " << mc_nue_counter << std::endl;

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

	const int flash_pe_threshold = 50;
	const double flash_time_start = 5;
	const double flash_time_end = 16;
	loop_flashes(f, optree, flash_pe_threshold, flash_time_start,
	             flash_time_end, passed_runs);

	//get vector with largest flashes y,z positions
	std::vector< std::vector< double> > * largest_flash_v_v = new std::vector < std::vector < double > >;
	SetXYflashVector(f, optree, largest_flash_v_v);
	std::cout << "Largest Flash Vector Size: " << largest_flash_v_v->size() << std::endl;

	int run_sum = 0;
	for(auto const run : * passed_runs) {run_sum = run_sum + run; }
	std::cout << "Passed Runs Size: " << run_sum << std::endl;

	int reco_nue_counter = 0;
	int reco_nue_counter_nue_cc = 0;
	int reco_nue_counter_cosmic = 0;
	int in_fv_counter = 0;
	int in_fv_counter_nue_cc = 0;
	int in_fv_counter_cosmic = 0;
	int vtx_flash_counter = 0;
	int vtx_flash_counter_nue_cc = 0;
	int vtx_flash_counter_cosmic = 0;
	int shwr_tpco_counter = 0;
	int shwr_tpco_counter_nue_cc = 0;
	int shwr_tpco_counter_cosmic = 0;
	int hit_threshold_counter = 0;
	int hit_threshold_counter_nue_cc = 0;
	int hit_threshold_counter_cosmic = 0;
	std::vector<int> tabulated_origins;

	//**********************************
	//now let's do the TPCO related cuts
	//*********************************
	for(int event = 0; event < total_entries; event++)
	{
		std::cout << "----------------------" << std::endl;
		std::cout << "[EVENT NUMBER] \t " << event << std::endl;
		std::cout << "----------------------" << std::endl;

		mytree->GetEntry(event);
		if(passed_runs->at(event) == 0)
		{
			std::cout << "[Failed In-Time Cut]" << std::endl;
			continue;
		}//false

		std::vector<std::string> *tpco_origin_v = new std::vector<std::string>;
		GetOrigins(tpc_object_container_v, tpco_origin_v);

		//XY Position of largest flash
		std::vector < double > largest_flash_v = largest_flash_v_v->at(event);

		//List of TPC Objects which pass the cuts
		std::vector<int> * passed_tpco = new std::vector<int>;
		passed_tpco->resize(tpc_object_container_v->size(), 1);

		//** start the cuts here **

		//reco nue cut
		HasNue(tpc_object_container_v, passed_tpco);
		if(ValidTPCObjects(passed_tpco) == false) {continue; }
		if(ValidTPCObjects(passed_tpco) == true ) {reco_nue_counter++; }
		tabulated_origins = TabulateOrigins(tpc_object_container_v, passed_tpco);
		reco_nue_counter_nue_cc = reco_nue_counter_nue_cc + tabulated_origins.at(0);
		reco_nue_counter_cosmic = reco_nue_counter_cosmic + tabulated_origins.at(1);

		//in fv cut
		const double _x1 = 0;
		const double _x2 = 0;
		const double _y1 = 0;
		const double _y2 = 0;
		const double _z1 = 0;
		const double _z2 = 0;
		fiducial_volume_cut(tpc_object_container_v, _x1, _x2, _y1, _y2, _z1, _z2, passed_tpco);
		if(ValidTPCObjects(passed_tpco) == false) {continue; }
		if(ValidTPCObjects(passed_tpco) == true ) {in_fv_counter++; }
		tabulated_origins = TabulateOrigins(tpc_object_container_v, passed_tpco);
		in_fv_counter_nue_cc = reco_nue_counter_nue_cc + tabulated_origins.at(0);
		in_fv_counter_cosmic = reco_nue_counter_cosmic + tabulated_origins.at(1);

		//vertex to flash
		const double tolerance = 100;//cm
		flashRecoVtxDist(largest_flash_v, tpc_object_container_v,
		                 tolerance, passed_tpco);
		if(ValidTPCObjects(passed_tpco) == false) {continue; }
		if(ValidTPCObjects(passed_tpco) == true) {vtx_flash_counter++; }
		tabulated_origins = TabulateOrigins(tpc_object_container_v, passed_tpco);
		vtx_flash_counter_nue_cc = reco_nue_counter_nue_cc + tabulated_origins.at(0);
		vtx_flash_counter_cosmic = reco_nue_counter_cosmic + tabulated_origins.at(1);

		//distance between pfp shower and nue object
		const double shwr_nue_tolerance = 50;//cm
		VtxNuDistance(tpc_object_container_v, shwr_nue_tolerance, passed_tpco);
		if(ValidTPCObjects(passed_tpco) == false) {continue; }
		if(ValidTPCObjects(passed_tpco) == true ) {shwr_tpco_counter++; }
		tabulated_origins = TabulateOrigins(tpc_object_container_v, passed_tpco);
		shwr_tpco_counter_nue_cc = reco_nue_counter_nue_cc + tabulated_origins.at(0);
		shwr_tpco_counter_cosmic = reco_nue_counter_cosmic + tabulated_origins.at(1);

		//hit threshold for showers
		const double shwr_hit_threshold = 50;//hits
		HitThreshold(tpc_object_container_v, shwr_hit_threshold, passed_tpco);
		if(ValidTPCObjects(passed_tpco) == false) {continue; }
		if(ValidTPCObjects(passed_tpco) == true ) {hit_threshold_counter++; }
		tabulated_origins = TabulateOrigins(tpc_object_container_v, passed_tpco);
		hit_threshold_counter_nue_cc = reco_nue_counter_nue_cc + tabulated_origins.at(0);
		hit_threshold_counter_cosmic = reco_nue_counter_cosmic + tabulated_origins.at(1);
	}
	std::cout << "------------------" << std::endl;
	std::cout << "End Selection" << std::endl;
	std::cout << "------------------" << std::endl;

	//we also want some metrics to print at the end
	PrintInfo( mc_nue_counter,
	           reco_nue_counter,
	           reco_nue_counter_nue_cc,
	           reco_nue_counter_cosmic,
	           "Reco Nue"
	           );
	PrintInfo( mc_nue_counter,
	           in_fv_counter,
	           in_fv_counter_nue_cc,
	           in_fv_counter_cosmic,
	           "In FV"
	           );
	PrintInfo( mc_nue_counter,
	           vtx_flash_counter,
	           vtx_flash_counter_nue_cc,
	           vtx_flash_counter_cosmic,
	           "Vtx-to-Flash"
	           );
	PrintInfo( mc_nue_counter,
	           shwr_tpco_counter,
	           shwr_tpco_counter_nue_cc,
	           shwr_tpco_counter_cosmic,
	           "Shower-to-TPCO"
	           );
	PrintInfo( mc_nue_counter,
	           hit_threshold_counter,
	           hit_threshold_counter_nue_cc,
	           hit_threshold_counter_cosmic,
	           "Hit Threshold"
	           );



	return 0;
}        //end selection


int main(){

	return selection();
}

#endif
