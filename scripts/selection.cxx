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
	for(auto const run : * _passed_runs) {std::cout << "In Time: " << run << std::endl; }

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
	std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
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
			std::cout << " \t " << i << " Passed!" << std::endl;
		}
		if(InFV == 0)//false
		{
			passed_tpco->at(i) = 0;
			std::cout << " \t " << i << " TPC Object Outside FV" << std::endl;
		}
	}//end tpc object loop
}

//***************************************************************************
//***************************************************************************

bool opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tollerance)
{
	const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
	if(distance <= tollerance) {return true; }
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
		largest_flash_v.push_back(fOpFlashCenterY);
		largest_flash_v.push_back(fOpFlashCenterZ);
		largest_flash_v_v->push_back(largest_flash_v);
		if(current_event != last_event) {largest_flash = 0; }
		largest_flash_v.clear();
	}
}


void flashRecoVtxDist(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                      double tollerance, std::vector<int> * passed_tpco)
{
	int n_tpc_obj = tpc_object_container_v->size();
	std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const bool is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), tollerance);
		if(is_close == 1)//true
		{
			passed_tpco->at(i) = 1;
			std::cout << " \t " << i << " Passed! " << std::endl;
		}
		if(is_close == 0)//false
		{
			passed_tpco->at(i) = 0;
			std::cout << " \t " << i << " TPC Object Vtx far from Flash" << std::endl;
		}
	}        //end tpc object loop
}//end flashRecoVtxDist

//***************************************************************************
//***************************************************************************

bool shwr_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z, double tollerance)
{
	const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );
	if(distance <= tollerance) {return true; }
	return false;
}

//***************************************************************************
//***************************************************************************

//this function wants to remove particles too far from the reconstructed neutrino vertex
void VtxNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                   double tollerance, std::vector<int> * passed_tpco)
{
	int n_tpc_obj = tpc_object_container_v->size();
	std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
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
			//check if at least one shower is within the tollerance to tpco vtx
			//if no shower within tollerance, discard tpco (i.e. it's not a cc nue)
			if(pfp_pdg == 11)        //if it's a shower
			{
				const double pfp_vtx_x = part.pfpVtxX();
				const double pfp_vtx_y = part.pfpVtxY();
				const double pfp_vtx_z = part.pfpVtxZ();
				close_shower = shwr_vtx_distance(tpc_vtx_x, tpc_vtx_y, tpc_vtx_z,
				                                 pfp_vtx_x, pfp_vtx_y, pfp_vtx_z, tollerance);
				if(close_shower == true)
				{
					std::cout << " \t " << i << " TPC Object Passed " << std::endl;
					passed_tpco->at(i) = 1;
					break;
				}
			}
		} //end loop pfps
		if(close_shower == 0)//false
		{
			std::cout << " \t " << i << " TPC Object Vtx far from All Shower Vtx" << std::endl;
			passed_tpco->at(i) = 0;
		}
	}//end tpc object loop
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
			std::cout << " \t " << i << " No Nue for TPC Object " << std::endl;
		}
		if(has_nue == true)
		{
			passed_tpco->at(i) = 1;
			std::cout << " \t " << i << " Has Nue for TPC Object" << std::endl;
		}
	}
}

//***************************************************************************
//***************************************************************************

int selection(){

	const char * _file1 = "../nue_xsec_extraction.root";
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
	mctree->SetBranchAddress("MC_PDG", &fMC_PDG);
	const int total_mc_entires = mctree->GetEntries();
	int mc_nue_counter = 0;
	for(int i = 0; i < total_mc_entires; i++)
	{
		if(fMC_PDG == 12) {mc_nue_counter++; }
	}
	std::cout << "MC Nue Counter: " << mc_nue_counter << std::endl;

	std::cout << "=====================" << std::endl;
	std::cout << "== Begin Selection ==" << std::endl;
	std::cout << "=====================" << std::endl;

	std::vector<int> * passed_runs = new std::vector<int>;
	std::cout << "Initialize Run Vector" << std::endl;

	//let's do the in-time cut as the very first thing
	std::cout << "=====================" << std::endl;
	std::cout << "==== In Time Cut ====" << std::endl;
	std::cout << "=====================" << std::endl;

	const int flash_pe_threshold = 50;
	const double flash_time_start = 5;
	const double flash_time_end = 16;
	loop_flashes(f, optree, flash_pe_threshold, flash_time_start,
	             flash_time_end, passed_runs);
	std::cout << "Passed Runs Size: " << passed_runs->size() << std::endl;

	//get vector with largest flashes y,z positions
	std::vector< std::vector< double> > * largest_flash_v_v = new std::vector < std::vector < double > >;
	SetXYflashVector(f, optree, largest_flash_v_v);
	std::cout << "Largest Flash Vector Size: " << largest_flash_v_v->size() << std::endl;

	//now let's do the TPCO related cuts

	const int total_entries = mytree->GetEntries();
	std::cout << "Total Events: " << total_entries << std::endl;
	for(int event = 0; event < total_entries; event++)
	{
		mytree->GetEntry(event);
		if(passed_runs->at(event) == 0) {continue; }//false

		//XY Position of largest flash
		std::vector < double > largest_flash_v = largest_flash_v_v->at(event);

		//List of TPC Objects which pass the cuts
		std::vector<int> * passed_tpco = new std::vector<int>;
		passed_tpco->resize(tpc_object_container_v->size(), 1);

		std::cout << "======================" << std::endl;
		std::cout << "==== Reco Nue Cut ====" << std::endl;
		std::cout << "======================" << std::endl;
		HasNue(tpc_object_container_v, passed_tpco);

		std::cout << "=======================" << std::endl;
		std::cout << "====== In FV Cut ======" << std::endl;
		std::cout << "=======================" << std::endl;

		const double _x1 = 0;
		const double _x2 = 0;
		const double _y1 = 0;
		const double _y2 = 0;
		const double _z1 = 0;
		const double _z2 = 0;
		fiducial_volume_cut(tpc_object_container_v, _x1, _x2, _y1, _y2, _z1, _z2, passed_tpco);

		std::cout << "============================" << std::endl;
		std::cout << "===== Vertex-To-Flash: =====" << std::endl;
		std::cout << "============================" << std::endl;
		const double tollerance = 100;//cm
		flashRecoVtxDist(largest_flash_v, tpc_object_container_v,
		                 tollerance, passed_tpco);


		std::cout << "============================" << std::endl;
		std::cout << "====== Shower-To-Nue: ======" << std::endl;
		std::cout << "============================" << std::endl;
		const double shwr_nue_tollerance = 50;//cm
		VtxNuDistance(tpc_object_container_v, shwr_nue_tollerance, passed_tpco);


		std::cout << "------------------" << std::endl;
		std::cout << "End Selection" << std::endl;
		std::cout << "------------------" << std::endl;
	}
	return 0;
}        //end selection


int main(){

	return selection();
}

#endif
