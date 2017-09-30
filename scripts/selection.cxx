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

bool flash_in_time(int flash_pe, int flash_pe_threshold, double flash_time, double flash_start, double flash_end)
{
	if(flash_time >= flash_start && flash_time <= flash_end)
	{
		if(flash_pe >= flash_pe_threshold)
		{
			return true;
		}
	}
	return false;
}
//***************************************************************************
//***************************************************************************

void loop_flashes(TFile * f, TTree * optical_tree, int flash_pe_threshold, double flash_time_start,
                  double flash_time_end, std::vector<bool> * passed_runs){

	optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;
	double fOpFlashWidthY = 0;
	double fOpFlashWidthZ = 0;
	double fOpFlashCenterY = 0;
	double fOpFlashCenterZ = 0;
	optical_tree->SetBranchAddress("Run", &fRun);
	optical_tree->SetBranchAddress("Event", &fEvent);
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE);
	optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime);
	optical_tree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY);
	optical_tree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ);
	optical_tree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY);
	optical_tree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ);

	const int optical_entries = optical_tree->GetEntries();

	bool in_time = false;
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
		bool in_time = false;
		if(passed_runs->size() != 0)
		{
			if(current_event == last_event && current_run == last_run)
			{
				if(passed_runs->back() == false)
				{
					passed_runs->pop_back();
				}
				if(passed_runs->back() == true)
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
		passed_runs->push_back(in_time);

	}//end loop optical entries
	for(auto const run : *passed_runs) {std::cout << "In Time: " << run << std::endl; }

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

void fiducial_volume_cut(TTree * mytree, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                         std::vector<bool> *passed_runs)
{
	const int total_entries = mytree->GetEntries();
	std::cout << "Total Events: " << total_entries << std::endl;
	for(int i = 0; i < total_entries; i++)
	{
		mytree->GetEntry(i);
		//if(passed_runs.at(i) == false) {continue; }
		int n_tpc_obj = tpc_object_container_v->size();
		std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
		for(int i = 0; i < n_tpc_obj; i++)
		{
			auto const tpc_obj = tpc_object_container_v->at(i);
			//this vertex (tpcobject container) is the vertex of the pfp neutrino
			const double tpc_vtx_x = tpc_obj.pfpVtxX();
			const double tpc_vtx_y = tpc_obj.pfpVtxY();
			const double tpc_vtx_z = tpc_obj.pfpVtxZ();
			const bool InFV = in_fv(tpc_vtx_x, tpc_vtx_y, tpc_vtx_z, _x1, _x2, _y1, _y2, _z1, _z2);
			//if(InFV == true) {std::cout << " \t Passed!" << std::endl; }
			if(InFV == false)
			{
				//std::cout << " \t Failed!" << std::endl;
				tpc_object_container_v->erase(tpc_object_container_v->begin() + i);
				n_tpc_obj = tpc_object_container_v->size();
			}
		}//end tpc object loop
		 //all of the tpc objects were outside the FV
		if(tpc_object_container_v->size() == 0)
		{
			bool *state = passed_runs->data()[i];
			//passed_runs[i] = false;
		}
	}//end events loop
}



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

	std::cout << "===================" << std::endl;
	std::cout << "Begin Selection" << std::endl;
	std::cout << "===================" << std::endl;


	std::cout << "===================" << std::endl;
	std::cout << "In Time Cut " << std::endl;
	std::cout << "===================" << std::endl;

	const int flash_pe_threshold = 50;
	const double flash_time_start = 5;
	const double flash_time_end = 16;
	std::vector<bool> *passed_runs;
	loop_flashes(f, optree, flash_pe_threshold, flash_time_start,
	             flash_time_end, passed_runs);
	std::cout << "Passed Runs Size: " << passed_runs->size() << std::endl;


	std::cout << "===================" << std::endl;
	std::cout << "In FV Cut " << std::endl;
	std::cout << "===================" << std::endl;

	const double _x1 = 0;
	const double _x2 = 0;
	const double _y1 = 0;
	const double _y2 = 0;
	const double _z1 = 0;
	const double _z2 = 0;
	fiducial_volume_cut(mytree, tpc_object_container_v, _x1, _x2, _y1, _y2, _z1, _z2, passed_runs);


	return 0;
}        //end selection


int main(){

	return selection();
}

#endif
