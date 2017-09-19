//#include "../xsecAna/TPCObject.h"
#include "../xsecAna/TpcObjectContainer.h"
#include "../xsecAna/ParticleContainer.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

int main()
{
	const char * _file1 = "../nue_xsec_extraction.root";
	std::cout << _file1 << std::endl;
	//first we need to open the root file
	TFile * f = new TFile(_file1);
	if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return 1; }
	TTree * mytree = (TTree*)f->Get("AnalyzeTPCO/tree");

	const std::vector<xsecAna::TPCObjectContainer> tpc_object_container_v;
	mytree->SetBranchAddress("TpcObject", &tpc_object_container_v);

	const int total_entries = mytree->GetEntries();

	for(int i = 0; i < total_entries; i++)
	{
		mytree->GetEntry(i);
		std::cout << tpc_object_container_v.size() << std::endl;
	}

	return 0;
}
