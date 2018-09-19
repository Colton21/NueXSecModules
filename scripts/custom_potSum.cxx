#include <iostream>

//Root Includes
#include "TFile.h"
#include "TTree.h"

int main(int argc, char * argv[])
{

	const char * _file1 = argv[1];
	std::cout << "File: " << _file1 << std::endl;
	//first we need to open the root file
	TFile * f = new TFile(_file1);
	if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return 1; }
	TTree * mytree = (TTree*)f->Get("AnalyzeTPCO/pot_tree");
        TTree * pot_tree2 = (TTree*)f->Get("AnalyzeTPCO/pottree");

	double pot_sum = 0;
	double pot;
	mytree->SetBranchAddress("pot", &pot);

	for(int i = 0; i < mytree->GetEntries(); i++)
	{
		mytree->GetEntry(i);
		//std::cout << pot << std::endl;
		pot_sum = pot_sum + pot;
	}

	double pot_sum2 = 0;
	double pot2;
	pot_tree2->SetBranchAddress("pot", &pot2);

	for(int i = 0; i < pot_tree2->GetEntries(); i++)
	{
		pot_tree2->GetEntry(i);
		pot_sum2 += pot2;
	}

	std::cout << "Total POT (Method 1): " << pot_sum  << std::endl;
	std::cout << "Total POT (Method 2): " << pot_sum2 << std::endl;
	return 0;
}
