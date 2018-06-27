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

	double pot_sum = 0;
	double pot;
	mytree->SetBranchAddress("pot", &pot);


	for(int i = 0; i < mytree->GetEntries(); i++)
	{
		mytree->GetEntry(i);
		//std::cout << pot << std::endl;
		pot_sum = pot_sum + pot;
	}

	double pot_sum_2 = 0;
	/*
	TTree * mytree_2 = (TTree*)f->Get("AnalyzeTPCO/pottree");
	double pot_2;
	mytree_2->SetBranchAddress("pot", &pot_2);
	for(int i = 0; i < mytree_2->GetEntries(); i++)
	{
		mytree_2->GetEntry(i);
		pot_sum_2 = pot_sum_2 + pot_2;
	}
	*/

	std::cout << "Total POT:     " << pot_sum << std::endl;
	std::cout << "Alt Total POT: " << pot_sum_2 << std::endl;
	return 0;
}
