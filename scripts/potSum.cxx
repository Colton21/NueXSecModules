#include <iostream>

//Root Includes
#include "TFile.h"
#include "TTree.h"

int main(int argc, char *argv[])
{

  const char * _file1 = argv[1];
  std::cout << "File: " << _file1 << std::endl;
  //first we need to open the root file
  TFile * f = new TFile(_file1);
  if(!f->IsOpen()){std::cout << "Could not open file!" << std::endl; return 1;}
  TTree * mytree = (TTree*)f->Get("POTFinder/pot_tree");

  double pot_sum = 0;
  double pot;
  mytree->SetBranchAddress("pot", &pot);


  for(int i = 0; i < mytree->GetEntries(); i++)
  {
    mytree->GetEntry(i);
    //std::cout << pot << std::endl;
    pot_sum = pot_sum + pot;
  }

  std::cout << "Total POT: " << pot_sum << std::endl;
  return 0;
}
