#include <iostream>

#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"

//let's find the integrated flux!

int flux_calc()
{
  //for printing the bins
  bool _verbose = true;

  std::cout << "Opening NuMIFlux.root" << std::endl;
  //first we need to open the root file
  TFile * f = new TFile("NuMIFlux.root");
  if(!f->IsOpen()){std::cout << "Could not open file!" << std::endl; return 1;}
  TH1D * h_nue_flux = (TH1D*)f->Get("nueFluxHisto");
  
  std::cout << "Opening argon_xsec_nue.root" << std::endl;
  //this is a genie file with the cross section
  TFile * xsec_f = new TFile("argon_xsec_nue.root");
  if(!xsec_f->IsOpen()){std::cout << "Could not open file!" << std::endl; return 1;}
  TGraph * g_nue_xsec = (TGraph*)xsec_f->Get("nu_e_Ar40/tot_cc");

  const int n_bins = h_nue_flux->GetNbinsX();
  double bin_flux_sum = 0;
  double bin_interaction_sum = 0;

  std::cout << "Loop over energy bins of 50 MeV" << std::endl;
  for(int bin = 0; bin < n_bins; bin++)
  {
    const double bin_flux = h_nue_flux->GetBinContent(bin);
    bin_flux_sum = bin_flux_sum + bin_flux;

    const double bin_energy = h_nue_flux->GetBinCenter(bin);
    const double bin_xsec_val = g_nue_xsec->Eval(bin_energy) * 1e-38 / 40;//this gets the units per nucleon per cm^2
    const double bin_interactions = bin_xsec_val * bin_flux;
    bin_interaction_sum = bin_interaction_sum + bin_interactions;

    if(_verbose == true && bin_flux > 0.0)
    {
      std::cout << "Cross Section: " << bin_xsec_val << std::endl;
      std::cout << "Bin Energy: " << bin_energy << '\t' << "Bin Flux: " << bin_flux << std::endl;
      std::cout << "===================" << std::endl;
      std::cout << "Interactions: " << bin_interactions << std::endl;
      std::cout << "===================" << std::endl;
    }
  }

  std::cout << "Total Neutrinos: " << bin_flux_sum << std::endl;
  std::cout << "Total Interactions: " << bin_interaction_sum << std::endl;
  //the website: https://cdcvs.fnal.gov/redmine/projects/ubooneoffline/wiki/NuMI_Flux_Histograms
  //states that the flux is scaled to 6e20 POT
  const double pot_used = 6 * pow(10, 20);
  std::cout << "Total POT used: " << pot_used << std::endl;

  std::cout << "Final Result is: " << bin_flux_sum / pot_used << " nues per POT per cm^2" << std::endl;
  std::cout << "SUM{Phi(E) * Sigma(E)} / SUM{Phi(E)} = " << bin_interaction_sum / bin_flux_sum << " cm^2" << std::endl;

  return 0;

}
