#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "nusimdata/SimulationBase/MCFlux.h"

#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

namespace mcfluxfinder {
  class MCFluxFinder;
}

class mcfluxfinder::MCFluxFinder : public art::EDAnalyzer {

public:
  explicit MCFluxFinder (fhicl::ParameterSet const & p);


  // Plugins should not be copied or assigned.
  MCFluxFinder(MCFluxFinder const &) = delete;
  MCFluxFinder(MCFluxFinder &&) = delete;
  MCFluxFinder & operator = (MCFluxFinder const &) = delete;
  MCFluxFinder & operator = (MCFluxFinder &&) = delete;

  //Required functions.
  void analyze(art::Event const & e) override;
  //void endSubRun(art::SubRun const & sr) override;

  void beginJob() override;
  void endJob() override;
private:

  const double dummy = 0;

  TH1D * h_theta = new TH1D ("h_theta", "h_theta", 40, 0, 80);
  TH1D * h_momentum = new TH1D ("h_momentum", "h_momentum", 40, 0, 20);
  TH2D * h_momentum_theta = new TH2D ("h_momentum_theta", "h_momentum_theta", 40, 0, 16, 40, 0, 80);
  TH2D * h_momentum_theta_piplus = new TH2D ("h_momentum_theta_pi+", "h_momentum_theta_pi+", 40, 0, 16, 40, 0, 80);
  TH2D * h_momentum_theta_K0L = new TH2D ("h_momentum_theta_K0L", "h_momentum_theta_K0L", 40, 0, 16, 40, 0, 80);
  TH2D * h_momentum_theta_p = new TH2D ("h_momentum_theta_p", "h_momentum_theta_p", 40, 0, 16, 40, 0, 80);
  TH1I * h_particle_type = new TH1I ("h_particle_type", "h_particle_type", 15, 0, 15);
  TH1I * h_particle_type_nue = new TH1I ("h_particle_type_nue", "h_particle_type_nue", 15, 0, 15);
  TH1I * h_particle_type_nuebar = new TH1I ("h_particle_type_nuebar", "h_particle_type_nuebar", 15, 0, 15);
  TH1I * h_particle_type_numu = new TH1I ("h_particle_type_numu", "h_particle_type_numu", 15, 0, 15);
  TH1I * h_particle_type_numubar = new TH1I ("h_particle_type_numubar", "h_particle_type_numubar", 15, 0, 15);

};//end class def

void mcfluxfinder::MCFluxFinder::beginJob()
{
  std::cout << "--- Begin Job --- " << dummy << std::endl;
  // Implementation of optional member function here.
  //art::ServiceHandle< art::TFileService > tfs;
  //define trees
  //pot_tree = tfs->make<TTree>("pot_tree", "pot_per_subrun");

  //pot_tree->Branch("pot", &pot, "pot/D");
}

mcfluxfinder::MCFluxFinder::MCFluxFinder(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
}

void mcfluxfinder::MCFluxFinder::analyze(art::Event const & e) {

  art::Handle < std::vector < simb::MCTruth > > MCTruthHandle;
  e.getByLabel("generator", MCTruthHandle);
  if(!MCTruthHandle.isValid()) {std::cout << "MCTruth Handle is not valid" << std::endl; exit(1); }

  bool nue_cc = false;
  bool nue_nc = false;
  bool nuebar_cc = false; 
  bool nuebar_nc = false;
  bool numu_cc = false;
  bool numu_nc = false;
  bool numubar_cc = false;
  bool numubar_nc = false;

  for(auto const & mctruth : (*MCTruthHandle) )
  {
    auto const mc_nu = mctruth.GetNeutrino().Nu();
    if(mc_nu.PdgCode() ==  12)
    {
      if(mctruth.GetNeutrino().CCNC() == 0){nue_cc = true;}
      if(mctruth.GetNeutrino().CCNC() == 1){nue_nc = true;}
    }
    if(mc_nu.PdgCode() == -12)
    {
      if(mctruth.GetNeutrino().CCNC() == 0){nuebar_cc = true;}
      if(mctruth.GetNeutrino().CCNC() == 1){nuebar_nc = true;}
    }
    if(mc_nu.PdgCode() == 14)
    {
      if(mctruth.GetNeutrino().CCNC() == 0){numu_cc = true;}
      if(mctruth.GetNeutrino().CCNC() == 1){numu_nc = true;}
    }
    if(mc_nu.PdgCode() == -14)
    {
      if(mctruth.GetNeutrino().CCNC() == 0){numubar_cc = true;}
      if(mctruth.GetNeutrino().CCNC() == 1){numubar_nc = true;}
    }
  }

  int mc_flux_pdg = 0;
  int mc_flux_pid = 0;
  double mc_flux_px = 0;
  double mc_flux_py = 0;
  double mc_flux_pz = 0;


  double angle = 0;
  double momentum = 0;

  art::Handle < std::vector < simb::MCFlux > > MCFluxHandle;
  e.getByLabel("generator", MCFluxHandle);
  if(!MCFluxHandle.isValid()) {std::cout << "MCFlux Handle is not valid!" << std::endl; exit(1); }
  for(auto const & mcflux : (*MCFluxHandle) )
  {
   //These pdg codes should correspond to the particles produced immediately after the proton
   //hits the target, not necissarily the nu-parent
   mc_flux_pdg = mcflux.ftgptype;

   //this is what I've seen so far...
   if(mc_flux_pdg == 13){mc_flux_pid = 1;} //muon
   if(mc_flux_pdg == -13){mc_flux_pid = 2;} //mu+
   if(mc_flux_pdg == 211){mc_flux_pid = 3;} //pi+
   if(mc_flux_pdg == -211){mc_flux_pid = 4;} //pi-
   if(mc_flux_pdg == 130){mc_flux_pid = 5;} //K0_L
   if(mc_flux_pdg == 310){mc_flux_pid = 6;} //K0_S
   if(mc_flux_pdg == 311){mc_flux_pid = 7;} //K0
   if(mc_flux_pdg == 321){mc_flux_pid = 8;} //K+
   if(mc_flux_pdg == -321){mc_flux_pid = 9;} //K-
   if(mc_flux_pdg == 2212){mc_flux_pid = 10;} //p
   if(mc_flux_pdg == 2112){mc_flux_pid = 11;} //n
   if(mc_flux_pdg == 3122){mc_flux_pid = 12;} //lambda
   h_particle_type->Fill(mc_flux_pid);

   if(nue_cc || nue_nc){h_particle_type_nue->Fill(mc_flux_pid);}
   if(nuebar_cc || nuebar_nc){h_particle_type_nuebar->Fill(mc_flux_pid);}
   if(numu_cc || numu_nc){h_particle_type_numu->Fill(mc_flux_pid);}
   if(numubar_cc || numu_nc){h_particle_type_numubar->Fill(mc_flux_pid);}


   //I think this momentum is with respect to the beam direction
   /*
   mc_flux_px = mcflux.ftgppx;
   mc_flux_py = mcflux.ftgppy;
   mc_flux_pz = mcflux.ftgppz;
   */
   //these variables are filled in flux file, above are not
   //they should be equivalent...
   mc_flux_px = mcflux.fppdxdz * mcflux.fpppz;
   mc_flux_py = mcflux.fppdydz * mcflux.fpppz;
   mc_flux_pz = mcflux.fpppz;
   momentum = sqrt((mc_flux_px * mc_flux_px) + (mc_flux_py * mc_flux_py) + (mc_flux_pz * mc_flux_pz));
   //create a unit vector in the Z direction - assuming z is the beam
   std::vector<double> mc_flux_unit;
   mc_flux_unit.push_back(0);
   mc_flux_unit.push_back(0);
   mc_flux_unit.push_back(1);

   //now calculate the angle between the unit vector and the momentum vector
   const double dot_prod = (mc_flux_unit.at(0) * mc_flux_px) + (mc_flux_unit.at(1) * mc_flux_py) + (mc_flux_unit.at(2) * mc_flux_pz);
   angle = acos(dot_prod / (1 * std::abs(momentum) ) ) * (180 / 3.141592);

   
   std::cout << "Particle Type: " << mc_flux_pdg << std::endl; 
   /*
   std::cout << "Momentum: " << mc_flux_px << ", " << mc_flux_py << ", " << mc_flux_pz << ", " << momentum << std::endl;
   std::cout << "Angle off-target: " << angle << std::endl;
   */

   h_momentum->Fill(std::abs(momentum));
   h_theta->Fill(angle);
   h_momentum_theta->Fill(std::abs(momentum), angle);
   if(mc_flux_pdg == 211){h_momentum_theta_piplus->Fill(std::abs(momentum), angle);}
   if(mc_flux_pdg == 2212){h_momentum_theta_p->Fill(std::abs(momentum), angle);}
   if(mc_flux_pdg == 130){h_momentum_theta_K0L->Fill(std::abs(momentum), angle);}

  }

  if(nue_cc == true){std::cout << "Nue CC" << std::endl;}
  if(nue_nc == true){std::cout << "Nue NC" << std::endl;}
  if(nuebar_cc == true){std::cout << "NueBar CC" << std::endl;}
  if(nuebar_nc == true){std::cout << "NueBar NC" << std::endl;}

  if(numu_cc == true){std::cout << "Numu CC" << std::endl;}
  if(numu_nc == true){std::cout << "Numu NC" << std::endl;}
  if(numubar_cc == true){std::cout << "NumuBar CC" << std::endl;}
  if(numubar_nc == true){std::cout << "NumuBar NC" << std::endl;}


}//end analyzer loop

void mcfluxfinder::MCFluxFinder::endJob()
{

  TCanvas * c1 = new TCanvas();
  c1->cd();
  h_momentum->Draw();
  c1->Print("plots/h_production_momentum.pdf");


  TCanvas * c2 = new TCanvas();
  c2->cd();
  h_theta->Draw();
  c2->Print("plots/h_production_theta.pdf");

  TCanvas * c3 = new TCanvas();
  c3->cd();
  h_momentum_theta->Draw("colz");
  h_momentum_theta->GetXaxis()->SetTitle("Momentum [GeV]");
  h_momentum_theta->GetYaxis()->SetTitle("#theta [Degrees]");
  c3->Print("plots/h_production_momentum_theta.pdf");

  TCanvas * c3b = new TCanvas();
  c3b->cd();
  h_momentum_theta_piplus->Draw("colz");
  h_momentum_theta_piplus->GetXaxis()->SetTitle("Momentum #pi+ [GeV]");
  h_momentum_theta_piplus->GetYaxis()->SetTitle("#theta #pi+ [Degrees]");
  c3b->Print("plots/h_production_momentum_theta_pi+.pdf");

  TCanvas * c3c = new TCanvas();
  c3c->cd();
  h_momentum_theta_K0L->Draw("colz");
  h_momentum_theta_K0L->GetXaxis()->SetTitle("Momentum K0L [GeV]");
  h_momentum_theta_K0L->GetYaxis()->SetTitle("#theta K0L [Degrees]");
  c3c->Print("plots/h_production_momentum_theta_K0L.pdf");

  TCanvas * c3d = new TCanvas();
  c3d->cd();
  h_momentum_theta_p->Draw("colz");
  h_momentum_theta_p->GetXaxis()->SetTitle("Momentum p [GeV]");
  h_momentum_theta_p->GetYaxis()->SetTitle("#theta p [Degrees]");
  c3d->Print("plots/h_production_momentum_theta_p.pdf");


  const char * str_particle[15] = {"0", "#mu-", "#mu+", "#pi+", "#pi-", "K0L", 
                                   "K0S", "K0", "K+", "K-", "p", "n", "#Lambda", "-", "-"};

  for(int i = 1; i <= 15; i++)
  {
    h_particle_type->GetXaxis()->SetBinLabel(i, str_particle[i-1]);
    h_particle_type_nue->GetXaxis()->SetBinLabel(i, str_particle[i-1]);
    h_particle_type_nuebar->GetXaxis()->SetBinLabel(i, str_particle[i-1]);
    h_particle_type_numu->GetXaxis()->SetBinLabel(i, str_particle[i-1]);
    h_particle_type_numubar->GetXaxis()->SetBinLabel(i, str_particle[i-1]);
  }

  TCanvas * c4 = new TCanvas();
  c4->cd();
  h_particle_type->Draw();
  c4->Print("plots/h_particle_type.pdf");

  TCanvas * c5 = new TCanvas();
  c5->cd();
  h_particle_type_nue->Draw();
  c5->Print("plots/h_particle_type_nue.pdf");

  TCanvas * c6 = new TCanvas();
  c6->cd();
  h_particle_type_nuebar->Draw();
  c6->Print("plots/h_particle_type_nuebar.pdf");

  TCanvas * c7 = new TCanvas();
  c7->cd();
  h_particle_type_numu->Draw();
  c7->Print("plots/h_particle_type_numu.pdf");

  TCanvas * c8 = new TCanvas();
  c8->cd();
  h_particle_type_numubar->Draw();
  c8->Print("plots/h_particle_type_numubar.pdf");

  std::cout << "End Job!" << std::endl;
}

DEFINE_ART_MODULE(mcfluxfinder::MCFluxFinder)
