#include "AnaHelper.h"
//#include "larcoreobj/POTSummary.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <iostream>

namespace potfinder {
  class POTFinder;
}

class potfinder::POTFinder : public art::EDAnalyzer {

public:
  explicit POTFinder (fhicl::ParameterSet const & p);


  // Plugins should not be copied or assigned.
  POTFinder(POTFinder const &) = delete;
  POTFinder(POTFinder &&) = delete;
  POTFinder & operator = (POTFinder const &) = delete;
  POTFinder & operator = (POTFinder &&) = delete;

  //Required functions.
  void analyze(art::Event const & e) override;
  void endSubRun(art::SubRun const & sr) override;

  void beginJob() override;

private:
  TTree * pot_tree;
  double pot = 0.0;

  int mc_numu_counter = 0;
  int mc_numu_bar_counter = 0;
  int mc_nue_counter = 0;
  int mc_nue_bar_counter = 0;
  int mc_numu_cc_qe_counter = 0;
  int mc_numu_bar_cc_qe_counter = 0;
  int mc_nue_cc_qe_counter = 0;
  int mc_nue_bar_cc_qe_counter = 0;

};//end class def

void potfinder::POTFinder::beginJob()
{
  std::cout << "--- Begin Job ---" << std::endl;
  // Implementation of optional member function here.
  art::ServiceHandle< art::TFileService > tfs;
  //define trees
  pot_tree = tfs->make<TTree>("pot_tree", "pot_per_subrun");

  pot_tree->Branch("pot", &pot, "pot/D");
}

potfinder::POTFinder::POTFinder(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{

}

void potfinder::POTFinder::analyze(art::Event const & e) {

  art::Handle < std::vector < simb::MCTruth > > MCTruthHandle;
  e.getByLabel("generator", MCTruthHandle);
  if(!MCTruthHandle.isValid()) {std::cout << "MCTruth Handle is not valid" << std::endl; exit(1); }
  for(auto const & mctruth : (*MCTruthHandle) )
  {
    auto const mc_nu = mctruth.GetNeutrino().Nu();
    if(mc_nu.PdgCode() ==  12)
    {
      mc_nue_counter++;
      if(mctruth.GetNeutrino().CCNC() == 0 && mctruth.GetNeutrino().Mode() == 0)
      //CCQE
      {mc_nue_cc_qe_counter++;}
    }
    if(mc_nu.PdgCode() == -12)
    {
      mc_nue_bar_counter++;
      if(mctruth.GetNeutrino().CCNC() == 0 && mctruth.GetNeutrino().Mode() == 0)
      //CCQE
      {mc_nue_bar_cc_qe_counter++;}
    }
    if(mc_nu.PdgCode() ==  14)
    {
      mc_numu_counter++;
      if(mctruth.GetNeutrino().CCNC() == 0 && mctruth.GetNeutrino().Mode() == 0)
      //CCQE
      {mc_numu_cc_qe_counter++;}
    }
    if(mc_nu.PdgCode() == -14)
    {
      mc_numu_bar_counter++;
      if(mctruth.GetNeutrino().CCNC() == 0 && mctruth.GetNeutrino().Mode() == 0)
      //CCQE
      {mc_numu_bar_cc_qe_counter++;}
    }
  }
}

void potfinder::POTFinder::endSubRun(art::SubRun const & sr) 
{ 
 
  auto const & POTSummaryHandle = sr.getValidHandle < sumdata::POTSummary >("generator");
  auto const & POTSummary(*POTSummaryHandle);
  const double total_pot = POTSummary.totpot;
  std::cout << "----------------------------" << std::endl;
  std::cout << "Total POT / subRun: " << total_pot << std::endl;
  std::cout << "----------------------------" << std::endl;

  pot = total_pot;
  pot_tree->Fill();

  std::cout << "MC Nue Counter           : " << mc_nue_counter << std::endl;
  std::cout << "MC Nue Bar Counter       : " << mc_nue_bar_counter << std::endl;
  std::cout << "MC Numu Counter          : " << mc_numu_counter << std::endl;
  std::cout << "MC Numu Bar Counter      : " << mc_numu_bar_counter << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << "MC Nue CCQE Counter      : " << mc_nue_cc_qe_counter << std::endl;
  std::cout << "MC Nue CCQE Bar Counter  : " << mc_nue_bar_cc_qe_counter << std::endl;
  std::cout << "MC Numu CCQE Counter     : " << mc_numu_cc_qe_counter << std::endl;
  std::cout << "MC Numu CCQE Bar Counter : " << mc_numu_bar_cc_qe_counter << std::endl;
  std::cout << '\n' << std::endl;

}//end endSubRun


DEFINE_ART_MODULE(potfinder::POTFinder)
