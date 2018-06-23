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
  for(auto const & mctruth : (*MCTruthHandle) )
  {
    auto const mc_nu = mctruth.GetNeutrino().Nu();
    if(mc_nu.PdgCode() ==  12)
    {
      if(mctruth.GetNeutrino().CCNC() == 0 && mctruth.GetNeutrino().Mode() == 0)
      {
        std::cout << "CCQE event!" << std::endl;
      }
      //CCQE
    }
  }

  art::Handle < std::vector < simb::MCFlux > > MCFluxHandle;
  e.getByLabel("generator", MCFluxHandle);
  if(!MCFluxHandle.isValid()) {std::cout << "MCFlux Handle is not valid!" << std::endl; exit(1); }
  for(auto const & mcflux : (*MCFluxHandle) )
  {
   std::cout << mcflux.ftgptype << std::endl; 
  }

}//end analyzer loop

void mcfluxfinder::MCFluxFinder::endJob()
{
  std::cout << "End Job!" << std::endl;
}

DEFINE_ART_MODULE(mcfluxfinder::MCFluxFinder)
