#include "AnaHelper.h"
//#include "larcoreobj/POTSummary.h"
#include "larcoreobj/SummaryData/POTSummary.h"

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

}//end endSubRun


DEFINE_ART_MODULE(potfinder::POTFinder)
