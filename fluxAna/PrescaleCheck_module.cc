#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

//#include "/uboone/app/users/chill2/uboonecode_v06_26_01_13/srcs/uboonecode/uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "nusimdata/SimulationBase/MCFlux.h"

#include <fstream>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

namespace prescalecheck {
  class PrescaleCheck;
}

class prescalecheck::PrescaleCheck : public art::EDAnalyzer {

public:
  explicit PrescaleCheck (fhicl::ParameterSet const & p);


  // Plugins should not be copied or assigned.
  PrescaleCheck(PrescaleCheck const &) = delete;
  PrescaleCheck(PrescaleCheck &&) = delete;
  PrescaleCheck & operator = (PrescaleCheck const &) = delete;
  PrescaleCheck & operator = (PrescaleCheck &&) = delete;

  //Required functions.
  void analyze(art::Event const & e) override;
  //void endSubRun(art::SubRun const & sr) override;

  void beginJob() override;
  void endJob() override;

private:

  TTree * prescale_tree;
  int fRun = 0;
  int fEvent = 0;
  int fSubRun = 0;
  double fPrescaleNuMI = 0;

  TH1D * h_numi_prescale = new TH1D("numi_prescale", "numi_prescale", 3500, 0, 8000);

};//end class def

void prescalecheck::PrescaleCheck::beginJob()
{
  std::cout << "--- Begin Job --- " << std::endl;
}

prescalecheck::PrescaleCheck::PrescaleCheck(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
  art::ServiceHandle<art::TFileService> fs;
  prescale_tree = fs->make<TTree>("prescale_tree", "");

  prescale_tree->Branch("Run", &fRun, "fRun/I");
  prescale_tree->Branch("SubRun", &fSubRun, "fSubRun/I");
  prescale_tree->Branch("Event", &fEvent, "fEvent/I");
  prescale_tree->Branch("PrescaleNuMI", &fPrescaleNuMI, "fPrescaleNuMI/D");
}

void prescalecheck::PrescaleCheck::analyze(art::Event const & e) {

  const int run = e.id().run();
  const int event = e.id().event();
  const int subRun = e.id().subRun();

  art::Handle<raw::ubdaqSoftwareTriggerData> SWTriggerHandle;
  e.getByLabel("daq", SWTriggerHandle);
  if(!SWTriggerHandle.isValid()) {std::cout << "SW Trigger Handle is not valid!" << std::endl; exit(1);}

  std::cout << "Number of algorithms = " << SWTriggerHandle->getNumberOfAlgorithms() << std::endl;	
  std::vector<std::string> algoNames = SWTriggerHandle->getListOfAlgorithms();	
  //for(unsigned int i = 0; i < algoNames.size(); i++)
  //{
  //  std::cout << algoNames.at(i) << std::endl;
  //}
  const float prescale_bnb = SWTriggerHandle->getPrescale("EXT_BNBwin_FEMBeamTriggerAlgo");
  const float prescale_numi = SWTriggerHandle->getPrescale("EXT_NUMIwin_FEMBeamTriggerAlgo");
  const float prescale_unbiased = SWTriggerHandle->getPrescale("EXT_unbiased_PrescaleAlgo");

  std::cout << "-------------" << std::endl;
  std::cout << "EXT BNB  : " << prescale_bnb << std::endl;
  std::cout << "EXT NuMI : " << prescale_numi << std::endl;
  std::cout << "EXT Unbiased Prescale: " << prescale_unbiased << std::endl;
  std::cout << "--------------" << std::endl;

  
  h_numi_prescale->Fill(run, 1./prescale_numi);

  fRun = run;
  fSubRun = subRun;
  fEvent = event;
  fPrescaleNuMI = 1./prescale_numi;

  prescale_tree->Fill();

}//end analyzer loop

void prescalecheck::PrescaleCheck::endJob()
{
  //TCanvas * c1 = new TCanvas();
  //c1->cd();
  //h_numi_prescale->GetXaxis()->SetTitle("Run Number");
  //h_numi_prescale->GetYaxis()->SetTitle("Prescale Factor");
  //h_numi_prescale->Draw();
  //c1->Print("plots/numi_prescale.pdf");
}

DEFINE_ART_MODULE(prescalecheck::PrescaleCheck)
