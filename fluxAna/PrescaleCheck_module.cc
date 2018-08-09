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

  double fPrescaleNuMI         = 0;
  double fPrescaleBNB          = 0;
  double fPrescaleUnbiased     = 0;
  double fPrescaleBeamNuMI     = 0;
  double fPrescaleNuMIUnbiased = 0;

  bool fPassedNuMI         = 0;
  bool fPassedBNB          = 0;
  bool fPassedUnbiased     = 0;
  bool fPassedBeamNuMI     = 0;
  bool fPassedNuMIUnbiased = 0;

  bool fPassedPrescaleNuMI         = 0;
  bool fPassedPrescaleBNB          = 0;
  bool fPassedPrescaleUnbiased     = 0;
  bool fPassedPrescaleBeamNuMI     = 0;
  bool fPassedPrescaleNuMIUnbiased = 0;

  int fPHMaxNuMI         = 0;
  int fPHMaxBNB          = 0;
  int fPHMaxUnbiased     = 0;
  int fPHMaxBeamNuMI     = 0;
  int fPHMaxNuMIUnbiased = 0;

  int fMultiplicityNuMI         = 0;
  int fMultiplicityBNB          = 0;
  int fMultiplicityUnbiased     = 0;
  int fMultiplicityBeamNuMI     = 0;
  int fMultiplicityNuMIUnbiased = 0;

  int fTriggerTickNuMI         = 0;
  int fTriggerTickBNB          = 0;
  int fTriggerTickUnbiased     = 0;
  int fTriggerTickBeamNuMI     = 0;
  int fTriggerTickNuMIUnbiased = 0;

  double fTimeSinceTriggerNuMI         = 0;
  double fTimeSinceTriggerBNB          = 0;
  double fTimeSinceTriggerUnbiased     = 0;
  double fTimeSinceTriggerBeamNuMI     = 0;
  double fTimeSinceTriggerNuMIUnbiased = 0;

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

  prescale_tree->Branch("Run",                    &fRun,                    "fRun/I");
  prescale_tree->Branch("SubRun",                 &fSubRun,                 "fSubRun/I");
  prescale_tree->Branch("Event",                  &fEvent,                  "fEvent/I");

  prescale_tree->Branch("PrescaleNuMI",           &fPrescaleNuMI,           "fPrescaleNuMI/D");
  prescale_tree->Branch("PrescaleBNB",            &fPrescaleBNB,            "fPrescaleBNB/D");
  prescale_tree->Branch("PrescaleUnbiased",       &fPrescaleUnbiased,       "fPrescaleUnbiased/D");
  prescale_tree->Branch("PrescaleBeamNuMI",       &fPrescaleBeamNuMI,       "fPrescaleBeamNuMI/D");
  prescale_tree->Branch("PrescaleNuMIUnbiased",   &fPrescaleNuMIUnbiased,   "fPrescaleNuMIUnbiased/D");

  prescale_tree->Branch("PassedNuMI",             &fPassedNuMI,             "fPassedNuMI/O");
  prescale_tree->Branch("PassedBNB",              &fPassedBNB,              "fPassedBNB/O");
  prescale_tree->Branch("PassedUnbiased",         &fPassedUnbiased,         "fPassedUnbiased/O");
  prescale_tree->Branch("PassedBeamNuMI",         &fPassedBeamNuMI,         "fPassedBeamNuMI/O");
  prescale_tree->Branch("PassedNuMIUnbiased",     &fPassedNuMIUnbiased,     "fPassedNuMIUnbiased/O");

  prescale_tree->Branch("PassedPrescaleNuMI",         &fPassedPrescaleNuMI,         "fPassedPrescaleNuMI/O");
  prescale_tree->Branch("PassedPrescaleBNB",          &fPassedPrescaleBNB,          "fPassedPrescaleBNB/O");
  prescale_tree->Branch("PassedPrescaleUnbiased",     &fPassedPrescaleUnbiased,     "fPassedPrescaleUnbiased/O");
  prescale_tree->Branch("PassedPrescaleBeamNuMI",     &fPassedPrescaleBeamNuMI,     "fPassedPrescaleBeamNuMI/O");
  prescale_tree->Branch("PassedPrescaleNuMIUnbiased", &fPassedPrescaleNuMIUnbiased, "fPassedPrescaleNuMIUnbiased/O");

  prescale_tree->Branch("PHMaxNuMI",                &fPHMaxNuMI,                "fPHMaxNuMI/I");
  prescale_tree->Branch("PHMaxBNB",                 &fPHMaxBNB,                 "fPHMaxBNB/I");
  prescale_tree->Branch("PHMaxUnbiased",            &fPHMaxUnbiased,            "fPHMaxUnbiased/I");
  prescale_tree->Branch("PHMaxBeamNuMI",            &fPHMaxBeamNuMI,            "fPHMaxBeamNuMI/I");
  prescale_tree->Branch("PHMaxNuMIUnbiased",        &fPHMaxNuMIUnbiased,        "fPHMaxNuMIUnbiased/I");

  prescale_tree->Branch("MultiplicityNuMI",         &fMultiplicityNuMI,         "fMultiplicityNuMI/I");
  prescale_tree->Branch("MultiplicityBNB",          &fMultiplicityBNB,          "fMultiplicityBNB/I");
  prescale_tree->Branch("MultiplicityUnbiased",     &fMultiplicityUnbiased,     "fMultiplicityUnbiased/I");
  prescale_tree->Branch("MultiplicityBeamNuMI",     &fMultiplicityBeamNuMI,     "fMultiplicityNuMIBeam/I");
  prescale_tree->Branch("MultplicityNuMIUnbiased",  &fMultiplicityNuMIUnbiased, "fMultiplicityNuMIUnbiased/I");

  prescale_tree->Branch("TriggerTickNuMI",          &fTriggerTickNuMI,          "fTriggerTickNuMI/I");
  prescale_tree->Branch("TriggerTickBNB",           &fTriggerTickBNB,           "fTriggerTickBNB/I");
  prescale_tree->Branch("TriggerTickUnbiased",      &fTriggerTickUnbiased,      "fTriggerTickUnbiased/I");
  prescale_tree->Branch("TriggerTickBeamNuMI",      &fTriggerTickBeamNuMI,      "fTriggerTickBeamNuMI/I");
  prescale_tree->Branch("TriggerTickNuMIUnbiased",  &fTriggerTickNuMIUnbiased,  "fTriggerTickNuMIUnbiased/I");

  prescale_tree->Branch("TimeSinceTriggerNuMI",         &fTimeSinceTriggerNuMI,         "fTimeSinceTriggerNuMI/D");
  prescale_tree->Branch("TimeSinceTriggerBNB",          &fTimeSinceTriggerBNB,          "fTimeSinceTriggerBNB/D");
  prescale_tree->Branch("TimeSinceTriggerUnbiased",     &fTimeSinceTriggerUnbiased,     "fTimeSinceTriggerUnbiased/D");
  prescale_tree->Branch("TimeSinceTriggerBeamNuMI",     &fTimeSinceTriggerBeamNuMI,     "fTimeSinceTriggerBeamNuMI/D");
  prescale_tree->Branch("TimeSinceTriggerNuMIUnbiased", &fTimeSinceTriggerNuMIUnbiased, "fTimeSinceTriggerNuMIUnbiased/D");

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
  //
  //NUMI_FEMBeamTriggerAlgo
  //NUMI_unbiased_PrescaleAlgo
  //
  const float prescale_bnb           = SWTriggerHandle->getPrescale("EXT_BNBwin_FEMBeamTriggerAlgo");
  const float prescale_numi          = SWTriggerHandle->getPrescale("EXT_NUMIwin_FEMBeamTriggerAlgo");
  const float prescale_unbiased      = SWTriggerHandle->getPrescale("EXT_unbiased_PrescaleAlgo");
  const float prescale_beam_numi     = SWTriggerHandle->getPrescale("NUMI_FEMBeamTriggerAlgo");
  const float prescale_numi_unbiased = SWTriggerHandle->getPrescale("NUMI_unbiased_PrescaleAlgo");

  const bool passed_bnb           = SWTriggerHandle->passedAlgo("EXT_BNBwin_FEMBeamTriggerAlgo");
  const bool passed_numi          = SWTriggerHandle->passedAlgo("EXT_NUMIwin_FEMBeamTriggerAlgo");
  const bool passed_unbiased      = SWTriggerHandle->passedAlgo("EXT_unbiased_PrescaleAlgo");
  const bool passed_beam_numi     = SWTriggerHandle->passedAlgo("NUMI_FEMBeamTriggerAlgo");
  const bool passed_numi_unbiased = SWTriggerHandle->passedAlgo("NUMI_unbiased_PrescaleAlgo");

  const bool passed_prescale_numi          = SWTriggerHandle->passedPrescaleAlgo("EXT_NUMIwin_FEMBeamTriggerAlgo");
  const bool passed_prescale_bnb           = SWTriggerHandle->passedPrescaleAlgo("EXT_BNBwin_FEMBeamTriggerAlgo");
  const bool passed_prescale_unbiased      = SWTriggerHandle->passedPrescaleAlgo("EXT_unbiased_PrescaleAlgo");
  const bool passed_prescale_beam_numi     = SWTriggerHandle->passedPrescaleAlgo("NUMI_FEMBeamTriggerAlgo");
  const bool passed_prescale_numi_unbiased = SWTriggerHandle->passedPrescaleAlgo("NUMI_unbiased_PrescaleAlgo");


  const int ph_max_numi          = SWTriggerHandle->getPhmax("EXT_NUMIwin_FEMBeamTriggerAlgo"); 
  const int ph_max_bnb           = SWTriggerHandle->getPhmax("EXT_BNBwin_FEMBeamTriggerAlgo"); 
  const int ph_max_unbiased      = SWTriggerHandle->getPhmax("EXT_unbiased_PrescaleAlgo");
  const int ph_max_beam_numi     = SWTriggerHandle->getPhmax("NUMI_FEMBeamTriggerAlgo");
  const int ph_max_numi_unbiased = SWTriggerHandle->getPhmax("NUMI_unbiased_PrescaleAlgo");


  const int mult_numi          = SWTriggerHandle->getMultiplicity("EXT_NUMIwin_FEMBeamTriggerAlgo");
  const int mult_bnb           = SWTriggerHandle->getMultiplicity("EXT_BNBwin_FEMBeamTriggerAlgo");
  const int mult_unbiased      = SWTriggerHandle->getMultiplicity("EXT_unbiased_PrescaleAlgo"); 
  const int mult_beam_numi     = SWTriggerHandle->getMultiplicity("NUMI_FEMBeamTriggerAlgo");
  const int mult_numi_unbiased = SWTriggerHandle->getMultiplicity("NUMI_unbiased_PrescaleAlgo");

  const int trigger_tick_numi          = SWTriggerHandle->getTriggerTick("EXT_NUMIwin_FEMBeamTriggerAlgo");
  const int trigger_tick_bnb           = SWTriggerHandle->getTriggerTick("EXT_BNBwin_FEMBeamTriggerAlgo");  
  const int trigger_tick_unbiased      = SWTriggerHandle->getTriggerTick("EXT_unbiased_PrescaleAlgo");
  const int trigger_tick_beam_numi     = SWTriggerHandle->getTriggerTick("NUMI_FEMBeamTriggerAlgo");
  const int trigger_tick_numi_unbiased = SWTriggerHandle->getTriggerTick("NUMI_unbiased_PrescaleAlgo");

  const double delta_trigger_numi          = SWTriggerHandle->getTimeSinceTrigger("EXT_NUMIwin_FEMBeamTriggerAlgo");  
  const double delta_trigger_bnb           = SWTriggerHandle->getTimeSinceTrigger("EXT_BNBwin_FEMBeamTriggerAlgo");
  const double delta_trigger_unbiased      = SWTriggerHandle->getTimeSinceTrigger("EXT_unbiased_PrescaleAlgo");
  const double delta_trigger_beam_numi     = SWTriggerHandle->getTimeSinceTrigger("NUMI_FEMBeamTriggerAlgo");
  const double delta_trigger_numi_unbiased = SWTriggerHandle->getTimeSinceTrigger("NUMI_unbiased_PrescaleAlgo");
  
  h_numi_prescale->Fill(run, 1./prescale_numi);

  fRun = run;
  fSubRun = subRun;
  fEvent = event;

  fPrescaleNuMI         = 1./prescale_numi;
  fPrescaleBNB          = 1./prescale_bnb;
  fPrescaleUnbiased     = 1./prescale_unbiased;
  fPrescaleBeamNuMI     = 1./prescale_beam_numi;
  fPrescaleNuMIUnbiased = 1./prescale_numi_unbiased;

  fPassedNuMI         = passed_numi;
  fPassedBNB          = passed_bnb;
  fPassedUnbiased     = passed_unbiased;
  fPassedBeamNuMI     = passed_beam_numi;
  fPassedNuMIUnbiased = passed_numi_unbiased;

  fPassedPrescaleNuMI         = passed_prescale_numi;
  fPassedPrescaleBNB          = passed_prescale_bnb;
  fPassedPrescaleUnbiased     = passed_prescale_unbiased;
  fPassedPrescaleBeamNuMI     = passed_prescale_beam_numi;
  fPassedPrescaleNuMIUnbiased = passed_prescale_numi_unbiased;

  fPHMaxNuMI         = ph_max_numi;
  fPHMaxBNB          = ph_max_bnb;
  fPHMaxUnbiased     = ph_max_unbiased;
  fPHMaxBeamNuMI     = ph_max_beam_numi;
  fPHMaxNuMIUnbiased = ph_max_numi_unbiased;

  fMultiplicityNuMI         = mult_numi;
  fMultiplicityBNB          = mult_bnb;
  fMultiplicityUnbiased     = mult_unbiased;
  fMultiplicityBeamNuMI     = mult_beam_numi;
  fMultiplicityNuMIUnbiased = mult_numi_unbiased;

  fTriggerTickNuMI         = trigger_tick_numi;
  fTriggerTickBNB          = trigger_tick_bnb;
  fTriggerTickUnbiased     = trigger_tick_unbiased;
  fTriggerTickBeamNuMI     = trigger_tick_beam_numi;
  fTriggerTickNuMIUnbiased = trigger_tick_numi_unbiased;

  fTimeSinceTriggerNuMI         = delta_trigger_numi;
  fTimeSinceTriggerBNB          = delta_trigger_bnb;
  fTimeSinceTriggerUnbiased     = delta_trigger_unbiased;
  fTimeSinceTriggerBeamNuMI     = delta_trigger_beam_numi;
  fTimeSinceTriggerNuMIUnbiased = delta_trigger_numi_unbiased;

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
