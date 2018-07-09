#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "uboone/EventWeight/MCEventWeight.h"

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

namespace beamweight {
  class BeamWeight;
}

class beamweight::BeamWeight : public art::EDAnalyzer {

public:
  explicit BeamWeight (fhicl::ParameterSet const & p);


  // Plugins should not be copied or assigned.
  BeamWeight(BeamWeight const &) = delete;
  BeamWeight(BeamWeight &&) = delete;
  BeamWeight & operator = (BeamWeight const &) = delete;
  BeamWeight & operator = (BeamWeight &&) = delete;

  //Required functions.
  void analyze(art::Event const & e) override;
  //void endSubRun(art::SubRun const & sr) override;

  void beginJob() override;
  void endJob() override;

private:

  std::string evtwght_tag = "eventweight";

  //systematic - universe 
  std::vector< std::vector< double > > Weights;

  std::vector< std::string > labels;

  TTree * prescale_tree;
  int fRun = 0;
  int fEvent = 0;
  int fSubRun = 0;

  TH1D * h_numi_prescale = new TH1D("numi_prescale", "numi_prescale", 3500, 0, 8000);

};//end class def

void beamweight::BeamWeight::beginJob()
{
  std::cout << "--- Begin Job --- " << std::endl;
}

beamweight::BeamWeight::BeamWeight(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
  art::ServiceHandle<art::TFileService> fs;
  prescale_tree = fs->make<TTree>("prescale_tree", "");

  prescale_tree->Branch("Run", &fRun, "fRun/I");
  prescale_tree->Branch("SubRun", &fSubRun, "fSubRun/I");
  prescale_tree->Branch("Event", &fEvent, "fEvent/I");
}

void beamweight::BeamWeight::analyze(art::Event const & e) {

  labels = {"FluxUnisim", "piplus", "piminus",
            "kplus", "kminus", "kzero","total"};
  Weights.resize(6);

  for(int i = 0; i < 6; i++)
  {
    Weights[i].resize(1000);
  }

  const int run = e.id().run();
  const int event = e.id().event();
  const int subRun = e.id().subRun();

  fRun = run;
  fSubRun = subRun;
  fEvent = event;

  art::Handle<std::vector<evwgh::MCEventWeight> > EventWeightHandle;
  e.getByLabel(evtwght_tag, EventWeightHandle);
  if(!EventWeightHandle.isValid()) {std::cout << "Event Weight Handle is not valid!" << std::endl; exit(1);}
  //auto const & EventWeightHandle = *e.getValidHandle<std::vector<evwgh::MCEventWeight> >(evtwght_tag);
  //const int event_weight_size = EventWeightHandle.size();

  for(auto const & event_weight : (*EventWeightHandle))
  //for(int i = 0; i < event_weight_size; i++)
  {
    //auto const event_weight = EventWeightHandle.at(i);
    /*
    for(auto last : evtwght.fWeight)
    {
      if(last.first.find("bnbcorrection") != std::string::npos)
      {
        bnbweight = last.second.at(0);
      }
    }
    */
    for(auto last : event_weight.fWeight)
    {
      for(int l = 0; l < int(labels.size())-1; l++)
      {
        if(last.first.find(labels[l].c_str()) != std::string::npos)
        {
          for(int i = 0; i < int(last.second.size()); i++)
          {         
              Weights[l][i] *= last.second.at(i);
              std::cout << Weights[l][i] << std::endl;
          }
        }
      }
    }
  }//end looping over size of event_weights (Handle)


  prescale_tree->Fill();

}//end analyzer loop

void beamweight::BeamWeight::endJob()
{
  //TCanvas * c1 = new TCanvas();
  //c1->cd();
  //h_numi_prescale->GetXaxis()->SetTitle("Run Number");
  //h_numi_prescale->GetYaxis()->SetTitle("Prescale Factor");
  //h_numi_prescale->Draw();
  //c1->Print("plots/numi_prescale.pdf");
}

DEFINE_ART_MODULE(beamweight::BeamWeight)
