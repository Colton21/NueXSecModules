////////////////////////////////////////////////////////////////////////
// Class:       NueXsecExtractor
// Plugin Type: analyzer (art v2_05_00)
// File:        NueXsecExtractor_module.cc
//
// Generated at Fri Jun  2 10:54:39 2017 by Colton Hill using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TTree.h"
#include "TFile.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


class FlashSaver;


class FlashSaver : public art::EDAnalyzer {
public:
explicit FlashSaver(fhicl::ParameterSet const & p);
// The compiler-generated destructor is fine for non-base
// classes without bare pointers or other resource use.

// Plugins should not be copied or assigned.
FlashSaver(FlashSaver const &) = delete;
FlashSaver(FlashSaver &&) = delete;
FlashSaver & operator = (FlashSaver const &) = delete;
FlashSaver & operator = (FlashSaver &&) = delete;

// Required functions.
void analyze(art::Event const & e) override;

// Selected optional functions.
void beginJob() override;
void endJob() override;
void reconfigure(fhicl::ParameterSet const & p) override;

private:

//
// Declare member data here.
//
bool _verbose = true;

TTree * event_tree;
TTree * optical_tree;

int fEvent_num = 0;
int fRun_num = 0;
int fSubrun_num = 0;

int fOpFlashPE = 0;
double fOpFlashTime = 0;
double fOpFlashWidthY = 0;
double fOpFlashWidthZ = 0;
double fOpFlashCenterY = 0;
double fOpFlashCenterZ = 0;
};


FlashSaver::FlashSaver(fhicl::ParameterSet const & p)
	:
	EDAnalyzer(p) // ,
	// More initializers here.
{
}

void FlashSaver::reconfigure(fhicl::ParameterSet const & p)
{
	// Implementation of optional member function here.
}

void FlashSaver::analyze(art::Event const & evt)
{
	fEvent_num = evt.id().event();
	fRun_num = evt.run();
	fSubrun_num = evt.subRun();

	std::string beam_flash_tag = "simpleFlashBeam";
	auto const & beam_opf = evt.getValidHandle<std::vector < recob::OpFlash> >(beam_flash_tag);
	auto const & beam_opflashes(*beam_opf);

	//double pe_threshold = 50;
	//double min_time = 5;
	//double max_time = 16;

	//simple flash beam
	for(auto const & opflsh : beam_opflashes)
	{
		//if(opflsh.Time() >= min_time && opflsh.Time() <= max_time)
		//{
			//if(opflsh.TotalPE() >= pe_threshold)
			//{
				fOpFlashPE = opflsh.TotalPE();
				fOpFlashTime = opflsh.Time();
				fOpFlashWidthY = opflsh.YWidth();
				fOpFlashWidthZ = opflsh.ZWidth();
				fOpFlashCenterY = opflsh.YCenter();
				fOpFlashCenterZ = opflsh.ZCenter();
				optical_tree->Fill();
			//}        //end if pe threshold
		//}        //end if in time
	}        //end loop flashes
	event_tree->Fill();
}

void FlashSaver::beginJob()
{
	std::cout << "--- Begin Job ---" << std::endl;
	// Implementation of optional member function here.
	art::ServiceHandle< art::TFileService > tfs;

	//define trees
	event_tree = tfs->make<TTree>("nue_tree", "data_products");
	optical_tree = tfs->make<TTree>("optical_tree", "optical_objects");

	//event tree
	event_tree->Branch("event", &fEvent_num, "fEvent_num/I");
	event_tree->Branch("Run", &fRun_num, "fRun_num/I");
	event_tree->Branch("Subrun", &fSubrun_num, "fSubrun_num/I");

        //optical tree
	optical_tree->Branch("event", &fEvent_num, "fEvent_num/I");
	optical_tree->Branch("Run", &fRun_num, "fRun_num/I");
	optical_tree->Branch("Subrun", &fSubrun_num, "fSubrun_num/I");
	optical_tree->Branch("OpFlashPE", &fOpFlashPE, "fOpFlashPe/I");
	optical_tree->Branch("OpFlashTime", &fOpFlashTime, "fOpFlashTime/D");
	optical_tree->Branch("OpFlashWidhtY", &fOpFlashWidthY, "fOpFlashWidthY/D");
	optical_tree->Branch("OpFlashWidthZ", &fOpFlashWidthZ, "fOpFlashWidhtZ/D");
	optical_tree->Branch("OpFlashCenterY", &fOpFlashCenterY, "fOpFlashCenterY/D");
	optical_tree->Branch("OpFlashCenterZ", &fOpFlashCenterZ, "fOpFlashCenterZ/D");

}

void FlashSaver::endJob()
{
  std::cout << "End Job" << std::endl;

}


DEFINE_ART_MODULE(FlashSaver)
