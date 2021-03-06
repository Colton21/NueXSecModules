////////////////////////////////////////////////////////////////////////
// Class:       RecoTrueMatching
// Plugin Type: producer (art v2_05_00)
// File:        RecoTrueMatching_module.cc
//
// Generated at Fri Aug 18 08:43:34 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

// Data product include
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTracker.h"
#include "MCGhost.h"

#include <memory>

// Algorithms include
#include "RecoTrueHelper.h"

namespace xsecAna {
class RecoTrueMatching;
}

class xsecAna::RecoTrueMatching : public art::EDProducer {
public:
explicit RecoTrueMatching(fhicl::ParameterSet const & p);
// The compiler-generated destructor is fine for non-base
// classes without bare pointers or other resource use.

// Plugins should not be copied or assigned.
RecoTrueMatching(RecoTrueMatching const &) = delete;
RecoTrueMatching(RecoTrueMatching &&) = delete;
RecoTrueMatching & operator = (RecoTrueMatching const &) = delete;
RecoTrueMatching & operator = (RecoTrueMatching &&) = delete;

// Required functions.
void produce(art::Event & e) override;

private:

std::string _pfp_producer;
std::string _spacepointLabel;
std::string _hitfinderLabel;
std::string _geantModuleLabel;

bool _is_data;
bool _debug;
bool _cosmic_only;
bool _use_premade_ass;
std::string _mcpHitAssLabel;
};


xsecAna::RecoTrueMatching::RecoTrueMatching(fhicl::ParameterSet const & p) {

	std::cout << "[RecoTrueMatching] Setting fcl Parameters" << std::endl;

	_pfp_producer                   = p.get<std::string>("PFParticleProducer");
	_hitfinderLabel                 = p.get<std::string>("HitProducer");
	_geantModuleLabel               = p.get<std::string>("GeantModule");
	_spacepointLabel                = p.get<std::string>("SpacePointProducer");

	_debug                          = p.get<bool>("Debug", true);
	_cosmic_only                    = p.get<bool>("CosmicOnly", false);
	_use_premade_ass                = p.get<bool>("UsePremadeAssociation", true);
	_mcpHitAssLabel                 = p.get<std::string>("MCPHitAssProducer", "pandoraCosmicHitRemoval");

	produces< std::vector<xsecAna::MCGhost> >();
	produces< art::Assns<simb::MCParticle, xsecAna::MCGhost> >();
	produces< art::Assns<recob::PFParticle, xsecAna::MCGhost> >();

	std::cout << "[RecoTrueMatching] End Setting fcl Parameters" << std::endl;

}

void xsecAna::RecoTrueMatching::produce(art::Event & e)
{
	nue_xsec::recotruehelper _recotruehelper_instance;

	if(_debug) std::cout << "[RecoTrueMatching] Starts" << std::endl;
	if(_debug) std::cout << "[RecoTrueMatching] event: " << e.id().event() << std::endl;


	// Instantiate the output
	std::unique_ptr< std::vector< xsecAna::MCGhost > >                 mcGhostVector   (new std::vector<xsecAna::MCGhost>);
	std::unique_ptr< art::Assns<simb::MCParticle, xsecAna::MCGhost> >  assnOutGhostMCP (new art::Assns<simb::MCParticle, xsecAna::MCGhost>);
	std::unique_ptr< art::Assns<recob::PFParticle, xsecAna::MCGhost> > assnOutGhostPFP (new art::Assns<recob::PFParticle, xsecAna::MCGhost>);


	if(_cosmic_only)
	{
		std::cout << "[RecoTrueMatching] Cosmic Only Configuration! - End Module" << std::endl;
		e.put(std::move(mcGhostVector));
		e.put(std::move(assnOutGhostMCP));
		e.put(std::move(assnOutGhostPFP));
		return;
	}

	_is_data = e.isRealData();

	if (_is_data) {
		std::cout << "[RecoTrueMatching] Running on a real data file. No MC-PFP matching will be attempted." << std::endl;
		e.put(std::move(mcGhostVector));
		e.put(std::move(assnOutGhostMCP));
		e.put(std::move(assnOutGhostPFP));
		return;
	}

	if(!_use_premade_ass) 
	{	
		std::cout << "[RecoTrueMatching] Constructing Associations w/ SimChannels for MCParticle<-->Hits " << std::endl;
		std::cout << "[RecoTrueMatching] Configuring associations " << std::endl;
		_recotruehelper_instance.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel); 
	}
	if(_use_premade_ass)  
	{
		std::cout << "[RecoTrueMatching] Using Premade Associations for MCParticle<-->Hits " << std::endl;
		std::cout << "[RecoTrueMatching] Configuring associations " << std::endl;
		_recotruehelper_instance.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel,
		                                                  _mcpHitAssLabel, lar_pandora::LArPandoraHelper::kAddDaughters); 
	}

	lar_pandora::MCParticlesToPFParticles matchedMCToPFParticles; // This is a map: MCParticle to matched PFParticle
	lar_pandora::MCParticlesToHits matchedParticleHits;

	_recotruehelper_instance.GetRecoToTrueMatches(matchedMCToPFParticles, matchedParticleHits);

	std::cout << "[RecoTrueMatching] Generating " << matchedMCToPFParticles.size() << " matched MC to Particle associations." << std::endl;

	for (auto const& iter : matchedMCToPFParticles) {

		art::Ptr<simb::MCParticle>  mc_par = iter.first;// The MCParticle
		art::Ptr<recob::PFParticle> pf_par = iter.second; // The matched PFParticle

		xsecAna::MCGhost mcGhost;
		mcGhost.SetMode("depEnergy");

		mcGhostVector->emplace_back(mcGhost);
		util::CreateAssn(*this, e, *mcGhostVector, pf_par, *assnOutGhostPFP);
		util::CreateAssn(*this, e, *mcGhostVector, mc_par, *assnOutGhostMCP);
	}

	e.put(std::move(mcGhostVector));
	e.put(std::move(assnOutGhostMCP));
	e.put(std::move(assnOutGhostPFP));
}

DEFINE_ART_MODULE(xsecAna::RecoTrueMatching)
