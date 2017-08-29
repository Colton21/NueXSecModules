#include "XSec.h"
#include "RecoTrueHelper.h"

namespace nue_xsec {

void XSec::reconfigure(fhicl::ParameterSet const &p)
{

	_pfp_producer                   = p.get<std::string>("PFParticleProducer");
	_hitfinderLabel                 = p.get<std::string>("HitProducer");
	_geantModuleLabel               = p.get<std::string>("GeantModule");
	_spacepointLabel                = p.get<std::string>("SpacePointProducer");
	_neutrino_flash_match_producer  = p.get<std::string>("NeutrinoFlashMatchProducer");
	_cosmic_flash_match_producer    = p.get<std::string>("CosmicFlashMatchProducer");
	_opflash_producer_beam          = p.get<std::string>("OpFlashBeamProducer");
	_acpt_producer                  = p.get<std::string>("ACPTProducer");
	_tpcobject_producer             = p.get<std::string>("TPCObjectProducer");
	_potsum_producer                = p.get<std::string>("POTSummaryProducer");
	_potsum_instance                = p.get<std::string>("POTSummaryInstance");
	_particle_id_producer           = p.get<std::string>("ParticleIDProducer");
	_mc_ghost_producer              = p.get<std::string>("MCGhostProducer");

	_useDaughterMCParticles         = p.get<std::string>("UseDaughterMCParticles");
	_useDaughterPFParticles         = p.get<std::string>("UseDaughterPFParticles");
	_addDaughterMCParticles         = p.get<std::string>("AddDaughterMCParticles");
	_addDaughterPFParticles         = p.get<std::string>("AddDaughterPFParticles");

	_use_genie_info                 = p.get<bool>("UseGENIEInfo", false);
	_minimumHitRequirement          = p.get<int>("MinimumHitRequirement", 3);

	_beam_spill_start               = p.get<double>("BeamSpillStart", 3.2);
	_beam_spill_end                 = p.get<double>("BeamSpillEnd",   4.8);

	_debug                          = p.get<std::string>("Debug", false);
	_verbose                        = p.get<std::string>("Verbose", false);

}


XSec::XSec: EDAnalyzer(p) {


	myTree->Branch("run", &run, "run/I");
	myTree->Branch("event", &event, "event/I");
	myTree->Branch("index", &index, "index/I");
	myTree->Branch("nMCParticles", &nMCParticles, "nMCParticles/I");
	myTree->Branch("nMCNeutrinos", &nMCNeutrinos, "nMCNeutrinos/I");
	myTree->Branch("nPFPNeutrinos", &nPFPNeutrinos, "nPFPNeutrinos/I");
	myTree->Branch("mcPdg", &mcPdg, "mcPdg/I");
	myTree->Branch("mcNuPdg", &mcNuPdg, "mcNuPdg/I");
	myTree->Branch("mcNuIndex", &mcNuIndex, "mcNuIndex/I");
	myTree->Branch("mcParentPdg", &mcParentPdg, "mcParentPdg/I");
	myTree->Branch("mcIsNeutirno", &mcIsNeutirno, "mcIsNeutirno/O");
	myTree->Branch("mcIsPrimary", &mcIsPrimary, "mcIsPrimary/O");
	myTree->Branch("mcMode", &mcMode, "mcMode/I");
	myTree->Branch("mcIsCC", &mcIsCC, "mcIsCC/O");
	myTree->Branch("pfpPdg", &pfpPdg, "pfpPdg/I");
	myTree->Branch("pfpNuPdg", &pfpNuPdg, "pfpNuPdg/I");
	myTree->Branch("pfpNuIndex", &pfpNuIndex, "pfpNuIndex/I");
	myTree->Branch("pfpParentPdg", &pfpParentPdg, "pfpParentPdg/I");
	myTree->Branch("pfpIsNeutrino", &pfpIsNeutrino, "pfpIsNeutrino/O");

	myTree->Branch("mcVtxX", &mcVtxX, "mcVtxX/D");
	myTree->Branch("mcVtxY", &mcVtxY, "mcVtxY/D");
	myTree->Branch("mcVtxZ", &mcVtxZ, "mcVtxZ/D");
	myTree->Branch("pfpVtxX", &pfpVtxX, "pfpVtxX/D");
	myTree->Branch("pfpVtxY", &pfpVtxY, "pfpVtxY/D");
	myTree->Branch("pfpVtxZ", &pfpVtxZ, "pfpVtxZ/D");

	myTree->Branch("mcDirX", &mcDirX, "mcDirX/D");
	myTree->Branch("mcDirY", &mcDirY, "mcDirY/D");
	myTree->Branch("mcDirZ", &mcDirZ, "mcDirZ/D");
	myTree->Branch("pfpDirX", &pfpDirX, "pfpDirX/D");
	myTree->Branch("pfpDirY", &pfpDirY, "pfpDirY/D");
	myTree->Branch("pfpDirZ", &pfpDirZ, "pfpDirZ/D");

	myTree->Branch("mcTheta", &mcTheta, "mcTheta/D");
	myTree->Branch("mcPhi", &mcPhi, "mcPhi/D");
	myTree->Branch("pfpTheta", &pfpTheta, "pfpTheta/D");
	myTree->Branch("pfpPhi", &pfpPhi, "pfpPhi/D");

	myTree->Branch("mcLength", &mcLength, "mcLength/D");
	myTree->Branch("pfpLength", &pfpLength, "pfpLength/D");

	myTree->Branch("mcEnergy", &mcEnergy, "mcEnergy/D");
	myTree->Branch("mcMomentum", &mcMomentum, "mcMomentum/D");
	myTree->Branch("pfpMomentum", &pfpMomentum, "pfpMomentum/D");

	myTree->Branch("completeness", &completeness, "completeness/D");
	myTree->Branch("purity", &purity, "purity/D");

	myTree->Branch("nMCHits",   &nMCHits, "mcHits/I");
	myTree->Branch("nMCHitsU",  &nMCHitsU, "mcHitsU/I");
	myTree->Branch("nMCHitsV",  &nMCHitsV, "mcHitsV/I");
	myTree->Branch("nMCHitsY",  &nMCHitsY, "mcHitsY/I");
	myTree->Branch("nPFPHits",  &nPFPHits, "pfpHits/I");
	myTree->Branch("nPFPHitsU", &nPFPHitsU, "pfpHitsU/I");
	myTree->Branch("nPFPHitsV", &nPFPHitsV, "pfpHitsV/I");
	myTree->Branch("nPFPHitsY", &nPFPHitsY, "pfpHitsY/I");

	myTree->Branch("mcOpenAngle", &mcOpenAngle, "mcOpenAngle/D");
	myTree->Branch("pfpOpenAngle", &pfpOpenAngle, "pfpOpenAngle");

}

void XSec::analyze(art::Event const & e) {

	//First thing is to zero all of the values
	run = -9999;
	event = -9999;
	index = -9999;
	nMCParticles = -9999;
	nMCNeutrinos = -9999;
	nPFPartilcles = -9999;
	nPFPNeutrinos = -9999;

	mcPdg = -9999;
	mcNuPdg = -9999;
	mcNuIndex = -9999;
	mcParentPdg = -9999;
	mcIsNeutirno = -9999;
	mcIsPrimary = -9999;
	mcMode = -9999;
	mcIsCC = -9999;
	pfpPdg = -9999;
	pfpNuPdg = -9999;
	pfpNuIndex = -9999;
	pfpParentPdg = -9999;
	pfpIsNeutrino = -9999;

	mcVtxX = -9999;
	mcVtxY = -9999;
	mcVtxZ = -9999;
	pfpVtxX = -9999;
	pfpVtxY = -9999;
	pfpVtxZ = -9999;

	mcDirX = -9999;
	mcDirY = -9999;
	mcDirZ = -9999;
	pfpDirX = -9999;
	pfpDirY = -9999;
	pfpDirZ = -9999;

	mcTheta = -9999;
	mcPhi = -9999;
	pfpTheta = -9999;
	pfpPhi = -9999;

	mcLength = -9999;
	pfpLength = -9999;

	mcEnergy = -9999;
	mcMomentum = -9999;
	//double pfpEnergy
	pfpMomentum = -9999;

	completeness = -9999;
	purity = -9999;

	nMCHits = -9999;
	nMCHitsU = -9999;
	nMCHitsV = -9999;
	nMCHitsY = -9999;
	nPFPHits = -9999;
	nPFPHitsU = -9999;
	nPFPHitsV = -9999;
	nPFPHitsY = -9999;
	nMatchedHits = -9999;
	nMatchedHitsU = -9999;
	nMatchedHitsV = -9999;
	nMatchedHitsY = -9999;

	mcOpenAngle = -9999;
	pfpOpenAngle = -9999;


	run = e.id().run();
	//int _subrun = e.id().subRun();
	event = e.id().event();
	//bool _is_data = e.isRealData();
	//bool _is_mc = !_is_data;

	//I want to move this all to a separate file, but I should test that this builds as is first!
	//performing reco-true matching
	//====================================================
	//lar_pandora::MCParticlesToPFParticles matchedMCToPFParticles; // This is a map: MCParticle to matched PFParticle
	//lar_pandora::MCParticlesToHits matchedParticleHits;

	// Collect Tracks and PFParticle <-> Track Associations
	// ====================================================
	lar_pandora::TrackVector recoTrackVector;
	lar_pandora::PFParticlesToTracks recoParticlesToTracks;
	LArPandoraHelper::CollectTracks(e, _pfp_producer, recoTrackVector, recoParticlesToTracks);

	lar_pandora::T0Vector t0Vector_trk;
	lar_pandora::TracksToT0s tracksToT0s;
	LArPandoraHelper::CollectT0s(e, _pfp_producer, t0Vector_trk, tracksToT0s);

	if (_verbose)
		std::cout << "  Tracks: " << recoTrackVector.size() << std::endl;

	// Collect Showers and PFParticle <-> Shower Associations
	// ====================================================
	lar_pandora::ShowerVector recoShowerVector;
	lar_pandora::PFParticlesToShowers recoParticlesToShowers;
	LArPandoraHelper::CollectShowers(e, _pfp_producer, recoShowerVector, recoParticlesToShowers);

	lar_pandora::T0Vector t0Vector_shwr;
	lar_pandora::ShowersToT0s showersToT0s;
	LArPandoraHelper::CollectT0s(e, _pfp_producer, t0Vector_shwr, showersToT0s);

	if (_verbose)
		std::cout << "  Showers: " << recoShowerVector.size() << std::endl;

	// Collect Vertices and PFParticle <-> Vertex Associations
	// =======================================================
	lar_pandora::VertexVector recoVertexVector;
	lar_pandora::PFParticlesToVertices recoParticlesToVertices;
	LArPandoraHelper::CollectVertices(e,  _pfp_producer, recoVertexVector, recoParticlesToVertices);

	if (_verbose)
		std::cout << "  Vertices: " << recoVertexVector.size() << std::endl;

	// Collect PFParticles and match Reco Particles to Hits
	// ====================================================
	lar_pandora::PFParticleVector recoParticleVector;
	lar_pandora::PFParticleVector recoNeutrinoVector;
	lar_pandora::PFParticlesToHits recoParticlesToHits;
	lar_pandora::HitsToPFParticles recoHitsToParticles;

	LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
	LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
	LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, recoParticlesToHits, recoHitsToParticles,
	                                         (_useDaughterPFParticles ?
	                                          (_addDaughterPFParticles ? LArPandoraHelper::kAddDaughters : LArPandoraHelper::kUseDaughters) : LArPandoraHelper::kIgnoreDaughters));

	if (_verbose)
		std::cout << "  RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;

	if (_verbose)
		std::cout << "  RecoParticles: " << recoParticleVector.size() << std::endl;

	// Collect Hits
	// ============
	lar_pandora::HitVector hitVector;
	LArPandoraHelper::CollectHits(evt, m_hitfinderLabel, hitVector);

	if (_verbose)
		std::cout << "  Hits: " << hitVector.size() << std::endl;

	// Collect SpacePoints and SpacePoint <-> Hit Associations
	// =======================================================
	lar_pandora::SpacePointVector spacePointVector;
	lar_pandora::SpacePointsToHits spacePointsToHits;
	lar_pandora::HitsToSpacePoints hitsToSpacePoints;
	LArPandoraHelper::CollectSpacePoints(evt, m_particleLabel, spacePointVector, spacePointsToHits, hitsToSpacePoints);

	if (_verbose)
		std::cout << "  SpacePoints: " << spacePointVector.size() << std::endl;

	// Collect MCParticles and match True Particles to Hits
	// ====================================================
	lar_pandora::MCParticleVector trueParticleVector;
	lar_pandora::MCTruthToMCParticles truthToParticles;
	lar_pandora::MCParticlesToMCTruth particlesToTruth;
	lar_pandora::MCParticlesToHits trueParticlesToHits;
	lar_pandora::HitsToMCParticles trueHitsToParticles;

	if (!e.isRealData())
	{
		LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, trueParticleVector);
		LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, truthToParticles, particlesToTruth);
		LArPandoraHelper::BuildMCParticleHitMaps(e, _geantModuleLabel, hitVector, trueParticlesToHits, trueHitsToParticles,
		                                         (_useDaughterMCParticles ?
		                                          (_addDaughterMCParticles ? LArPandoraHelper::kAddDaughters : LArPandoraHelper::kUseDaughters) : LArPandoraHelper::kIgnoreDaughters));
	}

	if (_verbose)
		std::cout << "  TrueParticles: " << particlesToTruth.size() << std::endl;

	if (_verbose)
		std::cout << "  TrueEvents: " << truthToParticles.size() << std::endl;

	if (trueParticlesToHits.empty())
	{
		myTree->Fill();
		return;
	}

	// Build Reco and True Particle Maps (for Parent/Daughter Navigation)
	// =================================================================
	lar_pandora::MCParticleMap trueParticleMap;
	lar_pandora::PFParticleMap recoParticleMap;

	myMatcher.BuildTrueParticleMap(trueParticleVector, trueParticleMap);
	myMatcher.BuildRecoParticleMap(recoParticleVector, recoParticleMap);

	//nMCParticles  = trueParticlesToHits.size();
	//nNeutrinoPfos = 0;
	//nPrimaryPfos  = 0;
	//nDaughterPfos = 0;

	// Count reconstructed particles
	// for (PFParticleVector::const_iterator iter = recoParticleVector.begin(), iterEnd = recoParticleVector.end(); iter != iterEnd; ++iter)
	// {
	//      const art::Ptr<recob::PFParticle> recoParticle = *iter;
	//
	//      if (LArPandoraHelper::IsNeutrino(recoParticle))
	//      {
	//              m_nNeutrinoPfos++;
	//      }
	//      else if (LArPandoraHelper::IsFinalState(recoParticleMap, recoParticle))
	//      {
	//              m_nPrimaryPfos++;
	//      }
	//      else
	//      {
	//              m_nDaughterPfos++;
	//      }
	// }

	// Match Reco Neutrinos to True Neutrinos
	// ======================================
	nue_xsec::recotruehelper myMatcher;
	lar_pandora::PFParticlesToHits recoNeutrinosToHits;
	lar_pandora::HitsToPFParticles recoHitsToNeutrinos;
	lar_pandora::HitsToMCTruth trueHitsToNeutrinos;
	lar_pandora::MCTruthToHits trueNeutrinosToHits;
	myMatcher.BuildRecoNeutrinoHitMaps(recoParticleMap, recoParticlesToHits, recoNeutrinosToHits, recoHitsToNeutrinos);
	myMatcher.BuildTrueNeutrinoHitMaps(truthToParticles, trueParticlesToHits, trueNeutrinosToHits, trueHitsToNeutrinos);

	lar_pandora::MCTruthToPFParticles matchedNeutrinos;
	lar_pandora::MCTruthToHits matchedNeutrinoHits;
	myMatcher.GetRecoToTrueMatches(recoNeutrinosToHits, trueHitsToNeutrinos, matchedNeutrinos, matchedNeutrinoHits);


	//Start to Fill the MC Truth information - First we do neutrinos
	//=================================================
	for (lar_pandora::MCTruthToHits::const_iterator iter = trueNeutrinosToHits.begin(), iterEnd = trueNeutrinosToHits.end(); iter != iterEnd; ++iter)
	{
		const art::Ptr<simb::MCTruth> trueEvent = iter->first;
		const lar_pandora::HitVector &trueHitVector = iter->second;

		if (trueHitVector.empty())
			continue;

		//**** WHAT DOES THIS DO??? ****
		if (!trueEvent->NeutrinoSet())
			continue;

		const simb::MCNeutrino trueNeutrino(trueEvent->GetNeutrino());
		const simb::MCParticle trueParticle(trueNeutrino.Nu());

		mcIsCC = ((simb::kCC == trueNeutrino.CCNC()) ? 1 : 0);
		mcPdg = trueParticle.PdgCode();
		mcNuPdg = trueNeutrino.Nu().PdgCode();
		mcIsNeutirno = true;
		mcParentPdg = 0;
		mcIsPrimary = 0;
		mcMode = trueNeutrino.Mode();

		mcVtxX = trueParticle.Vx();
		mcVtxY = trueParticle.Vy();
		mcVtxZ = trueParticle.Vz();

		mcDirX = trueParticle.Px() / trueParticle.P();
		mcDirY = trueParticle.Py() / trueParticle.P();
		mcDirZ = trueParticle.Pz() / trueParticle.P();
		mcTheta = acos(mcDirZ) * (180 / 3.1415);
		mcPhi = atan2(mcDirY, mcDirX) * (180 / 3.1415);

		mcLength = 0;

		mcEnergy = trueParticle.E();
		mcMomentum = trueParticle.P();

		completeness = 0.0;
		purity = 0.0;

		nMCHits = trueHitVector.size();
		nMCHitsU = myMatcher.CountHitsByType(geo::kU, trueHitVector);
		nMCHitsV = myMatcher.CountHitsByType(geo::kV, trueHitVector);
		nMCHitsY = myMatcher.CountHitsByType(geo::kW, trueHitVector);


		// Start Filling the PFP Neutrino information
		//==================================================
		lar_pandora::MCTruthToPFParticles::const_iterator pIter1 = matchedNeutrinos.find(trueEvent);
		if (matchedNeutrinos.end() != pIter1)
		{
			const art::Ptr<recob::PFParticle> recoParticle = pIter1->second;



			pfpPdg = recoParticle->PdgCode();
			pfpNuPdg = pfpPdg;
			pfpParentPdg = pfpPdg;
			pfpIsNeutrino = true;

			//these are not reconstructed for pfp neutrinos!
			pfpDirX = 0;
			pfpDirY = 0;
			pfpDirZ = 0;
			pfpTheta = 0;
			pfpPhi = 0;
			pfpLength = 0;
			pfpMomentum = 0;

			lar_pandoar::PFParticlesToHits::const_iterator pIter2 = recoNeutrinosToHits.find(recoParticle);
			if (recoNeutrinosToHits.end() == pIter2)
				throw cet::exception("LArPandora") << " PFParticleMonitoring::analyze --- Found a reco neutrino without any hits ";

			const lar_pandora::HitVector &recoHitVector = pIter2->second;

			//get the pfp hits
			lar_pandora::MCTruthToHits::const_iterator pIter3 = matchedNeutrinoHits.find(trueEvent);
			if (matchedNeutrinoHits.end() != pIter3)
			{
				const lar_pandora::HitVector &matchedHitVector = pIter3->second;
				nPFPHits = recoHitVector.size();
				nPFPHitsU = myMatcher.CountHitsByType(geo::kU, recoHitVector);
				nPFPHitsV = myMatcher.CountHitsByType(geo::kV, recoHitVector);
				nPFPHitsY = myMatcher.CountHitsByType(geo::kW, recoHitVector);

				nMatchedHits = matchedHitVector.size();
				nMatchedHitsU = myMatcher.CountHitsByType(geo::kU, matchedHitVector);
				nMatchedHitsV = myMatcher.CountHitsByType(geo::kV, matchedHitVector);
				nMatchedHitsY = myMatcher.CountHitsByType(geo::kW, matchedHitVector);
			}

			//get the reco vertex
			lar_pandora::PFParticlesToVertices::const_iterator pIter4 = recoParticlesToVertices.find(recoParticle);
			if (recoParticlesToVertices.end() != pIter4)
			{
				const lar_pandora::VertexVector &vertexVector = pIter4->second;
				if (!vertexVector.empty())
				{
					if (vertexVector.size() !=1 && _debug == true)
						std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;

					const art::Ptr<recob::Vertex> recoVertex = *(vertexVector.begin());
					double xyz[3] = {0.0, 0.0, 0.0};
					recoVertex->XYZ(xyz);
					pfpVtxX = xyz[0];
					pfpVtxY = xyz[1];
					pfpVtxZ = xyz[2];
				}
			}

		}
		purity = ((nPFPHits == 0) ? 0.0 : static_cast<double>(nMatchedHits) / static_cast<double>(nPFPHits));
		completeness = ((nPFPHits == 0) ? 0.0 : static_cast<double>(nMatchedHits) / static_cast<double>(nMCHits));

		//let's fill the tree with the neutrinos
		myTree->Fill();

	}//end looping neutrino map
	 //========================================================


	//start looping the track/shower maps
	//=======================================================
	lar_pandora::MCParticlesToPFParticles matchedParticles;
	lar_pandora::MCParticlesToHits matchedParticleHits;
	myMatcher.GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedParticleHits);

	// Compare true and reconstructed particles
	for (lar_pandora::MCParticlesToHits::const_iterator iter = trueParticlesToHits.begin(), iterEnd = trueParticlesToHits.end(); iter != iterEnd; ++iter)
	{
		const art::Ptr<simb::MCParticle> trueParticle = iter->first;
		const lar_pandora::HitVector &trueHitVector = iter->second;

		if (trueHitVector.empty())
			continue;

		mcPdg = trueParticle->PdgCode();
		mcIsNeutirno = false;

		try
		{
			int startT(-1);
			int endT(-1);
			myMatcher.GetStartAndEndPoints(trueParticle, startT, endT);

			mcVtxX = trueParticle->Vx(startT);
			mcVtxY = trueParticle->Vy(startT);
			mcVtxZ = trueParticle->Vz(startT);

			mcLength = myMatcher.GetLength(trueParticle, startT, endT);

			mcEnergy = trueParticle->E(startT);
			mcMomentum = trueParticle->P(startT);
			if(mcMomentum > 0.0) //don't want divide by 0 problems!
			{
				mcDirX = trueParticle->Px(startT) / mcMomentum;
				mcDirY = trueParticle->Py(startT) / mcMomentum;
				mcDirZ = trueParticle->Pz(startT) / mcMomentum;
				mcTheta = acos(mcDirZ) * (180 / 3.1415);
				mcPhi = atan2(mcDirY, mcDirX) * (180 / 3.1415);
			}
		}
		catch (cet::exception &e) {
		}

		// Get the true parent neutrino
		//=============================
		lar_pandora::MCParticlesToMCTruth::const_iterator nuIter = particlesToTruth.find(trueParticle);
		if (particlesToTruth.end() == nuIter)
			throw cet::exception("LArPandora") << " PFParticleMonitoring::analyze --- Found a true particle without any ancestry information ";

		const art::Ptr<simb::MCTruth> trueEvent = nuIter->second;
		if (trueEvent->NeutrinoSet())
		{
			const simb::MCNeutrino neutrino = trueEvent->GetNeutrino();
			mcNuPdg = neutrino.Nu().PdgCode();
			mcIsCC = ((simb::kCC == neutrino.CCNC()) ? 1 : 0);
			mcMode = neutrino.Mode();
		}

		// Get the true 'parent' and 'primary' particles
		//================================================
		try
		{
			const art::Ptr<simb::MCParticle> parentParticle(LArPandoraHelper::GetParentMCParticle(trueParticleMap, trueParticle));
			const art::Ptr<simb::MCParticle> primaryParticle(LArPandoraHelper::GetFinalStateMCParticle(trueParticleMap, trueParticle));
			mcParentPdg = ((parentParticle != trueParticle) ? parentParticle->PdgCode() : 0);
			//mcPrimaryPdg = primaryParticle->PdgCode();
			mcIsPrimary = (primaryParticle == trueParticle);//this returns the ID of the primary particle ****
		}
		catch (cet::exception &e) {
		}

		// Count number of available hits
		// Match true and reconstructed hits
		nMCHits = trueHitVector.size();
		nMCHitsU = myMatcher.CountHitsByType(geo::kU, trueHitVector);
		nMCHitsV = myMatcher.CountHitsByType(geo::kV, trueHitVector);
		nMCHitsY = myMatcher.CountHitsByType(geo::kW, trueHitVector);

		//Now we start working with the matched pfpartciles
		//========================================================
		lar_pandora::MCParticlesToPFParticles::const_iterator pIter1 = matchedParticles.find(trueParticle);
		if (matchedParticles.end() != pIter1)
		{
			const art::Ptr<recob::PFParticle> recoParticle = pIter1->second;
			pfpPdg = recoParticle->PdgCode();
			pfpNuPdg = LArPandoraHelper::GetParentNeutrino(recoParticleMap, recoParticle);
			if(pfpPdg == 12 || pfpPdg == 14) {pfpIsNeutrino = true; }
			//else{pfpIsNeutrino == false; }
			pfpIsPrimary = LArPandoraHelper::IsFinalState(recoParticleMap, recoParticle);

			const art::Ptr<recob::PFParticle> parentParticle = LArPandoraHelper::GetParentPFParticle(recoParticleMap, recoParticle);
			pfpParentPdg = parentParticle->PdgCode();

			// const art::Ptr<recob::PFParticle> primaryParticle = LArPandoraHelper::GetFinalStatePFParticle(recoParticleMap, recoParticle);
			// pfpPrimaryPdg = primaryParticle->PdgCode();

			//Find the matched hits!
			lar_pandora::PFParticlesToHits::const_iterator pIter2 = recoParticlesToHits.find(recoParticle);
			if (recoParticlesToHits.end() == pIter2)
				throw cet::exception("LArPandora") << " PFParticleMonitoring::analyze --- Found a reco particle without any hits ";

			const lar_pandora::HitVector &recoHitVector = pIter2->second;

			lar_pandora::MCParticlesToHits::const_iterator pIter3 = matchedParticleHits.find(trueParticle);
			if (matchedParticleHits.end() == pIter3)
				throw cet::exception("LArPandora") << " PFParticleMonitoring::analyze --- Found a matched true particle without matched hits ";

			const lar_pandora::HitVector &matchedHitVector = pIter3->second;

			nPFPHits = recoHitVector.size();
			nPFPHitsU = myMatcher.CountHitsByType(geo::kU, recoHitVector);
			nPFPHitsV = myMatcher.CountHitsByType(geo::kV, recoHitVector);
			nPFPHitsY = myMatcher.CountHitsByType(geo::kW, recoHitVector);

			nMatchedHits = matchedHitVector.size();
			nMatchedHitsU = myMatcher.CountHitsByType(geo::kU, matchedHitVector);
			nMatchedHitsV = myMatcher.CountHitsByType(geo::kV, matchedHitVector);
			nMatchedHitsY = myMatcher.CountHitsByType(geo::kW, matchedHitVector);

			//Find the reconstructed vertices!
			lar_pandora::PFParticlesToVertices::const_iterator pIter4 = recoParticlesToVertices.find(recoParticle);
			if (recoParticlesToVertices.end() != pIter4)
			{
				const lar_pandora::VertexVector &vertexVector = pIter4->second;
				if (!vertexVector.empty())
				{
					if (vertexVector.size() !=1 && _debug == true)
						std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;

					const art::Ptr<recob::Vertex> recoVertex = *(vertexVector.begin());
					double xyz[3] = {0.0, 0.0, 0.0};
					recoVertex->XYZ(xyz);
					pfpVtxX = xyz[0];
					pfpVtxY = xyz[1];
					pfpVtxZ = xyz[2];
				}
			}

			//stuff with directions - also gets track and shower specific!
			//===========================================================
			//tracks!
			//============================================================
			lar_pandora::PFParticlesToTracks::const_iterator pIter5 = recoParticlesToTracks.find(recoParticle);
			if (recoParticlesToTracks.end() != pIter5)
			{
				const lar_pandora::TrackVector &trackVector = pIter5->second;
				if (!trackVector.empty())
				{
					if (trackVector.size() !=1 && _debug == true)
						std::cout << " Warning: Found particle with more than one associated track " << std::endl;

					const art::Ptr<recob::Track> recoTrack = *(trackVector.begin());
					const TVector3 &vtxDirection = recoTrack->VertexDirection();

					pfpDirX = vtxDirection.x();
					pfpDirY = vtxDirection.y();
					pfpDirZ = vtxDirection.z();
					pfpTheta = acos(pfpDirZ) * (180 / 3.1415);
					pfpPhi = atan2(pfpDirY, pfpDirX) * (180 / 3.1415);
					pfpLength = recoTrack->Length();
					pfpMomentum = recoTrack->VertexMomentum();
				}
			}//end looping tracks
			 //showers!
			 //===================================================================================
			lar_pandora::PFParticlesToShowers::const_iterator pIter6 = recoParticlesToShowers.find(recoParticle);
			if (recoParticlesToShowers.end() != pIter6)
			{
				const lar_pandora::ShowerVector &showerVector = pIter6->second;
				if (!showerVector.empty())
				{
					if (showerVector.size() !=1 && _debug == true)
						std::cout << " Warning: Found particle with more than one associated shower " << std::endl;

					const art::Ptr<recob::Shower> recoShower = *(showerVector.begin());
					const TVector3 &vtxDirection = recoShower->Direction();

					pfpDirX = vtxDirection.x();
					pfpDirY = vtxDirection.y();
					pfpDirZ = vtxDirection.z();
					pfpTheta = acos(pfpDirZ) * (180 / 3.1415);
					pfpPhi = atan2(pfpDirY, pfpDirX) * (180 / 3.1415);
					pfpLength = recoShower->Length();
					const int bestplane = recoShower->best_plane();
					pfpMomentum = recoShower->Energy().at(bestplane);

					pfpOpenAngle = recoShower->OpenAngle();
				}
			}//end looping tracks
			purity = ((nPFPHits == 0) ? 0.0 : static_cast<double>(nMatchedHits) / static_cast<double>(nPFPHits));
			completeness = ((nPFPHits == 0) ? 0.0 : static_cast<double>(nMatchedHits) / static_cast<double>(nMCHits));
		}
		myTree->Fill();
	} //end looping track/shower map

//
// //==========================================================================
//      if(_is_mc == true)
//      {
//              matchinghelper.GetRecoToTrueMatches(e,
//                                                  _pfp_producer,
//                                                  _spacepointLabel,
//                                                  _geantModuleLabel,
//                                                  _hitfinderLabel,
//                                                  matchedMCToPFParticles,
//                                                  matchedParticleHits);
//      } //end if mc
//
// //loop over all matched particles
//      auto const &pfparticle_handle = e.getValidHandle<std::vector<recob::PFParticle> >(_pfp_producer);
//
//      for (auto const& iter : matchedMCToPFParticles)
//      {
//              art::Ptr<simb::MCParticle>  mc_part = iter.first;// The MCParticle
//              art::Ptr<recob::PFParticle> pf_part = iter.second; // The matched PFParticle
//
//              //fill MC info function
//
//              //fill Reco info function
//              //check if neutrino or track/shower
//              // int _pfpPdg = pf_pard->PdgCode()
//              //               if(_pfpPdg == 12 || _pfpPdg == 14)
//              // {
//              //
//              // }
//              // if(_pfpPdg == 11 || _pfpPdg == 13)
//              // {
//              //      //showers
//              //      if(_pfpPdg == 11)
//              //      {
//              //              art::FindOneP(recob::Shower) shower_for_pfp(pfparticle_handle, e, _pfp_producer);
//              //              auto const & shwr = shower_for_pfp.at(iter);
//              //      }
//              //      //tracks
//              //      if(_pfpPdg == 13)
//              //      {
//              //              art::FindOneP(recob::Track) track_for_pfp(pfparticle_handle, e, _pfp_producer);
//              //              auto const & trk = track_for_pfp.at(iter);
//              //      }
//              // }
//
//      }
//
	if(_is_data == true)
	{
		//I need to just fill all pfp information!
		std::cout << "We're looking at data!" << std::endl;
	}


}

void XSec::endSubRun(const art::SubRun& sr) {
	//probably want to fill the tree here
	std::cout << "[XSec_Module] End Running" << std::endl;

}

}//end namespace

DEFINE_ART_MODULE(nue_xsec::XSec)
