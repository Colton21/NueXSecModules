#include "TpcObjectAnalysis.h"

namespace xsec_ana {

void NueXSec::reconfigure(fhicl::ParameterSet const &p)
{
	_pfp_producer                   = p.get<std::string>("PFParticleProducer");
	_mc_ghost_producer              = p.get<std::string>("MCGhostProducer");
	_tpcobject_producer             = p.get<std::string>("TPCObjectProducer");

	_debug                          = p.get<bool>("Debug", false);
	_verbose                        = p.get<bool>("Verbose", false);
}
void NueXSecProducer::reconfigure(fhicl::ParameterSet const &p)
{
	_pfp_producer                   = p.get<std::string>("PFParticleProducer");
	_hitfinderLabel                 = p.get<std::string>("HitProducer");
	_geantModuleLabel               = p.get<std::string>("GeantModule");
	_spacepointLabel                = p.get<std::string>("SpacePointProducer");
	_neutrino_flash_match_producer  = p.get<std::string>("NeutrinoFlashMatchProducer");
	_cosmic_flash_match_producer    = p.get<std::string>("CosmicFlashMatchProducer");
	_opflash_producer_beam          = p.get<std::string>("OpFlashBeamProducer");
	_potsum_producer                = p.get<std::string>("POTSummaryProducer");
	_potsum_instance                = p.get<std::string>("POTSummaryInstance");
	_particle_id_producer           = p.get<std::string>("ParticleIDProducer");

	_vertexLabel                    =p.get<std::string>("VertexProducer");
	_trackLabel                     =p.get<std::string>("TrackProducer");
	_showerLabel                    =p.get<std::string>("ShowerProducer");

	_use_genie_info                 = p.get<bool>("UseGENIEInfo", "false");
	_minimumHitRequirement          = p.get<int>("MinimumHitRequirement", 3);

	_beam_spill_start               = p.get<double>("BeamSpillStart", 3.2);
	_beam_spill_end                 = p.get<double>("BeamSpillEnd",   4.8);

	_debug                          = p.get<bool>("Debug", false);
	_verbose                        = p.get<bool>("Verbose", false);
}


NueXSec::NueXSec(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
	myTree->Branch("TpcObjectContainer", &tpc_object_container_v, "tpc_object_container_v");
}

void NueXSecProducer::produce(art::Event & e){
	art::ServiceHandle<cheat::BackTracker> bt;
	nue_xsec::recotruehelper _recotruehelper_instance;
	xsec_ana::tpcobjecthelper _tpcobjecthelper_instance;

	run = e.id().run();
	event = e.id().event();
	bool _is_data = e.isRealData();
	bool _is_mc = !_is_data;

// Instantiate the output
	std::unique_ptr< std::vector< xsec_ana::TPCObject > >                 tpcObjectVector        (new std::vector<xsec_ana::TPCObject>);
	std::unique_ptr< art::Assns<xsec_ana::TPCObject, recob::Track> >      assnOutTPCObjectTrack  (new art::Assns<xsec_ana::TPCObject, recob::Track>      );
	std::unique_ptr< art::Assns<xsec_ana::TPCObject, recob::Shower> >     assnOutTPCObjectShower (new art::Assns<xsec_ana::TPCObject, recob::Shower>     );
	std::unique_ptr< art::Assns<xsec_ana::TPCObject, recob::PFParticle> > assnOutTPCObjectPFP    (new art::Assns<xsec_ana::TPCObject, recob::PFParticle> );


// Use LArPandoraHelper functions to collect Pandora information
	lar_pandora::PFParticleVector pfParticleList;        //vector of PFParticles
	lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, pfParticleList);

// Collect vertices, tracks and shower
	lar_pandora::VertexVector allPfParticleVertices;
	lar_pandora::PFParticlesToVertices pfParticleToVertexMap;
	lar_pandora::LArPandoraHelper::CollectVertices(e, _vertexLabel, allPfParticleVertices, pfParticleToVertexMap);
	lar_pandora::TrackVector allPfParticleTracks;
	lar_pandora::PFParticlesToTracks pfParticleToTrackMap;
	lar_pandora::LArPandoraHelper::CollectTracks(e, _trackLabel, allPfParticleTracks, pfParticleToTrackMap);
	lar_pandora::ShowerVector allPfParticleShowers;
	lar_pandora::PFParticlesToShowers pfParticleToShowerMap;
	lar_pandora::LArPandoraHelper::CollectShowers(e, _showerLabel, allPfParticleShowers, pfParticleToShowerMap);


	std::vector<lar_pandora::TrackVector     > track_v_v;
	std::vector<lar_pandora::ShowerVector    > shower_v_v;
	std::vector<lar_pandora::PFParticleVector> pfp_v_v;
	std::vector<int> p_v, t_v, s_v;

	_tpcobjecthelper_instance.tpcobjecthelper::GetTPCObjects(pfParticleList, pfParticleToTrackMap, pfParticleToShowerMap,
	                                                         pfParticleToVertexMap, pfp_v_v, track_v_v, shower_v_v, p_v, t_v, s_v);

	lar_pandora::MCParticlesToPFParticles matchedParticles;    // This is a map: MCParticle to matched PFParticle
	lar_pandora::MCParticlesToHits matchedParticleHits;
	if (_is_mc)
	{
		_recotruehelper_instance.GetRecoToTrueMatches(matchedParticles, matchedParticleHits);
	}

//reco true matching is performed here!
	std::vector<art::Ptr<recob::PFParticle> > neutrinoOriginPFP;
	std::vector<art::Ptr<recob::PFParticle> > cosmicOriginPFP;

	if(!neutrinoOriginPFP.empty()) {neutrinoOriginPFP.clear(); }
	if(!cosmicOriginPFP.empty()) {cosmicOriginPFP.clear(); }

	std::vector< std::pair< int, simb::Origin_t > > pfp_origin_v;
	if(!pfp_origin_v.empty()) {pfp_origin_v.clear(); }
	for (lar_pandora::MCParticlesToPFParticles::const_iterator iter1 = matchedParticles.begin(), iterEnd1 = matchedParticles.end();
	     iter1 != iterEnd1; ++iter1)
	{
		art::Ptr<simb::MCParticle>  mc_par = iter1->first;// The MCParticle
		art::Ptr<recob::PFParticle> pf_par = iter1->second; // The matched PFParticle

		const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());

		if (!mc_truth) {
			std::cerr << "[TPCObjectMaker] Problem with MCTruth pointer." << std::endl;
			continue;
		}
		std::pair < int, simb::Origin_t > pfp_origin (pf_par->Self(), mc_truth->Origin());
		pfp_origin_v.emplace_back(pfp_origin);
	} //end looping mc to pfp

//loop over TPC objects
	for (size_t pfparticle_vector = 0; pfparticle_vector < pfp_v_v.size(); pfparticle_vector++)
	{

		xsec_ana::TPCObject obj;

		// Set tracks
		std::vector<recob::Track> trk_v;
		//trk_v.clear();
		for (auto t : track_v_v[pfparticle_vector]) trk_v.emplace_back((*t));
		obj.SetTracks(trk_v);

		// Set showers
		std::vector<recob::Shower> shwr_v;
		//shwr_v.clear();
		for (auto s : shower_v_v[pfparticle_vector]) shwr_v.emplace_back((*s));
		obj.SetShowers(shwr_v);

		//set individual particle origins
		std::vector<simb::Origin_t> origin_v;
		// Set PFPs
		std::vector<recob::PFParticle> pfp_v;
		//pfp_v.clear();
		for (auto p : pfp_v_v[pfparticle_vector])
		{
			pfp_v.emplace_back((*p));
			const int pfp_id = p->Self();
			simb::Origin_t pfp_origin = simb::kUnknown; //this is for the case where the pfp is not matched
			bool matched = false;
			//const int pfp_pdg = p.PdgCode();
			for(auto const pp : pfp_origin_v)
			{
				if(pfp_id == pp.first)
				{
					matched = true;
					origin_v.emplace_back(pp.second);
				}
			}
			if(matched == false) {origin_v.emplace_back(pfp_origin); }
		}
		obj.SetPFPs(pfp_v);
		obj.SetParticleOrigins(origin_v);

		// Set vertex
		art::Ptr<recob::PFParticle> pfp = _tpcobjecthelper_instance.tpcobjecthelper::GetNuPFP(pfp_v_v[pfparticle_vector]);
		auto iter = pfParticleToVertexMap.find(pfp);
		if (iter != pfParticleToVertexMap.end()) {obj.SetVertex(*(iter->second[0])); }

		// Set origin
		//xsec_ana::TPCObjectOrigin origin = xsec_ana::kUnknown;
		simb::Origin_t origin = simb::kUnknown;
		if (_is_mc)
		{
			origin = _tpcobjecthelper_instance.tpcobjecthelper::GetSliceOrigin(neutrinoOriginPFP, cosmicOriginPFP, pfp_v_v[pfparticle_vector]);
		}
		obj.SetOrigin(origin);

		// Set Multiplicity
		obj.SetMultiplicity(p_v[pfparticle_vector], t_v[pfparticle_vector], s_v[pfparticle_vector]);

		tpcObjectVector->emplace_back(obj);
		util::CreateAssn(*this, e, *tpcObjectVector, track_v_v[pfparticle_vector],  *assnOutTPCObjectTrack);
		util::CreateAssn(*this, e, *tpcObjectVector, shower_v_v[pfparticle_vector], *assnOutTPCObjectShower);
		util::CreateAssn(*this, e, *tpcObjectVector, pfp_v_v[pfparticle_vector],    *assnOutTPCObjectPFP);

	} // end looping pfparticle_vectors (tpc objects)

// Put TPCObjects into the Event
	e.put(std::move(tpcObjectVector));
	e.put(std::move(assnOutTPCObjectTrack));
	e.put(std::move(assnOutTPCObjectShower));
	e.put(std::move(assnOutTPCObjectPFP));

	if (_debug) std::cout << "[TPCObjectMaker] Ends" << std::endl;

}
//***********************************
//finish constructing TPC Objects! **
//***********************************

//**************************
//start the analysis loops!
//*************************

void NueXSec::analyze(art::Event const & e)
{
	art::ServiceHandle<cheat::BackTracker> bt;

	run = e.id().run();
	event = e.id().event();
	//bool _is_data = e.isRealData();
	//bool _is_mc = !_is_data;


	// Get PFP
	art::Handle<std::vector<recob::PFParticle> > pfp_h;
	e.getByLabel(_pfp_producer,pfp_h);
	if(!pfp_h.isValid())
	{
		std::cout << "[UBXSec] PFP product " << _pfp_producer << " not found..." << std::endl;
		//throw std::exception();
	}
	if(pfp_h->empty())
	{
		std::cout << "[UBXSec] PFP " << _pfp_producer << " is empty." << std::endl;
	}
	art::FindManyP<recob::Track> tracks_from_pfp(pfp_h, e, _pfp_producer);
	art::FindManyP<recob::Shower> showers_from_pfp(pfp_h, e, _pfp_producer);

	// Get Ghosts
	art::Handle<std::vector<xsec_ana::MCGhost> > ghost_h;
	e.getByLabel(_mc_ghost_producer,ghost_h);
	if(!ghost_h.isValid()) {
		std::cout << "[UBXSec] MCGhost product " << _mc_ghost_producer << " not found..." << std::endl;
		//throw std::exception();
	}
	art::FindManyP<xsec_ana::MCGhost>   mcghost_from_pfp   (pfp_h,   e, _mc_ghost_producer);
	art::FindManyP<simb::MCParticle> mcpar_from_mcghost (ghost_h, e, _mc_ghost_producer);

	// Get TPCObjects from the Event
	art::Handle<std::vector<xsec_ana::TPCObject> > tpcobj_h;
	e.getByLabel(_tpcobject_producer, tpcobj_h);
	if (!tpcobj_h.isValid()) {
		std::cout << "[UBXSec] Cannote locate ubana::TPCObject." << std::endl;
	}
	//art::FindManyP<xsec_ana::FlashMatch> tpcobjToFlashMatchAssns(tpcobj_h, e, _neutrino_flash_match_producer);
	art::FindManyP<recob::Track>      tpcobjToTrackAssns(tpcobj_h, e, _tpcobject_producer);
	art::FindManyP<recob::Shower>     tpcobjToShowerAssns(tpcobj_h, e, _tpcobject_producer);
	art::FindManyP<recob::PFParticle> tpcobjToPFPAssns(tpcobj_h, e, _tpcobject_producer);

	std::cout << "TPC Objects in this Event: " << tpcobj_h->size() << std::endl;
	//loop over all of the tpc objects!

	int tpc_object_counter = 0;
	if(!tpc_object_container_v.empty()) {tpc_object_container_v.clear(); }

	for(size_t tpc_counter = 0; tpc_counter < tpcobj_h->size(); tpc_counter++)
	{
		const xsec_ana::TPCObject tpcobj = (*tpcobj_h)[tpc_counter];
		const int ntracks                   = tpcobj.GetNTracks();
		const int nshowers                  = tpcobj.GetNShowers();
		const int npfparticles              = tpcobj.GetNPFP();
		const simb::Origin_t tpcobj_origin  = tpcobj.GetOrigin();
		const std::vector<recob::PFParticle> pfp_v = tpcobj.GetPFPs();
		const std::vector<recob::Track>  track_v  = tpcobj.GetTracks();
		const std::vector<recob::Shower> shower_v = tpcobj.GetShowers();

		if(_verbose)
		{
			std::cout << "N PFPs   : " << npfparticles << std::endl;
			std::cout << "N Tracks : " << ntracks << std::endl;
			std::cout << "N Showers: " << nshowers << std::endl;
		}

		// Reco vertex
		double reco_nu_vtx[3];
		recob::Vertex tpcobj_nu_vtx = tpcobj.GetVertex();
		tpcobj_nu_vtx.XYZ(reco_nu_vtx);
		const double pfp_nu_vtx_x = reco_nu_vtx[0];
		const double pfp_nu_vtx_y = reco_nu_vtx[1];
		const double pfp_nu_vtx_z = reco_nu_vtx[2];
		//_slc_nuvtx_fv[slice] = (UBXSecHelper::InFV(reco_nu_vtx) ? 1 : 0);

		xsec_ana::pfp_container tpc_object_container;
		tpc_object_container.SetpfpVtxX(pfp_nu_vtx_x);
		tpc_object_container.SetpfpVtxY(pfp_nu_vtx_y);
		tpc_object_container.SetpfpVtxZ(pfp_nu_vtx_z);
		tpc_object_container.SetRunNumber(run);
		//tpc_object_container.SetSubRunNumber();
		tpc_object_container.SetEventNumber(event);
		tpc_object_container.SetIndex(tpc_object_counter);
		tpc_object_container.SetOrigin(tpcobj_origin);

		if(_verbose) {std::cout << "Number of PFP in this TPC Object: " << npfparticles << std::endl; }
		tpc_object_container.SetNumPFParticles(npfparticles);

		// Hits - we want the hits from both tracks and showers
		int nhits_u = 0;
		int nhits_v = 0;
		int nhits_w = 0;
		int total_nhits_u = 0;
		int total_nhits_v = 0;
		int total_nhits_w = 0;
		int total_nhits = 0;
		//need to sum all hits from both tracks and showers
		for(auto const track : track_v)
		{
			art::Ptr<recob::Track> & _track = track;
			xsec_ana::utility::GetNumberOfHitsPerPlane(e, _pfp_producer, _track, nhits_u, nhits_v, nhits_w);
			total_nhits_u += nhits_u;
			total_nhits_v += nhits_v;
			total_nhits_w += nhits_w;
		}
		for(auto const shower : shower_v)
		{
			art::Ptr<recob::Shower> & _shower = shower;
			xsec_ana::utility::GetNumberOfHitsPerPlane(e, _pfp_producer, _shower, nhits_u, nhits_v, nhits_w);
			total_nhits_u += nhits_u;
			total_nhits_v += nhits_v;
			total_nhits_w += nhits_w;
		}
		total_nhits = (total_nhits_u + total_nhits_v + total_nhits_w);

		tpc_object_container.SetNumPFPHits   (total_nhits);
		tpc_object_container.SetNumPFPHitsU  (total_nhits_u);
		tpc_object_container.SetNumPFPHitsV  (total_nhits_v);
		tpc_object_container.SetNumPFPHitsW  (total_nhits_w);

		int pfp_nu_counter = 0;
		//int pfpNuPdg = 0;//not set
		int pfpParentPdg = 0;
		int mode = -1;
		int ccnc = -1;
		bool is_neutrino = false;
		//bool is_primary = false; // not set

		//****************************************
		//loop over pfparticles in the tpc object
		//****************************************
		auto pfps_from_tpcobj = tpcobjToPFPAssns.at(tpc_object_counter);
		for(auto const pfp : pfps_from_tpcobj)
		//for(auto const pfp : pfp_v)
		{

			xsec_ana::ParticleContainer particle_container;

			int mcPdg = 0;
			//int mcNuPdg = 0; // not set
			//int mcParentPdg = 0; // not set
			double pfp_dir_x = 0;
			double pfp_dir_y = 0;
			double pfp_dir_z = 0;
			double pfp_theta = 0;
			double pfp_phi = 0;
			double pfp_length = 0;
			double pfp_momentum = 0;
			int pfp_hits = 0;
			int pfp_hits_u = 0;
			int pfp_hits_v = 0;
			int pfp_hits_w = 0;
			double pfp_open_angle = 0;

			double mc_vtx_x = 0;
			double mc_vtx_y = 0;
			double mc_vtx_z = 0;
			double mc_dir_x = 0;
			double mc_dir_y = 0;
			double mc_dir_z = 0;
			double mc_theta = 0;
			double mc_phi = 0;
			double mcLength = 0;
			double mcEnergy = 0;
			double mcMomentum = 0;
			//double mc_open_angle = 0; //unset

			const int pfpPdg = pfp->PdgCode();
			pfpParentPdg = pfp_v.at(pfp->Parent()).PdgCode();
			particle_container.SetpfpPdgCode(pfpPdg);
			particle_container.SetpfpNuPdgCode(pfpParentPdg);
			const int index = pfp->Self();
			particle_container.SetIndex(index);
			simb::Origin_t mcOrigin = simb::kUnknown;

			if(pfpPdg == 12 || pfpPdg == 14) {
				if(_verbose) {std::cout << "PFP Neutrino with PDG Code: " << pfpPdg << std::endl; }
				is_neutrino = true;
				pfp_nu_counter++;
			}
			particle_container.SetIsNeutrino(is_neutrino);

			// Reco vertex
			lar_pandora::VertexVector vertexVector;
			lar_pandora::PFParticlesToVertices particlesToVertices;
			lar_pandora::LArPandoraHelper::CollectVertices(e, _pfp_producer, vertexVector, particlesToVertices);

			auto iter = particlesToVertices.find(pfp);
			//auto iter = pfParticleToVertexMap.find(pfp);
			double pfp_vtx_x = 0;
			double pfp_vtx_y = 0;
			double pfp_vtx_z = 0;
			if (iter != particlesToVertices.end())
			//if(iter != pfParticleToVertexMap.end())
			{
				//lar_pandora::VertexVector vertex_v = pfParticleToVertexMap.find(pfp)->second;
				lar_pandora::VertexVector vertex_v = particlesToVertices.find(pfp)->second;
				double reco_vtx[3];
				vertex_v[0]->XYZ(reco_vtx);
				pfp_vtx_x = reco_vtx[0];
				pfp_vtx_y = reco_vtx[1];
				pfp_vtx_z = reco_vtx[2];
			}

			particle_container.SetpfpVtxX(pfp_vtx_x);
			particle_container.SetpfpVtxY(pfp_vtx_y);
			particle_container.SetpfpVtxZ(pfp_vtx_z);

			//mcghosts do accounting from pfp to mcghost to mc particle
			const std::vector<art::Ptr<MCGhost> > mcghost = mcghost_from_pfp.at(pfp.key());
			std::vector<art::Ptr<simb::MCParticle> > mcpart;
			if(mcghost.size() == 0) {std::cout << "No matched MC Ghost to PFP!" << std::endl; }
			if(mcghost.size() > 1) {std::cout << "Too many matched MC Ghost to PFP!" << std::endl; }
			if(mcghost.size() == 1)
			{
				mcpart = mcpar_from_mcghost.at(mcghost[0].key());
				const art::Ptr<simb::MCParticle> the_mcpart = mcpart.at(0);
				const art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(the_mcpart->TrackId());
				const simb::MCNeutrino mc_nu = mctruth->GetNeutrino();
				mode = mc_nu.Mode();
				ccnc = mc_nu.CCNC();

				mcOrigin = mctruth->Origin();
				mcPdg = the_mcpart->PdgCode();
				//mcNuPdg = the_mcpart->Mother();
				//mcParentPdg = the_mcpart->Mother();
				mc_vtx_x = the_mcpart->Vx();
				mc_vtx_y = the_mcpart->Vy();
				mc_vtx_z = the_mcpart->Vz();
				mc_dir_x = the_mcpart->Px();
				mc_dir_y = the_mcpart->Py();
				mc_dir_z = the_mcpart->Pz();
				mc_theta = acos(mc_dir_z) * (180 / 3.1415);
				mc_phi = atan2(mc_dir_y, mc_dir_x) * (180 / 3.1415);
				const double mc_length_x = the_mcpart->Position().X() - the_mcpart->EndPosition().X();
				const double mc_length_y = the_mcpart->Position().Y() - the_mcpart->EndPosition().Y();
				const double mc_length_z = the_mcpart->Position().Z() - the_mcpart->EndPosition().Z();
				mcLength = sqrt((mc_length_x * mc_length_x) + (mc_length_y * mc_length_y) + (mc_length_z * mc_length_z));
				mcEnergy = the_mcpart->E();
				mcMomentum = the_mcpart->P();

			}//end mcghost == 1
			particle_container.SetmcPdgCode(mcPdg);
			particle_container.SetOrigin(mcOrigin);
			particle_container.SetmcVtxX(mc_vtx_x);
			particle_container.SetmcVtxY(mc_vtx_y);
			particle_container.SetmcVtxZ(mc_vtx_z);
			particle_container.SetmcDirX(mc_dir_x);
			particle_container.SetmcDirY(mc_dir_y);
			particle_container.SetmcDirZ(mc_dir_z);
			particle_container.SetmcTheta(mc_theta);
			particle_container.SetmcPhi(mc_phi);
			particle_container.SetmcLength(mcLength);
			particle_container.SetmcEnergy(mcEnergy);
			particle_container.SetmcMomentum(mcMomentum);

			//pfp tracks
			if(pfpPdg == 13)
			{
				std::vector<art::Ptr<recob::Track> > tracks = tracks_from_pfp.at(pfp.key());
				std::cout << "[UBXSec] \t\t n tracks ass to this pfp: " << tracks.size() << std::endl;
				//we want to take the first association, right?
				const art::Ptr<recob::Track> this_track = tracks.at(0);
				pfp_dir_x = this_track->VertexDirection().X();
				pfp_dir_y = this_track->VertexDirection().Y();
				pfp_dir_z = this_track->VertexDirection().Z();
				pfp_theta = acos(pfp_dir_z) * (180 / 3.1415);
				pfp_phi = atan2(pfp_dir_y, pfp_dir_x) * (180 / 3.1415);
				pfp_length = this_track->Length();
				pfp_momentum = this_track->VertexMomentum();

				xsec_ana::utility::GetNumberOfHitsPerPlane(e, _pfp_producer, this_track, pfp_hits_u, pfp_hits_v, pfp_hits_w);
				// Check where the hit is coming from
				// for (unsigned int h = 0; h < hit_v.size(); h++) {
				//      if (hit_v[h]->View() == 0) pfp_hits_u++;
				//      if (hit_v[h]->View() == 1) pfp_hits_v++;
				//      if (hit_v[h]->View() == 2) pfp_hits_w++;
				// }
				pfp_hits = (pfp_hits_u + pfp_hits_v + pfp_hits_w);
			}//end pfp tracks
			 //pfp showers
			if(pfpPdg == 11)
			{
				std::vector<art::Ptr<recob::Shower> > showers = showers_from_pfp.at(pfp.key());
				std::cout << "[UBXSec] \t\t n showers ass to this pfp: " << showers.size() << std::endl;
				//we want to take the first association, right?
				const art::Ptr<recob::Shower> this_shower = showers.at(0);
				pfp_dir_x = this_shower->Direction().X();
				pfp_dir_y = this_shower->Direction().Y();
				pfp_dir_z = this_shower->Direction().Z();
				pfp_theta = acos(pfp_dir_z) * (180 / 3.1415);
				pfp_phi = atan2(pfp_dir_y, pfp_dir_x) * (180 / 3.1415);
				pfp_length = this_shower->Length();
				pfp_momentum = this_shower->Energy().at(this_shower->best_plane());
				pfp_open_angle = this_shower->OpenAngle();

				xsec_ana::utility::GetNumberOfHitsPerPlane(e, _pfp_producer, this_shower, pfp_hits_u, pfp_hits_v, pfp_hits_w);
				// Check where the hit is coming from
				// for (unsigned int h = 0; h < hit_v.size(); h++) {
				//      if (hit_v[h]->View() == 0) pfp_hits_u++;
				//      if (hit_v[h]->View() == 1) pfp_hits_v++;
				//      if (hit_v[h]->View() == 2) pfp_hits_w++;
				// }
				pfp_hits = (pfp_hits_u + pfp_hits_v + pfp_hits_w);
			}//end pfp showers

			particle_container.SetpfpDirX(pfp_dir_x);
			particle_container.SetpfpDirY(pfp_dir_y);
			particle_container.SetpfpDirZ(pfp_dir_z);
			particle_container.SetpfpTheta(pfp_theta);
			particle_container.SetpfpPhi(pfp_phi);
			particle_container.SetpfpLength(pfp_length);
			particle_container.SetpfpMomentum(pfp_momentum);
			particle_container.SetpfpOpenAngle(pfp_open_angle);
			particle_container.SetNumPFPHits(pfp_hits);
			particle_container.SetNumPFPHitsU(pfp_hits_u);
			particle_container.SetNumPFPHitsV(pfp_hits_v);
			particle_container.SetNumPFPHitsW(pfp_hits_w);

			tpc_object_container.AddParticle(particle_container);

		}//end loop over pfp in tpcobject
		tpc_object_container.SetNumPFPNeutrinos(pfp_nu_counter);
		tpc_object_container.SetMode(mode);
		tpc_object_container.SetIsCC(ccnc);

		tpc_object_container_v.emplace_back(tpc_object_container);
		tpc_object_counter++;
	}//end loop tpc objects

	//fill root tree per event
	myTree->Fill();

}//end analyze

void NueXSec::endSubRun(art::SubRun const & sr) {
	//probably want to fill the tree here
	std::cout << "[XSec_Module] End Running" << std::endl;

}

}//end namespace

DEFINE_ART_MODULE(xsec_ana::NueXSec)
