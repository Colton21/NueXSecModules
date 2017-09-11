#include "TpcObjectHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


namespace nue_xsec {

//_____________________________________________________________________________________
art::Ptr<recob::PFParticle> tpcobjecthelper::GetNuPFP(lar_pandora::PFParticleVector pfp_v){

	for (unsigned int pfp = 0; pfp < pfp_v.size(); pfp++) {

		if(lar_pandora::LArPandoraHelper::IsNeutrino(pfp_v.at(pfp))) {
			return pfp_v.at(pfp);
		}
	}
	std::cout << "[TPCObjectMaker] No neutrino PFP found." << std::endl;

	art::Ptr<recob::PFParticle> temp;
	return temp;

}

//___________________________________________________________________________________________________
void tpcobjecthelper::GetTPCObjects(lar_pandora::PFParticleVector pfParticleList,
                                    lar_pandora::PFParticlesToTracks pfParticleToTrackMap,
                                    lar_pandora::PFParticlesToShowers pfParticleToShowerMap,
                                    lar_pandora::PFParticlesToVertices pfParticleToVertexMap,
                                    std::vector<lar_pandora::PFParticleVector> & pfp_v_v,
                                    std::vector<lar_pandora::TrackVector> & track_v_v,
                                    std::vector<lar_pandora::ShowerVector> & shower_v_v,
                                    std::vector<int> & p_v, std::vector<int> & t_v, std::vector<int> & s_v) {

	track_v_v.clear();
	shower_v_v.clear();
	pfp_v_v.clear();
	p_v.clear();
	t_v.clear();
	s_v.clear();


	if (_debug) std::cout << "[TPCObjectMaker] Getting TPC Objects..." << std::endl;

	for (unsigned int n = 0; n < pfParticleList.size(); ++n) {
		const art::Ptr<recob::PFParticle> particle = pfParticleList.at(n);

		if(lar_pandora::LArPandoraHelper::IsNeutrino(particle)) {
			if (_debug) std::cout << "[TPCObjectMaker] \t Creating TPC Object " << track_v_v.size() << std::endl;

			lar_pandora::TrackVector track_v;
			lar_pandora::ShowerVector shower_v;
			lar_pandora::PFParticleVector pfp_v;
			int p, t, s;


			// Collect PFPs for this TPC object
			this->CollectPFP(pfParticleList, particle, pfp_v);

			// Collect Tracks and Showers for this TPC object
			this->CollectTracksAndShowers(pfParticleToTrackMap, pfParticleToShowerMap, pfp_v, // input
			                              track_v, shower_v);               // output

			// If filtering is on, filter the PFP for this TPC object
			lar_pandora::PFParticleVector filtered_pfp_v;
			if(_tpcobj_filter && _do_filter) {

				filtered_pfp_v = _tpcobj_filter->Filter(pfp_v, pfParticleToTrackMap, pfParticleToShowerMap, pfParticleToVertexMap);

				pfp_v = filtered_pfp_v;

				this->CollectTracksAndShowers(pfParticleToTrackMap, pfParticleToShowerMap, pfp_v, // input
				                              track_v, shower_v);         // output
			}

			// Calculate multiplicity for this TPC object
			this->GetMultiplicity(pfParticleList, pfp_v, particle, p, t, s);



			if (_debug) std::cout << "[TPCObjectMaker] \t Number of pfp for this TPC object: "    << pfp_v.size()   << std::endl;
			for (auto pfp : pfp_v) {
				if (_debug) std::cout << "[TPCObjectMaker] \t \t PFP " << pfp->Self() << " with pdg " << pfp->PdgCode();
				auto it = pfParticleToVertexMap.find(pfp);
				if (it == pfParticleToVertexMap.end()) {
					if (_debug) std::cout << " and vertex [vertex not available for this PFP]" << std::endl;
				} else {
					double xyz[3];
					(it->second)[0]->XYZ(xyz);
					if (_debug) std::cout << " and vertex " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
				}
			}
			if (_debug) {
				std::cout << "[TPCObjectMaker] \t Number of tracks for this TPC object:  " << track_v.size()  << std::endl;
				std::cout << "[TPCObjectMaker] \t Number of showers for this TPC object: " << shower_v.size() << std::endl;
				std::cout << "[TPCObjectMaker] \t Multiplicity (PFP) for this TPC object:    " << p << std::endl;
				std::cout << "[TPCObjectMaker] \t Multiplicity (Track) for this TPC object:  " << t << std::endl;
				std::cout << "[TPCObjectMaker] \t Multiplicity (Shower) for this TPC object: " << s << std::endl;
				std::cout << "[TPCObjectMaker]" << std::endl;
			}

			pfp_v_v.emplace_back(pfp_v);
			track_v_v.emplace_back(track_v);
			shower_v_v.emplace_back(shower_v);
			p_v.emplace_back(p);
			t_v.emplace_back(t);
			s_v.emplace_back(s);

		} // end if neutrino
	} // end pfp loop
}

//______________________________________________________________________________________________________________________________________
void tpcobjecthelper::CollectPFP(lar_pandora::PFParticleVector pfParticleList,
                                 art::Ptr<recob::PFParticle> particle,
                                 lar_pandora::PFParticleVector &pfp_v) {

	pfp_v.emplace_back(particle);

	// And their daughters
	const std::vector<size_t> &daughterIDs = particle->Daughters();
	if(daughterIDs.size() == 0) return;
	else {
		for (unsigned int m = 0; m < daughterIDs.size(); ++m) {
			const art::Ptr<recob::PFParticle> daughter = pfParticleList.at(daughterIDs.at(m));
			// Recursive call
			this->CollectPFP(pfParticleList, daughter, pfp_v);
		}
	}

}

//______________________________________________________________________________________________________________________________________
void tpcobjecthelper::CollectTracksAndShowers(lar_pandora::PFParticlesToTracks pfParticleToTrackMap,
                                              lar_pandora::PFParticlesToShowers pfParticleToShowerMap,
                                              lar_pandora::PFParticleVector pfp_v,
                                              lar_pandora::TrackVector &track_v,
                                              lar_pandora::ShowerVector &shower_v) {

	// Cleaning
	track_v.clear();
	shower_v.clear();

	// Loop over the PFPs
	for (auto pfp : pfp_v) {

		auto iter1 = pfParticleToTrackMap.find(pfp);
		if (iter1 != pfParticleToTrackMap.end()) {
			lar_pandora::TrackVector tracks = iter1->second;
			if (_debug) std::cout << "[TPCObjectMaker] \t PFP " << pfp->Self() << " has " << tracks.size() << " tracks ass." << std::endl;
			for (unsigned int trk = 0; trk < tracks.size(); trk++) {
				track_v.emplace_back(tracks[trk]);
			}
		}

		auto iter2 = pfParticleToShowerMap.find(pfp);
		if (iter2 != pfParticleToShowerMap.end()) {
			lar_pandora::ShowerVector showers = iter2->second;
			if (_debug) std::cout << "[TPCObjectMaker] \t PFP " << pfp->Self() << " has " << showers.size() << " showers ass." << std::endl;
			for (unsigned int s = 0; s < showers.size(); s++) {
				shower_v.emplace_back(showers[s]);
			}
		}
	}

}


//______________________________________________________________________________________________________________________________________
void tpcobjecthelper::GetMultiplicity(lar_pandora::PFParticleVector pfParticleList,
                                      lar_pandora::PFParticleVector pfp_v,
                                      art::Ptr<recob::PFParticle> particle,
                                      int & p,
                                      int & t,
                                      int & s) {

	// Input PFP has to be a neutrino
	if (!lar_pandora::LArPandoraHelper::IsNeutrino(particle)) {
		std::cerr << "[TPCObjectMaker] Using tpcobjecthelper::GetMultiplicity with a non neutrino PFP as input. Exiting now." << std::endl;
		throw std::exception();
	}

	// Initialize
	p = 0;
	t = 0;
	s = 0;

	const std::vector<size_t> &daughterIDs = particle->Daughters();
	if(daughterIDs.size() == 0) {
		if (_debug) std::cout << "[TPCObjectMaker] No daughters for this neutrino PFP." << std::endl;
		return;
	}
	else {
		for (unsigned int m = 0; m < daughterIDs.size(); ++m) {

			const art::Ptr<recob::PFParticle> daughter = pfParticleList.at(daughterIDs.at(m));

			bool found_in_tpcobj = false;
			for (auto pfp : pfp_v) {
				if (daughter == pfp)
					found_in_tpcobj = true;
			}
			if (!found_in_tpcobj) continue;

			p++;
			if (lar_pandora::LArPandoraHelper::IsTrack(daughter))  t++;
			if (lar_pandora::LArPandoraHelper::IsShower(daughter)) s++;
		}
	}

}

//__________________________________________________________________________
xsec_ana::TPCObjectOrigin tpcobjecthelper::GetSliceOrigin(std::vector<art::Ptr<recob::PFParticle> > neutrinoOriginPFP, std::vector<art::Ptr<recob::PFParticle> > cosmicOriginPFP, lar_pandora::PFParticleVector pfp_v) {

	xsec_ana::TPCObjectOrigin origin = xsec_ana::kUnknown;

	int nuOrigin     = 0;
	int cosmicOrigin = 0;

	// Loop over pfp in the slice
	for ( unsigned int i = 0; i < pfp_v.size(); i++) {

		// Loop over pfp from nu origin
		for ( unsigned int j = 0; j < neutrinoOriginPFP.size(); j++) {

			if (neutrinoOriginPFP[j] == pfp_v[i]) {
				nuOrigin++;
			}
		}

		// Loop over pfp from cosmic origin
		for ( unsigned int j = 0; j < cosmicOriginPFP.size(); j++) {

			if (cosmicOriginPFP[j] == pfp_v[i]) {
				cosmicOrigin++;
			}
		}
	}

	if (nuOrigin > 0  && cosmicOrigin == 0) origin = xsec_ana::kBeamNeutrino;
	if (nuOrigin == 0 && cosmicOrigin > 0 ) origin = xsec_ana::kCosmicRay;
	if (nuOrigin > 0  && cosmicOrigin > 0 ) origin = xsec_ana::kMixed;

	return origin;

}


}//end namepsace
