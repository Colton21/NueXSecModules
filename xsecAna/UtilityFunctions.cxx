#include "UtilityFunctions.h"

namespace xsec_ana {

//______________________________________________________________________________

void utility::GetNumberOfHitsPerPlane(art::Event const & e,
                                      std::string _particleLabel,
                                      lar_pandora::TrackVector track_v,
                                      int & nhits_u,
                                      int & nhits_v,
                                      int & nhits_w ) {

	nhits_u = 0;
	nhits_v = 0;
	nhits_w = 0;

	lar_pandora::TrackVector trackVector;
	lar_pandora::TracksToHits tracksToHits;
	lar_pandora::LArPandoraHelper::CollectTracks( e, _particleLabel, trackVector, tracksToHits );

	// Loop over the tracks in this TPC Object
	for (unsigned int t = 0; t < track_v.size(); t++) {

		// Get the hits associated with the track
		lar_pandora::HitVector hit_v = tracksToHits.at(track_v[t]);

		// Check where the hit is coming from
		for (unsigned int h = 0; h < hit_v.size(); h++) {

			if (hit_v[h]->View() == 0) nhits_u++;
			if (hit_v[h]->View() == 1) nhits_v++;
			if (hit_v[h]->View() == 2) nhits_w++;

		}
	}

}

//______________________________________________________________________________

void utility::GetNumberOfHitsPerPlane(art::Event const & e,
                                      std::string _particleLabel,
                                      lar_pandora::ShowerVector shower_v,
                                      int & nhits_u,
                                      int & nhits_v,
                                      int & nhits_w ) {

	nhits_u = 0;
	nhits_v = 0;
	nhits_w = 0;

	lar_pandora::ShowerVector showerVector;
	lar_pandora::ShowersToHits showersToHits;
	lar_pandora::LArPandoraHelper::CollectShowers( e, _particleLabel, showerVector, showersToHits );

	// Loop over the showers in this TPC Object
	for (unsigned int t = 0; t < shower_v.size(); t++) {

		// Get the hits associated with the track
		lar_pandora::HitVector hit_v = showersToHits.at(shower_v[t]);

		// Check where the hit is coming from
		for (unsigned int h = 0; h < hit_v.size(); h++) {

			if (hit_v[h]->View() == 0) nhits_u++;
			if (hit_v[h]->View() == 1) nhits_v++;
			if (hit_v[h]->View() == 2) nhits_w++;

		}
	}

}

//___________________________________________________________________________________________________
void utility::GetTrackPurityAndEfficiency( lar_pandora::HitVector recoHits, double & trackPurity, double & trackEfficiency ) {

	art::ServiceHandle<cheat::BackTracker> bt;

	// map from geant track id to true track deposited energy
	std::map<int,double> trkidToIDE;

	for(size_t h = 0; h < recoHits.size(); h++) {

		art::Ptr<recob::Hit> recoHit = recoHits[h];
		std::vector<sim::TrackIDE> eveIDs = bt->HitToEveID(recoHit);

		for(size_t e = 0; e < eveIDs.size(); ++e) {
			//std::cout<<"[Hit "<< h<<"] hit plane: "<<recoHit->View()<<" "<<e<<" "<<eveIDs[e].trackID<<" "<<eveIDs[e].energy<<" "<<eveIDs[e].energyFrac<<"   pdg "<< (bt->TrackIDToParticle(eveIDs[e].trackID))->PdgCode()<<"   process "<<(bt->TrackIDToParticle(eveIDs[e].trackID))->Process()<<std::endl;
			trkidToIDE[eveIDs[e].trackID] += eveIDs[e].energy;
		}
	}

	double maxe = -1;
	double tote = 0;
	int trackid;
	for(auto const& ii : trkidToIDE) {
		tote += ii.second;
		if ((ii.second)>maxe) {
			maxe = ii.second;
			trackid = ii.first;
		}
	}

	if (tote>0) {
		trackPurity = maxe/tote;
	}

	std::vector<sim::IDE> vide(bt->TrackIDToSimIDE(trackid));
	double totalEnergyFromMainTrack = 0;
	for (const sim::IDE& ide: vide) {
		totalEnergyFromMainTrack += ide.energy;
	}

	trackEfficiency = maxe/(totalEnergyFromMainTrack); //totalEnergyFromMainTrack includes both inductions and collection energies

	return;
}

}//end namespace
