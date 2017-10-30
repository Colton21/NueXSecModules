#include "UtilityFunctions.h"

namespace xsecAna {

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
                                      art::Ptr<recob::Track> track,
                                      int & nhits_u,
                                      int & nhits_v,
                                      int & nhits_w ) {

	nhits_u = 0;
	nhits_v = 0;
	nhits_w = 0;

	lar_pandora::TrackVector trackVector;
	lar_pandora::TracksToHits tracksToHits;
	lar_pandora::LArPandoraHelper::CollectTracks( e, _particleLabel, trackVector, tracksToHits );

	// Get the hits associated with the track
	lar_pandora::HitVector hit_v = tracksToHits.at(track);

	// Check where the hit is coming from
	for (unsigned int h = 0; h < hit_v.size(); h++) {

		if (hit_v[h]->View() == 0) nhits_u++;
		if (hit_v[h]->View() == 1) nhits_v++;
		if (hit_v[h]->View() == 2) nhits_w++;

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

//______________________________________________________________________________

void utility::GetNumberOfHitsPerPlane(art::Event const & e,
                                      std::string _particleLabel,
                                      art::Ptr<recob::Shower> shower,
                                      int & nhits_u,
                                      int & nhits_v,
                                      int & nhits_w ) {

	nhits_u = 0;
	nhits_v = 0;
	nhits_w = 0;

	lar_pandora::ShowerVector showerVector;
	lar_pandora::ShowersToHits showersToHits;
	lar_pandora::LArPandoraHelper::CollectShowers( e, _particleLabel, showerVector, showersToHits );

	// Get the hits associated with the track
	lar_pandora::HitVector hit_v = showersToHits.at(shower);

	// Check where the hit is coming from
	for (unsigned int h = 0; h < hit_v.size(); h++) {

		if (hit_v[h]->View() == 0) nhits_u++;
		if (hit_v[h]->View() == 1) nhits_v++;
		if (hit_v[h]->View() == 2) nhits_w++;

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
			//std::cout<<"[Hit "<< h<<"] hit plane: "<<recoHit->View()<<" "<<e<<" "<<eveIDs[e].trackID<<" "<<eveIDs[e].energy<<" "
			//<<eveIDs[e].energyFrac<<"   pdg "<< (bt->TrackIDToParticle(eveIDs[e].trackID))->PdgCode()<<"   process "
			//<<(bt->TrackIDToParticle(eveIDs[e].trackID))->Process()<<std::endl;
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

void utility::ConstructShowerdQdX(std::map <art::Ptr<recob::Cluster>, art::Ptr<std::vector<recob::Hit> > > ClusterToHitsMap,
                                  std::vector<art::Ptr<recob::Cluster> > clusters, const art::Ptr<recob::Shower> shower, bool _verbose)
{

	const double _gain = 0;

	detinfo::DetectorProperties const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
	const double drift = detprop->DriftVelocity() * 1e-3;

	// TODO Use variable from detector properties!
	// To get the time in ns -> 4.8 ms / 9600 ticks * 1e6 = 500
	// 0.3 wire spacing

	const double fromTickToNs = 4.8 / detprop->ReadOutWindowSize() * 1e6;
	const double wireSpacing = 0.3;

	const int n_clusters = clusters.size();

	if(_verbose) {std::cout << "[dQdx] Clusters size " << n_clusters << std::endl; }
	for(auto const cluster : clusters)
	{
		art::Ptr< std::vector < recob::Hits > > > hit_v = ClusterToHitsMap.find(cluster)->second;
		const int start_wire = cluster->StartWire();
		const double start_position = start_wire * wireSpacing;
		std::cout << start_wire << std::endl;
	}

}

}//end namespace
