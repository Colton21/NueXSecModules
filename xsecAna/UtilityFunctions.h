#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include "AnaHelper.h"

#include "GeometryHelper.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include <algorithm>

namespace xsecAna {

class utility {

public:

/**
 *  @brief Returns number of hits on each plane for a TPC obj
 *
 *  @param e the ART event
 *  @param _particleLabel the PFP procuder module
 *  @param track_v the TPC object (vector of tracks)
 *  @param nhits_u number of hits in the U plane
 *  @param nhits_v number of hits in the V plane
 *  @param nhits_w number of hits in the W plane  */
static void GetNumberOfHitsPerPlane(art::Event const & e, std::string _particleLabel, lar_pandora::TrackVector track_v, int & nhits_u, int & nhits_v, int & nhits_w );

/**
 *  @brief Returns number of hits on each plane for a TPC obj
 *
 *  @param e the ART event
 *  @param _particleLabel the PFP procuder module
 *  @param track_v the TPC object (vector of tracks)
 *  @param nhits_u number of hits in the U plane
 *  @param nhits_v number of hits in the V plane
 *  @param nhits_w number of hits in the W plane  */
static void GetNumberOfHitsPerPlane(art::Event const & e, std::string _particleLabel, art::Ptr<recob::Track> track, int & nhits_u, int & nhits_v, int & nhits_w );

/**
 *  @brief Returns number of hits on each plane for a TPC obj
 *
 *  @param e the ART event
 *  @param _particleLabel the PFP procuder module
 *  @param track_v the TPC object (vector of tracks)
 *  @param nhits_u number of hits in the U plane
 *  @param nhits_v number of hits in the V plane
 *  @param nhits_w number of hits in the W plane  */
static void GetNumberOfHitsPerPlane(art::Event const & e, std::string _particleLabel, lar_pandora::ShowerVector shower_v, int & nhits_u, int & nhits_v, int & nhits_w );

/**
 *  @brief Returns number of hits on each plane for a TPC obj
 *
 *  @param e the ART event
 *  @param _particleLabel the PFP procuder module
 *  @param track_v the TPC object (vector of tracks)
 *  @param nhits_u number of hits in the U plane
 *  @param nhits_v number of hits in the V plane
 *  @param nhits_w number of hits in the W plane  */
static void GetNumberOfHitsPerPlane(art::Event const & e, std::string _particleLabel, art::Ptr<recob::Shower> shower, int & nhits_u, int & nhits_v, int & nhits_w );

/**
 *  @brief Given a vector of hits associated to the PFP (or a track, indeed you are just passing the hits), returns the tracking eff and purity
 *
 *  @param recoHits the input vector of reconstructed hits
 *  @param trackPurity the output track (PFP, whatever) purity
 *  @param trackEfficiency the output track (PFP, whatever) efficiency
 */
static void GetTrackPurityAndEfficiency( lar_pandora::HitVector recoHits, double & trackPurity, double & trackEfficiency );


static void ConstructShowerdQdX(xsecAna::GeometryHelper geoHelper, bool is_data, std::map <art::Ptr<recob::Cluster>, std::vector<art::Ptr < recob::Hit> > > ClusterToHitsMap,
                                std::vector<art::Ptr<recob::Cluster> > clusters, double _dQdxRectangleLength, double _dQdxRectangleWidth,
				const art::Ptr<recob::Shower> shower, std::vector< std::vector < double > > shower_cluster_dqdx, bool _verbose);

};

}//end namespace

#endif
