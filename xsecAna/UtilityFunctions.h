#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include "AnaHelper.h"

#include "GeometryHelper.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibService.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

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

//this uses the calibration constants per plane as inputs and calculates the deposited energy for a track/shower based on the hit total ADCs
static void GetEnergyPerPlane(art::Event const & e,
                              std::string _particleLabel,
                              art::Ptr<recob::Shower> shower,
                              double & calibration_u,
                              double & calibration_v,
                              double & calibration_w,
                              double & energy_u,
                              double & energy_v,
                              double & energy_w );

static void GetEnergyPerPlane(art::Event const & e,
                              std::string _particleLabel,
                              art::Ptr<recob::Track> track,
                              double & calibration_u,
                              double & calibration_v,
                              double & calibration_w,
                              double & energy_u,
                              double & energy_v,
                              double & energy_w );

static void GetTrackPurityAndEfficiency( lar_pandora::HitVector recoHits, double & trackPurity, double & trackEfficiency );


static void ConstructShowerdQdX(xsecAna::GeometryHelper geoHelper, bool is_data, std::map <art::Ptr<recob::Cluster>, std::vector<art::Ptr < recob::Hit> > > ClusterToHitsMap,
                                std::vector<art::Ptr<recob::Cluster> > clusters, double _dQdxRectangleLength, double _dQdxRectangleWidth,
                                const art::Ptr<recob::Shower> shower, std::vector< std::vector < double > > & shower_cluster_dqdx,
                                std::vector< std::vector < double > > & shower_cluster_dq,
                                std::vector< std::vector < double > > & shower_cluster_dx,  bool _verbose);

static void ConvertdEdX(std::vector< std::vector < double > > & shower_cluster_dqdx, std::vector<double> & shower_dEdx);

static void ConstructShowerdQdXAlternative(xsecAna::GeometryHelper geoHelper, bool is_data, std::map <art::Ptr<recob::Cluster>, std::vector<art::Ptr< recob::Hit> > > ClusterToHitsMap,
                                           std::vector<art::Ptr<recob::Cluster> > clusters, double _dQdxRectangleLength, double _dQdxRectangleWidth,
                                           const art::Ptr<recob::Shower> shower, std::vector< std::vector < double > > & shower_cluster_dqdx,
                                           std::vector< std::vector < double > > & shower_cluster_dq, std::vector< std::vector < double > > & shower_cluster_dx,
                                           std::vector<double> & dqdx_cali,
                                           bool _verbose);

};

}//end namespace

#endif
