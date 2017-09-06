#ifndef TPCOBJECTHELPER_H
#define TPCOBJECTHELPER_H

#include "AnaHelper.h"

namespace xsec_ana {

class tpcobjecthelper {

public:

/**
 *  @brief Gets all the PFPs for a single Pandora slice, given the neutrino PFP as input
 *
 *  @param pfParticleList the list of PFP
 *  @param particle the neutrino PFP, input
 *  @param pfp_v output, a vector of PFP (the TPC object) */
void CollectPFP(lar_pandora::PFParticleVector pfParticleList, art::Ptr<recob::PFParticle> particle, lar_pandora::PFParticleVector &pfp_v);

/**
 *  @brief Gets all the tracks and PFP for a single Pandora slice
 *
 *  @param pfParticleToTrackMap map from PFP to tracks
 *  @param pfParticleToShowerMap map from PFP to showers
 *  @param pfp_v input, a vector of PFP (the TPC object)
 *  @param track_v output, a vector of tracks (the TPC object)
 *  @param shower_v output, a vector of showers (the TPC object) */
void CollectTracksAndShowers(lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticlesToShowers pfParticleToShowerMap, lar_pandora::PFParticleVector pfp_v, lar_pandora::TrackVector &track_v, lar_pandora::ShowerVector &shower_v);

/**
 *  @brief Gets the pfp, track and shower multiplicity for a neutrino PFP
 *
 *  @param pfParticleList the list of PFP
 *  @param particle the input neutrino PFP
 *  @param pfp_v input, a vector of PFP (the TPC object)
 *  @param p output, multiplicity in number of PFPs
 *  @param t output, multiplicity in number of tracks
 *  @param s output, multiplicity in number of showers */
void GetMultiplicity(lar_pandora::PFParticleVector pfParticleList, lar_pandora::PFParticleVector pfp_v, art::Ptr<recob::PFParticle> particle, int & p, int & t, int & s);

/**
 *  @brief Constructs TPC objects using Pandora PFP slices
 *
 *  @param pfParticleList the list of PFP
 *  @param pfParticleToTrackMap map from PFP to tracks
 *  @param pfParticleToShowerMap map from PFP to showers
 *  @param pfParticleToVertexMap map from PFP to vertices
 *  @param _pfp_producer the PFP producer module
 *  @param pfp_v_v output, a vector of vector of PFP (a vector of TPC objects)
 *  @param track_v_v output, a vector of vector of tracks (a vector of TPC objects)
 *  @param shower_v_v output, a vector of vector of showers (a vector of TPC objects)
 *  @param p_v output, multiplicity in number of PFPs
 *  @param t_v output, multiplicity in number of tracks
 *  @param s_v output, multiplicity in number of showers */
void GetTPCObjects(lar_pandora::PFParticleVector pfParticleList, lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticlesToShowers pfParticleToShowerMap, lar_pandora::PFParticlesToVertices pfParticleToVertexMap, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v, std::vector<lar_pandora::ShowerVector> & shower_v_v, std::vector<int> & p_v, std::vector<int> & t_v, std::vector<int> & s_v);


/**
 *  @brief Returns the nu PFP from a TPC object
 *
 *  @param pfp_v the TPC object (vector of PFP) */
art::Ptr<recob::PFParticle> GetNuPFP(lar_pandora::PFParticleVector pfp_v);

}; //end tpcobjecthelper

}//end namespace
