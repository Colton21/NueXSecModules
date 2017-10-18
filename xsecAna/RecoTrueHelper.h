#ifndef RECOTRUEHELPER_H
#define RECOTRUEHELPER_H


#include "AnaHelper.h"

typedef std::set<art::Ptr<recob::PFParticle> > PFParticleSet;
typedef std::set<art::Ptr<simb::MCParticle> > MCParticleSet;
typedef std::set< art::Ptr<simb::MCTruth> > MCTruthSet;


namespace nue_xsec {

class recotruehelper {

private:

lar_pandora::HitsToMCParticles _trueHitsToParticles;   ///< A map from recon hits to MCParticles
lar_pandora::PFParticlesToHits _recoParticlesToHits; ///< A map from PFParticles to recon hits

public:

/// Default constructor
//recotruehelper();

/// Default destructor
//~recotruehelper(){
//}

//I should export this as a fcl parameter in the future
bool _debug             = false;
bool _verbose           = false;
bool _recursiveMatching = false;

/// Configure function parameters
void Configure(art::Event const & e, std::string _pfp_producer,
               std::string _spacepointLabel, std::string _hitfinderLabel, std::string _geantModuleLabel);


/**
 *  @brief  Build mapping from true neutrinos to hits
 *
 *  @param truthToParticles  the input mapping from true event to true particles
 *  @param trueParticlesToHits  the input mapping from true particles to hits
 *  @param trueNeutrinosToHits  the output mapping from trues event to hits
 *  @param trueHitsToNeutrinos  the output mappign from hits to true events
 */
void BuildTrueNeutrinoHitMaps(const lar_pandora::MCTruthToMCParticles &truthToParticles, const lar_pandora::MCParticlesToHits &trueParticlesToHits,
                              lar_pandora::MCTruthToHits &trueNeutrinosToHits, lar_pandora::HitsToMCTruth &trueHitsToNeutrinos) const;

/**
 *  @brief  Build mapping from reconstructed neutrinos to hits
 *
 *  @param recoParticleMap  the input mapping from reconstructed particle and particle ID
 *  @param recoParticlesToHits  the input mapping from reconstructed particles to hits
 *  @param recoNeutrinosToHits  the output mapping from reconstructed particles to hits
 *  @param recoHitsToNeutrinos  the output mapping from reconstructed hits to particles
 */
void BuildRecoNeutrinoHitMaps(const lar_pandora::PFParticleMap &recoParticleMap, const lar_pandora::PFParticlesToHits &recoParticlesToHits,
                              lar_pandora::PFParticlesToHits &recoNeutrinosToHits, lar_pandora::HitsToPFParticles &recoHitsToNeutrinos) const;


/**
 *  @brief Returns matching between true and reconstructed particles
 *
 *  @param matchedParticles the output matches between reconstructed and true particles
 *  @param matchedHits the output matches between reconstructed particles and hits
 */
void GetRecoToTrueMatches(lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits);


/**
 *  @brief Perform matching between true and reconstructed neutrino events
 *
 *  @param recoNeutrinosToHits  the mapping from reconstructed neutrino events to hits
 *  @param trueHitsToNeutrinos  the mapping from hits to true neutrino events
 *  @param matchedNeutrinos  the output matches between reconstructed and true neutrinos
 *  @param matchedNeutrinoHits  the output matches between reconstructed neutrinos and hits
 */
void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoNeutrinosToHits, const lar_pandora::HitsToMCTruth &trueHitsToNeutrinos,
                          lar_pandora::MCTruthToPFParticles &matchedNeutrinos, lar_pandora::MCTruthToHits &matchedNeutrinoHits) const;

/**
 *  @brief Perform matching between true and reconstructed neutrino events
 *
 *  @param recoNeutrinosToHits  the mapping from reconstructed neutrino events to hits
 *  @param trueHitsToNeutrinos  the mapping from hits to true neutrino events
 *  @param matchedNeutrinos  the output matches between reconstructed and true neutrinos
 *  @param matchedNeutrinoHits  the output matches between reconstructed neutrinos and hits
 *  @param recoVeto  the veto list for reconstructed particles
 *  @param trueVeto  the veto list for true particles
 */
void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoNeutrinosToHits, const lar_pandora::HitsToMCTruth &trueHitsToNeutrinos,
                          lar_pandora::MCTruthToPFParticles &matchedNeutrinos, lar_pandora::MCTruthToHits &matchedNeutrinoHits,
                          PFParticleSet &recoVeto, MCTruthSet &trueVeto) const;

/**
 *  @brief Perform matching between true and reconstructed particles
 *
 *  @param recoParticlesToHits the mapping from reconstructed particles to hits
 *  @param trueHitsToParticles the mapping from hits to true particles
 *  @param matchedParticles the output matches between reconstructed and true particles
 *  @param matchedHits the output matches between reconstructed particles and hits
 */
void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                          lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits) const;

/**
 *  @brief Perform matching between true and reconstructed particles
 *
 *  @param recoParticlesToHits the mapping from reconstructed particles to hits
 *  @param trueHitsToParticles the mapping from hits to true particles
 *  @param matchedParticles the output matches between reconstructed and true particles
 *  @param matchedHits the output matches between reconstructed particles and hits
 *  @param recoVeto the veto list for reconstructed particles
 *  @param trueVeto the veto list for true particles
 */
void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                          lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits,
                          PFParticleSet &recoVeto, MCParticleSet &trueVeto) const;

/**
 *  @brief Build particle maps for reconstructed particles
 *
 *  @param particleVector the input vector of reconstructed particles
 *  @param particleMap the output mapping between reconstructed particles and particle ID
 */
void BuildRecoParticleMap(const lar_pandora::PFParticleVector &particleVector, lar_pandora::PFParticleMap &particleMap) const;

/**
 *  @brief Build particle maps for true particles
 *
 *  @param particleVector the input vector of true particles
 *  @param particleMap the output mapping between true particle and true track ID
 */
void BuildTrueParticleMap(const lar_pandora::MCParticleVector &particleVector, lar_pandora::MCParticleMap &particleMap) const;

/**
 *  @brief Count the number of reconstructed hits in a given wire plane
 *
 *  @param view the wire plane ID
 *  @param hitVector the input vector of reconstructed hits
 */
int CountHitsByType(const int view, const lar_pandora::HitVector &hitVector) const;

/**
 *  @brief Find the start and end points of the true particle in the active region of detector
 *
 *  @param trueParticle the input true particle
 *  @param startT  the true start point
 *  @param endT  the true end point
 */
void GetStartAndEndPoints(const art::Ptr<simb::MCParticle> trueParticle, int &startT, int &endT) const;

/**
 *  @brief Find the length of the true particle trajectory through the active region of the detector
 *
 *  @param trueParticle the input true particle
 *  @param startT  the true start point
 *  @param endT  the true end point
 */
double GetLength(const art::Ptr<simb::MCParticle> trueParticle, const int startT, const int endT) const;

};

}//end namespace lar_pandora

#endif
