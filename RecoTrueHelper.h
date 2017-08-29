#include "AnaHelper.h"

namespace nue_xsec {

class recotruehelper {

private:

typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
typedef std::set< art::Ptr<simb::MCParticle> > MCParticleSet;
typedef std::set< art::Ptr<simb::MCTruth> > MCTruthSet;

//I should export this as a fcl parameter in the future
bool _recursiveMatching = true;


/**
 *  @brief  Build mapping from true neutrinos to hits
 *
 *  @param truthToParticles  the input mapping from true event to true particles
 *  @param trueParticlesToHits  the input mapping from true particles to hits
 *  @param trueNeutrinosToHits  the output mapping from trues event to hits
 *  @param trueHitsToNeutrinos  the output mappign from hits to true events
 */
void BuildTrueNeutrinoHitMaps(const MCTruthToMCParticles &truthToParticles, const MCParticlesToHits &trueParticlesToHits,
                              MCTruthToHits &trueNeutrinosToHits, HitsToMCTruth &trueHitsToNeutrinos) const;

/**
 *  @brief  Build mapping from reconstructed neutrinos to hits
 *
 *  @param recoParticleMap  the input mapping from reconstructed particle and particle ID
 *  @param recoParticlesToHits  the input mapping from reconstructed particles to hits
 *  @param recoNeutrinosToHits  the output mapping from reconstructed particles to hits
 *  @param recoHitsToNeutrinos  the output mapping from reconstructed hits to particles
 */
void BuildRecoNeutrinoHitMaps(const PFParticleMap &recoParticleMap, const PFParticlesToHits &recoParticlesToHits,
                              PFParticlesToHits &recoNeutrinosToHits, HitsToPFParticles &recoHitsToNeutrinos) const;

/**
 *  @brief Perform matching between true and reconstructed neutrino events
 *
 *  @param recoNeutrinosToHits  the mapping from reconstructed neutrino events to hits
 *  @param trueHitsToNeutrinos  the mapping from hits to true neutrino events
 *  @param matchedNeutrinos  the output matches between reconstructed and true neutrinos
 *  @param matchedNeutrinoHits  the output matches between reconstructed neutrinos and hits
 */
void GetRecoToTrueMatches(const PFParticlesToHits &recoNeutrinosToHits, const HitsToMCTruth &trueHitsToNeutrinos,
                          MCTruthToPFParticles &matchedNeutrinos, MCTruthToHits &matchedNeutrinoHits) const;

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
void GetRecoToTrueMatches(const PFParticlesToHits &recoNeutrinosToHits, const HitsToMCTruth &trueHitsToNeutrinos,
                          MCTruthToPFParticles &matchedNeutrinos, MCTruthToHits &matchedNeutrinoHits, PFParticleSet &recoVeto, MCTruthSet &trueVeto) const;

/**
 *  @brief Perform matching between true and reconstructed particles
 *
 *  @param recoParticlesToHits the mapping from reconstructed particles to hits
 *  @param trueHitsToParticles the mapping from hits to true particles
 *  @param matchedParticles the output matches between reconstructed and true particles
 *  @param matchedHits the output matches between reconstructed particles and hits
 */
void GetRecoToTrueMatches(const PFParticlesToHits &recoParticlesToHits, const HitsToMCParticles &trueHitsToParticles,
                          MCParticlesToPFParticles &matchedParticles, MCParticlesToHits &matchedHits) const;

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
void GetRecoToTrueMatches(const PFParticlesToHits &recoParticlesToHits, const HitsToMCParticles &trueHitsToParticles,
                          MCParticlesToPFParticles &matchedParticles, MCParticlesToHits &matchedHits, PFParticleSet &recoVeto, MCParticleSet &trueVeto) const;

/**
 *  @brief Build particle maps for reconstructed particles
 *
 *  @param particleVector the input vector of reconstructed particles
 *  @param particleMap the output mapping between reconstructed particles and particle ID
 */
void BuildRecoParticleMap(const PFParticleVector &particleVector, PFParticleMap &particleMap) const;

/**
 *  @brief Build particle maps for true particles
 *
 *  @param particleVector the input vector of true particles
 *  @param particleMap the output mapping between true particle and true track ID
 */
void BuildTrueParticleMap(const MCParticleVector &particleVector, MCParticleMap &particleMap) const;

/**
 *  @brief Count the number of reconstructed hits in a given wire plane
 *
 *  @param view the wire plane ID
 *  @param hitVector the input vector of reconstructed hits
 */
int CountHitsByType(const int view, const HitVector &hitVector) const;

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
