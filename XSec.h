#include "AnaHelper.h"
#include "PostMCCorrections.h"
#include "RecoTrueHelper.h"

#ifndef XSec_H
#define XSec_H

class XSec;

class XSec : public art::EDAnalyzer {
public:
explicit XSec(fhicl::ParameterSet const & p);
// The compiler-generated destructor is fine for non-base
// classes without bare pointers or other resource use.

// Plugins should not be copied or assigned.
XSec(XSec const &) = delete;
XSec(XSec &&) = delete;
XSec & operator = (XSec const &) = delete;
XSec & operator = (XSec &&) = delete;

// Required functions.
void analyze(art::Event const & e) override;
void endSubRun(const art::SubRun &sr) override;

private:

TTree * myTree;
bool isMC;
bool isData;

//RecoTrueMatcher matchinghelper;

std::string _pfp_producer;
std::string _hitfinderLabe;
std::string _geantModuleLabel;
std::string _spacepointLabel;
std::string _neutrino_flash_match_producer;
std::string _cosmic_flash_match_producer;
std::string _opflash_producer_beam;
std::string _acpt_producer;
std::string _tpcobject_producer;
std::string _potsum_producer;
std::string _potsum_instance;
std::string _particle_id_producer;
std::string _mc_ghost_producer;

std::string _useDaughterPFParticles;
std::string _addDaughterPFParticles;

bool _use_genie_info;
int _minimumHitRequirement;

double _beam_spill_start;
double _beam_spill_end;

bool _debug;
bool _verbose;


int run;
int event;
int index;
int nMCParticles;
int nMCNeutrinos;
int nPFPartilcles;
int nPFPNeutrinos;
int mcPdg;
int mcNuPdg;
int mcNuIndex;
int mcParentPdg;
bool mcIsNeutirno;
bool mcIsPrimary;
int mcMode;
bool mcIsCC;
int pfpPdg;
int pfpNuPdg;
int pfpNuIndex;
int pfpParentPdg;
bool pfpIsNeutrino;

double mcVtxX;
double mcVtxY;
double mcVtxZ;
double pfpVtxX;
double pfpVtxY;
double pfpVtxZ;

double mcDirX;
double mcDirY;
double mcDirZ;
double pfpDirX;
double pfpDirY;
double pfpDirZ;

double mcTheta;
double mcPhi;
double pfpTheta;
double pfpPhi;

double mcLength;
double pfpLength;

double mcEnergy;
double mcMomentum;
//double pfpEnergy;
double pfpMomentum;

double completeness;
double purity;

int nMCHits;
int nMCHitsU;
int nMCHitsV;
int nMCHitsY;
int nPFPHits;
int nPFPHitsU;
int nPFPHitsV;
int nPFPHitsY;
int nMatchedHits;
int nMatchedHitsU;
int nMatchedHitsV;
int nMatchedHitsY;

double mcOpenAngle;
double pfpOpenAngle;

};

#endif
