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

RecoTrueMatcher matchinghelper;

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
int nMatchedHitsW;

double mcOpenAngle;
double pfpOpenAngle;

};

#endif
