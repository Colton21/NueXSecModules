#ifndef TPCOBJECTCONTAINER_H
#define TPCOBJECTCONTAINER_H

#include "AnaHelper.h"
#include "ParticleContainer.h"

namespace xsec_ana {

class pfp_container {

private:

std::vector<ParticleContainer> fParticleList;

int run;
int subrun;
int event;
int index;

int n_pfp;
int n_pfp_nu;

int mode;
int origin;
int is_cc;

int mc_vtx;
int mc_vtx_x;
int mc_vtx_y;
int mc_vtx_z;
int pfp_vtx_x;
int pfp_vtx_y;
int pfp_vtx_z;

double completeness;
double purity;

int n_mc_hits;
int n_mc_hits_u;
int n_mc_hits_v;
int n_mc_hits_w;
int n_pfp_hits;
int n_pfp_hits_u;
int n_pfp_hits_v;
int n_pfp_hits_w;

int n_matched_hits;
int n_matched_hits_u;
int n_matched_hits_v;
int n_matched_hits_w;

public:

pfp_container();

//destructor
virtual ~pfp_container();

void AddParticle ( ParticleContainer & particle);
const ParticleContainer & GetParticle (int i) const;

//----------------
//Setter Functions
//----------------

void SetRunNumber (int);
void SetSubRunNumber (int);
void SetEventNumber(int);   //event
void SetIndex (int); //tpc object num per event
//void SetNuMMCParticles (int);
//nMCNeutrinos;
void SetNumPFParticles (int);//number of pfp in the tpc object
void SetNumPFPNeutrinos (int); //number of pfp neutrinos in the tpc object
void SetMode (int);//interaction mode
//0 = CCQE-like
//1 = Resonant
//2 = Deep Inelastic Scattering
//3 = Coherent Scattering
void SetOrigin (int);//this is the origin coming from the tpcobject
void SetIsCC (int); //true = CC, false = NC

//the tpc object takes the neutrino vertex
void SetmcVtxX  (double);
void SetmcVtxY  (double);
void SetmcVtxZ  (double);
void SetpfpVtxX (double);
void SetpfpVtxY (double);
void SetpfpVtxZ (double);

void SetCompleteness (double);
//completeness is defined as:
void SetPurity (double);
//Purity is defined as:

//Hits are summed from the individual particles in the object
void SetNumMCHits   (int);
void SetNumMCHitsU  (int);
void SetNumMCHitsV  (int);
void SetNumMCHitsW  (int);
void SetNumPFPHits  (int);
void SetNumPFPHitsU (int);
void SetNumPFPHitsV (int);
void SetNumPFPHitsW (int);

void SetNumMatchedHits  (int);
void SetNumMatchedHitsU (int);
void SetNumMatchedHitsV (int);
void SetNumMatchedHitsW (int);

//----------------
//Getter Functions
//----------------

const int RunNumber ();
const int SubRunNumber ();
const int EventNumber();
const int Index ();

const int NumPFParticles ();
const int NumPFPNeutrinos ();
const int Mode ();
const int Origin ();
const int IsCC ();

const double mcVtxX  ();
const double mcVtxY  ();
const double mcVtxZ  ();
const double pfpVtxX ();
const double pfpVtxY ();
const double pfpVtxZ ();

const double Completeness ();
const double Purity ();

const int NumMCHits   ();
const int NumMCHitsU  ();
const int NumMCHitsV  ();
const int NumMCHitsW  ();
const int NumPFPHits  ();
const int NumPFPHitsU ();
const int NumPFPHitsV ();
const int NumPFPHitsW ();

const int NumMatchedHits  ();
const int NumMatchedHitsU ();
const int NumMatchedHitsV ();
const int NumMatchedHitsW ();

};

}//end namespace

#endif
