#ifndef TPCOBJECTCONTAINER_H
#define TPCOBJECTCONTAINER_H

//#include "AnaHelper.h"
#include "ParticleContainer.h"

#include <string>
#include <vector>

namespace xsecAna {

class TPCObjectContainer {

private:

std::vector<ParticleContainer> fParticleList;
//ParticleContainer leading_particle;

int is_data;
int run;
int subrun;
int event;
int index;

int n_pfp;
int n_pfp_nu;

int mode;
std::string origin;
int is_cc;
int mc_pdg_code;
int pfp_pdg_code;

int n_pfp_tracks;
int n_pfp_showers;

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

TPCObjectContainer()=default;

void AddParticle (ParticleContainer);
ParticleContainer GetParticle (int) const;

//----------------
//Setter Functions
//----------------

void SetIsData (int);//0 = false, 1 = true
void SetRunNumber (int);
void SetSubRunNumber (int);
void SetEventNumber(int);   //event
void SetIndex (int); //tpc object num per event
//void SetNuMMCParticles (int);
//nMCNeutrinos;
void SetNumPFParticles (int);//number of pfp in the tpc object
void SetNumPFPNeutrinos (int); //number of pfp neutrinos in the tpc object
void SetMode (int);//interaction mode
void SetmcPdgCode(int);
void SetpfpPdgCode(int);
//0 = CCQE-like
//1 = Resonant
//2 = Deep Inelastic Scattering
//3 = Coherent Scattering
void SetOrigin (std::string);//this is the origin coming from the tpcobject
void SetCCNC (int); //true = CC, false = NC

void SetNPfpTracks  (int);
void SetNPfpShowers (int);
//void SetLeadingParticle (ParticleContainer);

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

int IsData        () const;
int RunNumber     () const;
int SubRunNumber  () const;
int EventNumber   () const;
int Index         () const;

int NumPFParticles  () const;
int NumPFPNeutrinos () const;
int Mode            () const;
std::string Origin  () const;
int CCNC            () const;

int PFParticlePdgCode () const;
int MCParticlePdgCode () const;

int NPfpTracks  () const;
int NPfpShowers () const;
//ParticleContainer LeadingParticle() const;

double mcVtxX  () const;
double mcVtxY  () const;
double mcVtxZ  () const;
double pfpVtxX () const;
double pfpVtxY () const;
double pfpVtxZ () const;

double Completeness () const;
double Purity       () const;

int NumMCHits   () const;
int NumMCHitsU  () const;
int NumMCHitsV  () const;
int NumMCHitsW  () const;
int NumPFPHits  () const;
int NumPFPHitsU () const;
int NumPFPHitsV () const;
int NumPFPHitsW () const;

int NumMatchedHits  () const;
int NumMatchedHitsU () const;
int NumMatchedHitsV () const;
int NumMatchedHitsW () const;

};

}//end namespace

#endif
