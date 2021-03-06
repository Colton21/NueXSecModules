#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

//#include "AnaHelper.h"

#include <string>
#include <vector>

namespace xsecAna {

class ParticleContainer {

private:

int index;

std::string origin;
int mc_pdg;
int mc_nu_pdg;
int mc_parent_pdg;
bool is_neutrino;
bool is_primary;
int pfp_pdg;
int pfp_nu_pdg;
int pfp_parent_pdg;
int mode;
int is_cc;

int mc_vtx_x;
int mc_vtx_y;
int mc_vtx_z;
int pfp_vtx_x;
int pfp_vtx_y;
int pfp_vtx_z;

double mc_dir_x;
double mc_dir_y;
double mc_dir_z;
double pfp_dir_x;
double pfp_dir_y;
double pfp_dir_z;

double mc_theta;
double mc_phi;
double pfp_theta;
double pfp_phi;

double mc_length;
double pfp_length;

double mc_energy;
double mc_momentum;
double pfp_momentum;
double mc_neutrino_energy;

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
double pfp_energy_u;
double pfp_energy_v;
double pfp_energy_w;

int n_matched_hits;
int n_matched_hits_u;
int n_matched_hits_v;
int n_matched_hits_w;

double mc_open_angle;
double pfp_open_angle;

std::vector < std::vector < double > > pfp_cluster_dqdx;
std::vector < std::vector < double > > pfp_cluster_dq;
std::vector < std::vector < double > > pfp_cluster_dx;
std::vector<double> pfp_dedx;

std::vector < std::vector < double > > pfp_cluster_dqdx_alt;
std::vector < std::vector < double > > pfp_cluster_dq_alt;
std::vector < std::vector < double > > pfp_cluster_dx_alt;
std::vector<double> pfp_dedx_alt;

std::vector< std::vector < double > > pfp_cluster_dqdx_cali;
std::vector< std::vector < double > > pfp_cluster_dqdx_omit;
std::vector< std::vector < double > > pfp_cluster_dqdx_omit_cali;

std::vector<double> pfp_dedx_cali;
std::vector<double> pfp_dedx_omit;
std::vector<double> pfp_dedx_omit_cali;

public:

ParticleContainer()=default;

//-----------------------
//Setter functions
//------------------------

void SetIndex(int);

void SetOrigin(std::string);
void SetmcPdgCode(int);
void SetmcNuPdgCode(int);
void SetmcParentPdgCode(int);
void SetIsNeutrino(bool);
void SetmcIsPrimary(bool);
void SetpfpPdgCode(int);
void SetpfpNuPdgCode(int);
void SetpfpParentPdgCode(int);
void SetMode(int);
void SetCCNC(int);

void SetmcVtxX  (double);
void SetmcVtxY  (double);
void SetmcVtxZ  (double);
void SetpfpVtxX (double);
void SetpfpVtxY (double);
void SetpfpVtxZ (double);

void SetmcDirX  (double);
void SetmcDirY  (double);
void SetmcDirZ  (double);
void SetpfpDirX (double);
void SetpfpDirY (double);
void SetpfpDirZ (double);

void SetmcTheta (double);
void SetmcPhi (double);
void SetpfpTheta (double);
void SetpfpPhi (double);

void SetmcLength(double);
void SetpfpLength(double);

void SetmcEnergy(double);
void SetmcMomentum(double);
void SetpfpMomentum(double);
void SetmcNeutrinoEnergy(double);

void SetCompleteness(double);
void SetPurity(double);

void SetNumMCHits   (int);
void SetNumMCHitsU  (int);
void SetNumMCHitsV  (int);
void SetNumMCHitsW  (int);
void SetNumPFPHits  (int);
void SetNumPFPHitsU (int);
void SetNumPFPHitsV (int);
void SetNumPFPHitsW (int);
void SetpfpEnergyU  (double);
void SetpfpEnergyV  (double);
void SetpfpEnergyW  (double);

void SetNumMatchedHits  (int);
void SetNumMatchedHitsU (int);
void SetNumMatchedHitsV (int);
void SetNumMatchedHitsW (int);

void SetmcOpenAngle(double);
void SetpfpOpenAngle(double);

void SetPfpClusterdQdx(std::vector<std::vector<double> > );
void SetPfpClusterdQ(std::vector<std::vector<double> >);
void SetPfpClusterdX(std::vector<std::vector<double> >);
void SetPfpdEdx(std::vector<double>);

void SetPfpClusterdQdx_alt(std::vector<std::vector<double> > );
void SetPfpClusterdQ_alt(std::vector<std::vector<double> >);
void SetPfpClusterdX_alt(std::vector<std::vector<double> >);
void SetPfpdEdx_alt(std::vector<double>);

void SetPfpClusterdQdxCali(std::vector<std::vector<double>>);
void SetPfpClusterdQdxOmitFirst(std::vector<std::vector<double>>);
void SetPfpClusterdQdxOmitFirst_cali(std::vector<std::vector<double>>);

void SetPfpdEdxCali(std::vector<double>);
void SetPfpdEdxOmitFirst(std::vector<double>);
void SetPfpdEdxOmitFirst_cali(std::vector<double>);

//need to write the getter functions too!

int Index() const; //index is particle number in tpc object container

std::string Origin() const; //this is the particle origin

int MCPdgCode           () const;//true pdg code for the pfp
int MCNuPdgCode         () const;//true nu pdg code (parent) for the pfp
int MCParentPdg         () const;//true parent pdg of the pfp
bool IsNeutrino         () const;//is the pfparticle 12 / 14 pdg code
bool IsPrimary          () const;
int PFParticlePdgCode   () const;
int PFParticleNuPdgCode () const;
//const int PFParticleNuIndex();
int PFParticleParentPdgCode() const;
int Mode () const;
int CCNC () const;//0 is cc, 1 is nc

double mcVtxX () const;
double mcVtxY () const;
double mcVtxZ () const;
double pfpVtxX() const;
double pfpVtxY() const;
double pfpVtxZ() const;

double mcDirX   () const;
double mcDirY   () const;
double mcDirZ   () const;
double pfpDirX  () const;
double pfpDirY  () const;
double pfpDirZ  () const;

double mcTheta  () const;
double mcPhi    () const;
double pfpTheta () const;
double pfpPhi   () const;

double mcLength   () const;
double pfpLength  () const;

double mcEnergy         () const;
double mcMomentum       () const;
double pfpMomentum      () const;
double mcNeutrinoEnergy () const;

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
double PfpEnergyU () const;
double PfpEnergyV () const;
double PfpEnergyW () const;

int NumMatchedHits  () const;
int NumMatchedHitsU () const;
int NumMatchedHitsV () const;
int NumMatchedHitsW () const;

double mcOpenAngle  () const;
double pfpOpenAngle () const;

std::vector< std::vector< double > > PfpClusterdQdx () const;
std::vector<std::vector<double> > PfpClusterdQ () const;
std::vector<std::vector<double> > PfpClusterdX () const;
std::vector< double > PfpdEdx () const;

std::vector< std::vector< double > > PfpClusterdQdx_alt () const;
std::vector<std::vector<double> > PfpClusterdQ_alt () const;
std::vector<std::vector<double> > PfpClusterdX_alt () const;
std::vector< double > PfpdEdx_alt () const;

std::vector < std::vector < double > >  PfpClusterdQdxCali() const;
std::vector < std::vector < double > >  PfpClusterdQdxOmitFirst() const;
std::vector < std::vector < double > >  PfpClusterdQdxOmitFirst_cali() const;

std::vector < double > PfpdEdx_cali() const;
std::vector < double > PfpdEdxOmitFirst() const;
std::vector < double > PfpdEdxOmitFirst_cali() const;

};//end class

}//end namespace

#endif
