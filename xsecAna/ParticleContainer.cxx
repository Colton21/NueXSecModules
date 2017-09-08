#include "ParticleContainer.h"

namespace xsec_ana {

void SetIndex(int);

void SetOrigin(int _origin){
	origin = _origin;
}
void SetmcPdgCode(int _mc_pdg){
	mc_pdg = _mc_pdg;
}
void SetmcNuPdgCode(int _mc_nu_pdg){
	mc_nu_pdg = _mc_nu_pdg
}
void SetmcParentPdgCode(int _mc_parent_pdg){
	mc_parent_pdg = _mc_parent_pdg;
}
void SetIsNeutrino(bool _is_neutrino){
	is_neutrino = _is_neutrino;
}
void SetmcIsPrimary(bool _is_primary){
	is_primary = _is_primary;
}
void SetpfpPdgCode(int _pfp_pdg){
	pfp_pdg = _pfp_pdg;
}
void SetpfpNuPdgCode(int _pfp_nu_pdg){
	pfp_nu_pdg = _pfp_nu_pdg;
}
void SetpfpParentPdgCode(int _pfp_parent_pdg){
	pfp_parent_pdg = _pfp_parent_pdg;
}

void SetmcVtxX  (double _mc_vtx_x){
	mc_vtx_x = _mc_vtx_x;
}
void SetmcVtxY  (double _mc_vtx_y){
	mc_vtx_y = _mc_vtx_y;
}
void SetmcVtxZ  (double _mc_vtx_z){
	mc_vtx_z = _mc_vtx_z;
}
void SetpfpVtxX (double _pfp_vtx_x){
	pfp_vtx_x = _pfp_vtx_x;
}
void SetpfpVtxY (double _pfp_vtx_y){
	pfp_vtx_y = _pfp_vtx_y;
}
void SetpfpVtxZ (double _pfp_vtx_z){
	pfp_vtx_z = _pfp_vtx_z;
}

void SetmcDirX  (double _mc_dir_x){
	mc_dir_x = _mc_dir_x;
}
void SetmcDirY  (double _mc_dir_y){
	mc_dir_y = _mc_dir_y;
}
void SetmcDirZ  (double _mc_dir_z){
	mc_dir_z = _mc_dir_z;
}
void SetpfpDirX (double _pfp_dir_x){
	pfp_dir_x = _pfp_dir_x;
}
void SetpfpDirY (double _pfp_dir_y){
	pfp_dir_y = _pfp_dir_y;
}
void SetpfpDirZ (double _pfp_dir_z){
	pfp_dir_z = _pfp_dir_z;
}

void SetmcTheta (double _mc_theta){
	mc_theta = _mc_theta;
}
void SetmcPhi (double _mc_phi){
	mc_phi = _mc_phi;
}
void SetpfpTheta (double _pfp_theta){
	pfp_theta = _pfp_theta;
}
void SetpfpPhi (double _pfp_phi){
	pfp_phi = _pfp_phi;
}

void SetmcLength(double _mc_length){
	mc_length = _mc_length;
}
void SetpfpLength(double _pfp_length){
	pfp_length = _pfp_length;
}

void SetmcEnergy(double _mc_energy){
	mc_energy = _mc_energy;
}
void SetmcMomentum(double _mc_momentum){
	mc_momentum = _mc_momentum;
}
void SetpfpMomentum(double _pfp_momentum){
	pfp_momentum = _pfp_momentum;
}

void SetCompleteness(double _completeness){
	completeness = _completeness;
}
void SetPurity(double _purity){
	purity = _purity;
}

void SetNumMCHits   (int _n_mc_hits){
	n_mc_hits = _n_mc_hits;
}
void SetNumMCHitsU  (int _n_mc_hits_u){
	n_mc_hits_u = _n_mc_hits_u;
}
void SetNumMCHitsV  (int _n_mc_hits_v){
	n_mc_htis_v = _n_mc_hits_v;
}
void SetNumMCHitsW  (int _n_mc_hits_w){
	n_mc_hits_w = _n_mc_hits_w;
}
void SetNumPFPHits  (int _n_pfp_hits){
	n_pfp_hits = _n_pfp_hits;
}
void SetNumPFPHitsU (int _n_pfp_hits_u){
	n_pfp_hits_u = _n_pfp_hits_u;
}
void SetNumPFPHitsV (int _n_pfp_hits_v){
	n_pfp_hits_v = _n_pfp_hits_v;
}
void SetNumPFPHitsW (int _n_pfp_hits_w){
	n_pfp_hits_w = _n_pfp_hits_w;
}

void SetNumMatchedHits  (int _n_matched_hits){
	n_matched_hits = _n_matched_hits;
}
void SetNumMatchedHitsU (int _n_matched_hits_u){
	n_matched_hits_u = _n_matched_hits_u;
}
void SetNumMatchedHitsV (int _n_matched_hits_v){
	n_matched_hits_v = _n_matched_hits_v;
}
void SetNumMatchedHitsW (int _n_matched_hits_w){
	n_matched_hits_w = _n_matched_hits_w;
}

void SetmcOpenAngle(double _mc_open_angle){
	mc_open_angle = _mc_open_angle;
}
void SetpfpOpenAngle(double _pfp_open_angle){
	pfp_open_angle = _pfp_open_angle;
}

//need to write the getter functions too!

const int Index(){
	return index;
} //index is particle number in tpc object container

const int Origin(){
	return origin;
}  //this is the particle origin

const int MCPdgCode(){
	return mc_pdg;
}  //true pdg code for the pfp
const int MCNuPdgCode(){
	return mc_nu_pdg;
}  //true nu pdg code (parent) for the pfp
const int MCParentPdg(){
	return mc_parent_pdg;
}  //true parent pdg of the pfp
const bool IsNeutirno(){
	return is_neutirno;
}  //is the pfparticle 12 / 14 pdg code
const bool IsPrimary(){
	return is_primary;
}
const int PFParticlePdgCode(){
	return pfp_pdg;
}
const int PFParticleNuPdgCode(){
	return pfp_nu_pdg;
}
const int PFParticleParentPdgCode(){
	return pfp_parent_pdg;
}

const double mcVtxX(){
	return mc_vtx_x;
}
const double mcVtxY(){
	return mc_vtx_y;
}
const double mcVtxZ(){
	return mc_vtx_z;
}
const double pfpVtxX(){
	return pfp_vtx_x;
}
const double pfpVtxY(){
	return pfp_vtx_y;
}
const double pfpVtxZ(){
	return pfp_vtx_z;
}

const double mcDirX(){
	return mc_dir_x;
}
const double mcDirY(){
	return mc_dir_y;
}
const double mcDirZ(){
	return mc_dir_z;
}
const double pfpDirX(){
	return pfp_dir_x;
}
const double pfpDirY(){
	return pfp_dir_y;
}
const double pfpDirZ(){
	return pfp_dir_z;
}

const double mcTheta(){
	return mc_theta;
}
const double mcPhi(){
	return mc_phi;
}
const double pfpTheta(){
	return pfp_theta;
}
const double pfpPhi(){
	return pfp_phi;
}

const double mcLength(){
	return mc_length;
}
const double pfpLength(){
	return pfp_length;
}

const double mcEnergy(){
	return mc_energy;
}
const double mcMomentum(){
	return mc_momentum;
}
const double pfpMomentum(){
	return pfp_momentum;
}

const double completeness(){
	return completeness;
}
const double purity(){
	return purity;
}

const int NumMCHits   (){
	return n_mc_hits;
}
const int NumMCHitsU  (){
	return n_mc_hits_u;
}
const int NumMCHitsV  (){
	return n_mc_hits_v;
}
const int NumMCHitsW  (){
	return n_mc_hits_w;
}
const int NumPFPHits  (){
	return n_pfp_hits;
}
const int NumPFPHitsU (){
	return n_pfp_hits_u;
}
const int NumPFPHitsV (){
	return n_pfp_hits_v;
}
const int NumPFPHitsW (){
	return n_pfp_hits_w;
}

const int NumMatchedHits  (){
	return n_matched_hits;
}
const int NumMatchedHitsU (){
	return n_matched_hits_u;
}
const int NumMatchedHitsV (){
	return n_matched_hits_v;
}
const int NumMatchedHitsW (){
	return n_matched_hits_w;
}

const double mcOpenAngle(){
	return mc_open_angle;
}
const double pfpOpenAngle(){
	return pfp_open_angle;
}


}//end namespace
