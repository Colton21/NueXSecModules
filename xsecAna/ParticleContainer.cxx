#include "ParticleContainer.h"

namespace xsec_ana {

void ParticleContainer::SetIndex(int);

void ParticleContainer::SetOrigin(int _origin){
	origin = _origin;
}
void ParticleContainer::SetmcPdgCode(int _mc_pdg){
	mc_pdg = _mc_pdg;
}
void ParticleContainer::SetmcNuPdgCode(int _mc_nu_pdg){
	mc_nu_pdg = _mc_nu_pdg
}
void ParticleContainer::SetmcParentPdgCode(int _mc_parent_pdg){
	mc_parent_pdg = _mc_parent_pdg;
}
void ParticleContainer::SetIsNeutrino(bool _is_neutrino){
	is_neutrino = _is_neutrino;
}
void ParticleContainer::SetmcIsPrimary(bool _is_primary){
	is_primary = _is_primary;
}
void ParticleContainer::SetpfpPdgCode(int _pfp_pdg){
	pfp_pdg = _pfp_pdg;
}
void ParticleContainer::SetpfpNuPdgCode(int _pfp_nu_pdg){
	pfp_nu_pdg = _pfp_nu_pdg;
}
void ParticleContainer::SetpfpParentPdgCode(int _pfp_parent_pdg){
	pfp_parent_pdg = _pfp_parent_pdg;
}

void ParticleContainer::SetmcVtxX  (double _mc_vtx_x){
	mc_vtx_x = _mc_vtx_x;
}
void ParticleContainer::SetmcVtxY  (double _mc_vtx_y){
	mc_vtx_y = _mc_vtx_y;
}
void ParticleContainer::SetmcVtxZ  (double _mc_vtx_z){
	mc_vtx_z = _mc_vtx_z;
}
void ParticleContainer::SetpfpVtxX (double _pfp_vtx_x){
	pfp_vtx_x = _pfp_vtx_x;
}
void ParticleContainer::SetpfpVtxY (double _pfp_vtx_y){
	pfp_vtx_y = _pfp_vtx_y;
}
void ParticleContainer::SetpfpVtxZ (double _pfp_vtx_z){
	pfp_vtx_z = _pfp_vtx_z;
}

void ParticleContainer::SetmcDirX  (double _mc_dir_x){
	mc_dir_x = _mc_dir_x;
}
void ParticleContainer::SetmcDirY  (double _mc_dir_y){
	mc_dir_y = _mc_dir_y;
}
void ParticleContainer::SetmcDirZ  (double _mc_dir_z){
	mc_dir_z = _mc_dir_z;
}
void ParticleContainer::SetpfpDirX (double _pfp_dir_x){
	pfp_dir_x = _pfp_dir_x;
}
void ParticleContainer::SetpfpDirY (double _pfp_dir_y){
	pfp_dir_y = _pfp_dir_y;
}
void ParticleContainer::SetpfpDirZ (double _pfp_dir_z){
	pfp_dir_z = _pfp_dir_z;
}

void ParticleContainer::SetmcTheta (double _mc_theta){
	mc_theta = _mc_theta;
}
void ParticleContainer::SetmcPhi (double _mc_phi){
	mc_phi = _mc_phi;
}
void ParticleContainer::SetpfpTheta (double _pfp_theta){
	pfp_theta = _pfp_theta;
}
void ParticleContainer::SetpfpPhi (double _pfp_phi){
	pfp_phi = _pfp_phi;
}

void ParticleContainer::SetmcLength(double _mc_length){
	mc_length = _mc_length;
}
void ParticleContainer::SetpfpLength(double _pfp_length){
	pfp_length = _pfp_length;
}

void ParticleContainer::SetmcEnergy(double _mc_energy){
	mc_energy = _mc_energy;
}
void ParticleContainer::SetmcMomentum(double _mc_momentum){
	mc_momentum = _mc_momentum;
}
void ParticleContainer::SetpfpMomentum(double _pfp_momentum){
	pfp_momentum = _pfp_momentum;
}

void ParticleContainer::SetCompleteness(double _completeness){
	completeness = _completeness;
}
void ParticleContainer::SetPurity(double _purity){
	purity = _purity;
}

void ParticleContainer::SetNumMCHits   (int _n_mc_hits){
	n_mc_hits = _n_mc_hits;
}
void ParticleContainer::SetNumMCHitsU  (int _n_mc_hits_u){
	n_mc_hits_u = _n_mc_hits_u;
}
void ParticleContainer::SetNumMCHitsV  (int _n_mc_hits_v){
	n_mc_htis_v = _n_mc_hits_v;
}
void ParticleContainer::SetNumMCHitsW  (int _n_mc_hits_w){
	n_mc_hits_w = _n_mc_hits_w;
}
void ParticleContainer::SetNumPFPHits  (int _n_pfp_hits){
	n_pfp_hits = _n_pfp_hits;
}
void ParticleContainer::SetNumPFPHitsU (int _n_pfp_hits_u){
	n_pfp_hits_u = _n_pfp_hits_u;
}
void ParticleContainer::SetNumPFPHitsV (int _n_pfp_hits_v){
	n_pfp_hits_v = _n_pfp_hits_v;
}
void ParticleContainer::SetNumPFPHitsW (int _n_pfp_hits_w){
	n_pfp_hits_w = _n_pfp_hits_w;
}

void ParticleContainer::SetNumMatchedHits  (int _n_matched_hits){
	n_matched_hits = _n_matched_hits;
}
void ParticleContainer::SetNumMatchedHitsU (int _n_matched_hits_u){
	n_matched_hits_u = _n_matched_hits_u;
}
void ParticleContainer::SetNumMatchedHitsV (int _n_matched_hits_v){
	n_matched_hits_v = _n_matched_hits_v;
}
void ParticleContainer::SetNumMatchedHitsW (int _n_matched_hits_w){
	n_matched_hits_w = _n_matched_hits_w;
}

void ParticleContainer::SetmcOpenAngle(double _mc_open_angle){
	mc_open_angle = _mc_open_angle;
}
void ParticleContainer::SetpfpOpenAngle(double _pfp_open_angle){
	pfp_open_angle = _pfp_open_angle;
}

//need to write the getter functions too!

const int ParticleContainer::Index(){
	return index;
} //index is particle number in tpc object container

const int ParticleContainer::Origin(){
	return origin;
}  //this is the particle origin

const int ParticleContainer::MCPdgCode(){
	return mc_pdg;
}  //true pdg code for the pfp
const int ParticleContainer::MCNuPdgCode(){
	return mc_nu_pdg;
}  //true nu pdg code (parent) for the pfp
const int ParticleContainer::MCParentPdg(){
	return mc_parent_pdg;
}  //true parent pdg of the pfp
const bool ParticleContainer::IsNeutirno(){
	return is_neutirno;
}  //is the pfparticle 12 / 14 pdg code
const bool ParticleContainer::IsPrimary(){
	return is_primary;
}
const int ParticleContainer::PFParticlePdgCode(){
	return pfp_pdg;
}
const int ParticleContainer::PFParticleNuPdgCode(){
	return pfp_nu_pdg;
}
const int ParticleContainer::PFParticleParentPdgCode(){
	return pfp_parent_pdg;
}

const double ParticleContainer::mcVtxX(){
	return mc_vtx_x;
}
const double ParticleContainer::mcVtxY(){
	return mc_vtx_y;
}
const double ParticleContainer::mcVtxZ(){
	return mc_vtx_z;
}
const double ParticleContainer::pfpVtxX(){
	return pfp_vtx_x;
}
const double ParticleContainer::pfpVtxY(){
	return pfp_vtx_y;
}
const double ParticleContainer::pfpVtxZ(){
	return pfp_vtx_z;
}

const double ParticleContainer::mcDirX(){
	return mc_dir_x;
}
const double ParticleContainer::mcDirY(){
	return mc_dir_y;
}
const double ParticleContainer::mcDirZ(){
	return mc_dir_z;
}
const double ParticleContainer::pfpDirX(){
	return pfp_dir_x;
}
const double ParticleContainer::pfpDirY(){
	return pfp_dir_y;
}
const double ParticleContainer::pfpDirZ(){
	return pfp_dir_z;
}

const double ParticleContainer::mcTheta(){
	return mc_theta;
}
const double ParticleContainer::mcPhi(){
	return mc_phi;
}
const double ParticleContainer::pfpTheta(){
	return pfp_theta;
}
const double ParticleContainer::pfpPhi(){
	return pfp_phi;
}

const double ParticleContainer::mcLength(){
	return mc_length;
}
const double ParticleContainer::pfpLength(){
	return pfp_length;
}

const double ParticleContainer::mcEnergy(){
	return mc_energy;
}
const double ParticleContainer::mcMomentum(){
	return mc_momentum;
}
const double ParticleContainer::pfpMomentum(){
	return pfp_momentum;
}

const double ParticleContainer::completeness(){
	return completeness;
}
const double ParticleContainer::purity(){
	return purity;
}

const int ParticleContainer::NumMCHits   (){
	return n_mc_hits;
}
const int ParticleContainer::NumMCHitsU  (){
	return n_mc_hits_u;
}
const int ParticleContainer::NumMCHitsV  (){
	return n_mc_hits_v;
}
const int ParticleContainer::NumMCHitsW  (){
	return n_mc_hits_w;
}
const int ParticleContainer::NumPFPHits  (){
	return n_pfp_hits;
}
const int ParticleContainer::NumPFPHitsU (){
	return n_pfp_hits_u;
}
const int ParticleContainer::NumPFPHitsV (){
	return n_pfp_hits_v;
}
const int ParticleContainer::NumPFPHitsW (){
	return n_pfp_hits_w;
}

const int ParticleContainer::NumMatchedHits  (){
	return n_matched_hits;
}
const int ParticleContainer::NumMatchedHitsU (){
	return n_matched_hits_u;
}
const int ParticleContainer::NumMatchedHitsV (){
	return n_matched_hits_v;
}
const int ParticleContainer::NumMatchedHitsW (){
	return n_matched_hits_w;
}

const double ParticleContainer::mcOpenAngle(){
	return mc_open_angle;
}
const double ParticleContainer::pfpOpenAngle(){
	return pfp_open_angle;
}


}//end namespace
