#include "ParticleContainer.h"

namespace xsecAna {

// ParticleContainer::ParticleContainer(){
//
// }

void ParticleContainer::SetIndex(int _index){
	index = _index;
}
void ParticleContainer::SetOrigin(std::string _origin){
	origin = _origin;
}
void ParticleContainer::SetmcPdgCode(int _mc_pdg){
	mc_pdg = _mc_pdg;
}
void ParticleContainer::SetmcNuPdgCode(int _mc_nu_pdg){
	mc_nu_pdg = _mc_nu_pdg;
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
void ParticleContainer::SetMode(int _mode){
	mode = _mode;
}
void ParticleContainer::SetCCNC(int _is_cc){
	is_cc = _is_cc;
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
void ParticleContainer::SetmcNeutrinoEnergy(double _mc_neutrino_energy){
	mc_neutrino_energy = _mc_neutrino_energy;
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
	n_mc_hits_v = _n_mc_hits_v;
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
void ParticleContainer::SetpfpEnergyU (int _pfp_energy_u){
	pfp_energy_u = _pfp_energy_u;
}
void ParticleContainer::SetpfpEnergyV (int _pfp_energy_v){
	pfp_energy_v = _pfp_energy_v;
}
void ParticleContainer::SetpfpEnergyW (int _pfp_energy_w){
	pfp_energy_w = _pfp_energy_w;
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

void ParticleContainer::SetPfpClusterdQdx(std::vector<std::vector<double> > _pfp_cluster_dqdx){
	pfp_cluster_dqdx = _pfp_cluster_dqdx;
}

void ParticleContainer::SetPfpdEdx(std::vector<double> _pfp_dedx){
	pfp_dedx = _pfp_dedx;
}

//need to write the getter functions too!

int ParticleContainer::Index() const {
	return index;
} //index is particle number in tpc object container

std::string ParticleContainer::Origin() const {
	return origin;
}  //this is the particle origin

int ParticleContainer::MCPdgCode() const {
	return mc_pdg;
}  //true pdg code for the pfp
int ParticleContainer::MCNuPdgCode() const {
	return mc_nu_pdg;
}  //true nu pdg code (parent) for the pfp
int ParticleContainer::MCParentPdg() const {
	return mc_parent_pdg;
}  //true parent pdg of the pfp
bool ParticleContainer::IsNeutrino() const {
	return is_neutrino;
}  //is the pfparticle 12 / 14 pdg code
bool ParticleContainer::IsPrimary() const {
	return is_primary;
}
int ParticleContainer::PFParticlePdgCode() const {
	return pfp_pdg;
}
int ParticleContainer::PFParticleNuPdgCode() const {
	return pfp_nu_pdg;
}
int ParticleContainer::PFParticleParentPdgCode() const {
	return pfp_parent_pdg;
}
int ParticleContainer::Mode() const {
	return mode;
}
int ParticleContainer::CCNC() const {
	return is_cc;
}

double ParticleContainer::mcVtxX() const {
	return mc_vtx_x;
}
double ParticleContainer::mcVtxY() const {
	return mc_vtx_y;
}
double ParticleContainer::mcVtxZ() const {
	return mc_vtx_z;
}
double ParticleContainer::pfpVtxX() const {
	return pfp_vtx_x;
}
double ParticleContainer::pfpVtxY() const {
	return pfp_vtx_y;
}
double ParticleContainer::pfpVtxZ() const {
	return pfp_vtx_z;
}

double ParticleContainer::mcDirX() const {
	return mc_dir_x;
}
double ParticleContainer::mcDirY() const {
	return mc_dir_y;
}
double ParticleContainer::mcDirZ() const {
	return mc_dir_z;
}
double ParticleContainer::pfpDirX() const {
	return pfp_dir_x;
}
double ParticleContainer::pfpDirY() const {
	return pfp_dir_y;
}
double ParticleContainer::pfpDirZ() const {
	return pfp_dir_z;
}

double ParticleContainer::mcTheta() const {
	return mc_theta;
}
double ParticleContainer::mcPhi() const {
	return mc_phi;
}
double ParticleContainer::pfpTheta() const {
	return pfp_theta;
}
double ParticleContainer::pfpPhi() const {
	return pfp_phi;
}

double ParticleContainer::mcLength() const {
	return mc_length;
}
double ParticleContainer::pfpLength() const {
	return pfp_length;
}

double ParticleContainer::mcEnergy() const {
	return mc_energy;
}
double ParticleContainer::mcMomentum() const {
	return mc_momentum;
}
double ParticleContainer::pfpMomentum() const {
	return pfp_momentum;
}
double ParticleContainer::mcNeutrinoEnergy() const {
	return mc_neutrino_energy;
}

double ParticleContainer::Completeness() const {
	return completeness;
}
double ParticleContainer::Purity() const {
	return purity;
}

int ParticleContainer::NumMCHits   () const {
	return n_mc_hits;
}
int ParticleContainer::NumMCHitsU  () const {
	return n_mc_hits_u;
}
int ParticleContainer::NumMCHitsV  () const {
	return n_mc_hits_v;
}
int ParticleContainer::NumMCHitsW  () const {
	return n_mc_hits_w;
}
int ParticleContainer::NumPFPHits  () const {
	return n_pfp_hits;
}
int ParticleContainer::NumPFPHitsU () const {
	return n_pfp_hits_u;
}
int ParticleContainer::NumPFPHitsV () const {
	return n_pfp_hits_v;
}
int ParticleContainer::NumPFPHitsW () const {
	return n_pfp_hits_w;
}
double ParticleContainer::PfpEnergyU () const {
	return pfp_energy_u;
}
double ParticleContainer::PfpEnergyV () const {
	return pfp_energy_v;
}
double ParticleContainer::PfpEnergyW () const {
	return pfp_energy_w;
}

int ParticleContainer::NumMatchedHits  () const {
	return n_matched_hits;
}
int ParticleContainer::NumMatchedHitsU () const {
	return n_matched_hits_u;
}
int ParticleContainer::NumMatchedHitsV () const {
	return n_matched_hits_v;
}
int ParticleContainer::NumMatchedHitsW () const {
	return n_matched_hits_w;
}

double ParticleContainer::mcOpenAngle() const {
	return mc_open_angle;
}
double ParticleContainer::pfpOpenAngle() const {
	return pfp_open_angle;
}
std::vector < std::vector < double > > ParticleContainer::PfpClusterdQdx() const {
	return pfp_cluster_dqdx;
}

std::vector<double > ParticleContainer::PfpdEdx() const {
	return pfp_dedx;
}

}//end namespace
