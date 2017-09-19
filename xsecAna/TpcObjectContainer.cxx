#include "TpcObjectContainer.h"

namespace xsecAna {

TPCObjectContainer::TPCObjectContainer() {
}

TPCObjectContainer::~TPCObjectContainer(){
}

//

void TPCObjectContainer::AddParticle ( ParticleContainer _particle){
	fParticleList.push_back(_particle);
}

const ParticleContainer TPCObjectContainer::GetParticle(int _i) const {
	return fParticleList.at(_i);
}

//

void TPCObjectContainer::SetRunNumber (int run_number){
	run = run_number;
}
void TPCObjectContainer::SetSubRunNumber (int subrun_number){
	subrun = subrun_number;
}
void TPCObjectContainer::SetEventNumber(int event_number){
	event = event_number;
}
void TPCObjectContainer::SetIndex (int _index){
	index = _index;
}
void TPCObjectContainer::SetNumPFParticles (int _n_pfp){
	n_pfp = _n_pfp;
}
void TPCObjectContainer::SetNumPFPNeutrinos (int _n_pfp_nu){
	n_pfp_nu = _n_pfp_nu;
}
void TPCObjectContainer::SetMode (int _mode){
	mode = _mode;
}
void TPCObjectContainer::SetOrigin (std::string _origin){
	origin = _origin;
}
void TPCObjectContainer::SetIsCC (int _is_cc){
	is_cc = _is_cc;
}

//the tpc object takes the neutrino vertex
void TPCObjectContainer::SetmcVtxX  (double _mc_vtx_x){
	mc_vtx_x = _mc_vtx_x;
}
void TPCObjectContainer::SetmcVtxY  (double _mc_vtx_y){
	mc_vtx_y = _mc_vtx_y;
}
void TPCObjectContainer::SetmcVtxZ  (double _mc_vtx_z){
	mc_vtx_z = _mc_vtx_z;
}
void TPCObjectContainer::SetpfpVtxX (double _pfp_vtx_x){
	pfp_vtx_x = _pfp_vtx_x;
}
void TPCObjectContainer::SetpfpVtxY (double _pfp_vtx_y){
	pfp_vtx_y = _pfp_vtx_y;
}
void TPCObjectContainer::SetpfpVtxZ (double _pfp_vtx_z){
	pfp_vtx_z = _pfp_vtx_z;
}

void TPCObjectContainer::SetCompleteness (double _completeness){
	completeness = _completeness;
}
void TPCObjectContainer::SetPurity (double _purity){
	purity = _purity;
}

void TPCObjectContainer::SetNumMCHits   (int _n_mc_hits){
	n_mc_hits = _n_mc_hits;
}
void TPCObjectContainer::SetNumMCHitsU  (int _n_mc_hits_u){
	n_mc_hits_u = _n_mc_hits_u;
}
void TPCObjectContainer::SetNumMCHitsV  (int _n_mc_hits_v){
	n_mc_hits_v = _n_mc_hits_v;
}
void TPCObjectContainer::SetNumMCHitsW  (int _n_mc_hits_w){
	n_mc_hits_w = _n_mc_hits_w;
}
void TPCObjectContainer::SetNumPFPHits  (int _n_pfp_hits){
	n_pfp_hits = _n_pfp_hits;
}
void TPCObjectContainer::SetNumPFPHitsU (int _n_pfp_hits_u){
	n_pfp_hits_u = _n_pfp_hits_u;
}
void TPCObjectContainer::SetNumPFPHitsV (int _n_pfp_hits_v){
	n_pfp_hits_v = _n_pfp_hits_v;
}
void TPCObjectContainer::SetNumPFPHitsW (int _n_pfp_hits_w){
	n_pfp_hits_w = _n_pfp_hits_w;
}

void TPCObjectContainer::SetNumMatchedHits  (int _n_matched_hits){
	n_matched_hits = _n_matched_hits;
}
void TPCObjectContainer::SetNumMatchedHitsU (int _n_matched_hits_u){
	n_matched_hits_u = _n_matched_hits_u;
}
void TPCObjectContainer::SetNumMatchedHitsV (int _n_matched_hits_v){
	n_matched_hits_v = _n_matched_hits_v;
}
void TPCObjectContainer::SetNumMatchedHitsW (int _n_matched_hits_w){
	n_matched_hits_w = _n_matched_hits_w;
}

//----------------
//Getter Functions
//----------------

const int TPCObjectContainer::RunNumber () {
	return run;
}
const int TPCObjectContainer::SubRunNumber () {
	return subrun;
}
const int TPCObjectContainer::EventNumber() {
	return event;
}
const int TPCObjectContainer::Index () {
	return index;
}

const int TPCObjectContainer::NumPFParticles (){
	return n_pfp;
}
const int TPCObjectContainer::NumPFPNeutrinos (){
	return n_pfp_nu;
}
const int TPCObjectContainer::Mode (){
	return mode;
}
const std::string TPCObjectContainer::Origin (){
	return origin;
}
const int TPCObjectContainer::IsCC (){
	return is_cc;
}

const double TPCObjectContainer::mcVtxX  (){
	return mc_vtx_x;
}
const double TPCObjectContainer::mcVtxY  (){
	return mc_vtx_y;
}
const double TPCObjectContainer::mcVtxZ  (){
	return mc_vtx_z;
}
const double TPCObjectContainer::pfpVtxX (){
	return pfp_vtx_x;
}
const double TPCObjectContainer::pfpVtxY (){
	return pfp_vtx_y;
}
const double TPCObjectContainer::pfpVtxZ (){
	return pfp_vtx_z;
}

const double TPCObjectContainer::Completeness (){
	return completeness;
}
const double TPCObjectContainer::Purity (){
	return purity;
}

const int TPCObjectContainer::NumMCHits   () {
	return n_mc_hits;
}
const int TPCObjectContainer::NumMCHitsU  () {
	return n_mc_hits_u;
}
const int TPCObjectContainer::NumMCHitsV  () {
	return n_mc_hits_v;
}
const int TPCObjectContainer::NumMCHitsW  () {
	return n_mc_hits_w;
}
const int TPCObjectContainer::NumPFPHits  () {
	return n_pfp_hits;
}
const int TPCObjectContainer::NumPFPHitsU () {
	return n_pfp_hits_u;
}
const int TPCObjectContainer::NumPFPHitsV () {
	return n_pfp_hits_v;
}
const int TPCObjectContainer::NumPFPHitsW () {
	return n_pfp_hits_w;
}

const int TPCObjectContainer::NumMatchedHits  () {
	return n_matched_hits;
}
const int TPCObjectContainer::NumMatchedHitsU () {
	return n_matched_hits_u;
}
const int TPCObjectContainer::NumMatchedHitsV () {
	return n_matched_hits_v;
}
const int TPCObjectContainer::NumMatchedHitsW () {
	return n_matched_hits_w;
}



}//end namespace
