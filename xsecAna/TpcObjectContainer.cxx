#include "TpcObjectContainer.h"

namespace xsec_ana {

void pfp_container::SetRunNumber (int run_number){
	run = run_number;
}
void pfp_container::SetSubRunNumber (int subrun_number){
	subrun = subrun_number;
}
void pfp_container::SetEventNumber(int event_number){
	event = event_number;
}
void pfp_container::SetIndex (int _index){
	index = _index;
}
void pfp_container::SetNumPFParticles (int _n_pfp){
	n_pfp = _n_pfp;
}
void pfp_container::SetNumPFPNeutrinos (int _n_pfp_nu){
	n_pfp_nu = _n_pfp_nu;
}
void pfp_container::SetMode (int _mode){
	mode = _mode;
}
void pfp_container::SetOrigin (simb::Origin_t _origin){
	origin = _origin;
}
void pfp_container::SetIsCC (int _is_cc){
	is_cc = _is_cc;
}

//the tpc object takes the neutrino vertex
void pfp_container::SetmcVtxX  (double _mc_vtx_x){
	mc_vtx_x = _mc_vtx_x;
}
void pfp_container::SetmcVtxY  (double _mc_vtx_y){
	mc_vtx_y = _mc_vtx_y;
}
void pfp_container::SetmcVtxZ  (double _mc_vtx_z){
	mc_vtx_z = _mc_vtx_z;
}
void pfp_container::SetpfpVtxX (double _pfp_vtx_x){
	pfp_vtx_x = _pfp_vtx_x;
}
void pfp_container::SetpfpVtxY (double _pfp_vtx_y){
	pfp_vtx_y = _pfp_vtx_y;
}
void pfp_container::SetpfpVtxZ (double _pfp_vtx_z){
	pfp_vtx_z = _pfp_vtx_z;
}

void pfp_container::SetCompleteness (double _completeness){
	completeness = _completeness;
}
void pfp_container::SetPurity (double _purity){
	purity = _purity;
}

void pfp_container::SetNumMCHits   (int _n_mc_hits){
	n_mc_hits = _n_mc_hits;
}
void pfp_container::SetNumMCHitsU  (int _n_mc_hits_u){
	n_mc_hits_u = _n_mc_hits_u;
}
void pfp_container::SetNumMCHitsV  (int _n_mc_hits_v){
	n_mc_hits_v = _n_mc_hits_v;
}
void pfp_container::SetNumMCHitsW  (int _n_mc_hits_w){
	n_mc_hits_w = _n_mc_hits_w;
}
void pfp_container::SetNumPFPHits  (int _n_pfp_hits){
	n_pfp_hits = _n_pfp_hits;
}
void pfp_container::SetNumPFPHitsU (int _n_pfp_hits_u){
	n_pfp_hits_u = _n_pfp_hits_u;
}
void pfp_container::SetNumPFPHitsV (int _n_pfp_hits_v){
	n_pfp_hits_v = _n_pfp_hits_v;
}
void pfp_container::SetNumPFPHitsW (int _n_pfp_hits_w){
	n_pfp_hits_w = _n_pfp_hits_w;
}

void pfp_container::SetNumMatchedHits  (int _n_matched_hits){
	n_matched_hits = _n_matched_hits;
}
void pfp_container::SetNumMatchedHitsU (int _n_matched_hits_u){
	n_matched_hits_u = _n_matched_hits_u;
}
void pfp_container::SetNumMatchedHitsV (int _n_matched_hits_v){
	n_matched_hits_v = _n_matched_hits_v;
}
void pfp_container::SetNumMatchedHitsW (int _n_matched_hits_w){
	n_matched_hits_w = _n_matched_hits_w;
}

//----------------
//Getter Functions
//----------------

const int pfp_container::RunNumber () {
	return run;
}
const int pfp_container::SubRunNumber () {
	return subrun;
}
const int pfp_container::EventNumber() {
	return event;
}
const int pfp_container::Index () {
	return index;
}

const int pfp_container::NumPFParticles (){
	return n_pfp;
}
const int pfp_container::NumPFPNeutrinos (){
	return n_pfp_nu;
}
const int pfp_container::Mode (){
	return mode;
}
const simb::Origin_t pfp_container::Origin (){
	return origin;
}
const int pfp_container::IsCC (){
	return is_cc;
}

const double pfp_container::mcVtxX  (){
	return mc_vtx_x;
}
const double pfp_container::mcVtxY  (){
	return mc_vtx_y;
}
const double pfp_container::mcVtxZ  (){
	return mc_vtx_z;
}
const double pfp_container::pfpVtxX (){
	return pfp_vtx_x;
}
const double pfp_container::pfpVtxY (){
	return pfp_vtx_y;
}
const double pfp_container::pfpVtxZ (){
	return pfp_vtx_z;
}

const double pfp_container::Completeness (){
	return completeness;
}
const double pfp_container::Purity (){
	return purity;
}

const int pfp_container::NumMCHits   () {
	return n_mc_hits;
}
const int pfp_container::NumMCHitsU  () {
	return n_mc_hits_u;
}
const int pfp_container::NumMCHitsV  () {
	return n_mc_hits_v;
}
const int pfp_container::NumMCHitsW  () {
	return n_mc_hits_w;
}
const int pfp_container::NumPFPHits  () {
	return n_pfp_hits;
}
const int pfp_container::NumPFPHitsU () {
	return n_pfp_hits_u;
}
const int pfp_container::NumPFPHitsV () {
	return n_pfp_hits_v;
}
const int pfp_container::NumPFPHitsW () {
	return n_pfp_hits_w;
}

const int pfp_container::NumMatchedHits  () {
	return n_matched_hits;
}
const int pfp_container::NumMatchedHitsU () {
	return n_matched_hits_u;
}
const int pfp_container::NumMatchedHitsV () {
	return n_matched_hits_v;
}
const int pfp_container::NumMatchedHitsW () {
	return n_matched_hits_w;
}



}//end namespace
