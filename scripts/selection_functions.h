#ifndef SELECTION_FUNCTIONS_h
#define SELECTION_FUNCTIONS_h

#include "../xsecAna/TpcObjectContainer.h"
#include "../xsecAna/ParticleContainer.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMarker.h"

#include <iostream>
#include <vector>

#include "../xsecAna/LinkDef.h"




class selection_functions {

public:

selection_functions()=default;

//***************************************************************************
//***************************************************************************
int flash_in_time(int flash_pe, int flash_pe_threshold, double flash_time, double flash_start, double flash_end);
//***************************************************************************
//***************************************************************************
void loop_flashes(TFile * f, TTree * optical_tree, int flash_pe_threshold, double flash_time_start,
                  double flash_time_end, std::vector<int> * _passed_runs);
//***************************************************************************
//***************************************************************************
bool in_fv(double x, double y, double z,
           double x1, double x2, double y1,
           double y2, double z1, double z2);
//***************************************************************************
//***************************************************************************
void fiducial_volume_cut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                         std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
bool opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance);
//***************************************************************************
//***************************************************************************
bool opt_vtx_distance_width(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double flash_width_z, double tolerance);
//***************************************************************************
//***************************************************************************
void SetXYflashVector(TFile * f, TTree * optical_tree, std::vector< std::vector< double> > * largest_flash_v_v);
//***************************************************************************
//***************************************************************************
void flashRecoVtxDist(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                      double tolerance, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
bool shwr_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z, double tolerance);
//***************************************************************************
//***************************************************************************
//this function wants to remove particles too far from the reconstructed neutrino vertex
void VtxNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                   double tolerance, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
//this function wants to remove particles too far from the reconstructed neutrino vertex
void VtxTrackNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        double tolerance, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
//this function wants to remove particles too far from the reconstructed neutrino vertex
void HitThreshold(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                  double threshold, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
//this gives a list of all of the origins of the tpc objects
void GetOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::string> * tpco_origin_v);
//***************************************************************************
//***************************************************************************
//this function simply checks if the tpc object is a nue
void HasNue(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
void OpenAngleCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco, const double tolerance_open_angle,
                  const bool _verbose);
//***************************************************************************
//***************************************************************************
void dEdxCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
             const double tolerance_dedx_min, const double tolerance_dedx_max, const bool _verbose);
//***************************************************************************
//***************************************************************************
//this function just counts if at least 1 tpc object passes the cuts
bool ValidTPCObjects(std::vector<std::pair<int, std::string> > * passed_tpco);
//***************************************************************************
//***************************************************************************
std::vector<int> TabulateOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
                                 double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ);
//***************************************************************************
//***************************************************************************
void TotalOrigins(std::vector<int> tabulated_origins, std::vector<int> * total_cut_origins);
//***************************************************************************
//***************************************************************************
//modify this so it takes a string of the cut name so I only pass it a few variable at a time,
//then I can call this function several times later at the bottom
void PrintInfo(int mc_nue_cc_counter,
               std::vector<int> * counter_v,
               std::string cut_name);
//***************************************************************************
//***************************************************************************
double calcNumNucleons(double _x1, double _x2, double _y1,
                       double _y2, double _z1, double _z2);
//***************************************************************************
//***************************************************************************
void calcXSec(double _x1, double _x2, double _y1,
              double _y2, double _z1, double _z2,
              int n_total, int n_bkg, double flux, double efficiency, std::vector<double>  * xsec_cc);
//***************************************************************************
//***************************************************************************
void xsec_plot(bool _verbose, double genie_xsec, double xsec, double average_energy, double stat_error);
//***************************************************************************
//***************************************************************************
void PostCutPlots(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                  TH2I * h_tracks_showers, TH2I * h_tracks_showers_cosmic, TH2I * h_tracks_showers_numu,
                  TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_mixed,
                  TH1D * h_leading_shower_open_angle_numu_cc, TH1D * h_leading_shower_open_angle_numu_nc,
                  TH1D * h_leading_shower_open_angle_cosmic, TH1D * h_leading_shower_open_angle_nue_nc,
                  TH1D * h_leading_shower_open_angle_numu_cc_mixed, TH1D * h_leading_shower_open_angle_other_mixed,
                  TH1D * h_leading_shower_open_angle_unmatched,
                  TH1D * h_trk_vtx_dist_nue_cc, TH1D * h_trk_vtx_dist_nue_cc_mixed,
                  TH1D * h_trk_vtx_dist_numu_cc, TH1D * h_trk_vtx_dist_numu_nc,
                  TH1D * h_trk_vtx_dist_cosmic, TH1D * h_trk_vtx_dist_nue_nc,
                  TH1D * h_trk_vtx_dist_numu_cc_mixed, TH1D * h_trk_vtx_dist_other_mixed,
                  TH1D * h_trk_vtx_dist_unmatched);
//***************************************************************************
//***************************************************************************
void TopologyPlots(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
                   double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                   TH2D * h_pfp_track_shower_nue_cc_qe,
                   TH2D * h_pfp_track_shower_nue_cc_out_fv,
                   TH2D * h_pfp_track_shower_nue_cc_res,
                   TH2D * h_pfp_track_shower_nue_cc_dis,
                   TH2D * h_pfp_track_shower_nue_cc_coh,
                   TH2D * h_pfp_track_shower_nue_cc_mec,
                   TH2D * h_pfp_track_shower_nue_nc,
                   TH2D * h_pfp_track_shower_numu_cc_qe,
                   TH2D * h_pfp_track_shower_numu_cc_res,
                   TH2D * h_pfp_track_shower_numu_cc_dis,
                   TH2D * h_pfp_track_shower_numu_cc_coh,
                   TH2D * h_pfp_track_shower_numu_cc_mec,
                   TH2D * h_pfp_track_shower_numu_nc,
                   TH2D * h_pfp_track_shower_nue_cc_mixed,
                   TH2D * h_pfp_track_shower_numu_cc_mixed,
                   TH2D * h_pfp_track_shower_cosmic,
                   TH2D * h_pfp_track_shower_other_mixed,
                   TH2D * h_pfp_track_shower_unmatched,
                   TH2D * h_leading_shower_mc_pdg_nue_cc_qe,
                   TH2D * h_leading_shower_mc_pdg_nue_cc_out_fv,
                   TH2D * h_leading_shower_mc_pdg_nue_cc_res,
                   TH2D * h_leading_shower_mc_pdg_nue_cc_dis,
                   TH2D * h_leading_shower_mc_pdg_nue_cc_coh,
                   TH2D * h_leading_shower_mc_pdg_nue_cc_mec,
                   TH2D * h_leading_shower_mc_pdg_nue_nc,
                   TH2D * h_leading_shower_mc_pdg_numu_cc_qe,
                   TH2D * h_leading_shower_mc_pdg_numu_cc_res,
                   TH2D * h_leading_shower_mc_pdg_numu_cc_dis,
                   TH2D * h_leading_shower_mc_pdg_numu_cc_coh,
                   TH2D * h_leading_shower_mc_pdg_numu_cc_mec,
                   TH2D * h_leading_shower_mc_pdg_numu_nc,
                   TH2D * h_leading_shower_mc_pdg_nue_cc_mixed,
                   TH2D * h_leading_shower_mc_pdg_numu_cc_mixed,
                   TH2D * h_leading_shower_mc_pdg_cosmic,
                   TH2D * h_leading_shower_mc_pdg_other_mixed,
                   TH2D * h_leading_shower_mc_pdg_unmatched,
                   TH1D * h_pfp_track_nue_cc_qe,
                   TH1D * h_pfp_track_nue_cc_out_fv,
                   TH1D * h_pfp_track_nue_cc_res,
                   TH1D * h_pfp_track_nue_cc_dis,
                   TH1D * h_pfp_track_nue_cc_coh,
                   TH1D * h_pfp_track_nue_cc_mec,
                   TH1D * h_pfp_track_nue_nc,
                   TH1D * h_pfp_track_numu_cc_qe,
                   TH1D * h_pfp_track_numu_cc_res,
                   TH1D * h_pfp_track_numu_cc_dis,
                   TH1D * h_pfp_track_numu_cc_coh,
                   TH1D * h_pfp_track_numu_cc_mec,
                   TH1D * h_pfp_track_numu_nc,
                   TH1D * h_pfp_track_nue_cc_mixed,
                   TH1D * h_pfp_track_numu_cc_mixed,
                   TH1D * h_pfp_track_cosmic,
                   TH1D * h_pfp_track_other_mixed,
                   TH1D * h_pfp_track_unmatched,
                   TH1D * h_pfp_shower_nue_cc_qe,
                   TH1D * h_pfp_shower_nue_cc_out_fv,
                   TH1D * h_pfp_shower_nue_cc_res,
                   TH1D * h_pfp_shower_nue_cc_dis,
                   TH1D * h_pfp_shower_nue_cc_coh,
                   TH1D * h_pfp_shower_nue_cc_mec,
                   TH1D * h_pfp_shower_nue_nc,
                   TH1D * h_pfp_shower_numu_cc_qe,
                   TH1D * h_pfp_shower_numu_cc_res,
                   TH1D * h_pfp_shower_numu_cc_dis,
                   TH1D * h_pfp_shower_numu_cc_coh,
                   TH1D * h_pfp_shower_numu_cc_mec,
                   TH1D * h_pfp_shower_numu_nc,
                   TH1D * h_pfp_shower_nue_cc_mixed,
                   TH1D * h_pfp_shower_numu_cc_mixed,
                   TH1D * h_pfp_shower_cosmic,
                   TH1D * h_pfp_shower_other_mixed,
                   TH1D * h_pfp_shower_unmatched    );
//***************************************************************************
//***************************************************************************

};

#endif
