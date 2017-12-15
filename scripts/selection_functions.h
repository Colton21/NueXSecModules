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
bool flash_in_time(double flash_time, double flash_start, double flash_end);
//***************************************************************************
//***************************************************************************
bool flash_pe(int flash_pe, int flash_pe_threshold);
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
void SetXYflashVector(TFile * f, TTree * optical_tree, std::vector< std::vector< double> > * largest_flash_v_v,
                      double flash_time_start, double flash_time_end);
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
void OpenAngleCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
                  const std::vector<double> tolerance_open_angle, const bool _verbose);
//***************************************************************************
//***************************************************************************
void PostCutsdEdx(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                  double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                  double mc_nu_vtx_x, double mc_nu_vtx_y, double mc_nu_vtx_z,
                  TH1D * h_dedx_cuts_nue_cc,
                  TH1D * h_dedx_cuts_nue_cc_mixed,
                  TH1D * h_dedx_cuts_nue_cc_out_fv,
                  TH1D * h_dedx_cuts_numu_cc,
                  TH1D * h_dedx_cuts_nc,
                  TH1D * h_dedx_cuts_cosmic,
                  TH1D * h_dedx_cuts_nc_pi0,
                  TH1D * h_dedx_cuts_numu_cc_mixed,
                  TH1D * h_dedx_cuts_other_mixed,
                  TH1D * h_dedx_cuts_unmatched     );
//***************************************************************************
//***************************************************************************
void dEdxCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
             const double tolerance_dedx_min, const double tolerance_dedx_max, const bool _verbose);
//***************************************************************************
//***************************************************************************
void FillPostCutVector(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                       std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                       double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                       double mc_nu_vtx_x, double mc_nu_vtx_y, double mc_nu_vtx_z,
                       std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int> > * post_cuts_v);
//***************************************************************************
//***************************************************************************
void PrintPostCutVector(std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int> > * post_cuts_v, bool _post_cuts_verbose);
//***************************************************************************
//***************************************************************************
//this function just counts if at least 1 tpc object passes the cuts
bool ValidTPCObjects(std::vector<std::pair<int, std::string> > * passed_tpco);
//***************************************************************************
//***************************************************************************
std::vector<int> TabulateOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
                                 bool has_pi0, double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ);
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
std::pair<std::string, int> TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool has_pi0,
                                            double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ);
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
void PostCutOpenAngle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                      double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                      TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_mixed,
                      TH1D * h_leading_shower_open_angle_numu_cc, TH1D * h_leading_shower_open_angle_nc,
                      TH1D * h_leading_shower_open_angle_cosmic, TH1D * h_leading_shower_open_angle_nc_pi0,
                      TH1D * h_leading_shower_open_angle_numu_cc_mixed, TH1D * h_leading_shower_open_angle_other_mixed,
                      TH1D * h_leading_shower_open_angle_unmatched);
//***************************************************************************
//***************************************************************************
void PostCutTrkVtx(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                   bool has_pi0,
                   double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                   TH1D * h_trk_vtx_dist_nue_cc, TH1D * h_trk_vtx_dist_nue_cc_mixed,
                   TH1D * h_trk_vtx_dist_numu_cc, TH1D * h_trk_vtx_dist_nc,
                   TH1D * h_trk_vtx_dist_cosmic, TH1D * h_trk_vtx_dist_nc_pi0,
                   TH1D * h_trk_vtx_dist_numu_cc_mixed, TH1D * h_trk_vtx_dist_other_mixed,
                   TH1D * h_trk_vtx_dist_unmatched);
//***************************************************************************
//***************************************************************************
void TopologyPlots1(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                    std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                    double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                    TH2D * h_pfp_track_shower_nue_cc_qe,
                    TH2D * h_pfp_track_shower_nue_cc_out_fv,
                    TH2D * h_pfp_track_shower_nue_cc_res,
                    TH2D * h_pfp_track_shower_nue_cc_dis,
                    TH2D * h_pfp_track_shower_nue_cc_coh,
                    TH2D * h_pfp_track_shower_nue_cc_mec,
                    TH2D * h_pfp_track_shower_nc,
                    TH2D * h_pfp_track_shower_numu_cc_qe,
                    TH2D * h_pfp_track_shower_numu_cc_res,
                    TH2D * h_pfp_track_shower_numu_cc_dis,
                    TH2D * h_pfp_track_shower_numu_cc_coh,
                    TH2D * h_pfp_track_shower_numu_cc_mec,
                    TH2D * h_pfp_track_shower_nc_pi0,
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
                    TH2D * h_leading_shower_mc_pdg_nc,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_qe,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_res,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_dis,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_coh,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_mec,
                    TH2D * h_leading_shower_mc_pdg_nc_pi0,
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
                    TH1D * h_pfp_track_nc,
                    TH1D * h_pfp_track_numu_cc_qe,
                    TH1D * h_pfp_track_numu_cc_res,
                    TH1D * h_pfp_track_numu_cc_dis,
                    TH1D * h_pfp_track_numu_cc_coh,
                    TH1D * h_pfp_track_numu_cc_mec,
                    TH1D * h_pfp_track_nc_pi0,
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
                    TH1D * h_pfp_shower_nc,
                    TH1D * h_pfp_shower_numu_cc_qe,
                    TH1D * h_pfp_shower_numu_cc_res,
                    TH1D * h_pfp_shower_numu_cc_dis,
                    TH1D * h_pfp_shower_numu_cc_coh,
                    TH1D * h_pfp_shower_numu_cc_mec,
                    TH1D * h_pfp_shower_nc_pi0,
                    TH1D * h_pfp_shower_nue_cc_mixed,
                    TH1D * h_pfp_shower_numu_cc_mixed,
                    TH1D * h_pfp_shower_cosmic,
                    TH1D * h_pfp_shower_other_mixed,
                    TH1D * h_pfp_shower_unmatched    );
//***************************************************************************
//***************************************************************************
void TopologyPlots2(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                    std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                    double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                    TH2D * h_pfp_track_shower_nue_cc_qe,
                    TH2D * h_pfp_track_shower_nue_cc_out_fv,
                    TH2D * h_pfp_track_shower_nue_cc_res,
                    TH2D * h_pfp_track_shower_nue_cc_dis,
                    TH2D * h_pfp_track_shower_nue_cc_coh,
                    TH2D * h_pfp_track_shower_nue_cc_mec,
                    TH2D * h_pfp_track_shower_nc,
                    TH2D * h_pfp_track_shower_numu_cc_qe,
                    TH2D * h_pfp_track_shower_numu_cc_res,
                    TH2D * h_pfp_track_shower_numu_cc_dis,
                    TH2D * h_pfp_track_shower_numu_cc_coh,
                    TH2D * h_pfp_track_shower_numu_cc_mec,
                    TH2D * h_pfp_track_shower_nc_pi0,
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
                    TH2D * h_leading_shower_mc_pdg_nc,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_qe,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_res,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_dis,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_coh,
                    TH2D * h_leading_shower_mc_pdg_numu_cc_mec,
                    TH2D * h_leading_shower_mc_pdg_nc_pi0,
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
                    TH1D * h_pfp_track_nc,
                    TH1D * h_pfp_track_numu_cc_qe,
                    TH1D * h_pfp_track_numu_cc_res,
                    TH1D * h_pfp_track_numu_cc_dis,
                    TH1D * h_pfp_track_numu_cc_coh,
                    TH1D * h_pfp_track_numu_cc_mec,
                    TH1D * h_pfp_track_nc_pi0,
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
                    TH1D * h_pfp_shower_nc,
                    TH1D * h_pfp_shower_numu_cc_qe,
                    TH1D * h_pfp_shower_numu_cc_res,
                    TH1D * h_pfp_shower_numu_cc_dis,
                    TH1D * h_pfp_shower_numu_cc_coh,
                    TH1D * h_pfp_shower_numu_cc_mec,
                    TH1D * h_pfp_shower_nc_pi0,
                    TH1D * h_pfp_shower_nue_cc_mixed,
                    TH1D * h_pfp_shower_numu_cc_mixed,
                    TH1D * h_pfp_shower_cosmic,
                    TH1D * h_pfp_shower_other_mixed,
                    TH1D * h_pfp_shower_unmatched    );
//***************************************************************************
//***************************************************************************
void PostCutsVtxFlash(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                      double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                      TH1D * h_vtx_flash_nue_cc, TH1D * h_vtx_flash_nue_cc_mixed,
                      TH1D * h_vtx_flash_numu_cc, TH1D * h_vtx_flash_nc,
                      TH1D * h_vtx_flash_cosmic, TH1D * h_vtx_flash_nc_pi0,
                      TH1D * h_vtx_flash_numu_cc_mixed, TH1D * h_vtx_flash_other_mixed,
                      TH1D * h_vtx_flash_unmatched);
//***************************************************************************
//***************************************************************************
void PostCutsShwrVtx(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                     double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                     TH1D * h_shwr_vtx_dist_nue_cc,
                     TH1D * h_shwr_vtx_dist_nue_cc_mixed,
                     TH1D * h_shwr_vtx_dist_numu_cc,
                     TH1D * h_shwr_vtx_dist_nc,
                     TH1D * h_shwr_vtx_dist_cosmic,
                     TH1D * h_shwr_vtx_dist_nc_pi0,
                     TH1D * h_shwr_vtx_dist_numu_cc_mixed,
                     TH1D * h_shwr_vtx_dist_other_mixed,
                     TH1D * h_shwr_vtx_dist_unmatched     );
//***************************************************************************
//***************************************************************************
void PostCutHitThreshold(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                         double vtxX, double vtxY, double vtxZ,
                         double mc_nu_energy, double mc_ele_energy,
                         TH2D * h_shwr_hits_nu_eng, TH2D * h_shwr_hits_ele_eng);
//***************************************************************************
//***************************************************************************
void TopologyEfficiency(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                        double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                        double vtxX, double vtxY, double vtxZ,
                        std::vector<int> * no_track, std::vector<int> * has_track);
//***************************************************************************
//***************************************************************************
void SequentialTrueEnergyPlots(int mc_nu_id, double mc_nu_vtx_x, double mc_nu_vtx_y, double mc_nu_vtx_z,
                               double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                               std::vector<int> tabulated_origins, double mc_nu_energy,
                               double mc_ele_energy, TH1D * h_selected_nu_energy, TH1D * h_selected_ele_energy);
//***************************************************************************
//***************************************************************************
void ChargeShare(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                 double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                 double vtxX, double vtxY, double vtxZ, TH1D * h_charge_share_nue_cc_mixed);
//***************************************************************************
//***************************************************************************
void FlashTot0(std::vector< double> largest_flash_v, double mc_nu_time, int mc_nu_id, std::vector<int> tabulated_origins, TH1D * h_flash_t0_diff);
//***************************************************************************
//***************************************************************************
void dEdxVsOpenAngle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                     double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                     double vtxX, double vtxY, double vtxZ,
                     TH2D * h_dedx_open_angle_nue_cc,
                     TH2D * h_dedx_open_angle_nue_cc_out_fv,
                     TH2D * h_dedx_open_angle_nue_cc_mixed,
                     TH2D * h_dedx_open_angle_numu_cc,
                     TH2D * h_dedx_open_angle_numu_cc_mixed,
                     TH2D * h_dedx_open_angle_nc,
                     TH2D * h_dedx_open_angle_nc_pi0,
                     TH2D * h_dedx_open_angle_cosmic,
                     TH2D * h_dedx_open_angle_other_mixed,
                     TH2D * h_dedx_open_angle_unmatched);
//***************************************************************************
//***************************************************************************
void ShowerLengthvsHits(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                        double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                        double vtxX, double vtxY, double vtxZ,
                        TH2D * h_shwr_len_hits_nue_cc,
                        TH2D * h_shwr_len_hits_nue_cc_out_fv,
                        TH2D * h_shwr_len_hits_nue_cc_mixed,
                        TH2D * h_shwr_len_hits_numu_cc,
                        TH2D * h_shwr_len_hits_numu_cc_mixed,
                        TH2D * h_shwr_len_hits_nc,
                        TH2D * h_shwr_len_hits_nc_pi0,
                        TH2D * h_shwr_len_hits_cosmic,
                        TH2D * h_shwr_len_hits_other_mixed,
                        TH2D * h_shwr_len_hits_unmatched);
//***************************************************************************
//***************************************************************************
void SecondaryShowersDist(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                          double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                          double vtxX, double vtxY, double vtxZ,
                          TH1D * h_second_shwr_dist_nue_cc,
                          TH1D * h_second_shwr_dist_nue_cc_out_fv,
                          TH1D * h_second_shwr_dist_nue_cc_mixed,
                          TH1D * h_second_shwr_dist_numu_cc,
                          TH1D * h_second_shwr_dist_numu_cc_mixed,
                          TH1D * h_second_shwr_dist_nc,
                          TH1D * h_second_shwr_dist_nc_pi0,
                          TH1D * h_second_shwr_dist_cosmic,
                          TH1D * h_second_shwr_dist_other_mixed,
                          TH1D * h_second_shwr_dist_unmatched);
//***************************************************************************
//***************************************************************************
};


#endif
