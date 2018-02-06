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
void FillPostCutVector(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                       std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                       double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                       double mc_nu_vtx_x, double mc_nu_vtx_y, double mc_nu_vtx_z,
                       std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v);
//***************************************************************************
//***************************************************************************
void PrintPostCutVector(std::vector<std::tuple<int, int, double, double, double,
                                               std::string, std::string, int, int, double> > * post_cuts_v, bool _post_cuts_verbose);
//***************************************************************************
//***************************************************************************
void PostCutVectorPlots(std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v,
                        bool _post_cuts_verbose, TH1 * post_cuts_num_showers_purity_qe,
                        TH1 * post_cuts_num_showers_purity_res,
                        TH1 * post_cuts_num_showers_purity_dis,
                        TH1 * post_cuts_num_showers_purity_coh,
                        TH1 * post_cuts_num_showers_purity_mec);
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
void PostCutOpenAngle1Shower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                             double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                             TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_mixed,
                             TH1D * h_leading_shower_open_angle_numu_cc, TH1D * h_leading_shower_open_angle_nc,
                             TH1D * h_leading_shower_open_angle_cosmic, TH1D * h_leading_shower_open_angle_nc_pi0,
                             TH1D * h_leading_shower_open_angle_numu_cc_mixed, TH1D * h_leading_shower_open_angle_other_mixed,
                             TH1D * h_leading_shower_open_angle_unmatched);
//***************************************************************************
//***************************************************************************
void PostCutOpenAngle2PlusShower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
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
void NumShowersOpenAngle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                         double vtxX, double vtxY, double vtxZ,
                         TH1D * h_pfp_shower_open_angle_nue_cc_qe,
                         TH1D * h_pfp_shower_open_angle_nue_cc_out_fv,
                         TH1D * h_pfp_shower_open_angle_nue_cc_res,
                         TH1D * h_pfp_shower_open_angle_nue_cc_dis,
                         TH1D * h_pfp_shower_open_angle_nue_cc_coh,
                         TH1D * h_pfp_shower_open_angle_nue_cc_mec,
                         TH1D * h_pfp_shower_open_angle_nc,
                         TH1D * h_pfp_shower_open_angle_numu_cc_qe,
                         TH1D * h_pfp_shower_open_angle_numu_cc_res,
                         TH1D * h_pfp_shower_open_angle_numu_cc_dis,
                         TH1D * h_pfp_shower_open_angle_numu_cc_coh,
                         TH1D * h_pfp_shower_open_angle_numu_cc_mec,
                         TH1D * h_pfp_shower_open_angle_nc_pi0,
                         TH1D * h_pfp_shower_open_angle_nue_cc_mixed,
                         TH1D * h_pfp_shower_open_angle_numu_cc_mixed,
                         TH1D * h_pfp_shower_open_angle_cosmic,
                         TH1D * h_pfp_shower_open_angle_other_mixed,
                         TH1D * h_pfp_shower_open_angle_unmatched
                         );
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
                        std::vector<int> * no_track, std::vector<int> * has_track,
                        std::vector<int> * _1_shwr, std::vector<int> * _2_shwr,
                        std::vector<int> * _3_shwr, std::vector<int> * _4_shwr);
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
void FlashTot0(std::vector< double> largest_flash_v, double mc_nu_time, int mc_nu_id, std::vector<int> tabulated_origins,
               double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
               double vtxX, double vtxY, double vtxZ, TH1D * h_flash_t0_diff);
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
void HitLengthRatio(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                    double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                    double vtxX, double vtxY, double vtxZ,
                    TH1D * h_hit_length_ratio_nue_cc,
                    TH1D * h_hit_length_ratio_nue_cc_out_fv,
                    TH1D * h_hit_length_ratio_nue_cc_mixed,
                    TH1D * h_hit_length_ratio_numu_cc,
                    TH1D * h_hit_length_ratio_numu_cc_mixed,
                    TH1D * h_hit_length_ratio_nc,
                    TH1D * h_hit_length_ratio_nc_pi0,
                    TH1D * h_hit_length_ratio_cosmic,
                    TH1D * h_hit_length_ratio_other_mixed,
                    TH1D * h_hit_length_ratio_unmatched);
//***************************************************************************
//***************************************************************************
void TrackLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                 double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                 double vtxX, double vtxY, double vtxZ,
                 TH1D * h_trk_length_nue_cc,
                 TH1D * h_trk_length_nue_cc_out_fv,
                 TH1D * h_trk_length_nue_cc_mixed,
                 TH1D * h_trk_length_numu_cc,
                 TH1D * h_trk_length_numu_cc_mixed,
                 TH1D * h_trk_length_nc,
                 TH1D * h_trk_length_nc_pi0,
                 TH1D * h_trk_length_cosmic,
                 TH1D * h_trk_length_other_mixed,
                 TH1D * h_trk_length_unmatched);
//***************************************************************************
//***************************************************************************
void LongestTrackLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                        double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                        double vtxX, double vtxY, double vtxZ,
                        TH1D * h_trk_length_nue_cc,
                        TH1D * h_trk_length_nue_cc_out_fv,
                        TH1D * h_trk_length_nue_cc_mixed,
                        TH1D * h_trk_length_numu_cc,
                        TH1D * h_trk_length_numu_cc_mixed,
                        TH1D * h_trk_length_nc,
                        TH1D * h_trk_length_nc_pi0,
                        TH1D * h_trk_length_cosmic,
                        TH1D * h_trk_length_other_mixed,
                        TH1D * h_trk_length_unmatched);
//***************************************************************************
//***************************************************************************
void ShowerLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                  double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                  double vtxX, double vtxY, double vtxZ,
                  TH1D * h_shwr_length_nue_cc,
                  TH1D * h_shwr_length_nue_cc_out_fv,
                  TH1D * h_shwr_length_nue_cc_mixed,
                  TH1D * h_shwr_length_numu_cc,
                  TH1D * h_shwr_length_numu_cc_mixed,
                  TH1D * h_shwr_length_nc,
                  TH1D * h_shwr_length_nc_pi0,
                  TH1D * h_shwr_length_cosmic,
                  TH1D * h_shwr_length_other_mixed,
                  TH1D * h_shwr_length_unmatched);
//***************************************************************************
//***************************************************************************
void LongestShowerLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                         double vtxX, double vtxY, double vtxZ,
                         TH1D * h_shwr_length_nue_cc,
                         TH1D * h_shwr_length_nue_cc_out_fv,
                         TH1D * h_shwr_length_nue_cc_mixed,
                         TH1D * h_shwr_length_numu_cc,
                         TH1D * h_shwr_length_numu_cc_mixed,
                         TH1D * h_shwr_length_nc,
                         TH1D * h_shwr_length_nc_pi0,
                         TH1D * h_shwr_length_cosmic,
                         TH1D * h_shwr_length_other_mixed,
                         TH1D * h_shwr_length_unmatched);
//***************************************************************************
//***************************************************************************
void LeadingShowerLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                         double vtxX, double vtxY, double vtxZ,
                         TH1D * h_shwr_length_nue_cc,
                         TH1D * h_shwr_length_nue_cc_out_fv,
                         TH1D * h_shwr_length_nue_cc_mixed,
                         TH1D * h_shwr_length_numu_cc,
                         TH1D * h_shwr_length_numu_cc_mixed,
                         TH1D * h_shwr_length_nc,
                         TH1D * h_shwr_length_nc_pi0,
                         TH1D * h_shwr_length_cosmic,
                         TH1D * h_shwr_length_other_mixed,
                         TH1D * h_shwr_length_unmatched);
//***************************************************************************
//***************************************************************************
void LeadingShowerTrackLengths(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                               double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                               double vtxX, double vtxY, double vtxZ,
                               TH1D * h_shwr_trk_length_nue_cc,
                               TH1D * h_shwr_trk_length_nue_cc_out_fv,
                               TH1D * h_shwr_trk_length_nue_cc_mixed,
                               TH1D * h_shwr_trk_length_numu_cc,
                               TH1D * h_shwr_trk_length_numu_cc_mixed,
                               TH1D * h_shwr_trk_length_nc,
                               TH1D * h_shwr_trk_length_nc_pi0,
                               TH1D * h_shwr_trk_length_cosmic,
                               TH1D * h_shwr_trk_length_other_mixed,
                               TH1D * h_shwr_trk_length_unmatched);
//***************************************************************************
//***************************************************************************
void LongestShowerTrackLengths(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                               double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                               double vtxX, double vtxY, double vtxZ,
                               TH1D * h_shwr_trk_length_nue_cc,
                               TH1D * h_shwr_trk_length_nue_cc_out_fv,
                               TH1D * h_shwr_trk_length_nue_cc_mixed,
                               TH1D * h_shwr_trk_length_numu_cc,
                               TH1D * h_shwr_trk_length_numu_cc_mixed,
                               TH1D * h_shwr_trk_length_nc,
                               TH1D * h_shwr_trk_length_nc_pi0,
                               TH1D * h_shwr_trk_length_cosmic,
                               TH1D * h_shwr_trk_length_other_mixed,
                               TH1D * h_shwr_trk_length_unmatched);
//***************************************************************************
//***************************************************************************
void PlaneHitsComparisonShower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                               double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                               double vtxX, double vtxY, double vtxZ,
                               TH2D * h_collection_total_hits_shower_nue_cc,
                               TH2D * h_collection_total_hits_shower_nue_cc_out_fv,
                               TH2D * h_collection_total_hits_shower_nue_cc_mixed,
                               TH2D * h_collection_total_hits_shower_numu_cc,
                               TH2D * h_collection_total_hits_shower_numu_cc_mixed,
                               TH2D * h_collection_total_hits_shower_nc,
                               TH2D * h_collection_total_hits_shower_nc_pi0,
                               TH2D * h_collection_total_hits_shower_cosmic,
                               TH2D * h_collection_total_hits_shower_other_mixed,
                               TH2D * h_collection_total_hits_shower_unmatched);
//***************************************************************************
//***************************************************************************
void PlaneHitsComparisonLeadingShower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                      double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                      double vtxX, double vtxY, double vtxZ,
                                      TH2D * h_collection_total_hits_shower_nue_cc,
                                      TH2D * h_collection_total_hits_shower_nue_cc_out_fv,
                                      TH2D * h_collection_total_hits_shower_nue_cc_mixed,
                                      TH2D * h_collection_total_hits_shower_numu_cc,
                                      TH2D * h_collection_total_hits_shower_numu_cc_mixed,
                                      TH2D * h_collection_total_hits_shower_nc,
                                      TH2D * h_collection_total_hits_shower_nc_pi0,
                                      TH2D * h_collection_total_hits_shower_cosmic,
                                      TH2D * h_collection_total_hits_shower_other_mixed,
                                      TH2D * h_collection_total_hits_shower_unmatched);
//***************************************************************************
//***************************************************************************
void PlaneHitsComparisonTrack(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                              double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                              double vtxX, double vtxY, double vtxZ,
                              TH2D * h_collection_total_hits_track_nue_cc,
                              TH2D * h_collection_total_hits_track_nue_cc_out_fv,
                              TH2D * h_collection_total_hits_track_nue_cc_mixed,
                              TH2D * h_collection_total_hits_track_numu_cc,
                              TH2D * h_collection_total_hits_track_numu_cc_mixed,
                              TH2D * h_collection_total_hits_track_nc,
                              TH2D * h_collection_total_hits_track_nc_pi0,
                              TH2D * h_collection_total_hits_track_cosmic,
                              TH2D * h_collection_total_hits_track_other_mixed,
                              TH2D * h_collection_total_hits_track_unmatched);
//***************************************************************************
//***************************************************************************
void HitsPlots1D(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                 double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                 double vtxX, double vtxY, double vtxZ,
                 TH1D * h_collection_hits_track_nue_cc,
                 TH1D * h_collection_hits_track_nue_cc_out_fv,
                 TH1D * h_collection_hits_track_nue_cc_mixed,
                 TH1D * h_collection_hits_track_numu_cc,
                 TH1D * h_collection_hits_track_numu_cc_mixed,
                 TH1D * h_collection_hits_track_nc,
                 TH1D * h_collection_hits_track_nc_pi0,
                 TH1D * h_collection_hits_track_cosmic,
                 TH1D * h_collection_hits_track_other_mixed,
                 TH1D * h_collection_hits_track_unmatched,
                 TH1D * h_collection_hits_shower_nue_cc,
                 TH1D * h_collection_hits_shower_nue_cc_out_fv,
                 TH1D * h_collection_hits_shower_nue_cc_mixed,
                 TH1D * h_collection_hits_shower_numu_cc,
                 TH1D * h_collection_hits_shower_numu_cc_mixed,
                 TH1D * h_collection_hits_shower_nc,
                 TH1D * h_collection_hits_shower_nc_pi0,
                 TH1D * h_collection_hits_shower_cosmic,
                 TH1D * h_collection_hits_shower_other_mixed,
                 TH1D * h_collection_hits_shower_unmatched,
                 TH1D * h_collection_hits_leading_shower_nue_cc,
                 TH1D * h_collection_hits_leading_shower_nue_cc_out_fv,
                 TH1D * h_collection_hits_leading_shower_nue_cc_mixed,
                 TH1D * h_collection_hits_leading_shower_numu_cc,
                 TH1D * h_collection_hits_leading_shower_numu_cc_mixed,
                 TH1D * h_collection_hits_leading_shower_nc,
                 TH1D * h_collection_hits_leading_shower_nc_pi0,
                 TH1D * h_collection_hits_leading_shower_cosmic,
                 TH1D * h_collection_hits_leading_shower_other_mixed,
                 TH1D * h_collection_hits_leading_shower_unmatched,
                 TH1D * h_total_hits_leading_shower_nue_cc,
                 TH1D * h_total_hits_leading_shower_nue_cc_out_fv,
                 TH1D * h_total_hits_leading_shower_nue_cc_mixed,
                 TH1D * h_total_hits_leading_shower_numu_cc,
                 TH1D * h_total_hits_leading_shower_numu_cc_mixed,
                 TH1D * h_total_hits_leading_shower_nc,
                 TH1D * h_total_hits_leading_shower_nc_pi0,
                 TH1D * h_total_hits_leading_shower_cosmic,
                 TH1D * h_total_hits_leading_shower_other_mixed,
                 TH1D * h_total_hits_leading_shower_unmatched);
//***************************************************************************
//***************************************************************************
void EnergyHits(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0, bool _verbose,
                double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                double vtxX, double vtxY, double vtxZ, double mc_nu_energy, double mc_ele_energy,
                TH2D * h_ele_eng_total_hits, TH2D * h_ele_eng_colleciton_hits, TH2D * h_nu_eng_total_hits, TH2D * h_nu_eng_collection_hits);
//***************************************************************************
//***************************************************************************
int MapFailureCutToString(const std::string failure_cut);
//***************************************************************************
//***************************************************************************
void FailureReason(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                   std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0, bool _verbose,
                   double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                   double vtxX, double vtxY, double vtxZ,
                   TH1D * h_failure_reason_nue_cc,
                   TH1D * h_failure_reason_nue_cc_out_fv,
                   TH1D * h_failure_reason_nue_cc_mixed,
                   TH1D * h_failure_reason_numu_cc,
                   TH1D * h_failure_reason_numu_cc_mixed,
                   TH1D * h_failure_reason_nc,
                   TH1D * h_failure_reason_nc_pi0,
                   TH1D * h_failure_reason_cosmic,
                   TH1D * h_failure_reason_other_mixed,
                   TH1D * h_failure_reason_unmatched);
//***************************************************************************
//***************************************************************************
void LeadingCosTheta(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                     double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                     double vtxX, double vtxY, double vtxZ,
                     TH1D * h_ele_cos_theta_nue_cc,
                     TH1D * h_ele_cos_theta_nue_cc_out_fv,
                     TH1D * h_ele_cos_theta_nue_cc_mixed,
                     TH1D * h_ele_cos_theta_numu_cc,
                     TH1D * h_ele_cos_theta_numu_cc_mixed,
                     TH1D * h_ele_cos_theta_nc,
                     TH1D * h_ele_cos_theta_nc_pi0,
                     TH1D * h_ele_cos_theta_cosmic,
                     TH1D * h_ele_cos_theta_other_mixed,
                     TH1D * h_ele_cos_theta_unmatched);
//***************************************************************************
//***************************************************************************
//***************************************************************************
//***************************************************************************
};
#endif
