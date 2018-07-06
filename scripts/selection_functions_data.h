#ifndef SELECTION_FUNCTIONS_DATA_h
#define SELECTION_FUNCTIONS_DATA_h

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




class selection_functions_data {

public:

selection_functions_data()=default;

void FillPostCutVectorData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                           std::vector<std::pair<int, std::string> > * passed_tpco,
                           std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v);
//***************************************************************************
//***************************************************************************
void TabulateOriginsData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco,
                         std::vector<int> * tabulated_origins_data);
//***************************************************************************
//***************************************************************************
static void PrintInfoData(int counter, std::string cut_name);
//***************************************************************************
//***************************************************************************
std::pair<std::string, int> TPCO_Classifier_Data(xsecAna::TPCObjectContainer tpc_obj);
//***************************************************************************
//***************************************************************************
void PostCutsdEdxData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                      TH1D * h_dedx_cuts_data);
//***************************************************************************
//***************************************************************************
void PostCutOpenAngleData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                          TH1D * h_leading_shower_open_angle_data);
//***************************************************************************
//***************************************************************************
void PostCutTrkVtxData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                       TH1D * h_trk_vtx_dist_data);
//***************************************************************************
//***************************************************************************
void TopologyPlots1Data(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        std::vector<std::pair<int, std::string> > * passed_tpco,
                        TH2D * h_pfp_track_shower_data,
                        TH1D * h_pfp_track_data,
                        TH1D * h_pfp_shower_data
                        );
//***************************************************************************
//***************************************************************************
void TopologyPlots2Data(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        std::vector<std::pair<int, std::string> > * passed_tpco,
                        TH2D * h_pfp_track_shower_data,
                        TH1D * h_pfp_track_data,
                        TH1D * h_pfp_shower_data
                        );
//***************************************************************************
//***************************************************************************
void PostCutsVtxFlashData(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                          std::vector<std::pair<int, std::string> > * passed_tpco, TH1D * h_vtx_flash_data);
//***************************************************************************
//***************************************************************************
void PostCutsShwrVtxData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                         TH1D * h_shwr_vtx_dist_data);
//***************************************************************************
//***************************************************************************
void TopologyEfficiencyData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                            std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                            std::vector<int> * no_track, std::vector<int> * has_track,
                            std::vector<int> * _1_shwr, std::vector<int> * _2_shwr,
                            std::vector<int> * _3_shwr, std::vector<int> * _4_shwr);
//***************************************************************************
//***************************************************************************
void dEdxVsOpenAngleData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                         TH2D * h_dedx_open_angle_data);
//***************************************************************************
//***************************************************************************
void ShowerLengthvsHitsData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                            std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                            TH2D * h_shwr_len_hits_data);
//***************************************************************************
//***************************************************************************
void SecondaryShowersDistData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                              TH1D * h_second_shwr_dist_data);
//***************************************************************************
//***************************************************************************
void HitLengthRatioData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                        TH1D * h_hit_length_ratio_data);
//***************************************************************************
//***************************************************************************
int MapFailureCutToStringData(const std::string failure_cut);
//***************************************************************************
//***************************************************************************
void FailureReasonData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                       std::vector<std::pair<int, std::string> > * passed_tpco,
                       TH1D * h_failure_reason_data);
//***************************************************************************
//***************************************************************************
void LeadingCosThetaData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco, const double theta_translation, const double phi_translation,
                         bool _verbose, TH1D * h_ele_cos_theta_data);
//***************************************************************************
//***************************************************************************
void HitsPlots1DData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                     TH1D * h_collection_hits_track_data,
                     TH1D * h_collection_hits_shower_data,
                     TH1D * h_collection_hits_leading_shower_data,
                     TH1D * h_total_hits_leading_shower_data);
//***************************************************************************
//***************************************************************************
void NumShowersOpenAngleData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                             std::vector<std::pair<int, std::string> > * passed_tpco, TH1D * h_pfp_shower_open_angle_data);
//***************************************************************************
//***************************************************************************
void PostCutOpenAngle1ShowerData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                 TH1D * h_leading_shower_open_angle_data);
//***************************************************************************
//***************************************************************************
void PostCutOpenAngle2PlusShowerData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                     TH1D * h_leading_shower_open_angle_data);
//***************************************************************************
//***************************************************************************
void PlaneHitsComparisonTrackData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                  TH2D * h_collection_total_hits_track_data);
//***************************************************************************
//***************************************************************************
void PlaneHitsComparisonShowerData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                   TH2D * h_collection_total_hits_shower_data);
//***************************************************************************
//***************************************************************************
void PlaneHitsComparisonLeadingShowerData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                          TH2D * h_collection_total_hits_shower_data);
//***************************************************************************
//***************************************************************************
void LeadingPhiData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                    TH1D * h_ele_pfp_phi_data);
//***************************************************************************
//***************************************************************************
void LeadingThetaData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                      TH1D * h_ele_pfp_theta_data);
//***************************************************************************
//***************************************************************************
void LeadingMomentumData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, TH1D * h_ele_pfp_momentum_data);
//***************************************************************************
//***************************************************************************
void LeadingMomentumTrackTopologyData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                      TH1D * h_ele_pfp_momentum_no_track_data,
                                      TH1D * h_ele_pfp_momentum_has_track_data);
//***************************************************************************
//***************************************************************************
void XYZPositionData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, std::vector<TH1 * > * h_ele_pfp_xyz_data);
//***************************************************************************
//***************************************************************************
void EnergyCosThetaData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                        TH2 * h_ele_eng_costheta_data);
//***************************************************************************
//***************************************************************************
void EnergyCosThetaSlicesData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                              const double theta_translation, const double phi_translation,
                              TH1 * h_ele_eng_for_data,
                              TH1 * h_ele_eng_mid_data,
                              TH1 * h_ele_eng_back_data);
//***************************************************************************
//***************************************************************************
void LeadingShowerLengthData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                             TH1D * h_shwr_length_data);
//***************************************************************************
//***************************************************************************
void LeadingShowerTrackLengthsData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                   TH1D * h_shwr_trk_length_data);
//***************************************************************************
//***************************************************************************
void PostCutsVtxFlashUpstreamData(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                  std::vector<std::pair<int, std::string> > * passed_tpco, TH1D * h_vtx_flash_data);
//***************************************************************************
//***************************************************************************
void PostCutsVtxFlashDownstreamData(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                    std::vector<std::pair<int, std::string> > * passed_tpco, TH1D * h_vtx_flash_data);
//***************************************************************************
//***************************************************************************
void PostCutsLeadingMomentumData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                 TH1D * h_ele_momentum_data);
//***************************************************************************
//***************************************************************************

};


#endif
