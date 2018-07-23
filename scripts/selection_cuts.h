#ifndef SELECTION_CUTS_h
#define SELECTION_CUTS_h

#include "../xsecAna/TpcObjectContainer.h"
#include "../xsecAna/ParticleContainer.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <vector>

#include "../xsecAna/LinkDef.h"

class selection_cuts {

public:

selection_cuts()=default;


//***************************************************************************
bool flash_in_time(double flash_time, double flash_start, double flash_end);
//***************************************************************************
bool flash_pe(int flash_pe, int flash_pe_threshold);
//***************************************************************************
void loop_flashes(TFile * f, TTree * optical_tree, int flash_pe_threshold, double flash_time_start,
                  double flash_time_end, std::vector<int> * _passed_runs, std::vector<std::pair<double, int> > * flash_time, const int stream);
//***************************************************************************
bool in_fv(double x, double y, double z, std::vector<double> fv_boundary_v);
//***************************************************************************
void fiducial_volume_cut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                         std::vector<double> fv_boundary_v,
                         std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
bool opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance);
//***************************************************************************
bool opt_vtx_distance_width(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double flash_width_z, double tolerance);
//***************************************************************************
void SetXYflashVector(TFile * f, TTree * optical_tree, std::vector< std::vector< double> > * largest_flash_v_v,
                      double flash_time_start, double flash_time_end, double flash_pe_threshold);
//***************************************************************************
void flashRecoVtxDist(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                      double tolerance, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
bool shwr_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z, double tolerance);
//***************************************************************************
void VtxNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                   double tolerance, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
void VtxTrackNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        double tolerance, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
void HitThreshold(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                  double threshold, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
//this gives a list of all of the origins of the tpc objects
void GetOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::string> * tpco_origin_v);
//***************************************************************************
void HasNue(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
void OpenAngleCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
                  const std::vector<double> tolerance_open_angle, const bool _verbose);
//***************************************************************************
void dEdxCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
             const double tolerance_dedx_min, const double tolerance_dedx_max, const bool _verbose);
//***************************************************************************
void SecondaryShowersDistCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, const double dist_tolerance);
//***************************************************************************
void HitLengthRatioCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                       const double pfp_hits_length_tolerance);
//***************************************************************************
void HitThresholdCollection(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                            double threshold, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose);
//***************************************************************************
void LongestTrackLeadingShowerCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                  const double ratio_tolerance);
//***************************************************************************
bool IsContained(std::vector<double> track_start, std::vector<double> track_end, std::vector<double> fv_boundary_v);
//***************************************************************************
void ContainedTracksCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                        std::vector<double> fv_boundary_v, const bool enabled);
//***************************************************************************




};

#endif
