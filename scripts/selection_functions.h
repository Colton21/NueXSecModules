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
#include "TGraph.h"
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
                         std::vector<int> * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
bool opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance);
//***************************************************************************
//***************************************************************************
void SetXYflashVector(TFile * f, TTree * optical_tree, std::vector< std::vector< double> > * largest_flash_v_v);
//***************************************************************************
//***************************************************************************
void flashRecoVtxDist(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                      double tolerance, std::vector<int> * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
bool shwr_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z, double tolerance);
//***************************************************************************
//***************************************************************************
//this function wants to remove particles too far from the reconstructed neutrino vertex
void VtxNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                   double tolerance, std::vector<int> * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
//this function wants to remove particles too far from the reconstructed neutrino vertex
void HitThreshold(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                  double threshold, std::vector<int> * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
//this gives a list of all of the origins of the tpc objects
void GetOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::string> * tpco_origin_v);
//***************************************************************************
//***************************************************************************
//this function simply checks if the tpc object is a nue
void HasNue(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco, const bool _verbose);
//***************************************************************************
//***************************************************************************
//this function just counts if at least 1 tpc object passes the cuts
bool ValidTPCObjects(std::vector<int> * passed_tpco);
//***************************************************************************
//***************************************************************************
std::vector<int> TabulateOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco);
//***************************************************************************
//***************************************************************************
//modify this so it takes a string of the cut name so I only pass it a few variable at a time,
//then I can call this function several times later at the bottom
void PrintInfo(int mc_nue_cc_counter,
               int counter,
               int counter_nue_cc,
               int counter_nue_cc_mixed,
               int counter_cosmic,
               int counter_nue_nc,
               int counter_numu,
               int counter_unmatched,
               int counter_other_mixed,
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
void xsec_plot(bool _verbose, double genie_xsec, double xsec);
//***************************************************************************
//***************************************************************************
void PostCutPlots(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco, bool _verbose, TH2I * h_tracks_showers);

};

#endif