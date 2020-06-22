#ifndef VARIATION_OUTPUT_BKG_H
#define VARIATION_OUTPUT_BKG_H

// ------------------------------------------------
// Class description
// A class to get cut variables from the selection
// and add them to a root file for further inspection.
// The intended use of this is to investigate
// the effect of the various detector systematics.

// USAGE
// mode = "bnb" || "numi"
// plot_confg = "same" || "bnbnumi" || ""

// ./main.exe --var_mode_bkg <MC_Var_File> <mode> <plot_config>
// ------------------------------------------------

// Library and Other Includes
#include "../xsecAna/TpcObjectContainer.h"
#include "../xsecAna/ParticleContainer.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TKey.h"
#include "TLine.h"

// C++ includes
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <math.h>

#include "../xsecAna/LinkDef.h"



class variation_output_bkg {
    public:

    // Default constructor
    variation_output_bkg(){

        std::vector<double> bins = { 0, 15, 30, 45, 60, 75, 90}; 
        // std::vector<double> bins = { 0, 22.5, 45, 67.5, 90}; // Even smaller amount of bins
        double* bins_phi_wrapped = &bins[0];
        int n_bins = bins.size() - 1;

        TH1D_hist.resize(kTH1D_names_MAX);
        // Leading shower
        TH1D_hist.at(kldg_shwr_hits)             = new TH1D("h_ldg_shwr_hits","ldg_shwr_hits",                  50, 0, 600); 
        TH1D_hist.at(kldg_shwr_hits_WPlane)     = new TH1D("h_ldg_shwr_hits_WPlane","ldg_shwr_hits_WPlane",  30, 0, 1000);
        TH1D_hist.at(kldg_shwr_Open_Angle)        = new TH1D("h_ldg_shwr_Open_Angle","ldg_shwr_Open_Angle",      25, 0, 50);
        TH1D_hist.at(kldg_shwr_dEdx_WPlane)     = new TH1D("h_ldg_shwr_dEdx_WPlane","ldg_shwr_dEdx_WPlane",  40, 0, 10);
        TH1D_hist.at(kldg_shwr_HitPerLen)         = new TH1D("h_ldg_shwr_HitPerLen","ldg_shwr_HitPerLen",      20, 0, 20);
        TH1D_hist.at(kldg_shwr_Phi)             = new TH1D("h_ldg_shwr_Phi","ldg_shwr_Phi",                  12, -180, 180);
        TH1D_hist.at(kldg_shwr_Phi_wrapped)        = new TH1D("h_ldg_shwr_Phi_wrapped","ldg_shwr_Phi_wrapped",  n_bins, bins_phi_wrapped);
        TH1D_hist.at(kldg_shwr_Theta)             = new TH1D("h_ldg_shwr_Theta","ldg_shwr_Theta",              12, 0, 180);
        TH1D_hist.at(kldg_shwr_CTheta)             = new TH1D("h_ldg_shwr_CTheta","ldg_shwr_CTheta",              16, -1, 1);
        TH1D_hist.at(klong_Track_ldg_shwr)         = new TH1D("h_long_Track_ldg_shwr","long_Track_ldg_shwr",     20, 0, 3);
        
        // TPC
        TH1D_hist.at(ktpc_obj_vtx_x)             = new TH1D("h_tpc_obj_vtx_x","tpc_obj_vtx_x",                  20, 0, 260);
        TH1D_hist.at(ktpc_obj_vtx_y)             = new TH1D("h_tpc_obj_vtx_y","tpc_obj_vtx_y",                  20, -117, 117);
        TH1D_hist.at(ktpc_obj_vtx_z)             = new TH1D("h_tpc_obj_vtx_z","tpc_obj_vtx_z",                  40, 0, 1040);
        
        // Other
        TH1D_hist.at(ktotal_hits)                 = new TH1D("h_total_hits", "h_total_hits",                      50, 0, 600);
        TH1D_hist.at(kn_pfp)                     = new TH1D("h_n_pfp", "h_n_pfp",                              8, 0, 8);
        TH1D_hist.at(kn_pfp_50Hits)             = new TH1D("h_n_pfp_50Hits", "h_n_pfp_50Hits",                  8, 0, 8);
        TH1D_hist.at(kn_tracks)                 = new TH1D("h_n_tracks", "h_n_tracks",                          8, 0, 8);
        TH1D_hist.at(kn_tracks_50Hits)             = new TH1D("h_n_tracks_50Hits", "h_n_tracks_50Hits",          8, 0, 8);
        TH1D_hist.at(kn_showers)                 = new TH1D("h_n_showers", "h_n_showers",                      8, 0, 8);
        TH1D_hist.at(kn_showers_50Hits)         = new TH1D("h_n_showers_50Hits", "h_n_showers_50Hits",          8, 0, 8);
        TH1D_hist.at(ktrack_phi)                 = new TH1D("h_track_phi", "h_track_phi",                     12 , -180 ,180);
        TH1D_hist.at(kshower_phi)                 = new TH1D("h_shower_phi", "h_shower_phi",                     12 , -180 ,180); // Shower Phi
        TH1D_hist.at(kselected)                 = new TH1D("h_selected", "h_selected",                         1 , 0 , 1);

        TH1D_hist.at(kshower_Nu_vtx_Dist)        = new TH1D("h_shower_Nu_vtx_Dist","h_shower_Nu_vtx_Dist",     20, 0, 20);
        TH1D_hist.at(ktrack_Nu_vtx_Dist)        = new TH1D("h_track_Nu_vtx_Dist","h_track_Nu_vtx_Dist",          20, 0, 20);

        // Largest Flash
        TH1D_hist.at(klargest_flash_y)            = new TH1D("h_largest_flash_y", "h_largest_flash_y",         60, -40, 40);
        TH1D_hist.at(klargest_flash_z)            = new TH1D("h_largest_flash_z", "h_largest_flash_z",         125, 0, 1000);
        TH1D_hist.at(klargest_flash_time)        = new TH1D("h_largest_flash_time", "h_largest_flash_time",    50, 0, 20);
        TH1D_hist.at(klargest_flash_pe)            = new TH1D("h_largest_flash_pe", "h_largest_flash_pe",         30, 0, 6000);
        TH1D_hist.at(kFlash_TPCObj_Dist)        = new TH1D("h_Flash_TPCObj_Dist", "h_Flash_TPCObj_Dist",     50, 0, 200);     // Largest flash to TPC OBj Vtx Dist

        TH1D_hist.at(kshower_phi_pi0)             = new TH1D("h_shower_phi_pi0",        "h_shower_phi_pi0",        12 , -180 ,180); // Shower Phi pi0
        TH1D_hist.at(kshower_phi_bkg_cosmic)     = new TH1D("h_shower_phi_bkg_cosmic", "h_shower_phi_bkg_cosmic", 12 , -180 ,180); // Shower Phi bkg electrons
        TH1D_hist.at(kshower_phi_other)         = new TH1D("h_shower_phi_other",      "h_shower_phi_other",      12 , -180 ,180); // Shower Phi
        
        TH1D_hist.at(kshower_phi_pi0_wrapped)         = new TH1D("h_shower_phi_pi0_wrapped",        "h_shower_phi_pi0_wrapped",        n_bins, bins_phi_wrapped); // Shower Phi pi0 wrapped
        TH1D_hist.at(kshower_phi_bkg_cosmic_wrapped)  = new TH1D("h_shower_phi_bkg_cosmic_wrapped", "h_shower_phi_bkg_cosmic_wrapped", n_bins, bins_phi_wrapped); // Shower Phi bkg electrons wrapped
        TH1D_hist.at(kshower_phi_other_wrapped)       = new TH1D("h_shower_phi_other_wrapped",      "h_shower_phi_other_wrapped",      n_bins, bins_phi_wrapped); // Shower Phi wrapped

        TH1D_hist.at(kshower_E_pi0)               = new TH1D("h_shower_E_pi0",        "h_shower_E_pi0",         7 , 0, 3 ); // Shower E pi0
        TH1D_hist.at(kshower_E_bkg_cosmic)        = new TH1D("h_shower_E_bkg_cosmic", "h_shower_E_bkg_cosmic",  7 , 0, 5 ); // Shower E bkg cosmic
        TH1D_hist.at(kshower_E_other)             = new TH1D("h_shower_E_other",      "h_shower_E_other",       7 , 0, 3 ); // Shower E other
        TH1D_hist.at(kshower_E)                   = new TH1D("h_shower_E",            "h_shower_E",             7 , 0, 3 ); // Shower E

        TH1D_hist.at(kshower_Theta_pi0)               = new TH1D("h_shower_Theta_pi0",        "h_shower_Theta_pi0",         12, 0, 180 ); // Shower Theta pi0
        TH1D_hist.at(kshower_Theta_bkg_cosmic)        = new TH1D("h_shower_Theta_bkg_cosmic", "h_shower_Theta_bkg_cosmic",  12, 0, 180 ); // Shower Theta bkg cosmic
        TH1D_hist.at(kshower_Theta_other)             = new TH1D("h_shower_Theta_other",      "h_shower_Theta_other",       12, 0, 180 ); // Shower Theta other

        // The unselected histograms
        TH1D_hist.at(kshower_phi_pi0_unselected)           = new TH1D("h_shower_phi_pi0_unselected",        "h_shower_phi_pi0_unselected",        12 , -180 ,180); // Shower Phi pi0
        TH1D_hist.at(kshower_phi_bkg_cosmic_unselected)    = new TH1D("h_shower_phi_bkg_cosmic_unselected", "h_shower_phi_bkg_cosmic_unselected", 12 , -180 ,180); // Shower Phi bkg electrons
        TH1D_hist.at(kshower_phi_other_unselected)         = new TH1D("h_shower_phi_other_unselected",      "h_shower_phi_other_unselected",      12 , -180 ,180); // Shower Phi
        
        TH1D_hist.at(kshower_phi_pi0_wrapped_unselected)         = new TH1D("h_shower_phi_pi0_wrapped_unselected",        "h_shower_phi_pi0_wrapped_unselected",        n_bins, bins_phi_wrapped); // Shower Phi pi0 wrapped
        TH1D_hist.at(kshower_phi_bkg_cosmic_wrapped_unselected)  = new TH1D("h_shower_phi_bkg_cosmic_wrapped_unselected", "h_shower_phi_bkg_cosmic_wrapped_unselected", n_bins, bins_phi_wrapped); // Shower Phi bkg electrons wrapped
        TH1D_hist.at(kshower_phi_other_wrapped_unselected)       = new TH1D("h_shower_phi_other_wrapped_unselected",      "h_shower_phi_other_wrapped_unselected",      n_bins, bins_phi_wrapped); // Shower Phi wrapped

        TH1D_hist.at(kshower_E_pi0_unselected)               = new TH1D("h_shower_E_pi0_unselected",        "h_shower_E_pi0_unselected",         7 , 0, 3 ); // Shower E pi0
        TH1D_hist.at(kshower_E_bkg_cosmic_unselected)        = new TH1D("h_shower_E_bkg_cosmic_unselected", "h_shower_E_bkg_cosmic_unselected",  7 , 0, 5 ); // Shower E bkg cosmic
        TH1D_hist.at(kshower_E_other_unselected)             = new TH1D("h_shower_E_other_unselected",      "h_shower_E_other_unselected",       7 , 0, 3 ); // Shower E other
        TH1D_hist.at(kshower_E_unselected)                   = new TH1D("h_shower_E_unselected",            "h_shower_E_unselected",             7 , 0, 3 ); // Shower E

        TH1D_hist.at(kshower_Theta_pi0_unselected)           = new TH1D("h_shower_Theta_pi0_unselected",        "h_shower_Theta_pi0_unselected",         12, 0, 180 ); // Shower Theta pi0
        TH1D_hist.at(kshower_Theta_bkg_cosmic_unselected)    = new TH1D("h_shower_Theta_bkg_cosmic_unselected", "h_shower_Theta_bkg_cosmic_unselected",  12, 0, 180 ); // Shower Theta bkg cosmic
        TH1D_hist.at(kshower_Theta_other_unselected)         = new TH1D("h_shower_Theta_other_unselected",      "h_shower_Theta_other_unselected",       12, 0, 180 ); // Shower Theta other

        TH1D_hist.at(kldg_shwr_Phi_unselected)               = new TH1D("h_ldg_shwr_Phi_unselected","ldg_shwr_Phi_unselected",                  12, -180, 180);
        TH1D_hist.at(kldg_shwr_Phi_wrapped_unselected)       = new TH1D("h_ldg_shwr_Phi_wrapped_unselected","ldg_shwr_Phi_wrapped_unselected",  n_bins, bins_phi_wrapped);
        TH1D_hist.at(kldg_shwr_Theta_unselected)             = new TH1D("h_ldg_shwr_Theta_unselected","ldg_shwr_Theta_unselected",              12, 0, 180);

        TH1D_hist.at(kpre_hit_threshold_cut)                = new TH1D("h_shwr_hit_threshold_precut",            "h_shwr_hit_threshold_precut;Num MC Hits; Entries", 100, 0, 500);
        TH1D_hist.at(kpre_hit_threshold_collection_cut)     = new TH1D("h_shwr_hit_threshold_collection_precut", "h_shwr_hit_threshold_collection_precut;Num MC Hits Collection Plane; Entries", 100, 0, 500);
        TH1D_hist.at(kpre_open_angle_cut)                   = new TH1D("h_ldg_shwr_Open_Angle_precut",           "ldg_shwr_Open_Angle_precut;Open Angle [degrees]; Entries",      25, 0, 50);
        TH1D_hist.at(kpre_dedx_cut)                         = new TH1D("h_ldg_shwr_dEdx_WPlane_precut",          "ldg_shwr_dEdx_WPlane_precut;dEdx [MeV/cm]; Entries",  40, 0, 10);
        TH1D_hist.at(kpost_hit_threshold_cut)               = new TH1D("h_shwr_hit_threshold_postcut",           "h_shwr_hit_threshold_postcut;Num MC Hits; Entries", 100, 0 ,500);
        TH1D_hist.at(kpost_hit_threshold_collection_cut)    = new TH1D("h_shwr_hit_threshold_collection_postcut","h_shwr_hit_threshold_collection_postcut;Num MC Hits Collection Plane; Entries", 100, 0 ,500);
        TH1D_hist.at(kpost_open_angle_cut)                  = new TH1D("h_ldg_shwr_Open_Angle_postcut",          "ldg_shwr_Open_Angle_postcut;Open Angle [degrees]; Entries",      25, 0, 50);
        TH1D_hist.at(kpost_dedx_cut)                        = new TH1D("h_ldg_shwr_dEdx_WPlane_postcut",         "ldg_shwr_dEdx_WPlane_postcut;dEdx [MeV/cm]; Entries",  40, 0, 10);


        // Do the same but for the weighted histograms
        TH1D_hist_weighted.resize(kTH1D_names_weighted_MAX);

        // Loop over the variations an make weighted histograms for each variation
        for (int p=0; p < TH1D_hist_weighted.size(); p++){
            TH1D_hist_weighted.at(p).resize(bnbvars.size()); // There are 18 variations in total I think...
        }

        // A loop over the variations
        for (unsigned int k=0; k < TH1D_hist_weighted.at(0).size(); k++){
            // Leading shower
            TH1D_hist_weighted.at(kldg_shwr_Phi_w).at(k)                   = new TH1D(Form("h_ldg_shwr_Phi_%s_w", bnbvars.at(k).c_str()),                  Form("ldg_shwr_Phi %s", bnbvars.at(k).c_str()),                   12, -180, 180);
            TH1D_hist_weighted.at(kldg_shwr_Phi_wrapped_w).at(k)           = new TH1D(Form("h_ldg_shwr_Phi_wrapped_%s_w", bnbvars.at(k).c_str()),          Form("ldg_shwr_Phi_wrapped %s", bnbvars.at(k).c_str()),           n_bins, bins_phi_wrapped);
            TH1D_hist_weighted.at(kldg_shwr_Theta_w).at(k)                 = new TH1D(Form("h_ldg_shwr_Theta_%s_w", bnbvars.at(k).c_str()),                Form("ldg_shwr_Theta %s", bnbvars.at(k).c_str()),                 12, 0, 180);
            TH1D_hist_weighted.at(kselected_w).at(k)                       = new TH1D(Form("h_selected_%s_w", bnbvars.at(k).c_str()),                      Form("h_selected %s", bnbvars.at(k).c_str()),                     1 , 0 , 1);
            TH1D_hist_weighted.at(kshower_phi_pi0_w).at(k)                 = new TH1D(Form("h_shower_phi_pi0_%s_w", bnbvars.at(k).c_str()),                Form("h_shower_phi_pi0 %s", bnbvars.at(k).c_str()),               12 , -180 ,180); // Shower Phi pi0
            TH1D_hist_weighted.at(kshower_phi_bkg_cosmic_w).at(k)          = new TH1D(Form("h_shower_phi_bkg_cosmic_%s_w", bnbvars.at(k).c_str()),         Form("h_shower_phi_bkg_cosmic %s", bnbvars.at(k).c_str()),        12 , -180 ,180); // Shower Phi bkg electrons
            TH1D_hist_weighted.at(kshower_phi_other_w).at(k)               = new TH1D(Form("h_shower_phi_other_%s_w", bnbvars.at(k).c_str()),              Form("h_shower_phi_other %s", bnbvars.at(k).c_str()),             12 , -180 ,180); // Shower Phi
            TH1D_hist_weighted.at(kshower_phi_pi0_wrapped_w).at(k)         = new TH1D(Form("h_shower_phi_pi0_wrapped_%s_w", bnbvars.at(k).c_str()),        Form("h_shower_phi_pi0_wrapped %s", bnbvars.at(k).c_str()),       n_bins, bins_phi_wrapped); // Shower Phi pi0 wrapped
            TH1D_hist_weighted.at(kshower_phi_bkg_cosmic_wrapped_w).at(k)  = new TH1D(Form("h_shower_phi_bkg_cosmic_wrapped_%s_w", bnbvars.at(k).c_str()), Form("h_shower_phi_bkg_cosmic_wrapped %s", bnbvars.at(k).c_str()),n_bins, bins_phi_wrapped); // Shower Phi bkg electrons wrapped
            TH1D_hist_weighted.at(kshower_phi_other_wrapped_w).at(k)       = new TH1D(Form("h_shower_phi_other_wrapped_%s_w", bnbvars.at(k).c_str()),      Form("h_shower_phi_other_wrapped %s", bnbvars.at(k).c_str()),     n_bins, bins_phi_wrapped); // Shower Phi wrapped
            TH1D_hist_weighted.at(kshower_E_pi0_w).at(k)                   = new TH1D(Form("h_shower_E_pi0_%s_w", bnbvars.at(k).c_str()),                  Form("h_shower_E_pi0 %s", bnbvars.at(k).c_str()),                 7 , 0, 3 ); // Shower E pi0
            TH1D_hist_weighted.at(kshower_E_bkg_cosmic_w).at(k)            = new TH1D(Form("h_shower_E_bkg_cosmic_%s_w", bnbvars.at(k).c_str()),           Form("h_shower_E_bkg_cosmic %s", bnbvars.at(k).c_str()),          7 , 0, 5 ); // Shower E bkg cosmic
            TH1D_hist_weighted.at(kshower_E_other_w).at(k)                 = new TH1D(Form("h_shower_E_other_%s_w", bnbvars.at(k).c_str()),                Form("h_shower_E_other %s", bnbvars.at(k).c_str()),               7 , 0, 3 ); // Shower E other
            TH1D_hist_weighted.at(kshower_E_w).at(k)                       = new TH1D(Form("h_shower_E_%s_w", bnbvars.at(k).c_str()),                      Form("h_shower_E %s", bnbvars.at(k).c_str()),                     7 , 0, 3 ); // Shower E
            TH1D_hist_weighted.at(kshower_Theta_pi0_w).at(k)               = new TH1D(Form("h_shower_Theta_pi0_%s_w", bnbvars.at(k).c_str()),              Form("h_shower_Theta_pi0 %s", bnbvars.at(k).c_str()),             12, 0, 180 ); // Shower Theta pi0
            TH1D_hist_weighted.at(kshower_Theta_bkg_cosmic_w).at(k)        = new TH1D(Form("h_shower_Theta_bkg_cosmic_%s_w", bnbvars.at(k).c_str()),       Form("h_shower_Theta_bkg_cosmic %s", bnbvars.at(k).c_str()),      12, 0, 180 ); // Shower Theta bkg cosmic
            TH1D_hist_weighted.at(kshower_Theta_other_w).at(k)             = new TH1D(Form("h_shower_Theta_other_%s_w", bnbvars.at(k).c_str()),            Form("h_shower_Theta_other %s", bnbvars.at(k).c_str()),           12, 0, 180 ); // Shower Theta other

        }
        
        
        // 2D Histograms
        TH2D_hist.resize(kTH2D_names_MAX);
        TH2D_hist.at(kEBkg_Theta)       = new TH2D("h_EBkg_Theta",       "EBkg_Theta",       7, 0, 3, 12,    0, 180);
        TH2D_hist.at(kEBkg_Phi)         = new TH2D("h_EBkg_Phi",         "EBkg_Phi",         7, 0, 3, 12, -180, 180); 
    
        TH2D_hist.at(kEBkg_pi0_Theta)   = new TH2D("h_EBkg_pi0_Theta",   "EBkg_pi0_Theta",   7, 0, 3, 12,    0, 180);
        TH2D_hist.at(kEBkg_pi0_Phi)     = new TH2D("h_EBkg_pi0_Phi",     "EBkg_pi0_Phi",     7, 0, 3, 12, -180, 180);
        
        TH2D_hist.at(kEBkg_cosmic_Theta)= new TH2D("h_EBkg_cosmic_Theta","EBkg_cosmic_Theta",7, 0, 3, 12,    0, 180);
        TH2D_hist.at(kEBkg_cosmic_Phi)  = new TH2D("h_EBkg_cosmic_Phi",  "EBkg_cosmic_Phi",  7, 0, 3, 12, -180, 180);
        
        TH2D_hist.at(kEBkg_other_Theta) = new TH2D("h_EBkg_other_Theta", "EBkg_other_Theta", 7, 0, 3, 12,    0, 180);
        TH2D_hist.at(kEBkg_other_Phi)   = new TH2D("h_EBkg_other_Phi",   "EBkg_other_Phi",   7, 0, 3, 12, -180, 180);

        TH2D_hist.at(kEBkg_Phi_wrapped)         = new TH2D("h_EBkg_Phi_wrapped",         "EBkg_Phi_wrapped",         7, 0, 3, 5, 0, 90); 
        TH2D_hist.at(kEBkg_pi0_Phi_wrapped)     = new TH2D("h_EBkg_pi0_Phi_wrapped",     "EBkg_pi0_Phi_wrapped",     7, 0, 3, 5, 0, 90);
        TH2D_hist.at(kEBkg_cosmic_Phi_wrapped)  = new TH2D("h_EBkg_cosmic_Phi_wrapped",  "EBkg_cosmic_Phi_wrapped",  7, 0, 3, 5, 0, 90);
        TH2D_hist.at(kEBkg_other_Phi_wrapped)   = new TH2D("h_EBkg_other_Phi_wrapped",   "EBkg_other_Phi_wrapped",   7, 0, 3, 5, 0, 90);

        TH2D_hist.at(kThetaBkg_Phi_wrapped)         = new TH2D("h_ThetaBkg_Phi_wrapped",         "ThetaBkg_Phi_wrapped",        12, 0, 180, n_bins, bins_phi_wrapped); 
        TH2D_hist.at(kThetaBkg_pi0_Phi_wrapped)     = new TH2D("h_ThetaBkg_pi0_Phi_wrapped",     "ThetaBkg_pi0_Phi_wrapped",    12, 0, 180, n_bins, bins_phi_wrapped);
        TH2D_hist.at(kThetaBkg_cosmic_Phi_wrapped)  = new TH2D("h_ThetaBkg_cosmic_Phi_wrapped",  "ThetaBkg_cosmic_Phi_wrapped", 12, 0, 180, n_bins, bins_phi_wrapped);
        TH2D_hist.at(kThetaBkg_other_Phi_wrapped)   = new TH2D("h_ThetaBkg_other_Phi_wrapped",   "ThetaBkg_other_Phi_wrapped",  12, 0, 180, n_bins, bins_phi_wrapped);

        TH2D_hist.at(kThetaBkg_Phi)         = new TH2D("h_ThetaBkg_Phi",         "ThetaBkg_Phi",        12, 0, 180, 12, -180, 180); 
        TH2D_hist.at(kThetaBkg_pi0_Phi)     = new TH2D("h_ThetaBkg_pi0_Phi",     "ThetaBkg_pi0_Phi",    12, 0, 180, 12, -180, 180);
        TH2D_hist.at(kThetaBkg_cosmic_Phi)  = new TH2D("h_ThetaBkg_cosmic_Phi",  "ThetaBkg_cosmic_Phi", 12, 0, 180, 12, -180, 180);
        TH2D_hist.at(kThetaBkg_other_Phi)   = new TH2D("h_ThetaBkg_other_Phi",   "ThetaBkg_other_Phi",  12, 0, 180, 12, -180, 180);

        // unselected histograms
        TH2D_hist.at(kEBkg_Theta_unselected)       = new TH2D("h_EBkg_Theta_unselected",       "EBkg_Theta_unselected",       7, 0, 3, 12,    0, 180);
        TH2D_hist.at(kEBkg_Phi_unselected)         = new TH2D("h_EBkg_Phi_unselected",         "EBkg_Phi_unselected",         7, 0, 3, 12, -180, 180); 
    
        TH2D_hist.at(kEBkg_pi0_Theta_unselected)   = new TH2D("h_EBkg_pi0_Theta_unselected",   "EBkg_pi0_Theta_unselected",   7, 0, 3, 12,    0, 180);
        TH2D_hist.at(kEBkg_pi0_Phi_unselected)     = new TH2D("h_EBkg_pi0_Phi_unselected",     "EBkg_pi0_Phi_unselected",     7, 0, 3, 12, -180, 180);
        
        TH2D_hist.at(kEBkg_cosmic_Theta_unselected)= new TH2D("h_EBkg_cosmic_Theta_unselected","EBkg_cosmic_Theta_unselected",7, 0, 3, 12,    0, 180);
        TH2D_hist.at(kEBkg_cosmic_Phi_unselected)  = new TH2D("h_EBkg_cosmic_Phi_unselected",  "EBkg_cosmic_Phi_unselected",  7, 0, 3, 12, -180, 180);
        
        TH2D_hist.at(kEBkg_other_Theta_unselected) = new TH2D("h_EBkg_other_Theta_unselected", "EBkg_other_Theta_unselected", 7, 0, 3, 12,    0, 180);
        TH2D_hist.at(kEBkg_other_Phi_unselected)   = new TH2D("h_EBkg_other_Phi_unselected",   "EBkg_other_Phi_unselected",   7, 0, 3, 12, -180, 180);

        TH2D_hist.at(kEBkg_Phi_wrapped_unselected)         = new TH2D("h_EBkg_Phi_wrapped_unselected",         "EBkg_Phi_wrapped_unselected",         7, 0, 3, 5, 0, 90); 
        TH2D_hist.at(kEBkg_pi0_Phi_wrapped_unselected)     = new TH2D("h_EBkg_pi0_Phi_wrapped_unselected",     "EBkg_pi0_Phi_wrapped_unselected",     7, 0, 3, 5, 0, 90);
        TH2D_hist.at(kEBkg_cosmic_Phi_wrapped_unselected)  = new TH2D("h_EBkg_cosmic_Phi_wrapped_unselected",  "EBkg_cosmic_Phi_wrapped_unselected",  7, 0, 3, 5, 0, 90);
        TH2D_hist.at(kEBkg_other_Phi_wrapped_unselected)   = new TH2D("h_EBkg_other_Phi_wrapped_unselected",   "EBkg_other_Phi_wrapped_unselected",   7, 0, 3, 5, 0, 90);

        TH2D_hist.at(kThetaBkg_Phi_wrapped_unselected)         = new TH2D("h_ThetaBkg_Phi_wrapped_unselected",         "ThetaBkg_Phi_wrapped_unselected",        12, 0, 180, n_bins, bins_phi_wrapped); 
        TH2D_hist.at(kThetaBkg_pi0_Phi_wrapped_unselected)     = new TH2D("h_ThetaBkg_pi0_Phi_wrapped_unselected",     "ThetaBkg_pi0_Phi_wrapped_unselected",    12, 0, 180, n_bins, bins_phi_wrapped);
        TH2D_hist.at(kThetaBkg_cosmic_Phi_wrapped_unselected)  = new TH2D("h_ThetaBkg_cosmic_Phi_wrapped_unselected",  "ThetaBkg_cosmic_Phi_wrapped_unselected", 12, 0, 180, n_bins, bins_phi_wrapped);
        TH2D_hist.at(kThetaBkg_other_Phi_wrapped_unselected)   = new TH2D("h_ThetaBkg_other_Phi_wrapped_unselected",   "ThetaBkg_other_Phi_wrapped_unselected",  12, 0, 180, n_bins, bins_phi_wrapped);

        TH2D_hist.at(kThetaBkg_Phi_unselected)         = new TH2D("h_ThetaBkg_Phi_unselected",         "ThetaBkg_Phi_unselected",        12, 0, 180, 12, -180, 180); 
        TH2D_hist.at(kThetaBkg_pi0_Phi_unselected)     = new TH2D("h_ThetaBkg_pi0_Phi_unselected",     "ThetaBkg_pi0_Phi_unselected",    12, 0, 180, 12, -180, 180);
        TH2D_hist.at(kThetaBkg_cosmic_Phi_unselected)  = new TH2D("h_ThetaBkg_cosmic_Phi_unselected",  "ThetaBkg_cosmic_Phi_unselected", 12, 0, 180, 12, -180, 180);
        TH2D_hist.at(kThetaBkg_other_Phi_unselected)   = new TH2D("h_ThetaBkg_other_Phi_unselected",   "ThetaBkg_other_Phi_unselected",  12, 0, 180, 12, -180, 180);

        TBranch *bmc_phi       = VariableTree ->Branch("mc_phi",       &mc_Phi,       "mc_Phi");
        TBranch *bmc_phi_pi0   = VariableTree ->Branch("mc_phi_pi0",   &mc_Phi_pi0,   "mc_Phi_pi0");
        TBranch *bmc_phi_other = VariableTree ->Branch("mc_phi_other", &mc_Phi_other, "mc_Phi_other");
       
    }
    
    // ----------------------
    //   Main Function
    // ----------------------
    void run_var(const char * file1, TString mode, const std::vector<double> _config, TString plot_config);

    // ----------------------
    //   General Functions
    // ----------------------
    int GetLeadingShowerIndex(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj);     // Returns the index of the leading shower
    double GetLongestTrackLength(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj);  // Returns the length of the longest track
    void DrawTH1D(TH1D* h, double POT_Scaling);                                                         // Function that draws a TH1D histogram
    void DrawTH2D(TH2D* h, double POT_Scaling);                                                         // Function that draws a TH2D histogram
    double GetPOT(const char * _file1);                                                                 // Gets the POT stored in an external file
    void PlotVariatons(TFile* f_var_out, TString mode);                                             // Plots the variation files on the same plot
    void PlotVariatonsNuMIBNB();                                                                        // Plotting function to compare the NuMI and BNB CV variations
    std::vector<std::string> GrabDirs(TFile* f_var_out);                                                // Grabs the directories in the file
    void DrawTH1D_SAME(TH1D* hist, std::string variation, TLegend* legend, std::string histname);       // Function that draws a TH1D histogram for the same plot
    void DrawTH2D_SAME(TH2D* hist, std::string variation, std::string histname);                        // Function that draws a TH2D histogram for the same plot
    void DrawTH1D_Ratio(TH1D* hist,  std::string variation, TLegend* legend, std::string histname);                       // As above but for drawing the ratios
    
    void GetNumber_Track_Shower(const int n_pfp, int n_tpc_obj,                                         // Utility function to get the number of tracks and showers
                                     xsecAna::TPCObjectContainer tpc_obj, int &n_showers, int &n_tracks,
                                     int &n_pfp_50Hits, int &n_tracks_50Hits, int &n_showers_50Hits);
    
    double pfp_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,                       // Calculates the pfp to nu vertex distance
                                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z);
    
    bool in_fv(double x, double y, double z, std::vector<double> fv_boundary_v);
    std::pair<std::string, int> TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool true_in_tpc, bool has_pi0);
    std::string Background_Classifier(int mc_pdg, std::string tpc_obj_classfication);
    double WrapPhi(double phi);                                                                         // A function to wrap phi from 360 to 90 degrees to increse statistics
    void GenerateWeightHistograms(); // To generate a file with the histogram weights
    void WeightBNBVar(xsecAna::TPCObjectContainer tpc_obj, bool bool_sig, const int leading_shower_index, std::pair<std::string, int> tpc_classification, std::string dirname); // Function to weight the numi CV
    void CompareWeightedHistograms(); // Function to make the NuMI variaions and weighted NuMI CV variations
    void CompareWeightedDrawSpecs(TH1D* hist, std::string weighted_str, std::string variation, TLegend* legend, std::string histname); // Function to draw Weighted NuMI Histograms
    void GetBNBBkgWeight(double theta, double phi, double phi_wrapped, std::string bkg_class, double &weight_all, double &weight_indiv, std::string dirname ); // Reweight the backgrounds based on a 2D histogram in theta and phi
    void FillHistogramCuts( xsecAna::TPCObjectContainer tpc_obj, std::string histname, bool bool_sig); // Function to fill histograms before and after selection cuts

    // ----------------------
    // Flash Functions
    // ----------------------
    std::vector<std::vector<double>> GetLargestFlashVector(TFile* f, double flash_time_start, double flash_time_end); // Function to resize opical entries to same size of events and get largest flash vector
    bool flash_in_time(double flash_time, double flash_start, double flash_end);                                       // Decides whether flash is in time or not
    bool flash_pe(int flash_pe, int flash_pe_threshold);                                                               // Decides whether flash has sufficient PE
    double Flash_TPCObj_vtx_Dist(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z);         // Returns the 2D distance of the flash to TPC OBj Vertex
 
    // ----------------------
    //   Cut Functions
    // ----------------------
    void FlashinTime_FlashPE(TFile* f, double flash_start_time, double flash_end_time, std::vector<bool> &flash_cuts_pass_vec, TString mode );
    bool HasNue(xsecAna::TPCObjectContainer tpc_obj, const int n_pfp );
    bool opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance);
    bool flashRecoVtxDist(std::vector< double > largest_flash_v, double tolerance, const double tpc_vtx_x, const double tpc_vtx_y, const double tpc_vtx_z);
    bool VtxNuDistance(xsecAna::TPCObjectContainer tpc_obj,int pfp_pdg_type , double tolerance);
    bool HitThreshold(xsecAna::TPCObjectContainer tpc_obj, double threshold, bool useCollection);
    bool OpenAngleCut(xsecAna::TPCObjectContainer tpc_obj, const std::vector<double> tolerance_open_angle);
    bool dEdxCut( xsecAna::TPCObjectContainer tpc_obj, const double tolerance_dedx_min, const double tolerance_dedx_max);
    bool SecondaryShowersDistCut(xsecAna::TPCObjectContainer tpc_obj, const double dist_tolerance);
    bool HitLengthRatioCut(const double pfp_hits_length_tolerance, xsecAna::TPCObjectContainer tpc_obj);
    bool LongestTrackLeadingShowerCut(const double ratio_tolerance, xsecAna::TPCObjectContainer tpc_obj);
    bool IsContained(std::vector<double> track_start, std::vector<double> track_end, std::vector<double> fv_boundary_v);
    bool ContainedTracksCut(std::vector<double> fv_boundary_v, xsecAna::TPCObjectContainer tpc_obj);

    // ----------------------
    //     Reco Variables
    // ----------------------
    double ldg_shwr_hits{0};             // Leading shower hits (all planes)
    double ldg_shwr_hits_WPlane{0};      // Leading Shower Collection plane hits
    double lldg_shwr_Open_Angle{0};      // Leading Shower Opening angle of shower
    double ldg_shwr_dEdx_WPlane{0};      // Leading shower dE/dx at Collection Plane
    double ldg_shwr_HitPerLen{0};        // Leading shower hits per length
    double ldg_shwr_Phi{0};              // Leading shower phi
    double ldg_shwr_Theta{0};            // Leading shower theta
    double ldg_shwr_CTheta{0};           // Leading shower cos theta
    double long_Track_ldg_shwr{0};       // Longest track / leading shower length
    
    // Leading shower momentum
    int n_pfp{0};                       // Number of pfp in TPC Obj
    int n_pfp_50Hits{0};                // Number of pfp in TPC Obj with > 50 Hits
    int n_tracks{0};                    // Number of Tracks in TPC Obj
    int n_tracks_50Hits{0};             // Number of Tracks in TPC Obj with > 50 Hits
    int n_showers{0};                   // Number of Showers in TPC Obj
    int n_showers_50Hits{0};            // Number of showers in TPC Obj with > 50 Hits
    double track_phi{0};                // Track Phi
    int nue_cc_counter{0};              // Number of reco nue + nuebar
    int other_counter{0};               // Counter for not numu/nue
    int sig_counter{0};                 // Number reco signal events
    int bkg_counter{0};                 // Number reco background events
    int pi0_counter{0};                 // Number of pi0 background events that are selected
    int cosmic_counter{0};              // Number of cosmic background events that are selected
    int other_bkg_counter{0};           // Number of other background events that are selected
    
    
    double tpc_obj_vtx_x{0}, tpc_obj_vtx_y{0}, tpc_obj_vtx_z{0}; // TPCObj Vertex X, Y, Z

    // ----------------------
    //    Cut Variables
    // ----------------------
    double _x1;
    double _x2;
    double _y1;
    double _y2;
    double _z1;
    double _z2;
    double flash_pe_threshold;
    double flash_time_start;
    double flash_time_end;
    double tolerance;
    double shwr_nue_tolerance;
    double trk_nue_tolerance;
    double shwr_hit_threshold;
    double shwr_hit_threshold_collection;
    double tolerance_open_angle_min;
    double tolerance_open_angle_max;
    double tolerance_dedx_min;
    double tolerance_dedx_max;
    double dist_tolerance;
    double pfp_hits_length_tolerance;
    double ratio_tolerance;
    bool detector_variations;

    // ----------------------
    //    Cut Counters
    // ----------------------
    int counter_FlashinTime_FlashPE{0};
    int counter_HasNue{0};
    int counter_inFV{0};
    int counter_FlashRecoVtxDist{0};
    int counter_VtxNuDist{0};
    int counter_VtxTrackNuDist{0};
    int counter_HitThresh{0};
    int counter_HitThreshW{0};
    int counter_OpenAngle{0};
    int counter_dEdx{0};
    int counter_SecondaryShowers{0};
    int counter_HitLenghtRatio{0};
    int counter_LongestTrackLeadingShower{0};
    int counter_ContainedTrack{0};


    // ----------------------
    //  MC Tree Variables
    // ----------------------
    int mc_nue_cc_counter = 0;
    int mc_nue_nc_counter = 0;
    int mc_nue_cc_counter_bar = 0;
    int mc_nue_nc_counter_bar = 0;
    
    double mc_nu_energy = 0;
    int mc_nu_id = -1;
    
    double mc_nu_vtx_x = -999;
    double mc_nu_vtx_y = -999;
    double mc_nu_vtx_z = -999;
    
    double mc_nu_dir_x = -999;
    double mc_nu_dir_y = -999;
    double mc_nu_dir_z = -999;
    
    double mc_ele_dir_x = -999;
    double mc_ele_dir_y = -999;
    double mc_ele_dir_z = -999;
    
    double mc_ele_energy = 0;
    double mc_ele_momentum = 0;
    
    bool has_pi0 = false;
    double mc_nu_time = -1;
    int mc_nu_num_particles = 0;
    int mc_nu_num_charged_particles = 0;

    // Weighted backgrounds
    double weight_all{0.0}, weight_indiv{0.0};

    // ----------------------
    //     Other Variables
    // ----------------------
    // vector containing all the signal type events
    std::vector<std::string> signal_modes = {"nue_cc_qe",     "nue_cc_res",     "nue_cc_dis",     "nue_cc_coh",     "nue_cc_mec",
                                             "nue_bar_cc_qe", "nue_bar_cc_res", "nue_bar_cc_dis", "nue_bar_cc_coh", "nue_bar_cc_mec" };

    // ----------------------
    //      Histograms
    // ----------------------
    std::vector<TH1D*> TH1D_hist;
    enum TH1D_names{ ktotal_hits,             kldg_shwr_hits,             kldg_shwr_hits_WPlane,
                     kldg_shwr_Open_Angle,    kldg_shwr_dEdx_WPlane,      kldg_shwr_HitPerLen,
                     kldg_shwr_Phi,           kldg_shwr_Phi_wrapped,      kldg_shwr_Theta,               kldg_shwr_CTheta,
                     klong_Track_ldg_shwr,    ktpc_obj_vtx_x,             ktpc_obj_vtx_y,                ktpc_obj_vtx_z,
                     kn_pfp, kn_pfp_50Hits,   kn_tracks,                  kn_tracks_50Hits,              kn_showers,
                     kn_showers_50Hits,       ktrack_phi,                 kshower_phi, klargest_flash_y, klargest_flash_z,
                     klargest_flash_time,     klargest_flash_pe,          kFlash_TPCObj_Dist,
                     kshower_Nu_vtx_Dist,     ktrack_Nu_vtx_Dist,         kselected,
                     kshower_phi_pi0,         kshower_phi_bkg_cosmic,         kshower_phi_other,             kshower_E, 
                     kshower_phi_pi0_wrapped, kshower_phi_bkg_cosmic_wrapped, kshower_phi_other_wrapped,
                     kshower_E_pi0,           kshower_E_bkg_cosmic,           kshower_E_other,
                     kshower_Theta_pi0,       kshower_Theta_bkg_cosmic,       kshower_Theta_other,
                     kshower_phi_pi0_unselected,         kshower_phi_bkg_cosmic_unselected,         kshower_phi_other_unselected,             kshower_E_unselected, 
                     kshower_phi_pi0_wrapped_unselected, kshower_phi_bkg_cosmic_wrapped_unselected, kshower_phi_other_wrapped_unselected,
                     kshower_E_pi0_unselected,           kshower_E_bkg_cosmic_unselected,           kshower_E_other_unselected,
                     kshower_Theta_pi0_unselected,       kshower_Theta_bkg_cosmic_unselected,       kshower_Theta_other_unselected,
                     kldg_shwr_Phi_unselected, kldg_shwr_Phi_wrapped_unselected, kldg_shwr_Theta_unselected,
                     kpre_hit_threshold_cut, kpre_hit_threshold_collection_cut,kpre_open_angle_cut, kpre_dedx_cut,
                     kpost_hit_threshold_cut, kpost_hit_threshold_collection_cut, kpost_open_angle_cut, kpost_dedx_cut,
                     kTH1D_names_MAX};

    
    
    // 2D Histograms -- histograms now created in the default constructor
    std::vector<TH2D*> TH2D_hist;
    enum TH2D_names{ kEBkg_Theta,        kEBkg_Phi,       kEBkg_Phi_wrapped,        kThetaBkg_Phi_wrapped,        kThetaBkg_Phi, 
                     kEBkg_pi0_Theta,    kEBkg_pi0_Phi,   kEBkg_pi0_Phi_wrapped,    kThetaBkg_pi0_Phi_wrapped,    kThetaBkg_pi0_Phi, 
                     kEBkg_cosmic_Theta, kEBkg_cosmic_Phi,kEBkg_cosmic_Phi_wrapped, kThetaBkg_cosmic_Phi_wrapped, kThetaBkg_cosmic_Phi, 
                     kEBkg_other_Theta,  kEBkg_other_Phi, kEBkg_other_Phi_wrapped,  kThetaBkg_other_Phi_wrapped,  kThetaBkg_other_Phi, 
                     kEBkg_Theta_unselected,        kEBkg_Phi_unselected,       kEBkg_Phi_wrapped_unselected,        kThetaBkg_Phi_wrapped_unselected,        kThetaBkg_Phi_unselected, 
                     kEBkg_pi0_Theta_unselected,    kEBkg_pi0_Phi_unselected,   kEBkg_pi0_Phi_wrapped_unselected,    kThetaBkg_pi0_Phi_wrapped_unselected,    kThetaBkg_pi0_Phi_unselected, 
                     kEBkg_cosmic_Theta_unselected, kEBkg_cosmic_Phi_unselected,kEBkg_cosmic_Phi_wrapped_unselected, kThetaBkg_cosmic_Phi_wrapped_unselected, kThetaBkg_cosmic_Phi_unselected, 
                     kEBkg_other_Theta_unselected,  kEBkg_other_Phi_unselected, kEBkg_other_Phi_wrapped_unselected,  kThetaBkg_other_Phi_wrapped_unselected,  kThetaBkg_other_Phi_unselected, 
                     kTH2D_names_MAX};

    std::vector<std::vector<TH1D*>> TH1D_hist_weighted;
    enum TH1D_names_weighted{ 
                     kldg_shwr_Phi_w,           kldg_shwr_Phi_wrapped_w,      kldg_shwr_Theta_w, 
                     kselected_w,
                     kshower_phi_pi0_w,         kshower_phi_bkg_cosmic_w,         kshower_phi_other_w,             kshower_E_w, 
                     kshower_phi_pi0_wrapped_w, kshower_phi_bkg_cosmic_wrapped_w, kshower_phi_other_wrapped_w,
                     kshower_E_pi0_w,           kshower_E_bkg_cosmic_w,           kshower_E_other_w,
                     kshower_Theta_pi0_w,       kshower_Theta_bkg_cosmic_w,       kshower_Theta_other_w,
                     kTH1D_names_weighted_MAX};

    std::vector<std::string> bnbvars = {
                                        "BNBwithDIC",
                                        "BNBdataSCE",
                                        "BNBdeadSaturatedChannels",
                                        "BNBLArG4BugFix",
                                        "BNBDLup",
                                        "BNBDLdown",
                                        "BNBDTup",
                                        "BNBDTdown",
                                        "BNBnoiseAmpUp",
                                        "BNBnoiseAmpDown",
                                        "BNBaltDeadChannels",
                                        "BNBstretchResp",
                                        "BNBsqueezeResp",
                                        "BNBupPEnoise",
                                        "BNBEnhancedTPCVis",
                                        "BNBBirksRecomb",
                                        "BNBdownPEnoise"};

    // ----------------------
    //      Other
    // ----------------------
    TFile * f_var_out;
    TFile *fweight;

    double largest_flash_y{0};
    double largest_flash_z{0};
    double largest_flash_time{0};
    double largest_flash_pe{0};
    
    TTree* VariableTree = new TTree("VariableTree","VariableTree");
    double mc_Phi_pi0 = -999, mc_Phi_other = -999, mc_Phi = -999, mc_Theta = -999;
    double mc_Energy = -999;
    std::string bkg_class;

    std::string draw_mode = "";

    TTree *selected_tree = new TTree("selected_tree","selected_tree");

    private:

};

#endif
