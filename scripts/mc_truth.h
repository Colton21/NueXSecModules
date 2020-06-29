#ifndef MC_TRUTH_H
#define MC_TRUTH_H

// ------------------------------------------------
// Class description
// A class to get mc truth variables from the selecton
// so we can see the effect of event weights on them

// ./main.exe --mc_truth <MC_File>
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
#include "TVector3.h"

// C++ includes
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <math.h>

#include "../xsecAna/LinkDef.h"



class mc_truth {
    public:

    // ----------------------
    //   Main Function
    // ----------------------
    void run_var(const char * file1, const std::vector<double> _config);

    // vector containing all the signal type events
    std::vector<std::string> signal_modes = {"nue_cc_qe",     "nue_cc_res",     "nue_cc_dis",     "nue_cc_coh",     "nue_cc_mec",
                                             "nue_bar_cc_qe", "nue_bar_cc_res", "nue_bar_cc_dis", "nue_bar_cc_coh", "nue_bar_cc_mec" };

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

    void merge_weights_into_tree();

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
    
    
    double tpc_obj_vtx_x{0}, tpc_obj_vtx_y{0}, tpc_obj_vtx_z{0}; // TPCObj Vertex X, Y, Z

    std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
    std::vector<xsecAna::TPCObjectContainer> tpc_object_container_v_copy;

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
    int counter_ext{0};


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
    double nue_cc_counter{0.0};

    double largest_flash_y{0};
    double largest_flash_z{0};
    double largest_flash_time{0};
    double largest_flash_pe{0};

    int run{0}, subrun{0}, evt{0};

    std::string type = "mc"; // Use this to switch between MC and EXT

    TH2D* h_shr_dEdx_shr_dist_sig = new TH2D( "h_reco_shr_dEdx_shr_dist_sig",";Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
    TH2D* h_shr_dEdx_shr_dist_bkg = new TH2D( "h_reco_shr_dEdx_shr_dist_bkg",";Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
    TH2D* h_shr_dEdx_shr_dist_ext = new TH2D( "h_reco_shr_dEdx_shr_dist_ext",";Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);

    TH2D* h_shr_dEdx_shr_dist_sig_good = new TH2D( "h_reco_shr_dEdx_shr_dist_sig_good","0 < #theta < 60;Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 20, 0, 6, 10, 0, 20);
    TH2D* h_shr_dEdx_shr_dist_bkg_good = new TH2D( "h_reco_shr_dEdx_shr_dist_bkg_good","0 < #theta < 60;Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 20, 0, 6, 10, 0, 20);
    TH2D* h_shr_dEdx_shr_dist_ext_good = new TH2D( "h_reco_shr_dEdx_shr_dist_ext_good","0 < #theta < 60;Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 20, 0, 6, 10, 0, 20);

    
    
    bool fill_tree{false}; // Decide if to fill the tree
    bool write_tree{true};


    private:

};

#endif
