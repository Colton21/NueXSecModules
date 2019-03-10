#ifndef VARIATION_OUTPUT_H
#define VARIATION_OUTPUT_H

// ------------------------------------------------
// Class description
// A class to get cut variables from the selection
// and add them to a root file for further inspection.
// The intended use of this is to investigate
// the effect of the various detector systematics.

// USAGE
// mode = "same" || ""
// ./main.exe --var_mode <MC_Var_File> <mode> 
// ------------------------------------------------

// Library and Other Includes
#include "../xsecAna/TpcObjectContainer.h"
#include "../xsecAna/ParticleContainer.h"
// #include "histogram_functions.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

// C++ includes
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <math.h>

#include "../xsecAna/LinkDef.h"

class variation_output {
	public:
	variation_output()=default; //def_con
	
	// ----------------------
	//   Main Function
	// ----------------------
	void run_var(const char * file1, TString mode);

	// ----------------------
	//   Other Functions
	// ----------------------
	int GetLeadingShowerIndex(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj);     // Returns the index of the leading shower
	double GetLongestTrackLength(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj);  // Returns the length of the longest track
	void DrawTH1D(TH1D* h, double POT_Scaling); // Function that draws a TH1D histogram
	double GetPOT(const char * _file1); 		// Gets the POT stored in an external file
	void PlotVariatons(TFile* f_var_out); 		// Plots the variation files on the same plot
	std::vector<std::string> GrabDirs(TFile* f_var_out); // Grabs the directories in the file
	void DrawTH1D_SAME(TH1D* hist, std::string variation, TLegend* legend, std::string histname); 	// Function that draws a TH1D histogram for the same plot
	void GetNumber_Track_Shower(const int n_pfp, int n_tpc_obj,
									 xsecAna::TPCObjectContainer tpc_obj, int &n_showers, int &n_tracks,
									 int &n_pfp_50Hits, int &n_tracks_50Hits, int &n_showers_50Hits); // Utility function to get the number of tracks and showers
	double pfp_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z); // Calculates the pfp to nu vertex distance
	// Flash Functions
	std::vector<std::vector<double>> GetLargestFlashVector(TFile* f); 				// Function to resize opical entries to same size of events and get largest flash vector
	bool flash_in_time(double flash_time, double flash_start, double flash_end); 	// Decides whether flash is in time or not
	bool flash_pe(int flash_pe, int flash_pe_threshold); 							// Decides whether flash has sufficient PE
	double Flash_TPCObj_vtx_Dist(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z); // Returns the 2D distance of the flash to TPC OBj Vertex
	// ----------------------
	//      Variables
	// ----------------------
	double ldg_shwr_hits{0}; 			// Leading shower hits (all planes)
	double ldg_shwr_hits_WPlane{0}; 	// Leading Shower Collection plane hits
	double lldg_shwr_Open_Angle{0}; 	// Leading Shower Opening angle of shower
	double ldg_shwr_dEdx_WPlane{0};		// Leading shower dE/dx at Collection Plane
	double ldg_shwr_HitPerLen{0};		// Leading shower hits per length
	double ldg_shwr_Phi{0};				// Leading shower phi
	double ldg_shwr_Theta{0};			// Leading shower theta
	double ldg_shwr_CTheta{0};			// Leading shower cos theta
	double long_Track_ldg_shwr{0};		// Longest track / leading shower length
	// Leading shower momentum
	int n_pfp{0};					// Number of pfp in TPC Obj
	int n_pfp_50Hits{0};			// Number of pfp in TPC Obj with > 50 Hits
	int n_tracks{0};				// Number of Tracks in TPC Obj
	int n_tracks_50Hits{0};			// Number of Tracks in TPC Obj with > 50 Hits
	int n_showers{0};				// Number of Showers in TPC Obj
	int n_showers_50Hits{0};		// Number of showers in TPC Obj with > 50 Hits
	double track_phi{0};			// Track Phi
	
	
	double tpc_obj_vtx_x{0}, tpc_obj_vtx_y{0}, tpc_obj_vtx_z{0}; // TPCObj Vertex X, Y, Z
	// Distance of nue vertex and shower vertex
	// Distance of nue vertex and track vertex

	// ----------------------
	//      Histograms
	// ----------------------
	// Leading shower
	TH1D* h_ldg_shwr_hits 			= new TH1D("h_ldg_shwr_hits","ldg_shwr_hits", 				 50, 0, 600); 
	TH1D* h_ldg_shwr_hits_WPlane 	= new TH1D("h_ldg_shwr_hits_WPlane","ldg_shwr_hits_WPlane",  30, 0, 1000);
	TH1D* h_ldg_shwr_Open_Angle		= new TH1D("h_ldg_shwr_Open_Angle","ldg_shwr_Open_Angle", 	 25, 0, 50);
	TH1D* h_ldg_shwr_dEdx_WPlane 	= new TH1D("h_ldg_shwr_dEdx_WPlane","ldg_shwr_dEdx_WPlane",  40, 0, 10);
	TH1D* h_ldg_shwr_HitPerLen 		= new TH1D("h_ldg_shwr_HitPerLen","ldg_shwr_HitPerLen", 	 20, 0, 20);
	TH1D* h_ldg_shwr_Phi 			= new TH1D("h_ldg_shwr_Phi","ldg_shwr_Phi", 				 12, -180, 180);
	TH1D* h_ldg_shwr_Theta 			= new TH1D("h_ldg_shwr_Theta","ldg_shwr_Theta", 			 12, 0, 180);
	TH1D* h_ldg_shwr_CTheta 		= new TH1D("h_ldg_shwr_CTheta","ldg_shwr_CTheta", 			 16, -1, 1);
	TH1D* h_long_Track_ldg_shwr 	= new TH1D("h_long_Track_ldg_shwr","long_Track_ldg_shwr",	 20, 0, 3);
	
	// TPC
	TH1D* h_tpc_obj_vtx_x 			= new TH1D("h_tpc_obj_vtx_x","tpc_obj_vtx_x", 				 20, 0, 260);
	TH1D* h_tpc_obj_vtx_y 			= new TH1D("h_tpc_obj_vtx_y","tpc_obj_vtx_y", 				 20, -117, 117);
	TH1D* h_tpc_obj_vtx_z 			= new TH1D("h_tpc_obj_vtx_z","tpc_obj_vtx_z", 				 40, 0, 1040);
	
	// Other
	TH1D* h_total_hits 				= new TH1D("h_total_hits", "h_total_hits",	 				 50, 0, 600);
	TH1D* h_n_pfp 					= new TH1D("h_n_pfp", "h_n_pfp", 							 8, 0, 8);
	TH1D* h_n_pfp_50Hits 			= new TH1D("h_n_pfp_50Hits", "h_n_pfp_50Hits", 				 8, 0, 8);
	TH1D* h_n_tracks 				= new TH1D("h_n_tracks", "h_n_tracks", 						 8, 0, 8);
	TH1D* h_n_tracks_50Hits 		= new TH1D("h_n_tracks_50Hits", "h_n_tracks_50Hits", 		 8, 0, 8);
	TH1D* h_n_showers 				= new TH1D("h_n_showers", "h_n_showers", 					 8, 0, 8);
	TH1D* h_n_showers_50Hits 		= new TH1D("h_n_showers_50Hits", "h_n_showers_50Hits", 		 8, 0, 8);
	TH1D* h_track_phi 				= new TH1D("h_track_phi", "h_track_phi",					 12 , -180 ,180);
	TH1D* h_shower_phi 				= new TH1D("h_shower_phi", "h_shower_phi",					 12 , -180 ,180); // Shower Phi

	TH1D* h_shower_Nu_vtx_Dist		= new TH1D("h_shower_Nu_vtx_Dist","h_shower_Nu_vtx_Dist",	 20, 0, 20);
	TH1D* h_track_Nu_vtx_Dist		= new TH1D("h_track_Nu_vtx_Dist","h_track_Nu_vtx_Dist",	 	 20, 0, 20);

	// Largest Flash
	TH1D* h_largest_flash_y			= new TH1D("h_largest_flash_y", "h_largest_flash_y", 		60, -40, 40);
	TH1D* h_largest_flash_z			= new TH1D("h_largest_flash_z", "h_largest_flash_z", 		125, 0, 1000);
	TH1D* h_largest_flash_time		= new TH1D("h_largest_flash_time", "h_largest_flash_time",	50, 0, 20);
	TH1D* h_largest_flash_pe		= new TH1D("h_largest_flash_pe", "h_largest_flash_pe", 		30, 0, 6000);
	TH1D* h_Flash_TPCObj_Dist		= new TH1D("h_Flash_TPCObj_Dist", "h_Flash_TPCObj_Dist", 	50, 0, 200); // Largest flash to TPC OBj Vtx Dist
	// ----------------------
	//      Other
	// ----------------------
	TFile * f_var_out = new TFile("plots/variation_out.root","UPDATE");
	bool PlotVar{false};

	double largest_flash_y{0};
	double largest_flash_z{0};
	double largest_flash_time{0};
	double largest_flash_pe{0};



	private:

};

#endif
