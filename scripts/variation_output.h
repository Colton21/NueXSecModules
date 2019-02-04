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
	TH1D* h_tpc_obj_vtx_x 			= new TH1D("h_tpc_obj_vtx_x","tpc_obj_vtx_x", 				 21, 0, 260);
	TH1D* h_tpc_obj_vtx_y 			= new TH1D("h_tpc_obj_vtx_y","tpc_obj_vtx_y", 				 19, 0, 115);
	TH1D* h_tpc_obj_vtx_z 			= new TH1D("h_tpc_obj_vtx_z","tpc_obj_vtx_z", 				 39, 0, 1040);
	
	// Other
	TH1D* h_total_hits 				= new TH1D("h_total_hits", "h_total_hits",	 				 50, 0, 600);

    // ----------------------
    //      Other
    // ----------------------
	TFile * f_var_out = new TFile("plots/variation_out.root","UPDATE");
	bool PlotVar{false};

    private:

};

#endif
