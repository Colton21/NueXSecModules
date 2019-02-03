#include "variation_output.h"

//***************************************************************************
//***************************************************************************
int variation_output::GetLeadingShowerIndex(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj){
    int leading_index{0};
    int most_hits{0};
	
    // Loop over Particle Objects
    for (int j = 0; j < n_pfp; j++) {
        auto const pfp_obj = tpc_obj.GetParticle(j);
        const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
        
        if (pfp_pdg != 11) continue; // skip if not a shower

        const int n_pfp_hits = pfp_obj.NumPFPHits();
        
        // Compare for most hits
        if (n_pfp_hits > most_hits) {
            leading_index = j; 
            most_hits = n_pfp_hits; 
        }
    }

    return leading_index;
}
//***************************************************************************
//***************************************************************************
void variation_output::DrawTH1D(TH1D* h, double POT_Scaling){
	TCanvas* c = new TCanvas();
	c->cd();

	h->SetLineColor(kMagenta+3);
	h->SetLineWidth(2);
	h->SetLineStyle(1);

	h->Scale(POT_Scaling);
	h->Draw("");

}
//***************************************************************************
//***************************************************************************
void variation_output::DrawTH1D_SAME(TH1D* hist, std::string variation, TLegend* legend, std::string histname){

	
	// ----------------------
	//    Axis Specifiers
	// ----------------------
	if (histname == "h_total_hits"){
		hist->SetTitle("Total Hits (All Planes);PFP Total Hits;Entries");
		hist->GetYaxis()->SetRangeUser(0,4500);

	}
	else if (histname == "h_ldg_shwr_hits") {
		hist->SetTitle("Leading Shower Hits (All Planes); Leading Shower Hits;Entries");
		hist->GetYaxis()->SetRangeUser(0,750);
	}
	else if (histname == "h_ldg_shwr_hits_WPlane") {
		hist->SetTitle("Leading Shower Hits (Collection Plane);Leading Shower Collection Plane Hits;Entries");
		hist->GetYaxis()->SetRangeUser(0,5000);
	}
	else if (histname == "h_ldg_shwr_Open_Angle"){
		hist->SetTitle("Leading Shower Open Angle;Leading Shower Opening Angle [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,175);
	}
	else if (histname == "h_ldg_shwr_dEdx_WPlane"){
		hist->SetTitle("Leading Shower dEdx Collection Plane;Leading Shower dEdx Collection Plane [MeV/cm];Entries");
		hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_HitPerLen"){
		hist->SetTitle("Leading Shower Hits / Length;Leading Shower Hits / Length [cm^-1];Entries");
		hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_Phi"){
		hist->SetTitle("Leading Shower #phi;Leading Shower #phi [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_Theta"){
		hist->SetTitle("Leading Shower #theta;Leading Shower #theta [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_CTheta"){
		hist->SetTitle("Leading Shower cos(#theta);Leading Shower cos(#theta);Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	
	else return;

	// ----------------------
	//    Draw Specifiers
	// ----------------------
	if (variation == "CV"){
	hist->SetLineColor(kBlack);
	hist->SetLineWidth(2);
	hist->SetLineStyle(1);
	legend->AddEntry(hist, "CV", "l");
	hist->Draw("hist,same");
	} 
	else if  (variation == "withDIC"){
		hist->SetLineColor(kMagenta+2);
		hist->SetLineWidth(2);
		hist->SetLineStyle(1);
		legend->AddEntry(hist, "withDIC", "l");
		hist->Draw("hist,same");
	}
	else if  (variation == "EnhancedTPCVis" ){ 
		hist->SetLineColor(30);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "EnhancedTPCVis", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "altDeadChannels"){ 
		hist->SetLineColor(38);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "altDeadChannels", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "deadSaturatedChannels"){
		hist->SetLineColor(28);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "deadSaturatedChannels", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
		
	}
	else if  (variation == "stretchResp" ){
		hist->SetLineColor(36);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "stretchResp", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "squeezeResp"){
		hist->SetLineColor(1001);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "squeezeResp", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	// else if  (variation == "PPFXThinMeson"){
	// 	hist->SetLineColor(kBlue+1);
	// 	hist->SetLineWidth(2);
	// 	legend->AddEntry(hist, "Meson Incident.", "l");
	// 	hist->SetLineStyle(2);
	// 	hist->Draw("hist,same");
	// }
	// else if  (variation == "PPFXThinNeutron"){
	// 	hist->SetLineColor(42);
	// 	hist->SetLineWidth(2);
	// 	legend->AddEntry(hist, "nC #rightarrow #piX", "l");
	// 	hist->SetLineStyle(2);
	// 	hist->Draw("hist,same");
	// }
	else return;
}
//***************************************************************************
//***************************************************************************
double variation_output::GetPOT(const char * _file1){
	double POT{0};
	std::string line;

	std::string filename;
	std::string temp_filename = _file1; // cast to string
	filename = "File: " + temp_filename;

	std::cout << filename << std::endl;

	std::ifstream myfile ("POT_List_Varaitions.txt");
	int i_POT{0};
	
	if (myfile.is_open()) {

		// Loop over lines in file
		while ( getline (myfile,line) ) {

			if (i_POT == 1){
				line.erase(0, 21);
				POT = std::stod(line); // Convert string to double
				std::cout << "POT in File:\t"<< POT << std::endl;
				return POT;
			}

			// Found the correct variation file 
			if (line == filename){
				std::cout << "Found match for POT file"<< std::endl;
				i_POT++;
			}
			
		}

		myfile.close();
	}

	else std::cout << "Unable to open file" << std::endl; 
	return POT;
}
//***************************************************************************
//***************************************************************************
// Function that grabs the reweighted histogram names for plotting
std::vector<std::string> variation_output::GrabDirs(TFile* f_var_out) {
	std::vector<std::string> variations;

	f_var_out->cd();
	
	TKey *key;
	TIter nextkey(gDirectory->GetListOfKeys());

	TString same_plots = "SAME_Plots";

	std::cout << "\n=================================================" << std::endl;	
	std::cout << "Getting variation modes:" << std::endl;	
  	while ( ( key =  (TKey*)nextkey()) ) { // Extra brackets to omit a warning 
    	if (key->IsFolder()) {
			std::cout << key->GetName() << std::endl; // Print the variations
			if (key->GetName() == same_plots ) continue; // Skip this
			variations.push_back(key->GetName());
		}
	}
	std::cout << "=================================================\n" << std::endl;

	return (variations);
}
//***************************************************************************
//***************************************************************************
void variation_output::PlotVariatons(TFile* f_var_out){
	f_var_out->cd();
	
	// Grab the variation folders in the file
	std::vector<std::string> variations = variation_output::GrabDirs(f_var_out); 

	TH1D* hist;
	
	std::vector<std::string> histnames = {"h_total_hits","h_ldg_shwr_hits", "h_ldg_shwr_hits_WPlane",
										 "h_ldg_shwr_Open_Angle", "h_ldg_shwr_dEdx_WPlane", "h_ldg_shwr_HitPerLen",
										 "h_ldg_shwr_Phi", "h_ldg_shwr_Theta","h_ldg_shwr_CTheta"};

	// Loop over the histograms
	for (int j=0; j < histnames.size(); j++){
		
		// Canvas + Legend
		TCanvas* c = new TCanvas();
		// TLegend* legend = new TLegend(0.5, 0.65, 0.9, 0.9);
		TLegend* legend =new TLegend(0.75, 0.60, 0.95, 0.95);
		// legend->SetNColumns(2);
		// legend->SetBorderSize(0);
		// legend->SetFillStyle(0);

		// Loop over variation directories
		for (int i=0; i < variations.size(); i++){
			char name[500];
			snprintf(name, 500, "%s/%s", variations[i].c_str(),histnames[j].c_str() );

			hist = (TH1D*)f_var_out->Get(name);

			if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!" << std::endl;

			DrawTH1D_SAME(hist, variations[i], legend, histnames[j]);

		}

		// Print the Canvas
		char Canvas_name[500];
		snprintf(Canvas_name, 500, "plots/%s.pdf",histnames[j].c_str() ); 
		legend->Draw();
		c->Print(Canvas_name);

	}

	// Close the file
	f_var_out->Close();

}
//***************************************************************************
//***************************************************************************
// ----------------------
//         Main
// ----------------------
void variation_output::run_var(const char * _file1, TString mode) {
	// std::cout << "=================================================\n" << std::endl;
	// std::cout << "Warning, there are hardcoded values" << 
	// 	"in this script, grep for  \"HARDCODED\" for places\n" << std::endl;
	// std::cout << "=================================================\n" << std::endl;
	
	// create plots folder if it does not exist
	system("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi");

	gStyle->SetOptStat(0); // say no to stats box

	//***************************************************************************
	// if bool true just run this function
	if (mode == "same") PlotVar = true; 
	else PlotVar = false; 

	if (PlotVar) {
		PlotVariatons(f_var_out); 
		return; 
	}
	//***************************************************************************

	// Do some POT scaling
	std::cout << "=================================================\n" << std::endl;
	double CV_POT =  GetPOT("files/filter_CV.root");
	double POT_Scaling =  CV_POT / GetPOT(_file1);
	std::cout << "POT Scaling:\t" << POT_Scaling << std::endl;
	std::cout << "=================================================\n" << std::endl;
	
	// Check if the outfile opened successfully
	if ( f_var_out->IsOpen() ) std::cout << "Variation File opened successfully\n" << std::endl;
    
    // Get Variation file TPCObj Tree 
    TFile* inFile = new TFile(_file1);
    TTree* TPCObjTree = (TTree*) inFile->Get("AnalyzeTPCO/tree");
    
    std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
    TPCObjTree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

    // Num events in the tree
    const int tree_total_entries = TPCObjTree->GetEntries();
    std::cout << "Total Events: " << tree_total_entries << std::endl;

    // ----------------------
    //      Event loop
    // ----------------------
    std::cout << "Starting Eventloop..." << std::endl;
    for(int event = 0; event < tree_total_entries; event++){
        
        TPCObjTree->GetEntry(event);

        int n_tpc_obj = tpc_object_container_v->size();
        
        // Loop over TPCObj
        for(int i = 0; i < n_tpc_obj; i++){
            auto const tpc_obj = tpc_object_container_v->at(i);

            // TPC Obj vars
            tpc_obj_vtx_y = tpc_obj.pfpVtxY();
            tpc_obj_vtx_z = tpc_obj.pfpVtxZ();
            tpc_obj_vtx_x = tpc_obj.pfpVtxX();
           
            const int tpc_obj_mode = tpc_obj.Mode(); 
            const int n_pfp = tpc_obj.NumPFParticles();
            const int leading_shower_index = GetLeadingShowerIndex(n_pfp, n_tpc_obj, tpc_obj);

            // Loop over the Par Objects
            for (int j = 0; j < n_pfp ; j++){

                auto const pfp_obj = tpc_obj.GetParticle(j);

                // PFP vars                
                const std::string mc_origin = pfp_obj.Origin();
                const int  pfp_pdg          = pfp_obj.PFParticlePdgCode();
                const int  num_pfp_hits     = pfp_obj.NumPFPHits();
                const int  mc_parent_pdg    = pfp_obj.MCParentPdg();
                
                const double pff_length     = pfp_obj.pfpLength();
                const double pfp_open_angle = pfp_obj.pfpOpenAngle();

                // CC && BeamNu 
                if( pfp_obj.CCNC() == 0 && mc_origin == "kBeamNeutrino") {
                       
                    // Electron (Shower like) && (Nue || Nuebar)
                    if ( pfp_pdg == 11 && (mc_parent_pdg == 12 || mc_parent_pdg == -12) ) {
                        h_total_hits->Fill(num_pfp_hits);

						// Leading shower
						if (j == leading_shower_index){
							const double leading_shower_phi = atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
							const double leading_shower_theta = acos(pfp_obj.pfpDirZ()) * 180 / 3.1415;

							h_ldg_shwr_hits->Fill(num_pfp_hits);
							h_ldg_shwr_hits_WPlane->Fill(pfp_obj.NumPFPHitsW()); // W Plane
							h_ldg_shwr_Open_Angle->Fill(pfp_obj.pfpOpenAngle() * (180 / 3.1415) );
							h_ldg_shwr_dEdx_WPlane->Fill( pfp_obj.PfpdEdx().at(2) ); // W Plane
							h_ldg_shwr_HitPerLen->Fill( num_pfp_hits / pfp_obj.pfpLength() );
							h_ldg_shwr_Phi->Fill(leading_shower_phi);
							h_ldg_shwr_Theta->Fill(leading_shower_theta);
							h_ldg_shwr_CTheta->Fill( cos(leading_shower_theta) );
						}

                    }
                    
                }     

            } // END LOOP PAR OBJ

        } // END LOOP TPCO

    } // END EVENT LOOP
    std::cout << "Finished Eventloop..." << std::endl;

    
	// ----------------------
    //    Save to a file
    // ----------------------
	TDirectory* savedir = gDirectory; //  Create the directory
	TDirectory* subdir;

	// Get the directory for the file	
	std::string dirname = _file1;
	
	// Get the variation name by stripping the input file name -- HARDCODED
	// Required format "files/filter_<NAME>.root"
	dirname.erase(0,13); // - "files/filter_"
	dirname.erase(dirname.end()-5, dirname.end()); // - ".root"

	std::cout << "dirname:\t" << dirname << std::endl;

	f_var_out->cd();

	// If directory does not exist then make it
	savedir = (TDirectory*)f_var_out->Get(dirname.c_str());
	if (savedir == NULL ) {
		savedir = gDirectory;
		std::cout << dirname << " directory does not exist, creating..." << std::endl;
		subdir = savedir->mkdir(dirname.c_str() ) ;
	}
	else {
		std::cout << dirname <<" directory exists, overwriting..." << std::endl;
		subdir = savedir;
	}
	
	subdir->cd();

	// ----------------------
    //    Draw Histograms
    // ----------------------
    DrawTH1D(h_total_hits, POT_Scaling);
	DrawTH1D(h_ldg_shwr_hits, POT_Scaling);
	DrawTH1D(h_ldg_shwr_hits_WPlane, POT_Scaling);
	DrawTH1D(h_ldg_shwr_Open_Angle, POT_Scaling);
	DrawTH1D(h_ldg_shwr_dEdx_WPlane, POT_Scaling);
	DrawTH1D(h_ldg_shwr_HitPerLen, POT_Scaling);
	DrawTH1D(h_ldg_shwr_Phi, POT_Scaling);
	DrawTH1D(h_ldg_shwr_Theta, POT_Scaling);
	DrawTH1D(h_ldg_shwr_CTheta, POT_Scaling);
	// ----------------------
    //   Write and close
	//   the TFile to new/updated
	//   directory
    // ----------------------
	h_total_hits->Write("",TObject::kOverwrite);
	h_ldg_shwr_hits->Write("",TObject::kOverwrite);
	h_ldg_shwr_hits_WPlane-> Write("",TObject::kOverwrite);
	h_ldg_shwr_Open_Angle->Write("", TObject::kOverwrite);
	h_ldg_shwr_dEdx_WPlane->Write("", TObject::kOverwrite);
	h_ldg_shwr_HitPerLen->Write("", TObject::kOverwrite);
	h_ldg_shwr_Phi->Write("", TObject::kOverwrite);
	h_ldg_shwr_Theta->Write("", TObject::kOverwrite);
	h_ldg_shwr_CTheta->Write("", TObject::kOverwrite);
	
	f_var_out->Close(); 

} // END MAIN
