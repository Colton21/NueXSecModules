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
double variation_output::GetLongestTrackLength(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj){
	double longest_track = 0;

	// Loop over Particle Objects
	for(int i = 0; i < n_pfp; i++) {

		auto const pfp_obj = tpc_obj.GetParticle(i);
		const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
		
		// If it  is track like
		if(pfp_pdg == 13) {
			const double trk_length = pfp_obj.pfpLength();
			
			if(trk_length > longest_track) {

				longest_track = trk_length;
			}
		}
	}

	return longest_track;
}
//***************************************************************************
//***************************************************************************
void variation_output::GetNumber_Track_Shower(const int n_pfp, int n_tpc_obj,
									 xsecAna::TPCObjectContainer tpc_obj, int &n_showers, int &n_tracks,
									 int &n_pfp_50Hits, int &n_tracks_50Hits,int &n_showers_50Hits) { 
	// Loop over Particle Objects
	for (int j = 0; j < n_pfp; j++) {
		auto const pfp_obj = tpc_obj.GetParticle(j);
		const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
		
		// Don't count the neutrinos
		if (pfp_pdg != 12 || pfp_pdg != -12 || pfp_pdg != 14 || pfp_pdg != -14){
			// Cut with pfp with > 50 Hits
			if (pfp_obj.NumPFPHits() > 50){
				// For > 50 hits
				n_pfp_50Hits++; // Add one to n_pfp > 50 hits counter
				
				if (pfp_pdg == 11) 		n_showers_50Hits++; // Add to shower counter
				else if (pfp_pdg == 13) n_tracks_50Hits++;  // Add to track counter
				else std::cout << "Unknown pandora classification:\t" << pfp_pdg << std::endl;
					
			}

			// All PFP
			if (pfp_pdg == 11) 		n_showers++; // Add to shower counter
			else if (pfp_pdg == 13) n_tracks++;  // Add to track counter
			else return;
		}
		
	}
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
		hist->SetTitle(";PFP Total Hits;Entries");
		hist->GetYaxis()->SetRangeUser(0,4500);
	}
	else if (histname == "h_ldg_shwr_hits") {
		hist->SetTitle("; Leading Shower Hits (All Planes);Entries");
		hist->GetYaxis()->SetRangeUser(0,750);
	}
	else if (histname == "h_ldg_shwr_hits_WPlane") {
		hist->SetTitle(";Leading Shower Collection Plane Hits;Entries");
		hist->GetYaxis()->SetRangeUser(0,5000);
	}
	else if (histname == "h_ldg_shwr_Open_Angle"){
		hist->SetTitle(";Leading Shower Opening Angle [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,175);
	}
	else if (histname == "h_ldg_shwr_dEdx_WPlane"){
		hist->SetTitle(";Leading Shower dEdx Collection Plane [MeV/cm];Entries");
		hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_HitPerLen"){
		hist->SetTitle(";Leading Shower Hits / Length [ cm^{-1} ];Entries");
		hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_Phi"){
		hist->SetTitle(";Leading Shower #phi [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_Theta"){
		hist->SetTitle(";Leading Shower #theta [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_CTheta"){
		hist->SetTitle(";Leading Shower cos(#theta);Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_long_Track_ldg_shwr"){
		hist->SetTitle(";Longest Track Length / Leading Shower Length;Entries");
		hist->GetYaxis()->SetRangeUser(0,22000);
	}
	else if (histname == "h_tpc_obj_vtx_x"){
		hist->SetTitle(";TPC Object Vertex X [cm];Entries");
		hist->GetYaxis()->SetRangeUser(0,2500);
	}
	else if (histname == "h_tpc_obj_vtx_y"){
		hist->SetTitle(";TPC Object Vertex Y [cm];Entries");
		hist->GetYaxis()->SetRangeUser(0,2200);
	}
	else if (histname == "h_tpc_obj_vtx_z"){
		hist->SetTitle(";TPC Object Vertex Z [cm];Entries");
		hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_pfp"){
		hist->SetTitle(";Number of PFP;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_pfp_50Hits"){
		hist->SetTitle("; Number of PFP > 50 Hits;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_tracks"){
		hist->SetTitle(";Number of Tracks;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_tracks_50Hits"){
		hist->SetTitle("; Number of Tracks > 50 Hits;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_showers"){
		hist->SetTitle(";Number of Showers;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_showers_50Hits"){
		hist->SetTitle("; Number of Showers > 50 Hits;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_track_phi"){
		hist->SetTitle("; Track #phi [degrees];Entries");
		hist->GetYaxis()->SetRangeUser(0,6000);
	}
	else if (histname == "h_shower_phi"){
		hist->SetTitle("; Shower #phi [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}			
	else if (histname == "h_largest_flash_y"){
		hist->SetTitle(";Largest Flash Y [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_largest_flash_z"){
		hist->SetTitle("; Largest Flash Z [cm];Entries");
		hist->GetYaxis()->SetRangeUser(0,1000);
	}
	else if (histname == "h_largest_flash_time"){
		hist->SetTitle("; Largest Flash Time [#mus];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_largest_flash_pe"){
		hist->SetTitle("; Largest Flash PE;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_Flash_TPCObj_Dist"){
		hist->SetTitle("; 2D Distance from Largest Flash to Reco #nu Vertex [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
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
		legend->AddEntry(hist, "DIC", "l");
		hist->Draw("hist,same");
	}
	else if  (variation == "EnhancedTPCVis" ){ 
		hist->SetLineColor(30);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Enhanced TPC Vis.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "altDeadChannels"){ 
		hist->SetLineColor(38);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Alt. Dead Chan.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "deadSaturatedChannels"){
		hist->SetLineColor(28);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Dead Sat. Chan.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
		
	}
	else if  (variation == "stretchResp" ){
		hist->SetLineColor(36);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Stretch Resp.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "squeezeResp"){
		hist->SetLineColor(1001);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Squeeze Resp.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "upPEnoise"){
		hist->SetLineColor(kBlue+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "PE Noise Up", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "noiseAmpDown"){
		hist->SetLineColor(42);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Noise Amp. Down", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "downPEnoise"){
		hist->SetLineColor(50);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "PE Noise Down", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "noiseAmpUp"){
		hist->SetLineColor(kOrange+10);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Noise Amp. Up", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "DTdown"){
		hist->SetLineColor(kOrange+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "DT Down", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
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
	std::cout << "Could not find a match for POT file"<< std::endl;
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
bool variation_output::flash_in_time(double flash_time, double flash_start, double flash_end) {
	if(flash_time >= flash_start && flash_time <= flash_end) return true; // Pass in time
	else return false; // Fail in time
}
//***************************************************************************
//***************************************************************************
bool variation_output::flash_pe(int flash_pe, int flash_pe_threshold) {
	if (flash_pe >= flash_pe_threshold) return true; // Pass PE Thresh
	else return false; // Fail PE Thresh
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
										 "h_ldg_shwr_Phi", "h_ldg_shwr_Theta","h_ldg_shwr_CTheta",
										  "h_long_Track_ldg_shwr", "h_tpc_obj_vtx_x", "h_tpc_obj_vtx_y", "h_tpc_obj_vtx_z",
										  "h_n_pfp", "h_n_pfp_50Hits", "h_n_tracks", "h_n_tracks_50Hits", "h_n_showers",
										  "h_n_showers_50Hits", "h_track_phi", "h_shower_phi", "h_largest_flash_y", "h_largest_flash_z",
										  "h_largest_flash_time", "h_largest_flash_pe", "h_Flash_TPCObj_Dist" };

	// Loop over the histograms
	for (int j=0; j < histnames.size(); j++){
		
		// Canvas + Legend
		TCanvas* c = new TCanvas();
		TLegend* legend;
		if (histnames[j] == "h_ldg_shwr_CTheta") legend = new TLegend(0.15, 0.55, 0.37, 0.85); // Reposition
		else if (histnames[j] == "h_tpc_obj_vtx_x" || histnames[j] == "h_tpc_obj_vtx_y" || histnames[j] == "h_tpc_obj_vtx_z" )
				legend = new TLegend(0.15, 0.15, 0.37, 0.45); // Reposition
		else if (histnames[j] == "h_largest_flash_z" ) legend = new TLegend(0.35, 0.59, 0.57, 0.89);
		else legend = new TLegend(0.72, 0.59, 0.94, 0.89);

		legend->SetBorderSize(0);
		legend->SetFillStyle(0);

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
//************************** Flash Functions ********************************
std::vector<std::vector<double>> variation_output::GetLargestFlashVector(TFile* f ){

	// Get Optical Information from file
	TTree* optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

	std::vector<std::vector<double>> largest_flash_v_v;

	// ----------------------
	//    Optical Info
	// ----------------------
	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;
	double fOpFlashWidthY = 0;
	double fOpFlashWidthZ = 0;
	double fOpFlashCenterY = 0;
	double fOpFlashCenterZ = 0;

	optical_tree->SetBranchAddress("Run",              &fRun           ); // Run number of Flash
	optical_tree->SetBranchAddress("Event",            &fEvent         ); // Event number of Flash
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE     ); // PE of Flash
	optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime   ); // Time of flash
	optical_tree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY );
	optical_tree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ );
	optical_tree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY); // Flash Y center
	optical_tree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ); // Flash Z Center

	// Num events in the optical tree
	const int optical_entries = optical_tree->GetEntries();
	std::cout << "Total Optical Entries: " << optical_entries << std::endl;

	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;

	// Contains the entry number for a given OpFlash per event
	std::vector<int>					optical_list_pe;
	std::vector<std::vector<int> >		optical_list_pe_v;
	
	std::vector<double>					optical_list_flash_center_y; 
	std::vector<std::vector<double> >	optical_list_flash_center_y_v;
	
	std::vector<double>					optical_list_flash_center_z; 
	std::vector<std::vector<double> >	optical_list_flash_center_z_v;
	
	std::vector<double>					optical_list_flash_time;
	std::vector<std::vector<double> >	optical_list_flash_time_v;
	
	// Loop over the optical entries to get the largest flash vector
	
	// ----------------------
	// Resize the optical enties to be the same sizd as number of Events (TPC Obj)
	// ----------------------
	
	for(int i = 0; i < optical_entries; i++) {
		
		// Get the Optical entry
		optical_tree->GetEntry(i);

		current_run		= fRun;
		current_event 	= fEvent;

		// New event
		if(current_event != last_event) {
			optical_list_pe.clear();
			optical_list_flash_center_y.clear();
			optical_list_flash_center_z.clear();
			optical_list_flash_time.clear();

			optical_list_pe.push_back(fOpFlashPE);
			optical_list_flash_center_y.push_back(fOpFlashCenterY);
			optical_list_flash_center_z.push_back(fOpFlashCenterZ);
			optical_list_flash_time.push_back(fOpFlashTime);

		}
		// Same event
		if(current_event == last_event && current_run == last_run) {
			optical_list_pe_v.pop_back();
			optical_list_flash_center_y_v.pop_back();
			optical_list_flash_center_z_v.pop_back();
			optical_list_flash_time_v.pop_back();

			optical_list_pe.push_back(fOpFlashPE);
			optical_list_flash_center_y.push_back(fOpFlashCenterY);
			optical_list_flash_center_z.push_back(fOpFlashCenterZ);
			optical_list_flash_time.push_back(fOpFlashTime);

		}

		last_event = current_event;
		last_run   = current_run;

		optical_list_pe_v.push_back(optical_list_pe);
		optical_list_flash_center_y_v.push_back(optical_list_flash_center_y);
		optical_list_flash_center_z_v.push_back(optical_list_flash_center_z);
		optical_list_flash_time_v.push_back(optical_list_flash_time);

	}
	
	std::cout << "Resized Optical List Size: " << optical_list_pe_v.size() << std::endl;
	
	// Largest Flash Vector
	std::vector<double> largest_flash_v;
	largest_flash_v.resize(4, 0);
	
	// ----------------------
	//      Event loop
	// ----------------------
	for(int i = 0; i < optical_list_pe_v.size(); i++) {
		
		bool in_time 				= false;
		bool got_in_time 			= false;
		bool sufficient_flash 		= false;
		bool got_sufficient_flash 	= false;
		
		double largest_flash = 0.;
		double largest_center_y = 0;
		double largest_center_z = 0;
		double largest_flash_time = 0;
		
		// Cut Variables defined in main.h
		double flash_time_start = 5.5; 
		double flash_time_end = 16.0;
		int flash_pe_threshold = 50;

		// Loop through all flashes in event and find largest
		for(int j = 0; j < optical_list_pe_v.at(i).size(); j++) {
			
			auto const opt_time         = optical_list_flash_time_v.at(i).at(j);
			auto const opt_pe           = optical_list_pe_v.at(i).at(j);
			const double opt_center_y   = optical_list_flash_center_y_v.at(i).at(j);
			const double opt_center_z   = optical_list_flash_center_z_v.at(i).at(j);
			const double opt_flash_time = optical_list_flash_time_v.at(i).at(j);
			
			// See if flash was in time
			in_time = flash_in_time(opt_time, flash_time_start, flash_time_end); 
			if(in_time == true) got_in_time = true; 
			
			// See if flash meets the threshold requirements
			sufficient_flash = flash_pe(opt_pe, flash_pe_threshold); 
			if(sufficient_flash == true) got_sufficient_flash = true; 
			
			//Flash is both in time and over PE threshold
			if(in_time == true && sufficient_flash == true) {
				
				// Find the largest flash in this event
				if(opt_pe > largest_flash) {
					largest_flash      = opt_pe;
					largest_center_y   = opt_center_y;
					largest_center_z   = opt_center_z;
					largest_flash_time = opt_flash_time;
				}
			}
		}
		
		largest_flash_v.at(0) = largest_center_y;
		largest_flash_v.at(1) = largest_center_z;
		largest_flash_v.at(2) = largest_flash_time;
		largest_flash_v.at(3) = largest_flash;
		
		largest_flash_v_v.push_back(largest_flash_v);
		
	}

	return largest_flash_v_v;
}
//***************************************************************************
//***************************************************************************
double variation_output::Flash_TPCObj_vtx_Dist(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z) {
	const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
	return distance;
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

	//*************************** SAME PLOT *************************************
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
	
	// Get Variation file
	TFile* inFile = new TFile(_file1);
	
	// Get TPC Obj Information from file
	TTree* TPCObjTree = (TTree*) inFile->Get("AnalyzeTPCO/tree");
	
	std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
	TPCObjTree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

	// Num events in the tree
	const int tree_total_entries = TPCObjTree->GetEntries();
	std::cout << "Total Events: " << tree_total_entries << std::endl;

	// Get the largest Flash Vector of Vector	
	std::vector<std::vector<double>> largest_flash_v_v = GetLargestFlashVector(inFile);

	// ----------------------
	//      Event loop
	// ----------------------
	std::cout << "Starting Eventloop..." << std::endl;
	for(int event = 0; event < tree_total_entries; event++){
		
		TPCObjTree->GetEntry(event);

		int n_tpc_obj = tpc_object_container_v->size();


		// --------------- Flash Information ---------------
		std::vector<double> largest_flash_v 	= largest_flash_v_v.at(event); // Vec with the largest flash
		largest_flash_y 	= largest_flash_v.at(0);
		largest_flash_z 	= largest_flash_v.at(1);
		largest_flash_time 	= largest_flash_v.at(2);
		largest_flash_pe 	= largest_flash_v.at(3);

		// Only fill for matched events (where each entry is non-zero)
		if (largest_flash_y != 0 && largest_flash_z !=0 && largest_flash_time != 0 && largest_flash_pe != 0){
			h_largest_flash_y	->Fill(largest_flash_y);
			h_largest_flash_z	->Fill(largest_flash_z);
			h_largest_flash_time->Fill(largest_flash_time);
			h_largest_flash_pe	->Fill(largest_flash_pe);
		}
		
		// -------------------------------------------------

		// Loop over TPCObj
		for(int i = 0; i < n_tpc_obj; i++){
			auto const tpc_obj = tpc_object_container_v->at(i);


			// TPC Obj vars
			tpc_obj_vtx_x = tpc_obj.pfpVtxX();
			tpc_obj_vtx_y = tpc_obj.pfpVtxY();
			tpc_obj_vtx_z = tpc_obj.pfpVtxZ();
		   
			const int tpc_obj_mode = tpc_obj.Mode(); 
			n_pfp = tpc_obj.NumPFParticles();
			const int leading_shower_index = GetLeadingShowerIndex(n_pfp, n_tpc_obj, tpc_obj);

			// Other histograms
			n_tracks = 0; n_showers = 0; n_pfp_50Hits = 0; n_tracks_50Hits = 0; n_showers_50Hits = 0;
			GetNumber_Track_Shower(n_pfp, n_tpc_obj, tpc_obj, n_showers, n_tracks, n_pfp_50Hits, n_tracks_50Hits, n_showers_50Hits );// Get the number of tracks and showers
			h_n_pfp			  ->Fill(n_tracks + n_showers);
			h_n_pfp_50Hits	  ->Fill(n_pfp_50Hits);
			h_n_tracks		  ->Fill(n_tracks);
			h_n_tracks_50Hits ->Fill(n_tracks_50Hits);
			h_n_showers		  ->Fill(n_showers);
			h_n_showers_50Hits->Fill(n_showers_50Hits);

			const double Flash_TPCObj_Dist = Flash_TPCObj_vtx_Dist(tpc_obj_vtx_y, tpc_obj_vtx_z, largest_flash_y, largest_flash_z);
			h_Flash_TPCObj_Dist->Fill(Flash_TPCObj_Dist);

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

						const double shower_phi = atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
						h_shower_phi->Fill(shower_phi);

						//  ------------ Leading shower ------------
						if (j == leading_shower_index){
							const double leading_shower_phi 	= atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
							const double leading_shower_theta 	= acos(pfp_obj.pfpDirZ()) * 180 / 3.1415;
							
							const double leading_shower_length 	= pfp_obj.pfpLength();
							const double lonest_track_length 	= GetLongestTrackLength(n_pfp, n_tpc_obj, tpc_obj);

							h_ldg_shwr_hits			->Fill(num_pfp_hits);
							h_ldg_shwr_hits_WPlane	->Fill(pfp_obj.NumPFPHitsW()); 		// W Plane
							h_ldg_shwr_Open_Angle	->Fill(pfp_obj.pfpOpenAngle() * (180 / 3.1415) );
							h_ldg_shwr_dEdx_WPlane	->Fill(pfp_obj.PfpdEdx().at(2) ); 	// W Plane
							h_ldg_shwr_HitPerLen	->Fill(num_pfp_hits / pfp_obj.pfpLength() );
							h_ldg_shwr_Phi			->Fill(leading_shower_phi);
							h_ldg_shwr_Theta		->Fill(leading_shower_theta);
							h_ldg_shwr_CTheta		->Fill(cos(leading_shower_theta * 3.1414 / 180.));
							h_long_Track_ldg_shwr	->Fill(lonest_track_length / leading_shower_length);

							// Vertex Information - require a shower so fill once  when leading shower
							h_tpc_obj_vtx_x->Fill(tpc_obj_vtx_x);
							h_tpc_obj_vtx_y->Fill(tpc_obj_vtx_y);
							h_tpc_obj_vtx_z->Fill(tpc_obj_vtx_z);
						}
						
					}
					// Track like && (Nue || Nuebar)
					if ( pfp_pdg == 13 && (mc_parent_pdg == 12 || mc_parent_pdg == -12) ) {
						const double track_phi 	= atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
						h_track_phi->Fill(track_phi);
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
	DrawTH1D(h_long_Track_ldg_shwr, POT_Scaling);
	DrawTH1D(h_tpc_obj_vtx_x, POT_Scaling);
	DrawTH1D(h_tpc_obj_vtx_y, POT_Scaling);	
	DrawTH1D(h_tpc_obj_vtx_z, POT_Scaling);	
	DrawTH1D(h_n_pfp, POT_Scaling);
	DrawTH1D(h_n_pfp_50Hits, POT_Scaling);
	DrawTH1D(h_n_tracks, POT_Scaling);
	DrawTH1D(h_n_tracks_50Hits, POT_Scaling);
	DrawTH1D(h_n_showers, POT_Scaling);
	DrawTH1D(h_n_showers_50Hits, POT_Scaling);
	DrawTH1D(h_track_phi, POT_Scaling);
	DrawTH1D(h_shower_phi, POT_Scaling);
	DrawTH1D(h_largest_flash_y, POT_Scaling);
	DrawTH1D(h_largest_flash_z, POT_Scaling);
	DrawTH1D(h_largest_flash_time, POT_Scaling);
	DrawTH1D(h_largest_flash_pe, POT_Scaling);
	DrawTH1D(h_Flash_TPCObj_Dist, POT_Scaling);

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
	h_long_Track_ldg_shwr->Write("", TObject::kOverwrite);
	h_tpc_obj_vtx_x->Write("", TObject::kOverwrite);
	h_tpc_obj_vtx_y->Write("", TObject::kOverwrite);
	h_tpc_obj_vtx_z->Write("", TObject::kOverwrite);
	h_n_pfp->Write("", TObject::kOverwrite);
	h_n_pfp_50Hits->Write("", TObject::kOverwrite);
	h_n_tracks->Write("", TObject::kOverwrite);
	h_n_tracks_50Hits->Write("", TObject::kOverwrite);
	h_n_showers->Write("", TObject::kOverwrite);
	h_n_showers_50Hits->Write("", TObject::kOverwrite);
	h_track_phi->Write("", TObject::kOverwrite);
	h_shower_phi->Write("", TObject::kOverwrite);
	h_largest_flash_y->Write("", TObject::kOverwrite);
	h_largest_flash_z->Write("", TObject::kOverwrite);
	h_largest_flash_time->Write("", TObject::kOverwrite);
	h_largest_flash_pe->Write("", TObject::kOverwrite);
	h_Flash_TPCObj_Dist->Write("", TObject::kOverwrite);
	
	f_var_out->Close(); 

} // END MAIN


// Selection cuts--
// loop_flashes
// SetXYflashVector(

// largest_flash_v.at(0) = largest_center_y; const double flash_pe = largest_flash_v.at(0)...
// largest_flash_v.at(1) = largest_center_z;
// largest_flash_v.at(2) = current_event;
// largest_flash_v.at(3) = fOpFlashWidthZ;
// largest_flash_v.at(4) = largest_flash_time;
// largest_flash_v.at(5) = largest_flash;


