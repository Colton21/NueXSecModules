#include "selection_functions.h"

//***************************************************************************
//***************************************************************************

int selection_functions::flash_in_time(int flash_pe, int flash_pe_threshold, double flash_time, double flash_start, double flash_end)
{
	if(flash_time >= flash_start && flash_time <= flash_end)
	{
		if(flash_pe >= flash_pe_threshold)
		{
			return 1;//true
		}
	}
	return 0;//false
}

//***************************************************************************
//***************************************************************************

void selection_functions::loop_flashes(TFile * f, TTree * optical_tree, int flash_pe_threshold, double flash_time_start,
                                       double flash_time_end, std::vector<int> * _passed_runs)
{

	optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;

	optical_tree->SetBranchAddress("Run", &fRun);
	optical_tree->SetBranchAddress("Event", &fEvent);
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE);
	optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime);

	const int optical_entries = optical_tree->GetEntries();

	/// I want passed_runs to be the same size as the number of TPCObject_v
	/// This will give the vector a 1:1 mapping to the events

	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;

	std::cout << "Optical Entries: " << optical_entries << std::endl;

	for(int i = 0; i < optical_entries; i++)
	{
		optical_tree->GetEntry(i);

		//this function here is meant to construct a vector mapping to which
		//events successfully pass this cut
		current_run = fRun;
		current_event = fEvent;

		int in_time = 0;//false
		if(_passed_runs->size() != 0)
		{
			if(current_event == last_event && current_run == last_run)
			{
				if(_passed_runs->back() == 1)//true
				{
					last_event = current_event;
					last_run = current_run;
					continue;
				}
				if(_passed_runs->back() == 0)//false
				{
					_passed_runs->pop_back();
				}
			}
		}
		last_event = current_event;
		last_run = current_run;

		in_time = flash_in_time(fOpFlashPE, flash_pe_threshold, fOpFlashTime, flash_time_start, flash_time_end);
		_passed_runs->push_back(in_time);
	}//end loop optical entries

}//end optical loop function



//***************************************************************************
//***************************************************************************

bool selection_functions::in_fv(double x, double y, double z,
                                double x1, double x2, double y1,
                                double y2, double z1, double z2)
{
	const double det_x1 = 0;
	const double det_x2 = 256.35;
	const double det_y1 = -116.5;
	const double det_y2 = 116.5;
	const double det_z1 = 0;
	const double det_z2 = 1036.8;

	if(x <= det_x1 + x1 || x >= det_x2 - x2) {return false; }
	if(y <= det_y1 + y1 || y >= det_y2 - y2) {return false; }
	if(z <= det_z1 + z1 || z >= det_z2 - z2) {return false; }
	else{return true; }
}

//***************************************************************************
//***************************************************************************

void selection_functions::fiducial_volume_cut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                              std::vector<int> * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	if(passed_tpco->size() != n_tpc_obj) {std::cout << "Passed TPCO Vector Size != nTPCO!" << std::endl; }
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		//this vertex (tpcobject container) is the vertex of the pfp neutrino
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const bool InFV = in_fv(tpc_vtx_x, tpc_vtx_y, tpc_vtx_z, _x1, _x2, _y1, _y2, _z1, _z2);
		if(InFV == 1)//true
		{
			passed_tpco->at(i) = 1;
			if(_verbose) {std::cout << " \t " << i << "[Fid Volume Cut] \t Passed" << std::endl; }
		}
		if(InFV == 0)//false
		{
			passed_tpco->at(i) = 0;
		}
	}//end tpc object loop
}

//***************************************************************************
//***************************************************************************

bool selection_functions::opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance)
{
	const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
	if(distance <= tolerance) {return true; }
	return false;
}

//***************************************************************************
//***************************************************************************

void selection_functions::SetXYflashVector(TFile * f, TTree * optical_tree, std::vector< std::vector< double> > * largest_flash_v_v)
{
	optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");
	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashWidthY = 0;
	double fOpFlashWidthZ = 0;
	double fOpFlashCenterY = 0;
	double fOpFlashCenterZ = 0;

	optical_tree->SetBranchAddress("Run", &fRun);
	optical_tree->SetBranchAddress("Event", &fEvent);
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE);
	optical_tree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY);
	optical_tree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ);
	optical_tree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY);
	optical_tree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ);

	const int optical_entries = optical_tree->GetEntries();
	double largest_flash = 0;
	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;
	//since I only care about the largest flash in each event,
	//size(optical_entries)-->size(total_entries)

	for(int i = 0; i < optical_entries; i++)
	{
		optical_tree->GetEntry(i);
		current_run = fRun;
		current_event = fEvent;
		if(current_event != last_event) {largest_flash = 0; }
		double this_flash = fOpFlashPE;
		std::vector < double > largest_flash_v;//contains the y,z for largest flash

		if(this_flash <= largest_flash)
		{
			last_event = current_event;
			last_run = current_run;
			continue;
		}
		if(current_event == last_event && current_run == last_run)
		{
			if(this_flash > largest_flash)
			{
				largest_flash = this_flash;
				largest_flash_v_v->pop_back();
			}
		}
		last_event = current_event;
		last_run = current_run;
		largest_flash = this_flash;
		largest_flash_v.push_back(fOpFlashCenterY);
		largest_flash_v.push_back(fOpFlashCenterZ);
		largest_flash_v.push_back(current_event);
		largest_flash_v_v->push_back(largest_flash_v);
		largest_flash_v.clear();
	}
}

//***************************************************************************
//***************************************************************************

void selection_functions::flashRecoVtxDist(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                           double tolerance, std::vector<int> * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const bool is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), tolerance);
		if(is_close == 1)//true
		{
			passed_tpco->at(i) = 1;
			if(_verbose) std::cout << " \t " << i << "[Vertex-To-Flash] \t Passed " << std::endl;
		}
		if(is_close == 0)//false
		{
			passed_tpco->at(i) = 0;
		}
	}        //end tpc object loop
}//end flashRecoVtxDist

//***************************************************************************
//***************************************************************************

bool selection_functions::shwr_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                                            double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z, double tolerance)
{
	const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );
	if(distance <= tolerance) {return true; }
	return false;
}

//***************************************************************************
//***************************************************************************
//this function wants to remove particles too far from the reconstructed neutrino vertex
void selection_functions::VtxNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                        double tolerance, std::vector<int> * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		bool close_shower = false;
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			//check if at least one shower is within the tolerance to tpco vtx
			//if no shower within tolerance, discard tpco (i.e. it's not a cc nue)
			if(pfp_pdg == 11)        //if it's a shower
			{
				const double pfp_vtx_x = part.pfpVtxX();
				const double pfp_vtx_y = part.pfpVtxY();
				const double pfp_vtx_z = part.pfpVtxZ();
				close_shower = shwr_vtx_distance(tpc_vtx_x, tpc_vtx_y, tpc_vtx_z,
				                                 pfp_vtx_x, pfp_vtx_y, pfp_vtx_z, tolerance);
				if(close_shower == true)
				{
					if(_verbose) std::cout << " \t " << i << "[Shower-To-Nue] \t Passed " << std::endl;
					passed_tpco->at(i) = 1;
					break;
				}
			}
		} //end loop pfps
		if(close_shower == 0)//false
		{
			passed_tpco->at(i) = 0;
		}
	}//end tpc object loop
}

//***************************************************************************
//***************************************************************************
//this function wants to remove particles too far from the reconstructed neutrino vertex
void selection_functions::HitThreshold(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                       double threshold, std::vector<int> * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		bool over_threshold = false;
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_tpco_hits = tpc_obj.NumPFPHits();
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			//check if at least one shower is within the tolerance to tpco vtx
			//if no shower within tolerance, discard tpco (i.e. it's not a cc nue)
			if(pfp_pdg == 11)        //if it's a shower
			{
				const int n_pfp_hits = part.NumPFPHits();
				if(n_pfp_hits >= threshold) {over_threshold = true; }
				if(over_threshold == true)
				{
					if(_verbose) std::cout << " \t " << i << "[Hit Threshold] \t Passed " << std::endl;
					passed_tpco->at(i) = 1;
					break;
				}
			}
		} //end loop pfps
		if(over_threshold == 0)//false
		{
			passed_tpco->at(i) = 0;
		}
	}//end tpc object loop
}

//***************************************************************************
//***************************************************************************
//this gives a list of all of the origins of the tpc objects
void selection_functions::GetOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::string> * tpco_origin_v)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		auto const tpc_obj = tpc_object_container_v->at(i);
		const std::string tpc_obj_origin = tpc_obj.Origin();
		tpco_origin_v->push_back(tpc_obj_origin);
	}
}

//***************************************************************************
//***************************************************************************
//this function simply checks if the tpc object is a nue
void selection_functions::HasNue(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();

		bool has_nue = false;
		for(int j = 0; j <n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			if(pfp_pdg == 12)
			{
				has_nue = true;
				break;
			}
		}//end loop pfparticles
		if(has_nue == false)
		{
			passed_tpco->at(i) = 0;
		}
		if(has_nue == true)
		{
			passed_tpco->at(i) = 1;
			if(_verbose) std::cout << " \t " << i << "[Reco Nue Cut] \t Passed" << std::endl;
		}
	}
}

//***************************************************************************
//***************************************************************************
//this function just counts if at least 1 tpc object passes the cuts
bool selection_functions::ValidTPCObjects(std::vector<int> * passed_tpco)
{
	const int n_tpco = passed_tpco->size();
	int pass_sum = 0;
	for(int const tpco : * passed_tpco)
	{
		pass_sum = pass_sum + tpco;
	}
	//0 --> failed tpco, 1 --> passed tpco
	if(pass_sum == 0) {return false; }
	return true;
}
//***************************************************************************
//***************************************************************************

std::vector<int> selection_functions::TabulateOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco)
{
	int nue_cc = 0;
	int nue_cc_mixed = 0;
	int cosmic = 0;
	int nue_nc = 0;
	int numu   = 0;
	int unmatched = 0;
	int other_mixed = 0;
	int total = 0;
	std::vector<int> tabulated_origins;
	tabulated_origins.resize(8);

	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		int part_nue_cc = 0;
		int part_cosmic = 0;
		int part_nue_nc = 0;
		int part_numu   = 0;
		int part_unmatched = 0;

		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const std::string tpc_obj_origin = tpc_obj.Origin();
		if(tpc_obj_origin == "kCosmicRay") {cosmic++; }
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_cc++; }
			if(part.CCNC() == 1 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_nc++; }
			if(part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu++; }
			if(part.Origin() == "kCosmicRay") { part_cosmic++; }
			if(part.Origin() == "kUnknown")   { part_unmatched++; }
		}
		//now to catagorise the tpco
		if(part_cosmic > 0)
		{
			if(part_nue_cc > 0) {nue_cc_mixed++; continue; }
			if(part_nue_nc > 0 || part_numu > 0) {other_mixed++; continue; }
			cosmic++;
		}
		if(part_cosmic == 0)
		{
			if(part_nue_cc    > 0) {nue_cc++;    continue; }
			if(part_nue_nc    > 0) {nue_nc++;    continue; }
			if(part_numu      > 0) {numu++;      continue; }
			if(part_unmatched > 0) {unmatched++; continue; }
		}
	}

	total = nue_cc + nue_cc_mixed + cosmic + nue_nc + numu + unmatched + other_mixed;

	tabulated_origins.at(0) = nue_cc;
	tabulated_origins.at(1) = nue_cc_mixed;
	tabulated_origins.at(2) = cosmic;
	tabulated_origins.at(3) = nue_nc;
	tabulated_origins.at(4) = numu;
	tabulated_origins.at(5) = unmatched;
	tabulated_origins.at(6) = other_mixed;
	tabulated_origins.at(7) = total;
	return tabulated_origins;
}

//***************************************************************************
//***************************************************************************
//modify this so it takes a string of the cut name so I only pass it a few variable at a time,
//then I can call this function several times later at the bottom
void selection_functions::PrintInfo(int mc_nue_cc_counter,
                                    int counter,
                                    int counter_nue_cc,
                                    int counter_nue_cc_mixed,
                                    int counter_cosmic,
                                    int counter_nue_nc,
                                    int counter_numu,
                                    int counter_unmatched,
                                    int counter_other_mixed,
                                    std::string cut_name)
{
	std::cout << " <" << cut_name << "> " << std::endl;
	std::cout << " Total Candidate Nue    : " << counter << std::endl;
	std::cout << " Number of Nue CC       : " << counter_nue_cc << std::endl;
	std::cout << " Number of Nue CC Mixed : " << counter_nue_cc_mixed << std::endl;
	std::cout << " Number of Cosmic       : " << counter_cosmic << std::endl;
	std::cout << " Number of Nue NC       : " << counter_nue_nc << std::endl;
	std::cout << " Number of Numu         : " << counter_numu << std::endl;
	std::cout << " Number of Unmatched    : " << counter_unmatched << std::endl;
	std::cout << " Number of Other Mixed  : " << counter_other_mixed << std::endl;
	std::cout << "------------------------" << std::endl;
	const double efficiency = double(counter_nue_cc) / double(mc_nue_cc_counter);
	const double purity = double(counter_nue_cc) / double(counter);
	std::cout << " Efficiency       : " << efficiency << std::endl;
	std::cout << " Purity           : " << purity << std::endl;
	std::cout << "------------------------" << std::endl;
	std::cout << "------------------------" << std::endl;
}
//***************************************************************************
//***************************************************************************
double selection_functions::calcNumNucleons(double _x1, double _x2, double _y1,
                                            double _y2, double _z1, double _z2)
{
	const double det_x1 = 0;
	const double det_x2 = 256.35;
	const double det_y1 = -116.5;
	const double det_y2 = 116.5;
	const double det_z1 = 0;
	const double det_z2 = 1036.8;

	const double x1 = det_x1 + x1;
	const double x2 = det_x2 - x2;
	const double y1 = det_y1 + y1;
	const double y2 = det_y2 - y2;
	const double z1 = det_z1 + z1;
	const double z2 = det_z2 - z2;

	const double vol = (x2 - x1) * (y2 - y1) * (z2 - z1); //cm^3

	const double lar_density = 0.0013954; //kg/cm^3
	const double au = 1.67 * pow(10, -27); //kg
	const double n_target = vol * lar_density / au;
	return n_target;
}
//***************************************************************************
//***************************************************************************
void selection_functions::calcXSec(double _x1, double _x2, double _y1,
                                   double _y2, double _z1, double _z2,
                                   int n_total, int n_bkg, double flux, double efficiency, std::vector<double>  * xsec_cc)
{
	const int n_events = n_total - n_bkg;
	//scale_factor = 2.4 * math.pow(10, 17)  # POT / nue
	//calculate the number of nucleons based on the fiducial volume
	const double n_target = calcNumNucleons(_x1, _x2, _y1,
	                                        _y2, _z1, _z2);

	std::cout << "-------------------" << std::endl;
	std::cout << "N_total    :  " << n_total << std::endl;
	std::cout << "N_bkg      :  " << n_bkg << std::endl;
	std::cout << "N_target   :  " << n_target << std::endl;
	std::cout << "Flux       :  " << flux << std::endl;
	std::cout << "Efficiency :  " << efficiency << std::endl;
	std::cout << "-------------------" << std::endl;
	xsec_cc->push_back((n_events) / (flux * n_target * efficiency));
	const double n_error = (n_total / sqrt(n_total));
	xsec_cc->push_back((n_error) /  (flux * n_target * efficiency));
	const double sys_error = xsec_cc->at(0) * 0.30; //beam sys error of 30%
	xsec_cc->push_back(sys_error);
}

//***************************************************************************
//***************************************************************************

void selection_functions::xsec_plot(bool _verbose, double genie_xsec, double xsec)
{

	//setting verbose manually for this function...
	_verbose = false;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "Cross Section Plotting Function" << std::endl;
	std::cout << "-------------------------------" << std::endl;

	std::cout << "Opening NuMIFlux.root" << std::endl;
	//first we need to open the root file
	TFile * f = new TFile("../arxiv/NuMIFlux.root");
	if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; exit(1); }
	TH1D * h_nue_flux = (TH1D*)f->Get("nueFluxHisto");
	h_nue_flux->GetXaxis()->SetLimits(0.0, 5.0);//5 Gev

	std::cout << "Opening argon_xsec_nue.root" << std::endl;
	//this is a genie file with the cross section
	TFile * xsec_f = new TFile("../arxiv/argon_xsec_nue.root");
	if(!xsec_f->IsOpen()) {std::cout << "Could not open file!" << std::endl; exit(1); }
	TGraph * g_nue_xsec = (TGraph*)xsec_f->Get("nu_e_Ar40/tot_cc");
	g_nue_xsec->GetXaxis()->SetLimits(0.0, 5.0);//5 GeV
	g_nue_xsec->SetMinimum(0.0);
	g_nue_xsec->SetMaximum(50e-39);

	// double x[1] = {1.0};//take the average of the interacting energy or the flux energy?
	// double y[1] = {genie_xsec};
	// const int n = 1;
	//TGraph * g_genie_point = new TGraph(n, x, y);
	// g_genie_point->GetXaxis()->SetLimits(0.0, 5.0);//5 GeV
	// g_genie_point->SetMinimum(0.0);
	// g_genie_point->SetMaximum(50e-39);
	TMarker m_genie_point(1.0, genie_xsec, 1);
	TMarker m_my_point(1.0, xsec, 1);

	const int n_bins = h_nue_flux->GetNbinsX();
	double bin_flux_sum = 0;
	double bin_interaction_sum = 0;

	//std::cout << "Loop over energy bins of 50 MeV" << std::endl;
	for(int bin = 0; bin < n_bins; bin++)
	{
		const double bin_flux = h_nue_flux->GetBinContent(bin);
		bin_flux_sum = bin_flux_sum + bin_flux;

		const double bin_energy = h_nue_flux->GetBinCenter(bin);
		//g_nue_xsec->GetY()[bin] *= 1e-38 / 40;
		//g_nue_xsec->GetY() [bin] = g_nue_xsec->Eval(bin_energy) * 1e-38 / 40;
		//std::cout <<  bin_energy << " , " << g_nue_xsec->GetY() [bin] << std::endl;
		const double bin_xsec_val = g_nue_xsec->Eval(bin_energy) * 1e-38 / 40;//this gets the units per nucleon per cm^2
		const double bin_interactions = bin_xsec_val * bin_flux;
		bin_interaction_sum = bin_interaction_sum + bin_interactions;

		if(_verbose == true && bin_flux > 0.0)
		{
			std::cout << "Cross Section: " << bin_xsec_val << std::endl;
			std::cout << "Bin Energy: " << bin_energy << '\t' << "Bin Flux: " << bin_flux << std::endl;
			std::cout << "===================" << std::endl;
			std::cout << "Interactions: " << bin_interactions << std::endl;
			std::cout << "===================" << std::endl;
		}
	}
	for(int bin = 0; bin < g_nue_xsec->GetN(); bin++)
	{
		g_nue_xsec->GetY() [bin] *= 1e-38 / 40;
	}
	if(_verbose)
	{
		std::cout << "Total Neutrinos: " << bin_flux_sum << std::endl;
		std::cout << "Total Interactions: " << bin_interaction_sum << std::endl;
	}
	//the website: https://cdcvs.fnal.gov/redmine/projects/ubooneoffline/wiki/NuMI_Flux_Histograms
	//states that the flux is scaled to 6e20 POT
	const double pot_used = 6 * pow(10, 20);

	if(_verbose)
	{
		std::cout << "Total POT used: " << pot_used << std::endl;
		std::cout << "Final Result is: " << bin_flux_sum / pot_used << " nues per POT per cm^2" << std::endl;
		std::cout << "SUM{Phi(E) * Sigma(E)} / SUM{Phi(E)} = " << bin_interaction_sum / bin_flux_sum << " cm^2" << std::endl;
	}

	h_nue_flux->SetStats(kFALSE);
	h_nue_flux->SetTitle("");
	g_nue_xsec->SetTitle("");
	TCanvas * combined_c1 = new TCanvas();
	combined_c1->cd();
	TPad * pad1 = new TPad();
	TPad * pad2 = new TPad();
	pad2->SetFillStyle(4000); //will be transparent
	pad2->SetFrameFillStyle(0);
	pad1->Draw();
	pad1->cd();
	g_nue_xsec->Draw("");
	m_genie_point.Draw("SAME");
	m_my_point.Draw("SAME");
	pad2->Draw();
	pad2->cd();
	h_nue_flux->Draw("Y+");


	combined_c1->Print("combined_xsec.pdf");


	if(f->IsOpen()) {f->Close(); }
	if(f->IsOpen()) {xsec_f->Close(); }

}

void selection_functions::PostCutPlots(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                       std::vector<int> * passed_tpco, bool _verbose, TH2I * h_tracks_showers)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		int num_tracks = 0;
		int num_showers = 0;
		auto const tpc_obj = tpc_object_container_v->at(i);
		const std::string tpc_obj_origin = tpc_obj.Origin();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			if(part.PFParticlePdgCode() == 11) {num_showers++; }
			if(part.PFParticlePdgCode() == 13) {num_tracks++; }
		}
		h_tracks_showers->Fill(num_tracks, num_showers);
	}
}
