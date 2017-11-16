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

bool selection_functions::opt_vtx_distance_width(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double flash_width_z, double tolerance)
{
	const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) ) - flash_width_z;
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
		largest_flash_v.push_back(fOpFlashWidthZ);
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
void selection_functions::VtxTrackNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             double tolerance, std::vector<int> * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		bool close_track = false;
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		//loop over pfparticles in the TPCO
		//this ensures there is at least 1 track in the event we're looking at
		int n_tracks = 0;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			n_tracks++;
		}
		if(n_tracks == 0) {continue; }
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			//check if at least one shower is within the tolerance to tpco vtx
			//if no shower within tolerance, discard tpco (i.e. it's not a cc nue)
			if(pfp_pdg == 13)        //if it's a track
			{
				const double pfp_vtx_x = part.pfpVtxX();
				const double pfp_vtx_y = part.pfpVtxY();
				const double pfp_vtx_z = part.pfpVtxZ();
				close_track = shwr_vtx_distance(tpc_vtx_x, tpc_vtx_y, tpc_vtx_z,
				                                pfp_vtx_x, pfp_vtx_y, pfp_vtx_z, tolerance);
				if(close_track == true)
				{
					if(_verbose) std::cout << " \t " << i << "[Track-To-Nue] \t Passed " << std::endl;
					passed_tpco->at(i) = 1;
					break;
				}
			}
		} //end loop pfps
		if(close_track == 0)//false
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
void selection_functions::OpenAngleCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco,
                                       const double tolerance_open_angle, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		int leading_index = 0;
		int leading_hits  = 0;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			const int n_pfp_hits = part.NumPFPHits();
			if(pfp_pdg == 11 && n_pfp_hits > leading_hits)
			{
				leading_hits = n_pfp_hits;
				leading_index = j;
			}
		}//end loop pfparticles
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		if(leading_open_angle > tolerance_open_angle) {passed_tpco->at(i) = 0; }
	}//end loop tpc objects
}//end open angle cut

//***************************************************************************
//***************************************************************************
void selection_functions::dEdxCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco,
                                  const double tolerance_dedx_min, const double tolerance_dedx_max, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		int leading_index = 0;
		int leading_hits  = 0;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			const int n_pfp_hits = part.NumPFPHits();
			if(pfp_pdg == 11 && n_pfp_hits > leading_hits)
			{
				leading_hits = n_pfp_hits;
				leading_index = j;
			}
		}//end loop pfparticles
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		if(leading_dedx > tolerance_dedx_max || leading_dedx < tolerance_dedx_min) {passed_tpco->at(i) = 0; }
	}//end loop tpc objects
}//end dedx cut
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

std::vector<int> selection_functions::TabulateOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco,
                                                      double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ)
{
	int nue_cc        = 0;
	int nue_cc_qe     = 0;
	int nue_cc_res    = 0;
	int nue_cc_dis    = 0;
	int nue_cc_coh    = 0;
	int nue_cc_mec    = 0;
	int nue_cc_mixed  = 0;
	int nue_cc_out_fv = 0;
	int cosmic        = 0;
	int nue_nc        = 0;
	int numu_cc       = 0;
	int numu_cc_qe    = 0;
	int numu_cc_res   = 0;
	int numu_cc_dis   = 0;
	int numu_cc_coh   = 0;
	int numu_cc_mec   = 0;
	int numu_cc_mixed = 0;
	int numu_nc       = 0;
	int unmatched     = 0;
	int other_mixed   = 0;
	int total         = 0;
	int signal_tpco_num = -1;
	std::vector<int> tabulated_origins;
	tabulated_origins.resize(22);

	bool true_in_tpc = false;

	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		int part_nue_cc    = 0;
		int part_cosmic    = 0;
		int part_nue_nc    = 0;
		int part_numu_cc   = 0;
		int part_numu_nc   = 0;
		int part_unmatched = 0;

		if(passed_tpco->at(i) == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const std::string tpc_obj_origin = tpc_obj.Origin();
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_cc++; }
			if(part.CCNC() == 1 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_nc++; }
			if(part.Origin() == "kBeamNeutrino" && part.CCNC() == 0 && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu_cc++; }
			if(part.Origin() == "kBeamNeutrino" && part.CCNC() == 1 && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu_nc++; }
			if(part.Origin() == "kCosmicRay") { part_cosmic++; }
			if(part.Origin() == "kUnknown")   { part_unmatched++; }
		}
		//now to catagorise the tpco
		if(part_cosmic > 0)
		{
			if(part_nue_cc > 0)  {nue_cc_mixed++; continue; }
			if(part_numu_cc > 0) {numu_cc_mixed++; continue; }
			if(part_nue_nc > 0 || part_numu_nc > 0) {other_mixed++; continue; }
			cosmic++;
		}
		//this uses the true neutrino vertex for this specific event
		//not the true vtx per tpc object - maybe this can be fixed in the future...
		//but using the true nu vtx only matters for the pure signal events,
		//where the neutrino vertex IS the true tpc object vertex
		true_in_tpc = in_fv(vtxX, vtxY, vtxZ,
		                    _x1, _x2, _y1,
		                    _y2, _z1, _z2);
		if(part_cosmic == 0)
		{
			if(part_nue_cc    > 0 && true_in_tpc == false) { nue_cc_out_fv++; continue; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 0)   {nue_cc_qe++;    signal_tpco_num = i;         continue; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 1)   {nue_cc_res++;   signal_tpco_num = i;         continue; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 2)   {nue_cc_dis++;   signal_tpco_num = i;         continue; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 3)   {nue_cc_coh++;   signal_tpco_num = i;         continue; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 10)  {nue_cc_mec++;   signal_tpco_num = i;         continue; }
			if(part_nue_nc    > 0)                        {nue_nc++;           continue; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 0)   {numu_cc_qe++;       continue; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 1)   {numu_cc_res++;      continue; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 2)   {numu_cc_dis++;      continue; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 3)   {numu_cc_coh++;      continue; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 10)  {numu_cc_mec++;      continue; }
			if(part_numu_nc   > 0)                        {numu_nc++;          continue; }
			if(part_unmatched > 0)                        {unmatched++;        continue; }
		}
	}//end loop tpc objects

	nue_cc = nue_cc_qe + nue_cc_res + nue_cc_dis + nue_cc_coh + nue_cc_mec;
	numu_cc = numu_cc_qe + numu_cc_res + numu_cc_dis + numu_cc_coh + numu_cc_mec;
	total = nue_cc + nue_cc_mixed + nue_cc_out_fv + cosmic + nue_nc + numu_cc + numu_cc_mixed + numu_nc + unmatched + other_mixed;

	tabulated_origins.at(0)  = nue_cc;
	tabulated_origins.at(1)  = nue_cc_mixed;
	tabulated_origins.at(2)  = cosmic;
	tabulated_origins.at(3)  = nue_nc;
	tabulated_origins.at(4)  = numu_cc;
	tabulated_origins.at(5)  = unmatched;
	tabulated_origins.at(6)  = other_mixed;
	tabulated_origins.at(7)  = total;
	tabulated_origins.at(8)  = signal_tpco_num;
	tabulated_origins.at(9)  = nue_cc_out_fv;
	tabulated_origins.at(10) = numu_nc;
	tabulated_origins.at(11) = numu_cc_mixed;
	tabulated_origins.at(12) = nue_cc_qe;
	tabulated_origins.at(13) = nue_cc_res;
	tabulated_origins.at(14) = nue_cc_dis;
	tabulated_origins.at(15) = nue_cc_coh;
	tabulated_origins.at(16) = nue_cc_mec;
	tabulated_origins.at(17) = numu_cc_qe;
	tabulated_origins.at(18) = numu_cc_res;
	tabulated_origins.at(19) = numu_cc_dis;
	tabulated_origins.at(20) = numu_cc_coh;
	tabulated_origins.at(21) = numu_cc_mec;
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
                                    int counter_nue_cc_out_fv,
                                    int counter_cosmic,
                                    int counter_nue_nc,
                                    int counter_numu_cc,
                                    int counter_numu_cc_mixed,
                                    int counter_numu_nc,
                                    int counter_unmatched,
                                    int counter_other_mixed,
                                    int counter_nue_cc_qe,
                                    int counter_nue_cc_res,
                                    int counter_nue_cc_dis,
                                    int counter_nue_cc_coh,
                                    int counter_nue_cc_mec,
                                    int counter_numu_cc_qe,
                                    int counter_numu_cc_res,
                                    int counter_numu_cc_dis,
                                    int counter_numu_cc_coh,
                                    int counter_numu_cc_mec,
                                    std::string cut_name)
{
	std::cout << " <" << cut_name << "> " << std::endl;
	std::cout << " Total Candidate Nue     : " << counter << std::endl;
	std::cout << " Number of Nue CC        : " << counter_nue_cc << std::endl;
	std::cout << " Number of Nue CC Mixed  : " << counter_nue_cc_mixed << std::endl;
	std::cout << " Number of Nue CC out FV : " << counter_nue_cc_out_fv << std::endl;
	std::cout << " Number of Cosmic        : " << counter_cosmic << std::endl;
	std::cout << " Number of Nue NC        : " << counter_nue_nc << std::endl;
	std::cout << " Number of Numu CC       : " << counter_numu_cc << std::endl;
	std::cout << " Number of Numu CC Mixed : " << counter_numu_cc_mixed << std::endl;
	std::cout << " Number of Numu NC       : " << counter_numu_nc << std::endl;
	std::cout << " Number of Unmatched     : " << counter_unmatched << std::endl;
	std::cout << " Number of Other Mixed   : " << counter_other_mixed << std::endl;
	std::cout << "---------------------------" << std::endl;
	std::cout << " Nue CC QE               : " << counter_nue_cc_qe   << std::endl;
	std::cout << " Nue CC Res              : " << counter_nue_cc_res  << std::endl;
	std::cout << " Nue CC DIS              : " << counter_nue_cc_dis  << std::endl;
	std::cout << " Nue CC COH              : " << counter_nue_cc_coh  << std::endl;
	std::cout << " Nue CC MEC              : " << counter_nue_cc_mec  << std::endl;
	std::cout << " Numu CC QE              : " << counter_numu_cc_qe  << std::endl;
	std::cout << " Numu CC Res             : " << counter_numu_cc_res << std::endl;
	std::cout << " Numu CC DIS             : " << counter_numu_cc_dis << std::endl;
	std::cout << " Numu CC COH             : " << counter_numu_cc_coh << std::endl;
	std::cout << " Numu CC MEC             : " << counter_numu_cc_mec << std::endl;
	std::cout << "---------------------------" << std::endl;
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

	const double x1 = (det_x1 + _x1);
	const double x2 = (det_x2 - _x2);
	const double y1 = (det_y1 + _y1);
	const double y2 = (det_y2 - _y2);
	const double z1 = (det_z1 + _z1);
	const double z2 = (det_z2 - _z2);

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

void selection_functions::xsec_plot(bool _verbose, double genie_xsec, double xsec, double average_energy, double stat_error)
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

	double g_x[1] = {1.0};//take the average of the interacting energy or the flux energy?
	double g_y[1] = {genie_xsec};
	double g_ex[1] = {0.0};
	double g_ey[1] = {0.0};
	const int g_n = 1;
	TGraphErrors * g_genie_point = new TGraphErrors(g_n, g_x, g_y, g_ex, g_ey);
	g_genie_point->GetXaxis()->SetLimits(0.0, 5.0);//5 GeV
	g_genie_point->SetMinimum(0.0);
	g_genie_point->SetMaximum(50e-39);

	double x[1] = {average_energy};
	double y[1] = {xsec};
	double ex[1] = {stat_error};
	double ey[1] = {stat_error};
	const int n = 1;
	TGraphErrors * g_my_point = new TGraphErrors(n, x, y, ex, ey);
	g_my_point->GetXaxis()->SetLimits(0.0, 5.0);
	g_my_point->SetMinimum(0.0);
	g_my_point->SetMaximum(50e-39);

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
	g_nue_xsec->GetXaxis()->SetLabelOffset(999);
	g_nue_xsec->GetXaxis()->SetLabelSize(0);
	g_nue_xsec->GetYaxis()->SetTitle("#nu_{e} CC Cross Section cm^2");
	g_genie_point->GetXaxis()->SetLabelOffset(999);
	g_genie_point->GetXaxis()->SetLabelSize(0);
	g_genie_point->SetMarkerStyle(3);
	g_genie_point->SetMarkerSize(1.5);
	g_genie_point->Draw("SAMEP");
	g_my_point->GetXaxis()->SetLabelOffset(999);
	g_my_point->GetXaxis()->SetLabelSize(0);
	g_my_point->SetMarkerStyle(3);
	g_my_point->SetMarkerColor(30);
	g_my_point->SetMarkerSize(1.5);
	g_my_point->Draw("SAMEP");
	pad2->Draw();
	pad2->cd();
	h_nue_flux->Draw("Y+");


	combined_c1->Print("combined_xsec.pdf");


	if(f->IsOpen()) {f->Close(); }
	if(f->IsOpen()) {xsec_f->Close(); }

}

void selection_functions::PostCutPlots(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                       std::vector<int> * passed_tpco, bool _verbose,
                                       TH2I * h_tracks_showers, TH2I * h_tracks_showers_cosmic, TH2I * h_tracks_showers_numu,
                                       TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_mixed,
                                       TH1D * h_leading_shower_open_angle_numu_cc, TH1D * h_leading_shower_open_angle_numu_nc,
                                       TH1D * h_leading_shower_open_angle_cosmic, TH1D * h_leading_shower_open_angle_nue_nc,
                                       TH1D * h_leading_shower_open_angle_numu_cc_mixed, TH1D * h_leading_shower_open_angle_other_mixed,
                                       TH1D * h_leading_shower_open_angle_unmatched,
                                       TH1D * h_trk_vtx_dist_nue_cc, TH1D * h_trk_vtx_dist_nue_cc_mixed,
                                       TH1D * h_trk_vtx_dist_numu_cc, TH1D * h_trk_vtx_dist_numu_nc,
                                       TH1D * h_trk_vtx_dist_cosmic, TH1D * h_trk_vtx_dist_nue_nc,
                                       TH1D * h_trk_vtx_dist_numu_cc_mixed, TH1D * h_trk_vtx_dist_other_mixed,
                                       TH1D * h_trk_vtx_dist_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	int nue_cc        = 0;
	int nue_cc_mixed  = 0;
	int cosmic        = 0;
	int nue_nc        = 0;
	int numu_cc       = 0;
	int numu_cc_mixed = 0;
	int numu_nc       = 0;
	int unmatched     = 0;
	int other_mixed   = 0;
	int total = 0;
	int signal_tpco_num = -1;

	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i) == 0) {continue; }
		int part_nue_cc    = 0;
		int part_cosmic    = 0;
		int part_nue_nc    = 0;
		int part_numu_cc   = 0;
		int part_numu_nc   = 0;
		int part_unmatched = 0;

		int leading_index = 0;
		int most_hits = 0;

		int num_tracks = 0;
		int num_showers = 0;

		double smallest_trk_vtx_dist = 1000;
		bool has_track = false;

		auto const tpc_obj = tpc_object_container_v->at(i);
		const std::string tpc_obj_origin = tpc_obj.Origin();
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			const int n_pfp_hits = part.NumPFPHits();

			if(pfp_pdg == 13) //if it's a track
			{
				has_track = true;
				const double pfp_vtx_x = part.pfpVtxX();
				const double pfp_vtx_y = part.pfpVtxY();
				const double pfp_vtx_z = part.pfpVtxZ();
				const double trk_vtx_dist = sqrt(((tpc_vtx_x - pfp_vtx_x) * (tpc_vtx_x - pfp_vtx_x)) +
				                                 ((tpc_vtx_y - pfp_vtx_y) * (tpc_vtx_y - pfp_vtx_y)) +
				                                 ((tpc_vtx_z - pfp_vtx_z) * (tpc_vtx_z - pfp_vtx_z)));
				if(trk_vtx_dist < smallest_trk_vtx_dist) {smallest_trk_vtx_dist = trk_vtx_dist; }
				const double trk_length = part.pfpLength();
				//std::cout << trk_length << std::endl;
			}

			if(pfp_pdg == 11 && n_pfp_hits > most_hits) {leading_index = j; most_hits = n_pfp_hits; }
			if(part.PFParticlePdgCode() == 11) {num_showers++; }
			if(part.PFParticlePdgCode() == 13) {num_tracks++; }
			if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_cc++; }
			if(part.CCNC() == 1 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_nc++; }
			if(part.Origin() == "kBeamNeutrino" && part.CCNC() == 0 && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu_cc++; }
			if(part.Origin() == "kBeamNeutrino" && part.CCNC() == 1 && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu_nc++; }
			if(part.Origin() == "kCosmicRay") { part_cosmic++; }
			if(part.Origin() == "kUnknown")   { part_unmatched++; }
		}//end pfp loop
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);

		if(part_cosmic > 0)
		{
			if(part_nue_cc > 0)
			{
				nue_cc_mixed++;
				h_leading_shower_open_angle_nue_cc_mixed->Fill(leading_open_angle);
				if(has_track) {h_trk_vtx_dist_nue_cc_mixed->Fill(smallest_trk_vtx_dist); }
				continue;
			}
			if(part_numu_cc > 0)
			{
				numu_cc_mixed++;
				h_leading_shower_open_angle_numu_cc_mixed->Fill(leading_open_angle);
				if(has_track) {h_trk_vtx_dist_numu_cc_mixed->Fill(smallest_trk_vtx_dist); }
				continue;
			}
			if(part_nue_nc > 0 || part_numu_nc > 0)
			{
				other_mixed++;
				h_leading_shower_open_angle_other_mixed->Fill(leading_open_angle);
				if(has_track) {h_trk_vtx_dist_other_mixed->Fill(smallest_trk_vtx_dist); }
				continue;
			}
			cosmic++;
			h_tracks_showers_cosmic->Fill(num_tracks, num_showers);
			h_leading_shower_open_angle_cosmic->Fill(leading_open_angle);
			if(has_track) {h_trk_vtx_dist_cosmic->Fill(smallest_trk_vtx_dist); }
		}
		if(part_cosmic == 0)
		{
			if(part_nue_cc    > 0)
			{
				nue_cc++;
				h_tracks_showers->Fill(num_tracks, num_showers);
				h_leading_shower_open_angle_nue_cc->Fill(leading_open_angle);
				if(has_track) {h_trk_vtx_dist_nue_cc->Fill(smallest_trk_vtx_dist); }
				continue;
			}
			if(part_nue_nc    > 0)
			{
				nue_nc++;
				h_leading_shower_open_angle_nue_nc->Fill(leading_open_angle);
				if(has_track) {h_trk_vtx_dist_nue_nc->Fill(smallest_trk_vtx_dist); }
				continue;
			}
			if(part_numu_cc   > 0)
			{
				numu_cc++;
				h_tracks_showers_numu->Fill(num_tracks, num_showers);
				h_leading_shower_open_angle_numu_cc->Fill(leading_open_angle);
				if(has_track) {h_trk_vtx_dist_numu_cc->Fill(smallest_trk_vtx_dist); }
				continue;
			}
			if(part_numu_nc   > 0)
			{
				numu_nc++;
				h_leading_shower_open_angle_numu_nc->Fill(leading_open_angle);
				if(has_track) {h_trk_vtx_dist_numu_nc->Fill(smallest_trk_vtx_dist); }
				continue;
			}
			if(part_unmatched > 0)
			{
				unmatched++;
				h_leading_shower_open_angle_unmatched->Fill(leading_open_angle);
				if(has_track) {h_trk_vtx_dist_unmatched->Fill(smallest_trk_vtx_dist); }
				continue;
			}
		}
	}//end tpco loop
}

void selection_functions::TopologyPlots(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<int> * passed_tpco,
                                        double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                                        TH2D * h_pfp_track_shower_nue_cc_qe,
                                        TH2D * h_pfp_track_shower_nue_cc_out_fv,
                                        TH2D * h_pfp_track_shower_nue_cc_res,
                                        TH2D * h_pfp_track_shower_nue_cc_dis,
                                        TH2D * h_pfp_track_shower_nue_cc_coh,
                                        TH2D * h_pfp_track_shower_nue_cc_mec,
                                        TH2D * h_pfp_track_shower_nue_nc,
                                        TH2D * h_pfp_track_shower_numu_cc_qe,
                                        TH2D * h_pfp_track_shower_numu_cc_res,
                                        TH2D * h_pfp_track_shower_numu_cc_dis,
                                        TH2D * h_pfp_track_shower_numu_cc_coh,
                                        TH2D * h_pfp_track_shower_numu_cc_mec,
                                        TH2D * h_pfp_track_shower_numu_nc,
                                        TH2D * h_pfp_track_shower_nue_cc_mixed,
                                        TH2D * h_pfp_track_shower_numu_cc_mixed,
                                        TH2D * h_pfp_track_shower_cosmic,
                                        TH2D * h_pfp_track_shower_other_mixed,
                                        TH2D * h_pfp_track_shower_unmatched)
{
	bool true_in_tpc = false;

	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		int part_nue_cc    = 0;
		int part_cosmic    = 0;
		int part_nue_nc    = 0;
		int part_numu_cc   = 0;
		int part_numu_nc   = 0;
		int part_unmatched = 0;

		if(passed_tpco->at(i) == 0) {continue; }
		bool tpco_id_valid = false;
		auto const tpc_obj = tpc_object_container_v->at(i);
		const std::string tpc_obj_origin = tpc_obj.Origin();
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks  = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id;
		//loop over pfparticles in the TPCO
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_cc++; }
			if(part.CCNC() == 1 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_nc++; }
			if(part.Origin() == "kBeamNeutrino" && part.CCNC() == 0 && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu_cc++; }
			if(part.Origin() == "kBeamNeutrino" && part.CCNC() == 1 && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu_nc++; }
			if(part.Origin() == "kCosmicRay") { part_cosmic++; }
			if(part.Origin() == "kUnknown")   { part_unmatched++; }
		}
		//now to catagorise the tpco
		if(part_cosmic > 0)
		{
			if(part_nue_cc > 0)                                               { tpco_id = "nue_cc_mixed";  tpco_id_valid = true; }
			if(part_numu_cc > 0 && tpco_id_valid == false)                    { tpco_id = "numu_cc_mixed"; tpco_id_valid = true; }
			if((part_nue_nc > 0 || part_numu_nc > 0) && tpco_id_valid == false) { tpco_id = "other_mixed";   tpco_id_valid = true; }
			if(tpco_id_valid == false)                                        { tpco_id = "cosmic";        tpco_id_valid = true; }
		}
		//this uses the true neutrino vertex for this specific event
		//not the true vtx per tpc object - maybe this can be fixed in the future...
		//but using the true nu vtx only matters for the pure signal events,
		//where the neutrino vertex IS the true tpc object vertex
		true_in_tpc = in_fv(vtxX, vtxY, vtxZ,
		                    _x1, _x2, _y1,
		                    _y2, _z1, _z2);
		if(part_cosmic == 0)
		{
			if(part_nue_cc    > 0 && true_in_tpc == false && tpco_id_valid == false) { tpco_id = "nue_cc_out_fv";    tpco_id_valid = true; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 0    && tpco_id_valid == false) { tpco_id = "nue_cc_qe";        tpco_id_valid = true; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 1    && tpco_id_valid == false) { tpco_id = "nue_cc_res";       tpco_id_valid = true; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 2    && tpco_id_valid == false) { tpco_id = "nue_cc_dis";       tpco_id_valid = true; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 3    && tpco_id_valid == false) { tpco_id = "nue_cc_coh";       tpco_id_valid = true; }
			if(part_nue_cc    > 0 && tpc_obj_mode == 10   && tpco_id_valid == false) { tpco_id = "nue_cc_mec";       tpco_id_valid = true; }
			if(part_nue_nc    > 0                         && tpco_id_valid == false) { tpco_id = "nue_nc";           tpco_id_valid = true; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 0    && tpco_id_valid == false) { tpco_id = "numu_cc_qe";       tpco_id_valid = true; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 1    && tpco_id_valid == false) { tpco_id = "numu_cc_res";      tpco_id_valid = true; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 2    && tpco_id_valid == false) { tpco_id = "numu_cc_dis";      tpco_id_valid = true; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 3    && tpco_id_valid == false) { tpco_id = "numu_cc_coh";      tpco_id_valid = true; }
			if(part_numu_cc   > 0 && tpc_obj_mode == 10   && tpco_id_valid == false) { tpco_id = "numu_cc_mec";      tpco_id_valid = true; }
			if(part_numu_nc   > 0                         && tpco_id_valid == false) { tpco_id = "numu_nc";          tpco_id_valid = true; }
			if(part_unmatched > 0                         && tpco_id_valid == false) { tpco_id = "unmatched";        tpco_id_valid = true; }
		}
		if(tpco_id == "nue_cc_qe")        {h_pfp_track_shower_nue_cc_qe->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "nue_cc_out_fv")    {h_pfp_track_shower_nue_cc_out_fv->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "nue_cc_res")       {h_pfp_track_shower_nue_cc_res->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "nue_cc_dis")       {h_pfp_track_shower_nue_cc_dis->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "nue_cc_coh")       {h_pfp_track_shower_nue_cc_coh->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "nue_cc_mec")       {h_pfp_track_shower_nue_cc_mec->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "nue_nc")           {h_pfp_track_shower_nue_nc->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "numu_cc_qe")       {h_pfp_track_shower_numu_cc_qe->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "numu_cc_res")      {h_pfp_track_shower_numu_cc_res->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "numu_cc_dis")      {h_pfp_track_shower_numu_cc_dis->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "numu_cc_coh")      {h_pfp_track_shower_numu_cc_coh->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "numu_cc_mec")      {h_pfp_track_shower_numu_cc_mec->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "numu_nc")          {h_pfp_track_shower_numu_nc->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "nue_cc_mixed")     {h_pfp_track_shower_nue_cc_mixed->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "numu_cc_mixed")    {h_pfp_track_shower_numu_cc_mixed->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "cosmic")           {h_pfp_track_shower_cosmic->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "other_mixed")      {h_pfp_track_shower_other_mixed->Fill(n_pfp_tracks, n_pfp_showers); }
		if(tpco_id == "unmatched")        {h_pfp_track_shower_unmatched->Fill(n_pfp_tracks, n_pfp_showers); }
	}        //end loop tpc objects

}//end function
