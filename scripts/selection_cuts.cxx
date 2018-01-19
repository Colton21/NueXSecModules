#include "selection_cuts.h"

//***************************************************************************
bool selection_cuts::flash_in_time(double flash_time, double flash_start, double flash_end)
{
	if(flash_time >= flash_start && flash_time <= flash_end)
	{
		return true;//true
	}
	return false;//false
}
//***************************************************************************
bool selection_cuts::flash_pe(int flash_pe, int flash_pe_threshold)
{
	if(flash_pe >= flash_pe_threshold)
	{
		return true;//true
	}
	return false;//false
}
//***************************************************************************
void selection_cuts::loop_flashes(TFile * f, TTree * optical_tree, int flash_pe_threshold, double flash_time_start,
                                  double flash_time_end, std::vector<int> * _passed_runs)
{
	optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;

	optical_tree->SetBranchAddress("Run",              &fRun);
	optical_tree->SetBranchAddress("Event",            &fEvent);
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE);
	optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime);

	const int optical_entries = optical_tree->GetEntries();

	/// I want passed_runs to be the same size as the number of TPCObject_v
	/// This will give the vector a 1:1 mapping to the events

	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;

	//contains the entry number for a given OpFlash per event
	std::vector<int> optical_list_pe;
	std::vector<std::vector<int> > optical_list_pe_v;
	std::vector<double> optical_list_time;
	std::vector<std::vector<double> > optical_list_time_v;

	std::cout << "Optical Entries: " << optical_entries << std::endl;

	for(int i = 0; i < optical_entries; i++)
	{
		optical_tree->GetEntry(i);

		//this function here is meant to construct a vector mapping to which
		//events successfully pass this cut
		current_run = fRun;
		current_event = fEvent;

		//new event
		if(current_event != last_event)
		{
			optical_list_pe.clear();
			optical_list_time.clear();
			optical_list_pe.push_back(fOpFlashPE);
			optical_list_time.push_back(fOpFlashTime);
		}
		//same event
		if(current_event == last_event && current_run == last_run)
		{
			optical_list_pe_v.pop_back();
			optical_list_time_v.pop_back();
			optical_list_pe.push_back(fOpFlashPE);
			optical_list_time.push_back(fOpFlashTime);
		}
		last_event = current_event;
		last_run = current_run;
		optical_list_pe_v.push_back(optical_list_pe);
		optical_list_time_v.push_back(optical_list_time);
	}
	std::cout << "Optical List Vector Size: " << optical_list_pe_v.size() << std::endl;

	std::vector<int> optical_pass_list_v;
	int opt_list_counter = 0;
	for(int i = 0; i < optical_list_pe_v.size(); i++)
	{
		bool in_time = false;
		bool got_in_time = false;
		bool sufficient_flash = false;
		bool got_sufficient_flash = false;
		auto const opt_time_v = optical_list_time_v.at(i);
		auto const opt_pe_v = optical_list_pe_v.at(i);
		for(int j = 0; j < optical_list_pe_v.at(i).size(); j++)
		{
			auto const opt_time = opt_time_v.at(j);
			auto const opt_pe = opt_pe_v.at(j);
			in_time = flash_in_time(opt_time, flash_time_start, flash_time_end);
			if(in_time == true) {got_in_time = true; }
			sufficient_flash = flash_pe(opt_pe, flash_pe_threshold);
			if(sufficient_flash == true) {got_sufficient_flash = true; }
			//flash is both in time and over PE threshold
			if(in_time == true && sufficient_flash == true)
			{
				_passed_runs->at(i) = 1;
				break;
			}
		}
		//flash in time, PE too low
		if(got_sufficient_flash == false && got_in_time == true) {_passed_runs->at(i) = 2; }
		//no flash was in time
		if(got_in_time == false) {_passed_runs->at(i) = 0; }
	}
}//end optical loop function
//***************************************************************************
bool selection_cuts::in_fv(double x, double y, double z,
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


	//we also want to consider the dead region in the detector
	//let's expand our FV to exclude the region of wires 7136-7263
	//Since it's collection plane --> (7136-4800) * 0.3 = 700.8 cm start
	// (7263-4800) * 0.3 = 738.9 cm end
	// const double dead_z_start = 700.8;
	// const double dead_z_end = 738.9;
	// const double dead_tolerance = 10; //cm
	// if( z >= dead_z_start - dead_tolerance && z <= dead_z_end + dead_tolerance) {return false; }
	//did not have a positive effect ...

	else{return true; }
}
//***************************************************************************
void selection_cuts::fiducial_volume_cut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                         std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	if(passed_tpco->size() != n_tpc_obj) {std::cout << "Passed TPCO Vector Size != nTPCO!" << std::endl; }
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		//this vertex (tpcobject container) is the vertex of the pfp neutrino
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const bool InFV = in_fv(tpc_vtx_x, tpc_vtx_y, tpc_vtx_z, _x1, _x2, _y1, _y2, _z1, _z2);
		if(InFV == 1)//true
		{
			passed_tpco->at(i).first = 1;
			if(_verbose) {std::cout << " \t " << i << "[Fid Volume Cut] \t Passed" << std::endl; }
		}
		if(InFV == 0)//false
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "InFV";
		}
	}//end tpc object loop
}
//***************************************************************************
bool selection_cuts::opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance)
{
	const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
	//const double distance = tpc_vtx_z - flash_vtx_z;
	if(distance <= tolerance) {return true; }
	return false;
}
//***************************************************************************
bool selection_cuts::opt_vtx_distance_width(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z,
                                            double flash_width_z, double tolerance)
{
	const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) ) - flash_width_z;
	if(distance <= tolerance) {return true; }
	return false;
}
//***************************************************************************
void selection_cuts::SetXYflashVector(TFile * f, TTree * optical_tree, std::vector< std::vector< double> > * largest_flash_v_v,
                                      double flash_time_start, double flash_time_end)
{
	optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");
	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;
	double fOpFlashWidthY = 0;
	double fOpFlashWidthZ = 0;
	double fOpFlashCenterY = 0;
	double fOpFlashCenterZ = 0;

	optical_tree->SetBranchAddress("Run",              &fRun           );
	optical_tree->SetBranchAddress("Event",            &fEvent         );
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE     );
	optical_tree->SetBranchAddress("OpFlashTime",     &fOpFlashTime   );
	optical_tree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY );
	optical_tree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ );
	optical_tree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY);
	optical_tree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ);


	const int optical_entries = optical_tree->GetEntries();

	/// I want passed_runs to be the same size as the number of TPCObject_v
	/// This will give the vector a 1:1 mapping to the events

	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;
	int largest_flash = 0;

	//contains the entry number for a given OpFlash per event
	std::vector<int> largest_flash_list;
	std::vector<double> largest_flash_v;
	largest_flash_v.resize(6);

	for(int i = 0; i < optical_entries; i++)
	{
		optical_tree->GetEntry(i);
		bool in_time = false;
		if(fOpFlashTime >= flash_time_start && fOpFlashTime <= flash_time_end) {in_time = true; }

		//this function here is meant to construct a vector mapping to which
		//events successfully pass this cut
		current_run = fRun;
		current_event = fEvent;

		if(in_time == false)
		{
			largest_flash = 0;
			largest_flash_v.at(0) = fOpFlashCenterY;
			largest_flash_v.at(1) = fOpFlashCenterZ;
			largest_flash_v.at(2) = current_event;
			largest_flash_v.at(3) = fOpFlashWidthZ;
			largest_flash_v.at(4) = fOpFlashTime;
			largest_flash_v.at(5) = largest_flash;
			if(current_event != last_event) {largest_flash_v_v->push_back(largest_flash_v); }
		}

		//new event
		if(current_event != last_event && in_time == true)
		{
			largest_flash = fOpFlashPE;
			largest_flash_list.push_back(largest_flash);
			largest_flash_v.at(0) = fOpFlashCenterY;
			largest_flash_v.at(1) = fOpFlashCenterZ;
			largest_flash_v.at(2) = current_event;
			largest_flash_v.at(3) = fOpFlashWidthZ;
			largest_flash_v.at(4) = fOpFlashTime;
			largest_flash_v.at(5) = largest_flash;
			largest_flash_v_v->push_back(largest_flash_v);
		}
		//same event
		if(current_event == last_event && current_run == last_run)
		{
			if(fOpFlashPE > largest_flash && in_time == true)
			{
				largest_flash = fOpFlashPE;
				largest_flash_list.pop_back();
				largest_flash_v_v->pop_back();
				largest_flash_list.push_back(largest_flash);
				largest_flash_v.at(0) = fOpFlashCenterY;
				largest_flash_v.at(1) = fOpFlashCenterZ;
				largest_flash_v.at(2) = current_event;
				largest_flash_v.at(3) = fOpFlashWidthZ;
				largest_flash_v.at(4) = fOpFlashTime;
				largest_flash_v.at(5) = largest_flash;
				largest_flash_v_v->push_back(largest_flash_v);
			}
		}
		last_event = current_event;
		last_run = current_run;
	}//end loop optical entries
	int flash_counter = 0;
	for(auto const flash_v : * largest_flash_v_v)
	{
		//std::cout << "Largest Flash in this event: " << flash_v.at(5) << std::endl;
		if(flash_v.at(5) >= 50) {flash_counter++; }
	}
	std::cout << "Non-zero Largest Flashes: " << flash_counter << std::endl;
}
//***************************************************************************
void selection_cuts::flashRecoVtxDist(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      double tolerance, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const bool is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), tolerance);
		if(is_close == 1)//true
		{
			passed_tpco->at(i).first = 1;
			if(_verbose) std::cout << " \t " << i << "[Vertex-To-Flash] \t Passed " << std::endl;
		}
		if(is_close == 0)//false
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "FlashDist";
		}
	}        //end tpc object loop
}//end flashRecoVtxDist
//***************************************************************************
bool selection_cuts::shwr_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z, double tolerance)
{
	const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );
	if(distance <= tolerance) {return true; }
	return false;
}
//***************************************************************************
void selection_cuts::VtxNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                   double tolerance, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		bool close_shower = false;
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		//loop over pfparticles in the TPCO
		const int n_tracks = tpc_obj.NPfpTracks();
		//if(n_tracks == 0) {continue; } //this is here for testing how the cut responds when I only cut for !=0 tracks
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
					passed_tpco->at(i).first = 1;
					break;
				}
			}
		} //end loop pfps
		if(close_shower == 0)//false
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "ShwrVtx";
		}
	}//end tpc object loop
}
//***************************************************************************
void selection_cuts::VtxTrackNuDistance(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                        double tolerance, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		bool close_track = false;
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		//loop over pfparticles in the TPCO
		//this ensures there is at least 1 track in the event we're looking at
		const int n_tracks = tpc_obj.NPfpTracks();
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
					passed_tpco->at(i).first = 1;
					break;
				}
			}
		} //end loop pfps
		if(close_track == 0)//false
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "TrkVtx";
		}
	}//end tpc object loop
}
//***************************************************************************
void selection_cuts::HitThreshold(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                  double threshold, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
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
					passed_tpco->at(i).first = 1;
					break;
				}
			}
		} //end loop pfps
		if(over_threshold == 0)//false
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "HitThreshold";
		}
	}//end tpc object loop
}
//***************************************************************************
void selection_cuts::HasNue(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                            std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
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
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "HasNue";
		}
		if(has_nue == true)
		{
			//* this might cause some problems later - is it doing anything??? *//
			passed_tpco->at(i).first = 1;
			if(_verbose) std::cout << " \t " << i << "[Reco Nue Cut] \t Passed" << std::endl;
		}
	}
}
//***************************************************************************
void selection_cuts::OpenAngleCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
                                  const std::vector<double> tolerance_open_angle, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
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
		if(leading_open_angle > tolerance_open_angle.at(1) || leading_open_angle < tolerance_open_angle.at(0))
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "OpenAngle";
		}
	}//end loop tpc objects
}//end open angle cut
//***************************************************************************
void selection_cuts::dEdxCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
                             const double tolerance_dedx_min, const double tolerance_dedx_max, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
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
		if(leading_dedx > tolerance_dedx_max || leading_dedx < tolerance_dedx_min)
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "dEdX";
		}
	}//end loop tpc objects
}//end dedx cut
//***************************************************************************
void selection_cuts::SecondaryShowersDistCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, const double dist_tolerance)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		const double tpco_vtx_x = tpc_obj.pfpVtxX();
		const double tpco_vtx_y = tpc_obj.pfpVtxY();
		const double tpco_vtx_z = tpc_obj.pfpVtxZ();
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
		if(n_pfp_showers > 3)
		{
			for(int j = 0; j < n_pfp; j++)
			{
				if(j == leading_index) {continue; }//we assume leading shower == electron shower
				auto const part = tpc_obj.GetParticle(j);
				const int pfp_pdg = part.PFParticlePdgCode();
				const double pfp_vtx_x = part.pfpVtxX();
				const double pfp_vtx_y = part.pfpVtxY();
				const double pfp_vtx_z = part.pfpVtxZ();
				const double distance = sqrt(pow((pfp_vtx_x - tpco_vtx_x),2) +
				                             pow((pfp_vtx_y - tpco_vtx_y),2) +
				                             pow((pfp_vtx_z - tpco_vtx_z),2));
				if(pfp_pdg == 11)//22 cm is ~ 2 radiation lengths
				{
					if(distance > dist_tolerance)
					{
						passed_tpco->at(i).first = 0;
						passed_tpco->at(i).second = "SecondaryDist";
						if(_verbose) {std::cout << "[SecondaryDist] TPC Object Failed!" << std::endl; }
						break;
					}
				}
			}//end loop pfp
		}
	}
}
//***************************************************************************
void selection_cuts::HitLengthRatioCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                       const double pfp_hits_length_tolerance)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
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
		const int pfp_pdg = leading_shower.PFParticlePdgCode();
		const double pfp_hits = leading_shower.NumPFPHits();
		const double pfp_length = leading_shower.pfpLength();
		const double pfp_hits_length_ratio = (pfp_hits / pfp_length);

		if(pfp_pdg == 11)
		{
			if(pfp_hits_length_ratio < pfp_hits_length_tolerance)
			{
				passed_tpco->at(i).first = 0;
				passed_tpco->at(i).second = "HitLengthRatio";
				if(_verbose) {std::cout << "[HitLengthRatio] TPC Object Failed!" << std::endl; }
			}
		}
	}
}
//***************************************************************************
//this gives a list of all of the origins of the tpc objects
void selection_cuts::GetOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::string> * tpco_origin_v)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		auto const tpc_obj = tpc_object_container_v->at(i);
		const std::string tpc_obj_origin = tpc_obj.Origin();
		tpco_origin_v->push_back(tpc_obj_origin);
	}
}
