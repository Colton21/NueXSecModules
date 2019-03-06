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
// Driver function to sort the vector elements
// by second element of pairs
bool selection_cuts::sortbysec(const std::tuple<int,int,int> &a, const std::tuple<int,int,int> &b)
{
	return (std::get<0>(a) < std::get<0>(b));
}
//***************************************************************************
void selection_cuts::loop_flashes(TFile * f, TTree * optical_tree, int flash_pe_threshold, double flash_time_start,
                                  double flash_time_end, std::vector<int> * _passed_runs, std::vector<std::pair<double, int> > * flash_time, const int stream)
{
	optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

	int fRun = 0;
	int fSubrun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;

	optical_tree->SetBranchAddress("Run",              &fRun);
	optical_tree->SetBranchAddress("SubRun",           &fSubrun);
	optical_tree->SetBranchAddress("Event",            &fEvent);
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE);
	optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime);

	const int optical_entries = optical_tree->GetEntries();

	/// I want passed_runs to be the same size as the number of TPCObject_v
	/// This will give the vector a 1:1 mapping to the events

	int current_event = 0;
	int current_run = 0;
	int current_subrun = 0;
	int last_event = 0;
	int last_run = 0;
	int last_subrun = 0;

	//contains the entry number for a given OpFlash per event
	std::vector<int> optical_list_pe;
	std::vector<std::vector<int> > optical_list_pe_v;
	std::vector<double> optical_list_time;
	std::vector<std::vector<double> > optical_list_time_v;

	std::cout << "Optical Entries: " << optical_entries << std::endl;

	//std::vector<std::tuple<int, int, int> > temp_v;

	for(int i = 0; i < optical_entries; i++)
	{
		optical_tree->GetEntry(i);

		//this function here is meant to construct a vector mapping to which
		//events successfully pass this cut
		current_run = fRun;
		current_subrun = fSubrun;
		current_event = fEvent;
		double op_flash_time = fOpFlashTime;

		//const int stream values:
		// 0 = data (on-beam)
		// 1 = ext  (off-beam)
		// 2 = mc

		//Note: EXT and On-Beam triggers are shifted
		//this shifts the EXT to the On-Beam timing
		if(stream == 1) {op_flash_time = op_flash_time - 0.343; }
		//if(ext) {op_flash_time = op_flash_time - 0.406; }

		//in this case the MC is shifted compared to the on-beam data
		//this is an artefact of the MC simulation - should be corrected eventually
		//in the future.
		//based on comparisons, shift is ~1 us
		if(stream == 2) {op_flash_time = op_flash_time + 1.0; }

		auto const my_pair = std::make_pair(op_flash_time, current_run);
		flash_time->push_back(my_pair);
		// auto const tuple = std::make_tuple(current_run, current_subrun, current_event);
		// temp_v.push_back(tuple);

		///////////////////
		// Using sort() function to sort by 2nd element
		// of pair
		//std::sort(flash_time->begin(), flash_time->end(), sortbysec);

		//new event
		if(current_event != last_event)
		{
			optical_list_pe.clear();
			optical_list_time.clear();
			optical_list_pe.push_back(fOpFlashPE);
			optical_list_time.push_back(op_flash_time);
		}
		//same event
		if(current_event == last_event && current_run == last_run)
		{
			optical_list_pe_v.pop_back();
			optical_list_time_v.pop_back();
			optical_list_pe.push_back(fOpFlashPE);
			optical_list_time.push_back(op_flash_time);
		}
		last_event = current_event;
		last_run = current_run;
		last_subrun = current_subrun;

		optical_list_pe_v.push_back(optical_list_pe);
		optical_list_time_v.push_back(optical_list_time);
	}
	std::cout << "Optical List Vector Size: " << optical_list_pe_v.size() << std::endl;
	//
	// std::sort(temp_v.begin(), temp_v.end(), sortbysec);
	// const int last_run = 0;
	//
	// for(auto const tuple : temp_v)
	// {
	//      const int run = std::get<0>(tuple);
	//      const int subrun = std::get<1>(tuple);
	//      const int event = std::get<2>(tuple);
	//      //std::cout << "Run: " << run << " Subrun: " << subrun << " Event: " << event << std::endl;
	// }

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
bool selection_cuts::in_fv(double x, double y, double z, std::vector<double> fv_boundary_v)
{
	const double det_x1 = 0;
	const double det_x2 = 256.35;
	const double det_y1 = -116.5;
	const double det_y2 = 116.5;
	const double det_z1 = 0;
	const double det_z2 = 1036.8;

	const double x1 = fv_boundary_v.at(0);
	const double x2 = fv_boundary_v.at(1);
	const double y1 = fv_boundary_v.at(2);
	const double y2 = fv_boundary_v.at(3);
	const double z1 = fv_boundary_v.at(4);
	const double z2 = fv_boundary_v.at(5);

	if(x <= det_x1 + x1 || x >= det_x2 - x2) {return false; }
	if(y <= det_y1 + y1 || y >= det_y2 - y2) {return false; }
	if(z <= det_z1 + z1 || z >= det_z2 - z2) {return false; }
	return true;

	//we also want to consider the dead region in the detector
	//let's expand our FV to exclude the region of wires 7136-7263
	//Since it's collection plane --> (7136-4800) * 0.3 = 700.8 cm start
	// (7263-4800) * 0.3 = 738.9 cm end
	// const double dead_z_start = 700.8;
	// const double dead_z_end = 738.9;
	// const double dead_tolerance = 25; //cm
	// if( z >= dead_z_start - dead_tolerance && z <= dead_z_end + dead_tolerance) {return false; }
	//did not have a positive effect ...
}
//***************************************************************************
void selection_cuts::fiducial_volume_cut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                         std::vector<double> fv_boundary_v,
                                         std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		//this vertex (tpcobject container) is the vertex of the pfp neutrino
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		//std::cout << tpc_vtx_x << ", " << tpc_vtx_y << ", " << tpc_vtx_z << std::endl;
		const bool InFV = in_fv(tpc_vtx_x, tpc_vtx_y, tpc_vtx_z, fv_boundary_v);
		if(InFV == true) { if(_verbose) { std::cout << " \t " << i << "[Fid Volume Cut] \t Passed" << std::endl; } }
		if(InFV == false)
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
                                      double flash_time_start, double flash_time_end, double flash_pe_threshold)
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
	optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime   );
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

	//contains the entry number for a given OpFlash per event
	std::vector<int> optical_list_pe;
	std::vector<std::vector<int> > optical_list_pe_v;
	std::vector<double> optical_list_time;
	std::vector<std::vector<double> > optical_list_time_v;
	std::vector<double> optical_list_flash_center_y;
	std::vector<std::vector<double> > optical_list_flash_center_y_v;
	std::vector<double> optical_list_flash_center_z;
	std::vector<std::vector<double> > optical_list_flash_center_z_v;
	std::vector<std::vector<double> > optical_list_flash_time_v;
	std::vector<double> optical_list_flash_time;
	std::vector<double> optical_list_currentevent;
	std::vector<std::vector<double> > optical_list_currentevent_v;
	// largest_flash_v.at(3) = fOpFlashWidthZ;
	// largest_flash_v.at(4) = fOpFlashTime;
	// largest_flash_v.at(5) = largest_flash;

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
			optical_list_flash_center_y.clear();
			optical_list_flash_center_z.clear();
			optical_list_flash_time.clear();
			optical_list_currentevent.clear();
			optical_list_pe.push_back(fOpFlashPE);
			optical_list_time.push_back(fOpFlashTime);
			optical_list_flash_center_y.push_back(fOpFlashCenterY);
			optical_list_flash_center_z.push_back(fOpFlashCenterZ);
			optical_list_flash_time.push_back(fOpFlashTime);
			optical_list_currentevent.push_back(fEvent);
		}
		//same event
		if(current_event == last_event && current_run == last_run)
		{
			optical_list_pe_v.pop_back();
			optical_list_time_v.pop_back();
			optical_list_flash_center_y_v.pop_back();
			optical_list_flash_center_z_v.pop_back();
			optical_list_flash_time_v.pop_back();
			optical_list_currentevent_v.pop_back();
			optical_list_pe.push_back(fOpFlashPE);
			optical_list_time.push_back(fOpFlashTime);
			optical_list_flash_center_y.push_back(fOpFlashCenterY);
			optical_list_flash_center_z.push_back(fOpFlashCenterZ);
			optical_list_flash_time.push_back(fOpFlashTime);
			optical_list_currentevent.push_back(fEvent);
		}
		last_event = current_event;
		last_run = current_run;
		optical_list_pe_v.push_back(optical_list_pe);
		optical_list_time_v.push_back(optical_list_time);
		optical_list_flash_center_y_v.push_back(optical_list_flash_center_y);
		optical_list_flash_center_z_v.push_back(optical_list_flash_center_z);
		optical_list_flash_time_v.push_back(optical_list_flash_time);
		optical_list_currentevent_v.push_back(optical_list_currentevent);
	}
	std::cout << "Optical List Vector Size: " << optical_list_pe_v.size() << std::endl;

	std::vector<int> optical_pass_list_v;
	int opt_list_counter = 0;
	std::vector<double> largest_flash_v;
	largest_flash_v.resize(7, 0);
	//loop through all events
	for(int i = 0; i < optical_list_pe_v.size(); i++)
	{
		bool in_time = false;
		bool got_in_time = false;
		bool sufficient_flash = false;
		bool got_sufficient_flash = false;
		auto const opt_time_v = optical_list_time_v.at(i);
		auto const opt_pe_v = optical_list_pe_v.at(i);
		//loop through all flashes in event
		double largest_flash = 0.;
		double largest_center_y = 0;
		double largest_center_z = 0;
		double largest_flash_time = 0;

		for(int j = 0; j < optical_list_pe_v.at(i).size(); j++)
		{
			auto const opt_time         = opt_time_v.at(j);
			auto const opt_pe           = opt_pe_v.at(j);
			const double opt_center_y   = optical_list_flash_center_y_v.at(i).at(j);
			const double opt_center_z   = optical_list_flash_center_z_v.at(i).at(j);
			const double opt_flash_time = optical_list_flash_time_v.at(i).at(j);
			in_time = flash_in_time(opt_time, flash_time_start, flash_time_end);
			if(in_time == true) {got_in_time = true; }
			sufficient_flash = flash_pe(opt_pe, flash_pe_threshold);//update to flash_pe_threshold
			if(sufficient_flash == true) {got_sufficient_flash = true; }
			//flash is both in time and over PE threshold
			if(in_time == true && sufficient_flash == true)
			{
				//let's find the largest flash in this event
				if(opt_pe > largest_flash)
				{
					largest_flash      = opt_pe;
					largest_center_y   = opt_center_y;
					largest_center_z   = opt_center_z;
					largest_flash_time = opt_flash_time;
				}
			}
		}
		largest_flash_v.at(0) = largest_center_y;
		largest_flash_v.at(1) = largest_center_z;
		// largest_flash_v.at(2) = current_event;
		// largest_flash_v.at(3) = fOpFlashWidthZ;
		largest_flash_v.at(4) = largest_flash_time;
		largest_flash_v.at(5) = largest_flash; //PE
		largest_flash_v_v->push_back(largest_flash_v);
	}

	//
	// //contains the entry number for a given OpFlash per event
	// std::vector<int> largest_flash_list;
	//
	//
	// for(int i = 0; i < optical_entries; i++)
	// {
	//      optical_tree->GetEntry(i);
	//      bool in_time = false;
	//      if(fOpFlashTime >= flash_time_start && fOpFlashTime <= flash_time_end) {in_time = true; }
	//
	//      //this function here is meant to construct a vector mapping to which
	//      //events successfully pass this cut
	//      current_run = fRun;
	//      current_event = fEvent;
	//
	//      if(in_time == false)
	//      {
	//              largest_flash = 0;
	//              largest_flash_v.at(0) = fOpFlashCenterY;
	//              largest_flash_v.at(1) = fOpFlashCenterZ;
	//              largest_flash_v.at(2) = current_event;
	//              largest_flash_v.at(3) = fOpFlashWidthZ;
	//              largest_flash_v.at(4) = fOpFlashTime;
	//              largest_flash_v.at(5) = largest_flash;
	//              largest_flash_v.at(6) = 0; //not in time
	//              if(current_event != last_event) {largest_flash_v_v->push_back(largest_flash_v); }
	//      }
	//
	//      //new event
	//      if(current_event != last_event && in_time == true)
	//      {
	//              largest_flash = fOpFlashPE;
	//              largest_flash_list.push_back(largest_flash);
	//              largest_flash_v.at(0) = fOpFlashCenterY;
	//              largest_flash_v.at(1) = fOpFlashCenterZ;
	//              largest_flash_v.at(2) = current_event;
	//              largest_flash_v.at(3) = fOpFlashWidthZ;
	//              largest_flash_v.at(4) = fOpFlashTime;
	//              largest_flash_v.at(5) = largest_flash;
	//              largest_flash_v.at(6) = 1;
	//              largest_flash_v_v->push_back(largest_flash_v);
	//      }
	//      //same event
	//      if(current_event == last_event && current_run == last_run)
	//      {
	//              if(fOpFlashPE > largest_flash && in_time == true)
	//              {
	//                      largest_flash = fOpFlashPE;
	//                      largest_flash_list.pop_back();
	//                      largest_flash_v_v->pop_back();
	//                      largest_flash_list.push_back(largest_flash);
	//                      largest_flash_v.at(0) = fOpFlashCenterY;
	//                      largest_flash_v.at(1) = fOpFlashCenterZ;
	//                      largest_flash_v.at(2) = current_event;
	//                      largest_flash_v.at(3) = fOpFlashWidthZ;
	//                      largest_flash_v.at(4) = fOpFlashTime;
	//                      largest_flash_v.at(5) = largest_flash;
	//                      largest_flash_v.at(6) = 1;
	//                      largest_flash_v_v->push_back(largest_flash_v);
	//              }
	//      }
	//      last_event = current_event;
	//      last_run = current_run;
	// }//end loop optical entries
	int flash_counter = 0;
	int too_small_flash_counter = 0;
	int out_time_counter = 0;
	for(auto const flash_v : * largest_flash_v_v)
	{
		//std::cout << "Largest Flash in this event: " << flash_v.at(5) << std::endl;
		if(flash_v.at(5) >= 50) {flash_counter++; }
		if(flash_v.at(5) == 0) {out_time_counter++; }
	}
	// std::cout << "Largest Flashes >= 50 PE              : " << flash_counter << std::endl;
	// std::cout << "Largest Flashes Out-Time or Too Small : " << out_time_counter << std::endl;
	// std::cout << "Largest Flash Vector Size             : " << largest_flash_v_v->size() << std::endl;
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
		bool is_close;
		//flash is upstream
		if(tpc_vtx_z < largest_flash_v.at(1))
		{
			is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), tolerance);
		}
		//flash is downstream
		if(tpc_vtx_z >= largest_flash_v.at(1))
		{
			is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), (tolerance - 20));
		}
		if(is_close == 1)//true
		{
			//passed_tpco->at(i).first = 1;
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
					//passed_tpco->at(i).first = 1;
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
					//passed_tpco->at(i).first = 1;
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
					//passed_tpco->at(i).first = 1;
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
void selection_cuts::HitThresholdCollection(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                            double threshold, std::vector<std::pair<int, std::string> > * passed_tpco, const bool _verbose)
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
		const int n_pfp_hits_w = leading_shower.NumPFPHitsW();
		if(n_pfp_hits_w < threshold)
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "HitThresholdW";
		}
		if(n_pfp_hits_w >= threshold)
		{
			if(_verbose) {std::cout << "[Collection Hit Threshold]" << '\t' << "Passed" << std::endl; }
		}
	}

	//      bool over_threshold = false;
	//      auto const tpc_obj = tpc_object_container_v->at(i);
	//      const int n_pfp = tpc_obj.NumPFParticles();
	//      //loop over pfparticles in the TPCO
	//      for(int j = 0; j < n_pfp; j++)
	//      {
	//              auto const part = tpc_obj.GetParticle(j);
	//              const int pfp_pdg = part.PFParticlePdgCode();
	//              //check if at least one shower is within the tolerance to tpco vtx
	//              //if no shower within tolerance, discard tpco (i.e. it's not a cc nue)
	//              if(pfp_pdg == 11)        //if it's a shower
	//              {
	//                      const int n_pfp_hits_w = part.NumPFPHitsW();
	//                      if(n_pfp_hits_w >= threshold) {over_threshold = true; }
	//                      if(over_threshold == true)
	//                      {
	//                              if(_verbose) std::cout << " \t " << i << "[Hit Threshold] \t Passed " << std::endl;
	//                              //passed_tpco->at(i).first = 1;
	//                              break;
	//                      }
	//              }
	//      } //end loop pfps
	//      if(over_threshold == false)//false
	//      {
	//              passed_tpco->at(i).first = 0;
	//              passed_tpco->at(i).second = "HitThresholdW";
	//      }
	// }//end tpc object loop
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
		bool has_valid_shower = false;

		//test case -- difference is very small following all selection cuts
		//if(tpc_obj.NPfpShowers() != 0) {has_nue = true; has_valid_shower = true; }

		for(int j = 0; j <n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			const int pfp_hits = part.NumPFPHits();
			if(pfp_pdg == 11 && pfp_hits > 0) { has_valid_shower = true; }
			if(pfp_pdg == 12) { has_nue = true; }
		}//end loop pfparticles
		if(has_nue == false || has_valid_shower == false)
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "HasNue";
		}
		if(has_nue == true && has_valid_shower == true) { if(_verbose) std::cout << " \t " << i << "[Reco Nue Cut] \t Passed" << std::endl; }
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
		if(leading_open_angle <= tolerance_open_angle.at(1) && leading_open_angle >= tolerance_open_angle.at(0))
		{
			if(_verbose) {std::cout << "[Open Angle]" << '\t' << "Passed" << std::endl; }
		}
	}//end loop tpc objects
}//end open angle cut
//***************************************************************************
void selection_cuts::dEdxCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, std::vector<std::pair<int, std::string> > * passed_tpco,
                             const double tolerance_dedx_min, const double tolerance_dedx_max, const bool _verbose, const bool is_ext)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }

		auto const tpc_obj = tpc_object_container_v->at(i);
		const bool is_data = tpc_obj.IsData();
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
		double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		if(is_data && !is_ext) {leading_dedx = leading_dedx * (242.72 / 196.979); }
		if(leading_dedx > tolerance_dedx_max || leading_dedx < tolerance_dedx_min)
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "dEdX";
		}
		if(leading_dedx <= tolerance_dedx_max && leading_dedx >= tolerance_dedx_min)
		{
			if(_verbose) {std::cout << "[dE/dx]" << '\t' << "Passed" << std::endl; }
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
		if(n_pfp_showers <= 1) {continue; }
		//This cut does not target events with fewer than 2 showers
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
		for(int j = 0; j < n_pfp; j++)
		{
			if(j == leading_index) {continue; }        //we assume leading shower == electron shower
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			const double pfp_vtx_x = part.pfpVtxX();
			const double pfp_vtx_y = part.pfpVtxY();
			const double pfp_vtx_z = part.pfpVtxZ();
			const double distance = sqrt(pow((pfp_vtx_x - tpco_vtx_x),2) +
			                             pow((pfp_vtx_y - tpco_vtx_y),2) +
			                             pow((pfp_vtx_z - tpco_vtx_z),2));
			if(pfp_pdg == 11)        //22 cm is ~ 2 radiation lengths
			{
				if(distance > dist_tolerance)
				{
					passed_tpco->at(i).first = 0;
					passed_tpco->at(i).second = "SecondaryDist";
					if(_verbose) {std::cout << "[SecondaryDist] TPC Object Failed!" << std::endl; }
					break;
				}
				if(distance <= dist_tolerance)
				{
					if(_verbose) {std::cout << "[2nd Shower]" << '\t' << "Passed" << std::endl; }
				}
			}
		}        //end loop pfp
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
void selection_cuts::LongestTrackLeadingShowerCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                  const double ratio_tolerance)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		if(n_pfp_tracks == 0) {continue; }
		int leading_index = 0;
		int leading_hits  = 0;
		double longest_track = 0;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const pfp = tpc_obj.GetParticle(j);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			const int n_pfp_hits = pfp.NumPFPHits();
			if(pfp_pdg == 11 && n_pfp_hits > leading_hits)
			{
				leading_hits = n_pfp_hits;
				leading_index = j;
			}
			if(pfp_pdg == 13)
			{
				const double trk_length = pfp.pfpLength();
				if(trk_length > longest_track)
				{
					longest_track = trk_length;
				}
			}
		}//end loop pfparticles
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_length = leading_shower.pfpLength();
		const double longest_track_leading_shower_ratio = longest_track / leading_shower_length;

		//if the ratio is too large:
		if(longest_track_leading_shower_ratio > ratio_tolerance)
		{
			passed_tpco->at(i).first = 0;
			passed_tpco->at(i).second = "TrkShwrLenRatio";
			if(_verbose) {std::cout << "[TrkShwrLenhRatio] TPC Object Failed!" << std::endl; }
		}
	}
}
//***************************************************************************
bool selection_cuts::IsContained(std::vector<double> track_start, std::vector<double> track_end, std::vector<double> fv_boundary_v)
{
	if(in_fv(track_start.at(0), track_start.at(1), track_start.at(2), fv_boundary_v) == true
	   && in_fv(track_end.at(0), track_end.at(1), track_end.at(2), fv_boundary_v) == true)
	{
		return true;
	}
	return false;
}
//***************************************************************************
void selection_cuts::ContainedTracksCut(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                        std::vector<double> fv_boundary_v, const bool enabled)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(enabled == false) {continue; }
		if(passed_tpco->at(i).first == 0) {continue; }

		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		//this is normally enabled, but due to test cut below it is off
		if(n_pfp_tracks == 0) {continue; }

		for(int j = 0; j < n_pfp; j++)
		{
			auto const pfp = tpc_obj.GetParticle(j);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg == 13)
			{
				const double pfp_vtx_x = pfp.pfpVtxX();
				const double pfp_vtx_y = pfp.pfpVtxY();
				const double pfp_vtx_z = pfp.pfpVtxZ();
				const double pfp_dir_x = pfp.pfpDirX();
				const double pfp_dir_y = pfp.pfpDirY();
				const double pfp_dir_z = pfp.pfpDirZ();
				const double trk_length = pfp.pfpLength();
				const double pfp_end_x = (pfp.pfpVtxX() + (trk_length * pfp_dir_x));
				const double pfp_end_y = (pfp.pfpVtxY() + (trk_length * pfp_dir_y));
				const double pfp_end_z = (pfp.pfpVtxZ() + (trk_length * pfp_dir_z));

				std::vector<double> pfp_start_vtx {pfp_vtx_x, pfp_vtx_y, pfp_vtx_z};
				std::vector<double> pfp_end_vtx {pfp_end_x, pfp_end_y, pfp_end_z};

				const bool is_contained = IsContained(pfp_start_vtx, pfp_end_vtx, fv_boundary_v);

				//if not contained
				if(is_contained == false)
				{
					passed_tpco->at(i).first = 0;
					passed_tpco->at(i).second = "Contained";
					if(_verbose) {std::cout << "[Containment] TPC Object Failed!" << std::endl; }
					break;
				}
			}//end is track
		}//end loop pfparticles
		 //adding temporary cut on the shower direction here
		 // const int n_pfp_showers = tpc_obj.NPfpShowers();
		 // int leading_index = 0;
		 // int leading_hits  = 0;
		 // for(int j = 0; j < n_pfp; j++)
		 // {
		 //     auto const part = tpc_obj.GetParticle(j);
		 //     const int pfp_pdg = part.PFParticlePdgCode();
		 //     const int n_pfp_hits = part.NumPFPHits();
		 //     if(pfp_pdg == 11 && n_pfp_hits > leading_hits)
		 //     {
		 //             leading_hits = n_pfp_hits;
		 //             leading_index = j;
		 //     }
		 // }//end loop pfparticles
		 // auto const leading_shower = tpc_obj.GetParticle(leading_index);
		 // const double leading_theta = acos(leading_shower.pfpDirZ()) * 180 / 3.1415;
		 //
		 // if(leading_theta > 60)
		 // {
		 //     passed_tpco->at(i).first = 0;
		 //     passed_tpco->at(i).second = "Contained";
		 // }
	}//end loop tpc objects
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
//***************************************************************************
