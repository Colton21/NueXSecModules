#include "selection_functions_data.h"

//***************************************************************************
//***************************************************************************
void selection_functions_data::FillPostCutVectorData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                     std::vector<std::pair<int, std::string> > * passed_tpco,
                                                     std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int> > * post_cuts_v)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		const std::string reason = passed_tpco->at(i).second; //this should always return passed
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double pfp_vtx_x = tpc_obj.pfpVtxX();
		const double pfp_vtx_y = tpc_obj.pfpVtxY();
		const double pfp_vtx_z = tpc_obj.pfpVtxZ();
		const int run_num = tpc_obj.RunNumber();
		const int event_num = tpc_obj.EventNumber();
		const int num_tracks = tpc_obj.NPfpTracks();
		const int num_showers = tpc_obj.NPfpShowers();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;

		std::tuple<int, int, double, double, double, std::string, std::string, int, int> my_tuple =
		        std::make_tuple(event_num, run_num, pfp_vtx_x, pfp_vtx_y, pfp_vtx_z, reason, tpco_id, num_tracks, num_showers);
		post_cuts_v->push_back(my_tuple);
	}
}
//***************************************************************************
//***************************************************************************
std::vector<int> selection_functions_data::TabulateOriginsData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                               std::vector<std::pair<int, std::string> > * passed_tpco)
{
	int unmatched     = 0;
	int total         = 0;
	int signal_tpco_num = -1;
	std::vector<int> tabulated_origins;
	tabulated_origins.resize(22, 0);

	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		//loop over pfparticles in the TPCO
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
		if(tpco_id == "unmatched")     {unmatched++; }
	}//end loop tpc objects

	total = unmatched;

	tabulated_origins.at(5)  = unmatched;
	tabulated_origins.at(7)  = total;
	return tabulated_origins;
}
//***************************************************************************
//***************************************************************************
//modify this so it takes a string of the cut name so I only pass it a few variable at a time,
//then I can call this function several times later at the bottom
void selection_functions_data::PrintInfoData(std::vector<int> * counter_v, std::string cut_name)
{
	int counter = counter_v->at(7);
	std::cout << " <" << cut_name << "> " << std::endl;
	std::cout << " Total Candidate Nue     : " << counter << std::endl;
	std::cout << "------------------------" << std::endl;
	std::cout << "------------------------" << std::endl;
}
//***************************************************************************
//***************************************************************************
std::pair<std::string, int> selection_functions_data::TPCO_Classifier_Data(xsecAna::TPCObjectContainer tpc_obj)
{
	const int n_pfp = tpc_obj.NumPFParticles();
	int most_hits = 0;
	int leading_index = 0;
	for(int j = 0; j < n_pfp; j++)
	{
		auto const part = tpc_obj.GetParticle(j);
		const int n_pfp_hits = part.NumPFPHits();
		if(n_pfp_hits > most_hits) {leading_index = j; most_hits = n_pfp_hits; }
	}
	return std::make_pair("unmatched",     leading_index);
	//return the string for the tpco id
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutsdEdxData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                TH1D * h_dedx_cuts_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int most_hits = 0;
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		if(tpco_id == "unmatched")
		{
			h_dedx_cuts_data->Fill(leading_dedx);
		}
	}        //end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutOpenAngleData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                    TH1D * h_leading_shower_open_angle_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		if(tpco_id == "unmatched")
		{
			h_leading_shower_open_angle_data->Fill(leading_open_angle);
		}
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutTrkVtxData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                 TH1D * h_trk_vtx_dist_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		int * leading_index;
		int most_hits = 0;
		int num_tracks = 0;
		int num_showers = 0;

		double smallest_trk_vtx_dist = 1000;
		bool has_track = false;

		auto const tpc_obj = tpc_object_container_v->at(i);
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
		}//end pfp loop
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		if(tpco_id == "unmatched")
		{
			if(has_track) {h_trk_vtx_dist_data->Fill(smallest_trk_vtx_dist); }
		}
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::TopologyPlots1Data(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                  std::vector<std::pair<int, std::string> > * passed_tpco,
                                                  TH2D * h_pfp_track_shower_data,
                                                  TH1D * h_pfp_track_data,
                                                  TH1D * h_pfp_shower_data
                                                  )
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_tracks  = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		if(tpco_id == "unmatched")
		{
			h_pfp_track_shower_data->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_data->Fill(n_pfp_tracks);
			h_pfp_shower_data->Fill(n_pfp_showers);
		}
	}        //end loop tpc objects

}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::TopologyPlots2Data(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                  std::vector<std::pair<int, std::string> > * passed_tpco,
                                                  TH2D * h_pfp_track_shower_data,
                                                  TH1D * h_pfp_track_data,
                                                  TH1D * h_pfp_shower_data
                                                  )
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_tracks  = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		if(tpco_id == "unmatched")
		{
			h_pfp_track_shower_data->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_data->Fill(n_pfp_tracks);
			h_pfp_shower_data->Fill(n_pfp_showers);
		}
	}        //end loop tpc objects

}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutsVtxFlashData(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                    TH1D * h_vtx_flash_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const double flash_vtx_y = largest_flash_v.at(0);
		const double flash_vtx_z = largest_flash_v.at(1);
		const int n_pfp = tpc_obj.NumPFParticles();
		const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		if(tpco_id == "unmatched")
		{
			h_vtx_flash_data->Fill(distance);
		}
	}        //end loop tpc objects
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutsShwrVtxData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                   TH1D * h_shwr_vtx_dist_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_vtx_x = leading_shower.pfpVtxX();
		const double leading_vtx_y = leading_shower.pfpVtxY();
		const double leading_vtx_z = leading_shower.pfpVtxZ();
		const double distance = sqrt(pow((tpc_vtx_x - leading_vtx_x), 2) +
		                             pow((tpc_vtx_y - leading_vtx_y), 2) +
		                             pow((tpc_vtx_z - leading_vtx_z), 2));
		if(tpco_id == "unmatched")
		{
			h_shwr_vtx_dist_data->Fill(distance);
		}
	}        //end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::TopologyEfficiencyData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                      std::vector<int> * no_track_data, std::vector<int> * has_track_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();

		if(n_pfp_showers == 0)
		{
			std::cout << "Event with No Shower! WTF! " << std::endl;
			continue;
		}
		//signal
		if(n_pfp_tracks == 0) {no_track_data->at(0) += 1; }
		if(n_pfp_tracks >= 1) {has_track_data->at(0) += 1; }
	}//end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::dEdxVsOpenAngleData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                   TH2D * h_dedx_open_angle_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		if(tpco_id == "unmatched")
		{
			h_dedx_open_angle_data->Fill(leading_dedx, leading_open_angle);
		}
	}
}
//***************************************************************************
//***************************************************************************
//shower hits vs shower length, to see if we can better use the hit threshold cut to remove less signal?
void selection_functions_data::ShowerLengthvsHitsData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                      TH2D * h_shwr_len_hits_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int leading_hits = leading_shower.NumPFPHits();
		const double leading_length = leading_shower.pfpLength();
		if(tpco_id == "unmatched")
		{
			h_shwr_len_hits_data->Fill(leading_length, leading_hits);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::SecondaryShowersDistData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                        TH1D * h_second_shwr_dist_data)
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		int leading_index = tpco_class.second;
		//auto const leading_shower = tpc_obj.GetParticle(leading_index);
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
					if(tpco_id == "unmatched")
					{
						h_second_shwr_dist_data->Fill(distance);
					}
				}//end if reco shower
			}//end loop pfp
		}//end if tpco > 3 reco showers
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::HitLengthRatioData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                  TH1D * h_hit_length_ratio_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int pfp_pdg = leading_shower.PFParticlePdgCode();
		const double pfp_hits = leading_shower.NumPFPHits();
		const double pfp_length = leading_shower.pfpLength();
		const double pfp_hits_length_ratio = (pfp_hits / pfp_length);
		if(pfp_pdg == 11)
		{
			if(tpco_id == "unmatched")
			{
				h_hit_length_ratio_data->Fill(pfp_hits_length_ratio);
			}
		}                //end if reco shower
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
int selection_functions_data::MapFailureCutToStringData(const std::string failure_cut)
{
	int failure_reason = 10.0;
	if(failure_cut == "Passed")         {return failure_reason;   }
	if(failure_cut == "HasNue")         {failure_reason = 0.0;    }
	if(failure_cut == "InFV")           {failure_reason = 1.0;    }
	if(failure_cut == "FlashDist")      {failure_reason = 2.0;    }
	if(failure_cut == "ShwrVtx")        {failure_reason = 3.0;    }
	if(failure_cut == "TrkVtx")         {failure_reason = 4.0;    }
	if(failure_cut == "HitThreshold")   {failure_reason = 5.0;    }
	if(failure_cut == "OpenAngle")      {failure_reason = 6.0;    }
	if(failure_cut == "dEdX")           {failure_reason = 7.0;    }
	if(failure_cut == "SecondaryDist")  {failure_reason = 8.0;    }
	if(failure_cut == "HitLengthRatio") {failure_reason = 9.0;    }
	return failure_reason;
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::FailureReasonData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                 std::vector<std::pair<int, std::string> > * passed_tpco,
                                                 TH1D * h_failure_reason_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		//if(passed_tpco->at(i).first == 0) {continue; }
		const std::string failure_cut = passed_tpco->at(i).second;
		const int failure_reason = MapFailureCutToStringData(failure_cut);
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		if(tpco_id == "unmatched")
		{
			h_failure_reason_data->Fill(failure_reason);
		}
	}//end loop tpco
}
//***************************************************************************
//***************************************************************************

//end functions
