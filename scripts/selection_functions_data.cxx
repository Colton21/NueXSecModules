#include "selection_functions_data.h"

//***************************************************************************
//***************************************************************************
void selection_functions_data::FillPostCutVectorData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                     std::vector<std::pair<int, std::string> > * passed_tpco,
                                                     std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v)
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
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double opening_angle = leading_shower.pfpOpenAngle();

		std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> my_tuple =
		        std::make_tuple(event_num, run_num, pfp_vtx_x, pfp_vtx_y, pfp_vtx_z, reason, tpco_id, num_tracks, num_showers, opening_angle);
		post_cuts_v->push_back(my_tuple);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::TabulateOriginsData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco,
                                                   std::vector<int> * tabulated_origins_data)
{
	int num_data = 0;
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		num_data++;
	}
	tabulated_origins_data->at(0) = num_data;
}
//***************************************************************************
//***************************************************************************
//modify this so it takes a string of the cut name so I only pass it a few variable at a time,
//then I can call this function several times later at the bottom
void selection_functions_data::PrintInfoData(int counter, std::string cut_name)
{
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
		const int pfp_pdg = part.PFParticlePdgCode();
		if(pfp_pdg == 11)
		{
			if(n_pfp_hits > most_hits)
			{
				leading_index = j;
				most_hits = n_pfp_hits;
			}
		}
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
		h_dedx_cuts_data->Fill(leading_dedx * (242.72 / 196.979));
	}        //end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutsdEdxAltScaleData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, const double alt_scale,
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
		h_dedx_cuts_data->Fill(leading_dedx * (242.72 / 196.979) * alt_scale);
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
				//const double trk_length = part.pfpLength();
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
void selection_functions_data::PostCutOpenAngle1ShowerData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                           TH1D * h_leading_shower_open_angle_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }

		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		if(n_pfp_showers != 1) {continue; }
		const int n_pfp = tpc_obj.NumPFParticles();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		h_leading_shower_open_angle_data->Fill(leading_open_angle);
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutOpenAngle2PlusShowerData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                               TH1D * h_leading_shower_open_angle_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }

		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		if(n_pfp_showers  < 2) {continue; }
		const int n_pfp = tpc_obj.NumPFParticles();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		h_leading_shower_open_angle_data->Fill(leading_open_angle);
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingPhiData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                              TH1D * h_ele_pfp_phi_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_phi = atan2(leading_shower.pfpDirY(), leading_shower.pfpDirX()) * 180 / 3.1415;
		h_ele_pfp_phi_data->Fill(leading_shower_phi);
	}//end pfp loop
}
//leading shower phi
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingThetaData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                TH1D * h_ele_pfp_theta_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		//const double leading_shower_theta = acos(leading_shower.pfpDirZ()) * (180 / 3.1415);

		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, 0, 0);
		const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);


		h_ele_pfp_theta_data->Fill(leading_shower_theta);
	}//end pfp loop
}
//leading shower theta
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingMomentumData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, TH1D * h_ele_pfp_momentum_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_momentum = leading_shower.pfpMomentum();
		h_ele_pfp_momentum_data->Fill(leading_shower_momentum);
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingMomentumThetaSliceData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                             TH1D * h_ele_momentum_1_data,
                                                             TH1D * h_ele_momentum_2_data,
                                                             TH1D * h_ele_momentum_3_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_momentum = leading_shower.pfpMomentum();
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, 0, 0);
		const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);
		if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_data->Fill(leading_shower_momentum); }
		if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_data->Fill(leading_shower_momentum); }
		if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_data->Fill(leading_shower_momentum); }
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::dedxThetaSliceData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                  TH1D * h_dedx_1_data,
                                                  TH1D * h_dedx_2_data,
                                                  TH1D * h_dedx_3_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_dedx = leading_shower.PfpdEdx().at(2);
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, 0, 0);
		const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);
		if(leading_shower_theta >= 0 && leading_shower_theta < 60)    {h_dedx_1_data->Fill(leading_shower_dedx * (242.72 / 196.979)); }
		if(leading_shower_theta >= 60 && leading_shower_theta < 120)   {h_dedx_2_data->Fill(leading_shower_dedx * (242.72 / 196.979)); }
		if(leading_shower_theta >= 120 && leading_shower_theta <= 180) {h_dedx_3_data->Fill(leading_shower_dedx * (242.72 / 196.979)); }
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingMomentumTrackTopologyData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                                std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                                TH1D * h_ele_pfp_momentum_no_track_data,
                                                                TH1D * h_ele_pfp_momentum_has_track_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_momentum = leading_shower.pfpMomentum();
		if(n_pfp_tracks == 0)
		{
			h_ele_pfp_momentum_no_track_data->Fill(leading_shower_momentum);
		}
		if(n_pfp_tracks >= 1)
		{
			h_ele_pfp_momentum_has_track_data->Fill(leading_shower_momentum);
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingThetaTrackTopologyData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                             TH1D * h_ele_pfp_theta_no_track_data,
                                                             TH1D * h_ele_pfp_theta_has_track_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		//const double leading_shower_momentum = leading_shower.pfpMomentum();
		const double leading_shower_theta = acos(leading_shower.pfpDirZ()) * 180 / 3.1415;
		//const double leading_shower_phi = atan2(leading_shower.pfpDirY(), leading_shower.pfpDirX()) * 180 / 3.1415;
		if(n_pfp_tracks == 0)
		{
			h_ele_pfp_theta_no_track_data->Fill(leading_shower_theta);
		}
		if(n_pfp_tracks >= 1)
		{
			h_ele_pfp_theta_has_track_data->Fill(leading_shower_theta);
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingPhiTrackTopologyData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                           TH1D * h_ele_pfp_phi_no_track_data,
                                                           TH1D * h_ele_pfp_phi_has_track_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		//const double leading_shower_momentum = leading_shower.pfpMomentum();
		//const double leading_shower_theta = acos(leading_shower.pfpDirZ()) * 180 / 3.1415;
		const double leading_shower_phi = atan2(leading_shower.pfpDirY(), leading_shower.pfpDirX()) * 180 / 3.1415;
		if(n_pfp_tracks == 0)
		{
			h_ele_pfp_phi_no_track_data->Fill(leading_shower_phi);
		}
		if(n_pfp_tracks >= 1)
		{
			h_ele_pfp_phi_has_track_data->Fill(leading_shower_phi);
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutsVtxFlashData(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, TH1D * h_vtx_flash_data)
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
		const int n_pfp = tpc_obj.NumPFParticles();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;

		//let's plot the smallest distance in the case of multiple showers
		//that way we can see if the cut works properly - i.e. if there's a
		//smallest distance greater than cut value, can see it's not working
		double min_distance = 10000;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const pfp = tpc_obj.GetParticle(j);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg != 11) {continue; }
			const double pfp_vtx_x = pfp.pfpVtxX();
			const double pfp_vtx_y = pfp.pfpVtxY();
			const double pfp_vtx_z = pfp.pfpVtxZ();

			const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) +
			                             pow((tpc_vtx_y - pfp_vtx_y), 2) +
			                             pow((tpc_vtx_z - pfp_vtx_z), 2));
			if(distance < min_distance) {min_distance = distance; }
		}
		if(tpco_id == "unmatched")
		{
			h_shwr_vtx_dist_data->Fill(min_distance);
		}
	}        //end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::TopologyEfficiencyData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                      std::vector<int> * no_track, std::vector<int> * has_track,
                                                      std::vector<int> * _1_shwr, std::vector<int> * _2_shwr,
                                                      std::vector<int> * _3_shwr, std::vector<int> * _4_shwr)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		//signal
		if(n_pfp_showers == 1)
		{
			_1_shwr->at(1) += 1;
		}
		if(n_pfp_showers == 2)
		{
			_2_shwr->at(1) += 1;
		}
		if(n_pfp_showers == 3)
		{
			_3_shwr->at(1) += 1;
		}
		if(n_pfp_showers >= 4)
		{
			_4_shwr->at(1) += 1;
		}
		if(n_pfp_tracks == 0)
		{
			no_track->at(1) += 1;
		}
		if(n_pfp_tracks >= 1)
		{
			has_track->at(1) += 1;
		}
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
			h_dedx_open_angle_data->Fill(leading_dedx * (242.72 / 196.979), leading_open_angle);
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
		if(n_pfp_showers <= 1) {continue; }
		double max_distance = 0;
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
			if(pfp_pdg == 11 && distance > max_distance)        //22 cm is ~ 2 radiation lengths
			{max_distance = distance; }
		}
		if(tpco_id == "unmatched")
		{
			h_second_shwr_dist_data->Fill(max_distance);
		}
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingCosThetaData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco,
                                                   const double theta_translation, const double phi_translation, bool _verbose, TH1D * h_ele_cos_theta_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, theta_translation, phi_translation);
		double leading_shower_cos_theta = shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag());
		if(theta_translation == 0 && phi_translation == 0) {leading_shower_cos_theta = leading_shower_z; }
		h_ele_cos_theta_data->Fill(leading_shower_cos_theta);
	}//end pfp loop
}
//leading shower cos theta - before most selection cuts
//***************************************************************************
//***************************************************************************
void selection_functions_data::HitsPlots1DData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                               TH1D * h_collection_hits_track_data,
                                               TH1D * h_collection_hits_shower_data,
                                               TH1D * h_collection_hits_leading_shower_data,
                                               TH1D * h_total_hits_leading_shower_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int n_pfp_hits_w_leading_shower = leading_shower.NumPFPHitsW();
		const int n_pfp_hits_leading_shower = leading_shower.NumPFPHits();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		std::vector<int> n_pfp_hits_w_track_v;
		std::vector<int> n_pfp_hits_w_shower_v;
		for(int i = 0; i < n_pfp; i++)
		{
			auto const pfp = tpc_obj.GetParticle(i);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg == 13)
			{
				const int n_pfp_hits_w_track = pfp.NumPFPHitsW();
				n_pfp_hits_w_track_v.push_back(n_pfp_hits_w_track);
			}
			if(pfp_pdg == 11)
			{
				const int n_pfp_hits_w_shower = pfp.NumPFPHitsW();
				n_pfp_hits_w_shower_v.push_back(n_pfp_hits_w_shower);
			}
		}

		if(n_pfp_tracks != 0) {for(auto const n_pfp_hits_w_track : n_pfp_hits_w_track_v) {h_collection_hits_track_data->Fill(n_pfp_hits_w_track); }}
		for(auto const n_pfp_hits_w_shower : n_pfp_hits_w_shower_v) {h_collection_hits_shower_data->Fill(n_pfp_hits_w_shower); }
		h_collection_hits_leading_shower_data->Fill(n_pfp_hits_w_leading_shower);
		h_total_hits_leading_shower_data->Fill(n_pfp_hits_leading_shower);
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::NumShowersOpenAngleData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                       std::vector<std::pair<int, std::string> > * passed_tpco, TH1D * h_pfp_shower_open_angle_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		h_pfp_shower_open_angle_data->Fill(n_pfp_showers);
	}        //end loop tpc objects

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
void selection_functions_data::PlaneHitsComparisonTrackData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                            std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                            TH2D * h_collection_total_hits_track_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		if(n_pfp_tracks == 0) {continue; }
		int n_pfp_hits_w = 0;
		int n_pfp_hits = 0;
		for(int i = 0; i < n_pfp; i++)
		{
			auto const pfp = tpc_obj.GetParticle(i);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg == 13)
			{
				n_pfp_hits_w = pfp.NumPFPHitsW();
				n_pfp_hits = pfp.NumPFPHits();
			}
		}
		h_collection_total_hits_track_data->Fill(n_pfp_hits_w, n_pfp_hits);
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::PlaneHitsComparisonShowerData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                             TH2D * h_collection_total_hits_shower_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		int n_pfp_hits_w = 0;
		int n_pfp_hits = 0;
		for(int i = 0; i < n_pfp; i++)
		{
			auto const pfp = tpc_obj.GetParticle(i);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg == 11)
			{
				n_pfp_hits_w = pfp.NumPFPHitsW();
				n_pfp_hits = pfp.NumPFPHits();
			}
		}
		h_collection_total_hits_shower_data->Fill(n_pfp_hits_w, n_pfp_hits);
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::PlaneHitsComparisonLeadingShowerData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                                    TH2D * h_collection_total_hits_shower_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int n_pfp_hits_w = leading_shower.NumPFPHitsW();
		const int n_pfp_hits = leading_shower.NumPFPHits();
		h_collection_total_hits_shower_data->Fill(n_pfp_hits_w, n_pfp_hits);
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
void selection_functions_data::XYZPositionData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                               std::vector<TH1 * > * h_ele_pfp_xyz_data,
                                               TH2 * h_pfp_zy_vtx_data, double & xyz_near, double & xyz_far)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();

		if(tpc_vtx_x < 30 && tpc_vtx_y < -100 && tpc_vtx_z < 50) {xyz_near++; }
		if(tpc_vtx_x > 220 && tpc_vtx_y > 100 && tpc_vtx_z > 986) {xyz_far++; }

		h_ele_pfp_xyz_data->at(0)->Fill(tpc_vtx_x);
		h_ele_pfp_xyz_data->at(1)->Fill(tpc_vtx_y);
		h_ele_pfp_xyz_data->at(2)->Fill(tpc_vtx_z);
		h_pfp_zy_vtx_data->Fill(tpc_vtx_z, tpc_vtx_y);
	}
}
void selection_functions_data::XYZPositionData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                               std::vector<TH1 * > * h_ele_pfp_xyz_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		h_ele_pfp_xyz_data->at(0)->Fill(tpc_vtx_x);
		h_ele_pfp_xyz_data->at(1)->Fill(tpc_vtx_y);
		h_ele_pfp_xyz_data->at(2)->Fill(tpc_vtx_z);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::EnergyCosThetaData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                  TH2 * h_ele_eng_costheta_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double momentum = leading_shower.pfpMomentum();
		const double costheta = leading_shower.pfpDirZ();
		h_ele_eng_costheta_data->Fill(momentum, costheta);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::EnergyCosThetaSlicesData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                        const double theta_translation, const double phi_translation,
                                                        TH1 * h_ele_eng_for_data,
                                                        TH1 * h_ele_eng_mid_data,
                                                        TH1 * h_ele_eng_back_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double momentum = leading_shower.pfpMomentum();
		double costheta = leading_shower.pfpDirZ();
		double theta = (acos(costheta) * (180 / 3.1415)) + theta_translation;
		costheta = cos(theta);
		if(costheta >= 0.5) {h_ele_eng_for_data->Fill(momentum); }
		if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_data->Fill(momentum); }
		if(costheta <= -0.5) {h_ele_eng_back_data->Fill(momentum); }
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingShowerLengthData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                       TH1D * h_shwr_length_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shwr_length = leading_shower.pfpLength();
		h_shwr_length_data->Fill(leading_shwr_length);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingShowerTrackLengthsData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                             TH1D * h_shwr_trk_length_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		std::pair<std::string, int> tpco_class = TPCO_Classifier_Data(tpc_obj);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shwr_length = leading_shower.pfpLength();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		if(n_pfp_tracks == 0) {continue; }
		double longest_track = 0;
		for(int i = 0; i < n_pfp; i++)
		{
			auto const pfp = tpc_obj.GetParticle(i);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			const double trk_length = pfp.pfpLength();
			if(pfp_pdg == 13 && trk_length > longest_track)
			{
				longest_track = trk_length;
			}
		}
		const double longest_trk_leading_shwr_ratio = longest_track / leading_shwr_length;
		h_shwr_trk_length_data->Fill(longest_trk_leading_shwr_ratio);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutsVtxFlashUpstreamData(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                            std::vector<std::pair<int, std::string> > * passed_tpco, TH1D * h_vtx_flash_data)
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
		if(flash_vtx_z < tpc_vtx_z) {continue; }//skip this event if the flash is downstream
		const int n_pfp = tpc_obj.NumPFParticles();
		const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
		h_vtx_flash_data->Fill(distance);
	}        //end loop tpc objects
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutsVtxFlashDownstreamData(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                              std::vector<std::pair<int, std::string> > * passed_tpco, TH1D * h_vtx_flash_data)
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
		if(flash_vtx_z >= tpc_vtx_z) {continue; }//skip this event if the flash is upstream
		const int n_pfp = tpc_obj.NumPFParticles();
		const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
		h_vtx_flash_data->Fill(distance);
	}        //end loop tpc objects
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions_data::PostCutsLeadingMomentumData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                           TH1D * h_ele_momentum_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int most_hits = 0;
		int leading_index = 0;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int n_pfp_hits = part.NumPFPHits();
			const int pfp_pdg = part.PFParticlePdgCode();
			if(pfp_pdg == 11)
			{
				if(n_pfp_hits > most_hits)
				{
					leading_index = j;
					most_hits = n_pfp_hits;
				}
			}
		}
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_momentum = leading_shower.pfpMomentum();
		h_ele_momentum_data->Fill(leading_momentum);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::dEdxThetaData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                             TH2D * h_dedx_theta_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int most_hits = 0;
		int leading_index = 0;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int n_pfp_hits = part.NumPFPHits();
			const int pfp_pdg = part.PFParticlePdgCode();
			if(pfp_pdg == 11)
			{
				if(n_pfp_hits > most_hits)
				{
					leading_index = j;
					most_hits = n_pfp_hits;
				}
			}
		}
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_shower_theta = acos(leading_shower.pfpDirZ()) * 180 / 3.1415;
		h_dedx_theta_data->Fill(leading_dedx * (242.72 / 196.979), leading_shower_theta);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::EventMultiplicityData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                     TH1D * h_multiplicity_shower_data, TH1D * h_multiplicity_track_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{

		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		h_multiplicity_shower_data->Fill(n_pfp_showers);
		h_multiplicity_track_data->Fill(n_pfp_tracks);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::LeadingKinematicsShowerTopologyData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                                   TH1D * h_ele_pfp_momentum_1shwr_data,
                                                                   TH1D * h_ele_pfp_momentum_2shwr_data,
                                                                   TH1D * h_ele_pfp_theta_1shwr_data,
                                                                   TH1D * h_ele_pfp_theta_2shwr_data,
                                                                   TH1D * h_ele_pfp_phi_1shwr_data,
                                                                   TH1D * h_ele_pfp_phi_2shwr_data)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int most_hits = 0;
		int leading_index = 0;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int n_pfp_hits = part.NumPFPHits();
			const int pfp_pdg = part.PFParticlePdgCode();
			if(pfp_pdg == 11)
			{
				if(n_pfp_hits > most_hits)
				{
					leading_index = j;
					most_hits = n_pfp_hits;
				}
			}
		}
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_momentum = leading_shower.pfpMomentum();
		const double leading_theta = acos(leading_shower.pfpDirZ()) * 180 / 3.1415;
		const double leading_phi = atan2(leading_shower.pfpDirY(), leading_shower.pfpDirX()) * 180 / 3.1415;
		const int num_showers = tpc_obj.NPfpShowers();
		if(num_showers == 1)
		{
			h_ele_pfp_momentum_1shwr_data->Fill(leading_momentum);
			h_ele_pfp_theta_1shwr_data->Fill(leading_theta);
			h_ele_pfp_phi_1shwr_data->Fill(leading_phi);
		}
		if(num_showers >= 2)
		{
			h_ele_pfp_momentum_2shwr_data->Fill(leading_momentum);
			h_ele_pfp_theta_2shwr_data->Fill(leading_theta);
			h_ele_pfp_phi_2shwr_data->Fill(leading_phi);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::IsContainedPlotData(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                   std::vector<double> fv_boundary_v, TH1 * h_track_containment_data)
{
	int n_tpc_obj = tpc_object_container_v->size();

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

	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }

		auto const tpc_obj = tpc_object_container_v->at(i);

		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		if(n_pfp_tracks == 0) {continue; }

		int is_contained;

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

				if(pfp_vtx_x <= det_x1 + x1 || pfp_vtx_x >= det_x2 - x2) {is_contained = 0; break; }
				if(pfp_vtx_y <= det_y1 + y1 || pfp_vtx_y >= det_y2 - y2) {is_contained = 0; break; }
				if(pfp_vtx_z <= det_z1 + z1 || pfp_vtx_z >= det_z2 - z2) {is_contained = 0; break; }
				if(pfp_end_x <= det_x1 + x1 || pfp_end_x >= det_x2 - x2) {is_contained = 0; break; }
				if(pfp_end_y <= det_y1 + y1 || pfp_end_y >= det_y2 - y2) {is_contained = 0; break; }
				if(pfp_end_z <= det_z1 + z1 || pfp_end_z >= det_z2 - z2) {is_contained = 0; break; }
				is_contained = 1;
			}
		}
		h_track_containment_data->Fill(is_contained);
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions_data::EvaluatedEdxMethodData(
        std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
        std::vector<std::pair<int, std::string> > * passed_tpco,
        TH1D * h_dedx,
        TH1D * h_dedx_cali,
        TH1D * h_dedx_omit,
        TH1D * h_dedx_omit_cali,
        TH2D * dedx_yz_ratio_cali,
        TH2D * dedx_yz_ratio_omit,
        TH2D * dedx_yz_ratio_omit_cali
        )
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		int most_hits = 0;
		int leading_index = 0;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int n_pfp_hits = part.NumPFPHits();
			const int pfp_pdg = part.PFParticlePdgCode();
			if(pfp_pdg == 11)
			{
				if(n_pfp_hits > most_hits)
				{
					leading_index = j;
					most_hits = n_pfp_hits;
				}
			}
		}
		auto const leading_shower = tpc_obj.GetParticle(leading_index);

		//multiple methods of getting the dE/dx
		const double dedx = leading_shower.PfpdEdx().at(2);
		const double dedx_cali = leading_shower.PfpdEdx_cali().at(2);
		const double dedx_omit = leading_shower.PfpdEdxOmitFirst().at(2);
		const double dedx_omit_cali = leading_shower.PfpdEdxOmitFirst_cali().at(2);

		//const double leading_shower_x = leading_shower.pfpVtxX();
		const double leading_shower_y = leading_shower.pfpVtxY();
		const double leading_shower_z = leading_shower.pfpVtxZ();
		const double ratio_cali = (dedx / dedx_cali);
		const double ratio_omit = (dedx / dedx_omit);
		const double ratio_omit_cali = (dedx_omit / dedx_omit_cali);

		//fill these for all MC
		dedx_yz_ratio_cali->Fill(leading_shower_z, leading_shower_y, ratio_cali);
		dedx_yz_ratio_omit->Fill(leading_shower_z, leading_shower_y, ratio_omit);
		dedx_yz_ratio_omit_cali->Fill(leading_shower_z, leading_shower_y, ratio_omit_cali);

		h_dedx->Fill(dedx);
		h_dedx_cali->Fill(dedx_cali);
		h_dedx_omit->Fill(dedx_omit);
		h_dedx_omit_cali->Fill(dedx_omit_cali);
	} //end pfp loop
}
//***************************************************************************
//***************************************************************************


//end functions
