#include "selection_functions.h"
#include "selection_cuts.h"

//***************************************************************************
void selection_functions::PostCutsdEdx(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                       std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                       TH1D * h_dedx_cuts_nue_cc,
                                       TH1D * h_dedx_cuts_nue_cc_mixed,
                                       TH1D * h_dedx_cuts_nue_cc_out_fv,
                                       TH1D * h_dedx_cuts_numu_cc,
                                       TH1D * h_dedx_cuts_nc,
                                       TH1D * h_dedx_cuts_cosmic,
                                       TH1D * h_dedx_cuts_nc_pi0,
                                       TH1D * h_dedx_cuts_numu_cc_mixed,
                                       TH1D * h_dedx_cuts_other_mixed,
                                       TH1D * h_dedx_cuts_unmatched     )
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{

		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_dedx_cuts_nue_cc->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_dedx_cuts_nue_cc_out_fv->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_dedx_cuts_nue_cc->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_dedx_cuts_nue_cc->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_dedx_cuts_nue_cc->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_dedx_cuts_nue_cc->Fill(leading_dedx);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_dedx_cuts_numu_cc->Fill(leading_dedx);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_dedx_cuts_numu_cc->Fill(leading_dedx);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_dedx_cuts_numu_cc->Fill(leading_dedx);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_dedx_cuts_numu_cc->Fill(leading_dedx);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_dedx_cuts_numu_cc->Fill(leading_dedx);
		}
		if(tpco_id == "nc")
		{
			h_dedx_cuts_nc->Fill(leading_dedx);
		}
		if(tpco_id == "nc_pi0")
		{
			h_dedx_cuts_nc_pi0->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_dedx_cuts_nue_cc_mixed->Fill(leading_dedx);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_dedx_cuts_numu_cc_mixed->Fill(leading_dedx);
			h_dedx_cuts_numu_cc->Fill(leading_dedx);
		}
		if(tpco_id == "cosmic")
		{
			h_dedx_cuts_cosmic->Fill(leading_dedx);
		}
		if(tpco_id == "other_mixed")
		{
			h_dedx_cuts_other_mixed->Fill(leading_dedx);
		}
		if(tpco_id == "unmatched")
		{
			h_dedx_cuts_unmatched->Fill(leading_dedx);
		}
	}        //end loop tpc objects
}
//***************************************************************************
void selection_functions::PostCutsdEdxInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                             TH1D * h_dedx_cuts_intime)
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
		int most_hits = 0;
		int leading_index = 0;
		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			const int pfp_pdg = part.PFParticlePdgCode();
			if(pfp_pdg != 11) {continue; }
			const int n_pfp_hits = part.NumPFPHits();
			if(n_pfp_hits > most_hits) {leading_index = j; most_hits = n_pfp_hits; }
		}
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		h_dedx_cuts_intime->Fill(leading_dedx * (242.72 / 196.979));
	}        //end loop tpc objects
}
//***************************************************************************
void selection_functions::PostCutsdEdxTrueParticle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                   std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                   TH1D * h_dedx_cuts_electron,
                                                   TH1D * h_dedx_cuts_photon,
                                                   TH1D * h_dedx_cuts_proton,
                                                   TH1D * h_dedx_cuts_pion,
                                                   TH1D * h_dedx_cuts_muon,
                                                   TH1D * h_dedx_cuts_kaon,
                                                   TH1D * h_dedx_cuts_neutron,
                                                   TH1D * h_dedx_cuts_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_mc_pdg = leading_shower.MCPdgCode();
		bool good_id = false;
		if(leading_mc_pdg == 11 || leading_mc_pdg == -11)
		{
			h_dedx_cuts_electron->Fill(leading_dedx);
			good_id = true;
		}
		if(leading_mc_pdg == 13 || leading_mc_pdg == -13)
		{
			h_dedx_cuts_muon->Fill(leading_dedx);
			good_id = true;
		}
		if(leading_mc_pdg == 22)
		{
			h_dedx_cuts_photon->Fill(leading_dedx);
			good_id = true;
		}
		if(leading_mc_pdg == 2212)
		{
			h_dedx_cuts_proton->Fill(leading_dedx);
			good_id = true;
		}
		if(leading_mc_pdg == 211 || leading_mc_pdg == -211)
		{
			h_dedx_cuts_pion->Fill(leading_dedx);
			good_id = true;
		}
		if(leading_mc_pdg == 2112)
		{
			h_dedx_cuts_neutron->Fill(leading_dedx);
			good_id = true;
		}
		if(leading_mc_pdg == 130 || leading_mc_pdg == 310 || leading_mc_pdg == 311 || leading_mc_pdg == 321 || leading_mc_pdg == -321)
		{
			h_dedx_cuts_kaon->Fill(leading_dedx);
			good_id = true;
		}
		if(good_id == false )
		{
			h_dedx_cuts_unmatched->Fill(leading_dedx);
		}
	}        //end loop tpc objects
}
void selection_functions::PostCutsdEdxTrueParticleInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                         TH1 * h_dedx_cuts_ext_unmatched)
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
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		h_dedx_cuts_ext_unmatched->Fill(leading_dedx * (242.72 / 196.979));
	}//end loop tpc objects
}
//***************************************************************************
void selection_functions::PostCutsdEdxHitsTrueParticle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                       std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                       TH2D * h_dedx_cuts_electron,
                                                       TH2D * h_dedx_cuts_photon,
                                                       TH2D * h_dedx_cuts_proton,
                                                       TH2D * h_dedx_cuts_pion,
                                                       TH2D * h_dedx_cuts_muon,
                                                       TH2D * h_dedx_cuts_kaon,
                                                       TH2D * h_dedx_cuts_neutron,
                                                       TH2D * h_dedx_cuts_unmatched,
                                                       TH2D * h_dedx_cuts_collection_electron,
                                                       TH2D * h_dedx_cuts_collection_photon,
                                                       TH2D * h_dedx_cuts_collection_proton,
                                                       TH2D * h_dedx_cuts_collection_pion,
                                                       TH2D * h_dedx_cuts_collection_muon,
                                                       TH2D * h_dedx_cuts_collection_kaon,
                                                       TH2D * h_dedx_cuts_collection_neutron,
                                                       TH2D * h_dedx_cuts_collection_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_mc_pdg = leading_shower.MCPdgCode();
		const int leading_shower_hits = leading_shower.NumPFPHits();
		const int leading_shower_collection_hits = leading_shower.NumPFPHitsW();
		if(leading_mc_pdg == 11 || leading_mc_pdg == -11)
		{
			h_dedx_cuts_electron->Fill(leading_dedx, leading_shower_hits);
			h_dedx_cuts_collection_electron->Fill(leading_dedx, leading_shower_collection_hits);
		}
		if(leading_mc_pdg == 13 || leading_mc_pdg == -13)
		{
			h_dedx_cuts_muon->Fill(leading_dedx, leading_shower_hits);
			h_dedx_cuts_collection_muon->Fill(leading_dedx, leading_shower_collection_hits);
		}
		if(leading_mc_pdg == 22)
		{
			h_dedx_cuts_photon->Fill(leading_dedx, leading_shower_hits);
			h_dedx_cuts_collection_photon->Fill(leading_dedx, leading_shower_collection_hits);
		}
		if(leading_mc_pdg == 2212)
		{
			h_dedx_cuts_proton->Fill(leading_dedx, leading_shower_hits);
			h_dedx_cuts_collection_proton->Fill(leading_dedx, leading_shower_collection_hits);
		}
		if(leading_mc_pdg == 211 || leading_mc_pdg == -211)
		{
			h_dedx_cuts_pion->Fill(leading_dedx, leading_shower_hits);
			h_dedx_cuts_collection_pion->Fill(leading_dedx, leading_shower_collection_hits);
		}
		if(leading_mc_pdg == 2112)
		{
			h_dedx_cuts_neutron->Fill(leading_dedx, leading_shower_hits);
			h_dedx_cuts_collection_neutron->Fill(leading_dedx, leading_shower_collection_hits);
		}
		if(leading_mc_pdg == 130 || leading_mc_pdg == 310 || leading_mc_pdg == 311 || leading_mc_pdg == 321 || leading_mc_pdg == -321)
		{
			h_dedx_cuts_kaon->Fill(leading_dedx, leading_shower_hits);
			h_dedx_cuts_collection_kaon->Fill(leading_dedx, leading_shower_collection_hits);
		}
		if(leading_mc_pdg == 0)
		{
			h_dedx_cuts_unmatched->Fill(leading_dedx, leading_shower_hits);
			h_dedx_cuts_collection_unmatched->Fill(leading_dedx, leading_shower_collection_hits);
		}
	}        //end loop tpc objects
}
//***************************************************************************
void selection_functions::FillPostCutVector(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                            std::vector<std::pair<int, std::string> > * passed_tpco,
                                            std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                            std::vector<std::tuple<int, int, int, double, double, double,
                                                                   std::string, std::string, int, int, double> > * post_cuts_v)
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
		const int sub_run_num = tpc_obj.SubRunNumber();
		const int event_num = tpc_obj.EventNumber();
		const int num_tracks = tpc_obj.NPfpTracks();
		const int num_showers = tpc_obj.NPfpShowers();

		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();

		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double opening_angle = leading_shower.pfpOpenAngle();

		std::tuple<int, int, int, double, double, double, std::string, std::string, int, int, double> my_tuple =
		        std::make_tuple(event_num, run_num, sub_run_num, pfp_vtx_x, pfp_vtx_y, pfp_vtx_z, reason, tpco_id, num_tracks, num_showers, opening_angle);
		post_cuts_v->push_back(my_tuple);
	}
}
void selection_functions::FillPostCutVector(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                            std::vector<std::pair<int, std::string> > * passed_tpco,
                                            std::vector<std::tuple<int, int, int, double, double, double,
                                                                   std::string, std::string, int, int, double> > * post_cuts_v)
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
		const int sub_run_num = tpc_obj.SubRunNumber();
		const int event_num = tpc_obj.EventNumber();
		const int num_tracks = tpc_obj.NPfpTracks();
		const int num_showers = tpc_obj.NPfpShowers();

		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		std::string tpco_id = "Data";
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
		const double opening_angle = leading_shower.pfpOpenAngle();

		std::tuple<int, int, int, double, double, double, std::string, std::string, int, int, double> my_tuple =
		        std::make_tuple(event_num, run_num, sub_run_num, pfp_vtx_x, pfp_vtx_y, pfp_vtx_z, reason, tpco_id, num_tracks, num_showers, opening_angle);
		post_cuts_v->push_back(my_tuple);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::PrintPostCutVector(std::vector<std::tuple<int, int, int, double, double, double,
                                                                    std::string, std::string, int, int, double> > * post_cuts_v, bool _post_cuts_verbose)
{
	const int passed_events = post_cuts_v->size();
	std::cout << "* * * * * * * * * * * * * * * * *" << std::endl;
	std::cout << "Total Passed Events: " << passed_events << std::endl;
	std::cout << "* * * * * * * * * * * * * * * * *" << std::endl;
	for(auto const my_tuple: * post_cuts_v)
	{
		const int event_num = std::get<0>(my_tuple);
		const int run_num = std::get<1>(my_tuple);
		const int sub_run_num = std::get<2>(my_tuple);
		const double pfp_vtx_x = std::get<3>(my_tuple);
		const double pfp_vtx_y = std::get<4>(my_tuple);
		const double pfp_vtx_z = std::get<5>(my_tuple);
		const std::string reason = std::get<6>(my_tuple);
		const std::string event_type = std::get<7>(my_tuple);
		const int num_tracks = std::get<8>(my_tuple);
		const int num_showers = std::get<9>(my_tuple);
		const double opening_angle = std::get<10>(my_tuple);
		std::cout << "* * * * * * * * * * * * * * * * *" << std::endl;
		std::cout << "Event Type     : " << event_type << std::endl;
		std::cout << "Event Number   : " << event_num << std::endl;
		std::cout << "Run Number     : " << run_num << std::endl;
		std::cout << "SubRun Number  : " << sub_run_num << std::endl;
		std::cout << "Num PFP Tracks : " << num_tracks << std::endl;
		std::cout << "Num PFP Showers: " << num_showers << std::endl;
		std::cout << "Pfp Vtx X      : " << pfp_vtx_x << std::endl;
		std::cout << "Pfp Vtx Y      : " << pfp_vtx_y << std::endl;
		std::cout << "Pfp Vtx Z      : " << pfp_vtx_z << std::endl;
		std::cout << "Opening Angle  : " << opening_angle * (180 / 3.1415) << std::endl;
		std::cout << "TPCO Reason    : " << reason << std::endl;
		std::cout << "* * * * * * * * * * * * * * * * *" << std::endl;
	}
	std::cout << "   * * * END * * *   " << std::endl;
	std::cout << "* * * * * * * * * * * * * * * * *" << std::endl;
}
void selection_functions::PostCutVectorPlots(std::vector<std::tuple<int, int, int, double, double, double, std::string,
                                                                    std::string, int, int, double> > * post_cuts_v,
                                             bool _post_cuts_verbose, TH1 * post_cuts_num_showers_purity_qe,
                                             TH1 * post_cuts_num_showers_purity_res,
                                             TH1 * post_cuts_num_showers_purity_dis,
                                             TH1 * post_cuts_num_showers_purity_coh,
                                             TH1 * post_cuts_num_showers_purity_mec)
{
	int signal_events_1_qe = 0;
	int signal_events_2_qe = 0;
	int signal_events_3_qe = 0;
	int signal_events_4_qe = 0;

	int signal_events_1_res = 0;
	int signal_events_2_res = 0;
	int signal_events_3_res = 0;
	int signal_events_4_res = 0;

	int signal_events_1_dis = 0;
	int signal_events_2_dis = 0;
	int signal_events_3_dis = 0;
	int signal_events_4_dis = 0;

	int signal_events_1_coh = 0;
	int signal_events_2_coh = 0;
	int signal_events_3_coh = 0;
	int signal_events_4_coh = 0;

	int signal_events_1_mec = 0;
	int signal_events_2_mec = 0;
	int signal_events_3_mec = 0;
	int signal_events_4_mec = 0;

	int bkg_events_1 = 0;
	int bkg_events_2 = 0;
	int bkg_events_3 = 0;
	int bkg_events_4 = 0;
	//this loops through all events which passed the selection cuts
	for(auto const my_tuple : * post_cuts_v)
	{
		const int event_num = std::get<0>(my_tuple);
		const int run_num = std::get<1>(my_tuple);
		const int sub_run_num = std::get<2>(my_tuple);
		const double pfp_vtx_x = std::get<3>(my_tuple);
		const double pfp_vtx_y = std::get<4>(my_tuple);
		const double pfp_vtx_z = std::get<5>(my_tuple);
		const std::string reason = std::get<6>(my_tuple);
		const std::string event_type = std::get<7>(my_tuple);
		const int num_tracks = std::get<8>(my_tuple);
		const int num_showers = std::get<9>(my_tuple);
		const double opening_angle = std::get<10>(my_tuple);

		int signal_bkg = 0;
		if(event_type == "nue_cc_qe")   {signal_bkg = 1; }
		if(event_type == "nue_cc_res")  {signal_bkg = 2; }
		if(event_type == "nue_cc_dis")  {signal_bkg = 3; }
		if(event_type == "nue_cc_coh")  {signal_bkg = 4; }
		if(event_type == "nue_cc_mec")  {signal_bkg = 5; }
		if(signal_bkg == 1)
		{
			if(num_showers == 1) {signal_events_1_qe++; }
			if(num_showers == 2) {signal_events_2_qe++; }
			if(num_showers == 3) {signal_events_3_qe++; }
			if(num_showers >= 4) {signal_events_4_qe++; }
		}
		if(signal_bkg == 2)
		{
			if(num_showers == 1) {signal_events_1_res++; }
			if(num_showers == 2) {signal_events_2_res++; }
			if(num_showers == 3) {signal_events_3_res++; }
			if(num_showers >= 4) {signal_events_4_res++; }
		}
		if(signal_bkg == 3)
		{
			if(num_showers == 1) {signal_events_1_dis++; }
			if(num_showers == 2) {signal_events_2_dis++; }
			if(num_showers == 3) {signal_events_3_dis++; }
			if(num_showers >= 4) {signal_events_4_dis++; }
		}
		if(signal_bkg == 4)
		{
			if(num_showers == 1) {signal_events_1_coh++; }
			if(num_showers == 2) {signal_events_2_coh++; }
			if(num_showers == 3) {signal_events_3_coh++; }
			if(num_showers >= 4) {signal_events_4_coh++; }
		}
		if(signal_bkg == 5)
		{
			if(num_showers == 1) {signal_events_1_mec++; }
			if(num_showers == 2) {signal_events_2_mec++; }
			if(num_showers == 3) {signal_events_3_mec++; }
			if(num_showers >= 4) {signal_events_4_mec++; }
		}
		if(signal_bkg == 0)
		{
			if(num_showers == 1) {bkg_events_1++; }
			if(num_showers == 2) {bkg_events_2++; }
			if(num_showers == 3) {bkg_events_3++; }
			if(num_showers >= 4) {bkg_events_4++; }
		}
	}
	const double total_events_1 = signal_events_1_qe + signal_events_1_res + signal_events_1_dis + signal_events_1_coh + signal_events_1_mec + bkg_events_1;
	const double total_events_2 = signal_events_2_qe + signal_events_2_res + signal_events_2_dis + signal_events_2_coh + signal_events_2_mec + bkg_events_2;
	const double total_events_3 = signal_events_3_qe + signal_events_3_res + signal_events_3_dis + signal_events_3_coh + signal_events_3_mec + bkg_events_3;
	const double total_events_4 = signal_events_4_qe + signal_events_4_res + signal_events_4_dis + signal_events_4_coh + signal_events_4_mec + bkg_events_4;

	//std::cout << "4+ Shower Background Events: " << bkg_events_4 << std::endl;

	double purity_1_qe  = double(signal_events_1_qe)  / total_events_1;
	double purity_2_qe  = double(signal_events_2_qe)  / total_events_2;
	double purity_3_qe  = double(signal_events_3_qe)  / total_events_3;
	double purity_4_qe  = double(signal_events_4_qe)  / total_events_4;

	double purity_1_res = double(signal_events_1_res) / total_events_1;
	double purity_2_res = double(signal_events_2_res) / total_events_2;
	double purity_3_res = double(signal_events_3_res) / total_events_3;
	double purity_4_res = double(signal_events_4_res) / total_events_4;

	double purity_1_dis = double(signal_events_1_dis) / total_events_1;
	double purity_2_dis = double(signal_events_2_dis) / total_events_2;
	double purity_3_dis = double(signal_events_3_dis) / total_events_3;
	double purity_4_dis = double(signal_events_4_dis) / total_events_4;

	double purity_1_coh = double(signal_events_1_coh) / total_events_1;
	double purity_2_coh = double(signal_events_2_coh) / total_events_2;
	double purity_3_coh = double(signal_events_3_coh) / total_events_3;
	double purity_4_coh = double(signal_events_4_coh) / total_events_4;

	double purity_1_mec = double(signal_events_1_mec) / total_events_1;
	double purity_2_mec = double(signal_events_2_mec) / total_events_2;
	double purity_3_mec = double(signal_events_3_mec) / total_events_3;
	double purity_4_mec = double(signal_events_4_mec) / total_events_4;

	post_cuts_num_showers_purity_qe->SetBinContent(1, purity_1_qe);
	post_cuts_num_showers_purity_qe->SetBinContent(2, purity_2_qe);
	post_cuts_num_showers_purity_qe->SetBinContent(3, purity_3_qe);
	post_cuts_num_showers_purity_qe->SetBinContent(4, purity_4_qe);

	post_cuts_num_showers_purity_res->SetBinContent(1, purity_1_res);
	post_cuts_num_showers_purity_res->SetBinContent(2, purity_2_res);
	post_cuts_num_showers_purity_res->SetBinContent(3, purity_3_res);
	post_cuts_num_showers_purity_res->SetBinContent(4, purity_4_res);

	post_cuts_num_showers_purity_dis->SetBinContent(1, purity_1_dis);
	post_cuts_num_showers_purity_dis->SetBinContent(2, purity_2_dis);
	post_cuts_num_showers_purity_dis->SetBinContent(3, purity_3_dis);
	post_cuts_num_showers_purity_dis->SetBinContent(4, purity_4_dis);

	post_cuts_num_showers_purity_coh->SetBinContent(1, purity_1_coh);
	post_cuts_num_showers_purity_coh->SetBinContent(2, purity_2_coh);
	post_cuts_num_showers_purity_coh->SetBinContent(3, purity_3_coh);
	post_cuts_num_showers_purity_coh->SetBinContent(4, purity_4_coh);

	post_cuts_num_showers_purity_mec->SetBinContent(1, purity_1_mec);
	post_cuts_num_showers_purity_mec->SetBinContent(2, purity_2_mec);
	post_cuts_num_showers_purity_mec->SetBinContent(3, purity_3_mec);
	post_cuts_num_showers_purity_mec->SetBinContent(4, purity_4_mec);
}
//***************************************************************************
//this function just counts if at least 1 tpc object passes the cuts
bool selection_functions::ValidTPCObjects(std::vector<std::pair<int, std::string> > * passed_tpco)
{
	const int n_tpco = passed_tpco->size();
	int pass_sum = 0;
	for(auto const tpco : * passed_tpco)
	{
		pass_sum = pass_sum + tpco.first;
	}
	//0 --> failed tpco, 1 --> passed tpco
	if(pass_sum == 0) {return false; }
	return true;
}
//***************************************************************************
//***************************************************************************
void selection_functions::TabulateOriginsInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                std::vector<std::pair<int, std::string> > * passed_tpco,
                                                std::vector<int> * tabulated_origins_intime)
{
	int num_intime_cosmics = 0;
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		num_intime_cosmics++;
	}
	tabulated_origins_intime->at(0) = num_intime_cosmics;
}
//***************************************************************************
//***************************************************************************
void selection_functions::TabulateOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco,
                                          std::vector<int> * tabulated_origins,
                                          std::vector<std::pair<std::string, int> > * tpco_classifier_v)
{
	int nue_cc        = 0;

	int nue_cc_qe     = 0;
	int nue_cc_res    = 0;
	int nue_cc_dis    = 0;
	int nue_cc_coh    = 0;
	int nue_cc_mec    = 0;

	int nue_bar_cc_qe     = 0;
	int nue_bar_cc_res    = 0;
	int nue_bar_cc_dis    = 0;
	int nue_bar_cc_coh    = 0;
	int nue_bar_cc_mec    = 0;

	int total_nue_cc_qe     = 0;
	int total_nue_cc_res    = 0;
	int total_nue_cc_dis    = 0;
	int total_nue_cc_coh    = 0;
	int total_nue_cc_mec    = 0;

	int nue_cc_mixed  = 0;
	int nue_cc_out_fv = 0;
	int cosmic        = 0;
	int nc            = 0;
	int numu_cc       = 0;
	int numu_cc_qe    = 0;
	int numu_cc_res   = 0;
	int numu_cc_dis   = 0;
	int numu_cc_coh   = 0;
	int numu_cc_mec   = 0;
	int numu_cc_mixed = 0;
	int nc_pi0        = 0;
	int unmatched     = 0;
	int other_mixed   = 0;
	int total         = 0;
	int signal_tpco_num = -1;
	int only_nue_cc = 0;
	int only_nue_bar_cc = 0;

	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		//auto const tpc_obj = tpc_object_container_v->at(i);
		//int leading_index = tpco_class.second;
		std::string tpco_id = tpco_classifier_v->at(i).first;

		if(tpco_id == "nue_cc_qe")       {nue_cc_qe++;  }
		if(tpco_id == "nue_cc_res")      {nue_cc_res++; }
		if(tpco_id == "nue_cc_coh")      {nue_cc_coh++; }
		if(tpco_id == "nue_cc_dis")      {nue_cc_dis++; }
		if(tpco_id == "nue_cc_mec")      {nue_cc_mec++; }

		if(tpco_id == "nue_bar_cc_qe")   {nue_bar_cc_qe++;  }
		if(tpco_id == "nue_bar_cc_res")  {nue_bar_cc_res++; }
		if(tpco_id == "nue_bar_cc_coh")  {nue_bar_cc_coh++; }
		if(tpco_id == "nue_bar_cc_dis")  {nue_bar_cc_dis++; }
		if(tpco_id == "nue_bar_cc_mec")  {nue_bar_cc_mec++; }

		if(tpco_id == "nue_cc_out_fv")   {nue_cc_out_fv++; }
		if(tpco_id == "nue_cc_mixed")    {nue_cc_mixed++; }
		if(tpco_id == "nc")              {nc++; }
		if(tpco_id == "numu_cc_qe")      {numu_cc_qe++; }
		if(tpco_id == "numu_cc_res")     {numu_cc_res++; }
		if(tpco_id == "numu_cc_coh")     {numu_cc_coh++; }
		if(tpco_id == "numu_cc_dis")     {numu_cc_dis++; }
		if(tpco_id == "numu_cc_mec")     {numu_cc_mec++; }
		if(tpco_id == "numu_cc_mixed")   {numu_cc_mixed++; }
		if(tpco_id == "nc_pi0")          {nc_pi0++; }
		if(tpco_id == "cosmic")          {cosmic++; }
		if(tpco_id == "other_mixed")     {other_mixed++; }
		if(tpco_id == "unmatched")       {unmatched++; }
	}//end loop tpc objects

	only_nue_cc = nue_cc_qe + nue_cc_res + nue_cc_dis + nue_cc_coh + nue_cc_mec;
	only_nue_bar_cc = nue_bar_cc_qe + nue_bar_cc_res + nue_bar_cc_dis + nue_bar_cc_coh + nue_bar_cc_mec;

	total_nue_cc_qe  = nue_cc_qe  + nue_bar_cc_qe;
	total_nue_cc_res = nue_cc_res + nue_bar_cc_res;
	total_nue_cc_dis = nue_cc_dis + nue_bar_cc_dis;
	total_nue_cc_coh = nue_cc_coh + nue_bar_cc_coh;
	total_nue_cc_mec = nue_cc_mec + nue_bar_cc_mec;


	nue_cc = only_nue_cc + only_nue_bar_cc;
	numu_cc = numu_cc_qe + numu_cc_res + numu_cc_dis + numu_cc_coh + numu_cc_mec;
	total = nue_cc + nue_cc_mixed + nue_cc_out_fv + cosmic + nc + numu_cc + numu_cc_mixed + nc_pi0 + unmatched + other_mixed;

	tabulated_origins->at(0)  = nue_cc;//this is nue_cc + nue_bar_cc
	tabulated_origins->at(1)  = nue_cc_mixed;
	tabulated_origins->at(2)  = cosmic;
	tabulated_origins->at(3)  = nc;
	tabulated_origins->at(4)  = numu_cc;
	tabulated_origins->at(5)  = unmatched;
	tabulated_origins->at(6)  = other_mixed;
	tabulated_origins->at(7)  = total;
	tabulated_origins->at(8)  = signal_tpco_num;
	tabulated_origins->at(9)  = nue_cc_out_fv;
	tabulated_origins->at(10) = nc_pi0;
	tabulated_origins->at(11) = numu_cc_mixed;
	tabulated_origins->at(12) = total_nue_cc_qe;
	tabulated_origins->at(13) = total_nue_cc_res;
	tabulated_origins->at(14) = total_nue_cc_dis;
	tabulated_origins->at(15) = total_nue_cc_coh;
	tabulated_origins->at(16) = total_nue_cc_mec;
	tabulated_origins->at(17) = numu_cc_qe;
	tabulated_origins->at(18) = numu_cc_res;
	tabulated_origins->at(19) = numu_cc_dis;
	tabulated_origins->at(20) = numu_cc_coh;
	tabulated_origins->at(21) = numu_cc_mec;
	tabulated_origins->at(22) = only_nue_cc;
	tabulated_origins->at(23) = only_nue_bar_cc;
}
//***************************************************************************
//***************************************************************************
void selection_functions::TotalOrigins(std::vector<int> * tabulated_origins, std::vector<int> * total_cut_origins)
{
	total_cut_origins->at(0)  += tabulated_origins->at(0);
	total_cut_origins->at(1)  += tabulated_origins->at(1);
	total_cut_origins->at(9)  += tabulated_origins->at(9);
	total_cut_origins->at(2)  += tabulated_origins->at(2);
	total_cut_origins->at(3)  += tabulated_origins->at(3);
	total_cut_origins->at(4)  += tabulated_origins->at(4);
	total_cut_origins->at(11) += tabulated_origins->at(11);
	total_cut_origins->at(10) += tabulated_origins->at(10);
	total_cut_origins->at(5)  += tabulated_origins->at(5);
	total_cut_origins->at(6)  += tabulated_origins->at(6);
	total_cut_origins->at(7)  += tabulated_origins->at(7);
	total_cut_origins->at(12) += tabulated_origins->at(12);
	total_cut_origins->at(13) += tabulated_origins->at(13);
	total_cut_origins->at(14) += tabulated_origins->at(14);
	total_cut_origins->at(15) += tabulated_origins->at(15);
	total_cut_origins->at(16) += tabulated_origins->at(16);
	total_cut_origins->at(17) += tabulated_origins->at(17);
	total_cut_origins->at(18) += tabulated_origins->at(18);
	total_cut_origins->at(19) += tabulated_origins->at(19);
	total_cut_origins->at(20) += tabulated_origins->at(20);
	total_cut_origins->at(21) += tabulated_origins->at(21);
	total_cut_origins->at(22) += tabulated_origins->at(22);
	total_cut_origins->at(23) += tabulated_origins->at(23);
}
void selection_functions::TotalOriginsInTime(std::vector<int> * tabulated_origins, std::vector<int> * total_cut_origins)
{
	total_cut_origins->at(0)  += tabulated_origins->at(0);
}
//***************************************************************************
//***************************************************************************
//modify this so it takes a string of the cut name so I only pass it a few variable at a time,
//then I can call this function several times later at the bottom
void selection_functions::PrintInfo(int mc_nue_cc_counter, std::vector<int> * counter_v, int counter_intime_cosmics,
                                    double intime_scale_factor, double data_scale_factor, std::string cut_name)
{
	int counter                = counter_v->at(7);
	int counter_nue_cc         = counter_v->at(0);
	int counter_nue_cc_mixed   = counter_v->at(1);
	int counter_nue_cc_out_fv  = counter_v->at(9);
	int counter_cosmic         = counter_v->at(2);
	int counter_nc             = counter_v->at(3);
	int counter_numu_cc        = counter_v->at(4);
	int counter_numu_cc_mixed  = counter_v->at(11);
	int counter_nc_pi0         = counter_v->at(10);
	int counter_unmatched      = counter_v->at(5);
	int counter_other_mixed    = counter_v->at(6);
	int counter_nue_cc_qe      = counter_v->at(12);
	int counter_nue_cc_res     = counter_v->at(13);
	int counter_nue_cc_dis     = counter_v->at(14);
	int counter_nue_cc_coh     = counter_v->at(15);
	int counter_nue_cc_mec     = counter_v->at(16);
	int counter_numu_cc_qe     = counter_v->at(17);
	int counter_numu_cc_res    = counter_v->at(18);
	int counter_numu_cc_dis    = counter_v->at(19);
	int counter_numu_cc_coh    = counter_v->at(20);
	int counter_numu_cc_mec    = counter_v->at(21);

	counter = counter + (counter_intime_cosmics * (intime_scale_factor / data_scale_factor));

	std::cout << " <" << cut_name << "> " << std::endl;
	std::cout << " Total Candidate Nue     : " << counter                << "\t \t " << double(counter                * data_scale_factor  ) << std::endl;
	std::cout << " Number of Nue CC        : " << counter_nue_cc         << "\t \t " << double(counter_nue_cc         * data_scale_factor  ) << std::endl;
	std::cout << " Number of Nue CC Mixed  : " << counter_nue_cc_mixed   << "\t \t " << double(counter_nue_cc_mixed   * data_scale_factor  ) << std::endl;
	std::cout << " Number of Nue CC out FV : " << counter_nue_cc_out_fv  << "\t \t " << double(counter_nue_cc_out_fv  * data_scale_factor  ) << std::endl;
	std::cout << " Number of Cosmic        : " << counter_cosmic         << "\t \t " << double(counter_cosmic         * data_scale_factor  ) << std::endl;
	std::cout << " Number of Numu CC       : " << counter_numu_cc        << "\t \t " << double(counter_numu_cc        * data_scale_factor  ) << std::endl;
	std::cout << " Number of Numu CC Mixed : " << counter_numu_cc_mixed  << "\t \t " << double(counter_numu_cc_mixed  * data_scale_factor  ) << std::endl;
	std::cout << " Number of NC            : " << counter_nc             << "\t \t " << double(counter_nc             * data_scale_factor  ) << std::endl;
	std::cout << " Number of NC Pi0        : " << counter_nc_pi0         << "\t \t " << double(counter_nc_pi0         * data_scale_factor  ) << std::endl;
	std::cout << " Number of Unmatched     : " << counter_unmatched      << "\t \t " << double(counter_unmatched      * data_scale_factor  ) << std::endl;
	std::cout << " Number of Other Mixed   : " << counter_other_mixed    << "\t \t " << double(counter_other_mixed    * data_scale_factor  ) << std::endl;
	std::cout << " Number of InTime Cosmics: " << double(counter_intime_cosmics * (intime_scale_factor / data_scale_factor))
	          << "\t \t " << double(counter_intime_cosmics * intime_scale_factor) << std::endl;
	std::cout << "---------Unscaled----------" << std::endl;
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
	std::cout << " Efficiency       : " << "( " << counter_nue_cc << " / " << mc_nue_cc_counter << " ) = " << efficiency << std::endl;
	std::cout << " Purity           : " << purity << std::endl;
	std::cout << "------------------------" << std::endl;
	std::cout << "------------------------" << std::endl;
}
//***************************************************************************
//***************************************************************************
void selection_functions::ExportEfficiencyPurity(int mc_nue_cc_counter, std::vector<int> * counter_v, int counter_intime_cosmics,
                                                 double intime_scale_factor, double data_scale_factor, std::string cut_name,
                                                 std::vector<std::tuple< double, double, std::string> > * results_v)
{
	int counter           = counter_v->at(7);
	int counter_nue_cc    = counter_v->at(0);
	counter = counter + double(counter_intime_cosmics * (intime_scale_factor / data_scale_factor));
	const double efficiency = double(counter_nue_cc) / double(mc_nue_cc_counter);
	const double purity     = double(counter_nue_cc) / double(counter);

	auto result = std::make_tuple(efficiency, purity, cut_name);
	results_v->push_back(result);

}
//***************************************************************************
//***************************************************************************
void selection_functions::PrintTopologyPurity(std::vector<int> * no_track, std::vector<int> * has_track,
                                              std::vector<int> * _1_shwr, std::vector<int> * _2_shwr, std::vector<int> * _3_shwr, std::vector<int> * _4_shwr)
{
	std::cout << "---------------------" << std::endl;
	std::cout << "No Track Signal: " << no_track->at(0) << std::endl;
	std::cout << "No Track Bkg   : " << no_track->at(1) << std::endl;
	std::cout << "Purity         : " << double(no_track->at(0)) / double(no_track->at(0) + no_track->at(1)) << std::endl;
	std::cout << " ******************* " << std::endl;
	std::cout << "1+ Track Signal: " << has_track->at(0) << std::endl;
	std::cout << "1+ Track Bkg   : " << has_track->at(1) << std::endl;
	std::cout << "Purity         : " << double(has_track->at(0)) / double(has_track->at(0) + has_track->at(1)) << std::endl;
	std::cout << "---------------------" << std::endl;
	std::cout << "---------------------" << std::endl;
	std::cout << "1 Shower Signal : " << _1_shwr->at(0) << std::endl;
	std::cout << "1 Shower Bkg    : " << _1_shwr->at(1) << std::endl;
	std::cout << "Purity          : " << double(_1_shwr->at(0)) / double(_1_shwr->at(0) + _1_shwr->at(1)) << std::endl;
	std::cout << " ******************* " << std::endl;
	std::cout << "2 Shower Signal : " << _2_shwr->at(0) << std::endl;
	std::cout << "2 Shower Bkg    : " << _2_shwr->at(1) << std::endl;
	std::cout << "Purity          : " << double(_2_shwr->at(0)) / double(_2_shwr->at(0) + _2_shwr->at(1)) << std::endl;
	std::cout << " ******************* " << std::endl;
	std::cout << "3 Shower Signal : " << _3_shwr->at(0) << std::endl;
	std::cout << "3 Shower Bkg    : " << _3_shwr->at(1) << std::endl;
	std::cout << "Purity          : " << double(_3_shwr->at(0)) / double(_3_shwr->at(0) + _3_shwr->at(1)) << std::endl;
	std::cout << " ******************* " << std::endl;
	std::cout << "4+ Shower Signal: " << _4_shwr->at(0) << std::endl;
	std::cout << "4+ Shower Bkg   : " << _4_shwr->at(1) << std::endl;
	std::cout << "Purity          : " << double(_4_shwr->at(0)) / double(_4_shwr->at(0) + _4_shwr->at(1)) << std::endl;
	std::cout << "---------------------" << std::endl;
	std::cout << "---------------------" << std::endl;
}
//***************************************************************************
//***************************************************************************
//this function performs a classification on a tpc object basis,
//by looping over each pfparticle contained
std::pair<std::string, int> selection_functions::TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool has_pi0, bool true_in_tpc)
{
	int part_nue_cc     = 0;
	int part_nue_bar_cc = 0;
	int part_cosmic     = 0;
	int part_nc         = 0;
	int part_nc_pi0     = 0;
	int part_numu_cc    = 0;
	int part_unmatched  = 0;

	const int tpc_obj_mode = tpc_obj.Mode();
	const int n_pfp = tpc_obj.NumPFParticles();
	const int n_pfp_showers = tpc_obj.NPfpShowers();
	int most_hits = 0;
	int leading_index = -1;
	int leading_pdg = 0;
	int leading_mc_parent_pdg = 0;
	std::string leading_origin = "kNothing";

	for(int j = 0; j < n_pfp; j++)
	{
		auto const part = tpc_obj.GetParticle(j);
		const int n_pfp_hits = part.NumPFPHits();
		const int mc_parent_pdg = part.MCParentPdg();
		const int pfp_pdg = part.PFParticlePdgCode();
		if(pfp_pdg == 11)
		{
			if(n_pfp_hits > most_hits)
			{
				leading_index = j;
				most_hits = n_pfp_hits;
			}
		}
		//if(n_pfp_showers)
		if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && mc_parent_pdg == 12)  { part_nue_cc++; }
		if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && mc_parent_pdg == -12) { part_nue_bar_cc++; }
		if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && (mc_parent_pdg == 14 || mc_parent_pdg == -14)) { part_numu_cc++; }
		if(part.CCNC() == 1 && part.Origin() == "kBeamNeutrino")
		{
			if(has_pi0 == true)  {part_nc_pi0++; }
			if(has_pi0 == false) {part_nc++; }
		}
		if(part.Origin() == "kCosmicRay") { part_cosmic++;    }
		if(part.Origin() == "kUnknown"  ) { part_unmatched++; }
	}
	//some tpc objects actually have 0 hits - crazy!
	if(tpc_obj.NumPFPHits() == 0) {return std::make_pair("bad_reco", 0); }

	//currently, any tpc objects which only have a track end up with a leading_index of -1
	//this index will likely cause code to crash if called before the signal definition cuts

	//also some rare cases where nu_pfp = nue, and shower hits = 0 with track hits > 0 - how does this happen? (NC event?)

	//now to catagorise the tpco
	if(part_cosmic > 0)
	{
		if(part_nue_cc  > 0 || part_nue_bar_cc > 0)  { return std::make_pair("nue_cc_mixed",  leading_index); }
		if(part_numu_cc > 0 )                        { return std::make_pair("numu_cc_mixed", leading_index); }
		if(part_nc  > 0 || part_nc_pi0 > 0)          { return std::make_pair("other_mixed",   leading_index); }
		return std::make_pair("cosmic", leading_index);
	}
	//this uses the true neutrino vertex for this specific event
	//not the true vtx per tpc object - maybe this can be fixed in the future...
	//but using the true nu vtx only matters for the pure signal events,
	//where the neutrino vertex IS the true tpc object vertex
	if(part_cosmic == 0)
	{
		if(part_nue_cc      > 0 && true_in_tpc == false) { return std::make_pair("nue_cc_out_fv", leading_index);   }
		if(part_nue_bar_cc  > 0 && true_in_tpc == false) { return std::make_pair("nue_cc_out_fv", leading_index);   }

		if(part_nue_cc    > 0 && tpc_obj_mode == 0   ) { return std::make_pair("nue_cc_qe",     leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 1   ) { return std::make_pair("nue_cc_res",    leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 2   ) { return std::make_pair("nue_cc_dis",    leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 3   ) { return std::make_pair("nue_cc_coh",    leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 10  ) { return std::make_pair("nue_cc_mec",    leading_index);   }

		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 0   ) { return std::make_pair("nue_bar_cc_qe",     leading_index);   }
		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 1   ) { return std::make_pair("nue_bar_cc_res",    leading_index);   }
		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 2   ) { return std::make_pair("nue_bar_cc_dis",    leading_index);   }
		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 3   ) { return std::make_pair("nue_bar_cc_coh",    leading_index);   }
		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 10  ) { return std::make_pair("nue_bar_cc_mec",    leading_index);   }

		if(part_numu_cc     > 0 && tpc_obj_mode == 0   ) { return std::make_pair("numu_cc_qe",    leading_index);   }
		if(part_numu_cc     > 0 && tpc_obj_mode == 1   ) { return std::make_pair("numu_cc_res",   leading_index);   }
		if(part_numu_cc     > 0 && tpc_obj_mode == 2   ) { return std::make_pair("numu_cc_dis",   leading_index);   }
		if(part_numu_cc     > 0 && tpc_obj_mode == 3   ) { return std::make_pair("numu_cc_coh",   leading_index);   }
		if(part_numu_cc     > 0 && tpc_obj_mode == 10  ) { return std::make_pair("numu_cc_mec",   leading_index);   }
		if(part_nc          > 0                        ) { return std::make_pair("nc",            leading_index);   }
		if(part_nc_pi0      > 0                        ) { return std::make_pair("nc_pi0",        leading_index);   }
		if(part_unmatched   > 0                        ) { return std::make_pair("unmatched",     leading_index);   }
	}
	//this never happens :)
	std::cout << "HELP HELP HELP END OF TPCO CLASSIFIER AND NO CLASSIFICATION!" << std::endl;
	//return the string for the tpco id
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::FillTPCOClassV(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v, bool true_in_tpc, bool has_pi0,
                                         std::vector<std::pair<std::string, int> > * tpco_classifier_v)
{
	for(auto const tpc_obj : * tpc_object_container_v)
	{
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, true_in_tpc);
		//int leading_index = tpco_class.second;
		//std::string tpco_id = tpco_class.first;
		tpco_classifier_v->push_back(tpco_class);
	}
}
//***************************************************************************
//***************************************************************************
double selection_functions::calcNumNucleons(const bool is_data, double _x1, double _x2, double _y1,
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

	double lar_density = 0;
	if(!is_data) {lar_density = 1.3954; } // g/cm^3
	if(is_data)  {lar_density = 1.3836; } // g/cm^3 -- comes from Joseph's note
	const double avagadro_num = 6.022e23; // molecule per mol
	const double n_nucleon = 40.0; // argon has 40
	const double mol_mass = 39.95; // g/mol
	//const double au = 1.67 * pow(10, -27); //kg
	//const double n_target = vol * lar_density / au;

	const double n_target = lar_density * vol * avagadro_num * n_nucleon / mol_mass;

	std::cout << "--------------------" << std::endl;
	std::cout << "  Calculating N_t   " << std::endl;
	std::cout << "--- vol  : " << vol << " cm3" << std::endl;
	std::cout << "--- p_lar: " << lar_density << " g / cm3" << std::endl;
	//std::cout << "--- AU   : " << au << " kg " << std::endl;
	std::cout << "---------------------" << std::endl;


	return n_target;
}
//***************************************************************************
//***************************************************************************
void selection_functions::calcXSec(const bool is_data, double _x1, double _x2, double _y1,
                                   double _y2, double _z1, double _z2,
                                   double n_total, double n_bkg, double flux, double efficiency, std::vector<double>  * xsec_cc)
{
	const int n_events = n_total - n_bkg;
	//scale_factor = 2.4 * math.pow(10, 17)  # POT / nue
	//calculate the number of nucleons based on the fiducial volume
	const double n_target = selection_functions::calcNumNucleons(is_data, _x1, _x2, _y1,
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
void selection_functions::XSecWork(double final_counter, double final_counter_nue_cc, double final_counter_nue_bar_cc, double final_counter_nue_cc_mixed,
                                   double final_counter_nue_cc_out_fv, double final_counter_cosmic, double final_counter_nc, double final_counter_numu_cc,
                                   double final_counter_numu_cc_mixed, double final_counter_nc_pi0, double final_counter_unmatched,
                                   double final_counter_other_mixed, double final_counter_intime,
                                   double intime_scale_factor, double final_counter_data, double data_scale_factor,
                                   std::vector<double> fv_boundary_v, double flux_nue, double flux_nue_bar,
                                   std::vector<double> selected_energy_vector, double genie_xsec_nue, double genie_xsec_nue_bar,
                                   const int total_mc_entries_inFV_nue, const int total_mc_entries_inFV_nue_bar)
{
	std::vector<double> * xsec_cc = new std::vector<double>;
	const int n_bkg_mc = (final_counter_nue_cc_mixed + final_counter_nue_cc_out_fv + final_counter_cosmic + final_counter_nc
	                      + final_counter_numu_cc + final_counter_numu_cc_mixed + final_counter_nc_pi0 +
	                      final_counter_unmatched + final_counter_other_mixed);
	int n_bkg_intime = final_counter_intime;
	const double n_bkg = n_bkg_mc + double(n_bkg_intime * (intime_scale_factor / data_scale_factor));
	const double n_bkg_data = double(n_bkg_mc * data_scale_factor) + double(n_bkg_intime * intime_scale_factor);

	const double flux_data = (flux_nue + flux_nue_bar) * data_scale_factor;
	const double total_mc_entries_inFV = (total_mc_entries_inFV_nue + total_mc_entries_inFV_nue_bar);

	const double efficiency_nue = final_counter_nue_cc / double(total_mc_entries_inFV_nue);
	const double efficiency_stat_err_nue = (1 / sqrt(total_mc_entries_inFV_nue)) * sqrt(efficiency_nue * (1 - efficiency_nue));

	const double efficiency_nue_bar = final_counter_nue_bar_cc / double(total_mc_entries_inFV_nue_bar);
	const double efficiency_stat_err_nue_bar = (1 / sqrt(total_mc_entries_inFV_nue_bar)) * sqrt(efficiency_nue_bar * (1 - efficiency_nue_bar));

	const double efficiency = (final_counter_nue_cc + final_counter_nue_bar_cc) / double(total_mc_entries_inFV);
	const double efficiency_stat_err = (1 / sqrt(total_mc_entries_inFV)) * sqrt(efficiency * (1 - efficiency));

	const double flux_nue_nue_bar = (flux_nue + flux_nue_bar);
	const double efficiency_nue_nue_bar = double(final_counter_nue_cc + final_counter_nue_bar_cc) /
	                                      double(total_mc_entries_inFV_nue + total_mc_entries_inFV_nue_bar);

	const double _x1 = fv_boundary_v.at(0);
	const double _x2 = fv_boundary_v.at(1);
	const double _y1 = fv_boundary_v.at(2);
	const double _y2 = fv_boundary_v.at(3);
	const double _z1 = fv_boundary_v.at(4);
	const double _z2 = fv_boundary_v.at(5);

	std::cout << "--- Cross Section Calculations ---" << std::endl;
	//******************************
	//******** Data ****************
	//******************************

	const int n_total_data = final_counter_data;
	selection_functions::calcXSec(true, _x1, _x2, _y1, _y2, _z1, _z2,
	                              n_total_data, n_bkg_data, flux_data,
	                              efficiency, xsec_cc);
	double xsec_cc_data = xsec_cc->at(0);
	double xsec_cc_stat_data = xsec_cc->at(0) * (pow((sqrt(n_total_data) / n_total_data), 2) + pow((efficiency_stat_err / efficiency), 2));

	std::cout << "-------------------------" << std::endl;
	std::cout << " Cross Section Results (Data):  " << std::endl;
	std::cout << " " << xsec_cc_data << " +/- (stats) "
	          << xsec_cc_stat_data << " +/- (sys) "
	          << xsec_cc->at(2) << std::endl;
	std::cout << "-------------------------" << std::endl;
	xsec_cc->clear();

	//************************************
	//******** Monte Carlo ***************
	//************************************
	//nue
	const int n_total_mc = final_counter_nue_cc + final_counter_nue_bar_cc + n_bkg;
	selection_functions::calcXSec(false, _x1, _x2, _y1, _y2, _z1, _z2,
	                              final_counter_nue_cc, 0, flux_nue,
	                              efficiency_nue, xsec_cc);
	double xsec_cc_stat_mc_nue = xsec_cc->at(0) * (pow((sqrt(n_total_mc) / n_total_mc), 2) +
	                                               pow((efficiency_stat_err_nue / efficiency_nue), 2));

	std::cout << "-------------------------" << std::endl;
	std::cout << " Cross Section Results (MC - nue):  " << std::endl;
	std::cout << " " << xsec_cc->at(0) << " +/- (stats) "
	          << xsec_cc_stat_mc_nue << " +/- (sys) "
	          << xsec_cc->at(2) << std::endl;
	std::cout << "-------------------------" << std::endl;
	xsec_cc->clear();

	//nue_bar
	selection_functions::calcXSec(false, _x1, _x2, _y1, _y2, _z1, _z2,
	                              final_counter_nue_bar_cc, 0, flux_nue_bar,
	                              efficiency_nue_bar, xsec_cc);
	double xsec_cc_stat_mc_nue_bar = xsec_cc->at(0) * (pow((sqrt(n_total_mc) / n_total_mc), 2) +
	                                                   pow((efficiency_stat_err_nue_bar / efficiency_nue_bar), 2));

	std::cout << "-------------------------" << std::endl;
	std::cout << " Cross Section Results (MC - nuebar):  " << std::endl;
	std::cout << " " << xsec_cc->at(0) << " +/- (stats) "
	          << xsec_cc_stat_mc_nue_bar << " +/- (sys) "
	          << xsec_cc->at(2) << std::endl;
	std::cout << "-------------------------" << std::endl;
	xsec_cc->clear();

	//nue+nue_bar
	selection_functions::calcXSec(false, _x1, _x2, _y1, _y2, _z1, _z2,
	                              (final_counter_nue_cc + final_counter_nue_bar_cc), 0, flux_nue_nue_bar,
	                              efficiency_nue_nue_bar, xsec_cc);
	double xsec_cc_nue_nue_bar_mc = xsec_cc->at(0);

	std::cout << "-------------------------" << std::endl;
	std::cout << " Cross Section Results (MC - nue + nuebar):  " << std::endl;
	std::cout << " " << xsec_cc->at(0) << " +/- (stats) "
	          << xsec_cc_stat_mc_nue_bar << " +/- (sys) "
	          << xsec_cc->at(2) << std::endl;
	std::cout << "-------------------------" << std::endl;
	xsec_cc->clear();

	//************************************
	//******** True Level ****************
	//************************************
	selection_functions::calcXSec(false, _x1, _x2, _y1, _y2, _z1, _z2,
	                              total_mc_entries_inFV_nue, 0, flux_nue,
	                              1, xsec_cc);
	double xsec_cc_stat_truth = xsec_cc->at(1);

	std::cout << "-------------------------" << std::endl;
	std::cout << " Cross Section Results (Truth - nue):  " << std::endl;
	std::cout << " " << xsec_cc->at(0) << " +/- (stats) "
	          << xsec_cc_stat_truth << " +/- (sys) "
	          << xsec_cc->at(2) << std::endl;
	xsec_cc->clear();

	selection_functions::calcXSec(false, _x1, _x2, _y1, _y2, _z1, _z2,
	                              total_mc_entries_inFV_nue_bar, 0, flux_nue_bar,
	                              1, xsec_cc);
	xsec_cc_stat_truth = xsec_cc->at(1);

	std::cout << "-------------------------" << std::endl;
	std::cout << " Cross Section Results (Truth - nuebar):  " << std::endl;
	std::cout << " " << xsec_cc->at(0) << " +/- (stats) "
	          << xsec_cc_stat_truth << " +/- (sys) "
	          << xsec_cc->at(2) << std::endl;


	std::cout << "-------------------------" << std::endl;
	std::cout << "-------------------------" << std::endl;
	std::cout << " Genie value of Flux " << '\n' <<
	        " Integrated Xsec (nue)   :    " << genie_xsec_nue << '\n' <<
	        " Integrated Xsec (nuebar):    " << genie_xsec_nue_bar << std::endl;
	xsec_cc->clear();

	double all_energy = 0;
	for(auto const energy : selected_energy_vector) {all_energy += energy; }
	const double average_true_energy = all_energy / selected_energy_vector.size();
	bool _verbose = false;
	selection_functions::xsec_plot(_verbose, genie_xsec_nue, genie_xsec_nue_bar, xsec_cc_data, average_true_energy, xsec_cc_stat_data);

	const double data_xsec = xsec_cc_data;
	const double mc_xsec = xsec_cc_nue_nue_bar_mc;
	const double stat_err = xsec_cc_stat_data;
	const double sys_err = data_xsec * sqrt( pow(0.20, 2) + pow(0.20, 2));

	IntegratedXsecPlot(data_xsec, mc_xsec, stat_err, sys_err);

}
//***************************************************************************
//***************************************************************************
void selection_functions::xsec_plot(bool _verbose, double genie_xsec_nue, double genie_xsec_nue_bar, double xsec_nue, double average_energy, double stat_error)
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
	double g_y[1] = {genie_xsec_nue};
	double g_ex[1] = {0.0};
	double g_ey[1] = {0.0};
	const int g_n = 1;
	TGraphErrors * g_genie_point = new TGraphErrors(g_n, g_x, g_y, g_ex, g_ey);
	g_genie_point->GetXaxis()->SetLimits(0.0, 5.0);//5 GeV
	g_genie_point->SetMinimum(0.0);
	g_genie_point->SetMaximum(50e-39);

	double x[1] = {average_energy};
	double y[1] = {xsec_nue};
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


	combined_c1->Print("../scripts/plots/combined_xsec.pdf");


	if(f->IsOpen()) {f->Close(); }
	if(f->IsOpen()) {xsec_f->Close(); }
}
void selection_functions::IntegratedXsecPlot(const double data_xsec, const double mc_xsec, const double stat_err, const double sys_err)
{
	//const double full_error = sqrt( pow(stat_err, 2) + pow(sys_err, 2) );
	const double full_error = sqrt( pow(0, 2) + pow(sys_err, 2) );

	TCanvas * c_integrated = new TCanvas();
	c_integrated->cd();

	double x[1] = {1.0};
	double y[1] = {data_xsec};
	double ex[1] = {0.0};
	double ey[1] = {full_error};
	const int n = 1;

	TGraphErrors * xsec_graph = new TGraphErrors(n, x, y, ex, ey);
	xsec_graph->GetXaxis()->SetLimits(0.9, 1.1);//5 GeV
	xsec_graph->SetMinimum(3.6e-39);
	xsec_graph->SetMaximum(7.2e-39);
	xsec_graph->SetMarkerStyle(3);
	xsec_graph->Draw();

	TLine * line = new TLine(0.9, mc_xsec, 1.1, mc_xsec);
	line->SetLineColor(46);
	line->Draw("same");

	TLegend * leg = new TLegend();
	leg->AddEntry(xsec_graph, "Data Cross Section", "l");
	leg->AddEntry(line, "GENIE", "l");
	leg->Draw();

	c_integrated->Print("../scripts/plots/integrated_xsec.pdf");
}

//***************************************************************************
//***************************************************************************
void selection_functions::PostCutOpenAngle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                           std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                           TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_out_fv,
                                           TH1D * h_leading_shower_open_angle_nue_cc_mixed,
                                           TH1D * h_leading_shower_open_angle_numu_cc, TH1D * h_leading_shower_open_angle_nc,
                                           TH1D * h_leading_shower_open_angle_cosmic, TH1D * h_leading_shower_open_angle_nc_pi0,
                                           TH1D * h_leading_shower_open_angle_numu_cc_mixed, TH1D * h_leading_shower_open_angle_other_mixed,
                                           TH1D * h_leading_shower_open_angle_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		int leading_index = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);

		if(tpco_id == "nue_cc_mixed")  {h_leading_shower_open_angle_nue_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "nue_cc_out_fv") {h_leading_shower_open_angle_nue_cc_out_fv->Fill(leading_open_angle); }
		//if(tpco_id == "numu_cc_mixed") {h_leading_shower_open_angle_numu_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "numu_cc_mixed") {h_leading_shower_open_angle_numu_cc->Fill(leading_open_angle); }
		if(tpco_id == "other_mixed")   {h_leading_shower_open_angle_other_mixed->Fill(leading_open_angle); }
		if(tpco_id == "cosmic")        {h_leading_shower_open_angle_cosmic->Fill(leading_open_angle); }
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_mec" ||
		   tpco_id == "nue_bar_cc_qe" || tpco_id == "nue_bar_cc_res" || tpco_id == "nue_bar_cc_dis" || tpco_id == "nue_bar_cc_coh" || tpco_id == "nue_bar_cc_mec")
		{
			h_leading_shower_open_angle_nue_cc->Fill(leading_open_angle);
		}
		if(tpco_id == "numu_cc_qe" || tpco_id == "numu_cc_res" || tpco_id == "numu_cc_dis" || tpco_id == "numu_cc_coh" || tpco_id == "numu_cc_mec")
		{
			h_leading_shower_open_angle_numu_cc->Fill(leading_open_angle);
		}
		if(tpco_id == "nc")        {h_leading_shower_open_angle_nc->Fill(leading_open_angle); }
		if(tpco_id == "nc_pi0")    {h_leading_shower_open_angle_nc_pi0->Fill(leading_open_angle); }
		if(tpco_id == "unmatched") {h_leading_shower_open_angle_unmatched->Fill(leading_open_angle); }
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutOpenAngleInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                 TH1D * h_leading_shower_open_angle_intime)
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
			const int pfp_pdg = part.PFParticlePdgCode();
			if(pfp_pdg != 11) {continue; }
			const int n_pfp_hits = part.NumPFPHits();
			if(n_pfp_hits > most_hits) {leading_index = j; most_hits = n_pfp_hits; }
		}
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		h_leading_shower_open_angle_intime->Fill(leading_open_angle);
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutOpenAngle1Shower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                  std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                  TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_out_fv,
                                                  TH1D * h_leading_shower_open_angle_nue_cc_mixed,
                                                  TH1D * h_leading_shower_open_angle_numu_cc, TH1D * h_leading_shower_open_angle_nc,
                                                  TH1D * h_leading_shower_open_angle_cosmic, TH1D * h_leading_shower_open_angle_nc_pi0,
                                                  TH1D * h_leading_shower_open_angle_numu_cc_mixed, TH1D * h_leading_shower_open_angle_other_mixed,
                                                  TH1D * h_leading_shower_open_angle_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }

		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		if(n_pfp_showers != 1) {continue; }
		int leading_index = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);

		if(tpco_id == "nue_cc_mixed")  {h_leading_shower_open_angle_nue_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "nue_cc_out_fv") {h_leading_shower_open_angle_nue_cc_out_fv->Fill(leading_open_angle); }
		//if(tpco_id == "numu_cc_mixed") {h_leading_shower_open_angle_numu_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "numu_cc_mixed") {h_leading_shower_open_angle_numu_cc->Fill(leading_open_angle); }
		if(tpco_id == "other_mixed")   {h_leading_shower_open_angle_other_mixed->Fill(leading_open_angle); }
		if(tpco_id == "cosmic")        {h_leading_shower_open_angle_cosmic->Fill(leading_open_angle); }
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_mec" ||
		   tpco_id == "nue_bar_cc_qe" || tpco_id == "nue_bar_cc_res" || tpco_id == "nue_bar_cc_dis" || tpco_id == "nue_bar_cc_coh" || tpco_id == "nue_bar_cc_mec")
		{
			h_leading_shower_open_angle_nue_cc->Fill(leading_open_angle);
		}
		if(tpco_id == "numu_cc_qe" || tpco_id == "numu_cc_res" || tpco_id == "numu_cc_dis" || tpco_id == "numu_cc_coh" || tpco_id == "numu_cc_mec")
		{
			h_leading_shower_open_angle_numu_cc->Fill(leading_open_angle);
		}
		if(tpco_id == "nc")        {h_leading_shower_open_angle_nc->Fill(leading_open_angle); }
		if(tpco_id == "nc_pi0")    {h_leading_shower_open_angle_nc_pi0->Fill(leading_open_angle); }
		if(tpco_id == "unmatched") {h_leading_shower_open_angle_unmatched->Fill(leading_open_angle); }
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutOpenAngle1ShowerInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                        TH1D * h_leading_shower_open_angle_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }

		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		if(n_pfp_showers != 1) {continue; }
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
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		h_leading_shower_open_angle_intime->Fill(leading_open_angle);
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutOpenAngle2PlusShower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                      std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                      TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_out_fv,
                                                      TH1D * h_leading_shower_open_angle_nue_cc_mixed,
                                                      TH1D * h_leading_shower_open_angle_numu_cc, TH1D * h_leading_shower_open_angle_nc,
                                                      TH1D * h_leading_shower_open_angle_cosmic, TH1D * h_leading_shower_open_angle_nc_pi0,
                                                      TH1D * h_leading_shower_open_angle_numu_cc_mixed, TH1D * h_leading_shower_open_angle_other_mixed,
                                                      TH1D * h_leading_shower_open_angle_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }

		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		if(n_pfp_showers  < 2) {continue; }
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);

		if(tpco_id == "nue_cc_mixed")  {h_leading_shower_open_angle_nue_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "nue_cc_out_fv") {h_leading_shower_open_angle_nue_cc_out_fv->Fill(leading_open_angle); }
		//if(tpco_id == "numu_cc_mixed") {h_leading_shower_open_angle_numu_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "numu_cc_mixed") {h_leading_shower_open_angle_numu_cc->Fill(leading_open_angle); }
		if(tpco_id == "other_mixed")   {h_leading_shower_open_angle_other_mixed->Fill(leading_open_angle); }
		if(tpco_id == "cosmic")        {h_leading_shower_open_angle_cosmic->Fill(leading_open_angle); }
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_mec" ||
		   tpco_id == "nue_bar_cc_qe" || tpco_id == "nue_bar_cc_res" || tpco_id == "nue_bar_cc_dis" || tpco_id == "nue_bar_cc_coh" || tpco_id == "nue_bar_cc_mec")
		{
			h_leading_shower_open_angle_nue_cc->Fill(leading_open_angle);
		}
		if(tpco_id == "numu_cc_qe" || tpco_id == "numu_cc_res" || tpco_id == "numu_cc_dis" || tpco_id == "numu_cc_coh" || tpco_id == "numu_cc_mec")
		{
			h_leading_shower_open_angle_numu_cc->Fill(leading_open_angle);
		}
		if(tpco_id == "nc")        {h_leading_shower_open_angle_nc->Fill(leading_open_angle); }
		if(tpco_id == "nc_pi0")    {h_leading_shower_open_angle_nc_pi0->Fill(leading_open_angle); }
		if(tpco_id == "unmatched") {h_leading_shower_open_angle_unmatched->Fill(leading_open_angle); }
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutOpenAngle2PlusShowerInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                            std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                            TH1D * h_leading_shower_open_angle_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }

		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		if(n_pfp_showers  < 2) {continue; }
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
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		h_leading_shower_open_angle_intime->Fill(leading_open_angle);
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutTrkVtx(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                        std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                        TH1D * h_trk_vtx_dist_nue_cc, TH1D * h_trk_vtx_dist_nue_cc_out_fv,
                                        TH1D * h_trk_vtx_dist_nue_cc_mixed,
                                        TH1D * h_trk_vtx_dist_numu_cc, TH1D * h_trk_vtx_dist_nc,
                                        TH1D * h_trk_vtx_dist_cosmic, TH1D * h_trk_vtx_dist_nc_pi0,
                                        TH1D * h_trk_vtx_dist_numu_cc_mixed, TH1D * h_trk_vtx_dist_other_mixed,
                                        TH1D * h_trk_vtx_dist_unmatched)
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
			}
		}//end pfp loop
		std::string tpco_id = tpco_classifier_v->at(i).first;
		if(tpco_id == "nue_cc_mixed")
		{
			if(has_track) {h_trk_vtx_dist_nue_cc_mixed->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "nue_cc_out_fv") {if(has_track) {h_trk_vtx_dist_nue_cc_out_fv->Fill(smallest_trk_vtx_dist); }}
		if(tpco_id == "numu_cc_mixed")
		{
			//if(has_track) {h_trk_vtx_dist_numu_cc_mixed->Fill(smallest_trk_vtx_dist); }
			if(has_track) {h_trk_vtx_dist_numu_cc->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "other_mixed")
		{
			if(has_track) {h_trk_vtx_dist_other_mixed->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "cosmic")
		{
			if(has_track) {h_trk_vtx_dist_cosmic->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_mec" ||
		   tpco_id == "nue_bar_cc_qe" || tpco_id == "nue_bar_cc_res" || tpco_id == "nue_bar_cc_dis" || tpco_id == "nue_bar_cc_coh" || tpco_id == "nue_bar_cc_mec")
		{
			if(has_track) {h_trk_vtx_dist_nue_cc->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "numu_cc_qe" || tpco_id == "numu_cc_res" || tpco_id == "numu_cc_dis" || tpco_id == "numu_cc_coh" || tpco_id == "numu_cc_mec")
		{
			if(has_track) {h_trk_vtx_dist_numu_cc->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "nc")
		{
			if(has_track) {h_trk_vtx_dist_nc->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "nc_pi0")
		{
			if(has_track) {h_trk_vtx_dist_nc_pi0->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "unmatched")
		{
			if(has_track) {h_trk_vtx_dist_unmatched->Fill(smallest_trk_vtx_dist); }
		}
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutTrkVtxInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, TH1D * h_trk_vtx_dist_intime)
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
			}
		}//end pfp loop
		if(has_track) {h_trk_vtx_dist_intime->Fill(smallest_trk_vtx_dist); }
	}//end tpco loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::TopologyPlots1(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                         std::vector<std::pair<int, std::string> > * passed_tpco,
                                         std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                         TH2D * h_pfp_track_shower_nue_cc_qe,
                                         TH2D * h_pfp_track_shower_nue_cc_out_fv,
                                         TH2D * h_pfp_track_shower_nue_cc_res,
                                         TH2D * h_pfp_track_shower_nue_cc_dis,
                                         TH2D * h_pfp_track_shower_nue_cc_coh,
                                         TH2D * h_pfp_track_shower_nue_cc_mec,
                                         TH2D * h_pfp_track_shower_nc,
                                         TH2D * h_pfp_track_shower_numu_cc_qe,
                                         TH2D * h_pfp_track_shower_numu_cc_res,
                                         TH2D * h_pfp_track_shower_numu_cc_dis,
                                         TH2D * h_pfp_track_shower_numu_cc_coh,
                                         TH2D * h_pfp_track_shower_numu_cc_mec,
                                         TH2D * h_pfp_track_shower_nc_pi0,
                                         TH2D * h_pfp_track_shower_nue_cc_mixed,
                                         TH2D * h_pfp_track_shower_numu_cc_mixed,
                                         TH2D * h_pfp_track_shower_cosmic,
                                         TH2D * h_pfp_track_shower_other_mixed,
                                         TH2D * h_pfp_track_shower_unmatched,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_qe,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_out_fv,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_res,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_dis,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_coh,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_mec,
                                         TH2D * h_leading_shower_mc_pdg_nc,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_qe,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_res,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_dis,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_coh,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_mec,
                                         TH2D * h_leading_shower_mc_pdg_nc_pi0,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_mixed,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_mixed,
                                         TH2D * h_leading_shower_mc_pdg_cosmic,
                                         TH2D * h_leading_shower_mc_pdg_other_mixed,
                                         TH2D * h_leading_shower_mc_pdg_unmatched,
                                         TH1D * h_pfp_track_nue_cc_qe,
                                         TH1D * h_pfp_track_nue_cc_out_fv,
                                         TH1D * h_pfp_track_nue_cc_res,
                                         TH1D * h_pfp_track_nue_cc_dis,
                                         TH1D * h_pfp_track_nue_cc_coh,
                                         TH1D * h_pfp_track_nue_cc_mec,
                                         TH1D * h_pfp_track_nc,
                                         TH1D * h_pfp_track_numu_cc_qe,
                                         TH1D * h_pfp_track_numu_cc_res,
                                         TH1D * h_pfp_track_numu_cc_dis,
                                         TH1D * h_pfp_track_numu_cc_coh,
                                         TH1D * h_pfp_track_numu_cc_mec,
                                         TH1D * h_pfp_track_nc_pi0,
                                         TH1D * h_pfp_track_nue_cc_mixed,
                                         TH1D * h_pfp_track_numu_cc_mixed,
                                         TH1D * h_pfp_track_cosmic,
                                         TH1D * h_pfp_track_other_mixed,
                                         TH1D * h_pfp_track_unmatched,
                                         TH1D * h_pfp_shower_nue_cc_qe,
                                         TH1D * h_pfp_shower_nue_cc_out_fv,
                                         TH1D * h_pfp_shower_nue_cc_res,
                                         TH1D * h_pfp_shower_nue_cc_dis,
                                         TH1D * h_pfp_shower_nue_cc_coh,
                                         TH1D * h_pfp_shower_nue_cc_mec,
                                         TH1D * h_pfp_shower_nc,
                                         TH1D * h_pfp_shower_numu_cc_qe,
                                         TH1D * h_pfp_shower_numu_cc_res,
                                         TH1D * h_pfp_shower_numu_cc_dis,
                                         TH1D * h_pfp_shower_numu_cc_coh,
                                         TH1D * h_pfp_shower_numu_cc_mec,
                                         TH1D * h_pfp_shower_nc_pi0,
                                         TH1D * h_pfp_shower_nue_cc_mixed,
                                         TH1D * h_pfp_shower_numu_cc_mixed,
                                         TH1D * h_pfp_shower_cosmic,
                                         TH1D * h_pfp_shower_other_mixed,
                                         TH1D * h_pfp_shower_unmatched
                                         )
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{

		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks  = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const std::string leading_origin = leading_shower.Origin();
		double leading_origin_int;
		if(leading_origin == "kBeamNeutrino") {leading_origin_int = 0.0; }
		if(leading_origin == "kCosmicRay")    {leading_origin_int = 1.0; }
		if(leading_origin == "kUnknown")      {leading_origin_int = 2.0; }
		const int leading_mc_pdg = leading_shower.MCPdgCode();
		double leading_pdg_int;
		if(leading_mc_pdg == 11)                            {leading_pdg_int = 0.0;  }
		if(leading_mc_pdg == -11)                           {leading_pdg_int = 1.0;  }
		if(leading_mc_pdg == 13)                            {leading_pdg_int = 2.0;  }
		if(leading_mc_pdg == -13)                           {leading_pdg_int = 3.0;  }
		if(leading_mc_pdg == 22)                            {leading_pdg_int = 4.0;  }
		if(leading_mc_pdg == 211 ||
		   leading_mc_pdg == -211)                          {leading_pdg_int = 5.0;  }
		if(leading_mc_pdg == 2212)                          {leading_pdg_int = 6.0;  }
		if(leading_mc_pdg == 2112)                          {leading_pdg_int = 7.0;  }
		if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
		   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
		   leading_mc_pdg == -321)                          {leading_pdg_int = 8.0; }
		if(leading_mc_pdg == 0)                             {leading_pdg_int = 9.0; }

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_pfp_track_shower_nue_cc_qe->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_qe->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_qe->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_qe->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_pfp_track_shower_nue_cc_out_fv->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_out_fv->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_out_fv->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_out_fv->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_pfp_track_shower_nue_cc_res->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_res->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_res->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_res->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_pfp_track_shower_nue_cc_dis->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_dis->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_dis->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_dis->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_pfp_track_shower_nue_cc_coh->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_coh->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_coh->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_coh->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_pfp_track_shower_nue_cc_mec->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_mec->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_mec->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_mec->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_pfp_track_shower_numu_cc_qe->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_qe->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_qe->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_qe->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_pfp_track_shower_numu_cc_res->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_res->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_res->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_res->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_pfp_track_shower_numu_cc_dis->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_dis->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_dis->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_dis->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_pfp_track_shower_numu_cc_coh->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_coh->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_coh->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_coh->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_pfp_track_shower_numu_cc_mec->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_mec->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_mec->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_mec->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nc")
		{
			h_pfp_track_shower_nc->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nc->Fill(n_pfp_tracks);
			h_pfp_shower_nc->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nc->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nc_pi0")
		{
			h_pfp_track_shower_nc_pi0->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nc_pi0->Fill(n_pfp_tracks);
			h_pfp_shower_nc_pi0->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nc_pi0->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_pfp_track_shower_nue_cc_mixed->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_mixed->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_mixed->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_mixed->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			h_pfp_track_shower_numu_cc_mixed->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_mixed->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_mixed->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_mixed->Fill(leading_origin_int, leading_pdg_int);
			// h_pfp_track_shower_numu_cc->Fill(n_pfp_tracks, n_pfp_showers);
			// h_pfp_track_numu_cc->Fill(n_pfp_tracks);
			// h_pfp_shower_numu_cc->Fill(n_pfp_showers);
			// h_leading_shower_mc_pdg_numu_cc->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "cosmic")
		{
			h_pfp_track_shower_cosmic->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_cosmic->Fill(n_pfp_tracks);
			h_pfp_shower_cosmic->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_cosmic->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "other_mixed")
		{
			h_pfp_track_shower_other_mixed->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_other_mixed->Fill(n_pfp_tracks);
			h_pfp_shower_other_mixed->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_other_mixed->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "unmatched")
		{
			h_pfp_track_shower_unmatched->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_unmatched->Fill(n_pfp_tracks);
			h_pfp_shower_unmatched->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_unmatched->Fill(leading_origin_int, leading_pdg_int);
		}
	}        //end loop tpc objects

}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::NumShowersOpenAngle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco,
                                              std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                              TH1D * h_pfp_shower_open_angle_nue_cc_qe,
                                              TH1D * h_pfp_shower_open_angle_nue_cc_out_fv,
                                              TH1D * h_pfp_shower_open_angle_nue_cc_res,
                                              TH1D * h_pfp_shower_open_angle_nue_cc_dis,
                                              TH1D * h_pfp_shower_open_angle_nue_cc_coh,
                                              TH1D * h_pfp_shower_open_angle_nue_cc_mec,
                                              TH1D * h_pfp_shower_open_angle_nc,
                                              TH1D * h_pfp_shower_open_angle_numu_cc_qe,
                                              TH1D * h_pfp_shower_open_angle_numu_cc_res,
                                              TH1D * h_pfp_shower_open_angle_numu_cc_dis,
                                              TH1D * h_pfp_shower_open_angle_numu_cc_coh,
                                              TH1D * h_pfp_shower_open_angle_numu_cc_mec,
                                              TH1D * h_pfp_shower_open_angle_nc_pi0,
                                              TH1D * h_pfp_shower_open_angle_nue_cc_mixed,
                                              TH1D * h_pfp_shower_open_angle_numu_cc_mixed,
                                              TH1D * h_pfp_shower_open_angle_cosmic,
                                              TH1D * h_pfp_shower_open_angle_other_mixed,
                                              TH1D * h_pfp_shower_open_angle_unmatched
                                              )
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id = tpco_classifier_v->at(i).first;
		if(tpco_id == "nue_cc_qe")      {h_pfp_shower_open_angle_nue_cc_qe->Fill(n_pfp_showers);  }
		if(tpco_id == "nue_cc_res")     {h_pfp_shower_open_angle_nue_cc_res->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_dis")     {h_pfp_shower_open_angle_nue_cc_dis->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_coh")     {h_pfp_shower_open_angle_nue_cc_coh->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_mec")     {h_pfp_shower_open_angle_nue_cc_mec->Fill(n_pfp_showers); }
		if(tpco_id == "nue_bar_cc_qe")      {h_pfp_shower_open_angle_nue_cc_qe->Fill(n_pfp_showers);  }
		if(tpco_id == "nue_bar_cc_res")     {h_pfp_shower_open_angle_nue_cc_res->Fill(n_pfp_showers); }
		if(tpco_id == "nue_bar_cc_dis")     {h_pfp_shower_open_angle_nue_cc_dis->Fill(n_pfp_showers); }
		if(tpco_id == "nue_bar_cc_coh")     {h_pfp_shower_open_angle_nue_cc_coh->Fill(n_pfp_showers); }
		if(tpco_id == "nue_bar_cc_mec")     {h_pfp_shower_open_angle_nue_cc_mec->Fill(n_pfp_showers); }

		if(tpco_id == "nue_cc_out_fv")  {h_pfp_shower_open_angle_nue_cc_out_fv->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_qe")     {h_pfp_shower_open_angle_numu_cc_qe->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_res")    {h_pfp_shower_open_angle_numu_cc_res->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_dis")    {h_pfp_shower_open_angle_numu_cc_dis->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_coh")    {h_pfp_shower_open_angle_numu_cc_coh->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_mec")    {h_pfp_shower_open_angle_numu_cc_mec->Fill(n_pfp_showers); }
		if(tpco_id == "nc")             {h_pfp_shower_open_angle_nc->Fill(n_pfp_showers); }
		if(tpco_id == "nc_pi0")         {h_pfp_shower_open_angle_nc_pi0->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_mixed")   {h_pfp_shower_open_angle_nue_cc_mixed->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_mixed")  {h_pfp_shower_open_angle_numu_cc_mixed->Fill(n_pfp_showers); }
		//if(tpco_id == "numu_cc_mixed")  {h_pfp_shower_open_angle_numu_cc->Fill(n_pfp_showers); }
		if(tpco_id == "cosmic")         {h_pfp_shower_open_angle_cosmic->Fill(n_pfp_showers); }
		if(tpco_id == "other_mixed")    {h_pfp_shower_open_angle_other_mixed->Fill(n_pfp_showers); }
		if(tpco_id == "unmatched")      {h_pfp_shower_open_angle_unmatched->Fill(n_pfp_showers); }
	}        //end loop tpc objects

}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::NumShowersOpenAngleInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, TH1D * h_pfp_shower_open_angle_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		h_pfp_shower_open_angle_intime->Fill(n_pfp_showers);
	}        //end loop tpc objects

}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::TopologyPlots2(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                         std::vector<std::pair<int, std::string> > * passed_tpco,
                                         std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                         TH2D * h_pfp_track_shower_nue_cc_qe,
                                         TH2D * h_pfp_track_shower_nue_cc_out_fv,
                                         TH2D * h_pfp_track_shower_nue_cc_res,
                                         TH2D * h_pfp_track_shower_nue_cc_dis,
                                         TH2D * h_pfp_track_shower_nue_cc_coh,
                                         TH2D * h_pfp_track_shower_nue_cc_mec,
                                         TH2D * h_pfp_track_shower_nc,
                                         TH2D * h_pfp_track_shower_numu_cc_qe,
                                         TH2D * h_pfp_track_shower_numu_cc_res,
                                         TH2D * h_pfp_track_shower_numu_cc_dis,
                                         TH2D * h_pfp_track_shower_numu_cc_coh,
                                         TH2D * h_pfp_track_shower_numu_cc_mec,
                                         TH2D * h_pfp_track_shower_nc_pi0,
                                         TH2D * h_pfp_track_shower_nue_cc_mixed,
                                         TH2D * h_pfp_track_shower_numu_cc_mixed,
                                         TH2D * h_pfp_track_shower_cosmic,
                                         TH2D * h_pfp_track_shower_other_mixed,
                                         TH2D * h_pfp_track_shower_unmatched,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_qe,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_out_fv,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_res,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_dis,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_coh,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_mec,
                                         TH2D * h_leading_shower_mc_pdg_nc,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_qe,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_res,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_dis,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_coh,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_mec,
                                         TH2D * h_leading_shower_mc_pdg_nc_pi0,
                                         TH2D * h_leading_shower_mc_pdg_nue_cc_mixed,
                                         TH2D * h_leading_shower_mc_pdg_numu_cc_mixed,
                                         TH2D * h_leading_shower_mc_pdg_cosmic,
                                         TH2D * h_leading_shower_mc_pdg_other_mixed,
                                         TH2D * h_leading_shower_mc_pdg_unmatched,
                                         TH1D * h_pfp_track_nue_cc_qe,
                                         TH1D * h_pfp_track_nue_cc_out_fv,
                                         TH1D * h_pfp_track_nue_cc_res,
                                         TH1D * h_pfp_track_nue_cc_dis,
                                         TH1D * h_pfp_track_nue_cc_coh,
                                         TH1D * h_pfp_track_nue_cc_mec,
                                         TH1D * h_pfp_track_nc,
                                         TH1D * h_pfp_track_numu_cc_qe,
                                         TH1D * h_pfp_track_numu_cc_res,
                                         TH1D * h_pfp_track_numu_cc_dis,
                                         TH1D * h_pfp_track_numu_cc_coh,
                                         TH1D * h_pfp_track_numu_cc_mec,
                                         TH1D * h_pfp_track_nc_pi0,
                                         TH1D * h_pfp_track_nue_cc_mixed,
                                         TH1D * h_pfp_track_numu_cc_mixed,
                                         TH1D * h_pfp_track_cosmic,
                                         TH1D * h_pfp_track_other_mixed,
                                         TH1D * h_pfp_track_unmatched,
                                         TH1D * h_pfp_shower_nue_cc_qe,
                                         TH1D * h_pfp_shower_nue_cc_out_fv,
                                         TH1D * h_pfp_shower_nue_cc_res,
                                         TH1D * h_pfp_shower_nue_cc_dis,
                                         TH1D * h_pfp_shower_nue_cc_coh,
                                         TH1D * h_pfp_shower_nue_cc_mec,
                                         TH1D * h_pfp_shower_nc,
                                         TH1D * h_pfp_shower_numu_cc_qe,
                                         TH1D * h_pfp_shower_numu_cc_res,
                                         TH1D * h_pfp_shower_numu_cc_dis,
                                         TH1D * h_pfp_shower_numu_cc_coh,
                                         TH1D * h_pfp_shower_numu_cc_mec,
                                         TH1D * h_pfp_shower_nc_pi0,
                                         TH1D * h_pfp_shower_nue_cc_mixed,
                                         TH1D * h_pfp_shower_numu_cc_mixed,
                                         TH1D * h_pfp_shower_cosmic,
                                         TH1D * h_pfp_shower_other_mixed,
                                         TH1D * h_pfp_shower_unmatched
                                         )
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{

		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks  = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const std::string leading_origin = leading_shower.Origin();
		double leading_origin_int;
		if(leading_origin == "kBeamNeutrino") {leading_origin_int = 0.0; }
		if(leading_origin == "kCosmicRay")    {leading_origin_int = 1.0; }
		if(leading_origin == "kUnknown")      {leading_origin_int = 2.0; }
		const int leading_mc_pdg = leading_shower.MCPdgCode();
		double leading_pdg_int;
		if(leading_mc_pdg == 11)                            {leading_pdg_int = 0.0;  }
		if(leading_mc_pdg == -11)                           {leading_pdg_int = 1.0;  }
		if(leading_mc_pdg == 13)                            {leading_pdg_int = 2.0;  }
		if(leading_mc_pdg == -13)                           {leading_pdg_int = 3.0;  }
		if(leading_mc_pdg == 22)                            {leading_pdg_int = 4.0;  }
		if(leading_mc_pdg == 211 ||
		   leading_mc_pdg == -211)                          {leading_pdg_int = 5.0;  }
		if(leading_mc_pdg == 2212)                          {leading_pdg_int = 6.0;  }
		if(leading_mc_pdg == 2112)                          {leading_pdg_int = 7.0;  }
		if(leading_mc_pdg == 130 || leading_mc_pdg == 310 ||
		   leading_mc_pdg == 311 || leading_mc_pdg == 321 ||
		   leading_mc_pdg == -321)                          {leading_pdg_int = 8.0; }
		if(leading_mc_pdg == 0)                             {leading_pdg_int = 9.0; }

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_pfp_track_shower_nue_cc_qe->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_qe->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_qe->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_qe->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_pfp_track_shower_nue_cc_out_fv->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_out_fv->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_out_fv->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_out_fv->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_pfp_track_shower_nue_cc_res->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_res->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_res->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_res->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_pfp_track_shower_nue_cc_dis->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_dis->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_dis->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_dis->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_pfp_track_shower_nue_cc_coh->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_coh->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_coh->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_coh->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_pfp_track_shower_nue_cc_mec->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_mec->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_mec->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_mec->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_pfp_track_shower_numu_cc_qe->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_qe->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_qe->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_qe->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_pfp_track_shower_numu_cc_res->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_res->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_res->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_res->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_pfp_track_shower_numu_cc_dis->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_dis->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_dis->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_dis->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_pfp_track_shower_numu_cc_coh->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_coh->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_coh->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_coh->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_pfp_track_shower_numu_cc_mec->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_mec->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_mec->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_mec->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nc")
		{
			h_pfp_track_shower_nc->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nc->Fill(n_pfp_tracks);
			h_pfp_shower_nc->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nc->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nc_pi0")
		{
			h_pfp_track_shower_nc_pi0->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nc_pi0->Fill(n_pfp_tracks);
			h_pfp_shower_nc_pi0->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nc_pi0->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_pfp_track_shower_nue_cc_mixed->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_mixed->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_mixed->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_mixed->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			h_pfp_track_shower_numu_cc_mixed->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_numu_cc_mixed->Fill(n_pfp_tracks);
			h_pfp_shower_numu_cc_mixed->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_numu_cc_mixed->Fill(leading_origin_int, leading_pdg_int);
			// h_pfp_track_shower_numu_cc->Fill(n_pfp_tracks, n_pfp_showers);
			// h_pfp_track_numu_cc->Fill(n_pfp_tracks);
			// h_pfp_shower_numu_cc->Fill(n_pfp_showers);
			// h_leading_shower_mc_pdg_numu_cc->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "cosmic")
		{
			h_pfp_track_shower_cosmic->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_cosmic->Fill(n_pfp_tracks);
			h_pfp_shower_cosmic->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_cosmic->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "other_mixed")
		{
			h_pfp_track_shower_other_mixed->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_other_mixed->Fill(n_pfp_tracks);
			h_pfp_shower_other_mixed->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_other_mixed->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "unmatched")
		{
			h_pfp_track_shower_unmatched->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_unmatched->Fill(n_pfp_tracks);
			h_pfp_shower_unmatched->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_unmatched->Fill(leading_origin_int, leading_pdg_int);
		}
	}        //end loop tpc objects

}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsVtxFlash(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                           std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                           TH1D * h_vtx_flash_nue_cc, TH1D * h_vtx_flash_nue_cc_out_fv, TH1D * h_vtx_flash_nue_cc_mixed,
                                           TH1D * h_vtx_flash_numu_cc, TH1D * h_vtx_flash_nc,
                                           TH1D * h_vtx_flash_cosmic, TH1D * h_vtx_flash_nc_pi0,
                                           TH1D * h_vtx_flash_numu_cc_mixed, TH1D * h_vtx_flash_other_mixed,
                                           TH1D * h_vtx_flash_unmatched)
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
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_vtx_flash_nue_cc_out_fv->Fill(distance);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "nc")
		{
			h_vtx_flash_nc->Fill(distance);
		}
		if(tpco_id == "nc_pi0")
		{
			h_vtx_flash_nc_pi0->Fill(distance);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_vtx_flash_nue_cc_mixed->Fill(distance);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_vtx_flash_numu_cc_mixed->Fill(distance);
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "cosmic")
		{
			h_vtx_flash_cosmic->Fill(distance);
		}
		if(tpco_id == "other_mixed")
		{
			h_vtx_flash_other_mixed->Fill(distance);
		}
		if(tpco_id == "unmatched")
		{
			h_vtx_flash_unmatched->Fill(distance);
		}
	}        //end loop tpc objects
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsVtxFlashInTime(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, TH1D * h_vtx_flash_intime)
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
		h_vtx_flash_intime->Fill(distance);
	}        //end loop tpc objects
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsShwrVtx(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                          std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                          TH1D * h_shwr_vtx_dist_nue_cc,
                                          TH1D * h_shwr_vtx_dist_nue_cc_out_fv,
                                          TH1D * h_shwr_vtx_dist_nue_cc_mixed,
                                          TH1D * h_shwr_vtx_dist_numu_cc,
                                          TH1D * h_shwr_vtx_dist_nc,
                                          TH1D * h_shwr_vtx_dist_cosmic,
                                          TH1D * h_shwr_vtx_dist_nc_pi0,
                                          TH1D * h_shwr_vtx_dist_numu_cc_mixed,
                                          TH1D * h_shwr_vtx_dist_other_mixed,
                                          TH1D * h_shwr_vtx_dist_unmatched     )
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
		std::string tpco_id = tpco_classifier_v->at(i).first;

		// int leading_index   = tpco_classifier_v->at(i).second;
		// auto const leading_shower = tpc_obj.GetParticle(leading_index);
		// const double leading_vtx_x = leading_shower.pfpVtxX();
		// const double leading_vtx_y = leading_shower.pfpVtxY();
		// const double leading_vtx_z = leading_shower.pfpVtxZ();

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

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_shwr_vtx_dist_nue_cc->Fill(min_distance);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_shwr_vtx_dist_nue_cc_out_fv->Fill(min_distance);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_shwr_vtx_dist_nue_cc->Fill(min_distance);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_shwr_vtx_dist_nue_cc->Fill(min_distance);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_shwr_vtx_dist_nue_cc->Fill(min_distance);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_shwr_vtx_dist_nue_cc->Fill(min_distance);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_shwr_vtx_dist_numu_cc->Fill(min_distance);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_shwr_vtx_dist_numu_cc->Fill(min_distance);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_shwr_vtx_dist_numu_cc->Fill(min_distance);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_shwr_vtx_dist_numu_cc->Fill(min_distance);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_shwr_vtx_dist_numu_cc->Fill(min_distance);
		}
		if(tpco_id == "nc")
		{
			h_shwr_vtx_dist_nc->Fill(min_distance);
		}
		if(tpco_id == "nc_pi0")
		{
			h_shwr_vtx_dist_nc_pi0->Fill(min_distance);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_shwr_vtx_dist_nue_cc_mixed->Fill(min_distance);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_shwr_vtx_dist_numu_cc_mixed->Fill(distance);
			h_shwr_vtx_dist_numu_cc->Fill(min_distance);
		}
		if(tpco_id == "cosmic")
		{
			h_shwr_vtx_dist_cosmic->Fill(min_distance);
		}
		if(tpco_id == "other_mixed")
		{
			h_shwr_vtx_dist_other_mixed->Fill(min_distance);
		}
		if(tpco_id == "unmatched")
		{
			h_shwr_vtx_dist_unmatched->Fill(min_distance);
		}
	}        //end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsShwrVtxInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, TH1D * h_shwr_vtx_dist_intime)
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
		h_shwr_vtx_dist_intime->Fill(min_distance);
	}        //end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutHitThreshold(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                              std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                              double mc_nu_energy, double mc_ele_energy,
                                              TH2D * h_shwr_hits_nu_eng, TH2D * h_shwr_hits_ele_eng)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;

		int pfp_shower_hits = 0;
		if(tpco_id == "nue_cc_qe"  ||
		   tpco_id == "nue_cc_res" ||
		   tpco_id == "nue_cc_dis" ||
		   tpco_id == "nue_cc_mec" ||
		   tpco_id == "nue_cc_coh" ||
		   tpco_id == "nue_bar_cc_qe"  ||
		   tpco_id == "nue_bar_cc_res" ||
		   tpco_id == "nue_bar_cc_dis" ||
		   tpco_id == "nue_bar_cc_mec" ||
		   tpco_id == "nue_bar_cc_coh")
		{
			for(int j = 0; j < n_pfp; j++)
			{
				auto const part = tpc_obj.GetParticle(j);
				const int mc_pdg_code = part.MCPdgCode();
				if(mc_pdg_code == 11 || mc_pdg_code == -11) {pfp_shower_hits = part.NumPFPHits(); }
			}//end loop pfp objects
			h_shwr_hits_nu_eng->Fill(mc_nu_energy, pfp_shower_hits);
			h_shwr_hits_ele_eng->Fill(mc_ele_energy, pfp_shower_hits);
		}//if classified as pure nue cc in fv, i.e. signal
	} //end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions::TopologyEfficiency(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                             std::vector<std::pair<std::string, int> > * tpco_classifier_v,
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
		std::string tpco_id = tpco_classifier_v->at(i).first;
		//signal
		if(n_pfp_showers == 1)
		{
			if(tpco_id == "nue_cc_qe"  ||
			   tpco_id == "nue_cc_res" ||
			   tpco_id == "nue_cc_dis" ||
			   tpco_id == "nue_cc_mec" ||
			   tpco_id == "nue_cc_coh" ||
			   tpco_id == "nue_bar_cc_qe"  ||
			   tpco_id == "nue_bar_cc_res" ||
			   tpco_id == "nue_bar_cc_dis" ||
			   tpco_id == "nue_bar_cc_mec" ||
			   tpco_id == "nue_bar_cc_coh")
			{
				_1_shwr->at(0) += 1;
			}
			//not signal
			if(tpco_id != "nue_cc_qe"  &&
			   tpco_id != "nue_cc_res" &&
			   tpco_id != "nue_cc_dis" &&
			   tpco_id != "nue_cc_mec" &&
			   tpco_id != "nue_cc_coh" &&
			   tpco_id != "nue_bar_cc_qe"  &&
			   tpco_id != "nue_bar_cc_res" &&
			   tpco_id != "nue_bar_cc_dis" &&
			   tpco_id != "nue_bar_cc_mec" &&
			   tpco_id != "nue_bar_cc_coh" )
			{
				_1_shwr->at(1) += 1;
			}
		}
		if(n_pfp_showers == 2)
		{
			if(tpco_id == "nue_cc_qe"  ||
			   tpco_id == "nue_cc_res" ||
			   tpco_id == "nue_cc_dis" ||
			   tpco_id == "nue_cc_mec" ||
			   tpco_id == "nue_cc_coh" ||
			   tpco_id == "nue_bar_cc_qe"  ||
			   tpco_id == "nue_bar_cc_res" ||
			   tpco_id == "nue_bar_cc_dis" ||
			   tpco_id == "nue_bar_cc_mec" ||
			   tpco_id == "nue_bar_cc_coh")
			{
				_2_shwr->at(0) += 1;
			}
			//not signal
			if(tpco_id != "nue_cc_qe"  &&
			   tpco_id != "nue_cc_res" &&
			   tpco_id != "nue_cc_dis" &&
			   tpco_id != "nue_cc_mec" &&
			   tpco_id != "nue_cc_coh" &&
			   tpco_id != "nue_bar_cc_qe"  &&
			   tpco_id != "nue_bar_cc_res" &&
			   tpco_id != "nue_bar_cc_dis" &&
			   tpco_id != "nue_bar_cc_coh" &&
			   tpco_id != "nue_bar_cc_mec" )
			{
				_2_shwr->at(1) += 1;
			}
		}
		if(n_pfp_showers == 3)
		{
			if(tpco_id == "nue_cc_qe"  ||
			   tpco_id == "nue_cc_res" ||
			   tpco_id == "nue_cc_dis" ||
			   tpco_id == "nue_cc_mec" ||
			   tpco_id == "nue_cc_coh" ||
			   tpco_id == "nue_bar_cc_qe"  ||
			   tpco_id == "nue_bar_cc_res" ||
			   tpco_id == "nue_bar_cc_dis" ||
			   tpco_id == "nue_bar_cc_mec" ||
			   tpco_id == "nue_bar_cc_coh")
			{
				_3_shwr->at(0) += 1;
			}
			//not signal
			if(tpco_id != "nue_cc_qe"  &&
			   tpco_id != "nue_cc_res" &&
			   tpco_id != "nue_cc_dis" &&
			   tpco_id != "nue_cc_mec" &&
			   tpco_id != "nue_cc_coh" &&
			   tpco_id != "nue_bar_cc_qe"  &&
			   tpco_id != "nue_bar_cc_res" &&
			   tpco_id != "nue_bar_cc_dis" &&
			   tpco_id != "nue_bar_cc_mec" &&
			   tpco_id != "nue_bar_cc_coh" )
			{
				_3_shwr->at(1) += 1;
			}
		}
		if(n_pfp_showers >= 4)
		{
			if(tpco_id == "nue_cc_qe"  ||
			   tpco_id == "nue_cc_res" ||
			   tpco_id == "nue_cc_dis" ||
			   tpco_id == "nue_cc_mec" ||
			   tpco_id == "nue_cc_coh" ||
			   tpco_id == "nue_bar_cc_qe"  ||
			   tpco_id == "nue_bar_cc_res" ||
			   tpco_id == "nue_bar_cc_dis" ||
			   tpco_id == "nue_bar_cc_mec" ||
			   tpco_id == "nue_bar_cc_coh")
			{
				_4_shwr->at(0) += 1;
			}
			//not signal
			if(tpco_id != "nue_cc_qe"  &&
			   tpco_id != "nue_cc_res" &&
			   tpco_id != "nue_cc_dis" &&
			   tpco_id != "nue_cc_mec" &&
			   tpco_id != "nue_cc_coh" &&
			   tpco_id != "nue_bar_cc_qe"  &&
			   tpco_id != "nue_bar_cc_res" &&
			   tpco_id != "nue_bar_cc_dis" &&
			   tpco_id != "nue_bar_cc_mec" &&
			   tpco_id != "nue_bar_cc_coh" )
			{
				_4_shwr->at(1) += 1;
			}
		}
		if(n_pfp_tracks == 0)
		{
			if(tpco_id == "nue_cc_qe"  ||
			   tpco_id == "nue_cc_res" ||
			   tpco_id == "nue_cc_dis" ||
			   tpco_id == "nue_cc_mec" ||
			   tpco_id == "nue_cc_coh" ||
			   tpco_id == "nue_bar_cc_qe"  ||
			   tpco_id == "nue_bar_cc_res" ||
			   tpco_id == "nue_bar_cc_dis" ||
			   tpco_id == "nue_bar_cc_mec" ||
			   tpco_id == "nue_bar_cc_coh")
			{
				no_track->at(0) += 1;
			}
			//not signal
			if(tpco_id != "nue_cc_qe"  &&
			   tpco_id != "nue_cc_res" &&
			   tpco_id != "nue_cc_dis" &&
			   tpco_id != "nue_cc_mec" &&
			   tpco_id != "nue_cc_coh" &&
			   tpco_id != "nue_bar_cc_qe"  &&
			   tpco_id != "nue_bar_cc_res" &&
			   tpco_id != "nue_bar_cc_dis" &&
			   tpco_id != "nue_bar_cc_mec" &&
			   tpco_id != "nue_bar_cc_coh" )
			{
				no_track->at(1) += 1;
			}
		}
		if(n_pfp_tracks >= 1)
		{
			if(tpco_id == "nue_cc_qe"  ||
			   tpco_id == "nue_cc_res" ||
			   tpco_id == "nue_cc_dis" ||
			   tpco_id == "nue_cc_mec" ||
			   tpco_id == "nue_cc_coh" ||
			   tpco_id == "nue_bar_cc_qe"  ||
			   tpco_id == "nue_bar_cc_res" ||
			   tpco_id == "nue_bar_cc_dis" ||
			   tpco_id == "nue_bar_cc_mec" ||
			   tpco_id == "nue_bar_cc_coh")
			{
				has_track->at(0) += 1;
			}
			//not signal
			if(tpco_id != "nue_cc_qe"  &&
			   tpco_id != "nue_cc_res" &&
			   tpco_id != "nue_cc_dis" &&
			   tpco_id != "nue_cc_mec" &&
			   tpco_id != "nue_cc_coh" &&
			   tpco_id != "nue_bar_cc_qe"  &&
			   tpco_id != "nue_bar_cc_res" &&
			   tpco_id != "nue_bar_cc_dis" &&
			   tpco_id != "nue_bar_cc_mec" &&
			   tpco_id != "nue_bar_cc_coh" )
			{
				has_track->at(1) += 1;
			}
		}
	}//end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions::SequentialTrueEnergyPlots(int mc_nu_id, double mc_nu_vtx_x, double mc_nu_vtx_y, double mc_nu_vtx_z,
                                                    std::vector<double> fv_boundary_v,
                                                    std::vector<int> * tabulated_origins, double mc_nu_energy,
                                                    double mc_ele_energy, TH1D * h_selected_nu_energy, TH1D * h_selected_ele_energy)
{
	//this checks if there is a true nue/nue-bar CC event and a selected nue_cc signal event, true in FV
	selection_cuts _functions_instance;
	if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins->at(0) == 1) {
		if(_functions_instance.selection_cuts::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v) == true) {
			h_selected_nu_energy->Fill(mc_nu_energy);
			h_selected_ele_energy->Fill(mc_ele_energy);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::ChargeShare(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                      std::vector<std::pair<std::string, int> > * tpco_classifier_v, TH1D * h_charge_share_nue_cc_mixed)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id = tpco_classifier_v->at(i).first;
		if(tpco_id == "nue_cc_mixed")
		{
			int neutrino_charge = 0;
			int cosmic_charge = 0;
			for(int j = 0; j < n_pfp; j++)
			{
				auto const part = tpc_obj.GetParticle(j);
				const int n_pfp_hits = part.NumPFPHits();
				if(part.Origin() == "kBeamNeutrino") {neutrino_charge += n_pfp_hits; }
				if(part.Origin() == "kCosmicRay")    {cosmic_charge += n_pfp_hits; }
			}
			const double neutrino_charge_share = double(neutrino_charge) / double(neutrino_charge + cosmic_charge);
			h_charge_share_nue_cc_mixed->Fill(neutrino_charge_share);
		}//end if tpco == nue_cc_mixed
	}//end loop tpco
}
//***************************************************************************
//***************************************************************************
void selection_functions::FlashTot0(std::vector< double> largest_flash_v, double mc_nu_time, int mc_nu_id, std::vector<int> * tabulated_origins,
                                    std::vector<double> fv_boundary_v,
                                    double vtxX, double vtxY, double vtxZ, TH1D * h_flash_t0_diff)
{
	selection_cuts _functions_instance;
	if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins->at(0) == 1)
	{
		const bool InFV = _functions_instance.selection_cuts::in_fv(vtxX, vtxY, vtxZ, fv_boundary_v);
		if(InFV == true)
		{
			double largest_op_flash_time = largest_flash_v.at(4);
			double difference = largest_op_flash_time - (mc_nu_time / 1000);
			//std::cout << "Largest OpFlash Time: " << largest_op_flash_time << " , MC Nu Time: " << mc_nu_time / 1000 << std::endl;
			//std::cout << "Largest Flash Time - MC Nu Time: " << difference << std::endl;
			h_flash_t0_diff->Fill(difference);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::dEdxVsOpenAngle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                          std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                          TH2D * h_dedx_open_angle_nue_cc,
                                          TH2D * h_dedx_open_angle_nue_cc_out_fv,
                                          TH2D * h_dedx_open_angle_nue_cc_mixed,
                                          TH2D * h_dedx_open_angle_numu_cc,
                                          TH2D * h_dedx_open_angle_numu_cc_mixed,
                                          TH2D * h_dedx_open_angle_nc,
                                          TH2D * h_dedx_open_angle_nc_pi0,
                                          TH2D * h_dedx_open_angle_cosmic,
                                          TH2D * h_dedx_open_angle_other_mixed,
                                          TH2D * h_dedx_open_angle_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id = tpco_classifier_v->at(i).first;
		int leading_index   = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_dedx_open_angle_nue_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_dedx_open_angle_nue_cc_out_fv->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_dedx_open_angle_nue_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_dedx_open_angle_nue_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_dedx_open_angle_nue_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_dedx_open_angle_nue_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_dedx_open_angle_numu_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_dedx_open_angle_numu_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_dedx_open_angle_numu_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_dedx_open_angle_numu_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_dedx_open_angle_numu_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nc")
		{
			h_dedx_open_angle_nc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nc_pi0")
		{
			h_dedx_open_angle_nc_pi0->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_dedx_open_angle_nue_cc_mixed->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_dedx_open_angle_numu_cc_mixed->Fill(leading_dedx, leading_open_angle);
			h_dedx_open_angle_numu_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "cosmic")
		{
			h_dedx_open_angle_cosmic->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "other_mixed")
		{
			h_dedx_open_angle_other_mixed->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "unmatched")
		{
			h_dedx_open_angle_unmatched->Fill(leading_dedx, leading_open_angle);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::dEdxVsOpenAngleInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, TH2D * h_dedx_open_angle_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
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
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		h_dedx_open_angle_intime->Fill(leading_dedx * (242.72 / 196.979), leading_open_angle);
	}
}
//***************************************************************************
//***************************************************************************
//shower hits vs shower length, to see if we can better use the hit threshold cut to remove less signal?
void selection_functions::ShowerLengthvsHits(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                             std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                             TH2D * h_shwr_len_hits_nue_cc,
                                             TH2D * h_shwr_len_hits_nue_cc_out_fv,
                                             TH2D * h_shwr_len_hits_nue_cc_mixed,
                                             TH2D * h_shwr_len_hits_numu_cc,
                                             TH2D * h_shwr_len_hits_numu_cc_mixed,
                                             TH2D * h_shwr_len_hits_nc,
                                             TH2D * h_shwr_len_hits_nc_pi0,
                                             TH2D * h_shwr_len_hits_cosmic,
                                             TH2D * h_shwr_len_hits_other_mixed,
                                             TH2D * h_shwr_len_hits_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id = tpco_classifier_v->at(i).first;
		int leading_index   = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int leading_hits = leading_shower.NumPFPHits();
		const double leading_length = leading_shower.pfpLength();
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_shwr_len_hits_nue_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_shwr_len_hits_nue_cc_out_fv->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_shwr_len_hits_nue_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_shwr_len_hits_nue_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_shwr_len_hits_nue_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_shwr_len_hits_nue_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_shwr_len_hits_numu_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_shwr_len_hits_numu_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_shwr_len_hits_numu_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_shwr_len_hits_numu_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_shwr_len_hits_numu_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nc")
		{
			h_shwr_len_hits_nc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nc_pi0")
		{
			h_shwr_len_hits_nc_pi0->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_shwr_len_hits_nue_cc_mixed->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_shwr_len_hits_numu_cc_mixed->Fill(leading_length, leading_hits);
			h_shwr_len_hits_numu_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "cosmic")
		{
			h_shwr_len_hits_cosmic->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "other_mixed")
		{
			h_shwr_len_hits_other_mixed->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "unmatched")
		{
			h_shwr_len_hits_unmatched->Fill(leading_length, leading_hits);
		}
	}
}
//***************************************************************************
//***************************************************************************
//shower hits vs shower length, to see if we can better use the hit threshold cut to remove less signal?
void selection_functions::ShowerLengthvsHitsInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                   TH2D * h_shwr_len_hits_intime)
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
		const int leading_hits = leading_shower.NumPFPHits();
		const double leading_length = leading_shower.pfpLength();
		h_shwr_len_hits_intime->Fill(leading_length, leading_hits);
	}
}
//***************************************************************************
void selection_functions::SecondaryShowersDist(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                               std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                               TH1D * h_second_shwr_dist_nue_cc,
                                               TH1D * h_second_shwr_dist_nue_cc_out_fv,
                                               TH1D * h_second_shwr_dist_nue_cc_mixed,
                                               TH1D * h_second_shwr_dist_numu_cc,
                                               TH1D * h_second_shwr_dist_numu_cc_mixed,
                                               TH1D * h_second_shwr_dist_nc,
                                               TH1D * h_second_shwr_dist_nc_pi0,
                                               TH1D * h_second_shwr_dist_cosmic,
                                               TH1D * h_second_shwr_dist_other_mixed,
                                               TH1D * h_second_shwr_dist_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		if(n_pfp_showers <= 1) {continue; }
		const double tpco_vtx_x = tpc_obj.pfpVtxX();
		const double tpco_vtx_y = tpc_obj.pfpVtxY();
		const double tpco_vtx_z = tpc_obj.pfpVtxZ();
		std::string tpco_id = tpco_classifier_v->at(i).first;
		int leading_index   = tpco_classifier_v->at(i).second;
		//auto const leading_shower = tpc_obj.GetParticle(leading_index);
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

			if(pfp_pdg == 11 && distance > max_distance) {max_distance = distance; }
		}
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_second_shwr_dist_nue_cc->Fill(max_distance);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_second_shwr_dist_nue_cc_out_fv->Fill(max_distance);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_second_shwr_dist_nue_cc->Fill(max_distance);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_second_shwr_dist_nue_cc->Fill(max_distance);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_second_shwr_dist_nue_cc->Fill(max_distance);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_second_shwr_dist_nue_cc->Fill(max_distance);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_second_shwr_dist_numu_cc->Fill(max_distance);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_second_shwr_dist_numu_cc->Fill(max_distance);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_second_shwr_dist_numu_cc->Fill(max_distance);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_second_shwr_dist_numu_cc->Fill(max_distance);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_second_shwr_dist_numu_cc->Fill(max_distance);
		}
		if(tpco_id == "nc")
		{
			h_second_shwr_dist_nc->Fill(max_distance);
		}
		if(tpco_id == "nc_pi0")
		{
			h_second_shwr_dist_nc_pi0->Fill(max_distance);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_second_shwr_dist_nue_cc_mixed->Fill(max_distance);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_second_shwr_dist_numu_cc_mixed->Fill(distance);
			h_second_shwr_dist_numu_cc->Fill(max_distance);
		}
		if(tpco_id == "cosmic")
		{
			h_second_shwr_dist_cosmic->Fill(max_distance);
		}
		if(tpco_id == "other_mixed")
		{
			h_second_shwr_dist_other_mixed->Fill(max_distance);
		}
		if(tpco_id == "unmatched")
		{
			h_second_shwr_dist_unmatched->Fill(max_distance);
		}
	}//end loop tpco
}//end function
//***************************************************************************
void selection_functions::SecondaryShowersDistInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                     TH1D * h_second_shwr_dist_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		if(n_pfp_showers <= 1) {continue; }
		const double tpco_vtx_x = tpc_obj.pfpVtxX();
		const double tpco_vtx_y = tpc_obj.pfpVtxY();
		const double tpco_vtx_z = tpc_obj.pfpVtxZ();
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
		//auto const leading_shower = tpc_obj.GetParticle(leading_index);
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
			{
				max_distance = distance;
			}        //end if reco shower
		}        //end loop pfp
		h_second_shwr_dist_intime->Fill(max_distance);
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::HitLengthRatio(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                         std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                         TH1D * h_hit_length_ratio_nue_cc,
                                         TH1D * h_hit_length_ratio_nue_cc_out_fv,
                                         TH1D * h_hit_length_ratio_nue_cc_mixed,
                                         TH1D * h_hit_length_ratio_numu_cc,
                                         TH1D * h_hit_length_ratio_numu_cc_mixed,
                                         TH1D * h_hit_length_ratio_nc,
                                         TH1D * h_hit_length_ratio_nc_pi0,
                                         TH1D * h_hit_length_ratio_cosmic,
                                         TH1D * h_hit_length_ratio_other_mixed,
                                         TH1D * h_hit_length_ratio_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id = tpco_classifier_v->at(i).first;
		int leading_index   = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int pfp_pdg = leading_shower.PFParticlePdgCode();
		const double pfp_hits = leading_shower.NumPFPHits();
		const double pfp_length = leading_shower.pfpLength();
		const double pfp_hits_length_ratio = (pfp_hits / pfp_length);
		if(pfp_pdg == 11)
		{
			if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
			{
				h_hit_length_ratio_nue_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_out_fv")
			{
				h_hit_length_ratio_nue_cc_out_fv->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
			{
				h_hit_length_ratio_nue_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
			{
				h_hit_length_ratio_nue_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
			{
				h_hit_length_ratio_nue_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
			{
				h_hit_length_ratio_nue_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "numu_cc_qe")
			{
				h_hit_length_ratio_numu_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "numu_cc_res")
			{
				h_hit_length_ratio_numu_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "numu_cc_dis")
			{
				h_hit_length_ratio_numu_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "numu_cc_coh")
			{
				h_hit_length_ratio_numu_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "numu_cc_mec")
			{
				h_hit_length_ratio_numu_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nc")
			{
				h_hit_length_ratio_nc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nc_pi0")
			{
				h_hit_length_ratio_nc_pi0->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_mixed")
			{
				h_hit_length_ratio_nue_cc_mixed->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "numu_cc_mixed")
			{
				//h_hit_length_ratio_numu_cc_mixed->Fill(pfp_hits_length_ratio);
				h_hit_length_ratio_numu_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "cosmic")
			{
				h_hit_length_ratio_cosmic->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "other_mixed")
			{
				h_hit_length_ratio_other_mixed->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "unmatched")
			{
				h_hit_length_ratio_unmatched->Fill(pfp_hits_length_ratio);
			}
		}                //end if reco shower
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::HitLengthRatioInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                               TH1D * h_hit_length_ratio_intime)
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
		const int pfp_pdg = leading_shower.PFParticlePdgCode();
		const double pfp_hits = leading_shower.NumPFPHits();
		const double pfp_length = leading_shower.pfpLength();
		const double pfp_hits_length_ratio = (pfp_hits / pfp_length);
		if(pfp_pdg == 11)
		{
			h_hit_length_ratio_intime->Fill(pfp_hits_length_ratio);
		}                //end if reco shower
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::TrackLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                      std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                      TH1D * h_trk_length_nue_cc,
                                      TH1D * h_trk_length_nue_cc_out_fv,
                                      TH1D * h_trk_length_nue_cc_mixed,
                                      TH1D * h_trk_length_numu_cc,
                                      TH1D * h_trk_length_numu_cc_mixed,
                                      TH1D * h_trk_length_nc,
                                      TH1D * h_trk_length_nc_pi0,
                                      TH1D * h_trk_length_cosmic,
                                      TH1D * h_trk_length_other_mixed,
                                      TH1D * h_trk_length_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		if(n_pfp_tracks == 0) {continue; }
		std::vector< double > trk_length_v;
		for(int i = 0; i < n_pfp; i++)
		{
			auto const pfp = tpc_obj.GetParticle(i);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg == 13)
			{
				const double trk_length = pfp.pfpLength();
				trk_length_v.push_back(trk_length);
			}
		}
		std::string tpco_id = tpco_classifier_v->at(i).first;

		for(const double trk_length : trk_length_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_res")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_dis")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_coh")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_mec")     {h_trk_length_nue_cc->Fill(trk_length); }

			if(tpco_id == "nue_bar_cc_qe")      {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_bar_cc_res")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_bar_cc_dis")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_bar_cc_coh")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_bar_cc_mec")     {h_trk_length_nue_cc->Fill(trk_length); }

			if(tpco_id == "nue_cc_out_fv")  {h_trk_length_nue_cc_out_fv->Fill(trk_length); }
			if(tpco_id == "numu_cc_qe")     {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_res")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_dis")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_coh")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_mec")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "nc")             {h_trk_length_nc->Fill(trk_length); }
			if(tpco_id == "nc_pi0")         {h_trk_length_nc_pi0->Fill(trk_length); }
			if(tpco_id == "nue_cc_mixed")   {h_trk_length_nue_cc_mixed->Fill(trk_length); }
			//if(tpco_id == "numu_cc_mixed")  {h_trk_length_numu_cc_mixed->Fill(trk_length); }
			if(tpco_id == "numu_cc_mixed")  {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "cosmic")         {h_trk_length_cosmic->Fill(trk_length); }
			if(tpco_id == "other_mixed")    {h_trk_length_other_mixed->Fill(trk_length); }
			if(tpco_id == "unmatched")      {h_trk_length_unmatched->Fill(trk_length); }
		}
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::LongestTrackLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                             std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                             TH1D * h_trk_length_nue_cc,
                                             TH1D * h_trk_length_nue_cc_out_fv,
                                             TH1D * h_trk_length_nue_cc_mixed,
                                             TH1D * h_trk_length_numu_cc,
                                             TH1D * h_trk_length_numu_cc_mixed,
                                             TH1D * h_trk_length_nc,
                                             TH1D * h_trk_length_nc_pi0,
                                             TH1D * h_trk_length_cosmic,
                                             TH1D * h_trk_length_other_mixed,
                                             TH1D * h_trk_length_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		if(n_pfp_tracks == 0) {continue; }
		std::vector< double > trk_length_v;
		double longest_track = 0;
		for(int i = 0; i < n_pfp; i++)
		{
			auto const pfp = tpc_obj.GetParticle(i);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg == 13)
			{
				const double trk_length = pfp.pfpLength();
				if(trk_length > longest_track)
				{
					trk_length_v.clear();
					trk_length_v.push_back(trk_length);
					longest_track = trk_length;
				}
			}
		}
		std::string tpco_id = tpco_classifier_v->at(i).first;

		for(const double trk_length : trk_length_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_res")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_dis")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_coh")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_mec")     {h_trk_length_nue_cc->Fill(trk_length); }

			if(tpco_id == "nue_bar_cc_qe")      {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_bar_cc_res")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_bar_cc_dis")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_bar_cc_coh")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_bar_cc_mec")     {h_trk_length_nue_cc->Fill(trk_length); }

			if(tpco_id == "nue_cc_out_fv")  {h_trk_length_nue_cc_out_fv->Fill(trk_length); }
			if(tpco_id == "numu_cc_qe")     {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_res")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_dis")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_coh")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_mec")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "nc")             {h_trk_length_nc->Fill(trk_length); }
			if(tpco_id == "nc_pi0")         {h_trk_length_nc_pi0->Fill(trk_length); }
			if(tpco_id == "nue_cc_mixed")   {h_trk_length_nue_cc_mixed->Fill(trk_length); }
			//if(tpco_id == "numu_cc_mixed")  {h_trk_length_numu_cc_mixed->Fill(trk_length); }
			if(tpco_id == "numu_cc_mixed")  {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "cosmic")         {h_trk_length_cosmic->Fill(trk_length); }
			if(tpco_id == "other_mixed")    {h_trk_length_other_mixed->Fill(trk_length); }
			if(tpco_id == "unmatched")      {h_trk_length_unmatched->Fill(trk_length); }
		}
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::ShowerLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                       std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                       TH1D * h_shwr_length_nue_cc,
                                       TH1D * h_shwr_length_nue_cc_out_fv,
                                       TH1D * h_shwr_length_nue_cc_mixed,
                                       TH1D * h_shwr_length_numu_cc,
                                       TH1D * h_shwr_length_numu_cc_mixed,
                                       TH1D * h_shwr_length_nc,
                                       TH1D * h_shwr_length_nc_pi0,
                                       TH1D * h_shwr_length_cosmic,
                                       TH1D * h_shwr_length_other_mixed,
                                       TH1D * h_shwr_length_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::vector< double > shwr_length_v;
		for(int i = 0; i < n_pfp; i++)
		{
			auto const pfp = tpc_obj.GetParticle(i);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg == 11)
			{
				const double shwr_length = pfp.pfpLength();
				shwr_length_v.push_back(shwr_length);
			}
		}
		std::string tpco_id = tpco_classifier_v->at(i).first;

		for(const double shwr_length : shwr_length_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_res")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_dis")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_coh")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_mec")     {h_shwr_length_nue_cc->Fill(shwr_length); }

			if(tpco_id == "nue_bar_cc_qe")      {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_bar_cc_res")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_bar_cc_dis")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_bar_cc_coh")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_bar_cc_mec")     {h_shwr_length_nue_cc->Fill(shwr_length); }

			if(tpco_id == "nue_cc_out_fv")  {h_shwr_length_nue_cc_out_fv->Fill(shwr_length); }
			if(tpco_id == "numu_cc_qe")     {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_res")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_dis")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_coh")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_mec")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "nc")             {h_shwr_length_nc->Fill(shwr_length); }
			if(tpco_id == "nc_pi0")         {h_shwr_length_nc_pi0->Fill(shwr_length); }
			if(tpco_id == "nue_cc_mixed")   {h_shwr_length_nue_cc_mixed->Fill(shwr_length); }
			//if(tpco_id == "numu_cc_mixed")  {h_shwr_length_numu_cc_mixed->Fill(shwr_length); }
			if(tpco_id == "numu_cc_mixed")  {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "cosmic")         {h_shwr_length_cosmic->Fill(shwr_length); }
			if(tpco_id == "other_mixed")    {h_shwr_length_other_mixed->Fill(shwr_length); }
			if(tpco_id == "unmatched")      {h_shwr_length_unmatched->Fill(shwr_length); }
		}
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::LongestShowerLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                              std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                              TH1D * h_shwr_length_nue_cc,
                                              TH1D * h_shwr_length_nue_cc_out_fv,
                                              TH1D * h_shwr_length_nue_cc_mixed,
                                              TH1D * h_shwr_length_numu_cc,
                                              TH1D * h_shwr_length_numu_cc_mixed,
                                              TH1D * h_shwr_length_nc,
                                              TH1D * h_shwr_length_nc_pi0,
                                              TH1D * h_shwr_length_cosmic,
                                              TH1D * h_shwr_length_other_mixed,
                                              TH1D * h_shwr_length_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::vector< double > shwr_length_v;
		double longest_shower = 0;
		for(int i = 0; i < n_pfp; i++)
		{
			auto const pfp = tpc_obj.GetParticle(i);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg == 11)
			{
				const double shwr_length = pfp.pfpLength();
				if(shwr_length > longest_shower)
				{
					shwr_length_v.clear();
					shwr_length_v.push_back(shwr_length);
					longest_shower = shwr_length;
				}
			}
		}
		std::string tpco_id = tpco_classifier_v->at(i).first;

		for(const double shwr_length : shwr_length_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_res")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_dis")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_coh")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_mec")     {h_shwr_length_nue_cc->Fill(shwr_length); }

			if(tpco_id == "nue_bar_cc_qe")      {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_bar_cc_res")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_bar_cc_dis")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_bar_cc_coh")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_bar_cc_mec")     {h_shwr_length_nue_cc->Fill(shwr_length); }

			if(tpco_id == "nue_cc_out_fv")  {h_shwr_length_nue_cc_out_fv->Fill(shwr_length); }
			if(tpco_id == "numu_cc_qe")     {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_res")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_dis")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_coh")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_mec")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "nc")             {h_shwr_length_nc->Fill(shwr_length); }
			if(tpco_id == "nc_pi0")         {h_shwr_length_nc_pi0->Fill(shwr_length); }
			if(tpco_id == "nue_cc_mixed")   {h_shwr_length_nue_cc_mixed->Fill(shwr_length); }
			//if(tpco_id == "numu_cc_mixed")  {h_shwr_length_numu_cc_mixed->Fill(shwr_length); }
			if(tpco_id == "numu_cc_mixed")  {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "cosmic")         {h_shwr_length_cosmic->Fill(shwr_length); }
			if(tpco_id == "other_mixed")    {h_shwr_length_other_mixed->Fill(shwr_length); }
			if(tpco_id == "unmatched")      {h_shwr_length_unmatched->Fill(shwr_length); }
		}
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingShowerLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                              std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                              TH1D * h_shwr_length_nue_cc,
                                              TH1D * h_shwr_length_nue_cc_out_fv,
                                              TH1D * h_shwr_length_nue_cc_mixed,
                                              TH1D * h_shwr_length_numu_cc,
                                              TH1D * h_shwr_length_numu_cc_mixed,
                                              TH1D * h_shwr_length_nc,
                                              TH1D * h_shwr_length_nc_pi0,
                                              TH1D * h_shwr_length_cosmic,
                                              TH1D * h_shwr_length_other_mixed,
                                              TH1D * h_shwr_length_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shwr_length = leading_shower.pfpLength();

		if(tpco_id == "nue_cc_qe")      {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_res")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_dis")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_coh")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_mec")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }

		if(tpco_id == "nue_bar_cc_qe")      {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_bar_cc_res")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_bar_cc_dis")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_bar_cc_coh")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_bar_cc_mec")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }

		if(tpco_id == "nue_cc_out_fv")  {h_shwr_length_nue_cc_out_fv->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_qe")     {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_res")    {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_dis")    {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_coh")    {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_mec")    {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nc")             {h_shwr_length_nc->Fill(leading_shwr_length); }
		if(tpco_id == "nc_pi0")         {h_shwr_length_nc_pi0->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_mixed")   {h_shwr_length_nue_cc_mixed->Fill(leading_shwr_length); }
		//if(tpco_id == "numu_cc_mixed")  {h_shwr_length_numu_cc_mixed->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_mixed")  {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "cosmic")         {h_shwr_length_cosmic->Fill(leading_shwr_length); }
		if(tpco_id == "other_mixed")    {h_shwr_length_other_mixed->Fill(leading_shwr_length); }
		if(tpco_id == "unmatched")      {h_shwr_length_unmatched->Fill(leading_shwr_length); }
	}//end loop tpco
}//end function
void selection_functions::LeadingShowerLengthInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                    TH1D * h_shwr_length_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
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
		const double leading_shwr_length = leading_shower.pfpLength();
		h_shwr_length_intime->Fill(leading_shwr_length);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingShowerTrackLengths(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                    std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                    TH1D * h_shwr_trk_length_nue_cc,
                                                    TH1D * h_shwr_trk_length_nue_cc_out_fv,
                                                    TH1D * h_shwr_trk_length_nue_cc_mixed,
                                                    TH1D * h_shwr_trk_length_numu_cc,
                                                    TH1D * h_shwr_trk_length_numu_cc_mixed,
                                                    TH1D * h_shwr_trk_length_nc,
                                                    TH1D * h_shwr_trk_length_nc_pi0,
                                                    TH1D * h_shwr_trk_length_cosmic,
                                                    TH1D * h_shwr_trk_length_other_mixed,
                                                    TH1D * h_shwr_trk_length_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
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

		if(tpco_id == "nue_cc_qe")      {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_res")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_dis")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_coh")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_mec")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }

		if(tpco_id == "nue_bar_cc_qe")      {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_bar_cc_res")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_bar_cc_dis")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_bar_cc_coh")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_bar_cc_mec")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }

		if(tpco_id == "nue_cc_out_fv")  {h_shwr_trk_length_nue_cc_out_fv->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_qe")     {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_res")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_dis")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_coh")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_mec")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nc")             {h_shwr_trk_length_nc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nc_pi0")         {h_shwr_trk_length_nc_pi0->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_mixed")   {h_shwr_trk_length_nue_cc_mixed->Fill(longest_trk_leading_shwr_ratio); }
		//if(tpco_id == "numu_cc_mixed")  {h_shwr_trk_length_numu_cc_mixed->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_mixed")  {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "cosmic")         {h_shwr_trk_length_cosmic->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "other_mixed")    {h_shwr_trk_length_other_mixed->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "unmatched")      {h_shwr_trk_length_unmatched->Fill(longest_trk_leading_shwr_ratio); }
	}//end loop tpco
}//end function
void selection_functions::LeadingShowerTrackLengthsInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                          TH1D * h_shwr_trk_length_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
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
		h_shwr_trk_length_intime->Fill(longest_trk_leading_shwr_ratio);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::LongestShowerTrackLengths(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                    std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                    TH1D * h_shwr_trk_length_nue_cc,
                                                    TH1D * h_shwr_trk_length_nue_cc_out_fv,
                                                    TH1D * h_shwr_trk_length_nue_cc_mixed,
                                                    TH1D * h_shwr_trk_length_numu_cc,
                                                    TH1D * h_shwr_trk_length_numu_cc_mixed,
                                                    TH1D * h_shwr_trk_length_nc,
                                                    TH1D * h_shwr_trk_length_nc_pi0,
                                                    TH1D * h_shwr_trk_length_cosmic,
                                                    TH1D * h_shwr_trk_length_other_mixed,
                                                    TH1D * h_shwr_trk_length_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id = tpco_classifier_v->at(i).first;
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		if(n_pfp_tracks == 0) {continue; }
		std::vector< double > trk_length_v;
		std::vector< double > shwr_length_v;
		double longest_shower = 0;
		double longest_track = 0;
		for(int i = 0; i < n_pfp; i++)
		{
			auto const pfp = tpc_obj.GetParticle(i);
			const int pfp_pdg = pfp.PFParticlePdgCode();
			if(pfp_pdg == 13)
			{
				const double trk_length = pfp.pfpLength();
				if(trk_length > longest_track)
				{
					trk_length_v.clear();
					trk_length_v.push_back(trk_length);
					longest_track = trk_length;
				}
			}
			if(pfp_pdg == 11)
			{
				const double shwr_length = pfp.pfpLength();
				if(shwr_length > longest_shower)
				{
					shwr_length_v.clear();
					shwr_length_v.push_back(shwr_length);
					longest_shower = shwr_length;
				}
			}
		}
		const double longest_trk_longest_shwr_ratio = longest_track / longest_shower;

		if(tpco_id == "nue_cc_qe")      {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_res")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_dis")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_coh")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_mec")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }

		if(tpco_id == "nue_bar_cc_qe")      {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_bar_cc_res")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_bar_cc_dis")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_bar_cc_coh")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_bar_cc_mec")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }

		if(tpco_id == "nue_cc_out_fv")  {h_shwr_trk_length_nue_cc_out_fv->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_qe")     {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_res")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_dis")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_coh")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_mec")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nc")             {h_shwr_trk_length_nc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nc_pi0")         {h_shwr_trk_length_nc_pi0->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_mixed")   {h_shwr_trk_length_nue_cc_mixed->Fill(longest_trk_longest_shwr_ratio); }
		//if(tpco_id == "numu_cc_mixed")  {h_shwr_trk_length_numu_cc_mixed->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_mixed")  {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "cosmic")         {h_shwr_trk_length_cosmic->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "other_mixed")    {h_shwr_trk_length_other_mixed->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "unmatched")      {h_shwr_trk_length_unmatched->Fill(longest_trk_longest_shwr_ratio); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PlaneHitsComparisonShower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                    std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                    TH2D * h_collection_total_hits_shower_nue_cc,
                                                    TH2D * h_collection_total_hits_shower_nue_cc_out_fv,
                                                    TH2D * h_collection_total_hits_shower_nue_cc_mixed,
                                                    TH2D * h_collection_total_hits_shower_numu_cc,
                                                    TH2D * h_collection_total_hits_shower_numu_cc_mixed,
                                                    TH2D * h_collection_total_hits_shower_nc,
                                                    TH2D * h_collection_total_hits_shower_nc_pi0,
                                                    TH2D * h_collection_total_hits_shower_cosmic,
                                                    TH2D * h_collection_total_hits_shower_other_mixed,
                                                    TH2D * h_collection_total_hits_shower_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id = tpco_classifier_v->at(i).first;
		//const int n_pfp_tracks = tpc_obj.NPfpTracks();
		//if(n_pfp_tracks == 0) {continue; }
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

		if(tpco_id == "nue_cc_qe")      {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_res")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_dis")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_coh")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mec")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }

		if(tpco_id == "nue_bar_cc_qe")      {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_res")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_dis")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_coh")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_mec")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }

		if(tpco_id == "nue_cc_out_fv")  {h_collection_total_hits_shower_nue_cc_out_fv->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_qe")     {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_res")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_dis")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_coh")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mec")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc")             {h_collection_total_hits_shower_nc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc_pi0")         {h_collection_total_hits_shower_nc_pi0->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mixed")   {h_collection_total_hits_shower_nue_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		//if(tpco_id == "numu_cc_mixed")  {h_collection_total_hits_shower_numu_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mixed")  {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "cosmic")         {h_collection_total_hits_shower_cosmic->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "other_mixed")    {h_collection_total_hits_shower_other_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "unmatched")      {h_collection_total_hits_shower_unmatched->Fill(n_pfp_hits_w, n_pfp_hits); }
	}//end loop tpco
}//end function

//***************************************************************************
//***************************************************************************
void selection_functions::PlaneHitsComparisonShowerInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                          TH2D * h_collection_total_hits_shower_intime)
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
		h_collection_total_hits_shower_intime->Fill(n_pfp_hits_w, n_pfp_hits);
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PlaneHitsComparisonLeadingShower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                           std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                           TH2D * h_collection_total_hits_shower_nue_cc,
                                                           TH2D * h_collection_total_hits_shower_nue_cc_out_fv,
                                                           TH2D * h_collection_total_hits_shower_nue_cc_mixed,
                                                           TH2D * h_collection_total_hits_shower_numu_cc,
                                                           TH2D * h_collection_total_hits_shower_numu_cc_mixed,
                                                           TH2D * h_collection_total_hits_shower_nc,
                                                           TH2D * h_collection_total_hits_shower_nc_pi0,
                                                           TH2D * h_collection_total_hits_shower_cosmic,
                                                           TH2D * h_collection_total_hits_shower_other_mixed,
                                                           TH2D * h_collection_total_hits_shower_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		//const int n_pfp_tracks = tpc_obj.NPfpTracks();
		//if(n_pfp_tracks == 0) {continue; }
		const int n_pfp_hits_w = leading_shower.NumPFPHitsW();
		const int n_pfp_hits = leading_shower.NumPFPHits();

		if(tpco_id == "nue_cc_qe")      {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_res")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_dis")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_coh")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mec")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }

		if(tpco_id == "nue_bar_cc_qe")      {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_res")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_dis")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_coh")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_mec")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }

		if(tpco_id == "nue_cc_out_fv")  {h_collection_total_hits_shower_nue_cc_out_fv->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_qe")     {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_res")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_dis")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_coh")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mec")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc")             {h_collection_total_hits_shower_nc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc_pi0")         {h_collection_total_hits_shower_nc_pi0->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mixed")   {h_collection_total_hits_shower_nue_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		//if(tpco_id == "numu_cc_mixed")  {h_collection_total_hits_shower_numu_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mixed")  {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "cosmic")         {h_collection_total_hits_shower_cosmic->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "other_mixed")    {h_collection_total_hits_shower_other_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "unmatched")      {h_collection_total_hits_shower_unmatched->Fill(n_pfp_hits_w, n_pfp_hits); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PlaneHitsComparisonLeadingShowerInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                                 TH2D * h_collection_total_hits_shower_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
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
		const int n_pfp_hits_w = leading_shower.NumPFPHitsW();
		const int n_pfp_hits = leading_shower.NumPFPHits();
		h_collection_total_hits_shower_intime->Fill(n_pfp_hits_w, n_pfp_hits);
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PlaneHitsComparisonTrack(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                   std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                   TH2D * h_collection_total_hits_track_nue_cc,
                                                   TH2D * h_collection_total_hits_track_nue_cc_out_fv,
                                                   TH2D * h_collection_total_hits_track_nue_cc_mixed,
                                                   TH2D * h_collection_total_hits_track_numu_cc,
                                                   TH2D * h_collection_total_hits_track_numu_cc_mixed,
                                                   TH2D * h_collection_total_hits_track_nc,
                                                   TH2D * h_collection_total_hits_track_nc_pi0,
                                                   TH2D * h_collection_total_hits_track_cosmic,
                                                   TH2D * h_collection_total_hits_track_other_mixed,
                                                   TH2D * h_collection_total_hits_track_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id = tpco_classifier_v->at(i).first;
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

		if(tpco_id == "nue_cc_qe")      {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_res")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_dis")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_coh")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mec")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }

		if(tpco_id == "nue_bar_cc_qe")      {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_res")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_dis")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_coh")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_bar_cc_mec")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }

		if(tpco_id == "nue_cc_out_fv")  {h_collection_total_hits_track_nue_cc_out_fv->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_qe")     {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_res")    {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_dis")    {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_coh")    {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mec")    {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc")             {h_collection_total_hits_track_nc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc_pi0")         {h_collection_total_hits_track_nc_pi0->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mixed")   {h_collection_total_hits_track_nue_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		//if(tpco_id == "numu_cc_mixed")  {h_collection_total_hits_track_numu_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mixed")  {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "cosmic")         {h_collection_total_hits_track_cosmic->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "other_mixed")    {h_collection_total_hits_track_other_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "unmatched")      {h_collection_total_hits_track_unmatched->Fill(n_pfp_hits_w, n_pfp_hits); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PlaneHitsComparisonTrackInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                         TH2D * h_collection_total_hits_track_intime)
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
		h_collection_total_hits_track_intime->Fill(n_pfp_hits_w, n_pfp_hits);
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::HitsPlots1D(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                      std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                      TH1D * h_collection_hits_track_nue_cc,
                                      TH1D * h_collection_hits_track_nue_cc_out_fv,
                                      TH1D * h_collection_hits_track_nue_cc_mixed,
                                      TH1D * h_collection_hits_track_numu_cc,
                                      TH1D * h_collection_hits_track_numu_cc_mixed,
                                      TH1D * h_collection_hits_track_nc,
                                      TH1D * h_collection_hits_track_nc_pi0,
                                      TH1D * h_collection_hits_track_cosmic,
                                      TH1D * h_collection_hits_track_other_mixed,
                                      TH1D * h_collection_hits_track_unmatched,
                                      TH1D * h_collection_hits_shower_nue_cc,
                                      TH1D * h_collection_hits_shower_nue_cc_out_fv,
                                      TH1D * h_collection_hits_shower_nue_cc_mixed,
                                      TH1D * h_collection_hits_shower_numu_cc,
                                      TH1D * h_collection_hits_shower_numu_cc_mixed,
                                      TH1D * h_collection_hits_shower_nc,
                                      TH1D * h_collection_hits_shower_nc_pi0,
                                      TH1D * h_collection_hits_shower_cosmic,
                                      TH1D * h_collection_hits_shower_other_mixed,
                                      TH1D * h_collection_hits_shower_unmatched,
                                      TH1D * h_collection_hits_leading_shower_nue_cc,
                                      TH1D * h_collection_hits_leading_shower_nue_cc_out_fv,
                                      TH1D * h_collection_hits_leading_shower_nue_cc_mixed,
                                      TH1D * h_collection_hits_leading_shower_numu_cc,
                                      TH1D * h_collection_hits_leading_shower_numu_cc_mixed,
                                      TH1D * h_collection_hits_leading_shower_nc,
                                      TH1D * h_collection_hits_leading_shower_nc_pi0,
                                      TH1D * h_collection_hits_leading_shower_cosmic,
                                      TH1D * h_collection_hits_leading_shower_other_mixed,
                                      TH1D * h_collection_hits_leading_shower_unmatched,
                                      TH1D * h_total_hits_leading_shower_nue_cc,
                                      TH1D * h_total_hits_leading_shower_nue_cc_out_fv,
                                      TH1D * h_total_hits_leading_shower_nue_cc_mixed,
                                      TH1D * h_total_hits_leading_shower_numu_cc,
                                      TH1D * h_total_hits_leading_shower_numu_cc_mixed,
                                      TH1D * h_total_hits_leading_shower_nc,
                                      TH1D * h_total_hits_leading_shower_nc_pi0,
                                      TH1D * h_total_hits_leading_shower_cosmic,
                                      TH1D * h_total_hits_leading_shower_other_mixed,
                                      TH1D * h_total_hits_leading_shower_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
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

		if(n_pfp_tracks != 0) {
			for(auto const n_pfp_hits_w_track : n_pfp_hits_w_track_v)
			{
				if(tpco_id == "nue_cc_qe")      {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_res")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_dis")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_coh")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_mec")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }

				if(tpco_id == "nue_bar_cc_qe")      {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_bar_cc_res")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_bar_cc_dis")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_bar_cc_coh")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_bar_cc_mec")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }

				if(tpco_id == "nue_cc_out_fv")  {h_collection_hits_track_nue_cc_out_fv->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_qe")     {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_res")    {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_dis")    {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_coh")    {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_mec")    {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nc")             {h_collection_hits_track_nc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nc_pi0")         {h_collection_hits_track_nc_pi0->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_mixed")   {h_collection_hits_track_nue_cc_mixed->Fill(n_pfp_hits_w_track); }
				//if(tpco_id == "numu_cc_mixed")  {h_collection_hits_track_numu_cc_mixed->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_mixed")  {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "cosmic")         {h_collection_hits_track_cosmic->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "other_mixed")    {h_collection_hits_track_other_mixed->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "unmatched")      {h_collection_hits_track_unmatched->Fill(n_pfp_hits_w_track); }
			}
		}

		for(auto const n_pfp_hits_w_shower : n_pfp_hits_w_shower_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_res")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_dis")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_coh")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_mec")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }

			if(tpco_id == "nue_bar_cc_qe")      {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_bar_cc_res")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_bar_cc_dis")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_bar_cc_coh")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_bar_cc_mec")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }

			if(tpco_id == "nue_cc_out_fv")  {h_collection_hits_shower_nue_cc_out_fv->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_qe")     {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_res")    {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_dis")    {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_coh")    {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_mec")    {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nc")             {h_collection_hits_shower_nc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nc_pi0")         {h_collection_hits_shower_nc_pi0->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_mixed")   {h_collection_hits_shower_nue_cc_mixed->Fill(n_pfp_hits_w_shower); }
			//if(tpco_id == "numu_cc_mixed")  {h_collection_hits_shower_numu_cc_mixed->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_mixed")  {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "cosmic")         {h_collection_hits_shower_cosmic->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "other_mixed")    {h_collection_hits_shower_other_mixed->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "unmatched")      {h_collection_hits_shower_unmatched->Fill(n_pfp_hits_w_shower); }
		}

		if(tpco_id == "nue_cc_qe")      {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_res")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_dis")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_coh")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_mec")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }

		if(tpco_id == "nue_bar_cc_qe")      {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_bar_cc_res")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_bar_cc_dis")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_bar_cc_coh")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_bar_cc_mec")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }

		if(tpco_id == "nue_cc_out_fv")  {h_collection_hits_leading_shower_nue_cc_out_fv->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_qe")     {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_res")    {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_dis")    {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_coh")    {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_mec")    {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nc")             {h_collection_hits_leading_shower_nc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nc_pi0")         {h_collection_hits_leading_shower_nc_pi0->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_mixed")   {h_collection_hits_leading_shower_nue_cc_mixed->Fill(n_pfp_hits_w_leading_shower); }
		//if(tpco_id == "numu_cc_mixed")  {h_collection_hits_leading_shower_numu_cc_mixed->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_mixed")  {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "cosmic")         {h_collection_hits_leading_shower_cosmic->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "other_mixed")    {h_collection_hits_leading_shower_other_mixed->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "unmatched")      {h_collection_hits_leading_shower_unmatched->Fill(n_pfp_hits_w_leading_shower); }

		if(tpco_id == "nue_cc_qe")      {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_res")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_dis")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_coh")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_mec")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }

		if(tpco_id == "nue_bar_cc_qe")      {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_bar_cc_res")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_bar_cc_dis")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_bar_cc_coh")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_bar_cc_mec")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }

		if(tpco_id == "nue_cc_out_fv")  {h_total_hits_leading_shower_nue_cc_out_fv->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_qe")     {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_res")    {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_dis")    {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_coh")    {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_mec")    {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nc")             {h_total_hits_leading_shower_nc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nc_pi0")         {h_total_hits_leading_shower_nc_pi0->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_mixed")   {h_total_hits_leading_shower_nue_cc_mixed->Fill(n_pfp_hits_leading_shower); }
		//if(tpco_id == "numu_cc_mixed")  {h_total_hits_leading_shower_numu_cc_mixed->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_mixed")  {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "cosmic")         {h_total_hits_leading_shower_cosmic->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "other_mixed")    {h_total_hits_leading_shower_other_mixed->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "unmatched")      {h_total_hits_leading_shower_unmatched->Fill(n_pfp_hits_leading_shower); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::HitsPlots1DInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                            std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                            TH1D * h_collection_hits_track_intime,
                                            TH1D * h_collection_hits_shower_intime,
                                            TH1D * h_collection_hits_leading_shower_intime,
                                            TH1D * h_total_hits_leading_shower_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
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

		if(n_pfp_tracks != 0) {for(auto const n_pfp_hits_w_track : n_pfp_hits_w_track_v) {h_collection_hits_track_intime->Fill(n_pfp_hits_w_track); }}
		for(auto const n_pfp_hits_w_shower : n_pfp_hits_w_shower_v) {h_collection_hits_shower_intime->Fill(n_pfp_hits_w_shower); }
		h_collection_hits_leading_shower_intime->Fill(n_pfp_hits_w_leading_shower);
		h_total_hits_leading_shower_intime->Fill(n_pfp_hits_leading_shower);
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
int selection_functions::MapFailureCutToString(const std::string failure_cut)
{
	int failure_reason = -1;
	//if(failure_cut == "Passed")          {return failure_reason;   }
	if(failure_cut == "HasNue")          {failure_reason = 0.0;    }
	if(failure_cut == "InFV")            {failure_reason = 1.0;    }
	if(failure_cut == "FlashDist")       {failure_reason = 2.0;    }
	if(failure_cut == "ShwrVtx")         {failure_reason = 3.0;    }
	if(failure_cut == "TrkVtx")          {failure_reason = 4.0;    }
	if(failure_cut == "HitThreshold")    {failure_reason = 5.0;    }
	if(failure_cut == "OpenAngle")       {failure_reason = 6.0;    }
	if(failure_cut == "dEdX")            {failure_reason = 7.0;    }
	if(failure_cut == "SecondaryDist")   {failure_reason = 8.0;    }
	if(failure_cut == "HitLengthRatio")  {failure_reason = 9.0;    }
	if(failure_cut == "HitThresholdW")   {failure_reason = 10.0;   }
	if(failure_cut == "TrkShwrLenRatio") {failure_reason = 11.0;   }
	return failure_reason;
}
//***************************************************************************
//***************************************************************************
void selection_functions::EnergyHits(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                     std::vector<std::pair<std::string, int> > * tpco_classifier_v, double mc_nu_energy, double mc_ele_energy,
                                     TH2D * h_ele_eng_total_hits, TH2D * h_ele_eng_colleciton_hits, TH2D * h_nu_eng_total_hits, TH2D * h_nu_eng_collection_hits)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int n_pfp_hits_w  = leading_shower.NumPFPHitsW();
		const int n_pfp_hits    = leading_shower.NumPFPHits();
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_mec" ||
		   tpco_id == "nue_bar_cc_qe" || tpco_id == "nue_bar_cc_res" || tpco_id == "nue_bar_cc_coh" || tpco_id == "nue_bar_cc_dis" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_eng_total_hits->Fill(n_pfp_hits, mc_ele_energy);
			h_ele_eng_colleciton_hits->Fill(n_pfp_hits_w, mc_ele_energy);
			h_nu_eng_total_hits->Fill(n_pfp_hits, mc_nu_energy);
			h_nu_eng_collection_hits->Fill(n_pfp_hits_w, mc_nu_energy);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::FailureReason(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                        std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                        TH1D * h_failure_reason_nue_cc,
                                        TH1D * h_failure_reason_nue_cc_out_fv,
                                        TH1D * h_failure_reason_nue_cc_mixed,
                                        TH1D * h_failure_reason_numu_cc,
                                        TH1D * h_failure_reason_numu_cc_mixed,
                                        TH1D * h_failure_reason_nc,
                                        TH1D * h_failure_reason_nc_pi0,
                                        TH1D * h_failure_reason_cosmic,
                                        TH1D * h_failure_reason_other_mixed,
                                        TH1D * h_failure_reason_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		//if(passed_tpco->at(i).first == 0) {continue; }
		const std::string failure_cut = passed_tpco->at(i).second;
		const int failure_reason = MapFailureCutToString(failure_cut);
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id = tpco_classifier_v->at(i).first;

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_failure_reason_nue_cc->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_failure_reason_nue_cc_out_fv->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_failure_reason_nue_cc->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_failure_reason_nue_cc->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_failure_reason_nue_cc->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_failure_reason_nue_cc->Fill(failure_reason);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_failure_reason_numu_cc->Fill(failure_reason);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_failure_reason_numu_cc->Fill(failure_reason);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_failure_reason_numu_cc->Fill(failure_reason);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_failure_reason_numu_cc->Fill(failure_reason);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_failure_reason_numu_cc->Fill(failure_reason);
		}
		if(tpco_id == "nc")
		{
			h_failure_reason_nc->Fill(failure_reason);
		}
		if(tpco_id == "nc_pi0")
		{
			h_failure_reason_nc_pi0->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_failure_reason_nue_cc_mixed->Fill(failure_reason);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_failure_reason_numu_cc_mixed->Fill(failure_reason);
			h_failure_reason_numu_cc->Fill(failure_reason);
		}
		if(tpco_id == "cosmic")
		{
			h_failure_reason_cosmic->Fill(failure_reason);
		}
		if(tpco_id == "other_mixed")
		{
			h_failure_reason_other_mixed->Fill(failure_reason);
		}
		if(tpco_id == "unmatched")
		{
			h_failure_reason_unmatched->Fill(failure_reason);
		}
	}//end loop tpco
}
//leading shower cos theta - after selection cuts
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingCosTheta(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco,
                                          const double theta_translation, const double phi_translation, bool _verbose,
                                          std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                          TH1D * h_ele_cos_theta_nue_cc,
                                          TH1D * h_ele_cos_theta_nue_cc_out_fv,
                                          TH1D * h_ele_cos_theta_nue_cc_mixed,
                                          TH1D * h_ele_cos_theta_numu_cc,
                                          TH1D * h_ele_cos_theta_numu_cc_mixed,
                                          TH1D * h_ele_cos_theta_nc,
                                          TH1D * h_ele_cos_theta_nc_pi0,
                                          TH1D * h_ele_cos_theta_cosmic,
                                          TH1D * h_ele_cos_theta_other_mixed,
                                          TH1D * h_ele_cos_theta_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		double leading_shower_cos_theta = -99;
		if(theta_translation == 0 && phi_translation == 0) {leading_shower_cos_theta = leading_shower_z; }

		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, theta_translation, phi_translation);
		if(theta_translation != 0 && phi_translation != 0)
		{
			leading_shower_cos_theta = shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag());
		}

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_ele_cos_theta_nue_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_ele_cos_theta_nue_cc_out_fv->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_ele_cos_theta_nue_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_ele_cos_theta_nue_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_ele_cos_theta_nue_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_cos_theta_nue_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_ele_cos_theta_numu_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_ele_cos_theta_numu_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_ele_cos_theta_numu_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_ele_cos_theta_numu_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_ele_cos_theta_numu_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "nc")
		{
			h_ele_cos_theta_nc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "nc_pi0")
		{
			h_ele_cos_theta_nc_pi0->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_ele_cos_theta_nue_cc_mixed->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_cos_theta_numu_cc_mixed->Fill(leading_shower_cos_theta);
			h_ele_cos_theta_numu_cc->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "cosmic")
		{
			h_ele_cos_theta_cosmic->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "other_mixed")
		{
			h_ele_cos_theta_other_mixed->Fill(leading_shower_cos_theta);
		}
		if(tpco_id == "unmatched")
		{
			h_ele_cos_theta_unmatched->Fill(leading_shower_cos_theta);
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingCosThetaInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                std::vector<std::pair<int, std::string> > * passed_tpco,
                                                const double theta_translation, const double phi_translation, bool _verbose, TH1D * h_ele_cos_theta_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		int most_hits = 0;
		int leading_index = 0;
		const int n_pfp = tpc_obj.NumPFParticles();
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
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, theta_translation, phi_translation);
		double leading_shower_cos_theta = shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag());
		if(theta_translation == 0 && phi_translation == 0) {leading_shower_cos_theta = leading_shower_z; }
		h_ele_cos_theta_intime->Fill(leading_shower_cos_theta);
	}//end pfp loop
}
//leading shower cos theta - before most selection cuts
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingMomentum(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                          std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                          TH1D * h_ele_pfp_momentum_nue_cc,
                                          TH1D * h_ele_pfp_momentum_nue_cc_out_fv,
                                          TH1D * h_ele_pfp_momentum_nue_cc_mixed,
                                          TH1D * h_ele_pfp_momentum_numu_cc,
                                          TH1D * h_ele_pfp_momentum_numu_cc_mixed,
                                          TH1D * h_ele_pfp_momentum_nc,
                                          TH1D * h_ele_pfp_momentum_nc_pi0,
                                          TH1D * h_ele_pfp_momentum_cosmic,
                                          TH1D * h_ele_pfp_momentum_other_mixed,
                                          TH1D * h_ele_pfp_momentum_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_momentum = leading_shower.pfpMomentum();
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_ele_pfp_momentum_nue_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_ele_pfp_momentum_nue_cc_out_fv->Fill(leading_shower_momentum);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_ele_pfp_momentum_nue_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_ele_pfp_momentum_nue_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_ele_pfp_momentum_nue_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_pfp_momentum_nue_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_ele_pfp_momentum_numu_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_ele_pfp_momentum_numu_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_ele_pfp_momentum_numu_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_ele_pfp_momentum_numu_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_ele_pfp_momentum_numu_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "nc")
		{
			h_ele_pfp_momentum_nc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "nc_pi0")
		{
			h_ele_pfp_momentum_nc_pi0->Fill(leading_shower_momentum);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_ele_pfp_momentum_nue_cc_mixed->Fill(leading_shower_momentum);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_pfp_momentum_numu_cc_mixed->Fill(leading_shower_momentum);
			h_ele_pfp_momentum_numu_cc->Fill(leading_shower_momentum);
		}
		if(tpco_id == "cosmic")
		{
			h_ele_pfp_momentum_cosmic->Fill(leading_shower_momentum);
		}
		if(tpco_id == "other_mixed")
		{
			h_ele_pfp_momentum_other_mixed->Fill(leading_shower_momentum);
		}
		if(tpco_id == "unmatched")
		{
			h_ele_pfp_momentum_unmatched->Fill(leading_shower_momentum);
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingMomentumInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, TH1D * h_ele_pfp_momentum_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		int most_hits = 0;
		int leading_index = 0;
		const int n_pfp = tpc_obj.NumPFParticles();
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
		//const double leading_shower_cos_theta = leading_shower.pfpDirZ() / leading_shower.pfpMomentum();
		const double leading_shower_momentum = leading_shower.pfpMomentum();
		//std::cout << leading_shower.pfpDirZ() << " , " << leading_shower.pfpMomentum() << ", " << leading_shower_cos_theta << std::endl;
		h_ele_pfp_momentum_intime->Fill(leading_shower_momentum);
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingMomentumTrackTopology(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                       std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                       TH1D * h_ele_pfp_momentum_no_track_nue_cc,
                                                       TH1D * h_ele_pfp_momentum_no_track_nue_cc_out_fv,
                                                       TH1D * h_ele_pfp_momentum_no_track_nue_cc_mixed,
                                                       TH1D * h_ele_pfp_momentum_no_track_numu_cc,
                                                       TH1D * h_ele_pfp_momentum_no_track_numu_cc_mixed,
                                                       TH1D * h_ele_pfp_momentum_no_track_nc,
                                                       TH1D * h_ele_pfp_momentum_no_track_nc_pi0,
                                                       TH1D * h_ele_pfp_momentum_no_track_cosmic,
                                                       TH1D * h_ele_pfp_momentum_no_track_other_mixed,
                                                       TH1D * h_ele_pfp_momentum_no_track_unmatched,
                                                       TH1D * h_ele_pfp_momentum_has_track_nue_cc,
                                                       TH1D * h_ele_pfp_momentum_has_track_nue_cc_out_fv,
                                                       TH1D * h_ele_pfp_momentum_has_track_nue_cc_mixed,
                                                       TH1D * h_ele_pfp_momentum_has_track_numu_cc,
                                                       TH1D * h_ele_pfp_momentum_has_track_numu_cc_mixed,
                                                       TH1D * h_ele_pfp_momentum_has_track_nc,
                                                       TH1D * h_ele_pfp_momentum_has_track_nc_pi0,
                                                       TH1D * h_ele_pfp_momentum_has_track_cosmic,
                                                       TH1D * h_ele_pfp_momentum_has_track_other_mixed,
                                                       TH1D * h_ele_pfp_momentum_has_track_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_momentum = leading_shower.pfpMomentum();
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_nue_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_nue_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_nue_cc_out_fv->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_nue_cc_out_fv->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_has_track_nue_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_nue_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_nue_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_nue_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_nue_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_nue_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_nue_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_nue_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "numu_cc_qe")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_numu_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_numu_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "numu_cc_res")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_numu_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_numu_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "numu_cc_dis")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_numu_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_numu_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "numu_cc_coh")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_numu_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_numu_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "numu_cc_mec")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_numu_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_numu_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "nc")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_nc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_nc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "nc_pi0")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_nc_pi0->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_nc_pi0->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "nue_cc_mixed")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_nue_cc_mixed->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_nue_cc_mixed->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "numu_cc_mixed")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_numu_cc->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_numu_cc->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "cosmic")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_cosmic->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_cosmic->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "other_mixed")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_other_mixed->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_other_mixed->Fill(leading_shower_momentum);
			}
		}
		if(tpco_id == "unmatched")
		{
			if(n_pfp_tracks == 0)
			{
				h_ele_pfp_momentum_no_track_unmatched->Fill(leading_shower_momentum);
			}
			if(n_pfp_tracks >= 1)
			{
				h_ele_pfp_momentum_has_track_unmatched->Fill(leading_shower_momentum);
			}
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingMomentumTrackTopologyInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                             TH1D * h_ele_pfp_momentum_no_track_intime,
                                                             TH1D * h_ele_pfp_momentum_has_track_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		int most_hits = 0;
		int leading_index = 0;
		const int n_pfp = tpc_obj.NumPFParticles();
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
		//const double leading_shower_cos_theta = leading_shower.pfpDirZ() / leading_shower.pfpMomentum();
		const double leading_shower_momentum = leading_shower.pfpMomentum();
		if(n_pfp_tracks == 0)
		{
			h_ele_pfp_momentum_no_track_intime->Fill(leading_shower_momentum);
		}
		if(n_pfp_tracks >= 1)
		{
			h_ele_pfp_momentum_has_track_intime->Fill(leading_shower_momentum);
		}
	}//end pfp loop
}

//leading shower cos theta
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingTheta(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                       std::vector<std::pair<int, std::string> > * passed_tpco,
                                       const double theta_translation, const double phi_translation, bool _verbose,
                                       std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                       TH1D * h_ele_pfp_theta_nue_cc,
                                       TH1D * h_ele_pfp_theta_nue_cc_out_fv,
                                       TH1D * h_ele_pfp_theta_nue_cc_mixed,
                                       TH1D * h_ele_pfp_theta_numu_cc,
                                       TH1D * h_ele_pfp_theta_numu_cc_mixed,
                                       TH1D * h_ele_pfp_theta_nc,
                                       TH1D * h_ele_pfp_theta_nc_pi0,
                                       TH1D * h_ele_pfp_theta_cosmic,
                                       TH1D * h_ele_pfp_theta_other_mixed,
                                       TH1D * h_ele_pfp_theta_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		//const double leading_shower_theta = acos(leading_shower.pfpDirZ()) * (180 / 3.1415);

		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, theta_translation, phi_translation);
		const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_ele_pfp_theta_nue_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_ele_pfp_theta_nue_cc_out_fv->Fill(leading_shower_theta);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_ele_pfp_theta_nue_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_ele_pfp_theta_nue_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_ele_pfp_theta_nue_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_pfp_theta_nue_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_ele_pfp_theta_numu_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_ele_pfp_theta_numu_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_ele_pfp_theta_numu_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_ele_pfp_theta_numu_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_ele_pfp_theta_numu_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "nc")
		{
			h_ele_pfp_theta_nc->Fill(leading_shower_theta);
		}
		if(tpco_id == "nc_pi0")
		{
			h_ele_pfp_theta_nc_pi0->Fill(leading_shower_theta);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_ele_pfp_theta_nue_cc_mixed->Fill(leading_shower_theta);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_pfp_theta_numu_cc_mixed->Fill(leading_shower_theta);
			h_ele_pfp_theta_numu_cc->Fill(leading_shower_theta);
		}
		if(tpco_id == "cosmic")
		{
			h_ele_pfp_theta_cosmic->Fill(leading_shower_theta);
		}
		if(tpco_id == "other_mixed")
		{
			h_ele_pfp_theta_other_mixed->Fill(leading_shower_theta);
		}
		if(tpco_id == "unmatched")
		{
			h_ele_pfp_theta_unmatched->Fill(leading_shower_theta);
		}
	}//end pfp loop
}
//leading shower theta
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingThetaInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco,
                                             const double theta_translation, const double phi_translation, bool _verbose,
                                             TH1D * h_ele_pfp_theta_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		int most_hits = 0;
		int leading_index = 0;
		const int n_pfp = tpc_obj.NumPFParticles();
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
		//const double leading_shower_theta = acos(leading_shower.pfpDirZ()) * (180 / 3.1415);

		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, theta_translation, phi_translation);
		const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);

		h_ele_pfp_theta_intime->Fill(leading_shower_theta);
	}//end pfp loop
}
//leading shower theta
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingPhi(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                     std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                     TH1D * h_ele_pfp_phi_nue_cc,
                                     TH1D * h_ele_pfp_phi_nue_cc_out_fv,
                                     TH1D * h_ele_pfp_phi_nue_cc_mixed,
                                     TH1D * h_ele_pfp_phi_numu_cc,
                                     TH1D * h_ele_pfp_phi_numu_cc_mixed,
                                     TH1D * h_ele_pfp_phi_nc,
                                     TH1D * h_ele_pfp_phi_nc_pi0,
                                     TH1D * h_ele_pfp_phi_cosmic,
                                     TH1D * h_ele_pfp_phi_other_mixed,
                                     TH1D * h_ele_pfp_phi_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_phi = atan2(leading_shower.pfpDirY(), leading_shower.pfpDirX()) * 180 / 3.1415;
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_ele_pfp_phi_nue_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_ele_pfp_phi_nue_cc_out_fv->Fill(leading_shower_phi);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_ele_pfp_phi_nue_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_ele_pfp_phi_nue_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_ele_pfp_phi_nue_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_pfp_phi_nue_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_ele_pfp_phi_numu_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_ele_pfp_phi_numu_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_ele_pfp_phi_numu_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_ele_pfp_phi_numu_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_ele_pfp_phi_numu_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "nc")
		{
			h_ele_pfp_phi_nc->Fill(leading_shower_phi);
		}
		if(tpco_id == "nc_pi0")
		{
			h_ele_pfp_phi_nc_pi0->Fill(leading_shower_phi);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_ele_pfp_phi_nue_cc_mixed->Fill(leading_shower_phi);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_pfp_phi_numu_cc_mixed->Fill(leading_shower_phi);
			h_ele_pfp_phi_numu_cc->Fill(leading_shower_phi);
		}
		if(tpco_id == "cosmic")
		{
			h_ele_pfp_phi_cosmic->Fill(leading_shower_phi);
		}
		if(tpco_id == "other_mixed")
		{
			h_ele_pfp_phi_other_mixed->Fill(leading_shower_phi);
		}
		if(tpco_id == "unmatched")
		{
			h_ele_pfp_phi_unmatched->Fill(leading_shower_phi);
		}
	}//end pfp loop
}
//leading shower phi
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingPhiInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                           TH1D * h_ele_pfp_phi_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		int most_hits = 0;
		int leading_index = 0;
		const int n_pfp = tpc_obj.NumPFParticles();
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
		const double leading_shower_phi = atan2(leading_shower.pfpDirY(), leading_shower.pfpDirX()) * 180 / 3.1415;
		h_ele_pfp_phi_intime->Fill(leading_shower_phi);
	}//end pfp loop
}
//leading shower phi
//***************************************************************************
//***************************************************************************
void selection_functions::Leading1Shwr2Shwr(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                            std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                            std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                            TH1D * h_leading_shwr_length_1shwr_nue_cc,
                                            TH1D * h_leading_shwr_length_1shwr_nue_cc_out_fv,
                                            TH1D * h_leading_shwr_length_1shwr_nue_cc_mixed,
                                            TH1D * h_leading_shwr_length_1shwr_numu_cc,
                                            TH1D * h_leading_shwr_length_1shwr_numu_cc_mixed,
                                            TH1D * h_leading_shwr_length_1shwr_nc,
                                            TH1D * h_leading_shwr_length_1shwr_nc_pi0,
                                            TH1D * h_leading_shwr_length_1shwr_cosmic,
                                            TH1D * h_leading_shwr_length_1shwr_other_mixed,
                                            TH1D * h_leading_shwr_length_1shwr_unmatched,
                                            TH1D * h_leading_shwr_length_2shwr_nue_cc,
                                            TH1D * h_leading_shwr_length_2shwr_nue_cc_out_fv,
                                            TH1D * h_leading_shwr_length_2shwr_nue_cc_mixed,
                                            TH1D * h_leading_shwr_length_2shwr_numu_cc,
                                            TH1D * h_leading_shwr_length_2shwr_numu_cc_mixed,
                                            TH1D * h_leading_shwr_length_2shwr_nc,
                                            TH1D * h_leading_shwr_length_2shwr_nc_pi0,
                                            TH1D * h_leading_shwr_length_2shwr_cosmic,
                                            TH1D * h_leading_shwr_length_2shwr_other_mixed,
                                            TH1D * h_leading_shwr_length_2shwr_unmatched,
                                            TH1D * h_leading_shwr_hits_1shwr_nue_cc,
                                            TH1D * h_leading_shwr_hits_1shwr_nue_cc_out_fv,
                                            TH1D * h_leading_shwr_hits_1shwr_nue_cc_mixed,
                                            TH1D * h_leading_shwr_hits_1shwr_numu_cc,
                                            TH1D * h_leading_shwr_hits_1shwr_numu_cc_mixed,
                                            TH1D * h_leading_shwr_hits_1shwr_nc,
                                            TH1D * h_leading_shwr_hits_1shwr_nc_pi0,
                                            TH1D * h_leading_shwr_hits_1shwr_cosmic,
                                            TH1D * h_leading_shwr_hits_1shwr_other_mixed,
                                            TH1D * h_leading_shwr_hits_1shwr_unmatched,
                                            TH1D * h_leading_shwr_hits_2shwr_nue_cc,
                                            TH1D * h_leading_shwr_hits_2shwr_nue_cc_out_fv,
                                            TH1D * h_leading_shwr_hits_2shwr_nue_cc_mixed,
                                            TH1D * h_leading_shwr_hits_2shwr_numu_cc,
                                            TH1D * h_leading_shwr_hits_2shwr_numu_cc_mixed,
                                            TH1D * h_leading_shwr_hits_2shwr_nc,
                                            TH1D * h_leading_shwr_hits_2shwr_nc_pi0,
                                            TH1D * h_leading_shwr_hits_2shwr_cosmic,
                                            TH1D * h_leading_shwr_hits_2shwr_other_mixed,
                                            TH1D * h_leading_shwr_hits_2shwr_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		//const double leading_shower_phi = atan2(leading_shower.pfpDirY(), leading_shower.pfpDirX()) * 180 / 3.1415;
		const double leading_shower_length = leading_shower.pfpLength();
		const double leading_shower_hits = leading_shower.NumPFPHits();
		if(n_pfp_showers == 1)
		{
			if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
			{
				h_leading_shwr_length_1shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_out_fv")
			{
				h_leading_shwr_length_1shwr_nue_cc_out_fv->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_nue_cc_out_fv->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
			{
				h_leading_shwr_length_1shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
			{
				h_leading_shwr_length_1shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
			{
				h_leading_shwr_length_1shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
			{
				h_leading_shwr_length_1shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_qe")
			{
				h_leading_shwr_length_1shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_res")
			{
				h_leading_shwr_length_1shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_dis")
			{
				h_leading_shwr_length_1shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_coh")
			{
				h_leading_shwr_length_1shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_mec")
			{
				h_leading_shwr_length_1shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nc")
			{
				h_leading_shwr_length_1shwr_nc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_nc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nc_pi0")
			{
				h_leading_shwr_length_1shwr_nc_pi0->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_nc_pi0->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_mixed")
			{
				h_leading_shwr_length_1shwr_nue_cc_mixed->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_nue_cc_mixed->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_mixed")
			{
				// h_leading_shwr_length_1shwr_numu_cc_mixed->Fill(leading_shower_length);
				// h_leading_shwr_hits_1shwr_numu_cc_mixed->Fill(leading_shower_hits);
				h_leading_shwr_length_1shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "cosmic")
			{
				h_leading_shwr_length_1shwr_cosmic->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_cosmic->Fill(leading_shower_hits);
			}
			if(tpco_id == "other_mixed")
			{
				h_leading_shwr_length_1shwr_other_mixed->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_other_mixed->Fill(leading_shower_hits);
			}
			if(tpco_id == "unmatched")
			{
				h_leading_shwr_length_1shwr_unmatched->Fill(leading_shower_length);
				h_leading_shwr_hits_1shwr_unmatched->Fill(leading_shower_hits);
			}
		}
		if(n_pfp_showers >= 2)
		{
			if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
			{
				h_leading_shwr_length_2shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_out_fv")
			{
				h_leading_shwr_length_2shwr_nue_cc_out_fv->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_nue_cc_out_fv->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
			{
				h_leading_shwr_length_2shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
			{
				h_leading_shwr_length_2shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
			{
				h_leading_shwr_length_2shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
			{
				h_leading_shwr_length_2shwr_nue_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_nue_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_qe")
			{
				h_leading_shwr_length_2shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_res")
			{
				h_leading_shwr_length_2shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_dis")
			{
				h_leading_shwr_length_2shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_coh")
			{
				h_leading_shwr_length_2shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_mec")
			{
				h_leading_shwr_length_2shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nc")
			{
				h_leading_shwr_length_2shwr_nc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_nc->Fill(leading_shower_hits);
			}
			if(tpco_id == "nc_pi0")
			{
				h_leading_shwr_length_2shwr_nc_pi0->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_nc_pi0->Fill(leading_shower_hits);
			}
			if(tpco_id == "nue_cc_mixed")
			{
				h_leading_shwr_length_2shwr_nue_cc_mixed->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_nue_cc_mixed->Fill(leading_shower_hits);
			}
			if(tpco_id == "numu_cc_mixed")
			{
				// h_leading_shwr_length_2shwr_numu_cc_mixed->Fill(leading_shower_length);
				// h_leading_shwr_hits_2shwr_numu_cc_mixed->Fill(leading_shower_hits);
				h_leading_shwr_length_2shwr_numu_cc->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_numu_cc->Fill(leading_shower_hits);
			}
			if(tpco_id == "cosmic")
			{
				h_leading_shwr_length_2shwr_cosmic->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_cosmic->Fill(leading_shower_hits);
			}
			if(tpco_id == "other_mixed")
			{
				h_leading_shwr_length_2shwr_other_mixed->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_other_mixed->Fill(leading_shower_hits);
			}
			if(tpco_id == "unmatched")
			{
				h_leading_shwr_length_2shwr_unmatched->Fill(leading_shower_length);
				h_leading_shwr_hits_2shwr_unmatched->Fill(leading_shower_hits);
			}
		}
	}//end pfp loop
}
//leading shower phi
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutVector2DPlots(std::vector<std::tuple<int, int, int, double, double, double,
                                                                      std::string, std::string, int, int, double> > * post_cuts_v,
                                               bool _post_cuts_verbose, const double intime_scale_factor, const double data_scale_factor,
                                               TH2 * post_cuts_num_tracks_showers_purity_qe,
                                               TH2 * post_cuts_num_tracks_showers_purity_res,
                                               TH2 * post_cuts_num_tracks_showers_purity_dis,
                                               TH2 * post_cuts_num_tracks_showers_purity_coh,
                                               TH2 * post_cuts_num_tracks_showers_purity_mec,
                                               TH2 * post_cuts_num_tracks_showers_purity_total,
                                               TH2 * post_cuts_num_tracks_showers_signal_total,
                                               TH2 * post_cuts_num_tracks_showers_bkg_total,
                                               TH2 * post_cuts_num_tracks_showers_total_total)
{
	int signal_events_1_1_qe = 0;
	int signal_events_2_1_qe = 0;
	int signal_events_3_1_qe = 0;
	int signal_events_4_1_qe = 0;

	int signal_events_1_1_res = 0;
	int signal_events_2_1_res = 0;
	int signal_events_3_1_res = 0;
	int signal_events_4_1_res = 0;

	int signal_events_1_1_dis = 0;
	int signal_events_2_1_dis = 0;
	int signal_events_3_1_dis = 0;
	int signal_events_4_1_dis = 0;

	int signal_events_1_1_coh = 0;
	int signal_events_2_1_coh = 0;
	int signal_events_3_1_coh = 0;
	int signal_events_4_1_coh = 0;

	int signal_events_1_1_mec = 0;
	int signal_events_2_1_mec = 0;
	int signal_events_3_1_mec = 0;
	int signal_events_4_1_mec = 0;

	int signal_events_1_1_total = 0;
	int signal_events_2_1_total = 0;
	int signal_events_3_1_total = 0;
	int signal_events_4_1_total = 0;

	int bkg_events_1_1 = 0;
	int bkg_events_2_1 = 0;
	int bkg_events_3_1 = 0;
	int bkg_events_4_1 = 0;

	int intime_events_1_1 = 0;
	int intime_events_2_1 = 0;
	int intime_events_3_1 = 0;
	int intime_events_4_1 = 0;

	int signal_events_1_0_qe = 0;
	int signal_events_2_0_qe = 0;
	int signal_events_3_0_qe = 0;
	int signal_events_4_0_qe = 0;

	int signal_events_1_0_res = 0;
	int signal_events_2_0_res = 0;
	int signal_events_3_0_res = 0;
	int signal_events_4_0_res = 0;

	int signal_events_1_0_dis = 0;
	int signal_events_2_0_dis = 0;
	int signal_events_3_0_dis = 0;
	int signal_events_4_0_dis = 0;

	int signal_events_1_0_coh = 0;
	int signal_events_2_0_coh = 0;
	int signal_events_3_0_coh = 0;
	int signal_events_4_0_coh = 0;

	int signal_events_1_0_mec = 0;
	int signal_events_2_0_mec = 0;
	int signal_events_3_0_mec = 0;
	int signal_events_4_0_mec = 0;

	int signal_events_1_0_total = 0;
	int signal_events_2_0_total = 0;
	int signal_events_3_0_total = 0;
	int signal_events_4_0_total = 0;

	int bkg_events_1_0 = 0;
	int bkg_events_2_0 = 0;
	int bkg_events_3_0 = 0;
	int bkg_events_4_0 = 0;

	int intime_events_1_0 = 0;
	int intime_events_2_0 = 0;
	int intime_events_3_0 = 0;
	int intime_events_4_0 = 0;
	//this loops through all events which passed the selection cuts
	for(auto const my_tuple : * post_cuts_v)
	{
		const int event_num = std::get<0>(my_tuple);
		const int run_num = std::get<1>(my_tuple);
		const int sub_run_num = std::get<2>(my_tuple);
		const double pfp_vtx_x = std::get<3>(my_tuple);
		const double pfp_vtx_y = std::get<4>(my_tuple);
		const double pfp_vtx_z = std::get<5>(my_tuple);
		const std::string reason = std::get<6>(my_tuple);
		const std::string event_type = std::get<7>(my_tuple);
		const int num_tracks = std::get<8>(my_tuple);
		const int num_showers = std::get<9>(my_tuple);
		const double opening_angle = std::get<10>(my_tuple);

		int signal_bkg = 0;
		if(event_type == "nue_cc_qe"  || event_type == "nue_bar_cc_qe")   {signal_bkg = 1; }
		if(event_type == "nue_cc_res" || event_type == "nue_bar_cc_res")  {signal_bkg = 2; }
		if(event_type == "nue_cc_dis" || event_type == "nue_bar_cc_dis")  {signal_bkg = 3; }
		if(event_type == "nue_cc_coh" || event_type == "nue_bar_cc_coh")  {signal_bkg = 4; }
		if(event_type == "nue_cc_mec" || event_type == "nue_bar_cc_mec")  {signal_bkg = 5; }
		if(event_type == "Data")      {signal_bkg = 6; }
		//if(signal_bkg == 1)
		if(event_type == "nue_cc_qe" || event_type == "nue_bar_cc_qe")
		{
			if(num_tracks == 0)
			{
				if(num_showers == 1) {signal_events_1_0_qe++; }
				if(num_showers == 2) {signal_events_2_0_qe++; }
				if(num_showers == 3) {signal_events_3_0_qe++; }
				if(num_showers >= 4) {signal_events_4_0_qe++; }
			}
			if(num_tracks >= 1)
			{
				if(num_showers == 1) {signal_events_1_1_qe++; }
				if(num_showers == 2) {signal_events_2_1_qe++; }
				if(num_showers == 3) {signal_events_3_1_qe++; }
				if(num_showers >= 4) {signal_events_4_1_qe++; }
			}
		}
		//if(signal_bkg == 2)
		if(event_type == "nue_cc_res" || event_type == "nue_bar_cc_res")
		{
			if(num_tracks == 0)
			{
				if(num_showers == 1) {signal_events_1_0_res++; }
				if(num_showers == 2) {signal_events_2_0_res++; }
				if(num_showers == 3) {signal_events_3_0_res++; }
				if(num_showers >= 4) {signal_events_4_0_res++; }
			}
			if(num_tracks >= 1)
			{
				if(num_showers == 1) {signal_events_1_1_res++; }
				if(num_showers == 2) {signal_events_2_1_res++; }
				if(num_showers == 3) {signal_events_3_1_res++; }
				if(num_showers >= 4) {signal_events_4_1_res++; }
			}
		}
		//if(signal_bkg == 3)
		if(event_type == "nue_cc_dis" || event_type == "nue_bar_cc_dis")
		{
			if(num_tracks == 0)
			{
				if(num_showers == 1) {signal_events_1_0_dis++; }
				if(num_showers == 2) {signal_events_2_0_dis++; }
				if(num_showers == 3) {signal_events_3_0_dis++; }
				if(num_showers >= 4) {signal_events_4_0_dis++; }
			}
			if(num_tracks >= 1)
			{
				if(num_showers == 1) {signal_events_1_1_dis++; }
				if(num_showers == 2) {signal_events_2_1_dis++; }
				if(num_showers == 3) {signal_events_3_1_dis++; }
				if(num_showers >= 4) {signal_events_4_1_dis++; }
			}
		}
		//if(signal_bkg == 4)
		if(event_type == "nue_cc_coh" || event_type == "nue_bar_cc_coh")
		{
			if(num_tracks == 0)
			{
				if(num_showers == 1) {signal_events_1_0_coh++; }
				if(num_showers == 2) {signal_events_2_0_coh++; }
				if(num_showers == 3) {signal_events_3_0_coh++; }
				if(num_showers >= 4) {signal_events_4_0_coh++; }
			}
			if(num_tracks >= 1)
			{
				if(num_showers == 1) {signal_events_1_1_coh++; }
				if(num_showers == 2) {signal_events_2_1_coh++; }
				if(num_showers == 3) {signal_events_3_1_coh++; }
				if(num_showers >= 4) {signal_events_4_1_coh++; }
			}
		}
		//if(signal_bkg == 5)
		if(event_type == "nue_cc_mec" || event_type == "nue_bar_cc_mec")
		{
			if(num_tracks == 0)
			{
				if(num_showers == 1) {signal_events_1_0_mec++; }
				if(num_showers == 2) {signal_events_2_0_mec++; }
				if(num_showers == 3) {signal_events_3_0_mec++; }
				if(num_showers >= 4) {signal_events_4_0_mec++; }
			}
			if(num_tracks >= 1)
			{
				if(num_showers == 1) {signal_events_1_1_mec++; }
				if(num_showers == 2) {signal_events_2_1_mec++; }
				if(num_showers == 3) {signal_events_3_1_mec++; }
				if(num_showers >= 4) {signal_events_4_1_mec++; }
			}
		}
		if(signal_bkg == 0)
		{
			if(num_tracks == 0)
			{
				if(num_showers == 1) {bkg_events_1_0++; }
				if(num_showers == 2) {bkg_events_2_0++; }
				if(num_showers == 3) {bkg_events_3_0++; }
				if(num_showers >= 4) {bkg_events_4_0++; }
			}
			if(num_tracks >= 1)
			{
				if(num_showers == 1) {bkg_events_1_1++; }
				if(num_showers == 2) {bkg_events_2_1++; }
				if(num_showers == 3) {bkg_events_3_1++; }
				if(num_showers >= 4) {bkg_events_4_1++; }
			}
		}
		//for in-time, needs unique scaling
		if(signal_bkg == 6)
		{
			if(num_tracks == 0)
			{
				if(num_showers == 1) {intime_events_1_0++; }
				if(num_showers == 2) {intime_events_2_0++; }
				if(num_showers == 3) {intime_events_3_0++; }
				if(num_showers >= 4) {intime_events_4_0++; }
			}
			if(num_tracks >= 1)
			{
				if(num_showers == 1) {intime_events_1_1++; }
				if(num_showers == 2) {intime_events_2_1++; }
				if(num_showers == 3) {intime_events_3_1++; }
				if(num_showers >= 4) {intime_events_4_1++; }
			}
		}
	}
	const double total_events_1_0 = signal_events_1_0_qe + signal_events_1_0_res + signal_events_1_0_dis + signal_events_1_0_coh + signal_events_1_0_mec + bkg_events_1_0 + (intime_events_1_0 * intime_scale_factor / data_scale_factor);
	const double total_events_2_0 = signal_events_2_0_qe + signal_events_2_0_res + signal_events_2_0_dis + signal_events_2_0_coh + signal_events_2_0_mec + bkg_events_2_0 + (intime_events_2_0 * intime_scale_factor / data_scale_factor);
	const double total_events_3_0 = signal_events_3_0_qe + signal_events_3_0_res + signal_events_3_0_dis + signal_events_3_0_coh + signal_events_3_0_mec + bkg_events_3_0 + (intime_events_3_0 * intime_scale_factor / data_scale_factor);
	const double total_events_4_0 = signal_events_4_0_qe + signal_events_4_0_res + signal_events_4_0_dis + signal_events_4_0_coh + signal_events_4_0_mec + bkg_events_4_0 + (intime_events_4_0 * intime_scale_factor / data_scale_factor);

	const double total_signal_events_1_0 = signal_events_1_0_qe + signal_events_1_0_res + signal_events_1_0_dis + signal_events_1_0_coh + signal_events_1_0_mec;
	const double total_signal_events_2_0 = signal_events_2_0_qe + signal_events_2_0_res + signal_events_2_0_dis + signal_events_2_0_coh + signal_events_2_0_mec;
	const double total_signal_events_3_0 = signal_events_3_0_qe + signal_events_3_0_res + signal_events_3_0_dis + signal_events_3_0_coh + signal_events_3_0_mec;
	const double total_signal_events_4_0 = signal_events_4_0_qe + signal_events_4_0_res + signal_events_4_0_dis + signal_events_4_0_coh + signal_events_4_0_mec;

	const double total_events_1_1 = signal_events_1_1_qe + signal_events_1_1_res + signal_events_1_1_dis + signal_events_1_1_coh + signal_events_1_1_mec + bkg_events_1_1 + (intime_events_1_1 * intime_scale_factor / data_scale_factor);
	const double total_events_2_1 = signal_events_2_1_qe + signal_events_2_1_res + signal_events_2_1_dis + signal_events_2_1_coh + signal_events_2_1_mec + bkg_events_2_1 + (intime_events_2_1 * intime_scale_factor / data_scale_factor);
	const double total_events_3_1 = signal_events_3_1_qe + signal_events_3_1_res + signal_events_3_1_dis + signal_events_3_1_coh + signal_events_3_1_mec + bkg_events_3_1 + (intime_events_3_1 * intime_scale_factor / data_scale_factor);
	const double total_events_4_1 = signal_events_4_1_qe + signal_events_4_1_res + signal_events_4_1_dis + signal_events_4_1_coh + signal_events_4_1_mec + bkg_events_4_1 + (intime_events_4_1 * intime_scale_factor / data_scale_factor);

	const double total_signal_events_1_1 = signal_events_1_1_qe + signal_events_1_1_res + signal_events_1_1_dis + signal_events_1_1_coh + signal_events_1_1_mec;
	const double total_signal_events_2_1 = signal_events_2_1_qe + signal_events_2_1_res + signal_events_2_1_dis + signal_events_2_1_coh + signal_events_2_1_mec;
	const double total_signal_events_3_1 = signal_events_3_1_qe + signal_events_3_1_res + signal_events_3_1_dis + signal_events_3_1_coh + signal_events_3_1_mec;
	const double total_signal_events_4_1 = signal_events_4_1_qe + signal_events_4_1_res + signal_events_4_1_dis + signal_events_4_1_coh + signal_events_4_1_mec;

	std::cout << " ======================== " << std::endl;
	std::cout << "Total Purities: " << std::endl;
	std::cout << "1 Shower 0 Track: " << total_signal_events_1_0 << " , stat err from n signal: " << 1/sqrt(total_signal_events_1_0) << std::endl;
	std::cout << "2 Shower 0 Track: " << total_signal_events_2_0 << " , stat err from n signal: " << 1/sqrt(total_signal_events_2_0) << std::endl;
	std::cout << "3 Shower 0 Track: " << total_signal_events_3_0 << " , stat err from n signal: " << 1/sqrt(total_signal_events_3_0) << std::endl;
	std::cout << "4 Shower 0 Track: " << total_signal_events_4_0 << " , stat err from n signal: " << 1/sqrt(total_signal_events_4_0) << std::endl;
	std::cout << "1 Shower 1 Track: " << total_signal_events_1_1 << " , stat err from n signal: " << 1/sqrt(total_signal_events_1_1) << std::endl;
	std::cout << "2 Shower 1 Track: " << total_signal_events_2_1 << " , stat err from n signal: " << 1/sqrt(total_signal_events_2_1) << std::endl;
	std::cout << "3 Shower 1 Track: " << total_signal_events_3_1 << " , stat err from n signal: " << 1/sqrt(total_signal_events_3_1) << std::endl;
	std::cout << "4 Shower 1 Track: " << total_signal_events_4_1 << " , stat err from n signal: " << 1/sqrt(total_signal_events_4_1) << std::endl;
	std::cout << " ======================== " << std::endl;

	post_cuts_num_tracks_showers_signal_total->SetBinContent(1, 1, total_signal_events_1_0);
	post_cuts_num_tracks_showers_signal_total->SetBinContent(2, 1, total_signal_events_2_0);
	post_cuts_num_tracks_showers_signal_total->SetBinContent(3, 1, total_signal_events_3_0);
	post_cuts_num_tracks_showers_signal_total->SetBinContent(4, 1, total_signal_events_4_0);
	post_cuts_num_tracks_showers_signal_total->SetBinContent(1, 2, total_signal_events_1_1);
	post_cuts_num_tracks_showers_signal_total->SetBinContent(2, 2, total_signal_events_2_1);
	post_cuts_num_tracks_showers_signal_total->SetBinContent(3, 2, total_signal_events_3_1);
	post_cuts_num_tracks_showers_signal_total->SetBinContent(4, 2, total_signal_events_4_1);

	post_cuts_num_tracks_showers_bkg_total->SetBinContent(1, 1, bkg_events_1_0);
	post_cuts_num_tracks_showers_bkg_total->SetBinContent(2, 1, bkg_events_2_0);
	post_cuts_num_tracks_showers_bkg_total->SetBinContent(3, 1, bkg_events_3_0);
	post_cuts_num_tracks_showers_bkg_total->SetBinContent(4, 1, bkg_events_4_0);
	post_cuts_num_tracks_showers_bkg_total->SetBinContent(1, 2, bkg_events_1_1);
	post_cuts_num_tracks_showers_bkg_total->SetBinContent(2, 2, bkg_events_2_1);
	post_cuts_num_tracks_showers_bkg_total->SetBinContent(3, 2, bkg_events_3_1);
	post_cuts_num_tracks_showers_bkg_total->SetBinContent(4, 2, bkg_events_4_1);

	post_cuts_num_tracks_showers_total_total->SetBinContent(1, 1, total_events_1_0);
	post_cuts_num_tracks_showers_total_total->SetBinContent(2, 1, total_events_2_0);
	post_cuts_num_tracks_showers_total_total->SetBinContent(3, 1, total_events_3_0);
	post_cuts_num_tracks_showers_total_total->SetBinContent(4, 1, total_events_4_0);
	post_cuts_num_tracks_showers_total_total->SetBinContent(1, 2, total_events_1_1);
	post_cuts_num_tracks_showers_total_total->SetBinContent(2, 2, total_events_2_1);
	post_cuts_num_tracks_showers_total_total->SetBinContent(3, 2, total_events_3_1);
	post_cuts_num_tracks_showers_total_total->SetBinContent(4, 2, total_events_4_1);

	TCanvas * efficiency_c1 = new TCanvas();
	efficiency_c1->cd();

	TEfficiency * teff = new TEfficiency(*post_cuts_num_tracks_showers_signal_total, *post_cuts_num_tracks_showers_total_total);
	teff->SetTitle("Topology Purity");
	//teff->SetLineColor(kGreen+3);
	//teff->SetMarkerColor(kGreen+3);
	//teff->SetMarkerStyle(20);
	//teff->SetMarkerSize(0.5);
	teff->Draw("lego4 e");
	//teff->Draw("lego");
	efficiency_c1->Print("../scripts/plots/post_cuts_num_tracks_showers_purity_teff.pdf");


	double purity_1_0_qe  = double(signal_events_1_0_qe)  / total_events_1_0;
	double purity_2_0_qe  = double(signal_events_2_0_qe)  / total_events_2_0;
	double purity_3_0_qe  = double(signal_events_3_0_qe)  / total_events_3_0;
	double purity_4_0_qe  = double(signal_events_4_0_qe)  / total_events_4_0;
	double purity_1_1_qe  = double(signal_events_1_1_qe)  / total_events_1_1;
	double purity_2_1_qe  = double(signal_events_2_1_qe)  / total_events_2_1;
	double purity_3_1_qe  = double(signal_events_3_1_qe)  / total_events_3_1;
	double purity_4_1_qe  = double(signal_events_4_1_qe)  / total_events_4_1;

	double purity_1_0_res = double(signal_events_1_0_res) / total_events_1_0;
	double purity_2_0_res = double(signal_events_2_0_res) / total_events_2_0;
	double purity_3_0_res = double(signal_events_3_0_res) / total_events_3_0;
	double purity_4_0_res = double(signal_events_4_0_res) / total_events_4_0;
	double purity_1_1_res = double(signal_events_1_1_res) / total_events_1_1;
	double purity_2_1_res = double(signal_events_2_1_res) / total_events_2_1;
	double purity_3_1_res = double(signal_events_3_1_res) / total_events_3_1;
	double purity_4_1_res = double(signal_events_4_1_res) / total_events_4_1;

	double purity_1_0_dis = double(signal_events_1_0_dis) / total_events_1_0;
	double purity_2_0_dis = double(signal_events_2_0_dis) / total_events_2_0;
	double purity_3_0_dis = double(signal_events_3_0_dis) / total_events_3_0;
	double purity_4_0_dis = double(signal_events_4_0_dis) / total_events_4_0;
	double purity_1_1_dis = double(signal_events_1_1_dis) / total_events_1_1;
	double purity_2_1_dis = double(signal_events_2_1_dis) / total_events_2_1;
	double purity_3_1_dis = double(signal_events_3_1_dis) / total_events_3_1;
	double purity_4_1_dis = double(signal_events_4_1_dis) / total_events_4_1;

	double purity_1_0_coh = double(signal_events_1_0_coh) / total_events_1_0;
	double purity_2_0_coh = double(signal_events_2_0_coh) / total_events_2_0;
	double purity_3_0_coh = double(signal_events_3_0_coh) / total_events_3_0;
	double purity_4_0_coh = double(signal_events_4_0_coh) / total_events_4_0;
	double purity_1_1_coh = double(signal_events_1_1_coh) / total_events_1_1;
	double purity_2_1_coh = double(signal_events_2_1_coh) / total_events_2_1;
	double purity_3_1_coh = double(signal_events_3_1_coh) / total_events_3_1;
	double purity_4_1_coh = double(signal_events_4_1_coh) / total_events_4_1;

	double purity_1_0_mec = double(signal_events_1_0_mec) / total_events_1_0;
	double purity_2_0_mec = double(signal_events_2_0_mec) / total_events_2_0;
	double purity_3_0_mec = double(signal_events_3_0_mec) / total_events_3_0;
	double purity_4_0_mec = double(signal_events_4_0_mec) / total_events_4_0;
	double purity_1_1_mec = double(signal_events_1_1_mec) / total_events_1_1;
	double purity_2_1_mec = double(signal_events_2_1_mec) / total_events_2_1;
	double purity_3_1_mec = double(signal_events_3_1_mec) / total_events_3_1;
	double purity_4_1_mec = double(signal_events_4_1_mec) / total_events_4_1;

	// double purity_1_0_total = purity_1_0_qe + purity_1_0_res + purity_1_0_dis + purity_1_0_coh + purity_1_0_mec;
	// double purity_2_0_total = purity_2_0_qe + purity_2_0_res + purity_2_0_dis + purity_2_0_coh + purity_2_0_mec;
	// double purity_3_0_total = purity_3_0_qe + purity_3_0_res + purity_3_0_dis + purity_3_0_coh + purity_3_0_mec;
	// double purity_4_0_total = purity_4_0_qe + purity_4_0_res + purity_4_0_dis + purity_4_0_coh + purity_4_0_mec;
	// double purity_1_1_total = purity_1_1_qe + purity_1_1_res + purity_1_1_dis + purity_1_1_coh + purity_1_1_mec;
	// double purity_2_1_total = purity_2_1_qe + purity_2_1_res + purity_2_1_dis + purity_2_1_coh + purity_2_1_mec;
	// double purity_3_1_total = purity_3_1_qe + purity_3_1_res + purity_3_1_dis + purity_3_1_coh + purity_3_1_mec;
	// double purity_4_1_total = purity_4_1_qe + purity_4_1_res + purity_4_1_dis + purity_4_1_coh + purity_4_1_mec;

	double purity_1_0_total = double(signal_events_1_0_qe + signal_events_1_0_res + signal_events_1_0_dis + signal_events_1_0_coh + signal_events_1_0_mec) / total_events_1_0;
	double purity_2_0_total = double(signal_events_2_0_qe + signal_events_2_0_res + signal_events_2_0_dis + signal_events_2_0_coh + signal_events_2_0_mec) / total_events_2_0;
	double purity_3_0_total = double(signal_events_3_0_qe + signal_events_3_0_res + signal_events_3_0_dis + signal_events_3_0_coh + signal_events_3_0_mec) / total_events_3_0;
	double purity_4_0_total = double(signal_events_4_0_qe + signal_events_4_0_res + signal_events_4_0_dis + signal_events_4_0_coh + signal_events_4_0_mec) / total_events_4_0;
	double purity_1_1_total = double(signal_events_1_1_qe + signal_events_1_1_res + signal_events_1_1_dis + signal_events_1_1_coh + signal_events_1_1_mec) / total_events_1_1;
	double purity_2_1_total = double(signal_events_2_1_qe + signal_events_2_1_res + signal_events_2_1_dis + signal_events_2_1_coh + signal_events_2_1_mec) / total_events_2_1;
	double purity_3_1_total = double(signal_events_3_1_qe + signal_events_3_1_res + signal_events_3_1_dis + signal_events_3_1_coh + signal_events_3_1_mec) / total_events_3_1;
	double purity_4_1_total = double(signal_events_4_1_qe + signal_events_4_1_res + signal_events_4_1_dis + signal_events_4_1_coh + signal_events_4_1_mec) / total_events_4_1;

	post_cuts_num_tracks_showers_purity_qe->SetBinContent(1, 1, purity_1_0_qe);
	post_cuts_num_tracks_showers_purity_qe->SetBinContent(2, 1, purity_2_0_qe);
	post_cuts_num_tracks_showers_purity_qe->SetBinContent(3, 1, purity_3_0_qe);
	post_cuts_num_tracks_showers_purity_qe->SetBinContent(4, 1, purity_4_0_qe);
	post_cuts_num_tracks_showers_purity_qe->SetBinContent(1, 2, purity_1_1_qe);
	post_cuts_num_tracks_showers_purity_qe->SetBinContent(2, 2, purity_2_1_qe);
	post_cuts_num_tracks_showers_purity_qe->SetBinContent(3, 2, purity_3_1_qe);
	post_cuts_num_tracks_showers_purity_qe->SetBinContent(4, 2, purity_4_1_qe);

	post_cuts_num_tracks_showers_purity_res->SetBinContent(1, 1, purity_1_0_res);
	post_cuts_num_tracks_showers_purity_res->SetBinContent(2, 1, purity_2_0_res);
	post_cuts_num_tracks_showers_purity_res->SetBinContent(3, 1, purity_3_0_res);
	post_cuts_num_tracks_showers_purity_res->SetBinContent(4, 1, purity_4_0_res);
	post_cuts_num_tracks_showers_purity_res->SetBinContent(1, 2, purity_1_1_res);
	post_cuts_num_tracks_showers_purity_res->SetBinContent(2, 2, purity_2_1_res);
	post_cuts_num_tracks_showers_purity_res->SetBinContent(3, 2, purity_3_1_res);
	post_cuts_num_tracks_showers_purity_res->SetBinContent(4, 2, purity_4_1_res);

	post_cuts_num_tracks_showers_purity_dis->SetBinContent(1, 1, purity_1_0_dis);
	post_cuts_num_tracks_showers_purity_dis->SetBinContent(2, 1, purity_2_0_dis);
	post_cuts_num_tracks_showers_purity_dis->SetBinContent(3, 1, purity_3_0_dis);
	post_cuts_num_tracks_showers_purity_dis->SetBinContent(4, 1, purity_4_0_dis);
	post_cuts_num_tracks_showers_purity_dis->SetBinContent(1, 2, purity_1_1_dis);
	post_cuts_num_tracks_showers_purity_dis->SetBinContent(2, 2, purity_2_1_dis);
	post_cuts_num_tracks_showers_purity_dis->SetBinContent(3, 2, purity_3_1_dis);
	post_cuts_num_tracks_showers_purity_dis->SetBinContent(4, 2, purity_4_1_dis);

	post_cuts_num_tracks_showers_purity_coh->SetBinContent(1, 1, purity_1_0_coh);
	post_cuts_num_tracks_showers_purity_coh->SetBinContent(2, 1, purity_2_0_coh);
	post_cuts_num_tracks_showers_purity_coh->SetBinContent(3, 1, purity_3_0_coh);
	post_cuts_num_tracks_showers_purity_coh->SetBinContent(4, 1, purity_4_0_coh);
	post_cuts_num_tracks_showers_purity_coh->SetBinContent(1, 2, purity_1_1_coh);
	post_cuts_num_tracks_showers_purity_coh->SetBinContent(2, 2, purity_2_1_coh);
	post_cuts_num_tracks_showers_purity_coh->SetBinContent(3, 2, purity_3_1_coh);
	post_cuts_num_tracks_showers_purity_coh->SetBinContent(4, 2, purity_4_1_coh);

	post_cuts_num_tracks_showers_purity_mec->SetBinContent(1, 1, purity_1_0_mec);
	post_cuts_num_tracks_showers_purity_mec->SetBinContent(2, 1, purity_2_0_mec);
	post_cuts_num_tracks_showers_purity_mec->SetBinContent(3, 1, purity_3_0_mec);
	post_cuts_num_tracks_showers_purity_mec->SetBinContent(4, 1, purity_4_0_mec);
	post_cuts_num_tracks_showers_purity_mec->SetBinContent(1, 2, purity_1_1_mec);
	post_cuts_num_tracks_showers_purity_mec->SetBinContent(2, 2, purity_2_1_mec);
	post_cuts_num_tracks_showers_purity_mec->SetBinContent(3, 2, purity_3_1_mec);
	post_cuts_num_tracks_showers_purity_mec->SetBinContent(4, 2, purity_4_1_mec);

	post_cuts_num_tracks_showers_purity_total->SetBinContent(1, 1, purity_1_0_total);
	post_cuts_num_tracks_showers_purity_total->SetBinContent(2, 1, purity_2_0_total);
	post_cuts_num_tracks_showers_purity_total->SetBinContent(3, 1, purity_3_0_total);
	post_cuts_num_tracks_showers_purity_total->SetBinContent(4, 1, purity_4_0_total);
	post_cuts_num_tracks_showers_purity_total->SetBinContent(1, 2, purity_1_1_total);
	post_cuts_num_tracks_showers_purity_total->SetBinContent(2, 2, purity_2_1_total);
	post_cuts_num_tracks_showers_purity_total->SetBinContent(3, 2, purity_3_1_total);
	post_cuts_num_tracks_showers_purity_total->SetBinContent(4, 2, purity_4_1_total);
}
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingThetaPhi(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                          std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                          TH2D * h_ele_theta_phi_nue_cc,
                                          TH2D * h_ele_theta_phi_nue_cc_out_fv,
                                          TH2D * h_ele_theta_phi_nue_cc_mixed,
                                          TH2D * h_ele_theta_phi_numu_cc,
                                          TH2D * h_ele_theta_phi_numu_cc_mixed,
                                          TH2D * h_ele_theta_phi_nc,
                                          TH2D * h_ele_theta_phi_nc_pi0,
                                          TH2D * h_ele_theta_phi_cosmic,
                                          TH2D * h_ele_theta_phi_other_mixed,
                                          TH2D * h_ele_theta_phi_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shower_theta = acos(leading_shower.pfpDirZ()) * 180 / 3.1415;
		const double leading_shower_phi = atan2(leading_shower.pfpDirY(), leading_shower.pfpDirX()) * 180 / 3.1415;
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_ele_theta_phi_nue_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_ele_theta_phi_nue_cc_out_fv->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_ele_theta_phi_nue_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_ele_theta_phi_nue_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_ele_theta_phi_nue_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_theta_phi_nue_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_ele_theta_phi_numu_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_ele_theta_phi_numu_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_ele_theta_phi_numu_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_ele_theta_phi_numu_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_ele_theta_phi_numu_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "nc")
		{
			h_ele_theta_phi_nc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "nc_pi0")
		{
			h_ele_theta_phi_nc_pi0->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_ele_theta_phi_nue_cc_mixed->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_theta_phi_numu_cc_mixed->Fill(leading_shower_phi, leading_shower_theta);
			h_ele_theta_phi_numu_cc->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "cosmic")
		{
			h_ele_theta_phi_cosmic->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "other_mixed")
		{
			h_ele_theta_phi_other_mixed->Fill(leading_shower_phi, leading_shower_theta);
		}
		if(tpco_id == "unmatched")
		{
			h_ele_theta_phi_unmatched->Fill(leading_shower_phi, leading_shower_theta);
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingThetaPhiInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                TH2D * h_ele_theta_phi_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		//if(passed_tpco->at(i).first == 0) {continue; }
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
		const double leading_shower_theta = acos(leading_shower.pfpDirZ()) * 180 / 3.1415;
		const double leading_shower_phi = atan2(leading_shower.pfpDirY(), leading_shower.pfpDirX()) * 180 / 3.1415;
		h_ele_theta_phi_intime->Fill(leading_shower_phi, leading_shower_theta);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::XYZPosition(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                      std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                      const double mc_nu_vtx_x, const double mc_nu_vtx_y, const double mc_nu_vtx_z,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nue_cc,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nue_cc_out_fv,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nue_cc_mixed,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_numu_cc,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_numu_cc_mixed,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nc,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nc_pi0,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_cosmic,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_other_mixed,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_unmatched)
{
	TH2D * dummy = new TH2D();
	XYZPosition(tpc_object_container_v,
	            passed_tpco, _verbose,
	            tpco_classifier_v,
	            mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z,
	            h_ele_pfp_xyz_nue_cc,
	            h_ele_pfp_xyz_nue_cc_out_fv,
	            h_ele_pfp_xyz_nue_cc_mixed,
	            h_ele_pfp_xyz_numu_cc,
	            h_ele_pfp_xyz_numu_cc_mixed,
	            h_ele_pfp_xyz_nc,
	            h_ele_pfp_xyz_nc_pi0,
	            h_ele_pfp_xyz_cosmic,
	            h_ele_pfp_xyz_other_mixed,
	            h_ele_pfp_xyz_unmatched,
	            dummy,
	            dummy);
}
//***************************************************************************
//***************************************************************************
void selection_functions::XYZPosition(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                      std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                      const double mc_nu_vtx_x, const double mc_nu_vtx_y, const double mc_nu_vtx_z,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nue_cc,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nue_cc_out_fv,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nue_cc_mixed,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_numu_cc,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_numu_cc_mixed,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nc,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nc_pi0,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_cosmic,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_other_mixed,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_unmatched,
                                      TH2 * h_pfp_zy_vtx_nue_cc,
                                      TH2 * h_pfp_zy_vtx_all)

{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;

		const double tpco_vtx_x = tpc_obj.pfpVtxX();
		const double tpco_vtx_y = tpc_obj.pfpVtxY();
		const double tpco_vtx_z = tpc_obj.pfpVtxZ();

		h_pfp_zy_vtx_all->Fill(tpco_vtx_z, tpco_vtx_y);

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
			h_pfp_zy_vtx_nue_cc->Fill(tpco_vtx_z, tpco_vtx_y);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_ele_pfp_xyz_nue_cc_out_fv->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc_out_fv->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc_out_fv->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
			h_pfp_zy_vtx_nue_cc->Fill(tpco_vtx_z, tpco_vtx_y);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
			h_pfp_zy_vtx_nue_cc->Fill(tpco_vtx_z, tpco_vtx_y);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
			h_pfp_zy_vtx_nue_cc->Fill(tpco_vtx_z, tpco_vtx_y);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
			h_pfp_zy_vtx_nue_cc->Fill(tpco_vtx_z, tpco_vtx_y);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nc")
		{
			h_ele_pfp_xyz_nc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nc_pi0")
		{
			h_ele_pfp_xyz_nc_pi0->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nc_pi0->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nc_pi0->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_ele_pfp_xyz_nue_cc_mixed->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc_mixed->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc_mixed->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			// h_ele_pfp_xyz_numu_cc_mixed->at(0)->Fill(tpco_vtx_x);
			// h_ele_pfp_xyz_numu_cc_mixed->at(1)->Fill(tpco_vtx_y);
			// h_ele_pfp_xyz_numu_cc_mixed->at(2)->Fill(tpco_vtx_z);
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "cosmic")
		{
			h_ele_pfp_xyz_cosmic->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_cosmic->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_cosmic->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "other_mixed")
		{
			h_ele_pfp_xyz_other_mixed->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_other_mixed->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_other_mixed->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "unmatched")
		{
			h_ele_pfp_xyz_unmatched->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_unmatched->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_unmatched->at(2)->Fill(tpco_vtx_z);
		}
	}                                        //end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::XYZPosition(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                      std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                      const double mc_nu_vtx_x, const double mc_nu_vtx_y, const double mc_nu_vtx_z,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nue_cc,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nue_cc_out_fv,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nue_cc_mixed,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_numu_cc,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_numu_cc_mixed,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nc,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_nc_pi0,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_cosmic,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_other_mixed,
                                      std::vector<TH1 *> * h_ele_pfp_xyz_unmatched,
                                      TH2 * h_mc_vtx_xy_nue_cc,
                                      TH2 * h_mc_vtx_xz_nue_cc,
                                      TH2 * h_mc_vtx_yz_nue_cc,
                                      TH2 * h_reco_vtx_xy_nue_cc,
                                      TH2 * h_reco_vtx_xz_nue_cc,
                                      TH2 * h_reco_vtx_yz_nue_cc,
                                      TH2 * h_mc_vtx_xy_nue_cc_out_fv,
                                      TH2 * h_mc_vtx_xz_nue_cc_out_fv,
                                      TH2 * h_mc_vtx_yz_nue_cc_out_fv,
                                      TH2 * h_reco_vtx_xy_nue_cc_out_fv,
                                      TH2 * h_reco_vtx_xz_nue_cc_out_fv,
                                      TH2 * h_reco_vtx_yz_nue_cc_out_fv,
                                      TH2 * h_mc_reco_vtx_x_nue_cc,
                                      TH2 * h_mc_reco_vtx_y_nue_cc,
                                      TH2 * h_mc_reco_vtx_z_nue_cc,
                                      TH2 * h_mc_reco_vtx_x_nue_cc_out_fv,
                                      TH2 * h_mc_reco_vtx_y_nue_cc_out_fv,
                                      TH2 * h_mc_reco_vtx_z_nue_cc_out_fv)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;

		const double tpco_vtx_x = tpc_obj.pfpVtxX();
		const double tpco_vtx_y = tpc_obj.pfpVtxY();
		const double tpco_vtx_z = tpc_obj.pfpVtxZ();

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
			h_mc_vtx_xy_nue_cc->Fill(mc_nu_vtx_x, mc_nu_vtx_y);
			h_mc_vtx_xz_nue_cc->Fill(mc_nu_vtx_x, mc_nu_vtx_z);
			h_mc_vtx_yz_nue_cc->Fill(mc_nu_vtx_z, mc_nu_vtx_y);
			h_reco_vtx_xy_nue_cc->Fill(tpco_vtx_x, tpco_vtx_y);
			h_reco_vtx_xz_nue_cc->Fill(tpco_vtx_x, tpco_vtx_z);
			h_reco_vtx_yz_nue_cc->Fill(tpco_vtx_z, tpco_vtx_y);
			h_mc_reco_vtx_x_nue_cc->Fill(mc_nu_vtx_x, tpco_vtx_x);
			h_mc_reco_vtx_y_nue_cc->Fill(mc_nu_vtx_y, tpco_vtx_y);
			h_mc_reco_vtx_z_nue_cc->Fill(mc_nu_vtx_z, tpco_vtx_z);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_ele_pfp_xyz_nue_cc_out_fv->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc_out_fv->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc_out_fv->at(2)->Fill(tpco_vtx_z);
			h_mc_vtx_xy_nue_cc_out_fv->Fill(mc_nu_vtx_x, mc_nu_vtx_y);
			h_mc_vtx_xz_nue_cc_out_fv->Fill(mc_nu_vtx_x, mc_nu_vtx_z);
			h_mc_vtx_yz_nue_cc_out_fv->Fill(mc_nu_vtx_z, mc_nu_vtx_y);
			h_reco_vtx_xy_nue_cc_out_fv->Fill(tpco_vtx_x, tpco_vtx_y);
			h_reco_vtx_xz_nue_cc_out_fv->Fill(tpco_vtx_x, tpco_vtx_z);
			h_reco_vtx_yz_nue_cc_out_fv->Fill(tpco_vtx_z, tpco_vtx_y);
			h_mc_reco_vtx_x_nue_cc_out_fv->Fill(mc_nu_vtx_x, tpco_vtx_x);
			h_mc_reco_vtx_y_nue_cc_out_fv->Fill(mc_nu_vtx_y, tpco_vtx_y);
			h_mc_reco_vtx_z_nue_cc_out_fv->Fill(mc_nu_vtx_z, tpco_vtx_z);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_pfp_xyz_nue_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nc")
		{
			h_ele_pfp_xyz_nc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nc_pi0")
		{
			h_ele_pfp_xyz_nc_pi0->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nc_pi0->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nc_pi0->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_ele_pfp_xyz_nue_cc_mixed->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_nue_cc_mixed->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_nue_cc_mixed->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			// h_ele_pfp_xyz_numu_cc_mixed->at(0)->Fill(tpco_vtx_x);
			// h_ele_pfp_xyz_numu_cc_mixed->at(1)->Fill(tpco_vtx_y);
			// h_ele_pfp_xyz_numu_cc_mixed->at(2)->Fill(tpco_vtx_z);
			h_ele_pfp_xyz_numu_cc->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_numu_cc->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_numu_cc->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "cosmic")
		{
			h_ele_pfp_xyz_cosmic->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_cosmic->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_cosmic->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "other_mixed")
		{
			h_ele_pfp_xyz_other_mixed->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_other_mixed->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_other_mixed->at(2)->Fill(tpco_vtx_z);
		}
		if(tpco_id == "unmatched")
		{
			h_ele_pfp_xyz_unmatched->at(0)->Fill(tpco_vtx_x);
			h_ele_pfp_xyz_unmatched->at(1)->Fill(tpco_vtx_y);
			h_ele_pfp_xyz_unmatched->at(2)->Fill(tpco_vtx_z);
		}
	}//end pfp loop
}
void selection_functions::XYZPositionInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                            std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                            std::vector<TH1 *> * h_ele_pfp_xyz_intime,
                                            TH2 * h_pfp_zy_vtx_ext)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);

		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		h_ele_pfp_xyz_intime->at(0)->Fill(tpc_vtx_x);
		h_ele_pfp_xyz_intime->at(1)->Fill(tpc_vtx_y);
		h_ele_pfp_xyz_intime->at(2)->Fill(tpc_vtx_z);
		h_pfp_zy_vtx_ext->Fill(tpc_vtx_z, tpc_vtx_y);
	}
}
void selection_functions::XYZPositionInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                            std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                            std::vector<TH1 *> * h_ele_pfp_xyz_intime)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);

		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		h_ele_pfp_xyz_intime->at(0)->Fill(tpc_vtx_x);
		h_ele_pfp_xyz_intime->at(1)->Fill(tpc_vtx_y);
		h_ele_pfp_xyz_intime->at(2)->Fill(tpc_vtx_z);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::EnergyCosTheta(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                         std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                         TH2 * h_ele_eng_costheta_nue_cc,
                                         TH2 * h_ele_eng_costheta_nue_cc_out_fv,
                                         TH2 * h_ele_eng_costheta_nue_cc_mixed,
                                         TH2 * h_ele_eng_costheta_numu_cc,
                                         TH2 * h_ele_eng_costheta_numu_cc_mixed,
                                         TH2 * h_ele_eng_costheta_nc,
                                         TH2 * h_ele_eng_costheta_nc_pi0,
                                         TH2 * h_ele_eng_costheta_cosmic,
                                         TH2 * h_ele_eng_costheta_other_mixed,
                                         TH2 * h_ele_eng_costheta_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double momentum = leading_shower.pfpMomentum();
		const double costheta = leading_shower.pfpDirZ();

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_ele_eng_costheta_nue_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_ele_eng_costheta_nue_cc_out_fv->Fill(momentum, costheta);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_ele_eng_costheta_nue_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_ele_eng_costheta_nue_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_ele_eng_costheta_nue_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_eng_costheta_nue_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_ele_eng_costheta_numu_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_ele_eng_costheta_numu_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_ele_eng_costheta_numu_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_ele_eng_costheta_numu_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_ele_eng_costheta_numu_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "nc")
		{
			h_ele_eng_costheta_nc->Fill(momentum, costheta);
		}
		if(tpco_id == "nc_pi0")
		{
			h_ele_eng_costheta_nc_pi0->Fill(momentum, costheta);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_ele_eng_costheta_nue_cc_mixed->Fill(momentum, costheta);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_eng_costheta_numu_cc_mixed->Fill(momentum, costheta);
			h_ele_eng_costheta_numu_cc->Fill(momentum, costheta);
		}
		if(tpco_id == "cosmic")
		{
			h_ele_eng_costheta_cosmic->Fill(momentum, costheta);
		}
		if(tpco_id == "other_mixed")
		{
			h_ele_eng_costheta_other_mixed->Fill(momentum, costheta);
		}
		if(tpco_id == "unmatched")
		{
			h_ele_eng_costheta_unmatched->Fill(momentum, costheta);
		}
	}//end pfp loop
}
void selection_functions::EnergyCosThetaInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                               TH2 * h_ele_eng_costheta_intime)
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
		const double momentum = leading_shower.pfpMomentum();
		const double costheta = leading_shower.pfpDirZ();
		h_ele_eng_costheta_intime->Fill(momentum, costheta);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::TrueRecoEle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                      std::vector<std::pair<std::string, int> > * tpco_classifier_v, double mc_ele_momentum, double mc_ele_cos_theta,
                                      TH2D * h_true_reco_ele_momentum, TH2D * h_true_reco_ele_costheta, TH1D * h_true_num_e)
{
	int n_tpc_obj = tpc_object_container_v->size();
	int this_event_num_e = 0;
	int mc_pdg_code = 0;
	double pfp_ele_momentum = 0;
	double pfp_ele_costheta = 0;
	//this asks if there is a second true electron in the event with lower energy
	bool second_lower = false;
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		std::string tpco_id = tpco_classifier_v->at(i).first;

		for(int j = 0; j < n_pfp; j++)
		{
			auto const part = tpc_obj.GetParticle(j);
			mc_pdg_code = part.MCPdgCode();
			if(mc_pdg_code == 11 || mc_pdg_code == -11)
			{
				if(this_event_num_e > 0)
				{
					if(part.pfpMomentum() < pfp_ele_momentum) {second_lower = true; }
					if(part.pfpMomentum() > pfp_ele_momentum) {second_lower = false; }
				}
				if(second_lower == false)
				{
					pfp_ele_momentum = part.pfpMomentum();
					pfp_ele_costheta = part.pfpDirZ();
				}
				this_event_num_e++;
			}
		}//end loop pfps
		 //check if true nue cc event
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_mec" ||
		   tpco_id == "nue_bar_cc_qe" || tpco_id == "nue_bar_cc_res" || tpco_id == "nue_bar_cc_coh" || tpco_id == "nue_bar_cc_dis" || tpco_id == "nue_bar_cc_mec")
		{
			//when this_event_num_e == 0, we had a true nue cc interaction, and a shower was reco, but it wasn't the true e
			if(this_event_num_e != 0) {h_true_reco_ele_momentum->Fill(mc_ele_momentum, pfp_ele_momentum); }
			if(this_event_num_e != 0) {h_true_reco_ele_costheta->Fill(mc_ele_cos_theta, pfp_ele_costheta); }
			//this one is the number of reconstructed true electrons
			h_true_num_e->Fill(this_event_num_e);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::EnergyCosThetaSlices(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                               const double theta_translation, const double phi_translation,
                                               std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                               TH1 * h_ele_eng_for_nue_cc,
                                               TH1 * h_ele_eng_for_nue_cc_out_fv,
                                               TH1 * h_ele_eng_for_nue_cc_mixed,
                                               TH1 * h_ele_eng_for_numu_cc,
                                               TH1 * h_ele_eng_for_numu_cc_mixed,
                                               TH1 * h_ele_eng_for_nc,
                                               TH1 * h_ele_eng_for_nc_pi0,
                                               TH1 * h_ele_eng_for_cosmic,
                                               TH1 * h_ele_eng_for_other_mixed,
                                               TH1 * h_ele_eng_for_unmatched,
                                               TH1 * h_ele_eng_mid_nue_cc,
                                               TH1 * h_ele_eng_mid_nue_cc_out_fv,
                                               TH1 * h_ele_eng_mid_nue_cc_mixed,
                                               TH1 * h_ele_eng_mid_numu_cc,
                                               TH1 * h_ele_eng_mid_numu_cc_mixed,
                                               TH1 * h_ele_eng_mid_nc,
                                               TH1 * h_ele_eng_mid_nc_pi0,
                                               TH1 * h_ele_eng_mid_cosmic,
                                               TH1 * h_ele_eng_mid_other_mixed,
                                               TH1 * h_ele_eng_mid_unmatched,
                                               TH1 * h_ele_eng_back_nue_cc,
                                               TH1 * h_ele_eng_back_nue_cc_out_fv,
                                               TH1 * h_ele_eng_back_nue_cc_mixed,
                                               TH1 * h_ele_eng_back_numu_cc,
                                               TH1 * h_ele_eng_back_numu_cc_mixed,
                                               TH1 * h_ele_eng_back_nc,
                                               TH1 * h_ele_eng_back_nc_pi0,
                                               TH1 * h_ele_eng_back_cosmic,
                                               TH1 * h_ele_eng_back_other_mixed,
                                               TH1 * h_ele_eng_back_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double momentum = leading_shower.pfpMomentum();
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, theta_translation, phi_translation);
		const double costheta = shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag());

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			if(costheta >= 0.5) {h_ele_eng_for_nue_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_nue_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_nue_cc->Fill(momentum); }
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			if(costheta >= 0.5) {h_ele_eng_for_nue_cc_out_fv->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_nue_cc_out_fv->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_nue_cc_out_fv->Fill(momentum); }
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			if(costheta >= 0.5) {h_ele_eng_for_nue_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_nue_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_nue_cc->Fill(momentum); }
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			if(costheta >= 0.5) {h_ele_eng_for_nue_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_nue_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_nue_cc->Fill(momentum); }
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			if(costheta >= 0.5) {h_ele_eng_for_nue_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_nue_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_nue_cc->Fill(momentum); }
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			if(costheta >= 0.5) {h_ele_eng_for_nue_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_nue_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_nue_cc->Fill(momentum); }
		}
		if(tpco_id == "numu_cc_qe")
		{
			if(costheta >= 0.5) {h_ele_eng_for_numu_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_numu_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_numu_cc->Fill(momentum); }
		}
		if(tpco_id == "numu_cc_res")
		{
			if(costheta >= 0.5) {h_ele_eng_for_numu_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_numu_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_numu_cc->Fill(momentum); }
		}
		if(tpco_id == "numu_cc_dis")
		{
			if(costheta >= 0.5) {h_ele_eng_for_numu_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_numu_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_numu_cc->Fill(momentum); }
		}
		if(tpco_id == "numu_cc_coh")
		{
			if(costheta >= 0.5) {h_ele_eng_for_numu_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_numu_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_numu_cc->Fill(momentum); }
		}
		if(tpco_id == "numu_cc_mec")
		{
			if(costheta >= 0.5) {h_ele_eng_for_numu_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_numu_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_numu_cc->Fill(momentum); }
		}
		if(tpco_id == "nc")
		{
			if(costheta >= 0.5) {h_ele_eng_for_nc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_nc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_nc->Fill(momentum); }
		}
		if(tpco_id == "nc_pi0")
		{
			if(costheta >= 0.5) {h_ele_eng_for_nc_pi0->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_nc_pi0->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_nc_pi0->Fill(momentum); }
		}
		if(tpco_id == "nue_cc_mixed")
		{
			if(costheta >= 0.5) {h_ele_eng_for_nue_cc_mixed->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_nue_cc_mixed->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_nc->Fill(momentum); }
		}
		if(tpco_id == "numu_cc_mixed")
		{
			// if(costheta >= 0.5) {h_ele_eng_for_numu_cc_mixed->Fill(momentum); }
			// if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_numu_cc_mixed->Fill(momentum); }
			// if(costheta <= -0.5) {h_ele_eng_back_numu_cc_mixed->Fill(momentum); }
			if(costheta >= 0.5) {h_ele_eng_for_numu_cc->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_numu_cc->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_numu_cc->Fill(momentum); }
		}
		if(tpco_id == "cosmic")
		{
			if(costheta >= 0.5) {h_ele_eng_for_cosmic->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_cosmic->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_cosmic->Fill(momentum); }
		}
		if(tpco_id == "other_mixed")
		{
			if(costheta >= 0.5) {h_ele_eng_for_other_mixed->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_other_mixed->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_other_mixed->Fill(momentum); }
		}
		if(tpco_id == "unmatched")
		{
			if(costheta >= 0.5) {h_ele_eng_for_unmatched->Fill(momentum); }
			if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_unmatched->Fill(momentum); }
			if(costheta <= -0.5) {h_ele_eng_back_unmatched->Fill(momentum); }
		}
	}//end pfp loop
}
void selection_functions::EnergyCosThetaSlicesInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                     const double theta_translation, const double phi_translation,
                                                     TH1 * h_ele_eng_for_intime,
                                                     TH1 * h_ele_eng_mid_intime,
                                                     TH1 * h_ele_eng_back_intime)
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
		const double momentum = leading_shower.pfpMomentum();
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, theta_translation, phi_translation);
		const double costheta = shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag());
		if(costheta >= 0.5) {h_ele_eng_for_intime->Fill(momentum); }
		if(costheta < 0.5 && costheta > -0.5) {h_ele_eng_mid_intime->Fill(momentum); }
		if(costheta <= -0.5) {h_ele_eng_back_intime->Fill(momentum); }
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::DifferentialEnergySlices(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                   std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                   TH1 * h_low_true_momentum, TH1 * h_med_true_momentum, TH1 * h_high_true_momentum)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj      = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_mec" ||
		   tpco_id == "nue_bar_cc_qe" || tpco_id == "nue_bar_cc_res" || tpco_id == "nue_bar_cc_dis" || tpco_id == "nue_bar_cc_coh" || tpco_id == "nue_bar_cc_mec")
		{
			for(int j = 0; j < tpc_obj.NumPFParticles(); j++)
			{
				auto const pfp = tpc_obj.GetParticle(j);
				const int mc_pdg = pfp.MCPdgCode();
				if(mc_pdg == 11)
				{
					const double momentum = pfp.pfpMomentum();
					const double mc_momentum = pfp.mcMomentum();
					if(momentum <= 0.5) {h_low_true_momentum->Fill(mc_momentum); }
					if(momentum <= 1.0) {h_med_true_momentum->Fill(mc_momentum); }
					if(momentum <= 1.5) {h_high_true_momentum->Fill(mc_momentum); }
				}
			}
		}
	}
}
//***************************************************************************
void selection_functions::IsContainedPlot(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                          std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                          std::vector<double> fv_boundary_v,
                                          TH1 * h_track_containment_nue_cc,
                                          TH1 * h_track_containment_nue_cc_out_fv,
                                          TH1 * h_track_containment_nue_cc_mixed,
                                          TH1 * h_track_containment_numu_cc,
                                          TH1 * h_track_containment_numu_cc_mixed,
                                          TH1 * h_track_containment_nc,
                                          TH1 * h_track_containment_nc_pi0,
                                          TH1 * h_track_containment_cosmic,
                                          TH1 * h_track_containment_other_mixed,
                                          TH1 * h_track_containment_unmatched)
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
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;

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

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_track_containment_nue_cc->Fill(is_contained);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_track_containment_nue_cc_out_fv->Fill(is_contained);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_track_containment_nue_cc->Fill(is_contained);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_track_containment_nue_cc->Fill(is_contained);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_track_containment_nue_cc->Fill(is_contained);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_track_containment_nue_cc->Fill(is_contained);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_track_containment_numu_cc->Fill(is_contained);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_track_containment_numu_cc->Fill(is_contained);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_track_containment_numu_cc->Fill(is_contained);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_track_containment_numu_cc->Fill(is_contained);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_track_containment_numu_cc->Fill(is_contained);
		}
		if(tpco_id == "nc")
		{
			h_track_containment_nc->Fill(is_contained);
		}
		if(tpco_id == "nc_pi0")
		{
			h_track_containment_nc_pi0->Fill(is_contained);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_track_containment_nue_cc_mixed->Fill(is_contained);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_eng_costheta_numu_cc_mixed->Fill(momentum, costheta);
			h_track_containment_numu_cc->Fill(is_contained);
		}
		if(tpco_id == "cosmic")
		{
			h_track_containment_cosmic->Fill(is_contained);
		}
		if(tpco_id == "other_mixed")
		{
			h_track_containment_other_mixed->Fill(is_contained);
		}
		if(tpco_id == "unmatched")
		{
			h_track_containment_unmatched->Fill(is_contained);
		}
	}//end pfp loop
}
//***************************************************************************
void selection_functions::IsContainedPlotInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                std::vector<double> fv_boundary_v, TH1 * h_track_containment_intime)
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
		h_track_containment_intime->Fill(is_contained);
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::dEdxCollectionAngle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                              std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                              TH2D * h_dedx_collection_angle_nue_cc,
                                              TH2D * h_dedx_collection_angle_nue_cc_mixed,
                                              TH2D * h_dedx_collection_angle_nue_cc_out_fv,
                                              TH2D * h_dedx_collection_angle_numu_cc,
                                              TH2D * h_dedx_collection_angle_nc,
                                              TH2D * h_dedx_collection_angle_cosmic,
                                              TH2D * h_dedx_collection_angle_nc_pi0,
                                              TH2D * h_dedx_collection_angle_numu_cc_mixed,
                                              TH2D * h_dedx_collection_angle_other_mixed,
                                              TH2D * h_dedx_collection_angle_unmatched     )
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{

		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_dir_y = leading_shower.pfpDirY();
		const double leading_dir_z = leading_shower.pfpDirZ();
		const double leading_angle_collection = atan2(leading_dir_y, leading_dir_z);
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_dedx_collection_angle_nue_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_dedx_collection_angle_nue_cc_out_fv->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_dedx_collection_angle_nue_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_dedx_collection_angle_nue_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_dedx_collection_angle_nue_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_dedx_collection_angle_nue_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_dedx_collection_angle_numu_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_dedx_collection_angle_numu_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_dedx_collection_angle_numu_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_dedx_collection_angle_numu_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_dedx_collection_angle_numu_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "nc")
		{
			h_dedx_collection_angle_nc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "nc_pi0")
		{
			h_dedx_collection_angle_nc_pi0->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_dedx_collection_angle_nue_cc_mixed->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_dedx_cuts_numu_cc_mixed->Fill(leading_dedx);
			h_dedx_collection_angle_numu_cc->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "cosmic")
		{
			h_dedx_collection_angle_cosmic->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "other_mixed")
		{
			h_dedx_collection_angle_other_mixed->Fill(leading_dedx, leading_angle_collection);
		}
		if(tpco_id == "unmatched")
		{
			h_dedx_collection_angle_unmatched->Fill(leading_dedx, leading_angle_collection);
		}
	}        //end loop tpc objects
}

void selection_functions::dEdxCollectionAngleInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                    TH2D * h_dedx_collection_angle_intime)
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
		const double leading_dir_y = leading_shower.pfpDirY();
		const double leading_dir_z = leading_shower.pfpDirZ();
		const double leading_angle_collection = atan2(leading_dir_y, leading_dir_z);
		h_dedx_collection_angle_intime->Fill(leading_dedx * (242.72 / 196.979), leading_angle_collection);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsVtxFlashUpstream(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                   std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                   TH1D * h_vtx_flash_nue_cc, TH1D * h_vtx_flash_nue_cc_out_fv, TH1D * h_vtx_flash_nue_cc_mixed,
                                                   TH1D * h_vtx_flash_numu_cc, TH1D * h_vtx_flash_nc,
                                                   TH1D * h_vtx_flash_cosmic, TH1D * h_vtx_flash_nc_pi0,
                                                   TH1D * h_vtx_flash_numu_cc_mixed, TH1D * h_vtx_flash_other_mixed,
                                                   TH1D * h_vtx_flash_unmatched)
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
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_vtx_flash_nue_cc_out_fv->Fill(distance);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "nc")
		{
			h_vtx_flash_nc->Fill(distance);
		}
		if(tpco_id == "nc_pi0")
		{
			h_vtx_flash_nc_pi0->Fill(distance);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_vtx_flash_nue_cc_mixed->Fill(distance);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_vtx_flash_numu_cc_mixed->Fill(distance);
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "cosmic")
		{
			h_vtx_flash_cosmic->Fill(distance);
		}
		if(tpco_id == "other_mixed")
		{
			h_vtx_flash_other_mixed->Fill(distance);
		}
		if(tpco_id == "unmatched")
		{
			h_vtx_flash_unmatched->Fill(distance);
		}
	}        //end loop tpc objects
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsVtxFlashUpstreamInTime(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, TH1D * h_vtx_flash_intime)
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
		h_vtx_flash_intime->Fill(distance);
	}        //end loop tpc objects
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsVtxFlashDownstream(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                     std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                     std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                     TH1D * h_vtx_flash_nue_cc, TH1D * h_vtx_flash_nue_cc_out_fv, TH1D * h_vtx_flash_nue_cc_mixed,
                                                     TH1D * h_vtx_flash_numu_cc, TH1D * h_vtx_flash_nc,
                                                     TH1D * h_vtx_flash_cosmic, TH1D * h_vtx_flash_nc_pi0,
                                                     TH1D * h_vtx_flash_numu_cc_mixed, TH1D * h_vtx_flash_other_mixed,
                                                     TH1D * h_vtx_flash_unmatched)
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
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_vtx_flash_nue_cc_out_fv->Fill(distance);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "nc")
		{
			h_vtx_flash_nc->Fill(distance);
		}
		if(tpco_id == "nc_pi0")
		{
			h_vtx_flash_nc_pi0->Fill(distance);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_vtx_flash_nue_cc_mixed->Fill(distance);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_vtx_flash_numu_cc_mixed->Fill(distance);
			h_vtx_flash_numu_cc->Fill(distance);
		}
		if(tpco_id == "cosmic")
		{
			h_vtx_flash_cosmic->Fill(distance);
		}
		if(tpco_id == "other_mixed")
		{
			h_vtx_flash_other_mixed->Fill(distance);
		}
		if(tpco_id == "unmatched")
		{
			h_vtx_flash_unmatched->Fill(distance);
		}
	}        //end loop tpc objects
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsVtxFlashDownstreamInTime(std::vector< double > largest_flash_v, std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, TH1D * h_vtx_flash_intime)
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
		if(flash_vtx_z > tpc_vtx_z) {continue; }//skip this event if the flash is upstream
		const int n_pfp = tpc_obj.NumPFParticles();
		const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
		h_vtx_flash_intime->Fill(distance);
	}        //end loop tpc objects
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::FillTrueRecoEnergy(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco,
                                             std::vector<std::pair<std::string, int> > * tpco_classifier_v, double mc_ele_energy,
                                             TH1D * h_mc_ele_energy,
                                             TH1D * h_reco_ele_energy,
                                             TH2D * h_mc_reco_ele_energy)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double reco_ele_momentum = leading_shower.pfpMomentum();
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_mc_ele_energy->Fill(mc_ele_energy);
			h_reco_ele_energy->Fill(reco_ele_momentum);
			h_mc_reco_ele_energy->Fill(mc_ele_energy, reco_ele_momentum);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_mc_ele_energy->Fill(mc_ele_energy);
			h_reco_ele_energy->Fill(reco_ele_momentum);
			h_mc_reco_ele_energy->Fill(mc_ele_energy, reco_ele_momentum);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_mc_ele_energy->Fill(mc_ele_energy);
			h_reco_ele_energy->Fill(reco_ele_momentum);
			h_mc_reco_ele_energy->Fill(mc_ele_energy, reco_ele_momentum);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_mc_ele_energy->Fill(mc_ele_energy);
			h_reco_ele_energy->Fill(reco_ele_momentum);
			h_mc_reco_ele_energy->Fill(mc_ele_energy, reco_ele_momentum);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_mc_ele_energy->Fill(mc_ele_energy);
			h_reco_ele_energy->Fill(reco_ele_momentum);
			h_mc_reco_ele_energy->Fill(mc_ele_energy, reco_ele_momentum);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsLeadingMomentumTrueParticle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                              std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                              TH1D * h_leading_momentum_electron,
                                                              TH1D * h_leading_momentum_photon,
                                                              TH1D * h_leading_momentum_proton,
                                                              TH1D * h_leading_momentum_pion,
                                                              TH1D * h_leading_momentum_muon,
                                                              TH1D * h_leading_momentum_kaon,
                                                              TH1D * h_leading_momentum_neutron,
                                                              TH1D * h_leading_momentum_mc_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_momentum = leading_shower.pfpMomentum();
		const int leading_mc_pdg = leading_shower.MCPdgCode();
		bool good_id = false;
		if(leading_mc_pdg == 11 || leading_mc_pdg == -11)
		{
			h_leading_momentum_electron->Fill(leading_momentum);
			good_id = true;
		}
		if(leading_mc_pdg == 13 || leading_mc_pdg == -13)
		{
			h_leading_momentum_muon->Fill(leading_momentum);
			good_id = true;
		}
		if(leading_mc_pdg == 22)
		{
			h_leading_momentum_photon->Fill(leading_momentum);
			good_id = true;
		}
		if(leading_mc_pdg == 2212)
		{
			h_leading_momentum_proton->Fill(leading_momentum);
			good_id = true;
		}
		if(leading_mc_pdg == 211 || leading_mc_pdg == -211)
		{
			h_leading_momentum_pion->Fill(leading_momentum);
			good_id = true;
		}
		if(leading_mc_pdg == 2112)
		{
			h_leading_momentum_neutron->Fill(leading_momentum);
			good_id = true;
		}
		if(leading_mc_pdg == 130 || leading_mc_pdg == 310 || leading_mc_pdg == 311 || leading_mc_pdg == 321 || leading_mc_pdg == -321)
		{
			h_leading_momentum_kaon->Fill(leading_momentum);
			good_id = true;
		}
		if(good_id == false )
		{
			h_leading_momentum_mc_unmatched->Fill(leading_momentum);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsLeadingMomentumTrueParticleInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                                    TH1D * h_leading_momentum_type_cosmic)
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
		h_leading_momentum_type_cosmic->Fill(leading_momentum);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsLeadingMomentum(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                  std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                  TH1D * h_ele_momentum_nue_cc,
                                                  TH1D * h_ele_momentum_nue_cc_out_fv,
                                                  TH1D * h_ele_momentum_nue_cc_mixed,
                                                  TH1D * h_ele_momentum_numu_cc,
                                                  TH1D * h_ele_momentum_numu_cc_mixed,
                                                  TH1D * h_ele_momentum_nc,
                                                  TH1D * h_ele_momentum_nc_pi0,
                                                  TH1D * h_ele_momentum_cosmic,
                                                  TH1D * h_ele_momentum_other_mixed,
                                                  TH1D * h_ele_momentum_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_momentum = leading_shower.pfpMomentum();
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_ele_momentum_nue_cc->Fill(leading_momentum);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_ele_momentum_nue_cc_out_fv->Fill(leading_momentum);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_ele_momentum_nue_cc->Fill(leading_momentum);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_ele_momentum_nue_cc->Fill(leading_momentum);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_ele_momentum_nue_cc->Fill(leading_momentum);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_ele_momentum_nue_cc->Fill(leading_momentum);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_ele_momentum_numu_cc->Fill(leading_momentum);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_ele_momentum_numu_cc->Fill(leading_momentum);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_ele_momentum_numu_cc->Fill(leading_momentum);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_ele_momentum_numu_cc->Fill(leading_momentum);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_ele_momentum_numu_cc->Fill(leading_momentum);
		}
		if(tpco_id == "nc")
		{
			h_ele_momentum_nc->Fill(leading_momentum);
		}
		if(tpco_id == "nc_pi0")
		{
			h_ele_momentum_nc_pi0->Fill(leading_momentum);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_ele_momentum_nue_cc_mixed->Fill(leading_momentum);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_pfp_phi_numu_cc_mixed->Fill(leading_shower_phi);
			h_ele_momentum_numu_cc->Fill(leading_momentum);
		}
		if(tpco_id == "cosmic")
		{
			h_ele_momentum_cosmic->Fill(leading_momentum);
		}
		if(tpco_id == "other_mixed")
		{
			h_ele_momentum_other_mixed->Fill(leading_momentum);
		}
		if(tpco_id == "unmatched")
		{
			h_ele_momentum_unmatched->Fill(leading_momentum);
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsLeadingMomentumInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                        TH1D * h_ele_momentum_intime)
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
		h_ele_momentum_intime->Fill(leading_momentum);
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsLeadingMomentumThetaSlice(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                            std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                            std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                            TH1D * h_ele_momentum_1_nue_cc,
                                                            TH1D * h_ele_momentum_1_nue_cc_out_fv,
                                                            TH1D * h_ele_momentum_1_nue_cc_mixed,
                                                            TH1D * h_ele_momentum_1_numu_cc,
                                                            TH1D * h_ele_momentum_1_numu_cc_mixed,
                                                            TH1D * h_ele_momentum_1_nc,
                                                            TH1D * h_ele_momentum_1_nc_pi0,
                                                            TH1D * h_ele_momentum_1_cosmic,
                                                            TH1D * h_ele_momentum_1_other_mixed,
                                                            TH1D * h_ele_momentum_1_unmatched,
                                                            TH1D * h_ele_momentum_2_nue_cc,
                                                            TH1D * h_ele_momentum_2_nue_cc_out_fv,
                                                            TH1D * h_ele_momentum_2_nue_cc_mixed,
                                                            TH1D * h_ele_momentum_2_numu_cc,
                                                            TH1D * h_ele_momentum_2_numu_cc_mixed,
                                                            TH1D * h_ele_momentum_2_nc,
                                                            TH1D * h_ele_momentum_2_nc_pi0,
                                                            TH1D * h_ele_momentum_2_cosmic,
                                                            TH1D * h_ele_momentum_2_other_mixed,
                                                            TH1D * h_ele_momentum_2_unmatched,
                                                            TH1D * h_ele_momentum_3_nue_cc,
                                                            TH1D * h_ele_momentum_3_nue_cc_out_fv,
                                                            TH1D * h_ele_momentum_3_nue_cc_mixed,
                                                            TH1D * h_ele_momentum_3_numu_cc,
                                                            TH1D * h_ele_momentum_3_numu_cc_mixed,
                                                            TH1D * h_ele_momentum_3_nc,
                                                            TH1D * h_ele_momentum_3_nc_pi0,
                                                            TH1D * h_ele_momentum_3_cosmic,
                                                            TH1D * h_ele_momentum_3_other_mixed,
                                                            TH1D * h_ele_momentum_3_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_momentum = leading_shower.pfpMomentum();
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, 0, 0);
		const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_nue_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_nue_cc_out_fv->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_nue_cc_out_fv->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_nue_cc_out_fv->Fill(leading_momentum); }
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_nue_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_nue_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_nue_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_nue_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_nue_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "numu_cc_qe")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_numu_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "numu_cc_res")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_numu_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "numu_cc_dis")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_numu_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "numu_cc_coh")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_numu_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "numu_cc_mec")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_numu_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "nc")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_nc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_nc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_nc->Fill(leading_momentum); }
		}
		if(tpco_id == "nc_pi0")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_nc_pi0->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_nc_pi0->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_nc_pi0->Fill(leading_momentum); }
		}
		if(tpco_id == "nue_cc_mixed")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_nue_cc_mixed->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_nue_cc_mixed->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_nue_cc_mixed->Fill(leading_momentum); }
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_pfp_phi_numu_cc_mixed->Fill(leading_shower_phi);
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_numu_cc->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_numu_cc->Fill(leading_momentum); }
		}
		if(tpco_id == "cosmic")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_cosmic->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_cosmic->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_cosmic->Fill(leading_momentum); }
		}
		if(tpco_id == "other_mixed")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_other_mixed->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_other_mixed->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_other_mixed->Fill(leading_momentum); }
		}
		if(tpco_id == "unmatched")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_unmatched->Fill(leading_momentum); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_unmatched->Fill(leading_momentum); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_unmatched->Fill(leading_momentum); }
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsLeadingMomentumThetaSliceInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                                  TH1D * h_ele_momentum_1_intime,
                                                                  TH1D * h_ele_momentum_2_intime,
                                                                  TH1D * h_ele_momentum_3_intime)
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
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, 0, 0);
		const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);
		if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_ele_momentum_1_intime->Fill(leading_momentum); }
		if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_ele_momentum_2_intime->Fill(leading_momentum); }
		if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_ele_momentum_3_intime->Fill(leading_momentum); }
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsdedxThetaSlice(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                 std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                 std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                                 TH1D * h_dedx_1_nue_cc,
                                                 TH1D * h_dedx_1_nue_cc_out_fv,
                                                 TH1D * h_dedx_1_nue_cc_mixed,
                                                 TH1D * h_dedx_1_numu_cc,
                                                 TH1D * h_dedx_1_numu_cc_mixed,
                                                 TH1D * h_dedx_1_nc,
                                                 TH1D * h_dedx_1_nc_pi0,
                                                 TH1D * h_dedx_1_cosmic,
                                                 TH1D * h_dedx_1_other_mixed,
                                                 TH1D * h_dedx_1_unmatched,
                                                 TH1D * h_dedx_2_nue_cc,
                                                 TH1D * h_dedx_2_nue_cc_out_fv,
                                                 TH1D * h_dedx_2_nue_cc_mixed,
                                                 TH1D * h_dedx_2_numu_cc,
                                                 TH1D * h_dedx_2_numu_cc_mixed,
                                                 TH1D * h_dedx_2_nc,
                                                 TH1D * h_dedx_2_nc_pi0,
                                                 TH1D * h_dedx_2_cosmic,
                                                 TH1D * h_dedx_2_other_mixed,
                                                 TH1D * h_dedx_2_unmatched,
                                                 TH1D * h_dedx_3_nue_cc,
                                                 TH1D * h_dedx_3_nue_cc_out_fv,
                                                 TH1D * h_dedx_3_nue_cc_mixed,
                                                 TH1D * h_dedx_3_numu_cc,
                                                 TH1D * h_dedx_3_numu_cc_mixed,
                                                 TH1D * h_dedx_3_nc,
                                                 TH1D * h_dedx_3_nc_pi0,
                                                 TH1D * h_dedx_3_cosmic,
                                                 TH1D * h_dedx_3_other_mixed,
                                                 TH1D * h_dedx_3_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		std::string tpco_id     = tpco_classifier_v->at(i).first;
		const int leading_index = tpco_classifier_v->at(i).second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, 0, 0);
		const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);

		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_nue_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_nue_cc_out_fv->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_nue_cc_out_fv->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_nue_cc_out_fv->Fill(leading_dedx); }
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_nue_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_nue_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_nue_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_nue_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_nue_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "numu_cc_qe")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_numu_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "numu_cc_res")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_numu_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "numu_cc_dis")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_numu_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "numu_cc_coh")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_numu_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "numu_cc_mec")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_numu_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "nc")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_nc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_nc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_nc->Fill(leading_dedx); }
		}
		if(tpco_id == "nc_pi0")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_nc_pi0->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_nc_pi0->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_nc_pi0->Fill(leading_dedx); }
		}
		if(tpco_id == "nue_cc_mixed")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_nue_cc_mixed->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_nue_cc_mixed->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_nue_cc_mixed->Fill(leading_dedx); }
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_ele_pfp_phi_numu_cc_mixed->Fill(leading_shower_phi);
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_numu_cc->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_numu_cc->Fill(leading_dedx); }
		}
		if(tpco_id == "cosmic")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_cosmic->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_cosmic->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_cosmic->Fill(leading_dedx); }
		}
		if(tpco_id == "other_mixed")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_other_mixed->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_other_mixed->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_other_mixed->Fill(leading_dedx); }
		}
		if(tpco_id == "unmatched")
		{
			if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_unmatched->Fill(leading_dedx); }
			if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_unmatched->Fill(leading_dedx); }
			if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_unmatched->Fill(leading_dedx); }
		}
	}//end pfp loop
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutsdedxThetaSliceInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                                       TH1D * h_dedx_1_intime,
                                                       TH1D * h_dedx_2_intime,
                                                       TH1D * h_dedx_3_intime)
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
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_shower_z = leading_shower.pfpDirZ();
		const double leading_shower_y = leading_shower.pfpDirY();
		const double leading_shower_x = leading_shower.pfpDirX();
		TVector3 shower_vector(leading_shower_x, leading_shower_y, leading_shower_z);
		TVector3 numi_vector;
		numi_vector.SetMagThetaPhi(1, 0, 0);
		const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);
		if(leading_shower_theta >= 0 && leading_shower_theta < 40)    {h_dedx_1_intime->Fill(leading_dedx * (242.72 / 196.979)); }
		if(leading_shower_theta >= 40 && leading_shower_theta < 90)   {h_dedx_2_intime->Fill(leading_dedx * (242.72 / 196.979)); }
		if(leading_shower_theta >= 90 && leading_shower_theta <= 180) {h_dedx_3_intime->Fill(leading_dedx * (242.72 / 196.979)); }
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::dEdxTheta(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                    std::vector<std::pair<std::string, int> > * tpco_classifier_v,
                                    TH2D * h_dedx_theta_nue_cc,
                                    TH2D * h_dedx_theta_nue_cc_mixed,
                                    TH2D * h_dedx_theta_nue_cc_out_fv,
                                    TH2D * h_dedx_theta_numu_cc,
                                    TH2D * h_dedx_theta_nc,
                                    TH2D * h_dedx_theta_cosmic,
                                    TH2D * h_dedx_theta_nc_pi0,
                                    TH2D * h_dedx_theta_numu_cc_mixed,
                                    TH2D * h_dedx_theta_other_mixed,
                                    TH2D * h_dedx_theta_unmatched)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{

		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int leading_index   = tpco_classifier_v->at(i).second;
		std::string tpco_id = tpco_classifier_v->at(i).first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_shower_theta = acos(leading_shower.pfpDirZ()) * 180 / 3.1415;
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_bar_cc_qe")
		{
			h_dedx_theta_nue_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_dedx_theta_nue_cc_out_fv->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res")
		{
			h_dedx_theta_nue_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis")
		{
			h_dedx_theta_nue_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh")
		{
			h_dedx_theta_nue_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec")
		{
			h_dedx_theta_nue_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_dedx_theta_numu_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_dedx_theta_numu_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_dedx_theta_numu_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_dedx_theta_numu_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_dedx_theta_numu_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "nc")
		{
			h_dedx_theta_nc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "nc_pi0")
		{
			h_dedx_theta_nc_pi0->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_dedx_theta_nue_cc_mixed->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			//h_dedx_cuts_numu_cc_mixed->Fill(leading_dedx);
			h_dedx_theta_numu_cc->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "cosmic")
		{
			h_dedx_theta_cosmic->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "other_mixed")
		{
			h_dedx_theta_other_mixed->Fill(leading_dedx, leading_shower_theta);
		}
		if(tpco_id == "unmatched")
		{
			h_dedx_theta_unmatched->Fill(leading_dedx, leading_shower_theta);
		}
	}        //end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions::dEdxThetaInTime(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                          TH2D * h_dedx_theta_intime)
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
		h_dedx_theta_intime->Fill(leading_dedx * (242.72 / 196.979), leading_shower_theta);
	}
}
//***************************************************************************
//***************************************************************************
//end functions
