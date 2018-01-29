#include "selection_functions.h"
#include "selection_cuts.h"

//***************************************************************************
void selection_functions::PostCutsdEdx(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose,
                                       bool has_pi0,
                                       double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                       double vtxX, double vtxY, double vtxZ,
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
		const double tpc_vtx_x = tpc_obj.pfpVtxX();
		const double tpc_vtx_y = tpc_obj.pfpVtxY();
		const double tpc_vtx_z = tpc_obj.pfpVtxZ();
		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();
		//loop over pfparticles in the TPCO
		int most_hits = 0;
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		if(tpco_id == "nue_cc_qe")
		{
			h_dedx_cuts_nue_cc->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_dedx_cuts_nue_cc_out_fv->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_res")
		{
			h_dedx_cuts_nue_cc->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_dis")
		{
			h_dedx_cuts_nue_cc->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_coh")
		{
			h_dedx_cuts_nue_cc->Fill(leading_dedx);
		}
		if(tpco_id == "nue_cc_mec")
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
			h_dedx_cuts_numu_cc_mixed->Fill(leading_dedx);
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
void selection_functions::FillPostCutVector(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                            std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                                            double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                            double vtxX, double vtxY, double vtxZ,
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

		const int tpc_obj_mode = tpc_obj.Mode();
		const int n_pfp = tpc_obj.NumPFParticles();

		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
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
void selection_functions::PrintPostCutVector(std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v,
                                             bool _post_cuts_verbose)
{
	const int passed_events = post_cuts_v->size();
	std::cout << "* * * * * * * * * * * * * * * * *" << std::endl;
	std::cout << "Total Passed Events: " << passed_events << std::endl;
	std::cout << "* * * * * * * * * * * * * * * * *" << std::endl;
	for(auto const my_tuple: * post_cuts_v)
	{
		const int event_num = std::get<0>(my_tuple);
		const int run_num = std::get<1>(my_tuple);
		const double pfp_vtx_x = std::get<2>(my_tuple);
		const double pfp_vtx_y = std::get<3>(my_tuple);
		const double pfp_vtx_z = std::get<4>(my_tuple);
		const std::string reason = std::get<5>(my_tuple);
		const std::string event_type = std::get<6>(my_tuple);
		const int num_tracks = std::get<7>(my_tuple);
		const int num_showers = std::get<8>(my_tuple);
		const double opening_angle = std::get<9>(my_tuple);
		std::cout << "* * * * * * * * * * * * * * * * *" << std::endl;
		std::cout << "Event Type     : " << event_type << std::endl;
		std::cout << "Event Number   : " << event_num << std::endl;
		std::cout << "Run Number     : " << run_num << std::endl;
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
void selection_functions::PostCutVectorPlots(std::vector<std::tuple<int, int, double, double, double, std::string, std::string, int, int, double> > * post_cuts_v,
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
		const double pfp_vtx_x = std::get<2>(my_tuple);
		const double pfp_vtx_y = std::get<3>(my_tuple);
		const double pfp_vtx_z = std::get<4>(my_tuple);
		const std::string reason = std::get<5>(my_tuple);
		const std::string event_type = std::get<6>(my_tuple);
		const int num_tracks = std::get<7>(my_tuple);
		const int num_showers = std::get<8>(my_tuple);
		const double opening_angle = std::get<9>(my_tuple);

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

	std::cout << "4+ Shower Background Events: " << bkg_events_4 << std::endl;

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
std::vector<int> selection_functions::TabulateOrigins(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                                                      double _x1, double _x2, double _y1, double _y2,
                                                      double _z1, double _z2, double vtxX, double vtxY, double vtxZ)
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
	std::vector<int> tabulated_origins;
	tabulated_origins.resize(22);

	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		//loop over pfparticles in the TPCO
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;

		if(tpco_id == "nue_cc_qe")     {nue_cc_qe++; }
		if(tpco_id == "nue_cc_out_fv") {nue_cc_out_fv++; }
		if(tpco_id == "nue_cc_res")    {nue_cc_res++; }
		if(tpco_id == "nue_cc_coh")    {nue_cc_coh++; }
		if(tpco_id == "nue_cc_dis")    {nue_cc_dis++; }
		if(tpco_id == "nue_cc_mec")    {nue_cc_mec++; }
		if(tpco_id == "nue_cc_mixed")  {nue_cc_mixed++; }
		if(tpco_id == "nc")            {nc++; }
		if(tpco_id == "numu_cc_qe")    {numu_cc_qe++; }
		if(tpco_id == "numu_cc_res")   {numu_cc_res++; }
		if(tpco_id == "numu_cc_coh")   {numu_cc_coh++; }
		if(tpco_id == "numu_cc_dis")   {numu_cc_dis++; }
		if(tpco_id == "numu_cc_mec")   {numu_cc_mec++; }
		if(tpco_id == "numu_cc_mixed") {numu_cc_mixed++; }
		if(tpco_id == "nc_pi0")        {nc_pi0++; }
		if(tpco_id == "cosmic")        {cosmic++; }
		if(tpco_id == "other_mixed")   {other_mixed++; }
		if(tpco_id == "unmatched")     {unmatched++; }
	}//end loop tpc objects

	nue_cc = nue_cc_qe + nue_cc_res + nue_cc_dis + nue_cc_coh + nue_cc_mec;
	numu_cc = numu_cc_qe + numu_cc_res + numu_cc_dis + numu_cc_coh + numu_cc_mec;
	total = nue_cc + nue_cc_mixed + nue_cc_out_fv + cosmic + nc + numu_cc + numu_cc_mixed + nc_pi0 + unmatched + other_mixed;

	tabulated_origins.at(0)  = nue_cc;
	tabulated_origins.at(1)  = nue_cc_mixed;
	tabulated_origins.at(2)  = cosmic;
	tabulated_origins.at(3)  = nc;
	tabulated_origins.at(4)  = numu_cc;
	tabulated_origins.at(5)  = unmatched;
	tabulated_origins.at(6)  = other_mixed;
	tabulated_origins.at(7)  = total;
	tabulated_origins.at(8)  = signal_tpco_num;
	tabulated_origins.at(9)  = nue_cc_out_fv;
	tabulated_origins.at(10) = nc_pi0;
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
void selection_functions::TotalOrigins(std::vector<int> tabulated_origins, std::vector<int> * total_cut_origins)
{
	total_cut_origins->at(0)  += tabulated_origins.at(0);
	total_cut_origins->at(1)  += tabulated_origins.at(1);
	total_cut_origins->at(9)  += tabulated_origins.at(9);
	total_cut_origins->at(2)  += tabulated_origins.at(2);
	total_cut_origins->at(3)  += tabulated_origins.at(3);
	total_cut_origins->at(4)  += tabulated_origins.at(4);
	total_cut_origins->at(11) += tabulated_origins.at(11);
	total_cut_origins->at(10) += tabulated_origins.at(10);
	total_cut_origins->at(5)  += tabulated_origins.at(5);
	total_cut_origins->at(6)  += tabulated_origins.at(6);
	total_cut_origins->at(7)  += tabulated_origins.at(7);
	total_cut_origins->at(12) += tabulated_origins.at(12);
	total_cut_origins->at(13) += tabulated_origins.at(13);
	total_cut_origins->at(14) += tabulated_origins.at(14);
	total_cut_origins->at(15) += tabulated_origins.at(15);
	total_cut_origins->at(16) += tabulated_origins.at(16);
	total_cut_origins->at(17) += tabulated_origins.at(17);
	total_cut_origins->at(18) += tabulated_origins.at(18);
	total_cut_origins->at(19) += tabulated_origins.at(19);
	total_cut_origins->at(20) += tabulated_origins.at(20);
	total_cut_origins->at(21) += tabulated_origins.at(21);
}
//***************************************************************************
//***************************************************************************
//modify this so it takes a string of the cut name so I only pass it a few variable at a time,
//then I can call this function several times later at the bottom
void selection_functions::PrintInfo(int mc_nue_cc_counter,
                                    std::vector<int> * counter_v,
                                    std::string cut_name)
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

	std::cout << " <" << cut_name << "> " << std::endl;
	std::cout << " Total Candidate Nue     : " << counter << std::endl;
	std::cout << " Number of Nue CC        : " << counter_nue_cc << std::endl;
	std::cout << " Number of Nue CC Mixed  : " << counter_nue_cc_mixed << std::endl;
	std::cout << " Number of Nue CC out FV : " << counter_nue_cc_out_fv << std::endl;
	std::cout << " Number of Cosmic        : " << counter_cosmic << std::endl;
	std::cout << " Number of Numu CC       : " << counter_numu_cc << std::endl;
	std::cout << " Number of Numu CC Mixed : " << counter_numu_cc_mixed << std::endl;
	std::cout << " Number of NC            : " << counter_nc << std::endl;
	std::cout << " Number of NC Pi0        : " << counter_nc_pi0 << std::endl;
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
std::pair<std::string, int> selection_functions::TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool has_pi0,
                                                                 double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                                                 double vtxX, double vtxY, double vtxZ)
{
	int part_nue_cc    = 0;
	int part_cosmic    = 0;
	int part_nc        = 0;
	int part_nc_pi0    = 0;
	int part_numu_cc   = 0;
	int part_unmatched = 0;
	bool true_in_tpc = false;
	selection_cuts _functions_instance;

	const int tpc_obj_mode = tpc_obj.Mode();
	const int n_pfp = tpc_obj.NumPFParticles();
	int most_hits = 0;
	int leading_index = 0;
	for(int j = 0; j < n_pfp; j++)
	{
		auto const part = tpc_obj.GetParticle(j);
		const int n_pfp_hits = part.NumPFPHits();
		if(n_pfp_hits > most_hits) {leading_index = j; most_hits = n_pfp_hits; }
		if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && (part.MCParentPdg() == 12 || part.MCParentPdg() == -12)) { part_nue_cc++; }
		if(part.Origin() == "kBeamNeutrino" && part.CCNC() == 0 && (part.MCParentPdg() == 14 || part.MCParentPdg() == -14)) { part_numu_cc++; }
		if(part.CCNC() == 1 && part.Origin() == "kBeamNeutrino")
		{
			if(has_pi0 == true)  {part_nc_pi0++; }
			if(has_pi0 == false) {part_nc++; }
		}
		if(part.Origin() == "kCosmicRay") { part_cosmic++; }
		if(part.Origin() == "kUnknown")   { part_unmatched++; }
	}

	//now to catagorise the tpco
	if(part_cosmic > 0)
	{
		if(part_nue_cc  > 0 )                        { return std::make_pair("nue_cc_mixed", leading_index);  }
		if(part_numu_cc > 0 )                        { return std::make_pair("numu_cc_mixed", leading_index); }
		if(part_nc  > 0 || part_nc_pi0 > 0)          { return std::make_pair("other_mixed", leading_index);   }
		return std::make_pair("cosmic", leading_index);
	}
	//this uses the true neutrino vertex for this specific event
	//not the true vtx per tpc object - maybe this can be fixed in the future...
	//but using the true nu vtx only matters for the pure signal events,
	//where the neutrino vertex IS the true tpc object vertex
	true_in_tpc = _functions_instance.selection_cuts::in_fv(vtxX, vtxY, vtxZ,
	                                                        _x1, _x2, _y1,
	                                                        _y2, _z1, _z2);
	if(part_cosmic == 0)
	{
		if(part_nue_cc    > 0 && true_in_tpc == false) { return std::make_pair("nue_cc_out_fv", leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 0   ) { return std::make_pair("nue_cc_qe",     leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 1   ) { return std::make_pair("nue_cc_res",    leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 2   ) { return std::make_pair("nue_cc_dis",    leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 3   ) { return std::make_pair("nue_cc_coh",    leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 10  ) { return std::make_pair("nue_cc_mec",    leading_index);   }
		if(part_numu_cc   > 0 && tpc_obj_mode == 0   ) { return std::make_pair("numu_cc_qe",    leading_index);   }
		if(part_numu_cc   > 0 && tpc_obj_mode == 1   ) { return std::make_pair("numu_cc_res",   leading_index);   }
		if(part_numu_cc   > 0 && tpc_obj_mode == 2   ) { return std::make_pair("numu_cc_dis",   leading_index);   }
		if(part_numu_cc   > 0 && tpc_obj_mode == 3   ) { return std::make_pair("numu_cc_coh",   leading_index);   }
		if(part_numu_cc   > 0 && tpc_obj_mode == 10  ) { return std::make_pair("numu_cc_mec",   leading_index);   }
		if(part_nc        > 0                        ) { return std::make_pair("nc",            leading_index);   }
		if(part_nc_pi0    > 0                        ) { return std::make_pair("nc_pi0",        leading_index);   }
		if(part_unmatched > 0                        ) { return std::make_pair("unmatched",     leading_index);   }
	}
	//return the string for the tpco id
}//end function
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
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutOpenAngle(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                           double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                                           TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_mixed,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);

		if(tpco_id == "nue_cc_mixed")  {h_leading_shower_open_angle_nue_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "numu_cc_mixed") {h_leading_shower_open_angle_numu_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "other_mixed")   {h_leading_shower_open_angle_other_mixed->Fill(leading_open_angle); }
		if(tpco_id == "cosmic")        {h_leading_shower_open_angle_cosmic->Fill(leading_open_angle); }
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_mec")
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
void selection_functions::PostCutOpenAngle1Shower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                  std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                                  double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                                                  TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_mixed,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);

		if(tpco_id == "nue_cc_mixed")  {h_leading_shower_open_angle_nue_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "numu_cc_mixed") {h_leading_shower_open_angle_numu_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "other_mixed")   {h_leading_shower_open_angle_other_mixed->Fill(leading_open_angle); }
		if(tpco_id == "cosmic")        {h_leading_shower_open_angle_cosmic->Fill(leading_open_angle); }
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_mec")
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
void selection_functions::PostCutOpenAngle2PlusShower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                                      double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                                                      TH1D * h_leading_shower_open_angle_nue_cc, TH1D * h_leading_shower_open_angle_nue_cc_mixed,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);

		if(tpco_id == "nue_cc_mixed")  {h_leading_shower_open_angle_nue_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "numu_cc_mixed") {h_leading_shower_open_angle_numu_cc_mixed->Fill(leading_open_angle); }
		if(tpco_id == "other_mixed")   {h_leading_shower_open_angle_other_mixed->Fill(leading_open_angle); }
		if(tpco_id == "cosmic")        {h_leading_shower_open_angle_cosmic->Fill(leading_open_angle); }
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_mec")
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
void selection_functions::PostCutTrkVtx(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                        double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                                        TH1D * h_trk_vtx_dist_nue_cc, TH1D * h_trk_vtx_dist_nue_cc_mixed,
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
				//std::cout << trk_length << std::endl;
			}
		}//end pfp loop
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0,
		                                                         _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		if(tpco_id == "nue_cc_mixed")
		{
			if(has_track) {h_trk_vtx_dist_nue_cc_mixed->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "numu_cc_mixed")
		{
			if(has_track) {h_trk_vtx_dist_numu_cc_mixed->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "other_mixed")
		{
			if(has_track) {h_trk_vtx_dist_other_mixed->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "cosmic")
		{
			if(has_track) {h_trk_vtx_dist_cosmic->Fill(smallest_trk_vtx_dist); }
		}
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_mec")
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
void selection_functions::TopologyPlots1(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                         std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                         double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
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

		if(tpco_id == "nue_cc_qe")
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
		if(tpco_id == "nue_cc_res")
		{
			h_pfp_track_shower_nue_cc_res->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_res->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_res->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_res->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_dis")
		{
			h_pfp_track_shower_nue_cc_dis->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_dis->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_dis->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_dis->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_coh")
		{
			h_pfp_track_shower_nue_cc_coh->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_coh->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_coh->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_coh->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_mec")
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
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                                              double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                              double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		if(tpco_id == "nue_cc_qe")      {h_pfp_shower_open_angle_nue_cc_qe->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_out_fv")  {h_pfp_shower_open_angle_nue_cc_out_fv->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_res")     {h_pfp_shower_open_angle_nue_cc_res->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_dis")     {h_pfp_shower_open_angle_nue_cc_dis->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_coh")     {h_pfp_shower_open_angle_nue_cc_coh->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_mec")     {h_pfp_shower_open_angle_nue_cc_mec->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_qe")     {h_pfp_shower_open_angle_numu_cc_qe->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_res")    {h_pfp_shower_open_angle_numu_cc_res->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_dis")    {h_pfp_shower_open_angle_numu_cc_dis->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_coh")    {h_pfp_shower_open_angle_numu_cc_coh->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_mec")    {h_pfp_shower_open_angle_numu_cc_mec->Fill(n_pfp_showers); }
		if(tpco_id == "nc")             {h_pfp_shower_open_angle_nc->Fill(n_pfp_showers); }
		if(tpco_id == "nc_pi0")         {h_pfp_shower_open_angle_nc_pi0->Fill(n_pfp_showers); }
		if(tpco_id == "nue_cc_mixed")   {h_pfp_shower_open_angle_nue_cc_mixed->Fill(n_pfp_showers); }
		if(tpco_id == "numu_cc_mixed")  {h_pfp_shower_open_angle_numu_cc_mixed->Fill(n_pfp_showers); }
		if(tpco_id == "cosmic")         {h_pfp_shower_open_angle_cosmic->Fill(n_pfp_showers); }
		if(tpco_id == "other_mixed")    {h_pfp_shower_open_angle_other_mixed->Fill(n_pfp_showers); }
		if(tpco_id == "unmatched")      {h_pfp_shower_open_angle_unmatched->Fill(n_pfp_showers); }
	}        //end loop tpc objects

}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::TopologyPlots2(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                         std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0,
                                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                         double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
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

		if(tpco_id == "nue_cc_qe")
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
		if(tpco_id == "nue_cc_res")
		{
			h_pfp_track_shower_nue_cc_res->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_res->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_res->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_res->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_dis")
		{
			h_pfp_track_shower_nue_cc_dis->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_dis->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_dis->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_dis->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_coh")
		{
			h_pfp_track_shower_nue_cc_coh->Fill(n_pfp_tracks, n_pfp_showers);
			h_pfp_track_nue_cc_coh->Fill(n_pfp_tracks);
			h_pfp_shower_nue_cc_coh->Fill(n_pfp_showers);
			h_leading_shower_mc_pdg_nue_cc_coh->Fill(leading_origin_int, leading_pdg_int);
		}
		if(tpco_id == "nue_cc_mec")
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
                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                           double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                                           TH1D * h_vtx_flash_nue_cc, TH1D * h_vtx_flash_nue_cc_mixed,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;

		if(tpco_id == "nue_cc_qe")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		//if(tpco_id == "nue_cc_out_fv")
		//{
		//	h_vtx_flash_nue_cc_out_fv->Fill(distance);
		//}
		if(tpco_id == "nue_cc_res")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_dis")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_coh")
		{
			h_vtx_flash_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_mec")
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
			h_vtx_flash_numu_cc_mixed->Fill(distance);
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
void selection_functions::PostCutsShwrVtx(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                          double _x1, double _x2, double _y1, double _y2, double _z1, double _z2, double vtxX, double vtxY, double vtxZ,
                                          TH1D * h_shwr_vtx_dist_nue_cc,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_vtx_x = leading_shower.pfpVtxX();
		const double leading_vtx_y = leading_shower.pfpVtxY();
		const double leading_vtx_z = leading_shower.pfpVtxZ();
		const double distance = sqrt(pow((tpc_vtx_x - leading_vtx_x), 2) +
		                             pow((tpc_vtx_y - leading_vtx_y), 2) +
		                             pow((tpc_vtx_z - leading_vtx_z), 2));

		if(tpco_id == "nue_cc_qe")
		{
			h_shwr_vtx_dist_nue_cc->Fill(distance);
		}
		//if(tpco_id == "nue_cc_out_fv")
		//{
		//	h_vtx_flash_nue_cc_out_fv->Fill(distance);
		//}
		if(tpco_id == "nue_cc_res")
		{
			h_shwr_vtx_dist_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_dis")
		{
			h_shwr_vtx_dist_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_coh")
		{
			h_shwr_vtx_dist_nue_cc->Fill(distance);
		}
		if(tpco_id == "nue_cc_mec")
		{
			h_shwr_vtx_dist_nue_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_qe")
		{
			h_shwr_vtx_dist_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_res")
		{
			h_shwr_vtx_dist_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_dis")
		{
			h_shwr_vtx_dist_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_coh")
		{
			h_shwr_vtx_dist_numu_cc->Fill(distance);
		}
		if(tpco_id == "numu_cc_mec")
		{
			h_shwr_vtx_dist_numu_cc->Fill(distance);
		}
		if(tpco_id == "nc")
		{
			h_shwr_vtx_dist_nc->Fill(distance);
		}
		if(tpco_id == "nc_pi0")
		{
			h_shwr_vtx_dist_nc_pi0->Fill(distance);
		}
		if(tpco_id == "nue_cc_mixed")
		{
			h_shwr_vtx_dist_nue_cc_mixed->Fill(distance);
		}
		if(tpco_id == "numu_cc_mixed")
		{
			h_shwr_vtx_dist_numu_cc_mixed->Fill(distance);
		}
		if(tpco_id == "cosmic")
		{
			h_shwr_vtx_dist_cosmic->Fill(distance);
		}
		if(tpco_id == "other_mixed")
		{
			h_shwr_vtx_dist_other_mixed->Fill(distance);
		}
		if(tpco_id == "unmatched")
		{
			h_shwr_vtx_dist_unmatched->Fill(distance);
		}
	}        //end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions::PostCutHitThreshold(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                              double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                              double vtxX, double vtxY, double vtxZ,
                                              double mc_nu_energy, double mc_ele_energy,
                                              TH2D * h_shwr_hits_nu_eng, TH2D * h_shwr_hits_ele_eng)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		int leading_index = tpco_class.second;
		std::string tpco_id = tpco_class.first;

		int pfp_shower_hits = 0;
		if(tpco_id == "nue_cc_qe"  ||
		   tpco_id == "nue_cc_res" ||
		   tpco_id == "nue_cc_dis" ||
		   tpco_id == "nue_cc_mec" ||
		   tpco_id == "nue_cc_coh")
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
                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                             double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                             double vtxX, double vtxY, double vtxZ,
                                             std::vector<int> * no_track, std::vector<int> * has_track)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		//signal
		if(n_pfp_tracks == 0)
		{
			if(tpco_id == "nue_cc_qe"  ||
			   tpco_id == "nue_cc_res" ||
			   tpco_id == "nue_cc_dis" ||
			   tpco_id == "nue_cc_mec" ||
			   tpco_id == "nue_cc_coh")
			{
				no_track->at(0) += 1;
			}
			//not signal
			if(tpco_id != "nue_cc_qe"  &&
			   tpco_id != "nue_cc_res" &&
			   tpco_id != "nue_cc_dis" &&
			   tpco_id != "nue_cc_mec" &&
			   tpco_id != "nue_cc_coh")
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
			   tpco_id == "nue_cc_coh")
			{
				has_track->at(0) += 1;
			}
			//not signal
			if(tpco_id != "nue_cc_qe"  &&
			   tpco_id != "nue_cc_res" &&
			   tpco_id != "nue_cc_dis" &&
			   tpco_id != "nue_cc_mec" &&
			   tpco_id != "nue_cc_coh")
			{
				has_track->at(1) += 1;
			}
		}
	}//end loop tpc objects
}
//***************************************************************************
//***************************************************************************
void selection_functions::SequentialTrueEnergyPlots(int mc_nu_id, double mc_nu_vtx_x, double mc_nu_vtx_y, double mc_nu_vtx_z,
                                                    double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                                    std::vector<int> tabulated_origins, double mc_nu_energy,
                                                    double mc_ele_energy, TH1D * h_selected_nu_energy, TH1D * h_selected_ele_energy)
{
	//this checks if there is a true nue/nue-bar CC event and a selected nue_cc signal event, true in FV
	selection_cuts _functions_instance;
	if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins.at(0) == 1) {
		if(_functions_instance.selection_cuts::in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, _x1, _x2, _y1, _y2, _z1, _z2) == true) {
			h_selected_nu_energy->Fill(mc_nu_energy);
			h_selected_ele_energy->Fill(mc_ele_energy);
		}
	}
}
//***************************************************************************
//***************************************************************************
void selection_functions::ChargeShare(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                      double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                      double vtxX, double vtxY, double vtxZ, TH1D * h_charge_share_nue_cc_mixed)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		const int n_pfp_tracks = tpc_obj.NPfpTracks();
		const int n_pfp_showers = tpc_obj.NPfpShowers();
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
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
void selection_functions::FlashTot0(std::vector< double> largest_flash_v, double mc_nu_time, int mc_nu_id, std::vector<int> tabulated_origins,
                                    double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                    double vtxX, double vtxY, double vtxZ, TH1D * h_flash_t0_diff)
{
	selection_cuts _functions_instance;
	if((mc_nu_id == 1 || mc_nu_id == 5) && tabulated_origins.at(0) == 1)
	{
		const bool InFV = _functions_instance.selection_cuts::in_fv(vtxX, vtxY, vtxZ, _x1, _x2, _y1, _y2, _z1, _z2);
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
                                          std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                          double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                          double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
		const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);
		if(tpco_id == "nue_cc_qe")
		{
			h_dedx_open_angle_nue_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_dedx_open_angle_nue_cc_out_fv->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_res")
		{
			h_dedx_open_angle_nue_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_dis")
		{
			h_dedx_open_angle_nue_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_coh")
		{
			h_dedx_open_angle_nue_cc->Fill(leading_dedx, leading_open_angle);
		}
		if(tpco_id == "nue_cc_mec")
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
			h_dedx_open_angle_numu_cc_mixed->Fill(leading_dedx, leading_open_angle);
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
//shower hits vs shower length, to see if we can better use the hit threshold cut to remove less signal?
void selection_functions::ShowerLengthvsHits(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                             double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                             double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int leading_hits = leading_shower.NumPFPHits();
		const double leading_length = leading_shower.pfpLength();
		if(tpco_id == "nue_cc_qe")
		{
			h_shwr_len_hits_nue_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_shwr_len_hits_nue_cc_out_fv->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_res")
		{
			h_shwr_len_hits_nue_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_dis")
		{
			h_shwr_len_hits_nue_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_coh")
		{
			h_shwr_len_hits_nue_cc->Fill(leading_length, leading_hits);
		}
		if(tpco_id == "nue_cc_mec")
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
			h_shwr_len_hits_numu_cc_mixed->Fill(leading_length, leading_hits);
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
void selection_functions::SecondaryShowersDist(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                               std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                               double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                               double vtxX, double vtxY, double vtxZ,
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
		if(n_pfp_showers <= 3) {continue; }
		const double tpco_vtx_x = tpc_obj.pfpVtxX();
		const double tpco_vtx_y = tpc_obj.pfpVtxY();
		const double tpco_vtx_z = tpc_obj.pfpVtxZ();
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		int leading_index = tpco_class.second;
		//auto const leading_shower = tpc_obj.GetParticle(leading_index);
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
				if(tpco_id == "nue_cc_qe")
				{
					h_second_shwr_dist_nue_cc->Fill(distance);
				}
				if(tpco_id == "nue_cc_out_fv")
				{
					h_second_shwr_dist_nue_cc_out_fv->Fill(distance);
				}
				if(tpco_id == "nue_cc_res")
				{
					h_second_shwr_dist_nue_cc->Fill(distance);
				}
				if(tpco_id == "nue_cc_dis")
				{
					h_second_shwr_dist_nue_cc->Fill(distance);
				}
				if(tpco_id == "nue_cc_coh")
				{
					h_second_shwr_dist_nue_cc->Fill(distance);
				}
				if(tpco_id == "nue_cc_mec")
				{
					h_second_shwr_dist_nue_cc->Fill(distance);
				}
				if(tpco_id == "numu_cc_qe")
				{
					h_second_shwr_dist_numu_cc->Fill(distance);
				}
				if(tpco_id == "numu_cc_res")
				{
					h_second_shwr_dist_numu_cc->Fill(distance);
				}
				if(tpco_id == "numu_cc_dis")
				{
					h_second_shwr_dist_numu_cc->Fill(distance);
				}
				if(tpco_id == "numu_cc_coh")
				{
					h_second_shwr_dist_numu_cc->Fill(distance);
				}
				if(tpco_id == "numu_cc_mec")
				{
					h_second_shwr_dist_numu_cc->Fill(distance);
				}
				if(tpco_id == "nc")
				{
					h_second_shwr_dist_nc->Fill(distance);
				}
				if(tpco_id == "nc_pi0")
				{
					h_second_shwr_dist_nc_pi0->Fill(distance);
				}
				if(tpco_id == "nue_cc_mixed")
				{
					h_second_shwr_dist_nue_cc_mixed->Fill(distance);
				}
				if(tpco_id == "numu_cc_mixed")
				{
					h_second_shwr_dist_numu_cc_mixed->Fill(distance);
				}
				if(tpco_id == "cosmic")
				{
					h_second_shwr_dist_cosmic->Fill(distance);
				}
				if(tpco_id == "other_mixed")
				{
					h_second_shwr_dist_other_mixed->Fill(distance);
				}
				if(tpco_id == "unmatched")
				{
					h_second_shwr_dist_unmatched->Fill(distance);
				}
			}        //end if reco shower
		}        //end loop pfp
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::HitLengthRatio(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                         std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                         double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                         double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int pfp_pdg = leading_shower.PFParticlePdgCode();
		const double pfp_hits = leading_shower.NumPFPHits();
		const double pfp_length = leading_shower.pfpLength();
		const double pfp_hits_length_ratio = (pfp_hits / pfp_length);
		if(pfp_pdg == 11)
		{
			if(tpco_id == "nue_cc_qe")
			{
				h_hit_length_ratio_nue_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_out_fv")
			{
				h_hit_length_ratio_nue_cc_out_fv->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_res")
			{
				h_hit_length_ratio_nue_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_dis")
			{
				h_hit_length_ratio_nue_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_coh")
			{
				h_hit_length_ratio_nue_cc->Fill(pfp_hits_length_ratio);
			}
			if(tpco_id == "nue_cc_mec")
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
				h_hit_length_ratio_numu_cc_mixed->Fill(pfp_hits_length_ratio);
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
void selection_functions::TrackLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                      double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                      double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;

		for(const double trk_length : trk_length_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_out_fv")  {h_trk_length_nue_cc_out_fv->Fill(trk_length); }
			if(tpco_id == "nue_cc_res")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_dis")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_coh")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_mec")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_qe")     {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_res")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_dis")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_coh")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_mec")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "nc")             {h_trk_length_nc->Fill(trk_length); }
			if(tpco_id == "nc_pi0")         {h_trk_length_nc_pi0->Fill(trk_length); }
			if(tpco_id == "nue_cc_mixed")   {h_trk_length_nue_cc_mixed->Fill(trk_length); }
			if(tpco_id == "numu_cc_mixed")  {h_trk_length_numu_cc_mixed->Fill(trk_length); }
			if(tpco_id == "cosmic")         {h_trk_length_cosmic->Fill(trk_length); }
			if(tpco_id == "other_mixed")    {h_trk_length_other_mixed->Fill(trk_length); }
			if(tpco_id == "unmatched")      {h_trk_length_unmatched->Fill(trk_length); }
		}
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::LongestTrackLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                             std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                             double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                             double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;

		for(const double trk_length : trk_length_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_out_fv")  {h_trk_length_nue_cc_out_fv->Fill(trk_length); }
			if(tpco_id == "nue_cc_res")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_dis")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_coh")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "nue_cc_mec")     {h_trk_length_nue_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_qe")     {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_res")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_dis")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_coh")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "numu_cc_mec")    {h_trk_length_numu_cc->Fill(trk_length); }
			if(tpco_id == "nc")             {h_trk_length_nc->Fill(trk_length); }
			if(tpco_id == "nc_pi0")         {h_trk_length_nc_pi0->Fill(trk_length); }
			if(tpco_id == "nue_cc_mixed")   {h_trk_length_nue_cc_mixed->Fill(trk_length); }
			if(tpco_id == "numu_cc_mixed")  {h_trk_length_numu_cc_mixed->Fill(trk_length); }
			if(tpco_id == "cosmic")         {h_trk_length_cosmic->Fill(trk_length); }
			if(tpco_id == "other_mixed")    {h_trk_length_other_mixed->Fill(trk_length); }
			if(tpco_id == "unmatched")      {h_trk_length_unmatched->Fill(trk_length); }
		}
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::ShowerLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                       std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                       double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                       double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;

		for(const double shwr_length : shwr_length_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_out_fv")  {h_shwr_length_nue_cc_out_fv->Fill(shwr_length); }
			if(tpco_id == "nue_cc_res")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_dis")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_coh")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_mec")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_qe")     {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_res")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_dis")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_coh")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_mec")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "nc")             {h_shwr_length_nc->Fill(shwr_length); }
			if(tpco_id == "nc_pi0")         {h_shwr_length_nc_pi0->Fill(shwr_length); }
			if(tpco_id == "nue_cc_mixed")   {h_shwr_length_nue_cc_mixed->Fill(shwr_length); }
			if(tpco_id == "numu_cc_mixed")  {h_shwr_length_numu_cc_mixed->Fill(shwr_length); }
			if(tpco_id == "cosmic")         {h_shwr_length_cosmic->Fill(shwr_length); }
			if(tpco_id == "other_mixed")    {h_shwr_length_other_mixed->Fill(shwr_length); }
			if(tpco_id == "unmatched")      {h_shwr_length_unmatched->Fill(shwr_length); }
		}
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::LongestShowerLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                              double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                              double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;

		for(const double shwr_length : shwr_length_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_out_fv")  {h_shwr_length_nue_cc_out_fv->Fill(shwr_length); }
			if(tpco_id == "nue_cc_res")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_dis")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_coh")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "nue_cc_mec")     {h_shwr_length_nue_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_qe")     {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_res")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_dis")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_coh")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "numu_cc_mec")    {h_shwr_length_numu_cc->Fill(shwr_length); }
			if(tpco_id == "nc")             {h_shwr_length_nc->Fill(shwr_length); }
			if(tpco_id == "nc_pi0")         {h_shwr_length_nc_pi0->Fill(shwr_length); }
			if(tpco_id == "nue_cc_mixed")   {h_shwr_length_nue_cc_mixed->Fill(shwr_length); }
			if(tpco_id == "numu_cc_mixed")  {h_shwr_length_numu_cc_mixed->Fill(shwr_length); }
			if(tpco_id == "cosmic")         {h_shwr_length_cosmic->Fill(shwr_length); }
			if(tpco_id == "other_mixed")    {h_shwr_length_other_mixed->Fill(shwr_length); }
			if(tpco_id == "unmatched")      {h_shwr_length_unmatched->Fill(shwr_length); }
		}
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingShowerLength(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                              std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                              double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                              double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shwr_length = leading_shower.pfpLength();

		if(tpco_id == "nue_cc_qe")      {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_out_fv")  {h_shwr_length_nue_cc_out_fv->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_res")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_dis")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_coh")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_mec")     {h_shwr_length_nue_cc->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_qe")     {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_res")    {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_dis")    {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_coh")    {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_mec")    {h_shwr_length_numu_cc->Fill(leading_shwr_length); }
		if(tpco_id == "nc")             {h_shwr_length_nc->Fill(leading_shwr_length); }
		if(tpco_id == "nc_pi0")         {h_shwr_length_nc_pi0->Fill(leading_shwr_length); }
		if(tpco_id == "nue_cc_mixed")   {h_shwr_length_nue_cc_mixed->Fill(leading_shwr_length); }
		if(tpco_id == "numu_cc_mixed")  {h_shwr_length_numu_cc_mixed->Fill(leading_shwr_length); }
		if(tpco_id == "cosmic")         {h_shwr_length_cosmic->Fill(leading_shwr_length); }
		if(tpco_id == "other_mixed")    {h_shwr_length_other_mixed->Fill(leading_shwr_length); }
		if(tpco_id == "unmatched")      {h_shwr_length_unmatched->Fill(leading_shwr_length); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::LeadingShowerTrackLengths(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                                    double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                                    double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const double leading_shwr_length = leading_shower.pfpLength();
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
		const double longest_trk_leading_shwr_ratio = longest_track / leading_shwr_length;

		if(tpco_id == "nue_cc_qe")      {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_out_fv")  {h_shwr_trk_length_nue_cc_out_fv->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_res")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_dis")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_coh")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_mec")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_qe")     {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_res")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_dis")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_coh")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_mec")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nc")             {h_shwr_trk_length_nc->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nc_pi0")         {h_shwr_trk_length_nc_pi0->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "nue_cc_mixed")   {h_shwr_trk_length_nue_cc_mixed->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "numu_cc_mixed")  {h_shwr_trk_length_numu_cc_mixed->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "cosmic")         {h_shwr_trk_length_cosmic->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "other_mixed")    {h_shwr_trk_length_other_mixed->Fill(longest_trk_leading_shwr_ratio); }
		if(tpco_id == "unmatched")      {h_shwr_trk_length_unmatched->Fill(longest_trk_leading_shwr_ratio); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::LongestShowerTrackLengths(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                                    double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                                    double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
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
		if(tpco_id == "nue_cc_out_fv")  {h_shwr_trk_length_nue_cc_out_fv->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_res")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_dis")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_coh")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_mec")     {h_shwr_trk_length_nue_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_qe")     {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_res")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_dis")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_coh")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_mec")    {h_shwr_trk_length_numu_cc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nc")             {h_shwr_trk_length_nc->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nc_pi0")         {h_shwr_trk_length_nc_pi0->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "nue_cc_mixed")   {h_shwr_trk_length_nue_cc_mixed->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "numu_cc_mixed")  {h_shwr_trk_length_numu_cc_mixed->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "cosmic")         {h_shwr_trk_length_cosmic->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "other_mixed")    {h_shwr_trk_length_other_mixed->Fill(longest_trk_longest_shwr_ratio); }
		if(tpco_id == "unmatched")      {h_shwr_trk_length_unmatched->Fill(longest_trk_longest_shwr_ratio); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PlaneHitsComparisonShower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                    std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                                    double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                                    double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
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
		if(tpco_id == "nue_cc_out_fv")  {h_collection_total_hits_shower_nue_cc_out_fv->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_res")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_dis")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_coh")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mec")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_qe")     {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_res")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_dis")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_coh")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mec")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc")             {h_collection_total_hits_shower_nc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc_pi0")         {h_collection_total_hits_shower_nc_pi0->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mixed")   {h_collection_total_hits_shower_nue_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mixed")  {h_collection_total_hits_shower_numu_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "cosmic")         {h_collection_total_hits_shower_cosmic->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "other_mixed")    {h_collection_total_hits_shower_other_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "unmatched")      {h_collection_total_hits_shower_unmatched->Fill(n_pfp_hits_w, n_pfp_hits); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PlaneHitsComparisonLeadingShower(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                           std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                                           double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                                           double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		//const int n_pfp_tracks = tpc_obj.NPfpTracks();
		//if(n_pfp_tracks == 0) {continue; }
		const int n_pfp_hits_w = leading_shower.NumPFPHitsW();
		const int n_pfp_hits = leading_shower.NumPFPHits();

		if(tpco_id == "nue_cc_qe")      {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_out_fv")  {h_collection_total_hits_shower_nue_cc_out_fv->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_res")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_dis")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_coh")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mec")     {h_collection_total_hits_shower_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_qe")     {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_res")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_dis")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_coh")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mec")    {h_collection_total_hits_shower_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc")             {h_collection_total_hits_shower_nc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc_pi0")         {h_collection_total_hits_shower_nc_pi0->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mixed")   {h_collection_total_hits_shower_nue_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mixed")  {h_collection_total_hits_shower_numu_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "cosmic")         {h_collection_total_hits_shower_cosmic->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "other_mixed")    {h_collection_total_hits_shower_other_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "unmatched")      {h_collection_total_hits_shower_unmatched->Fill(n_pfp_hits_w, n_pfp_hits); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::PlaneHitsComparisonTrack(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                                   std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                                   double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                                   double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
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
		if(tpco_id == "nue_cc_out_fv")  {h_collection_total_hits_track_nue_cc_out_fv->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_res")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_dis")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_coh")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mec")     {h_collection_total_hits_track_nue_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_qe")     {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_res")    {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_dis")    {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_coh")    {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mec")    {h_collection_total_hits_track_numu_cc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc")             {h_collection_total_hits_track_nc->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nc_pi0")         {h_collection_total_hits_track_nc_pi0->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "nue_cc_mixed")   {h_collection_total_hits_track_nue_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "numu_cc_mixed")  {h_collection_total_hits_track_numu_cc_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "cosmic")         {h_collection_total_hits_track_cosmic->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "other_mixed")    {h_collection_total_hits_track_other_mixed->Fill(n_pfp_hits_w, n_pfp_hits); }
		if(tpco_id == "unmatched")      {h_collection_total_hits_track_unmatched->Fill(n_pfp_hits_w, n_pfp_hits); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
void selection_functions::HitsPlots1D(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                      std::vector<std::pair<int, std::string> > * passed_tpco, bool _verbose, bool has_pi0,
                                      double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                      double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
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

		if(n_pfp_tracks != 0) {
			for(auto const n_pfp_hits_w_track : n_pfp_hits_w_track_v)
			{
				if(tpco_id == "nue_cc_qe")      {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_out_fv")  {h_collection_hits_track_nue_cc_out_fv->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_res")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_dis")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_coh")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_mec")     {h_collection_hits_track_nue_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_qe")     {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_res")    {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_dis")    {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_coh")    {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_mec")    {h_collection_hits_track_numu_cc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nc")             {h_collection_hits_track_nc->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nc_pi0")         {h_collection_hits_track_nc_pi0->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "nue_cc_mixed")   {h_collection_hits_track_nue_cc_mixed->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "numu_cc_mixed")  {h_collection_hits_track_numu_cc_mixed->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "cosmic")         {h_collection_hits_track_cosmic->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "other_mixed")    {h_collection_hits_track_other_mixed->Fill(n_pfp_hits_w_track); }
				if(tpco_id == "unmatched")      {h_collection_hits_track_unmatched->Fill(n_pfp_hits_w_track); }
			}
		}

		for(auto const n_pfp_hits_w_shower : n_pfp_hits_w_shower_v)
		{
			if(tpco_id == "nue_cc_qe")      {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_out_fv")  {h_collection_hits_shower_nue_cc_out_fv->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_res")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_dis")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_coh")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_mec")     {h_collection_hits_shower_nue_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_qe")     {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_res")    {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_dis")    {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_coh")    {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_mec")    {h_collection_hits_shower_numu_cc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nc")             {h_collection_hits_shower_nc->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nc_pi0")         {h_collection_hits_shower_nc_pi0->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "nue_cc_mixed")   {h_collection_hits_shower_nue_cc_mixed->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "numu_cc_mixed")  {h_collection_hits_shower_numu_cc_mixed->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "cosmic")         {h_collection_hits_shower_cosmic->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "other_mixed")    {h_collection_hits_shower_other_mixed->Fill(n_pfp_hits_w_shower); }
			if(tpco_id == "unmatched")      {h_collection_hits_shower_unmatched->Fill(n_pfp_hits_w_shower); }
		}

		if(tpco_id == "nue_cc_qe")      {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_out_fv")  {h_collection_hits_leading_shower_nue_cc_out_fv->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_res")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_dis")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_coh")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_mec")     {h_collection_hits_leading_shower_nue_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_qe")     {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_res")    {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_dis")    {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_coh")    {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_mec")    {h_collection_hits_leading_shower_numu_cc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nc")             {h_collection_hits_leading_shower_nc->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nc_pi0")         {h_collection_hits_leading_shower_nc_pi0->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "nue_cc_mixed")   {h_collection_hits_leading_shower_nue_cc_mixed->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "numu_cc_mixed")  {h_collection_hits_leading_shower_numu_cc_mixed->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "cosmic")         {h_collection_hits_leading_shower_cosmic->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "other_mixed")    {h_collection_hits_leading_shower_other_mixed->Fill(n_pfp_hits_w_leading_shower); }
		if(tpco_id == "unmatched")      {h_collection_hits_leading_shower_unmatched->Fill(n_pfp_hits_w_leading_shower); }

		if(tpco_id == "nue_cc_qe")      {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_out_fv")  {h_total_hits_leading_shower_nue_cc_out_fv->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_res")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_dis")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_coh")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_mec")     {h_total_hits_leading_shower_nue_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_qe")     {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_res")    {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_dis")    {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_coh")    {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_mec")    {h_total_hits_leading_shower_numu_cc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nc")             {h_total_hits_leading_shower_nc->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nc_pi0")         {h_total_hits_leading_shower_nc_pi0->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "nue_cc_mixed")   {h_total_hits_leading_shower_nue_cc_mixed->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "numu_cc_mixed")  {h_total_hits_leading_shower_numu_cc_mixed->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "cosmic")         {h_total_hits_leading_shower_cosmic->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "other_mixed")    {h_total_hits_leading_shower_other_mixed->Fill(n_pfp_hits_leading_shower); }
		if(tpco_id == "unmatched")      {h_total_hits_leading_shower_unmatched->Fill(n_pfp_hits_leading_shower); }
	}//end loop tpco
}//end function
//***************************************************************************
//***************************************************************************
int selection_functions::MapFailureCutToString(const std::string failure_cut)
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
void selection_functions::EnergyHits(std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v,
                                     std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0, bool _verbose,
                                     double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                     double vtxX, double vtxY, double vtxZ, double mc_nu_energy, double mc_ele_energy,
                                     TH2D * h_ele_eng_total_hits, TH2D * h_ele_eng_colleciton_hits, TH2D * h_nu_eng_total_hits, TH2D * h_nu_eng_collection_hits)
{
	int n_tpc_obj = tpc_object_container_v->size();
	for(int i = 0; i < n_tpc_obj; i++)
	{
		if(passed_tpco->at(i).first == 0) {continue; }
		auto const tpc_obj = tpc_object_container_v->at(i);
		const int n_pfp = tpc_obj.NumPFParticles();
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;
		const int leading_index = tpco_class.second;
		auto const leading_shower = tpc_obj.GetParticle(leading_index);
		const int n_pfp_hits_w  = leading_shower.NumPFPHitsW();
		const int n_pfp_hits    = leading_shower.NumPFPHits();
		if(tpco_id == "nue_cc_qe" || tpco_id == "nue_cc_res" || tpco_id == "nue_cc_coh" || tpco_id == "nue_cc_dis" || tpco_id == "nue_cc_mec")
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
                                        std::vector<std::pair<int, std::string> > * passed_tpco, bool has_pi0, bool _verbose,
                                        double _x1, double _x2, double _y1, double _y2, double _z1, double _z2,
                                        double vtxX, double vtxY, double vtxZ,
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
		std::pair<std::string, int> tpco_class = TPCO_Classifier(tpc_obj, has_pi0, _x1, _x2, _y1, _y2, _z1, _z2, vtxX, vtxY, vtxZ);
		std::string tpco_id = tpco_class.first;

		if(tpco_id == "nue_cc_qe")
		{
			h_failure_reason_nue_cc->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_out_fv")
		{
			h_failure_reason_nue_cc_out_fv->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_res")
		{
			h_failure_reason_nue_cc->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_dis")
		{
			h_failure_reason_nue_cc->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_coh")
		{
			h_failure_reason_nue_cc->Fill(failure_reason);
		}
		if(tpco_id == "nue_cc_mec")
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
			h_failure_reason_numu_cc_mixed->Fill(failure_reason);
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
//***************************************************************************
//***************************************************************************

//end functions
