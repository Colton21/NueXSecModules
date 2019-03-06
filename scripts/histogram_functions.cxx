#include "histogram_functions.h"

//void histogram_functions::Plot1DHistogram (TH1 * histogram, const char * x_axis_name, const char * print_name)
void histogram_functions::Plot1DHistogram (TH1 * histogram, const char * x_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();

	histogram->GetXaxis()->SetTitle(x_axis_name);
	histogram->Draw();

	c1->Print(print_name);
}

void histogram_functions::Plot1DHistogramGausFit (TH1 * histogram, const char * x_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();

	histogram->GetXaxis()->SetTitle(x_axis_name);
	histogram->Draw();

	histogram->Fit("gaus");

	c1->Print(print_name);
}

void histogram_functions::PlotTEfficiency (TH1 *h_num, TH1 *h_den, const char * title, const char * print_name)
{
	PlotTEfficiency(h_num, h_den, false, title, print_name);
}

void histogram_functions::PlotTEfficiency (TH1 *h_num, TH1 *h_den, const bool rebin, const char * title, const char * print_name)
{
	TCanvas * efficiency_c1 = new TCanvas();
	efficiency_c1->cd();
	TH1* h_num_clone = (TH1*)h_num->Clone("h_num_clone");
	TH1* h_den_clone = (TH1*)h_den->Clone("h_den_clone");

	if(rebin)
	{
		double new_bins [7] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 4.0};
		TH1* h_num_clone_rebin = (TH1*)h_num_clone->Rebin(6, "h_num_clone_rebin", new_bins);
		TH1* h_den_clone_rebin = (TH1*)h_den_clone->Rebin(6, "h_num_clone_rebin", new_bins);

		TEfficiency * teff = new TEfficiency(*h_num_clone_rebin, *h_den_clone_rebin);
		teff->SetTitle(title);
		teff->SetLineColor(kGreen+3);
		teff->SetMarkerColor(kGreen+3);
		teff->SetMarkerStyle(20);
		teff->SetMarkerSize(0.5);
		teff->Draw("AP");
		//std::cout << "----------" << std::endl;
		// for(int i = 1; i < h_num_clone_rebin->GetNbinsX()+1; i++)
		// {
		//      std::cout << "Bin: " << i << " Value: " << h_num_clone_rebin->GetBinContent(i) / h_den_clone_rebin->GetBinContent(i) << std::endl;
		// }
		// std::cout << "----------" << std::endl;

	}
	if(!rebin)
	{
		// std::cout << h_num_clone->GetNbinsX()<< ", " << h_den_clone->GetNbinsX() << std::endl;
		// // num
		// std::cout << "[";
		// for (int i(0); i< h_num_clone->GetSize(); i++) {std::cout << h_num_clone->GetBinContent(i) << ", "; }
		// std::cout << "]" << std::endl;
		// std::cout << "[";
		// for (int i(0); i< h_den_clone->GetSize(); i++) {std::cout << h_den_clone->GetBinContent(i) << ", "; }
		// std::cout << "]" << std::endl;
		TEfficiency * teff = new TEfficiency(*h_num_clone, *h_den_clone);
		teff->SetTitle(title);
		teff->SetLineColor(kGreen+3);
		teff->SetMarkerColor(kGreen+3);
		teff->SetMarkerStyle(20);
		teff->SetMarkerSize(0.5);
		teff->Draw("AP");
	}
	efficiency_c1->Print(print_name);
	std::cout << "Print Name: " << print_name << std::endl;
}
void histogram_functions::PlotTEfficiencyOverlay(TH1 * h_num,
                                                 TH1 * h_intime,
                                                 TH1 * h_pe,
                                                 TH1 * h_reco_nue,
                                                 TH1 * h_in_fv,
                                                 TH1 * h_vtx_flash,
                                                 TH1 * h_shwr_vtx,
                                                 TH1 * h_trk_vtx,
                                                 TH1 * h_hit,
                                                 TH1 * h_yhit,
                                                 TH1 * h_open_angle,
                                                 TH1 * h_dedx,
                                                 TH1 * h_2shwr,
                                                 TH1 * h_hit_len,
                                                 TH1 * h_trk_shwr,
                                                 TH1 * h_contain,
                                                 TH1 * h_den, const bool rebin, const char * title, const char * print_name)
{
	TCanvas * efficiency_c1 = new TCanvas();
	efficiency_c1->cd();
	TH1 * h_num_clone                    = (TH1*)h_num->Clone("h_num_clone");
	TH1 * h_nue_eng_eff_intime_clone     = (TH1*)h_intime->Clone("h_intime_clone");
	TH1 * h_nue_eng_eff_pe_clone         = (TH1*)h_pe->Clone("h_pe_clone");
	TH1 * h_nue_eng_eff_reco_nue_clone   = (TH1*)h_reco_nue->Clone("h_reco_nue_clone");
	TH1 * h_nue_eng_eff_in_fv_clone      = (TH1*)h_in_fv->Clone("h_in_fv_clone");
	TH1 * h_nue_eng_eff_vtx_flash_clone  = (TH1*)h_vtx_flash->Clone("h_vtx_flash_clone");
	TH1 * h_nue_eng_eff_shwr_vtx_clone   = (TH1*)h_shwr_vtx->Clone("h_shwr_vtx_clone");
	TH1 * h_nue_eng_eff_trk_vtx_clone    = (TH1*)h_trk_vtx->Clone("h_trk_vtx_clone");
	TH1 * h_nue_eng_eff_hit_clone        = (TH1*)h_hit->Clone("h_hit_clone");
	TH1 * h_nue_eng_eff_yhit_clone       = (TH1*)h_yhit->Clone("h_yhit_clone");
	TH1 * h_nue_eng_eff_open_angle_clone = (TH1*)h_open_angle->Clone("h_open_angle_clone");
	TH1 * h_nue_eng_eff_dedx_clone       = (TH1*)h_dedx->Clone("h_dedx_clone");
	TH1 * h_nue_eng_eff_2shwr_clone      = (TH1*)h_2shwr->Clone("h_2shwr_clone");
	TH1 * h_nue_eng_eff_hit_len_clone    = (TH1*)h_hit_len->Clone("h_hit_len_clone");
	TH1 * h_nue_eng_eff_trk_shwr_clone   = (TH1*)h_trk_shwr->Clone("h_trk_shwr_clone");
	TH1 * h_nue_eng_eff_contain_clone    = (TH1*)h_contain->Clone("h_contain_clone");
	TH1 * h_den_clone                    = (TH1*)h_den->Clone("h_den_clone");

	// std::cout << "Integral: " << h_num_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_reco_nue_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_in_fv_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_vtx_flash_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_shwr_vtx_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_trk_vtx_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_hit_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_yhit_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_open_angle_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_dedx_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_2shwr_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_hit_len_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_trk_shwr_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_nue_eng_eff_contain_clone->Integral() << std::endl;
	// std::cout << "Integral: " << h_den_clone->Integral() << std::endl;



	if(rebin)
	{
		double new_bins [7] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 4.0};
		TH1* h_num_clone_rebin = (TH1*)h_num_clone->Rebin(6, "h_num_clone_rebin", new_bins);
		TH1* h_den_clone_rebin = (TH1*)h_den_clone->Rebin(6, "h_num_clone_rebin", new_bins);
		TH1 * h_nue_eng_eff_reco_nue_clone_rebin   = (TH1*)h_nue_eng_eff_reco_nue_clone->Rebin(6,   "h_nue_eng_eff_reco_nue_clone_rebin",   new_bins);
		TH1 * h_nue_eng_eff_intime_clone_rebin     = (TH1*)h_nue_eng_eff_intime_clone->Rebin(6,     "h_nue_eng_eff_intime_clone_rebin",     new_bins);
		TH1 * h_nue_eng_eff_pe_clone_rebin         = (TH1*)h_nue_eng_eff_pe_clone->Rebin(6,         "h_nue_eng_eff_pe_clone_rebin",         new_bins);
		TH1 * h_nue_eng_eff_in_fv_clone_rebin      = (TH1*)h_nue_eng_eff_in_fv_clone->Rebin(6,      "h_nue_eng_eff_in_fv_clone_rebin",      new_bins);
		TH1 * h_nue_eng_eff_vtx_flash_clone_rebin  = (TH1*)h_nue_eng_eff_vtx_flash_clone->Rebin(6,  "h_nue_eng_eff_vtx_flash_clone_rebin",  new_bins);
		TH1 * h_nue_eng_eff_shwr_vtx_clone_rebin   = (TH1*)h_nue_eng_eff_shwr_vtx_clone->Rebin(6,   "h_nue_eng_eff_shwr_vtx_clone_rebin",   new_bins);
		TH1 * h_nue_eng_eff_trk_vtx_clone_rebin    = (TH1*)h_nue_eng_eff_trk_vtx_clone->Rebin(6,    "h_nue_eng_eff_trk_vtx_clone_rebin",    new_bins);
		TH1 * h_nue_eng_eff_hit_clone_rebin        = (TH1*)h_nue_eng_eff_hit_clone->Rebin(6,        "h_nue_eng_eff_hit_clone_rebin",        new_bins);
		TH1 * h_nue_eng_eff_yhit_clone_rebin       = (TH1*)h_nue_eng_eff_yhit_clone->Rebin(6,       "h_nue_eng_eff_yhit_clone_rebin",       new_bins);
		TH1 * h_nue_eng_eff_open_angle_clone_rebin = (TH1*)h_nue_eng_eff_open_angle_clone->Rebin(6, "h_nue_eng_eff_open_angle_clone_rebin", new_bins);
		TH1 * h_nue_eng_eff_dedx_clone_rebin       = (TH1*)h_nue_eng_eff_dedx_clone->Rebin(6,       "h_nue_eng_eff_dedx_clone_rebin",       new_bins);
		TH1 * h_nue_eng_eff_2shwr_clone_rebin      = (TH1*)h_nue_eng_eff_2shwr_clone->Rebin(6,      "h_nue_eng_eff_2shwr_clone_rebin",      new_bins);
		TH1 * h_nue_eng_eff_hit_len_clone_rebin    = (TH1*)h_nue_eng_eff_hit_len_clone->Rebin(6,    "h_nue_eng_eff_hit_len_clone_rebin",    new_bins);
		TH1 * h_nue_eng_eff_trk_shwr_clone_rebin   = (TH1*)h_nue_eng_eff_trk_shwr_clone->Rebin(6,   "h_nue_eng_eff_trk_shwr_clone_rebin",   new_bins);
		TH1 * h_nue_eng_eff_contain_clone_rebin    = (TH1*)h_nue_eng_eff_contain_clone->Rebin(6,    "h_nue_eng_eff_contain_clone_rebin",    new_bins);

		TEfficiency * teff            = new TEfficiency(*h_num_clone_rebin,                     *h_den_clone_rebin);
		TEfficiency * teff_reco_nue   = new TEfficiency(*h_nue_eng_eff_reco_nue_clone_rebin,    *h_den_clone_rebin);
		TEfficiency * teff_intime     = new TEfficiency(*h_nue_eng_eff_intime_clone_rebin,      *h_den_clone_rebin);
		TEfficiency * teff_pe         = new TEfficiency(*h_nue_eng_eff_pe_clone_rebin,          *h_den_clone_rebin);
		TEfficiency * teff_in_fv      = new TEfficiency(*h_nue_eng_eff_in_fv_clone_rebin,       *h_den_clone_rebin);
		TEfficiency * teff_vtx_flash  = new TEfficiency(*h_nue_eng_eff_vtx_flash_clone_rebin,   *h_den_clone_rebin);
		TEfficiency * teff_shwr_vtx   = new TEfficiency(*h_nue_eng_eff_shwr_vtx_clone_rebin,    *h_den_clone_rebin);
		TEfficiency * teff_trk_vtx    = new TEfficiency(*h_nue_eng_eff_trk_vtx_clone_rebin,     *h_den_clone_rebin);
		TEfficiency * teff_hit        = new TEfficiency(*h_nue_eng_eff_hit_clone_rebin,         *h_den_clone_rebin);
		TEfficiency * teff_yhit       = new TEfficiency(*h_nue_eng_eff_yhit_clone_rebin,        *h_den_clone_rebin);
		TEfficiency * teff_open_angle = new TEfficiency(*h_nue_eng_eff_open_angle_clone_rebin,  *h_den_clone_rebin);
		TEfficiency * teff_dedx       = new TEfficiency(*h_nue_eng_eff_dedx_clone_rebin,        *h_den_clone_rebin);
		TEfficiency * teff_2shwr      = new TEfficiency(*h_nue_eng_eff_2shwr_clone_rebin,       *h_den_clone_rebin);
		TEfficiency * teff_hit_len    = new TEfficiency(*h_nue_eng_eff_hit_len_clone_rebin,     *h_den_clone_rebin);
		TEfficiency * teff_trk_shwr   = new TEfficiency(*h_nue_eng_eff_trk_shwr_clone_rebin,    *h_den_clone_rebin);
		TEfficiency * teff_contain    = new TEfficiency(*h_nue_eng_eff_contain_clone_rebin,     *h_den_clone_rebin);

		teff->SetTitle(title);
		//teff->SetLineColor(kGreen+3);
		//teff->SetMarkerColor(kGreen+3);
		teff->SetMarkerStyle(20);
		teff->SetMarkerSize(0.5);
		teff_in_fv->SetMarkerStyle(20);
		teff_vtx_flash->SetMarkerStyle(20);
		teff_shwr_vtx->SetMarkerStyle(20);
		teff_trk_vtx->SetMarkerStyle(20);
		teff_hit->SetMarkerStyle(20);
		teff_yhit->SetMarkerStyle(20);
		teff_open_angle->SetMarkerStyle(20);
		teff_dedx->SetMarkerStyle(20);
		teff_2shwr->SetMarkerStyle(20);
		teff_hit_len->SetMarkerStyle(20);
		teff_trk_shwr->SetMarkerStyle(20);
		teff_contain->SetMarkerStyle(20);


		teff_in_fv->SetMarkerSize(0.5);
		teff_vtx_flash->SetMarkerSize(0.5);
		teff_shwr_vtx->SetMarkerSize(0.5);
		teff_trk_vtx->SetMarkerSize(0.5);
		teff_hit->SetMarkerSize(0.5);
		teff_yhit->SetMarkerSize(0.5);
		teff_open_angle->SetMarkerSize(0.5);
		teff_dedx->SetMarkerSize(0.5);
		teff_2shwr->SetMarkerSize(0.5);
		teff_hit_len->SetMarkerSize(0.5);
		teff_trk_shwr->SetMarkerSize(0.5);
		teff_contain->SetMarkerSize(0.5);

		teff->SetMarkerColor(30);
		teff_in_fv->SetMarkerColor(38);
		teff_vtx_flash->SetMarkerColor(43);
		teff_shwr_vtx->SetMarkerColor(9);
		teff_trk_vtx->SetMarkerColor(34);
		teff_hit->SetMarkerColor(8);
		teff_yhit->SetMarkerColor(49);
		teff_open_angle->SetMarkerColor(7);
		teff_dedx->SetMarkerColor(6);
		teff_2shwr->SetMarkerColor(41);
		teff_hit_len->SetMarkerColor(40);
		teff_trk_shwr->SetMarkerColor(46);
		teff_contain->SetMarkerColor(1);

		teff->SetLineColor(30);
		teff_in_fv->SetLineColor(38);
		teff_vtx_flash->SetLineColor(43);
		teff_shwr_vtx->SetLineColor(9);
		teff_trk_vtx->SetLineColor(34);
		teff_hit->SetLineColor(8);
		teff_yhit->SetLineColor(49);
		teff_open_angle->SetLineColor(7);
		teff_dedx->SetLineColor(6);
		teff_2shwr->SetLineColor(41);
		teff_hit_len->SetLineColor(40);
		teff_trk_shwr->SetLineColor(46);
		teff_contain->SetLineColor(1);

		//teff->Draw("AP");
		teff_in_fv->Draw("AP");
		teff_vtx_flash->Draw("PSAME");
		//teff_shwr_vtx->Draw("PSAME");
		teff_trk_vtx->Draw("PSAME");
		//teff_hit->Draw("PSAME");
		teff_yhit->Draw("PSAME");
		//teff_open_angle->Draw("PSAME");
		teff_dedx->Draw("PSAME");
		// teff_2shwr->Draw("PSAME");
		// teff_hit_len->Draw("PSAME");
		// teff_trk_shwr->Draw("PSAME");
		teff_contain->Draw("PSAME");

		// gPad->Update();
		// auto graph = teff->GetPaintedGraph();
		// graph->SetMinimum(0.0);
		// graph->SetMaximum(0.80);
		// gPad->Update();
		efficiency_c1->Update();
		teff_in_fv->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 0.80);
		teff_in_fv->GetPaintedGraph()->GetXaxis()->SetLabelSize(12);
		// teff->GetPaintedGraph()->GetXaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
		// teff->GetPaintedGraph()->GetYaxis()->SetLabelSize(11);
		// teff->GetPaintedGraph()->GetYaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
		// teff->GetPaintedGraph()->GetXaxis()->SetTitleOffset(3.6);
		teff_in_fv->GetPaintedGraph()->GetXaxis()->SetTitleSize(16);
		teff_in_fv->GetPaintedGraph()->GetXaxis()->SetTitleFont(46);

		teff_in_fv->GetPaintedGraph()->GetYaxis()->SetTitleSize(16);
		teff_in_fv->GetPaintedGraph()->GetYaxis()->SetTitleFont(46);
		efficiency_c1->Update();

		TLegend * leg1 = new TLegend(0.75, 0.98, 0.98, 0.60);
		leg1->AddEntry(teff_in_fv,           "Reco #nu_{e} in Fid. Vol.", "l");
		//leg1->AddEntry(teff_reco_nue,  "Reco #nu_{e} InFV",  "l");
		leg1->AddEntry(teff_vtx_flash, "Vertex-Flash",      "l");
		leg1->AddEntry(teff_trk_vtx,   "Shower/Track to #nu", "l");
		//leg1->AddEntry(teff_hit,       "Hit",           "l");
		leg1->AddEntry(teff_yhit,      "Hit Threshold",          "l");
		leg1->AddEntry(teff_dedx,       "dE/dx",        "l");
		leg1->AddEntry(teff_contain,    "Track Contain",   "l");
		leg1->Draw();

	}
	efficiency_c1->Print(print_name);
}

void histogram_functions::PlotTEfficiency (TH2 *h_num, TH2 *h_den, const char * title, const char * print_name)
{
	TCanvas * efficiency_c1 = new TCanvas();
	efficiency_c1->cd();
	TH2* h_num_clone = (TH2*)h_num->Clone("h_num_clone");
	TH2* h_den_clone = (TH2*)h_den->Clone("h_den_clone");

	//std::cout << h_num_clone->GetNbinsX()<< ", " << h_den_clone->GetNbinsX() << std::endl;
	// num
	// std::cout << "[";
	// for (int i(0); i< h_num_clone->GetSize(); i++) {std::cout << h_num_clone->GetBinContent(i) << ", "; }
	// std::cout << "]" << std::endl;
	// std::cout << "[";
	// for (int i(0); i< h_den_clone->GetSize(); i++) {std::cout << h_den_clone->GetBinContent(i) << ", "; }
	// std::cout << "]" << std::endl;
	TEfficiency * teff = new TEfficiency(*h_num_clone, *h_den_clone);
	teff->SetTitle(title);
	teff->SetLineColor(kGreen+3);
	teff->SetMarkerColor(kGreen+3);
	teff->SetMarkerStyle(20);
	teff->SetMarkerSize(0.5);
	teff->Draw("AP");

	efficiency_c1->Print(print_name);
	std::cout << "Print Name: " << print_name << std::endl;
}

void histogram_functions::Plot2DHistogram (TH2 * histogram, const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	histogram->GetXaxis()->SetTitle(x_axis_name);
	histogram->GetYaxis()->SetTitle(y_axis_name);
	histogram->SetTitle(title);
	histogram->SetStats(kFALSE);
	histogram->Draw("colz");
	c1->Print(print_name);
}

void histogram_functions::Plot2DHistogramNormZ (TH2 * histogram_1, TH2 * histogram_2, const char * title_1, const char * title_2,
                                                const char * x_axis_name, const char * y_axis_name, const char * print_name_1, const char * print_name_2)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	histogram_1->GetXaxis()->SetTitle(x_axis_name);
	histogram_1->GetYaxis()->SetTitle(y_axis_name);
	histogram_1->SetTitle(title_1);
	histogram_1->SetStats(kFALSE);
	histogram_1->Draw("colz");

	TCanvas * c2 = new TCanvas();
	c2->cd();
	histogram_2->GetXaxis()->SetTitle(x_axis_name);
	histogram_2->GetYaxis()->SetTitle(y_axis_name);
	histogram_2->SetTitle(title_2);
	histogram_2->SetStats(kFALSE);
	histogram_2->Draw("colz");

	histogram_1->GetZaxis()->SetRangeUser(histogram_2->GetMinimum(),//GetBinLowEdge(1),
	                                      histogram_2->GetMaximum());//GetBinLowEdge(histogram_2->GetNbinsZ()+1));
	std::cout << "\t" << histogram_2->GetZaxis()->GetBinLowEdge(1) << ", " << histogram_2->GetZaxis()->GetBinLowEdge(histogram_2->GetNbinsZ()+1) << std::endl;
	c1->cd();
	histogram_1->Draw("colz");
	c2->cd();
	histogram_2->Draw("colz");

	c1->Print(print_name_1);
	c2->Print(print_name_2);
}

void histogram_functions::TimingHistograms(TH1 * histogram_1, TH1 * histogram_2, TH1 * histogram_3, TH1 * histogram_4,
                                           const double data_scale_factor, const double intime_scale_factor, const double dirt_scale_factor,
                                           const char * x_axis_name, const char * print_name, const char * print_name2)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();

	TH1 * h_1_clone = (TH1*)histogram_1->Clone("h_1_clone");
	TH1 * h_2_clone = (TH1*)histogram_2->Clone("h_2_clone");
	TH1 * h_3_clone = (TH1*)histogram_3->Clone("h_3_clone");
	TH1 * h_4_clone = (TH1*)histogram_4->Clone("h_4_clone");

	h_1_clone->Sumw2();
	h_2_clone->Sumw2();
	h_3_clone->Sumw2();
	h_4_clone->Sumw2();

	h_1_clone->SetFillColor(30);
	h_2_clone->SetFillColor(9);
	h_4_clone->SetFillColor(46);

	h_1_clone->Scale(data_scale_factor);
	h_2_clone->Scale(intime_scale_factor);
	TH1 * h_3_clone_clone = (TH1*)h_3_clone->Clone("h_3_clone_clone");
	h_3_clone->Add(h_2_clone, -1);
	h_4_clone->Scale(dirt_scale_factor);

	h_3_clone->GetXaxis()->SetTitle(x_axis_name);
	h_3_clone->GetYaxis()->SetRangeUser(-3000, 3000);

	h_3_clone->Draw("hist");
	h_1_clone->SetLineColor(46);
	h_1_clone->Draw("hist same");

	c1->Print(print_name);

	TCanvas * c2 = new TCanvas();
	c2->cd();

	THStack * stack = new THStack();

	stack->Add(h_2_clone);
	stack->Add(h_4_clone);
	stack->Add(h_1_clone);

	stack->Draw("hist");

	stack->GetXaxis()->SetTitle(x_axis_name);
	stack->Draw("hist");

	//we want this before the off-beam is subtracted
	h_3_clone_clone->Draw("same");

	c2->Print(print_name2);

}

//used to overlay both the intime and data flash histograms
void histogram_functions::TimingHistogramsOverlay(std::vector<std::pair<double, int> > * data_flash_time, TH1 * histogram_1, TH1 * histogram_2,
                                                  const double intime_scale_factor, const char * x_axis_name,
                                                  const char * print_name1, const char * print_name2)
{
	TCanvas * c1a = new TCanvas();
	c1a->cd();

	TH1 * h_1_clone = (TH1*)histogram_1->Clone("h_1_clone");
	TH1 * h_2_clone = (TH1*)histogram_2->Clone("h_2_clone");

	h_1_clone->GetXaxis()->SetTitle(x_axis_name);
	h_1_clone->SetLineColor(46);
	h_1_clone->Sumw2();
	h_2_clone->Sumw2();
	h_1_clone->Scale(intime_scale_factor);

	h_1_clone->SetStats(kFALSE);
	h_2_clone->SetStats(kFALSE);
	h_1_clone->SetTitle("");
	h_1_clone->Draw("e");
	h_2_clone->Draw("e same");

	TLegend * leg1 = new TLegend(0.7, 0.85, 0.98, 0.98);
	leg1->AddEntry(h_1_clone,    "NuMI EXT Run 1",      "l");
	leg1->AddEntry(h_2_clone,    "NuMI On-Beam Run 1",  "l");
	//leg1->SetTextFont(132);
	leg1->Draw();

	c1a->Print(print_name1);

	TCanvas * c1b = new TCanvas();
	c1b->cd();

	TH1 * h_divide_clone = (TH1*)h_2_clone->Clone("h_divide_clone");
	h_divide_clone->Divide(h_1_clone);
	h_divide_clone->Draw("e");
	c1b->Print("../scripts/plots/numi_timing_on_off_divide.pdf");

	TCanvas * c2 = new TCanvas();
	c2->cd();
	TH1 * h_3_clone = (TH1*)h_2_clone->Clone("h_3_clone");
	h_3_clone->Add(h_1_clone, -1);
	h_3_clone->SetStats(kFALSE);
	//h_3_clone->SetBins(80, 0, 20);
	h_3_clone->Sumw2();
	h_3_clone->SetTitle("");
	h_3_clone->Draw("e");

	TLegend * leg2 = new TLegend(0.7, 0.85, 0.98, 0.98);
	leg2->AddEntry(h_3_clone,    "NuMI On-Beam - Off-Beam Run 1",  "l");
	//leg2->SetTextFont(132);
	leg2->Draw();

	c2->Print(print_name2);

	TH1D * h_flash_time_data_first_half   = new TH1D ("h_flash_time_data_first_half",    "h_flash_time_data_first_half",   80, 0, 20);
	TH1D * h_flash_time_data_second_half  = new TH1D ("h_flash_time_data_second_half",   "h_flash_time_data_second_half",  80, 0, 20);

	h_flash_time_data_first_half->Sumw2();
	h_flash_time_data_second_half->Sumw2();

	for(auto const flash_timing : * data_flash_time)
	{
		const double run_number = flash_timing.second;
		if(run_number <= 6450) {h_flash_time_data_first_half->Fill(flash_timing.first); }
		if(run_number > 6450) {h_flash_time_data_second_half->Fill(flash_timing.first); }
	}

	Plot1DHistogram(h_flash_time_data_first_half,   "Flash Time [#mus]", "../scripts/plots/flash_time_data_first_half.pdf");
	Plot1DHistogram(h_flash_time_data_second_half,  "Flash Time [#mus]", "../scripts/plots/flash_time_data_second_half.pdf");

	TCanvas * c3a = new TCanvas();
	c3a->cd();
	const double integral_1 = h_flash_time_data_first_half->Integral();
	const double integral_2 = h_flash_time_data_second_half->Integral();
	//h_flash_time_data_first_half->Scale(1./integral);
	h_flash_time_data_second_half->Scale(integral_1 / integral_2);
	//h_flash_time_data_second_half->Integral();


	h_flash_time_data_first_half->Sumw2();
	h_flash_time_data_second_half->Sumw2();
	h_flash_time_data_first_half->SetStats(kFALSE);
	h_flash_time_data_second_half->SetStats(kFALSE);

	h_flash_time_data_first_half->GetXaxis()->SetTitle("Flash Time [#mus]");
	h_flash_time_data_first_half->SetTitle("");
	h_flash_time_data_first_half->Draw("p e");
	h_flash_time_data_second_half->SetLineColor(46);
	h_flash_time_data_second_half->Draw("p e same");

	TLegend * leg3a = new TLegend(0.7, 0.85, 0.98, 0.98);
	leg3a->AddEntry(h_flash_time_data_first_half,  "NuMI On-Beam Run 1 (<= 6450)", "l");
	leg3a->AddEntry(h_flash_time_data_second_half, "NuMI On-Beam Run 1 (> 6450)", "l");
	//leg3a->SetTextFont(132);
	leg3a->Draw();

	c3a->Print("../scripts/plots/flash_time_data_first_second_half.pdf");


	TCanvas * c3b = new TCanvas();
	c3b->cd();
	TH1 * h_flash_time_data_divide = (TH1*)h_flash_time_data_first_half->Clone("h_flash_time_data_divide");
	h_flash_time_data_divide->Divide(h_flash_time_data_second_half);
	h_flash_time_data_divide->GetYaxis()->SetRangeUser(0.8, 1.3);
	h_flash_time_data_divide->GetYaxis()->SetTitleSize(17);
	h_flash_time_data_divide->GetYaxis()->SetTitleOffset(4);
	h_flash_time_data_divide->GetXaxis()->SetTitleSize(17);
	h_flash_time_data_divide->GetXaxis()->SetTitleOffset(10);
	h_flash_time_data_divide->SetTitleSize(17);
	h_flash_time_data_divide->SetTitleOffset(10);
	h_flash_time_data_divide->GetXaxis()->SetTitle("Flash Time [#mus]");
	h_flash_time_data_divide->SetTitle("NuMI Run 1 On-Beam First Half / Second Half");
	h_flash_time_data_divide->GetYaxis()->SetTitle("Ratio Flashes 1st Half / 2nd Half");
	h_flash_time_data_divide->SetStats(kFALSE);
	h_flash_time_data_divide->Sumw2();
	h_flash_time_data_divide->SetTitle("");
	h_flash_time_data_divide->Draw("e");

	TLine * line = new TLine(0, 1, 20, 1);
	line->SetLineColor(46);
	line->Draw("same");

	TLegend * leg3b = new TLegend(0.7, 0.85, 0.98, 0.98);
	leg3b->AddEntry(h_flash_time_data_divide, "NuMI Run 1 1st-Half / 2nd-Half", "l");
	//leg3b->SetTextFont(132);
	leg3b->Draw();

	c3b->Print("../scripts/plots/flash_time_data_divide.pdf");

}
void histogram_functions::PlotFlashInfo(TH1 * h_flash_mc, TH1 * h_flash_intime, TH1 * h_flash_data, TH1 * h_flash_dirt,
                                        const double intime_scale_factor, const double data_scale_factor, const double dirt_scale_factor,
                                        const char * x_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	THStack * stack = new THStack();

	h_flash_mc->SetStats(kFALSE);
	h_flash_intime->SetStats(kFALSE);
	h_flash_data->SetStats(kFALSE);
	h_flash_dirt->SetStats(kFALSE);

	h_flash_mc->SetFillColor(49);
	h_flash_intime->SetFillColor(41);
	h_flash_intime->SetFillStyle(3345);
	h_flash_dirt->SetFillColor(30);

	TH1 * h_mc_clone       = (TH1*)h_flash_mc->Clone("h_mc_clone");
	TH1 * h_intime_clone   = (TH1*)h_flash_intime->Clone("h_intime_clone");
	TH1 * h_dirt_clone     = (TH1*)h_flash_dirt->Clone("h_dirt_clone");

	const double y_maximum = std::max(h_flash_data->GetMaximum(), stack->GetMaximum());
	stack->SetMaximum(y_maximum * 1.2);

	h_mc_clone->Sumw2();
	h_intime_clone->Sumw2();
	h_dirt_clone->Sumw2();

	h_mc_clone->Scale(data_scale_factor);
	h_intime_clone->Scale(intime_scale_factor);
	h_dirt_clone->Scale(dirt_scale_factor);

	h_flash_data->SetMarkerStyle(20);
	h_flash_data->SetMarkerSize(0.5);
	h_flash_data->Sumw2();

	stack->Add(h_mc_clone);
	stack->Add(h_intime_clone);
	stack->Add(h_dirt_clone);

	stack->Draw("hist");
	stack->GetXaxis()->SetTitle(x_axis_name);
	h_flash_data->Draw("same PE");

	TH1 * h_error_hist = (TH1*)h_mc_clone->Clone("h_error_hist");
	h_error_hist->Add(h_intime_clone, 1);
	h_error_hist->Add(h_dirt_clone, 1);

	h_error_hist->SetFillColorAlpha(12, 0.15);
	h_error_hist->Draw("e2 hist same");

	TLegend * leg_stack = new TLegend(0.85,0.85,0.95,0.95);
	leg_stack->AddEntry(h_flash_mc,      "MC",   "f");
	leg_stack->AddEntry(h_flash_intime,  "EXT",  "f");
	leg_stack->AddEntry(h_flash_dirt,    "Dirt", "f");
	leg_stack->Draw();
	c1->Print(print_name);

}


void histogram_functions::LegoStackData(TH2 * h_nue_cc, TH2 * h_nue_cc_mixed, TH2 * h_nue_cc_out_fv, TH2 * h_numu_cc, TH2 * h_cosmic, TH2 * h_nc,
                                        TH2 * h_nc_pi0, TH2 * h_other_mixed, TH2 * h_unmatched, TH2 * h_intime, const double intime_scale_factor,
                                        TH2 * h_data, const double data_scale_factor,
                                        const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                        const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	THStack * stack = new THStack();
	h_nue_cc->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_numu_cc->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_intime->SetStats(kFALSE);
	h_data->SetStats(kFALSE);
	h_nue_cc->SetFillColor(30);
	h_nue_cc_mixed->SetFillColor(38);
	h_numu_cc->SetFillColor(28);
	h_nc_pi0->SetFillColor(36);
	h_cosmic->SetFillStyle(1001);
	h_cosmic->SetFillColor(kBlack);
	h_nc->SetFillColor(46);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	h_intime->SetFillColor(41);
	h_intime->SetFillStyle(3345);
	//h_intime->SetFillStyle(3354);
	TH2 * h_intime_clone        = (TH2*)h_intime->Clone("h_intime_clone");
	TH2 * h_nue_cc_clone        = (TH2*)h_nue_cc->Clone("h_nue_cc_clone");
	TH2 * h_nue_cc_mixed_clone  = (TH2*)h_nue_cc_mixed->Clone("h_nue_cc_mixed_clone");
	TH2 * h_nue_cc_out_fv_clone = (TH2*)h_nue_cc_out_fv->Clone("h_nue_cc_out_fv_clone");
	TH2 * h_cosmic_clone        = (TH2*)h_cosmic->Clone("h_cosmic_clone");
	TH2 * h_numu_cc_clone       = (TH2*)h_numu_cc->Clone("h_numu_cc_clone");
	TH2 * h_nc_clone            = (TH2*)h_nc->Clone("h_nc_clone");
	TH2 * h_nc_pi0_clone        = (TH2*)h_nc_pi0->Clone("h_nc_pi0_clone");
	TH2 * h_other_mixed_clone   = (TH2*)h_other_mixed->Clone("h_other_mixed_clone");
	TH2 * h_unmatched_clone     = (TH2*)h_unmatched->Clone("h_unmatched_clone");

	h_nue_cc_clone->Sumw2();
	h_nue_cc_mixed_clone->Sumw2();
	h_nue_cc_out_fv_clone->Sumw2();
	h_cosmic_clone->Sumw2();
	h_numu_cc_clone->Sumw2();
	h_nc_clone->Sumw2();
	h_nc_pi0_clone->Sumw2();
	h_other_mixed_clone->Sumw2();
	h_unmatched_clone->Sumw2();
	h_intime_clone->Sumw2();

	h_nue_cc_clone->Scale(data_scale_factor);
	h_nue_cc_mixed_clone->Scale(data_scale_factor);
	h_nue_cc_out_fv_clone->Scale(data_scale_factor);
	h_cosmic_clone->Scale(data_scale_factor);
	h_numu_cc_clone->Scale(data_scale_factor);
	h_nc_clone->Scale(data_scale_factor);
	h_nc_pi0_clone->Scale(data_scale_factor);
	h_other_mixed_clone->Scale(data_scale_factor);
	h_unmatched_clone->Scale(data_scale_factor);
	h_intime_clone->Scale(intime_scale_factor);

	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(0.5);
	h_data->Sumw2();

	stack->Add(h_nue_cc_clone);
	stack->Add(h_nue_cc_mixed_clone);
	stack->Add(h_cosmic_clone);
	stack->Add(h_numu_cc_clone);
	stack->Add(h_nc_clone);
	stack->Add(h_nc_pi0_clone);
	stack->Add(h_other_mixed_clone);
	stack->Add(h_unmatched_clone);
	stack->Add(h_intime_clone);

	stack->Draw("lego4 0");
	stack->GetXaxis()->SetTitle(x_axis_name);
	//h_data->Draw("lego1 0 same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "Nue CC",        "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "Nue CC Mixed",  "f");
	leg_stack->AddEntry(h_nue_cc_out_fv,   "Nue CC Out FV", "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic",        "f");
	leg_stack->AddEntry(h_numu_cc,         "Numu CC",       "f");
	leg_stack->AddEntry(h_nc,              "NC",            "f");
	leg_stack->AddEntry(h_nc_pi0,          "NC Pi0",        "f");
	leg_stack->AddEntry(h_other_mixed,     "Other Mixed",   "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",     "f");
	leg_stack->AddEntry(h_intime,          "InTime",        "f");
	leg_stack->Draw();
	c1->Print(print_name);
}
void histogram_functions::Plot2DHistogram (TH2 * histogram, const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name,
                                           const char * draw_option)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	histogram->GetXaxis()->SetTitle(x_axis_name);
	histogram->GetYaxis()->SetTitle(y_axis_name);
	histogram->SetTitle(title);
	histogram->SetStats(kFALSE);
	histogram->Draw(draw_option);
	c1->Print(print_name);
}
void histogram_functions::Plot2DHistogramLogZ (TH2 * histogram, const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name,
                                               const char * draw_option)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	histogram->GetXaxis()->SetTitle(x_axis_name);
	histogram->GetYaxis()->SetTitle(y_axis_name);
	c1->SetLogz();
	histogram->SetTitle(title);
	histogram->SetStats(kFALSE);
	histogram->Draw(draw_option);
	c1->Print(print_name);
}
void histogram_functions::Plot2DHistogram (TH2 * histogram, const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name,
                                           const char * draw_option, const int x_divisions, const int y_divisions)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	histogram->GetXaxis()->SetTitle(x_axis_name);
	histogram->GetXaxis()->SetNdivisions(x_divisions);
	histogram->GetYaxis()->SetTitle(y_axis_name);
	histogram->GetYaxis()->SetNdivisions(y_divisions);
	histogram->SetTitle(title);
	histogram->SetStats(kFALSE);
	histogram->Draw(draw_option);
	c1->Print(print_name);
}

void histogram_functions::Plot2DHistogramOffSet (TH2 * histogram, const double label_offset, const double title_offset, const char * title,
                                                 const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	histogram->Draw("colz");
	histogram->GetYaxis()->SetLabelOffset(label_offset);
	histogram->GetYaxis()->SetTitleOffset(title_offset);
	histogram->GetYaxis()->SetTitle(x_axis_name);
	histogram->GetXaxis()->SetTitle(y_axis_name);
	histogram->SetTitle(title);
	histogram->SetStats(kFALSE);
	c1->Print(print_name);
}
void histogram_functions::PlotSimpleStack (TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed,
                                           TH1 * h_cosmic, TH1 * h_nc,
                                           TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, const char * title,
                                           const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	//default legend position
	//0.75,0.75,0.95,0.95
	PlotSimpleStack(h_nue_cc, h_nue_cc_mixed, h_nue_cc_out_fv, h_numu_cc, h_numu_cc_mixed, h_cosmic, h_nc,
	                h_nc_pi0, h_other_mixed, h_unmatched,
	                0.75, 0.95, 0.75, 0.95,
	                title, x_axis_name, y_axis_name, print_name);
}
void histogram_functions::PlotSimpleStackInTime (TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                                 TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                                 TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime,
                                                 const double intime_scale_factor, const double data_scale_factor,
                                                 const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	//default legend position
	//0.75,0.75,0.95,0.95
	PlotSimpleStackInTime(h_nue_cc, h_nue_cc_mixed, h_nue_cc_out_fv, h_numu_cc, h_numu_cc_mixed, h_cosmic, h_nc,
	                      h_nc_pi0, h_other_mixed, h_unmatched, h_intime, intime_scale_factor, data_scale_factor,
	                      0.80, 0.95, 0.75, 0.95,
	                      title, x_axis_name, y_axis_name, print_name);
}
void histogram_functions::PlotSimpleStackData (TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                               TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                               TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                               TH1 * h_data, const double data_scale_factor, TH1 * h_dirt, const double dirt_scale_factor,
                                               const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	//default legend position
	//0.75,0.75,0.95,0.95

	//I recommend xlow = 0.5, xhigh = 0.85, yhigh = 0.85, and ylow = 0.85-(0.05*nlabels)
	PlotSimpleStackData(h_nue_cc, h_nue_cc_mixed, h_nue_cc_out_fv, h_numu_cc, h_numu_cc_mixed, h_cosmic, h_nc,
	                    h_nc_pi0, h_other_mixed, h_unmatched, h_intime, intime_scale_factor,
	                    h_data, data_scale_factor, h_dirt, dirt_scale_factor,
	                    0.75, 0.98, 0.98, 0.50, false, false, 1.2,
	                    title, x_axis_name, y_axis_name, print_name);
}
void histogram_functions::PlotSimpleStackData (TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                               TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                               TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                               TH1 * h_data, const double data_scale_factor, TH1 * h_dirt, const double dirt_scale_factor,
                                               const double x_min, const double x_max, const double y_min, const double y_max,
                                               const bool logy, const bool area_norm,
                                               const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	//default legend position
	//0.75,0.75,0.95,0.95

	//I recommend xlow = 0.5, xhigh = 0.85, yhigh = 0.85, and ylow = 0.85-(0.05*nlabels)
	PlotSimpleStackData(h_nue_cc, h_nue_cc_mixed, h_nue_cc_out_fv, h_numu_cc, h_numu_cc_mixed, h_cosmic, h_nc,
	                    h_nc_pi0, h_other_mixed, h_unmatched, h_intime, intime_scale_factor,
	                    h_data, data_scale_factor, h_dirt, dirt_scale_factor,
	                    x_min, x_max, y_min, y_max, logy, area_norm, 1.2,
	                    title, x_axis_name, y_axis_name, print_name);
}

void histogram_functions::PlotSimpleStackData (TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                               TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                               TH1 * h_data, const double data_scale_factor, TH1 * h_dirt, const double dirt_scale_factor, const double y_scale_factor,
                                               const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	PlotSimpleStackData(h_nue_cc, h_nue_cc_mixed, h_nue_cc_out_fv, h_numu_cc, h_numu_cc_mixed, h_cosmic, h_nc,
	                    h_nc_pi0, h_other_mixed, h_unmatched, h_intime, intime_scale_factor,
	                    h_data, data_scale_factor, h_dirt, dirt_scale_factor,
	                    0.745, 0.98, 0.98, 0.50, false, false, y_scale_factor,
	                    title, x_axis_name, y_axis_name, print_name);
}

void histogram_functions::PlotSimpleStackParticle(TH1 * h_electron, TH1 * h_proton, TH1 * h_photon, TH1 * h_pion,
                                                  TH1 * h_kaon, TH1 * h_muon, TH1 * h_neutron, TH1 * h_unmatched,
                                                  TH1 * h_unmatched_ext, const double intime_scale_factor, const double data_scale_factor,
                                                  const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                                  const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	THStack * stack = new THStack();
	h_electron->SetStats(kFALSE);
	h_proton->SetStats(kFALSE);
	h_photon->SetStats(kFALSE);
	h_pion->SetStats(kFALSE);
	h_kaon->SetStats(kFALSE);
	h_muon->SetStats(kFALSE);
	h_neutron->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_unmatched_ext->SetStats(kFALSE);
	h_electron->SetFillColor(30);
	h_proton->SetFillColor(38);
	h_photon->SetFillColor(28);
	h_pion->SetFillColor(36);
	h_kaon->SetFillColor(1);
	h_muon->SetFillColor(46);
	h_neutron->SetFillColor(kTeal);
	h_unmatched->SetFillColor(12);
	h_unmatched_ext->SetFillColor(41);
	h_unmatched_ext->SetFillStyle(3345);

	TH1 * h_electron_clone      = (TH1*)h_electron->Clone("h_electron_clone");
	TH1 * h_proton_clone        = (TH1*)h_proton->Clone("h_proton_clone");
	TH1 * h_photon_clone        = (TH1*)h_photon->Clone("h_photon_clone");
	TH1 * h_pion_clone          = (TH1*)h_pion->Clone("h_pion_clone");
	TH1 * h_kaon_clone          = (TH1*)h_kaon->Clone("h_kaon_clone");
	TH1 * h_muon_clone          = (TH1*)h_muon->Clone("h_muon_clone");
	TH1 * h_neutron_clone       = (TH1*)h_neutron->Clone("h_neutron_clone");
	TH1 * h_unmatched_clone     = (TH1*)h_unmatched->Clone("h_unmatched_clone");
	TH1 * h_unmatched_ext_clone = (TH1*)h_unmatched_ext->Clone("h_unmatched_ext_clone");

	h_electron_clone->Sumw2();
	h_proton_clone->Sumw2();
	h_photon_clone->Sumw2();
	h_pion_clone->Sumw2();
	h_kaon_clone->Sumw2();
	h_muon_clone->Sumw2();
	h_neutron_clone->Sumw2();
	h_unmatched_clone->Sumw2();
	h_unmatched_ext_clone->Sumw2();

	h_electron_clone->Scale(data_scale_factor);
	h_proton_clone->Scale(data_scale_factor);
	h_photon_clone->Scale(data_scale_factor);
	h_pion_clone->Scale(data_scale_factor);
	h_kaon_clone->Scale(data_scale_factor);
	h_muon_clone->Scale(data_scale_factor);
	h_neutron_clone->Scale(data_scale_factor);
	h_unmatched_clone->Scale(data_scale_factor);
	h_unmatched_ext_clone->Scale(intime_scale_factor);

	stack->Add(h_electron_clone);
	stack->Add(h_proton_clone);
	stack->Add(h_photon_clone);
	stack->Add(h_pion_clone);
	stack->Add(h_kaon_clone);
	stack->Add(h_muon_clone);
	stack->Add(h_neutron_clone);
	stack->Add(h_unmatched_clone);
	stack->Add(h_unmatched_ext_clone);
	stack->Draw("hist");
	stack->GetXaxis()->SetTitle(x_axis_name);

	TH1 * h_error_hist = (TH1*)h_electron_clone->Clone("h_error_hist");
	h_error_hist->Add(h_proton_clone, 1);
	h_error_hist->Add(h_photon_clone, 1);
	h_error_hist->Add(h_pion_clone, 1);
	h_error_hist->Add(h_kaon_clone, 1);
	h_error_hist->Add(h_muon_clone, 1);
	h_error_hist->Add(h_neutron_clone, 1);
	h_error_hist->Add(h_unmatched_clone, 1);
	h_error_hist->Add(h_unmatched_ext_clone, 1);

	h_error_hist->SetFillColorAlpha(12, 0.15);
	h_error_hist->Draw("e2 hist same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_electron,        "Electron",   "f");
	leg_stack->AddEntry(h_proton,          "Proton",     "f");
	leg_stack->AddEntry(h_photon,          "Photon",     "f");
	leg_stack->AddEntry(h_pion,            "Pion",       "f");
	leg_stack->AddEntry(h_kaon,            "Kaon",       "f");
	leg_stack->AddEntry(h_muon,            "Muon",       "f");
	leg_stack->AddEntry(h_neutron,         "Neutron",    "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",  "f");
	leg_stack->AddEntry(h_unmatched_ext,   "EXT",        "f");
	leg_stack->Draw();
	c1->Print(print_name);

	delete h_electron_clone;
	delete h_proton_clone;
	delete h_photon_clone;
	delete h_pion_clone;
	delete h_kaon_clone;
	delete h_muon_clone;
	delete h_neutron_clone;
	delete h_unmatched_clone;
	delete h_unmatched_ext_clone;
	delete stack;
}

void histogram_functions::PlotSimpleStack(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                          TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                          TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched,
                                          const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                          const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	THStack * stack = new THStack();
	h_nue_cc->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_numu_cc->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_numu_cc_mixed->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_nue_cc_out_fv->SetStats(kFALSE);
	h_nue_cc->SetFillColor(30);
	h_nue_cc_mixed->SetFillColor(38);
	h_numu_cc->SetFillColor(28);
	h_nc_pi0->SetFillColor(36);
	h_cosmic->SetFillColor(1);
	h_nc->SetFillColor(46);
	h_nue_cc_out_fv->SetFillColor(kTeal);
	//h_numu_cc_mixed->SetFillColor(28);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "#nu_{e} CC",       "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "#nu_{e} CC Mixed", "f");
	leg_stack->AddEntry(h_nue_cc_out_fv,   "#nu_{e} CC OutFV", "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic",           "f");
	leg_stack->AddEntry(h_numu_cc,         "#nu_{#mu} CC",     "f");
	//leg_stack->AddEntry(h_numu_cc_mixed,   "Numu CC Mixed",  "f");
	leg_stack->AddEntry(h_nc,              "NC",               "f");
	leg_stack->AddEntry(h_nc_pi0,          "NC #pi^{0}",       "f");
	leg_stack->AddEntry(h_other_mixed,     "NC Mixed",         "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",        "f");
	leg_stack->Draw();
	c1->Print(print_name);

	delete c1;
	delete stack;
}
void histogram_functions::PlotSimpleStackInTime(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                                TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                                TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime,
                                                const double intime_scale_factor, const double data_scale_factor,
                                                const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                                const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	THStack * stack = new THStack();
	h_nue_cc->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_numu_cc->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_numu_cc_mixed->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_intime->SetStats(kFALSE);
	h_nue_cc_out_fv->SetStats(kFALSE);
	h_nue_cc->SetFillColor(30);
	h_nue_cc_mixed->SetFillColor(38);
	h_numu_cc->SetFillColor(28);
	h_nc_pi0->SetFillColor(36);
	h_cosmic->SetFillColor(1);
	h_nc->SetFillColor(46);
	h_nue_cc_out_fv->SetFillColor(kTeal);
	//h_numu_cc_mixed->SetFillColor(28);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	h_intime->SetFillColor(41);
	h_intime->SetFillStyle(3345);
	//h_intime->SetFillStyle(3354);
	TH1 * h_intime_clone = (TH1*)h_intime->Clone("h_intime_clone");
	TH1 * h_nue_cc_clone = (TH1*)h_nue_cc->Clone("h_nue_cc_clone");
	TH1 * h_nue_cc_mixed_clone = (TH1*)h_nue_cc_mixed->Clone("h_nue_cc_mixed_clone");
	TH1 * h_nue_cc_out_fv_clone = (TH1*)h_nue_cc_out_fv->Clone("h_nue_cc_out_fv_clone");
	TH1 * h_cosmic_clone = (TH1*)h_cosmic->Clone("h_cosmic_clone");
	TH1 * h_numu_cc_clone = (TH1*)h_numu_cc->Clone("h_numu_cc_clone");
	TH1 * h_nc_clone = (TH1*)h_nc->Clone("h_nc_clone");
	TH1 * h_nc_pi0_clone = (TH1*)h_nc_pi0->Clone("h_nc_pi0_clone");
	TH1 * h_other_mixed_clone = (TH1*)h_other_mixed->Clone("h_other_mixed_clone");
	TH1 * h_unmatched_clone = (TH1*)h_unmatched->Clone("h_unmatched_clone");

	h_nue_cc_clone->Sumw2();
	h_nue_cc_mixed_clone->Sumw2();
	h_nue_cc_out_fv_clone->Sumw2();
	h_cosmic_clone->Sumw2();
	h_numu_cc_clone->Sumw2();
	h_nc_clone->Sumw2();
	h_nc_pi0_clone->Sumw2();
	h_other_mixed_clone->Sumw2();
	h_unmatched_clone->Sumw2();
	h_intime_clone->Sumw2();

	h_nue_cc_clone->Scale(data_scale_factor);
	h_nue_cc_mixed_clone->Scale(data_scale_factor);
	h_nue_cc_out_fv_clone->Scale(data_scale_factor);
	h_cosmic_clone->Scale(data_scale_factor);
	h_numu_cc_clone->Scale(data_scale_factor);
	h_nc_clone->Scale(data_scale_factor);
	h_nc_pi0_clone->Scale(data_scale_factor);
	h_other_mixed_clone->Scale(data_scale_factor);
	h_unmatched_clone->Scale(data_scale_factor);
	h_intime_clone->Scale(intime_scale_factor);

	stack->Add(h_nue_cc_clone);
	stack->Add(h_nue_cc_mixed_clone);
	stack->Add(h_nue_cc_out_fv_clone);
	stack->Add(h_cosmic_clone);
	stack->Add(h_numu_cc_clone);
	stack->Add(h_nc_clone);
	stack->Add(h_nc_pi0_clone);
	stack->Add(h_other_mixed_clone);
	stack->Add(h_unmatched_clone);
	stack->Add(h_intime_clone);
	stack->Draw("hist");
	stack->GetXaxis()->SetTitle(x_axis_name);

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "Nue CC",        "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "Nue CC Mixed",  "f");
	leg_stack->AddEntry(h_nue_cc_out_fv,   "Nue CC OutFV",  "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic",        "f");
	leg_stack->AddEntry(h_numu_cc,         "Numu CC",       "f");
	//leg_stack->AddEntry(h_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack->AddEntry(h_nc,              "NC",            "f");
	leg_stack->AddEntry(h_nc_pi0,          "NC Pi0",        "f");
	leg_stack->AddEntry(h_other_mixed,     "NC Mixed",      "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",     "f");
	leg_stack->AddEntry(h_intime,          "InTime",        "f");
	leg_stack->Draw();
	c1->Print(print_name);

	delete h_nue_cc_clone;
	delete h_nue_cc_mixed_clone;
	delete h_nue_cc_out_fv_clone;
	delete h_cosmic_clone;
	delete h_numu_cc_clone;
	delete h_nc_clone;
	delete h_nc_pi0_clone;
	delete h_other_mixed_clone;
	delete h_unmatched_clone;
	delete h_intime_clone;
	delete stack;
	delete c1;
}
void histogram_functions::PlotSimpleStackData(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                              TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                              TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                              TH1 * h_data, const double data_scale_factor, TH1 * h_dirt, const double dirt_scale_factor,
                                              const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                              const bool logy, const bool area_norm,
                                              const double y_scale_factor,
                                              const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{

	const bool p_value = false;

	TCanvas * c1 = new TCanvas(title, title, 500, 500);
	c1->cd();

	//TPad *pad2 = new TPad("pad2", "pad2", 0, 0.3, 1, 1.0);
	//TPad *pad2_2 = new TPad("pad2_2", "pad2_2", 0, 0.05, 1, 0.3);

	//TPad *topPad = new TPad("topPad", "", 0.005, 0.32, 0.995, 0.995);
	//TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.28);
	TPad * topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
	TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
	//was 0.01
	topPad->SetBottomMargin(0.05);
	//bottomPad->SetTopMargin(0.0);
	//bottomPad->SetBottomMargin(0.12);
	bottomPad->SetTopMargin(0.04);
	bottomPad->SetBottomMargin(0.25);
	bottomPad->SetGridy();
	topPad->Draw();
	bottomPad->Draw();
	topPad->cd();

	THStack * stack = new THStack();
	h_nue_cc->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_numu_cc->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_intime->SetStats(kFALSE);
	h_data->SetStats(kFALSE);
	h_nue_cc_out_fv->SetStats(kFALSE);
	h_nue_cc->SetFillColor(30);
	h_nue_cc_mixed->SetFillColor(38);
	h_numu_cc->SetFillColor(28);
	h_nc_pi0->SetFillColor(36);
	h_cosmic->SetFillColor(1);
	h_nc->SetFillColor(46);
	h_nue_cc_out_fv->SetFillColor(kTeal);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	h_intime->SetFillColor(41);
	h_intime->SetFillStyle(3345);

	h_dirt->SetFillColor(2);
	h_dirt->SetFillStyle(3354);

	//h_data->Scale(data_scale_factor);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(0.5);
	//h_data->Sumw2();

	TH1 * h_nue_cc_clone        = (TH1*)h_nue_cc->Clone("h_nue_cc_clone");
	TH1 * h_nue_cc_mixed_clone  = (TH1*)h_nue_cc_mixed->Clone("h_nue_cc_mixed_clone");
	TH1 * h_nue_cc_out_fv_clone = (TH1*)h_nue_cc_out_fv->Clone("h_nue_cc_out_fv_clone");
	TH1 * h_cosmic_clone        = (TH1*)h_cosmic->Clone("h_cosmic_clone");
	TH1 * h_numu_cc_clone       = (TH1*)h_numu_cc->Clone("h_numu_cc_clone");
	TH1 * h_nc_clone            = (TH1*)h_nc->Clone("h_nc_clone");
	TH1 * h_nc_pi0_clone        = (TH1*)h_nc_pi0->Clone("h_nc_pi0_clone");
	TH1 * h_other_mixed_clone   = (TH1*)h_other_mixed->Clone("h_other_mixed_clone");
	TH1 * h_unmatched_clone     = (TH1*)h_unmatched->Clone("h_unmatched_clone");
	TH1 * h_intime_clone        = (TH1*)h_intime->Clone("h_intime_clone");
	TH1 * h_dirt_clone          = (TH1*)h_dirt->Clone("h_dirt_clone");
	TH1 * h_data_clone          = (TH1*)h_data->Clone("h_data_clone");

	h_nue_cc_clone->Sumw2();
	h_nue_cc_mixed_clone->Sumw2();
	h_nue_cc_out_fv_clone->Sumw2();
	h_numu_cc_clone->Sumw2();
	h_nc_pi0_clone->Sumw2();
	h_cosmic_clone->Sumw2();
	h_nc_clone->Sumw2();
	h_other_mixed_clone->Sumw2();
	h_unmatched_clone->Sumw2();
	h_intime_clone->Sumw2();
	h_dirt_clone->Sumw2();
	h_data_clone->Sumw2();

	h_nue_cc_clone->Scale(data_scale_factor);
	h_nue_cc_mixed_clone->Scale(data_scale_factor);
	h_nue_cc_out_fv_clone->Scale(data_scale_factor);
	h_cosmic_clone->Scale(data_scale_factor);
	h_numu_cc_clone->Scale(data_scale_factor);
	h_nc_clone->Scale(data_scale_factor);
	h_nc_pi0_clone->Scale(data_scale_factor);
	h_other_mixed_clone->Scale(data_scale_factor);
	h_unmatched_clone->Scale(data_scale_factor);
	h_intime_clone->Scale(intime_scale_factor);
	h_dirt_clone->Scale(dirt_scale_factor);

	double integral_data = h_data_clone->Integral();

	if(area_norm)
	{
		double integral_mc_ext = h_nue_cc_clone->Integral() +
		                         h_nue_cc_mixed_clone->Integral() +
		                         h_nue_cc_out_fv_clone->Integral() +
		                         h_cosmic_clone->Integral() +
		                         h_numu_cc_clone->Integral() +
		                         h_nc_clone->Integral() +
		                         h_nc_pi0_clone->Integral() +
		                         h_other_mixed_clone->Integral() +
		                         h_unmatched_clone->Integral() +
		                         h_dirt_clone->Integral(); // +
		//h_intime_clone->Integral();

		TH1 * h_data_scaling_clone = (TH1*)h_data->Clone("h_data_scaling_clone");
		h_data_scaling_clone->Add(h_intime_clone, -1);
		const double integral_on_minus_off = h_data_scaling_clone->Integral();

		h_nue_cc_clone->Scale(integral_on_minus_off / integral_mc_ext);
		h_nue_cc_mixed_clone->Scale(integral_on_minus_off / integral_mc_ext);
		h_nue_cc_out_fv_clone->Scale(integral_on_minus_off / integral_mc_ext);
		h_cosmic_clone->Scale(integral_on_minus_off / integral_mc_ext);
		h_numu_cc_clone->Scale(integral_on_minus_off / integral_mc_ext);
		h_nc_clone->Scale(integral_on_minus_off / integral_mc_ext);
		h_nc_pi0_clone->Scale(integral_on_minus_off / integral_mc_ext);
		h_other_mixed_clone->Scale(integral_on_minus_off / integral_mc_ext);
		h_unmatched_clone->Scale(integral_on_minus_off / integral_mc_ext);
		h_dirt_clone->Scale(integral_on_minus_off / integral_mc_ext);

		// h_nue_cc_clone->Scale(integral_data / integral_mc_ext);
		// h_nue_cc_mixed_clone->Scale(integral_data / integral_mc_ext);
		// h_nue_cc_out_fv_clone->Scale(integral_data / integral_mc_ext);
		// h_cosmic_clone->Scale(integral_data / integral_mc_ext);
		// h_numu_cc_clone->Scale(integral_data / integral_mc_ext);
		// h_nc_clone->Scale(integral_data / integral_mc_ext);
		// h_nc_pi0_clone->Scale(integral_data / integral_mc_ext);
		// h_other_mixed_clone->Scale(integral_data / integral_mc_ext);
		// h_unmatched_clone->Scale(integral_data / integral_mc_ext);
		//h_intime_clone->Scale(integral_data / integral_mc_ext);

		h_nue_cc_clone->Scale(1. / integral_data);
		h_nue_cc_mixed_clone->Scale(1. / integral_data);
		h_nue_cc_out_fv_clone->Scale(1. / integral_data);
		h_cosmic_clone->Scale(1. / integral_data);
		h_numu_cc_clone->Scale(1. / integral_data);
		h_nc_clone->Scale(1. / integral_data);
		h_nc_pi0_clone->Scale(1. / integral_data);
		h_other_mixed_clone->Scale(1. / integral_data);
		h_unmatched_clone->Scale(1. / integral_data);
		h_intime_clone->Scale(1. / integral_data);
		h_dirt_clone->Scale(1. / integral_data);
		h_data_clone->Scale(1. / integral_data);

	}


	stack->Add(h_nue_cc_clone);
	stack->Add(h_nue_cc_mixed_clone);
	stack->Add(h_nue_cc_out_fv_clone);
	stack->Add(h_cosmic_clone);
	stack->Add(h_numu_cc_clone);
	stack->Add(h_nc_clone);
	stack->Add(h_nc_pi0_clone);
	stack->Add(h_other_mixed_clone);
	stack->Add(h_unmatched_clone);
	stack->Add(h_dirt_clone);
	stack->Add(h_intime_clone);

	const double y_maximum = std::max(h_data_clone->GetMaximum(), stack->GetMaximum());
	//stack->SetMaximum(y_maximum * y_scale_factor);

	if(logy == true)
	{
		TH1 * h_scale_axes = (TH1*)h_data_clone->Clone("h_scale_axes");
		//h_scale_axes->GetYaxis()->SetRangeUser(0.1, y_maximum * (y_scale_factor * 100));
		if(h_nue_cc_clone->GetMinimum() != 0.0) {h_scale_axes->SetMinimum(h_nue_cc_clone->GetMinimum() / 2.); }
		if(h_nue_cc_clone->GetMinimum() == 0.0) {h_scale_axes->SetMinimum(h_nue_cc_clone->GetMinimum() + 0.0001 / 2.); }
		h_scale_axes->SetMaximum(y_maximum * (y_scale_factor * 500));
		h_scale_axes->SetLineColor(0);
		h_scale_axes->SetFillColor(0);
		h_scale_axes->GetYaxis()->SetTitle("Entries [A.U.]");
		h_scale_axes->SetTitle(" ");
		h_scale_axes->GetXaxis()->SetTitle(" ");
		h_scale_axes->GetXaxis()->SetLabelSize(0);
		h_scale_axes->GetXaxis()->SetLabelFont(0); // Absolute font size in pixel (precision 3)
		h_scale_axes->Draw();
		stack->Draw("same hist");

		h_scale_axes->GetYaxis()->SetTitle("Entries [A.U.]");
	}

	if(!logy)
	{
		stack->SetMinimum(0);
		stack->SetMaximum(y_maximum * y_scale_factor);
		stack->Draw("hist");
	}
	//stack->GetXaxis()->SetTitle(x_axis_name);
	if(!area_norm) {stack->GetYaxis()->SetTitle("Entries"); }
	if(area_norm) {stack->GetYaxis()->SetTitle("Entries [A.U.]"); }

	stack->GetYaxis()->SetTitleFont(46);
	stack->GetYaxis()->SetTitleSize(17);
	stack->GetXaxis()->SetLabelOffset(10);
	h_data_clone->Draw("same PE");

	TH1 * h_error_hist = (TH1*)h_nue_cc_clone->Clone("h_error_hist");
	h_error_hist->Add(h_nue_cc_mixed_clone,  1);
	h_error_hist->Add(h_nue_cc_out_fv_clone, 1);
	h_error_hist->Add(h_numu_cc_clone,       1);
	h_error_hist->Add(h_nc_pi0_clone,        1);
	h_error_hist->Add(h_nc_clone,            1);
	h_error_hist->Add(h_other_mixed_clone,   1);
	h_error_hist->Add(h_cosmic_clone,        1);
	h_error_hist->Add(h_unmatched_clone,     1);
	h_error_hist->Add(h_dirt_clone,          1);
	h_error_hist->Add(h_intime_clone,        1);

	h_error_hist->SetFillColorAlpha(12, 0.15);
	h_error_hist->Draw("e2 hist same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "#nu_{e} CC",        "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "#nu_{e} CC Mixed",  "f");
	leg_stack->AddEntry(h_nue_cc_out_fv,   "#nu_{e} CC OutFV",  "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic",            "f");
	leg_stack->AddEntry(h_numu_cc,         "#nu_{#mu} CC",      "f");
	leg_stack->AddEntry(h_nc,              "NC",                "f");
	leg_stack->AddEntry(h_nc_pi0,          "NC #pi^{0}",        "f");
	leg_stack->AddEntry(h_other_mixed,     "NC Mixed",          "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",         "f");
	leg_stack->AddEntry(h_dirt,            "Dirt",              "f");
	leg_stack->AddEntry(h_intime,          "InTime",            "f");
	leg_stack->Draw();

	if(!logy) {stack->GetYaxis()->SetRangeUser(0, y_maximum * y_scale_factor); }
	if(logy == true)
	{
		topPad->SetLogy();
		//stack->GetYaxis()->SetRangeUser(h_nue_cc_clone->GetMinimum(), y_maximum * y_scale_factor);
	}

	TH1 * h_last = (TH1*) stack->GetStack()->Last();
	std::vector <double> chi2  = Chi2Calc(h_last, h_data_clone, area_norm, integral_data);
	//chi2 : chi2/ndf, mc+ext, data

	//x_min, y_min, x_max, y_max
	//reduced chi2
	TPaveText * pt = new TPaveText(.46,.80,.73,1.06, "NBNDC");
	std::ostringstream o_string;
	o_string.precision(3);
	o_string << std::fixed;
	o_string << float(chi2.at(0));
	std::string convert_string = o_string.str();
	std::string chi2_string = "#chi_{Stat}^{2}/DOF=" + convert_string;
	pt->AddText(chi2_string.c_str());
	pt->SetFillStyle(0);
	pt->SetBorderSize(0);
	//pt->Draw();

	//num events
	TPaveText * pt2 = new TPaveText(.13,.80,.46,1.06, "NBNDC");
	std::ostringstream o_string2a;
	std::ostringstream o_string2b;
	if(!area_norm)
	{
		o_string2a << int(chi2.at(2));
		o_string2b << int(chi2.at(1));
	}
	if(area_norm)
	{
		o_string2a << int(chi2.at(2) * integral_data);
		o_string2b << int(chi2.at(1) * integral_data);
	}
	std::string convert_string2a = o_string2a.str();
	std::string convert_string2b = o_string2b.str();
	std::string chi2_string2 = "Data: " + convert_string2a + "|MC+EXT:" + convert_string2b;
	pt2->AddText(chi2_string2.c_str());
	pt2->SetFillStyle(0);
	pt2->SetBorderSize(0);
	// this is removed for public distributions
	//pt2->Draw();

	//num bins
	TPaveText * pt3 = new TPaveText(.60,.80,.73,.973, "NBNDC");
	std::ostringstream o_string3;
	o_string3 << int(chi2.at(3));
	std::string convert_string3 = o_string3.str();
	std::string ndf_string = "DOF=" + convert_string3;
	pt3->AddText(ndf_string.c_str());
	pt3->SetFillStyle(0);
	pt3->SetBorderSize(0);
	//pt3->Draw();

	//p value
	//optional
	TPaveText * pt4 = new TPaveText(.45,.80,.60,.973, "NBNDC");
	std::ostringstream o_string4;
	o_string4.precision(4);
	o_string4 << std::fixed;
	o_string4 << chi2.at(4);
	std::string convert_string4 = o_string4.str();
	std::string p_string = "P=" + convert_string4;
	pt4->AddText(p_string.c_str());
	pt4->SetFillStyle(0);
	pt4->SetBorderSize(0);
	if(p_value) {pt4->Draw(); }

	bottomPad->cd();
	TH1 * ratioPlot = (TH1*)h_data_clone->Clone("ratioPlot");
	// ratioPlot->Add(h_nue_cc_clone,        -1);
	// ratioPlot->Add(h_nue_cc_mixed_clone,  -1);
	// ratioPlot->Add(h_nue_cc_out_fv_clone, -1);
	// ratioPlot->Add(h_cosmic_clone,        -1);
	// ratioPlot->Add(h_numu_cc_clone,       -1);
	// ratioPlot->Add(h_nc_clone,            -1);
	// ratioPlot->Add(h_nc_pi0_clone,        -1);
	// ratioPlot->Add(h_other_mixed_clone,   -1);
	// ratioPlot->Add(h_unmatched_clone,     -1);
	// ratioPlot->Add(h_intime_clone,        -1);
	//ratioPlot->Divide(h_data);
	TH1 * h_mc_ext_sum = (TH1*)h_nue_cc_clone->Clone("h_mc_ext_sum");
	//h_mc_ext_sum->Add(h_nue_cc_clone,        1);
	h_mc_ext_sum->Add(h_nue_cc_mixed_clone,  1);
	h_mc_ext_sum->Add(h_nue_cc_out_fv_clone, 1);
	h_mc_ext_sum->Add(h_cosmic_clone,        1);
	h_mc_ext_sum->Add(h_numu_cc_clone,       1);
	h_mc_ext_sum->Add(h_nc_clone,            1);
	h_mc_ext_sum->Add(h_nc_pi0_clone,        1);
	h_mc_ext_sum->Add(h_other_mixed_clone,   1);
	h_mc_ext_sum->Add(h_unmatched_clone,     1);
	h_mc_ext_sum->Add(h_dirt_clone,          1);
	h_mc_ext_sum->Add(h_intime_clone,        1);

	ratioPlot->GetXaxis()->SetLabelSize(12);
	ratioPlot->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	ratioPlot->GetYaxis()->SetLabelSize(11);
	ratioPlot->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	ratioPlot->GetXaxis()->SetTitleOffset(3.0);
	ratioPlot->GetXaxis()->SetTitleSize(16);
	ratioPlot->GetXaxis()->SetTitleFont(45);

	ratioPlot->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

	ratioPlot->Add(h_mc_ext_sum, -1);
	ratioPlot->Divide(h_mc_ext_sum);
	ratioPlot->GetYaxis()->SetRangeUser(-1,1);
	ratioPlot->GetXaxis()->SetTitle(x_axis_name);
	ratioPlot->GetYaxis()->SetTitle("(Data - MC) / MC ");
	ratioPlot->GetYaxis()->SetTitleSize(17);
	ratioPlot->GetYaxis()->SetTitleFont(46);
	ratioPlot->GetYaxis()->SetTitleOffset(2);
	ratioPlot->SetTitle(" ");
	ratioPlot->Draw();

	//now doing this stuff on the bottom pad

	//x_min, y_min, x_max, y_max
	//reduced chi2
	TPaveText * pt_bottom = new TPaveText(.12, .80, .30, .96, "NBNDC");
	std::ostringstream o_string_bottom;
	o_string_bottom.precision(3);
	o_string_bottom << std::fixed;
	o_string_bottom << float(chi2.at(0) * chi2.at(3));
	std::string convert_string_bottom = o_string_bottom.str();

	std::ostringstream o_string3_bottom;
	o_string3_bottom << int(chi2.at(3));
	std::string convert_string3_bottom = o_string3_bottom.str();

	std::string chi2_string_bottom = "#chi_{Stat}^{2}/DOF=(" + convert_string_bottom + "/" + convert_string3_bottom + ")";
	pt_bottom->AddText(chi2_string_bottom.c_str());
	pt_bottom->SetFillStyle(0);
	pt_bottom->SetBorderSize(0);
	pt_bottom->Draw();


	//num bins
	// TPaveText * pt3 = new TPaveText(.60,.80,.73,.973, "NBNDC");
	// std::ostringstream o_string3;
	// o_string3 << int(chi2.at(3));
	// std::string convert_string3 = o_string3.str();
	// std::string ndf_string = "DOF=" + convert_string3;
	// pt3->AddText(ndf_string.c_str());
	// pt3->SetFillStyle(0);
	// pt3->SetBorderSize(0);
	// pt3->Draw();

	std::cout << print_name << std::endl;
	c1->Print(print_name);

	delete h_nue_cc_clone;
	delete h_nue_cc_mixed_clone;
	delete h_nue_cc_out_fv_clone;
	delete h_cosmic_clone;
	delete h_numu_cc_clone;
	delete h_nc_clone;
	delete h_nc_pi0_clone;
	delete h_other_mixed_clone;
	delete h_unmatched_clone;
	delete h_intime_clone;
	delete h_dirt_clone;
	delete h_data_clone;

	delete stack;
	delete leg_stack;
	delete pt;
	delete pt2;
	delete pt3;
	delete pt4;
	delete pt_bottom;

	delete h_error_hist;
	delete ratioPlot;
	delete h_mc_ext_sum;
	delete topPad;
	delete bottomPad;
	delete c1;
}

void histogram_functions::PlotSimpleStackDataMomentumRebin(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                                           TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                                           TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime,
                                                           const double intime_scale_factor, TH1 * h_data, const double data_scale_factor,
                                                           TH1 * h_dirt, const double dirt_scale_factor,
                                                           const bool area_norm,
                                                           const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{

	const bool p_value = false;

	TCanvas * c1 = new TCanvas(title, title, 500, 500);
	c1->cd();

	//TPad *pad2 = new TPad("pad2", "pad2", 0, 0.3, 1, 1.0);
	//TPad *pad2_2 = new TPad("pad2_2", "pad2_2", 0, 0.05, 1, 0.3);

	//TPad *topPad = new TPad("topPad", "", 0.005, 0.32, 0.995, 0.995);
	//TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.28);
	TPad * topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
	TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
	//was 0.01
	topPad->SetBottomMargin(0.05);
	//bottomPad->SetTopMargin(0.0);
	//bottomPad->SetBottomMargin(0.12);
	bottomPad->SetTopMargin(0.04);
	bottomPad->SetBottomMargin(0.25);
	bottomPad->SetGridy();
	topPad->Draw();
	bottomPad->Draw();
	topPad->cd();


	THStack * stack = new THStack();
	h_nue_cc->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_numu_cc->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_intime->SetStats(kFALSE);
	h_data->SetStats(kFALSE);
	h_nue_cc_out_fv->SetStats(kFALSE);
	h_nue_cc->SetFillColor(30);
	h_nue_cc_mixed->SetFillColor(38);
	h_numu_cc->SetFillColor(28);
	h_nc_pi0->SetFillColor(36);
	h_cosmic->SetFillColor(1);
	h_nc->SetFillColor(46);
	h_nue_cc_out_fv->SetFillColor(kTeal);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	h_intime->SetFillColor(41);
	h_intime->SetFillStyle(3345);

	h_dirt->SetFillColor(2);
	h_dirt->SetFillStyle(3354);

	//h_data->Scale(data_scale_factor);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(0.5);
	h_data->Sumw2();

	// TH1 * h_nue_cc_clone        = (TH1*)h_nue_cc->Clone("h_nue_cc_clone");
	// TH1 * h_nue_cc_mixed_clone  = (TH1*)h_nue_cc_mixed->Clone("h_nue_cc_mixed_clone");
	// TH1 * h_nue_cc_out_fv_clone = (TH1*)h_nue_cc_out_fv->Clone("h_nue_cc_out_fv_clone");
	// TH1 * h_cosmic_clone        = (TH1*)h_cosmic->Clone("h_cosmic_clone");
	// TH1 * h_numu_cc_clone       = (TH1*)h_numu_cc->Clone("h_numu_cc_clone");
	// TH1 * h_nc_clone            = (TH1*)h_nc->Clone("h_nc_clone");
	// TH1 * h_nc_pi0_clone        = (TH1*)h_nc_pi0->Clone("h_nc_pi0_clone");
	// TH1 * h_other_mixed_clone   = (TH1*)h_other_mixed->Clone("h_other_mixed_clone");
	// TH1 * h_unmatched_clone     = (TH1*)h_unmatched->Clone("h_unmatched_clone");
	// TH1 * h_intime_clone        = (TH1*)h_intime->Clone("h_intime_clone");

	double new_bins [5] = {0.0, 0.20, 0.40, 0.8, 4.0};
	const double num_bins = 4; //new_bins size - 1

	TH1 * h_nue_cc_rebin        = (TH1*)h_nue_cc->Rebin(num_bins, "h_nue_cc_rebin", new_bins);
	TH1 * h_nue_cc_mixed_rebin  = (TH1*)h_nue_cc_mixed->Rebin(num_bins, "h_nue_cc_mixed_rebin", new_bins);
	TH1 * h_nue_cc_out_fv_rebin = (TH1*)h_nue_cc_out_fv->Rebin(num_bins, "h_nue_cc_out_fv_rebin", new_bins);
	TH1 * h_cosmic_rebin        = (TH1*)h_cosmic->Rebin(num_bins, "h_cosmic_rebin", new_bins);
	TH1 * h_numu_cc_rebin       = (TH1*)h_numu_cc->Rebin(num_bins, "h_numu_cc_rebin", new_bins);
	TH1 * h_nc_rebin            = (TH1*)h_nc->Rebin(num_bins, "h_nc_rebin", new_bins);
	TH1 * h_nc_pi0_rebin        = (TH1*)h_nc_pi0->Rebin(num_bins, "h_nc_pi0_rebin", new_bins);
	TH1 * h_other_mixed_rebin   = (TH1*)h_other_mixed->Rebin(num_bins, "h_other_mixed_rebin", new_bins);
	TH1 * h_unmatched_rebin     = (TH1*)h_unmatched->Rebin(num_bins, "h_unmatched_rebin", new_bins);
	TH1 * h_intime_rebin        = (TH1*)h_intime->Rebin(num_bins, "h_intime_rebin", new_bins);
	TH1 * h_dirt_rebin          = (TH1*)h_dirt->Rebin(num_bins, "h_dirt_rebin", new_bins);
	TH1 * h_data_rebin          = (TH1*)h_data->Rebin(num_bins, "h_data_rebin", new_bins);

	h_nue_cc_rebin->Sumw2();
	h_nue_cc_mixed_rebin->Sumw2();
	h_nue_cc_out_fv_rebin->Sumw2();
	h_numu_cc_rebin->Sumw2();
	h_nc_pi0_rebin->Sumw2();
	h_cosmic_rebin->Sumw2();
	h_nc_rebin->Sumw2();
	h_other_mixed_rebin->Sumw2();
	h_unmatched_rebin->Sumw2();
	h_intime_rebin->Sumw2();
	h_dirt_rebin->Sumw2();
	h_data_rebin->Sumw2();

	h_nue_cc_rebin->Scale(data_scale_factor);
	h_nue_cc_mixed_rebin->Scale(data_scale_factor);
	h_nue_cc_out_fv_rebin->Scale(data_scale_factor);
	h_cosmic_rebin->Scale(data_scale_factor);
	h_numu_cc_rebin->Scale(data_scale_factor);
	h_nc_rebin->Scale(data_scale_factor);
	h_nc_pi0_rebin->Scale(data_scale_factor);
	h_other_mixed_rebin->Scale(data_scale_factor);
	h_unmatched_rebin->Scale(data_scale_factor);
	h_intime_rebin->Scale(intime_scale_factor);
	h_dirt_rebin->Scale(dirt_scale_factor);

	const double integral_data = h_data_rebin->Integral();

	if(area_norm)
	{
		const double integral_mc_ext = h_nue_cc_rebin->Integral() +
		                               h_nue_cc_mixed_rebin->Integral() +
		                               h_nue_cc_out_fv_rebin->Integral() +
		                               h_cosmic_rebin->Integral() +
		                               h_numu_cc_rebin->Integral() +
		                               h_nc_rebin->Integral() +
		                               h_nc_pi0_rebin->Integral() +
		                               h_other_mixed_rebin->Integral() +
		                               h_unmatched_rebin->Integral() +
		                               h_dirt_rebin->Integral();                                                                                          // +
		//h_intime_rebin->Integral();

		TH1 * h_data_scaling_clone = (TH1*)h_data_rebin->Clone("h_data_scaling_clone");
		h_data_scaling_clone->Add(h_intime_rebin, -1);
		const double integral_on_minus_off = h_data_scaling_clone->Integral();

		h_nue_cc_rebin->Scale(integral_on_minus_off / integral_mc_ext);
		h_nue_cc_mixed_rebin->Scale(integral_on_minus_off / integral_mc_ext);
		h_nue_cc_out_fv_rebin->Scale(integral_on_minus_off / integral_mc_ext);
		h_cosmic_rebin->Scale(integral_on_minus_off / integral_mc_ext);
		h_numu_cc_rebin->Scale(integral_on_minus_off / integral_mc_ext);
		h_nc_rebin->Scale(integral_on_minus_off / integral_mc_ext);
		h_nc_pi0_rebin->Scale(integral_on_minus_off / integral_mc_ext);
		h_other_mixed_rebin->Scale(integral_on_minus_off / integral_mc_ext);
		h_unmatched_rebin->Scale(integral_on_minus_off / integral_mc_ext);
		h_dirt_rebin->Scale(integral_on_minus_off / integral_mc_ext);

		// h_nue_cc_rebin->Scale(integral_data / integral_mc_ext);
		// h_nue_cc_mixed_rebin->Scale(integral_data / integral_mc_ext);
		// h_nue_cc_out_fv_rebin->Scale(integral_data / integral_mc_ext);
		// h_cosmic_rebin->Scale(integral_data / integral_mc_ext);
		// h_numu_cc_rebin->Scale(integral_data / integral_mc_ext);
		// h_nc_rebin->Scale(integral_data / integral_mc_ext);
		// h_nc_pi0_rebin->Scale(integral_data / integral_mc_ext);
		// h_other_mixed_rebin->Scale(integral_data / integral_mc_ext);
		// h_unmatched_rebin->Scale(integral_data / integral_mc_ext);
		//h_intime_rebin->Scale(integral_data / integral_mc_ext);

		h_nue_cc_rebin->Scale(1. / integral_data);
		h_nue_cc_mixed_rebin->Scale(1. / integral_data);
		h_nue_cc_out_fv_rebin->Scale(1. / integral_data);
		h_cosmic_rebin->Scale(1. / integral_data);
		h_numu_cc_rebin->Scale(1. / integral_data);
		h_nc_rebin->Scale(1. / integral_data);
		h_nc_pi0_rebin->Scale(1. / integral_data);
		h_other_mixed_rebin->Scale(1. / integral_data);
		h_unmatched_rebin->Scale(1. / integral_data);
		h_intime_rebin->Scale(1. / integral_data);
		h_data_rebin->Scale(1. / integral_data);
		h_dirt_rebin->Scale(1. / integral_data);

	}

	stack->Add(h_nue_cc_rebin);
	stack->Add(h_nue_cc_mixed_rebin);
	stack->Add(h_nue_cc_out_fv_rebin);
	stack->Add(h_cosmic_rebin);
	stack->Add(h_numu_cc_rebin);
	stack->Add(h_nc_rebin);
	stack->Add(h_nc_pi0_rebin);
	stack->Add(h_other_mixed_rebin);
	stack->Add(h_unmatched_rebin);
	stack->Add(h_dirt_rebin);
	stack->Add(h_intime_rebin);

	const double y_maximum = std::max(h_data_rebin->GetMaximum(), stack->GetMaximum());
	stack->SetMaximum(y_maximum * 1.2);

	stack->Draw("hist");
	//stack->GetYaxis()->SetTitle(y_axis_name);
	if(!area_norm) {stack->GetYaxis()->SetTitle("Entries"); }
	if(area_norm) {stack->GetYaxis()->SetTitle("Entries [A.U.]"); }
	stack->GetXaxis()->SetLabelOffset(10);
	stack->GetYaxis()->SetTitleSize(17);
	stack->GetYaxis()->SetTitleFont(46);
	stack->GetYaxis()->SetTitleOffset(4);

	h_data_rebin->Draw("same PE");

	TH1 * h_error_hist = (TH1*)h_nue_cc_rebin->Clone("h_error_hist");
	h_error_hist->Add(h_nue_cc_mixed_rebin,  1);
	h_error_hist->Add(h_nue_cc_out_fv_rebin, 1);
	h_error_hist->Add(h_numu_cc_rebin,       1);
	h_error_hist->Add(h_nc_pi0_rebin,        1);
	h_error_hist->Add(h_nc_rebin,            1);
	h_error_hist->Add(h_other_mixed_rebin,   1);
	h_error_hist->Add(h_cosmic_rebin,        1);
	h_error_hist->Add(h_unmatched_rebin,     1);
	h_error_hist->Add(h_dirt_rebin,          1);
	h_error_hist->Add(h_intime_rebin,        1);

	h_error_hist->SetFillColorAlpha(12, 0.15);
	h_error_hist->Draw("e2 hist same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(0.75, 0.98, 0.98, 0.50);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "#nu_{e} CC",        "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "#nu_{e} CC Mixed",  "f");
	leg_stack->AddEntry(h_nue_cc_out_fv,   "#nu_{e} CC OutFV",  "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic",            "f");
	leg_stack->AddEntry(h_numu_cc,         "#nu_{#mu} CC",      "f");
	leg_stack->AddEntry(h_nc,              "NC",                "f");
	leg_stack->AddEntry(h_nc_pi0,          "NC #pi^{0}",        "f");
	leg_stack->AddEntry(h_other_mixed,     "NC Mixed",          "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",         "f");
	leg_stack->AddEntry(h_dirt,            "Dirt",              "f");
	leg_stack->AddEntry(h_intime,          "InTime",            "f");
	leg_stack->Draw();

	TH1 * h_last = (TH1*) stack->GetStack()->Last();
	std::vector <double> chi2  = Chi2Calc(h_last, h_data_rebin, area_norm, integral_data);
	//chi2 : chi2/ndf, mc+ext, data

	//x_min, y_min, x_max, y_max
	//reduced chi2
	TPaveText * pt = new TPaveText(.46,.80,.73,1.06, "NBNDC");
	std::ostringstream o_string;
	o_string << float(chi2.at(0));
	std::string convert_string = o_string.str();
	std::string chi2_string = "#chi_{Stat}^{2}/DOF=" + convert_string;
	pt->AddText(chi2_string.c_str());
	pt->SetFillStyle(0);
	pt->SetBorderSize(0);
	//pt->Draw();

	//num events
	TPaveText * pt2 = new TPaveText(.13,.80,.46,1.06, "NBNDC");
	std::ostringstream o_string2a;
	std::ostringstream o_string2b;
	if(!area_norm)
	{
		o_string2a << int(chi2.at(2));
		o_string2b << int(chi2.at(1));
	}
	if(area_norm)
	{
		o_string2a << int(chi2.at(2) * integral_data);
		o_string2b << int(chi2.at(1) * integral_data);
	}
	std::string convert_string2a = o_string2a.str();
	std::string convert_string2b = o_string2b.str();
	std::string chi2_string2 = "Data: " + convert_string2a + "|MC+EXT:" + convert_string2b;
	pt2->AddText(chi2_string2.c_str());
	pt2->SetFillStyle(0);
	pt2->SetBorderSize(0);
	//pt2->Draw();

	//degrees of freedom
	TPaveText * pt3 = new TPaveText(.60,.80,.73,.973, "NBNDC");
	std::ostringstream o_string3;
	o_string3 << int(chi2.at(3));
	std::string convert_string3 = o_string3.str();
	std::string ndf_string = "DOF=" + convert_string3;
	pt3->AddText(ndf_string.c_str());
	pt3->SetFillStyle(0);
	pt3->SetBorderSize(0);
	//pt3->Draw();

	//p value
	//optional
	TPaveText * pt4 = new TPaveText(.45,.80,.60,.973, "NBNDC");
	std::ostringstream o_string4;
	o_string4.precision(4);
	o_string4 << std::fixed;
	o_string4 << chi2.at(4);
	std::string convert_string4 = o_string4.str();
	std::string p_string = "P=" + convert_string4;
	pt4->AddText(p_string.c_str());
	pt4->SetFillStyle(0);
	pt4->SetBorderSize(0);
	if(p_value) {pt4->Draw(); }

	bottomPad->cd();
	TH1 * ratioPlot = (TH1*)h_data_rebin->Clone("ratioPlot");
	// ratioPlot->Add(h_nue_cc_clone,        -1);
	// ratioPlot->Add(h_nue_cc_mixed_clone,  -1);
	// ratioPlot->Add(h_nue_cc_out_fv_clone, -1);
	// ratioPlot->Add(h_cosmic_clone,        -1);
	// ratioPlot->Add(h_numu_cc_clone,       -1);
	// ratioPlot->Add(h_nc_clone,            -1);
	// ratioPlot->Add(h_nc_pi0_clone,        -1);
	// ratioPlot->Add(h_other_mixed_clone,   -1);
	// ratioPlot->Add(h_unmatched_clone,     -1);
	// ratioPlot->Add(h_intime_clone,        -1);
	//ratioPlot->Divide(h_data);
	TH1 * h_mc_ext_sum = (TH1*)h_nue_cc_rebin->Clone("h_mc_ext_sum");
	//h_mc_ext_sum->Add(h_nue_cc_clone,        1);
	h_mc_ext_sum->Add(h_nue_cc_mixed_rebin,  1);
	h_mc_ext_sum->Add(h_nue_cc_out_fv_rebin, 1);
	h_mc_ext_sum->Add(h_cosmic_rebin,        1);
	h_mc_ext_sum->Add(h_numu_cc_rebin,       1);
	h_mc_ext_sum->Add(h_nc_rebin,            1);
	h_mc_ext_sum->Add(h_nc_pi0_rebin,        1);
	h_mc_ext_sum->Add(h_other_mixed_rebin,   1);
	h_mc_ext_sum->Add(h_unmatched_rebin,     1);
	h_mc_ext_sum->Add(h_dirt_rebin,          1);
	h_mc_ext_sum->Add(h_intime_rebin,        1);

	ratioPlot->GetXaxis()->SetLabelSize(12);
	ratioPlot->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	ratioPlot->GetYaxis()->SetLabelSize(11);
	ratioPlot->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	ratioPlot->GetXaxis()->SetTitleOffset(3.6);
	ratioPlot->GetXaxis()->SetTitleSize(15);
	ratioPlot->GetXaxis()->SetTitleFont(43);

	ratioPlot->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

	ratioPlot->Add(h_mc_ext_sum, -1);
	ratioPlot->Divide(h_mc_ext_sum);
	ratioPlot->GetYaxis()->SetRangeUser(-1,1);
	ratioPlot->GetXaxis()->SetTitle(x_axis_name);
	ratioPlot->GetYaxis()->SetTitle("(Data - MC) / MC ");
	ratioPlot->GetYaxis()->SetTitleSize(17);
	ratioPlot->GetYaxis()->SetTitleFont(46);
	ratioPlot->GetYaxis()->SetTitleOffset(2);
	ratioPlot->SetTitle(" ");
	ratioPlot->Draw();

	//now doing this stuff on the bottom pad

	//x_min, y_min, x_max, y_max
	//reduced chi2
	TPaveText * pt_bottom = new TPaveText(.12, .80, .30, .96, "NBNDC");
	std::ostringstream o_string_bottom;
	o_string_bottom.precision(3);
	o_string_bottom << std::fixed;
	o_string_bottom << float(chi2.at(0) * chi2.at(3));
	std::string convert_string_bottom = o_string_bottom.str();

	std::ostringstream o_string3_bottom;
	o_string3_bottom << int(chi2.at(3));
	std::string convert_string3_bottom = o_string3_bottom.str();

	std::string chi2_string_bottom = "#chi_{Stat}^{2}/DOF=(" + convert_string_bottom + "/" + convert_string3_bottom + ")";
	pt_bottom->AddText(chi2_string_bottom.c_str());
	pt_bottom->SetFillStyle(0);
	pt_bottom->SetBorderSize(0);
	pt_bottom->Draw();

	c1->Print(print_name);

	delete h_nue_cc_rebin;
	delete h_nue_cc_mixed_rebin;
	delete h_nue_cc_out_fv_rebin;
	delete h_cosmic_rebin;
	delete h_numu_cc_rebin;
	delete h_nc_rebin;
	delete h_nc_pi0_rebin;
	delete h_other_mixed_rebin;
	delete h_unmatched_rebin;
	delete h_intime_rebin;
	delete h_dirt_rebin;
	delete h_data_rebin;

	delete h_error_hist;
	delete h_last;

	delete ratioPlot;
	delete h_mc_ext_sum;

	delete c1;
	delete pt_bottom;
}

void histogram_functions::PlotdEdxTheta(
        TH2 * h_nue_cc, TH2 * h_nue_cc_mixed, TH2 * h_nue_cc_out_fv, TH2 * h_numu_cc,
        TH2 * h_numu_cc_mixed, TH2 * h_cosmic, TH1 * h_nc,
        TH2 * h_nc_pi0, TH2 * h_other_mixed, TH2 * h_unmatched, TH2 * h_intime, const double intime_scale_factor,
        TH2 * h_data, const double data_scale_factor,
        const char * title, const char * x_axis_name, const char * y_axis_name,
        const char * print_name1, const char * print_name2, const char * print_name3
        )
{

	TCanvas * c1 = new TCanvas();
	c1->cd();

	h_nue_cc->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_numu_cc->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_numu_cc_mixed->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_intime->SetStats(kFALSE);
	h_data->SetStats(kFALSE);
	h_nue_cc_out_fv->SetStats(kFALSE);

	TH2 * h_nue_cc_clone = (TH2*)h_nue_cc->Clone("h_nue_cc_clone");
	TH2 * h_nue_cc_mixed_clone = (TH2*)h_nue_cc_mixed->Clone("h_nue_cc_mixed_clone");
	TH2 * h_nue_cc_out_fv_clone = (TH2*)h_nue_cc_out_fv->Clone("h_nue_cc_out_fv_clone");
	TH2 * h_cosmic_clone = (TH2*)h_cosmic->Clone("h_cosmic_clone");
	TH2 * h_numu_cc_clone = (TH2*)h_numu_cc->Clone("h_numu_cc_clone");
	TH2 * h_nc_clone = (TH2*)h_nc->Clone("h_nc_clone");
	TH2 * h_nc_pi0_clone = (TH2*)h_nc_pi0->Clone("h_nc_pi0_clone");
	TH2 * h_other_mixed_clone = (TH2*)h_other_mixed->Clone("h_other_mixed_clone");
	TH2 * h_unmatched_clone = (TH2*)h_unmatched->Clone("h_unmatched_clone");
	TH2 * h_intime_clone = (TH2*)h_intime->Clone("h_intime_clone");

	h_nue_cc_clone->Scale(data_scale_factor);


	TH2 * h_dummy_mc = (TH2*)h_nue_cc_clone->Clone("h_dummy_mc");

	h_dummy_mc->Add(h_nue_cc_mixed_clone,  data_scale_factor);
	h_dummy_mc->Add(h_nue_cc_out_fv_clone, data_scale_factor);
	h_dummy_mc->Add(h_cosmic_clone,        data_scale_factor);
	h_dummy_mc->Add(h_numu_cc_clone,       data_scale_factor);
	h_dummy_mc->Add(h_nc_clone,            data_scale_factor);
	h_dummy_mc->Add(h_nc_pi0_clone,        data_scale_factor);
	h_dummy_mc->Add(h_other_mixed_clone,   data_scale_factor);
	h_dummy_mc->Add(h_unmatched_clone,     data_scale_factor);
	h_dummy_mc->Add(h_intime_clone,        intime_scale_factor);

	h_dummy_mc->GetXaxis()->SetTitle("Leading Shower dE/dx [MeV/cm]");
	h_dummy_mc->GetYaxis()->SetTitle("Leading Shower Theta [Degrees]");
	h_dummy_mc->Draw("colz");
	c1->Print(print_name1);

	TCanvas * c2 = new TCanvas();
	c2->cd();
	h_data->SetStats(kFALSE);
	TH2 * h_data_clone = (TH2*)h_data->Clone("h_data_clone");
	h_data_clone->GetXaxis()->SetTitle("Leading Shower dE/dx [MeV/cm]");
	h_data_clone->GetYaxis()->SetTitle("Leading Shower Theta [Degrees]");
	h_data_clone->Draw("colz");
	c2->Print(print_name2);

	TCanvas * c3 = new TCanvas();
	c3->cd();
	TH2 * h_division = (TH2*)h_data_clone->Clone("h_division");
	h_division->Divide(h_dummy_mc);
	h_division->Draw("colz");
	c3->Print(print_name3);

	delete h_nue_cc_clone;
	delete h_nue_cc_mixed_clone;
	delete h_nue_cc_out_fv_clone;
	delete h_cosmic_clone;
	delete h_numu_cc_clone;
	delete h_nc_clone;
	delete h_nc_pi0_clone;
	delete h_other_mixed_clone;
	delete h_unmatched_clone;
	delete h_intime_clone;

	delete c1;
	delete c2;
	delete c3;

}
void histogram_functions::PlotDetailStack(TH1 * h_nue_cc_qe,
                                          TH1 * h_nue_cc_out_fv,
                                          TH1 * h_nue_cc_res,
                                          TH1 * h_nue_cc_dis,
                                          TH1 * h_nue_cc_coh,
                                          TH1 * h_nue_cc_mec,
                                          TH1 * h_nue_cc_mixed,
                                          TH1 * h_numu_cc_qe,
                                          TH1 * h_numu_cc_res,
                                          TH1 * h_numu_cc_dis,
                                          TH1 * h_numu_cc_coh,
                                          TH1 * h_numu_cc_mec,
                                          TH1 * h_numu_cc_mixed,
                                          TH1 * h_cosmic,
                                          TH1 * h_nc,
                                          TH1 * h_nc_pi0,
                                          TH1 * h_other_mixed,
                                          TH1 * h_unmatched,
                                          const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	THStack * stack = new THStack();
	h_nue_cc_qe->SetStats(kFALSE);
	h_nue_cc_out_fv->SetStats(kFALSE);
	h_nue_cc_res->SetStats(kFALSE);
	h_nue_cc_dis->SetStats(kFALSE);
	h_nue_cc_coh->SetStats(kFALSE);
	h_nue_cc_mec->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_numu_cc_qe->SetStats(kFALSE);
	h_numu_cc_res->SetStats(kFALSE);
	h_numu_cc_dis->SetStats(kFALSE);
	h_numu_cc_coh->SetStats(kFALSE);
	h_numu_cc_mec->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_numu_cc_mixed->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_nue_cc_qe->SetFillColor(30);
	h_nue_cc_out_fv->SetFillColor(45);
	h_nue_cc_res->SetFillColor(31);
	h_nue_cc_dis->SetFillColor(32);
	h_nue_cc_coh->SetFillColor(33);
	h_nue_cc_mec->SetFillColor(34);
	h_nc->SetFillColor(46);
	h_numu_cc_qe->SetFillColor(28);
	h_numu_cc_res->SetFillColor(27);
	h_numu_cc_dis->SetFillColor(26);
	h_numu_cc_coh->SetFillColor(23);
	h_numu_cc_mec->SetFillColor(22);
	h_nc_pi0->SetFillColor(36);
	h_nue_cc_mixed->SetFillColor(38);
	h_numu_cc_mixed->SetFillColor(25);
	h_cosmic->SetFillColor(39);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	stack->Add(h_nue_cc_qe    );
	stack->Add(h_nue_cc_out_fv);
	stack->Add(h_nue_cc_res   );
	stack->Add(h_nue_cc_dis   );
	stack->Add(h_nue_cc_coh   );
	stack->Add(h_nue_cc_mec   );
	stack->Add(h_nc       );
	stack->Add(h_numu_cc_qe   );
	stack->Add(h_numu_cc_res  );
	stack->Add(h_numu_cc_dis  );
	stack->Add(h_numu_cc_coh  );
	stack->Add(h_numu_cc_mec  );
	stack->Add(h_nc_pi0      );
	stack->Add(h_nue_cc_mixed );
	stack->Add(h_numu_cc_mixed);
	stack->Add(h_cosmic       );
	stack->Add(h_other_mixed  );
	stack->Add(h_unmatched    );
	stack->Draw();
	stack->GetXaxis()->SetTitle(x_axis_name);
	TLegend * leg_track_stack_l1 = new TLegend(0.75, 0.50, 0.95, 0.95);
	leg_track_stack_l1->AddEntry(h_nue_cc_qe,      "Nue CC QE", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_out_fv,  "Nue CC Out FV", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_res,     "Nue CC Res", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_dis,     "Nue CC DIS", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_coh,     "Nue CC Coh", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_mec,     "Nue CC MEC", "f");
	leg_track_stack_l1->AddEntry(h_nc,         "NC", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_qe,     "Numu CC QE", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_res,    "Numu CC Res", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_dis,    "Numu CC DIS", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_coh,    "Numu CC Coh", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_mec,    "Numu CC MEC", "f");
	leg_track_stack_l1->AddEntry(h_nc_pi0,        "NC Pi0", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_mixed,   "Nue CC Mixed", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_mixed,  "Numu CC Mixed", "f");
	leg_track_stack_l1->AddEntry(h_cosmic,         "Cosmic", "f");
	leg_track_stack_l1->AddEntry(h_other_mixed,    "Other Mixed", "f");
	leg_track_stack_l1->AddEntry(h_unmatched,      "Unmatched", "f");
	leg_track_stack_l1->Draw();
	c1->Print(print_name);

	delete stack;
	delete c1;

}

void histogram_functions::EnergyOverlay(
        TH1 * h_0,
        TH1 * h_1,
        TH1 * h_2,
        TH1 * h_3,
        TH1 * h_4,
        TH1 * h_5,
        TH1 * h_6,
        TH1 * h_7,
        TH1 * h_8,
        TH1 * h_9,
        TH1 * h_10,
        TH1 * h_11,
        TH1 * h_12,
        TH1 * h_13,
        const char * x_axis_name, const char * y_axis_name, const char * print_name
        )
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	h_0->SetStats(kFALSE);
	h_0->SetFillColor(1);
	h_0->GetXaxis()->SetTitle(x_axis_name);
	const double nu_energy_no_cut_integral = h_0->Integral();
	h_0->Scale(1./nu_energy_no_cut_integral);
	h_0->Draw("hist");
	h_1->SetStats(kFALSE);
	h_1->SetFillColor(29);
	h_1->Scale(1./nu_energy_no_cut_integral);
	h_1->Draw("hist same");
	h_2->SetFillColor(30);
	h_2->SetStats(kFALSE);
	h_2->Scale(1./nu_energy_no_cut_integral);
	h_2->Draw("hist same");
	h_3->SetFillColor(45);
	h_3->SetStats(kFALSE);
	h_3->Scale(1./nu_energy_no_cut_integral);
	h_3->Draw("hist same");
	h_4->SetFillColor(28);
	h_4->SetStats(kFALSE);
	h_4->Scale(1./nu_energy_no_cut_integral);
	h_4->Draw("hist same");
	h_5->SetFillColor(26);
	h_5->SetStats(kFALSE);
	h_5->Scale(1./nu_energy_no_cut_integral);
	h_5->Draw("hist same");
	h_6->SetFillColor(36);
	h_6->SetStats(kFALSE);
	h_6->Scale(1./nu_energy_no_cut_integral);
	h_6->Draw("hist same");
	h_7->SetFillColor(39);
	h_7->SetStats(kFALSE);
	h_7->Scale(1./nu_energy_no_cut_integral);
	h_7->Draw("hist same");
	h_8->SetFillColor(42);
	h_8->SetStats(kFALSE);
	h_8->Scale(1./nu_energy_no_cut_integral);
	h_8->Draw("hist same");
	h_9->SetFillColor(12);
	h_9->SetStats(kFALSE);
	h_9->Scale(1./nu_energy_no_cut_integral);
	h_9->Draw("hist same");
	h_10->SetFillColor(1);
	h_10->SetStats(kFALSE);
	h_10->Scale(1./nu_energy_no_cut_integral);
	h_10->Draw("hist same");
	h_11->SetFillColor(46);
	h_11->SetStats(kFALSE);
	h_11->Scale(1./nu_energy_no_cut_integral);
	h_11->Draw("hist same");
	h_12->SetFillColor(19);
	h_12->SetStats(kFALSE);
	h_12->Scale(1./nu_energy_no_cut_integral);
	h_12->Draw("hist same");
	h_13->SetFillColor(41);
	h_13->SetStats(kFALSE);
	h_13->Scale(1./nu_energy_no_cut_integral);
	h_13->Draw("hist same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_sequential1 = new TLegend(0.65,0.65,0.85,0.85);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_sequential1->AddEntry(h_0,   "No Cut",           "f");
	leg_sequential1->AddEntry(h_1,   "Reco Nue",         "f");
	leg_sequential1->AddEntry(h_2,   "Fiducial Volume",  "f");
	leg_sequential1->AddEntry(h_3,   "Vtx To Flash",     "f");
	leg_sequential1->AddEntry(h_4,   "Shower Vtx",       "f");
	leg_sequential1->AddEntry(h_5,   "Track Vtx",        "f");
	leg_sequential1->AddEntry(h_6,   "Total Hits",       "f");
	leg_sequential1->AddEntry(h_7,   "Collection Hits",  "f");
	leg_sequential1->AddEntry(h_8,   "Open Angle",       "f");
	leg_sequential1->AddEntry(h_9,   "dE/dx",            "f");
	leg_sequential1->AddEntry(h_10,  "Secondary Shower", "f");
	leg_sequential1->AddEntry(h_11,  "Hits/Length",      "f");
	leg_sequential1->AddEntry(h_12,  "TrkShwr Length",   "f");
	leg_sequential1->AddEntry(h_13,  "Trk Containment",  "f");
	leg_sequential1->Draw();
	c1->Print(print_name);
}


void histogram_functions::PostHistogramOverlay(TH1 * h_no_cut, TH1 * h_reco_nue,
                                               TH1 * h_in_fv, TH1 * h_vtx_flash, TH1 * h_shwr_vtx,
                                               TH1 * h_trk_vtx, TH1 * h_hit_threshold,
                                               TH1 * h_open_angle, TH1 * h_dedx,
                                               const char * x_axis_name, const char * y_axis_name, const char * print_name)
{

	TCanvas * c1 = new TCanvas();
	c1->cd();
	h_no_cut->SetStats(kFALSE);
	h_no_cut->SetFillColor(29);
	h_no_cut->GetXaxis()->SetTitle(x_axis_name);
	const double nu_energy_no_cut_integral = h_no_cut->Integral();
	h_no_cut->Scale(1./nu_energy_no_cut_integral);
	h_no_cut->Draw("hist");
	h_reco_nue->SetFillColor(30);
	h_reco_nue->SetStats(kFALSE);
	h_reco_nue->Scale(1./nu_energy_no_cut_integral);
	h_reco_nue->Draw("hist same");
	h_in_fv->SetFillColor(45);
	h_in_fv->SetStats(kFALSE);
	h_in_fv->Scale(1./nu_energy_no_cut_integral);
	h_in_fv->Draw("hist same");
	h_vtx_flash->SetFillColor(28);
	h_vtx_flash->SetStats(kFALSE);
	h_vtx_flash->Scale(1./nu_energy_no_cut_integral);
	h_vtx_flash->Draw("hist same");
	h_shwr_vtx->SetFillColor(26);
	h_shwr_vtx->SetStats(kFALSE);
	h_shwr_vtx->Scale(1./nu_energy_no_cut_integral);
	h_shwr_vtx->Draw("hist same");
	h_trk_vtx->SetFillColor(36);
	h_trk_vtx->SetStats(kFALSE);
	h_trk_vtx->Scale(1./nu_energy_no_cut_integral);
	h_trk_vtx->Draw("hist same");
	h_hit_threshold->SetFillColor(39);
	h_hit_threshold->SetStats(kFALSE);
	h_hit_threshold->Scale(1./nu_energy_no_cut_integral);
	h_hit_threshold->Draw("hist same");
	h_open_angle->SetFillColor(42);
	h_open_angle->SetStats(kFALSE);
	h_open_angle->Scale(1./nu_energy_no_cut_integral);
	h_open_angle->Draw("hist same");
	h_dedx->SetFillColor(12);
	h_dedx->SetStats(kFALSE);
	h_dedx->Scale(1./nu_energy_no_cut_integral);
	h_dedx->Draw("hist same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_sequential1 = new TLegend(0.65,0.65,0.85,0.85);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_sequential1->AddEntry(h_no_cut,         "No Cuts",  "f");
	leg_sequential1->AddEntry(h_reco_nue,       "Reco Nue", "f");
	leg_sequential1->AddEntry(h_in_fv,          "Fiducial Volume", "f");
	leg_sequential1->AddEntry(h_vtx_flash,      "Vtx To Flash", "f");
	leg_sequential1->AddEntry(h_shwr_vtx,       "Shower Vtx", "f");
	leg_sequential1->AddEntry(h_trk_vtx,        "Track Vtx", "f");
	leg_sequential1->AddEntry(h_hit_threshold,  "Hits", "f");
	leg_sequential1->AddEntry(h_open_angle,     "Open Angle", "f");
	leg_sequential1->AddEntry(h_dedx,           "dE/dx", "f");
	leg_sequential1->Draw();
	c1->Print(print_name);
}

void histogram_functions::PlotDataMC(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                     TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_data, const double data_mc_scale_factor,
                                     const char * x_axis_name, const char * y_axis_name, const char * print_name, const char * data_print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	THStack * stack = new THStack();
	h_nue_cc->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_numu_cc->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_numu_cc_mixed->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_nue_cc->SetFillColor(30);
	h_nue_cc_mixed->SetFillColor(38);
	h_numu_cc->SetFillColor(28);
	h_nc_pi0->SetFillColor(36);
	h_cosmic->SetFillColor(1);
	h_nc->SetFillColor(46);
	h_numu_cc_mixed->SetFillColor(kTeal);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	stack->Add(h_nue_cc);
	stack->Add(h_nue_cc_mixed);
	stack->Add(h_cosmic);
	stack->Add(h_numu_cc);
	stack->Add(h_numu_cc_mixed);
	stack->Add(h_nc);
	stack->Add(h_nc_pi0);
	stack->Add(h_other_mixed);
	stack->Add(h_unmatched);
	stack->Draw();
	stack->GetXaxis()->SetTitle(x_axis_name);

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "Nue CC", "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic", "f");
	leg_stack->AddEntry(h_numu_cc,         "Numu CC", "f");
	leg_stack->AddEntry(h_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack->AddEntry(h_nc,              "NC", "f");
	leg_stack->AddEntry(h_nc_pi0,          "NC Pi0", "f");
	leg_stack->AddEntry(h_other_mixed,     "Other Mixed", "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched", "f");
	leg_stack->Draw();
	c1->Print(print_name);

	double integral_nue_cc          = h_nue_cc->Integral();
	double integral_nue_cc_mixed    = h_nue_cc_mixed->Integral();
	double integral_cosmic          = h_cosmic->Integral();
	double integral_numu_cc         = h_numu_cc->Integral();
	double integral_numu_cc_mixed   = h_numu_cc_mixed->Integral();
	double integral_nc              = h_nc->Integral();
	double integral_nc_pi0          = h_nc_pi0->Integral();
	double integral_other_mixed     = h_other_mixed->Integral();
	double integral_unmatched       = h_unmatched->Integral();
	double integral                 = integral_nue_cc
	                                  + integral_nue_cc_mixed
	                                  + integral_cosmic
	                                  + integral_numu_cc
	                                  + integral_numu_cc_mixed
	                                  + integral_nc
	                                  + integral_nc_pi0
	                                  + integral_other_mixed
	                                  + integral_unmatched;
	//h_leading_shower_open_angle_data->Scale(1.29916); //- relative scaling (on-beam - off-beam) and POT
	//h_leading_shower_open_angle_data->Scale(5.619); //- POT scaling
	h_data->Scale((integral / h_data->Integral()) * data_mc_scale_factor);
	h_data->Draw("P E1 same");
	c1->Print(data_print_name);
}

void histogram_functions::PlotDetailDataMCStack(TH1 * h_nue_cc_qe,
                                                TH1 * h_nue_cc_out_fv,
                                                TH1 * h_nue_cc_res,
                                                TH1 * h_nue_cc_dis,
                                                TH1 * h_nue_cc_coh,
                                                TH1 * h_nue_cc_mec,
                                                TH1 * h_nue_cc_mixed,
                                                TH1 * h_numu_cc_qe,
                                                TH1 * h_numu_cc_res,
                                                TH1 * h_numu_cc_dis,
                                                TH1 * h_numu_cc_coh,
                                                TH1 * h_numu_cc_mec,
                                                TH1 * h_numu_cc_mixed,
                                                TH1 * h_cosmic,
                                                TH1 * h_nc,
                                                TH1 * h_nc_pi0,
                                                TH1 * h_other_mixed,
                                                TH1 * h_unmatched,
                                                TH1 * h_data,
                                                const double data_mc_scale_factor,
                                                const char * x_axis_name, const char * y_axis_name, const char * print_name, const char * data_print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	THStack * stack = new THStack();
	h_nue_cc_qe->SetStats(kFALSE);
	h_nue_cc_out_fv->SetStats(kFALSE);
	h_nue_cc_res->SetStats(kFALSE);
	h_nue_cc_dis->SetStats(kFALSE);
	h_nue_cc_coh->SetStats(kFALSE);
	h_nue_cc_mec->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_numu_cc_qe->SetStats(kFALSE);
	h_numu_cc_res->SetStats(kFALSE);
	h_numu_cc_dis->SetStats(kFALSE);
	h_numu_cc_coh->SetStats(kFALSE);
	h_numu_cc_mec->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_numu_cc_mixed->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	h_nue_cc_qe->SetFillColor(30);
	h_nue_cc_out_fv->SetFillColor(45);
	h_nue_cc_res->SetFillColor(31);
	h_nue_cc_dis->SetFillColor(32);
	h_nue_cc_coh->SetFillColor(33);
	h_nue_cc_mec->SetFillColor(34);
	h_nc->SetFillColor(46);
	h_numu_cc_qe->SetFillColor(28);
	h_numu_cc_res->SetFillColor(27);
	h_numu_cc_dis->SetFillColor(26);
	h_numu_cc_coh->SetFillColor(23);
	h_numu_cc_mec->SetFillColor(22);
	h_nc_pi0->SetFillColor(36);
	h_nue_cc_mixed->SetFillColor(38);
	h_numu_cc_mixed->SetFillColor(25);
	h_cosmic->SetFillColor(39);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	stack->Add(h_nue_cc_qe    );
	stack->Add(h_nue_cc_out_fv);
	stack->Add(h_nue_cc_res   );
	stack->Add(h_nue_cc_dis   );
	stack->Add(h_nue_cc_coh   );
	stack->Add(h_nue_cc_mec   );
	stack->Add(h_nc       );
	stack->Add(h_numu_cc_qe   );
	stack->Add(h_numu_cc_res  );
	stack->Add(h_numu_cc_dis  );
	stack->Add(h_numu_cc_coh  );
	stack->Add(h_numu_cc_mec  );
	stack->Add(h_nc_pi0      );
	stack->Add(h_nue_cc_mixed );
	stack->Add(h_numu_cc_mixed);
	stack->Add(h_cosmic       );
	stack->Add(h_other_mixed  );
	stack->Add(h_unmatched    );
	stack->Draw();
	stack->GetXaxis()->SetTitle(x_axis_name);
	TLegend * leg_track_stack_l1 = new TLegend(0.75, 0.50, 0.95, 0.95);
	leg_track_stack_l1->AddEntry(h_nue_cc_qe,      "Nue CC QE", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_out_fv,  "Nue CC Out FV", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_res,     "Nue CC Res", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_dis,     "Nue CC DIS", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_coh,     "Nue CC Coh", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_mec,     "Nue CC MEC", "f");
	leg_track_stack_l1->AddEntry(h_nc,         "NC", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_qe,     "Numu CC QE", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_res,    "Numu CC Res", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_dis,    "Numu CC DIS", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_coh,    "Numu CC Coh", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_mec,    "Numu CC MEC", "f");
	leg_track_stack_l1->AddEntry(h_nc_pi0,        "NC Pi0", "f");
	leg_track_stack_l1->AddEntry(h_nue_cc_mixed,   "Nue CC Mixed", "f");
	leg_track_stack_l1->AddEntry(h_numu_cc_mixed,  "Numu CC Mixed", "f");
	leg_track_stack_l1->AddEntry(h_cosmic,         "Cosmic", "f");
	leg_track_stack_l1->AddEntry(h_other_mixed,    "Other Mixed", "f");
	leg_track_stack_l1->AddEntry(h_unmatched,      "Unmatched", "f");
	leg_track_stack_l1->Draw();
	c1->Print(print_name);

	//h_pfp_track_data->Scale(1.2991626);
	h_data->Scale(data_mc_scale_factor);
	h_data->Draw("P E1 same");
	c1->Print(data_print_name);
}
void histogram_functions::PurityStack(TH1 * h_qe, TH1 * h_res, TH1 * h_dis, TH1 * h_coh, TH1 * h_mec, const char * x_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	THStack * stack = new THStack();
	h_qe->SetStats(kFALSE);
	h_res->SetStats(kFALSE);
	h_dis->SetStats(kFALSE);
	h_coh->SetStats(kFALSE);
	h_mec->SetStats(kFALSE);
	h_qe->SetFillColor(30);
	h_res->SetFillColor(38);
	h_dis->SetFillColor(28);
	h_coh->SetFillColor(1);
	h_mec->SetFillColor(36);
	stack->Add(h_qe);
	stack->Add(h_res);
	stack->Add(h_dis);
	stack->Add(h_coh);
	stack->Add(h_mec);
	stack->Draw();
	stack->GetXaxis()->SetTitle(x_axis_name);

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(0.75,0.75,0.95,0.95);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_qe,          "QE", "f");
	leg_stack->AddEntry(h_res,         "Res", "f");
	leg_stack->AddEntry(h_dis,         "DIS", "f");
	leg_stack->AddEntry(h_coh,         "Coh", "f");
	leg_stack->AddEntry(h_mec,         "MEC", "f");
	leg_stack->Draw();
	c1->Print(print_name);

	delete stack;
	delete c1;

}
void histogram_functions::OverlayScatter(TH2 * h_nue_cc, TH2 * h_nue_cc_mixed, TH2 * h_nue_cc_out_fv, TH2 * h_numu_cc,
                                         TH2 * h_numu_cc_mixed, TH2 * h_cosmic, TH2 * h_nc,
                                         TH2 * h_nc_pi0, TH2 * h_other_mixed, TH2 * h_unmatched,
                                         const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                         const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	h_nue_cc->SetStats(kFALSE);
	h_nue_cc_mixed->SetStats(kFALSE);
	h_nue_cc_out_fv->SetStats(kFALSE);
	h_numu_cc->SetStats(kFALSE);
	h_nc_pi0->SetStats(kFALSE);
	h_cosmic->SetStats(kFALSE);
	h_nc->SetStats(kFALSE);
	h_numu_cc_mixed->SetStats(kFALSE);
	h_other_mixed->SetStats(kFALSE);
	h_unmatched->SetStats(kFALSE);
	// h_nue_cc->SetFillColor(30);
	// h_nue_cc_mixed->SetFillColor(38);
	// h_numu_cc->SetFillColor(28);
	// h_nc_pi0->SetFillColor(36);
	// h_cosmic->SetFillColor(1);
	// h_nc->SetFillColor(46);
	// h_numu_cc_mixed->SetFillColor(kTeal);
	// h_other_mixed->SetFillColor(42);
	// h_unmatched->SetFillColor(12);
	h_nue_cc->SetMarkerColor(30);
	h_nue_cc_mixed->SetMarkerColor(38);
	h_numu_cc->SetMarkerColor(28);
	h_nc_pi0->SetMarkerColor(36);
	h_cosmic->SetMarkerColor(1);
	h_nc->SetMarkerColor(46);
	h_nue_cc_out_fv->SetMarkerColor(kTeal);
	h_numu_cc_mixed->SetMarkerColor(28);
	h_other_mixed->SetMarkerColor(42);
	h_unmatched->SetMarkerColor(12);
	h_nue_cc->SetMarkerStyle(8);
	h_nue_cc_mixed->SetMarkerStyle(8);
	h_nue_cc_out_fv->SetMarkerStyle(8);
	h_numu_cc->SetMarkerStyle(8);
	h_nc_pi0->SetMarkerStyle(8);
	h_cosmic->SetMarkerStyle(8);
	h_nc->SetMarkerStyle(8);
	h_numu_cc_mixed->SetMarkerStyle(8);
	h_other_mixed->SetMarkerStyle(8);
	h_unmatched->SetMarkerStyle(8);
	h_nue_cc->Draw();
	h_nue_cc->GetXaxis()->SetTitle(x_axis_name);
	h_nue_cc->GetYaxis()->SetTitle(y_axis_name);
	h_nue_cc_mixed->Draw("SAME");
	h_nue_cc_out_fv->Draw("SAME");
	h_cosmic->Draw("SAME");
	h_numu_cc->Draw("SAME");
	h_numu_cc_mixed->Draw("SAME");
	h_nc->Draw("SAME");
	h_nc_pi0->Draw("SAME");
	h_other_mixed->Draw("SAME");
	h_unmatched->Draw("SAME");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "Nue CC",        "p");
	leg_stack->AddEntry(h_nue_cc_mixed,    "Nue CC Mixed",  "p");
	leg_stack->AddEntry(h_nue_cc_out_fv,   "Nue CC OutFV",  "p");
	leg_stack->AddEntry(h_cosmic,          "Cosmic",        "p");
	leg_stack->AddEntry(h_numu_cc,         "Numu CC",       "p");
	//leg_stack->AddEntry(h_numu_cc_mixed,   "Numu CC Mixed", "p");
	leg_stack->AddEntry(h_nc,              "NC",            "p");
	leg_stack->AddEntry(h_nc_pi0,          "NC Pi0",        "p");
	leg_stack->AddEntry(h_other_mixed,     "NC Mixed",   "p");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",     "p");
	leg_stack->Draw();
	c1->Print(print_name);
}

void histogram_functions::Plot2DdEdxMap(TH2 * histogram, const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();
	histogram->GetXaxis()->SetTitle(x_axis_name);
	histogram->GetYaxis()->SetTitle(y_axis_name);
	histogram->SetTitle(title);
	histogram->SetStats(kFALSE);
	histogram->Scale( 1. / histogram->Integral());
	histogram->Draw("colz");
	c1->Print(print_name);
}

std::vector <double> histogram_functions::Chi2Calc(TH1 * h_mc_ext, TH1 * h_data, const bool area_norm, const double return_norm){
	const int n_bins = h_mc_ext->GetNbinsX();

	const double f_1 = h_mc_ext->Integral();
	const double f_2 = h_data->Integral();

	//area normalised?
	TH1 * h_mc_ext_clone = (TH1*)h_mc_ext->Clone("h_mc_ext_clone");
	TH1 * h_data_clone = (TH1*)h_data->Clone("h_data_clone");
	if(!area_norm) {h_mc_ext_clone->Scale(f_2/f_1); }
	if(area_norm)
	{
		//this keeps them area normalised,
		//but at the original values, not 0->1
		//which messes with the chi2 and p calc
		h_mc_ext_clone->Scale(return_norm);
		h_data_clone->Scale(return_norm);

		const double f_1_adj = h_mc_ext->Integral();
		const double f_2_adj = h_data->Integral();
		h_mc_ext_clone->Scale(f_2_adj/f_1_adj);

	}
	//h_data_clone->Scale(1./f_2);

	std::vector <double> chi2;
	double chi2_val = 0;
	double n_mc_ext_val = 0;
	double n_data_val = 0;
	for( int i = 1; i < n_bins; i++)
	{
		const double n_mc_ext = h_mc_ext_clone->GetBinContent(i);
		const double n_data   = h_data_clone->GetBinContent(i);

		//don't calculate chi2 for bins where no comparison possible
		if(n_data == 0 || n_mc_ext == 0) { continue; }

		//chi2_val += (pow((n_mc_ext - n_data),2) / n_mc_ext);
		//chi2_val += (pow((n_data - n_mc_ext),2)) / (((n_data * f_2) / pow(f_2, 2)) + ((n_mc_ext * f_1) / pow(f_1, 2)));
		chi2_val += 2 * (n_mc_ext - n_data + (n_data * TMath::Log(n_data/n_mc_ext)));

		n_mc_ext_val += n_mc_ext;
		n_data_val += n_data;
	}
	const double reduced_chi2 = chi2_val / (n_bins - 1);
	const double p = TMath::Prob(chi2_val, n_bins);

	chi2.push_back(reduced_chi2);
	//correct this value back to the un-normalised
	//chi2.push_back(n_mc_ext_val * f_1);
	//chi2.push_back(n_data_val * f_2);
	chi2.push_back(n_mc_ext_val * (f_1 / f_2));
	chi2.push_back(n_data_val);
	chi2.push_back(n_bins - 1);
	chi2.push_back(p);
	return chi2;
}
