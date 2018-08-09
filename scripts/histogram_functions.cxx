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
		std::cout << "----------" << std::endl;
		for(int i = 1; i < h_num_clone_rebin->GetNbinsX()+1; i++)
		{
			std::cout << "Bin: " << i << " Value: " << h_num_clone_rebin->GetBinContent(i) / h_den_clone_rebin->GetBinContent(i) << std::endl;
		}
		std::cout << "----------" << std::endl;

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
void histogram_functions::TimingHistograms(TH1 * histogram_1, TH1 * histogram_2, TH1 * histogram_3,
                                           const double data_scale_factor, const double intime_scale_factor,
                                           const char * x_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();

	TH1 * h_1_clone = (TH1*)histogram_1->Clone("h_1_clone");
	TH1 * h_2_clone = (TH1*)histogram_2->Clone("h_2_clone");
	TH1 * h_3_clone = (TH1*)histogram_3->Clone("h_3_clone");

	h_1_clone->Sumw2();
	h_2_clone->Sumw2();
	h_3_clone->Sumw2();

	h_1_clone->Scale(data_scale_factor);
	h_2_clone->Scale(intime_scale_factor);
	h_3_clone->Add(h_2_clone, -1);

	h_3_clone->GetXaxis()->SetTitle(x_axis_name);
	h_3_clone->GetYaxis()->SetRangeUser(-3000, 3000);

	h_3_clone->Draw("hist");
	h_1_clone->SetLineColor(46);
	h_1_clone->Draw("hist same");

	c1->Print(print_name);

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
	h_flash_time_data_divide->GetYaxis()->SetRangeUser(0.8, 1.5);
	h_flash_time_data_divide->GetXaxis()->SetTitle("Flash Time [#mus]");
	h_flash_time_data_divide->GetYaxis()->SetTitle("NuMI Run 1 On-Beam First Half / Second Half");
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
	//h_intime->SetFillSytle(3354);
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
	                      0.75, 0.95, 0.75, 0.95,
	                      title, x_axis_name, y_axis_name, print_name);
}
void histogram_functions::PlotSimpleStackData (TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                               TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                               TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                               TH1 * h_data, const double data_scale_factor,
                                               const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	//default legend position
	//0.75,0.75,0.95,0.95
	PlotSimpleStackData(h_nue_cc, h_nue_cc_mixed, h_nue_cc_out_fv, h_numu_cc, h_numu_cc_mixed, h_cosmic, h_nc,
	                    h_nc_pi0, h_other_mixed, h_unmatched, h_intime, intime_scale_factor,
	                    h_data, data_scale_factor,
	                    0.75, 0.95, 0.75, 0.95,
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
	h_neutron->SetFillColor(20);
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
	h_nue_cc_out_fv->SetFillColor(20);
	//h_numu_cc_mixed->SetFillColor(28);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	stack->Add(h_nue_cc);
	stack->Add(h_nue_cc_mixed);
	stack->Add(h_nue_cc_out_fv);
	stack->Add(h_cosmic);
	//h_numu_cc->Add(h_numu_cc_mixed, 1);
	stack->Add(h_numu_cc);
	//stack->Add(h_numu_cc_mixed);
	stack->Add(h_nc);
	stack->Add(h_nc_pi0);
	stack->Add(h_other_mixed);
	stack->Add(h_unmatched);
	stack->Draw();
	stack->GetXaxis()->SetTitle(x_axis_name);

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "Nue CC",       "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack->AddEntry(h_nue_cc_out_fv,   "Nue CC OutFV", "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic",       "f");
	leg_stack->AddEntry(h_numu_cc,         "Numu CC",      "f");
	//leg_stack->AddEntry(h_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack->AddEntry(h_nc,              "NC",           "f");
	leg_stack->AddEntry(h_nc_pi0,          "NC Pi0",       "f");
	leg_stack->AddEntry(h_other_mixed,     "NC Mixed",     "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",    "f");
	leg_stack->Draw();
	c1->Print(print_name);
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
	h_nue_cc_out_fv->SetFillColor(20);
	//h_numu_cc_mixed->SetFillColor(28);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	h_intime->SetFillColor(41);
	h_intime->SetFillStyle(3345);
	//h_intime->SetFillSytle(3354);
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
}
void histogram_functions::PlotSimpleStackData(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                              TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                              TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                              TH1 * h_data, const double data_scale_factor,
                                              const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                              const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas(title, title, 500, 500);
	c1->cd();

	TPad *topPad = new TPad("topPad", "", 0.005, 0.32, 0.995, 0.995);
	TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.28);
	topPad->SetBottomMargin(0.08);
	bottomPad->SetTopMargin(0.0);
	bottomPad->SetBottomMargin(0.12);
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
	h_nue_cc_out_fv->SetFillColor(20);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	h_intime->SetFillColor(41);
	h_intime->SetFillStyle(3345);
	//h_intime->SetFillSytle(3354);

	//h_data->Scale(data_scale_factor);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(0.5);
	h_data->Sumw2();

	TH1 * h_nue_cc_clone        = (TH1*)h_nue_cc->Clone("h_nue_cc_clone");
	TH1 * h_nue_cc_mixed_clone  = (TH1*)h_nue_cc_mixed->Clone("h_nue_cc_mixed_clone");
	TH1 * h_nue_cc_out_fv_clone = (TH1*)h_nue_cc_out_fv->Clone("h_nue_cc_out_fv_clone");
	TH1 * h_cosmic_clone        = (TH1*)h_cosmic->Clone("h_cosmic_clone");
	TH1 * h_numu_cc_clone       = (TH1*)h_numu_cc->Clone("h_numu_cc_clone");
	TH1 * h_nc_clone            = (TH1*)h_nc->Clone("h_nc_clone");
	TH1 * h_nc_pi0_clone        = (TH1*)h_nc_pi0->Clone("h_nc_pi0_clone");
	TH1 * h_other_mixed_clone   = (TH1*)h_other_mixed->Clone("h_other_mixed_clone");
	TH1 * h_unmatched_clone     = (TH1*)h_unmatched->Clone("h_unmatched_clone");
	TH1 * h_intime_clone = (TH1*)h_intime->Clone("h_intime_clone");

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
	h_data->Sumw2();

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

	const double y_maximum = std::max(h_data->GetMaximum(), stack->GetMaximum());
	stack->SetMaximum(y_maximum * 1.2);

	stack->Draw("hist");
	stack->GetXaxis()->SetTitle(x_axis_name);
	h_data->Draw("same PE");

	TH1 * h_error_hist = (TH1*)h_nue_cc_clone->Clone("h_error_hist");
	h_error_hist->Add(h_nue_cc_mixed_clone, 1);
	h_error_hist->Add(h_nue_cc_out_fv_clone, 1);
	h_error_hist->Add(h_numu_cc_clone, 1);
	h_error_hist->Add(h_nc_pi0_clone, 1);
	h_error_hist->Add(h_nc_clone, 1);
	h_error_hist->Add(h_other_mixed_clone, 1);
	h_error_hist->Add(h_cosmic_clone, 1);
	h_error_hist->Add(h_unmatched_clone, 1);
	h_error_hist->Add(h_intime_clone, 1);

	h_error_hist->SetFillColorAlpha(12, 0.15);
	h_error_hist->Draw("e2 hist same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "Nue CC",        "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "Nue CC Mixed",  "f");
	leg_stack->AddEntry(h_nue_cc_out_fv,   "Nue CC OutFV",  "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic",        "f");
	leg_stack->AddEntry(h_numu_cc,         "Numu CC",       "f");
	leg_stack->AddEntry(h_nc,              "NC",            "f");
	leg_stack->AddEntry(h_nc_pi0,          "NC Pi0",        "f");
	leg_stack->AddEntry(h_other_mixed,     "NC Mixed",      "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",     "f");
	leg_stack->AddEntry(h_intime,          "InTime",        "f");
	leg_stack->Draw();

	bottomPad->cd();
	TH1 * ratioPlot = (TH1*)h_data->Clone("ratioPlot");
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
	h_mc_ext_sum->Add(h_intime_clone,        1);

	ratioPlot->Add(h_mc_ext_sum, -1);
	ratioPlot->Divide(h_mc_ext_sum);
	ratioPlot->GetYaxis()->SetRangeUser(-1,1);
	ratioPlot->Draw();


	c1->Print(print_name);
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
	h_numu_cc_mixed->SetFillColor(20);
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
	// h_numu_cc_mixed->SetFillColor(20);
	// h_other_mixed->SetFillColor(42);
	// h_unmatched->SetFillColor(12);
	h_nue_cc->SetMarkerColor(30);
	h_nue_cc_mixed->SetMarkerColor(38);
	h_numu_cc->SetMarkerColor(28);
	h_nc_pi0->SetMarkerColor(36);
	h_cosmic->SetMarkerColor(1);
	h_nc->SetMarkerColor(46);
	h_nue_cc_out_fv->SetMarkerColor(20);
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
