#include "histogram_functions.h"

void histogram_functions::Plot1DHistogram (TH1 * histogram, const char * x_axis_name, const char * print_name)
{
	TCanvas * c1 = new TCanvas();
	c1->cd();

	histogram->GetXaxis()->SetTitle(x_axis_name);
	histogram->Draw();

	c1->Print(print_name);
}

void histogram_functions::PlotTEfficiency (TH1 *h_num, TH1 *h_den, const char * title, const char * print_name)
{

	TCanvas * efficiency_c1 = new TCanvas();
	efficiency_c1->cd();

	TEfficiency * teff = new TEfficiency(*h_num, *h_den);
	teff->SetTitle(title);
	teff->SetLineColor(kGreen+3);
	teff->SetMarkerColor(kGreen+3);
	teff->SetMarkerStyle(20);
	teff->SetMarkerSize(0.5);
	teff->Draw("AP");
	efficiency_c1->Print(print_name);
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
void histogram_functions::LegoStackData(TH2 * h_nue_cc, TH2 * h_nue_cc_mixed, TH2 * h_numu_cc, TH2 * h_numu_cc_mixed, TH2 * h_cosmic, TH2 * h_nc,
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
	h_numu_cc_mixed->SetStats(kFALSE);
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
	h_numu_cc_mixed->SetFillColor(20);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	h_intime->SetFillColor(41);
	h_intime->SetFillStyle(3345);
	//h_intime->SetFillSytle(3354);
	h_intime->Scale(intime_scale_factor);
	h_data->Scale(data_scale_factor);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(0.5);
	h_data->Sumw2();
	h_data->Scale(1./h_data->Integral());

	const double integral = h_nue_cc->Integral() +
	                        h_nue_cc_mixed->Integral() +
	                        h_cosmic->Integral() +
	                        h_numu_cc->Integral() +
	                        h_numu_cc_mixed->Integral() +
	                        h_nc->Integral() +
	                        h_nc_pi0->Integral() +
	                        h_other_mixed->Integral() +
	                        h_unmatched->Integral() +
	                        h_intime->Integral();

	h_nue_cc->Scale(1./integral);
	h_nue_cc_mixed->Scale(1./integral);
	h_cosmic->Scale(1./integral);
	h_numu_cc->Scale(1./integral);
	h_numu_cc_mixed->Scale(1./integral);
	h_nc->Scale(1./integral);
	h_nc_pi0->Scale(1./integral);
	h_other_mixed->Scale(1./integral);
	h_unmatched->Scale(1./integral);
	h_intime->Scale(1./integral);
	stack->Add(h_nue_cc);
	stack->Add(h_nue_cc_mixed);
	stack->Add(h_cosmic);
	stack->Add(h_numu_cc);
	stack->Add(h_numu_cc_mixed);
	stack->Add(h_nc);
	stack->Add(h_nc_pi0);
	stack->Add(h_other_mixed);
	stack->Add(h_unmatched);
	stack->Add(h_intime);

	// const double y_maximum = std::max(h_data->GetMaximum(), stack->GetMaximum());
	// stack->SetMaximum(y_maximum * 1.2);

	stack->Draw("lego4 0");
	stack->GetXaxis()->SetTitle(x_axis_name);
	//h_data->Draw("lego1 0 same");

	//gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	TLegend * leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	//leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg_stack->AddEntry(h_nue_cc,          "Nue CC",        "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "Nue CC Mixed",  "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic",        "f");
	leg_stack->AddEntry(h_numu_cc,         "Numu CC",       "f");
	leg_stack->AddEntry(h_numu_cc_mixed,   "Numu CC Mixed", "f");
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
                                                 TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                                 const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
{
	//default legend position
	//0.75,0.75,0.95,0.95
	PlotSimpleStackInTime(h_nue_cc, h_nue_cc_mixed, h_nue_cc_out_fv, h_numu_cc, h_numu_cc_mixed, h_cosmic, h_nc,
	                      h_nc_pi0, h_other_mixed, h_unmatched, h_intime, intime_scale_factor,
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
	h_numu_cc_mixed->SetFillColor(28);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	stack->Add(h_nue_cc);
	stack->Add(h_nue_cc_mixed);
	stack->Add(h_nue_cc_out_fv);
	stack->Add(h_cosmic);
	h_numu_cc->Add(h_numu_cc_mixed, 1);
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
	leg_stack->AddEntry(h_nue_cc,          "Nue CC", "f");
	leg_stack->AddEntry(h_nue_cc_mixed,    "Nue CC Mixed", "f");
	leg_stack->AddEntry(h_nue_cc_out_fv,   "Nue CC OutFV", "f");
	leg_stack->AddEntry(h_cosmic,          "Cosmic", "f");
	leg_stack->AddEntry(h_numu_cc,         "Numu CC", "f");
	//leg_stack->AddEntry(h_numu_cc_mixed,   "Numu CC Mixed", "f");
	leg_stack->AddEntry(h_nc,              "NC", "f");
	leg_stack->AddEntry(h_nc_pi0,          "NC Pi0", "f");
	leg_stack->AddEntry(h_other_mixed,     "Other Mixed", "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched", "f");
	leg_stack->Draw();
	c1->Print(print_name);
}
void histogram_functions::PlotSimpleStackInTime(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc,
                                                TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                                TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
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
	h_numu_cc_mixed->SetFillColor(28);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	h_intime->SetFillColor(41);
	h_intime->SetFillStyle(3345);
	//h_intime->SetFillSytle(3354);
	h_intime->Scale(intime_scale_factor);
	stack->Add(h_nue_cc);
	stack->Add(h_nue_cc_mixed);
	stack->Add(h_nue_cc_out_fv);
	stack->Add(h_cosmic);
	h_numu_cc->Add(h_numu_cc_mixed, 1);
	stack->Add(h_numu_cc);
	//stack->Add(h_numu_cc_mixed);
	stack->Add(h_nc);
	stack->Add(h_nc_pi0);
	stack->Add(h_other_mixed);
	stack->Add(h_unmatched);
	stack->Add(h_intime);
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
	leg_stack->AddEntry(h_other_mixed,     "Other Mixed",   "f");
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
	h_data->SetStats(kFALSE);
	h_nue_cc_out_fv->SetStats(kFALSE);
	h_nue_cc->SetFillColor(30);
	h_nue_cc_mixed->SetFillColor(38);
	h_numu_cc->SetFillColor(28);
	h_nc_pi0->SetFillColor(36);
	h_cosmic->SetFillColor(1);
	h_nc->SetFillColor(46);
	h_nue_cc_out_fv->SetFillColor(20);
	h_numu_cc_mixed->SetFillColor(28);
	h_other_mixed->SetFillColor(42);
	h_unmatched->SetFillColor(12);
	h_intime->SetFillColor(41);
	h_intime->SetFillStyle(3345);
	//h_intime->SetFillSytle(3354);
	h_intime->Scale(intime_scale_factor);
	h_data->Scale(data_scale_factor);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(0.5);
	h_data->Sumw2();
	h_data->Scale(1./h_data->Integral());
	//h_data->GetXaxis()->SetTitle(x_axis_name);

	const double integral = h_nue_cc->Integral() +
	                        h_nue_cc_mixed->Integral() +
	                        h_cosmic->Integral() +
	                        h_numu_cc->Integral() +
	                        h_numu_cc_mixed->Integral() +
	                        h_nc->Integral() +
	                        h_nc_pi0->Integral() +
	                        h_other_mixed->Integral() +
	                        h_unmatched->Integral() +
	                        h_intime->Integral() +
	                        h_nue_cc_out_fv->Integral();

	h_nue_cc->Scale(1./integral);
	h_nue_cc_mixed->Scale(1./integral);
	h_cosmic->Scale(1./integral);
	h_numu_cc->Scale(1./integral);
	h_numu_cc_mixed->Scale(1./integral);
	h_nc->Scale(1./integral);
	h_nc_pi0->Scale(1./integral);
	h_other_mixed->Scale(1./integral);
	h_unmatched->Scale(1./integral);
	h_intime->Scale(1./integral);
	h_nue_cc_out_fv->Scale(1./integral);
	stack->Add(h_nue_cc);
	stack->Add(h_nue_cc_mixed);
	stack->Add(h_nue_cc_out_fv);
	stack->Add(h_cosmic);
	h_numu_cc->Add(h_numu_cc_mixed, 1);
	stack->Add(h_numu_cc);
	//stack->Add(h_numu_cc_mixed);
	stack->Add(h_nc);
	stack->Add(h_nc_pi0);
	stack->Add(h_other_mixed);
	stack->Add(h_unmatched);
	stack->Add(h_intime);

	const double y_maximum = std::max(h_data->GetMaximum(), stack->GetMaximum());
	stack->SetMaximum(y_maximum * 1.2);

	//h_data->Draw();
	stack->Draw("hist");
	stack->GetXaxis()->SetTitle(x_axis_name);
	h_data->Draw("same");

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
	leg_stack->AddEntry(h_other_mixed,     "Other Mixed",   "f");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",     "f");
	leg_stack->AddEntry(h_intime,          "InTime",        "f");
	leg_stack->Draw();
	c1->Print(print_name);
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
void histogram_functions::OverlayScatter(TH2 * h_nue_cc, TH2 * h_nue_cc_mixed, TH2 * h_numu_cc, TH2 * h_numu_cc_mixed, TH2 * h_cosmic, TH2 * h_nc,
                                         TH2 * h_nc_pi0, TH2 * h_other_mixed, TH2 * h_unmatched,
                                         const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                         const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name)
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
	h_numu_cc_mixed->SetMarkerColor(20);
	h_other_mixed->SetMarkerColor(42);
	h_unmatched->SetMarkerColor(12);
	h_nue_cc->SetMarkerStyle(8);
	h_nue_cc_mixed->SetMarkerStyle(8);
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
	leg_stack->AddEntry(h_cosmic,          "Cosmic",        "p");
	leg_stack->AddEntry(h_numu_cc,         "Numu CC",       "p");
	leg_stack->AddEntry(h_numu_cc_mixed,   "Numu CC Mixed", "p");
	leg_stack->AddEntry(h_nc,              "NC",            "p");
	leg_stack->AddEntry(h_nc_pi0,          "NC Pi0",        "p");
	leg_stack->AddEntry(h_other_mixed,     "Other Mixed",   "p");
	leg_stack->AddEntry(h_unmatched,       "Unmatched",     "p");
	leg_stack->Draw();
	c1->Print(print_name);
}
