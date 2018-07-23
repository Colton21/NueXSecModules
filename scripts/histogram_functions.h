#ifndef HISTOGRAM_FUNCTIONS_h
#define HISTOGRAM_FUNCTIONS_h

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMarker.h"
#include "TLine.h"

#include <iostream>

class histogram_functions {

public:

histogram_functions()=default;

static void Plot1DHistogram (TH1 * histogram, const char * x_axis_name, const char * print_name);
static void Plot1DHistogramGausFit (TH1 * histogram, const char * x_axis_name, const char * print_name);
static void PlotTEfficiency (TH1 *h_num, TH1 *h_den, const char * title, const char * print_name);
static void PlotTEfficiency (TH1 *h_num, TH1 *h_den, const bool rebin, const char * title, const char * print_name);
static void Plot2DHistogram (TH2 * histogram, const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name);
static void Plot2DHistogram (TH2 * histogram, const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name,
                             const char * draw_option);
static void Plot2DHistogram (TH2 * histogram, const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name,
                             const char * draw_option, const int x_divisions, const int y_divisions);
static void Plot2DHistogramOffSet (TH2 * histogram, const double label_offset, const double title_offset, const char * title,
                                   const char * x_axis_name, const char * y_axis_name, const char * print_name);
static void PlotSimpleStackParticle(TH1 * h_electron, TH1 * h_proton, TH1 * h_photon, TH1 * h_pion,
                                    TH1 * h_kaon, TH1 * h_muon, TH1 * h_neutron, TH1 * h_unmatched, TH1 * h_ext_unmatched,
                                    const double intime_scale_factor, const double data_scale_factor,
                                    const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                    const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name);
static void PlotSimpleStack (TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                             TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, const char * title,
                             const char * x_axis_name, const char * y_axis_name, const char * print_name);
static void PlotSimpleStackInTime (TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed, TH1 * h_cosmic,
                                   TH1 * h_nc, TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                   const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name);
static void PlotSimpleStackData (TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                 TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                 TH1 * h_data, const double data_scale_factor,
                                 const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name);
static void TimingHistograms(TH1 * histogram_1, TH1 * histogram_2, TH1 * histogram_3,
                             const double data_scale_factor, const double intime_scale_factor,
                             const char * x_axis_name, const char* print_name);
static void TimingHistogramsOverlay(std::vector<std::pair<double, int> > * data_flash_time, TH1 * histogram_1, TH1 * histogram_2,
                                    const double intime_scale_factor, const char * x_axis_name,
                                    const char * print_name1, const char * print_name2);
static void LegoStackData(TH2 * h_nue_cc, TH2 * h_nue_cc_mixed, TH2 * h_numu_cc, TH2 * h_numu_cc_mixed, TH2 * h_cosmic, TH2 * h_nc,
                          TH2 * h_nc_pi0, TH2 * h_other_mixed, TH2 * h_unmatched, TH2 * h_intime, const double intime_scale_factor,
                          TH2 * h_data, const double data_scale_factor,
                          const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                          const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name);
static void PlotSimpleStack(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                            TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched,
                            const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                            const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name);
static void PlotSimpleStackInTime(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                  TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                  const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                  const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name);
static void PlotSimpleStackData(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_nue_cc_out_fv, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                                TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_intime, const double intime_scale_factor,
                                TH1 * h_data, const double data_scale_factor,
                                const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                                const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name);

static void PlotDetailStack(TH1 * h_nue_cc_qe,
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
                            const char * x_axis_name, const char * y_axis_name, const char * print_name);

static void PostHistogramOverlay(TH1 * h_no_cut, TH1 * h_reco_nue,
                                 TH1 * h_in_fv, TH1 * h_vtx_flash, TH1 * h_shwr_vtx,
                                 TH1 * h_trk_vtx, TH1 * h_hit_threshold,
                                 TH1 * h_open_angle, TH1 * h_dedx,
                                 const char * x_axis_name, const char * y_axis_name, const char * print_name);

static void PlotDataMC(TH1 * h_nue_cc, TH1 * h_nue_cc_mixed, TH1 * h_numu_cc, TH1 * h_numu_cc_mixed, TH1 * h_cosmic, TH1 * h_nc,
                       TH1 * h_nc_pi0, TH1 * h_other_mixed, TH1 * h_unmatched, TH1 * h_data, const double data_mc_scale_factor,
                       const char * x_axis_name, const char * y_axis_name, const char * print_name, const char * data_print_name);

static void PlotDetailDataMCStack(TH1 * h_nue_cc_qe,
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
                                  const char * x_axis_name, const char * y_axis_name, const char * print_name, const char * data_print_name);

static void PurityStack(TH1 * h_qe, TH1 * h_res, TH1 * h_dis, TH1 * h_coh, TH1 * h_mec, const char * x_axis_name, const char * print_name);

static void OverlayScatter(TH2 * h_nue_cc, TH2 * h_nue_cc_mixed, TH2 * h_nue_cc_out_fv, TH2 * h_numu_cc, TH2 * h_numu_cc_mixed, TH2 * h_cosmic, TH2 * h_nc,
                           TH2 * h_nc_pi0, TH2 * h_other_mixed, TH2 * h_unmatched,
                           const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                           const char * title, const char * x_axis_name, const char * y_axis_name, const char * print_name);

static void EnergyOverlay(
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
        );
static void PlotdEdxTheta(
        TH2 * h_nue_cc, TH2 * h_nue_cc_mixed, TH2 * h_nue_cc_out_fv, TH2 * h_numu_cc,
        TH2 * h_numu_cc_mixed, TH2 * h_cosmic, TH1 * h_nc,
        TH2 * h_nc_pi0, TH2 * h_other_mixed, TH2 * h_unmatched, TH2 * h_intime, const double intime_scale_factor,
        TH2 * h_data, const double data_scale_factor,
        const char * title, const char * x_axis_name, const char * y_axis_name,
        const char * print_name1, const char * print_name2, const char * print_name3
        );
};
#endif
