// Script to get the 2d plot of dedx and shower distance and plot it




void plot_2d_dedx_shrdist(){

    gStyle->SetOptStat(0);

    TFile *filein = new TFile("files/mc_truth_out.root", "READ");

    TH2D* h_sig = (TH2D*)filein->Get("h_reco_shr_dEdx_shr_dist_sig_good");
    h_sig->SetFillColorAlpha(kGreen+2,0.8);
    h_sig->Scale(0.1301);
    
    TH2D* h_bkg = (TH2D*)filein->Get("h_reco_shr_dEdx_shr_dist_bkg_good");
    h_bkg->Scale(0.1301);
    
    TH2D* h_ext = (TH2D*)filein->Get("h_reco_shr_dEdx_shr_dist_ext_good");
    h_ext->Scale(1.0154);
    
    
    h_bkg->Add(h_ext);
    
    TCanvas *c = new TCanvas();

    // h_bkg->Draw("box");
    h_sig->Draw("box,same");
    h_bkg ->SetFillColorAlpha(kRed+1, 0.4);
    h_bkg->Draw("box,same");

}


