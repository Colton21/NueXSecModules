#include "variation_output_bkg.h"

//***************************************************************************
//***************************************************************************

/*
 * Get File Name from a Path with or without extension
 */
std::string getFileName(std::string filePath, bool withExtension = true, char seperator = '/') {
    // Get last dot position
    std::size_t dotPos = filePath.rfind('.');
    std::size_t sepPos = filePath.rfind(seperator);
 
    if(sepPos != std::string::npos)
    {
        return filePath.substr(sepPos + 1, filePath.size() - (withExtension || dotPos != std::string::npos ? 1 : dotPos) );
    }
    return "";
}
//***************************************************************************
//***************************************************************************

// ----------------------
//         Main
// ----------------------
void variation_output_bkg::run_var(const char * _file1, TString mode, const std::vector<double> _config, TString plot_config) {
    // std::cout << "=================================================\n" << std::endl;
    // std::cout << "Warning, there are hardcoded values" << 
    //     "in this script, grep for  \"HARDCODED\" for places\n" << std::endl;
    // std::cout << "=================================================\n" << std::endl;
    
    // create plots folder if it does not exist
    system("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi");

    gStyle->SetOptStat(0); // say no to stats box

    // Choose which file 
    if (mode == "bnb"){
        f_var_out = new TFile("plots/variation_out_bnb_bkg.root","UPDATE");
        draw_mode = "bnb";
    }
    else if (mode == "numi"){
        f_var_out = new TFile("plots/variation_out_numi_bkg.root","UPDATE");
    }
    else return;

    //*************************** Configure the cut parameters *************************************
    std::cout << "\n --- Configuring Parameters --- \n" << std::endl;
    _x1 = _config[0];
    _x2 = _config[1];
    _y1 = _config[2];
    _y2 = _config[3];
    _z1 = _config[4];
    _z2 = _config[5];
    flash_pe_threshold = _config[6];
    
    if (mode == "bnb"){
        std::cout << "Using BNB Params" << std::endl;
        // flash_time_start = 3.19; // Manually Override for BNB [us]
        // flash_time_end   = 4.87; // Manually Override for BNB [us]
        flash_time_start = 3.125; // Manually Override for BNB [us]
        flash_time_end   = 4.725; // Manually Override for BNB [us]

    }
    else if (mode == "numi") {
        std::cout << "Using NuMI Params" << std::endl;
        flash_time_start = _config[7]; // Use default numi config
        flash_time_end   = _config[8]; // Use default numi config
    }
    else return;
    std::cout << "flash_time_start:\t" << flash_time_start << "   flash_time_end:\t" << flash_time_end << std::endl;
    
    tolerance                     = _config[9];
    shwr_nue_tolerance            = _config[10];
    trk_nue_tolerance             = _config[11];
    shwr_hit_threshold            = _config[12];
    shwr_hit_threshold_collection = _config[13];
    tolerance_open_angle_min      = _config[14];
    tolerance_open_angle_max      = _config[15];
    tolerance_dedx_min            = _config[16];
    tolerance_dedx_max            = _config[17];
    dist_tolerance                = _config[18];
    pfp_hits_length_tolerance     = _config[19];
    ratio_tolerance               = _config[20];
    detector_variations           = _config[21];
    const std::vector<double> tolerance_open_angle {tolerance_open_angle_min, tolerance_open_angle_max};

    // Get the directory for the file    
    std::string dirname;
    // Get File name with extension from file path
    dirname = getFileName(_file1, false);
    
    dirname.erase(0, 7); // - "filter_"
    dirname.erase(dirname.end()-5, dirname.end()); // - ".root"

    std::cout << "dirname:\t" << dirname << std::endl;

    //*************************** RUN PLOTTING FUNCTION PLOT***********************
    // if bool true just run this function
    if (plot_config == "same") {
        PlotVariatons(f_var_out, mode); 
        return; 
    }
    if (plot_config == "bnbnumi"){
        PlotVariatonsNuMIBNB(); 
        return; 
    }
    if (plot_config == "makeweights"){
        GenerateWeightHistograms(); // Generate the weight histograms if the BNB CV is run
        return;
    }
    if(plot_config == "weightednumiplots"){
        CompareWeightedHistograms();
        return;
    }
    //*************************** POT Scaling *************************************
    std::cout << "=================================================\n" << std::endl;
    
    double CV_POT;
    
    if (mode == "bnb") 
        CV_POT =  GetPOT("/uboone/data/users/kmistry/work/NueXSection_Outputs/detector_variations/bnb_det_var_cv/filter_BNBCV.root");
    else
        CV_POT = GetPOT("/uboone/data/users/kmistry/work/NueXSection_Outputs/detector_variations/numi_det_var/filter_NuMICV.root");
        
    double POT_Scaling;
    // if (mode == "numi") POT_Scaling =  1.0;
    
    POT_Scaling =  CV_POT / GetPOT(_file1);
    
    std::cout << "POT Scaling:\t" << POT_Scaling << std::endl;
    std::cout << "=================================================\n" << std::endl;
    
    // Get the weighted histogram file
    fweight = TFile::Open("plots/variation_weights.root");

    //*********** Open Variation File and get branches ******************************
    // Check if the outfile opened successfully
    if ( f_var_out->IsOpen() ) std::cout << "Variation File opened successfully\n" << std::endl;
    
    // Get Variation file
    TFile* inFile = new TFile(_file1);
    
    // Get TPC Obj Information from file
    TTree* TPCObjTree = (TTree*) inFile->Get("AnalyzeTPCO/tree");    
    std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
    TPCObjTree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

    // Num events in the tree
    const int tree_total_entries = TPCObjTree->GetEntries();
    std::cout << "Total Events: " << tree_total_entries << std::endl;

    // Get the largest Flash Vector of Vector    
    std::vector<std::vector<double>> largest_flash_v_v = GetLargestFlashVector(inFile, flash_time_start, flash_time_end);

    //************* Get list of flashes that pass flash in time and flash pe cut ************************
    std::vector<bool> flash_cuts_pass_vec;
    FlashinTime_FlashPE(inFile, flash_time_start, flash_time_end, flash_cuts_pass_vec, mode );
    //***************************************************************************************************

    // MCTruth Counting Tree
    TTree * mctruth_counter_tree = (TTree*)inFile->Get("AnalyzeTPCO/mctruth_counter_tree"); 
    mctruth_counter_tree->SetBranchAddress("fMCNuEnergy", &mc_nu_energy);
    mctruth_counter_tree->SetBranchAddress("fMCNuID", &mc_nu_id);
    mctruth_counter_tree->SetBranchAddress("fMCNuVtxX", &mc_nu_vtx_x);
    mctruth_counter_tree->SetBranchAddress("fMCNuVtxY", &mc_nu_vtx_y);
    mctruth_counter_tree->SetBranchAddress("fMCNuVtxZ", &mc_nu_vtx_z);
    mctruth_counter_tree->SetBranchAddress("fMCNuDirX", &mc_nu_dir_x);
    mctruth_counter_tree->SetBranchAddress("fMCNuDirY", &mc_nu_dir_y);
    mctruth_counter_tree->SetBranchAddress("fMCNuDirZ", &mc_nu_dir_z);
    mctruth_counter_tree->SetBranchAddress("fMCNumParticles", &mc_nu_num_particles);
    mctruth_counter_tree->SetBranchAddress("fMCNumChargedParticles", &mc_nu_num_charged_particles);
    mctruth_counter_tree->SetBranchAddress("fMCEleDirX", &mc_ele_dir_x);
    mctruth_counter_tree->SetBranchAddress("fMCEleDirY", &mc_ele_dir_y);
    mctruth_counter_tree->SetBranchAddress("fMCEleDirZ", &mc_ele_dir_z);
    mctruth_counter_tree->SetBranchAddress("fMCEleEnergy", &mc_ele_energy);
    mctruth_counter_tree->SetBranchAddress("fMCEleMomentum", &mc_ele_momentum);
    mctruth_counter_tree->SetBranchAddress("has_pi0", &has_pi0);
    mctruth_counter_tree->SetBranchAddress("fMCNuTime", &mc_nu_time);

    // Make branches for the selected tree
    selected_tree->Branch("mc_Phi",    &mc_Phi,   "mc_Phi/D");
    selected_tree->Branch("mc_Theta",  &mc_Theta, "mc_Theta/D");
    selected_tree->Branch("mc_Energy", &mc_Energy,"mc_Energy/D");
    selected_tree->Branch("bkg_class", &bkg_class);


    // Define the FV
    std::vector<double> fv_boundary_v = {_x1, _x2, _y1, _y2, _z1, _z2};

    std::cout << "Total Events (MC): " << mctruth_counter_tree->GetEntries() << std::endl;

    // ----------------------
    //      Event loop
    // ----------------------
    std::cout << "Starting Eventloop..." << std::endl;
    for(int event = 0; event < tree_total_entries; event++){

        if (event % 100000 == 0) std::cout << "On entry " << event/100000.0 <<"00k" << std::endl;
        
        TPCObjTree->GetEntry(event);
        mctruth_counter_tree->GetEntry(event);

        int n_tpc_obj = tpc_object_container_v->size();

        // --------------- MC Counters ---------------
        if(mc_nu_id == 1) mc_nue_cc_counter++;
        if(mc_nu_id == 3) mc_nue_nc_counter++;
        if(mc_nu_id == 5) mc_nue_cc_counter_bar++;
        if(mc_nu_id == 7) mc_nue_nc_counter_bar++;
        
        // Filter for truth Nues
        // if (mc_nu_id == 2 || mc_nu_id == 4 || mc_nu_id == 6 || mc_nu_id == 8) continue;

        // --------------- Flash Information ---------------
        std::vector<double> largest_flash_v     = largest_flash_v_v.at(event); // Vec with the largest flash
        largest_flash_y     = largest_flash_v.at(0);
        largest_flash_z     = largest_flash_v.at(1);
        largest_flash_time     = largest_flash_v.at(2);
        largest_flash_pe     = largest_flash_v.at(3);

        // Only fill for matched events (where each entry is non-zero)
        if (largest_flash_y != 0 && largest_flash_z !=0 && largest_flash_time != 0 && largest_flash_pe != 0){
            TH1D_hist.at(klargest_flash_y)      ->Fill(largest_flash_y);
            TH1D_hist.at(klargest_flash_z)      ->Fill(largest_flash_z);
            TH1D_hist.at(klargest_flash_time) ->Fill(largest_flash_time);
            TH1D_hist.at(klargest_flash_pe)      ->Fill(largest_flash_pe);
        }
        
        // -------------------------------------------------

        // Loop over TPCObj
        for (int i = 0; i < n_tpc_obj; i++){
            auto const tpc_obj = tpc_object_container_v->at(i);

            // Check to see in truth the vertex was inside the tpc FV
            const bool true_in_tpc = in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v);

            // Classify the event  
            std::pair<std::string, int> tpc_classification = TPCO_Classifier(tpc_obj, true_in_tpc, has_pi0);

            // Checks the classification for Signal
            bool bool_sig{false}; 
            for (unsigned int k=0; k < signal_modes.size(); k++){
                if (signal_modes.at(k).find(tpc_classification.first) != std::string::npos) {
                    bool_sig = true;
                }
            }

            // TPC Obj vars
            tpc_obj_vtx_x = tpc_obj.pfpVtxX();
            tpc_obj_vtx_y = tpc_obj.pfpVtxY();
            tpc_obj_vtx_z = tpc_obj.pfpVtxZ();
           
            const int tpc_obj_mode = tpc_obj.Mode(); 
            n_pfp = tpc_obj.NumPFParticles();
            const int leading_shower_index = GetLeadingShowerIndex(n_pfp, n_tpc_obj, tpc_obj);

            // Other histograms
            n_tracks = 0; n_showers = 0; n_pfp_50Hits = 0; n_tracks_50Hits = 0; n_showers_50Hits = 0;
            GetNumber_Track_Shower(n_pfp, n_tpc_obj, tpc_obj, n_showers, n_tracks, n_pfp_50Hits, n_tracks_50Hits, n_showers_50Hits );// Get the number of tracks and showers
            TH1D_hist.at(kn_pfp)            ->Fill(n_tracks + n_showers);
            TH1D_hist.at(kn_pfp_50Hits)     ->Fill(n_pfp_50Hits);
            TH1D_hist.at(kn_tracks)         ->Fill(n_tracks);
            TH1D_hist.at(kn_tracks_50Hits)  ->Fill(n_tracks_50Hits);
            TH1D_hist.at(kn_showers)        ->Fill(n_showers);
            TH1D_hist.at(kn_showers_50Hits) ->Fill(n_showers_50Hits);

            const double Flash_TPCObj_Dist = Flash_TPCObj_vtx_Dist(tpc_obj_vtx_y, tpc_obj_vtx_z, largest_flash_y, largest_flash_z);
            TH1D_hist.at(kFlash_TPCObj_Dist) ->Fill(Flash_TPCObj_Dist);

            if (bool_sig) nue_cc_counter++;

            // Fill the unselected histograms
            if (!bool_sig && (dirname == "BNBCV" || dirname == "NuMICV")) {

                    // Loop over the Par Objects
                for (int j = 0; j < n_pfp ; j++){

                    auto const pfp_obj = tpc_obj.GetParticle(j);
      
                    const double mc_Theta  = pfp_obj.mcTheta();
                    mc_Phi    = pfp_obj.mcPhi();
                    mc_Energy = pfp_obj.mcEnergy();
                    const double mc_pdg    = pfp_obj.MCPdgCode();
                    const double mc_pdg_parent = pfp_obj.MCParentPdg();
                    
                    // Background events
                    if (!bool_sig) {

                        //  ------------ Leading shower ------------
                        if (j == leading_shower_index){

                            bkg_class = Background_Classifier(mc_pdg, tpc_classification.first);

                            TH1D_hist.at(kldg_shwr_Phi_unselected)             ->Fill(mc_Phi);
                            TH1D_hist.at(kldg_shwr_Theta_unselected)         ->Fill(mc_Theta);
                           
                            // Wrap phi from 0 to 90
                            double leading_shower_phi_wrapped   = WrapPhi(mc_Phi);
                            double mc_phi_wrapped               = WrapPhi(mc_Phi);
                            TH1D_hist.at(kldg_shwr_Phi_wrapped_unselected) ->Fill(mc_phi_wrapped);
                            TH1D_hist.at(kshower_E_unselected)             ->Fill(mc_Energy);

                            // Fill the phi distribtuons for the bkgs
                            if (bkg_class == "pi0_gamma") {
                                mc_Phi_pi0 = mc_Phi;
                                TH1D_hist.at(kshower_phi_pi0_unselected)          ->Fill(mc_Phi);
                                TH1D_hist.at(kshower_E_pi0_unselected)            ->Fill(mc_Energy);
                                TH1D_hist.at(kshower_phi_pi0_wrapped_unselected)  ->Fill(mc_phi_wrapped);
                                TH1D_hist.at(kshower_Theta_pi0_unselected)        ->Fill(mc_Theta);
                                TH2D_hist.at(kEBkg_pi0_Theta_unselected)          ->Fill(mc_Energy, mc_Theta);
                                TH2D_hist.at(kEBkg_pi0_Phi_unselected)            ->Fill(mc_Energy, mc_Phi);
                                TH2D_hist.at(kEBkg_pi0_Phi_wrapped_unselected)    ->Fill(mc_Energy, mc_phi_wrapped);
                                TH2D_hist.at(kThetaBkg_pi0_Phi_wrapped_unselected)->Fill(mc_Theta, mc_phi_wrapped);
                                TH2D_hist.at(kThetaBkg_pi0_Phi_unselected)        ->Fill(mc_Theta, mc_Phi);
                            }
                            
                            if (bkg_class == "cosmic") {
                                TH1D_hist.at(kshower_phi_bkg_cosmic_unselected)         ->Fill(mc_Phi);
                                TH1D_hist.at(kshower_E_bkg_cosmic_unselected)           ->Fill(mc_Energy);
                                TH1D_hist.at(kshower_phi_bkg_cosmic_wrapped_unselected) ->Fill(mc_phi_wrapped);
                                TH1D_hist.at(kshower_Theta_bkg_cosmic_unselected)        ->Fill(mc_Theta);
                                TH2D_hist.at(kEBkg_cosmic_Theta_unselected)             ->Fill(mc_Energy, mc_Theta);
                                TH2D_hist.at(kEBkg_cosmic_Phi_unselected)               ->Fill(mc_Energy, mc_Phi);
                                TH2D_hist.at(kEBkg_cosmic_Phi_wrapped_unselected)       ->Fill(mc_Energy, mc_phi_wrapped);
                                TH2D_hist.at(kThetaBkg_cosmic_Phi_wrapped_unselected)   ->Fill(mc_Theta, mc_phi_wrapped);
                                TH2D_hist.at(kThetaBkg_cosmic_Phi_unselected)           ->Fill(mc_Theta, mc_Phi);
                            }
                            
                            if (bkg_class == "other_bkg") {
                                mc_Phi_other = mc_Phi;
                                TH1D_hist.at(kshower_phi_other_unselected)          ->Fill(mc_Phi);
                                TH1D_hist.at(kshower_E_other_unselected)            ->Fill(mc_Energy);
                                TH1D_hist.at(kshower_phi_other_wrapped_unselected)  ->Fill(mc_phi_wrapped);
                                TH1D_hist.at(kshower_Theta_other_unselected)        ->Fill(mc_Theta);
                                TH2D_hist.at(kEBkg_other_Theta_unselected)          ->Fill(mc_Energy, mc_Theta);
                                TH2D_hist.at(kEBkg_other_Phi_unselected)            ->Fill(mc_Energy, mc_Phi);
                                TH2D_hist.at(kEBkg_other_Phi_wrapped_unselected)    ->Fill(mc_Energy, mc_phi_wrapped);
                                TH2D_hist.at(kThetaBkg_other_Phi_wrapped_unselected)->Fill(mc_Theta, mc_phi_wrapped);
                                TH2D_hist.at(kThetaBkg_other_Phi_unselected)        ->Fill(mc_Theta, mc_Phi);
                            }
                        
                            // 2D histos of the energy of the background particle vs angle
                            TH2D_hist.at(kEBkg_Theta_unselected)         ->Fill(mc_Energy, mc_Theta);
                            TH2D_hist.at(kEBkg_Phi_unselected)           ->Fill(mc_Energy, mc_Phi);
                            TH2D_hist.at(kEBkg_Phi_wrapped_unselected)   ->Fill(mc_Energy, mc_phi_wrapped);
                            
                            TH2D_hist.at(kThetaBkg_Phi_wrapped_unselected)   ->Fill(mc_Theta, mc_phi_wrapped);
                            TH2D_hist.at(kThetaBkg_Phi_unselected)           ->Fill(mc_Theta, mc_Phi);


                        }
                    
                    }
                

                } // END LOOP PAR OBJ

            }

            //************************ Apply Flash in time and Flash PE cut *************************************
            if (flash_cuts_pass_vec.at(event) == false) continue;
            if (bool_sig) counter_FlashinTime_FlashPE++;
            //***************************************************************************************************
            
            //****************************** Apply Pandora Reco Nue cut *****************************************
            bool bool_HasNue = HasNue(tpc_obj, n_pfp );
            if ( bool_HasNue == false ) continue;
            if (bool_sig) counter_HasNue++;
            //***************************************************************************************************
            
            //****************************** Apply In FV cut ****************************************************
            bool bool_inFV = in_fv(tpc_obj_vtx_x, tpc_obj_vtx_y, tpc_obj_vtx_z, fv_boundary_v);
            if ( bool_inFV == false ) continue;
            if (bool_sig) counter_inFV++;
            //***************************************************************************************************
            
            //****************************** Apply flash vtx cut ************************************************
            bool bool_flashvtx = flashRecoVtxDist(largest_flash_v, tolerance, tpc_obj_vtx_x, tpc_obj_vtx_y,  tpc_obj_vtx_z);
            if ( bool_flashvtx == false ) continue;
            if (bool_sig) counter_FlashRecoVtxDist++;
            //***************************************************************************************************     

            //****************************** Apply vtx nu distance cut ******************************************
            bool bool_vtxnudist = VtxNuDistance( tpc_obj, 11, shwr_nue_tolerance);
            if ( bool_vtxnudist == false ) continue;
            if (bool_sig) counter_VtxNuDist++;
            //***************************************************************************************************
        
            //****************************** Apply track vtx nu distance cut ************************************
            bool bool_trackvtxnudist = VtxNuDistance( tpc_obj, 13, trk_nue_tolerance);
            if ( bool_trackvtxnudist == false ) continue;
            if (bool_sig) counter_VtxTrackNuDist++;
            //***************************************************************************************************
            std::string hist_fill_name = "kpre_hit_threshold_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            //****************************** Apply Hit threshold cut**** cut ************************************
            bool bool_hitTh = HitThreshold(tpc_obj, shwr_hit_threshold, false);
            if ( bool_hitTh == false ) continue;
            if (bool_sig) counter_HitThresh++;
            //***************************************************************************************************
            hist_fill_name = "kpost_hit_threshold_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            hist_fill_name = "kpre_hit_threshold_collection_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            //****************************** Apply Hit threshold collection cut *********************************
            bool bool_hitThW = HitThreshold(tpc_obj, shwr_hit_threshold_collection, true);
            if ( bool_hitThW == false ) continue;
            if (bool_sig) counter_HitThreshW++;
            //***************************************************************************************************
            hist_fill_name = "kpost_hit_threshold_collection_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            hist_fill_name = "kpre_open_angle_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            //****************************** Apply Open Angle cut ***********************************************
            bool bool_OpenAngle = OpenAngleCut(tpc_obj, tolerance_open_angle);
            if ( bool_OpenAngle == false ) continue;
            if (bool_sig) counter_OpenAngle++;
            //***************************************************************************************************
            hist_fill_name = "kpost_open_angle_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            hist_fill_name = "kpre_dedx_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            //****************************** Apply dEdx cut *****************************************************
            bool bool_dEdx = dEdxCut(tpc_obj, tolerance_dedx_min, tolerance_dedx_max);
            if ( bool_dEdx == false ) continue;
            if (bool_sig) counter_dEdx++;
            //***************************************************************************************************
            hist_fill_name = "kpost_dedx_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            //************************* Apply Secondary shower dist cut *****************************************
            bool bool_sshwrcut = SecondaryShowersDistCut(tpc_obj, dist_tolerance);
            if ( bool_sshwrcut == false ) continue;
            if (bool_sig) counter_SecondaryShowers++;
            //***************************************************************************************************

            //************************* Apply hit per lengh ratio cut *******************************************
            bool bool_hitlengthratio = HitLengthRatioCut( pfp_hits_length_tolerance, tpc_obj);
            if ( bool_hitlengthratio == false ) continue;
            if (bool_sig) counter_HitLenghtRatio++;
            //***************************************************************************************************

            //************************* Apply Longest Track Leading Shower cut **********************************
            bool bool_longtrackleadingshwr = LongestTrackLeadingShowerCut(ratio_tolerance, tpc_obj);
            if ( bool_longtrackleadingshwr == false ) continue;
            if (bool_sig) counter_LongestTrackLeadingShower++;
            //***************************************************************************************************
            
            //************************* Apply Contained Track Cut ***********************************************
            bool bool_contTrack = ContainedTracksCut(fv_boundary_v, tpc_obj);
            if ( bool_contTrack  == false ) continue;
            if (bool_sig) counter_ContainedTrack++;
            //***************************************************************************************************
            
            //************************* Skip unmatched and other mixed events ***********************************
            if (tpc_classification.first == "unmatched" || tpc_classification.first == "other_mixed") continue;
            //***************************************************************************************************

            if (bool_sig) sig_counter++;
            if (!bool_sig) bkg_counter++;

            //if (!bool_sig ) std::cout << tpc_classification.first << std::endl;

            // Weight the BNB Variations
            WeightBNBVar(tpc_obj, bool_sig, leading_shower_index, tpc_classification, dirname);

            // Loop over the Par Objects
            for (int j = 0; j < n_pfp ; j++){

                auto const pfp_obj = tpc_obj.GetParticle(j);

                // PFP vars                
                const std::string mc_origin = pfp_obj.Origin();
                const int  pfp_pdg          = pfp_obj.PFParticlePdgCode();
                const int  num_pfp_hits     = pfp_obj.NumPFPHits();
                const int  mc_parent_pdg    = pfp_obj.MCParentPdg();
                const int  n_pfp_hits_w     = pfp_obj.NumPFPHitsW(); // Collection plane hits
                
                const double pfp_length     = pfp_obj.pfpLength();
                const double pfp_open_angle = pfp_obj.pfpOpenAngle();
                const double leading_dedx   = pfp_obj.PfpdEdx().at(2);//just the collection plane!

                const double pfp_vtx_x =  pfp_obj.pfpVtxX();
                const double pfp_vtx_y =  pfp_obj.pfpVtxY();
                const double pfp_vtx_z =  pfp_obj.pfpVtxZ();

                const double pfp_dir_x = pfp_obj.pfpDirX();
                const double pfp_dir_y = pfp_obj.pfpDirY();
                const double pfp_dir_z = pfp_obj.pfpDirZ();

                mc_Theta  = pfp_obj.mcTheta();
                mc_Phi    = pfp_obj.mcPhi();
                mc_Energy = pfp_obj.mcEnergy();
                const double mc_pdg    = pfp_obj.MCPdgCode();
                const double mc_pdg_parent = pfp_obj.MCParentPdg();

                const double pfp_end_x = (pfp_obj.pfpVtxX() + (pfp_length * pfp_dir_x));
                const double pfp_end_y = (pfp_obj.pfpVtxY() + (pfp_length * pfp_dir_y));
                const double pfp_end_z = (pfp_obj.pfpVtxZ() + (pfp_length * pfp_dir_z));

                std::vector<double> pfp_start_vtx {pfp_vtx_x, pfp_vtx_y, pfp_vtx_z};
                std::vector<double> pfp_end_vtx {pfp_end_x, pfp_end_y, pfp_end_z};

                const double pfp_Nu_vtx_Dist =  pfp_vtx_distance(tpc_obj_vtx_x, tpc_obj_vtx_y, tpc_obj_vtx_z, pfp_vtx_x, pfp_vtx_y, pfp_vtx_z);
                
                // Background events
                if (!bool_sig) {

                    // std::cout << tpc_classification.first << std::endl;
                    
                    TH1D_hist.at(ktotal_hits) ->Fill(num_pfp_hits);

                    const double shower_phi = atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
                    TH1D_hist.at(kshower_phi) ->Fill(shower_phi);

                    TH1D_hist.at(kshower_Nu_vtx_Dist) ->Fill(pfp_Nu_vtx_Dist);

                    //  ------------ Leading shower ------------
                    if (j == leading_shower_index){

                        bkg_class = Background_Classifier(mc_pdg, tpc_classification.first);

                        // std::cout << "mc_pdg: " << mc_pdg << " parent pdg:  " << mc_parent_pdg <<"  Classifier: " << tpc_classification.first << "  bkg class:  " << bkg_class << std::endl;

                        TH1D_hist.at(kselected) ->Fill(0);

                        double leading_shower_phi     = atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
                        const double leading_shower_theta     = acos(pfp_obj.pfpDirZ()) * 180 / 3.1415;
                        
                        const double leading_shower_length     = pfp_obj.pfpLength();
                        const double longest_track_length      = GetLongestTrackLength(n_pfp, n_tpc_obj, tpc_obj);

                        TH1D_hist.at(kldg_shwr_hits)          ->Fill(num_pfp_hits);
                        TH1D_hist.at(kldg_shwr_hits_WPlane)   ->Fill(pfp_obj.NumPFPHitsW());         // W Plane
                        TH1D_hist.at(kldg_shwr_Open_Angle)    ->Fill(pfp_obj.pfpOpenAngle() * (180 / 3.1415) );
                        TH1D_hist.at(kldg_shwr_dEdx_WPlane)   ->Fill(pfp_obj.PfpdEdx().at(2) );     // W Plane
                        TH1D_hist.at(kldg_shwr_HitPerLen)     ->Fill(num_pfp_hits / pfp_obj.pfpLength() );
                        TH1D_hist.at(kldg_shwr_Phi)           ->Fill(mc_Phi);
                        TH1D_hist.at(kldg_shwr_Theta)         ->Fill(mc_Theta);
                        TH1D_hist.at(kldg_shwr_CTheta)        ->Fill(cos(leading_shower_theta * 3.1414 / 180.));
                        TH1D_hist.at(klong_Track_ldg_shwr)    ->Fill(longest_track_length / leading_shower_length);

                        // Wrap phi from 0 to 90
                        double leading_shower_phi_wrapped   = WrapPhi(mc_Phi);
                        double mc_phi_wrapped               = WrapPhi(mc_Phi);
                        TH1D_hist.at(kldg_shwr_Phi_wrapped) ->Fill(mc_phi_wrapped);
                        TH1D_hist.at(kshower_E)             ->Fill(mc_Energy);

                        // Weight the backgrounds
                        GetBNBBkgWeight(mc_Theta, mc_Phi, mc_phi_wrapped, bkg_class, weight_all, weight_indiv, dirname );
                        // std::cout << "weight all:  " << weight_all << "  weight_indiv:   " << weight_indiv << "  bkg CV: " << bkg_counter << std::endl;

                        // Fill the phi distribtuons for the bkgs
                        if (bkg_class == "pi0_gamma") {
                            mc_Phi_pi0 = mc_Phi;
                            TH1D_hist.at(kshower_phi_pi0)          ->Fill(mc_Phi);
                            TH1D_hist.at(kshower_E_pi0)            ->Fill(mc_Energy);
                            TH1D_hist.at(kshower_phi_pi0_wrapped)  ->Fill(mc_phi_wrapped);
                            TH1D_hist.at(kshower_Theta_pi0)        ->Fill(mc_Theta);
                            TH2D_hist.at(kEBkg_pi0_Theta)          ->Fill(mc_Energy, mc_Theta);
                            TH2D_hist.at(kEBkg_pi0_Phi)            ->Fill(mc_Energy, mc_Phi);
                            TH2D_hist.at(kEBkg_pi0_Phi_wrapped)    ->Fill(mc_Energy, mc_phi_wrapped);
                            TH2D_hist.at(kThetaBkg_pi0_Phi_wrapped)->Fill(mc_Theta, mc_phi_wrapped);
                            TH2D_hist.at(kThetaBkg_pi0_Phi)        ->Fill(mc_Theta, mc_Phi);
                            pi0_counter++;
                        }
                        
                        if (bkg_class == "cosmic") {
                            TH1D_hist.at(kshower_phi_bkg_cosmic)         ->Fill(mc_Phi);
                            TH1D_hist.at(kshower_E_bkg_cosmic)           ->Fill(mc_Energy);
                            TH1D_hist.at(kshower_phi_bkg_cosmic_wrapped) ->Fill(mc_phi_wrapped);
                            TH1D_hist.at(kshower_Theta_bkg_cosmic)       ->Fill(mc_Theta);
                            TH2D_hist.at(kEBkg_cosmic_Theta)             ->Fill(mc_Energy, mc_Theta);
                            TH2D_hist.at(kEBkg_cosmic_Phi)               ->Fill(mc_Energy, mc_Phi);
                            TH2D_hist.at(kEBkg_cosmic_Phi_wrapped)       ->Fill(mc_Energy, mc_phi_wrapped);
                            TH2D_hist.at(kThetaBkg_cosmic_Phi_wrapped)   ->Fill(mc_Theta, mc_phi_wrapped);
                            TH2D_hist.at(kThetaBkg_cosmic_Phi)           ->Fill(mc_Theta, mc_Phi);
                            cosmic_counter++;
                        }
                        
                        if (bkg_class == "other_bkg") {
                            mc_Phi_other = mc_Phi;
                            TH1D_hist.at(kshower_phi_other)          ->Fill(mc_Phi);
                            TH1D_hist.at(kshower_E_other)            ->Fill(mc_Energy);
                            TH1D_hist.at(kshower_phi_other_wrapped)  ->Fill(mc_phi_wrapped);
                            TH1D_hist.at(kshower_Theta_other)        ->Fill(mc_Theta);
                            TH2D_hist.at(kEBkg_other_Theta)          ->Fill(mc_Energy, mc_Theta);
                            TH2D_hist.at(kEBkg_other_Phi)            ->Fill(mc_Energy, mc_Phi);
                            TH2D_hist.at(kEBkg_other_Phi_wrapped)    ->Fill(mc_Energy, mc_phi_wrapped);
                            TH2D_hist.at(kThetaBkg_other_Phi_wrapped)->Fill(mc_Theta, mc_phi_wrapped);
                            TH2D_hist.at(kThetaBkg_other_Phi)        ->Fill(mc_Theta, mc_Phi);
                            other_bkg_counter++;
                        }

                        // Vertex Information - require a shower so fill once  when leading shower
                        TH1D_hist.at(ktpc_obj_vtx_x) ->Fill(tpc_obj_vtx_x);
                        TH1D_hist.at(ktpc_obj_vtx_y) ->Fill(tpc_obj_vtx_y);
                        TH1D_hist.at(ktpc_obj_vtx_z) ->Fill(tpc_obj_vtx_z);
                    
                        // 2D histos of the energy of the background particle vs angle
                        TH2D_hist.at(kEBkg_Theta)         ->Fill(mc_Energy, mc_Theta);
                        TH2D_hist.at(kEBkg_Phi)           ->Fill(mc_Energy, mc_Phi);
                        TH2D_hist.at(kEBkg_Phi_wrapped)   ->Fill(mc_Energy, mc_phi_wrapped);
                        
                        TH2D_hist.at(kThetaBkg_Phi_wrapped)   ->Fill(mc_Theta, mc_phi_wrapped);
                        TH2D_hist.at(kThetaBkg_Phi)           ->Fill(mc_Theta, mc_Phi);

                        VariableTree->Fill();
                        selected_tree->Fill();

                    }
                
                }
                
                
                // Track like
                if ( pfp_pdg == 13) {
                    // Background events
                    if (!bool_sig) {
                        const double track_phi     = atan2(pfp_obj.    pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
                        TH1D_hist.at(ktrack_phi) ->Fill(track_phi);
                        TH1D_hist.at(ktrack_Nu_vtx_Dist) ->Fill(pfp_Nu_vtx_Dist);
                    }
                }

            } // END LOOP PAR OBJ

        } // END LOOP TPCO

    } // END EVENT LOOP
    std::cout << "Finished Eventloop..." << std::endl;
        
    std::cout << "--------------- MC Truth COUNTERS -----------------" << std::endl;
    std::cout << "MC Nue CC Counter      --- " << mc_nue_cc_counter << std::endl;
    std::cout << "MC Nue NC Counter      --- " << mc_nue_nc_counter << std::endl;
    std::cout << "MC Nue CC Counter Bar  --- " << mc_nue_cc_counter_bar << std::endl;
    std::cout << "MC Nue NC Counter Bar  --- " << mc_nue_nc_counter_bar << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "--------------- Cut COUNTERS ---------------------" << std::endl;
    std::cout << "Flash in Time, Flash PE      --- " << counter_FlashinTime_FlashPE << std::endl;
    std::cout << "Reco Nue                     --- " << counter_HasNue << std::endl;
    std::cout << "In FV                        --- " << counter_inFV << std::endl;
    std::cout << "Flash Reco Vtx Dist          --- " << counter_FlashRecoVtxDist << std::endl;
    std::cout << "Vtx Nu Dist                  --- " << counter_VtxNuDist << std::endl;
    std::cout << "Vtx Track Nu Dist            --- " << counter_VtxTrackNuDist << std::endl;
    std::cout << "Hit Threshold                --- " << counter_HitThresh << std::endl;
    std::cout << "Hit Threshold Collection     --- " << counter_HitThreshW << std::endl;
    std::cout << "Open Angle                   --- " << counter_OpenAngle << std::endl;
    std::cout << "dEdx                         --- " << counter_dEdx << std::endl;
    std::cout << "Secondary Showers            --- " << counter_SecondaryShowers << std::endl;
    std::cout << "Hit Length Ratio             --- " << counter_HitLenghtRatio << std::endl;
    std::cout << "Longest Track Leading Shower --- " << counter_LongestTrackLeadingShower << std::endl;
    std::cout << "Contained Tracks             --- " << counter_ContainedTrack << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "--------------- Reco COUNTERS ---------------------" << std::endl;
    std::cout << "True Nue CC in FV --- " << nue_cc_counter << std::endl;
    std::cout << "RECO Signal       --- " << sig_counter << std::endl;
    std::cout << "RECO Background   --- " << bkg_counter << std::endl;
    std::cout << "RECO Bkg W All    --- " << weight_all << std::endl;
    std::cout << "RECO Bkg W Indiv  --- " << weight_indiv << std::endl;
    std::cout << "Pi0 Counter       --- " << pi0_counter << std::endl;
    std::cout << "Cosmic Counter    --- " << cosmic_counter << std::endl;
    std::cout << "Other Counter     --- " << other_bkg_counter << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    
    // ----------------------
    //    Save to a file
    // ----------------------
    TDirectory* savedir = gDirectory; //  Create the directory
    TDirectory* subdir;

    f_var_out->cd();

    // If directory does not exist then make it
    savedir = (TDirectory*)f_var_out->Get(dirname.c_str());
    if (savedir == NULL ) {
        savedir = gDirectory;
        std::cout << dirname << " directory does not exist, creating..." << std::endl;
        subdir = savedir->mkdir(dirname.c_str() ) ;
    }
    else {
        std::cout << dirname <<" directory exists, overwriting..." << std::endl;
        subdir = savedir;
    }
    
    subdir->cd();

    // ----------------------
    //    Draw Histograms and write to file
    // ----------------------
    for (int i=0; i < TH1D_hist.size(); i++){
        DrawTH1D(TH1D_hist.at(i), POT_Scaling);
        TH1D_hist.at(i)->Write("",TObject::kOverwrite);
    } 

    for (int i=0; i < TH2D_hist.size(); i++){
        DrawTH2D(TH2D_hist.at(i), POT_Scaling);
        TH2D_hist.at(i)->Write("",TObject::kOverwrite);
    }

    VariableTree->Write("",TObject::kOverwrite);
    selected_tree->Write("",TObject::kOverwrite);
    f_var_out->Close();

    // -------------------------------------------------------------------------

    TFile *fnumi = new TFile("plots/variation_out_numi_bkg.root","UPDATE");
    fnumi->cd();

    // Now we want to overwrite the weighted histograms
    // We only want to save when the specific variaion is being ran
    std::vector<std::string> variations = {
                                    "BNBwithDIC",
                                    "BNBdataSCE",
                                    "BNBdeadSaturatedChannels",
                                    "BNBLArG4BugFix",
                                    "BNBDLup",
                                    "BNBDLdown",
                                    "BNBDTup",
                                    "BNBDTdown",
                                    "BNBnoiseAmpUp",
                                    "BNBnoiseAmpDown",
                                    "BNBaltDeadChannels",
                                    "BNBstretchResp",
                                    "BNBsqueezeResp",
                                    "BNBupPEnoise",
                                    "BNBEnhancedTPCVis",
                                    "BNBBirksRecomb",
                                    "BNBdownPEnoise"};

    int var_index = -1; // Index of variation

    // Get the index of variation to weight
    for (unsigned int j=0; j< variations.size(); j++){
        if (dirname == variations.at(j)) var_index = j;
    }

    // Fill weighted histograms, only do this for matched variations
    if (var_index != -1){
        dirname = "BNBWeighted_Variations";

        // If directory does not exist then make it
        savedir = (TDirectory*)fnumi->Get(dirname.c_str());
        if (savedir == NULL ) {
            savedir = gDirectory;
            std::cout << dirname << " directory does not exist, creating..." << std::endl;
            subdir = savedir->mkdir(dirname.c_str() ) ;
        }
        else {
            std::cout << dirname <<" directory exists, overwriting..." << std::endl;
            subdir = savedir;
        }
        
        subdir->cd();

        // Do the weighted histograms
        for (int i=0; i < TH1D_hist_weighted.size() ; i++){

            TCanvas* c = new TCanvas();
            TLegend* dummy_leg = new TLegend();
            c->cd();
            TH1D_hist_weighted.at(i).at(var_index)->SetLineWidth(2);
            TH1D_hist_weighted.at(i).at(var_index)->SetLineStyle(1);
            TH1D_hist_weighted.at(i).at(var_index)->Scale(POT_Scaling);
            TH1D_hist_weighted.at(i).at(var_index)->SetOption("hist, E");
            DrawTH1D_SAME(TH1D_hist_weighted.at(i).at(var_index), bnbvars.at(i), dummy_leg, "dummy");

            TH1D_hist_weighted.at(i).at(var_index)->Write("",TObject::kOverwrite);

        } 
    }

} // END MAIN
//***************************************************************************
//***************************************************************************
void variation_output_bkg::PlotVariatons(TFile* f_var_out, TString mode){
    f_var_out->cd();

    system("if [ ! -d \"plots/bnb/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/bnb/; fi");
    system("if [ ! -d \"plots/numi/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/numi/; fi");

    system("if [ ! -d \"plots_png/bnb/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots_png/bnb/; fi");
    system("if [ ! -d \"plots_png/numi/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots_png/numi/; fi");

    std::string mode_str;

    if (mode == "bnb")
        mode_str = "bnb";
    else 
        mode_str = "numi";
    
    // Grab the variation folders in the file
    std::vector<std::string> variations = variation_output_bkg::GrabDirs(f_var_out); 

     // ************************** 1D Histograms *******************************
    std::vector<std::string> histnames = {"h_total_hits",             "h_ldg_shwr_hits",             "h_ldg_shwr_hits_WPlane",
                                          "h_ldg_shwr_Open_Angle",    "h_ldg_shwr_dEdx_WPlane",      "h_ldg_shwr_HitPerLen",
                                          "h_ldg_shwr_Phi",           "h_ldg_shwr_Phi_wrapped",      "h_ldg_shwr_Theta",       "h_ldg_shwr_CTheta",
                                          "h_long_Track_ldg_shwr",    "h_tpc_obj_vtx_x",             "h_tpc_obj_vtx_y",        "h_tpc_obj_vtx_z",
                                          "h_n_pfp",                  "h_n_pfp_50Hits",              "h_n_tracks",             "h_n_tracks_50Hits",  "h_n_showers",
                                          "h_n_showers_50Hits",       "h_track_phi",                 "h_shower_phi",           "h_largest_flash_y",  "h_largest_flash_z",
                                          "h_largest_flash_time",     "h_largest_flash_pe",          "h_Flash_TPCObj_Dist",
                                          "h_shower_Nu_vtx_Dist",     "h_track_Nu_vtx_Dist",         "h_selected",
                                          "h_shower_phi_pi0",         "h_shower_phi_bkg_cosmic",         "h_shower_phi_other",
                                          "h_shower_phi_pi0_wrapped", "h_shower_phi_bkg_cosmic_wrapped", "h_shower_phi_other_wrapped",
                                          "h_shower_E_pi0",           "h_shower_E_bkg_cosmic",           "h_shower_E_other",       "h_shower_E",
                                          "h_shower_Theta_pi0",       "h_shower_Theta_bkg_cosmic",       "h_shower_Theta_other", };

    // Loop over the histograms
    for (int j=0; j < histnames.size(); j++){
        TH1D* hist;
        
        // Canvas + Legend
        TCanvas* c = new TCanvas();
        TLegend* legend;
        if (histnames[j] == "h_ldg_shwr_CTheta") legend = new TLegend(0.15, 0.55, 0.37, 0.85); // Reposition
        else if (histnames[j] == "h_tpc_obj_vtx_x" || histnames[j] == "h_tpc_obj_vtx_y" || histnames[j] == "h_tpc_obj_vtx_z" )
                legend = new TLegend(0.15, 0.15, 0.37, 0.45); // Reposition
        else if (histnames[j] == "h_largest_flash_z" ) legend = new TLegend(0.35, 0.59, 0.57, 0.89);
        else legend = new TLegend(0.72, 0.59, 0.94, 0.89);

        legend->SetBorderSize(0);
        legend->SetFillStyle(0);

        // Loop over variation directories
        for (int i=0; i < variations.size(); i++){
            char name[500];
            snprintf(name, 500, "%s/%s", variations[i].c_str(),histnames[j].c_str() );

            hist = (TH1D*)f_var_out->Get(name);
            if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!" << std::endl;

            DrawTH1D_SAME(hist, variations[i], legend, histnames[j]);


        }

        // Print the Canvas
        char Canvas_name[500];
        snprintf(Canvas_name, 500, "plots/%s/%s.pdf", mode_str.c_str(), histnames[j].c_str() ); 
        legend->Draw();
        c->Print(Canvas_name);
        snprintf(Canvas_name, 500, "plots_png/%s/%s.png", mode_str.c_str(), histnames[j].c_str() );
        c->Print(Canvas_name);

    }
    
    // ************************ 1D Histograms ratio ****************************
    std::string histname_ratio;
    std::vector<std::string> histnames_phi_theta = {"h_ldg_shwr_Phi", "h_ldg_shwr_Theta", "h_ldg_shwr_Phi_wrapped", "h_selected"};
    
    for (int j=0; j < histnames_phi_theta.size(); j++){
        TH1D* hist;
        TH1D* hist_CV;
        
        // Canvas + Legend
        TCanvas* c = new TCanvas();
        TLegend* legend = new TLegend(0.72, 0.59, 0.94, 0.89);

        legend->SetBorderSize(0);
        legend->SetFillStyle(0);

        // Loop over variation directories
        for (int i=0; i < variations.size(); i++){
            char name[500];
            snprintf(name, 500, "%s/%s", variations[i].c_str(), histnames_phi_theta[j].c_str() );
            hist = (TH1D*)f_var_out->Get(name);
            if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!" << std::endl;

            char name_CV[500];
            
            if (mode_str == "bnb")
                snprintf(name_CV, 500, "BNBCV/%s", histnames_phi_theta[j].c_str() );
            else
                snprintf(name_CV, 500, "NuMICV/%s", histnames_phi_theta[j].c_str() );
            
            
            hist_CV = (TH1D*)f_var_out->Get(name_CV);
            if (hist_CV == NULL ) std::cout << "ERROR: Can't get CV Histogram!" << std::endl;

            histname_ratio = histnames_phi_theta[j] + "_ratio";
            
            // hist->Sumw2();
            // hist_CV->Sumw2();

            TH1D* hist_divide = (TH1D*) hist->Clone("hist_divide");
            hist_divide->Sumw2();
            
            hist_divide->Divide(hist_CV);

            if (variations[i] == "BNBCV" || variations[i] == "NuMICV"){
                for (unsigned int k=1; k < hist_divide->GetNbinsX()+1;     k++ ){
                    hist_divide->SetBinContent(k, 1);
                }
            }
            
            DrawTH1D_SAME(hist_divide, variations[i], legend, histname_ratio);
            
        }
        
        // Print the Canvas
        char Canvas_name[500];
        snprintf(Canvas_name, 500, "plots/%s/%s.pdf",mode_str.c_str(), histname_ratio.c_str() ); 
        legend->Draw();
        c->Print(Canvas_name);

        snprintf(Canvas_name, 500, "plots_png/%s/%s.png",mode_str.c_str(), histname_ratio.c_str() ); 
        c->Print(Canvas_name);

    }
    // ************************** 2D Histograms ********************************
    std::vector<std::string> histnames_2D = {"h_EBkg_Theta",       "h_EBkg_Phi",       "h_EBkg_Phi_wrapped",        "h_ThetaBkg_Phi_wrapped",        "h_ThetaBkg_Phi", 
                                             "h_EBkg_pi0_Theta",   "h_EBkg_pi0_Phi",   "h_EBkg_pi0_Phi_wrapped",    "h_ThetaBkg_pi0_Phi_wrapped",    "h_ThetaBkg_pi0_Phi", 
                                             "h_EBkg_cosmic_Theta","h_EBkg_cosmic_Phi","h_EBkg_cosmic_Phi_wrapped", "h_ThetaBkg_cosmic_Phi_wrapped", "h_ThetaBkg_cosmic_Phi", 
                                             "h_EBkg_other_Theta", "h_EBkg_other_Phi", "h_EBkg_other_Phi_wrapped",  "h_ThetaBkg_other_Phi_wrapped",  "h_ThetaBkg_other_Phi"};
    

    // Loop over the histograms
    for (int j=0; j < histnames_2D.size(); j++){
        TH2D* hist;
        
        // Loop over variation directories
        for (int i=0; i < variations.size(); i++){
            // Canvas + Legend
            TCanvas* c = new TCanvas();

            char name[500];
            snprintf(name, 500, "%s/%s", variations[i].c_str(),histnames_2D[j].c_str() );

            hist = (TH2D*)f_var_out->Get(name);
            if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!" << std::endl;

            DrawTH2D_SAME(hist, variations[i], histnames_2D[j]);

            // Print the Canvas
            char Canvas_name[1000];
            snprintf(Canvas_name, 1000, "plots/%s/%s_%s.pdf",mode_str.c_str(),histnames_2D[j].c_str(), variations[i].c_str() ); 
            c->Print(Canvas_name);
            snprintf(Canvas_name, 1000, "plots_png/%s/%s_%s.png",mode_str.c_str(),histnames_2D[j].c_str(), variations[i].c_str() ); 
            c->Print(Canvas_name);
            c->Close();
        }
    }
    
    // Close the file
    f_var_out->Close();

}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::PlotVariatonsNuMIBNB(){
    
    // Get the files
    TFile *f_bnb, *f_numi;

    // create plots folder if it does not exist
    system("if [ ! -d \"plots/numi_bnb_comparisons/CV/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/numi_bnb_comparisons/CV; fi");
    system("if [ ! -d \"plots_png/numi_bnb_comparisons/CV/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots_png/numi_bnb_comparisons/CV/; fi");

    system("if [ ! -d \"plots/numi_bnb_comparisons/DIC/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/numi_bnb_comparisons/DIC; fi");
    system("if [ ! -d \"plots_png/numi_bnb_comparisons/DIC/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots_png/numi_bnb_comparisons/DIC/; fi");

    system("if [ ! -d \"plots/numi_bnb_comparisons/LArG4BugFix/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/numi_bnb_comparisons/LArG4BugFix; fi");
    system("if [ ! -d \"plots_png/numi_bnb_comparisons/LArG4BugFix/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots_png/numi_bnb_comparisons/LArG4BugFix/; fi");
    
    f_bnb  = new TFile("plots/variation_out_bnb_bkg.root" ,"UPDATE");
    f_numi = new TFile("plots/variation_out_numi_bkg.root","UPDATE");
   
    
    // Grab the variation folders in the file
    std::vector<std::string> variations = {"NuMICV", "BNBCV", "NuMIwithDIC", "BNBwithDIC", "NuMILArG4BugFix", "BNBLArG4BugFix"};

     // ************************** 1D Histograms *******************************
    std::vector<std::string> histnames = {"h_total_hits",             "h_ldg_shwr_hits",             "h_ldg_shwr_hits_WPlane",
                                         "h_ldg_shwr_Open_Angle",     "h_ldg_shwr_dEdx_WPlane",      "h_ldg_shwr_HitPerLen",
                                         "h_ldg_shwr_Phi",            "h_ldg_shwr_Phi_wrapped",      "h_ldg_shwr_Theta",     "h_ldg_shwr_CTheta",
                                          "h_long_Track_ldg_shwr",    "h_tpc_obj_vtx_x",             "h_tpc_obj_vtx_y",      "h_tpc_obj_vtx_z",
                                          "h_n_pfp",                  "h_n_pfp_50Hits",              "h_n_tracks",           "h_n_tracks_50Hits", "h_n_showers",
                                          "h_n_showers_50Hits",       "h_track_phi",                 "h_shower_phi",         "h_largest_flash_y", "h_largest_flash_z",
                                          "h_largest_flash_time",     "h_largest_flash_pe",          "h_Flash_TPCObj_Dist",
                                          "h_shower_Nu_vtx_Dist",     "h_track_Nu_vtx_Dist",         "h_selected",
                                          "h_shower_phi_pi0",         "h_shower_phi_other",
                                          "h_shower_phi_pi0_wrapped", "h_shower_phi_other_wrapped",
                                          "h_shower_E_pi0",           "h_shower_E_other",            "h_shower_E",
                                          "h_shower_phi_bkg_cosmic", "h_shower_phi_bkg_cosmic_wrapped", "h_shower_E_bkg_cosmic"
                                          };

    std::vector<std::string> histnames_ratio = {
        "h_ldg_shwr_Phi",           "h_ldg_shwr_Phi_wrapped",      "h_ldg_shwr_Theta",
        "h_shower_phi_pi0_wrapped", "h_shower_phi_other_wrapped",  "h_shower_phi_bkg_cosmic_wrapped", 
        "h_shower_E_pi0",           "h_shower_E_bkg_cosmic",       "h_shower_E_other", "h_shower_E"
    };

    // Loop over the histograms
    for (int j=0; j < histnames.size(); j++){
        std::vector<TH1D*> hist(variations.size());
        
        // Canvas + Legend
        TCanvas* c = new TCanvas();
        TLegend* legend = new TLegend(0.72, 0.59, 0.94, 0.89);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);

        TLegend* legend_ratio = new TLegend(0.72, 0.59, 0.94, 0.89);
        legend_ratio->SetBorderSize(0);
        legend_ratio->SetFillStyle(0);

        // Loop over variation directories and get the histograms
        for (int i=0; i < variations.size(); i++){
            char name[500];
            snprintf(name, 500, "%s/%s", variations[i].c_str(),histnames[j].c_str() );

            hist.at(i) = (TH1D*)f_bnb->Get(name);
            if (hist.at(i) == NULL ) hist.at(i) = (TH1D*)f_numi->Get(name); // If no histogram in bnb file, try the numi one
            if (hist.at(i) == NULL ) std::cout << "ERROR: Can't get Histogram! " << name << std::endl;

        }

        // Get the histogram areas
        std::vector<double> hist_areas(variations.size());
        for (int i=0; i < variations.size(); i++){
           hist_areas.at(i) = hist.at(i)->Integral(0, -1);
        }

        // Now we want to normalise the histograms by area -- scale the pairs of variations to the numi
        for (int i=1; i < variations.size(); i+=2){
           hist.at(i)->Scale(hist_areas.at(i - 1) / hist_areas.at(i));
        }

        // Look to see if we want to make a ratio plot
        bool bool_string{false};
        for (int i = 0; i < histnames_ratio.size(); i++){
        
            size_t found = histnames.at(j).find(histnames_ratio.at(i)); 
            
            // Got the variation match
            if (found != std::string::npos) {
                bool_string = true;
            }
        }

        // Divde out the bnb to numi histgram
        TH1D *hist_divide = (TH1D*) hist.at(1)->Clone("hist_divide"); // BNB CV
        if (bool_string) {
            hist_divide->Sumw2();
            hist_divide->Divide(hist.at(0)); // NuMI CV
        }

        TH1D *hist_divide_withDIC = (TH1D*) hist.at(3)->Clone("hist_divide_withDIC"); // BNB withDIC
        if (bool_string) {
            hist_divide_withDIC->Sumw2();
            hist_divide_withDIC->Divide(hist.at(2)); // NuMI withDIC
        }

        TH1D *hist_divide_LArG4BugFix = (TH1D*) hist.at(5)->Clone("hist_divide_LArG4BugFix"); // BNB LArG4BugFix
        if (bool_string) {
            hist_divide_LArG4BugFix->Sumw2();
            hist_divide_LArG4BugFix->Divide(hist.at(4)); // NuMI LArG4BugFix
        }


        char Canvas_name[500];

        // Now draw the histograms, one plot for each pair of variaions
        for (int i=0; i < variations.size(); i+=2){
            c->Clear();
            c = new TCanvas();

            legend->Clear();
            legend = new TLegend(0.72, 0.59, 0.94, 0.89);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            
            DrawTH1D_SAME(hist.at(i),   variations[i],   legend, histnames[j]);
            DrawTH1D_SAME(hist.at(i+1), variations[i+1], legend, histnames[j]);
        
            std::string var_str;

            if (i == 0)      var_str = "CV";
            else if (i == 2) var_str = "DIC";
            else             var_str = "LArG4BugFix";

            // Print the Canvas
            snprintf(Canvas_name, 500, "plots/numi_bnb_comparisons/%s/%s_%s.pdf",var_str.c_str(), histnames[j].c_str(), var_str.c_str() ); 
            legend->Draw();
            c->Print(Canvas_name);

            snprintf(Canvas_name, 500, "plots_png/numi_bnb_comparisons/%s/%s_%s.png",var_str.c_str(), histnames[j].c_str(), var_str.c_str() ); 
            c->Print(Canvas_name);
        
        }

        // Plot the ratios, 1 plot for all variations
        if (bool_string) {
            c->Clear();
            c = new TCanvas();
            DrawTH1D_Ratio(hist_divide,             "CV",          legend_ratio, histnames[j]);
            DrawTH1D_Ratio(hist_divide_withDIC,     "withDIC",     legend_ratio, histnames[j]);
            DrawTH1D_Ratio(hist_divide_LArG4BugFix, "LArG4BugFix", legend_ratio, histnames[j]);
            legend_ratio->Draw();
            snprintf(Canvas_name, 500, "plots/numi_bnb_comparisons/%s_ratio.pdf",histnames[j].c_str() );
            c->Print(Canvas_name);
            snprintf(Canvas_name, 500, "plots_png/numi_bnb_comparisons/%s_ratio.png",histnames[j].c_str() );
            c->Print(Canvas_name);
        }
    }
    
    // ************************** 2D Histograms ********************************
    std::vector<std::string> histnames_2D = {"h_EBkg_Theta",       "h_EBkg_Phi",       "h_EBkg_Phi_wrapped",        "h_ThetaBkg_Phi_wrapped",        "h_ThetaBkg_Phi", 
                                             "h_EBkg_pi0_Theta",   "h_EBkg_pi0_Phi",   "h_EBkg_pi0_Phi_wrapped",    "h_ThetaBkg_pi0_Phi_wrapped",    "h_ThetaBkg_pi0_Phi", 
                                             "h_EBkg_cosmic_Theta","h_EBkg_cosmic_Phi","h_EBkg_cosmic_Phi_wrapped", "h_ThetaBkg_cosmic_Phi_wrapped", "h_ThetaBkg_cosmic_Phi", 
                                             "h_EBkg_other_Theta", "h_EBkg_other_Phi", "h_EBkg_other_Phi_wrapped",  "h_ThetaBkg_other_Phi_wrapped",  "h_ThetaBkg_other_Phi"};

    // Loop over the histograms
    for (int j=0; j < histnames_2D.size(); j++){
        TH2D* hist;
        
        // Loop over variation directories
        for (int i=0; i < variations.size(); i++){
            // Canvas + Legend
            TCanvas* c = new TCanvas();

            char name[500];
            snprintf(name, 500, "%s/%s", variations[i].c_str(),histnames_2D[j].c_str() );

            hist = (TH2D*)f_bnb->Get(name);
            if (hist == NULL ) hist = (TH2D*)f_numi->Get(name); // If no histogram in bnb file, try the numi one
            if (hist == NULL ) std::cout << "ERROR: Can't get Histogram! " << name << std::endl;

            DrawTH2D_SAME(hist, variations[i], histnames_2D[j]);

            std::string var_str;

            if (i == 0)      var_str = "CV";
            else if (i == 2) var_str = "DIC";
            else             var_str = "LArG4BugFix";

            // Print the Canvas
            char Canvas_name[1000];
            snprintf(Canvas_name, 1000, "plots/numi_bnb_comparisons/%s/%s_%s.pdf",var_str.c_str(), histnames_2D[j].c_str(), variations[i].c_str() ); 
            c->Print(Canvas_name);
            snprintf(Canvas_name, 1000, "plots_png/numi_bnb_comparisons/%s/%s_%s.png",var_str.c_str(), histnames_2D[j].c_str(), variations[i].c_str() ); 
            c->Print(Canvas_name);
            c->Close();
        }
    }
    
    // Close the file
    f_bnb ->Close();
    f_numi->Close();

}
//***************************************************************************
//***************************************************************************
int variation_output_bkg::GetLeadingShowerIndex(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj){
    int leading_index{0};
    int most_hits{0};
    
    // Loop over Particle Objects
    for (int j = 0; j < n_pfp; j++) {
        auto const pfp_obj = tpc_obj.GetParticle(j);
        const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
        
        if (pfp_pdg != 11) continue; // skip if not a shower

        const int n_pfp_hits = pfp_obj.NumPFPHits();
        
        // Compare for most hits
        if (n_pfp_hits > most_hits) {
            leading_index = j; 
            most_hits = n_pfp_hits; 
        }
    }

    return leading_index;
}
//***************************************************************************
//***************************************************************************
double variation_output_bkg::GetLongestTrackLength(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj){
    double longest_track = 0;

    // Loop over Particle Objects
    for(int i = 0; i < n_pfp; i++) {

        auto const pfp_obj = tpc_obj.GetParticle(i);
        const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
        
        // If it  is track like
        if(pfp_pdg == 13) {
            const double trk_length = pfp_obj.pfpLength();
            
            if(trk_length > longest_track) {

                longest_track = trk_length;
            }
        }
    }

    return longest_track;
}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::GetNumber_Track_Shower(const int n_pfp, int n_tpc_obj,
                                     xsecAna::TPCObjectContainer tpc_obj, int &n_showers, int &n_tracks,
                                     int &n_pfp_50Hits, int &n_tracks_50Hits,int &n_showers_50Hits) { 
    // Loop over Particle Objects
    for (int j = 0; j < n_pfp; j++) {
        auto const pfp_obj = tpc_obj.GetParticle(j);
        const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
        
        // Don't count the neutrinos
        if (pfp_pdg != 12 || pfp_pdg != -12 || pfp_pdg != 14 || pfp_pdg != -14){
            // Cut with pfp with > 50 Hits
            if (pfp_obj.NumPFPHits() > 50){
                // For > 50 hits
                n_pfp_50Hits++; // Add one to n_pfp > 50 hits counter
                
                if (pfp_pdg == 11)         n_showers_50Hits++; // Add to shower counter
                else if (pfp_pdg == 13) n_tracks_50Hits++;  // Add to track counter
                else std::cout << "Unknown pandora classification:\t" << pfp_pdg << std::endl;
                    
            }

            // All PFP
            if (pfp_pdg == 11)         n_showers++; // Add to shower counter
            else if (pfp_pdg == 13) n_tracks++;  // Add to track counter
            else return;
        }
        
    }
}
//***************************************************************************
//***************************************************************************
double variation_output_bkg::pfp_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z) {
    const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );

    return distance;
}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::DrawTH1D(TH1D* h, double POT_Scaling){
    TCanvas* c = new TCanvas();
    c->cd();

    h->SetLineColor(kMagenta+3);
    h->SetLineWidth(2);
    h->SetLineStyle(1);
    h->Scale(POT_Scaling);
    h->SetOption("hist, E");
    h->Draw("hist, E");

}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::DrawTH2D(TH2D* h, double POT_Scaling){
    TCanvas* c = new TCanvas();
    c->cd();

    h->Scale(POT_Scaling);
    h->Draw("colz");
    c->Close();

}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::DrawTH1D_SAME(TH1D* hist, std::string variation, TLegend* legend, std::string histname){
    
    // ---------------------- 
    //    Axis Specifiers
    // ----------------------
    if (histname == "h_total_hits"){
        hist->SetTitle(";PFP Total Hits;Entries");
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_ldg_shwr_hits") {
        hist->SetTitle("; Leading Shower Hits (All Planes);Entries");
        // hist->GetYaxis()->SetRangeUser(0,750);
    }
    else if (histname == "h_ldg_shwr_hits_WPlane") {
        hist->SetTitle(";Leading Shower Collection Plane Hits;Entries");
        // hist->GetYaxis()->SetRangeUser(0,5000);
    }
    else if (histname == "h_ldg_shwr_Open_Angle"){
        hist->SetTitle(";Leading Shower Opening Angle [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,175);
    }
    else if (histname == "h_ldg_shwr_dEdx_WPlane"){
        hist->SetTitle(";Leading Shower dEdx Collection Plane [MeV/cm];Entries");
        hist->GetYaxis()->SetRangeUser(0,135);
    }
    else if (histname == "h_ldg_shwr_HitPerLen"){
        hist->SetTitle(";Leading Shower Hits / Length [ cm^{-1} ];Entries");
        // hist->GetYaxis()->SetRangeUser(0,7000);
    }
    else if (histname == "h_ldg_shwr_Phi"){
        hist->SetTitle(";Leading Shower #phi [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,250);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,50);
    }
    else if (histname == "h_ldg_shwr_Theta"){
        hist->SetTitle(";Leading Shower #theta [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,200);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,60);
    }
    else if (histname == "h_ldg_shwr_Phi_ratio"){
        hist->SetTitle(";Leading Shower #phi [degrees];Ratio");
        hist->GetYaxis()->SetRangeUser(-2,6);
    }
    else if (histname == "h_ldg_shwr_Phi_wrapped"){
        hist->SetTitle(";Leading Shower #phi (wrapped) [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,400);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,60);
    }
    else if (histname == "h_ldg_shwr_Phi_wrapped_ratio"){
        hist->SetTitle(";Leading Shower #phi (wrapped) ratio [degrees];Ratio");
        hist->GetYaxis()->SetRangeUser(-0.2,3);
    }
    else if (histname == "h_ldg_shwr_Theta_ratio"){
        hist->SetTitle(";Leading Shower #theta [degrees];Ratio");
        hist->GetYaxis()->SetRangeUser(-2,5);
    }
    else if (histname == "h_ldg_shwr_CTheta"){
        hist->SetTitle(";Leading Shower cos(#theta);Entries");
        hist->GetYaxis()->SetRangeUser(0,300);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,60);
    }
    else if (histname == "h_long_Track_ldg_shwr"){
        hist->SetTitle(";Longest Track Length / Leading Shower Length;Entries");
        // hist->GetYaxis()->SetRangeUser(0,22000);
    }
    else if (histname == "h_tpc_obj_vtx_x"){
        hist->SetTitle(";TPC Object Vertex X [cm];Entries");
        // hist->GetYaxis()->SetRangeUser(0,2500);
    }
    else if (histname == "h_tpc_obj_vtx_y"){
        hist->SetTitle(";TPC Object Vertex Y [cm];Entries");
        // hist->GetYaxis()->SetRangeUser(0,2200);
    }
    else if (histname == "h_tpc_obj_vtx_z"){
        hist->SetTitle(";TPC Object Vertex Z [cm];Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_n_pfp"){
        hist->SetTitle(";Number of PFP;Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_n_pfp_50Hits"){
        hist->SetTitle("; Number of PFP > 50 Hits;Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_n_tracks"){
        hist->SetTitle(";Number of Tracks;Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_n_tracks_50Hits"){
        hist->SetTitle("; Number of Tracks > 50 Hits;Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_n_showers"){
        hist->SetTitle(";Number of Showers;Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_n_showers_50Hits"){
        hist->SetTitle("; Number of Showers > 50 Hits;Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_track_phi"){
        hist->SetTitle("; Track #phi [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,60);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,20);
    }
    else if (histname == "h_shower_phi"){
        hist->SetTitle("; Shower #phi [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,800);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,300);
    }            
    else if (histname == "h_largest_flash_y"){
        hist->SetTitle(";Largest Flash Y [cm];Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_largest_flash_z"){
        hist->SetTitle("; Largest Flash Z [cm];Entries");
        // hist->GetYaxis()->SetRangeUser(0,1000);
    }
    else if (histname == "h_largest_flash_time"){
        hist->SetTitle("; Largest Flash Time [#mus];Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_largest_flash_pe"){
        hist->SetTitle("; Largest Flash PE;Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_Flash_TPCObj_Dist"){
        hist->SetTitle("; 2D Distance from Largest Flash to #nu Vertex [cm];Entries");
        // hist->GetYaxis()->SetRangeUser(0,1200);
    }
    else if (histname == "h_shower_Nu_vtx_Dist"){
        hist->SetTitle("; 3D Distance of Shower to #nu Vertex [cm];Entries");
        hist->GetYaxis()->SetRangeUser(0,1000);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,300);
    }
    else if (histname == "h_track_Nu_vtx_Dist"){
        hist->SetTitle("; 3D Distance of Track to #nu Vertex [cm];Entries");
        hist->GetYaxis()->SetRangeUser(0,150);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,50);
    }
    else if (histname == "h_selected"){
        hist->SetTitle("; Num Background Selected;Entries");
        hist->GetXaxis()->SetLabelOffset(999);
        hist->GetXaxis()->SetLabelSize(0);
        hist->GetYaxis()->SetRangeUser(400,700);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,200); 
    }
    else if (histname == "h_selected_ratio"){
        hist->SetTitle("; Num Background Selected;Ratio");
        hist->GetXaxis()->SetLabelOffset(999);
        hist->GetXaxis()->SetLabelSize(0);
        hist->GetYaxis()->SetRangeUser(0,2);
    }
    else if (histname ==  "h_shower_phi_pi0"){
        hist->SetTitle("; Leading Shower Phi #pi^{0} Bkg [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,50);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,30);
    }
    else if (histname ==  "h_shower_phi_bkg_cosmic"){
        hist->SetTitle("; Leading Shower Phi cosmic Bkg [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,50);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,20);
    }
    else if (histname ==  "h_shower_phi_other"){
        hist->SetTitle("; Leading Shower Phi Other Bkg [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,50);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,10);
    }
    else if (histname ==  "h_shower_phi_pi0_wrapped"){
        hist->SetTitle("; Leading Shower Phi Wrapped #pi^{0} Bkg [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,220);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,60);
    }
    else if (histname ==  "h_shower_phi_bkg_cosmic_wrapped"){
        hist->SetTitle("; Leading Shower Phi Wrapped cosmic Bkg [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,40);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,15);
    }
    else if (histname ==  "h_shower_phi_other_wrapped"){
        hist->SetTitle("; Leading Shower Phi Wrapped Other Bkg [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,100);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,20);
    }
    else if (histname ==  "h_shower_E_pi0"){
        hist->SetTitle("; Leading Shower Energy #pi^{0} Bkg [GeV];Entries");
        hist->GetYaxis()->SetRangeUser(0,300);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,80);
    }
    else if (histname ==  "h_shower_E_bkg_cosmic"){
        hist->SetTitle("; Leading Shower Energy cosmic Bkg [GeV];Entries");
        // hist->GetYaxis()->SetRangeUser(0,10);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,10);
    }
    else if (histname ==  "h_shower_E_other"){
        hist->SetTitle("; Leading Shower Energy Other Bkg [GeV];Entries");
        hist->GetYaxis()->SetRangeUser(0,100);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,30);
    }
    else if (histname ==  "h_shower_E"){
        hist->SetTitle("; Leading Shower Energy Bkg [GeV];Entries");
        hist->GetYaxis()->SetRangeUser(0,300);

        if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,110);
    }
    else if (histname ==  "h_shower_Theta_pi0"){
        hist->SetTitle("; Leading Shower Theta #pi^{0} Bkg [Degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,300);

        // if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,80);
    }
    else if (histname ==  "h_shower_Theta_bkg_cosmic"){
        hist->SetTitle("; Leading Shower Theta cosmic Bkg [Degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,10);

        // if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,10);
    }
    else if (histname ==  "h_shower_Theta_other"){
        hist->SetTitle("; Leading Shower Theta Other Bkg [Degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,100);

        // if (draw_mode == "bnb") hist->GetYaxis()->SetRangeUser(0,30);
    }
    else return;

    std::string draw_spec = "hist,E, same";
    hist->SetOption("hist, E");

    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);

    if ((histname == "h_ldg_shwr_Theta_ratio" || histname == "h_ldg_shwr_Phi_ratio" || histname == "h_ldg_shwr_Phi_wrapped_ratio") && variation != "BNBCV"  )
        draw_spec = "E, same";

    // ----------------------
    //    Draw Specifiers
    // ----------------------
    if (variation == "BNBCV" || variation == "NuMICV"){
        hist->SetLineColor(kBlack);
        hist->SetLineWidth(2);
        hist->SetLineStyle(1);
        

        if (variation == "NuMICV"){
            hist->SetLineColor(kAzure-6);
            legend->AddEntry(hist, "NuMI CV", "l");
        }
        else
            legend->AddEntry(hist, "BNB CV", "l");

            hist->Draw(draw_spec.c_str());
    } 
    else if  (variation == "BNBwithDIC"         || variation == "NuMIwithDIC"){
        hist->SetLineColor(kMagenta+2);
        hist->SetLineWidth(2);
        hist->SetLineStyle(1);
        legend->AddEntry(hist, "DIC", "l");
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBEnhancedTPCVis"  || variation == "NuMIEnhancedTPCVis" ){ 
        hist->SetLineColor(30);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Enhanced TPC Vis.", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBaltDeadChannels" || variation == "NuMIaltDeadChannels"){ 
        hist->SetLineColor(38);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Alt. Dead Chan.", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBdeadSaturatedChannels" || variation == "NuMIdeadSaturatedChannels"){
        hist->SetLineColor(28);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Dead Sat. Chan.", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
        
    }
    else if  (variation == "BNBstretchResp"  || variation == "NuMIstretchResp" ){
        hist->SetLineColor(36);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Stretch Resp.", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBsqueezeResp"  || variation == "NuMIsqueezeResp"){
        hist->SetLineColor(1001);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Squeeze Resp.", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBupPEnoise"    || variation == "NuMIupPEnoise"){
        hist->SetLineColor(kBlue+1);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "PE Noise Up", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBnoiseAmpDown" || variation == "NuMInoiseAmpDown"){
        hist->SetLineColor(42);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Noise Amp. Down", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBdownPEnoise"  || variation == "NuMIdownPEnoise"){
        hist->SetLineColor(50);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "PE Noise Down", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBnoiseAmpUp"   || variation == "NuMInoiseAmpUp"){
        hist->SetLineColor(kOrange+10);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Noise Amp. Up", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBDTdown"       || variation == "NuMIDTdown"){
        hist->SetLineColor(kOrange+1);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "DT Down", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBDTup"         || variation == "NuMIDTup"){
        hist->SetLineColor(kMagenta-10);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "DT Up", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBDLup"         || variation == "NuMIDLup"){
        hist->SetLineColor(kMagenta);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "DL Up", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBDLdown"       || variation == "NuMIDLdown"){
        hist->SetLineColor(kTeal+6);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "DL Down", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBdataSCE"      || variation == "NuMIdataSCE"){
        hist->SetLineColor(kAzure-9);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "SCE", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBLArG4BugFix"  || variation == "NuMILArG4BugFix"){
        hist->SetLineColor(kSpring-7);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "LArG4BugFix", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBBirksRecomb"  || variation == "NuMIBirksRecomb"){
        hist->SetLineColor(kRed+1);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Birks Recomb.","l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else {
        std::cout << "Error! Could not match varations:\t" << variation << std::endl;
        return;
    }

}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::DrawTH1D_Ratio(TH1D* hist, std::string variation, TLegend* legend, std::string histname){

    TLine *line;

    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);

    if (histname == "h_ldg_shwr_Phi"){
        hist->SetTitle(";Leading Shower #phi [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,40);
        line = new TLine(-180,1,180,1); 
    }
    else if (histname == "h_ldg_shwr_Theta"){
        hist->SetTitle(";Leading Shower #theta [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,7000);
        line = new TLine(0,1,180,1); 
    }
    else if (histname == "h_ldg_shwr_Phi_wrapped"){
        hist->SetTitle(";Leading Shower #phi (wrapped) [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,40);
        line = new TLine(0,1,90,1); 
    }
    else if (histname ==  "h_shower_phi_pi0"){
        hist->SetTitle("; Leading Shower Phi #pi^{0} Bkg [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,50);
        line = new TLine(-180,1,180,1); 
    }
    else if (histname ==  "h_shower_phi_bkg_cosmic"){
        hist->SetTitle("; Leading Shower Phi cosmic Bkg [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,10);
        line = new TLine(-180,1,180,1); 
    }
    else if (histname ==  "h_shower_phi_other"){
        hist->SetTitle("; Leading Shower Phi Other Bkg [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,30);
        line = new TLine(-180,1,180,1);  
    }
    else if (histname ==  "h_shower_phi_pi0_wrapped"){
        hist->SetTitle("; Leading Shower Phi Wrapped #pi^{0} Bkg [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,120);
        line = new TLine(0,1,90,1); 
    }
    else if (histname ==  "h_shower_phi_bkg_cosmic_wrapped"){
        hist->SetTitle("; Leading Shower Phi Wrapped cosmic Bkg [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,15);
        line = new TLine(0,1,90,1); 
    }
    else if (histname ==  "h_shower_phi_other_wrapped"){
        hist->SetTitle("; Leading Shower Phi Wrapped Other Bkg [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,70);
        line = new TLine(0,1,90,1); 
    }
    else if (histname ==  "h_shower_E_pi0"){
        hist->SetTitle("; Leading Shower Energy #pi^{0} Bkg [GeV];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,50);
        line = new TLine(0,1,3,1); 
    }
    else if (histname ==  "h_shower_E_bkg_cosmic"){
        hist->SetTitle("; Leading Shower Energy cosmic Bkg [GeV];Ratio BNB / NuMI");
        hist->GetYaxis()->SetRangeUser(0,10);
        line = new TLine(0,1,5,1);  
    }
    else if (histname ==  "h_shower_E_other"){
        hist->SetTitle("; Leading Shower Energy Other Bkg [GeV];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,30);
        line = new TLine(0,1,3,1); 
    }
    else if (histname ==  "h_shower_E"){
        hist->SetTitle("; Leading Shower Bkg Energy [GeV];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,30);
        line = new TLine(0,1,3,1); 
    }
    else if (histname ==  "h_shower_Theta_pi0"){
        hist->SetTitle("; Leading Shower Theta #pi^{0} Bkg [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,50);
        line = new TLine(0,1,3,1); 
    }
    else if (histname ==  "h_shower_Theta_bkg_cosmic"){
        hist->SetTitle("; Leading Shower Theta cosmic Bkg [degrees];Ratio BNB / NuMI");
        hist->GetYaxis()->SetRangeUser(0,10);
        line = new TLine(0,1,5,1);  
    }
    else if (histname ==  "h_shower_Theta_other"){
        hist->SetTitle("; Leading Shower Theta Other Bkg [degrees];Ratio BNB / NuMI");
        // hist->GetYaxis()->SetRangeUser(0,30);
        line = new TLine(0,1,3,1); 
    }
    


    // ----------------------
    //    Draw Specifiers
    // ----------------------
    if (variation == "CV"){
        hist->SetLineColor(kBlack);
        hist->SetLineWidth(2);
        hist->SetLineStyle(1);
        legend->AddEntry(hist, "CV", "l");
        
    } 
    else if  (variation == "withDIC"){
        hist->SetLineColor(kMagenta+2);
        hist->SetLineWidth(2);
        hist->SetLineStyle(1);
        legend->AddEntry(hist, "DIC", "l");

    }
    else if  (variation == "LArG4BugFix"){
        hist->SetLineColor(kSpring-7);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "LArG4BugFix", "l");
        hist->SetLineStyle(1);
    }
    else {
        std::cout << "Error! Could not match varations:\t" << variation << std::endl;
        return;
    }

    hist->GetYaxis()->SetRangeUser(-2,5);

    // hist->SetOption("E");
    hist->Draw("E, same");

    
    line->SetLineStyle(7);
    line->Draw();

}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::DrawTH2D_SAME(TH2D* hist, std::string variation, std::string histname){
    
    // IncreaseLabelSize(hist);

    std::string draw_spec;

    // Find the ratio in the string
    bool bool_ratio{false};
    size_t found = histname.find("ratio"); 
    
    // Got a ratio plot instead
    if (found != std::string::npos) bool_ratio = true;

    draw_spec = "colz";

    std::string variation_name = " ";

    // ----------------------
    //    Draw Specifiers
    // ----------------------
    if (variation == "BNBCV"){
        variation_name = "BNB CV";
    } 
    else if (variation == "NuMICV"){
         variation_name = "NuMI CV";
    }
    else if  (variation == "BNBwithDIC" || variation == "NuMIwithDIC"){
        variation_name = "DIC";
    }
    else if  (variation == "BNBEnhancedTPCVis" || variation == "NuMIEnhancedTPCVis"){ 
        variation_name = "Enhanced TPC Vis.";
    }
    else if  (variation == "BNBaltDeadChannels" || variation == "NuMIaltDeadChannels"){ 
        variation_name = "Alt. Dead Chan.";
    }
    else if  (variation == "BNBdeadSaturatedChannels" || variation == "NuMIdeadSaturatedChannels"){
        variation_name = "Dead Sat. Chan.";
    }
    else if  (variation == "BNBstretchResp" || variation == "NuMIstretchResp"){
        variation_name = "Stretch Resp.";
    }
    else if  (variation == "BNBsqueezeResp" || variation == "NuMIsqueezeResp"){
        variation_name = "Squeeze Resp.";
    }
    else if  (variation == "BNBupPEnoise" || variation == "NuMIupPEnoise"){
        variation_name = "PE Noise Up";
    }
    else if  (variation == "BNBnoiseAmpDown" || variation == "NuMInoiseAmpDown"){
        variation_name = "Noise Amp. Down";
    }
    else if  (variation == "BNBdownPEnoise" || variation == "NuMIdownPEnoise"){
        variation_name = "PE Noise Down";
    }
    else if  (variation == "BNBnoiseAmpUp" || variation == "NuMInoiseAmpUp"){
        variation_name = "Noise Amp. Up";
    }
    else if  (variation == "BNBDTdown" || variation == "NuMIDTdown"){
        variation_name = "DT Down";
    }
    else if  (variation == "BNBDTup" || variation == "NuMIDTup"){
        variation_name = "DT Up";
    }
    else if  (variation == "BNBDLup" || variation == "NuMIDLup"){
        variation_name = "DL Up";
    }
    else if  (variation == "BNBDLdown" || variation == "NuMIDLdown"){
        variation_name = "DL Down";
    }
    else if  (variation == "BNBdataSCE" || variation == "NuMIdataSCE"){
        variation_name = "SCE";
    }
    else if  (variation == "BNBLArG4BugFix" || variation == "NuMILArG4BugFix"){
        variation_name = "LArG4BugFix";
    }
    else if  (variation == "BNBBirksRecomb" || variation == "NuMIBirksRecomb"){
        variation_name = "Birks Recomb.";
    }
    else {
        std::cout << "Error! Could not match varations:\t" << variation << std::endl;
        return;
    }

    // ----------------------
    //    Axis Specifiers
    // ----------------------

    // All Bkg
    if (histname == "h_EBkg_Theta" || histname == "h_EBkg_Theta_ratio"){
        hist->SetTitle(Form("%s;Leading Shower MC Energy [GeV]; Leading Shower MC Theta [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower MC Energy [GeV]; Leading Shower MC Theta [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_EBkg_Phi" || histname == "h_EBkg_Phi_ratio"){
        hist->SetTitle(Form("%s;Leading Shower MC Energy [GeV]; Leading Shower MC Phi [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower MC Energy [GeV]; Leading Shower MC Phi [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_EBkg_Phi_wrapped" || histname == "h_EBkg_Phi_wrapped_ratio"){
        hist->SetTitle(Form("%s;Leading Shower MC Energy [GeV]; Leading Shower MC Phi Wrapped [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower MC Energy [GeV]; Leading Shower MC Phi Wrapped [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_ThetaBkg_Phi_wrapped" || histname == "h_ThetaBkg_Phi_wrapped_ratio"){
        hist->SetTitle(Form("%s;Leading Shower MC Theta [degrees]; Leading Shower MC Phi Wrapped [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower MC Theta [degrees]; Leading Shower MC Phi Wrapped [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_ThetaBkg_Phi" || histname == "h_ThetaBkg_Phi_ratio"){
        hist->SetTitle(Form("%s;Leading Shower MC Theta [degrees]; Leading Shower MC Phi [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower MC Theta [degrees]; Leading Shower MC Phi [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }

    // Pi0
    else if (histname == "h_EBkg_pi0_Theta" || histname == "h_EBkg_pi0_Theta_ratio"){
        hist->SetTitle(Form("%s;Leading Shower (#pi^{0}) MC Energy [GeV]; Leading Shower (#pi^{0}) MC Theta [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower (#pi^{0}) MC Energy [GeV]; Leading Shower (#pi^{0}) MC Theta [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_EBkg_pi0_Phi" || histname == "h_EBkg_pi0_Phi_ratio"){
        hist->SetTitle(Form("%s;Leading Shower (#pi^{0}) MC Energy [GeV]; Leading Shower (#pi^{0}) MC Phi [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower (#pi^{0}) MC Energy [GeV]; Leading Shower (#pi^{0}) MC Phi [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_EBkg_pi0_Phi_wrapped" || histname == "h_EBkg_pi0_Phi_wrapped_ratio"){
        hist->SetTitle(Form("%s;Leading Shower (#pi^{0}) MC Energy [GeV]; Leading Shower (#pi^{0}) MC Phi Wrapped [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower (#pi^{0}) MC Energy [GeV]; Leading Shower (#pi^{0}) MC Phi Wrapped [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_ThetaBkg_pi0_Phi_wrapped" || histname == "h_ThetaBkg_pi0_Phi_wrapped_ratio"){
        hist->SetTitle(Form("%s;Leading Shower (#pi^{0}) MC Theta [degrees]; Leading Shower (#pi^{0}) MC Phi Wrapped [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower (#pi^{0}) MC Theta [degrees]; Leading Shower (#pi^{0}) MC Phi Wrapped [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_ThetaBkg_pi0_Phi" || histname == "h_ThetaBkg_pi0_Phi_ratio"){
        hist->SetTitle(Form("%s;Leading Shower (#pi^{0}) MC Theta [degrees]; Leading Shower (#pi^{0}) MC Phi [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower (#pi^{0}) MC Theta [degrees]; Leading Shower (#pi^{0}) MC Phi [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
   
   // Cosmic
    else if (histname == "h_EBkg_cosmic_Theta" || histname == "h_EBkg_cosmic_Theta_ratio"){
        hist->SetTitle(Form("%s;Leading Shower Bkg cosmic MC Energy [GeV]; Leading Shower Bkg cosmic MC Theta [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower Bkg cosmic MC Energy [GeV]; Leading Shower Bkg cosmic MC Theta [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_EBkg_cosmic_Phi" || histname == "h_EBkg_cosmic_Phi_ratio"){
        hist->SetTitle(Form("%s;Leading Shower Bkg cosmic MC Energy [GeV]; Leading Shower Bkg cosmic MC Phi [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower Bkg cosmic MC Energy [GeV]; Leading Shower Bkg cosmic MC Phi [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_EBkg_cosmic_Phi_wrapped" || histname == "h_EBkg_cosmic_Phi_wrapped_ratio"){
        hist->SetTitle(Form("%s;Leading Shower Bkg cosmic MC Energy [GeV]; Leading Shower Bkg cosmic MC Phi Wrapped [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower Bkg cosmic MC Energy [GeV]; Leading Shower Bkg cosmic MC Phi Wrapped [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_ThetaBkg_cosmic_Phi_wrapped" || histname == "h_ThetaBkg_cosmic_Phi_wrapped_ratio"){
        hist->SetTitle(Form("%s;Leading Shower cosmic MC Theta [degrees]; Leading Shower cosmic MC Phi Wrapped [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower cosmic MC Theta [degrees]; Leading Shower cosmic MC Phi Wrapped [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_ThetaBkg_cosmic_Phi" || histname == "h_ThetaBkg_cosmic_Phi_ratio"){
        hist->SetTitle(Form("%s;Leading Shower cosmic MC Theta [degrees]; Leading Shower cosmic MC Phi [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower cosmic MC Theta [degrees]; Leading Shower cosmic MC Phi [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
   
   // Other
    else if (histname == "h_EBkg_other_Theta" || histname == "h_EBkg_other_Theta_ratio"){
        hist->SetTitle(Form("%s;Leading Shower Bkg Other MC Energy [GeV]; Leading Shower Bkg Other MC Theta [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower Bkg Other MC Energy [GeV]; Leading Shower Bkg Other MC Theta [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_EBkg_other_Phi" || histname == "h_EBkg_other_Phi_ratio"){
        hist->SetTitle(Form("%s;Leading Shower Bkg Other MC Energy [GeV]; Leading Shower Bkg Other MC Phi [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower Bkg Other MC Energy [GeV]; Leading Shower Bkg Other MC Phi [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_EBkg_other_Phi_wrapped" || histname == "h_EBkg_other_Phi_wrapped_ratio"){
        hist->SetTitle(Form("%s;Leading Shower Bkg Other MC Energy [GeV]; Leading Shower Bkg Other MC Phi Wrapped [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower Bkg Other MC Energy [GeV]; Leading Shower Bkg Other MC Phi Wrapped [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_ThetaBkg_other_Phi_wrapped" || histname == "h_ThetaBkg_other_Phi_wrapped_ratio"){
        hist->SetTitle(Form("%s;Leading Shower Other MC Theta [degrees]; Leading Shower Other MC Phi Wrapped [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower Other MC Theta [degrees]; Leading Shower Other MC Phi Wrapped [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }
    else if (histname == "h_ThetaBkg_other_Phi" || histname == "h_ThetaBkg_other_Phi_ratio"){
        hist->SetTitle(Form("%s;Leading Shower Other MC Theta [degrees]; Leading Shower Other MC Phi [deg]", variation_name.c_str()));
        if (bool_ratio) hist->SetTitle(Form("%s Ratio to CV;Leading Shower Other MC Theta [degrees]; Leading Shower Other MC Phi [deg]", variation_name.c_str()));
        // hist->GetYaxis()->SetRangeUser(0,4500);
    }


    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetZaxis()->SetLabelSize(0.05);
    hist->GetZaxis()->SetTitleSize(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.2);
    gPad->SetBottomMargin(0.13);
    hist->SetMarkerSize(1.8);
   
    
    // gPad->SetLogz();
    hist->Draw(draw_spec.c_str());

}
//***************************************************************************
//***************************************************************************
double variation_output_bkg::GetPOT(const char * _file1){
    double POT{0};
    std::string line;

    std::string filename;
    std::string temp_filename = getFileName(_file1, false);
    filename = "File: " + temp_filename;

    std::cout << filename << std::endl;

    std::ifstream myfile ("POT_List_Varaitions.txt");
    int i_POT{0};
    
    if (myfile.is_open()) {

        // Loop over lines in file
        while ( getline (myfile,line) ) {

            if (i_POT == 1){
                line.erase(0, 21);
                POT = std::stod(line); // Convert string to double
                std::cout << "POT in File:\t"<< POT << std::endl;
                return POT;
            }

            // Found the correct variation file 
            if (line == filename){
                std::cout << "Found match for POT file"<< std::endl;
                i_POT++;
            }
            
        }

        myfile.close();
    }

    else std::cout << "Unable to open file" << std::endl; 
    std::cout << "Could not find a match for POT file"<< std::endl;
    return POT;
}
//***************************************************************************
//***************************************************************************
// Function that grabs the reweighted histogram names for plotting
std::vector<std::string> variation_output_bkg::GrabDirs(TFile* f_var_out) {
    std::vector<std::string> variations;

    f_var_out->cd();
    
    TKey *key;
    TIter nextkey(gDirectory->GetListOfKeys());

    TString same_plots = "SAME_Plots";

    std::cout << "\n=================================================" << std::endl;    
    std::cout << "Getting variation modes:" << std::endl;    
      while ( ( key =  (TKey*)nextkey()) ) { // Extra brackets to omit a warning 
        if (key->IsFolder()) {
            std::cout << key->GetName() << std::endl; // Print the variations
            if (key->GetName() == same_plots ) continue; // Skip this
            variations.push_back(key->GetName());
        }
    }
    std::cout << "=================================================\n" << std::endl;

    return (variations);
}
//***************************************************************************
//************************** Flash Functions ********************************
std::vector<std::vector<double>> variation_output_bkg::GetLargestFlashVector(TFile* f, double flash_start_time, double flash_end_time ){

    // Get Optical Information from file
    TTree* optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

    std::vector<std::vector<double>> largest_flash_v_v;

    // ----------------------
    //    Optical Info
    // ----------------------
    int fRun = 0;
    int fEvent = 0;
    int fOpFlashPE = 0;
    double fOpFlashTime = 0;
    double fOpFlashWidthY = 0;
    double fOpFlashWidthZ = 0;
    double fOpFlashCenterY = 0;
    double fOpFlashCenterZ = 0;

    optical_tree->SetBranchAddress("Run",              &fRun           ); // Run number of Flash
    optical_tree->SetBranchAddress("Event",            &fEvent         ); // Event number of Flash
    optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE     ); // PE of Flash
    optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime   ); // Time of flash
    optical_tree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY );
    optical_tree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ );
    optical_tree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY); // Flash Y center
    optical_tree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ); // Flash Z Center

    // Num events in the optical tree
    const int optical_entries = optical_tree->GetEntries();
    std::cout << "Total Optical Entries: " << optical_entries << std::endl;

    int current_event = 0;
    int current_run = 0;
    int last_event = 0;
    int last_run = 0;

    // Contains the entry number for a given OpFlash per event
    std::vector<int>                    optical_list_pe;
    std::vector<std::vector<int> >        optical_list_pe_v;
    
    std::vector<double>                    optical_list_flash_center_y; 
    std::vector<std::vector<double> >    optical_list_flash_center_y_v;
    
    std::vector<double>                    optical_list_flash_center_z; 
    std::vector<std::vector<double> >    optical_list_flash_center_z_v;
    
    std::vector<double>                    optical_list_flash_time;
    std::vector<std::vector<double> >    optical_list_flash_time_v;
    
    // Loop over the optical entries to get the largest flash vector
    
    // ----------------------
    // Resize the optical enties to be the same sizd as number of Events (TPC Obj)
    // ----------------------
    
    for(int i = 0; i < optical_entries; i++) {
        
        // Get the Optical entry
        optical_tree->GetEntry(i);

        current_run        = fRun;
        current_event     = fEvent;

        // New event
        if(current_event != last_event) {
            optical_list_pe.clear();
            optical_list_flash_center_y.clear();
            optical_list_flash_center_z.clear();
            optical_list_flash_time.clear();

            optical_list_pe.push_back(fOpFlashPE);
            optical_list_flash_center_y.push_back(fOpFlashCenterY);
            optical_list_flash_center_z.push_back(fOpFlashCenterZ);
            optical_list_flash_time.push_back(fOpFlashTime);

        }
        // Same event
        if(current_event == last_event && current_run == last_run) {
            optical_list_pe_v.pop_back();
            optical_list_flash_center_y_v.pop_back();
            optical_list_flash_center_z_v.pop_back();
            optical_list_flash_time_v.pop_back();

            optical_list_pe.push_back(fOpFlashPE);
            optical_list_flash_center_y.push_back(fOpFlashCenterY);
            optical_list_flash_center_z.push_back(fOpFlashCenterZ);
            optical_list_flash_time.push_back(fOpFlashTime);

        }

        last_event = current_event;
        last_run   = current_run;

        optical_list_pe_v.push_back(optical_list_pe);
        optical_list_flash_center_y_v.push_back(optical_list_flash_center_y);
        optical_list_flash_center_z_v.push_back(optical_list_flash_center_z);
        optical_list_flash_time_v.push_back(optical_list_flash_time);

    }
    
    std::cout << "Resized Optical List Size: " << optical_list_pe_v.size() << std::endl;
    
    // Largest Flash Vector
    std::vector<double> largest_flash_v;
    largest_flash_v.resize(4, 0);
    
    // ----------------------
    //      Event loop
    // ----------------------
    for(int i = 0; i < optical_list_pe_v.size(); i++) {
        
        bool in_time                 = false;
        bool got_in_time             = false;
        bool sufficient_flash         = false;
        bool got_sufficient_flash     = false;
        
        double largest_flash = 0.;
        double largest_center_y = 0;
        double largest_center_z = 0;
        double largest_flash_time = 0;
        
        // Cut Variables defined in main.h
        int flash_pe_threshold = 50;

        // Loop through all flashes in event and find largest
        for(int j = 0; j < optical_list_pe_v.at(i).size(); j++) {
            
            auto const opt_time         = optical_list_flash_time_v.at(i).at(j) ; // shift due to MC and offset
            auto const opt_pe           = optical_list_pe_v.at(i).at(j);
            const double opt_center_y   = optical_list_flash_center_y_v.at(i).at(j);
            const double opt_center_z   = optical_list_flash_center_z_v.at(i).at(j);
            const double opt_flash_time = optical_list_flash_time_v.at(i).at(j); // shift due to MC and offset
            
            // See if flash was in time
            in_time = flash_in_time(opt_time, flash_time_start, flash_time_end); 
            if(in_time == true) got_in_time = true; 
            
            // See if flash meets the threshold requirements
            sufficient_flash = flash_pe(opt_pe, flash_pe_threshold); 
            if(sufficient_flash == true) got_sufficient_flash = true; 
            
            //Flash is both in time and over PE threshold
            if(in_time == true && sufficient_flash == true) {
                
                // Find the largest flash in this event
                if(opt_pe > largest_flash) {
                    largest_flash      = opt_pe;
                    largest_center_y   = opt_center_y;
                    largest_center_z   = opt_center_z;
                    largest_flash_time = opt_flash_time;
                }
            }
        }
        
        largest_flash_v.at(0) = largest_center_y;
        largest_flash_v.at(1) = largest_center_z;
        largest_flash_v.at(2) = largest_flash_time;
        largest_flash_v.at(3) = largest_flash;
        
        largest_flash_v_v.push_back(largest_flash_v);
        
    }

    return largest_flash_v_v;
}
//***************************************************************************
//***************************************************************************
double variation_output_bkg::Flash_TPCObj_vtx_Dist(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z) {
    const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
    return distance;
}
//***************************************************************************
//***************************************************************************
std::pair<std::string, int> variation_output_bkg::TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool true_in_tpc, bool has_pi0) {
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
    return std::make_pair("non_match", 0);
}
//*********************** Selection Cuts ************************************
//***************************************************************************
//IN FV
bool variation_output_bkg::in_fv(double x, double y, double z, std::vector<double> fv_boundary_v) {
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
}
//***************************************************************************
//***************************************************************************
// Flash in time
bool variation_output_bkg::flash_in_time(double flash_time, double flash_start, double flash_end) {
    if(flash_time >= flash_start && flash_time <= flash_end) return true; // Pass in time
    else return false; // Fail in time
}
//***************************************************************************
//***************************************************************************
// Flash PE
bool variation_output_bkg::flash_pe(int flash_pe, int flash_pe_threshold) {
    if (flash_pe >= flash_pe_threshold) return true; // Pass PE Thresh
    else return false; // Fail PE Thresh
}
//***************************************************************************
//***************************************************************************
// Get the vector of flashes with true/false on whether it passed the selection or not
void variation_output_bkg::FlashinTime_FlashPE(TFile* f, double flash_start_time, double flash_end_time, std::vector<bool> &flash_cuts_pass_vec, TString mode ){
    // Get Optical Information from file
    TTree* optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

    std::vector<std::vector<double>> largest_flash_v_v;

    // ----------------------
    //    Optical Info
    // ----------------------
    int fRun = 0;
    int fEvent = 0;
    int fOpFlashPE = 0;
    double fOpFlashTime = 0;
    double fOpFlashWidthY = 0;
    double fOpFlashWidthZ = 0;
    double fOpFlashCenterY = 0;
    double fOpFlashCenterZ = 0;

    optical_tree->SetBranchAddress("Run",              &fRun           ); // Run number of Flash
    optical_tree->SetBranchAddress("Event",            &fEvent         ); // Event number of Flash
    optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE     ); // PE of Flash
    optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime   ); // Time of flash
    optical_tree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY );
    optical_tree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ );
    optical_tree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY); // Flash Y center
    optical_tree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ); // Flash Z Center

    // Num events in the optical tree
    const int optical_entries = optical_tree->GetEntries();
    std::cout << "Total Optical Entries: " << optical_entries << std::endl;

    int current_event = 0;
    int current_run = 0;
    int last_event = 0;
    int last_run = 0;

    // Contains the entry number for a given OpFlash per event
    std::vector<int>                    optical_list_pe;
    std::vector<std::vector<int> >        optical_list_pe_v;
    
    std::vector<double>                    optical_list_flash_center_y; 
    std::vector<std::vector<double> >    optical_list_flash_center_y_v;
    
    std::vector<double>                    optical_list_flash_center_z; 
    std::vector<std::vector<double> >    optical_list_flash_center_z_v;
    
    std::vector<double>                    optical_list_flash_time;
    std::vector<std::vector<double> >    optical_list_flash_time_v;
    
    // Loop over the optical entries to get the largest flash vector
    
    // ----------------------
    // Resize the optical enties to be the same sizd as number of Events (TPC Obj)
    // ----------------------
    
    for(int i = 0; i < optical_entries; i++) {
        
        // Get the Optical entry
        optical_tree->GetEntry(i);

        current_run        = fRun;
        current_event     = fEvent;

        // New event
        if(current_event != last_event) {
            optical_list_pe.clear();
            optical_list_flash_center_y.clear();
            optical_list_flash_center_z.clear();
            optical_list_flash_time.clear();

            optical_list_pe.push_back(fOpFlashPE);
            optical_list_flash_center_y.push_back(fOpFlashCenterY);
            optical_list_flash_center_z.push_back(fOpFlashCenterZ);
            optical_list_flash_time.push_back(fOpFlashTime);

        }
        // Same event
        if(current_event == last_event && current_run == last_run) {
            optical_list_pe_v.pop_back();
            optical_list_flash_center_y_v.pop_back();
            optical_list_flash_center_z_v.pop_back();
            optical_list_flash_time_v.pop_back();

            optical_list_pe.push_back(fOpFlashPE);
            optical_list_flash_center_y.push_back(fOpFlashCenterY);
            optical_list_flash_center_z.push_back(fOpFlashCenterZ);
            optical_list_flash_time.push_back(fOpFlashTime);

        }

        last_event = current_event;
        last_run   = current_run;

        optical_list_pe_v.push_back(optical_list_pe);
        optical_list_flash_center_y_v.push_back(optical_list_flash_center_y);
        optical_list_flash_center_z_v.push_back(optical_list_flash_center_z);
        optical_list_flash_time_v.push_back(optical_list_flash_time);

    }
    
    std::cout << "Resized Optical List Size: " << optical_list_pe_v.size() << std::endl;
    
    flash_cuts_pass_vec.resize(optical_list_pe_v.size());

    // Loop over the optical list
    for(int i = 0; i < optical_list_pe_v.size(); i++) {
        
        bool in_time = false;
        bool sufficient_flash = false;
        
        auto const opt_time_v = optical_list_flash_time_v.at(i);
        auto const opt_pe_v   = optical_list_pe_v.at(i);
        
        // Loop over the optical list vec
        for (int j = 0; j < optical_list_pe_v.at(i).size(); j++) {
            
            auto opt_time =  opt_time_v.at(j);
            if (mode == "numi") opt_time = opt_time_v.at(j) + 1.0;
            else opt_time = opt_time_v.at(j);
            
            auto const opt_pe = opt_pe_v.at(j);
            
            in_time          = flash_in_time(opt_time, flash_start_time, flash_end_time);
            sufficient_flash = flash_pe(opt_pe, flash_pe_threshold);
            
            // Flash is both in time and over PE threshold
            if(in_time == true && sufficient_flash == true){
                flash_cuts_pass_vec.at(i) = true;
                break; // once pased we are done, so dont loop any more otherwise we may overwrite this

            }
        }

        if (in_time == false && sufficient_flash == false) {
            flash_cuts_pass_vec.at(i) = false;
        }

    }

}
//***************************************************************************
//***************************************************************************
// Check for valid reconstructed nue
bool variation_output_bkg::HasNue(xsecAna::TPCObjectContainer tpc_obj, const int n_pfp ) {

    bool has_nue = false;
    bool has_valid_shower = false;

    for(int j = 0; j < n_pfp; j++) {
            auto const part     = tpc_obj.GetParticle(j);
            const int  pfp_pdg  = part.PFParticlePdgCode();
            const int  pfp_hits = part.NumPFPHits();
            
            if(pfp_pdg == 11 && pfp_hits > 0) has_valid_shower = true; 
            
            if(pfp_pdg == 12) has_nue = true; 
    }

    
    if(has_nue == true && has_valid_shower == true)
        return true; 
    else 
        return false;
    
}
//***************************************************************************
//***************************************************************************
// Flash reco Vertex Distance 
bool variation_output_bkg::opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance) {
    const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
    
    if(distance <= tolerance) return true;
    return false;
}
bool variation_output_bkg::flashRecoVtxDist(std::vector< double > largest_flash_v, double tolerance, const double tpc_vtx_x, const double tpc_vtx_y, const double tpc_vtx_z) {
    
    bool is_close;
    
    //flash is upstream
    if(tpc_vtx_z < largest_flash_v.at(1)) 
        is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), tolerance);
    
    //flash is downstream
    if(tpc_vtx_z >= largest_flash_v.at(1)) 
        is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), (tolerance - 20));
    
    if (is_close == true )
        return true;
    else
        return false;
    
}
//***************************************************************************
//***************************************************************************
// Vertex to shower distance
bool variation_output_bkg::VtxNuDistance(xsecAna::TPCObjectContainer tpc_obj,int pfp_pdg_type , double tolerance){
    
    const int n_pfp = tpc_obj.NumPFParticles();
    const double tpc_vtx_x = tpc_obj.pfpVtxX();
    const double tpc_vtx_y = tpc_obj.pfpVtxY();
    const double tpc_vtx_z = tpc_obj.pfpVtxZ();

    const int n_tracks = tpc_obj.NPfpTracks();
    if (n_tracks == 0 && pfp_pdg_type == 13 ) return true; 

    for (int j = 0; j < n_pfp; j++) {

        auto const part   = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();

        if (pfp_pdg == pfp_pdg_type) {

            const double pfp_vtx_x = part.pfpVtxX();
            const double pfp_vtx_y = part.pfpVtxY();
            const double pfp_vtx_z = part.pfpVtxZ();

            const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );

            if (distance <= tolerance) return true; 
        }

    }
    return false;

}
//***************************************************************************
//***************************************************************************
// Hit thresholds all planes
bool variation_output_bkg::HitThreshold(xsecAna::TPCObjectContainer tpc_obj, double threshold, bool useCollection){

    const int n_pfp = tpc_obj.NumPFParticles();

    for (int j = 0; j < n_pfp; j++) {

        auto const pfp_obj = tpc_obj.GetParticle(j);
        int  num_pfp_hits  = pfp_obj.NumPFPHits();
        const int  pfp_pdg = pfp_obj.PFParticlePdgCode();

        if (useCollection) num_pfp_hits = pfp_obj.NumPFPHitsW(); // Collection plane hits

        if (pfp_pdg == 11 && num_pfp_hits >= threshold) return true;
    }
    
    return false;
    

}
//***************************************************************************
//***************************************************************************
// Leading shower open angle
bool variation_output_bkg::OpenAngleCut(xsecAna::TPCObjectContainer tpc_obj, const std::vector<double> tolerance_open_angle){

    const int n_pfp = tpc_obj.NumPFParticles();

    int leading_index = 0;
    int leading_hits  = 0;
    
    for (int j = 0; j < n_pfp; j++) {
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const int n_pfp_hits = part.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
    }
    
    auto const leading_shower       = tpc_obj.GetParticle(leading_index);
    const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);

    if (leading_open_angle <= tolerance_open_angle.at(1) && leading_open_angle >= tolerance_open_angle.at(0))
        return true;
    else 
        return false;

}
//***************************************************************************
//***************************************************************************
// leading shower dEdx cut
bool variation_output_bkg::dEdxCut( xsecAna::TPCObjectContainer tpc_obj, const double tolerance_dedx_min, const double tolerance_dedx_max){

    const int n_pfp = tpc_obj.NumPFParticles();
    int leading_index = 0;
    int leading_hits  = 0;

    for (int j = 0; j < n_pfp; j++) {
        
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const int n_pfp_hits = part.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
    } 
    
    auto const leading_shower = tpc_obj.GetParticle(leading_index);
    double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
    leading_dedx = leading_dedx * (196.979 /242.72);

    if (leading_dedx <= tolerance_dedx_max && leading_dedx >= tolerance_dedx_min) return true;
     
    return false;

}
//***************************************************************************
//***************************************************************************
// Secondar Showers Cut
bool variation_output_bkg::SecondaryShowersDistCut(xsecAna::TPCObjectContainer tpc_obj, const double dist_tolerance){

    const int n_pfp = tpc_obj.NumPFParticles();
    const int n_pfp_showers = tpc_obj.NPfpShowers();
    
    // This cut does not target events with fewer than 2 showers
    if (n_pfp_showers <= 1) return true; 
    
    const double tpco_vtx_x = tpc_obj.pfpVtxX();
    const double tpco_vtx_y = tpc_obj.pfpVtxY();
    const double tpco_vtx_z = tpc_obj.pfpVtxZ();
    int leading_index = 0;
    int leading_hits  = 0;
    
    for (int j = 0; j < n_pfp; j++) {
        
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const int n_pfp_hits = part.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
    }
    
    for (int j = 0; j < n_pfp; j++) {
        
        if (j == leading_index) continue; // We assume leading shower == electron shower
        
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const double pfp_vtx_x = part.pfpVtxX();
        const double pfp_vtx_y = part.pfpVtxY();
        const double pfp_vtx_z = part.pfpVtxZ();
        
        const double distance = sqrt(pow((pfp_vtx_x - tpco_vtx_x),2) + pow((pfp_vtx_y - tpco_vtx_y),2) + pow((pfp_vtx_z - tpco_vtx_z),2));
        
        if (pfp_pdg == 11) {
            if (distance > dist_tolerance) return false;
            // if (distance <= dist_tolerance) return true;
                
        }
    }
    
    return true;
    
}
//***************************************************************************
//***************************************************************************
// Hits per length cut
bool variation_output_bkg::HitLengthRatioCut(const double pfp_hits_length_tolerance, xsecAna::TPCObjectContainer tpc_obj){
    
    const int n_pfp = tpc_obj.NumPFParticles();
    const int n_pfp_showers = tpc_obj.NPfpShowers();
    int leading_index = 0;
    int leading_hits  = 0;
    
    for (int j = 0; j < n_pfp; j++) {
        
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const int n_pfp_hits = part.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
    }
    
    auto const leading_shower = tpc_obj.GetParticle(leading_index);
    const int pfp_pdg = leading_shower.PFParticlePdgCode();
    const double pfp_hits = leading_shower.NumPFPHits();
    const double pfp_length = leading_shower.pfpLength();
    const double pfp_hits_length_ratio = (pfp_hits / pfp_length);

    if (pfp_pdg == 11 && pfp_hits_length_ratio > pfp_hits_length_tolerance ) return true;

    return false;
    
}
//***************************************************************************
//***************************************************************************
// Longest Track Leading Shower Cut
bool variation_output_bkg::LongestTrackLeadingShowerCut(const double ratio_tolerance, xsecAna::TPCObjectContainer tpc_obj){

    const int n_pfp = tpc_obj.NumPFParticles();
    const int n_pfp_tracks = tpc_obj.NPfpTracks();
    
    if (n_pfp_tracks == 0) return true;
    
    int leading_index = 0;
    int leading_hits  = 0;
    double longest_track = 0;
    
    for (int j = 0; j < n_pfp; j++) {
        
        auto const pfp = tpc_obj.GetParticle(j);
        const int pfp_pdg = pfp.PFParticlePdgCode();
        const int n_pfp_hits = pfp.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
        
        if(pfp_pdg == 13) {
            
            const double trk_length = pfp.pfpLength();
            
            if (trk_length > longest_track) longest_track = trk_length;
            
        }
    
    } //end loop pfparticles
    
    auto const leading_shower = tpc_obj.GetParticle(leading_index);
    const double leading_shower_length = leading_shower.pfpLength();
    const double longest_track_leading_shower_ratio = longest_track / leading_shower_length;

    //if the ratio is too large:
    if (longest_track_leading_shower_ratio > ratio_tolerance) return false;

    return true;

}
//***************************************************************************
//***************************************************************************
// Contained Tracks Cut
bool variation_output_bkg::IsContained(std::vector<double> track_start, std::vector<double> track_end, std::vector<double> fv_boundary_v) {
    
    if(in_fv(track_start.at(0), track_start.at(1), track_start.at(2), fv_boundary_v) == true
       && in_fv(track_end.at(0), track_end.at(1), track_end.at(2), fv_boundary_v) == true) {
        return true;
    }
    else 
        return false;
}
bool variation_output_bkg::ContainedTracksCut(std::vector<double> fv_boundary_v, xsecAna::TPCObjectContainer tpc_obj){

    const int n_pfp = tpc_obj.NumPFParticles();
    const int n_pfp_tracks = tpc_obj.NPfpTracks();
    
    // This is normally enabled, but due to test cut below it is off
    if (n_pfp_tracks == 0) return true;

    for (int j = 0; j < n_pfp; j++) {
        
        auto const pfp = tpc_obj.GetParticle(j);
        const int pfp_pdg = pfp.PFParticlePdgCode();
        
        if (pfp_pdg == 13) {
            
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
            if(is_contained == false) return false;
        
        } // end is track
    
    } // end loop pfprticles
    return true;    
}
//***************************************************************************
//***************************************************************************
std::string variation_output_bkg::Background_Classifier(int mc_pdg, std::string tpc_obj_classification){
    std::string classification = "NULL";

    // Due to low statistics we remove this catogory and absorb it into "other"
    // if (mc_pdg == 11)
    //     classification = "bkg_e";                                // out of fv, mixed, cosmic (see two nc_pi0 ??)
    
    if ( (mc_pdg == 22 ) && tpc_obj_classification != "cosmic"){ // photons from pi 0
        classification = "pi0_gamma";
        classification = "other_bkg"; // put this here to absorb these backgrounds into this category
    }
    
    else if (tpc_obj_classification == "cosmic"){
        classification = "cosmic";
    }
    // else if (tpc_obj_classification == "nue_cc_out_fv" || tpc_obj_classification == "nue_cc_mixed" || tpc_obj_classification == "nc") {
    //     classification = "e_bkg";
    // }
    else{
        classification = "other_bkg";
        
        // std::cout << "bkg class:\t" << tpc_obj_classification << std::endl;
    }

    return classification; 
}
//***************************************************************************
//***************************************************************************
double variation_output_bkg::WrapPhi(double phi){
    double phi_wrapped{0};


    // >90 to 180
    if      (phi > 90   && phi <= 180) phi_wrapped = 180 - phi;
    
    // <0 to -90
    else if (phi >= -90  && phi < 0)   phi_wrapped = phi * -1;
    
    // <-90 to -180
    else if (phi >= -180 && phi < -90) phi_wrapped = 180 - (phi * -1);
    
    else phi_wrapped = phi;

    return phi_wrapped;
}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::GenerateWeightHistograms(){

    std::cout << "Generating the weight histograms" << std::endl;

    // First open the BNB histogram file
    TFile *fweight = new TFile("plots/variation_weights.root","UPDATE");
    TFile *fnumi = TFile::Open("plots/variation_out_numi_bkg.root"); // Get the numi file 

    // Loop over all folders (loop over variations)
    std::string histname_ratio;
    std::vector<std::string> histnames = { "h_total_hits",             "h_ldg_shwr_hits",             "h_ldg_shwr_hits_WPlane",
                                          "h_ldg_shwr_Open_Angle",    "h_ldg_shwr_dEdx_WPlane",      "h_ldg_shwr_HitPerLen",
                                          "h_ldg_shwr_CTheta",
                                          "h_long_Track_ldg_shwr",    "h_tpc_obj_vtx_x",             "h_tpc_obj_vtx_y",        "h_tpc_obj_vtx_z",
                                          "h_n_pfp",                  "h_n_pfp_50Hits",              "h_n_tracks",             "h_n_tracks_50Hits",  "h_n_showers",
                                          "h_n_showers_50Hits",       "h_track_phi",                 "h_shower_phi",           "h_largest_flash_y",  "h_largest_flash_z",
                                          "h_largest_flash_time",     "h_largest_flash_pe",          "h_Flash_TPCObj_Dist",
                                          "h_shower_Nu_vtx_Dist",     "h_track_Nu_vtx_Dist",
                                          "h_ldg_shwr_Phi",           "h_ldg_shwr_Phi_wrapped",          "h_ldg_shwr_Theta",
                                          "h_shower_phi_pi0",         "h_shower_phi_bkg_cosmic",         "h_shower_phi_other",
                                          "h_shower_phi_pi0_wrapped", "h_shower_phi_bkg_cosmic_wrapped", "h_shower_phi_other_wrapped",
                                          "h_shower_E_pi0",           "h_shower_E_bkg_cosmic",           "h_shower_E_other", "h_shower_E", "h_selected",
                                          "h_shower_Theta_pi0",       "h_shower_Theta_bkg_cosmic",       "h_shower_Theta_other"
                                        //   "h_ldg_shwr_Phi_unselected",           "h_ldg_shwr_Phi_wrapped_unselected",          "h_ldg_shwr_Theta_unselected",
                                        //   "h_shower_phi_pi0_unselected",         "h_shower_phi_bkg_cosmic_unselected",         "h_shower_phi_other_unselected",
                                        //   "h_shower_phi_pi0_wrapped_unselected", "h_shower_phi_bkg_cosmic_wrapped_unselected", "h_shower_phi_other_wrapped_unselected",
                                        //   "h_shower_E_pi0_unselected",           "h_shower_E_bkg_cosmic_unselected",           "h_shower_E_other_unselected", "h_shower_E_unselected",
                                        //   "h_shower_Theta_pi0_unselected",       "h_shower_Theta_bkg_cosmic_unselected",       "h_shower_Theta_other_unselected"
                                          };
    
    // Loop over the histograms
    for (int j=0; j < histnames.size(); j++){
        TH1D* hist_NuMI;
        TH1D* hist_BNB;

        fnumi->cd();
        
        char name_NuMI[500];
        snprintf(name_NuMI, 500, "NuMICV/%s", histnames[j].c_str() );
        hist_NuMI = (TH1D*)fnumi->Get(name_NuMI);
        if (hist_NuMI == NULL ) std::cout << "ERROR: Can't get Histogram!   " << name_NuMI <<  std::endl;

        char name_BNB[500];
        snprintf(name_BNB, 500, "BNBCV/%s", histnames[j].c_str() );
        hist_BNB = (TH1D*)f_var_out->Get(name_BNB);
        if (hist_BNB == NULL ) std::cout << "ERROR: Can't get CV Histogram!    " << name_BNB << std::endl;

        histname_ratio = histnames[j] + "_ratio";
        
        TH1D* hist_divide = (TH1D*) hist_NuMI->Clone(Form("%s",histname_ratio.c_str()));
        hist_divide->Sumw2();
        
        hist_divide->Divide(hist_BNB);
        hist_divide->SetOption("E");

        fweight->cd();
        hist_divide->Write("",TObject::kOverwrite);
        f_var_out->cd();

    }

    // Now lets save the 2D histograms for weighting
    histnames.clear();
    histnames = { "h_ThetaBkg_Phi_wrapped", "h_ThetaBkg_pi0_Phi_wrapped", "h_ThetaBkg_cosmic_Phi_wrapped",
                  "h_ThetaBkg_other_Phi_wrapped", "h_ThetaBkg_Phi", "h_ThetaBkg_pi0_Phi", "h_ThetaBkg_cosmic_Phi", "h_ThetaBkg_other_Phi"
                //   "h_ThetaBkg_Phi_wrapped_unselected", "h_ThetaBkg_pi0_Phi_wrapped_unselected", "h_ThetaBkg_cosmic_Phi_wrapped_unselected",
                //   "h_ThetaBkg_other_Phi_wrapped_unselected", "h_ThetaBkg_Phi_unselected", "h_ThetaBkg_pi0_Phi_unselected", "h_ThetaBkg_cosmic_Phi_unselected", "h_ThetaBkg_other_Phi_unselected"
                   };

    for (int j=0; j < histnames.size(); j++){
        TH2D* hist_NuMI;
        TH2D* hist_BNB;

        fnumi->cd();
        
        char name_NuMI[500];
        snprintf(name_NuMI, 500, "NuMICV/%s", histnames[j].c_str() );
        hist_NuMI = (TH2D*)fnumi->Get(name_NuMI);
        if (hist_NuMI == NULL ) std::cout << "ERROR: Can't get Histogram!   " << name_NuMI <<  std::endl;

        char name_BNB[500];
        snprintf(name_BNB, 500, "BNBCV/%s", histnames[j].c_str() );
        hist_BNB = (TH2D*)f_var_out->Get(name_BNB);
        if (hist_BNB == NULL ) std::cout << "ERROR: Can't get CV Histogram!    " << name_BNB << std::endl;

        histname_ratio = histnames[j] + "_ratio";
        
        TH2D* hist_divide = (TH2D*) hist_NuMI->Clone(Form("%s",histname_ratio.c_str()));
        hist_divide->Sumw2();
        
        hist_divide->Divide(hist_BNB);
        hist_divide->SetOption("colz");

        fweight->cd();
        hist_divide->Write("",TObject::kOverwrite);
        f_var_out->cd();

    }

    fweight->Close();
}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::WeightBNBVar(xsecAna::TPCObjectContainer tpc_obj ,bool bool_sig, const int leading_shower_index, std::pair<std::string, int> tpc_classification, std::string dirname){

    // Grab the variation folders in the file
    std::vector<std::string> variations = {
                                        "BNBwithDIC",
                                        "BNBdataSCE",
                                        "BNBdeadSaturatedChannels",
                                        "BNBLArG4BugFix",
                                        "BNBDLup",
                                        "BNBDLdown",
                                        "BNBDTup",
                                        "BNBDTdown",
                                        "BNBnoiseAmpUp",
                                        "BNBnoiseAmpDown",
                                        "BNBaltDeadChannels",
                                        "BNBstretchResp",
                                        "BNBsqueezeResp",
                                        "BNBupPEnoise",
                                        "BNBEnhancedTPCVis",
                                        "BNBBirksRecomb",
                                        "BNBdownPEnoise"};

    int p = -1; // Index of variation

    // Get the index of variation to weight
    for (unsigned int j=0; j< variations.size(); j++){
        if (dirname == variations.at(j)) p = j;
    }

    // No variation was matched so we just return
    if (p == -1) return;


    n_pfp = tpc_obj.NumPFParticles();

    // Now we loop over the particles in the tpc object and weight these
    for (int j = 0; j < n_pfp ; j++){

        auto const pfp_obj = tpc_obj.GetParticle(j);

        const double mc_Theta  = pfp_obj.mcTheta();
        mc_Phi    = pfp_obj.mcPhi();
        mc_Energy = pfp_obj.mcEnergy();
        const double mc_pdg    = pfp_obj.MCPdgCode();
        
        // Background events
        if (!bool_sig) {

            const double shower_phi = atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;

            //  ------------ Leading shower ------------
            if (j == leading_shower_index){

                bkg_class = Background_Classifier(mc_pdg, tpc_classification.first);

                double mc_phi_wrapped = WrapPhi(mc_Phi);
            
                std::string histname = "h_shower";

                TH1D *hist; // temp hist for getting weighted histograms
                char name[500];
                int xbin{0};
                double weight{0};

                // Fill the phi distribtuons for the bkgs
                if (bkg_class == "pi0_gamma") {
                    std::string histname_phi_pi0 = histname + "_phi_pi0_ratio";
                    snprintf(name, 500, "%s", histname_phi_pi0.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name<<  std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_Phi); // Get the bin at Phi
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_phi_pi0_w).at(p)->Fill(mc_Phi, weight);

                    std::string histname_phi_pi0_wrapped = histname + "_phi_pi0_wrapped_ratio";
                    snprintf(name, 500, "%s", histname_phi_pi0_wrapped.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_phi_wrapped); // Get the bin at Phi
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_phi_pi0_wrapped_w).at(p)->Fill(mc_phi_wrapped, weight);


                    std::string histname_E_pi0 = histname + "_E_pi0_ratio";
                    snprintf(name, 500, "%s", histname_E_pi0.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_Energy); // Get the bin at Phi
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_E_pi0_w).at(p)->Fill(mc_Energy, weight);

                    std::string histname_Theta_pi0 = histname + "_Theta_pi0_ratio";
                    snprintf(name, 500, "%s", histname_Theta_pi0.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name<<  std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_Theta); // Get the bin at Theta
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_Theta_pi0_w).at(p)->Fill(mc_Theta, weight);
                    
                }
                
                if (bkg_class == "cosmic") {
                    std::string histname_phi_cosmic = histname + "_phi_bkg_cosmic_ratio";
                    snprintf(name, 500, "%s", histname_phi_cosmic.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: "<< name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_Phi); // Get the bin at Phi
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_phi_bkg_cosmic_w).at(p)->Fill(mc_Phi, weight);


                    std::string histname_phi_cosmic_wrapped = histname + "_phi_bkg_cosmic_wrapped_ratio";
                    snprintf(name, 500, "%s", histname_phi_cosmic_wrapped.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_phi_wrapped); // Get the bin at Phi
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_phi_bkg_cosmic_wrapped_w).at(p)->Fill(mc_phi_wrapped, weight);


                    std::string histname_E_cosmic = histname + "_E_bkg_cosmic_ratio";
                    snprintf(name, 500, "%s", histname_E_cosmic.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_Energy); // Get the bin at Phi
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_E_bkg_cosmic_w).at(p)->Fill(mc_Energy, weight);

                    std::string histname_Theta_cosmic = histname + "_Theta_bkg_cosmic_ratio";
                    snprintf(name, 500, "%s", histname_Theta_cosmic.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: "<< name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_Theta); // Get the bin at Theta
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_Theta_bkg_cosmic_w).at(p)->Fill(mc_Theta, weight);

                }
                if (bkg_class == "other_bkg") {
                    std::string histname_phi_other = histname + "_phi_other_ratio";
                    snprintf(name, 500, "%s", histname_phi_other.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_Phi); // Get the bin at Phi
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_phi_other_w).at(p)->Fill(mc_Phi, weight);


                    std::string histname_phi_other_wrapped = histname + "_phi_other_wrapped_ratio";
                    snprintf(name, 500, "%s", histname_phi_other_wrapped.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_phi_wrapped); // Get the bin at Phi
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_phi_other_wrapped_w).at(p)->Fill(mc_phi_wrapped, weight);

                    std::string histname_E_other = histname + "_E_other_ratio";
                    snprintf(name, 500, "%s", histname_E_other.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_Energy); // Get the bin at Phi
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_E_other_w).at(p)->Fill(mc_Energy, weight);

                    std::string histname_Theta_other = histname + "_Theta_other_ratio";
                    snprintf(name, 500, "%s", histname_Theta_other.c_str() );
                    hist = (TH1D*)fweight->Get(name); 
                    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                    xbin = hist->GetXaxis()->FindBin(mc_Theta); // Get the bin at Theta
                    weight = hist->GetBinContent(xbin); // Get the weight at that bin

                    // Here we would fill the histogram
                    if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                    TH1D_hist_weighted.at(kshower_Theta_other_w).at(p)->Fill(mc_Theta, weight);

                } // End if background 


                // Phi
                std::string histname_phi ="h_ldg_shwr_Phi_ratio";
                snprintf(name, 500, "%s", histname_phi.c_str() );
                hist = (TH1D*)fweight->Get(name); 
                if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                xbin = hist->GetXaxis()->FindBin(mc_Phi); // Get the bin at Phi
                weight = hist->GetBinContent(xbin); // Get the weight at that bin

                // Here we would fill the histogram
                if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                TH1D_hist_weighted.at(kldg_shwr_Phi_w).at(p)->Fill(mc_Phi, weight);

                // Phi Wrapped
                std::string histname_phi_wrapped = "h_ldg_shwr_Phi_wrapped_ratio";
                snprintf(name, 500, "%s", histname_phi_wrapped.c_str() );
                hist = (TH1D*)fweight->Get(name); 
                if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                xbin = hist->GetXaxis()->FindBin(mc_phi_wrapped); // Get the bin at Phi
                weight = hist->GetBinContent(xbin); // Get the weight at that bin

                // Here we would fill the histogram
                if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                TH1D_hist_weighted.at(kldg_shwr_Phi_wrapped_w).at(p)->Fill(mc_phi_wrapped, weight);                

                // Energy
                std::string histname_E = histname + "_E_ratio";
                snprintf(name, 500, "%s", histname_E.c_str() );
                hist = (TH1D*)fweight->Get(name); 
                if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                xbin = hist->GetXaxis()->FindBin(mc_Energy); // Get the bin at Phi
                weight = hist->GetBinContent(xbin); // Get the weight at that bin

                // Here we would fill the histogram
                if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                TH1D_hist_weighted.at(kshower_E_w).at(p)->Fill(mc_Energy, weight);

                // Theta
                std::string histname_Theta = "h_ldg_shwr_Theta_ratio";
                snprintf(name, 500, "%s", histname_Theta.c_str() );
                hist = (TH1D*)fweight->Get(name); 
                if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                xbin = hist->GetXaxis()->FindBin(mc_Theta); // Get the bin at Phi
                weight = hist->GetBinContent(xbin); // Get the weight at that bin

                // Here we would fill the histogram
                if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                TH1D_hist_weighted.at(kldg_shwr_Theta_w).at(p)->Fill(mc_Theta, weight);
                
                // Selected
                std::string histname_Selected = "h_selected_ratio";
                snprintf(name, 500, "%s", histname_Selected.c_str() );
                hist = (TH1D*)fweight->Get(name); 
                if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
                xbin = hist->GetXaxis()->FindBin(0.5); // Get the bin at Phi
                weight = hist->GetBinContent(xbin); // Get the weight at that bin

                // Here we would fill the histogram
                if (weight == 0) weight = 1; // We shouldt fill with empty weights because of lack of stats
                TH1D_hist_weighted.at(kselected_w).at(p)->Fill(1, weight);
                
            }
        }
    }
}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::CompareWeightedHistograms(){

    system("if [ ! -d \"plots/weighted_numi/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/weighted_numi/; fi");
    system("if [ ! -d \"plots_png/weighted_numi/\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots_png/weighted_numi/; fi");

    std::vector<std::string> variations= {"NuMILArG4BugFix", "NuMIwithDIC", "NuMIdeadSaturatedChannels"};
    std::vector<std::string> variations_BNB= {"BNBLArG4BugFix", "BNBwithDIC", "BNBdeadSaturatedChannels"};

    std::vector<std::string> histnames = {
                                          "h_ldg_shwr_Phi",           "h_ldg_shwr_Phi_wrapped",          "h_ldg_shwr_Theta",
                                          "h_shower_phi_pi0",         "h_shower_phi_bkg_cosmic",         "h_shower_phi_other",
                                          "h_shower_phi_pi0_wrapped", "h_shower_phi_bkg_cosmic_wrapped", "h_shower_phi_other_wrapped",
                                          "h_shower_E_pi0",           "h_shower_E_bkg_cosmic",           "h_shower_E_other", "h_shower_E", "h_selected",
                                          "h_shower_Theta_pi0",       "h_shower_Theta_bkg_cosmic",       "h_shower_Theta_other",};

    char name[500];

    for (unsigned int i = 0; i<variations.size(); i++){

        for (unsigned int k=0; k< histnames.size(); k++){

            TCanvas* c = new TCanvas();
            TLegend* legend = new TLegend(0.60, 0.80, 0.90, 0.89);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);

            TH1D *hist_NuMI; // Thea actual numi variation
            std::string histname_NuMI = variations.at(i) + "/" + histnames.at(k);
            snprintf(name, 500, "%s", histname_NuMI.c_str() );
            hist_NuMI = (TH1D*)f_var_out->Get(name); 
            if (hist_NuMI == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;

            TH1D *hist_NuMI_Weighted;  // The NuMI CV weighted by the BNB
            std::string histname_NuMI_Weighted = "BNBWeighted_Variations/"+ histnames.at(k)+"_"+variations_BNB.at(i)+"_w";
            snprintf(name, 500, "%s", histname_NuMI_Weighted.c_str() );
            hist_NuMI_Weighted = (TH1D*)f_var_out->Get(name);
            if (hist_NuMI_Weighted == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;

            CompareWeightedDrawSpecs(hist_NuMI_Weighted, "BNBWeighted_Variations", variations.at(i), legend, histnames.at(k));
            CompareWeightedDrawSpecs(hist_NuMI, "", variations.at(i), legend, histnames.at(k)); 
            
            legend->Draw();

            char Canvas_name[1000];
            snprintf(Canvas_name, 1000, "plots/weighted_numi/%s_%s.pdf", histnames.at(k).c_str(), variations.at(i).c_str()); 
            c->Print(Canvas_name);
            snprintf(Canvas_name, 1000, "plots_png/weighted_numi/%s_%s.png",histnames.at(k).c_str(), variations.at(i).c_str() ); 
            c->Print(Canvas_name);
            c->Close();

        }

    }

}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::CompareWeightedDrawSpecs(TH1D* hist, std::string weighted_str,  std::string variation, TLegend* legend, std::string histname){


    if (histname == "h_ldg_shwr_Phi"){
        hist->SetTitle(";Leading Shower #phi [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,250);
    }
    else if (histname == "h_ldg_shwr_Theta"){
        hist->SetTitle(";Leading Shower #theta [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,200);

    }
    else if (histname == "h_ldg_shwr_Phi_wrapped"){
        hist->SetTitle(";Leading Shower #phi (wrapped) [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,900);

    }
    else if (histname == "h_selected"){
        hist->SetTitle("; Num Background Selected;Entries");
        hist->GetXaxis()->SetLabelOffset(999);
        hist->GetXaxis()->SetLabelSize(0);
        // hist->GetYaxis()->SetRangeUser(400,700);

    }
    else if (histname ==  "h_shower_phi_pi0"){
        hist->SetTitle("; Leading Shower Phi #pi^{0} Bkg [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,50);

    }
    else if (histname ==  "h_shower_phi_bkg_cosmic"){
        hist->SetTitle("; Leading Shower Phi cosmic Bkg [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,50);
    }
    else if (histname ==  "h_shower_phi_other"){
        hist->SetTitle("; Leading Shower Phi Other Bkg [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,170);

    }
    else if (histname ==  "h_shower_phi_pi0_wrapped"){
        hist->SetTitle("; Leading Shower Phi Wrapped #pi^{0} Bkg [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,600);

    }
    else if (histname ==  "h_shower_phi_bkg_cosmic_wrapped"){
        hist->SetTitle("; Leading Shower Phi Wrapped cosmic Bkg [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,40);

    }
    else if (histname ==  "h_shower_phi_other_wrapped"){
        hist->SetTitle("; Leading Shower Phi Wrapped Other Bkg [degrees];Entries");
        hist->GetYaxis()->SetRangeUser(0,250);

    }
    else if (histname ==  "h_shower_E_pi0"){
        hist->SetTitle("; Leading Shower Energy #pi^{0} Bkg [GeV];Entries");
        // hist->GetYaxis()->SetRangeUser(0,300);

    }
    else if (histname ==  "h_shower_E_bkg_cosmic"){
        hist->SetTitle("; Leading Shower Energy cosmic Bkg [GeV];Entries");
        // hist->GetYaxis()->SetRangeUser(0,10);

    }
    else if (histname ==  "h_shower_E_other"){
        hist->SetTitle("; Leading Shower Energy Other Bkg [GeV];Entries");
        // hist->GetYaxis()->SetRangeUser(0,200);

    }
    else if (histname ==  "h_shower_E"){
        hist->SetTitle("; Leading Shower Energy Bkg [GeV];Entries");
        hist->GetYaxis()->SetRangeUser(0,1000);

    }
    else if (histname ==  "h_shower_Theta_pi0"){
        hist->SetTitle("; Leading Shower Theta #pi^{0} Bkg [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,600);

    }
    else if (histname ==  "h_shower_Theta_bkg_cosmic"){
        hist->SetTitle("; Leading Shower Theta cosmic Bkg [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,40);

    }
    else if (histname ==  "h_shower_Theta_other"){
        hist->SetTitle("; Leading Shower Theta Other Bkg [degrees];Entries");
        // hist->GetYaxis()->SetRangeUser(0,250);

    }
    else return;

    std::string draw_spec = "hist,E, same";
    hist->SetOption("hist, E");

    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);


    if (variation == "BNBCV" || variation == "NuMICV"){
        hist->SetLineColor(kBlack);
        hist->SetLineWidth(2);
        hist->SetLineStyle(1);
        legend->AddEntry(hist, "BNB CV", "l");
        hist->Draw(draw_spec.c_str());
    } 
    else if  (variation == "BNBwithDIC"         || variation == "NuMIwithDIC"){
        hist->SetLineWidth(2);
        hist->SetLineStyle(1);
        
        if (weighted_str == "BNBWeighted_Variations"){
            hist->SetLineColor(kMagenta+2);
            legend->AddEntry(hist, "Weighted BNB DIC", "l");
        } 
        else{
            hist->SetLineColor(kBlack);
            legend->AddEntry(hist, "NuMI DIC", "l");
        }
        
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBEnhancedTPCVis"  || variation == "NuMIEnhancedTPCVis" ){ 
        hist->SetLineColor(30);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Enhanced TPC Vis.", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBaltDeadChannels" || variation == "NuMIaltDeadChannels"){ 
        hist->SetLineColor(38);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Alt. Dead Chan.", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBdeadSaturatedChannels" || variation == "NuMIdeadSaturatedChannels"){
        hist->SetLineWidth(2);
        hist->SetLineStyle(1);
    
        if (weighted_str == "BNBWeighted_Variations"){
            hist->SetLineColor(28);
            legend->AddEntry(hist, "Weighted BNB Dead  Sat. Chan.", "l");
        } 
        else{
            hist->SetLineColor(kBlack);
            legend->AddEntry(hist, "NuMI Dead Sat. Chan.", "l");
        }

        hist->Draw(draw_spec.c_str());
        
    }
    else if  (variation == "BNBstretchResp"  || variation == "NuMIstretchResp" ){
        hist->SetLineColor(36);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Stretch Resp.", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBsqueezeResp"  || variation == "NuMIsqueezeResp"){
        hist->SetLineColor(1001);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Squeeze Resp.", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBupPEnoise"    || variation == "NuMIupPEnoise"){
        hist->SetLineColor(kBlue+1);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "PE Noise Up", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBnoiseAmpDown" || variation == "NuMInoiseAmpDown"){
        hist->SetLineColor(42);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Noise Amp. Down", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBdownPEnoise"  || variation == "NuMIdownPEnoise"){
        hist->SetLineColor(50);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "PE Noise Down", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBnoiseAmpUp"   || variation == "NuMInoiseAmpUp"){
        hist->SetLineColor(kOrange+10);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Noise Amp. Up", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBDTdown"       || variation == "NuMIDTdown"){
        hist->SetLineColor(kOrange+1);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "DT Down", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBDTup"         || variation == "NuMIDTup"){
        hist->SetLineColor(kMagenta-10);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "DT Up", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBDLup"         || variation == "NuMIDLup"){
        hist->SetLineColor(kMagenta);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "DL Up", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBDLdown"       || variation == "NuMIDLdown"){
        hist->SetLineColor(kTeal+6);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "DL Down", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBdataSCE"      || variation == "NuMIdataSCE"){
        hist->SetLineColor(kAzure-9);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "SCE", "l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else if  (variation == "BNBLArG4BugFix"  || variation == "NuMILArG4BugFix"){
        
        hist->SetLineWidth(2);
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());

        if (weighted_str == "BNBWeighted_Variations"){
            hist->SetLineColor(kSpring-7);
            legend->AddEntry(hist, "Weighted BNB LArG4BugFix", "l");
        } 
        else{
            hist->SetLineColor(kBlack);
            legend->AddEntry(hist, "NuMI LArG4BugFix", "l");
        }
    }
    else if  (variation == "BNBBirksRecomb"  || variation == "NuMIBirksRecomb"){
        hist->SetLineColor(kRed+1);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Birks Recomb.","l");
        hist->SetLineStyle(1);
        hist->Draw(draw_spec.c_str());
    }
    else {
        std::cout << "Error! Could not match varations:\t" << variation << std::endl;
        return;
    }


}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::GetBNBBkgWeight(double theta, double phi, double phi_wrapped, std::string bkg_class, double &weight_all, double &weight_indiv, std::string dirname ){
    
    // Grab the variation folders in the file
    std::vector<std::string> variations = {
                                        "BNBCV",
                                        "BNBwithDIC",
                                        "BNBdataSCE",
                                        "BNBdeadSaturatedChannels",
                                        "BNBLArG4BugFix",
                                        "BNBDLup",
                                        "BNBDLdown",
                                        "BNBDTup",
                                        "BNBDTdown",
                                        "BNBnoiseAmpUp",
                                        "BNBnoiseAmpDown",
                                        "BNBaltDeadChannels",
                                        "BNBstretchResp",
                                        "BNBsqueezeResp",
                                        "BNBupPEnoise",
                                        "BNBEnhancedTPCVis",
                                        "BNBBirksRecomb",
                                        "BNBdownPEnoise"};

    int p = -1; // Index of variation

    // Get the index of variation to weight
    for (unsigned int j=0; j< variations.size(); j++){
        if (dirname == variations.at(j)) p = j;
    }

    // No variation was matched so we just return
    if (p == -1) return;

    double xbin{1.0}, ybin{1.0};
    TH2D* hist;
    double weight{0.0};

    double mc_phi_wrapped = WrapPhi(phi);

    char name[500];

    // First weight by the total histograms
    std::string histname_theta_phi = "h_ThetaBkg_Phi_wrapped_ratio";
    snprintf(name, 500, "%s", histname_theta_phi.c_str() );
    hist = (TH2D*)fweight->Get(name);
    if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
    xbin = hist->GetXaxis()->FindBin(theta);
    ybin = hist->GetYaxis()->FindBin(mc_phi_wrapped);
    weight = hist->GetBinContent(xbin, ybin);

    if (weight == 0) weight_all+=1.0;
    else weight_all+=weight;


    double cosmic_sf = 181.0/153.0; // Ratio of cosmic events to bnb events selected

    // // Now try to weight split using background categories
    // if (bkg_class == "pi0_gamma") histname_theta_phi = "h_ThetaBkg_pi0_Phi_ratio";
    // if (bkg_class == "cosmic")    histname_theta_phi = "h_ThetaBkg_cosmic_Phi_ratio";
    // if (bkg_class == "other_bkg") histname_theta_phi = "h_ThetaBkg_other_Phi_ratio";

    // snprintf(name, 500, "%s", histname_theta_phi.c_str() );
    // hist = (TH2D*)fweight->Get(name);
    // if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
    // xbin = hist->GetXaxis()->FindBin(theta);
    
    // if (bkg_class != "cosmic") ybin = hist->GetYaxis()->FindBin(phi); // not cosmic
    // else ybin = hist->GetYaxis()->FindBin(phi); //cosmic
    
    // weight = hist->GetBinContent(xbin, ybin);

    // // if (bkg_class == "cosmic") weight = cosmic_sf;

    // if (weight == 0) weight=1.0;

    // // if (weight > 9) weight = 1.0;
    
    // weight_indiv+=weight;


    // Here we weight in only phi
    TH1D* hist_1D;
    std::string histname_phi;
    
    if (bkg_class == "pi0_gamma") histname_phi = "h_shower_phi_pi0_ratio";
    if (bkg_class == "cosmic")    histname_phi = "h_shower_phi_bkg_cosmic_ratio";
    if (bkg_class == "other_bkg") histname_phi = "h_shower_phi_other_ratio";

    snprintf(name, 500, "%s", histname_phi.c_str() );
    hist_1D = (TH1D*)fweight->Get(name);
    
    if (hist_1D == NULL ) std::cout << "ERROR: Can't get Histogram!: " << name << std::endl;
    
    xbin = hist_1D->GetXaxis()->FindBin(phi);
    
    weight = hist_1D->GetBinContent(xbin);

    // if (bkg_class == "cosmic") weight = cosmic_sf;

    if (weight == 0) weight=1.0;

    // if (weight > 9) weight = 1.0;
    
    weight_indiv+=weight;

    std::cout << "Weight: " << weight << "  " << weight_indiv <<"  " << mc_phi_wrapped<< "  " << theta << "  " <<bkg_class <<   std::endl;

}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::FillHistogramCuts( xsecAna::TPCObjectContainer tpc_obj, std::string histname, bool bool_sig){

    n_pfp = tpc_obj.NumPFParticles();

    const int leading_shower_index = GetLeadingShowerIndex(n_pfp, 0, tpc_obj);

    for (int j = 0; j < n_pfp ; j++){

        auto const pfp_obj = tpc_obj.GetParticle(j);

        // PFP vars                
        const int  num_pfp_hits       = pfp_obj.NumPFPHits();
        const int  num_pfp_hits_w     = pfp_obj.NumPFPHitsW(); // Collection plane hits
        
        const double mc_open_angle  = pfp_obj.mcOpenAngle();
        const double leading_dedx   = pfp_obj.PfpdEdx().at(2);//just the collection plane!

        // Background events
        if (!bool_sig) {

            if (j == leading_shower_index){

                if (histname == "kpre_hit_threshold_cut")             TH1D_hist.at(kpre_hit_threshold_cut)            ->Fill(num_pfp_hits);
                if (histname == "kpost_hit_threshold_cut")            TH1D_hist.at(kpost_hit_threshold_cut)           ->Fill(num_pfp_hits);
                if (histname == "kpre_hit_threshold_collection_cut")  TH1D_hist.at(kpre_hit_threshold_collection_cut) ->Fill(num_pfp_hits_w);
                if (histname == "kpost_hit_threshold_collection_cut") TH1D_hist.at(kpost_hit_threshold_collection_cut)->Fill(num_pfp_hits_w);
                if (histname == "kpre_open_angle_cut")                TH1D_hist.at(kpre_open_angle_cut)               ->Fill(mc_open_angle * 180/3.14159);
                if (histname == "kpost_open_angle_cut")               TH1D_hist.at(kpost_open_angle_cut)              ->Fill(mc_open_angle * 180/3.14159);
                if (histname == "kpre_dedx_cut")                      TH1D_hist.at(kpre_dedx_cut)                     ->Fill(leading_dedx);
                if (histname == "kpost_dedx_cut")                     TH1D_hist.at(kpost_dedx_cut)                    ->Fill(leading_dedx);
            }

        }

    }

}
//***************************************************************************
//***************************************************************************