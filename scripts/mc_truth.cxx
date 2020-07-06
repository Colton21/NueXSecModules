#include "mc_truth.h"

//***************************************************************************
//***************************************************************************

// ----------------------
//         Main
// ----------------------
void mc_truth::run_var(const char * _file1, const std::vector<double> _config) {
    
    // create plots folder if it does not exist
    system("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi");

    // Set this bool to call function to merge weights into main tree
    if (false){
        merge_weights_into_tree();
        return;
    }

    // Set this to true to make the histograms 
    if (true){
        runweightmode();
        return;
    }

    gStyle->SetOptStat(0); // say no to stats box

    TString mode = "numi";

    std::string infile = std::string(_file1);

    if (infile == "files/filter_ext_v7_full.root") type = "ext";
    else type = "mc";

    if (infile == "files/mc_truth_tree_out.root") write_tree = false;

    std::cout << "Running in mode: "<< type << std::endl;

    //*************************** Configure the cut parameters *************************************
    std::cout << "\n --- Configuring Parameters --- \n" << std::endl;
    _x1 = _config[0];
    _x2 = _config[1];
    _y1 = _config[2];
    _y2 = _config[3];
    _z1 = _config[4];
    _z2 = _config[5];
    flash_pe_threshold = _config[6];
    flash_time_start = _config[7]; // Use default numi config
    flash_time_end   = _config[8]; // Use default numi config
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

    // Create output file
    TFile *outFile_root;
    TTree *mc_tree_out;
    if (write_tree){
        outFile_root = TFile::Open("files/mc_truth_tree_out.root", "UPDATE");
        mc_tree_out = new TTree("mc_tree_out", "mc_tree_out");
    } 
    else {
        outFile_root = TFile::Open("files/mc_truth_tree_out_dummy.root", "UPDATE");
        mc_tree_out = new TTree("mc_tree_out_dummy", "mc_tree_out_dummy");
    }
    

    // Create output file
    TFile *outFile = TFile::Open("files/mc_truth_out.root", "UPDATE");
    
    // Get input file
    TFile* inFile = TFile::Open(_file1, "READ");
    inFile->cd();

    // Force the tree to be associated with the correct output file
    mc_tree_out->SetDirectory(outFile_root);

    
    // Get TPC Obj Information from file
    TTree* TPCObjTree;
    
    if (infile == "files/filter_numi_cosmic_v4.root") TPCObjTree = (TTree*) inFile->Get("AnalyzeTPCO/tree");    
    if (infile == "files/mc_truth_tree_out.root") TPCObjTree = (TTree*) inFile->Get("mc_tree_out");
    TPCObjTree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

    // Num events in the tree
    const int tree_total_entries = TPCObjTree->GetEntries();
    std::cout << "Total Events: " << tree_total_entries << std::endl;

    // Get the largest Flash Vector of Vector    
    std::vector<std::vector<double>> largest_flash_v_v;
    
    if (write_tree) largest_flash_v_v = GetLargestFlashVector(inFile, flash_time_start, flash_time_end);

    //************* Get list of flashes that pass flash in time and flash pe cut ************************
    std::vector<bool> flash_cuts_pass_vec;
    if (write_tree){
        
        if (type == "mc") FlashinTime_FlashPE(inFile, flash_time_start, flash_time_end, flash_cuts_pass_vec, mode );
        else  FlashinTime_FlashPE(inFile, flash_time_start+0.343, flash_time_end+0.343, flash_cuts_pass_vec, mode );
    }
    //***************************************************************************************************

    // MCTruth Counting Tree
    TTree * mctruth_counter_tree;
    if (type == "mc"){
        
        if (infile == "files/filter_numi_cosmic_v4.root") mctruth_counter_tree = (TTree*)inFile->Get("AnalyzeTPCO/mctruth_counter_tree"); 
        if (infile == "files/mc_truth_tree_out.root") mctruth_counter_tree = (TTree*) inFile->Get("mc_tree_out");
        
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
        if (infile == "files/mc_truth_tree_out.root") {
            mctruth_counter_tree->SetBranchAddress("ccqe_weights", &ccqe_weights);
            mctruth_counter_tree->SetBranchAddress("ccmec_weights", &ccmec_weights);
        }

        // Output tree stuff
        mc_tree_out->Branch("run",    &run);
        mc_tree_out->Branch("subrun", &subrun);
        mc_tree_out->Branch("evt",    &evt);
        mc_tree_out->Branch("fMCNuEnergy", &mc_nu_energy);
        mc_tree_out->Branch("fMCNuID",     &mc_nu_id);
        mc_tree_out->Branch("fMCNuVtxX",   &mc_nu_vtx_x);
        mc_tree_out->Branch("fMCNuVtxY",   &mc_nu_vtx_y);
        mc_tree_out->Branch("fMCNuVtxZ",   &mc_nu_vtx_z);
        mc_tree_out->Branch("fMCNuDirX",   &mc_nu_dir_x);
        mc_tree_out->Branch("fMCNuDirY",   &mc_nu_dir_y);
        mc_tree_out->Branch("fMCNuDirZ",   &mc_nu_dir_z);
        mc_tree_out->Branch("fMCNumParticles", &mc_nu_num_particles);
        mc_tree_out->Branch("fMCNumChargedParticles", &mc_nu_num_charged_particles);
        mc_tree_out->Branch("fMCEleDirX", &mc_ele_dir_x);
        mc_tree_out->Branch("fMCEleDirY", &mc_ele_dir_y);
        mc_tree_out->Branch("fMCEleDirZ", &mc_ele_dir_z);
        mc_tree_out->Branch("fMCEleEnergy",   &mc_ele_energy);
        mc_tree_out->Branch("fMCEleMomentum", &mc_ele_momentum);
        mc_tree_out->Branch("has_pi0",   &has_pi0);
        mc_tree_out->Branch("fMCNuTime", &mc_nu_time);
        mc_tree_out->Branch("TpcObjectContainerV", &tpc_object_container_v_copy);
        mc_tree_out->Branch("passed_selection", &passed_selection);
        mc_tree_out->Branch("tpc_object_pass", &tpc_object_pass);
    }

    // Define the FV
    std::vector<double> fv_boundary_v = {_x1, _x2, _y1, _y2, _z1, _z2};

    if (type == "mc") std::cout << "Total Events (MC): " << mctruth_counter_tree->GetEntries() << std::endl;

    // ----------------------
    //      Event loop
    // ----------------------
    std::cout << "Starting Eventloop..." << std::endl;
    for(int event = 0; event < tree_total_entries; event++){

        if (gen_event || passed_selection){
            mc_tree_out->Fill();
        }

        passed_selection = false;

        if (event % 100000 == 0) std::cout << "On entry " << event/100000.0 <<"00k" << std::endl;
        
        TPCObjTree->GetEntry(event);
        if (type == "mc") mctruth_counter_tree->GetEntry(event);

        tpc_object_container_v_copy = *tpc_object_container_v;

        int n_tpc_obj = tpc_object_container_v->size();

        // Skip numu events
        // if (mc_nu_id == 2 || mc_nu_id == 4 || mc_nu_id == 6 || mc_nu_id == 8) continue;

        // Check to see in truth the vertex was inside the tpc FV
        bool true_in_tpc;
        if (type == "mc") true_in_tpc = in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v);

        if(true_in_tpc == true && (mc_nu_id == 1 || mc_nu_id == 5)) {
            if (tpc_object_container_v->size() != 0){
                
                auto const tpc_object_container = tpc_object_container_v->at(0);
                evt  = tpc_object_container.EventNumber();
                run    = tpc_object_container.RunNumber();
                subrun = tpc_object_container.SubRunNumber();
                gen_event = true;

            }
        }
        else {
            gen_event = false;
        }

        // --------------- MC Counters ---------------
        if (type == "mc"){
            if (mc_nu_id == 1 && true_in_tpc) mc_nue_cc_counter++;
            if (mc_nu_id == 3 && true_in_tpc) mc_nue_nc_counter++;
            if (mc_nu_id == 5 && true_in_tpc) mc_nue_cc_counter_bar++;
            if (mc_nu_id == 7 && true_in_tpc) mc_nue_nc_counter_bar++;
        
        
            // Filter for truth Nues
            // if (mc_nu_id == 2 || mc_nu_id == 4 || mc_nu_id == 6 || mc_nu_id == 8) continue;
        }

        // --------------- Flash Information ---------------
        std::vector<double> largest_flash_v;
        
        if (write_tree){
            largest_flash_v = largest_flash_v_v.at(event); // Vec with the largest flash
            largest_flash_y     = largest_flash_v.at(0);
            largest_flash_z     = largest_flash_v.at(1);
            largest_flash_time     = largest_flash_v.at(2);
            largest_flash_pe     = largest_flash_v.at(3);
        }
        
        // -------------------------------------------------

        tpc_object_pass.clear();
        tpc_object_pass.resize(n_tpc_obj, 0);

        // Loop over TPCObj
        for (int i = 0; i < n_tpc_obj; i++){
            auto const tpc_obj = tpc_object_container_v->at(i);

            // Classify the event  
            std::pair<std::string, int> tpc_classification;
            if (type == "mc") tpc_classification = TPCO_Classifier(tpc_obj, true_in_tpc, has_pi0);
            else tpc_classification = std::make_pair("ext", 0);

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

            double Flash_TPCObj_Dist{0.0};
            if (write_tree) Flash_TPCObj_Dist = Flash_TPCObj_vtx_Dist(tpc_obj_vtx_y, tpc_obj_vtx_z, largest_flash_y, largest_flash_z);

            if (bool_sig) nue_cc_counter++;

            // Loop over the Par Objects
            for (int j = 0; j < n_pfp ; j++){

                auto const pfp_obj = tpc_obj.GetParticle(j);

                // PFP vars                
                const std::string mc_origin = pfp_obj.Origin();
                const int  pfp_pdg          = pfp_obj.PFParticlePdgCode();
                const int  num_pfp_hits     = pfp_obj.NumPFPHits();
                const int  n_pfp_hits_w     = pfp_obj.NumPFPHitsW(); // Collection plane hits
                
                const double pfp_length     = pfp_obj.pfpLength();
                const double pfp_open_angle = pfp_obj.pfpOpenAngle();
                double leading_dedx   = pfp_obj.PfpdEdx().at(2);//just the collection plane!
                if (type == "mc") leading_dedx = leading_dedx * (196.979 /242.72);

                const double pfp_vtx_x =  pfp_obj.pfpVtxX();
                const double pfp_vtx_y =  pfp_obj.pfpVtxY();
                const double pfp_vtx_z =  pfp_obj.pfpVtxZ();

                const double pfp_dir_x = pfp_obj.pfpDirX();
                const double pfp_dir_y = pfp_obj.pfpDirY();
                const double pfp_dir_z = pfp_obj.pfpDirZ();

                double mc_Theta{-1.0}, mc_Phi{-1.0}, mc_Energy{-1.0}, mc_pdg{-1.0}, mc_pdg_parent{-1.0};
                
                if (type == "mc"){
                    mc_Theta  = pfp_obj.mcTheta();
                    mc_Phi    = pfp_obj.mcPhi();
                    mc_Energy = pfp_obj.mcEnergy();
                    mc_pdg    = pfp_obj.MCPdgCode();
                    mc_pdg_parent = pfp_obj.MCParentPdg();
                }

                const double pfp_end_x = (pfp_obj.pfpVtxX() + (pfp_length * pfp_dir_x));
                const double pfp_end_y = (pfp_obj.pfpVtxY() + (pfp_length * pfp_dir_y));
                const double pfp_end_z = (pfp_obj.pfpVtxZ() + (pfp_length * pfp_dir_z));

                std::vector<double> pfp_start_vtx {pfp_vtx_x, pfp_vtx_y, pfp_vtx_z};
                std::vector<double> pfp_end_vtx {pfp_end_x, pfp_end_y, pfp_end_z};

                // // Take a look at protons
                // if (mc_pdg == 2212) {
                //     // std::cout << "proton: " << mc_Energy << std::endl;
                // }

                // // Take a look at electrons
                // if (mc_pdg == 11){
                //     // std::cout << "electron: " << mc_Energy << std::endl;
                // }

                const double pfp_Nu_vtx_Dist =  pfp_vtx_distance(tpc_obj_vtx_x, tpc_obj_vtx_y, tpc_obj_vtx_z, pfp_vtx_x, pfp_vtx_y, pfp_vtx_z);
                

            } // END LOOP PAR OBJ



            //************************ Apply Flash in time and Flash PE cut *************************************
            if (write_tree){
                if (flash_cuts_pass_vec.at(event) == false) continue;
                if (bool_sig) counter_FlashinTime_FlashPE++;
            }
            
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
            if (write_tree){
                bool bool_flashvtx = flashRecoVtxDist(largest_flash_v, tolerance, tpc_obj_vtx_x, tpc_obj_vtx_y,  tpc_obj_vtx_z);
                if ( bool_flashvtx == false ) continue;
                if (bool_sig) counter_FlashRecoVtxDist++;
            }
            
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
            // std::string hist_fill_name = "kpre_hit_threshold_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            //****************************** Apply Hit threshold cut**** cut ************************************
            bool bool_hitTh = HitThreshold(tpc_obj, shwr_hit_threshold, false);
            if ( bool_hitTh == false ) continue;
            if (bool_sig) counter_HitThresh++;
            //***************************************************************************************************
            // hist_fill_name = "kpost_hit_threshold_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            // hist_fill_name = "kpre_hit_threshold_collection_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            //****************************** Apply Hit threshold collection cut *********************************
            bool bool_hitThW = HitThreshold(tpc_obj, shwr_hit_threshold_collection, true);
            if ( bool_hitThW == false ) continue;
            if (bool_sig) counter_HitThreshW++;
            //***************************************************************************************************
            // hist_fill_name = "kpost_hit_threshold_collection_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            // hist_fill_name = "kpre_open_angle_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            //****************************** Apply Open Angle cut ***********************************************
            bool bool_OpenAngle = OpenAngleCut(tpc_obj, tolerance_open_angle);
            if ( bool_OpenAngle == false ) continue;
            if (bool_sig) counter_OpenAngle++;
            //***************************************************************************************************
            // hist_fill_name = "kpost_open_angle_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            // hist_fill_name = "kpre_dedx_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
            //****************************** Apply dEdx cut *****************************************************
            bool bool_dEdx = dEdxCut(tpc_obj, tolerance_dedx_min, tolerance_dedx_max);
            if ( bool_dEdx == false ) continue;
            if (bool_sig) counter_dEdx++;
            //***************************************************************************************************
            // hist_fill_name = "kpost_dedx_cut"; FillHistogramCuts( tpc_obj, hist_fill_name, bool_sig);
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

            passed_selection = true;
            tpc_object_pass.at(i) = 1;

            if (type != "mc") counter_ext++;

            // Loop over the Par Objects
            for (int j = 0; j < n_pfp ; j++){

                auto const pfp_obj = tpc_obj.GetParticle(j);

                // PFP vars                
                const std::string mc_origin = pfp_obj.Origin();
                const int  pfp_pdg          = pfp_obj.PFParticlePdgCode();
                const int  num_pfp_hits     = pfp_obj.NumPFPHits();
                const int  n_pfp_hits_w     = pfp_obj.NumPFPHitsW(); // Collection plane hits
                
                const double pfp_length     = pfp_obj.pfpLength();
                const double pfp_open_angle = pfp_obj.pfpOpenAngle();
                double leading_dedx   = pfp_obj.PfpdEdx().at(2);//just the collection plane!
                if (type == "mc") leading_dedx = leading_dedx * (196.979 /242.72);

                const double pfp_vtx_x =  pfp_obj.pfpVtxX();
                const double pfp_vtx_y =  pfp_obj.pfpVtxY();
                const double pfp_vtx_z =  pfp_obj.pfpVtxZ();

                const double pfp_dir_x = pfp_obj.pfpDirX();
                const double pfp_dir_y = pfp_obj.pfpDirY();
                const double pfp_dir_z = pfp_obj.pfpDirZ();

                

                double mc_Theta{-1.0}, mc_Phi{-1.0}, mc_Energy{-1.0}, mc_pdg{-1.0}, mc_pdg_parent{-1.0};
                
                if (type == "mc"){
                    mc_Theta  = pfp_obj.mcTheta();
                    mc_Phi    = pfp_obj.mcPhi();
                    mc_Energy = pfp_obj.mcEnergy();
                    mc_pdg    = pfp_obj.MCPdgCode();
                    mc_pdg_parent = pfp_obj.MCParentPdg();
                }

                const double pfp_end_x = (pfp_obj.pfpVtxX() + (pfp_length * pfp_dir_x));
                const double pfp_end_y = (pfp_obj.pfpVtxY() + (pfp_length * pfp_dir_y));
                const double pfp_end_z = (pfp_obj.pfpVtxZ() + (pfp_length * pfp_dir_z));

                std::vector<double> pfp_start_vtx {pfp_vtx_x, pfp_vtx_y, pfp_vtx_z};
                std::vector<double> pfp_end_vtx {pfp_end_x, pfp_end_y, pfp_end_z};

                const double shr_vtx_dist =  pfp_vtx_distance(tpc_obj_vtx_x, tpc_obj_vtx_y, tpc_obj_vtx_z, pfp_vtx_x, pfp_vtx_y, pfp_vtx_z);

                // Leading shower candidate
                if (j == leading_shower_index){


                    TVector3 shower_vector(pfp_dir_x, pfp_dir_y, pfp_dir_z);
                    TVector3 numi_vector;
                    numi_vector.SetMagThetaPhi(1, 0, 0);
                    const double leading_shower_theta = acos(shower_vector.Dot(numi_vector) / (shower_vector.Mag() * numi_vector.Mag())) * (180/3.1415);

                    if (type == "mc"){
                        if (bool_sig) h_shr_dEdx_shr_dist_sig->Fill(leading_dedx, shr_vtx_dist);
                        if (!bool_sig) h_shr_dEdx_shr_dist_bkg->Fill(leading_dedx, shr_vtx_dist);

                        if (leading_shower_theta>=0 && leading_shower_theta <=60){
                            if (bool_sig) h_shr_dEdx_shr_dist_sig_good->Fill(leading_dedx, shr_vtx_dist);
                            if (!bool_sig) h_shr_dEdx_shr_dist_bkg_good->Fill(leading_dedx, shr_vtx_dist);
                        }
                        
                    }
                    else {
                        h_shr_dEdx_shr_dist_ext->Fill(leading_dedx, shr_vtx_dist);
                        
                        if (leading_shower_theta>=0 && leading_shower_theta <= 60){
                            h_shr_dEdx_shr_dist_ext_good->Fill(leading_dedx, shr_vtx_dist);
                        }
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
    std::cout << "EXT                          --- " << counter_ext << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "--------------- Reco COUNTERS ---------------------" << std::endl;
    std::cout << "True Nue CC in FV --- " << nue_cc_counter << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;


    // Write histograms to file
    
    outFile->cd();

    if (type == "mc"){
        h_shr_dEdx_shr_dist_sig->SetOption("box");
        h_shr_dEdx_shr_dist_sig->SetFillColor(30);
        h_shr_dEdx_shr_dist_sig->Write("",TObject::kOverwrite);

        h_shr_dEdx_shr_dist_bkg->SetOption("box");
        h_shr_dEdx_shr_dist_bkg->SetFillColor(46);
        h_shr_dEdx_shr_dist_bkg->Write("",TObject::kOverwrite);

        h_shr_dEdx_shr_dist_sig_good->SetOption("box");
        h_shr_dEdx_shr_dist_sig_good->SetFillColor(30);
        h_shr_dEdx_shr_dist_sig_good->Write("",TObject::kOverwrite);

        h_shr_dEdx_shr_dist_bkg_good->SetOption("box");
        h_shr_dEdx_shr_dist_bkg_good->SetFillColor(46);
        h_shr_dEdx_shr_dist_bkg_good->Write("",TObject::kOverwrite);
    }
    else {
        h_shr_dEdx_shr_dist_ext->SetOption("box");
        h_shr_dEdx_shr_dist_ext->SetFillColor(46);
        h_shr_dEdx_shr_dist_ext->Write("",TObject::kOverwrite);

        h_shr_dEdx_shr_dist_ext_good->SetOption("box");
        h_shr_dEdx_shr_dist_ext_good->SetFillColor(46);
        h_shr_dEdx_shr_dist_ext_good->Write("",TObject::kOverwrite);

    }

    // Here we create another root file where we can save a slimmed version of the selection for re-weighting
    if (write_tree){
        outFile_root->cd();
        mc_tree_out->Write("",TObject::kOverwrite);
    }
    
} // END MAIN
//***************************************************************************
//***************************************************************************
int mc_truth::GetLeadingShowerIndex(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj){
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
double mc_truth::GetLongestTrackLength(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj){
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
void mc_truth::GetNumber_Track_Shower(const int n_pfp, int n_tpc_obj,
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
double mc_truth::pfp_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z) {
    const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );

    return distance;
}
//***************************************************************************
//************************** Flash Functions ********************************
std::vector<std::vector<double>> mc_truth::GetLargestFlashVector(TFile* f, double flash_start_time, double flash_end_time ){

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
double mc_truth::Flash_TPCObj_vtx_Dist(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z) {
    const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
    return distance;
}
//***************************************************************************
//***************************************************************************
std::pair<std::string, int> mc_truth::TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool true_in_tpc, bool has_pi0) {
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
        if(part_nue_cc  > 0 || part_nue_bar_cc > 0)  { return std::make_pair("nue_cc_qe",  leading_index); }
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
bool mc_truth::in_fv(double x, double y, double z, std::vector<double> fv_boundary_v) {
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
bool mc_truth::flash_in_time(double flash_time, double flash_start, double flash_end) {
    if(flash_time >= flash_start && flash_time <= flash_end) return true; // Pass in time
    else return false; // Fail in time
}
//***************************************************************************
//***************************************************************************
// Flash PE
bool mc_truth::flash_pe(int flash_pe, int flash_pe_threshold) {
    if (flash_pe >= flash_pe_threshold) return true; // Pass PE Thresh
    else return false; // Fail PE Thresh
}
//***************************************************************************
//***************************************************************************
// Get the vector of flashes with true/false on whether it passed the selection or not
void mc_truth::FlashinTime_FlashPE(TFile* f, double flash_start_time, double flash_end_time, std::vector<bool> &flash_cuts_pass_vec, TString mode ){
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
bool mc_truth::HasNue(xsecAna::TPCObjectContainer tpc_obj, const int n_pfp ) {

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
bool mc_truth::opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance) {
    const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
    
    if(distance <= tolerance) return true;
    return false;
}
bool mc_truth::flashRecoVtxDist(std::vector< double > largest_flash_v, double tolerance, const double tpc_vtx_x, const double tpc_vtx_y, const double tpc_vtx_z) {
    
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
bool mc_truth::VtxNuDistance(xsecAna::TPCObjectContainer tpc_obj,int pfp_pdg_type , double tolerance){
    
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
bool mc_truth::HitThreshold(xsecAna::TPCObjectContainer tpc_obj, double threshold, bool useCollection){

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
bool mc_truth::OpenAngleCut(xsecAna::TPCObjectContainer tpc_obj, const std::vector<double> tolerance_open_angle){

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
bool mc_truth::dEdxCut( xsecAna::TPCObjectContainer tpc_obj, const double tolerance_dedx_min, const double tolerance_dedx_max){

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
    if (type == "mc") leading_dedx = leading_dedx * (196.979 /242.72);

    if (leading_dedx <= tolerance_dedx_max && leading_dedx >= tolerance_dedx_min) return true;
     
    return false;

}
//***************************************************************************
//***************************************************************************
// Secondar Showers Cut
bool mc_truth::SecondaryShowersDistCut(xsecAna::TPCObjectContainer tpc_obj, const double dist_tolerance){

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
bool mc_truth::HitLengthRatioCut(const double pfp_hits_length_tolerance, xsecAna::TPCObjectContainer tpc_obj){
    
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
bool mc_truth::LongestTrackLeadingShowerCut(const double ratio_tolerance, xsecAna::TPCObjectContainer tpc_obj){

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
bool mc_truth::IsContained(std::vector<double> track_start, std::vector<double> track_end, std::vector<double> fv_boundary_v) {
    
    if(in_fv(track_start.at(0), track_start.at(1), track_start.at(2), fv_boundary_v) == true
       && in_fv(track_end.at(0), track_end.at(1), track_end.at(2), fv_boundary_v) == true) {
        return true;
    }
    else 
        return false;
}
bool mc_truth::ContainedTracksCut(std::vector<double> fv_boundary_v, xsecAna::TPCObjectContainer tpc_obj){

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
void mc_truth::merge_weights_into_tree(){

    // Get the file with the event weights in
    TFile* weight_file = TFile::Open("files/NuMIEventWeight_CCQE_CCMEC.root");
    TTree* weight_tree = (TTree*) weight_file->Get("eventweightreader/tree");    
    std::vector<double> * ccqe_weights = nullptr;
    std::vector<double> * ccmec_weights = nullptr;
    int run_w, subrun_w, event_w;
    weight_tree->SetBranchAddress("run", &run_w);
    weight_tree->SetBranchAddress("subrun", &subrun_w);
    weight_tree->SetBranchAddress("event", &event_w);
    weight_tree->SetBranchAddress("ccqe", &ccqe_weights);
    weight_tree->SetBranchAddress("ccmec", &ccmec_weights);
    const int weighttree_total_entries = weight_tree->GetEntries();
    std::cout << "Total Weighted Events: " << weighttree_total_entries << std::endl;
    weight_tree->SetDirectory(weight_file);


    // Get the main file
    TFile *outFile_root = TFile::Open("files/mc_truth_tree_out.root", "UPDATE");
    TTree *mctruth_counter_tree = (TTree*) outFile_root->Get("mc_tree_out");
    mctruth_counter_tree->SetDirectory(outFile_root);

    mctruth_counter_tree->SetBranchAddress("run", &run);
    mctruth_counter_tree->SetBranchAddress("subrun", &subrun);
    mctruth_counter_tree->SetBranchAddress("evt", &evt);

    TBranch *b_ccqe  = mctruth_counter_tree->Branch("ccqe_weights", &ccqe_weights);
    TBranch *b_ccmec = mctruth_counter_tree->Branch("ccmec_weights",&ccmec_weights);

    const int tree_total_entries = mctruth_counter_tree->GetEntries();

    for (int event = 0; event < tree_total_entries; event++){

        if (event % 1000 == 0) std::cout << "On entry " << event/1000.0 <<"k" << std::endl;

        mctruth_counter_tree->GetEntry(event);

        // Loop over the weight tree and get the weight vectors for this specific run subrun event
        for(int evt_w = 0; evt_w < weighttree_total_entries; evt_w++){

            weight_tree->GetEntry(evt_w);

            if (event_w == evt && subrun_w == subrun && run_w == run){
                // std::cout << "Found match of subruns" << std::endl;
                b_ccqe->Fill();
                b_ccmec->Fill();
            }
        }

    }
    
    outFile_root->cd();
    mctruth_counter_tree->Write("",TObject::kOverwrite);

}
//***************************************************************************
//***************************************************************************
void mc_truth::runweightmode(){

    gStyle->SetOptStat(0);

    // Get the main file
    TFile *outFile_root = TFile::Open("files/mc_truth_tree_out.root", "READ");
    TTree *mctruth_counter_tree = (TTree*) outFile_root->Get("mc_tree_out");
    mctruth_counter_tree->SetDirectory(outFile_root);

    mctruth_counter_tree->SetBranchAddress("run", &run);
    mctruth_counter_tree->SetBranchAddress("subrun", &subrun);
    mctruth_counter_tree->SetBranchAddress("evt", &evt);
    mctruth_counter_tree->SetBranchAddress("ccqe_weights", &ccqe_weights);
    mctruth_counter_tree->SetBranchAddress("ccmec_weights", &ccmec_weights);
    mctruth_counter_tree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);
    mctruth_counter_tree->SetBranchAddress("passed_selection", &passed_selection);
    mctruth_counter_tree->SetBranchAddress("fMCNuID", &mc_nu_id);
    mctruth_counter_tree->SetBranchAddress("fMCNuEnergy", &mc_nu_energy);
    mctruth_counter_tree->SetBranchAddress("fMCEleEnergy",   &mc_ele_energy);
    mctruth_counter_tree->SetBranchAddress("fMCNumChargedParticles", &mc_nu_num_charged_particles);
    mctruth_counter_tree->SetBranchAddress("fMCNuDirX",   &mc_nu_dir_x);
    mctruth_counter_tree->SetBranchAddress("fMCNuDirY",   &mc_nu_dir_y);
    mctruth_counter_tree->SetBranchAddress("fMCNuDirZ",   &mc_nu_dir_z);
    mctruth_counter_tree->SetBranchAddress("fMCEleDirX", &mc_ele_dir_x);
    mctruth_counter_tree->SetBranchAddress("fMCEleDirY", &mc_ele_dir_y);
    mctruth_counter_tree->SetBranchAddress("fMCEleDirZ", &mc_ele_dir_z);
    mctruth_counter_tree->SetBranchAddress("tpc_object_pass", &tpc_object_pass_temp);

    const int tree_total_entries = mctruth_counter_tree->GetEntries();

    // Lets define our histograms first index of vector is reweightor type
    std::vector<std::vector<TH1D*>> h_proton_num; // proton efficiency numerator
    std::vector<std::vector<TH1D*>> h_proton_den; // proton efficiency denomiator

    std::vector<std::vector<TH1D*>> h_electron_num; // electron efficiency numerator
    std::vector<std::vector<TH1D*>> h_electron_den; // electron efficiency denomiator

    TH1D* h_proton_num_cv = new TH1D("h_proton_num_CV", "",10, 0, 1.2 );
    TH1D* h_proton_den_cv = new TH1D("h_proton_den_CV", "",10, 0, 1.2 );

    TH1D* h_electron_num_cv = new TH1D("h_electron_num_CV", "",10, 0, 2 );
    TH1D* h_electron_den_cv = new TH1D("h_electron_den_CV", "",10, 0, 2 );

    TH1D* h_E_transfer_CV_den = new TH1D("h_E_transfer_CV_den", ";Energy Transfer [GeV]; Entries",10, 0, 2 );
    TH1D* h_P_transfer_CV_den = new TH1D("h_P_transfer_CV_den", ";Momentum Transfer [GeV/c]; Entries",10, 0, 2 );
    TH1D* h_NumCharged_CV_den = new TH1D("h_NumCharged_CV_den", ";Number of Charged Particles; Entries",10, 0, 10 );

    TH1D* h_E_transfer_CV_num = new TH1D("h_E_transfer_CV_num", ";Energy Transfer [GeV]; Entries",10, 0, 2 );
    TH1D* h_P_transfer_CV_num = new TH1D("h_P_transfer_CV_num", ";Momentum Transfer [GeV/c]; Entries",10, 0, 2 );
    TH1D* h_NumCharged_CV_num = new TH1D("h_NumCharged_CV_num", ";Number of Charged Particles; Entries",10, 0, 10 );

    // Energy/P Transfer and Num charged particles
    std::vector<std::vector<TH1D*>> h_E_transfer_den_v;
    std::vector<std::vector<TH1D*>> h_P_transfer_den_v;
    std::vector<std::vector<TH1D*>> h_NumCharged_den_v;

    h_proton_num.resize(2);
    h_proton_den.resize(2);
    h_electron_num.resize(2);
    h_electron_den.resize(2);
    h_E_transfer_den_v.resize(2);
    h_P_transfer_den_v.resize(2);
    h_NumCharged_den_v.resize(2);

    // Resize to the number of universes
    for (unsigned int k = 0; k < h_proton_num.size(); k++){
        h_proton_num.at(k).resize(1000);
        h_proton_den.at(k).resize(1000);
        h_electron_num.at(k).resize(1000);
        h_electron_den.at(k).resize(1000);
        h_E_transfer_den_v.at(k).resize(1000);
        h_P_transfer_den_v.at(k).resize(1000);
        h_NumCharged_den_v.at(k).resize(1000);
    }

    // Loop over the histograms and create the histograms
    for (unsigned int l = 0; l < h_proton_num.size(); l++){
        for (unsigned int p = 0; p< h_proton_num.at(l).size(); p++){
            h_proton_num.at(l).at(p) = new TH1D(Form("h_proton_num_%i_%i",l, p), ";Proton Momentum [GeV/c];Entries", 10, 0, 1.2);
            h_proton_den.at(l).at(p) = new TH1D(Form("h_proton_den_%i_%i",l, p), ";Proton Momentum [GeV/c];Entries", 10, 0, 1.2);
            h_electron_num.at(l).at(p) = new TH1D(Form("h_electron_num_%i_%i",l, p), ";Electron Energy [GeV];Entries", 10, 0, 2);
            h_electron_den.at(l).at(p) = new TH1D(Form("h_electron_den_%i_%i",l, p), ";Electron Energy [GeV];Entries", 10, 0, 2);
            h_E_transfer_den_v.at(l).at(p) = new TH1D(Form("h_E_transfer_den_v_%i_%i",l, p), ";Energy Transfer [GeV];Entries", 10, 0, 2);
            h_P_transfer_den_v.at(l).at(p) = new TH1D(Form("h_P_transfer_den_v_%i_%i",l, p), ";Momentum Transfer [GeV/c];Entries", 10, 0, 2);
            h_NumCharged_den_v.at(l).at(p) = new TH1D(Form("h_NumCharged_den_v_%i_%i",l, p), ";Number of CHarged Particles;Entries", 10, 0, 10);
        }
    }


    for (int event = 0; event < tree_total_entries; event++){

        if (event % 1000 == 0) std::cout << "On entry " << event/1000.0 <<"k" << std::endl;

        mctruth_counter_tree->GetEntry(event);

        // Only looking at signal events for now, so skip the numus
        if (mc_nu_id == 2 || mc_nu_id == 4 || mc_nu_id == 6 || mc_nu_id == 8) continue;


        int n_tpc_obj = tpc_object_container_v->size();


        // Fill the MC Truth histograms
        if (mc_nu_id == 1 || mc_nu_id == 5){

            // Calculate the momentum transfer
            TVector3 v_nu(mc_nu_dir_x, mc_nu_dir_y, mc_nu_dir_z);
            TVector3 v_ele(mc_ele_dir_x, mc_ele_dir_y, mc_ele_dir_z);

            double angle = v_nu.Dot(v_ele);

            double p_ele = std::sqrt(mc_ele_energy*mc_ele_energy - 0.511e-3*0.511e-3 );

            double p_transfer = std::sqrt( mc_nu_energy*mc_nu_energy + p_ele*p_ele - 2 * mc_nu_energy * p_ele * std::cos(angle) );

            h_E_transfer_CV_den->Fill(mc_nu_energy - mc_ele_energy);
            h_P_transfer_CV_den->Fill(p_transfer);
            h_NumCharged_CV_den->Fill(mc_nu_num_charged_particles);

            if (passed_selection){
                h_E_transfer_CV_num->Fill(mc_nu_energy - mc_ele_energy);
                h_P_transfer_CV_num->Fill(p_transfer);
                h_NumCharged_CV_num->Fill(mc_nu_num_charged_particles);
            }

            for (unsigned int t =0 ; t < ccqe_weights->size(); t++){
                h_E_transfer_den_v.at(0).at(t)->Fill(mc_nu_energy - mc_ele_energy, ccqe_weights->at(t));
                h_E_transfer_den_v.at(1).at(t)->Fill(mc_nu_energy - mc_ele_energy, ccmec_weights->at(t));

                h_P_transfer_den_v.at(0).at(t)->Fill(p_transfer, ccqe_weights->at(t));
                h_P_transfer_den_v.at(1).at(t)->Fill(p_transfer, ccmec_weights->at(t));

                h_NumCharged_den_v.at(0).at(t)->Fill(mc_nu_num_charged_particles, ccqe_weights->at(t));
                h_NumCharged_den_v.at(1).at(t)->Fill(mc_nu_num_charged_particles, ccmec_weights->at(t));
            }


        }

        // Loop over TPCObj
        for (int i = 0; i < n_tpc_obj; i++){
            // if (tpc_object_pass_temp->at(i) == 0) continue;

            auto const tpc_obj = tpc_object_container_v->at(i);

            n_pfp = tpc_obj.NumPFParticles();

            int leading_proton_index{0};
            int leading_electron_index{0};

            double proton_mom{0.0};
            double electron_energy{0.0};

            // Loop over the Par Objects -- get the leading proton and leading electron
            // Leading defined as the the pfp with the most energy
            for (int j = 0; j < n_pfp ; j++){

                auto const pfp_obj = tpc_obj.GetParticle(j);

                int mc_pdg       = pfp_obj.MCPdgCode();
                double mc_Energy = pfp_obj.mcEnergy();
                double mc_Mom    = pfp_obj.mcMomentum();

                // Take a look at protons
                if (mc_pdg == 2212) {
                    if (mc_Mom > proton_mom){
                        proton_mom = mc_Mom;
                        leading_proton_index = j;
                        

                    } 
                }

                // Take a look at electrons
                if (mc_pdg == 11){
                    if (mc_Energy > electron_energy){
                        electron_energy = mc_Energy;
                        leading_electron_index =j;
                    } 
                }

            } // END loop over pfp


            // std::cout << proton_energy<< std::endl;

            for (int j = 0; j < n_pfp ; j++){
                auto const pfp_obj = tpc_obj.GetParticle(j);


                // Leading proton
                if (j == leading_proton_index){

                    if (proton_mom == 0.0) continue;
                    
                    // Passed the selection so numerator
                    if (tpc_object_pass_temp->at(i) == 1){
                        for (unsigned int t =0 ; t < ccqe_weights->size(); t++){
                            h_proton_num.at(0).at(t)->Fill(proton_mom, ccqe_weights->at(t));
                            h_proton_num.at(1).at(t)->Fill(proton_mom, ccmec_weights->at(t));
                            h_proton_num_cv->Fill(proton_mom);
                        }

                    }
                    // Failed the selection so denominator
                    else {
                        for (unsigned int t =0 ; t < ccqe_weights->size(); t++){
                            h_proton_den.at(0).at(t)->Fill(proton_mom, ccqe_weights->at(t));
                            h_proton_den.at(1).at(t)->Fill(proton_mom, ccmec_weights->at(t));
                            h_proton_den_cv->Fill(proton_mom);
                        }

                    }
                } // end if proton

                // Leading electron
                if (j == leading_electron_index){

                    if (electron_energy == 0.0) continue;

                    // Passed the selection so numerator
                    if (tpc_object_pass_temp->at(i) == 1){
                        for (unsigned int t =0 ; t < ccqe_weights->size(); t++){
                            h_electron_num.at(0).at(t)->Fill(electron_energy, ccqe_weights->at(t));
                            h_electron_num.at(1).at(t)->Fill(electron_energy, ccmec_weights->at(t));
                            h_electron_num_cv->Fill(electron_energy);
                        }

                    }
                    // Failed the selection so denominator
                    else {
                        for (unsigned int t =0 ; t < ccqe_weights->size(); t++){
                            h_electron_den.at(0).at(t)->Fill(electron_energy, ccqe_weights->at(t));
                            h_electron_den.at(1).at(t)->Fill(electron_energy, ccmec_weights->at(t));
                            h_electron_den_cv->Fill(electron_energy);
                        }

                    }

                } // end if electron

            } // END loop over pfp

        
        }// End loop over tpc objects


    } // END loop over the entires

   draw_weight_hists(h_proton_num.at(0), 0, 1500, "plots/h_ccqe_proton_num.pdf");
   draw_weight_hists(h_proton_num.at(1), 0, 1500, "plots/h_ccmec_proton_num.pdf");

   draw_weight_hists(h_proton_den.at(0), 0, 1500, "plots/h_ccqe_proton_den.pdf");
   draw_weight_hists(h_proton_den.at(1), 0, 1500, "plots/h_ccmec_proton_den.pdf");

   draw_weight_hists(h_electron_num.at(0), 0, 2000, "plots/h_ccqe_electron_num.pdf");
   draw_weight_hists(h_electron_num.at(1), 0, 2000, "plots/h_ccmec_electron_num.pdf");

   draw_weight_hists(h_electron_den.at(0), 0, 2000, "plots/h_ccqe_electron_den.pdf");
   draw_weight_hists(h_electron_den.at(1), 0, 2000, "plots/h_ccmec_electron_den.pdf");

    // Now lets take the ratio and plot that
    for (unsigned int t =0 ; t < ccqe_weights->size(); t++){
        h_proton_num.at(0).at(t)->Divide(h_proton_den.at(0).at(t));
        h_proton_num.at(1).at(t)->Divide(h_proton_den.at(1).at(t));
        h_electron_num.at(0).at(t)->Divide(h_electron_den.at(0).at(t));
        h_electron_num.at(1).at(t)->Divide(h_electron_den.at(1).at(t));

        h_proton_num.at(0).at(t)->SetTitle(";Proton Momentum [GeV/c];Efficiency");
        h_proton_num.at(1).at(t)->SetTitle(";Proton Momentum [GeV/c];Efficiency");
        
        h_electron_num.at(0).at(t)->SetTitle(";Electron Energy [GeV];Efficiency");
        h_electron_num.at(1).at(t)->SetTitle(";Electron Energy [GeV];Efficiency");
    }

    // h_proton_num_cv  ->Divide(h_proton_cv_den);
    // h_electron_num_cv->Divide(h_electron_cv_den);


    draw_weight_hists(h_proton_num.at(0), 0, 0.3, "plots/h_ccqe_proton_eff.pdf");
    draw_weight_hists(h_proton_num.at(1), 0, 0.3, "plots/h_ccmec_proton_eff.pdf");

    draw_weight_hists(h_electron_num.at(0), 0, 0.2, "plots/h_ccqe_electron_eff.pdf");
    draw_weight_hists(h_electron_num.at(1), 0, 0.2, "plots/h_ccmec_electron_eff.pdf");


    // Divide the by the denominator to get the CV efficiency
    h_E_transfer_CV_num->Divide(h_E_transfer_CV_den);
    h_P_transfer_CV_num->Divide(h_P_transfer_CV_den);
    h_NumCharged_CV_num->Divide(h_NumCharged_CV_den);


    // For the energy transfer histograms we need to normalise them to the CV first
    double int_E_transfer_CV = h_E_transfer_CV_den->Integral();
    double int_P_transfer_CV = h_P_transfer_CV_den->Integral();
    double int_NumCharged_CV = h_NumCharged_CV_den->Integral();
    for (unsigned int t =0 ; t < ccqe_weights->size(); t++){
        
        double E_transfer_v_0_sf = int_E_transfer_CV / h_E_transfer_den_v.at(0).at(t)->Integral();
        double E_transfer_v_1_sf = int_E_transfer_CV / h_E_transfer_den_v.at(1).at(t)->Integral();

        double P_transfer_v_0_sf = int_P_transfer_CV / h_P_transfer_den_v.at(0).at(t)->Integral();
        double P_transfer_v_1_sf = int_P_transfer_CV / h_P_transfer_den_v.at(1).at(t)->Integral();

        double NumCharged_v_0_sf = int_NumCharged_CV / h_NumCharged_den_v.at(0).at(t)->Integral();
        double NumCharged_v_1_sf = int_NumCharged_CV / h_NumCharged_den_v.at(1).at(t)->Integral();
        
        h_E_transfer_den_v.at(0).at(t)->Scale(E_transfer_v_0_sf);
        h_E_transfer_den_v.at(1).at(t)->Scale(E_transfer_v_1_sf);

        h_P_transfer_den_v.at(0).at(t)->Scale(P_transfer_v_0_sf);
        h_P_transfer_den_v.at(1).at(t)->Scale(P_transfer_v_1_sf);

        h_NumCharged_den_v.at(0).at(t)->Scale(NumCharged_v_0_sf);
        h_NumCharged_den_v.at(1).at(t)->Scale(NumCharged_v_1_sf);
    }

    std::vector<std::vector<double>> E_transfer_std;
    std::vector<std::vector<double>> P_transfer_std;
    std::vector<std::vector<double>> NumCharged_std;
    E_transfer_std.resize(2);
    P_transfer_std.resize(2);
    NumCharged_std.resize(2);

    E_transfer_std.at(0).resize(h_E_transfer_CV_den->GetNbinsX()+1,0);
    E_transfer_std.at(1).resize(h_E_transfer_CV_den->GetNbinsX()+1,0);

    P_transfer_std.at(0).resize(h_P_transfer_CV_den->GetNbinsX()+1,0);
    P_transfer_std.at(1).resize(h_P_transfer_CV_den->GetNbinsX()+1,0);

    NumCharged_std.at(0).resize(h_NumCharged_CV_den->GetNbinsX()+1,0);
    NumCharged_std.at(1).resize(h_NumCharged_CV_den->GetNbinsX()+1,0);




    for (unsigned int t =0 ; t < ccqe_weights->size(); t++){
        
        // Energy Transfer
        for (unsigned int j = 1; j < h_E_transfer_CV_den->GetNbinsX()+1; j++){
            E_transfer_std.at(0).at(j - 1) += ( (h_E_transfer_den_v.at(0).at(t)->GetBinContent(j) - h_E_transfer_CV_den->GetBinContent(j)) * (h_E_transfer_den_v.at(0).at(t)->GetBinContent(j) - h_E_transfer_CV_den->GetBinContent(j)) );
            E_transfer_std.at(1).at(j - 1) += ( (h_E_transfer_den_v.at(1).at(t)->GetBinContent(j) - h_E_transfer_CV_den->GetBinContent(j)) * (h_E_transfer_den_v.at(1).at(t)->GetBinContent(j) - h_E_transfer_CV_den->GetBinContent(j)) );

        }

        // P Transfer
        for (unsigned int j = 1; j < h_P_transfer_CV_den->GetNbinsX()+1; j++){
            P_transfer_std.at(0).at(j - 1) += ( (h_P_transfer_den_v.at(0).at(t)->GetBinContent(j) - h_P_transfer_CV_den->GetBinContent(j)) * (h_P_transfer_den_v.at(0).at(t)->GetBinContent(j) - h_P_transfer_CV_den->GetBinContent(j)) );
            P_transfer_std.at(1).at(j - 1) += ( (h_P_transfer_den_v.at(1).at(t)->GetBinContent(j) - h_P_transfer_CV_den->GetBinContent(j)) * (h_P_transfer_den_v.at(1).at(t)->GetBinContent(j) - h_P_transfer_CV_den->GetBinContent(j)) );

        }    

        // Num Charged
        for (unsigned int j = 1; j < h_NumCharged_CV_den->GetNbinsX()+1; j++){
            NumCharged_std.at(0).at(j - 1) += ( (h_NumCharged_den_v.at(0).at(t)->GetBinContent(j) - h_NumCharged_CV_den->GetBinContent(j)) * (h_NumCharged_den_v.at(0).at(t)->GetBinContent(j) - h_NumCharged_CV_den->GetBinContent(j)) );
            NumCharged_std.at(1).at(j - 1) += ( (h_NumCharged_den_v.at(1).at(t)->GetBinContent(j) - h_NumCharged_CV_den->GetBinContent(j)) * (h_NumCharged_den_v.at(1).at(t)->GetBinContent(j) - h_NumCharged_CV_den->GetBinContent(j)) );

        }    
    }

    // Energy Transfer
    for (unsigned int j = 1; j < h_E_transfer_CV_den->GetNbinsX()+1; j++){
        E_transfer_std.at(0).at(j - 1) = std::sqrt(E_transfer_std.at(0).at(j - 1)/1000);
        E_transfer_std.at(1).at(j - 1) = std::sqrt(E_transfer_std.at(1).at(j - 1)/1000);

        h_E_transfer_CV_den->SetBinError(j, E_transfer_std.at(0).at(j - 1) );

    }

    

    // P Transfer
    for (unsigned int j = 1; j < h_P_transfer_CV_den->GetNbinsX()+1; j++){
        P_transfer_std.at(0).at(j - 1) = std::sqrt(P_transfer_std.at(0).at(j - 1)/1000);
        P_transfer_std.at(1).at(j - 1) = std::sqrt(P_transfer_std.at(1).at(j - 1)/1000);

        h_P_transfer_CV_den->SetBinError(j, P_transfer_std.at(0).at(j - 1) );

    } 

    // Energy Transfer
    for (unsigned int j = 1; j < h_NumCharged_CV_den->GetNbinsX()+1; j++){
        NumCharged_std.at(0).at(j - 1) = std::sqrt(NumCharged_std.at(0).at(j - 1)/1000);
        NumCharged_std.at(1).at(j - 1) = std::sqrt(NumCharged_std.at(1).at(j - 1)/1000);

        h_NumCharged_CV_den->SetBinError(j, NumCharged_std.at(0).at(j - 1) );

    } 


    draw_standard_hists(h_E_transfer_CV_den, h_E_transfer_CV_num, "plots/h_E_transfer_CV_den_ccqe.pdf");
    draw_standard_hists(h_P_transfer_CV_den, h_P_transfer_CV_num, "plots/h_P_transfer_CV_den_ccqe.pdf");
    draw_standard_hists(h_NumCharged_CV_den, h_NumCharged_CV_num, "plots/h_NumCharged_CV_den_ccqe.pdf");


    for (unsigned int j = 1; j < h_E_transfer_CV_den->GetNbinsX()+1; j++){
        h_E_transfer_CV_den->SetBinError(j, E_transfer_std.at(1).at(j - 1) );
    }

    for (unsigned int j = 1; j < h_P_transfer_CV_den->GetNbinsX()+1; j++){
        h_P_transfer_CV_den->SetBinError(j, P_transfer_std.at(1).at(j - 1) );
    }

    for (unsigned int j = 1; j < h_NumCharged_CV_den->GetNbinsX()+1; j++){
        h_NumCharged_CV_den->SetBinError(j, NumCharged_std.at(1).at(j - 1) );
    }

    draw_standard_hists(h_E_transfer_CV_den, h_E_transfer_CV_num, "plots/h_E_transfer_CV_den_ccmec.pdf");
    draw_standard_hists(h_P_transfer_CV_den, h_P_transfer_CV_num, "plots/h_P_transfer_CV_den_ccmec.pdf");
    draw_standard_hists(h_NumCharged_CV_den, h_NumCharged_CV_num, "plots/h_NumCharged_CV_den_ccmec.pdf");

}
//***************************************************************************
//***************************************************************************
void mc_truth::draw_weight_hists(std::vector<TH1D*> hist_v, double y1, double y2, const char*  printname){


     TCanvas *c = new TCanvas();

    for (unsigned int t =0 ; t < hist_v.size(); t++){
        hist_v.at(t)->GetYaxis()->SetRangeUser(y1, y2);
        hist_v.at(t)->Draw("hist,same");
    }

    c->Print(printname);

}
//***************************************************************************
//***************************************************************************
void mc_truth::draw_standard_hists(TH1D* hist, TH1D* h_eff, const char*  printname){

    TH1D *h_eff_clone = (TH1D *)h_eff->Clone("h_clone");


    TCanvas *c = new TCanvas();
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(2);
    hist->Draw("hist,E");

    // Now draw the efficiency on top
    gPad->SetRightMargin(0.17);

    c->Update();
    
    Float_t rightmax = 1.1 * h_eff_clone->GetMaximum();
    Float_t scale = gPad->GetUymax() / rightmax;
    h_eff_clone->SetLineColor(kAzure - 6);
    h_eff_clone->SetLineWidth(2);
    h_eff_clone->Scale(scale);
    h_eff_clone->Draw("hist,same");

    

    TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, rightmax, 510, "+L");
    axis->SetTitle("Efficiency");
    axis->SetTitleOffset(1.1);
    axis->SetLineColor(kAzure - 6);
    axis->SetLabelColor(kAzure - 6);
    axis->SetTitleColor(kAzure - 6);
    axis->SetTextFont(42);
    axis->SetLabelFont(42);

    axis->Draw();

    hist->Draw("hist,E,same");
    
    
    
    c->Print(printname);

    delete c;
    delete h_eff_clone;

}