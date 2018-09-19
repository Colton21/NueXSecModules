#include "test.h"
//This way you can be lazy
using namespace art;
using namespace std;


int main(int argv, char** argc) {
  
  double fDetLength = 1036.8; 
  double fDetHalfWidth = 128.175;
  double fDetHalfHeight = 116.5;

  std::pair<float, float>  _xRange;
  std::pair<float, float>  _yRange;
  std::pair<float, float>  _zRange;

  _xRange.first  = 0;
  _xRange.second = 2*fDetHalfWidth;
  _yRange.first  = -1*fDetHalfHeight;
  _yRange.second = fDetHalfHeight;
  _zRange.first =  0;
  _zRange.second = fDetLength;

  //geoalgo::GeoAlgo const _geo_algo_instance;

  //geoalgo::AABox volAVTPC( _xRange.first, _yRange.first, _zRange.first, _xRange.second, _yRange.second, _zRange.second);

  /* Here you name the thing that "produced" the data product that you want to look at.
     In our event dump we see two things:
     GenieGen.. | corsika............. | .... | std::vector<simb::MCTruth>.... | ....1
     GenieGen.. | generator........... | .... | std::vector<simb::MCTruth>.... | ....1

     This means that if you want to look at the simb::MCTruth data product there are two "producers"
     We want to look at neutrinos so we choose "generator" if you wanted to look at cosmics you
     could pick "coriska"
  */
  InputTag mctruths_tag { "flux" };
  InputTag  evtwght_tag { "eventweight" };
 
  //We have passed the input file as an argument to the function 
  vector<string> filename;
  for (int i = 1; i < argv; ++i) { 
    std::cout << "FILE : " << argc[i] << std::endl; 
    filename.push_back(string(argc[i]));
  }

  //Define the histograms you want here

  //flavors
  std::vector< TH1D* > Enu_CV_Window;
  std::vector< TH1D* > Enu_CV_AV_TPC;
  //flavors - systematic - universe 
  std::vector< std::vector < std::vector < TH1D* > > > Enu_Syst_Window;
  std::vector< std::vector < std::vector < TH1D* > > > Enu_Syst_AV_TPC;

  std::vector< double > TotWeight;
  TotWeight.resize(1000);

  //systematic - universe 
  std::vector< std::vector< double > > Weights;   
  Weights.resize(6);
  for(int i = 0; i < 6; i++){
    Weights[i].resize(1000);
  }

  //
   
  std::vector< string > flav;
  flav = {"numu","nue","numubar","nuebar"};
  
  std::vector< string > labels;
  labels = {"FluxUnisim", "piplus", "piminus",
            "kplus", "kminus", "kzero","total"};



  Enu_CV_Window.resize(4);
  Enu_CV_AV_TPC.resize(4);
  Enu_Syst_Window.resize(4);
  Enu_Syst_AV_TPC.resize(4);

  for(int i = 0; i < int(flav.size()); i++){//flav

    Enu_CV_Window[i] = new TH1D(Form("%s_CV_Window",flav[i].c_str()),"", 200, 0, 10);
    Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"", 200, 0, 10);

    Enu_Syst_Window[i].resize(labels.size());
    Enu_Syst_AV_TPC[i].resize(labels.size());

    for(int j = 0; j < int(labels.size()); j++){//labels
      
      Enu_Syst_Window[i][j].resize(1000);
      Enu_Syst_AV_TPC[i][j].resize(1000);
     
      
      for(int k = 0; k < 1000; k++){//unis
	Enu_Syst_Window[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_Window",flav[i].c_str(), labels[j].c_str(), k),"", 200, 0, 10);
	Enu_Syst_AV_TPC[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_AV_TPC",flav[i].c_str(), labels[j].c_str(), k),"", 200, 0, 10);
      }//iterate over universes
    }//iterate over systematics
  }//iterate over flavors 

  //Let's Do Science! 
      int y = 0;
      int n = 0;

  for (gallery::Event ev(filename) ; !ev.atEnd(); ev.next()) {

    //Next we want to grab from the event the data-product that you want
    auto const& mctruths = *ev.getValidHandle<vector<simb::MCTruth>>(mctruths_tag);
    auto const& evtwghts = *ev.getValidHandle<vector<evwgh::MCEventWeight>>(evtwght_tag);

    if (mctruths.empty() || evtwghts.empty()) continue;
    
    
    //Now we'll iterate through these 
    for (size_t i = 0; i < mctruths.size(); i++) {
      auto const& mctruth = mctruths.at(i);
      auto const& evtwght = evtwghts.at(i);


      //mctruth.GetNeutrino().Nu().PdgCode()
      //mctruth.GetNeutrino().Nu().E()

      //mctruth.GetNeutrino().Nu().Vx()
      //mctruth.GetNeutrino().Nu().Vy()
      //mctruth.GetNeutrino().Nu().Vz()

      //mctruth.GetNeutrino().Nu().Px()
      //mctruth.GetNeutrino().Nu().Py()
      //mctruth.GetNeutrino().Nu().Pz()

      int pdg = 9; 
      if(mctruth.GetNeutrino().Nu().PdgCode() == 14) {pdg = 0;}
      else if(mctruth.GetNeutrino().Nu().PdgCode() == 12) {pdg = 1;}
      else if(mctruth.GetNeutrino().Nu().PdgCode() ==-14) {pdg = 2;}
      else if(mctruth.GetNeutrino().Nu().PdgCode() ==-12) {pdg = 3;}
      else {std::cout << "shits cray " << mctruth.GetNeutrino().Nu().PdgCode() << std::endl;}

      /*
      geoalgo::HalfLine ray(mctruth.GetNeutrino().Nu().Vx()*100,
			    mctruth.GetNeutrino().Nu().Vy()*100,
			    mctruth.GetNeutrino().Nu().Vz()*100,
			    mctruth.GetNeutrino().Nu().Px(),
			    mctruth.GetNeutrino().Nu().Py(),
			    mctruth.GetNeutrino().Nu().Pz());
      */      

      //geoalgo::HalfLine ray(_xRange.second/2, fDetHalfHeight-10, -1000, 0, 0, -1);

      //auto vec = _geo_algo_instance.Intersection(volAVTPC, ray);
     
      //      std::cout << "# of intersections : " << vec.size()  << std::endl;

      //bool intercept = false;

      //if(vec.size() == 0){ intercept = false; }
      //if(vec.size() == 2){ intercept = true; }
      //if(vec.size() != 2 && vec.size() != 0){ std::cout << "you dum" << std::endl;}

      double bnbweight = 1; 
      for(int l = 0; l < int(labels.size())-1; l++){
	std::fill(Weights[l].begin(), Weights[l].end(), 1);
      }

      //      std::cout << "Getting weights!" << std::endl; 
      for(auto last : evtwght.fWeight){
        if(last.first.find("bnbcorrection") != std::string::npos)
	  bnbweight = last.second.at(0);       
      }
      //      std::cout << "BNB Correction Weight : " << bnbweight << std::endl;
/*
      for(auto last : evtwght.fWeight){
        for(int l = 0; l < int(labels.size())-1; l++){
          if(last.first.find(labels[l].c_str()) != std::string::npos){
            for(int i = 0; i < int(last.second.size()); i++){	      
              Weights[l][i] *= last.second.at(i);
	    }
	  }
	}
      }
*/
      //      std::cout << "got all my weights! " << std::endl; 

      Enu_CV_Window[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), bnbweight);
      /*
      if(intercept){
	Enu_CV_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), bnbweight);
      }
      */
      std::fill(TotWeight.begin(), TotWeight.end(), 1);

      /*

      for(int l = 0; l < int(labels.size())-1; l++){
        for(unsigned int i = 0; i < Weights[l].size(); i++){
	  TotWeight[i] *= Weights[l][i];
	  Enu_Syst_Window[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(), Weights[l][i]*bnbweight);

	  if(intercept){	  
	    Enu_Syst_AV_TPC[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(), Weights[l][i]*bnbweight);
	  }
	}
      }

      int full_lable = labels.size()-1;

      for(unsigned int i = 0; i < Weights[0].size(); i++){
	Enu_Syst_Window[pdg][full_lable][i]->Fill(mctruth.GetNeutrino().Nu().E(), TotWeight[i]*bnbweight);
	
	if(intercept){	  
	  Enu_Syst_AV_TPC[pdg][full_lable][i]->Fill(mctruth.GetNeutrino().Nu().E(), TotWeight[i]*bnbweight);
	}	
      }
      */   
         
    }//Iterate through neutrino interactions
  }// Iterate through events

  //Drawing the histograms don't work like they do with macros

  /*

  //We can also write the output to a output root file like this:
  TFile* output = new TFile("output.root","RECREATE");

  TDirectory *savdir = gDirectory;

  std::vector< std::vector< std::vector< TDirectory* > > > subdir; //flav //syst //cont 
  subdir.resize(4);
  for(int i = 0; i < 4; i++){
    subdir[i].resize(8);
    for(int j = 0; j < 8; j++){
      subdir[i][j].resize(3);
    }
  }
  */

  //std::vector<string> cont; cont = {"Window","Active_TPC_Volume",""};

  /*for(int f = 0; f < 4; f++){
    std::cout << flav[f] << std::endl;
    subdir[f][0][0] = savdir->mkdir(Form("%s",flav[f].c_str()));
    subdir[f][0][0]->cd();
    
    Enu_CV_Window[f]->Write();      
    Enu_CV_AV_TPC[f]->Write();

    for(int s = 1; s < 8; s++){
      std::cout << labels[s-1] << std::endl;
      subdir[f][s][0] = subdir[f][0][0]->mkdir(Form("%s",labels[s-1].c_str()));
      subdir[f][s][0]->cd();
      for(int c = 1; c < 3; c++){
	std::cout << cont[c-1] << std::endl;
	subdir[f][s][c] = subdir[f][s][0]->mkdir(Form("%s",cont[c-1].c_str()));
	subdir[f][s][c]->cd();
	
	if(c == 1){
	  for(int i = 0; i < 1000; i++){
	    Enu_Syst_Window[f][s-1][i]->Write();
	  }	
	}

	if(c == 2){
	  for(int i = 0; i < 1000; i++){
	    Enu_Syst_AV_TPC[f][s-1][i]->Write();
	  }	
	}

      }//cont
      
    }//systs
    savdir->cd();
  }//flavs

  output->Close();
  */
  return 1;

    
}//main 

    /*
  TDirectory *numudir = savdir->mkdir("numu");
  numudir->cd();

  Enu_CV_Window[0]->Write();
  Enu_CV_AV_TPC[0]->Write();

  TDirectory *numusystdir = numudir->mkdir("syst");
  numusystdir->cd();

  TDirectory *numusyst1dir = numusystdir->mkdir(Form("numu_%s",labels[0].c_str()));
  numusyst1dir->cd();

  for(int i = 0; i < Enu_Syst_Window[pdg][l][i])

  numusystdir->cd();
  TDirectory *numusyst2dir = numusystdir->mkdir(Form("numu_%s",labels[1].c_str()));
  numusyst2dir->cd();

  numusystdir->cd();
  TDirectory *numusyst3dir = numusystdir->mkdir(Form("numu_%s",labels[2].c_str()));
  numusyst3dir->cd();

  numusystdir->cd();
  TDirectory *numusyst4dir = numusystdir->mkdir(Form("numu_%s",labels[3].c_str()));
  numusyst4dir->cd();

  numusystdir->cd();
  TDirectory *numusyst5dir = numusystdir->mkdir(Form("numu_%s",labels[4].c_str()));
  numusyst5dir->cd();

  numusystdir->cd();
  TDirectory *numusyst6dir = numusystdir->mkdir(Form("numu_%s",labels[5].c_str()));
  numusyst6dir->cd();




  savdir->cd();
  TDirectory *nuedir = savdir->mkdir("nue");
  nuedir->cd();

  Enu_CV_Window[1]->Write();
  Enu_CV_AV_TPC[1]->Write();

  TDirectory *nuesystdir = nuedir->mkdir("syst");
  nuesystdir->cd();


  savdir->cd();
  TDirectory *numubardir = savdir->mkdir("numubar");
  numubardir->cd();

  Enu_CV_Window[2]->Write();
  Enu_CV_AV_TPC[2]->Write();

  TDirectory *numubarsystdir = numubardir->mkdir("syst");
  numubarsystdir->cd();

  savdir->cd();
  TDirectory *nuebardir = savdir->mkdir("nuebar");
  nuebardir->cd();

  Enu_CV_Window[3]->Write();
  Enu_CV_AV_TPC[3]->Write();

  TDirectory *nuebarsystdir = nuebardir->mkdir("syst");
  nuebarsystdir->cd();
    

  output->Close();
  
  return 1;
  }*/
