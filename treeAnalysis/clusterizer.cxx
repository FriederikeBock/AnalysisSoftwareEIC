#include <algorithm>
// ANCHOR debug output verbosity
int verbosityCLS = 1;

// ANCHOR define global variables
bool _do_V1clusterizer = false;
bool _do_V3clusterizer = false;
bool _do_XNclusterizer = false;
bool _do_NxNclusterizer = false;

bool _do_V1clusterizerFEMC = false;
bool _do_V3clusterizerFEMC = false;
bool _do_XNclusterizerFEMC = false;
bool _do_NxNclusterizerFEMC = false;

TH2F*  h_clusterizer_matched_dx_dy[5][4]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_all_dx_dy[5][4]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspec_matched_E_eta[5][4]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_clusterizer_clsspec_PDG[5][4]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_clusterizer_clsspec_matched_PDG[5][4]; // [calorimeter_enum][algorithm_enum]

TString str_CLSRZR_mcparticles[5] = {"electron", "cpion", "proton", "ckaon", "muon"};
int int_CLSRZR_mcparticles_PDG[5] = {11 /*e*/,211  /*pi*/,  2212/*p*/,  321/*K*/,  13/*mu*/};

TH2F*  h_clusterizer_clsspec_E_eta[5][4]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecMC_E_eta[5][4]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspec_particle_E_eta[5][4][5]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecMC_particle_E_eta[5][4][5]; // [calorimeter_enum][algorithm_enum]

bool _doCalibration = false;

float calibrationFEMC(float clusterE){
  float paramsECALcalib[5]= {3.15612e+00,1.57658e-01,1.66441e+00,-2.50979e+02,3.81511e+02}; //noshift
  return (( paramsECALcalib[0] + paramsECALcalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramsECALcalib[2] * TMath::Exp( ( clusterE - paramsECALcalib[3] ) / paramsECALcalib[4] ) ) ));
}

float calibrationFHCAL(float clusterE){
  float paramsFHCALcalib[5]= {1.75773e+00, 3.87923e-01, 1.20862e+00, -1.25194e+02, 7.54438e+01}; //noshift
  return (( paramsFHCALcalib[0] + paramsFHCALcalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramsFHCALcalib[2] * TMath::Exp( ( clusterE - paramsFHCALcalib[3] ) / paramsFHCALcalib[4] ) ) ));
}

bool isClusterMatched(int clsID, float matchingwindow, float* clusters_X, float* clusters_Y, float* clusters_E, int caloEnum, int clusterizerEnum, bool useProjection = true);

// ANCHOR main function to be called in event loop
void runclusterizer(
  int clusterizerEnum,
  int caloEnum,
  float seedE,
  float aggE,
  int &nclusters,
  float* &clusters_E,
  float* &clusters_Eta,
  float* &clusters_Phi,
  float* &clusters_M02,
  float* &clusters_M20,
  bool* &clusters_isMatched,
  int* &clusters_NTower,
  int* &clusters_trueID,
  int* &clusters_NtrueID,
  float* &clusters_X,
  float* &clusters_Y,
  float* &clusters_Z
){



  nclusters = 0;
  if(verbosityCLS>1)cout << "XN: new event" << endl;

  for(int icalo=0;icalo<5;icalo++){
    for(int ialgo=0;ialgo<4;ialgo++){
      if(!h_clusterizer_clsspec_PDG[icalo][ialgo])h_clusterizer_clsspec_PDG[icalo][ialgo] 	= new TH1F(Form("h_clusterizer_clsspec_PDG_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 2500,0,2500);
      if(!h_clusterizer_clsspec_matched_PDG[icalo][ialgo])h_clusterizer_clsspec_matched_PDG[icalo][ialgo] 	= new TH1F(Form("h_clusterizer_clsspec_matched_PDG_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 2500,0,2500);
      if(!h_clusterizer_clsspec_matched_E_eta[icalo][ialgo])h_clusterizer_clsspec_matched_E_eta[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_clsspec_matched_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, 80, -4.0,4.0);

      if(!h_clusterizer_clsspec_E_eta[icalo][ialgo])h_clusterizer_clsspec_E_eta[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_clsspec_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      if(!h_clusterizer_clsspecMC_E_eta[icalo][ialgo])h_clusterizer_clsspecMC_E_eta[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_clsspecMC_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      for(int ipart=0;ipart<5;ipart++){
        if(!h_clusterizer_clsspec_particle_E_eta[icalo][ialgo][ipart])h_clusterizer_clsspec_particle_E_eta[icalo][ialgo][ipart] 	= new TH2F(Form("h_clusterizer_clsspec_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CLSRZR_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
        if(!h_clusterizer_clsspecMC_particle_E_eta[icalo][ialgo][ipart])h_clusterizer_clsspecMC_particle_E_eta[icalo][ialgo][ipart] 	= new TH2F(Form("h_clusterizer_clsspecMC_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CLSRZR_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      }
    }
  }

  // vector of towers used in the clusterizer
  std::vector<towersStrct> input_towers;
  // vector of towers within the currently found cluster
  std::vector<towersStrct> cluster_towers;

  // fill vector with towers for clusterization above aggregation threshold
  if(caloEnum==kFHCAL){
    for(int itow=0; itow<_nTowers_FHCAL; itow++){
      if(_tower_FHCAL_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E = _tower_FHCAL_E[itow];
        tempstructT.tower_iEta = _tower_FHCAL_iEta[itow];
        tempstructT.tower_iPhi = _tower_FHCAL_iPhi[itow];
        tempstructT.tower_trueID = _tower_FHCAL_trueID[itow];
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kFEMC){
    for(int itow=0; itow<_nTowers_FEMC; itow++){
      if(_tower_FEMC_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E = _tower_FEMC_E[itow];
        tempstructT.tower_iEta = _tower_FEMC_iEta[itow];
        tempstructT.tower_iPhi = _tower_FEMC_iPhi[itow];
        tempstructT.tower_trueID = _tower_FEMC_trueID[itow];
        input_towers.push_back(tempstructT);
      }
    }
  } else {
      cout << "Incorrect calorimeter selected! Enum " << clusterizerEnum << " not defined!" << endl;
      return;
    }
  // sort vector in descending energy order
  std::sort(input_towers.begin(), input_towers.end(), &acompare);
  std::vector<int> clslabels;
  while (!input_towers.empty()) {
    cluster_towers.clear();
    clslabels.clear();
    // always start with highest energetic tower
    if(input_towers.at(0).tower_E > seedE){
      // fill seed cell information into current cluster
      clusters_E[nclusters] = input_towers.at(0).tower_E;
      clusters_NTower[nclusters] = 1;
      clusters_NtrueID[nclusters] = 1;
      clusters_trueID[nclusters] = input_towers.at(0).tower_trueID; // TODO save all MC labels?
      cluster_towers.push_back(input_towers.at(0));
      clslabels.push_back(input_towers.at(0).tower_trueID);

      // ANCHOR XN clusterizer logic
      if(clusterizerEnum==kXN){
        for (int tit = 1; tit < input_towers.size(); tit++){
          // towers must be within 3x3 matrix (delta Eta and delta Phi < 2)
          if( ( std::abs(input_towers.at(tit).tower_iEta-input_towers.at(0).tower_iEta) + std::abs(input_towers.at(tit).tower_iPhi-input_towers.at(0).tower_iPhi) ) <= 2 ){
              clusters_E[nclusters]+=input_towers.at(tit).tower_E;
              clusters_NTower[nclusters]++;
              cluster_towers.push_back(input_towers.at(tit));
              if(!(std::find(clslabels.begin(), clslabels.end(), input_towers.at(tit).tower_trueID) != clslabels.end())){
                clusters_NtrueID[nclusters]++;
                clslabels.push_back(input_towers.at(tit).tower_trueID);
              }
              input_towers.erase(input_towers.begin()+tit);
          }
        }
      }
      // ANCHOR NxN clusterizer logic
      else if(clusterizerEnum==kNxN){
        for (int tit = 1; tit < input_towers.size(); tit++){
          // towers must be within 3x3 matrix (delta Eta and delta Phi < 2)
          if(std::abs(input_towers.at(tit).tower_iEta-input_towers.at(0).tower_iEta)<2){
            if(std::abs(input_towers.at(tit).tower_iPhi-input_towers.at(0).tower_iPhi)<2){
              clusters_E[nclusters]+=input_towers.at(tit).tower_E;
              clusters_NTower[nclusters]++;
              cluster_towers.push_back(input_towers.at(tit));
              if(!(std::find(clslabels.begin(), clslabels.end(), input_towers.at(tit).tower_trueID) != clslabels.end())){
                clusters_NtrueID[nclusters]++;
                clslabels.push_back(input_towers.at(tit).tower_trueID);
              }
              input_towers.erase(input_towers.begin()+tit);
            }
          }
        }
      }

      // ANCHOR V3 clusterizer logic
      else if(clusterizerEnum==kV3){
        // remove seed tower from sample
        input_towers.erase(input_towers.begin());
        for (int tit = 0; tit < cluster_towers.size(); tit++){
          // Now go recursively to the next 4 neighbours and add them to the cluster if they fulfill the conditions
          int iEtaTwr = cluster_towers.at(tit).tower_iEta;
          int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
          for (int ait = 0; ait < input_towers.size(); ait++){
            int iEtaTwrAgg = input_towers.at(ait).tower_iEta;
            int iPhiTwrAgg = input_towers.at(ait).tower_iPhi;
            if( (TMath::Abs(iEtaTwrAgg-iEtaTwr)+TMath::Abs(iPhiTwrAgg-iPhiTwr)) == 1){
              // only aggregate towers with lower energy than current tower
              if(input_towers.at(ait).tower_E >= cluster_towers.at(tit).tower_E) continue;
              clusters_E[nclusters]+=input_towers.at(ait).tower_E;
              clusters_NTower[nclusters]++;
              cluster_towers.push_back(input_towers.at(ait));
              if(!(std::find(clslabels.begin(), clslabels.end(), input_towers.at(ait).tower_trueID) != clslabels.end())){
                clusters_NtrueID[nclusters]++;
                clslabels.push_back(input_towers.at(ait).tower_trueID);
              }
              input_towers.erase(input_towers.begin()+ait);
            }
          }
        }
      }
      else {
        cout << "incorrect clusterizer selected! Enum " << clusterizerEnum << " not defined!" << endl;
        return;
      }

      // determine remaining cluster properties from its towers
      float* showershape_eta_phi = CalculateM02andWeightedPosition(cluster_towers, weightM02, clusters_E[nclusters], caloEnum, kFALSE);
      clusters_M02[nclusters] = showershape_eta_phi[0];
      clusters_M20[nclusters] = showershape_eta_phi[1];
      clusters_Eta[nclusters] = showershape_eta_phi[2];
      clusters_Phi[nclusters] = showershape_eta_phi[3];
      clusters_X[nclusters] = showershape_eta_phi[4];
      clusters_Y[nclusters] = showershape_eta_phi[5];
      clusters_Z[nclusters] = showershape_eta_phi[6];
      clusters_isMatched[nclusters] = isClusterMatched(nclusters, 20,clusters_X,clusters_Y,clusters_E, caloEnum, clusterizerEnum, true);
     // if(verbosityCLS>1) cout << "\tXN cluster with E = " << clusters_E[nclusters] << "\tEta: " << clusters_Eta[nclusters]<< "\tPhi: " << clusters_Phi[nclusters]<< "\tntowers: " << clusters_NTower[nclusters] << "\ttrueID: " << clusters_trueID[nclusters] << endl;
      // remove clusterized towers
      if(!(clusterizerEnum==kV3)){
        input_towers.erase(input_towers.begin());
      }

      // apply calibration if desired
      if(_doCalibration){
        if(caloEnum==kFHCAL)
          clusters_E[nclusters]/=calibrationFHCAL(clusters_E[nclusters]);
        else if(caloEnum==kFEMC)
          clusters_E[nclusters]/=calibrationFEMC(clusters_E[nclusters]);
      }

      // fill inputs for efficiency
      h_clusterizer_clsspec_E_eta[caloEnum][clusterizerEnum]->Fill(clusters_E[nclusters], clusters_Eta[nclusters]);
      h_clusterizer_clsspecMC_E_eta[caloEnum][clusterizerEnum]->Fill(_mcpart_E[clusters_trueID[nclusters]-1], _mcpart_Eta[clusters_trueID[nclusters]-1]);

      for(int ipart=0;ipart<5;ipart++){
        if(int_CLSRZR_mcparticles_PDG[ipart]==_mcpart_PDG[clusters_trueID[nclusters]-1]){
          h_clusterizer_clsspec_particle_E_eta[caloEnum][clusterizerEnum][ipart]->Fill(clusters_E[nclusters], clusters_Eta[nclusters]);
          h_clusterizer_clsspecMC_particle_E_eta[caloEnum][clusterizerEnum][ipart]->Fill(_mcpart_E[clusters_trueID[nclusters]-1], _mcpart_Eta[clusters_trueID[nclusters]-1]);
        }
      }

      if(clusters_isMatched[nclusters]) h_clusterizer_clsspec_matched_E_eta[caloEnum][clusterizerEnum]->Fill(clusters_E[nclusters], clusters_Eta[nclusters]);
      h_clusterizer_clsspec_PDG[caloEnum][clusterizerEnum]->Fill(_mcpart_PDG[clusters_trueID[nclusters]-1]);
      if(clusters_isMatched[nclusters]) h_clusterizer_clsspec_matched_PDG[caloEnum][clusterizerEnum]->Fill(_mcpart_PDG[clusters_trueID[nclusters]-1]);

      nclusters++;
    } else {
      input_towers.clear();
    }
  }
}


// ANCHOR track/projection matching function
// TODO at some point we might want E/p
bool isClusterMatched(int clsID, float matchingwindow, float* clusters_X, float* clusters_Y, float* clusters_E, int caloEnum, int clusterizerEnum, bool useProjection = true){
  for(int icalo=0;icalo<5;icalo++){
    for(int ialgo=0;ialgo<4;ialgo++){
      if(!h_clusterizer_all_dx_dy[icalo][ialgo])h_clusterizer_all_dx_dy[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_all_dx_dy_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 100,-100,100, 100,-100,100);
      if(!h_clusterizer_matched_dx_dy[icalo][ialgo])h_clusterizer_matched_dx_dy[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_matched_dx_dy_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 100,-100,100, 100,-100,100);
    }
  }
  if(useProjection){
    int projectionlayer = 0;
    if(caloEnum==kFHCAL)projectionlayer = 2;
    else if(caloEnum==kFEMC)projectionlayer = 1;
    else {
      cout << "isClusterMatched: caloEnum noch defined, returning FALSE" << endl;
      return false;
    }
    for(Int_t iproj=0; iproj<_nProjections; iproj++){
      if(_track_ProjLayer[iproj]!=projectionlayer) continue;
      if(_track_Proj_t[iproj]<-9000) continue;
      if(_track_Proj_x[iproj]==0 && _track_Proj_y[iproj]==0) continue;
      if(clusters_X[iproj]==0 && clusters_Y[iproj]==0) continue;
      // check eta difference
      h_clusterizer_all_dx_dy[caloEnum][clusterizerEnum]->Fill(_track_Proj_x[iproj]-clusters_X[clsID],_track_Proj_y[iproj]-clusters_Y[clsID]);

      if(abs(_track_Proj_x[iproj]-clusters_X[clsID]) < matchingwindow){
        if(verbosityCLS>1) cout << "\tProj matched in x! prj: " << _track_Proj_x[iproj] << "\tcls: " << clusters_X[clsID] << "\tdelta: " << abs(_track_Proj_x[iproj]-clusters_X[clsID]) << endl;
        // check phi difference
        if(abs(_track_Proj_y[iproj]-clusters_Y[clsID]) < matchingwindow){
          if(verbosityCLS>1) cout << "\tProj matched in y! prj: " << _track_Proj_y[iproj] << "\tcls: " << clusters_Y[clsID] << "\tdelta: " << abs(_track_Proj_y[iproj]-clusters_Y[clsID]) << endl;
          h_clusterizer_matched_dx_dy[caloEnum][clusterizerEnum]->Fill(_track_Proj_x[iproj]-clusters_X[clsID],_track_Proj_y[iproj]-clusters_Y[clsID]);
          return true;
        }
      }
    }
  } else {
    float zHC = 400;
    if(caloEnum==kFHCAL)zHC=350;
    if(caloEnum==kFEMC)zHC=291;
    for(Int_t itrk=0; itrk<_nTracks; itrk++){
      // check eta difference
      TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
      trackvec*=(zHC/trackvec.Z());

      h_clusterizer_all_dx_dy[caloEnum][clusterizerEnum]->Fill(trackvec.X()-clusters_X[clsID],trackvec.Y()-clusters_Y[clsID]);

      if(abs(trackvec.X()-clusters_X[clsID]) < matchingwindow){
        // if(verbosityCLS>1) cout << "\tTrck matched in x! trk: " << trackvec.X() << "\tcls: " << clusters_X[clsID] << "\tdelta: " << abs(trackvec.X()-clusters_X[clsID]) << endl;
        float trkphi = (trackvec.Phi()<0 ? trackvec.Phi()+TMath::Pi() : trackvec.Phi()-TMath::Pi());
        // check phi difference
        if(abs(trackvec.Y()-clusters_Y[clsID]) < matchingwindow){
          h_clusterizer_matched_dx_dy[caloEnum][clusterizerEnum]->Fill(trackvec.X()-clusters_X[clsID],trackvec.Y()-clusters_Y[clsID]);
          // if(verbosityCLS>1)  cout << "\tTrck matched in y! prj: " << trackvec.Y() << "\tcls: " << clusters_Y[clsID] << "\tdelta: " << abs(trackvec.Y()-clusters_Y[clsID]) << endl;
          return true;
        }
      }
    }
  }
  return false;
}


void clusterizerSave(){
// make output directory
    gSystem->Exec("mkdir -p "+outputDir + "/clusterizer");

  // define output file
  TFile* fileOutput = new TFile(Form("%s/clusterizer/output_CLSIZER.root",outputDir.Data()),"RECREATE");

  for(int icalo=0;icalo<5;icalo++){
    for(int ialgo=0;ialgo<4;ialgo++){
      if(h_clusterizer_all_dx_dy[icalo][ialgo]) h_clusterizer_all_dx_dy[icalo][ialgo]->Write();
      if(h_clusterizer_matched_dx_dy[icalo][ialgo]) h_clusterizer_matched_dx_dy[icalo][ialgo]->Write();
      if(h_clusterizer_clsspec_matched_E_eta[icalo][ialgo]) h_clusterizer_clsspec_matched_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspec_PDG[icalo][ialgo]) h_clusterizer_clsspec_PDG[icalo][ialgo]->Write();
      if(h_clusterizer_clsspec_matched_PDG[icalo][ialgo]) h_clusterizer_clsspec_matched_PDG[icalo][ialgo]->Write();
      if(h_clusterizer_clsspec_E_eta[icalo][ialgo]) h_clusterizer_clsspec_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspecMC_E_eta[icalo][ialgo]) h_clusterizer_clsspecMC_E_eta[icalo][ialgo]->Write();
      for(int ipart=0;ipart<5;ipart++){
        if(h_clusterizer_clsspec_particle_E_eta[icalo][ialgo][ipart])h_clusterizer_clsspec_particle_E_eta[icalo][ialgo][ipart]->Write();
        if(h_clusterizer_clsspecMC_particle_E_eta[icalo][ialgo][ipart])h_clusterizer_clsspecMC_particle_E_eta[icalo][ialgo][ipart]->Write();
      }
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
