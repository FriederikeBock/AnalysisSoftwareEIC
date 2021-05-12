#include <algorithm>
// ANCHOR debug output verbosity
int verbosityCLS = 1;

// ANCHOR define global variables
bool _do_V1clusterizer = false;
bool _do_V3clusterizer = false;
float aggregation_margin_V3 = 0.03;
bool _do_MAclusterizer = false;
float aggregation_margin_MA = 0.03;
bool _do_C3clusterizer = false;
bool _do_C5clusterizer = false;
bool _do_3x3clusterizer = false;
bool _do_5x5clusterizer = false;

bool _do_V1clusterizerFEMC = false;
bool _do_V3clusterizerFEMC = false;
bool _do_MAclusterizerFEMC = false;
bool _do_C3clusterizerFEMC = false;
bool _do_C5clusterizerFEMC = false;
bool _do_3x3clusterizerFEMC = false;
bool _do_5x5clusterizerFEMC = false;

TH2F*  h_clusterizer_nonagg_towers[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_matched_dx_dy[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_all_dx_dy[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
bool _doClusterECalibration = true;

bool isClusterMatched(int clsID, float matchingwindow, float* clusters_X, float* clusters_Y, float* clusters_E, int caloEnum, int clusterizerEnum, unsigned short primaryTrackSource, bool useProjection = true);

// ANCHOR main function to be called in event loop
void runclusterizer(
  int clusterizerEnum,
  int caloEnum,
  float seedE,
  float aggE,
  unsigned short primaryTrackSource,
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
  if(verbosityCLS>1)cout << "C3: new event" << endl;

  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!h_clusterizer_nonagg_towers[icalo][ialgo])h_clusterizer_nonagg_towers[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_nonagg_towers%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 50,0,50,120,0,0.6);

    }
  }

  // vector of towers used in the clusterizer
  std::vector<towersStrct> input_towers;
  std::vector<towersStrct> cluster_1Cell_towers;
  // vector of towers within the currently found cluster
  std::vector<towersStrct> cluster_towers;

  // fill vector with towers for clusterization above aggregation threshold
  float E_scaling = 2; // = 2 would be a first calibration fix
  if(caloEnum==kFHCAL){
    for(int itow=0; itow<_nTowers_FHCAL; itow++){
      if(_tower_FHCAL_E[itow]>(aggE/E_scaling)){
        towersStrct tempstructT;
        tempstructT.tower_E = _tower_FHCAL_E[itow]*E_scaling;
        tempstructT.tower_iEta = _tower_FHCAL_iEta[itow];
        tempstructT.tower_iPhi = _tower_FHCAL_iPhi[itow];
        tempstructT.tower_trueID = GetCorrectMCArrayEntry(_tower_FHCAL_trueID[itow]);
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
        tempstructT.tower_trueID = GetCorrectMCArrayEntry(_tower_FEMC_trueID[itow]);
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

      // ANCHOR C3 clusterizer logic
      if(clusterizerEnum==kC3){
//         cout << "running C3" << endl;
        for (int tit = 1; tit < input_towers.size(); tit++){
          // towers must be within cross of delta Eta and delta Phi <= 1
          if( ( std::abs(input_towers.at(tit).tower_iEta-input_towers.at(0).tower_iEta) + std::abs(input_towers.at(tit).tower_iPhi-input_towers.at(0).tower_iPhi) ) <= 1 ){
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
      // ANCHOR C5 clusterizer logic
      else if(clusterizerEnum==kC5){
//         cout << "running C5" << endl;
        for (int tit = 1; tit < input_towers.size(); tit++){
          // towers must be within cross of delta Eta and delta Phi <= 2 (-> 3x3 plus 4 additional towers)
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
      // ANCHOR 3x3 clusterizer logic
      else if(clusterizerEnum==k3x3){
//         cout << "running 3x3" << endl;
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
      // ANCHOR 5x5 clusterizer logic
      else if(clusterizerEnum==k5x5){
//         cout << "running 5x5" << endl;
        for (int tit = 1; tit < input_towers.size(); tit++){
          // towers must be within 5x5 matrix (delta Eta and delta Phi < 2)
          if(std::abs(input_towers.at(tit).tower_iEta-input_towers.at(0).tower_iEta)<3){
            if(std::abs(input_towers.at(tit).tower_iPhi-input_towers.at(0).tower_iPhi)<3){
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
//         cout << "running V3" << endl;
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
              if(input_towers.at(ait).tower_E >= (cluster_towers.at(tit).tower_E + aggregation_margin_V3)) continue;
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

      // ANCHOR MA clusterizer logic
      else if(clusterizerEnum==kMA){
//         cout << "running MA" << endl;
        // remove seed tower from sample
        input_towers.erase(input_towers.begin());
        for (int tit = 0; tit < cluster_towers.size(); tit++){
          // Now go recursively to all neighbours and add them to the cluster if they fulfill the conditions
          int iEtaTwr = cluster_towers.at(tit).tower_iEta;
          int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
          for (int ait = 0; ait < input_towers.size(); ait++){
            int iEtaTwrAgg = input_towers.at(ait).tower_iEta;
            int iPhiTwrAgg = input_towers.at(ait).tower_iPhi;
            // first condition asks for V3-like neighbors, while second condition also checks diagonally attached towers
            if( ((TMath::Abs(iEtaTwrAgg-iEtaTwr)+TMath::Abs(iPhiTwrAgg-iPhiTwr)) == 1) || ((TMath::Abs(iEtaTwrAgg-iEtaTwr)==1) && (TMath::Abs(iPhiTwrAgg-iPhiTwr)==1))){
              // only aggregate towers with lower energy than current tower
              if(input_towers.at(ait).tower_E >= (cluster_towers.at(tit).tower_E + aggregation_margin_V3)) continue;
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
      // ANCHOR V1 clusterizer logic
      else if(clusterizerEnum==kV1){
//         cout << "running MA" << endl;
        // remove seed tower from sample
        input_towers.erase(input_towers.begin());
        for (int tit = 0; tit < cluster_towers.size(); tit++){
          // Now go recursively to all neighbours and add them to the cluster if they fulfill the conditions
          int iEtaTwr = cluster_towers.at(tit).tower_iEta;
          int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
          for (int ait = 0; ait < input_towers.size(); ait++){
            int iEtaTwrAgg = input_towers.at(ait).tower_iEta;
            int iPhiTwrAgg = input_towers.at(ait).tower_iPhi;
            // first condition asks for V3-like neighbors, while second condition also checks diagonally attached towers
            if( ((TMath::Abs(iEtaTwrAgg-iEtaTwr)+TMath::Abs(iPhiTwrAgg-iPhiTwr)) == 1) || ((TMath::Abs(iEtaTwrAgg-iEtaTwr)==1) && (TMath::Abs(iPhiTwrAgg-iPhiTwr)==1))){
              // only aggregate towers with lower energy than current tower
              // if(input_towers.at(ait).tower_E >= (cluster_towers.at(tit).tower_E + aggregation_margin_V3)) continue;
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
      float* showershape_eta_phi = CalculateM02andWeightedPosition(cluster_towers, weightM02, clusters_E[nclusters], caloEnum, false);//(clusterizerEnum==kMA && caloEnum==kFHCAL) ? true : false);
      clusters_M02[nclusters] = showershape_eta_phi[0];
      clusters_M20[nclusters] = showershape_eta_phi[1];
      clusters_Eta[nclusters] = showershape_eta_phi[2];
      clusters_Phi[nclusters] = showershape_eta_phi[3];
      clusters_X[nclusters] = showershape_eta_phi[4];
      clusters_Y[nclusters] = showershape_eta_phi[5];
      clusters_Z[nclusters] = showershape_eta_phi[6];
      clusters_isMatched[nclusters] = isClusterMatched(nclusters, 20,clusters_X,clusters_Y,clusters_E, caloEnum, clusterizerEnum, primaryTrackSource, true);
      if(verbosityCLS>1) cout << clusterizerEnum << "\t" << nclusters << "\tC3 cluster with E = " << clusters_E[nclusters] << "\tEta: " << clusters_Eta[nclusters]<< "\tPhi: " << clusters_Phi[nclusters]<< "\tntowers: " << clusters_NTower[nclusters] << "\ttrueID: " << clusters_trueID[nclusters] << endl;
      // remove clusterized towers
      if(!(clusterizerEnum==kV3) && !(clusterizerEnum==kMA)){
        input_towers.erase(input_towers.begin());
      }

      // apply calibration if desired
      if(_doClusterECalibration){
          clusters_E[nclusters]/=getCalibrationValue(clusters_E[nclusters], caloEnum, clusterizerEnum);
          clusters_E[nclusters]*=getEnergySmearing( caloEnum, clusterizerEnum);
      }
      nclusters++;
    } else {
      for (int ait = 0; ait < input_towers.size(); ait++){
        h_clusterizer_nonagg_towers[caloEnum][clusterizerEnum]->Fill(input_towers.size(),input_towers.at(ait).tower_E);
      }
      input_towers.clear();
    }
  }
}


// ANCHOR track/projection matching function
// TODO at some point we might want E/p
bool isClusterMatched(int clsID, float matchingwindow, float* clusters_X, float* clusters_Y, float* clusters_E, int caloEnum, int clusterizerEnum, unsigned short primaryTrackSource, bool useProjection = true){
  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!h_clusterizer_all_dx_dy[icalo][ialgo])h_clusterizer_all_dx_dy[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_all_dx_dy_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 100,-100,100, 100,-100,100);
      if(!h_clusterizer_matched_dx_dy[icalo][ialgo])h_clusterizer_matched_dx_dy[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_matched_dx_dy_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 100,-100,100, 100,-100,100);
    }
  }
  if(useProjection){
    int projectionlayer = 0;
    if(caloEnum==kFHCAL)projectionlayer = 5;
    else if(caloEnum==kFEMC)projectionlayer = 6;
    else {
      cout << "isClusterMatched: caloEnum not defined, returning FALSE" << endl;
      return false;
    }
    for(Int_t iproj=0; iproj<_nProjections; iproj++){
      if(_track_ProjLayer[iproj]!=projectionlayer) continue;
      if(_track_Proj_t[iproj]<-9000) continue;
      if(_track_Proj_x[iproj]==0 && _track_Proj_y[iproj]==0) continue;
      if(clusters_X[clsID]==0 && clusters_Y[clsID]==0) continue;
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
      // Select track source
      if (_track_source[itrk] != primaryTrackSource) { continue; }
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

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_CLSIZER.root",outputDir.Data()),"RECREATE");

  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(h_clusterizer_nonagg_towers[icalo][ialgo]) h_clusterizer_nonagg_towers[icalo][ialgo]->Write();
      if(h_clusterizer_all_dx_dy[icalo][ialgo]) h_clusterizer_all_dx_dy[icalo][ialgo]->Write();
      if(h_clusterizer_matched_dx_dy[icalo][ialgo]) h_clusterizer_matched_dx_dy[icalo][ialgo]->Write();
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
