#include <algorithm>
// ANCHOR debug output verbosity
int verbosityCLS = 1;

// ANCHOR define global variables
float aggregation_margin_V3 = 0.03;
float aggregation_margin_MA = 0.03;

TH2F*  h_clusterizer_nonagg_towers[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_matched_2D_delta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_all_2D_delta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
bool _doClusterECalibration = true;


// ANCHOR track/projection matching function
// TODO at some point we might want E/p
bool isClusterMatched(  int clsID, 
                        float* clusters_X, 
                        float* clusters_Y, 
                        float* clusters_Eta, 
                        float* clusters_Phi, 
                        float* clusters_E, 
                        int caloEnum, 
                        int clusterizerEnum, 
                        unsigned short primaryTrackSource, 
                        bool useProjection = true
                     ){
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    bool isFwd = IsForwardCalorimeter(icalo);
    float nbins2Ddelta = isFwd ? 78 : 80;
    float min2Ddeltahist = isFwd ? -39 : -0.2;
    float max2Ddeltahist = isFwd ? 39 : 0.2;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!h_clusterizer_all_2D_delta[icalo][ialgo])h_clusterizer_all_2D_delta[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_all_2D_delta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
      if(!h_clusterizer_matched_2D_delta[icalo][ialgo])h_clusterizer_matched_2D_delta[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_matched_2D_delta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
    }
  }
  // in case the cluster position is 0,0 -> no matching to be done
  if(clusters_X[clsID]==0 && clusters_Y[clsID]==0) return false;

  float matchingwindow = ReturnTrackMatchingWindowForCalo(caloEnum);
  bool isFwd = IsForwardCalorimeter(caloEnum);

  if(useProjection){
    int projectionlayer = ReturnProjectionIndexForCalorimeter(caloEnum);
    if(projectionlayer==-1) return false;
    for(Int_t iproj=0; iproj<_nProjections; iproj++){
      // check for correct projection layer
      if(_track_ProjLayer[iproj]!=projectionlayer) continue;
      // bad timing = bad projection
      if(_track_Proj_t[iproj]<-9000) continue;
      // no projection should end up on the beampipe
      if(isFwd && _track_Proj_x[iproj]==0 && _track_Proj_y[iproj]==0) continue;

      // for the barrel region we need eta/phi instead of x/y
      float projeta = -20;
      float projphi = -20;
      if(!isFwd){
        TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
        projeta = projvec.Eta();
        projphi = projvec.Phi();
        if((projphi==0 && projeta==0) || (projphi==-20 && projeta==-20)) continue;
      }

      // calculate delta x/y or delta eta/phi between projection and cluster
      float delta_1 = isFwd ? _track_Proj_x[iproj]-clusters_X[clsID] : projeta-clusters_Eta[clsID];
      float delta_2 = isFwd ? _track_Proj_y[iproj]-clusters_Y[clsID] : projphi-clusters_Phi[clsID];

      // fill delta histogram
      h_clusterizer_all_2D_delta[caloEnum][clusterizerEnum]->Fill(delta_1,delta_2);

      // check x or eta difference
      if(abs(delta_1) < matchingwindow){
        // check y or phi difference
        if(abs(delta_2) < matchingwindow){
          h_clusterizer_matched_2D_delta[caloEnum][clusterizerEnum]->Fill(delta_1,delta_2);
          return true; // cluster is matched to projection!
        }
      }
    }
  } else {
    float zHC = ReturnFwdCalorimeterPosition(caloEnum);
    for(Int_t itrk=0; itrk<_nTracks; itrk++){
      // Select track source
      if (_track_source[itrk] != primaryTrackSource) { continue; }
      // check eta difference
      TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
      float tracketa,trackphi = -20;
      if(isFwd){
        trackvec*=(zHC/trackvec.Z());
      } else {
          tracketa = trackvec.Eta();
          trackphi = trackvec.Phi();
          if((tracketa==0 && trackphi==0) || (tracketa==-20 && trackphi==-20)) continue;
      }

      // calculate delta x/y or delta eta/phi between projection and cluster
      float delta_1 = isFwd ? trackvec.X()-clusters_X[clsID] : tracketa-clusters_Eta[clsID];
      float delta_2 = isFwd ? trackvec.Y()-clusters_Y[clsID] : trackphi-clusters_Phi[clsID];

      // fill delta histogram
      h_clusterizer_all_2D_delta[caloEnum][clusterizerEnum]->Fill(delta_1,delta_2);

      // check x or eta difference
      if(abs(delta_1) < matchingwindow){
        // check y or phi difference
        if(abs(delta_2) < matchingwindow){
          h_clusterizer_matched_2D_delta[caloEnum][clusterizerEnum]->Fill(delta_1,delta_2);
          return true; // cluster is matched to track!
        }
      }
    }
  }
  return false;
}


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
  if(verbosityCLS>1)cout << "runclusterizer: new event" << endl;

  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
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
        tempstructT.tower_E       = _tower_FHCAL_E[itow]*E_scaling;
        tempstructT.tower_iEta    = _tower_FHCAL_iEta[itow];
        tempstructT.tower_iPhi    = _tower_FHCAL_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_FHCAL_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kFEMC){
    for(int itow=0; itow<_nTowers_FEMC; itow++){
      if(_tower_FEMC_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = 1.25*_tower_FEMC_E[itow]; // temp fix for ROS setting
        tempstructT.tower_iEta    = _tower_FEMC_iEta[itow];
        tempstructT.tower_iPhi    = _tower_FEMC_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_FEMC_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kCEMC){
    for(int itow=0; itow<_nTowers_CEMC; itow++){
      if(_tower_CEMC_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_CEMC_E[itow];
        tempstructT.tower_iEta    = _tower_CEMC_iEta[itow];
        tempstructT.tower_iPhi    = _tower_CEMC_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_CEMC_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kHCALIN){
    for(int itow=0; itow<_nTowers_HCALIN; itow++){
      if(_tower_HCALIN_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_HCALIN_E[itow];
        tempstructT.tower_iEta    = _tower_HCALIN_iEta[itow];
        tempstructT.tower_iPhi    = _tower_HCALIN_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_HCALIN_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kHCALOUT){
    for(int itow=0; itow<_nTowers_HCALOUT; itow++){
      if(_tower_HCALOUT_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_HCALOUT_E[itow];
        tempstructT.tower_iEta    = _tower_HCALOUT_iEta[itow];
        tempstructT.tower_iPhi    = _tower_HCALOUT_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_HCALOUT_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kLFHCAL){
    for(int itow=0; itow<_nTowers_LFHCAL; itow++){
      if(_tower_LFHCAL_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = 0.55*_tower_LFHCAL_E[itow]; // dirty fix
        tempstructT.tower_iEta    = _tower_LFHCAL_iEta[itow];
        tempstructT.tower_iPhi    = _tower_LFHCAL_iPhi[itow];
        tempstructT.tower_iL      = _tower_LFHCAL_iL[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_LFHCAL_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kDRCALO){
    for(int itow=0; itow<_nTowers_DRCALO; itow++){
      if(_tower_DRCALO_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_DRCALO_E[itow];
        tempstructT.tower_iEta    =   _tower_DRCALO_iEta[itow];
        tempstructT.tower_iPhi    = _tower_DRCALO_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_DRCALO_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kEHCAL){
    for(int itow=0; itow<_nTowers_EHCAL; itow++){
      if(_tower_EHCAL_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_EHCAL_E[itow];
        tempstructT.tower_iEta    = _tower_EHCAL_iEta[itow];
        tempstructT.tower_iPhi    = _tower_EHCAL_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_EHCAL_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kBECAL){
    for(int itow=0; itow<_nTowers_BECAL; itow++){
      if(_tower_BECAL_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_BECAL_E[itow];
        tempstructT.tower_iEta    = _tower_BECAL_iEta[itow];
        tempstructT.tower_iPhi    = _tower_BECAL_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_BECAL_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kEEMC){
    for(int itow=0; itow<_nTowers_EEMC; itow++){
      if(_tower_EEMC_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_EEMC_E[itow];
        tempstructT.tower_iEta    = _tower_EEMC_iEta[itow];
        tempstructT.tower_iPhi    = _tower_EEMC_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_EEMC_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kEEMCG){
    for(int itow=0; itow<_nTowers_EEMCG; itow++){
      if(_tower_EEMCG_E[itow]>aggE){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_EEMCG_E[itow];
        tempstructT.tower_iEta    = _tower_EEMCG_iEta[itow];
        tempstructT.tower_iPhi    = _tower_EEMCG_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_EEMCG_trueID[itow]);
        input_towers.push_back(tempstructT);
      }
    }
  } else {
    cout << "Incorrect calorimeter selected! Enum " << caloEnum << " not defined!" << endl;
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
      if (verbosityCLS > 2) std::cout << "new cluster" << std::endl;
      // fill seed cell information into current cluster
      clusters_E[nclusters] = input_towers.at(0).tower_E;
      clusters_NTower[nclusters] = 1;
      clusters_NtrueID[nclusters] = 1;
      clusters_trueID[nclusters] = input_towers.at(0).tower_trueID; // TODO save all MC labels?
      cluster_towers.push_back(input_towers.at(0));
      clslabels.push_back(input_towers.at(0).tower_trueID);
      if (verbosityCLS > 2) std::cout << "seed: "<<  input_towers.at(0).tower_iEta << "\t" << input_towers.at(0).tower_iPhi << "\t" << input_towers.at(0).tower_iL << "\t E:"<< clusters_E[nclusters] << std::endl;
      
      // ANCHOR C3 clusterizer logic
      if(clusterizerEnum==kC3){
//         cout << "running C3" << endl;
        for (int tit = 1; tit < (int)input_towers.size(); tit++){
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
            tit--;
          }
        }
      }
      // ANCHOR C5 clusterizer logic
      else if(clusterizerEnum==kC5){
//         cout << "running C5" << endl;
        for (int tit = 1; tit < (int)input_towers.size(); tit++){
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
            tit--;
          }
        }
      }
      // ANCHOR 3x3 clusterizer logic
      else if(clusterizerEnum==k3x3){
//         cout << "running 3x3" << endl;
        for (int tit = 1; tit < (int)input_towers.size(); tit++){
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
              tit--;
            }
          }
        }
      }
      // ANCHOR 5x5 clusterizer logic
      else if(clusterizerEnum==k5x5){
//         cout << "running 5x5" << endl;
        for (int tit = 1; tit < (int)input_towers.size(); tit++){
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
              tit--;
            }
          }
        }
      }

      // ANCHOR V3 clusterizer logic
      else if(clusterizerEnum==kV3){
//         cout << "running V3" << endl;
        // remove seed tower from sample
        input_towers.erase(input_towers.begin());
        for (int tit = 0; tit < (int)cluster_towers.size(); tit++){
          // Now go recursively to the next 4 neighbours and add them to the cluster if they fulfill the conditions
          int iEtaTwr = cluster_towers.at(tit).tower_iEta;
          int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
          for (int ait = 0; ait < (int)input_towers.size(); ait++){
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
              ait--;
            }
          }
        }
      }

      // ANCHOR MA clusterizer logic
      else if(clusterizerEnum==kMA){
//         cout << "running MA" << endl;
        // remove seed tower from sample
        input_towers.erase(input_towers.begin());
        for (int tit = 0; tit < (int)cluster_towers.size(); tit++){
          // Now go recursively to all neighbours and add them to the cluster if they fulfill the conditions
          int iEtaTwr = cluster_towers.at(tit).tower_iEta;
          int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
          int iLTwr   = cluster_towers.at(tit).tower_iL;
          int refC = 0;
          for (int ait = 0; ait < (int)input_towers.size(); ait++){
            int iEtaTwrAgg = input_towers.at(ait).tower_iEta;
            int iPhiTwrAgg = input_towers.at(ait).tower_iPhi;
            int iLTwrAgg   = input_towers.at(ait).tower_iL;
            int deltaL    = TMath::Abs(iLTwrAgg-iLTwr) ;
            int deltaPhi  = TMath::Abs(iPhiTwrAgg-iPhiTwr) ;
            int deltaEta  = TMath::Abs(iEtaTwrAgg-iEtaTwr) ;
            bool neighbor = (deltaL+deltaPhi+deltaEta == 1);
            bool corner2D = (deltaL == 0 && deltaPhi == 1 && deltaEta == 1) || (deltaL == 1 && deltaPhi == 0 && deltaEta == 1) || (deltaL == 1 && deltaPhi == 1 && deltaEta == 0);          
            // first condition asks for V3-like neighbors, while second condition also checks diagonally attached towers
            if(neighbor || corner2D ){

              // only aggregate towers with lower energy than current tower
              if(caloEnum != kLFHCAL){
                if(input_towers.at(ait).tower_E >= (cluster_towers.at(tit).tower_E + aggregation_margin_MA)) continue;
              } 
              clusters_E[nclusters]+=input_towers.at(ait).tower_E;
              clusters_NTower[nclusters]++;
              cluster_towers.push_back(input_towers.at(ait));
              if(!(std::find(clslabels.begin(), clslabels.end(), input_towers.at(ait).tower_trueID) != clslabels.end())){
                clusters_NtrueID[nclusters]++;
                clslabels.push_back(input_towers.at(ait).tower_trueID);
              }
              if (verbosityCLS > 2) std::cout << "aggregated: "<< iEtaTwrAgg << "\t" << iPhiTwrAgg << "\t" << iLTwrAgg << "\t E:" << input_towers.at(ait).tower_E << "\t reference: "<< refC << "\t"<< iEtaTwr << "\t" << iPhiTwr << "\t" << iLTwr << "\t cond.: \t"<< neighbor << "\t" << corner2D << "\t  diffs: " << deltaEta << "\t" << deltaPhi << "\t" << deltaL<< std::endl;

              input_towers.erase(input_towers.begin()+ait);
              ait--;
              refC++;
            }
          }
        }
      }
      // ANCHOR V1 clusterizer logic
      else if(clusterizerEnum==kV1){
//         cout << "running MA" << endl;
        // remove seed tower from sample
        input_towers.erase(input_towers.begin());
        for (int tit = 0; tit < (int)cluster_towers.size(); tit++){
          // Now go recursively to all neighbours and add them to the cluster if they fulfill the conditions
          int iEtaTwr = cluster_towers.at(tit).tower_iEta;
          int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
          for (int ait = 0; ait < (int)input_towers.size(); ait++){
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
              ait--;
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
      if (tracksEnabled){
        clusters_isMatched[nclusters] = isClusterMatched(nclusters,clusters_X,clusters_Y,clusters_Eta,clusters_Phi,clusters_E, caloEnum, clusterizerEnum, primaryTrackSource, true);
      } else {
        clusters_isMatched[nclusters] = false;
      }
        if(verbosityCLS>1) cout << clusterizerEnum << "\t" << nclusters << "\tcluster with E = " << clusters_E[nclusters] << "\tEta: " << clusters_Eta[nclusters]<< "\tPhi: " << clusters_Phi[nclusters]<< "\tX: " << clusters_X[nclusters]<< "\tX: " << clusters_X[nclusters]<< "\tZ: " << clusters_Z[nclusters]<< "\tntowers: " << clusters_NTower[nclusters] << "\ttrueID: " << clusters_trueID[nclusters] << endl;
      // remove clusterized towers
      if(!(clusterizerEnum==kV3) && !(clusterizerEnum==kMA) && !(clusterizerEnum==kV1)){
        input_towers.erase(input_towers.begin());
      }

      // apply calibration if desired
      if(_doClusterECalibration){
          clusters_E[nclusters]/=getCalibrationValue(clusters_E[nclusters], caloEnum, clusterizerEnum);
          clusters_E[nclusters]*=getEnergySmearing( caloEnum, clusterizerEnum);
      }
      nclusters++;
    } else {
      for (int ait = 0; ait < (int)input_towers.size(); ait++){
        h_clusterizer_nonagg_towers[caloEnum][clusterizerEnum]->Fill(input_towers.size(),input_towers.at(ait).tower_E);
      }
      input_towers.clear();
    }
  }
}


void clusterizerSave(){

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_CLSIZER.root",outputDir.Data()),"RECREATE");

  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(h_clusterizer_nonagg_towers[icalo][ialgo]) h_clusterizer_nonagg_towers[icalo][ialgo]->Write();
      if(h_clusterizer_all_2D_delta[icalo][ialgo]) h_clusterizer_all_2D_delta[icalo][ialgo]->Write();
      if(h_clusterizer_matched_2D_delta[icalo][ialgo]) h_clusterizer_matched_2D_delta[icalo][ialgo]->Write();
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
