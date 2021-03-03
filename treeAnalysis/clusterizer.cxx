#include <algorithm>
// ANCHOR debug output verbosity
int verbosityCLS = 1;

// ANCHOR define global variables
enum clusterizertype {
    kV1         = 0,
    kV3         = 1,
    kNxN        = 2,
    kXN         = 3
};

bool _do_V1clusterizer = false;
bool _do_V3clusterizer = false;
bool _do_XNclusterizer = false;
bool _do_NxNclusterizer = false;

bool _do_V1clusterizerFEMC = false;
bool _do_V3clusterizerFEMC = false;
bool _do_XNclusterizerFEMC = false;
bool _do_NxNclusterizerFEMC = false;


bool isClusterMatched(int clsID, float matchingwindow, float* clusters_Eta, float* clusters_Phi, float* clusters_E, bool useProjection = false, int projectionlayer = 0);

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
  int* &clusters_NtrueID
){
  nclusters = 0;
  if(verbosityCLS>1)cout << "XN: new event" << endl;

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
        tempstructT.tower_trueID = _tower_FHCAL_trueID[itow]-1;
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
        tempstructT.tower_trueID = _tower_FEMC_trueID[itow]-1;
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
      clusters_isMatched[nclusters] = isClusterMatched(nclusters, 0.3,clusters_Eta,clusters_Phi,clusters_E);
     // if(verbosityCLS>1) cout << "\tXN cluster with E = " << clusters_E[nclusters] << "\tEta: " << clusters_Eta[nclusters]<< "\tPhi: " << clusters_Phi[nclusters]<< "\tntowers: " << clusters_NTower[nclusters] << "\ttrueID: " << clusters_trueID[nclusters] << endl;
      // remove clusterized towers
      if(!(clusterizerEnum==kV3)){
        input_towers.erase(input_towers.begin());
      }
      nclusters++;
    } else {
      input_towers.clear();
    }
  }
}


// ANCHOR track/projection matching function
// TODO at some point we might want E/p
bool isClusterMatched(int clsID, float matchingwindow, float* clusters_Eta, float* clusters_Phi, float* clusters_E, bool useProjection = false, int projectionlayer /*FTTL2*/){
  if(useProjection){
    for(Int_t iproj=0; iproj<_nProjections; iproj++){
      if(_track_ProjLayer[iproj]!=projectionlayer) continue;
      // check eta difference
      TVector3 trackvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
      if(abs(trackvec.Eta()-clusters_Eta[clsID]) < matchingwindow){
        if(verbosityCLS>1) cout << "\tProj matched in eta! trk: " << trackvec.Eta() << "\tcls: " << clusters_Eta[clsID] << "\tdelta: " << abs(trackvec.Eta()-clusters_Eta[clsID]) << endl;
        float trkphi = (trackvec.Phi()<0 ? trackvec.Phi()+TMath::Pi() : trackvec.Phi()-TMath::Pi());
        // check phi difference
        if(abs(trkphi-clusters_Phi[clsID]) < matchingwindow){
          if(verbosityCLS>1) cout << "\tProj matched in phi! trk: " << trackvec.Mag() << "\tcls: " << clusters_E[clsID] << "\tdelta: " << abs(trkphi-clusters_Phi[clsID]) << endl;
          return true;
        }
      }
    }
  } else {
    for(Int_t itrk=0; itrk<_nTracks; itrk++){
      // check eta difference
      TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
      if(abs(trackvec.Eta()-clusters_Eta[clsID]) < matchingwindow){
        if(verbosityCLS>1) cout << "\tTrack matched in eta! trk: " << trackvec.Mag() << "\tcls: " << clusters_E[clsID] << "\tdelta: " << abs(trackvec.Eta()-clusters_Eta[clsID]) << endl;
        float trkphi = (trackvec.Phi()<0 ? trackvec.Phi()+TMath::Pi() : trackvec.Phi()-TMath::Pi());
        // check phi difference
        if(abs(trkphi-clusters_Phi[clsID]) < matchingwindow){
          if(verbosityCLS>1) cout << "\tTrack matched in phi! trk: " << trackvec.Mag() << "\tcls: " << clusters_E[clsID] << "\tdelta: " << abs(trkphi-clusters_Phi[clsID]) << endl;
          return true;
        }
      }
    }
  }
  return false;
}
