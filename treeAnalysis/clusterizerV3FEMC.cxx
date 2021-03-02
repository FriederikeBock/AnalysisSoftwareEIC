#include <algorithm>
// ANCHOR debug output verbosity
int verbosityV3FEMC = 1;

bool _do_V3clusterizerFEMC = false;

// ANCHOR define global variables
int _nclusters_V3_FEMC = 0;
float* _clusters_V3_FEMC_E            = new float[_maxNclusters];
float* _clusters_V3_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_V3_FEMC_Phi         = new float[_maxNclusters];
float* _clusters_V3_FEMC_M02         = new float[_maxNclusters];
float* _clusters_V3_FEMC_M20         = new float[_maxNclusters];
bool* _clusters_V3_FEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_V3_FEMC_NTower         = new int[_maxNclusters];
int* _clusters_V3_FEMC_trueID       = new int[_maxNclusters];

bool isClusterMatchedV3FEMC(int clsID, float matchingwindow = 0.3, bool useProjection = false, int projectionlayer = 0);


// ANCHOR main function to be called in event loop
void runV3clusterizerFEMC(float seedE_V3 = 0.3, float aggE_V3 = 0.1){
  _nclusters_V3_FEMC = 0;
  if(verbosityV3FEMC>1)cout << "V3: new event" << endl;

  // vector of towers used in the clusterizer
  std::vector<towersStrct> input_towers;
  // vector of towers within the currently found cluster
  std::vector<towersStrct> cluster_towers;

  // fill vector with towers for clusterization above aggregation threshold
  for(int itow=0; itow<_nTowers_FEMC; itow++){
    if(_tower_FEMC_E[itow]>aggE_V3){
      towersStrct tempstructT;
      tempstructT.tower_E = _tower_FEMC_E[itow];
      tempstructT.tower_iEta = _tower_FEMC_iEta[itow];
      tempstructT.tower_iPhi = _tower_FEMC_iPhi[itow];
      tempstructT.tower_trueID = _tower_FEMC_trueID[itow];
      input_towers.push_back(tempstructT);
    }
  }
  // sort vector in descending energy order
  std::sort(input_towers.begin(), input_towers.end(), &acompare);

  while (!input_towers.empty()) {
    cluster_towers.clear();
    // always start with highest energetic tower
    if(input_towers.at(0).tower_E > seedE_V3){
      // fill seed cell information into current cluster
      _clusters_V3_FEMC_E[_nclusters_V3_FEMC] = input_towers.at(0).tower_E;
      _clusters_V3_FEMC_NTower[_nclusters_V3_FEMC] = 1;
      _clusters_V3_FEMC_trueID[_nclusters_V3_FEMC] = input_towers.at(0).tower_trueID; // TODO save all MC labels?
      cluster_towers.push_back(input_towers.at(0));

      // remove seed tower from sample
      input_towers.erase(input_towers.begin());

      // loop over all aggregated towers
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
              _clusters_V3_FEMC_E[_nclusters_V3_FEMC]+=input_towers.at(ait).tower_E;
              _clusters_V3_FEMC_NTower[_nclusters_V3_FEMC]++;
              cluster_towers.push_back(input_towers.at(ait));
              input_towers.erase(input_towers.begin()+ait);
            }
          }
      }

      // determine remaining cluster properties from its towers
      float* showershape_eta_phi = CalculateM02andWeightedPosition(cluster_towers, weightM02, _clusters_V3_FEMC_E[_nclusters_V3_FEMC], true, kFALSE);
      _clusters_V3_FEMC_M02[_nclusters_V3_FEMC] = showershape_eta_phi[0];
      _clusters_V3_FEMC_M20[_nclusters_V3_FEMC] = showershape_eta_phi[1];
      _clusters_V3_FEMC_Eta[_nclusters_V3_FEMC] = showershape_eta_phi[2];
      _clusters_V3_FEMC_Phi[_nclusters_V3_FEMC] = showershape_eta_phi[3];
      _clusters_V3_FEMC_isMatched[_nclusters_V3_FEMC] = isClusterMatchedV3FEMC(_nclusters_V3_FEMC);
      // if(verbosityV3FEMC>1) cout << "\tV3 cluster with E = " << _clusters_V3_FEMC_E[_nclusters_V3_FEMC] << "\tEta: " << _clusters_V3_FEMC_Eta[_nclusters_V3_FEMC]<< "\tPhi: " << _clusters_V3_FEMC_Phi[_nclusters_V3_FEMC]<< "\tntowers: " << _clusters_V3_FEMC_NTower[_nclusters_V3_FEMC] << "\ttrueID: " << _clusters_V3_FEMC_trueID[_nclusters_V3_FEMC] << endl;
      // remove clusterized towers

      _nclusters_V3_FEMC++;
    } else {
      input_towers.clear();
    }
  }
}


// ANCHOR track/projection matching function
// TODO at some point we might want E/p
bool isClusterMatchedV3FEMC(int clsID, float matchingwindow = 0.3, bool useProjection = false, int projectionlayer = 2 /*FTTL2*/){
  if(useProjection){
    for(Int_t iproj=0; iproj<_nProjections; iproj++){
      if(_track_ProjLayer[iproj]!=projectionlayer) continue;
      // check eta difference
      TVector3 trackvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
      if(abs(trackvec.Eta()-_clusters_V3_FEMC_Eta[clsID]) < matchingwindow){
        if(verbosityV3FEMC>1) cout << "\tProj matched in eta! trk: " << trackvec.Eta() << "\tcls: " << _clusters_V3_FEMC_Eta[clsID] << "\tdelta: " << abs(trackvec.Eta()-_clusters_V3_FEMC_Eta[clsID]) << endl;
        float trkphi = (trackvec.Phi()<0 ? trackvec.Phi()+TMath::Pi() : trackvec.Phi()-TMath::Pi());
        // check phi difference
        if(abs(trkphi-_clusters_V3_FEMC_Phi[clsID]) < matchingwindow){
          if(verbosityV3FEMC>1) cout << "\tProj matched in phi! trk: " << trackvec.Mag() << "\tcls: " << _clusters_V3_FEMC_E[clsID] << "\tdelta: " << abs(trkphi-_clusters_V3_FEMC_Phi[clsID]) << endl;
          return true;
        }
      }
    }
  } else {
    for(Int_t itrk=0; itrk<_nTracks; itrk++){
      // check eta difference
      TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
      if(abs(trackvec.Eta()-_clusters_V3_FEMC_Eta[clsID]) < matchingwindow){
        if(verbosityV3FEMC>1) cout << "\tTrack matched in eta! trk: " << trackvec.Mag() << "\tcls: " << _clusters_V3_FEMC_E[clsID] << "\tdelta: " << abs(trackvec.Eta()-_clusters_V3_FEMC_Eta[clsID]) << endl;
        float trkphi = (trackvec.Phi()<0 ? trackvec.Phi()+TMath::Pi() : trackvec.Phi()-TMath::Pi());
        // check phi difference
        if(abs(trkphi-_clusters_V3_FEMC_Phi[clsID]) < matchingwindow){
          if(verbosityV3FEMC>1) cout << "\tTrack matched in phi! trk: " << trackvec.Mag() << "\tcls: " << _clusters_V3_FEMC_E[clsID] << "\tdelta: " << abs(trkphi-_clusters_V3_FEMC_Phi[clsID]) << endl;
          return true;
        }
      }
    }
  }
  return false;
}


// ANCHOR save function after event loop
void saveHistosV3FEMCclusterizer(){
  // make output directory
  gSystem->Exec("mkdir -p treeProcessing/V3clusterizer");
  // define output file
  TFile* fileOutput = new TFile("treeProcessing/V3clusterizer/output_V3.root","RECREATE");

  // write histograms (e.g. cluster spectrum, M02 distribution, etc)
  // h_pion_FEMC_E->Write();

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
