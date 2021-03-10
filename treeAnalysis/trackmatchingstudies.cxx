#include <algorithm>
// ANCHOR debug output verbosity
int verbosityTM = 1;


TH2F*  h_TMstudies_dx_dy[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]

// ANCHOR track/projection matching function
// TODO at some point we might want E/p
bool trackmatchingstudies(int enum_clusterizer, int enum_calorimeter, bool useProjection = false){

  int _nclusters_TM = 0;
  float* clusters_TM_E            = new float[_maxNclusters];
  float* clusters_TM_Eta         = new float[_maxNclusters];
  float* clusters_TM_Phi         = new float[_maxNclusters];
  float* clusters_TM_M02         = new float[_maxNclusters];
  float* clusters_TM_M20         = new float[_maxNclusters];
  bool* clusters_TM_isMatched         = new bool[_maxNclusters];
  int* clusters_TM_NTower         = new int[_maxNclusters];
  int* clusters_TM_trueID       = new int[_maxNclusters];
  int* clusters_TM_NtrueID       = new int[_maxNclusters];
  float* clusters_TM_X         = new float[_maxNclusters];
  float* clusters_TM_Y         = new float[_maxNclusters];
  float* clusters_TM_Z         = new float[_maxNclusters];
  loadClusterizerInput(
    enum_clusterizer,
    enum_calorimeter,
    _nclusters_TM,
    clusters_TM_E,
    clusters_TM_Eta,
    clusters_TM_Phi,
    clusters_TM_M02,
    clusters_TM_M20,
    clusters_TM_isMatched,
    clusters_TM_NTower,
    clusters_TM_trueID,
    clusters_TM_NtrueID,
    clusters_TM_X,
    clusters_TM_Y,
    clusters_TM_Z);

  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!h_TMstudies_dx_dy[icalo][ialgo])h_TMstudies_dx_dy[icalo][ialgo] 	= new TH2F(Form("h_TMstudies_dx_dy_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 100,-100,100, 100,-100,100);
    }
  }

  int projectionlayer = -1;
  if(enum_calorimeter==kFHCAL) projectionlayer = 2;
  if(enum_calorimeter==kFEMC) projectionlayer = 1;

  for(Int_t iclus=0; iclus<_nclusters_TM; iclus++){
    // if(!(clusters_TM_M02[iclus]<0.1))continue;
    if(!(clusters_TM_Eta[iclus]>3.0))continue;
    // if(!(clusters_TM_NTower[iclus]<3))continue;
    // bool matchfound = false;
    if(useProjection){
      for(Int_t iproj=0; iproj<_nProjections; iproj++){
        // if(matchfound) continue;
        if(_track_ProjLayer[iproj]!=projectionlayer) continue;
        // if((abs(_track_Proj_x[iclus])<15 && abs(_track_Proj_y[iclus])<15))continue;
        if(_track_Proj_t[iproj]<-9000) continue;
        if(_track_Proj_x[iproj]==0 && _track_Proj_y[iproj]==0) continue;
        if(clusters_TM_X[iproj]==0 && clusters_TM_Y[iproj]==0) continue;
        // check eta difference
        h_TMstudies_dx_dy[enum_calorimeter][enum_clusterizer]->Fill(_track_Proj_x[iproj]-clusters_TM_X[iclus],_track_Proj_y[iproj]-clusters_TM_Y[iclus]);
        // if(TMath::Sqrt(TMath::Power(_track_Proj_x[iproj]-clusters_TM_X[iclus],2)+TMath::Power(_track_Proj_y[iproj]-clusters_TM_Y[iclus],2)) <15)matchfound=true;
        // TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
        // float projeta = projvec.Eta();
        // float projphi = (projvec.Phi()<0 ? projvec.Phi()+TMath::Pi() : projvec.Phi()-TMath::Pi());
        // h_TMstudies_dx_dy[enum_calorimeter][enum_clusterizer]->Fill(projeta-clusters_TM_Eta[iclus],projphi-clusters_TM_Phi[iclus]);
      }
    } else {
      for(Int_t itrk=0; itrk<_nTracks; itrk++){
        // check eta difference
        h_TMstudies_dx_dy[enum_calorimeter][enum_clusterizer]->Fill(_track_px[itrk]-clusters_TM_X[iclus],_track_py[itrk]-clusters_TM_Y[iclus]);
        // TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
        // float trketa = trackvec.Eta();
        // float trkphi = (trackvec.Phi()<0 ? trackvec.Phi()+TMath::Pi() : trackvec.Phi()-TMath::Pi());
        // h_TMstudies_dx_dy[enum_calorimeter][enum_clusterizer]->Fill(trketa-clusters_TM_Eta[iclus],trkphi-clusters_TM_Phi[iclus]);
        // if(abs(trackvec.Eta()-clusters_TM_Eta[iclus]) < matchingwindow){
          // if(verbosityTM>1) cout << "\tTrack matched in eta! trk: " << trackvec.Mag() << "\tcls: " << clusters_TM_E[iclus] << "\tdelta: " << abs(trackvec.Eta()-clusters_TM_Eta[iclus]) << endl;
          // float trkphi = (trackvec.Phi()<0 ? trackvec.Phi()+TMath::Pi() : trackvec.Phi()-TMath::Pi());
          // // check phi difference
          // if(abs(trkphi-clusters_TM_Phi[iclus]) < matchingwindow){
          //   if(verbosityTM>1) cout << "\tTrack matched in phi! trk: " << trackvec.Mag() << "\tcls: " << clusters_TM_E[iclus] << "\tdelta: " << abs(trkphi-clusters_TM_Phi[iclus]) << endl;
          //   return true;
          // }
        // }
      }
    }
  }
  return false;
}

void trackmatchingstudiesSave(){
// make output directory
    gSystem->Exec("mkdir -p "+outputDir + "/trackmatchingstudies");

  // define output file
  TFile* fileOutput = new TFile(Form("%s/trackmatchingstudies/output_TMSTUD.root",outputDir.Data()),"RECREATE");

  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(h_TMstudies_dx_dy[icalo][ialgo]) h_TMstudies_dx_dy[icalo][ialgo]->Write();
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}