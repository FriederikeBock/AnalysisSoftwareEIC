#include <algorithm>
// ANCHOR debug output verbosity
int verbosityTM = 1;


TH2F*  h_TMstudies_x_y_clus[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_x_y_track; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_x_y_projection[_active_calo]; // [calorimeter_enum][algorithm_enum]

TH2F*  h_TMstudies_dx_dy_track[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_dy_projection[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_isMatched[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_isMatched_E[_active_calo]; // [calorimeter_enum][algorithm_enum]

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

// ANCHOR track/projection matching function
// TODO at some point we might want E/p
bool trackmatchingstudies(){
  if(!h_TMstudies_x_y_track)h_TMstudies_x_y_track 	= new TH2F(Form("h_TMstudies_x_y_track"), "", 250,-250,250, 250,-250,250);
  bool b_TMstudies_x_y_track_filled = false;
  for(int icalo=0;icalo<_active_calo;icalo++){
      if(!h_TMstudies_x_y_projection[icalo])h_TMstudies_x_y_projection[icalo] 	= new TH2F(Form("h_TMstudies_x_y_projection_%s",str_calorimeter[icalo].Data()), "", 250,-250,250, 250,-250,250);
  }
  for(int icalo=0;icalo<_active_calo;icalo++){
      bool b_TMstudies_x_y_projection_filled = false;

      if(!h_TMstudies_isMatched[icalo]){
        h_TMstudies_isMatched[icalo] 	= new TH1F(Form("h_TMstudies_isMatched_%s",str_calorimeter[icalo].Data()), "", _active_algo*3,0,_active_algo*3);
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          h_TMstudies_isMatched[icalo]->GetXaxis()->SetBinLabel(ialgo*3+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
          h_TMstudies_isMatched[icalo]->GetXaxis()->SetBinLabel(ialgo*3+2, "track matched");
          h_TMstudies_isMatched[icalo]->GetXaxis()->SetBinLabel(ialgo*3+3, "proj matched");
        }
      }
      if(!h_TMstudies_isMatched_E[icalo]){
        h_TMstudies_isMatched_E[icalo] 	= new TH2F(Form("h_TMstudies_isMatched_E_%s",str_calorimeter[icalo].Data()), "", _active_algo*3,0,_active_algo*3,100,0,100);
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          h_TMstudies_isMatched_E[icalo]->GetXaxis()->SetBinLabel(ialgo*3+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
          h_TMstudies_isMatched_E[icalo]->GetXaxis()->SetBinLabel(ialgo*3+2, "track matched");
          h_TMstudies_isMatched_E[icalo]->GetXaxis()->SetBinLabel(ialgo*3+3, "proj matched");
        }
      }
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!h_TMstudies_dx_dy_track[icalo][ialgo])h_TMstudies_dx_dy_track[icalo][ialgo] 	= new TH2F(Form("h_TMstudies_dx_dy_track_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 100,-100,100, 100,-100,100);
      if(!h_TMstudies_dx_dy_projection[icalo][ialgo])h_TMstudies_dx_dy_projection[icalo][ialgo] 	= new TH2F(Form("h_TMstudies_dx_dy_projection_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 100,-100,100, 100,-100,100);
      if(!h_TMstudies_x_y_clus[icalo][ialgo])h_TMstudies_x_y_clus[icalo][ialgo] 	= new TH2F(Form("h_TMstudies_x_y_clus_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 250,-250,250, 250,-250,250);

      if(!loadClusterizerInput(
        ialgo,
        icalo,
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
        clusters_TM_Z))continue;

      if(_nclusters_TM==0) continue;
      int projectionlayer = -1;
      if(icalo==kFHCAL) projectionlayer = 5;
      if(icalo==kFEMC) projectionlayer = 6;
      float zHC = 400;
      if(icalo==kFHCAL)zHC=350;
      if(icalo==kFEMC)zHC=291;

      float matchingwdw = 20;

      for(Int_t iclus=0; iclus<_nclusters_TM; iclus++){
        bool ismatched_track = false;
        bool ismatched_projection = false;
        // if(!(clusters_TM_M02[iclus]<0.1))continue;
        // if(!(clusters_TM_Eta[iclus]>3.0))continue;
        // if(!(clusters_TM_NTower[iclus]<3))continue;
        // bool matchfound = false;
        h_TMstudies_isMatched[icalo]->Fill(3*ialgo);
        h_TMstudies_isMatched_E[icalo]->Fill(3*ialgo,clusters_TM_E[iclus] );
        h_TMstudies_x_y_clus[icalo][ialgo]->Fill(clusters_TM_X[iclus],clusters_TM_Y[iclus]);
        for(Int_t iproj=0; iproj<_nProjections; iproj++){
          // if(matchfound) continue;
          // if((abs(_track_Proj_x[iclus])<15 && abs(_track_Proj_y[iclus])<15))continue;
          if(_track_Proj_t[iproj]<-9000) continue;
          if(_track_Proj_x[iproj]==0 && _track_Proj_y[iproj]==0) continue;
          if(!b_TMstudies_x_y_projection_filled){
            if(_track_ProjLayer[iproj]==kFHCAL)h_TMstudies_x_y_projection[kFHCAL]->Fill(_track_Proj_x[iproj],_track_Proj_y[iproj]);
            if(_track_ProjLayer[iproj]==kFEMC)h_TMstudies_x_y_projection[kFEMC]->Fill(_track_Proj_x[iproj],_track_Proj_y[iproj]);
          }
          if(_track_ProjLayer[iproj]!=projectionlayer) continue;
          if(clusters_TM_X[iclus]==0 && clusters_TM_Y[iclus]==0) continue;
          // check eta difference
          h_TMstudies_dx_dy_projection[icalo][ialgo]->Fill(_track_Proj_x[iproj]-clusters_TM_X[iclus],_track_Proj_y[iproj]-clusters_TM_Y[iclus]);
          // if(TMath::Sqrt(TMath::Power(_track_Proj_x[iproj]-clusters_TM_X[iclus],2)+TMath::Power(_track_Proj_y[iproj]-clusters_TM_Y[iclus],2)) <15)matchfound=true;
          // TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
          // float projeta = projvec.Eta();
          // float projphi = (projvec.Phi()<0 ? projvec.Phi()+TMath::Pi() : projvec.Phi()-TMath::Pi());
          // h_TMstudies_dx_dy_projection[icalo][ialgo]->Fill(projeta-clusters_TM_Eta[iclus],projphi-clusters_TM_Phi[iclus]);
          if(abs(_track_Proj_x[iproj]-clusters_TM_X[iclus]) < matchingwdw){
            if(abs(_track_Proj_y[iproj]-clusters_TM_Y[iclus]) < matchingwdw){
              ismatched_projection = true;
            }
          }
        }
        b_TMstudies_x_y_projection_filled = true;

        for(Int_t itrk=0; itrk<_nTracks; itrk++){
          // check eta difference
          TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
          trackvec*=(zHC/trackvec.Z());
          if(!b_TMstudies_x_y_track_filled)h_TMstudies_x_y_track->Fill(trackvec.X(),trackvec.Y());

          if(clusters_TM_X[iclus]==0 && clusters_TM_Y[iclus]==0) continue;

          h_TMstudies_dx_dy_track[icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[iclus],trackvec.Y()-clusters_TM_Y[iclus]);
          if(abs(trackvec.X()-clusters_TM_X[iclus]) < matchingwdw){
            if(abs(trackvec.Y()-clusters_TM_Y[iclus]) < matchingwdw){
              ismatched_track = true;
            }
          }
        }
        b_TMstudies_x_y_track_filled = true;

        if(ismatched_track){
          h_TMstudies_isMatched[icalo]->Fill(3*ialgo+1);
          h_TMstudies_isMatched_E[icalo]->Fill(3*ialgo+1,clusters_TM_E[iclus] );
        }
        if(ismatched_projection){
          h_TMstudies_isMatched[icalo]->Fill(3*ialgo+2);
          h_TMstudies_isMatched_E[icalo]->Fill(3*ialgo+2,clusters_TM_E[iclus] );
        }
      }
    }
  }
  return false;
}

void trackmatchingstudiesSave(){

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_TMSTUD.root",outputDir.Data()),"RECREATE");
  if(h_TMstudies_x_y_track)h_TMstudies_x_y_track->Write();
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(h_TMstudies_isMatched[icalo])h_TMstudies_isMatched[icalo]->Write();
    if(h_TMstudies_isMatched_E[icalo])h_TMstudies_isMatched_E[icalo]->Write();
    if(h_TMstudies_x_y_projection[icalo])h_TMstudies_x_y_projection[icalo]->Write();
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(h_TMstudies_x_y_clus[icalo][ialgo]) h_TMstudies_x_y_clus[icalo][ialgo]->Write();
      if(h_TMstudies_dx_dy_track[icalo][ialgo]) h_TMstudies_dx_dy_track[icalo][ialgo]->Write();
      if(h_TMstudies_dx_dy_projection[icalo][ialgo]) h_TMstudies_dx_dy_projection[icalo][ialgo]->Write();
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
