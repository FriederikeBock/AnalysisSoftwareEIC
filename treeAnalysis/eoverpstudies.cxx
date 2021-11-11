#include <algorithm>
// ANCHOR debug output verbosity
int verbosityEOP = 1;
bool initializeHistsEOP = false;
//************************************************************************************************************
//************************************************************************************************************
// declaration of histogramas and global variables
//************************************************************************************************************
//************************************************************************************************************

TH1F*  h_EoverP_trackSpec_pions[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_EoverP_trackSpec_electrons[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_EoverP_pions_EoverP_E_trackbased[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_EoverP_electrons_EoverP_E_trackbased[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_EoverP_pions_EoverP_E_projectionbased[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_EoverP_electrons_EoverP_E_projectionbased[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_pions_EoverP_E_Eta_trackbased[_active_calo][2]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_electrons_EoverP_E_Eta_trackbased[_active_calo][2]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_pions_EoverP_E_Eta_projectionbased[_active_calo][2]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_electrons_EoverP_E_Eta_projectionbased[_active_calo][2]; // [calorimeter_enum][algorithm_enum]

TString str_EoverP_FillTrue[2] = {"recP","trueP"};
//************************************************************************************************************
//************************************************************************************************************
// ANCHOR track/projection matching function
// TODO at some point we might want E/p
//************************************************************************************************************
//************************************************************************************************************
bool eoverpstudies( int primaryTrackSource = 0, bool runSpecialCuts = false){

  
  // boolean for single event histo?
  bool b_TMstudies_2D_track_filled[_active_calo] = {false};
  
  //************************************************************************************************************
  //*********************************  create histograms  ******************************************************
  //************************************************************************************************************
  if (!initializeHistsEOP){
    for(int icalo=0;icalo<_active_calo;icalo++){
      if (!caloEnabled[icalo]) continue;
      int nbinsEStud = 500;
      int nbinsEoP = 200;
      if(!h_EoverP_trackSpec_pions[icalo])
        h_EoverP_trackSpec_pions[icalo]     = new TH1F(Form("h_EoverP_trackSpec_pions_%s",str_calorimeter[icalo].Data()), "", nbinsEStud,0,50);
      if(!h_EoverP_trackSpec_electrons[icalo])
        h_EoverP_trackSpec_electrons[icalo]     = new TH1F(Form("h_EoverP_trackSpec_electrons_%s",str_calorimeter[icalo].Data()), "", nbinsEStud,0,50);
      if(!h_EoverP_pions_EoverP_E_trackbased[icalo])
        h_EoverP_pions_EoverP_E_trackbased[icalo]        = new TH2F(Form("h_EoverP_pions_EoverP_E_trackbased_%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
      if(!h_EoverP_electrons_EoverP_E_trackbased[icalo])
        h_EoverP_electrons_EoverP_E_trackbased[icalo]        = new TH2F(Form("h_EoverP_electrons_EoverP_E_trackbased_%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
      if(!h_EoverP_pions_EoverP_E_projectionbased[icalo])
        h_EoverP_pions_EoverP_E_projectionbased[icalo]        = new TH2F(Form("h_EoverP_pions_EoverP_E_projectionbased_%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
      if(!h_EoverP_electrons_EoverP_E_projectionbased[icalo])
        h_EoverP_electrons_EoverP_E_projectionbased[icalo]        = new TH2F(Form("h_EoverP_electrons_EoverP_E_projectionbased_%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
        // versus Eta
      for(int ifilltrue=0;ifilltrue<2;ifilltrue++){
        if(!h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][ifilltrue])
          h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_pions_EoverP_E_Eta_trackbased_%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50, 160, -4, 4);
        if(!h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][ifilltrue])
          h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_electrons_EoverP_E_Eta_trackbased_%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50, 160, -4, 4);
        if(!h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][ifilltrue])
          h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_pions_EoverP_E_Eta_projectionbased_%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50, 160, -4, 4);
        if(!h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][ifilltrue])
          h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_electrons_EoverP_E_Eta_projectionbased_%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50, 160, -4, 4);
      }
    }
    initializeHistsEOP = true;
  }
  
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    Int_t caloDir           = GetCaloDirection(icalo);
    Float_t minEtaCurCalo   = partEtaCalo[minEtaBinCalo[caloDir]];
    Float_t maxEtaCurCalo   = partEtaCalo[maxEtaBinCalo[caloDir]];
    bool isFwd = IsForwardCalorimeter(icalo);
    bool b_TMstudies_2D_projection_filled = false;
    if(icalo==kFEMC){
      minEtaCurCalo = nominalEtaRegion[2][0];
      maxEtaCurCalo = nominalEtaRegion[2][1];
    } else if(icalo==kBECAL){
      minEtaCurCalo = nominalEtaRegion[1][0];
      maxEtaCurCalo = nominalEtaRegion[1][1];
    } else if(icalo==kEEMC){
      minEtaCurCalo = nominalEtaRegion[0][0];
      maxEtaCurCalo = nominalEtaRegion[0][1];
    }

    // loop over all different algorithms
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!loadClusterizerInput( ialgo, icalo)) continue;

      int projectionlayer = ReturnProjectionIndexForCalorimeter(icalo, false);
      float zHC = isFwd ? ReturnFwdCalorimeterPosition(icalo) : 1;
      float matchingwdw = ReturnTrackMatchingWindowForCalo(icalo);

      for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus++){
        bool ismatched_track      = false;
        bool ismatched_projection = false;
        int mcID                = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID;
        
        // if(((_clusters_calo[ialgo][icalo].at(iclus)).cluster_M02<0.1))continue;
          if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta < minEtaCurCalo || (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta > maxEtaCurCalo) continue;
        
        // if(!((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta>3.0))continue;
        // if(((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers<2))continue;
        // if(((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta>3.2))continue;

        // find largest cluster for current MC label and make matching plots for only the largest cluster

        for(Int_t iproj=0; iproj<_nProjections; iproj++){
          // perform matching only if projection layer is the correct one for the current calorimeter
          if(_track_ProjLayer[iproj]!=projectionlayer) continue;
          if (_track_source[(int)_track_ProjTrackID[iproj]] != primaryTrackSource) { continue; }
          // if timing is -9000 then projection failed
          if(_track_Proj_t[iproj]<-9000) continue;
          // if x and y in forward direction is 0, then projection can not be used
          if(isFwd && _track_Proj_true_x[iproj]==0 && _track_Proj_true_y[iproj]==0) continue;
          // create projection vector from x,y,z positions
          TVector3 projvec(_track_Proj_true_x[iproj],_track_Proj_true_y[iproj],_track_Proj_true_z[iproj]);
          float projeta = projvec.Eta();
          float projphi = projvec.Phi();
          if(!isFwd && projphi==0 && projeta==0) continue;
          // this condition is called only once to fill 2D histograms with all projections

          // if a cluster is found at x,y=0 (inside the beampipe), then something went wrong
          if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_X==0 && (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y==0) continue;
          // check deta/dphi or dx/dy difference
          float delta_1 = isFwd ? _track_Proj_true_x[iproj]-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_X : projeta-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta;
          float delta_2 = isFwd ? _track_Proj_true_y[iproj]-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y : projphi-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi;

          // if(TMath::Sqrt(TMath::Power(delta_1,2)+TMath::Power(_track_Proj_true_y[iproj]-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y,2)) <15)matchfound=true;
          // h_TMstudies_2D_delta_projection[0][icalo][ialgo]->Fill(projeta-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,projphi-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi);
          if(abs(delta_1) < matchingwdw){
            if(abs(delta_2) < matchingwdw){
              ismatched_projection = true;
            }
          }
          if(ialgo==kMA){
            TLorentzVector PiGen(_mcpart_px[mcID], _mcpart_py[mcID], _mcpart_pz[mcID], _mcpart_E[mcID]);
            // // TODO condition 1
            // if((Int_t)_clusters_calo[ialgo][icalo].size()==1){ // only one cluster per event
              if(ReturnCaloFwdBarrelBckw(icalo)==0 && _track_pz[(int)_track_ProjTrackID[iproj]]<0) continue;
              if(ReturnCaloFwdBarrelBckw(icalo)==2 && _track_pz[(int)_track_ProjTrackID[iproj]]>0) continue;
              if (_track_source[(int)_track_ProjTrackID[iproj]] != primaryTrackSource) { continue; }
              float trackmom = sqrt(pow(_track_px[(int)_track_ProjTrackID[iproj]], 2) + pow(_track_py[(int)_track_ProjTrackID[iproj]], 2) + pow(_track_pz[(int)_track_ProjTrackID[iproj]], 2));
              // // TODO condition 2
              // if( ((_clusters_calo[ialgo][icalo].at(0)).cluster_E) / trackmom > 1.05 ) { continue; }
              if(abs(_mcpart_PDG[(int)_track_trueID[(int)_track_ProjTrackID[iproj]]])==211){
                if(ismatched_projection){
                  h_EoverP_pions_EoverP_E_projectionbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,PiGen.P());
                  h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][0]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,PiGen.P(),PiGen.Eta());
                  h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][1]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/PiGen.P(),PiGen.P(),PiGen.Eta());
                  // h_EoverP_pions_EoverP_E_projectionbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,trackmom);
                  // h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,trackmom,projeta);
                }
              }
              if(abs(_mcpart_PDG[(int)_track_trueID[(int)_track_ProjTrackID[iproj]]])==11){
                if(ismatched_projection){
                  h_EoverP_electrons_EoverP_E_projectionbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,PiGen.P());
                  h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][0]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,PiGen.P(),PiGen.Eta());
                  h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][1]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/PiGen.P(),PiGen.P(),PiGen.Eta());
                  // h_EoverP_electrons_EoverP_E_projectionbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,trackmom);
                  // h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,trackmom,projeta);
                }
              }
            // }
          }
        }
        b_TMstudies_2D_projection_filled = true;

        for(Int_t itrk=0; itrk<_nTracks; itrk++){
          // don't match backward tracks to forward calo and vise versa
          if(ReturnCaloFwdBarrelBckw(icalo)==0 && _track_pz[itrk]<0) continue;
          if(ReturnCaloFwdBarrelBckw(icalo)==2 && _track_pz[itrk]>0) continue;
          if (_track_source[itrk] != primaryTrackSource) { continue; }

          // check eta difference
          TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
          if(isFwd)trackvec*=abs(zHC/trackvec.Z());
          float tracketa = trackvec.Eta();
          float trackphi = trackvec.Phi(); //(trackvec.Phi()<0 ? trackvec.Phi()+TMath::Pi() : trackvec.Phi()-TMath::Pi());

          if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_X==0 && (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y==0) continue;
          float delta_1 = isFwd ? trackvec.X()-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_X : tracketa-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta;
          float delta_2 = isFwd ? trackvec.Y()-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y : trackphi-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi;
          if(abs(delta_1) < matchingwdw){
            if(abs(delta_2) < matchingwdw){
              ismatched_track = true;
            }
          }
          if(ialgo==kMA){
            TLorentzVector PiGen(_mcpart_px[mcID], _mcpart_py[mcID], _mcpart_pz[mcID], _mcpart_E[mcID]);
            // // TODO condition 1
            // if((Int_t)_clusters_calo[ialgo][icalo].size()==1){ // only one cluster per event
              float trackmom = sqrt(pow(_track_px[itrk], 2) + pow(_track_py[itrk], 2) + pow(_track_pz[itrk], 2));
              // // TODO condition 2
              // if( ((_clusters_calo[ialgo][icalo].at(0)).cluster_E) / trackmom > 1.05 ) { continue; }
              if(abs(_mcpart_PDG[(int)_track_trueID[itrk]])==211){
                if(ismatched_track){
                  h_EoverP_pions_EoverP_E_trackbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,PiGen.P());
                  h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][0]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,PiGen.P(),PiGen.Eta());
                  h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][1]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/PiGen.P(),PiGen.P(),PiGen.Eta());
                  // h_EoverP_pions_EoverP_E_trackbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,trackmom);
                  // h_EoverP_pions_EoverP_E_Eta_trackbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,trackmom,tracketa);
                }
                h_EoverP_trackSpec_pions[icalo]->Fill(trackmom);
              }
              if(abs(_mcpart_PDG[(int)_track_trueID[itrk]])==11){
                if(ismatched_track){
                  h_EoverP_electrons_EoverP_E_trackbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,PiGen.P());
                  h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][0]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,PiGen.P(),PiGen.Eta());
                  h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][1]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/PiGen.P(),PiGen.P(),PiGen.Eta());
                  // h_EoverP_electrons_EoverP_E_trackbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,trackmom);
                  // h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackmom,trackmom,tracketa);
                }
                h_EoverP_trackSpec_electrons[icalo]->Fill(trackmom);
              }
            // }
          }
        }
        b_TMstudies_2D_track_filled[icalo] = true;
      }
    }
  }

  return false;
}

//************************************************************************************************************
//************************************************************************************************************
// Function to save the histograms for the track matching
//************************************************************************************************************
//************************************************************************************************************
void eoverpstudiesSave(){

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_EOPSTUD.root",outputDir.Data()),"RECREATE");
  // fileOutput->mkdir("etabins");
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;

    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    if(h_EoverP_pions_EoverP_E_trackbased[icalo])h_EoverP_pions_EoverP_E_trackbased[icalo]->Write();
    if(h_EoverP_electrons_EoverP_E_trackbased[icalo])h_EoverP_electrons_EoverP_E_trackbased[icalo]->Write();
    if(h_EoverP_pions_EoverP_E_projectionbased[icalo])h_EoverP_pions_EoverP_E_projectionbased[icalo]->Write();
    if(h_EoverP_electrons_EoverP_E_projectionbased[icalo])h_EoverP_electrons_EoverP_E_projectionbased[icalo]->Write();
    for(int ifilltrue=0;ifilltrue<2;ifilltrue++){
      if(h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][ifilltrue])h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][ifilltrue]->Write();
      if(h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][ifilltrue])h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][ifilltrue]->Write();
      if(h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][ifilltrue])h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][ifilltrue]->Write();
      if(h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][ifilltrue])h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][ifilltrue]->Write();
    }
    if(h_EoverP_trackSpec_pions[icalo])h_EoverP_trackSpec_pions[icalo]->Write();
    if(h_EoverP_trackSpec_electrons[icalo])h_EoverP_trackSpec_electrons[icalo]->Write();
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
