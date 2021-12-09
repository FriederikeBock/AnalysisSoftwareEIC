#include <algorithm>
// ANCHOR debug output verbosity
int verbosityEOP = 1;
bool initializeHistsEOP = false;
bool b_EoverP_dotracks = false;
//************************************************************************************************************
//************************************************************************************************************
// declaration of histogramas and global variables
//************************************************************************************************************
//************************************************************************************************************

TH1F*  h_EoverP_trackSpec_pions[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_EoverP_trackSpec_electrons[_active_calo]; // [calorimeter_enum][algorithm_enum]
// track based hists
TH2F*  h_EoverP_pions_EoverP_E_trackbased[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_EoverP_electrons_EoverP_E_trackbased[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_pions_EoverP_E_Eta_trackbased[_active_calo][2]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_electrons_EoverP_E_Eta_trackbased[_active_calo][2]; // [calorimeter_enum][algorithm_enum]
// projection based hists
TH2F*  h_EoverP_pions_EoverP_E_projectionbased[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_EoverP_electrons_EoverP_E_projectionbased[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_pions_EoverP_E_Eta_projectionbased[_active_calo][2]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_electrons_EoverP_E_Eta_projectionbased[_active_calo][2]; // [calorimeter_enum][algorithm_enum]
// projection based with M02
TH2F*  h_EoverP_pions_EoverP_E_projectionbasedWithM02[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_EoverP_electrons_EoverP_E_projectionbasedWithM02[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_pions_EoverP_E_Eta_projectionbasedWithM02[_active_calo][2]; // [calorimeter_enum][algorithm_enum]
TH3F*  h_EoverP_electrons_EoverP_E_Eta_projectionbasedWithM02[_active_calo][2]; // [calorimeter_enum][algorithm_enum]

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
      if (IsHCALCalorimeter(icalo)) continue;
      int nbinsEStud = 500;
      int nbinsEoP = 200;
      if(!h_EoverP_trackSpec_pions[icalo])
        h_EoverP_trackSpec_pions[icalo]     = new TH1F(Form("h_EoverP_trackSpec_pions_%s",str_calorimeter[icalo].Data()), "", nbinsEStud,0,50);
      if(!h_EoverP_trackSpec_electrons[icalo])
        h_EoverP_trackSpec_electrons[icalo]     = new TH1F(Form("h_EoverP_trackSpec_electrons_%s",str_calorimeter[icalo].Data()), "", nbinsEStud,0,50);
      if(b_EoverP_dotracks){
        if(!h_EoverP_pions_EoverP_E_trackbased[icalo])
          h_EoverP_pions_EoverP_E_trackbased[icalo]        = new TH2F(Form("h_EoverP_pions_EoverP_E_trackbased_%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
        if(!h_EoverP_electrons_EoverP_E_trackbased[icalo])
          h_EoverP_electrons_EoverP_E_trackbased[icalo]        = new TH2F(Form("h_EoverP_electrons_EoverP_E_trackbased_%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
      }
      if(!h_EoverP_pions_EoverP_E_projectionbased[icalo])
        h_EoverP_pions_EoverP_E_projectionbased[icalo]        = new TH2F(Form("h_EoverP_pions_EoverP_E_projectionbased_%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
      if(!h_EoverP_electrons_EoverP_E_projectionbased[icalo])
        h_EoverP_electrons_EoverP_E_projectionbased[icalo]        = new TH2F(Form("h_EoverP_electrons_EoverP_E_projectionbased_%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
      if(!h_EoverP_pions_EoverP_E_projectionbasedWithM02[icalo])
        h_EoverP_pions_EoverP_E_projectionbasedWithM02[icalo]        = new TH2F(Form("h_EoverP_pions_EoverP_E_projectionbasedWithM02%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
      if(!h_EoverP_electrons_EoverP_E_projectionbasedWithM02[icalo])
        h_EoverP_electrons_EoverP_E_projectionbasedWithM02[icalo]        = new TH2F(Form("h_EoverP_electrons_EoverP_E_projectionbasedWithM02%s",str_calorimeter[icalo].Data()), "", nbinsEoP, 0, 4, nbinsEStud, 0, 50);
        // versus Eta
      for(int ifilltrue=0;ifilltrue<2;ifilltrue++){
        if(b_EoverP_dotracks){
          if(!h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][ifilltrue])
            h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_pions_EoverP_E_Eta_trackbased_%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 2, nbinsEStud, 0, 50, 160, -4, 4);
          if(!h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][ifilltrue])
            h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_electrons_EoverP_E_Eta_trackbased_%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 2, nbinsEStud, 0, 50, 160, -4, 4);
        }
        if(!h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][ifilltrue])
          h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_pions_EoverP_E_Eta_projectionbased_%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 2, nbinsEStud, 0, 50, 160, -4, 4);
        if(!h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][ifilltrue])
          h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_electrons_EoverP_E_Eta_projectionbased_%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 2, nbinsEStud, 0, 50, 160, -4, 4);
        if(!h_EoverP_pions_EoverP_E_Eta_projectionbasedWithM02[icalo][ifilltrue])
          h_EoverP_pions_EoverP_E_Eta_projectionbasedWithM02[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_pions_EoverP_E_Eta_projectionbasedWithM02%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 2, nbinsEStud, 0, 50, 160, -4, 4);
        if(!h_EoverP_electrons_EoverP_E_Eta_projectionbasedWithM02[icalo][ifilltrue])
          h_EoverP_electrons_EoverP_E_Eta_projectionbasedWithM02[icalo][ifilltrue]        = new TH3F(Form("h_EoverP_electrons_EoverP_E_Eta_projectionbasedWithM02%s_%s",str_calorimeter[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data()), "", nbinsEoP, 0, 2, nbinsEStud, 0, 50, 160, -4, 4);
      }
    }
    initializeHistsEOP = true;
  }
  
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    if (IsHCALCalorimeter(icalo)) continue;
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
        bool isHighestForPart = true;

        // if(((_clusters_calo[ialgo][icalo].at(iclus)).cluster_M02<0.1))continue;
        if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta < minEtaCurCalo || (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta > maxEtaCurCalo) continue;
        if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers < 2) continue;

        for(Int_t iclus2=0; iclus2<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus2++){
          if( iclus == iclus2) continue;
          if ( (_clusters_calo[ialgo][icalo].at(iclus2)).cluster_NTowers<2) continue;
          if ( mcID == (_clusters_calo[ialgo][icalo].at(iclus2)).cluster_trueID  ) {
            if ( (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E < (_clusters_calo[ialgo][icalo].at(iclus2)).cluster_E  ){
              isHighestForPart = false;
            }
          }
        }
        if(!isHighestForPart) continue;
        Double_t trackP = 0;
        Double_t currR = 10000;
        if (ialgo==kMA && (_clusters_calo[ialgo][icalo].at(iclus)).cluster_isMatched){
          int trackIDMatched = -1;
          for (Int_t k = 0; k < (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).size(); k++){
            double matching_R = 10000;
            if(isFwd) matching_R = sqrt(pow((((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(k)).dX,2)+pow((((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(k)).dY,2));
            else      matching_R = sqrt(pow((((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(k)).dEta,2)+pow((((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(k)).dPhi,2));
            if(matching_R < currR){
              currR = matching_R;
              trackIDMatched = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(k)).id;
            }
          }
          if(currR==10000) cout << "something is wrong" << endl;
          TVector3 recpartvec(_track_px[trackIDMatched],_track_py[trackIDMatched],_track_pz[trackIDMatched]);
          trackP = recpartvec.Mag();

          TLorentzVector PiGen(_mcpart_px[mcID], _mcpart_py[mcID], _mcpart_pz[mcID], _mcpart_E[mcID]);
          if (verbosityEOP > 4)
            std::cout << _mcpart_E[mcID] << "\t" << icalo << "\t" << GetM02CutValueForEMC(icalo, _mcpart_E[mcID]) << std::endl;
          
          if(abs(_mcpart_PDG[(int)_track_trueID[trackIDMatched]])==211){
            h_EoverP_pions_EoverP_E_projectionbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackP,PiGen.P());
            h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][0]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackP,PiGen.P(),PiGen.Eta());
            h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][1]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/PiGen.P(),PiGen.P(),PiGen.Eta());
            if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_M02 < GetM02CutValueForEMC(icalo, _mcpart_E[mcID])){ 
              h_EoverP_pions_EoverP_E_projectionbasedWithM02[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackP,PiGen.P());
              h_EoverP_pions_EoverP_E_Eta_projectionbasedWithM02[icalo][0]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackP,PiGen.P(),PiGen.Eta());
              h_EoverP_pions_EoverP_E_Eta_projectionbasedWithM02[icalo][1]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/PiGen.P(),PiGen.P(),PiGen.Eta());
            }
          }
          if(abs(_mcpart_PDG[(int)_track_trueID[trackIDMatched]])==11){
            h_EoverP_electrons_EoverP_E_projectionbased[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackP,PiGen.P());
            h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][0]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackP,PiGen.P(),PiGen.Eta());
            h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][1]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/PiGen.P(),PiGen.P(),PiGen.Eta());
            if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_M02 < GetM02CutValueForEMC(icalo, _mcpart_E[mcID])){ 
              h_EoverP_electrons_EoverP_E_projectionbasedWithM02[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackP,PiGen.P());
              h_EoverP_electrons_EoverP_E_Eta_projectionbasedWithM02[icalo][0]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/trackP,PiGen.P(),PiGen.Eta());
              h_EoverP_electrons_EoverP_E_Eta_projectionbasedWithM02[icalo][1]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/PiGen.P(),PiGen.P(),PiGen.Eta());
            }
          }
        }

        b_TMstudies_2D_projection_filled = true;
        if(b_EoverP_dotracks){
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
    if (IsHCALCalorimeter(icalo)) continue;
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    if(b_EoverP_dotracks){
      if(h_EoverP_pions_EoverP_E_trackbased[icalo])h_EoverP_pions_EoverP_E_trackbased[icalo]->Write();
      if(h_EoverP_electrons_EoverP_E_trackbased[icalo])h_EoverP_electrons_EoverP_E_trackbased[icalo]->Write();
    }
    if(h_EoverP_pions_EoverP_E_projectionbased[icalo])h_EoverP_pions_EoverP_E_projectionbased[icalo]->Write();
    if(h_EoverP_electrons_EoverP_E_projectionbased[icalo])h_EoverP_electrons_EoverP_E_projectionbased[icalo]->Write();
    if(h_EoverP_pions_EoverP_E_projectionbasedWithM02[icalo])h_EoverP_pions_EoverP_E_projectionbasedWithM02[icalo]->Write();
    if(h_EoverP_electrons_EoverP_E_projectionbasedWithM02[icalo])h_EoverP_electrons_EoverP_E_projectionbasedWithM02[icalo]->Write();
    for(int ifilltrue=0;ifilltrue<2;ifilltrue++){
      if(b_EoverP_dotracks){
        if(h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][ifilltrue])h_EoverP_pions_EoverP_E_Eta_trackbased[icalo][ifilltrue]->Write();
        if(h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][ifilltrue])h_EoverP_electrons_EoverP_E_Eta_trackbased[icalo][ifilltrue]->Write();
      }
      if(h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][ifilltrue])h_EoverP_pions_EoverP_E_Eta_projectionbased[icalo][ifilltrue]->Write();
      if(h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][ifilltrue])h_EoverP_electrons_EoverP_E_Eta_projectionbased[icalo][ifilltrue]->Write();
      if(h_EoverP_pions_EoverP_E_Eta_projectionbasedWithM02[icalo][ifilltrue])h_EoverP_pions_EoverP_E_Eta_projectionbasedWithM02[icalo][ifilltrue]->Write();
      if(h_EoverP_electrons_EoverP_E_Eta_projectionbasedWithM02[icalo][ifilltrue])h_EoverP_electrons_EoverP_E_Eta_projectionbasedWithM02[icalo][ifilltrue]->Write();
    }
    if(h_EoverP_trackSpec_pions[icalo])h_EoverP_trackSpec_pions[icalo]->Write();
    if(h_EoverP_trackSpec_electrons[icalo])h_EoverP_trackSpec_electrons[icalo]->Write();
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
