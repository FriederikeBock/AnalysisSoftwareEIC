#include <algorithm>
// ANCHOR debug output verbosity
int verbosityTM = 1;
bool initializeHistsTMS = false;
//************************************************************************************************************
//************************************************************************************************************
// declaration of histogramas and global variables
//************************************************************************************************************
//************************************************************************************************************
TH2F*  h_TMstudies_2D_clus[_active_calo][_active_algo];                         // [calorimeter_enum][algorithm_enum]

TH1F*  h_TMstudies_clusSpec[_active_calo][_active_algo];                        // [calorimeter_enum][algorithm_enum]-
TH1F*  h_TMstudies_clusSpec_Eta[_active_calo][_active_algo][nEta+1];            // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_Matched[_active_calo][_active_algo];                // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_Matched_Eta[_active_calo][_active_algo][nEta+1];    // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_MatchedCalo[_active_calo][_active_algo];            // [calorimeter_enum][algorithm_enum]

TH2F*  h_TMstudies_2D_delta[_active_calo][_active_algo];                        // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsEta[_active_calo][_active_algo];                        // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dy_vsEta[_active_calo][_active_algo];                        // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsClsE[_active_calo][_active_algo][nEta+1];               // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dy_vsClsE[_active_calo][_active_algo][nEta+1];               // [calorimeter_enum][algorithm_enum]

TH2F*  h_TMstudiesCalo_2D_delta[_active_calo][_active_calo][_active_algo];            // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudiesCalo_dEta_vsEta[_active_calo][_active_calo][_active_algo];          // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudiesCalo_dPhi_vsEta[_active_calo][_active_calo][_active_algo];          // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudiesCalo_dEta_vsClsE[_active_calo][_active_calo][_active_algo][nEta+1]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudiesCalo_dPhi_vsClsE[_active_calo][_active_calo][_active_algo][nEta+1]; // [calorimeter_enum][algorithm_enum]

//************************************************************************************************************
//************************************************************************************************************
// ANCHOR track/projection matching function
//************************************************************************************************************
//************************************************************************************************************
bool trackmatchingstudies( int primaryTrackSource = 0){

  //************************************************************************************************************
  //*********************************  create histograms  ******************************************************
  //************************************************************************************************************
  if (!initializeHistsTMS){
    for(int icalo=0;icalo<_active_calo;icalo++){
      if (!caloEnabled[icalo]) continue;
      
      Int_t caloDir           = GetCaloDirection(icalo);
      Float_t minEtaCurCalo   = partEtaCalo[minEtaBinCalo[caloDir]];
      Float_t maxEtaCurCalo   = partEtaCalo[maxEtaBinCalo[caloDir]];
        
      bool isFwd              = IsForwardCalorimeter(icalo);
      float nbins2D           = isFwd ? 299 : 300;
      float min2Dhist         = isFwd ? -299 : -2.0;
      float max2Dhist         = isFwd ? 299 : 2.0;
      float min2Dhist2        = isFwd ? -299 : -TMath::Pi();
      float max2Dhist2        = isFwd ? 299 : TMath::Pi();
      TString str2DhistoName  = isFwd ? "x_y" : "eta_phi";
      float nbins2Ddelta      = isFwd ? 78*4 : 200;
      float min2Ddeltahist    = isFwd ? -39 : -0.25;
      float max2Ddeltahist    = isFwd ? 39 : 0.25;
      
      for(int ialgo=0;ialgo<_active_algo;ialgo++){
  
        h_TMstudies_clusSpec[icalo][ialgo]             = new TH1F(Form("h_TMstudies_clusSpec_%s_%s",str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 120,0,60);
        h_TMstudies_clusSpec_Matched[icalo][ialgo]     = new TH1F(Form("h_TMstudies_clusSpec_Matched_%s_%s",str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 120,0,60);

        h_TMstudies_2D_delta[icalo][ialgo]            = new TH2F(Form("h_TMstudies_2D_delta_%s_%s",str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 
                                                                                  nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
        h_TMstudies_dx_vsEta[icalo][ialgo]            = new TH2F(Form("h_TMstudies_dx_vsEta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 
                                                                                nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);
        h_TMstudies_dy_vsEta[icalo][ialgo]            = new TH2F(Form("h_TMstudies_dy_vsEta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 
                                                                                nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);

        h_TMstudies_2D_clus[icalo][ialgo]       = new TH2F(Form("h_TMstudies_%s_clus_%s_%s",str2DhistoName.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nbins2D, min2Dhist, max2Dhist, nbins2D, min2Dhist2, max2Dhist2);
        
        // add histos for HCal cluster to ECal cluster matching
        if (IsHCALCalorimeter(icalo)){
          h_TMstudies_clusSpec_MatchedCalo[icalo][ialgo] = new TH1F(Form("h_TMstudies_clusSpec_MatchedCalo_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), 
                                                                    "", 120,0,60);
          
          if (_combCalo[icalo] != -1){
            if (caloEnabled[_combCalo[icalo]]){
              h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]][ialgo]    = new TH2F(Form("h_TMstudiesCalo_2D_delta_%s_%s_%s",str_calorimeter[icalo].Data(),
                                                                                     str_calorimeter[_combCalo[icalo]].Data(), str_clusterizer[ialgo].Data()), "", 
                                                                                     nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
              h_TMstudiesCalo_dEta_vsEta[icalo][_combCalo[icalo]][ialgo]  = new TH2F(Form("h_TMstudiesCalo_dEta_vsEta_%s_%s_%s", str_calorimeter[icalo].Data(),
                                                                                    str_calorimeter[_combCalo[icalo]].Data(), str_clusterizer[ialgo].Data()), "", 
                                                                                    200, -0.25, 0.25, 
                                                                                    (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);
              h_TMstudiesCalo_dPhi_vsEta[icalo][_combCalo[icalo]][ialgo]  = new TH2F(Form("h_TMstudiesCalo_dPhi_vsEta_%s_%s_%s", str_calorimeter[icalo].Data(),
                                                                                    str_calorimeter[_combCalo[icalo]].Data(), str_clusterizer[ialgo].Data()), "", 
                                                                                    200, -0.25, 0.25, 
                                                                                    (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);
            }
          }
          
          if (_combCalo2[icalo] != -1){
            if (caloEnabled[_combCalo2[icalo]]){
              h_TMstudiesCalo_2D_delta[icalo][_combCalo2[icalo]][ialgo]   = new TH2F(Form("h_TMstudiesCalo_2D_delta_%s_%s_%s",str_calorimeter[icalo].Data(),
                                                                                     str_calorimeter[_combCalo2[icalo]].Data(), str_clusterizer[ialgo].Data()), "", 
                                                                                     200, -0.25, 0.25, 200, -0.25, 0.25);
              h_TMstudiesCalo_dEta_vsEta[icalo][_combCalo2[icalo]][ialgo] = new TH2F(Form("h_TMstudiesCalo_dEta_vsEta_%s_%s_%s", str_calorimeter[icalo].Data(),
                                                                                    str_calorimeter[_combCalo2[icalo]].Data(), str_clusterizer[ialgo].Data()), "", 
                                                                                    200, -0.25, 0.25, 
                                                                                    (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);
              h_TMstudiesCalo_dPhi_vsEta[icalo][_combCalo2[icalo]][ialgo] = new TH2F(Form("h_TMstudiesCalo_dPhi_vsEta_%s_%s_%s", str_calorimeter[icalo].Data(),
                                                                                    str_calorimeter[_combCalo2[icalo]].Data(), str_clusterizer[ialgo].Data()), "", 
                                                                                    200, -0.25, 0.25, 
                                                                                    (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);
            }
          }          
        }
        
        for (Int_t et = minEtaBinCalo[caloDir]; et<maxEtaBinCalo[caloDir]; et++){
          Double_t etaMin = partEtaCalo[0];
          Double_t etaMax = partEtaCalo[nEta];
          if (et < nEta){
            etaMin = partEtaCalo[et];
            etaMax = partEtaCalo[et+1];
          }
          // reconstructed  particles
          h_TMstudies_clusSpec_Eta[icalo][ialgo][et]          = new TH1F(Form("h_TMstudies_clusSpec_Eta_%1.1f_%1.1f_%s_%s", etaMin, etaMax, str_calorimeter[icalo].Data(),
                                                                               str_clusterizer[ialgo].Data()), "", 120, 0, 60);
          // reconstructed & matched particles
          h_TMstudies_clusSpec_Matched_Eta[icalo][ialgo][et]  = new TH1F(Form("h_TMstudies_clusSpec_Matched_Eta_%1.1f_%1.1f_%s_%s", etaMin, etaMax, str_calorimeter[icalo].Data(),
                                                                                str_clusterizer[ialgo].Data()), "", 120, 0, 60);
          h_TMstudies_dx_vsClsE[icalo][ialgo][et]           = new TH2F(Form("h_TMstudies_dx_vsClsE_Eta_%1.1f_%1.1f_%s_%s",etaMin,etaMax, str_calorimeter[icalo].Data(),
                                                                  str_clusterizer[ialgo].Data()), "", nBinsP, binningP,
                                                                  nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
          h_TMstudies_dy_vsClsE[icalo][ialgo][et]           = new TH2F(Form("h_TMstudies_dy_vsClsE_Eta_%1.1f_%1.1f_%s_%s",etaMin,etaMax, str_calorimeter[icalo].Data(),
                                                                  str_clusterizer[ialgo].Data()), "", nBinsP, binningP,
                                                                  nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
          if (IsHCALCalorimeter(icalo)){
            if (_combCalo[icalo] != -1){
              if (caloEnabled[_combCalo[icalo]]){
                h_TMstudiesCalo_dEta_vsClsE[icalo][_combCalo[icalo]][ialgo][et] = new TH2F(Form("h_TMstudiesCalo_dEta_vsClsE_Eta_%1.1f_%1.1f_%s_%s_%s", etaMin, etaMax,
                                                                                          str_calorimeter[icalo].Data(), str_calorimeter[_combCalo[icalo]].Data(),
                                                                                          str_clusterizer[ialgo].Data()), "", nBinsP, binningP,
                                                                                          200, -0.25, 0.25);
                h_TMstudiesCalo_dPhi_vsClsE[icalo][_combCalo[icalo]][ialgo][et] = new TH2F(Form("h_TMstudiesCalo_dPhi_vsClsE_Eta_%1.1f_%1.1f_%s_%s_%s", etaMin, etaMax,
                                                                                          str_calorimeter[icalo].Data(), str_calorimeter[_combCalo[icalo]].Data(),
                                                                                          str_clusterizer[ialgo].Data()), "", nBinsP, binningP,
                                                                                          200, -0.25, 0.25);
              }
            }
            if (_combCalo2[icalo] != -1){
              if (caloEnabled[_combCalo2[icalo]]){
                h_TMstudiesCalo_dEta_vsClsE[icalo][_combCalo2[icalo]][ialgo][et]= new TH2F(Form("h_TMstudiesCalo_dEta_vsClsE_Eta_%1.1f_%1.1f_%s_%s_%s", etaMin, etaMax,
                                                                                          str_calorimeter[icalo].Data(), str_calorimeter[_combCalo2[icalo]].Data(),
                                                                                          str_clusterizer[ialgo].Data()), "", nBinsP, binningP,
                                                                                          200, -0.25, 0.25);
                h_TMstudiesCalo_dPhi_vsClsE[icalo][_combCalo2[icalo]][ialgo][et]= new TH2F(Form("h_TMstudiesCalo_dPhi_vsClsE_Eta_%1.1f_%1.1f_%s_%s_%s", etaMin, etaMax,
                                                                                          str_calorimeter[icalo].Data(), str_calorimeter[_combCalo2[icalo]].Data(),
                                                                                          str_clusterizer[ialgo].Data()), "", nBinsP, binningP,
                                                                                          200, -0.25, 0.25);
              }
            }
          }
        }
      }
    }
    initializeHistsTMS = true;
  }
  
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    Int_t caloDir           = GetCaloDirection(icalo);
    Float_t minEtaCurCalo   = partEtaCalo[minEtaBinCalo[caloDir]];
    Float_t maxEtaCurCalo   = partEtaCalo[maxEtaBinCalo[caloDir]];
    bool isFwd = IsForwardCalorimeter(icalo);
  
    if(verbosityTM >2)std::cout << str_calorimeter[icalo].Data() << endl;
    // loop over all different algorithms
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!loadClusterizerInput( ialgo, icalo)) continue;

      for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus++){
        int mcID                = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID;
        if(verbosityTM >2) std::cout << "\t " << iclus << "\t" << (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E << std::endl;
        
        h_TMstudies_clusSpec[icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
        
        // determine eta bin
        Int_t et = 0;
        while (partEtaCalo[et+1] < _mcpart_Eta[mcID] && et < nEta) et++;

        if(et >= minEtaBinCalo[caloDir] && et<maxEtaBinCalo[caloDir]){
          h_TMstudies_clusSpec_Eta[icalo][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
        }
        h_TMstudies_2D_clus[icalo][ialgo]->Fill(isFwd ? (_clusters_calo[ialgo][icalo].at(iclus)).cluster_X : (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,
                                                isFwd ? (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y : (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi);
        
        if(verbosityTM >2) std::cout << "\t\t tracks: " << (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks.size()) << std::endl;
        for (int itr = 0; itr < (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks.size()); itr++){
          int trackID = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(itr)).id;
          float dPhi  = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(itr)).dPhi;
          float dEta  = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(itr)).dEta;
          float dX    = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(itr)).dX;
          float dY    = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(itr)).dY;
          
          // check deta/dphi or dx/dy difference
          float delta_1 = isFwd ? dX : dEta;
          float delta_2 = isFwd ? dY : dPhi;

          h_TMstudies_2D_delta[icalo][ialgo]->Fill(delta_1,delta_2);
          h_TMstudies_dx_vsEta[icalo][ialgo]->Fill(delta_1, _mcpart_Eta[mcID]);
          h_TMstudies_dy_vsEta[icalo][ialgo]->Fill(delta_2, _mcpart_Eta[mcID]);
          
          if(h_TMstudies_dx_vsClsE[icalo][ialgo][et]) h_TMstudies_dx_vsClsE[icalo][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, delta_1);
          if(h_TMstudies_dy_vsClsE[icalo][ialgo][et]) h_TMstudies_dy_vsClsE[icalo][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, delta_2);
        }
        if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks.size() > 0){
          h_TMstudies_clusSpec_Matched[icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
          if(et >= minEtaBinCalo[caloDir] && et<maxEtaBinCalo[caloDir]){
            h_TMstudies_clusSpec_Matched_Eta[icalo][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
          }
        }
        
        if( IsHCALCalorimeter(icalo)){
          if(verbosityTM >2) std::cout << "\t\t ecals: " << (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals.size()) << std::endl;
          for (int iCl2 = 0; iCl2 < (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals.size()); iCl2++){
            int clID = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).id;
            int calID = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).caloid;
            float dPhi  = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).dPhi;
            float dEta  = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).dEta;
            if(verbosityTM >2) std::cout << "\t\t\t"<< clID << "\t" << calID << "\t" << dPhi << "\t" << dEta << std::endl;
//             float dX    = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).dX;
//             float dY    = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).dY;
            if(h_TMstudiesCalo_2D_delta[icalo][calID][ialgo]) h_TMstudiesCalo_2D_delta[icalo][calID][ialgo]->Fill(dEta, dPhi);
            if(h_TMstudiesCalo_dPhi_vsEta[icalo][calID][ialgo]) h_TMstudiesCalo_dPhi_vsEta[icalo][calID][ialgo]->Fill(dPhi, _mcpart_Eta[mcID]);
            if(h_TMstudiesCalo_dEta_vsEta[icalo][calID][ialgo]) h_TMstudiesCalo_dEta_vsEta[icalo][calID][ialgo]->Fill(dEta, _mcpart_Eta[mcID]);
            if(h_TMstudiesCalo_dPhi_vsClsE[icalo][calID][ialgo][et]) h_TMstudiesCalo_dPhi_vsClsE[icalo][calID][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, dPhi);
            if(h_TMstudiesCalo_dEta_vsClsE[icalo][calID][ialgo][et]) h_TMstudiesCalo_dEta_vsClsE[icalo][calID][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, dEta);
          }
          if(verbosityTM >2) std::cout << "\t\t hcals: "  << "\t"<< (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals.size()) << std::endl;
          for (int iCl2 = 0; iCl2 < (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals.size()); iCl2++){
            int clID = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl2)).id;
            int calID = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl2)).caloid;
            float dPhi  = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl2)).dPhi;
            float dEta  = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl2)).dEta;
            if(verbosityTM >2) std::cout << "\t\t\t"<< clID << "\t" << calID << "\t" << dPhi << "\t" << dEta << std::endl;
//             float dX    = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl2)).dX;
//             float dY    = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl2)).dY;
            if(h_TMstudiesCalo_2D_delta[icalo][calID][ialgo]) h_TMstudiesCalo_2D_delta[icalo][calID][ialgo]->Fill(dEta, dPhi);
            if(h_TMstudiesCalo_dPhi_vsEta[icalo][calID][ialgo]) h_TMstudiesCalo_dPhi_vsEta[icalo][calID][ialgo]->Fill(dPhi, _mcpart_Eta[mcID]);
            if(h_TMstudiesCalo_dEta_vsEta[icalo][calID][ialgo]) h_TMstudiesCalo_dEta_vsEta[icalo][calID][ialgo]->Fill(dEta, _mcpart_Eta[mcID]);
            if(h_TMstudiesCalo_dPhi_vsClsE[icalo][calID][ialgo][et]) h_TMstudiesCalo_dPhi_vsClsE[icalo][calID][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, dPhi);
            if(h_TMstudiesCalo_dEta_vsClsE[icalo][calID][ialgo][et]) h_TMstudiesCalo_dEta_vsClsE[icalo][calID][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, dEta);
          }
          if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals.size() > 0){
            h_TMstudies_clusSpec_MatchedCalo[icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
          }
        }
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
void trackmatchingstudiesSave(){

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_TMSTUD.root",outputDir.Data()),"RECREATE");
  // fileOutput->mkdir("etabins");
  for(int icalo=0;icalo<_active_calo;icalo++){
    
    if (!caloEnabled[icalo]) continue;
    Int_t caloDir           = GetCaloDirection(icalo);
    
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(h_TMstudies_clusSpec[icalo][ialgo])h_TMstudies_clusSpec[icalo][ialgo]->Write();
      if(h_TMstudies_clusSpec_Matched[icalo][ialgo])h_TMstudies_clusSpec_Matched[icalo][ialgo]->Write();
      if(h_TMstudies_clusSpec_MatchedCalo[icalo][ialgo])h_TMstudies_clusSpec_MatchedCalo[icalo][ialgo]->Write();
      if(h_TMstudies_2D_clus[icalo][ialgo]) h_TMstudies_2D_clus[icalo][ialgo]->Write();
      
      if(h_TMstudies_2D_delta[icalo][ialgo]) h_TMstudies_2D_delta[icalo][ialgo]->Write();
      if(h_TMstudies_dx_vsEta[icalo][ialgo]) h_TMstudies_dx_vsEta[icalo][ialgo]->Write();
      if(h_TMstudies_dy_vsEta[icalo][ialgo]) h_TMstudies_dy_vsEta[icalo][ialgo]->Write();
      
      for (int jcalo=0; jcalo<_active_calo; jcalo++){
        if (!caloEnabled[jcalo]) continue;
        if (h_TMstudiesCalo_dPhi_vsEta[icalo][jcalo][ialgo]) h_TMstudiesCalo_dPhi_vsEta[icalo][jcalo][ialgo]->Write();
        if (h_TMstudiesCalo_dEta_vsEta[icalo][jcalo][ialgo]) h_TMstudiesCalo_dEta_vsEta[icalo][jcalo][ialgo]->Write();
        if (h_TMstudiesCalo_2D_delta[icalo][jcalo][ialgo])   h_TMstudiesCalo_2D_delta[icalo][jcalo][ialgo]->Write();
      }
        
      for (Int_t et = minEtaBinCalo[caloDir]; et<maxEtaBinCalo[caloDir]; et++){
        if(h_TMstudies_clusSpec_Eta[icalo][ialgo][et])        h_TMstudies_clusSpec_Eta[icalo][ialgo][et]->Write();
        if(h_TMstudies_clusSpec_Matched_Eta[icalo][ialgo][et])h_TMstudies_clusSpec_Matched_Eta[icalo][ialgo][et]->Write();

        if(h_TMstudies_dx_vsClsE[icalo][ialgo][et])           h_TMstudies_dx_vsClsE[icalo][ialgo][et]->Write();
        if(h_TMstudies_dy_vsClsE[icalo][ialgo][et])           h_TMstudies_dy_vsClsE[icalo][ialgo][et]->Write();
        for (int jcalo=0; jcalo<_active_calo; jcalo++){
          if (h_TMstudiesCalo_dPhi_vsClsE[icalo][jcalo][ialgo][et]) h_TMstudiesCalo_dPhi_vsClsE[icalo][jcalo][ialgo][et]->Write();
          if (h_TMstudiesCalo_dEta_vsClsE[icalo][jcalo][ialgo][et]) h_TMstudiesCalo_dEta_vsClsE[icalo][jcalo][ialgo][et]->Write();
        } 
      }
    }  
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
