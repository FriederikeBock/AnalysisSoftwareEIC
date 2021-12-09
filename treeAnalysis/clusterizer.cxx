#include <algorithm>
#include "caloheader.h"
// ANCHOR debug output verbosity
int verbosityCLS = 0;

// ANCHOR define global variables
float aggregation_margin_V3 = 0.03;
//                                       FHCAL   FEMC    DRCALO  EEMC   CEMC   EHCAL HCALIN  HCALOUT LFHCAL  EEMCG BECAL
float aggregation_margin_MA[maxcalo] = { 0.05,   0.03,   0.03,   0.075, 0.05,  0.1,  0.1,    0.1,    0.03,   0.05,  0.03   };

TH2F*  h_clusterizer_nonagg_towers[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_matched_2D_delta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_all_2D_delta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_all_2D_deltaCalo[_active_calo][_active_calo][_active_algo]; // [calorimeter_enum][calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_calomatch_2D_deltaCalo[_active_calo][_active_calo][_active_algo]; // [calorimeter_enum][calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_calomatchtrue_2D_deltaCalo[_active_calo][_active_calo][_active_algo]; // [calorimeter_enum][calorimeter_enum][algorithm_enum]

TH2F*  h_clusterizer_isMatched_E[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_isMatchedCalo_E[_active_calo][_active_calo][_active_algo]; // [calorimeter_enum][calorimeter_enum][algorithm_enum]


bool _doClusterECalibration = true;

//**************************************************************************************************************
//**************************************************************************************************************
// ANCHOR track/projection matching function
//**************************************************************************************************************
//**************************************************************************************************************
std::vector<matchingStrct> isClusterMatched(  clustersStrct tempcluster, 
                        int caloEnum, 
                        int clusterizerEnum, 
                        unsigned short primaryTrackSource, 
                        bool useProjection = true
                     ){
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    bool isFwd = IsForwardCalorimeter(icalo);
    float nbins2Ddelta = isFwd ? 78*4 : 80*2;
    float min2Ddeltahist = isFwd ? -39 : -0.2;
    float max2Ddeltahist = isFwd ? 39 : 0.2;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!h_clusterizer_all_2D_delta[icalo][ialgo])
        h_clusterizer_all_2D_delta[icalo][ialgo]         = new TH2F(Form("h_clusterizer_all_2D_delta_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 
                                                                    nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
      if(!h_clusterizer_matched_2D_delta[icalo][ialgo])
        h_clusterizer_matched_2D_delta[icalo][ialgo]     = new TH2F(Form("h_clusterizer_matched_2D_delta_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 
                                                                    nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
    }
  }

  std::vector<matchingStrct> matching_trackIDs;

  // in case the cluster position is 0,0 -> no matching to be done
  if(tempcluster.cluster_X==0 && tempcluster.cluster_Y==0) return matching_trackIDs;

  matchingStrct tempMatch;
  float matchingwindow = ReturnTrackMatchingWindowForCalo(caloEnum);
  bool isFwd = IsForwardCalorimeter(caloEnum);

  if(useProjection){
    int projectionlayer = ReturnProjectionIndexForCalorimeter(caloEnum, _useAlternateForProjections);    // use the projections to the last TTL layers
//     std::cout << caloEnum << "\t" << projectionlayer << endl;
    if(projectionlayer==-1) return matching_trackIDs;
    for(Int_t iproj=0; iproj<_nProjections; iproj++){
      // check for correct projection layer
      if(_track_ProjLayer[iproj]!=projectionlayer) continue;
      // bad timing = bad projection
      if(_track_Proj_t[iproj] < 0) continue;
      if (_track_source[(int)_track_ProjTrackID[iproj]] != primaryTrackSource) { continue; }
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

      float dPhi = projphi-tempcluster.cluster_Phi;
      if (dPhi >= TMath::Pi()) dPhi = dPhi-2*TMath::Pi();
      if (dPhi <= -TMath::Pi()) dPhi = dPhi+2*TMath::Pi();

      
      // calculate delta x/y or delta eta/phi between projection and cluster
      float delta_1 = isFwd ? _track_Proj_x[iproj]-tempcluster.cluster_X : projeta-tempcluster.cluster_Eta;
      float delta_2 = isFwd ? _track_Proj_y[iproj]-tempcluster.cluster_Y : dPhi;

      // fill delta histogram
      h_clusterizer_all_2D_delta[caloEnum][clusterizerEnum]->Fill(delta_1,delta_2);

      // check x or eta difference
      if(abs(delta_1) < matchingwindow){
        // check y or phi difference
        if(abs(delta_2) < matchingwindow){
          tempMatch.dX = _track_Proj_x[iproj]-tempcluster.cluster_X;
          tempMatch.dY = _track_Proj_y[iproj]-tempcluster.cluster_Y;
          tempMatch.dPhi = dPhi;
          tempMatch.dEta = projeta-tempcluster.cluster_Eta;
          tempMatch.id = (int)_track_ProjTrackID[iproj];
          h_clusterizer_matched_2D_delta[caloEnum][clusterizerEnum]->Fill(delta_1,delta_2);
          matching_trackIDs.push_back(tempMatch);
          // reset temp struct 
          tempMatch = matchingStrct();
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
      tracketa = trackvec.Eta();
      trackphi = trackvec.Phi();
      if(isFwd){
        trackvec*=(zHC/abs(trackvec.Z()));
      } else {
          if((tracketa==0 && trackphi==0) || (tracketa==-20 && trackphi==-20)) continue;
      }

      float dPhi = trackphi-tempcluster.cluster_Phi;
      if (dPhi >= TMath::Pi()) dPhi = dPhi-2*TMath::Pi();
      if (dPhi <= -TMath::Pi()) dPhi = dPhi+2*TMath::Pi();

      // calculate delta x/y or delta eta/phi between projection and cluster
      float delta_1 = isFwd ? trackvec.X()-tempcluster.cluster_X : tracketa-tempcluster.cluster_Eta;
      float delta_2 = isFwd ? trackvec.Y()-tempcluster.cluster_Y : trackphi-tempcluster.cluster_Phi;

//				cout << "CLUST  \t TrackID " << itrk << " Calo " << caloEnum << " dR = " << sqrt(delta_1*delta_1 + delta_2*delta_2) << endl;
      // fill delta histogram
      h_clusterizer_all_2D_delta[caloEnum][clusterizerEnum]->Fill(delta_1,delta_2);

      // check x or eta difference
      if(abs(delta_1) < matchingwindow){
        // check y or phi difference
        if(abs(delta_2) < matchingwindow){
          tempMatch.dX    = trackvec.X()-tempcluster.cluster_X;
          tempMatch.dY    = trackvec.Y()-tempcluster.cluster_Y;
          tempMatch.dPhi  = dPhi;
          tempMatch.dEta  = tracketa-tempcluster.cluster_Eta;
          tempMatch.id    = itrk;
          h_clusterizer_matched_2D_delta[caloEnum][clusterizerEnum]->Fill(delta_1,delta_2);
          matching_trackIDs.push_back(tempMatch);
          // reset temp struct 
          tempMatch = matchingStrct();
        }
      }
    }
  }
  return matching_trackIDs;
}

//**************************************************************************************************************
//**************************************************************************************************************
// associate clusters in vectors to track arrays
//**************************************************************************************************************
//**************************************************************************************************************
void SetClustersMatchedToTracks(){
  for (int caloEnum = 0; caloEnum < maxcalo; caloEnum++){
    if (!caloEnabled[caloEnum]) continue;
    if (caloEnum == kHCALIN) continue; // skip HCAL in for now in track association
    for (int ialgo = 0; ialgo < _active_algo; ialgo++){
      // create bookkeeping histo
      if(!h_clusterizer_isMatched_E[caloEnum][ialgo]){
        h_clusterizer_isMatched_E[caloEnum][ialgo]             = new TH2F(Form("h_clusterizer_isMatched_E_%s_%s", str_calorimeter[caloEnum].Data(), str_clusterizer[ialgo].Data()), "; condition; E (GeV)", 3, 0, 3,100,0,100);
        h_clusterizer_isMatched_E[caloEnum][ialgo]->GetXaxis()->SetBinLabel(1, "clusters");
        h_clusterizer_isMatched_E[caloEnum][ialgo]->GetXaxis()->SetBinLabel(2, "clusters 1 match");
        h_clusterizer_isMatched_E[caloEnum][ialgo]->GetXaxis()->SetBinLabel(3, "clusters > 1 match");
      }
      
      for (int iCl = 0; iCl < (int)_clusters_calo[ialgo][caloEnum].size(); iCl++){
        h_clusterizer_isMatched_E[caloEnum][ialgo]->Fill(0.5,(_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_E );
        
        if ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_isMatched){
          if ((((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).size()) == 1)
            h_clusterizer_isMatched_E[caloEnum][ialgo]->Fill(1.5,(_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_E );
          else 
            h_clusterizer_isMatched_E[caloEnum][ialgo]->Fill(2.5,(_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_E );
            
          for (int iT = 0; iT < (int)(((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).size()); iT++){
            int trkID = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).id;
            // do separate checks for HCal and ECal
            if (IsHCALCalorimeter(caloEnum)){
              // is this the first match?
              if (_track_matchHCal[trkID].caloid == -1){
                _track_matchHCal[trkID].caloid  = caloEnum;
                _track_matchHCal[trkID].id      = iCl;
                _track_matchHCal[trkID].dPhi    = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dPhi;
                _track_matchHCal[trkID].dEta    = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dEta;
                _track_matchHCal[trkID].dX      = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dX;
                _track_matchHCal[trkID].dY      = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dY;
              } else {
                // calc 2D distance for current cluster
                float dRC = TMath::Sqrt(TMath::Power(((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dPhi,2)+TMath::Power(((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dEta,2));
                // calc 2D distance for already matched cluster
                float dRP = TMath::Sqrt(TMath::Power(_track_matchHCal[trkID].dPhi,2)+ TMath::Power(_track_matchHCal[trkID].dEta,2));
                // update track info if cluster better matched
                if (dRC < dRP){
                  _track_matchHCal[trkID].caloid  = caloEnum;
                  _track_matchHCal[trkID].id      = iCl;
                  _track_matchHCal[trkID].dPhi    = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dPhi;
                  _track_matchHCal[trkID].dEta    = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dEta;
                  _track_matchHCal[trkID].dX      = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dX;
                  _track_matchHCal[trkID].dY      = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dY;
                }
              }
            } else {
              // is this the first match?
              if (_track_matchECal[trkID].caloid == -1){
                _track_matchECal[trkID].caloid  = caloEnum;
                _track_matchECal[trkID].id      = iCl;
                _track_matchECal[trkID].dPhi    = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dPhi;
                _track_matchECal[trkID].dEta    = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dEta;
                _track_matchECal[trkID].dX      = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dX;
                _track_matchECal[trkID].dY      = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dY;
              } else {
                // calc 2D distance for current cluster
                float dRC = TMath::Sqrt(TMath::Power(((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dPhi,2)+ TMath::Power(((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dEta,2));
                // calc 2D distance for already matched cluster
                float dRP = TMath::Sqrt(TMath::Power(_track_matchECal[trkID].dPhi,2)+ TMath::Power(_track_matchECal[trkID].dEta,2));
                // update track info if cluster better matched
                if (dRC < dRP){
                  _track_matchECal[trkID].caloid  = caloEnum;
                  _track_matchECal[trkID].id      = iCl;
                  _track_matchECal[trkID].dPhi    = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dPhi;
                  _track_matchECal[trkID].dEta    = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dEta;
                  _track_matchECal[trkID].dX      = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dX;
                  _track_matchECal[trkID].dY      = ((_clusters_calo[ialgo][caloEnum].at(iCl)).cluster_matchedTracks).at(iT).dY;
                }
              }
            }
          }
        }
      }
    }
  }
}


//**************************************************************************************************************
//**************************************************************************************************************
// associate HCal and Ecal clusters with each other
//**************************************************************************************************************
//**************************************************************************************************************
void MatchClustersWithOtherCalos( int caloEnum, int clusterizerEnum ){
  
  if (!caloEnabled[caloEnum]) return;
  if( !(_combCalo[caloEnum] != -1  && IsHCALCalorimeter(caloEnum))) return;
  if (!caloEnabled[_combCalo[caloEnum]]) return;
  int calID2 = _combCalo[caloEnum];
  
//   std::cout <<  "starting matching for "<< str_calorimeter[caloEnum].Data() << " and " << str_calorimeter[calID2].Data() << std::endl;
  
  if (!h_clusterizer_all_2D_deltaCalo[caloEnum][calID2][clusterizerEnum])
    h_clusterizer_all_2D_deltaCalo[caloEnum][calID2][clusterizerEnum]           = new TH2F(Form("h_clusterizer_all_2D_deltaCalo_%s_%s_%s",str_calorimeter[caloEnum].Data(),str_calorimeter[calID2].Data(), str_clusterizer[clusterizerEnum].Data()), "; dEta; dPhi", 200, -0.25, 0.25, 200, -0.25, 0.25);
  if (!h_clusterizer_calomatch_2D_deltaCalo[caloEnum][calID2][clusterizerEnum])
    h_clusterizer_calomatch_2D_deltaCalo[caloEnum][calID2][clusterizerEnum]           = new TH2F(Form("h_clusterizer_calomatch_2D_deltaCalo_%s_%s_%s",str_calorimeter[caloEnum].Data(),str_calorimeter[calID2].Data(), str_clusterizer[clusterizerEnum].Data()), "; dEta; dPhi", 200, -0.25, 0.25, 200, -0.25, 0.25);
  if (!h_clusterizer_calomatchtrue_2D_deltaCalo[caloEnum][calID2][clusterizerEnum])
    h_clusterizer_calomatchtrue_2D_deltaCalo[caloEnum][calID2][clusterizerEnum]           = new TH2F(Form("h_clusterizer_calomatch_true_2D_deltaCalo_%s_%s_%s",str_calorimeter[caloEnum].Data(),str_calorimeter[calID2].Data(), str_clusterizer[clusterizerEnum].Data()), "; dEta; dPhi", 200, -0.25, 0.25, 200, -0.25, 0.25);

  // create bookkeeping hist
  if(!h_clusterizer_isMatchedCalo_E[caloEnum][calID2][clusterizerEnum]){
    h_clusterizer_isMatchedCalo_E[caloEnum][calID2][clusterizerEnum]             = new TH2F(Form("h_clusterizer_isMatchedCalo_E_%s_%s_%s", str_calorimeter[caloEnum].Data(), str_calorimeter[calID2].Data(),  str_clusterizer[clusterizerEnum].Data()), "; condition; E (GeV)", 3, 0, 3,100,0,100);  
    h_clusterizer_isMatchedCalo_E[caloEnum][calID2][clusterizerEnum]->GetXaxis()->SetBinLabel(1, "clusters");
    h_clusterizer_isMatchedCalo_E[caloEnum][calID2][clusterizerEnum]->GetXaxis()->SetBinLabel(2, "clusters 1 match");
    h_clusterizer_isMatchedCalo_E[caloEnum][calID2][clusterizerEnum]->GetXaxis()->SetBinLabel(3, "clusters > 1 match");
  }

  for (int iCl = 0; iCl < (int)_clusters_calo[clusterizerEnum][caloEnum].size(); iCl++){
    float etaCl1 = (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_Eta;
    float phiCl1 = (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_Phi;
    float xCl1 = (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_X;
    float yCl1 = (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_Y;
    
    if (xCl1 == 0 && yCl1 == 0) continue;
    
    h_clusterizer_isMatchedCalo_E[caloEnum][calID2][clusterizerEnum]->Fill(0.5,(_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_E );
    int nMatches = 0;
    for (int jCl = 0; jCl < (int)_clusters_calo[clusterizerEnum][calID2].size(); jCl++){
      float etaCl2 = (_clusters_calo[clusterizerEnum][calID2].at(jCl)).cluster_Eta;
      float phiCl2 = (_clusters_calo[clusterizerEnum][calID2].at(jCl)).cluster_Phi;
      float xCl2 = (_clusters_calo[clusterizerEnum][calID2].at(jCl)).cluster_X;
      float yCl2 = (_clusters_calo[clusterizerEnum][calID2].at(jCl)).cluster_Y;
      if (xCl2 == 0 && yCl2 == 0) continue;
        
      float dEta  = etaCl1 -etaCl2;      
      float dPhi  = phiCl1 -phiCl2;  // phi is periodic (correct for that)
      if (dPhi >= TMath::Pi()) dPhi = dPhi-2*TMath::Pi();
      if (dPhi <= -TMath::Pi()) dPhi = dPhi+2*TMath::Pi();
      float dX    = xCl1 -xCl2;
      float dY    = yCl1 -yCl2;
    
      h_clusterizer_all_2D_deltaCalo[caloEnum][calID2][clusterizerEnum]->Fill(dEta, dPhi);
      
      if (TMath::Abs(dPhi) > ReturnPhiCaloMatching(caloEnum, calID2)) continue;
      if (TMath::Abs(dEta) > ReturnEtaCaloMatching(caloEnum, calID2)) continue;
      
      h_clusterizer_calomatch_2D_deltaCalo[caloEnum][calID2][clusterizerEnum]->Fill(dEta, dPhi);
      if ((_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_trueID == (_clusters_calo[clusterizerEnum][calID2].at(jCl)).cluster_trueID)
        h_clusterizer_calomatchtrue_2D_deltaCalo[caloEnum][calID2][clusterizerEnum]->Fill(dEta, dPhi);
      
      nMatches++;
      
      if (caloEnum == kEHCAL || caloEnum == kLFHCAL){
        matchingCalStrct tempMatch = matchingCalStrct();
        tempMatch.caloid  = calID2;
        tempMatch.id      = jCl;
        tempMatch.dPhi    = dPhi;
        tempMatch.dEta    = dEta;
        tempMatch.dX      = dX;
        tempMatch.dY      = dY;
        (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_matchedECals.push_back(tempMatch);
        tempMatch.caloid  = caloEnum;
        tempMatch.id      = iCl;
        (_clusters_calo[clusterizerEnum][calID2].at(jCl)).cluster_matchedHCals.push_back(tempMatch);
      }
      
      if (caloEnum == kHCALOUT){
        matchingCalStrct tempMatch = matchingCalStrct();
        tempMatch.caloid  = calID2;
        tempMatch.id      = jCl;
        tempMatch.dPhi    = dPhi;
        tempMatch.dEta    = dEta;
        tempMatch.dX      = dX;
        tempMatch.dY      = dY;
        (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_matchedECals.push_back(tempMatch);
        tempMatch.caloid  = caloEnum;
        tempMatch.id      = iCl;
        (_clusters_calo[clusterizerEnum][calID2].at(jCl)).cluster_matchedHCals.push_back(tempMatch);

      }

      if (caloEnum == kHCALIN){
        matchingCalStrct tempMatch = matchingCalStrct();
        tempMatch.caloid  = calID2;
        tempMatch.id      = jCl;
        tempMatch.dPhi    = dPhi;
        tempMatch.dEta    = dEta;
        tempMatch.dX      = dX;
        tempMatch.dY      = dY;
        (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_matchedECals.push_back(tempMatch);
        tempMatch.caloid  = caloEnum;
        tempMatch.id      = iCl;
        (_clusters_calo[clusterizerEnum][calID2].at(jCl)).cluster_matchedECals.push_back(tempMatch);
      } 
    }
    if (nMatches == 1)
      h_clusterizer_isMatchedCalo_E[caloEnum][calID2][clusterizerEnum]->Fill(1.5,(_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_E );
    else if (nMatches > 1)
      h_clusterizer_isMatchedCalo_E[caloEnum][calID2][clusterizerEnum]->Fill(2.5,(_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_E );
  }

  if (caloEnum == kHCALOUT){
    if( _combCalo2[caloEnum] == -1 ) return;
    if (!caloEnabled[_combCalo2[caloEnum]]) return;
    int calID3 = _combCalo2[caloEnum];
    
    if (!h_clusterizer_all_2D_deltaCalo[caloEnum][calID3][clusterizerEnum])
      h_clusterizer_all_2D_deltaCalo[caloEnum][calID3][clusterizerEnum]           = new TH2F(Form("h_clusterizer_all_2D_deltaCalo_%s_%s_%s",str_calorimeter[caloEnum].Data(),str_calorimeter[calID3].Data(), str_clusterizer[clusterizerEnum].Data()), "; dEta; dPhi", 200, -0.25, 0.25, 200, -0.25, 0.25);
    if (!h_clusterizer_calomatch_2D_deltaCalo[caloEnum][calID3][clusterizerEnum])
      h_clusterizer_calomatch_2D_deltaCalo[caloEnum][calID3][clusterizerEnum]           = new TH2F(Form("h_clusterizer_calomatch_2D_deltaCalo_%s_%s_%s",str_calorimeter[caloEnum].Data(),str_calorimeter[calID3].Data(), str_clusterizer[clusterizerEnum].Data()), "; dEta; dPhi", 200, -0.25, 0.25, 200, -0.25, 0.25);
    if (!h_clusterizer_calomatchtrue_2D_deltaCalo[caloEnum][calID3][clusterizerEnum])
      h_clusterizer_calomatchtrue_2D_deltaCalo[caloEnum][calID3][clusterizerEnum]           = new TH2F(Form("h_clusterizer_calomatch_true_2D_deltaCalo_%s_%s_%s",str_calorimeter[caloEnum].Data(),str_calorimeter[calID3].Data(), str_clusterizer[clusterizerEnum].Data()), "; dEta; dPhi", 200, -0.25, 0.25, 200, -0.25, 0.25);

      // create bookkeeping hist
    if(!h_clusterizer_isMatchedCalo_E[caloEnum][calID3][clusterizerEnum]){
      h_clusterizer_isMatchedCalo_E[caloEnum][calID3][clusterizerEnum]             = new TH2F(Form("h_clusterizer_isMatchedCalo_E_%s_%s_%s", str_calorimeter[caloEnum].Data(), str_calorimeter[calID3].Data(),  str_clusterizer[clusterizerEnum].Data()), "; condition; E (GeV)", 3, 0, 3,100,0,100);  
      h_clusterizer_isMatchedCalo_E[caloEnum][calID3][clusterizerEnum]->GetXaxis()->SetBinLabel(1, "clusters");
      h_clusterizer_isMatchedCalo_E[caloEnum][calID3][clusterizerEnum]->GetXaxis()->SetBinLabel(2, "clusters 1 match");
      h_clusterizer_isMatchedCalo_E[caloEnum][calID3][clusterizerEnum]->GetXaxis()->SetBinLabel(3, "clusters > 1 match");
    }

//     std::cout <<  "starting matching for "<< str_calorimeter[caloEnum].Data() << " and " << str_calorimeter[calID3].Data() << std::endl;
    
    for (int iCl = 0; iCl < (int)_clusters_calo[clusterizerEnum][caloEnum].size(); iCl++){
      float etaCl1 = (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_Eta;
      float phiCl1 = (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_Phi;
      float xCl1 = (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_X;
      float yCl1 = (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_Y;
      
      h_clusterizer_isMatchedCalo_E[caloEnum][calID3][clusterizerEnum]->Fill(0.5,(_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_E );
      int nMatches = 0;
      
      if (xCl1 == 0 && yCl1 == 0) continue;
      
      for (int jCl = 0; jCl < (int)_clusters_calo[clusterizerEnum][calID3].size(); jCl++){
        float etaCl2 = (_clusters_calo[clusterizerEnum][calID3].at(jCl)).cluster_Eta;
        float phiCl2 = (_clusters_calo[clusterizerEnum][calID3].at(jCl)).cluster_Phi;
        float xCl2 = (_clusters_calo[clusterizerEnum][calID3].at(jCl)).cluster_X;
        float yCl2 = (_clusters_calo[clusterizerEnum][calID3].at(jCl)).cluster_Y;
        if (xCl2 == 0 && yCl2 == 0) continue;
          
        float dEta  = etaCl1 -etaCl2;      
        float dPhi  = phiCl1 -phiCl2;  // phi is periodic (correct for that)
        if (dPhi >= TMath::Pi()) dPhi = dPhi-2*TMath::Pi();
        if (dPhi <= -TMath::Pi()) dPhi = dPhi+2*TMath::Pi();
        float dX    = xCl1 -xCl2;
        float dY    = yCl1 -yCl2;
      
        h_clusterizer_all_2D_deltaCalo[caloEnum][calID3][clusterizerEnum]->Fill(dEta, dPhi);
      
        if (TMath::Abs(dPhi) > ReturnPhiCaloMatching(caloEnum, calID3)) continue;
        if (TMath::Abs(dEta) > ReturnEtaCaloMatching(caloEnum, calID3)) continue;
        nMatches++;
        h_clusterizer_calomatch_2D_deltaCalo[caloEnum][calID3][clusterizerEnum]->Fill(dEta, dPhi);
        if ((_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_trueID == (_clusters_calo[clusterizerEnum][calID3].at(jCl)).cluster_trueID)
          h_clusterizer_calomatchtrue_2D_deltaCalo[caloEnum][calID3][clusterizerEnum]->Fill(dEta, dPhi);

        matchingCalStrct tempMatch = matchingCalStrct();
        tempMatch.caloid  = calID3;
        tempMatch.id      = jCl;
        tempMatch.dEta    = dEta;
        tempMatch.dPhi    = dPhi;
        tempMatch.dX      = dX;
        tempMatch.dY      = dY;
        (_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_matchedHCals.push_back(tempMatch);
        tempMatch.caloid  = caloEnum;
        tempMatch.id      = iCl;
        (_clusters_calo[clusterizerEnum][calID3].at(jCl)).cluster_matchedHCals.push_back(tempMatch);      }
      if (nMatches == 1)
        h_clusterizer_isMatchedCalo_E[caloEnum][calID3][clusterizerEnum]->Fill(1.5,(_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_E );
      else if (nMatches > 1)
        h_clusterizer_isMatchedCalo_E[caloEnum][calID3][clusterizerEnum]->Fill(2.5,(_clusters_calo[clusterizerEnum][caloEnum].at(iCl)).cluster_E );

    }
  }
}


//**************************************************************************************************************
//**************************************************************************************************************
// put calorimeter towers in temporary structure for clusterization
//**************************************************************************************************************
//**************************************************************************************************************
std::vector<towersStrct> readTowersForCalo( int caloEnum, float aggE, float escaling = 1){
  std::vector<towersStrct> input_towers_temp;
  if(caloEnum==kFHCAL){
    if (verbosityCLS > 1) std::cout << "FHCAL: "<<_nTowers_FHCAL << std::endl;
    for(int itow=0; itow<_nTowers_FHCAL; itow++){
      if(_tower_FHCAL_E[itow]>(aggE/escaling)){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_FHCAL_E[itow]*escaling;
        tempstructT.tower_iEta    = _tower_FHCAL_iEta[itow];
        tempstructT.tower_iPhi    = _tower_FHCAL_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_FHCAL_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kFEMC){
    if (verbosityCLS > 1) std::cout << "FEMC: "<<_nTowers_FEMC << std::endl;
    for(int itow=0; itow<_nTowers_FEMC; itow++){
      if(_tower_FEMC_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_FEMC_E[itow]; 
        tempstructT.tower_iEta    = _tower_FEMC_iEta[itow];
        tempstructT.tower_iPhi    = _tower_FEMC_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_FEMC_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kCEMC){
    if (verbosityCLS > 1) std::cout << "CEMC: "<< _nTowers_CEMC << std::endl;
    for(int itow=0; itow<_nTowers_CEMC; itow++){
      if(_tower_CEMC_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_CEMC_E[itow];
        tempstructT.tower_iEta    = _tower_CEMC_iEta[itow];
        tempstructT.tower_iPhi    = _tower_CEMC_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_CEMC_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kHCALIN){
    if (verbosityCLS > 1) std::cout << "HCALIN: "<< _nTowers_HCALIN << std::endl;
    for(int itow=0; itow<_nTowers_HCALIN; itow++){
      if(_tower_HCALIN_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_HCALIN_E[itow];
        tempstructT.tower_iEta    = _tower_HCALIN_iEta[itow];
        tempstructT.tower_iPhi    = _tower_HCALIN_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_HCALIN_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kHCALOUT){
    if (verbosityCLS > 1) std::cout << "HCALOUT: "<< _nTowers_HCALOUT << std::endl;
    for(int itow=0; itow<_nTowers_HCALOUT; itow++){
      if(_tower_HCALOUT_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_HCALOUT_E[itow];
        tempstructT.tower_iEta    = _tower_HCALOUT_iEta[itow];
        tempstructT.tower_iPhi    = _tower_HCALOUT_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_HCALOUT_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kLFHCAL){
    if (verbosityCLS > 1) std::cout << "LFHCAL: "<< _nTowers_LFHCAL << std::endl;
    for(int itow=0; itow<_nTowers_LFHCAL; itow++){
      if(_tower_LFHCAL_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_LFHCAL_E[itow]; 
        tempstructT.tower_iEta    = _tower_LFHCAL_iEta[itow];
        tempstructT.tower_iPhi    = _tower_LFHCAL_iPhi[itow];
        tempstructT.tower_iL      = _tower_LFHCAL_iL[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_LFHCAL_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kDRCALO){
    if (verbosityCLS > 1) std::cout << "DRCALO: "<< _nTowers_DRCALO << std::endl;
    for(int itow=0; itow<_nTowers_DRCALO; itow++){
      if(_tower_DRCALO_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_DRCALO_E[itow];
        tempstructT.tower_iEta    = _tower_DRCALO_iEta[itow];
        tempstructT.tower_iPhi    = _tower_DRCALO_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_DRCALO_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kEHCAL){
    if (verbosityCLS > 1) std::cout << "EHCAL: "<< _nTowers_EHCAL << std::endl;
    for(int itow=0; itow<_nTowers_EHCAL; itow++){
      if(_tower_EHCAL_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_EHCAL_E[itow];
        tempstructT.tower_iEta    = _tower_EHCAL_iEta[itow];
        tempstructT.tower_iPhi    = _tower_EHCAL_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_EHCAL_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kBECAL){
    if (verbosityCLS > 1) std::cout << "BECAL: "<< _nTowers_BECAL << std::endl;
    for(int itow=0; itow<_nTowers_BECAL; itow++){
      if(_tower_BECAL_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_BECAL_E[itow];
        tempstructT.tower_iEta    = _tower_BECAL_iEta[itow];
        tempstructT.tower_iPhi    = _tower_BECAL_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_BECAL_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kEEMC){
    if (verbosityCLS > 1) std::cout << "EEMC: "<< _nTowers_EEMC << std::endl;
    for(int itow=0; itow<_nTowers_EEMC; itow++){
      if(_tower_EEMC_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_EEMC_E[itow];
        tempstructT.tower_iEta    = _tower_EEMC_iEta[itow];
        tempstructT.tower_iPhi    = _tower_EEMC_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_EEMC_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else if(caloEnum==kEEMCG){
    if (verbosityCLS > 1) std::cout << "EEMCG: "<< _nTowers_EEMCG << std::endl;
    for(int itow=0; itow<_nTowers_EEMCG; itow++){
      if(_tower_EEMCG_E[itow]>aggE/escaling){
        towersStrct tempstructT;
        tempstructT.tower_E       = _tower_EEMCG_E[itow];
        tempstructT.tower_iEta    = _tower_EEMCG_iEta[itow];
        tempstructT.tower_iPhi    = _tower_EEMCG_iPhi[itow];
        tempstructT.tower_trueID  = GetCorrectMCArrayEntry(_tower_EEMCG_trueID[itow]);
        input_towers_temp.push_back(tempstructT);
      }
    }
  } else {
    std::cout << "Incorrect calorimeter selected! Enum " << caloEnum << " not defined!" << std::endl;
  }
  return input_towers_temp;
}

//**************************************************************************************************************
//**************************************************************************************************************
// find clusters in a circle of X-cells around leading tower
//**************************************************************************************************************
//**************************************************************************************************************
clustersStrct findCircularCluster(
                                    int size,                                       // size of circle in cells
                                    float seed,                                     // minimum seed energy
                                    int caloEnum,                                   // calorimeter type
                                    std::vector<towersStrct> &input_towers_temp,    // temporary full tower array
                                    std::vector<towersStrct> &cluster_towers_temp,  // towers associated to cluster
                                    std::vector<int> clslabels_temp                 // MC labels in cluster
                                 ){
  int maxdiff = 0;
  if (size == 3)
    maxdiff = 1;
  else if (size == 5)
    maxdiff = 2;

  clustersStrct tempstructC;
  
  if(input_towers_temp.at(0).tower_E > seed){
    if (verbosityCLS > 2) std::cout << "new cluster" << std::endl;
    // fill seed cell information into current cluster
    tempstructC.cluster_E       = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_seed    = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_NTowers = 1;
    tempstructC.cluster_NtrueID = 1;
    tempstructC.cluster_trueID = input_towers_temp.at(0).tower_trueID; // TODO save all MC labels?
    cluster_towers_temp.push_back(input_towers_temp.at(0));
    clslabels_temp.push_back(input_towers_temp.at(0).tower_trueID);
    if (verbosityCLS > 2) std::cout << "seed: "<<  input_towers_temp.at(0).tower_iEta << "\t" << input_towers_temp.at(0).tower_iPhi << "\t" << input_towers_temp.at(0).tower_iL << "\t E:"<< tempstructC.cluster_E << std::endl;
    
    
    for (int tit = 1; tit < (int)input_towers_temp.size(); tit++){
      // towers must be within cross of delta Eta and delta Phi <= 1
      int iPhiRef = input_towers_temp.at(0).tower_iPhi;
      int iEtaRef = input_towers_temp.at(0).tower_iEta;
      int iPhiAgg = input_towers_temp.at(tit).tower_iPhi;
      int iEtaAgg = input_towers_temp.at(tit).tower_iEta;
      // ensure clusters aren't split in barrel calorimeters at pi/-pi
      if (!IsForwardCalorimeter(caloEnum)) {
        if (iPhiRef < 5 && iPhiAgg > _caloTowersPhi[caloEnum]-5){
          iPhiAgg= iPhiAgg-_caloTowersPhi[caloEnum];
        }
        if (iPhiRef > _caloTowersPhi[caloEnum]-5 && iPhiAgg < 5){
          iPhiRef= iPhiRef-_caloTowersPhi[caloEnum];
        }
      }
      int deltaEta = std::abs(iEtaAgg-iEtaRef);
      int deltaPhi = std::abs(iPhiAgg-iPhiRef);
      if( ( deltaEta + deltaPhi ) <= maxdiff ){
        tempstructC.cluster_E+=input_towers_temp.at(tit).tower_E;
        tempstructC.cluster_NTowers++;
        cluster_towers_temp.push_back(input_towers_temp.at(tit));
        if(!(std::find(clslabels_temp.begin(), clslabels_temp.end(), input_towers_temp.at(tit).tower_trueID) != clslabels_temp.end())){
          tempstructC.cluster_NtrueID++;
          clslabels_temp.push_back(input_towers_temp.at(tit).tower_trueID);
        }
        input_towers_temp.erase(input_towers_temp.begin()+tit);
        tit--;
      }
    }
  } 
  input_towers_temp.erase(input_towers_temp.begin());
  return tempstructC;
}

//**************************************************************************************************************
//**************************************************************************************************************
// find clusters in a square of X-cells around leading tower
//**************************************************************************************************************
//**************************************************************************************************************
clustersStrct findSquareCluster(
                                  int size,                                       // size of square in cells
                                  float seed,                                     // minimum seed energy
                                  int caloEnum,                                   // calorimeter type
                                  std::vector<towersStrct> &input_towers_temp,    // temporary full tower array
                                  std::vector<towersStrct> &cluster_towers_temp,  // towers associated to cluster
                                  std::vector<int> clslabels_temp                 // MC labels in cluster
                               ){
  int maxdiff = 0;
  if (size == 3)
    maxdiff = 2;
  else if (size == 5)
    maxdiff = 3;

  clustersStrct tempstructC;
  
  if(input_towers_temp.at(0).tower_E > seed){
    if (verbosityCLS > 2) std::cout << "new cluster" << std::endl;
    // fill seed cell information into current cluster
    tempstructC.cluster_E       = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_seed    = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_NTowers = 1;
    tempstructC.cluster_NtrueID = 1;
    tempstructC.cluster_trueID = input_towers_temp.at(0).tower_trueID; // TODO save all MC labels?
    cluster_towers_temp.push_back(input_towers_temp.at(0));
    clslabels_temp.push_back(input_towers_temp.at(0).tower_trueID);
    if (verbosityCLS > 2) std::cout << "seed: "<<  input_towers_temp.at(0).tower_iEta << "\t" << input_towers_temp.at(0).tower_iPhi << "\t" << input_towers_temp.at(0).tower_iL << "\t E:"<< tempstructC.cluster_E << std::endl;
    
    
    for (int tit = 1; tit < (int)input_towers_temp.size(); tit++){
      // towers must be within cross of delta Eta and delta Phi <= 1
      int iPhiRef = input_towers_temp.at(0).tower_iPhi;
      int iEtaRef = input_towers_temp.at(0).tower_iEta;
      int iPhiAgg = input_towers_temp.at(tit).tower_iPhi;
      int iEtaAgg = input_towers_temp.at(tit).tower_iEta;
      // ensure clusters aren't split in barrel calorimeters at pi/-pi
      if (!IsForwardCalorimeter(caloEnum)) {
        if (iPhiRef < 5 && iPhiAgg > _caloTowersPhi[caloEnum]-5){
          iPhiAgg= iPhiAgg-_caloTowersPhi[caloEnum];
        }
        if (iPhiRef > _caloTowersPhi[caloEnum]-5 && iPhiAgg < 5){
          iPhiRef= iPhiRef-_caloTowersPhi[caloEnum];
        }
      }
      int deltaEta = std::abs(iEtaAgg-iEtaRef);
      int deltaPhi = std::abs(iPhiAgg-iPhiRef);
      if( deltaEta < maxdiff ){
        if( deltaPhi < maxdiff ){
          tempstructC.cluster_E+=input_towers_temp.at(tit).tower_E;
          tempstructC.cluster_NTowers++;
          cluster_towers_temp.push_back(input_towers_temp.at(tit));
          if(!(std::find(clslabels_temp.begin(), clslabels_temp.end(), input_towers_temp.at(tit).tower_trueID) != clslabels_temp.end())){
            tempstructC.cluster_NtrueID++;
            clslabels_temp.push_back(input_towers_temp.at(tit).tower_trueID);
          }
          input_towers_temp.erase(input_towers_temp.begin()+tit);
          tit--;
        }
      }
    }
  }
  input_towers_temp.erase(input_towers_temp.begin());
  return tempstructC;
}

//**************************************************************************************************************
//**************************************************************************************************************
// find clusters with continuous energy distribution and common edges
//**************************************************************************************************************
//**************************************************************************************************************
clustersStrct findV1Cluster(
                              float seed,                                     // minimum seed energy
                              int caloEnum,                                   // calorimeter type
                              std::vector<towersStrct> &input_towers_temp,    // temporary full tower array
                              std::vector<towersStrct> &cluster_towers_temp,  // towers associated to cluster
                              std::vector<int> clslabels_temp                 // MC labels in cluster
                           ){

  clustersStrct tempstructC;
  
  if(input_towers_temp.at(0).tower_E > seed){
    if (verbosityCLS > 2) std::cout << "new cluster" << std::endl;
    // fill seed cell information into current cluster
    tempstructC.cluster_E       = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_seed    = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_NTowers = 1;
    tempstructC.cluster_NtrueID = 1;
    tempstructC.cluster_trueID = input_towers_temp.at(0).tower_trueID; // TODO save all MC labels?
    cluster_towers_temp.push_back(input_towers_temp.at(0));
    clslabels_temp.push_back(input_towers_temp.at(0).tower_trueID);
    if (verbosityCLS > 2) std::cout << "seed: "<<  input_towers_temp.at(0).tower_iEta << "\t" << input_towers_temp.at(0).tower_iPhi << "\t" << input_towers_temp.at(0).tower_iL << "\t E:"<< tempstructC.cluster_E << std::endl;
    
    
    // remove seed tower from sample
    input_towers_temp.erase(input_towers_temp.begin());
    for (int tit = 0; tit < (int)cluster_towers_temp.size(); tit++){
      // Now go recursively to all neighbours and add them to the cluster if they fulfill the conditions
      int iEtaTwr = cluster_towers_temp.at(tit).tower_iEta;
      int iPhiTwr = cluster_towers_temp.at(tit).tower_iPhi;
      for (int ait = 0; ait < (int)input_towers_temp.size(); ait++){
        int iEtaTwrAgg = input_towers_temp.at(ait).tower_iEta;
        int iPhiTwrAgg = input_towers_temp.at(ait).tower_iPhi;
        
        if (!IsForwardCalorimeter(caloEnum)) {
          if (iPhiTwr < 5 && iPhiTwrAgg > _caloTowersPhi[caloEnum]-5){
            iPhiTwrAgg= iPhiTwrAgg-_caloTowersPhi[caloEnum];
          }
          if (iPhiTwr > _caloTowersPhi[caloEnum]-5 && iPhiTwrAgg < 5){
            iPhiTwr= iPhiTwr-_caloTowersPhi[caloEnum];
          }
        }
        int deltaEta = std::abs(iEtaTwrAgg-iEtaTwr);
        int deltaPhi = std::abs(iPhiTwrAgg-iPhiTwr);
        // first condition asks for V3-like neighbors, while second condition also checks diagonally attached towers
        if( ((deltaEta+deltaPhi) == 1) || (deltaEta==1 && deltaPhi==1)){
          // only aggregate towers with lower energy than current tower
          tempstructC.cluster_E+=input_towers_temp.at(ait).tower_E;
          tempstructC.cluster_NTowers++;
          cluster_towers_temp.push_back(input_towers_temp.at(ait));
          if(!(std::find(clslabels_temp.begin(), clslabels_temp.end(), input_towers_temp.at(ait).tower_trueID) != clslabels_temp.end())){
            tempstructC.cluster_NtrueID++;
            clslabels_temp.push_back(input_towers_temp.at(ait).tower_trueID);
          }
          input_towers_temp.erase(input_towers_temp.begin()+ait);
          ait--;
        }
      }
    }
  } 
  return tempstructC;
}

//**************************************************************************************************************
//**************************************************************************************************************
// find clusters with common edges, separate if energy increases in neighboring cell
//**************************************************************************************************************
//**************************************************************************************************************
clustersStrct findV3Cluster(  
                              float seed,                                     // minimum seed energy
                              int caloEnum,                                   // calorimeter type
                              std::vector<towersStrct> &input_towers_temp,    // temporary full tower array
                              std::vector<towersStrct> &cluster_towers_temp,  // towers associated to cluster
                              std::vector<int> clslabels_temp                 // MC labels in cluster
                           ){

  clustersStrct tempstructC;
  
  if(input_towers_temp.at(0).tower_E > seed){
    if (verbosityCLS > 2) std::cout << "new cluster" << std::endl;
    // fill seed cell information into current cluster
    tempstructC.cluster_E       = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_seed    = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_NTowers = 1;
    tempstructC.cluster_NtrueID = 1;
    tempstructC.cluster_trueID = input_towers_temp.at(0).tower_trueID; // TODO save all MC labels?
    cluster_towers_temp.push_back(input_towers_temp.at(0));
    clslabels_temp.push_back(input_towers_temp.at(0).tower_trueID);
    if (verbosityCLS > 2) std::cout << "seed: "<<  input_towers_temp.at(0).tower_iEta << "\t" << input_towers_temp.at(0).tower_iPhi << "\t" << input_towers_temp.at(0).tower_iL << "\t E:"<< tempstructC.cluster_E << std::endl;
    
    // remove seed tower from sample
    input_towers_temp.erase(input_towers_temp.begin());    
    for (int tit = 0; tit < (int)cluster_towers_temp.size(); tit++){
      // Now go recursively to the next 4 neighbours and add them to the cluster if they fulfill the conditions
      int iEtaTwr = cluster_towers_temp.at(tit).tower_iEta;
      int iPhiTwr = cluster_towers_temp.at(tit).tower_iPhi;
      for (int ait = 0; ait < (int)input_towers_temp.size(); ait++){
        int iEtaTwrAgg = input_towers_temp.at(ait).tower_iEta;
        int iPhiTwrAgg = input_towers_temp.at(ait).tower_iPhi;
        
        if (!IsForwardCalorimeter(caloEnum)) {
          if (iPhiTwr < 5 && iPhiTwrAgg > _caloTowersPhi[caloEnum]-5){
            iPhiTwrAgg= iPhiTwrAgg-_caloTowersPhi[caloEnum];
          }
          if (iPhiTwr > _caloTowersPhi[caloEnum]-5 && iPhiTwrAgg < 5){
            iPhiTwr= iPhiTwr-_caloTowersPhi[caloEnum];
          }
        }
        int deltaEta = std::abs(iEtaTwrAgg-iEtaTwr);
        int deltaPhi = std::abs(iPhiTwrAgg-iPhiTwr);

        if( (deltaEta+deltaPhi) == 1){
          // only aggregate towers with lower energy than current tower
          if(input_towers_temp.at(ait).tower_E >= (cluster_towers_temp.at(tit).tower_E + aggregation_margin_V3)) continue;
          tempstructC.cluster_E+=input_towers_temp.at(ait).tower_E;
          tempstructC.cluster_NTowers++;
          cluster_towers_temp.push_back(input_towers_temp.at(ait));
          if(!(std::find(clslabels_temp.begin(), clslabels_temp.end(), input_towers_temp.at(ait).tower_trueID) != clslabels_temp.end())){
            tempstructC.cluster_NtrueID++;
            clslabels_temp.push_back(input_towers_temp.at(ait).tower_trueID);
          }
          input_towers_temp.erase(input_towers_temp.begin()+ait);
          ait--;
        }
      }
    }
  } 
  return tempstructC;
}

//**************************************************************************************************************
//**************************************************************************************************************
// find clusters with common edges or corners, separate if energy increases in neighboring cell
//**************************************************************************************************************
//**************************************************************************************************************
clustersStrct findMACluster(
                                float seed,                                     // minimum seed energy
                              int caloEnum,                                   // calorimeter type
                              std::vector<towersStrct> &input_towers_temp,    // temporary full tower array
                              std::vector<towersStrct> &cluster_towers_temp,  // towers associated to cluster
                              std::vector<int> clslabels_temp                 // MC labels in cluster
                            ){
  clustersStrct tempstructC;
  
  if(input_towers_temp.at(0).tower_E > seed){
    if (verbosityCLS > 2) std::cout << "new cluster" << std::endl;
    // fill seed cell information into current cluster
    tempstructC.cluster_E       = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_seed    = input_towers_temp.at(0).tower_E;
    tempstructC.cluster_NTowers = 1;
    tempstructC.cluster_NtrueID = 1;
    tempstructC.cluster_trueID = input_towers_temp.at(0).tower_trueID; // TODO save all MC labels?
    cluster_towers_temp.push_back(input_towers_temp.at(0));
    clslabels_temp.push_back(input_towers_temp.at(0).tower_trueID);
    if (verbosityCLS > 2) std::cout << "seed: "<<  input_towers_temp.at(0).tower_iEta << "\t" << input_towers_temp.at(0).tower_iPhi << "\t" << input_towers_temp.at(0).tower_iL << "\t E:"<< tempstructC.cluster_E << std::endl;
    
    
    // remove seed tower from sample
    input_towers_temp.erase(input_towers_temp.begin());
    for (int tit = 0; tit < (int)cluster_towers_temp.size(); tit++){
      // Now go recursively to all neighbours and add them to the cluster if they fulfill the conditions
      int iEtaTwr = cluster_towers_temp.at(tit).tower_iEta;
      int iPhiTwr = cluster_towers_temp.at(tit).tower_iPhi;
      int iLTwr   = cluster_towers_temp.at(tit).tower_iL;
      int refC = 0;
      for (int ait = 0; ait < (int)input_towers_temp.size(); ait++){
        int iEtaTwrAgg = input_towers_temp.at(ait).tower_iEta;
        int iPhiTwrAgg = input_towers_temp.at(ait).tower_iPhi;
        int iLTwrAgg   = input_towers_temp.at(ait).tower_iL;
        
        if (!IsForwardCalorimeter(caloEnum)) {
          if (iPhiTwr < 5 && iPhiTwrAgg > _caloTowersPhi[caloEnum]-5){
            iPhiTwrAgg= iPhiTwrAgg-_caloTowersPhi[caloEnum];
          }
          if (iPhiTwr > _caloTowersPhi[caloEnum]-5 && iPhiTwrAgg < 5){
            iPhiTwr= iPhiTwr-_caloTowersPhi[caloEnum];
          }
        }
        
        int deltaL    = TMath::Abs(iLTwrAgg-iLTwr) ;
        int deltaPhi  = TMath::Abs(iPhiTwrAgg-iPhiTwr) ;
        int deltaEta  = TMath::Abs(iEtaTwrAgg-iEtaTwr) ;
        bool neighbor = (deltaL+deltaPhi+deltaEta == 1);
        bool corner2D = (deltaL == 0 && deltaPhi == 1 && deltaEta == 1) || (deltaL == 1 && deltaPhi == 0 && deltaEta == 1) || (deltaL == 1 && deltaPhi == 1 && deltaEta == 0);          
//         first condition asks for V3-like neighbors, while second condition also checks diagonally attached towers
        if(neighbor || corner2D ){

          // only aggregate towers with lower energy than current tower
          if(caloEnum != kLFHCAL){
            if(input_towers_temp.at(ait).tower_E >= (cluster_towers_temp.at(tit).tower_E + aggregation_margin_MA[caloEnum])) continue;
          } 
          tempstructC.cluster_E+=input_towers_temp.at(ait).tower_E;
          tempstructC.cluster_NTowers++;
          cluster_towers_temp.push_back(input_towers_temp.at(ait));
          if(!(std::find(clslabels_temp.begin(), clslabels_temp.end(), input_towers_temp.at(ait).tower_trueID) != clslabels_temp.end())){
            tempstructC.cluster_NtrueID++;
            clslabels_temp.push_back(input_towers_temp.at(ait).tower_trueID);
          }
          if (verbosityCLS > 2) std::cout << "aggregated: "<< iEtaTwrAgg << "\t" << iPhiTwrAgg << "\t" << iLTwrAgg << "\t E:" << input_towers_temp.at(ait).tower_E << "\t reference: "<< refC << "\t"<< iEtaTwr << "\t" << iPhiTwr << "\t" << iLTwr << "\t cond.: \t"<< neighbor << "\t" << corner2D << "\t  diffs: " << deltaEta << "\t" << deltaPhi << "\t" << deltaL<< std::endl;

          input_towers_temp.erase(input_towers_temp.begin()+ait);
          ait--;
          refC++;
        }
      }
    }
  } 
  return tempstructC;
}

//**************************************************************************************************************
//**************************************************************************************************************
// ANCHOR main function to be called in event loop
// - run clusterizer and create new cluster list
//**************************************************************************************************************
//**************************************************************************************************************
void runclusterizer(
  int clusterizerEnum,                // which clusterizer are you running
  int caloEnum,                       // which calo are you evaluating
  float seedE,                        // what is the minimum energy of the leading cell in the cluster
  float aggE,                         // what is the minimum energy of a cell in the clusters
  unsigned short primaryTrackSource   // which track source are you matching the cluster to
){
  
  if (caloEnum == kLFHCAL && clusterizerEnum > 0) return;
  int nclusters = 0;
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if (icalo == kLFHCAL && ialgo > 0) continue;
      if(!h_clusterizer_nonagg_towers[icalo][ialgo])h_clusterizer_nonagg_towers[icalo][ialgo] 	= new TH2F(Form("h_clusterizer_nonagg_towers%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 50,0,50,120,0,0.6);

    }
  }

  // vector of towers used in the clusterizer
  std::vector<towersStrct> input_towers = readTowersForCalo(caloEnum, aggE);
  // fill vector with towers for clusterization above aggregation threshold
  if (caloEnum == kFHCAL){
    input_towers = readTowersForCalo(caloEnum, aggE, 2.);
  } else {
    input_towers = readTowersForCalo(caloEnum, aggE, 1.);
  }
  // vector of towers within the currently found cluster
  std::vector<towersStrct> cluster_towers;

  // sort vector in descending energy order
  std::sort(input_towers.begin(), input_towers.end(), &acompare);
  std::vector<int> clslabels;
  if (verbosityCLS > 3) std::cout<< "running :" << str_clusterizer[clusterizerEnum].Data() << std::endl;
  while (!input_towers.empty()) {
    cluster_towers.clear();
    clslabels.clear();
    clustersStrct tempstructC;
    // always start with highest energetic tower
    if(input_towers.at(0).tower_E > seedE){
      if (verbosityCLS > 3) std::cout<< "seed: " << input_towers.at(0).tower_E << "\t" << input_towers.at(0).tower_iEta <<  "\t" << input_towers.at(0).tower_iPhi<< std::endl;
      if(clusterizerEnum==kC3){
        tempstructC = findCircularCluster(3, seedE, caloEnum, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==kC5){
        tempstructC = findCircularCluster(5, seedE, caloEnum, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==k3x3){
        tempstructC = findSquareCluster(3, seedE, caloEnum, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==k5x5){
        tempstructC = findSquareCluster(5, seedE, caloEnum, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==kV1){  
        tempstructC = findV1Cluster(seedE, caloEnum, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==kV3){  
        tempstructC = findV3Cluster(seedE, caloEnum, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==kMA){  
        tempstructC = findMACluster(seedE, caloEnum, input_towers, cluster_towers, clslabels);
      } else {
        std::cout << "incorrect clusterizer selected! Enum " << clusterizerEnum << " not defined!" << std::endl;
        return;
      }

      // determine remaining cluster properties from its towers
      float* showershape_eta_phi = CalculateM02andWeightedPosition(cluster_towers, tempstructC.cluster_E, caloEnum, false);
      tempstructC.cluster_M02 = showershape_eta_phi[0];
      tempstructC.cluster_M20 = showershape_eta_phi[1];
      tempstructC.cluster_Eta = showershape_eta_phi[2];
      tempstructC.cluster_Phi = showershape_eta_phi[3];
      tempstructC.cluster_X = showershape_eta_phi[4];
      tempstructC.cluster_Y = showershape_eta_phi[5];
      tempstructC.cluster_Z = showershape_eta_phi[6];
      if (tracksEnabled){
        tempstructC.cluster_matchedTracks = isClusterMatched(tempstructC, caloEnum, clusterizerEnum, primaryTrackSource, true);
        if(tempstructC.cluster_matchedTracks.size() > 0) {
                tempstructC.cluster_isMatched = true;
        } else {
          tempstructC.cluster_isMatched = false;
        }
      } else {
        tempstructC.cluster_isMatched = false;
      }
      if(verbosityCLS>1) std::cout << clusterizerEnum << "\t" << nclusters << "\tcluster with E = " << tempstructC.cluster_E << "\tEta: " << tempstructC.cluster_Eta<< "\tPhi: " << tempstructC.cluster_Phi
                              << "\tX: " << tempstructC.cluster_X<< "\tY: " << tempstructC.cluster_Y<< "\tZ: " << tempstructC.cluster_Z<< "\tntowers: " << tempstructC.cluster_NTowers 
                              << "\ttrueID: " << tempstructC.cluster_trueID << std::endl;

      // apply calibration if desired
      if(_doClusterECalibration){
          tempstructC.cluster_E/=getCalibrationValue(tempstructC.cluster_E, caloEnum, clusterizerEnum, tempstructC.cluster_trueID);
          tempstructC.cluster_E*=getEnergySmearing( caloEnum, clusterizerEnum);
      }
      
      _clusters_calo[clusterizerEnum][caloEnum].push_back(tempstructC);
      
      nclusters++;
    } else {
      if (verbosityCLS > 3) std::cout<< "remaing: "<< (int)input_towers.size() << " largest:" << input_towers.at(0).tower_E << "\t" << input_towers.at(0).tower_iEta <<  "\t" << input_towers.at(0).tower_iPhi<< std::endl;
      for (int ait = 0; ait < (int)input_towers.size(); ait++){
        if (verbosityCLS > 3) std::cout<< input_towers.at(ait).tower_E << "\t" << input_towers.at(ait).tower_iEta <<  "\t" << input_towers.at(ait).tower_iPhi<< std::endl;
        h_clusterizer_nonagg_towers[caloEnum][clusterizerEnum]->Fill(input_towers.size(),input_towers.at(ait).tower_E);
      }
      input_towers.clear();
      
    }
  }
  
  if (verbosityCLS > 3) std::cout<< "finished this event for " << str_clusterizer[clusterizerEnum].Data() << std::endl;
  
  std::sort(_clusters_calo[clusterizerEnum][caloEnum].begin(), _clusters_calo[clusterizerEnum][caloEnum].end(), &acompareCl);    

  return; 
}

//**************************************************************************************************************
//**************************************************************************************************************
// store diagnostics histograms for clusterizer
//**************************************************************************************************************
//**************************************************************************************************************
void clusterizerSave(){

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_CLSIZER.root",outputDir.Data()),"RECREATE");

  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if (h_clusterizer_nonagg_towers[icalo][ialgo])    h_clusterizer_nonagg_towers[icalo][ialgo]->Write();
      if (h_clusterizer_all_2D_delta[icalo][ialgo])     h_clusterizer_all_2D_delta[icalo][ialgo]->Write();
      if (h_clusterizer_matched_2D_delta[icalo][ialgo]) h_clusterizer_matched_2D_delta[icalo][ialgo]->Write();
      if (h_clusterizer_isMatched_E[icalo][ialgo])      h_clusterizer_isMatched_E[icalo][ialgo]->Write();
      for (int jcalo=0; jcalo < _active_calo; jcalo++){
        if (h_clusterizer_all_2D_deltaCalo[icalo][jcalo][ialgo])            h_clusterizer_all_2D_deltaCalo[icalo][jcalo][ialgo]->Write();
        if (h_clusterizer_calomatch_2D_deltaCalo[icalo][jcalo][ialgo])      h_clusterizer_calomatch_2D_deltaCalo[icalo][jcalo][ialgo]->Write();
        if (h_clusterizer_calomatchtrue_2D_deltaCalo[icalo][jcalo][ialgo])  h_clusterizer_calomatchtrue_2D_deltaCalo[icalo][jcalo][ialgo]->Write();
        if (h_clusterizer_isMatchedCalo_E[icalo][jcalo][ialgo])             h_clusterizer_isMatchedCalo_E[icalo][jcalo][ialgo]->Write();
      }
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}

//**************************************************************************************************************
//**************************************************************************************************************
// clear all cluster vectors properly after each event
//**************************************************************************************************************
//**************************************************************************************************************
void clearClusterVectors(){
  for (int icalo = 0; icalo < maxcalo; icalo++){
    for (int ialgo = 0; ialgo < maxAlgo; ialgo++){
      if (!_clusters_calo[ialgo][icalo].empty() && caloEnabled[icalo] ){
        if(verbosityCLS > 3) std::cout << "contained " <<  _clusters_calo[ialgo][icalo].size() << " for calo " << icalo << " \t clusterizer " << ialgo << std::endl ;
        if(verbosityCLS > 3) std::cout << "clearing ...." << std::endl ;
        _clusters_calo[ialgo][icalo].clear();
      }
    }
  }
  
}
