
// ANCHOR debug output verbosity
Int_t verbosityHITS =0;

// ANCHOR create histograms globally
bool histsCreatedHITS = false;
TH2F*  h_hits_layer_xy[_maxProjectionLayers]                        = {NULL};
TH2F*  h_hits_layer_etaphi[_maxProjectionLayers]                    = {NULL};
TH2F*  h_trackProj_layer_xy[_maxProjectionLayers]                   = {NULL};
TH2F*  h_trackProj_layer_etaphi[_maxProjectionLayers]               = {NULL};
TH2F*  h_hitslayer_vs_tracks[_maxProjectionLayers]                  = {NULL};
TH2F*  h_hitslayer_vs_clusters[_maxProjectionLayers][_active_calo]  = {{NULL}};
TH2F*  h_hitslayer_vs_towers[_maxProjectionLayers][_active_calo]    = {{NULL}};
TH2F*  h_fwdtracks_vs_clusters[_active_calo]                        = {NULL};

// *****************************************************************************************************
// *****************************************************************************************************
// ANCHOR main function to be called in event loop
// *****************************************************************************************************
// *****************************************************************************************************
void hitstudies(unsigned short primaryTrackSource){
  
  // *****************************************************************************************************
  // ********************** Create histograms ************************************************************
  // *****************************************************************************************************
  if (! histsCreatedHITS){
    for(int il=0;il<_maxProjectionLayers;il++){
      if (!enableLayerHist[il]) continue;
      if(!h_hits_layer_etaphi[il])h_hits_layer_etaphi[il]             = new TH2F(Form("h_hits_layer_etaphi_%s", GetProjectionNameFromIndex(layerIndexHist[il]).Data()), "", 450, -4.5 , 4.5,300, -TMath::Pi(),TMath::Pi());
      if(!h_trackProj_layer_etaphi[il])h_trackProj_layer_etaphi[il]   = new TH2F(Form("h_trackProj_layer_etaphi_%s", GetProjectionNameFromIndex(layerIndexHist[il]).Data()), "", 450, -4.5 , 4.5,300, -TMath::Pi(),TMath::Pi());
      if(!h_hitslayer_vs_tracks[il])h_hitslayer_vs_tracks[il]         = new TH2F(Form("h_hitslayer_vs_tracks_%s", GetProjectionNameFromIndex(layerIndexHist[il]).Data()), "", 100, 0 , 100,100, 0 , 100);
      if (GetRegionFromIndex(layerIndexHist[il]) != 1){
        if(!h_hits_layer_xy[il])h_hits_layer_xy[il]                   = new TH2F(Form("h_hits_layer_xy_%s", GetProjectionNameFromIndex(layerIndexHist[il]).Data()), "", 300, -200,200,300, -200,200);
        if(!h_trackProj_layer_xy[il])h_trackProj_layer_xy[il]         = new TH2F(Form("h_trackProj_layer_xy_%s", GetProjectionNameFromIndex(layerIndexHist[il]).Data()), "", 300, -200,200,300, -200,200);
      }
      if (IsCaloProjection(layerIndexHist[il]) || HasTimingLayer(layerIndexHist[il]) ){
        for(int icalo=0;icalo<_active_calo;icalo++){
          if (!caloEnabled[icalo]) continue;
          if(!h_hitslayer_vs_clusters[il][icalo])h_hitslayer_vs_clusters[il][icalo]   = new TH2F(Form("h_hitslayer_vs_clusters_%s_%s", GetProjectionNameFromIndex(layerIndexHist[il]).Data(),str_calorimeter[icalo].Data()), "", 100, 0 , 100,30, 0 , 30);
          if(!h_hitslayer_vs_towers[il][icalo])h_hitslayer_vs_towers[il][icalo]       = new TH2F(Form("h_hitslayer_vs_towers_%s_%s", GetProjectionNameFromIndex(layerIndexHist[il]).Data(),str_calorimeter[icalo].Data()), "", 100, 0 , 100,300, 0 , 300);
          if(!h_fwdtracks_vs_clusters[icalo])h_fwdtracks_vs_clusters[icalo]           = new TH2F(Form("h_fwdtracks_vs_clusters_%s",str_calorimeter[icalo].Data()), "", 100, 0 , 100,30, 0 , 30);
        }
      }
    }
    histsCreatedHITS = true;
  }

  // *****************************************************************************************************
  // ANCHOR Hits, fill eta-phi and x-y distribitions for hits
  // *****************************************************************************************************
  int ihitslayer[_maxProjectionLayers] = {0};
  for(Int_t ihit=0; ihit<_nHitsLayers; ihit++){
    if(verbosityHITS>1) 
      std::cout << "\tHIT: hit " << ihit<<"/"<< _nHitsLayers << "\tin layer " << _hits_layerID[ihit] << "\twith X = " << _hits_x[ihit] << " cm" << std::endl;

    TVector3 hitvec(_hits_x[ihit],_hits_y[ihit],_hits_z[ihit]);
    float hiteta = hitvec.Eta();
    float hitphi = (hitvec.Phi()<0 ? hitvec.Phi()+TMath::Pi() : hitvec.Phi()-TMath::Pi());

    Int_t histIndex = ReturnIndexForwardLayer(_hits_layerID[ihit]);
    if (histIndex == -1) continue;
    ihitslayer[histIndex]++;
    if(verbosityHITS>1) std::cout << "layer index: \t" << histIndex << "\t" << _hits_layerID[ihit] << std::endl;
    if( histIndex<_maxProjectionLayers){
      if (GetRegionFromIndex(_hits_layerID[ihit]) != 1) 
        if (h_hits_layer_xy[histIndex])h_hits_layer_xy[histIndex]->Fill(_hits_x[ihit],_hits_y[ihit]);
      if (h_hits_layer_etaphi[histIndex])h_hits_layer_etaphi[histIndex]->Fill(hiteta,hitphi);
    }
  }
  // *****************************************************************************************************
  // count tracks in different eta regions
  // *****************************************************************************************************
  int itracks[3] = {0, 0, 0};
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    // Select track source
    if (_track_source[itrk] != primaryTrackSource) { continue; }

    TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    float trketa = trackvec.Eta();
    if(trketa>1.2 && trketa<4)            itracks[2]++;
    else if(trketa>-4 && trketa<-1.7)     itracks[0]++;
    else if(trketa>=-1.7 && trketa<=1.2)  itracks[1]++;
  }
  // *****************************************************************************************************
  // correlate hits in different layers with tracks/cluster 
  // *****************************************************************************************************
  for(int ik = 0; ik < (int)_maxProjectionLayers; ik++){
    if (!enableLayerHist[ik]) continue;
    h_hitslayer_vs_tracks[ik]->Fill(ihitslayer[ik], itracks[GetRegionFromIndex(ik)]);
    if (IsCaloProjection(layerIndexHist[ik]) || HasTimingLayer(layerIndexHist[ik]) ){
      for(int icalo=0;icalo<_active_calo;icalo++){
        if (!caloEnabled[icalo]) continue;
        h_hitslayer_vs_clusters[ik][icalo]->Fill(ihitslayer[ik], _clusters_calo[kMA][icalo].size() );
        h_hitslayer_vs_towers[ik][icalo]->Fill(ihitslayer[ik], ReturnMaxTowerCalo(icalo));
      }
    }
  }
  // *****************************************************************************************************
  // ANCHOR Track projections
  // *****************************************************************************************************
  for(Int_t iproj=0; iproj<_nProjections; iproj++){
//     if(_track_Proj_x[iproj]==0 && _track_Proj_y[iproj]==0) continue;
    if(TMath::Abs(_track_Proj_t[iproj])< 2e-20) continue;
    if(verbosityHITS>1) 
      std::cout << "\tProjection: proj " << iproj << "\tin layer " << _track_ProjLayer[iproj] << "\t t = " << _track_Proj_t[iproj] << " ps" << "\t X = " << _track_Proj_x[iproj] << " cm" << "\t Y = " << _track_Proj_y[iproj] << " cm" << "\t Z = " << _track_Proj_z[iproj] << " cm"<< std::endl;
    TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
    float projeta = projvec.Eta();
    float projphi = (projvec.Phi()<0 ? projvec.Phi()+TMath::Pi() : projvec.Phi()-TMath::Pi());

    Int_t histIndex = ReturnIndexForwardLayer(_track_ProjLayer[iproj]);
    if (histIndex == -1) continue;
    if( histIndex<_maxProjectionLayers){
      if (GetRegionFromIndex(_track_ProjLayer[iproj]) != 1)       
        if (h_trackProj_layer_xy[histIndex]) h_trackProj_layer_xy[histIndex]->Fill(_track_Proj_x[iproj],_track_Proj_y[iproj]);
      if (h_trackProj_layer_etaphi[histIndex]) h_trackProj_layer_etaphi[histIndex]->Fill(projeta,projphi);
    }
  }

}
// *****************************************************************************************************
// *****************************************************************************************************
// ANCHOR save function after event loop
// *****************************************************************************************************
// *****************************************************************************************************
void hitstudiesSave(){
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_HITS.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int il=0;il<_maxProjectionLayers;il++){
    if (!enableLayerHist[il]) continue;    
    fileOutput->mkdir(GetProjectionNameFromIndex(layerIndexHist[il]).Data());
    fileOutput->cd(GetProjectionNameFromIndex(layerIndexHist[il]).Data());
    if (h_hits_layer_xy[il]) h_hits_layer_xy[il]->Scale(1/_nEventsTree);
    if (h_hits_layer_xy[il]) h_hits_layer_xy[il]->Write();
    if (h_hits_layer_etaphi[il]) h_hits_layer_etaphi[il]->Scale(1/_nEventsTree);
    if (h_hits_layer_etaphi[il])h_hits_layer_etaphi[il]->Write();

    if (h_trackProj_layer_xy[il]) h_trackProj_layer_xy[il]->Scale(1/_nEventsTree);
    if (h_trackProj_layer_xy[il]) h_trackProj_layer_xy[il]->Write();
    if (h_trackProj_layer_etaphi[il]) h_trackProj_layer_etaphi[il]->Scale(1/_nEventsTree);
    if (h_trackProj_layer_etaphi[il]) h_trackProj_layer_etaphi[il]->Write();

    if (h_hitslayer_vs_tracks[il]) h_hitslayer_vs_tracks[il]->Scale(1./h_hitslayer_vs_tracks[il]->GetEntries());
    if (h_hitslayer_vs_tracks[il]) h_hitslayer_vs_tracks[il]->Write();
    for(int icalo=0;icalo<_active_calo;icalo++){
      if (h_hitslayer_vs_clusters[il][icalo]) h_hitslayer_vs_clusters[il][icalo]->Scale(1./h_hitslayer_vs_clusters[il][icalo]->GetEntries());
      if (h_hitslayer_vs_clusters[il][icalo]) h_hitslayer_vs_clusters[il][icalo]->Write();
      if (h_hitslayer_vs_towers[il][icalo]) h_hitslayer_vs_towers[il][icalo]->Scale(1./h_hitslayer_vs_towers[il][icalo]->GetEntries());
      if (h_hitslayer_vs_towers[il][icalo]) h_hitslayer_vs_towers[il][icalo]->Write();
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
