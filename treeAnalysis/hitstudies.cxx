
// ANCHOR debug output verbosity
Int_t verbosityHITS = 0;

// ANCHOR create histograms globally
// TString str_TRKEFF_mcparticles[5] = {"electron", "cpion", "proton", "ckaon", "muon"};
// int int_TRKEFF_mcparticles_PDG[5] = {11 /*e*/,211  /*pi*/,  2212/*p*/,  321/*K*/,  13/*mu*/};

TH2F*  h_hits_layer_xy[_maxProjectionLayers] 	= {NULL};
TH2F*  h_hits_layer_etaphi[_maxProjectionLayers] 	= {NULL};

TH2F*  h_trackProj_layer_xy[_maxProjectionLayers] 	= {NULL};
TH2F*  h_trackProj_layer_etaphi[_maxProjectionLayers] 	= {NULL};

TH2F*  h_hitslayer_vs_tracks[_maxProjectionLayers] 	= {NULL};
TH2F*  h_hitslayer_vs_clusters[_maxProjectionLayers][_active_calo] 	= {{NULL}};
TH2F*  h_hitslayer_vs_towers[_maxProjectionLayers][_active_calo] 	= {{NULL}};
TH2F*  h_fwdtracks_vs_clusters[_active_calo] 	= {NULL};


// ANCHOR main function to be called in event loop
void hitstudies(){
  for(int il=0;il<_maxProjectionLayers;il++){
    if(!h_hits_layer_xy[il])h_hits_layer_xy[il] = new TH2F(Form("h_hits_layer_xy_%s", GetProjectionNameFromIndex(il).Data()), "", 300, -200,200,300, -200,200);
    if(!h_hits_layer_etaphi[il])h_hits_layer_etaphi[il]  	= new TH2F(Form("h_hits_layer_etaphi_%s", GetProjectionNameFromIndex(il).Data()), "", 300, 1 , 4,300, -TMath::Pi(),TMath::Pi());
    if(!h_trackProj_layer_xy[il])h_trackProj_layer_xy[il] 	= new TH2F(Form("h_trackProj_layer_xy_%s", GetProjectionNameFromIndex(il).Data()), "", 300, -200,200,300, -200,200);
    if(!h_trackProj_layer_etaphi[il])h_trackProj_layer_etaphi[il] 	= new TH2F(Form("h_trackProj_layer_etaphi_%s", GetProjectionNameFromIndex(il).Data()), "", 300, 1 , 4,300, -TMath::Pi(),TMath::Pi());

    if(!h_hitslayer_vs_tracks[il])h_hitslayer_vs_tracks[il] 	= new TH2F(Form("h_hitslayer_vs_tracks_%s", GetProjectionNameFromIndex(il).Data()), "", 100, 0 , 100,100, 0 , 100);
    for(int icalo=0;icalo<_active_calo;icalo++){
      if(!h_hitslayer_vs_clusters[il][icalo])h_hitslayer_vs_clusters[il][icalo] 	= new TH2F(Form("h_hitslayer_vs_clusters_%s_%s", GetProjectionNameFromIndex(il).Data(),str_calorimeter[icalo].Data()), "", 100, 0 , 100,30, 0 , 30);
      if(!h_hitslayer_vs_towers[il][icalo])h_hitslayer_vs_towers[il][icalo] 	= new TH2F(Form("h_hitslayer_vs_towers_%s_%s", GetProjectionNameFromIndex(il).Data(),str_calorimeter[icalo].Data()), "", 100, 0 , 100,300, 0 , 300);
      if(!h_fwdtracks_vs_clusters[icalo])h_fwdtracks_vs_clusters[icalo] 	= new TH2F(Form("h_fwdtracks_vs_clusters_%s",str_calorimeter[icalo].Data()), "", 100, 0 , 100,30, 0 , 30);

    }

  }
  // ANCHOR Hits
  int ihitslayer[10] = {0};
  for(Int_t ihit=0; ihit<_nHitsLayers; ihit++){
    if(verbosityHITS>1) cout << "\tHIT: hit " << ihit << "\tin layer " << _hits_layerID[ihit] << "\twith X = " << _hits_x[ihit] << " cm" << endl;

    TVector3 hitvec(_hits_x[ihit],_hits_y[ihit],_hits_z[ihit]);
    float hiteta = hitvec.Eta();
    float hitphi = (hitvec.Phi()<0 ? hitvec.Phi()+TMath::Pi() : hitvec.Phi()-TMath::Pi());

    ihitslayer[_hits_layerID[ihit]]++;

    if(_hits_layerID[ihit]<_maxProjectionLayers){
      h_hits_layer_xy[_hits_layerID[ihit]]->Fill(_hits_x[ihit],_hits_y[ihit]);
      h_hits_layer_etaphi[_hits_layerID[ihit]]->Fill(hiteta,hitphi);
    }
  }
  int itracksfwd = 0;
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    float trketa = trackvec.Eta();
    if(trketa>1.05)itracksfwd++;
  }
  for(int il=0;il<_maxProjectionLayers;il++){
    h_hitslayer_vs_tracks[il]->Fill(ihitslayer[il],itracksfwd);
    for(int icalo=0;icalo<_active_calo;icalo++){
      h_hitslayer_vs_clusters[il][icalo]->Fill(ihitslayer[il],icalo==0 ? _nclusters_V3_FHCAL : _nclusters_V3_FEMC);
      h_hitslayer_vs_towers[il][icalo]->Fill(ihitslayer[il],icalo==0 ? _nTowers_FHCAL : _nTowers_FEMC);
    }
  }

  // ANCHOR Track projections
  for(Int_t iproj=0; iproj<_nProjections; iproj++){
    if(verbosityHITS>1) cout << "\tProjection: proj " << iproj << "\tin layer " << _track_ProjLayer[iproj] << "\twith X = " << _track_Proj_x[iproj] << " cm" << endl;
    if(_track_Proj_x[iproj]==0 && _track_Proj_y[iproj]==0) continue;
    if(_track_Proj_t[iproj]<-9000) continue;
    TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
    float projeta = projvec.Eta();
    float projphi = (projvec.Phi()<0 ? projvec.Phi()+TMath::Pi() : projvec.Phi()-TMath::Pi());

    if(_track_ProjLayer[iproj]<_maxProjectionLayers){
      h_trackProj_layer_xy[_track_ProjLayer[iproj]]->Fill(_track_Proj_x[iproj],_track_Proj_y[iproj]);
      h_trackProj_layer_etaphi[_track_ProjLayer[iproj]]->Fill(projeta,projphi);
    }
  }

}

// ANCHOR save function after event loop
void hitstudiesSave(){
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_HITS.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int il=0;il<_maxProjectionLayers;il++){

    h_hits_layer_xy[il]->Scale(1/_nEventsTree);
    h_hits_layer_xy[il]->Write();
    h_hits_layer_etaphi[il]->Scale(1/_nEventsTree);
    h_hits_layer_etaphi[il]->Write();

    h_trackProj_layer_xy[il]->Scale(1/_nEventsTree);
    h_trackProj_layer_xy[il]->Write();
    h_trackProj_layer_etaphi[il]->Scale(1/_nEventsTree);
    h_trackProj_layer_etaphi[il]->Write();

    h_hitslayer_vs_tracks[il]->Scale(1./h_hitslayer_vs_tracks[il]->GetEntries());
    h_hitslayer_vs_tracks[il]->Write();
    for(int icalo=0;icalo<_active_calo;icalo++){
      h_hitslayer_vs_clusters[il][icalo]->Scale(1./h_hitslayer_vs_clusters[il][icalo]->GetEntries());
      h_hitslayer_vs_clusters[il][icalo]->Write();
      h_hitslayer_vs_towers[il][icalo]->Scale(1./h_hitslayer_vs_towers[il][icalo]->GetEntries());
      h_hitslayer_vs_towers[il][icalo]->Write();
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
