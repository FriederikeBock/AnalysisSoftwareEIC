
// ANCHOR debug output verbosity
Int_t verbosityHITS = 0;

// ANCHOR create histograms globally
// TString str_TRKEFF_mcparticles[5] = {"electron", "cpion", "proton", "ckaon", "muon"};
// int int_TRKEFF_mcparticles_PDG[5] = {11 /*e*/,211  /*pi*/,  2212/*p*/,  321/*K*/,  13/*mu*/};

TH2F*  h_hits_FTTL1_xy 	= new TH2F("h_hits_FTTL1_xy", "", 300, -200,200,300, -200,200);
TH2F*  h_hits_FTTL1_etaphi 	= new TH2F("h_hits_FTTL1_etaphi", "", 300, 1 , 4,300, -TMath::Pi(),TMath::Pi());
TH2F*  h_hits_FTTL2_xy 	= new TH2F("h_hits_FTTL2_xy", "", 300, -200,200,300, -200,200);
TH2F*  h_hits_FTTL2_etaphi 	= new TH2F("h_hits_FTTL2_etaphi", "", 300, 1 , 4,300, -TMath::Pi(),TMath::Pi());
TH2F*  h_hits_CTTL0_etaphi 	= new TH2F("h_hits_CTTL0_etaphi", "", 300, -1.5 , 1.5,300, -TMath::Pi(),TMath::Pi());
TH2F*  h_hits_CTTL1_etaphi 	= new TH2F("h_hits_CTTL1_etaphi", "", 300, -1.5 , 1.5,300, -TMath::Pi(),TMath::Pi());

TH2F*  h_trackProj_FTTL1_xy 	= new TH2F("h_trackProj_FTTL1_xy", "", 300, -200,200,300, -200,200);
TH2F*  h_trackProj_FTTL1_etaphi 	= new TH2F("h_trackProj_FTTL1_etaphi", "", 300, 1 , 4,300, -TMath::Pi(),TMath::Pi());
TH2F*  h_trackProj_FTTL2_xy 	= new TH2F("h_trackProj_FTTL2_xy", "", 300, -200,200,300, -200,200);
TH2F*  h_trackProj_FTTL2_etaphi 	= new TH2F("h_trackProj_FTTL2_etaphi", "", 300, 1 , 4,300, -TMath::Pi(),TMath::Pi());
TH2F*  h_trackProj_CTTL0_etaphi 	= new TH2F("h_trackProj_CTTL0_etaphi", "", 300, -1.5 , 1.5,300, -TMath::Pi(),TMath::Pi());
TH2F*  h_trackProj_CTTL1_etaphi 	= new TH2F("h_trackProj_CTTL1_etaphi", "", 300, -1.5 , 1.5,300, -TMath::Pi(),TMath::Pi());

// ANCHOR main function to be called in event loop
void hitstudies(){

  // ANCHOR Hits
  for(Int_t ihit=0; ihit<_nHitsLayers; ihit++){
    if(verbosityHITS>1) cout << "\tHIT: hit " << ihit << "\tin layer " << _hits_layerID[ihit] << "\twith X = " << _hits_x[ihit] << " cm" << endl;

    TVector3 hitvec(_hits_x[ihit],_hits_y[ihit],_hits_z[ihit]);
    float hiteta = hitvec.Eta();
    float hitphi = (hitvec.Phi()<0 ? hitvec.Phi()+TMath::Pi() : hitvec.Phi()-TMath::Pi());
    // FTTL1
    if(_hits_layerID[ihit]==1){
      h_hits_FTTL1_xy->Fill(_hits_x[ihit],_hits_y[ihit]);
      h_hits_FTTL1_etaphi->Fill(hiteta,hitphi);
    }
    // FTTL2
    if(_hits_layerID[ihit]==2){
      h_hits_FTTL2_xy->Fill(_hits_x[ihit],_hits_y[ihit]);
      h_hits_FTTL2_etaphi->Fill(hiteta,hitphi);
    }
    // CTTL0
    if(_hits_layerID[ihit]==7){
      h_hits_CTTL0_etaphi->Fill(hiteta,hitphi);
    }
    // CTTL1
    if(_hits_layerID[ihit]==8){
      h_hits_CTTL1_etaphi->Fill(hiteta,hitphi);
    }
  }

  // ANCHOR Track projections
  for(Int_t iproj=0; iproj<_nProjections; iproj++){
    if(verbosityHITS>1) cout << "\tProjection: proj " << iproj << "\tin layer " << _track_ProjLayer[iproj] << "\twith X = " << _track_Proj_x[iproj] << " cm" << endl;
    TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
    float projeta = projvec.Eta();
    float projphi = (projvec.Phi()<0 ? projvec.Phi()+TMath::Pi() : projvec.Phi()-TMath::Pi());
    // FTTL1
    if(_track_ProjLayer[iproj]==1){
      h_trackProj_FTTL1_xy->Fill(_track_Proj_x[iproj],_track_Proj_y[iproj]);
      h_trackProj_FTTL1_etaphi->Fill(projeta,projphi);
    }
    // FTTL2
    if(_track_ProjLayer[iproj]==2){
      h_trackProj_FTTL2_xy->Fill(_track_Proj_x[iproj],_track_Proj_y[iproj]);
      h_trackProj_FTTL2_etaphi->Fill(projeta,projphi);
    }
    // CTTL0
    if(_track_ProjLayer[iproj]==7){
      h_trackProj_CTTL0_etaphi->Fill(projeta,projphi);
    }
    // CTTL1
    if(_track_ProjLayer[iproj]==8){
      h_trackProj_CTTL1_etaphi->Fill(projeta,projphi);
    }
  }

}

// ANCHOR save function after event loop
void hitstudiesSave(){
  // make output directory
    gSystem->Exec("mkdir -p "+outputDir + "/hitstudies");

  // define output file
  TFile* fileOutput = new TFile(Form("%s/hitstudies/output_HITS.root",outputDir.Data()),"RECREATE");

  // write histograms
  if(h_hits_FTTL1_xy){
    h_hits_FTTL1_xy->Scale(1/_nEventsTree);
    h_hits_FTTL1_xy->Write();
  }
  if(h_hits_FTTL1_etaphi){
    h_hits_FTTL1_etaphi->Scale(1/_nEventsTree);
    h_hits_FTTL1_etaphi->Write();
  }
  if(h_hits_FTTL2_xy){
    h_hits_FTTL2_xy->Scale(1/_nEventsTree);
    h_hits_FTTL2_xy->Write();
  }
  if(h_hits_FTTL2_etaphi){
    h_hits_FTTL2_etaphi->Scale(1/_nEventsTree);
    h_hits_FTTL2_etaphi->Write();
  }
  if(h_hits_CTTL0_etaphi){
    h_hits_CTTL0_etaphi->Scale(1/_nEventsTree);
    h_hits_CTTL0_etaphi->Write();
  }
  if(h_hits_CTTL1_etaphi){
    h_hits_CTTL1_etaphi->Scale(1/_nEventsTree);
    h_hits_CTTL1_etaphi->Write();
  }

  if(h_trackProj_FTTL1_xy){
    h_trackProj_FTTL1_xy->Scale(1/_nEventsTree);
    h_trackProj_FTTL1_xy->Write();
  }
  if(h_trackProj_FTTL1_etaphi){
    h_trackProj_FTTL1_etaphi->Scale(1/_nEventsTree);
    h_trackProj_FTTL1_etaphi->Write();
  }
  if(h_trackProj_FTTL2_xy){
    h_trackProj_FTTL2_xy->Scale(1/_nEventsTree);
    h_trackProj_FTTL2_xy->Write();
  }
  if(h_trackProj_FTTL2_etaphi){
    h_trackProj_FTTL2_etaphi->Scale(1/_nEventsTree);
    h_trackProj_FTTL2_etaphi->Write();
  }
  if(h_trackProj_CTTL0_etaphi){
    h_trackProj_CTTL0_etaphi->Scale(1/_nEventsTree);
    h_trackProj_CTTL0_etaphi->Write();
  }
  if(h_trackProj_CTTL1_etaphi){
    h_trackProj_CTTL1_etaphi->Scale(1/_nEventsTree);
    h_trackProj_CTTL1_etaphi->Write();
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}