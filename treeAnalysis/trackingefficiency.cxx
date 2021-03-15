
// ANCHOR debug output verbosity
Int_t verbosityTRKEFF = 0;

// ANCHOR create histograms globally
const int nPart_TRKEFF = 5;
TString str_TRKEFF_mcparticles[nPart_TRKEFF+1] = {"electron", "cpion", "proton", "ckaon", "muon", "all"};
int int_TRKEFF_mcparticles_PDG[nPart_TRKEFF] = {11 /*e*/,211  /*pi*/,  2212/*p*/,  321/*K*/,  13/*mu*/};

// **********************************************************************************************
// ********************************** Set up histos efficiency histos ***************************
// **********************************************************************************************
TH2F*  h_chargedpart_MC_pT        = new TH2F("h_chargedpart_MC_pT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_particle_MC_pT[nPart_TRKEFF];
TH2F*  h_chargedpart_rec_pT       = new TH2F("h_chargedpart_rec_pT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_chargedpart_rec_truepT   = new TH2F("h_chargedpart_rec_truepT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_particle_rec_pT[nPart_TRKEFF];
TH2F*  h_particle_rec_truepT[nPart_TRKEFF];

TH2F*  h_chargedpart_MC_p         = new TH2F("h_chargedpart_MC_p", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_particle_MC_p[nPart_TRKEFF];
TH2F*  h_chargedpart_rec_p        = new TH2F("h_chargedpart_rec_p", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_chargedpart_rec_truep    = new TH2F("h_chargedpart_rec_truep", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_chargedpart_rec_trueE    = new TH2F("h_chargedpart_rec_trueE", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_particle_rec_p[nPart_TRKEFF];
TH2F*  h_particle_rec_truep[nPart_TRKEFF];
TH2F*  h_particle_rec_trueE[nPart_TRKEFF];

TH2F*  h_chargedpart_MC_E         = new TH2F("h_chargedpart_MC_E", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_particle_MC_E[nPart_TRKEFF];

// **********************************************************************************************
// ********************************** Set up histos resolution histos ***************************
// **********************************************************************************************
TH1D* hNEvents                                = new TH1D("nEvents","",1,0.5,1.5);
TH1D* hNTracks                                = new TH1D("nTracks","",200,-0.5,199.5);
TH2F* h_tracks_reso_pT[nResoSt][nEta+1][nPart_TRKEFF+1] = {{{NULL}}};
TH2F* h_tracks_reso_p[nResoSt][nEta+1][nPart_TRKEFF+1]  = {{NULL}};
TH2F* h_tracks_resoEta_pT[nResoSt][nEta+1]              = {{NULL}};
TH2F* h_tracks_resoPhi_pT[nResoSt][nEta+1]              = {{NULL}};
TH2F* h_tracksTrue_Eta_pT[nPart_TRKEFF+1][nCuts]        = {{NULL}};
TH2F* h_tracksRec_Eta_pT[nPart_TRKEFF+1][nCuts]         = {{NULL}};
TH2F* h_tracksTrue_Eta_p[nPart_TRKEFF+1][nCuts]         = {{NULL}};
TH2F* h_tracksRec_Eta_p[nPart_TRKEFF+1][nCuts]          = {{NULL}};
Bool_t initResoHist                                     = kFALSE;

// **********************************************************************************************
// **************************** tracking efficiency processing **********************************
// **********************************************************************************************
void trackingefficiency(){
  // cout << "new evt" << endl;
  for(int imc=0; imc<_nMCPart; imc++){
    TVector3 mcpartvec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
    float trueeta = mcpartvec.Eta();
    float truept  = mcpartvec.Pt();
    float truep   = mcpartvec.Mag();

    // cout << "\tMCID " << imc << "\tPDG: " << _mcpart_PDG[imc] << "\tpT: " << TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)) << "\tEta " << trueeta << endl;
    for(int ipart=0;ipart<nPart_TRKEFF;ipart++){
      if(!h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart] 	= new TH2F(Form("h_%s_MC_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);
      if(!h_particle_rec_pT[ipart]) h_particle_rec_pT[ipart] 	= new TH2F(Form("h_%s_rec_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);
      if(!h_particle_rec_truepT[ipart]) h_particle_rec_truepT[ipart] 	= new TH2F(Form("h_%s_rec_truepT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);

      if(!h_particle_MC_p[ipart]) h_particle_MC_p[ipart] 	= new TH2F(Form("h_%s_MC_p",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      if(!h_particle_rec_p[ipart]) h_particle_rec_p[ipart] 	= new TH2F(Form("h_%s_rec_p",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      if(!h_particle_rec_truep[ipart]) h_particle_rec_truep[ipart] 	= new TH2F(Form("h_%s_rec_truep",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      if(!h_particle_rec_trueE[ipart]) h_particle_rec_trueE[ipart] 	= new TH2F(Form("h_%s_rec_trueE",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      
      
      if(!h_particle_MC_E[ipart]) h_particle_MC_E[ipart] 	= new TH2F(Form("h_%s_MC_E",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      
      if(abs(_mcpart_PDG[imc]) == int_TRKEFF_mcparticles_PDG[ipart]){
        h_particle_MC_pT[ipart]->Fill(truept,trueeta);
        h_chargedpart_MC_pT->Fill(truept,trueeta);

        h_particle_MC_p[ipart]->Fill(truep,trueeta);
        h_chargedpart_MC_p->Fill(truep,trueeta);
        h_particle_MC_E[ipart]->Fill(_mcpart_E[imc],trueeta);
        h_chargedpart_MC_E->Fill(_mcpart_E[imc],trueeta);
      }
    }
  }

  // see if true reco. track was found
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]-1],_mcpart_py[(int)_track_trueID[itrk]-1],_mcpart_pz[(int)_track_trueID[itrk]-1]);
    // cout << itrk << endl;
    float receta = recpartvec.Eta();
    float trueeta = mcpartvec.Eta();
    // cout << "\tTRKTRUEID " << _track_trueID[itrk]-1 << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]-1] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << endl;
    float pt = recpartvec.Pt();
    float truept = mcpartvec.Pt();
    float pmom = recpartvec.Mag();
    float truepmom = mcpartvec.Mag();
    TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1]));
    Double_t mass = part->Mass();
    float trueE     = TMath::Sqrt(truepmom*truepmom+mass*mass);
    
    for(int ipart=0;ipart<nPart_TRKEFF;ipart++){
      if(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1])==int_TRKEFF_mcparticles_PDG[ipart]){
        h_particle_rec_pT[ipart]->Fill(pt,trueeta);
        h_particle_rec_truepT[ipart]->Fill(truept,trueeta);
        h_chargedpart_rec_pT->Fill(pt,trueeta);
        h_chargedpart_rec_truepT->Fill(truept,trueeta);

        h_particle_rec_p[ipart]->Fill(pmom,trueeta);
        h_particle_rec_truep[ipart]->Fill(truepmom,trueeta);
        h_particle_rec_trueE[ipart]->Fill(trueE,trueeta);
        h_chargedpart_rec_p->Fill(pmom,trueeta);
        h_chargedpart_rec_truep->Fill(truepmom,trueeta);
        h_chargedpart_rec_trueE->Fill(trueE,trueeta);
      }
    }
  }
}


// **********************************************************************************************
// **************************** tracking resolution processing **********************************
// **********************************************************************************************
void trackingresolution(){
  // =======================================================================
  // **************** initializing hists ***********************************
  // =======================================================================
  if (!initResoHist){
    for (Int_t ipart = 0; ipart < nPart_TRKEFF+1; ipart++){
      for (Int_t k = 0; k < 12; k++){
//         std::cout << ipart << "\t" << str_TRKEFF_mcparticles[ipart].Data() << "\t" << k << "\t" << nameCuts[k].Data()  << "\t" 
//                   << Form("h_tracks_%s_%s_True_Eta_pT", str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data()) ;
        
        h_tracksTrue_Eta_pT[ipart][k]         = new TH2F(Form("h_tracks_%s_%s_True_Eta_pT", str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data() ), "; #it{p}_{T,MC} (GeV/c); #eta_{MC}", 200, 0, 20, 200, -4, 4.);
        h_tracksRec_Eta_pT[ipart][k]          = new TH2F(Form("h_tracks_%s_%s_Rec_Eta_pT", str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data() ), "; #it{p}_{T,rec} (GeV/c); #eta_{rec}", 200, 0, 20, 200, -4, 4.);
        h_tracksTrue_Eta_p[ipart][k]          = new TH2F(Form("h_tracks_%s_%s_True_Eta_p", str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data() ), "; #it{p}_{MC} (GeV/c); #eta_{MC}", nBinsP, binningP,  200, -4, 4.);
        h_tracksRec_Eta_p[ipart][k]           = new TH2F(Form("h_tracks_%s_%s_Rec_Eta_p", str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data() ), "; #it{p}_{rc} (GeV/c); #eta_{rec}", nBinsP, binningP, 200, -4, 4.);
//         std::cout << h_tracksTrue_Eta_pT[ipart][k]->GetName() << "\t" << h_tracksRec_Eta_pT[ipart][k]->GetName()  << "\t" << h_tracksTrue_Eta_p[ipart][k]->GetName() << "\t" << h_tracksRec_Eta_p[ipart][k]->GetName()<< std::endl;
      }
    }
    for (Int_t et = 0; et<nEta+1; et++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (et < nEta){
        etaMin = partEta[et];
        etaMax = partEta[et+1];
      }      

      for (Int_t i = 0; i< 5; i++){
        for (Int_t ipart = 0; ipart < 6; ipart++){
          h_tracks_reso_pT[i][et][ipart]     = new TH2F(Form("h_tracks_reso_pT_%s_%s_%d",nameResoAdd[i].Data(), str_TRKEFF_mcparticles[ipart].Data(),et), 
                                                Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#it{p}_{T,rec}-#it{p}_{T,MC})/#it{p}_{T,MC}",etaMin, etaMax), 
                                                200, 0, 20, 250, -0.5, 0.5);
          h_tracks_reso_p[i][et][ipart]         = new TH2F(Form("h_tracks_reso_p_%s_%s_%d",nameResoAdd[i].Data(), str_TRKEFF_mcparticles[ipart].Data(),et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (#it{p}_{rec}-#it{p}_{MC})/#it{p}_{MC}", etaMin, etaMax), 
                                                  nBinsP, binningP, 250, -0.5, 0.5);
        }
        h_tracks_resoEta_pT[i][et]     = new TH2F(Form("h_tracks_reso_Eta_pT_%s_%d",nameResoAdd[i].Data(),et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#eta_{rec}-#eta_{MC})/#eta_{MC}",
                                                                                        etaMin, etaMax), 200, 0, 20, 500, -0.1, 0.1);
        h_tracks_resoPhi_pT[i][et]     = new TH2F(Form("h_tracks_reso_Phi_pT_%s_%d",nameResoAdd[i].Data(),et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#phi_{rec}-#phi_{MC})/#phi_{MC}",
                                                                                        etaMin, etaMax), 200, 0, 20, 500, -0.5, 0.5);
      }
    }
    initResoHist = kTRUE;
  }
  hNEvents->Fill(1);
  hNTracks->Fill(_nTracks);
  
  // =======================================================================
  // ******************* filling resolution histos *************************
  // =======================================================================
  Int_t nCurrProj = 0; 
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]-1],_mcpart_py[(int)_track_trueID[itrk]-1],_mcpart_pz[(int)_track_trueID[itrk]-1]);
    float receta    = recpartvec.Eta();
    float trueeta   = mcpartvec.Eta();
    float recphi    = recpartvec.Phi();
    float truephi   = mcpartvec.Phi();
//     cout << "\tTRKTRUEID " << _track_trueID[itrk]-1 << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]-1] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << endl;
    float pt        = recpartvec.Pt();
    float truept    = mcpartvec.Pt();
    float pmom      = recpartvec.Mag();
    float truepmom  = mcpartvec.Mag();
    
    Int_t parIdx = 5;
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1]) == 11)   parIdx = 0;
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1]) == 211)  parIdx = 1;
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1]) == 2212) parIdx = 2;          
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1]) == 321)  parIdx = 3;
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1]) == 13)   parIdx = 4;
    
    TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1]));
    Double_t mass = part->Mass();
    float trueE     = TMath::Sqrt(truepmom*truepmom+mass*mass);
    
    Bool_t hasTL      = kFALSE;
    Bool_t hasTLAE    = kFALSE;
    Bool_t hasFTrL    = kFALSE;
    Bool_t hasTrL     = kFALSE;
    Int_t nTrL        = 0;
    Int_t nTL         = 0;
    Int_t nTrT        = 0;
    
    for(Int_t iproj=nCurrProj; iproj<_nProjections; iproj++){
      if (itrk != _track_ProjTrackID[iproj]) continue;
      if(_track_Proj_t[iproj]==0.) continue;
//       cout << "\t layer: "<< _track_ProjLayer[iproj] << "\t" << _track_Proj_x[iproj] << "\t" << _track_Proj_y[iproj] << "\t" << _track_Proj_z[iproj] << endl;
      hasTL     = (hasTL || HasTimingLayer(_track_ProjLayer[iproj]));
      if (HasTimingLayer(_track_ProjLayer[iproj])) nTL++;
        
      hasTLAE   = (hasTLAE || HasTimingLayerAfterECal(_track_ProjLayer[iproj]));
      hasFTrL   = (hasFTrL || HasFirstTwoLayers(_track_ProjLayer[iproj]));
      hasTrL    = (hasTrL || IsTrackerLayer(_track_ProjLayer[iproj]));
      if (IsTrackerLayer(_track_ProjLayer[iproj])) nTrL++;
      
      nCurrProj = iproj;
    }
    nTrT = nTrL + nTL;
//     cout << "\t summary:" << hasTL << "\t" << hasTLAE << "\t"<< nTL << "\t" << hasTrL << "\t" << hasFTrL << "\t" << nTrL << "\t tot: " << nTrT << endl;
      

    // determine eta bin
    Int_t et = 0;
    while (partEta[et+1] < trueeta && et < nEta) et++;
    
    
    h_tracks_reso_pT[0][et][parIdx]->Fill(truept,(pt-truept)/truept);
    h_tracks_reso_p[0][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);
    h_tracks_resoEta_pT[0][et]->Fill(truept,(receta-trueeta));
    h_tracks_resoPhi_pT[0][et]->Fill(truept,(recphi-truephi));
    if (nTL > 0){
      h_tracks_reso_pT[2][et][parIdx]->Fill(truept,(pt-truept)/truept);
      h_tracks_reso_p[2][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);
      h_tracks_resoEta_pT[2][et]->Fill(truept,(receta-trueeta));
      h_tracks_resoPhi_pT[2][et]->Fill(truept,(recphi-truephi));
    } else {
      h_tracks_reso_pT[1][et][parIdx]->Fill(truept,(pt-truept)/truept);
      h_tracks_reso_p[1][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);        
      h_tracks_resoEta_pT[1][et]->Fill(truept,(receta-trueeta));
      h_tracks_resoPhi_pT[1][et]->Fill(truept,(recphi-truephi));
    }
    if (nTrL >= 2){
      h_tracks_reso_pT[3][et][parIdx]->Fill(truept,(pt-truept)/truept);
      h_tracks_reso_p[3][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);        
      h_tracks_resoEta_pT[3][et]->Fill(truept,(receta-trueeta));
      h_tracks_resoPhi_pT[3][et]->Fill(truept,(recphi-truephi));
    }
    if (nTrL >= 3){
      h_tracks_reso_pT[4][et][parIdx]->Fill(truept,(pt-truept)/truept);
      h_tracks_reso_p[4][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);        
      h_tracks_resoEta_pT[4][et]->Fill(truept,(receta-trueeta));
      h_tracks_resoPhi_pT[4][et]->Fill(truept,(recphi-truephi));
    }
        
    // all tracks "N"
    h_tracksTrue_Eta_pT[parIdx][0]->Fill(truept, trueeta);
    h_tracksRec_Eta_pT[parIdx][0]->Fill(pt, receta);
    h_tracksTrue_Eta_p[parIdx][0]->Fill(truepmom, trueeta);
    h_tracksRec_Eta_p[parIdx][0]->Fill(pmom, receta);
    // at least 3 layers
    if (nTrT >= 3){
      h_tracksTrue_Eta_pT[parIdx][1]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[parIdx][1]->Fill(pt, receta);
      h_tracksTrue_Eta_p[parIdx][1]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[parIdx][1]->Fill(pmom, receta);
      // at least 3 layers & inner most tracking layer fw
      if (hasFTrL){ 
        h_tracksTrue_Eta_pT[parIdx][2]->Fill(truept, trueeta);
        h_tracksRec_Eta_pT[parIdx][2]->Fill(pt, receta);
        h_tracksTrue_Eta_p[parIdx][2]->Fill(truepmom, trueeta);
        h_tracksRec_Eta_p[parIdx][2]->Fill(pmom, receta);
      }
    }
    // only tracker layers
    if (nTL < 1){
      h_tracksTrue_Eta_pT[parIdx][3]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[parIdx][3]->Fill(pt, receta);
      h_tracksTrue_Eta_p[parIdx][3]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[parIdx][3]->Fill(pmom, receta);
      // only tracker layers & inner most tracking layer fw
      if (hasFTrL){ 
        h_tracksTrue_Eta_pT[parIdx][4]->Fill(truept, trueeta);
        h_tracksRec_Eta_pT[parIdx][4]->Fill(pt, receta);
        h_tracksTrue_Eta_p[parIdx][4]->Fill(truepmom, trueeta);
        h_tracksRec_Eta_p[parIdx][4]->Fill(pmom, receta);
      }        
    }
    // timing layer hit before ECal but not after
    if (hasTL && !hasTLAE){
      h_tracksTrue_Eta_pT[parIdx][5]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[parIdx][5]->Fill(pt, receta);
      h_tracksTrue_Eta_p[parIdx][5]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[parIdx][5]->Fill(pmom, receta);
      // timing layer hit before ECal but not after & inner most tracking layer fw
      if (hasFTrL){ 
        h_tracksTrue_Eta_pT[parIdx][6]->Fill(truept, trueeta);
        h_tracksRec_Eta_pT[parIdx][6]->Fill(pt, receta);
        h_tracksTrue_Eta_p[parIdx][6]->Fill(truepmom, trueeta);
        h_tracksRec_Eta_p[parIdx][6]->Fill(pmom, receta);
      }        
    }
    // timing layer hit after ECal 
    if (hasTLAE){
      h_tracksTrue_Eta_pT[parIdx][7]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[parIdx][7]->Fill(pt, receta);
      h_tracksTrue_Eta_p[parIdx][7]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[parIdx][7]->Fill(pmom, receta);
      // timing layer hit after ECal & inner most tracking layer fw
      if (hasFTrL){ 
        h_tracksTrue_Eta_pT[parIdx][8]->Fill(truept, trueeta);
        h_tracksRec_Eta_pT[parIdx][8]->Fill(pt, receta);
        h_tracksTrue_Eta_p[parIdx][8]->Fill(truepmom, trueeta);
        h_tracksRec_Eta_p[parIdx][8]->Fill(pmom, receta);
      }        
    }
  
    if (nTL > 0){
      h_tracksTrue_Eta_pT[parIdx][9]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[parIdx][9]->Fill(pt, receta);
      h_tracksTrue_Eta_p[parIdx][9]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[parIdx][9]->Fill(pmom, receta);
    }
    if (nTrL >= 2){
      h_tracksTrue_Eta_pT[parIdx][10]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[parIdx][10]->Fill(pt, receta);
      h_tracksTrue_Eta_p[parIdx][10]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[parIdx][10]->Fill(pmom, receta);
    }
    if (nTrL >= 3){
      h_tracksTrue_Eta_pT[parIdx][11]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[parIdx][11]->Fill(pt, receta);
      h_tracksTrue_Eta_p[parIdx][11]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[parIdx][11]->Fill(pmom, receta);
    }  
  }
}

// **********************************************************************************************
// **************************** saving efficiency histos ****************************************
// **********************************************************************************************
// ANCHOR save function after event loop
void trackingefficiencyhistosSave(){

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_TRKEFF.root",outputDir.Data()),"RECREATE");

  // write histograms
  if(h_chargedpart_MC_pT) h_chargedpart_MC_pT->Write();
  if(h_chargedpart_rec_pT) h_chargedpart_rec_pT->Write();
  if(h_chargedpart_rec_truepT) h_chargedpart_rec_truepT->Write();
  if(h_chargedpart_MC_pT && h_chargedpart_rec_pT && h_chargedpart_rec_truepT) {
    TH2F*  h_chargedpart_eff_pT = (TH2F*)h_chargedpart_rec_pT->Clone("h_chargedpart_eff_pT");
    h_chargedpart_eff_pT->Divide(h_chargedpart_MC_pT);
    h_chargedpart_eff_pT->Write("h_chargedpart_eff_pT", TObject::kOverwrite);
    TH2F*  h_chargedpart_eff_truepT = (TH2F*)h_chargedpart_rec_truepT->Clone("h_chargedpart_eff_truepT");
    h_chargedpart_eff_truepT->Divide(h_chargedpart_MC_pT);
    h_chargedpart_eff_truepT->Write("h_chargedpart_eff_truepT", TObject::kOverwrite);
  }
  if(h_chargedpart_MC_p) h_chargedpart_MC_p->Write();
  if(h_chargedpart_rec_p) h_chargedpart_rec_p->Write();
  if(h_chargedpart_rec_truep) h_chargedpart_rec_truep->Write();
  if(h_chargedpart_rec_trueE) h_chargedpart_rec_trueE->Write();
  if(h_chargedpart_MC_p && h_chargedpart_rec_p && h_chargedpart_rec_truep) {
    TH2F*  h_chargedpart_eff_p = (TH2F*)h_chargedpart_rec_p->Clone("h_chargedpart_eff_p");
    h_chargedpart_eff_p->Divide(h_chargedpart_MC_p);
    h_chargedpart_eff_p->Write("h_chargedpart_eff_p",TObject::kOverwrite);
    TH2F*  h_chargedpart_eff_truep = (TH2F*)h_chargedpart_rec_truep->Clone("h_chargedpart_eff_truep");
    h_chargedpart_eff_truep->Divide(h_chargedpart_MC_p);
    h_chargedpart_eff_truep->Write("h_chargedpart_eff_truep",TObject::kOverwrite);
  }
  if(h_chargedpart_MC_E) h_chargedpart_MC_E->Write();
  TH2F*  h_particle_eff_pT[nPart_TRKEFF];
  TH2F*  h_particle_eff_truepT[nPart_TRKEFF];
  TH2F*  h_particle_eff_p[nPart_TRKEFF];
  TH2F*  h_particle_eff_truep[nPart_TRKEFF];
  for(int ipart=0;ipart<nPart_TRKEFF;ipart++){
    if(h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart]->Write();
    if(h_particle_rec_pT[ipart]) h_particle_rec_pT[ipart]->Write();
    if(h_particle_rec_truepT[ipart]) h_particle_rec_truepT[ipart]->Write();
    if(h_particle_MC_pT[ipart] && h_particle_rec_pT[ipart] && h_particle_rec_truepT[ipart]) {
      h_particle_eff_pT[ipart] = (TH2F*)h_particle_rec_pT[ipart]->Clone(Form("h_%s_eff_pT",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_pT[ipart]->Divide(h_particle_MC_pT[ipart]);
      h_particle_eff_pT[ipart]->Write(Form("h_%s_eff_pT",str_TRKEFF_mcparticles[ipart].Data()),TObject::kOverwrite);
      h_particle_eff_truepT[ipart] = (TH2F*)h_particle_rec_truepT[ipart]->Clone(Form("h_%s_eff_truepT",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_truepT[ipart]->Divide(h_particle_MC_pT[ipart]);
      h_particle_eff_truepT[ipart]->Write(Form("h_%s_eff_truepT",str_TRKEFF_mcparticles[ipart].Data()),TObject::kOverwrite);
    }
    if(h_particle_MC_p[ipart]) h_particle_MC_p[ipart]->Write();
    if(h_particle_rec_p[ipart]) h_particle_rec_p[ipart]->Write();
    if(h_particle_rec_truep[ipart]) h_particle_rec_truep[ipart]->Write();
    if(h_particle_rec_trueE[ipart]) h_particle_rec_trueE[ipart]->Write();
    if(h_particle_MC_p[ipart] && h_particle_rec_p[ipart] && h_particle_rec_truep[ipart]) {
      h_particle_eff_p[ipart] = (TH2F*)h_particle_rec_p[ipart]->Clone(Form("h_%s_eff_p",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_p[ipart]->Divide(h_particle_MC_p[ipart]);
      h_particle_eff_p[ipart]->Write(Form("h_%s_eff_p",str_TRKEFF_mcparticles[ipart].Data()),TObject::kOverwrite);
      h_particle_eff_truep[ipart] = (TH2F*)h_particle_rec_truep[ipart]->Clone(Form("h_%s_eff_truep",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_truep[ipart]->Divide(h_particle_MC_p[ipart]);
      h_particle_eff_truep[ipart]->Write(Form("h_%s_eff_truep",str_TRKEFF_mcparticles[ipart].Data()),TObject::kOverwrite);
    }
    if(h_particle_MC_E[ipart]) h_particle_MC_E[ipart]->Write();
    
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}

// **********************************************************************************************
// **************************** saving resolution histos ****************************************
// **********************************************************************************************
// ANCHOR save function after event loop
void trackingresolutionhistosSave(){

  // ***********************************************************************************************
  // ************************** sum up hists correctly *********************************************
  // ***********************************************************************************************
  for (Int_t id = 0; id < nPart_TRKEFF; id++){
    for (Int_t k = 0; k < nCuts; k++){
      if (id == 1){
        h_tracksTrue_Eta_pT[nPart_TRKEFF][k]->Sumw2();
        h_tracksRec_Eta_pT[nPart_TRKEFF][k]->Sumw2();
        h_tracksTrue_Eta_p[nPart_TRKEFF][k]->Sumw2();
        h_tracksRec_Eta_p[nPart_TRKEFF][k]->Sumw2();
      }
//       if (k == 0) cout << id << "\t" << h_tracksTrue_Eta_pT[nPart_TRKEFF][k]->GetEntries() << "\t" << h_tracksTrue_Eta_pT[id][k]->GetEntries() << endl;
      h_tracksTrue_Eta_pT[nPart_TRKEFF][k]->Add(h_tracksTrue_Eta_pT[id][k]);
      h_tracksRec_Eta_pT[nPart_TRKEFF][k]->Add(h_tracksRec_Eta_pT[id][k]);
      h_tracksTrue_Eta_p[nPart_TRKEFF][k]->Add(h_tracksTrue_Eta_p[id][k]);
      h_tracksRec_Eta_p[nPart_TRKEFF][k]->Add(h_tracksRec_Eta_p[id][k]);     
    }
    for (Int_t i = 0; i < nResoSt; i++){
      for (Int_t et = 0; et < nEta; et++){
        if (id == 0 && et == 0){
          h_tracks_reso_pT[i][et][nPart_TRKEFF]->Sumw2();
          h_tracks_reso_pT[i][nEta][nPart_TRKEFF]->Sumw2();
          h_tracks_reso_pT[i][nEta][id]->Sumw2();
          h_tracks_reso_p[i][et][nPart_TRKEFF]->Sumw2();
          h_tracks_reso_p[i][nEta][nPart_TRKEFF]->Sumw2();
          h_tracks_reso_p[i][nEta][id]->Sumw2();
        }
        h_tracks_reso_pT[i][et][nPart_TRKEFF]->Add(h_tracks_reso_pT[i][et][id]);
        h_tracks_reso_pT[i][nEta][nPart_TRKEFF]->Add(h_tracks_reso_pT[i][et][id]);
        h_tracks_reso_pT[i][nEta][id]->Add(h_tracks_reso_pT[i][et][id]);
        h_tracks_reso_p[i][et][nPart_TRKEFF]->Add(h_tracks_reso_p[i][et][id]);
        h_tracks_reso_p[i][nEta][nPart_TRKEFF]->Add(h_tracks_reso_p[i][et][id]);
        h_tracks_reso_p[i][nEta][id]->Add(h_tracks_reso_p[i][et][id]);
      }
    }
  }
  for (Int_t i = 0; i < nResoSt; i++){
    for (Int_t et = 0; et < nEta; et++){
      if (et == 0){
        h_tracks_resoEta_pT[i][nEta]->Sumw2();
        h_tracks_resoPhi_pT[i][nEta]->Sumw2();
      }
      h_tracks_resoEta_pT[i][nEta]->Add(h_tracks_resoEta_pT[i][et]);
      h_tracks_resoPhi_pT[i][nEta]->Add(h_tracks_resoPhi_pT[i][et]);
    }
  }

  
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_TRKRS.root",outputDir.Data()),"RECREATE");
    hNEvents->Write();
    hNTracks->Write();
    for (Int_t et = 0; et < nEta+1; et++){
      for (Int_t i = 0; i< nResoSt; i++){
        for (Int_t id = 0; id < nPart_TRKEFF+1 ; id++) {
          if(h_tracks_reso_pT[i][et][id]) h_tracks_reso_pT[i][et][id]->Write();
          if(h_tracks_reso_p[i][et][id]) h_tracks_reso_p[i][et][id]->Write();
        }
        if (h_tracks_resoEta_pT[i][et]) h_tracks_resoEta_pT[i][et]->Write();
        if (h_tracks_resoPhi_pT[i][et]) h_tracks_resoPhi_pT[i][et]->Write();
      }
    }
    
    for (Int_t id = 0; id < nPart_TRKEFF+1; id++){
      for (Int_t k = 0; k < nCuts; k++){
        if(h_tracksTrue_Eta_pT[id][k]) h_tracksTrue_Eta_pT[id][k]->Write();
        if(h_tracksRec_Eta_pT[id][k]) h_tracksRec_Eta_pT[id][k]->Write();
        if(h_tracksTrue_Eta_p[id][k]) h_tracksTrue_Eta_p[id][k]->Write();
        if(h_tracksRec_Eta_p[id][k]) h_tracksRec_Eta_p[id][k]->Write();
      }
    }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
