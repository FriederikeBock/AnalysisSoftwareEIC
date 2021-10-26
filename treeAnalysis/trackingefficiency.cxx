
// ANCHOR debug output verbosity
Int_t verbosityTRKEFF = 0;

// ANCHOR create histograms globally
const int nPart_TRKEFF = 5;
const unsigned int nTrackSources = 3;
TString str_TRKEFF_mcparticles[nPart_TRKEFF+1] = {"electron", "cpion", "proton", "ckaon", "muon", "all"};
int int_TRKEFF_mcparticles_PDG[nPart_TRKEFF] = {11 /*e*/,211  /*pi*/,  2212/*p*/,  321/*K*/,  13/*mu*/};

// **********************************************************************************************
// ********************************** Set up histos efficiency histos ***************************
// **********************************************************************************************
TH2F*  h_chargedpart_MC_pT        = new TH2F("h_chargedpart_MC_pT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_particle_MC_pT[nPart_TRKEFF] = {nullptr};
TH2F*  h_chargedpart_rec_pT[nTrackSources]       = {nullptr};
TH2F*  h_chargedpart_rec_truepT[nTrackSources]   = {nullptr};
TH2F*  h_particle_rec_pT[nTrackSources][nPart_TRKEFF] = {{nullptr}};
TH2F*  h_particle_rec_truepT[nTrackSources][nPart_TRKEFF] = {{nullptr}};

TH2F*  h_chargedpart_MC_p         = new TH2F("h_chargedpart_MC_p", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_particle_MC_p[nPart_TRKEFF] = {nullptr};
TH2F*  h_chargedpart_rec_p[nTrackSources]        = {nullptr};
TH2F*  h_chargedpart_rec_truep[nTrackSources]    = {nullptr};
TH2F*  h_chargedpart_rec_trueE[nTrackSources]    = {nullptr};
TH2F*  h_particle_rec_p[nTrackSources][nPart_TRKEFF] = {{nullptr}};
TH2F*  h_particle_rec_truep[nTrackSources][nPart_TRKEFF] = {{nullptr}};
TH2F*  h_particle_rec_trueE[nTrackSources][nPart_TRKEFF] = {{nullptr}};

TH2F*  h_chargedpart_MC_E         = new TH2F("h_chargedpart_MC_E", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_particle_MC_E[nPart_TRKEFF];

// **********************************************************************************************
// ********************************** Set up histos resolution histos ***************************
// **********************************************************************************************
TH1D* hNEvents                                = new TH1D("nEvents","",1,0.5,1.5);
TH1D* hNTracks[nTrackSources]                                           = {nullptr};
TH2F* h_tracks_reso_pT[nTrackSources][nResoSt][nEta+1][nPart_TRKEFF+1] = {{{{nullptr}}}};
TH2F* h_tracks_reso_p[nTrackSources][nResoSt][nEta+1][nPart_TRKEFF+1]  = {{{{nullptr}}}};
TH2F* h_tracks_resoEta_pT[nTrackSources][nResoSt][nEta+1]              = {{{nullptr}}};
TH2F* h_tracks_resoPhi_pT[nTrackSources][nResoSt][nEta+1]              = {{{nullptr}}};
TH2F* h_tracksTrue_Eta_pT[nTrackSources][nPart_TRKEFF+1][nCuts]        = {{{nullptr}}};
TH2F* h_tracksRec_Eta_pT[nTrackSources][nPart_TRKEFF+1][nCuts]         = {{{nullptr}}};
TH2F* h_tracksTrue_Eta_p[nTrackSources][nPart_TRKEFF+1][nCuts]         = {{{nullptr}}};
TH2F* h_tracksRec_Eta_p[nTrackSources][nPart_TRKEFF+1][nCuts]          = {{{nullptr}}};
Bool_t initResoHist                                     = kFALSE;

// **********************************************************************************************
// ********************************* Set up tracking comparison histos **************************
// **********************************************************************************************
// NOTE: We'll only fill in half of these, but it's easier to create the fully array, and then
//       we only initialize + write what we need.
TH2F* h_trackingComparison_p_eta[nTrackSources][nTrackSources] = {{nullptr}};
TH2F* h_trackingComparison_p_phi[nTrackSources][nTrackSources] = {{nullptr}};
TH2F* h_trackingComparison_p_pt[nTrackSources][nTrackSources] = {{nullptr}};
TH2F* h_trackingComparison_eta_p[nTrackSources][nTrackSources] = {{nullptr}};
TH2F* h_trackingComparison_eta_phi[nTrackSources][nTrackSources] = {{nullptr}};
TH2F* h_trackingComparison_eta_pt[nTrackSources][nTrackSources] = {{nullptr}};
bool initTrackingComparionHists = false;

// **********************************************************************************************
// ****************** create vectors for matching diff rec tracks to MC part ********************
// **********************************************************************************************
void prepareMCMatchInfo(){
  for(Int_t itrk=0; itrk<(Int_t)_nTracks; itrk++){
    if (verbosityTRKEFF > 2) cout << "processing track: " << itrk << endl;
    if (_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].size() > 0){
      if (verbosityTRKEFF > 2)
        cout << "adding rec track: \t" << itrk << "\t track source: " << _track_source[itrk] << " MC: "<< (int)_track_trueID[itrk] << endl;
    } else {
      if (verbosityTRKEFF > 2)
        cout << "found rec track: \t" << itrk << "\t track source: " << _track_source[itrk] << " MC: "<< (int)_track_trueID[itrk] << endl;
    }
    _mcpart_RecTrackIDs[(int)_track_trueID[itrk]].push_back(itrk);
  }
  Int_t nCurrProj = 0;
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    if (verbosityTRKEFF > 2) cout << "current track: " << itrk <<endl;
    unsigned short trackSource = _track_source[itrk];
    for(Int_t iproj=nCurrProj; iproj<_nProjections; iproj++){
      if (itrk != _track_ProjTrackID[iproj]) continue;
      if(_track_Proj_t[iproj]==0.) continue;
      if (verbosityTRKEFF > 2) cout << "\t prjection layer: "<< _track_ProjLayer[iproj] << "\t" << _track_Proj_x[iproj] << "\t" << _track_Proj_y[iproj] << "\t" << _track_Proj_z[iproj] << endl;
      _track_hasTTL[itrk]     = (_track_hasTTL[itrk] || HasTimingLayer(_track_ProjLayer[iproj]));
      if (HasTimingLayer(_track_ProjLayer[iproj])) _track_nTTL[itrk]++;
      if (verbosityTRKEFF > 2) cout << "timing layer count: " << _track_nTTL[itrk] << endl;
//       if (IsTrackerLayer(_track_ProjLayer[iproj])) nTrL++;
      nCurrProj = iproj;
    }

    if (trackSource == 1 )
      _track_hasOL[itrk] = (_track_hasOL[itrk] || true);
    if (trackSource == 2 )
      _track_hasIL[itrk]   = (_track_hasIL[itrk] || true);

    if ((int)_track_trueID[itrk] < 0) continue;
    if (_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].size() > 1){
      if (verbosityTRKEFF > 2) cout << "identified multiple tracks for: " << (int)_track_trueID[itrk] << " current track id: " << itrk << endl;
      for (int rtr = 0; rtr < _mcpart_RecTrackIDs[(int)_track_trueID[itrk]].size(); rtr++){
        if (verbosityTRKEFF > 2) cout << rtr << "\t MC id " << _mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr) << "\t track source\t " << _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] << endl;
        if ( _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] > trackSource &&  _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] == 1 )
          _track_hasOL[itrk] = (_track_hasOL[itrk] || true);
        if ( _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] > trackSource  &&  _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] == 2 )
          _track_hasIL[itrk]   = (_track_hasIL[itrk] || true);
      }
    } else {
      if (verbosityTRKEFF > 2) cout << "only one track found for: " << (int)_track_trueID[itrk] << "\t source: " << trackSource<< endl;
    }


  }
}

// **********************************************************************************************
// **************************** tracking efficiency processing **********************************
// **********************************************************************************************
void trackingefficiency(){
  // cout << "new evt" << endl;
  // Setup hists if necessary
  for (unsigned int iTrkSource = 0; iTrkSource < nTrackSources; iTrkSource++) {
    if (!h_chargedpart_rec_pT[iTrkSource]) {
      h_chargedpart_rec_pT[iTrkSource]       = new TH2F(TString::Format("h_chargedpart_rec_pT_%u", iTrkSource), "", 60, 0,30,80, -4.0,4.0);
    }
    if (!h_chargedpart_rec_truepT[iTrkSource]) {
      h_chargedpart_rec_truepT[iTrkSource]   = new TH2F(TString::Format("h_chargedpart_rec_truepT_%u", iTrkSource), "", 60, 0,30,80, -4.0,4.0);
    }
    if (!h_chargedpart_rec_p[iTrkSource]) {
      h_chargedpart_rec_p[iTrkSource]        = new TH2F(TString::Format("h_chargedpart_rec_p_%u", iTrkSource), "", nBinsP, binningP, 80, -4.0,4.0);
    }
    if (!h_chargedpart_rec_truep[iTrkSource]) {
      h_chargedpart_rec_truep[iTrkSource]    = new TH2F(TString::Format("h_chargedpart_rec_truep_%u", iTrkSource), "", nBinsP, binningP, 80, -4.0,4.0);
    }
    if (!h_chargedpart_rec_trueE[iTrkSource]) {
      h_chargedpart_rec_trueE[iTrkSource]    = new TH2F(TString::Format("h_chargedpart_rec_trueE_%u", iTrkSource), "", nBinsP, binningP, 80, -4.0,4.0);
    }
  }

  for(int imc=0; imc<_nMCPart; imc++){
    TVector3 mcpartvec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
    float trueeta = _mcpart_Eta[imc];
    float truept  = mcpartvec.Pt();
    float truep   = mcpartvec.Mag();

    if (verbosityTRKEFF) cout << "\tMCID " << imc << "\tiMC id: " << _mcpart_ID[imc] << "\tPDG: " << _mcpart_PDG[imc] << "\tpT: " << TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)) << "\tEta " << trueeta << endl;
    for(int ipart=0;ipart<nPart_TRKEFF;ipart++){
      // MC
      if(!h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart] 	= new TH2F(Form("h_%s_MC_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);
      if(!h_particle_MC_p[ipart]) h_particle_MC_p[ipart] 	= new TH2F(Form("h_%s_MC_p",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      if(!h_particle_MC_E[ipart]) h_particle_MC_E[ipart] 	= new TH2F(Form("h_%s_MC_E",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);

      // Rec
      for (unsigned int i = 0; i < nTrackSources; i++) {
        // pt
        if(!h_particle_rec_pT[i][ipart]) {
          h_particle_rec_pT[i][ipart] 	= new TH2F(Form("h_%s_rec_pT_%u",str_TRKEFF_mcparticles[ipart].Data(), i), "", 60, 0,30,80, -4.0,4.0);
        }
        if(!h_particle_rec_truepT[i][ipart]) {
          h_particle_rec_truepT[i][ipart] 	= new TH2F(Form("h_%s_rec_truepT_%u",str_TRKEFF_mcparticles[ipart].Data(), i), "", 60, 0,30,80, -4.0,4.0);
        }
        // p
        if(!h_particle_rec_p[i][ipart]) {
          h_particle_rec_p[i][ipart] 	= new TH2F(Form("h_%s_rec_p_%u",str_TRKEFF_mcparticles[ipart].Data(), i), "", nBinsP, binningP, 80, -4.0,4.0);
        }
        if(!h_particle_rec_truep[i][ipart]) {
          h_particle_rec_truep[i][ipart] 	= new TH2F(Form("h_%s_rec_truep_%u",str_TRKEFF_mcparticles[ipart].Data(), i), "", nBinsP, binningP, 80, -4.0,4.0);
        }
        // E
        if(!h_particle_rec_trueE[i][ipart]) {
          h_particle_rec_trueE[i][ipart] 	= new TH2F(Form("h_%s_rec_trueE_%u",str_TRKEFF_mcparticles[ipart].Data(), i), "", nBinsP, binningP, 80, -4.0,4.0);
        }
      }

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
  for(Int_t itrk=0; itrk<(Int_t)_nTracks; itrk++){
    unsigned short trackSource = _track_source[itrk];
    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
    // cout << itrk << endl;
    float receta = recpartvec.Eta();
    float trueeta = mcpartvec.Eta();
    if (verbosityTRKEFF) cout << "\tTRKTRUEID " << (int)_track_trueID[itrk] << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << endl;
    float pt = recpartvec.Pt();
    float truept = mcpartvec.Pt();
    float pmom = recpartvec.Mag();
    float truepmom = mcpartvec.Mag();
    TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(abs(_mcpart_PDG[(int)_track_trueID[itrk]]));
    Double_t mass = part->Mass();
    float trueE     = TMath::Sqrt(truepmom*truepmom+mass*mass);
    for(int ipart=0;ipart<(int)nPart_TRKEFF;ipart++){
      if(abs(_mcpart_PDG[(int)_track_trueID[itrk]])==int_TRKEFF_mcparticles_PDG[ipart]){
        h_particle_rec_pT[trackSource][ipart]->Fill(pt,trueeta);
        h_particle_rec_truepT[trackSource][ipart]->Fill(truept,trueeta);
        h_chargedpart_rec_pT[trackSource]->Fill(pt,trueeta);
        h_chargedpart_rec_truepT[trackSource]->Fill(truept,trueeta);

        h_particle_rec_p[trackSource][ipart]->Fill(pmom,trueeta);
        h_particle_rec_truep[trackSource][ipart]->Fill(truepmom,trueeta);
        h_particle_rec_trueE[trackSource][ipart]->Fill(trueE,trueeta);
        h_chargedpart_rec_p[trackSource]->Fill(pmom,trueeta);
        h_chargedpart_rec_truep[trackSource]->Fill(truepmom,trueeta);
        h_chargedpart_rec_trueE[trackSource]->Fill(trueE,trueeta);
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
    for (unsigned int iTrkSource = 0; iTrkSource < nTrackSources; iTrkSource++) {
      for (Int_t ipart = 0; ipart < nPart_TRKEFF+1; ipart++){
        for (Int_t k = 0; k < nCuts; k++){
//           std::cout << ipart << "\t" << str_TRKEFF_mcparticles[ipart].Data() << "\t" << k << "\t" << nameCuts[k].Data()  << "\t"
//                     << Form("h_tracks_%s_%s_True_Eta_pT", str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data()) ;

          h_tracksTrue_Eta_pT[iTrkSource][ipart][k]         = new TH2F(Form("h_tracks_%u_%s_%s_True_Eta_pT", iTrkSource, str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data() ), "; #it{p}_{T,MC} (GeV/c); #eta_{MC}", 200, 0, 20, 200, -4, 4.);
          h_tracksRec_Eta_pT[iTrkSource][ipart][k]          = new TH2F(Form("h_tracks_%u_%s_%s_Rec_Eta_pT", iTrkSource, str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data() ), "; #it{p}_{T,rec} (GeV/c); #eta_{rec}", 200, 0, 20, 200, -4, 4.);
          h_tracksTrue_Eta_p[iTrkSource][ipart][k]          = new TH2F(Form("h_tracks_%u_%s_%s_True_Eta_p", iTrkSource, str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data() ), "; #it{p}_{MC} (GeV/c); #eta_{MC}", nBinsP, binningP,  200, -4, 4.);
          h_tracksRec_Eta_p[iTrkSource][ipart][k]           = new TH2F(Form("h_tracks_%u_%s_%s_Rec_Eta_p", iTrkSource, str_TRKEFF_mcparticles[ipart].Data(), nameCuts[k].Data() ), "; #it{p}_{rc} (GeV/c); #eta_{rec}", nBinsP, binningP, 200, -4, 4.);
//           std::cout << h_tracksTrue_Eta_pT[ipart][k]->GetName() << "\t" << h_tracksRec_Eta_pT[ipart][k]->GetName()  << "\t" << h_tracksTrue_Eta_p[ipart][k]->GetName() << "\t" << h_tracksRec_Eta_p[ipart][k]->GetName()<< std::endl;
        }
      }
      for (Int_t et = 0; et<nEta+1; et++){
        Double_t etaMin = partEta[0];
        Double_t etaMax = partEta[nEta];
        if (et < nEta){
          etaMin = partEta[et];
          etaMax = partEta[et+1];
        }

        for (Int_t i = 0; i< nResoSt; i++){
          for (Int_t ipart = 0; ipart < nPart_TRKEFF+1; ipart++){
            h_tracks_reso_pT[iTrkSource][i][et][ipart]     = new TH2F(Form("h_tracks_reso_pT_%u_%s_%s_%d", iTrkSource, nameResoAdd[i].Data(), str_TRKEFF_mcparticles[ipart].Data(),et),
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#it{p}_{T,rec}-#it{p}_{T,MC})/#it{p}_{T,MC}",etaMin, etaMax),
                                                  200, 0, 20, 250, -0.5, 0.5);
            h_tracks_reso_p[iTrkSource][i][et][ipart]         = new TH2F(Form("h_tracks_reso_p_%u_%s_%s_%d", iTrkSource, nameResoAdd[i].Data(), str_TRKEFF_mcparticles[ipart].Data(),et),
                                                    Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (#it{p}_{rec}-#it{p}_{MC})/#it{p}_{MC}", etaMin, etaMax),
                                                    nBinsP, binningP, 250, -0.5, 0.5);
          }
          h_tracks_resoEta_pT[iTrkSource][i][et]     = new TH2F(Form("h_tracks_reso_Eta_pT_%u_%s_%d", iTrkSource, nameResoAdd[i].Data(),et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#eta_{rec}-#eta_{MC})/#eta_{MC}",
                                                                                          etaMin, etaMax), 200, 0, 20, 500, -0.1, 0.1);
          h_tracks_resoPhi_pT[iTrkSource][i][et]     = new TH2F(Form("h_tracks_reso_Phi_pT_%u_%s_%d", iTrkSource, nameResoAdd[i].Data(),et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#phi_{rec}-#phi_{MC})/#phi_{MC}",
                                                                                          etaMin, etaMax), 200, 0, 20, 500, -0.5, 0.5);
        }
      }
      hNTracks[iTrkSource] = new TH1D(TString::Format("nTracks_%u", iTrkSource),"",200,-0.5,199.5);
    }
    initResoHist = kTRUE;
  }
  // Fill n_events and n_tracks
  hNEvents->Fill(1);
  std::vector<unsigned int> nTracks(nTrackSources, 0);
  for (unsigned int i = 0; i < (unsigned int)_nTracks; i++) {
      ++nTracks[_track_source[i]];
  }
  // Fill each track source
  for (unsigned int i = 0; i < nTracks.size(); i++) {
    hNTracks[i]->Fill(nTracks[i]);
  }

  // =======================================================================
  // ******************* filling resolution histos *************************
  // =======================================================================
  Int_t nCurrProj = 0;
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
    float receta    = recpartvec.Eta();
    float trueeta   = mcpartvec.Eta();
    float recphi    = recpartvec.Phi();
    float truephi   = mcpartvec.Phi();
    if (verbosityTRKEFF > 2) cout << "\tTRKTRUEID " << (int)_track_trueID[itrk] << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << endl;
    float pt        = recpartvec.Pt();
    float truept    = mcpartvec.Pt();
    float pmom      = recpartvec.Mag();
    float truepmom  = mcpartvec.Mag();
    unsigned short trackSource = _track_source[itrk];

    Int_t parIdx = 5;
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]]) == 11)   parIdx = 0;
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]]) == 211)  parIdx = 1;
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]]) == 2212) parIdx = 2;
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]]) == 321)  parIdx = 3;
    if(abs(_mcpart_PDG[(int)_track_trueID[itrk]]) == 13)   parIdx = 4;

    TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(abs(_mcpart_PDG[(int)_track_trueID[itrk]]));
    Double_t mass = part->Mass();
    float trueE     = TMath::Sqrt(truepmom*truepmom+mass*mass);

    Bool_t hasTL      = _track_hasTTL[itrk];
    Bool_t hasFTrL    = _track_hasIL[itrk];
    Bool_t hasTrL     = _track_hasOL[itrk];
    Int_t nTrL        = 0;
    Int_t nTL         = _track_nTTL[itrk];
    Int_t nTrT        = 0;

    for(Int_t iproj=nCurrProj; iproj<_nProjections; iproj++){
//       hasFTrL   = (hasFTrL || HasFirstTwoLayers(_track_ProjLayer[iproj]));
//       hasTrL    = (hasTrL || IsTrackerLayer(_track_ProjLayer[iproj]));
      if (IsTrackerLayer(_track_ProjLayer[iproj])) nTrL++;
      nCurrProj = iproj;
    }
    nTrT = nTrL + nTL;
    if (verbosityTRKEFF > 2) cout << "\t summary: TL " << hasTL  << "\t TL hits "<< nTL << "\t tr " << hasTrL << "\t tr silicon " << hasFTrL << "\t tr hits " << nTrL << "\t tot: " << nTrT << endl;


    // determine eta bin
    Int_t et = 0;
    while (partEta[et+1] < trueeta && et < nEta) et++;


    h_tracks_reso_pT[trackSource][0][et][parIdx]->Fill(truept,(pt-truept)/truept);
    h_tracks_reso_p[trackSource][0][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);
    h_tracks_resoEta_pT[trackSource][0][et]->Fill(truept,(receta-trueeta));
    h_tracks_resoPhi_pT[trackSource][0][et]->Fill(truept,(recphi-truephi));
    if (nTL > 0){
      h_tracks_reso_pT[trackSource][2][et][parIdx]->Fill(truept,(pt-truept)/truept);
      h_tracks_reso_p[trackSource][2][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);
      h_tracks_resoEta_pT[trackSource][2][et]->Fill(truept,(receta-trueeta));
      h_tracks_resoPhi_pT[trackSource][2][et]->Fill(truept,(recphi-truephi));
    } else {
      h_tracks_reso_pT[trackSource][1][et][parIdx]->Fill(truept,(pt-truept)/truept);
      h_tracks_reso_p[trackSource][1][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);
      h_tracks_resoEta_pT[trackSource][1][et]->Fill(truept,(receta-trueeta));
      h_tracks_resoPhi_pT[trackSource][1][et]->Fill(truept,(recphi-truephi));
    }
//     if (nTrL >= 2){
    if (hasTrL){
      h_tracks_reso_pT[trackSource][3][et][parIdx]->Fill(truept,(pt-truept)/truept);
      h_tracks_reso_p[trackSource][3][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);
      h_tracks_resoEta_pT[trackSource][3][et]->Fill(truept,(receta-trueeta));
      h_tracks_resoPhi_pT[trackSource][3][et]->Fill(truept,(recphi-truephi));
    }
//     if (nTrL >= 3){
    if (hasFTrL){
      h_tracks_reso_pT[trackSource][4][et][parIdx]->Fill(truept,(pt-truept)/truept);
      h_tracks_reso_p[trackSource][4][et][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom);
      h_tracks_resoEta_pT[trackSource][4][et]->Fill(truept,(receta-trueeta));
      h_tracks_resoPhi_pT[trackSource][4][et]->Fill(truept,(recphi-truephi));
    }

    // all tracks "N"
    h_tracksTrue_Eta_pT[trackSource][parIdx][0]->Fill(truept, trueeta);
    h_tracksRec_Eta_pT[trackSource][parIdx][0]->Fill(pt, receta);
    h_tracksTrue_Eta_p[trackSource][parIdx][0]->Fill(truepmom, trueeta);
    h_tracksRec_Eta_p[trackSource][parIdx][0]->Fill(pmom, receta);
    // at least 3 layers
    if (hasTrL){
      h_tracksTrue_Eta_pT[trackSource][parIdx][1]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[trackSource][parIdx][1]->Fill(pt, receta);
      h_tracksTrue_Eta_p[trackSource][parIdx][1]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[trackSource][parIdx][1]->Fill(pmom, receta);
      // at least 3 layers & inner most tracking layer fw
      if (hasFTrL){
        h_tracksTrue_Eta_pT[trackSource][parIdx][2]->Fill(truept, trueeta);
        h_tracksRec_Eta_pT[trackSource][parIdx][2]->Fill(pt, receta);
        h_tracksTrue_Eta_p[trackSource][parIdx][2]->Fill(truepmom, trueeta);
        h_tracksRec_Eta_p[trackSource][parIdx][2]->Fill(pmom, receta);
      }
    }
    // only tracker layers
    if (nTL < 1){
      h_tracksTrue_Eta_pT[trackSource][parIdx][3]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[trackSource][parIdx][3]->Fill(pt, receta);
      h_tracksTrue_Eta_p[trackSource][parIdx][3]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[trackSource][parIdx][3]->Fill(pmom, receta);
      // only tracker layers & inner most tracking layer fw
      if (hasFTrL){
        h_tracksTrue_Eta_pT[trackSource][parIdx][4]->Fill(truept, trueeta);
        h_tracksRec_Eta_pT[trackSource][parIdx][4]->Fill(pt, receta);
        h_tracksTrue_Eta_p[trackSource][parIdx][4]->Fill(truepmom, trueeta);
        h_tracksRec_Eta_p[trackSource][parIdx][4]->Fill(pmom, receta);
      }
    }
//     // timing layer hit before ECal but not after
//     if (hasTL && !hasTLAE){
//       h_tracksTrue_Eta_pT[trackSource][parIdx][5]->Fill(truept, trueeta);
//       h_tracksRec_Eta_pT[trackSource][parIdx][5]->Fill(pt, receta);
//       h_tracksTrue_Eta_p[trackSource][parIdx][5]->Fill(truepmom, trueeta);
//       h_tracksRec_Eta_p[trackSource][parIdx][5]->Fill(pmom, receta);
//       // timing layer hit before ECal but not after & inner most tracking layer fw
//       if (hasFTrL){
//         h_tracksTrue_Eta_pT[trackSource][parIdx][6]->Fill(truept, trueeta);
//         h_tracksRec_Eta_pT[trackSource][parIdx][6]->Fill(pt, receta);
//         h_tracksTrue_Eta_p[trackSource][parIdx][6]->Fill(truepmom, trueeta);
//         h_tracksRec_Eta_p[trackSource][parIdx][6]->Fill(pmom, receta);
//       }
//     }
//     // timing layer hit after ECal
//     if (hasTLAE){
//       h_tracksTrue_Eta_pT[trackSource][parIdx][7]->Fill(truept, trueeta);
//       h_tracksRec_Eta_pT[trackSource][parIdx][7]->Fill(pt, receta);
//       h_tracksTrue_Eta_p[trackSource][parIdx][7]->Fill(truepmom, trueeta);
//       h_tracksRec_Eta_p[trackSource][parIdx][7]->Fill(pmom, receta);
//       // timing layer hit after ECal & inner most tracking layer fw
//       if (hasFTrL){
//         h_tracksTrue_Eta_pT[trackSource][parIdx][8]->Fill(truept, trueeta);
//         h_tracksRec_Eta_pT[trackSource][parIdx][8]->Fill(pt, receta);
//         h_tracksTrue_Eta_p[trackSource][parIdx][8]->Fill(truepmom, trueeta);
//         h_tracksRec_Eta_p[trackSource][parIdx][8]->Fill(pmom, receta);
//       }
//     }

    if (nTL > 0){
      h_tracksTrue_Eta_pT[trackSource][parIdx][9]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[trackSource][parIdx][9]->Fill(pt, receta);
      h_tracksTrue_Eta_p[trackSource][parIdx][9]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[trackSource][parIdx][9]->Fill(pmom, receta);
    }
    if (hasTrL){
      h_tracksTrue_Eta_pT[trackSource][parIdx][10]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[trackSource][parIdx][10]->Fill(pt, receta);
      h_tracksTrue_Eta_p[trackSource][parIdx][10]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[trackSource][parIdx][10]->Fill(pmom, receta);
    }
    if (hasFTrL){
      h_tracksTrue_Eta_pT[trackSource][parIdx][11]->Fill(truept, trueeta);
      h_tracksRec_Eta_pT[trackSource][parIdx][11]->Fill(pt, receta);
      h_tracksTrue_Eta_p[trackSource][parIdx][11]->Fill(truepmom, trueeta);
      h_tracksRec_Eta_p[trackSource][parIdx][11]->Fill(pmom, receta);
    }
  }
}

// **********************************************************************************************
// ************************ compare tracking from different sources *****************************
// **********************************************************************************************
void trackingcomparison() {
  if (verbosityTRKEFF > 3) cout << "start comparison" << endl;
  std::map<unsigned int, std::string> trackMap = {
    {0, "all"},
    {1, "inner"},
    {2, "silicon"}
  };
  if (!initTrackingComparionHists) {
    for (unsigned int iOuterTrkSource = 0; iOuterTrkSource < nTrackSources; ++iOuterTrkSource) {
      for (unsigned int iInnerTrkSource = 0; iInnerTrkSource < nTrackSources; ++iInnerTrkSource) {
        if (iOuterTrkSource >= iInnerTrkSource) { continue; }
        h_trackingComparison_p_eta[iOuterTrkSource][iInnerTrkSource] = new TH2F(TString::Format("h_trackingComparison_p_eta_%u_%u", iOuterTrkSource, iInnerTrkSource), TString::Format(";p (GeV/c);#eta_{%s} - #eta_{%s}", trackMap[iOuterTrkSource].c_str(), trackMap[iInnerTrkSource].c_str()), nBinsP, binningP, 500, -0.2, 0.2);
        h_trackingComparison_p_phi[iOuterTrkSource][iInnerTrkSource] = new TH2F(TString::Format("h_trackingComparison_p_phi_%u_%u", iOuterTrkSource, iInnerTrkSource), TString::Format(";p (GeV/c);#varphi_{%s} - #varphi_{%s}", trackMap[iOuterTrkSource].c_str(), trackMap[iInnerTrkSource].c_str()), nBinsP, binningP, 500, -0.2, 0.2);
        h_trackingComparison_p_pt[iOuterTrkSource][iInnerTrkSource] = new TH2F(TString::Format("h_trackingComparison_p_pt_%u_%u", iOuterTrkSource, iInnerTrkSource), TString::Format(";p (GeV/c);p_{T,%s} - p_{T,%s} (GeV/c)", trackMap[iOuterTrkSource].c_str(), trackMap[iInnerTrkSource].c_str()), nBinsP, binningP, 500, -0.5, 0.5);
        h_trackingComparison_eta_p[iOuterTrkSource][iInnerTrkSource] = new TH2F(TString::Format("h_trackingComparison_eta_p_%u_%u", iOuterTrkSource, iInnerTrkSource), TString::Format(";#eta;p_{%s} - p_{%s} (GeV/c)", trackMap[iOuterTrkSource].c_str(), trackMap[iInnerTrkSource].c_str()), nEta, partEta, 500, -0.5, 0.5);
        h_trackingComparison_eta_phi[iOuterTrkSource][iInnerTrkSource] = new TH2F(TString::Format("h_trackingComparison_eta_phi_%u_%u", iOuterTrkSource, iInnerTrkSource), TString::Format(";#eta;#varphi_{%s} - #varphi_{%s}", trackMap[iOuterTrkSource].c_str(), trackMap[iInnerTrkSource].c_str()), nEta, partEta, 500, -0.5, 0.5);
        h_trackingComparison_eta_pt[iOuterTrkSource][iInnerTrkSource] = new TH2F(TString::Format("h_trackingComparison_eta_pt_%u_%u", iOuterTrkSource, iInnerTrkSource), TString::Format(";#eta;p_{T,%s} - p_{T,%s} (GeV/c)", trackMap[iOuterTrkSource].c_str(), trackMap[iInnerTrkSource].c_str()), nEta, partEta, 500, -0.5, 0.5);
      }
    }
    initTrackingComparionHists = true;
  }

  // First, iterate once over track sources and track index
  for (unsigned int iOuterTrk = 0; iOuterTrk < (unsigned int)_nTracks; ++iOuterTrk) {
    int iOuterTrkSource = _track_source[iOuterTrk];
    if (verbosityTRKEFF > 3)cout << "outer track: "<< iOuterTrk << endl;
    if ((int)_track_trueID[iOuterTrk] < 0) continue;
    if (_mcpart_RecTrackIDs[(int)_track_trueID[iOuterTrk]].size() > 1){
      for (int rtr = 0; rtr < _mcpart_RecTrackIDs[(int)_track_trueID[iOuterTrk]].size(); rtr++){
        int iInnerTrk = _mcpart_RecTrackIDs[(int)_track_trueID[iOuterTrk]].at(rtr);
        int iInnerTrkSource = _track_source[iInnerTrk];
        if ( _track_source[iOuterTrk] >= iInnerTrkSource ) continue;

        TVector3 outerVec(_track_px[iOuterTrk], _track_py[iOuterTrk], _track_pz[iOuterTrk]);
        TVector3 innerVec(_track_px[iInnerTrk], _track_py[iInnerTrk], _track_pz[iInnerTrk]);
        if (verbosityTRKEFF > 3) cout <<"O source: " << iOuterTrkSource<< " I: " <<  iInnerTrk<< " source I: " << iInnerTrkSource << " eta diff: " << outerVec.Eta() - innerVec.Eta() << " pt diff " << outerVec.Pt() - innerVec.Pt() << endl;
        h_trackingComparison_p_eta[iOuterTrkSource][iInnerTrkSource]->Fill(outerVec.Mag(), outerVec.Eta() - innerVec.Eta());
        h_trackingComparison_p_phi[iOuterTrkSource][iInnerTrkSource]->Fill(outerVec.Mag(), outerVec.Phi() - innerVec.Phi());
        h_trackingComparison_p_pt[iOuterTrkSource][iInnerTrkSource]->Fill(outerVec.Mag(), outerVec.Pt() - innerVec.Pt());
        h_trackingComparison_eta_p[iOuterTrkSource][iInnerTrkSource]->Fill(outerVec.Eta(), outerVec.Mag() - innerVec.Mag());
        h_trackingComparison_eta_phi[iOuterTrkSource][iInnerTrkSource]->Fill(outerVec.Eta(), outerVec.Phi() - innerVec.Phi());
        h_trackingComparison_eta_pt[iOuterTrkSource][iInnerTrkSource]->Fill(outerVec.Eta(), outerVec.Pt() - innerVec.Pt());
      }
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
  for (unsigned int i = 0; i < nTrackSources; i++) {
    if(h_chargedpart_rec_pT[i]) h_chargedpart_rec_pT[i]->Write();
    if(h_chargedpart_rec_truepT[i]) h_chargedpart_rec_truepT[i]->Write();
    if(h_chargedpart_MC_pT && h_chargedpart_rec_pT[i] && h_chargedpart_rec_truepT[i]) {
      std::string name = "h_chargedpart_eff_pT_" + std::to_string(i);
      TH2F*  h_chargedpart_eff_pT = (TH2F*)h_chargedpart_rec_pT[i]->Clone(name.c_str());
      h_chargedpart_eff_pT->Divide(h_chargedpart_MC_pT);
      h_chargedpart_eff_pT->Write(name.c_str(), TObject::kOverwrite);
      name = "h_chargedpart_eff_truepT_" + std::to_string(i);
      TH2F*  h_chargedpart_eff_truepT = (TH2F*)h_chargedpart_rec_truepT[i]->Clone(name.c_str());
      h_chargedpart_eff_truepT->Divide(h_chargedpart_MC_pT);
      h_chargedpart_eff_truepT->Write(name.c_str(), TObject::kOverwrite);
    }
  }
  if(h_chargedpart_MC_p) h_chargedpart_MC_p->Write();
  for (unsigned int i = 0; i < nTrackSources; i++) {
    if(h_chargedpart_rec_p[i]) h_chargedpart_rec_p[i]->Write();
    if(h_chargedpart_rec_truep[i]) h_chargedpart_rec_truep[i]->Write();
    if(h_chargedpart_rec_trueE[i]) h_chargedpart_rec_trueE[i]->Write();
    if(h_chargedpart_MC_p && h_chargedpart_rec_p[i] && h_chargedpart_rec_truep[i]) {
      std::string name = "h_chargedpart_eff_p_" + std::to_string(i);
      TH2F*  h_chargedpart_eff_p = (TH2F*)h_chargedpart_rec_p[i]->Clone(name.c_str());
      h_chargedpart_eff_p->Divide(h_chargedpart_MC_p);
      h_chargedpart_eff_p->Write(name.c_str(),TObject::kOverwrite);
      name = "h_chargedpart_eff_truep_" + std::to_string(i);
      TH2F*  h_chargedpart_eff_truep = (TH2F*)h_chargedpart_rec_truep[i]->Clone(name.c_str());
      h_chargedpart_eff_truep->Divide(h_chargedpart_MC_p);
      h_chargedpart_eff_truep->Write(name.c_str(),TObject::kOverwrite);
    }
  }
  if(h_chargedpart_MC_E) h_chargedpart_MC_E->Write();
  TH2F*  h_particle_eff_pT[nTrackSources][nPart_TRKEFF] = {{nullptr}};
  TH2F*  h_particle_eff_truepT[nTrackSources][nPart_TRKEFF] = {{nullptr}};
  TH2F*  h_particle_eff_p[nTrackSources][nPart_TRKEFF] = {{nullptr}};
  TH2F*  h_particle_eff_truep[nTrackSources][nPart_TRKEFF] = {{nullptr}};
  for(int ipart=0;ipart<nPart_TRKEFF;ipart++){
    if(h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart]->Write();
    for (unsigned int i = 0; i < nTrackSources; i++) {
      if(h_particle_rec_pT[i][ipart]) h_particle_rec_pT[i][ipart]->Write();
      if(h_particle_rec_truepT[i][ipart]) h_particle_rec_truepT[i][ipart]->Write();
      if(h_particle_MC_pT[ipart] && h_particle_rec_pT[i][ipart] && h_particle_rec_truepT[ipart]) {
        std::string name = std::string("h_") + str_TRKEFF_mcparticles[ipart].Data() + "_eff_pT" + std::to_string(i);
        h_particle_eff_pT[i][ipart] = (TH2F*)h_particle_rec_pT[i][ipart]->Clone(name.c_str());
        h_particle_eff_pT[i][ipart]->Divide(h_particle_MC_pT[ipart]);
        h_particle_eff_pT[i][ipart]->Write(name.c_str(), TObject::kOverwrite);
        name = std::string("h_") + str_TRKEFF_mcparticles[ipart].Data() + "eff_truepT" + std::to_string(i);
        h_particle_eff_truepT[i][ipart] = (TH2F*)h_particle_rec_truepT[i][ipart]->Clone(name.c_str());
        h_particle_eff_truepT[i][ipart]->Divide(h_particle_MC_pT[ipart]);
        h_particle_eff_truepT[i][ipart]->Write(name.c_str(), TObject::kOverwrite);
      }
    }
    if(h_particle_MC_p[ipart]) h_particle_MC_p[ipart]->Write();
    for (unsigned int i = 0; i < nTrackSources; i++) {
      if(h_particle_rec_p[i][ipart]) h_particle_rec_p[i][ipart]->Write();
      if(h_particle_rec_truep[i][ipart]) h_particle_rec_truep[i][ipart]->Write();
      if(h_particle_rec_trueE[i][ipart]) h_particle_rec_trueE[i][ipart]->Write();
      if(h_particle_MC_p[ipart] && h_particle_rec_p[i][ipart] && h_particle_rec_truep[i][ipart]) {
        std::string name = std::string("h_") + str_TRKEFF_mcparticles[ipart].Data() + "_eff_p" + std::to_string(i);
        h_particle_eff_p[i][ipart] = (TH2F*)h_particle_rec_p[i][ipart]->Clone(name.c_str());
        h_particle_eff_p[i][ipart]->Divide(h_particle_MC_p[ipart]);
        h_particle_eff_p[i][ipart]->Write(name.c_str(), TObject::kOverwrite);
        name = std::string("h_") + str_TRKEFF_mcparticles[ipart].Data() + "_eff_truep" + std::to_string(i);
        h_particle_eff_truep[i][ipart] = (TH2F*)h_particle_rec_truep[i][ipart]->Clone(name.c_str());
        h_particle_eff_truep[i][ipart]->Divide(h_particle_MC_p[ipart]);
        h_particle_eff_truep[i][ipart]->Write(name.c_str(), TObject::kOverwrite);
      }
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
  for (unsigned int iTrkSource = 0; iTrkSource < nTrackSources; iTrkSource++) {
    for (Int_t id = 0; id < nPart_TRKEFF; id++){
      for (Int_t k = 0; k < nCuts; k++){
        if (id == 1){
          if(h_tracksTrue_Eta_pT[iTrkSource][nPart_TRKEFF][k]) h_tracksTrue_Eta_pT[iTrkSource][nPart_TRKEFF][k]->Sumw2();
          if(h_tracksRec_Eta_pT[iTrkSource][nPart_TRKEFF][k]) h_tracksRec_Eta_pT[iTrkSource][nPart_TRKEFF][k]->Sumw2();
          if(h_tracksTrue_Eta_p[iTrkSource][nPart_TRKEFF][k]) h_tracksTrue_Eta_p[iTrkSource][nPart_TRKEFF][k]->Sumw2();
          if(h_tracksRec_Eta_p[iTrkSource][nPart_TRKEFF][k]) h_tracksRec_Eta_p[iTrkSource][nPart_TRKEFF][k]->Sumw2();
        }
//         if (k == 0) cout << id << "\t" << h_tracksTrue_Eta_pT[nPart_TRKEFF][k]->GetEntries() << "\t" << h_tracksTrue_Eta_pT[id][k]->GetEntries() << endl;
        if(h_tracksTrue_Eta_pT[iTrkSource][nPart_TRKEFF][k]) h_tracksTrue_Eta_pT[iTrkSource][nPart_TRKEFF][k]->Add(h_tracksTrue_Eta_pT[iTrkSource][id][k]);
        if(h_tracksRec_Eta_pT[iTrkSource][nPart_TRKEFF][k]) h_tracksRec_Eta_pT[iTrkSource][nPart_TRKEFF][k]->Add(h_tracksRec_Eta_pT[iTrkSource][id][k]);
        if(h_tracksTrue_Eta_p[iTrkSource][nPart_TRKEFF][k]) h_tracksTrue_Eta_p[iTrkSource][nPart_TRKEFF][k]->Add(h_tracksTrue_Eta_p[iTrkSource][id][k]);
        if(h_tracksRec_Eta_p[iTrkSource][nPart_TRKEFF][k]) h_tracksRec_Eta_p[iTrkSource][nPart_TRKEFF][k]->Add(h_tracksRec_Eta_p[iTrkSource][id][k]);
      }
      for (Int_t i = 0; i < nResoSt; i++){
        for (Int_t et = 0; et < nEta; et++){
          if (id == 0 && et == 0){
            if(h_tracks_reso_pT[iTrkSource][i][et][nPart_TRKEFF]) h_tracks_reso_pT[iTrkSource][i][et][nPart_TRKEFF]->Sumw2();
            if(h_tracks_reso_pT[iTrkSource][i][nEta][nPart_TRKEFF]) h_tracks_reso_pT[iTrkSource][i][nEta][nPart_TRKEFF]->Sumw2();
            if(h_tracks_reso_pT[iTrkSource][i][nEta][id]) h_tracks_reso_pT[iTrkSource][i][nEta][id]->Sumw2();
            if(h_tracks_reso_p[iTrkSource][i][et][nPart_TRKEFF]) h_tracks_reso_p[iTrkSource][i][et][nPart_TRKEFF]->Sumw2();
            if(h_tracks_reso_p[iTrkSource][i][nEta][nPart_TRKEFF]) h_tracks_reso_p[iTrkSource][i][nEta][nPart_TRKEFF]->Sumw2();
            if(h_tracks_reso_p[iTrkSource][i][nEta][id]) h_tracks_reso_p[iTrkSource][i][nEta][id]->Sumw2();
          }
          if(h_tracks_reso_pT[iTrkSource][i][et][id]){
            if(h_tracks_reso_pT[iTrkSource][i][et][nPart_TRKEFF]) h_tracks_reso_pT[iTrkSource][i][et][nPart_TRKEFF]->Add(h_tracks_reso_pT[iTrkSource][i][et][id]);
            if(h_tracks_reso_pT[iTrkSource][i][nEta][nPart_TRKEFF]) h_tracks_reso_pT[iTrkSource][i][nEta][nPart_TRKEFF]->Add(h_tracks_reso_pT[iTrkSource][i][et][id]);
            if(h_tracks_reso_pT[iTrkSource][i][nEta][id]) h_tracks_reso_pT[iTrkSource][i][nEta][id]->Add(h_tracks_reso_pT[iTrkSource][i][et][id]);
          }
          if(h_tracks_reso_p[iTrkSource][i][et][id]){
            if(h_tracks_reso_p[iTrkSource][i][et][nPart_TRKEFF]) h_tracks_reso_p[iTrkSource][i][et][nPart_TRKEFF]->Add(h_tracks_reso_p[iTrkSource][i][et][id]);
            if(h_tracks_reso_p[iTrkSource][i][nEta][nPart_TRKEFF]) h_tracks_reso_p[iTrkSource][i][nEta][nPart_TRKEFF]->Add(h_tracks_reso_p[iTrkSource][i][et][id]);
            if(h_tracks_reso_p[iTrkSource][i][nEta][id]) h_tracks_reso_p[iTrkSource][i][nEta][id]->Add(h_tracks_reso_p[iTrkSource][i][et][id]);
          }
        }
      }
    }
    for (Int_t i = 0; i < nResoSt; i++){
      for (Int_t et = 0; et < nEta; et++){
        if (et == 0){
          if(h_tracks_resoEta_pT[iTrkSource][i][nEta]) h_tracks_resoEta_pT[iTrkSource][i][nEta]->Sumw2();
          if(h_tracks_resoPhi_pT[iTrkSource][i][nEta]) h_tracks_resoPhi_pT[iTrkSource][i][nEta]->Sumw2();
        }
        if(h_tracks_resoEta_pT[iTrkSource][i][nEta]) h_tracks_resoEta_pT[iTrkSource][i][nEta]->Add(h_tracks_resoEta_pT[iTrkSource][i][et]);
        if(h_tracks_resoPhi_pT[iTrkSource][i][nEta]) h_tracks_resoPhi_pT[iTrkSource][i][nEta]->Add(h_tracks_resoPhi_pT[iTrkSource][i][et]);
      }
    }
  }

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_TRKRS.root",outputDir.Data()),"RECREATE");
  hNEvents->Write();
  for (unsigned int iTrkSource = 0; iTrkSource < nTrackSources; iTrkSource++) {
    fileOutput->mkdir(Form("TrackSource_%d", iTrkSource));
    fileOutput->cd(Form("TrackSource_%d", iTrkSource));
    if(hNTracks[iTrkSource])hNTracks[iTrkSource]->Write();
    for (Int_t et = 0; et < nEta+1; et++){
      for (Int_t i = 0; i< nResoSt; i++){
        for (Int_t id = 0; id < nPart_TRKEFF+1 ; id++) {
          if(h_tracks_reso_pT[iTrkSource][i][et][id]) h_tracks_reso_pT[iTrkSource][i][et][id]->Write();
          if(h_tracks_reso_p[iTrkSource][i][et][id]) h_tracks_reso_p[iTrkSource][i][et][id]->Write();
        }
        if (h_tracks_resoEta_pT[iTrkSource][i][et]) h_tracks_resoEta_pT[iTrkSource][i][et]->Write();
        if (h_tracks_resoPhi_pT[iTrkSource][i][et]) h_tracks_resoPhi_pT[iTrkSource][i][et]->Write();
      }
    }

    for (Int_t id = 0; id < nPart_TRKEFF+1; id++){
      for (Int_t k = 0; k < nCuts; k++){
        if(h_tracksTrue_Eta_pT[iTrkSource][id][k]) h_tracksTrue_Eta_pT[iTrkSource][id][k]->Write();
        if(h_tracksRec_Eta_pT[iTrkSource][id][k]) h_tracksRec_Eta_pT[iTrkSource][id][k]->Write();
        if(h_tracksTrue_Eta_p[iTrkSource][id][k]) h_tracksTrue_Eta_p[iTrkSource][id][k]->Write();
        if(h_tracksRec_Eta_p[iTrkSource][id][k]) h_tracksRec_Eta_p[iTrkSource][id][k]->Write();
      }
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}


// **********************************************************************************************
// ****************** saving comparison of tracking from different sources **********************
// **********************************************************************************************
void trackingcomparisonhistosSave() {

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_TRKS_Comparison.root", outputDir.Data()),"RECREATE");

  for (unsigned int iOuterTrkSource = 0; iOuterTrkSource < nTrackSources; ++iOuterTrkSource) {
    for (unsigned int iInnerTrkSource = 0; iInnerTrkSource < nTrackSources; ++iInnerTrkSource) {
      if (iOuterTrkSource >= iInnerTrkSource) { continue; }
      if (h_trackingComparison_p_eta[iOuterTrkSource][iInnerTrkSource]) { h_trackingComparison_p_eta[iOuterTrkSource][iInnerTrkSource]->Write(); }
      if (h_trackingComparison_p_phi[iOuterTrkSource][iInnerTrkSource]) { h_trackingComparison_p_phi[iOuterTrkSource][iInnerTrkSource]->Write(); }
      if (h_trackingComparison_p_pt[iOuterTrkSource][iInnerTrkSource]) { h_trackingComparison_p_pt[iOuterTrkSource][iInnerTrkSource]->Write(); }
      if (h_trackingComparison_eta_p[iOuterTrkSource][iInnerTrkSource]) { h_trackingComparison_eta_p[iOuterTrkSource][iInnerTrkSource]->Write(); }
      if (h_trackingComparison_eta_phi[iOuterTrkSource][iInnerTrkSource]) { h_trackingComparison_eta_phi[iOuterTrkSource][iInnerTrkSource]->Write(); }
      if (h_trackingComparison_eta_pt[iOuterTrkSource][iInnerTrkSource]) { h_trackingComparison_eta_pt[iOuterTrkSource][iInnerTrkSource]->Write(); }
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}


void clearMCRecMatchVectors(){
  if (verbosityTRKEFF > 2) cout << "clearing MC vectors" << endl;
  for(int imc=0; imc<_maxNMCPart; imc++){
    if (_mcpart_RecTrackIDs[imc].size() > 0){
      _mcpart_RecTrackIDs[imc].clear();
      _mcpart_RecTrackIDs[imc].resize(0);
    }
  }

  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    _track_hasTTL[itrk]  = 0;
    _track_nTTL[itrk]  = 0;
    _track_hasIL[itrk]   = 0;
    _track_hasOL[itrk]   = 0;
  }
}
