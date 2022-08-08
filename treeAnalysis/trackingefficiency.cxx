
// ANCHOR debug output verbosity
Int_t verbosityTRKEFF = 0;

// ANCHOR create histograms globally
const int nPart_TRKEFF = 5;
const unsigned int nTrackSources = 4;
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
TH3F* h_tracks_reso_pT[nTrackSources][nPart_TRKEFF+1]                   = {{nullptr}};
TH3F* h_tracks_reso_p[nTrackSources][nPart_TRKEFF+1]                    = {{nullptr}};
TH3F* h_tracks_resoEta_pT[nTrackSources]                                = {nullptr};
TH3F* h_tracks_resoPhi_pT[nTrackSources]                                = {nullptr};
TH3F* h_tracks_dca2d_pT[nTrackSources]                                  = {nullptr};
TH3F* h_tracks_dca2d_p[nTrackSources]                                   = {nullptr};
TH3F* h_tracks_chi2_p[nTrackSources]                                   = {nullptr};
TH3F* h_tracks_ndf_p[nTrackSources]                                   = {nullptr};
TH2F* h_tracksTrue_Eta_pT[nTrackSources][nPart_TRKEFF+1]                = {{nullptr}};
TH2F* h_tracksRec_Eta_pT[nTrackSources][nPart_TRKEFF+1]                 = {{nullptr}};
TH2F* h_tracksTrue_Eta_p[nTrackSources][nPart_TRKEFF+1]                 = {{nullptr}};
TH2F* h_tracksRec_Eta_p[nTrackSources][nPart_TRKEFF+1]                  = {{nullptr}};
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
// **************************** tracking efficiency processing **********************************
// **********************************************************************************************
void trackingefficiency(){
  // std::cout << "new evt" << std::endl;
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

    if (verbosityTRKEFF) std::cout << "\tMCID " << imc << "\tiMC id: " << _mcpart_ID[imc] << "\tPDG: " << _mcpart_PDG[imc] << "\tpT: " << TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)) << "\tEta " << trueeta << std::endl;
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
    if (trackSource >= (unsigned short)nTrackSources) continue;
    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
    // std::cout << itrk << std::endl;
    float receta = recpartvec.Eta();
    float trueeta = mcpartvec.Eta();
    if (verbosityTRKEFF) std::cout << "\tTRKTRUEID " << (int)_track_trueID[itrk] << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << "\t track source\t" << trackSource << std::endl;
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
    Double_t binningReso [251];
    for (Int_t kt = 0; kt < 250; kt++){
      binningReso[kt] = -0.5 + kt*(0.5-(-0.5))/250;  
      if (kt == 249)
          binningReso[kt+1] = 0.5;
    }
//     for (Int_t kt = 0; kt < 250+1; kt++) std::cout << binningReso[kt] <<","; 
//     std::cout << std::endl;
      
    Double_t binningEtaC [81];
    for (Int_t kt = 0; kt < 80; kt++){
      binningEtaC[kt] = -4. + kt*(4.-(-4.))/80;  
      if (kt == 79)
          binningEtaC[kt+1] = 4;
        
    }
//     for (Int_t kt = 0; kt < 80+1; kt++) std::cout << binningEtaC[kt] <<",";
//     std::cout << std::endl;
    Double_t minDCA = -500;
    Double_t maxDCA = 500;
    Double_t binningdca2D [1001];
    for (Int_t kt = 0; kt < 1000; kt++){
      binningdca2D[kt] = minDCA + kt*(maxDCA-minDCA)/1000;  
      if (kt == 1000-1)
          binningdca2D[kt+1] = maxDCA;   
    }
    Double_t minNDF = 0;
    Double_t maxNDF = 20;
    Double_t binningNDF [21];
    for (Int_t kt = 0; kt < 20; kt++){
      binningNDF[kt] = minNDF + kt*(maxNDF-minNDF)/20;  
      if (kt == 20-1)
          binningNDF[kt+1] = maxNDF;   
    }
//     for (Int_t kt = 0; kt < 20+1; kt++) std::cout << binningNDF[kt] <<",";
//     std::cout << std::endl;

    Double_t minChi2 = 0;
    Double_t maxChi2 = 50;
    Double_t binningChi2 [101];
    for (Int_t kt = 0; kt < 100; kt++){
      binningChi2[kt] = minChi2 + kt*(maxChi2-minChi2)/100;  
      if (kt == 100-1)
          binningChi2[kt+1] = maxChi2;   
    }
//     for (Int_t kt = 0; kt < 100+1; kt++) std::cout << binningChi2[kt] <<",";
//     std::cout << std::endl;

    for (unsigned int iTrkSource = 0; iTrkSource < nTrackSources; iTrkSource++) {
      for (Int_t ipart = 0; ipart < nPart_TRKEFF+1; ipart++){
        h_tracksTrue_Eta_pT[iTrkSource][ipart]         = new TH2F(Form("h_tracks_%u_%s_True_Eta_pT", iTrkSource, str_TRKEFF_mcparticles[ipart].Data() ), 
                                                                      "; #it{p}_{T,MC} (GeV/c); #eta_{MC}", 200, 0, 20, 200, -4, 4.);
        h_tracksRec_Eta_pT[iTrkSource][ipart]          = new TH2F(Form("h_tracks_%u_%s_Rec_Eta_pT", iTrkSource, str_TRKEFF_mcparticles[ipart].Data() ), 
                                                                      "; #it{p}_{T,rec} (GeV/c); #eta_{rec}", 200, 0, 20, 200, -4, 4.);
        h_tracksTrue_Eta_p[iTrkSource][ipart]          = new TH2F(Form("h_tracks_%u_%s_True_Eta_p", iTrkSource, str_TRKEFF_mcparticles[ipart].Data() ), 
                                                                      "; #it{p}_{MC} (GeV/c); #eta_{MC}", nBinsP, binningP,  200, -4, 4.);
        h_tracksRec_Eta_p[iTrkSource][ipart]           = new TH2F(Form("h_tracks_%u_%s_Rec_Eta_p", iTrkSource, str_TRKEFF_mcparticles[ipart].Data() ), 
                                                                       "; #it{p}_{rc} (GeV/c); #eta_{rec}", nBinsP, binningP, 200, -4, 4.);
      }
      for (Int_t ipart = 0; ipart < nPart_TRKEFF+1; ipart++){
        h_tracks_reso_pT[iTrkSource][ipart]   = new TH3F(Form("h_tracks_reso_pT_%u_%s", iTrkSource, str_TRKEFF_mcparticles[ipart].Data()),
                                                          "; #it{p}_{T,MC} (GeV/c); (#it{p}_{T,rec}-#it{p}_{T,MC})/#it{p}_{T,MC}; #eta",
                                                          200, 0, 20, 250, -0.5, 0.5, 80, -4, 4);
        h_tracks_reso_p[iTrkSource][ipart]    = new TH3F(Form("h_tracks_reso_p_%u_%s", iTrkSource, str_TRKEFF_mcparticles[ipart].Data()),
                                                          "; #it{p}_{MC} (GeV/c); (#it{p}_{rec}-#it{p}_{MC})/#it{p}_{MC}; #eta",
                                                          nBinsP, binningP, 250, binningReso, 80, binningEtaC);
      }
      h_tracks_resoEta_pT[iTrkSource]         = new TH3F(Form("h_tracks_reso_Eta_pT_%u", iTrkSource), "; #it{p}_{T,MC} (GeV/c); (#eta_{rec}-#eta_{MC})/#eta_{MC}; #eta", 
                                                         200, 0, 20, 500, -0.1, 0.1, 80, -4, 4);
      h_tracks_resoPhi_pT[iTrkSource]         = new TH3F(Form("h_tracks_reso_Phi_pT_%u", iTrkSource), "; #it{p}_{T,MC} (GeV/c); (#phi_{rec}-#phi_{MC})/#phi_{MC}; #eta",
                                                         200, 0, 20, 500, -0.5, 0.5, 80, -4, 4);
      if (iTrkSource == 0 || iTrkSource == 4 ){
        h_tracks_dca2d_pT[iTrkSource]           = new TH3F(Form("h_tracks_dca_pT_%u", iTrkSource), "; #it{p}_{T,MC} (GeV/c); dca_{2D} (#mum); #eta", 
                                                          200, 0, 20, 1000, minDCA, maxDCA, 80, -4, 4);
        h_tracks_dca2d_p[iTrkSource]            = new TH3F(Form("h_tracks_dca_p_%u", iTrkSource), "; #it{p}_{MC} (GeV/c); dca_{2D} (#mum); #eta",
                                                          nBinsP, binningP, 1000, binningdca2D, 80, binningEtaC);
        h_tracks_ndf_p[iTrkSource]              = new TH3F(Form("h_tracks_ndf_p_%u", iTrkSource), "; #it{p}_{MC} (GeV/c); ndf; #eta",
                                                          nBinsP, binningP, 20, binningNDF, 80, binningEtaC);
        h_tracks_chi2_p[iTrkSource]             = new TH3F(Form("h_tracks_chi2_p_%u", iTrkSource), "; #it{p}_{MC} (GeV/c); #chi^{2}; #eta",
                                                          nBinsP, binningP, 100, binningChi2, 80, binningEtaC);
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
    unsigned short trackSource = _track_source[itrk];
    if (trackSource >= (unsigned short)nTrackSources) continue;

    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
    float receta    = recpartvec.Eta();
    float trueeta   = mcpartvec.Eta();
    float recphi    = recpartvec.Phi();
    float truephi   = mcpartvec.Phi();
    if (verbosityTRKEFF > 2) std::cout << "\tTRKTRUEID " << (int)_track_trueID[itrk] << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << std::endl;
    float pt        = recpartvec.Pt();
    float truept    = mcpartvec.Pt();
    float pmom      = recpartvec.Mag();
    float truepmom  = mcpartvec.Mag();

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
    Int_t nTrL        = _track_nTrL[itrk];
    Int_t nTL         = _track_nTTL[itrk];
    Int_t nTrT        = _track_nTrL[itrk]-_track_nTTL[itrk];

    if (verbosityTRKEFF > 2) std::cout << "\t summary: TL " << hasTL  << "\t TL hits "<< nTL << "\t tr " << hasTrL << "\t tr silicon " << hasFTrL << "\t tr hits " << nTrL << "\t tot: " << nTrT << std::endl;
    if (verbosityTRKEFF > 2) std::cout << "\t summary: chi2 " << _track_chi2[itrk] << "\t ndf " << _track_ndf[itrk]<< std::endl;

    // determine eta bin
    h_tracks_reso_pT[trackSource][parIdx]->Fill(truept,(pt-truept)/truept,trueeta);
    h_tracks_reso_p[trackSource][parIdx]->Fill(truepmom,(pmom-truepmom)/truepmom,trueeta);
    h_tracks_resoEta_pT[trackSource]->Fill(truept,(receta-trueeta),trueeta);
    h_tracks_resoPhi_pT[trackSource]->Fill(truept,(recphi-truephi),trueeta);

    // dca distributions
    Double_t dca2Dmum = _track_dca2D[itrk]*1e4;
//     std::cout << dca2Dmum << std::endl;
    if (h_tracks_dca2d_pT[trackSource]) h_tracks_dca2d_pT[trackSource]->Fill(truept,dca2Dmum,trueeta);
    if (h_tracks_dca2d_p[trackSource]) h_tracks_dca2d_p[trackSource]->Fill(truepmom,dca2Dmum,trueeta);
    if (h_tracks_ndf_p[trackSource]) h_tracks_ndf_p[trackSource]->Fill(truepmom,_track_ndf[itrk],trueeta);
    if (h_tracks_chi2_p[trackSource]) h_tracks_chi2_p[trackSource]->Fill(truepmom,_track_chi2[itrk],trueeta);

    // all tracks "N"
    h_tracksTrue_Eta_pT[trackSource][parIdx]->Fill(truept, trueeta);
    h_tracksRec_Eta_pT[trackSource][parIdx]->Fill(pt, receta);
    h_tracksTrue_Eta_p[trackSource][parIdx]->Fill(truepmom, trueeta);
    h_tracksRec_Eta_p[trackSource][parIdx]->Fill(pmom, receta);
  }
}

// **********************************************************************************************
// ************************ compare tracking from different sources *****************************
// **********************************************************************************************
void trackingcomparison() {
  if (verbosityTRKEFF > 3) std::cout << "start comparison" << std::endl;
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
    if (iOuterTrkSource >= (int)nTrackSources) continue;
    if (verbosityTRKEFF > 3)std::cout << "outer track: "<< iOuterTrk << std::endl;
    if ((int)_track_trueID[iOuterTrk] < 0) continue;
    if (_mcpart_RecTrackIDs[(int)_track_trueID[iOuterTrk]].size() > 1){
      for (int rtr = 0; rtr < (int)_mcpart_RecTrackIDs[(int)_track_trueID[iOuterTrk]].size(); rtr++){
        int iInnerTrk = _mcpart_RecTrackIDs[(int)_track_trueID[iOuterTrk]].at(rtr);
        int iInnerTrkSource = _track_source[iInnerTrk];
        if (iInnerTrkSource >= (int)nTrackSources) continue;
        if ( _track_source[iOuterTrk] >= iInnerTrkSource ) continue;

        TVector3 outerVec(_track_px[iOuterTrk], _track_py[iOuterTrk], _track_pz[iOuterTrk]);
        TVector3 innerVec(_track_px[iInnerTrk], _track_py[iInnerTrk], _track_pz[iInnerTrk]);
        if (verbosityTRKEFF > 3) std::cout <<"O source: " << iOuterTrkSource<< " I: " <<  iInnerTrk<< " source I: " << iInnerTrkSource << " eta diff: " << outerVec.Eta() - innerVec.Eta() << " pt diff " << outerVec.Pt() - innerVec.Pt() << std::endl;
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
//     if(h_chargedpart_MC_pT && h_chargedpart_rec_pT[i] && h_chargedpart_rec_truepT[i]) {
//       std::string name = "h_chargedpart_eff_pT_" + std::to_string(i);
//       TH2F*  h_chargedpart_eff_pT = (TH2F*)h_chargedpart_rec_pT[i]->Clone(name.c_str());
//       h_chargedpart_eff_pT->Divide(h_chargedpart_MC_pT);
//       h_chargedpart_eff_pT->Write(name.c_str(), TObject::kOverwrite);
//       name = "h_chargedpart_eff_truepT_" + std::to_string(i);
//       TH2F*  h_chargedpart_eff_truepT = (TH2F*)h_chargedpart_rec_truepT[i]->Clone(name.c_str());
//       h_chargedpart_eff_truepT->Divide(h_chargedpart_MC_pT);
//       h_chargedpart_eff_truepT->Write(name.c_str(), TObject::kOverwrite);
//     }
  }
  if(h_chargedpart_MC_p) h_chargedpart_MC_p->Write();
  for (unsigned int i = 0; i < nTrackSources; i++) {
    if(h_chargedpart_rec_p[i]) h_chargedpart_rec_p[i]->Write();
    if(h_chargedpart_rec_truep[i]) h_chargedpart_rec_truep[i]->Write();
    if(h_chargedpart_rec_trueE[i]) h_chargedpart_rec_trueE[i]->Write();
//     if(h_chargedpart_MC_p && h_chargedpart_rec_p[i] && h_chargedpart_rec_truep[i]) {
//       std::string name = "h_chargedpart_eff_p_" + std::to_string(i);
//       TH2F*  h_chargedpart_eff_p = (TH2F*)h_chargedpart_rec_p[i]->Clone(name.c_str());
//       h_chargedpart_eff_p->Divide(h_chargedpart_MC_p);
//       h_chargedpart_eff_p->Write(name.c_str(),TObject::kOverwrite);
//       name = "h_chargedpart_eff_truep_" + std::to_string(i);
//       TH2F*  h_chargedpart_eff_truep = (TH2F*)h_chargedpart_rec_truep[i]->Clone(name.c_str());
//       h_chargedpart_eff_truep->Divide(h_chargedpart_MC_p);
//       h_chargedpart_eff_truep->Write(name.c_str(),TObject::kOverwrite);
//     }
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
//       if(h_particle_MC_pT[ipart] && h_particle_rec_pT[i][ipart] && h_particle_rec_truepT[ipart]) {
//         std::string name = std::string("h_") + str_TRKEFF_mcparticles[ipart].Data() + "_eff_pT" + std::to_string(i);
//         h_particle_eff_pT[i][ipart] = (TH2F*)h_particle_rec_pT[i][ipart]->Clone(name.c_str());
//         h_particle_eff_pT[i][ipart]->Divide(h_particle_MC_pT[ipart]);
//         h_particle_eff_pT[i][ipart]->Write(name.c_str(), TObject::kOverwrite);
//         name = std::string("h_") + str_TRKEFF_mcparticles[ipart].Data() + "eff_truepT" + std::to_string(i);
//         h_particle_eff_truepT[i][ipart] = (TH2F*)h_particle_rec_truepT[i][ipart]->Clone(name.c_str());
//         h_particle_eff_truepT[i][ipart]->Divide(h_particle_MC_pT[ipart]);
//         h_particle_eff_truepT[i][ipart]->Write(name.c_str(), TObject::kOverwrite);
//       }
    }
    if(h_particle_MC_p[ipart]) h_particle_MC_p[ipart]->Write();
    for (unsigned int i = 0; i < nTrackSources; i++) {
      if(h_particle_rec_p[i][ipart]) h_particle_rec_p[i][ipart]->Write();
      if(h_particle_rec_truep[i][ipart]) h_particle_rec_truep[i][ipart]->Write();
      if(h_particle_rec_trueE[i][ipart]) h_particle_rec_trueE[i][ipart]->Write();
//       if(h_particle_MC_p[ipart] && h_particle_rec_p[i][ipart] && h_particle_rec_truep[i][ipart]) {
//         std::string name = std::string("h_") + str_TRKEFF_mcparticles[ipart].Data() + "_eff_p" + std::to_string(i);
//         h_particle_eff_p[i][ipart] = (TH2F*)h_particle_rec_p[i][ipart]->Clone(name.c_str());
//         h_particle_eff_p[i][ipart]->Divide(h_particle_MC_p[ipart]);
//         h_particle_eff_p[i][ipart]->Write(name.c_str(), TObject::kOverwrite);
//         name = std::string("h_") + str_TRKEFF_mcparticles[ipart].Data() + "_eff_truep" + std::to_string(i);
//         h_particle_eff_truep[i][ipart] = (TH2F*)h_particle_rec_truep[i][ipart]->Clone(name.c_str());
//         h_particle_eff_truep[i][ipart]->Divide(h_particle_MC_p[ipart]);
//         h_particle_eff_truep[i][ipart]->Write(name.c_str(), TObject::kOverwrite);
//       }
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
      if (id == 1){
        if(h_tracksTrue_Eta_pT[iTrkSource][nPart_TRKEFF]) h_tracksTrue_Eta_pT[iTrkSource][nPart_TRKEFF]->Sumw2();
        if(h_tracksRec_Eta_pT[iTrkSource][nPart_TRKEFF]) h_tracksRec_Eta_pT[iTrkSource][nPart_TRKEFF]->Sumw2();
        if(h_tracksTrue_Eta_p[iTrkSource][nPart_TRKEFF]) h_tracksTrue_Eta_p[iTrkSource][nPart_TRKEFF]->Sumw2();
        if(h_tracksRec_Eta_p[iTrkSource][nPart_TRKEFF]) h_tracksRec_Eta_p[iTrkSource][nPart_TRKEFF]->Sumw2();
      }
//         if (k == 0) std::cout << id << "\t" << h_tracksTrue_Eta_pT[nPart_TRKEFF]->GetEntries() << "\t" << h_tracksTrue_Eta_pT[id]->GetEntries() << std::endl;
      if(h_tracksTrue_Eta_pT[iTrkSource][nPart_TRKEFF]) h_tracksTrue_Eta_pT[iTrkSource][nPart_TRKEFF]->Add(h_tracksTrue_Eta_pT[iTrkSource][id]);
      if(h_tracksRec_Eta_pT[iTrkSource][nPart_TRKEFF]) h_tracksRec_Eta_pT[iTrkSource][nPart_TRKEFF]->Add(h_tracksRec_Eta_pT[iTrkSource][id]);
      if(h_tracksTrue_Eta_p[iTrkSource][nPart_TRKEFF]) h_tracksTrue_Eta_p[iTrkSource][nPart_TRKEFF]->Add(h_tracksTrue_Eta_p[iTrkSource][id]);
      if(h_tracksRec_Eta_p[iTrkSource][nPart_TRKEFF]) h_tracksRec_Eta_p[iTrkSource][nPart_TRKEFF]->Add(h_tracksRec_Eta_p[iTrkSource][id]);
      if (id == 0 ){
        if(h_tracks_reso_pT[iTrkSource][nPart_TRKEFF]) h_tracks_reso_pT[iTrkSource][nPart_TRKEFF]->Sumw2();
        if(h_tracks_reso_p[iTrkSource][nPart_TRKEFF]) h_tracks_reso_p[iTrkSource][nPart_TRKEFF]->Sumw2();
      }
      if(h_tracks_reso_pT[iTrkSource][id]){
        if(h_tracks_reso_pT[iTrkSource][nPart_TRKEFF]) h_tracks_reso_pT[iTrkSource][nPart_TRKEFF]->Add(h_tracks_reso_pT[iTrkSource][id]);
      }
      if(h_tracks_reso_p[iTrkSource][id]){
        if(h_tracks_reso_p[iTrkSource][nPart_TRKEFF]) h_tracks_reso_p[iTrkSource][nPart_TRKEFF]->Add(h_tracks_reso_p[iTrkSource][id]);
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
    for (Int_t id = 0; id < nPart_TRKEFF+1 ; id++) {
      if(h_tracks_reso_pT[iTrkSource][id]) h_tracks_reso_pT[iTrkSource][id]->Write();
      if(h_tracks_reso_p[iTrkSource][id]) h_tracks_reso_p[iTrkSource][id]->Write();
    }
    if (h_tracks_resoEta_pT[iTrkSource]) h_tracks_resoEta_pT[iTrkSource]->Write();
    if (h_tracks_resoPhi_pT[iTrkSource]) h_tracks_resoPhi_pT[iTrkSource]->Write();
    if (h_tracks_dca2d_pT[iTrkSource]) h_tracks_dca2d_pT[iTrkSource]->Write();
    if (h_tracks_dca2d_p[iTrkSource]) h_tracks_dca2d_p[iTrkSource]->Write();
    if (h_tracks_ndf_p[iTrkSource]) h_tracks_ndf_p[iTrkSource]->Write();
    if (h_tracks_chi2_p[iTrkSource]) h_tracks_chi2_p[iTrkSource]->Write();

    for (Int_t id = 0; id < nPart_TRKEFF+1; id++){
      if(h_tracksTrue_Eta_pT[iTrkSource][id]) h_tracksTrue_Eta_pT[iTrkSource][id]->Write();
      if(h_tracksRec_Eta_pT[iTrkSource][id]) h_tracksRec_Eta_pT[iTrkSource][id]->Write();
      if(h_tracksTrue_Eta_p[iTrkSource][id]) h_tracksTrue_Eta_p[iTrkSource][id]->Write();
      if(h_tracksRec_Eta_p[iTrkSource][id]) h_tracksRec_Eta_p[iTrkSource][id]->Write();
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
