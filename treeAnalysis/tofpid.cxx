// TOF pid analysis module

// steering 
int verbosityTOF = 0;

//----------------------------------------------------------------------------
// define histos
//----------------------------------------------------------------------------
bool bHistTOFCreated    = false;
TH2F* h_TPH_InvGenBeta_Eta_p[4][nPID][nEta+1]             = {{{NULL}}};
TH2F* h_TPH_InvBeta_Eta_p[4][nPID][nEta+1]                = {{{NULL}}};
TH2F* h_TPH_InvSmearBeta_Eta_p[4][nPID][nEta+1]           = {{{NULL}}};
TH2F* h_TPH_InvSmearAndT0Beta_Eta_p[4][nPID][nEta+1]      = {{{NULL}}};
TH2F* h_TPH_Res_InvBeta_Eta_p[4][nPID][nEta+1]            = {{{NULL}}};
TH2F* h_TPH_Res_InvSmearBeta_Eta_p[4][nPID][nEta+1]       = {{{NULL}}};
TH2F* h_TPH_Res_InvSmearAndT0Beta_Eta_p[4][nPID][nEta+1]  = {{{NULL}}};

// const Double_t mPtCutForT0[3] = {0.6, 2, 4};
const Double_t mPtCutForT0[3] = {6, 6, 6};
const Double_t tDiffAbort = 1.e-6;
const   Int_t nSces       = 2;
TString sceName[nSces]    = {"wScatE", "woScatE"};
TString nameLayerTTL[3]   = {"ETTL", "CTTL", "FTTL"};

TH2D *h_TPH_InitT0VsNtrk[nSces];
TH2D *h_TPH_InitT0VsNttlhit[nSces];
TH2D *h_TPH_NttlhitVsNtrk_wInitT0[nSces];
TH2D *h_TPH_NttlhitVsNtrk_wFinalT0[nSces];
TH2D *h_TPH_T0VsNtrk[nSces];
TH2D *h_TPH_T0VsNttlhit[nSces];
TH2D *h_TPH_InitDbetaVsP_pi[nSces];
TH2D *h_TPH_InitDbetaVsP_k[nSces];
TH2D *h_TPH_InitDbetaVsP_p[nSces];
TH2D *h_TPH_DbetaVsP_pi[nSces];
TH2D *h_TPH_DbetaVsP_k[nSces];
TH2D *h_TPH_DbetaVsP_p[nSces];
TH2D *h_TPH_BetaVsP_T0[3]                           = {NULL};
TH2D *h_TPH_DbetaVsP_wT0[nPID][3]                   = {NULL};
TH2D *h_TPH_DbetaVsP_woT0[nPID][3]                  = {NULL};
TH2D *h_TPH_DgenbetaVsGenP[nPID][3]                 = {NULL};
TH2D *h_TPH_DbetaVsP_woT0_test[nPID][3]             = {NULL};

// define c
const double c                  = 29.9792458; // cm/ns

const Int_t nIterT = 4; //4
TH2D *h_TPH_T0VsIter;
TH2D *h_TPH_UPartIter;


//___________________________________________________________________________
Int_t GetParIndexFromPDG(Int_t pdg){
  Int_t parIdx = -1;
  if(TMath::Abs(pdg) == 11)   parIdx = 1;
  if(TMath::Abs(pdg) == 13)   parIdx = 2;
  if(TMath::Abs(pdg) == 211)  parIdx = 3;
  if(TMath::Abs(pdg) == 321)  parIdx = 4;
  if(TMath::Abs(pdg) == 2212) parIdx = 5;
  return parIdx;
}  


//----------------------------------------------------------------------------
// main routine for TOF pid extractions
//----------------------------------------------------------------------------
void tofpidhistos(){

  //----------------------------------------------------------------------------
  // create histos if not yet available
  //----------------------------------------------------------------------------
  if (!bHistTOFCreated){
    testBinning(binningPFine, nBinsPFine);
    
    h_TPH_T0VsIter      = new TH2D("h_T0VsIter", "; Iter; T_{0} (ns)", nSces*nIterT, -0.5, nSces*nIterT-0.5, 500, -0.5, 0.5);
    h_TPH_UPartIter     = new TH2D("h_TPH_UsedParticlesVsIter", "; Iter; N_{particles}", nSces*nIterT, -0.5, nSces*nIterT-0.5, 50, 0, 50);
    for(int sc=0; sc<nSces; sc++){
      h_TPH_InitT0VsNtrk[sc]    = new TH2D(Form("h_InitT0VsNtrk_%s", sceName[sc].Data()), "; N_{trk}; Initial T_{0} (ns)", 100, 0, 100, 500, -0.5, 0.5);
      h_TPH_InitT0VsNttlhit[sc] = new TH2D(Form("h_InitT0VsNttlhit_%s", sceName[sc].Data()), ";  N_{trk w/ TTL}; Initial T_{0} (ns)", 100, 0, 100, 500, -0.5, 0.5);
      
      
      
      h_TPH_NttlhitVsNtrk_wInitT0[sc]  = new TH2D(Form("h_NttlhitVsNtrk_wInitT0_%s", sceName[sc].Data()), "; N_{trk}; N_{trk w/ TTL}", 100, 0, 100, 100, 0, 100);
      h_TPH_NttlhitVsNtrk_wFinalT0[sc] = new TH2D(Form("h_NttlhitVsNtrk_wFinalT0_%s", sceName[sc].Data()), "; N_{trk}; N_{trk w/ TTL}", 100, 0, 100, 100, 0, 100);

      h_TPH_T0VsNtrk[sc]          = new TH2D(Form("h_T0VsNtrk_%s", sceName[sc].Data()), "; N_{trk}; Final T_{0} (ns)", 100, 0, 100, 500, -0.5, 0.5);
      h_TPH_T0VsNttlhit[sc]       = new TH2D(Form("h_T0VsNttlhit_%s", sceName[sc].Data()), "; N_{trk w/ TTL}; Final T_{0} (ns)", 100, 0, 100, 500, -0.5, 0.5);

      h_TPH_InitDbetaVsP_pi[sc]   = new TH2D(Form("h_InitDbetaVsP_pi_%s", sceName[sc].Data()), "#pi; p (GeV); 1/#beta - 1/#beta_{#pi}", 300, 0, 15, 500, -0.5, 0.5);
      h_TPH_InitDbetaVsP_k[sc]    = new TH2D(Form("h_InitDbetaVsP_k_%s", sceName[sc].Data()), "K; p (GeV); 1/#beta - 1/#beta_{k}", 300, 0, 15, 500, -0.5, 0.5);
      h_TPH_InitDbetaVsP_p[sc]    = new TH2D(Form("h_InitDbetaVsP_p_%s", sceName[sc].Data()), "P; p (GeV); 1/#beta - 1/#beta_{p}", 300, 0, 15, 500, -0.5, 0.5);

      h_TPH_DbetaVsP_pi[sc]       = new TH2D(Form("h_DbetaVsP_pi_%s", sceName[sc].Data()), "#pi; p (GeV); 1/#beta - 1/#beta_{#pi}", 300, 0, 15, 500, -0.5, 0.5);
      h_TPH_DbetaVsP_k[sc]        = new TH2D(Form("h_DbetaVsP_k_%s", sceName[sc].Data()), "k; p (GeV); 1/#beta - 1/#beta_{k}", 300, 0, 15, 500, -0.5, 0.5);
      h_TPH_DbetaVsP_p[sc]        = new TH2D(Form("h_DbetaVsP_p_%s", sceName[sc].Data()), "p; p (GeV); 1/#beta - 1/#beta_{p}", 300, 0, 15, 500, -0.5, 0.5);
    }
    
    for (Int_t dir = 0; dir < 3; dir++){
      h_TPH_BetaVsP_T0[dir] = new TH2D(Form("h_BetaVsP_T0_%s", nameLayerTTL[dir].Data()), Form ("%s; p (GeV); 1/#beta", nameLayerTTL[dir].Data()), nBinsPFine, binningPFine, 600, 0.9, 1.5);
      for (Int_t id = 0; id < nPID; id++){
        h_TPH_DbetaVsP_wT0[id][dir]    = new TH2D(Form("h_DbetaVsP_wT0_%s_%s", nameLayerTTL[dir].Data(), partName[id].Data()), 
                                                    Form("%s %s; p_{MC} (GeV); 1/#beta - 1/#beta_{ident}", nameLayerTTL[dir].Data(), partName[id].Data()),  nBinsPFine, binningPFine, 500, -0.1, 0.1);
        h_TPH_DbetaVsP_woT0[id][dir]   = new TH2D(Form("h_DbetaVsP_woT0_%s_%s", nameLayerTTL[dir].Data(), partName[id].Data()), 
                                                    Form("%s %s; p_{MC} (GeV); 1/#beta - 1/#beta_{ident}", nameLayerTTL[dir].Data(), partName[id].Data()),  nBinsPFine, binningPFine, 500, -0.1, 0.1);
        h_TPH_DgenbetaVsGenP[id][dir]  = new TH2D(Form("h_DgenbetaVsGenP_%s_%s", nameLayerTTL[dir].Data(), partName[id].Data()), 
                                                    Form("%s %s; p_{MC} (GeV); 1/#beta - 1/#beta_{ident}", nameLayerTTL[dir].Data(), partName[id].Data()),  nBinsPFine, binningPFine, 600, -0.3, 0.3);
        h_TPH_DbetaVsP_woT0_test[id][dir]  = new TH2D(Form("h_DbetaVsP_woT0_test_%s_%s", nameLayerTTL[dir].Data(), partName[id].Data()), 
                                                    Form("%s %s; p_{MC} (GeV); 1/#beta - 1/#beta_{ident}", nameLayerTTL[dir].Data(), partName[id].Data()),  nBinsPFine, binningPFine, 500, -0.1, 0.1);
      }
    }
    
    for (Int_t et = 0; et<nEta+1; et++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (et < nEta){
        etaMin = partEta[et];
        etaMax = partEta[et+1];
      }      
      
      
      for (Int_t e = 0; e< 4; e++){
        for (Int_t id = 0; id < nPID; id++){
          h_TPH_InvGenBeta_Eta_p[e][id][et]      = new TH2F(Form("h_InvGenBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(), partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{MC}", etaMin, etaMax), 
                                                  nBinsPFine, binningPFine, 300, 0.9, 1.5);

          h_TPH_InvBeta_Eta_p[e][id][et]         = new TH2F(Form("h_InvBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(), partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{rec}", etaMin, etaMax), 
                                                  nBinsPFine, binningPFine, 300, 0.9, 1.5);
          h_TPH_InvSmearBeta_Eta_p[e][id][et]    = new TH2F(Form("h_InvSmearBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{smear}", etaMin, etaMax), 
                                                  nBinsPFine, binningPFine, 300, 0.9, 1.5);
          h_TPH_InvSmearAndT0Beta_Eta_p[e][id][et]    = new TH2F(Form("h_InvSmearAndT0Beta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{smear}", etaMin, etaMax), 
                                                  nBinsPFine, binningPFine, 300, 0.9, 1.5);
          h_TPH_Res_InvBeta_Eta_p[e][id][et]     = new TH2F(Form("h_Res_InvBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (1/#beta_{rec} - 1/#beta_{MC})/(1/#beta_{MC}) ", etaMin, etaMax), 
                                                  nBinsPFine, binningPFine, 500, -0.1, 0.1);
          h_TPH_Res_InvSmearBeta_Eta_p[e][id][et]= new TH2F(Form("h_Res_InvSmearBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (1/#beta_{smear} - 1/#beta_{MC})/(1/#beta_{MC}) ", etaMin, etaMax), 
                                                  nBinsPFine, binningPFine, 500, -0.1, 0.1);
          h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][id][et]= new TH2F(Form("h_Res_InvSmearAndT0Beta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (1/#beta_{smear} - 1/#beta_{MC})/(1/#beta_{MC}) ", etaMin, etaMax), 
                                                  nBinsPFine, binningPFine, 500, -0.1, 0.1);
        }
      }
    }
    bHistTOFCreated = true;
  }
  float mT0[nIterT];
  for(int it=0; it<nIterT; it++) mT0[it] = -9999;
  
  // fixed mass assumptions
  float mElectron = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  float mPion     = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  float mProton   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  float mKaon     = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  //----------------------------------------------------------------------------
  // determine start time based on scattered electron
  //----------------------------------------------------------------------------
  int nScatElect  = 0;            // number of scattered electron candidates
  int scEIdx      = -1;           // what kind of determination
  Int_t nTracksWTime = 0;         // number of track with timing
  Int_t nTracksF     = 0;         // number of track with timing
  for(Int_t itrk=0; itrk<(Int_t)_nTracks; itrk++){
    if (_track_source[itrk] != 0) continue;
    if (_track_trueID[itrk] < 0) continue;
    nTracksF++;
    if (!_track_hasTTL[itrk]) continue;
    // count tracks with timing
    nTracksWTime++; // count tracks with timing
  
    // abort if not a probable scattered electron
    if (!_track_TOFmeas[itrk].at(0).isProbEScat) continue;
    
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
    float trueP         = mcpartvec.Mag();
    // generated beta for electrons
    float betaGen = trueP/sqrt(trueP*trueP + mElectron*mElectron);
    // calculate average start time for single electron track
    float tstart  = 0;
    for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
      float tflight  = _track_TOFmeas[itrk].at(tl).pathlength/(betaGen*c);
      tstart  += _track_TOFmeas[itrk].at(tl).recTime-tflight;
    }
    tstart  = tstart/_track_TOFmeas[itrk].size();
    if (nScatElect == 0)  mT0[0]  = tstart;
    else                  mT0[0]  += tstart;
    nScatElect++;
  }
//   cout << __LINE__ << endl;
  //----------------------------------------------------------------------------
  // Did we find the scattered electron?
  //----------------------------------------------------------------------------
  if (nScatElect > 0){ 
    //----------------------------------------------------------------------------
    // YES! average start time if its more than 1
    //----------------------------------------------------------------------------
    scEIdx  = 0;
    mT0[0]  = mT0[0]/nScatElect;
    h_TPH_UPartIter->Fill(0., nScatElect);
  } else {
    //----------------------------------------------------------------------------
    // NO! assume all particles are pions and determine starttime from that
    //----------------------------------------------------------------------------
    for(Int_t itrk=0; itrk<(Int_t)_nTracks; itrk++){
      if (_track_source[itrk] != 0) continue;
      if (_track_trueID[itrk] < 0) continue;
      if (! _track_hasTTL[itrk]) continue;
      TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
      TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
      float recP          = recpartvec.Mag();
      float trueP         = mcpartvec.Mag();
      // generted beta for pion
      float betaGen = trueP/sqrt(trueP*trueP + mPion* mPion);
      // calculate average start time for single electron track
      float tstart  = 0;
      for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
        float tflight  = _track_TOFmeas[itrk].at(tl).pathlength/(betaGen*c);
        tstart  += _track_TOFmeas[itrk].at(tl).recTime-tflight;
      }
      tstart  = tstart/_track_TOFmeas[itrk].size();
      if (itrk == 0)  mT0[0]  = tstart;
      else            mT0[0]  += tstart;
    } 
    // make average over all particles
    if (nTracksWTime > 1){
      scEIdx = 1;
      mT0[0] = mT0[0] / nTracksWTime;
      h_TPH_UPartIter->Fill(nIterT, nTracksWTime);
    }
  }

//   cout << __LINE__ << endl;
  // fill histos for 0th iteration
  if(scEIdx > -1){
    h_TPH_InitT0VsNtrk[scEIdx]->Fill(nTracksF, mT0[0]);
    h_TPH_InitT0VsNttlhit[scEIdx]->Fill(nTracksWTime, mT0[0]);
  }
  //----------------------------------------------------------------------------
  // try to improve start time reso with iterative PID
  //----------------------------------------------------------------------------
  for(int it=0; it<nIterT; it++){
    // abort if previous iteration was good enough
    if(fabs(mT0[it-1] + 9999) < tDiffAbort) break;
    Int_t  nUsedTTLHits = 0;
    for(Int_t itrk=0; itrk<(Int_t)_nTracks; itrk++){
      if (_track_source[itrk] != 0) continue;
      if (_track_trueID[itrk] < 0) continue;
      if (! _track_hasTTL[itrk]) continue;
      TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
      TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
      float recP          = recpartvec.Mag();
      float trueP         = mcpartvec.Mag();
      int pid             = _mcpart_PDG[(int)_track_trueID[itrk]];
      Double_t mass       = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
    
      Bool_t isPID = kFALSE;

      float InvBetaWT0  = 0;
      float InvBetaW0T0 = 0;
      float InvBetaPath = 0;
      for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
        InvBetaPath += c*(_track_TOFmeas[itrk].at(tl).genTime)/_track_TOFmeas[itrk].at(tl).pathlength; 
        InvBetaWT0  += c*(_track_TOFmeas[itrk].at(tl).recTime -mT0[it-1] )/_track_TOFmeas[itrk].at(tl).pathlength; 
        InvBetaW0T0 += c*(_track_TOFmeas[itrk].at(tl).recTime)/_track_TOFmeas[itrk].at(tl).pathlength; 
      }
      InvBetaPath = InvBetaPath/_track_nTTL[itrk];
      InvBetaWT0  = InvBetaWT0/_track_nTTL[itrk];
      InvBetaW0T0 = InvBetaW0T0/_track_nTTL[itrk];

      float beta_e  = recP / sqrt(recP*recP + mElectron*mElectron);
      float beta_pi = recP / sqrt(recP*recP + mPion*mPion);
      float beta_k  = recP / sqrt(recP*recP + mKaon*mKaon);
      float beta_p  = recP / sqrt(recP*recP + mProton*mProton);
      float beta_id  = recP / sqrt(recP*recP + mass*mass);
      
      float dBeta_pi = InvBetaWT0 - 1/beta_pi;
      float dBeta_k  = InvBetaWT0 - 1/beta_k;
      float dBeta_p  = InvBetaWT0 - 1/beta_p;

      float genBeta_pi = trueP / sqrt(trueP*trueP + mPion*mPion);
      float genBeta_k  = trueP / sqrt(trueP*trueP + mKaon*mKaon);
      float genBeta_p  = trueP / sqrt(trueP*trueP + mProton*mProton);
      float genBeta_id  = trueP / sqrt(trueP*trueP + mass*mass);
      
      Int_t dir = GetRegionFromEta(mcpartvec.Eta());
      Int_t parIdx = GetParIndexFromPDG(pid);
      
      // show the starting point for the T0 improvements
      if(!_track_TOFmeas[itrk].at(0).isProbEScat && it==1){
        h_TPH_InitDbetaVsP_pi[scEIdx]->Fill(trueP, dBeta_pi);
        h_TPH_InitDbetaVsP_k[scEIdx]->Fill(trueP, dBeta_k);
        h_TPH_InitDbetaVsP_p[scEIdx]->Fill(trueP, dBeta_p);
        
        if (dir != -1 && parIdx != -1){
          h_TPH_DgenbetaVsGenP[parIdx][dir]->Fill(trueP, InvBetaPath - 1/genBeta_id);
          h_TPH_DbetaVsP_woT0_test[parIdx][dir]->Fill(trueP, InvBetaW0T0 - 1/beta_id);
        }
      }
      dBeta_pi  = dBeta_pi/(1/beta_pi);
      dBeta_k   = dBeta_k/(1/beta_k);
      dBeta_p   = dBeta_p/(1/beta_p);
      
      //----------------------------------------------------------------------------
      // use most probable PID for t0 determination
      //----------------------------------------------------------------------------
      float tstart  = 0;
      // scattered electrons
      if (_track_TOFmeas[itrk].at(0).isProbEScat){    
        isPID = kTRUE;  
        for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
          float tflight  = _track_TOFmeas[itrk].at(tl).pathlength/(beta_e*c);
          tstart  += _track_TOFmeas[itrk].at(tl).recTime-tflight;
        }
        tstart  = tstart/_track_TOFmeas[itrk].size();
        nUsedTTLHits++;
      // pions
      } else if (recP < mPtCutForT0[0] && dBeta_pi < 0.01 && dBeta_pi > -0.01){
        isPID = kTRUE;
        for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
          float tflight  = _track_TOFmeas[itrk].at(tl).pathlength/(beta_pi*c);
          tstart  += _track_TOFmeas[itrk].at(tl).recTime-tflight;
        }
        tstart  = tstart/_track_TOFmeas[itrk].size();
        nUsedTTLHits++;
      // kaons  
      } else if ( recP < mPtCutForT0[1] && fabs(dBeta_k) < 0.01 ){
        isPID = kTRUE;
        for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
          float tflight  = _track_TOFmeas[itrk].at(tl).pathlength/(beta_k*c);
          tstart  += _track_TOFmeas[itrk].at(tl).recTime-tflight;
        }
        tstart  = tstart/_track_TOFmeas[itrk].size();
        nUsedTTLHits++;
      // protons  
      } else if (recP < mPtCutForT0[2] && fabs(dBeta_p) < 0.01){ 
        isPID = kTRUE;
        for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
          float tflight  = _track_TOFmeas[itrk].at(tl).pathlength/(beta_p*c);
          tstart  += _track_TOFmeas[itrk].at(tl).recTime-tflight;
        }
        tstart  = tstart/_track_TOFmeas[itrk].size();
        nUsedTTLHits++;      
      } else if ( recP > 6){
//       } else if ( recP > 15){
        isPID = kTRUE;  
        for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
          float tflight  = _track_TOFmeas[itrk].at(tl).pathlength/(beta_e*c);
          tstart  += _track_TOFmeas[itrk].at(tl).recTime-tflight;
        }
        tstart  = tstart/_track_TOFmeas[itrk].size();
        nUsedTTLHits++;     
      }
      // add up start times
      if(isPID){
        if(fabs(mT0[it] + 9999) < tDiffAbort) mT0[it] = tstart;
        else mT0[it] += tstart;
      }
    }
    // calculate new average over event
    if(nUsedTTLHits > 0){
      mT0[it] = mT0[it]/nUsedTTLHits;
      h_TPH_UPartIter->Fill(scEIdx*nIterT+it, nUsedTTLHits);
    }
  }
//   cout << __LINE__ << endl;
  //----------------------------------------------------------------------------
  // determine final T0
  //----------------------------------------------------------------------------
  float mFinalT0 = -9999;
  // for the events with the scattering electron detected, the initial T0 could be used as the final T0
  if(nScatElect > 0){
    for(Int_t it=0; it<nIterT; it++){ 
      if(mT0[it] > -9999) mFinalT0 = mT0[it];
    }
  } else {
    if (nIterT == 1){
      if(mT0[0] > -9999) mFinalT0 = mT0[0];
    } else {
      for(Int_t it=1; it<nIterT; it++){
        if(mT0[it] > -9999) mFinalT0 = mT0[it];
      }
    }
  }
  for (Int_t it=0; it<nIterT; it++){
    if(mT0[it] > -9999) h_TPH_T0VsIter->Fill(scEIdx*nIterT+it, mT0[it]);
  }
  // fill controll hists for T0 vs track numbers 
  if(mT0[0] > -9999){
    h_TPH_NttlhitVsNtrk_wInitT0[scEIdx]->Fill(nTracksF, nTracksWTime);
    if(mFinalT0 > -9999) h_TPH_NttlhitVsNtrk_wFinalT0[scEIdx]->Fill(nTracksF, nTracksWTime);
  }
  if(mFinalT0 > -9999){
    h_TPH_T0VsNtrk[scEIdx]->Fill(nTracksF, mFinalT0);
    h_TPH_T0VsNttlhit[scEIdx]->Fill(nTracksWTime, mFinalT0);
  }
  
//   cout << __LINE__ << endl;
  //----------------------------------------------------------------------------
  // run over tracks
  //----------------------------------------------------------------------------
  for(Int_t itrk=0; itrk<(Int_t)_nTracks; itrk++){
    if (_track_source[itrk] != 0) continue;
    if (_track_trueID[itrk] < 0) continue;
    if (! _track_hasTTL[itrk]) continue;
    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
    // cout << itrk << endl;
    float recEta        = recpartvec.Eta();
    float trueEta       = mcpartvec.Eta();
    float recP          = recpartvec.Mag();
    float trueP         = mcpartvec.Mag();
    int pid             = _mcpart_PDG[(int)_track_trueID[itrk]];
    Double_t mass       = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
    Float_t beta[4]     = {0.};       // 1/beta: 
                                      //  0 - generated
                                      //  1 - path length included
                                      //  2 - path length included + smeared time
                                      //  3 - path length included + smeared time + start time
    
    Int_t dir     = GetRegionFromEta(mcpartvec.Eta());
    Int_t parIdx  = GetParIndexFromPDG(pid);
    if (parIdx == -1) continue;
    if (verbosityTOF > 2) cout << "\tTRKTRUEID " << (int)_track_trueID[itrk] << "\tPDG: " << pid << "\tpT: " << recP << "\tEta " << recEta << endl;
    
    if (verbosityTOF){
      std::cout << "Track p " << _track_px[itrk] << "\t" << _track_py[itrk] << "\t" << _track_pz[itrk] << "\t eta " << recEta <<  std::endl;
      std::cout << "\t TTL: " << _track_hasTTL[itrk] << "\t n layer TTL " << _track_nTTL[itrk] << "\t Tracker: " << _track_hasIL[itrk] << "\t" << _track_hasOL[itrk] << "\t n Tr. layers " << _track_nTrL[itrk] << std::endl;
      for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
        std::cout << "\t\t" << tl << "\t" << tl << "\t pathlengh " << _track_TOFmeas[itrk].at(tl).pathlength << "\t gen time: "  << _track_TOFmeas[itrk].at(tl).genTime << "\t rec t \t"<<  _track_TOFmeas[itrk].at(tl).recTime << std::endl;
      }
    }
//     cout << __LINE__ << endl;
  
    //*************************************************************
    // calculate average 1/\beta for multiple measurements per track
    //*************************************************************
    for (int tl = 0; tl < int(_track_TOFmeas[itrk].size()); tl++){
      beta[1] += c*_track_TOFmeas[itrk].at(tl).genTime/_track_TOFmeas[itrk].at(tl).pathlength; 
      beta[2] += c*_track_TOFmeas[itrk].at(tl).recTime/_track_TOFmeas[itrk].at(tl).pathlength; 
      if (mFinalT0 > -9999) beta[3] += c*(_track_TOFmeas[itrk].at(tl).recTime-mFinalT0)/_track_TOFmeas[itrk].at(tl).pathlength; 
    }
    float betaRec = sqrt(recP*recP + mass*mass)/recP;
    beta[0] = sqrt(trueP*trueP + mass*mass)/trueP;
    beta[1] = beta[1]/_track_nTTL[itrk];
    beta[2] = beta[2]/_track_nTTL[itrk];
    if (mFinalT0 > -9999) beta[3] = beta[3]/_track_nTTL[itrk];
    
//     cout << __LINE__ << endl;
    // determine eta bin
    Int_t et = 0;
    while (partEta[et+1] < trueEta && et < nEta) et++;
//     cout << __LINE__ << endl;
    //*************************************************************
    // fill beta hists
    //*************************************************************
    if (mFinalT0 > -9999){
//       cout << __LINE__ << endl;
      if (dir != -1 ) h_TPH_BetaVsP_T0[dir]->Fill(trueP, beta[3]);
      float dBeta_pi = beta[3] - 1/(recP / sqrt(recP*recP + mPion*mPion));
      float dBeta_k  = beta[3] - 1/(recP / sqrt(recP*recP + mKaon*mKaon));
      float dBeta_p  = beta[3] - 1/(recP / sqrt(recP*recP + mProton*mProton));
      if (scEIdx > -1){
//         cout << __LINE__ << "\t" << scEIdx << "\t" << h_TPH_DbetaVsP_pi[scEIdx] << "\t" <<  h_TPH_DbetaVsP_k[scEIdx] <<"\t" << h_TPH_DbetaVsP_p[scEIdx]<< endl;
        h_TPH_DbetaVsP_pi[scEIdx]->Fill(trueP, dBeta_pi);
        h_TPH_DbetaVsP_k[scEIdx]->Fill(trueP, dBeta_k);
        h_TPH_DbetaVsP_p[scEIdx]->Fill(trueP, dBeta_p);
      }
    }
//     cout << __LINE__ << endl;
    h_TPH_InvGenBeta_Eta_p[0][parIdx][et]->Fill(trueP, beta[0]);
    h_TPH_InvBeta_Eta_p[0][parIdx][et]->Fill(trueP, beta[1]);
    h_TPH_InvSmearBeta_Eta_p[0][parIdx][et]->Fill(trueP, beta[2]);
    h_TPH_Res_InvBeta_Eta_p[0][parIdx][et]->Fill(trueP, (beta[1]-beta[0])/beta[0]);
    h_TPH_Res_InvSmearBeta_Eta_p[0][parIdx][et]->Fill(trueP, (beta[2]-beta[0])/beta[0]);
    if (mFinalT0 > -9999){
      h_TPH_InvSmearAndT0Beta_Eta_p[0][parIdx][et]->Fill(trueP, beta[3]);
      h_TPH_Res_InvSmearAndT0Beta_Eta_p[0][parIdx][et]->Fill(trueP, (beta[3]-beta[0])/beta[0]);
    }
//     cout << __LINE__ << endl;
    if ( _track_nTrL[itrk] < 3){
      h_TPH_InvGenBeta_Eta_p[1][parIdx][et]->Fill(trueP, beta[0]);
      h_TPH_InvBeta_Eta_p[1][parIdx][et]->Fill(trueP, beta[1]);
      h_TPH_InvSmearBeta_Eta_p[1][parIdx][et]->Fill(trueP, beta[2]);
      h_TPH_Res_InvBeta_Eta_p[1][parIdx][et]->Fill(trueP, (beta[1]-beta[0])/beta[0]);
      h_TPH_Res_InvSmearBeta_Eta_p[1][parIdx][et]->Fill(trueP, (beta[2]-beta[0])/beta[0]);
      if (mFinalT0 > -9999){
        h_TPH_InvSmearAndT0Beta_Eta_p[1][parIdx][et]->Fill(trueP, beta[3]);
        h_TPH_Res_InvSmearAndT0Beta_Eta_p[1][parIdx][et]->Fill(trueP, (beta[3]-beta[0])/beta[0]);
      } 
    }
//     cout << __LINE__ << endl;
    if ((_track_nTTL[itrk]+_track_nTrL[itrk]) >= 3){
      h_TPH_InvGenBeta_Eta_p[2][parIdx][et]->Fill(trueP, beta[0]);
      h_TPH_InvBeta_Eta_p[2][parIdx][et]->Fill(trueP, beta[1]);
      h_TPH_InvSmearBeta_Eta_p[2][parIdx][et]->Fill(trueP, beta[2]);
      h_TPH_Res_InvBeta_Eta_p[2][parIdx][et]->Fill(trueP, (beta[1]-beta[0])/beta[0]);
      h_TPH_Res_InvSmearBeta_Eta_p[2][parIdx][et]->Fill(trueP, (beta[2]-beta[0])/beta[0]);
      if (mFinalT0 > -9999){
        h_TPH_InvSmearAndT0Beta_Eta_p[2][parIdx][et]->Fill(trueP, beta[3]);
        h_TPH_Res_InvSmearAndT0Beta_Eta_p[2][parIdx][et]->Fill(trueP, (beta[3]-beta[0])/beta[0]);
      }
    }
    if ((_track_nTTL[itrk]+_track_nTrL[itrk]) >= 3 && _track_hasIL[itrk]){
      h_TPH_InvGenBeta_Eta_p[3][parIdx][et]->Fill(trueP, beta[0]);
      h_TPH_InvBeta_Eta_p[3][parIdx][et]->Fill(trueP, beta[1]);
      h_TPH_InvSmearBeta_Eta_p[3][parIdx][et]->Fill(trueP, beta[2]);
      h_TPH_Res_InvBeta_Eta_p[3][parIdx][et]->Fill(trueP, (beta[1]-beta[0])/beta[0]);
      h_TPH_Res_InvSmearBeta_Eta_p[3][parIdx][et]->Fill(trueP, (beta[2]-beta[0])/beta[0]);
      if (mFinalT0 > -9999){
        h_TPH_InvSmearAndT0Beta_Eta_p[2][parIdx][et]->Fill(trueP, beta[3]);
        h_TPH_Res_InvSmearAndT0Beta_Eta_p[2][parIdx][et]->Fill(trueP, (beta[3]-beta[0])/beta[0]);
      }
    }
//     cout << __LINE__ << endl;
    if (dir != -1){
      if (mFinalT0 > -9999) h_TPH_DbetaVsP_wT0[parIdx][dir]->Fill(trueP, beta[3]-betaRec);
      h_TPH_DbetaVsP_woT0[parIdx][dir]->Fill(trueP, beta[2]-betaRec);
    }
  }
//   cout << __LINE__ << endl;
}

//----------------------------------------------------------------------------
// save TOF PID hists
//----------------------------------------------------------------------------
void tofpidhistosSave(){
  
  // ***********************************************************************************************
  // ************************** sum up hists correctly *********************************************
  // ***********************************************************************************************
  for (Int_t id = 1; id < 6; id++){
    for (Int_t e = 0; e < 4; e++){
      for (Int_t et = 0; et < nEta; et++){
        if (et == 0 && id == 0){
          h_TPH_InvGenBeta_Eta_p[e][0][et]->Sumw2();
          h_TPH_InvGenBeta_Eta_p[e][0][nEta]->Sumw2();
          h_TPH_InvGenBeta_Eta_p[e][id][nEta]->Sumw2();
          h_TPH_InvBeta_Eta_p[e][0][et]->Sumw2();
          h_TPH_InvBeta_Eta_p[e][0][nEta]->Sumw2();
          h_TPH_InvBeta_Eta_p[e][id][nEta]->Sumw2();
          h_TPH_InvSmearBeta_Eta_p[e][0][et]->Sumw2();
          h_TPH_InvSmearBeta_Eta_p[e][0][nEta]->Sumw2();
          h_TPH_InvSmearBeta_Eta_p[e][id][nEta]->Sumw2();
          h_TPH_Res_InvBeta_Eta_p[e][0][et]->Sumw2();
          h_TPH_Res_InvBeta_Eta_p[e][0][nEta]->Sumw2();
          h_TPH_Res_InvBeta_Eta_p[e][id][nEta]->Sumw2();
          h_TPH_Res_InvSmearBeta_Eta_p[e][0][et]->Sumw2();
          h_TPH_Res_InvSmearBeta_Eta_p[e][0][nEta]->Sumw2();
          h_TPH_Res_InvSmearBeta_Eta_p[e][id][nEta]->Sumw2();
          h_TPH_InvSmearAndT0Beta_Eta_p[e][0][et]->Sumw2();
          h_TPH_InvSmearAndT0Beta_Eta_p[e][0][nEta]->Sumw2();
          h_TPH_InvSmearAndT0Beta_Eta_p[e][id][nEta]->Sumw2();
          h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][0][et]->Sumw2();
          h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][0][nEta]->Sumw2();
          h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][id][nEta]->Sumw2();
        }
        h_TPH_InvGenBeta_Eta_p[e][0][et]->Add(h_TPH_InvGenBeta_Eta_p[e][id][et]);
        h_TPH_InvGenBeta_Eta_p[e][0][nEta]->Add(h_TPH_InvGenBeta_Eta_p[e][id][et]);
        h_TPH_InvGenBeta_Eta_p[e][id][nEta]->Add(h_TPH_InvGenBeta_Eta_p[e][id][et]);
        h_TPH_InvBeta_Eta_p[e][0][et]->Add(h_TPH_InvBeta_Eta_p[e][id][et]);
        h_TPH_InvBeta_Eta_p[e][0][nEta]->Add(h_TPH_InvBeta_Eta_p[e][id][et]);
        h_TPH_InvBeta_Eta_p[e][id][nEta]->Add(h_TPH_InvBeta_Eta_p[e][id][et]);
        h_TPH_InvSmearBeta_Eta_p[e][0][et]->Add(h_TPH_InvSmearBeta_Eta_p[e][id][et]);
        h_TPH_InvSmearBeta_Eta_p[e][0][nEta]->Add(h_TPH_InvSmearBeta_Eta_p[e][id][et]);
        h_TPH_InvSmearBeta_Eta_p[e][id][nEta]->Add(h_TPH_InvSmearBeta_Eta_p[e][id][et]);
        h_TPH_Res_InvBeta_Eta_p[e][0][et]->Add(h_TPH_Res_InvBeta_Eta_p[e][id][et]);
        h_TPH_Res_InvBeta_Eta_p[e][0][nEta]->Add(h_TPH_Res_InvBeta_Eta_p[e][id][et]);
        h_TPH_Res_InvBeta_Eta_p[e][id][nEta]->Add(h_TPH_Res_InvBeta_Eta_p[e][id][et]);
        h_TPH_Res_InvSmearBeta_Eta_p[e][0][et]->Add(h_TPH_Res_InvSmearBeta_Eta_p[e][id][et]);
        h_TPH_Res_InvSmearBeta_Eta_p[e][0][nEta]->Add(h_TPH_Res_InvSmearBeta_Eta_p[e][id][et]);
        h_TPH_Res_InvSmearBeta_Eta_p[e][id][nEta]->Add(h_TPH_Res_InvSmearBeta_Eta_p[e][id][et]);
        h_TPH_InvSmearAndT0Beta_Eta_p[e][0][et]->Add(h_TPH_InvSmearAndT0Beta_Eta_p[e][id][et]);
        h_TPH_InvSmearAndT0Beta_Eta_p[e][0][nEta]->Add(h_TPH_InvSmearAndT0Beta_Eta_p[e][id][et]);
        h_TPH_InvSmearAndT0Beta_Eta_p[e][id][nEta]->Add(h_TPH_InvSmearAndT0Beta_Eta_p[e][id][et]);
        h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][0][et]->Add(h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][id][et]);
        h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][0][nEta]->Add(h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][id][et]);
        h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][id][nEta]->Add(h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][id][et]);
      }
    }
  }

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_TOF.root",outputDir.Data()),"RECREATE");

  h_TPH_T0VsIter->Write();
  h_TPH_UPartIter->Write();
  for (int sc = 0; sc < nSces; sc++){
    if (h_TPH_InitT0VsNtrk[sc]) h_TPH_InitT0VsNtrk[sc]->Write();
    if (h_TPH_InitT0VsNttlhit[sc]) h_TPH_InitT0VsNttlhit[sc]->Write();
    if (h_TPH_NttlhitVsNtrk_wInitT0[sc]) h_TPH_NttlhitVsNtrk_wInitT0[sc]->Write();
    if (h_TPH_NttlhitVsNtrk_wFinalT0[sc]) h_TPH_NttlhitVsNtrk_wFinalT0[sc]->Write();
    if (h_TPH_T0VsNtrk[sc]) h_TPH_T0VsNtrk[sc]->Write();
    if (h_TPH_T0VsNttlhit[sc]) h_TPH_T0VsNttlhit[sc]->Write();
    if (h_TPH_InitDbetaVsP_pi[sc]) h_TPH_InitDbetaVsP_pi[sc]->Write();
    if (h_TPH_InitDbetaVsP_k[sc]) h_TPH_InitDbetaVsP_k[sc]->Write();
    if (h_TPH_InitDbetaVsP_p[sc]) h_TPH_InitDbetaVsP_p[sc]->Write();
    if (h_TPH_DbetaVsP_pi[sc]) h_TPH_DbetaVsP_pi[sc]->Write();
    if (h_TPH_DbetaVsP_k[sc]) h_TPH_DbetaVsP_k[sc]->Write();
    if (h_TPH_DbetaVsP_p[sc]) h_TPH_DbetaVsP_p[sc]->Write();
  }
    
  for (Int_t dir = 0; dir < 3; dir++){
    if (h_TPH_BetaVsP_T0[dir])h_TPH_BetaVsP_T0[dir]->Write();
    for (Int_t id = 0; id < nPID; id++){
      if (h_TPH_DbetaVsP_wT0[id][dir]){
        if (h_TPH_DbetaVsP_wT0[id][dir]->GetEntries() > 0)
          h_TPH_DbetaVsP_wT0[id][dir]->Write();
      }
      if (h_TPH_DbetaVsP_woT0[id][dir]){
        if (h_TPH_DbetaVsP_woT0[id][dir]->GetEntries() > 0)
          h_TPH_DbetaVsP_woT0[id][dir]->Write();
      }
      if (h_TPH_DgenbetaVsGenP[id][dir]){
        if (h_TPH_DgenbetaVsGenP[id][dir]->GetEntries() > 0)
          h_TPH_DgenbetaVsGenP[id][dir]->Write();
      }
      if (h_TPH_DbetaVsP_woT0_test[id][dir]){
        if (h_TPH_DbetaVsP_woT0_test[id][dir]->GetEntries() > 0)
          h_TPH_DbetaVsP_woT0_test[id][dir]->Write();
      }
    }
  }

  for (Int_t et = 0; et<nEta+1; et++){
    for (Int_t e = 0; e< 4; e++){
      for (Int_t id = 0; id < nPID; id++){
        if (h_TPH_InvGenBeta_Eta_p[e][id][et]){
          if (h_TPH_InvGenBeta_Eta_p[e][id][et]->GetEntries() > 0)
          h_TPH_InvGenBeta_Eta_p[e][id][et]->Write();
        }
        if (h_TPH_InvBeta_Eta_p[e][id][et]){
          if (h_TPH_InvBeta_Eta_p[e][id][et]->GetEntries() > 0)
            h_TPH_InvBeta_Eta_p[e][id][et]->Write();
        }
        if (h_TPH_InvSmearBeta_Eta_p[e][id][et]){
          if (h_TPH_InvSmearBeta_Eta_p[e][id][et]->GetEntries() > 0)
            h_TPH_InvSmearBeta_Eta_p[e][id][et]->Write();
        }
        if (h_TPH_InvSmearAndT0Beta_Eta_p[e][id][et]){
          if (h_TPH_InvSmearAndT0Beta_Eta_p[e][id][et]->GetEntries() > 0)
            h_TPH_InvSmearAndT0Beta_Eta_p[e][id][et]->Write();
        }
        if (h_TPH_Res_InvBeta_Eta_p[e][id][et]){
          if (h_TPH_Res_InvBeta_Eta_p[e][id][et]->GetEntries() > 0)
            h_TPH_Res_InvBeta_Eta_p[e][id][et]->Write();
        }
        if (h_TPH_Res_InvSmearBeta_Eta_p[e][id][et]){
          if (h_TPH_Res_InvSmearBeta_Eta_p[e][id][et]->GetEntries() > 0)
            h_TPH_Res_InvSmearBeta_Eta_p[e][id][et]->Write();
        }
        if (h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][id][et]){
          if (h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][id][et]->GetEntries() > 0)
            h_TPH_Res_InvSmearAndT0Beta_Eta_p[e][id][et]->Write();
        }
      }
    }
  }
}

