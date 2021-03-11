#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
#include <fstream>

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);



void analyseTreeForTrackingResolAndPID(
                                        TString studyname           = "defaultLBL",
                                        TString singleFileName      = "/home/fbock/eic/Singularity/Fun4All_G4_LBLDetector_Org/10000/G4LBLDetector_g4tracking_eval.root",
                                        TString inputFileNameList   = "",
                                        Bool_t isLANL               = kFALSE,
                                        Bool_t hasTiming            = kTRUE
                                      ){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput           = ReturnDateStringForOutput();

  TString outputDir               = "rootfiles/";
  gSystem->Exec("mkdir -p "+outputDir);

  // **********************************************************************************************
  // ***************************** general variable setups ****************************************
  // **********************************************************************************************
  // random number generator
  TRandom3 r3;
  const double c                  = 29.9792458; // cm/ns
  const double sigmat             = 20e-3; // ns

  
  
  // regions: 0 = Backwards (electron going), 1 = central barrel, 2 = forward (hadron going)
  Int_t maxLayer[3]               = {5, 6, 5};
  Int_t maxLayerTTL[3]            = {2, 2, 3};
  Int_t layerOffsetLBL[3]         = {30, 10, 20 };
  Int_t layerTTLBEC[3]            = {2, 1, 2 };
  
  if (studyname.Contains("FTTLS2LC")){
    maxLayerTTL[2] = 2;
  } else if( studyname.Contains("FTTLSE2LC")){
    maxLayerTTL[2] = 2;
    layerTTLBEC[2] = 1;
  } else if (studyname.Contains("FTTLSE1") ){
    maxLayerTTL[2] = 1;
    layerTTLBEC[2] = 1;
  }
  
  if (studyname.Contains("ETTLSE1")){
    maxLayerTTL[0] = 1;
    layerTTLBEC[0] = 1;
  }
  if (studyname.Contains("CTTLSE1")){
    maxLayerTTL[1] = 1;
    layerTTLBEC[1] = 1;
  }
    
  if (isLANL){
    maxLayer[0]   = 0;
  }
    
  TString nameLayerLANL[3]        = {"", "BARREL", "FST" };  
  TString nameLayerLBL[3]         = {"BACKWARD", "CENTRAL", "FORWARD" };
  TString nameLayerTTL[3]         = {"ETTL", "CTTL", "FTTL" };
  
  // **********************************************************************************************
  // ********************************** Set up histos *********************************************
  // **********************************************************************************************
  TH1D* hNEvents                                = new TH1D("nEvents","",1,0,1);
  TH1D* hNTracks                                = new TH1D("nTracks","",200,-0.5,199.5);
  TH2F* h_tracks_reso_pT[5][nEta+1][6]          = {{{NULL}}};
  TH2F* h_tracks_reso_p[5][nEta+1][6]           = {{NULL}};
  TH2F* h_tracks_resoEta_pT[5][nEta+1]          = {{NULL}};
  TH2F* h_tracks_resoPhi_pT[5][nEta+1]          = {{NULL}};
  TH2F* h_tracksTrue_Eta_pT[6][12]              = {{NULL}};
  TH2F* h_tracksRec_Eta_pT[6][12]               = {{NULL}};
  TH2F* h_tracksTrue_Eta_p[6][12]               = {{NULL}};
  TH2F* h_tracksRec_Eta_p[6][12]                = {{NULL}};
  TH2F* h_InvGenBeta_Eta_p[4][6][nEta+1]        = {{NULL}};
  TH2F* h_InvBeta_Eta_p[4][6][nEta+1]           = {{NULL}};
  TH2F* h_InvSmearBeta_Eta_p[4][6][nEta+1]      = {{NULL}};
  TH2F* h_Res_InvBeta_Eta_p[4][6][nEta+1]       = {{NULL}};
  TH2F* h_Res_InvSmearBeta_Eta_p[4][6][nEta+1]  = {{NULL}};
  
  for (Int_t id = 0; id < 6; id++){
    for (Int_t k = 0; k < 12; k++){
      std::cout << id << "\t" << partName[id].Data() << "\t" << k << "\t" << nameCuts[k].Data()  << "\t" 
                << Form("h_tracks_%s_%s_True_Eta_pT", partName[id].Data(), nameCuts[k].Data()) ;
      h_tracksTrue_Eta_pT[id][k]         = new TH2F(Form("h_tracks_%s_%s_True_Eta_pT", partName[id].Data(), nameCuts[k].Data() ), "", 200, 0, 20, 200, -4, 4.);
//       h_tracksTrue_Eta_pT[id][k]->Sumw2();
      h_tracksRec_Eta_pT[id][k]          = new TH2F(Form("h_tracks_%s_%s_Rec_Eta_pT", partName[id].Data(), nameCuts[k].Data() ), "", 200, 0, 20, 200, -4, 4.);
//       h_tracksRec_Eta_pT[id][k]->Sumw2();
      h_tracksTrue_Eta_p[id][k]          = new TH2F(Form("h_tracks_%s_%s_True_Eta_p", partName[id].Data(), nameCuts[k].Data() ), "", nBinsP, binningP,  200, -4, 4.);
//       h_tracksTrue_Eta_p[id][k]->Sumw2();
      h_tracksRec_Eta_p[id][k]           = new TH2F(Form("h_tracks_%s_%s_Rec_Eta_p", partName[id].Data(), nameCuts[k].Data() ), "", nBinsP, binningP, 200, -4, 4.);
//       h_tracksRec_Eta_p[id][k]->Sumw2();
      
      std::cout << h_tracksTrue_Eta_pT[id][k]->GetName() << "\t" << h_tracksRec_Eta_pT[id][k]->GetName()  << "\t" << h_tracksTrue_Eta_p[id][k]->GetName() << "\t" << h_tracksRec_Eta_p[id][k]->GetName()<< std::endl;
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
      for (Int_t id = 0; id < 6; id++){
        h_tracks_reso_pT[i][et][id]     = new TH2F(Form("h_tracks_reso_pT_%s_%s_%d",nameResoAdd[i].Data(), partName[id].Data(),et), 
                                              Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#it{p}_{T,rec}-#it{p}_{T,MC})/#it{p}_{T,MC}",etaMin, etaMax), 
                                              200, 0, 20, 250, -0.5, 0.5);
        h_tracks_reso_p[i][et][id]         = new TH2F(Form("h_tracks_reso_p_%s_%s_%d",nameResoAdd[i].Data(), partName[id].Data(),et), 
                                                Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (#it{p}_{rec}-#it{p}_{MC})/#it{p}_{MC}", etaMin, etaMax), 
                                                nBinsP, binningP, 250, -0.5, 0.5);
      }
      h_tracks_resoEta_pT[i][et]     = new TH2F(Form("h_tracks_reso_Eta_pT_%s_%d",nameResoAdd[i].Data(),et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#eta_{rec}-#eta_{MC})/#eta_{MC}",
                                                                                      etaMin, etaMax), 200, 0, 20, 500, -0.1, 0.1);
      h_tracks_resoPhi_pT[i][et]     = new TH2F(Form("h_tracks_reso_Phi_pT_%s_%d",nameResoAdd[i].Data(),et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#phi_{rec}-#phi_{MC})/#phi_{MC}",
                                                                                      etaMin, etaMax), 200, 0, 20, 500, -0.5, 0.5);
    }
    for (Int_t e = 0; e< 4; e++){
      for (Int_t id = 0; id < 6; id++){
        
        h_InvGenBeta_Eta_p[e][id][et]      = new TH2F(Form("h_InvGenBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(), partName[id].Data(), et), 
                                                Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{MC}", etaMin, etaMax), 
                                                 nBinsPLow, binningP, 300, 0.9, 1.5);

        h_InvBeta_Eta_p[e][id][et]         = new TH2F(Form("h_InvBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(), partName[id].Data(), et), 
                                                Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{rec}", etaMin, etaMax), 
                                                nBinsPLow, binningP, 300, 0.9, 1.5);
        h_InvSmearBeta_Eta_p[e][id][et]    = new TH2F(Form("h_InvSmearBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{smear}", etaMin, etaMax), 
                                                nBinsPLow, binningP, 300, 0.9, 1.5);
        h_Res_InvBeta_Eta_p[e][id][et]     = new TH2F(Form("h_Res_InvBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (1/#beta_{rec} - 1/#beta_{MC})/(1/#beta_{MC}) ", etaMin, etaMax), 
                                                nBinsPLow, binningP, 500, -0.1, 0.1);
        h_Res_InvSmearBeta_Eta_p[e][id][et]= new TH2F(Form("h_Res_InvSmearBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (1/#beta_{smear} - 1/#beta_{MC})/(1/#beta_{MC}) ", etaMin, etaMax), 
                                                nBinsPLow, binningP, 500, -0.1, 0.1);
      }
    }
  }
  
  
  // **********************************************************************************************
  //************************** Read data **********************************************************
  // **********************************************************************************************
  TChain* t_tracks = new TChain("tracks", "tracks");
  if (singleFileName.CompareTo("") != 0){
    t_tracks->Add(singleFileName.Data());
  } else {
    ifstream in(inputFileNameList.Data());
    std::cout<<"File: "<<inputFileNameList.Data()<<std::endl;
    std::cout<<"reading input files "<<std::endl;
    // general number of triggers set
    while(!in.eof()){
      TString tempFilename = "";
      in >> tempFilename;
      if (tempFilename.CompareTo("") != 0 ){
        std::cout << "adding " << tempFilename.Data() << std::endl;
        t_tracks->Add(tempFilename.Data());
      }
    }
  }

  if(t_tracks){
    Int_t eventID;
    Float_t trk_px,trk_py,trk_pz,trk_px_true,trk_py_true,trk_pz_true;
    Int_t pid;
    // recHit and trueHit arrays as follows: [etaRegion][layer][dimension] 
    // currently maximum number of layers = 6 (central barrel) organized inner to outer
    Float_t recHit[3][6][4];
    Float_t trueHit[3][6][4];
    // recHit and trueHit arrays as follows: [etaRegion][layer][dimension]
    // currently maximum number of layers = 3 (forward) organized inner to outer
    Float_t recHitTTL[3][3][4];
    Float_t trueHitTTL[3][3][4];
    // smeared time for layers [etaRegion][layer]
    Float_t timeSmeared[3][3];
    // 0 true, 1 rec, 2 smeared
    Float_t beta[3];
    
    t_tracks->SetBranchAddress("event", &eventID);
    t_tracks->SetBranchAddress("gflavor", &pid);
    t_tracks->SetBranchAddress("px", &trk_px);
    t_tracks->SetBranchAddress("py", &trk_py);
    t_tracks->SetBranchAddress("pz", &trk_pz);
    t_tracks->SetBranchAddress("gpx", &trk_px_true);
    t_tracks->SetBranchAddress("gpy", &trk_py_true);
    t_tracks->SetBranchAddress("gpz", &trk_pz_true);
    if (isLANL){
      for (Int_t eR = 0; eR < 3; eR++){
        for (Int_t l = 0; l < maxLayer[eR]; l++){
          TString nameCurrBase = Form("%s_%d", nameLayerLANL[eR].Data(), l);
          std::cout << eR << "\t" << l << "\t" << nameCurrBase.Data() << std::endl;
          t_tracks->SetBranchAddress(Form("%s_proj_x", nameCurrBase.Data()), &trueHit[eR][l][0]);
          t_tracks->SetBranchAddress(Form("%s_proj_y", nameCurrBase.Data()), &trueHit[eR][l][1]);
          t_tracks->SetBranchAddress(Form("%s_proj_z", nameCurrBase.Data()), &trueHit[eR][l][2]);
          t_tracks->SetBranchAddress(Form("%s_proj_t", nameCurrBase.Data()), &trueHit[eR][l][3]);
          t_tracks->SetBranchAddress(Form("%s_x", nameCurrBase.Data()), &recHit[eR][l][0]);
          t_tracks->SetBranchAddress(Form("%s_y", nameCurrBase.Data()), &recHit[eR][l][1]);
          t_tracks->SetBranchAddress(Form("%s_z", nameCurrBase.Data()), &recHit[eR][l][2]);
          t_tracks->SetBranchAddress(Form("%s_t", nameCurrBase.Data()), &recHit[eR][l][3]);
        }
      }
    } else {
      for (Int_t eR = 0; eR < 3; eR++){
        for (Int_t l = 0; l < maxLayer[eR]; l++){
          
          TString nameCurrBase = Form("LBLVTX_%s_%d", nameLayerLBL[eR].Data(), layerOffsetLBL[eR]+l);
          std::cout << eR << "\t" << l << "\t" << nameCurrBase.Data() << std::endl;
          t_tracks->SetBranchAddress(Form("%s_proj_x", nameCurrBase.Data()), &trueHit[eR][l][0]);
          t_tracks->SetBranchAddress(Form("%s_proj_y", nameCurrBase.Data()), &trueHit[eR][l][1]);
          t_tracks->SetBranchAddress(Form("%s_proj_z", nameCurrBase.Data()), &trueHit[eR][l][2]);
          t_tracks->SetBranchAddress(Form("%s_proj_t", nameCurrBase.Data()), &trueHit[eR][l][3]);
          t_tracks->SetBranchAddress(Form("%s_x", nameCurrBase.Data()), &recHit[eR][l][0]);
          t_tracks->SetBranchAddress(Form("%s_y", nameCurrBase.Data()), &recHit[eR][l][1]);
          t_tracks->SetBranchAddress(Form("%s_z", nameCurrBase.Data()), &recHit[eR][l][2]);
          t_tracks->SetBranchAddress(Form("%s_t", nameCurrBase.Data()), &recHit[eR][l][3]);      
        }
      }
    }

    if (hasTiming){
      for (Int_t eR = 0; eR < 3; eR++){
        for (Int_t l = 0;  l < maxLayerTTL[eR]; l++){
          TString nameCurrBase = Form("%s_%d", nameLayerTTL[eR].Data(), l);
          std::cout << eR << "\t" << l << "\t" << nameCurrBase.Data() << std::endl;
          t_tracks->SetBranchAddress(Form("%s_proj_x", nameCurrBase.Data()), &trueHitTTL[eR][l][0]);
          t_tracks->SetBranchAddress(Form("%s_proj_y", nameCurrBase.Data()), &trueHitTTL[eR][l][1]);
          t_tracks->SetBranchAddress(Form("%s_proj_z", nameCurrBase.Data()), &trueHitTTL[eR][l][2]);
          t_tracks->SetBranchAddress(Form("%s_proj_t", nameCurrBase.Data()), &trueHitTTL[eR][l][3]);
          t_tracks->SetBranchAddress(Form("%s_x", nameCurrBase.Data()), &recHitTTL[eR][l][0]);
          t_tracks->SetBranchAddress(Form("%s_y", nameCurrBase.Data()), &recHitTTL[eR][l][1]);
          t_tracks->SetBranchAddress(Form("%s_z", nameCurrBase.Data()), &recHitTTL[eR][l][2]);
          t_tracks->SetBranchAddress(Form("%s_t", nameCurrBase.Data()), &recHitTTL[eR][l][3]);
        }
      }
    }
    Int_t nentries = Int_t(t_tracks->GetEntries());
    Float_t lastEvt = 0;
    Int_t currEventID = -1;
    Int_t nTracksCurr = 0;
    for (Long_t i = 0; i < nentries; ++i) {
      if (t_tracks->LoadTree(i) < 0)
        break;
//       std::cout << i << std::endl;
      t_tracks->GetEntry(i);
      
      if (eventID != currEventID){
        hNEvents->AddBinContent(1,1);
        if (i > 0) hNTracks->Fill(nTracksCurr);
        nTracksCurr = 1;
        if (i == nentries-1) hNTracks->Fill(nTracksCurr);
        currEventID = eventID;
      } else if (i == nentries-1 ){
        hNTracks->Fill(nTracksCurr);
      } else {
        nTracksCurr++; 
      }
      TVector3 vec_mom_true(trk_px, trk_py, trk_pz);
      TVector3 vec_mom(trk_px_true, trk_py_true, trk_pz_true);

      Int_t parIdx = -1;
      if(TMath::Abs(pid) == 11)   parIdx = 1;
      if(TMath::Abs(pid) == 13)   parIdx = 2;
      if(TMath::Abs(pid) == 211)  parIdx = 3;
      if(TMath::Abs(pid) == 321)  parIdx = 4;
      if(TMath::Abs(pid) == 2212) parIdx = 5;

      
      TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(pid);
      Double_t mass = part->Mass();

      Double_t truePt     = TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2));
      Double_t trueP      = TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)+TMath::Power(trk_pz_true,2));
      Double_t trueEta    = vec_mom_true.Eta();
      Double_t truePhi    = vec_mom_true.Phi();
      Double_t recPt      = TMath::Sqrt(TMath::Power(trk_px,2)+TMath::Power(trk_py,2));
      Double_t recP       = TMath::Sqrt(TMath::Power(trk_px,2)+TMath::Power(trk_py,2)+TMath::Power(trk_pz,2));
      Double_t recEta     = vec_mom.Eta();
      Double_t recPhi     = vec_mom.Phi();
     
      Int_t trk_hits[3]         = {0, 0, 0};
      Int_t trk_hitsTime[3]     = {0, 0, 0};
      Int_t trk_hitsTracker[3]  = {0, 0, 0};
      Int_t totTimeHits     = 0;
      Int_t totTimeTracker  = 0;
      Bool_t kHitFL         = kFALSE;
      Bool_t kHitBEMC       = kFALSE;
      Bool_t kHitAEMC       = kFALSE;

      beta[0] = 0.;
      beta[1] = 0.;
      beta[2] = 0.;
      
      for (Int_t eR = 0; eR < 3; eR++){
        for (Int_t l = 0; l < maxLayer[eR]; l++){
          if (recHit[eR][l][3] > -9999){
            trk_hits[eR]++;  
            trk_hitsTracker[eR]++;  
          }
          if (l < 2 && recHit[eR][l][3] > -9999 ) kHitFL = kTRUE;
        }

        if (hasTiming){
          for (Int_t lt = 0; lt < maxLayerTTL[eR]; lt++){
            if (recHitTTL[eR][lt][3] > -9999){ 
              trk_hits[eR]++;  
              trk_hitsTime[eR]++;
              timeSmeared[eR][lt] = r3.Gaus(recHitTTL[eR][lt][3], sigmat);
              beta[1] += c*recHitTTL[eR][lt][3]/trueHitTTL[eR][lt][3]; 
              beta[2] += c*timeSmeared[eR][lt]/trueHitTTL[eR][lt][3]; 
            }
            if (lt < layerTTLBEC[eR] && recHitTTL[eR][lt][3] > -9999 ) kHitBEMC = kTRUE;
            if (lt >= layerTTLBEC[eR] && recHitTTL[eR][lt][3] > -9999 ) kHitAEMC = kTRUE; 
          }
        }
      }
      for (Int_t eR = 0; eR < 3; eR++){
        totTimeHits += trk_hitsTime[eR];
        totTimeTracker += trk_hitsTracker[eR];
      }
      
      if (hasTiming){
        beta[0] = sqrt(trueP*trueP + mass*mass)/trueP;
        if (totTimeHits > 0){
          beta[1] = beta[1]/totTimeHits;
          beta[2] = beta[2]/totTimeHits;
  //         std::cout << beta[0] << "\t" << beta[1] << "\t" << beta[2] << std::endl;
        } else {
          beta[1] = -10000;
          beta[2] = -10000;
        }
      }
      
      // determine eta bin
      Int_t et = 0;
      while (partEta[et+1] < trueEta && et < nEta) et++;
      
      
      h_tracks_reso_pT[0][et][parIdx]->Fill(truePt,(recPt-truePt)/truePt);
      h_tracks_reso_p[0][et][parIdx]->Fill(trueP,(recP-trueP)/trueP);
      h_tracks_resoEta_pT[0][et]->Fill(truePt,(recEta-trueEta));
      h_tracks_resoPhi_pT[0][et]->Fill(truePt,(recPhi-truePhi));
      if (totTimeHits > 0){
        h_tracks_reso_pT[2][et][parIdx]->Fill(truePt,(recPt-truePt)/truePt);
        h_tracks_reso_p[2][et][parIdx]->Fill(trueP,(recP-trueP)/trueP);
        h_tracks_resoEta_pT[2][et]->Fill(truePt,(recEta-trueEta));
        h_tracks_resoPhi_pT[2][et]->Fill(truePt,(recPhi-truePhi));
      } else {
        h_tracks_reso_pT[1][et][parIdx]->Fill(truePt,(recPt-truePt)/truePt);
        h_tracks_reso_p[1][et][parIdx]->Fill(trueP,(recP-trueP)/trueP);        
        h_tracks_resoEta_pT[1][et]->Fill(truePt,(recEta-trueEta));
        h_tracks_resoPhi_pT[1][et]->Fill(truePt,(recPhi-truePhi));
      }
      if (totTimeTracker >= 2){
        h_tracks_reso_pT[3][et][parIdx]->Fill(truePt,(recPt-truePt)/truePt);
        h_tracks_reso_p[3][et][parIdx]->Fill(trueP,(recP-trueP)/trueP);        
        h_tracks_resoEta_pT[3][et]->Fill(truePt,(recEta-trueEta));
        h_tracks_resoPhi_pT[3][et]->Fill(truePt,(recPhi-truePhi));
      }
      if (totTimeTracker >= 3){
        h_tracks_reso_pT[4][et][parIdx]->Fill(truePt,(recPt-truePt)/truePt);
        h_tracks_reso_p[4][et][parIdx]->Fill(trueP,(recP-trueP)/trueP);        
        h_tracks_resoEta_pT[4][et]->Fill(truePt,(recEta-trueEta));
        h_tracks_resoPhi_pT[4][et]->Fill(truePt,(recPhi-truePhi));
      }
            
      if (hasTiming && totTimeHits > 0){
        // fill beta hists
        h_InvGenBeta_Eta_p[0][parIdx][et]->Fill(trueP, beta[0]);
        h_InvBeta_Eta_p[0][parIdx][et]->Fill(trueP, beta[1]);
        h_InvSmearBeta_Eta_p[0][parIdx][et]->Fill(trueP, beta[2]);
        h_Res_InvBeta_Eta_p[0][parIdx][et]->Fill(trueP, (beta[1]-beta[0])/beta[0]);
        h_Res_InvSmearBeta_Eta_p[0][parIdx][et]->Fill(trueP, (beta[2]-beta[0])/beta[0]);
        // fill beta hists after EMC
        if ( kHitAEMC){
          h_InvGenBeta_Eta_p[1][parIdx][et]->Fill(trueP, beta[0]);
          h_InvBeta_Eta_p[1][parIdx][et]->Fill(trueP, beta[1]);
          h_InvSmearBeta_Eta_p[1][parIdx][et]->Fill(trueP, beta[2]);
          h_Res_InvBeta_Eta_p[1][parIdx][et]->Fill(trueP, (beta[1]-beta[0])/beta[0]);
          h_Res_InvSmearBeta_Eta_p[1][parIdx][et]->Fill(trueP, (beta[2]-beta[0])/beta[0]);
        }
        if (totTimeTracker >= 3){
          h_InvGenBeta_Eta_p[2][parIdx][et]->Fill(trueP, beta[0]);
          h_InvBeta_Eta_p[2][parIdx][et]->Fill(trueP, beta[1]);
          h_InvSmearBeta_Eta_p[2][parIdx][et]->Fill(trueP, beta[2]);
          h_Res_InvBeta_Eta_p[2][parIdx][et]->Fill(trueP, (beta[1]-beta[0])/beta[0]);
          h_Res_InvSmearBeta_Eta_p[2][parIdx][et]->Fill(trueP, (beta[2]-beta[0])/beta[0]);
        }
        if (totTimeTracker >= 3 && kHitAEMC){
          h_InvGenBeta_Eta_p[3][parIdx][et]->Fill(trueP, beta[0]);
          h_InvBeta_Eta_p[3][parIdx][et]->Fill(trueP, beta[1]);
          h_InvSmearBeta_Eta_p[3][parIdx][et]->Fill(trueP, beta[2]);
          h_Res_InvBeta_Eta_p[3][parIdx][et]->Fill(trueP, (beta[1]-beta[0])/beta[0]);
          h_Res_InvSmearBeta_Eta_p[3][parIdx][et]->Fill(trueP, (beta[2]-beta[0])/beta[0]);
        }
      }

      // all tracks "N"
      h_tracksTrue_Eta_pT[parIdx][0]->Fill(truePt, trueEta);
      h_tracksRec_Eta_pT[parIdx][0]->Fill(recPt, recEta);
      h_tracksTrue_Eta_p[parIdx][0]->Fill(trueP, trueEta);
      h_tracksRec_Eta_p[parIdx][0]->Fill(recP, recEta);
      // at least 3 layers
      if ((trk_hits[0]+trk_hits[1]+trk_hits[2]) >= 3){
        h_tracksTrue_Eta_pT[parIdx][1]->Fill(truePt, trueEta);
        h_tracksRec_Eta_pT[parIdx][1]->Fill(recPt, recEta);
        h_tracksTrue_Eta_p[parIdx][1]->Fill(trueP, trueEta);
        h_tracksRec_Eta_p[parIdx][1]->Fill(recP, recEta);
        // at least 3 layers & inner most tracking layer fw
        if (kHitFL){ 
          h_tracksTrue_Eta_pT[parIdx][2]->Fill(truePt, trueEta);
          h_tracksRec_Eta_pT[parIdx][2]->Fill(recPt, recEta);
          h_tracksTrue_Eta_p[parIdx][2]->Fill(trueP, trueEta);
          h_tracksRec_Eta_p[parIdx][2]->Fill(recP, recEta);
        }
      }
      // only tracker layers
      if (totTimeHits < 1){
        h_tracksTrue_Eta_pT[parIdx][3]->Fill(truePt, trueEta);
        h_tracksRec_Eta_pT[parIdx][3]->Fill(recPt, recEta);
        h_tracksTrue_Eta_p[parIdx][3]->Fill(trueP, trueEta);
        h_tracksRec_Eta_p[parIdx][3]->Fill(recP, recEta);
        // only tracker layers & inner most tracking layer fw
        if (kHitFL){ 
          h_tracksTrue_Eta_pT[parIdx][4]->Fill(truePt, trueEta);
          h_tracksRec_Eta_pT[parIdx][4]->Fill(recPt, recEta);
          h_tracksTrue_Eta_p[parIdx][4]->Fill(trueP, trueEta);
          h_tracksRec_Eta_p[parIdx][4]->Fill(recP, recEta);
        }        
      }
      // timing layer hit before ECal but not after
      if (kHitBEMC && !kHitAEMC){
        h_tracksTrue_Eta_pT[parIdx][5]->Fill(truePt, trueEta);
        h_tracksRec_Eta_pT[parIdx][5]->Fill(recPt, recEta);
        h_tracksTrue_Eta_p[parIdx][5]->Fill(trueP, trueEta);
        h_tracksRec_Eta_p[parIdx][5]->Fill(recP, recEta);
        // timing layer hit before ECal but not after & inner most tracking layer fw
        if (kHitFL){ 
          h_tracksTrue_Eta_pT[parIdx][6]->Fill(truePt, trueEta);
          h_tracksRec_Eta_pT[parIdx][6]->Fill(recPt, recEta);
          h_tracksTrue_Eta_p[parIdx][6]->Fill(trueP, trueEta);
          h_tracksRec_Eta_p[parIdx][6]->Fill(recP, recEta);
        }        
      }
      // timing layer hit after ECal 
      if (kHitAEMC){
        h_tracksTrue_Eta_pT[parIdx][7]->Fill(truePt, trueEta);
        h_tracksRec_Eta_pT[parIdx][7]->Fill(recPt, recEta);
        h_tracksTrue_Eta_p[parIdx][7]->Fill(trueP, trueEta);
        h_tracksRec_Eta_p[parIdx][7]->Fill(recP, recEta);
        // timing layer hit after ECal & inner most tracking layer fw
        if (kHitFL){ 
          h_tracksTrue_Eta_pT[parIdx][8]->Fill(truePt, trueEta);
          h_tracksRec_Eta_pT[parIdx][8]->Fill(recPt, recEta);
          h_tracksTrue_Eta_p[parIdx][8]->Fill(trueP, trueEta);
          h_tracksRec_Eta_p[parIdx][8]->Fill(recP, recEta);
        }        
      }
   
      if (totTimeHits > 0){
        h_tracksTrue_Eta_pT[parIdx][9]->Fill(truePt, trueEta);
        h_tracksRec_Eta_pT[parIdx][9]->Fill(recPt, recEta);
        h_tracksTrue_Eta_p[parIdx][9]->Fill(trueP, trueEta);
        h_tracksRec_Eta_p[parIdx][9]->Fill(recP, recEta);
      }
      if (totTimeTracker >= 2){
        h_tracksTrue_Eta_pT[parIdx][10]->Fill(truePt, trueEta);
        h_tracksRec_Eta_pT[parIdx][10]->Fill(recPt, recEta);
        h_tracksTrue_Eta_p[parIdx][10]->Fill(trueP, trueEta);
        h_tracksRec_Eta_p[parIdx][10]->Fill(recP, recEta);
      }
      if (totTimeTracker >= 3){
        h_tracksTrue_Eta_pT[parIdx][11]->Fill(truePt, trueEta);
        h_tracksRec_Eta_pT[parIdx][11]->Fill(recPt, recEta);
        h_tracksTrue_Eta_p[parIdx][11]->Fill(trueP, trueEta);
        h_tracksRec_Eta_p[parIdx][11]->Fill(recP, recEta);
      }
    }
  }

  // ***********************************************************************************************
  // ************************** sum up hists correctly *********************************************
  // ***********************************************************************************************
  for (Int_t id = 1; id < 6; id++){
    for (Int_t k = 0; k < 12; k++){
      if (id == 1){
        h_tracksTrue_Eta_pT[0][k]->Sumw2();
        h_tracksRec_Eta_pT[0][k]->Sumw2();
        h_tracksTrue_Eta_p[0][k]->Sumw2();
        h_tracksRec_Eta_p[0][k]->Sumw2();
      }
      h_tracksTrue_Eta_pT[0][k]->Add(h_tracksTrue_Eta_pT[id][k]);
      h_tracksRec_Eta_pT[0][k]->Add(h_tracksRec_Eta_pT[id][k]);
      h_tracksTrue_Eta_p[0][k]->Add(h_tracksTrue_Eta_p[id][k]);
      h_tracksRec_Eta_p[0][k]->Add(h_tracksRec_Eta_p[id][k]);     
    }
    for (Int_t i = 0; i < 5; i++){
      for (Int_t et = 0; et < nEta; et++){
        if (id == 1 && et == 0){
          h_tracks_reso_pT[i][et][0]->Sumw2();
          h_tracks_reso_pT[i][nEta][0]->Sumw2();
          h_tracks_reso_pT[i][nEta][id]->Sumw2();
          h_tracks_reso_p[i][et][0]->Sumw2();
          h_tracks_reso_p[i][nEta][0]->Sumw2();
          h_tracks_reso_p[i][nEta][id]->Sumw2();
        }
        h_tracks_reso_pT[i][et][0]->Add(h_tracks_reso_pT[i][et][id]);
        h_tracks_reso_pT[i][nEta][0]->Add(h_tracks_reso_pT[i][et][id]);
        h_tracks_reso_pT[i][nEta][id]->Add(h_tracks_reso_pT[i][et][id]);
        h_tracks_reso_p[i][et][0]->Add(h_tracks_reso_p[i][et][id]);
        h_tracks_reso_p[i][nEta][0]->Add(h_tracks_reso_p[i][et][id]);
        h_tracks_reso_p[i][nEta][id]->Add(h_tracks_reso_p[i][et][id]);
      }
    }
    for (Int_t e = 0; e < 4; e++){
      for (Int_t et = 0; et < nEta; et++){
        if (et == 0 && id == 0){
          h_InvGenBeta_Eta_p[e][0][et]->Sumw2();
          h_InvGenBeta_Eta_p[e][0][nEta]->Sumw2();
          h_InvGenBeta_Eta_p[e][id][nEta]->Sumw2();
          h_InvBeta_Eta_p[e][0][et]->Sumw2();
          h_InvBeta_Eta_p[e][0][nEta]->Sumw2();
          h_InvBeta_Eta_p[e][id][nEta]->Sumw2();
          h_InvSmearBeta_Eta_p[e][0][et]->Sumw2();
          h_InvSmearBeta_Eta_p[e][0][nEta]->Sumw2();
          h_InvSmearBeta_Eta_p[e][id][nEta]->Sumw2();
          h_Res_InvBeta_Eta_p[e][0][et]->Sumw2();
          h_Res_InvBeta_Eta_p[e][0][nEta]->Sumw2();
          h_Res_InvBeta_Eta_p[e][id][nEta]->Sumw2();
          h_Res_InvSmearBeta_Eta_p[e][0][et]->Sumw2();
          h_Res_InvSmearBeta_Eta_p[e][0][nEta]->Sumw2();
          h_Res_InvSmearBeta_Eta_p[e][id][nEta]->Sumw2();
          
        }
        h_InvGenBeta_Eta_p[e][0][et]->Add(h_InvGenBeta_Eta_p[e][id][et]);
        h_InvGenBeta_Eta_p[e][0][nEta]->Add(h_InvGenBeta_Eta_p[e][id][et]);
        h_InvGenBeta_Eta_p[e][id][nEta]->Add(h_InvGenBeta_Eta_p[e][id][et]);
        h_InvBeta_Eta_p[e][0][et]->Add(h_InvBeta_Eta_p[e][id][et]);
        h_InvBeta_Eta_p[e][0][nEta]->Add(h_InvBeta_Eta_p[e][id][et]);
        h_InvBeta_Eta_p[e][id][nEta]->Add(h_InvBeta_Eta_p[e][id][et]);
        h_InvSmearBeta_Eta_p[e][0][et]->Add(h_InvSmearBeta_Eta_p[e][id][et]);
        h_InvSmearBeta_Eta_p[e][0][nEta]->Add(h_InvSmearBeta_Eta_p[e][id][et]);
        h_InvSmearBeta_Eta_p[e][id][nEta]->Add(h_InvSmearBeta_Eta_p[e][id][et]);
        h_Res_InvBeta_Eta_p[e][0][et]->Add(h_Res_InvBeta_Eta_p[e][id][et]);
        h_Res_InvBeta_Eta_p[e][0][nEta]->Add(h_Res_InvBeta_Eta_p[e][id][et]);
        h_Res_InvBeta_Eta_p[e][id][nEta]->Add(h_Res_InvBeta_Eta_p[e][id][et]);
        h_Res_InvSmearBeta_Eta_p[e][0][et]->Add(h_Res_InvSmearBeta_Eta_p[e][id][et]);
        h_Res_InvSmearBeta_Eta_p[e][0][nEta]->Add(h_Res_InvSmearBeta_Eta_p[e][id][et]);
        h_Res_InvSmearBeta_Eta_p[e][id][nEta]->Add(h_Res_InvSmearBeta_Eta_p[e][id][et]);
      }
    }
  }
  for (Int_t i = 0; i < 5; i++){
    for (Int_t et = 0; et < nEta; et++){
      if (et == 0){
        h_tracks_resoEta_pT[i][nEta]->Sumw2();
        h_tracks_resoPhi_pT[i][nEta]->Sumw2();
      }
      h_tracks_resoEta_pT[i][nEta]->Add(h_tracks_resoEta_pT[i][et]);
      h_tracks_resoPhi_pT[i][nEta]->Add(h_tracks_resoPhi_pT[i][et]);
    }
  }
  
  // ***********************************************************************************************
  // **************************** Write outputs to file ********************************************
  // ***********************************************************************************************
  cout << "writing" << endl;
  TString nameOutput      = Form("%s/histosTrackingPerformance_%s.root",outputDir.Data(), studyname.Data());
  TFile* fOutput          = new TFile(nameOutput,"RECREATE");
    hNEvents->Write();
    hNTracks->Write();
    for (Int_t et = 0; et < nEta+1; et++){
      for (Int_t i = 0; i< 5; i++){
        for (Int_t id = 0; id < 6; id++) {
          h_tracks_reso_pT[i][et][id]->Write();
          h_tracks_reso_p[i][et][id]->Write();
        }
        h_tracks_resoEta_pT[i][et]->Write();
        h_tracks_resoPhi_pT[i][et]->Write();
      }
      if (hasTiming){
        for (Int_t id = 0; id < 6; id++) {
          for (Int_t e = 0; e < 4; e++){
            h_InvBeta_Eta_p[e][id][et]->Write();
            h_InvGenBeta_Eta_p[e][id][et]->Write();
            h_InvSmearBeta_Eta_p[e][id][et]->Write();
            h_Res_InvBeta_Eta_p[e][id][et]->Write();
            h_Res_InvSmearBeta_Eta_p[e][id][et]->Write();
          }
        }
      }
    }
    
    for (Int_t id = 0; id < 6; id++){
      for (Int_t k = 0; k < 12; k++){
        std::cout << id << "\t" << k << "\t" << h_tracksTrue_Eta_pT[id][k]->GetName() << "\t" << h_tracksRec_Eta_pT[id][k]->GetName()  << "\t" << h_tracksTrue_Eta_p[id][k]->GetName() << "\t" << h_tracksRec_Eta_p[id][k]->GetName()<< std::endl;
        h_tracksTrue_Eta_pT[id][k]->Write();
        h_tracksRec_Eta_pT[id][k]->Write();
        h_tracksTrue_Eta_p[id][k]->Write();
        h_tracksRec_Eta_p[id][k]->Write();
      }
    }
  fOutput->Write();
  fOutput->Close();
  delete fOutput;
}

