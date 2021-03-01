#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))


void pidreso_Pythia(
                            TString inputFileName   = "file.root",
                            TString suffix          = "pdf",
                            TString addLabel        = "",
                            Bool_t properFit        = kTRUE,
                            Int_t trackCuts         = 0
                            
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  SetPlotStyle();

  TString dateForOutput             = ReturnDateStringForOutput();
  TString readTrackClass            = "";
  TString writeLabel                = "";
  TString labelPlotCuts             = "";
  if (trackCuts == 2){
    readTrackClass                  = "LI3";
    addLabel                        = addLabel+"-LI3";
    writeLabel                      = "LI3";
    labelPlotCuts                   = "#geq 3 tracker hits";
  }

  TString outputDir                 = Form("plots/%s/%s",dateForOutput.Data(),addLabel.Data());
  TString outputDirBetaRes          = Form("plots/%s/%s/BetaRes",dateForOutput.Data(),addLabel.Data());
  TString outputDirBeta             = Form("plots/%s/%s/Beta",dateForOutput.Data(),addLabel.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  gSystem->Exec("mkdir -p "+outputDirBetaRes);
  gSystem->Exec("mkdir -p "+outputDirBeta);
  TString detLabel = "";
  
  TString collisionSystem = "Pythia 6, e+p, 10+250 GeV";
  TString pTHard = "";
  Int_t nLinesCol = 1;
  if (addLabel.Contains("pTHard5GeV")) {
    pTHard = "#it{p}_{T}^{hard} #geq 5 GeV/#it{c}";
    nLinesCol++;
  }


  Bool_t enablePlot[20]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  Int_t nActiveEta            = 8;
  Double_t maxBetaSigma       = 0.05;
  if (addLabel.Contains("LBLwithLGAD") ){
    detLabel  = "LBL+TTL";
    enablePlot[19] = 0;
    enablePlot[18] = 0;
    enablePlot[0] = 0;
    enablePlot[1] = 0;
    enablePlot[2] = 0;
    nActiveEta    = 16;
  } else if (addLabel.Contains("LBLwithACLGAD") ){
    detLabel  = "LBL+TTL(AC-LGAD)";
    enablePlot[19] = 0;
    enablePlot[18] = 0;
    enablePlot[0] = 0;
    enablePlot[1] = 0;
    enablePlot[2] = 0;
    nActiveEta    = 16;
  } else if (addLabel.Contains("LBLwithFTTLS2LC-ETTL-CTTL") ){
    detLabel  = "LBL+TTL(2 l's)";
    enablePlot[19] = 0;
    enablePlot[18] = 0;
    enablePlot[0] = 0;
    enablePlot[1] = 0;
    enablePlot[2] = 0;
    nActiveEta    = 16;
  } else if (addLabel.Contains("LBLwithFTTLS3LVC-ETTLLC-CTTLLC") ){
    detLabel  = "LBL+TTL(1.3mm)";
    enablePlot[19] = 0;
    enablePlot[18] = 0;
    enablePlot[0] = 0;
    enablePlot[1] = 0;
    enablePlot[2] = 0;
    nActiveEta    = 16;
  } else if (addLabel.Contains("LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1") ){
    detLabel  = "LBL+TTL(1 l b. ECal)";
    enablePlot[19] = 0;
    enablePlot[18] = 0;
    enablePlot[0] = 0;
    enablePlot[1] = 0;
    enablePlot[2] = 0;
    nActiveEta    = 16;
  } else if (addLabel.Contains("LBLwithFTTLSE2LC-ETTL-CTTLSE1") ){
    detLabel  = "LBL+TTL(1|2 l b. ECal)";
    enablePlot[19] = 0;
    enablePlot[18] = 0;
    enablePlot[0] = 0;
    enablePlot[1] = 0;
    enablePlot[2] = 0;
    nActiveEta    = 16;
  } else if (addLabel.Contains("LANLwithLGAD") ){
    detLabel  = "LANL+TTL";
    enablePlot[19] = 0;
    enablePlot[18] = 0;
    enablePlot[0] = 0;
    enablePlot[1] = 0;
    enablePlot[2] = 0;
    nActiveEta    = 16;
  } else if (addLabel.Contains("LANLwithACLGAD") ){
    detLabel  = "LANL+TTL(AC-LGAD)";
    enablePlot[19] = 0;
    enablePlot[18] = 0;
    enablePlot[0] = 0;
    enablePlot[1] = 0;
    enablePlot[2] = 0;
    nActiveEta    = 16;
  } else if (addLabel.Contains("LANLwithFTTLS3LVC-ETTLLC-CTTLLC") ){
    detLabel  = "LANL+TTL(1.3mm)";
    enablePlot[19] = 0;
    enablePlot[18] = 0;
    enablePlot[0] = 0;
    enablePlot[1] = 0;
    enablePlot[2] = 0;
    nActiveEta    = 16;
  }
  if (addLabel.Contains("Hist")){
    maxBetaSigma         = 0.05;
  }
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;

  TString partName[6]                   = {"All", "Pion", "Kaon", "Proton", "Electron", "Muon"};
  TString partLabel[6]                  = {"(h/e)^{#pm}", "#pi^{#pm}", "K^{#pm}", "p/#bar{p}", "e^{#pm}", "#mu^{#pm}"};
  //************************** Read data **************************************************
  const Int_t nPt                  = 13;
  const static Double_t partPt[]     = {0., 0.5, 1.0,  2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 
                                      9.0, 10.0, 15.0, 20.0};
//   const Int_t nPt                  = 20;
//   const static Double_t partE[]     = {0., 0.5, 1., 1.5, 2., 2.5, 3.0, 3.5, 4.0, 4.5, 
//                                       5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
//                                       10.};
  const Int_t nEta                = 19;                                        
  Double_t partEta[20]            = { -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.2, -1.0, -0.6, -0.2, 
                                       0.2, 0.6, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};

  const Int_t rebinEta_BetaResol[20]   = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 
                                         2, 2, 2, 2, 2, 2, 2, 2, 4,  1};
  Float_t betaResolFit[20]    = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};                                    
  
  Color_t colorEta[20]         = {kBlue+1, kBlue-6, kViolet+2, kViolet-5, kMagenta+1, kMagenta-6, kPink-5, kPink-9, kRed+1, kRed-7, 
                                kOrange+7, kOrange, kYellow-6, kSpring+5, kGreen+1, kGreen-5, kCyan+1, kCyan+3, kAzure+2, kBlack };
  Style_t markerStyleEta[20]   = {24, 25, 27, 28, 30, 42, 46, 24, 25, 27, 
                                 28, 30, 42, 46, 24, 25, 27, 28, 30, 20};
  Size_t markerSizeEta[20]     = {1.5, 1.4, 1.9, 1.5, 1.8, 1.8, 1.5, 1.5, 1.4, 1.9,
                                 1.5, 1.8, 1.8, 1.5, 1.5, 1.4, 1.9, 1.5, 1.8, 1.5 };
  Color_t colorPID[6]          = {kBlack, kRed+1, kGreen+2, kCyan+2, kBlue+1, kOrange};
  Style_t markerStylePID[6]    = {20, 24, 25, 27, 28, 30};
  Size_t markerSizePID[6]      = {1.5, 1.4, 1.5, 1.9, 1.6, 1.8};
  
                                 
                                 
  TH1D* h_mean_p_betaReso[6][nEta+1]            = {{NULL}};
  TH1D* h_sigma_p_betaReso[6][nEta+1]           = {{NULL}};  
  TH1D* h_beta_reso_bins[6][nEta+1][nPt]        = {{{NULL}}};
  TF1* fit_tracks_reso_bins[6][nEta+1][nPt]     = {{{NULL}}};

  TH1D* h_mean_p_betaResoAEMC[6][nEta+1]        = {{NULL}};
  TH1D* h_sigma_p_betaResoAEMC[6][nEta+1]       = {{NULL}};  
  TH1D* h_beta_reso_binsAEMC[6][nEta+1][nPt]    = {{{NULL}}};
  TF1* fit_tracks_reso_binsAEMC[6][nEta+1][nPt] = {{{NULL}}};
  
  TH2F* h_beta_reso[6][nEta+1]            = {{NULL}};
  TH2F* h_beta_resoAEMC[6][nEta+1]        = {{NULL}};
  TH2F* h_beta_p[6][nEta+1]               = {{NULL}};
  TH2F* h_betaGen_p[6][nEta+1]            = {{NULL}};
  TH2F* h_betaTrack_p[6][nEta+1]          = {{NULL}};
  TH2F* h_beta_p_AEMC[6][nEta+1]          = {{NULL}};
  TH2F* h_betaGen_p_AEMC[6][nEta+1]       = {{NULL}};
  TH2F* h_betaTrack_p_AEMC[6][nEta+1]     = {{NULL}};
  
  TFile* inputFile  = new TFile(inputFileName.Data());
  TH1D* histNEvents = (TH1D*)inputFile->Get("nEvents");
  Long_t nEvents    = histNEvents->GetBinContent(1);
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    for (Int_t pid = 0; pid < 6; pid++){    
      h_beta_reso[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_Res_InvSmearBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_beta_reso[pid][iEta]->Sumw2();
      h_mean_p_betaReso[pid][iEta] = new TH1D(Form("histBetaResol%s_%s_mean_%d", readTrackClass.Data(), partName[pid].Data(), iEta), 
                                            ";#it{p}^{MC} (GeV/#it{c}); #LT (1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC} #GT",
                                            nPt, partPt);
      h_sigma_p_betaReso[pid][iEta] = new TH1D(Form("histBetaResol%s_%s_sigma_%d", readTrackClass.Data(), partName[pid].Data(), iEta), 
                                              ";#it{p}^{MC} (GeV/#it{c}); #sigma( (1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC} )", 
                                              nPt, partPt);
      h_beta_resoAEMC[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_Res_InvSmearBeta%sAEMC_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_beta_resoAEMC[pid][iEta]->Sumw2();
      h_mean_p_betaResoAEMC[pid][iEta] = new TH1D(Form("histBetaResol%sAEMC_%s_mean_%d", readTrackClass.Data(), partName[pid].Data(), iEta), 
                                            ";#it{p}^{MC} (GeV/#it{c}); #LT (1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC} #GT",
                                            nPt, partPt);
      h_sigma_p_betaResoAEMC[pid][iEta] = new TH1D(Form("histBetaResol%sAEMC_%s_sigma_%d", readTrackClass.Data(), partName[pid].Data(), iEta), 
                                              ";#it{p}^{MC} (GeV/#it{c}); #sigma( (1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC} )", 
                                              nPt, partPt);
      h_betaTrack_p[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_betaTrack_p[pid][iEta]->Sumw2();
      h_betaGen_p[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvGenBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_betaGen_p[pid][iEta]->Sumw2();
      h_beta_p[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvSmearBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_beta_p[pid][iEta]->Sumw2();
      h_betaTrack_p_AEMC[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvBeta%sAEMC_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_betaTrack_p_AEMC[pid][iEta]->Sumw2();
      h_betaGen_p_AEMC[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvGenBeta%sAEMC_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_betaGen_p_AEMC[pid][iEta]->Sumw2();
      h_beta_p_AEMC[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvSmearBeta%sAEMC_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_beta_p_AEMC[pid][iEta]->Sumw2();
    }
  }
  
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput)std::cout << Form("%1.1f<#eta<%1.1f",etaMin,etaMax)  << std::endl;
    for (Int_t iPt = 0; iPt < nPt; iPt++){
      for (Int_t pid = 0; pid < 6; pid++){
  //       std::cout << partPt[iPt] << "\t" << h_beta_reso[iEta]->GetXaxis()->FindBin(partPt[iPt]) << std::endl;
        h_beta_reso_bins[pid][iEta][iPt] = (TH1D*)h_beta_reso[pid][iEta]->ProjectionY(Form("projection_%s_dBeta_%d_%d",partName[pid].Data(),iPt, iEta), 
                                                                          h_beta_reso[pid][iEta]->GetXaxis()->FindBin(partPt[iPt]), h_beta_reso[pid][iEta]->GetXaxis()->FindBin(partPt[iPt+1]),"e"); 
        h_beta_reso_bins[pid][iEta][iPt]->Rebin(rebinEta_BetaResol[iEta]);
        h_beta_reso_bins[pid][iEta][iPt]->GetXaxis()->SetRangeUser(-0.1,0.1);
        if(h_beta_reso_bins[pid][iEta][iPt]->GetEntries() > 10 && properFit){
            fit_tracks_reso_bins[pid][iEta][iPt] = new TF1(Form("fitGauss%s%d_%d",partName[pid].Data() ,iPt, iEta),  "gaus", -betaResolFit[iEta],betaResolFit[iEta]);
            fit_tracks_reso_bins[pid][iEta][iPt]->SetParameter(0,h_beta_reso_bins[pid][iEta][iPt]->GetMaximum() );
            fit_tracks_reso_bins[pid][iEta][iPt]->SetParLimits(0,0.9*h_beta_reso_bins[pid][iEta][iPt]->GetMaximum(),4*h_beta_reso_bins[pid][iEta][iPt]->GetMaximum() );
            h_beta_reso_bins[pid][iEta][iPt]->Fit(fit_tracks_reso_bins[pid][iEta][iPt],"QNRME+","",-betaResolFit[iEta],betaResolFit[iEta]);
        }
        h_beta_reso_binsAEMC[pid][iEta][iPt] = (TH1D*)h_beta_resoAEMC[pid][iEta]->ProjectionY(Form("projectionAEMC_%s_dBeta_%d_%d",partName[pid].Data(),iPt, iEta), 
                                                                          h_beta_resoAEMC[pid][iEta]->GetXaxis()->FindBin(partPt[iPt]), h_beta_resoAEMC[pid][iEta]->GetXaxis()->FindBin(partPt[iPt+1]),"e"); 
        h_beta_reso_binsAEMC[pid][iEta][iPt]->Rebin(rebinEta_BetaResol[iEta]);
        h_beta_reso_binsAEMC[pid][iEta][iPt]->GetXaxis()->SetRangeUser(-0.1,0.1);
        if(h_beta_reso_binsAEMC[pid][iEta][iPt]->GetEntries() > 10 && properFit){
            fit_tracks_reso_binsAEMC[pid][iEta][iPt] = new TF1(Form("fitGaussAEMC%s%d_%d",partName[pid].Data() ,iPt, iEta),  "gaus", -betaResolFit[iEta],betaResolFit[iEta]);
            fit_tracks_reso_binsAEMC[pid][iEta][iPt]->SetParameter(0,h_beta_reso_binsAEMC[pid][iEta][iPt]->GetMaximum() );
            fit_tracks_reso_binsAEMC[pid][iEta][iPt]->SetParLimits(0,0.9*h_beta_reso_binsAEMC[pid][iEta][iPt]->GetMaximum(),4*h_beta_reso_bins[pid][iEta][iPt]->GetMaximum() );
            h_beta_reso_binsAEMC[pid][iEta][iPt]->Fit(fit_tracks_reso_binsAEMC[pid][iEta][iPt],"QNRME+","",-betaResolFit[iEta],betaResolFit[iEta]);
        }
      }
      
      if(iEta == 0 && debugOutput){
        std::cout <<h_mean_p_betaReso[0][iEta]->GetBinCenter(iPt+1) << "\t"<<iPt+1<<  "\t"<< h_beta_reso_bins[0][iEta][iPt]->GetMean() <<  "\t"<< h_beta_reso_bins[0][iEta][iPt]->GetRMS() << std::endl;
        if (fit_tracks_reso_bins[0][iEta][iPt]) std::cout<< fit_tracks_reso_bins[0][iEta][iPt]->GetParameter(1) << "\t" << fit_tracks_reso_bins[0][iEta][iPt]->GetParameter(2) << std::endl;
      }
      if (properFit){
        for (Int_t pid = 0; pid < 6; pid++){
          if (fit_tracks_reso_bins[pid][iEta][iPt]){
            h_mean_p_betaReso[pid][iEta]->SetBinContent(iPt+1,fit_tracks_reso_bins[pid][iEta][iPt]->GetParameter(1));
            h_mean_p_betaReso[pid][iEta]->SetBinError(iPt+1,fit_tracks_reso_bins[pid][iEta][iPt]->GetParError(1));
            h_sigma_p_betaReso[pid][iEta]->SetBinContent(iPt+1,fit_tracks_reso_bins[pid][iEta][iPt]->GetParameter(2));
            h_sigma_p_betaReso[pid][iEta]->SetBinError(iPt+1,fit_tracks_reso_bins[pid][iEta][iPt]->GetParError(2));
          } else {
            h_mean_p_betaReso[pid][iEta]->SetBinContent(iPt+1,-100000);
            h_mean_p_betaReso[pid][iEta]->SetBinError(iPt+1,0);
            h_sigma_p_betaReso[pid][iEta]->SetBinContent(iPt+1,-1000);
            h_sigma_p_betaReso[pid][iEta]->SetBinError(iPt+1,0);         
          }
          if (fit_tracks_reso_binsAEMC[pid][iEta][iPt]){
            h_mean_p_betaResoAEMC[pid][iEta]->SetBinContent(iPt+1,fit_tracks_reso_binsAEMC[pid][iEta][iPt]->GetParameter(1));
            h_mean_p_betaResoAEMC[pid][iEta]->SetBinError(iPt+1,fit_tracks_reso_binsAEMC[pid][iEta][iPt]->GetParError(1));
            h_sigma_p_betaResoAEMC[pid][iEta]->SetBinContent(iPt+1,fit_tracks_reso_binsAEMC[pid][iEta][iPt]->GetParameter(2));
            h_sigma_p_betaResoAEMC[pid][iEta]->SetBinError(iPt+1,fit_tracks_reso_binsAEMC[pid][iEta][iPt]->GetParError(2));
          } else {
            h_mean_p_betaResoAEMC[pid][iEta]->SetBinContent(iPt+1,-100000);
            h_mean_p_betaResoAEMC[pid][iEta]->SetBinError(iPt+1,0);
            h_sigma_p_betaResoAEMC[pid][iEta]->SetBinContent(iPt+1,-1000);
            h_sigma_p_betaResoAEMC[pid][iEta]->SetBinError(iPt+1,0);         
          }
        }
      } else {
        for (Int_t pid = 0; pid < 6; pid++){
          if(h_beta_reso_bins[pid][iEta][iPt]->GetEntries() > 10){
            h_mean_p_betaReso[pid][iEta]->SetBinContent(iPt+1,h_beta_reso_bins[pid][iEta][iPt]->GetMean());
            h_mean_p_betaReso[pid][iEta]->SetBinError(iPt+1,h_beta_reso_bins[pid][iEta][iPt]->GetMeanError());
            h_sigma_p_betaReso[pid][iEta]->SetBinContent(iPt+1,h_beta_reso_bins[pid][iEta][iPt]->GetRMS());
            h_sigma_p_betaReso[pid][iEta]->SetBinError(iPt+1,h_beta_reso_bins[pid][iEta][iPt]->GetRMSError());      
//       h_sigma_p_betaReso[iEta]->SetBinContent(iPt+1,h_beta_reso_bins[iEta][iPt]->GetStdDev());
//       h_sigma_p_betaReso[iEta]->SetBinError(iPt+1,h_beta_reso_bins[iEta][iPt]->GetStdDevError());
          } else {
            h_mean_p_betaReso[pid][iEta]->SetBinContent(iPt+1,-100000);
            h_mean_p_betaReso[pid][iEta]->SetBinError(iPt+1,0);
            h_sigma_p_betaReso[pid][iEta]->SetBinContent(iPt+1,-1000);
            h_sigma_p_betaReso[pid][iEta]->SetBinError(iPt+1,0);         
          }
          if(h_beta_reso_binsAEMC[pid][iEta][iPt]->GetEntries() > 10){
            h_mean_p_betaResoAEMC[pid][iEta]->SetBinContent(iPt+1,h_beta_reso_binsAEMC[pid][iEta][iPt]->GetMean());
            h_mean_p_betaResoAEMC[pid][iEta]->SetBinError(iPt+1,h_beta_reso_binsAEMC[pid][iEta][iPt]->GetMeanError());
            h_sigma_p_betaResoAEMC[pid][iEta]->SetBinContent(iPt+1,h_beta_reso_binsAEMC[pid][iEta][iPt]->GetRMS());
            h_sigma_p_betaResoAEMC[pid][iEta]->SetBinError(iPt+1,h_beta_reso_binsAEMC[pid][iEta][iPt]->GetRMSError());      
//       h_sigma_p_betaReso[iEta]->SetBinContent(iPt+1,h_beta_reso_bins[iEta][iPt]->GetStdDev());
//       h_sigma_p_betaReso[iEta]->SetBinError(iPt+1,h_beta_reso_bins[iEta][iPt]->GetStdDevError());
          } else {
            h_mean_p_betaResoAEMC[pid][iEta]->SetBinContent(iPt+1,-100000);
            h_mean_p_betaResoAEMC[pid][iEta]->SetBinError(iPt+1,0);
            h_sigma_p_betaResoAEMC[pid][iEta]->SetBinContent(iPt+1,-1000);
            h_sigma_p_betaResoAEMC[pid][iEta]->SetBinError(iPt+1,0);         
          }
        }  
      }
    }
  }

  
  
  
  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
  split_canvas(cPNG, "cPNG1", nPt);
  
  TCanvas* cSingle2D = new TCanvas("cSingle2D","",0,0,1000,800);
  DrawGammaCanvasSettings( cSingle2D, 0.1, 0.12, 0.025, 0.105);
  cSingle2D->SetLogz();
  
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    Int_t padnum =0;
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput) std::cout << Form("%1.1f<#eta<%1.1f",etaMin,etaMax)  << std::endl;
      
    for (Int_t pid = 0; pid < 6; pid++){
      padnum =0;
      for(Int_t iPt=0; iPt< nPt; iPt++){
        cPNG->cd(padnum+1);
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNG->cd(padnum+1)->SetLogy(1);
        SetStyleHistoTH1ForGraphs(h_beta_reso_bins[pid][iEta][iPt], "(1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC}", "counts",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        h_beta_reso_bins[pid][iEta][iPt]->GetYaxis()->SetRangeUser(0.8,h_beta_reso_bins[pid][iEta][iPt]->GetMaximum()*2);
        
        DrawGammaSetMarker(h_beta_reso_bins[pid][iEta][iPt], 20, 1.5, kBlack, kBlack);      
        h_beta_reso_bins[pid][iEta][iPt]->Draw("p,e");
        
        if (fit_tracks_reso_bins[pid][iEta][iPt] && properFit){
          DrawGammaSetMarkerTF1( fit_tracks_reso_bins[pid][iEta][iPt], 7, 2, kBlue+1);
          fit_tracks_reso_bins[pid][iEta][iPt]->Draw("same");
        }
        
        if (padnum == 0){
          drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        
        drawLatexAdd(Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f ",partPt[iPt],partPt[iPt+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(iEta == 0 && debugOutput) std::cout << Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f \t",partPt[iPt],partPt[iPt+1]) << h_mean_p_betaReso[pid][iEta]->GetBinCenter(iPt+1)<< "\t"<< h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1) << std::endl;
        DrawGammaLines(h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1), 
                      h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1), 
                      0., 0.05*h_beta_reso_bins[pid][iEta][iPt]->GetMaximum(), 1, kGray+2, 2);
        DrawGammaLines(h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1)-h_sigma_p_betaReso[pid][iEta]->GetBinContent(iPt+1), 
                      h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1)-h_sigma_p_betaReso[pid][iEta]->GetBinContent(iPt+1), 
                      0., 0.05*h_beta_reso_bins[pid][iEta][iPt]->GetMaximum(), 1, kRed+2, 2);
        DrawGammaLines(h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1)+h_sigma_p_betaReso[pid][iEta]->GetBinContent(iPt+1), 
                      h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1)+h_sigma_p_betaReso[pid][iEta]->GetBinContent(iPt+1), 
                      0., 0.05*h_beta_reso_bins[pid][iEta][iPt]->GetMaximum(), 1, kRed+2, 2);

        padnum++;
      }
      cPNG->Print(Form("%s/Beta_%s_Resolution_%d_%d.%s", outputDirBetaRes.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      
      padnum =0;
      for(Int_t iPt=0; iPt< nPt; iPt++){
        cPNG->cd(padnum+1);
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNG->cd(padnum+1)->SetLogy(1);
        SetStyleHistoTH1ForGraphs(h_beta_reso_binsAEMC[pid][iEta][iPt], "(1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC}", "counts",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        h_beta_reso_binsAEMC[pid][iEta][iPt]->GetYaxis()->SetRangeUser(0.8,h_beta_reso_bins[pid][iEta][iPt]->GetMaximum()*2);
        
        DrawGammaSetMarker(h_beta_reso_binsAEMC[pid][iEta][iPt], 20, 1.5, kBlack, kBlack);      
        h_beta_reso_binsAEMC[pid][iEta][iPt]->Draw("p,e");
        
        if (fit_tracks_reso_binsAEMC[pid][iEta][iPt] && properFit){
          DrawGammaSetMarkerTF1( fit_tracks_reso_binsAEMC[pid][iEta][iPt], 7, 2, kBlue+1);
          fit_tracks_reso_binsAEMC[pid][iEta][iPt]->Draw("same");
        }
        
        if (padnum == 0){
          drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        
        drawLatexAdd(Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f ",partPt[iPt],partPt[iPt+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(iEta == 0 && debugOutput) std::cout << Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f \t",partPt[iPt],partPt[iPt+1]) << h_mean_p_betaReso[pid][iEta]->GetBinCenter(iPt+1)<< "\t"<< h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1) << std::endl;
        DrawGammaLines(h_mean_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1), 
                      h_mean_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1), 
                      0., 0.05*h_beta_reso_binsAEMC[pid][iEta][iPt]->GetMaximum(), 1, kGray+2, 2);
        DrawGammaLines(h_mean_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1)-h_sigma_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1), 
                      h_mean_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1)-h_sigma_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1), 
                      0., 0.05*h_beta_reso_binsAEMC[pid][iEta][iPt]->GetMaximum(), 1, kRed+2, 2);
        DrawGammaLines(h_mean_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1)+h_sigma_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1), 
                      h_mean_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1)+h_sigma_p_betaResoAEMC[pid][iEta]->GetBinContent(iPt+1), 
                      0., 0.05*h_beta_reso_binsAEMC[pid][iEta][iPt]->GetMaximum(), 1, kRed+2, 2);

        padnum++;
      }
      cPNG->Print(Form("%s/Beta_%s_ResolutionWithHitAfterECal_%d_%d.%s", outputDirBetaRes.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
  }
  split_canvas(cPNG, "cPNG1", nEta+1);
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.11, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      SetStyleHistoTH2ForGraphs(h_beta_reso[pid][iEta], "#it{p} (GeV/#it{c})", "(1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      h_beta_reso[pid][iEta]->Scale(1./nEvents);
      h_beta_reso[pid][iEta]->GetYaxis()->SetRangeUser(-0.02,0.02);
      h_beta_reso[pid][iEta]->GetZaxis()->SetRangeUser(0.9/nEvents,1./nEvents*1e4);
      h_beta_reso[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
    }
    cPNG->Print(Form("%s/BetaP_%s_Resolution.%s", outputDirBetaRes.Data(), partName[pid].Data(), suffix.Data()));
  }
  
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.11, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      SetStyleHistoTH2ForGraphs(h_beta_resoAEMC[pid][iEta], "#it{p} (GeV/#it{c})", "(1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      h_beta_resoAEMC[pid][iEta]->Scale(1./nEvents);
      h_beta_resoAEMC[pid][iEta]->GetYaxis()->SetRangeUser(-0.02,0.02);
      h_beta_resoAEMC[pid][iEta]->GetZaxis()->SetRangeUser(0.9/nEvents,1./nEvents*1e4);
      h_beta_resoAEMC[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iEta == 0)  drawLatexAdd("w/ hit after ECal",0.16,0.86-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
    }
    cPNG->Print(Form("%s/BetaP_%s_ResolutionWithHitAfterECal.%s", outputDirBetaRes.Data(), partName[pid].Data(), suffix.Data()));
  }
  
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      SetStyleHistoTH2ForGraphs(h_betaGen_p[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{MC}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      h_betaGen_p[pid][iEta]->Scale(1./nEvents);
//       h_betaGen_p[pid][iEta]->GetYaxis()->SetRangeUser(-0.02,0.02);
      h_betaGen_p[pid][iEta]->GetZaxis()->SetRangeUser(0.9/nEvents,1./nEvents*1e4);
      h_betaGen_p[pid][iEta]->GetXaxis()->SetRangeUser(0,15.);
      h_betaGen_p[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
      
      cSingle2D->cd();
        h_betaGen_p[pid][iEta]->Draw("colz");
  
  
        drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(collisionSystem,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.91-(nLinesCol*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/GenBetaP_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
    cPNG->Print(Form("%s/GenBetaP_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
    
    
    
  }
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      SetStyleHistoTH2ForGraphs(h_betaGen_p_AEMC[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{MC}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      h_betaGen_p_AEMC[pid][iEta]->Scale(1./nEvents);
      h_betaGen_p_AEMC[pid][iEta]->GetZaxis()->SetRangeUser(0.9/nEvents,1./nEvents*1e4);
      h_betaGen_p_AEMC[pid][iEta]->GetXaxis()->SetRangeUser(0,15.);
      h_betaGen_p_AEMC[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iEta == 0)  drawLatexAdd("w/ hit after ECal",0.16,0.86-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
      
      cSingle2D->cd();
        h_betaGen_p_AEMC[pid][iEta]->Draw("colz");
        
        drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(collisionSystem,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.91-(nLinesCol*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("w/ hit after ECal",0.16,0.91-((nLinesCol+1)*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/GenBetaPWithHitAfterECal_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
    cPNG->Print(Form("%s/GenBetaPWithHitAfterECal_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
  }
 
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      SetStyleHistoTH2ForGraphs(h_betaTrack_p[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{path}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      h_betaTrack_p[pid][iEta]->Scale(1./nEvents);
      h_betaTrack_p[pid][iEta]->GetXaxis()->SetRangeUser(0,15.);
      h_betaTrack_p[pid][iEta]->GetZaxis()->SetRangeUser(0.9/nEvents,1./nEvents*1e4);
      h_betaTrack_p[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;

      cSingle2D->cd();
        h_betaTrack_p[pid][iEta]->Draw("colz");
        
        drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(collisionSystem,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.91-(nLinesCol*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/TrackBetaP_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    }
    cPNG->Print(Form("%s/TrackBetaP_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
  }
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      SetStyleHistoTH2ForGraphs(h_betaTrack_p_AEMC[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{path}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      h_betaTrack_p_AEMC[pid][iEta]->Scale(1./nEvents);
      h_betaTrack_p_AEMC[pid][iEta]->GetXaxis()->SetRangeUser(0,15.);
      h_betaTrack_p_AEMC[pid][iEta]->GetZaxis()->SetRangeUser(0.9/nEvents,1./nEvents*1e4);
      h_betaTrack_p_AEMC[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iEta == 0)  drawLatexAdd("w/ hit after ECal",0.16,0.86-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
      
      cSingle2D->cd();
        h_betaTrack_p_AEMC[pid][iEta]->Draw("colz");

        drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(collisionSystem,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.91-(nLinesCol*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("w/ hit after ECal",0.16,0.91-((nLinesCol+1)*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/TrackBetaPWithHitAfterECal_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
    cPNG->Print(Form("%s/TrackBetaPWithHitAfterECal_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
  }

  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      SetStyleHistoTH2ForGraphs(h_beta_p[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{rec}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      h_beta_p[pid][iEta]->Scale(1./nEvents);
      h_beta_p[pid][iEta]->GetXaxis()->SetRangeUser(0,15.);
      h_beta_p[pid][iEta]->GetZaxis()->SetRangeUser(0.9/nEvents,1./nEvents*1e4);
      h_beta_p[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;

      cSingle2D->cd();
        h_beta_p[pid][iEta]->Draw("colz");
        
        drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(collisionSystem,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.91-(nLinesCol*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/BetaP_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
    cPNG->Print(Form("%s/BetaP_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
  }
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      SetStyleHistoTH2ForGraphs(h_beta_p_AEMC[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{rec}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      h_beta_p_AEMC[pid][iEta]->Scale(1./nEvents);
      h_beta_p_AEMC[pid][iEta]->GetXaxis()->SetRangeUser(0,15.);
      h_beta_p_AEMC[pid][iEta]->GetZaxis()->SetRangeUser(0.9/nEvents,1./nEvents*1e4);
      h_beta_p_AEMC[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iEta == 0)  drawLatexAdd("w/ hit after ECal",0.16,0.86-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
      
      cSingle2D->cd();
        h_beta_p_AEMC[pid][iEta]->Draw("colz");
        
        drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(collisionSystem,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.85, 0.91-(nLinesCol*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("w/ hit after ECal",0.16,0.91-((nLinesCol+1)*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/BetaPWithHitAfterECal_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
    cPNG->Print(Form("%s/BetaPWithHitAfterECal_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
  }
 
  // 1D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,800);
  DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.045, 0.105);
  // cReso->SetLogz();

  TH2F* histoDummyBetaResMean   = new TH2F("histoDummyBetaResMean","histoDummyBetaResMean",1000,0, 20,1000,-0.01, 0.025);
  SetStyleHistoTH2ForGraphs(histoDummyBetaResMean, "#it{p}^{MC} (GeV/#it{c})","#LT (1/#beta^{rec} - 1/#beta^{MC}) / 1/#beta^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyBetaResMean->GetXaxis()->SetNoExponent();
  histoDummyBetaResMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyBetaResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyBetaResSigma   = new TH2F("histoDummyBetaResSigma","histoDummyBetaResSigma",1000,0, 20,1000,-0.0, maxBetaSigma);
  SetStyleHistoTH2ForGraphs(histoDummyBetaResSigma, "#it{p}^{MC} (GeV/#it{c})","#sigma((1/#beta^{rec} - 1/#beta^{MC}) / 1/#beta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyBetaResSigma->GetXaxis()->SetNoExponent();
  histoDummyBetaResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyBetaResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);

  TLegend* legendPtResM  = GetAndSetLegend2(0.14, 0.91-(nActiveEta/3*0.75*textSizeLabelsRel), 0.6, 0.91,0.75*textSizeLabelsPixel, 3, "", 43, 0.2);
  TLegend* legendPtResPID  = GetAndSetLegend2(0.14, 0.91-(3*0.85*textSizeLabelsRel), 0.35, 0.91,0.85*textSizeLabelsPixel, 2, "", 43, 0.25);
  
  for (Int_t pid = 0; pid < 6; pid++){
    histoDummyBetaResMean->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_mean_p_betaReso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_mean_p_betaReso[pid][iEta]->Draw("same,p");
      if (pid == 0){
        if (iEta == nEta )
          legendPtResM->AddEntry(h_mean_p_betaReso[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
        else 
          legendPtResM->AddEntry(h_mean_p_betaReso[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
      }
    }
    legendPtResM->Draw();
    DrawGammaLines(0, 20, 0., 0., 2, kGray+2, 7);
   
    drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s",partLabel[pid].Data(), detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolution_%s_Mean_p.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
    
    histoDummyBetaResSigma->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_sigma_p_betaReso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_sigma_p_betaReso[pid][iEta]->Draw("same,p");
    }
    legendPtResM->Draw();
    drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s",partLabel[pid].Data(), detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolution_%s_Sigma_p.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
    
    histoDummyBetaResMean->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_mean_p_betaResoAEMC[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_mean_p_betaResoAEMC[pid][iEta]->Draw("same,p");
    }
    legendPtResM->Draw();
    DrawGammaLines(0, 20, 0., 0., 2, kGray+2, 7);

    drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s",partLabel[pid].Data(), detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("w/ hit after ECal",0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionWithHitAfterECal_%s_Mean_p.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
    
    histoDummyBetaResSigma->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_sigma_p_betaResoAEMC[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_sigma_p_betaResoAEMC[pid][iEta]->Draw("same,p");
    }
    legendPtResM->Draw();
    drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s",partLabel[pid].Data(), detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("w/ hit after ECal",0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionWithHitAfterECal_%s_Sigma_p.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
  }
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    histoDummyBetaResSigma->Draw();
    legendPtResPID->Clear();
    for (Int_t pid =0; pid < 6; pid++){
      
      DrawGammaSetMarker(h_sigma_p_betaReso[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_sigma_p_betaReso[pid][iEta]->Draw("same,p");      
      legendPtResPID->AddEntry(h_sigma_p_betaReso[pid][iEta],partLabel[pid].Data(),"p");
    }
    legendPtResPID->Draw();

    drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s",etaMin,etaMax, detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionPID_Sigma_p_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
    histoDummyBetaResSigma->Draw();
    for (Int_t pid =0; pid < 6; pid++){
      DrawGammaSetMarker(h_sigma_p_betaResoAEMC[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_sigma_p_betaResoAEMC[pid][iEta]->Draw("same,p");      
    }
    legendPtResPID->Draw();
    drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s",etaMin,etaMax, detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("w/ hit after ECal",0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionWithHitAfterECalPID_Sigma_p_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    histoDummyBetaResMean->Draw();
    for (Int_t pid =0; pid < 6; pid++){
      
      DrawGammaSetMarker(h_mean_p_betaReso[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_mean_p_betaReso[pid][iEta]->Draw("same,p");      
    }
    legendPtResPID->Draw();

    drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s",etaMin,etaMax, detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionPID_Mean_p_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
    histoDummyBetaResMean->Draw();
    for (Int_t pid =0; pid < 6; pid++){
      DrawGammaSetMarker(h_mean_p_betaResoAEMC[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_mean_p_betaResoAEMC[pid][iEta]->Draw("same,p");      
    }
    legendPtResPID->Draw();
    drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s",etaMin,etaMax, detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("w/ hit after ECal",0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionWithHitAfterECalPID_Mean_p_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  }
  
  
  TFile* outputFile  = new TFile(inputFileName.Data(),"UPDATE");
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    if (properFit){
      for (Int_t pid = 0; pid < 6; pid++){
        h_mean_p_betaReso[pid][iEta]->Write(Form("histBetaResol%s_%s_FitMean_%d", readTrackClass.Data(), partName[pid].Data(),iEta),TObject::kOverwrite);
        h_sigma_p_betaReso[pid][iEta]->Write(Form("histBetaResol%s_%s_FitSigma_%d",readTrackClass.Data(),  partName[pid].Data(), iEta),TObject::kOverwrite);
        h_mean_p_betaResoAEMC[pid][iEta]->Write(Form("histBetaResol%sAEMC_%s_FitMean_%d", readTrackClass.Data(), partName[pid].Data(),iEta),TObject::kOverwrite);
        h_sigma_p_betaResoAEMC[pid][iEta]->Write(Form("histBetaResol%sAEMC_%s_FitSigma_%d", readTrackClass.Data(),  partName[pid].Data(), iEta),TObject::kOverwrite);
      }
    } else {
      for (Int_t pid = 0; pid < 6; pid++){
        h_mean_p_betaReso[pid][iEta]->Write(Form("histBetaResol%s_%s_mean_%d", readTrackClass.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
        h_sigma_p_betaReso[pid][iEta]->Write(Form("histBetaResol%s_%s_sigma_%d", readTrackClass.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);      
        h_mean_p_betaResoAEMC[pid][iEta]->Write(Form("histBetaResol%sAEMC_%s_mean_%d", readTrackClass.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
        h_sigma_p_betaResoAEMC[pid][iEta]->Write(Form("histBetaResol%sAEMC_%s_sigma_%d", readTrackClass.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);      
      }
    }
  }
  outputFile->Write();
  outputFile->Close();
}
